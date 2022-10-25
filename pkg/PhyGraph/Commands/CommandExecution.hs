{- |
Module      :  CommandExecution.hs
Description :  Module to coordinate command execution
Copyright   :  (c) 2021 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
License     :

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

Maintainer  :  Ward Wheeler <wheeler@amnh.org>
Stability   :  unstable
Portability :  portable (I hope)

-}

{-# LANGUAGE BangPatterns #-}

module Commands.CommandExecution
  ( executeCommands
  , executeRenameReblockCommands
  , getDataListList
  ) where

import qualified Data.CSV               as CSV
import qualified Data.List              as L
import           Data.Maybe
import           Text.Read
import qualified Data.Text.Lazy         as T
import qualified Data.Vector            as V
import           Debug.Trace
import           GeneralUtilities
import           System.IO
import           Types.Types
import qualified Utilities.Utilities    as U
import qualified Data.Char              as C
import qualified Search.Build           as B
import qualified Reconciliation.ReconcileGraphs as R
import qualified Search.Refinement as REF
import qualified Search.Search as S
import qualified Commands.Transform as TRANS
import           System.Timing
import qualified Support.Support as SUP
import           Data.Char
import Graphs.GraphOperations as GO
import qualified GraphOptimization.Traversals as TRAV
import qualified Commands.Verify             as VER
import qualified Data.InfList                as IL
import           Commands.CommandUtilities



-- | executeCommands reads input files and returns raw data
-- need to close files after read
executeCommands :: GlobalSettings -> Int -> String -> ProcessedData -> ProcessedData -> [PhylogeneticGraph] -> [[VertexCost]] -> [Int] -> [PhylogeneticGraph] -> [Command] -> IO ([PhylogeneticGraph], GlobalSettings, [Int], [PhylogeneticGraph])
executeCommands globalSettings numInputFiles crossReferenceString origProcessedData processedData curGraphs pairwiseDist seedList supportGraphList commandList = do
    if null commandList then return (curGraphs, globalSettings, seedList, supportGraphList)
    else do
        let (firstOption, firstArgs) = head commandList
        
        -- skip "Read" and "Rename "commands already processed
        if firstOption == Read then error ("Read command should already have been processed: " ++ show (firstOption, firstArgs))
        else if firstOption == Rename then error ("Rename command should already have been processed: " ++ show (firstOption, firstArgs))
        else if firstOption == Reblock then error ("Reblock command should already have been processed: " ++ show (firstOption, firstArgs))
        else if firstOption == Run then error ("Run command should already have been processed: " ++ show (firstOption, firstArgs))
        
        -- other commands
        else if firstOption == Build then do
            (elapsedSeconds, newGraphList) <- timeOp $ pure $ B.buildGraph firstArgs globalSettings processedData pairwiseDist (head seedList) 
                
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : (searchData globalSettings)
            
            executeCommands (globalSettings {searchData = newSearchData}) numInputFiles crossReferenceString origProcessedData processedData (curGraphs ++ newGraphList) pairwiseDist (tail seedList) supportGraphList (tail commandList)
        
        else if firstOption == Refine then do 
            (elapsedSeconds, newGraphList) <- timeOp $ pure $ REF.refineGraph firstArgs globalSettings processedData (head seedList) curGraphs
            
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : (searchData globalSettings)   
            
            executeCommands (globalSettings {searchData = newSearchData}) numInputFiles crossReferenceString origProcessedData processedData newGraphList pairwiseDist (tail seedList) supportGraphList (tail commandList)
        
        else if firstOption == Fuse then do
            (elapsedSeconds, newGraphList) <- timeOp $ pure $ REF.fuseGraphs firstArgs globalSettings processedData (head seedList) curGraphs
            
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : (searchData globalSettings)   
            
            executeCommands (globalSettings {searchData = newSearchData}) numInputFiles crossReferenceString origProcessedData processedData newGraphList pairwiseDist (tail seedList) supportGraphList (tail commandList)
        
        else if firstOption == Report then do
            
            let doDotPDF = any (=="dotpdf") $ fmap (fmap toLower . fst) firstArgs
            let collapse' = any (=="collapse") $ fmap (fmap toLower . fst) firstArgs
            let noCollapse' = any (=="nocollapse") $ fmap (fmap toLower . fst) firstArgs

            -- set default collapse for dotPDF to True, False otherwise
            let collapse = if collapse' then True
                           else if noCollapse' then False
                           else if doDotPDF then True
                           else False

            let curGraphs' = if (not collapse) then curGraphs
                             else fmap U.collapseGraph curGraphs
            
            -- use 'temp' updated graphs s don't repeatedly add model and root complexityies
            -- reporting collapsed 
            -- reverse sorting graphs by cost
            let graphsWithUpdatedCosts = reverse (L.sortOn snd6 $ fmap (TRAV.updateGraphCostsComplexities globalSettings) curGraphs')
                reportStuff@(reportString, outFile, writeMode) = reportCommand globalSettings firstArgs numInputFiles crossReferenceString  processedData graphsWithUpdatedCosts supportGraphList pairwiseDist
            
            if null reportString then do
                executeCommands globalSettings numInputFiles crossReferenceString origProcessedData processedData curGraphs pairwiseDist seedList supportGraphList (tail commandList)
            else  do
                hPutStrLn stderr ("Report writing to " ++ outFile)
                
                if doDotPDF then do
                    let reportString' = changeDotPreamble "digraph {" "digraph G {\n\trankdir = LR;\tnode [ shape = rect];\n" reportString
                    printGraphVizDot reportString' outFile
                    executeCommands globalSettings numInputFiles crossReferenceString origProcessedData processedData curGraphs pairwiseDist seedList supportGraphList (tail commandList)

                else do
                    if outFile == "stderr" then hPutStr stderr reportString
                    else if outFile == "stdout" then putStr reportString
                    else if writeMode == "overwrite" then writeFile outFile reportString
                    else if writeMode == "append" then appendFile outFile reportString
                    else error ("Error 'read' command not properly formatted" ++ show reportStuff)
                    executeCommands globalSettings numInputFiles crossReferenceString origProcessedData processedData curGraphs pairwiseDist seedList supportGraphList (tail commandList)
        
        else if firstOption == Search then do
            (elapsedSeconds, output) <- timeOp $ S.search firstArgs globalSettings processedData pairwiseDist (head seedList) curGraphs
                --in pure result
            -- (newGraphList, serchInfoList) <- S.search firstArgs globalSettings origProcessedData processedData pairwiseDist (head seedList) curGraphs
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs (fst output) (fromIntegral $ toMilliseconds elapsedSeconds) (concat $ fmap (L.intercalate "\n") $ snd output)
            let newSearchData = searchInfo : (searchData globalSettings)
            executeCommands (globalSettings {searchData = newSearchData})  numInputFiles crossReferenceString origProcessedData processedData  (fst output) pairwiseDist (tail seedList) supportGraphList (tail commandList)
        
        else if firstOption == Select then do
            (elapsedSeconds, newGraphList) <- timeOp $ pure $ GO.selectPhylogeneticGraph firstArgs (head seedList) VER.selectArgList curGraphs
                
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : (searchData globalSettings)   
            
            executeCommands (globalSettings {searchData = newSearchData}) numInputFiles crossReferenceString origProcessedData processedData newGraphList pairwiseDist (tail seedList) supportGraphList (tail commandList)
        
        else if firstOption == Set then 
            -- if set changes graph aspects--may nned to reoptimize
            let (newGlobalSettings, newProcessedData, seedList') = setCommand firstArgs globalSettings processedData seedList
                newGraphList = if not (requireReoptimization globalSettings newGlobalSettings) then curGraphs
                               else trace ("Reoptimizing gaphs") fmap (TRAV.multiTraverseFullyLabelGraph newGlobalSettings newProcessedData True True Nothing) (fmap fst6 curGraphs)
                
                searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList 0 "No Comment"
                newSearchData = searchInfo : (searchData newGlobalSettings)   
            in
            
            executeCommands (newGlobalSettings {searchData = newSearchData}) numInputFiles crossReferenceString origProcessedData processedData newGraphList pairwiseDist seedList' supportGraphList (tail commandList)
        
        else if firstOption == Swap then do
            (elapsedSeconds, newGraphList) <- timeOp $ pure $ REF.swapMaster firstArgs globalSettings processedData (head seedList)  curGraphs
                
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : (searchData globalSettings)   
            
            executeCommands (globalSettings {searchData = newSearchData}) numInputFiles crossReferenceString origProcessedData processedData newGraphList pairwiseDist (tail seedList) supportGraphList (tail commandList)
        
        else if firstOption == Support then do
            (elapsedSeconds, newSupportGraphList) <- timeOp $ pure $ SUP.supportGraph firstArgs globalSettings processedData (head seedList)  curGraphs
                
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newSupportGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : (searchData globalSettings)   
            
            executeCommands (globalSettings {searchData = newSearchData}) numInputFiles crossReferenceString origProcessedData processedData curGraphs pairwiseDist (tail seedList) (supportGraphList ++ newSupportGraphList) (tail commandList)

        else if firstOption == Transform then do
            (elapsedSeconds, (newGS, newOrigData, newProcessedData, newGraphs)) <- timeOp $ pure $ TRANS.transform firstArgs globalSettings origProcessedData processedData (head seedList) curGraphs
                
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphs (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : (searchData globalSettings)   
            
            executeCommands (newGS {searchData = newSearchData}) numInputFiles crossReferenceString newOrigData newProcessedData newGraphs pairwiseDist (tail seedList) supportGraphList (tail commandList)

        else error ("Command " ++ (show firstOption) ++ " not recognized/implemented")

-- | makeSearchRecord take sbefore and after data of a commend and returns SearchData record
makeSearchRecord :: Instruction -> [Argument] -> [PhylogeneticGraph] -> [PhylogeneticGraph] -> Int -> String -> SearchData
makeSearchRecord firstOption firstArgs curGraphs newGraphList elapsedTime comment =
    SearchData { instruction = firstOption
               , arguments = firstArgs
               , minGraphCostIn = if null curGraphs then infinity 
                                  else minimum $ fmap snd6 curGraphs
               , maxGraphCostIn = if null curGraphs then infinity 
                                  else maximum $ fmap snd6 curGraphs
               , numGraphsIn = length curGraphs
               , minGraphCostOut = if null newGraphList then infinity 
                                   else minimum $ fmap snd6 newGraphList
               , maxGraphCostOut = if null newGraphList then infinity 
                                   else maximum $ fmap snd6 newGraphList
               , numGraphsOut = length newGraphList
               , commentString = comment
               , duration = elapsedTime
               }


-- | setCommand takes arguments to change globalSettings and multiple data aspects (e.g. 'blocks')
-- needs to be abtracted--too long
setCommand :: [Argument] -> GlobalSettings -> ProcessedData -> [Int] -> (GlobalSettings, ProcessedData, [Int])
setCommand argList globalSettings processedData inSeedList =
    let commandList = fmap (fmap C.toLower) $ filter (/= "") $ fmap fst argList
        optionList = fmap (fmap C.toLower) $ filter (/= "") $ fmap snd argList
        checkCommandList = checkCommandArgs "set" commandList VER.setArgList
        leafNameVect = fst3 processedData

    in
    if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'set': " ++ show argList)

    -- this could be changed later
    else if length commandList > 1 || length optionList > 1 then errorWithoutStackTrace ("Set option error: can only have one set argument for each command: " ++ (show (commandList,optionList))) 

    -- early extraction of partition character and bc2-gt64 follows from null inputs
    -- this due to not having all info required for all global settings, so optoin resitrcted and repeated
    else if (null inSeedList) then 
        if head commandList == "partitioncharacter"  then
            let localPartitionChar = head optionList
            in
            if length localPartitionChar /= 1 then errorWithoutStackTrace ("Error in 'set' command. Partitioncharacter '" ++ (show localPartitionChar) ++ "' must be a single character")
            else 
                trace ("PartitionCharacter set to '" ++ (head optionList) ++ "'")
                (globalSettings {partitionCharacter = localPartitionChar}, processedData, inSeedList)

        else if head commandList == "bc2"  then
            let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                noChangeValue = readMaybe noChangeString :: Maybe Double
                changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                changeValue = readMaybe changeString :: Maybe Double
            in
            if length commandList /= length optionList then errorWithoutStackTrace ("Set option error: number of values and options do not match: " ++ (show (commandList,optionList)))
            else if (null . head) optionList then errorWithoutStackTrace ("Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): no values found ")
            else if (',' `notElem` (head optionList)) then errorWithoutStackTrace ("Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): no comma found ")
            else if isNothing noChangeValue then errorWithoutStackTrace ("Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): " ++ (head optionList))
            else if isNothing changeValue then errorWithoutStackTrace ("Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): " ++ (head optionList))
            else 
                if (bc2 globalSettings) /= (fromJust noChangeValue, fromJust changeValue) then 
                    trace ("bit cost 2 state set to " ++ (show (fromJust noChangeValue, fromJust changeValue)))
                    (globalSettings {bc2 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                else (globalSettings {bc2 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
        
        
        else if head commandList == "bc4"  then
            let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                noChangeValue = readMaybe noChangeString :: Maybe Double
                changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                changeValue = readMaybe changeString :: Maybe Double
            in
            if (null . head) optionList then errorWithoutStackTrace ("Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): no values found ")
            else if (',' `notElem` (head optionList)) then errorWithoutStackTrace ("Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): no comma found ")
            else if isNothing noChangeValue then errorWithoutStackTrace ("Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): " ++ (head optionList))
            else if isNothing changeValue then errorWithoutStackTrace ("Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): " ++ (head optionList))
            else 
                if (bc4 globalSettings) /= (fromJust noChangeValue, fromJust changeValue) then 
                    trace ("bit cost 4 state set to " ++ (show (fromJust noChangeValue, fromJust changeValue)))
                    (globalSettings {bc4 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                else (globalSettings {bc4 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
        
        
       else if head commandList == "bc5"  then
            let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                noChangeValue = readMaybe noChangeString :: Maybe Double
                changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                changeValue = readMaybe changeString :: Maybe Double
            in
            if (null . head) optionList then errorWithoutStackTrace ("Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): no values found ")
            else if (',' `notElem` (head optionList)) then errorWithoutStackTrace ("Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): no comma found ")
            else if isNothing noChangeValue then errorWithoutStackTrace ("Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): " ++ (head optionList))
            else if isNothing changeValue then errorWithoutStackTrace ("Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): " ++ (head optionList))
            else 
                if (bc5 globalSettings) /= (fromJust noChangeValue, fromJust changeValue) then 
                    trace ("bit cost 5 state set to " ++ (show (fromJust noChangeValue, fromJust changeValue)))
                    (globalSettings {bc5 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                else (globalSettings {bc5 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
        
        
         else if head commandList == "bc8"  then
            let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                noChangeValue = readMaybe noChangeString :: Maybe Double
                changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                changeValue = readMaybe changeString :: Maybe Double
            in
            if (null . head) optionList then errorWithoutStackTrace ("Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): no values found ")
            else if (',' `notElem` (head optionList)) then errorWithoutStackTrace ("Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): no comma found ")
            else if isNothing noChangeValue then errorWithoutStackTrace ("Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): " ++ (head optionList))
            else if isNothing changeValue then errorWithoutStackTrace ("Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): " ++ (head optionList))
            else 
                if (bc8 globalSettings) /= (fromJust noChangeValue, fromJust changeValue) then 
                    trace ("bit cost 8 state set to " ++ (show (fromJust noChangeValue, fromJust changeValue)))
                    (globalSettings {bc8 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                else (globalSettings {bc8 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
        
        
         else if head commandList == "bc64"  then
            let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                noChangeValue = readMaybe noChangeString :: Maybe Double
                changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                changeValue = readMaybe changeString :: Maybe Double
            in
            if (null . head) optionList then errorWithoutStackTrace ("Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): no values found ")
            else if (',' `notElem` (head optionList)) then errorWithoutStackTrace ("Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): no comma found ")
            else if isNothing noChangeValue then errorWithoutStackTrace ("Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): " ++ (head optionList))
            else if isNothing changeValue then errorWithoutStackTrace ("Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): " ++ (head optionList))
            else 
                if (bc64 globalSettings) /= (fromJust noChangeValue, fromJust changeValue) then 
                    trace ("bit cost 64 state set to " ++ (show (fromJust noChangeValue, fromJust changeValue)))
                    (globalSettings {bc64 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                else (globalSettings {bc64 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
        
        
         else if head commandList == "bcgt64"  then
            let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                noChangeValue = readMaybe noChangeString :: Maybe Double
                changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                changeValue = readMaybe changeString :: Maybe Double
            in
            if (null . head) optionList then errorWithoutStackTrace ("Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): no values found ")
            else if (',' `notElem` (head optionList)) then errorWithoutStackTrace ("Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): no comma found ")
            else if isNothing noChangeValue then errorWithoutStackTrace ("Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): " ++ (head optionList))
            else if isNothing changeValue then errorWithoutStackTrace ("Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): " ++ (head optionList))
            else 
                if (bcgt64 globalSettings) /= (fromJust noChangeValue, fromJust changeValue) then 
                    trace ("bit cost > 64 state set to " ++ (show (fromJust noChangeValue, fromJust changeValue)))
                    (globalSettings {bcgt64 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                else (globalSettings {bcgt64 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
        
        
         -- partition character to reset
        else 
            -- trace ("PartitionCharacter set to '" ++ (partitionCharacter globalSettings) ++ "'")
            (globalSettings, processedData, inSeedList)
    -- regular command stuff not initial at start   
    else
        if head commandList == "bc2"  then
            let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                noChangeValue = readMaybe noChangeString :: Maybe Double
                changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                changeValue = readMaybe changeString :: Maybe Double
            in
            if length commandList /= length optionList then errorWithoutStackTrace ("Set option error: number of values and options do not match: " ++ (show (commandList,optionList)))
            else if (null . head) optionList then errorWithoutStackTrace ("Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): no values found ")
            else if (',' `notElem` (head optionList)) then errorWithoutStackTrace ("Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): no comma found ")
            else if isNothing noChangeValue then errorWithoutStackTrace ("Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): " ++ (head optionList))
            else if isNothing changeValue then errorWithoutStackTrace ("Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): " ++ (head optionList))
            else 
                if (bc2 globalSettings) /= (fromJust noChangeValue, fromJust changeValue) then 
                    trace ("bit cost 2 state set to " ++ (show (fromJust noChangeValue, fromJust changeValue)))
                    (globalSettings {bc2 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                else (globalSettings {bc2 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
        
        
        else if head commandList == "bc4"  then
            let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                noChangeValue = readMaybe noChangeString :: Maybe Double
                changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                changeValue = readMaybe changeString :: Maybe Double
            in
            if (null . head) optionList then errorWithoutStackTrace ("Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): no values found ")
            else if (',' `notElem` (head optionList)) then errorWithoutStackTrace ("Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): no comma found ")
            else if isNothing noChangeValue then errorWithoutStackTrace ("Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): " ++ (head optionList))
            else if isNothing changeValue then errorWithoutStackTrace ("Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): " ++ (head optionList))
            else 
                if (bc4 globalSettings) /= (fromJust noChangeValue, fromJust changeValue) then 
                    trace ("bit cost 4 state set to " ++ (show (fromJust noChangeValue, fromJust changeValue)))
                    (globalSettings {bc4 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                else (globalSettings {bc4 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
        
        
       else if head commandList == "bc5"  then
            let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                noChangeValue = readMaybe noChangeString :: Maybe Double
                changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                changeValue = readMaybe changeString :: Maybe Double
            in
            if (null . head) optionList then errorWithoutStackTrace ("Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): no values found ")
            else if (',' `notElem` (head optionList)) then errorWithoutStackTrace ("Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): no comma found ")
            else if isNothing noChangeValue then errorWithoutStackTrace ("Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): " ++ (head optionList))
            else if isNothing changeValue then errorWithoutStackTrace ("Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): " ++ (head optionList))
            else 
                if (bc5 globalSettings) /= (fromJust noChangeValue, fromJust changeValue) then 
                    trace ("bit cost 5 state set to " ++ (show (fromJust noChangeValue, fromJust changeValue)))
                    (globalSettings {bc5 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                else (globalSettings {bc5 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
        
        
         else if head commandList == "bc8"  then
            let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                noChangeValue = readMaybe noChangeString :: Maybe Double
                changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                changeValue = readMaybe changeString :: Maybe Double
            in
            if (null . head) optionList then errorWithoutStackTrace ("Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): no values found ")
            else if (',' `notElem` (head optionList)) then errorWithoutStackTrace ("Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): no comma found ")
            else if isNothing noChangeValue then errorWithoutStackTrace ("Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): " ++ (head optionList))
            else if isNothing changeValue then errorWithoutStackTrace ("Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): " ++ (head optionList))
            else 
                if (bc8 globalSettings) /= (fromJust noChangeValue, fromJust changeValue) then 
                    trace ("bit cost 8 state set to " ++ (show (fromJust noChangeValue, fromJust changeValue)))
                    (globalSettings {bc8 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                else (globalSettings {bc8 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
        
        
         else if head commandList == "bc64"  then
            let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                noChangeValue = readMaybe noChangeString :: Maybe Double
                changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                changeValue = readMaybe changeString :: Maybe Double
            in
            if (null . head) optionList then errorWithoutStackTrace ("Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): no values found ")
            else if (',' `notElem` (head optionList)) then errorWithoutStackTrace ("Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): no comma found ")
            else if isNothing noChangeValue then errorWithoutStackTrace ("Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): " ++ (head optionList))
            else if isNothing changeValue then errorWithoutStackTrace ("Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): " ++ (head optionList))
            else 
                if (bc64 globalSettings) /= (fromJust noChangeValue, fromJust changeValue) then 
                    trace ("bit cost 64 state set to " ++ (show (fromJust noChangeValue, fromJust changeValue)))
                    (globalSettings {bc64 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                else (globalSettings {bc64 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
        
        
         else if head commandList == "bcgt64"  then
            let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                noChangeValue = readMaybe noChangeString :: Maybe Double
                changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                changeValue = readMaybe changeString :: Maybe Double
            in
            if (null . head) optionList then errorWithoutStackTrace ("Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): no values found ")
            else if (',' `notElem` (head optionList)) then errorWithoutStackTrace ("Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): no comma found ")
            else if isNothing noChangeValue then errorWithoutStackTrace ("Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): " ++ (head optionList))
            else if isNothing changeValue then errorWithoutStackTrace ("Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): " ++ (head optionList))
            else 
                if (bcgt64 globalSettings) /= (fromJust noChangeValue, fromJust changeValue) then 
                    trace ("bit cost > 64 state set to " ++ (show (fromJust noChangeValue, fromJust changeValue)))
                    (globalSettings {bcgt64 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                else (globalSettings {bcgt64 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
        
         --  modified criterion causes changes in graphfactor and root cost 
         else if head commandList == "criterion"  then
            let localCriterion
                  | (head optionList == "parsimony") = Parsimony
                  | (head optionList == "pmdl") = PMDL
                  | (head optionList == "ml") = Likelihood
                  | otherwise = errorWithoutStackTrace ("Error in 'set' command. Criterion '" ++ (head optionList) ++ "' is not 'parsimony', 'ml', or 'pmdl'")

                -- create lazy list of graph complexity indexed by number of network nodes--need leaf number for base tree complexity
                lGraphComplexityList = if localCriterion == Parsimony then IL.repeat (0.0, 0.0)
                                       else if localCriterion `elem` [PMDL, Likelihood] then U.calculateGraphComplexity processedData
                                       else errorWithoutStackTrace ("Optimality criterion not recognized: " ++ (show localCriterion))

                lRootComplexity = if localCriterion == Parsimony then 0.0
                                 else if localCriterion `elem`  [PMDL, Likelihood] then U.calculateW15RootCost processedData
                                 else error ("Optimality criterion not recognized: " ++ (show localCriterion))

                lGraphFactor = if localCriterion `elem` [PMDL, Likelihood] then PMDLGraph
                               else graphFactor globalSettings
            in
            trace ("Optimality criterion set to " ++ (show localCriterion) ++ " Tree Complexity = " ++ (show $ fst $ IL.head lGraphComplexityList) ++ " bits")
            (globalSettings {optimalityCriterion = localCriterion, graphComplexityList = lGraphComplexityList, rootComplexity = lRootComplexity, graphFactor = lGraphFactor}, processedData, inSeedList)

        else if head commandList == "compressresolutions"  then
            let localCriterion
                  | (head optionList == "true") = True
                  | (head optionList == "false") = False
                  | otherwise = errorWithoutStackTrace ("Error in 'set' command. CompressResolutions '" ++ (head optionList) ++ "' is not 'true' or 'false'")
            in
            trace ("CompressResolutions set to " ++ head optionList)
            (globalSettings {compressResolutions = localCriterion}, processedData, inSeedList)

        -- this not intended to be for users
        else if head commandList == "dynamicepsilon"  then
            let localValue = readMaybe (head optionList) :: Maybe Double
            in
            if localValue == Nothing then error ("Set option 'dynamicEpsilon' must be set to an double value >= 0.0 (e.g. dynamicepsilon:0.02): " ++ (head optionList))
            else if (fromJust localValue) < 0.0 then errorWithoutStackTrace ("Set option 'dynamicEpsilon' must be set to an double value >= 0.0 (e.g. dynamicepsilon:0.02): " ++ (head optionList))
            else 
                trace ("Dynamic Epsilon factor set to " ++ head optionList)
                (globalSettings {dynamicEpsilon = 1.0 + ((fromJust localValue) * (fractionDynamic globalSettings))}, processedData, inSeedList)

        else if head commandList == "finalassignment"  then
            let localMethod
                  | ((head optionList == "do") || (head optionList == "directoptimization")) = DirectOptimization
                  | ((head optionList == "ia") || (head optionList == "impliedalignment")) = ImpliedAlignment
                  | otherwise = errorWithoutStackTrace ("Error in 'set' command. FinalAssignment  '" ++ (head optionList) ++ "' is not 'DirectOptimization (DO)' or 'ImpliedAlignment (IA)'")
            in
            if (graphType globalSettings) == Tree then
                trace ("FinalAssignment set to " ++ head optionList)
                (globalSettings {finalAssignment = localMethod}, processedData, inSeedList)
            else if localMethod == DirectOptimization then
                (globalSettings {finalAssignment = localMethod}, processedData, inSeedList)
            else 
                trace ("FinalAssignment set to DO (ignoring IA option) for non-Tree graphs")
                (globalSettings {finalAssignment = DirectOptimization}, processedData, inSeedList)

        else if head commandList == "graphfactor"  then
            let localMethod
                  | (head optionList == "nopenalty") = NoNetworkPenalty
                  | (head optionList == "w15") = Wheeler2015Network
                  | (head optionList == "w23") = Wheeler2023Network
                  | (head optionList == "pmdl") = PMDLGraph
                  | otherwise = errorWithoutStackTrace ("Error in 'set' command. GraphFactor  '" ++ (head optionList) ++ "' is not 'NoPenalty', 'W15', 'W23', or 'PMDL'")
            in
            trace ("GraphFactor set to " ++ (show localMethod))
            (globalSettings {graphFactor = localMethod}, processedData, inSeedList)

        else if head commandList == "graphssteepest"  then
            let localValue = readMaybe (head optionList) :: Maybe Int
            in
            if localValue == Nothing then error ("Set option 'graphsSteepest' must be set to an integer value (e.g. graphsSteepest:5): " ++ (head optionList))
            else 
                trace ("GraphsStreepest set to " ++ head optionList)
                (globalSettings {graphsSteepest = (fromJust localValue)}, processedData, inSeedList)

        else if head commandList == "graphtype"  then
            let localGraphType
                  | (head optionList == "tree") = Tree
                  | (head optionList == "softwired") = SoftWired
                  | (head optionList == "hardwired") = HardWired
                  | otherwise = errorWithoutStackTrace ("Error in 'set' command. Graphtype '" ++ (head optionList) ++ "' is not 'tree', 'hardwired', or 'softwired'")
            in
            if localGraphType /= Tree then 
                trace ("Graphtype set to " ++ (head optionList) ++ " and final assignment to DO")
                (globalSettings {graphType = localGraphType, finalAssignment = DirectOptimization}, processedData, inSeedList)
            else 
                trace ("Graphtype set to " ++ head optionList)
                (globalSettings {graphType = localGraphType}, processedData, inSeedList)

        else if head commandList == "modelcomplexity"  then
            let localValue = readMaybe (head optionList) :: Maybe Double
            in
            if localValue == Nothing then error ("Set option 'modelComplexity' must be set to a double value (e.g. modelComplexity:123.456): " ++ (head optionList))
            else 
                trace ("Model Complexity set to " ++ head optionList)
                (globalSettings {modelComplexity = (fromJust localValue)}, processedData, inSeedList)
        
        else if head commandList == "outgroup"  then
            let outTaxonName = T.pack $ filter (/= '"') $ head $ filter (/= "") $ fmap snd argList
                outTaxonIndex = V.elemIndex outTaxonName leafNameVect

            in
            if isNothing outTaxonIndex then errorWithoutStackTrace ("Error in 'set' command. Out-taxon " ++ T.unpack outTaxonName ++ " not found in input leaf list" ++ show (fmap (T.unpack) leafNameVect))
            else trace ("Outgroup set to " ++ T.unpack outTaxonName) (globalSettings {outgroupIndex = fromJust outTaxonIndex, outGroupName = outTaxonName}, processedData, inSeedList)

        else if head commandList == "partitioncharacter"  then
            let localPartitionChar = head optionList
            in
            if length localPartitionChar /= 1 then errorWithoutStackTrace ("Error in 'set' command. Partitioncharacter '" ++ (show localPartitionChar) ++ "' must be a single character")
            else 
                if (localPartitionChar /= partitionCharacter globalSettings) then
                    trace ("PartitionCharacter set to '" ++ (head optionList) ++ "'")
                    (globalSettings {partitionCharacter = localPartitionChar}, processedData, inSeedList)
                else 
                    (globalSettings, processedData, inSeedList)
        
        else if head commandList == "rootcost"  then
            let localMethod
                  | (head optionList == "norootcost") = NoRootCost
                  | (head optionList == "w15") = Wheeler2015Root
                  | (head optionList == "pmdl") = PMDLRoot
                  | (head optionList == "ml") = MLRoot
                  | otherwise = errorWithoutStackTrace ("Error in 'set' command. RootCost  '" ++ (head optionList) ++ "' is not 'NoRootCost', 'W15', or 'PMDL'")

                lRootComplexity = if localMethod == NoRootCost then 0.0
                                 else if localMethod `elem` [Wheeler2015Root, PMDLRoot, MLRoot] then U.calculateW15RootCost processedData
                                 else error ("Root cost method not recognized: " ++ (show localMethod))
            in
            trace ("RootCost set to " ++ (show localMethod) ++ " " ++ (show lRootComplexity) ++ " bits")
            (globalSettings {rootCost = localMethod, rootComplexity = lRootComplexity}, processedData, inSeedList)

        else if head commandList == "seed"  then
            let localValue = readMaybe (head optionList) :: Maybe Int
            in
            if localValue == Nothing then error ("Set option 'seed' must be set to an integer value (e.g. seed:123): " ++ (head optionList))
            else 
                trace ("Random Seed set to " ++ head optionList)
                (globalSettings {seed = (fromJust localValue)}, processedData, randomIntList (fromJust localValue))

        else trace ("Warning--unrecognized/missing 'set' option in " ++ show argList) (globalSettings, processedData, inSeedList)



-- | reportCommand takes report options, current data and graphs and returns a
-- (potentially large) String to print and the channel to print it to
-- and write mode overwrite/append
reportCommand :: GlobalSettings -> [Argument] -> Int -> String -> ProcessedData -> [PhylogeneticGraph] -> [PhylogeneticGraph] -> [[VertexCost]] -> (String, String, String)
reportCommand globalSettings argList numInputFiles crossReferenceString processedData curGraphs supportGraphs pairwiseDistanceMatrix =
    let argListWithoutReconcileCommands = filter ((`notElem` VER.reconcileArgList) .fst) argList
        --check for balances double quotes and only one pair
        outFileNameList = filter (/= "") $ fmap snd argListWithoutReconcileCommands --argList
        commandList = fmap (fmap C.toLower) $ filter (/= "") $ fmap fst argListWithoutReconcileCommands
        -- reconcileList = filter (/= "") $ fmap fst argList
    in
    if length outFileNameList > 1 then errorWithoutStackTrace ("Report can only have one file name: " ++ (show outFileNameList) ++ " " ++ (show argList))
    else
        let checkCommandList = checkCommandArgs "report" commandList VER.reportArgList
            outfileName = if null outFileNameList then "stderr"
                          else tail $ L.init $ head outFileNameList
            writeMode = if "overwrite" `elem` commandList then "overwrite"
                        else "append"

        in
        -- error too harsh, lose everything else
        --if (null $ filter (/= "overwrite") $ filter (/= "append") commandList) then errorWithoutStackTrace ("Error: Missing 'report' option in " ++ show commandList)
        --else
        if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in report: " ++ show argList)
        else
            -- This for reconciled data
            if "crossrefs" `elem` commandList then
                let dataString = crossReferenceString
                in
                (dataString, outfileName, writeMode)

            else if "data" `elem` commandList then
                let dataString = phyloDataToString 0 $ thd3 processedData
                    baseData = ("There were " ++ (show numInputFiles) ++ " input data files with " ++ show (length $ thd3 processedData) ++ " blocks and " ++ (show ((length dataString) - 1)) ++ " total characters\n")
                    charInfoFields = ["Index", "Block", "Name", "Type", "Activity", "Weight", "Prealigned", "Alphabet", "TCM"]
                in
                (baseData ++ CSV.genCsvFile (charInfoFields : dataString), outfileName, writeMode)

            else if "diagnosis" `elem` commandList then
                let dataString = CSV.genCsvFile $ concatMap (getGraphDiagnosis processedData) (zip curGraphs [0.. ((length curGraphs) - 1)])
                in
                if null curGraphs then 
                    trace ("No graphs to diagnose")
                    ("No graphs to diagnose", outfileName, writeMode)
                else 
                    trace ("Diagnosing " ++ (show $ length curGraphs) ++ " graphs at minimum cost " ++ (show $ minimum $ fmap snd6 curGraphs))
                    (dataString, outfileName, writeMode)

            else if "displaytrees" `elem` commandList then
                -- need to specify -O option for multiple graphs
                let inputDisplayVVList = fmap fth6 curGraphs
                    costList = fmap snd6 curGraphs
                    displayCostListList = fmap GO.getDisplayTreeCostList curGraphs
                    displayInfoString = ("DisplayTree costs : " ++ (show (fmap sum $ fmap fst displayCostListList, displayCostListList))) 
                    treeIndexStringList = fmap ((++ "\n") . ("Canonical Tree " ++)) (fmap show [0..(length inputDisplayVVList - 1)])
                    canonicalGraphPairList = zip treeIndexStringList inputDisplayVVList
                    blockStringList = concatMap (++ "\n") (fmap (outputBlockTrees commandList costList (outgroupIndex globalSettings)) canonicalGraphPairList)
                    -- graphString = outputGraphString commandList (outgroupIndex globalSettings) (fmap thd6 curGraphs) (fmap snd6 curGraphs)
                in
                if null curGraphs || (graphType globalSettings) /= SoftWired then 
                    trace ("No soft-wired graphs to report display trees")
                    ("No soft-wired graphs to report display trees", outfileName, writeMode)
                else 
                    (displayInfoString ++ "\n" ++ blockStringList, outfileName, writeMode)
                

            else if "graphs" `elem` commandList then
            --else if (not .null) (L.intersect ["graphs", "newick", "dot", "dotpdf"] commandList) then
                let 
                    graphString = outputGraphString commandList (outgroupIndex globalSettings) (fmap thd6 curGraphs) (fmap snd6 curGraphs)
                in
                if null curGraphs then 
                    trace ("No graphs to report")
                    ("No graphs to report", outfileName, writeMode)
                else 
                    trace ("Reporting " ++ (show $ length curGraphs) ++ " graph(s) at minimum cost " ++ (show $ minimum $ fmap snd6 curGraphs))
                    (graphString, outfileName, writeMode)

            else if "ia" `elem` commandList || "impliedalignment" `elem` commandList then
                if null curGraphs then 
                    trace ("No graphs to create implied alignments")
                    ("No impliedAlgnments to report", outfileName, writeMode)
                else
                    let includeMissing = any (=="includemissing") commandList
                        concatSeqs = any (=="concatenate") commandList
                        iaContentList = zipWith (getImpliedAlignmentString globalSettings (includeMissing || concatSeqs) concatSeqs processedData) curGraphs [0.. (length curGraphs - 1)]
                    in
                    trace ("\tWarning--prealigned sequence data with non-additive type costs (all change values equal) have been recoded to non-additve characters and will not appear in implied alignment output.")
                    (concat iaContentList, outfileName, writeMode)

            else if "pairdist" `elem` commandList then
                let nameData = L.intercalate "," (V.toList $ fmap T.unpack $ fst3 processedData) ++ "\n"
                    dataString = CSV.genCsvFile $ fmap (fmap show) pairwiseDistanceMatrix
                in
                (nameData ++ dataString, outfileName, writeMode)

            else if "reconcile" `elem` commandList then
                let (reconcileString, _) = R.makeReconcileGraph VER.reconcileArgList argList (fmap fst6 curGraphs)
                in
                if null curGraphs then 
                    trace ("No graphs to reconcile")
                    ([], outfileName, writeMode)
                else 
                    (reconcileString, outfileName, writeMode)

            else if "search" `elem` commandList then
                let dataString' = fmap showSearchFields $ reverse $ searchData globalSettings
                    -- reformat the "search" command fields a bit 
                    dataString = processSearchFields dataString'
                    baseData = ("SearchData\nRandom seed " ++ (show $ seed globalSettings))
                    charInfoFields = ["Command", "Arguments", "Min cost in", "Max cost in", "Num graphs in", "Min cost out", "Max cost out", "Num graphs out", "CPU time (secs)", "Comment"]
                in
                (baseData ++ CSV.genCsvFile (charInfoFields : dataString), outfileName, writeMode)

            else if "support" `elem` commandList then
                let graphString = outputGraphStringSimple commandList (outgroupIndex globalSettings) (fmap fst6 supportGraphs) (fmap snd6 supportGraphs)
                in
                -- trace ("Rep Sup: " ++ (LG.prettify $ fst6 $ head supportGraphs)) (
                if null supportGraphs then 
                    trace ("\tNo support graphs to report")
                    ([], outfileName, writeMode)
                else 
                trace ("Reporting " ++ (show $ length curGraphs) ++ " support graph(s)")
                (graphString, outfileName, writeMode)
                -- )

            else if "tnt" `elem` commandList then 
                if null curGraphs then 
                    trace ("No graphs to create implied alignments for TNT output")
                    ("No impliedAlgnments for TNT to report", outfileName, writeMode)
                else
                    let tntContentList = zipWith (getTNTString globalSettings processedData) curGraphs [0.. (length curGraphs - 1)]
                    in
                    (concat tntContentList, outfileName, writeMode)
           
            else 
                trace ("\nWarning--unrecognized/missing report option in " ++ (show commandList) ++ " defaulting to 'graphs'") (
                let graphString = outputGraphString commandList (outgroupIndex globalSettings) (fmap thd6 curGraphs) (fmap snd6 curGraphs)
                in
                if null curGraphs then 
                    trace ("No graphs to report")
                    ("No graphs to report", outfileName, writeMode)
                else 
                    trace ("Reporting " ++ (show $ length curGraphs) ++ " graph(s) at minimum cost " ++ (show $ minimum $ fmap snd6 curGraphs) ++"\n")
                    (graphString, outfileName, writeMode)
                )

