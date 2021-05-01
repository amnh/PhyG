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

{-
Functions to manage command execution after data have been read and processed

Need top add:
    Report 
        reconciled graphs ala EUN2 codes
    Select() 
        Trees, taxa and blocks or characters to use or not use

    Set()
        Outgroup
        Optimality criterion (Pars/PMDL)
        Blocks  
            create new based on character names
            filter any empty blocks
            Reorganize optimimized data
                Consolidating chracter types, packing etc
    Build
        Trees
        Networks
        Distance trees/Networks
            based on "Block" trees and reconciled via EUN2 code

    Swap 
        NNI, SPR,TBR with distnce from orginal placement
            Edges spiral out from original location
            For networks based on inflectgion points
            Forests etc
    Refine
        Other forms of refinement
        GA (Fusing), Perturbation, Similated Annealing
        Thomsons algorithm and timerd search

    Transformed data 
        Implied alignment
        Search till no improvement

    Support 
        Bremer (vis swap)
        Bootstrap and JackKnife resampling with search options

-}

module CommandExecution ( executeCommands
                        , executeRenameCommands) where

import           Types
import           Debug.Trace
import           GeneralUtilities
import           System.IO
import           GraphFormatUtilities
import qualified LocalGraph as LG
import           Data.List
import qualified Data.Text.Lazy as T
import qualified Data.CSV as CSV
import qualified Data.Vector as V
import qualified Data.Text.Short as ST


-- | executeReadCommands reads iput files and returns raw data 
-- need to close files after read
executeCommands :: [RawData] -> ProcessedData -> [SimpleGraph] -> [[VertexCost]] -> [Command] -> IO [SimpleGraph]
executeCommands rawData processedData curGraphs pairwiseDist commandList = do
    if null commandList then return curGraphs
    else do
        let (firstOption, firstArgs) = head commandList

        -- skip "Read" and "Rename "commands already processed
        if firstOption == Read then error ("Read command should already have been processed: " ++ show (firstOption, firstArgs))
        else if firstOption == Rename then error ("Rename command should already have been processed: " ++ show (firstOption, firstArgs))
        -- report command    
        else if firstOption == Report then do
            let reportStuff@(reportString, outFile, writeMode) = reportCommand firstArgs rawData processedData curGraphs pairwiseDist
            hPutStrLn stderr ("Report writing to " ++ outFile)
            if outFile == "stderr" then hPutStr stderr reportString
            else if outFile == "stdout" then hPutStr stdout reportString
            else if writeMode == "overwrite" then writeFile outFile reportString
            else if writeMode == "append" then appendFile outFile reportString
            else error ("Error 'read' command not properly formatted" ++ (show reportStuff))
            executeCommands rawData processedData curGraphs pairwiseDist (tail commandList)
        else 
            executeCommands rawData processedData curGraphs pairwiseDist (tail commandList)
            
-- | reportArgList contains valid report arguments
reportArgList :: [String]
reportArgList = ["all", "data", "graphs", "overwrite", "append", "dot", "newick", "ascii", "crossrefs", "pairdist"]

-- | checkReportCommands takes commands and verifies that they are in list
checkReportCommands :: [String] -> [String] -> Bool
checkReportCommands commandList permittedList =
    if null commandList then True
    else 
        let firstCommand = head commandList
            foundCommand = firstCommand `elem` permittedList
        in
        if foundCommand then checkReportCommands (tail commandList) permittedList
        else 
            let errorMatch = snd $ getBestMatch (maxBound :: Int ,"no suggestion") permittedList firstCommand
            in
            errorWithoutStackTrace ("\nError: Unrecognized 'report' option. By \'" ++ firstCommand ++ "\' did you mean \'" ++ errorMatch ++ "\'?\n") 

-- | reportCommand takes report options, current data and graphs and returns a 
-- (potentially large) String to print and the channel to print it to 
-- and write mode overwrite/append
reportCommand :: [Argument] -> [RawData] -> ProcessedData -> [SimpleGraph] -> [[VertexCost]] -> (String, String, String)
reportCommand argList rawData processedData curGraphs pairwiseDistanceMatrix =
    let outFileNameList = filter (/= "") $ fmap snd argList
        commandList = filter (/= "") $ fmap fst argList
    in
    if length outFileNameList > 1 then errorWithoutStackTrace ("Report can only have one file name: " ++ show outFileNameList)
    else 
        let checkCommandList = checkReportCommands commandList reportArgList
            outfileName = if null outFileNameList then "stderr" 
                          else tail $ init $ head outFileNameList
            writeMode = if "overwrite" `elem` commandList then "overwrite"
                        else "append"
            
        in
        -- error too harsh, lose everything else
        --if (null $ filter (/= "overwrite") $ filter (/= "append") commandList) then errorWithoutStackTrace ("Error: Missing 'report' option in " ++ show commandList)
        --else 
        if checkCommandList == False then errorWithoutStackTrace ("Unrecognized command in report: " ++ (show argList))
        else 
            -- This for reconciled data
            if "crossrefs" `elem` commandList then
                let dataString = CSV.genCsvFile $ getDataListList rawData processedData
                in 
                (dataString, outfileName, writeMode)
            else if "pairdist" `elem` commandList then
                let nameData = (intercalate "," $ V.toList $ fmap T.unpack $ fst processedData) ++ "\n"
                    dataString = CSV.genCsvFile $ fmap (fmap show) pairwiseDistanceMatrix
                in 
                (nameData ++ dataString, outfileName, writeMode)
            else if "data" `elem` commandList then 
                let dataString = phyloDataToString 0 $ snd processedData
                    baseData = ("There were " ++ (show $ length rawData) ++ " input data files with " ++ (show $ length $ snd processedData) ++ " blocks and " ++ (show $ ((length dataString) - 1)) ++ " total characters\n")
                    charInfoFields = ["Index", "Block", "Name", "Type", "Activity", "Weight", "Prealigned", "Alphabet"]
                in
                (baseData ++ (CSV.genCsvFile $ charInfoFields : dataString), outfileName, writeMode)
            else if "graphs" `elem` commandList then 
                -- need to specify -O option for multiple graphs
                if "dot" `elem` commandList then
                    let graphString = concat $ intersperse "\n" $ fmap fgl2DotString curGraphs
                    in 
                    (graphString, outfileName, writeMode)
                else if "newick" `elem` commandList then
                    let graphString = fglList2ForestEnhancedNewickString curGraphs  True True
                    in 
                    (graphString, outfileName, writeMode)
                else if "ascii" `elem` commandList then
                    let graphString = concat $ fmap LG.fglToPrettyString curGraphs
                    in 
                    (graphString, outfileName, writeMode)
                else -- "dot" as default
                    let graphString = concat $ fmap fgl2DotString curGraphs
                    in 
                    (graphString, outfileName, writeMode)
            else trace ("Warning--unrecognized/missing report option in " ++ show commandList) ("No report specified", outfileName, writeMode)

-- | getDataListList returns a list of lists of Strings for data output as csv
-- for row is source file names, suubsequent rows by taxon with +/- for present absent taxon in 
-- input file
getDataListList :: [RawData] -> ProcessedData -> [[String]]
getDataListList inDataList processedData = 
    if null inDataList then []
    else 
        let fileNames = " " : (fmap (takeWhile (/= ':')) $ fmap T.unpack $ fmap name $ fmap head $ fmap snd inDataList)
            fullTaxList = V.toList $ fst  processedData
            presenceAbsenceList = fmap (isThere inDataList) fullTaxList
            fullMatrix = zipWith (:) (fmap T.unpack fullTaxList) presenceAbsenceList
        in
        --trace (show fileNames)
        fileNames : fullMatrix

-- | isThere takes a list of Rawdata and reurns a String of + - 
isThere :: [RawData] -> T.Text -> [String]
isThere inData inName =
    if null inData then []
    else
        let firstTaxList = fmap fst $ fst $ head inData
        in
        if inName `elem` firstTaxList then "+" : isThere (tail inData) inName
        else  "-" : isThere (tail inData) inName

-- | phyloDataToString converts RawData type to String
-- for additive chars--multiply states by weight is < 1 when outputtting due to conversion on input
phyloDataToString :: Int -> V.Vector BlockData -> [[String]]
phyloDataToString charIndexStart inDataVect = 
    if V.null inDataVect then []
    else 
        let (blockName, _, charInfoVect) = V.head inDataVect
            charStrings = zipWith (:) (replicate (V.length charInfoVect) (T.unpack blockName)) (fmap getCharInfoStrings $ V.toList charInfoVect)
            charNumberString = fmap show [charIndexStart..(charIndexStart + (length charStrings) - 1)] 
            fullMatrix = zipWith (:) charNumberString charStrings
        in 
        fullMatrix ++ phyloDataToString (charIndexStart + (length charStrings)) (V.tail inDataVect)
    
-- | getCharInfoStrings takes charInfo and returns list of Strings of fields
getCharInfoStrings :: CharInfo -> [String]
getCharInfoStrings inChar =
    let activityString = if (activity inChar) == True then "active"
                         else "inactive"
        prealignedString = if (prealigned inChar) == True then "prealigned"
                         else "unaligned"
    in
    [T.unpack $ name inChar, show $ charType inChar, activityString, show $ weight inChar, prealignedString] ++ (fmap ST.toString $ alphabet inChar)

-- | executeRenameCommands takes all the "Rename commands" pairs and 
-- creates a list of pairs of new name and list of old names to be converted
-- as Text
executeRenameCommands :: [(T.Text, T.Text)] -> [Command] -> IO [(T.Text, T.Text)]
executeRenameCommands curPairs commandList  =
    if null commandList then return curPairs
    else do
        let (firstOption, firstArgs) = head commandList

        -- skip "Read" and "Rename "commands already processed
        if firstOption /= Rename then executeRenameCommands curPairs (tail commandList)
        else 
            let newName = T.pack $ snd $ head firstArgs
                newNameList = replicate (length $ tail firstArgs) newName
                newPairs = zip newNameList (fmap T.pack $ fmap snd $ tail firstArgs)
            in
            executeRenameCommands (curPairs ++ newPairs) (tail commandList)

-- | executeSet processes the "set" command
-- set command very general can set outgroup, optimality criterion, blocks
-- executeSet :: 
