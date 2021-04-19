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


-- | executeReadCommands reads iput files and returns raw data 
-- need to close files after read
executeCommands :: [RawData] -> ProcessedData -> [SimpleGraph] -> [Command] -> IO [SimpleGraph]
executeCommands rawData processedData curGraphs commandList = do
    if null commandList then return curGraphs
    else do
        let (firstOption, firstArgs) = head commandList

        -- skip "Read" and "Rename "commands already processed
        if firstOption == Read then error ("Read command should already have been processed: " ++ show (firstOption, firstArgs))
        else if firstOption == Rename then error ("Rename command should already have been processed: " ++ show (firstOption, firstArgs))
        -- report command    
        else if firstOption == Report then do
            let reportStuff@(reportString, outFile, writeMode) = reportCommand firstArgs rawData processedData curGraphs
            hPutStrLn stderr ("Report writing to " ++ outFile)
            if outFile == "stderr" then hPutStr stderr reportString
            else if outFile == "stdout" then hPutStr stdout reportString
            else if writeMode == "overwrite" then writeFile outFile reportString
            else if writeMode == "append" then appendFile outFile reportString
            else error ("Error 'read' command not properly formatted" ++ (show reportStuff))
            executeCommands rawData processedData curGraphs (tail commandList)
        else 
            executeCommands rawData processedData curGraphs (tail commandList)
            
-- | reportArgList contains valid report arguments
reportArgList :: [String]
reportArgList = ["all", "data","graphs", "overwrite", "append", "dot", reverse "newick", "ascii", "crossrefs"]

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
reportCommand :: [Argument] -> [RawData] -> ProcessedData -> [SimpleGraph] -> (String, String, String)
reportCommand argList rawData processedData curGraphs =
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
        if checkCommandList == False then errorWithoutStackTrace ("Unrecognized command in report: " ++ (show argList))
        else 
            -- This for reconciled data
            if "data" `elem` commandList then 
                let baseData = ("There were " ++ (show $ length rawData) ++ " input data files and " ++ (show $ length curGraphs) ++ " graphs\n")
                    dataString = phyloDataToString rawData
                in
                (baseData ++ dataString, outfileName, writeMode)
            else if "graphs" `elem` commandList then 
                -- need to specify -O option for multiple graphs
                if "dot" `elem` commandList then
                    let graphString = concat $ intersperse "\n" $ fmap fgl2DotString curGraphs
                    in 
                    (graphString, outfileName, writeMode)
                else if (reverse "newick") `elem` (take 6 $ fmap reverse commandList) then
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
            else ("Blah", outfileName, writeMode)
            
-- | phyloDataToString converts RawData type to String
-- for additive chars--multiply states by weight is < 1 when outputtting due to conversion on input
phyloDataToString :: [RawData] -> String
phyloDataToString inData = show inData

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
