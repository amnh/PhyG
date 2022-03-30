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
  ) where

import Data.Foldable
import qualified Data.CSV               as CSV
import qualified Data.List              as L
import           Data.Maybe
import           Text.Read
import qualified Data.Text.Lazy         as T
import qualified Data.Text.Short        as ST
import qualified Data.Vector            as V
import qualified Data.Vector.Storable   as SV
import qualified Data.Vector.Unboxed    as UV
import           Data.Time
import           Data.Time.Clock.POSIX
import           Debug.Trace
import           GeneralUtilities
import           GraphFormatUtilities
import qualified Graphs.GraphOperations as GO
import           System.IO
import           Types.Types
import qualified Utilities.LocalGraph   as LG
import qualified Utilities.Utilities    as U
import qualified Data.Char              as C
import qualified Search.Build           as B
import qualified Reconciliation.ReconcileGraphs as R
import qualified GraphOptimization.Traversals as TRA
import qualified Search.Refinement as REF
import qualified Search.Search as S
import qualified Commands.Transform as TRANS
import System.Timing
import Control.DeepSeq
import           Control.Concurrent
import qualified Support.Support as SUP
import           Data.Char
import qualified Data.List.Split as SL 
import Graphs.GraphOperations as GO
import GraphOptimization.Traversals as TRAV
import System.Info
import System.Process
import System.Directory
import qualified Commands.Transform as TRANS
import qualified SymMatrix                   as S
import           Data.Alphabet
import Data.Bits


-- | setArgLIst contains valid 'set' arguments
setArgList :: [String]
setArgList = ["outgroup", "criterion", "graphtype", "compressresolutions", "finalassignment", "graphfactor", "rootcost", "seed"]

-- | reportArgList contains valid 'report' arguments
reportArgList :: [String]
reportArgList = ["all", "data", "search", "graphs", "overwrite", "append", "dot", "dotpdf", "newick", "ascii", "crossrefs", "pairdist", "diagnosis","displaytrees", "reconcile", "support", "ia", "impliedalignment", "tnt"]


-- | buildArgList is the list of valid build arguments
selectArgList :: [String]
selectArgList = ["best", "all", "unique", "atrandom"]

-- | executeCommands reads input files and returns raw data
-- need to close files after read
executeCommands :: GlobalSettings -> [RawData] -> ProcessedData -> ProcessedData -> [PhylogeneticGraph] -> [[VertexCost]] -> [Int] -> [PhylogeneticGraph] -> [Command] -> IO ([PhylogeneticGraph], GlobalSettings, [Int], [PhylogeneticGraph])
executeCommands globalSettings rawData origProcessedData processedData curGraphs pairwiseDist seedList supportGraphList commandList = do
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
            
            executeCommands (globalSettings {searchData = newSearchData}) rawData origProcessedData processedData (curGraphs ++ newGraphList) pairwiseDist (tail seedList) supportGraphList (tail commandList)
        
        else if firstOption == Refine then do 
            (elapsedSeconds, newGraphList) <- timeOp $ pure $ REF.refineGraph firstArgs globalSettings processedData (head seedList) curGraphs
            
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : (searchData globalSettings)   
            
            executeCommands (globalSettings {searchData = newSearchData}) rawData origProcessedData processedData newGraphList pairwiseDist (tail seedList) supportGraphList (tail commandList)
        
        else if firstOption == Fuse then do
            (elapsedSeconds, newGraphList) <- timeOp $ pure $ REF.fuseGraphs firstArgs globalSettings processedData (head seedList) curGraphs
            
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : (searchData globalSettings)   
            
            executeCommands (globalSettings {searchData = newSearchData}) rawData origProcessedData processedData newGraphList pairwiseDist (tail seedList) supportGraphList (tail commandList)
        
        else if firstOption == Report then do
            let reportStuff@(reportString, outFile, writeMode) = reportCommand globalSettings firstArgs rawData  processedData curGraphs supportGraphList pairwiseDist
            let doDotPDF = any (=="dotpdf") $ fmap (fmap toLower . fst) firstArgs

            if null reportString then do
                executeCommands globalSettings rawData origProcessedData processedData curGraphs pairwiseDist seedList supportGraphList (tail commandList)
            else  do
                hPutStrLn stderr ("Report writing to " ++ outFile)

                if doDotPDF then do
                    let reportString' = changeDotPreamble "digraph {" "digraph G {\n\trankdir = LR;\tnode [ shape = rect];\n" reportString
                    printGraphVizDot reportString' outFile
                    executeCommands globalSettings rawData origProcessedData processedData curGraphs pairwiseDist seedList supportGraphList (tail commandList)

                else do
                    if outFile == "stderr" then hPutStr stderr reportString
                    else if outFile == "stdout" then putStr reportString
                    else if writeMode == "overwrite" then writeFile outFile reportString
                    else if writeMode == "append" then appendFile outFile reportString
                    else error ("Error 'read' command not properly formatted" ++ show reportStuff)
                    executeCommands globalSettings rawData origProcessedData processedData curGraphs pairwiseDist seedList supportGraphList (tail commandList)
        
        else if firstOption == Search then do
            (elapsedSeconds, output) <- timeOp $ S.search firstArgs globalSettings processedData pairwiseDist (head seedList) curGraphs
                --in pure result
            -- (newGraphList, serchInfoList) <- S.search firstArgs globalSettings origProcessedData processedData pairwiseDist (head seedList) curGraphs
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs (fst output) (fromIntegral $ toMilliseconds elapsedSeconds) (concat $ fmap (L.intercalate "\n") $ snd output)
            let newSearchData = searchInfo : (searchData globalSettings)
            executeCommands (globalSettings {searchData = newSearchData})  rawData origProcessedData processedData  (fst output) pairwiseDist (tail seedList) supportGraphList (tail commandList)
        
        else if firstOption == Select then do
            (elapsedSeconds, newGraphList) <- timeOp $ pure $ GO.selectPhylogeneticGraph firstArgs (head seedList) selectArgList curGraphs
                
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : (searchData globalSettings)   
            
            executeCommands (globalSettings {searchData = newSearchData}) rawData origProcessedData processedData newGraphList pairwiseDist (tail seedList) supportGraphList (tail commandList)
        
        else if firstOption == Set then 
            -- if set changes graph aspects--may nned to reoptimize
            let (newGlobalSettings, newProcessedData, seedList') = setCommand firstArgs globalSettings processedData seedList
                newGraphList = if not (requireReoptimization globalSettings newGlobalSettings) then curGraphs
                               else trace ("Reoptimizing gaphs") fmap (TRA.multiTraverseFullyLabelGraph newGlobalSettings newProcessedData True True Nothing) (fmap fst6 curGraphs)
                
                searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList 0 "No Comment"
                newSearchData = searchInfo : (searchData newGlobalSettings)   
            in
            
            executeCommands (newGlobalSettings {searchData = newSearchData}) rawData origProcessedData processedData newGraphList pairwiseDist seedList' supportGraphList (tail commandList)
        
        else if firstOption == Swap then do
            (elapsedSeconds, newGraphList) <- timeOp $ pure $ REF.swapMaster firstArgs globalSettings processedData (head seedList)  curGraphs
                
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : (searchData globalSettings)   
            
            executeCommands (globalSettings {searchData = newSearchData}) rawData origProcessedData processedData newGraphList pairwiseDist (tail seedList) supportGraphList (tail commandList)
        
        else if firstOption == Support then do
            (elapsedSeconds, newSupportGraphList) <- timeOp $ pure $ SUP.supportGraph firstArgs globalSettings processedData (head seedList)  curGraphs
                
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newSupportGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : (searchData globalSettings)   
            
            executeCommands (globalSettings {searchData = newSearchData}) rawData origProcessedData processedData curGraphs pairwiseDist (tail seedList) (supportGraphList ++ newSupportGraphList) (tail commandList)

        else if firstOption == Transform then do
            (elapsedSeconds, (newGS, newProcessedData, newGraphs)) <- timeOp $ pure $ TRANS.transform firstArgs globalSettings origProcessedData processedData (head seedList) curGraphs
                
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphs (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : (searchData globalSettings)   
            
            executeCommands (newGS {searchData = newSearchData}) rawData origProcessedData newProcessedData newGraphs pairwiseDist (tail seedList) supportGraphList (tail commandList)

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
setCommand :: [Argument] -> GlobalSettings -> ProcessedData -> [Int] -> (GlobalSettings, ProcessedData, [Int])
setCommand argList globalSettings processedData inSeedList =
    let commandList = fmap (fmap C.toLower) $ filter (/= "") $ fmap fst argList
        optionList = fmap (fmap C.toLower) $ filter (/= "") $ fmap snd argList
        checkCommandList = checkCommandArgs "set" commandList setArgList
        leafNameVect = fst3 processedData

    in
    if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'set': " ++ show argList)
    else
        if head commandList == "outgroup"  then
            let outTaxonName = T.pack $ filter (/= '"') $ head $ filter (/= "") $ fmap snd argList
                outTaxonIndex = V.elemIndex outTaxonName leafNameVect

            in
            if isNothing outTaxonIndex then errorWithoutStackTrace ("Error in 'set' command. Out-taxon " ++ T.unpack outTaxonName ++ " not found in input leaf list" ++ show (fmap (T.unpack) leafNameVect))
            else trace ("Outgroup set to " ++ T.unpack outTaxonName) (globalSettings {outgroupIndex = fromJust outTaxonIndex, outGroupName = outTaxonName}, processedData, inSeedList)
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
        else if head commandList == "criterion"  then
            let localCriterion
                  | (head optionList == "parsimony") = Parsimony
                  | (head optionList == "pmdl") = PMDL
                  | otherwise = errorWithoutStackTrace ("Error in 'set' command. Criterion '" ++ (head optionList) ++ "' is not 'parsimony' or 'pmdl'")
            in
            trace ("Optimality criterion set to " ++ head optionList)
            (globalSettings {optimalityCriterion = localCriterion}, processedData, inSeedList)
        else if head commandList == "compressresolutions"  then
            let localCriterion
                  | (head optionList == "true") = True
                  | (head optionList == "false") = False
                  | otherwise = errorWithoutStackTrace ("Error in 'set' command. CompressResolutions '" ++ (head optionList) ++ "' is not 'true' or 'false'")
            in
            trace ("CompressResolutions set to " ++ head optionList)
            (globalSettings {compressResolutions = localCriterion}, processedData, inSeedList)
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
                  | (head optionList == "pmdl") = PMDLGraph
                  | otherwise = errorWithoutStackTrace ("Error in 'set' command. GraphFactor  '" ++ (head optionList) ++ "' is not 'NoPenalty', 'W15', or 'PMDL'")
            in
            trace ("GraphFactor set to " ++ head optionList)
            (globalSettings {graphFactor = localMethod}, processedData, inSeedList)
        else if head commandList == "rootcost"  then
            let localMethod
                  | (head optionList == "norootcost") = NoRootCost
                  | (head optionList == "w15") = Wheeler2015Root
                  | (head optionList == "pmdl") = PMDLRoot
                  | otherwise = errorWithoutStackTrace ("Error in 'set' command. RootCost  '" ++ (head optionList) ++ "' is not 'NoRootCost', 'W15', or 'PMDL'")
            in
            trace ("RootCost set to " ++ head optionList)
            (globalSettings {rootCost = localMethod}, processedData, inSeedList)
        else if head commandList == "seed"  then
            let localValue = readMaybe (head optionList) :: Maybe Int
            in
            if localValue == Nothing then error ("Set option 'seed' muts be set to an integer value (e.g. seed:123): " ++ (head optionList))
            else 
                trace ("Random Seed set to " ++ head optionList)
                (globalSettings {seed = (fromJust localValue)}, processedData, randomIntList (fromJust localValue))
        else trace ("Warning--unrecognized/missing 'set' option in " ++ show argList) (globalSettings, processedData, inSeedList)



-- | reportCommand takes report options, current data and graphs and returns a
-- (potentially large) String to print and the channel to print it to
-- and write mode overwrite/append
reportCommand :: GlobalSettings -> [Argument] -> [RawData] -> ProcessedData -> [PhylogeneticGraph] -> [PhylogeneticGraph] -> [[VertexCost]] -> (String, String, String)
reportCommand globalSettings argList rawData processedData curGraphs supportGraphs pairwiseDistanceMatrix =
    let argListWithoutReconcileCommands = filter ((`notElem` R.reconcileCommandList) .fst) argList
        outFileNameList = filter (/= "") $ fmap snd argListWithoutReconcileCommands --argList
        commandList = fmap (fmap C.toLower) $ filter (/= "") $ fmap fst argListWithoutReconcileCommands
        -- reconcileList = filter (/= "") $ fmap fst argList
    in
    if length outFileNameList > 1 then errorWithoutStackTrace ("Report can only have one file name: " ++ (show outFileNameList) ++ " " ++ (show argList))
    else
        let checkCommandList = checkCommandArgs "report" commandList reportArgList
            outfileName = if null outFileNameList then "stderr"
                          else tail $ init $ head outFileNameList
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
                let dataString = CSV.genCsvFile $ getDataListList rawData processedData
                in
                (dataString, outfileName, writeMode)

            else if "data" `elem` commandList then
                let dataString = phyloDataToString 0 $ thd3 processedData
                    baseData = ("There were " ++ show (length rawData) ++ " input data files with " ++ show (length $ thd3 processedData) ++ " blocks and " ++ (show ((length dataString) - 1)) ++ " total characters\n")
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
                    treeIndexStringList = fmap ((++ "\n") . ("Canonical Tree " ++)) (fmap show [0..(length inputDisplayVVList - 1)])
                    canonicalGraphPairList = zip treeIndexStringList inputDisplayVVList
                    blockStringList = concatMap (++ "\n") (fmap (outputBlockTrees commandList (outgroupIndex globalSettings)) canonicalGraphPairList)
                    -- graphString = outputGraphString commandList (outgroupIndex globalSettings) (fmap thd6 curGraphs) (fmap snd6 curGraphs)
                in
                if null curGraphs || (graphType globalSettings) /= SoftWired then 
                    trace ("No soft-wired graphs to report display trees")
                    ("No soft-wired graphs to report display trees", outfileName, writeMode)
                else 
                    (blockStringList, outfileName, writeMode)

            else if "graphs" `elem` commandList then
                let graphString = outputGraphString commandList (outgroupIndex globalSettings) (fmap thd6 curGraphs) (fmap snd6 curGraphs)
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
                    let iaContentList = zipWith (getImpliedAlignmentString globalSettings processedData) curGraphs [0.. (length curGraphs - 1)]
                    in
                    (concat iaContentList, outfileName, writeMode)

            else if "pairdist" `elem` commandList then
                let nameData = L.intercalate "," (V.toList $ fmap T.unpack $ fst3 processedData) ++ "\n"
                    dataString = CSV.genCsvFile $ fmap (fmap show) pairwiseDistanceMatrix
                in
                (nameData ++ dataString, outfileName, writeMode)

            else if "reconcile" `elem` commandList then
                let (reconcileString, _) = R.makeReconcileGraph R.reconcileCommandList argList (fmap fst6 curGraphs)
                in
                if null curGraphs then 
                    trace ("No graphs to reconcile")
                    ("No graphs to reconcile", outfileName, writeMode)
                else 
                    (reconcileString, outfileName, writeMode)

            else if "search" `elem` commandList then
                let dataString = fmap showSearchFields $ reverse $ searchData globalSettings
                    baseData = ("SearchData\n")
                    charInfoFields = ["Command", "Arguments", "Min cost in", "Max cost in", "Num graphs in", "Min cost out", "Max cost out", "Num graphs out", "Duration (secs)", "Comment"]
                in
                (baseData ++ CSV.genCsvFile (charInfoFields : dataString), outfileName, writeMode)

            else if "support" `elem` commandList then
                let graphString = outputGraphStringSimple commandList (outgroupIndex globalSettings) (fmap fst6 supportGraphs) (fmap snd6 supportGraphs)
                in
                -- trace ("Rep Sup: " ++ (LG.prettify $ fst6 $ head supportGraphs)) (
                if null supportGraphs then 
                    trace ("No support graphs to report")
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
           
            else trace ("Warning--unrecognized/missing report option in " ++ show commandList) ("No report specified", outfileName, writeMode)

-- changeDotPreamble takes an input string to search for and a new one to add in its place
-- searches through dot file (can have multipl graphs) replacing teh search string each time.
changeDotPreamble :: String -> String -> String -> String
changeDotPreamble findString newString inDotString = 
    if null inDotString then []
    else 
        changePreamble' findString newString [] (lines inDotString)

-- changeDotPreamble' internal process for changeDotPreamble
changePreamble' :: String -> String -> [String] -> [String] -> String
changePreamble' findString newString accumList inLineList =
    if null inLineList then unlines $ reverse accumList
    else 
        -- trace ("CP':" ++ (head inLineList) ++ " " ++  findString ++ " " ++ newString) (
        let firstLine = head inLineList
        in
        if firstLine == findString then changePreamble' findString newString (newString : accumList) (tail inLineList)
        else changePreamble' findString newString (firstLine : accumList) (tail inLineList)
        -- )

--printGraph graphviz simple dot file of graph
--execute with "dot -Tps test.dot -o test.ps"
--need to add output to argument filename and call 
--graphviz via System.Process.runprocess
--also, reorder GenForest so smalles (num leaves) is either first or
--last so can print small to large all the way so easier to read
printGraphVizDot :: String -> String -> IO ()
printGraphVizDot graphDotString dotFile =
    if null graphDotString then error "No graph to report"
    else do
        myHandle <- openFile dotFile WriteMode
        if os /= "darwin" then hPutStrLn  stderr ("\tOutputting graphviz to " ++ dotFile ++ ".pdf.")
        else hPutStrLn  stderr ("\tOutputting graphviz to " ++ dotFile ++ ".eps.")
        let outputType = if os == "darwin" then "-Teps"
                         else "-Tpdf"
        --hPutStrLn myHandle "digraph G {"
        --hPutStrLn myHandle "\trankdir = LR;"
        --hPutStrLn myHandle "\tnode [ shape = rect];"
        --hPutStr myHandle $ (unlines . tail . lines) graphDotString
        hPutStr myHandle graphDotString
        -- hPutStrLn myHandle "}"
        hClose myHandle
        pCode <- findExecutable "dot" --system "dot" --check for Graphviz
        {-
        hPutStrLn stderr
            (if isJust pCode then --pCode /= Nothing then
                "executed dot " ++ outputType ++ dotFile ++ " -O " else
                "Graphviz call failed (not installed or found).  Dot file still created. Dot can be obtained from https://graphviz.org/download")
        -}
        if isJust pCode then do
            cpResult  <- createProcess (proc "dot" [outputType, dotFile, "-O"])
            hPutStrLn stderr ("\tExecuted dot " ++ outputType ++ " " ++ dotFile ++ " -O ") 
        else 
            hPutStrLn stderr "\tGraphviz call failed (not installed or found).  Dot file still created. Dot can be obtained from https://graphviz.org/download"


-- | showSearchFields cretes a String list for SearchData Fields
showSearchFields :: SearchData -> [String]
showSearchFields sD =
    [show $ instruction sD, concat $ fmap showArg $ arguments sD, show $ minGraphCostIn sD, show $ maxGraphCostIn sD, show $ numGraphsIn sD, show $ minGraphCostOut sD, show $ maxGraphCostOut sD, show $ numGraphsOut sD, 
    show $ (fromIntegral $ duration sD) / 1000, commentString sD]
    where showArg a = "(" ++ (fst a) ++ "," ++ (snd a) ++ ")"

-- | requireReoptimization checks if set command in globalk settings requires reoptimization of graphs due to change in
-- graph type, optimality criterion etc.
requireReoptimization :: GlobalSettings -> GlobalSettings -> Bool 
requireReoptimization gsOld gsNew =
    if graphType gsOld /= graphType gsNew then True
    else if optimalityCriterion gsOld /= optimalityCriterion gsNew then True
    else if finalAssignment gsOld /= finalAssignment gsNew then True
    else if graphFactor gsOld /= graphFactor gsNew then True
    else if rootCost gsOld /= rootCost gsNew then True
    else False

-- | outputBlockTrees takes a PhyloGeneticTree and outputs BlockTrees
outputBlockTrees :: [String] -> Int -> (String , V.Vector [BlockDisplayForest]) -> String
outputBlockTrees commandList lOutgroupIndex (labelString, graphLV) =
    let blockIndexStringList = fmap ((++ "\n") . ("Block " ++)) (fmap show [0..((V.length graphLV) - 1)])
        blockStrings = concatMap (++ "\n") (fmap (makeBlockGraphStrings commandList lOutgroupIndex ) $ zip blockIndexStringList (V.toList graphLV))
    in
    labelString ++ blockStrings

-- | makeBlockGraphStrings makes individual block display trees--potentially multiple
makeBlockGraphStrings :: [String] -> Int -> (String ,[BlockDisplayForest]) -> String
makeBlockGraphStrings commandList lOutgroupIndex (labelString, graphL) =
    let diplayIndexString =("Display Tree(s): " ++ show (length graphL) ++ "\n")
        displayString = (++ "\n") $ outputDisplayString commandList lOutgroupIndex graphL
    in
    labelString ++ diplayIndexString ++ displayString

-- | outputDisplayString is a wrapper around graph output functions--but without cost list
outputDisplayString :: [String] -> Int -> [DecoratedGraph] -> String
outputDisplayString commandList lOutgroupIndex graphList
  | "dot" `elem` commandList = makeDotList lOutgroupIndex (fmap GO.convertDecoratedToSimpleGraph graphList)
  | "newick" `elem` commandList = makeNewickList lOutgroupIndex (fmap GO.convertDecoratedToSimpleGraph graphList) (replicate (length graphList) 0.0)
  | "ascii" `elem` commandList = makeAsciiList lOutgroupIndex (fmap GO.convertDecoratedToSimpleGraph graphList)
  | otherwise = -- "dot" as default
    makeDotList lOutgroupIndex (fmap GO.convertDecoratedToSimpleGraph graphList)

-- | outputGraphString is a wrapper arounf graph output functions
outputGraphString :: [String] -> Int -> [DecoratedGraph] ->  [VertexCost] -> String
outputGraphString commandList lOutgroupIndex graphList costList
  | "dot" `elem` commandList = makeDotList lOutgroupIndex (fmap GO.convertDecoratedToSimpleGraph graphList)
  | "newick" `elem` commandList = makeNewickList lOutgroupIndex (fmap GO.convertDecoratedToSimpleGraph graphList) costList
  | "ascii" `elem` commandList = makeAsciiList lOutgroupIndex (fmap GO.convertDecoratedToSimpleGraph graphList)
  | otherwise = -- "dot" as default
    makeDotList lOutgroupIndex (fmap GO.convertDecoratedToSimpleGraph graphList)

-- | outputGraphStringSimple is a wrapper arounf graph output functions
outputGraphStringSimple :: [String] -> Int -> [SimpleGraph] ->  [VertexCost] -> String
outputGraphStringSimple commandList lOutgroupIndex graphList costList
  | "dot" `elem` commandList = makeDotList lOutgroupIndex graphList
  | "newick" `elem` commandList = makeNewickList lOutgroupIndex graphList costList
  | "ascii" `elem` commandList = makeAsciiList lOutgroupIndex graphList
  | otherwise = -- "dot" as default
    makeDotList lOutgroupIndex graphList


-- | makeDotList takes a list of fgl trees and outputs a single String cointaining the graphs in Dot format
-- need to specify -O option for multiple graph(outgroupIndex globalSettings)s
makeDotList :: Int -> [SimpleGraph] -> String
makeDotList rootIndex graphList =
     L.intercalate "\n" (fmap fgl2DotString $ fmap (GO.rerootTree rootIndex) graphList)

-- | makeNewickList takes a list of fgl trees and outputs a single String cointaining the graphs in Newick format
makeNewickList ::  Int -> [SimpleGraph] -> [VertexCost] -> String
makeNewickList rootIndex graphList costList =
    let graphString = fglList2ForestEnhancedNewickString (fmap (GO.rerootTree rootIndex) graphList)  True True
        newickStringList = fmap init $ filter (not . null) $ lines graphString
        costStringList  = fmap (('[' :) . (++ "];\n")) (fmap show costList)
        graphStringCost = concat $ zipWith (++) newickStringList costStringList
    in
    graphStringCost

-- | makeAsciiList takes a list of fgl trees and outputs a single String cointaining the graphs in ascii format
makeAsciiList :: Int -> [SimpleGraph] -> String
makeAsciiList rootIndex graphList =
    concatMap LG.prettify (fmap (GO.rerootTree rootIndex) graphList)

-- | getDataListList returns a list of lists of Strings for data output as csv
-- for row is source file names, suubsequent rows by taxon with +/- for present absent taxon in
-- input file
getDataListList :: [RawData] -> ProcessedData -> [[String]]
getDataListList inDataList processedData =
    if null inDataList then []
    else
        let fileNames = " " : fmap (takeWhile (/= ':')) (fmap T.unpack $ fmap name $ fmap head $ fmap snd inDataList)
            fullTaxList = V.toList $ fst3  processedData
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
            charStrings = zipWith (:) (replicate (V.length charInfoVect) (T.unpack blockName)) (getCharInfoStrings <$> V.toList charInfoVect)
            charNumberString = fmap show [charIndexStart..(charIndexStart + length charStrings - 1)]
            fullMatrix = zipWith (:) charNumberString charStrings
        in
        fullMatrix ++ phyloDataToString (charIndexStart + length charStrings) (V.tail inDataVect)

-- | getCharInfoStrings takes charInfo and returns list of Strings of fields
getCharInfoStrings :: CharInfo -> [String]
getCharInfoStrings inChar =
    let activityString = if (activity inChar) then "active"
                         else "inactive"
        prealignedString = if (prealigned inChar) then "prealigned"
                         else "unaligned"
    in
    [T.unpack $ name inChar, show $ charType inChar, activityString, show $ weight inChar, prealignedString] <> (fmap ST.toString . toList $ alphabet inChar) <> [show $ costMatrix inChar]

-- | executeRenameReblockCommands takes all the "Rename commands" pairs and
-- creates a list of pairs of new name and list of old names to be converted
-- as Text
executeRenameReblockCommands :: [(T.Text, T.Text)] -> [Command] -> IO [(T.Text, T.Text)]
executeRenameReblockCommands curPairs commandList  =
    if null commandList then return curPairs
    else do
        let (firstOption, firstArgs) = head commandList

        -- skip "Read" and "Rename "commands already processed
        if (firstOption /= Rename) && (firstOption /= Reblock) then executeRenameReblockCommands curPairs (tail commandList)
        else
            let newName = T.filter C.isPrint $ T.filter (/= '"') $ T.pack $ snd $ head firstArgs
                newNameList = replicate (length $ tail firstArgs) newName
                oldNameList = (fmap (T.filter (/= '"') . T.pack) (fmap snd $ tail firstArgs))
                newPairs = zip newNameList oldNameList
            in
            executeRenameReblockCommands (curPairs ++ newPairs) (tail commandList)

-- | getGraphDiagnosis creates basic for CSV of graph vertex and node information
-- nodes first then vertices
getGraphDiagnosis :: ProcessedData -> (PhylogeneticGraph, Int) -> [[String]]
getGraphDiagnosis inData (inGraph, graphIndex) =
    let decGraph = thd6 inGraph
    in
    if LG.isEmpty decGraph then []
    else
        let vertexList = LG.labNodes decGraph
            edgeList = LG.labEdges decGraph
            topHeaderList  = ["Graph Index", "Vertex Index", "Vertex Name", "Vertex Type", "Child Vertices", "Parent Vertices", "Data Block", "Character Name", "Character Type", "Preliminary State", "Final State", "Local Cost"]
            vertexInfoList =  concatMap (getVertexCharInfo (thd3 inData) (fst6 inGraph) (six6 inGraph)) vertexList
            edgeHeaderList = [[" "],[" ", "Edge Head Vertex", "Edge Tail Vertex", "Edge Type", "Minimum Length", "Maximum Length", "MidRange Length"]]
            edgeInfoList = fmap getEdgeInfo edgeList
        in
        [topHeaderList, [show graphIndex]] ++ vertexInfoList ++ edgeHeaderList ++ edgeInfoList

-- | getVertexCharInfo returns a list of list of Strings of vertex infomation
-- one list for each character at the vertex
getVertexCharInfo :: V.Vector BlockData -> SimpleGraph -> V.Vector (V.Vector CharInfo) -> LG.LNode VertexInfo -> [[String]]
getVertexCharInfo blockDataVect inGraph charInfoVectVect inVert =
    let leafParents = LG.parents inGraph (fst inVert)
        parentNodes
          | nodeType  (snd inVert) == RootNode = "None"
          | nodeType  (snd inVert) == LeafNode = show leafParents
          | otherwise = show $  parents  (snd inVert)
        childNodes = if nodeType  (snd inVert) == LeafNode then "None" else show $  children  (snd inVert)
        basicInfoList = [" ", show $ fst inVert, T.unpack $ vertName (snd inVert), show $ nodeType  (snd inVert), childNodes, parentNodes, " ", " ", " ", " ", " ", show $ vertexCost (snd inVert)]
        blockCharVect = V.zip3  (V.map fst3 blockDataVect)  (vertData  (snd inVert)) charInfoVectVect
        blockInfoList = concat $ V.toList $ V.map getBlockList blockCharVect
    in
    basicInfoList : blockInfoList

-- | getBlockList takes a pair of Vector of chardata and vector of charInfo and returns Strings
getBlockList :: (NameText, V.Vector CharacterData, V.Vector CharInfo) -> [[String]]
getBlockList (blockName, blockDataVect, charInfoVect) =
    let firstLine = [" ", " ", " ", " ", " ", " ", T.unpack blockName]
        charlines = V.toList $ V.map makeCharLine (V.zip blockDataVect charInfoVect)
    in
    firstLine : charlines

-- | makeCharLine takes character data
-- will be less legible for optimized data--so should use a diagnosis
-- based on "naive" data for human legible output
-- need to add back-converting to observed states using alphabet in charInfo
makeCharLine :: (CharacterData, CharInfo) -> [String]
makeCharLine (blockDatum, charInfo) =
    let localType = charType charInfo
        localAlphabet = fmap ST.toString $ alphabet charInfo
        isPrealigned = if prealigned charInfo == True then "Prealigned "
                       else ""
        enhancedCharType = if localType `elem` sequenceCharacterTypes then (isPrealigned ++ (show localType))
                        else if localType `elem`  [Add, NonAdd, Matrix] then (show localType)
                        else error ("Character Type :" ++ (show localType) ++ "unrecogniized or not implemented")

        (stringPrelim, stringFinal) = if localType == Add then (show $ snd3 $ rangePrelim blockDatum, show $ rangeFinal blockDatum)
                                      else if localType == NonAdd then (concat $ V.map (U.bitVectToCharState localAlphabet) $ snd3 $ stateBVPrelim blockDatum, concat $ V.map (U.bitVectToCharState localAlphabet) $ stateBVFinal blockDatum)
                                      else if localType == Matrix then (show $ matrixStatesPrelim blockDatum, show $ fmap (fmap fst3) $ matrixStatesFinal blockDatum)
                                      else if localType `elem` sequenceCharacterTypes
                                      then case localType of
                                             x | x `elem` [SlimSeq, NucSeq  ] -> (SV.foldMap (U.bitVectToCharState localAlphabet) $ slimPrelim blockDatum, SV.foldMap (U.bitVectToCharState localAlphabet) $ slimFinal blockDatum)
                                             x | x `elem` [WideSeq, AminoSeq] -> (UV.foldMap (U.bitVectToCharState localAlphabet) $ widePrelim blockDatum, UV.foldMap (U.bitVectToCharState localAlphabet) $ wideFinal blockDatum)
                                             x | x `elem` [HugeSeq]           -> (   foldMap (U.bitVectToCharState localAlphabet) $ hugePrelim blockDatum,    foldMap (U.bitVectToCharState localAlphabet) $ hugeFinal blockDatum)
                                             x | x `elem` [AlignedSlim]       -> (SV.foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ alignedSlimPrelim blockDatum, SV.foldMap (U.bitVectToCharState localAlphabet) $ alignedSlimFinal blockDatum)
                                             x | x `elem` [AlignedWide]       -> (UV.foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ alignedWidePrelim blockDatum, UV.foldMap (U.bitVectToCharState localAlphabet) $ alignedWideFinal blockDatum)
                                             x | x `elem` [AlignedHuge]       -> (   foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ alignedHugePrelim blockDatum,    foldMap (U.bitVectToCharState localAlphabet) $ alignedHugeFinal blockDatum)
                                             
                                             _                                -> error ("Un-implemented data type " ++ show localType)
                                      else error ("Un-implemented data type " ++ show localType)
        in
        [" ", " ", " ", " ", " ", " ", " ", T.unpack $ name charInfo, enhancedCharType, stringPrelim, stringFinal, show $ localCost blockDatum]


-- | getEdgeInfo returns a list of Strings of edge infomation
getEdgeInfo :: LG.LEdge EdgeInfo -> [String]
getEdgeInfo inEdge =
    [" ", show $ fst3 inEdge, show $ snd3 inEdge, show $ edgeType (thd3 inEdge), show $ minLength (thd3 inEdge), show $ maxLength (thd3 inEdge), show $ midRangeLength (thd3 inEdge)]

-- | TNT report functions

-- | getTNTStrings returns as a single String the implied alignments of all sequence characters
-- softwired use display trees, hardWired transform to softwired then proceed with display trees
-- key to keep cost matrices and weights
getTNTString :: GlobalSettings -> ProcessedData -> PhylogeneticGraph -> Int -> String
getTNTString inGS inData inGraph graphNumber = 
    if LG.isEmpty (fst6 inGraph) then error "No graphs for create TNT data for in getTNTString"
    else
        let leafList = snd4 $ LG.splitVertexList (thd6 inGraph)
            leafNameList = fmap (++ "\t") $ fmap T.unpack $ fmap (vertName . snd) leafList
            headerString = "xread\n'TNT data for Graph " ++ (show graphNumber) ++ " generated by PhylogeneticGraph (PhyG)\n\tSource characters:\n"
            finalString = "proc/;\n\n"
            numTaxa = V.length $ fst3 inData
            charInfoVV = six6 inGraph

            -- get character information in 3-tuples and list of lengths--to match lengths
            ccCodeInfo = getCharacterInfo charInfoVV

        in
            
        if graphType inGS == Tree then 
            let leafDataList = V.fromList $ fmap (vertData . snd) leafList

                -- get character strings
                taxonCharacterStringList = V.toList $ fmap (++ "\n") $ fmap (getTaxonCharString charInfoVV) leafDataList 
                nameCharStringList = concat $ zipWith (++) leafNameList taxonCharacterStringList

                -- length information for cc code extents
                charLengthList = concat $ V.toList $ V.zipWith getBlockLength (V.head leafDataList) charInfoVV

                -- Block/Character names for use in comment to show sources of new characters
                charNameList = concat $ V.toList $ fmap getBlockNames charInfoVV

                nameLengthPairList = zip charNameList charLengthList
                nameLengthString = concat $ pairListToStringList nameLengthPairList 0

                -- merge lengths and cc codes
                ccCodeString = mergeCharInfoCharLength ccCodeInfo charLengthList 0
            in
            -- trace ("GTNTS:" ++ (show charLengthList)) 
            headerString  ++ nameLengthString ++ "'\n" ++ (show $ sum charLengthList) ++ " " ++ (show numTaxa) ++ "\n" 
                ++ nameCharStringList ++ ";\n" ++ ccCodeString ++ finalString
            

        -- for softwired networks--use display trees
        else if graphType inGS == SoftWired then 

            -- get display trees for each data block-- takes first of potentially multiple
            let blockDisplayList = fmap GO.convertDecoratedToSimpleGraph $ fmap head $ fth6 inGraph

                -- create seprate processed data for each block
                blockProcessedDataList = fmap (makeBlockData (fst3 inData) (snd3 inData)) (thd3 inData)

                -- Perform full optimizations on display trees (as trees) with single block data (blockProcessedDataList) to creeate IAs
                decoratedBlockTreeList = V.zipWith (TRAV.multiTraverseFullyLabelGraph' (inGS {graphType = Tree}) False False Nothing) blockProcessedDataList blockDisplayList

                -- create leaf data by merging display graph block data (each one a phylogentic graph)
                (leafDataList, mergedCharInfoVV) = mergeDataBlocks (V.toList decoratedBlockTreeList) [] []

                -- get character strings
                taxonCharacterStringList = V.toList $ fmap (++ "\n") $ fmap (getTaxonCharString mergedCharInfoVV) leafDataList 
                nameCharStringList = concat $ zipWith (++) leafNameList taxonCharacterStringList

                -- length information for cc code extents
                charLengthList = concat $ V.toList $ V.zipWith getBlockLength (V.head leafDataList) mergedCharInfoVV

                -- Block/Character names for use in comment to show sources of new characters
                charNameList = concat $ V.toList $ fmap getBlockNames charInfoVV
                
                nameLengthPairList = zip charNameList charLengthList
                nameLengthString = concat $ pairListToStringList nameLengthPairList 0

                -- merge lengths and cc codes
                ccCodeString = mergeCharInfoCharLength ccCodeInfo charLengthList 0
            in
            headerString ++ nameLengthString ++ "'\n" ++ (show $ sum charLengthList) ++ " " ++ (show numTaxa) ++ "\n" 
                ++ nameCharStringList ++ ";\n" ++ ccCodeString ++ finalString     

        else 
            trace ("TNT  not yet implemented for graphtype " ++ show (graphType inGS))
            "There is no implied alignment for hard-wired graphs--at least not yet.\n\tCould transform graph to softwired and generate TNT text that way"

-- | pairListToStringList takes  alist of (String, Int) and a starting index and returns scope of charcter for leading comment
pairListToStringList :: [(String, Int)] -> Int -> [String]
pairListToStringList pairList startIndex =
    if null pairList then []
    else
        let (a, b) = head pairList
        in
        ("\t\t" ++ (show startIndex) ++ "-" ++ (show $ b + startIndex - 1) ++ " : " ++ a ++ "\n") : pairListToStringList (tail pairList) (startIndex + b)

-- | mergeDataBlocks takes a list of Phylogenetic Graphs (Trees) and merges the data blocks (each graph should have only 1)
-- and merges the charInfo Vectors returning data and charInfo
mergeDataBlocks :: [PhylogeneticGraph] -> [[(V.Vector CharacterData)]] -> [V.Vector CharInfo] -> (V.Vector (V.Vector (V.Vector CharacterData)), V.Vector (V.Vector CharInfo))
mergeDataBlocks inGraphList curDataList curInfoList = 
    if null inGraphList then (V.fromList $ fmap V.fromList $ fmap reverse curDataList, V.fromList $ reverse curInfoList)
    else 
        let firstGraph = head inGraphList
            firstTree = thd6 firstGraph
            firstCharInfo = V.head $ six6 firstGraph
            leafList = snd4 $ LG.splitVertexList firstTree

            -- since each graph has a single block--take head to get vector of characters
            leafCharacterList = V.toList $ fmap V.head $ fmap (vertData . snd) (V.fromList leafList)

            -- zip data for each taxon
            newDataList = if null curDataList then fmap (:[]) $ leafCharacterList
                          else zipWith (:) leafCharacterList curDataList
        in
        mergeDataBlocks (tail inGraphList) newDataList (firstCharInfo : curInfoList)

-- | getTaxonCharString returns the total character string for a taxon
-- length and zipping for missing data
getTaxonCharString ::  V.Vector (V.Vector CharInfo) -> VertexBlockData -> String
getTaxonCharString charInfoVV charDataVV =
    let lengthBlock = maximum $ V.zipWith U.getCharacterLength (V.head charDataVV) (V.head charInfoVV)
    in
    concat $ V.zipWith (getBlockString lengthBlock) charInfoVV charDataVV

-- | getBlockString returns the String for a character block
-- returns all '?' if missing
getBlockString :: Int -> V.Vector CharInfo -> V.Vector CharacterData -> String
getBlockString lengthBlock charInfoV charDataV =
    -- this to deal with missing characters
    -- trace ("GBS: " ++ (show $ V.length charDataV)) (
    if V.null charDataV then L.replicate lengthBlock '?' 
    else concat $ V.zipWith getCharacterString  charDataV charInfoV
    -- )



-- | mergeCharInfoCharLength merges cc code char info and char lengths for scope
mergeCharInfoCharLength :: [(String, String, String)] -> [Int] -> Int -> String
mergeCharInfoCharLength codeList lengthList charIndex =      
    if null codeList then []
    else 
        let (ccCodeString, costsString, weightString) = head codeList
            charLength = head lengthList
            startScope = show charIndex
            endScope = show (charIndex + charLength - 1)
            scope = startScope ++ "." ++ endScope
            weightString' = if null weightString then []
                            else "cc /" ++ weightString ++ scope ++ ";\n"
            costsString' = if null costsString then []
                           else "costs " ++ scope ++ " = " ++ costsString ++ ";\n"
            ccCodeString' = "cc " ++ ccCodeString ++ " " ++ scope ++ ";\n"
        in
        (ccCodeString' ++ weightString' ++ costsString') ++ (mergeCharInfoCharLength (tail codeList) (tail lengthList) (charIndex + charLength))
        


-- | getCharacterInfo takes charInfo vect vect and reiurns triples of ccCode, costs, and weight values 
-- for each character
getCharacterInfo :: V.Vector (V.Vector CharInfo) -> [(String, String, String)]
getCharacterInfo inCharInfoVV = 
   concat $ V.toList $ fmap getBlockInfo inCharInfoVV

-- | getBlockInfo gets character code info for a block
getBlockInfo :: V.Vector CharInfo -> [(String, String, String)]
getBlockInfo inCharInfoV = V.toList $ fmap getCharCodeInfo inCharInfoV

-- | getCharCodeInfo extracts 3-tuple of cc code, costs and weight as strings
-- from charInfo
getCharCodeInfo :: CharInfo -> (String, String, String)
getCharCodeInfo inCharInfo =
    let inCharType = charType inCharInfo
        charWeightString = if weight inCharInfo ==  1 then ""
                           else show $ weight inCharInfo
        inAlph  = alphabet inCharInfo
        inMatrix = costMatrix inCharInfo
        (costMatrixType, _) = TRANS.getRecodingType inMatrix
        matrixString = if costMatrixType == "nonAdd" then ""
                       else makeMatrixString inAlph inMatrix
    in
    let codeTriple = case inCharType of
                      x | x `elem` [Add              ] -> ("+", "", charWeightString)
                      x | x `elem` [NonAdd           ] -> ("-", "", charWeightString)
                      x | x `elem` [Matrix           ] -> ("(", matrixString, charWeightString)
                      x | x `elem` sequenceCharacterTypes -> if costMatrixType == "nonAdd" then ("-", "", charWeightString) 
                                                             else ("(", matrixString, charWeightString)
                      _                                -> error ("Un-implemented data type " ++ show inCharType)
    in
    codeTriple

-- | makeMatrixString takes alphabet and input cost matrix and creates TNT 
-- matrix cost line
-- could be lesser but might not be symmetrical
makeMatrixString :: Alphabet ST.ShortText -> S.Matrix Int -> String
makeMatrixString inAlphabet inMatrix = 
    let elementList =  fmap ST.toString $ toList inAlphabet

        -- get element pairs (Strings)
        elementPairList = filter notDiag $ getListPairs elementList

        -- index pairs for accessing matrix
        elementIndexPairList = filter notDiag $ getListPairs [0 .. (length elementList - 1)]
        elementPairCosts = fmap (inMatrix S.!) elementIndexPairList

        -- make strings of form state_i / state_j cost ...
        costString = makeCostString elementPairList elementPairCosts 

    in
    costString
    where notDiag (a,b) = if a < b then True else False

-- | makeCostString takes list of state pairs and list of costs and creates tnt cost string
makeCostString :: [(String, String)] -> [Int] -> String
makeCostString namePairList costList =
    if null namePairList then []
    else 
        let (a,b) = head namePairList
            c = head costList
        in
        (a ++ "/" ++ b ++ " " ++ (show c) ++ " ") ++ (makeCostString (tail namePairList) (tail costList))

-- | getBlockLength returns a list of the lengths of all characters in a blocks
getBlockLength :: V.Vector CharacterData -> V.Vector CharInfo -> [Int]
getBlockLength inCharDataV inCharInfoV =
    -- trace ("GBL:" ++ (show $ V.zipWith U.getCharacterLength inCharDataV inCharInfoV))
    V.toList $ V.zipWith U.getCharacterLength inCharDataV inCharInfoV

-- | getBlockNames returns a list of the lengths of all characters in a blocks
getBlockNames :: V.Vector CharInfo -> [String]
getBlockNames inCharInfoV =
    -- trace ("GBL:" ++ (show $ V.zipWith U.getCharacterLength inCharDataV inCharInfoV))
    V.toList $ fmap T.unpack $ fmap name inCharInfoV


-- | getCharacterString returns a string of character states
-- need to add splace between (large alphabets etc)
-- local alphabet for charactes where that is input.  MAytrix and additive are integers
getCharacterString :: CharacterData -> CharInfo -> String
getCharacterString inCharData inCharInfo = 
    let inCharType = charType inCharInfo
        localAlphabet = if inCharType /= NonAdd then fmap ST.toString $ alphabet inCharInfo
                        else fmap ST.toString discreteAlphabet
    in
    let charString = case inCharType of
                      x | x `elem` [NonAdd           ] ->    foldMap (bitVectToCharStringTNT localAlphabet) $ snd3 $ stateBVPrelim inCharData
                      x | x `elem` [Add              ] ->    foldMap  U.additivStateToString $ snd3 $ rangePrelim inCharData
                      x | x `elem` [Matrix           ] ->    foldMap  U.matrixStateToString  $ matrixStatesPrelim inCharData
                      x | x `elem` [SlimSeq, NucSeq  ] -> SV.foldMap (bitVectToCharStringTNT localAlphabet) $ snd3 $ slimAlignment inCharData
                      x | x `elem` [WideSeq, AminoSeq] -> UV.foldMap (bitVectToCharStringTNT localAlphabet) $ snd3 $ wideAlignment inCharData
                      x | x `elem` [HugeSeq]           ->    foldMap (bitVectToCharStringTNT localAlphabet) $ snd3 $ hugeAlignment inCharData
                      x | x `elem` [AlignedSlim]       -> SV.foldMap (bitVectToCharStringTNT localAlphabet) $ snd3 $ alignedSlimPrelim inCharData
                      x | x `elem` [AlignedWide]       -> UV.foldMap (bitVectToCharStringTNT localAlphabet) $ snd3 $ alignedWidePrelim inCharData
                      x | x `elem` [AlignedHuge]       ->    foldMap (bitVectToCharStringTNT localAlphabet) $ snd3 $ alignedHugePrelim inCharData 
                      _                                -> error ("Un-implemented data type " ++ show inCharType)
    in
    charString
        
-- | bitVectToCharStringTNT wraps '[]' around ambiguous states and removes commas between states
bitVectToCharStringTNT ::  Bits b => Alphabet String -> b -> String
bitVectToCharStringTNT localAlphabet bitValue = 
    let stateString = U.bitVectToCharState localAlphabet bitValue
    in
    if length stateString > 1 then "[" ++ (filter (/=',') stateString) ++ "]"
    else stateString

-- | bitVectToCharNumTNT wraps '[]' around ambiguous states and removes commas between states
bitVectToCharNumTNT ::  Bits b => Alphabet String -> b -> String
bitVectToCharNumTNT localAlphabet bitValue = 
    let stateString = U.bitVectToCharState localAlphabet bitValue
    in
    if length stateString > 1 then "[" ++ (filter (/=',') stateString) ++ "]"
    else stateString

-- | Implied Alignment report functions

-- | getImpliedAlignmentString returns as a single String the implied alignments of all sequence characters
-- softwired use display trees, hardWired transform to softwired then proceed with display trees
getImpliedAlignmentString :: GlobalSettings -> ProcessedData -> PhylogeneticGraph -> Int -> String
getImpliedAlignmentString inGS inData inGraph graphNumber =
    if LG.isEmpty (fst6 inGraph) then error "No graphs for create IAs for in getImpliedAlignmentStrings"
    else
        let headerString = "Implied Alignments for Graph " ++ (show graphNumber) ++ "\n"
        in
        if graphType inGS == Tree then headerString ++ (getTreeIAString inGraph)

        -- for softwired networks--use display trees
        else if graphType inGS == SoftWired then 
            -- get display trees for each data block-- takes first of potentially multiple
            let blockDisplayList = fmap GO.convertDecoratedToSimpleGraph $ fmap head $ fth6 inGraph

                -- create seprate processed data for each block
                blockProcessedDataList = fmap (makeBlockData (fst3 inData) (snd3 inData)) (thd3 inData)

                -- Perform full optimizations on display trees (as trees) with single block data (blockProcessedDataList) to creeate IAs
                decoratedBlockTreeList = V.zipWith (TRAV.multiTraverseFullyLabelGraph' (inGS {graphType = Tree}) False False Nothing) blockProcessedDataList blockDisplayList

                -- extract IA strings as if mutiple graphs
                diplayIAStringList = fmap getTreeIAString $ V.toList decoratedBlockTreeList

            in
            concat diplayIAStringList

        -- There is no IA for Hardwired at least as of yet 
        else 
            trace ("IA  not yet implemented for graphtype " ++ show (graphType inGS))
            "There is no implied alignment for hard-wired graphs--at least not yet.\n\tCould transform graph to softwired and generate an implied alignment that way"

-- | getTreeIAString takes a Tree Decorated Graph and returns Implied ALignmentString
getTreeIAString :: PhylogeneticGraph -> String
getTreeIAString inGraph =
    let leafList = snd4 $ LG.splitVertexList (thd6 inGraph)
        leafNameList = fmap (vertName . snd) leafList
        leafDataList = V.fromList $ fmap (vertData . snd) leafList
        charInfoVV = six6 inGraph
        characterStringList = makeFullIAStrings charInfoVV leafNameList leafDataList
    in
    concat characterStringList

-- | makeBlockData cretes new single block processed data
makeBlockData :: V.Vector NameText-> V.Vector NameBV -> BlockData -> ProcessedData
makeBlockData a b c = (a, b, V.singleton c)

-- | makeFullIAStrings goes block by block, creating fasta strings for each
makeFullIAStrings ::  V.Vector (V.Vector CharInfo) -> [NameText] -> V.Vector VertexBlockData -> [String]
makeFullIAStrings charInfoVV leafNameList leafDataList = 
    let numBlocks = V.length charInfoVV
    in
    concat $ fmap (makeBlockIAStrings leafNameList leafDataList charInfoVV) (V.fromList [0.. numBlocks - 1])

-- | makeBlockIAStrings extracts data for a block (via index) and calls funciton to make iaStrings for each character
makeBlockIAStrings :: [NameText] -> V.Vector (V.Vector (V.Vector CharacterData)) -> V.Vector (V.Vector CharInfo) -> Int -> [String]
makeBlockIAStrings leafNameList leafDataList charInfoVV blockIndex =
    let thisBlockCharInfo = charInfoVV V.! blockIndex
        numChars = V.length thisBlockCharInfo
        thisBlockCharData = fmap (V.! blockIndex) leafDataList
        blockCharacterStringList = V.zipWith (makeBlockCharacterString leafNameList thisBlockCharData) thisBlockCharInfo (V.fromList [0 .. (numChars - 1)])
    in
    filter (/= []) $ V.toList blockCharacterStringList

-- | makeBlockCharacterString creates implied alignmennt string for sequnce charactes and null if not
makeBlockCharacterString :: [NameText] -> V.Vector (V.Vector CharacterData) -> CharInfo -> Int -> String
makeBlockCharacterString leafNameList leafDataVV thisCharInfo charIndex =
    -- check if sequence type character
    let thisCharType = charType thisCharInfo
        thisCharName = name thisCharInfo
    in
    if thisCharType `notElem` sequenceCharacterTypes then []
    else 
        let -- thisCharData = fmap (V.! charIndex) leafDataVV
            thisCharData = getTaxDataOrMissing leafDataVV charIndex 0 []
            nameDataPairList = zip leafNameList thisCharData 
            fastaString = pairList2Fasta thisCharInfo nameDataPairList
        in 
        -- trace ("MBCS: " ++ (show $ length leafNameList) ++ " " ++ (show $ V.length thisCharData) ++ "\n" ++ (show leafDataVV))
        "\nSequence character " ++ (T.unpack thisCharName) ++ "\n" ++ fastaString ++ "\n"

{-
-- | getCharacterDataOrMissing takes a vector of vector of charcter data and returns list
-- of taxa for a given sequnce character.  If there are no data for that character for a taxon
getCharacterDataOrMissing :: V.Vector (V.Vector CharacterData) -> Int -> [[CharacterData]] -> [[CharacterData]]
getCharacterDataOrMissing leafDataVV charIndex newCharList =
    if charIndex == V.length leafDataVV then reverse newCharList
    else 
        let firstCharData = getTaxDataOrMissing leafDataVV charIndex 0 []
        in
        getCharacterDataOrMissing leafDataVV (charIndex + 1) (firstCharData : newCharList)
-}

-- | getTaxDataOrMissing gets teh index character if data not null, empty character if not
getTaxDataOrMissing :: V.Vector (V.Vector CharacterData) -> Int -> Int -> [CharacterData] -> [CharacterData]  
getTaxDataOrMissing charDataV charIndex taxonIndex newTaxList =
    if taxonIndex == V.length charDataV then reverse newTaxList
    else if V.null (charDataV V.! taxonIndex) then getTaxDataOrMissing charDataV charIndex (taxonIndex + 1) (emptyCharacter : newTaxList)
    else getTaxDataOrMissing charDataV charIndex (taxonIndex + 1) (((charDataV V.! taxonIndex) V.! charIndex) : newTaxList)

-- | pairList2Fasta takes a character type and list of pairs of taxon names (as T.Text) 
-- and character data and returns fasta formated string
pairList2Fasta :: CharInfo -> [(NameText, CharacterData)] -> String
pairList2Fasta inCharInfo nameDataPairList = 
    if null nameDataPairList then []
    else 
        let (firstName, blockDatum) = head nameDataPairList
            inCharType = charType inCharInfo
            localAlphabet = fmap ST.toString $ alphabet inCharInfo
            sequenceString = case inCharType of
                               x | x `elem` [SlimSeq, NucSeq  ] -> SV.foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ slimAlignment blockDatum
                               x | x `elem` [WideSeq, AminoSeq] -> UV.foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ wideAlignment blockDatum
                               x | x `elem` [HugeSeq]           ->    foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ hugeAlignment blockDatum
                               x | x `elem` [AlignedSlim]       -> SV.foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ alignedSlimPrelim blockDatum
                               x | x `elem` [AlignedWide]       -> UV.foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ alignedWidePrelim blockDatum
                               x | x `elem` [AlignedHuge]       ->    foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ alignedHugePrelim blockDatum 
                               _                                -> error ("Un-implemented data type " ++ show inCharType)

            sequenceChunks = fmap (++ "\n") $ SL.chunksOf 50 sequenceString

        in
        if blockDatum == emptyCharacter then (pairList2Fasta inCharInfo (tail nameDataPairList))
        else (concat $ (('>' : (T.unpack firstName)) ++ "\n") : sequenceChunks) ++ (pairList2Fasta inCharInfo (tail nameDataPairList))

