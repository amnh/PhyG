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
import System.Timing
import Control.DeepSeq
import           Control.Concurrent
import qualified Support.Support as SUP


-- | setArgLIst contains valid 'set' arguments
setArgList :: [String]
setArgList = ["outgroup", "criterion", "graphtype", "compressresolutions", "finalassignment", "graphfactor", "rootcost", "seed"]

-- | reportArgList contains valid 'report' arguments
reportArgList :: [String]
reportArgList = ["all", "data", "search", "graphs", "overwrite", "append", "dot", "newick", "ascii", "crossrefs", "pairdist", "diagnosis","displaytrees", "reconcile"]


-- | reconcileCommandList list of allowable commands
reconcileCommandList :: [String]
reconcileCommandList = ["method", "compare", "threshold", "outformat", "outfile", "connect", "edgelabel", "vertexlabel"]

-- | buildArgList is the list of valid build arguments
selectArgList :: [String]
selectArgList = ["best", "all", "unique", "atrandom"]

-- | executeCommands reads input files and returns raw data
-- need to close files after read
executeCommands :: GlobalSettings -> [RawData] -> ProcessedData -> [PhylogeneticGraph] -> [[VertexCost]] -> [Int] -> [Command] -> IO ([PhylogeneticGraph], GlobalSettings, [Int])
executeCommands globalSettings rawData processedData curGraphs pairwiseDist seedList commandList = do
    if null commandList then return (curGraphs, globalSettings, seedList)
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
            
            executeCommands (globalSettings {searchData = newSearchData}) rawData processedData (curGraphs ++ newGraphList) pairwiseDist (tail seedList) (tail commandList)
        
        else if firstOption == Refine then do 
            (elapsedSeconds, newGraphList) <- timeOp $ pure $ REF.refineGraph firstArgs globalSettings processedData (head seedList) curGraphs
            
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : (searchData globalSettings)   
            
            executeCommands (globalSettings {searchData = newSearchData}) rawData processedData newGraphList pairwiseDist (tail seedList) (tail commandList)
        
        else if firstOption == Fuse then do
            (elapsedSeconds, newGraphList) <- timeOp $ pure $ REF.fuseGraphs firstArgs globalSettings processedData (head seedList) curGraphs
            
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : (searchData globalSettings)   
            
            executeCommands (globalSettings {searchData = newSearchData}) rawData processedData newGraphList pairwiseDist (tail seedList) (tail commandList)
        
        else if firstOption == Report then do
            let reportStuff@(reportString, outFile, writeMode) = reportCommand globalSettings firstArgs rawData processedData curGraphs pairwiseDist
            hPutStrLn stderr ("Report writing to " ++ outFile)
            if outFile == "stderr" then hPutStr stderr reportString
            else if outFile == "stdout" then putStr reportString
            else if writeMode == "overwrite" then writeFile outFile reportString
            else if writeMode == "append" then appendFile outFile reportString
            else error ("Error 'read' command not properly formatted" ++ show reportStuff)
            executeCommands globalSettings rawData processedData curGraphs pairwiseDist seedList (tail commandList)
        
        else if firstOption == Search then do
            (elapsedSeconds, output) <- timeOp $ S.search firstArgs globalSettings processedData pairwiseDist (head seedList) curGraphs
                --in pure result
            -- (newGraphList, serchInfoList) <- S.search firstArgs globalSettings processedData pairwiseDist (head seedList) curGraphs
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs (fst output) (fromIntegral $ toMilliseconds elapsedSeconds) (concat $ fmap (L.intercalate "\n") $ snd output)
            let newSearchData = searchInfo : (searchData globalSettings)
            executeCommands (globalSettings {searchData = newSearchData})  rawData processedData  (fst output) pairwiseDist (tail seedList) (tail commandList)
        
        else if firstOption == Select then do
            (elapsedSeconds, newGraphList) <- timeOp $ pure $ GO.selectPhylogeneticGraph firstArgs (head seedList) selectArgList curGraphs
                
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : (searchData globalSettings)   
            
            executeCommands (globalSettings {searchData = newSearchData}) rawData processedData newGraphList pairwiseDist (tail seedList) (tail commandList)
        
        else if firstOption == Set then 
            -- if set changes graph aspects--may nned to reoptimize
            let (newGlobalSettings, newProcessedData, seedList') = setCommand firstArgs globalSettings processedData seedList
                newGraphList = if not (requireReoptimization globalSettings newGlobalSettings) then curGraphs
                               else trace ("Reoptimizing gaphs") fmap (TRA.multiTraverseFullyLabelGraph newGlobalSettings newProcessedData True True Nothing) (fmap fst6 curGraphs)
                
                searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList 0 "No Comment"
                newSearchData = searchInfo : (searchData newGlobalSettings)   
            in
            
            executeCommands (newGlobalSettings {searchData = newSearchData}) rawData newProcessedData newGraphList pairwiseDist seedList' (tail commandList)
        
        else if firstOption == Swap then do
            (elapsedSeconds, newGraphList) <- timeOp $ pure $ REF.swapMaster firstArgs globalSettings processedData (head seedList)  curGraphs
                
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : (searchData globalSettings)   
            
            executeCommands (globalSettings {searchData = newSearchData}) rawData processedData newGraphList pairwiseDist (tail seedList) (tail commandList)
        
        else if firstOption == Support then do
            (elapsedSeconds, newGraphList) <- timeOp $ pure $ SUP.supportGraph firstArgs globalSettings processedData (head seedList)  curGraphs
                
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : (searchData globalSettings)   
            
            executeCommands (globalSettings {searchData = newSearchData}) rawData processedData newGraphList pairwiseDist (tail seedList) (tail commandList)

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
        checkCommandList = U.checkCommandArgs "set" commandList setArgList
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
reportCommand :: GlobalSettings -> [Argument] -> [RawData] -> ProcessedData -> [PhylogeneticGraph] -> [[VertexCost]] -> (String, String, String)
reportCommand globalSettings argList rawData processedData curGraphs pairwiseDistanceMatrix =
    let argListWithoutReconcileCommands = filter ((`notElem` reconcileCommandList) .fst) argList
        outFileNameList = filter (/= "") $ fmap snd argListWithoutReconcileCommands --argList
        commandList = filter (/= "") $ fmap fst argListWithoutReconcileCommands
        -- reconcileList = filter (/= "") $ fmap fst argList
    in
    if length outFileNameList > 1 then errorWithoutStackTrace ("Report can only have one file name: " ++ (show outFileNameList) ++ " " ++ (show argList))
    else
        let checkCommandList = U.checkCommandArgs "report" commandList reportArgList
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
                    trace ("Reporting " ++ (show $ length curGraphs) ++ " graphs at minimum cost " ++ (show $ minimum $ fmap snd6 curGraphs))
                    (graphString, outfileName, writeMode)

           else if "pairdist" `elem` commandList then
                let nameData = L.intercalate "," (V.toList $ fmap T.unpack $ fst3 processedData) ++ "\n"
                    dataString = CSV.genCsvFile $ fmap (fmap show) pairwiseDistanceMatrix
                in
                (nameData ++ dataString, outfileName, writeMode)

            else if "reconcile" `elem` commandList then
                let (reconcileString, _) = R.makeReconcileGraph reconcileCommandList argList (fmap fst6 curGraphs)
                in
                if null curGraphs then 
                    trace ("No graphs to reconcile")
                    ("No graphs to reconcile", outfileName, writeMode)
                else 
                    (reconcileString, outfileName, writeMode)

            else if "search" `elem` commandList then
                let dataString = fmap showSearchFields $ reverse $ searchData globalSettings
                    baseData = ("SearchData\n")
                    charInfoFields = ["Command", "Arguments", "Min cost in", "Max cost in", "Num graphs in", "Min cost out", "Max cost out", "Num graphs out", "Duration (msecs)", "Comment"]
                in
                (baseData ++ CSV.genCsvFile (charInfoFields : dataString), outfileName, writeMode)


            else trace ("Warning--unrecognized/missing report option in " ++ show commandList) ("No report specified", outfileName, writeMode)

-- | showSearchFields cretes a String list for SearchData Fields
showSearchFields :: SearchData -> [String]
showSearchFields sD =
    [show $ instruction sD, concat $ fmap showArg $ arguments sD, show $ minGraphCostIn sD, show $ maxGraphCostIn sD, show $ numGraphsIn sD, show $ minGraphCostOut sD, show $ maxGraphCostOut sD, show $ numGraphsOut sD, 
    show $ duration sD, commentString sD]
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
  | "dot" `elem` commandList = makeDotList lOutgroupIndex graphList
  | "newick" `elem` commandList = makeNewickList lOutgroupIndex graphList (replicate (length graphList) 0.0)
  | "ascii" `elem` commandList = makeAsciiList lOutgroupIndex graphList
  | otherwise = -- "dot" as default
    makeDotList lOutgroupIndex graphList

-- | outputGraphString is a wrapper arounf graph output functions
outputGraphString :: [String] -> Int -> [DecoratedGraph] ->  [VertexCost] -> String
outputGraphString commandList lOutgroupIndex graphList costList
  | "dot" `elem` commandList = makeDotList lOutgroupIndex graphList
  | "newick" `elem` commandList = makeNewickList lOutgroupIndex graphList costList
  | "ascii" `elem` commandList = makeAsciiList lOutgroupIndex graphList
  | otherwise = -- "dot" as default
    makeDotList lOutgroupIndex graphList

-- | makeDotList takes a list of fgl trees and outputs a single String cointaining the graphs in Dot format
-- need to specify -O option for multiple graph(outgroupIndex globalSettings)s
makeDotList :: Int -> [DecoratedGraph] -> String
makeDotList rootIndex graphList =
     L.intercalate "\n" (fmap fgl2DotString $ fmap (GO.rerootTree rootIndex) $ fmap GO.convertDecoratedToSimpleGraph graphList)

-- | makeNewickList takes a list of fgl trees and outputs a single String cointaining the graphs in Newick format
makeNewickList ::  Int -> [DecoratedGraph] -> [VertexCost] -> String
makeNewickList rootIndex graphList costList =
    let graphString = fglList2ForestEnhancedNewickString (fmap (GO.rerootTree rootIndex . GO.convertDecoratedToSimpleGraph) graphList)  True True
        newickStringList = fmap init $ filter (not . null) $ lines graphString
        costStringList  = fmap (('[' :) . (++ "];\n")) (fmap show costList)
        graphStringCost = concat $ zipWith (++) newickStringList costStringList
    in
    graphStringCost

-- | makeAsciiList takes a list of fgl trees and outputs a single String cointaining the graphs in ascii format
makeAsciiList :: Int -> [DecoratedGraph] -> String
makeAsciiList rootIndex graphList =
    concatMap LG.prettify (fmap (GO.rerootTree rootIndex) $ fmap GO.convertDecoratedToSimpleGraph graphList)

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
        enhancedCharType = if localType `elem`  [SlimSeq, WideSeq, NucSeq, AminoSeq, HugeSeq] then (isPrealigned ++ (show localType))
                        else if localType `elem`  [Add, NonAdd, Matrix] then (show localType)
                        else error ("Character Type :" ++ (show localType) ++ "unrecogniized or not implemented")

        (stringPrelim, stringFinal) = if localType == Add then (show $ snd3 $ rangePrelim blockDatum, show $ rangeFinal blockDatum)
                                      else if localType == NonAdd then (concat $ V.map (U.bitVectToCharState localAlphabet) $ snd3 $ stateBVPrelim blockDatum, concat $ V.map (U.bitVectToCharState localAlphabet) $ stateBVFinal blockDatum)
                                      else if localType == Matrix then (show $ matrixStatesPrelim blockDatum, show $ fmap (fmap fst3) $ matrixStatesFinal blockDatum)
                                      else if localType `elem` [SlimSeq, WideSeq, NucSeq, AminoSeq, HugeSeq]
                                      then case localType of
                                             x | x `elem` [SlimSeq, NucSeq  ] -> (SV.foldMap (U.bitVectToCharState localAlphabet) $ slimPrelim blockDatum, SV.foldMap (U.bitVectToCharState localAlphabet) $ slimFinal blockDatum)
                                             x | x `elem` [WideSeq, AminoSeq] -> (UV.foldMap (U.bitVectToCharState localAlphabet) $ widePrelim blockDatum, UV.foldMap (U.bitVectToCharState localAlphabet) $ wideFinal blockDatum)
                                             x | x `elem` [HugeSeq]           -> (   foldMap (U.bitVectToCharState localAlphabet) $ hugePrelim blockDatum,    foldMap (U.bitVectToCharState localAlphabet) $ hugeFinal blockDatum)
                                             _                                -> error ("Un-implemented data type " ++ show localType)
                                      else error ("Un-implemented data type " ++ show localType)
        in
        [" ", " ", " ", " ", " ", " ", " ", T.unpack $ name charInfo, enhancedCharType, stringPrelim, stringFinal, show $ localCost blockDatum]


-- | getEdgeInfo returns a list of Strings of edge infomation
getEdgeInfo :: LG.LEdge EdgeInfo -> [String]
getEdgeInfo inEdge =
    [" ", show $ fst3 inEdge, show $ snd3 inEdge, show $ edgeType (thd3 inEdge), show $ minLength (thd3 inEdge), show $ maxLength (thd3 inEdge), show $ midRangeLength (thd3 inEdge)]


