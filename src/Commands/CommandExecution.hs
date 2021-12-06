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

module Commands.CommandExecution
  ( executeCommands
  , executeRenameReblockCommands
  ) where

import Data.Foldable
import qualified Data.CSV               as CSV
import qualified Data.List              as L
import           Data.Maybe
-- import           Text.Read
import qualified Data.Text.Lazy         as T
import qualified Data.Text.Short        as ST
import qualified Data.Vector            as V
import qualified Data.Vector.Storable   as SV
import qualified Data.Vector.Unboxed    as UV
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
import qualified Search.Swap            as SW
import qualified Reconciliation.ReconcileGraphs as R


-- | setArgLIst contains valid 'set' arguments
setArgList :: [String]
setArgList = ["outgroup", "criterion", "graphtype","compressresolutions", "finalassignment", "graphfactor", "rootcost"]

-- | reportArgList contains valid 'report' arguments
reportArgList :: [String]
reportArgList = ["all", "data", "graphs", "overwrite", "append", "dot", "newick", "ascii", "crossrefs", "pairdist", "diagnosis","displaytrees", "reconcile"]


-- | reconcileCommandList list of allowable commands
reconcileCommandList :: [String]
reconcileCommandList = ["method", "compare", "threshold", "outformat", "outfile", "connect", "edgelabel", "vertexlabel"]

-- | buildArgList is the list of valid build arguments
selectArgList :: [String]
selectArgList = ["best", "all", "unique", "random"]


-- | executeCommands reads input files and returns raw data
-- need to close files after read
executeCommands :: GlobalSettings -> [RawData] -> ProcessedData -> [PhylogeneticGraph] -> [[VertexCost]] -> [Int] -> [Command] -> IO ([PhylogeneticGraph], GlobalSettings)
executeCommands globalSettings rawData processedData curGraphs pairwiseDist seedList commandList = do
    if null commandList then return (curGraphs, globalSettings)
    else do
        let (firstOption, firstArgs) = head commandList

        -- skip "Read" and "Rename "commands already processed
        if firstOption == Read then error ("Read command should already have been processed: " ++ show (firstOption, firstArgs))
        else if firstOption == Rename then error ("Rename command should already have been processed: " ++ show (firstOption, firstArgs))
        -- report command
        else if firstOption == Report then do
            let reportStuff@(reportString, outFile, writeMode) = reportCommand globalSettings firstArgs rawData processedData curGraphs pairwiseDist
            hPutStrLn stderr ("Report writing to " ++ outFile)
            if outFile == "stderr" then hPutStr stderr reportString
            else if outFile == "stdout" then putStr reportString
            else if writeMode == "overwrite" then writeFile outFile reportString
            else if writeMode == "append" then appendFile outFile reportString
            else error ("Error 'read' command not properly formatted" ++ show reportStuff)
            executeCommands globalSettings rawData processedData curGraphs pairwiseDist seedList (tail commandList)
        else if firstOption == Set then
            let (newGlobalSettings, newProcessedData) = setCommand firstArgs globalSettings processedData
            in
            executeCommands newGlobalSettings rawData newProcessedData curGraphs pairwiseDist seedList (tail commandList)
        else if firstOption == Build then
            let newGraphList = B.buildGraph firstArgs globalSettings processedData pairwiseDist (head seedList)
            in
            executeCommands globalSettings rawData processedData (curGraphs ++ newGraphList) pairwiseDist (tail seedList) (tail commandList)
        else if firstOption == Select then
            let newGraphList = GO.selectPhylogeneticGraph firstArgs (head seedList) selectArgList curGraphs
            in
            executeCommands globalSettings rawData processedData newGraphList pairwiseDist (tail seedList) (tail commandList)
        else if firstOption == Swap then
            let newGraphList = SW.swapMaster firstArgs globalSettings processedData (head seedList)  curGraphs
            in
            executeCommands globalSettings rawData processedData newGraphList pairwiseDist (tail seedList) (tail commandList)
        else error ("Command " ++ (show firstOption) ++ " not recognized/implemented")

-- | setCommand takes arguments to change globalSettings and multiple data aspects (e.g. 'blocks')
setCommand :: [Argument] -> GlobalSettings -> ProcessedData -> (GlobalSettings, ProcessedData)
setCommand argList globalSettings processedData =
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
            else trace ("Outgroup set to " ++ T.unpack outTaxonName) (globalSettings {outgroupIndex = fromJust outTaxonIndex, outGroupName = outTaxonName}, processedData)
        else if head commandList == "graphtype"  then
            let localGraphType
                  | (head optionList == "tree") = Tree
                  | (head optionList == "softwired") = SoftWired
                  | (head optionList == "hardwired") = HardWired
                  | otherwise = errorWithoutStackTrace ("Error in 'set' command. Graphtype '" ++ (head optionList) ++ "' is not 'tree', 'hardwired', or 'softwired'")
            in
            trace ("Graphtype set to " ++ head optionList)
            (globalSettings {graphType = localGraphType}, processedData)
        else if head commandList == "criterion"  then
            let localCriterion
                  | (head optionList == "parsimony") = Parsimony
                  | (head optionList == "pmdl") = PMDL
                  | otherwise = errorWithoutStackTrace ("Error in 'set' command. Criterion '" ++ (head optionList) ++ "' is not 'parsimony' or 'pmdl'")
            in
            trace ("Optimality criterion set to " ++ head optionList)
            (globalSettings {optimalityCriterion = localCriterion}, processedData)
        else if head commandList == "compressresolutions"  then
            let localCriterion
                  | (head optionList == "true") = True
                  | (head optionList == "false") = False
                  | otherwise = errorWithoutStackTrace ("Error in 'set' command. CompressResolutions '" ++ (head optionList) ++ "' is not 'true' or 'false'")
            in
            trace ("CompressResolutions set to " ++ head optionList)
            (globalSettings {compressResolutions = localCriterion}, processedData)
        else if head commandList == "finalassignment"  then
            let localMethod
                  | ((head optionList == "do") || (head optionList == "directoptimization")) = DirectOptimization
                  | ((head optionList == "ia") || (head optionList == "impliedalignment")) = ImpliedAlignment
                  | otherwise = errorWithoutStackTrace ("Error in 'set' command. FinalAssignment  '" ++ (head optionList) ++ "' is not 'DirectOptimization (DO)' or 'ImpliedAlignment (IA)'")
            in
            trace ("FinalAssignment set to " ++ head optionList)
            (globalSettings {finalAssignment = localMethod}, processedData)
        else if head commandList == "graphfactor"  then
            let localMethod
                  | (head optionList == "nopenalty") = NoNetworkPenalty
                  | (head optionList == "w15") = Wheeler2015Network
                  | (head optionList == "pmdl") = PMDLGraph
                  | otherwise = errorWithoutStackTrace ("Error in 'set' command. GraphFactor  '" ++ (head optionList) ++ "' is not 'NoPenalty', 'W15', or 'PMDL'")
            in
            trace ("GraphFactor set to " ++ head optionList)
            (globalSettings {graphFactor = localMethod}, processedData)
        else if head commandList == "rootcost"  then
            let localMethod
                  | (head optionList == "norootcost") = NoRootCost
                  | (head optionList == "w15") = Wheeler2015Root
                  | (head optionList == "pmdl") = PMDLRoot
                  | otherwise = errorWithoutStackTrace ("Error in 'set' command. RootCost  '" ++ (head optionList) ++ "' is not 'NoRootCost', 'W15', or 'PMDL'")
            in
            trace ("RootCost set to " ++ head optionList)
            (globalSettings {rootCost = localMethod}, processedData)
        else trace ("Warning--unrecognized/missing 'set' option in " ++ show argList) (globalSettings, processedData)



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
                (dataString, outfileName, writeMode)

            else if "graphs" `elem` commandList then
                let graphString = outputGraphString commandList (outgroupIndex globalSettings) (fmap thd6 curGraphs) (fmap snd6 curGraphs)
                in
                trace ("Reporting " ++ (show $ length curGraphs) ++ " graphs at minimum cost " ++ (show $ minimum $ fmap snd6 curGraphs))
                (graphString, outfileName, writeMode)

            else if "displaytrees" `elem` commandList then
                -- need to specify -O option for multiple graphs
                let inputDisplayVVList = fmap fth6 curGraphs
                    treeIndexStringList = fmap ((++ "\n") . ("Canonical Tree " ++)) (fmap show [0..(length inputDisplayVVList - 1)])
                    canonicalGraphPairList = zip treeIndexStringList inputDisplayVVList
                    blockStringList = concatMap (++ "\n") (fmap (outputBlockTrees commandList (outgroupIndex globalSettings)) canonicalGraphPairList)
                    -- graphString = outputGraphString commandList (outgroupIndex globalSettings) (fmap thd6 curGraphs) (fmap snd6 curGraphs)
                in
                (blockStringList, outfileName, writeMode)

            else if "pairdist" `elem` commandList then
                let nameData = L.intercalate "," (V.toList $ fmap T.unpack $ fst3 processedData) ++ "\n"
                    dataString = CSV.genCsvFile $ fmap (fmap show) pairwiseDistanceMatrix
                in
                (nameData ++ dataString, outfileName, writeMode)

            else if "reconcile" `elem` commandList then
                let (reconcileString, _) = R.makeReconcileGraph reconcileCommandList argList (fmap fst6 curGraphs)
                in
                (reconcileString, outfileName, writeMode)


            else trace ("Warning--unrecognized/missing report option in " ++ show commandList) ("No report specified", outfileName, writeMode)

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

        (stringPrelim, stringFinal) = if localType == Add then (show $ rangePrelim blockDatum, show $ rangeFinal blockDatum)
                                      else if localType == NonAdd then (concat $ V.map (U.bitVectToCharState localAlphabet) $ fst3 $ stateBVPrelim blockDatum, concat $ V.map (U.bitVectToCharState localAlphabet) $ stateBVFinal blockDatum)
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


