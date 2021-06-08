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

module Commands.CommandExecution ( executeCommands
                                 , executeRenameCommands) where

import           Types.Types
import           Debug.Trace
import           GeneralUtilities
import           System.IO
import           Data.Maybe
import           Data.List
import qualified Data.Text.Lazy as T
import qualified Data.CSV as CSV
import qualified Data.Vector as V
import qualified Data.Text.Short as ST
import qualified Graphs.GraphOperations as GO
import           GraphFormatUtilities
import qualified Utilities.LocalGraph as LG
import qualified Data.BitVector.LittleEndian as BV

-- | executeCommands reads iput files and returns raw data 
-- need to close files after read
executeCommands :: GlobalSettings -> [RawData] -> ProcessedData -> [PhylogeneticGraph] -> [[VertexCost]] -> [Command] -> IO ([PhylogeneticGraph], GlobalSettings)
executeCommands globalSettings rawData processedData curGraphs pairwiseDist commandList = do
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
            else if outFile == "stdout" then hPutStr stdout reportString
            else if writeMode == "overwrite" then writeFile outFile reportString
            else if writeMode == "append" then appendFile outFile reportString
            else error ("Error 'read' command not properly formatted" ++ (show reportStuff))
            executeCommands globalSettings rawData processedData curGraphs pairwiseDist (tail commandList)
        else if firstOption == Set then
            let (newGlobalSettings, newProcessedData) = setCommand firstArgs globalSettings processedData
            in
            executeCommands newGlobalSettings rawData newProcessedData curGraphs pairwiseDist (tail commandList)
        else 
            executeCommands globalSettings rawData processedData curGraphs pairwiseDist (tail commandList)

-- | setArgLIst contains valid 'set' arguments
setArgList :: [String]
setArgList = ["outgroup", "criterion", "graphtype","block"]

-- | setCommand takes arguments to change globalSettings and multiple data aspects (e.g. 'blocks')
setCommand :: [Argument] -> GlobalSettings -> ProcessedData -> (GlobalSettings, ProcessedData)
setCommand argList globalSettings processedData =
    let commandList = filter (/= "") $ fmap fst argList
        optionList = filter (/= "") $ fmap snd argList
        checkCommandList = checkCommandArgs "set" commandList setArgList
        leafNameVect = fst3 processedData

    in
    if checkCommandList == False then errorWithoutStackTrace ("Unrecognized command in 'set': " ++ (show argList))
    else 
        if head commandList == "outgroup"  then 
            let outTaxonName = T.pack $ filter (/= '"') $ head optionList
                outTaxonIndex = V.elemIndex outTaxonName leafNameVect

            in
            if outTaxonIndex == Nothing then errorWithoutStackTrace ("Error in 'set' command. Out-taxon " ++ (T.unpack outTaxonName) ++ " not found in input leaf list" ++ (show $ fmap (T.unpack) leafNameVect))
            else trace ("Outgroup set to " ++ T.unpack outTaxonName) (globalSettings {outgroupIndex = fromJust outTaxonIndex, outGroupName = outTaxonName}, processedData)
        else if head commandList == "graphtype"  then 
            let localGraphType = if (head optionList == "tree") then Tree
                                 else if (head optionList == "softwired") then SoftWired
                                 else if (head optionList == "hardwired") then HardWired
                                 else errorWithoutStackTrace ("Error in 'set' command. Graphtype '" ++ (head optionList) ++ "' is not 'tree', 'hardwired', or 'softwired'")
            in
            trace ("Graphtype set to " ++ (head optionList))
            (globalSettings {graphType = localGraphType}, processedData)
        else if head commandList == "criterion"  then 
            let localCriterion = if (head optionList == "parsimony") then Parsimony
                                 else if (head optionList == "pmdl") then PMDL
                                 else errorWithoutStackTrace ("Error in 'set' command. Criterion '" ++ (head optionList) ++ "' is not 'parsimony' or 'pmdl'")
            in
            trace ("Optimality criterion set to " ++ (head optionList))
            (globalSettings {optimalityCriterion = localCriterion}, processedData)
        else trace ("Warning--unrecognized/missing 'set' option in " ++ show argList) (globalSettings, processedData)



            
-- | reportArgList contains valid 'report' arguments
reportArgList :: [String]
reportArgList = ["all", "data", "graphs", "overwrite", "append", "dot", "newick", "ascii", "crossrefs", "pairdist", "diagnosis"]

-- | checkCommandArgs takes comamnd and args and verifies that they are in list
checkCommandArgs :: String -> [String] -> [String] -> Bool
checkCommandArgs commandString commandList permittedList =
    if null commandList then True
    else 
        let firstCommand = head commandList
            foundCommand = firstCommand `elem` permittedList
        in
        if foundCommand then checkCommandArgs commandString (tail commandList) permittedList
        else 
            let errorMatch = snd $ getBestMatch (maxBound :: Int ,"no suggestion") permittedList firstCommand
            in
            errorWithoutStackTrace ("\nError: Unrecognized '"++ commandString ++"' option. By \'" ++ firstCommand ++ "\' did you mean \'" ++ errorMatch ++ "\'?\n") 

-- | reportCommand takes report options, current data and graphs and returns a 
-- (potentially large) String to print and the channel to print it to 
-- and write mode overwrite/append
reportCommand :: GlobalSettings -> [Argument] -> [RawData] -> ProcessedData -> [PhylogeneticGraph] -> [[VertexCost]] -> (String, String, String)
reportCommand globalSettings argList rawData processedData curGraphs pairwiseDistanceMatrix =
    let outFileNameList = filter (/= "") $ fmap snd argList
        commandList = filter (/= "") $ fmap fst argList
    in
    if length outFileNameList > 1 then errorWithoutStackTrace ("Report can only have one file name: " ++ show outFileNameList)
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
        if checkCommandList == False then errorWithoutStackTrace ("Unrecognized command in report: " ++ (show argList))
        else 
            -- This for reconciled data
            if "crossrefs" `elem` commandList then
                let dataString = CSV.genCsvFile $ getDataListList rawData processedData
                in 
                (dataString, outfileName, writeMode)
            else if "data" `elem` commandList then 
                let dataString = phyloDataToString 0 $ thd3 processedData
                    baseData = ("There were " ++ (show $ length rawData) ++ " input data files with " ++ (show $ length $ thd3 processedData) ++ " blocks and " ++ (show $ ((length dataString) - 1)) ++ " total characters\n")
                    charInfoFields = ["Index", "Block", "Name", "Type", "Activity", "Weight", "Prealigned", "Alphabet", "TCM"]
                in
                (baseData ++ (CSV.genCsvFile $ charInfoFields : dataString), outfileName, writeMode)
            else if "diagnosis" `elem` commandList then
                let dataString = CSV.genCsvFile $ concat $ fmap (getGraphDiagnosis processedData) $ zip curGraphs [0.. ((length curGraphs) - 1)]
                in 
                (dataString, outfileName, writeMode)
            else if "graphs" `elem` commandList then 
                -- need to specify -O option for multiple graphs
                if "dot" `elem` commandList then
                    let graphString = concat $ intersperse "\n" $ fmap fgl2DotString $ fmap (GO.rerootGraph (outgroupIndex globalSettings)) $ fmap fst6 curGraphs
                    in 
                    (graphString, outfileName, writeMode)
                else if "newick" `elem` commandList then
                    let graphString = fglList2ForestEnhancedNewickString (fmap (GO.rerootGraph (outgroupIndex globalSettings)) (fmap GO.convertDecoratedToSimpleGraph $ fmap thd6 curGraphs))  True True
                        --graphString = fglList2ForestEnhancedNewickString (fmap (GO.rerootGraph (outgroupIndex globalSettings)) (fmap fst6 curGraphs))  True True
                        newickStringList = fmap init $ lines graphString
                        costStringList  = fmap ('[' :) $ fmap (++ "];\n") $ fmap show $ fmap snd6 curGraphs
                        graphStringCost = concat $ zipWith (++) newickStringList costStringList
                    in
                    (graphStringCost, outfileName, writeMode)
                else if "ascii" `elem` commandList then
                    let graphString = concat $ fmap LG.prettify  $ fmap (GO.rerootGraph (outgroupIndex globalSettings)) $ fmap GO.convertDecoratedToSimpleGraph $ fmap thd6 curGraphs 
                    in 
                    (graphString, outfileName, writeMode)
                else -- "dot" as default
                    let graphString = concat $ fmap fgl2DotString $ fmap (GO.rerootGraph (outgroupIndex globalSettings)) $ fmap fst6 curGraphs
                    in 
                    (graphString, outfileName, writeMode)
            else if "pairdist" `elem` commandList then
                let nameData = (intercalate "," $ V.toList $ fmap T.unpack $ fst3 processedData) ++ "\n"
                    dataString = CSV.genCsvFile $ fmap (fmap show) pairwiseDistanceMatrix
                in 
                (nameData ++ dataString, outfileName, writeMode)
            else trace ("Warning--unrecognized/missing report option in " ++ show commandList) ("No report specified", outfileName, writeMode)

-- | getDataListList returns a list of lists of Strings for data output as csv
-- for row is source file names, suubsequent rows by taxon with +/- for present absent taxon in 
-- input file
getDataListList :: [RawData] -> ProcessedData -> [[String]]
getDataListList inDataList processedData = 
    if null inDataList then []
    else 
        let fileNames = " " : (fmap (takeWhile (/= ':')) $ fmap T.unpack $ fmap name $ fmap head $ fmap snd inDataList)
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
    [T.unpack $ name inChar, show $ charType inChar, activityString, show $ weight inChar, prealignedString] ++ (fmap ST.toString $ alphabet inChar) ++ [show $ costMatrix inChar]

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
            let newName = T.filter (/= '"') $ T.pack $ snd $ head firstArgs
                newNameList = replicate (length $ tail firstArgs) newName
                newPairs = zip newNameList (fmap (T.filter (/= '"')) $ fmap T.pack $ fmap snd $ tail firstArgs)
            in
            executeRenameCommands (curPairs ++ newPairs) (tail commandList)

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
            vertexInfoList =  concat $ fmap (getVertexCharInfo (thd3 inData) (fst6 inGraph) (six6 inGraph)) vertexList
            edgeHeaderList = [[""],["", "Edge Head Vertex", "Edge Tail Vertex", "Edge Type", "Minimum Length", "Maximum Length"]]
            edgeInfoList = fmap getEdgeInfo edgeList
        in
        [topHeaderList, [show graphIndex]] ++ vertexInfoList ++ edgeHeaderList ++ edgeInfoList

-- | getVertexCharInfo returns a list of list of Strings of vertex infomation
-- one list for each character at the vertex
getVertexCharInfo :: V.Vector BlockData -> SimpleGraph -> V.Vector (V.Vector CharInfo) -> LG.LNode VertexInfo -> [[String]]
getVertexCharInfo blockDataVect inGraph charInfoVectVect inVert = 
    let leafParents = LG.parents inGraph (fst inVert)
        parentNodes = if nodeType  (snd inVert) == RootNode then "None" 
                     else if nodeType  (snd inVert) == LeafNode then show leafParents 
                     else show $  parents  (snd inVert)
        childNodes = if nodeType  (snd inVert) == LeafNode then "None" else show $  children  (snd inVert)
        basicInfoList = ["", show $ fst inVert, T.unpack $ vertName (snd inVert), show $ nodeType  (snd inVert), childNodes, parentNodes, "", "", "", "", "", show $ vertexCost (snd inVert)]
        blockCharVect = V.zip3  (V.map fst3 blockDataVect)  (vertData  (snd inVert)) charInfoVectVect 
        blockInfoList = concat $ V.toList $ V.map getBlockList blockCharVect
    in
    basicInfoList : blockInfoList

-- | getBlockList takes a pair of Vector of chardata and vector of charInfo and returns Strings
getBlockList :: (NameText, V.Vector CharacterData, V.Vector CharInfo) -> [[String]]
getBlockList (blockName, blockDataVect, charInfoVect) = 
    let firstLine = ["", "", "", "", "", "", T.unpack blockName]
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
        enhancedCharType = if localType `elem`  [SmallAlphSeq, NucSeq, AminoSeq, GenSeq] then (isPrealigned ++ (show localType))
                           else (show localType)
        (stringPrelim, stringFinal) = if localType == Add then (show $ rangePrelim blockDatum, show $ rangeFinal blockDatum) 
                                      else if localType == NonAdd then (concat $ V.map (bitVectToCharState localAlphabet) $ stateBVPrelim blockDatum, concat $ V.map (bitVectToCharState localAlphabet) $ stateBVFinal blockDatum)
                                      else if localType == Matrix then (show $ matrixStatesPrelim blockDatum, show $ matrixStatesFinal blockDatum)
                                      else if localType `elem` [SmallAlphSeq, NucSeq, AminoSeq, GenSeq] then (show $ sequencePrelim blockDatum, show $ sequenceFinal blockDatum)
                                      else error ("Un-implemented data type " ++ show localType)
        in
        ["", "", "", "", "", "", "", T.unpack $ name charInfo, enhancedCharType, stringPrelim, stringFinal, show $ localCost blockDatum]


-- | getEdgeInfo returns a list of Strings of edge infomation
getEdgeInfo :: LG.LEdge EdgeInfo -> [String]
getEdgeInfo inEdge = 
    ["", show $ fst3 inEdge, show $ snd3 inEdge, show $ edgeType (thd3 inEdge), show $ minLength (thd3 inEdge), show $ maxLength (thd3 inEdge)]


-- | executeSet processes the "set" command
-- set command very general can set outgroup, optimality criterion, blocks
-- executeSet :: 

-- bitVectToCharState  takes a bit vector representation and returns a list states as integers
bitVectToCharState :: [String] -> BV.BitVector -> String
bitVectToCharState localAlphabet inBit = 
    let bitBoolPairList = zip (BV.toBits inBit) localAlphabet
        (_, stateList) = unzip $ filter ((==True).fst) bitBoolPairList
        in
        intercalate "," stateList



