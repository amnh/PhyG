{- |
Module      :  CommandUtilities.hs
Description :  Module helper fiunctions for command rpocessing
Copyright   :  (c) 2021-2022 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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



module Commands.CommandUtilities where


import           Data.Alphabet
import           Data.Bits
import qualified Data.Char                    as C
import           Data.Foldable
import qualified Data.List                    as L
import qualified Data.List.Split              as LS
import qualified Data.List.Split              as SL
import           Data.Maybe
import qualified Data.Text.Lazy               as T
import qualified Data.Text.Short              as ST
import qualified Data.Vector                  as V
import qualified Data.Vector.Storable         as SV
import qualified Data.Vector.Unboxed          as UV
import           Debug.Trace
import           GeneralUtilities
import           GraphFormatUtilities
import qualified GraphOptimization.Traversals as TRAV
import           Graphs.GraphOperations       as GO
import qualified Input.Reorganize             as IR
import qualified SymMatrix                    as S
import           System.Directory
import           System.IO
import           System.Info
import           System.Process
import           Types.Types
import qualified Utilities.LocalGraph         as LG
import qualified Utilities.Utilities          as U
-- import qualified Commands.Transform          as DT
import qualified Data.Set                     as SET


-- | processSearchFields takes a [String] and reformats the String associated with the
-- "search" commands and especially Thompson sampling data,
-- otherwise leaves the list unchanged
processSearchFields :: [[String]] -> [[String]]
processSearchFields inStringListList =
    if null inStringListList then []
    else
        let firstList = head inStringListList
        in
        if head firstList /= "Search" then firstList : processSearchFields (tail inStringListList)
        else
            let newHeader = ["Iteration","Search Type", "Delta", "Min Cost out", "CPU time (secs)"]
                instanceSplitList = LS.splitOn "*" (L.last firstList)
                (instanceStringListList, searchBanditListList) = unzip $ fmap processSearchInstance instanceSplitList -- (L.last firstList)
            in
            -- trace ("GSI: " ++ (show firstList))
            [L.init firstList] ++ [newHeader ++ head searchBanditListList ++ ["Arguments"]] ++ concat instanceStringListList ++ processSearchFields (tail inStringListList)

-- processSearchInstance takes the String of instance information and
-- returns appropriate [[String]] for pretty csv output
processSearchInstance :: String -> ([[String]], [String])
processSearchInstance inString =
    if null inString then ([], [])
    else
        let tempList = getSearchIterations inString
            iterationList = L.init tempList
            iterationCounterList = fmap ((:[]) . show) [0..(length iterationList - 1)]
            searchBanditList = getBanditNames $ tail $ dropWhile (/= '[') $ head iterationList
            preArgStringList = fmap getPreArgString iterationList
            searchArgStringList = fmap getSearchArgString iterationList
            searchBanditProbsList = fmap (getBanditProbs . tail . (dropWhile (/= '['))) iterationList
            processedSearchList = L.zipWith4 concat4 iterationCounterList preArgStringList  searchBanditProbsList searchArgStringList
            instanceStringList = processedSearchList ++ [LS.splitOn "," $ L.last tempList]
        in
        (instanceStringList, searchBanditList)
        where concat4 a b c d = a ++ b ++ c ++ d

-- | getBanditProbs parses bandit prob line for probabilities
getBanditProbs :: String -> [String]
getBanditProbs  inString =
    if null inString then []
    else
        let stuff = dropWhile (/= ',') inString
            stuff2 = tail $ takeWhile (/=  ')') stuff
            remainder = tail $ dropWhile (/=  ')') stuff
            remainder' = if length remainder > 2 then drop 2 remainder
                        else []
        in
        stuff2 : getBanditProbs remainder'

-- | getSearchArgString get seach iteration arguments and concats, removing ','
getSearchArgString :: String -> [String]
getSearchArgString inString =
    [L.intercalate " " $ drop 4 $ tail $ LS.splitOn "," $ takeWhile (/= '[') inString | not (null inString)]

-- getPreArgString gets search srtategy fields (type, delta, min cost, CPUtime)
-- befroe arguments fields
getPreArgString :: String -> [String]
getPreArgString inString =
    if null inString then []
    else
        take 4 $ tail $ LS.splitOn "," inString

-- | getBanditNames extracts the names of search bandits from comment list
-- already first part filtered out so only pairs in "(,)"
getBanditNames :: String -> [String]
getBanditNames  inString =
    if null inString then []
    else
        let firstBanditName = takeWhile (/= ',') $ tail inString
            remainder = dropWhile (/= '(') $ tail inString
        in
         firstBanditName : getBanditNames remainder

-- | getSearchIterations breaks up comment feild into individual iteration lines
getSearchIterations :: String -> [String]
getSearchIterations inList =
    if null inList then []
    else
        let commentField = filter (/= '"') inList
            commentLines = LS.splitOn "]" commentField
        in
        commentLines

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
--execute with "dot -Teps test.dot -o test.eps"
--need to add output to argument filename and call
--graphviz via System.Process.runprocess
--also, reorder GenForest so smalles (num leaves) is either first or
--last so can print small to large all the way so easier to read
-- eps on OSX because ps gets cutt off for some reason and no pdf onOSX
-- -O foir multiple graphs I htink
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
            _  <- createProcess (proc "dot" [outputType, dotFile, "-O"])
            hPutStrLn stderr ("\tExecuted dot " ++ outputType ++ " " ++ dotFile ++ " -O ")
        else
            hPutStrLn stderr "\tGraphviz call failed (not installed or found).  Dot file still created. Dot can be obtained from https://graphviz.org/download"


-- | showSearchFields cretes a String list for SearchData Fields
-- special processing for notACommand that contains info for initial data and graph processing
showSearchFields :: SearchData -> [String]
showSearchFields sD =
    let inInstruction = instruction sD
        (instructionString, durationString, commentString') = if inInstruction /= NotACommand then
                                                (show $ instruction sD, show ((fromIntegral $ duration sD) / 1000 :: Double), commentString sD)
                                              else (commentString sD, show ((fromIntegral $ duration sD) / 1000000000000 :: Double), "No Comment")

    in
    [instructionString, unwords $ (showArg <$> arguments sD), show $ minGraphCostIn sD, show $ maxGraphCostIn sD, show $ numGraphsIn sD, show $ minGraphCostOut sD, show $ maxGraphCostOut sD, show $ numGraphsOut sD,
        durationString, commentString']
    where showArg a = fst a ++ ":" ++ snd a

-- | requireReoptimization checks if set command in globalk settings requires reoptimization of graphs due to change in
-- graph type, optimality criterion etc.
requireReoptimization :: GlobalSettings -> GlobalSettings -> Bool
requireReoptimization gsOld gsNew
  | graphType gsOld /= graphType gsNew = True
  | optimalityCriterion gsOld /= optimalityCriterion gsNew = True
  | finalAssignment gsOld /= finalAssignment gsNew = True
  | graphFactor gsOld /= graphFactor gsNew = True
  | rootCost gsOld /= rootCost gsNew = True
  | otherwise = False

-- | outputBlockTrees takes a PhyloGeneticTree and outputs BlockTrees
outputBlockTrees :: [String] -> [VertexCost] -> Int -> (String , V.Vector [DecoratedGraph]) -> String
outputBlockTrees commandList costList lOutgroupIndex (labelString, graphLV) =
    let blockIndexStringList = fmap (((++ "\n") . ("//Block " ++)) . show) [0..((V.length graphLV) - 1)]
        blockStrings = unlines (makeBlockGraphStrings commandList costList lOutgroupIndex <$> zip blockIndexStringList (V.toList graphLV))
    in
    labelString ++ blockStrings

-- | makeBlockGraphStrings makes individual block display trees--potentially multiple
makeBlockGraphStrings :: [String] -> [VertexCost] -> Int -> (String ,[DecoratedGraph]) -> String
makeBlockGraphStrings commandList costList lOutgroupIndex (labelString, graphL) =
    let diplayIndexString =("//Display Tree(s): " ++ show (length graphL) ++ "\n")
        displayString = (++ "\n") $ outputDisplayString commandList costList lOutgroupIndex graphL
    in
    labelString ++ diplayIndexString ++ displayString

-- | outputDisplayString is a wrapper around graph output functions--but without cost list
outputDisplayString :: [String] -> [VertexCost] -> Int -> [DecoratedGraph] -> String
outputDisplayString commandList costList lOutgroupIndex graphList
  | "dot" `elem` commandList = makeDotList costList lOutgroupIndex (fmap GO.convertDecoratedToSimpleGraph graphList)
  | "newick" `elem` commandList = GO.makeNewickList (not (elem "nobranchlengths" commandList)) (not (elem "nohtulabels" commandList)) lOutgroupIndex (fmap GO.convertDecoratedToSimpleGraph graphList) (replicate (length graphList) 0.0)
  | "ascii" `elem` commandList = makeAsciiList lOutgroupIndex (fmap GO.convertDecoratedToSimpleGraph graphList)
  | otherwise = -- "dot" as default
    makeDotList costList lOutgroupIndex (fmap GO.convertDecoratedToSimpleGraph graphList)

-- | outputGraphString is a wrapper around graph output functions
outputGraphString :: [String] -> Int -> [DecoratedGraph] ->  [VertexCost] -> String
outputGraphString commandList lOutgroupIndex graphList costList
  | "dot" `elem` commandList = makeDotList costList lOutgroupIndex (fmap GO.convertDecoratedToSimpleGraph graphList)
  | "newick" `elem` commandList = GO.makeNewickList (not (elem "nobranchlengths" commandList)) (not (elem "nohtulabels" commandList)) lOutgroupIndex (fmap GO.convertDecoratedToSimpleGraph graphList) costList
  | "ascii" `elem` commandList = makeAsciiList lOutgroupIndex (fmap GO.convertDecoratedToSimpleGraph graphList)
  | otherwise = -- "dot" as default
    makeDotList costList lOutgroupIndex (fmap GO.convertDecoratedToSimpleGraph graphList)

-- | outputGraphStringSimple is a wrapper around graph output functions
outputGraphStringSimple :: [String] -> Int -> [SimpleGraph] ->  [VertexCost] -> String
outputGraphStringSimple commandList lOutgroupIndex graphList costList
  | "dot" `elem` commandList = makeDotList costList lOutgroupIndex graphList
  | "newick" `elem` commandList = GO.makeNewickList True True lOutgroupIndex graphList costList
  | "ascii" `elem` commandList = makeAsciiList lOutgroupIndex graphList
  | otherwise = -- "dot" as default
    makeDotList costList lOutgroupIndex graphList


-- | makeDotList takes a list of fgl trees and outputs a single String cointaining the graphs in Dot format
-- need to specify -O option for multiple graph(outgroupIndex globalSettings)s
makeDotList :: [VertexCost] -> Int -> [SimpleGraph] -> String
makeDotList costList rootIndex graphList =
    let graphStringList = fmap (fgl2DotString . LG.rerootTree rootIndex) graphList
        costStringList = fmap (("\n//" ++) . show) costList
    in
    L.intercalate "\n" (zipWith (++) graphStringList costStringList)

-- | makeAsciiList takes a list of fgl trees and outputs a single String cointaining the graphs in ascii format
makeAsciiList :: Int -> [SimpleGraph] -> String
makeAsciiList rootIndex graphList =
    concatMap LG.prettify (fmap (LG.rerootTree rootIndex) graphList)

{- Older version wiht more data dependenncy
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
-}

-- | getDataListList returns a list of lists of Strings for data output as csv
-- for row is source file names, subsequent rows by taxon with +/- for present absent taxon in
-- input file
-- different from getDataListList in removeal or processed data requiremenrt replaced with taxan name list
getDataListList :: [RawData] -> [T.Text] -> [[String]]
getDataListList inDataList fullTaxList =
    if null inDataList then []
    else
        let fileNames = " " : fmap (takeWhile (/= ':') . T.unpack) (fmap name $ fmap head $ fmap snd inDataList)
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
    let activityString = if activity inChar then "active"
                         else "inactive"
        prealignedString = if prealigned inChar then "prealigned"
                         else "unaligned"
    in
    [T.unpack $ name inChar, show $ charType inChar, activityString, show $ weight inChar, prealignedString] ++ [init $ concat $ fmap (++ ",") $ fmap ST.toString . toList $ alphabet inChar] ++ [show $ costMatrix inChar]

-- | executeRenameReblockCommands takes all the "Rename commands" pairs and
-- creates a list of pairs of new name and list of old names to be converted
-- as Text
executeRenameReblockCommands :: Instruction -> [(T.Text, T.Text)] -> [Command] -> IO [(T.Text, T.Text)]
executeRenameReblockCommands thisInStruction curPairs commandList  =
    if null commandList then return curPairs
    else do
        let (firstOption, firstArgs) = head commandList

        -- skip "Read" and "Rename "commands already processed
        if firstOption /= thisInStruction then executeRenameReblockCommands thisInStruction curPairs (tail commandList)
        else
            let newName = T.filter C.isPrint $ T.filter (/= '"') $ T.pack $ snd $ head firstArgs
                newNameList = replicate (length $ tail firstArgs) newName
                oldNameList = fmap (T.filter (/= '"') . T.pack) (fmap snd $ tail firstArgs)
                newPairs = zip newNameList oldNameList
            in
            executeRenameReblockCommands thisInStruction (curPairs ++ newPairs) (tail commandList)

-- | getGraphDiagnosis creates basic for CSV of graph vertex and node information
-- nodes first then vertices
getGraphDiagnosis :: GlobalSettings -> ProcessedData -> (PhylogeneticGraph, Int) -> [[String]]
getGraphDiagnosis _ inData (inGraph, graphIndex) =
    let decGraph = thd6 inGraph
    in
    if LG.isEmpty decGraph then []
    else
        let useIA = True
            useDO = False
            vertexList = LG.labNodes decGraph
            edgeList = LG.labEdges decGraph

            -- Vertex final character states for currect node
            vertexTitle = ["Vertex Character State Information"]
            -- topHeaderList  = ["Graph Index", "Vertex Index", "Vertex Name", "Vertex Type", "Child Vertices", "Parent Vertices", "Data Block", "Character Name", "Character Type", "Preliminary State", "Final State", "Local Cost"]
            topHeaderList  = ["Graph Index", "Vertex Index", "Vertex Name", "Vertex Type", "Child Vertices", "Parent Vertices", "Data Block", "Character Name", "Character Type", "Final State"]
            vertexInfoList =  concatMap (getVertexCharInfo useDO (thd3 inData) (fst6 inGraph) (six6 inGraph)) vertexList

            -- this using IA fields to get changes
            vertexInfoListChanges =  concatMap (getVertexCharInfo useIA (thd3 inData) (fst6 inGraph) (six6 inGraph)) vertexList

            -- Edge length information
            edgeTitle = [[" "],["Edge Weight/Length Information"]]
            edgeHeaderList = [[" ", "Edge Head Vertex", "Edge Tail Vertex", "Edge Type", "Minimum Length", "Maximum Length", "MidRange Length"]]
            edgeInfoList = fmap getEdgeInfo edgeList

            vertexChangeTitle = [[" "],["Vertex Character Changes"], ["Graph Index", "Vertex Index", "Vertex Name", "Vertex Type", "Child Vertices", "Parent Vertices", "Data Block", "Character Name", "Character Type", "Parent Final State", "Node Final State", "Sequence Changes (position, parent final state, node final state)"]]

            vertexParentStateList =  fmap ((:[]) . last) (concatMap (getVertexAndParentCharInfo useIA (thd3 inData) (fst6 inGraph) (six6 inGraph) (V.fromList vertexList)) vertexList)
            vertexStateList = fmap (drop 9) vertexInfoListChanges

            -- process to change to lines of individual changes--basically a transpose
            vertexChangeListByPosition = fmap (getAlignmentBasedChanges' 0) (zip vertexParentStateList vertexStateList)

            -- putting parent states before current state
            vertexStateInfoList = fmap (take 9) vertexInfoListChanges

            vertexChangeList = L.zipWith4 concat4 vertexStateInfoList vertexParentStateList vertexStateList vertexChangeListByPosition

            -- filter out those that are the same states
            differenceList = removeNoChangeLines vertexChangeList
        in
        -- trace ("GGD: " ++ (show $ snd6 staticGraph))
        [vertexTitle, topHeaderList, [show graphIndex]] ++ vertexInfoList ++ edgeTitle ++ edgeHeaderList ++ edgeInfoList ++ vertexChangeTitle ++ differenceList
        where concat4 a b c d = a ++ b ++ c ++ d

-- | getAlignmentBasedChanges' takes two equal length implied Alignments and outputs list of element changes between the two
-- assumes single String in each list
getAlignmentBasedChanges' :: Int -> ([String], [String]) -> [String]
getAlignmentBasedChanges' index (a, b)
  | length a > 1 || length b < 1 = error ("Should only have length 1 lists here: " ++ (show (length a, length b)))
  | null a = []
  | otherwise = 
        -- empty spaces sometimes
        let stringList1 = filter (not . null) $ LS.splitOn (" ") (head a)
            stringList2 = filter (not . null) $ LS.splitOn (" ") (head b)
        in

        -- this so returns empty for non--sequence characters
        if null stringList1 then []

        else if length stringList1 /= length stringList2 then
                 error ("Unequal characters in parent and node state lists in getAlignmentBasedChanges': "
                 ++ (show (length stringList1, length stringList2) ++ "\n" ++ (show stringList1) ++ "\n" ++ (show stringList2)))
        else getAlignmentBasedChangesGuts index stringList1 stringList2


-- | getAlignmentBasedChangesGuts takes processed element lists and cretes string of changes
getAlignmentBasedChangesGuts :: Int -> [String] -> [String] -> [String]
getAlignmentBasedChangesGuts index a b
  | null a = []
  | (head a) == (head b) = getAlignmentBasedChangesGuts (index + 1) (tail a) (tail b)
  | otherwise = ((show index) ++ ":" ++ (head a) ++ "," ++ (head b)) : getAlignmentBasedChangesGuts (index + 1) (tail a) (tail b)


-- | getAlignmentBasedChanges takes two equal length implied Alignments and outputs list of element changes between the two
-- only working for nucleotide prealigned or not
getAlignmentBasedChanges :: Int -> ([String], [String]) -> [String]
getAlignmentBasedChanges index (a, b) =
    if null a then []
    else 
        -- empty spaces sometimes
        let string1 = if (head $ head a) == ' ' then tail $ head a
                      else head a
            string2 = if (head $ head b) == ' ' then tail $ head b
                      else head b
        in
        if null string1 then []
            
        -- this so returns empty for non--sequence characters
        else if (length string1 < 2) || (length string2 < 2) then []
        else if length string1 /= length string2 then []
        else 
            -- index stuff in case gets out of sync in list (missing data, no changes etc)
            if (take 1 string1) == (take 1 string2) then 
               getAlignmentBasedChanges (index + 1) ([tail string1], [tail string2])
            else 
                ((show index) ++ ":" ++ (take 1 string1) ++ "," ++ (take 1 string2)) : getAlignmentBasedChanges (index + 1) ([tail string1], [tail string2])

-- | removeNoChangeLines takes lines of vertex changes and removes lines where parent and child startes are the same
-- so missing or ambiguous in one and not the other will be maintained
removeNoChangeLines :: [[String]] -> [[String]]
removeNoChangeLines inStringList =
    if null inStringList then []
    else
        let parentState = head inStringList !! 9
            childState  = head inStringList !! 10
        in
        (if (parentState == " ") || (parentState /= childState) then (head inStringList) : removeNoChangeLines (tail inStringList) else removeNoChangeLines (tail inStringList))

-- | getVertexCharInfo returns a list of list of Strings of vertex information
-- one list for each character at the vertex
getVertexCharInfo :: Bool -> V.Vector BlockData -> SimpleGraph -> V.Vector (V.Vector CharInfo) -> LG.LNode VertexInfo -> [[String]]
getVertexCharInfo useIA blockDataVect inGraph charInfoVectVect inVert =
    let nodeParents = LG.parents inGraph (fst inVert)
        parentNodes
          | nodeType  (snd inVert) == RootNode = "None"
          | nodeType  (snd inVert) == LeafNode = show nodeParents
          | otherwise = show $  parents  (snd inVert)
        childNodes = if nodeType  (snd inVert) == LeafNode then "None" else show $  children  (snd inVert)
        basicInfoList = [" ", show $ fst inVert, T.unpack $ vertName (snd inVert), show $ nodeType  (snd inVert), childNodes, parentNodes, " ", " ", " ", " ", " "]
        blockCharVect = V.zip3  (V.map fst3 blockDataVect)  (vertData  (snd inVert)) charInfoVectVect
        blockInfoList = concat $ V.toList $ V.map (getBlockList useIA) blockCharVect
    in
    basicInfoList : blockInfoList


-- | getVertexAndParentCharInfo returns a list of list of Strings of vertex information
-- for child and its parent
getVertexAndParentCharInfo :: Bool -> V.Vector BlockData -> SimpleGraph -> V.Vector (V.Vector CharInfo) -> V.Vector (LG.LNode VertexInfo) -> LG.LNode VertexInfo -> [[String]]
getVertexAndParentCharInfo useIA blockDataVect inGraph charInfoVectVect allVertVect inVert =
    let nodeParents = LG.parents inGraph (fst inVert)
        parentNodes
          | nodeType  (snd inVert) == RootNode = "None"
          | nodeType  (snd inVert) == LeafNode = show nodeParents
          | otherwise = show $  parents  (snd inVert)
        childNodes = if nodeType  (snd inVert) == LeafNode then "None" else show $  children  (snd inVert)
        basicInfoList = [" ", show $ fst inVert, T.unpack $ vertName (snd inVert), show $ nodeType  (snd inVert), childNodes, parentNodes, " ",
         " ", " ", " ", " "]
        blockCharVectNode = V.zip3  (V.map fst3 blockDataVect)  (vertData  (snd inVert)) charInfoVectVect

        -- for root--gets its own values as parent--filtered out in diff list later
        blockCharVectParent = if parentNodes == "None" then blockCharVectNode
                              else V.zip3  (V.map fst3 blockDataVect)  (vertData  (snd $ allVertVect V.! head nodeParents)) charInfoVectVect
        blockInfoListParent = concat $ V.toList $ V.map (getBlockList useIA) blockCharVectParent
    in
    basicInfoList : blockInfoListParent


-- | getBlockList takes a pair of Vector of chardata and vector of charInfo and returns Strings
getBlockList :: Bool -> (NameText, V.Vector CharacterData, V.Vector CharInfo) -> [[String]]
getBlockList useIA (blockName, blockDataVect, charInfoVect) =
    let firstLine = [" ", " ", " ", " ", " ", " ", T.unpack blockName," ", " ", " ", " ", " "]
        charlines = V.toList $ V.map (makeCharLine useIA) (V.zip blockDataVect charInfoVect)
    in
    firstLine : charlines

-- | makeCharLine takes character data
-- will be less legible for optimized data--so should use a diagnosis
-- based on "naive" data for human legible output
-- need to add back-converting to observed states using alphabet in charInfo
-- nothing here for packed since not "entered"
-- useIA for using alignment fields for changes in diagnosis
-- is useIA == False then just printing final seqeunce and removes spaces
-- for single character sequences (e.g. DNA/Protein)
makeCharLine :: Bool -> (CharacterData, CharInfo) -> [String]
makeCharLine useIA (blockDatum, charInfo) =
    let localType = charType charInfo
        localAlphabet = (ST.toString <$> alphabet charInfo)
        isPrealigned = if prealigned charInfo then "Prealigned "
                       else ""
        enhancedCharType
          | localType `elem` sequenceCharacterTypes = (isPrealigned ++ (show localType))
          | localType `elem` exactCharacterTypes = (show localType)
          | otherwise = error ("Character Type :" ++ (show localType) ++ "unrecogniized or not implemented")

        -- set where get string from, IA for change lists
        (slimField, wideField, hugeField) = if useIA then (slimIAFinal blockDatum, wideIAFinal blockDatum, hugeIAFinal blockDatum)
                                            else (slimFinal blockDatum, wideFinal blockDatum, hugeFinal blockDatum)


        -- (stringPrelim, stringFinal) = if localType == Add then (show $ snd3 $ rangePrelim blockDatum, show $ rangeFinal blockDatum)
        stringFinal
          | localType == Add = (show $ rangeFinal blockDatum)
          | localType == NonAdd = (concat $ V.map (U.bitVectToCharStateQual localAlphabet) $ stateBVFinal blockDatum)
          | localType `elem` packedNonAddTypes = (UV.foldMap (U.bitVectToCharState localAlphabet) $ packedNonAddFinal blockDatum)
          | localType == Matrix = (show $ fmap (fmap fst3) $ matrixStatesFinal blockDatum)
          | localType `elem` sequenceCharacterTypes = case localType of
                                    x | x `elem` [NucSeq  ] -> (SV.foldMap (U.bitVectToCharState' localAlphabet) slimField)
                                    x | x `elem` [SlimSeq ] -> (SV.foldMap (U.bitVectToCharState localAlphabet) slimField)
                                    x | x `elem` [WideSeq] -> (UV.foldMap (U.bitVectToCharState localAlphabet) wideField)
                                    x | x `elem` [AminoSeq] -> (UV.foldMap (U.bitVectToCharState' localAlphabet) wideField)
                                    x | x `elem` [HugeSeq]           -> (foldMap (U.bitVectToCharState localAlphabet) hugeField)
                                    x | x `elem` [AlignedSlim]       -> if (isAlphabetDna $ alphabet charInfo) || (isAlphabetRna $ alphabet charInfo) then
                                                                            (SV.foldMap (U.bitVectToCharState' localAlphabet) $ alignedSlimFinal blockDatum)
                                                                        else
                                                                            (SV.foldMap (U.bitVectToCharState localAlphabet) $ alignedSlimFinal blockDatum)
                                    x | x `elem` [AlignedWide]       -> if (isAlphabetAminoAcid $ alphabet charInfo) then
                                                                            (UV.foldMap (U.bitVectToCharState' localAlphabet) $ alignedWideFinal blockDatum)
                                                                        else
                                                                            (UV.foldMap (U.bitVectToCharState localAlphabet) $ alignedWideFinal blockDatum)
                                    x | x `elem` [AlignedHuge]       -> (foldMap (U.bitVectToCharState localAlphabet) $ alignedHugeFinal blockDatum)

                                    _ -> error ("Un-implemented data type " ++ show localType)
          | otherwise = error ("Un-implemented data type " ++ show localType)

        -- this removes ' ' between elements if sequence elements are a single character (e.g. DNA)
        stringFinal'
          | useIA = stringFinal
          | localType `elem` [NucSeq,  AminoSeq] = filter (/= ' ') stringFinal
          | otherwise = let maxSymbolLength = maximum $ fmap length $ SET.toList (alphabetSymbols localAlphabet)
                        in
                        if maxSymbolLength > 1 then removeRecurrrentSpaces $ fmap nothingToNothing stringFinal
                        else filter (/= ' ') stringFinal

        in

        -- trace ("MCL:" ++ (show localType) ++ " " ++ stringFinal)
        -- [" ", " ", " ", " ", " ", " ", " ", T.unpack $ name charInfo, enhancedCharType, stringPrelim, stringFinal, show $ localCost blockDatum]
        [" ", " ", " ", " ", " ", " ", " ", T.unpack $ name charInfo, enhancedCharType, stringFinal']
        where nothingToNothing a = if a == '\8220' then '\00'
                                   else a

-- | removeRecurrrentSpaces removes spaces that are followed by spaces
removeRecurrrentSpaces :: String -> String
removeRecurrrentSpaces inString
  | null inString = []
  | length inString == 1 = inString
  | head inString == ' ' = if (head $ tail inString) == ' ' then removeRecurrrentSpaces (tail inString)
        else (head inString) : removeRecurrrentSpaces (tail inString)
  | otherwise = (head inString) : removeRecurrrentSpaces (tail inString)

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
            leafNameList = fmap ((++ "\t") . T.unpack) (fmap (vertName . snd) leafList)
            headerString = "xread\n'TNT data for Graph " ++ show graphNumber ++ " generated by PhylogeneticGraph (PhyG)\n\tSource characters:\n"
            finalString = "proc/;\n\n"
            numTaxa = V.length $ fst3 inData
            charInfoVV = six6 inGraph

            -- get character information in 3-tuples and list of lengths--to match lengths
            ccCodeInfo = getCharacterInfo charInfoVV

        in

        if graphType inGS == Tree then
            let leafDataList = V.fromList $ fmap (vertData . snd) leafList

                -- get character strings
                taxonCharacterStringList = V.toList $ fmap ((++ "\n") . getTaxonCharString charInfoVV) leafDataList
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
            headerString  ++ nameLengthString ++ "'\n" ++ show (sum charLengthList) ++ " " ++ show numTaxa ++ "\n"
                ++ nameCharStringList ++ ";\n" ++ ccCodeString ++ finalString


        -- for softwired networks--use display trees
        else if graphType inGS == SoftWired then

            -- get display trees for each data block-- takes first of potentially multiple
            let middleStuffString = createDisplayTreeTNT inGS inData inGraph
            in
            headerString ++ middleStuffString ++ finalString

        -- for hard-wired networks--transfoirm to softwired and use display trees
        else if graphType inGS == HardWired then
            let newGS = inGS {graphType = SoftWired}

                pruneEdges = False
                warnPruneEdges = False
                startVertex = Nothing

                newGraph = TRAV.multiTraverseFullyLabelGraph newGS inData pruneEdges warnPruneEdges startVertex (fst6 newGraph)

                middleStuffString = createDisplayTreeTNT inGS inData newGraph

            in
            trace "There is no implied alignment for hard-wired graphs--at least not yet. Ggenerating TNT text via softwired transformation"
            --headerString ++ nameLengthString ++ "'\n" ++ (show $ sum charLengthList) ++ " " ++ (show numTaxa) ++ "\n"
            --    ++ nameCharStringList ++ ";\n" ++ ccCodeString ++ finalString
            headerString ++ middleStuffString ++ finalString

        else
            trace ("TNT  not yet implemented for graphtype " ++ show (graphType inGS))
            ("There is no implied alignment for "++ show (graphType inGS))

-- | createDisplayTreeTNT take a softwired graph and creates TNT data string
createDisplayTreeTNT :: GlobalSettings -> ProcessedData -> PhylogeneticGraph -> String
createDisplayTreeTNT inGS inData inGraph =
    let leafList = snd4 $ LG.splitVertexList (thd6 inGraph)
        leafNameList = fmap ((++ "\t") . T.unpack) (fmap (vertName . snd) leafList)
        charInfoVV = six6 inGraph
        numTaxa = V.length $ fst3 inData
        ccCodeInfo = getCharacterInfo charInfoVV
        blockDisplayList = fmap (GO.convertDecoratedToSimpleGraph . head) (fth6 inGraph)

        -- create seprate processed data for each block
        blockProcessedDataList = fmap (makeBlockData (fst3 inData) (snd3 inData)) (thd3 inData)

        -- Perform full optimizations on display trees (as trees) with single block data (blockProcessedDataList) to creeate IAs
        decoratedBlockTreeList = V.zipWith (TRAV.multiTraverseFullyLabelGraph' (inGS {graphType = Tree}) False False Nothing) blockProcessedDataList blockDisplayList

        -- create leaf data by merging display graph block data (each one a phylogentic graph)
        (leafDataList, mergedCharInfoVV) = mergeDataBlocks (V.toList decoratedBlockTreeList) [] []

        -- get character strings
        taxonCharacterStringList = V.toList $ fmap ((++ "\n") . getTaxonCharString mergedCharInfoVV) leafDataList
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
    nameLengthString ++ "'\n" ++ show (sum charLengthList) ++ " " ++ show numTaxa ++ "\n"
                ++ nameCharStringList ++ ";\n" ++ ccCodeString




-- | pairListToStringList takes  alist of (String, Int) and a starting index and returns scope of charcter for leading comment
pairListToStringList :: [(String, Int)] -> Int -> [String]
pairListToStringList pairList startIndex =
    if null pairList then []
    else
        let (a, b) = head pairList
        in
        ("\t\t" ++ show startIndex ++ "-" ++ show (b + startIndex - 1) ++ " : " ++ a ++ "\n") : pairListToStringList (tail pairList) (startIndex + b)

-- | mergeDataBlocks takes a list of Phylogenetic Graphs (Trees) and merges the data blocks (each graph should have only 1)
-- and merges the charInfo Vectors returning data and charInfo
mergeDataBlocks :: [PhylogeneticGraph] -> [[V.Vector CharacterData]] -> [V.Vector CharInfo] -> (V.Vector (V.Vector (V.Vector CharacterData)), V.Vector (V.Vector CharInfo))
mergeDataBlocks inGraphList curDataList curInfoList =
    if null inGraphList then (V.fromList $ fmap (V.fromList . reverse) curDataList, V.fromList $ reverse curInfoList)
    else
        let firstGraph = head inGraphList
            firstTree = thd6 firstGraph
            firstCharInfo = V.head $ six6 firstGraph
            leafList = snd4 $ LG.splitVertexList firstTree

            -- since each graph has a single block--take head to get vector of characters
            leafCharacterList = V.toList $ fmap (V.head . vertData . snd) (V.fromList leafList)

            -- zip data for each taxon
            newDataList = if null curDataList then (:[]) <$> leafCharacterList
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
        (ccCodeString' ++ weightString' ++ costsString') ++ mergeCharInfoCharLength (tail codeList) (tail lengthList) (charIndex + charLength)



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
        (costMatrixType, _) = IR.getRecodingType inMatrix
        matrixString = if costMatrixType == "nonAdd" then ""
                       else makeMatrixString inAlph inMatrix
    in
    let codeTriple = case inCharType of
                      x | x == Add -> ("+", "", charWeightString)
                      x | x == NonAdd -> ("-", "", charWeightString)
                      x | x `elem` packedNonAddTypes   -> ("-", "", charWeightString)
                      x | x == Matrix -> ("(", matrixString, charWeightString)
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
    let elementList =  (ST.toString <$> toList inAlphabet)

        -- get element pairs (Strings)
        elementPairList = filter notDiag $ getListPairs elementList

        -- index pairs for accessing matrix
        elementIndexPairList = filter notDiag $ getListPairs [0 .. (length elementList - 1)]
        elementPairCosts = fmap (inMatrix S.!) elementIndexPairList

        -- make strings of form state_i / state_j cost ...
        costString = makeCostString elementPairList elementPairCosts

    in
    costString
    where notDiag (a,b) = a < b

-- | makeCostString takes list of state pairs and list of costs and creates tnt cost string
makeCostString :: [(String, String)] -> [Int] -> String
makeCostString namePairList costList =
    if null namePairList then []
    else
        let (a,b) = head namePairList
            c = head costList
        in
        (a ++ "/" ++ b ++ " " ++ show c ++ " ") ++ makeCostString (tail namePairList) (tail costList)

-- | getBlockLength returns a list of the lengths of all characters in a blocks
getBlockLength :: V.Vector CharacterData -> V.Vector CharInfo -> [Int]
getBlockLength inCharDataV inCharInfoV =
    -- trace ("GBL:" ++ (show $ V.zipWith U.getCharacterLength inCharDataV inCharInfoV))
    V.toList $ V.zipWith U.getCharacterLength inCharDataV inCharInfoV

-- | getBlockNames returns a list of the lengths of all characters in a blocks
getBlockNames :: V.Vector CharInfo -> [String]
getBlockNames inCharInfoV =
    -- trace ("GBL:" ++ (show $ V.zipWith U.getCharacterLength inCharDataV inCharInfoV))
    V.toList $ fmap (T.unpack . name) inCharInfoV


-- | getCharacterString returns a string of character states
-- need to add splace between (large alphabets etc)
-- local alphabet for charactes where that is input.  MAytrix and additive are integers
getCharacterString :: CharacterData -> CharInfo -> String
getCharacterString inCharData inCharInfo =
    let inCharType = charType inCharInfo
        localAlphabet = if inCharType /= NonAdd then ST.toString <$> alphabet inCharInfo
                        else fmap ST.toString discreteAlphabet
    in
    let charString = case inCharType of
                      x | x == NonAdd ->    foldMap (bitVectToCharStringTNT localAlphabet) $ snd3 $ stateBVPrelim inCharData
                      x | x `elem` packedNonAddTypes   -> UV.foldMap (bitVectToCharStringTNT localAlphabet) $ snd3 $ packedNonAddPrelim inCharData
                      x | x == Add ->    foldMap  U.additivStateToString $ snd3 $ rangePrelim inCharData
                      x | x == Matrix ->    foldMap  U.matrixStateToString  $ matrixStatesPrelim inCharData
                      x | x `elem` [SlimSeq, NucSeq  ] -> SV.foldMap (bitVectToCharStringTNT localAlphabet) $ snd3 $ slimAlignment inCharData
                      x | x `elem` [WideSeq, AminoSeq] -> UV.foldMap (bitVectToCharStringTNT localAlphabet) $ snd3 $ wideAlignment inCharData
                      x | x == HugeSeq           ->    foldMap (bitVectToCharStringTNT localAlphabet) $ snd3 $ hugeAlignment inCharData
                      x | x == AlignedSlim       -> SV.foldMap (bitVectToCharStringTNT localAlphabet) $ snd3 $ alignedSlimPrelim inCharData
                      x | x == AlignedWide       -> UV.foldMap (bitVectToCharStringTNT localAlphabet) $ snd3 $ alignedWidePrelim inCharData
                      x | x == AlignedHuge       ->    foldMap (bitVectToCharStringTNT localAlphabet) $ snd3 $ alignedHugePrelim inCharData
                      _                                -> error ("Un-implemented data type " ++ show inCharType)
        allMissing = not (any (/= '-') charString)
    in
    if not allMissing then charString
    else replicate (length charString) '?'

-- | bitVectToCharStringTNT wraps '[]' around ambiguous states and removes commas between states
bitVectToCharStringTNT ::  (FiniteBits b, Bits b) => Alphabet String -> b -> String
bitVectToCharStringTNT localAlphabet bitValue =
    let stateString = U.bitVectToCharState localAlphabet bitValue
    in
    if length stateString > 1 then "[" ++ filter (/=',') stateString ++ "]"
    else stateString

-- | Implied Alignment report functions

-- | getImpliedAlignmentString returns as a single String the implied alignments of all sequence characters
-- softwired use display trees, hardWired transform to softwired then proceed with display trees
getImpliedAlignmentString :: GlobalSettings 
                          -> Bool 
                          -> Bool 
                          -> ProcessedData 
                          -> PhylogeneticGraph 
                          -> Int 
                          -> String
getImpliedAlignmentString inGS includeMissing concatSeqs inData inGraph graphNumber =
    if LG.isEmpty (fst6 inGraph) then error "No graphs for create IAs for in getImpliedAlignmentStrings"
    else
        let headerString = "Implied Alignments for Graph " ++ show graphNumber ++ "\n"
        in
        if graphType inGS == Tree then
            if not concatSeqs then headerString ++ getTreeIAString includeMissing inGraph
            else headerString ++ U.concatFastas (getTreeIAString includeMissing inGraph)

        -- for softwired networks--use display trees
        else if graphType inGS == SoftWired then
            -- get display trees for each data block-- takes first of potentially multiple
            let blockDisplayList = fmap (GO.contractIn1Out1EdgesRename . GO.convertDecoratedToSimpleGraph . head) (fth6 inGraph)

                -- create seprate processed data for each block
                blockProcessedDataList = fmap (makeBlockData (fst3 inData) (snd3 inData)) (thd3 inData)

                -- Perform full optimizations on display trees (as trees) with single block data (blockProcessedDataList) to create IAs
                decoratedBlockTreeList = V.zipWith (TRAV.multiTraverseFullyLabelGraph' (inGS {graphType = Tree}) False False Nothing) blockProcessedDataList blockDisplayList

                -- extract IA strings as if mutiple graphs
                diplayIAStringList = (getTreeIAString includeMissing <$> V.toList decoratedBlockTreeList)

            in
            if not concatSeqs then headerString ++ concat diplayIAStringList
            else headerString ++ U.concatFastas (concat diplayIAStringList)

        -- There is no IA for Hardwired at least as of yet
        -- so convert to softwired and use display trees
        else if graphType inGS == HardWired then
            let newGS = inGS {graphType = SoftWired}

                pruneEdges = False
                warnPruneEdges = False
                startVertex = Nothing

                newGraph = TRAV.multiTraverseFullyLabelGraph newGS inData pruneEdges warnPruneEdges startVertex (fst6 newGraph)

                blockDisplayList = fmap (GO.convertDecoratedToSimpleGraph . head) (fth6 newGraph)

                -- create seprate processed data for each block
                blockProcessedDataList = fmap (makeBlockData (fst3 inData) (snd3 inData)) (thd3 inData)

                -- Perform full optimizations on display trees (as trees) with single block data (blockProcessedDataList) to creeate IAs
                decoratedBlockTreeList = V.zipWith (TRAV.multiTraverseFullyLabelGraph' (inGS {graphType = Tree}) False False Nothing) blockProcessedDataList blockDisplayList

                -- extract IA strings as if mutiple graphs
                diplayIAStringList = (getTreeIAString includeMissing <$> V.toList decoratedBlockTreeList)
            in
            trace "There is no implied alignment for hard-wired graphs--at least not yet. Transfroming to softwired and generate an implied alignment that way" (
            if not concatSeqs then headerString ++ concat diplayIAStringList
            else headerString ++ U.concatFastas (concat diplayIAStringList)
            )

        else
            trace ("IA  not yet implemented for graphtype " ++ show (graphType inGS))
            ("There is no implied alignment for " ++  show (graphType inGS))

-- | getTreeIAString takes a Tree Decorated Graph and returns Implied AlignmentString
getTreeIAString :: Bool -> PhylogeneticGraph -> String
getTreeIAString includeMissing inGraph =
    let leafList = snd4 $ LG.splitVertexList (thd6 inGraph)
        leafNameList = fmap (vertName . snd) leafList
        leafDataList = V.fromList $ fmap (vertData . snd) leafList
        charInfoVV = six6 inGraph
        characterStringList = makeFullIAStrings includeMissing charInfoVV leafNameList leafDataList
    in
    concat characterStringList

-- | makeBlockData cretes new single block processed data
makeBlockData :: V.Vector NameText-> V.Vector NameBV -> BlockData -> ProcessedData
makeBlockData a b c = (a, b, V.singleton c)

-- | makeFullIAStrings goes block by block, creating fasta strings for each
makeFullIAStrings ::  Bool -> V.Vector (V.Vector CharInfo) -> [NameText] -> V.Vector VertexBlockData -> [String]
makeFullIAStrings includeMissing charInfoVV leafNameList leafDataList =
    let numBlocks = V.length charInfoVV
    in
    concatMap (makeBlockIAStrings includeMissing leafNameList leafDataList charInfoVV) (V.fromList [0.. numBlocks - 1])

-- | makeBlockIAStrings extracts data for a block (via index) and calls function to make iaStrings for each character
makeBlockIAStrings :: Bool -> [NameText] -> V.Vector (V.Vector (V.Vector CharacterData)) -> V.Vector (V.Vector CharInfo) -> Int -> [String]
makeBlockIAStrings includeMissing leafNameList leafDataList charInfoVV blockIndex =
    let thisBlockCharInfo = charInfoVV V.! blockIndex
        numChars = V.length thisBlockCharInfo
        thisBlockCharData = fmap (V.! blockIndex) leafDataList
        blockCharacterStringList = V.zipWith (makeBlockCharacterString includeMissing leafNameList thisBlockCharData) thisBlockCharInfo (V.fromList [0 .. (numChars - 1)])
    in
    V.toList blockCharacterStringList

-- | isAllGaps checks whether a sequence is all gap charcaters '-'
isAllGaps :: String -> Bool
isAllGaps inSeq
  | null inSeq = True
  | length (filter (`notElem` ['-', '\n']) inSeq) == 0 = True
  | otherwise = False

-- | makeBlockCharacterString creates implied alignmennt string for sequnce charactes and null if not
makeBlockCharacterString :: Bool -> [NameText] -> V.Vector (V.Vector CharacterData) -> CharInfo -> Int -> String
makeBlockCharacterString includeMissing leafNameList leafDataVV thisCharInfo charIndex =
    -- check if sequence type character
    let thisCharType = charType thisCharInfo
        thisCharName = name thisCharInfo
    in
    if thisCharType `notElem` sequenceCharacterTypes then []
    else
        let -- thisCharData = fmap (V.! charIndex) leafDataVV
            thisCharData = getTaxDataOrMissing leafDataVV charIndex 0 []
            nameDataPairList = zip leafNameList thisCharData
            fastaString = pairList2Fasta includeMissing thisCharInfo nameDataPairList
        in
        -- trace ("MBCS: " ++ (show $ length leafNameList) ++ " " ++ (show $ V.length thisCharData) ++ "\n" ++ (show leafDataVV))
        "\nSequence character " ++ T.unpack thisCharName ++ "\n" ++ fastaString ++ "\n"

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
getTaxDataOrMissing charDataV charIndex taxonIndex newTaxList
  | taxonIndex == V.length charDataV = reverse newTaxList
  | V.null (charDataV V.! taxonIndex) = getTaxDataOrMissing charDataV charIndex (taxonIndex + 1) (emptyCharacter : newTaxList)
  | otherwise = getTaxDataOrMissing charDataV charIndex (taxonIndex + 1) (((charDataV V.! taxonIndex) V.! charIndex) : newTaxList)

-- | pairList2Fasta takes a character type and list of pairs of taxon names (as T.Text)
-- and character data and returns fasta formated string
pairList2Fasta :: Bool -> CharInfo -> [(NameText, CharacterData)] -> String
pairList2Fasta includeMissing inCharInfo nameDataPairList =
    if null nameDataPairList then []
    else
        let (firstName, blockDatum) = head nameDataPairList
            inCharType = charType inCharInfo
            localAlphabet = (ST.toString <$> alphabet inCharInfo)

            sequenceString = case inCharType of
                               -- x | x `elem` [SlimSeq, NucSeq  ] -> SV.foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ slimAlignment blockDatum
                               x | x == SlimSeq -> SV.foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ slimAlignment blockDatum
                               x | x == NucSeq -> SV.foldMap (U.bitVectToCharState' localAlphabet) $ snd3 $ slimAlignment blockDatum
                               -- x | x `elem` [WideSeq, AminoSeq] -> UV.foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ wideAlignment blockDatum
                               x | x == WideSeq -> UV.foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ wideAlignment blockDatum
                               x | x == AminoSeq -> UV.foldMap (U.bitVectToCharState' localAlphabet) $ snd3 $ wideAlignment blockDatum
                               x | x == HugeSeq           ->    foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ hugeAlignment blockDatum
                               x | x == AlignedSlim       -> if isAlphabetDna (alphabet inCharInfo) || isAlphabetRna (alphabet inCharInfo) then SV.foldMap (U.bitVectToCharState' localAlphabet) $ snd3 $ alignedSlimPrelim blockDatum
                                                                   else SV.foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ alignedSlimPrelim blockDatum
                               x | x == AlignedWide       -> if isAlphabetAminoAcid $ alphabet inCharInfo then UV.foldMap (U.bitVectToCharState' localAlphabet) $ snd3 $ alignedWidePrelim blockDatum
                                                                   else UV.foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ alignedWidePrelim blockDatum
                               x | x == AlignedHuge       ->    foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ alignedHugePrelim blockDatum
                               _                                -> error ("Un-implemented data type " ++ show inCharType)

            sequenceString' = if not (any (/= '-') sequenceString) then
                                replicate (length sequenceString) '?'
                              else sequenceString
            sequenceChunks = ((++ "\n") <$> SL.chunksOf 50 sequenceString')

        in
        (if ((not includeMissing) && (isAllGaps sequenceString)) || (blockDatum == emptyCharacter) then pairList2Fasta includeMissing inCharInfo (tail nameDataPairList) else (concat $ (('>' : (T.unpack firstName)) ++ "\n") : sequenceChunks) ++ (pairList2Fasta includeMissing inCharInfo (tail nameDataPairList)))


