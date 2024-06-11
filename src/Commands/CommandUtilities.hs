{- |
Module exposing helper functions used for command processing.
-}
module Commands.CommandUtilities where

import Bio.DynamicCharacter
import Bio.DynamicCharacter.Element
import Control.Monad (when)
import Control.Monad.IO.Class (MonadIO (..))
import Control.Parallel.Strategies
import Data.Alphabet
import Data.Alphabet.Codec (decodeState)
import Data.Alphabet.Special
import Data.BitVector.LittleEndian qualified as BV
import Data.Bits
import Data.Char qualified as C
import Complexity.Utilities qualified as CU
import Data.Foldable
import Data.List qualified as L
import Data.List.NonEmpty qualified as NE
import Data.List.Split qualified as LS
import Data.Maybe
import Data.MetricRepresentation qualified as MR
import Data.Set qualified as SET
import Data.Text.Lazy qualified as T
import Data.Text.Short qualified as ST
import Data.Vector qualified as V
import Data.Vector.Generic qualified as GV
import Data.Vector.Storable qualified as SV
import Data.Vector.Unboxed qualified as UV
import Debug.Trace
import DirectOptimization.Pairwise
import GeneralUtilities
import GraphFormatUtilities
import GraphOptimization.Medians qualified as M
import GraphOptimization.PreOrderFunctions qualified as PRE
import GraphOptimization.Traversals qualified as TRAV
import Graphs.GraphOperations qualified as GO
import Input.Reorganize qualified as IR
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import PHANE.Evaluation.Verbosity (Verbosity (..))
import SymMatrix qualified as S
import System.Directory
import System.IO
import System.Info
import System.Process
import Types.Types
import Utilities.LocalGraph qualified as LG
import Utilities.Utilities qualified as U


{- | processSearchFields takes a [String] and reformats the String associated with the
"search" commands and especially Thompson sampling data,
otherwise leaves the list unchanged
-}
processSearchFields ∷ [[String]] → [[String]]
processSearchFields inStringListList =
    if null inStringListList
        then []
        else
            let firstList = head inStringListList
            in  if head firstList /= "Search"
                    then firstList : processSearchFields (tail inStringListList)
                    else
                        let newHeader = ["Iteration", "Search Type", "Delta", "Min Cost out", "CPU time (secs)"]
                            instanceSplitList = LS.splitOn "*" (L.last firstList)
                            hitsMinimum = filter (/= '*') $ last $ LS.splitOn "," (L.last firstList)
                            (instanceStringListList, searchBanditListList) = unzip $ fmap processSearchInstance instanceSplitList -- (L.last firstList)
                        in  -- trace ("GSI: " <> (show firstList) <> "\nLF: " <> (hitsMinimum) <> "\nILL: " <> (show instanceStringListList) <> "\nSB: " <> (show searchBanditListList))
                            fmap (fmap (filter (/= '\n'))) $
                                [L.init firstList]
                                    <> [newHeader <> (head searchBanditListList) <> ["Arguments"]]
                                    <> concat instanceStringListList
                                    <> [[hitsMinimum]]
                                    <> processSearchFields (tail inStringListList)


-- processSearchInstance takes the String of instance information and
-- returns appropriate [[String]] for pretty csv output
processSearchInstance ∷ String → ([[String]], [String])
processSearchInstance inString =
    if null inString
        then ([], [])
        else
            let tempList = getSearchIterations inString
                iterationList = L.init tempList
                iterationCounterList = fmap ((: []) . show) [0 .. (length iterationList - 1)]
                searchBanditList' = getBanditNames $ tail $ dropWhile (/= '[') $ head iterationList
                preArgStringList = fmap getPreArgString iterationList
                searchArgStringList = fmap getSearchArgString iterationList
                searchBanditProbsList' = fmap (getBanditProbs . tail . (dropWhile (/= '['))) iterationList

                -- add columns between graph and search bandits
                searchBanditList = (take 3 searchBanditList') <> [" "] <> (drop 3 searchBanditList')
                searchBanditProbsList = fmap (addColumn 3) searchBanditProbsList'

                processedSearchList = L.zipWith4 concat4 iterationCounterList preArgStringList searchBanditProbsList searchArgStringList
                finalString = addColumn 8 $ LS.splitOn "," $ L.last tempList
                instanceStringList = processedSearchList <> [finalString]
            in  if null iterationList
                    then ([], [])
                    else (instanceStringList, searchBanditList)
    where
        concat4 ∷ ∀ {a}. (Semigroup a) ⇒ a → a → a → a → a
        concat4 a b c d = a <> b <> c <> d


-- | addColumn takea list of strings and adds a list of empty strings after third string in each
addColumn ∷ Int → [String] → [String]
addColumn index inList =
    if null inList
        then []
        else (take index inList) <> [" "] <> (drop index inList)


-- | getBanditProbs parses bandit prob line for probabilities
getBanditProbs ∷ String → [String]
getBanditProbs inString =
    if null inString
        then []
        else
            let stuff = dropWhile (/= ',') inString
                stuff2 = tail $ takeWhile (/= ')') stuff
                remainder = tail $ dropWhile (/= ')') stuff
                remainder' =
                    if length remainder > 2
                        then drop 2 remainder
                        else []
            in  stuff2 : getBanditProbs remainder'


-- | getSearchArgString get seach iteration arguments and concats, removing ','
getSearchArgString ∷ String → [String]
getSearchArgString inString =
    [L.intercalate " " $ drop 4 $ tail $ LS.splitOn "," $ takeWhile (/= '[') inString | not (null inString)]


-- getPreArgString gets search srtategy fields (type, delta, min cost, CPUtime)
-- befroe arguments fields
getPreArgString ∷ String → [String]
getPreArgString inString =
    if null inString
        then []
        else take 4 $ tail $ LS.splitOn "," inString


{- | getBanditNames extracts the names of search bandits from comment list
already first part filtered out so only pairs in "(,)"
-}
getBanditNames ∷ String → [String]
getBanditNames inString =
    if null inString
        then []
        else
            let firstBanditName = takeWhile (/= ',') $ tail inString
                remainder = dropWhile (/= '(') $ tail inString
            in  firstBanditName : getBanditNames remainder


-- | getSearchIterations breaks up comment feild into individual iteration lines
getSearchIterations ∷ String → [String]
getSearchIterations inList =
    if null inList
        then []
        else
            let commentField = filter (/= '"') inList
                commentLines = LS.splitOn "]" commentField
            in  commentLines


-- changeDotPreamble takes an input string to search for and a new one to add in its place
-- searches through dot file (can have multipl graphs) replacing teh search string each time.
changeDotPreamble ∷ String → String → String → String
changeDotPreamble findString newString inDotString =
    if null inDotString
        then []
        else changePreamble' findString newString [] (lines inDotString)


-- changeDotPreamble' internal process for changeDotPreamble
changePreamble' ∷ String → String → [String] → [String] → String
changePreamble' findString newString accumList inLineList =
    if null inLineList
        then unlines $ reverse accumList
        else -- trace ("CP':" <> (head inLineList) <> " " <>  findString <> " " <> newString) (

            let firstLine = head inLineList
            in  if firstLine == findString
                    then changePreamble' findString newString (newString : accumList) (tail inLineList)
                    else changePreamble' findString newString (firstLine : accumList) (tail inLineList)


-- )

-- printGraph graphviz simple dot file of graph
-- execute with "dot -Teps test.dot -o test.eps"
-- need to add output to argument filename and call
-- graphviz via System.Process.runprocess
-- also, reorder GenForest so smalles (num leaves) is either first or
-- last so can print small to large all the way so easier to read
-- eps on OSX because ps gets cutt off for some reason and no pdf onOSX
-- -O foir multiple graphs I htink
printGraphVizDot ∷ String → String → PhyG ()
printGraphVizDot graphDotString dotFile =
    if null graphDotString
        then do failWithPhase Outputting "No graph to report"
        else do
            myHandle ← liftIO $ openFile dotFile WriteMode
            if os /= "darwin"
                then do logWith LogInfo ("\tOutputting graphviz to " <> dotFile <> ".pdf.\n")
                else do logWith LogInfo ("\tOutputting graphviz to " <> dotFile <> ".eps.\n")
            let outputType =
                    if os == "darwin"
                        then "-Teps"
                        else "-Tpdf"
            -- hPutStrLn myHandle "digraph G {"
            -- hPutStrLn myHandle "\trankdir = LR;"
            -- hPutStrLn myHandle "\tnode [ shape = rect];"
            -- hPutStr myHandle $ (unlines . tail . lines) graphDotString
            liftIO $ hPutStr myHandle graphDotString
            -- hPutStrLn myHandle "}"
            liftIO $ hClose myHandle
            pCode ← liftIO $ findExecutable "dot" -- system "dot" --check for Graphviz
            {-
            hPutStrLn stderr
                (if isJust pCode then --pCode /= Nothing then
                    "executed dot " <> outputType <> dotFile <> " -O " else
                    "Graphviz call failed (not installed or found).  Dot file still created. Dot can be obtained from https://graphviz.org/download")
            -}
            if isJust pCode
                then do
                    liftIO $ createProcess (proc "dot" [outputType, dotFile, "-O"])
                    logWith LogInfo ("\tExecuted dot " <> outputType <> " " <> dotFile <> " -O \n")
                else
                    logWith
                        LogInfo
                        "\tGraphviz call failed (not installed or found).  Dot file still created. Dot can be obtained from https://graphviz.org/download\n"


{- | showSearchFields creates a String list for SearchData Fields
special processing for notACommand that contains info for initial data and graph processing
-}
showSearchFields ∷ SearchData → [String]
showSearchFields sD =
    let inInstruction = instruction sD
        (instructionString, durationString, commentString') =
            if inInstruction /= NotACommand
                then (show $ instruction sD, show ((fromIntegral $ duration sD) / 1000 ∷ Double), commentString sD)
                else (commentString sD, show ((fromIntegral $ duration sD) / 1000000000000 ∷ Double), "No Comment")
    in  [ instructionString
        , unwords $ (showArg <$> arguments sD)
        , show $ minGraphCostIn sD
        , show $ maxGraphCostIn sD
        , show $ numGraphsIn sD
        , show $ minGraphCostOut sD
        , show $ maxGraphCostOut sD
        , show $ numGraphsOut sD
        , durationString
        , commentString'
        ]
    where
        showArg a = fst a <> ":" <> snd a


{- | requireReoptimization checks if set command in globalk settings requires reoptimization of graphs due to change in
graph type, optimality criterion etc.
-}
requireReoptimization ∷ GlobalSettings → GlobalSettings → Bool
requireReoptimization gsOld gsNew
    | graphType gsOld /= graphType gsNew = True
    | optimalityCriterion gsOld /= optimalityCriterion gsNew = True
    | finalAssignment gsOld /= finalAssignment gsNew = True
    | graphFactor gsOld /= graphFactor gsNew = True
    | rootCost gsOld /= rootCost gsNew = True
    | otherwise = False


-- | outputBlockTrees takes a PhyloGeneticTree and outputs BlockTrees
outputBlockTrees ∷ [String] → [VertexCost] → Int → (String, V.Vector [DecoratedGraph]) → String
outputBlockTrees commandList costList lOutgroupIndex (labelString, graphLV) =
    let blockIndexStringList =
            if ("dot" `elem` commandList) || ("dotpdf" `elem` commandList)
                then -- dot comments
                    fmap (((<> "\n") . ("//Block " <>)) . show) [0 .. ((V.length graphLV) - 1)]
                else -- newick
                    fmap (((<> "]\n") . ("[Block " <>)) . show) [0 .. ((V.length graphLV) - 1)]
        blockStrings = unlines (makeBlockGraphStrings commandList costList lOutgroupIndex <$> zip blockIndexStringList (V.toList graphLV))
    in  labelString <> blockStrings


-- | makeBlockGraphStrings makes individual block display trees--potentially multiple
makeBlockGraphStrings ∷ [String] → [VertexCost] → Int → (String, [DecoratedGraph]) → String
makeBlockGraphStrings commandList costList lOutgroupIndex (labelString, graphL) =
    let diplayIndexString =
            if ("dot" `elem` commandList) || ("dotpdf" `elem` commandList)
                then ("//Display Tree(s): " <> show (length graphL) <> "\n")
                else ("[Display Tree[s]: " <> show (length graphL) <> "]\n")
        displayString = (<> "\n") $ outputDisplayString commandList costList lOutgroupIndex graphL
    in  labelString <> diplayIndexString <> displayString


-- | outputDisplayString is a wrapper around graph output functions--but without cost list
outputDisplayString ∷ [String] → [VertexCost] → Int → [DecoratedGraph] → String
outputDisplayString commandList costList lOutgroupIndex graphList
    | "dot" `elem` commandList =
        makeDotList
            (not (elem "nobranchlengths" commandList))
            (not (elem "nohtulabels" commandList))
            costList
            lOutgroupIndex
            (fmap GO.convertDecoratedToSimpleGraph graphList)
    | "newick" `elem` commandList =
        GO.makeNewickList
            False
            (not (elem "nobranchlengths" commandList))
            (not (elem "nohtulabels" commandList))
            lOutgroupIndex
            (fmap GO.convertDecoratedToSimpleGraph graphList)
            (replicate (length graphList) 0.0)
    | "ascii" `elem` commandList = makeAsciiList lOutgroupIndex (fmap GO.convertDecoratedToSimpleGraph graphList)
    | otherwise -- "dot" as default
        =
        makeDotList
            (not (elem "nobranchlengths" commandList))
            (not (elem "nohtulabels" commandList))
            costList
            lOutgroupIndex
            (fmap GO.convertDecoratedToSimpleGraph graphList)


-- | outputGraphString is a wrapper around graph output functions
outputGraphString ∷ [String] → Int → [DecoratedGraph] → [VertexCost] → String
outputGraphString commandList lOutgroupIndex graphList costList
    | "dot" `elem` commandList =
        makeDotList
            (not (elem "nobranchlengths" commandList))
            (not (elem "nohtulabels" commandList))
            costList
            lOutgroupIndex
            (fmap GO.convertDecoratedToSimpleGraph graphList)
    | "newick" `elem` commandList =
        GO.makeNewickList
            False
            (not (elem "nobranchlengths" commandList))
            (not (elem "nohtulabels" commandList))
            lOutgroupIndex
            (fmap GO.convertDecoratedToSimpleGraph graphList)
            costList
    | "ascii" `elem` commandList = makeAsciiList lOutgroupIndex (fmap GO.convertDecoratedToSimpleGraph graphList)
    | otherwise -- "dot" as default
        =
        makeDotList
            (not (elem "nobranchlengths" commandList))
            (not (elem "nohtulabels" commandList))
            costList
            lOutgroupIndex
            (fmap GO.convertDecoratedToSimpleGraph graphList)


-- | outputGraphStringSimple is a wrapper around graph output functions
outputGraphStringSimple ∷ [String] → Int → [SimpleGraph] → [VertexCost] → String
outputGraphStringSimple commandList lOutgroupIndex graphList costList
    | "dot" `elem` commandList =
        makeDotList (not (elem "nobranchlengths" commandList)) (not (elem "nohtulabels" commandList)) costList lOutgroupIndex graphList
    | "newick" `elem` commandList = GO.makeNewickList False True True lOutgroupIndex graphList costList
    | "ascii" `elem` commandList = makeAsciiList lOutgroupIndex graphList
    | otherwise -- "dot" as default
        =
        makeDotList (not (elem "nobranchlengths" commandList)) (not (elem "nohtulabels" commandList)) costList lOutgroupIndex graphList


{- | makeDotList takes a list of fgl trees and outputs a single String cointaining the graphs in Dot format
need to specify -O option for multiple graph(outgroupIndex globalSettings)s
-}
makeDotList ∷ Bool → Bool → [VertexCost] → Int → [SimpleGraph] → String
makeDotList writeEdgeWeight writeNodeLabel costList rootIndex graphList =
    let graphStringList' = fmap (fgl2DotString . LG.rerootTree rootIndex) graphList
        graphStringList = fmap (stripDotLabels writeEdgeWeight writeNodeLabel) graphStringList'
        costStringList = fmap (("\n//" <>) . show) costList
    in  L.intercalate "\n" (zipWith (<>) graphStringList costStringList)


-- | stripDotLabels strips away edge and HTU labels from dot files
stripDotLabels ∷ Bool → Bool → String → String
stripDotLabels writeEdgeWeight writeNodeLabel inGraphString =
    if null inGraphString
        then inGraphString
        else
            if writeNodeLabel && writeEdgeWeight
                then inGraphString
                else
                    let firstString =
                            if writeEdgeWeight
                                then inGraphString
                                else stripEdge inGraphString

                        secondString =
                            if writeNodeLabel
                                then firstString
                                else stripNode firstString
                    in  secondString


-- | stripEdge removes edge labels from HTUs in graphviz format string
stripEdge ∷ String → String
stripEdge inString =
    if null inString
        then inString
        else
            let lineStringList = lines inString
                newLines = fmap makeNewLine lineStringList
            in  unlines newLines
    where
        makeNewLine a =
            if (null $ L.intersect "->" a)
                then a
                else
                    let b = words a
                    in  "    " <> (concat [b !! 0, " ", b !! 1, " ", b !! 2]) <> " [];"


-- | stripNode removes edge labels from HTUs in graphviz format string
stripNode ∷ String → String
stripNode inString =
    if null inString
        then inString
        else
            let lineStringList = lines inString
                newLines = fmap makeNewLine lineStringList
            in  unlines newLines
    where
        makeNewLine a =
            if (null $ L.intersect "HTU" a)
                then a
                else
                    let b = words a
                        c = take 10 $ b !! 1
                        newLine =
                            if c == "[label=HTU"
                                then "    " <> (head b) <> " [];"
                                else a
                    in  newLine


-- | makeAsciiList takes a list of fgl trees and outputs a single String cointaining the graphs in ascii format
makeAsciiList ∷ Int → [SimpleGraph] → String
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

{- | getDataListList returns a list of lists of Strings for data output as csv
for row is source file names, subsequent rows by taxon with +/- for present absent taxon in
input file
different from getDataListList in removeal or processed data requiremenrt replaced with taxan name list
-}
getDataListList ∷ (Foldable f) ⇒ [RawData] → f T.Text → [[String]]
getDataListList inDataList comprehensiveTaxaSet
    | null inDataList = []
    | otherwise =
        let taxaList = toList comprehensiveTaxaSet
            fileNames = " " : fmap (takeWhile (/= ':') . T.unpack) (fmap name $ fmap head $ fmap snd inDataList)
            presenceAbsenceList = fmap (isThere inDataList) taxaList
            fullMatrix = zipWith (:) (fmap T.unpack taxaList) presenceAbsenceList
        in  -- trace (show fileNames)
            fileNames : fullMatrix


-- | isThere takes a list of Rawdata and reurns a String of + -
isThere ∷ [RawData] → T.Text → [String]
isThere inData inName =
    if null inData
        then []
        else
            let firstTaxList = fmap fst $ fst $ head inData
            in  if inName `elem` firstTaxList
                    then "+" : isThere (tail inData) inName
                    else "-" : isThere (tail inData) inName


{- | phyloDataToString converts RawData type to String
for additive chars--multiply states by weight is < 1 when outputtting due to conversion on input
-}
phyloDataToString ∷ Int → V.Vector BlockData → [[String]]
phyloDataToString charIndexStart inDataVect =
    if V.null inDataVect
        then []
        else
            let (blockName, blockData, charInfoVect) = V.head inDataVect
                charStrings = zipWith (:) (replicate (V.length charInfoVect) (T.unpack blockName)) (getCharInfoStrings <$> V.toList charInfoVect)
                charNumberString = fmap show [charIndexStart .. (charIndexStart + length charStrings - 1)]
                fullMatrix = zipWith (:) charNumberString charStrings
            in  fullMatrix <> phyloDataToString (charIndexStart + length charStrings) (V.tail inDataVect)


-- | extractAlphabetStrings takes vector of block data and returns block x character string list
extractAlphabetStrings ∷ V.Vector BlockData → [[[String]]]
extractAlphabetStrings inBlockDataV =
    let result = V.toList $ fmap extractAlphabetStringsBlock inBlockDataV
    in  trace
            ("EAS: " <> (show result))
            result


-- | extractAlphabetStringsBlock takes block data and returns character string list
extractAlphabetStringsBlock ∷ BlockData → [[String]]
extractAlphabetStringsBlock inBlockData = V.toList $ fmap extractAlphabetStringsChar (thd3 inBlockData)


-- | extractAlphabetStringsChar takes char data and returns character string list
extractAlphabetStringsChar ∷ CharInfo → [String]
extractAlphabetStringsChar inCharInfo = (alphabetSymbols $ ST.toString <$> alphabet inCharInfo)


{- | getDataElementTransformations takes alphabet, parent and child final states
and calculates and formats the transition matrix in frequency and raw numbers
this for list of blocks each with list of characters
over all edges
need to ddeal woith missing data
-}
getDataElementTransformations ∷ [[[String]]] → [[String]] → [[String]] → [[String]] → [[String]]
getDataElementTransformations alphabetStrings parentLL childLL parentChildLL =
    if null childLL
        then []
        else -- convert to block x character (Parent String, Child String) for counting exercise

            let numBlocks = length alphabetStrings
                numCharsList = fmap length alphabetStrings

                numParent = length parentLL
                numChild = length childLL
            in  -- dataByBlock = reorderBySection numBlocks diffLL []
                -- dataByBlockChar = zipeWith reorderBySection numCharsList dataByBlock

                trace
                    ("GDET: " <> (show (numBlocks, numCharsList, numParent, numChild)) <> " " <> (show alphabetStrings)) -- <>
                    -- (show parentLL) <> "\n\n" <> (show childLL)) $
                    parentChildLL


-- | reorderBySection takes a list of list of Strings and reorders klists into first diominant (as in blocks here)
reorderBySection ∷ Int → [[String]] → [[String]] → [[String]]
reorderBySection numSections inData reorgData =
    if null inData
        then reorgData
        else
            if length inData < numSections
                then error ("Incorrect input list size: " <> (show numSections) <> " sections and " <> (show $ length inData) <> " pieces")
                else
                    let firstGroup = filter (not . null) $ take numSections inData
                        newData =
                            if null reorgData
                                then firstGroup
                                else zipWith (<>) reorgData firstGroup
                    in  trace ("RBS: " <> (show firstGroup)) $
                            reorderBySection numSections (drop numSections inData) newData


{- | getBlockElementTransformations
need to deal with missing data
-}
getBlockElementTransformations ∷ [[String]] → [String] → [String]
getBlockElementTransformations alphabetStringL diffL =
    if null diffL
        then []
        else
            trace ("GBET: " <> (show alphabetStringL) <> " " <> (concat diffL)) $
                zipWith getCharElementTransformations alphabetStringL diffL


{- | getCharElementTransformations
need to deal with missing data
-}
getCharElementTransformations ∷ [String] → String → String
getCharElementTransformations alphabetString diffState =
    trace
        ("GCET: " <> (show alphabetString) <> " " <> diffState)
        []


-- | getDataElementFrequencies takes a vecor of block data and returns character element frequencies
-- as tiple (in String) of element, frequncy and number of elements
-- 
getDataElementFrequencies ∷ Bool → V.Vector BlockData → [[String]]
getDataElementFrequencies useIA inBlockDataV =
    if V.null inBlockDataV
        then []
        else
            let blockCharFreqLLL = V.toList $ V.zipWith (getBlockCharElementFrequencies useIA) (fmap snd3 inBlockDataV) (fmap thd3 inBlockDataV)
            in  
            -- this for triple info (element, frequency, and number)
            -- fmap (fmap show) blockCharFreqLLL
            -- this for frequency only
            fmap (fmap show) $ fmap (fmap (fmap snd3)) blockCharFreqLLL


{- | getBlockElementFrequencies gets the element grequencies for each character in a block
if an unaligned sequence type infers based on length differences of inputs using number of gaps to make
inputs square (so minimum number of gaps)
-}
getBlockCharElementFrequencies ∷ Bool → V.Vector (V.Vector CharacterData) → V.Vector CharInfo → [[(String, Double, Int)]]
getBlockCharElementFrequencies useIA charDataV charInfoV =
    if V.null charDataV
        then []
        else -- set to preliminary states or IA states

            let dataV =
                    if not useIA
                        then PRE.setFinalToPreliminaryStates charDataV
                        else charDataV
            in  V.toList $ V.zipWith (getCharElementFrequencies useIA) (U.transposeVector dataV) charInfoV


{- | getCharElementFrequencies gets element frequencies for a character
if an unaligned sequence type infers based on length differences of inputs using number of gaps to make
inputs square (so minimum number of gaps)
returns frequencies and number for each alphabet element
ignores ambiguities/polymorphism
-}
getCharElementFrequencies ∷ Bool → V.Vector CharacterData → CharInfo → [(String, Double, Int)]
getCharElementFrequencies useIA charData charInfo =
    if V.null charData
        then []
        else
            let -- alphabet element strings
                alphabetElementStrings = (alphabetSymbols $ ST.toString <$> alphabet charInfo)

                charInfoV = V.replicate (V.length charData) charInfo
                charPairV = V.zip charData charInfoV
                taxonElementList = V.toList $ fmap L.last $ fmap (makeCharLine useIA) charPairV

                -- get implicit gaps from unaligned seqs (will end up zero if all equal in length)
                numExtraGaps =
                    if "-" `notElem` alphabetElementStrings
                        then 0
                        else
                            let maxLength = maximum $ fmap length taxonElementList
                            in  getMinGapNumber maxLength 0 taxonElementList

                totalElementList = concat $ (L.replicate numExtraGaps '-') : taxonElementList
                elementsGroups = L.group $ L.sort totalElementList
                elementList = fmap (: []) $ fmap head elementsGroups
                doubleList = zip elementList (fmap length elementsGroups)
                elementNumberList = reorderAsInAlphabet alphabetElementStrings doubleList
                numElements = fromIntegral $ sum elementNumberList
                elementfreqList = fmap (/ numElements) $ fmap fromIntegral elementNumberList
            in  -- trace ("GCEF: " <> totalElementList ) $ -- <> " -> " <> (concat $ concat $ fmap (makeCharLine usaIA) charPairV))
                zip3 alphabetElementStrings elementfreqList elementNumberList


-- | getMinGapNumber gets implicit gap number by summing length differneces among sequences in order
getMinGapNumber ∷ Int → Int → [String] → Int
getMinGapNumber maxLength curNum inElementLL =
    if null inElementLL
        then curNum
        else
            let increment = maxLength - (length $ head inElementLL)
            in  getMinGapNumber maxLength (curNum + increment) (tail inElementLL)


-- | reorderAsInAlphabet
reorderAsInAlphabet ∷ [String] → [(String, Int)] → [Int]
reorderAsInAlphabet inAlphList inDoubleList =
    -- trace ("RASIA: " <> (show inAlphList) <> " " <> (show inDoubleList)) $
    if null inDoubleList
        then []
        else fmap (findDouble inDoubleList) inAlphList


-- | findDouble takes list of strings then pulls the tipple wthat matches the string
findDouble ∷ [(String, Int)] → String → Int
findDouble inDoubleList matchString =
    if null inDoubleList
        then 0
        else
            let (alphElement, number) = head inDoubleList
            in  if alphElement == matchString
                    then number
                    else findDouble (tail inDoubleList) matchString


-- | getCharInfoStrings takes charInfo and returns list of Strings of fields
getCharInfoStrings ∷ CharInfo → [String]
getCharInfoStrings inChar =
    let activityString =
            if activity inChar
                then "active"
                else "inactive"
        prealignedString =
            if prealigned inChar
                then "prealigned"
                else "unaligned"
    in  [T.unpack $ name inChar, show $ charType inChar, activityString, show $ weight inChar, prealignedString]
            <> [init $ concat $ fmap (<> ",") $ fmap ST.toString . toList $ alphabet inChar]
            <> [show $ costMatrix inChar]


{- | executeRenameReblockCommands takes all the "Rename commands" pairs and
creates a list of pairs of new name and list of old names to be converted
as Text
-}
executeRenameReblockCommands ∷ Instruction → [(T.Text, T.Text)] → [Command] → IO [(T.Text, T.Text)]
executeRenameReblockCommands thisInStruction curPairs commandList =
    if null commandList
        then return curPairs
        else do
            let (firstOption, firstArgs) = head commandList

            -- skip "Read" and "Rename "commands already processed
            if firstOption /= thisInStruction
                then executeRenameReblockCommands thisInStruction curPairs (tail commandList)
                else
                    let newName = T.filter C.isPrint $ T.filter (/= '"') $ T.pack $ snd $ head firstArgs
                        newNameList = replicate (length $ tail firstArgs) newName
                        oldNameList = fmap (T.filter (/= '"') . T.pack) (fmap snd $ tail firstArgs)
                        newPairs = zip newNameList oldNameList
                        changeNamePairs = filter areDiff (curPairs <> newPairs)
                        newChangeNames = fmap fst changeNamePairs
                        origChangeNames = fmap snd changeNamePairs
                        intersectionChangeNames = L.intersect newChangeNames origChangeNames
                    in  if (not $ null intersectionChangeNames)
                            then errorWithoutStackTrace ("Renaming of " <> (show intersectionChangeNames) <> " as both a new name and to be renamed")
                            else executeRenameReblockCommands thisInStruction (curPairs <> newPairs) (tail commandList)
    where
        areDiff (a, b) =
            if a /= b
                then True
                else False

{- | getGraphParameters estimates basic parameters--element frequencies and transformation frequencies
    From uniquely optimized (ie not 'R') edges and veritces
-}
getGraphParameters ∷ GlobalSettings → ProcessedData → (ReducedPhylogeneticGraph, Int) → PhyG [[String]]
getGraphParameters _ inData (inGraph, graphIndex) =
    let decGraph = thd5 inGraph
    in  
    if LG.isEmpty decGraph then pure []
    else do
                    let useIA = True
                    let useDO = False
                    let vertexList = LG.labNodes decGraph
                    let edgeList = LG.labEdges decGraph
                    
                    vertexInfo <-
                            let action :: LG.LNode VertexInfo → [[String]]
                                action = getVertexCharInfo useDO (thd3 inData) (fst5 inGraph) (fft5 inGraph) 
                            in do
                                diagPar ← getParallelChunkMap
                                let result = diagPar action vertexList
                                pure result 
                  
                    -- False for Use IA here
                    
                    let alphabetInfo = getDataElementFrequencies False (thd3 inData)

                    -- Get changes based on edges
                    let edgeIndexList = fmap LG.toEdge edgeList

                    -- should be parallel
                    let vertexVect = V.fromList $ vertexList
                    edgeTransformationTotalList <- 
                            let action :: (LG.Node, LG.Node) → [[(String, [(String, Int)], [[Int]])]]
                                action = getEdgeTransformations vertexVect (fft5 inGraph)
                            in do
                                actionPar <- getParallelChunkMap
                                let result = actionPar action edgeIndexList
                                pure result

                    -- extract relevent info
                    let edgeTransformationList = fmap (fmap (fmap fst3)) edgeTransformationTotalList
                    let transformationNumbersLLL = fmap (fmap (fmap thd3)) edgeTransformationTotalList

                    let overallElementTransformations = sumEdgeTransformationLists [] transformationNumbersLLL
                    overallElementTransformationsFreq <-
                            let action :: [[[Int]]] → [[[Double]]]
                                action = fmap U.normalizeMatrix
                            in do
                                actionPar <- getParallelChunkMap
                                let result = actionPar action overallElementTransformations
                                pure result

                    let vertexHeader = fmap (fmap (take 9)) vertexInfo

                    let edgeListLists = knitTitlesChangeInfo vertexHeader edgeTransformationList

                    let transformationHeader = fmap (drop 5) $ tail $ head vertexHeader

                    let statsListList =
                            addBlockCharStrings
                                transformationHeader
                                alphabetInfo
                                (fmap (fmap show) overallElementTransformationsFreq)
                                (fmap (fmap show) overallElementTransformations)

                    let edgeTransformationTitle =
                            [ [" "]
                            , ["Edge Transformation Statistics"]
                            ,
                                [ ""
                                , "Data Block"
                                , "Character Name"
                                , "Character Type"
                                , "Element Frequencies"
                                , "Transformation Frequencies"
                                ]
                            ]
                    pure $  edgeTransformationTitle
                            <> statsListList
                    where
                        -- <> fmap show overallElementTransformations

                        concat4 ∷ ∀ {a}. (Semigroup a) ⇒ a → a → a → a → a
                        concat4 a b c d = a <> b <> c <> d
                        concat3 ∷ ∀ {a}. (Semigroup a) ⇒ a → a → a → a
                        concat3 a b c = a <> b <> c



{- | getGraphMetaData creates basic strings for CSV of graph vertex and node information
    nodes first then vertices
    Now--spawned from getDiagnosis, but assuming they will diverge and want metData to 
    be stable to visualization codes
-}
getGraphMetaData ∷ GlobalSettings → ProcessedData → (ReducedPhylogeneticGraph, Int) → PhyG [[String]]
getGraphMetaData _ inData (inGraph, graphIndex) =
    let decGraph = thd5 inGraph
    in  
    if LG.isEmpty decGraph then pure []
    else do
                    let useIA = True
                    let useDO = False
                    let vertexList = LG.labNodes decGraph
                    let edgeList = LG.labEdges decGraph
                    let numLeaves = V.length $ fst3 inData

                        -- Vertex final character states for currect node
                    let vertexTitle = ["Vertex Character State Information"]
                        -- topHeaderList  = ["Graph Index", "Vertex Index", "Vertex Name", "Vertex Type", "Child Vertices", "Parent Vertices", "Data Block", "Character Name", "Character Type", "Preliminary State", "Final State", "Local Cost"]
                    let topHeaderList =
                                [ "Graph Index"
                                , "Vertex Index"
                                , "Vertex Name"
                                , "Vertex Type"
                                , "Child Vertices"
                                , "Parent Vertices"
                                , "Data Block"
                                , "Character Name"
                                , "Character Type"
                                , "Final State"
                                ]

                    -- let vertexInfo = fmap (getVertexCharInfo useDO (thd3 inData) (fst5 inGraph) (fft5 inGraph)) vertexList
                    vertexInfo <-
                            let action :: LG.LNode VertexInfo → [[String]]
                                action = getVertexCharInfo useDO (thd3 inData) (fst5 inGraph) (fft5 inGraph) 
                            in do
                                diagPar ← getParallelChunkMap
                                let result = diagPar action vertexList
                                pure result 
                    let vertexInfoList = concat vertexInfo

                    -- this using IA fields to get changes
                    let vertexInfoListChanges = vertexInfoList -- concatMap (getVertexCharInfo useIA (thd3 inData) (fst5 inGraph) (fft5 inGraph)) vertexList

                    -- Edge length information
                    let edgeTitle = [[" "], ["Edge Weight/Length Information"]]
                    let edgeHeaderList = [[" ", "Edge Head Vertex", "Edge Tail Vertex", "Edge Type", "Minimum Length", "Maximum Length", "MidRange Length"]]
                    
                    -- this is basically a Show
                    let edgeInfoList = fmap U.getEdgeInfo edgeList

                    -- Alphabet element numbers
                    let alphabetTitle = [["Alphabet (element, frequency, number) Gap, if estimated from unaligned sequences, is a minimum"]]
                    -- False for Use IA here
                    let alphabetInfo = getDataElementFrequencies False (thd3 inData)

                    let alphbetStringLL = extractAlphabetStrings (thd3 inData)

                    let vertexChangeTitle =
                            [ [" "]
                            , ["Vertex Character Changes"]
                            ,
                                [ "Graph Index"
                                , "Vertex Index"
                                , "Vertex Name"
                                , "Vertex Type"
                                , "Child Vertices"
                                , "Parent Vertices"
                                , "Data Block"
                                , "Character Name"
                                , "Character Type"
                                , "Parent Final State"
                                , "Node Final State"
                                -- , "Sequence Changes (position, parent final state, node final state)"
                                ]
                            ]

                    -- let vertexParentStateList = fmap ((: []) . last) (concatMap (getVertexAndParentCharInfo useDO (thd3 inData) (fst5 inGraph) (fft5 inGraph) (V.fromList vertexList)) vertexList)
                    vertParCharInfoList <-
                            let action :: LG.LNode VertexInfo → [[String]]
                                action = getVertexAndParentCharInfo useDO (thd3 inData) (fst5 inGraph) (fft5 inGraph) (V.fromList vertexList)
                            in do
                                actionPar <- getParallelChunkMap
                                let result = actionPar action vertexList
                                pure result

                    let vertexParentStateList = fmap ((: []) . last) $ concat vertParCharInfoList
                    let vertexStateList = fmap (drop 9) vertexInfoListChanges

                    -- process to change to lines of individual changes--basically a transpose
                    -- True to only report diffs
                    -- vertexChangeListByPosition = fmap (getAlignmentBasedChanges' True 0) (zip vertexParentStateList vertexStateList)
                    -- let parentChildStatesList = fmap (getAlignmentBasedChanges' False 0) (zip vertexParentStateList vertexStateList)
                    parentChildStatesList <- 
                            let action :: ([String], [String]) → [String]
                                action = getAlignmentBasedChanges' False 0
                            in do
                                actionPar <- getParallelChunkMap
                                let result = actionPar action (zip vertexParentStateList vertexStateList)
                                pure result

                    -- putting parent states before current state
                    let vertexStateInfoList = fmap (take 9) vertexInfoListChanges

                    let vertexChangeList = L.zipWith3 concat3 vertexStateInfoList vertexParentStateList vertexStateList -- 
                    -- filter out those that are the same states
                    let differenceList = removeNoChangeLines vertexChangeList

                    -- element transformation numbers
                    let elementTransformationTitle = [["Element Transformations (element<->element, frequency, number) based in Implied Alignment for unaligned sequences"]]
                    -- get element transfomation by re-parsing formated results--uses teh character states strings this way
                    let elementTransformationInfo = getDataElementTransformations alphbetStringLL vertexParentStateList vertexStateList parentChildStatesList

                    -- Get changes based on edges
                    let edgeIndexList = fmap LG.toEdge edgeList

                    -- should be parallel
                    let vertexVect = V.fromList $ vertexList
                    -- let edgeTransformationTotalList = fmap (getEdgeTransformations vertexVect (fft5 inGraph)) edgeIndexList
                    edgeTransformationTotalList <- 
                            let action :: (LG.Node, LG.Node) → [[(String, [(String, Int)], [[Int]])]]
                                action = getEdgeTransformations vertexVect (fft5 inGraph)
                            in do
                                actionPar <- getParallelChunkMap
                                let result = actionPar action edgeIndexList
                                pure result

                    -- extract relevent info
                    let edgeTransformationList = fmap (fmap (fmap fst3)) edgeTransformationTotalList
                    -- elementNumbersLLL = fmap (fmap (fmap snd3)) edgeTransformationTotalList
                    let transformationNumbersLLL = fmap (fmap (fmap thd3)) edgeTransformationTotalList

                    -- overallElementNumbers = sumEdgeElementLists [] $ take numLeaves elementNumbersLLL
                    let overallElementTransformations = sumEdgeTransformationLists [] transformationNumbersLLL
                    -- let overallElementTransformationsFreq = fmap (fmap U.normalizeMatrix) overallElementTransformations
                    overallElementTransformationsFreq <-
                            let action :: [[[Int]]] → [[[Double]]]
                                action = fmap U.normalizeMatrix
                            in do
                                actionPar <- getParallelChunkMap
                                let result = actionPar action overallElementTransformations
                                pure result

                    let vertexHeader = fmap (fmap (take 9)) vertexInfo

                    let edgeListLists = knitTitlesChangeInfo vertexHeader edgeTransformationList

                    let transformationHeader = fmap (drop 5) $ tail $ head vertexHeader

                    
                    let vertexChangeTitleNew =
                            [ [" "]
                            , ["Vertex Character Changes"]
                            ,
                                [ ""
                                , "Vertex Index"
                                , "Vertex Name"
                                , "Vertex Type"
                                , "Child Vertices"
                                , "Parent Vertices"
                                , "Data Block"
                                , "Character Name"
                                , "Character Type"
                                , "Parent Final State"
                                , "Node Final State"
                                , "Unambiguous Transformations"
                                ]
                            ]

                    pure $  [vertexTitle, topHeaderList, [show graphIndex]]
                            <> vertexInfoList
                            <> edgeTitle
                            <> edgeHeaderList
                            <> edgeInfoList
                            -- <> vertexChangeTitle
                            -- <> differenceList
                            <> vertexChangeTitleNew
                            <> edgeListLists
                            -- <> elementTransformationTitle
                            -- <> elementTransformationInfo
                            -- <> alphabetTitle
                            -- <> alphabetInfo
                    where
                        -- <> fmap show overallElementTransformations

                        concat4 ∷ ∀ {a}. (Semigroup a) ⇒ a → a → a → a → a
                        concat4 a b c d = a <> b <> c <> d
                        concat3 ∷ ∀ {a}. (Semigroup a) ⇒ a → a → a → a
                        concat3 a b c = a <> b <> c

{- | getGraphDiagnosis creates basic strings for CSV of graph vertex and node information
nodes first then vertices
-}
getGraphDiagnosis ∷ GlobalSettings → ProcessedData → (ReducedPhylogeneticGraph, Int) → PhyG [[String]]
getGraphDiagnosis inGS inData (inGraph, graphIndex) =
    let decGraph = thd5 inGraph
    in  
    if LG.isEmpty decGraph then pure []
    else do
                    let useIA = True
                    let useDO = False
                    let vertexList = LG.labNodes decGraph
                    let edgeList = LG.labEdges decGraph
                    let numLeaves = V.length $ fst3 inData

                        -- Vertex final character states for currect node
                    let vertexTitle = ["Vertex Character State Information"]
                        -- topHeaderList  = ["Graph Index", "Vertex Index", "Vertex Name", "Vertex Type", "Child Vertices", "Parent Vertices", "Data Block", "Character Name", "Character Type", "Preliminary State", "Final State", "Local Cost"]
                    let topHeaderList =
                                [ "Graph Index"
                                , "Vertex Index"
                                , "Vertex Name"
                                , "Vertex Type"
                                , "Child Vertices"
                                , "Parent Vertices"
                                , "Data Block"
                                , "Character Name"
                                , "Character Type"
                                , "Final State"
                                ]

                    -- let vertexInfo = fmap (getVertexCharInfo useDO (thd3 inData) (fst5 inGraph) (fft5 inGraph)) vertexList
                    vertexInfo <-
                            let action :: LG.LNode VertexInfo → [[String]]
                                action = getVertexCharInfo useDO (thd3 inData) (fst5 inGraph) (fft5 inGraph) 
                            in do
                                diagPar ← getParallelChunkMap
                                let result = diagPar action vertexList
                                pure result 
                    let vertexInfoList = concat vertexInfo

                    -- this using IA fields to get changes
                    let vertexInfoListChanges = vertexInfoList -- concatMap (getVertexCharInfo useIA (thd3 inData) (fst5 inGraph) (fft5 inGraph)) vertexList

                    
                    let edgeInfoList = fmap U.getEdgeInfo edgeList

                    -- Get complexity information--empty of not complexity
                    edgeTripleList <- U.getEdgeComplexityFactors inGS inData vertexList edgeList 
                    let edgeComplexityFactor = fmap thd3 edgeTripleList

                    let (vertexComplexityLabel, vertexComplexityList) = if (optimalityCriterion inGS `elem` [PMDL, SI]) then
                                                                ("Complexity Factor", fmap (:[]) $ fmap show edgeComplexityFactor)
                                               else ("", fmap (:[]) $ replicate (length edgeList) "")

                    -- Edge length information
                    let edgeTitle = [[" "], ["Edge Weight/Length Information"]]
                    let edgeHeaderList = [[" ", "Edge Head Vertex", "Edge Tail Vertex", "Edge Type", "Minimum Length", "Maximum Length", "MidRange Length", vertexComplexityLabel]]

                    -- Alphabet element numbers
                    let alphabetTitle = [["Alphabet (element, frequency, number) Gap, if estimated from unaligned sequences, is a minimum"]]
                    -- False for Use IA here
                    let alphabetInfo = getDataElementFrequencies False (thd3 inData)

                    let alphbetStringLL = extractAlphabetStrings (thd3 inData)

                    let vertexChangeTitle =
                            [ [" "]
                            , ["Vertex Character Changes"]
                            ,
                                [ "Graph Index"
                                , "Vertex Index"
                                , "Vertex Name"
                                , "Vertex Type"
                                , "Child Vertices"
                                , "Parent Vertices"
                                , "Data Block"
                                , "Character Name"
                                , "Character Type"
                                , "Parent Final State"
                                , "Node Final State"
                                -- , "Sequence Changes (position, parent final state, node final state)"
                                ]
                            ]

                    -- let vertexParentStateList = fmap ((: []) . last) (concatMap (getVertexAndParentCharInfo useDO (thd3 inData) (fst5 inGraph) (fft5 inGraph) (V.fromList vertexList)) vertexList)
                    vertParCharInfoList <-
                            let action :: LG.LNode VertexInfo → [[String]]
                                action = getVertexAndParentCharInfo useDO (thd3 inData) (fst5 inGraph) (fft5 inGraph) (V.fromList vertexList)
                            in do
                                actionPar <- getParallelChunkMap
                                let result = actionPar action vertexList
                                pure result

                    let vertexParentStateList = fmap ((: []) . last) $ concat vertParCharInfoList
                    let vertexStateList = fmap (drop 9) vertexInfoListChanges

                    -- process to change to lines of individual changes--basically a transpose
                    -- True to only report diffs
                    -- vertexChangeListByPosition = fmap (getAlignmentBasedChanges' True 0) (zip vertexParentStateList vertexStateList)
                    -- let parentChildStatesList = fmap (getAlignmentBasedChanges' False 0) (zip vertexParentStateList vertexStateList)
                    parentChildStatesList <- 
                            let action :: ([String], [String]) → [String]
                                action = getAlignmentBasedChanges' False 0
                            in do
                                actionPar <- getParallelChunkMap
                                let result = actionPar action (zip vertexParentStateList vertexStateList)
                                pure result

                    -- putting parent states before current state
                    let vertexStateInfoList = fmap (take 9) vertexInfoListChanges

                    let vertexChangeList = L.zipWith3 concat3 vertexStateInfoList vertexParentStateList vertexStateList -- 
                    -- filter out those that are the same states
                    let differenceList = removeNoChangeLines vertexChangeList

                    -- element transformation numbers
                    let elementTransformationTitle = [["Element Transformations (element<->element, frequency, number) based in Implied Alignment for unaligned sequences"]]
                    -- get element transfomation by re-parsing formated results--uses teh character states strings this way
                    let elementTransformationInfo = getDataElementTransformations alphbetStringLL vertexParentStateList vertexStateList parentChildStatesList

                    -- Get changes based on edges
                    let edgeIndexList = fmap LG.toEdge edgeList

                    -- should be parallel
                    let vertexVect = V.fromList $ vertexList
                    -- let edgeTransformationTotalList = fmap (getEdgeTransformations vertexVect (fft5 inGraph)) edgeIndexList
                    edgeTransformationTotalList <- 
                            let action :: (LG.Node, LG.Node) → [[(String, [(String, Int)], [[Int]])]]
                                action = getEdgeTransformations vertexVect (fft5 inGraph)
                            in do
                                actionPar <- getParallelChunkMap
                                let result = actionPar action edgeIndexList
                                pure result

                    -- extract relevent info
                    let edgeTransformationList = fmap (fmap (fmap fst3)) edgeTransformationTotalList
                    -- elementNumbersLLL = fmap (fmap (fmap snd3)) edgeTransformationTotalList
                    let transformationNumbersLLL = fmap (fmap (fmap thd3)) edgeTransformationTotalList

                    -- overallElementNumbers = sumEdgeElementLists [] $ take numLeaves elementNumbersLLL
                    let overallElementTransformations = sumEdgeTransformationLists [] transformationNumbersLLL
                    -- let overallElementTransformationsFreq = fmap (fmap U.normalizeMatrix) overallElementTransformations
                    overallElementTransformationsFreq <-
                            let action :: [[[Int]]] → [[[Double]]]
                                action = fmap U.normalizeMatrix
                            in do
                                actionPar <- getParallelChunkMap
                                let result = actionPar action overallElementTransformations
                                pure result

                    let vertexHeader = fmap (fmap (take 9)) vertexInfo

                    let edgeListLists = knitTitlesChangeInfo vertexHeader edgeTransformationList

                    -- test code to get conditional info of child given parent K(c|p) = K(c<>p) - K(p), K(c)/K(c|p) randomnes index 
                    let characterParentFields = fmap (getIndex 9) edgeListLists
                    let characterChildFields = fmap (getIndex 10) edgeListLists
                    let characterEdgeInformationContent = zipWith getEdgeInformationContent characterChildFields characterParentFields

                    let transformationHeader = fmap (drop 5) $ tail $ head vertexHeader

                    
                    let vertexChangeTitleNew =
                            [ [" "]
                            , ["Vertex Character Changes"]
                            ,
                                [ ""
                                , "Vertex Index"
                                , "Vertex Name"
                                , "Vertex Type"
                                , "Child Vertices"
                                , "Parent Vertices"
                                , "Data Block"
                                , "Character Name"
                                , "Character Type"
                                , "Parent Final State"
                                , "Node Final State"
                                , "Unambiguous Transformations"
                                ]
                            ]

                    pure $  [vertexTitle, topHeaderList, [show graphIndex]]
                            <> vertexInfoList
                            <> edgeTitle
                            <> edgeHeaderList
                            <> zipWith (<>) edgeInfoList vertexComplexityList
                            -- <> vertexChangeTitle
                            -- <> differenceList
                            <> vertexChangeTitleNew
                            <> edgeListLists
                            --compression lists 
                            -- <> vertexComplexityList
                            -- <> zipWith (<>) edgeListLists (fmap (:[]) $ fmap showInfo characterEdgeInformationContent)
                            -- <> elementTransformationTitle
                            -- <> elementTransformationInfo
                            -- <> alphabetTitle
                            -- <> alphabetInfo
                    where
                        -- <> fmap show overallElementTransformations

                        concat4 ∷ ∀ {a}. (Semigroup a) ⇒ a → a → a → a → a
                        concat4 a b c d = a <> b <> c <> d
                        concat3 ∷ ∀ {a}. (Semigroup a) ⇒ a → a → a → a
                        concat3 a b c = a <> b <> c

                        getIndex a s = if (length s > a) then s !! a
                                       else "" 

                        showInfo a = if snd5 a == -1.0 then ""
                                     else show a

{- U.getEdgeInformationContent--Experimental
    returns conditional complexity of child state given parent--K(c|p) and
    ration of child complexity to child conditinal complexity-- K(c) / K(c|p)
    compression based (zip)
-}
getEdgeInformationContent :: String -> String -> (Double, Double,Double, Double, Double)
getEdgeInformationContent childString' parentString' =
    -- filter out solitary gaps--put in by IA for transformation counting
    let childString = filter (/= '-') childString'
        parentString = filter (/= '-') parentString'
    in
    -- empty line
    if null childString && null parentString then (0,-1.0, 0,0,0)
    -- child empty or missing
    else if null childString then (0,-1.0, 0,0,0)
    -- parent empty or missing
    else if null parentString then 
        -- only K(c)
        let (_, _, _, kc)  = CU.getInformationContent childString
        in (kc, 1.0, 0,0,0)
    else 
        let (_, _, _, kc)  = CU.getInformationContent childString
            (_, _, _, kcGp)  = CU.getInformationContent  (childString <> parentString)
            (_, _, _, kp)  = CU.getInformationContent parentString
        in
        (kcGp - kp,  kcGp / kc, kc, kp, kcGp)


{- | sumEdgeElementLists takes list of edges of block of characters transomftion matricwes and
sums over edge (outermost) list for eventual csv output
-}
sumEdgeElementLists ∷ [[[(String, Int)]]] → [[[[(String, Int)]]]] → [[[(String, Int)]]]
sumEdgeElementLists curSum inEdgeList =
    if null inEdgeList
        then curSum
        else
            let firstEdge = head inEdgeList
                newSum =
                    if null curSum
                        then firstEdge
                        else addBlockCharacterElements curSum firstEdge
            in  sumEdgeElementLists newSum (tail inEdgeList)


-- | addBlockCharacterElements adds two block lists of element lists
addBlockCharacterElements ∷ [[[(String, Int)]]] → [[[(String, Int)]]] → [[[(String, Int)]]]
addBlockCharacterElements inBlockEdgeL1 inBlockEdgeL2 =
    if null inBlockEdgeL1
        then inBlockEdgeL2
        else
            if null inBlockEdgeL2
                then inBlockEdgeL1
                else zipWith addCharacterElementLists inBlockEdgeL1 inBlockEdgeL2


-- | addCharacterElementLists adds two lists of charcter element lists
addCharacterElementLists ∷ [[(String, Int)]] → [[(String, Int)]] → [[(String, Int)]]
addCharacterElementLists inCharList1 inCharList2 =
    if null inCharList1
        then inCharList2
        else
            if null inCharList2
                then inCharList1
                else zipWith addElementLists inCharList1 inCharList2


{- | addElementLists adds two lists of elements and numbers
    Assumes lists are in same element order
-}
addElementLists ∷ [(String, Int)] → [(String, Int)] → [(String, Int)]
addElementLists elementL1 elementL2 =
    if null elementL1
        then elementL2
        else
            if null elementL2
                then elementL1
                else
                    if (fmap fst elementL1) /= (fmap fst elementL2)
                        then error ("Element lists do not match: " <> (show $ fmap fst elementL1) <> " versus " <> (show $ fmap fst elementL2))
                        else
                            let numberL1 = fmap snd elementL1
                                numberL2 = fmap snd elementL2
                                newNUmberL = zipWith (+) numberL1 numberL2
                            in  zip (fmap fst elementL1) newNUmberL


{- | sumEdgeTransformationLists takes list of edges of block of characters transomftion matricwes and
sums over edge (outermost) list for eventual csv output
-}
sumEdgeTransformationLists ∷ [[[[Int]]]] → [[[[[Int]]]]] → [[[[Int]]]]
sumEdgeTransformationLists curSum inEdgeList =
    if null inEdgeList
        then curSum
        else
            let firstEdge = head inEdgeList
                newSum =
                    if null curSum
                        then firstEdge
                        else addBlockCharacterMatrix curSum firstEdge
            in  sumEdgeTransformationLists newSum (tail inEdgeList)


-- | addBlockCharacterMatrix adds two block lists of character lists of transrmation matrices
addBlockCharacterMatrix ∷ [[[[Int]]]] → [[[[Int]]]] → [[[[Int]]]]
addBlockCharacterMatrix inBlockEdgeL1 inBlockEdgeL2 =
    if null inBlockEdgeL1
        then inBlockEdgeL2
        else
            if null inBlockEdgeL2
                then inBlockEdgeL1
                else zipWith addCharacterLists inBlockEdgeL1 inBlockEdgeL2


-- | addCharacterLists adds two lists of charcter transfomratin matrices
addCharacterLists ∷ [[[Int]]] → [[[Int]]] → [[[Int]]]
addCharacterLists inCharList1 inCharList2 =
    if null inCharList1
        then inCharList2
        else
            if null inCharList2
                then inCharList1
                else zipWith (U.combineMatrices (+)) inCharList1 inCharList2


{- | addBlockCharStrings adds block and character strings to transformation info
    for CSV output
-}
addBlockCharStrings ∷ [[String]] → [[String]] → [[String]] → [[String]] → [[String]]
addBlockCharStrings labelStringList elementStringList matrixStringList matrixStringList2 =
    if null labelStringList || null matrixStringList
        then []
        else
            let blockTitle = head labelStringList
                charMatrixList = head matrixStringList
                charMatrixList2 = head matrixStringList2
                elementList = head elementStringList
                charNameList = take (length charMatrixList) (tail labelStringList)
            in  -- first doesn't ahve numbers of transforms
                -- (blockTitle : (zipWith3 (concat3) charNameList (fmap (:[]) elementList) (fmap (:[]) charMatrixList))) <> addBlockCharStrings (drop (1 + length charMatrixList) labelStringList) (tail elementStringList) (tail matrixStringList)
                ( blockTitle
                    : (L.zipWith4 (concat4) charNameList (fmap (: []) elementList) (fmap (: []) charMatrixList) (fmap (: []) charMatrixList2))
                )
                    <> addBlockCharStrings
                        (drop (1 + length charMatrixList) labelStringList)
                        (tail elementStringList)
                        (tail matrixStringList)
                        (tail matrixStringList2)
    where
        -- concat3 a b c = a <> b <> c
        concat4 a b c d = a <> b <> c <> d


{- | knitTitlesChangeInfo takes [[[String]]] of title info and knits with [[[String]]] of character change info
    into [[String]] for CSV output
    Edges x Blocks x characters (final is transformation info)
-}
knitTitlesChangeInfo ∷ [[[String]]] → [[[String]]] → [[String]]
knitTitlesChangeInfo titlesLLL transLLL =
    if null titlesLLL || null transLLL
        then []
        else -- trace ("KTCI: " <> (show titlesLLL)) $
        -- trace ("KTIC2: " <> (show transLLL))
            concat $ zipWith formatEdge titlesLLL transLLL


{- | formatEdge formats edges string for CSV
    divides list into blocks by lengfth olf charcters and block titles haveing an extra line
-}
formatEdge ∷ [[String]] → [[String]] → [[String]]
formatEdge titleLL transLL =
    if null titleLL || null transLL
        then []
        else
            let edgeTitle = head titleLL
                numCharsL = fmap length transLL
                blockTitlesLLL = U.divideList (fmap (+ 1) numCharsL) (tail titleLL)
                blockTransLLL = transLL -- U.divideList numCharsL transLL
            in  -- trace ("FE: " <> (show edgeTitle)) $
                edgeTitle : (concat $ zipWith (formatBlock numCharsL) blockTitlesLLL blockTransLLL)


-- | formatBlocks formats block string for CSV
formatBlock ∷ [Int] → [[String]] → [String] → [[String]]
formatBlock charLengthL titleLL transL =
    if null titleLL || null transL
        then []
        else
            let blockTitle = head titleLL
                charTitleLL = tail titleLL
            in  -- trace ("FB: " <> (show blockTitle)) $
                blockTitle : (zipWith formatCharacter charTitleLL transL)


-- | formatCharacter formats charcater string for CSV
formatCharacter ∷ [String] → String → [String]
formatCharacter titleLine transS =
    if null titleLine || null transS
        then []
        else -- trace ("FC: " <> (show titleLine)) $
            titleLine <> (LS.splitOn "," transS)


{- | getEdgeTransformations get tranformations for an edge by block and character
    changes by block then character
-}
getEdgeTransformations
    ∷ V.Vector (LG.LNode VertexInfo) → V.Vector (V.Vector CharInfo) → (LG.Node, LG.Node) → [[(String, [(String, Int)], [[Int]])]]
getEdgeTransformations nodeVect charInfoVV (parentIndex, childIndex) =
    let parentNodeLabel = snd $ nodeVect V.! parentIndex
        childNodeLabel = snd $ nodeVect V.! childIndex
        parentBlockData = vertData parentNodeLabel
        childBlockData = vertData childNodeLabel
    in  V.toList $ V.zipWith3 getEdgeBlockChanges parentBlockData childBlockData charInfoVV


-- | getEdgeBlockChanges takes VertexBlockData from parent and child node and gets transformations by character
getEdgeBlockChanges
    ∷ V.Vector CharacterData → V.Vector CharacterData → V.Vector CharInfo → [(String, [(String, Int)], [[Int]])]
getEdgeBlockChanges parentBlockData childBlockData charInfoV =
    V.toList $ V.zipWith3 getCharacterChanges parentBlockData childBlockData charInfoV


{- | getCharacterChanges takes strings of character pairs and outputs differences as pairs,
    for unaligned sequences performs a DO and uses aligned states
-}
getCharacterChanges ∷ CharacterData → CharacterData → CharInfo → (String, [(String, Int)], [[Int]])
getCharacterChanges parentChar nodeChar charInfo =
    let localType = charType charInfo
        localAlphabet = (ST.toString <$> alphabet charInfo)

        -- this to avoid recalculations and list access issues
        lANES = (fromJust $ NE.nonEmpty $ alphabetSymbols localAlphabet)
        lAVect = V.fromList $ NE.toList $ lANES

        -- getCharState
        --    ∷ ∀ {b}
        --     . (Show b, Bits b)
        --    ⇒ b
        --    → String
        getCharState a = U.bitVectToCharState localAlphabet lANES lAVect a -- concat $ NE.toList $ decodeState localAlphabet a --
        (slimParent, slimChild) =
            if localType `notElem` [NucSeq, SlimSeq]
                then ([], [])
                else
                    let (_, r) =
                            slimPairwiseDO
                                (slimTCM charInfo)
                                (M.makeDynamicCharacterFromSingleVector $ slimFinal parentChar)
                                (M.makeDynamicCharacterFromSingleVector $ slimFinal nodeChar)
                    in  -- trace (show (r, length $ SV.foldMap getCharState $ extractMediansLeftGapped r, length $ SV.foldMap getCharState $ extractMediansRightGapped r))
                        (SV.foldMap getCharState $ extractMediansLeftGapped r, SV.foldMap getCharState $ extractMediansRightGapped r)
        (wideParent, wideChild) =
            if localType `notElem` [WideSeq, AminoSeq]
                then ([], [])
                else
                    let coefficient = MR.minInDelCost (wideTCM charInfo)
                        (_, r) =
                            widePairwiseDO
                                coefficient
                                (MR.retreivePairwiseTCM $ wideTCM charInfo)
                                (M.makeDynamicCharacterFromSingleVector $ wideFinal parentChar)
                                (M.makeDynamicCharacterFromSingleVector $ wideFinal nodeChar)
                    in  (UV.foldMap getCharState $ extractMediansLeftGapped r, UV.foldMap getCharState $ extractMediansRightGapped r)
        (hugeParent, hugeChild) =
            if localType `notElem` [HugeSeq]
                then ([], [])
                else
                    let coefficient = MR.minInDelCost (hugeTCM charInfo)
                        (_, r) =
                            hugePairwiseDO
                                coefficient
                                (MR.retreivePairwiseTCM $ hugeTCM charInfo)
                                (M.makeDynamicCharacterFromSingleVector $ hugeFinal parentChar)
                                (M.makeDynamicCharacterFromSingleVector $ hugeFinal nodeChar)
                    in  (foldMap getCharState $ extractMediansLeftGapped r, foldMap getCharState $ extractMediansRightGapped r)

        (parentState, nodeState)
            | localType == Add = (show $ rangeFinal parentChar, show $ rangeFinal nodeChar)
            | localType == NonAdd =
                ( concat $ V.map (U.bitVectToCharStateQual localAlphabet) $ stateBVFinal parentChar
                , concat $ V.map (U.bitVectToCharStateQual localAlphabet) $ stateBVFinal nodeChar
                )
            | localType `elem` packedNonAddTypes =
                (UV.foldMap getCharState $ packedNonAddFinal parentChar, UV.foldMap getCharState $ packedNonAddFinal nodeChar)
            | localType == Matrix =
                (show $ fmap (fmap fst3) $ matrixStatesFinal parentChar, UV.foldMap getCharState $ packedNonAddFinal nodeChar)
            | localType `elem` sequenceCharacterTypes = case localType of
                x | x `elem` [NucSeq, SlimSeq] → (slimParent, slimChild)
                x | x `elem` [WideSeq, AminoSeq] → (wideParent, wideChild)
                x | x `elem` [HugeSeq] → (hugeParent, hugeChild)
                x
                    | x `elem` [AlignedSlim] →
                        (SV.foldMap getCharState $ alignedSlimFinal parentChar, SV.foldMap getCharState $ alignedSlimFinal nodeChar)
                x
                    | x `elem` [AlignedWide] →
                        (UV.foldMap getCharState $ alignedWideFinal parentChar, UV.foldMap getCharState $ alignedWideFinal nodeChar)
                x
                    | x `elem` [AlignedHuge] →
                        (foldMap getCharState $ alignedHugeFinal parentChar, foldMap getCharState $ alignedHugeFinal nodeChar)
                _ → error ("Un-implemented data type " <> show localType)
            | otherwise = error ("Un-implemented data type " <> show localType)

        -- character states and transformation numbers (only counting unambiguous)
        elementsParent = getElementNumbers (alphabetSymbols localAlphabet) parentState
        elementsChild = getElementNumbers (alphabetSymbols localAlphabet) nodeState
        elementsCombined = elementsChild -- zip (alphabetSymbols localAlphabet)  (zipWith (+) (fmap snd elementsParent) (fmap snd elementsChild))
        elementsCombinedString = L.intercalate "," $ fmap makeElementString elementsCombined

        emptyMatrix = replicate (length localAlphabet) (replicate (length localAlphabet) 0)
        elementTransformations = getTransformations (alphabetSymbols localAlphabet) parentState nodeState emptyMatrix
    in  -- convert String to pair list
        -- (parentState <> "," <> nodeState <> "," <> elementsCombinedString <> "," <> (replaceComma $ show elementTransformations), elementsCombined, elementTransformations)
        ( parentState <> "," <> nodeState <> "," <> (replaceComma $ show elementTransformations)
        , elementsCombined
        , elementTransformations
        )
    where
        -- ((replaceComma $ show elementTransformations), elementsCombined, elementTransformations)

        makeElementString (a, b) = a <> " " <> (show b)
        replaceComma a =
            if null a
                then []
                else
                    if head a == ','
                        then ' ' : replaceComma (tail a)
                        else (head a) : replaceComma (tail a)


{- | getTransformations takes alphabet and parent and child SINGE CHARACXTRE STATES
and returns matrix of numbers of changes
This an absurdely slow implementation n^2 at least
-}
getTransformations ∷ [String] → String → String → [[Int]] → [[Int]]
getTransformations alphabet parentCharList childCharList curMatrix =
    if null parentCharList || null childCharList
        then curMatrix
        else
            let pChar = head parentCharList
                cChar = head childCharList
            in  -- trace ("GT: " <> (show curMatrix)) $
                if pChar == cChar
                    then getTransformations alphabet (tail parentCharList) (tail childCharList) curMatrix
                    else
                        if (pChar : []) `elem` alphabet && (cChar : []) `elem` alphabet
                            then
                                let pIndex = fromJust $ L.elemIndex (pChar : []) alphabet
                                    cIndex = fromJust $ L.elemIndex (cChar : []) alphabet
                                    newMatrix = incrementMatrix curMatrix pIndex cIndex
                                in  getTransformations alphabet (tail parentCharList) (tail childCharList) newMatrix
                            else getTransformations alphabet (tail parentCharList) (tail childCharList) curMatrix


-- | incrementMatrix updates matrix very stupidly
incrementMatrix ∷ [[Int]] → Int → Int → [[Int]]
incrementMatrix inMatrix p c =
    if null inMatrix
        then []
        else
            let firstRows = take p inMatrix
                lastRows = drop (p + 1) inMatrix
                updateRow = inMatrix !! p
                firstPart = take c updateRow
                lastPart = drop (c + 1) updateRow
                newRow = firstPart <> [1 + (updateRow !! c)] <> lastPart
            in  firstRows <> [newRow] <> lastRows


-- | getElementNumbers takes a String of SINGLECHARACTRE states and checks versus alphabet for unambiguous states numbers
getElementNumbers ∷ [String] → String → [(String, Int)]
getElementNumbers alphabet charList =
    if null charList
        then []
        else
            if null alphabet
                then []
                else
                    let stringList = fmap (: []) charList
                        alphNumber = length $ L.findIndices (== (head alphabet)) stringList
                    in  ((head alphabet), alphNumber) : getElementNumbers (tail alphabet) charList


{- | getAlignmentBasedChanges' takes two equal length implied Alignments and outputs list of element changes between the two
assumes single String in each list
-}
getAlignmentBasedChanges' ∷ Bool → Int → ([String], [String]) → [String]
getAlignmentBasedChanges' onlyDiffs index (a, b)
    | length a > 1 || length b < 1 = error ("Should only have length 1 lists here: " <> (show (length a, length b)))
    | null a = []
    | otherwise =
        -- empty spaces sometimes
        let stringList1 = filter (not . null) $ LS.splitOn (" ") (head a)
            stringList2 = filter (not . null) $ LS.splitOn (" ") (head b)
        in  -- this so returns empty for non--sequence characters
            if null stringList1
                then []
                else -- this shouldn't happen but there seem to be extraneous '-' at time in wide and huge IAs
                -- since all zip functions as well--changed getAlignmentBasedChangesGuts to check for either null stringList1 or stringList2
                -- to avoid issue-a bit of a hack
                -- else if length stringList1 /= length stringList2 then
                --         error ("Unequal characters in parent and node state lists in getAlignmentBasedChanges': "
                --         <> (show (length stringList1, length stringList2) <> "\n" <> (show stringList1) <> "\n" <> (show stringList2)))
                    getAlignmentBasedChangesGuts onlyDiffs index stringList1 stringList2


-- | getAlignmentBasedChangesGuts takes processed element lists and creates string of changes
getAlignmentBasedChangesGuts ∷ Bool → Int → [String] → [String] → [String]
getAlignmentBasedChangesGuts onlyDiffs index a b
    | null a || null b = []
    | (head a) == (head b) && onlyDiffs = getAlignmentBasedChangesGuts onlyDiffs (index + 1) (tail a) (tail b)
    | otherwise =
        ((show index) <> ":" <> (head a) <> "," <> (head b)) : getAlignmentBasedChangesGuts onlyDiffs (index + 1) (tail a) (tail b)


{- | getAlignmentBasedChanges takes two equal length implied Alignments and outputs list of element changes between the two
only working for nucleotide prealigned or not
-}
getAlignmentBasedChanges ∷ Int → ([String], [String]) → [String]
getAlignmentBasedChanges index (a, b) =
    if null a
        then []
        else -- empty spaces sometimes

            let string1 =
                    if (head $ head a) == ' '
                        then tail $ head a
                        else head a
                string2 =
                    if (head $ head b) == ' '
                        then tail $ head b
                        else head b
            in  if null string1
                    then []
                    else -- this so returns empty for non--sequence characters

                        if (length string1 < 2) || (length string2 < 2)
                            then []
                            else
                                if length string1 /= length string2
                                    then []
                                    else -- index stuff in case gets out of sync in list (missing data, no changes etc)

                                        if (take 1 string1) == (take 1 string2)
                                            then getAlignmentBasedChanges (index + 1) ([tail string1], [tail string2])
                                            else
                                                ((show index) <> ":" <> (take 1 string1) <> "," <> (take 1 string2))
                                                    : getAlignmentBasedChanges (index + 1) ([tail string1], [tail string2])


{- | removeNoChangeLines takes lines of vertex changes and removes lines where parent and child startes are the same
so missing or ambiguous in one and not the other will be maintained
-}
removeNoChangeLines ∷ [[String]] → [[String]]
removeNoChangeLines inStringList =
    if null inStringList
        then []
        else
            let parentState = head inStringList !! 9
                childState = head inStringList !! 10
            in  ( if (parentState == " ") || (parentState /= childState)
                    then (head inStringList) : removeNoChangeLines (tail inStringList)
                    else removeNoChangeLines (tail inStringList)
                )


{- | getVertexCharInfo returns a list of list of Strings of vertex information
one list for each character at the vertex
-}
getVertexCharInfo
    ∷ Bool → V.Vector BlockData → SimpleGraph → V.Vector (V.Vector CharInfo) → LG.LNode VertexInfo → [[String]]
getVertexCharInfo useIA blockDataVect inGraph charInfoVectVect inVert =
    let nodeParents = LG.parents inGraph (fst inVert)
        parentNodes
            | nodeType (snd inVert) == RootNode = "None"
            | nodeType (snd inVert) == LeafNode = show nodeParents
            | otherwise = show $ parents (snd inVert)
        childNodes = if nodeType (snd inVert) == LeafNode then "None" else show $ children (snd inVert)
        basicInfoList =
            [ " "
            , show $ fst inVert
            , T.unpack $ vertName (snd inVert)
            , show $ nodeType (snd inVert)
            , childNodes
            , parentNodes
            , " "
            , " "
            , " "
            , " "
            , " "
            ]
        blockCharVect = V.zip3 (V.map fst3 blockDataVect) (vertData (snd inVert)) charInfoVectVect
        blockInfoList = concat $ V.toList $ V.map (getBlockList useIA) blockCharVect
    in  basicInfoList : blockInfoList


{- | getVertexAndParentCharInfo returns a list of list of Strings of vertex information
for child and its parent
-}
getVertexAndParentCharInfo
    ∷ Bool
    → V.Vector BlockData
    → SimpleGraph
    → V.Vector (V.Vector CharInfo)
    → V.Vector (LG.LNode VertexInfo)
    → LG.LNode VertexInfo
    → [[String]]
getVertexAndParentCharInfo useIA blockDataVect inGraph charInfoVectVect allVertVect inVert =
    let nodeParents = LG.parents inGraph (fst inVert)
        parentNodes
            | nodeType (snd inVert) == RootNode = "None"
            | nodeType (snd inVert) == LeafNode = show nodeParents
            | otherwise = show $ parents (snd inVert)
        childNodes = if nodeType (snd inVert) == LeafNode then "None" else show $ children (snd inVert)
        basicInfoList =
            [ " "
            , show $ fst inVert
            , T.unpack $ vertName (snd inVert)
            , show $ nodeType (snd inVert)
            , childNodes
            , parentNodes
            , " "
            , " "
            , " "
            , " "
            , " "
            ]
        blockCharVectNode = V.zip3 (V.map fst3 blockDataVect) (vertData (snd inVert)) charInfoVectVect

        -- for root--gets its own values as parent--filtered out in diff list later
        blockCharVectParent =
            if parentNodes == "None"
                then blockCharVectNode
                else V.zip3 (V.map fst3 blockDataVect) (vertData (snd $ allVertVect V.! head nodeParents)) charInfoVectVect
        blockInfoListParent = concat $ V.toList $ V.map (getBlockList useIA) blockCharVectParent
    in  basicInfoList : blockInfoListParent


-- | getBlockList takes a pair of Vector of chardata and vector of charInfo and returns Strings
getBlockList ∷ Bool → (NameText, V.Vector CharacterData, V.Vector CharInfo) → [[String]]
getBlockList useIA (blockName, blockDataVect, charInfoVect) =
    let firstLine = [" ", " ", " ", " ", " ", " ", T.unpack blockName, " ", " ", " ", " ", " "]
        charlines = V.toList $ V.map (makeCharLine useIA) (V.zip blockDataVect charInfoVect)
    in  firstLine : charlines


{- | makeCharLine takes character data
will be less legible for optimized data--so should use a diagnosis
based on "naive" data for human legible output
need to add back-converting to observed states using alphabet in charInfo
nothing here for packed since not "entered"
useIA for using alignment fields for changes in diagnosis
is useIA == False then just printing final sequence and removes spaces
for single character sequences (e.g. DNA/Protein)
-}
makeCharLine ∷ Bool → (CharacterData, CharInfo) → [String]
makeCharLine useIA (blockDatum, charInfo) =
    let localType = charType charInfo
        localAlphabet = (ST.toString <$> alphabet charInfo)
        isPrealigned =
            if prealigned charInfo
                then "Prealigned "
                else ""
        enhancedCharType
            | localType `elem` sequenceCharacterTypes = (isPrealigned <> (show localType))
            | localType `elem` exactCharacterTypes = (show localType)
            | otherwise = error ("Character Type :" <> (show localType) <> "unrecogniized or not implemented")

        -- set where get string from, IA for change lists
        (slimField, wideField, hugeField) =
            if useIA
                then (slimIAFinal blockDatum, wideIAFinal blockDatum, hugeIAFinal blockDatum)
                else (slimFinal blockDatum, wideFinal blockDatum, hugeFinal blockDatum)

        -- this to avoid recalculations and list access issues
        lANES = (fromJust $ NE.nonEmpty $ alphabetSymbols localAlphabet)
        lAVect = V.fromList $ NE.toList $ lANES

        getCharState
            ∷ ∀ {b}
             . (Show b, Bits b)
            ⇒ b
            → String
        getCharState a = U.bitVectToCharState localAlphabet lANES lAVect a

        -- (stringPrelim, stringFinal) = if localType == Add then (show $ snd3 $ rangePrelim blockDatum, show $ rangeFinal blockDatum)
        stringFinal
            | localType == Add = (show $ rangeFinal blockDatum)
            | localType == NonAdd = (concat $ V.map (U.bitVectToCharStateQual localAlphabet) $ stateBVFinal blockDatum)
            | localType `elem` packedNonAddTypes = (UV.foldMap getCharState $ packedNonAddFinal blockDatum)
            | localType == Matrix = (show $ fmap (fmap fst3) $ matrixStatesFinal blockDatum)
            | localType `elem` sequenceCharacterTypes = case localType of
                x | x `elem` [NucSeq] → (SV.foldMap getCharState slimField)
                x | x `elem` [SlimSeq] → (SV.foldMap getCharState slimField)
                x | x `elem` [WideSeq] → (UV.foldMap getCharState wideField)
                x | x `elem` [AminoSeq] → (UV.foldMap getCharState wideField)
                x | x `elem` [HugeSeq] → (foldMap getCharState hugeField)
                x | x `elem` [AlignedSlim] → (SV.foldMap getCharState $ alignedSlimFinal blockDatum)
                x | x `elem` [AlignedWide] → (UV.foldMap getCharState $ alignedWideFinal blockDatum)
                x | x `elem` [AlignedHuge] → (foldMap getCharState $ alignedHugeFinal blockDatum)
                _ → error ("Un-implemented data type " <> show localType)
            | otherwise = error ("Un-implemented data type " <> show localType)

        -- this removes ' ' between elements if sequence elements are a single character (e.g. DNA)
        stringFinal'
            | useIA = stringFinal
            | localType `elem` [NucSeq, AminoSeq] = filter (/= ' ') stringFinal
            | otherwise =
                let maxSymbolLength = maximum $ fmap length $ SET.toList (alphabetSymbols localAlphabet)
                in  if maxSymbolLength > 1
                        then removeRecurrrentSpaces $ fmap nothingToNothing stringFinal
                        else filter (/= ' ') stringFinal
    in  -- trace ("MCL:" <> (show localType) <> " " <> stringFinal)
        -- [" ", " ", " ", " ", " ", " ", " ", T.unpack $ name charInfo, enhancedCharType, stringPrelim, stringFinal, show $ localCost blockDatum]
        [" ", " ", " ", " ", " ", " ", " ", T.unpack $ name charInfo, enhancedCharType, stringFinal']
    where
        nothingToNothing a =
            if a == '\8220'
                then '\00'
                else a


-- | removeRecurrrentSpaces removes spaces that are followed by spaces
removeRecurrrentSpaces ∷ String → String
removeRecurrrentSpaces inString
    | null inString = []
    | length inString == 1 = inString
    | head inString == ' ' =
        if (head $ tail inString) == ' '
            then removeRecurrrentSpaces (tail inString)
            else (head inString) : removeRecurrrentSpaces (tail inString)
    | otherwise = (head inString) : removeRecurrrentSpaces (tail inString)


-- | TNT report functions

{- | getTNTStrings returns as a set of interleaved b;ocks==one for each "character"  so not mix numerical and
sequence characters
softwired use display trees, hardWired transform to softwired then proceed with display trees
key to keep cost matrices and weights
Uses Phylogenetic graph to noit repeat functions for display and charcter trees
-}
getTNTString ∷ GlobalSettings → ProcessedData → (PhylogeneticGraph, Int) → PhyG String
getTNTString inGS inData (inGraph, graphNumber) =
    if LG.isEmpty (fst6 inGraph)
        then error "No graphs for create TNT data for in getTNTString"
        else
            let leafList = snd4 $ LG.splitVertexList (thd6 inGraph)
                leafNameList = fmap ((<> "\t") . T.unpack) (fmap (vertName . snd) leafList)
                headerString =
                    "taxname + 90;\nxread\n'TNT data for Graph "
                        <> show graphNumber
                        <> " generated by PhylogeneticGraph (PhyG)\n\tSource characters:\n"
                finalString = "proc/;\nlength;\n\n"
                numTaxa = V.length $ fst3 inData
                charInfoVV = six6 inGraph

                -- get character information in 3-tuples and list of lengths--to match lengths
                ccCodeInfo = getCharacterInfo charInfoVV

                -- only first of block list type--so assumes single sequence per block and only that type--no exact types
                charTypeList = V.toList $ fmap charType $ fmap V.head charInfoVV

                -- Tree in TNT format
                tntTreeString =
                    "tread 'Trees from PhyG'\n"
                        <> (GO.makeNewickList True False False (outgroupIndex inGS) [GO.convertDecoratedToSimpleGraph $ thd6 inGraph] [snd6 inGraph])
            in  if graphType inGS == Tree
                    then
                        let leafDataList = fmap (vertData . snd) leafList
                        in  do
                                -- get character strings
                                blockStringListstList ← mapM (getTaxonCharStringListPair charInfoVV) (zip leafDataList leafNameList)
                                -- let blockStringListstList = V.toList blockStringListstList'
                                interleavedBlocks' ← makePairInterleave blockStringListstList charTypeList
                                let interleavedBlocks = concat interleavedBlocks'
                                -- taxonCharacterStringList = V.toList $ fmap ((<> "\n") . getTaxonCharString charInfoVV) leafDataList
                                -- nameCharStringList = concat $ zipWith (<>) leafNameList taxonCharacterStringList

                                -- length information for cc code extents
                                let charLengthList = concat $ V.toList $ V.zipWith getBlockLength (head leafDataList) charInfoVV

                                -- Block/Character names for use in comment to show sources of new characters
                                let charNameList = concat $ V.toList $ fmap getBlockNames charInfoVV

                                let nameLengthPairList = zip charNameList charLengthList
                                let nameLengthString = concat $ pairListToStringList nameLengthPairList 0

                                -- merge lengths and cc codes
                                let ccCodeString = mergeCharInfoCharLength ccCodeInfo charLengthList 0

                                -- trace ("GTNTS:" <> (show charLengthList))
                                pure $
                                    headerString
                                        <> nameLengthString
                                        <> "'\n"
                                        <> show (sum charLengthList)
                                        <> " "
                                        <> show numTaxa
                                        <> "\n"
                                        <> interleavedBlocks
                                        <> ";\n"
                                        <> ccCodeString
                                        <> tntTreeString
                                        <> finalString
                    else -- for softwired networks--use display trees

                        if graphType inGS == SoftWired
                            then do
                                -- get display trees for each data block-- takes first of potentially multiple
                                middleStuffString ← createDisplayTreeTNT inGS inData inGraph

                                pure $ headerString <> middleStuffString <> tntTreeString <> finalString
                            else -- for hard-wired networks--transfoirm to softwired and use display trees

                                if graphType inGS == HardWired
                                    then
                                        let newGS = inGS{graphType = SoftWired}

                                            pruneEdges = False
                                            warnPruneEdges = False

                                            startVertex ∷ ∀ {a}. Maybe a
                                            startVertex = Nothing
                                        in  do
                                                newGraph ← TRAV.multiTraverseFullyLabelGraph newGS inData pruneEdges warnPruneEdges startVertex (fst6 inGraph)

                                                middleStuffString ← createDisplayTreeTNT inGS inData newGraph

                                                logWith
                                                    LogWarn
                                                    "There is no implied alignment for hard-wired graphs--at least not yet. Ggenerating TNT text via softwired transformation\n"
                                                -- headerString <> nameLengthString <> "'\n" <> (show $ sum charLengthList) <> " " <> (show numTaxa) <> "\n"
                                                --    <> nameCharStringList <> ";\n" <> ccCodeString <> finalString
                                                pure $ headerString <> middleStuffString <> tntTreeString <> finalString
                                    else do
                                        logWith LogWarn ("TNT  not yet implemented for graphtype " <> show (graphType inGS) <> "\n")
                                        pure $ ("There is no implied alignment for " <> show (graphType inGS))


{- | makePairInterleave creates interleave block strings from character name block list
list of taxa and blocck with taxon name
-}
makePairInterleave ∷ [[String]] → [CharType] → PhyG [String]
makePairInterleave inTaxCharStringList alphabetType =
    -- concat $ fmap concat inTaxCharStringList
    if null inTaxCharStringList
        then do pure []
        else
            if null $ head inTaxCharStringList
                then do pure []
                else
                    let firstInterleave = concatMap (<> "\n") $ fmap head inTaxCharStringList
                        remainder = fmap tail inTaxCharStringList
                        interleaveMarkerString
                            | head alphabetType `elem` exactCharacterTypes = "& [num]"
                            | head alphabetType `elem` [NucSeq, AlignedSlim] = "& [dna gaps]"
                            | head alphabetType `elem` [AminoSeq, AlignedWide] = "& [prot gaps]"
                            | otherwise = "& [other]"
                    in  do
                            endString ← makePairInterleave remainder (tail alphabetType)
                            let returnString = (interleaveMarkerString <> "\n") : (firstInterleave : endString)
                            if (head alphabetType `notElem` (exactCharacterTypes <> [NucSeq, AlignedSlim] <> [AminoSeq, AlignedWide]))
                                then do
                                    logWith LogWarn ("Warning--sequence data in tnt ouput not of type TNT accepts" <> "\n")
                                    pure returnString
                                else do
                                    pure returnString


-- | createDisplayTreeTNT take a softwired graph and creates TNT data string
createDisplayTreeTNT ∷ GlobalSettings → ProcessedData → PhylogeneticGraph → PhyG String
createDisplayTreeTNT inGS inData inGraph =
    let leafList = snd4 $ LG.splitVertexList (thd6 inGraph)
        leafNameList = fmap ((<> "\t") . T.unpack) (fmap (vertName . snd) leafList)
        charInfoVV = six6 inGraph

        -- only first of block list type--so assumes single sequence per block and only that type--no exact types
        charTypeList = V.toList $ fmap charType $ fmap V.head charInfoVV

        numTaxa = V.length $ fst3 inData
        ccCodeInfo = getCharacterInfo charInfoVV

        -- parallel stuff
        contract ∷ [DecoratedGraph] → SimpleGraph
        contract = GO.contractIn1Out1EdgesRename . GO.convertDecoratedToSimpleGraph . head

        block ∷ BlockData → ProcessedData
        block = makeBlockData (fst3 inData) (snd3 inData)

        traverseAction ∷ (ProcessedData, SimpleGraph) → PhyG PhylogeneticGraph
        traverseAction = TRAV.multiTraverseFullyLabelGraphPair (inGS{graphType = Tree}) False False Nothing

        taxonString ∷ (VertexBlockData, String) → PhyG [String]
        taxonString = getTaxonCharStringListPair charInfoVV
    in  do
            contractPar ← getParallelChunkMap
            let blockDisplayList = contractPar contract (V.toList $ fth6 inGraph)

            -- blockDisplayList = fmap (GO.contractIn1Out1EdgesRename . GO.convertDecoratedToSimpleGraph . head) (V.toList $ fth6 inGraph) `using` PU.myParListChunkRDS
            -- blockDisplayList = PU.seqParMap PU.myStrategyHighLevel  (GO.contractIn1Out1EdgesRename . GO.convertDecoratedToSimpleGraph . head) (V.toList $ fth6 inGraph)

            -- create separate processed data for each block
            blockPar ← getParallelChunkMap
            let blockProcessedDataList = blockPar block (V.toList $ thd3 inData)
            -- blockProcessedDataList = PU.seqParMap PU.myStrategyHighLevel (makeBlockData (fst3 inData) (snd3 inData)) (thd3 inData)

            -- Perform full optimizations on display trees (as trees) with single block data (blockProcessedDataList) to creeate IAs
            decoratedBlockTreeList ←
                getParallelChunkTraverse >>= \pTraverse →
                    traverseAction `pTraverse` zip blockProcessedDataList blockDisplayList

            -- create leaf data by merging display graph block data (each one a phylogentic graph)
            let (leafDataList, mergedCharInfoVV) = mergeDataBlocks decoratedBlockTreeList [] []

            -- get character block strings as interleaved groups
            blockStringListstList ←
                getParallelChunkTraverse >>= \pTraverse →
                    taxonString `pTraverse` zip (V.toList leafDataList) leafNameList

            interleavedBlocks ← fold <$> makePairInterleave blockStringListstList charTypeList

            -- taxonCharacterStringList = V.toList $ fmap ((<> "\n") . getTaxonCharString mergedCharInfoVV) leafDataList
            -- nameCharStringList = concat $ zipWith (<>) leafNameList taxonCharacterStringList

            -- length information for cc code extents
            let charLengthList = concat $ V.toList $ V.zipWith getBlockLength (V.head leafDataList) mergedCharInfoVV

            -- Block/Character names for use in comment to show sources of new characters
            let charNameList = concat $ V.toList $ fmap getBlockNames charInfoVV

            let nameLengthPairList = zip charNameList charLengthList
            let nameLengthString = concat $ pairListToStringList nameLengthPairList 0

            -- merge lengths and cc codes
            let ccCodeString = mergeCharInfoCharLength ccCodeInfo charLengthList 0
            pure $
                nameLengthString
                    <> "'\n"
                    <> show (sum charLengthList)
                    <> " "
                    <> show numTaxa
                    <> "\n"
                    <> interleavedBlocks
                    <> ";\n"
                    <> ccCodeString


-- | pairListToStringList takes  alist of (String, Int) and a starting index and returns scope of charcter for leading comment
pairListToStringList ∷ [(String, Int)] → Int → [String]
pairListToStringList pairList startIndex =
    if null pairList
        then []
        else
            let (a, b) = head pairList
            in  ("\t\t" <> show startIndex <> "-" <> show (b + startIndex - 1) <> " : " <> a <> "\n")
                    : pairListToStringList (tail pairList) (startIndex + b)


{- | mergeDataBlocks takes a list of Phylogenetic Graphs (Trees) and merges the data blocks (each graph should have only 1)
and merges the charInfo Vectors returning data and charInfo
-}
mergeDataBlocks
    ∷ [PhylogeneticGraph]
    → [[V.Vector CharacterData]]
    → [V.Vector CharInfo]
    → (V.Vector (V.Vector (V.Vector CharacterData)), V.Vector (V.Vector CharInfo))
mergeDataBlocks inGraphList curDataList curInfoList =
    if null inGraphList
        then (V.fromList $ fmap (V.fromList . reverse) curDataList, V.fromList $ reverse curInfoList)
        else
            let firstGraph = head inGraphList
                firstTree = thd6 firstGraph
                firstCharInfo = V.head $ six6 firstGraph
                leafList = snd4 $ LG.splitVertexList firstTree

                -- since each graph has a single block--take head to get vector of characters
                leafCharacterList = V.toList $ fmap (V.head . vertData . snd) (V.fromList leafList)

                -- zip data for each taxon
                newDataList =
                    if null curDataList
                        then (: []) <$> leafCharacterList
                        else zipWith (:) leafCharacterList curDataList
            in  mergeDataBlocks (tail inGraphList) newDataList (firstCharInfo : curInfoList)


{- | getTaxonCharString returns the total character string for a taxon
length and zipping for missing data
-}
getTaxonCharString ∷ V.Vector (V.Vector CharInfo) → VertexBlockData → PhyG String
getTaxonCharString charInfoVV charDataVV =
    -- False for not use IA field
    let lengthBlock = maximum $ V.zipWith (U.getCharacterLength False) (V.head charDataVV) (V.head charInfoVV)
        -- parallel stuff
        action ∷ (V.Vector CharInfo, V.Vector CharacterData) → PhyG String
        action = getBlockStringPair lengthBlock
    in  getParallelChunkTraverse >>= \pTraverse →
            fmap fold . pTraverse action . zip (V.toList charInfoVV) $ V.toList charDataVV


-- concat (zipWith (getBlockString lengthBlock) (V.toList charInfoVV) (V.toList charDataVV) `using` PU.myParListChunkRDS)

-- | getTaxonCharStringListPair wrapper for getTaxonCharStringList with differnet parameter passing
getTaxonCharStringListPair ∷ V.Vector (V.Vector CharInfo) → (VertexBlockData, String) → PhyG [String]
getTaxonCharStringListPair charInfoVV (charDataVV, leafName) = getTaxonCharStringList charInfoVV charDataVV leafName


{- | getTaxonCharStringList returns the total character string list (over blocks) for a taxon
length and zipping for missing data
-}
getTaxonCharStringList ∷ V.Vector (V.Vector CharInfo) → VertexBlockData → String → PhyG [String]
getTaxonCharStringList charInfoVV charDataVV leafName =
    let -- False for not use IA field
        lengthBlock = maximum $ V.zipWith (U.getCharacterLength False) (V.head charDataVV) (V.head charInfoVV)
        -- parallel stuff
        action ∷ (V.Vector CharInfo, V.Vector CharacterData) → PhyG String
        action = getBlockStringPair lengthBlock
        prefix = (fmap (leafName <>))
    in  getParallelChunkTraverse >>= \pTraverse →
            fmap prefix . pTraverse action . zip (V.toList charInfoVV) $ V.toList charDataVV


-- fmap (leafName <>) $ (zipWith (getBlockString lengthBlock) (V.toList charInfoVV) (V.toList charDataVV) `using` PU.myParListChunkRDS)

-- | getBlockStringPair wrapper around getBlockString
getBlockStringPair ∷ Int → (V.Vector CharInfo, V.Vector CharacterData) → PhyG String
getBlockStringPair lengthBlock (charInfoV, charDataV) = getBlockString lengthBlock charInfoV charDataV


{- | getBlockString returns the String for a character block
returns all '?' if missing
-}
getBlockString ∷ Int → V.Vector CharInfo → V.Vector CharacterData → PhyG String
getBlockString lengthBlock charInfoV charDataV =
    -- this to deal with missing characters
    -- trace ("GBS: " <> (show $ V.length charDataV)) (
    if V.null charDataV
        then pure $ L.replicate lengthBlock '?'
        else -- parallel stuff

            let action ∷ (CharacterData, CharInfo) → String
                action = getCharacterStringPair
            in  do
                    actionPar ← getParallelChunkMap
                    let result = actionPar action (zip (V.toList charDataV) (V.toList charInfoV))
                    pure $ concat result


-- concat (zipWith getCharacterString  (V.toList charDataV) (V.toList charInfoV) `using` PU.myParListChunkRDS)
-- )

-- | mergeCharInfoCharLength merges cc code char info and char lengths for scope
mergeCharInfoCharLength ∷ [(String, String, String)] → [Int] → Int → String
mergeCharInfoCharLength codeList lengthList charIndex =
    if null codeList
        then []
        else
            let (ccCodeString, costsString, weightString) = head codeList
                charLength = head lengthList
                startScope = show charIndex
                endScope = show (charIndex + charLength - 1)
                scope = startScope <> "." <> endScope
                weightString' =
                    if null weightString
                        then []
                        else "cc /" <> weightString <> scope <> ";\n"
                costsString' =
                    if null costsString
                        then []
                        else "costs " <> scope <> " = " <> costsString <> ";\n"
                ccCodeString' = "cc " <> ccCodeString <> " " <> scope <> ";\n"
            in  (ccCodeString' <> weightString' <> costsString')
                    <> mergeCharInfoCharLength (tail codeList) (tail lengthList) (charIndex + charLength)


{- | getCharacterInfo takes charInfo vect vect and reiurns triples of ccCode, costs, and weight values
for each character
-}
getCharacterInfo ∷ V.Vector (V.Vector CharInfo) → [(String, String, String)]
getCharacterInfo inCharInfoVV =
    concat $ V.toList $ fmap getBlockInfo inCharInfoVV


-- | getBlockInfo gets character code info for a block
getBlockInfo ∷ V.Vector CharInfo → [(String, String, String)]
getBlockInfo inCharInfoV = V.toList $ fmap getCharCodeInfo inCharInfoV


{- | getCharCodeInfo extracts 3-tuple of cc code, costs and weight as strings
from charInfo
-}
getCharCodeInfo ∷ CharInfo → (String, String, String)
getCharCodeInfo inCharInfo =
    let inCharType = charType inCharInfo
        charWeightString =
            if weight inCharInfo == 1
                then ""
                else show $ weight inCharInfo
        inAlph = alphabet inCharInfo
        inMatrix = costMatrix inCharInfo
        (costMatrixType, _) = IR.getRecodingType inMatrix
        matrixString =
            if costMatrixType == "nonAdd"
                then ""
                else makeMatrixString inAlph inMatrix
    in  let codeTriple = case inCharType of
                x | x == Add → ("+", "", charWeightString)
                x | x == NonAdd → ("-", "", charWeightString)
                x | x `elem` packedNonAddTypes → ("-", "", charWeightString)
                x | x == Matrix → ("(", matrixString, charWeightString)
                x
                    | x `elem` sequenceCharacterTypes →
                        if costMatrixType == "nonAdd"
                            then ("-", "", charWeightString)
                            else ("(", matrixString, charWeightString)
                _ → error ("Un-implemented data type " <> show inCharType)
        in  codeTriple


{- | makeMatrixString takes alphabet and input cost matrix and creates TNT
matrix cost line
could be lesser but might not be symmetrical
-}
makeMatrixString ∷ Alphabet ST.ShortText → S.Matrix Int → String
makeMatrixString inAlphabet inMatrix =
    let elementList = (ST.toString <$> toList inAlphabet)

        -- get element pairs (Strings)
        elementPairList = filter notDiag $ getListPairs elementList

        -- index pairs for accessing matrix
        elementIndexPairList = filter notDiag $ getListPairs [0 .. (length elementList - 1)]
        elementPairCosts = fmap (inMatrix S.!) elementIndexPairList

        -- make strings of form state_i / state_j cost ...
        costString = makeCostString elementPairList elementPairCosts
    in  costString
    where
        notDiag ∷ ∀ {a}. (Ord a) ⇒ (a, a) → Bool
        notDiag (a, b) = a < b


-- | makeCostString takes list of state pairs and list of costs and creates tnt cost string
makeCostString ∷ [(String, String)] → [Int] → String
makeCostString namePairList costList =
    if null namePairList
        then []
        else
            let (a, b) = head namePairList
                c = head costList
            in  (a <> "/" <> b <> " " <> show c <> " ") <> makeCostString (tail namePairList) (tail costList)


-- | getBlockLength returns a list of the lengths of all characters in a blocks
getBlockLength ∷ V.Vector CharacterData → V.Vector CharInfo → [Int]
getBlockLength inCharDataV inCharInfoV =
    -- trace ("GBL:" <> (show $ V.zipWith U.getCharacterLength inCharDataV inCharInfoV))
    -- False so not use IA field
    V.toList $ V.zipWith (U.getCharacterLength False) inCharDataV inCharInfoV


-- | getBlockNames returns a list of the lengths of all characters in a blocks
getBlockNames ∷ V.Vector CharInfo → [String]
getBlockNames inCharInfoV =
    -- trace ("GBL:" <> (show $ V.zipWith U.getCharacterLength inCharDataV inCharInfoV))
    V.toList $ fmap (T.unpack . name) inCharInfoV


-- | getCharacterStringPair is a wrapper around  getCharacterString
getCharacterStringPair ∷ (CharacterData, CharInfo) → String
getCharacterStringPair (inCharData, inCharInfo) = getCharacterString inCharData inCharInfo


{- | getCharacterString returns a string of character states
need to add space between (large alphabets etc)
local alphabet for characters where that is input.  Matrix and additive are integers
-}
getCharacterString ∷ CharacterData → CharInfo → String
getCharacterString inCharData inCharInfo =
    let inCharType = charType inCharInfo
        localAlphabet =
            if inCharType /= NonAdd
                then ST.toString <$> alphabet inCharInfo
                else fmap ST.toString discreteAlphabet
        -- alphSize =  S.rows $ costMatrix inCharInfo

        -- this to avoid recalculations and list access issues
        lANES = (fromJust $ NE.nonEmpty $ alphabetSymbols localAlphabet)
        lAVect = V.fromList $ NE.toList $ lANES

        getCharState
            ∷ ∀ {b}
             . (Show b, Bits b)
            ⇒ b
            → String
        getCharState a = U.bitVectToCharState localAlphabet lANES lAVect a
    in  let charString = case inCharType of
                x
                    | x == NonAdd →
                        filter (/= ' ') $ foldMap (U.bitVectToCharStateNonAdd localAlphabet) $ snd3 $ stateBVPrelim inCharData
                x | x `elem` packedNonAddTypes → UV.foldMap (U.bitVectToCharStateQual localAlphabet) $ snd3 $ packedNonAddPrelim inCharData
                x | x == Add → filter (/= ' ') $ foldMap (U.additivStateToString lAVect) $ snd3 $ rangePrelim inCharData
                x | x == Matrix → filter (/= ' ') $ foldMap U.matrixStateToString $ matrixStatesPrelim inCharData
                x | x `elem` [SlimSeq, NucSeq] → filter (/= ' ') $ SV.foldMap getCharState $ snd3 $ slimAlignment inCharData
                x | x `elem` [WideSeq, AminoSeq] → filter (/= ' ') $ UV.foldMap getCharState $ snd3 $ wideAlignment inCharData
                x | x == HugeSeq → foldMap getCharState $ snd3 $ hugeAlignment inCharData
                x | x == AlignedSlim → filter (/= ' ') $ SV.foldMap getCharState $ snd3 $ alignedSlimPrelim inCharData
                x | x == AlignedWide → filter (/= ' ') $ UV.foldMap getCharState $ snd3 $ alignedWidePrelim inCharData
                x | x == AlignedHuge → foldMap getCharState $ snd3 $ alignedHugePrelim inCharData
                _ → error ("Un-implemented data type " <> show inCharType)
        in  if not (isAllGaps charString)
                then charString
                else fmap replaceDashWithQuest charString
    where
        replaceDashWithQuest s =
            if s == '-'
                then '?'
                else s


-- nothingToNothing a = if a == '\8220' then '\00'
--                     else a
{-
-- | bitVectToCharStringTNT wraps '[]' around ambiguous states and removes commas between states
bitVectToCharStringTNT ::  (Show b, FiniteBits b, Bits b) =>  Alphabet String -> b -> String
bitVectToCharStringTNT localAlphabet bitValue =
   -- let stateString = U.bitVectToCharState'' (fromJust $ NE.nonEmpty $ alphabetSymbols localAlphabet) bitValue
   let stateString = U.bitVectToCharState' localAlphabet bitValue
   in
   stateString
-}

-- | Implied Alignment report functions

{- | getImpliedAlignmentString returns as a single String the implied alignments of all sequence characters
softwired use display trees, hardWired transform to softwired then proceed with display trees
-}
getImpliedAlignmentString
    ∷ GlobalSettings
    → Bool
    → Bool
    → ProcessedData
    → (ReducedPhylogeneticGraph, Int)
    → PhyG String
getImpliedAlignmentString inGS includeMissing concatSeqs inData (inReducedGraph, graphNumber) =
    if LG.isEmpty (fst5 inReducedGraph)
        then error "No graphs for create IAs for in getImpliedAlignmentStrings"
        else
            let headerString = "Implied Alignments for Graph " <> show graphNumber <> "\n"
                inGraph = GO.convertReduced2PhylogeneticGraph inReducedGraph

                -- parallel stuff
                reoptimize ∷ (ProcessedData, SimpleGraph) → PhyG PhylogeneticGraph
                reoptimize = TRAV.multiTraverseFullyLabelGraphPair (inGS{graphType = Tree}) False False Nothing

                getIAAction ∷ PhylogeneticGraph → PhyG String
                getIAAction = getTreeIAString includeMissing
            in  if graphType inGS == Tree
                    then do
                        resultIA ← getTreeIAString includeMissing inGraph
                        if not concatSeqs
                            then do
                                pure $ headerString <> resultIA
                            else do
                                pure $ headerString <> U.concatFastas resultIA -- (getTreeIAString includeMissing inGraph)
                    else -- for softwired networks--use display trees

                        if graphType inGS == SoftWired
                            then -- get display trees for each data block-- takes first of potentially multiple

                                let blockDisplayList = fmap (GO.contractIn1Out1EdgesRename . GO.convertDecoratedToSimpleGraph . head) (fth6 inGraph)

                                    -- create seprate processed data for each block
                                    blockProcessedDataList = fmap (makeBlockData (fst3 inData) (snd3 inData)) (thd3 inData)
                                in  do
                                        decoratedBlockTreeList' ←
                                            getParallelChunkTraverse >>= \pTraverse →
                                                pTraverse reoptimize . zip (V.toList blockProcessedDataList) $ V.toList blockDisplayList
                                        -- Perform full optimizations on display trees (as trees) with single block data (blockProcessedDataList) to create IAs
                                        let decoratedBlockTreeList = V.fromList decoratedBlockTreeList'

                                        -- extract IA strings as if multiple graphs
                                        diplayIAStringList ← mapM (getTreeIAString includeMissing) (V.toList decoratedBlockTreeList)

                                        if not concatSeqs
                                            then do
                                                pure $ headerString <> concat diplayIAStringList
                                            else do
                                                pure $ headerString <> U.concatFastas (concat diplayIAStringList)
                            else -- There is no IA for Hardwired at least as of yet
                            -- so convert to softwired and use display trees

                                if graphType inGS == HardWired
                                    then
                                        let newGS = inGS{graphType = SoftWired}

                                            pruneEdges = False
                                            warnPruneEdges = False

                                            startVertex ∷ ∀ {a}. Maybe a
                                            startVertex = Nothing
                                        in  do
                                                -- create seprate processed data for each block
                                                newGraph ← TRAV.multiTraverseFullyLabelGraph newGS inData pruneEdges warnPruneEdges startVertex (fst6 inGraph)

                                                let blockDisplayList = fmap (GO.convertDecoratedToSimpleGraph . head) (fth6 newGraph)

                                                let blockProcessedDataList = fmap (makeBlockData (fst3 inData) (snd3 inData)) (thd3 inData)

                                                decoratedBlockTreeList' ←
                                                    getParallelChunkTraverse >>= \pTraverse →
                                                        pTraverse reoptimize . zip (V.toList blockProcessedDataList) $ V.toList blockDisplayList
                                                -- Perform full optimizations on display trees (as trees) with single block data (blockProcessedDataList) to creeate IAs
                                                let decoratedBlockTreeList = V.fromList decoratedBlockTreeList'

                                                -- extract IA strings as if mutiple graphs
                                                diplayIAStringList ←
                                                    getParallelChunkTraverse >>= \pTraverse →
                                                        getIAAction `pTraverse` V.toList decoratedBlockTreeList

                                                logWith
                                                    LogWarn
                                                    "There is no implied alignment for hard-wired graphs--at least not yet. Transfroming to softwired and generate an implied alignment that way\n"
                                                if not concatSeqs
                                                    then do
                                                        pure $ headerString <> concat diplayIAStringList
                                                    else do
                                                        pure $ headerString <> U.concatFastas (concat diplayIAStringList)
                                    else do
                                        logWith LogWarn ("IA  not yet implemented for graphtype " <> show (graphType inGS) <> "\n")
                                        pure $ "There is no implied alignment for " <> show (graphType inGS)


-- | getTreeIAString takes a Tree Decorated Graph and returns Implied AlignmentString
getTreeIAString ∷ Bool → PhylogeneticGraph → PhyG String
getTreeIAString includeMissing inGraph =
    let leafList = snd4 $ LG.splitVertexList (thd6 inGraph)
        leafNameList = fmap (vertName . snd) leafList
        leafDataList = V.fromList $ fmap (vertData . snd) leafList
        charInfoVV = six6 inGraph
    in  do
            characterStringList ← makeFullIAStrings includeMissing charInfoVV leafNameList leafDataList

            pure $ concat characterStringList


-- | makeBlockData cretes new single block processed data
makeBlockData ∷ V.Vector NameText → V.Vector NameBV → BlockData → ProcessedData
makeBlockData a b c = (a, b, V.singleton c)


-- | makeFullIAStrings goes block by block, creating fasta strings for each
makeFullIAStrings ∷ Bool → V.Vector (V.Vector CharInfo) → [NameText] → V.Vector VertexBlockData → PhyG [String]
makeFullIAStrings includeMissing charInfoVV leafNameList leafDataList =
    let numBlocks = V.length charInfoVV
        -- parallel stuff
        action ∷ Int → PhyG [String]
        action = makeBlockIAStrings includeMissing leafNameList leafDataList charInfoVV
    in  getParallelChunkTraverse >>= \pTraverse →
            fold <$> pTraverse action [0 .. numBlocks - 1]


-- | makeBlockIAStrings extracts data for a block (via index) and calls function to make iaStrings for each character
makeBlockIAStrings
    ∷ Bool → [NameText] → V.Vector (V.Vector (V.Vector CharacterData)) → V.Vector (V.Vector CharInfo) → Int → PhyG [String]
makeBlockIAStrings includeMissing leafNameList leafDataList charInfoVV blockIndex =
    let thisBlockCharInfo = V.toList $ charInfoVV V.! blockIndex
        numChars = length thisBlockCharInfo
        thisBlockCharData = fmap (V.! blockIndex) leafDataList

        -- parallel stuff
        action ∷ (CharInfo, Int) → String
        action = makeBlockCharacterStringPair includeMissing leafNameList thisBlockCharData
    in  do
            actionPar ← getParallelChunkMap
            let blockCharacterStringList = actionPar action (zip thisBlockCharInfo [0 .. (numChars - 1)])
            pure blockCharacterStringList


-- blockCharacterStringList = zipWith (makeBlockCharacterStringPair includeMissing leafNameList thisBlockCharData) thisBlockCharInfo [0 .. (numChars - 1)] `using` PU.myParListChunkRDS

-- | isAllGaps checks whether a sequence is all gap charcaters '-'
isAllGaps ∷ String → Bool
isAllGaps inSeq
    | null inSeq = True
    | length (filter (`notElem` [' ', '-', '\n']) inSeq) == 0 = True
    | otherwise = False


-- | makeBlockCharacterStringPair is a wrapper for makeBlockCharacterString
makeBlockCharacterStringPair ∷ Bool → [NameText] → V.Vector (V.Vector CharacterData) → (CharInfo, Int) → String
makeBlockCharacterStringPair includeMissing leafNameList leafDataVV (thisCharInfo, charIndex) =
    makeBlockCharacterString includeMissing leafNameList leafDataVV thisCharInfo charIndex


-- | makeBlockCharacterString creates implied alignmennt string for sequence charactes and null if not
makeBlockCharacterString ∷ Bool → [NameText] → V.Vector (V.Vector CharacterData) → CharInfo → Int → String
makeBlockCharacterString includeMissing leafNameList leafDataVV thisCharInfo charIndex =
    -- check if sequence type character
    let thisCharType = charType thisCharInfo
        thisCharName = name thisCharInfo
    in  if thisCharType `notElem` sequenceCharacterTypes
            then []
            else
                let -- thisCharData = fmap (V.! charIndex) leafDataVV
                    thisCharData = getTaxDataOrMissing leafDataVV charIndex 0 []
                    nameDataPairList = zip leafNameList thisCharData
                    fastaString = pairList2Fasta includeMissing thisCharInfo nameDataPairList
                in  -- trace ("MBCS: " <> (show $ length leafNameList) <> " " <> (show $ V.length thisCharData) <> "\n" <> (show leafDataVV))
                    "\nSequence character " <> T.unpack thisCharName <> "\n" <> fastaString <> "\n"


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

-- | getTaxDataOrMissing gets the index character if data not null, empty character if not
getTaxDataOrMissing ∷ V.Vector (V.Vector CharacterData) → Int → Int → [CharacterData] → [CharacterData]
getTaxDataOrMissing charDataV charIndex taxonIndex newTaxList
    | taxonIndex == V.length charDataV = reverse newTaxList
    | V.null (charDataV V.! taxonIndex) = getTaxDataOrMissing charDataV charIndex (taxonIndex + 1) (emptyCharacter : newTaxList)
    | otherwise = getTaxDataOrMissing charDataV charIndex (taxonIndex + 1) (((charDataV V.! taxonIndex) V.! charIndex) : newTaxList)


{- | pairList2Fasta takes a character type and list of pairs of taxon names (as T.Text)
and character data and returns fasta formated string
-}
pairList2Fasta ∷ Bool → CharInfo → [(NameText, CharacterData)] → String
pairList2Fasta includeMissing inCharInfo nameDataPairList =
    if null nameDataPairList
        then []
        else
            let (firstName, blockDatum) = head nameDataPairList
                inCharType = charType inCharInfo

                -- this to avoid recalculations and list access issues
                localAlphabet = (ST.toString <$> alphabet inCharInfo)
                lANES = (fromJust $ NE.nonEmpty $ alphabetSymbols localAlphabet)
                lAVect = V.fromList $ NE.toList $ lANES

                getCharState
                    ∷ ∀ {b}
                     . (Show b, Bits b)
                    ⇒ b
                    → String
                getCharState a = U.bitVectToCharState localAlphabet lANES lAVect a

                sequenceString = case inCharType of
                    -- x | x `elem` [SlimSeq, NucSeq  ] -> SV.foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ slimAlignment blockDatum
                    x | x == SlimSeq → SV.foldMap getCharState $ snd3 $ slimAlignment blockDatum
                    x | x == NucSeq → SV.foldMap getCharState $ snd3 $ slimAlignment blockDatum
                    -- x | x `elem` [WideSeq, AminoSeq] -> UV.foldMap (U.bitVectToCharState localAlphabet) $ snd3 $ wideAlignment blockDatum
                    x | x == WideSeq → UV.foldMap getCharState $ snd3 $ wideAlignment blockDatum
                    x | x == AminoSeq → UV.foldMap getCharState $ snd3 $ wideAlignment blockDatum
                    x | x == HugeSeq → foldMap getCharState $ snd3 $ hugeAlignment blockDatum
                    x | x == AlignedSlim → SV.foldMap getCharState $ snd3 $ alignedSlimPrelim blockDatum
                    x | x == AlignedWide → UV.foldMap getCharState $ snd3 $ alignedWidePrelim blockDatum
                    x | x == AlignedHuge → foldMap getCharState $ snd3 $ alignedHugePrelim blockDatum
                    _ → error ("Un-implemented data type " <> show inCharType)

                -- If all gaps then change to all question marks since missing
                sequenceString' =
                    if isAllGaps sequenceString
                        then fmap replaceDashWithQuest sequenceString
                        else sequenceString

                -- Make lines 50 chars long
                sequenceChunks = ((<> "\n") <$> LS.chunksOf 50 sequenceString')
            in  if ((not includeMissing) && (isAllGaps sequenceString)) || (blockDatum == emptyCharacter)
                    then pairList2Fasta includeMissing inCharInfo (tail nameDataPairList)
                    else
                        (concat $ (('>' : (T.unpack firstName)) <> "\n") : sequenceChunks)
                            <> (pairList2Fasta includeMissing inCharInfo (tail nameDataPairList))
    where
        replaceDashWithQuest s =
            if s == '-'
                then '?'
                else s
