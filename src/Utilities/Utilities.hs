{- |
Module specifying utility functions for use with PhyGraph
-}
module Utilities.Utilities where

import Bio.DynamicCharacter.Element (SlimState, WideState, HugeState)
import Complexity.Graphs qualified as GC
import Complexity.Utilities qualified as GCU
import Control.DeepSeq (force)
import Control.Monad (replicateM)
import Control.Monad.IO.Class (MonadIO (..))
import Control.Monad.Random.Class
import Data.Alphabet
import Data.Alphabet.IUPAC
import Data.Alphabet.Special
-- import Data.Alphabet.Special
import Data.Bimap qualified as BM
import Data.BitVector.LittleEndian qualified as BV
import Data.Bits
import Data.Foldable
import Data.Functor ((<&>))
import Data.InfList qualified as IL
import Data.List qualified as L
import Data.List.NonEmpty (NonEmpty (..))
import Data.List.NonEmpty qualified as NE
import Data.List.Split qualified as SL
import Data.Maybe
import Data.Set qualified as SET
import Data.Text.Lazy qualified as T
import Data.Text.Short qualified as ST
import Data.Vector qualified as V
import Data.Vector.Storable qualified as SV
import Data.Vector.Unboxed qualified as UV
import Debug.Trace
import GeneralUtilities
import GeneralUtilities qualified as GU
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Verbosity (Verbosity (..))
import SymMatrix qualified as S
import Types.Types
import Utilities.LocalGraph qualified as LG
import GraphOptimization.Medians (get2WaySlim, get2WayWideHuge,)


{- | strict2of5 esures parallelism and demands strict return of 2nd of 5 tuple elements
    this is used in lazy-ish parallel evalution functions in PHANE evaluation
-}
strict2of5 :: ReducedPhylogeneticGraph -> ReducedPhylogeneticGraph
strict2of5 b@(_,!x,_,_,_) =
        force x `seq` b

{- | needTwoEdgeNoCostAdjust checks global data for PMDL or SI
and whether the required median is a distance (ie single edge)
or two edge median (as in creating a vertex for post order traversal)
and returns a boolean is adjustment is required.
this to add in the extra "noChange" costs when required
-}
needTwoEdgeNoCostAdjust ∷ GlobalSettings → Bool → Bool
needTwoEdgeNoCostAdjust inGS isTwoEdgeMedian =
    if not isTwoEdgeMedian
        then False
        else
            if optimalityCriterion inGS `notElem` [SI, PMDL, MAPA]
                then False
                else True


{- | collapseGraph collapses zero-length edges in 3rd field of a phylogenetic graph
does not affect display trees or character graphs
fst6 and thd6 (both) are modified since this is used for output
Collapsed frst can no longer be used for graph optimization since non-binary
this is done by removing the end vertex 'v' of min length 0 edge (u->v)
this removes (u,v) and the two (v, x) and (v, y) out edges from v.  New edges are created
(u, x) and (u,y) with length labels of (v,x) and (v,y)
assumes all indexing is the same between the simple and decorated graph
done recusively until no minLength == zero edges so edges renumbered properly
network edges, pendant edges and root edges, are not collapsed
this wierd thype is to allow for polymorphism in graph type--basically a general phylogenetic graph
-}
collapseGraph
    ∷ GenPhyloGraph a b
    → GenPhyloGraph a b
collapseGraph inPhylograph@(inSimple, inC, inDecorated, inD, inE, inF) =
    if LG.isEmpty inSimple
        then inPhylograph
        else
            let inDecTreeEdgeList = filter (not . LG.isNetworkLabEdge inDecorated) $ LG.labEdges inDecorated
                zeroEdgeList' = filter ((== 0.0) . minLength . thd3) inDecTreeEdgeList

                -- remove cases of pendent edges--don't remove those either
                leafIndexList = fmap fst $ snd4 $ LG.splitVertexList inSimple
                rootChildIndexList = fmap snd3 $ concatMap (LG.out inSimple) $ fmap fst $ fst4 $ LG.splitVertexList inSimple
                zeroEdgeList = filter ((`notElem` (leafIndexList <> rootChildIndexList)) . snd3) zeroEdgeList'
            in  if null zeroEdgeList
                    then inPhylograph
                    else -- get node to be deleted and its out edges

                        let nodeToDelete = snd3 $ head zeroEdgeList
                            sourceEdgeToDelete = fst3 $ head zeroEdgeList

                            -- new dec edges
                            firstOutEdgeDec = head $ LG.out inDecorated nodeToDelete
                            secondOutEdgeDec = last $ LG.out inDecorated nodeToDelete

                            newFirstEdgeDec = (sourceEdgeToDelete, snd3 firstOutEdgeDec, thd3 firstOutEdgeDec)
                            newSecondEdgeDec = (sourceEdgeToDelete, snd3 secondOutEdgeDec, thd3 secondOutEdgeDec)

                            -- new simple edges
                            firstOutEdgeSimple = head $ LG.out inSimple nodeToDelete
                            secondOutEdgeSimple = last $ LG.out inSimple nodeToDelete

                            newFirstEdgeSimple = (sourceEdgeToDelete, snd3 firstOutEdgeSimple, thd3 firstOutEdgeSimple)
                            newSecondEdgeSimple = (sourceEdgeToDelete, snd3 secondOutEdgeSimple, thd3 secondOutEdgeSimple)

                            -- make new decorated--deleting node removes all incident edges
                            newDecGraph = LG.insEdges [newFirstEdgeDec, newSecondEdgeDec] $ LG.delNode nodeToDelete inDecorated

                            -- make new simple--deleting node removes all incident edges
                            newSimpleGraph = LG.insEdges [newFirstEdgeSimple, newSecondEdgeSimple] $ LG.delNode nodeToDelete inSimple
                        in  (newSimpleGraph, inC, newDecGraph, inD, inE, inF)


-- | collapseReducedGraph is a wrpper to collapseGraph
collapseReducedGraph ∷ ReducedPhylogeneticGraph → ReducedPhylogeneticGraph
collapseReducedGraph (inSimple, inC, inDecorated, inD, inF) =
    let (newSimpleGraph, _, newDecGraph, _, _, _) = collapseGraph (inSimple, inC, inDecorated, mempty, mempty, inF)
    in  (newSimpleGraph, inC, newDecGraph, inD, inF)


{- | calculateGraphComplexity returns an infinite list of graph complexities indexed by
number of network nodes-assumes for now--single component graph not forest
first in pair is softwired complexity, second hardwired complexity
could coppy to vector for say 100 or so and offset infinite list after
to reduce lisyt access facvtor
-}
calculateGraphComplexity ∷ ProcessedData → IL.InfList (VertexCost, VertexCost)
calculateGraphComplexity (nameVect, _, blockDatVect) =
    let numNetNodesList = IL.fromList [(0 ∷ Int) ..]
        numRoots = 1
        numBlocks = V.length blockDatVect
        graphComplexity = IL.map (getGraphComplexity (V.length nameVect) numRoots numBlocks) numNetNodesList
    in  graphComplexity


{- | old version with direct bit calcualtions
| getGraphComplexity' takes the number of leaves and number of
network nodes and calculates the graph complexity in bits
tree num edges (2n-2) n leaves * 2 nodes for each edge * (log 2n -1 vertices-- min specify)
-}
getGraphComplexity' ∷ Int → Int → Int → (VertexCost, VertexCost)
getGraphComplexity' numLeaves numRoots numNetNodes =
    -- place holder for now
    let nodeComplexity = logBase 2.0 (fromIntegral $ (2 * numLeaves) - 1 + numNetNodes) -- bits to specify each vertex
        treeEdges = (2 * numLeaves) - 2
        extraRootEdges = 2 * (numRoots - 1)
        baseTreeComplexity = nodeComplexity * fromIntegral (2 * (treeEdges - extraRootEdges))
        numDisplayTrees = 2.0 ** fromIntegral numNetNodes
        harWiredEdges = (treeEdges - extraRootEdges) + (3 * numNetNodes)
        hardwiredAddComplexity = nodeComplexity * fromIntegral (2 * harWiredEdges)
    in  -- maybe softwired is numDisplatTrees * harWired since have those edges in input
        (baseTreeComplexity * numDisplayTrees, hardwiredAddComplexity)


{- | getGraphComplexity takes the number of leaves and number of
network nodes and calculates the algorithmic graph complexity in bits
-}
getGraphComplexity ∷ Int → Int → Int → Int → (VertexCost, VertexCost)
getGraphComplexity numLeaves numRoots numBlocks numNetNodes =
    -- place holder for now
    let graphProgram = GC.makeProgramStringGraph numLeaves 0 numRoots numNetNodes
        (_, _, _, gzipGraph) = GCU.getInformationContent graphProgram

        graphDisplayProgram = GC.makeDisplayGraphString numLeaves 0 numRoots numNetNodes
        (_, _, _, gzipDisplay) = GCU.getInformationContent graphDisplayProgram

        displayTreeSwitchingComplexity = fromIntegral numNetNodes
        marginalDisplayComplexity = gzipDisplay - gzipGraph -- graphDisplayShannonBits - graphShannonBits

        -- cost of swithing (speciying) 1 bit per netNode then minimum of blocks as duspolay tree number since only have a few block usually
        softWiredFactor = displayTreeSwitchingComplexity + ((min (2 ** fromIntegral numNetNodes) (fromIntegral numBlocks)) * marginalDisplayComplexity)
    in  (gzipGraph + softWiredFactor, gzipGraph)


{- | calculateMAPARootCost-- for  now used NCM--but better to reflect empirical Pi (frequency) values
won't affect the search choice since a constant factor
-}
calculateMAPARootCost ∷ ProcessedData → VertexCost
calculateMAPARootCost = calculateNCMRootCost


{- | calculateNCMRootCost calcuates the contant fact of SUM -log 10 1/r over all characters
approximate for packed data (based on alphabet size for packed)
-}
calculateNCMRootCost ∷ ProcessedData → VertexCost
calculateNCMRootCost (_, _, blockDataV) =
    if V.null blockDataV
        then 0.0
        else V.sum $ fmap getBlockNCMRootCost blockDataV


-- | getBlockNCMRootCost gets NCM root cost for character block
getBlockNCMRootCost ∷ BlockData → VertexCost
getBlockNCMRootCost (_, charDataVV, charInfoV) =
    if V.null charDataVV || V.null charInfoV
        then 0
        else -- get length of each characters
        -- this for prealigned and non-aligned sequences mainly
        -- but if data are reorganized and packed--all data

            let numChars = V.length charInfoV
                leafCharListV = fmap (charDataVV V.!) [0 .. numChars - 1]
                -- False fo not use IA field
                maxCharLengthList = zipWith (getMaxCharacterLength False) (V.toList charInfoV) (fmap V.toList leafCharListV)
                weightList = fmap weight (V.toList charInfoV)
                rootCostList = zipWith (*) weightList (fmap fromIntegral maxCharLengthList)
            in  -- trace ("GNCMR: " <> (show (numChars, maxCharLengthList, weightList, rootCostList))) $
                sum rootCostList


{- | calculatePMDLVertexComplexity creates a vertex cost as either the 'insertion' of character data
or as probbaility of seqeunce in bit based on element frequencies.  
-}
calculatePMDLVertexComplexity ∷ Bool -> Maybe DecoratedGraph -> ProcessedData -> Maybe Int -> VertexCost
calculatePMDLVertexComplexity useLogPiValues decGraph (nameVect, _, blockDataV) index =
    --trace ("In CPMDLRC") $
    if useLogPiValues then 
        -- root complexity base on log2 Pis
        --trace ("New way: " <> (show (getLogPiRootCost blockDataV))) $
        getLogPiRootCost blockDataV

    else if isNothing index then --False so insert--if index /= Nothing 
        -- average of leaves
        -- use insert based (but pi of '-' prob way underestimated)
        let numLeaves = V.length nameVect
            insertDataCost = V.sum $ fmap getblockInsertDataCost blockDataV
        in  
        insertDataCost / fromIntegral numLeaves
    else -- optimized graph specific vertex cost-- not averager
        error "Optimized data complexity noy yet implemented"

{- | calculatePMDLRootCost creates a root cost as either the 'insertion' of character data
or as probbaility of seqeunce in bit based on element frequencies.  
For sequence data averaged over leaf taxa
so root independent

Wrapper around more flexible calculatePMDLVertexComplexity
-}
calculatePMDLRootCost ∷ Bool -> ProcessedData -> VertexCost
calculatePMDLRootCost useLogPiValues inData =
    calculatePMDLVertexComplexity useLogPiValues Nothing inData Nothing
   


{- | getblockInsertDataCost gets the total cost of 'inserting' the data in a block
this most easily done before bit packing since won't vary anyway.
then store value in Global Settings
-}
getblockInsertDataCost ∷ BlockData → Double
getblockInsertDataCost (_, characterDataVV, charInfoV) =
    V.sum $ fmap (getLeafInsertCost charInfoV) characterDataVV


{- | getLeafInsertCost is the cost or originating or 'inserting' leaf data
for all characters in a block
-}
getLeafInsertCost ∷ V.Vector CharInfo → V.Vector CharacterData → Double
getLeafInsertCost charInfoV charDataV =
    V.sum $ V.zipWith getCharacterInsertCost charDataV charInfoV

{- | getLogPiRootCost gets log 2 ofelement dreqnecies over all data blocks
-}
getLogPiRootCost :: V.Vector BlockData -> Double 
getLogPiRootCost inBlockDataV = 
    --trace ("In GLPRC: " <> (show (V.sum $ fmap getBlockLogPiCost inBlockDataV))) $
    V.sum $ fmap getBlockLogPiCost inBlockDataV

{- | getBlockLogPiCost sums up the character root complexity of a block
relies on pulling charcters out of block
-}
getBlockLogPiCost :: BlockData -> Double
getBlockLogPiCost (_, leafCharVV, charInfoV) = 
    if V.null charInfoV then 0.0
    else
        -- create leaf by single character structure to be mapped over 
        let charNumList = [0.. (V.length $ V.head leafCharVV) - 1]
            charVect = fmap (leafCharVV V.!) (V.fromList charNumList)
        in
        --trace ("In GBLPC: " <> (show ((V.length $ V.head leafCharVV),  V.sum $ V.zipWith getCharacterLogPiCost charVect charInfoV))) $
        V.sum $ V.zipWith getCharacterLogPiCost charVect charInfoV

{- | getCharacterLogPiCost takes a characterdata (vector of terminals for single character) 
and returns root complexity cost based on frequency of 
sequence character elements and state number for no-sequence 
-}
getCharacterLogPiCost :: V.Vector CharacterData -> CharInfo -> Double
getCharacterLogPiCost charDataV charInfo =
    if V.null charDataV then 0.0
    else 
        let nonSequenceCharRootCost = getNonSequenceCharacterLogPiCost (V.head charDataV) charInfo
            sequenceCharRootCost = getSequenceCharacterLogPiCost charDataV charInfo
        in
        --trace ("In GCLPC: " <> (show (nonSequenceCharRootCost, sequenceCharRootCost))) $
        nonSequenceCharRootCost + sequenceCharRootCost

{- | getNonSequenceCharacterLogPiCost creates root cost for all add/nonAdd/Matrix characters
based number of character states and the log2 of that same for each leaf (since static)
so no need to average--can use a single leaf
-}
getNonSequenceCharacterLogPiCost ∷ CharacterData → CharInfo → Double
getNonSequenceCharacterLogPiCost inChar charInfo =
    -- trace ("In GCIC: " <> (show $ charType charInfo)) $
    let localCharType = charType charInfo
        
        numStates =
            if localCharType == NonAdd
                then BV.dimension $ (V.head $ GU.snd3 $ stateBVPrelim inChar)
                else
                    if localCharType == Add
                        then
                            let (a, b) = V.head $ GU.snd3 $ rangePrelim inChar
                            in  toEnum (b - a)
                        else
                            if localCharType == Matrix
                                then toEnum $ V.length $ V.head $ matrixStatesPrelim inChar
                                else 0 ∷ Word
        alphabetWeight = logBase 2.0 $ (fromIntegral $ numStates ∷ Double)

        
        insertCost  
            | localCharType == Add                      = alphabetWeight * fromIntegral (V.length $ GU.snd3 $ rangePrelim inChar)
            | localCharType == NonAdd                   = alphabetWeight * fromIntegral (V.length $ GU.snd3 $ stateBVPrelim inChar)
            | localCharType `elem` packedNonAddTypes    = fromIntegral (UV.length $ GU.snd3 $ packedNonAddPrelim inChar)
            | localCharType == Matrix                   = alphabetWeight * fromIntegral (V.length $ matrixStatesPrelim inChar)
            | otherwise = error ("Non static charcter type: " <> (show localCharType))

    in  --trace ("GCIC: " <> (show (insertCost))) $
    if localCharType `notElem` ([Add, NonAdd, Matrix] <> packedNonAddTypes) then 0.0
    else insertCost

{- | getSequenceCharacterLogPiCost creates root cost for all characters as with add/nonAdd
based on the frequencies of character states and the log2 of that frequency for each state or sequence element
average cost over leaves so root placement independent
-}
getSequenceCharacterLogPiCost :: V.Vector CharacterData -> CharInfo -> Double
getSequenceCharacterLogPiCost charDataV charInfo =
    --trace ("In GSCLPC") $
    if V.null charDataV then error "No data to calculate root cost"
    else if charType charInfo `notElem` sequenceCharacterTypes then 0.0
    else 
        let isNeyman = checkNeyman (costMatrix charInfo)
            numLeaves = V.length charDataV
            leafList = V.toList $ fmap (getSequenceLeafCharStateNumbers charInfo) charDataV
            elementList = fmap fst $ head leafList
            numberList = fmap (getElementNumberList (concat leafList)) elementList
            elementNumList = zip elementList numberList
            totalElements = fromIntegral $ sum $ fmap snd elementNumList
            elemFreqList = fmap (/ totalElements) $ fmap fromIntegral $ fmap snd elementNumList
            bitList = if isNeyman then
                        replicate (length $ alphabet charInfo) $ logBase 2.0 $ (fromIntegral $ (length $ alphabet charInfo) ∷ Double)
                      else fmap (hereLogBase2) elemFreqList
        in
        -- trace ("GSCLPC: " <> (show (isNeyman,elemFreqList, numberList, bitList))) $
        abs $ (sum $ zipWith (*) (fmap fromIntegral numberList) bitList) / (fromIntegral numLeaves)

        where hereLogBase2 a = if a < epsilon then 0.0
                               else logBase 2.0 a

{- | checkNeyman take a matrix (type from PHANE) and returns True if
    the matrix has all non-diagonal values the same
    otherwise False
-}
checkNeyman :: S.Matrix Int -> Bool
checkNeyman inMatrix =
    if S.null inMatrix then False
    else if S.rows inMatrix < 2 then False
    else 
        let nonDiagVal = inMatrix S.! (0,1)
        in
        checkMatrixVals inMatrix nonDiagVal 0 0 

{- | checkMatrixVals checks of all non-adiagnonal values equal input value
-} 
checkMatrixVals :: S.Matrix Int -> Int -> Int -> Int -> Bool
checkMatrixVals inMatrix nonDiagVal rowIndex columnIndex =
    --trace ("Check Matrix: " <> (show (rowIndex, columnIndex, S.rows inMatrix, S.cols inMatrix))) $
    if rowIndex >= S.rows inMatrix then True
    else if columnIndex >= S.cols inMatrix then checkMatrixVals inMatrix nonDiagVal (rowIndex + 1) 0 
    else if rowIndex == columnIndex then checkMatrixVals inMatrix nonDiagVal rowIndex (columnIndex + 1)
    else if inMatrix S.! (rowIndex, columnIndex) /= nonDiagVal then False
    else checkMatrixVals inMatrix nonDiagVal rowIndex (columnIndex + 1)

{- | getCharacterElementFreq gets the alphabet elements of sequence characters in a 
specific character over all leaves
-}
getCharacterElementFreq :: V.Vector CharacterData -> CharInfo -> [(String, Double)]
getCharacterElementFreq charDataV charInfo =
    if V.null charDataV then []
    else 
        let leafList = V.toList $ fmap (getSequenceLeafCharStateNumbers charInfo) charDataV
            elementList = fmap fst $ head leafList
            elementNumList = zip elementList (fmap (getElementNumberList (concat leafList)) elementList)
            totalElements = fromIntegral $ sum $ fmap snd elementNumList
            elemFreqList = fmap (/ totalElements) $ fmap fromIntegral $ fmap snd elementNumList
        in
        zip elementList elemFreqList

{- |  getElementNumberList takes a list of (String, Int) and returns 
the summed Ints based on lisyt of String
-}
getElementNumberList :: [(String, Int)] -> String -> Int 
getElementNumberList pairList matchSymbol =
    if null pairList then 0
    else 
        if (fst $ head pairList) == matchSymbol then
            (snd $ head pairList) + getElementNumberList (tail pairList) matchSymbol
        else getElementNumberList (tail pairList) matchSymbol

{- | getSequenceLeafCharStateNumbers gets the numebr of each charcter state as in alphabet
-}
getSequenceLeafCharStateNumbers :: CharInfo -> CharacterData → [(String, Int)]
getSequenceLeafCharStateNumbers charInfo blockDatum =
    let localType = charType charInfo
        localAlphabet = (ST.toString <$> alphabet charInfo)

        -- this to avoid recalculations and list access issues
        lANES = (fromJust $ NE.nonEmpty $ alphabetSymbols localAlphabet)
        lAList = NE.toList $ lANES
        lAVect = V.fromList $ NE.toList $ lANES

        getCharState
            ∷ ∀ {b}
             . (Show b, Bits b)
            ⇒ b
            → String
        getCharState a = (bitVectToCharState localAlphabet lANES lAVect a) <> ","

        -- get strings (with ",") of elements
        stringStates
            | localType `elem` sequenceCharacterTypes = case localType of
                x | x `elem` [NucSeq] → (SV.foldMap getCharState (slimPrelim blockDatum))
                x | x `elem` [SlimSeq] → (SV.foldMap getCharState (slimPrelim blockDatum))
                x | x `elem` [WideSeq] → (UV.foldMap getCharState (widePrelim blockDatum))
                x | x `elem` [AminoSeq] → (UV.foldMap getCharState (widePrelim blockDatum))
                x | x `elem` [HugeSeq] → (foldMap getCharState (hugePrelim blockDatum))
                x | x `elem` [AlignedSlim] → (SV.foldMap getCharState $ snd3 $ alignedSlimPrelim blockDatum)
                x | x `elem` [AlignedWide] → (UV.foldMap getCharState $ snd3 $ alignedWidePrelim blockDatum)
                x | x `elem` [AlignedHuge] → (foldMap getCharState $ snd3 $ alignedHugePrelim blockDatum)
                _ → error ("Un-implemented sequence data type " <> show localType)
            | otherwise = error ("Non-sequence data type " <> show localType)

        -- filter out ambiguities (ie not in element list)
        elementListList = filter (`elem` lAList) $  SL.splitOn "," stringStates
        elementList = L.group $ L.sort elementListList
        elementNumberList = fmap (getElementNumber elementList) lAList

    in
    zip lAList elementNumberList


{- | getElementNumber takes a list of list of elements that has been grouped and checkes
if matches String and if such returns length else goes on to next input String
-}
getElementNumber :: [[String]] -> String -> Int
getElementNumber symbolLL matchSymbol =
    if null symbolLL then 0
    else 
        if (head $ head symbolLL) == matchSymbol then length $ head symbolLL
        else getElementNumber (tail symbolLL) matchSymbol 


{- | getCharacterInsertCost takes a character and characterInfo and returns origination/insert cost for the character
for PMDL of add, non add matrix log2 of alphabet size
This uses actual inserts of different elements as oposedd to average, although averaged over leaves so not 
root dependant
-}
getCharacterInsertCost ∷ CharacterData → CharInfo → Double
getCharacterInsertCost inChar charInfo =
    --trace ("In GCIC: " <> (show $ charType charInfo)) $
    let localCharType = charType charInfo
        thisWeight = weight charInfo
        
        numStates =
            if localCharType == NonAdd
                then BV.dimension $ (V.head $ GU.snd3 $ stateBVPrelim inChar)
                else
                    if localCharType == Add
                        then
                            let (a, b) = V.head $ GU.snd3 $ rangePrelim inChar
                            in  toEnum (b - a)
                        else
                            if localCharType == Matrix
                                then toEnum $ V.length $ V.head $ matrixStatesPrelim inChar
                                else 0 ∷ Word
        alphabetWeight = logBase 2.0 $ (fromIntegral $ numStates ∷ Double)

        slimLength = max (SV.length $ slimPrelim inChar) (SV.length $ snd3 $ alignedSlimPrelim inChar)
        wideLength = max (UV.length $ widePrelim inChar) (UV.length $ snd3 $ alignedWidePrelim inChar)

        allGapSlim = SV.replicate slimLength $ (0 ∷ SlimState) `setBit` fromEnum gapIndex
        allGapWide = UV.replicate wideLength $ (0 ∷ WideState) `setBit` fromEnum gapIndex

        hugeGapChar = ((V.head $ hugePrelim inChar) `xor` (V.head $ hugePrelim inChar)) `setBit` fromEnum gapIndex 
        allGapHuge = V.replicate (V.length $ hugePrelim inChar) hugeGapChar

        hugeGapCharAligned = ((V.head $ snd3 $ alignedHugePrelim inChar) `xor` (V.head $ snd3 $ alignedHugePrelim inChar)) `setBit` fromEnum gapIndex 
        allGapHugeAligned = V.replicate (V.length $ snd3 $ alignedHugePrelim inChar) hugeGapCharAligned
     
        -- don't need to worry about noCost Adjustment and medioan issue since "-" -> "-" is zero
        -- otherwise would hahve to subtract off extra no-change cost (gap to gap)
        insertCost  
            | localCharType == Add                      = alphabetWeight * fromIntegral (V.length $ GU.snd3 $ rangePrelim inChar)
            | localCharType == NonAdd                   = alphabetWeight * fromIntegral (V.length $ GU.snd3 $ stateBVPrelim inChar)
            | localCharType `elem` packedNonAddTypes    = fromIntegral (UV.length $ GU.snd3 $ packedNonAddPrelim inChar)
            | localCharType == Matrix                   = alphabetWeight * fromIntegral (V.length $ matrixStatesPrelim inChar)
            | localCharType == SlimSeq                  = fromIntegral $ snd $ get2WaySlim (slimTCM charInfo) allGapSlim (slimPrelim inChar)
            | localCharType == NucSeq                   = fromIntegral $ snd $ get2WaySlim (slimTCM charInfo) allGapSlim (slimPrelim inChar) 
            | localCharType == WideSeq                  = fromIntegral $ snd $ get2WayWideHuge (wideTCM charInfo) allGapWide (widePrelim inChar) 
            | localCharType == AminoSeq                 = fromIntegral $ snd $ get2WayWideHuge (wideTCM charInfo) allGapWide (widePrelim inChar)
            | localCharType == HugeSeq                  = fromIntegral $ snd $ get2WayWideHuge (hugeTCM charInfo) allGapHuge (hugePrelim inChar)
            | localCharType == AlignedSlim              = fromIntegral $ snd $ get2WaySlim (slimTCM charInfo) allGapSlim (snd3 $ alignedSlimPrelim inChar)
            | localCharType == AlignedWide              = fromIntegral $ snd $ get2WayWideHuge (wideTCM charInfo) allGapWide (snd3 $ alignedWidePrelim inChar)
            | localCharType == AlignedHuge              = fromIntegral $ snd $ get2WayWideHuge (hugeTCM charInfo) allGapHugeAligned (snd3 $ alignedHugePrelim inChar)
            | otherwise = error ("Character type unimplemented : " <> show localCharType)

    in  --trace ("GCIC: " <> (show (insertCost, thisWeight * insertCost, allGapSlim, snd3 $ alignedSlimPrelim inChar))) $
        thisWeight * insertCost

  
{- | splitSequence takes a ShortText divider and splits a list of ShortText on
that ShortText divider consuming it akin to Text.splitOn
-}
splitSequence ∷ ST.ShortText → [ST.ShortText] → [[ST.ShortText]]
splitSequence partitionST stList =
    if null stList
        then []
        else
            let firstPart = takeWhile (/= partitionST) stList
                restList = dropWhile (/= partitionST) stList
            in  if restList == [partitionST]
                    then firstPart : [[ST.fromString "#"]]
                    else
                        if not $ null restList
                            then firstPart : splitSequence partitionST (tail restList)
                            else [firstPart]


-- See Bio.DynamicCharacter.decodeState for a better implementation for dynamic character elements
bitVectToCharStateQual ∷ (Show b, FiniteBits b, Bits b) ⇒ Alphabet String → b → String
bitVectToCharStateQual localAlphabet bitValue =
    let charString = L.intercalate "," $ foldr pollSymbol mempty indices
    in  if popCount bitValue == finiteBitSize bitValue
            then "?"
            else
                if popCount bitValue > 1
                    then "[" <> charString <> "]"
                    else charString
    where
        indices = [0 .. len - 1]
        len = length vec
        -- this is a hack--the alphabets for non-additive charcaters gets truncated to binary at some point earlier
        vec = V.fromList $ fmap show [0 .. finiteBitSize bitValue - 1]
        pollSymbol i polled
            | bitValue `testBit` i = (vec V.! i) : polled
            | otherwise = polled


-- See Bio.DynamicCharacter.decodeState for a better implementation for dynamic character elements
-- this for TNT output of qualitative characters
bitVectToCharStateNonAdd ∷ (Show b, FiniteBits b, Bits b) ⇒ Alphabet String → b → String
bitVectToCharStateNonAdd localAlphabet bitValue =
    let stateList = [0 .. (finiteBitSize bitValue) - 1]
        stateCharList = fmap (: []) $ ['0' .. '9'] <> ['A' .. 'Z'] <> ['a' .. 'z']
        bitOnList = fmap (testBit bitValue) stateList
        statesON = fmap fst $ filter ((== True) . snd) $ zip stateCharList bitOnList
        charString = concat statesON
    in  -- trace ("BVNA: " <> (show (bitValue, bitOnList, charString))) $
        if popCount bitValue == finiteBitSize bitValue
            then "?"
            else
                if popCount bitValue > 1
                    then "[" <> charString <> "]"
                    else charString


-- See Bio.DynamicCharacter.decodeState for a better implementation for dynamic character elements
bitVectToCharState' ∷ (FiniteBits b, Bits b) ⇒ Alphabet String → b → String
bitVectToCharState' localAlphabet bitValue =
    -- check for symbol length > 1 then add space (since sorted last element longest)
    let -- maxSymbolLength = maximum (length <$> SET.toList (alphabetSymbols localAlphabet))
        charString = foldr pollSymbol mempty indices
        charString' = L.intercalate "," $ filter (/= "\8220") charString
    in  -- trace ("BV2CSA:" <> (show (maxSymbolLength, SET.toList (alphabetSymbols localAlphabet) ))) (
        if popCount bitValue == finiteBitSize bitValue
            then "?"
            else
                if popCount bitValue > 1
                    then "[" <> charString' <> "]" <> " "
                    else charString' <> " "
    where
        -- )

        indices = [0 .. len - 1]
        len = length vec
        vec = alphabetSymbols localAlphabet
        pollSymbol i polled
            | bitValue `testBit` i = (vec V.! i) : polled
            | otherwise = polled


bitVectToCharState ∷ (Show b, Bits b) ⇒ Alphabet String → NonEmpty String → V.Vector String → b → String
bitVectToCharState localAlphabet localAlphabetNEString localAlphabetVect bitValue =
    -- trace ("BVCA': " <> (show $ isAlphabetAminoAcid localAlphabet) <> " " <> (show $ isAlphabetDna localAlphabet) <> " " <> (show localAlphabet)) (
    let stringVal' = foldr pollSymbol mempty indices
        stringVal = concat stringVal'
    in  {-
        if length stringVal == 1 then L.intercalate "," stringVal' <> " "
        else

          -- if isAlphabetAminoAcid localAlphabet then
          if SET.size (alphabetSymbols localAlphabet) > 5 then

              if stringVal == "DN" then "B" <> " "
              else if stringVal == "EQ" then "Z" <> " "
              else if stringVal == "ACDEFGHIKLMNPQRSTVWY" then "X" <> " "
              else if stringVal == "-ACDEFGHIKLMNPQRSTVWY" then "?" <> " "
               -- amino acid polymorphisms without ambiguity codes
              else "[" <> stringVal <> "]" <> " "

          -- Nucleotide IUPAC
          -- hack until fix isDNA and RNA alphabet
          -- else if ((show localAlphabet) == ("Alphabet: {\"-\", \"A\", \"C\", \"G\", \"T\"}")) || ((show localAlphabet) == ("Alphabet: {\"-\", \"A\", \"C\", \"G\", \"U\"}")) then
          -}

        if (isAlphabetDna localAlphabet || isAlphabetRna localAlphabet) && (SET.size (alphabetSymbols localAlphabet) == 5)
            then
                if stringVal `elem` ["", "-"]
                    then "-"
                    else
                        if length stringVal == 1
                            then stringVal
                            else
                                if stringVal == "AG"
                                    then "R"
                                    else
                                        if stringVal == "CT"
                                            then "Y"
                                            else
                                                if stringVal == "CG"
                                                    then "S"
                                                    else
                                                        if stringVal == "AT"
                                                            then "W"
                                                            else
                                                                if stringVal == "GT"
                                                                    then "K"
                                                                    else
                                                                        if stringVal == "AC"
                                                                            then "M"
                                                                            else
                                                                                if stringVal == "CGT"
                                                                                    then "B"
                                                                                    else
                                                                                        if stringVal == "AGT"
                                                                                            then "D"
                                                                                            else
                                                                                                if stringVal == "ACT"
                                                                                                    then "H"
                                                                                                    else
                                                                                                        if stringVal == "ACG"
                                                                                                            then "V"
                                                                                                            else
                                                                                                                if stringVal == "ACGT"
                                                                                                                    then "N"
                                                                                                                    else
                                                                                                                        if stringVal == "-ACGT"
                                                                                                                            then "?"
                                                                                                                            else -- ours for gap chars and nuc

                                                                                                                                if stringVal == "-A"
                                                                                                                                    then "a"
                                                                                                                                    else
                                                                                                                                        if stringVal == "-C"
                                                                                                                                            then "c"
                                                                                                                                            else
                                                                                                                                                if stringVal == "-G"
                                                                                                                                                    then "g"
                                                                                                                                                    else
                                                                                                                                                        if stringVal == "-T"
                                                                                                                                                            then "t"
                                                                                                                                                            else
                                                                                                                                                                if stringVal == "-AG"
                                                                                                                                                                    then "r"
                                                                                                                                                                    else
                                                                                                                                                                        if stringVal == "-CT"
                                                                                                                                                                            then "y"
                                                                                                                                                                            else
                                                                                                                                                                                if stringVal == "-CG"
                                                                                                                                                                                    then "s"
                                                                                                                                                                                    else
                                                                                                                                                                                        if stringVal == "-AT"
                                                                                                                                                                                            then "w"
                                                                                                                                                                                            else
                                                                                                                                                                                                if stringVal == "-GT"
                                                                                                                                                                                                    then "k"
                                                                                                                                                                                                    else
                                                                                                                                                                                                        if stringVal == "-AC"
                                                                                                                                                                                                            then "m"
                                                                                                                                                                                                            else
                                                                                                                                                                                                                if stringVal == "-CGT"
                                                                                                                                                                                                                    then "b"
                                                                                                                                                                                                                    else
                                                                                                                                                                                                                        if stringVal == "-AGT"
                                                                                                                                                                                                                            then "d"
                                                                                                                                                                                                                            else
                                                                                                                                                                                                                                if stringVal == "-ACT"
                                                                                                                                                                                                                                    then "h"
                                                                                                                                                                                                                                    else
                                                                                                                                                                                                                                        if stringVal == "-ACG"
                                                                                                                                                                                                                                            then "v"
                                                                                                                                                                                                                                            else "Unrecognized nucleic acid ambiguity code : " <> "|" <> stringVal <> "|"
            else -- AA IUPAC

                if isAlphabetAminoAcid localAlphabet && (SET.size (alphabetSymbols localAlphabet) > 5)
                    then
                        if length stringVal == 1
                            then stringVal <> " "
                            else
                                if stringVal == "DN"
                                    then "B"
                                    else
                                        if stringVal == "EQ"
                                            then "Z"
                                            else
                                                if stringVal == "ACDEFGHIKLMNPQRSTVWY"
                                                    then "X"
                                                    else
                                                        if stringVal == "-ACDEFGHIKLMNPQRSTVWY"
                                                            then "?"
                                                            else -- amino acid polymorphisms without ambiguity codes
                                                                "[" <> stringVal <> "]" <> " "
                                                                
                    else (bitVectToCharState'' localAlphabetNEString localAlphabetVect bitValue) <> " "
    where
        indices = [0 .. len - 1]
        len = length vec
        vec = localAlphabetVect -- alphabetSymbols localAlphabet
        pollSymbol i polled
            | bitValue `testBit` i = (vec V.! i) : polled
            | otherwise = polled


-- bitVectToCharState''  takes a bit vector representation and returns a list states as integers
bitVectToCharState'' ∷ (Bits b) ⇒ NonEmpty String → V.Vector String → b → String
bitVectToCharState'' localAlphabet localAlphabetVect bitValue
    | isAlphabetDna hereAlphabet = fold $ iupacToDna BM.!> (NE.fromList observedSymbols)
    | isAlphabetAminoAcid hereAlphabet = fold $ iupacToAminoAcid BM.!> (NE.fromList observedSymbols)
    | otherwise = -- L.intercalate "," $ toList observedSymbols
        let symbolList = toList observedSymbols 
            symbolComma = L.intercalate "," observedSymbols
        in
        if length observedSymbols == 1 then symbolComma
        else if length observedSymbols == 0 then "-"
        else ('[' : symbolComma) <> "]"
    where
        hereAlphabet = fromSymbols localAlphabet
        symbolCountH = length localAlphabet
        observedSymbols =
            -- NE.fromList $
                foldMap
                    -- (\ i -> [localAlphabet NE.!! i | bitValue `testBit` i])
                    (\i → [localAlphabetVect V.! i | bitValue `testBit` i])
                    [0 .. symbolCountH - 1]


-- | matrixStateToStringtakes a matrix state and returns a string representation
matrixStateToString ∷ V.Vector MatrixTriple → String
matrixStateToString inStateVect =
    let minCost = V.minimum $ fmap fst3 inStateVect
        minCostStates = V.toList $ V.filter ((== minCost) . fst3) inStateVect
        statesStringList = fmap show minCostStates
    in  if length statesStringList == 1
            then head statesStringList
            else "[" <> unwords statesStringList <> "]"


{- | additivStateToString take an additive range and prints single state if range equal or
[ab] if not
[a-b] causes problems with TNT
-}
additivStateToString ∷ V.Vector String → (Int, Int) → String
additivStateToString localAlphabet (a, b) =
    if a == b
        then show a
        else
            if (show a == V.head localAlphabet) && (show b == V.last localAlphabet)
                then "?"
                else "[" <> show a <> show b <> "]"


{- | filledDataFields takes rawData and checks taxon to see what percent
"characters" are found.
call with (0,0)
-}
filledDataFields ∷ (Int, Int) → TermData → (NameText, Int, Int)
filledDataFields (hasData, totalData) (taxName, taxData)
    | null taxData = (taxName, hasData, totalData)
    | ST.length (head taxData) == 0 = filledDataFields (hasData, 1 + totalData) (taxName, tail taxData)
    | otherwise = filledDataFields (1 + hasData, 1 + totalData) (taxName, tail taxData)


{- | stripComments removes all lines that being with "--" haskell stype comments
needs to be reversed on return to maintain order.
-}
stripComments ∷ [String] → [String]
stripComments inStringList =
    if null inStringList
        then []
        else
            let strippedLine = GU.stripString $ head inStringList
            in  if null strippedLine
                    then stripComments $ tail inStringList
                    else
                        if length strippedLine < 2
                            then strippedLine : stripComments (tail inStringList)
                            else
                                if "--" == take 2 strippedLine
                                    then stripComments $ tail inStringList
                                    else strippedLine : stripComments (tail inStringList)


-- | getDecoratedGraphBlockCharInformation takes decorated graph and reports number of blosk and size of each
getDecoratedGraphBlockCharInformation ∷ DecoratedGraph → ((Int, Int), [V.Vector Int])
getDecoratedGraphBlockCharInformation inGraph =
    if LG.isEmpty inGraph
        then ((0, 0), [])
        else -- get a vertices from graph and take their information

            let inVertDataList = fmap (vertData . snd) (LG.labNodes inGraph)
                blockNumMax = maximum $ fmap length inVertDataList
                blocknumMin = minimum $ fmap length inVertDataList
                blockLengthList = fmap (fmap length) inVertDataList
            in  ((blockNumMax, blocknumMin), blockLengthList)


{- | vectMaybeHead takes a vector and returns JUst V.head if not V.empty
Nothing otherwise
-}
vectMaybeHead ∷ V.Vector a → Maybe a
vectMaybeHead inVect =
    if V.null inVect
        then Nothing
        else Just (V.head inVect)


-- vectResolveMaybe takes a Vector of Maybe a
-- and returns Just a or V.empty
vectResolveMaybe ∷ V.Vector (Maybe a) → V.Vector a
vectResolveMaybe inVect =
    -- trace ("VRM " <> show (length inVect)) $
    if isNothing (V.head inVect)
        then V.empty
        else V.singleton $ fromJust $ V.head inVect


{- | getNumberPrealignedCharacters takes processed data and returns the number of prealigned sequence characters
used to special case procedurs with prealigned sequences
-}
getNumberPrealignedCharacters ∷ V.Vector BlockData → Int
getNumberPrealignedCharacters blockDataVect =
    if V.null blockDataVect
        then 0
        else
            let firstBlock = GU.thd3 $ V.head blockDataVect
                characterTypes = V.map charType firstBlock
                sequenceChars = length $ V.filter id $ V.map (`elem` prealignedCharacterTypes) characterTypes
            in  sequenceChars + getNumberPrealignedCharacters (V.tail blockDataVect)


{- | getNumberNonExactCharacters takes processed data and returns the number of non-exact (unaligned seqeunce) characters
used to special case procedures with unaligned sequences
-}
getNumberNonExactCharacters ∷ V.Vector BlockData → Int
getNumberNonExactCharacters blockDataVect =
    if V.null blockDataVect
        then 0
        else
            let firstBlock = GU.thd3 $ V.head blockDataVect
                characterTypes = V.map charType firstBlock
                sequenceChars = length $ V.filter id $ V.map (`elem` nonExactCharacterTypes) characterTypes
            in  sequenceChars + getNumberNonExactCharacters (V.tail blockDataVect)


{- | getNumberSequenceCharacters takes processed data and returns the number of non-exact (= sequence) characters
utilized to special case datasets with limited non-exact characters
-}
getNumberSequenceCharacters ∷ V.Vector BlockData → Int
getNumberSequenceCharacters blockDataVect =
    if V.null blockDataVect
        then 0
        else
            let firstBlock = GU.thd3 $ V.head blockDataVect
                characterTypes = V.map charType firstBlock
                sequenceChars = length $ V.filter id $ V.map (`elem` sequenceCharacterTypes) characterTypes
            in  sequenceChars + getNumberSequenceCharacters (V.tail blockDataVect)


{- | getNumber4864PackedChars takes processed data and returns the number of
packed characters with states  in 4, 8, and 64
this for NCM since weightiong may be apperoximate and needs to be rediagnosed
-}
getNumber4864PackedChars ∷ V.Vector BlockData → Int
getNumber4864PackedChars blockDataVect =
    if V.null blockDataVect
        then 0
        else
            let firstBlock = GU.thd3 $ V.head blockDataVect
                characterTypes = V.map charType firstBlock
                packedChars = length $ V.filter id $ V.map (`elem` [Packed4, Packed8, Packed64]) characterTypes
            in  packedChars + getNumber4864PackedChars (V.tail blockDataVect)


{- | has4864PackedChars takes processed data and if has
packed characters with states  in 4, 8, and 64
this for NCM since weightiong may be apperoximate and needs to be rediagnosed
-}
has4864PackedChars ∷ V.Vector BlockData → Bool
has4864PackedChars blockDataVect =
    not (V.null blockDataVect)
        && ( let firstBlock = GU.thd3 $ V.head blockDataVect
                 characterTypes = V.map charType firstBlock
                 packedChars = length $ V.filter id $ V.map (`elem` [Packed4, Packed8, Packed64]) characterTypes
             in  ((packedChars > 0) || has4864PackedChars (V.tail blockDataVect))
           )


{- | getLengthSequenceCharacters takes processed data and returns the total length (maximum) of non-exact (= sequence) characters
utilised to get rough estimate of fraction of non-exact characters for
dynamic epsilon adjustment to data types
maximum since that ius th eminimum length of optimized HTU sequences
-}
getLengthSequenceCharacters ∷ V.Vector BlockData → Int
getLengthSequenceCharacters blockDataVect =
    if V.null blockDataVect
        then 0
        else
            let -- get character info
                firstBlock = GU.thd3 $ V.head blockDataVect
                characterTypes = V.map charType firstBlock

                -- get sequences in block
                firstBlockCharacters = GU.snd3 $ V.head blockDataVect
                (sequenceCharVect, _) = V.unzip $ V.filter snd (V.zip firstBlockCharacters (V.map (`elem` sequenceCharacterTypes) characterTypes))

                -- get max length sequence data
                sequenceCharsLength = V.sum $ fmap (V.maximum . fmap getMaxCharLength) sequenceCharVect
            in  -- trace ("GLSC: " <> (show (sequenceCharsLength, V.length firstBlockCharacters, fmap V.length firstBlockCharacters)))
                sequenceCharsLength + getLengthSequenceCharacters (V.tail blockDataVect)


-- | getMaxCharLength takes characterData and returns the length of the longest character field from preliminary
getMaxCharLength ∷ CharacterData → Int
getMaxCharLength inChardata =
    let nonAdd = (V.length . snd3) $ stateBVPrelim inChardata
        add = (V.length . snd3) $ rangePrelim inChardata
        matrix = V.length $ matrixStatesPrelim inChardata
        slim = (SV.length . snd3) $ slimGapped inChardata
        wide = (UV.length . snd3) $ wideGapped inChardata
        huge = (V.length . snd3) $ hugeGapped inChardata
        aSlim = (SV.length . snd3) $ alignedSlimPrelim inChardata
        aWide = (UV.length . snd3) $ alignedWidePrelim inChardata
        aHuge = (V.length . snd3) $ alignedHugePrelim inChardata
        packed = 32 * ((UV.length . snd3) $ packedNonAddPrelim inChardata)
        aBitChar = (UV.head . snd3) $ packedNonAddPrelim inChardata
        missingBitChar = complement $ aBitChar `xor` aBitChar
        packedNotMissingChars = 32 * (UV.length $ UV.filter (/= missingBitChar) $ snd3 $ packedNonAddPrelim inChardata)
    in  -- trace ("GMCL: " <> (show [nonAdd, add, matrix, slim, wide, huge, aSlim, aWide, aHuge, packed]))
        -- trace ("GMCL: " <> show (nonAdd, add, matrix, packed, packedNotMissingChars)) $
        maximum [nonAdd, add, matrix, slim, wide, huge, aSlim, aWide, aHuge, packedNotMissingChars]


{- | getNumberExactCharacters takes processed data and returns the number of non-exact characters
ised to special case datasets with limited non-exact characters
-}
getNumberExactCharacters ∷ V.Vector BlockData → Int
getNumberExactCharacters blockDataVect =
    if V.null blockDataVect
        then 0
        else
            let firstBlock = GU.thd3 $ V.head blockDataVect
                characterTypes = V.map charType firstBlock
                exactChars = length $ V.filter id $ V.map (`elem` exactCharacterTypes) characterTypes
            in  exactChars + getNumberExactCharacters (V.tail blockDataVect)


{- | getPairwiseObservationsGraph gets the observations between a pairs of vertices on a graph
that are non-missing
used to normalize Wagner deltas
-}
getPairwiseObservationsGraph ∷ VertexBlockData → VertexBlockData → VertexCost
getPairwiseObservationsGraph vertexBlockDataI vertexBlockDataJ =
    let blockObsV = V.zipWith getPairedBlockObs vertexBlockDataI vertexBlockDataJ
    in  if V.null vertexBlockDataI || V.null vertexBlockDataJ
            then 0.0
            else fromIntegral $ V.sum blockObsV


-- | getPairedBlockObs takes zipped block of chardata vectors and returns block paired observaiton that are non-missing
getPairedBlockObs ∷ V.Vector CharacterData → V.Vector CharacterData → Int
getPairedBlockObs charVectI charVectJ =
    if V.null charVectI || V.null charVectJ
        then 0
        else V.sum $ V.zipWith getPairCharNonMissing charVectI charVectJ


-- | getPairCharNonMissing gets non-missing character numbers from a pair of characters
getPairCharNonMissing ∷ CharacterData → CharacterData → Int
getPairCharNonMissing iTaxon jTaxon =
    let obsI = getMaxCharLength iTaxon
        obsJ = getMaxCharLength jTaxon
    in  -- trace ("GPCL: " <> (show (iIndex, jIndex, V.length charDataV))) $
        if obsI == 0 || obsJ == 0
            then 0
            else max obsI obsJ


{- | getPairwiseObservations gets the observations between a pairs of leaves that are non-missing
used to normalize distances
-}
getPairwiseObservations ∷ V.Vector BlockData → (Int, Int) → VertexCost
getPairwiseObservations blocKDataV pairTax =
    if V.null blocKDataV
        then 0
        else fromIntegral $ V.sum (fmap (getPairBlockObs pairTax) blocKDataV)


-- | getMaxBlockObs gets the supremum over taxa number of characters in a block of data
getPairBlockObs ∷ (Int, Int) → BlockData → Int
getPairBlockObs pairTax (_, charDataVV, _) =
    if V.null charDataVV
        then 0
        else
            let newListList = L.transpose $ V.toList $ fmap V.toList charDataVV
                charTaxVect = V.fromList $ fmap V.fromList newListList
            in  -- trace ("GPBO: " <> (show (V.length charDataVV, V.length charTaxVect, fmap V.length charTaxVect))) $
                V.sum (fmap (getPairCharLength pairTax) charTaxVect)


{- | getPairBlockObs get non-missing observations between taxa
NB--does not go into qualitative or packed charcters and check for missing values
other than "all missing"  packed
-}
getPairCharLength ∷ (Int, Int) → V.Vector CharacterData → Int
getPairCharLength (iIndex, jIndex) charDataV =
    if V.null charDataV
        then 0
        else
            let iTaxon = charDataV V.! iIndex
                jTaxon = charDataV V.! jIndex
                obsI = getMaxCharLength iTaxon
                obsJ = getMaxCharLength jTaxon
            in  -- trace ("GPCL: " <> (show (iIndex, jIndex, V.length charDataV))) $
                if obsI == 0 || obsJ == 0
                    then 0
                    else max obsI obsJ


{- | getMaxNumberObservations takes data set and returns the supremum of character numbers from all
taxa over all charcaters (sequence and qiualitative)
used for various normalizations
-}
getMaxNumberObservations ∷ V.Vector BlockData → PhyG VertexCost
getMaxNumberObservations blocKDataV
    | V.null blocKDataV = pure 0
    | otherwise =
        getParallelChunkTraverse >>= \pTraverse →
            fmap (fromIntegral . sum) . pTraverse getMaxBlockObs $ V.toList blocKDataV


-- | getMaxBlockObs gets the supremum over taxa number of characters in a block of data
getMaxBlockObs ∷ BlockData → PhyG Int
getMaxBlockObs (_, charDataVV, _)
    | V.null charDataVV = pure 0
    | otherwise =
        let newListList = L.transpose $ V.toList $ fmap V.toList charDataVV
        in  getParallelChunkTraverse >>= \pTraverse →
                sum <$> pTraverse getSupCharLength newListList


{- | getMaxCharLength takes a vector of charcters and returns the supremum of observations for that character
over all taxa
-}
getSupCharLength ∷ [CharacterData] → PhyG Int
getSupCharLength charDataV
    | null charDataV = pure 0
    | otherwise =
        getParallelChunkMap <&> \pMap →
            maximum $ getMaxCharLength `pMap` charDataV


-- getFractionDynamic returns fraction (really of length) of dynamic charcters for adjustment to dynamicEpsilon
getFractionDynamic ∷ ProcessedData → Double
getFractionDynamic inData =
    let numStaticCharacters = getNumberExactCharacters $ thd3 inData
        lengthDynamicCharacters = getLengthSequenceCharacters $ thd3 inData
    in  fromIntegral lengthDynamicCharacters / fromIntegral (lengthDynamicCharacters + numStaticCharacters)


{- | splitBlockCharacters takes a block of characters (vector) and splits into two partitions of exact (Add, NonAdd, Matrix) and sequence characters
(= nonExact) using accumulators
-}
splitBlockCharacters
    ∷ V.Vector (V.Vector CharacterData)
    → V.Vector CharInfo
    → Int
    → [([CharacterData], CharInfo)]
    → [([CharacterData], CharInfo)]
    → (BlockData, BlockData)
splitBlockCharacters inDataVV inCharInfoV localIndex exactCharPairList seqCharPairList =
    if localIndex == V.length inCharInfoV
        then
            let (exactDataList, exactCharInfoList) = unzip exactCharPairList
                (sequenceDataList, sequenceCharInfoList) = unzip seqCharPairList
                newExactCharInfoVect = V.fromList $ reverse exactCharInfoList
                newSeqCharCharInfoVect = V.fromList $ reverse sequenceCharInfoList
                newExactData = V.fromList $ fmap (V.fromList . reverse) (L.transpose exactDataList)
                newSeqCharData = V.fromList $ fmap (V.fromList . reverse) (L.transpose sequenceDataList)
            in  ( (T.pack "ExactCharacters", newExactData, newExactCharInfoVect)
                , (T.pack "Non-ExactCharacters", newSeqCharData, newSeqCharCharInfoVect)
                )
        else
            let localCharacterType = charType (inCharInfoV V.! localIndex)
                thisCharacterData = V.toList $ fmap (V.! localIndex) inDataVV
                newPair = (thisCharacterData, inCharInfoV V.! localIndex)
            in  if localCharacterType `elem` exactCharacterTypes
                    then splitBlockCharacters inDataVV inCharInfoV (localIndex + 1) (newPair : exactCharPairList) seqCharPairList
                    else
                        if localCharacterType `elem` sequenceCharacterTypes
                            then splitBlockCharacters inDataVV inCharInfoV (localIndex + 1) exactCharPairList (newPair : seqCharPairList)
                            else error ("Unrecongized/implemented character type: " <> show localCharacterType)


-- | safeVectorHead safe vector head--throws error if null
safeVectorHead ∷ V.Vector a → a
safeVectorHead inVect =
    if V.null inVect
        then error "Empty vector in safeVectorHead"
        else V.head inVect


{- | get leftRightChildLabelBV takes a pair of vertex labels and returns left and right
based on their bitvector representation.  This ensures left/right consistancey in
pre and postoder passes, and with bitvectors of leaves determined by data hash,
ensures label invariance with repect to leaves
larger bitvector goes second (bigge)
-}
leftRightChildLabelBV ∷ (VertexInfo, VertexInfo) → (VertexInfo, VertexInfo)
leftRightChildLabelBV inPair@(firstNode, secondNode) =
    let firstLabel = bvLabel firstNode
        secondLabel = bvLabel secondNode
    in  if firstLabel > secondLabel
            then (secondNode, firstNode)
            else inPair


{- | get leftRightChildLabelBVNode takes a pair ofnodes and returns left and right
based on their bitvector representation.  This ensures left/right consistancey in
pre and postoder passes, and with bitvectors of leaves determined by data hash,
ensures label invariance with repect to leaves
larger bitvector goes second (bigge)
-}
leftRightChildLabelBVNode ∷ (LG.LNode VertexInfo, LG.LNode VertexInfo) → (LG.LNode VertexInfo, LG.LNode VertexInfo)
leftRightChildLabelBVNode inPair@(firstNode, secondNode) =
    let firstLabel = bvLabel $ snd firstNode
        secondLabel = bvLabel $ snd secondNode
    in  if firstLabel > secondLabel
            then (secondNode, firstNode)
            else inPair


{- | prettyPrintVertexInfo returns a string with formated version of
vertex info
-}
prettyPrintVertexInfo ∷ VertexInfo → String
prettyPrintVertexInfo inVertData =
    let zerothPart = "Vertex name " <> T.unpack (vertName inVertData) <> " Index " <> show (index inVertData)
        firstPart = "\n\tBitVector (as number) " <> show ((BV.toUnsignedNumber $ bvLabel inVertData) ∷ Int)
        secondPart = "\n\tParents " <> show (parents inVertData) <> " Children " <> show (children inVertData)
        thirdPart =
            "\n\tType "
                <> show (nodeType inVertData)
                <> " Local Cost "
                <> show (vertexCost inVertData)
                <> " SubGraph Cost "
                <> show (subGraphCost inVertData)
        fourthPart =
            "\n\tData Blocks: "
                <> show (V.length $ vertData inVertData)
                <> " Characters (by block) "
                <> show (V.length <$> vertData inVertData)
        fifthPart = "\n\t" <> show (vertData inVertData)
    in  zerothPart <> firstPart <> secondPart <> thirdPart <> fourthPart <> fifthPart


-- | add3 adds three values
add3 ∷ (Num a) ⇒ a → a → a → a
add3 x y z = x + y + z


{- | getProcessDataByBlock takes ProcessData and returns a list of Processed data with one block
per processed data element
argument to filter terminals with missing taxa
wraps around getProcessDataByBlock' with counter
-}
getProcessDataByBlock ∷ Bool → ProcessedData → [ProcessedData]
getProcessDataByBlock filterMissing (nameVect, nameBVVect, blockDataVect) = reverse $ getProcessDataByBlock' filterMissing 0 (nameVect, nameBVVect, blockDataVect)


{- | getProcessDataByBlock' called by getProcessDataByBlock with counter
and later reversed
-}
getProcessDataByBlock' ∷ Bool → Int → ProcessedData → [ProcessedData]
getProcessDataByBlock' filterMissing counter (nameVect, nameBVVect, blockDataVect)
    | V.null blockDataVect = []
    | counter == V.length blockDataVect = []
    | otherwise =
        let thisBlockData = blockDataVect V.! counter
        in  if not filterMissing
                then
                    (nameVect, nameBVVect, V.singleton thisBlockData)
                        : getProcessDataByBlock' filterMissing (counter + 1) (nameVect, nameBVVect, blockDataVect)
                else
                    let (blockName, charDataLeafVect, blockCharInfo) = thisBlockData
                        isMissingVect = V.map V.null charDataLeafVect
                        (nonMissingNameVect, nonMissingBVVect, nonMissingLeafData, _) = V.unzip4 $ V.filter (not . GU.fth4) (V.zip4 nameVect nameBVVect charDataLeafVect isMissingVect)
                        nonMissingBlockData = (blockName, nonMissingLeafData, blockCharInfo)
                    in  (nonMissingNameVect, nonMissingBVVect, V.singleton nonMissingBlockData)
                            : getProcessDataByBlock' filterMissing (counter + 1) (nameVect, nameBVVect, blockDataVect)


{- | copyToNothing takes VertexBlockData and copies to VertexBlockDataMaybe
data as nothing
-}
copyToNothing ∷ VertexBlockData → VertexBlockDataMaybe
copyToNothing = fmap setNothing
    where
        setNothing a = V.replicate (V.length a) Nothing


{- | copyToJust takes VertexBlockData and copies to VertexBlockDataMaybe
data as Just CharacterData
-}
copyToJust ∷ VertexBlockData → VertexBlockDataMaybe
copyToJust = fmap (fmap Just)


{- | simAnnealAccept takes simulated annealing parameters, current best graph (e) cost,
candidate graph cost (e') and a uniform random integer and returns a Bool to accept or reject
the candidate solution
the basic method is
 1) accepts if current is better
 2) Otherwise prob accept = exp(-(e' -e)/T)
where T is a step from max to min
maxT and minT can probbaly be set to 100 and 1 or something but leaving some flexibility
curStep == 0 random walk (always accept)
curStep == (numSteps -1) greedy False is not better
-}
simAnnealAccept ∷ Maybe SAParams → VertexCost → VertexCost → PhyG (Bool, Maybe SAParams)
simAnnealAccept inParams curBestCost candCost = case inParams of
    Nothing → error "simAnnealAccept Simulated anneling parameters = Nothing"
    Just simAnealVals → case method simAnealVals of
        -- drifting probs
        Drift → driftAccept inParams curBestCost candCost
        _ →
            -- simulated annealing probs
            let numSteps = numberSteps simAnealVals
                curStep = currentStep simAnealVals

                -- stepFactor =  (fromIntegral $ numSteps - curStep) / (fromIntegral numSteps)
                -- tempFactor = curBestCost  * stepFactor

                candCost'
                    | curBestCost == candCost = candCost + 1
                    | otherwise = candCost

                -- factors here for tweaking
                energyFactor = 10.0 * (100 * (curBestCost - candCost') / curBestCost)
                tempFactor' = 10.0 * fromIntegral (numSteps - curStep) / fromIntegral numSteps

                -- flipped order - (e' -e)
                -- probAcceptance = exp ((curBestCost - candCost) / ((maxTemp - minTemp) * tempFactor))
                -- probAcceptance' = exp ( (fromIntegral (curStep + 1)) * (curBestCost - candCost') / tempFactor)

                probAcceptance = exp (energyFactor / tempFactor')

                -- multiplier for resolution 1000, 100 prob be ok
                randMultiplier ∷ Word
                randMultiplier = 1000
                intAccept = floor $ fromIntegral randMultiplier * probAcceptance

                nextSAParams = simAnealVals{currentStep = curStep + 1}

                withUpdatedParams ∷ Bool → (Bool, Maybe SAParams)
                withUpdatedParams b = (b, Just nextSAParams)

                costCheck
                    -- lowest cost-- greedy
                    -- but increment this if using heuristic costs
                    | candCost < curBestCost = pure True
                    -- not better and at lowest temp
                    | curStep >= numSteps - 1 = pure False
                    -- test for non-lowest temp conditions
                    | otherwise =
                        -- use remainder for testing
                        (< intAccept) . snd . (`divMod` randMultiplier) . abs <$> getRandom
            in  withUpdatedParams <$> costCheck


-- | incrementSimAnnealParams increments the step number by 1 but returns all other the same
incrementSimAnnealParams ∷ Maybe SAParams → Maybe SAParams
incrementSimAnnealParams =
    let incrementor params = case method params of
            SimAnneal → params{currentStep = currentStep params + 1}
            _ → params{driftChanges = driftChanges params + 1}
    in  fmap incrementor


-- | generateRandLists generates n random lists from seed
generateRandIntLists ∷ Int → PhyG [[Int]]
generateRandIntLists count = replicateM count getRandoms


{- | generateUniqueRandList take a int and simulated anealing parameter slist and creates
a list of SA paramter values with unique rnandomInt lists
sets current step to 0
-}
generateUniqueRandList ∷ Int → Maybe SAParams → [Maybe SAParams]
generateUniqueRandList number inParams = replicate number inParams


{- | driftAccept takes SAParams, currrent best cost, and candidate cost
and returns a Boolean and an incremented set of params
this based on a percentage of diffference in graph cost
-}
driftAccept ∷ Maybe SAParams → VertexCost → VertexCost → PhyG (Bool, Maybe SAParams)
driftAccept simAnealVals curBestCost candCost = case simAnealVals of
    Nothing → error "Nothing value in driftAccept"
    Just params →
        let -- prob acceptance for better, same, and worse costs
            probAcceptance
                | candCost < curBestCost = 1.0
                | candCost == curBestCost = driftAcceptEqual params
                | otherwise = 1.0 / (driftAcceptWorse params + (100.0 * (candCost - curBestCost) / curBestCost))

            -- multiplier for resolution 1000, 100 prob be ok
            randMultiplier ∷ Word
            randMultiplier = 1000
            intAccept = floor $ fromIntegral randMultiplier * probAcceptance

            -- not always incrementing becasue may not result in changes
            nextSAParams = Just $ params{driftChanges = driftChanges params + 1}
            nextSAPAramsNoChange = simAnealVals

            resultParams
                -- only increment numberof changes for True values
                -- but increment this if using heuristic costs
                | candCost < curBestCost = pure (True, nextSAParams)
                | otherwise = do
                    -- use remainder for testing--passing infinite list and take head
                    intRandVal ← snd . (`divMod` randMultiplier) . abs <$> getRandom
                    pure $
                        if intRandVal < intAccept
                            then -- trace ("Drift T: " <> (show (curNumChanges, candCost, curBestCost, probAcceptance, intAccept, intRandVal)) <> " True")
                                (True, nextSAParams)
                            else -- trace ("Drift F: " <> (show (curNumChanges, candCost, curBestCost, probAcceptance, intAccept, intRandVal)) <> " False")
                                (False, nextSAPAramsNoChange)
        in  resultParams


-- | getTraversalCosts takes a Phylogenetic Graph and returns costs of traversal trees
getTraversalCosts ∷ PhylogeneticGraph → [VertexCost]
getTraversalCosts inGraph =
    let traversalTrees = V.toList (V.toList <$> fft6 inGraph)
        traversalRoots = fmap (head . LG.getRoots) (concat traversalTrees)
        traversalRootCosts = fmap (subGraphCost . snd) traversalRoots
    in  traversalRootCosts


{- SAme as getCharacterLength
-- | getSequenceCharacterLengths returns a the length of block characters
getSequenceCharacterLengths :: CharacterData -> CharInfo -> Int
getSequenceCharacterLengths inCharData inCharInfo =
    let inCharType = charType inCharInfo
    in
    -- trace ("GCL:" <> (show inCharType) <> " " <> (show $ snd3 $ stateBVPrelim inCharData)) (
    case inCharType of
      x | x == NonAdd -> 0 -- V.length  $ snd3 $ stateBVPrelim inCharData
      x | x `elem` packedNonAddTypes   -> 0 -- UV.length  $ snd3 $ packedNonAddPrelim inCharData
      x | x == Add -> 0 -- V.length  $ snd3 $ rangePrelim inCharData
      x | x == Matrix -> 0 -- V.length  $ matrixStatesPrelim inCharData
      x | x `elem` [SlimSeq, NucSeq  ] -> SV.length $ snd3 $ slimAlignment inCharData
      x | x `elem` [WideSeq, AminoSeq] -> UV.length $ snd3 $ wideAlignment inCharData
      x | x == HugeSeq           -> V.length  $ snd3 $ hugeAlignment inCharData
      x | x == AlignedSlim       -> SV.length $ snd3 $ alignedSlimPrelim inCharData
      x | x == AlignedWide       -> UV.length $ snd3 $ alignedWidePrelim inCharData
      x | x == AlignedHuge       -> V.length  $ snd3 $ alignedHugePrelim inCharData
      _                                -> error ("Un-implemented data type " <> show inCharType)
      -- )
-}

-- | getCharacterLengths returns a the length of block characters
getCharacterLength ∷ Bool → CharacterData → CharInfo → Int
getCharacterLength useIA inCharData inCharInfo =
    let inCharType = charType inCharInfo
    in  -- trace ("GCL:" <> (show inCharType) <> " " <> (show $ snd3 $ stateBVPrelim inCharData)) (
        case inCharType of
            x | x == NonAdd → V.length $ snd3 $ stateBVPrelim inCharData
            x | x `elem` packedNonAddTypes → UV.length $ snd3 $ packedNonAddPrelim inCharData
            x | x == Add → V.length $ snd3 $ rangePrelim inCharData
            x | x == Matrix → V.length $ matrixStatesPrelim inCharData
            x | x `elem` [SlimSeq, NucSeq] && useIA → SV.length $ snd3 $ slimAlignment inCharData -- slimAlignment inCharData
            x | x `elem` [SlimSeq, NucSeq] → SV.length $ snd3 $ slimGapped inCharData -- slimAlignment inCharData
            x | x `elem` [WideSeq, AminoSeq] && useIA → UV.length $ snd3 $ wideAlignment inCharData --  wideAlignment inCharData
            x | x `elem` [WideSeq, AminoSeq] → UV.length $ snd3 $ wideGapped inCharData --  wideAlignment inCharData
            x | x == HugeSeq && useIA → V.length $ snd3 $ hugeAlignment inCharData --  hugeAlignment inCharData
            x | x == HugeSeq → V.length $ snd3 $ hugeGapped inCharData --  hugeAlignment inCharData
            x | x == AlignedSlim → SV.length $ snd3 $ alignedSlimPrelim inCharData
            x | x == AlignedWide → UV.length $ snd3 $ alignedWidePrelim inCharData
            x | x == AlignedHuge → V.length $ snd3 $ alignedHugePrelim inCharData
            _ → error ("Un-implemented data type " <> show inCharType)


-- )

-- | getCharacterLengths' flipped arg version of getCharacterLength
getCharacterLength' ∷ Bool → CharInfo → CharacterData → Int
getCharacterLength' useIA inCharInfo inCharData = getCharacterLength useIA inCharData inCharInfo


-- | getMaxCharacterLengths get maximum charcter legnth from a list
getMaxCharacterLength ∷ Bool → CharInfo → [CharacterData] → Int
getMaxCharacterLength useIA inCharInfo inCharDataList = maximum $ fmap (getCharacterLength' useIA inCharInfo) inCharDataList


-- | getSingleTaxon takes a taxa x characters block and an index and returns the character vector for that index
getSingleTaxon ∷ V.Vector (V.Vector CharacterData) → Int → V.Vector CharacterData
getSingleTaxon singleCharVect taxonIndex = fmap (V.! taxonIndex) singleCharVect


{- | glueBackTaxChar takes single chartacter taxon vectors and glues them back inot multiple characters for each
taxon as expected in Blockdata.  Like a transpose.  FIlters out zero length characters
-}
glueBackTaxChar ∷ V.Vector (V.Vector CharacterData) → V.Vector (V.Vector CharacterData)
glueBackTaxChar singleCharVect =
    let numTaxa = V.length $ V.head singleCharVect
        multiCharVect = fmap (getSingleTaxon singleCharVect) (V.fromList [0 .. numTaxa - 1])
    in  multiCharVect

{- | concatFastas is used by "report(ia)" to create a single concatenated fasta
for use by programs such as RAxML andf TNT
takes a single string of multiple fasta output, each entity to be concated (by taxon name)
separate by a line starting with "Seqeunce character" from ia output
there must be the same number of taxa and names in each element to be concatenated
-}
concatFastas ∷ String → String
concatFastas inMultFastaString =
    if null inMultFastaString
        then []
        else -- split on "Sequence character"

            let fastaFileStringList = filter (not . null) $ spitIntoFastas inMultFastaString

                -- make pairs from each "file"
                fastaTaxDataPairLL = fmap fasta2PairList fastaFileStringList

                fastTaxNameList = (fst <$> head fastaTaxDataPairLL)

                -- merge sequences with (<>)
                fastaDataLL = fmap (fmap snd) fastaTaxDataPairLL
                mergedData = mergeFastaData fastaDataLL

                -- create new tuples
                newFastaPairs = zipWith (<>) fastTaxNameList mergedData
            in  -- trace ("CF:" <> (show $ length fastaFileStringList) <> " " <> (show $ fmap length fastaFileStringList))
                concat newFastaPairs


{- | spitIntoFastas takes String generated by reprot ia functions,
splits on "Sequence character" line, creating separate fastas
-}
spitIntoFastas ∷ String → [String]
spitIntoFastas inString =
    if null inString
        then []
        else
            let linesList = splitFastaLines [] [] $ filter (not . null) $ lines inString
                fastaStringList = fmap unlines linesList
            in  fastaStringList


-- | splitFastaLines splits a list of lines into lists based on the line "Sequence character"
splitFastaLines ∷ [String] → [[String]] → [String] → [[String]]
splitFastaLines curFasta curFastaList inLineList =
    if null inLineList
        then reverse (reverse curFasta : curFastaList)
        else
            let firstLine = head inLineList
                firstWord = if (not $ null $ words firstLine) then head $ words firstLine
                            else ""
            in  if null firstLine
                    then splitFastaLines curFasta curFastaList (tail inLineList)
                    else
                        if firstWord == "Sequence"
                            then splitFastaLines [] (reverse curFasta : curFastaList) (tail inLineList)
                            else splitFastaLines (firstLine : curFasta) curFastaList (tail inLineList)


-- | mergeFastaData takes list of list of Strings and merges into a single list of Strings
mergeFastaData ∷ [[String]] → [String]
mergeFastaData inDataLL =
    if null $ head inDataLL
        then []
        else
            let firstData = concatMap head inDataLL
                formattedData = unlines $ SL.chunksOf 50 firstData
            in  formattedData : mergeFastaData (fmap tail inDataLL)


{- | fasta2PairList takes an individual fasta file as single string and returns
pairs data of taxon name and data
assumes single character alphabet
deletes '-' (unless "prealigned"), and spaces
-}
fasta2PairList ∷ String → [(String, String)]
fasta2PairList fileContents' =
    if null fileContents'
        then []
        else
            let fileContents = unlines $ filter (not . null) $ lines fileContents'
            in  if null fileContents
                    then []
                    else
                        if head fileContents /= '>'
                            then errorWithoutStackTrace "\n\n'Read' command error: fasta file must start with '>'"
                            else
                                let terminalSplits = SL.splitWhen (== '>') fileContents

                                    -- tail because initial split will an empty text
                                    pairData = getRawDataPairs (tail terminalSplits)
                                in  -- trace ("FSPL: " <> (show $ length pairData) <> " " <> (show $ fmap fst pairData))
                                    pairData


-- | getRawDataPairstakes splits of Text and returns terminalName, Data pairs--minimal error checking
getRawDataPairs ∷ [String] → [(String, String)]
getRawDataPairs inList =
    if null inList
        then []
        else
            let firstText = head inList
                firstName = head $ lines firstText
                firstData = concat $ tail $ lines firstText
            in  ('>' : firstName <> "\n", firstData) : getRawDataPairs (tail inList)


{- | hasResolutionDuplicateEdges checks resolution subtree in verteexInfo for duplicate d edges in softwired
subtree field
-}
hasResolutionDuplicateEdges ∷ ResolutionData → Bool
hasResolutionDuplicateEdges inResData =
    let edgeList = snd $ displaySubGraph inResData
        edgeDupList = length $ fmap LG.toEdge edgeList L.\\ L.nub (fmap LG.toEdge edgeList)
    in  edgeDupList > 0


-- | getUnionFieldsNode reutnrs String of union fields nicely
getUnionFieldsNode ∷ VertexBlockData → String
getUnionFieldsNode inVertData =
    "UnionFields " <> show inVertData


-- | transposeVector transposes a vrctor via conversino to lists transposing those and converting back
transposeVector ∷ V.Vector (V.Vector a) → V.Vector (V.Vector a)
transposeVector inVect =
    if V.null inVect
        then V.empty
        else
            let newListList = L.transpose $ V.toList $ fmap V.toList inVect
            in  V.fromList $ fmap V.fromList newListList


{- | combineMatrices takes two matrices [[a]] and applied function to be ziped and cretes new matrix
better be small becasue of list access would be n^3
-}
combineMatrices ∷ (a → a → a) → [[a]] → [[a]] → [[a]]
combineMatrices f m1 m2
    | null m1 = error "Null matrix 1 in combineMatrices"
    | null m2 = error "Null matrix 2 in combineMatrices"
    | (fmap length m1) /= (fmap length m2) =
        error ("Cannot combine matrices with unequal dimensions " <> (show (fmap length m1) <> " " <> show (fmap length m2)))
    | otherwise = zipWith (zipWith f) m1 m2


-- | normalizeMatrix takes a [[Int]] and returns normalized, symmetrical frequencies
normalizeMatrix ∷ [[Int]] → [[Double]]
normalizeMatrix inMatrix =
    if null inMatrix
        then []
        else
            let totalNumber = (sum $ (fmap sum) inMatrix) ∷ Int
                transposeMatrix = L.transpose inMatrix
                symMatrix = zipWith (zipWith (+)) inMatrix transposeMatrix
                newMatrix = fmap (fmap (/ (fromIntegral totalNumber))) $ fmap (fmap fromIntegral) symMatrix
            in  newMatrix


{- | divideList takes a [a] and returns [[a]]
    by dividing according to lengths in input list to lengths
-}
divideList ∷ [Int] → [a] → [[a]]
divideList lengthList inList =
    if null lengthList
        then []
        else
            if null inList
                then []
                else
                    if sum lengthList /= length inList
                        then error ("List division and list length do not match: " <> (show (sum lengthList)) <> " versus " <> (show $ length inList))
                        else
                            let firstLength = head lengthList
                            in  (take firstLength inList) : divideList (tail lengthList) (drop firstLength inList)
