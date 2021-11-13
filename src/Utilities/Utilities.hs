{- |
Module      :  Utilities.hs
Description :  Module specifying utility functions for use with PhyGraph
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

module Utilities.Utilities  where

import qualified Data.Vector as V
import           Data.Maybe
import qualified Data.List as L
import           Data.BitVector.LittleEndian (BitVector)
import qualified Data.BitVector.LittleEndian as BV
import Data.List.NonEmpty (NonEmpty(..))
import qualified Data.List.NonEmpty as NE
import Bio.Character.Encodable.Dynamic  -- this has DynamicCharacter reference
import Data.Alphabet
import Data.Alphabet.IUPAC
-- import Data.Alphabet.Special
import Data.Foldable
import Data.Bits
import Data.Set (Set)
import qualified Data.Set as Set
import qualified Data.Bimap as BM
import Types.Types
import           Data.Text.Short (ShortText)
import qualified Data.Text.Short             as ST
import qualified GeneralUtilities as GU
import Debug.Trace
import qualified Utilities.LocalGraph as LG
import qualified Data.Text.Lazy  as T
import qualified Data.Vector.Storable        as SV
import qualified Data.Vector.Unboxed         as UV
import qualified SymMatrix                   as S

-- | getblockInsertDataCost gets teh total cost of 'inserting' the data in a block
getblockInsertDataCost :: BlockData -> Double
getblockInsertDataCost (_, characterDataVV, charInfoV) =
    V.sum $ fmap (getLeafInsertCost charInfoV) characterDataVV

-- | getLeafInsertCost is the cost or ortiginating or 'inserting' leaf data
-- for all characters in a block
getLeafInsertCost :: V.Vector CharInfo -> V.Vector CharacterData -> Double
getLeafInsertCost charInfoV charDataV = 
    V.sum $ V.zipWith getCharacterInsertCost charDataV charInfoV

-- | getCharacterInsertCost takes a character and characterInfo and retujrns origination/insert cost for the character
getCharacterInsertCost :: CharacterData -> CharInfo -> Double
getCharacterInsertCost inChar charInfo =
    let localCharType = charType charInfo
        thisWeight = weight charInfo
        inDelCost = (costMatrix charInfo) S.! (0, (length (alphabet charInfo) - 1))
    in
    if localCharType == Add then thisWeight * (fromIntegral $ V.length $ GU.fst3 $ rangePrelim inChar)
    else if localCharType == NonAdd then thisWeight * (fromIntegral $ V.length $ GU.fst3 $ stateBVPrelim inChar) 
    else if localCharType == Matrix then thisWeight * (fromIntegral $ V.length $ matrixStatesPrelim inChar)
    else if localCharType == SlimSeq || localCharType == NucSeq then thisWeight * (fromIntegral inDelCost) * (fromIntegral $ SV.length $ slimPrelim inChar)
    else if localCharType == WideSeq || localCharType ==  AminoSeq then thisWeight * (fromIntegral inDelCost) * (fromIntegral $ UV.length $ widePrelim inChar)
    else if localCharType == HugeSeq then thisWeight * (fromIntegral inDelCost) * (fromIntegral $ V.length $ hugePrelim inChar)
    else error ("Character type unimplemented : " ++ show localCharType)


-- | splitSequence takes a ShortText divider and splits a list of ShortText on 
-- that ShortText divider consuming it akin to Text.splitOn
splitSequence :: ST.ShortText -> [ST.ShortText] -> [[ST.ShortText]]
splitSequence partitionST stList =
    if null stList then []
    else
        let firstPart = takeWhile (/= partitionST) stList
            restList = dropWhile (/= partitionST) stList
        in
        if restList == [partitionST] then firstPart : [[ST.fromString "#"]]
        else if not $ null restList then firstPart : splitSequence partitionST (tail restList)
        else [firstPart]



-- | dynamicCharacterTo3Vector takes a DYnamicCharacter and returns three Vectors
dynamicCharacterTo3Vector :: DynamicCharacter -> (Word, V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector)
dynamicCharacterTo3Vector (Missing x) = (x, V.empty, V.empty, V.empty)
dynamicCharacterTo3Vector (DC x) =
    let neVect = V.fromList $ toList x
        (a,b,c) = V.unzip3 neVect
    in
    (0 :: Word, a, b, c)


convertVectorToDynamicCharacter :: V.Vector BV.BitVector -> DynamicCharacter
convertVectorToDynamicCharacter inVector =
    let lenAlph = BV.dimension $ V.head inVector
        arbitraryAlphabet = fromSymbols $ show <$> 0 :| [1 .. lenAlph - 1]
    in
    encodeStream arbitraryAlphabet $ fmap (NE.fromList . f 0 . BV.toBits) . NE.fromList $ toList  inVector
    where
        f :: Word -> [Bool] -> [String]
        f _ [] = []
        f n (x:xs)
            | x = show n : f (n+1) xs
            | otherwise = f (n+1) xs


-- See Bio.DynamicCharacter.decodeState for a better implementation for dynamic character elements
bitVectToCharState :: Bits b => Alphabet String -> b -> String
bitVectToCharState localAlphabet bitValue = L.intercalate "," $ foldr pollSymbol mempty indices 
  where
    indices = [ 0 .. len - 1 ]
    len = length vec
    vec = alphabetSymbols localAlphabet
    pollSymbol i polled
      | bitValue `testBit` i = (vec V.! i) : polled
      | otherwise         = polled


-- bitVectToCharState  takes a bit vector representation and returns a list states as integers
bitVectToCharState' :: (Bits b) => [String] -> b -> String
bitVectToCharState' localAlphabet bitValue
  | isAlphabetDna       hereAlphabet = fold $ iupacToDna       BM.!> observedSymbols
  | isAlphabetAminoAcid hereAlphabet = fold $ iupacToAminoAcid BM.!> observedSymbols
  | otherwise = L.intercalate "," $ toList observedSymbols
  where
      hereAlphabet = fromSymbols localAlphabet
      symbolCountH = length localAlphabet
      observedSymbols
        = NE.fromList
            $ foldMap
                (\ i -> [localAlphabet !! i | bitValue `testBit` i])
                [0 .. symbolCountH - 1]

{-}
    --if DNA use IUPAC ambiguity codes
    if localAlphabet == ["A","C","G","T","-"] then 
        let numBit = (BV.toUnsignedNumber inBit) 
        in 
        if numBit == 1 then "A"
        else if numBit == 2 then "C"
        else if numBit == 4 then "G"
        else if numBit == 8 then "T"
        else if numBit == 16 then "-"
        -- IUPAC ambiguity
        else if numBit == 5 then "R"
        else if numBit == 10 then "Y"
        else if numBit == 3 then "M"
        else if numBit == 9 then "W"
        else if numBit == 6 then "S"
        else if numBit == 12 then "K"
        else if numBit == 14 then "B"
        else if numBit == 13 then "D"
        else if numBit == 11 then "H"
        else if numBit == 7 then "V"
        else if numBit == 15 then "N"
        -- indel ambiguity
        else if numBit == 17 then "A|"
        else if numBit == 18 then "C|"
        else if numBit == 20 then "G|"
        else if numBit == 24 then "T|"
        else if numBit == 21 then "R|"
        else if numBit == 26 then "Y|"
        else if numBit == 19 then "M|"
        else if numBit == 25 then "W|"
        else if numBit == 22 then "S|"
        else if numBit == 28 then "K|"
        else if numBit == 30 then "B|"
        else if numBit == 29 then "D|"
        else if numBit == 27 then "H|"
        else if numBit == 23 then "V|"
        else if numBit == 31 then "?"
        else error ("Unrecognized bit coded ambiguity group " ++ show numBit) 
    else     
        let bitBoolPairList = zip (BV.toBits inBit) localAlphabet
            (_, stateList) = unzip $ filter ((==True).fst) bitBoolPairList
            in
            intercalate "," stateList


-}

-- | filledDataFields takes rawData and checks taxon to see what percent
-- "characters" are found.
-- call with (0,0)
filledDataFields :: (Int, Int) -> TermData -> (NameText, Int, Int)
filledDataFields (hasData, totalData) (taxName, taxData)
  | null taxData = (taxName, hasData, totalData)
  | ST.length (head taxData) == 0 = filledDataFields (hasData, 1 + totalData) (taxName, tail taxData)
  | otherwise = filledDataFields (1 + hasData, 1 + totalData) (taxName, tail taxData)

-- | stripComments removes all lines that being with "--" haskell stype comments
-- needs to be reversed on return to maintain order.
stripComments :: [String] -> [String]
stripComments inStringList =
    if null inStringList then []
    else
        let strippedLine = GU.stripString $ head inStringList
        in
        if null strippedLine then stripComments $ tail inStringList
        else if length strippedLine < 2 then strippedLine : stripComments (tail inStringList)
        else
            if "--" == take 2 strippedLine then stripComments $ tail inStringList
            else strippedLine : stripComments (tail inStringList)

-- | getDecoratedGraphBlockCharInformation takes decorated graph and reports number of blosk and size of each
getDecoratedGraphBlockCharInformation :: DecoratedGraph -> ((Int, Int), [V.Vector Int])
getDecoratedGraphBlockCharInformation inGraph =
    if LG.isEmpty inGraph then ((0,0), [])
    else
        -- get a vertices from graph and take their information
        let inVertDataList = fmap (vertData . snd) (LG.labNodes inGraph)
            blockNumMax = maximum $ fmap length inVertDataList
            blocknumMin = minimum $ fmap length inVertDataList
            blockLengthList = fmap (fmap length) inVertDataList
        in
        ((blockNumMax, blocknumMin), blockLengthList)

-- | vectMaybeHead takes a vector and returns JUst V.head if not V.empty
-- Nothing otherwise
vectMaybeHead :: V.Vector a -> Maybe a
vectMaybeHead inVect =
    if V.null inVect then Nothing
    else Just (V.head inVect)

-- vectResolveMaybe takes a Vector of Maybe a 
-- and returns Just a or V.empty
vectResolveMaybe :: V.Vector (Maybe a) -> V.Vector a
vectResolveMaybe inVect =
    trace ("VRM " ++ show (length inVect)) (
    if isNothing (V.head inVect) then V.empty
    else V.singleton $ fromJust $ V.head inVect
    )

-- | getNumberNonExactCharacters takes processed data and returns the number of non-exact characters
-- ised to special case datasets with limited non-exact characters
getNumberNonExactCharacters :: V.Vector BlockData -> Int
getNumberNonExactCharacters blockDataVect =
    if V.null blockDataVect then 0
    else
        let firstBlock = GU.thd3 $ V.head blockDataVect
            characterTypes = V.map charType firstBlock
            nonExactChars = length $ V.filter (== True) $ V.map (`elem` nonExactCharacterTypes) characterTypes
        in
        nonExactChars + getNumberNonExactCharacters (V.tail blockDataVect)

-- | getNumberExactCharacters takes processed data and returns the number of non-exact characters
-- ised to special case datasets with limited non-exact characters
getNumberExactCharacters :: V.Vector BlockData -> Int
getNumberExactCharacters blockDataVect =
    if V.null blockDataVect then 0
    else
        let firstBlock = GU.thd3 $ V.head blockDataVect
            characterTypes = V.map charType firstBlock
            nonExactChars = length $ V.filter (== True) $ V.map (`elem` exactCharacterTypes) characterTypes
        in
        nonExactChars + getNumberExactCharacters (V.tail blockDataVect)

-- | splitBlockCharacters takes a block of characters (vector) and splits into two partitions of exact and non-exact characters
-- using accumulators
splitBlockCharacters :: V.Vector (V.Vector CharacterData)
                     -> V.Vector CharInfo
                     -> Int
                     -> [([CharacterData], CharInfo)]
                     -> [([CharacterData], CharInfo)]
                     -> (BlockData, BlockData)
splitBlockCharacters inDataVV inCharInfoV localIndex exactCharPairList nonExactCharPairList =
    if localIndex == V.length inCharInfoV then
        let (exactDataList, exactCharInfoList) = unzip exactCharPairList
            (nonExactDataList, nonExactCharInfoList) = unzip nonExactCharPairList
            newExactCharInfoVect = V.fromList $ reverse exactCharInfoList
            newNonExactCharInfoVect = V.fromList $ reverse nonExactCharInfoList
            newExactData = V.fromList $ fmap (V.fromList . reverse) (L.transpose exactDataList)
            newNonExactData = V.fromList $ fmap (V.fromList . reverse) (L.transpose nonExactDataList)
        in
        ((T.pack "ExactCharacters", newExactData, newExactCharInfoVect), (T.pack "Non-ExactCharacters", newNonExactData, newNonExactCharInfoVect))
    else
        let localCharacterType = charType (inCharInfoV V.! localIndex)
            thisCharacterData = V.toList $ fmap (V.! localIndex) inDataVV
            newPair = (thisCharacterData, inCharInfoV V.! localIndex)
        in
        if localCharacterType `elem` exactCharacterTypes then
            splitBlockCharacters inDataVV inCharInfoV (localIndex + 1) (newPair : exactCharPairList) nonExactCharPairList
        else if localCharacterType `elem` nonExactCharacterTypes then
            splitBlockCharacters inDataVV inCharInfoV (localIndex + 1) exactCharPairList (newPair : nonExactCharPairList)
        else error ("Unrecongized/implemented character type: " ++ show localCharacterType)



-- | safeVectorHead safe vector head--throws error if null
safeVectorHead :: V.Vector a -> a
safeVectorHead inVect =
    if V.null inVect then error "Empty vector in safeVectorHead"
    else V.head inVect


-- | checkCommandArgs takes comamnd and args and verifies that they are in list
checkCommandArgs :: String -> [String] -> [String] -> Bool
checkCommandArgs commandString commandList permittedList =
    null commandList || (
    let firstCommand = head commandList
        foundCommand = firstCommand `elem` permittedList
    in
    if foundCommand then checkCommandArgs commandString (tail commandList) permittedList
    else
        let errorMatch = snd $ GU.getBestMatch (maxBound :: Int ,"no suggestion") permittedList firstCommand
        in
        errorWithoutStackTrace ("\nError: Unrecognized '"++ commandString ++"' option. By \'" ++ firstCommand ++ "\' did you mean \'" ++ errorMatch ++ "\'?\n"))

-- | get leftRightChilLabelBV takes a pair of vertex labels and returns left and right
-- based on their bitvector representation.  This ensures left/right consistancey in
-- pre and postoder passes, and with bitvectors of leaves determined by data hash,
-- ensures label invariance with repect to leaves 
-- larger bitvector goes second (bigge)
leftRightChildLabelBV :: (VertexInfo, VertexInfo) -> (VertexInfo, VertexInfo)
leftRightChildLabelBV inPair@(firstNode, secondNode) =
    let firstLabel  = bvLabel firstNode
        secondLabel = bvLabel secondNode
    in
    if firstLabel > secondLabel then (secondNode, firstNode)
    else inPair

-- | get leftRightChildLabelBVNode takes a pair ofnodes and returns left and right
-- based on their bitvector representation.  This ensures left/right consistancey in
-- pre and postoder passes, and with bitvectors of leaves determined by data hash,
-- ensures label invariance with repect to leaves 
-- larger bitvector goes second (bigge)
leftRightChildLabelBVNode :: (LG.LNode VertexInfo, LG.LNode VertexInfo) -> (LG.LNode VertexInfo, LG.LNode VertexInfo)
leftRightChildLabelBVNode inPair@(firstNode, secondNode) =
    let firstLabel  = bvLabel $ snd firstNode
        secondLabel = bvLabel $ snd secondNode
    in
    if firstLabel > secondLabel then (secondNode, firstNode)
    else inPair


{-
-- | getBridgeList takes a Decorated graph and returns the list of edges that, if removed, will
-- split the graph, that is increase the number of unconnected components
-- this is useful in graph rearrangement
-- root connected edges are not returned-- but in a general graph context can be bridges
getBridgeList :: DecoratedGraph -> [LG.LEdge EdgeInfo]
getBridgeList inGraph =
    if LG.isEmpty inGraph then []
    else
        let vertexList = LG.labNodes inGraph
            labEdgeList = LG.labEdges inGraph
            vertBVList = fmap (bvLabel . snd) networkVertexList
            vertPairVect = V.fromList $ L.sortOn fst $ zip (fmap fst vertexList) vertBVList
            (_, _, _, networkVertexList) = LG.splitVertexList inGraph
            netVertBVList = fmap (bvLabel . snd) networkVertexList
            -- netVertPairList = zip (fmap fst networkVertexList) netVertBVList
            bridgeList = getBridgeList' vertPairVect netVertBVList labEdgeList

        in
        if null networkVertexList then labEdgeList
        else bridgeList

-- | getBridgeList takes a vector of (vertex, bitvector label) pairs, a list of network 
-- vertex bitvector labels, and a list of labelled edge and returns a list of bridge edges 
-- checks whether for each edge (u,v), the bitvector labels of all the network nodes are 
-- `compatible' with bit vector of vertex v.  a, and b are compatible if a .&. b = a, b, or empty
getBridgeList' :: V.Vector (Int, BV.BitVector) -> [BV.BitVector] -> [LG.LEdge EdgeInfo] -> [LG.LEdge EdgeInfo]
getBridgeList' vertexPairVect netVertBVList inEdgeList =
    if null inEdgeList then []
    else
        let firstEdge@(_, vVertex, _) = head inEdgeList
            isBridge = GU.isBVCompatible (snd $ vertexPairVect V.! vVertex) netVertBVList
        in
        if isBridge then firstEdge : getBridgeList' vertexPairVect netVertBVList (tail inEdgeList)
        else getBridgeList' vertexPairVect netVertBVList (tail inEdgeList)
-}

-- | prettyPrintVertexInfo returns a string with formated version of 
-- vertex info
prettyPrintVertexInfo :: VertexInfo -> String
prettyPrintVertexInfo inVertData =
    let zerothPart = "Vertex name " ++ T.unpack (vertName inVertData) ++ " Index " ++ show (index inVertData)
        firstPart = "\n\tBitVector (as number) " ++ show (BV.toUnsignedNumber $ bvLabel inVertData)
        secondPart = "\n\tParents " ++ show (parents inVertData) ++ " Children " ++ show (children inVertData)
        thirdPart = "\n\tType " ++ show (nodeType inVertData) ++ " Local Cost " ++ show (vertexCost inVertData) ++ " SubGraph Cost " ++ show (subGraphCost inVertData)
        fourthPart = "\n\tData Blocks: " ++ show (V.length $ vertData inVertData) ++ " Characters (by block) " ++ show (V.length <$> vertData inVertData)
        fifthPart = "\n\t" ++ show (vertData inVertData)
    in
    zerothPart ++ firstPart ++ secondPart ++ thirdPart ++ fourthPart ++ fifthPart


{-# INLINEABLE encodeState #-}
{-# SPECIALISE encodeState :: (Bits b, Foldable f, Ord s) => Alphabet s         -> (Word -> b        ) -> f s           -> b         #-}
{-# SPECIALISE encodeState :: (Bits b,             Ord s) => Alphabet s         -> (Word -> b        ) -> Set s         -> b         #-}
{-# SPECIALISE encodeState ::  Bits b                     => Alphabet String    -> (Word -> b        ) -> Set String    -> b         #-}
{-# SPECIALISE encodeState ::  Bits b                     => Alphabet ShortText -> (Word -> b        ) -> Set ShortText -> b         #-}
{-# SPECIALISE encodeState :: (        Foldable f, Ord s) => Alphabet s         -> (Word -> BitVector) -> f s           -> BitVector #-}
{-# SPECIALISE encodeState ::                      Ord s  => Alphabet s         -> (Word -> BitVector) -> Set s         -> BitVector #-}
{-# SPECIALISE encodeState ::                                Alphabet String    -> (Word -> BitVector) -> Set String    -> BitVector #-}
{-# SPECIALISE encodeState ::                                Alphabet ShortText -> (Word -> BitVector) -> Set ShortText -> BitVector #-}
encodeState
  :: ( Bits e
     , Foldable f
     , Ord s
     )
  => Alphabet s  -- ^ Alphabet of symbols
  -> (Word -> e) -- ^ Constructor for an empty element, taking the alphabet size
  -> f s         -- ^ ambiguity groups of symbols
  -> e           -- ^ Encoded dynamic character element
encodeState alphabet' f symbols = getSubsetIndex alphabet' symbolsSet emptyElement
  where
    emptyElement = f . toEnum $ length alphabet'
    symbolsSet   = Set.fromList $ toList symbols


-- | add3 adds three values
add3 :: (Num a) => a -> a -> a -> a
add3 x y z = x + y + z
