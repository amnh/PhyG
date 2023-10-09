{- |
Module      :  Medians.hs
Description :  Module specifying data type medians
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

{--
TODO:

  Parallelize  median2Vect
--}

module GraphOptimization.Medians  ( median2
                                  , median2Single
                                  , median2NonExact
                                  -- , median2SingleNonExact
                                  , median2StaticIA
                                  , makeIAUnionPrelimLeaf
                                  , makeIAPrelimCharacter
                                  , makeIAFinalCharacter
                                  -- , createUngappedMedianSequence
                                  , intervalAdd
                                  , interUnion
                                  , getNewRange
                                  , addMatrix
                                  , getDOMedian
                                  , getPreAligned2Median
                                  , getDOMedianCharInfo
                                  , getNoGapPrelimContext
                                  , pairwiseDO
                                  , makeDynamicCharacterFromSingleVector
                                  , makeEdgeData
                                  , createEdgeUnionOverBlocks
                                  , get2WaySlim
                                  , get2WayWideHuge
                                  , getFinal3WaySlim
                                  , getFinal3WayWideHuge
                                  , generalSequenceDiff
                                  , union2Single
                                  , distance2Unions
                                  ) where

import Bio.DynamicCharacter
import Bio.DynamicCharacter.Element
import Data.Alphabet
import Data.BitVector.LittleEndian qualified as BV
import Data.Bits
import Data.Foldable
import Data.List qualified as L
import Data.Maybe
import Data.MetricRepresentation qualified  as MR
import Data.TCM.Dense qualified as TCMD
import Data.Vector qualified as V
import Data.Vector.Generic qualified  as GV
import Data.Vector.Storable qualified as SV
import Data.Word
import DirectOptimization.Pairwise
import Foreign.C.Types (CUInt)
import GeneralUtilities
import Input.BitPack qualified as BP
import SymMatrix qualified as S
import Types.Types
import Utilities.LocalGraph qualified as LG
--import           Debug.Trace


-- | makeDynamicCharacterFromSingleVector takes a single vector (usually a 'final' state)
-- and returns a dynamic character that canbe used with other functions
makeDynamicCharacterFromSingleVector :: (GV.Vector v a) => v a -> (v a, v a, v a)
makeDynamicCharacterFromSingleVector dc = unsafeCharacterBuiltByST (toEnum $ GV.length dc) $ \dc' -> GV.imapM_ (\k v -> setAlign dc' k v v v) dc

-- | median2 takes the vectors of characters and applies median2Single to each
-- character
-- for parallel fmap over all then parallelized by type and sequences
-- used for distances and post-order assignments
median2 ::   V.Vector CharacterData -> V.Vector CharacterData -> V.Vector CharInfo -> V.Vector (CharacterData, VertexCost)
median2 = V.zipWith3 (median2Single False)


-- | median2NonExact takes the vectors of characters and applies median2NonExact to each
-- character for parallel fmap over all then parallelized by type and sequences
-- this only reoptimized the nonexact characters (sequence characters for now, perhpas otehrs later)
-- and takes the existing optimization for exact (Add, NonAdd, Matrix) for the others.
median2NonExact :: V.Vector CharacterData -> V.Vector CharacterData -> V.Vector CharInfo -> V.Vector (CharacterData, VertexCost)
median2NonExact = V.zipWith3 median2SingleNonExact

-- | median2StaticIA takes the vectors of characters and applies median2SingleStaticIA to each
-- character for parallel fmap over all then parallelized by type and sequences
-- this reoptimized only IA fields for the nonexact characters (sequence characters for now, perhpas others later)
-- and takes the existing optimization for exact (Add, NonAdd, Matrix) for the others.
median2StaticIA :: V.Vector CharacterData -> V.Vector CharacterData -> V.Vector CharInfo -> V.Vector (CharacterData, VertexCost)
median2StaticIA = V.zipWith3 (median2Single True)

-- | median2Single takes character data and returns median character and cost
-- median2single assumes that the character vectors in the various states are the same length
-- that is--all leaves (hence other vertices later) have the same number of each type of character
-- used for post-order assignments
-- this is from preliminary states
-- staticIA for dynm,aic assumes all same length
-- PMDL costs are claculate by type--additive by conversion to non-additive --but if states> 129 worn't do it so warning in docs
-- bp2,4,5,8,64, nonadd are by weights vis set command, matrix, sequence are set by tcm with non-zero diagnonal
median2Single :: Bool -> CharacterData -> CharacterData -> CharInfo -> (CharacterData, VertexCost)
median2Single staticIA firstVertChar secondVertChar inCharInfo =
    let thisType    = charType inCharInfo
        thisWeight  = weight inCharInfo
        thisMatrix  = costMatrix inCharInfo
        thisSlimTCM = slimTCM inCharInfo
        thisWideTCM = wideTCM inCharInfo
        thisHugeTCM = hugeTCM inCharInfo
        thisActive  = activity inCharInfo
        thisNoChangeCost = noChangeCost inCharInfo
        thisChangeCost = changeCost inCharInfo
    in
    if not thisActive then (firstVertChar, 0)
    else if thisType == Add then
        let newCharVect = intervalAdd thisWeight firstVertChar secondVertChar
        in
        (newCharVect, localCost  newCharVect)

    else if thisType == NonAdd then
        let newCharVect = interUnion thisWeight (thisNoChangeCost, thisChangeCost) firstVertChar secondVertChar
        in
        (newCharVect, localCost  newCharVect)

    else if thisType `elem` packedNonAddTypes then
        --assumes all weight 1
        let newCharVect = BP.median2Packed thisType thisWeight (thisNoChangeCost, thisChangeCost) firstVertChar secondVertChar
        in
        (newCharVect, localCost  newCharVect)

    else if thisType == Matrix then
      let newCharVect = addMatrix thisWeight thisMatrix firstVertChar secondVertChar
        in
        (newCharVect, localCost  newCharVect)

    else if thisType `elem` prealignedCharacterTypes then
      let newCharVect = getPreAligned2Median inCharInfo emptyCharacter firstVertChar secondVertChar
      in
      (newCharVect, localCost  newCharVect)

    else if thisType `elem` nonExactCharacterTypes then
      let newCharVect = if staticIA then makeIAPrelimCharacter inCharInfo emptyCharacter firstVertChar secondVertChar
                        else getDOMedian thisWeight thisMatrix thisSlimTCM thisWideTCM thisHugeTCM thisType firstVertChar secondVertChar
      in
      -- trace ("M2S: " <> (show $ localCost  newCharVect))
      (newCharVect, localCost  newCharVect)

    else error ("Character type " <> show thisType <> " unrecognized/not implemented")

-- | median2SingleNonExact takes character data and returns median character and cost
-- median2single assumes that the character vectors in the various states are the same length
-- that is--all leaves (hencee other vertices later) have the same number of each type of character
-- this only reoptimized the nonexact characters (sequence characters for now, perhpas otehrs later)
-- and skips optimization placing a dummy value exact (Add, NonAdd, Matrix) for the others.
median2SingleNonExact :: CharacterData -> CharacterData -> CharInfo -> (CharacterData, VertexCost)
median2SingleNonExact firstVertChar secondVertChar inCharInfo =
    let thisType    = charType inCharInfo
        thisWeight  = weight inCharInfo
        thisMatrix  = costMatrix inCharInfo
        thisSlimTCM = slimTCM inCharInfo
        thisWideTCM = wideTCM inCharInfo
        thisHugeTCM = hugeTCM inCharInfo
        thisActive  = activity inCharInfo
        dummyStaticCharacter = emptyCharacter
    in
    if not thisActive || (thisType `elem` exactCharacterTypes) then (dummyStaticCharacter, 0)
    else if thisType `elem` prealignedCharacterTypes then
         let newCharVect = getPreAligned2Median inCharInfo dummyStaticCharacter firstVertChar secondVertChar
         in
         -- trace ("M2S:" <> (show $ localCost  newCharVect) <> (show (firstVertChar, secondVertChar)))
         -- trace ("M2SNEP: " <> (show thisType))
         (newCharVect, localCost  newCharVect)

    else if thisType `elem` nonExactCharacterTypes then
         let newCharVect = getDOMedian thisWeight thisMatrix thisSlimTCM thisWideTCM thisHugeTCM thisType firstVertChar secondVertChar
         in
         --trace ("M2SNE: " <> (show thisType) <> (show $ localCost  newCharVect))
         (newCharVect, localCost  newCharVect)

    else error ("Character type " <> show thisType <> " unrecognized/not implemented")

-- | distance2Unions is a block wrapper around  distance2UnionsBlock
-- really only to get union distance--keeping union states in there in csae needed later
distance2Unions ::VertexBlockData ->  VertexBlockData -> V.Vector (V.Vector CharInfo) -> (VertexBlockData, VertexCost)
distance2Unions firstBlock secondBlock charInfoVV =
    let (newBlockV, newCostV) = V.unzip $ V.zipWith3 distance2UnionsBlock firstBlock secondBlock charInfoVV
    in
    (newBlockV, V.sum newCostV)

-- | distance2UnionsBlock is a block wrapper around  distance2UnionsCharacter
distance2UnionsBlock :: V.Vector CharacterData ->  V.Vector CharacterData -> V.Vector CharInfo -> (V.Vector CharacterData, VertexCost)
distance2UnionsBlock firstBlock secondBlock charInfoV =
    let (newBlockV, newCostV) = V.unzip $ V.zipWith3 distance2UnionsCharacter firstBlock secondBlock charInfoV
    in
    (newBlockV, V.sum newCostV)

-- | distance2UnionsCharacter takes character data and returns median of union characters and cost
-- median2Unions assumes that the character vectors in the various states are the same length
-- this is from union states
-- assumes all same length
-- PMDL costs are claculate by type--additive by conversion to non-additive --but if states> 129 won't do it so warning in docs
-- bp2,4,5,8,64, nonadd are by weights vis set command, matrix, sequence are set by tcm with non-zero diagnonal

-- wrong for prealigned--needs to DO
    -- 1) Check length/composition union fields
    -- 2) check cost for each type
    -- 3) use DO for alignments even on unoin fields (triple)
-- as a result can use non-alignred dynamic characters as unions for comparison
distance2UnionsCharacter :: CharacterData -> CharacterData -> CharInfo -> (CharacterData, VertexCost)
distance2UnionsCharacter firstVertChar secondVertChar inCharInfo =
    let thisType    = charType inCharInfo
        thisWeight  = weight inCharInfo
        thisMatrix  = costMatrix inCharInfo
        thisActive  = activity inCharInfo
        thisNoChangeCost = noChangeCost inCharInfo
        thisChangeCost = changeCost inCharInfo
    in
    if not thisActive then (firstVertChar, 0)
    else if thisType == Add then
        let newCharVect = intervalAddUnionField thisWeight firstVertChar secondVertChar
        in
        (newCharVect, localCost  newCharVect)

    else if thisType == NonAdd then
        let newCharVect = interUnionUnionField thisWeight (thisNoChangeCost, thisChangeCost) firstVertChar secondVertChar
        in
        (newCharVect, localCost  newCharVect)

    else if thisType `elem` packedNonAddTypes then
        --assumes all weight 1
        let newCharVect = BP.median2PackedUnionField thisType thisWeight (thisNoChangeCost, thisChangeCost) firstVertChar secondVertChar
        in
        (newCharVect, localCost  newCharVect)

    else if thisType == Matrix then
      let newCharVect = addMatrixUnionField thisWeight thisMatrix firstVertChar secondVertChar
        in
        (newCharVect, localCost  newCharVect)

    else if thisType `elem` prealignedCharacterTypes then
      let newCharVect = getPreAligned2MedianUnionFields inCharInfo emptyCharacter firstVertChar secondVertChar
      in
      (newCharVect, localCost  newCharVect)

    else if thisType `elem` nonExactCharacterTypes then
      let newCharVect = getDOMedianCharInfoUnion inCharInfo firstVertChar secondVertChar
      -- let newCharVect = getNonExactUnionFields inCharInfo emptyCharacter firstVertChar secondVertChar
      in
      -- trace ("M2S: " <> (show $ localCost  newCharVect))
     (newCharVect, localCost  newCharVect)

    else error ("Character type " <> show thisType <> " unrecognized/not implemented")

-- | localOr wrapper for BV.or for vector elements
localOr :: BV.BitVector -> BV.BitVector -> BV.BitVector
localOr lBV rBV = lBV .|. rBV

-- | interUnionBV takes two bitvectors and returns new state, nochange number (1 or 0), change number (0 or 1)
interUnionBV :: BV.BitVector -> BV.BitVector -> (BV.BitVector, Int, Int)
interUnionBV leftBV rightBV =
    if BV.isZeroVector (leftBV .&. rightBV)  then (leftBV .|. rightBV, 0, 1)
    else (leftBV .&. rightBV, 1, 0)

-- | interUnion takes two non-additive chars and creates newCharcter as 2-median
-- in post-order pass to create preliminary states assignment
-- assumes a single weight for all
-- performs two passes though chars to get cost of assignments
-- snd3 $ rangePrelim left/rightChar due to triple in prelim
-- could have done noChnageCoast/Chaneg cost with length subtraction but very small issue in real use since
-- only for nonadd characters with > 64 states.
interUnion :: Double -> (Double, Double) -> CharacterData -> CharacterData -> CharacterData
interUnion thisWeight (lNoChangeCost, lChangeCost) leftChar rightChar =
    let (newStateVect, noChangeCostVect, changeCostVect) = V.unzip3 $ V.zipWith interUnionBV (snd3 $ stateBVPrelim leftChar) (snd3 $ stateBVPrelim rightChar)
        newCost = thisWeight * ((lNoChangeCost * fromIntegral (V.sum noChangeCostVect)) + (lChangeCost * fromIntegral (V.sum changeCostVect)))
        newCharacter = emptyCharacter { stateBVPrelim = (snd3 $ stateBVPrelim leftChar, newStateVect, snd3 $ stateBVPrelim rightChar)
                                      , localCost = newCost
                                      , globalCost = newCost + globalCost leftChar + globalCost rightChar
                                      }
    in
    --trace ("NonAdditive: " <> (show numUnions) <> " " <> (show newCost) <> "\t" <> (show $ stateBVPrelim leftChar) <> "\t" <> (show $ stateBVPrelim rightChar) <> "\t"
    --   <> (show intersectVect) <> "\t" <> (show unionVect) <> "\t" <> (show newStateVect))
    newCharacter

-- | interUnionUnionField takes two non-additive chars and creates newCharcter as 2-median
-- in post-order pass to create union states assignment and cost
-- assumes a single weight for all
-- performs two passes though chars to get cost of assignments
-- snd3 $ rangePrelim left/rightChar due to triple in prelim
-- could have done noChnageCoast/Chaneg cost with length subtraction but very small issue in real use since
-- only for nonadd characters with > 64 states.
interUnionUnionField :: Double -> (Double, Double) -> CharacterData -> CharacterData -> CharacterData
interUnionUnionField thisWeight (lNoChangeCost, lChangeCost) leftChar rightChar =
    let (newStateVect, noChangeCostVect, changeCostVect) = V.unzip3 $ V.zipWith interUnionBV (stateBVUnion leftChar) (stateBVUnion rightChar)
        newCost = thisWeight * ((lNoChangeCost * fromIntegral (V.sum noChangeCostVect)) + (lChangeCost * fromIntegral (V.sum changeCostVect)))
        newCharacter = emptyCharacter { stateBVUnion = newStateVect
                                      , localCost = newCost
                                      , globalCost = newCost + globalCost leftChar + globalCost rightChar
                                      }
    in
    --trace ("NonAdditive: " <> (show numUnions) <> " " <> (show newCost) <> "\t" <> (show $ stateBVPrelim leftChar) <> "\t" <> (show $ stateBVPrelim rightChar) <> "\t"
    --   <> (show intersectVect) <> "\t" <> (show unionVect) <> "\t" <> (show newStateVect))
    newCharacter


-- | localUnion takes two non-additive chars and creates newCharcter as 2-union/or
-- assumes a single weight for all
-- performs single
-- bsaed on final states
localUnion :: CharacterData -> CharacterData -> CharacterData
localUnion leftChar rightChar =
    let unionVect = V.zipWith localOr (stateBVFinal leftChar) (stateBVFinal rightChar)
        newCharacter = emptyCharacter { stateBVPrelim = (unionVect, unionVect, unionVect)
                                      , stateBVFinal = unionVect
                                      }
    in
    --trace ("NonAdditive: " <> (show numUnions) <> " " <> (show newCost) <> "\t" <> (show $ stateBVPrelim leftChar) <> "\t" <> (show $ stateBVPrelim rightChar) <> "\t"
    --   <> (show intersectVect) <> "\t" <> (show unionVect) <> "\t" <> (show newStateVect))
    newCharacter

-- | getNewRange takes min and max range of two additive characters and returns
-- a triple of (newMin, newMax, Cost)
getNewRange :: Int-> Int -> Int -> Int -> (Int, Int, Int)
getNewRange lMin lMax rMin rMax =
    -- subset
    if (rMin >= lMin) && (rMax <= lMax) then (rMin, rMax, 0)
    else if  (lMin >= rMin) && (lMax <= rMax) then (lMin, lMax, 0)
    -- overlaps
    else if (rMin >= lMin) && (rMax >= lMax) && (rMin <= lMax) then (rMin, lMax,0)
    else if (lMin >= rMin) && (lMax >= rMax) && (lMin <= rMax) then (lMin, rMax,0)
    -- newInterval
    else if lMax <= rMin then (lMax, rMin, rMin - lMax)
    else if rMax <= lMin then (rMax, lMin, lMin - rMax)
    else error ("This can't happen " <> show (lMin, lMax, rMin, rMax))

-- | intervalAdd takes two additive chars and creates newCharcter as 2-median
-- in post-order pass to create preliminary states assignment
-- assumes a single weight for all
-- snd3 $ rangePrelim left/rightChar due to triple in prelim
intervalAdd :: Double -> CharacterData -> CharacterData -> CharacterData
intervalAdd thisWeight leftChar rightChar =
    let newRangeCosts = V.zipWith4 getNewRange (V.map fst $ snd3 $ rangePrelim leftChar) (V.map snd $ snd3 $ rangePrelim leftChar) (V.map fst $ snd3 $ rangePrelim rightChar) (V.map snd $ snd3 $ rangePrelim rightChar)
        newMinRange = V.map fst3 newRangeCosts
        newMaxRange = V.map snd3 newRangeCosts
        newCost = thisWeight * fromIntegral (V.sum $ V.map thd3 newRangeCosts)
        newCharacter = emptyCharacter { rangePrelim = (snd3 $ rangePrelim leftChar, V.zip newMinRange newMaxRange, snd3 $ rangePrelim rightChar)
                                      , localCost = newCost
                                      , globalCost = newCost + globalCost leftChar + globalCost rightChar
                                      }
    in
    newCharacter


-- | intervalAddUnionField takes two additive chars and creates newCharcter as 2-median
-- in post-order pass to create union states assignment and cost
-- assumes a single weight for all
intervalAddUnionField :: Double -> CharacterData -> CharacterData -> CharacterData
intervalAddUnionField thisWeight leftChar rightChar =
    let newRangeCosts = V.zipWith4 getNewRange (V.map fst $ rangeUnion leftChar) (V.map snd $ rangeUnion leftChar) (V.map fst $ rangeUnion rightChar) (V.map snd $ rangeUnion rightChar)
        newMinRange = V.map fst3 newRangeCosts
        newMaxRange = V.map snd3 newRangeCosts
        newCost = thisWeight * fromIntegral (V.sum $ V.map thd3 newRangeCosts)
        newCharacter = emptyCharacter { rangeUnion = V.zip newMinRange newMaxRange
                                      , localCost = newCost
                                      , globalCost = newCost + globalCost leftChar + globalCost rightChar
                                      }
    in
    newCharacter

-- | getUnionRange takes min and max range of two additive characters and returns
-- a pair of (newMin, newMax)
getUnionRange :: Int -> Int -> Int -> Int -> (Int, Int)
getUnionRange lMin lMax rMin rMax =
    (min lMin rMin, max lMax rMax)

-- | intervalUnion takes two additive chars and creates newCharcter as 2-union
-- min of all lower, max of all higher
-- final states used and assigned to obthe prelim and final for use in swap delta
intervalUnion :: CharacterData -> CharacterData -> CharacterData
intervalUnion leftChar rightChar =
    let newRangeCosts = V.zipWith4 getUnionRange (V.map fst $ rangeFinal leftChar) (V.map snd $ rangeFinal leftChar) (V.map fst $ rangeFinal rightChar) (V.map snd $ rangeFinal rightChar)
        newMinRange = V.map fst newRangeCosts
        newMaxRange = V.map snd newRangeCosts
        newRange = V.zip newMinRange newMaxRange
        newCharacter = emptyCharacter { rangePrelim = (newRange, newRange, newRange)
                                      , rangeFinal = newRange
                                      }
    in
    --trace ("Additive: " <> (show newCost) <> "\t" <> (show $ rangeFinal leftChar) <> "\t" <> (show $ rangeFinal rightChar)
     --   <> (show newRangeCosts))
    newCharacter

-- | getMinCostStates takes cost matrix and vector of states (cost, _, _) and retuns a list of (totalCost, best child state)
getMinCostStates :: S.Matrix Int -> V.Vector MatrixTriple -> Int -> Int -> Int -> [(Int, ChildStateIndex)]-> Int -> [(Int, ChildStateIndex)]
getMinCostStates thisMatrix childVect bestCost numStates childState currentBestStates stateIndex =
   --trace (show thisMatrix <> "\n" <> (show  childVect) <> "\n" <> show (numStates, childState, stateIndex)) (
   if V.null childVect then reverse (filter ((== bestCost).fst) currentBestStates)
   else
      let (childCost, _, _)  = V.head childVect
          childStateCost = if childCost /= (maxBound :: Int) then childCost + (thisMatrix S.! (childState, stateIndex))
                           else (maxBound :: Int)
      in
      if childStateCost > bestCost then getMinCostStates thisMatrix (V.tail childVect) bestCost numStates (childState + 1) currentBestStates stateIndex
      else if childStateCost == bestCost then getMinCostStates thisMatrix (V.tail childVect) bestCost numStates (childState + 1) ((childStateCost, childState) : currentBestStates) stateIndex
      else getMinCostStates thisMatrix (V.tail childVect) childStateCost numStates (childState + 1) [(childStateCost, childState)] stateIndex
    -- )


-- | getNewVector takes the vector of states and costs from the child nodes and the
-- cost matrix and calculates a new vector n^2 in states
getNewVector :: S.Matrix Int -> Int -> (V.Vector MatrixTriple, V.Vector MatrixTriple) -> V.Vector MatrixTriple
getNewVector thisMatrix  numStates (lChild, rChild) =
  let newStates = [0..(numStates -1)]
      leftPairs = fmap (getMinCostStates thisMatrix lChild (maxBound :: Int) numStates 0 []) newStates
      rightPairs = fmap (getMinCostStates thisMatrix rChild (maxBound :: Int) numStates 0 []) newStates
      stateCosts = zipWith (+) (fmap (fst . head) leftPairs) (fmap (fst . head) rightPairs)
      newStateTripleList = zip3 stateCosts (fmap (fmap snd) leftPairs) (fmap (fmap snd) rightPairs)
  in
  V.fromList newStateTripleList

-- | addMatrix thisWeight thisMatrix firstVertChar secondVertChar matrix character
-- assumes each character has same cost matrix
-- Need to add approximation ala DO tcm lookup later
-- Local and global costs are based on current not necessaril;y optimal minimum cost states
addMatrix :: Double -> S.Matrix Int -> CharacterData -> CharacterData -> CharacterData
addMatrix thisWeight thisMatrix firstVertChar secondVertChar =
  if null thisMatrix then error "Null cost matrix in addMatrix"
  else
    let numStates = length thisMatrix
        initialMatrixVector = getNewVector thisMatrix numStates <$> V.zip (matrixStatesPrelim firstVertChar) (matrixStatesPrelim secondVertChar)
        initialCostVector = fmap (V.minimum . fmap fst3) initialMatrixVector
        newCost = thisWeight * fromIntegral (V.sum initialCostVector)
        newCharacter = emptyCharacter { matrixStatesPrelim = initialMatrixVector
                                      , localCost = newCost  - globalCost firstVertChar - globalCost secondVertChar
                                      , globalCost = newCost
                                      }
        in
        --trace ("Matrix: " <> (show newCost) <> "\n\t" <> (show $ matrixStatesPrelim firstVertChar)  <> "\n\t" <> (show $ matrixStatesPrelim secondVertChar) <>
        --  "\n\t" <> (show initialMatrixVector) <> "\n\t" <> (show initialCostVector))
        newCharacter

-- | addMatrixUnionField thisWeight thisMatrix firstVertChar secondVertChar matrix character
-- uinion fields
-- assumes each character has same cost matrix
-- Need to add approximation ala DO tcm lookup later
-- Local and global costs are based on current not necessarily optimal minimum cost states
addMatrixUnionField :: Double -> S.Matrix Int -> CharacterData -> CharacterData -> CharacterData
addMatrixUnionField thisWeight thisMatrix firstVertChar secondVertChar =
  if null thisMatrix then error "Null cost matrix in addMatrix"
  else
    let numStates = length thisMatrix
        initialMatrixVector = getNewVector thisMatrix numStates <$> V.zip (matrixStatesUnion firstVertChar) (matrixStatesUnion secondVertChar)
        initialCostVector = fmap (V.minimum . fmap fst3) initialMatrixVector
        newCost = thisWeight * fromIntegral (V.sum initialCostVector)
        newCharacter = emptyCharacter { matrixStatesUnion = initialMatrixVector
                                      , localCost = newCost  - globalCost firstVertChar - globalCost secondVertChar
                                      , globalCost = newCost
                                      }
        in
        newCharacter

-- | getUnionVector takes the vector of states and costs from two nodes
-- and sets the states with min cost in the two vertices and maxBound in other states
getUnionVector :: S.Matrix Int -> Int -> (V.Vector MatrixTriple, V.Vector MatrixTriple) -> V.Vector MatrixTriple
getUnionVector thisMatrix  numStates (lChild, rChild) =
  let newStates = [0..(numStates -1)]
      leftPairs = fmap (getMinCostStates thisMatrix lChild (maxBound :: Int) numStates 0 []) newStates
      rightPairs = fmap (getMinCostStates thisMatrix rChild (maxBound :: Int) numStates 0 []) newStates
      stateCosts = zipWith (+) (fmap (fst . head) leftPairs) (fmap (fst . head) rightPairs)
      minStateCost = minimum stateCosts
      stateCosts' = fmap (minOrMax minStateCost) stateCosts
      newStateTripleList = zip3 stateCosts' (fmap (fmap snd) leftPairs) (fmap (fmap snd) rightPairs)
  in
  V.fromList newStateTripleList
    where minOrMax minVal curVal = if curVal == minVal then minVal
                                    else maxBound :: Int

-- | unionMatrix  thisMatrix firstVertChar secondVertChar matrix character
-- assumes each character has same cost matrix
-- Need to add approximation ala DO tcm lookup later
-- Local and global costs are based on current not necessaril;y optimal minimum cost states
unionMatrix :: S.Matrix Int -> CharacterData -> CharacterData -> CharacterData
unionMatrix thisMatrix firstVertChar secondVertChar =
  if null thisMatrix then error "Null cost matrix in addMatrix"
  else
    let numStates = length thisMatrix
        initialMatrixVector = getUnionVector thisMatrix numStates <$> V.zip (matrixStatesFinal firstVertChar) (matrixStatesFinal secondVertChar)
        newCharacter = emptyCharacter { matrixStatesPrelim = initialMatrixVector
                                      , matrixStatesFinal = initialMatrixVector
                                      }
        in
        --trace ("Matrix: " <> (show newCost) <> "\n\t" <> (show $ matrixStatesPrelim firstVertChar)  <> "\n\t" <> (show $ matrixStatesPrelim secondVertChar) <>
        --  "\n\t" <> (show initialMatrixVector) <> "\n\t" <> (show initialCostVector))
        newCharacter

-- | pairwiseDO is a wrapper around slim/wise/hugeParwiseDO to allow direct call and return of
-- DO medians and cost.  This is used in final state assignment
pairwiseDO :: CharInfo
           -> (SlimDynamicCharacter, WideDynamicCharacter, HugeDynamicCharacter)
           -> (SlimDynamicCharacter, WideDynamicCharacter, HugeDynamicCharacter)
           -> (SlimDynamicCharacter, WideDynamicCharacter, HugeDynamicCharacter, Double)
pairwiseDO charInfo (slim1, wide1, huge1) (slim2, wide2, huge2) =
    let thisType = charType charInfo
    in
    if thisType `elem` [SlimSeq,   NucSeq]      then
        let (cost, r) = slimPairwiseDO (slimTCM charInfo) slim1 slim2
        in
        -- trace ("pDO:" <> (show (GV.length $ fst3 slim1, GV.length $ snd3 slim1)) <> " " <> (show (GV.length $ fst3 slim2, GV.length $ snd3 slim2)))
        (r, mempty, mempty, weight charInfo * fromIntegral cost)

    else if thisType `elem` [WideSeq, AminoSeq] then
        let coefficient = MR.minInDelCost (wideTCM charInfo)
            (cost, r) = widePairwiseDO coefficient (MR.retreivePairwiseTCM $ wideTCM charInfo) wide1 wide2
        in
        (mempty, r, mempty, weight charInfo * fromIntegral cost)

    else if thisType == HugeSeq           then
        let coefficient = MR.minInDelCost (hugeTCM charInfo)
            (cost, r) = hugePairwiseDO coefficient (MR.retreivePairwiseTCM $ hugeTCM charInfo) huge1 huge2
        in
        (mempty, mempty, r, weight charInfo * fromIntegral cost)

    else error $ fold ["Unrecognised character type '", show thisType, "'in a DYNAMIC character branch" ]

-- | getDOMedianCharInfoUnion  is a wrapper around getDOMedian with CharInfo-based interface
-- for union fields.
-- Strips out solo gapos (0/1) before DO step
getDOMedianCharInfoUnion :: CharInfo -> CharacterData -> CharacterData -> CharacterData
getDOMedianCharInfoUnion charInfo = getDOMedianUnion (weight charInfo) (costMatrix charInfo) (slimTCM charInfo) (wideTCM charInfo) (hugeTCM charInfo) (charType charInfo)

-- | getDOMedianUnion calls appropriate pairwise DO to create sequence median after some type wrangling
-- works on union states
-- filters out gaps (0/1) values before DO (>1)
getDOMedianUnion
  :: Double
  -> S.Matrix Int
  -> TCMD.DenseTransitionCostMatrix
  -> MR.MetricRepresentation WideState
  -> MR.MetricRepresentation BV.BitVector
  -> CharType
  -> CharacterData
  -> CharacterData
  -> CharacterData
getDOMedianUnion thisWeight thisMatrix thisSlimTCM thisWideTCM thisHugeTCM thisType leftChar rightChar
  | null thisMatrix = error "Null cost matrix in getDOMedian"
  | thisType `elem` [SlimSeq,   NucSeq] = newSlimCharacterData
  | thisType `elem` [WideSeq, AminoSeq] = newWideCharacterData
  | thisType == HugeSeq           = newHugeCharacterData
  | otherwise = error $ fold ["Unrecognised character type '", show thisType, "'in a DYNAMIC character branch" ]
  where
    blankCharacterData = emptyCharacter

    newSlimCharacterData =
        let newCost     = thisWeight * fromIntegral cost
            subtreeCost = sum [ newCost, globalCost leftChar, globalCost rightChar]
            slimIAUnionNoGapsLeft = extractMediansSingle $ slimIAUnion leftChar
            slimIAUnionNoGapsRight = extractMediansSingle $ slimIAUnion rightChar
            (cost, r)   = slimPairwiseDO
                thisSlimTCM (makeDynamicCharacterFromSingleVector slimIAUnionNoGapsLeft) (makeDynamicCharacterFromSingleVector slimIAUnionNoGapsRight)
        in
        --trace ("GDOMU:" <> show (cost, extractMedians r, slimIAUnionNoGapsLeft, slimIAUnionNoGapsRight)) $
        blankCharacterData
              { slimIAUnion   = extractMedians r
              , localCostVect = V.singleton $ fromIntegral cost
              , localCost     = newCost
              , globalCost    = subtreeCost
              }

    newWideCharacterData =
        let newCost     = thisWeight * fromIntegral cost
            coefficient = MR.minInDelCost thisWideTCM
            subtreeCost = sum [ newCost, globalCost leftChar, globalCost rightChar]
            wideIAUnionNoGapsLeft = extractMediansSingle $ wideIAUnion leftChar
            wideIAUnionNoGapsRight = extractMediansSingle $ wideIAUnion rightChar
            (cost, r)   = widePairwiseDO
                coefficient
                (MR.retreivePairwiseTCM thisWideTCM)
                (makeDynamicCharacterFromSingleVector wideIAUnionNoGapsLeft) (makeDynamicCharacterFromSingleVector wideIAUnionNoGapsRight)
        in  blankCharacterData
              { wideIAUnion    = extractMedians r
              , localCostVect = V.singleton $ fromIntegral cost
              , localCost     = newCost
              , globalCost    = subtreeCost
              }

    newHugeCharacterData =
        let newCost     = thisWeight * fromIntegral cost
            coefficient = MR.minInDelCost thisHugeTCM
            subtreeCost = newCost + globalCost leftChar + globalCost rightChar
            hugeIAUnionNoGapsLeft = extractMediansSingle $ hugeIAUnion leftChar
            hugeIAUnionNoGapsRight = extractMediansSingle $ hugeIAUnion rightChar
            (cost, r)   = hugePairwiseDO
                coefficient
                (MR.retreivePairwiseTCM thisHugeTCM)
                (makeDynamicCharacterFromSingleVector hugeIAUnionNoGapsLeft) (makeDynamicCharacterFromSingleVector hugeIAUnionNoGapsRight)
        in blankCharacterData
              { hugeIAUnion = extractMedians r
              , localCostVect = V.singleton $ fromIntegral cost
              , localCost  = newCost
              , globalCost = subtreeCost
              }


-- | getDOMedianCharInfo  is a wrapper around getDOMedian with CharInfo-based interface
getDOMedianCharInfo :: CharInfo -> CharacterData -> CharacterData -> CharacterData
getDOMedianCharInfo charInfo = getDOMedian (weight charInfo) (costMatrix charInfo) (slimTCM charInfo) (wideTCM charInfo) (hugeTCM charInfo) (charType charInfo)

-- | getDOMedian calls appropriate pairwise DO to create sequence median after some type wrangling
-- works on preliminary states
getDOMedian
  :: Double
  -> S.Matrix Int
  -> TCMD.DenseTransitionCostMatrix
  -> MR.MetricRepresentation WideState
  -> MR.MetricRepresentation BV.BitVector
  -> CharType
  -> CharacterData
  -> CharacterData
  -> CharacterData
getDOMedian thisWeight thisMatrix thisSlimTCM thisWideTCM thisHugeTCM thisType leftChar rightChar
  | null thisMatrix = error "Null cost matrix in getDOMedian"
  | thisType `elem` [SlimSeq,   NucSeq] = newSlimCharacterData
  | thisType `elem` [WideSeq, AminoSeq] = newWideCharacterData
  | thisType == HugeSeq           = newHugeCharacterData
  | otherwise = error $ fold ["Unrecognised character type '", show thisType, "'in a DYNAMIC character branch" ]
  where
    blankCharacterData = emptyCharacter

    newSlimCharacterData =
        let newCost     = thisWeight * fromIntegral cost
            subtreeCost = sum [ newCost, globalCost leftChar, globalCost rightChar]
            (cost, r)   = slimPairwiseDO
                thisSlimTCM (slimGapped leftChar) (slimGapped rightChar)
        in  blankCharacterData
              { slimPrelim    = extractMedians r
              , slimGapped    = r
              , localCostVect = V.singleton $ fromIntegral cost
              , localCost     = newCost
              , globalCost    = subtreeCost
              }

    newWideCharacterData =
        let newCost     = thisWeight * fromIntegral cost
            coefficient = MR.minInDelCost thisWideTCM
            subtreeCost = sum [ newCost, globalCost leftChar, globalCost rightChar]
            (cost, r)   = widePairwiseDO
                coefficient
                (MR.retreivePairwiseTCM thisWideTCM)
                (wideGapped leftChar) (wideGapped rightChar)
        in  blankCharacterData
              { widePrelim    = extractMedians r
              , wideGapped    = r
              , localCostVect = V.singleton $ fromIntegral cost
              , localCost     = newCost
              , globalCost    = subtreeCost
              }

    newHugeCharacterData =
        let newCost     = thisWeight * fromIntegral cost
            coefficient = MR.minInDelCost thisHugeTCM
            subtreeCost = newCost + globalCost leftChar + globalCost rightChar
            (cost, r)   = hugePairwiseDO
                coefficient
                (MR.retreivePairwiseTCM thisHugeTCM)
                (hugeGapped leftChar) (hugeGapped rightChar)
        in blankCharacterData
              { hugePrelim = extractMedians r
              , hugeGapped = r
              , localCostVect = V.singleton $ fromIntegral cost
              , localCost  = newCost
              , globalCost = subtreeCost
              }

-- | getPrealignedUnion calls appropriate pairwise function to create sequence union of final states
-- for prealigned states
getPrealignedUnion :: CharType
                   -> CharacterData
                   -> CharacterData
                   -> CharacterData
getPrealignedUnion thisType leftChar rightChar =
    let blankCharacterData = emptyCharacter
    in
    if thisType == AlignedSlim then
        let finalUnion = GV.zipWith (.|.) (alignedSlimFinal leftChar) (alignedSlimFinal rightChar)
            prelimState = (finalUnion, finalUnion, finalUnion)
        in
        blankCharacterData { alignedSlimPrelim  = prelimState
                           , alignedSlimFinal = finalUnion
                           }
    else if thisType == AlignedWide then
        let finalUnion = GV.zipWith (.|.) (alignedWideFinal leftChar) (alignedWideFinal rightChar)
            prelimState = (finalUnion, finalUnion, finalUnion)
        in
        blankCharacterData { alignedWidePrelim  = prelimState
                           , alignedWideFinal = finalUnion
                           }

    else if thisType == AlignedHuge then
        let finalUnion = GV.zipWith (.|.) (alignedHugeFinal leftChar) (alignedHugeFinal rightChar)
            prelimState = (finalUnion, finalUnion, finalUnion)
        in
        blankCharacterData { alignedHugePrelim  = prelimState
                           , alignedHugeFinal = finalUnion
                           }

    else error ("Unrecognised character type '" <> show thisType <> "' in getPrealignedUnion")


-- | getDynamicUnion calls appropriate pairwise function to create sequence median after some type wrangling
-- if using IA--takes IAFInal for each node, creates union of IAFinals states
-- if DO then calculated DO medians and takes union of left and right states
-- gaps need to be fitered if DO used later (as in Wagner), or as in SPR/TBR rearragement
-- sets final and IA states for Swap delta heuristics
getDynamicUnion :: Bool
                -> Bool
                -> CharType
                -> CharacterData
                -> CharacterData
                -> TCMD.DenseTransitionCostMatrix
                -> MR.MetricRepresentation WideState
                -> MR.MetricRepresentation BV.BitVector
                -> CharacterData
getDynamicUnion useIA filterGaps thisType leftChar rightChar thisSlimTCM thisWideTCM thisHugeTCM
  | thisType `elem` [SlimSeq,   NucSeq] = newSlimCharacterData
  | thisType `elem` [WideSeq, AminoSeq] = newWideCharacterData
  | thisType == HugeSeq           = newHugeCharacterData
  | otherwise = error $ fold ["Unrecognised character type '", show thisType, "'in a DYNAMIC character branch" ]
  where
    blankCharacterData = emptyCharacter

    newSlimCharacterData =
        let r  = if useIA then GV.zipWith (.|.) (slimIAFinal leftChar) (slimIAFinal rightChar)
                 else
                    let (_, (lG, _, rG)) = slimPairwiseDO thisSlimTCM (makeDynamicCharacterFromSingleVector $ slimFinal leftChar) (makeDynamicCharacterFromSingleVector $ slimFinal rightChar)
                    in
                    GV.zipWith (.|.) lG rG

            r' = if filterGaps then extractMediansSingle r
                 else r

        in  blankCharacterData { slimPrelim = r'
                               , slimGapped = (r', r',r')
                               , slimFinal = r'
                               , slimIAPrelim = (r, r, r)
                               , slimIAFinal = r
                               }

    newWideCharacterData =
        let r  = if useIA then GV.zipWith (.|.) (wideIAFinal leftChar) (wideIAFinal rightChar)
                 else
                    let coefficient = MR.minInDelCost thisWideTCM
                        (_, (lG, _, rG)) = widePairwiseDO coefficient (MR.retreivePairwiseTCM thisWideTCM) (makeDynamicCharacterFromSingleVector $ wideFinal leftChar) (makeDynamicCharacterFromSingleVector $ wideFinal rightChar)
                    in
                    GV.zipWith (.|.) lG rG
         -- r   = GV.zipWith (.|.) (wideIAFinal leftChar) (wideIAFinal rightChar)

            r' = if filterGaps then extractMediansSingle  r
                 else r

        in  blankCharacterData { widePrelim = r'
                               , wideGapped = (r', r',r')
                               , wideFinal = r'
                               , wideIAPrelim = (r, r, r)
                               , wideIAFinal = r
                               }

    newHugeCharacterData =
        let r  = if useIA then GV.zipWith (.|.) (hugeIAFinal leftChar) (hugeIAFinal rightChar)
                 else
                    let coefficient = MR.minInDelCost thisHugeTCM
                        (_, (lG, _, rG)) = hugePairwiseDO coefficient (MR.retreivePairwiseTCM thisHugeTCM) (makeDynamicCharacterFromSingleVector $ hugeFinal leftChar) (makeDynamicCharacterFromSingleVector $ hugeFinal rightChar)
                    in
                    GV.zipWith (.|.) lG rG

         -- r   = GV.zipWith (.|.) (hugeIAFinal leftChar) (hugeIAFinal rightChar)

            r' = if filterGaps then extractMediansSingle  r
                 else r

        in  blankCharacterData { hugePrelim = r'
                               , hugeGapped = (r', r',r')
                               , hugeFinal = r'
                               , hugeIAPrelim = (r, r, r)
                               , hugeIAFinal = r
                               }

-- | union2 takes the vectors of characters and applies union2Single to each character
-- used for edge states in buikd and rearrangement
union2 ::  Bool -> Bool -> V.Vector CharacterData -> V.Vector CharacterData -> V.Vector CharInfo -> V.Vector CharacterData
union2 useIA filterGaps = V.zipWith3 (union2Single useIA filterGaps)

-- | union2Single takes character data and returns union character data
-- union2Single assumes that the character vectors in the various states are the same length
-- that is--all leaves (hence other vertices later) have the same number of each type of character
-- used IAFinal states for dynamic characters
-- used in heurstic graph build and rearrangement
union2Single :: Bool -> Bool -> CharacterData -> CharacterData -> CharInfo -> CharacterData
union2Single useIA filterGaps firstVertChar secondVertChar inCharInfo =
    let thisType    = charType inCharInfo
        thisMatrix  = costMatrix inCharInfo
        thisActive  = activity inCharInfo
        thisSlimTCM = slimTCM inCharInfo
        thisWideTCM = wideTCM inCharInfo
        thisHugeTCM = hugeTCM inCharInfo
    in
    if not thisActive then firstVertChar
    else if thisType == Add then
        intervalUnion firstVertChar secondVertChar

    else if thisType == NonAdd then
        localUnion firstVertChar secondVertChar

    else if thisType `elem` packedNonAddTypes then
        BP.unionPacked firstVertChar secondVertChar

    else if thisType == Matrix then
        unionMatrix thisMatrix firstVertChar secondVertChar

    else if thisType `elem` prealignedCharacterTypes then
        getPrealignedUnion thisType firstVertChar secondVertChar

    else if thisType `elem` nonExactCharacterTypes then
        getDynamicUnion useIA filterGaps thisType firstVertChar secondVertChar thisSlimTCM thisWideTCM thisHugeTCM

    else error ("Character type " <> show thisType <> " unrecognized/not implemented")

-- | makeEdgeData takes and edge and makes the VertData for the edge from the union of the two vertices
-- using IA assignments not so great for search deltas
makeEdgeData :: Bool -> Bool -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.LEdge b -> VertexBlockData
makeEdgeData useIA filterGaps inGraph charInfoVV (eNode, vNode, _) =
   let eNodeVertData = vertData $ fromJust $ LG.lab inGraph eNode
       vNodeVertData = vertData $ fromJust $ LG.lab inGraph vNode
   in
   createEdgeUnionOverBlocks useIA filterGaps eNodeVertData vNodeVertData charInfoVV []

-- | createEdgeUnionOverBlocks creates the union of the final states characters on an edge
-- The function takes data in blocks and block vector of char info and
-- extracts the triple for each block and creates new block data
-- this is used for delta's in edge invastion in Wagner and SPR/TBR
-- filter gaps for using with DO (flterGaps = True) or IA (filterGaps = False)
createEdgeUnionOverBlocks :: Bool
                          -> Bool
                          -> VertexBlockData
                          -> VertexBlockData
                          -> V.Vector (V.Vector CharInfo)
                          -> [V.Vector CharacterData]
                          -> V.Vector (V.Vector CharacterData)
createEdgeUnionOverBlocks useIA filterGaps leftBlockData rightBlockData blockCharInfoVect curBlockData =
    if V.null leftBlockData then
        --trace ("Blocks: " <> (show $ length curBlockData) <> " Chars  B0: " <> (show $ V.map snd $ head curBlockData))
        V.fromList $ reverse curBlockData
    else
        let leftBlockLength = length $ V.head leftBlockData
            rightBlockLength =  length $ V.head rightBlockData
            -- firstBlock = V.zip3 (V.head leftBlockData) (V.head rightBlockData) (V.head blockCharInfoVect)

            -- missing data cases first or zip defaults to zero length
            firstBlockMedian
              | (leftBlockLength == 0) = V.head rightBlockData
              | (rightBlockLength == 0) = V.head leftBlockData
              | otherwise = union2 useIA filterGaps (V.head leftBlockData) (V.head rightBlockData) (V.head blockCharInfoVect)
        in
        createEdgeUnionOverBlocks useIA filterGaps (V.tail leftBlockData) (V.tail rightBlockData) (V.tail blockCharInfoVect) (firstBlockMedian : curBlockData)

-- | getPreAligned2Median takes prealigned character types (AlignedSlim, AlignedWide, AlignedHuge) and returns 2-median and cost
-- uses IA-type functions for slim/wide/huge
getPreAligned2Median :: CharInfo -> CharacterData -> CharacterData -> CharacterData -> CharacterData
getPreAligned2Median charInfo nodeChar leftChar rightChar =
    let setCost cVal r = r
            { localCost = weight charInfo * fromIntegral cVal
            , globalCost = sum [ weight charInfo * fromIntegral cVal, globalCost leftChar, globalCost rightChar]
            }

        setSlimPrelim v r = r { alignedSlimPrelim = v }
        setWidePrelim v r = r { alignedWidePrelim = v }
        setHugePrelim v r = r { alignedHugePrelim = v }

        getCharL f = extractMediansGapped $ f leftChar
        getCharR f = extractMediansGapped $ f rightChar

        (setter, cost) = case charType charInfo of
            AlignedSlim ->
                let cL = getCharL alignedSlimPrelim
                    cR = getCharR alignedSlimPrelim
                    (cM, score) = get2WaySlim (slimTCM charInfo) cL cR
                in  (setSlimPrelim (cL, cM, cR), score)

            AlignedWide ->
                let cL = getCharL alignedWidePrelim
                    cR = getCharR alignedWidePrelim
                    (cM, score) = get2WayWideHuge (wideTCM charInfo) cL cR
                in  (setWidePrelim (cL, cM, cR), score)

            AlignedHuge ->
                let cL = getCharL alignedHugePrelim
                    cR = getCharR alignedHugePrelim
                    (cM, score) = get2WayWideHuge (hugeTCM charInfo) cL cR
                in  (setHugePrelim (cL, cM, cR), score)

            other -> error $ "Unrecognized character type: " <> show other

    in  setter $ setCost cost nodeChar
{-
    let characterType = charType charInfo
    in
    if characterType == AlignedSlim then
        let (prelimChar, cost) = get2WaySlim (slimTCM charInfo) (extractMediansGapped $ alignedSlimPrelim leftChar) (extractMediansGapped $ alignedSlimPrelim rightChar)
        in
        -- trace ("GPA2M: " <> (show $ GV.length prelimChar))
        nodeChar { alignedSlimPrelim = (extractMediansGapped $ alignedSlimPrelim leftChar, prelimChar,  extractMediansGapped $ alignedSlimPrelim rightChar)
                 , localCost = weight charInfo * fromIntegral cost
                 , globalCost = sum [ weight charInfo * fromIntegral cost, globalCost leftChar, globalCost rightChar]
                 }

    else if characterType == AlignedWide then
        let (prelimChar, cost) = get2WayWideHuge (wideTCM charInfo) (extractMediansGapped $ alignedWidePrelim leftChar) (extractMediansGapped $ alignedWidePrelim rightChar)
        in
        nodeChar { alignedWidePrelim = (extractMediansGapped $ alignedWidePrelim leftChar, prelimChar,  extractMediansGapped $ alignedWidePrelim rightChar)
                 , localCost = weight charInfo * fromIntegral cost
                 , globalCost = sum [ weight charInfo * fromIntegral cost, globalCost leftChar, globalCost rightChar]
                 }

    else if characterType == AlignedHuge then
        let (prelimChar, cost) = get2WayWideHuge (hugeTCM charInfo) (extractMediansGapped $ alignedHugePrelim leftChar) (extractMediansGapped $ alignedHugePrelim rightChar)
        in
        nodeChar { alignedHugePrelim = (extractMediansGapped $ alignedHugePrelim leftChar, prelimChar,  extractMediansGapped $ alignedHugePrelim rightChar)
                 , localCost = weight charInfo * fromIntegral cost
                 , globalCost = sum [ weight charInfo * fromIntegral cost, globalCost leftChar, globalCost rightChar]
                 }

    else error ("Unrecognized character type " <> show characterType)
-}


-- | getPreAligned2MedianUnionFields takes prealigned character types (AlignedSlim, AlignedWide, AlignedHuge) and returns 2-median and cost
-- uses IA-type functions for slim/wide/huge-- based on union fields
getPreAligned2MedianUnionFields :: CharInfo -> CharacterData -> CharacterData -> CharacterData -> CharacterData
getPreAligned2MedianUnionFields charInfo nodeChar leftChar rightChar =
    let characterType = charType charInfo
    in
    if characterType == AlignedSlim then
        let (prelimChar, cost) = get2WaySlim (slimTCM charInfo) (alignedSlimUnion leftChar) (alignedSlimUnion rightChar)
        in
        -- trace ("GPA2M-slim: " <> (show (GV.length prelimChar, GV.length $ alignedSlimUnion leftChar, GV.length $ alignedSlimUnion rightChar)))
        nodeChar { alignedSlimUnion = prelimChar
                 , localCost = weight charInfo * fromIntegral cost
                 , globalCost = sum [ weight charInfo * fromIntegral cost, globalCost leftChar, globalCost rightChar]
                 }

    else if characterType == AlignedWide then
        let (prelimChar, cost) = get2WayWideHuge (wideTCM charInfo) (alignedWideUnion leftChar) (alignedWideUnion rightChar)
        in
        --trace ("GPA2M-wide: " <> (show $ GV.length prelimChar))
        nodeChar { alignedWideUnion = prelimChar
                 , localCost = weight charInfo * fromIntegral cost
                 , globalCost = sum [ weight charInfo * fromIntegral cost, globalCost leftChar, globalCost rightChar]
                 }

    else if characterType == AlignedHuge then
        let (prelimChar, cost) = get2WayWideHuge (hugeTCM charInfo) (alignedHugeUnion leftChar) (alignedHugeUnion rightChar)
        in
        --trace ("GPA2M-huge: " <> (show $ GV.length prelimChar))
        nodeChar { alignedHugeUnion = prelimChar
                 , localCost = weight charInfo * fromIntegral cost
                 , globalCost = sum [ weight charInfo * fromIntegral cost, globalCost leftChar, globalCost rightChar]
                 }

    else error ("Unrecognized character type " <> show characterType)



-- | makeIAUnionPrelimLeaf makes union and sets for leaf characters--leaves alignment fields unchanged
makeIAUnionPrelimLeaf :: CharInfo -> CharacterData -> CharacterData
makeIAUnionPrelimLeaf charInfo nodeChar  =
    --traceNoLF ("In makeIAUnionPrelimLeaf") (
    let characterType = charType charInfo
    in
    if characterType == NonAdd then
        let prelimState = snd3 $ stateBVPrelim nodeChar
        in
        nodeChar {stateBVUnion = prelimState}

    else if characterType == Add then
        let prelimState = snd3 $ rangePrelim nodeChar
        in
        nodeChar {rangeUnion = prelimState}

    else if characterType == Matrix then
        let prelimState = matrixStatesPrelim nodeChar
        in
        nodeChar {matrixStatesUnion = prelimState}

    -- the checking for null in prelimstate (slim/wide/huge sequences) comes form case of splitting graph
    -- and reoptimizing in swap and fuse but single leaf is split off

    else if characterType `elem` [SlimSeq, NucSeq] then
        let prelimState = extractMediansGapped $ slimAlignment nodeChar
            unionState = unionBVOrMissing prelimState (length $ alphabet charInfo) (slimGapped nodeChar)
        in
        nodeChar { slimIAPrelim = makeDynamicCharacterFromSingleVector unionState -- slimAlignment nodeChar
                 , slimIAFinal = unionState -- prelimState
                 , slimIAUnion = unionState
                 }

    else if characterType `elem` [WideSeq, AminoSeq] then
        let prelimState = extractMediansGapped $ wideAlignment nodeChar
            unionState = unionBVOrMissing prelimState (length $ alphabet charInfo) (wideGapped nodeChar)
        in
        nodeChar { wideIAPrelim = makeDynamicCharacterFromSingleVector unionState -- ideAlignment nodeChar
                 , wideIAFinal = unionState -- prelimState
                 , wideIAUnion = unionState
                 }

    else if characterType == HugeSeq then
        let prelimState = extractMediansGapped $ hugeAlignment nodeChar
            unionState = unionBVOrMissing prelimState (length $ alphabet charInfo) (hugeGapped nodeChar)
        in
        nodeChar { hugeIAPrelim = makeDynamicCharacterFromSingleVector unionState -- hugeAlignment nodeChar
                 , hugeIAFinal = unionState -- prelimState
                 , hugeIAUnion = unionState
                 }

    else if characterType == AlignedSlim then
        let prelimState = snd3 $ alignedSlimPrelim nodeChar
        in
        nodeChar {alignedSlimUnion = prelimState}

    else if characterType == AlignedWide then
        let prelimState = snd3 $ alignedWidePrelim nodeChar
        in
        nodeChar {alignedWideUnion = prelimState}

    else if characterType == AlignedHuge then
        let prelimState = snd3 $ alignedHugePrelim nodeChar
        in
        nodeChar {alignedHugeUnion = prelimState}

    else if characterType `elem` [Packed2, Packed4, Packed5, Packed8, Packed64] then
        let prelimState = snd3 $ packedNonAddPrelim nodeChar
        in
        nodeChar {packedNonAddUnion = prelimState}

    else error ("Unrecognized character type " <> show characterType)
    --)

-- | unionBVOrMissing returns leaf algnment state, if all '-' converts to missing characters
--  so BV unions and IAFinal get a zero cost as opposed to checking against all '-' seqquence
unionBVOrMissing :: (FiniteBits e, GV.Vector v e) => v e -> Int -> (v e, v e, v e) -> v e
unionBVOrMissing prelimState alphSize nodeGapped =
    if not $ GV.null prelimState then
        if GV.null $ extractMediansSingle prelimState then
        --trace ("MIAUPL: " <> (show $ convertIfAllGapsToAllBitsOn (length $ alphabet charInfo) prelimState))
            convertAllGapsToAllBitsOn alphSize prelimState
        else prelimState
    else extractMedians nodeGapped

-- | convertAllGapsToAllBitsOn takes a single fields of a dynamic character and
-- converts replaces
-- 'gaps' with all bits on--in essence '?' or missing element
convertAllGapsToAllBitsOn :: (FiniteBits e, GV.Vector v e) => Int -> v e -> v e
convertAllGapsToAllBitsOn alphSize inVect =
    if GV.null inVect then inVect
    else
        let numElements = GV.length inVect
            onBitList =  fmap (setBit (inVect GV.! 0)) [0..alphSize - 1]
            onBits = L.foldl1' (.|.) onBitList
        in
        GV.replicate numElements onBits

-- | allMissingBits test if all bits in alphabet size are ON
allMissingBits :: (FiniteBits e, GV.Vector v e) => Int -> v e -> Bool
allMissingBits alphSize inVect =
    not (GV.null inVect) && (let onBitList =  fmap (setBit (inVect GV.! 0)) [0..alphSize - 1]
                                 missingBits = L.foldl1' (.|.) onBitList
                                 offBits = GV.filter (/= missingBits) inVect
                             in
                             GV.null offBits)
{-
-- | getNonExactUnionFields takes two non-exact characters and union field assignment
-- based on character type and nodeChar
-- modifies both IA and union fields
-- union fields are unions of parent states (aligned, or IA, or static)
getNonExactUnionFields :: CharInfo -> CharacterData -> CharacterData -> CharacterData -> CharacterData
getNonExactUnionFields charInfo nodeChar leftChar rightChar =
    let characterType = charType charInfo
    in
    if characterType `elem` [SlimSeq, NucSeq] then
        let (prelimChar, cost) = get2WaySlim (slimTCM charInfo) (slimIAUnion leftChar) (slimIAUnion rightChar)
        in
        nodeChar { slimIAUnion = prelimChar
                   , localCost = (weight charInfo) * (fromIntegral cost)
                   , globalCost = sum [ (weight charInfo) * (fromIntegral cost), globalCost leftChar, globalCost rightChar]
                    }
    else if characterType `elem` [WideSeq, AminoSeq] then
        let (prelimChar, minCost) = get2WayWideHuge (wideTCM charInfo) (wideIAUnion leftChar) (wideIAUnion rightChar)
        in
        nodeChar { wideIAUnion = prelimChar
                   , localCost = (weight charInfo) * (fromIntegral minCost)
                   , globalCost = sum [ (weight charInfo) * (fromIntegral minCost), globalCost leftChar, globalCost rightChar]
                   }
    else if characterType == HugeSeq then
        let (prelimChar, minCost) = get2WayWideHuge (hugeTCM charInfo) (hugeIAUnion leftChar) (hugeIAUnion rightChar)
        in
        nodeChar {hugeIAUnion = prelimChar
                , localCost = (weight charInfo) * (fromIntegral minCost)
                , globalCost = sum [ (weight charInfo) * (fromIntegral minCost), globalCost leftChar, globalCost rightChar]
                }
    else error ("Unrecognized character type " <> show characterType)
-}


-- | makeIAPrelimCharacter takes two characters and performs 2-way assignment
-- based on character type and nodeChar
-- modifies both IA and union fields
-- union fields are unions of parent states (aligned, or IA, or static)
makeIAPrelimCharacter :: CharInfo -> CharacterData -> CharacterData -> CharacterData -> CharacterData
makeIAPrelimCharacter charInfo nodeChar leftChar rightChar =
     --traceNoLF ("In makeIAPrelimCharacter") (
     let characterType = charType charInfo
     in
     if characterType == NonAdd then
        let leftState = stateBVUnion leftChar
            rightState = stateBVUnion rightChar
            unionState = V.zipWith (.|.) leftState rightState
        in
        nodeChar {stateBVUnion = unionState}

    else if characterType == Add then
        let prelimState = V.zipWith4 getUnionRange (V.map fst $ rangeUnion leftChar) (V.map snd $ rangeUnion leftChar) (V.map fst $ rangeUnion rightChar) (V.map snd $ rangeUnion rightChar)
        in
        nodeChar {rangeUnion = prelimState}

    else if characterType == Matrix then
        let numStates = length $ costMatrix charInfo
            newMatrixVector = getUnionVector (costMatrix charInfo) numStates <$> V.zip (matrixStatesUnion leftChar) (matrixStatesUnion rightChar)
        in
        nodeChar {matrixStatesUnion = newMatrixVector}
    else if characterType == AlignedSlim then
        let prelimState = GV.zipWith (.|.) (alignedSlimUnion leftChar) (alignedSlimUnion rightChar)
        in
        nodeChar {alignedSlimUnion = prelimState}

    else if characterType == AlignedWide then
        let prelimState = GV.zipWith (.|.) (alignedWideUnion leftChar) (alignedWideUnion rightChar)
        in
        nodeChar {alignedWideUnion = prelimState}

    else if characterType == AlignedHuge then
        let prelimState = GV.zipWith (.|.) (alignedHugeUnion leftChar) (alignedHugeUnion rightChar)
        in
        nodeChar {alignedHugeUnion = prelimState}

    else if characterType `elem` [Packed2, Packed4, Packed5, Packed8, Packed64] then
        let prelimState = GV.zipWith (.|.) (packedNonAddUnion leftChar) (packedNonAddUnion leftChar)
        in
        nodeChar {packedNonAddUnion = prelimState}

   else if characterType `elem` [SlimSeq, NucSeq] then
        let (prelimChar, cost) = get2WaySlim (slimTCM charInfo) (extractMediansGapped $ slimIAPrelim leftChar) (extractMediansGapped $ slimIAPrelim rightChar)
        in
        -- trace ("MPC: " <> (show prelimChar) <> "\nleft: " <> (show $ extractMediansGapped $ slimIAPrelim leftChar) <> "\nright: " <> (show $ extractMediansGapped $ slimIAPrelim rightChar))
        -- trace ("MIAUP-C: " <> (show $ GV.length $ GV.zipWith (.|.) (slimIAUnion leftChar) (slimIAUnion rightChar)))

        -- the check for all missing basically creates an intersection result so all missing isn't porpagated post-order
        nodeChar {slimIAPrelim = (extractMediansGapped $ slimIAPrelim leftChar
                , prelimChar,  extractMediansGapped $ slimIAPrelim rightChar)
                , slimIAUnion = orBVOrMissingIntersection (length $ alphabet charInfo) (slimIAUnion leftChar) (slimIAUnion rightChar)
                , localCost = if GV.null (slimIAUnion leftChar) || GV.null (slimIAUnion rightChar) then 0
                              else weight charInfo * fromIntegral cost
                , globalCost = sum [ weight charInfo * fromIntegral cost, globalCost leftChar, globalCost rightChar]
                }
    else if characterType `elem` [WideSeq, AminoSeq] then
        let (prelimChar, minCost)  = get2WayWideHuge (wideTCM charInfo) (extractMediansGapped $ wideIAPrelim leftChar) (extractMediansGapped $ wideIAPrelim rightChar)
        in
        nodeChar {wideIAPrelim = (extractMediansGapped $ wideIAPrelim leftChar
                , prelimChar, extractMediansGapped $ wideIAPrelim rightChar)
                , wideIAUnion = orBVOrMissingIntersection (length $ alphabet charInfo) (wideIAUnion leftChar) (wideIAUnion rightChar)
                , localCost = weight charInfo * fromIntegral minCost
                , globalCost = sum [ weight charInfo * fromIntegral minCost, globalCost leftChar, globalCost rightChar]
                }
    else if characterType == HugeSeq then
        let (prelimChar, minCost)  = get2WayWideHuge (hugeTCM charInfo) (extractMediansGapped $ hugeIAPrelim leftChar) (extractMediansGapped $ hugeIAPrelim rightChar)
        in
        nodeChar {hugeIAPrelim = (extractMediansGapped $ hugeIAPrelim leftChar
                , prelimChar, extractMediansGapped $ hugeIAPrelim rightChar)
                , hugeIAUnion = orBVOrMissingIntersection (length $ alphabet charInfo) (hugeIAUnion leftChar) (hugeIAUnion rightChar)
                , localCost = weight charInfo * fromIntegral minCost
                , globalCost = sum [ weight charInfo * fromIntegral minCost, globalCost leftChar, globalCost rightChar]
                }
    else nodeChar --error ("Unrecognized character type " <> show characterType)
    -- )

-- | orBVOrMissingIntersection takes two uninBV seqs and returns union or intersectino if one/both missing
orBVOrMissingIntersection :: (FiniteBits e, GV.Vector v e) => Int -> v e -> v e -> v e
orBVOrMissingIntersection alphSize unionIALeft unionIARight
  | GV.null unionIALeft || allMissingBits alphSize unionIALeft = unionIARight
  | GV.null unionIARight  || allMissingBits alphSize unionIARight = unionIALeft
  | otherwise = GV.zipWith (.|.) unionIALeft unionIARight

-- | makeIAFinalCharacterStaticIA takes two characters and performs 2-way assignment
-- based on character type and nodeChar--only IA fields are modified
-- this pulls from current node for left and right states
-- skips other than unaligned sequence characters
makeIAFinalCharacter :: AssignmentMethod -> CharInfo -> CharacterData -> CharacterData-> CharacterData
makeIAFinalCharacter finalMethod charInfo nodeChar parentChar  =
     let characterType = charType charInfo
     in
     if characterType `elem` [SlimSeq, NucSeq] then
        let finalIAChar = getFinal3WaySlim (slimTCM charInfo) (slimIAFinal parentChar) (extractMediansLeftGapped $ slimIAPrelim nodeChar) (extractMediansRightGapped $ slimIAPrelim nodeChar)
            finalChar =  if finalMethod == ImpliedAlignment then extractMedians $ makeDynamicCharacterFromSingleVector finalIAChar
                         else slimFinal nodeChar
        in
        nodeChar { slimIAFinal = finalIAChar
                 , slimFinal = finalChar
                 }
     else if characterType `elem` [WideSeq, AminoSeq] then
        let finalIAChar = getFinal3WayWideHuge (wideTCM charInfo) (wideIAFinal parentChar) (extractMediansLeftGapped $ wideIAPrelim nodeChar) (extractMediansRightGapped $ wideIAPrelim nodeChar)
            finalChar = if finalMethod == ImpliedAlignment then extractMedians $ makeDynamicCharacterFromSingleVector finalIAChar
                        else wideFinal nodeChar
        in
        nodeChar { wideIAFinal = finalIAChar
                 , wideFinal = finalChar
                 }
     else if characterType == HugeSeq then
        let finalIAChar = getFinal3WayWideHuge (hugeTCM charInfo) (hugeIAFinal parentChar) (extractMediansLeftGapped $ hugeIAPrelim nodeChar) (extractMediansRightGapped $ hugeIAPrelim nodeChar)
            finalChar = if finalMethod == ImpliedAlignment then extractMedians $ makeDynamicCharacterFromSingleVector finalIAChar
                        else hugeFinal nodeChar
        in
        nodeChar { hugeIAFinal = finalIAChar
                 , hugeFinal = finalChar
                 }
     else nodeChar -- error ("Unrecognized character type " <> show characterType)

-- | get2WaySlim takes two slim vectors an produces a preliminary median
get2WayGeneric :: (FiniteBits e, GV.Vector v e) => (e -> e -> (e, Word)) -> v e -> v e -> (v e, Word)
get2WayGeneric tcm descendantLeftPrelim descendantRightPrelim =
   let -- this should not be needed some problems at times with IA
       -- len   = GV.length descendantLeftPrelim
       len   = min (GV.length descendantLeftPrelim) (GV.length descendantRightPrelim)
       vt    = V.generate len $ \i -> tcm (descendantLeftPrelim GV.! i) (descendantRightPrelim GV.! i) -- :: V.Vector (SlimState, Word)
       gen v = let med i = fst $ v V.! i in GV.generate len med
       add   = V.foldl' (\x e -> x + snd e) 0
   in  (,) <$> gen <*> add $ vt


-- | get2WaySlim takes two slim vectors an produces a preliminary median
get2WaySlim :: TCMD.DenseTransitionCostMatrix -> SV.Vector SlimState -> SV.Vector SlimState -> (SV.Vector SlimState, Word)
get2WaySlim lSlimTCM = get2WayGeneric (TCMD.lookupPairwise lSlimTCM)


-- | get2WayWideHuge like get2WaySlim but for wide and huge characters
get2WayWideHuge :: (FiniteBits a, GV.Vector v a) => MR.MetricRepresentation a -> v a -> v a -> (v a, Word)
get2WayWideHuge whTCM = get2WayGeneric (MR.retreivePairwiseTCM whTCM)

-- | getFinal3Way takes parent final assignment (including indel characters) and descendent
-- preliminary gapped assingment from postorder and creates a gapped final assignment based on
-- minimum cost median for the three inputs.  THis is done to preserve the ImpliedAlignment
-- information to create a final assingment with out an additional DO call to keep the
-- creation linear in sequence length.  Since gaps remain--they must be filtered when output or
-- used as true final sequence assignments using M.createUngappedMedianSequence
getFinal3WaySlim :: TCMD.DenseTransitionCostMatrix -> SV.Vector SlimState -> SV.Vector SlimState -> SV.Vector SlimState -> SV.Vector SlimState
getFinal3WaySlim lSlimTCM parentFinal descendantLeftPrelim descendantRightPrelim =
   let newFinal = removeGapAndNil $ SV.zipWith3 (local3WaySlim lSlimTCM) parentFinal descendantLeftPrelim descendantRightPrelim
   in
   newFinal

-- | getFinal3WayWideHuge like getFinal3WaySlim but for wide and huge characters
getFinal3WayWideHuge :: (FiniteBits a, GV.Vector v a) => MR.MetricRepresentation a -> v a -> v a -> v a -> v a
getFinal3WayWideHuge whTCM parentFinal descendantLeftPrelim descendantRightPrelim =
   let newFinal = removeGapAndNil $ GV.zipWith3 (local3WayWideHuge whTCM) parentFinal descendantLeftPrelim descendantRightPrelim
   in
   newFinal

-- | local3WayWideHuge takes tripples for wide and huge sequence types and returns median
local3WayWideHuge :: (FiniteBits a) => MR.MetricRepresentation a -> a -> a -> a -> a
local3WayWideHuge lWideTCM b c d =
   let  -- b' = if b == zeroBits then gap else b
        -- c' = if c == zeroBits then gap else c
        -- d' = if d == zeroBits then gap else d
        (median, _) = MR.retreiveThreewayTCM lWideTCM b c d
   in
   -- trace ((show b) <> " " <> (show c) <> " " <> (show d) <> " => " <> (show median))
   median

-- | local3WaySlim takes triple of SlimState and retuns median
local3WaySlim :: TCMD.DenseTransitionCostMatrix -> SlimState -> SlimState -> SlimState -> SlimState
local3WaySlim lSlimTCM b c d =
 -- trace ("L3WS: " <> (show (b,c,d))) (
 let  -- b' = if b == zeroBits then gap else b
      -- c' = if c == zeroBits then gap else c
      -- d' = if d == zeroBits then gap else d
 in
 -- trace ("L3WS: " <> (show (b',c',d'))) (
 let (median, _) = TCMD.lookupThreeway lSlimTCM b c d
 in
 -- trace ("3way: " <> (show b) <> " " <> (show c) <> " " <> (show d) <> " => " <> (show (median, cost)))
 median
 -- )

 -- | generalSequenceDiff  takes two sequence elemental bit types and retuns min and max integer
-- cost differences using matrix values
-- if value has no bits on--it is set to 0th bit on for GAP
generalSequenceDiff :: (Show a, FiniteBits a) => S.Matrix Int -> Int -> a -> a -> (Int, Int)
generalSequenceDiff thisMatrix numStates uState vState =
    -- trace ("GSD: " <> (show (numStates, uState, vState))) (
    let gapIfNil x
            | popCount x == 0 = (x `xor` x) `setBit` fromEnum gapIndex
            | otherwise = x
        uState' = gapIfNil uState
        vState' = gapIfNil vState
        uStateList = fmap snd $ filter fst $ zip (fmap (testBit uState') [0.. numStates - 1]) [0.. numStates - 1]
        vStateList = fmap snd $ filter fst $ zip (fmap (testBit vState') [0.. numStates - 1]) [0.. numStates - 1]
        uvCombinations = cartProd uStateList vStateList
        costOfPairs = fmap (thisMatrix S.!) uvCombinations
    in
    -- trace ("GSD: " <> (show uStateList) <> " " <> (show vStateList) <> " min " <> (show $ minimum costOfPairs) <> " max " <> (show $  maximum costOfPairs))
    (minimum costOfPairs, maximum costOfPairs)
    -- )


-- | getNoGapPrelimContext takes gaps and nils out of left, median, and right of gapped structure
getNoGapPrelimContext
  :: ( FiniteBits e
     , GV.Vector v e
     )
  => OpenDynamicCharacter v e
  -> OpenDynamicCharacter v e
getNoGapPrelimContext prelimContext =
   let lhs = extractMediansLeft prelimContext
       med = extractMedians prelimContext
       rhs = extractMediansRight prelimContext
   in  (lhs, med, rhs)
