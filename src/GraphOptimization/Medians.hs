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
                                  , median2SingleNonExact
                                  ) where

import           Types.Types
import           Debug.Trace
import qualified Data.Vector as V
import qualified Data.BitVector.LittleEndian as BV
--import qualified Data.BitVector as BV
import           GeneralUtilities

import qualified ParallelUtilities as P
import qualified SymMatrix as S 
import qualified DirectOptimization.DOWrapper as DOW
import qualified Data.TCM.Dense as TCMD
import qualified Data.MetricRepresentation as MR
import qualified Bio.Character.Encodable.Dynamic.AmbiguityGroup as AG
import Data.Bits ((.&.), (.|.))

--import qualified Data.Alphabet as DALPH

-- | median2 takes the vectors of characters and applies media2 to each 
-- character 
-- for parallel fmap over all then parallelized by type and seqeunces
median2 :: V.Vector (CharacterData, CharacterData, CharInfo) -> V.Vector (CharacterData, VertexCost)
median2 inData = V.map median2Single inData


-- | median2NonExact takes the vectors of characters and applies median2NonExact to each 
-- character for parallel fmap over all then parallelized by type and seqeunces
-- this only reoptimized the nonexact characters (sequence characters for now, perhpas otehrs later)
-- and takes the existing optimization for exact (Add, NonAdd, Matrix) for the others. 
median2NonExact :: V.Vector (CharacterData, CharacterData, CharInfo) -> V.Vector (CharacterData, VertexCost)
median2NonExact inData  = V.map median2SingleNonExact inData 
    

-- | median2Single takes character data and returns median character and cost
-- median2single assumes that the character vectors in the various states are the same length
-- that is--all leaves (hencee other vertices later) have the same number of each type of character
median2Single :: (CharacterData, CharacterData, CharInfo) -> (CharacterData, VertexCost)
median2Single (firstVertChar, secondVertChar, inCharInfo) = 
    let thisType = charType inCharInfo
        thisWeight = weight inCharInfo
        thisMatrix = costMatrix inCharInfo
        thisDenseMatrix = denseTCM inCharInfo
        thisTCMMemo = memoTCM inCharInfo
        thisActive = activity inCharInfo
    in
    if thisActive == False then (firstVertChar, 0)
    else if thisType == Add then 
        let newCharVect = intervalAdd thisWeight firstVertChar secondVertChar 
        in
        (newCharVect, localCost  newCharVect)

    else if thisType == NonAdd then 
        let newCharVect = interUnion thisWeight firstVertChar secondVertChar 
        in
        (newCharVect, localCost  newCharVect)

    else if thisType == Matrix then 
      let newCharVect = addMatrix thisWeight thisMatrix firstVertChar secondVertChar 
        in
        --trace (show $ alphabet inCharInfo)
        (newCharVect, localCost  newCharVect)

    else if thisType `elem` [SmallAlphSeq, NucSeq] then 
      -- ffi to POY-C/PCG code
      let newCharVect = getDOMedian thisWeight thisMatrix thisDenseMatrix thisTCMMemo thisType firstVertChar secondVertChar
      in
      (newCharVect, localCost  newCharVect)

    else if thisType `elem` [AminoSeq, GenSeq] then 
      let newCharVect = getDOMedian thisWeight thisMatrix thisDenseMatrix thisTCMMemo thisType firstVertChar secondVertChar
      in
      (newCharVect, localCost  newCharVect)

    else error ("Character type " ++ show thisType ++ " unrecongized/not implemented")


-- | median2SingleNonExact takes character data and returns median character and cost
-- median2single assumes that the character vectors in the various states are the same length
-- that is--all leaves (hencee other vertices later) have the same number of each type of character
-- this only reoptimized the nonexact characters (sequence characters for now, perhpas otehrs later)
-- and skips optimization placing a dummy value exact (Add, NonAdd, Matrix) for the others. 
median2SingleNonExact :: (CharacterData, CharacterData, CharInfo) -> (CharacterData, VertexCost)
median2SingleNonExact (firstVertChar, secondVertChar, inCharInfo) = 
    let thisType = charType inCharInfo
        thisWeight = weight inCharInfo
        thisMatrix = costMatrix inCharInfo
        thisDenseMatrix = denseTCM inCharInfo
        thisTCMMemo = memoTCM inCharInfo
        thisActive = activity inCharInfo
        dummyStaticCharacter = CharacterData {  stateBVPrelim = V.empty
                                      , stateBVFinal = V.empty
                                      , rangePrelim = V.empty
                                      , rangeFinal = V.empty
                                      , matrixStatesPrelim = V.empty
                                      , matrixStatesFinal = V.empty
                                      , sequencePrelim = V.empty
                                      , sequenceGapped = (V.empty, V.empty, V.empty)
                                      , sequenceFinal = V.empty
                                      , localCostVect = V.singleton 0
                                      , localCost = 0.0
                                      , globalCost = 0.0
                                      }
    in
    if thisActive == False then (dummyStaticCharacter, 0)
      
    else if thisType == Add then (dummyStaticCharacter, 0)

    else if thisType == NonAdd then (dummyStaticCharacter, 0)

    else if thisType == Matrix then (dummyStaticCharacter, 0)

    else if thisType `elem` [SmallAlphSeq, NucSeq] then 
      -- ffi to POY-C/PCG code
      let newCharVect = getDOMedian thisWeight thisMatrix thisDenseMatrix thisTCMMemo thisType firstVertChar secondVertChar
      in
      (newCharVect, localCost  newCharVect)

    else if thisType `elem` [AminoSeq, GenSeq] then 
      let newCharVect = getDOMedian thisWeight thisMatrix thisDenseMatrix thisTCMMemo thisType firstVertChar secondVertChar
      in
      (newCharVect, localCost  newCharVect)

    else error ("Character type " ++ show thisType ++ " unrecongized/not implemented")


-- | localOr wrapper for BV.or for vector elements
localOr :: BV.BitVector -> BV.BitVector -> BV.BitVector
localOr lBV rBV = lBV .|. rBV

-- | localAnd wrapper for BV.and for vector elements
localAnd :: BV.BitVector -> BV.BitVector -> BV.BitVector
localAnd lBV rBV = lBV .&. rBV

-- | localAndOr takes the intesection vect and union vect elements
-- and return intersection is /= 0 otherwise union
localAndOr ::BV.BitVector -> BV.BitVector -> BV.BitVector
localAndOr interBV unionBV = if (BV.isZeroVector interBV) then unionBV else interBV


-- | interUnion takes two non-additive chars and creates newCharcter as 2-median
-- in post-order pass to create preliminary states assignment
-- assumes a single weight for all
-- performs two passes though chars to get cost of assignments
interUnion :: Double -> CharacterData -> CharacterData -> CharacterData
interUnion thisWeight leftChar rightChar =
    let intersectVect =  V.zipWith localAnd (stateBVPrelim leftChar) (stateBVPrelim rightChar)
        unionVect = V.zipWith localOr (stateBVPrelim leftChar) (stateBVPrelim rightChar)
        numUnions = V.length $ V.filter BV.isZeroVector intersectVect
        newCost = thisWeight * (fromIntegral numUnions)
        newStateVect = V.zipWith localAndOr intersectVect unionVect
        newCharcater = CharacterData {  stateBVPrelim = newStateVect
                                      , stateBVFinal = V.singleton (BV.fromBits [False])
                                      , rangePrelim = V.empty
                                      , rangeFinal = V.empty
                                      , matrixStatesPrelim = V.empty
                                      , matrixStatesFinal = V.empty
                                      , sequencePrelim = V.empty
                                      , sequenceGapped = (V.empty, V.empty, V.empty)
                                      , sequenceFinal = V.empty
                                      , localCostVect = V.singleton 0
                                      , localCost = newCost
                                      , globalCost = newCost + (globalCost leftChar) + (globalCost rightChar) 
                                      }
    in 
    
    --trace ("NonAdditive: " ++ (show numUnions) ++ " " ++ (show newCost) ++ "\t" ++ (show $ stateBVPrelim leftChar) ++ "\t" ++ (show $ stateBVPrelim rightChar) ++ "\t"
    --   ++ (show intersectVect) ++ "\t" ++ (show unionVect) ++ "\t" ++ (show newStateVect))
    
    newCharcater

-- | getNewRange takes min and max range of two additive charcaters and returns 
-- a triple of (newMin, newMax, Cost)
getNewRange :: (Int, Int, Int, Int) -> (Int, Int, Int)
getNewRange inStuff@(lMin, lMax, rMin, rMax) = 
    -- subset
    if (rMin >= lMin) && (rMax <= lMax) then (rMin, rMax, 0)
    else if  (lMin >= rMin) && (lMax <= rMax) then (lMin, lMax, 0)
    -- overlaps
    else if (rMin >= lMin) && (rMax >= lMax) && (rMin <= lMax) then (rMin, lMax,0)
    else if (lMin >= rMin) && (lMax >= rMax) && (lMin <= rMax) then (lMin, rMax,0)
    -- newInterval
    else if (lMax <= rMin) then (lMax, rMin, rMin - lMax)
    else if (rMax <= lMin) then (rMax, lMin, lMin - rMax)
    else error ("This can't happen " ++ show inStuff)

-- | intervalAdd takes two additive chars and creates newCharcter as 2-median
-- in post-order pass to create preliminary states assignment
-- assumes a single weight for all
intervalAdd :: Double -> CharacterData -> CharacterData -> CharacterData
intervalAdd thisWeight leftChar rightChar =
    let newRangeCosts = V.map getNewRange $ V.zip4 (V.map fst $ rangePrelim leftChar) (V.map snd $ rangePrelim leftChar) (V.map fst $ rangePrelim rightChar) (V.map snd $ rangePrelim rightChar)
        newMinRange = V.map fst3 newRangeCosts
        newMaxRange = V.map snd3 newRangeCosts
        newCost = thisWeight * (fromIntegral $ V.sum $ V.map thd3 newRangeCosts)
        newCharcater = CharacterData {  stateBVPrelim = V.empty
                                      , stateBVFinal = V.empty
                                      , rangePrelim = V.zip newMinRange newMaxRange
                                      , rangeFinal = V.empty
                                      , matrixStatesPrelim = V.empty
                                      , matrixStatesFinal = V.empty
                                      , sequencePrelim = V.empty
                                      , sequenceGapped = (V.empty, V.empty, V.empty)
                                      , sequenceFinal = V.empty
                                      , localCostVect = V.singleton 0
                                      , localCost = newCost
                                      , globalCost = newCost + (globalCost leftChar) + (globalCost rightChar) 
                                      }
    in 
    
    --trace ("Additive: " ++ (show newCost) ++ "\t" ++ (show $ rangePrelim leftChar) ++ "\t" ++ (show $ rangePrelim rightChar) 
     --   ++ (show newRangeCosts))
    
    newCharcater

-- | getMinCostStates takes cost matrix and vector of states (cost, _, _) and retuns a list of (toitalCost, best child state) 
getMinCostStates :: S.Matrix Int -> V.Vector MatrixTriple -> Int -> Int -> Int -> [(Int, ChildStateIndex)]-> Int -> [(Int, ChildStateIndex)]
getMinCostStates thisMatrix childVect bestCost numStates childState currentBestStates stateIndex = 
   --trace (show thisMatrix ++ "\n" ++ (show  childVect) ++ "\n" ++ show (numStates, childState, stateIndex)) (
   if childState == numStates then reverse (filter ((== bestCost).fst) currentBestStates) 
   else
      let (childCost, _, _)  = V.head childVect
          childStateCost = if childCost /= (maxBound :: Int) then childCost + (thisMatrix S.! (childState, stateIndex))
                           else (maxBound :: Int) 
      in
      if childStateCost > bestCost then getMinCostStates thisMatrix (V.tail childVect) bestCost numStates (childState + 1) currentBestStates stateIndex
      else if childStateCost == bestCost then getMinCostStates thisMatrix (V.tail childVect) bestCost numStates (childState + 1) ((childStateCost, childState) : currentBestStates) stateIndex
      else getMinCostStates thisMatrix (V.tail childVect) childStateCost numStates (childState + 1) [(childStateCost, childState)] stateIndex
    --)


-- | getNewVector takes the vector of states and costs from teh child nodes and the 
-- cost matrix and calculates a new verctor n^2 in states
getNewVector :: S.Matrix Int -> Int -> (V.Vector MatrixTriple, V.Vector MatrixTriple) -> V.Vector MatrixTriple
getNewVector thisMatrix  numStates (lChild, rChild) =
  let newStates = [0..(numStates -1)] 
      leftPairs = fmap (getMinCostStates thisMatrix lChild (maxBound :: Int) numStates 0 []) newStates
      rightPairs = fmap (getMinCostStates thisMatrix rChild (maxBound :: Int) numStates 0 []) newStates
      stateCosts = zipWith (+) (fmap fst $ fmap head leftPairs) (fmap fst $ fmap head rightPairs) 
      newStateTripleList = zip3 stateCosts (fmap (fmap snd) leftPairs) (fmap (fmap snd) rightPairs)
  in
  V.fromList newStateTripleList

-- | addMatrix thisWeight thisMatrix firstVertChar secondVertChar 
-- assumes each character has asme cost matrix 
-- Need to add approximation ala DO tcm lookup later
-- Local and global costs are based on current not necessaril;y optimal minimum cost states
addMatrix :: Double -> S.Matrix Int -> CharacterData -> CharacterData -> CharacterData
addMatrix thisWeight thisMatrix firstVertChar secondVertChar =
  if null thisMatrix then error "Null cost matrix in addMatrix"
  else 
    let numStates = length thisMatrix 
        initialMatrixVector = V.map (getNewVector thisMatrix numStates) $ V.zip (matrixStatesPrelim firstVertChar) (matrixStatesPrelim secondVertChar) 
        initialCostVector = V.map V.minimum $ V.map (V.map fst3) initialMatrixVector
        newCost = thisWeight * (fromIntegral $ V.sum initialCostVector)
        newCharcater = CharacterData {  stateBVPrelim = V.empty  
                                      , stateBVFinal = V.empty
                                      , rangePrelim = V.empty
                                      , rangeFinal = V.empty
                                      , matrixStatesPrelim = initialMatrixVector
                                      , matrixStatesFinal = V.empty
                                      , sequencePrelim = V.empty
                                      , sequenceGapped = (V.empty, V.empty, V.empty)
                                      , sequenceFinal = V.empty
                                      , localCostVect = initialCostVector
                                      , localCost = newCost  - (globalCost firstVertChar) - (globalCost secondVertChar)
                                      , globalCost = newCost 
                                      }
        in 
        --trace ("Matrix: " ++ (show newCost) ++ "\n\t" ++ (show $ matrixStatesPrelim firstVertChar)  ++ "\n\t" ++ (show $ matrixStatesPrelim secondVertChar) ++
        --  "\n\t" ++ (show initialMatrixVector) ++ "\n\t" ++ (show initialCostVector))

        newCharcater

-- | getDOMedian calls PCG/POY/C ffi to create sequcne median after some type wrangling
getDOMedian ::  Double -> S.Matrix Int -> TCMD.DenseTransitionCostMatrix -> MR.MetricRepresentation BV.BitVector
              -> CharType -> CharacterData -> CharacterData -> CharacterData
getDOMedian thisWeight thisMatrix thisDenseMatrix thisMemoTCM thisType leftChar rightChar =
  if null thisMatrix then error "Null cost matrix in addMatrix"
  else 
    let resultVect = DOW.getDOMedian thisMatrix thisDenseMatrix thisMemoTCM thisType ((sequencePrelim leftChar), (sequencePrelim rightChar))
        newStateVect = fst resultVect
        newCostVect = snd resultVect
        newCost = (thisWeight *) (fromIntegral newCostVect)
        newCharcater = CharacterData {  stateBVPrelim = V.empty  
                                      , stateBVFinal = V.empty
                                      , rangePrelim = V.empty
                                      , rangeFinal = V.empty
                                      , matrixStatesPrelim = V.empty
                                      , matrixStatesFinal = V.empty
                                      , sequencePrelim = newStateVect
                                      , sequenceGapped = (V.empty, V.empty, V.empty)
                                      , sequenceFinal = V.empty
                                      , localCostVect = V.singleton newCostVect
                                      , localCost = newCost
                                      , globalCost = newCost + (globalCost leftChar) + (globalCost rightChar) 
                                      }
    in 
    --trace ("Sequence: " ++ (show newStateVect) ++ " " ++ (show medianCost) ++ "\n" ++ (show $ stateBVPrelim leftChar) ++ "\n" ++ (show $ stateBVPrelim rightChar))
    newCharcater
    
  