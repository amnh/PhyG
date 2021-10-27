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
                                  , createUngappedMedianSequence
                                  , intervalAdd
                                  , interUnion
                                  , addMatrix
                                  , getDOMedian
                                  , getDOMedianCharInfo
                                  , pairwiseDO
                                  ) where

import           Bio.DynamicCharacter
import           Data.Bits
import qualified Data.BitVector.LittleEndian                                 as BV
import           Data.Foldable
import qualified Data.MetricRepresentation                                   as MR
import qualified Data.TCM.Dense                                              as TCMD
import qualified Data.Vector                                                 as V
import qualified Data.Vector.Generic                                         as GV
import           Data.Word
import           DirectOptimization.Pairwise
import           GeneralUtilities
import qualified SymMatrix                                                   as S
import           Types.Types
import           Debug.Trace


--import qualified Data.Alphabet as DALPH

-- | median2 takes the vectors of characters and applies media2 to each
-- character
-- for parallel fmap over all then parallelized by type and sequences
-- used for distances and post-order assignments
median2 ::   V.Vector CharacterData -> V.Vector CharacterData -> V.Vector CharInfo -> V.Vector (CharacterData, VertexCost)
median2 = V.zipWith3 median2Single


-- | median2NonExact takes the vectors of characters and applies median2NonExact to each
-- character for parallel fmap over all then parallelized by type and sequences
-- this only reoptimized the nonexact characters (sequence characters for now, perhpas otehrs later)
-- and takes the existing optimization for exact (Add, NonAdd, Matrix) for the others.
median2NonExact :: V.Vector CharacterData -> V.Vector CharacterData -> V.Vector CharInfo -> V.Vector (CharacterData, VertexCost)
median2NonExact = V.zipWith3 median2SingleNonExact


-- | median2Single takes character data and returns median character and cost
-- median2single assumes that the character vectors in the various states are the same length
-- that is--all leaves (hencee other vertices later) have the same number of each type of character
-- used for post-order assignments
median2Single :: CharacterData -> CharacterData -> CharInfo -> (CharacterData, VertexCost)
median2Single firstVertChar secondVertChar inCharInfo =
    let thisType    = charType inCharInfo
        thisWeight  = weight inCharInfo
        thisMatrix  = costMatrix inCharInfo
        thisSlimTCM = slimTCM inCharInfo
        thisWideTCM = wideTCM inCharInfo
        thisHugeTCM = hugeTCM inCharInfo
        thisActive  = activity inCharInfo
    in
    if not thisActive then (firstVertChar, 0)
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

    else if thisType `elem` [SlimSeq, NucSeq, AminoSeq, WideSeq, HugeSeq] then
      let newCharVect = getDOMedian thisWeight thisMatrix thisSlimTCM thisWideTCM thisHugeTCM thisType firstVertChar secondVertChar
      in
      (newCharVect, localCost  newCharVect)

    else error ("Character type " ++ show thisType ++ " unrecongized/not implemented")


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
    if ((not thisActive || (thisType == Add)) || (thisType == NonAdd)) || (thisType == Matrix) then (dummyStaticCharacter, 0) else (if thisType `elem` [SlimSeq, NucSeq, AminoSeq, WideSeq, HugeSeq] then
                                                                                                                                 let newCharVect = getDOMedian thisWeight thisMatrix thisSlimTCM thisWideTCM thisHugeTCM thisType firstVertChar secondVertChar
                                                                                                                                 in
                                                                                                                                 (newCharVect, localCost  newCharVect)

                                                                                                                               else error ("Character type " ++ show thisType ++ " unrecongized/not implemented"))


-- | localOr wrapper for BV.or for vector elements
localOr :: BV.BitVector -> BV.BitVector -> BV.BitVector
localOr lBV rBV = lBV .|. rBV

-- | localAnd wrapper for BV.and for vector elements
localAnd :: BV.BitVector -> BV.BitVector -> BV.BitVector
localAnd lBV rBV = lBV .&. rBV

-- | localAndOr takes the intesection vect and union vect elements
-- and return intersection is /= 0 otherwise union
localAndOr ::BV.BitVector -> BV.BitVector -> BV.BitVector
localAndOr interBV unionBV = if BV.isZeroVector interBV then unionBV else interBV


-- | interUnion takes two non-additive chars and creates newCharcter as 2-median
-- in post-order pass to create preliminary states assignment
-- assumes a single weight for all
-- performs two passes though chars to get cost of assignments
-- fst3 $ rangePrelim left/rightChar due to triple in prelim
interUnion :: Double -> CharacterData -> CharacterData -> CharacterData
interUnion thisWeight leftChar rightChar =
    let intersectVect =  V.zipWith localAnd (fst3 $ stateBVPrelim leftChar) (fst3 $ stateBVPrelim rightChar)
        unionVect = V.zipWith localOr (fst3 $ stateBVPrelim leftChar) (fst3 $ stateBVPrelim rightChar)
        numUnions = V.length $ V.filter BV.isZeroVector intersectVect
        newCost = thisWeight * fromIntegral numUnions
        newStateVect = V.zipWith localAndOr intersectVect unionVect
        newCharcater = emptyCharacter { stateBVPrelim = (newStateVect, fst3 $ stateBVPrelim leftChar, fst3 $ stateBVPrelim rightChar)
                                      , localCost = newCost
                                      , globalCost = newCost + globalCost leftChar + globalCost rightChar
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
    else if lMax <= rMin then (lMax, rMin, rMin - lMax)
    else if rMax <= lMin then (rMax, lMin, lMin - rMax)
    else error ("This can't happen " ++ show inStuff)

-- | intervalAdd takes two additive chars and creates newCharcter as 2-median
-- in post-order pass to create preliminary states assignment
-- assumes a single weight for all
-- fst3 $ rangePrelim left/rightChar due to triple in prelim
intervalAdd :: Double -> CharacterData -> CharacterData -> CharacterData
intervalAdd thisWeight leftChar rightChar =
    let newRangeCosts = V.map getNewRange $ V.zip4 (V.map fst $ fst3 $ rangePrelim leftChar) (V.map snd $ fst3 $ rangePrelim leftChar) (V.map fst $ fst3 $ rangePrelim rightChar) (V.map snd $ fst3 $ rangePrelim rightChar)
        newMinRange = V.map fst3 newRangeCosts
        newMaxRange = V.map snd3 newRangeCosts
        newCost = thisWeight * fromIntegral (V.sum $ V.map thd3 newRangeCosts)
        newCharcater = emptyCharacter { rangePrelim = (V.zip newMinRange newMaxRange, fst3 $ rangePrelim leftChar, fst3 $ rangePrelim rightChar)
                                      , localCost = newCost
                                      , globalCost = newCost + globalCost leftChar + globalCost rightChar
                                      }
    in

    --trace ("Additive: " ++ (show newCost) ++ "\t" ++ (show $ rangePrelim leftChar) ++ "\t" ++ (show $ rangePrelim rightChar)
     --   ++ (show newRangeCosts))

    newCharcater

-- | getMinCostStates takes cost matrix and vector of states (cost, _, _) and retuns a list of (totalCost, best child state)
getMinCostStates :: S.Matrix Int -> V.Vector MatrixTriple -> Int -> Int -> Int -> [(Int, ChildStateIndex)]-> Int -> [(Int, ChildStateIndex)]
getMinCostStates thisMatrix childVect bestCost numStates childState currentBestStates stateIndex =
   --trace (show thisMatrix ++ "\n" ++ (show  childVect) ++ "\n" ++ show (numStates, childState, stateIndex)) (
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
-- cost matrix and calculates a new verctor n^2 in states
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
        --trace ("Matrix: " ++ (show newCost) ++ "\n\t" ++ (show $ matrixStatesPrelim firstVertChar)  ++ "\n\t" ++ (show $ matrixStatesPrelim secondVertChar) ++
        --  "\n\t" ++ (show initialMatrixVector) ++ "\n\t" ++ (show initialCostVector))

        newCharacter

-- | pairwiseDO is a wrapper around slim/wise/hugeParwiseDO to allow direct call and return of 
-- DO medians and cost.  This is used in final state assignment
pairwiseDO :: CharInfo
           -> (SlimDynamicCharacter, WideDynamicCharacter, HugeDynamicCharacter)
           -> (SlimDynamicCharacter, WideDynamicCharacter, HugeDynamicCharacter)
           -> (SlimDynamicCharacter, WideDynamicCharacter, HugeDynamicCharacter, Double)
pairwiseDO charInfo (slim1, wide1, huge1) (slim2, wide2, huge2) =
    let thisType = charType charInfo
        symbolCount = length $ costMatrix charInfo
    in
    if thisType `elem` [SlimSeq,   NucSeq]      then
        let (cost, r) = slimPairwiseDO (slimTCM charInfo) slim1 slim2
        in
        --trace ("pDO:" ++ (show (GV.length $ fst3 slim1)) ++ " " ++ (show (GV.length $ fst3 slim2))) 
        (r, mempty, mempty, weight charInfo * fromIntegral cost)

    else if thisType `elem` [WideSeq, AminoSeq] then
        let coefficient = MR.minInDelCost (wideTCM charInfo)
            gapState    = bit $ symbolCount - 1
            (cost, r) = widePairwiseDO coefficient gapState (MR.retreivePairwiseTCM $ wideTCM charInfo) wide1 wide2

        in
        (mempty, r, mempty, weight charInfo * fromIntegral cost)

    else if thisType == HugeSeq           then
        let coefficient = MR.minInDelCost (hugeTCM charInfo)
            gapState    = bit $ symbolCount - 1
            (cost, r) = hugePairwiseDO coefficient gapState (MR.retreivePairwiseTCM $ hugeTCM charInfo) huge1 huge2

        in
        (mempty, mempty, r, weight charInfo * fromIntegral cost)

    else error $ fold ["Unrecognised character type '", show thisType, "'in a DYNAMIC character branch" ]




-- | getDOMedianCharInfo  is a wrapper around getDOMedian with CharInfo-based interface
getDOMedianCharInfo :: CharInfo -> CharacterData -> CharacterData -> CharacterData
getDOMedianCharInfo charInfo = getDOMedian (weight charInfo) (costMatrix charInfo) (slimTCM charInfo) (wideTCM charInfo) (hugeTCM charInfo) (charType charInfo)

-- | getDOMedian calls PCG/POY/C ffi to create sequence median after some type wrangling
getDOMedian
  :: Double
  -> S.Matrix Int
  -> TCMD.DenseTransitionCostMatrix
  -> MR.MetricRepresentation Word64
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

    symbolCount = length thisMatrix
    newSlimCharacterData =
        let newCost     = thisWeight * fromIntegral cost
            subtreeCost = sum [ newCost, globalCost leftChar, globalCost rightChar]
            (cost, r)   = slimPairwiseDO
                thisSlimTCM  -- (slimGapped leftChar) (slimGapped rightChar)
                (slimPrelim  leftChar, slimPrelim  leftChar, slimPrelim  leftChar)
                (slimPrelim rightChar, slimPrelim rightChar, slimPrelim rightChar)
            --gapChar = getGapBV symbolCount
        in  blankCharacterData
              { --slimPrelim    = GV.filter (notGapNought gapChar) medians
                slimPrelim    = createUngappedMedianSequence symbolCount r
              , slimGapped    = r
              , localCostVect = V.singleton $ fromIntegral cost
              , localCost     = newCost
              , globalCost    = subtreeCost
              }

    newWideCharacterData =
        let newCost     = thisWeight * fromIntegral cost
            coefficient = MR.minInDelCost thisWideTCM
            gapState    = bit $ symbolCount - 1
            subtreeCost = sum [ newCost, globalCost leftChar, globalCost rightChar]
            (cost, r)   = widePairwiseDO
                coefficient
                gapState
                (MR.retreivePairwiseTCM thisWideTCM)
                (wideGapped leftChar) (wideGapped rightChar)
                -- (widePrelim  leftChar, widePrelim  leftChar, widePrelim  leftChar)
                -- (widePrelim rightChar, widePrelim rightChar, widePrelim rightChar)
            --gapChar = getGapBV symbolCount
        in  blankCharacterData
              { --widePrelim    = GV.filter (notGapNought gapChar) medians
                widePrelim    = createUngappedMedianSequence symbolCount r
              , wideGapped    = r
              , localCostVect = V.singleton $ fromIntegral cost
              , localCost     = newCost
              , globalCost    = subtreeCost
              }

    newHugeCharacterData =
        let newCost     = thisWeight * fromIntegral cost
            coefficient = MR.minInDelCost thisHugeTCM
            gapState    = bit $ symbolCount - 1
            subtreeCost = newCost + globalCost leftChar + globalCost rightChar
            (cost, r)   = hugePairwiseDO
            -- (cost, r)   = hugePairwiseDO
                coefficient
                gapState
                (MR.retreivePairwiseTCM thisHugeTCM)
                (hugeGapped leftChar) (hugeGapped rightChar)
                -- (hugePrelim  leftChar, hugePrelim  leftChar, hugePrelim  leftChar)
                -- (hugePrelim rightChar, hugePrelim rightChar, hugePrelim rightChar)
            --gapChar = getGap symbolCount
        in  blankCharacterData
              { --hugePrelim = GV.filter (notGapNought gapChar) medians -- 
                hugePrelim = createUngappedMedianSequence symbolCount r
              , hugeGapped = r
              , localCostVect = V.singleton $ fromIntegral cost
              , localCost  = newCost
              , globalCost = subtreeCost
              }

-- |
-- createUngappedMedianSequence enter `Symbol Count` (symbols from alphabet) and context
createUngappedMedianSequence :: (FiniteBits a, GV.Vector v a) => Int -> (v a, v a, v a) -> v a
createUngappedMedianSequence symbols v@(m,_,_) = GV.ifilter f m
  where
    gap = bit $ symbols - 1
    f i e = e /= gap && not (isGapped v i)

