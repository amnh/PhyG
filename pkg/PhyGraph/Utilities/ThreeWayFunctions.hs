{- |
Module      :  ThreeWayFunctions.hs
Description :  Module specifying three way optimization functions for use in pre-order
               optimization of HardWired graphs and iterative pass-type optimization for Trees
Copyright   :  (c) 2022 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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

{-
ToDo:
   Add parallel optimization overblocks and characters?
-}


module Utilities.ThreeWayFunctions  ( threeMedianFinal
                                    , addGapsToChildren
                                    ) where

import           Data.Alphabet
import           Data.Bits
import qualified Data.BitVector.LittleEndian as BV
import qualified Data.List                   as L
import qualified Data.MetricRepresentation   as MR
import qualified Data.TCM.Dense              as TCMD
import qualified Data.Vector                 as V
import qualified Data.Vector.Generic         as GV
import qualified Data.Vector.Storable        as SV
import qualified Data.Vector.Unboxed         as UV
import           Data.Word
import           Foreign.C.Types             (CUInt)
import           GeneralUtilities
import qualified GraphOptimization.Medians   as M
import qualified Input.BitPack               as BP
import qualified SymMatrix                   as S
import           Types.Types
import qualified Utilities.LocalGraph        as LG
--import           Debug.Trace

-- | threeMedianFinal calculates a 3-median for data types in a single character
-- for dynamic characters this is done by 3 min-trees
-- taking best result.  Later true 3-way via ffi can be incorporated
-- first type of operation two parents and current node-- since prelim
-- has (left child preliminary, node preliminary, right child preliminary)
-- that information can be used if needed
-- since assumes in 2 out 1 only the node preliminary field is used
threeMedianFinal :: GlobalSettings -> AssignmentMethod -> CharInfo -> CharacterData -> CharacterData -> CharacterData -> CharacterData
threeMedianFinal inGS finalMethod charInfo parent1 parent2 curNode =
   let localCharType = charType charInfo

   in
   if localCharType == Add then
      let threeFinal = V.zipWith3 threeWayAdditive (rangeFinal parent1) (rangeFinal parent2) (snd3 $ rangePrelim curNode)
      in curNode {rangeFinal = threeFinal}

   else if localCharType == NonAdd then
      let threeFinal = V.zipWith3 threeWayNonAdditive (stateBVFinal parent1) (stateBVFinal parent2) (snd3 $ stateBVPrelim curNode)
      in curNode {stateBVFinal = threeFinal}

   else if localCharType `elem` packedNonAddTypes then
      let threeFinal = BP.threeWayPacked localCharType (packedNonAddFinal parent1) (packedNonAddFinal parent2) (snd3 $ packedNonAddPrelim curNode)
      in curNode {packedNonAddFinal = threeFinal}

   else if localCharType == Matrix then
      let threeFinal = V.zipWith3 (threeWayMatrix (costMatrix charInfo)) (matrixStatesFinal parent1) (matrixStatesFinal parent2) (matrixStatesPrelim curNode)
      in curNode {matrixStatesFinal = threeFinal}

   else if localCharType == AlignedSlim then
      let threeFinal = M.getFinal3WaySlim (slimTCM charInfo) (alignedSlimFinal parent1) (alignedSlimFinal parent2) (snd3 $ alignedSlimPrelim curNode)
      in curNode {alignedSlimFinal = threeFinal}

   else if localCharType == AlignedWide then
      let threeFinal = M.getFinal3WayWideHuge (wideTCM charInfo) (alignedWideFinal parent1) (alignedWideFinal parent2) (snd3 $ alignedWidePrelim curNode)
      in curNode {alignedWideFinal = threeFinal}

   else if localCharType == AlignedHuge then
      let threeFinal = M.getFinal3WayWideHuge  (hugeTCM charInfo) (alignedHugeFinal parent1) (alignedHugeFinal parent2) (snd3 $ alignedHugePrelim curNode)
      in curNode {alignedHugeFinal = threeFinal}

   else if (localCharType == SlimSeq) || (localCharType == NucSeq) then
      let threeFinal = threeWaySlim charInfo parent1 parent2 curNode
      in
      curNode { slimFinal = threeFinal
              }

   else if (localCharType == WideSeq) || (localCharType == AminoSeq) then
      let threeFinal = threeWayWide charInfo parent1 parent2 curNode
      in
      curNode { wideFinal = threeFinal
              }

   else if localCharType == HugeSeq then
      let threeFinal = threeWayHuge charInfo parent1 parent2 curNode
      in
      curNode { hugeFinal = threeFinal
              }

   else error ("Unrecognized/implemented character type: " ++ show localCharType)


-- | threeWayNonAdditive takes the union/intersection operation over 3 non additive states
threeWayNonAdditive :: BV.BitVector -> BV.BitVector -> BV.BitVector -> BV.BitVector
threeWayNonAdditive inA inB inC =
   let intersection3  = (inA .&. inB) .&. inC
       intersectionAB = inA .&. inB
       intersectionAC = inA .&. inC
       intersectionBC = inB .&. inC
       union3 = (inA .|. inB) .|. inC
   in
   if not (BV.isZeroVector intersection3) then intersection3
   else if not (BV.isZeroVector intersectionAB) then intersectionAB .|. inC
   else if not (BV.isZeroVector intersectionAC) then intersectionAC .|. inB
   else if not (BV.isZeroVector intersectionBC) then intersectionBC .|. inA
   else union3

-- | threeWayAdditive take three additive states and returns median
-- the idea is the interval between the minimum of all three maximum (min maxA, maxB, maxC)
-- and the maximum of all three minima (max minA, minB, minC)
-- ordered such that the fst of pair not greater than second
threeWayAdditive :: (Int, Int) -> (Int, Int) -> (Int, Int) -> (Int, Int)
threeWayAdditive inA@(minA, maxA) inB@(minB, maxB) inC@(minC, maxC) =
   let minOfMaxs = minimum [maxA, maxB, maxC]
       maxOfMins = maximum [minA, minB, minC]
   in
   if maxOfMins > minOfMaxs then (maxOfMins, minOfMaxs)
   else (minOfMaxs, maxOfMins)

-- | threeWayMatrix creates median best state vector from a traceback, since parents could conflict
-- on traceback does a minimization.
-- The final states of parents will have non-maximum costs and these compared
-- to the the child states with pointers to their children are set for
-- traceback from current node to child(ren) from preliminary assignment
-- since type assumes two children--they are both set to same value so if either left or right
-- is set later the process will be correct
threeWayMatrix :: S.Matrix Int -> V.Vector MatrixTriple -> V.Vector MatrixTriple -> V.Vector MatrixTriple -> V.Vector MatrixTriple
threeWayMatrix costMatrix parent1 parent2 curNode =
   let numStates = S.rows costMatrix

       -- get the costs of each state for each node, for prents non-maximal cost will be final states
       parent1StatesCost = fmap fst3 parent1
       parent2StatesCost = fmap fst3 parent2
       curNodeStatesCost = fmap fst3 curNode

       -- get the minimum cost for each state given combinations of all three nodes and the min cost child state
       minCost3States =  getMinStatePair costMatrix (maxBound :: StateCost) numStates parent1StatesCost parent2StatesCost curNodeStatesCost
       -- minStateCost = V.minimum $ fmap fst minCost3States
       -- finalStatesTriple = fmap (assignMinMaxCost minStateCost (maxBound :: StateCost)) minCost3States
   in
   minCost3States

-- | getMinStatePair takes cost matrix and state costs (vector of Int) and returns best median cost state of child for that best cost
-- if either parent or child has maxbound cost then that state get max bound cost
getMinStatePair :: S.Matrix Int -> StateCost -> Int -> V.Vector StateCost -> V.Vector StateCost -> V.Vector StateCost -> V.Vector (StateCost, [ChildStateIndex], [ChildStateIndex])
getMinStatePair costMatrix maxCost numStates p1CostV p2CostV curCostV =
   let -- get costs to parents-- will assume parent costs are 0 or max
       bestMedianP1Cost = fmap (getBestPairCost costMatrix maxCost numStates p1CostV) [0..(numStates - 1)]
       bestMedianP2Cost = fmap (getBestPairCost costMatrix maxCost numStates p1CostV) [0..(numStates - 1)]


       -- get costs to single child via preliminary states
       medianChildCostPairVect = fmap (getBestPairCostAndState costMatrix maxCost numStates curCostV) [0..(numStates - 1)]

       -- get 3 sum costs and best state value
       threeWayStateCostList = zipWith3 f bestMedianP1Cost bestMedianP2Cost (fmap fst medianChildCostPairVect)
       minThreeWayCost = minimum threeWayStateCostList

       finalStateCostL = zipWith (assignBestMax minThreeWayCost maxCost) threeWayStateCostList medianChildCostPairVect

   in
   V.fromList finalStateCostL
   where f a b c = a + b + c

-- | assignBestMax checks 3-way median state cost and if minimum sets to that otherwise sets to max
-- double 2nd field for 2-child type asumption
assignBestMax :: StateCost -> StateCost -> StateCost -> (StateCost, [ChildStateIndex]) -> (StateCost, [ChildStateIndex], [ChildStateIndex])
assignBestMax minCost maxCost stateCost (_, stateChildList) =
   if stateCost == minCost then (minCost, stateChildList, stateChildList)
   else (maxCost, stateChildList, stateChildList)

-- | getBestPairCost gets the baest cost for a state to each of parent states--does not keep parent state
getBestPairCost :: S.Matrix Int -> StateCost -> Int -> V.Vector StateCost -> Int -> StateCost
getBestPairCost costMatrix maxCost numStates parentStateCostV medianStateIndex =
   let stateCost = V.minimum $ V.zipWith (g costMatrix maxCost medianStateIndex)parentStateCostV (V.fromList [0..(numStates - 1)])
   in
   stateCost

   where g cM mC mS pC pS = if pC == mC then mC
                               else cM S.! (mS, pS)

-- | getBestPairCostAndState gets best pair of median state and chikd states based on preliminarr states of node
getBestPairCostAndState :: S.Matrix Int -> StateCost -> Int -> V.Vector StateCost -> Int -> (StateCost, [ChildStateIndex])
getBestPairCostAndState inCostMatrix maxCost numStates childStateCostV medianStateIndex =
   let statecostV = V.zipWith (g inCostMatrix maxCost medianStateIndex) childStateCostV (V.fromList [0..(numStates - 1)])
       minStateCost = V.minimum $ fmap fst statecostV
       bestPairs = V.filter ((== minStateCost) . fst) statecostV
       bestChildStates = V.toList $ fmap snd bestPairs
   in
   (minStateCost, L.sort bestChildStates)

   where g cM mC mS pC pS = if pC == mC then (mC, pS)
                               else (cM S.! (mS, pS), pS)

-- | threeWaySlim take charInfo, 2 parents, and curNOde and creates 3 median via
-- 1) 3 DO medians (choosing lowest cost median) ((p1,p2), cn), ((cn,p1), p2), and ((cn,p2), p1)
-- 2) inserting gaps to make all 3 line up
-- 3) creating 3-medians
-- 4) choosing lowest cost median
threeWaySlim ::  CharInfo -> CharacterData -> CharacterData -> CharacterData -> SV.Vector CUInt
threeWaySlim charInfo parent1 parent2 curNode =
   let -- pairwise median structures
       p1p2 = M.getDOMedianCharInfo charInfo parent1 parent2
       p1cN = M.getDOMedianCharInfo charInfo parent1 curNode
       p2cN = M.getDOMedianCharInfo charInfo parent2 curNode

       -- get 3rd to pairwise
       p1p2cN = M.getDOMedianCharInfo charInfo p1p2 curNode
       p1cNp2 = M.getDOMedianCharInfo charInfo p1cN parent2
       p2cNp1 = M.getDOMedianCharInfo charInfo p2cN parent1

       (a1,b1,c1) = addGapsToChildren (slimGapped p1p2cN) (slimGapped p1p2)
       (median1, cost1) =  get3WayGeneric (TCMD.lookupThreeway (slimTCM charInfo)) a1 b1 c1

       (a2,b2,c2) = addGapsToChildren (slimGapped p1cNp2) (slimGapped p1cN)
       (median2, cost2) =  get3WayGeneric (TCMD.lookupThreeway (slimTCM charInfo)) a2 b2 c2

       (a3,b3,c3) = addGapsToChildren (slimGapped p2cNp1) (slimGapped p2cN)
       (median3, cost3) =  get3WayGeneric (TCMD.lookupThreeway (slimTCM charInfo)) a3 b3 c3

       minCost = minimum [cost1, cost2, cost3]

   in
   if cost1 == minCost then median1
   else if cost2 == minCost then median2
   else median3



-- | threeWayWide take charInfo, 2 parents, and curNOde and creates 3 median via
-- 1) 3 DO medians (choosing lowest cost median) ((p1,p2), cn), ((cn,p1), p2), and ((cn,p2), p1)
-- 2) inserting gaps to make all 3 line up
-- 3) creating 3-medians
-- 4) choosing lowest cost median
threeWayWide ::  CharInfo -> CharacterData -> CharacterData -> CharacterData -> UV.Vector Word64
threeWayWide charInfo parent1 parent2 curNode =
   let -- pairwise median structures
       p1p2 = M.getDOMedianCharInfo charInfo parent1 parent2
       p1cN = M.getDOMedianCharInfo charInfo parent1 curNode
       p2cN = M.getDOMedianCharInfo charInfo parent2 curNode

       -- get 3rd to pairwise
       p1p2cN = M.getDOMedianCharInfo charInfo p1p2 curNode
       p1cNp2 = M.getDOMedianCharInfo charInfo p1cN parent2
       p2cNp1 = M.getDOMedianCharInfo charInfo p2cN parent1

       (a1,b1,c1) = addGapsToChildren (wideGapped p1p2cN) (wideGapped p1p2)
       (median1, cost1) =  get3WayGeneric (MR.retreiveThreewayTCM (wideTCM charInfo)) a1 b1 c1

       (a2,b2,c2) = addGapsToChildren (wideGapped p1cNp2) (wideGapped p1cN)
       (median2, cost2) =  get3WayGeneric (MR.retreiveThreewayTCM (wideTCM charInfo)) a2 b2 c2

       (a3,b3,c3)= addGapsToChildren (wideGapped p2cNp1) (wideGapped p2cN)
       (median3, cost3) =  get3WayGeneric (MR.retreiveThreewayTCM (wideTCM charInfo)) a3 b3 c3

       minCost = minimum [cost1, cost2, cost3]

   in
   if cost1 == minCost then median1
   else if cost2 == minCost then median2
   else median3

-- | threeWayHuge take charInfo, 2 parents, and curNOde and creates 3 median via
-- 1) 3 DO medians (choosing lowest cost median) ((p1,p2), cn), ((cn,p1), p2), and ((cn,p2), p1)
-- 2) inserting gaps to make all 3 line up
-- 3) creating 3-medians
-- 4) choosing lowest cost median
threeWayHuge ::  CharInfo -> CharacterData -> CharacterData -> CharacterData -> V.Vector BV.BitVector
threeWayHuge charInfo parent1 parent2 curNode =
   let -- pairwise median structures
       p1p2 = M.getDOMedianCharInfo charInfo parent1 parent2
       p1cN = M.getDOMedianCharInfo charInfo parent1 curNode
       p2cN = M.getDOMedianCharInfo charInfo parent2 curNode

       -- get 3rd to pairwise
       p1p2cN = M.getDOMedianCharInfo charInfo p1p2 curNode
       p1cNp2 = M.getDOMedianCharInfo charInfo p1cN parent2
       p2cNp1 = M.getDOMedianCharInfo charInfo p2cN parent1

       (a1,b1,c1) = addGapsToChildren (hugeGapped p1p2cN) (hugeGapped p1p2)
       (median1, cost1) =  get3WayGeneric (MR.retreiveThreewayTCM (hugeTCM charInfo)) a1 b1 c1

       (a2,b2,c2) = addGapsToChildren (hugeGapped p1cNp2) (hugeGapped p1cN)
       (median2, cost2) =  get3WayGeneric (MR.retreiveThreewayTCM (hugeTCM charInfo)) a2 b2 c2

       (a3,b3,c3) = addGapsToChildren (hugeGapped p2cNp1) (hugeGapped p2cN)
       (median3, cost3) =  get3WayGeneric (MR.retreiveThreewayTCM (hugeTCM charInfo)) a3 b3 c3

       minCost = minimum [cost1, cost2, cost3]

   in
   if cost1 == minCost then median1
   else if cost2 == minCost then median2
   else median3



-- | addGapsToChildren pads out "new" gaps based on identity--if not identical--adds a gap based on cost matrix size
addGapsToChildren :: (FiniteBits a, GV.Vector v a) => (v a, v a, v a) -> (v a, v a, v a) -> (v a, v a, v a)
addGapsToChildren  (reGappedParentFinal, _, reGappedNodePrelim) (gappedLeftChild, gappedNodePrelim, gappedRightChild) =
   let (reGappedLeft, reGappedRight) = slideRegap reGappedNodePrelim gappedNodePrelim gappedLeftChild gappedRightChild mempty mempty
   in
   if (GV.length reGappedParentFinal /= GV.length reGappedLeft) || (GV.length reGappedParentFinal /= GV.length reGappedRight) then error ("Vectors not same length "
      ++ show (GV.length reGappedParentFinal, GV.length reGappedLeft, GV.length reGappedRight))
   else (reGappedParentFinal, reGappedLeft, reGappedRight)

-- | slideRegap takes two version of same vectors (1st and snd) one with additional gaps and if the two aren't equal then adds gaps
-- to the 3rd and 4th input vectors
slideRegap :: (FiniteBits a, GV.Vector v a) => v a -> v a -> v a -> v a -> [a] -> [a] -> (v a, v a)
slideRegap reGappedNode gappedNode gappedLeft gappedRight newLeftList newRightList =
   -- trace ("SRG: " ++ (show (GV.length reGappedNode, GV.length gappedNode))) (
   if GV.null reGappedNode then (GV.fromList $ reverse newLeftList, GV.fromList $ reverse newRightList)
   else
      let firstRGN = GV.head reGappedNode
          firstGN = GV.head  gappedNode
      in

      -- gap in reGappedNode, null gappedNode is gap at end of reGappedNode
      -- can copmplete the remainder of the slide as gaps only
      if GV.null gappedNode then
        let gapList = replicate  (GV.length reGappedNode) (bit gapIndex)
        in
        (GV.fromList $ reverse (gapList ++ newLeftList), GV.fromList $ reverse (gapList ++ newRightList))

      else if firstRGN /= firstGN then
         let gap = bit gapIndex
         in
         slideRegap (GV.tail reGappedNode) gappedNode gappedLeft gappedRight (gap : newLeftList) (gap : newRightList)

      -- no "new gap"
      else -- if firstRGN == firstGN then
         slideRegap (GV.tail reGappedNode) (GV.tail gappedNode) (GV.tail gappedLeft) (GV.tail gappedRight) (GV.head gappedLeft : newLeftList) (GV.head gappedRight : newRightList)
    -- )

-- | get3WayGeneric takes thee vectors and produces a (median, cost) pair
get3WayGeneric :: (FiniteBits e, GV.Vector v e) => (e -> e -> e -> (e, Word)) -> v e -> v e -> v e -> (v e, Word)
get3WayGeneric tcm in1 in2 in3 =
   let len   = GV.length in1
       vt    = V.generate len $ \i -> tcm (in1 GV.! i) (in2 GV.! i) (in3 GV.! i)
       gen v = let med i = fst $ v V.! i in GV.generate len med
       add   = V.foldl' (\x e -> x + snd e) 0
   in  (,) <$> gen <*> add $ vt


{-Not using this now--but could.  Would need to add Aligned Types-}
-- | threeWayGeneric take charInfo, 2 parents, and curNOde and creates 3 median via
-- 1) 3 DO medians (choosing lowest cost median) ((p1,p2), cn), ((cn,p1), p2), and ((cn,p2), p1)
-- 2) inserting gaps to make all 3 line up
-- 3) creating 3-medians
-- 4) choosing lowest cost median
threeWayGeneric :: CharInfo -> CharacterData -> CharacterData -> CharacterData -> CharacterData
threeWayGeneric charInfo parent1 parent2 curNode =
   let localCharType = charType charInfo
      -- pairwise medina structures
       p1p2 = M.getDOMedianCharInfo charInfo parent1 parent2
       p1cN = M.getDOMedianCharInfo charInfo parent1 curNode
       p2cN = M.getDOMedianCharInfo charInfo parent2 curNode

       -- get 3rd to pairwise
       p1p2cN = M.getDOMedianCharInfo charInfo p1p2 curNode
       p1cNp2 = M.getDOMedianCharInfo charInfo p1cN parent2
       p2cNp1 = M.getDOMedianCharInfo charInfo p2cN parent1

       (median1Slim, median1Wide, median1Huge, cost1) =  if localCharType `elem` [SlimSeq, NucSeq]  then
                                                            let (a,b,c) = addGapsToChildren (slimGapped p1p2cN) (slimGapped p1p2)
                                                                (median, cost) = get3WayGeneric (TCMD.lookupThreeway (slimTCM charInfo)) a b c
                                                            in
                                                            (median, mempty, mempty, cost)
                                                         else if localCharType `elem` [AminoSeq, WideSeq] then
                                                            let (a,b,c) = addGapsToChildren (wideGapped p1p2cN) (wideGapped p1p2)
                                                                (median, cost) = get3WayGeneric (MR.retreiveThreewayTCM (wideTCM charInfo)) a b c
                                                            in
                                                            (mempty, median, mempty, cost)
                                                         else if localCharType == HugeSeq then
                                                            let (a,b,c) = addGapsToChildren (hugeGapped p1p2cN) (hugeGapped p1p2)
                                                                (median, cost) = get3WayGeneric (MR.retreiveThreewayTCM (hugeTCM charInfo)) a b c
                                                            in
                                                            (mempty, mempty, median, cost)
                                                         else error ("Unrecognized character type: " ++ show localCharType)

       (median2Slim, median2Wide, median2Huge, cost2) =  if localCharType `elem` [SlimSeq, NucSeq]  then
                                                            let (a,b,c) = addGapsToChildren (slimGapped p1cNp2) (slimGapped p1cN)
                                                                (median, cost) = get3WayGeneric (TCMD.lookupThreeway (slimTCM charInfo)) a b c
                                                            in
                                                            (median, mempty, mempty, cost)
                                                         else if localCharType `elem` [AminoSeq, WideSeq] then
                                                            let (a,b,c) = addGapsToChildren (wideGapped p1cNp2) (wideGapped p1cN)
                                                                (median, cost) = get3WayGeneric (MR.retreiveThreewayTCM (wideTCM charInfo)) a b c
                                                            in
                                                            (mempty, median, mempty, cost)
                                                         else if localCharType == HugeSeq then
                                                            let (a,b,c) = addGapsToChildren (hugeGapped p1cNp2) (hugeGapped p1cN)
                                                                (median, cost) = get3WayGeneric (MR.retreiveThreewayTCM (hugeTCM charInfo)) a b c
                                                            in
                                                            (mempty, mempty, median, cost)
                                                         else error ("Unrecognized character type: " ++ show localCharType)

       (median3Slim, median3Wide, median3Huge, cost3) =  if localCharType `elem` [SlimSeq, NucSeq]  then
                                                            let (a,b,c) = addGapsToChildren (slimGapped p2cNp1) (slimGapped p2cN)
                                                                (median, cost) = get3WayGeneric (TCMD.lookupThreeway (slimTCM charInfo)) a b c
                                                            in
                                                            (median, mempty, mempty, cost)
                                                         else if localCharType `elem` [AminoSeq, WideSeq] then
                                                            let (a,b,c) = addGapsToChildren (wideGapped p2cNp1) (wideGapped p2cN)
                                                                (median, cost) = get3WayGeneric (MR.retreiveThreewayTCM (wideTCM charInfo)) a b c
                                                            in
                                                             (mempty, median, mempty, cost)
                                                         else if localCharType == HugeSeq then
                                                            let (a,b,c) = addGapsToChildren (hugeGapped p2cNp1) (hugeGapped p2cN)
                                                                (median, cost) = get3WayGeneric (MR.retreiveThreewayTCM (hugeTCM charInfo)) a b c
                                                            in
                                                            (mempty, mempty, median, cost)
                                                         else error ("Unrecognized character type: " ++ show localCharType)


       minCost = minimum [cost1, cost2, cost3]
       (medianBestSlim, medianBestWide, medianBestHuge) =  if cost1 == minCost then (median1Slim, median1Wide, median1Huge)
                                                           else if cost2 == minCost then (median2Slim, median2Wide, median2Huge)
                                                           else (median3Slim, median3Wide, median3Huge)


   in
   -- set for correct data type
   if localCharType `elem` [SlimSeq, NucSeq]  then emptyCharacter {slimFinal = medianBestSlim}

   else if localCharType `elem` [AminoSeq, WideSeq] then emptyCharacter {wideFinal = medianBestWide}

   else if localCharType == HugeSeq then emptyCharacter {hugeFinal = medianBestHuge}

   else error ("Unrecognized character type: " ++ show localCharType)
