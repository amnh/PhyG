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
                                    ) where

import           Bio.DynamicCharacter
import           Data.Alphabet
import           Data.Bits
import qualified Data.BitVector.LittleEndian as BV
import qualified Data.List                   as L
import           Data.Maybe
import qualified Data.Vector                 as V
import qualified Data.Vector.Generic         as GV
import           Debug.Trace
import qualified DirectOptimization.PreOrder as DOP
import           GeneralUtilities
import qualified GraphOptimization.Medians   as M
import           Types.Types
import qualified Utilities.LocalGraph        as LG
import qualified SymMatrix                   as S

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
      
   else if localCharType == Matrix then 
      let threeFinal = V.zipWith3 (threeWayMatrix (costMatrix charInfo)) (matrixStatesFinal parent1) (matrixStatesFinal parent2) (matrixStatesPrelim curNode) 
      in curNode {matrixStatesFinal = threeFinal}
      
   else if (localCharType == SlimSeq) || (localCharType == NucSeq) then
      curNode
         
   else if (localCharType == WideSeq) || (localCharType == AminoSeq) then
      curNode
      
   else if localCharType == HugeSeq then
      curNode
      
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
       threeWayStateCostList = zipWith3 a3 bestMedianP1Cost bestMedianP2Cost (fmap fst medianChildCostPairVect)
       minThreeWayCost = minimum threeWayStateCostList

       finalStateCostL = zipWith (assignBestMax minThreeWayCost maxCost) threeWayStateCostList medianChildCostPairVect

   in 
   V.fromList finalStateCostL
   where a3 a b c = a + b + c

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
getBestPairCostAndState costMatrix maxCost numStates childStateCostV medianStateIndex =
   let statecostV = V.zipWith (g costMatrix maxCost medianStateIndex) childStateCostV (V.fromList [0..(numStates - 1)])
       minStateCost = V.minimum $ fmap fst statecostV
       bestPairs = V.filter ((== minStateCost) . fst) statecostV
       bestChildStates = V.toList $ fmap snd bestPairs
   in
   (minStateCost, L.sort bestChildStates)

   where g cM mC mS pC pS = if pC == mC then (mC, pS)
                               else (cM S.! (mS, pS), pS)


