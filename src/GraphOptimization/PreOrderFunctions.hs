{- |
Module      :  PostOrderFunctions.hs
Description :  Module specifying pre-order graph functions
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

{-
ToDo:
   Add parallel optimization overblocks and characters?
-}


module GraphOptimization.PreOrderFunctions  ( createFinalAssignmentOverBlocks
                                            ) where

import           Types.Types
import qualified Utilities.LocalGraph as LG
import GeneralUtilities
import qualified Debug.Debug as D
import qualified DirectOptimization.PreOrder as DOP
import qualified GraphOptimization.Medians as M
import qualified Data.Vector.Generic  as GV
import qualified Data.Vector                                                 as V
import qualified Data.BitVector.LittleEndian as BV
import           Data.Bits                 ((.&.), (.|.))
import           Data.Maybe
import qualified Data.List as L
import Debug.Trace



-- | createFinalAssignment takes vertex data (child or current vertex) and creates the final 
-- assignment from parent (if not root or leaf) and 'child' ie current vertex
-- if root or leaf preliminary is assigned to final
   -- need to watch zipping for missing sequence data
createFinalAssignmentOverBlocks :: NodeType
                                -> VertexBlockData 
                                -> VertexBlockData 
                                -> V.Vector (V.Vector CharInfo) 
                                -> Bool
                                -> VertexBlockData
createFinalAssignmentOverBlocks childType childBlockData parentBlockData blockCharInfoVect isLeft =
   -- if root or leaf final assignment <- preliminary asssignment
   let childParentBlockCharInfoTriple = D.debugVectorZip3 childBlockData parentBlockData blockCharInfoVect
       rootBlockPair = D.debugVectorZip childBlockData blockCharInfoVect
   in
   fmap (assignFinal childType isLeft) childParentBlockCharInfoTriple

  -- | assignFinal takes a vertex type and single block of zip3 of child info, parent info, and character type 
-- to create pre-order assignments
assignFinal :: NodeType -> Bool -> (V.Vector CharacterData, V.Vector CharacterData, V.Vector CharInfo)-> V.Vector CharacterData
assignFinal childType isLeft inTriple@(childCharacterVect, parentCharacterVect, charInfoVect) =
   let childCharInfoTripleList = D.debugVectorZip3 childCharacterVect parentCharacterVect charInfoVect
   in
   fmap (setFinal childType isLeft) childCharInfoTripleList


-- | setFinal takes a vertex type and single character of zip3 of child info, parent info, and character type 
-- to create pre-order assignments
   -- | setFinalHTU takes a single character and its parent and sets the final state to prelim based 
-- on character info. 
-- non exact charcaters are vectors of characters of same type
-- this does the same things for seqeunce types, but also 
-- performs preorder logic for exact characters
setFinal :: NodeType -> Bool -> (CharacterData, CharacterData, CharInfo) -> CharacterData
setFinal childType isLeft inData@(childChar, parentChar, charInfo) =
   let localCharType = charType charInfo
       symbolCount = toEnum $ length $ costMatrix charInfo
   in

   -- Three cases, Root, leaf, HTU
   if childType == RootNode then 

      if localCharType == Add then 
         --trace ("Root " ++ show (fst3 $ rangePrelim childChar))
         childChar {rangeFinal = fst3 $ rangePrelim childChar}

      else if localCharType == NonAdd then childChar {stateBVFinal = fst3 $ stateBVPrelim childChar}

      else if localCharType == Matrix then childChar {matrixStatesFinal = setMinCostStatesMatrix (fromEnum symbolCount) (localCostVect childChar) (matrixStatesPrelim childChar)}

      -- need to set both final and alignment for sequence characters
      else if (localCharType == SlimSeq) || (localCharType == NucSeq) then 
         trace ("root " ++ show (slimPrelim childChar, slimPrelim childChar, slimGapped childChar)) 
         childChar {slimFinal = slimPrelim childChar, slimAlignment = slimGapped childChar}
         
      else if (localCharType == WideSeq) || (localCharType == AminoSeq) then 
         childChar {wideFinal = widePrelim childChar, wideAlignment = wideGapped childChar}
         
      else if localCharType == HugeSeq then 
         childChar {hugeFinal = hugePrelim childChar, hugeAlignment = hugeGapped childChar}
         
      else error ("Unrecognized/implemented character type: " ++ show localCharType)

   else if childType == LeafNode then 

      if localCharType == Add then childChar {rangeFinal = fst3 $ rangePrelim childChar}

      else if localCharType == NonAdd then childChar {stateBVFinal = fst3 $ stateBVPrelim childChar}

      else if localCharType == Matrix then 
         childChar {matrixStatesFinal = setMinCostStatesMatrix (fromEnum symbolCount) (V.replicate  (fromEnum symbolCount) 0) (matrixStatesPrelim childChar)}

      -- need to set both final and alignment for sequence characters
      else if (localCharType == SlimSeq) || (localCharType == NucSeq) then 
         let -- finalAlignment = DOP.preOrderLogic symbolCount isLeft (slimAlignment parentChar) (slimGapped parentChar) (slimGapped childChar)
             finalAlignment = DOP.preOrderLogic symbolCount isLeft (slimAlignment parentChar) (slimGapped parentChar) (slimGapped childChar)
             finalNoGaps = M.createUngappedMedianSequence (fromEnum symbolCount) finalAlignment
         in
         trace ("Leaf " ++ show (slimPrelim childChar, finalNoGaps, finalAlignment, slimGapped childChar, slimAlignment parentChar))
         childChar {slimFinal = finalNoGaps, slimAlignment = finalAlignment}
         
      else if (localCharType == WideSeq) || (localCharType == AminoSeq) then 
         let finalAlignment = DOP.preOrderLogic symbolCount isLeft (wideAlignment parentChar) (wideGapped parentChar) (wideGapped childChar)
             finalNoGaps = M.createUngappedMedianSequence (fromEnum symbolCount) finalAlignment
         in
         childChar {wideFinal = finalNoGaps, wideAlignment = finalAlignment}
         
      else if localCharType == HugeSeq then 
         let finalAlignment = DOP.preOrderLogic symbolCount isLeft (hugeAlignment parentChar) (hugeGapped parentChar) (hugeGapped childChar)
             finalNoGaps = M.createUngappedMedianSequence (fromEnum symbolCount) finalAlignment
         in
         childChar {hugeFinal = finalNoGaps, hugeAlignment = finalAlignment}
         
      else error ("Unrecognized/implemented character type: " ++ show localCharType)

   else if childType == TreeNode then
      
      if localCharType == Add then 
         -- add logic for pre-order
         let finalAssignment = additivePreorder (rangePrelim childChar) (rangeFinal parentChar)
         in
         childChar {rangeFinal = finalAssignment}

      else if localCharType == NonAdd then 
         -- add logic for pre-order
         let finalAssignment = nonAdditivePreorder (stateBVPrelim childChar) (stateBVFinal parentChar)
         in
         childChar {stateBVFinal = finalAssignment}

      else if localCharType == Matrix then 
         -- add logic for pre-order
         let finalAssignment = matrixPreorder isLeft (matrixStatesPrelim childChar) (matrixStatesFinal parentChar)
         in
         childChar {matrixStatesFinal = finalAssignment}

      -- need to set both final and alignment for sequence characters
      else if (localCharType == SlimSeq) || (localCharType == NucSeq) then 
         let finalGapped = DOP.preOrderLogic symbolCount isLeft (slimAlignment parentChar) (slimGapped parentChar) (slimGapped childChar)
             finalNoGaps = M.createUngappedMedianSequence (fromEnum symbolCount) finalGapped
         in 
         trace ("HTU " ++ show (slimPrelim childChar, finalNoGaps, finalGapped, slimGapped childChar, slimAlignment parentChar)) 
         childChar {slimFinal = finalNoGaps, slimAlignment = finalGapped}
         
      else if (localCharType == WideSeq) || (localCharType == AminoSeq) then 
         let finalGapped = DOP.preOrderLogic symbolCount isLeft (wideAlignment parentChar) (wideGapped parentChar) (wideGapped childChar)
             finalNoGaps = M.createUngappedMedianSequence (fromEnum symbolCount) finalGapped
         in
         childChar {wideFinal = finalNoGaps, wideAlignment = finalGapped}
         
      else if localCharType == HugeSeq then 
         let finalGapped@(finalGField, _, _) = DOP.preOrderLogic symbolCount isLeft (hugeAlignment parentChar) (hugeGapped parentChar) (hugeGapped childChar)
             --  should be like slim and wide--but useing second becuae of error on post order for huge characters
             --gapChar = M.getGapBV (fromEnum symbolCount)
             --finalNoGaps =  GV.filter (M.notGapNought gapChar) finalGField
             finalNoGaps = M.createUngappedMedianSequence (fromEnum symbolCount) finalGapped 

         in
         childChar {hugeFinal = finalNoGaps, hugeAlignment = finalGapped}
         
      else error ("Unrecognized/implemented character type: " ++ show localCharType)

   else error ("Node type should not be here (pre-order on tree node only): " ++ show  childType)


-- |  additivePreorder assignment takes preliminary triple of child (= current vertex) and
-- final states of parent to create preorder final states of child
additivePreorder :: (V.Vector (Int, Int), V.Vector (Int, Int), V.Vector (Int, Int)) -> V.Vector (Int, Int) ->  V.Vector (Int, Int) 
additivePreorder vertexData@(nodePrelim, leftChild, rightChild) parentFinal =
   if null nodePrelim then mempty
   else 
      let allFour = D.debugVectorZip4 nodePrelim leftChild rightChild parentFinal
      in
      fmap makeAdditiveCharacterFinal allFour

-- |  nonAdditivePreorder assignment takes preliminary triple of child (= current vertex) and
-- final states of parent to create preorder final states of child
nonAdditivePreorder :: (V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector) -> V.Vector BV.BitVector -> V.Vector BV.BitVector
nonAdditivePreorder vertexData@(nodePrelim, leftChild, rightChild) parentFinal =
   if null nodePrelim then mempty
   else 
      let allFour = D.debugVectorZip4 nodePrelim leftChild rightChild parentFinal
      in
      fmap makeNonAdditiveCharacterFinal allFour


-- | matrixPreorder assigment akes preliminary matrix states of child (= current vertex) and
-- final states of parent to create preorder final states of child
-- th eboolean says whether the node is a 'left' node or right based on bitvetor label
matrixPreorder :: Bool -> V.Vector (V.Vector MatrixTriple) -> V.Vector (V.Vector MatrixTriple) -> V.Vector (V.Vector MatrixTriple)
matrixPreorder isLeft nodePrelim parentFinal =
   if null nodePrelim then mempty
   else 
      let bothTwo = D.debugVectorZip nodePrelim parentFinal
      in
      fmap (makeMatrixCharacterFinal isLeft) bothTwo


-- | makeAdditiveCharacterFinal takes vertex preliminary, and child preliminary states with well as parent final state
-- and constructs final state assignment
makeAdditiveCharacterFinal :: ((Int, Int), (Int, Int), (Int, Int), (Int, Int)) -> (Int, Int)
makeAdditiveCharacterFinal inData@(nodePrelim, leftChild, rightChild, parentFinal) = 
   -- From Wheeler (20012) after Goloboff (1993) 
   let interNodeParent = intervalIntersection nodePrelim parentFinal
   in
   -- trace (show inData) (
   -- Rule 1
   if (interNodeParent /= Nothing) && (fromJust interNodeParent == parentFinal) then 
      -- trace ("R1 " ++ show parentFinal) 
      parentFinal
   -- Rule 2
   else if ((leftChild `intervalUnion` rightChild) `intervalIntersection` parentFinal) /= Nothing then
      let xFactor = ((leftChild `intervalUnion` rightChild) `intervalUnion` nodePrelim) `intervalIntersection` parentFinal
      in
      if xFactor == Nothing then error ("I don't think this should happen in makeAdditiveCharacterFinal" ++ show inData)
      else 
         if (fromJust xFactor) `intervalIntersection` nodePrelim /= Nothing then 
            -- trace ("R2a " ++ show (fromJust xFactor)) 
            fromJust xFactor
         else 
            -- trace ("Rb " ++ show (lciClosest (fromJust xFactor) nodePrelim)) 
            lciClosest (fromJust xFactor) nodePrelim

   -- Rule 3
   else 
      let unionLR = leftChild `intervalUnion` rightChild
          closestPtoA = stateFirstClosestToSecond nodePrelim parentFinal
          closestULRtoA = stateFirstClosestToSecond unionLR parentFinal
      in
      -- trace ("R3 " ++ show (min closestPtoA closestULRtoA, max closestPtoA closestULRtoA)) 
      (min closestPtoA closestULRtoA, max closestPtoA closestULRtoA)
   -- )

-- | stateFirstClosestToSecond takes teh states of the first interval and finds the state wiht smallest distance 
-- to either state in the second
 -- assumes a <= b, x<= y
stateFirstClosestToSecond :: (Int, Int) -> (Int, Int) -> Int
stateFirstClosestToSecond firstInt@(a,b) secondInt@(x,y) =
   let distASecond = if x > b then x - a
                     else if y < a then a - y
                     else error ("I don't think this should happen in makeAdditiveCharacterFinal" ++ show (a,b,x,y))
       distBSecond = if x > b then x - b
                     else if y < a then b - y
                     else error ("I don't think this should happen in makeAdditiveCharacterFinal" ++ show (a,b,x,y))
   in
   if distASecond <= distBSecond then a
   else b

-- | lciClosest returns thhe "largest closed interval" between the first interval
-- and the closest state in the second interval
 -- assumes a <= b, x<= y
lciClosest :: (Int, Int) -> (Int, Int) -> (Int, Int)
lciClosest xFactor@(a,b) nodePrelim@(x,y) = 
   if x > b then (a,x)
   else if y < a then (y,b)
   else error ("I don't think this should happen in lciClosest" ++ show (a,b,x,y))

 -- | intervalIntersection is bit-analogue intersection for additive character operations
 -- takes two intervals and returnas intesection
 -- Nothing signifies an empty intersection
 -- assumes a <= b, x<= y
intervalIntersection :: (Int, Int) -> (Int, Int) -> Maybe (Int, Int)
intervalIntersection (a,b) (x,y) = 
   let newPair = (max a x, min b y)
   in
   if max a x > min b y then Nothing
   else Just newPair

  
-- | intervalUnion is bit-analogue union for additive character operations
-- takes two intervals and returnas union
 -- assumes a <= b, x<= y
intervalUnion :: (Int, Int) -> (Int, Int) -> (Int, Int)
intervalUnion (a,b) (x,y) = (min a x, max b y)

-- | largestClosedInterval is the maximum interval created from two intervals
largestClosedInterval :: (Int, Int) -> (Int, Int) -> (Int, Int)
largestClosedInterval = intervalUnion

-- | makeNonAdditiveCharacterFinal takes vertex preliminary, and child preliminary states with well as parent final state
-- and constructs final state assignment
makeNonAdditiveCharacterFinal :: (BV.BitVector, BV.BitVector, BV.BitVector, BV.BitVector) -> BV.BitVector
makeNonAdditiveCharacterFinal inData@(nodePrelim, leftChild, rightChild, parentFinal) = 
   -- From Wheeler (2012) after Fitch (1971)
   trace (show inData) (
   if (nodePrelim .&. parentFinal) == parentFinal then 
      --trace ("R1 " ++ show parentFinal) 
      parentFinal
   else if (BV.isZeroVector (leftChild .&. rightChild)) && (leftChild .|. rightChild) == nodePrelim then 
      --trace ("R2 " ++ show (nodePrelim .|. parentFinal)) 
      nodePrelim .|. parentFinal 
   else 
      -- trace ("R3 " ++ show (nodePrelim .|.  (leftChild .&. parentFinal) .|. (rightChild .&. parentFinal))) 
      nodePrelim .|.  (leftChild .&. parentFinal) .|. (rightChild .&. parentFinal)
   )

-- | makeMatrixCharacterFinal vertex preliminaryavnd parent final state
-- and constructs final state assignment
-- really just tracks the states on a traceback and sets the cost to maxBound:: Int for states not in the traceback
-- path
-- Bool for left right node
makeMatrixCharacterFinal :: Bool -> (V.Vector MatrixTriple, V.Vector MatrixTriple) -> V.Vector MatrixTriple
makeMatrixCharacterFinal isLeft inData@(nodePrelim, parentFinal) = 
   let numStates = length nodePrelim
       stateIndexList = V.fromList [0..(numStates - 1)]
       (stateCostList, stateLeftChildList, stateRightChildList) = V.unzip3 parentFinal
       (prelimStateCostList, prelimStateLeftChildList, prelimStateRightChildList) = V.unzip3 nodePrelim
       allThree = if isLeft then D.debugVectorZip3 stateCostList stateLeftChildList stateIndexList
                  else D.debugVectorZip3 stateCostList stateRightChildList stateIndexList
       bestParentThree = V.filter ((/= (maxBound :: StateCost)). fst3) allThree
       bestPrelimStates = L.sort $ L.nub $ concat $ fmap snd3 bestParentThree
       allFour = D.debugVectorZip4 prelimStateCostList  prelimStateLeftChildList prelimStateRightChildList stateIndexList
       finalBestTriple = V.filter ((/= (maxBound :: StateCost)).fst3) $ fmap (setCostsAndStates bestPrelimStates) allFour
   in
   finalBestTriple

-- | setCostsAndStates takes a list of states that are in teh set of 'best' and a four-tuple
-- of a matrix triple annd a fourth field of the state index
-- if the state is in the list of `best' indices it is kept and not if it isn't
setCostsAndStates :: [Int] -> (StateCost, [ChildStateIndex], [ChildStateIndex], Int) -> (StateCost, [ChildStateIndex], [ChildStateIndex])
setCostsAndStates bestPrelimStates inQuad@(cost, leftChildState, rightChildStates, stateIndex) = 
   if stateIndex `elem` bestPrelimStates then (stateIndex, leftChildState, rightChildStates)
   else (maxBound :: StateCost, leftChildState, rightChildStates)



-- | setMinCostStatesMatrix  sets the cost of non-minimal cost states to maxBounnd :: StateCost (Int) 
setMinCostStatesMatrix ::  Int -> V.Vector StateCost -> V.Vector (V.Vector MatrixTriple) ->  V.Vector (V.Vector MatrixTriple)
setMinCostStatesMatrix numStates inCostVect inStateVect = 
    fmap (V.filter ((/= (maxBound :: StateCost)).fst3)) $ fmap (nonMinCostStatesToMaxCost (V.fromList [0.. (numStates - 1)])) $ D.debugVectorZip inCostVect inStateVect 

-- | nonMinCostStatesToMaxCost takes an individual pair of minimum state cost and matrix character triple 
-- retiurns a new character with the states cost either the minium value or maxBound iof not
-- this only really useful at root--other vertices minimu costs may not be paert of the
-- miniumm cost assignment, but may be useful heuristically
nonMinCostStatesToMaxCost :: V.Vector StateCost -> (StateCost, V.Vector MatrixTriple) -> V.Vector MatrixTriple
nonMinCostStatesToMaxCost stateIndexList (minStateCost, tripleVect) = 
   let result = fmap (modifyStateCost minStateCost) $ V.zip tripleVect stateIndexList
   in 
   -- trace ((show stateIndexList) ++ " " ++ (show $ V.zip tripleVect stateIndexList))
   result
      where
         modifyStateCost d ((a,b,c), e) = if a == d then (e,b,c)
                                          else (maxBound :: StateCost ,b,c)


