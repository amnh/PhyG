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
                                            , get2WaySlim
                                            , get2WayWideHuge
                                            , getFinal3WaySlim
                                            , getFinal3WayWideHuge
                                            ) where

import           Types.Types
import GeneralUtilities
import qualified DirectOptimization.PreOrder as DOP
import qualified Data.Vector                                                 as V
import qualified Data.BitVector.LittleEndian as BV
import           Data.Maybe
import qualified Data.List as L
import Debug.Trace
import qualified Data.TCM.Dense as TCMD
import qualified Data.MetricRepresentation as MR
import qualified Data.Vector.Generic                                         as GV
import Data.Bits
import qualified Data.Vector.Storable         as SV
import           Foreign.C.Types             (CUInt)
import qualified GraphOptimization.Medians as M


-- | createFinalAssignment takes vertex data (child or current vertex) and creates the final 
-- assignment from parent (if not root or leaf) and 'child' ie current vertex
-- if root or leaf preliminary is assigned to final
   -- need to watch zipping for missing sequence data
-- this creates the IA during preorder from which final assignments are contructed
-- via addition post and preorder passes on IA fields.
createFinalAssignmentOverBlocks :: AssignmentMethod 
                                -> NodeType
                                -> VertexBlockData 
                                -> VertexBlockData 
                                -> CharInfo 
                                -> Bool
                                -> Bool
                                -> VertexBlockData
createFinalAssignmentOverBlocks finalMethod childType childBlockData parentBlockData charInfo isLeft isOutDegree1 =
   -- if root or leaf final assignment <- preliminary asssignment
   V.zipWith (assignFinal finalMethod childType isLeft charInfo isOutDegree1) childBlockData parentBlockData

  -- | assignFinal takes a vertex type and single block of zip3 of child info, parent info, and character type 
-- to create pre-order assignments
assignFinal :: AssignmentMethod -> NodeType -> Bool -> CharInfo -> Bool -> V.Vector CharacterData -> V.Vector CharacterData -> V.Vector CharacterData
assignFinal finalMethod childType isLeft charInfo isOutDegree1 childCharacterVect parentCharacterVect =
   V.zipWith (setFinal finalMethod childType isLeft charInfo isOutDegree1) childCharacterVect parentCharacterVect

-- | setFinal takes a vertex type and single character of zip3 of child info, parent info, and character type 
-- to create pre-order assignments
   -- | setFinalHTU takes a single character and its parent and sets the final state to prelim based 
-- on character info. 
-- non exact charcaters are vectors of characters of same type
-- this does the same things for seqeunce types, but also 
-- performs preorder logic for exact characters
setFinal :: AssignmentMethod -> NodeType -> Bool -> CharInfo -> Bool -> CharacterData-> CharacterData -> CharacterData
setFinal finalMethod childType isLeft charInfo isOutDegree1 childChar parentChar =
   let localCharType = charType charInfo
       symbolCount = toEnum $ length $ costMatrix charInfo
   in
   -- Three cases, Root, leaf, HTU
   --trace ("set final:" ++ (show isLeft) ++ " " ++ (show isOutDegree1) ++ " " ++ (show $ slimAlignment parentChar) ++ " " 
   --   ++ (show $ slimGapped parentChar) ++ " " ++ (show $ slimGapped childChar)) (
   if childType == RootNode then 

      if localCharType == Add then 
         childChar {rangeFinal = fst3 $ rangePrelim childChar}

      else if localCharType == NonAdd then childChar {stateBVFinal = fst3 $ stateBVPrelim childChar}

      else if localCharType == Matrix then childChar {matrixStatesFinal = setMinCostStatesMatrix (fromEnum symbolCount) (localCostVect childChar) (matrixStatesPrelim childChar)}

      -- need to set both final and alignment for sequence characters
      else if (localCharType == SlimSeq) || (localCharType == NucSeq) then 
         let finalAssignment' = M.createUngappedMedianSequence (fromEnum symbolCount) (slimGapped childChar)
         in
         childChar {slimFinal = finalAssignment', slimAlignment = slimGapped childChar}
         
      else if (localCharType == WideSeq) || (localCharType == AminoSeq) then 
         let finalAssignment' = M.createUngappedMedianSequence (fromEnum symbolCount) (wideGapped childChar)
         in
         childChar {wideFinal = finalAssignment', wideAlignment = wideGapped childChar}
         
      else if localCharType == HugeSeq then 
         let finalAssignment' = M.createUngappedMedianSequence (fromEnum symbolCount) (hugeGapped childChar)
         in
         childChar {hugeFinal = finalAssignment', hugeAlignment = hugeGapped childChar}
         
      else error ("Unrecognized/implemented character type: " ++ show localCharType)

   else if childType == LeafNode then 
      -- since leaf no neeed to precess final alignment fields for sequence characters
      if localCharType == Add then childChar {rangeFinal = fst3 $ rangePrelim childChar}

      else if localCharType == NonAdd then childChar {stateBVFinal = fst3 $ stateBVPrelim childChar}

      else if localCharType == Matrix then 
         childChar {matrixStatesFinal = setMinCostStatesMatrix (fromEnum symbolCount) (V.replicate  (fromEnum symbolCount) 0) (matrixStatesPrelim childChar)}

      -- need to set both final and alignment for sequence characters
      else if (localCharType == SlimSeq) || (localCharType == NucSeq) then 
         let finalAlignment = DOP.preOrderLogic symbolCount isLeft (slimAlignment parentChar) (slimGapped parentChar) (slimGapped childChar)
             finalAssignment' = M.createUngappedMedianSequence (fromEnum symbolCount) (slimGapped childChar)
         in
         --trace ("Leaf " ++ show (slimPrelim childChar, slimPrelim childChar, finalAlignment, slimGapped childChar, slimAlignment parentChar))
         childChar {slimFinal = finalAssignment', slimAlignment = finalAlignment}
         
      else if (localCharType == WideSeq) || (localCharType == AminoSeq) then 
         let finalAlignment = DOP.preOrderLogic symbolCount isLeft (wideAlignment parentChar) (wideGapped parentChar) (wideGapped childChar)
             finalAssignment' = M.createUngappedMedianSequence (fromEnum symbolCount) (wideGapped childChar)
         in
         childChar {wideFinal = finalAssignment', wideAlignment = finalAlignment}
         
      else if localCharType == HugeSeq then 
         let finalAlignment = DOP.preOrderLogic symbolCount isLeft (hugeAlignment parentChar) (hugeGapped parentChar) (hugeGapped childChar)
             finalAssignment' = M.createUngappedMedianSequence (fromEnum symbolCount) (hugeGapped childChar)
         in
         childChar {hugeFinal = finalAssignment', hugeAlignment = finalAlignment}
         
      else error ("Unrecognized/implemented character type: " ++ show localCharType)

   else if (childType == TreeNode && not isOutDegree1) then
      
      if localCharType == Add then 
         -- add logic for pre-order
         let finalAssignment' = additivePreorder (rangePrelim childChar) (rangeFinal parentChar)
         in
         childChar {rangeFinal = finalAssignment'}

      else if localCharType == NonAdd then 
         -- add logic for pre-order
         let finalAssignment' = nonAdditivePreorder (stateBVPrelim childChar) (stateBVFinal parentChar)
         in
         childChar {stateBVFinal = finalAssignment'}

      else if localCharType == Matrix then 
         -- add logic for pre-order
         let finalAssignment' = matrixPreorder isLeft (matrixStatesPrelim childChar) (matrixStatesFinal parentChar)
         in
         childChar {matrixStatesFinal = finalAssignment'}

      -- need to set both final and alignment for sequence characters
      else if (localCharType == SlimSeq) || (localCharType == NucSeq) then 
         let finalGapped = DOP.preOrderLogic symbolCount isLeft (slimAlignment parentChar) (slimGapped parentChar) (slimGapped childChar)
             finalAssignmentDOGapped = slimFinal $ getDOFinal parentChar childChar charInfo
             finalAssignmentDO = if finalMethod == DirectOptimization then M.createUngappedMedianSequence (fromEnum symbolCount) (finalAssignmentDOGapped, mempty, mempty)
                                 else mempty
         in 
         childChar {slimFinal = finalAssignmentDO, slimAlignment = finalGapped}
         
      else if (localCharType == WideSeq) || (localCharType == AminoSeq) then 
         let finalGapped = DOP.preOrderLogic symbolCount isLeft (wideAlignment parentChar) (wideGapped parentChar) (wideGapped childChar)
             finalAssignmentDOGapped = wideFinal $ getDOFinal parentChar childChar charInfo
             finalAssignmentDO = if finalMethod == DirectOptimization then M.createUngappedMedianSequence (fromEnum symbolCount) (finalAssignmentDOGapped, mempty, mempty)
                                 else mempty
         in 
         childChar {wideFinal = finalAssignmentDO, wideAlignment = finalGapped}
         
      else if localCharType == HugeSeq then 
         let finalGapped = DOP.preOrderLogic symbolCount isLeft (hugeAlignment parentChar) (hugeGapped parentChar) (hugeGapped childChar)
             finalAssignmentDOGapped = hugeFinal $ getDOFinal parentChar childChar charInfo
             finalAssignmentDO = if finalMethod == DirectOptimization then M.createUngappedMedianSequence (fromEnum symbolCount) (finalAssignmentDOGapped, mempty, mempty)
                                 else mempty
         in 
         childChar {hugeFinal = finalAssignmentDO, hugeAlignment = finalGapped}
         
      else error ("Unrecognized/implemented character type: " ++ show localCharType)
   
   -- display tree indegree=outdegree=1
   -- since display trees here--indegree should be one as well
   else if isOutDegree1 then 
      -- trace ("InOut1 preorder") (
      if localCharType == Add then 
         -- add logic for pre-order
         let lFinalAssignment = (rangeFinal parentChar)
         in
         childChar {rangeFinal = lFinalAssignment}

      else if localCharType == NonAdd then 
         -- add logic for pre-order
         let lFinalAssignment = (stateBVFinal parentChar)
         in
         childChar {stateBVFinal = lFinalAssignment}

      else if localCharType == Matrix then 
         -- add logic for pre-order
         let lFinalAssignment = (matrixStatesFinal parentChar)
         in
         childChar {matrixStatesFinal = lFinalAssignment}

      -- need to set both final and alignment for sequence characters
      else if (localCharType == SlimSeq) || (localCharType == NucSeq) then 
         childChar { slimFinal = slimFinal parentChar
                   , slimAlignment = slimAlignment parentChar
                   , slimIAFinal = slimFinal parentChar}
         
      else if (localCharType == WideSeq) || (localCharType == AminoSeq) then 
         childChar { wideFinal = wideFinal parentChar
                   , wideAlignment = wideAlignment parentChar
                   , wideIAFinal = wideFinal parentChar}
       
      else if localCharType == HugeSeq then 
         childChar { hugeFinal = hugeFinal parentChar
                   , hugeAlignment =  hugeAlignment parentChar
                   , hugeIAFinal = hugeFinal parentChar}
       
      else error ("Unrecognized/implemented character type: " ++ show localCharType)
      -- )
      
   else error ("Node type should not be here (pre-order on tree node only): " ++ show  childType)
   -- )

-- | getDOFinal takes parent final, and node gapped (including its parent gapped) and performs a DO median
-- to get the finla state.  This takes place in several steps
--    1) align (DOMedian) parent final with node gapped (ie node preliminary)
--    2) propagate new gaps in naligned node prelimiinary to child gapped in node tripe (snd and thd)
--       creating a 3-way alignment with parent final and child preliminary
--    3) apply approprouae get3way for the structure
-- The final is then returned--with gaps to be filtered afterwards
-- getDOFinal :: (FiniteBits a, GV.Vector v a) => v a -> (v a, v a, v a) -> CharInfo -> v a
getDOFinal :: CharacterData -> CharacterData -> CharInfo -> CharacterData
getDOFinal parentData nodeData charInfo =
   let parentNodeChar = M.getDOMedianCharInfo charInfo parentData nodeData

       -- put "new" gaps into 2nd and thd gapped fileds of appropriate seqeunce type
       gappedLeftRightChar = makeGappedLeftRight parentNodeChar nodeData charInfo
   in
   gappedLeftRightChar

-- | makeGappedLeftRight takes an alignment parent charcater and original node character and inserts "new" gaps into nodeCharcater
makeGappedLeftRight :: CharacterData -> CharacterData -> CharInfo -> CharacterData
makeGappedLeftRight gappedLeftRight nodeChar charInfo =
   let localCharType = charType charInfo
   in
   if localCharType `elem` [SlimSeq, NucSeq] then
      let (parentGapped, leftChildGapped, rightChildGapped) = addGapsToChildren  (length $ costMatrix charInfo) (slimGapped gappedLeftRight) (slimGapped nodeChar)
          newFinalGapped = getFinal3WaySlim (slimTCM charInfo) (length $ costMatrix charInfo) parentGapped leftChildGapped rightChildGapped
      in
      nodeChar {slimFinal = newFinalGapped}

   else if localCharType `elem` [AminoSeq, WideSeq] then
      emptyCharacter

   else if localCharType == HugeSeq then 
      emptyCharacter

   else error ("Unrecognized character type: " ++ (show localCharType))

-- | addGapsToChildren pads out "new" gaps based on identity--if not identical--adds a gap based on cost matrix size
addGapsToChildren :: (FiniteBits a, GV.Vector v a) => Int -> (v a, v a, v a) -> (v a, v a, v a) -> (v a, v a, v a)
addGapsToChildren symbols (_, reGappedParentFinal, reGappedNodePrelim) (gappedNodePrelim, gappedLeftChild, gappedRightChild) =
   let (reGappedLeft, reGappedRight) = slideRegap symbols reGappedNodePrelim gappedNodePrelim gappedLeftChild gappedRightChild mempty mempty
   in
   if (GV.length reGappedParentFinal /= GV.length reGappedLeft) || (GV.length reGappedParentFinal /= GV.length reGappedRight) then error ("Vectors not smae length " 
      ++ (show (GV.length reGappedParentFinal, GV.length reGappedLeft, GV.length reGappedRight)))
   else (reGappedParentFinal, reGappedLeft, reGappedRight)

-- | slideRegap takes two version of same vectors (1st and snd) one with additional gaps and if the two aren't equal then adds gaps 
-- to the 3rd and 4th input vectors
slideRegap :: (FiniteBits a, GV.Vector v a) => Int -> v a -> v a -> v a -> v a -> [a] -> [a] -> (v a, v a)
slideRegap symbols reGappedNode gappedNode gappedLeft gappedRight newLeftList newRightList =
   if GV.null reGappedNode then (GV.fromList $ reverse newLeftList, GV.fromList $ reverse newRightList)
   else 
      let firstRGN = GV.head reGappedNode
          firstGN = GV.head  gappedNode
      in
      -- no "new gap"
      if firstRGN == firstGN then 
         slideRegap symbols (GV.tail reGappedNode) (GV.tail gappedNode) (GV.tail gappedLeft) (GV.tail gappedRight) (GV.head gappedLeft : newLeftList) (GV.head gappedRight : newRightList)
      else 
         let gap = bit $ symbols - 1
         in
         slideRegap symbols (GV.tail reGappedNode) gappedNode gappedLeft gappedRight (gap : newLeftList) (gap : newRightList)


-- | getFinal3Way takes parent final assignment (including indel characters) and descendent
-- preliminary gapped assingment from postorder and creates a gapped final assignment based on 
-- minimum cost median for the three inputs.  THis is done to preserve the ImpliedAlignment
-- information to create a final assingment with out an additional DO call to keep the 
-- creation linear in sequence length.  Since gaps remain--they must be filtered when output or 
-- used as true final sequence assignments using M.createUngappedMedianSequence
getFinal3WaySlim :: TCMD.DenseTransitionCostMatrix -> Int -> SV.Vector CUInt -> SV.Vector CUInt -> SV.Vector CUInt -> SV.Vector CUInt
getFinal3WaySlim lSlimTCM symbolCount parentFinal descendantLeftPrelim descendantRightPrelim =
   let gap = bit $ symbolCount - 1 
       newFinal = SV.zipWith3 (local3WaySlim lSlimTCM gap) parentFinal descendantLeftPrelim descendantRightPrelim
   in
   newFinal

-- | getFinal3WayWideHuge like getFinal3WaySlim but for wide and huge characters
getFinal3WayWideHuge :: (FiniteBits a, GV.Vector v a) => MR.MetricRepresentation a -> Int -> v a -> v a -> v a -> v a
getFinal3WayWideHuge whTCM symbolCount parentFinal descendantLeftPrelim descendantRightPrelim =
   let gap = bit $ symbolCount - 1 
       newFinal = GV.zipWith3 (local3WayWideHuge whTCM gap) parentFinal descendantLeftPrelim descendantRightPrelim
   in
   newFinal

-- | local3WayWideHuge takes tripples for wide and huge sequence types and returns median
local3WayWideHuge :: (FiniteBits a) => MR.MetricRepresentation a -> a-> a -> a -> a -> a
local3WayWideHuge lWideTCM gap b c d =
   let  b' = if popCount b == 0 then gap else b
        c' = if popCount c == 0 then gap else c
        d' = if popCount d == 0 then gap else d
        (median, _) = MR.retreiveThreewayTCM lWideTCM b' c' d'
   in
   -- trace ((show b) ++ " " ++ (show c) ++ " " ++ (show d) ++ " => " ++ (show median))
   median

-- | local3WaySlim takes triple of CUInt and retuns median
local3WaySlim :: TCMD.DenseTransitionCostMatrix -> CUInt -> CUInt -> CUInt -> CUInt -> CUInt
local3WaySlim lSlimTCM gap b c d =
 let  b' = if popCount b == 0 then gap else b
      c' = if popCount c == 0 then gap else c
      d' = if popCount d == 0 then gap else d
      (median, _) = TCMD.lookupThreeway lSlimTCM b' c' d'
 in
 -- trace ((show b) ++ " " ++ (show c) ++ " " ++ (show d) ++ " => " ++ (show median))
 median

-- | get2WaySlim takes two slim vectors an produces a preliminary median
get2WaySlim :: TCMD.DenseTransitionCostMatrix -> Int -> SV.Vector CUInt -> SV.Vector CUInt -> SV.Vector CUInt
get2WaySlim lSlimTCM symbolCount descendantLeftPrelim descendantRightPrelim =
   let gap = bit $ symbolCount - 1 
       median = SV.zipWith (local2WaySlim lSlimTCM gap) descendantLeftPrelim descendantRightPrelim
   in
   median

-- | local2WaySlim takes pair of CUInt and retuns median
local2WaySlim :: TCMD.DenseTransitionCostMatrix -> CUInt -> CUInt -> CUInt -> CUInt
local2WaySlim lSlimTCM gap b c =
 let  b' = if popCount b == 0 then gap else b
      c' = if popCount c == 0 then gap else c
      (median, _) = TCMD.lookupPairwise lSlimTCM b' c'
 in
 -- trace ((show b) ++ " " ++ (show c) ++ " " ++ (show d) ++ " => " ++ (show median))
 median

-- | get2WayWideHuge like get2WaySlim but for wide and huge characters
get2WayWideHuge :: (FiniteBits a, GV.Vector v a) => MR.MetricRepresentation a -> Int -> v a -> v a -> v a
get2WayWideHuge whTCM symbolCount descendantLeftPrelim descendantRightPrelim =
   let gap = bit $ symbolCount - 1 
       median = GV.zipWith (local2WayWideHuge whTCM gap) descendantLeftPrelim descendantRightPrelim
   in
   median

-- | local3WayWideHuge takes tripples for wide and huge sequence types and returns median
local2WayWideHuge :: (FiniteBits a) => MR.MetricRepresentation a -> a -> a -> a -> a
local2WayWideHuge lWideTCM gap b c =
   let  b' = if popCount b == 0 then gap else b
        c' = if popCount c == 0 then gap else c
        (median, _) = MR.retreivePairwiseTCM lWideTCM b' c'
   in
   -- trace ((show b) ++ " " ++ (show c) ++ " " ++ (show d) ++ " => " ++ (show median))
   median


-- |  additivePreorder assignment takes preliminary triple of child (= current vertex) and
-- final states of parent to create preorder final states of child
additivePreorder :: (V.Vector (Int, Int), V.Vector (Int, Int), V.Vector (Int, Int)) -> V.Vector (Int, Int) ->  V.Vector (Int, Int) 
additivePreorder (nodePrelim, leftChild, rightChild) parentFinal =
   if null nodePrelim then mempty
   else 
      V.zipWith4 makeAdditiveCharacterFinal nodePrelim leftChild rightChild parentFinal

-- |  nonAdditivePreorder assignment takes preliminary triple of child (= current vertex) and
-- final states of parent to create preorder final states of child
nonAdditivePreorder :: (V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector) -> V.Vector BV.BitVector -> V.Vector BV.BitVector
nonAdditivePreorder (nodePrelim, leftChild, rightChild) parentFinal =
   if null nodePrelim then mempty
   else 
      V.zipWith4 makeNonAdditiveCharacterFinal nodePrelim leftChild rightChild parentFinal


-- | matrixPreorder assigment akes preliminary matrix states of child (= current vertex) and
-- final states of parent to create preorder final states of child
-- th eboolean says whether the node is a 'left' node or right based on bitvetor label
matrixPreorder :: Bool -> V.Vector (V.Vector MatrixTriple) -> V.Vector (V.Vector MatrixTriple) -> V.Vector (V.Vector MatrixTriple)
matrixPreorder isLeft nodePrelim parentFinal =
   if null nodePrelim then mempty
   else 
      V.zipWith (makeMatrixCharacterFinal isLeft)  nodePrelim parentFinal


-- | makeAdditiveCharacterFinal takes vertex preliminary, and child preliminary states with well as parent final state
-- and constructs final state assignment
makeAdditiveCharacterFinal :: (Int, Int) -> (Int, Int) -> (Int, Int) -> (Int, Int) -> (Int, Int)
makeAdditiveCharacterFinal nodePrelim leftChild rightChild parentFinal = 
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
      if xFactor == Nothing then error ("I don't think this should happen in makeAdditiveCharacterFinal" ++ (show (nodePrelim, leftChild, rightChild, parentFinal)))
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
stateFirstClosestToSecond (a,b) (x,y) =
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
lciClosest (a,b) (x,y) = 
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

{-
-- | largestClosedInterval is the maximum interval created from two intervals
largestClosedInterval :: (Int, Int) -> (Int, Int) -> (Int, Int)
largestClosedInterval = intervalUnion
-}

-- | makeNonAdditiveCharacterFinal takes vertex preliminary, and child preliminary states with well as parent final state
-- and constructs final state assignment
makeNonAdditiveCharacterFinal :: BV.BitVector -> BV.BitVector-> BV.BitVector-> BV.BitVector -> BV.BitVector
makeNonAdditiveCharacterFinal nodePrelim leftChild rightChild parentFinal = 
   -- From Wheeler (2012) after Fitch (1971)
   -- trace (show inData) (
   if (nodePrelim .&. parentFinal) == parentFinal then 
      --trace ("R1 " ++ show parentFinal) 
      parentFinal
   else if (BV.isZeroVector (leftChild .&. rightChild)) && (leftChild .|. rightChild) == nodePrelim then 
      --trace ("R2 " ++ show (nodePrelim .|. parentFinal)) 
      nodePrelim .|. parentFinal 
   else 
      -- trace ("R3 " ++ show (nodePrelim .|.  (leftChild .&. parentFinal) .|. (rightChild .&. parentFinal))) 
      nodePrelim .|.  (leftChild .&. parentFinal) .|. (rightChild .&. parentFinal)
   -- )

-- | makeMatrixCharacterFinal vertex preliminaryavnd parent final state
-- and constructs final state assignment
-- really just tracks the states on a traceback and sets the cost to maxBound:: Int for states not in the traceback
-- path
-- Bool for left right node
makeMatrixCharacterFinal :: Bool -> V.Vector MatrixTriple -> V.Vector MatrixTriple -> V.Vector MatrixTriple
makeMatrixCharacterFinal isLeft nodePrelim parentFinal = 
   let numStates = length nodePrelim
       stateIndexList = V.fromList [0..(numStates - 1)]
       (stateCostList, stateLeftChildList, stateRightChildList) = V.unzip3 parentFinal
       (_, prelimStateLeftChildList, prelimStateRightChildList) = V.unzip3 nodePrelim
       allThree = if isLeft then V.zip3 stateCostList stateLeftChildList stateIndexList
                  else V.zip3 stateCostList stateRightChildList stateIndexList
       bestParentThree = V.filter ((/= (maxBound :: StateCost)). fst3) allThree
       bestPrelimStates = L.sort $ L.nub $ concat $ fmap snd3 bestParentThree
       allFour = V.zipWith3 (setCostsAndStates bestPrelimStates)  prelimStateLeftChildList prelimStateRightChildList stateIndexList
       finalBestTriple = V.filter ((/= (maxBound :: StateCost)).fst3) allFour
   in
   finalBestTriple

-- | setCostsAndStates takes a list of states that are in teh set of 'best' and a four-tuple
-- of a matrix triple annd a fourth field of the state index
-- if the state is in the list of `best' indices it is kept and not if it isn't
setCostsAndStates :: [Int] -> [ChildStateIndex] -> [ChildStateIndex] -> Int -> (StateCost, [ChildStateIndex], [ChildStateIndex])
setCostsAndStates bestPrelimStates leftChildState rightChildStates stateIndex = 
   if stateIndex `elem` bestPrelimStates then (stateIndex, leftChildState, rightChildStates)
   else (maxBound :: StateCost, leftChildState, rightChildStates)


-- | setMinCostStatesMatrix  sets the cost of non-minimal cost states to maxBounnd :: StateCost (Int) 
setMinCostStatesMatrix ::  Int -> V.Vector StateCost -> V.Vector (V.Vector MatrixTriple) ->  V.Vector (V.Vector MatrixTriple)
setMinCostStatesMatrix numStates inCostVect inStateVect = 
    fmap (V.filter ((/= (maxBound :: StateCost)).fst3)) $ V.zipWith (nonMinCostStatesToMaxCost (V.fromList [0.. (numStates - 1)])) inCostVect inStateVect
    
-- | nonMinCostStatesToMaxCost takes an individual pair of minimum state cost and matrix character triple 
-- retiurns a new character with the states cost either the minium value or maxBound iof not
-- this only really useful at root--other vertices minimu costs may not be paert of the
-- miniumm cost assignment, but may be useful heuristically
nonMinCostStatesToMaxCost :: V.Vector StateCost -> StateCost -> V.Vector MatrixTriple -> V.Vector MatrixTriple
nonMinCostStatesToMaxCost stateIndexList minStateCost tripleVect = 
   let result = V.zipWith (modifyStateCost minStateCost) tripleVect stateIndexList
   in 
   -- trace ((show stateIndexList) ++ " " ++ (show $ V.zip tripleVect stateIndexList))
   result
      where
         modifyStateCost d (a,b,c) e = if a == d then (e,b,c)
                                          else (maxBound :: StateCost ,b,c)


