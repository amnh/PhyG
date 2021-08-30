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



-- | createFinalAssignment takes vertex data (child or current vertex) and creates the final 
-- assignment from parent (if not root or leaf) and 'child' ie current vertex
-- if root or leaf preliminary is assigned to final
   -- need to watch zipping for missing sequence data
createFinalAssignmentOverBlocks :: NodeType
                                -> VertexBlockData 
                                -> VertexBlockData 
                                -> V.Vector (V.Vector CharInfo) 
                                -> VertexBlockData
createFinalAssignmentOverBlocks childType childBlockData parentBlockData blockCharInfoVect =
   -- if root or leaf final assignment <- preliminary asssignment
   let childParentBlockCharInfoTriple = D.debugVectorZip3 childBlockData parentBlockData blockCharInfoVect
       rootBlockPair = D.debugVectorZip childBlockData blockCharInfoVect
   in
   -- root node 
   if (childType == RootNode) then
      fmap assignPrelimToFinalRootBlock rootBlockPair
   
   -- leaf node triple for this because of implied alignment
   else if (childType == LeafNode) then 
      fmap assignPrelimToFinalLeafBlock childParentBlockCharInfoTriple

   -- internal node
   else 
      fmap assignFinalHTUBlock childParentBlockCharInfoTriple

-- | assignPrelimToFinalBlock takes a pair of block data (child, character info)
-- and assigns the prliminary state to the final based on character type
-- then returns the final CharacterData 
-- this is on;y used for root due to issues contructing the implied alignment for leaves and internal vertices
assignPrelimToFinalRootBlock :: (V.Vector CharacterData, V.Vector CharInfo)-> V.Vector CharacterData
assignPrelimToFinalRootBlock inPair@(childCharacterVect, charInfoVect) =
   let childCharInfoPairList = D.debugVectorZip childCharacterVect charInfoVect
   in
   fmap setFinalPrelimCharRoot childCharInfoPairList

-- | assignPrelimToFinalLeafBlock takes a triple of block data (child, parent, character info)
-- and assigns the preliminary state to the final based on character type
-- then returns the final CharacterData 
-- this is only used for leaves due to issues contructing the implied alignment for leaves and internal vertices
assignPrelimToFinalLeafBlock :: (V.Vector CharacterData, V.Vector CharacterData, V.Vector CharInfo)-> V.Vector CharacterData
assignPrelimToFinalLeafBlock inPair@(childCharacterVect, parentCharacterVect, charInfoVect) =
   let childCharInfoTripleList = D.debugVectorZip3 childCharacterVect parentCharacterVect charInfoVect
   in
   fmap setFinalPrelimCharLeaf childCharInfoTripleList

-- | assignFinalHTUBlock takes a triple of block data (child, parent, character info)
-- and assigns the preliminary state to the final based on character type
-- then returns the final CharacterData 
-- this is only used for HTUs due to issues contructing the implied alignment for leaves and internal vertices
assignFinalHTUBlock :: (V.Vector CharacterData, V.Vector CharacterData, V.Vector CharInfo)-> V.Vector CharacterData
assignFinalHTUBlock inPair@(childCharacterVect, parentCharacterVect, charInfoVect) =
   let childCharInfoTripleList = D.debugVectorZip3 childCharacterVect parentCharacterVect charInfoVect
   in
   fmap setFinalHTU childCharInfoTripleList



-- | setFinalPrelimCharRoot takes a single character and sets the final state to prelim based 
-- on character info. Alignment fields are set to gapped field
-- non exact charcaters are vectors of characters of same type
setFinalPrelimCharRoot :: (CharacterData, CharInfo) -> CharacterData
setFinalPrelimCharRoot inData@(childChar, charInfo) =
   let localCharType = charType charInfo
       nonAddPrelim = stateBVPrelim childChar
       addPrelim = rangePrelim childChar
       matrixPrelim = matrixStatesPrelim childChar
       localSlimPrelim = slimPrelim childChar
       localSlimAlignment = slimGapped childChar
       localWidePrelim = widePrelim childChar
       localWideAlignment = wideGapped childChar
       localHugePrelim = hugePrelim childChar 
       localHugeAlignment = hugeGapped childChar
   in

   if localCharType == Add then childChar {stateBVFinal = nonAddPrelim}

   else if localCharType == NonAdd then childChar {rangeFinal = addPrelim}

   else if localCharType == Matrix then childChar {matrixStatesFinal = matrixPrelim}

   -- need to set both final and alignment for sequence characters
   else if (localCharType == SlimSeq) || (localCharType == NucSeq) then childChar {slimFinal = localSlimPrelim, slimAlignment = localSlimAlignment}
      
   else if (localCharType == WideSeq) || (localCharType == AminoSeq) then childChar {wideFinal = localWidePrelim, wideAlignment = localWideAlignment}
      
   else if localCharType == HugeSeq then childChar {hugeFinal = localHugePrelim, hugeAlignment = localHugeAlignment}
      
   else error ("UNrecognized/implemented charcater type: " ++ show localCharType)

-- | setFinalPrelimCharLeaf takes a single character and its parent and sets the final state to prelim based 
-- on character info. 
-- non exact charcaters are vectors of characters of same type
-- this unlike root due to the Implied Alignment of leaves as
-- well as internal vertices
setFinalPrelimCharLeaf :: (CharacterData, CharacterData, CharInfo) -> CharacterData
setFinalPrelimCharLeaf inData@(childChar, parentChar, charInfo) =
   let localCharType = charType charInfo
       symbolCount = toEnum $ length $ costMatrix charInfo
       nonAddPrelim = stateBVPrelim childChar
       addPrelim = rangePrelim childChar
       matrixPrelim = matrixStatesPrelim childChar
       localSlimPrelim = slimPrelim childChar
       localSlimAlignment = slimGapped childChar
       localWidePrelim = widePrelim childChar
       localWideAlignment = wideGapped childChar
       localHugePrelim = hugePrelim childChar 
       localHugeAlignment = hugeGapped childChar
   in

   if localCharType == Add then childChar {stateBVFinal = nonAddPrelim}

   else if localCharType == NonAdd then childChar {rangeFinal = addPrelim}

   else if localCharType == Matrix then childChar {matrixStatesFinal = matrixPrelim}

   -- need to set both final and alignment for sequence characters
   else if (localCharType == SlimSeq) || (localCharType == NucSeq) then 
      let finalGapped = DOP.preOrderLogic symbolCount True (slimAlignment parentChar) (slimGapped parentChar) (slimGapped childChar)
          finalNoGaps = M.createUngappedMedianSequence (fromEnum symbolCount) finalGapped
      in
      childChar {slimFinal = finalNoGaps, slimAlignment = finalGapped}
      
   else if (localCharType == WideSeq) || (localCharType == AminoSeq) then 
      let finalGapped = DOP.preOrderLogic symbolCount True (wideAlignment parentChar) (wideGapped parentChar) (wideGapped childChar)
          finalNoGaps = M.createUngappedMedianSequence (fromEnum symbolCount) finalGapped
      in
      childChar {wideFinal = finalNoGaps, wideAlignment = finalGapped}
      
   else if localCharType == HugeSeq then 
      let finalGapped@(finalGField, _, _) = DOP.preOrderLogic symbolCount True (hugeAlignment parentChar) (hugeGapped parentChar) (hugeGapped childChar)
          --  should be like slim and wide--but useing second becuae of error on post order for huge characters
          -- finalNoGaps = M.createUngappedMedianSequence symbolCount finalGapped
          gapChar = M.getGapBV (fromEnum symbolCount)
          finalNoGaps = GV.filter (M.notGapNought gapChar) finalGField

      in
      childChar {hugeFinal = finalNoGaps, hugeAlignment = finalGapped}
      
   else error ("UNrecognized/implemented charcater type: " ++ show localCharType)


   -- | setFinalHTU takes a single character and its parent and sets the final state to prelim based 
-- on character info. 
-- non exact charcaters are vectors of characters of same type
-- this does the same things for seqeunce types, but also 
-- performs preorder logic for exact characters
setFinalHTU :: (CharacterData, CharacterData, CharInfo) -> CharacterData
setFinalHTU inData@(childChar, parentChar, charInfo) =
   let localCharType = charType charInfo
       symbolCount = toEnum $ length $ costMatrix charInfo
       nonAddPrelim = stateBVPrelim childChar
       addPrelim = rangePrelim childChar
       matrixPrelim = matrixStatesPrelim childChar
       localSlimPrelim = slimPrelim childChar
       localSlimAlignment = slimGapped childChar
       localWidePrelim = widePrelim childChar
       localWideAlignment = wideGapped childChar
       localHugePrelim = hugePrelim childChar 
       localHugeAlignment = hugeGapped childChar
   in

   if localCharType == Add then 
      -- add logic for pre-order
      childChar {stateBVFinal = nonAddPrelim}

   else if localCharType == NonAdd then 
      -- add logic for pre-order
      childChar {rangeFinal = addPrelim}

   else if localCharType == Matrix then 
      -- add logic for pre-order
      childChar {matrixStatesFinal = matrixPrelim}

   -- need to set both final and alignment for sequence characters
   else if (localCharType == SlimSeq) || (localCharType == NucSeq) then 
      let finalGapped = DOP.preOrderLogic symbolCount True (slimAlignment parentChar) (slimGapped parentChar) (slimGapped childChar)
          finalNoGaps = M.createUngappedMedianSequence (fromEnum symbolCount) finalGapped
      in
      childChar {slimFinal = finalNoGaps, slimAlignment = finalGapped}
      
   else if (localCharType == WideSeq) || (localCharType == AminoSeq) then 
      let finalGapped = DOP.preOrderLogic symbolCount True (wideAlignment parentChar) (wideGapped parentChar) (wideGapped childChar)
          finalNoGaps = M.createUngappedMedianSequence (fromEnum symbolCount) finalGapped
      in
      childChar {wideFinal = finalNoGaps, wideAlignment = finalGapped}
      
   else if localCharType == HugeSeq then 
      let finalGapped@(finalGField, _, _) = DOP.preOrderLogic symbolCount True (hugeAlignment parentChar) (hugeGapped parentChar) (hugeGapped childChar)
          --  should be like slim and wide--but useing second becuae of error on post order for huge characters
          gapChar = M.getGapBV (fromEnum symbolCount)
          finalNoGaps =  GV.filter (M.notGapNought gapChar) finalGField
          --finalNoGaps = M.createUngappedMedianSequence (fromEnum symbolCount) finalGapped 

      in
      childChar {hugeFinal = finalNoGaps, hugeAlignment = finalGapped}
      
   else error ("UNrecognized/implemented charcater type: " ++ show localCharType)