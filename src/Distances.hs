{- |
Module      :  Distances.hs
Description :  Module specifying data types
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

module Distances (getPairwiseDistances
                , getBlockDistance
                , getPairwiseBlocDistance
                ) where

import           Types
import           Debug.Trace
import qualified Medians as M
import qualified SymMatrix as S
import qualified Data.Vector as V
import           Data.List
import           GeneralUtilities
import qualified ParallelUtilities as P


-- | getPairwiseDistances takes Processed data
-- and retuns a matrix (list of lists of Double) of pairwise
-- distances among vertices in data set over blocks ans all character types
-- sums over blocks
getPairwiseDistances :: ProcessedData ->  [[VertexCost]]
getPairwiseDistances (nameVect, blockDataVect) =
    if V.null nameVect then error "Null name vector in getPairwiseDistances"
    else if V.null blockDataVect then error "Null Block Data vector in getPairwiseDistances"
    else 
        let blockDistancesList =  V.toList $ V.map (getPairwiseBlocDistance (V.length nameVect)) blockDataVect 
            summedBlock = foldl' (S.zipWith (+)) (head blockDistancesList) (tail blockDistancesList)
        in
        trace ("Generating pairwise distances") 
        S.toFullLists summedBlock 



-- | getBlockDistance takes Block data and returns distance between
-- vertices based on block data  
-- this can be done for ;leaves only or all via the input processed
-- data leaves are first--then HTUs follow
getBlockDistance :: BlockData -> (Int, Int) -> VertexCost
getBlockDistance (_, blockCharData, blockCharInfo) (firstIndex, secondIndex) =
    let vertData = V.map snd blockCharData
        pairCost = V.sum $ V.map snd $ M.median2 $ V.zip3 (vertData V.! firstIndex) (vertData V.! secondIndex) blockCharInfo
    in
    pairCost

-- | getPairwiseBlocDistance returns pairwsie ditances among vertices for 
-- a block of data
-- this can be done for ;leaves only or all via the input processed
-- data leaves are first--then HTUs follow
getPairwiseBlocDistance :: Int -> BlockData-> S.Matrix VertexCost
getPairwiseBlocDistance  numVerts inData =
    let pairList = makeIndexPairs numVerts numVerts 0 0
        initialPairMatrix = S.fromLists $ replicate numVerts $ replicate numVerts 0.0
        pairListCosts = P.seqParMap P.myStrategyRPAR (getBlockDistance inData) pairList
        (iLst, jList) = unzip pairList
        threeList = zip3 iLst jList pairListCosts
        newMatrix = S.updateMatrix initialPairMatrix threeList
    in 
    --trace ("NM:\n" ++ (show $ S.toFullLists newMatrix))
    newMatrix
    
