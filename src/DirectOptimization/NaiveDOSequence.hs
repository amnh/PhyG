{- |
Module      :  NaiveDOSequence 
Description :  Functions for parsimony DO Naive version using bitvectors (so large alphabet) 
               without Ukkonen space/time saving 
               Uses "Sequence" type to avoid O(n) of Vector.cons
               Prototype non-optimized, restricted Haskell functisns
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

module NaiveDOSequence
( naiveDO
) where

import Debug.Trace
import qualified Data.Vector as V
--import qualified Data.BitVector as BV
import Data.Bits
import qualified Data.BitVector.LittleEndian as BV
import qualified LocalSequence as LS
import qualified Data.BitVector.LittleEndian as BL
import Data.Bits ((.&.), (.|.))

data Direction = LeftDir | DownDir | DiagDir
    deriving (Read, Show, Eq)


-- | NaiveDo wraps around naive_do_BV 
-- unlimited alphabet size no space saving
naiveDO :: V.Vector BV.BitVector -> V.Vector BV.BitVector -> Int -> (V.Vector BV.BitVector, Int)
naiveDO lBV rBV inDelCost =
    -- missing data inputs
    if V.null lBV then (rBV, 0)
    else if V.null rBV then (lBV, 0)
    else 
        -- not missing
        -- get inDelCost 
        let bvLength = fromIntegral $ BV.dimension (V.head lBV) 
            --setting left most bit to 1 same purpose as inDelBit for Ukkonen
            bv1 = BV.fromBits (True : (replicate (bvLength -1) False))
            inDelBitBV = shiftL bv1 (fromIntegral $ bvLength - 1)
            (newMedianLarge, medianCostLarge) = naive_do_BV lBV rBV inDelCost inDelBitBV
        in
        --trace ("DO: " ++ (show inDelCost) ++ " " ++ (show $ V.head $ V.last thisMatrix)) (
        (newMedianLarge, medianCostLarge) 
        --)

-- | naive_do takes two input sequences and returns median sequence and cost 
--based on charInfo-1:1 for now
--to do:
--      different costs
--      left/right based on name
--      Ukkonnen
--      C via FFI
--      Affine
naive_do_BV :: V.Vector BV.BitVector  -> V.Vector BV.BitVector  -> Int -> BV.BitVector  ->  (V.Vector BV.BitVector , Int)
naive_do_BV lSeq rSeq inDelCost inDelBitBV =
    if V.null lSeq then (rSeq, 0)
    else if V.null rSeq then (lSeq, 0)
    else 
        let subCost = 1
            --this for left right constant--want longer in left for Ukkonnen
            --(lSeq, lLength, rSeq, rLength) = setLeftRight inlSeq inrSeq
            --lSeq = max inlSeq inrSeq
            --rSeq = min inlSeq inrSeq
            lLength = V.length lSeq
            rLength = V.length rSeq
            nwMatrix = LS.cons firstRow (getRows lSeq rSeq inDelCost subCost 1 firstRow inDelBitBV)
            firstRow = getFirstRow inDelCost lLength 0 0 lSeq inDelBitBV
            (cost, _, _) = (nwMatrix LS.! rLength) LS.! lLength
            median = traceback nwMatrix (V.length rSeq) (V.length lSeq) inDelBitBV
            --revNWMatrix = LS.reverse nwMatrix
            --nwMatrix = LS.snocFlip firstRow (getRows lSeq rSeq inDelCost subCost 1 firstRow inDelBitBV)
            --(cost, _, _) = (nwMatrix LS.! 0) LS.! lLength
            --median = tracebackReverse nwMatrix 0  (V.length lSeq) inDelBitBV (V.length rSeq) (V.length lSeq)
            
        in
        --trace ("NW: " ++ show cost ++ "\n" ++ (show $ fmap (fmap fst3) nwMatrix))
        (LS.toVector median, cost)

-- | traceback creates REVERSE mediian from nwMatrix, reverse to make tail
--recusive
traceback :: LS.Seq (LS.Seq (Int, BV.BitVector, Direction)) -> Int -> Int -> BV.BitVector -> LS.Seq BV.BitVector
traceback nwMatrix posL posR inDelBitBV =
        --trace ("psLR " ++ show posL ++ " " ++ show posR) (
        if posL == 0 && posR == 0 then LS.empty
        else 
            let (_, state, direction) = (nwMatrix LS.! posL ) LS.! posR
            in
                if (state /= inDelBitBV) then
                    if direction == LeftDir then LS.snocFlip  state (traceback nwMatrix (posL) (posR - 1) inDelBitBV) 
                    else if direction == DownDir then LS.snocFlip state (traceback nwMatrix (posL - 1) (posR) inDelBitBV)  
                    else LS.snocFlip state (traceback nwMatrix (posL - 1) (posR - 1) inDelBitBV)
                else 
                    if direction == LeftDir then (traceback nwMatrix (posL) (posR - 1) inDelBitBV) 
                    else if direction == DownDir then (traceback nwMatrix (posL - 1) (posR) inDelBitBV)  
                    else (traceback nwMatrix (posL - 1) (posR - 1) inDelBitBV) 

-- | getFirstRow initializes foirst row of NW matrix
getFirstRow :: Int -> Int -> Int -> Int -> V.Vector BV.BitVector  -> BV.BitVector -> LS.Seq (Int, BV.BitVector, Direction)
getFirstRow indelCost rowLength position prevCost lSeq inDelBitBV = 
    if position == (rowLength + 1) then LS.empty
    else
        if position == 0 then LS.cons (0, inDelBitBV, DiagDir) (getFirstRow indelCost rowLength (position + 1) 0 lSeq inDelBitBV) 
        else
            let newCost = prevCost + indelCost
                newState = getUnionIntersectionStateBV inDelBitBV (lSeq V.! (position - 1))
            in
            --trace ("FRC " ++ show newCost)
            if (newState /= inDelBitBV) then --if there was no inDel overlap between states
                LS.cons (newCost, newState, LeftDir) 
                    (getFirstRow  indelCost rowLength (position + 1) newCost lSeq inDelBitBV)
            else                           --indel in both states so no cost
                LS.cons (prevCost, newState, LeftDir)
                    (getFirstRow  indelCost rowLength (position + 1) prevCost lSeq inDelBitBV)

-- | getRow starts at second row (=1) and cretes each row in turn
getRows :: V.Vector BV.BitVector  -> V.Vector BV.BitVector  -> Int -> Int -> Int -> LS.Seq (Int, BV.BitVector, Direction) -> BV.BitVector -> LS.Seq (LS.Seq (Int, BV.BitVector, Direction))
getRows lSeq rSeq indelCost subCost rowNum prevRow inDelBitBV =
    if rowNum == ((V.length rSeq) + 1) then LS.empty
    else 
        let thisRow =  getThisRow lSeq rSeq indelCost subCost rowNum prevRow 0 (V.length lSeq) 0 inDelBitBV
        in
        --trace ("Row " ++ show rowNum)
        --LS.snocFlip thisRow (getRows lSeq rSeq indelCost subCost (rowNum + 1) thisRow inDelBitBV) 
        LS.cons thisRow (getRows lSeq rSeq indelCost subCost (rowNum + 1) thisRow inDelBitBV) 

-- | getThisRow takes sequences and parameters with row number and make a non-first
--row
getThisRow :: V.Vector BV.BitVector  -> V.Vector BV.BitVector  -> Int -> Int -> Int ->  LS.Seq (Int, BV.BitVector, Direction) -> Int -> Int -> Int -> BV.BitVector -> LS.Seq (Int, BV.BitVector, Direction)
getThisRow lSeq rSeq indelCost subCost rowNum prevRow position rowLength prevCost inDelBitBV =
    if position == (rowLength + 1) then LS.empty
    else if position == 0 then
        let newState = getUnionIntersectionStateBV inDelBitBV (rSeq V.! (rowNum - 1))
            (upValue, _, _) = prevRow LS.! position
        in --following in case overlap of inDelBit in leading gaps
        if (newState /= inDelBitBV) then
            LS.cons (upValue + indelCost, newState, DownDir) 
                (getThisRow lSeq rSeq indelCost subCost rowNum prevRow (position + 1) rowLength (upValue + indelCost) inDelBitBV)
        else 
            LS.cons (upValue, newState, DownDir)
                (getThisRow lSeq rSeq indelCost subCost rowNum prevRow (position + 1) rowLength upValue inDelBitBV)
    else 
        let lSeqPos = position - 1 --since first is '-' the index is row/pos - 1
            rSeqRow = rowNum - 1 --since first is '-' the index is row/pos - 1
            leftCost = getOverlapCostBV prevCost indelCost (lSeq V.! lSeqPos) inDelBitBV --need to check for overlap
            (upValue, _, _) = prevRow LS.! position 
            downCost = getOverlapCostBV upValue indelCost (rSeq V.! rSeqRow) inDelBitBV--need to check for overlap
            (diagValue, _, _) = prevRow LS.! (position - 1)
            --intersection = BV.and [(lSeq V.! lSeqPos), (rSeq V.! rSeqRow)]
            --unionLocal = BV.or [(lSeq V.! lSeqPos), (rSeq V.! rSeqRow)]
            intersection = (lSeq V.! lSeqPos) .&. (rSeq V.! rSeqRow)
            unionLocal = (lSeq V.! lSeqPos) .|. (rSeq V.! rSeqRow)
            (diagCost, diagState) = getDiagDirCostBV diagValue intersection unionLocal subCost
            (minCost, minState, minDir) = getMinCostDirBV leftCost downCost diagCost diagState 
                (getUnionIntersectionStateBV inDelBitBV (lSeq V.! lSeqPos)) (getUnionIntersectionStateBV inDelBitBV (rSeq V.! rSeqRow)) 
        in
        --trace ("preRow " ++ show prevRow ++ "row " ++ show rowNum ++ " col " ++ show position 
        --    ++ " trip " ++ show minCost ++ " " ++ show minState ++ " " ++ show minDir)
        LS.cons (minCost, minState, minDir) (getThisRow lSeq rSeq indelCost subCost rowNum prevRow (position + 1) rowLength minCost inDelBitBV)        


-- | getDiagDirCostBV takes union intersection and state to get diagonla sub or no-sub
--cost
getDiagDirCostBV :: Int -> BV.BitVector -> BV.BitVector -> Int -> (Int, BV.BitVector)
getDiagDirCostBV upLeftDirCost intersection unionLocal subCost =
    --trace ("DiagCost " ++ show upLeftDirCost ++ " int " ++ show intersection ++ " union " ++ show union) (
    if not $ BV.isZeroVector intersection then (upLeftDirCost, intersection)
    else (upLeftDirCost + subCost, unionLocal)

-- | getUnionIntersectionBV
getUnionIntersectionStateBV :: BV.BitVector -> BV.BitVector -> BV.BitVector
getUnionIntersectionStateBV lState rState =
    --let intersection = BV.and [lState, rState]
    let intersection = lState .&. rState
    in
    if not $ BV.isZeroVector intersection then intersection
    --else BV.or [lState, rState]
    else lState .|. rState


-- | getMinCostDirBV takes costs and states of three directins and returns min cost,
--directin, and state
--ORDER diag, down, left so same as POY4-5.
getMinCostDirBV :: Int -> Int -> Int -> BV.BitVector -> BV.BitVector -> BV.BitVector -> (Int, BV.BitVector, Direction)
getMinCostDirBV leftCost downCost diagCost diagState leftState downState =
    let minValue = minimum [leftCost, downCost, diagCost]
    in
    --trace ("costs " ++ show leftCost ++ " " ++ show downCost ++ " " ++ show diagCost ++ " -> " ++ show minValue) (
    if diagCost == minValue then (diagCost, diagState, DiagDir)
    else if downCost == minValue then (downCost, downState, DownDir)
    else (leftCost, leftState, LeftDir)
    
-- | getOverlapCostBV cheks for ovelap in gap so if indel, but opossite a gap
--ambiguity--there is no cost
getOverlapCostBV :: Int -> Int -> BV.BitVector -> BV.BitVector -> Int
getOverlapCostBV preCost indelCost oppositeState inDelBitBV =
    --trace("bits " ++ show oppositeState ++ " overAND " ++ show ((.&.) oppositeState inDelBit) ++ " control " ++ show ((.&.) inDelBit inDelBit)) ( 
    --if preCost == barrierCost then barrierCost 
    --else 
    --if BV.and [oppositeState, inDelBitBV] == 0 then preCost + indelCost
    if BV.isZeroVector (oppositeState .&. inDelBitBV) then preCost + indelCost
    else preCost
    --)

