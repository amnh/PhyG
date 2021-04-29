{- |
Module      :  NaiveDOSequence 
Description :  Functions for parsimony DO Naive version using bitvectors (so large alphabet) 
               without Ukkonen space/time saving 
               Uses "Sequence" type to avoid O(n) of Vector.cons
               Prototype non-optimized, restricted Haskell functisns
Copyright   :  (c) 2014 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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
import Data.Int
import Data.Bits
import qualified Data.Vector as V
import Types
import qualified Data.BitVector as BV
import qualified LocalSequence as LS
import GeneralUtilities

data Direction = LeftDir | DownDir | DiagDir
    deriving (Read, Show, Eq)


-- | NaiveDo wraps around naive_do_BV 
-- unlimited alphabet size no space saving
naiveDO :: V.Vector BV.BV -> V.Vector BV.BV -> Int -> (V.Vector BV.BV, Int)
naiveDO lBV rBV inDelCost =
    -- missing data inputs
    if V.null lBV then (rBV, 0)
    else if V.null rBV then (lBV, 0)
    else 
        -- not missing
        -- get inDelCost 
        let bvLength = BV.size (V.head lBV) 
            --setting left most bit to 1 same purpose as inDelBit for Ukkonen
            bv1 = BV.bitVec bvLength (1 :: Integer)
            inDelBitBV = bv1 BV.<<.(BV.bitVec bvLength (bvLength - 1))
            (newMedianLarge, medianCostLarge) = naive_do_BV lBV rBV inDelCost inDelBitBV
        in
        --trace ("DO: " ++ (show inDelCost) ++ " " ++ (show $ V.head $ V.last thisMatrix)) (
        (newMedianLarge, medianCostLarge) 
        --)
        
-- | setLeftRight returns sequence that is longer first,
--shorter second
setLeftRight :: (Ord a) => V.Vector a -> V.Vector a  -> (V.Vector a, Int, V.Vector a, Int)
setLeftRight inL inR = 
        if V.length inL < V.length inR then (inR, V.length inR, inL, V.length inL)
        else if V.length inL > V.length inR then (inL, V.length inL, inR, V.length inR)
        else 
            let outL = max inL inR
                outR = min inL inR 
            in
            (outL, V.length outL, outR, V.length outR)


-- | naive_do takes two input sequences and returns median sequence and cost 
--based on charInfo-1:1 for now
--to do:
--      different costs
--      left/right based on name
--      Ukkonnen
--      C via FFI
--      Affine
naive_do_BV :: V.Vector BV.BV  -> V.Vector BV.BV  -> Int -> BV.BV  ->  (V.Vector BV.BV , Int)
naive_do_BV inlSeq inrSeq inDelCost inDelBitBV =
    if V.null inlSeq then (inrSeq, 0)
    else if V.null inrSeq then (inlSeq, 0)
    else 
        let subCost = 1
            --this for left right constant--want longer in left for Ukkonnen
            (lSeq, lLength, rSeq, rLength) = setLeftRight inlSeq inrSeq
            --lSeq = max inlSeq inrSeq
            --rSeq = min inlSeq inrSeq
            --lLength = V.length lSeq
            --rLength = V.length rSeq
            nwMatrix = LS.cons firstRow (getRows lSeq rSeq inDelCost subCost 1 firstRow inDelBitBV)
            firstRow = getFirstRow inDelCost lLength 0 0 lSeq inDelBitBV
            (cost, _, _) = (nwMatrix LS.! rLength) LS.! lLength
            median = LS.reverse (traceback nwMatrix (V.length rSeq) (V.length lSeq) inDelBitBV)
            --revNWMatrix = LS.reverse nwMatrix
            --nwMatrix = LS.snocFlip firstRow (getRows lSeq rSeq inDelCost subCost 1 firstRow inDelBitBV)
            --(cost, _, _) = (nwMatrix LS.! 0) LS.! lLength
            --median = tracebackReverse nwMatrix 0  (V.length lSeq) inDelBitBV (V.length rSeq) (V.length lSeq)
            
        in
        --trace ("NW: " ++ show cost ++ "\n" ++ (show $ fmap (fmap fst3) nwMatrix))
        (LS.toVector median, cost)
  
-- | tracebackReverse uses revese matrix (ie last elemnt 0,0) 
tracebackReverse :: LS.Seq (LS.Seq (Int, BV.BV, Direction)) -> Int -> Int -> BV.BV -> Int -> Int -> LS.Seq BV.BV
tracebackReverse nwMatrix posL posR inDelBitBV maxL maxR =
        --trace ((show (maxL, maxR)) ++ " psLR " ++ (show (posL,posR)) ++ " " ++ (show ((nwMatrix LS.! posL) LS.! posR))) (
        if posL == maxL && posR == 0 then LS.empty
        else 
            let (_, state, direction) = (nwMatrix LS.! 0 ) LS.! posR
            in
                if (state /= inDelBitBV) then
                    if direction == LeftDir then LS.cons  state (tracebackReverse nwMatrix (posL) (posR - 1) inDelBitBV  maxL maxR ) 
                    else if direction == DownDir then LS.cons state (tracebackReverse (LS.tail nwMatrix) (posL + 1) (posR) inDelBitBV  maxL maxR )  
                    else LS.cons state (tracebackReverse (LS.tail nwMatrix) (posL + 1) (posR - 1) inDelBitBV  maxL maxR )
                else 
                    if direction == LeftDir then (tracebackReverse nwMatrix (posL) (posR - 1) inDelBitBV  maxL maxR ) 
                    else if direction == DownDir then (tracebackReverse (LS.tail nwMatrix) (posL + 1) (posR) inDelBitBV  maxL maxR )  
                    else (tracebackReverse (LS.tail nwMatrix) (posL + 1) (posR - 1) inDelBitBV  maxL maxR ) 
                    

-- | traceback creates REVERSE mediian from nwMatrix, reverse to make tail
--recusive
traceback :: LS.Seq (LS.Seq (Int, BV.BV, Direction)) -> Int -> Int -> BV.BV -> LS.Seq BV.BV
traceback nwMatrix posL posR inDelBitBV =
        --trace ("psLR " ++ show posL ++ " " ++ show posR) (
        if posL == 0 && posR == 0 then LS.empty
        else 
            let (_, state, direction) = (nwMatrix LS.! posL ) LS.! posR
            in
                if (state /= inDelBitBV) then
                    if direction == LeftDir then LS.cons  state (traceback nwMatrix (posL) (posR - 1) inDelBitBV) 
                    else if direction == DownDir then LS.cons state (traceback nwMatrix (posL - 1) (posR) inDelBitBV)  
                    else LS.cons state (traceback nwMatrix (posL - 1) (posR - 1) inDelBitBV)
                else 
                    if direction == LeftDir then (traceback nwMatrix (posL) (posR - 1) inDelBitBV) 
                    else if direction == DownDir then (traceback nwMatrix (posL - 1) (posR) inDelBitBV)  
                    else (traceback nwMatrix (posL - 1) (posR - 1) inDelBitBV) 

-- | getFirstRow initializes foirst row of NW matrix
getFirstRow :: Int -> Int -> Int -> Int -> V.Vector BV.BV  -> BV.BV -> LS.Seq (Int, BV.BV, Direction)
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
getRows :: V.Vector BV.BV  -> V.Vector BV.BV  -> Int -> Int -> Int -> LS.Seq (Int, BV.BV, Direction) -> BV.BV -> LS.Seq (LS.Seq (Int, BV.BV, Direction))
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
getThisRow :: V.Vector BV.BV  -> V.Vector BV.BV  -> Int -> Int -> Int ->  LS.Seq (Int, BV.BV, Direction) -> Int -> Int -> Int -> BV.BV -> LS.Seq (Int, BV.BV, Direction)
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
            intersection = BV.and [(lSeq V.! lSeqPos), (rSeq V.! rSeqRow)]
            unionLocal = BV.or [(lSeq V.! lSeqPos), (rSeq V.! rSeqRow)]
            (diagCost, diagState) = getDiagDirCostBV diagValue intersection unionLocal subCost
            (minCost, minState, minDir) = getMinCostDirBV leftCost downCost diagCost diagState 
                (getUnionIntersectionStateBV inDelBitBV (lSeq V.! lSeqPos)) (getUnionIntersectionStateBV inDelBitBV (rSeq V.! rSeqRow)) 
        in
        --trace ("preRow " ++ show prevRow ++ "row " ++ show rowNum ++ " col " ++ show position 
        --    ++ " trip " ++ show minCost ++ " " ++ show minState ++ " " ++ show minDir)
        LS.cons (minCost, minState, minDir) (getThisRow lSeq rSeq indelCost subCost rowNum prevRow (position + 1) rowLength minCost inDelBitBV)        


-- | getDiagDirCostBV takes union intersection and state to get diagonla sub or no-sub
--cost
getDiagDirCostBV :: Int -> BV.BV -> BV.BV -> Int -> (Int, BV.BV)
getDiagDirCostBV upLeftDirCost intersection unionLocal subCost =
    --trace ("DiagCost " ++ show upLeftDirCost ++ " int " ++ show intersection ++ " union " ++ show union) (
    if intersection /= 0 then (upLeftDirCost, intersection)
    else (upLeftDirCost + subCost, unionLocal)

-- | getUnionIntersectionBV
getUnionIntersectionStateBV :: BV.BV -> BV.BV -> BV.BV
getUnionIntersectionStateBV lState rState =
    let intersection = BV.and [lState, rState]
    in
    if intersection /= 0 then intersection
    else BV.or [lState, rState]


-- | getMinCostDirBV takes costs and states of three directins and returns min cost,
--directin, and state
--ORDER diag, down, left so same as POY4-5.
getMinCostDirBV :: Int -> Int -> Int -> BV.BV -> BV.BV -> BV.BV -> (Int, BV.BV, Direction)
getMinCostDirBV leftCost downCost diagCost diagState leftState downState =
    let minValue = minimum [leftCost, downCost, diagCost]
    in
    --trace ("costs " ++ show leftCost ++ " " ++ show downCost ++ " " ++ show diagCost ++ " -> " ++ show minValue) (
    if diagCost == minValue then (diagCost, diagState, DiagDir)
    else if downCost == minValue then (downCost, downState, DownDir)
    else (leftCost, leftState, LeftDir)
    
-- | getOverlapCostBV cheks for ovelap in gap so if indel, but opossite a gap
--ambiguity--there is no cost
getOverlapCostBV :: Int -> Int -> BV.BV -> BV.BV -> Int
getOverlapCostBV preCost indelCost oppositeState inDelBitBV =
    --trace("bits " ++ show oppositeState ++ " overAND " ++ show ((.&.) oppositeState inDelBit) ++ " control " ++ show ((.&.) inDelBit inDelBit)) ( 
    --if preCost == barrierCost then barrierCost 
    --else 
    if BV.and [oppositeState, inDelBitBV] == 0 then preCost + indelCost
    else preCost
    --)

