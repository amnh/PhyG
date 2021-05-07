{- |
Module      :  DOUkkonnenSequence 
Description :  Functions for parsimony DO optimization functions. using Sequence type
               mainly to reduce vector cons cost
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

{-Improvements
Lots of cons  O(n) stuff--could be improved
-}

-- {-# LANGUAGE DeriveGeneric, DerivingVia, UndecidableInstances #-}

module DOUkkonnenSequence
( ukkonenDO
) where

import Debug.Trace
import Data.Bits
--import qualified Data.BitVector as BV
import qualified Data.BitVector.LittleEndian as BV
import qualified Data.Vector  as V
import qualified LocalSequence as LS
import Data.Bits ((.&.), (.|.))


data Direction = LeftDir | DownDir | DiagDir
  --deriving Eq 
  --deriving Enum
  --deriving V.Unboxable via V.Enum Direction
  deriving (Read, Show, Eq)


type NWElement = (Int, BV.BitVector, Direction)

type BaseChar = V.Vector BV.BitVector 



{-
--Old stuff later bitV.Vector s
inDelBitBV :: BV.BitVector
inDelBitBV = (bit 63) :: BV.BitVector --(bit 63) :: Int --set inDelBitBV to 64th bit in Int
-}

-- | barrierCost really should ne Infnity or somehting like that
barrierCost:: Int
barrierCost = (bit 60) :: Int --really big Int--asssumes 64 bit at least, but can be added to without rolling over.

barrierBit :: BV.BitVector
barrierBit = BV.fromBits [False] --(bit 63) :: BV.BitVector

--FOR both DO's  lseq is a row, acrosss so num columns = length of lseq
--There are rseq rows
-- | UkkonenDO takes two input sequences and returns median sequence and cost
--only 1:1 for now. Uses Ukkonen's space/time saving algorithm
--need to make sure Left/Right and diag/ins/del orders consistent and with
--POY4/5
--lseq > rseq appeard more efficient--could be wrong
--move to C via FFI
--Still occasional error in cost and median (disagreement) show in Chel.seq
ukkonenDO :: BaseChar -> BaseChar -> Int -> (BaseChar, Int)
ukkonenDO lSeq rSeq inDelCost =
    if V.null lSeq then (rSeq, 0)
    else if V.null rSeq then (lSeq, 0)
    else 
        let subCost = 1
            --this for left right constant--want longer in left for Ukkonnen
            lLength = V.length lSeq
            rLength = V.length rSeq

            maxGap = 1 + lLength - rLength  --10000 :: Int --holder lseq - rSeq + 1
            bvLength = BV.dimension (V.head lSeq) 
            --setting left most bit to 1 same purpose as inDelBitBV for Ukkonen
            --bv1 = BV.bitVec bvLength (1 :: Integer)
            --inDelBitBV = bv1 BV.<<.(BV.bitVec bvLength (bvLength - 1))
            bv1 = BV.fromBits [True]
            inDelBitBV = shiftR bv1 (fromIntegral $ bvLength - 1)
            (median, cost) = ukkonenCore lSeq lLength rSeq rLength maxGap inDelCost subCost inDelBitBV
            --firstRow = getFirstRowUkkonen indelCost lLength 0 0 lSeq maxGap
            --nwMatrix = V.cons firstRow (getRowsUkkonen lSeq rSeq indelCost subCost 1 firstRow maxGap)
            --(cost, _, _) = (nwMatrix V.! rLength) V.! (transformFullYShortY lLength rLength  maxGap) --fix for offset
            --median = trace ("getting median " ++ show cost) V.reverse (fromMaybe V.empty (tracebackUkkonen nwMatrix rLength lLength maxGap 0 0))
        in
        --trace ("\nCost/median " ++ show cost ++ "->" ++ show median)
        --(median, fromIntegral cost)
        (LS.toVector $ median, cost)


-- | transformFullYShortY take full Y value (if did entire NW matrix) and returns
--short (Ukkonnen Y) given Y, Y length and row numbera
--remove error when working--overhead
transformFullYShortY :: Int -> Int -> Int -> Int
transformFullYShortY currentY rowNumber maxGap =
    let transformY = currentY - (max 0 (rowNumber - maxGap - 1))
    in
    if (transformY < 0) then error (show currentY ++ " " ++ show rowNumber ++ " " ++ show maxGap ++ " Impossible negative value for transfomred Y")
    else transformY

-- | ukkonenCore core functions of Ukkonen to allow for recursing with maxGap
--doubled if not large enough (returns Nothing)  
ukkonenCore :: BaseChar -> Int -> BaseChar -> Int -> Int -> Int -> Int -> BV.BitVector -> (LS.Seq BV.BitVector, Int)
ukkonenCore lSeq lLength rSeq rLength maxGap indelCost subCost inDelBitBV =  
    let firstRow = getFirstRowUkkonen indelCost lLength 0 0 lSeq maxGap inDelBitBV 
        nwMatrix = LS.cons firstRow (getRowsUkkonen lSeq rSeq indelCost subCost 1 firstRow maxGap inDelBitBV )
        (cost, _, _) = LS.last (LS.last nwMatrix) -- V.! rLength) --V.! (transformFullYShortY lLength rLength  maxGap) --fix for offset
        medianTest = tracebackUkkonen nwMatrix rLength lLength maxGap 0 0  inDelBitBV 
    in
    --trace ("mh " ++ show (V.head medianTest)) (
    if not $ BV.isZeroVector (LS.last medianTest) then -- (0 :: Int) then 
        --trace (show nwMatrix) 
        (medianTest, fromIntegral cost)
    else --trace ("Going back!! " ++ show cost) 
        ukkonenCore lSeq lLength rSeq rLength (2 * maxGap) indelCost subCost  inDelBitBV 
    


-- | tracebackUkkonen creates REVERSE mediian from nwMatrix, reverse to make tail
--recusive, for Ukkonen space/time saving offsets
--need to count gaps in traceback for threshold/barrier stuff
--CHANGE TO MAYBE (V.Vector  Int) FOR BARRIER CHECK
tracebackUkkonen :: LS.Seq  (LS.Seq NWElement) -> Int -> Int -> Int -> Int -> Int -> BV.BitVector -> LS.Seq  BV.BitVector
tracebackUkkonen nwMatrix posR posL maxGap rInDel lInDel  inDelBitBV =
        --trace ("psLR " ++ show posR ++ " " ++ show posL ++ " Left " ++ show lInDel ++ " Right " ++ show rInDel ++ " maxGap " ++ show maxGap) (
        if (rInDel  > (maxGap - 2)) || (lInDel > (maxGap - 2)) then LS.singleton $ BV.fromBits [False] --
        else if posL == 0 && posR == 0 then LS.empty
        else 
            --trace ("TB " ++ show posL ++ " " ++ show posR) (
            let (_, state, direction) = (nwMatrix LS.! posR) LS.! (transformFullYShortY posL posR  maxGap) --(transformFullYShortY posL posR maxGap)
            in
                --trace ("state " ++ show state ++ " dir " ++ show direction) (
                if (state /= inDelBitBV) then
                    if direction == LeftDir then LS.snocFlip  state (tracebackUkkonen nwMatrix (posR) (posL - 1) maxGap rInDel (lInDel + 1) inDelBitBV) 
                    else if direction == DownDir then LS.snocFlip state (tracebackUkkonen nwMatrix (posR - 1) (posL) maxGap (rInDel + 1) lInDel inDelBitBV)  
                    else LS.snocFlip state (tracebackUkkonen nwMatrix (posR - 1) (posL - 1) maxGap rInDel lInDel inDelBitBV)
                else 
                    if direction == LeftDir then tracebackUkkonen nwMatrix (posR) (posL - 1) maxGap rInDel (lInDel + 1) inDelBitBV 
                    else if direction == DownDir then tracebackUkkonen nwMatrix (posR - 1) (posL) maxGap (rInDel + 1) lInDel inDelBitBV  
                    else tracebackUkkonen nwMatrix (posR - 1) (posL - 1) maxGap rInDel lInDel inDelBitBV 
            --)--)

-- | getFirstRowUkkonen initializes first row of NW-Ukkonen matrix
getFirstRowUkkonen :: Int -> Int -> Int -> Int -> BaseChar -> Int -> BV.BitVector -> LS.Seq NWElement
getFirstRowUkkonen indelCost rowLength position prevCost lSeq  maxGap inDelBitBV = 
    --trace ("row 0 pos " ++ show position ++ "/" ++ show (maxShortY rowLength 0 maxGap) ++ " rowLength " ++ show rowLength ++ " maxGap " ++ show maxGap ++ " lseq " ++ show lSeq) (
    if position == rowLength  + 1 then LS.empty
    else if position == (maxGap + 1) then LS.singleton (barrierCost, barrierBit, LeftDir) 
    else
        if position == 0 then LS.cons (0, inDelBitBV, DiagDir) (getFirstRowUkkonen indelCost rowLength (position + 1) 0 lSeq maxGap inDelBitBV) 
        else
            let newCost = prevCost + indelCost
                newState = getUnionIntersectionState inDelBitBV (lSeq V.! (position - 1))
            in
            --trace ("FRC " ++ show newCost)
            if (newState /= inDelBitBV) then --if there was no inDel overlap between states
                LS.cons (newCost, newState, LeftDir) 
                    (getFirstRowUkkonen  indelCost rowLength (position + 1) newCost lSeq maxGap inDelBitBV)
            else                           --indel in both states so no cost
                LS.cons (prevCost, newState, LeftDir)
                    (getFirstRowUkkonen  indelCost rowLength (position + 1) prevCost lSeq maxGap inDelBitBV)
   --) 

-- | getRowUkkonen starts at second row (=1) and creates each row in turn--Ukkonen
getRowsUkkonen :: BaseChar -> BaseChar -> Int -> Int -> Int -> LS.Seq NWElement -> Int -> BV.BitVector -> LS.Seq  (LS.Seq NWElement)
getRowsUkkonen lSeq rSeq indelCost subCost rowNum prevRow maxGap  inDelBitBV=
    --trace ("Row " ++ show rowNum) (
    if rowNum == ((V.length rSeq) + 1) then LS.empty
    else 
        let startPosition = max 0 (rowNum - maxGap) --check for left barriers 
            thisRowZero =  getThisRowUkkonen lSeq rSeq indelCost subCost rowNum prevRow startPosition (V.length lSeq) 0 maxGap inDelBitBV
            thisRowNonZero = LS.cons (barrierCost, barrierBit, DownDir) (getThisRowUkkonen lSeq rSeq indelCost subCost rowNum prevRow startPosition  (V.length lSeq) barrierCost maxGap  inDelBitBV)
        in
        if startPosition == 0 then 
            --trace ("Row " ++ show rowNum ++ " of " ++ show (V.length rSeq) ++ " starts " ++ show startPosition ++ ":" ++ show thisRowZero) (
                LS.cons thisRowZero (getRowsUkkonen lSeq rSeq indelCost subCost (rowNum + 1) thisRowZero maxGap inDelBitBV) 
                --)
        else
            --trace ("Row " ++ show rowNum ++ " of " ++ show (V.length rSeq) ++" starts " ++ show startPosition ++ ":" ++ show thisRowNonZero) (
            LS.cons thisRowNonZero (getRowsUkkonen lSeq rSeq indelCost subCost (rowNum + 1) thisRowNonZero maxGap inDelBitBV)
            --)

-- | getThisRowUkkonen takes sequences and parameters with row number and make a non-first
--row--Ukkonen
getThisRowUkkonen :: BaseChar -> BaseChar -> Int -> Int -> Int ->  LS.Seq NWElement -> Int -> Int -> Int -> Int -> BV.BitVector -> LS.Seq NWElement
getThisRowUkkonen lSeq rSeq indelCost subCost rowNum prevRow position rowLength prevCost maxGap inDelBitBV=
    --trace ("Column " ++ show position) (
    if  position ==  rowLength  + 1 then LS.empty
    else if position == (rowNum + maxGap + 1) then LS.singleton (barrierCost, barrierBit, LeftDir)
    else if position == 0 then
        let newState = getUnionIntersectionState inDelBitBV (rSeq V.! (rowNum - 1))
            (upValue, _, _) = prevRow LS.! position 
        in --following in case overlap of inDelBitBV in leading gaps
        if (newState /= inDelBitBV) then
            LS.cons (upValue + indelCost, newState, DownDir) 
                (getThisRowUkkonen lSeq rSeq indelCost subCost rowNum prevRow (position + 1) rowLength (upValue + indelCost) maxGap inDelBitBV)
        else 
            LS.cons (upValue, newState, DownDir)
                (getThisRowUkkonen lSeq rSeq indelCost subCost rowNum prevRow (position + 1) rowLength upValue maxGap inDelBitBV)
    else 
        let lSeqPos = position - 1 --since first is '-' the index is row/pos - 1
            rSeqRow = rowNum - 1 --since first is '-' the index is row/pos - 1
            leftCost = getOverlapCost prevCost indelCost (lSeq V.! lSeqPos) inDelBitBV --need to check for overlap
            (upValue, _, _) = prevRow LS.! (transformFullYShortY  position (rowNum - 1) maxGap)
            downCost = getOverlapCost upValue indelCost (rSeq V.! rSeqRow)  inDelBitBV --need to check for overlap
            (diagValue, _, _) = prevRow LS.! (transformFullYShortY  (position - 1) (rowNum - 1) maxGap)
            --intersection = BV.and [(lSeq V.! lSeqPos) , (rSeq V.! rSeqRow)]
            --unionLocal = BV.or [(lSeq V.! lSeqPos) , (rSeq V.! rSeqRow)]
            intersection = (lSeq V.! lSeqPos) .&. (rSeq V.! rSeqRow)
            unionLocal = (lSeq V.! lSeqPos) .|. (rSeq V.! rSeqRow)
            (diagCost, diagState) = getDiagDirCost diagValue intersection unionLocal subCost
            (minCost, minState, minDir) = getMinCostDir leftCost downCost diagCost diagState 
                (getUnionIntersectionState inDelBitBV (lSeq V.! lSeqPos)) (getUnionIntersectionState inDelBitBV (rSeq V.! rSeqRow)) 
        in
        --trace ("preRow " ++ show prevRow ) -- ++ "row " ++ show rowNum ++ " col " ++ show position 
        --    ++ " trip " ++ show minCost ++ " " ++ show minState ++ " " ++ show minDir)
        LS.cons (minCost, minState, minDir) (getThisRowUkkonen lSeq rSeq indelCost subCost rowNum prevRow (position + 1) rowLength minCost maxGap inDelBitBV)        
        --)

-- | getMinCostDir takes costs and states of three directins and returns min cost,
--directin, and state
--ORDER diag, down, left so same as POY4-5.
getMinCostDir :: Int -> Int -> Int -> BV.BitVector -> BV.BitVector -> BV.BitVector -> NWElement
getMinCostDir leftCost downCost diagCost diagState leftState downState =
    let minValue = minimum [leftCost, downCost, diagCost]
    in
    --trace ("costs " ++ show leftCost ++ " " ++ show downCost ++ " " ++ show diagCost ++ " -> " ++ show minValue) (
    if diagCost == minValue then (diagCost, diagState, DiagDir)
    else if downCost == minValue then (downCost, downState, DownDir)
    else (leftCost, leftState, LeftDir)
    --else error "Error in getMinCost"
    {-else (downCost, downState, DownDir)
    else (leftCost, leftState, LeftDir)
    if leftCost == minValue then (leftCost, leftState, LeftDir)
    else if downCost == minValue then (downCost, downState, DownDir)
    else (diagCost, diagState, DiagDir)
    -}

-- | getDiagDirCost takes union intersection and state to get diagonla sub or no-sub
--cost
getDiagDirCost :: Int -> BV.BitVector -> BV.BitVector -> Int -> (Int, BV.BitVector)
getDiagDirCost upLeftDirCost intersection unionLocal subCost =
    --trace ("DiagCost " ++ show upLeftDirCost ++ " int " ++ show intersection ++ " union " ++ show union) (
    if not $ BV.isZeroVector intersection then (upLeftDirCost, intersection)
    else (upLeftDirCost + subCost, unionLocal)

-- | getUnionIntersection
getUnionIntersectionState :: BV.BitVector -> BV.BitVector -> BV.BitVector
getUnionIntersectionState lState rState =
    --let intersection = BV.and [lState , rState] 
    let intersection = lState .&. rState 
    in
    if not $ BV.isZeroVector intersection then intersection
    --else BV.or [lState , rState] 
    else lState .|. rState

-- | getOverlapCost cheks for ovelap in gap so if indel, but opossite a gap
--ambiguity--there is no cost
getOverlapCost :: Int -> Int -> BV.BitVector -> BV.BitVector -> Int
getOverlapCost preCost indelCost oppositeState inDelBitBV =
    --trace("bits " ++ show oppositeState ++ " overAND " ++ show ((.&.) oppositeState inDelBitBV) ++ " control " ++ show ((.&.) inDelBitBV inDelBitBV)) ( 
    --if preCost == barrierCost then barrierCost 
    --else 
    if BV.isZeroVector (oppositeState .&. inDelBitBV) then preCost + indelCost
    else preCost
    --)

