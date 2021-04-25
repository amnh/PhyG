{- |
Module      :  DOPrototype 
Description :  Functions for parsimony DO optimization functions
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

{-# LANGUAGE DeriveGeneric, DerivingVia, UndecidableInstances #-}

module DOPrototype
( getDOMedian
) where

import Debug.Trace
import Data.Int
import Data.Bits
import qualified Data.Vector as VB
-- import qualified Data.Vector.Unboxed as V
import qualified Data.Vector.Generic as VG (convert)
import Types
import qualified Data.BitVector as BV
import qualified NaiveDO as NDO
import qualified Data.Vector.Unboxing as V --Unboxing
import GHC.Generics


data Direction = LeftDir | DownDir | DiagDir
  deriving Eq 
  deriving Enum
  deriving V.Unboxable via V.Enum Direction
    --deriving (Read, Show, Eq)
    
type NWElement = (Int, Int64, Direction)

type VectorNWElement = V.Vector NWElement


data NWElementStruct = NWElementStruct {fstS :: (Int, Int64, Direction)}
  deriving Generic
  deriving V.Unboxable via V.Generics NWElementStruct

type BaseChar = V.Vector Int64 


--  Wrappers to interface with PhyGraph Types
--      This is a temporary step until PCG/POY code 
--      Imported vi FFIs
--      UkkonenDO is for small alphabets keeping bit64
--      NaiveDO is for large alphabets migrated to general bitvector

{-
-- | convert between boxed and unboxed vector
boxedToUnboxed :: (VU.Unbox a) => V.Vector a -> VU.Vector a
boxedToUnboxed inVect = VG.convert inVect

unbBoxedtoBoxed :: (VU.Unbox a) => VU.Vector a -> V.Vector a
unbBoxedtoBoxed inVect = VG.convert inVect
-}

-- | getDOMedian wraps around getPrelim and getPrelim3
-- changing types and checking for missing cases
getDOMedian :: VB.Vector BV.BV -> VB.Vector BV.BV -> VB.Vector (VB.Vector Int) -> CharType -> (VB.Vector BV.BV, Int)
getDOMedian lBV rBV thisMatrix thisType =
    -- missing data inputs
    if VB.null lBV then (rBV, 0)
    else if VB.null rBV then (lBV, 0)
    else 
        -- not missing
        -- get inDelCost 
        let bvLength = BV.size (VB.head lBV) 
            inDelCost = VB.head (VB.last thisMatrix)
            leftChar64 = VG.convert $ VB.map convertBVTo64 lBV
            rightChar64 = VG.convert $ VB.map convertBVTo64 rBV
            (newMedianSmall, medianCostSmall) = ukkonenDO leftChar64 rightChar64 inDelCost
            newMedianSmallBV = VB.map (convert64ToBV bvLength) $ VG.convert newMedianSmall
            --setting left most bit to 1 same purpose as inDelBit for Ukkonen
            (newMedianLarge, medianCostLarge) = NDO.naiveDO lBV rBV inDelCost
        in
        --trace ("DO: " ++ (s(VectorNWElement)how inDelCost) ++ " " ++ (show $ V.head $ V.last thisMatrix)) (
        if thisType == NucSeq then (newMedianSmallBV, medianCostSmall)
        else if thisType == GenSeq then (newMedianLarge, medianCostLarge) 
        else error "Unrecognized/Not implemented character type"
        --)
        

-- | convertBVTo64 converts bitvector type to bit64 Int64 type 
convertBVTo64 :: BV.BV -> Int64
convertBVTo64 inBV = fromIntegral (BV.nat inBV) 

-- | convert64ToBV converts bitvector type to bit64 Int64 type 
convert64ToBV :: Int -> Int64 -> BV.BV
convert64ToBV bvLength in64  =  BV.bitVec  bvLength in64

--Old stuff later bitvectors
inDelBit :: Int64
inDelBit = (bit 63) :: Int64 --(bit 63) :: Int64 --set indelBit to 64th bit in Int64

barrierCost:: Int
barrierCost = (bit 60) :: Int --really big Int--asssumes 64 bit at least, but can be added to without rolling over.

barrierBit :: Int64
barrierBit = (bit 63) :: Int64

-- | transformFullYShortY take full Y value (if did entire NW matrix) and returns
--short (Ukkonnen Y) given Y, Y length and row numbera
--remove error when working--overhead
transformFullYShortY :: Int -> Int -> Int -> Int
transformFullYShortY currentY rowNumber maxGap =
    let transformY = currentY - (max 0 (rowNumber - maxGap - 1))
    in
    if (transformY < 0) then error (show currentY ++ " " ++ show rowNumber ++ " " ++ show maxGap ++ " Impossible negative value for transfomred Y")
    else transformY

-- | setLeftRight returns sequence that is longer first,
--shorter second
setLeftRight :: BaseChar -> BaseChar  -> (BaseChar, Int, BaseChar, Int)
setLeftRight inL inR = 
        if V.length inL < V.length inR then (inR, V.length inR, inL, V.length inL)
        else if V.length inL > V.length inR then (inL, V.length inL, inR, V.length inR)
        else 
            let outL = max inL inR
                outR = min inL inR 
            in
            (outL, V.length outL, outR, V.length outR)

-- | ukkonenCore core functions of Ukkonen to allow for recursing with maxGap
--doubled if not large enough (returns Nothing)  
ukkonenCore :: BaseChar -> Int -> BaseChar -> Int -> Int -> Int -> Int -> (BaseChar, Int)
ukkonenCore lSeq lLength rSeq rLength maxGap indelCost subCost =  
    let firstRow = getFirstRowUkkonen indelCost lLength 0 0 lSeq maxGap
        nwMatrix = V.cons firstRow (getRowsUkkonen lSeq rSeq indelCost subCost 1 firstRow maxGap)
        (cost, _, _) = V.last (V.last nwMatrix) -- V.! rLength) --V.! (transformFullYShortY lLength rLength  maxGap) --fix for offset
        medianTest = tracebackUkkonen nwMatrix rLength lLength maxGap 0 0
    in
    --trace ("mh " ++ show (V.head medianTest)) (
    if  (V.last medianTest) /= 0 then -- (0 :: Int64) then 
        --trace (show nwMatrix) 
        (V.reverse medianTest, fromIntegral cost)
    else --trace ("Going back!! " ++ show cost) 
        ukkonenCore lSeq lLength rSeq rLength (2 * maxGap) indelCost subCost
    

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
ukkonenDO inlSeq inrSeq inDelCost =
    if V.null inlSeq then (inrSeq, 0)
    else if V.null inrSeq then (inlSeq, 0)
    else 
        let subCost = 1
            --this for left right constant--want longer in left for Ukkonnen
            (lSeq, lLength, rSeq, rLength) = setLeftRight inlSeq inrSeq
            maxGap = 1 + lLength - rLength  --10000 :: Int --holder lseq - rSeq + 1
            (median, cost) = ukkonenCore lSeq lLength rSeq rLength maxGap inDelCost subCost
            --firstRow = getFirstRowUkkonen indelCost lLength 0 0 lSeq maxGap
            --nwMatrix = V.cons firstRow (getRowsUkkonen lSeq rSeq indelCost subCost 1 firstRow maxGap)
            --(cost, _, _) = (nwMatrix V.! rLength) V.! (transformFullYShortY lLength rLength  maxGap) --fix for offset
            --median = trace ("getting median " ++ show cost) V.reverse (fromMaybe V.empty (tracebackUkkonen nwMatrix rLength lLength maxGap 0 0))
        in
        --trace ("\nCost/median " ++ show cost ++ "->" ++ show median)
        --(median, fromIntegral cost)
        (median, cost)

-- | tracebackUkkonen creates REVERSE mediian from nwMatrix, reverse to make tail
--recusive, for Ukkonen space/time saving offsets
--need to count gaps in traceback for threshold/barrier stuff
--CHANGE TO MAYBE (V.Vector Int64) FOR BARRIER CHECK
tracebackUkkonen :: V.Vector (VectorNWElement) -> Int -> Int -> Int -> Int -> Int -> V.Vector Int64
tracebackUkkonen nwMatrix posR posL maxGap rInDel lInDel =
        --trace ("psLR " ++ show posR ++ " " ++ show posL ++ " Left " ++ show lInDel ++ " Right " ++ show rInDel ++ " maxGap " ++ show maxGap) (
        if (rInDel  > (maxGap - 2)) || (lInDel > (maxGap - 2)) then V.singleton 0 --
        else if posL == 0 && posR == 0 then V.empty
        else 
            --trace ("TB " ++ show posL ++ " " ++ show posR) (
            let (_, state, direction) = (nwMatrix V.! posR) V.! (transformFullYShortY posL posR  maxGap) --(transformFullYShortY posL posR maxGap)
            in
                --trace ("state " ++ show state ++ " dir " ++ show direction) (
                if (state /= inDelBit) then
                    if direction == LeftDir then (V.cons  state (tracebackUkkonen nwMatrix (posR) (posL - 1) maxGap rInDel (lInDel + 1))) 
                    else if direction == DownDir then (V.cons state (tracebackUkkonen nwMatrix (posR - 1) (posL) maxGap (rInDel + 1) lInDel))  
                    else (V.cons state (tracebackUkkonen nwMatrix (posR - 1) (posL - 1) maxGap rInDel lInDel))
                else 
                    if direction == LeftDir then (tracebackUkkonen nwMatrix (posR) (posL - 1) maxGap rInDel (lInDel + 1)) 
                    else if direction == DownDir then (tracebackUkkonen nwMatrix (posR - 1) (posL) maxGap (rInDel + 1) lInDel)  
                    else (tracebackUkkonen nwMatrix (posR - 1) (posL - 1) maxGap rInDel lInDel) 
            --)--)

-- | getFirstRowUkkonen initializes first row of NW-Ukkonen matrix
getFirstRowUkkonen :: Int -> Int -> Int -> Int -> BaseChar -> Int -> VectorNWElement
getFirstRowUkkonen indelCost rowLength position prevCost lSeq  maxGap = 
    --trace ("row 0 pos " ++ show position ++ "/" ++ show (maxShortY rowLength 0 maxGap) ++ " rowLength " ++ show rowLength ++ " maxGap " ++ show maxGap ++ " lseq " ++ show lSeq) (
    if position == rowLength  + 1 then V.empty
    else if position == (maxGap + 1) then V.singleton (barrierCost, barrierBit, LeftDir) 
    else
        if position == 0 then V.cons (0, inDelBit, DiagDir) (getFirstRowUkkonen indelCost rowLength (position + 1) 0 lSeq maxGap) 
        else
            let newCost = prevCost + indelCost
                newState = getUnionIntersectionState inDelBit (lSeq V.! (position - 1))
            in
            --trace ("FRC " ++ show newCost)
            if (newState /= inDelBit) then --if there was no inDel overlap between states
                V.cons (newCost, newState, LeftDir) 
                    (getFirstRowUkkonen  indelCost rowLength (position + 1) newCost lSeq maxGap)
            else                           --indel in both states so no cost
                V.cons (prevCost, newState, LeftDir)
                    (getFirstRowUkkonen  indelCost rowLength (position + 1) prevCost lSeq maxGap)
   --) 

-- | getRowUkkonen starts at second row (=1) and creates each row in turn--Ukkonen
getRowsUkkonen :: BaseChar -> BaseChar -> Int -> Int -> Int -> VectorNWElement -> Int -> V.Vector (VectorNWElement)
getRowsUkkonen lSeq rSeq indelCost subCost rowNum prevRow maxGap =
    if rowNum == ((V.length rSeq) + 1) then V.empty
    else 
        let startPosition = max 0 (rowNum - maxGap) --check for left barriers 
            thisRowZero =  getThisRowUkkonen lSeq rSeq indelCost subCost rowNum prevRow startPosition (V.length lSeq) 0 maxGap
            thisRowNonZero = V.cons (barrierCost, barrierBit, DownDir) (getThisRowUkkonen lSeq rSeq indelCost subCost rowNum prevRow startPosition  (V.length lSeq) barrierCost maxGap )
        in
        if startPosition == 0 then 
            --trace ("Row " ++ show rowNum ++ " of " ++ show (V.length rSeq) ++ " starts " ++ show startPosition ++ ":" ++ show thisRowZero) (
                V.cons thisRowZero (getRowsUkkonen lSeq rSeq indelCost subCost (rowNum + 1) thisRowZero maxGap) 
                --)
        else
            --trace ("Row " ++ show rowNum ++ " of " ++ show (V.length rSeq) ++" starts " ++ show startPosition ++ ":" ++ show thisRowNonZero) (
            V.cons thisRowNonZero (getRowsUkkonen lSeq rSeq indelCost subCost (rowNum + 1) thisRowNonZero maxGap)
            --)

-- | getThisRowUkkonen takes sequences and parameters with row number and make a non-first
--row--Ukkonen
getThisRowUkkonen :: BaseChar -> BaseChar -> Int -> Int -> Int ->  VectorNWElement -> Int -> Int -> Int -> Int -> VectorNWElement
getThisRowUkkonen lSeq rSeq indelCost subCost rowNum prevRow position rowLength prevCost maxGap =
    if  position ==  rowLength  + 1 then V.empty
    else if position == (rowNum + maxGap + 1) then V.singleton (barrierCost, barrierBit, LeftDir)
    else if position == 0 then
        let newState = getUnionIntersectionState inDelBit (rSeq V.! (rowNum - 1))
            (upValue, _, _) = prevRow V.! position 
        in --following in case overlap of inDelBit in leading gaps
        if (newState /= inDelBit) then
            V.cons (upValue + indelCost, newState, DownDir) 
                (getThisRowUkkonen lSeq rSeq indelCost subCost rowNum prevRow (position + 1) rowLength (upValue + indelCost) maxGap)
        else 
            V.cons (upValue, newState, DownDir)
                (getThisRowUkkonen lSeq rSeq indelCost subCost rowNum prevRow (position + 1) rowLength upValue maxGap)
    else 
        let lSeqPos = position - 1 --since first is '-' the index is row/pos - 1
            rSeqRow = rowNum - 1 --since first is '-' the index is row/pos - 1
            leftCost = getOverlapCost prevCost indelCost (lSeq V.! lSeqPos) --need to check for overlap
            (upValue, _, _) = prevRow V.! (transformFullYShortY  position (rowNum - 1) maxGap)
            downCost = getOverlapCost upValue indelCost (rSeq V.! rSeqRow) --need to check for overlap
            (diagValue, _, _) = prevRow V.! (transformFullYShortY  (position - 1) (rowNum - 1) maxGap)
            intersection = (lSeq V.! lSeqPos) .&. (rSeq V.! rSeqRow)
            unionLocal = (lSeq V.! lSeqPos) .|. (rSeq V.! rSeqRow)
            (diagCost, diagState) = getDiagDirCost diagValue intersection unionLocal subCost
            (minCost, minState, minDir) = getMinCostDir leftCost downCost diagCost diagState 
                (getUnionIntersectionState inDelBit (lSeq V.! lSeqPos)) (getUnionIntersectionState inDelBit (rSeq V.! rSeqRow)) 
        in
        --trace ("preRow " ++ show prevRow ) -- ++ "row " ++ show rowNum ++ " col " ++ show position 
        --    ++ " trip " ++ show minCost ++ " " ++ show minState ++ " " ++ show minDir)
        V.cons (minCost, minState, minDir) (getThisRowUkkonen lSeq rSeq indelCost subCost rowNum prevRow (position + 1) rowLength minCost maxGap)        

-- | getMinCostDir takes costs and states of three directins and returns min cost,
--directin, and state
--ORDER diag, down, left so same as POY4-5.
getMinCostDir :: Int -> Int -> Int -> Int64 -> Int64 -> Int64 -> NWElement
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
getDiagDirCost :: Int -> Int64 -> Int64 -> Int -> (Int, Int64)
getDiagDirCost upLeftDirCost intersection unionLocal subCost =
    --trace ("DiagCost " ++ show upLeftDirCost ++ " int " ++ show intersection ++ " union " ++ show union) (
    if intersection /= 0 then (upLeftDirCost, intersection)
    else (upLeftDirCost + subCost, unionLocal)

-- | getUnionIntersection
getUnionIntersectionState :: Int64 -> Int64 -> Int64
getUnionIntersectionState lState rState =
    let intersection = (lState .&. rState) :: Int64
    in
    if intersection /= 0 then intersection
    else (lState .|. rState) :: Int64

-- | getOverlapCost cheks for ovelap in gap so if indel, but opossite a gap
--ambiguity--there is no cost
getOverlapCost :: Int -> Int -> Int64 -> Int
getOverlapCost preCost indelCost oppositeState =
    --trace("bits " ++ show oppositeState ++ " overAND " ++ show ((.&.) oppositeState inDelBit) ++ " control " ++ show ((.&.) inDelBit inDelBit)) ( 
    --if preCost == barrierCost then barrierCost 
    --else 
    if (.&.) oppositeState inDelBit == (0 :: Int64) then preCost + indelCost
    else preCost
    --)

