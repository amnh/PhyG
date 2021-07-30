{- |
Module      :  DOWrapper
Description :  Controles the various DO functions for differnet situations (alphabet size, tcm etc)
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

{-Improvements
Lots of cons  O(n) stuff--could be improved
-}

-- {-# LANGUAGE DeriveGeneric, DerivingVia, UndecidableInstances #-}

module DirectOptimization.DOWrapper
( getDOMedian
) where

import           Data.Int
import           Debug.Trace
--import qualified DirectOptimization.DOSmallFFI as DOSmallFFI
--import qualified DirectOptimization.DOLargeFFI as DOLargeFFI

import qualified Data.Vector                                    as V
import qualified Utilities.Utilities                            as U
--import qualified Data.BitVector  as BV
import qualified Bio.Character.Encodable.Dynamic.AmbiguityGroup as AG
import           Data.Alphabet
import qualified Data.BitVector.LittleEndian                    as BV
import           Data.Bits                                      (shiftL, shiftR, (.&.), (.|.))
import           Data.Foldable
import           Data.List.NonEmpty                             (NonEmpty(..))
import           Data.MetricRepresentation
import qualified Data.MetricRepresentation                      as MR
import qualified Data.TCM                                       as TCM
import qualified Data.TCM.Dense                                 as TCMD
import qualified Data.TCM.Memoized                              as TCMM
import qualified DirectOptimization.DOUkkonenSequence           as DKS
import qualified DirectOptimization.DOUkkonenSequenceInt64      as DKS64
import qualified DirectOptimization.NaiveDOSequence             as NS
import qualified SymMatrix                                      as SM
import           Types.Types

import qualified DirectOptimization.DOHuge                      as HUGE
import qualified DirectOptimization.DOSlim                      as SLIM
import qualified DirectOptimization.DOWide                      as WIDE



-- | thesholdUKLength sets threshold for where its worth doing Ukkonen stuff
-- short seqeuneces not worth it.  This should be tested empirically
thesholdUKLength :: Int
thesholdUKLength = 15

-- | getDOMedian wraps around getPrelim and getPrelim3
-- changing types and checking for missing cases
getDOMedian :: V.Vector  (V.Vector  Int)
            -> TCMD.DenseTransitionCostMatrix
            -> MR.MetricRepresentation BV.BitVector
            -> CharType -> (V.Vector  BV.BitVector, V.Vector  BV.BitVector) -> (V.Vector  BV.BitVector, Int)
getDOMedian thisMatrix tcmDense tcmMemo thisType (origLBV, origRBV) =
    -- missing data inputs
    if V.null origLBV then (origRBV, 0)
    else if V.null origRBV then (origLBV, 0)
    else
        -- not missing
        -- get inDelCost
        --trace ("2 in :\n" ++ (concat $ V.map (U.bitVectToCharState ["A","C","G","T","-"]) origLBV) ++ "\n" ++ (concat $ V.map (U.bitVectToCharState ["A","C","G","T","-"]) origRBV)) (
        let inDelCost = V.head (V.last thisMatrix)
            (lBV, rBV) = setLeftRight origLBV origRBV
            -- Int64 version faster for small aphabets
            bvLength = BV.dimension (V.head lBV)
            bv1 = BV.fromBits (True : (replicate ((fromIntegral bvLength) -1) False))
            gapCharBit =  shiftL bv1 ((fromIntegral bvLength) - 1)
            leftChar64 = V.map convertBVTo64 lBV
            rightChar64 = V.map convertBVTo64 rBV
            (newMedian64, medianCost64) = DKS64.ukkonenDO leftChar64 rightChar64 inDelCost
            newMedianBV = V.map (convert64ToBV bvLength) newMedian64

            --setting left most bit to 1 same purpose as inDelBit for Ukkonen--I thnik small error in this
            -- after bitvector change but Naive OK
            --(newMedianSmall, medianCostSmall) = DKS.ukkonenDO lBV rBV inDelCost
            --(newMedianLarge, medianCostLarge) = NS.naiveDO lBV rBV inDelCost
            (newMedianSmall, medianCostSmall) = SLIM.wrapperSlimDO lBV rBV thisMatrix
            (newMedianMedium, medianCostMedium) = WIDE.wrapperWideDO lBV rBV tcmMemo
            (newMedianLarge, medianCostLarge) = HUGE.wrapperHugeDO lBV rBV tcmMemo

            --(mediansFFI, costFFI) = DOSmallFFI.wrapperPCG_DO_Small_FFI lBV rBV thisMatrix tcmDense
            --mediansFFI' =  V.filter (/= gapCharBit) mediansFFI

            -- Problems with tcmMemo FFI calls--erratic/inconsistent behavior
            --(mediansLargeFFI, costLargeFFI) = DOLargeFFI.wrapperPCG_DO_Large lBV rBV thisMatrix tcmMemo
        in

        if (thisType == NucSeq || thisType == SlimSeq) then(newMedianSmall, medianCostSmall) -- (newMedianLarge, medianCostLarge) --
        else if (thisType == AminoSeq || thisType == WideSeq) then (newMedianMedium, medianCostMedium)
        else if (thisType == HugeSeq) then  (newMedianLarge, medianCostLarge)
        else error "Unrecognized/Not implemented character type"
        --)


-- | convertBVTo64 converts bitV.Vector  type to bit64 Int type
convertBVTo64 :: BV.BitVector -> Int64
convertBVTo64 inBV = fromIntegral (BV.toUnsignedNumber inBV)

-- | convert64ToBV converts bitV.Vector  type to bit64 Int type
convert64ToBV :: Word -> Int64 -> BV.BitVector
convert64ToBV bvLength in64  =  BV.fromNumber bvLength in64

-- | setLeftRight returns sequence that is longer first,
--shorter second
setLeftRight :: V.Vector  BV.BitVector -> V.Vector  BV.BitVector -> (V.Vector  BV.BitVector, V.Vector  BV.BitVector)
setLeftRight inL inR =
        if V.length inL < V.length inR then (inR, inL)
        else if V.length inL > V.length inR then (inL, inR)
        else
            let outL = max inL inR
                outR = min inL inR
            in
            (outL, outR)
