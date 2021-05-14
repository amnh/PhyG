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

module DOWrapper
( getDOMedian
) where

import Debug.Trace
import Data.Int
import qualified DOSmallFFI as DOSmallFFI
import qualified DOLargeFFI as DOLargeFFI
import qualified Data.Vector  as V
--import qualified Data.BitVector  as BV
import qualified Data.BitVector.LittleEndian as BV
import qualified NaiveDOSequence as NS
import qualified DOUkkonnenSequence as DKS
import qualified DOUkkonnenSequenceInt64 as DKS64
import qualified Data.TCM.Dense as TCMD
import qualified SymMatrix as SM
import qualified Data.TCM as TCM
import qualified Data.MetricRepresentation as MR
import qualified Bio.Character.Encodable.Dynamic.AmbiguityGroup as AG
import qualified Data.TCM.Memoized as TCMM
import Data.Alphabet
import Data.Foldable
import Data.MetricRepresentation
import Data.List.NonEmpty (NonEmpty(..))


import Types

-- | thesholdUKLength sets threshold for where its worth doing Ukkonen stuff
-- short seqeuneces not worth it.  This should be tested empirically
thesholdUKLength :: Int 
thesholdUKLength = 1

-- | getDOMedian wraps around getPrelim and getPrelim3
-- changing types and checking for missing cases
getDOMedian :: V.Vector  BV.BitVector -> V.Vector  BV.BitVector -> V.Vector  (V.Vector  Int) 
            -> TCMD.DenseTransitionCostMatrix 
            -> MR.MetricRepresentation (AG.AmbiguityGroup -> AG.AmbiguityGroup -> (AG.AmbiguityGroup, Word)) 
            -> CharType -> (V.Vector  BV.BitVector, Int)
getDOMedian origLBV origRBV thisMatrix tcmDense tcmMemo thisType =
    -- missing data inputs
    if V.null origLBV then (origRBV, 0)
    else if V.null origRBV then (origLBV, 0)
    else 
        -- not missing
        -- get inDelCost 
        let inDelCost = V.head (V.last thisMatrix)
            (lBV, rBV) = setLeftRight origLBV origRBV
            -- Int64 version faster for small aphabets
            bvLength = BV.dimension (V.head lBV) 
            leftChar64 = V.map convertBVTo64 lBV
            rightChar64 = V.map convertBVTo64 rBV
            (newMedian64, medianCost64) = DKS64.ukkonenDO leftChar64 rightChar64 inDelCost
            newMedianBV = V.map (convert64ToBV bvLength) newMedian64

            {-
            tcmMemo' = 
                let sigma i j       = toEnum . fromEnum $ tcm TCM.! (fromEnum i, fromEnum j)
                    memoMatrixValue = TCMM.generateMemoizedTransitionCostMatrix (toEnum $ length arbitraryAlphabet) sigma
                in  ExplicitLayout tcm (TCMM.getMedianAndCost2D memoMatrixValue)
        
            (weight, tcm) = TCM.fromRows $ SM.getFullVects thisMatrix

            arbitraryAlphabet :: Alphabet String
            arbitraryAlphabet = fromSymbols $ show <$> 0 :| [1 .. length thisMatrix - 1]
            -}

            --setting left most bit to 1 same purpose as inDelBit for Ukkonen
            (newMedianSmall, medianCostSmall) = DKS.ukkonenDO lBV rBV inDelCost
            (newMedianLarge, medianCostLarge) = NS.naiveDO lBV rBV inDelCost

            (mediansFFI, costFFI) = DOSmallFFI.wrapperPCG_DO_Small_FFI lBV rBV thisMatrix tcmDense
            (mediansLargeFFI, costLargeFFI) = DOLargeFFI.wrapperPCG_DO_Large lBV rBV thisMatrix tcmMemo
        in
        --trace ("DO: " ++ (s(V.Vector (Int, Int, Int))how inDelCost) ++ " " ++ (show $ V.head $ V.last thisMatrix)) (
        -- Naive (ie no Ukkonene if short sequneces
        if (thisType == NucSeq || thisType == SmallAlphSeq) then (mediansFFI, costFFI)
        --else if (thisType == GenSeq || thisType == AminoSeq) then  (mediansLargeFFI, costLargeFFI) 
        else if (thisType == GenSeq || thisType == AminoSeq) then  (newMedianLarge, medianCostLarge)
        --if (min (V.length lBV)  (V.length rBV)) < thesholdUKLength then (newMedianLarge, medianCostLarge)
        --else if thisType == NucSeq then (newMedianBV, medianCost64)
        --else if thisType == GenSeq then (newMedianSmall, medianCostSmall)
        else error "Unrecognized/Not implemented character type"
       
        
-- | convertBVTo64 converts bitV.Vector  type to bit64 Int type 
convertBVTo64 :: BV.BitVector -> Int64
convertBVTo64 inBV = fromIntegral (BV.toSignedNumber inBV) 

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