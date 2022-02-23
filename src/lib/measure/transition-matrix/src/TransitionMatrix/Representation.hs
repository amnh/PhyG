-----------------------------------------------------------------------------
-- |
-- Module      :  TransitionMatrix.Representation
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE DeriveAnyClass        #-}
{-# LANGUAGE DeriveGeneric         #-}
{-# LANGUAGE DerivingStrategies    #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE LambdaCase            #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE Strict                #-}
{-# LANGUAGE UnboxedSums           #-}

module TransitionMatrix.Representation
  ( -- * Specialized Representation
    TransitionMatrix(..)
{-
    -- * Measures
    -- * Special Constructors
  , discreteCrossGap
  , discreteMetric
  , linearNorm
    -- * General Constructors
  , fromList
  , fromColumns
  , fromRows
  , generate
-}
  ) where

import           Control.DeepSeq
import           Data.Bits
import           Data.Coerce
import           Data.Hashable
import           Data.Word                                          (Word64)
import           Foreign.C.Types                                    (CUInt)
import           GHC.Generics
import           Measure.Transition
import           Measure.Range
import           Layout.Special.Type
import           Layout.Memoize.States
import           Layout.Compact.States
import           Measure.Unit.SymbolCount
import           Measure.Unit.SymbolIndex


-- |
-- Represents the metric for some discrete characters or dynamic characters.
-- The representation is highly optimized for both time and space efficiency.
-- Because of this, the term "compact" is used to describe the representation.
--
-- If any measure has an 'iota' symbol count, i.e. a symbol count of less than
-- or equal to 'infimumSymbolLimit', the compact representation will pre-compute
-- the entire transition cost matrix in a structure which is 'Storable' and
-- inter-operable with the C FFI. If and only if the measure is specified for an
-- 'iota' symbol count will 'getCompactPairwise' return a @Just@ value.
--
-- Additionally, the compact representation notes if the discrete metric, the
-- discrete metric adjoined by the gap symbol, or the L1 norm are the specified
-- measure. If any of these metrics are specified, specialized functions which
-- are more efficient will be returned when calling 'symbolDistances', 'stateTransitionPairwiseDispersion', and
-- 'getStateTransitionCubeλ'.
--
-- Notably, if it is the case that /both/ the measure has an 'iota' symbol count
-- /and/ the measure is a specialized metric described above, then calling 'getCompactPairwise'
-- returns a /compile-time/ pre-computed transition cost matrix. Multiple
-- "constructions" of the same metric will result in the same "singlton" compact
-- representation.
--
-- Finally, if a specified measure has more than an 'iota' symbol count /and/ is
-- not one of the specialized metrics described above, then a /sparse/ and
-- /memoized/ representation is stored.
--
-- It is important to use this type in the metadata decorations rather than store
-- a function because a function cannot be contained in a compact region.
--
-- Use the 'symbolDistances', 'stateTransitionPairwiseDispersion', 'getStateTransitionCubeλ', and 'getCompactPairwise' to the retrieve the
-- desired functions.
data  TransitionMatrix a
    = IotaMatrix      {-# UNPACK #-} (Maybe SpecializableMetric) {-# UNPACK #-} StateTransitionsCompact
    | VastMatrix      {-# UNPACK #-} (Sparse a)
    | VastSpecialized {-# UNPACK #-} SpecializableMetric
    deriving stock    (Generic)
    deriving anyclass (NFData)


instance Eq (TransitionMatrix a) where

    (==) (VastMatrix            x   ) (VastMatrix            x'   ) = x == x'
    (==) (VastSpecialized       x   ) (VastSpecialized       x'   ) = x == x'
    (==) (IotaMatrix      (Just x) _) (IotaMatrix      (Just x') _) = x == x'
    (==) (IotaMatrix       Nothing x) (IotaMatrix       Nothing x') = x == x'
    (==) _ _ = False
    

instance HasEditExtrema (TransitionMatrix a) where

    {-# INLINEABLE maxEdit #-}
    maxEdit      (IotaMatrix  _ tcm) = maxEdit      tcm
    maxEdit      (VastMatrix    sdm) = maxEdit      sdm
    maxEdit      (VastSpecialized m) = maxEdit        m

    maxDeletion  (IotaMatrix  _ tcm) = maxDeletion  tcm
    maxDeletion  (VastMatrix    sdm) = maxDeletion  sdm
    maxDeletion  (VastSpecialized m) = maxDeletion    m

    maxInsertion (IotaMatrix  _ tcm) = maxInsertion tcm
    maxInsertion (VastMatrix    sdm) = maxInsertion sdm
    maxInsertion (VastSpecialized m) = maxInsertion   m

    {-# INLINEABLE minEdit #-}
    minEdit      (IotaMatrix  _ tcm) = minEdit      tcm
    minEdit      (VastMatrix    sdm) = minEdit      sdm
    minEdit      (VastSpecialized m) = minEdit        m

    minDeletion  (IotaMatrix  _ tcm) = minDeletion  tcm
    minDeletion  (VastMatrix    sdm) = minDeletion  sdm
    minDeletion  (VastSpecialized m) = minDeletion    m

    minInsertion (IotaMatrix  _ tcm) = minInsertion tcm
    minInsertion (VastMatrix    sdm) = minInsertion sdm
    minInsertion (VastSpecialized m) = minInsertion   m


instance HasSymbolDistances (TransitionMatrix a) where

    symbolDistances (IotaMatrix  _ tcm) = symbolDistances tcm
    symbolDistances (VastMatrix    sdm) = symbolDistances sdm
    symbolDistances (VastSpecialized m) = symbolDistances   m


instance HasSymbolCount (TransitionMatrix a) where

    symbolCount (IotaMatrix  _ tcm) = coerce $ matrixSize tcm
    symbolCount (VastMatrix    sdm) = symbolCount sdm
    symbolCount (VastSpecialized m) = symbolCount   m


instance HasStateTransitionsCompact (TransitionMatrix a) where

    getCompactPairwise (IotaMatrix _ tcm) = getCompactPairwise tcm
    getCompactPairwise  _                 = Nothing

    getCompactThreeway (IotaMatrix _ tcm) = getCompactThreeway tcm
    getCompactThreeway  _                 = Nothing



-- |
-- Generally less desireable instance.
instance (FiniteBits a, Hashable a) => HasStateTransitions (TransitionMatrix a) a where

    -- |
    -- /O(1)/
    stateTransitionPairwiseDispersion (IotaMatrix  _ tcm) = stateTransitionPairwiseDispersion tcm
    stateTransitionPairwiseDispersion (VastMatrix    sdm) = stateTransitionPairwiseDispersion sdm
    stateTransitionPairwiseDispersion (VastSpecialized m) = stateTransitionPairwiseDispersion   m

    -- |
    -- /O(1)/
    stateTransitionPairwiseDistance   (IotaMatrix  _ tcm) = stateTransitionPairwiseDistance   tcm
    stateTransitionPairwiseDistance   (VastMatrix    sdm) = stateTransitionPairwiseDistance   sdm
    stateTransitionPairwiseDistance   (VastSpecialized m) = stateTransitionPairwiseDistance     m

    -- |
    -- /O(1)/
    stateTransitionPairwiseCentroid   (IotaMatrix  _ tcm) = stateTransitionPairwiseCentroid   tcm
    stateTransitionPairwiseCentroid   (VastMatrix    sdm) = stateTransitionPairwiseCentroid   sdm
    stateTransitionPairwiseCentroid   (VastSpecialized m) = stateTransitionPairwiseCentroid     m

    -- |
    -- /O(1)/
    stateTransitionThreewayDispersion (IotaMatrix  _ tcm) = stateTransitionThreewayDispersion tcm
    stateTransitionThreewayDispersion (VastMatrix    sdm) = stateTransitionThreewayDispersion sdm
    stateTransitionThreewayDispersion (VastSpecialized m) = stateTransitionThreewayDispersion   m

    -- |
    -- /O(1)/
    stateTransitionThreewayDistance   (IotaMatrix  _ tcm) = stateTransitionThreewayDistance   tcm
    stateTransitionThreewayDistance   (VastMatrix    sdm) = stateTransitionThreewayDistance   sdm
    stateTransitionThreewayDistance   (VastSpecialized m) = stateTransitionThreewayDistance     m

    -- |
    -- /O(1)/
    stateTransitionThreewayCentroid   (IotaMatrix  _ tcm) = stateTransitionThreewayCentroid   tcm
    stateTransitionThreewayCentroid   (VastMatrix    sdm) = stateTransitionThreewayCentroid   sdm
    stateTransitionThreewayCentroid   (VastSpecialized m) = stateTransitionThreewayCentroid     m



instance {-# OVERLAPPING #-} HasStateTransitions (TransitionMatrix CUInt) CUInt where

    -- |
    -- /O(1)/
    stateTransitionPairwiseDispersion (IotaMatrix  _ tcm) = stateTransitionPairwiseDispersion tcm
    stateTransitionPairwiseDispersion (VastMatrix    sdm) = stateTransitionPairwiseDispersion sdm
    stateTransitionPairwiseDispersion (VastSpecialized m) = stateTransitionPairwiseDispersion   m

    -- |
    -- /O(1)/
    stateTransitionPairwiseDistance   (IotaMatrix  _ tcm) = stateTransitionPairwiseDistance   tcm
    stateTransitionPairwiseDistance   (VastMatrix    sdm) = stateTransitionPairwiseDistance   sdm
    stateTransitionPairwiseDistance   (VastSpecialized m) = stateTransitionPairwiseDistance     m

    -- |
    -- /O(1)/
    stateTransitionPairwiseCentroid   (IotaMatrix  _ tcm) = stateTransitionPairwiseCentroid   tcm
    stateTransitionPairwiseCentroid   (VastMatrix    sdm) = stateTransitionPairwiseCentroid   sdm
    stateTransitionPairwiseCentroid   (VastSpecialized m) = stateTransitionPairwiseCentroid     m

    -- |
    -- /O(1)/
    stateTransitionThreewayDispersion (IotaMatrix  _ tcm) = stateTransitionThreewayDispersion tcm
    stateTransitionThreewayDispersion (VastMatrix    sdm) = stateTransitionThreewayDispersion sdm
    stateTransitionThreewayDispersion (VastSpecialized m) = stateTransitionThreewayDispersion   m

    -- |
    -- /O(1)/
    stateTransitionThreewayDistance   (IotaMatrix  _ tcm) = stateTransitionThreewayDistance   tcm
    stateTransitionThreewayDistance   (VastMatrix    sdm) = stateTransitionThreewayDistance   sdm
    stateTransitionThreewayDistance   (VastSpecialized m) = stateTransitionThreewayDistance     m

    -- |
    -- /O(1)/
    stateTransitionThreewayCentroid   (IotaMatrix  _ tcm) = stateTransitionThreewayCentroid   tcm
    stateTransitionThreewayCentroid   (VastMatrix    sdm) = stateTransitionThreewayCentroid   sdm
    stateTransitionThreewayCentroid   (VastSpecialized m) = stateTransitionThreewayCentroid     m


instance {-# OVERLAPPING #-} HasStateTransitions (TransitionMatrix Word64) Word64 where

    -- |
    -- /O(1)/
    stateTransitionPairwiseDispersion (IotaMatrix  _ tcm) = stateTransitionPairwiseDispersion tcm
    stateTransitionPairwiseDispersion (VastMatrix    sdm) = stateTransitionPairwiseDispersion sdm
    stateTransitionPairwiseDispersion (VastSpecialized m) = stateTransitionPairwiseDispersion   m

    -- |
    -- /O(1)/
    stateTransitionPairwiseDistance   (IotaMatrix  _ tcm) = stateTransitionPairwiseDistance   tcm
    stateTransitionPairwiseDistance   (VastMatrix    sdm) = stateTransitionPairwiseDistance   sdm
    stateTransitionPairwiseDistance   (VastSpecialized m) = stateTransitionPairwiseDistance     m

    -- |
    -- /O(1)/
    stateTransitionPairwiseCentroid   (IotaMatrix  _ tcm) = stateTransitionPairwiseCentroid   tcm
    stateTransitionPairwiseCentroid   (VastMatrix    sdm) = stateTransitionPairwiseCentroid   sdm
    stateTransitionPairwiseCentroid   (VastSpecialized m) = stateTransitionPairwiseCentroid     m

    -- |
    -- /O(1)/
    stateTransitionThreewayDispersion (IotaMatrix  _ tcm) = stateTransitionThreewayDispersion tcm
    stateTransitionThreewayDispersion (VastMatrix    sdm) = stateTransitionThreewayDispersion sdm
    stateTransitionThreewayDispersion (VastSpecialized m) = stateTransitionThreewayDispersion   m

    -- |
    -- /O(1)/
    stateTransitionThreewayDistance   (IotaMatrix  _ tcm) = stateTransitionThreewayDistance   tcm
    stateTransitionThreewayDistance   (VastMatrix    sdm) = stateTransitionThreewayDistance   sdm
    stateTransitionThreewayDistance   (VastSpecialized m) = stateTransitionThreewayDistance     m

    -- |
    -- /O(1)/
    stateTransitionThreewayCentroid   (IotaMatrix  _ tcm) = stateTransitionThreewayCentroid   tcm
    stateTransitionThreewayCentroid   (VastMatrix    sdm) = stateTransitionThreewayCentroid   sdm
    stateTransitionThreewayCentroid   (VastSpecialized m) = stateTransitionThreewayCentroid     m


instance MeasurableRange (TransitionMatrix a) SymbolIndex where

    measureRange (IotaMatrix  _ tcm) = (coerce (0:: Word), coerce . pred $ matrixSize tcm)
    measureRange (VastMatrix    sdm) = measureRange sdm
    measureRange (VastSpecialized m) = measureRange   m


instance Show (TransitionMatrix a) where

    show =
        \case
            IotaMatrix (Just sm) tcm ->
                let dimension = fromEnum $ matrixSize tcm
                in  renderWithBytes dimension sm
            IotaMatrix  _ tcm -> renderSummary tcm
            VastMatrix    sdm -> show sdm
            VastSpecialized m -> show   m
