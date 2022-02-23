-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Transition.States
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE MonoLocalBinds        #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE RankNTypes            #-}
{-# LANGUAGE Strict                #-}
{-# LANGUAGE TypeSynonymInstances  #-}

module Measure.Transition.States
  ( HasStateTransitions(..)
  ) where

import Data.Bits
import Data.Word                      (Word64)
import Foreign.C.Types                (CUInt)
import Layout.Compact.States.Indexing
import Measure.Transition.Types


-- |
-- Any structural representation which can produce a functions measuring the
-- 'Measure.Dispersion.Dispersion' between pairs and triples of states, i.e.
-- produce 2D and 3D "Transition Cost Matrices."
--
-- /NOTE:/ Measurability of 'Measure.Dispersion.Dispersion' implies measurability
-- of 'Measure.Centroid.Centroid' and 'Measure.Distance.Distance'.
class HasStateTransitions a e where

    {-# MINIMAL stateTransitionPairwiseDispersion, stateTransitionThreewayDispersion #-}

    stateTransitionPairwiseDispersion :: a -> StateTransitionPairwiseDispersionλ e

    stateTransitionPairwiseDistance   :: a -> StateTransitionPairwiseDistanceλ   e
    stateTransitionPairwiseDistance   = ((fst .) .) . stateTransitionPairwiseDispersion

    stateTransitionPairwiseCentroid   :: a -> StateTransitionPairwiseCentroidλ     e
    stateTransitionPairwiseCentroid   = ((snd .) .) . stateTransitionPairwiseDispersion

    stateTransitionThreewayDispersion :: a -> StateTransitionThreewayDispersionλ e

    stateTransitionThreewayDistance   :: a -> StateTransitionThreewayDistanceλ   e
    stateTransitionThreewayDistance   = (((fst .) .) .) . stateTransitionThreewayDispersion

    stateTransitionThreewayCentroid   :: a -> StateTransitionThreewayCentroidλ     e
    stateTransitionThreewayCentroid   = (((snd .) .) .) . stateTransitionThreewayDispersion


-- |
-- Generally less desireable instance.
instance Bits b => HasStateTransitions StateTransitionsCompact b where

    -- |
    -- /O(1)/
    stateTransitionPairwiseDispersion = both2D'

    -- |
    -- /O(1)/
    stateTransitionPairwiseDistance   = cost2D'

    -- |
    -- /O(1)/
    stateTransitionPairwiseCentroid   = mean2D'

    -- |
    -- /O(1)/
    stateTransitionThreewayDispersion = both3D'

    -- |
    -- /O(1)/
    stateTransitionThreewayDistance   = cost3D'

    -- |
    -- /O(1)/
    stateTransitionThreewayCentroid   = mean3D'


instance {-# OVERLAPPING #-} HasStateTransitions StateTransitionsCompact CUInt where

    -- |
    -- /O(1)/
    {-# SPECIALISE INLINE stateTransitionPairwiseDispersion :: StateTransitionsCompact -> StateTransitionPairwiseDispersionλ CUInt #-}
    stateTransitionPairwiseDispersion = both2D

    -- |
    -- /O(1)/
    {-# SPECIALISE INLINE stateTransitionPairwiseDistance   :: StateTransitionsCompact -> StateTransitionPairwiseDistanceλ   CUInt #-}
    stateTransitionPairwiseDistance   = cost2D

    -- |
    -- /O(1)/
    {-# SPECIALISE INLINE stateTransitionPairwiseCentroid   :: StateTransitionsCompact -> StateTransitionPairwiseCentroidλ   CUInt #-}
    stateTransitionPairwiseCentroid   = mean2D

    -- |
    -- /O(1)/
    {-# SPECIALISE INLINE stateTransitionThreewayDispersion :: StateTransitionsCompact -> StateTransitionThreewayDispersionλ CUInt #-}
    stateTransitionThreewayDispersion = both3D

    -- |
    -- /O(1)/
    {-# SPECIALISE INLINE stateTransitionThreewayDistance   :: StateTransitionsCompact -> StateTransitionThreewayDistanceλ   CUInt #-}
    stateTransitionThreewayDistance   = cost3D

    -- |
    -- /O(1)/
    {-# SPECIALISE INLINE stateTransitionThreewayCentroid   :: StateTransitionsCompact -> StateTransitionThreewayCentroidλ   CUInt #-}
    stateTransitionThreewayCentroid   = mean3D


instance {-# OVERLAPPING #-} HasStateTransitions StateTransitionsCompact Word64 where

    -- |
    -- /O(1)/
    {-# SPECIALISE INLINE stateTransitionPairwiseDispersion :: StateTransitionsCompact -> StateTransitionPairwiseDispersionλ Word64 #-}
    stateTransitionPairwiseDispersion = both2D

    -- |
    -- /O(1)/
    {-# SPECIALISE INLINE stateTransitionPairwiseDistance   :: StateTransitionsCompact -> StateTransitionPairwiseDistanceλ   Word64 #-}
    stateTransitionPairwiseDistance   = cost2D

    -- |
    -- /O(1)/
    {-# SPECIALISE INLINE stateTransitionPairwiseCentroid   :: StateTransitionsCompact -> StateTransitionPairwiseCentroidλ   Word64 #-}
    stateTransitionPairwiseCentroid   = mean2D

    -- |
    -- /O(1)/
    {-# SPECIALISE INLINE stateTransitionThreewayDispersion :: StateTransitionsCompact -> StateTransitionThreewayDispersionλ Word64 #-}
    stateTransitionThreewayDispersion = both3D

    -- |
    -- /O(1)/
    {-# SPECIALISE INLINE stateTransitionThreewayDistance   :: StateTransitionsCompact -> StateTransitionThreewayDistanceλ   Word64 #-}
    stateTransitionThreewayDistance   = cost3D

    -- |
    -- /O(1)/
    {-# SPECIALISE INLINE stateTransitionThreewayCentroid   :: StateTransitionsCompact -> StateTransitionThreewayCentroidλ   Word64 #-}
    stateTransitionThreewayCentroid   = mean3D
