-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Transition
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Measure.Transition
  ( -- * Core "Transition Matrix" Type-class
    HasTransitionMatrix
    -- * Symbol Measure
  , SymbolDistanceλ
    -- * State Pairwise Measures
  , StateTransitionPairwiseCentroidλ
  , StateTransitionPairwiseDispersionλ
  , StateTransitionPairwiseDistanceλ
  -- * State Threeway Measures
  , StateTransitionThreewayCentroidλ
  , StateTransitionThreewayDispersionλ
  , StateTransitionThreewayDistanceλ
  -- * Synonyms
  , SDMλ
  , TCM2Dλ
  , TCM3Dλ
  -- * Class
  , HasStateTransitions(..)
  , HasStateTransitionsCompact(..)
  , HasSymbolDistances(..)
  , HasEditExtrema(..)
  ) where

import Measure.Range
import Measure.Transition.Compact
import Measure.Transition.Edits
import Measure.Transition.States
import Measure.Transition.Symbols
import Measure.Transition.Types
import Measure.Unit.SymbolCount
import Measure.Unit.SymbolIndex


-- |
-- Collection of measure functionality required for facilitating string alignment.
--
-- Supports the following measures:
--
--   * "Dimensionality"
--     - 'symbolCount'
--
--   * 'Measure.Range.InclusiveRange' of "Symbols"
--     - 'Measure.Range.measureRange'
--
--   * 'Measure.Dispersion.Dispersion'
--     - 'stateTransitionPairwiseDispersion' : state ⨉ state         → (ℕ , state)
--     - 'stateTransitionThreewayDispersion' : state ⨉ state ⨉ state → (ℕ , state)
--
--   * 'Measure.Centroid.Centroid'
--     - 'stateTransitionPairwiseCentroid'   : state ⨉ state         → state
--     - 'stateTransitionThreewayCentroid'   : state ⨉ state ⨉ state → state
--
--   * 'Measure.Distance.Distance'
--     - 'symbolDistances'                   : symbol ⨉ symbol         → ℕ
--     - 'stateTransitionPairwiseDistance'   : state  ⨉ state          → ℕ
--     - 'stateTransitionThreewayDistance'   : state  ⨉ state  ⨉ state → ℕ
--
--   * "Symbol Edit Distance"
--     - 'maxEdit'      : ℕ
--     - 'maxDeletion'  : ℕ
--     - 'maxInsertion' : ℕ
--     - 'minEdit'      : ℕ
--     - 'minDeletion'  : ℕ
--     - 'minInsertion' : ℕ
--
class
    ( HasEditExtrema             a
    , HasStateTransitions        a b
    , HasStateTransitionsCompact a
    , HasSymbolCount             a
    , HasSymbolDistances         a
    , MeasurableRange            a SymbolIndex
    ) => HasTransitionMatrix a b where
{-
    transitionList

    transitionColumns

    transitionRows

    transitionGeneration
-}
