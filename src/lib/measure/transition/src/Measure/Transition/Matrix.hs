-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Transition.Matrix
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

module Measure.Transition.Matrix
  ( -- * Comparators
    SCMλ
  , TCM2Dλ
  , TCM3Dλ
    -- * Comparator Type-classes
  , HasDenseMatrix(..)
  , HasSymbolChangeMatrix(..)
  , HasTransitionCostMatrix(..)
  ) where

import Measure.Transition.Matrix.Types
import Measure.Transition.Matrix.Class
