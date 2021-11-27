-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Matrix
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

module Measure.Matrix
  ( -- * Comparators
    SCMλ
  , TCM2Dλ
  , TCM3Dλ
    -- * Comparator Type-classes
  , HasDenseMatrix(..)
  , HasSymbolChangeMatrix(..)
  , HasTransitionCostMatrix(..)
  ) where

import Measure.Matrix.Types
import Measure.Matrix.Class
