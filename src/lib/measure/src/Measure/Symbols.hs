-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Symbols
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-- A matrix of transition costs between alphabet symbols.
-- Exposes row-major monomorphic maps, folds, and traversals.
-----------------------------------------------------------------------------

module Measure.Symbols
  ( TCM()
  , DiagnosisOfSCM(..)
  , StructureOfSCM(..)
    -- * Construction
  , fromCols
  , fromList
  , fromRows
  , generate
    -- * Indexing
  , (!)
    -- * Queries
  , size
    -- * Specialization Utility
  , diagnoseTcm
  ) where

import Measure.Symbols.Dense.Internal
