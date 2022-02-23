-----------------------------------------------------------------------------
-- |
-- Module      :  TransitionMatrix
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

module TransitionMatrix
  ( -- * Specialized Representation
    TransitionMatrix()
  , TransitionMeasureDiagnosis(..)
    -- * Special Constructors
  , discreteCrossGap
  , discreteMetric
  , linearNorm
    -- * General Constructors
  , fromList
  , fromColumns
  , fromRows
  , generate
  ) where

import TransitionMatrix.Diagnosis
