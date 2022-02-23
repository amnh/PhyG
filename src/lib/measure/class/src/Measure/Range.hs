-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Range
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE MultiParamTypeClasses #-}

module Measure.Range
  ( MeasurableRange(..)
  , InclusiveRange
  ) where

import Data.Ix


-- |
-- The inclusive lower and upper bounds of a range, respectively.
type InclusiveRange i = (i, i)


-- |
-- Structure which measures a contiguous 'InclusiveRange'.
class Ix i => MeasurableRange a i where

    measureRange :: a -> InclusiveRange i
