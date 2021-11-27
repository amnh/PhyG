-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Unit.Distance
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE DeriveDataTypeable         #-}
{-# LANGUAGE DeriveGeneric              #-}
{-# LANGUAGE DerivingStrategies         #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Measure.Unit.Distance
  ( Distance(..)
  , HasEditCosts(..)
  ) where

import Control.DeepSeq
import Data.Data
import GHC.Generics


-- |
-- The distance between two measurable elements.
newtype Distance = Distance Word
    deriving newtype (Bounded, Enum, Eq, Integral, NFData, Num, Ord, Read, Real, Show, Typeable)
    deriving stock   (Data, Generic)


-- |
-- Any structural representation which can produce a Symbol Change Matrix.
class HasEditCosts a where

    {-# MINIMAL maxDeletion, maxInsertion, minDeletion, minInsertion #-}

    {-# INLINEABLE maxEdit #-}
    maxEdit :: a -> Distance
    maxEdit a = max (maxDeletion a) $ maxInsertion a

    maxDeletion  :: a -> Distance

    maxInsertion :: a -> Distance

    {-# INLINEABLE minEdit #-}
    minEdit :: a -> Distance
    minEdit a = min (minDeletion a) $ minInsertion a

    minDeletion  :: a -> Distance

    minInsertion :: a -> Distance
