-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Unit.SymbolChangeCost
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

module Measure.Unit.SymbolChangeCost
  ( SymbolChangeCost(..)
--  , HasEditCosts(..)
  ) where

import Control.DeepSeq
import Data.Data
import Data.Hashable
import GHC.Generics


-- |
-- The distance between two measurable elements.
newtype SymbolChangeCost = SymbolChangeCost Word
    deriving newtype (Bounded, Enum, Eq, Hashable, Integral, NFData, Num, Ord, Read, Real, Show, Typeable)
    deriving stock   (Data, Generic)


{-
-- |
-- Any structural representation which can produce a Symbol Change Matrix.
class HasEditCosts a where

    {-# MINIMAL maxDeletion, maxInsertion, minDeletion, minInsertion #-}

    {-# INLINEABLE maxEdit #-}
    maxEdit :: a -> SymbolChangeCost
    maxEdit a = max (maxDeletion a) $ maxInsertion a

    maxDeletion  :: a -> SymbolChangeCost

    maxInsertion :: a -> SymbolChangeCost

    {-# INLINEABLE minEdit #-}
    minEdit :: a -> SymbolChangeCost
    minEdit a = min (minDeletion a) $ minInsertion a

    minDeletion  :: a -> SymbolChangeCost

    minInsertion :: a -> SymbolChangeCost
-}
