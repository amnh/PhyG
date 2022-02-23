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
  ) where

import Control.DeepSeq
import Data.Data
import Data.Hashable
import Foreign.Storable
import GHC.Generics


-- |
-- The distance between two measurable elements.
newtype SymbolChangeCost = SymbolChangeCost Word
    deriving newtype (Bounded, Enum, Eq, Hashable, Integral, NFData, Num, Ord, Read, Real, Show, Storable, Typeable)
    deriving stock   (Data, Generic)

