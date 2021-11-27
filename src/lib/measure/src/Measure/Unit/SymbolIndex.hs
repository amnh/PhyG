-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Unit.SymbolIndex
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

module Measure.Unit.SymbolIndex
  ( SymbolIndex(..)
  ) where

import Control.DeepSeq
import Data.Data
import GHC.Generics


-- |
-- The index of a symbol in an alphabet.
newtype SymbolIndex = SymbolIndex Word
    deriving newtype (Bounded, Enum, Eq, Integral, NFData, Num, Ord, Read, Real, Show, Typeable)
    deriving stock   (Data, Generic)
