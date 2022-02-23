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
  , atSymbolIndex
  ) where

import Control.DeepSeq
import Data.Data
import Data.Hashable
import Data.Int
import Data.Ix
import Data.Word
import Foreign.C.Types (CUInt)
import Foreign.Storable
import GHC.Generics


-- |
-- The index of a symbol in an alphabet.
newtype SymbolIndex = SymbolIndex Word
    deriving newtype (Enum, Eq, Hashable, Ix, NFData, Ord, Read, Show, Storable, Typeable)
    deriving stock   (Data, Generic)


-- |
-- Use a 'SymbolIndex' in place of an 'Integral' value for index related operations.
{-# SCC        atSymbolIndex #-}
{-# INLINEABLE atSymbolIndex #-}
{-# SPECIALISE atSymbolIndex :: SymbolIndex -> CUInt  #-}
{-# SPECIALISE atSymbolIndex :: SymbolIndex -> Int    #-}
{-# SPECIALISE atSymbolIndex :: SymbolIndex -> Int8   #-}
{-# SPECIALISE atSymbolIndex :: SymbolIndex -> Int16  #-}
{-# SPECIALISE atSymbolIndex :: SymbolIndex -> Int32  #-}
{-# SPECIALISE atSymbolIndex :: SymbolIndex -> Int64  #-}
{-# SPECIALISE atSymbolIndex :: SymbolIndex -> Word   #-}
{-# SPECIALISE atSymbolIndex :: SymbolIndex -> Word8  #-}
{-# SPECIALISE atSymbolIndex :: SymbolIndex -> Word16 #-}
{-# SPECIALISE atSymbolIndex :: SymbolIndex -> Word32 #-}
{-# SPECIALISE atSymbolIndex :: SymbolIndex -> Word64 #-}
atSymbolIndex :: Integral i => SymbolIndex -> i
atSymbolIndex (SymbolIndex i) = fromIntegral i
