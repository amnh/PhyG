-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Unit.SymbolCount
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
{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Measure.Unit.SymbolCount
  ( SymbolCount(..)
  , HasSymbolCount(..)
  , iota
  , infimumSymbolLimit
  ) where

import Control.DeepSeq
import Data.Coerce
import Data.Data
import Data.Int
import Data.Word
import Foreign.C.Types
import GHC.Generics
import GHC.Natural


-- |
-- The index of a symbol in an alphabet.
newtype SymbolCount = SymbolCount Word
    deriving newtype (Bounded, Enum, Eq, Integral, NFData, Num, Ord, Read, Real, Show, Typeable)
    deriving stock   (Data, Generic)


-- |
-- A structure which can derive the number of alphabet symbols associated with it.
class HasSymbolCount a where

    symbolCount :: a -> SymbolCount


instance HasSymbolCount CChar where

    symbolCount = symbolCount . (coerce :: CChar -> Int8)


instance HasSymbolCount CSChar where

    symbolCount = symbolCount . (coerce :: CSChar -> Int8)


instance HasSymbolCount CUChar where

    symbolCount = symbolCount . (coerce :: CUChar -> Word8)


instance HasSymbolCount CShort where

    symbolCount = symbolCount . (coerce :: CShort -> Int16)


instance HasSymbolCount CUShort where

    symbolCount = symbolCount . (coerce :: CUShort -> Word16)


instance HasSymbolCount CInt where

    symbolCount = symbolCount . (coerce :: CInt -> Int32)


instance HasSymbolCount CUInt where

    symbolCount = symbolCount . (coerce :: CUInt -> Word32)


instance HasSymbolCount CLong where

    symbolCount = symbolCount . (coerce :: CLong -> Int64)


instance HasSymbolCount CULong where

    symbolCount = symbolCount . (coerce :: CULong -> Word64)


instance HasSymbolCount CSize where

    symbolCount = symbolCount . (coerce :: CSize -> Word64)


instance HasSymbolCount CLLong where

    symbolCount = symbolCount . (coerce :: CLLong -> Int64)


instance HasSymbolCount CULLong where

    symbolCount = symbolCount . (coerce :: CULLong -> Word64)


instance HasSymbolCount Int where

    symbolCount = clamp


instance HasSymbolCount Int8 where

    symbolCount = clamp 


instance HasSymbolCount Int16 where

    symbolCount = clamp 


instance HasSymbolCount Int32 where

    symbolCount = clamp 


instance HasSymbolCount Int64 where

    symbolCount = clamp 


instance HasSymbolCount Integer where

    symbolCount = clamp 


instance HasSymbolCount Natural where

    symbolCount = SymbolCount . naturalToWord


instance HasSymbolCount SymbolCount where

    symbolCount = id


instance HasSymbolCount Word where

    symbolCount = SymbolCount


instance HasSymbolCount Word8 where

    symbolCount = SymbolCount . fromIntegral 


instance HasSymbolCount Word16 where

    symbolCount = SymbolCount . fromIntegral 


instance HasSymbolCount Word32 where

    symbolCount = SymbolCount . fromIntegral 


instance HasSymbolCount Word64 where

    symbolCount = SymbolCount . fromIntegral 


-- |
-- The largest 'SymbolCount' value for which the predicate 'iota' holds.
--
-- Useful for partitioning a collection of symbols based on whether it is too
-- large for the C FFI.
infimumSymbolLimit :: SymbolCount
infimumSymbolLimit = 8


-- |
-- Predicate to deciding if a 'SymbolCount' is small enough to be compatible with
-- the C FFI.
iota :: HasSymbolCount a => a -> Bool
iota = (<= infimumSymbolLimit) . symbolCount


clamp :: Integral i => i -> SymbolCount
clamp i
  | i <= 0    = SymbolCount 0
  | otherwise = fromIntegral i
