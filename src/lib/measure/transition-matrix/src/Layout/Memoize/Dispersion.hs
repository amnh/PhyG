-----------------------------------------------------------------------------
-- |
-- Module      :  Layout.Memoize.Dispersion
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MonoLocalBinds   #-}
{-# LANGUAGE RankNTypes       #-}
{-# LANGUAGE Strict           #-}

module Layout.Memoize.Dispersion
  ( bitDispersion
  , bitDispersionPairwise
  , bitDispersionThreeway
  ) where

import Data.Bits
import Data.Foldable
import Data.Ix
import Data.List.NonEmpty            (NonEmpty(..))
import Data.Semigroup
import Data.Semigroup.Foldable
import Data.Word
import Foreign.C.Types               (CUInt)
import Measure.Dispersion
import Measure.Distance
import Measure.Unit.SymbolChangeCost
import Measure.Unit.SymbolIndex


data  Bounds b i
    = Bounds
    { _lBound :: i
    , _uBound :: i
    , _bValue :: b
    }


-- |
-- Takes one or more elements of 'FiniteBits' and a symbol change cost function
-- and returns a tuple of a new character, along with the cost of obtaining that
-- character. The return character may be (or is even likely to be) ambiguous.
-- Will attempt to intersect the two characters, but will union them if that is
-- not possible, based on the symbol change cost function.
--
-- To clarify, the return character is an intersection of all possible least-cost
-- combinations, so for instance, if @ char1 == A,T @ and @ char2 == G,C @, and
-- the two (non-overlapping) least cost pairs are A,C and T,G, then the return
-- value is A,C,G,T.
{-# INLINEABLE bitDispersion #-}
{-# SPECIALISE bitDispersion :: (Bounded c, Enum i, Integral c, Ix i, Num d) => (i, i) -> Distance c i -> Dispersion d CUInt  #-}
{-# SPECIALISE bitDispersion :: (Bounded c, Enum i, Integral c, Ix i, Num d) => (i, i) -> Distance c i -> Dispersion d Word   #-}
{-# SPECIALISE bitDispersion :: (Bounded c, Enum i, Integral c, Ix i, Num d) => (i, i) -> Distance c i -> Dispersion d Word8  #-}
{-# SPECIALISE bitDispersion :: (Bounded c, Enum i, Integral c, Ix i, Num d) => (i, i) -> Distance c i -> Dispersion d Word16 #-}
{-# SPECIALISE bitDispersion :: (Bounded c, Enum i, Integral c, Ix i, Num d) => (i, i) -> Distance c i -> Dispersion d Word32 #-}
{-# SPECIALISE bitDispersion :: (Bounded c, Enum i, Integral c, Ix i, Num d) => (i, i) -> Distance c i -> Dispersion d Word64 #-}
{-# SPECIALISE bitDispersion :: (Bounded c, Integral c, Num d) => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> Dispersion d CUInt  #-}
{-# SPECIALISE bitDispersion :: (Bounded c, Integral c, Num d) => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> Dispersion d Word   #-}
{-# SPECIALISE bitDispersion :: (Bounded c, Integral c, Num d) => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> Dispersion d Word8  #-}
{-# SPECIALISE bitDispersion :: (Bounded c, Integral c, Num d) => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> Dispersion d Word16 #-}
{-# SPECIALISE bitDispersion :: (Bounded c, Integral c, Num d) => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> Dispersion d Word32 #-}
{-# SPECIALISE bitDispersion :: (Bounded c, Integral c, Num d) => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> Dispersion d Word64 #-}
{-# SPECIALISE bitDispersion :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> Dispersion SymbolChangeCost CUInt  #-}
{-# SPECIALISE bitDispersion :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> Dispersion SymbolChangeCost Word   #-}
{-# SPECIALISE bitDispersion :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> Dispersion SymbolChangeCost Word8  #-}
{-# SPECIALISE bitDispersion :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> Dispersion SymbolChangeCost Word16 #-}
{-# SPECIALISE bitDispersion :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> Dispersion SymbolChangeCost Word32 #-}
{-# SPECIALISE bitDispersion :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> Dispersion SymbolChangeCost Word64 #-}
{-# SPECIALISE bitDispersion :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> Dispersion SymbolChangeCost CUInt  #-}
{-# SPECIALISE bitDispersion :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> Dispersion SymbolChangeCost Word   #-}
{-# SPECIALISE bitDispersion :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> Dispersion SymbolChangeCost Word8  #-}
{-# SPECIALISE bitDispersion :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> Dispersion SymbolChangeCost Word16 #-}
{-# SPECIALISE bitDispersion :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> Dispersion SymbolChangeCost Word32 #-}
{-# SPECIALISE bitDispersion :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> Dispersion SymbolChangeCost Word64 #-}
bitDispersion
  :: ( Bounded c
     , Enum i
     , FiniteBits b
     , Integral c
     , Ix i
     , Num d
     )
  => (i, i)         -- ^ Magnitude
  -> Distance   c i -- ^ Symbol change matrix (SDM) to determine cost
  -> Dispersion d b -- ^ List of elements for of which to find the k-median and cost
bitDispersion ixBounds sigma xs = (fromIntegral distance, medoid)
  where
    (distance, medoid) = foldl' processRange (maxBound, zero) $ range ixBounds
    withBounds = getBitBounds <$> toList xs
    wlog  = getFirst $ foldMap1 First xs
    zero  = wlog `xor` wlog

    processRange (oldCost, bits) i =
        let newCost = foldl' (+) 0 $ getDistance i <$> withBounds
        in  case oldCost `compare` newCost of
              LT -> (oldCost, bits                    )
              EQ -> (oldCost, bits `setBit` fromEnum i)
              GT -> (newCost, zero `setBit` fromEnum i)

    getDistance i (Bounds lo hi b) = foldl' processSubrange maxBound $ range (lo, hi)
      where
        processSubrange cost j
          | b `testBit` fromEnum j = min cost $ sigma i j
          | otherwise              = cost



-- |
-- Calculate the median between /two/ states.
{-# INLINEABLE bitDispersionPairwise #-}
{-# SPECIALISE bitDispersionPairwise :: (Bounded c, Enum i, Integral c, Ix i, Num d) => (i, i) -> Distance c i -> DispersionPairwise d CUInt  #-}
{-# SPECIALISE bitDispersionPairwise :: (Bounded c, Enum i, Integral c, Ix i, Num d) => (i, i) -> Distance c i -> DispersionPairwise d Word   #-}
{-# SPECIALISE bitDispersionPairwise :: (Bounded c, Enum i, Integral c, Ix i, Num d) => (i, i) -> Distance c i -> DispersionPairwise d Word8  #-}
{-# SPECIALISE bitDispersionPairwise :: (Bounded c, Enum i, Integral c, Ix i, Num d) => (i, i) -> Distance c i -> DispersionPairwise d Word16 #-}
{-# SPECIALISE bitDispersionPairwise :: (Bounded c, Enum i, Integral c, Ix i, Num d) => (i, i) -> Distance c i -> DispersionPairwise d Word32 #-}
{-# SPECIALISE bitDispersionPairwise :: (Bounded c, Enum i, Integral c, Ix i, Num d) => (i, i) -> Distance c i -> DispersionPairwise d Word64 #-}
{-# SPECIALISE bitDispersionPairwise :: (Bounded c, Integral c, Num d) => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionPairwise d CUInt  #-}
{-# SPECIALISE bitDispersionPairwise :: (Bounded c, Integral c, Num d) => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionPairwise d Word   #-}
{-# SPECIALISE bitDispersionPairwise :: (Bounded c, Integral c, Num d) => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionPairwise d Word8  #-}
{-# SPECIALISE bitDispersionPairwise :: (Bounded c, Integral c, Num d) => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionPairwise d Word16 #-}
{-# SPECIALISE bitDispersionPairwise :: (Bounded c, Integral c, Num d) => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionPairwise d Word32 #-}
{-# SPECIALISE bitDispersionPairwise :: (Bounded c, Integral c, Num d) => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionPairwise d Word64 #-}
{-# SPECIALISE bitDispersionPairwise :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionPairwise SymbolChangeCost CUInt  #-}
{-# SPECIALISE bitDispersionPairwise :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionPairwise SymbolChangeCost Word   #-}
{-# SPECIALISE bitDispersionPairwise :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionPairwise SymbolChangeCost Word8  #-}
{-# SPECIALISE bitDispersionPairwise :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionPairwise SymbolChangeCost Word16 #-}
{-# SPECIALISE bitDispersionPairwise :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionPairwise SymbolChangeCost Word32 #-}
{-# SPECIALISE bitDispersionPairwise :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionPairwise SymbolChangeCost Word64 #-}
{-# SPECIALISE bitDispersionPairwise :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionPairwise SymbolChangeCost CUInt  #-}
{-# SPECIALISE bitDispersionPairwise :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionPairwise SymbolChangeCost Word   #-}
{-# SPECIALISE bitDispersionPairwise :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionPairwise SymbolChangeCost Word8  #-}
{-# SPECIALISE bitDispersionPairwise :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionPairwise SymbolChangeCost Word16 #-}
{-# SPECIALISE bitDispersionPairwise :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionPairwise SymbolChangeCost Word32 #-}
{-# SPECIALISE bitDispersionPairwise :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionPairwise SymbolChangeCost Word64 #-}
bitDispersionPairwise
  :: ( Bounded c
     , Enum i
     , FiniteBits b
     , Integral c
     , Ix i
     , Num d
     )
  => (i, i)
  -> Distance           c i
  -> DispersionPairwise d b
bitDispersionPairwise size sigma char1 char2 = bitDispersion size sigma $ char1 :| [char2]


-- |
-- Calculate the median between /three/ states.
{-# INLINEABLE bitDispersionThreeway #-}
{-# SPECIALISE bitDispersionThreeway :: (Bounded c, Enum i, Integral c, Ix i, Num d) => (i, i) -> Distance c i -> DispersionThreeway d CUInt  #-}
{-# SPECIALISE bitDispersionThreeway :: (Bounded c, Enum i, Integral c, Ix i, Num d) => (i, i) -> Distance c i -> DispersionThreeway d Word   #-}
{-# SPECIALISE bitDispersionThreeway :: (Bounded c, Enum i, Integral c, Ix i, Num d) => (i, i) -> Distance c i -> DispersionThreeway d Word8  #-}
{-# SPECIALISE bitDispersionThreeway :: (Bounded c, Enum i, Integral c, Ix i, Num d) => (i, i) -> Distance c i -> DispersionThreeway d Word16 #-}
{-# SPECIALISE bitDispersionThreeway :: (Bounded c, Enum i, Integral c, Ix i, Num d) => (i, i) -> Distance c i -> DispersionThreeway d Word32 #-}
{-# SPECIALISE bitDispersionThreeway :: (Bounded c, Enum i, Integral c, Ix i, Num d) => (i, i) -> Distance c i -> DispersionThreeway d Word64 #-}
{-# SPECIALISE bitDispersionThreeway :: (Bounded c, Integral c, Num d) => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionThreeway d CUInt  #-}
{-# SPECIALISE bitDispersionThreeway :: (Bounded c, Integral c, Num d) => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionThreeway d Word   #-}
{-# SPECIALISE bitDispersionThreeway :: (Bounded c, Integral c, Num d) => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionThreeway d Word8  #-}
{-# SPECIALISE bitDispersionThreeway :: (Bounded c, Integral c, Num d) => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionThreeway d Word16 #-}
{-# SPECIALISE bitDispersionThreeway :: (Bounded c, Integral c, Num d) => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionThreeway d Word32 #-}
{-# SPECIALISE bitDispersionThreeway :: (Bounded c, Integral c, Num d) => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionThreeway d Word64 #-}
{-# SPECIALISE bitDispersionThreeway :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionThreeway SymbolChangeCost CUInt  #-}
{-# SPECIALISE bitDispersionThreeway :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionThreeway SymbolChangeCost Word   #-}
{-# SPECIALISE bitDispersionThreeway :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionThreeway SymbolChangeCost Word8  #-}
{-# SPECIALISE bitDispersionThreeway :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionThreeway SymbolChangeCost Word16 #-}
{-# SPECIALISE bitDispersionThreeway :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionThreeway SymbolChangeCost Word32 #-}
{-# SPECIALISE bitDispersionThreeway :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionThreeway SymbolChangeCost Word64 #-}
{-# SPECIALISE bitDispersionThreeway :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionThreeway SymbolChangeCost CUInt  #-}
{-# SPECIALISE bitDispersionThreeway :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionThreeway SymbolChangeCost Word   #-}
{-# SPECIALISE bitDispersionThreeway :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionThreeway SymbolChangeCost Word8  #-}
{-# SPECIALISE bitDispersionThreeway :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionThreeway SymbolChangeCost Word16 #-}
{-# SPECIALISE bitDispersionThreeway :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionThreeway SymbolChangeCost Word32 #-}
{-# SPECIALISE bitDispersionThreeway :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionThreeway SymbolChangeCost Word64 #-}
bitDispersionThreeway
  :: ( Bounded c
     , Enum i
     , FiniteBits b
     , Integral c
     , Ix i
     , Num d
     )
  => (i, i)
  -> Distance           c i
  -> DispersionThreeway d b
bitDispersionThreeway size sigma char1 char2 char3 = bitDispersion size sigma $ char1 :| [char2, char3]


-- |
-- Gets the lowest set bit and the highest set bit in the collection.
getBitBounds
  ::
  ( Enum i
  , FiniteBits b
  )
  => b
  -> Bounds b i
getBitBounds b =
    let bitZero   = (b `xor` b) `setBit` 0
        bigEndian = countLeadingZeros bitZero > 0 -- Check the endianness

        (f,g) | bigEndian = (countTrailingZeros, countLeadingZeros )
              | otherwise = (countLeadingZeros , countTrailingZeros)

        lZeroes = f b
        uZeroes = g b
        lower   = toEnum $ lZeroes
        upper   = toEnum . max 0 $ finiteBitSize b - uZeroes - 1
    in  Bounds lower upper b
