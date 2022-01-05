-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Dispersion.General
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

module Measure.Dispersion.General
  ( dispersion
  , dispersionPairwise
  , dispersionThreeway
  ) where

import Data.Bits
import Data.Foldable
import Data.Ix
import Data.List.NonEmpty (NonEmpty(..))
import Data.Semigroup
import Data.Semigroup.Foldable
import Measure.Distance
import Measure.Dispersion.Types
import Data.Word
import Foreign.C.Types         (CUInt)
import Measure.Unit.SymbolIndex
import Measure.Unit.SymbolChangeCost


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
-- the two (non-dispersionping) least cost pairs are A,C and T,G, then the return
-- value is A,C,G,T.
{-# INLINEABLE dispersion #-}
{-# SPECIALISE dispersion :: (Bounded c, Enum i, Ix i) => (i, i) -> Distance c i -> Dispersion c CUInt  #-}
{-# SPECIALISE dispersion :: (Bounded c, Enum i, Ix i) => (i, i) -> Distance c i -> Dispersion c Word   #-}
{-# SPECIALISE dispersion :: (Bounded c, Enum i, Ix i) => (i, i) -> Distance c i -> Dispersion c Word8  #-}
{-# SPECIALISE dispersion :: (Bounded c, Enum i, Ix i) => (i, i) -> Distance c i -> Dispersion c Word16 #-}
{-# SPECIALISE dispersion :: (Bounded c, Enum i, Ix i) => (i, i) -> Distance c i -> Dispersion c Word32 #-}
{-# SPECIALISE dispersion :: (Bounded c, Enum i, Ix i) => (i, i) -> Distance c i -> Dispersion c Word64 #-}
{-# SPECIALISE dispersion :: Bounded c => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> Dispersion c CUInt  #-}
{-# SPECIALISE dispersion :: Bounded c => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> Dispersion c Word   #-}
{-# SPECIALISE dispersion :: Bounded c => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> Dispersion c Word8  #-}
{-# SPECIALISE dispersion :: Bounded c => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> Dispersion c Word16 #-}
{-# SPECIALISE dispersion :: Bounded c => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> Dispersion c Word32 #-}
{-# SPECIALISE dispersion :: Bounded c => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> Dispersion c Word64 #-}
{-# SPECIALISE dispersion :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> Dispersion SymbolChangeCost CUInt  #-}
{-# SPECIALISE dispersion :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> Dispersion SymbolChangeCost Word   #-}
{-# SPECIALISE dispersion :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> Dispersion SymbolChangeCost Word8  #-}
{-# SPECIALISE dispersion :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> Dispersion SymbolChangeCost Word16 #-}
{-# SPECIALISE dispersion :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> Dispersion SymbolChangeCost Word32 #-}
{-# SPECIALISE dispersion :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> Dispersion SymbolChangeCost Word64 #-}
{-# SPECIALISE dispersion :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> Dispersion SymbolChangeCost CUInt  #-}
{-# SPECIALISE dispersion :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> Dispersion SymbolChangeCost Word   #-}
{-# SPECIALISE dispersion :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> Dispersion SymbolChangeCost Word8  #-}
{-# SPECIALISE dispersion :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> Dispersion SymbolChangeCost Word16 #-}
{-# SPECIALISE dispersion :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> Dispersion SymbolChangeCost Word32 #-}
{-# SPECIALISE dispersion :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> Dispersion SymbolChangeCost Word64 #-}
dispersion
  :: ( Bounded c
     , Enum i
     , FiniteBits b
     , Ix i
     )
  => (i, i)         -- ^ Magnitude
  -> Distance   c i -- ^ Symbol change matrix (SCM) to determine cost
  -> Dispersion c b -- ^ List of elements for of which to find the k-median and cost
dispersion ixBounds sigma xs = foldl' processRange (maxBound, zero) $ range ixBounds
  where
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
{-# INLINEABLE dispersionPairwise #-}
{-# SPECIALISE dispersionPairwise :: (Bounded c, Enum i, Ix i) => (i, i) -> Distance c i -> DispersionPairwise c CUInt  #-}
{-# SPECIALISE dispersionPairwise :: (Bounded c, Enum i, Ix i) => (i, i) -> Distance c i -> DispersionPairwise c Word   #-}
{-# SPECIALISE dispersionPairwise :: (Bounded c, Enum i, Ix i) => (i, i) -> Distance c i -> DispersionPairwise c Word8  #-}
{-# SPECIALISE dispersionPairwise :: (Bounded c, Enum i, Ix i) => (i, i) -> Distance c i -> DispersionPairwise c Word16 #-}
{-# SPECIALISE dispersionPairwise :: (Bounded c, Enum i, Ix i) => (i, i) -> Distance c i -> DispersionPairwise c Word32 #-}
{-# SPECIALISE dispersionPairwise :: (Bounded c, Enum i, Ix i) => (i, i) -> Distance c i -> DispersionPairwise c Word64 #-}
{-# SPECIALISE dispersionPairwise :: Bounded c => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionPairwise c CUInt  #-}
{-# SPECIALISE dispersionPairwise :: Bounded c => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionPairwise c Word   #-}
{-# SPECIALISE dispersionPairwise :: Bounded c => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionPairwise c Word8  #-}
{-# SPECIALISE dispersionPairwise :: Bounded c => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionPairwise c Word16 #-}
{-# SPECIALISE dispersionPairwise :: Bounded c => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionPairwise c Word32 #-}
{-# SPECIALISE dispersionPairwise :: Bounded c => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionPairwise c Word64 #-}
{-# SPECIALISE dispersionPairwise :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionPairwise SymbolChangeCost CUInt  #-}
{-# SPECIALISE dispersionPairwise :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionPairwise SymbolChangeCost Word   #-}
{-# SPECIALISE dispersionPairwise :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionPairwise SymbolChangeCost Word8  #-}
{-# SPECIALISE dispersionPairwise :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionPairwise SymbolChangeCost Word16 #-}
{-# SPECIALISE dispersionPairwise :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionPairwise SymbolChangeCost Word32 #-}
{-# SPECIALISE dispersionPairwise :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionPairwise SymbolChangeCost Word64 #-}
{-# SPECIALISE dispersionPairwise :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionPairwise SymbolChangeCost CUInt  #-}
{-# SPECIALISE dispersionPairwise :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionPairwise SymbolChangeCost Word   #-}
{-# SPECIALISE dispersionPairwise :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionPairwise SymbolChangeCost Word8  #-}
{-# SPECIALISE dispersionPairwise :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionPairwise SymbolChangeCost Word16 #-}
{-# SPECIALISE dispersionPairwise :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionPairwise SymbolChangeCost Word32 #-}
{-# SPECIALISE dispersionPairwise :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionPairwise SymbolChangeCost Word64 #-}
dispersionPairwise
  :: ( Bounded c
     , Enum i
     , FiniteBits b
     , Ix i
     )
  => (i, i)
  -> Distance           c i
  -> DispersionPairwise c b
dispersionPairwise size sigma char1 char2 = dispersion size sigma $ char1 :| [char2]


-- |
-- Calculate the median between /three/ states.
{-# INLINEABLE dispersionThreeway #-}
{-# SPECIALISE dispersionThreeway :: (Bounded c, Enum i, Ix i) => (i, i) -> Distance c i -> DispersionThreeway c CUInt  #-}
{-# SPECIALISE dispersionThreeway :: (Bounded c, Enum i, Ix i) => (i, i) -> Distance c i -> DispersionThreeway c Word   #-}
{-# SPECIALISE dispersionThreeway :: (Bounded c, Enum i, Ix i) => (i, i) -> Distance c i -> DispersionThreeway c Word8  #-}
{-# SPECIALISE dispersionThreeway :: (Bounded c, Enum i, Ix i) => (i, i) -> Distance c i -> DispersionThreeway c Word16 #-}
{-# SPECIALISE dispersionThreeway :: (Bounded c, Enum i, Ix i) => (i, i) -> Distance c i -> DispersionThreeway c Word32 #-}
{-# SPECIALISE dispersionThreeway :: (Bounded c, Enum i, Ix i) => (i, i) -> Distance c i -> DispersionThreeway c Word64 #-}
{-# SPECIALISE dispersionThreeway :: Bounded c => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionThreeway c CUInt  #-}
{-# SPECIALISE dispersionThreeway :: Bounded c => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionThreeway c Word   #-}
{-# SPECIALISE dispersionThreeway :: Bounded c => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionThreeway c Word8  #-}
{-# SPECIALISE dispersionThreeway :: Bounded c => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionThreeway c Word16 #-}
{-# SPECIALISE dispersionThreeway :: Bounded c => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionThreeway c Word32 #-}
{-# SPECIALISE dispersionThreeway :: Bounded c => (SymbolIndex, SymbolIndex) -> Distance c SymbolIndex -> DispersionThreeway c Word64 #-}
{-# SPECIALISE dispersionThreeway :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionThreeway SymbolChangeCost CUInt  #-}
{-# SPECIALISE dispersionThreeway :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionThreeway SymbolChangeCost Word   #-}
{-# SPECIALISE dispersionThreeway :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionThreeway SymbolChangeCost Word8  #-}
{-# SPECIALISE dispersionThreeway :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionThreeway SymbolChangeCost Word16 #-}
{-# SPECIALISE dispersionThreeway :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionThreeway SymbolChangeCost Word32 #-}
{-# SPECIALISE dispersionThreeway :: (Enum i, Ix i) => (i, i) -> Distance SymbolChangeCost i -> DispersionThreeway SymbolChangeCost Word64 #-}
{-# SPECIALISE dispersionThreeway :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionThreeway SymbolChangeCost CUInt  #-}
{-# SPECIALISE dispersionThreeway :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionThreeway SymbolChangeCost Word   #-}
{-# SPECIALISE dispersionThreeway :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionThreeway SymbolChangeCost Word8  #-}
{-# SPECIALISE dispersionThreeway :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionThreeway SymbolChangeCost Word16 #-}
{-# SPECIALISE dispersionThreeway :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionThreeway SymbolChangeCost Word32 #-}
{-# SPECIALISE dispersionThreeway :: (SymbolIndex, SymbolIndex) -> Distance SymbolChangeCost SymbolIndex -> DispersionThreeway SymbolChangeCost Word64 #-}
dispersionThreeway
  :: ( Bounded c
     , Enum i
     , FiniteBits b
     , Ix i
     )
  => (i, i)
  -> Distance           c i
  -> DispersionThreeway c b
dispersionThreeway size sigma char1 char2 char3 = dispersion size sigma $ char1 :| [char2, char3]


-- |
-- Gets the lowest set bit and the highest set bit in the collection.
{-# INLINEABLE getBitBounds #-}
{-# SPECIALISE getBitBounds :: Enum i => CUInt  -> Bounds CUInt  i #-}
{-# SPECIALISE getBitBounds :: Enum i => Word   -> Bounds Word   i #-}
{-# SPECIALISE getBitBounds :: Enum i => Word8  -> Bounds Word8  i #-}
{-# SPECIALISE getBitBounds :: Enum i => Word16 -> Bounds Word16 i #-}
{-# SPECIALISE getBitBounds :: Enum i => Word32 -> Bounds Word32 i #-}
{-# SPECIALISE getBitBounds :: Enum i => Word64 -> Bounds Word64 i #-}
{-# SPECIALISE getBitBounds :: CUInt  -> Bounds CUInt  SymbolIndex #-}
{-# SPECIALISE getBitBounds :: Word   -> Bounds Word   SymbolIndex #-}
{-# SPECIALISE getBitBounds :: Word8  -> Bounds Word8  SymbolIndex #-}
{-# SPECIALISE getBitBounds :: Word16 -> Bounds Word16 SymbolIndex #-}
{-# SPECIALISE getBitBounds :: Word32 -> Bounds Word32 SymbolIndex #-}
{-# SPECIALISE getBitBounds :: Word64 -> Bounds Word64 SymbolIndex #-}
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
