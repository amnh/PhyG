-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Compact.Overlap
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE Strict #-}

module Measure.Compact.Overlap
  ( overlap
  , overlap2
  , overlap3
  ) where

import Data.Bits
import Data.Foldable
import Data.List.NonEmpty      (NonEmpty(..))
import Data.Semigroup
import Data.Semigroup.Foldable
import Data.Word
import Foreign.C.Types         (CUInt)
import Measure.Matrix
import Measure.Unit.Distance
import Measure.Unit.SymbolCount
import Measure.Unit.SymbolIndex


data  Bounds b
    = Bounds
    { _lBound :: {-# UNPACK #-} SymbolIndex
    , _uBound :: {-# UNPACK #-} SymbolIndex
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
{-# INLINEABLE overlap #-}
{-# SPECIALISE overlap :: (Foldable1 f, Functor f) => SymbolCount -> SCMλ -> f        CUInt  -> (Distance,CUInt ) #-}
{-# SPECIALISE overlap :: (Foldable1 f, Functor f) => SymbolCount -> SCMλ -> f        Word   -> (Distance,Word  ) #-}
{-# SPECIALISE overlap :: (Foldable1 f, Functor f) => SymbolCount -> SCMλ -> f        Word8  -> (Distance,Word8 ) #-}
{-# SPECIALISE overlap :: (Foldable1 f, Functor f) => SymbolCount -> SCMλ -> f        Word16 -> (Distance,Word16) #-}
{-# SPECIALISE overlap :: (Foldable1 f, Functor f) => SymbolCount -> SCMλ -> f        Word32 -> (Distance,Word32) #-}
{-# SPECIALISE overlap :: (Foldable1 f, Functor f) => SymbolCount -> SCMλ -> f        Word64 -> (Distance,Word64) #-}
{-# SPECIALISE overlap :: FiniteBits b             => SymbolCount -> SCMλ -> NonEmpty b      -> (Distance,b     ) #-}
{-# SPECIALISE overlap ::                             SymbolCount -> SCMλ -> NonEmpty CUInt  -> (Distance,CUInt ) #-}
{-# SPECIALISE overlap ::                             SymbolCount -> SCMλ -> NonEmpty Word   -> (Distance,Word  ) #-}
{-# SPECIALISE overlap ::                             SymbolCount -> SCMλ -> NonEmpty Word8  -> (Distance,Word8 ) #-}
{-# SPECIALISE overlap ::                             SymbolCount -> SCMλ -> NonEmpty Word16 -> (Distance,Word16) #-}
{-# SPECIALISE overlap ::                             SymbolCount -> SCMλ -> NonEmpty Word32 -> (Distance,Word32) #-}
{-# SPECIALISE overlap ::                             SymbolCount -> SCMλ -> NonEmpty Word64 -> (Distance,Word64) #-}
overlap
  :: ( FiniteBits b
     , Foldable1 f
     , Functor f
     )
  => SymbolCount   -- ^ Alphabet size
  -> SCMλ          -- ^ Symbol change matrix (SCM) to determine cost
  -> f b           -- ^ List of elements for of which to find the k-median and cost
  -> (Distance, b) -- ^ K-median and cost
overlap size sigma xs = go (fromIntegral size) maxBound zero
  where
    withBounds = getBitBounds <$> xs
    wlog  = getFirst $ foldMap1 First xs
    zero  = wlog `xor` wlog

    go 0 theCost bits = (theCost, bits)
    go i oldCost bits =
        let i' = i - 1
            newCost = foldl' (+) 0 $ getDistance i' <$> withBounds
            (minCost, bits') = case oldCost `compare` newCost of
                                 LT -> (oldCost, bits                     )
                                 EQ -> (oldCost, bits `setBit` fromEnum i')
                                 GT -> (newCost, zero `setBit` fromEnum i')
        in go i' minCost bits'

    getDistance :: FiniteBits b => SymbolIndex -> Bounds b -> Distance
    getDistance i (Bounds lo hi b) = go' (hi+1) (maxBound :: Distance)
      where
        go' :: SymbolIndex -> Distance -> Distance
        go' j a | j <= lo = a
        go' j a =
          let j' = j - 1
              a' | b `testBit` fromEnum j' = min a $ sigma i j'
                 | otherwise               = a
          in  go' j' a'


-- |
-- Calculate the median between /two/ states.
{-# INLINEABLE overlap2 #-}
{-# SPECIALISE overlap2 :: SymbolCount -> SCMλ -> TCM2Dλ CUInt  #-}
{-# SPECIALISE overlap2 :: SymbolCount -> SCMλ -> TCM2Dλ Word   #-}
{-# SPECIALISE overlap2 :: SymbolCount -> SCMλ -> TCM2Dλ Word8  #-}
{-# SPECIALISE overlap2 :: SymbolCount -> SCMλ -> TCM2Dλ Word16 #-}
{-# SPECIALISE overlap2 :: SymbolCount -> SCMλ -> TCM2Dλ Word32 #-}
{-# SPECIALISE overlap2 :: SymbolCount -> SCMλ -> TCM2Dλ Word64 #-}
overlap2
  :: FiniteBits b
  => SymbolCount
  -> SCMλ
  -> TCM2Dλ b
overlap2 size sigma char1 char2 = overlap size sigma $ char1 :| [char2]


-- |
-- Calculate the median between /three/ states.
{-# INLINE     overlap3 #-}
{-# SPECIALISE overlap3 :: SymbolCount -> SCMλ -> TCM3Dλ CUInt  #-}
{-# SPECIALISE overlap3 :: SymbolCount -> SCMλ -> TCM3Dλ Word   #-}
{-# SPECIALISE overlap3 :: SymbolCount -> SCMλ -> TCM3Dλ Word8  #-}
{-# SPECIALISE overlap3 :: SymbolCount -> SCMλ -> TCM3Dλ Word16 #-}
{-# SPECIALISE overlap3 :: SymbolCount -> SCMλ -> TCM3Dλ Word32 #-}
{-# SPECIALISE overlap3 :: SymbolCount -> SCMλ -> TCM3Dλ Word64 #-}
overlap3
  :: FiniteBits b
  => SymbolCount
  -> SCMλ
  -> TCM3Dλ b
overlap3 size sigma char1 char2 char3 = overlap size sigma $ char1 :| [char2, char3]


-- |
-- Gets the lowest set bit and the highest set bit in the collection.
{-# INLINEABLE getBitBounds #-}
{-# SPECIALISE getBitBounds :: FiniteBits b => b      -> Bounds b      #-}
{-# SPECIALISE getBitBounds ::                 CUInt  -> Bounds CUInt  #-}
{-# SPECIALISE getBitBounds ::                 Word   -> Bounds Word   #-}
{-# SPECIALISE getBitBounds ::                 Word8  -> Bounds Word8  #-}
{-# SPECIALISE getBitBounds ::                 Word16 -> Bounds Word16 #-}
{-# SPECIALISE getBitBounds ::                 Word32 -> Bounds Word32 #-}
{-# SPECIALISE getBitBounds ::                 Word64 -> Bounds Word64 #-}
getBitBounds
  :: FiniteBits b
  => b
  -> Bounds b
getBitBounds b =
    let bitZero   = (b `xor` b) `setBit` 0
        bigEndian = countLeadingZeros bitZero > 0 -- Check the endianness

        (f,g) | bigEndian = (countTrailingZeros, countLeadingZeros )
              | otherwise = (countLeadingZeros , countTrailingZeros)

        lZeroes = f b
        uZeroes = g b
        lower   = toEnum lZeroes
        upper   = toEnum . max 0 $ finiteBitSize b - uZeroes - 1
    in  Bounds lower upper b
