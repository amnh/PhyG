-----------------------------------------------------------------------------
-- |
-- Module      :  Data.TCM.Overlap
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE Strict #-}

module Data.TCM.Overlap
  ( getBitBounds
  , overlap
  , overlap2
  , overlap3
  ) where

import Data.Bits
import Data.Foldable
import Data.List.NonEmpty      (NonEmpty(..))
import Data.Semigroup
import Data.Semigroup.Foldable
import Data.Word


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
{-# SPECIALISE overlap :: (Foldable1 f, Functor f) => (Word -> Word -> Word) -> f        Word   -> (Word  , Word) #-}
{-# SPECIALISE overlap :: (Foldable1 f, Functor f) => (Word -> Word -> Word) -> f        Word8  -> (Word8 , Word) #-}
{-# SPECIALISE overlap :: (Foldable1 f, Functor f) => (Word -> Word -> Word) -> f        Word16 -> (Word16, Word) #-}
{-# SPECIALISE overlap :: (Foldable1 f, Functor f) => (Word -> Word -> Word) -> f        Word32 -> (Word32, Word) #-}
{-# SPECIALISE overlap :: (Foldable1 f, Functor f) => (Word -> Word -> Word) -> f        Word64 -> (Word64, Word) #-}
{-# SPECIALISE overlap :: FiniteBits b             => (Word -> Word -> Word) -> NonEmpty b      -> (b     , Word) #-}
{-# SPECIALISE overlap ::                             (Word -> Word -> Word) -> NonEmpty Word   -> (Word  , Word) #-}
{-# SPECIALISE overlap ::                             (Word -> Word -> Word) -> NonEmpty Word8  -> (Word8 , Word) #-}
{-# SPECIALISE overlap ::                             (Word -> Word -> Word) -> NonEmpty Word16 -> (Word16, Word) #-}
{-# SPECIALISE overlap ::                             (Word -> Word -> Word) -> NonEmpty Word32 -> (Word32, Word) #-}
{-# SPECIALISE overlap ::                             (Word -> Word -> Word) -> NonEmpty Word64 -> (Word64, Word) #-}
overlap
  :: ( FiniteBits b
     , Foldable1 f
     , Functor f
     )
  => (Word -> Word -> Word) -- ^ Symbol change matrix (SCM) to determine cost
  -> f b                    -- ^ List of elements for of which to find the k-median and cost
  -> (b, Word)              -- ^ K-median and cost
overlap sigma xs = go (upper+1) maxBound zero
  where
    zero = let wlog = getFirst $ foldMap1 First xs
           in  wlog `xor` wlog
    (lower, upper) = getBitBounds xs

    go i theCost bits | i <= lower = (bits, theCost)
    go i oldCost bits =
        let i' = i - 1
            newCost = foldl' (+) 0 $ getDistance i' <$> xs
            (minCost, bits') = case oldCost `compare` newCost of
                                 EQ -> (oldCost, bits `setBit` fromEnum i')
                                 LT -> (oldCost, bits                     )
                                 GT -> (newCost, zero `setBit` fromEnum i')
        in go i' minCost bits'

    getDistance :: FiniteBits b => Word -> b -> Word
    getDistance i b = go' (start+1) (maxBound :: Word)
      where
        (end, start) = getBitBounds $ b:|[]
        go' :: Word -> Word -> Word
        go' j a | j <= end = a
        go' j a =
          let j' = j - 1
              a' | b `testBit` fromEnum j' = min a $ sigma i j'
                 | otherwise               = a
          in  go' j' a'


-- |
-- Calculate the median between /two/ states.
{-# INLINEABLE overlap2 #-}
{-# SPECIALISE overlap2 :: (Word -> Word -> Word) -> Word   -> Word   -> (Word  , Word) #-}
{-# SPECIALISE overlap2 :: (Word -> Word -> Word) -> Word8  -> Word8  -> (Word8 , Word) #-}
{-# SPECIALISE overlap2 :: (Word -> Word -> Word) -> Word16 -> Word16 -> (Word16, Word) #-}
{-# SPECIALISE overlap2 :: (Word -> Word -> Word) -> Word32 -> Word32 -> (Word32, Word) #-}
{-# SPECIALISE overlap2 :: (Word -> Word -> Word) -> Word64 -> Word64 -> (Word64, Word) #-}
overlap2
  :: FiniteBits b
  => (Word -> Word -> Word)
  -> b
  -> b
  -> (b, Word)
overlap2 sigma char1 char2 = overlap sigma $ char1 :| [char2]


-- |
-- Calculate the median between /three/ states.
{-# INLINE     overlap3 #-}
{-# SPECIALISE overlap3 :: (Word -> Word -> Word) -> Word   -> Word   -> Word   -> (Word  , Word) #-}
{-# SPECIALISE overlap3 :: (Word -> Word -> Word) -> Word8  -> Word8  -> Word8  -> (Word8 , Word) #-}
{-# SPECIALISE overlap3 :: (Word -> Word -> Word) -> Word16 -> Word16 -> Word16 -> (Word16, Word) #-}
{-# SPECIALISE overlap3 :: (Word -> Word -> Word) -> Word32 -> Word32 -> Word32 -> (Word32, Word) #-}
{-# SPECIALISE overlap3 :: (Word -> Word -> Word) -> Word64 -> Word64 -> Word64 -> (Word64, Word) #-}
overlap3
  :: FiniteBits b
  => (Word -> Word -> Word)
  -> b
  -> b
  -> b
  -> (b, Word)
overlap3 sigma char1 char2 char3 = overlap sigma $ char1 :| [char2, char3]


-- |
-- Gets the lowest set bit and the highest set bit in the collection.
{-# INLINEABLE getBitBounds #-}
{-# SPECIALISE getBitBounds :: Foldable1  f => f        Word   -> (Word, Word) #-}
{-# SPECIALISE getBitBounds :: Foldable1  f => f        Word8  -> (Word, Word) #-}
{-# SPECIALISE getBitBounds :: Foldable1  f => f        Word16 -> (Word, Word) #-}
{-# SPECIALISE getBitBounds :: Foldable1  f => f        Word32 -> (Word, Word) #-}
{-# SPECIALISE getBitBounds :: Foldable1  f => f        Word64 -> (Word, Word) #-}
{-# SPECIALISE getBitBounds :: FiniteBits b => NonEmpty b      -> (Word, Word) #-}
{-# SPECIALISE getBitBounds ::                 NonEmpty Word   -> (Word, Word) #-}
{-# SPECIALISE getBitBounds ::                 NonEmpty Word8  -> (Word, Word) #-}
{-# SPECIALISE getBitBounds ::                 NonEmpty Word16 -> (Word, Word) #-}
{-# SPECIALISE getBitBounds ::                 NonEmpty Word32 -> (Word, Word) #-}
{-# SPECIALISE getBitBounds ::                 NonEmpty Word64 -> (Word, Word) #-}
getBitBounds
  :: ( FiniteBits b
     , Foldable1 f
     )
  => f b
  -> (Word, Word)
getBitBounds xs =
    let wlog      = getFirst $ foldMap1 First xs
        bitZero   = (wlog `xor` wlog) `setBit` 0
        bigEndian = countLeadingZeros bitZero > 0 -- Check the endianness

        (f,g) | bigEndian = (countTrailingZeros, countLeadingZeros )
              | otherwise = (countLeadingZeros , countTrailingZeros)

        lZeroes = getMin   $ foldMap1 (Min . f) xs
        uZeroes = getMin   $ foldMap1 (Min . g) xs
        lower   = toEnum     lZeroes
        upper   = toEnum   $ finiteBitSize wlog - uZeroes - 1
    in  (lower, upper)
