{-# LANGUAGE Safe #-}
{-# LANGUAGE Strict #-}

module Data.TCM.Overlap (
    overlap,
    overlap2,
    overlap3,
) where

import Data.Bits
import Data.Foldable
import Data.List.NonEmpty (NonEmpty (..))
import Data.Semigroup
import Data.Semigroup.Foldable
import Data.Word
import Foreign.C.Types (CUInt)


data Bounds b = Bounds
    { _lBound ∷ {-# UNPACK #-} Word
    , _uBound ∷ {-# UNPACK #-} Word
    , _bValue ∷ b
    }


{- |
Takes one or more elements of 'FiniteBits' and a symbol change cost function
and returns a tuple of a new character, along with the cost of obtaining that
character. The return character may be (or is even likely to be) ambiguous.
Will attempt to intersect the two characters, but will union them if that is
not possible, based on the symbol change cost function.

To clarify, the return character is an intersection of all possible least-cost
combinations, so for instance, if @ char1 == A,T @ and @ char2 == G,C @, and
the two (non-overlapping) least cost pairs are A,C and T,G, then the return
value is A,C,G,T.
-}
{-# INLINEABLE overlap #-}
{-# SPECIALIZE overlap ∷ (Foldable1 f, Functor f) ⇒ Word → (Word → Word → Word) → f CUInt → (CUInt, Word) #-}
{-# SPECIALIZE overlap ∷ (Foldable1 f, Functor f) ⇒ Word → (Word → Word → Word) → f Word → (Word, Word) #-}
{-# SPECIALIZE overlap ∷ (Foldable1 f, Functor f) ⇒ Word → (Word → Word → Word) → f Word8 → (Word8, Word) #-}
{-# SPECIALIZE overlap ∷ (Foldable1 f, Functor f) ⇒ Word → (Word → Word → Word) → f Word16 → (Word16, Word) #-}
{-# SPECIALIZE overlap ∷ (Foldable1 f, Functor f) ⇒ Word → (Word → Word → Word) → f Word32 → (Word32, Word) #-}
{-# SPECIALIZE overlap ∷ (Foldable1 f, Functor f) ⇒ Word → (Word → Word → Word) → f Word64 → (Word64, Word) #-}
{-# SPECIALIZE overlap ∷ (FiniteBits b) ⇒ Word → (Word → Word → Word) → NonEmpty b → (b, Word) #-}
{-# SPECIALIZE overlap ∷ Word → (Word → Word → Word) → NonEmpty CUInt → (CUInt, Word) #-}
{-# SPECIALIZE overlap ∷ Word → (Word → Word → Word) → NonEmpty Word → (Word, Word) #-}
{-# SPECIALIZE overlap ∷ Word → (Word → Word → Word) → NonEmpty Word8 → (Word8, Word) #-}
{-# SPECIALIZE overlap ∷ Word → (Word → Word → Word) → NonEmpty Word16 → (Word16, Word) #-}
{-# SPECIALIZE overlap ∷ Word → (Word → Word → Word) → NonEmpty Word32 → (Word32, Word) #-}
{-# SPECIALIZE overlap ∷ Word → (Word → Word → Word) → NonEmpty Word64 → (Word64, Word) #-}
overlap
    ∷ ( FiniteBits b
      , Foldable1 f
      , Functor f
      )
    ⇒ Word
    -- ^ Alphabet size
    → (Word → Word → Word)
    -- ^ Symbol change matrix (SCM) to determine cost
    → f b
    -- ^ List of elements for of which to find the k-median and cost
    → (b, Word)
    -- ^ K-median and cost
overlap size sigma xs = go size maxBound zero
    where
        withBounds = getBitBounds <$> xs
        wlog = getFirst $ foldMap1 First xs
        zero = wlog `xor` wlog

        go 0 theCost bits = (bits, theCost)
        go i oldCost bits =
            let i' = i - 1
                newCost = foldl' (+) 0 $ getDistance i' <$> withBounds
                (minCost, bits') = case oldCost `compare` newCost of
                    LT → (oldCost, bits)
                    EQ → (oldCost, bits `setBit` fromEnum i')
                    GT → (newCost, zero `setBit` fromEnum i')
            in  go i' minCost bits'

        getDistance ∷ (FiniteBits b) ⇒ Word → Bounds b → Word
        getDistance i (Bounds lo hi b) = go' (hi + 1) (maxBound ∷ Word)
            where
                go' ∷ Word → Word → Word
                go' j a | j <= lo = a
                go' j a =
                    let j' = j - 1
                        a'
                            | b `testBit` fromEnum j' = min a $ sigma i j'
                            | otherwise = a
                    in  go' j' a'


{- |
Calculate the median between /two/ states.
-}
{-# INLINEABLE overlap2 #-}
{-# SPECIALIZE overlap2 ∷ Word → (Word → Word → Word) → CUInt → CUInt → (CUInt, Word) #-}
{-# SPECIALIZE overlap2 ∷ Word → (Word → Word → Word) → Word → Word → (Word, Word) #-}
{-# SPECIALIZE overlap2 ∷ Word → (Word → Word → Word) → Word8 → Word8 → (Word8, Word) #-}
{-# SPECIALIZE overlap2 ∷ Word → (Word → Word → Word) → Word16 → Word16 → (Word16, Word) #-}
{-# SPECIALIZE overlap2 ∷ Word → (Word → Word → Word) → Word32 → Word32 → (Word32, Word) #-}
{-# SPECIALIZE overlap2 ∷ Word → (Word → Word → Word) → Word64 → Word64 → (Word64, Word) #-}
overlap2
    ∷ (FiniteBits b)
    ⇒ Word
    → (Word → Word → Word)
    → b
    → b
    → (b, Word)
overlap2 size sigma char1 char2 = overlap size sigma $ char1 :| [char2]


{- |
Calculate the median between /three/ states.
-}
{-# INLINE overlap3 #-}
{-# SPECIALIZE overlap3 ∷ Word → (Word → Word → Word) → CUInt → CUInt → CUInt → (CUInt, Word) #-}
{-# SPECIALIZE overlap3 ∷ Word → (Word → Word → Word) → Word → Word → Word → (Word, Word) #-}
{-# SPECIALIZE overlap3 ∷ Word → (Word → Word → Word) → Word8 → Word8 → Word8 → (Word8, Word) #-}
{-# SPECIALIZE overlap3 ∷ Word → (Word → Word → Word) → Word16 → Word16 → Word16 → (Word16, Word) #-}
{-# SPECIALIZE overlap3 ∷ Word → (Word → Word → Word) → Word32 → Word32 → Word32 → (Word32, Word) #-}
{-# SPECIALIZE overlap3 ∷ Word → (Word → Word → Word) → Word64 → Word64 → Word64 → (Word64, Word) #-}
overlap3
    ∷ (FiniteBits b)
    ⇒ Word
    → (Word → Word → Word)
    → b
    → b
    → b
    → (b, Word)
overlap3 size sigma char1 char2 char3 = overlap size sigma $ char1 :| [char2, char3]


{- |
Gets the lowest set bit and the highest set bit in the collection.
-}
{-# INLINEABLE getBitBounds #-}
{-# SPECIALIZE getBitBounds ∷ CUInt → Bounds CUInt #-}
{-# SPECIALIZE getBitBounds ∷ Word → Bounds Word #-}
{-# SPECIALIZE getBitBounds ∷ Word8 → Bounds Word8 #-}
{-# SPECIALIZE getBitBounds ∷ Word16 → Bounds Word16 #-}
{-# SPECIALIZE getBitBounds ∷ Word32 → Bounds Word32 #-}
{-# SPECIALIZE getBitBounds ∷ Word64 → Bounds Word64 #-}
{-# SPECIALIZE getBitBounds ∷ (FiniteBits b) ⇒ b → Bounds b #-}
{-# SPECIALIZE getBitBounds ∷ CUInt → Bounds CUInt #-}
{-# SPECIALIZE getBitBounds ∷ Word → Bounds Word #-}
{-# SPECIALIZE getBitBounds ∷ Word8 → Bounds Word8 #-}
{-# SPECIALIZE getBitBounds ∷ Word16 → Bounds Word16 #-}
{-# SPECIALIZE getBitBounds ∷ Word32 → Bounds Word32 #-}
{-# SPECIALIZE getBitBounds ∷ Word64 → Bounds Word64 #-}
getBitBounds
    ∷ (FiniteBits b)
    ⇒ b
    → Bounds b
getBitBounds b =
    let bitZero = (b `xor` b) `setBit` 0
        bigEndian = countLeadingZeros bitZero > 0 -- Check the endianness
        (f, g)
            | bigEndian = (countTrailingZeros, countLeadingZeros)
            | otherwise = (countLeadingZeros, countTrailingZeros)

        lZeroes = f b
        uZeroes = g b
        lower = toEnum lZeroes
        upper = toEnum . max 0 $ finiteBitSize b - uZeroes - 1
    in  Bounds lower upper b
