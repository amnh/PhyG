-----------------------------------------------------------------------------
-- |
-- Module      :  DirectOptimization.Pairwise.Test
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-- Test suite for dynamic characters
--
-----------------------------------------------------------------------------

{-# LANGUAGE DerivingStrategies         #-}
{-# LANGUAGE FlexibleContexts           #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE TypeFamilies               #-}

module Test.QuickCheck.Instances.DynamicCharacter
  ( DNA(..)
  , DyadDNA(..)
  , SnippedDNA(..)
  , Nucleotide(..)
  , buildDNA
  , nucleotideAlphabet
  , nucleotideGap
  ) where

import           Bio.DynamicCharacter
import           Control.Applicative
import           Data.Alphabet
import           Data.Alphabet.Codec
import           Data.Alphabet.IUPAC       (iupacToDna)
import qualified Data.Bimap                as BM
import           Data.Bits
import           Data.Foldable
import           Data.List                 (intercalate)
import           Data.List.NonEmpty        (NonEmpty(..))
import qualified Data.List.NonEmpty        as NE
import           Data.MetricRepresentation
import           Data.Set                  (Set)
import qualified Data.Set                  as Set
import qualified Data.Vector               as V
import           Data.Vector.Generic       (Vector)
import qualified Data.Vector.Generic       as GV
import qualified Data.Vector.Storable      as SV
import           Foreign.C.Types           (CUInt(..))
import           Test.Tasty.QuickCheck     hiding ((.&.))


data DyadDNA = !DNA :×: !DNA
    deriving stock (Eq, Ord)


newtype DNA = DNA SlimDynamicCharacter
    deriving newtype (Eq, Ord)


newtype SnippedDNA = Snip { getSnippedDNA :: SV.Vector SlimState }
    deriving newtype (Eq, Ord)


newtype Nucleotide = N { getNucleotide :: SlimState }
    deriving newtype (Eq, Ord)


newtype NucleotideContext a = NC { getContext :: (a, a, a) }
    deriving newtype (Eq, Ord, Show)


newtype TestDynamicCharacter v a = DC (OpenDynamicCharacter v a)
    deriving newtype (Eq, Ord)


instance Arbitrary DyadDNA where

    arbitrary = do
        lhs <- arbitrary :: Gen DNA
        -- Define a predicate to ensure that the Ukkonen calls
        -- actually use Ukkonen's method rather than defaulting
        -- to the full direction matrix with swapping cost vector.
        let predicate =
              case sizeOfDNA lhs of
              -- If the character is empty, don;t worry about Ukkonen
              0 -> const True
              -- Otherwise, ensure that the size disparity of DNA
              -- is not large enough to trigger "No Gain From Ukkonen."
              n -> \x -> let m = sizeOfDNA x in 2 * max m n < 3 * min m n
        rhs <- arbitrary `suchThat` predicate
        pure $ lhs :×: rhs

    shrink (lhs :×: rhs) = uncurry (:×:) <$> shrink (lhs, rhs)


instance Arbitrary DNA where

    arbitrary = do
        DC (lc,mc,rc) <- arbitrary :: Gen (TestDynamicCharacter V.Vector Nucleotide)
        let lc' = toList $ getNucleotide <$> lc
        let mc' = toList $ getNucleotide <$> mc
        let rc' = toList $ getNucleotide <$> rc
        pure $ DNA ( SV.fromList lc', SV.fromList mc', SV.fromList rc' )

    shrink =
        let makeShrinkable
              :: DNA
              -> [NucleotideContext SlimState]
            makeShrinkable (DNA (x,y,z)) = NC <$> zip3 (SV.toList x) (SV.toList y) (SV.toList z)

            fromShrinkable
              :: [NucleotideContext SlimState]
              -> DNA
            fromShrinkable = DNA . tripleVector . unzip3 . fmap getContext

            tripleVector
              :: ([SlimState], [SlimState], [SlimState])
              -> (SV.Vector SlimState, SV.Vector SlimState, SV.Vector SlimState)
            tripleVector (x,y,z) = (SV.fromList x, SV.fromList y, SV.fromList z)

        in  fmap fromShrinkable . shrink . makeShrinkable


instance (Arbitrary a, FiniteBits a, Vector v a) => Arbitrary (TestDynamicCharacter v a) where

    arbitrary = do
        makeFilled <- chooseEnum (0, 9 :: Word)
        case makeFilled of
          0 -> pure $ DC (GV.empty, GV.empty, GV.empty)
          _ -> do
            randLength <- chooseInt (8, 24 :: Int)
            (lc,mc,rc) <- unzip3 . fmap getContext <$> vectorOf randLength arbitrary
            pure $ DC (GV.fromList lc, GV.fromList mc, GV.fromList rc)


instance (Arbitrary a, Bits a, Eq a, FiniteBits a) => Arbitrary (NucleotideContext a) where

    -- 50% »»» 'A' Alignment
    -- 25% »»» 'D' Deletion
    -- 25% »»» 'I' Insertion
    arbitrary = do
        ctx <- chooseEnum (0, 4 :: Word)
        tmp <- arbitrary
        let nil = tmp `xor` tmp  -- Create a value with zero set bits
        let gap = nil `setBit` 0 -- Create a gap by setting just least significant bit
        let constVoid = elements [nil]
        (x,y) <- case ctx of
                   0 -> liftA2 (,) arbitrary constVoid -- Insert
                   1 -> liftA2 (,) constVoid arbitrary -- Delete
                   _ -> liftA2 (,) arbitrary arbitrary -- Align
        let med = generateMedian nil gap x y
        pure $ NC (x, med, y)
      where
        generateMedian n g x y =
          let tcm = retreivePairwiseTCM discreteMetric
              f v | v == n    = g
                  | otherwise = v
          in  fst $ tcm (f x) (f y)


instance Arbitrary Nucleotide where

    -- 10% »»» '-' Gap
    -- 15% »»» 'A' Adenine
    -- 15% »»» 'C' Cytosine
    -- 15% »»» 'G' Guanine
    -- 15% »»» 'T' Thymine
    --  5% »»» '?' Anything/Missing/Unknown
    -- 25% »»»     Other Random Ambiguity
    arbitrary = do
        n <- chooseEnum (0, 3 :: Word)
        case n of
          0 -> generateUncommon
          _ -> generateCommon

      where
        makeNuke :: Foldable f => f String -> Gen Nucleotide
        makeNuke = pure . N . encodeState nucleotideAlphabet (const 0)
        gap      = makeNuke ["-"]
        adenine  = makeNuke ["A"]
        cytosine = makeNuke ["C"]
        guanine  = makeNuke ["G"]
        thymine  = makeNuke ["T"]
        unknown  = makeNuke symbols
        symbols  = alphabetSymbols nucleotideAlphabet :: Set String

        -- 13.3̅% (2/15) »»» '-' Gap
        -- 20.0% (3/15) »»» 'A' Adenine
        -- 20.0% (3/15) »»» 'C' Cytosine
        -- 20.0% (3/15) »»» 'G' Guanine
        -- 20.0% (3/15) »»» 'T' Thymine
        --  6.6̅% (1/15) »»» '?' Anything/Missing/Unknown
        generateCommon :: Gen Nucleotide
        generateCommon = do
            m <- chooseEnum (0, 4 :: Word)
            case m of
              0 -> adenine
              1 -> cytosine
              2 -> guanine
              3 -> thymine
              _ -> do
                n <- chooseEnum (0, 2 :: Word)
                case n of
                  0 -> unknown
                  _ -> gap

        -- Generate an ambiguity group other than the total ambiguity
        generateUncommon :: Gen Nucleotide
        generateUncommon = do
            -- We want to choose the remaining ambiguous states
            -- uniformly at random!
            -- To do so we will sample a random variable "n"
            -- with the uniform distribution between 0 and 4, inclusively.
            n <- chooseEnum (0, 4 :: Word)
            -- We observe the distribution of our remaining
            -- ambiguous states, by considering set bits:
            -- 5 choose 2 = 10
            -- 5 choose 3 = 10
            -- 5 choose 4 =  5
            -- So we divide the random variable "n's" outcomes as so:
            -- 4 & 3 -> 5 choose 2
            -- 2 & 1 -> 5 choose 3
            -- 0     -> 5 choose 4
            w <- elements $ toList symbols
            let s' = Set.delete w symbols
            x <- elements . toList $ s'
            if n >= 3
            then makeNuke [w, x]
            else do
              let s'' = Set.delete x s'
              y <- elements . toList $ s''
              if n >= 1
              then makeNuke [w, x, y]
              else do
                let s''' = Set.delete y s''
                z <- elements . toList $ s'''
                makeNuke [w, x, y, z]


-- Helper function for Bits instance of the Nucleotide type
binOp :: (SlimState -> SlimState -> SlimState) -> Nucleotide -> Nucleotide -> Nucleotide
binOp op (N x) (N y) = N $ x `op` y


idxOp :: (SlimState -> Int -> SlimState) -> Nucleotide -> Int -> Nucleotide
idxOp op n@(N x) i
  | i > length nucleotideAlphabet = n
  | otherwise = N $ x `op` i


instance Bits Nucleotide where

    {-# INLINE (.&.) #-}
    (.&.) = binOp (.&.)

    {-# INLINE (.|.) #-}
    (.|.) = binOp (.|.)

    {-# INLINE xor #-}
    xor = binOp xor

    {-# INLINE complement #-}
    complement = N . complement . getNucleotide

    {-# INLINE zeroBits #-}
    zeroBits = N 0

    {-# INLINE bit #-}
    bit i
      | i > length nucleotideAlphabet = N 0
      | otherwise = N $ bit i

    {-# INLINE clearBit #-}
    clearBit = idxOp clearBit

    {-# INLINE setBit #-}
    setBit = idxOp setBit

    {-# INLINE testBit #-}
    testBit (N x) i
      | i > length nucleotideAlphabet = False
      | otherwise = x `testBit` i

    bitSize = const $ length nucleotideAlphabet

    {-# INLINE bitSizeMaybe #-}
    bitSizeMaybe = const . Just $ length nucleotideAlphabet

    {-# INLINE isSigned #-}
    isSigned = const False

    {-# INLINE shiftL #-}
    shiftL (N x) k =
        let mask = bit (length nucleotideAlphabet) - 1
        in  N $ mask .&. (x `shiftL` k)

    {-# INLINE shiftR #-}
    shiftR (N x) = N . shiftR x

    {-# INLINE rotateL #-}
    rotateL (N x) k =
        let rot     = k `mod` length nucleotideAlphabet
            maskTop = (bit rot - 1) `shiftL` (length nucleotideAlphabet - rot)
            maskAll = bit (length nucleotideAlphabet) - 1
            mask    = maskTop .&. maskAll
            top = (x .&. mask) `shiftR` (length nucleotideAlphabet - rot)
            bot =  x           `shiftL` rot
        in  N $ top .|. bot

    {-# INLINE rotateR #-}
    rotateR (N x) k =
        let rot  = k `mod` length nucleotideAlphabet
            mask = bit rot - 1
            top  = (x .&. mask) `shiftL` (length nucleotideAlphabet - rot)
            bot  =  x           `shiftR` rot
        in  N $ top .|. bot

    {-# INLINE popCount #-}
    popCount = popCount . getNucleotide


instance FiniteBits Nucleotide where

    {-# INLINE finiteBitSize #-}
    finiteBitSize = const $ length nucleotideAlphabet

    {-# INLINE countTrailingZeros #-}
    countTrailingZeros = countTrailingZeros . getNucleotide

    {-# INLINE countLeadingZeros #-}
    countLeadingZeros = countLeadingZeros . getNucleotide


instance Show DyadDNA where

    show (lhs :×: rhs) = unlines
        [ ""
        , boxTop
        , boxDNA lhs
        , boxMid
        , boxDNA rhs
        , boxBot
        ]
      where
        padDNA = ("  │" <>) . (<> "│")
        boxDNA = intercalate "\n" . fmap padDNA . lines . show
        boxTop = "  ┌" <> replicate lhsLen '─' <> "┐"
        boxBot = "  └" <> replicate rhsLen '─' <> "┘"
        boxMid =
          let (mid, end) =
                case compare lhsLen rhsLen of
                  LT -> ("┴","┐")
                  EQ -> ( "","┤")
                  GT -> ("┬","┘")
          in  fold [ "  ├", replicate minLen '─', mid, replicate remLen '─', end ]

        remLen = maxLen - minLen - 1
        maxLen = max lhsLen rhsLen
        minLen = min lhsLen rhsLen
        lhsLen = strLen lhs
        rhsLen = strLen rhs
        strLen x =
          case sizeOfDNA x of
            0 -> length $ show x
            n -> n


instance Show DNA where

    show dna@(DNA (lc,mc,rc)) =
        case sizeOfDNA dna of
          0 -> "<∅,∅,∅>"
          _ ->
            let render = show . Snip
            in  unlines
                [ render lc
                , render mc
                , render rc
                ]


instance Show SnippedDNA where

    show = foldMap (show . N) . SV.toList . getSnippedDNA


instance Show Nucleotide where

    show (N state) =
        case symbols of
          [] | state == 0 -> " "
          []              -> "▒"
          x:xs            -> pure . head . iupac $ x:|xs
      where
         symbols = decodeState nucleotideAlphabet state :: [String]
         iupac :: NonEmpty String -> String
         iupac s =
           let value = BM.lookupR s iupacToDna :: Maybe (NonEmpty String)
           in  maybe " " NE.head value


sizeOfDNA :: DNA -> Int
sizeOfDNA (DNA (_,mc,_)) = GV.length mc


nucleotideAlphabet :: Alphabet String
nucleotideAlphabet = fromSymbols . NE.fromList $ pure <$> "-ACGT"


nucleotideGap :: Nucleotide
nucleotideGap = N $ encodeState nucleotideAlphabet (const 0) ["-"]


buildDNA :: String -> String -> String -> DNA
buildDNA lStr mStr rStr = DNA (lVec, mVec, rVec)
  where
    lBit = encode <$> lStr
    mBit = encode <$> mStr
    rBit = encode <$> rStr

    lVec = SV.fromList lBit
    mVec = SV.fromList mBit
    rVec = SV.fromList rBit

    encode :: Char -> SlimState
    encode ' ' = 0
    encode chr = valMay . iupac $ [chr] :| []

    valMay = maybe 0 (encodeState nucleotideAlphabet (const 0))

    iupac :: NonEmpty String -> Maybe (NonEmpty String)
    iupac s = BM.lookup s iupacToDna
