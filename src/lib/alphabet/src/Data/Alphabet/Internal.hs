-----------------------------------------------------------------------------
-- |
-- Module      :  Data.Alphabet.Internal
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-- We must ensure that missing and gap are appropriately
-- code as "-" & "?", respectively, before this module is used, i.e., as output
-- from either parsers or in unification step.
--
-----------------------------------------------------------------------------

{-# LANGUAGE BangPatterns       #-}
{-# LANGUAGE DeriveAnyClass     #-}
{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE DeriveFunctor      #-}
{-# LANGUAGE DeriveGeneric      #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE FlexibleContexts   #-}
{-# LANGUAGE FlexibleInstances  #-}
{-# LANGUAGE OverloadedStrings  #-}
{-# LANGUAGE TypeFamilies       #-}

module Data.Alphabet.Internal
  ( Alphabet(..)
  , alphabetStateNames
  , alphabetSymbols
  , fromSymbols
  , fromSymbolsWithStateNames
  , fromSymbolsWithStateNamesAndTCM
  , fromSymbolsWithTCM
  , gapIndex
  , gapSymbol
  , symbolSet
  , truncateAtSymbol
  , truncateAtMaxSymbol
  , getSubsetIndex
  , getSubsetIndices
  ) where

import           Control.Arrow
import           Control.DeepSeq            (NFData)
import           Control.Monad.State.Strict
import           Data.Bifunctor             (bimap)
import           Data.Binary                (Binary)
import           Data.Bits
import           Data.Data
import           Data.Foldable
import           Data.IntSet                (IntSet)
import qualified Data.IntSet                as Int
import           Data.Key
import           Data.List                  (elemIndex, intercalate, sort)
import           Data.List.NonEmpty         (NonEmpty(..), unzip)
import qualified Data.List.NonEmpty         as NE
import qualified Data.Map                   as Map
import           Data.Matrix.NotStupid      (Matrix, matrix)
import           Data.Maybe
import           Data.Monoid
import           Data.Ord
import           Data.Semigroup.Foldable
import           Data.Set                   (Set)
import qualified Data.Set                   as Set
import           Data.String
import           Data.Text.Short            (ShortText)
import qualified Data.Vector                as V
import           Data.Vector.NonEmpty       (Vector)
import qualified Data.Vector.NonEmpty       as NEV
import           Data.Word
import           GHC.Exts                   (IsList(fromList), Item)
import           GHC.Generics               (Generic)
import           Measure.Unit.SymbolCount
import           Measure.Unit.SymbolIndex
import           Numeric.Natural
import           Prelude                    hiding (lookup, unzip, zip)
import           Test.QuickCheck
--import           Test.QuickCheck.Arbitrary.Instances ()


-- |
-- The index of the vector where the gap state is stored.
--
-- /NOTE:/ This index value is very important for many gap-related operations,
-- both internally for the 'Alphabet' module/linbrary and externally in general.
gapIndex :: SymbolIndex
gapIndex =  SymbolIndex 0


-- |
-- A collection of symbols and optional corresponding state names.
data Alphabet a =
     Alphabet
     { isSorted     :: !Bool
     , symbolVector :: {-# UNPACK #-} !(Vector a)
     , stateNames   :: [a]
     }
    deriving anyclass (Binary, NFData)
    deriving stock    (Data, Generic, Functor, Typeable)


type instance Key Alphabet = SymbolIndex


-- Newtypes for corecing and consolidation of alphabet input processing logic
newtype AlphabetInputSingle a = ASI  { toSingle ::  a    }
    deriving anyclass (NFData)
    deriving stock    (Data, Eq, Ord, Generic)


newtype AlphabetInputTuple  a = ASNI { toTuple  :: (a,a) }
    deriving anyclass (NFData)
    deriving stock    (Data, Eq, Ord, Generic)


-- |
-- -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
-- -   Supporting code and data structures:
-- -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
--

newtype UnnamedSymbol a = Unnamed  a
    deriving anyclass (NFData)
    deriving stock    (Generic)


newtype NamedSymbol   a = Named (a,a)
    deriving anyclass (NFData)
    deriving stock    (Generic)


class InternalClass a where

    gapSymbol'        :: a

    isGapSymboled     :: a -> Bool

    isMissingSymboled :: a -> Bool


instance (Ord a, IsString a) => Arbitrary (Alphabet a) where

    arbitrary = do
        n <- (arbitrary :: Gen Int) `suchThat` (\x -> 0 < x && x <= 62)
        pure . fromSymbols $ take n symbolSpace
      where
        -- We do this to simplify Alphabet generation, ensuring that there is at least one non gap symbol.
        symbolSpace = fromString . pure <$> "-" <>['0'..'9'] <> ['A'..'Z'] <> ['a'..'z'] <> "?"


-- |
-- \( \mathcal{O} \left( n * \log_2 n \right) \)
instance Ord a => Eq (Alphabet a) where

    lhs == rhs =
      (length lhs == length rhs) &&
        (  (isSorted lhs && isSorted rhs && symbolVector lhs == symbolVector rhs)
        ||  sort (toList lhs) == sort (toList rhs)
        )


instance Foldable Alphabet where

    {-# INLINE toList #-}
    toList = toList . symbolVector

    {-# INLINE foldMap #-}
    foldMap f = foldMap f . symbolVector

    {-# INLINE foldr #-}
    foldr  f e = foldr  f e . symbolVector

    {-# INLINE foldl #-}
    foldl  f e = foldl  f e . symbolVector

    {-# INLINE foldr1 #-}
    foldr1 f   = foldr1 f   . symbolVector

    {-# INLINE foldl1 #-}
    foldl1 f   = foldl1 f   . symbolVector

    {-# INLINE length #-}
    length = length . symbolVector


instance Foldable1 Alphabet where

    {-# INLINE fold1 #-}
    fold1      = fold1 . symbolVector

    {-# INLINE foldMap1 #-}
    foldMap1 f = foldMap1 f . symbolVector

    {-# INLINE toNonEmpty #-}
    toNonEmpty = toNonEmpty . symbolVector


instance FoldableWithKey Alphabet where

    {-# INLINE foldrWithKey #-}
    foldrWithKey f e = foldrWithKey (f . fromSymbolIndex) e . symbolVector

    {-# INLINE foldlWithKey #-}
    foldlWithKey f e = foldlWithKey (flip (flip f . fromSymbolIndex)) e . symbolVector


instance FoldableWithKey1 Alphabet where

    {-# INLINE foldMapWithKey1 #-}
    foldMapWithKey1 f = foldMapWithKey1 (f . fromSymbolIndex) . symbolVector


instance HasSymbolCount (Alphabet a) where

    symbolCount = SymbolCount . toEnum . length . symbolVector


instance Indexable Alphabet where

    {-# INLINE index #-}
    index a i = fromMaybe raiseError $ i `lookup` a
      where
        raiseError = error $ fold
            ["Error indexing Alphabet at location "
            , show i
            , ", valid inclusive index range is [0, "
            , show $ length a - 1
            , "]."
            ]


instance (Eq a, IsString a) => InternalClass (AlphabetInputSingle a) where

    gapSymbol'        = ASI $ fromString "-"
    isGapSymboled     = (gapSymbol' ==)
    isMissingSymboled = (ASI (fromString "?") ==)


instance (Eq a, IsString a) => InternalClass (AlphabetInputTuple a) where

    gapSymbol'                     = ASNI (fromString "-", fromString "-")
    isGapSymboled     (ASNI (x,_)) = x == fromString "-"
    isMissingSymboled (ASNI (x,_)) = x == fromString "?"


instance Lookup Alphabet where

    {-# INLINE lookup #-}
    lookup i = lookup (atSymbolIndex i) . symbolVector


instance Ord a => Ord (Alphabet a) where

    compare = comparing symbolVector


instance Show a => Show (Alphabet a) where

    show x = fold
        [ "Alphabet: {"
        , intercalate ", " $ show <$> toList x
        , "}"
        ]


-- |
-- \( \mathcal{O} \left( n \right) \)
--
-- Retrieves the state names for the symbols of the 'Alphabet'.
--
-- If the symbols of the 'Alphabet' were not given state names during
-- construction then an empty list is returned.
alphabetStateNames :: (IsList (f a), Item (f a) ~ a) => Alphabet a -> f a
alphabetStateNames = fromList . toList . stateNames


-- |
-- \( \mathcal{O} \left( n \right) \)
--
-- Retrieves the symbols of the 'Alphabet'. Synonym for 'toList'.
{-# INLINE[1] alphabetSymbols #-}
{-# RULES
"alphabetSymbols/Set"    forall (x :: Ord a => Alphabet a). alphabetSymbols x = symbolSet x
"alphabetSymbols/Vector" alphabetSymbols = (NEV.unwrap :: Vector a -> V.Vector a) . symbolVector
#-}
alphabetSymbols :: (IsList (f a), Item (f a) ~ a) => Alphabet a -> f a
alphabetSymbols = fromList . toList


-- |
-- \( \mathcal{O} \left( 1 \right) \)
--
-- Retrieves the "gap character" from the alphabet.
gapSymbol :: Alphabet a -> a
gapSymbol alphabet = alphabet ! gapIndex


-- |
-- \( \mathcal{O} \left( n \right) \) for a sorted alphabet.
--
-- \( \mathcal{O} \left( n log n \right) \) for a /non-sorted/ alphabet.
--
-- Retrieves the set of all symbols from the alphabet.
symbolSet :: Ord a => Alphabet a -> Set a
symbolSet a =
    case toList $ symbolVector a of
      []   -> mempty
      -- First element is 'gap'
      g:xs -> Set.insert g $ f xs
  where
    f | isSorted a = Set.fromDistinctAscList
      | otherwise  = Set.fromList


-- |
-- For a given subset of symbols, this function returns a positive 'Natural' number
-- in the range @[0, 2^|alphabet| - 1]@.
-- This number is the unique index of the given subset in the powerset of the alphabet.
{-# INLINE      getSubsetIndices #-}
{-# SPECIALISE  getSubsetIndices :: Alphabet String    -> Set String    -> IntSet #-}
{-# SPECIALISE  getSubsetIndices :: Alphabet ShortText -> Set ShortText -> IntSet #-}
getSubsetIndices :: Ord a => Alphabet a -> Set a -> IntSet
getSubsetIndices a s
  | isSorted a = produceSet . go low $ consumeSet s -- /O(log a + n)/, a >= n
  | otherwise  = produceSet . mo low $ consumeSet s -- /O(a)/
  where
    vec = symbolVector a
    gap = gapSymbol a
    idx = atSymbolIndex gapIndex
    low = idx + 1

    inputHadGap = Set.member gap s
    consumeSet  = Set.toAscList . Set.delete gap
    produceSet  = addGapVal . Int.fromDistinctAscList

    addGapVal
      | inputHadGap = Int.insert idx
      | otherwise   = id

    -- Faster binary search for a sorted alphabet
    go   _    []  = []
    go !lo (x:xs) =
      case withinVec vec x lo of
        Right i -> i : go (i+1) xs
        Left  i ->     go  i    xs

    -- Slower version for an unsorted alphabet
    mo _    []  = []
    mo i (x:xs)
      | i > length vec = []
      | x == (vec ! i) = i : mo 0 xs
      | otherwise      =     mo (i+1) (x:xs)


-- |
-- For a given subset of symbols, this function returns a positive 'Natural' number
-- in the range @[0, 2^|alphabet| - 1]@.
-- This number is the unique index of the given subset in the powerset of the alphabet.
{-# INLINE      getSubsetIndex #-}
{-# SPECIALISE  getSubsetIndex :: Alphabet String    -> Set String    -> Word    -> Word    #-}
{-# SPECIALISE  getSubsetIndex :: Alphabet ShortText -> Set ShortText -> Word    -> Word    #-}
{-# SPECIALISE  getSubsetIndex :: Alphabet String    -> Set String    -> Word32  -> Word32  #-}
{-# SPECIALISE  getSubsetIndex :: Alphabet ShortText -> Set ShortText -> Word32  -> Word32  #-}
{-# SPECIALISE  getSubsetIndex :: Alphabet String    -> Set String    -> Word64  -> Word64  #-}
{-# SPECIALISE  getSubsetIndex :: Alphabet ShortText -> Set ShortText -> Word64  -> Word64  #-}
{-# SPECIALISE  getSubsetIndex :: Alphabet String    -> Set String    -> Natural -> Natural #-}
{-# SPECIALISE  getSubsetIndex :: Alphabet ShortText -> Set ShortText -> Natural -> Natural #-}
getSubsetIndex :: (Bits b, Ord a) => Alphabet a -> Set a -> b -> b
getSubsetIndex a s zero
  | isSorted a = addGapVal . go zero low $ consumeSet s -- /O(log a + n)/, a >= n
  | otherwise  = addGapVal . mo zero low $ consumeSet s -- /O(a)/
  where
    vec = symbolVector a
    gap = gapSymbol a
    idx = atSymbolIndex gapIndex
    low = idx + 1

    consumeSet  = Set.toAscList . Set.delete gap
    inputHadGap = Set.member gap s

    addGapVal
      | inputHadGap = (`setBit` idx)
      | otherwise   = id

    -- Faster binary search for a sorted alphabet
    go !bits   _    []  = bits
    go !bits !lo (x:xs) =
      case withinVec vec x lo of
        Right i -> go (bits .|. bit i) (i+1) xs
        Left  i -> go  bits             i    xs

    -- Slower version for an unsorted alphabet
    mo bits  _    []    = bits
    mo bits i (x:xs)
      | i > length vec  = bits
      | x == (vec ! i)  = mo (bits .|. bit i) low     xs
      | otherwise       = mo  bits           (i+1) (x:xs)


-- |
-- \( \mathcal{O} \left( n * \log_2 n \right) \)
--
-- Constructs an 'Alphabet' from a 'Foldable' structure of symbols which are
-- 'IsString' values.
{-# INLINE[1]  fromSymbols #-}
{-# SPECIALISE fromSymbols :: Foldable t => t String    -> Alphabet String    #-}
{-# SPECIALISE fromSymbols :: Foldable t => t ShortText -> Alphabet ShortText #-}
{-# SPECIALISE fromSymbols ::               [String]    -> Alphabet String    #-}
{-# SPECIALISE fromSymbols ::               [ShortText] -> Alphabet ShortText #-}
{-# RULES
"fromSymbols/Set" forall (s :: (IsString x, Ord x) => Set x). fromSymbols s =
    let g = fromString "-"
        x = g : toList (Set.delete g s)
        v = NEV.fromNonEmpty $ NE.fromList x
    in  Alphabet True v []
#-}
fromSymbols :: (Ord a, IsString a, Foldable t) => t a -> Alphabet a
fromSymbols inputSymbols = Alphabet sorted symbols []
  where
    symbols = NEV.fromNonEmpty . fmap toSingle . alphabetPreprocessing . fmap fromSingle $ toList inputSymbols

    sorted =
        -- Coerce to a plain Vector, drop the last (gap) element
        let v = init $ toList symbols
        -- Zip each element with the next element,
        -- and assert that all pairs are less-then-equal
        in all (uncurry (<=)) . zip v $ tail v


-- |
-- \( \mathcal{O} \left( n * \log_2 n \right) \)
--
-- Constructs an 'Alphabet' from a 'Foldable' structure of symbols and
-- corresponding state names, both of which are 'IsString' values.
--
-- The input ordering is preserved.
{-# SPECIALISE fromSymbolsWithStateNames :: Foldable t => t (   String,    String)  -> Alphabet String    #-}
{-# SPECIALISE fromSymbolsWithStateNames :: Foldable t => t (ShortText, ShortText)  -> Alphabet ShortText #-}
{-# SPECIALISE fromSymbolsWithStateNames ::                [(   String,    String)] -> Alphabet String    #-}
{-# SPECIALISE fromSymbolsWithStateNames ::                [(ShortText, ShortText)] -> Alphabet ShortText #-}
fromSymbolsWithStateNames :: (Ord a, IsString a, Foldable t) => t (a, a) -> Alphabet a
fromSymbolsWithStateNames inputSymbols = Alphabet False symbols names
  where
    (symbols, names) = bimap NEV.fromNonEmpty toList . unzip . fmap toTuple . alphabetPreprocessing . fmap fromTuple $ toList inputSymbols


-- |
-- Constructs an 'Alphabet' with a corresponding TCM. Permutes TCM rows and
-- columns as the 'Alphabet' is reordered. Deletes TCM rows and columns where
-- 'Alphabet' symbols are eliminated.
--
-- If the alphabet has been permuted the corresponding TCM needs to be permuted in the same manner.
--
-- /O(n*log(n) + n^2)/
{-# SPECIALISE fromSymbolsWithTCM :: FoldableWithKey t => t String    -> Matrix b -> (Alphabet String   , Matrix b) #-}
{-# SPECIALISE fromSymbolsWithTCM :: FoldableWithKey t => t ShortText -> Matrix b -> (Alphabet ShortText, Matrix b) #-}
{-# SPECIALISE fromSymbolsWithTCM ::                       [String]    -> Matrix b -> (Alphabet String   , Matrix b) #-}
{-# SPECIALISE fromSymbolsWithTCM ::                       [ShortText] -> Matrix b -> (Alphabet ShortText, Matrix b) #-}
fromSymbolsWithTCM :: (Ord a, IsString a, FoldableWithKey t) => t a -> Matrix b -> (Alphabet a, Matrix b)
fromSymbolsWithTCM symbols originalTcm = (alphabet, permutedTcm)
  where
    (uniqueSymbols, permuted) = removeSpecialSymbolsAndDuplicates symbols
    alphabet    = Alphabet True (NEV.fromNonEmpty uniqueSymbols) []
    permutedTcm = getPermutionContext alphabet originalTcm permuted

    removeSpecialSymbolsAndDuplicates xs = (uniqueVals, permutedKeys)
      where
        count = length uniqueSymbols - 1
        uniqueVals = NE.fromList . fmap toSingle . ([gapSymbol']<>) . Set.toAscList $ Map.keysSet uniques
        permutedKeys =
          case gapMay of
            Nothing -> NE.fromList $ toList uniques <> [count]
            Just i  -> NE.fromList $ toList uniques <> [i]

        (uniques, gapMay) = foldrWithKey go (mempty, Nothing) . fmap fromSingle $ toList xs
        go i s acc@(m, g)
          | isGapSymboled     s = (m, Just i)
          | isMissingSymboled s = acc
          | otherwise           = (Map.insert s i m, g)


-- |
-- Constructs an 'Alphabet' with a corresponding TCM. Permutes TCM rows and
-- columns as the 'Alphabet' is reordered. Deletes TCM rows and columns where
-- 'Alphabet' symbols are eliminated.
--
-- If the alphabet has been permuted the corresponding TCM needs to be permuted in the same manner.
--
-- /O(n*log(n) + n^2)/
{-# SPECIALISE fromSymbolsWithStateNamesAndTCM :: FoldableWithKey t => t (   String,    String)  -> Matrix b -> (Alphabet String   , Matrix b) #-}
{-# SPECIALISE fromSymbolsWithStateNamesAndTCM :: FoldableWithKey t => t (ShortText, ShortText)  -> Matrix b -> (Alphabet ShortText, Matrix b) #-}
{-# SPECIALISE fromSymbolsWithStateNamesAndTCM ::                       [(   String,    String)] -> Matrix b -> (Alphabet String   , Matrix b) #-}
{-# SPECIALISE fromSymbolsWithStateNamesAndTCM ::                       [(ShortText, ShortText)] -> Matrix b -> (Alphabet ShortText, Matrix b) #-}
fromSymbolsWithStateNamesAndTCM :: (Ord a, IsString a, FoldableWithKey t) => t (a, a) -> Matrix b -> (Alphabet a, Matrix b)
fromSymbolsWithStateNamesAndTCM symbols originalTcm = (alphabet, permutedTcm)
  where
    (uniqueSymbols, uniqueStates, permuted) = removeSpecialSymbolsAndDuplicates symbols
    alphabet    = Alphabet True (NEV.fromNonEmpty uniqueSymbols) uniqueStates
    permutedTcm = getPermutionContext alphabet originalTcm permuted

    removeSpecialSymbolsAndDuplicates xs = (uniqueVals, uniqueNames, permutedKeys)
      where
        count = length uniqueSymbols - 1
        (uniqueVals, uniqueNames) = first NE.fromList . unzip . fmap toTuple . ([gapSymbol']<>) . Set.toAscList $ Map.keysSet uniques
        permutedKeys =
          case gapMay of
            Nothing -> NE.fromList $ toList uniques <> [count]
            Just i  -> NE.fromList $ toList uniques <> [i]

        (uniques, gapMay) = foldrWithKey go (mempty, Nothing) . fmap fromTuple $ toList xs
        go i s acc@(m, g)
          | isGapSymboled     s = (m, Just i)
          | isMissingSymboled s = acc
          | otherwise           = (Map.insert s i m, g)


getPermutionContext :: (Foldable1 f, Foldable t) => t a -> Matrix b -> f Int -> Matrix b
getPermutionContext alphabet originalTcm permuted
  | isPermuted = matrix len len f
  | otherwise  = originalTcm
  where
    len         = length alphabet
    oldOrdering = NEV.fromNonEmpty permuted
    isPermuted  = oldOrdering /= NEV.generate len id
    f (i,j)     =  originalTcm ! (oldOrdering ! i, oldOrdering ! j)


-- |
-- \( \mathcal{O} \left( n * \log_2 n \right) \)
--
-- Attempts to find the symbol in the 'Alphabet'.
-- If the symbol exists, returns an alphabet that includes all symbols occurring
-- before the supplied symbol and excludes all symbols occurring after the
-- supplied symbol. The gap character is preserved in the alphabet
-- regardless of the supplied symbol.
--
-- The resulting alphabet /includes/ the input symbol.
truncateAtSymbol :: (Ord a, IsString a) => a -> Alphabet a -> Alphabet a
truncateAtSymbol symbol alphabet =
    case elemIndex symbol $ toList alphabet of
      Nothing -> alphabet
      Just i  ->
        case alphabetStateNames alphabet of
          [] -> fromSymbols               . take (i + 1) $      alphabetSymbols alphabet
          xs -> fromSymbolsWithStateNames . take (i + 1) $ zip (alphabetSymbols alphabet) xs


-- |
-- \( \mathcal{O} \left( n * \log_2 n \right) \)
--
-- Attempts to find the maximum provided symbol in the Alphabet.
-- If the any of the provided symbols exists, returns an alphabet including all
-- the symbols occurring before the maximum provided symbol and excluding all symbols
-- occurring after the maximum supplied symbol. The gap character is
-- preserved in the alphabet regardless of the supplied symbol.
--
-- The resulting alphabet /includes/ the input symbol.
truncateAtMaxSymbol :: (Foldable t, Ord a, IsString a) => t a -> Alphabet a -> Alphabet a
truncateAtMaxSymbol symbols alphabet =
    case maxIndex of
      Nothing -> alphabet
      Just i  ->
        case alphabetStateNames alphabet of
          [] -> fromSymbols               . take (atSymbolIndex i + 1) $      alphabetSymbols alphabet
          xs -> fromSymbolsWithStateNames . take (atSymbolIndex i + 1) $ zip (alphabetSymbols alphabet) xs
  where
    maxIndex = foldlWithKey' f Nothing alphabet
    f e k v
      | v `notElem` symbols = e
      | otherwise =
        case e of
          Nothing -> Just k
          Just  i -> Just $ max k i


-- |
-- \( \mathcal{O} \left( n * \log_2 n \right) \)
alphabetPreprocessing :: (Ord a, InternalClass a, Foldable t) => t a -> NonEmpty a
alphabetPreprocessing = prependGapSymbol . sort . removeSpecialSymbolsAndDuplicates . toList
  where
    prependGapSymbol xs =
        case xs of
          []   -> gapSymbol':|[]
          y:ys -> gapSymbol':|(y:ys)


    removeSpecialSymbolsAndDuplicates = (`evalState` mempty) . filterM f
      where
        f x
          | isGapSymboled     x = pure False
          | isMissingSymboled x = pure False
          | otherwise           = do
              seenSet <- get
              _       <- put $ x `Set.insert` seenSet
              pure $ x `notElem` seenSet


fromSingle :: a -> AlphabetInputSingle a
fromSingle = ASI


fromTuple :: (a, a) -> AlphabetInputTuple a
fromTuple  = ASNI


{-# INLINE withinVec #-}
{-# SPECIALISE withinVec :: Vector String    -> String    -> Int -> Either Int Int #-}
-- {-# SPECIALISE withinVec :: Vector ShortText -> ShortText -> Int -> Either Int Int #-}
withinVec :: Ord a => Vector a -> a -> Int -> Either Int Int
withinVec v e m
  | e == gap  = Right $ atSymbolIndex gapIndex
  | otherwise = go m  $ length v - 1
  where
    gap = v ! atSymbolIndex gapIndex
    -- Perform a binary search on the unboxed vector
    -- to determine if a symbol is present.
    --
    -- Equally fast, and uses less memory than a Set.
    {-# INLINE go #-}
    go !lo !hi
      | lo > hi   = Left hi
      | otherwise = let !md = (hi + lo) `div` 2
                        !z  = v ! md
                    in  case z `compare` e of
                          EQ -> Right md
                          LT -> go    (md + 1) hi
                          GT -> go lo (md - 1)


fromSymbolIndex :: Int -> SymbolIndex
fromSymbolIndex = SymbolIndex . toEnum
