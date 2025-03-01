-----------------------------------------------------------------------------
-- |
-- Module      :  Data.Vector.NonEmpty
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE DeriveTraversable #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE RoleAnnotations #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE TypeFamilies #-}

{-# OPTIONS_GHC -fno-warn-unused-imports #-}

module Data.Vector.NonEmpty
  ( Vector(..)
  -- * Construction
  , fromNonEmpty
  , generate
  , generateM
  , singleton
  , unfoldr
  -- * Conversion
  , toVector
  , fromVector
  , unsafeFromVector
  -- * Deconstruction
  , uncons
  -- * Evaluation
  , force
  ) where


import Control.DeepSeq hiding (force)
import Control.Foldl qualified as L
import Data.Binary
import Data.Coerce
import Data.Data
import Data.Foldable
import Data.Functor.Alt
import Data.Functor.Bind
import Data.Functor.Classes
import Data.Functor.Extend
import Data.Hashable
import Data.Key
import Data.List.NonEmpty (NonEmpty((:|)))
import Data.List.NonEmpty qualified as NE
import Data.Pointed
import Data.Semigroup.Foldable
import Data.Semigroup.Traversable
import Data.Vector qualified as V
import Data.Vector.Binary ()
import Data.Vector.Instances ()
import Test.QuickCheck hiding (generate)
import Text.Read


-- |
-- A sequence of values that are repeated multiple times in contiguous blocks.
newtype Vector a = NEV { unwrap :: V.Vector a }
  deriving stock   (Data, Eq, Functor, Ord, Foldable, Traversable)
  deriving newtype ( Adjustable
                   , Applicative
                   , Apply
                   , Binary
                   , Bind
                   , Eq1
                   , Extend
                   , FoldableWithKey
                   , Hashable
                   , Indexable
                   , Keyed
                   , Lookup
                   , Monad
                   , NFData
                   , Ord1
                   , Pointed
                   , Read
                   , Read1
                   , Semigroup
                   , Zip
                   , ZipWithKey
                   )

type role Vector representational


-- |
-- Generation biases towards medium length.
instance Arbitrary a => Arbitrary (Vector a) where

    arbitrary = do
      list   <- arbitrary
      values <- case list of
                  [] -> pure <$> arbitrary
                  xs -> pure xs
      pure . NEV $ V.fromList values


instance Alt Vector where

    (<!>) = (<>)


instance Foldable1 Vector where

    {-# INLINE fold1 #-}
    fold1 = fold1 . toNonEmpty

    {-# INLINE foldMap1 #-}
    foldMap1 f = foldMap1 f . toNonEmpty

    {-# INLINE toNonEmpty #-}
    toNonEmpty = NE.fromList . toList . unwrap


instance FoldableWithKey1 Vector where

    foldMapWithKey1 f = foldMapWithKey1 f . toNonEmpty


type instance Key Vector = Int


instance Traversable1 Vector where

    {-# INLINE traverse1 #-}
    traverse1 f = fmap fromNonEmpty . traverse1 f . toNonEmpty


instance TraversableWithKey Vector where

    {-# INLINE traverseWithKey #-}
    traverseWithKey f = fmap NEV . traverseWithKey f . unwrap


instance TraversableWithKey1 Vector where

    {-# INLINE traverseWithKey1 #-}
    traverseWithKey1 f = fmap fromNonEmpty . traverseWithKey1 f . toNonEmpty


{-
instance (Read a) => Read (Vector a) where

    readPrec = parens $ prec 10 $ do
        x:xs <- readPrec
        pure . fromNonEmpty $ x:|xs

    readListPrec = readListPrecDefault
-}


instance Show a => Show (Vector a) where

    show = show . unwrap


-- |
-- /O(1)/
--
-- A synonym for 'point'.
{-# INLINE singleton #-}
singleton :: a -> Vector a
singleton = NEV . V.singleton


-- |
-- /O(n)/
--
-- Construct a 'Vector' from a non-empty structure.
{-# INLINE fromNonEmpty #-}
fromNonEmpty :: Foldable1 f => f a -> Vector a
fromNonEmpty xs =
    let !n = length xs
    in  NEV . V.fromListN n . NE.toList $ toNonEmpty xs


-- |
-- /O(n)/
--
-- Construct a vector by repeatedly applying the generator function to a seed.
-- The generator function always yields the next element and either @ Just @ the
-- new seed or 'Nothing' if there are no more elements to be generated.
--
-- > unfoldr (\n -> (n, if n == 0 then Nothing else Just (n-1))) 10
-- >  = <10,9,8,7,6,5,4,3,2,1>
{-# INLINE unfoldr #-}
unfoldr :: (b -> (a, Maybe b)) -> b -> Vector a
unfoldr f = NEV . uncurry V.fromListN . go 0
  where
--  go :: Int -> b -> (Int, [a])
    go n b =
         let (v, mb) = f b
         in  (v:) <$> maybe (n, []) (go (n + 1)) mb


-- |
-- /O(n)/
--
-- Construct a vector of the given length by applying the function to each index
generate :: Int -> (Int -> a) -> Vector a
generate n f
  | n < 1     = error $ "Called Vector.Nonempty.generate on a non-positive dimension " <> show n
  | otherwise = NEV $ V.generate n f


-- |
-- /O(n)/
--
-- Construct a vector of the given length by applying the monadic function to each index
generateM :: Monad m => Int -> (Int -> m a) -> m (Vector a)
generateM n f
  | n < 1     = error $ "Called Vector.Nonempty.generateM on a non-positive dimension " <> show n
  | otherwise = NEV <$> V.generateM n f


-- |
-- /O(1)/
--
-- Get the underlying 'V.Vector'.
toVector :: Vector a -> V.Vector a
toVector = unwrap


-- |
-- /O(1)/
--
-- Attempt to convert a 'V.Vector' to a non-empty 'Vector'.
fromVector :: V.Vector a -> Maybe (Vector a)
fromVector v
  | V.null v  = Nothing
  | otherwise = Just $ NEV v


-- |
-- /O(1)/
--
-- Attempt to convert a 'V.Vector' to a non-empty 'Vector' throwing an
-- error if the vector received is empty.
unsafeFromVector :: V.Vector a -> Vector a
unsafeFromVector v
  | V.null v  = error "NonEmpty.unsafeFromVector: empty vector"
  | otherwise = NEV v


-- | /O(n)/
--
-- 'uncons' produces both the first element of the 'Vector' and a
-- 'Vector' of the remaining elements, if any.
uncons :: Vector a -> (a, Maybe (Vector a))
uncons (NEV v) = (first, stream)
  where
    stream
      | len == 1  = Nothing
      | otherwise = Just . NEV $ V.slice 1 (len-1) v
    first = v ! 0
    len   = length v


-- |
-- Fully evaluates the 'Vector'.
force :: forall a . Vector a -> Vector a
force = coerce $ V.force @a
