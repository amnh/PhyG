-----------------------------------------------------------------------------
-- |
-- Module      :  Control.Monad.Trans.Validation
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-- The 'Validation' type's monad transformer definition.
--
-- The transformer will collect failure up until the first monadic bind,
-- after which the presence of errors on the left hand side of the bind
-- will cause the computation to short-circuit.
--
-----------------------------------------------------------------------------

{-# Language DeriveGeneric #-}
{-# Language DerivingStrategies #-}
{-# Language FlexibleContexts #-}
{-# Language FlexibleInstances #-}
{-# Language MultiParamTypeClasses #-}

module Control.Monad.Trans.Validation
    ( ValidationT (..)
    , emap
    , invalid
    ) where

import Control.Applicative
import Control.DeepSeq
import Control.Monad.Fix (MonadFix(..))
import Control.Monad.IO.Class
import Control.Monad.Trans.Class
import Control.Monad.Zip (MonadZip(..))
import Data.Bifunctor
import Data.Functor.Alt (Alt(..))
import Data.Functor.Apply (Apply(..))
import Data.Functor.Bind (Bind(..))
import Data.Functor.Classes (Eq1(..), Ord1(..), Show1(..))
import Data.String
import Data.Validation
import GHC.Generics
import Test.QuickCheck hiding (Failure, Success)


-- |
-- A monad transformer of 'Data.Validation.Validation'.
newtype ValidationT e m a
    = ValidationT { runValidationT :: m (Validation e a) }
    -- | Run the 'Control.Monad.Trans.Validation.ValidationT' monad transformer
    deriving stock (Generic)


instance Alt m => Alt (ValidationT e m)  where

    {-# INLINABLE (<!>) #-}

    (<!>) x y = ValidationT $ runValidationT x <!> runValidationT y


instance  (Monad m, Semigroup e) => Apply (ValidationT e m)  where

    {-# INLINABLE (<.>) #-}
    {-# INLINE (.>) #-}

    (<.>) f v = ValidationT $ do
        x <- runValidationT f
        case x of
            Failure e -> pure $ Failure e
            Success y -> fmap y <$> runValidationT v

    (.>) x y = ValidationT $ do
        z <- runValidationT x
        case z of
            Failure e -> pure $ Failure e
            Success _ -> runValidationT y


instance (Monad m, Semigroup e) => Applicative (ValidationT e m) where

    {-# INLINABLE (<*>) #-}
    {-# INLINABLE liftA2 #-}
    {-# INLINE (*>) #-}
    {-# INLINE (<*) #-}
    {-# INLINE pure #-}

    liftA2 f x = ValidationT . liftA2 (liftA2 f) (runValidationT x) . runValidationT

    pure  = ValidationT . pure . Success

    (<*>) = (<.>)

    (*>)  = (.>)

    (<*)  = (<.)


instance (Arbitrary a, Arbitrary e, Arbitrary1 m) => Arbitrary (ValidationT e m a) where

    {-# INLINE arbitrary #-}

    arbitrary = ValidationT <$> liftArbitrary genValidation
        where
            genValidation =
                arbitrary >>= \success -> if success then Success <$> arbitrary else Failure <$> arbitrary


instance (Apply m, Monad m, Semigroup e) => Bind (ValidationT e m) where

    {-# INLINABLE (>>-) #-}

    (>>-) v f = ValidationT $ do
        x <- runValidationT v
        case x of
            Failure e -> pure $ Failure e
            Success a -> runValidationT $ f a


instance (Eq a, Eq e, Eq1 m) => Eq (ValidationT e m a) where

    {-# INLINE (==) #-}

    (==) x = liftEq (==) (runValidationT x) . runValidationT


{-
instance Eq1 m => Eq1 (ValidationT e m) where

    {-# INLINE liftEq #-}

    liftEq f lhs = liftEq (liftEq f) (runValidationT lhs) . runValidationT
-}


instance Foldable m => Foldable (ValidationT e m) where

    {-# INLINABLE foldMap #-}

    foldMap f = foldMap (foldMap f) . runValidationT


instance Functor m => Functor (ValidationT e m) where

    {-# INLINABLE fmap #-}

    fmap f = ValidationT . fmap (fmap f) . runValidationT


instance (Monad m, NFData a, NFData e) => NFData (ValidationT e m a) where

    {-# INLINE rnf #-}

    rnf (ValidationT x) = (force <$> x) `seq` ()


instance (Monad m, Semigroup e) => Monad (ValidationT e m) where

    {-# INLINABLE (>>=) #-}
    {-# INLINE (>>) #-}
    {-# INLINE return #-}

    (>>=) v f = ValidationT $ do
        x <- runValidationT v
        case x of
            Failure e -> pure $ Failure e
            Success a -> runValidationT $ f a

    (>>)   = (*>)

    return = pure


instance (IsString e, Monad m, Semigroup e) => MonadFail (ValidationT e m) where

    {-# INLINE fail #-}

    fail = ValidationT . pure . Failure . fromString


instance (Monad m, Semigroup e) => MonadFix (ValidationT e m) where

    mfix f = let a = a >>= f in a


instance (MonadIO m, Semigroup e) => MonadIO (ValidationT e m) where

    {-# INLINE liftIO #-}

    liftIO = lift . liftIO


instance MonadTrans (ValidationT e) where

    {-# INLINE lift #-}

    lift = ValidationT . fmap Success


instance (Monad m, Semigroup e) => MonadZip (ValidationT e m) where

    {-# INLINABLE mzip #-}
    {-# INLINABLE munzip #-}
    {-# INLINE mzipWith #-}

    mzip     = liftA2 (,)

    mzipWith = liftA2

    munzip x =
        let extracting
                :: Functor f
                => ((Validation err a, Validation err b) -> Validation e c)
                -> f (Validation err (a, b))
                -> ValidationT e f c
            extracting t = ValidationT . fmap (t . vunzip)

            valuation = runValidationT x

            vunzip :: Validation err (a, b) -> (Validation err a, Validation err b)
            vunzip (Failure e     ) = (Failure e, Failure e)
            vunzip (Success (a, b)) = (Success a, Success b)
        in  (extracting fst valuation, extracting snd valuation)


instance (Ord a, Ord e, Ord1 m) => Ord (ValidationT e m a) where

    {-# INLINE compare #-}

    compare x = liftCompare compare (runValidationT x) . runValidationT


{-
instance (Ord a, Ord e, Ord1 m) => Ord (ValidationT e m a) where

    {-# INLINE compare #-}

    compare = liftCompare compare


instance Ord1 m => Ord1 (ValidationT e m) where

    {-# INLINE liftCompare #-}

    liftCompare cmp lhs = liftCompare (liftCompare cmp) (runValidationT lhs) . runValidationT
-}


instance (Applicative m, Semigroup e) => Semigroup (ValidationT e m a) where

    {-# INLINE (<>) #-}

    x <> y = ValidationT $ liftA2 (<>) (runValidationT x) (runValidationT y)


instance (Show a, Show e, Show1 m) => Show (ValidationT e m a) where

    showsPrec n = liftShowsPrec showsPrec showList n . runValidationT


instance Traversable m => Traversable (ValidationT e m) where

    {-# INLINABLE traverse #-}

    traverse f = fmap ValidationT . traverse (traverse f) . runValidationT


-- |
-- Map over the error value.
{-# INLINABLE emap #-}
emap :: Functor f => (e -> b) -> ValidationT e f a -> ValidationT b f a
emap f = ValidationT . fmap (first f) . runValidationT


-- |
-- Place an error value into the Monad transformer.
{-# INLINABLE invalid #-}
invalid :: Applicative f => e -> ValidationT e f a
invalid = ValidationT . pure . Failure
