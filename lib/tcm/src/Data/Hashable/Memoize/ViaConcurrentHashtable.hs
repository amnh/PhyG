{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

-- {-# OPTIONS_GHC -fno-full-laziness #-}

{- |
Exposes memoization combinators. Assumes that the supplied functions are
side effect free. If this assumption is violated, undefined and unexpected
behavior may result.
-}
module Data.Hashable.Memoize.ViaConcurrentHashtable (
    memoize,
) where

import Control.DeepSeq
import Control.Monad.IO.Class (MonadIO(..))
import Data.Functor (($>))
import Data.HashTable
import Data.Hashable
import Data.IORef
import GHC.Conc.Sync
import System.IO
import System.IO.Unsafe
import Prelude hiding (lookup)


{- |
/O(1)/
-}
{-# NOINLINE memoize #-}
memoize ∷ ∀ a b m. (Hashable a, MonadIO m, NFData b) ⇒ (a → b) → m (a → b)
memoize f = liftIO $ do
    let initialSize = 2 ^ (16 ∷ Word)

    !htRef ← (newWithDefaults initialSize ∷ IO (HashTable a b)) >>= newIORef
    pure $ \k → unsafeDupablePerformIO $ do
        ht ← readIORef htRef
        result ← ht `lookup` k
        -- Here we check if the memoized value exists
        case result of
            -- If the value exists return it
            Just v → pure v
            -- If the value doesn't exist:
            Nothing →
                -- Perform the expensive calculation to determine the value
                -- associated with the key, fully evaluated.
                let v = force $ f k
                in  -- we want to perform the following modification atomically.
                    -- Insert the key-value pair into the HashTable
                    -- Then we return the value associated with the key
                    insert ht k v $> v
