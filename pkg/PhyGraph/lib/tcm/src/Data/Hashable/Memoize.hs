-----------------------------------------------------------------------------
-- |
-- Module      :  Data.Hashable.Memoize
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-- Exposes memoization combinators. Assumes that the supplied functions are
-- side effect free. If this assumption is violated, undefined and unexpected
-- behavior may result.
-----------------------------------------------------------------------------

{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE CPP                 #-}
{-# LANGUAGE ScopedTypeVariables #-}

#define USING_TVAR 0
#define USING_CONC 0
#define SCREAM_ON_ACCESS 0

module Data.Hashable.Memoize
  ( memoize
  , memoize2
  , memoize3
  ) where


import Control.Concurrent.STM
import Control.Concurrent.STM.TVar
import Control.DeepSeq
import Control.Monad          (join)
--import Control.Monad.ST
import Data.Functor (($>))
import Data.Hashable
--import qualified Data.HashTable as HT
#if USING_CONC == 1
import Data.HashTable
#else
import Data.HashTable.IO
#endif
--import Data.HashTable.ST.Basic
--import qualified Data.HashTable as CHT
import Data.IORef
import Prelude           hiding (lookup)
import System.IO
import System.IO.Unsafe


-- |
-- /O(1)/
--
-- Takes a function with a hashable and equatable first argument and returns a
-- memoized function with constant time access for already computed values.
--
-- __Note:__ This does /not/ memoize recursively defined functions.
--
-- To memoize a recursively defined function, you must redefine the function's
-- recursive calls to internally call a memoized definition in a mutually recursive
-- manner.
--
-- === Example: Does /not/ memoize the recursive definitions
--
-- > fib 0 = 0
-- > fib 1 = 1
-- > fib x = fib (x-1) + fib (x-2)
--
-- >>> let memo = memoize fib in memo 10000
--
--
-- === Example: Does memoize the recursive definitions
--
-- > fibM = f
-- >   where
-- >     f 0 = 0
-- >     f 1 = 1
-- >     f x = g (x-1) + g (x-2)
-- >     g = memoize f
--
-- >>> fibM 10000
--
{-# NOINLINE memoize #-}
memoize :: forall a b. (Eq a, Hashable a, NFData b) => (a -> b) -> a -> b
memoize =
#if USING_CONC == 1
    memoize_Conc
#elif USING_TVAR == 1
    memoize_TVar
#else
    memoize_IO
#endif


#if USING_CONC == 1
{-# NOINLINE memoize_Conc #-}
memoize_Conc :: forall a b. (Eq a, Hashable a, NFData b) => (a -> b) -> a -> b
memoize_Conc f = unsafePerformIO $ do

{-
    modifyIORef memoEntries succ
    entries <- readIORef memoEntries
    if entries `mod` 50 == 0 then print entries else pure ()
-}
    let initialSize = 2 ^ (16 :: Word)

    !htRef <- (newWithDefaults initialSize :: IO (HashTable a b)) >>= newIORef
    pure $ \k -> unsafeDupablePerformIO $ do
        ht <- readIORef htRef
        result <- ht `lookup` k
        -- Here we check if the memoized value exists
        case result of
          -- If the value exists return it
          Just v  -> pure v
          -- If the value doesn't exist:
          Nothing ->
            -- Perform the expensive calculation to determine the value
            -- associated with the key, fully evaluated.
            let v = force $ f k
            -- we want to perform the following modification atomically.
            in do insert ht k v       -- Insert the key-value pair into the HashTable
                  writeIORef htRef ht -- Place the updated hashtable back in the IO-Ref
                  -- After performing the update side effects,
                  -- we return the value associated with the key
                  pure v
#else


{-# NOINLINE memoize_IO #-}
memoize_IO :: forall a b. (Eq a, Hashable a, NFData b) => (a -> b) -> a -> b
memoize_IO f = unsafePerformIO $ do

    let initialSize = 2 ^ (16 :: Word)

    !htRef <- (newSized initialSize :: IO (BasicHashTable a b)) >>= newIORef
    pure $ \k -> unsafeDupablePerformIO $ do
#if SCREAM_ON_ACCESS == 1
        -- Increment access counter
        modifyIORef memoEntries succ
        entries <- readIORef memoEntries
        if entries `mod` 50 == 0 then hPutStrLn stderr (show entries) else pure ()
#endif
        ht <- readIORef htRef
        result <- ht `lookup` k
        -- Here we check if the memoized value exists
        case result of
          -- If the value exists return it
          Just v  -> pure v
          -- If the value doesn't exist:
          Nothing ->
            -- Perform the expensive calculation to determine the value
            -- associated with the key, fully evaluated.
            let v = force $ f k
            -- we want to perform the following modification atomically.
            in do insert ht k v       -- Insert the key-value pair into the HashTable
                  writeIORef htRef ht -- Place the updated hashtable back in the IO-Ref
                  -- After performing the update side effects,
                  -- we return the value associated with the key
                  pure v


{-# NOINLINE memoize_TVar #-}
memoize_TVar :: forall a b. (Eq a, Hashable a, NFData b) => (a -> b) -> a -> b
memoize_TVar f = unsafePerformIO $ do

    let initialSize = 2 ^ (16 :: Word)

    -- Create a TVar which holds the ST state and the HashTable
    !htRef <- newTVarIO (newSized initialSize :: IO (BasicHashTable a b))
    -- This is the returned closure of a memozized f
    -- The closure captures the "mutable" reference to the hashtable above
    -- through the TVar.
    --
    -- Once the mutable hashtable reference is escaped from the IO monad,
    -- this creates a new memoized reference to f.
    -- The technique should be safe for all pure functions, probably, I think.
    pure $ \k -> unsafePerformIO $ do
        -- Read the TVar, we use IO since it is the outer monad
        -- and the documentation says that this doesn't perform a complete transaction,
        -- it just reads the current value from the TVar
        ht <- join $ readTVarIO htRef
        -- We use the HashTable to try and lookup the memoized value
        result <- ht `lookup` k
        -- Here we check if the memoized value exists
        case result of
          -- If the value exists return it
          Just v  -> pure v
          -- If the value doesn't exist:
          Nothing ->
            -- Perform the expensive calculation to determine the value
            -- associated with the key, fully evaluated.
            let v = force $ f k
            -- we want to perform the following modification atomically.
            in  atomically $
                  -- Don't use writeTVar or use a reference to the HashTable from above.
                  -- It may have been concurrently modified before reaching this point!
                  -- We *atomically* insert the new key-value pair into the existing
                  -- HashTable behind the TVar, modifying the results of the TVar.
                  modifyTVar' htRef
                    (\st -> st                 -- Get the ST state from the TVar
                        >>= (\ht' ->           -- Bind the hashtable in the state to x
                                insert ht' k v -- Insert the key-value pair into the HashTable
                                $> ht'         -- Return the HashTable as the value in ST state
                            )
                    )
                  -- After performing the update side effects,
                  -- we return the value associated with the key
                  $> v
#endif


#if SCREAM_ON_ACCESS == 1
{-# NOINLINE memoEntries #-}
memoEntries :: IORef Word
memoEntries = unsafePerformIO $ newIORef 0
#endif


-- |
-- A memoizing combinator similar to 'memoize' except that it that acts on a
-- function of two inputs rather than one.
{-# NOINLINE memoize2 #-}
memoize2 :: (Eq a, Eq b, Hashable a, Hashable b, NFData c) => (a -> b -> c) -> a -> b -> c
memoize2 f = let f' = memoize (uncurry f)
             in curry f'


-- |
-- A memoizing combinator similar to 'memoize' except that it that acts on a
-- function of two inputs rather than one.
{-# NOINLINE memoize3 #-}
memoize3
  :: ( Eq a
     , Eq b
     , Eq c
     , Hashable a
     , Hashable b
     , Hashable c
     , NFData d
     )
  => (a -> b -> c -> d)
  -> a
  -> b
  -> c
  -> d
memoize3 f = let f' = memoize (uncurry3 f)
             in curry3 f'
  where
    curry3   g  x y z  = g (x,y,z)
    uncurry3 g (x,y,z) = g x y z



{-
-- These are included for haddock generation
fib 0 = 0
fib 1 = 1
fib x = fib (x-1) + fib (x-2)

fibM :: Integer -> Integer
fibM = f
  where
    f 0 = 0
    f 1 = 1
    f x = g (x-1) + g (x-2)
    g = memoize f
-}
