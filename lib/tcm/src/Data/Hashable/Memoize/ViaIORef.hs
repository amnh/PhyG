{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

-- {-# OPTIONS_GHC -fno-full-laziness #-}

{- |
Exposes memoization combinators. Assumes that the supplied functions are
side effect free. If this assumption is violated, undefined and unexpected
behavior may result.
-}
module Data.Hashable.Memoize.ViaIORef (
    memoize,
) where

import Control.DeepSeq
import Data.HashTable.IO
import Data.Hashable
import Data.IORef
import System.IO
import System.IO.Unsafe
import Prelude hiding (lookup)


{- |
/O(1)/

Takes a function with a hashable and equatable first argument and returns a
memoized function with constant time access for already computed values.

__Note:__ This does /not/ memoize recursively defined functions.

To memoize a recursively defined function, you must redefine the function's
recursive calls to internally call a memoized definition in a mutually recursive
manner.

=== Example: Does /not/ memoize the recursive definitions

> fib 0 = 0
> fib 1 = 1
> fib x = fib (x-1) + fib (x-2)

>>> let memo = memoize fib in memo 10000


=== Example: Does memoize the recursive definitions

> fibM = f
>   where
>     f 0 = 0
>     f 1 = 1
>     f x = g (x-1) + g (x-2)
>     g = memoize f

>>> fibM 10000
-}
{-# NOINLINE memoize #-}
memoize ∷ ∀ a b. (Hashable a, NFData b) ⇒ (a → b) → a → b
memoize f = unsafePerformIO $ do
    let initialSize = 2 ^ (16 ∷ Word)

    !htRef ← (newSized initialSize ∷ IO (BasicHashTable a b)) >>= newIORef
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
                    do
                        insert ht k v -- Insert the key-value pair into the HashTable
                        writeIORef htRef ht -- Place the updated hashtable back in the IO-Ref
                        -- After performing the update side effects,
                        -- we return the value associated with the key
                        pure v
