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
    memoize2,
    memoize3,
) where

import Control.DeepSeq
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
memoize ∷ ∀ a b. (Eq a, Hashable a, NFData b) ⇒ (a → b) → a → b
memoize f = unsafePerformIO $ do
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


{- |
A memoizing combinator similar to 'memoize' except that it that acts on a
function of two inputs rather than one.
-}
{-# NOINLINE memoize2 #-}
memoize2 ∷ (Hashable a, Hashable b, NFData c) ⇒ (a → b → c) → a → b → c
memoize2 f =
    let f' = memoize (uncurry f)
    in  curry f'


{- |
A memoizing combinator similar to 'memoize' except that it that acts on a
function of two inputs rather than one.
-}
{-# NOINLINE memoize3 #-}
memoize3
    ∷ ( Hashable a
      , Hashable b
      , Hashable c
      , NFData d
      )
    ⇒ (a → b → c → d)
    → a
    → b
    → c
    → d
memoize3 f =
    let curry3 ∷ ((a, b, c) → t) → a → b → c → t
        curry3 g x y z = g (x, y, z)

        uncurry3 ∷ (t1 → t2 → t3 → t4) → (t1, t2, t3) → t4
        uncurry3 g (x, y, z) = g x y z

        f' = memoize (uncurry3 f)
    in  curry3 f'

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
