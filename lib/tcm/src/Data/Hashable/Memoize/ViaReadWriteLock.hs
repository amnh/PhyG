{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE CPP #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE RoleAnnotations #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE StrictData #-}

-- {-# OPTIONS_GHC -fno-full-laziness #-}

{- |
Exposes memoization combinators. Assumes that the supplied functions are
side effect free. If this assumption is violated, undefined and unexpected
behavior may result.
-}
module Data.Hashable.Memoize.ViaReadWriteLock (
    memoize,
    memoize2,
    memoize3,
) where

import Control.Concurrent.ReadWriteVar qualified as RWLock
import Control.DeepSeq
import Data.Functor (($>))
import Data.HashTable.IO
import Data.Hashable
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

    -- Create a RWVar which holds the HashTable
    tableRef ← RWLock.new =<< (newSized initialSize ∷ IO (BasicHashTable a b))

    -- This is the returned closure of a memozized f
    -- The closure captures the "mutable" reference to the hashtable above
    -- through the TVar.
    --
    -- Once the mutable hashtable reference is escaped from the IO monad,
    -- this creates a new memoized reference to f.
    -- The technique should be safe for all pure functions, probably, I think.
    pure $ \k → unsafePerformIO $ do
        result ← RWLock.with tableRef (`lookup` k)
        case result of
            Just v → {-# SCC memoize_Lock_GET #-} pure v
            Nothing →
                {-# SCC memoize_Lock_PUT #-}
                let v = force $ f k
                in  RWLock.modify tableRef $ \t → insert t k v $> (t, v)


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