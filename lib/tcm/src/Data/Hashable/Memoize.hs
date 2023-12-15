{-# LANGUAGE CPP #-}

-- {-# OPTIONS_GHC -fno-full-laziness #-}

{- |
Exposes memoization combinators. Assumes that the supplied functions are
side effect free. If this assumption is violated, undefined and unexpected
behavior may result.
-}
module Data.Hashable.Memoize (
    memoize,
    memoize2,
    memoize3,
) where

import Control.DeepSeq (NFData)
import Control.Monad.IO.Class (MonadIO(..))
import Data.Hashable (Hashable)
#if defined (Memoize_Via_ConcurrentHashtable)
import Data.Hashable.Memoize.ViaConcurrentHashtable qualified as Memo (memoize)
#elif defined (Memoize_Via_IORef)
import Data.Hashable.Memoize.ViaIORef qualified as Memo (memoize)
#elif defined (Memoize_Via_ManualLock)
import Data.Hashable.Memoize.ViaManualLock qualified as Memo (memoize)
#elif defined (Memoize_Via_ReadWriteLock)
import Data.Hashable.Memoize.ViaReadWriteLock qualified as Memo (memoize)
#elif defined (Memoize_Via_Semaphore)
import Data.Hashable.Memoize.ViaSemaphore qualified as Memo (memoize)
#elif defined (Memoize_Via_TVar)
import Data.Hashable.Memoize.ViaTVar qualified as Memo (memoize)
#endif


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
memoize ∷ ∀ a b m. (Eq a, Hashable a, MonadIO m, NFData b) ⇒ (a → b) → m (a → b)
memoize = Memo.memoize


{- |
A memoizing combinator similar to 'memoize' except that it that acts on a
function of two inputs rather than one.
-}
{-# NOINLINE memoize2 #-}
memoize2 ∷ (Hashable a, Hashable b, MonadIO m, NFData c) ⇒ (a → b → c) → m (a → b → c)
memoize2 f =
    let f' = memoize (uncurry f)
    in  curry <$> f'


{- |
A memoizing combinator similar to 'memoize' except that it that acts on a
function of two inputs rather than one.
-}
{-# NOINLINE memoize3 #-}
memoize3
    ∷ ( Hashable a
      , Hashable b
      , Hashable c
      , MonadIO m
      , NFData d
      )
    ⇒ (a → b → c → d)
    → m (a → b → c → d)
memoize3 f =
    let curry3 ∷ ((a, b, c) → t) → a → b → c → t
        curry3 g x y z = g (x, y, z)

        uncurry3 ∷ (t1 → t2 → t3 → t4) → (t1, t2, t3) → t4
        uncurry3 g (x, y, z) = g x y z

        f' = memoize (uncurry3 f)
    in  curry3 <$> f'

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
