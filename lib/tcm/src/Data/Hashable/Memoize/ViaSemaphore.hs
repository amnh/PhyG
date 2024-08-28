{-# LANGUAGE RoleAnnotations #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE StrictData #-}

-- {-# OPTIONS_GHC -fno-full-laziness #-}

{- |
Exposes memoization combinators. Assumes that the supplied functions are
side effect free. If this assumption is violated, undefined and unexpected
behavior may result.
-}
module Data.Hashable.Memoize.ViaSemaphore (
    memoize,
) where

import Control.Concurrent.QSemN
import Control.DeepSeq
import Control.Exception (bracket_, mask_)
import Data.Functor (($>))
import Data.HashTable.IO
import Data.Hashable
import Data.IORef
import Data.Int (Int16, Int8)
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

    -- Create a TVar which holds the HashTable
    queueSem ← newQSemN quota
    tableRef ← newIORef =<< (newSized initialSize ∷ IO (BasicHashTable a b))

    let accessUsing ∷ Int → IO c → IO c
        accessUsing quanta op =
            bracket_ (waitQSemN queueSem quanta) (signalQSemN queueSem quanta) $ mask_ op
    -- This is the returned closure of a memozized f
    -- The closure captures the "mutable" reference to the hashtable above
    -- through the TVar.
    --
    -- Once the mutable hashtable reference is escaped from the IO monad,
    -- this creates a new memoized reference to f.
    -- The technique should be safe for all pure functions, probably, I think.
    pure $ \k → unsafePerformIO $ do
        result ← accessUsing quantaRead $ readIORef tableRef >>= (`lookup` k)
        case result of
            Just v → {-# SCC memoize_Lock_GET #-} pure v
            Nothing →
                {-# SCC memoize_Lock_PUT #-}
                let v = force $ f k
                in  --                    in  (takeHashTableAccess tabRef >>= giveHashTableAccess tabRef k v) $> v
                    accessUsing quantaWrite $ readIORef tableRef >>= \t → (insert t k v $> v)


--                    in  (takeHashTableAccess tabRef >>= giveHashTableAccess_old tabRef k v) $> v

{-
{-
-=-=-=-=-=-=-=-
    T O D O
-=-=-=-=-=-=-=-
Consider another implementation which only locks on *resize,*
permitting truly concurrent reads and writes in all but an infintesimal number of cases.
-}
data Lockbox k v = Lockbox
    {-# UNPACK #-} QSemN -- ^ The "Read Queue" of the nlock box
    {-# UNPACK #-} (IORef (BasicHashTable k v)) -- ^ The Hashtable reference

type role Lockbox representational representational
-}

-- Implement a Read/Write lock.

quota ∷ Int
quota = fromIntegral (maxBound ∷ Int16) + fromIntegral (minBound ∷ Int8)


quantaRead ∷ Int
quantaRead = 1


quantaWrite ∷ Int
quantaWrite = quota


{-
operation :: Int -> QSemN -> IO a -> IO a
operation quanta semaphore =
    bracket_
        (waitQSemN semaphore quanta)
        (signalQSemN semaphore quanta)
        . mask_

safeRead :: QSemN -> IO a -> IO a
safeRead = operation quantaRead

safeWrite :: QSemN -> IO a -> IO a
safeWrite = operation quantaWrite

{-# NOINLINE initializeSafeHashTable #-}
initializeSafeHashTable :: Int -> IO (Lockbox k v)
initializeSafeHashTable size = do
    queue <- newQSemN quota
    table <- newIORef =<< (newSized size :: IO (BasicHashTable a b))
    pure $ Lockbox queue table

{-# NOINLINE safelyReadValue #-}
safelyReadValue :: Hashable k => Lockbox k v -> k -> IO (Maybe v)
safelyReadValue lockbox@(Lockbox queue table) key = {-# SCC readHashTableAccess #-} safeRead queue $ do
    tab <- {-# SCC readHashTableAccess_Atomic_Block #-} readTVarIO table
    tab `lookup` key

{-# NOINLINE safelyWriteValue #-}
safelyWriteValue :: Hashable k => Lockbox k v -> k -> v -> IO ()
safelyWriteValue lockbox@(Lockbox queue table) key val = {-# SCC updateHashTableAccess #-} safeWrite queue $ do
    tab <- {-# SCC readHashTableAccess_Atomic_Block #-} readTIORef table
    atomicModifyIORef' $ insert tab key val
-}
