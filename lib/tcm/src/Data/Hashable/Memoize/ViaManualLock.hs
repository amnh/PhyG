{-# LANGUAGE BangPatterns #-}
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
module Data.Hashable.Memoize.ViaManualLock (
    memoize,
) where

import Control.Concurrent.STM
import Control.DeepSeq
import Control.Exception (bracket, bracket_)
import Data.Functor (($>))
import Data.HashTable.IO
import Data.Hashable
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
    --    let initialSize = 2 ^ (3 :: Word)

    -- Create a TVar which holds the HashTable
    !tabRef ← newHashTableAccess initialSize
    -- This is the returned closure of a memozized f
    -- The closure captures the "mutable" reference to the hashtable above
    -- through the TVar.
    --
    -- Once the mutable hashtable reference is escaped from the IO monad,
    -- this creates a new memoized reference to f.
    -- The technique should be safe for all pure functions, probably, I think.
    pure $ \k → unsafeDupablePerformIO $ do
        readHashTableAccess tabRef k
            >>= {-# SCC memoize_Lock_CASE #-}
            \case
                Just v → {-# SCC memoize_Lock_GET #-} pure v
                Nothing →
                    {-# SCC memoize_Lock_PUT #-}
                    let v = force $ f k
                    in  --                in  (takeHashTableAccess tabRef >>= giveHashTableAccess tabRef k v) $> v
                        updateHashTableAccess tabRef k v $> v


--                in  (takeHashTableAccess tabRef >>= giveHashTableAccess_old tabRef k v) $> v

{-
-=-=-=-=-=-=-=-
    T O D O
-=-=-=-=-=-=-=-
Consider another implementation which only locks on *resize,*
permitting truly concurrent reads and writes in all but an infintesimal number of cases.
-}
data HashTableAccess k v
    = Access
        {-# UNPACK #-} Bool
        -- ^ Can Read the hash table?
        {-# UNPACK #-} (BasicHashTable k v)
        -- ^ Hash table reference


data Lockbox k v = Lockbox
    { readQueue ∷ TQueue ()
    , lockToken ∷ TMVar ()
    , memoTable ∷ TVar (HashTableAccess k v)
    }


type role HashTableAccess representational representational


type role Lockbox representational representational


{-# NOINLINE newHashTableAccess #-}
newHashTableAccess ∷ Int → IO (Lockbox k v)
newHashTableAccess size = do
    token ← newTMVarIO ()
    table ← newTVarIO =<< (Access True <$> (newSized size ∷ IO (BasicHashTable a b)))
    queue ← newTQueueIO
    pure $
        Lockbox
            { readQueue = queue
            , lockToken = token
            , memoTable = table
            }


forbidReadAccess ∷ HashTableAccess k v → HashTableAccess k v
forbidReadAccess (Access _ table) = table `seq` Access False table


permitReadAccess ∷ HashTableAccess k v → HashTableAccess k v
permitReadAccess (Access _ table) = table `seq` Access True table


{-
forbidWriteAccess :: HashTableAccess k v -> HashTableAccess k v
forbidWriteAccess (Access r _ table) = table `seq` Access r False table

permitWriteAccess :: HashTableAccess k v -> HashTableAccess k v
permitWriteAccess (Access r _ table) = table `seq` Access r True table
-}

{-# NOINLINE readHashTableAccess #-}
readHashTableAccess ∷ (Hashable k) ⇒ Lockbox k v → k → IO (Maybe v)
readHashTableAccess lockbox@(Lockbox queue _ _) key =
    {-# SCC readHashTableAccess #-}
    let request = getReadableTable lockbox
        release ∷ b → IO ()
        release = const . atomically $ readTQueue queue
    in  bracket request release $ (`lookup` key)


{-# NOINLINE updateHashTableAccess #-}
updateHashTableAccess ∷ (Hashable k) ⇒ Lockbox k v → k → v → IO ()
updateHashTableAccess lockbox key val =
    {-# SCC updateHashTableAccess #-}
    bracket_ (markWriteableTable lockbox) (freeWriteableTable lockbox) $ do
        tab ← gainWriteableTable lockbox
        insert tab key val


{-# NOINLINE getReadableTable #-}
getReadableTable ∷ Lockbox k v → IO (BasicHashTable k v)
getReadableTable (Lockbox queue _ table) = atomically $ do
    Access readable memo ← readTVar table
    check readable
    writeTQueue queue ()
    pure memo


{-# NOINLINE markWriteableTable #-}
markWriteableTable ∷ Lockbox k v → IO ()
markWriteableTable (Lockbox _ token table) = atomically $ do
    takeTMVar token
    modifyTVar' table forbidReadAccess


{-# NOINLINE gainWriteableTable #-}
gainWriteableTable ∷ Lockbox k v → IO (BasicHashTable k v)
gainWriteableTable (Lockbox queue _ table) = atomically $ do
    check =<< isEmptyTQueue queue
    Access _ memo ← readTVar table
    pure memo


{-# NOINLINE freeWriteableTable #-}
freeWriteableTable ∷ Lockbox k v → IO ()
freeWriteableTable (Lockbox _ token table) = atomically $ do
    modifyTVar' table permitReadAccess
    putTMVar token ()
