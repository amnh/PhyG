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

{-# Language BangPatterns #-}
{-# Language CPP #-}
{-# Language LambdaCase #-}
{-# LANGUAGE RoleAnnotations #-}
{-# Language ScopedTypeVariables #-}
{-# Language StrictData #-}

-- {-# OPTIONS_GHC -fno-full-laziness #-}

#define SCREAM_ON_ACCESS 0
#define USING_CONC 0
#define USING_IO   0
#define USING_LOCK 1
#define USING_SEM  0
#define USING_TVAR 0

module Data.Hashable.Memoize
  ( memoize
  , memoize2
  , memoize3
  ) where


import GHC.Conc.Sync
--import Control.Concurrent.STM
--import Control.Concurrent.STM.TMVar
--import Control.Concurrent.STM.TVar
import Control.DeepSeq
import Control.Monad          (join)
--import Control.Monad.ST
import Data.Bits (Bits(popCount))
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
#if USING_LOCK == 1
import Control.Concurrent.STM
import Control.Concurrent.STM.TMVar
import Control.Concurrent.STM.TQueue
import Control.Exception (bracket, bracket_, mask_)
import GHC.Conc (unsafeIOToSTM)
#endif
import System.IO
import System.IO.Unsafe
#if USING_SEM == 1
import Control.Exception (bracket_, mask_)
import Control.Concurrent.QSemN
import Data.Int (Int8, Int16)
#endif


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
#if   USING_CONC == 1
{-# NOINLINE memoize #-}
memoize :: forall a b. (Eq a, Hashable a, NFData b) => (a -> b) -> a -> b
memoize = memoize_Conc
#elif USING_IO   == 1
{-# NOINLINE memoize #-}
memoize :: forall a b. (Eq a, Hashable a, NFData b) => (a -> b) -> a -> b
memoize = memoize_IO
#elif USING_LOCK == 1
{-# NOINLINE memoize #-}
memoize :: forall a b. (Hashable a, NFData b) => (a -> b) -> a -> b
memoize = memoize_Lock
#elif USING_SEM == 1
{-# NOINLINE memoize #-}
memoize :: forall a b. (Eq a, Hashable a, NFData b) => (a -> b) -> a -> b
memoize = memoize_Sem
#elif USING_TVAR == 1
{-# NOINLINE memoize #-}
memoize :: forall a b. (Eq a, Hashable a, NFData b) => (a -> b) -> a -> b
memoize = memoize_TVar
#else
{-# NOINLINE memoize #-}
memoize :: forall a b. (Eq a, Hashable a, NFData b) => (a -> b) -> a -> b
memoize = error "No memoization option specified"
#endif


#if USING_SEM == 1
{-
-=-=-=-=-=-=-=-
    T O D O
-=-=-=-=-=-=-=-
Consider another implementation which only locks on *resize,*
permitting truly concurrent reads and writes in all but an infintesimal number of cases.
-}
data Lockbox k v = Lockbox
    { readQueue :: QSemN
    , memoTable :: IORef (BasicHashTable k v)
    }


type role Lockbox representational representational


-- Implement a Read/Write lock.

quota :: Int
quota = fromIntegral (maxBound :: Int16) + fromIntegral (minBound :: Int8)


quantaRead :: Int
quantaRead = 1


quantaWrite :: Int
quantaWrite = quota


operation :: Int -> QSemN -> IO a -> IO a
operation quanta semaphore = bracket_
    (waitQSemN semaphore quanta)
    (signalQSemN semaphore quanta) . mask_


safeRead :: QSemN -> IO a -> IO a
safeRead = operation quantaRead 


safeWrite :: QSemN -> IO a -> IO a
safeWrite = operation quantaWrite


{-# NOINLINE initializeSafeHashTable #-}
initializeSafeHashTable :: Int -> IO (Lockbox k v)
initializeSafeHashTable size = do
    queue <- newQSemN quota
    table <- newIORef =<< (newSized size :: IO (BasicHashTable a b))
    pure $ Lockbox
        { readQueue = queue
        , memoTable = table
        }

{-
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

{-# NOINLINE memoize_Sem #-}
memoize_Sem :: forall a b. (Hashable a, NFData b) => (a -> b) -> a -> b
memoize_Sem f = unsafePerformIO $ do

    let initialSize = 2 ^ (16 :: Word)

    -- Create a TVar which holds the HashTable
    queueSem <- newQSemN quota
    tableRef <- newIORef =<< (newSized initialSize :: IO (BasicHashTable a b))

    let accessUsing :: Int -> IO c -> IO c
        accessUsing quanta op =
          bracket_ (waitQSemN queueSem quanta) (signalQSemN queueSem quanta) $ mask_ op
    -- This is the returned closure of a memozized f
    -- The closure captures the "mutable" reference to the hashtable above
    -- through the TVar.
    --
    -- Once the mutable hashtable reference is escaped from the IO monad,
    -- this creates a new memoized reference to f.
    -- The technique should be safe for all pure functions, probably, I think.
    pure $ \k -> unsafePerformIO $ do
            result <- accessUsing quantaRead $ readIORef tableRef >>= (`lookup` k)
            case result of
                Just v  -> {-# SCC memoize_Lock_GET #-} pure v
                Nothing -> {-# SCC memoize_Lock_PUT #-}
                    let v = force $ f k
--                    in  (takeHashTableAccess tabRef >>= giveHashTableAccess tabRef k v) $> v
                    in  accessUsing quantaWrite $ readIORef tableRef >>= \t -> (insert t k v $> v)
--                    in  (takeHashTableAccess tabRef >>= giveHashTableAccess_old tabRef k v) $> v
#endif


#if USING_LOCK == 1
{-
-=-=-=-=-=-=-=-
    T O D O
-=-=-=-=-=-=-=-
Consider another implementation which only locks on *resize,*
permitting truly concurrent reads and writes in all but an infintesimal number of cases.
-}
data HashTableAccess k v = Access
    { canRead :: {-# UNPACK #-} Bool
    , hashTable :: {-# UNPACK #-} BasicHashTable k v
    }

data Lockbox k v = Lockbox
    { readQueue :: TQueue ()
    , lockToken :: TMVar ()
    , memoTable :: TVar (HashTableAccess k v)
    }


type role HashTableAccess representational representational


{-# NOINLINE newHashTableAccess #-}
newHashTableAccess :: Int -> IO (Lockbox k v)
newHashTableAccess size = do
    token <- newTMVarIO ()
    table <- newTVarIO =<< (Access True <$> (newSized size :: IO (BasicHashTable a b)))
    queue <- newTQueueIO
    pure $ Lockbox
        { readQueue = queue
        , lockToken = token
        , memoTable = table
        }


forbidReadAccess :: HashTableAccess k v -> HashTableAccess k v
forbidReadAccess (Access _ table) = table `seq` Access False table


permitReadAccess :: HashTableAccess k v -> HashTableAccess k v
permitReadAccess (Access _ table) = table `seq` Access True table


{-
forbidWriteAccess :: HashTableAccess k v -> HashTableAccess k v
forbidWriteAccess (Access r _ table) = table `seq` Access r False table


permitWriteAccess :: HashTableAccess k v -> HashTableAccess k v
permitWriteAccess (Access r _ table) = table `seq` Access r True table
-}


{-# NOINLINE readHashTableAccess #-}
readHashTableAccess :: Hashable k => Lockbox k v -> k -> IO (Maybe v)
readHashTableAccess lockbox@(Lockbox queue _ _) key = {-# SCC readHashTableAccess #-}
    let request = getReadableTable lockbox
        release = const . atomically $ readTQueue queue
    in  bracket request release $ (`lookup` key)


{-# NOINLINE updateHashTableAccess #-}
updateHashTableAccess :: Hashable k => Lockbox k v -> k -> v -> IO ()
updateHashTableAccess lockbox key val = {-# SCC updateHashTableAccess #-}
    bracket_ (markWriteableTable lockbox) (freeWriteableTable lockbox) $ do
        tab <- gainWriteableTable lockbox
        insert tab key val


{-# NOINLINE getReadableTable #-}
getReadableTable :: Lockbox k v -> IO (BasicHashTable k v)
getReadableTable (Lockbox queue token table) = atomically $ do
    Access readable memo <- readTVar table
    check readable
    writeTQueue queue ()
    pure memo


{-# NOINLINE markWriteableTable #-}
markWriteableTable :: Lockbox k v -> IO ()
markWriteableTable (Lockbox queue token table) = atomically $ do
    takeTMVar token
    modifyTVar' table forbidReadAccess


{-# NOINLINE gainWriteableTable #-}
gainWriteableTable :: Lockbox k v -> IO (BasicHashTable k v)
gainWriteableTable (Lockbox queue _ table) = atomically $ do
    check =<< isEmptyTQueue queue
    Access _ memo <- readTVar table
    pure memo


{-# NOINLINE freeWriteableTable #-}
freeWriteableTable :: Lockbox k v -> IO ()
freeWriteableTable (Lockbox _ token table) = atomically $ do
    modifyTVar' table permitReadAccess
    putTMVar token ()


{-# NOINLINE memoize_Lock #-}
memoize_Lock :: forall a b. (Hashable a, NFData b) => (a -> b) -> a -> b
memoize_Lock f = unsafePerformIO $ do

    let initialSize = 2 ^ (16 :: Word)
--    let initialSize = 2 ^ (3 :: Word)

    -- Create a TVar which holds the HashTable
    !tabRef <- newHashTableAccess initialSize
    -- This is the returned closure of a memozized f
    -- The closure captures the "mutable" reference to the hashtable above
    -- through the TVar.
    --
    -- Once the mutable hashtable reference is escaped from the IO monad,
    -- this creates a new memoized reference to f.
    -- The technique should be safe for all pure functions, probably, I think.
    pure $ \k -> unsafeDupablePerformIO $ do
        readHashTableAccess tabRef k >>= {-# SCC memoize_Lock_CASE #-} \case
              Just v  -> {-# SCC memoize_Lock_GET #-} pure v
              Nothing -> {-# SCC memoize_Lock_PUT #-}
                let v = force $ f k
--                in  (takeHashTableAccess tabRef >>= giveHashTableAccess tabRef k v) $> v
                in  updateHashTableAccess tabRef k v $> v
--                in  (takeHashTableAccess tabRef >>= giveHashTableAccess_old tabRef k v) $> v
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
#endif


#if USING_IO == 1
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
        if popCount entries == 1 || popCount == 2 then hPutStrLn stderr (show entries) else pure ()
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
#endif


#if USING_TVAR == 1
{-# NOINLINE memoize_TVar #-}
memoize_TVar :: forall a b. (Eq a, Hashable a, NFData b) => (a -> b) -> a -> b
memoize_TVar f = unsafePerformIO $ do

    let initialSize = 2 ^ (16 :: Word)

    -- Create a TVar which holds the ST state and the HashTable
    !htRef <- ((newSized initialSize) >>= newTVarIO) :: IO (TVar (BasicHashTable a b))
    -- This is the returned closure of a memozized f
    -- The closure captures the "mutable" reference to the hashtable above
    -- through the TVar.
    --
    -- Once the mutable hashtable reference is escaped from the IO monad,
    -- this creates a new memoized reference to f.
    -- The technique should be safe for all pure functions, probably, I think.
    pure $ \k -> unsafeDupablePerformIO $ do
        -- Read the TVar, we use IO since it is the outer monad
        -- and the documentation says that this doesn't perform a complete transaction,
        -- it just reads the current value from the TVar
--        ht <- readTVarIO htRef
        -- We use the HashTable to try and lookup the memoized value
        result <- readTVarIO >>= (`lookup` k)
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
memoize2 :: (Hashable a, Hashable b, NFData c) => (a -> b -> c) -> a -> b -> c
memoize2 f = let f' = memoize (uncurry f)
             in curry f'


-- |
-- A memoizing combinator similar to 'memoize' except that it that acts on a
-- function of two inputs rather than one.
{-# NOINLINE memoize3 #-}
memoize3
  :: ( Hashable a
     , Hashable b
     , Hashable c
     , NFData d
     )
  => (a -> b -> c -> d)
  -> a
  -> b
  -> c
  -> d
memoize3 f =
    let curry3 :: ((a, b, c) -> t) -> a -> b -> c -> t
        curry3 g x y z = g (x,y,z)

        uncurry3 :: (t1 -> t2 -> t3 -> t4) -> (t1, t2, t3) -> t4
        uncurry3 g (x,y,z) = g x y z

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
