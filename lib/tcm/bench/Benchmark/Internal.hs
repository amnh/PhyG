module Benchmark.Internal (
    makeMemoizedHashTableIO,
    measureMemoizedHashTable,
    measureMemoizedHashTableParallel,
) where

import Control.Concurrent.Async (mapConcurrently)
import Control.Exception
import Control.DeepSeq
import Control.Monad.State.Strict
import Criterion.Main
import Data.Bits ((.&.))
import Data.Bifunctor (first)
import Data.TCM.Overlap (overlap3)
import System.Random (StdGen, genWord32R, mkStdGen)
import Data.Vector (Vector)
import Data.Vector qualified as V
import Data.Word
import GHC.IO (evaluate)

type HashKey = Word32


type HashVal = (Word8, Word)


maxKey ∷ HashKey
maxKey = pred $ 2 ^ 24


randomKey ∷ State StdGen HashKey
randomKey = state $ first succ . genWord32R (maxKey - 1)


measureKey ∷ HashKey → HashVal
measureKey key =
    let x = α 0xFF0000
        y = α 0x00FF00
        z = α 0x0000FF
        α = fromIntegral . (key .&.)
        σ :: Word -> Word -> Word
        σ p q = max p q - min p q
    in  overlap3 8 σ x y z


seed :: Int
seed = 1029110096060097500


createHashTable ∷ Word → Word → ((HashKey → HashVal) → HashKey → HashVal) → (Vector HashKey, HashKey → HashVal)
createHashTable valueCount queryCount memoizer =
    let memo = memoizer measureKey
        rgen = mkStdGen seed
        note = V.replicateM (fromEnum queryCount) (randomKey)
        load = V.replicateM (fromEnum valueCount) (memo <$> randomKey)
        (keys, vals) = liftA2 (,) note load `evalState` rgen
    in  (keys, force vals `seq` memo)


createHashTableFour ∷ Word → Word → ((HashKey → HashVal) → HashKey → HashVal) → (Vector HashKey, Vector HashKey, Vector HashKey, Vector HashKey, HashKey → HashVal)
createHashTableFour valueCount queryCount memoizer =
    let memo = memoizer measureKey
        rgen = mkStdGen seed
        key1 = V.replicateM (fromEnum queryCount) (randomKey)
        key2 = V.replicateM (fromEnum queryCount) (randomKey)
        key3 = V.replicateM (fromEnum queryCount) (randomKey)
        key4 = V.replicateM (fromEnum queryCount) (randomKey)
        load = V.replicateM (fromEnum valueCount) (memo <$> randomKey)
        comp = (,,,,) <$> key1 <*> key2 <*> key3 <*> key4 <*> load
        vals :: Vector HashVal -> HashKey -> HashVal
        vals x = force x `seq` memo
    in  vals <$> comp `evalState` rgen


{- |
Abstract utility function for measuring the run time of the supplied
file parser on the specified input file.
-}
measureMemoizedHashTable
    ∷ Word
    → Word
    → ((HashKey → HashVal) → HashKey → HashVal)
    → String
    → Benchmark
measureMemoizedHashTable values queries memoizer name =
    let (keys, memo) = createHashTable values queries memoizer
    in  bench name $ nf (fmap memo) keys


makeMemoizedHashTableIO
  ∷ Word → Word → ((HashKey → HashVal) → IO (HashKey → HashVal)) → IO (Vector HashKey, HashKey → HashVal)
makeMemoizedHashTableIO valueCount queryCount memoizer =
    let rgen = mkStdGen seed
        note = V.replicateM (fromEnum queryCount) (randomKey)
    in  do memo <- memoizer measureKey
           let load = V.replicateM (fromEnum valueCount) (memo <$> randomKey)
           let (keys, vals) = liftA2 (,) note load `evalState` rgen
           pure (keys, force vals `seq` memo)


{- |
Abstract utility function for measuring the run time of the supplied
file parser on the specified input file.
-}
measureMemoizedHashTableParallel
    ∷ Word
    → Word
    → ((HashKey → HashVal) → HashKey → HashVal)
    → String
    → Benchmark
measureMemoizedHashTableParallel values queries memoizer name =
    let (keys1, keys2, keys3, keys4, memo) = createHashTableFour values queries memoizer
        evaluator = mapConcurrently $ \xs -> onException (evaluate . force $ memo <$> xs) (evaluate 0)
    in  bench name $ nfAppIO evaluator [ keys1, keys2, keys3, keys4 ]
