{- |
Benchmarks for different memoized hashtable implmentations.
-}
module Main (main) where


import Benchmark.Internal
import Criterion.Main
import Data.Hashable.Memoize.ViaConcurrentHashtable qualified as Conc
import Data.Hashable.Memoize.ViaIORef qualified as IORf
import Data.Hashable.Memoize.ViaManualLock qualified as Lock
import Data.Hashable.Memoize.ViaReadWriteLock qualified as RWLk
import Data.Hashable.Memoize.ViaSemaphore qualified as Semp
import Data.Hashable.Memoize.ViaTVar qualified as TVar
import Data.Word


-- |
-- Entry point for the run time performance benchmark suite /all/ the file parsers.
main :: IO ()
main = defaultMain
    [ benchAsSequential
    , benchAsParallel
    ]



benchAsSequential :: Benchmark
benchAsSequential =
    let basicMeasure :: ((Word32 -> (Word8, Word)) -> Word32 -> (Word8, Word)) -> String -> Benchmark
        basicMeasure = measureMemoizedHashTable (2 ^ 10) (2 ^ 20)
    in  bgroup "Sequential"
            [ basicMeasure Lock.memoize "Custom locking definition"
            , basicMeasure IORf.memoize "Manual access through - IORef"
            , basicMeasure TVar.memoize "Manual access through - TVar"
            , basicMeasure Semp.memoize "Manual access through - Semaphore"
            , basicMeasure Conc.memoize "Package `concurrent-hashtables` Hash-table"
            , basicMeasure RWLk.memoize "Package `concurrent-extra`      Read/Write Lock"
            ]


benchAsParallel :: Benchmark
benchAsParallel = 
    let basicMeasure :: ((Word32 -> (Word8, Word)) -> Word32 -> (Word8, Word)) -> String -> Benchmark
        basicMeasure = measureMemoizedHashTableParallel (2 ^ 10) (2 ^ 20)
    in  bgroup "Parallel"
            [ basicMeasure Conc.memoize "Package `concurrent-hashtables` Hash-table"
            , basicMeasure RWLk.memoize "Package `concurrent-extra`      Read/Write Lock"
--            , basicMeasure TVar.memoize "Manual access through - TVar"
--            , basicMeasure IORf.memoize "Manual access through - IORef"
            , basicMeasure Lock.memoize "Custom locking definition"
            , basicMeasure Semp.memoize "Manual access through - Semaphore"
            ]
