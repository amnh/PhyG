------------------------------------------------------------------------------
-- |
-- Module      :  Data.MutualExclusionSet.Bench
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE BangPatterns #-}

module Data.MutualExclusionSet.Bench (benchmarks) where

import Control.DeepSeq
import Criterion.Main
import Data.Bits
import Data.MutualExclusionSet


-- |
-- Com plet set of benchmarks for 'MutualExclusionSet'.
benchmarks :: Benchmark
benchmarks = bgroup "MutualExclusionSet"
    [ singletonBench
    , invertBench
    , isCoherentBench
    , isExcludedBench
    , isIncludedBench
    , excludedLookupBench
    , includedLookupBench
    , excludedSetBench
    , includedSetBench
    , mutualExclusivePairsBench
    , mergeBench
    , isPermissibleBench
    ]


-- |
-- Benchmarks the construction of a singleton 'MutualExclusionSet'.
singletonBench :: Benchmark
singletonBench = bench "MutualExclusionSet singleton is constant-time construction" . nf (singleton 42) $ (1 :: Int)


-- |
-- Benchmarks the 'invert' transformation of a 'MutualExclusionSet'.
invertBench :: Benchmark
invertBench = bench "MutualExclusionSet invert is constant-time" . whnf invert $ force (ofSize 50)


-- |
-- Benchmarks the 'isCoherent' check of a 'MutualExclusionSet'.
isCoherentBench :: Benchmark
isCoherentBench = bench "MutualExclusionSet isCoherent is constant-time" . whnf isCoherent $ force (ofSize 50)


-- |
-- Benchmarks the 'isPermissible' check of a 'MutualExclusionSet'.
isPermissibleBench :: Benchmark
isPermissibleBench = linearBenchmark "MutualExclusionSet isPermissible log-access" (force . ofSize) (const isPermissible)


-- |
-- Benchmarks the 'isExcluded' query of a 'MutualExclusionSet'.
isExcludedBench :: Benchmark
isExcludedBench = logBenchmark "MutualExclusionSet isExcluded log-access" ofSize f
  where
    -- We negate i to consider both the included and excluded cases
    f i xs = (i `isExcluded` xs) `seq` negate i `isExcluded` xs


-- |
-- Benchmarks the 'isIncluded' query of a 'MutualExclusionSet'.
isIncludedBench :: Benchmark
isIncludedBench = logBenchmark "MutualExclusionSet isIncluded log-access" ofSize f
  where
    -- We negate i to consider both the included and excluded cases
    f i xs = (i `isIncluded` xs) `seq` negate i `isIncluded` xs


excludedLookupBench :: Benchmark
excludedLookupBench = logBenchmark "MutualExclusionSet excludedLookup log-access" ofSize f
  where
    -- We negate i to consider both the included and excluded cases
    f i xs = (i `excludedLookup` xs) `seq` negate i `excludedLookup` xs


includedLookupBench :: Benchmark
includedLookupBench = logBenchmark "MutualExclusionSet includedLookup log-access" ofSize f
  where
    -- We negate i to consider both the included and excluded cases
    f i xs = (i `includedLookup` xs) `seq` negate i `includedLookup` xs


excludedSetBench :: Benchmark
excludedSetBench = linearBenchmark "MutualExclusionSet excludedSet linear access" (force . ofSize) (const excludedSet)


includedSetBench :: Benchmark
includedSetBench = linearBenchmark "MutualExclusionSet includedSet linear access" (force . ofSize) (const includedSet)


mutualExclusivePairsBench :: Benchmark
mutualExclusivePairsBench = linearBenchmark "MutualExclusionSet mutuallyExclusivePairs linear access" (force . ofSize) (const mutuallyExclusivePairs)


mergeBench :: Benchmark
mergeBench = linearBenchmark2 "merge (<>) is linear" (force . ofSize) (force . ofSizeEven) (const (<>))


linearBenchmark :: (NFData a, NFData b) => String -> (Int -> a) -> (Int -> a -> b) -> Benchmark
linearBenchmark  label f g = bgroup label $ generateBenchmark <$> [0 .. 9]
  where
    generateBenchmark expVal = bench (show domainSize) $ nf app target
      where
        !target    = force $ f domainSize
        !app       = g expVal
        domainSize = 10 * (expVal + 1)


linearBenchmark2 :: (NFData a, NFData b, NFData c) => String -> (Int -> a) -> (Int -> b) -> (Int -> a -> b -> c) -> Benchmark
linearBenchmark2  label f g h = bgroup label $ generateBenchmark <$> [0 .. 9]
  where
    generateBenchmark expVal = bench (show domainSize) $ nf app rhs
      where
        !lhs       = force $ f domainSize
        !rhs       = force $ g domainSize
        !app       = h expVal lhs
        domainSize = 10 * (expVal + 1)


logBenchmark :: String -> (Int -> a) -> (Int -> a -> b) -> Benchmark
logBenchmark label f g = bgroup label $ generateBenchmark <$> [0 .. 9]
  where
    generateBenchmark expVal = bench (show domainSize) $ whnf app target
      where
        !app       = g indexpValrod
        !target    = f domainSize
        indexpValrod  = product [1..expVal] `mod` domainSize
        domainSize = 2 `shiftL` expVal


ofSize :: Int -> MutualExclusionSet Int
ofSize n = unsafeFromList $ (\x -> (x, negate x)) <$> [1..n]


ofSizeEven :: Int -> MutualExclusionSet Int
ofSizeEven n = unsafeFromList $ (\x -> (x, negate x)) <$> [2,4..2*n]
