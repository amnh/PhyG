-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Metricity
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE DeriveAnyClass     #-}
{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE DeriveGeneric      #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE Strict             #-}
{-# LANGUAGE TypeFamilies       #-}

module Measure.Metricity
  ( -- * Diagnosable Metrics
    Metricity(..)
    -- * Diagnoses
  , metricityOfSCMρ
  , metricityOfSCMλ
  ) where

import Control.DeepSeq
import Data.Bits
import Data.Data
import Data.Foldable
import GHC.Generics
import Measure.Symbols.Internal
import Measure.Matrix
import Measure.Unit.SymbolCount
import Measure.Unit.SymbolIndex
import Measure.Unit.Distance


-- |
-- There is a heirachichal nature if a Symbol Change Matrix's metricity.
--
-- @
--       Non-metric
--
--           |
--
--         Metric
--
--         /    \\
--
-- Ultrametric  L1 Norm
--
--      |
--
-- Discrete metric ⨉ Gap
--
--      |
--
-- Discrete metric
--
-- @
--
-- Each SCM structure has certain properties that must hold and allows for space
-- & time optimizations.
--
--   * Non-metric: Generate a warning for the user!
--
--   * Metric:
--
--       * /σ(i,j) = 0 iff i = j/
--
--       * /σ(i,j) = σ(j,i)/
--
--       * /σ(i,k) ≤ σ(i,j) + σ(j,k)/
--
--  * Ultrametric: Allows for some runtime optimizations
--
--       * /σ(i,j) = 0 iff i = j/
--
--       * /σ(i,j) = σ(j,i)/
--
--       * /σ(i,k) ≤ max { σ(i,j),  σ(j,k) }/
--
--   * Discrete metric ⨉ Gap: Allows for /no/ space allocation and runtime optimizations
--
--       * /σ(i,j) = 0 iff i = j/
--
--       * /σ(0,j) = c iff i ≠ j/
--
--       * /σ(i,0) = c iff i ≠ j/
--
--       * /σ(i,j) = 1 iff i ≠ j, i > 0, j > 0/
--
--   * Discrete metric: Allows for /no/ space allocation and runtime optimizations
--
--       * /σ(i,j) = 0 iff i = j/
--
--       * /σ(i,j) = 1 iff i ≠ j/
--
--   * L1 Norm: Allows for /no/ space allocation and runtime optimizations
--
--       * /σ(i,j) = max(i,j) - min(i,j)/
--
data  Metricity
    = NonMetric
    | Metric
    | Ultrametric
    | DiscreteCrossGap {-# UNPACK #-} !Word {-# UNPACK #-} !Word
    | DiscreteMetric
    | L1Norm
    deriving stock    (Data, Eq, Generic, Show, Typeable)
    deriving anyclass (NFData)


metricityOfSCMρ :: SCMρ -> Metricity
metricityOfSCMρ scmρ =
    let dimension = sizeOfSCM scmρ
        indexSCMρ = index scmρ
    in  metricityOfSCMλ dimension indexSCMρ


metricityOfSCMλ :: SymbolCount -> SCMλ -> Metricity
metricityOfSCMλ dim scmλ =
  let range  = [ 0 .. fromIntegral dim - 1 ]
  in  if nonZeroDiagonal scmλ range
      then  NonMetric
      else  case pass2D scmλ range of
              Just m -> m
              -- We know that the matrix is symetric,
              -- With zeroes along the diagonal,
              -- but it is not the Distcrete metric,
              -- nor the Discrete metric plus Gap
              -- nor the 1st linear norm.
              --
              -- Hence we check for metricity ancd ultrametricity
              Nothing ->
                let points =
                      [ (scmλ i k, scmλ i j, scmλ j k)
                      | i <- range
                      , j <- range
                      , i < j
                      , k <- range
                      , j <= k
                      ]

                    triangleInequality ~(x,y,z) = x <= y + z
                    strongerInequality ~(x,y,z) = x <= max y z

                    metricCheck  result             [] = result
                    metricCheck  (False,     _)     _  = (False, False)
                    metricCheck  (_    , False) (p:ps) = metricCheck (triangleInequality p, False) ps
                    metricCheck              _  (p:ps) = metricCheck (triangleInequality p, strongerInequality p) ps

                in  case metricCheck (True, True) points of
                      (_, True) -> Ultrametric
                      (True, _) -> Metric
                      _         -> NonMetric


-- |
-- An internal helper function used in both 'isMetric' & 'isUltraMetric' exported functions.
nonZeroDiagonal :: Foldable f => SCMλ -> f SymbolIndex -> Bool
nonZeroDiagonal scmλ = any (\i -> scmλ i i /= 0)


-- |
-- Just one pass over the upper triangluar of the matrix
pass2D :: Foldable f => SCMλ -> f SymbolIndex -> Maybe Metricity
pass2D scmλ range
  | resultBits `testBit` discrete = Just DiscreteMetric
  | resultBits `testBit` subindel = Just gapSubMetric
  | resultBits `testBit` linear1N = Just L1Norm
  | resultBits `testBit` symetric = Nothing
  | otherwise = Just NonMetric
  where
    resultBits = foldr check initialBits pairs

    gapSubMetric = DiscreteCrossGap (fromIntegral gDist) $ fromIntegral sDist

    dim = length range

    discrete = 0 :: Int
    subindel = 1
    linear1N = 2
    symetric = 3

    isDiscrete (_,_,x,_) = x == 1
    isInDelSub (i,j,x,_) = x == if i == 0 || j == 0 then gDist else sDist
    isLinear1N (i,j,x,_) = x == fromIntegral (max i j - min i j)

    bitSet = [ (discrete, isDiscrete)
             , (subindel, isInDelSub)
             , (linear1N, isLinear1N)
             , (symetric, const True)
             ]

    gDist, sDist :: Distance
    gDist = scmλ 0 1
    sDist | dim > 2   = scmλ 1 2
          | otherwise = maxBound

    check p@(_, _, x, y) bits
      | x /= y    = zeroBits
      | otherwise =
          let f b (i, g)
                | b `testBit` i && not (g p) = b `clearBit` i
                | otherwise = b
          in  foldl' f bits $ init bitSet

    initialBits :: Word
    initialBits =
        let b = foldl' setBit zeroBits $ fst <$> bitSet
        in  if dim > 2
            then b
            else b `clearBit` subindel

    pairs = [ (i, j, scmλ i j, scmλ j i) | i <- toList range, j <- toList range, i < j ]

