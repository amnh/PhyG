-----------------------------------------------------------------------------
-- |
-- Module      :  Bio.MetricRepresentation
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE BangPatterns       #-}
{-# LANGUAGE DeriveAnyClass     #-}
{-# LANGUAGE DeriveGeneric      #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE StrictData         #-}
{-# LANGUAGE UnboxedSums        #-}

module Data.MetricRepresentation
  ( -- * Smart Constructors
    MetricRepresentation()
  , discreteMetric
  , linearNorm
  , metricRepresentation
    -- * Accessors
  , minInDelCost
  , maxInDelCost
  , retreiveSCM
  , retreivePairwiseTCM
  , retreiveThreewayTCM
    -- * Metric Definitions
  , discreteMetricPairwiseLogic
  ) where

import Control.DeepSeq
import Data.Binary
import Data.Bits
import Data.Hashable
import Data.Hashable.Memoize
import Data.TCM              (TCM, (!), size)
import Data.TCM.Overlap
import GHC.Generics


-- |
-- Represents the metric for some discrete characters or dynamic characters.
-- The representation notes if the discrete metric or the L1 norm are the
-- specified metric for the character. If either of these metrics are specified,
-- specialized functions which are more efficient will be returned when
-- retrieving the pairwise of threeway transition cost matrix.
--
-- It is important to use this type in the metadata decorations rather than store
-- a function because a function cannot be contained in a compact region.
--
-- Use the elimination functions 'retreiveSCM', 'retreivePairwiseTCM', and 'retreiveThreewayTCM'
-- to the retrieve the desired functions.
data  MetricRepresentation a
    = ExplicitLayout {-# UNPACK #-} !TCM {-# UNPACK #-} !Word {-# UNPACK #-} !Word !(a -> a -> (a, Word)) !(a -> a -> a -> (a, Word))
    | DiscreteMetric
    | LinearNorm
    deriving stock    (Generic)
    deriving anyclass (NFData)


instance Eq (MetricRepresentation a) where

    (==)  DiscreteMetric             DiscreteMetric            = True
    (==)  LinearNorm                 LinearNorm                = True
    (==) (ExplicitLayout x _ _ _ _) (ExplicitLayout y _ _ _ _) = x == y
    (==) _ _                                                   = False

instance Show (MetricRepresentation a) where

    show  DiscreteMetric            = "Discrete-Metric"
    show  LinearNorm                = "1st-Linear-Norm"
    show (ExplicitLayout _ _ _ _ _) = "General-Metric"


-- |
-- Nullary constructor for the <https://en.wikipedia.org/wiki/Discrete_space discrete metric>.
discreteMetric :: MetricRepresentation a
discreteMetric = DiscreteMetric


-- |
-- Nullary constructor for the <https://en.wikipedia.org/wiki/Lp_space 1st linear norm>.
linearNorm :: MetricRepresentation a
linearNorm = LinearNorm


-- |
-- General constructor for an arbitrary metric.
--
-- Performs memoization so repeated value queries are not recomputed.
metricRepresentation
  :: ( FiniteBits a
     , Hashable a
     , NFData a
     )
  => TCM
  -> MetricRepresentation a
metricRepresentation tcm =
    let scm = makeSCM tcm
    in  ExplicitLayout tcm minInDel maxInDel
          (memoize2 (overlap2 scm))
          (memoize3 (overlap3 scm))
  where
    -- /O(2*(a - 1))/
    --
    -- This was taken from Ukkonen's original 1985 paper wherein the coefficient
    -- delta @(Δ)@ was defined by the minimum transition cost from any symbol in
    -- the alphabet @(Σ)@ to the gap symbol @'-'@.
    --
    -- If there is any transition to a gap from a non-gap for which the cost is
    -- zero, then this coefficient will be zero. This leaves us with no way to
    -- determine if optimality is preserved, and the Ukkonen algorithm will hang.
    -- Consequently, we do not perform Ukkonen's algorithm if the coefficient is
    -- zero.
    minInDel       = toEnum . fromEnum . minimum $ inDelCost min <$> nonGapElements
    maxInDel       = toEnum . fromEnum . maximum $ inDelCost max <$> nonGapElements
    alphabetSize   = size tcm
    gap            = alphabetSize - 1
    nonGapElements = [ 0 .. alphabetSize - 2 ]                                                                    
    inDelCost f i  = f (tcm ! (i  , gap))
                       (tcm ! (gap,   i))


-- |
-- /O(1)/
--
-- Extract the /minimum/ cost between a non-gap symbol and the gap symbol.
minInDelCost :: MetricRepresentation a -> Word
minInDelCost  DiscreteMetric = 1
minInDelCost  LinearNorm     = 1
minInDelCost (ExplicitLayout _ x _ _ _) = x


-- |
-- /O(1)/
--
-- Extract the /maximum/ cost between a non-gap symbol and the gap symbol.
maxInDelCost :: MetricRepresentation a -> Word
maxInDelCost  DiscreteMetric = 1
maxInDelCost  LinearNorm     = 1
maxInDelCost (ExplicitLayout _ _ x _ _) = x


-- |
-- Extract the "symbol change matrix" from a 'MetricRepresentation'.
retreiveSCM :: MetricRepresentation a -> Word -> Word -> Word
retreiveSCM (ExplicitLayout tcm _ _ _ _) = makeSCM tcm
retreiveSCM DiscreteMetric               = \i j -> if i == j then 0 else 1
retreiveSCM LinearNorm                   = l1normMetric


-- |
-- Extract the "transition cost matrix" from a 'MetricRepresentation',
-- using the elimination function.
retreivePairwiseTCM
  :: FiniteBits a
  => MetricRepresentation a
  -> a
  -> a
  -> (a, Word)
retreivePairwiseTCM (ExplicitLayout _ _ _ f _) = f
retreivePairwiseTCM DiscreteMetric             =  discreteMetricPairwiseLogic
retreivePairwiseTCM LinearNorm                 = firstLinearNormPairwiseLogic


-- |
-- Extract the threeway "transition cost matrix" from a 'MetricRepresentation',
-- using the elimination function.
retreiveThreewayTCM
  :: FiniteBits a
  => MetricRepresentation a
  -> a
  -> a
  -> a
  -> (a, Word)
retreiveThreewayTCM (ExplicitLayout _ _ _ _ f) = f
retreiveThreewayTCM DiscreteMetric             =  discreteMetricThreewayLogic
retreiveThreewayTCM LinearNorm                 = firstLinearNormThreewayLogic


-- |
-- Definition of the discrete metric.
{-# SCC        discreteMetricPairwiseLogic #-}
{-# INLINE     discreteMetricPairwiseLogic #-}
{-# SPECIALISE discreteMetricPairwiseLogic :: Bits b => b      -> b      -> (b     , Word) #-}
{-# SPECIALISE discreteMetricPairwiseLogic ::           Int    -> Int    -> (Int   , Word) #-}
{-# SPECIALISE discreteMetricPairwiseLogic ::           Word   -> Word   -> (Word  , Word) #-}
{-# SPECIALISE discreteMetricPairwiseLogic ::           Word8  -> Word8  -> (Word8 , Word) #-}
{-# SPECIALISE discreteMetricPairwiseLogic ::           Word16 -> Word16 -> (Word16, Word) #-}
{-# SPECIALISE discreteMetricPairwiseLogic ::           Word32 -> Word32 -> (Word32, Word) #-}
{-# SPECIALISE discreteMetricPairwiseLogic ::           Word64 -> Word64 -> (Word64, Word) #-}
discreteMetricPairwiseLogic
  :: ( Bits b
     , Num c
     )
  => b
  -> b
  -> (b, c)
discreteMetricPairwiseLogic !lhs !rhs
  | popCount intersect > 0 = (intersect, 0)
  | otherwise              = (  unioned, 1)
  where
    !intersect = lhs .&. rhs
    !unioned   = lhs .|. rhs


-- |
-- if           x    ⋂    y    ⋂    z    ≠ Ø ⮕  (    x    ⋂    y    ⋂    z    , 0)
--
-- else if   (x ⋂ y) ⋃ (x ⋂ z) ⋃ (y ⋂ z) ≠ Ø ⮕  ( (x ⋂ y) ⋃ (x ⋂ z) ⋃ (y ⋂ z) , 1)
--
-- otherwise                                 ⮕  (    x    ⋃    y    ⋃    z    , 2)
--
--
{-# SCC        discreteMetricThreewayLogic #-}
{-# INLINE     discreteMetricThreewayLogic #-}
{-# SPECIALISE discreteMetricThreewayLogic :: Bits b => b      -> b      -> b      -> (b     , Word) #-}
{-# SPECIALISE discreteMetricThreewayLogic ::           Int    -> Int    -> Int    -> (Int   , Word) #-}
{-# SPECIALISE discreteMetricThreewayLogic ::           Word   -> Word   -> Word   -> (Word  , Word) #-}
{-# SPECIALISE discreteMetricThreewayLogic ::           Word8  -> Word8  -> Word8  -> (Word8 , Word) #-}
{-# SPECIALISE discreteMetricThreewayLogic ::           Word16 -> Word16 -> Word16 -> (Word16, Word) #-}
{-# SPECIALISE discreteMetricThreewayLogic ::           Word32 -> Word32 -> Word32 -> (Word32, Word) #-}
{-# SPECIALISE discreteMetricThreewayLogic ::           Word64 -> Word64 -> Word64 -> (Word64, Word) #-}
discreteMetricThreewayLogic
  :: ( Bits b
     , Num c
     )
  => b
  -> b
  -> b
  -> (b, c)
discreteMetricThreewayLogic !x !y !z
  | popCount fullIntersection > 0 = (fullIntersection, 0)
  | popCount joinIntersection > 0 = (joinIntersection, 1)
  | otherwise                     = (fullUnion,        2)
  where
    !fullIntersection =  x        .&.  y        .&.  z
    !joinIntersection = (x .&. y) .|. (y .&. z) .|. (z .&. x)
    !fullUnion        =  x        .|.  y        .|.  z


-- |
-- Definition of the L1 norm metric.
{-# SCC        firstLinearNormPairwiseLogic #-}
{-# INLINEABLE firstLinearNormPairwiseLogic #-}
{-# SPECIALISE firstLinearNormPairwiseLogic :: FiniteBits b => b      -> b      -> (b     , Word) #-}
{-# SPECIALISE firstLinearNormPairwiseLogic ::                 Int    -> Int    -> (Int   , Word) #-}
{-# SPECIALISE firstLinearNormPairwiseLogic ::                 Word   -> Word   -> (Word  , Word) #-}
{-# SPECIALISE firstLinearNormPairwiseLogic ::                 Word8  -> Word8  -> (Word8 , Word) #-}
{-# SPECIALISE firstLinearNormPairwiseLogic ::                 Word16 -> Word16 -> (Word16, Word) #-}
{-# SPECIALISE firstLinearNormPairwiseLogic ::                 Word32 -> Word32 -> (Word32, Word) #-}
{-# SPECIALISE firstLinearNormPairwiseLogic ::                 Word64 -> Word64 -> (Word64, Word) #-}
firstLinearNormPairwiseLogic
  :: FiniteBits b
  => b
  -> b
  -> (b, Word)
firstLinearNormPairwiseLogic !lhs !rhs = overlap2 l1normMetric lhs rhs
{-
firstLinearNormPairwiseLogic !lhs !rhs
  | popCount intersect > 0 = (intersect, 0)
  | otherwise              = overlap2 l1normMetric lhs rhs
  where
    !intersect = lhs .&. rhs
-}


-- |
-- Definition of the L1 norm metric in three dimensions.
{-# SCC        firstLinearNormThreewayLogic #-}
{-# INLINEABLE firstLinearNormThreewayLogic #-}
{-# SPECIALISE firstLinearNormThreewayLogic :: FiniteBits b => b      -> b      -> b      -> (b     , Word) #-}
{-# SPECIALISE firstLinearNormThreewayLogic ::                 Int    -> Int    -> Int    -> (Int   , Word) #-}
{-# SPECIALISE firstLinearNormThreewayLogic ::                 Word   -> Word   -> Word   -> (Word  , Word) #-}
{-# SPECIALISE firstLinearNormThreewayLogic ::                 Word8  -> Word8  -> Word8  -> (Word8 , Word) #-}
{-# SPECIALISE firstLinearNormThreewayLogic ::                 Word16 -> Word16 -> Word16 -> (Word16, Word) #-}
{-# SPECIALISE firstLinearNormThreewayLogic ::                 Word32 -> Word32 -> Word32 -> (Word32, Word) #-}
{-# SPECIALISE firstLinearNormThreewayLogic ::                 Word64 -> Word64 -> Word64 -> (Word64, Word) #-}
firstLinearNormThreewayLogic
  :: FiniteBits b
  => b
  -> b
  -> b
  -> (b, Word)
firstLinearNormThreewayLogic !x !y !z
  | popCount intersect > 0 = (intersect, 0)
  | otherwise              = overlap3 l1normMetric x y z
  where
    !intersect = x .&. y .&. z


-- |
-- The L1 Norm.
-- See:
--   https://en.wikipedia.org/wiki/Lp_space
l1normMetric :: Word -> Word -> Word
l1normMetric i j = max i j - min i j


-- |
-- Create a symbol change matrix from a 'TCM'.
makeSCM :: TCM -> Word -> Word -> Word
makeSCM tcm i j = toEnum . fromEnum $ tcm ! (i, j)
