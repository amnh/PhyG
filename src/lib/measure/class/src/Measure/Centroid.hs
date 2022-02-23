-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Centroid
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE RankNTypes            #-}

module Measure.Centroid
  ( Centroid
  , CentroidPairwise
  , CentroidThreeway
  , MeasurableCentroid(..)
  , MeasurableCentroidPairwise(..)
  , MeasurableCentroidThreeway(..)
  ) where

import Data.Semigroup.Foldable (Foldable1)


-- |
-- Abstract function computing a [Centroid](https://en.wikipedia.org/wiki/Centroid).
type Centroid e = forall f. Foldable1 f => f e -> e


-- |
-- Compute the [Centroid](https://en.wikipedia.org/wiki/Centroid) of a pair.
type CentroidPairwise e = e -> e -> e


-- |
-- Compute the [Centroid](https://en.wikipedia.org/wiki/Centroid) of a triple.
type CentroidThreeway e = e -> e -> e -> e


-- |
-- Structure which can compute a [Centroid](https://en.wikipedia.org/wiki/Centroid).
class MeasurableCentroid a e where

    measureCentroid :: a -> Centroid e


-- |
-- Structure which can compute the [Centroid](https://en.wikipedia.org/wiki/Centroid) of a pair.
class MeasurableCentroidPairwise a e where

    measureCentroidPairwise :: a -> CentroidPairwise e


-- |
-- Structure which can compute the [Centroid](https://en.wikipedia.org/wiki/Centroid) of a triple.
class MeasurableCentroidThreeway a e where

    measureCentroidThreeway :: a -> CentroidThreeway e
