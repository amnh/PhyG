-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Dispersion
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

module Measure.Dispersion
  ( Dispersion
  , DispersionPairwise
  , DispersionThreeway
  , MeasurableDispersion(..)
  , MeasurableDispersionPairwise(..)
  , MeasurableDispersionThreeway(..)
  ) where

import Data.Semigroup.Foldable (Foldable1)


-- |
-- Abstract function computing the [Dispersion](https://en.wikipedia.org/wiki/Statistical_dispersion)
-- of a collection of points, returning the dispersion 'Measure.Distance.Distance'
-- and 'Measure.Centroid.Centroid'.
type Dispersion c e = forall f. Foldable1 f => f e -> (c, e)


-- |
-- Abstract function computing the [Dispersion](https://en.wikipedia.org/wiki/Statistical_dispersion)
-- of a pair.
type DispersionPairwise c e = e -> e -> (c, e)


-- |
-- Abstract function computing the [Dispersion](https://en.wikipedia.org/wiki/Statistical_dispersion)
-- of a triple.
type DispersionThreeway c e = e -> e -> e -> (c, e)


-- |
-- Structure which can compute a [Dispersion](https://en.wikipedia.org/wiki/Statistical_dispersion).
class MeasurableDispersion a c e where

    measureDispersion :: a -> Dispersion c e


-- |
-- Structure which can compute the [Dispersion](https://en.wikipedia.org/wiki/Statistical_dispersion) of a pair.
class MeasurableDispersionPairwise a c e where

    measureDispersionPairwise :: a -> DispersionPairwise c e


-- |
-- Structure which can compute the [Dispersion](https://en.wikipedia.org/wiki/Statistical_dispersion) of a triple.
class MeasurableDispersionThreeway a c e where

    measureDispersionThreeway :: a -> DispersionThreeway c e
