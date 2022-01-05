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

{-# LANGUAGE RankNTypes            #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Measure.Centroid
  ( Centroid
  , CentroidPairwise
  , CentroidThreeway
  , MeasurableCentroid(..)
  , MeasurableCentroidPairwise(..)
  , MeasurableCentroidThreeway(..)
  ) where

import Data.Semigroup.Foldable (Foldable1)


type Centroid e = forall f. Foldable1 f => f e -> e


type CentroidPairwise e = e -> e -> e


type CentroidThreeway e = e -> e -> e -> e


class MeasurableCentroid a e where

    measureCentroid :: a -> Centroid e


class MeasurableCentroidPairwise a e where

    measureCentroidPairwise :: a -> CentroidPairwise e


class MeasurableCentroidThreeway a e where

    measureCentroidThreeway :: a -> CentroidThreeway e
