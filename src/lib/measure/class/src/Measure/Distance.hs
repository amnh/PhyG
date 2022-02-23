-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Distance
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

module Measure.Distance
  ( Distance
  , MeasurableDistance(..)
  ) where


-- |
-- Abstract computation of the [Distance](https://en.wikipedia.org/wiki/Distance) between two "points."
type Distance c e = e -> e -> c


-- |
-- Structure which computes the [Distance](https://en.wikipedia.org/wiki/Distance) between two "points."
class MeasurableDistance a c e where

    measureDistance :: a -> Distance c e

