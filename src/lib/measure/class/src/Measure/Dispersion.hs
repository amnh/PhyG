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
  , dispersion
  , dispersionPairwise
  , dispersionThreeway
  ) where

import Measure.Dispersion.General
import Measure.Dispersion.Types


class MeasurableDispersion a c e where

    measureDispersion :: a -> Dispersion c e


class MeasurableDispersionPairwise a c e where

    measureDispersionPairwise :: a -> DispersionPairwise c e


class MeasurableDispersionThreeway a c e where

    measureDispersionThreeway :: a -> DispersionThreeway c e
