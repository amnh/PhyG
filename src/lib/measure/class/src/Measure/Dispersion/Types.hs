-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Dispersion.Types
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE RankNTypes            #-}

module Measure.Dispersion.Types
  ( Dispersion
  , DispersionPairwise
  , DispersionThreeway
  ) where

import Data.Semigroup.Foldable (Foldable1)


type Dispersion c e = forall f. (Foldable1 f, Integral c) => f e -> (c, e)


type DispersionPairwise c e = Integral c => e -> e -> (c, e)


type DispersionThreeway c e = Integral c => e -> e -> e -> (c, e)
