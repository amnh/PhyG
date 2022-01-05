-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Distance.Types
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE RankNTypes #-}

module Measure.Distance.Types
  ( Distance
  ) where


type Distance c e = Integral c => e -> e -> c
