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

{-# LANGUAGE RankNTypes            #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Measure.Distance
  ( Distance
  , MeasurableDistance(..)
  ) where

import Measure.Distance.Types


class MeasurableDistance a c e where

    measureDistance :: a -> Distance c e
