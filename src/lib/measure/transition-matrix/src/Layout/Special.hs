-----------------------------------------------------------------------------
-- |
-- Module      :  Layout.Special
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
{-# LANGUAGE LambdaCase         #-}
{-# LANGUAGE MultiParamTypeClasses #-} 
{-# LANGUAGE Strict             #-}
{-# LANGUAGE UnboxedSums           #-}

module Layout.Special
  ( -- * Specializable Metrics
    SpecializableMetric()
  ) where

import Layout.Special.Type
