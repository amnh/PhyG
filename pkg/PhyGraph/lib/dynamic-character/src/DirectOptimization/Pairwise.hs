-----------------------------------------------------------------------------
-- |
-- Module      :  DirectOptimization.Pairwise
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-- Pairwise direct optimization alignment functions using a variety of techniques.
--
-----------------------------------------------------------------------------

module DirectOptimization.Pairwise
  (
    -- * Slim characters
    SlimDynamicCharacter
  , SlimState
  , slimPairwiseDO
    -- * Wide characters
  , WideDynamicCharacter
  , WideState
  , widePairwiseDO
    -- * Huge characters
  , HugeDynamicCharacter
  , HugeState
  , hugePairwiseDO
  ) where

import DirectOptimization.Pairwise.Huge
import DirectOptimization.Pairwise.Slim
import DirectOptimization.Pairwise.Wide

