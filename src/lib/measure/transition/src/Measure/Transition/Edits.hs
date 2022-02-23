-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Transition.Edits
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE MultiParamTypeClasses #-}

module Measure.Transition.Edits
  ( HasEditExtrema(..)
  ) where

import Data.Coerce
import Layout.Compact.States
import Measure.Unit.SymbolChangeCost


-- |
-- Memoized/pre-computed access to the extrema costs for each type of edit.
class HasEditExtrema a where

    {-# MINIMAL maxDeletion, maxInsertion, minDeletion, minInsertion #-}

    {-# INLINEABLE maxEdit #-}
    maxEdit :: a -> SymbolChangeCost
    maxEdit a = max (maxDeletion a) $ maxInsertion a

    maxDeletion  :: a -> SymbolChangeCost

    maxInsertion :: a -> SymbolChangeCost

    {-# INLINEABLE minEdit #-}
    minEdit :: a -> SymbolChangeCost
    minEdit a = min (minDeletion a) $ minInsertion a

    minDeletion  :: a -> SymbolChangeCost

    minInsertion :: a -> SymbolChangeCost


instance HasEditExtrema StateTransitionsCompact where

    {-# INLINE maxDeletion #-}
    maxDeletion  = (coerce :: Word -> SymbolChangeCost) . maxDelCost

    {-# INLINE maxInsertion #-}
    maxInsertion = (coerce :: Word -> SymbolChangeCost) . maxInsCost

    {-# INLINE minDeletion #-}
    minDeletion  = (coerce :: Word -> SymbolChangeCost) . minDelCost

    {-# INLINE minInsertion #-}
    minInsertion = (coerce :: Word -> SymbolChangeCost) . minInsCost
