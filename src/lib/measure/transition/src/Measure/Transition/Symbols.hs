-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Transition.Symbols
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

module Measure.Transition.Symbols
  ( HasSymbolDistances(..)
  ) where

import Data.Coerce
import Layout.Compact.States.Indexing
import Measure.Transition.Types
import Measure.Unit.SymbolChangeCost
import Measure.Unit.SymbolIndex


-- |
-- Any structural representation which can produce a function measuring the
-- 'Measure.Distance.Distance' between two symbols, i.e. produce a "Symbol
-- Change Matrix."
class HasSymbolDistances a where

    symbolDistances :: a -> SymbolDistanceλ


instance HasSymbolDistances StateTransitionsCompact where

    symbolDistances tcm i =
        let f = coerce :: SymbolIndex -> Word
            g = coerce :: Word -> SymbolChangeCost
        in  g . costSymbol tcm (f i) . f
