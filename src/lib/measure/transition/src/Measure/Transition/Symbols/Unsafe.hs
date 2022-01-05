-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Transition.Symbols.Unsafe
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE Strict #-}

module Measure.Transition.Symbols.Unsafe
  ( -- * Unsafe
    unsafeFromSCMλ
  ) where

import Data.Vector.Storable     (force, generate)
import Measure.Transition.Symbols.Internal (SCMρ(..))
import Measure.Transition.Matrix
import Measure.Unit.SymbolCount


unsafeFromSCMλ :: SCMλ -> SymbolCount -> SCMρ
unsafeFromSCMλ scmλ n =
  let dim = fromIntegral $ n * n
      row = fromIntegral n
      g i = fromIntegral . uncurry scmλ $ fromIntegral i `quotRem` row
  in  SCMρ n . force $ generate dim g
