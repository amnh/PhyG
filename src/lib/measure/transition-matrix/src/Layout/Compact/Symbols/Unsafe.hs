-----------------------------------------------------------------------------
-- |
-- Module      :  Layout.Compact.Symbols.Unsafe
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RankNTypes       #-}
{-# LANGUAGE Strict           #-}

module Layout.Compact.Symbols.Unsafe
  ( -- * Unsafe
    unsafeFromSDMλSquare
  , unsafeFromSDMλTriangular
  , unsafeFromVectorSquare
  , unsafeFromVectorTriangular
  , unsafeCompactStateFromSDMS
  ) where

import Data.Coerce
import Data.Vector.Storable              (Vector, force, generate)
import Data.Word (Word16)
import Layout.Compact.States             (StateTransitionsCompact, initialize)
import Layout.Compact.Symbols.Internal   (SymbolDistanceMatrix(..))
import Layout.Compact.Symbols.Square     (SymbolDistanceMatrixSquare(..), rowMajorVector)
import Layout.Compact.Symbols.Triangular (SymbolDistanceMatrixTriangular(..))
import Measure.Transition
import Measure.Unit.SymbolChangeCost
import Measure.Unit.SymbolCount
import Measure.Unit.SymbolIndex


-- |
-- /ϴ(a⁵)/ where /a ≤ 8/ is the size of the character alphabet.
--
-- Generate the 2D and 3D compact state transiton cost matricies from the supplied
-- symbol distance matrix with linear dimensions of the alphabet symbol count.
unsafeCompactStateFromSDMS 
  :: SymbolChangeCost           -- ^ The gap open cost. A zero value indicates non-affine alignment context
  -> SymbolDistanceMatrixSquare -- ^ The dense, pre-computed matrix of costs to shift between symbols.
  -> StateTransitionsCompact
unsafeCompactStateFromSDMS penalty sdms =
    let dimension  = coerce $ symbolCount sdms
        gapSeqCost = coerce penalty
        costVector = coerce $ rowMajorVector sdms
    in  initialize dimension gapSeqCost costVector


unsafeFromSDMλSquare :: SDMλ -> SymbolCount -> SymbolDistanceMatrixSquare
unsafeFromSDMλSquare sdmλ sc@(SymbolCount n) =
    let dim = fromEnum $ n * n
        g i = let (q,r) = toEnum i `quotRem` n
              in  fromIntegral . sdmλ (SymbolIndex q) $ SymbolIndex r
    in  coerce . SymbolDistanceMatrix sc . force $ generate dim g


unsafeFromSDMλTriangular :: SDMλ -> SymbolCount -> SymbolDistanceMatrixTriangular
unsafeFromSDMλTriangular sdmλ sc@(SymbolCount n) =
    let dim = fromEnum $ n * n
        g i = let (q,r) = toEnum i `quotRem` n
              in  fromIntegral . sdmλ (SymbolIndex q) $ SymbolIndex r
    in  coerce . SymbolDistanceMatrix sc . force $ generate dim g


unsafeFromVectorSquare :: SymbolCount -> Vector Word16 -> SymbolDistanceMatrixSquare
unsafeFromVectorSquare = (coerce .) . (. force) . SymbolDistanceMatrix

unsafeFromVectorTriangular :: SymbolCount -> Vector Word16 -> SymbolDistanceMatrix
unsafeFromVectorTriangular = (coerce .) . (. force) . SymbolDistanceMatrix
