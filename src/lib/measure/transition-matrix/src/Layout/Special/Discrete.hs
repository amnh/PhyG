-----------------------------------------------------------------------------
-- |
-- Module      :  Layout.Special.Discrete
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

module Layout.Special.Discrete
  ( -- * Specialized SDMs
    sdmλ
  , sdmρ2
  , sdmρ3
  , sdmρ4
  , sdmρ5
  , sdmρ6
  , sdmρ7
  , sdmρ8
    -- * Specialized TCMs
  , tcm2Dλ
  , tcm3Dλ
  , tcmρ2
  , tcmρ3
  , tcmρ4
  , tcmρ5
  , tcmρ6
  , tcmρ7
  , tcmρ8
  ) where

import Data.Bits
import Data.Word
import Foreign.C.Types                   (CUInt)
import Measure.Transition
import Layout.Compact.States         (StateTransitionsCompact)
import Layout.Compact.Symbols        (SymbolDistanceMatrixSquare)
import Layout.Compact.Symbols.Unsafe (unsafeCompactStateFromSDMS, unsafeFromSDMλSquare)
import Measure.Unit.SymbolCount


sdmλ :: SDMλ
sdmλ i j
  | i == j    = 0
  | otherwise = 1


sdmρ :: Word -> SymbolDistanceMatrixSquare
sdmρ = unsafeFromSDMλSquare sdmλ . SymbolCount


sdmρ2, sdmρ3, sdmρ4, sdmρ5, sdmρ6, sdmρ7, sdmρ8 :: SymbolDistanceMatrixSquare
sdmρ2 = sdmρ 2
sdmρ3 = sdmρ 3
sdmρ4 = sdmρ 4
sdmρ5 = sdmρ 5
sdmρ6 = sdmρ 6
sdmρ7 = sdmρ 7
sdmρ8 = sdmρ 8


-- |
-- Definition of the discrete metric.
{-# SCC        tcm2Dλ #-}
{-# INLINE     tcm2Dλ #-}
{-# SPECIALISE tcm2Dλ :: Bits b => TCM2Dλ b      #-}
{-# SPECIALISE tcm2Dλ ::           TCM2Dλ Int    #-}
{-# SPECIALISE tcm2Dλ ::           TCM2Dλ CUInt  #-}
{-# SPECIALISE tcm2Dλ ::           TCM2Dλ Word   #-}
{-# SPECIALISE tcm2Dλ ::           TCM2Dλ Word8  #-}
{-# SPECIALISE tcm2Dλ ::           TCM2Dλ Word16 #-}
{-# SPECIALISE tcm2Dλ ::           TCM2Dλ Word32 #-}
{-# SPECIALISE tcm2Dλ ::           TCM2Dλ Word64 #-}
tcm2Dλ
  :: Bits b
  => TCM2Dλ b
tcm2Dλ lhs rhs
  | popCount intersect > 0 = (0, intersect)
  | otherwise              = (1,   unioned)
  where
    intersect = lhs .&. rhs
    unioned   = lhs .|. rhs


-- |
-- if           x    ⋂    y    ⋂    z    ≠ Ø ⮕  (    x    ⋂    y    ⋂    z    , 0)
--
-- else if   (x ⋂ y) ⋃ (x ⋂ z) ⋃ (y ⋂ z) ≠ Ø ⮕  ( (x ⋂ y) ⋃ (x ⋂ z) ⋃ (y ⋂ z) , 1)
--
-- otherwise                                 ⮕  (    x    ⋃    y    ⋃    z    , 2)
--
--
{-# SCC        tcm3Dλ #-}
{-# INLINE     tcm3Dλ #-}
{-# SPECIALISE tcm3Dλ :: Bits b => TCM3Dλ b      #-}
{-# SPECIALISE tcm3Dλ ::           TCM3Dλ Int    #-}
{-# SPECIALISE tcm3Dλ ::           TCM3Dλ CUInt  #-}
{-# SPECIALISE tcm3Dλ ::           TCM3Dλ Word   #-}
{-# SPECIALISE tcm3Dλ ::           TCM3Dλ Word8  #-}
{-# SPECIALISE tcm3Dλ ::           TCM3Dλ Word16 #-}
{-# SPECIALISE tcm3Dλ ::           TCM3Dλ Word32 #-}
{-# SPECIALISE tcm3Dλ ::           TCM3Dλ Word64 #-}
tcm3Dλ
  :: Bits b
  => TCM3Dλ b
tcm3Dλ x y z
  | popCount fullIntersection > 0 = (0, fullIntersection)
  | popCount joinIntersection > 0 = (1, joinIntersection)
  | otherwise                     = (2, fullUnion)
  where
    fullIntersection =  x        .&.  y        .&.  z
    joinIntersection = (x .&. y) .|. (y .&. z) .|. (z .&. x)
    fullUnion        =  x        .|.  y        .|.  z


tcmρ2, tcmρ3, tcmρ4, tcmρ5, tcmρ6, tcmρ7, tcmρ8 :: StateTransitionsCompact
tcmρ2 = unsafeCompactStateFromSDMS 0 sdmρ2
tcmρ3 = unsafeCompactStateFromSDMS 0 sdmρ3
tcmρ4 = unsafeCompactStateFromSDMS 0 sdmρ4
tcmρ5 = unsafeCompactStateFromSDMS 0 sdmρ5
tcmρ6 = unsafeCompactStateFromSDMS 0 sdmρ6
tcmρ7 = unsafeCompactStateFromSDMS 0 sdmρ7
tcmρ8 = unsafeCompactStateFromSDMS 0 sdmρ8


