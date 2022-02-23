-----------------------------------------------------------------------------
-- |
-- Module      :  Layout.Special.DiscreteCrossGap
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

module Layout.Special.DiscreteCrossGap
  ( -- * General functions
    sdmλ
  , tcmρ
  , tcm2Dλ
  , tcm3Dλ
    -- * Specialized SDMs
  , sdm12λ
  , sdm12ρ2
  , sdm12ρ3
  , sdm12ρ4
  , sdm12ρ5
  , sdm12ρ6
  , sdm12ρ7
  , sdm12ρ8
    -- * Specialized TCMs
  , tcm2D12λ
  , tcm3D12λ
  , tcm12ρ2
  , tcm12ρ3
  , tcm12ρ4
  , tcm12ρ5
  , tcm12ρ6
  , tcm12ρ7
  , tcm12ρ8
  ) where

import Data.Bits
import Data.Word
import Foreign.C.Types                          (CUInt)
import Measure.Transition
import Layout.Memoize.Dispersion
import Layout.Compact.States         (StateTransitionsCompact)
import Layout.Compact.Symbols        (SymbolDistanceMatrixSquare)
import Layout.Compact.Symbols.Unsafe (unsafeCompactStateFromSDMS, unsafeFromSDMλSquare)
import Measure.Unit.SymbolChangeCost
import Measure.Unit.SymbolCount
import Measure.Unit.SymbolIndex


-- |
-- The Discrete metric crossed with the gap element where the gap element where
-- the weighting ratio is 2:1
sdmλ :: SymbolChangeCost -> SymbolChangeCost -> SDMλ
sdmλ g s i j
  | i == j             = 0
  | i == SymbolIndex 0 = g
  | j == SymbolIndex 0 = g
  | otherwise          = s


-- |
-- The Discrete metric crossed with the gap element where the gap element where
-- the weighting ratio is 2:1
sdm12λ :: SDMλ
sdm12λ i j
  | i == j             = 0
  | i == SymbolIndex 0 = 2
  | j == SymbolIndex 0 = 2
  | otherwise          = 1



tcmρ
  :: SymbolCount
  -> SymbolChangeCost
  -> SymbolChangeCost
  -> StateTransitionsCompact
tcmρ n g s = unsafeCompactStateFromSDMS 0 $ unsafeFromSDMλSquare (sdmλ g s) n


sdm12ρ :: Word -> SymbolDistanceMatrixSquare
sdm12ρ = unsafeFromSDMλSquare sdm12λ . SymbolCount


sdm12ρ2, sdm12ρ3, sdm12ρ4, sdm12ρ5, sdm12ρ6, sdm12ρ7, sdm12ρ8 :: SymbolDistanceMatrixSquare
sdm12ρ2 = sdm12ρ 2
sdm12ρ3 = sdm12ρ 3
sdm12ρ4 = sdm12ρ 4
sdm12ρ5 = sdm12ρ 5
sdm12ρ6 = sdm12ρ 6
sdm12ρ7 = sdm12ρ 7
sdm12ρ8 = sdm12ρ 8


-- |
-- Definition of the L1 norm metric.
{-# SCC        tcm2Dλ #-}
{-# INLINEABLE tcm2Dλ #-}
{-# SPECIALISE tcm2Dλ :: FiniteBits b => SymbolCount -> SymbolChangeCost -> SymbolChangeCost -> TCM2Dλ b      #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> SymbolChangeCost -> SymbolChangeCost -> TCM2Dλ Int    #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> SymbolChangeCost -> SymbolChangeCost -> TCM2Dλ CUInt  #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> SymbolChangeCost -> SymbolChangeCost -> TCM2Dλ Word   #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> SymbolChangeCost -> SymbolChangeCost -> TCM2Dλ Word8  #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> SymbolChangeCost -> SymbolChangeCost -> TCM2Dλ Word16 #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> SymbolChangeCost -> SymbolChangeCost -> TCM2Dλ Word32 #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> SymbolChangeCost -> SymbolChangeCost -> TCM2Dλ Word64 #-}
tcm2Dλ
  :: FiniteBits b
  => SymbolCount
  -> SymbolChangeCost
  -> SymbolChangeCost
  -> TCM2Dλ b
tcm2Dλ elementBitWidth gapCost subCost =
    bitDispersionPairwise (symbolBounds elementBitWidth) (sdmλ gapCost subCost)


-- |
-- Definition of the L1 norm metric.
{-# SCC        tcm2D12λ #-}
{-# INLINEABLE tcm2D12λ #-}
{-# SPECIALISE tcm2D12λ :: FiniteBits b => SymbolCount -> TCM2Dλ b      #-}
{-# SPECIALISE tcm2D12λ ::                 SymbolCount -> TCM2Dλ Int    #-}
{-# SPECIALISE tcm2D12λ ::                 SymbolCount -> TCM2Dλ CUInt  #-}
{-# SPECIALISE tcm2D12λ ::                 SymbolCount -> TCM2Dλ Word   #-}
{-# SPECIALISE tcm2D12λ ::                 SymbolCount -> TCM2Dλ Word8  #-}
{-# SPECIALISE tcm2D12λ ::                 SymbolCount -> TCM2Dλ Word16 #-}
{-# SPECIALISE tcm2D12λ ::                 SymbolCount -> TCM2Dλ Word32 #-}
{-# SPECIALISE tcm2D12λ ::                 SymbolCount -> TCM2Dλ Word64 #-}
tcm2D12λ
  :: FiniteBits b
  => SymbolCount
  -> TCM2Dλ b
tcm2D12λ elementBitWidth =
    bitDispersionPairwise (symbolBounds elementBitWidth) sdm12λ


-- |
-- Definition of the L1 norm metric.
{-# SCC        tcm3Dλ #-}
{-# INLINEABLE tcm3Dλ #-}
{-# SPECIALISE tcm3Dλ :: FiniteBits b => SymbolCount -> SymbolChangeCost -> SymbolChangeCost -> TCM3Dλ b      #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> SymbolChangeCost -> SymbolChangeCost -> TCM3Dλ Int    #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> SymbolChangeCost -> SymbolChangeCost -> TCM3Dλ CUInt  #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> SymbolChangeCost -> SymbolChangeCost -> TCM3Dλ Word   #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> SymbolChangeCost -> SymbolChangeCost -> TCM3Dλ Word8  #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> SymbolChangeCost -> SymbolChangeCost -> TCM3Dλ Word16 #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> SymbolChangeCost -> SymbolChangeCost -> TCM3Dλ Word32 #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> SymbolChangeCost -> SymbolChangeCost -> TCM3Dλ Word64 #-}
tcm3Dλ
  :: FiniteBits b
  => SymbolCount
  -> SymbolChangeCost
  -> SymbolChangeCost
  -> TCM3Dλ b
tcm3Dλ elementBitWidth gapCost subCost =
    bitDispersionThreeway (symbolBounds elementBitWidth) (sdmλ gapCost subCost)


-- |
-- Definition of the L1 norm metric in three dimensions.
{-# SCC        tcm3D12λ #-}
{-# INLINEABLE tcm3D12λ #-}
{-# SPECIALISE tcm3D12λ :: FiniteBits b => SymbolCount -> TCM3Dλ b      #-}
{-# SPECIALISE tcm3D12λ ::                 SymbolCount -> TCM3Dλ Int    #-}
{-# SPECIALISE tcm3D12λ ::                 SymbolCount -> TCM3Dλ CUInt  #-}
{-# SPECIALISE tcm3D12λ ::                 SymbolCount -> TCM3Dλ Word   #-}
{-# SPECIALISE tcm3D12λ ::                 SymbolCount -> TCM3Dλ Word8  #-}
{-# SPECIALISE tcm3D12λ ::                 SymbolCount -> TCM3Dλ Word16 #-}
{-# SPECIALISE tcm3D12λ ::                 SymbolCount -> TCM3Dλ Word32 #-}
{-# SPECIALISE tcm3D12λ ::                 SymbolCount -> TCM3Dλ Word64 #-}
tcm3D12λ
  :: FiniteBits b
  => SymbolCount
  -> TCM3Dλ b
tcm3D12λ elementBitWidth =
    bitDispersionThreeway (symbolBounds elementBitWidth) sdm12λ


tcm12ρ2, tcm12ρ3, tcm12ρ4, tcm12ρ5, tcm12ρ6, tcm12ρ7, tcm12ρ8 :: StateTransitionsCompact
tcm12ρ2 = unsafeCompactStateFromSDMS 0 sdm12ρ2
tcm12ρ3 = unsafeCompactStateFromSDMS 0 sdm12ρ3
tcm12ρ4 = unsafeCompactStateFromSDMS 0 sdm12ρ4
tcm12ρ5 = unsafeCompactStateFromSDMS 0 sdm12ρ5
tcm12ρ6 = unsafeCompactStateFromSDMS 0 sdm12ρ6
tcm12ρ7 = unsafeCompactStateFromSDMS 0 sdm12ρ7
tcm12ρ8 = unsafeCompactStateFromSDMS 0 sdm12ρ8
