-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Compact.DiscreteCrossGap
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE Strict #-}

module Measure.Compact.DiscreteCrossGap
  ( -- * General functions
    scmλ
  , tcmρ
  , tcm2Dλ
  , tcm3Dλ 
    -- * Specialized SCMs
  , scm12λ
  , scm12ρ2
  , scm12ρ3
  , scm12ρ4
  , scm12ρ5
  , scm12ρ6
  , scm12ρ7
  , scm12ρ8
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
import Foreign.C.Types        (CUInt)
import Measure.Compact.Overlap
import Measure.States.Dense   (TCMρ, fromSCMρ, fromSCMλ)
import Measure.Symbols.Dense  (SCMρ)
import Measure.Symbols.Unsafe (unsafeFromSCMλ)
import Measure.Unit.SymbolCount
import Measure.Unit.Distance
import Measure.Matrix


-- |
-- The Discrete metric crossed with the gap element where the gap element where
-- the weighting ratio is 2:1
scmλ :: Distance -> Distance -> SCMλ
scmλ g s i j
  | i == j    = 0
  | i == 0    = g
  | j == 0    = g
  | otherwise = s


-- |
-- The Discrete metric crossed with the gap element where the gap element where
-- the weighting ratio is 2:1
scm12λ :: SCMλ
scm12λ i j
  | i == j    = 0
  | i == 0    = 2
  | j == 0    = 2
  | otherwise = 1



tcmρ
  :: SymbolCount
  -> Distance
  -> Distance
  -> TCMρ
tcmρ n g s = fromSCMλ n 0 $ scmλ g s


scm12ρ :: SymbolCount -> SCMρ
scm12ρ = unsafeFromSCMλ scm12λ


scm12ρ2, scm12ρ3, scm12ρ4, scm12ρ5, scm12ρ6, scm12ρ7, scm12ρ8 :: SCMρ
scm12ρ2 = scm12ρ 2
scm12ρ3 = scm12ρ 3
scm12ρ4 = scm12ρ 4
scm12ρ5 = scm12ρ 5
scm12ρ6 = scm12ρ 6
scm12ρ7 = scm12ρ 7
scm12ρ8 = scm12ρ 8


-- |
-- Definition of the L1 norm metric.
{-# SCC        tcm2Dλ #-}
{-# INLINEABLE tcm2Dλ #-}
{-# SPECIALISE tcm2Dλ :: FiniteBits b => SymbolCount -> Distance -> Distance -> TCM2Dλ b      #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> Distance -> Distance -> TCM2Dλ Int    #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> Distance -> Distance -> TCM2Dλ CUInt  #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> Distance -> Distance -> TCM2Dλ Word   #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> Distance -> Distance -> TCM2Dλ Word8  #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> Distance -> Distance -> TCM2Dλ Word16 #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> Distance -> Distance -> TCM2Dλ Word32 #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> Distance -> Distance -> TCM2Dλ Word64 #-}
tcm2Dλ
  :: FiniteBits b
  => SymbolCount
  -> Distance
  -> Distance
  -> TCM2Dλ b
tcm2Dλ elementBitWidth gapCost subCost lhs rhs = overlap2 elementBitWidth (scmλ gapCost subCost) lhs rhs


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
tcm2D12λ elementBitWidth lhs rhs = overlap2 elementBitWidth scm12λ lhs rhs


-- |
-- Definition of the L1 norm metric.
{-# SCC        tcm3Dλ #-}
{-# INLINEABLE tcm3Dλ #-}
{-# SPECIALISE tcm3Dλ :: FiniteBits b => SymbolCount -> Distance -> Distance -> TCM3Dλ b      #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> Distance -> Distance -> TCM3Dλ Int    #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> Distance -> Distance -> TCM3Dλ CUInt  #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> Distance -> Distance -> TCM3Dλ Word   #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> Distance -> Distance -> TCM3Dλ Word8  #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> Distance -> Distance -> TCM3Dλ Word16 #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> Distance -> Distance -> TCM3Dλ Word32 #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> Distance -> Distance -> TCM3Dλ Word64 #-}
tcm3Dλ
  :: FiniteBits b
  => SymbolCount
  -> Distance
  -> Distance
  -> TCM3Dλ b
tcm3Dλ elementBitWidth gapCost subCost x y z = overlap3 elementBitWidth (scmλ gapCost subCost) x y z


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
tcm3D12λ elementBitWidth x y z = overlap3 elementBitWidth scm12λ x y z


tcm12ρ2, tcm12ρ3, tcm12ρ4, tcm12ρ5, tcm12ρ6, tcm12ρ7, tcm12ρ8 :: TCMρ
tcm12ρ2 = fromSCMρ 0 scm12ρ2
tcm12ρ3 = fromSCMρ 0 scm12ρ3
tcm12ρ4 = fromSCMρ 0 scm12ρ4
tcm12ρ5 = fromSCMρ 0 scm12ρ5
tcm12ρ6 = fromSCMρ 0 scm12ρ6
tcm12ρ7 = fromSCMρ 0 scm12ρ7
tcm12ρ8 = fromSCMρ 0 scm12ρ8
