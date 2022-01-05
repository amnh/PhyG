-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Transition.Compact.L1Norm
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE Strict #-}

module Measure.Transition.Compact.L1Norm
  ( -- * Specialized SCMs
    scmλ
  , scmρ2
  , scmρ3
  , scmρ4
  , scmρ5
  , scmρ6
  , scmρ7
  , scmρ8
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
import Foreign.C.Types        (CUInt)
import Measure.Transition.Compact.Overlap
import Measure.Transition.States.Dense   (TCMρ, fromSCMρ)
import Measure.Transition.Symbols.Dense  (SCMρ)
import Measure.Transition.Symbols.Unsafe (unsafeFromSCMλ)
import Measure.Unit.SymbolCount
import Measure.Transition.Matrix


-- |
-- The L1 Norm.
-- See:
--   https://en.wikipedia.org/wiki/Lp_space
scmλ :: SCMλ
scmλ i j = fromIntegral $ max i j - min i j


scmρ :: SymbolCount -> SCMρ
scmρ = unsafeFromSCMλ scmλ


scmρ2, scmρ3, scmρ4, scmρ5, scmρ6, scmρ7, scmρ8 :: SCMρ
scmρ2 = scmρ 2
scmρ3 = scmρ 3
scmρ4 = scmρ 4
scmρ5 = scmρ 5
scmρ6 = scmρ 6
scmρ7 = scmρ 7
scmρ8 = scmρ 8


-- |
-- Definition of the L1 norm metric.
{-# SCC        tcm2Dλ #-}
{-# INLINEABLE tcm2Dλ #-}
{-# SPECIALISE tcm2Dλ :: FiniteBits b => SymbolCount -> TCM2Dλ b      #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> TCM2Dλ Int    #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> TCM2Dλ CUInt  #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> TCM2Dλ Word   #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> TCM2Dλ Word8  #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> TCM2Dλ Word16 #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> TCM2Dλ Word32 #-}
{-# SPECIALISE tcm2Dλ ::                 SymbolCount -> TCM2Dλ Word64 #-}
tcm2Dλ
  :: FiniteBits b
  => SymbolCount
  -> TCM2Dλ b
tcm2Dλ elementBitWidth lhs rhs = overlap2 elementBitWidth scmλ lhs rhs
{-
tcm2Dλ lhs rhs
  | popCount intersect > 0 = (intersect, 0)
  | otherwise              = overlap2 scmλ lhs rhs
  where
    intersect = lhs .&. rhs
-}


-- |
-- Definition of the L1 norm metric in three dimensions.
{-# SCC        tcm3Dλ #-}
{-# INLINEABLE tcm3Dλ #-}
{-# SPECIALISE tcm3Dλ :: FiniteBits b => SymbolCount -> TCM3Dλ b      #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> TCM3Dλ Int    #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> TCM3Dλ CUInt  #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> TCM3Dλ Word   #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> TCM3Dλ Word8  #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> TCM3Dλ Word16 #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> TCM3Dλ Word32 #-}
{-# SPECIALISE tcm3Dλ ::                 SymbolCount -> TCM3Dλ Word64 #-}
tcm3Dλ
  :: FiniteBits b
  => SymbolCount
  -> TCM3Dλ b
tcm3Dλ elementBitWidth x y z = overlap3 elementBitWidth scmλ x y z


tcmρ2, tcmρ3, tcmρ4, tcmρ5, tcmρ6, tcmρ7, tcmρ8 :: TCMρ
tcmρ2 = fromSCMρ 0 scmρ2
tcmρ3 = fromSCMρ 0 scmρ3
tcmρ4 = fromSCMρ 0 scmρ4
tcmρ5 = fromSCMρ 0 scmρ5
tcmρ6 = fromSCMρ 0 scmρ6
tcmρ7 = fromSCMρ 0 scmρ7
tcmρ8 = fromSCMρ 0 scmρ8
