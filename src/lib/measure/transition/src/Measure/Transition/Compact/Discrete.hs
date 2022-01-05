-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Transition.Compact.Discrete
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE Strict #-}

module Measure.Transition.Compact.Discrete
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
import Measure.Transition.States.Dense   (TCMρ, fromSCMρ)
import Measure.Transition.Symbols.Dense  (SCMρ)
import Measure.Transition.Symbols.Unsafe (unsafeFromSCMλ)
import Measure.Unit.SymbolCount
import Measure.Transition.Matrix


scmλ :: SCMλ
scmλ i j
  | i == j    = 0
  | otherwise = 1


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


tcmρ2, tcmρ3, tcmρ4, tcmρ5, tcmρ6, tcmρ7, tcmρ8 :: TCMρ
tcmρ2 = fromSCMρ 0 scmρ2
tcmρ3 = fromSCMρ 0 scmρ3
tcmρ4 = fromSCMρ 0 scmρ4
tcmρ5 = fromSCMρ 0 scmρ5
tcmρ6 = fromSCMρ 0 scmρ6
tcmρ7 = fromSCMρ 0 scmρ7
tcmρ8 = fromSCMρ 0 scmρ8


