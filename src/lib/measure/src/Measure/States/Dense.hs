-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.States.Dense
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-- Construction of a dense transition cost matrix 'TCMρ'
--
-- /O(a^5)/ where /a/ is the size of the character alphabet
-----------------------------------------------------------------------------

module Measure.States.Dense
  ( -- * Dense TCM Representation
    TCMρ()
    -- * FFI TCM Components
  , FFI2D()
  , FFI3D()
    -- * FFI Alignment Method
  , AlignmentStrategy(..)
    -- * Accessor Functions
  , getAlignmentStrategy
  , getFFI2D
  , getFFI3D
    -- * Construction
  , fromSCMρ
  , fromSCMλ
  ) where

import Foreign                  (Ptr)
import Measure.States.Dense.Allocation
import Measure.States.Dense.Structure


-- |
-- /ϴ(1)/
--
-- Extract sub-component required for C FFI to perform alignment between /two/
-- sequences.
getFFI2D :: TCMρ -> Ptr FFI2D
getFFI2D = matrix2D


-- |
-- /ϴ(1)/
--
-- Extract sub-component required for C FFI to perform alignment between /three/
-- sequences.
getFFI3D :: TCMρ -> Ptr FFI3D
getFFI3D = matrix3D
