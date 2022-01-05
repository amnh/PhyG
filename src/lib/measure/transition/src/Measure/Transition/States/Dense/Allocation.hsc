-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Transition.States.Dense.Allocation
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-- Generate the 2D and 3D dense TCM matricies used for FFI calls.
--
-- For notes on usage, data construction and external see referenced C
-- compilation units, and also driver.c, which is not imported, but is
-- included indirectory for reference.
-----------------------------------------------------------------------------

{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE Strict                   #-}

#include "c_code_alloc_setup.h"

module Measure.Transition.States.Dense.Allocation
  ( -- * Construction
    fromSCMρ
  , fromSCMλ
  ) where

import Data.Coerce
import Data.Vector.Storable (Vector)
import qualified Data.Vector.Storable as V
import Foreign
import Foreign.C.Types
import Measure.Transition.Symbols.Internal
import Measure.Transition.States.Dense.Structure
import Measure.Transition.Symbols.Unsafe
import Measure.Transition.Matrix
import Measure.Unit.Distance
import Measure.Unit.SymbolCount
import System.IO.Unsafe (unsafePerformIO)


-- |
-- Create and allocate cost matrix first argument, TCM, is only for non-ambiguous
-- nucleotides, and it used to generate the entire cost matrix, which includes ambiguous elements.
-- TCM is row-major, with each row being the left character element.
-- It is therefore indexed not by powers of two, but by cardinal integer.
foreign import ccall unsafe "c_code_alloc_setup.h setUp2dCostMtx"

    initializeCostMatrix2D_FFI
      :: Ptr FFI2D
      -> Ptr CUShort -- ^ tcm (The SCM row-major vector)
      -> CUShort     -- ^ gap_open_cost
      -> CSize       -- ^ alphSize
      -> IO ()


foreign import ccall unsafe "c_code_alloc_setup.h setUp3dCostMtx"

    initializeCostMatrix3D_FFI
      :: Ptr FFI3D
      -> Ptr CUShort  -- ^ tcm (The SCM row-major vector)
      -> CUShort      -- ^ gap_open_cost
      -> CSize        -- ^ alphSize
      -> IO ()


-- |
-- /ϴ(a⁵)/ where /a ≤ 8/ is the size of the character alphabet.
--
-- Generate the 2D and 3D dense transiton cost matricies ('TCMρ') from the
-- supplied symbol change matrix function ('SCMρ') with linear dimensions of
-- the alphabet symbol count.
--
-- Prefer 'fromSCMρ' to 'fromSCMλ', as this function performs /ϴ(a²)/ less
-- allocation than that function.
fromSCMρ
  :: Distance -- ^ The gap open cost. A zero value indicates non-affine alignment context
  -> SCMρ     -- ^ The dense, pre-computed matrix of costs to shift between symbols.
  -> TCMρ
fromSCMρ penalty scmρ =
    let (clip,size) = clampSize $ sizeOfSCM scmρ
        scmVector   = clip $ rowMajorVector scmρ
    in  initialize size penalty scmVector


-- |
-- /ϴ(a⁵)/ where /a ≤ 8/ is the size of the character alphabet.
--
-- Generate the 2D and 3D dense transiton cost matricies ('TCMρ') from the
-- supplied symbol change matrix function ('SCMλ') with linear dimensions of
-- the alphabet symbol count.
--
-- Prefer 'fromSCMρ' to 'fromSCMλ', as that function performs /ϴ(a²)/ less
-- allocation than this function.
fromSCMλ
  :: SymbolCount -- ^ The character alphabet size
  -> Distance    -- ^ The gap open cost. A zero value indicates non-affine alignment context
  -> SCMλ        -- ^ The function defining the costs to shift between symbols.
  -> TCMρ
fromSCMλ alphabetSize openningCost scmλ =
    let (_, size) = clampSize alphabetSize
        scmVector = coerce . rowMajorVector $ unsafeFromSCMλ scmλ size
    in  initialize size openningCost scmVector


initialize
  :: SymbolCount -- ^ Number of symbols to allocate for
  -> Distance    -- ^ Penalty cost to begin a gap sequence
  -> Vector CUShort -> TCMρ
initialize size penalty scmVector =
    let openPenalty = fromIntegral penalty
        dimension   = fromIntegral size
        rowLen      = fromIntegral size
        firstRow    = V.slice  1 (rowLen - 1) scmVector
        firstCol    = V.generate (rowLen - 1) $ \i -> scmVector V.! ((i+1) * rowLen)
        maxDel      = V.maximum firstCol
        maxIns      = V.maximum firstRow
        minDel      = V.minimum firstCol
        minIns      = V.minimum firstRow
    in  unsafePerformIO . V.unsafeWith scmVector $ \arr -> do
          cm2D <- malloc :: IO (Ptr FFI2D)
          cm3D <- malloc :: IO (Ptr FFI3D)
          _    <- initializeCostMatrix2D_FFI cm2D arr openPenalty dimension
          _    <- initializeCostMatrix3D_FFI cm3D arr openPenalty dimension
          pure TCMρ
            { gapPenalty = penalty
            , maxDelCost = fromIntegral maxDel
            , maxInsCost = fromIntegral maxIns
            , minDelCost = fromIntegral minDel
            , minInsCost = fromIntegral minIns
            , matrixSize = size
            , matrix2D   = cm2D
            , matrix3D   = cm3D
            }


-- |
-- Ensure that the size does not exceed 'maximumSize'.
clampSize :: SymbolCount -> (Vector a -> Vector CUShort, SymbolCount)
clampSize n =
  let truncateVector :: Vector CUShort -> Vector CUShort
      truncateVector = let x = fromEnum infimumSymbolLimit in V.slice 0 (x*x)
      
      transform :: Vector CUShort -> Vector CUShort
      transform
        | iota n    = id
        | otherwise = truncateVector

  in  (transform . coerce, min infimumSymbolLimit n)
