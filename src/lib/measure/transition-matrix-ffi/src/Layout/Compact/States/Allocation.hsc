-----------------------------------------------------------------------------
-- |
-- Module      :  Layout.Compact.States.Allocation
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-- Generate the 2D and 3D compact TCM matricies used for FFI calls.
--
-- For notes on usage, data construction and external see referenced C
-- compilation units, and also driver.c, which is not imported, but is
-- included indirectory for reference.
-----------------------------------------------------------------------------

{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE Strict                   #-}

#include "c_code_alloc_setup.h"

module Layout.Compact.States.Allocation
  ( -- * Construction
    initialize
--    fromSDMρ
--  , fromSDMλ
  ) where

import Data.Coerce
import Data.Vector.Storable (Vector)
import qualified Data.Vector.Storable as V
import Foreign
import Foreign.C.Types
import Layout.Compact.States.Structure
import System.IO.Unsafe (unsafePerformIO)


-- |
-- Create and allocate cost matrix first argument, TCM, is only for non-ambiguous
-- nucleotides, and it used to generate the entire cost matrix, which includes ambiguous elements.
-- TCM is row-major, with each row being the left character element.
-- It is therefore indexed not by powers of two, but by cardinal integer.
foreign import ccall unsafe "c_code_alloc_setup.h setUp2dCostMtx"

    initializeCostMatrix2D_FFI
      :: Ptr FFI2D
      -> Ptr CUShort -- ^ tcm (The SDM row-major vector)
      -> CUShort     -- ^ gap_open_cost
      -> CSize       -- ^ alphSize
      -> IO ()


foreign import ccall unsafe "c_code_alloc_setup.h setUp3dCostMtx"

    initializeCostMatrix3D_FFI
      :: Ptr FFI3D
      -> Ptr CUShort  -- ^ tcm (The SDM row-major vector)
      -> CUShort      -- ^ gap_open_cost
      -> CSize        -- ^ alphSize
      -> IO ()


{-
-- |
-- /ϴ(a⁵)/ where /a ≤ 8/ is the size of the character alphabet.
--
-- Generate the 2D and 3D dense transiton cost matricies ('TCMρ') from the
-- supplied symbol change matrix function ('SDMρ') with linear dimensions of
-- the alphabet symbol count.
--
-- Prefer 'fromSDMρ' to 'fromSDMλ', as this function performs /ϴ(a²)/ less
-- allocation than that function.
fromSDMρ
  :: SymbolChangeCost           -- ^ The gap open cost. A zero value indicates non-affine alignment context
  -> SymbolDistanceMatrixSquare -- ^ The dense, pre-computed matrix of costs to shift between symbols.
  -> TCMρ
fromSDMρ penalty scmρ =
    let (clip,size) = clampSize $ sizeOfSDM scmρ
        scmVector   = clip $ rowMajorVector scmρ
    in  initialize size penalty scmVector


-- |
-- /ϴ(a⁵)/ where /a ≤ 8/ is the size of the character alphabet.
--
-- Generate the 2D and 3D dense transiton cost matricies ('TCMρ') from the
-- supplied symbol change matrix function ('SDMλ') with linear dimensions of
-- the alphabet symbol count.
--
-- Prefer 'fromSDMρ' to 'fromSDMλ', as that function performs /ϴ(a²)/ less
-- allocation than this function.
fromSDMλ
  :: SymbolCount      -- ^ The character alphabet size
  -> SymbolChangeCost -- ^ The gap open cost. A zero value indicates non-affine alignment context
  -> SDMλ             -- ^ The function defining the costs to shift between symbols.
  -> TCMρ
fromSDMλ alphabetSize openningCost scmλ =
    let (_, size) = clampSize alphabetSize
        scmVector = coerce . rowMajorVector $ unsafeFromSDMλ scmλ size
    in  initialize size openningCost scmVector
-}


-- |
-- /ϴ(a⁵)/ where /a ≤ 8/ is the number of symbols in the character alphabet.
--
-- Generate the 2D and 3D compact state transition matricies from the supplied
-- row-major vecotr of symbol transtion distances linear dimensions of the
-- symbol count.
initialize
  :: Word -- ^ Number of symbols to allocate for
  -> Word -- ^ Penalty cost to begin a gap sequence
  -> Vector CUShort
  -> TCMρ
initialize dim penalty inputVector =
    let (clip,size) = clampSize dim
        safeVector  = clip inputVector
        openPenalty = fromIntegral penalty
        dimension   = fromIntegral size
        rowLen      = fromEnum     size
        firstRow    = V.slice  1 (rowLen - 1) safeVector
        firstCol    = V.generate (rowLen - 1) $ \i -> safeVector V.! ((i+1) * rowLen)
        maxDel      = V.maximum firstCol
        maxIns      = V.maximum firstRow
        minDel      = V.minimum firstCol
        minIns      = V.minimum firstRow
    in  unsafePerformIO . V.unsafeWith safeVector $ \arr -> do
          cm2D <- malloc :: IO (Ptr FFI2D)
          cm3D <- malloc :: IO (Ptr FFI3D)
          _    <- initializeCostMatrix2D_FFI cm2D arr openPenalty dimension
          _    <- initializeCostMatrix3D_FFI cm3D arr openPenalty dimension
          pure StateTransitionsCompact
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
-- Ensure that the size does not exceed 'maximumDimension'.
clampSize :: Word -> (Vector a -> Vector CUShort, Word)
clampSize n =
  let truncateVector :: Vector CUShort -> Vector CUShort
      truncateVector =
        let x = fromEnum maximumDimension
        in  V.slice 0 (x*x)
      
      transform :: Vector CUShort -> Vector CUShort
      transform
        | n <= maximumDimension = id
        | otherwise             = truncateVector

  in  (transform . coerce, min maximumDimension n)
