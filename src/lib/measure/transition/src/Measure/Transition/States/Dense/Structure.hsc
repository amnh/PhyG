{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE DeriveAnyClass           #-}
{-# LANGUAGE DeriveGeneric            #-}
{-# LANGUAGE DerivingStrategies       #-}
{-# LANGUAGE Strict                   #-}

#include "costMatrix.h"
#include "alignmentMatrices.h"

module Measure.Transition.States.Dense.Structure
  ( AlignmentStrategy(..)
  , TCMρ(..)
  , FFI2D()
  , FFI3D()
  -- * Accessor
  , getAlignmentStrategy
  -- * Indexing
  , getTCMρ2D
  , getTCMρ3D
  ) where

import Control.DeepSeq
import Data.Foldable
import Foreign
import Foreign.C.Types
import GHC.Generics     (Generic)
import System.IO.Unsafe (unsafePerformIO)
import Measure.Unit.SymbolChangeCost
import Measure.Unit.SymbolCount


-- |
-- Specify which alignment to perform
data  AlignmentStrategy = Linear | Affine
    deriving stock    (Eq, Generic, Show)
    deriving anyclass (NFData)


-- |
-- Holds single cost matrix, which contains costs and medians for all possible
-- character elements. It is completely filled using a TCM.
data  FFI2D
    = FFI2D
    { alphSize2D            :: CInt
    -- ^ alphabet size including gap, and including ambiguities if
    , costMatrixDimension2D :: CInt
    -- ^ ceiling of log_2 (alphSize)
    , gapChar2D             :: CInt
    -- ^ gap value (least significant bit set)
    , costModelType2D       :: CInt
    {- ^ The type of cost model to be used in the alignment, i.e. affine or not.
         Based on cost_matrix.ml, values are:
           • linear == 0
           • no_alignment == 2
           • affine == 3
         but I updated it. See costMatrix.h.
    -}
    , include_ambiguities2D :: CInt
    {- ^ This is a flag set to true if we are going to accept all possible
         combinations of the elements in the alphabet in the alignments. This is
         not true for protein characters for example, where the number of elements
         of the alphabet is already too big to build all the possible combinations.
    -}
    , gapOpenCost2D         :: CInt
    {- ^ The cost of opening a gap. This is only useful in certain cost_model_types
         (type 3: affine, based on my reading of ML code).
    -}
    , isMetric2D            :: CInt
    -- ^ if tcm is metric
    , allElems2D            :: CInt
    -- ^ total number of elements
    , bestCost2D            :: Ptr CUInt
    {- ^ The transformation cost matrix, including ambiguities, storing the **best**
         cost for each ambiguity pair
    -}
    , medians2D             :: Ptr CUInt
    {- ^ The matrix of possible medians between elements in the alphabet. The best
         possible medians according to the cost matrix.
    -}
    , worstCost2D           :: Ptr CInt
    {- ^ The transformation cost matrix, including ambiguities, storing the **worst**
         cost for each ambiguity pair
    -}
    , prependCost2D         :: Ptr CInt
    {- ^ The cost of going from gap -> each base. For ambiguities, use best cost.
         Set up as num_elements x num_elements matrix, but seemingly only first row is used.
    -}
    , tailCost2D            :: Ptr CInt
    {- ^ As prepend_cost, but with reverse directionality, so base -> gap.
         As with prepend_cost, seems to be allocated as too large.
    -}
    }
    deriving stock    (Eq, Generic)
    deriving anyclass (NFData)


-- |
-- A representation of the 3D cost matrix structure used on the C side.
data  FFI3D
    = FFI3D          -- See FFI2D datatype for field description
    { alphSize3D            :: CInt
    , costMatrixDimension3D :: CInt
    , gapChar3D             :: CInt
    , costModelType3D       :: CInt
    , include_ambiguities3D :: CInt
    , gapOpenCost3D         :: CInt
    , allElems3D            :: CInt
    , bestCost3D            :: Ptr CUInt
    , medians3D             :: Ptr CUInt
    }
    deriving stock    (Eq, Generic)
    deriving anyclass (NFData)


-- |
-- Exposed wrapper for C allocated cost matrix structures.
--
-- The 'matrix2D' and 'matrix3D' components are required for string alignment
-- utilizing the C FFI.
--
-- Use 'getSCMλ', 'getTCM2Dλ', and 'getTCM3Dλ' to extract symbol ⨉ symbol,
-- state ⨉ state, and state ⨉ state ⨉ state measures, respectively.
data  TCMρ
    = TCMρ
    { gapPenalty :: {-# UNPACK #-} SymbolChangeCost
    , maxDelCost :: {-# UNPACK #-} SymbolChangeCost
    , maxInsCost :: {-# UNPACK #-} SymbolChangeCost
    , minDelCost :: {-# UNPACK #-} SymbolChangeCost
    , minInsCost :: {-# UNPACK #-} SymbolChangeCost
    , matrixSize :: {-# UNPACK #-} SymbolCount
    , matrix2D   :: {-# UNPACK #-} (Ptr FFI2D)
    , matrix3D   :: {-# UNPACK #-} (Ptr FFI3D)
    }
    deriving stock    (Eq, Generic)
    deriving anyclass (NFData)


instance Enum AlignmentStrategy where

    fromEnum Linear = 0
    fromEnum Affine = 3

    toEnum 3 = Affine
    toEnum _ = Linear


{-
instance HasEditCosts TCMρ where

    maxDeletion  = maxDelCost

    maxInsertion = maxInsCost

    minDeletion  = minDelCost

    minInsertion = minInsCost
-}


instance HasSymbolCount TCMρ where

    symbolCount = matrixSize


instance Show FFI2D where

    show = unlines . (fieldRendering <*>)  . pure
      where
        fieldRendering =
            [ show .            alphSize2D
            , show . costMatrixDimension2D
            , show .             gapChar2D
            , show .       costModelType2D
            , show . include_ambiguities2D
            , show .         gapOpenCost2D
            , show .            isMetric2D
            , show .            allElems2D
            , show .            bestCost2D
            , show .             medians2D
            , show .           worstCost2D
            , show .         prependCost2D
            , show .            tailCost2D
            ]


instance Show FFI3D where

    show = unlines . (fieldRendering <*>)  . pure
      where
        fieldRendering =
            [ show .            alphSize3D
            , show . costMatrixDimension3D
            , show .             gapChar3D
            , show .       costModelType3D
            , show . include_ambiguities3D
            , show .         gapOpenCost3D
            , show .            allElems3D
            , show .            bestCost3D
            , show .             medians3D
            ]


instance Show TCMρ where

    show m =
        let n = matrixSize m
        in  fold [ "TCMρ "
                 ,  show n
                 , "X"
                 ,  show n
                 , " Matrix"
                 ]


instance Storable FFI2D where

    sizeOf = const $ #size struct cost_matrices_2d_t

    alignment = sizeOf -- alignment (undefined :: StablePtr FFI2D)

    peek ptr = FFI2D
        <$> (#peek struct cost_matrices_2d_t, alphSize           ) ptr
        <*> (#peek struct cost_matrices_2d_t, costMatrixDimension) ptr
        <*> (#peek struct cost_matrices_2d_t, gap_char           ) ptr
        <*> (#peek struct cost_matrices_2d_t, cost_model_type    ) ptr
        <*> (#peek struct cost_matrices_2d_t, include_ambiguities) ptr
        <*> (#peek struct cost_matrices_2d_t, gap_open_cost      ) ptr
        <*> (#peek struct cost_matrices_2d_t, is_metric          ) ptr
        <*> (#peek struct cost_matrices_2d_t, num_elements       ) ptr
        <*> (#peek struct cost_matrices_2d_t, cost               ) ptr
        <*> (#peek struct cost_matrices_2d_t, median             ) ptr
        <*> (#peek struct cost_matrices_2d_t, worst              ) ptr
        <*> (#peek struct cost_matrices_2d_t, prepend_cost       ) ptr
        <*> (#peek struct cost_matrices_2d_t, tail_cost          ) ptr

    poke ptr cm2D = do -- to modify values in the C app
        (#poke struct cost_matrices_2d_t, alphSize           ) ptr $            alphSize2D cm2D
        (#poke struct cost_matrices_2d_t, costMatrixDimension) ptr $ costMatrixDimension2D cm2D
        (#poke struct cost_matrices_2d_t, gap_char           ) ptr $             gapChar2D cm2D
        (#poke struct cost_matrices_2d_t, cost_model_type    ) ptr $       costModelType2D cm2D
        (#poke struct cost_matrices_2d_t, include_ambiguities) ptr $ include_ambiguities2D cm2D
        (#poke struct cost_matrices_2d_t, gap_open_cost      ) ptr $         gapOpenCost2D cm2D
        (#poke struct cost_matrices_2d_t, is_metric          ) ptr $            isMetric2D cm2D
        (#poke struct cost_matrices_2d_t, num_elements       ) ptr $             medians2D cm2D
        (#poke struct cost_matrices_2d_t, cost               ) ptr $            bestCost2D cm2D
        (#poke struct cost_matrices_2d_t, median             ) ptr $             medians2D cm2D
        (#poke struct cost_matrices_2d_t, worst              ) ptr $           worstCost2D cm2D
        (#poke struct cost_matrices_2d_t, prepend_cost       ) ptr $         prependCost2D cm2D
        (#poke struct cost_matrices_2d_t, tail_cost          ) ptr $            tailCost2D cm2D


instance Storable FFI3D where

    sizeOf = const $ #size struct cost_matrices_2d_t

    alignment = sizeOf -- alignment (undefined :: StablePtr FFI2D)

    peek ptr = FFI3D
        <$> (#peek struct cost_matrices_3d_t, alphSize           ) ptr
        <*> (#peek struct cost_matrices_3d_t, costMatrixDimension) ptr
        <*> (#peek struct cost_matrices_3d_t, gap_char           ) ptr
        <*> (#peek struct cost_matrices_3d_t, cost_model_type    ) ptr
        <*> (#peek struct cost_matrices_3d_t, include_ambiguities) ptr
        <*> (#peek struct cost_matrices_3d_t, gap_open_cost      ) ptr
        <*> (#peek struct cost_matrices_3d_t, num_elements       ) ptr
        <*> (#peek struct cost_matrices_3d_t, cost               ) ptr
        <*> (#peek struct cost_matrices_3d_t, median             ) ptr

    poke ptr cm3D = do -- to modify values in the C app
        (#poke struct cost_matrices_3d_t, alphSize           ) ptr $            alphSize3D cm3D
        (#poke struct cost_matrices_3d_t, costMatrixDimension) ptr $ costMatrixDimension3D cm3D 
        (#poke struct cost_matrices_3d_t, gap_char           ) ptr $             gapChar3D cm3D
        (#poke struct cost_matrices_3d_t, cost_model_type    ) ptr $       costModelType3D cm3D
        (#poke struct cost_matrices_3d_t, include_ambiguities) ptr $ include_ambiguities3D cm3D
        (#poke struct cost_matrices_3d_t, gap_open_cost      ) ptr $         gapOpenCost3D cm3D
        (#poke struct cost_matrices_3d_t, num_elements       ) ptr $            allElems3D cm3D
        (#poke struct cost_matrices_3d_t, cost               ) ptr $            bestCost3D cm3D
        (#poke struct cost_matrices_3d_t, median             ) ptr $             medians3D cm3D


-- |
-- /ϴ(1)/
--
-- Determine the alignment strategy encoded for the matrix.
getAlignmentStrategy :: FFI2D -> AlignmentStrategy
getAlignmentStrategy = toEnum . fromEnum . costModelType2D


-- |
-- /O(1)/
{-# INLINEABLE getTCMρ2D #-}
{-# SPECIALISE INLINE getTCMρ2D :: TCMρ -> CUInt -> CUInt -> (SymbolChangeCost, CUInt) #-}
getTCMρ2D
  :: ( Enum a
     , Enum b
     , Enum c
     , Num d
     )
  => TCMρ
  -> a
  -> b
  -> (d, c)
getTCMρ2D tcmρ e1 e2 = unsafePerformIO $ do
    cm2D <- peek $ matrix2D tcmρ
    let dim = fromEnum $ alphSize2D cm2D
    let off = (fromEnum e1 `shiftL` dim) + fromEnum e2
    let get = peek . flip advancePtr off
    cost <- get $ bestCost2D cm2D
    med  <- get $  medians2D cm2D
    pure (fromIntegral cost, toEnum $ fromEnum med)


-- |
-- /O(1)/
{-# INLINEABLE getTCMρ3D #-}
{-# SPECIALISE INLINE getTCMρ3D :: TCMρ -> CUInt -> CUInt -> CUInt -> (SymbolChangeCost, CUInt) #-}
getTCMρ3D
  :: ( Enum a
     , Enum b
     , Enum c
     , Enum e
     , Num d
     )
  => TCMρ
  -> a
  -> b
  -> c
  -> (d, e)
getTCMρ3D tcmρ e1 e2 e3 = unsafePerformIO $ do
    cm3D <- peek $ matrix3D tcmρ
    let dim = fromEnum $ alphSize3D cm3D
    let off = (((fromEnum e1 `shiftL` dim) + fromEnum e2) `shiftL` dim) + fromEnum e3
    let get = peek . flip advancePtr off
    cost <- get $ bestCost3D cm3D
    med  <- get $  medians3D cm3D
    pure (fromIntegral cost, toEnum $ fromEnum med)

