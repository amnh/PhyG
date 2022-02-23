-----------------------------------------------------------------------------
-- |
-- Module      :  Layout.Memoize.States
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE DeriveAnyClass        #-}
{-# LANGUAGE DeriveGeneric         #-}
{-# LANGUAGE DerivingStrategies    #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE Strict                #-}

module Layout.Memoize.States
  ( Sparse(..)
  , initialize
  ) where

import           Control.DeepSeq
import           Data.Bits
import           Data.Foldable (fold)
import           Data.Hashable
import           Data.Word                                          (Word64)
import           Foreign.C.Types (CUInt)
import           Foreign.Ptr (Ptr)
import           Foreign.Storable (sizeOf)
import           GHC.Generics
import           Layout.Compact.Symbols.Square     (SymbolDistanceMatrixSquare    , bytesSizeMatrixSquare    )
import           Layout.Compact.Symbols.Triangular (SymbolDistanceMatrixTriangular, bytesSizeMatrixTriangular)
import           Layout.Memoize.Dispersion
import           Layout.Memoize.Hashable
import           Measure.Range
import           Measure.Transition
import           Measure.Unit.SymbolChangeCost
import           Measure.Unit.SymbolCount
import           Measure.Unit.SymbolIndex


data  Sparse a
    = Sparse
    { maxDelCost  :: {-# UNPACK #-} !SymbolChangeCost
    , maxInsCost  :: {-# UNPACK #-} !SymbolChangeCost
    , minDelCost  :: {-# UNPACK #-} !SymbolChangeCost
    , minInsCost  :: {-# UNPACK #-} !SymbolChangeCost 
    , matrixForm  :: {-# UNPACK #-} !(Either SymbolDistanceMatrixTriangular SymbolDistanceMatrixSquare)
    , memoized2Dλ :: !(TCM2Dλ a)
    , memoized3Dλ :: !(TCM3Dλ a)
    }
    deriving stock    (Generic)
    deriving anyclass (NFData)


instance Eq (Sparse a) where

    (==) lhs rhs = and
        [ maxDelCost lhs == maxDelCost rhs
        , maxInsCost lhs == maxInsCost rhs
        , minDelCost lhs == minDelCost rhs
        , minInsCost lhs == minInsCost rhs
        , matrixForm lhs == matrixForm rhs
        ]


instance HasEditExtrema (Sparse a) where

    maxDeletion  = maxDelCost

    maxInsertion = maxInsCost

    minDeletion  = minDelCost

    minInsertion = minInsCost


instance HasSymbolCount (Sparse a) where

    symbolCount = either symbolCount symbolCount . matrixForm


instance HasSymbolDistances (Sparse a) where

    symbolDistances = either symbolDistances symbolDistances . matrixForm


instance HasStateTransitions (Sparse a) a where

    {-# INLINE            stateTransitionPairwiseDispersion #-}
    {-# SPECIALISE INLINE stateTransitionPairwiseDispersion :: Sparse CUInt  -> TCM2Dλ CUInt  #-}
    {-# SPECIALISE INLINE stateTransitionPairwiseDispersion :: Sparse Word64 -> TCM2Dλ Word64 #-}
    stateTransitionPairwiseDispersion = memoized2Dλ

    {-# INLINE            stateTransitionThreewayDispersion #-}
    {-# SPECIALISE INLINE stateTransitionThreewayDispersion :: Sparse CUInt  -> TCM3Dλ CUInt  #-}
    {-# SPECIALISE INLINE stateTransitionThreewayDispersion :: Sparse Word64 -> TCM3Dλ Word64 #-}
    stateTransitionThreewayDispersion = memoized3Dλ


instance HasStateTransitionsCompact (Sparse a) where

    {-# INLINE getCompactPairwise #-}
    getCompactPairwise = const Nothing

    {-# INLINE getCompactThreeway #-}
    getCompactThreeway = const Nothing


instance MeasurableRange (Sparse a) SymbolIndex where

    {-# INLINE measureRange #-}
    measureRange = either measureRange measureRange . matrixForm


instance Show (Sparse a) where

    show = renderSummary


{-# INLINEABLE        initialize #-}
{-# SPECIALISE INLINE initialize :: Either SymbolDistanceMatrixTriangular SymbolDistanceMatrixSquare -> Sparse Word64 #-}
initialize
  :: ( FiniteBits b
     , Hashable b
     , NFData b
     )
  => Either SymbolDistanceMatrixTriangular SymbolDistanceMatrixSquare
  -> Sparse b
initialize eMatrix =
    let range = either measureRange    measureRange    eMatrix
        sdmλ  = either symbolDistances symbolDistances eMatrix
        tcm2D = memoize2 $ bitDispersionPairwise range sdmλ
        tcm3D = memoize3 $ bitDispersionThreeway range sdmλ
    in  Sparse
        { maxDelCost  = either maxDeletion  maxDeletion  eMatrix
        , maxInsCost  = either maxInsertion maxInsertion eMatrix
        , minDelCost  = either minDeletion  minDeletion  eMatrix
        , minInsCost  = either minInsertion minInsertion eMatrix
        , matrixForm  = eMatrix
        , memoized2Dλ = tcm2D
        , memoized3Dλ = tcm3D
        }


renderSummary :: Sparse a -> String
renderSummary tcm =
    let b = bytesSizeOfSparse tcm
        n = symbolCount tcm
    in  unlines
        [ "General Metric (Memoized) with"
        , fold [ "  allocated "     , show b, "bytes + memoized value space" ]
        , fold [ "  dimension "     , show n ]
        ]


bytesSizeOfSparse :: Sparse a -> Int
bytesSizeOfSparse tms =
    let bytesPtr    = sizeOf (undefined :: Ptr Word)
        bytesWord   = sizeOf (undefined :: Word)
        bytesMatrix = either bytesSizeMatrixTriangular bytesSizeMatrixSquare $ matrixForm tms
    in  sum
        [ 4 * bytesWord       -- 4 unpacked fields
        , bytesPtr            -- Point to Either
        , bytesPtr + bytesPtr -- Either points to Left / Right
        , bytesMatrix         -- Space behind Either's pointer
        , bytesPtr            -- Point to memoized 2D
        , bytesPtr            -- Point to memoized 3D
        ]
