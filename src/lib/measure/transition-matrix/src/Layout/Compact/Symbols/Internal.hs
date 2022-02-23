-----------------------------------------------------------------------------
-- |
-- Module      :  Layout.Compact.Symbols.Internal
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE DeriveAnyClass        #-}
{-# LANGUAGE DeriveDataTypeable    #-}
{-# LANGUAGE DeriveGeneric         #-}
{-# LANGUAGE DerivingStrategies    #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE RankNTypes            #-}
{-# LANGUAGE Strict                #-}
{-# LANGUAGE TypeFamilies          #-}

module Layout.Compact.Symbols.Internal
  ( -- * Representational Type
    SymbolDistanceMatrix(..)
    -- * Queries
  , index
  , rowMajorVector
  , bytesSizeSymbolMatrix
  ) where

import           Control.DeepSeq
import           Data.Data
import           Data.Vector.Storable            (Vector)
import qualified Data.Vector.Storable            as V
import           Data.Word
import           Foreign.Storable
import           Foreign.Ptr                     (Ptr)
import           GHC.Generics
import           Measure.Range
import           Measure.Unit.SymbolChangeCost
import           Measure.Unit.SymbolCount
import           Measure.Unit.SymbolIndex


-- |
-- General symbol distance matrix layout.
data  SymbolDistanceMatrix
    = SymbolDistanceMatrix
        {-# UNPACK #-} !SymbolCount     -- ^ Matrix dimension
        {-# UNPACK #-} !(Vector Word16) -- ^ Contiguous array of matrix values
    deriving stock    (Data, Eq, Generic, Typeable)
    deriving anyclass (NFData)


instance MeasurableRange SymbolDistanceMatrix SymbolIndex where

    measureRange = symbolBounds . symbolCount


-- |
-- A structure which can derive the number of alphabet symbols associated with it.
instance HasSymbolCount SymbolDistanceMatrix where

    {-# INLINE symbolCount #-}
    symbolCount (SymbolDistanceMatrix n _) = n


-- |
-- /O(1)/
--
-- Indexing without bounds checking.
{-# INLINE index #-}
{-# SPECIALISE INLINE index :: SymbolDistanceMatrix -> Int -> SymbolChangeCost #-}
index :: Integral c => SymbolDistanceMatrix -> Int -> c
index (SymbolDistanceMatrix _ v) i =
    fromIntegral $ v `V.unsafeIndex` i


-- |
-- Deconstructs the 'SymbolDistanceMatrix' to expose the underlying unboxed 'Vector'.
{-# INLINE rowMajorVector #-}
rowMajorVector :: SymbolDistanceMatrix -> Vector Word16
rowMajorVector (SymbolDistanceMatrix _ v) = v


-- |
-- Computes the number of bytes used to store the symbol distance matrix.
{-# INLINE bytesSizeSymbolMatrix #-}
bytesSizeSymbolMatrix :: SymbolDistanceMatrix -> Int
bytesSizeSymbolMatrix (SymbolDistanceMatrix _ vec) =
    let c = sizeOf (undefined :: SymbolCount)
        p = sizeOf (undefined :: Ptr Word)
        i = sizeOf (undefined :: Ptr Int )
        e = sizeOf (undefined :: Word16)
        n = V.length vec
    in  sum
        [ c
        , p + p
        , i
        , e * n
        ]
       
