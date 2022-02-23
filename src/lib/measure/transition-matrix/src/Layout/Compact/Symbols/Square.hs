-----------------------------------------------------------------------------
-- |
-- Module      :  Layout.Compact.Symbols.Square
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
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE RankNTypes            #-}
{-# LANGUAGE Strict                #-}
{-# LANGUAGE TypeFamilies          #-}

module Layout.Compact.Symbols.Square
  ( -- * Representational Type
    SymbolDistanceMatrixSquare(..)
  , SDMS
    -- * Queries
  , bytesSizeMatrixSquare
  , rowMajorVector
  ) where

import           Control.DeepSeq
import           Data.Coerce
import           Data.Data
import           Data.Foldable
import           Data.Ix                         (Ix(range))
import           Data.Ord
import           Data.Vector.Storable            (Vector)
import qualified Data.Vector.Storable            as V
import           Data.Word
import           GHC.Generics
import           Measure.Range
import           Measure.Transition
import           Measure.Unit.SymbolCount
import           Measure.Unit.SymbolIndex
import           Layout.Compact.Symbols.Internal (SymbolDistanceMatrix(..))
import qualified Layout.Compact.Symbols.Internal as SDM


-- |
-- A data structure for storing a square array of sizeality
-- greater than or equal to two, with positive cost values at the array indices.
-- Values are stored in an unboxed structure for cache efficiency.
--
-- A 'SDMS' can be constructed by calling one of the following functions:
--
--   * 'fromList'
--
--   * 'fromCols'
--
--   * 'fromRows'
--
--   * 'generate'
--
-- Attempts to construct an empty or singleton 'SDMS' through the above
-- constructors will result in a runtime exception.
newtype SymbolDistanceMatrixSquare = SDMS { sdm :: SymbolDistanceMatrix }
    deriving newtype (Eq, NFData, Typeable)
    deriving stock   (Data, Generic)


type SDMS = SymbolDistanceMatrixSquare


-- |
-- Any structural representation which can produce a Symbol Change Matrix.
instance HasEditExtrema SymbolDistanceMatrixSquare where

    {-# INLINEABLE maxEdit #-}
    maxEdit a = max (maxInsertion a) $ maxDeletion a

    {-# INLINEABLE minEdit #-}
    minEdit a = min (minInsertion a) $ minDeletion a

    maxDeletion  = firstColExtrama maximumBy

    minDeletion  = firstColExtrama minimumBy

    maxInsertion = firstRowExtrema V.maximum

    minInsertion = firstRowExtrema V.minimum


-- |
-- Any structural representation which can produce a Symbol Change Matrix.
instance HasSymbolDistances SymbolDistanceMatrixSquare where

    {-# INLINE symbolDistances #-}
    symbolDistances = symbolIndexing


-- |
-- A structure which can derive the number of alphabet symbols associated with it.
instance HasSymbolCount SymbolDistanceMatrixSquare where

    {-# INLINE symbolCount #-}
    symbolCount = symbolCount . sdm


instance MeasurableRange SymbolDistanceMatrixSquare SymbolIndex where

    {-# INLINE measureRange #-}
    measureRange = measureRange . sdm


-- |
-- A pretty printed custom show instance for Symbol Distance Matrix.
instance Show SymbolDistanceMatrixSquare where

    show sdms = headerLine <> matrixLines
      where
        renderRow i = ("  "<>) . unwords $ renderValue <$> [ symbolIndexing sdms i j | j <- rangeValues ]
        matrixLines = unlines $ renderRow <$> rangeValues
        rangeValues = range $ measureRange sdms
        headerLine  = '\n' : unwords [ "SDM:", show $ symbolCount sdms, "x", show $ symbolCount sdms, "\n"]
        maxValue    = V.maximum $ rowMajorVector sdms
        padSpacing  = length $ show maxValue
        renderValue x = pad <> shown
          where
            shown = show x
            pad   = (padSpacing - length shown) `replicate` ' '


-- |
-- /O(1)/
--
-- Indexing without bounds checking.
{-# INLINE symbolIndexing #-}
symbolIndexing :: SymbolDistanceMatrixSquare -> SymbolDistanceλ
symbolIndexing sdms i j =
    let m  = sdm sdms
        n' = coerce (symbolCount m) :: Word
        i' = coerce i :: Word
        j' = coerce j :: Word
        p  = fromEnum $ i' * n' + j'
    in  SDM.index m p


-- |
-- /O(1)/
--
-- Computes the number of bytes used to store the 'SymbolDistanceMatrixSquare'.
{-# INLINE bytesSizeMatrixSquare #-}
bytesSizeMatrixSquare :: SymbolDistanceMatrixSquare -> Int
bytesSizeMatrixSquare = SDM.bytesSizeSymbolMatrix . sdm


-- |
-- Deconstructs the 'SymbolDistanceMatrixSquare' to expose the underlying unboxed 'Vector'.
{-# INLINE rowMajorVector #-}
rowMajorVector :: SymbolDistanceMatrixSquare -> Vector Word16
rowMajorVector = SDM.rowMajorVector . sdm


firstRowExtrema :: (Integral b, Num c) => (Vector Word16 -> b) -> SymbolDistanceMatrixSquare -> c
firstRowExtrema f (SDMS (SymbolDistanceMatrix (SymbolCount n) v)) =
    fromIntegral . f $ V.slice 0 (fromEnum n) v


firstColExtrama :: (Integral a, Num b) => ((Int -> Int -> Ordering) -> [Int] -> a) -> SymbolDistanceMatrixSquare -> b
firstColExtrama f (SDMS (SymbolDistanceMatrix (SymbolCount n) v)) =
    let r = fromEnum n
        g i = v V.! (i*r)
    in  fromIntegral $ f (comparing g) [ 0 .. r - 1 ]
