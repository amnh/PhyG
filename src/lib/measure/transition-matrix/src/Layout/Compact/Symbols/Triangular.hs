-----------------------------------------------------------------------------
-- |
-- Module      :  Layout.Compact.Symbols.Triangular
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

module Layout.Compact.Symbols.Triangular
  ( -- * Representational Type
    SymbolDistanceMatrixTriangular(..)
  , SDMT
    -- * Queries
  , bytesSizeMatrixTriangular
  , rowMajorVector
    -- * Constructors
  , fromList
  , fromCols
  , fromRows
  , generate
  , lowerTriangleOfSquare
  ) where

import           Control.Arrow                   ((***))
import           Control.DeepSeq
import           Data.Bits (shiftR)
import           Data.Coerce
import           Data.Data
import           Data.Foldable
import           Data.Ix                         (Ix(range))
import           Data.List                       (transpose)
import           Data.List.Utility               (equalityOf, occurrences)
import           Data.Map                        (delete, findMax, keys)
import qualified Data.Map                        as Map (fromList)
import           Data.Ord
import           Data.Ratio
import           Data.Vector.Storable            (Vector)
import qualified Data.Vector.Storable            as V
import           Data.Word
import           GHC.Generics
import           Measure.Range
import           Measure.Transition
import           Measure.Unit.SymbolCount
import           Measure.Unit.SymbolIndex
import           Layout.Compact.Symbols.Square   (SymbolDistanceMatrixSquare)
import qualified Layout.Compact.Symbols.Square   as SDMS
import           Layout.Compact.Symbols.Internal (SymbolDistanceMatrix(..))
import qualified Layout.Compact.Symbols.Internal as SDM


-- |
-- A data structure for storing a square array of sizeality
-- greater than or equal to two, with positive cost values at the array indices.
-- Values are stored in an unboxed structure for cache efficiency.
--
-- A 'SDMT' can be constructed by calling one of the following functions:
--
--   * 'fromList'
--
--   * 'fromCols'
--
--   * 'fromRows'
--
--   * 'generate'
--
-- Attempts to construct an empty or singleton 'SDMT' through the above
-- constructors will result in a runtime exception.
newtype SymbolDistanceMatrixTriangular = SDMT { sdm :: SymbolDistanceMatrix }
    deriving newtype (Eq, NFData, Typeable)
    deriving stock   (Data, Generic)


type SDMT = SymbolDistanceMatrixTriangular


-- |
-- Any structural representation which can produce a Symbol Change Matrix.
instance HasEditExtrema SymbolDistanceMatrixTriangular where

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
instance HasSymbolDistances SymbolDistanceMatrixTriangular where

    {-# INLINE symbolDistances #-}
    symbolDistances = symbolIndexing


-- |
-- A structure which can derive the number of alphabet symbols associated with it.
instance HasSymbolCount SymbolDistanceMatrixTriangular where

    {-# INLINE symbolCount #-}
    symbolCount = symbolCount . sdm


instance MeasurableRange SymbolDistanceMatrixTriangular SymbolIndex where

    {-# INLINE measureRange #-}
    measureRange = measureRange . sdm


-- |
-- A pretty printed custom show instance for Symbol Distance Matrix.
instance Show SymbolDistanceMatrixTriangular where

    show sdmt = headerLine <> matrixLines
      where
        renderRow i = ("  "<>) . unwords $ renderValue <$> [ symbolIndexing sdmt i j | j <- rangeValues ]
        matrixLines = unlines $ renderRow <$> rangeValues
        rangeValues = range $ measureRange sdmt
        headerLine  = '\n' : unwords [ "SDM:", show $ symbolCount sdmt, "x", show $ symbolCount sdmt, "\n"]
        maxValue    = V.maximum $ rowMajorVector sdmt
        padSpacing  = length $ show maxValue
        renderValue x = pad <> shown
          where
            shown = show x
            pad   = (padSpacing - length shown) `replicate` ' '


lowerTriangleOfSquare :: SymbolDistanceMatrixSquare -> SymbolDistanceMatrixTriangular
lowerTriangleOfSquare square =
    let (SymbolDistanceMatrix sc@(SymbolCount w) vec) = coerce square
        dim = fromEnum w
        len = dim * (dim + 1) `shiftR` 1
        arr =   let r = [ 0 .. dim - 1 ]
                    p = [ V.unsafeIndex vec (dim * i + j) | i <- r, j <- r, i >= j ] 
                in  V.fromListN len p
    in  SDMT $ SymbolDistanceMatrix sc arr


-- |
-- /O(1)/
--
-- Indexing without bounds checking.
{-# INLINE symbolIndexing #-}
symbolIndexing :: SymbolDistanceMatrixTriangular -> SymbolDistanceλ
symbolIndexing sdmt i j | i < j = symbolIndexing sdmt j i
symbolIndexing sdmt i j =
    let i' = coerce i :: Word
        j' = coerce j :: Word
        t  = i' * (i' + 1) `shiftR` 1
        p  = fromEnum $ t + j'
    in  SDM.index (sdm sdmt) p


-- |
-- /O(1)/
--
-- Computes the number of bytes used to store the 'SymbolDistanceMatrixTriangular'.
{-# INLINE bytesSizeMatrixTriangular #-}
bytesSizeMatrixTriangular :: SymbolDistanceMatrixTriangular -> Int
bytesSizeMatrixTriangular = SDM.bytesSizeSymbolMatrix . sdm


-- |
-- Deconstructs the 'SymbolDistanceMatrixTriangular' to expose the underlying unboxed 'Vector'.
{-# INLINE rowMajorVector #-}
rowMajorVector :: SymbolDistanceMatrixTriangular -> Vector Word16
rowMajorVector = SDM.rowMajorVector . sdm


-- |
-- /O(n*n)/
--
-- Construct a 'SDM' from a list of elements in row major order.
--
-- ==== __Examples__
--
-- >>> fromList [1..9]
-- SDM: 3 x 3
--   1 2 3
--   4 5 6
--   7 8 9
--
-- >>> fromList []
-- *** Exception: fromList: An empty structure was supplied. Cannot construct an empty SDM!
--
-- >>> fromList [42]
-- *** Exception: fromList: A singleton structure was supplied. Cannot construct a SDM with size of 1, must have size of 2 or greater.
--
-- >>> fromList [1..12]
-- *** Exception: fromList: The number of element (12) is not a square number. Cannot construct an non-square SDM! The number of elements (12) lies between the valid square numbers (9) and (16).
--
fromList :: (Foldable t, Real a) => t a -> (Rational, SymbolDistanceMatrixTriangular)
fromList xs
  | null xs       = error $ nullarySizeMessage "fromList"
  | notTriangularList = error notTriangularErrorMsg
  | size < 2      = error $ singletonSizeMessage "fromList"
  | otherwise     = fromListUnsafe xs
  where
    len           = length xs
    size          = floor $ sqrt (fromIntegral len :: Double)
    notTriangularList = square size /= len
    square x      = x*x
    notTriangularErrorMsg = unwords
      [ "fromList: The number of elements"
      , parenShow len
      , "is not a square number."
      , "Cannot construct an non-square SDM! "
      , "The number of elements"
      , parenShow len
      , "lies between the valid square numbers"
      , parenShow $ square size
      , "and"
      , parenShow . square $ size + 1
      , ")."
      ]


-- |
-- /O(n*n)/
--
-- Construct a 'SDM' from a list of columns.
--
-- ==== __Examples__
--
-- >>> fromCols [[1,2,3],[4,5,6],[7,8,9]]
-- SDM: 3 x 3
--   1 4 7
--   2 5 8
--   3 6 9
--
fromCols :: (Foldable t, Foldable t', Real a) => t (t' a) -> (Rational, SymbolDistanceMatrixTriangular)
fromCols xs
  | null xs          = error $ nullarySizeMessage "fromCols"
  | hasJaggedCols    = error jaggedColsErrorMsg
  | width /= height  = error notTriangularErrorMsg
  | height < 2       = error $ singletonSizeMessage "fromCols"
  | otherwise        = fromListUnsafe . fold . transpose $ toList <$> toList xs
  where
    width            = length xs
    height           = length . head $ toList xs
    hasJaggedCols    = not . equalityOf length $ toList xs

    jaggedColsErrorMsg = unwords
        [ "fromCols: All the columns did not have the same height!"
        , "Expected modal height of"
        , parenShow mode
        , "but found other heights of"
        , parenShow otherLengths
        ]
      where
        (mode, otherLengths) = modeAndOutlierLengths xs

    notTriangularErrorMsg = unwords
        [ "fromRows: The number of rows"
        , parenShow height
        ,"did not match the number of columns"
        , parenShow width
        , "!"
        ]


-- |
-- /O(n*n)/
--
-- Construct a 'SDM' from a list of rows.
--
-- ==== __Examples__
--
-- >>> fromRows [[1,2,3],[4,5,6],[7,8,9]]
-- SDM: 3 x 3
--   1 2 3
--   4 5 6
--   7 8 9
--
fromRows :: (Foldable t, Foldable t', Real a) => t (t' a) -> (Rational, SymbolDistanceMatrixTriangular)
fromRows xs
  | null xs          = error $ nullarySizeMessage "fromRows"
  | hasJaggedRows    = error jaggedRowsErrorMsg
  | width /= height  = error notTriangularErrorMsg
  | height < 2       = error $ singletonSizeMessage "fromRows"
  | otherwise        = fromListUnsafe . foldMap toList $ toList xs
  where
    height           = length xs
    width            = length . head $  toList xs
    hasJaggedRows    = not $ equalityOf length xs

    jaggedRowsErrorMsg = unwords
        [ "fromRows: All the rows did not have the same width!"
        , "Expected modal width of"
        , parenShow mode
        , "but found other widths of"
        , show otherLengths
        ]
      where
        (mode, otherLengths) = modeAndOutlierLengths xs

    notTriangularErrorMsg = unwords
        [ "fromRows: The number of rows"
        , parenShow height
        ,"did not match the number of columns"
        , parenShow width
        , "!"
        ]


-- |
-- /O(n*n)/
--
-- A generating function for a 'SDM'. Efficiently constructs a 'SDM' of the
-- specified size with each value defined by the result of the supplied
-- function.
--
-- ==== __Examples__
--
-- >>> generate 5 $ const 5
-- SDM: 5 x 5
--   5 5 5 5 5
--   5 5 5 5 5
--   5 5 5 5 5
--   5 5 5 5 5
--   5 5 5 5 5
--
-- >>> generate 4 $ \(i,j) -> abs (i - j)
-- SDM: 4 x 4
--   0 1 2 3
--   1 0 1 2
--   2 1 0 1
--   3 2 1 0
--
-- >>> generate 8 $ \(i,j) -> if i == j || i + j == 6 then 0 else 1
-- SDM: 8 x 8
--   0 1 1 1 1 1 0 1
--   1 0 1 1 1 0 1 1
--   1 1 0 1 0 1 1 1
--   1 1 1 0 1 1 1 1
--   1 1 0 1 0 1 1 1
--   1 0 1 1 1 0 1 1
--   0 1 1 1 1 1 0 1
--   1 1 1 1 1 1 1 0
--
generate :: ( Integral c
            , Integral d
            , Integral i
            )
         => d            -- ^ Number of rows & columns in the SDM.
         -> ((i,i) -> c) -- ^ Function to determine the value of a given index.
         -> SymbolDistanceMatrixTriangular
generate n f
  | n <  0    = error negativeErrorMessage
  | n == 0    = error $ nullarySizeMessage "generate 0 f"
  | n == 1    = error $ singletonSizeMessage "generate 1 f"
  | otherwise = SDMT $ SymbolDistanceMatrix dim resultVector
  where
    dim = SymbolCount $ fromIntegral n
    r   = fromIntegral n :: Int
    resultVector = V.generate (r*r) g
      where

        g :: Int -> Word16
        g i = fromIntegral . f . (fromIntegral *** fromIntegral) $ (i `divMod` r)

    negativeErrorMessage = fold
      [ "The call to 'generate ", show r, " f' is malformed, "
      , "the size ", parenShow r, " is a negative number. "
      , "Cannot construct a SDM with a negative size!"
      ]


-- Un-exported Functionality
--------------------------------------------------------------------------------


-- |
-- Precondition of the list having a square number of elements was already
-- checked in the function wrapping calls to this method.
--
-- This method take a list of values coercible to 'Rational' values via the
-- 'Real' type-class and produces an integral valued SDM with a rational weight.
fromListUnsafe :: (Foldable t, Real a) => t a -> (Rational, SymbolDistanceMatrixTriangular)
fromListUnsafe xs
  | not $ null negativeValues = error $ "The following negative values were found in the SDM: " <> show negativeValues
  | not $ null overflowValues = error $ "The following values are either too small or two large for the SDM's 16-bit precision: " <> show overflowValues
  | otherwise                 = ( 1 % coefficient, SDMT $ SymbolDistanceMatrix dim coercedVector)
  where
    dim               = SymbolCount . floor $ sqrt (fromIntegral (length xs) :: Double)
    coercedVector     = V.fromList $ toEnum . fromEnum <$> prospectiveValues
    negativeValues    = filter (< 0) rationalValues
    overflowValues    = fmap fst . filter (\(_,y) -> y > toRational (maxBound :: Word16)) $ zip rationalValues prospectiveValues
    prospectiveValues = ((coefficient % 1) *) <$> rationalValues
    coefficient       = foldl1 lcm $ abs . denominator <$> rationalValues
    rationalValues    = toRational <$> toList  xs


-- |
-- Determines the mode length and the other lengths of a nested foldable
-- structure.
modeAndOutlierLengths :: (Foldable t, Foldable t') => t (t' a) -> (Int, [Int])
modeAndOutlierLengths xs = (mode, otherLengths)
  where
    occuranceMap = Map.fromList . occurrences $ length <$> toList xs
    (mode,_)     = findMax occuranceMap
    otherLengths = keys  $ mode `delete` occuranceMap


parenShow :: Show a => a -> String
parenShow x = fold ["(", show x, ")"]


nullarySizeMessage :: String -> String
nullarySizeMessage = messagePrefix "An empty structure was supplied. Cannot construct an empty SDM!"

singletonSizeMessage :: String -> String
singletonSizeMessage = messagePrefix "A singleton structure was supplied. Cannot construct a SDM with size of 1, must have size of 2 or greater."


messagePrefix :: String -> String -> String
messagePrefix s x = fold ["Layout.Compact.Symbols.", x, ": ", s]


firstRowExtrema
  :: ( Integral b
     , Num c
     )
  => (Vector Word16 -> b)
  -> SymbolDistanceMatrixTriangular
  -> c
firstRowExtrema f (SDMT (SymbolDistanceMatrix (SymbolCount n) v)) =
    fromIntegral . f $ V.slice 0 (fromEnum n) v


firstColExtrama
  :: ( Integral a
     , Num b
     )
  => ((Int -> Int -> Ordering) -> [Int] -> a)
  -> SymbolDistanceMatrixTriangular
  -> b
firstColExtrama f (SDMT (SymbolDistanceMatrix (SymbolCount n) v)) =
    let r = fromEnum n
        g i = v V.! (i*r)
    in  fromIntegral $ f (comparing g) [ 0 .. r - 1 ]


