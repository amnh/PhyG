-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Symbols.Internal
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE DeriveAnyClass     #-}
{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE DeriveGeneric      #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE Strict             #-}
{-# LANGUAGE TypeFamilies       #-}

module Measure.Symbols.Internal
  ( -- * Representational Type
    SCMρ(..)
    -- * Queries
  , index
  , (!)
  , sizeOfSCM
  , rowMajorVector
    -- * Constructors
  , fromList
  , fromCols
  , fromRows
  , generate
  ) where

import           Control.Arrow        ((***))
import           Control.DeepSeq
import           Data.Data
import           Data.Foldable
import           Data.List            (transpose)
import           Data.List.Utility    (equalityOf, occurrences)
import           Data.Map             (delete, findMax, keys)
import qualified Data.Map             as Map (fromList)
import           Data.Ord
import           Data.Ratio
import           Data.Vector.Storable (Vector)
import qualified Data.Vector.Storable as V
import           Data.Word
import           GHC.Generics
import           Measure.Unit
import           Measure.Matrix


-- |
-- A data structure for storing a square array of sizeality
-- greater than or equal to two, with positive cost values at the array indices.
-- Values are stored in an unboxed structure for cache efficiency.
--
-- A 'SCM' can be constructed by calling one of the following functions:
--
--   * 'fromList'
--
--   * 'fromCols'
--
--   * 'fromRows'
--
--   * 'generate'
--
-- Attempts to construct an empty or singleton 'SCM' through the above
-- constructors will result in a runtime exception.
data  SCMρ
    = SCMρ {-# UNPACK #-} !SymbolCount {-# UNPACK #-} !(Vector Word16)
    deriving stock    (Data, Eq, Generic, Typeable)
    deriving anyclass (NFData)


-- |
-- Any structural representation which can produce a Symbol Change Matrix.
instance HasEditCosts SCMρ where

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
instance HasSymbolChangeMatrix SCMρ where

    {-# INLINE getSCMλ #-}
    getSCMλ = index


-- |
-- A structure which can derive the number of alphabet symbols associated with it.
instance HasSymbolCount SCMρ where

    {-# INLINE symbolCount #-}
    symbolCount = sizeOfSCM


-- |
-- A pretty printed custom show instance for 'SCM'.
instance Show SCMρ where

    show scm = headerLine <> matrixLines
      where
        renderRow i = ("  "<>) . unwords $ renderValue <$> [ scm ! (i,j) | j <- rangeValues ]
        matrixLines = unlines $ renderRow   <$> rangeValues
        rangeValues = [0 .. symbolCount scm - 1]
        headerLine  = '\n' : unwords [ "SCM:", show $ symbolCount scm, "x", show $ symbolCount scm, "\n"]
        maxValue    = V.maximum $ rowMajorVector scm
        padSpacing  = length $ show maxValue
        renderValue x = pad <> shown
          where
            shown = show (x :: Word)
            pad   = (padSpacing - length shown) `replicate` ' '


-- |
-- /O(1)/
--
-- Indexing without bounds checking.
{-# INLINE index #-}
index :: (Integral c, Integral i) => SCMρ -> i -> i -> c
index (SCMρ n v) i j = fromIntegral $ v `V.unsafeIndex` (fromIntegral i * fromIntegral n + fromIntegral j)


-- |
-- /O(1)/
--
-- Indexing without bounds checking.
{-# INLINE (!) #-}
(!) :: (Integral c, Integral i) => SCMρ -> (i,i) -> c
(!) scm = uncurry $ index scm


-- |
-- /O(1)/
--
-- Size of the symbol change matrix.
sizeOfSCM :: SCMρ -> SymbolCount
sizeOfSCM (SCMρ n _) = n


-- |
-- Deconstructs the 'SCMρ' to expose the underlying unboxed 'Vector'.
{-# INLINE rowMajorVector #-}
rowMajorVector :: SCMρ -> Vector Word16
rowMajorVector (SCMρ _ v) = v


-- |
-- /O(n*n)/
--
-- Construct a 'SCM' from a list of elements in row major order.
--
-- ==== __Examples__
--
-- >>> fromList [1..9]
-- SCM: 3 x 3
--   1 2 3
--   4 5 6
--   7 8 9
--
-- >>> fromList []
-- *** Exception: fromList: An empty structure was supplied. Cannot construct an empty SCM!
--
-- >>> fromList [42]
-- *** Exception: fromList: A singleton structure was supplied. Cannot construct a SCM with size of 1, must have size of 2 or greater.
--
-- >>> fromList [1..12]
-- *** Exception: fromList: The number of element (12) is not a square number. Cannot construct an non-square SCM! The number of elements (12) lies between the valid square numbers (9) and (16).
--
fromList :: (Foldable t, Real a) => t a -> (Rational, SCMρ)
fromList xs
  | null xs       = error $ nullarySizeMessage "fromList"
  | notSquareList = error notSquareErrorMsg
  | size < 2      = error $ singletonSizeMessage "fromList"
  | otherwise     = fromListUnsafe xs
  where
    len           = length xs
    size          = floor $ sqrt (fromIntegral len :: Double)
    notSquareList = square size /= len
    square x      = x*x
    notSquareErrorMsg = unwords
      [ "fromList: The number of elements"
      , parenShow len
      , "is not a square number."
      , "Cannot construct an non-square SCM! "
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
-- Construct a 'SCM' from a list of columns.
--
-- ==== __Examples__
--
-- >>> fromCols [[1,2,3],[4,5,6],[7,8,9]]
-- SCM: 3 x 3
--   1 4 7
--   2 5 8
--   3 6 9
--
fromCols :: (Foldable t, Foldable t', Real a) => t (t' a) -> (Rational, SCMρ)
fromCols xs
  | null xs          = error $ nullarySizeMessage "fromCols"
  | hasJaggedCols    = error jaggedColsErrorMsg
  | width /= height  = error notSquareErrorMsg
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

    notSquareErrorMsg = unwords
        [ "fromRows: The number of rows"
        , parenShow height
        ,"did not match the number of columns"
        , parenShow width
        , "!"
        ]


-- |
-- /O(n*n)/
--
-- Construct a 'SCM' from a list of rows.
--
-- ==== __Examples__
--
-- >>> fromRows [[1,2,3],[4,5,6],[7,8,9]]
-- SCM: 3 x 3
--   1 2 3
--   4 5 6
--   7 8 9
--
fromRows :: (Foldable t, Foldable t', Real a) => t (t' a) -> (Rational, SCMρ)
fromRows xs
  | null xs          = error $ nullarySizeMessage "fromRows"
  | hasJaggedRows    = error jaggedRowsErrorMsg
  | width /= height  = error notSquareErrorMsg
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

    notSquareErrorMsg = unwords
        [ "fromRows: The number of rows"
        , parenShow height
        ,"did not match the number of columns"
        , parenShow width
        , "!"
        ]


-- |
-- /O(n*n)/
--
-- A generating function for a 'SCM'. Efficiently constructs a 'SCM' of the
-- specified size with each value defined by the result of the supplied
-- function.
--
-- ==== __Examples__
--
-- >>> generate 5 $ const 5
-- SCM: 5 x 5
--   5 5 5 5 5
--   5 5 5 5 5
--   5 5 5 5 5
--   5 5 5 5 5
--   5 5 5 5 5
--
-- >>> generate 4 $ \(i,j) -> abs (i - j)
-- SCM: 4 x 4
--   0 1 2 3
--   1 0 1 2
--   2 1 0 1
--   3 2 1 0
--
-- >>> generate 8 $ \(i,j) -> if i == j || i + j == 6 then 0 else 1
-- SCM: 8 x 8
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
         => d            -- ^ Number of rows & columns in the SCM.
         -> ((i,i) -> c) -- ^ Function to determine the value of a given index.
         -> SCMρ
generate n f
  | n <  0    = error negativeErrorMessage
  | n == 0    = error $ nullarySizeMessage "generate 0 f"
  | n == 1    = error $ singletonSizeMessage "generate 1 f"
  | otherwise = SCMρ (fromIntegral n) resultVector
  where
    r = fromIntegral n :: Int
    resultVector = V.generate (r*r) g
      where

        g :: Int -> Word16
        g i = fromIntegral . f . (fromIntegral *** fromIntegral) $ (i `divMod` r)

    negativeErrorMessage = fold
      [ "The call to 'generate ", show r, " f' is malformed, "
      , "the size ", parenShow r, " is a negative number. "
      , "Cannot construct a SCM with a negative size!"
      ]


-- Un-exported Functionality
--------------------------------------------------------------------------------


-- |
-- Precondition of the list having a square number of elements was already
-- checked in the function wrapping calls to this method.
--
-- This method take a list of values coercible to 'Rational' values via the
-- 'Real' type-class and produces an integral valued SCM with a rational weight.
fromListUnsafe :: (Foldable t, Real a) => t a -> (Rational, SCMρ)
fromListUnsafe xs
  | not $ null negativeValues = error $ "The following negative values were found in the SCM: " <> show negativeValues
  | not $ null overflowValues = error $ "The following values are either too small or two large for the SCM's 16-bit precision: " <> show overflowValues
  | otherwise                 = ( 1 % coefficient, SCMρ dim coercedVector)
  where
    dim               = floor $ sqrt (fromIntegral (length xs) :: Double)
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
nullarySizeMessage = messagePrefix "An empty structure was supplied. Cannot construct an empty SCM!"

singletonSizeMessage :: String -> String
singletonSizeMessage = messagePrefix "A singleton structure was supplied. Cannot construct a SCM with size of 1, must have size of 2 or greater."


messagePrefix :: String -> String -> String
messagePrefix s x = fold ["Measure.Symbols.", x, ": ", s]


firstRowExtrema :: (Integral b, Num c) => (Vector Word16 -> b) -> SCMρ -> c
firstRowExtrema f (SCMρ n v) = fromIntegral . f $ V.slice 0 (fromEnum n) v


firstColExtrama :: (Integral a, Num b) => ((Int -> Int -> Ordering) -> [Int] -> a) -> SCMρ -> b
firstColExtrama f (SCMρ n v) =
    let r = fromEnum n
        g i = v V.! (i*r)
    in  fromIntegral $ f (comparing g) [ 0 .. r - 1 ]


