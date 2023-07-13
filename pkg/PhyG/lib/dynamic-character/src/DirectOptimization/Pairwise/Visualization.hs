-----------------------------------------------------------------------------
-- |
-- Module      :  DirectOptimization.Pairwise.Visualization
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-- Direct optimization pairwise alignment using the Needleman-Wunsch algorithm.
-- These functions will allocate an M * N matrix.
--
-----------------------------------------------------------------------------

{-# LANGUAGE ApplicativeDo         #-}
{-# LANGUAGE ConstraintKinds       #-}
{-# LANGUAGE DerivingStrategies    #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE Strict                #-}
{-# LANGUAGE TypeFamilies          #-}
{-# LANGUAGE UnboxedTuples         #-}

{-# OPTIONS_GHC -Wno-incomplete-patterns #-}
{-# OPTIONS_GHC -Wno-incomplete-uni-patterns #-}

module DirectOptimization.Pairwise.Visualization
  ( Direction()
    -- * Operational
  , directOptimizationDiffDirectionMatricies
    -- * Rendering
  , renderMatrix
  , renderDirectionMatrix
  , diffDirectionMatrix
  ) where

import Bio.DynamicCharacter
import Bio.DynamicCharacter.Measure
import Data.Bits
import Data.Foldable                         (fold)
import Data.Matrix.Class                     (Matrix, dim, toLists, unsafeIndex)
import Data.Set                              (Set, fromDistinctAscList, member)
import Data.Vector.Generic                   (Vector, basicLength, toList)
import DirectOptimization.Pairwise.Direction


directOptimizationDiffDirectionMatricies
  :: ( FiniteBits e
     , Matrix m t Direction
     , Ord (v e)
     , Vector v e
     )
  => (v e -> v e -> (Word, m t Direction))
  -> (v e -> v e -> (Word, m t Direction))
  -> OpenDynamicCharacter v e
  -> OpenDynamicCharacter v e
  -> String
directOptimizationDiffDirectionMatricies matrixGenerator1 matrixGenerator2 lhs rhs =
    let -- Remove gaps from the inputs and measure the results to determine
        -- which ungapped character is longer and which is shorter.
        -- Always pass the shorter character into alignment functions first!
        ~(_, _, _, lesser, longer) = measureCharactersWithoutGaps lhs rhs
        lesserMeds = extractMediansGapped lesser
        longerMeds = extractMediansGapped longer
    in  case basicLength lesserMeds of
          -- Neither character was Missing, but one or both are empty when gaps are removed
          0 ->  "One character was all gaps"
          -- Both have some non-gap elements, perform string alignment
          _ ->  let dm1 = snd $ matrixGenerator1 lesserMeds longerMeds
                    dm2 = snd $ matrixGenerator2 lesserMeds longerMeds
                in  diffDirectionMatrix lesserMeds longerMeds dm1 dm2


-- |
-- Serializes an alignment matrix to a 'String'. Uses input characters for row
-- and column labelings.
--
-- Useful for debugging purposes.
renderMatrix
  :: ( Matrix m t x
     , Vector v e
     , Show x
     )
  => v e   -- ^ Shorter vector of elements
  -> v e   -- ^ Longer  vector of elements
  -> m t x -- ^ Matrix of cells
  -> String
renderMatrix lesser longer mtx = unlines
    [ dimensionPrefix
    , headerRow
    , barRow
    , renderedRows
    ]
  where
    (colCount, rowCount, lesserTokens, longerTokens, maxPrefixWidth, maxColumnWidth, [matrixTokens]) =
        getMatrixConstants lesser longer [mtx]
{-
    toShownIntegers   = fmap (show . showBitsValue) . otoList

    showBitsValue :: FiniteBits b => b -> Word
    showBitsValue b = go (finiteBitSize b) 0
      where
        go 0 v = v
        go i v = let i' = i-1
                     v' | b `testBit` i' = v + bit i'
                        | otherwise      = v
                 in  go i' v'
-}

    dimensionPrefix  = " " <> unwords
        [ "Dimensions:"
        , show rowCount
        , "X"
        , show colCount
        ]

    headerRow = fold
        [ " "
        , pad maxPrefixWidth "\\"
        , "| "
        , pad maxColumnWidth "*"
        , concatMap (pad maxColumnWidth) longerTokens
        ]

    barRow = fold
        [ " "
        , bar maxPrefixWidth
        , "+"
        , concatMap (const (bar maxColumnWidth)) $ undefined : longerTokens
        ]
      where
        bar n = replicate (n+1) '-'

    renderedRows = unlines $ zipWith renderRow ("*":lesserTokens) matrixTokens
      where
        renderRow e cs = prefix <> suffix
          where
            prefix = fold [" ", pad maxPrefixWidth e, "| "]
            suffix = concatMap (pad maxColumnWidth) cs

    pad :: Int -> String -> String
    pad n e = replicate (n - len) ' ' <> e <> " "
      where
        len = length e


-- |
-- Serializes an alignment matrix to a 'String'. Uses input characters for row
-- and column labelings.
--
-- Useful for debugging purposes.
renderDirectionMatrix
  :: ( Matrix m t Direction
     , Vector v e
     )
  => v e
  -> v e
  -> m t Direction
  -> String
renderDirectionMatrix lesser longer mtx = unlines
    [ dimensionPrefix
    , headerRow
    , barRow
    , renderedRows
    ]
  where
    (colCount, rowCount, lesserTokens, longerTokens, maxPrefixWidth, maxColumnWidth, [matrixTokens]) =
        getMatrixConstants lesser longer [mtx]

    tracebackCells = getTracebackIndices mtx

    dimensionPrefix  = " " <> unwords
        [ "Dimensions:"
        , show rowCount
        , "X"
        , show colCount
        ]

    headerRow = fold
        [ " "
        , pad maxPrefixWidth "\\"
        , "| "
        , pad maxColumnWidth "*"
        , concatMap (pad maxColumnWidth) longerTokens
        ]

    barRow = fold
        [ " "
        , bar maxPrefixWidth
        , "+"
        , concatMap (const (bar maxColumnWidth)) $ undefined : longerTokens
        ]
      where
        bar n = replicate (n+1) '-'

    renderedRows = unlines $ zipWith3 renderRow [0..] ("*":lesserTokens) matrixTokens
      where
        renderRow i e cs = prefix <> suffix
          where
            prefix = fold [" ", pad maxPrefixWidth e, "| "]
            suffix = fold $ zipWith (renderCell i) [0..] cs

        renderCell i j = fmap f . pad maxColumnWidth
          where
            f | (i,j) `elem` tracebackCells = boldDirection
              | otherwise = id

    pad :: Int -> String -> String
    pad n e = replicate (n - len) ' ' <> e <> " "
      where
        len = length e


-- |
-- Serializes an alignment matrix to a 'String'. Uses input characters for row
-- and column labelings.
--
-- Useful for debugging purposes.
diffDirectionMatrix
  :: ( Matrix m t Direction
     , Vector v e
     )
  => v e
  -> v e
  -> m t Direction
  -> m t Direction
  -> String
diffDirectionMatrix lesser longer mtx1 mtx2 = unlines
    [ dimensionPrefix
    , headerRow
    , barRow
    , renderedRows
    ]
  where
    (colCount, rowCount, lesserTokens, longerTokens, maxPrefixWidth, maxColumnWidth, [tokMtx1,tokMtx2]) =
        getMatrixConstants lesser longer [mtx1, mtx2]

    tracebackCells1 = getTracebackIndices mtx1
    tracebackCells2 = getTracebackIndices mtx2

    dimensionPrefix  = " " <> unwords
        [ "Dimensions:"
        , show rowCount
        , "X"
        , show colCount
        ]

    headerRow = fold
        [ " "
        , pad maxPrefixWidth "\\"
        , "| "
        , pad maxColumnWidth "*"
        , concatMap (pad maxColumnWidth) longerTokens
        ]

    barRow = fold
        [ " "
        , bar maxPrefixWidth
        , "+"
        , concatMap (const (bar maxColumnWidth)) $ undefined : longerTokens
        ]
      where
        bar n = replicate (n+1) '-'

    renderedRows = unlines $ zipWith3 renderRow ("*":lesserTokens) strMtx1 strMtx2
      where
        renderRow e xs ys = prefix <> suffix
          where
            prefix = fold [" ", pad maxPrefixWidth e, "| "]
            suffix = fold $ zipWith renderCell xs ys

        renderCell x y
          | x' == y'  = x
          | otherwise = replicate (max (length x) (length y)) ' '
          where
            x' = boldDirection <$> x
            y' = boldDirection <$> y

    strMtx1 = tok2Str tracebackCells1 tokMtx1
    strMtx2 = tok2Str tracebackCells2 tokMtx2

    tok2Str s = zipWith f [0..]
      where
        f i   = zipWith (g i) [0..]
        g i j = fmap h . pad maxColumnWidth
          where
            h | (i,j) `member` s = boldDirection
              | otherwise = id


    pad :: Int -> String -> String
    pad n e = replicate (n - len) ' ' <> e <> " "
      where
        len = length e


-- |
-- Get the indices of the traceback route.
getTracebackIndices
  :: Matrix m t Direction
  => m t Direction
  -> Set (Int, Int)
getTracebackIndices mtx = fromDistinctAscList $ go (# m - 1, n - 1 #)
  where
    getDirection = curry $ unsafeIndex mtx
    (m,n) = dim mtx
    go (# i, j #)
      | i < 0 || j < 0 = []
      | (i,j) == (0,0) = [(0,0)]
      | otherwise =
          (i,j) : case getDirection i j of
                    LeftArrow -> go (# i    , j - 1 #)
                    DiagArrow -> go (# i - 1, j - 1 #)
                    UpArrow   -> go (# i - 1, j     #)


characterVectorToIndices :: Vector v e => v e -> [String]
characterVectorToIndices =
    let numbers = tail $ pure <$> cycle ['0'..'9']
    in  zipWith const numbers . toList


tokenizeMatrix :: (Matrix m t x, Show x) => m t x -> [[String]]
tokenizeMatrix = fmap (fmap show) . toLists


maxLengthOfGrid :: (Foldable g, Foldable r, Foldable f, Functor g, Functor r) => g (r (f a)) -> Int
maxLengthOfGrid = maximum . fmap maxLengthOfRow


maxLengthOfRow :: (Foldable r, Foldable f, Functor r) => r (f a) -> Int
maxLengthOfRow = maximum . fmap length


getMatrixConstants
  :: ( Matrix m t x
     , Show x
     , Vector v e
     )
  => v e
  -> v e
  -> [m t x]
  -> (Int, Int, [String], [String], Int, Int, [[[String]]])
getMatrixConstants lesser longer matrices =
    (colCount, rowCount, lesserTokens, longerTokens, maxPrefixWidth, maxColumnWidth, matrixTokens)
  where
    colCount       = basicLength longer + 1
    rowCount       = basicLength lesser + 1
    lesserTokens   = characterVectorToIndices lesser
    longerTokens   = characterVectorToIndices longer
    maxPrefixWidth = maxLengthOfRow lesserTokens
    maxHeaderWidth = maxLengthOfRow longerTokens

    matrixTokens   = tokenizeMatrix <$> matrices
    maxColumnWidth = maximum $
        maxHeaderWidth : (maxLengthOfGrid <$> matrixTokens)


