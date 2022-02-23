-----------------------------------------------------------------------------
-- |
-- Module      :  Layout.Compact.Symbols
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

module Layout.Compact.Symbols
  ( -- * Compact Representation
    SymbolDistanceMatrixSquare()
    -- * Constructors
  , fromList
  , fromCols
  , fromRows
  , generate
    -- * Denstructors
  , rowMajorVector
    -- * Indexing
  , (!)
  ) where

import           Layout.Compact.Symbols.Square (SymbolDistanceMatrixSquare, rowMajorVector, (!))
import qualified Layout.Compact.Symbols.Square as SDM
import           Measure.Transition.Diagnosis


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
-- *** Exception: fromList: A singleton structure was supplied. Cannot construct a SDM with dimension of 1, must have dimension of 2 or greater.
--
-- >>> fromList [1..12]
-- *** Exception: fromList: The number of element (12) is not a square number. Cannot construct an non-square SDM! The number of elements (12) lies between the valid square numbers (9) and (16).
--
fromList :: (Foldable t, Real a) => t a -> DiagnosisOfSymbolDistance
fromList = diagnosisBuild SDM.fromList


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
fromCols :: (Foldable t, Foldable t', Real a) => t (t' a) -> DiagnosisOfSymbolDistance
fromCols = diagnosisBuild SDM.fromCols


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
fromRows :: (Foldable t, Foldable t', Real a) => t (t' a) -> DiagnosisOfSymbolDistance
fromRows = diagnosisBuild SDM.fromRows


-- |
-- /O(n*n)/
--
-- A generating function for a 'SDM'. Efficiently constructs a 'SDM' of the
-- specified dimensions with each value defined by the result of the supplied
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
         -> DiagnosisOfSymbolDistance
generate n =
    let g s = (1,s)
    in  diagnosisBuild (g . SDM.generate n)


diagnosisBuild :: (t -> (Rational, SymbolDistanceMatrixSquare)) -> t -> DiagnosisOfSymbolDistance
diagnosisBuild f xs =
    let (factor, sdm) = f xs
        diagnosedSDM  = diagnoseSymbolDistanceMatrixSquare sdm
        coefficient'  = factoredCoefficient diagnosedSDM * factor
    in  diagnosedSDM { factoredCoefficient = coefficient' }


