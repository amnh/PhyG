-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Symbols.Dense
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

module Measure.Symbols.Dense
  ( -- * Dense Representation
    SCMρ()
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

import           Measure.Diagnosis
import           Measure.Symbols.Internal (SCMρ, rowMajorVector, (!))
import qualified Measure.Symbols.Internal as SCM


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
-- *** Exception: fromList: A singleton structure was supplied. Cannot construct a SCM with dimension of 1, must have dimension of 2 or greater.
--
-- >>> fromList [1..12]
-- *** Exception: fromList: The number of element (12) is not a square number. Cannot construct an non-square SCM! The number of elements (12) lies between the valid square numbers (9) and (16).
--
fromList :: (Foldable t, Real a) => t a -> DiagnosisOfSCM SCMρ
fromList = diagnosisBuild SCM.fromList


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
fromCols :: (Foldable t, Foldable t', Real a) => t (t' a) -> DiagnosisOfSCM SCMρ
fromCols = diagnosisBuild SCM.fromCols


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
fromRows :: (Foldable t, Foldable t', Real a) => t (t' a) -> DiagnosisOfSCM SCMρ
fromRows = diagnosisBuild SCM.fromRows


-- |
-- /O(n*n)/
--
-- A generating function for a 'SCM'. Efficiently constructs a 'SCM' of the
-- specified dimensions with each value defined by the result of the supplied
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
         -> DiagnosisOfSCM SCMρ
generate n =
    let g s = (1,s)
    in  diagnosisBuild (g . SCM.generate n)


diagnosisBuild :: (t -> (Rational, SCMρ)) -> t -> DiagnosisOfSCM SCMρ
diagnosisBuild f xs =
    let (factor, scm) = f xs
        diagnosedSCM  = diagnoseSCMρ scm
        coefficient'  = coefficient diagnosedSCM * factor
    in  diagnosedSCM { coefficient = coefficient' }


