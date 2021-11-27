-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Matrix.Types
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

module Measure.Matrix.Types
  ( -- * Different 
    SCMλ
  , TCM2Dλ
  , TCM3Dλ
  ) where

import Measure.Unit.Distance
import Measure.Unit.SymbolIndex


-- |
-- An abstract representation of the /distance/ bewteen two symbols in an alphabet.
-- Given the indicies of two symbols, the distance between the symbols is returned.
type SCMλ = SymbolIndex -> SymbolIndex -> Distance


-- |
-- An abstract representation of the /distance/ bewteen /two/ states.
-- Given two states, the distance and median between states is returned.
--
-- Abstractly, this function is the cross-product between the collection of all
-- possible states, forming a matrix of each pair's diatance and median. This is
-- the Transition Cost Matrix (TCM) for the collection of possible states.
--
-- Naturally, matrices are 2-dimensional. However, there exists a similar three-way
-- distance and median calculation. Rather than call this a Transition Cost Cube
-- (TCC), we instead call these 'TCM2Dλ' and 'TCM3Dλ'.
type TCM2Dλ e = e -> e -> (Distance, e)


-- |
-- An abstract representation of the /distance/ bewteen /three/ states.
-- Given three states, the distance and median between states is returned.
--
-- A higher dimensional version of 'TCM2Dλ'.
type TCM3Dλ e = e -> e -> e -> (Distance, e)
