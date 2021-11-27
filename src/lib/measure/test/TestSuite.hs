------------------------------------------------------------------------------
-- |
-- Module      :  Main
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

module Main where

import Measure.SymbolChangeMatrix.Dense.Test
import Test.Tasty
import Test.Tasty.Ingredients.Rerun          (rerunningTests)


-- |
-- The entry point for the 'Measure.SymbolChangeMatrix.Dense.TCM' test-suite.
main :: IO ()
main =
    defaultMainWithIngredients
    [ rerunningTests defaultIngredients ]
    testSuite
