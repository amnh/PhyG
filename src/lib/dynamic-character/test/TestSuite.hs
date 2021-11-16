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

module Main
  ( main
  ) where

import qualified DirectOptimization.Pairwise.Test as Pairwise
import           Test.Tasty


-- |
-- Entry point for the test-suite of the "dynamic-character" library.
main :: IO ()
main = defaultMain testSuite


testSuite :: TestTree
testSuite = testGroup "Dynamic Character Test-Suite"
    [ Pairwise.testSuite
    ]
