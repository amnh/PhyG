{-# LANGUAGE ImportQualifiedPost #-}

{- |
Test-suite for the dynamic characters.
-}
module Main (
    main,
) where

import DirectOptimization.Pairwise.Test qualified as Pairwise
import Test.Tasty


{- |
Entry point for the test-suite of the "dynamic-character" library.
-}
main ∷ IO ()
main = defaultMain testSuite


testSuite ∷ TestTree
testSuite =
    testGroup
        "Dynamic Character Test-Suite"
        [ Pairwise.testSuite
        ]
