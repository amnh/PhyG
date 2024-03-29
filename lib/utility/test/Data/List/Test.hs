------------------------------------------------------------------------------
-- |
-- Module      :  Data.List.Test
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

module Data.List.Test
  ( testSuite
  ) where

import Data.List             (nub, sort)
import Data.List.Utility
import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck as QC


-- |
-- The test-suite for custom "list-like" utilities.
testSuite :: TestTree
testSuite = testGroup "List Tests"
    [ testExampleCases
    , testInvariantProperties
    ]


testExampleCases :: TestTree
testExampleCases = testGroup "Example cases from documentation"
    [ transposeCases
    , isSingletonCases
    , duplicatesCases
    , mostCommonCases
    , occurancesCases
    , chunksOfCases
    , subsetOfCases
    , equalityOfCases
    , invariantTransformationCases
    , transitivePropertyHoldsCases
    , pairwiseSequenceCases
    , maximaByCases
    , minimaByCases
    ]


testInvariantProperties :: TestTree
testInvariantProperties = testGroup "Invariant properties"
    [ transposeProperties
    , isSingletonProperties
    , duplicatesProperties
    , mostCommonProperties
    , occurancesProperties
    , chunksOfProperties
    , subsetOfProperties
    , equalityOfProperties
    , transitivePropertyHoldsProperties
    , pairwiseSequenceProperties
    , maximaByProperties
    , minimaByProperties
    ]


transposeProperties :: TestTree
transposeProperties = testGroup "Properties of transpose"
  [
  ]

isSingletonProperties :: TestTree
isSingletonProperties = testGroup "Properties of isSingleton"
  [ QC.testProperty "Length of singleton list is 1" singletonLength
  ]

  where
    singletonLength :: [()] -> Property
    singletonLength us =
      if isSingleton us
      then length us === 1
      else property True

duplicatesProperties :: TestTree
duplicatesProperties = testGroup "Properties of duplicates"
  [QC.testProperty "duplicates of repeated list contains all original elements" doubleList
  ]

  where
    doubleList :: [Int] -> Property
    doubleList ns
      = (sort . nub . duplicates $ (ns <> ns))
      === (sort . nub $ ns)

mostCommonProperties :: TestTree
mostCommonProperties = testGroup "Properties of mostCommon"
  [ QC.testProperty "mostCommon xs == mostCommon (xs ++ xs)" doubleList
  ]

  where
    doubleList :: [Int] -> Property
    doubleList xs
      = mostCommon xs
      === mostCommon (xs <> xs)

occurancesProperties :: TestTree
occurancesProperties = testGroup "Properties of occurrences"
  [
  ]

chunksOfProperties :: TestTree
chunksOfProperties = testGroup "Properties of chunksOf"
  [
  ]

subsetOfProperties :: TestTree
subsetOfProperties  = testGroup "Properties of subsetOf"
  [
  ]

equalityOfProperties :: TestTree
equalityOfProperties = testGroup "Properties of equalityOf"
  [
  ]

transitivePropertyHoldsProperties :: TestTree
transitivePropertyHoldsProperties = testGroup "Properties of transitiveProperty"
  [
  ]

pairwiseSequenceProperties :: TestTree
pairwiseSequenceProperties = testGroup "Properties of pairwiseSequence"
  [
  ]

maximaByProperties :: TestTree
maximaByProperties = testGroup "Properties of maximaBy"
  [
  ]

minimaByProperties :: TestTree
minimaByProperties = testGroup "Properties of minimaBy"
  [
  ]


-- Unit tests to verify each of the examples in the documentation

transposeCases :: TestTree
transposeCases = testGroup "Cases of transpose"
  [ testCase "transpose [] == [[]]" ex1
  , testCase "transpose [[1]] == [[1]]" ex2
  , testCase "transpose [[1,2], [3,4]] == [[1, 3], [2,4]]" ex3
  , testCase "transpose [[1,2,3],[4,5,6],[7,8,9]] == [[1,4,7],[2,5,8],[3,6,9]]" ex4
  , testCase "transpose [[1,2,3,0,0],[4,5,6,0],[7,8,9]] == [[1,4,7],[2,5,8],[3,6,9]]" ex5
  ]
  where
    ex1, ex2, ex3, ex4, ex5 :: Assertion
    ex1 = transpose []                              @?= ([[]]                      :: [[Int]])
    ex2 = transpose [[1]]                           @?= ([[1]]                     :: [[Int]])
    ex3 = transpose [[1,2], [3,4]]                  @?= ([[1,3], [2,4]]            :: [[Int]])
    ex4 = transpose [[1,2,3],[4,5,6],[7,8,9]]       @?= ([[1,4,7],[2,5,8],[3,6,9]] :: [[Int]])
    ex5 = transpose [[1,2,3,0,0],[4,5,6,0],[7,8,9]] @?= ([[1,4,7],[2,5,8],[3,6,9]] :: [[Int]])


isSingletonCases :: TestTree
isSingletonCases = testGroup "Cases of isSingleton"
  [ testCase "isSingleton [] == False" ex1
  , testCase "isSingleton [()] == True" ex2
  , testCase "isSingleton [(), ()] == False" ex3
  ]

  where
    ex1, ex2, ex3 :: Assertion
    ex1 = isSingleton []       @?= False
    ex2 = isSingleton [()]     @?= True
    ex3 = isSingleton [(), ()] @?= False

duplicatesCases :: TestTree
duplicatesCases = testGroup "Cases of duplicates"
  [ testCase "duplicates \"duplicate string\" == \"it\"" ex1
  , testCase "duplicates \"GATACACATCAGATT\" == \"ACGT\"" ex2
  , testCase "duplicates [ \'A\' .. \'Z\'] == []" ex3
  ]

  where
    ex1, ex2, ex3 :: Assertion
    ex1 = duplicates "duplicate string" @?= "it"
    ex2 = duplicates "GATACACATCAGATT"  @?= "ACGT"
    ex3 = duplicates ['A'..'Z']         @?= []

mostCommonCases :: TestTree
mostCommonCases = testGroup "Cases of mostCommon"
  [ testCase "mostCommon \"GATACACATCAGATT\" == Just 'A'" ex1
  , testCase "mostCommon \"AABCDDDEFGGT\" == Just 'D'" ex2
  ]

  where
    ex1, ex2 :: Assertion
    ex1 = mostCommon "GATACACATCAGATT" @?= Just 'A'
    ex2 = mostCommon "AABCDDDEFGGT"    @?= Just 'D'

occurancesCases :: TestTree
occurancesCases = testGroup "Cases of occurrences"
  [ testCase "occurrences \"GATACACATCAGATT\" == [('A',6),('T',4),('C',3),('G',2)]" ex1
  , testCase
      ( unlines
      [ "occurrences \"AABCDDDEFGGT\""
      , "== [('D',3),('A',2),('G',2),('B',1),('C',1),('E',1),('F',1),('T',1)"
      ]
      )
      ex2
  ]

  where
    ex1, ex2 :: Assertion
    ex1 = occurrences "GATACACATCAGATT"
      @?= [('A',6),('T',4),('C',3),('G',2)]
    ex2 = occurrences "AABCDDDEFGGT"
      @?= [('D',3),('A',2),('G',2),('B',1),('C',1),('E',1),('F',1),('T',1)]

chunksOfCases :: TestTree
chunksOfCases = testGroup "Cases of chunksOf"
  [ testCase "chunksOf 3 [1..13] == [[1,2,3],[4,5,6],[7,8,9],[10,11,12],[13]]" ex1
  , testCase "chunksOf 5 [1..13] == [[1,2,3,4,5],[6,7,8,9,10],[11,12,13]]" ex2
  ]

  where
    ex1, ex2 :: Assertion
    ex1 = chunksOf 3 [1..13 :: Int] @?= [[1,2,3],[4,5,6],[7,8,9],[10,11,12],[13]]
    ex2 = chunksOf 5 [1..13 :: Int] @?= [[1,2,3,4,5],[6,7,8,9,10],[11,12,13]]

subsetOfCases :: TestTree
subsetOfCases  = testGroup "Cases of subsetsOf"
  [ testCase "([5..10] `subsetOf` [1..13]) == True" ex1
  , testCase "([11..15] `subsetOf` [1..13] == False" ex2
  ]

  where
    ex1, ex2 :: Assertion
    ex1 = [5..10]  `subsetOf` [1..13 :: Int] @?= True
    ex2 = [11..15] `subsetOf` [1..13 :: Int] @?= False

equalityOfCases :: TestTree
equalityOfCases = testGroup "Cases of equalityOf"
  [ testCase "equalityOf (`mod` 10) [9,19,29,39,49] == True" ex1
  , testCase "equalityOf (`mod` 7) [9,19,29,39,49] == False" ex2
  ]

  where
    ex1, ex2 :: Assertion
    ex1 = equalityOf (`mod` 10) [9,19..49 :: Int] @?= True
    ex2 = equalityOf (`mod` 7)  [9,19..49 :: Int] @?= False

invariantTransformationCases :: TestTree
invariantTransformationCases = testGroup "Cases of invariantTransformation"
  [ testCase "invariantTransformation (`mod` 10) [9,19,29,39,49] == Just 9" ex1
  , testCase "invariantTransformation (`mod`  7) [9,19,29,39,49] == Nothing" ex2
  ]

  where
    ex1, ex2 :: Assertion
    ex1 = invariantTransformation (`mod` 10) [9,19..49 :: Int]  @?= Just 9
    ex2 = invariantTransformation (`mod` 7 ) [9,19..49 :: Int]  @?= Nothing

transitivePropertyHoldsCases :: TestTree
transitivePropertyHoldsCases = testGroup "cases of transitivePropertyHolds"
  [ testCase
      ("transitivePropertyHolds (\\ x y -> snd x >= fst y)"
      <> "[(9,9), (8,7), (6,6), (6,5), (3,4), (3,0) ]"
      )
      ex1
  ]

  where
    ex1 :: Assertion
    ex1 =
       transitivePropertyHolds (\x y -> snd x >= fst y)
         [(9,9), (8,7), (6,6), (6,5), (3,4), (3,0) :: (Int,Int)]
       @?= True

pairwiseSequenceCases :: TestTree
pairwiseSequenceCases = testGroup "Cases of pairwiseSequence"
  [ testCase
      ( unlines
      [ "pairwiseSequence (\\ x y -> snd x /= snd y)"
      , "[[('A',1),('B',2)],[('X',1),('Y',2),('Z',3)],[('I',1),('J',2),('K',3),('L',4)]]"
      , "=="
      , "[ [('A',1),('Y',2),('K',3)]"
      , ", [('A',1),('Y',2),('L',4)]"
      , ", [('A',1),('Z',3),('J',2)]"
      , ", [('A',1),('Z',3),('L',4)]"
      , ", [('B',2),('X',1),('K',3)]"
      , ", [('B',2),('X',1),('L',4)]"
      , ", [('B',2),('Z',3),('I',1)]"
      , ", [('B',2),('Z',3),('L',4)]"
      , "]"
      ]
      )
      ex1
  ]

  where
    ex1 :: Assertion
    ex1 =
      pairwiseSequence
      (\x y -> snd x /= snd y)
      [ [('A',1),('B',2)]
      , [('X',1),('Y',2),('Z',3)]
      , [('I',1),('J',2),('K',3),('L',4) :: (Char,Int)]
      ]
      @?=
      [ [('A',1),('Y',2),('K',3)]
      , [('A',1),('Y',2),('L',4)]
      , [('A',1),('Z',3),('J',2)]
      , [('A',1),('Z',3),('L',4)]
      , [('B',2),('X',1),('K',3)]
      , [('B',2),('X',1),('L',4)]
      , [('B',2),('Z',3),('I',1)]
      , [('B',2),('Z',3),('L',4)]
      ]

maximaByCases :: TestTree
maximaByCases = testGroup "Cases of maximaBy"
  [
  ]

minimaByCases :: TestTree
minimaByCases = testGroup "Cases of minimaBy"
  [
  ]
