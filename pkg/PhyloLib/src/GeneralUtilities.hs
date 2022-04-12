{- |
Module      :  GeneralUtilities.hs
Description :  Module with useful functions
Copyright   :  (c) 2021 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
License     :

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

Maintainer  :  Ward Wheeler <wheeler@amnh.org>
Stability   :  unstable
Portability :  portable (I hope)

-}

{-# Language BangPatterns #-}
{-# Language ImportQualifiedPost #-}
{-# Language ScopedTypeVariables #-}

module GeneralUtilities
    ( module GeneralUtilities
    ) where

import Data.Array
import Data.Text qualified  as T
import Data.Text.Lazy qualified as TL
import System.Random
import Data.Array.IO
import Data.Foldable
import Control.Monad
import Control.DeepSeq
import Data.Time
import Data.Time.Clock.POSIX
import System.IO.Unsafe
import Text.Read
import Data.Maybe
import Data.Bits
import Data.BitVector.LittleEndian qualified as BV
import Data.List qualified as L


-- | functions for triples, quadruples
fst3 :: (a,b,c) -> a
fst3 (d,_,_) = d

snd3 :: (a,b,c) -> b
snd3 (_,e,_) = e

thd3 :: (a,b,c) -> c
thd3 (_,_,e) = e

fst4 :: (a,b,c,d) -> a
fst4 (e,_,_,_) = e

snd4 :: (a,b,c,d) -> b
snd4 (_,e,_,_) = e

thd4 :: (a,b,c,d) -> c
thd4 (_,_,e,_) = e

fth4 :: (a,b,c,d) -> d
fth4 (_,_,_,f) = f

fst5 :: (a,b,c,d,e) -> a
fst5 (e,_,_,_,_) = e

snd5 :: (a,b,c,d,e) -> b
snd5 (_,e,_,_,_) = e

thd5 :: (a,b,c,d,e)-> c
thd5 (_,_,e,_,_) = e

fth5 :: (a,b,c,d,e) -> d
fth5 (_,_,_,e,_) = e

fft5 :: (a,b,c,d,e) -> e
fft5 (_,_,_,_,e) = e

fst6 :: (a,b,c,d,e,f) -> a
fst6 (e,_,_,_,_,_) = e

snd6 :: (a,b,c,d,e,f) -> b
snd6 (_,e,_,_,_,_) = e

thd6 :: (a,b,c,d,e,f)-> c
thd6 (_,_,e,_,_,_) = e

fth6 :: (a,b,c,d,e,f) -> d
fth6 (_,_,_,e,_,_) = e

fft6 :: (a,b,c,d,e,f) -> e
fft6 (_,_,_,_,e,_) = e

six6 :: (a,b,c,d,e,f) -> f
six6 (_,_,_,_,_,e) = e

-- | doubleAsInt takes floor and ceil of Double and retuns Maybe Int
-- nothing if not, Just Int if it is
doubleAsInt :: Double -> Maybe Int
doubleAsInt inDouble =
    if ceiling inDouble /= floor inDouble then Nothing
    else Just (floor inDouble :: Int)

-- | editDistance is a naive edit distance between two lists
-- takes two  lists and returns edit distance
--- from  https://wiki.haskell.org/Edit_distance
editDistance :: Eq a => [a] -> [a] -> Int
editDistance xs ys = table ! (m,n)
    where
    (m,n) = (length xs, length ys)
    x     = array (1,m) (zip [1..] xs)
    y     = array (1,n) (zip [1..] ys)

    table :: Array (Int,Int) Int
    table = array bnds [(ij, dist ij) | ij <- range bnds]
    bnds  = ((0,0),(m,n))

    dist (0,j) = j
    dist (i,0) = i
    dist (i,j) = minimum [table ! (i-1,j) + 1, table ! (i,j-1) + 1,
        if x ! i == y ! j then table ! (i-1,j-1) else 1 + table ! (i-1,j-1)]

-- | checkCommandArgs takes comamnd and args and verifies that they are in list
checkCommandArgs :: String -> [String] -> [String] -> Bool
checkCommandArgs commandString commandList permittedList =
    null commandList || (
    let firstCommand = head commandList
        foundCommand = firstCommand `elem` permittedList
    in
    if foundCommand then checkCommandArgs commandString (tail commandList) permittedList
    else
        let errorMatch = snd $ getBestMatch (maxBound :: Int ,"no suggestion") permittedList firstCommand
        in  errorWithoutStackTrace $ fold
              [ "\nError: Unrecognized '"
              , commandString
              , "' option. By '"
              , firstCommand
              , "' did you mean '"
              , errorMatch
              , "'?\n"
              ] )


-- | getBestMatch compares input to allowable commands and checks if in list and if not outputs
-- closest match
-- call with (maxBound :: Int ,"no suggestion") commandList inString
getBestMatch :: (Int, String) -> [String] -> String -> (Int, String)
getBestMatch curBest@(minDist, _) allowedStrings inString =
    if null allowedStrings then curBest
    else
        let candidate =  head allowedStrings
            candidateEditCost = editDistance candidate inString
        in
        if candidateEditCost == 0 then (0, candidate)
        else if candidateEditCost < minDist then getBestMatch (candidateEditCost, candidate) (tail allowedStrings) inString
        else getBestMatch curBest (tail allowedStrings) inString

-- | getCommandErrorString takes list of non zero edits to allowed commands and reurns meaningful error string
getCommandErrorString :: [(Int, String, String)] -> String
getCommandErrorString noMatchList =
    if null noMatchList then ""
    else
        let (_, firstCommand, firstMatch) = head noMatchList
            firstError = "\tBy \'" ++ firstCommand ++ "\' did you mean \'" ++ firstMatch ++ "\'?\n"
        in
        firstError ++ getCommandErrorString (tail noMatchList)

-- | isSequentialSubsequence takes two lists and determines if the first List is
-- a subsequence of the second but the elements must be sequencetial unlike
-- isSubsequenceOf in Data.List
-- Uses Text.filter to see if there is a match
--isSequentialSubsequence :: (Eq a) => [a] -> [a] -> Bool
isSequentialSubsequence :: String -> String -> Bool
isSequentialSubsequence firstL secondL
  | null firstL = False
  | length firstL > length secondL = False
  | otherwise =
    let foundNumber = T.count  (T.pack firstL) (T.pack secondL)
    in
    foundNumber /= 0

-- | shuffle Randomly shuffles a list
--   /O(N)/
-- from https://wiki.haskell.org/Random_shuffle
shuffleIO :: [a] -> IO [a]
shuffleIO xs = do
        ar <- newArrayLocal  n xs
        forM [1..n] $ \i -> do
            j <- randomRIO (i,n)
            vi <- readArray ar i
            vj <- readArray ar j
            writeArray ar j vi
            return vj
  where
    n = length xs
    newArrayLocal :: Int -> [a] -> IO (IOArray Int a)
    newArrayLocal  nL xsL =  newListArray (1,nL) xsL

-- | shuffleInt takes a seed, number of replicates and a list of Ints and 
-- repeately shuffles the order 
shuffleInt :: Int -> Int -> [Int] -> [[Int]]
shuffleInt seed numReplicates inIntList =
    if null inIntList then []
    else if numReplicates < 1 then []
    else 
        let randList = take (length inIntList) $ randomIntList seed 
            pairList = L.sortOn fst $ zip randList inIntList
            (_, newList) = unzip pairList
        in
        newList : shuffleInt (seed + 1) (numReplicates - 1) inIntList


{-# NOINLINE randomList #-}
-- | randomList generates an infinite random list from a seed--no IO or ST monad 
-- but needs a good seed--perhaps system tiem
-- can cast to to other types like :: [Int]
randomList :: Int -> [Double]
randomList seed = randoms (mkStdGen seed) :: [Double]

{-# NOINLINE randomIntList #-}
-- | randomIntList generates an infinite random list of Ints 
randomIntList :: Int -> [Int]
randomIntList seed = randoms (mkStdGen seed) :: [Int]

{-# NOINLINE permuteList #-}
-- | permuteList ranomzes list order with seed 
permuteList :: Int -> [a] -> [a]
permuteList rSeed inList =
    if null inList then []
    else if length inList == 1 then inList
    else 
        fst $ unzip $ L.sortOn snd $ zip inList (randomIntList rSeed)

-- | takeRandom premutes a list and takes a number b ased on sed aned number to take
takeRandom :: Int -> Int -> [a] -> [a]
takeRandom rSeed number inList =
    if null inList then []
    else if number >= length inList then inList
    else
        L.take number $ permuteList rSeed inList

-- | takeNth takes n elments (each nth) of a list of length m
takeNth :: Int -> [a] -> [a]
takeNth number inList =
    if null inList then []
    else if number == 0 then []
    else if number == 1 then [head inList]
    else if number >= length inList then inList
    else 
        let (value, _) = divMod (length inList) number
            indexList = [0..(length inList - 1)]
            (_, remList) = unzip $ zipWith divMod indexList (L.replicate (length inList) value) 
            (outList, _) = unzip $ filter ((== 1) . snd) $ zip inList remList
        in
        take number outList

-- | getRandomElement returns the nth random element uniformly
-- at random
getRandomElement :: Int -> [a] -> a
getRandomElement rVal inList = 
    if null inList then error "Null list in getRandomElement"
    else if length inList == 1 then head inList
    else 
        let (_, idx) = divMod (abs rVal) (length inList)
        in
        inList !! idx

-- | selectListCostPairs is general to list of (a, Double)
-- but here used for graph sorting and selecting)takes a pair of graph representation (such as String or fgl graph), and
-- a Double cost and returns the whole of number of 'best', 'unique' or  'random' cost
-- need an Eq function such as '==' for Strings or equal for fgl
-- assumes options are all lower case
-- options are pairs of String and number for number or graphs to keeep, if number is set to (-1) then all are kept
-- if the numToKeep to return graphs is lower than number of graphs, the "best" number are returned
-- except for random.
selectListCostPairs :: forall a . (a -> a -> Bool) -> [(a, Double)] -> [String] -> Int -> Int -> [(a, Double)] 
selectListCostPairs compFun pairList optionList numToKeep seed = 
  if null optionList then error "No options specified for selectGraphCostPairs"
  else if null pairList then []
  else 
    let firstPass =
          let compFunPair :: forall b c. (a, b) -> (a, c) -> Bool
              compFunPair x = compFun (fst x) . fst
          in  if ("unique" `elem` optionList)
              then L.nubBy compFunPair pairList
              else pairList
        secondPass
          | ("best" `elem` optionList) = reverse $ L.sortOn snd firstPass
          | otherwise = firstPass
    in
    if ("random" `notElem` optionList) then take numToKeep secondPass 
    else -- shuffling with hash of structure as seed (not the best but simple for here)
      let randList = randomList seed
          pairListWRand = zip randList secondPass
          thirdPass = fmap snd $ L.sortOn fst pairListWRand

      in  take numToKeep thirdPass


-- | getSystemTimeSeconds gets teh syste time and returns IO Int
getSystemTimeSeconds :: IO Int
getSystemTimeSeconds = do
    systemTime <- getCurrentTime  
    let timeD = (round $ utcTimeToPOSIXSeconds systemTime) :: Int
    return timeD

{-# NOINLINE getSystemTimeNDT #-}
-- | getSystemTimeNDT gets the syste time and returns IO NominalDiffTime
getSystemTimeNDT :: IO NominalDiffTime
getSystemTimeNDT = do
    systemTime <- getCurrentTime  
    let !timeD = utcTimeToPOSIXSeconds systemTime
    return timeD

{-# NOINLINE getSystemTimeNDTUnsafe #-}
-- | getSystemTimeNDTUnsafe gets the syste time and returns IO NominalDiffTime
getSystemTimeNDTUnsafe :: NominalDiffTime
getSystemTimeNDTUnsafe = unsafePerformIO getSystemTimeNDT


{-# NOINLINE getSystemTimeSecondsUnsafe #-}
-- | getSystemTimeSecondsUnsafe gets the system time and returns Int via unsafePerformIO
-- without the NOINLINE the function would probbaly be comverted to a 
-- constant which would be "safe" and OK as a random seed or if only called once
getSystemTimeSecondsUnsafe :: Int
getSystemTimeSecondsUnsafe = unsafePerformIO $ force <$> getSystemTimeSeconds
    
-- | stringToInt converts a String to an Int
stringToInt :: String -> String -> Int
stringToInt fileName inStr = 
  let result = readMaybe inStr :: Maybe Int
  in
  if result == Nothing then errorWithoutStackTrace ("\n\n'Read' 'tcm' format error non-Integer value " ++ inStr ++ " in " ++ fileName)
  else fromJust result

-- | makeIndexPairs takes n and creates upper triangular matrix pairs (0,m)
makeIndexPairs :: Bool -> Int -> Int -> Int -> Int -> [(Int, Int)]
makeIndexPairs doDiagValues numI numJ indexI indexJ =
    if indexI == numI then []
    else if indexJ == numJ then makeIndexPairs doDiagValues numI numJ (indexI + 1) 0
    else 
        if doDiagValues && (indexI == indexJ) then (indexI, indexJ) : makeIndexPairs doDiagValues numI numJ indexI (indexJ + 1)
        else if (indexI < indexJ) then (indexI, indexJ) : makeIndexPairs doDiagValues numI numJ indexI (indexJ + 1)
        else makeIndexPairs doDiagValues numI numJ indexI (indexJ + 1)

-- | stripString  removes leading and trailing spaces from String 
-- akin to Text 'strip'
stripString :: String -> String
stripString inString = 
    if null inString then inString
    else 
        let firstS = dropWhile (== ' ') inString
            secondS = dropWhile (== ' ') $ reverse firstS
        in
        reverse secondS

-- | replaceVal replaces first value with second value e.g.  carriage return '\r' with line newlinme '\n'
-- call with [] accumulator
replaceVal :: (Eq a) => a -> a -> [a] -> [a] -> [a]  
replaceVal target replacement inList curList =
    if null inList then reverse curList
    else 
        let firstVal = head inList
        in
        if firstVal == target then replaceVal target replacement (tail inList) (replacement : curList)
        else replaceVal target replacement (tail inList) (firstVal : curList)


-- | cartProd takes two lists and retuns carteian product as list of pairs
cartProd :: [a] -> [b] -> [(a,b)]
cartProd xs ys = [(x,y) | x <- xs, y <- ys] 


-- | cartProdPair takes a pair of lists and retuns carteian product as list of pairs
cartProdPair :: ([a], [b]) -> [(a,b)]
cartProdPair (xs, ys) = [(x,y) | x <- xs, y <- ys] 


-- | isCompatible takes a bit vector and a list of bit vectors
-- and returns True if the fist bit vector is compatible will all in the list
isBVCompatible :: BV.BitVector -> [BV.BitVector] -> Bool
isBVCompatible inBV bvList =
    if null bvList then True
    else 
        let firstBV = head bvList
            bvVal = inBV .&. firstBV
        in
        if bvVal == inBV then isBVCompatible inBV (tail bvList)
        else if bvVal == firstBV then isBVCompatible inBV (tail bvList)
        else False

-- | textMatchWildcards takes two Text's first may have wildcards and second without
-- return True if they match, False otherwise.
textMatchWildcards :: TL.Text -> TL.Text -> Bool
textMatchWildcards straightText wildText =
    if TL.null wildText && TL.null straightText then 
        True
    else if TL.null wildText then
        False
    else if TL.null straightText then 
        False
    else if ((TL.head wildText == '*') &&  ((TL.length $ TL.dropWhile (== '*') wildText ) > 0)) && (TL.null straightText) then 
        False
    else if (TL.head wildText == '?') || (TL.head wildText == TL.head straightText) then 
        textMatchWildcards (TL.tail straightText) (TL.tail wildText)
    else if (TL.head wildText == '*') then 
        (textMatchWildcards (TL.tail straightText) wildText) || (textMatchWildcards straightText (TL.tail wildText))
    else 
        False

-- | elemWildards checks if a Text matches (without wildcards) at least one element of a List of Wildcard Text
elemWildcards :: TL.Text -> [TL.Text] -> Bool
elemWildcards straightText wildTextList =
    if null wildTextList then False
    else 
        if textMatchWildcards straightText (head wildTextList) then True
        else elemWildcards straightText (tail wildTextList)
    

-- | notElemWildcards checks if a Text matches (without wildcards) no elements of a List of Wildcard Text
notElemWildcards :: TL.Text -> [TL.Text] -> Bool
notElemWildcards straightText wildTextList =
    if null wildTextList then True
    else 
        if textMatchWildcards straightText (head wildTextList) then False
        else notElemWildcards straightText (tail wildTextList)
        
-- | getListPairs takes a list and returns all unique pairs of elements
-- order is (first found in list, second found in list)
getListPairs :: [a] -> [(a,a)]
getListPairs inList =
    if null inList then []
    else
        let firstElem = head inList
            firstPairs = zip (replicate (length $ tail inList) firstElem) (tail inList)
        in
        firstPairs ++ (getListPairs (tail inList)) 

