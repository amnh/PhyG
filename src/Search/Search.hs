{- |
Module      :  Search.hs
Description :  Module controlling timed randomized search functions
Copyright   :  (c) 2022 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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

module Search.Search  (search
                      ) where

import           Control.Concurrent
import Types.Types
import qualified ParallelUtilities       as PU
import Control.Parallel.Strategies
import Control.DeepSeq
import Debug.Trace
import GeneralUtilities
import qualified Graphs.GraphOperations  as GO
import Utilities.Utilities               as U
import Data.Maybe
import           Data.Char
import           Text.Read
import qualified Search.Fuse as F
import qualified Search.GeneticAlgorithm as GA
import qualified Search.Build as B
import qualified Search.Swap as S
import qualified Search.Refinement as R
import System.Timing
import Data.Bifunctor (bimap)
import Data.Traversable
import Data.Foldable
import Control.Monad
import Data.Vector (fromListN)


traverseInParallelWith :: (Foldable f, NFData c) => (a -> b -> IO c) -> f a -> f b -> IO [c]
traverseInParallelWith f xs ys = sequenceA $ parMap strat (uncurry f) $ zip (toList xs) (toList ys)


strat :: NFData a => Strategy (IO a)
strat x = pure $ do 
   r <- x
   rnf r `seq` pure r


-- | search arguments
searchArgList :: [String]
searchArgList = ["days", "hours", "minutes", "seconds"]


-- | search timed randomized search
search :: [Argument] -> GlobalSettings -> ProcessedData -> [[VertexCost]] -> Int -> [PhylogeneticGraph] -> IO [PhylogeneticGraph]
search inArgs inGS inData pairwiseDistances rSeed inGraphList = 
   let (searchTime, keepNum) = getSearchParams inArgs
       threshold   = fromSeconds . fromIntegral $ (9 * searchTime) `div` 10 
   in  do  threadCount <- (max 1) <$> getNumCapabilities
           putStrLn $ unwords [ "threadCount:", show threadCount ]
           let seedList  = take threadCount $ randomIntList rSeed
           let graphList = replicate threadCount (inGraphList, [])
           resultList <- traverseInParallelWith (searchForDuration inArgs inGS inData pairwiseDistances keepNum threshold) seedList graphList
           let (newGraphList, commentList) = unzip resultList
           pure . take keepNum . GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] $ (inGraphList ++ concat newGraphList)


searchForDuration :: [Argument] -> GlobalSettings -> ProcessedData -> [[VertexCost]] -> Int -> CPUTime -> Int -> ([PhylogeneticGraph], [String]) -> IO ([PhylogeneticGraph], [String])
searchForDuration inArgs inGS inData pairwiseDistances keepNum allotedSeconds rSeed input@(inGraphList, infoStringList) = do
   (elapsedSeconds, output) <- timeOp $ 
       let result = force $ performSearch inArgs inGS inData pairwiseDistances keepNum rSeed input
       in  pure result
   let remainingTime = allotedSeconds `timeDifference` elapsedSeconds
   putStrLn $ unlines ["Alloted: " <> show allotedSeconds, "Ellapsed: " <> show elapsedSeconds, "Remaining: " <> show remainingTime]
   if   elapsedSeconds >= allotedSeconds
   then pure output
   else searchForDuration inArgs inGS inData pairwiseDistances keepNum remainingTime rSeed $ bimap (inGraphList <>) (infoStringList <>) output 


-- | perform search takes in put graphs and performs randomized build and search with time limit
-- if run completres before 90% of time limit then will keep going
performSearch :: [Argument] -> GlobalSettings -> ProcessedData -> [[VertexCost]] -> Int -> Int -> ([PhylogeneticGraph], [String]) -> ([PhylogeneticGraph], [String])
performSearch inArgs inGS inData pairwiseDistances keepNum rSeed (inGraphList, infoStringList) = 
      -- for use later
      let randIntList = randomIntList rSeed
          buildType = getRandomElement (randIntList !! 0) ["distance", "character"]
          buildMethod = getRandomElement (randIntList !! 1) ["unitary"] -- , "block"]
             
          -- general build options
          numToCharBuild = 5
          numToDistBuild = 100
          numToKeep = numToCharBuild


          -- build block options
          reconciliationMethod = getRandomElement (randIntList !! 2) ["eun", "cun"]

          distOptions = if buildType == "distance" then [("replicates", show numToDistBuild), ("best", show numToKeep), ("rdwag", "")]
                        else [("replicates", show numToCharBuild)]

          blockOptions = if buildMethod == "block" then [("block", ""),("atRandom", ""),("displaytrees", show numToCharBuild),(reconciliationMethod, "")]
                         else []

          buildArgs = [(buildType, "")] ++ distOptions ++ blockOptions

          --swap options 
          swapType =  getRandomElement (randIntList !! 3) ["nni", "spr", "tbr"]
          swapKeep = 10

          swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep)]
      in

      -- no input graphs so must build and refine a bit to start
      if null inGraphList then 
         
         -- do build
         let buildGraphs = B.buildGraph buildArgs inGS inData pairwiseDistances (randIntList !! 4)
             uniqueBuildGraphs = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] buildGraphs
             swapGraphs = R.swapMaster swapArgs inGS inData (randIntList !! 5) uniqueBuildGraphs
             uniqueSwapGraphs = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] swapGraphs
             searchString = "Build " ++ show buildArgs ++ " Swap " ++ show swapArgs
         in  (uniqueSwapGraphs, [searchString])
         -- performSearch inArgs inGS inData pairwiseDistances keepNum startTime searchTime (searchString : infoStringList) (randIntList !! 6) uniqueSwapGraphs


      -- already have some input gaphs
      else 
         let operation = getRandomElement (randIntList !! 7) ["build"]
         in  if operation == "build" then
               let buildGraphs = B.buildGraph buildArgs inGS inData pairwiseDistances (randIntList !! 8)
                   uniqueBuildGraphs = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] buildGraphs
                   searchString = "Build " ++ show buildArgs
               in  (uniqueBuildGraphs, [searchString])
             -- performSearch inArgs inGS inData pairwiseDistances keepNum startTime searchTime (searchString : infoStringList) (randIntList !! 9) uniqueBuildGraphs
              
         else (inGraphList, infoStringList)
      

-- | getSearchParams takes arguments and returns search params
getSearchParams :: [Argument] -> (Int, Int)
getSearchParams inArgs = 
   let fstArgList = fmap (fmap toLower . fst) inArgs
       sndArgList = fmap (fmap toLower . snd) inArgs
       lcArgList = zip fstArgList sndArgList
       checkCommandList = U.checkCommandArgs "search" fstArgList searchArgList
   in
   -- check for valid command options
   if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'search': " ++ show inArgs)
   else
      let keepList = filter ((=="keep").fst) lcArgList
          keepNum
            | length keepList > 1 =
              errorWithoutStackTrace ("Multiple 'keep' number specifications in search command--can have only one: " ++ show inArgs)
            | null keepList = Just 10
            | otherwise = readMaybe (snd $ head keepList) :: Maybe Int

          daysList = filter ((=="days").fst) lcArgList
          days
            | length daysList > 1 =
              errorWithoutStackTrace ("Multiple 'days' number specifications in search command--can have only one: " ++ show inArgs)
            | null daysList = Just 0
            | otherwise = readMaybe (snd $ head daysList) :: Maybe Int

          hoursList = filter ((=="hours").fst) lcArgList
          hours
            | length hoursList > 1 =
              errorWithoutStackTrace ("Multiple 'hours' number specifications in search command--can have only one: " ++ show inArgs)
            | null hoursList = Just 0
            | otherwise = readMaybe (snd $ head hoursList) :: Maybe Int

          minutesList = filter ((=="minutes").fst) lcArgList
          minutes
            | length minutesList > 1 =
              errorWithoutStackTrace ("Multiple 'minutes' number specifications in search command--can have only one: " ++ show inArgs)
            | null minutesList = Just 0
            | otherwise = readMaybe (snd $ head minutesList) :: Maybe Int

          secondsList = filter ((=="seconds").fst) lcArgList
          seconds
            | length secondsList > 1 =
              errorWithoutStackTrace ("Multiple 'seconds' number specifications in search command--can have only one: " ++ show inArgs)
            | null secondsList = Just 30
            | otherwise = readMaybe (snd $ head secondsList) :: Maybe Int

      in 
      if isNothing keepNum then errorWithoutStackTrace ("Keep specification not an integer in search: "  ++ show (head keepList))
      else if isNothing days then errorWithoutStackTrace ("Days specification not an integer in search: "  ++ show (head daysList))
      else if isNothing hours then errorWithoutStackTrace ("Hours specification not an integer in search: "  ++ show (head hoursList))
      else if isNothing minutes then errorWithoutStackTrace ("Minutes specification not an integer in search: "  ++ show (head minutesList))
      else if isNothing seconds then errorWithoutStackTrace ("seconds factor specification not an integer in search: "  ++ show (head secondsList))

      else 
         let seconds' = if ((fromJust minutes > 0) || (fromJust hours > 0) || (fromJust days > 0)) && (null secondsList) then Just 0
                        else seconds
             searchTime = (fromJust seconds') + (60 * (fromJust minutes)) + (3600 * (fromJust hours))
         in
         (searchTime, fromJust keepNum)
             