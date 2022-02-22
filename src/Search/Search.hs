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

import Types.Types
import qualified ParallelUtilities       as PU
import Control.Parallel.Strategies
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


-- | search arguments
searchArgList :: [String]
searchArgList = ["days", "hours", "minutes", "seconds"]

-- | search timed randomized search
search :: [Argument] -> GlobalSettings -> ProcessedData -> [[VertexCost]] -> Int -> [PhylogeneticGraph] -> [PhylogeneticGraph]
search inArgs inGS inData pairwiseDistances rSeed inGraphList = 
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
         -- spawns as many 'search's as parallel threads
         let numThreads = if PU.getNumThreads > 0 then  PU.getNumThreads
                          else 1
             randIntList = randomIntList rSeed
             seconds' = if ((fromJust minutes > 0) || (fromJust hours > 0) || (fromJust days > 0)) && (null secondsList) then Just 0
                        else seconds
             searchTime = (fromJust seconds') + (60 * (fromJust minutes)) + (3600 * (fromJust hours))
             startTime = getSystemTimeSecondsUnsafe
             newGraphList = zipWith (performSearch inArgs inGS inData pairwiseDistances (fromJust keepNum) startTime searchTime) randIntList (replicate numThreads inGraphList) `using` PU.myParListChunkRDS
         in
         GO.selectPhylogeneticGraph [("unique", (show $ fromJust keepNum))] 0 ["unique"] $ (inGraphList ++ concat newGraphList)

-- | perform search takes in put graphs and performs randomized build and search with time limit
-- if run completres before 90% of time limit then will keep going
performSearch :: [Argument] -> GlobalSettings -> ProcessedData -> [[VertexCost]] -> Int -> Int -> Int -> Int -> [PhylogeneticGraph] -> [PhylogeneticGraph]
performSearch inArgs inGS inData pairwiseDistances keepNum startTime searchTime rSeed inGraphList = 
   let currentTime = getSystemTimeSecondsUnsafe
   in
   --check time to see if should stop
   if (currentTime - startTime) > floor (0.9 * (fromIntegral searchTime)) then inGraphList

   -- time left go on
   else 
      trace ("Time left: " ++ (show (searchTime - (currentTime - startTime))) ++ " seconds") (
      
      -- for use later
      let randIntList = randomIntList rSeed
          buildType = getRandomElement (randIntList !! 0) ["distance", "character"]
          buildMethod = getRandomElement (randIntList !! 1) ["unitary", "block"]
             
          -- general build options
          numToCharBuild = 10
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
         in
         performSearch inArgs inGS inData pairwiseDistances keepNum startTime searchTime (randIntList !! 6) uniqueSwapGraphs


      -- already have some input gaphs
      else 
         inGraphList
      )