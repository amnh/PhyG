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
import Control.DeepSeq
import GeneralUtilities
import qualified Graphs.GraphOperations  as GO
import Utilities.Utilities               as U
import Data.Maybe
import           Data.Char
import           Text.Read
import qualified Search.Build as B
import qualified Search.Refinement as R
import System.Timing
import Data.Bifunctor (bimap)
import Data.Foldable
import Control.Concurrent.Async


-- | A strict, three-way version of 'uncurry'.
uncurry3' :: (Functor f, NFData d) => (a -> b -> c -> f d) -> (a, b, c) -> f d
uncurry3' f (a, b, c) = force <$> f a b c


-- | search arguments
searchArgList :: [String]
searchArgList = ["days", "hours", "minutes", "seconds", "instances"]


-- | search timed randomized search returns graph list and comment list with info String for each serch instance
search :: [Argument] -> GlobalSettings -> ProcessedData -> [[VertexCost]] -> Int -> [PhylogeneticGraph] -> IO ([PhylogeneticGraph], [[String]])
search inArgs inGS inData pairwiseDistances rSeed inGraphList = 
   let (searchTime, keepNum, instances) = getSearchParams inArgs
       threshold   = fromSeconds . fromIntegral $ (9 * searchTime) `div` 10
       searchTimed = uncurry3' $ searchForDuration inGS inData pairwiseDistances keepNum threshold
       infoIndices = [1..]
       seadStreams = randomIntList <$> randomIntList rSeed        
   in  do  --  threadCount <- (max 1) <$> getNumCapabilities
           let threadCount = instances -- <- (max 1) <$> getNumCapabilities
           let startGraphs = replicate threadCount (inGraphList, mempty)
           let threadInits = zip3 infoIndices seadStreams startGraphs
           resultList <- mapConcurrently searchTimed threadInits
           pure $
               let (newGraphList, commentList) = unzip resultList
                   completeGraphList = inGraphList <> fold newGraphList
                   filteredGraphList = GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] completeGraphList
                   selectedGraphList = take keepNum filteredGraphList
               in  (selectedGraphList, commentList)


searchForDuration :: GlobalSettings -> ProcessedData -> [[VertexCost]] -> Int -> CPUTime -> Int -> [Int] -> ([PhylogeneticGraph], [String]) -> IO ([PhylogeneticGraph], [String])
searchForDuration inGS inData pairwiseDistances keepNum allotedSeconds refIndex seedList input@(inGraphList, infoStringList) = do
   (elapsedSeconds, output) <- timeOp $ 
       let result = force $ performSearch inGS inData pairwiseDistances keepNum (head seedList) input
       in  pure result
   let remainingTime = allotedSeconds `timeDifference` elapsedSeconds
   putStrLn $ unlines [ "Thread   \t" <> show refIndex
                      , "Alloted  \t" <> show allotedSeconds
                      , "Ellapsed \t" <> show elapsedSeconds
                      , "Remaining\t" <> show remainingTime
                      ]
   if   elapsedSeconds >= allotedSeconds
   then pure output
   else searchForDuration inGS inData pairwiseDistances keepNum remainingTime refIndex (tail seedList) $ bimap (inGraphList <>) (infoStringList <>) output 


-- | perform search takes in put graphs and performs randomized build and search with time limit
-- if run completres before 90% of time limit then will keep going
performSearch :: GlobalSettings -> ProcessedData -> [[VertexCost]] -> Int -> Int -> ([PhylogeneticGraph], [String]) -> ([PhylogeneticGraph], [String])
performSearch inGS inData pairwiseDistances keepNum rSeed (inGraphList, infoStringList) = 
      -- for use later
      let randIntList = randomIntList rSeed
          buildType = getRandomElement (randIntList !! 0) ["distance", "character"]
          buildMethod = getRandomElement (randIntList !! 1) ["unitary", "block"]
             
          -- general build options
          numToCharBuild = 10  :: Int
          numToDistBuild = 100 :: Int

          -- build block options
          reconciliationMethod = getRandomElement (randIntList !! 2) ["eun", "cun"]

          distOptions = if buildType == "distance" then 
                           if buildMethod == "block" then [("replicates", show numToDistBuild), ("rdwag", ""), ("best", show (1 :: Int))]
                           else  [("replicates", show numToDistBuild), ("rdwag", ""), ("best", show numToCharBuild)]
                        else 
                           if buildMethod == "block" then [("replicates", show numToCharBuild)]
                           else [("replicates", show (1 :: Int))]

          blockOptions = if buildMethod == "block" then [("block", ""),("atRandom", ""),("displaytrees", show numToCharBuild),(reconciliationMethod, "")]
                         else []

          buildArgs = [(buildType, "")] ++ distOptions ++ blockOptions

          --swap options 
          swapType =  getRandomElement (randIntList !! 3) ["nni", "spr", "tbr"]
          swapKeep = 10 :: Int

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
         let operation = getRandomElement (randIntList !! 7) ["buildSwap"]
         in  
         if operation == "buildSwap" then
            let buildGraphs = B.buildGraph buildArgs inGS inData pairwiseDistances (randIntList !! 4)
                uniqueBuildGraphs = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] buildGraphs
                swapGraphs = R.swapMaster swapArgs inGS inData (randIntList !! 5) uniqueBuildGraphs
                uniqueSwapGraphs = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] swapGraphs
                searchString = "Build " ++ show buildArgs ++ " Swap " ++ show swapArgs
         in  
         (uniqueSwapGraphs, [searchString])

         
              
         else (inGraphList, infoStringList)
      

-- | getSearchParams takes arguments and returns search params
getSearchParams :: [Argument] -> (Int, Int, Int)
getSearchParams inArgs = 
   let fstArgList = fmap (fmap toLower . fst) inArgs
       sndArgList = fmap (fmap toLower . snd) inArgs
       lcArgList = zip fstArgList sndArgList
       checkCommandList = U.checkCommandArgs "search" fstArgList searchArgList
   in
   -- check for valid command options
   if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'search': " ++ show inArgs)
   else
      let instancesList = filter ((=="instances").fst) lcArgList
          instances
            | length instancesList > 1 =
              errorWithoutStackTrace ("Multiple 'keep' number specifications in search command--can have only one: " ++ show inArgs)
            | null instancesList = Just 1
            | otherwise = readMaybe (snd $ head instancesList) :: Maybe Int

          keepList = filter ((=="keep").fst) lcArgList
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
      else if isNothing instances then errorWithoutStackTrace ("Instnaces specification not an integer in search: "  ++ show (head instancesList))
      else if isNothing days then errorWithoutStackTrace ("Days specification not an integer in search: "  ++ show (head daysList))
      else if isNothing hours then errorWithoutStackTrace ("Hours specification not an integer in search: "  ++ show (head hoursList))
      else if isNothing minutes then errorWithoutStackTrace ("Minutes specification not an integer in search: "  ++ show (head minutesList))
      else if isNothing seconds then errorWithoutStackTrace ("seconds factor specification not an integer in search: "  ++ show (head secondsList))

      else 
         let seconds' = if ((fromJust minutes > 0) || (fromJust hours > 0) || (fromJust days > 0)) && (null secondsList) then Just 0
                        else seconds
             searchTime = (fromJust seconds') + (60 * (fromJust minutes)) + (3600 * (fromJust hours))
         in
         (searchTime, fromJust keepNum, fromJust instances)
             
