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

module Search.Search  ( search
                      ) where

import qualified Commands.Transform           as TRANS
import qualified Commands.Verify              as VER
import           Control.Concurrent.Async
import           Control.DeepSeq
import           Data.Bifunctor               (bimap)
import           Data.Char
import           Data.Foldable
import qualified Data.List                    as L
import           Data.Maybe
import           GeneralUtilities
import qualified GraphOptimization.Traversals as T
import qualified Graphs.GraphOperations       as GO
import qualified Search.Build                 as B
import qualified Search.Refinement            as R
import           System.Timing
import           Text.Read
import           Types.Types
import           Debug.Trace
import qualified ParallelUtilities            as PU

-- | A strict, three-way version of 'uncurry'.
uncurry3' :: (Functor f, NFData d) => (a -> b -> c -> f d) -> (a, b, c) -> f d
uncurry3' f (a, b, c) = force <$> f a b c

-- | search timed randomized search returns graph list and comment list with info String for each search instance
search :: [Argument] -> GlobalSettings -> ProcessedData -> [[VertexCost]] -> Int -> [PhylogeneticGraph] -> IO ([PhylogeneticGraph], [[String]])
search inArgs inGS inData pairwiseDistances rSeed inGraphList =
   let (searchTime, keepNum, instances) = getSearchParams inArgs
       threshold   = fromSeconds . fromIntegral $ (9 * searchTime) `div` 10
       searchTimed = uncurry3' $ searchForDuration inGS inData pairwiseDistances keepNum threshold []
       infoIndices = [1..]
       seadStreams = randomIntList <$> randomIntList rSeed
   in 
   trace ("Randomized seach for " ++ (show searchTime) ++ " seconds with " ++ (show instances) ++ " instances keeping at most " ++ (show keepNum) ++ " graphs") (
            
        do  --  threadCount <- (max 1) <$> getNumCapabilities
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
    )

-- this CPUtime is total over all threads--not wall clock
searchForDuration :: GlobalSettings -> ProcessedData -> [[VertexCost]] -> Int -> CPUTime -> [String] -> Int -> [Int] -> ([PhylogeneticGraph], [String]) -> IO ([PhylogeneticGraph], [String])
searchForDuration inGS inData pairwiseDistances keepNum allotedSeconds inCommentList refIndex seedList input@(inGraphList, infoStringList) = do
   (elapsedSeconds, output) <- timeOp $
       let result = force $ performSearch inGS inData pairwiseDistances keepNum (head seedList) input
       in  pure result
   let elapsedSecondsThreadAdjusted = fromSeconds $ fromIntegral (fst $ divMod (read (show $ toSeconds elapsedSeconds) :: Int) (read (show PU.getNumThreads) :: Int)) 
   -- let remainingTime = allotedSeconds `timeDifference` elapsedSeconds
   let remainingTime = allotedSeconds `timeDifference` elapsedSecondsThreadAdjusted
   putStrLn $ unlines [ "Thread   \t" <> show refIndex
                      , "Alloted  \t" <> show allotedSeconds
                      , "Ellapsed \t" <> show elapsedSeconds
                      , "Remaining\t" <> show remainingTime
                      ]
   if elapsedSeconds >= allotedSeconds
   then pure output
   else searchForDuration inGS inData pairwiseDistances keepNum remainingTime (inCommentList ++ (snd output)) refIndex (tail seedList) $ bimap (inGraphList <>) (infoStringList <>) output


-- | perform search takes in put graphs and performs randomized build and search with time limit
-- if run completres before 90% of time limit then will keep going
performSearch :: GlobalSettings -> ProcessedData -> [[VertexCost]] -> Int -> Int -> ([PhylogeneticGraph], [String]) -> ([PhylogeneticGraph], [String])
performSearch inGS' inData' pairwiseDistances keepNum rSeed (inGraphList', _) =
      -- set up basic parameters for search/refine methods
      let randIntList = randomIntList rSeed
          buildType = getRandomElement (randIntList !! 0) ["distance", "character"]
          buildMethod = getRandomElement (randIntList !! 1) ["unitary", "block"]

          -- build options
          numToCharBuild = keepNum
          numToDistBuild = (100 :: Int)

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

          -- swap options
          swapType = getRandomElement (randIntList !! 3) ["nni", "spr", "tbr"]
          swapKeep = keepNum
          swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep)]


          -- fuse options
          fuseSwap = getRandomElement (randIntList !! 8) ["nni", "spr", "tbr"]
          fusePairs = getRandomElement (randIntList !! 9) ["20", "40", "100"]
          fuseKeep = 2 * keepNum
          fuseArgs = [(fuseSwap, ""), ("steepest",""), ("unique",""), ("atrandom", ""), ("pairs", fusePairs), ("keep", show fuseKeep)]

          -- net edit options
          netGenArgs = [("keep", show keepNum), ("steepest", ""), ("atRandom", "")]

          -- net move options
          netMoveArgs = ("netMove", "") : netGenArgs

          -- net add options
          netAddArgs = ("netAdd", "") : netGenArgs

          -- net delete options
          netDelArgs = ("netDel", "") : netGenArgs


          -- Genetic Algorithm Arguments
          popSize = getRandomElement (randIntList !! 11) ["10", "20", "40"]
          generations = getRandomElement (randIntList !! 12) ["1"] -- , "2" , "4"]
          severity = getRandomElement (randIntList !! 13) ["0.0", "1.0", "2.0"]
          recombinations = getRandomElement (randIntList !! 14) ["20", "40", "100"]

          gaArgs = [("popsize", popSize), ("generations", generations), ("severity", severity), ("recombinations", recombinations)]

          -- drift options
          maxChanges = getRandomElement (randIntList !! 15) ["5", "10", "15"]
          acceptEqual = getRandomElement (randIntList !! 16) ["0.1", "0.5", "0.75"]
          acceptWorse = getRandomElement (randIntList !! 17) ["0.0", "2.0", "10.0"]

          drfitArgs = [("drift", ""),("maxChanges", maxChanges), ("acceptEqual", acceptEqual), ("acceptWorse", acceptWorse)]

          -- simulated annealing options
          tempSteps = getRandomElement (randIntList !! 18) ["5", "10", "15"]

          simulatedAnnealArgs = [("annealing", ""),("steps", tempSteps)]

          -- swap with drift arguments
          swapDriftArgs = swapArgs ++ drfitArgs

          -- swap with simulated anneling options
          swapAnnealArgs = swapArgs ++ simulatedAnnealArgs

          -- choose staticApproximation or not
          transformToStaticApproximation = (not . null) inGraphList' && getRandomElement (randIntList !! 19) [True, False, False]
          ((inGS, _, inData, inGraphList), staticApproxString) = if transformToStaticApproximation then
                                                                    (TRANS.transform [("staticapprox",[])] inGS' inData' inData' 0 inGraphList', "StaticApprox ")
                                                                 else ((inGS', inData', inData', inGraphList'), "")


      in

      -- no input graphs so must build and refine a bit to start
      if null inGraphList then

         -- do build + Net add (if network) + swap to crete initial solutions
         let buildGraphs = B.buildGraph buildArgs inGS inData pairwiseDistances (randIntList !! 4)
             uniqueBuildGraphs = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] buildGraphs

             netAddGraphs = if (graphType inGS == Tree) then uniqueBuildGraphs
                            else R.netEdgeMaster netAddArgs inGS inData (randIntList !! 10) uniqueBuildGraphs

             swapGraphs = R.swapMaster swapArgs inGS inData (randIntList !! 5) netAddGraphs

             uniqueSwapGraphs = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] swapGraphs
             searchString = if (graphType inGS == Tree) then "Build " ++ (L.intercalate "," $ fmap showArg buildArgs) ++ " Swap " ++ (L.intercalate "," $ fmap showArg  swapArgs)
                            else "Build " ++ (L.intercalate "," $ fmap showArg buildArgs) ++ " Net Add " ++ (L.intercalate "," $ fmap showArg netAddArgs) ++ " Swap " ++ (L.intercalate "," $ fmap showArg swapArgs)
         in  (uniqueSwapGraphs, [searchString])
         -- performSearch inArgs inGS inData pairwiseDistances keepNum startTime searchTime (searchString : infoStringList) (randIntList !! 6) uniqueSwapGraphs


      -- already have some input graphs
      -- choose a method and paramteres at random
      else
         let operation = if (graphType inGS == Tree) then getRandomElement (randIntList !! 7) ["buildSwap","fuse", "GeneticAlgorithm", "swapAnneal", "swapDrift"]
                         else getRandomElement (randIntList !! 7) ["buildSwap","fuse", "GeneticAlgorithm", "swapAnneal", "swapDrift", "netAdd", "netDelete"] -- , "netMove"] -- add/del/move edges with and without drifting

             -- this so 1/2 time annealing
             saDrift = getRandomElement (randIntList !! 19) ["noSA", "noSA", "drift", "anneal"]

             -- for rediagnosis after static approx
             -- sequential rediagnosis since assuming that if in parallel teh search operations are running in parallel
             pruneEdges = False
             warnPruneEdges = False
             startVertex = Nothing
         in
         if operation == "buildSwap" then
            let buildGraphs = B.buildGraph buildArgs inGS inData pairwiseDistances (randIntList !! 4)
                uniqueBuildGraphs = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] buildGraphs
                swapGraphs = R.swapMaster swapArgs inGS inData (randIntList !! 5) uniqueBuildGraphs
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (swapGraphs ++ inGraphList)
                uniqueGraphs = if not transformToStaticApproximation then uniqueGraphs'
                               else fmap (T.multiTraverseFullyLabelGraph inGS' inData' pruneEdges warnPruneEdges startVertex) (fmap fst6 uniqueGraphs')
                searchString = staticApproxString ++ "Build " ++ (L.intercalate "," $ fmap showArg buildArgs) ++ " Swap " ++ (L.intercalate "," $ fmap showArg  swapArgs)
            in
            (uniqueGraphs, [searchString])

         else if operation == "fuse" then
            let fuseGraphs = R.fuseGraphs fuseArgs inGS inData (randIntList !! 10) inGraphList
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (fuseGraphs ++ inGraphList)
                uniqueGraphs = if not transformToStaticApproximation then uniqueGraphs'
                               else fmap (T.multiTraverseFullyLabelGraph inGS' inData' pruneEdges warnPruneEdges startVertex) (fmap fst6 uniqueGraphs')
                searchString = staticApproxString ++ "Fuse " ++ (L.intercalate "," $ fmap showArg  fuseArgs)
            in
            (uniqueGraphs, [searchString])

         else if operation == "GeneticAlgorithm" then
            let gaGraphs = R.geneticAlgorithmMaster gaArgs inGS inData (randIntList !! 10) inGraphList
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (gaGraphs ++ inGraphList)
                uniqueGraphs = if not transformToStaticApproximation then uniqueGraphs'
                               else fmap (T.multiTraverseFullyLabelGraph inGS' inData' pruneEdges warnPruneEdges startVertex) (fmap fst6 uniqueGraphs')
                searchString = staticApproxString ++ "Genetic Algorithm " ++ (L.intercalate "," $ fmap showArg  gaArgs)
            in
            (uniqueGraphs, [searchString])

         else if operation == "swapDrift" then
            let swapDriftGraphs = R.swapMaster swapDriftArgs inGS inData (randIntList !! 10) inGraphList
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (swapDriftGraphs ++ inGraphList)
                uniqueGraphs = if not transformToStaticApproximation then uniqueGraphs'
                               else fmap (T.multiTraverseFullyLabelGraph inGS' inData' pruneEdges warnPruneEdges startVertex) (fmap fst6 uniqueGraphs')
                searchString = staticApproxString ++ "SwapDrift " ++ (L.intercalate "," $ fmap showArg  swapDriftArgs)
            in
            (uniqueGraphs, [searchString])

         else if operation == "swapAnneal" then
            let swapAnnealGraphs = R.swapMaster swapAnnealArgs inGS inData (randIntList !! 10) inGraphList
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (swapAnnealGraphs ++ inGraphList)
                uniqueGraphs = if not transformToStaticApproximation then uniqueGraphs'
                               else fmap (T.multiTraverseFullyLabelGraph inGS' inData' pruneEdges warnPruneEdges startVertex) (fmap fst6 uniqueGraphs')
                searchString = staticApproxString ++ "SwapAnneal " ++ (L.intercalate "," $ fmap showArg  swapAnnealArgs)
            in
            (uniqueGraphs, [searchString])

         else if operation == "netMove" then
            let netMoveArgs' = if saDrift == "noSA" then netMoveArgs
                               else if saDrift == "drift" then netMoveArgs ++ drfitArgs
                               else netMoveArgs ++ simulatedAnnealArgs
                netMoveGraphs = R.netEdgeMaster netMoveArgs' inGS inData (randIntList !! 10) inGraphList
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (netMoveGraphs ++ inGraphList)
                uniqueGraphs = if not transformToStaticApproximation then uniqueGraphs'
                               else fmap (T.multiTraverseFullyLabelGraph inGS' inData' pruneEdges warnPruneEdges startVertex) (fmap fst6 uniqueGraphs')
                searchString = staticApproxString ++ "NetMove " ++ (L.intercalate "," $ fmap showArg  netMoveArgs')
            in
            (uniqueGraphs, [searchString])

         else if operation == "netAdd" then
            let netAddArgs' = if saDrift == "noSA" then netAddArgs
                              else if saDrift == "drift" then netAddArgs ++ drfitArgs
                              else netAddArgs ++ simulatedAnnealArgs
                netAddGraphs = R.netEdgeMaster netAddArgs' inGS inData (randIntList !! 10) inGraphList
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (netAddGraphs ++ inGraphList)
                uniqueGraphs = if not transformToStaticApproximation then uniqueGraphs'
                               else fmap (T.multiTraverseFullyLabelGraph inGS' inData' pruneEdges warnPruneEdges startVertex) (fmap fst6 uniqueGraphs')
                searchString = staticApproxString ++ "NetAdd " ++ (L.intercalate "," $ fmap showArg  netAddArgs')
            in
            (uniqueGraphs, [searchString])

         else if operation == "netDelete" then
            let netDelArgs' = if saDrift == "noSA" then netDelArgs
                              else if saDrift == "drift" then netDelArgs ++ drfitArgs
                              else netDelArgs ++ simulatedAnnealArgs
                netDelGraphs = R.netEdgeMaster netDelArgs' inGS inData (randIntList !! 10) inGraphList
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (netDelGraphs ++ inGraphList)
                uniqueGraphs = if not transformToStaticApproximation then uniqueGraphs'
                               else fmap (T.multiTraverseFullyLabelGraph inGS' inData' pruneEdges warnPruneEdges startVertex) (fmap fst6 uniqueGraphs')
                searchString = staticApproxString ++ "netDelete " ++ (L.intercalate "," $ fmap showArg  netDelArgs')
            in
            (uniqueGraphs, [searchString])


         else error ("Unknown/unimplemented method in search: " ++ operation)
      where showArg a = "(" ++ (fst a) ++ "," ++ (snd a) ++ ")"

-- | getSearchParams takes arguments and returns search params
getSearchParams :: [Argument] -> (Int, Int, Int)
getSearchParams inArgs =
   let fstArgList = fmap (fmap toLower . fst) inArgs
       sndArgList = fmap (fmap toLower . snd) inArgs
       lcArgList = zip fstArgList sndArgList
       checkCommandList = checkCommandArgs "search" fstArgList VER.searchArgList
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

