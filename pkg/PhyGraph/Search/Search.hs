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

import qualified Commands.Transform       as TRANS
import qualified Commands.Verify          as VER
import           Control.Concurrent.Async
import           Control.DeepSeq
import           Data.Bifunctor           (bimap)
import           Data.Char
import           Data.Foldable
import qualified Data.List                as L
import           Data.Maybe
import           GeneralUtilities
-- import qualified GraphOptimization.Traversals as T
import qualified Data.Vector              as V
import           Debug.Trace
import qualified Graphs.GraphOperations   as GO
import qualified Search.Build             as B
import qualified Search.Refinement        as R
import           System.IO
import           System.Random
import           System.Timing
import           Text.Read
import           Types.Types
import qualified Utilities.Utilities      as U



-- Add non reroot thing liike IA for faster

-- | treeBanditList is list of search types to be chosen from if graphType is tree
treeBanditList :: [String]
treeBanditList = [
                 "buildCharacter", "buildDistance", -- "buildSPR", "buildAlternate",
                 "swapSPR", "swapAlternate",
                 "fuse", "fuseSPR", "fuseTBR",
                 "driftSPR", "driftAlternate", "annealSPR", "annealAlternate",
                 "geneticAlgorithm"
                 ]

-- | netWorkBanditList is list of search types unique to graphType network
netWorkBanditList :: [String]
netWorkBanditList = ["networkAdd", "networkDelete", "networkAddDelete", "driftNetwork", "annealNetwork", "networkMove"]

-- | fullBanditList is list of search types to be chosen from if Network
fullBanditList :: [String]
fullBanditList = treeBanditList ++ netWorkBanditList

-- | A strict, three-way version of 'uncurry'.
uncurry3' :: (Functor f, NFData d) => (a -> b -> c -> f d) -> (a, b, c) -> f d
uncurry3' f (a, b, c) = force <$> f a b c

-- | search timed randomized search returns graph list and comment list with info String for each search instance
search :: [Argument] -> GlobalSettings -> ProcessedData -> [[VertexCost]] -> Int -> [PhylogeneticGraph] -> IO ([PhylogeneticGraph], [[String]])
search inArgs inGS inData pairwiseDistances rSeed inGraphList =
   let (searchTime, keepNum, instances, thompsonSample, mFactor, mFunction, maxNetEdges, stopNum) = getSearchParams inArgs

       -- flatThetaList is the initial prior list (flat) of search (bandit) choices
       -- can also be used in search for non-Thomspon search
       flatThetaList = if graphType inGS == Tree then
                            zip treeBanditList (L.replicate (length treeBanditList) (1.0 / (fromIntegral $ length treeBanditList)))
                       else
                            zip fullBanditList (L.replicate (length fullBanditList) (1.0 / (fromIntegral $ length fullBanditList)))

       threshold   = fromSeconds . fromIntegral $ (95 * searchTime) `div` 100
       initialSeconds = fromSeconds . fromIntegral $ (0 :: Int)
       searchTimed = uncurry3' $ searchForDuration inGS inData pairwiseDistances keepNum thompsonSample mFactor mFunction flatThetaList 1 maxNetEdges initialSeconds threshold 0 stopNum
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
                   newCostList = L.group $ L.sort $ fmap getMinGraphListCost newGraphList
                   iterationHitString = ("Hit minimum cost in " ++ (show $ length $ head newCostList) ++ " of " ++ (show $ length newGraphList) ++ " iterations")
                   completeGraphList = inGraphList <> fold newGraphList
                   filteredGraphList = GO.selectGraphs Unique (maxBound::Int) 0.0 (-1) completeGraphList
                   selectedGraphList = take keepNum filteredGraphList
               in  (selectedGraphList, commentList ++ [[iterationHitString]])
    )
    where getMinGraphListCost a = minimum $ fmap snd6 a

-- this CPUtime is total over all threads--not wall clock
-- so changed to crappier getCurrentTime in System.Timing to
-- get wall clock-like ellapsed time
searchForDuration :: GlobalSettings
                  -> ProcessedData
                  -> [[VertexCost]]
                  -> Int
                  -> Bool
                  -> Int
                  -> String
                  -> [(String, Double)]
                  -> Int
                  -> Int
                  -> CPUTime
                  -> CPUTime
                  -> Int
                  -> Int
                  -> Int
                  -> [Int]
                  -> ([PhylogeneticGraph], [String])
                  -> IO ([PhylogeneticGraph], [String])
searchForDuration inGS inData pairwiseDistances keepNum thompsonSample mFactor mFunction thetaList counter maxNetEdges inTotalSeconds allotedSeconds stopCount stopNum refIndex seedList (inGraphList, infoStringList) = do
   -- (elapsedSeconds, output) <- timeOpUT $
   (elapsedSeconds, elapsedSecondsCPU, output) <- timeOpCPUWall $
       let  -- this line to keep control of graph number
            inGraphList' = take keepNum $ GO.selectGraphs Unique (maxBound::Int) 0.0 (-1) inGraphList
            result = force $ performSearch inGS inData pairwiseDistances keepNum thompsonSample thetaList maxNetEdges (head seedList) inTotalSeconds (inGraphList', infoStringList)
       in  pure result

   -- update theta list based on performance
   let outTotalSeconds = timeSum inTotalSeconds elapsedSecondsCPU
   let finalTimeString = ",Final Values,,," ++ (show $ toSeconds outTotalSeconds)
   -- passing time as CPU time not wall clock so parallel timings change to elapsedSeconds for wall clock
   let (updatedThetaList, newStopCount) = updateTheta thompsonSample mFactor mFunction counter (snd output) thetaList elapsedSecondsCPU outTotalSeconds stopCount stopNum
   let thetaString = L.intercalate "," $ fmap (show . snd) updatedThetaList

   let remainingTime = allotedSeconds `timeLeft` elapsedSeconds
   hPutStrLn stderr $ unlines [ "Thread   \t" <> show refIndex
                      , "Alloted  \t" <> show allotedSeconds
                      , "Ellapsed \t" <> show elapsedSeconds
                      , "Remaining\t" <> show remainingTime
                      ]
   if (toPicoseconds remainingTime) == 0 || newStopCount >= stopNum
   then pure (fst output, infoStringList ++ (snd output) ++ [finalTimeString ++ "," ++ thetaString ++ "," ++ "*"]) -- output with strings correctly added together
   else searchForDuration inGS inData pairwiseDistances keepNum thompsonSample mFactor mFunction updatedThetaList (counter + 1) maxNetEdges outTotalSeconds remainingTime newStopCount stopNum refIndex (tail seedList) $ bimap (inGraphList <>) (infoStringList <>) output

-- | updateTheta updates the expected success parameters for the bandit search list
updateTheta :: Bool -> Int -> String -> Int -> [String] -> [(String, Double)] -> CPUTime -> CPUTime -> Int -> Int -> ([(String, Double)], Int)
updateTheta thompsonSample mFactor mFunction counter infoStringList inPairList elapsedSeconds totalSeconds stopCount stopNum =
    if null inPairList then ([], stopCount)
    else if not thompsonSample then
        let searchBandit = takeWhile (/= ',') (tail $ head infoStringList)
            searchDeltaString = takeWhile (/= ',') $ drop (1 + length searchBandit) (tail $ head infoStringList)
            searchDelta = read searchDeltaString :: Double
            newStopCount = if searchDelta <= 0.0 then stopCount + 1
                           else 0
            stopString = if newStopCount >= stopNum then ("\n\tSearch iterations have not improved for " ++ (show newStopCount) ++ " iterations--terminating this search command")
                         else ""
        in
        trace  stopString
        (inPairList, newStopCount)
    else
        -- trace ("UT1: " ++ (show infoStringList)) (
        -- update via results, previous history, memory \factor and type of memory "loss"
        let -- get search info for last search iteration
            -- thjis tail/head stuff due to leading ',' in comment field
            searchBandit = takeWhile (/= ',') (tail $ head infoStringList)
            searchDeltaString = takeWhile (/= ',') $ drop (1 + length searchBandit) (tail $ head infoStringList)
            searchDelta = read searchDeltaString :: Double

            -- get timing and benefit accounting for 0's
            -- average time ratio in time factor for benefit adjustment
            totalTime = (fromIntegral $ toSeconds totalSeconds) / (fromIntegral counter) :: Double
            timeFactor = if counter == 1 then 1.0
                         else if toSeconds elapsedSeconds == 0 then 0.1
                         else (fromIntegral $ toSeconds elapsedSeconds) / totalTime

            durationTime = if toSeconds elapsedSeconds <= 1 then 1
                           else toSeconds elapsedSeconds

            -- get bandit index and orignial theta value from pair list
            indexBandit' = L.elemIndex searchBandit $ fmap fst inPairList
            indexBandit = fromJust  indexBandit'

            inThetaBandit = snd $ inPairList !! indexBandit

            newStopCount = if searchDelta <= 0.0 then stopCount + 1
                           else 0

            stopString = if newStopCount >= stopNum then ("\n\tSearch iterations have not improved for " ++ (show newStopCount) ++ " iterations--terminating this search command")
                         else ""
        in
        -- check error
        if isNothing indexBandit' then error ("Bandit index not found: " ++ searchBandit ++ " in " ++ (show inPairList))

        -- "simple" for testing, sucess=1, no time factor, full average overall interations
        else if mFunction == "simple" then
            -- first Simple update based on 0 = no better. 1 == any better irrespective of magnitude or time, full memory
            -- all thetas (second field) sum to one.
            -- if didn't do anything but increment everyone else with 1/(num bandits - 1), then renormalize
            -- dowenweights unsiucessful bandit
            -- need to add more downweight if search long
            let previousSuccessList = fmap (* (fromIntegral counter)) $ fmap snd inPairList

                benefit = if searchDelta <= 0.0 then 0.0
                          else searchDelta / (fromIntegral durationTime)

                -- see if bandit was sucessful and if so set increment
                incrementBandit = if benefit == 0.0 then 0.0 else 1.0
                newBanditVal = incrementBandit + ((fromIntegral counter) * inThetaBandit)

                -- for nonBandit if not successful increment them (basically to downweight bandit), other wise not
                incrementNonBanditVals = if benefit == 0.0 then 1.0 / ((fromIntegral $ length inPairList) - 1.0)
                                         else 0.0
                updatedSuccessList = fmap (+ incrementNonBanditVals) previousSuccessList

                -- uopdate the bandit from list by splitting and rejoining
                firstBanditPart = take indexBandit updatedSuccessList
                thirdBanditPart = drop (indexBandit + 1) updatedSuccessList
                newSuccessList = firstBanditPart ++ (newBanditVal : thirdBanditPart)
                totalTheta = sum newSuccessList

                newThetaList = fmap (/ totalTheta) newSuccessList
            in
            -- trace ("UT2: " ++ searchBandit ++ " index " ++ (show indexBandit) ++ " total time: " ++ (show $ toSeconds totalSeconds) ++ " elapsed time: " ++ (show $ toSeconds elapsedSeconds) ++ " -> " ++ (show (searchDelta, toSeconds elapsedSeconds, benefit)) ++ "\n" ++ (show $ fmap snd inPairList) ++ "\n" ++ (show newThetaList) ++ "\n" ++ (head infoStringList)) (
            trace  stopString
            (zip (fmap fst inPairList) newThetaList, newStopCount)
            -- )

        -- more complex 'recency' options
        else if mFunction `elem` ["linear","exponential"] then
            let -- weight factors for previous theta (wN-1)* \theta
                -- maxed ot counter so averaging not wierd for early iterations
                mFactor' = min (fromIntegral counter :: Double) (fromIntegral mFactor :: Double)
                (wN_1, wN) = if mFunction == "linear" then
                                    (mFactor' / (mFactor' + 1), 1.0 / (mFactor' + 1))
                             else if mFunction == "exponential" then
                                    (1.0 - (1.0 / (2.0 ** mFactor')), (1.0 / (2.0 ** mFactor')))
                             else error ("Thompson search option " ++ mFunction ++ " not recognized " ++ (show ["simple", "linear","exponential"]))

                -- simple "success-based" benefit, scaled to average time of search iteration
                searchBenefit = if searchDelta <= 0.0 then 0.0
                                else searchDelta / timeFactor

                previousSuccessList = fmap (* (fromIntegral counter)) $ fmap snd inPairList

                -- average of new am]nd previous if m=1, no memory if m=0, longer memory with larger m
                -- m should be limited to counter
                -- linear ((m * previous value) + new value) / m+1
                -- exponential ((2^(-m)  * previous value) + new value) / (2^(-m) + 1)
                newBanditVal = if searchBenefit > inThetaBandit then
                                    (wN * searchBenefit) + (wN_1 * inThetaBandit)
                               else
                                    (wN_1 * inThetaBandit) + (wN * (inThetaBandit + searchBenefit))

                -- for nonBandit if not successful increment them (basically to downweight bandit), other wise not
                incrementNonBanditVals = if searchDelta <= 0.0 then 1.0 / ((fromIntegral $ length inPairList) - 1.0)
                                         else 0.0
                updatedSuccessList = fmap (+ incrementNonBanditVals) previousSuccessList

                -- uopdate the bandit from list by splitting and rejoining, then normalizing to 1.0
                firstBanditPart = take indexBandit updatedSuccessList
                thirdBanditPart = drop (indexBandit + 1) updatedSuccessList
                newSuccessList = firstBanditPart ++ (newBanditVal : thirdBanditPart)
                totalTheta = sum newSuccessList

                newThetaList = fmap (/ totalTheta) newSuccessList

            in
            -- trace ("Not simple: " ++ mFunction ++ " search benefit " ++ (show searchBenefit) ++ " " ++ searchBandit ++ " index " ++ (show indexBandit) ++ " total time: " ++ (show $ toSeconds totalSeconds) ++ " elapsed time: " ++ (show $ toSeconds elapsedSeconds) ++ " -> " ++ (show (searchDelta, toSeconds elapsedSeconds)) ++ "\n" ++ (show $ fmap snd inPairList) ++ "\n" ++ (show newThetaList) ++ "\n" ++ (head infoStringList)) (
            trace  stopString
            (zip (fmap fst inPairList) newThetaList, newStopCount)
            -- )


        else errorWithoutStackTrace ("Thompson search option " ++ mFunction ++ " not recognized " ++ (show ["simple", "linear","exponential"]))
        -- )


-- | performSearch takes in put graphs and performs randomized build and search with time limit
-- Thompson sampling and mFactor to pick strategy from updated theta success values
-- the random calls return the tail of the input list to avoid long list access--can do vector since infinite
performSearch :: GlobalSettings
               -> ProcessedData
               -> [[VertexCost]]
               -> Int
               -> Bool
               -> [(String, Double)]
               -> Int
               -> Int
               -> CPUTime
               -> ([PhylogeneticGraph], [String])
               -> ([PhylogeneticGraph], [String])
performSearch inGS' inData' pairwiseDistances keepNum _ thetaList maxNetEdges rSeed inTime (inGraphList', _) =
      -- set up basic parameters for search/refine methods
      let -- set up log for sample
          thompsonString = "," ++ (show thetaList)

          -- get infinite lists if integers and doubles
          randIntList = randomIntList rSeed
          randDoubleList = randoms (mkStdGen rSeed) :: [Double]

          -- this for constant access to random doubles need take for infinite list
          -- need to update as more random doubles are needed
          randDoubleVect = V.fromList $ take 20 randDoubleList

          -- choose search type from list with frequencies as input from searchForDuration
          searchBandit = chooseElementAtRandomPair (randDoubleVect V.! 0) thetaList

          -- common build arguments including block and distance
          buildMethod  = chooseElementAtRandomPair (randDoubleVect V.! 10) [("unitary", 0.8), ("block", 0.2)]
          buildType = if searchBandit == "buildCharacter" then "character"
                      else if searchBandit == "buildDistance" then "distance"
                      else chooseElementAtRandomPair (randDoubleVect V.! 11) [("distance", 0.5), ("character", 0.5)]

          numToCharBuild = (10 :: Int)
          numToDistBuild = (1000 :: Int)
          numDistToKeep = (100 :: Int)

          -- to resolve block build graphs
          reconciliationMethod = chooseElementAtRandomPair (randDoubleVect V.! 12) [("eun", 0.5), ("cun", 0.5)]

          wagnerOptions = if buildType == "distance" then
                            if buildMethod == "block" then [("replicates", show numToCharBuild), ("rdwag", ""), ("best", show (1 :: Int))]
                            else  [("replicates", show numToDistBuild), ("rdwag", ""), ("best", show numDistToKeep), ("return", show numToCharBuild)]
                          else if buildType == "character" then
                             if buildMethod == "block" then [("replicates", show (1 :: Int))]
                             else [("replicates", show numToCharBuild)]
                          else []

          blockOptions = if buildMethod == "block" then
                            [("block", ""),("atRandom", ""),("displaytrees", show numToCharBuild),(reconciliationMethod, "")]
                         else []

          -- common swap arguments
          swapKeep = min keepNum (chooseElementAtRandomPair (randDoubleVect V.! 15)  [(1, 0.50), (2, 0.33), (4, 0.17)])

          -- common drift arguments
          maxChanges = chooseElementAtRandomPair (randDoubleVect V.! 1) [("5", 0.33), ("10", 0.34), ("20", 0.33)]
          acceptEqual = chooseElementAtRandomPair (randDoubleVect V.! 2) [("0.1", 0.5), ("0.5", 0.5)]
          acceptWorse = chooseElementAtRandomPair (randDoubleVect V.! 3) [("10.0", 0.33), ("20.0", 0.34), ("40", 0.33)]
          driftArgs = [("drift", ""),("maxChanges", maxChanges), ("acceptEqual", acceptEqual), ("acceptWorse", acceptWorse)]

          -- common annealing arguments
          tempSteps = chooseElementAtRandomPair (randDoubleVect V.! 4)  [("5", 0.33), ("10", 0.34), ("20", 0.33)]
          annealArgs = [("annealing", ""),("steps", tempSteps)]

          -- common fuse options
          fusePairs = chooseElementAtRandomPair (randDoubleVect V.! 5) [("20", 0.45), ("40", 0.45), ("100", 0.1)]
          fuseKeep = 2 * keepNum

          -- network edit options
          netGeneralArgs = [("keep", show keepNum), ("steepest", ""), ("atRandom", ""), ("maxnetedges", show maxNetEdges)]
          netMoveArgs = ("netMove", "") : netGeneralArgs
          netAddArgs = ("netAdd", "") : netGeneralArgs
          netDelArgs = ("netDel", "") : netGeneralArgs
          netAddDelArgs = ("netadddel", "") : netGeneralArgs


          -- Genetic Algorithm Arguments
          -- stops after 2 rounds with no improvement (if generations > 2)
          popSize           = chooseElementAtRandomPair (randDoubleVect V.! 6)  [("10", 0.50), ("20", 0.25), ("40", 0.25)]
          generations       = chooseElementAtRandomPair (randDoubleVect V.! 7)  [("1", 1.0)] -- , "2" , "4"]
          severity          = chooseElementAtRandomPair (randDoubleVect V.! 8)  [("0.0", 0.33), ("1.0", 0.34), ("2.0", 0.33)]
          recombinations    = chooseElementAtRandomPair (randDoubleVect V.! 9)  [("20", 0.45), ("40", 0.45), ("100", 0.1)]

          gaArgs = [("popsize", popSize), ("generations", generations), ("severity", severity), ("recombinations", recombinations), ("stop", "2")]

      in

      -- no input graphs so must build to start
      if null inGraphList' then

         let buildArgs = [(buildType, "")] ++ wagnerOptions ++ blockOptions


             buildGraphs = B.buildGraph buildArgs inGS' inData' pairwiseDistances (randIntList !! 0)
             uniqueGraphs = take keepNum $ GO.selectGraphs Unique (maxBound::Int) 0.0 (-1) buildGraphs

             buildString = if searchBandit `elem` ["buildCharacter", "buildDistance"] then searchBandit
                           else if buildType == "character" then "buildCharacter"
                           else "buildDistance"

             -- string of delta cost of graphs
             deltaString = if null inGraphList' then "10.0"
                           else show ((minimum $ fmap snd6 inGraphList') - (minimum $ fmap snd6 uniqueGraphs))

             currentBestString = show $ minimum $ fmap snd6 uniqueGraphs

             searchString = "," ++ buildString ++ "," ++ deltaString ++ "," ++ currentBestString ++ "," ++ (show $ toSeconds inTime) ++ "," ++ (L.intercalate "," $ fmap showArg buildArgs)


         in  (uniqueGraphs, [searchString ++ thompsonString])


      -- already have some input graphs
      -- choose a method and parameters at random
      -- fuse on single graph will build a couple more first
      else
         let -- choose staticApproximation or not
             -- up top here because used by other non-build options
             -- if happens--need to transform  back before returning
             -- checks if there are non-exact characters to transform
             transformToStaticApproximation = if U.getNumberNonExactCharacters (thd3 inData') == 0 then False
                                              else chooseElementAtRandomPair (randDoubleVect V.! 13) [(True, 0.25), (False, 0.75)]

             -- choose to togle multitraverse (n-speed up compared to default multi-traverse)
             -- if already False--don't change
             transformMultiTraverse = if transformToStaticApproximation then False
                                      else if not (multiTraverseCharacters inGS') then False
                                      else chooseElementAtRandomPair (randDoubleVect V.! 14) [(True, 0.66), (False, 0.34)]

             -- Can't do both static approx and multitraverse:False
             ((inGS, origData, inData, inGraphList), transformString) = if transformToStaticApproximation && (useIA inGS') then
                                                                            (TRANS.transform [("staticapprox",[])] inGS' inData' inData' 0 inGraphList', ",StaticApprox")
                                                                        else if transformMultiTraverse then
                                                                            (TRANS.transform [("multitraverse","false")] inGS' inData' inData' 0 inGraphList', ",MultiTraverse:False")
                                                                        else
                                                                            ((inGS', inData', inData', inGraphList'), "")




         in
         -- bandit list with search arguments set
         -- primes (') for build to start with untransformed data
         let (searchGraphs, searchArgs) = if searchBandit == "buildCharacter" then
                                            let -- build options
                                                buildArgs = [(buildType, "")] ++ wagnerOptions ++ blockOptions
                                            in
                                            -- search
                                            (B.buildGraph buildArgs inGS' inData' pairwiseDistances (randIntList !! 0), buildArgs)

                                          else if searchBandit == "buildDistance" then
                                            let -- build options
                                                buildArgs = [(buildType, "")] ++ wagnerOptions ++ blockOptions
                                            in
                                            -- search for dist builds 1000, keeps 10 best distance then selects 10 best after rediagnosis
                                            -- this line in here to allow for returning lots of rediagnosed distance trees, then
                                            -- reducing to unique best cost trees--but is a memory pig
                                            let graphList = B.buildGraph buildArgs inGS' inData' pairwiseDistances (randIntList !! 0)
                                            in
                                            (take numToCharBuild $ GO.selectGraphs Unique (maxBound::Int) 0.0 (-1) graphList, buildArgs)

                                          else if searchBandit == "buildSPR" then
                                            let -- build part
                                                buildArgs = [(buildType, "")] ++ wagnerOptions ++ blockOptions
                                                buildGraphs = B.buildGraph buildArgs inGS' inData' pairwiseDistances (randIntList !! 0)
                                                buildGraphs' = GO.selectGraphs Unique (maxBound::Int) 0.0 (-1) buildGraphs

                                                -- swap options
                                                swapType = "spr"
                                                swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep), ("atrandom","")]
                                            in
                                            -- search
                                            (R.swapMaster swapArgs inGS inData (randIntList !! 1) buildGraphs', buildArgs ++ swapArgs)

                                          else if searchBandit == "buildAlternate" then
                                            let -- build part
                                                buildArgs = [(buildType, "")] ++ wagnerOptions ++ blockOptions
                                                buildGraphs = B.buildGraph buildArgs inGS' inData' pairwiseDistances (randIntList !! 0)
                                                buildGraphs' = GO.selectGraphs Unique (maxBound::Int) 0.0 (-1) buildGraphs

                                                -- swap options
                                                swapType = "alternate" --default anyway
                                                swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep), ("atrandom","")]
                                            in
                                            -- search
                                            (R.swapMaster swapArgs inGS inData (randIntList !! 1) buildGraphs', buildArgs ++ swapArgs)

                                          else if searchBandit == "swapSPR" then
                                            let -- swap options
                                                swapType = "spr"
                                                swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep), ("atrandom","")]
                                            in
                                            -- search
                                            (R.swapMaster swapArgs inGS inData (randIntList !! 1) inGraphList, swapArgs)

                                          else if searchBandit == "swapAlternate" then
                                            let -- swap options
                                                swapType = "alternate" --default anyway
                                                swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep), ("atrandom","")]
                                            in
                                            -- search
                                            (R.swapMaster swapArgs inGS inData (randIntList !! 1) inGraphList, swapArgs)

                                          else if searchBandit == "driftSPR" then
                                            let -- swap args
                                                swapType = "spr"
                                                swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep), ("atrandom","")]

                                                -- swap with drift (common) arguments
                                                swapDriftArgs = swapArgs ++ driftArgs
                                            in
                                            -- perform search
                                            (R.swapMaster swapDriftArgs inGS inData (randIntList !! 1) inGraphList, swapArgs)

                                          else if searchBandit == "driftAlternate" then
                                            let -- swap args
                                                swapType = "alternate"
                                                swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep), ("atrandom","")]

                                                -- swap with drift (common) arguments
                                                swapDriftArgs = swapArgs ++ driftArgs
                                            in
                                            -- perform search
                                            (R.swapMaster swapDriftArgs inGS inData (randIntList !! 1) inGraphList, swapDriftArgs)

                                          else if searchBandit == "annealSPR" then
                                            let -- swap args
                                                swapType = "spr"
                                                swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep), ("atrandom","")]

                                                -- swap with anneal (common) arguments
                                                swapAnnealArgs = swapArgs ++ annealArgs
                                            in
                                            -- perform search
                                            (R.swapMaster swapAnnealArgs inGS inData (randIntList !! 1) inGraphList, swapAnnealArgs)

                                          else if searchBandit == "annealAlternate" then
                                            let -- swap args
                                                swapType = "alternate"
                                                swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep), ("atrandom","")]

                                                -- swap with anneal (common) arguments
                                                swapAnnealArgs = swapArgs ++ annealArgs
                                            in
                                            -- perform search
                                            (R.swapMaster swapAnnealArgs inGS inData (randIntList !! 1) inGraphList, swapAnnealArgs)

                                          else if searchBandit == "geneticAlgorithm" then
                                            -- args from above
                                            -- perform search
                                            (R.geneticAlgorithmMaster gaArgs inGS inData (randIntList !! 1) inGraphList, gaArgs)

                                          else if searchBandit == "fuse" then
                                            -- should more graphs be added if only one?  Would downweight fuse perhpas too much
                                            let -- fuse arguments
                                                fuseArgs = [("none",""), ("all",""), ("unique",""), ("atrandom", ""), ("pairs", fusePairs), ("keep", show fuseKeep)]
                                            in
                                            -- perform search
                                            (R.fuseGraphs fuseArgs inGS inData(randIntList !! 1) inGraphList, fuseArgs)

                                          else if searchBandit == "fuseSPR" then
                                            let -- fuse arguments
                                                fuseArgs = [("spr", ""), ("all",""), ("unique",""), ("atrandom", ""), ("pairs", fusePairs), ("keep", show fuseKeep)]
                                            in
                                            -- perform search
                                            (R.fuseGraphs fuseArgs inGS inData (randIntList !! 1) inGraphList, fuseArgs)

                                          else if searchBandit == "fuseTBR" then
                                            let -- fuse arguments
                                                fuseArgs = [("tbr", ""), ("all",""), ("unique",""), ("atrandom", ""), ("pairs", fusePairs), ("keep", show fuseKeep)]
                                            in
                                            -- perform search
                                            (R.fuseGraphs fuseArgs inGS inData (randIntList !! 1) inGraphList, fuseArgs)

                                          else if searchBandit == "networkAdd" then
                                            let -- network add args
                                                netEditArgs = netAddArgs
                                            in
                                            -- perform search
                                            (R.netEdgeMaster netEditArgs inGS inData (randIntList !! 1) inGraphList, netEditArgs)

                                          else if searchBandit == "networkDelete" then
                                            let -- network delete args
                                                netEditArgs = netDelArgs
                                            in
                                            -- perform search
                                            (R.netEdgeMaster netEditArgs inGS inData (randIntList !! 1) inGraphList, netEditArgs)

                                          else if searchBandit == "networkAddDelete" then
                                            let -- network add/delete args
                                                netEditArgs = netAddDelArgs
                                            in
                                            -- perform search
                                            (R.netEdgeMaster netEditArgs inGS inData (randIntList !! 1) inGraphList, netEditArgs)

                                          else if searchBandit == "networkAddDelete" then
                                            let -- network add/delete args
                                                netEditArgs = netAddDelArgs
                                            in
                                            -- perform search
                                            (R.netEdgeMaster netEditArgs inGS inData (randIntList !! 1) inGraphList, netEditArgs)

                                          else if searchBandit == "networkMove" then
                                            let -- network move args
                                                netEditArgs = netMoveArgs
                                            in
                                            -- perform search
                                            (R.netEdgeMaster netEditArgs inGS inData (randIntList !! 1) inGraphList, netEditArgs)

                                          else if searchBandit == "driftNetwork" then
                                            let -- network add/delete  + drift args
                                                netEditArgs = netAddDelArgs ++ driftArgs
                                            in
                                            -- perform search
                                            (R.netEdgeMaster netEditArgs inGS inData (randIntList !! 1) inGraphList, netEditArgs)

                                          else if searchBandit == "annealNetwork" then
                                            let -- network add/delete  + annealing  args
                                                netEditArgs = netAddDelArgs ++ annealArgs
                                            in
                                            -- perform search
                                            (R.netEdgeMaster netEditArgs inGS inData (randIntList !! 1) inGraphList, netEditArgs)


                                          else error ("Unknown/unimplemented method in search: " ++ searchBandit)

             -- process
             uniqueGraphs' = take keepNum $ GO.selectGraphs Unique (maxBound::Int) 0.0 (-1) (searchGraphs ++ inGraphList)
             (uniqueGraphs, transString) = if (not transformToStaticApproximation && not transformMultiTraverse) then (uniqueGraphs', "")
                                           else if transformToStaticApproximation then
                                                (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', ",Dynamic")
                                           else
                                                (fth4 $ TRANS.transform [("multiTraverse","true")] inGS origData inData 0 uniqueGraphs', ",MultiTraverse:True")

             -- string of delta and cost of graphs
             deltaString = if null inGraphList' then "10.0,"
                           else show ((minimum $ fmap snd6 inGraphList') - (minimum $ fmap snd6 uniqueGraphs)) -- ++ ","

             currentBestString = show $ minimum $ fmap snd6 uniqueGraphs

             -- create string for search stats
             searchString = "," ++ searchBandit ++ "," ++ deltaString ++ "," ++ currentBestString ++ "," ++ (show $ toSeconds inTime) ++ "," ++ (L.intercalate "," $ fmap showArg searchArgs) ++ transformString ++ transString

         in
         (uniqueGraphs, [searchString ++ thompsonString])
         where showArg a = (fst a) ++ ":" ++ (snd a)


-- | getSearchParams takes arguments and returns search params
getSearchParams :: [Argument] -> (Int, Int, Int, Bool, Int, String, Int, Int)
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
            | null minutesList = Just 1
            | otherwise = readMaybe (snd $ head minutesList) :: Maybe Int

          secondsList = filter ((=="seconds").fst) lcArgList
          seconds
            | length secondsList > 1 =
              errorWithoutStackTrace ("Multiple 'seconds' number specifications in search command--can have only one: " ++ show inArgs)
            | null secondsList = Just 0
            | otherwise = readMaybe (snd $ head secondsList) :: Maybe Int

          maxNetEdgesList = filter ((=="maxnetedges").fst) lcArgList
          maxNetEdges
              | length maxNetEdgesList > 1 =
                errorWithoutStackTrace ("Multiple 'maxNetEdges' number specifications in netEdge command--can have only one: " ++ show inArgs)
              | null maxNetEdgesList = Just 10
              | otherwise = readMaybe (snd $ head maxNetEdgesList) :: Maybe Int

          thompsonList = filter ((=="thompson").fst) lcArgList
          mFactor
            | length thompsonList > 1 =
              errorWithoutStackTrace ("Multiple 'Thompson' number specifications in search command--can have only one: " ++ show inArgs)
            | null thompsonList = Just 1
            | otherwise = readMaybe (snd $ head thompsonList) :: Maybe Int

          stopList = filter ((=="stop").fst) lcArgList
          stopNum
            | length stopList > 1 =
              errorWithoutStackTrace ("Multiple 'stop' number specifications in search command--can have only one: " ++ show inArgs)
            | null stopList = Just (maxBound :: Int)
            | otherwise = readMaybe (snd $ head stopList) :: Maybe Int

          thompson = any ((=="thompson").fst) lcArgList
          mLinear = any ((=="linear").fst) lcArgList
          mExponential = any ((=="exponential").fst) lcArgList
          mSimple = any ((=="simple").fst) lcArgList

          mFunction = if mLinear && mExponential then
                            trace ("Thompson recency function specification has both 'linear' and 'exponential', defaulting to 'linear'")
                            "linear"
                      else if mLinear then "linear"
                      else if mExponential then "exponential"
                      else if mSimple then "simple"
                      else "linear"

      in
      if isNothing keepNum then errorWithoutStackTrace ("Keep specification not an integer in search: "  ++ show (head keepList))
      else if isNothing instances then errorWithoutStackTrace ("Instances specification not an integer in search: "  ++ show (head instancesList))
      else if isNothing days then errorWithoutStackTrace ("Days specification not an integer in search: "  ++ show (head daysList))
      else if isNothing hours then errorWithoutStackTrace ("Hours specification not an integer in search: "  ++ show (head hoursList))
      else if isNothing minutes then errorWithoutStackTrace ("Minutes specification not an integer in search: "  ++ show (head minutesList))
      else if isNothing seconds then errorWithoutStackTrace ("Seconds factor specification not an integer in search: "  ++ show (head secondsList))
      else if isNothing mFactor then errorWithoutStackTrace ("Thompson mFactor specification not an integer or not found in search (e.g. Thompson:1) "  ++ show (head thompsonList))
      else if isNothing maxNetEdges then errorWithoutStackTrace ("Search 'maxNetEdges' specification not an integer or not found (e.g. maxNetEdges:8): "  ++ show (snd $ head maxNetEdgesList))
      else if isNothing stopNum then errorWithoutStackTrace ("Search stop specification not an integer or not found in search (e.g. stop:10) "  ++ show (head stopList))

      else
         let seconds' = if ((fromJust minutes > 0) || (fromJust hours > 0) || (fromJust days > 0)) && (null secondsList) then Just 0
                        else seconds
             searchTime = (fromJust seconds') + (60 * (fromJust minutes)) + (3600 * (fromJust hours))
         in
         (searchTime, fromJust keepNum, fromJust instances, thompson, fromJust mFactor, mFunction, fromJust maxNetEdges, fromJust stopNum)



