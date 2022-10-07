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
-- import qualified GraphOptimization.Traversals as T
import qualified Graphs.GraphOperations       as GO
import qualified Search.Build                 as B
import qualified Search.Refinement            as R
import           System.Timing
import           Text.Read
import           Types.Types
import           Debug.Trace
import           System.Random
import qualified Data.Vector                  as V

-- | treeBanditList is list of search types to be chosen from if graphType is tree
treeBanditList :: [String]
treeBanditList = ["buildCharacter", "buildDistance", "swapSPR", "swapAlternate", "driftSPR", 
                 "driftAlternate", "annealSPR", "annealAlternate", "geneticAlgorithm", 
                 "fuse", "fuseSPR", "fuseAlternate"] 

-- | netWorkBanditList is list of search types unique to graphType network
netWorkBanditList :: [String]
netWorkBanditList = ["networkAdd", "networkDelete", "networkAddDelete", "driftNetwork", "annealNetwork"] -- "networkMove"

-- | fullBanditList is list of search types to be chosen from if Network
fullBanditList :: [String]
fullBanditList = treeBanditList ++ netWorkBanditList

-- | A strict, three-way version of 'uncurry'.
uncurry3' :: (Functor f, NFData d) => (a -> b -> c -> f d) -> (a, b, c) -> f d
uncurry3' f (a, b, c) = force <$> f a b c

-- | search timed randomized search returns graph list and comment list with info String for each search instance
search :: [Argument] -> GlobalSettings -> ProcessedData -> [[VertexCost]] -> Int -> [PhylogeneticGraph] -> IO ([PhylogeneticGraph], [[String]])
search inArgs inGS inData pairwiseDistances rSeed inGraphList =
   let (searchTime, keepNum, instances, thompsonSample, mFactor, mFunction, maxNetEdges) = getSearchParams inArgs
       
       -- flatThetaList is the initial prior list (flat) of search (bandit) choices
       -- can also be used in search for non-Thomspon search
       flatThetaList = if graphType inGS == Tree then 
                            zip treeBanditList (L.replicate (length treeBanditList) (1.0 / (fromIntegral $ length treeBanditList)))
                       else 
                            zip fullBanditList (L.replicate (length fullBanditList) (1.0 / (fromIntegral $ length fullBanditList)))
       
       threshold   = fromSeconds . fromIntegral $ (95 * searchTime) `div` 100
       searchTimed = uncurry3' $ searchForDuration inGS inData pairwiseDistances keepNum thompsonSample mFactor mFunction flatThetaList maxNetEdges 0 threshold []
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
                  -> [String] 
                  -> Int 
                  -> [Int] 
                  -> ([PhylogeneticGraph], [String]) 
                  -> IO ([PhylogeneticGraph], [String])
searchForDuration inGS inData pairwiseDistances keepNum thompsonSample mFactor mFunction thetaList counter maxNetEdges allotedSeconds inCommentList refIndex seedList input@(inGraphList, infoStringList) = do
   (elapsedSeconds, output) <- timeOpUT $
       let result = force $ performSearch' inGS inData pairwiseDistances keepNum thompsonSample thetaList maxNetEdges (head seedList) input
       in  pure result

   -- update theta list based on performance    
   let updatedThetaList = updateTheta mFactor mFunction counter (snd output) thetaList 

   let remainingTime = allotedSeconds `timeDifference` elapsedSeconds
   putStrLn $ unlines [ "Thread   \t" <> show refIndex
                      , "Alloted  \t" <> show allotedSeconds
                      , "Ellapsed \t" <> show elapsedSeconds
                      , "Remaining\t" <> show remainingTime
                      ]
   if elapsedSeconds >= allotedSeconds
   then pure output
   else searchForDuration inGS inData pairwiseDistances keepNum thompsonSample mFactor mFunction updatedThetaList (counter + 1) maxNetEdges remainingTime (inCommentList ++ (snd output)) refIndex (tail seedList) $ bimap (inGraphList <>) (infoStringList <>) output

-- | updateTheta updates the expected success parameters for the bandit search list
updateTheta :: Int -> String -> Int -> [String] -> [(String, Double)] -> [(String, Double)]
updateTheta  mFactor mFunction counter infoStringList inPairList =
    if null inPairList then []
    else 
        inPairList


-- | performSearch' takes in put graphs and performs randomized build and search with time limit
-- Thompson sampling and mFactor to pick strategy from updated theta success values
-- the random calls return the tail of the input list to avoid long list access--can do vector since infinite
performSearch' :: GlobalSettings 
               -> ProcessedData 
               -> [[VertexCost]] 
               -> Int 
               -> Bool 
               -> [(String, Double)] 
               -> Int
               -> Int 
               -> ([PhylogeneticGraph], [String]) 
               -> ([PhylogeneticGraph], [String])
performSearch' inGS' inData' pairwiseDistances keepNum thompsonSample thetaList maxNetEdges rSeed (inGraphList', _) =
      -- set up basic parameters for search/refine methods
      let -- set up log for sample
          thompsonString = if not thompsonSample then ","
                           else "," ++ (show thetaList) 

          -- get infinite lists if integers and doubles
          randIntList = randomIntList rSeed
          randDoubleList = randoms (mkStdGen rSeed) :: [Double]

          -- this for constant access to random doubles need take for infinite list
          -- need to update as more random doubles are needed
          randDoubleVect = V.fromList $ take 20 randDoubleList
          
          -- choose search type from list with frequencies as input from searchForDuration
          searchBandit = chooseElementAtRandomPair (randDoubleVect V.! 0) thetaList

          -- common swap arguments
          swapKeep = keepNum

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

          -- network edit options 'move' disabled for now
          netGeneralArgs = [("keep", show keepNum), ("steepest", ""), ("atRandom", ""), ("maxnetedges", show maxNetEdges)]
          -- netMoveArgs = ("netMove", "") : netGeneralArgs until fix network move
          netAddArgs = ("netAdd", "") : netGeneralArgs
          netDelArgs = ("netDel", "") : netGeneralArgs
          netAddDelArgs = ("netadddel", "") : netGeneralArgs


          -- Genetic Algorithm Arguments
          popSize           = chooseElementAtRandomPair (randDoubleVect V.! 6)  [("10", 0.50), ("20", 0.25), ("40", 0.25)]
          generations       = chooseElementAtRandomPair (randDoubleVect V.! 7)  [("1", 1.0)] -- , "2" , "4"]
          severity          = chooseElementAtRandomPair (randDoubleVect V.! 8)  [("0.0", 0.33), ("1.0", 0.34), ("2.0", 0.33)]
          recombinations    = chooseElementAtRandomPair (randDoubleVect V.! 9)  [("20", 0.45), ("40", 0.45), ("100", 0.1)]

          gaArgs = [("popsize", popSize), ("generations", generations), ("severity", severity), ("recombinations", recombinations)]
         
      in

      -- no input graphs so must build to start--or if buildX chosen as searchBandit 
      if null inGraphList' || searchBandit `elem` ["buildCharacter", "buildDistance"] then

         let -- build options including block and distance
             -- primes (') in data and graphlist since not reseat by potnatila statric apporx transformation
             buildMethod  = chooseElementAtRandomPair (randDoubleVect V.! 10) [("unitary", 0.9), ("block", 0.1)]
             buildType = if searchBandit == "buildCharacter" then "character"
                         else if searchBandit == "buildDistance" then "distance"
                         else 
                            chooseElementAtRandomPair (randDoubleVect V.! 11) [("distance", 0.5), ("character", 0.5)]
                            
             numToCharBuild = (10 :: Int)
             numToDistBuild = (100 :: Int)

             -- tail here in case both called to increment random value
             reconciliationMethod = chooseElementAtRandomPair (randDoubleVect V.! 12) [("eun", 0.5), ("cun", 0.5)]

             wagnerOptions = if buildType == "distance" then
                                if buildMethod == "block" then [("replicates", show numToDistBuild), ("rdwag", ""), ("best", show (1 :: Int))]
                                else  [("replicates", show numToDistBuild), ("rdwag", ""), ("best", show numToCharBuild)]
                             else if buildType == "character" then
                                 if buildMethod == "block" then [("replicates", show (1 :: Int))]
                                 else [("replicates", show numToCharBuild)] 
                             else []

             blockOptions = if buildMethod == "block" then 
                                [("block", ""),("atRandom", ""),("displaytrees", show numToCharBuild),(reconciliationMethod, "")]
                            else []

             buildArgs = [(buildType, "")] ++ wagnerOptions ++ blockOptions


             buildGraphs = B.buildGraph buildArgs inGS' inData' pairwiseDistances (randIntList !! 0)
             uniqueBuildGraphs = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] buildGraphs

             searchString = "Build " ++ (L.intercalate "," $ fmap showArg buildArgs)
                            
         in  (uniqueBuildGraphs, [searchString ++ thompsonString])


      -- already have some input graphs
      -- choose a method and parameters at random
      -- fuse on single graph will build a couple more first
      else
         let -- choose staticApproximation or not
             -- up top here because used by other non-build options
             -- if happens--need to rtansfomr back before returning
             transformToStaticApproximation = chooseElementAtRandomPair (randDoubleVect V.! 13) [(True, 0.67), (False, 0.33)]
             ((inGS, origData, inData, inGraphList), staticApproxString) = if transformToStaticApproximation then
                                                                            (TRANS.transform [("staticapprox",[])] inGS' inData' inData' 0 inGraphList', "StaticApprox ")
                                                                           else ((inGS', inData', inData', inGraphList'), "")
         in
         if searchBandit == "swapSPR" then 
            let -- swap options
                swapType = "spr"
                swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep)]
                
                -- search
                searchGraphs = R.swapMaster swapArgs inGS inData (head randIntList)inGraphList

                -- process
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (searchGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")

                -- create string for search stats
                searchString = staticApproxString ++ searchBandit ++ (L.intercalate "," $ fmap showArg swapArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])


         else if searchBandit == "swapAlternate" then 
            let -- swap options
                swapType = "alternate" --default anyway
                swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep)]
                
                -- search
                searchGraphs = R.swapMaster swapArgs inGS inData (head randIntList)inGraphList
                
                -- process
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (searchGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")
                
                -- create string for search stats
                searchString = staticApproxString ++ searchBandit ++ (L.intercalate "," $ fmap showArg swapArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])

         else if searchBandit == "driftSPR" then
            let 
                -- swap args
                swapType = "spr"
                swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep)]
                
                -- swap with drift (common) arguments
                swapDriftArgs = swapArgs ++ driftArgs

                -- perform search
                searchGraphs = R.swapMaster swapDriftArgs inGS inData (head randIntList)inGraphList
                
                -- process
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (searchGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")

                -- create string for search stats
                searchString = staticApproxString ++ searchBandit ++ (L.intercalate "," $ fmap showArg  swapDriftArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])
            
         else if searchBandit == "driftAlternate" then
            let 
                -- swap args
                swapType = "alternate"
                swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep)]
                
                -- swap with drift (common) arguments
                swapDriftArgs = swapArgs ++ driftArgs

                -- perform search
                searchGraphs = R.swapMaster swapDriftArgs inGS inData (head randIntList)inGraphList
                
                -- process
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (searchGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")

                -- create string for search stats
                searchString = staticApproxString ++ searchBandit ++ (L.intercalate "," $ fmap showArg  swapDriftArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])

            
         else if searchBandit == "annealSPR" then
            let 
                -- swap args
                swapType = "spr"
                swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep)]
                
                -- swap with anneal (common) arguments
                swapAnnealArgs = swapArgs ++ annealArgs

                -- perform search
                searchGraphs = R.swapMaster swapAnnealArgs inGS inData (head randIntList)inGraphList
                
                -- process
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (searchGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")

                -- create string for search stats
                searchString = staticApproxString ++ searchBandit ++ (L.intercalate "," $ fmap showArg  swapAnnealArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])
            
         else if searchBandit == "annealAlternate" then
            let 
                -- swap args
                swapType = "alternate"
                swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep)]
                
                -- swap with anneal (common) arguments
                swapAnnealArgs = swapArgs ++ annealArgs

                -- perform search
                searchGraphs = R.swapMaster swapAnnealArgs inGS inData (head randIntList)inGraphList
                
                -- process
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (searchGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")

                -- create string for search stats
                searchString = staticApproxString ++ searchBandit ++ (L.intercalate "," $ fmap showArg  swapAnnealArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])

            
         else if searchBandit == "geneticAlgorithm" then
            let -- args from above
                -- perform search
                searchGraphs = R.geneticAlgorithmMaster gaArgs inGS inData (head randIntList)inGraphList
                
                -- process
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (searchGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")

                -- create string for search stats
                searchString = staticApproxString ++ searchBandit ++ (L.intercalate "," $ fmap showArg  gaArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])

         else if searchBandit == "fuse" then
            -- should more graphs be added if only one?  Would downweight fuse perhpas too much
            let -- fuse arguments
                fuseArgs = [("steepest",""), ("unique",""), ("atrandom", ""), ("pairs", fusePairs), ("keep", show fuseKeep)]

                -- perform search
                searchGraphs = R.fuseGraphs fuseArgs inGS inData (head randIntList)inGraphList
                
                -- process
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (searchGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")
                
                -- create string for search stats
                searchString = staticApproxString ++ searchBandit ++ (L.intercalate "," $ fmap showArg  fuseArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])

         else if searchBandit == "fuseSPR" then
            let -- fuse arguments
                fuseArgs = [("spr", ""), ("steepest",""), ("unique",""), ("atrandom", ""), ("pairs", fusePairs), ("keep", show fuseKeep)]

                -- perform search
                searchGraphs = R.fuseGraphs fuseArgs inGS inData (head randIntList)inGraphList
                
                -- process
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (searchGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")
                
                -- create string for search stats
                searchString = staticApproxString ++ searchBandit ++ (L.intercalate "," $ fmap showArg  fuseArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])
            
         else if searchBandit == "fuseAlternate" then
            let -- fuse arguments
                fuseArgs = [("alternate", ""), ("steepest",""), ("unique",""), ("atrandom", ""), ("pairs", fusePairs), ("keep", show fuseKeep)]

                -- perform search
                searchGraphs = R.fuseGraphs fuseArgs inGS inData (head randIntList)inGraphList
                
                -- process
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (searchGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")
                
                -- create string for search stats
                searchString = staticApproxString ++ searchBandit ++ (L.intercalate "," $ fmap showArg  fuseArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])
            
         else if searchBandit == "networkAdd" then
            let -- network add args
                netEditArgs = netAddArgs

                -- perform search
                searchGraphs = R.netEdgeMaster netEditArgs inGS inData (head randIntList) inGraphList

                -- process
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (searchGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")

                -- create string for search stats
                searchString = staticApproxString ++ searchBandit ++ (L.intercalate "," $ fmap showArg  netEditArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])
            
         else if searchBandit == "networkDelete" then
            let -- network delete args
                netEditArgs = netDelArgs

                -- perform search
                searchGraphs = R.netEdgeMaster netEditArgs inGS inData (head randIntList) inGraphList

                -- process
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (searchGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")

                -- create string for search stats
                searchString = staticApproxString ++ searchBandit ++ (L.intercalate "," $ fmap showArg  netEditArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])
            
         else if searchBandit == "networkAddDelete" then
            let -- network add/delete args
                netEditArgs = netAddDelArgs

                -- perform search
                searchGraphs = R.netEdgeMaster netEditArgs inGS inData (head randIntList) inGraphList

                -- process
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (searchGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")

                -- create string for search stats
                searchString = staticApproxString ++ searchBandit ++ (L.intercalate "," $ fmap showArg  netEditArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])
            
         {-Inactive till fix net move
         else if searchBandit == "networkMove" then
            let -- network move args
                netEditArgs = netMoveArgs

                -- perform search
                searchGraphs = R.netEdgeMaster netEditArgs inGS inData (head randIntList) inGraphList

                -- process
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (searchGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")

                -- create string for search stats
                searchString = staticApproxString ++ searchBandit ++ (L.intercalate "," $ fmap showArg  netEditArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])
            -}
            
         else if searchBandit == "driftNetwork" then
            let -- network add/delete  + drift args
                netEditArgs = netAddDelArgs ++ driftArgs

                -- perform search
                searchGraphs = R.netEdgeMaster netEditArgs inGS inData (head randIntList) inGraphList

                -- process
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (searchGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")

                -- create string for search stats
                searchString = staticApproxString ++ searchBandit ++ (L.intercalate "," $ fmap showArg  netEditArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])

         else if searchBandit == "annealNetwork" then
            let -- network add/delete  + annealing  args
                netEditArgs = netAddDelArgs ++ annealArgs

                -- perform search
                searchGraphs = R.netEdgeMaster netEditArgs inGS inData (head randIntList) inGraphList

                -- process
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (searchGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")

                -- create string for search stats
                searchString = staticApproxString ++ searchBandit ++ (L.intercalate "," $ fmap showArg  netEditArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])

         else error ("Unknown/unimplemented method in search: " ++ searchBandit)
      where showArg a = "(" ++ (fst a) ++ "," ++ (snd a) ++ ")"



-- | perform search takes in put graphs and performs randomized build and search with time limit
-- Thompson sampling and mFactor to pick strategy from updated theta success values
performSearch :: GlobalSettings -> ProcessedData -> [[VertexCost]] -> Int -> Bool -> [Double] -> Int -> ([PhylogeneticGraph], [String]) -> ([PhylogeneticGraph], [String])
performSearch inGS' inData' pairwiseDistances keepNum thompsonSample thetaList rSeed (inGraphList', _) =
      -- set up basic parameters for search/refine methods
      let thompsonString = if not thompsonSample then ","
                           else if graphType inGS == Tree then 
                                "," ++ (show (zip treeBanditList thetaList)) 
                           else 
                                "," ++ (show (zip fullBanditList thetaList)) 

          randIntList = randomIntList rSeed
          buildType = getRandomElement (head randIntList)["distance", "character"]
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
          acceptWorse = getRandomElement (randIntList !! 17) ["10.0", "20.0", "40"]

          driftArgs = [("drift", ""),("maxChanges", maxChanges), ("acceptEqual", acceptEqual), ("acceptWorse", acceptWorse)]

          -- simulated annealing options
          tempSteps = getRandomElement (randIntList !! 18) ["5", "10", "15"]

          simulatedAnnealArgs = [("annealing", ""),("steps", tempSteps)]

          -- swap with drift arguments
          swapDriftArgs = swapArgs ++ driftArgs

          -- swap with simulated anneling options
          swapAnnealArgs = swapArgs ++ simulatedAnnealArgs

          -- choose staticApproximation or not
          transformToStaticApproximation = (not . null) inGraphList' && getRandomElement (randIntList !! 19) [True, False, False]
          ((inGS, origData, inData, inGraphList), staticApproxString) = if transformToStaticApproximation then
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
         in  (uniqueSwapGraphs, [searchString ++ thompsonString])
         -- performSearch inArgs inGS inData pairwiseDistances keepNum startTime searchTime (searchString : infoStringList) (randIntList !! 6) uniqueSwapGraphs


      -- already have some input graphs
      -- choose a method and paramteres at random
      else
         let operation = if (graphType inGS == Tree) then getRandomElement (randIntList !! 7) ["buildSwap","fuse", "GeneticAlgorithm", "swapAnneal", "swapDrift"]
                         else getRandomElement (randIntList !! 7) ["buildSwap","fuse", "GeneticAlgorithm", "swapAnneal", "swapDrift", "netAdd", "netDelete"] -- , "netMove"] -- add/del/move edges with and without drifting

             -- this so 1/2 time annealing
             saDrift = getRandomElement (randIntList !! 19) ["noSA", "noSA", "drift", "anneal"]

             {-
             -- for rediagnosis after static approx
             -- sequential rediagnosis since assuming that if in parallel teh search operations are running in parallel
             pruneEdges = False
             warnPruneEdges = False
             startVertex = Nothing
             -}
         in
         if operation == "buildSwap" then
            let buildGraphs = B.buildGraph buildArgs inGS inData pairwiseDistances (randIntList !! 4)
                uniqueBuildGraphs = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] buildGraphs
                swapGraphs = R.swapMaster swapArgs inGS inData (randIntList !! 5) uniqueBuildGraphs
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (swapGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")
                                              -- else fmap (T.multiTraverseFullyLabelGraph inGS' inData' pruneEdges warnPruneEdges startVertex) (fmap fst6 uniqueGraphs')
                searchString = staticApproxString ++ "Build " ++ (L.intercalate "," $ fmap showArg buildArgs) ++ " Swap " ++ (L.intercalate "," $ fmap showArg  swapArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])

         else if operation == "fuse" then
            let fuseGraphs = R.fuseGraphs fuseArgs inGS inData (randIntList !! 10) inGraphList
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (fuseGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")
                                              -- uniqueGraphs = if not transformToStaticApproximation then uniqueGraphs'
                                              -- else fmap (T.multiTraverseFullyLabelGraph inGS' inData' pruneEdges warnPruneEdges startVertex) (fmap fst6 uniqueGraphs')
                searchString = staticApproxString ++ "Fuse " ++ (L.intercalate "," $ fmap showArg  fuseArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])

         else if operation == "GeneticAlgorithm" then
            let gaGraphs = R.geneticAlgorithmMaster gaArgs inGS inData (randIntList !! 10) inGraphList
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (gaGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")
                                              -- if not transformToStaticApproximation then uniqueGraphs'
                                              -- else fmap (T.multiTraverseFullyLabelGraph inGS' inData' pruneEdges warnPruneEdges startVertex) (fmap fst6 uniqueGraphs')
                searchString = staticApproxString ++ "Genetic Algorithm " ++ (L.intercalate "," $ fmap showArg  gaArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])

         else if operation == "swapDrift" then
            let swapDriftGraphs = R.swapMaster swapDriftArgs inGS inData (randIntList !! 10) inGraphList
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (swapDriftGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")
                                              -- if not transformToStaticApproximation then uniqueGraphs'
                                              -- else fmap (T.multiTraverseFullyLabelGraph inGS' inData' pruneEdges warnPruneEdges startVertex) (fmap fst6 uniqueGraphs')
                searchString = staticApproxString ++ "SwapDrift " ++ (L.intercalate "," $ fmap showArg  swapDriftArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])

         else if operation == "swapAnneal" then
            let swapAnnealGraphs = R.swapMaster swapAnnealArgs inGS inData (randIntList !! 10) inGraphList
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (swapAnnealGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")
                                              -- if not transformToStaticApproximation then uniqueGraphs'
                                              -- else fmap (T.multiTraverseFullyLabelGraph inGS' inData' pruneEdges warnPruneEdges startVertex) (fmap fst6 uniqueGraphs')
                searchString = staticApproxString ++ "SwapAnneal " ++ (L.intercalate "," $ fmap showArg  swapAnnealArgs) ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])

         else if operation == "netMove" then
            let netMoveArgs' = if saDrift == "noSA" then netMoveArgs
                               else if saDrift == "drift" then netMoveArgs ++ driftArgs
                               else netMoveArgs ++ simulatedAnnealArgs
                netMoveGraphs = R.netEdgeMaster netMoveArgs' inGS inData (randIntList !! 10) inGraphList
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (netMoveGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")
                                              -- if not transformToStaticApproximation then uniqueGraphs'
                                              -- else fmap (T.multiTraverseFullyLabelGraph inGS' inData' pruneEdges warnPruneEdges startVertex) (fmap fst6 uniqueGraphs')
                searchString = staticApproxString ++ "NetMove " ++ (L.intercalate "," $ fmap showArg  netMoveArgs') ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])

         else if operation == "netAdd" then
            let netAddArgs' = if saDrift == "noSA" then netAddArgs
                              else if saDrift == "drift" then netAddArgs ++ driftArgs
                              else netAddArgs ++ simulatedAnnealArgs
                netAddGraphs = R.netEdgeMaster netAddArgs' inGS inData (randIntList !! 10) inGraphList
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (netAddGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")
                                              -- if not transformToStaticApproximation then uniqueGraphs'
                                              -- else fmap (T.multiTraverseFullyLabelGraph inGS' inData' pruneEdges warnPruneEdges startVertex) (fmap fst6 uniqueGraphs')
                searchString = staticApproxString ++ "NetAdd " ++ (L.intercalate "," $ fmap showArg  netAddArgs') ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])

         else if operation == "netDelete" then
            let netDelArgs' = if saDrift == "noSA" then netDelArgs
                              else if saDrift == "drift" then netDelArgs ++ driftArgs
                              else netDelArgs ++ simulatedAnnealArgs
                netDelGraphs = R.netEdgeMaster netDelArgs' inGS inData (randIntList !! 10) inGraphList
                uniqueGraphs' = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (netDelGraphs ++ inGraphList)
                (uniqueGraphs, transString) = if not transformToStaticApproximation then (uniqueGraphs', "")
                                              else (fth4 $ TRANS.transform [("dynamic",[])] inGS' origData inData 0 uniqueGraphs', " Dynamic, ")
                                              -- if not transformToStaticApproximation then uniqueGraphs'
                                              -- else fmap (T.multiTraverseFullyLabelGraph inGS' inData' pruneEdges warnPruneEdges startVertex) (fmap fst6 uniqueGraphs')
                searchString = staticApproxString ++ "netDelete " ++ (L.intercalate "," $ fmap showArg  netDelArgs') ++ transString
            in
            (uniqueGraphs, [searchString ++ thompsonString])


         else error ("Unknown/unimplemented method in search: " ++ operation)
      where showArg a = "(" ++ (fst a) ++ "," ++ (snd a) ++ ")"

-- | getSearchParams takes arguments and returns search params
getSearchParams :: [Argument] -> (Int, Int, Int, Bool, Int, String, Int)
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

          thompson = any ((=="thompson").fst) lcArgList
          mLinear = any ((=="linear").fst) lcArgList
          mExponential = any ((=="exponential").fst) lcArgList

          mFunction = if mLinear && mExponential then
                            trace ("Thompson recency function specification has both 'linear' and 'exponential', defaulting to 'linear'")
                            "linear"
                      else if mLinear then "linear"
                      else if mExponential then "exponential"
                      else "linear"

      in
      if isNothing keepNum then errorWithoutStackTrace ("Keep specification not an integer in search: "  ++ show (head keepList))
      else if isNothing instances then errorWithoutStackTrace ("Instnaces specification not an integer in search: "  ++ show (head instancesList))
      else if isNothing days then errorWithoutStackTrace ("Days specification not an integer in search: "  ++ show (head daysList))
      else if isNothing hours then errorWithoutStackTrace ("Hours specification not an integer in search: "  ++ show (head hoursList))
      else if isNothing minutes then errorWithoutStackTrace ("Minutes specification not an integer in search: "  ++ show (head minutesList))
      else if isNothing seconds then errorWithoutStackTrace ("seconds factor specification not an integer in search: "  ++ show (head secondsList))
      else if isNothing mFactor then errorWithoutStackTrace ("Thompson mFactor specification not an integer in search: "  ++ show (head secondsList))
      else if isNothing maxNetEdges then errorWithoutStackTrace ("Search 'maxNetEdges' specification not an integer (e.g. maxNetEdges:8): "  ++ show (snd $ head maxNetEdgesList))
         
      else
         let seconds' = if ((fromJust minutes > 0) || (fromJust hours > 0) || (fromJust days > 0)) && (null secondsList) then Just 0
                        else seconds
             searchTime = (fromJust seconds') + (60 * (fromJust minutes)) + (3600 * (fromJust hours))
         in
         (searchTime, fromJust keepNum, fromJust instances, thompson, fromJust mFactor, mFunction, fromJust maxNetEdges)

