{- |
Module      :  NetworkAddDelete.hs
Description :  Module specifying graph egde adding and deleting functions
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

module Search.NetworkAddDelete  ( deleteAllNetEdges
                                , insertAllNetEdges
                                , moveAllNetEdges
                                , deltaPenaltyAdjustment
                                , deleteNetEdge
                                , deleteOneNetAddAll
                                , addDeleteNetEdges
                                , getCharacterDelta
                                , getBlockDelta
                                -- these are not used but to quiet warnings
                                , heuristicDeleteDelta
                                , heuristicAddDelta
                                , heuristicAddDelta'
                                ) where

import           Data.Bits
import qualified Data.InfList                                     as IL
import           Data.Maybe
import qualified Data.Text.Lazy                                   as TL
import qualified Data.Vector                                      as V
import           Debug.Trace
import           GeneralUtilities
import qualified GraphOptimization.Medians                        as M
import qualified GraphOptimization.PostOrderSoftWiredFunctions    as POSW
import qualified GraphOptimization.PostOrderSoftWiredFunctionsNew as NEW
import qualified GraphOptimization.PreOrderFunctions              as PRE
import qualified GraphOptimization.Traversals                     as T
import qualified Graphs.GraphOperations                           as GO
import qualified ParallelUtilities                                as PU
import           Types.Types
import qualified Utilities.LocalGraph                             as LG
import qualified Utilities.Utilities                              as U
-- import qualified Data.List                            as L


-- | addDeleteNetEdges is a wrapper for addDeleteNetEdges' allowing for multiple simulated annealing rounds
addDeleteNetEdges :: GlobalSettings
                  -> ProcessedData
                  -> Int
                  -> Int
                  -> Int
                  -> Int
                  -> Int
                  -> Bool
                  -> Bool
                  -> Bool
                  -> ([PhylogeneticGraph], VertexCost)
                  -> (Maybe SAParams, [PhylogeneticGraph])
                  -> ([PhylogeneticGraph], Int)
addDeleteNetEdges inGS inData rSeed maxNetEdges numToKeep maxRounds counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
   if isNothing inSimAnnealParams then addDeleteNetEdges' inGS inData rSeed maxNetEdges numToKeep maxRounds counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) Nothing inPhyloGraphList
   else
      let -- create list of params with unique list of random values for rounds of annealing
         annealingRounds = rounds $ fromJust inSimAnnealParams
         saPAramList = (U.generateUniqueRandList annealingRounds inSimAnnealParams) -- (replicate annealingRounds inPhyloGraphList)

         (annealRoundsList, counterList) = unzip (PU.seqParMap (parStrategy $ lazyParStrat inGS) (addDeleteNetEdges'' inGS inData rSeed maxNetEdges numToKeep maxRounds counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost)) (zip saPAramList (replicate annealingRounds inPhyloGraphList)))

      in
      (GO.selectGraphs Best numToKeep 0.0 (-1) (concat annealRoundsList) , sum counterList)

-- | addDeleteNetEdges'' is wrapper around addDeleteNetEdges' to use parmap
addDeleteNetEdges'' :: GlobalSettings
                    -> ProcessedData
                    -> Int
                    -> Int
                    -> Int
                    -> Int
                    -> Int
                    -> Bool
                    -> Bool
                    -> Bool
                    -> ([PhylogeneticGraph], VertexCost)
                    -> (Maybe SAParams, [PhylogeneticGraph])
                    -> ([PhylogeneticGraph], Int)
addDeleteNetEdges'' inGS inData rSeed maxNetEdges numToKeep maxRounds counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
   addDeleteNetEdges' inGS inData rSeed maxNetEdges numToKeep maxRounds counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) inSimAnnealParams inPhyloGraphList

-- | addDeleteNetEdges' removes each edge and adds an edge to all possible places (or steepest) each round
-- until no better or additional graphs are found (or max rounds met)
-- call with ([], infinity) [single input graph]
-- doesn't have to be random, but likely to converge quickly if not
addDeleteNetEdges' :: GlobalSettings
                   -> ProcessedData
                   -> Int
                   -> Int
                   -> Int
                   -> Int
                   -> Int
                   -> Bool
                   -> Bool
                   -> Bool
                   -> ([PhylogeneticGraph], VertexCost)
                   -> Maybe SAParams
                   -> [PhylogeneticGraph]
                   -> ([PhylogeneticGraph], Int)
addDeleteNetEdges' inGS inData rSeed maxNetEdges numToKeep maxRounds counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) inSimAnnealParams inPhyloGraphList =
   if null inPhyloGraphList then (take numToKeep curBestGraphList, counter)

   -- if hit maxmimum rounds then return
   else if counter == maxRounds then (take numToKeep curBestGraphList, counter)

   -- other wise add/delete
   else
      -- trace ("\tRound " ++ (show counter)) (
      let -- insert edges first
          randIntList = randomIntList rSeed
          (insertGraphList, _) = insertAllNetEdges' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) randIntList inSimAnnealParams inPhyloGraphList

          -- this to update randlists in SAPArams for subsequent calls
          updatedSAParamList = if isJust inSimAnnealParams then U.generateUniqueRandList 2 inSimAnnealParams
                               else [Nothing, Nothing]

         -- if no better--take input for delte phase
          randIntList2 = randomIntList (randIntList !! 1)
          (insertGraphList', insertGraphCost, toDeleteList) = if null insertGraphList then (curBestGraphList, curBestGraphCost, inPhyloGraphList)
                                                              else
                                                                  let newList = GO.selectGraphs Best (maxBound::Int) 0.0 (-1) insertGraphList
                                                                  in
                                                                  (newList, snd6 $ head newList, newList)

         -- delete edges
          (deleteGraphList, _) =  deleteAllNetEdges' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (insertGraphList', insertGraphCost) randIntList2 (head updatedSAParamList) toDeleteList

         -- gather beter if any
          (newBestGraphList, newBestGraphCost, graphsToDoNext) = if null deleteGraphList then (curBestGraphList, curBestGraphCost, inPhyloGraphList)
                                                                 else
                                                                      let newDeleteGraphs = GO.selectGraphs Best (maxBound::Int) 0.0 (-1) deleteGraphList
                                                                      in
                                                                      (newDeleteGraphs, snd6 $ head newDeleteGraphs, newDeleteGraphs)
      in
      -- check is same then return
      if newBestGraphCost == curBestGraphCost then
         (take numToKeep curBestGraphList, counter)

      -- if better (or nothing) keep going
      else
         addDeleteNetEdges' inGS inData (randIntList !! 2) maxNetEdges numToKeep maxRounds (counter + 1) returnMutated doSteepest doRandomOrder (newBestGraphList, newBestGraphCost) (last updatedSAParamList) graphsToDoNext
      -- )

-- | moveAllNetEdges is a wrapper for moveAllNetEdges' allowing for multiple simulated annealing rounds
moveAllNetEdges :: GlobalSettings
                -> ProcessedData
                -> Int
                -> Int
                -> Int
                -> Int
                -> Bool
                -> Bool
                -> Bool
                -> ([PhylogeneticGraph], VertexCost)
                -> (Maybe SAParams, [PhylogeneticGraph])
                -> ([PhylogeneticGraph], Int)
moveAllNetEdges inGS inData rSeed maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
   if isNothing inSimAnnealParams then moveAllNetEdges' inGS inData rSeed maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) Nothing inPhyloGraphList
   else
      let -- create list of params with unique list of random values for rounds of annealing
         annealingRounds = rounds $ fromJust inSimAnnealParams
         saPAramList = (U.generateUniqueRandList annealingRounds inSimAnnealParams) -- (replicate annealingRounds inPhyloGraphList)

         (annealRoundsList, counterList) = unzip (PU.seqParMap (parStrategy $ lazyParStrat inGS) (moveAllNetEdges'' inGS inData rSeed maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost)) (zip saPAramList (replicate annealingRounds inPhyloGraphList)))

      in
      (GO.selectGraphs Best numToKeep 0.0 (-1) (concat annealRoundsList) , sum counterList)

-- | moveAllNetEdges'' is wrapper around moveAllNetEdges' to use parmap
moveAllNetEdges'' :: GlobalSettings
                  -> ProcessedData
                  -> Int
                  -> Int
                  -> Int
                  -> Int
                  -> Bool
                  -> Bool
                  -> Bool
                  -> ([PhylogeneticGraph], VertexCost)
                  -> (Maybe SAParams, [PhylogeneticGraph])
                  -> ([PhylogeneticGraph], Int)
moveAllNetEdges'' inGS inData rSeed maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
   moveAllNetEdges' inGS inData rSeed maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) inSimAnnealParams inPhyloGraphList

-- | moveAllNetEdges' removes each edge and adds an edge to all possible places (or steepest) each round
-- until no better or additional graphs are found
-- call with ([], infinity) [single input graph]
moveAllNetEdges' :: GlobalSettings
                 -> ProcessedData
                 -> Int
                 -> Int
                 -> Int
                 -> Int
                 -> Bool
                 -> Bool
                 -> Bool
                 -> ([PhylogeneticGraph], VertexCost)
                 -> Maybe SAParams
                 -> [PhylogeneticGraph]
                 -> ([PhylogeneticGraph], Int)
moveAllNetEdges' inGS inData rSeed maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) inSimAnnealParams inPhyloGraphList =
   if null inPhyloGraphList then (take numToKeep curBestGraphList, counter)
   else
      let firstPhyloGraph = head inPhyloGraphList
          currentCost = min curBestGraphCost (snd6 firstPhyloGraph)

          -- randomize order of edges to try moving
          netEdgeList = if not doRandomOrder then
                           LG.labNetEdges (thd6 firstPhyloGraph)
                        else
                           permuteList rSeed $ LG.labNetEdges (thd6 firstPhyloGraph)


          newGraphList' =  deleteOneNetAddAll inGS inData maxNetEdges numToKeep doSteepest doRandomOrder firstPhyloGraph (fmap LG.toEdge netEdgeList) rSeed inSimAnnealParams
          newGraphList = GO.selectGraphs Best numToKeep 0.0 (-1) newGraphList'
          newGraphCost = if (not . null) newGraphList' then snd6 $ head newGraphList
                         else infinity
      in

      -- if graph is a tree no edges to delete
      {-
      if LG.isTree (fst6 firstPhyloGraph) then
         trace ("\tGraph in move network edges is tree--skipping")
         moveAllNetEdges' inGS inData rSeed maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (firstPhyloGraph : curBestGraphList, currentCost) inSimAnnealParams (tail inPhyloGraphList)
      -}

      -- if graph is a tree no edges to delete
      if null netEdgeList then trace ("\t\tGraph in move has no network edges to move--skipping") (inPhyloGraphList, counter)

      -- regular move keeping best
      else if isNothing inSimAnnealParams then
         if newGraphCost > currentCost then
            -- trace ("\t MANE : Worse")
            moveAllNetEdges' inGS inData maxNetEdges rSeed numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (firstPhyloGraph : curBestGraphList, currentCost) inSimAnnealParams (tail inPhyloGraphList)
         else if newGraphCost < currentCost then
            -- trace ("\tMANE-> " ++ (show newGraphCost)) (
            if doSteepest then
               moveAllNetEdges' inGS inData rSeed maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (newGraphList, newGraphCost) inSimAnnealParams newGraphList
            else
               moveAllNetEdges' inGS inData rSeed maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (newGraphList, newGraphCost) inSimAnnealParams (newGraphList ++ (tail inPhyloGraphList))
            -- )

         else
            -- new graph list contains the input graph if equal and filterd unique already in moveAllNetEdges
            let newCurSameBestList = GO.selectGraphs Unique numToKeep 0.0 (-1) (curBestGraphList ++ newGraphList)
            in
            -- trace ("\t MANE : Equal")
            moveAllNetEdges' inGS inData rSeed maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (newCurSameBestList, currentCost) inSimAnnealParams (tail inPhyloGraphList)

      -- sim anneal choice
      else if True then errorWithoutStackTrace "Simulated Annealing/Drift not implemented for Network Move"
      else
            let -- abstract stopping criterion to continue
               numDone = if (method $ fromJust inSimAnnealParams) == SimAnneal then currentStep $ fromJust inSimAnnealParams
                         else driftChanges $ fromJust inSimAnnealParams
               numMax  = if (method $ fromJust inSimAnnealParams) == SimAnneal then numberSteps $ fromJust inSimAnnealParams
                         else driftMaxChanges $ fromJust inSimAnnealParams

               -- get acceptance based on heuristic costs
               uniqueGraphList = GO.selectGraphs Unique numToKeep 0.0 (-1) newGraphList'
               annealBestCost = if (not . null) uniqueGraphList then min curBestGraphCost (snd6 $ head uniqueGraphList)
                                 else curBestGraphCost
               (acceptFirstGraph, newSAParams) = if (not . null) uniqueGraphList then U.simAnnealAccept inSimAnnealParams annealBestCost (snd6 $ head uniqueGraphList)
                                                  else (False, U.incrementSimAnnealParams inSimAnnealParams)
            in
            -- trace ("ACG" ++ (show acceptFirstGraph) ++ " " ++ (show $ snd6 $ head uniqueGraphList)) (
            if (numDone < numMax) then
               -- this fixes tail fail
               let nextUniqueList = if (not . null) uniqueGraphList then tail uniqueGraphList
                                    else []
               in

               if acceptFirstGraph then
                  moveAllNetEdges' inGS inData rSeed maxNetEdges  numToKeep (counter + 1) returnMutated doSteepest doRandomOrder ((head uniqueGraphList) :  curBestGraphList, annealBestCost) newSAParams (nextUniqueList ++ (tail inPhyloGraphList))
               else
                  moveAllNetEdges' inGS inData rSeed maxNetEdges  numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (curBestGraphList, annealBestCost) newSAParams (nextUniqueList ++ (tail inPhyloGraphList))

            -- if want non-optimized list for GA or whatever
            else if returnMutated then (take numToKeep curBestGraphList, counter)

            -- optimize list and return
            else
               let (bestMoveList', counter') =  moveAllNetEdges' inGS inData rSeed maxNetEdges numToKeep (counter + 1) False doSteepest doRandomOrder ([], annealBestCost) Nothing (take numToKeep curBestGraphList)
                   bestMoveList = GO.selectGraphs Best numToKeep 0.0 (-1) bestMoveList'
               in
               --trace ("BM: " ++ (show $ snd6 $ head  bestMoveList))
               (take numToKeep bestMoveList, counter')
      -- )

-- | (curBestGraphList, annealBestCost) is a wrapper for moveAllNetEdges' allowing for multiple simulated annealing rounds
insertAllNetEdges :: GlobalSettings
                  -> ProcessedData
                  -> Int
                  -> Int
                  -> Int
                  -> Int
                  -> Int
                  -> Bool
                  -> Bool
                  -> Bool
                  -> ([PhylogeneticGraph], VertexCost)
                  -> (Maybe SAParams, [PhylogeneticGraph])
                  -> ([PhylogeneticGraph], Int)
insertAllNetEdges inGS inData rSeed maxNetEdges numToKeep maxRounds counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
   if isNothing inSimAnnealParams then

      -- check for multiple rounds of addition--if > 1 then need to randomize order
      if maxRounds == 1 then
            insertAllNetEdges' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (randomIntList rSeed) Nothing inPhyloGraphList
      else
         -- need to concat and send different randomization lists for each "round"
         let randSeedList = take maxRounds (randomIntList rSeed)
             randIntListList = fmap randomIntList randSeedList
             (insertGraphList, counterList) = unzip $ PU.seqParMap (parStrategy $ lazyParStrat inGS) (insertAllNetEdgesRand inGS inData maxNetEdges numToKeep counter returnMutated doSteepest (curBestGraphList, curBestGraphCost) Nothing inPhyloGraphList) randIntListList
         in
         -- insert functions take care of returning "better" or empty
         -- should be empty if nothing better
         (GO.selectGraphs Best numToKeep 0.0 (-1) (concat insertGraphList), sum counterList)
   else
      let -- create list of params with unique list of random values for rounds of annealing
          annealingRounds = rounds $ fromJust inSimAnnealParams
          annealParamGraphList = U.generateUniqueRandList annealingRounds inSimAnnealParams
          replicateRandIntList = fmap randomIntList (take annealingRounds (randomIntList rSeed))

          (annealRoundsList, counterList) = unzip (PU.seqParMap (parStrategy $ lazyParStrat inGS) (insertAllNetEdges'' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost)) (zip3 replicateRandIntList annealParamGraphList (replicate annealingRounds inPhyloGraphList)))
      in
      if (not returnMutated) || isNothing inSimAnnealParams then
         (GO.selectGraphs Best numToKeep 0.0 (-1) (concat annealRoundsList) , sum counterList)
      else
         (GO.selectGraphs Unique numToKeep 0.0 (-1) (concat annealRoundsList) , sum counterList)

-- | insertAllNetEdgesRand is a wrapper around insertAllNetEdges'' to pass unique randomLists to insertAllNetEdges'
insertAllNetEdgesRand :: GlobalSettings
                      -> ProcessedData
                      -> Int
                      -> Int
                      -> Int
                      -> Bool
                      -> Bool
                      -> ([PhylogeneticGraph], VertexCost)
                      -> Maybe SAParams
                      -> [PhylogeneticGraph]
                      -> [Int]
                      -> ([PhylogeneticGraph], Int)
insertAllNetEdgesRand inGS inData maxNetEdges numToKeep counter returnMutated doSteepest (curBestGraphList, curBestGraphCost) inSimAnnealParams inPhyloGraphList randIntList =
   insertAllNetEdges' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest True (curBestGraphList, curBestGraphCost) randIntList inSimAnnealParams inPhyloGraphList

-- | insertAllNetEdges'' is a wrapper around insertAllNetEdges' to allow for seqParMap
insertAllNetEdges'' :: GlobalSettings
                    -> ProcessedData
                    -> Int
                    -> Int
                    -> Int
                    -> Bool
                    -> Bool
                    -> Bool
                    -> ([PhylogeneticGraph], VertexCost)
                    -> ([Int], Maybe SAParams, [PhylogeneticGraph])
                    -> ([PhylogeneticGraph], Int)
insertAllNetEdges'' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (randIntList, inSimAnnealParams, inPhyloGraphList) =
      insertAllNetEdges' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) randIntList inSimAnnealParams inPhyloGraphList


-- | insertAllNetEdges' adds network edges one each each round until no better or additional
-- graphs are found
-- call with ([], infinity) [single input graph]
insertAllNetEdges' :: GlobalSettings
                   -> ProcessedData
                   -> Int
                   -> Int
                   -> Int
                   -> Bool
                   -> Bool
                   -> Bool
                   -> ([PhylogeneticGraph], VertexCost)
                   -> [Int]
                   -> Maybe SAParams
                   -> [PhylogeneticGraph]
                   -> ([PhylogeneticGraph], Int)
insertAllNetEdges' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) randIntList inSimAnnealParams inPhyloGraphList =
   if null inPhyloGraphList then
      -- this logic so don't return mutated if finish insertion before hitting other stopping points
      -- and don't want mutated--straight deletion on current best graphs
      if isNothing inSimAnnealParams || returnMutated then (take numToKeep curBestGraphList, counter)
      else
         deleteAllNetEdges' inGS inData maxNetEdges numToKeep counter False doSteepest doRandomOrder ([], curBestGraphCost) (randomIntList $ head randIntList) Nothing (take numToKeep curBestGraphList)
   else
      let currentCost = min curBestGraphCost (snd6 $ head inPhyloGraphList)

           --check for max net edges
          (_, _, _, netNodes) = LG.splitVertexList (thd6 $ head inPhyloGraphList)

          (newGraphList, _, newSAParams) = insertEachNetEdge inGS inData (head randIntList) maxNetEdges numToKeep doSteepest doRandomOrder Nothing inSimAnnealParams (head inPhyloGraphList)


          bestNewGraphList = GO.selectGraphs Best numToKeep 0.0 (-1) newGraphList
          newGraphCost = if (not . null) bestNewGraphList then snd6 $ head bestNewGraphList
                         else infinity

      in
      trace ("\t\tNumber of network edges: " ++ (show $ length netNodes)) (
      -- trace ("IANE: " ++ (show $ length netNodes)) (
      if length netNodes >= maxNetEdges then
         trace ("Maximum number of network edges reached: " ++ (show $ length netNodes))
         (take numToKeep curBestGraphList, counter)


      else if null newGraphList then
         trace ("\t\tNumber of network edges: " ++ (show $ length netNodes))
         (take numToKeep curBestGraphList, counter)

      -- regular insert keeping best
      else if isNothing inSimAnnealParams then
         postProcessNetworkAdd inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (newGraphList, newGraphCost) (tail randIntList) inSimAnnealParams netNodes currentCost (tail inPhyloGraphList)

      -- simulated annealing--needs new SAParams
      else
         -- trace ("IANE: " ++ (show $ method $ fromJust inSimAnnealParams)) (
         let (saGraphs, saCounter) = postProcessNetworkAddSA inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (newGraphList, newGraphCost) (tail randIntList) newSAParams netNodes currentCost (tail inPhyloGraphList)
         in

         --if want mutated then return
         if returnMutated then (saGraphs, saCounter)

         -- delete non-minimal edges if any
         -- not sim anneal/drift regular optimal searching
         else
            let annealBestCost = minimum $ fmap snd6 saGraphs
                (bestList', counter') =  deleteAllNetEdges' inGS inData maxNetEdges numToKeep saCounter False doSteepest doRandomOrder (saGraphs, annealBestCost) (randomIntList $ head randIntList) Nothing saGraphs
                bestList = GO.selectGraphs Best numToKeep 0.0 (-1) bestList'
            in
            (bestList, counter')
         -- )
      )


-- | postProcessNetworkAddSA processes simaneal/drift
-- assumes SAParams are updated during return of graph list above
postProcessNetworkAddSA :: GlobalSettings
                      -> ProcessedData
                      -> Int
                      -> Int
                      -> Int
                      -> Bool
                      -> Bool
                      -> Bool
                      -> ([PhylogeneticGraph], VertexCost)
                      -> ([PhylogeneticGraph], VertexCost)
                      -> [Int]
                      -> Maybe SAParams
                      -> [LG.LNode VertexInfo]
                      -> VertexCost
                      -> [PhylogeneticGraph]
                      -> ([PhylogeneticGraph], Int)
postProcessNetworkAddSA inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (newGraphList, newGraphCost) randIntList inSimAnnealParams _ _ inPhyloGraphList =

   -- trace ("\t\tNumber of network edges: " ++ (show $ length netNodes)) (
   -- this to deal with empty list issues if nothing found
   let (nextNewGraphList, firstNewGraphList) = if (not . null) newGraphList then (tail newGraphList, [head newGraphList])
                                               else ([], [])
       graphsToInsert = if doSteepest then newGraphList
                        else take numToKeep $ newGraphList ++ inPhyloGraphList
   in

   -- always accept if found better
   if newGraphCost < curBestGraphCost then
      trace ("\t-> " ++ (show newGraphCost))
      insertAllNetEdges' inGS inData maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (newGraphList, newGraphCost) randIntList inSimAnnealParams graphsToInsert

   -- check if hit max change/ cooling steps
   else if ((currentStep $ fromJust inSimAnnealParams) >= (numberSteps $ fromJust inSimAnnealParams)) || ((driftChanges $ fromJust inSimAnnealParams) >= (driftMaxChanges $ fromJust inSimAnnealParams)) then
         --trace ("PPA return: " ++ (show (newMinCost, curBestCost)))
         (GO.selectGraphs Unique numToKeep 0.0 (-1) (newGraphList ++ curBestGraphList), counter)


    -- more to do
    else
       let annealBestCost = min curBestGraphCost newGraphCost
       in
       insertAllNetEdges' inGS inData maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (firstNewGraphList ++  curBestGraphList, annealBestCost) randIntList inSimAnnealParams (nextNewGraphList ++ inPhyloGraphList)

       -- )

-- | postProcessNetworkAdd prcesses non-simaneal/drift--so no updating of SAParams
postProcessNetworkAdd :: GlobalSettings
                      -> ProcessedData
                      -> Int
                      -> Int
                      -> Int
                      -> Bool
                      -> Bool
                      -> Bool
                      -> ([PhylogeneticGraph], VertexCost)
                      -> ([PhylogeneticGraph], VertexCost)
                      -> [Int]
                      -> Maybe SAParams
                      -> [LG.LNode VertexInfo]
                      -> VertexCost
                      -> [PhylogeneticGraph]
                      -> ([PhylogeneticGraph], Int)
postProcessNetworkAdd inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, _) (newGraphList, newGraphCost) randIntList inSimAnnealParams _ currentCost inPhyloGraphList =

   -- "steepest style descent" abandons existing list if better cost found
   -- trace ("\t\tNumber of network edges: " ++ (show $ length netNodes)) (
   if newGraphCost < currentCost then
      -- check if graph OK--done in insert function
      let --- isCyclicList = filter (== True) $ fmap LG.cyclic $ fmap thd6 newGraphList
          --- hasDupEdges = filter (== True) $ fmap LG.hasDuplicateEdge $ fmap thd6 newGraphList
          graphsToInsert = if doSteepest then newGraphList
                           else take numToKeep $ newGraphList ++ inPhyloGraphList
      in
      trace ("\t-> " ++ (show newGraphCost))
      insertAllNetEdges' inGS inData maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (newGraphList, newGraphCost) randIntList inSimAnnealParams graphsToInsert

   -- worse graphs found--go on
   else if  newGraphCost > currentCost then
      -- trace ("IANE: Worse")
      insertAllNetEdges' inGS inData maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (curBestGraphList, currentCost) randIntList inSimAnnealParams inPhyloGraphList

   -- equal cost
   -- not sure if should add new graphs to queue to do edge deletion again
   else
      -- new graph list contains the input graph if equal and filterd unique already in insertAllNetEdges
      let newCurSameBestList = GO.selectGraphs Unique numToKeep 0.0 (-1) (curBestGraphList ++ newGraphList)
      in
      -- trace ("IANE: same " ++ (show $ length (tail inPhyloGraphList)))
      insertAllNetEdges' inGS inData maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (newCurSameBestList, currentCost) randIntList inSimAnnealParams inPhyloGraphList
      -- )


-- | insertEachNetEdge takes a phylogenetic graph and inserts all permissible network edges one at time
-- and returns unique list of new Phylogenetic Graphs and cost
-- even if worse--could be used for simulated annealing later
-- if equal returns unique graph list
insertEachNetEdge :: GlobalSettings
                  -> ProcessedData
                  -> Int
                  -> Int
                  -> Int
                  -> Bool
                  -> Bool
                  -> Maybe VertexCost
                  -> Maybe SAParams
                  -> PhylogeneticGraph
                  -> ([PhylogeneticGraph], VertexCost, Maybe SAParams)
insertEachNetEdge inGS inData rSeed maxNetEdges numToKeep doSteepest doRandomOrder preDeleteCost inSimAnnealParams inPhyloGraph =
   if LG.isEmpty $ fst6 inPhyloGraph then error "Empty input insertEachNetEdge graph in deleteAllNetEdges"
   else
      let currentCost = if isNothing preDeleteCost then snd6 inPhyloGraph
                        else fromJust preDeleteCost

          (_, _, _, netNodes) = LG.splitVertexList (thd6 inPhyloGraph)

          candidateNetworkEdgeList' = getPermissibleEdgePairs inGS (thd6 inPhyloGraph)

          -- radomize pair list
          rSeedList = randomIntList rSeed
          candidateNetworkEdgeList = if doRandomOrder then permuteList (head rSeedList) candidateNetworkEdgeList'
                                     else candidateNetworkEdgeList'

          -- newGraphList = concat (fmap (insertNetEdgeBothDirections inGS inData inPhyloGraph) candidateNetworkEdgeList `using`  PU.myParListChunkRDS)
          (newGraphList, newSAParams) = if not doSteepest then
                                             let genNewSimAnnealParams = if isNothing inSimAnnealParams then Nothing
                                                                         else U.incrementSimAnnealParams inSimAnnealParams
                                             in
                                             (filter (/= emptyPhylogeneticGraph) (PU.seqParMap (parStrategy $ lazyParStrat inGS) (insertNetEdge inGS inData inPhyloGraph preDeleteCost) candidateNetworkEdgeList), genNewSimAnnealParams)
                                        else insertNetEdgeRecursive inGS inData (tail rSeedList) maxNetEdges doSteepest doRandomOrder inPhyloGraph preDeleteCost inSimAnnealParams candidateNetworkEdgeList

          minCost = if null candidateNetworkEdgeList || null newGraphList then infinity
                    else minimum $ fmap snd6 newGraphList
      in
      trace ("\tExamining at most " ++ (show $ length candidateNetworkEdgeList) ++ " candidate edge pairs") (

      -- no network edges to insert
      -- trace ("IENE: " ++ (show minCost)) (
      if length netNodes >= maxNetEdges then
            trace ("Maximum number of network edges reached: " ++ (show $ length netNodes))
            ([inPhyloGraph], snd6 inPhyloGraph, inSimAnnealParams)

      -- no edges to add
      else if null candidateNetworkEdgeList then
         -- trace ("IENE num cand edges:" ++ (show $ length candidateNetworkEdgeList))
         ([inPhyloGraph], currentCost, newSAParams)

      -- single if steepest so no need to unique
      else if doSteepest then
         --  trace ("IENE: All " ++ (show minCost))
         (GO.selectGraphs Best numToKeep 0.0 (-1) $ newGraphList, minCost, newSAParams)

      -- "all" option needs to recurse since does all available edges at each step
      -- logic is here since not in the deleteNetEdge function
      else if isNothing inSimAnnealParams then
         if minCost < currentCost then
            let  annealParamList = U.generateUniqueRandList (length newGraphList) newSAParams
                 allRandIntList = take (length newGraphList) (randomIntList (rSeedList !! 1))
                 (allGraphListList, costList, allSAParamList) = unzip3 $ PU.seqParMap (parStrategy $ lazyParStrat inGS)  (insertEachNetEdge' inGS inData maxNetEdges numToKeep doSteepest doRandomOrder preDeleteCost) (zip3 allRandIntList annealParamList newGraphList)
                 (allMinCost, allMinCostGraphs)  = if (not . null . concat) allGraphListList then
                                                  (minimum costList, GO.selectGraphs Unique numToKeep 0.0 (-1) $ concat allGraphListList)
                                                   else (infinity, [])
            in
            (allMinCostGraphs, allMinCost, U.incrementSimAnnealParams $ head allSAParamList)

         else
            (GO.selectGraphs Unique numToKeep 0.0 (-1) $ newGraphList, minCost, newSAParams)

      -- SA anneal/Drift
      else
         -- always take better
         if minCost < currentCost then
               (newGraphList, minCost, newSAParams)

         -- check if hit step limit--more for SA than drift
         else if ((currentStep $ fromJust inSimAnnealParams) >= (numberSteps $ fromJust inSimAnnealParams)) || ((driftChanges $ fromJust inSimAnnealParams) >= (driftMaxChanges $ fromJust inSimAnnealParams)) then
               ([inPhyloGraph], snd6 inPhyloGraph, inSimAnnealParams)

         -- otherwise do the anneal/Drift accept, or keep going on input graph
         else
            let (acceptGraph, nextSAParams) = U.simAnnealAccept inSimAnnealParams currentCost minCost
            in
            if acceptGraph then
               (newGraphList, minCost, newSAParams)
            else
               insertEachNetEdge inGS inData (head $ randomIntList rSeed) maxNetEdges numToKeep doSteepest doRandomOrder preDeleteCost nextSAParams inPhyloGraph
      ) -- )

-- | insertEachNetEdge' is a wrapper around insertEachNetEdge to allow for parmapping with multiple parameters
insertEachNetEdge'   :: GlobalSettings
                     -> ProcessedData
                     -> Int
                     -> Int
                     -> Bool
                     -> Bool
                     -> Maybe VertexCost
                     -> (Int, Maybe SAParams, PhylogeneticGraph)
                     -> ([PhylogeneticGraph], VertexCost, Maybe SAParams)
insertEachNetEdge' inGS inData maxNetEdges numToKeep doSteepest doRandomOrder preDeleteCost (rSeed, inSimAnnealParams, inPhyloGraph) =
   insertEachNetEdge inGS inData rSeed maxNetEdges numToKeep doSteepest doRandomOrder preDeleteCost inSimAnnealParams inPhyloGraph

-- | insertNetEdgeRecursive recursively inserts edges and returns new graph only if better
-- if parallel evaluated numThreads each time (steepest scenario)
insertNetEdgeRecursive :: GlobalSettings
                       -> ProcessedData
                       -> [Int]
                       -> Int
                       -> Bool
                       -> Bool
                       ->  PhylogeneticGraph
                       -> Maybe VertexCost
                       -> Maybe SAParams
                       -> [(LG.LEdge EdgeInfo, LG.LEdge EdgeInfo)]
                       -> ([PhylogeneticGraph], Maybe SAParams)
insertNetEdgeRecursive inGS inData rSeedList maxNetEdges doSteepest doRandomOrder inPhyloGraph preDeleteCost inSimAnnealParams inEdgePairList =
   -- trace ("Edges pairs to go : " ++ (show $ length edgePairList)) (
   if null inEdgePairList then ([inPhyloGraph], inSimAnnealParams)
   else
      -- don't want to over saturate the parallel thread system
      let {-saRounds = if isNothing inSimAnnealParams then 1
                     else rounds $ fromJust inSimAnnealParams
          (numGraphsToExamine, _) = divMod PU.getNumThreads saRounds -- this may not "drift" if finds alot better, but that's how its supposed to work
          -}

          numGraphsToExamine = min (graphsSteepest inGS) PU.getNumThreads
          -- firstEdgePair = head edgePairList
          edgePairList = take numGraphsToExamine inEdgePairList

          --check for max net edges
          (_, _, _, netNodes) = LG.splitVertexList (thd6 inPhyloGraph)

          -- need to check display/character trees not conical graph
          -- newGraph = insertNetEdge inGS inData leafGraph inPhyloGraph preDeleteCost firstEdgePair
          -- these graph costs are "exact" or at least non-heuristic--needs to be updated when get a good heuristic
          newGraphList'' = PU.seqParMap (parStrategy $ lazyParStrat inGS) (insertNetEdge inGS inData inPhyloGraph preDeleteCost) edgePairList
          newGraphList' = filter (/= emptyPhylogeneticGraph) newGraphList''
          newGraphList = GO.selectGraphs Best (maxBound::Int) 0.0 (-1) newGraphList'
          newGraphCost = snd6 $ head newGraphList



      in
      -- traceNoLF ("*")  (      -- trace ("INER: " ++ (show $ snd6 newGraph) ++ " " ++ (show preDeleteCost)) (
      if length netNodes >= maxNetEdges then
         trace ("Maximum number of network edges reached: " ++ (show $ length netNodes))
         ([inPhyloGraph], inSimAnnealParams)

      -- malformed graph--returns nothing for either regular or simAnneal/drift
      else if null newGraphList' then
         -- trace ("INER: Empty more to go : " ++ (show $ length $ tail edgePairList))
         insertNetEdgeRecursive inGS inData rSeedList maxNetEdges doSteepest doRandomOrder inPhyloGraph preDeleteCost inSimAnnealParams (drop numGraphsToExamine inEdgePairList)

      -- "regular" insert, within steepest
      else if isNothing inSimAnnealParams then
         -- better cost
         if newGraphCost < snd6 inPhyloGraph then
            -- cyclic check in insert edge function
            -- trace ("INER: Better -> " ++ (show $ snd6 newGraph))
            (newGraphList, inSimAnnealParams)

         -- not better
         else -- trace ("INER: Really Not Better")
            insertNetEdgeRecursive inGS inData rSeedList maxNetEdges doSteepest doRandomOrder inPhyloGraph preDeleteCost inSimAnnealParams (drop numGraphsToExamine inEdgePairList)

      -- sim annealing/drift
      else
         -- trace ("IENR:" ++ (show (newGraphCost, snd6 inPhyloGraph)) ++ " params: " ++ (show (currentStep $ fromJust inSimAnnealParams, numberSteps $ fromJust inSimAnnealParams, driftChanges $ fromJust inSimAnnealParams, driftMaxChanges $ fromJust inSimAnnealParams))) (
         -- if better always accept
         if newGraphCost < snd6 inPhyloGraph then
            -- cyclic check in insert edge function
            -- trace ("INER: Better -> " ++ (show $ snd6 newGraph))
            -- these graph costs are "exact" or at least non-heuristic--needs to be updated when get a good heuristic
            let (_, nextSAParams) = U.simAnnealAccept inSimAnnealParams (snd6 inPhyloGraph) newGraphCost
            in
            (newGraphList, nextSAParams)

         -- check if hit step limit--more for SA than drift
         else if ((currentStep $ fromJust inSimAnnealParams) >= (numberSteps $ fromJust inSimAnnealParams)) || ((driftChanges $ fromJust inSimAnnealParams) >= (driftMaxChanges $ fromJust inSimAnnealParams)) then
                  ([inPhyloGraph], inSimAnnealParams)

         -- otherwise do the anneal/Drift accept
         else
            let (acceptGraph, nextSAParams) = U.simAnnealAccept inSimAnnealParams (snd6 inPhyloGraph) newGraphCost
            in
            if acceptGraph then (newGraphList, nextSAParams)
            else insertNetEdgeRecursive inGS inData rSeedList maxNetEdges doSteepest doRandomOrder inPhyloGraph preDeleteCost nextSAParams (drop numGraphsToExamine inEdgePairList)

      -- )
      -- )

-- | insertNetEdge inserts an edge between two other edges, creating 2 new nodes and rediagnoses graph
-- contacts deletes 2 orginal edges and adds 2 nodes and 5 new edges
-- does not check any edge reasonable-ness properties
-- new edge directed from first to second edge
-- naive for now
-- predeletecost of edge move
-- no choice of graph--just makes and returns
insertNetEdge :: GlobalSettings
              -> ProcessedData
              -> PhylogeneticGraph
              -> Maybe VertexCost
              -> (LG.LEdge b, LG.LEdge b)
              -> PhylogeneticGraph
insertNetEdge inGS inData inPhyloGraph _ edgePair@((u,v, _), (u',v', _)) =
   -- trace ("InsertEdge between: " ++ (show ((u,v), (u',v'))) )( -- ++ " into:\n " ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph $ thd6 inPhyloGraph)) (
   if LG.isEmpty $ thd6 inPhyloGraph then error "Empty input phylogenetic graph in insNetEdge"

   else
       let inSimple = fst6 inPhyloGraph

           -- get children of u' to make sure no net children--moved to permissiable edges
           --u'ChildrenNetNodes = filter (== True) $ fmap (LG.isNetworkNode inSimple) $ LG.descendants inSimple u'

           numNodes = length $ LG.nodes inSimple
           newNodeOne = (numNodes, TL.pack ("HTU" ++ (show numNodes)))
           newNodeTwo = (numNodes + 1, TL.pack ("HTU" ++ (show $ numNodes + 1)))
           newEdgeList = [(u, fst newNodeOne, 0.0),(fst newNodeOne, v, 0.0),(u', fst newNodeTwo, 0.0),(fst newNodeTwo, v', 0.0),(fst newNodeOne, fst newNodeTwo, 0.0)]
           edgesToDelete = [(u,v), (u',v')]
           newSimple = LG.insEdges newEdgeList $ LG.delEdges edgesToDelete $ LG.insNodes [newNodeOne, newNodeTwo] inSimple


           -- do not prune other edges if now unused
           pruneEdges = False

           -- don't warn that edges are being pruned
           warnPruneEdges = False

           -- graph optimization from root
           startVertex = Nothing

           -- full two-pass optimization
           leafGraph = LG.extractLeafGraph $ thd6 inPhyloGraph

           newPhyloGraph = T.multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex newSimple

           -- calculates heursitic graph delta
           -- (heuristicDelta, _, _, _, _)  = heuristicAddDelta inGS inPhyloGraph edgePair (fst newNodeOne) (fst newNodeTwo)
           heuristicDelta' = heuristicAddDelta' inGS inPhyloGraph edgePair


           edgeAddDelta = deltaPenaltyAdjustment inGS inPhyloGraph "add"

           heuristicFactor = (heuristicDelta' + edgeAddDelta) / edgeAddDelta

           -- use or not Net add heuristics
           metHeuristicThreshold = if useNetAddHeuristic inGS then heuristicFactor < (2/3)
                                   else True

       in

       -- remove these checks when working
       if not $ LG.isPhylogeneticGraph newSimple then
         emptyPhylogeneticGraph
       

       else
         -- need heuristics in here
         -- if (heuristicDelta + edgeAddDelta) < 0 then newPhyloGraph
         -- let oldPhyloGraph = T.multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex inSimple
         -- in
         -- trace ("INE: OK " ++ (show (numNodes, newNodeOne, newNodeTwo, newEdgeList, edgesToDelete, snd6 oldPhyloGraph)) ++ "\nOrig\n" ++ (-- if (heuristicDelta + edgeAddDelta) < 0 then newPhyloGraph
         -- trace ("INE: " ++ (show (heuristicDelta, heuristicDelta', edgeAddDelta, snd6 inPhyloGraph)) ++ " -> " ++ (show (heuristicDelta + edgeAddDelta + (snd6 inPhyloGraph), heuristicDelta' + edgeAddDelta + (snd6 inPhyloGraph), snd6 newPhyloGraph))) $
         if metHeuristicThreshold then 
            if (snd6 newPhyloGraph <= snd6 inPhyloGraph) then newPhyloGraph
            else emptyPhylogeneticGraph
         else emptyPhylogeneticGraph
         -- )


-- | (curBestGraphList, annealBestCost) is a wrapper for moveAllNetEdges' allowing for multiple simulated annealing rounds
deleteAllNetEdges :: GlobalSettings
                  -> ProcessedData
                  -> Int
                  -> Int
                  -> Int
                  -> Int
                  -> Bool
                  -> Bool
                  -> Bool
                  -> ([PhylogeneticGraph], VertexCost)
                  -> (Maybe SAParams, [PhylogeneticGraph])
                  -> ([PhylogeneticGraph], Int)
deleteAllNetEdges inGS inData rSeed maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
   if isNothing inSimAnnealParams then
      deleteAllNetEdges' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (randomIntList rSeed) inSimAnnealParams inPhyloGraphList
   else
      let -- create list of params with unique list of random values for rounds of annealing
          annealingRounds = rounds $ fromJust inSimAnnealParams
          annealParamGraphList = U.generateUniqueRandList annealingRounds inSimAnnealParams
          replicateRandIntList = fmap randomIntList (take annealingRounds (randomIntList rSeed))

          -- (annealRoundsList, counterList) = unzip (zipWith3 (deleteAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost)) replicateRandIntList annealParamGraphList (replicate annealingRounds inPhyloGraphList) `using` PU.myParListChunkRDS)
          (annealRoundsList, counterList) = unzip (PU.seqParMap (parStrategy $ lazyParStrat inGS) (deleteAllNetEdges'' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost)) (zip3 replicateRandIntList annealParamGraphList (replicate annealingRounds inPhyloGraphList)))
      in
      (GO.selectGraphs Best numToKeep 0.0 (-1) (concat annealRoundsList) , sum counterList)


-- | deleteAllNetEdges'' is a wrapper around deleteAllNetEdges' to allow use of seqParMap
deleteAllNetEdges'' :: GlobalSettings
                    -> ProcessedData
                    -> Int
                    -> Int
                    -> Int
                    -> Bool
                    -> Bool
                    -> Bool
                    -> ([PhylogeneticGraph], VertexCost)
                    -> ([Int], Maybe SAParams, [PhylogeneticGraph])
                    -> ([PhylogeneticGraph], Int)
deleteAllNetEdges'' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (randIntList, inSimAnnealParams, inPhyloGraphList) =
      deleteAllNetEdges' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) randIntList inSimAnnealParams inPhyloGraphList


-- | deleteAllNetEdges deletes network edges one each each round until no better or additional
-- graphs are found
-- call with ([], infinity) [single input graph]
deleteAllNetEdges' :: GlobalSettings
                   -> ProcessedData
                   -> Int
                   -> Int
                   -> Int
                   -> Bool
                   -> Bool
                   -> Bool
                   -> ([PhylogeneticGraph], VertexCost)
                   -> [Int]
                   -> Maybe SAParams
                   -> [PhylogeneticGraph]
                   -> ([PhylogeneticGraph], Int)
deleteAllNetEdges' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) randIntList inSimAnnealParams inPhyloGraphList =
   -- trace ("In deleteAllNetEdges " ++ (show $ length inPhyloGraphList)) (
   if null inPhyloGraphList then (take numToKeep curBestGraphList, counter)
   else
      let currentCost = min curBestGraphCost (snd6 $ head inPhyloGraphList)

          (newGraphList', _, newSAParams) = deleteEachNetEdge inGS inData (head randIntList) numToKeep doSteepest doRandomOrder False inSimAnnealParams (head inPhyloGraphList)

          newGraphList = GO.selectGraphs Best numToKeep 0.0 (-1) newGraphList'
          newGraphCost = if (not . null) newGraphList then snd6 $ head newGraphList
                         else infinity

      in
      -- trace ("DANE: " ++ (show (newGraphCost, length newGraphList))) (
      -- if graph is a tree no edges to delete
      if LG.isTree (fst6 $ head inPhyloGraphList) then
         -- let (a,b,c,d) = LG.splitVertexList (fst6 $ head inPhyloGraphList)
         -- in
         trace ("\tGraph in delete network edges is tree--skipping") --  :" ++ (show $ (snd6 $ head inPhyloGraphList, length a, length b, length c, length d)))
         deleteAllNetEdges' inGS inData maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder ((head inPhyloGraphList) : curBestGraphList, currentCost) (tail randIntList) inSimAnnealParams (tail inPhyloGraphList)

      -- is this an issue for SA?
      else if null newGraphList then (take numToKeep curBestGraphList, counter + 1)

      -- regular delete wihtout simulated annealing
      else if isNothing inSimAnnealParams then
            postProcessNetworkDelete inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (tail randIntList) inSimAnnealParams inPhyloGraphList newGraphList newGraphCost currentCost

      -- simulated annealing
      else
         let (saGraphs, saCounter) = postProcessNetworkDeleteSA inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (tail randIntList) newSAParams inPhyloGraphList newGraphList newGraphCost currentCost
         in

         --if want mutated then return
         if returnMutated then (saGraphs, saCounter)

         -- insert non-minimal edges if any
         -- not sim anneal/drift regular optimal searching
         else
            let annealBestCost = minimum $ fmap snd6 saGraphs
                (bestList', counter') =  insertAllNetEdges' inGS inData maxNetEdges numToKeep saCounter False doSteepest doRandomOrder (saGraphs, annealBestCost) (randomIntList $ head randIntList) Nothing saGraphs
                bestList = GO.selectGraphs Best numToKeep 0.0 (-1) bestList'
            in
            (bestList, counter')


-- | postProcessNetworkDeleteSA postprocesses results from delete actions for non-annealing/Drift network delete operations
-- assumes SAParams are updated during return of graph list above
postProcessNetworkDeleteSA :: GlobalSettings
                         -> ProcessedData
                         -> Int
                         -> Int
                         -> Int
                         -> Bool
                         -> Bool
                         -> Bool
                         -> ([PhylogeneticGraph], VertexCost)
                         -> [Int]
                         -> Maybe SAParams
                         -> [PhylogeneticGraph]
                         -> [PhylogeneticGraph]
                         -> VertexCost
                         -> VertexCost
                         -> ([PhylogeneticGraph], Int)
postProcessNetworkDeleteSA inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) randIntList inSimAnnealParams inPhyloGraphList newGraphList newGraphCost currentCost =

   -- this to deal with empty list issues if nothing found
   let (nextNewGraphList, firstNewGraphList) = if (not . null) newGraphList then (tail newGraphList, [head newGraphList])
                                               else ([], [])
       graphsToDelete = if doSteepest then newGraphList
                        else take numToKeep $ newGraphList ++ inPhyloGraphList
   in

   -- always accept if found better
   if newGraphCost < currentCost then
      trace ("\t-> " ++ (show newGraphCost))
      deleteAllNetEdges' inGS inData maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (newGraphList, newGraphCost) randIntList inSimAnnealParams graphsToDelete

   -- check if hit max change/ cooling steps
   else if ((currentStep $ fromJust inSimAnnealParams) >= (numberSteps $ fromJust inSimAnnealParams)) || ((driftChanges $ fromJust inSimAnnealParams) >= (driftMaxChanges $ fromJust inSimAnnealParams)) then
         --trace ("PPA return: " ++ (show (newMinCost, curBestCost)))
         (GO.selectGraphs Unique numToKeep 0.0 (-1) (newGraphList ++ curBestGraphList), counter)

   -- more to do
    else
       let annealBestCost = min curBestGraphCost newGraphCost
       in
       deleteAllNetEdges' inGS inData maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (firstNewGraphList ++  curBestGraphList, annealBestCost) randIntList inSimAnnealParams (nextNewGraphList ++ inPhyloGraphList)


-- | postProcessNetworkDelete postprocesses results from delete actions for "regular" ie non-annealing/Drift network delete operations
postProcessNetworkDelete :: GlobalSettings
                         -> ProcessedData
                         -> Int
                         -> Int
                         -> Int
                         -> Bool
                         -> Bool
                         -> Bool
                         -> ([PhylogeneticGraph], VertexCost)
                         -> [Int]
                         -> Maybe SAParams
                         -> [PhylogeneticGraph]
                         -> [PhylogeneticGraph]
                         -> VertexCost
                         -> VertexCost
                         -> ([PhylogeneticGraph], Int)
postProcessNetworkDelete inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, _) randIntList inSimAnnealParams inPhyloGraphList newGraphList newGraphCost currentCost =

   -- worse graphs found--go on
   if  newGraphCost > currentCost then deleteAllNetEdges' inGS inData maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder ((head inPhyloGraphList) : curBestGraphList, currentCost) randIntList inSimAnnealParams (tail inPhyloGraphList)

   -- "steepest style descent" abandons existing list if better cost found
   else if newGraphCost < currentCost then
      trace ("\t-> " ++ (show newGraphCost)) (
      if doSteepest then
         deleteAllNetEdges' inGS inData maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (newGraphList, newGraphCost) randIntList inSimAnnealParams newGraphList
      else
         deleteAllNetEdges' inGS inData maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (newGraphList, newGraphCost) randIntList inSimAnnealParams (newGraphList ++ (tail inPhyloGraphList))
      )

   -- equal cost
   -- not sure if should add new graphs to queue to do edge deletion again
   else
      -- new graph list contains the input graph if equal and filterd unique already in deleteEachNetEdge
      let newCurSameBestList =  GO.selectGraphs Unique numToKeep 0.0 (-1) (curBestGraphList ++ newGraphList)
      in
      deleteAllNetEdges' inGS inData maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (newCurSameBestList, currentCost) randIntList inSimAnnealParams (tail inPhyloGraphList)

-- | deleteOneNetAddAll version deletes net edges in turn and readds-based on original cost
-- but this cost in graph (really not correct) but allows logic of insert edge to function better
-- unlike deleteOneNetAddAll' only deals with single edge deletion at a time
deleteOneNetAddAll :: GlobalSettings
                   -> ProcessedData
                   -> Int
                   -> Int
                   -> Bool
                   -> Bool
                   -> PhylogeneticGraph
                   -> [LG.Edge]
                   -> Int
                   -> Maybe SAParams
                   -> [PhylogeneticGraph]
deleteOneNetAddAll inGS inData maxNetEdges numToKeep doSteepest doRandomOrder inPhyloGraph edgeToDeleteList rSeed inSimAnnealParams =
   if null edgeToDeleteList then
      -- trace ("\tGraph has no edges to move---skipping")
      [inPhyloGraph]
   else if LG.isEmpty $ thd6 inPhyloGraph then error "Empty graph in deleteOneNetAddAll"
   else
      -- trace ("DONAA-New: " ++ (show $ snd6 inPhyloGraph) ++ " Steepest:" ++ (show doSteepest)) (
      trace ("Moving " ++ (show $ length edgeToDeleteList) ++ " network edges, current best cost: " ++ (show $ snd6 inPhyloGraph)) (
      -- start with initial graph cost
      let inGraphCost = snd6 inPhyloGraph

          -- get deleted simple graphs and bool for changed
          delGraphBoolPair = deleteNetworkEdge (fst6 inPhyloGraph) (head edgeToDeleteList)

      in

      -- no change in network structure
      if snd delGraphBoolPair == False then
         deleteOneNetAddAll inGS inData maxNetEdges numToKeep doSteepest doRandomOrder inPhyloGraph (tail edgeToDeleteList) rSeed inSimAnnealParams

      else
          let simpleGraphToInsert = fst delGraphBoolPair

              (_, _, _, curNetNodes) = LG.splitVertexList simpleGraphToInsert
              curNumNetNodes = length curNetNodes

              -- optimize deleted graph and update cost with input cost
              leafGraph = LG.extractLeafGraph $ thd6 inPhyloGraph
   
              graphToInsert = T.multiTraverseFullyLabelSoftWired inGS inData False False leafGraph Nothing simpleGraphToInsert -- `using` PU.myParListChunkRDS

              -- keep same cost and just keep better--check if better than original later
              graphToInsert' = T.updatePhylogeneticGraphCost graphToInsert inGraphCost


              insertedGraphTripleList = insertEachNetEdge inGS inData rSeed (curNumNetNodes + 1) numToKeep doSteepest doRandomOrder Nothing inSimAnnealParams graphToInsert'

              newMinimumCost = snd3 insertedGraphTripleList

              newBestGraphs = filter ((== newMinimumCost) . snd6) $ fst3 insertedGraphTripleList

          in
          -- trace ("DONAA-New: " ++ (show (inGraphCost, fmap snd6 graphsToInsert, fmap snd6 graphsToInsert', newMinimumCost))) (
          if newMinimumCost < inGraphCost then
             -- trace ("DONA-> ")
             newBestGraphs

          else deleteOneNetAddAll inGS inData maxNetEdges numToKeep doSteepest doRandomOrder inPhyloGraph (tail edgeToDeleteList) rSeed inSimAnnealParams

          ) -- )


-- | getPermissibleEdgePairs takes a DecoratedGraph and returns the list of all pairs
-- of edges that can be joined by a network edge and meet all necessary conditions

-- add in other conditions
--   reproducable--ie not tree noide with two net node children--other stuff
getPermissibleEdgePairs ::  GlobalSettings -> DecoratedGraph -> [(LG.LEdge EdgeInfo, LG.LEdge EdgeInfo)]
getPermissibleEdgePairs inGS inGraph =
   if LG.isEmpty inGraph then error "Empty input graph in isEdgePairPermissible"
   else
       let edgeList = LG.labEdges inGraph

           -- edges to potentially conenct
           edgePairs = cartProd edgeList edgeList

           -- get coeval node pairs in existing grap
           coevalNodeConstraintList = LG.coevalNodePairs inGraph
           coevalNodeConstraintList' = PU.seqParMap (parStrategy $ lazyParStrat inGS) (LG.addBeforeAfterToPair inGraph) coevalNodeConstraintList -- `using`  PU.myParListChunkRDS

           edgeTestList = PU.seqParMap (parStrategy $ lazyParStrat inGS) (isEdgePairPermissible inGraph coevalNodeConstraintList') edgePairs -- `using`  PU.myParListChunkRDS
           pairList = fmap fst $ filter ((== True) . snd) $ zip edgePairs edgeTestList
       in
       -- trace ("Edge Pair list :" ++ (show $ fmap f pairList) ++ "\n"
       --  ++ "GPEP\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inGraph))
       pairList
       -- where f (a, b) = (LG.toEdge a, LG.toEdge b)

-- | isEdgePairPermissible takes a graph and two edges, coeval contraints, and tests whether a
-- pair of edges can be linked by a new edge and satify three consitions:
--    1) neither edge is a network edge
--    2) one edge cannot be "before" while the other is "after" in any of the constraint pairs
--    3) neither edge is an ancestor or descendent edge of the other (tested via bv of nodes)
-- the result should apply to a new edge in either direction
-- new edge to be creted is edge1 -> ege2
-- Could change to LG.isPhylogeneticGraph
isEdgePairPermissible :: DecoratedGraph -> [(LG.LNode a, LG.LNode a, [LG.LNode a], [LG.LNode a], [LG.LNode a], [LG.LNode a])] -> (LG.LEdge EdgeInfo, LG.LEdge EdgeInfo) -> Bool
isEdgePairPermissible inGraph constraintList (edge1@(u,v,_), edge2@(u',v',_)) =
   if LG.isEmpty inGraph then error "Empty input graph in isEdgePairPermissible"
   else
       if u == u' then False
       else if v == v' then False
       -- equality implied in above two
       -- else if LG.toEdge edge1 == LG.toEdge edge2 then False
       else if (LG.isNetworkNode inGraph u) || (LG.isNetworkNode inGraph u') then False
       else if (LG.isNetworkLabEdge inGraph edge1) || (LG.isNetworkLabEdge inGraph edge2) then False
       else if not (LG.meetsAllCoevalConstraintsNodes (fmap removeNodeLabels constraintList) edge1 edge2) then False
       else if (isAncDescEdge inGraph edge1 edge2) then False
       -- get children of u' to make sure no net children
       else if (not . null) $ filter (== True) $ fmap (LG.isNetworkNode inGraph) $ LG.descendants inGraph u' then False
       else True
   where removeNodeLabels (a,b,c,d,e,f) = (LG.toNode a, LG.toNode b, fmap LG.toNode c, fmap LG.toNode d, fmap LG.toNode e, fmap LG.toNode f)

-- | isAncDescEdge takes a graph and two edges and examines whethe either edge is the ancestor or descendent of the other
-- this is done via examination of teh bitvector fields of the node
isAncDescEdge :: DecoratedGraph ->  LG.LEdge EdgeInfo -> LG.LEdge EdgeInfo -> Bool
isAncDescEdge inGraph (a,_,_) (b, _, _) =
   if LG.isEmpty inGraph then error "Empty input graph in isAncDescEdge"
   else
      let aBV = bvLabel $ fromJust $ LG.lab inGraph a
          bBV = bvLabel $ fromJust $ LG.lab inGraph b
      in
      -- trace ("IADE: " ++ (show (a, aBV, b, bBV, aBV .&. bBV))) (
      if aBV .&. bBV == aBV then True
      else if aBV .&. bBV == bBV then True
      else False
      --- )


{- These heuristics do not seem tom work well at all-}

-- | heuristic add delta' based on new display tree and delta from existing costs by block--assumming < 0
-- original edges subtree1 ((u,l),(u,v)) and subtree2 ((u',v'),(u',l')) create a directed edge from
-- subtree 1 to subtree 2 via
-- 1) Add node x and y, delete edges (u,v) and (u'v') and create edges (u,x), (x,v), (u',y), and (y,v')
-- 2) real cost is the sum of block costs that are lower for new graph versus older
-- 3) heuristic is when new subtree is lower than existing block by block
--    so calculate d(u,v) + d(u',v') [existing display tree cost estimate] compared to
--    d((union u,v), v') - d(u'.v') [New display tree cost estimate] over blocks
--    so blockDelta = if d((union u,v), v') - d(u'.v') < d(u,v) + d(u',v') then d((union u,v), v') - d(u'.v')
--                     else 0 [existing better]
--    graphDelta = egdeAddDelta (separately calculated) + sum [blockDelta]
--    Compare to real delta to check behavior
-- original subtrees u -> (a,v) and u' -> (v',b)
heuristicAddDelta' :: GlobalSettings -> PhylogeneticGraph -> (LG.LEdge b, LG.LEdge b) -> VertexCost
heuristicAddDelta' _ inPhyloGraph ((u,v, _), (u',v', _)) =
  if LG.isEmpty (fst6 inPhyloGraph) then error "Empty graph in heuristicAddDelta"
  else
      let a = head $ filter (/= v) $ LG.descendants (fst6 inPhyloGraph) u
          b = head $ filter (/= v') $ LG.descendants (fst6 inPhyloGraph) u'
          blockDeltaV = V.zipWith (getBlockDelta (u,v,u',v',a,b)) (fft6 inPhyloGraph) (six6 inPhyloGraph)
      in
      V.sum blockDeltaV

-- | getBlockDelta determines the network add delta for each block (vector of characters)
-- if existing is lower then zero, else (existing - new)
getBlockDelta :: (LG.Node, LG.Node, LG.Node, LG.Node, LG.Node, LG.Node) -> V.Vector DecoratedGraph -> V.Vector CharInfo -> VertexCost
getBlockDelta (u,v,u',v',a,b) inCharV charInfoV =
   if V.null inCharV then error "Empty charcter tree vector in getBlockDelta"
   else
      let (charNewV, charExistingV) = V.unzip $ V.zipWith (getCharacterDelta (u,v,u',v',a,b)) inCharV charInfoV
          newCost = V.sum charNewV
          existingCost = V.sum charExistingV
      in
      -- trace ("GBD: " ++ (show (newCost, existingCost))) (
      if (newCost < existingCost) then newCost - existingCost
      else 0.0
      -- )

-- | getCharacterDelta determines the network add delta for each block (vector of characters)
-- if existing is lower then zero, else (existing - new)
--  calculate d(u,v) + d(u',v') [existing display tree cost estimate] compared to
--  d((union u,v), v') - d(u'.v')
-- need to use final assignemnts--so set prelim to final first
getCharacterDelta :: (LG.Node, LG.Node, LG.Node, LG.Node, LG.Node, LG.Node) -> DecoratedGraph -> CharInfo -> (VertexCost, VertexCost)
getCharacterDelta (_,v,_,v',a,b) inCharTree charInfo =
-- getCharacterDelta (u,v,u',v',a,b) inCharTree charInfo =
   let doIA = False
       -- filterGaps = True
       -- uData = V.head $ V.head $ vertData $ fromJust $ LG.lab inCharTree u
       vData = V.head $ V.head $ vertData $ fromJust $ LG.lab inCharTree v
       vFinalData = V.head $ V.head $ PRE.setPreliminaryToFinalStates  $ vertData $ fromJust $ LG.lab inCharTree v
       -- u'Data = V.head $ V.head $ vertData $ fromJust $ LG.lab inCharTree u'
       v'Data = V.head $ V.head $ vertData $ fromJust $ LG.lab inCharTree v'
       v'FinalData = V.head $ V.head $ PRE.setPreliminaryToFinalStates  $ vertData $ fromJust $ LG.lab inCharTree v'
       aData = V.head $ V.head $ vertData $ fromJust $ LG.lab inCharTree a
       aFinalData = V.head $ V.head $ PRE.setPreliminaryToFinalStates  $ vertData $ fromJust $ LG.lab inCharTree a
       bData = V.head $ V.head $ vertData $ fromJust $ LG.lab inCharTree b

       -- unionUV = M.union2Single doIA filterGaps uData vData charInfo
       -- (_,dUV) =  M.median2Single doIA uData vData charInfo
       -- dUV = vertexCost $ fromJust $ LG.lab inCharTree u
       -- dU'V' = vertexCost $ fromJust $ LG.lab inCharTree u'
       -- (_, dUnionUVV') = M.median2Single doIA unionUV v'Data charInfo

       (newX, dVV') = M.median2Single doIA vFinalData v'FinalData charInfo
       (_, dAX) = M.median2Single doIA aFinalData newX charInfo
       (_, dAV) = M.median2Single doIA aData vData charInfo
       (_, dV'B) = M.median2Single doIA v'Data bData charInfo
   in
   -- trace ("GCD: " ++ (show (dVV' + dAX, dAV + dV'B))) (
   (dVV' + dAX, dAV + dV'B)
   -- if dUnionUVV' - dU'V' < dU'V' then dUnionUVV' - dU'V'
   -- else 0.0
   -- )


-- | heuristicAddDelta takes the existing graph, edge pair, and new nodes to create and makes
-- the new nodes and reoptimizes starting nodes of two edges.  Returns cost delta based on
-- previous and new node resolution caches
-- returns cost delta and the reoptimized nodes for use in incremental optimization
-- original edges (to be deleted) (u,v) and (u',v'), n1 inserted in (u,v) and n2 inserted into (u',v')
-- creates (n1, n2), (u,n1), (n1,v), (u',n2), (n2, v')
heuristicAddDelta :: GlobalSettings 
                  -> PhylogeneticGraph 
                  -> (LG.LEdge b, LG.LEdge b) 
                  -> LG.Node 
                  -> LG.Node 
                  -> (VertexCost, LG.LNode VertexInfo, LG.LNode VertexInfo, LG.LNode VertexInfo, LG.LNode VertexInfo)
heuristicAddDelta inGS inPhyloGraph ((u,v, _), (u',v', _)) n1 n2 =
  if LG.isEmpty (fst6 inPhyloGraph) then error "Empty graph in heuristicAddDelta"
  else if graphType inGS == HardWired then
      let uvVertData = M.makeEdgeData  False True (thd6 inPhyloGraph) (six6 inPhyloGraph) (u, v, dummyEdge)
          uvPrimeData =  M.makeEdgeData  False True (thd6 inPhyloGraph) (six6 inPhyloGraph) (u', v', dummyEdge)
          hardDelta = V.sum $ fmap V.sum $ fmap (fmap snd) $ POSW.createVertexDataOverBlocks uvVertData uvPrimeData (six6 inPhyloGraph) []
      in
      (hardDelta, dummyNode, dummyNode, dummyNode, dummyNode)

  -- softwired
  else
      let uLab =      fromJust $ LG.lab (thd6 inPhyloGraph) u
          uPrimeLab = fromJust $ LG.lab (thd6 inPhyloGraph) u'
          vLab =      fromJust $ LG.lab (thd6 inPhyloGraph) v
          vPrimeLab = fromJust $ LG.lab (thd6 inPhyloGraph) v'
          uPrimeOtherChild = head $ filter ((/= v') . fst) $ LG.labDescendants (thd6 inPhyloGraph) (u', uPrimeLab)
          uOtherChild      = head $ filter ((/= v) . fst) $ LG.labDescendants (thd6 inPhyloGraph) (u, uLab)

          -- direction first edge to second so n2 is outdegree 1 to v'
          n2Lab          = NEW.getOutDegree1VertexSoftWired n2 vPrimeLab (thd6 inPhyloGraph) [n2]
          uPrimeLabAfter = NEW.getOutDegree2VertexSoftWired inGS (six6 inPhyloGraph) u' (n2, n2Lab) uPrimeOtherChild (thd6 inPhyloGraph)
          n1Lab          = NEW.getOutDegree2VertexSoftWired inGS (six6 inPhyloGraph) n1 (v, vLab) (n2, n2Lab) (thd6 inPhyloGraph)
          uLabAfter      = NEW.getOutDegree2VertexSoftWired inGS (six6 inPhyloGraph) u uOtherChild (n1, n1Lab) (thd6 inPhyloGraph)

          -- cost of resolutions
          (_, uCostBefore) = NEW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uLab)
          (_, uPrimeCostBefore) = NEW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uPrimeLab)
          (_, uCostAfter) = NEW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uLabAfter)
          (_, uPrimeCostAfter) = NEW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uPrimeLabAfter)

          addNetDelta = (uCostAfter - uCostBefore) +  (uPrimeCostAfter - uPrimeCostBefore)


      in
      -- trace ("HAD: " ++ (show (uCostAfter, uCostBefore, uPrimeCostAfter, uPrimeCostBefore)) ++ " -> " ++ (show addNetDelta)) $
      if null (filter ((/= v') . fst) $ LG.labDescendants (thd6 inPhyloGraph) (u', uPrimeLab)) || null (filter ((/= v) . fst) $ LG.labDescendants (thd6 inPhyloGraph) (u, uLab)) then (infinity, dummyNode, dummyNode, dummyNode, dummyNode)
      -- this should not happen--should try to create new edges from children of net edges
      else if (length $ LG.descendants (thd6 inPhyloGraph) u) < 2 ||  (length $ LG.descendants (thd6 inPhyloGraph) u') < 2 then error ("Outdegree 1 nodes in heuristicAddDelta")
      else
         (addNetDelta, (u, uLabAfter), (u', uPrimeLabAfter), (n1, n1Lab), (n2, n2Lab))
      


-- | deltaPenaltyAdjustment takes number of leaves and Phylogenetic graph and returns a heuristic graph penalty for adding a single network edge
-- if Wheeler2015Network, this is based on all changes affecting a single block (most permissive) and Wheeler 2015 calculation of penalty
-- if PMDLGraph -- KMDL not yet implemented
-- if NoNetworkPenalty then 0
-- modification "add" or subtrct to calculate delta
-- always delta is positive--whether neg or pos is deltermined when used
deltaPenaltyAdjustment :: GlobalSettings 
                       -> PhylogeneticGraph 
                       -> String 
                       -> VertexCost
deltaPenaltyAdjustment inGS inGraph modification =
   -- trace ("DPA: entering: " ++ (show $ graphFactor inGS)) (
   let numLeaves = numDataLeaves inGS
       edgeCostModel = graphFactor inGS
       (_, _, _, networkNodeList) = LG.splitVertexList (fst6 inGraph)
   in
   if edgeCostModel == NoNetworkPenalty then
      -- trace ("DPA: No penalty")
      0.0

   -- else if length networkNodeList == 0 then
   -- trace ("DPA: No cost")
   --   0.0

   else if edgeCostModel == Wheeler2015Network then
      -- trace  ("DPW: In Wheeler2015Network") (
      let graphCost = snd6 inGraph -- this includes any existing penalties--would be better not to include
          -- numBlocks = V.length $ fth6 inGraph
      in
      -- if (graphType inGS) == HardWired then 0.0
      -- trace ("DPA Value: " ++ (show $ graphCost / (fromIntegral $ numBlocks * 2 * ((2 * numLeaves) - 2))))
      -- else
      -- graphCost / (fromIntegral $ numBlocks * 2 * ((2 * numLeaves) - 2))
      -- this one--1/2 total edge number for factor
      graphCost / (fromIntegral $ 2 * ((2 * numLeaves) - 2) + (2 * (length networkNodeList)))
      -- )

   else if edgeCostModel == PMDLGraph then
      -- trace  ("DPW: In PMDLGraph") (
      if graphType inGS == Tree then
         fst $ IL.head (graphComplexityList inGS)

      else if graphType inGS == SoftWired then
         let currentComplexity = fst $ (graphComplexityList inGS) IL.!!! (length networkNodeList)
             nextComplexity    = if modification == "add" then fst $ (graphComplexityList inGS) IL.!!! ((length networkNodeList) + 1)
                                 else if modification == "delete" then fst $ (graphComplexityList inGS) IL.!!!  ((length networkNodeList) - 1)
                                 else error ("SoftWired deltaPenaltyAdjustment modification not recognized: " ++ modification)
         in
         abs (currentComplexity - nextComplexity)

      else if graphType inGS == HardWired then
         let currentComplexity = snd $ (graphComplexityList inGS) IL.!!! (length networkNodeList)
             nextComplexity    = if modification == "add" then snd $ (graphComplexityList inGS) IL.!!! ((length networkNodeList) + 1)
                                 else if modification == "delete" then snd $ (graphComplexityList inGS) IL.!!!  ((length networkNodeList) - 1)
                                 else error ("HardWired deltaPenaltyAdjustment modification not recognized: " ++ modification)
         in
         abs (currentComplexity - nextComplexity)


      else error ("Graph type not yet implemented: " ++ (show $ graphType inGS))
      -- )

   else error ("Network edge cost model not yet implemented: " ++ (show edgeCostModel))
   -- )



-- | deleteEachNetEdge takes a phylogenetic graph and deletes all network edges one at time
-- and returns best list of new Phylogenetic Graphs and cost
-- even if worse--could be used for simulated annealing later
-- if equal returns unique graph list
deleteEachNetEdge :: GlobalSettings
                  -> ProcessedData
                  -> Int
                  -> Int
                  -> Bool
                  -> Bool
                  -> Bool
                  -> Maybe SAParams
                  -> PhylogeneticGraph
                  -> ([PhylogeneticGraph], VertexCost, Maybe SAParams)
deleteEachNetEdge inGS inData rSeed numToKeep doSteepest doRandomOrder force inSimAnnealParams inPhyloGraph =
   -- trace ("DENE start") (
   if LG.isEmpty $ thd6 inPhyloGraph then ([], infinity, inSimAnnealParams) -- error "Empty input phylogenetic graph in deleteAllNetEdges"
   else
      let currentCost = snd6 inPhyloGraph

          -- potentially randomize order of list
          networkEdgeList' = LG.netEdges $ thd6 inPhyloGraph
          networkEdgeList = if not doRandomOrder then networkEdgeList'
                            else permuteList rSeed networkEdgeList'


          (newGraphList, newSAParams) = if not doSteepest then
                                          (PU.seqParMap (parStrategy $ lazyParStrat inGS) (deleteNetEdge inGS inData inPhyloGraph force) networkEdgeList, U.incrementSimAnnealParams inSimAnnealParams)
                                        else
                                          deleteNetEdgeRecursive inGS inData inPhyloGraph force inSimAnnealParams networkEdgeList

          bestCostGraphList = filter ((/= infinity) . snd6) $ GO.selectGraphs Best numToKeep 0.0 (-1) newGraphList
          minCost = if null bestCostGraphList then infinity
                    else minimum $ fmap snd6 bestCostGraphList
      in
      -- no network edges to delete
      if null networkEdgeList then trace ("\tNo network edges to delete") ([inPhyloGraph], currentCost, inSimAnnealParams)

      else
         -- single if steepest so no neeed to unique--and have run through all options (including SA stuff) via recursive call
         if doSteepest then
            (newGraphList, minCost, newSAParams)

      -- "all" option needs to recurse since does all available edges at each step
      -- logic is here since not in the deleteNetEdge function
      else if isNothing inSimAnnealParams then
         if minCost < currentCost then
            -- trace ("DENE--Delete net edge return:" ++ (show (minCost,length uniqueCostGraphList))) (
            let newRandIntList = take (length bestCostGraphList) (randomIntList rSeed)
                annealParamList = U.generateUniqueRandList (length bestCostGraphList) newSAParams
                nextGraphTripleList = PU.seqParMap (parStrategy $ lazyParStrat inGS) (deleteEachNetEdge' inGS inData numToKeep doSteepest doRandomOrder force) (zip3 annealParamList newRandIntList bestCostGraphList)

                newMinCost = minimum $ fmap snd3 nextGraphTripleList
                newGraphListBetter = filter ((== newMinCost) . snd6) $ concatMap fst3 nextGraphTripleList
            in
            (GO.selectGraphs Unique numToKeep 0.0 (-1) $ newGraphListBetter, newMinCost, newSAParams)

         else
            (bestCostGraphList, currentCost, newSAParams)

      -- SA anneal/Drift
      else
            -- always take better
            if minCost < currentCost then
               (bestCostGraphList, minCost, newSAParams)

            -- check if hit step limit--more for SA than drift
            else if ((currentStep $ fromJust inSimAnnealParams) >= (numberSteps $ fromJust inSimAnnealParams)) || ((driftChanges $ fromJust inSimAnnealParams) >= (driftMaxChanges $ fromJust inSimAnnealParams)) then
                  ([inPhyloGraph], snd6 inPhyloGraph, inSimAnnealParams)

            -- otherwise do the anneal/Drift accept, or keep going on input graph
            else
               let (acceptGraph, nextSAParams) = U.simAnnealAccept inSimAnnealParams currentCost minCost
               in
               if acceptGraph then
                  (bestCostGraphList, minCost, newSAParams)
               else
                  deleteEachNetEdge inGS inData (head $ randomIntList rSeed) numToKeep doSteepest doRandomOrder force nextSAParams inPhyloGraph

-- | deleteEachNetEdge' is a wrapper around deleteEachNetEdge to allow for zipping new random seeds for each
-- replicate
deleteEachNetEdge' :: GlobalSettings
                  -> ProcessedData
                  -> Int
                  -> Bool
                  -> Bool
                  -> Bool
                  -> (Maybe SAParams, Int, PhylogeneticGraph)
                  -> ([PhylogeneticGraph], VertexCost, Maybe SAParams)
deleteEachNetEdge' inGS inData numToKeep doSteepest doRandomOrder force (inSimAnnealParams, rSeed, inPhyloGraph) =
   deleteEachNetEdge inGS inData rSeed numToKeep doSteepest doRandomOrder force inSimAnnealParams inPhyloGraph


-- | deleteNetEdgeRecursive like deleteEdge, deletes an edge (checking if network) and rediagnoses graph
-- contacts in=out=1 edges and removes node, reindexing nodes and edges
-- except returns on first better (as opposed to do all deletes first)
-- or sim annleal/drift
deleteNetEdgeRecursive :: GlobalSettings
                       -> ProcessedData
                       -> PhylogeneticGraph
                       -> Bool
                       -> Maybe SAParams
                       -> [LG.Edge]
                       -> ([PhylogeneticGraph], Maybe SAParams)
deleteNetEdgeRecursive inGS inData inPhyloGraph force inSimAnnealParams inEdgeToDeleteList =
   if null inEdgeToDeleteList then ([], inSimAnnealParams)
   else
       let {- Unclear if should adjust to number of rounds if already limiting to graphsSteepest value
            saRounds = if isNothing inSimAnnealParams then 1
                      else rounds $ fromJust inSimAnnealParams

            (numGraphsToExamine, _) = divMod PU.getNumThreads saRounds -- this may not "drift" if finds alot better, but that's how its supposed to work
           -}
           numGraphsToExamine = min (graphsSteepest inGS) PU.getNumThreads
           -- edgeToDelete = head inEdgeToDeleteList
           edgeToDeleteList = take numGraphsToExamine inEdgeToDeleteList

           -- calls general funtion to remove network graph edge
           -- (delSimple, wasModified) = deleteNetworkEdge (fst6 inPhyloGraph) edgeToDelete
           simpleGraphList = fmap fst $ filter ((== True) . snd) $ PU.seqParMap (parStrategy $ lazyParStrat inGS) (deleteNetworkEdge (fst6 inPhyloGraph)) edgeToDeleteList

           -- delSimple = GO.contractIn1Out1EdgesRename $ LG.delEdge edgeToDelete $ fst6 inPhyloGraph

           -- prune other edges if now unused
           pruneEdges = False

           -- don't warn that edges are being pruned
           warnPruneEdges = False

           -- graph optimization from root
           startVertex = Nothing

           -- (heuristicDelta, _, _) = heuristicDeleteDelta inGS inPhyloGraph edgeToDelete
           -- heuristicDelta = 0.0

           -- can treat as negative for delete
           -- edgeAddDelta = deltaPenaltyAdjustment inGS inPhyloGraph "delete"

           -- full two-pass optimization
           leafGraph = LG.extractLeafGraph $ thd6 inPhyloGraph
                                 
           newPhyloGraphList' = if (graphType inGS == SoftWired) then
                                  PU.seqParMap (parStrategy $ lazyParStrat inGS) (T.multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex) simpleGraphList
                                 else if (graphType inGS == HardWired) then
                                    PU.seqParMap (parStrategy $ lazyParStrat inGS) (T.multiTraverseFullyLabelHardWired inGS inData leafGraph startVertex) simpleGraphList
                                 else error "Unsupported graph type in deleteNetEdge.  Must be soft or hard wired"

           newPhyloGraphList = GO.selectGraphs Best (maxBound::Int) 0.0 (-1) newPhyloGraphList'

       in
       -- if not modified return original graph
       -- This check seems to be issue with delete not functioning properly
       if null simpleGraphList then ([inPhyloGraph], inSimAnnealParams)

       -- forcing delete for move
       else if force then
         -- trace ("DNERec forced")
         (newPhyloGraphList, inSimAnnealParams)

       -- regular search not sim anneal/drift
       else if (isNothing inSimAnnealParams) then

          --return if better
          if (snd6 $ head newPhyloGraphList) < (snd6 inPhyloGraph) then
            -- trace  ("DNERec Better -> " ++ (show $ snd6 newPhyloGraph))
            (newPhyloGraphList, inSimAnnealParams)

          else
            -- need to update edge list for new graph
            -- potentially randomize order of list
            deleteNetEdgeRecursive inGS inData inPhyloGraph force inSimAnnealParams (drop numGraphsToExamine inEdgeToDeleteList)

      -- sim anneal/drift
      else

         -- if better always accept
         if (snd6 $ head newPhyloGraphList) < (snd6 inPhyloGraph) then
            -- these graph costs are "exact" or at least non-heuristic--needs to be updated when get a good heuristic
            let (_, nextSAParams) = U.simAnnealAccept inSimAnnealParams (snd6 inPhyloGraph) (snd6 $ head newPhyloGraphList)
            in
            (newPhyloGraphList, nextSAParams)


         -- check if hit step limit--more for SA than drift
         else if ((currentStep $ fromJust inSimAnnealParams) >= (numberSteps $ fromJust inSimAnnealParams)) || ((driftChanges $ fromJust inSimAnnealParams) >= (driftMaxChanges $ fromJust inSimAnnealParams)) then
                  ([inPhyloGraph], inSimAnnealParams)


          -- otherwise do the anneal/Drift accept
         else
            let (acceptGraph, nextSAParams) = U.simAnnealAccept inSimAnnealParams (snd6 inPhyloGraph) (snd6 $ head newPhyloGraphList)
            in
            if acceptGraph then
                (newPhyloGraphList, nextSAParams)

            else
               deleteNetEdgeRecursive inGS inData inPhyloGraph force nextSAParams (drop numGraphsToExamine inEdgeToDeleteList)

-- | deleteEdge deletes an edge (checking if network) and rediagnoses graph
-- contacts in=out=1 edgfes and removes node, reindexing nodes and edges
-- naive for now
-- force requires reoptimization no matter what--used for net move
-- skipping heuristics for now--awful
-- calls deleteNetworkEdge that has various graph checks
deleteNetEdge :: GlobalSettings 
              -> ProcessedData 
              -> PhylogeneticGraph 
              -> Bool 
              -> LG.Edge 
              -> PhylogeneticGraph
deleteNetEdge inGS inData inPhyloGraph force edgeToDelete =
   if LG.isEmpty $ thd6 inPhyloGraph then error "Empty input phylogenetic graph in deleteNetEdge"
   else if not (LG.isNetworkEdge (fst6 inPhyloGraph) edgeToDelete) then error ("Edge to delete: " ++ (show edgeToDelete) ++ " not in graph:\n" ++ (LG.prettify $ fst6 inPhyloGraph))
   else
       -- trace ("DNE: " ++ (show edgeToDelete)) (
       let (delSimple, wasModified)  = deleteNetworkEdge (fst6 inPhyloGraph) edgeToDelete

           -- delSimple = GO.contractIn1Out1EdgesRename $ LG.delEdge edgeToDelete $ fst6 inPhyloGraph

           -- prune other edges if now unused
           pruneEdges = False

           -- don't warn that edges are being pruned
           warnPruneEdges = False

           -- graph optimization from root
           startVertex = Nothing

           -- (heuristicDelta, _, _) = heuristicDeleteDelta inGS inPhyloGraph edgeToDelete


           -- edgeAddDelta = deltaPenaltyAdjustment inGS inPhyloGraph "delete"


           -- full two-pass optimization--cycles checked in edge deletion function
           leafGraph = LG.extractLeafGraph $ thd6 inPhyloGraph
                                 
           newPhyloGraph = if (graphType inGS == SoftWired) then
                              T.multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex delSimple
                           else if (graphType inGS == HardWired) then
                              T.multiTraverseFullyLabelHardWired inGS inData leafGraph startVertex delSimple
                           else error "Unsupported graph type in deleteNetEdge.  Must be soft or hard wired"
       in
       --check if deletino modified graph
       if not wasModified then inPhyloGraph

       -- else if force || (graphType inGS) == HardWired then
       else if force then
         -- trace ("DNE forced")
         newPhyloGraph
       else -- if (heuristicDelta / (dynamicEpsilon inGS)) - edgeAddDelta < 0 then newPhyloGraph
         if (snd6 newPhyloGraph) < (snd6 inPhyloGraph) then
            -- trace ("DNE Better: " ++ (show $ snd6 newPhyloGraph))
            newPhyloGraph
         else
            -- trace ("DNE Not Better: " ++ (show $ snd6 newPhyloGraph))
            inPhyloGraph
      -- )

-- | deleteNetworkEdge deletes a network edges from a simple graph
-- retuns newGraph if can be modified or input graph with Boolean to tell if modified
-- and contracts, reindexes/names internaledges/veritices around deletion
-- can't raise to general graph level due to vertex info
-- in edges (b,a) (c,a) (a,d), deleting (a,b) deletes node a, inserts edge (b,d)
-- contacts node c since  now in1out1 vertex
-- checks for chained network edges--can be created by progressive deletion
-- checks for cycles now
-- shouldn't need for check for creating a node with children that are both network nodes
-- since that would require that condition coming in and shodl be there--ie checked earlier in addition and input
deleteNetworkEdge :: SimpleGraph -> LG.Edge -> (SimpleGraph, Bool)
deleteNetworkEdge inGraph inEdge@(p1, nodeToDelete) =
   if LG.isEmpty inGraph then error ("Cannot delete edge from empty graph")
   else
      let childrenNodeToDelete = LG.descendants inGraph nodeToDelete
          parentsNodeToDelete = LG.parents inGraph nodeToDelete
          --parentNodeToKeep = head $ filter (/= p1) parentsNodeToDelete
          --newEdge = (parentNodeToKeep, head childrenNodeToDelete, 0.0)
          -- newGraph = LG.insEdge newEdge $ LG.delNode nodeToDelete inGraph
          newGraph = LG.delEdge inEdge inGraph
          -- newGraph' = GO.contractIn1Out1EdgesRename newGraph

          -- conversion as if input--see if affects length
          -- newGraph'' = GO.convertGeneralGraphToPhylogeneticGraph False newGraph
          newGraph'' = GO.contractIn1Out1EdgesRename newGraph

      in
      -- error conditions and creation of chained network edges (forbidden in phylogenetic graph--causes resolutoin cache issues)
      if length childrenNodeToDelete /= 1 then error ("Cannot delete non-network edge in deleteNetworkEdge: (1)" ++ (show inEdge) ++ "\n" ++ (LG.prettyIndices inGraph))
      else if length parentsNodeToDelete /= 2 then error ("Cannot delete non-network edge in deleteNetworkEdge (2): " ++ (show inEdge) ++ "\n" ++ (LG.prettyIndices inGraph))

      -- warning if chained on input, skip if chained net edges in output
      else if (LG.isNetworkNode inGraph p1) then
         -- error ("Error: Chained network nodes in deleteNetworkEdge : " ++ (show inEdge) ++ "\n" ++ (LG.prettyIndices inGraph) ++ " skipping")
         trace ("\tWarning: Chained network nodes in deleteNetworkEdge skipping deletion")
         (LG.empty, False)

      else if LG.hasChainedNetworkNodes newGraph'' then
         trace ("\tWarning: Chained network nodes in deleteNetworkEdge skipping deletion (2)")
         (LG.empty, False)

      else if LG.isEmpty newGraph'' then (LG.empty, False)

      else
         {-trace ("DNE: Edge to delete " ++ (show inEdge) ++ " cnd " ++ (show childrenNodeToDelete) ++ " pnd " ++ (show parentsNodeToDelete) ++ " pntk " ++ (show parentNodeToKeep)
            ++ " ne " ++ (show newEdge) ++ "\nInGraph: " ++ (LG.prettyIndices inGraph) ++ "\nNewGraph: " ++ (LG.prettyIndices newGraph) ++ "\nNewNewGraph: "
            ++ (LG.prettyIndices newGraph')) -}
         (newGraph'', True)


-- | heuristicDeleteDelta takes the existing graph, edge to delete,
-- reoptimizes starting nodes of two created edges.  Returns cost delta based on
-- previous and new node resolution caches
-- delete n1 -> n2, create u -> v, u' -> v'
-- assumes original is edge n1 -> n2, u' -> (n2, X), n1 -> (n2,v), u (n1,Y)
heuristicDeleteDelta :: GlobalSettings 
                     -> PhylogeneticGraph 
                     -> LG.Edge 
                     -> (VertexCost, LG.LNode VertexInfo, LG.LNode VertexInfo)
heuristicDeleteDelta inGS inPhyloGraph (n1, n2) =
  if LG.isEmpty (fst6 inPhyloGraph) then error "Empty graph in heuristicDeleteDelta"
  else if graphType inGS == HardWired then
      -- ensures delete--will always be lower or equakl cost if delete edge from HardWired
      (-1, dummyNode, dummyNode)
  else
      let inGraph = thd6 inPhyloGraph
          u  = head $ LG.parents inGraph n1
          u' = head $ filter (/= n1) $ LG.parents inGraph n2
          v' = head $ LG.descendants inGraph n2
          v  = head $ filter (/= n2) $ LG.descendants inGraph n1

          uLab      = fromJust $ LG.lab inGraph u
          uPrimeLab = fromJust $ LG.lab inGraph u'
          vLab =      fromJust $ LG.lab inGraph v
          vPrimeLab = fromJust $ LG.lab inGraph v'

          uOtherChild      = head $ filter ((/= n1) . fst) $ LG.labDescendants inGraph (u, uLab)
          uPrimeOtherChild = head $ filter ((/= n2) . fst) $ LG.labDescendants inGraph (u', uPrimeLab)

          -- skip over netnodes
          uLabAfter      = NEW.getOutDegree2VertexSoftWired inGS (six6 inPhyloGraph) u (v, vLab) uOtherChild inGraph
          uPrimeLabAfter = NEW.getOutDegree2VertexSoftWired inGS (six6 inPhyloGraph) u' (v', vPrimeLab) uPrimeOtherChild inGraph

          -- cost of resolutions
          (_, uCostBefore) = NEW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uLab)
          (_, uPrimeCostBefore) = NEW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uPrimeLab)
          (_, uCostAfter) = NEW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uLabAfter)
          (_, uPrimeCostAfter) = NEW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uPrimeLabAfter)

          addNetDelta = uCostAfter - uCostBefore +  uPrimeCostAfter - uPrimeCostBefore


      in
      -- this should not happen--should try to crete new edges from children of net edges
      if null (LG.parents inGraph n1) || null (filter (/= n1) $ LG.parents inGraph n2) || null (LG.descendants inGraph n2) || null (filter (/= n2) $ LG.descendants inGraph n1) || null (filter ((/= n2) . fst) $ LG.labDescendants inGraph (u', uPrimeLab)) || null (filter ((/= n1) . fst) $ LG.labDescendants inGraph (u, uLab)) then (infinity, dummyNode, dummyNode)
      -- this should not happen--should try to crete new edges from children of net edges
      else if (length (LG.parents inGraph n1) /= 1) || (length (LG.parents inGraph n2) /= 2) || (length (LG.descendants inGraph n2) /= 1) || (length (LG.descendants inGraph n1) /= 2) then error ("Graph malformation in numbers of parents and children in heuristicDeleteDelta")
      else
         (addNetDelta, (u, uLabAfter), (u', uPrimeLabAfter))

{-
-- | insertNetEdgeBothDirections calls insertNetEdge for both u -> v and v -> u new edge orientations
insertNetEdgeBothDirections :: GlobalSettings -> ProcessedData -> PhylogeneticGraph ->  Maybe VertexCost -> (LG.LEdge b, LG.LEdge b) -> [PhylogeneticGraph]
insertNetEdgeBothDirections inGS inData inPhyloGraph preDeleteCost (u,v) = fmap (insertNetEdge inGS inData inPhyloGraph preDeleteCost) [(u,v), (v,u)]
-}
