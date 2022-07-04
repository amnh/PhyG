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
                                ) where

import           Control.Parallel.Strategies
import           Data.Bits
import           Data.Maybe
import qualified Data.Text.Lazy                       as TL
import qualified Data.Vector                          as V
import           Debug.Trace
import           GeneralUtilities
import qualified GraphOptimization.Medians            as M
import qualified GraphOptimization.PostOrderFunctions as POS
import qualified GraphOptimization.Traversals         as T
import qualified Graphs.GraphOperations               as GO
import qualified ParallelUtilities                    as PU
import           Types.Types
import qualified Utilities.LocalGraph                 as LG
import qualified Utilities.Utilities                  as U
import qualified Data.InfList                         as IL
import qualified GraphOptimization.PostOrderSoftWiredFunctions as POSW
import qualified GraphOptimization.PreOrderFunctions as PRE
-- import qualified Data.List                            as L

-- Need a "steepest" that takes first better netowrk for add and delete.
-- Choose order in max branch length, root-end, leaf-end, and at random

-- | moveAllNetEdges is a wrapper for moveAllNetEdges' allowing for multiple simulated annealing rounds
moveAllNetEdges :: GlobalSettings -> ProcessedData -> Int -> Int -> Int -> Int -> Bool -> Bool -> Bool -> ([PhylogeneticGraph], VertexCost) -> (Maybe SAParams, [PhylogeneticGraph]) -> ([PhylogeneticGraph], Int)
moveAllNetEdges inGS inData rSeed maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
   if inSimAnnealParams == Nothing then moveAllNetEdges' inGS inData rSeed maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) Nothing inPhyloGraphList
   else
      let -- create list of params with unique list of random values for rounds of annealing
         annealingRounds = rounds $ fromJust inSimAnnealParams
         saPAramList = (U.generateUniqueRandList annealingRounds inSimAnnealParams) -- (replicate annealingRounds inPhyloGraphList)

         (annealRoundsList, counterList) = unzip (zipWith (moveAllNetEdges' inGS inData rSeed maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost)) saPAramList (replicate annealingRounds inPhyloGraphList)`using` PU.myParListChunkRDS)
      in
      (take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] (concat annealRoundsList) , sum counterList)

-- | moveAllNetEdges' removes each edge and adds an edge to all possible places (or steepest) each round
-- until no better or additional graphs are found
-- call with ([], infinity) [single input graph]
moveAllNetEdges' :: GlobalSettings -> ProcessedData -> Int -> Int -> Int -> Int -> Bool -> Bool -> Bool -> ([PhylogeneticGraph], VertexCost) -> Maybe SAParams -> [PhylogeneticGraph] -> ([PhylogeneticGraph], Int)
moveAllNetEdges' inGS inData rSeed maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) inSimAnnealParams inPhyloGraphList =
   if null inPhyloGraphList then (take numToKeep curBestGraphList, counter)
   else
      let firstPhyloGraph = head inPhyloGraphList
          leafGraph = LG.extractLeafGraph $ thd6 firstPhyloGraph
          currentCost = min curBestGraphCost (snd6 firstPhyloGraph)
          netEdgeList = LG.labNetEdges (thd6 firstPhyloGraph)
          -- newGraphList' = concat (zipWith3 (deleteOneNetAddAll inGS inData leafGraph maxNetEdges numToKeep doSteepest doRandomOrder firstPhyloGraph) (randomIntList rSeed) (U.generateUniqueRandList (length netEdgeList) inSimAnnealParams) (fmap LG.toEdge netEdgeList) `using` PU.myParListChunkRDS)
          newGraphList' =  deleteOneNetAddAll inGS inData leafGraph maxNetEdges numToKeep doSteepest doRandomOrder firstPhyloGraph (fmap LG.toEdge netEdgeList) rSeed inSimAnnealParams
          newGraphList = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] newGraphList'
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
      else if inSimAnnealParams == Nothing then
         if newGraphCost > currentCost then 
            -- trace ("\t MANE : Worse") 
            moveAllNetEdges' inGS inData maxNetEdges rSeed numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (firstPhyloGraph : curBestGraphList, currentCost) inSimAnnealParams (tail inPhyloGraphList)
         else if newGraphCost < currentCost then
            trace ("\t-> " ++ (show newGraphCost))
            moveAllNetEdges' inGS inData rSeed maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (newGraphList, newGraphCost) inSimAnnealParams (newGraphList ++ (tail inPhyloGraphList))

         else
            -- new graph list contains the input graph if equal and filterd unique already in moveAllNetEdges
            let newCurSameBestList = take numToKeep $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (curBestGraphList ++ newGraphList)
            in
            -- trace ("\t MANE : Equal") 
            moveAllNetEdges' inGS inData rSeed maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (newCurSameBestList, currentCost) inSimAnnealParams (tail inPhyloGraphList)

      -- sim anneal choice
      else
            let -- abstract stopping criterion to continue
               numDone = if (method $ fromJust inSimAnnealParams) == SimAnneal then currentStep $ fromJust inSimAnnealParams
                         else driftChanges $ fromJust inSimAnnealParams
               numMax  = if (method $ fromJust inSimAnnealParams) == SimAnneal then numberSteps $ fromJust inSimAnnealParams
                         else driftMaxChanges $ fromJust inSimAnnealParams

               -- get acceptance based on heuristic costs
               uniqueGraphList = take numToKeep $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] newGraphList'
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
                   bestMoveList = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] bestMoveList'
               in
               --trace ("BM: " ++ (show $ snd6 $ head  bestMoveList))
               (take numToKeep bestMoveList, counter')
      -- )


-- | (curBestGraphList, annealBestCost) is a wrapper for moveAllNetEdges' allowing for multiple simulated annealing rounds
insertAllNetEdges :: GlobalSettings -> ProcessedData -> Int -> Int -> Int -> Int -> Bool -> Bool -> Bool -> ([PhylogeneticGraph], VertexCost) -> (Maybe SAParams, [PhylogeneticGraph]) -> ([PhylogeneticGraph], Int)
insertAllNetEdges inGS inData rSeed maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
   let leafGraph = LG.extractLeafGraph $ (thd6 . head) inPhyloGraphList
   in
   if inSimAnnealParams == Nothing then
      insertAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (randomIntList rSeed) Nothing inPhyloGraphList
   else
      let -- create list of params with unique list of random values for rounds of annealing
          annealingRounds = rounds $ fromJust inSimAnnealParams
          annealParamGraphList = U.generateUniqueRandList annealingRounds inSimAnnealParams
          replicateRandIntList = fmap randomIntList (take annealingRounds (randomIntList rSeed))

          (annealRoundsList, counterList) = unzip (zipWith3 (insertAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost)) replicateRandIntList annealParamGraphList (replicate annealingRounds inPhyloGraphList) `using` PU.myParListChunkRDS)
      in
      (take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] (concat annealRoundsList) , sum counterList)

-- | insertAllNetEdges' adds network edges one each each round until no better or additional
-- graphs are found
-- call with ([], infinity) [single input graph]
insertAllNetEdges' :: GlobalSettings -> ProcessedData -> DecoratedGraph -> Int -> Int -> Int -> Bool -> Bool -> Bool -> ([PhylogeneticGraph], VertexCost) ->  [Int] -> Maybe SAParams -> [PhylogeneticGraph] -> ([PhylogeneticGraph], Int)
insertAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) randIntList inSimAnnealParams inPhyloGraphList =
   if null inPhyloGraphList then (take numToKeep curBestGraphList, counter)
   else
      let currentCost = min curBestGraphCost (snd6 $ head inPhyloGraphList)

           --check for max net edges
          (_, _, _, netNodes) = LG.splitVertexList (thd6 $ head inPhyloGraphList)

          (newGraphList, _) = insertEachNetEdge inGS inData leafGraph (head randIntList) maxNetEdges numToKeep doSteepest doRandomOrder Nothing inSimAnnealParams (head inPhyloGraphList)


          bestNewGraphList = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] newGraphList
          newGraphCost = if (not . null) bestNewGraphList then snd6 $ head bestNewGraphList
                         else infinity

          -- this to deal with potential tail of empty list
          nextNewGraphList = if (not . null) newGraphList then tail newGraphList
                             else []

      in
      -- trace ("IANE: " ++ (show $ length netNodes)) (
      if length netNodes >= maxNetEdges then
         trace ("Maximum number of network edges reached: " ++ (show $ length netNodes))
         (take numToKeep curBestGraphList, counter)


      else if null newGraphList then (take numToKeep curBestGraphList, counter)

      -- regular insert keeping best
      else if inSimAnnealParams == Nothing then
         -- "steepest style descent" abandons existing list if better cost found
         if newGraphCost < currentCost then
            -- check if graph OK--done in insert function
            let --- isCyclicList = filter (== True) $ fmap LG.cyclic $ fmap thd6 newGraphList
                --- hasDupEdges = filter (== True) $ fmap LG.hasDuplicateEdge $ fmap thd6 newGraphList

                graphsToInsert = if doSteepest then newGraphList
                                 else take numToKeep $ newGraphList ++ (tail inPhyloGraphList)
            in
            trace ("\t-> " ++ (show newGraphCost)) 
               insertAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (newGraphList, newGraphCost) (tail randIntList) inSimAnnealParams graphsToInsert
            {-
            if (null isCyclicList && null hasDupEdges) then 
               trace ("\t-> " ++ (show newGraphCost)) 
               insertAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (newGraphList, newGraphCost) (tail randIntList) inSimAnnealParams graphsToInsert
            else 
               trace ("Cycle " ++ (show newGraphCost) ++ " Duplicate Edges " ++ (show $ hasDupEdges)) 
               insertAllNetEdges' inGS inData maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (curBestGraphList, currentCost) (tail randIntList) inSimAnnealParams (tail inPhyloGraphList)
            -}

         -- worse graphs found--go on
         else if  newGraphCost > currentCost then
            -- trace ("IANE: Worse")
            insertAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (curBestGraphList, currentCost) (tail randIntList) inSimAnnealParams (tail inPhyloGraphList)

         -- equal cost
         -- not sure if should add new graphs to queue to do edge deletion again
         else
            -- new graph list contains the input graph if equal and filterd unique already in insertAllNetEdges
            let newCurSameBestList =  GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"]  $ take numToKeep (curBestGraphList ++ newGraphList)
            in
            -- trace ("IANE: same " ++ (show $ length (tail inPhyloGraphList)))
            insertAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (newCurSameBestList, currentCost) (tail randIntList) inSimAnnealParams (tail inPhyloGraphList)
            

      -- simulated annealing
      else
            -- if steepest -- simAnneal/Drift stuff done during net add function so jusyt take results and move on
            -- create new
            if doSteepest then
               let annealBestCost = min curBestGraphCost newGraphCost
                   newSAParams = Just $ (fromJust (U.incrementSimAnnealParams inSimAnnealParams)) {currentStep = 0, driftChanges = 0}
               in
               insertAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder ((head newGraphList) :  curBestGraphList, annealBestCost) (tail randIntList) newSAParams (nextNewGraphList ++ (tail inPhyloGraphList))

            -- simmAnneal on "all list"
            else
               -- sim anneal choice
               let numDone = if (method $ fromJust inSimAnnealParams) == SimAnneal then currentStep $ fromJust inSimAnnealParams
                            else driftChanges $ fromJust inSimAnnealParams
                   numMax  = if (method $ fromJust inSimAnnealParams) == SimAnneal then numberSteps $ fromJust inSimAnnealParams
                            else driftMaxChanges $ fromJust inSimAnnealParams
                   annealBestCost = min curBestGraphCost (snd6 $ head newGraphList)
                   (acceptFirstGraph, newSAParams) = U.simAnnealAccept inSimAnnealParams annealBestCost (snd6 $ head newGraphList)
               in
               -- trace ("ACG" ++ (show acceptFirstGraph) ++ " " ++ (show $ snd6 $ head uniqueGraphList)) (
               if (numDone < numMax) then
                  -- this fixes tail fail
                  if acceptFirstGraph then
                     insertAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder ((head newGraphList) :  curBestGraphList, annealBestCost) (tail randIntList) newSAParams (nextNewGraphList ++ (tail inPhyloGraphList))
                  else
                     insertAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (curBestGraphList, annealBestCost) (tail randIntList) newSAParams (nextNewGraphList ++ (tail inPhyloGraphList))

               -- returns non-optimized list for GA or whatever
               else if returnMutated then (take numToKeep curBestGraphList, counter)

               -- run net delete regular to get back to optimized edges
               else
                  let (bestList', counter') =  deleteAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep (counter + 1) False doSteepest doRandomOrder ([], annealBestCost) (tail randIntList) Nothing (take numToKeep curBestGraphList)
                      bestList = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] bestList'
               in
               --trace ("BM: " ++ (show $ snd6 $ head  bestMoveList))
               (take numToKeep bestList, counter')
               -- )


-- | insertEachNetEdge takes a phylogenetic graph and inserts all permissible network edges one at time
-- and returns unique list of new Phylogenetic Graphs and cost
-- even if worse--could be used for simulated annealing later
-- if equal returns unique graph list
insertEachNetEdge :: GlobalSettings 
                  -> ProcessedData 
                  -> DecoratedGraph 
                  -> Int 
                  -> Int 
                  -> Int 
                  -> Bool 
                  -> Bool 
                  -> Maybe VertexCost 
                  -> Maybe SAParams 
                  -> PhylogeneticGraph 
                  -> ([PhylogeneticGraph], VertexCost)
insertEachNetEdge inGS inData leafGraph  rSeed maxNetEdges numToKeep doSteepest doRandomOrder preDeleteCost inSimAnnealParams inPhyloGraph =
   if LG.isEmpty $ fst6 inPhyloGraph then error "Empty input insertEachNetEdge graph in deleteAllNetEdges"
   else
      let currentCost = if preDeleteCost == Nothing then snd6 inPhyloGraph
                        else fromJust preDeleteCost

          candidateNetworkEdgeList' = getPermissibleEdgePairs (thd6 inPhyloGraph)

          -- radomize pair list
          rSeedList = randomIntList rSeed
          candidateNetworkEdgeList = if not doSteepest then candidateNetworkEdgeList'
                                     else if doRandomOrder then permuteList (head rSeedList) candidateNetworkEdgeList'
                                     else candidateNetworkEdgeList'

          -- newGraphList = concat (fmap (insertNetEdgeBothDirections inGS inData inPhyloGraph) candidateNetworkEdgeList `using`  PU.myParListChunkRDS)
          newGraphList = if not doSteepest then filter (/= emptyPhylogeneticGraph) (fmap (insertNetEdge inGS inData leafGraph inPhyloGraph preDeleteCost) candidateNetworkEdgeList `using`  PU.myParListChunkRDS)
                         else insertNetEdgeRecursive inGS inData leafGraph rSeedList maxNetEdges doSteepest doRandomOrder inPhyloGraph preDeleteCost inSimAnnealParams candidateNetworkEdgeList

          --minCostGraphList = GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] newGraphList
          minCost = if null candidateNetworkEdgeList || null newGraphList then infinity
                    else minimum $ fmap snd6 newGraphList
      in
      trace ("\tExamining at most " ++ (show $ length candidateNetworkEdgeList) ++ " candidate edge pairs") (

      -- no network edges to insert
      -- trace ("IENE: " ++ (show minCost)) (
      if null candidateNetworkEdgeList then 
         -- trace ("IENE num cand edges:" ++ (show $ length candidateNetworkEdgeList)) 
         ([inPhyloGraph], currentCost)
      -- filter later
      -- single if steepest so no neeed to unique
      else if doSteepest then 
         -- trace ("IENE: Steepest " ++ (show minCost))
         (newGraphList, minCost)

      else 
         --  trace ("IENE: All " ++ (show minCost))
         (take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] $ newGraphList, minCost)
      ) -- )

-- | insertNetEdgeRecursive recursively inserts edges and returns new graph only if better
insertNetEdgeRecursive :: GlobalSettings 
                       -> ProcessedData 
                       -> DecoratedGraph
                       -> [Int] 
                       -> Int
                       -> Bool 
                       -> Bool 
                       ->  PhylogeneticGraph 
                       -> Maybe VertexCost 
                       -> Maybe SAParams 
                       -> [(LG.LEdge EdgeInfo, LG.LEdge EdgeInfo)] 
                       -> [PhylogeneticGraph]
insertNetEdgeRecursive inGS inData leafGraph rSeedList maxNetEdges doSteepest doRandomOrder inPhyloGraph preDeleteCost inSimAnnealParams edgePairList =
   -- trace ("Edges pairs to go : " ++ (show $ length edgePairList)) (
   if null edgePairList then [inPhyloGraph]
   else
      let 
          firstEdgePair = head edgePairList

          --check for max net edges
          (_, _, _, netNodes) = LG.splitVertexList (thd6 inPhyloGraph)

          -- needf to check disapaly/charcter trees not conical graph
          newGraph = insertNetEdge inGS inData leafGraph inPhyloGraph preDeleteCost firstEdgePair

          
      in
      traceNoLF ("*")  (      -- trace ("INER: " ++ (show $ snd6 newGraph) ++ " " ++ (show preDeleteCost)) (
      if length netNodes >= maxNetEdges then
         trace ("Maximum number of network edges reached: " ++ (show $ length netNodes))
         [inPhyloGraph]

      -- malformed graph
      else if newGraph == emptyPhylogeneticGraph then 
         -- trace ("INER: Empty more to go : " ++ (show $ length $ tail edgePairList)) 
         insertNetEdgeRecursive inGS inData leafGraph rSeedList maxNetEdges doSteepest doRandomOrder inPhyloGraph preDeleteCost inSimAnnealParams (tail edgePairList)

      else if (inSimAnnealParams == Nothing) then
         -- better cost
         if snd6 newGraph < snd6 inPhyloGraph then 
            -- cyclic check in insert edge function
            -- trace ("INER: Better -> " ++ (show $ snd6 newGraph)) 
            [newGraph]

            {-
            let isCyclic = V.filter (== True) $ fmap LG.cyclic $ fmap head $ fth6 newGraph
                dupList = V.filter (== True) $ fmap LG.hasDuplicateEdge $ fmap head $ fth6 newGraph
            in
            if (null isCyclic && null dupList) then 
               -- trace ("INER: Better -> " ++ (show $ snd6 newGraph)) 
               -- trace  ("INER:" ++ (LG.prettyIndices $ thd6 newGraph)) 
               [newGraph]
            else 
               trace ("Cycle " ++ (show isCyclic) ++ " Duplicate Edges " ++ (show dupList)) 
               insertNetEdgeRecursive inGS inData rSeedList maxNetEdges doSteepest doRandomOrder inPhyloGraph preDeleteCost inSimAnnealParams (tail edgePairList)
            -}

         -- not better
         else -- trace ("INER: Really Not Better")  
            insertNetEdgeRecursive inGS inData leafGraph rSeedList maxNetEdges doSteepest doRandomOrder inPhyloGraph preDeleteCost inSimAnnealParams (tail edgePairList)

      -- sim annealing/drift
      else
         let numDone = if (method $ fromJust inSimAnnealParams) == SimAnneal then currentStep $ fromJust inSimAnnealParams
                       else driftChanges $ fromJust inSimAnnealParams
             numMax  = if (method $ fromJust inSimAnnealParams) == SimAnneal then numberSteps $ fromJust inSimAnnealParams
                       else driftMaxChanges $ fromJust inSimAnnealParams
             (acceptGraph, nextSAParams) = U.simAnnealAccept inSimAnnealParams (snd6 inPhyloGraph) (snd6 newGraph)
         in
         if (numDone < numMax) then
            if acceptGraph then [newGraph]
            else insertNetEdgeRecursive inGS inData leafGraph rSeedList maxNetEdges doSteepest doRandomOrder inPhyloGraph preDeleteCost nextSAParams (tail edgePairList)

         -- hit end of SA/Drift
         else [inPhyloGraph]
      )
      -- )

-- | insertNetEdge inserts an edge between two other edges, creating 2 new nodes and rediagnoses graph
-- contacts deletes 2 orginal edges and adds 2 nodes and 5 new edges
-- does not check any edge reasonable-ness properties
-- new edge directed from first to second edge
-- naive for now
-- predeletecost of edge move
insertNetEdge :: GlobalSettings -> ProcessedData -> DecoratedGraph -> PhylogeneticGraph -> Maybe VertexCost -> (LG.LEdge b, LG.LEdge b) -> PhylogeneticGraph
insertNetEdge inGS inData leafGraph inPhyloGraph preDeleteCost edgePair@((u,v, _), (u',v', _)) =
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

           -- conversion as if input--see if affects length
               -- removed check after checks moved to permissible edges
               -- can add back if there are malformed graphs being generated
           -- newSimple' = GO.convertGeneralGraphToPhylogeneticGraph newSimple
               -- permissibale not catching timeconsistency issues with edges
           newSimple' = GO.makeGraphTimeConsistent newSimple


           -- full two-pass optimization
           newPhyloGraph = T.multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex newSimple'

           -- calculates heursitic graph delta
           -- (heuristicDelta, _, _, _, _)  = heuristicAddDelta inGS inPhyloGraph edgePair (fst newNodeOne) (fst newNodeTwo)
           -- heuristicDelta = heuristicAddDelta' inGS inPhyloGraph edgePair


           -- edgeAddDelta = deltaPenaltyAdjustment inGS inPhyloGraph "add"

       in

       if (graphType inGS) == HardWired then newPhyloGraph

       else 
         -- need heuristics in here
         -- if (heuristicDelta + edgeAddDelta) < 0 then newPhyloGraph
         
         {-
         if True then 
            trace ("INE: " ++ (show (heuristicDelta, edgeAddDelta, snd6 inPhyloGraph)) ++ " -> " ++ (show (heuristicDelta + edgeAddDelta + (snd6 inPhyloGraph), snd6 inPhyloGraph)))
            newPhyloGraph
         -}
         -- if (heuristicDelta + edgeAddDelta) < 0 then newPhyloGraph
         if True then newPhyloGraph
         else emptyPhylogeneticGraph

       {-
       -- preDelete cost changes criterion for edge move
       else if isNothing preDeleteCost then
         -- trace ("INE: " ++ (show (heuristicDelta, edgeAddDelta, (snd6 inPhyloGraph), (snd6 newPhyloGraph)))) ( 
         if True then -- heuristicDelta + edgeAddDelta < 0 then 
            if (snd6 newPhyloGraph) < (snd6 inPhyloGraph) then 
               -- trace ("INE: Better")
               newPhyloGraph
            else 
               -- trace ("INE: Worse")
               emptyPhylogeneticGraph
         else 
            -- trace ("INE: Worse")
            emptyPhylogeneticGraph
         -- )

       else
         -- trace ("INE: " ++ (show $ LG.isEmpty newSimple) ++ " " ++ (show $ LG.isEmpty $ thd6 newPhyloGraph) ++ " " ++ (show $ isNothing preDeleteCost) ++ " " ++ (show $ fromJust preDeleteCost) ++ " versus " ++ (show $ snd6 newPhyloGraph) ++ " incost " ++ (show $ snd6 inPhyloGraph)) ( --  ++ "\nNewGraph:" ++ (LG.prettyIndices $ thd6 inPhyloGraph)) (
         -- no net add cost because the number of net nodes is unchanged in add/delete when preDelete cost /= Noting
         --if heuristicDelta + (snd6 inPhyloGraph) <= fromJust preDeleteCost then newPhyloGraph
         if (snd6 newPhyloGraph) < fromJust preDeleteCost then newPhyloGraph
         else emptyPhylogeneticGraph
         -- )
       -- )
       -}


-- | (curBestGraphList, annealBestCost) is a wrapper for moveAllNetEdges' allowing for multiple simulated annealing rounds
deleteAllNetEdges :: GlobalSettings -> ProcessedData -> Int -> Int -> Int -> Int -> Bool -> Bool -> Bool -> ([PhylogeneticGraph], VertexCost) -> (Maybe SAParams, [PhylogeneticGraph]) -> ([PhylogeneticGraph], Int)
deleteAllNetEdges inGS inData rSeed maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
   let leafGraph = LG.extractLeafGraph $ (thd6 . head) inPhyloGraphList
   in
   if inSimAnnealParams == Nothing then
      deleteAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (randomIntList rSeed) inSimAnnealParams inPhyloGraphList
   else
      let -- create list of params with unique list of random values for rounds of annealing
          annealingRounds = rounds $ fromJust inSimAnnealParams
          annealParamGraphList = U.generateUniqueRandList annealingRounds inSimAnnealParams
          replicateRandIntList = fmap randomIntList (take annealingRounds (randomIntList rSeed))

          (annealRoundsList, counterList) = unzip (zipWith3 (deleteAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost)) replicateRandIntList annealParamGraphList (replicate annealingRounds inPhyloGraphList) `using` PU.myParListChunkRDS)
      in
      (take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] (concat annealRoundsList) , sum counterList)

-- | deleteAllNetEdges deletes network edges one each each round until no better or additional
-- graphs are found
-- call with ([], infinity) [single input graph]
deleteAllNetEdges' :: GlobalSettings -> ProcessedData -> DecoratedGraph -> Int ->Int ->  Int -> Bool -> Bool -> Bool -> ([PhylogeneticGraph], VertexCost) -> [Int] -> Maybe SAParams -> [PhylogeneticGraph]-> ([PhylogeneticGraph], Int)
deleteAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) randIntList inSimAnnealParams inPhyloGraphList =
   -- trace ("In deleteAllNetEdges " ++ (show $ length inPhyloGraphList)) ( 
   if null inPhyloGraphList then (take numToKeep curBestGraphList, counter)
   else
      let currentCost = min curBestGraphCost (snd6 $ head inPhyloGraphList)

          (newGraphList', _) = deleteEachNetEdge inGS inData leafGraph (head randIntList) numToKeep doSteepest doRandomOrder False inSimAnnealParams (head inPhyloGraphList)

          newGraphList = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] newGraphList'
          newGraphCost = if (not . null) newGraphList then snd6 $ head newGraphList
                         else infinity

      in
      -- trace ("DANE: " ++ (show (newGraphCost, length newGraphList))) (
      -- if graph is a tree no edges to delete
      if LG.isTree (fst6 $ head inPhyloGraphList) then
         let (a,b,c,d) = LG.splitVertexList (fst6 $ head inPhyloGraphList)
         in
         trace ("\tGraph in delete network edges is tree--skipping :" ++ (show $ (snd6 $ head inPhyloGraphList, length a, length b, length c, length d)))
         deleteAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder ((head inPhyloGraphList) : curBestGraphList, currentCost) (tail randIntList) inSimAnnealParams (tail inPhyloGraphList)

      else if null newGraphList then (take numToKeep curBestGraphList, counter + 1)

         -- regular delte wihtout simulated annealing
         -- worse graphs found--go on
      else if inSimAnnealParams == Nothing then
         if  newGraphCost > currentCost then deleteAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder ((head inPhyloGraphList) : curBestGraphList, currentCost) (tail randIntList) inSimAnnealParams (tail inPhyloGraphList)

         -- "steepest style descent" abandons existing list if better cost found
         else if newGraphCost < currentCost then
            trace ("\t-> " ++ (show newGraphCost))
            deleteAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (newGraphList, newGraphCost) (tail randIntList) inSimAnnealParams (newGraphList ++ (tail inPhyloGraphList))

         -- equal cost
         -- not sure if should add new graphs to queue to do edge deletion again
         else
            -- new grapjh list contains the input graph if equal and filterd unique already in deleteEachNetEdge
            let newCurSameBestList =  take numToKeep $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (curBestGraphList ++ newGraphList)
            in
            deleteAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (newCurSameBestList, currentCost) (tail randIntList) inSimAnnealParams (tail inPhyloGraphList)

      -- simulated annealing
      else
         -- if steepest -- simAnneal/Drift stuff done during net delete function so jusyt take results and move on
            -- create new
            if doSteepest then
               let annealBestCost = min curBestGraphCost newGraphCost
                   newSAParams = Just $ (fromJust (U.incrementSimAnnealParams newSAParams)) {currentStep = 0, driftChanges = 0}
               in
               deleteAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder ((head newGraphList') :  curBestGraphList, annealBestCost) (tail randIntList) newSAParams (newGraphList' ++ (tail inPhyloGraphList))

            -- simmAnneal on "all list"
            else
            -- sim anneal choice
            let numDone = if (method $ fromJust inSimAnnealParams) == SimAnneal then currentStep $ fromJust inSimAnnealParams
                         else driftChanges $ fromJust inSimAnnealParams
                numMax  = if (method $ fromJust inSimAnnealParams) == SimAnneal then numberSteps $ fromJust inSimAnnealParams
                         else driftMaxChanges $ fromJust inSimAnnealParams
                annealBestCost = min curBestGraphCost (snd6 $ head newGraphList')
                (acceptFirstGraph, nextSAParams) = U.simAnnealAccept inSimAnnealParams annealBestCost (snd6 $ head newGraphList')
            in
            -- trace ("ACG" ++ (show acceptFirstGraph) ++ " " ++ (show $ snd6 $ head uniqueGraphList)) (
            if (numDone < numMax) then
               let nextNewGraphList = if (not . null) newGraphList' then tail newGraphList'
                                    else []
               in

               if acceptFirstGraph then
                  deleteAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder ((head newGraphList') :  curBestGraphList, annealBestCost) (tail randIntList) nextSAParams (nextNewGraphList ++ (tail inPhyloGraphList))
               else
                  deleteAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep (counter + 1) returnMutated doSteepest doRandomOrder (curBestGraphList, annealBestCost)  (tail randIntList) nextSAParams (nextNewGraphList ++ (tail inPhyloGraphList))

            -- if want non-optimized list for GA or whatever
            else if returnMutated then (take numToKeep curBestGraphList, counter)

            -- optimize with net insert to add back edges to optimiality
            else
               let (bestList', counter') =  insertAllNetEdges' inGS inData leafGraph maxNetEdges numToKeep (counter + 1) False doSteepest doRandomOrder ([], annealBestCost) (tail randIntList) Nothing (take numToKeep curBestGraphList)
                   bestList = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] bestList'
               in
               --trace ("BM: " ++ (show $ snd6 $ head  bestMoveList))
               (take numToKeep bestList, counter')
            -- ) )


-- | deleteOneNetAddAll new version deletes net edges in turn and readds-based on original cost
-- but his cost in graph (really not correct) but allows logic of insert edge to function better
-- seems like after delete--or insert--the graphs are i nimporpper condition and returnin infinite cost due to that
-- seems likely true fir deleteOneNetAddAll' as well
-- this has no SA in it
deleteOneNetAddAll :: GlobalSettings 
                   -> ProcessedData 
                   -> DecoratedGraph 
                   -> Int 
                   -> Int 
                   -> Bool 
                   -> Bool 
                   -> PhylogeneticGraph 
                   -> [LG.Edge] 
                   -> Int 
                   -> Maybe SAParams 
                   -> [PhylogeneticGraph]
deleteOneNetAddAll inGS inData leafGraph maxNetEdges numToKeep doSteepest doRandomOrder inPhyloGraph edgeToDeleteList rSeed inSimAnnealParams =
   if null edgeToDeleteList then 
      -- trace ("\tGraph has no edges to move---skipping") 
      [inPhyloGraph]
   else if LG.isEmpty $ thd6 inPhyloGraph then error "Empty graph in deleteOneNetAddAll"
   else
      -- trace ("DONAA-New: " ++ (show $ snd6 inPhyloGraph) ++ " Steepest:" ++ (show doSteepest)) (
      trace ("Moving " ++ (show $ length edgeToDeleteList) ++ " network edges, curent best cost: " ++ (show $ snd6 inPhyloGraph)) (
      -- start with initial graph cost
      let inGraphCost = snd6 inPhyloGraph

          -- get deleted simple graphs and bool for chenged
          delGraphBoolPairList = fmap (flip deleteNetworkEdge (fst6 inPhyloGraph)) edgeToDeleteList
          (simpleGraphsToInsert, _) = unzip $ filter ((== True ) . snd) delGraphBoolPairList

          -- check for cycles -- already done
          -- simpleGraphsToInsert' = filter ((== False) . LG.cyclic) simpleGraphsToInsert

          -- optimize deleted graph and update cost with input cost
          graphsToInsert = fmap (T.multiTraverseFullyLabelSoftWired inGS inData False False leafGraph Nothing) simpleGraphsToInsert -- `using` PU.myParListChunkRDS

          -- keep same cost and just keep better--check if better than original later
          graphsToInsert' = fmap (flip T.updatePhylogeneticGraphCost inGraphCost) graphsToInsert


          insertedGraphPairList = if not doSteepest then fmap (insertEachNetEdge inGS inData leafGraph rSeed maxNetEdges numToKeep doSteepest doRandomOrder Nothing inSimAnnealParams) graphsToInsert' -- `using` PU.myParListChunkRDS
                                  else 
                                       let -- potentially randomize order of list
                                           graphsToInsert'' = if not doRandomOrder then graphsToInsert'
                                                              else permuteList rSeed graphsToInsert'
                                       in
                                       [insertEachNetEdgeRecursive inGS inData leafGraph maxNetEdges numToKeep doSteepest doRandomOrder (head $ randomIntList rSeed)inSimAnnealParams graphsToInsert'']


          newMinimumCost = minimum $ fmap snd insertedGraphPairList

          newBestGraphs = filter ((== newMinimumCost) . snd6) $ concat $ fmap fst insertedGraphPairList

      in
      trace ("DONAA-New: " ++ (show (inGraphCost, fmap snd6 graphsToInsert, fmap snd6 graphsToInsert', newMinimumCost))) (
      if null simpleGraphsToInsert then [emptyPhylogeneticGraph]

      else if newMinimumCost < inGraphCost then
         trace ("-> ")
         newBestGraphs

      else [inPhyloGraph]

     ) )


-- | insertEachNetEdgeRecursive is a wrapper arounf insertEachNet each for edge move heuristic and reutnrs better graph when found immediately
insertEachNetEdgeRecursive  :: GlobalSettings 
                            -> ProcessedData 
                            -> DecoratedGraph 
                            -> Int 
                            -> Int 
                            -> Bool 
                            -> Bool 
                            -> Int 
                            -> Maybe SAParams 
                            -> [PhylogeneticGraph] 
                            -> ([PhylogeneticGraph], VertexCost)
insertEachNetEdgeRecursive inGS inData leafGraph maxNetEdges numToKeep doSteepest doRandomOrder rSeed inSimAnnealParams inPhyloGraphList =
   if null inPhyloGraphList then ([],infinity)
   else 
         let firstGraph = head inPhyloGraphList
             (newGraphList, newCost) = insertEachNetEdge inGS inData leafGraph rSeed maxNetEdges numToKeep doSteepest doRandomOrder Nothing inSimAnnealParams firstGraph

             minCost = if null newGraphList then infinity
                       else newCost
         in
         -- return immediately 
         if minCost < (snd6 firstGraph) then 
            (take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] $ newGraphList, minCost)
         else 
            insertEachNetEdgeRecursive inGS inData leafGraph maxNetEdges numToKeep doSteepest doRandomOrder rSeed inSimAnnealParams  (tail inPhyloGraphList) 



-- | deleteOneNetAddAll deletes the specified edge from a graph--creating a fully optimized new one--then re-adds
-- and keeps best based on delta, reoptimizes those and compares to the oringal cost
-- if better or ssame keeps as per usual
deleteOneNetAddAll' :: GlobalSettings 
                    -> ProcessedData 
                    -> DecoratedGraph 
                    -> Int 
                    -> Int 
                    -> Bool 
                    -> Bool 
                    -> PhylogeneticGraph 
                    -> LG.Edge 
                    -> Int 
                    -> Maybe SAParams 
                    -> [PhylogeneticGraph]
deleteOneNetAddAll' inGS inData leafGraph maxNetEdges numToKeep doSteepest doRandomOrder inPhyloGraph edgeToDelete rSeed inSimAnnealParams  =
   if LG.isEmpty $ thd6 inPhyloGraph then error "Empty graph in deleteOneNetAddAll"
   else
      trace ("DONAA: moving " ++ (show edgeToDelete) ++ " Steepest:" ++ (show doSteepest)) (
      -- True to force reoptimization of delete
      let inGraphCost = snd6 inPhyloGraph
          
          -- deletedEdgeGraph = deleteNetEdge inGS inData inPhyloGraph True edgeToDelete
          
          -- generate only simple graph--doni't need costs or anything
          -- pass minimal phylogenetic graph
          (deletedGraphSimple, wasModified)  = deleteNetworkEdge edgeToDelete $ fst6 inPhyloGraph
          delPhyloGraph =  if (graphType inGS == SoftWired) then T.multiTraverseFullyLabelSoftWired inGS inData False False leafGraph Nothing deletedGraphSimple
                           else if (graphType inGS == HardWired) then T.multiTraverseFullyLabelHardWired inGS inData leafGraph Nothing deletedGraphSimple
                           else error "Unsupported graph type in deleteNetEdge.  Must be soft or hard wired"deletedGraphSimple

          (insertedGraphList, _) = insertEachNetEdge inGS inData leafGraph rSeed maxNetEdges numToKeep doSteepest doRandomOrder (Just inGraphCost) inSimAnnealParams delPhyloGraph
      in
      -- this checks to see if input graph edge was deleted--it might not be due to phylogenetic graph constraints
      -- if graph not modified just returns the original--does not insert new edges
      if not wasModified then trace ("DONAA: skipping") [emptyPhylogeneticGraph]

      else
         trace ("DONAA returning: " ++ (show $ fmap snd6 insertedGraphList))
         filter ((/= infinity) . snd6) $ take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] insertedGraphList
      )


-- | getPermissibleEdgePairs takes a DecoratedGraph and returns the list of all pairs
-- of edges that can be joined by a network edge and meet all necessary conditions

-- add in other conditions
--   reproducable--ie not tree noide with two net node children--other stuff
getPermissibleEdgePairs :: DecoratedGraph -> [(LG.LEdge EdgeInfo, LG.LEdge EdgeInfo)]
getPermissibleEdgePairs inGraph =
   if LG.isEmpty inGraph then error "Empty input graph in isEdgePairPermissible"
   else
       let edgeList = LG.labEdges inGraph

           -- edges to potentially conenct
           edgePairs = cartProd edgeList edgeList
           
           -- get coeval node pairs in existing grap
           coevalNodeConstraintList = LG.coevalNodePairs inGraph
           coevalNodeConstraintList' = fmap (GO.addBeforeAfterToPair inGraph) coevalNodeConstraintList `using`  PU.myParListChunkRDS
           
           edgeTestList = fmap (isEdgePairPermissible inGraph coevalNodeConstraintList') edgePairs `using`  PU.myParListChunkRDS
           pairList = fmap fst $ filter ((== True) . snd) $ zip edgePairs edgeTestList
       in
       -- trace ("Edge Pair list :" ++ (show $ fmap f pairList) ++ "\n"
       --  ++ "GPEP\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inGraph))
       pairList
       -- where f (a, b) = (LG.toEdge a, LG.toEdge b)

getPermissibleEdgePairs' :: DecoratedGraph -> [(LG.LEdge EdgeInfo, LG.LEdge EdgeInfo)]
getPermissibleEdgePairs' inGraph =
   if LG.isEmpty inGraph then error "Empty input graph in isEdgePairPermissible"
   else
       let edgeList = LG.labEdges inGraph
           edgePairs = cartProd edgeList edgeList
           contraintList = LG.getGraphCoevalConstraints inGraph
           edgeTestList = fmap (isEdgePairPermissible' inGraph contraintList) edgePairs `using`  PU.myParListChunkRDS
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
--    3) neither neither edge is an ancestor or descndent edge of the other (tested via bv of nodes)
-- the result should apply to a new edge in either direction
-- new edge to be creted is edge1 -> ege2
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


-- | isEdgePairPermissible takes a graph and two edges, coeval contraints, and tests whether a
-- pair of edges can be linked by a new edge and satify three consitions:
--    1) neither edge is a network edge
--    2) one edge cannot be "before" while the other is "after" in any of the constraint pairs
--    3) neither neither edge is an ancestor or descndent edge of the other (tested via bv of nodes)
-- the result should apply to a new edge in either direction
-- new edge to be creted is edge1 -> ege2
isEdgePairPermissible' :: DecoratedGraph -> [([LG.LEdge EdgeInfo],[LG.LEdge EdgeInfo])] -> (LG.LEdge EdgeInfo, LG.LEdge EdgeInfo) -> Bool
isEdgePairPermissible' inGraph constraintList (edge1@(u,v,_), edge2@(u',v',_)) =
   if LG.isEmpty inGraph then error "Empty input graph in isEdgePairPermissible"
   else
       if u == u' then False
       else if v == v' then False
       -- equality implied in above two
       -- else if LG.toEdge edge1 == LG.toEdge edge2 then False
       else if (LG.isNetworkNode inGraph u) || (LG.isNetworkNode inGraph u') then False
       else if (LG.isNetworkLabEdge inGraph edge1) || (LG.isNetworkLabEdge inGraph edge2) then False
       else if not (LG.meetsAllCoevalConstraintsEdges constraintList edge1 edge2) then False
       else if (isAncDescEdge inGraph edge1 edge2) then False
       -- get children of u' to make sure no net children
       else if (not . null) $ filter (== True) $ fmap (LG.isNetworkNode inGraph) $ LG.descendants inGraph u' then False
       else True



-- | isAncDescEdge takes a graph and two edges and examines whethe either edge is the ancestor or descendent of the other
-- this is done via examination of teh bitvector fields of the node
isAncDescEdge :: DecoratedGraph ->  LG.LEdge EdgeInfo -> LG.LEdge EdgeInfo -> Bool
isAncDescEdge inGraph (a,_,_) (b, _, _) =
   if LG.isEmpty inGraph then error "Empty input graph in isAncDescEdge"
   else
      let aBV = bvLabel $ fromJust $ LG.lab inGraph a
          bBV = bvLabel $ fromJust $ LG.lab inGraph b
      in
      if aBV .&. bBV == aBV then True
      else if aBV .&. bBV == bBV then True
      else False

{-
-- | insertNetEdgeBothDirections calls insertNetEdge for both u -> v and v -> u new edge orientations
insertNetEdgeBothDirections :: GlobalSettings -> ProcessedData -> PhylogeneticGraph ->  Maybe VertexCost -> (LG.LEdge b, LG.LEdge b) -> [PhylogeneticGraph]
insertNetEdgeBothDirections inGS inData inPhyloGraph preDeleteCost (u,v) = fmap (insertNetEdge inGS inData inPhyloGraph preDeleteCost) [(u,v), (v,u)]
-}

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
--    graphDelta = egdeAddDelta (separately calcualated) + sum [blockDelta]
--    Compare to real delta to check behavior
-- original subtrees u -> (a,v) and u' -> (v',b)
heuristicAddDelta' :: GlobalSettings -> PhylogeneticGraph -> (LG.LEdge b, LG.LEdge b) -> VertexCost
heuristicAddDelta' inGS inPhyloGraph ((u,v, _), (u',v', _)) = 
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
getCharacterDelta (u,v,u',v',a,b) inCharTree charInfo =
   let doIA = False
       filterGaps = True
       uData = V.head $ V.head $ vertData $ fromJust $ LG.lab inCharTree u
       vData = V.head $ V.head $ vertData $ fromJust $ LG.lab inCharTree v
       vFinalData = V.head $ V.head $ PRE.setPreliminaryToFinalStates  $ vertData $ fromJust $ LG.lab inCharTree v
       u'Data = V.head $ V.head $ vertData $ fromJust $ LG.lab inCharTree u'
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
heuristicAddDelta :: GlobalSettings -> PhylogeneticGraph -> (LG.LEdge b, LG.LEdge b) -> LG.Node -> LG.Node -> (VertexCost, LG.LNode VertexInfo, LG.LNode VertexInfo, LG.LNode VertexInfo, LG.LNode VertexInfo)
heuristicAddDelta inGS inPhyloGraph ((u,v, _), (u',v', _)) n1 n2 =
  if LG.isEmpty (fst6 inPhyloGraph) then error "Empty graph in heuristicAddDelta"
  else if graphType inGS == HardWired then
      let uvVertData = M.makeEdgeData  False (thd6 inPhyloGraph) (six6 inPhyloGraph) (u, v, dummyEdge)
          uvPrimeData =  M.makeEdgeData  False (thd6 inPhyloGraph) (six6 inPhyloGraph) (u', v', dummyEdge)
          hardDelta = V.sum $ fmap V.sum $ fmap (fmap snd) $ POS.createVertexDataOverBlocks uvVertData uvPrimeData (six6 inPhyloGraph) []
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
          n2Lab          = POSW.getOutDegree1VertexSoftWired n2 vPrimeLab (thd6 inPhyloGraph) [n2]
          uPrimeLabAfter = POSW.getOutDegree2VertexSoftWired inGS (six6 inPhyloGraph) u' (n2, n2Lab) uPrimeOtherChild (thd6 inPhyloGraph)
          n1Lab          = POSW.getOutDegree2VertexSoftWired inGS (six6 inPhyloGraph) n1 (v, vLab) (n2, n2Lab) (thd6 inPhyloGraph)
          uLabAfter      = POSW.getOutDegree2VertexSoftWired inGS (six6 inPhyloGraph) u uOtherChild (n1, n1Lab) (thd6 inPhyloGraph)

          -- cost of resolutions
          (_, uCostBefore) = POSW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uLab)
          (_, uPrimeCostBefore) = POSW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uPrimeLab)
          (_, uCostAfter) = POSW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uLabAfter)
          (_, uPrimeCostAfter) = POSW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uPrimeLabAfter)

          addNetDelta = uCostAfter - uCostBefore +  uPrimeCostAfter - uPrimeCostBefore


      in
      trace ("HAD: " ++ (show (uCostAfter, uCostBefore, uPrimeCostAfter, uPrimeCostBefore) ++ " -> " ++ (show $ uCostAfter - uCostBefore +  uPrimeCostAfter - uPrimeCostBefore))) (
      if null (filter ((/= v') . fst) $ LG.labDescendants (thd6 inPhyloGraph) (u', uPrimeLab)) || null (filter ((/= v) . fst) $ LG.labDescendants (thd6 inPhyloGraph) (u, uLab)) then (infinity, dummyNode, dummyNode, dummyNode, dummyNode)
      -- this should not happen--should try to crete new edges from children of net edges
      else if (length $ LG.descendants (thd6 inPhyloGraph) u) < 2 ||  (length $ LG.descendants (thd6 inPhyloGraph) u') < 2 then error ("Outdegree 1 nodes in heuristicAddDelta")
      else
         (addNetDelta, (u, uLabAfter), (u', uPrimeLabAfter), (n1, n1Lab), (n2, n2Lab))
      )



-- | deltaPenaltyAdjustment takes number of leaves and Phylogenetic graph and returns a heuristic graph penalty for adding a single network edge
-- if Wheeler2015Network, this is based on a all changes affecting a single block (most permissive)  and Wheeler 2015 calculation of penalty
-- if PMDLGraph -- KMDL not yet implemented
-- if NoNetworkPenalty then 0
-- modification "add" or subtrct to calculate delta
-- always delta is positive--whether neg or pos is deltermined when used
deltaPenaltyAdjustment :: GlobalSettings -> PhylogeneticGraph -> String -> VertexCost
deltaPenaltyAdjustment inGS inGraph modification =
   -- trace ("DPA: entering: " ++ (show $ graphFactor inGS)) (
   let numLeaves = numDataLeaves inGS
       edgeCostModel = graphFactor inGS
       (_, _, _, networkNodeList) = LG.splitVertexList (fst6 inGraph)
   in
   if edgeCostModel == NoNetworkPenalty then 
      -- trace ("DPA: No penalty") 
      0.0

   else if length networkNodeList == 0 then 
      -- trace ("DPA: No cost") 
      0.0      
   
   else if edgeCostModel == Wheeler2015Network then
      -- trace  ("DPW: In Wheeler2015Network") (
      let graphCost = snd6 inGraph -- this includes any existing penalties--would be better not to include
          numBlocks = V.length $ fth6 inGraph
      in
      -- trace ("DPA Value: " ++ (show $ graphCost / (fromIntegral $ numBlocks * 2 * ((2 * numLeaves) - 2)))) 
      graphCost / (fromIntegral $ numBlocks * 2 * ((2 * numLeaves) - 2))
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
deleteEachNetEdge :: GlobalSettings -> ProcessedData -> DecoratedGraph -> Int -> Int -> Bool ->  Bool -> Bool-> Maybe SAParams -> PhylogeneticGraph -> ([PhylogeneticGraph], VertexCost)
deleteEachNetEdge inGS inData leafGraph rSeed numToKeep doSteepest doRandomOrder force inSimAnnealParams inPhyloGraph =
   -- trace ("DENE start") (
   if LG.isEmpty $ thd6 inPhyloGraph then ([], infinity) -- error "Empty input phylogenetic graph in deleteAllNetEdges"
   else
      let currentCost = snd6 inPhyloGraph

          -- potentially randomize order of list
          networkEdgeList' = LG.netEdges $ thd6 inPhyloGraph
          networkEdgeList = if not doRandomOrder then networkEdgeList'
                            else permuteList rSeed networkEdgeList'


          newGraphList = if not doSteepest then fmap (deleteNetEdge inGS inData leafGraph inPhyloGraph force) networkEdgeList `using`  PU.myParListChunkRDS
                         else deleteNetEdgeRecursive inGS inData leafGraph inPhyloGraph force inSimAnnealParams networkEdgeList

          uniqueCostGraphList = filter ((/= infinity) . snd6) $ take numToKeep $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] newGraphList
          minCost = if null uniqueCostGraphList then infinity
                    else minimum $ fmap snd6 uniqueCostGraphList
      in
      -- no network edges to delete
      if null networkEdgeList then trace ("\tNo network edges to delete") ([inPhyloGraph], currentCost)

      else
         -- single if steepest so no neeed to unique
         if doSteepest then (newGraphList, minCost)

         else
            -- trace ("DENE end: " ++ (show (currentCost, minCost, length newGraphList, length uniqueCostGraphList))) (
            if minCost < currentCost then 
               -- trace ("DENE--Delete net edge return:" ++ (show (minCost,length uniqueCostGraphList))) (
               let nextGraphPairList = fmap (deleteEachNetEdge inGS inData leafGraph rSeed numToKeep doSteepest doRandomOrder force inSimAnnealParams) (filter ((== minCost) .snd6) uniqueCostGraphList) `using`  PU.myParListChunkRDS
                   newMinCost = minimum $ fmap snd nextGraphPairList
                   newGraphListBetter = filter ((== newMinCost) . snd6) $ concatMap fst nextGraphPairList
               in
               (take numToKeep $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] $ newGraphListBetter, newMinCost)

               -- )

            else 
               -- filter later
               -- trace ("DENE returning same")
               (take numToKeep $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] $ uniqueCostGraphList, currentCost)

-- | deleteNetworkEdge deletes a network edges from a simple graph
-- retuns newGraph if can be modified or input graph with Boolean to tell if modified
-- and contracts, reindexes/names internaledges/veritices around deletion 
-- can't raise to general graph level due to vertex info 
-- in edges (b,a) (c,a) (a,d), deleting (a,b) deletes node a, inserts edge (b,d)
-- contacts node c since  now in1out1 vertex
-- checks for chained network edges--can be created by progressive deletion
-- checks for cycles now
-- shouldn't need for check for creting a node with children that are both network nodes
-- sine that woudl require that condition coming in and shodl be there--ie checked earlier in addition and input
deleteNetworkEdge :: LG.Edge -> SimpleGraph -> (SimpleGraph, Bool)
deleteNetworkEdge inEdge@(p1, nodeToDelete) inGraph =
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
          -- newGraph'' = GO.convertGeneralGraphToPhylogeneticGraph newGraph
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
      
      else 
         {-trace ("DNE: Edge to delete " ++ (show inEdge) ++ " cnd " ++ (show childrenNodeToDelete) ++ " pnd " ++ (show parentsNodeToDelete) ++ " pntk " ++ (show parentNodeToKeep) 
            ++ " ne " ++ (show newEdge) ++ "\nInGraph: " ++ (LG.prettyIndices inGraph) ++ "\nNewGraph: " ++ (LG.prettyIndices newGraph) ++ "\nNewNewGraph: " 
            ++ (LG.prettyIndices newGraph')) -}
         (newGraph'', True)


-- | deleteNetEdgeRecursive like deleteEdge, deletes an edge (checking if network) and rediagnoses graph
-- contacts in=out=1 edges and removes node, reindexing nodes and edges
-- except returns on first better (as opposed to do all deletes first) 
-- or sim annleal/drift
deleteNetEdgeRecursive :: GlobalSettings -> ProcessedData -> DecoratedGraph -> PhylogeneticGraph -> Bool -> Maybe SAParams -> [LG.Edge] -> [PhylogeneticGraph]
deleteNetEdgeRecursive inGS inData leafGraph inPhyloGraph force inSimAnnealParams edgeToDeleteList =
   if null edgeToDeleteList then []
   else
       let edgeToDelete = head edgeToDeleteList

           -- calls general funtion to remove network graph edge
           (delSimple, wasModified) = deleteNetworkEdge edgeToDelete $ fst6 inPhyloGraph

           -- delSimple = GO.contractIn1Out1EdgesRename $ LG.delEdge edgeToDelete $ fst6 inPhyloGraph
   
           -- prune other edges if now unused
           pruneEdges = False

           -- don't warn that edges are being pruned
           warnPruneEdges = False

           -- graph optimization from root
           startVertex = Nothing

           (heuristicDelta, _, _) = heuristicDeleteDelta inGS inPhyloGraph edgeToDelete


           edgeAddDelta = deltaPenaltyAdjustment inGS inPhyloGraph "delete"


           -- full two-pass optimization
           newPhyloGraph = if (graphType inGS == SoftWired) then T.multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex delSimple
                           else if (graphType inGS == HardWired) then
                              T.multiTraverseFullyLabelHardWired inGS inData leafGraph startVertex delSimple
                           else error "Unsupported graph type in deleteNetEdge.  Must be soft or hard wired"

       in

       -- trace  ("DNERec: " ++ (show edgeToDelete) ++ " at " ++ (show $ snd6 newPhyloGraph)) (
       -- if not modified retiurn original graph
       
      -- This check seems to be issue with delete not functioning properly
       if not wasModified then [inPhyloGraph]

       -- forcing delete for move
       else if force then 
         -- trace ("DNERec forced") 
         [newPhyloGraph]

       -- regular search not sim anneal/drift
       else if (inSimAnnealParams == Nothing) then

          -- harwired must use fiull optimization cost
          if (graphType inGS) == HardWired then
            -- better
            if (snd6 newPhyloGraph < snd6 inPhyloGraph) then [newPhyloGraph]
            -- not better continue
            else deleteNetEdgeRecursive inGS inData leafGraph inPhyloGraph force inSimAnnealParams (tail edgeToDeleteList)

          -- better for tree / softwired
          -- heuristic no good so just checing result
          else if (snd6 newPhyloGraph < snd6 inPhyloGraph) then 
            -- trace  ("DNERec Better -> " ++ (show $ snd6 newPhyloGraph))
            [newPhyloGraph]
         
          else
            -- need to update edge list for new graph
            -- potentially randomize order of list
            deleteNetEdgeRecursive inGS inData leafGraph inPhyloGraph force inSimAnnealParams (tail edgeToDeleteList)

      -- sim anneal/drift
      else
         let numDone = if (method $ fromJust inSimAnnealParams) == SimAnneal then currentStep $ fromJust inSimAnnealParams
                       else driftChanges $ fromJust inSimAnnealParams
             numMax  = if (method $ fromJust inSimAnnealParams) == SimAnneal then numberSteps $ fromJust inSimAnnealParams
                       else driftMaxChanges $ fromJust inSimAnnealParams

             candidateGraphCost = if (graphType inGS) == HardWired then snd6 newPhyloGraph
                                  else ((snd6 inPhyloGraph) + (heuristicDelta - edgeAddDelta))

             (acceptGraph, nextSAParams) = U.simAnnealAccept inSimAnnealParams (snd6 inPhyloGraph) candidateGraphCost
         in
         if (numDone < numMax) then
            if acceptGraph then [newPhyloGraph]
            else deleteNetEdgeRecursive inGS inData leafGraph newPhyloGraph force nextSAParams (tail edgeToDeleteList)

         -- hit end of SA/Drift
         else [inPhyloGraph]
         -- )

-- | deleteEdge deletes an edge (checking if network) and rediagnoses graph
-- contacts in=out=1 edgfes and removes node, reindexing nodes and edges
-- naive for now
-- force requires reoptimization no matter what--used for net move
-- slipping heuristics for now--awful
deleteNetEdge :: GlobalSettings -> ProcessedData -> DecoratedGraph -> PhylogeneticGraph -> Bool -> LG.Edge -> PhylogeneticGraph
deleteNetEdge inGS inData leafGraph inPhyloGraph force edgeToDelete =
   if LG.isEmpty $ thd6 inPhyloGraph then error "Empty input phylogenetic graph in deleteNetEdge"
   else if not (LG.isNetworkEdge (fst6 inPhyloGraph) edgeToDelete) then error ("Edge to delete: " ++ (show edgeToDelete) ++ " not in graph:\n" ++ (LG.prettify $ fst6 inPhyloGraph))
   else
       -- trace ("DNE: " ++ (show edgeToDelete)) (
       let (delSimple, wasModified)  = deleteNetworkEdge edgeToDelete $ fst6 inPhyloGraph

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
           newPhyloGraph = if (graphType inGS == SoftWired) then T.multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex delSimple
                           else if (graphType inGS == HardWired) then
                              T.multiTraverseFullyLabelHardWired inGS inData leafGraph startVertex delSimple
                           else error "Unsupported graph type in deleteNetEdge.  Must be soft or hard wired"
       in
       --check if deletino modified graph
       if not wasModified then inPhyloGraph

       else if force || (graphType inGS) == HardWired then 
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




-- | heuristicDeleteDelta takes the existing graph, edge to delete,
-- reoptimizes starting nodes of two created edges.  Returns cost delta based on
-- previous and new node resolution caches
-- delete n1 -> n2, create u -> v, u' -> v'
-- assumes original is edge n1 -> n2, u' -> (n2, X), n1 -> (n2,v), u (n1,Y)
heuristicDeleteDelta :: GlobalSettings -> PhylogeneticGraph -> LG.Edge -> (VertexCost, LG.LNode VertexInfo, LG.LNode VertexInfo)
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
          uLabAfter      = POSW.getOutDegree2VertexSoftWired inGS (six6 inPhyloGraph) u (v, vLab) uOtherChild inGraph
          uPrimeLabAfter = POSW.getOutDegree2VertexSoftWired inGS (six6 inPhyloGraph) u' (v', vPrimeLab) uPrimeOtherChild inGraph

          -- cost of resolutions
          (_, uCostBefore) = POSW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uLab)
          (_, uPrimeCostBefore) = POSW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uPrimeLab)
          (_, uCostAfter) = POSW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uLabAfter)
          (_, uPrimeCostAfter) = POSW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uPrimeLabAfter)

          addNetDelta = uCostAfter - uCostBefore +  uPrimeCostAfter - uPrimeCostBefore


      in
      -- this should not happen--should try to crete new edges from children of net edges
      if null (LG.parents inGraph n1) || null (filter (/= n1) $ LG.parents inGraph n2) || null (LG.descendants inGraph n2) || null (filter (/= n2) $ LG.descendants inGraph n1) || null (filter ((/= n2) . fst) $ LG.labDescendants inGraph (u', uPrimeLab)) || null (filter ((/= n1) . fst) $ LG.labDescendants inGraph (u, uLab)) then (infinity, dummyNode, dummyNode)
      -- this should not happen--should try to crete new edges from children of net edges
      else if (length (LG.parents inGraph n1) /= 1) || (length (LG.parents inGraph n2) /= 2) || (length (LG.descendants inGraph n2) /= 1) || (length (LG.descendants inGraph n1) /= 2) then error ("Graph malformation in numbers of parents and children in heuristicDeleteDelta")
      else
         (addNetDelta, (u, uLabAfter), (u', uPrimeLabAfter))

