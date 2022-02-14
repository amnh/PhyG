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
                                ) where

import Types.Types
import qualified ParallelUtilities       as PU
import Control.Parallel.Strategies
import qualified Utilities.LocalGraph    as LG
import Debug.Trace
import qualified GraphOptimization.Traversals as T
import qualified GraphOptimization.PreOrderFunctions as PRE
import qualified GraphOptimization.PostOrderFunctions as POS
import qualified Data.List as L
import qualified Data.Text.Lazy              as TL
import GeneralUtilities
import qualified Graphs.GraphOperations  as GO
import Data.Bits
import Data.Maybe
import qualified Data.Vector as V
import qualified GraphOptimization.Medians as M
import qualified Utilities.Utilities as U

-- Need a "steepest" that takes first better netowrk for add and delete.
-- Choose order in max branch length, root-end, leaf-end, and at random

-- moveAllNetEdges is a wrapper for moveAllNetEdges' allowing for multiple simulated annealing rounds
moveAllNetEdges :: GlobalSettings -> ProcessedData -> Int -> Int -> Bool -> ([PhylogeneticGraph], VertexCost) -> (Maybe SAParams, [PhylogeneticGraph]) -> ([PhylogeneticGraph], Int) 
moveAllNetEdges inGS inData numToKeep counter returnMutated (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
   if inSimAnnealParams == Nothing then moveAllNetEdges' inGS inData numToKeep counter returnMutated (curBestGraphList, curBestGraphCost) (Nothing, inPhyloGraphList)
   else 
      let -- create list of params with unique list of random values for rounds of annealing
         annealingRounds = rounds $ fromJust inSimAnnealParams
         annealParamGraphListPair = zip (U.generateUniqueRandList annealingRounds inSimAnnealParams) (replicate annealingRounds inPhyloGraphList)

         (annealRoundsList, counterList) = unzip (fmap (moveAllNetEdges' inGS inData numToKeep counter returnMutated (curBestGraphList, curBestGraphCost)) annealParamGraphListPair `using` PU.myParListChunkRDS)
      in
      (GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] (concat annealRoundsList) , sum counterList)

-- | moveAllNetEdges' removes each edge and adds an edge to all possible plasses each round 
-- until no better or additional graphs are found
-- call with ([], infinity) [single input graph]
moveAllNetEdges' :: GlobalSettings -> ProcessedData -> Int -> Int -> Bool -> ([PhylogeneticGraph], VertexCost) -> (Maybe SAParams, [PhylogeneticGraph]) -> ([PhylogeneticGraph], Int)
moveAllNetEdges' inGS inData numToKeep counter returnMutated (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
   if null inPhyloGraphList then (take numToKeep curBestGraphList, counter)
   else 
      let firstPhyloGraph = head inPhyloGraphList
          currentCost = min curBestGraphCost (snd6 firstPhyloGraph)
          netEdgeList = LG.labNetEdges (thd6 firstPhyloGraph)
          newGraphList' = concatMap (deleteOneNetAddAll inGS inData numToKeep currentCost firstPhyloGraph) netEdgeList -- parallelizes later
          newGraphList = GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] newGraphList'
          newGraphCost = if (not . null) newGraphList' then snd6 $ head newGraphList
                         else infinity
      in
      if null netEdgeList then trace ("\tNo network edges to move") (inPhyloGraphList, counter)

      -- regular move keeping best
      else if inSimAnnealParams == Nothing then
         if newGraphCost > currentCost then moveAllNetEdges' inGS inData numToKeep (counter + 1) returnMutated (firstPhyloGraph : curBestGraphList, currentCost) (inSimAnnealParams, (tail inPhyloGraphList))
         else if newGraphCost < currentCost then
            -- trace ("Found better move")
            moveAllNetEdges' inGS inData numToKeep (counter + 1) returnMutated (newGraphList, newGraphCost) (inSimAnnealParams, (newGraphList ++ (tail inPhyloGraphList)))

         else 
            -- new graph list contains the input graph if equal and filterd unique already in moveAllNetEdges
            let newCurSameBestList = GO.selectPhylogeneticGraph [("unique", (show numToKeep))] 0 ["unique"] (curBestGraphList ++ newGraphList)
            in
            -- trace (show $ length newCurSameBestList) 
            moveAllNetEdges' inGS inData numToKeep (counter + 1) returnMutated (newCurSameBestList, currentCost) (inSimAnnealParams, (tail inPhyloGraphList))

      -- sim anneal choice
      else 
            let numSteps = numberSteps $ fromJust inSimAnnealParams
                curNumSteps = currentStep $ fromJust inSimAnnealParams
                uniqueGraphList = GO.selectPhylogeneticGraph [("unique", (show numToKeep))] 0 ["unique"] newGraphList'
                annealBestCost = min curBestGraphCost (snd6 $ head uniqueGraphList)
                acceptFirstGraph = U.simAnnealAccept (fromJust inSimAnnealParams) annealBestCost (snd6 $ head uniqueGraphList)
            in
            -- trace ("ACG" ++ (show acceptFirstGraph) ++ " " ++ (show $ snd6 $ head uniqueGraphList)) (
            if (curNumSteps < numSteps) then 
               if acceptFirstGraph then
                  moveAllNetEdges' inGS inData numToKeep (counter + 1) returnMutated ((head uniqueGraphList) :  curBestGraphList, annealBestCost) (U.incrementSimAnnealParams   inSimAnnealParams, ((tail uniqueGraphList) ++ (tail inPhyloGraphList)))
               else 
                  moveAllNetEdges' inGS inData numToKeep (counter + 1) returnMutated (curBestGraphList, annealBestCost) (U.incrementSimAnnealParams inSimAnnealParams, ((tail uniqueGraphList) ++ (tail inPhyloGraphList)))

            -- if want non-optimized list for GA or whatever
            else if returnMutated then (take numToKeep curBestGraphList, counter)

            -- optimize list and return
            else 
               let (bestMoveList', counter') =  moveAllNetEdges' inGS inData numToKeep (counter + 1) returnMutated ([], annealBestCost) (Nothing, (take numToKeep curBestGraphList))
                   bestMoveList = GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] bestMoveList'
               in
               --trace ("BM: " ++ (show $ snd6 $ head  bestMoveList))
               (take numToKeep bestMoveList, counter')
      -- )
      

-- | (curBestGraphList, annealBestCost) is a wrapper for moveAllNetEdges' allowing for multiple simulated annealing rounds
insertAllNetEdges :: GlobalSettings -> ProcessedData -> Int -> Int -> Bool -> ([PhylogeneticGraph], VertexCost) -> (Maybe SAParams, [PhylogeneticGraph]) -> ([PhylogeneticGraph], Int)
insertAllNetEdges inGS inData numToKeep counter returnMutated (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
   if inSimAnnealParams == Nothing then 
      insertAllNetEdges' inGS inData numToKeep counter returnMutated (curBestGraphList, curBestGraphCost) (Nothing, inPhyloGraphList)
   else 
      let -- create list of params with unique list of random values for rounds of annealing
         annealingRounds = rounds $ fromJust inSimAnnealParams
         annealParamGraphListPair = zip (U.generateUniqueRandList annealingRounds inSimAnnealParams) (replicate annealingRounds inPhyloGraphList)

         (annealRoundsList, counterList) = unzip (fmap (insertAllNetEdges' inGS inData numToKeep counter returnMutated (curBestGraphList, curBestGraphCost)) annealParamGraphListPair `using` PU.myParListChunkRDS)
      in
      (GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] (concat annealRoundsList) , sum counterList)

-- | insertAllNetEdges' adds network edges one each each round until no better or additional 
-- graphs are found
-- call with ([], infinity) [single input graph]
insertAllNetEdges' :: GlobalSettings -> ProcessedData -> Int -> Int -> Bool -> ([PhylogeneticGraph], VertexCost) -> (Maybe SAParams, [PhylogeneticGraph]) -> ([PhylogeneticGraph], Int)
insertAllNetEdges' inGS inData numToKeep counter returnMutated (curBestGraphList, curBestGraphCost) (inSimAnnealParams,inPhyloGraphList) =
   if null inPhyloGraphList then (take numToKeep curBestGraphList, counter)
   else
      let currentCost = min curBestGraphCost (snd6 $ head inPhyloGraphList)
          (newGraphList, _) = insertEachNetEdge inGS inData numToKeep Nothing (head inPhyloGraphList)
          bestNewGraphList = GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] newGraphList
          newGraphCost = if (not . null) bestNewGraphList then snd6 $ head bestNewGraphList
                         else infinity

      in
      if null newGraphList then (take numToKeep curBestGraphList, counter)

      -- regular insert keeping best
      else if inSimAnnealParams == Nothing then

         -- worse graphs found--go on 
         if  newGraphCost > currentCost then 
            if length curBestGraphList < numToKeep then 
               insertAllNetEdges' inGS inData numToKeep (counter + 1) returnMutated ((head inPhyloGraphList) : curBestGraphList, currentCost) (inSimAnnealParams, (tail inPhyloGraphList))
            else 
               insertAllNetEdges' inGS inData numToKeep (counter + 1) returnMutated (curBestGraphList, currentCost) (inSimAnnealParams, (tail inPhyloGraphList))

         -- "steepest style descent" abandons existing list if better cost found
         else if newGraphCost < currentCost then 
            -- trace ("Found a new edge to insert")
            insertAllNetEdges' inGS inData numToKeep (counter + 1) returnMutated (newGraphList, newGraphCost) (inSimAnnealParams, (newGraphList ++ (tail inPhyloGraphList)))

         -- equal cost
         -- not sure if should add new graphs to queue to do edge deletion again
         else 
            -- new grapjh list contains the input graph if equal and filterd unique already in insertAllNetEdges
            if length curBestGraphList < numToKeep then 
               let newCurSameBestList =  take numToKeep $  GO.selectPhylogeneticGraph [("unique", (show numToKeep))] 0 ["unique"] (curBestGraphList ++ newGraphList)
               in
               insertAllNetEdges' inGS inData numToKeep (counter + 1) returnMutated (newCurSameBestList, currentCost) (inSimAnnealParams, (tail inPhyloGraphList))
            else insertAllNetEdges' inGS inData numToKeep (counter + 1) returnMutated (curBestGraphList, curBestGraphCost) (inSimAnnealParams, (tail inPhyloGraphList))

      -- simulated annealing 
      else 
            -- sim anneal choice
            let numSteps = numberSteps $ fromJust inSimAnnealParams
                curNumSteps = currentStep $ fromJust inSimAnnealParams
                annealBestCost = min curBestGraphCost (snd6 $ head newGraphList)
                acceptFirstGraph = U.simAnnealAccept (fromJust inSimAnnealParams) annealBestCost (snd6 $ head newGraphList)
            in
            -- trace ("ACG" ++ (show acceptFirstGraph) ++ " " ++ (show $ snd6 $ head uniqueGraphList)) (
            if (curNumSteps < numSteps) then 
               if acceptFirstGraph then
                  insertAllNetEdges' inGS inData numToKeep (counter + 1) returnMutated ((head newGraphList) :  curBestGraphList, annealBestCost) (U.incrementSimAnnealParams inSimAnnealParams, ((tail newGraphList) ++ (tail inPhyloGraphList)))
               else 
                  insertAllNetEdges' inGS inData numToKeep (counter + 1) returnMutated (curBestGraphList, annealBestCost) (U.incrementSimAnnealParams inSimAnnealParams, ((tail newGraphList) ++ (tail inPhyloGraphList)))

            -- if want non-optimized list for GA or whatever
            else if returnMutated then (take numToKeep curBestGraphList, counter)

            -- optimize list and return
            else 
               let (bestInsertList', counter') =  insertAllNetEdges' inGS inData numToKeep counter returnMutated ([], annealBestCost) (Nothing, curBestGraphList)
                   bestInsertList = GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] bestInsertList'
               in
               --trace ("BM: " ++ (show $ snd6 $ head  bestMoveList))
               (take numToKeep bestInsertList, counter')
      



-- | (curBestGraphList, annealBestCost) is a wrapper for moveAllNetEdges' allowing for multiple simulated annealing rounds
deleteAllNetEdges :: GlobalSettings -> ProcessedData -> Int -> Int -> Bool -> ([PhylogeneticGraph], VertexCost) -> (Maybe SAParams, [PhylogeneticGraph]) -> ([PhylogeneticGraph], Int)
deleteAllNetEdges inGS inData numToKeep counter returnMutated (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
   if inSimAnnealParams == Nothing then 
      deleteAllNetEdges' inGS inData numToKeep counter returnMutated (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList)
   else 
      let -- create list of params with unique list of random values for rounds of annealing
         annealingRounds = rounds $ fromJust inSimAnnealParams
         annealParamGraphListPair = zip (U.generateUniqueRandList annealingRounds inSimAnnealParams) (replicate annealingRounds inPhyloGraphList)

         (annealRoundsList, counterList) = unzip (fmap (deleteAllNetEdges' inGS inData numToKeep counter returnMutated (curBestGraphList, curBestGraphCost)) annealParamGraphListPair `using` PU.myParListChunkRDS)
      in
      (GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] (concat annealRoundsList) , sum counterList)

-- | deleteAllNetEdges deletes network edges one each each round until no better or additional 
-- graphs are found
-- call with ([], infinity) [single input graph]
deleteAllNetEdges' :: GlobalSettings -> ProcessedData -> Int -> Int -> Bool -> ([PhylogeneticGraph], VertexCost) -> (Maybe SAParams, [PhylogeneticGraph])-> ([PhylogeneticGraph], Int)
deleteAllNetEdges' inGS inData numToKeep counter returnMutated (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
   if null inPhyloGraphList then (take numToKeep curBestGraphList, counter)
   else
      let currentCost = min curBestGraphCost (snd6 $ head inPhyloGraphList)

          (newGraphList', newGraphCost') = deleteEachNetEdge inGS inData numToKeep False (head inPhyloGraphList)
          newGraphList = GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] newGraphList'
          newGraphCost = if (not . null) newGraphList then snd6 $ head newGraphList
                         else infinity

      in
      if null newGraphList then (take numToKeep curBestGraphList, counter)

         -- regular delte wihtout simulated annealing
         -- worse graphs found--go on 
      else if inSimAnnealParams == Nothing then
         if  newGraphCost > currentCost then deleteAllNetEdges' inGS inData numToKeep (counter + 1) returnMutated ((head inPhyloGraphList) : curBestGraphList, currentCost) (inSimAnnealParams, (tail inPhyloGraphList)) 

         -- "steepest style descent" abandons existing list if better cost found
         else if newGraphCost < currentCost then deleteAllNetEdges' inGS inData numToKeep (counter + 1) returnMutated (newGraphList, newGraphCost) (inSimAnnealParams, (newGraphList ++ (tail inPhyloGraphList)))

         -- equal cost
         -- not sure if should add new graphs to queue to do edge deletion again
         else 
            -- new grapjh list contains the input graph if equal and filterd unique already in deleteEachNetEdge
            let newCurSameBestList =  GO.selectPhylogeneticGraph [("unique", (show numToKeep))] 0 ["unique"] (curBestGraphList ++ newGraphList)
            in
            deleteAllNetEdges inGS inData numToKeep  (counter + 1) returnMutated (newCurSameBestList, currentCost) (inSimAnnealParams, (tail inPhyloGraphList))

      -- simulated annealing
      else  
         -- sim anneal choice
            let numSteps = numberSteps $ fromJust inSimAnnealParams
                curNumSteps = currentStep $ fromJust inSimAnnealParams
                annealBestCost = min curBestGraphCost (snd6 $ head newGraphList')
                acceptFirstGraph = U.simAnnealAccept (fromJust inSimAnnealParams) annealBestCost (snd6 $ head newGraphList')
            in
            -- trace ("ACG" ++ (show acceptFirstGraph) ++ " " ++ (show $ snd6 $ head uniqueGraphList)) (
            if (curNumSteps < numSteps) then 
               if acceptFirstGraph then
                  deleteAllNetEdges' inGS inData numToKeep (counter + 1) returnMutated ((head newGraphList') :  curBestGraphList, annealBestCost) (U.incrementSimAnnealParams inSimAnnealParams, ((tail newGraphList') ++ (tail inPhyloGraphList)))
               else 
                  deleteAllNetEdges' inGS inData numToKeep (counter + 1) returnMutated (curBestGraphList, annealBestCost) (U.incrementSimAnnealParams inSimAnnealParams, ((tail newGraphList') ++ (tail inPhyloGraphList)))

            -- if want non-optimized list for GA or whatever
            else if returnMutated then (take numToKeep curBestGraphList, counter)

            -- optimize list and return
            else 
               let (besDeleteList', counter') =  deleteAllNetEdges' inGS inData numToKeep counter returnMutated ([], annealBestCost) (Nothing, curBestGraphList)
                   bestDeleteList = GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] besDeleteList'
               in
               --trace ("BM: " ++ (show $ snd6 $ head  bestMoveList))
               (take numToKeep bestDeleteList, counter')

      
     
-- | deleteOneNetAddAll deletes the specified edge from a graph--creating a fully optimized new one--then readds
-- and keeps best based on delta, reoptimizes those and compares to the oringal cost
-- if better or ssame keeps as per usual
deleteOneNetAddAll :: GlobalSettings -> ProcessedData -> Int -> VertexCost -> PhylogeneticGraph -> LG.LEdge EdgeInfo -> [PhylogeneticGraph]
deleteOneNetAddAll inGS inData numToKeep currentCost inPhyloGraph edgeToDelete = 
   if LG.isEmpty $ thd6 inPhyloGraph then error "Empty graph in deleteOneNetAddAll"
   else
      -- True to force reoptimization of delete
      let deletedEdgeGraph = deleteNetEdge inGS inData inPhyloGraph True (LG.toEdge edgeToDelete)
          (insertedGraphList, minNewCost) = insertEachNetEdge inGS inData numToKeep (Just $ snd6 inPhyloGraph) deletedEdgeGraph
      in
      -- if minNewCost <= currentCost then GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] insertedGraphList
      -- else [] 
      --return all unique-- filtered later
      filter ((/= infinity) . snd6) $ GO.selectPhylogeneticGraph [("unique", (show numToKeep))] 0 ["unique"] insertedGraphList



-- | insertEachNetEdge takes a phylogenetic graph and inserts all permissible network edges one at time
-- and returns unique list of new Phylogenetic Graphs and cost
-- even if worse--could be used for simulated annealing later
-- if equal returns unique graph list
insertEachNetEdge :: GlobalSettings -> ProcessedData -> Int -> Maybe VertexCost -> PhylogeneticGraph -> ([PhylogeneticGraph], VertexCost)
insertEachNetEdge inGS inData numToKeep preDeleteCost inPhyloGraph =
   if LG.isEmpty $ thd6 inPhyloGraph then error "Empty input insertEachNetEdge graph in deleteAllNetEdges"
   else
      let currentCost = if preDeleteCost == Nothing then snd6 inPhyloGraph
                        else fromJust preDeleteCost

          candidateNetworkEdgeList = getPermissibleEdgePairs (thd6 inPhyloGraph)

          -- newGraphList = concat (fmap (insertNetEdgeBothDirections inGS inData inPhyloGraph) candidateNetworkEdgeList `using`  PU.myParListChunkRDS)
          newGraphList = filter (/= emptyPhylogeneticGraph) (fmap (insertNetEdge inGS inData inPhyloGraph preDeleteCost) candidateNetworkEdgeList `using`  PU.myParListChunkRDS)

          --minCostGraphList = GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] newGraphList
          minCost = if null candidateNetworkEdgeList || null newGraphList then infinity
                    else minimum $ fmap snd6 newGraphList
      in
      trace ("\tExamining " ++ (show $ length candidateNetworkEdgeList) ++ " candidate edge pairs") (
      
      -- no network edges to insert
      if null candidateNetworkEdgeList then ([inPhyloGraph], currentCost)
      -- filter later 
      else
         -- (minCostGraphList, minCost)
         (GO.selectPhylogeneticGraph [("unique", (show numToKeep))] 0 ["unique"] $ newGraphList, minCost)
      )

-- | getPermissibleEdgePairs takes a DecoratedGraph and returns the list of all pairs
-- of edges that can be joined by a network edge and meet all necessary conditions
getPermissibleEdgePairs :: DecoratedGraph -> [(LG.LEdge EdgeInfo, LG.LEdge EdgeInfo)]
getPermissibleEdgePairs inGraph = 
   if LG.isEmpty inGraph then error "Empty input graph in isEdgePairPermissible"
   else 
       let edgeList = LG.labEdges inGraph
           edgePairs = cartProd edgeList edgeList
           contraintList = LG.getGraphCoevalConstraints inGraph
           edgeTestList = fmap (isEdgePairPermissible inGraph contraintList) edgePairs `using`  PU.myParListChunkRDS
           pairList = fmap fst $ filter ((== True) . snd) $ zip edgePairs edgeTestList
       in
       -- trace ("Edge Pair list :" ++ (show $ fmap f pairList) ++ "\n"
       --  ++ "GPEP\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inGraph))
       pairList
       where f (a, b) = (LG.toEdge a, LG.toEdge b)

-- | isEdgePairPermissible takes a graph and two edges, coeval contraints, and tests whether a
-- pair of edges can be linked by a new edge and satify three consitions:
--    1) neither edge is a network edge
--    2) one edge cannot be "before" while the other is "after" in any of the constraint pairs
--    3) neither neither edge is an ancestor or descndent edge of the other (tested via bv of nodes)
-- the result should apply to a new edge in either direction
isEdgePairPermissible :: DecoratedGraph -> [([LG.LEdge EdgeInfo],[LG.LEdge EdgeInfo])] -> (LG.LEdge EdgeInfo, LG.LEdge EdgeInfo) -> Bool
isEdgePairPermissible inGraph constraintList (edge1@(u,v,_), edge2@(u',v',_)) =
   if LG.isEmpty inGraph then error "Empty input graph in isEdgePairPermissible"
   else 
       if u == u' then False
       else if v == v' then False
       -- equality implied in above two 
       -- else if LG.toEdge edge1 == LG.toEdge edge2 then False
       else if (LG.isNetworkNode inGraph u) || (LG.isNetworkNode inGraph u') then False
       else if (LG.isNetworkLabEdge inGraph edge1) || (LG.isNetworkLabEdge inGraph edge2) then False
       else if not (LG.meetsAllCoevalConstraints constraintList edge1 edge2) then False
       else if (isAncDescEdge inGraph edge1 edge2) then False
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

-- | insertNetEdgeBothDirections calls insertNetEdge for both u -> v and v -> u new edge orientations
insertNetEdgeBothDirections :: GlobalSettings -> ProcessedData -> PhylogeneticGraph ->  Maybe VertexCost -> (LG.LEdge b, LG.LEdge b) -> [PhylogeneticGraph]
insertNetEdgeBothDirections inGS inData inPhyloGraph preDeleteCost (u,v) = fmap (insertNetEdge inGS inData inPhyloGraph preDeleteCost) [(u,v), (v,u)]


-- | heuristicAddDelta takes teh existing graph, edge pair, and new nodes to create and makes
-- the new nodes and reoprtimizes starting nodes of two edges.  Returns cost delta based on 
-- previous and new node resolution caches
-- returns cost delta and the reoptimized nodes for use in incremental optimization
-- creates (v', n2), (n2, X)u' new, (v, n2)n1  (n1,Y) new  
-- results of new edge n1 -> n2, u' -> (n2, X), n1 -> (n2,v), u (n1,Y)
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
          n2Lab          = POS.getOutDegree1VertexSoftWired n2 vPrimeLab (thd6 inPhyloGraph) [n2]
          uPrimeLabAfter = POS.getOutDegree2VertexSoftWired inGS (six6 inPhyloGraph) u' (n2, n2Lab) uPrimeOtherChild (thd6 inPhyloGraph)
          n1Lab          = POS.getOutDegree2VertexSoftWired inGS (six6 inPhyloGraph) n1 (v, vLab) (n2, n2Lab) (thd6 inPhyloGraph)
          uLabAfter      = POS.getOutDegree2VertexSoftWired inGS (six6 inPhyloGraph) u uOtherChild (n1, n1Lab) (thd6 inPhyloGraph)

          -- cost of resolutions
          (_, uCostBefore) = POS.extractDisplayTrees (Just (-1)) False (vertexResolutionData uLab)
          (_, uPrimeCostBefore) = POS.extractDisplayTrees (Just (-1)) False (vertexResolutionData uPrimeLab)
          (_, uCostAfter) = POS.extractDisplayTrees (Just (-1)) False (vertexResolutionData uLabAfter)
          (_, uPrimeCostAfter) = POS.extractDisplayTrees (Just (-1)) False (vertexResolutionData uPrimeLabAfter)

          addNetDelta = uCostAfter - uCostBefore +  uPrimeCostAfter - uPrimeCostBefore 


      in
      if null (filter ((/= v') . fst) $ LG.labDescendants (thd6 inPhyloGraph) (u', uPrimeLab)) || null (filter ((/= v) . fst) $ LG.labDescendants (thd6 inPhyloGraph) (u, uLab)) then (infinity, dummyNode, dummyNode, dummyNode, dummyNode)
      -- this should not happen--should try to crete new edges from children of net edges
      else if (length $ LG.descendants (thd6 inPhyloGraph) u) < 2 ||  (length $ LG.descendants (thd6 inPhyloGraph) u') < 2 then error ("Outdegree 1 nodes in heuristicAddDelta")
      else 
         (addNetDelta, (u, uLabAfter), (u', uPrimeLabAfter), (n1, n1Lab), (n2, n2Lab))
           


-- | deltaPenaltyAdjustment takes number of leaves and Phylogenetic graph and returns a heuristic graph penalty for adding a single network edge
-- if Wheeler2015Network, this is based on a all changes affecting a single block (most permissive)  and Wheeler 2015 calcualtion of penalty
-- if PMDLGraph -- KMDL not yet implemented
-- if NoNetworkPenalty then 0
deltaPenaltyAdjustment :: GraphFactor -> Int -> PhylogeneticGraph -> VertexCost
deltaPenaltyAdjustment edgeCostModel numLeaves inGraph =
   if edgeCostModel == NoNetworkPenalty then 0.0
   else if edgeCostModel == Wheeler2015Network then
      let graphCost = snd6 inGraph -- this includes any existing penalties--would be better not to include
          numBlocks = V.length $ fth6 inGraph
      in
      graphCost / (fromIntegral $ numBlocks * 2 * ((2 * numLeaves) - 2))

   else if edgeCostModel == PMDLGraph then error ("PMDLGraph  not yet implemented")
   else error ("Network edge cost model not yet implemented: " ++ (show edgeCostModel))



-- | deleteEachNetEdge takes a phylogenetic graph and deletes all network edges one at time
-- and returns best list of new Phylogenetic Graphs and cost
-- even if worse--could be used for simulated annealing later
-- if equal returns unique graph list
deleteEachNetEdge :: GlobalSettings -> ProcessedData -> Int -> Bool ->  PhylogeneticGraph -> ([PhylogeneticGraph], VertexCost)
deleteEachNetEdge inGS inData numToKeep force inPhyloGraph =
   if LG.isEmpty $ thd6 inPhyloGraph then ([], infinity) -- error "Empty input phylogenetic graph in deleteAllNetEdges"
   else
      let currentCost = snd6 inPhyloGraph
          networkEdgeList = LG.netEdges $ thd6 inPhyloGraph
          newGraphList = fmap (deleteNetEdge inGS inData inPhyloGraph force) networkEdgeList `using`  PU.myParListChunkRDS
          uniqueCostGraphList = filter ((/= infinity) . snd6) $ GO.selectPhylogeneticGraph [("unique", (show numToKeep))] 0 ["unique"] newGraphList
          minCost = if null uniqueCostGraphList then infinity 
                    else minimum $ fmap snd6 uniqueCostGraphList
      in
      -- no network edges to delete
      if null networkEdgeList then trace ("\tNo network edges to delete") ([inPhyloGraph], currentCost)
      else
         -- filter later
         (GO.selectPhylogeneticGraph [("unique", (show numToKeep))] 0 ["unique"] $ uniqueCostGraphList, currentCost)



-- | deleteEdge deletes an edge (checking if network) and rediagnoses graph
-- contacts in=out=1 edgfes and removes node, reindexing nodes and edges
-- but only returns simple graph field 
deleteNetEdgeSimple :: GlobalSettings -> ProcessedData -> PhylogeneticGraph -> Bool -> LG.Edge -> PhylogeneticGraph
deleteNetEdgeSimple inGS inData inPhyloGraph force edgeToDelete =
 if LG.isEmpty $ thd6 inPhyloGraph then error "Empty input phylogenetic graph in deleteNetEdge"
   else if not (LG.isNetworkEdge (fst6 inPhyloGraph) edgeToDelete) then error ("Edge to delete: " ++ (show edgeToDelete) ++ " not in graph:\n" ++ (LG.prettify $ fst6 inPhyloGraph))
   else 
       let delSimple = GO.contractIn1Out1EdgesRename $ LG.delEdge edgeToDelete $ fst6 inPhyloGraph 
       in
       -- 3rd field out of sync byut may be needed for leaf graph later
       (delSimple, snd6 inPhyloGraph, thd6 inPhyloGraph, fth6 inPhyloGraph, fft6 inPhyloGraph, six6 inPhyloGraph)


-- | insertNetEdge inserts an edge between two other edges, creating 2 new nodes and rediagnoses graph
-- contacts deletes 2 orginal edges and adds 2 nodes and 5 new edges
-- does not check any edge reasonable-ness properties
-- new edge directed from first to second edge
-- naive for now
-- predeletecost ofr edge move
insertNetEdge :: GlobalSettings -> ProcessedData -> PhylogeneticGraph -> Maybe VertexCost -> (LG.LEdge b, LG.LEdge b) -> PhylogeneticGraph
insertNetEdge inGS inData inPhyloGraph preDeleteCost edgePair@((u,v, _), (u',v', _)) =
   -- trace ("InsertEdge " ++ (show ((u,v), (u',v'))) ++ " into:\n " ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph $ thd6 inPhyloGraph)) (
   if LG.isEmpty $ thd6 inPhyloGraph then error "Empty input phylogenetic graph in insNetEdge"
   else 
       let inSimple = fst6 inPhyloGraph
           numNodes = length $ LG.nodes inSimple
           newNodeOne = (numNodes, TL.pack ("HTU" ++ (show numNodes)))
           newNodeTwo = (numNodes + 1, TL.pack ("HTU" ++ (show $ numNodes + 1)))
           newEdgeList = [(u, fst newNodeOne, 0.0),(fst newNodeOne, v, 0.0),(u', fst newNodeTwo, 0.0),(fst newNodeTwo, v', 0.0),(fst newNodeOne, fst newNodeTwo, 0.0)]
           newSimple = LG.insEdges newEdgeList $ LG.delEdges [(u,v), (u',v')] $ LG.insNodes [newNodeOne, newNodeTwo] inSimple
           leafGraph = LG.extractLeafGraph $ thd6 inPhyloGraph

           -- do not prune other edges if now unused
           pruneEdges = False

           -- don't warn that edges are being pruned
           warnPruneEdges = False

           -- graph optimization from root
           startVertex = Nothing

           
           -- full two-pass optimization
           newPhyloGraph = T.multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex newSimple
            
           -- calculates heursitic graph delta
           (heuristicDelta, _, _, _, _)  = heuristicAddDelta inGS inPhyloGraph edgePair (fst newNodeOne) (fst newNodeTwo) 
                                          

           edgeAddDelta = deltaPenaltyAdjustment (graphFactor inGS) (V.length $ fst3 inData) inPhyloGraph
                          


       in
       -- trace ("INE Deltas: " ++ (show (heuristicDelta, edgeAddDelta)) ++ " preDelete " ++ (show preDeleteCost)
       --  ++ "New Nodes " ++ (show [newNodeOne, newNodeTwo]) ++ " delete edges " ++ (show [(u,v), (u',v')]) ++ " New edges " ++ (show newEdgeList)
       --  ++ "\nInGraph:\n" ++ (LG.prettify inSimple) ++ "\nNewGraph:\n" ++ (LG.prettify newSimple) ) (

       -- Check for cycles
       if LG.cyclic newSimple then 
         -- trace ("INE cyclic") 
         emptyPhylogeneticGraph

       else if GO.parentInChain newSimple then 
         --trace ("INE PIN") 
         emptyPhylogeneticGraph

       -- force evaluation for HardWired
       else if (graphType inGS) == HardWired then newPhyloGraph
       
       -- preDelete cost changes criterion for edge move
       else if preDeleteCost == Nothing then 
          if heuristicDelta + edgeAddDelta < 0 then newPhyloGraph
          else emptyPhylogeneticGraph
          
       else 
         -- no net add cost becasue the numbe rof net nodes is unchaned in add/delete when preDelete cost /= Noting
          if heuristicDelta + (snd6 inPhyloGraph) <= fromJust preDeleteCost then newPhyloGraph
          else emptyPhylogeneticGraph
       --   )

-- | deleteEdge deletes an edge (checking if network) and rediagnoses graph
-- contacts in=out=1 edgfes and removes node, reindexing nodes and edges
-- naive for now
-- force requires reoptimization no matter what--used for net move
deleteNetEdge :: GlobalSettings -> ProcessedData -> PhylogeneticGraph -> Bool -> LG.Edge -> PhylogeneticGraph
deleteNetEdge inGS inData inPhyloGraph force edgeToDelete =
   if LG.isEmpty $ thd6 inPhyloGraph then error "Empty input phylogenetic graph in deleteNetEdge"
   else if not (LG.isNetworkEdge (fst6 inPhyloGraph) edgeToDelete) then error ("Edge to delete: " ++ (show edgeToDelete) ++ " not in graph:\n" ++ (LG.prettify $ fst6 inPhyloGraph))
   else 
       let delSimple = GO.contractIn1Out1EdgesRename $ LG.delEdge edgeToDelete $ fst6 inPhyloGraph
           leafGraph = LG.extractLeafGraph $ thd6 inPhyloGraph

           -- prune other edges if now unused
           pruneEdges = True

           -- don't warn that edges are being pruned
           warnPruneEdges = False

           -- graph optimization from root
           startVertex = Nothing

           (heuristicDelta, _, _) = heuristicDeleteDelta inGS inPhyloGraph edgeToDelete
                                    

           edgeAddDelta = deltaPenaltyAdjustment (graphFactor inGS) (V.length $ fst3 inData) inPhyloGraph
                          
           
           -- full two-pass optimization
           newPhyloGraph = if (graphType inGS == SoftWired) then T.multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex delSimple
                           else if (graphType inGS == HardWired) then 
                              -- trace ("Delete delSimple\n:" ++ (LG.prettify delSimple)) 
                              if (not . LG.cyclic) delSimple then T.multiTraverseFullyLabelHardWired inGS inData leafGraph startVertex delSimple
                              else emptyPhylogeneticGraph
                           else error "Unsupported graph type in deleteNetEdge.  Must be soft or hard wired"
       in
       if force || (graphType inGS) == HardWired then newPhyloGraph
       else if heuristicDelta - edgeAddDelta < 0 then newPhyloGraph
       else emptyPhylogeneticGraph
       

-- | heuristicDeleteDelta takes the existing graph, edge to delete, 
-- reoptimizes starting nodes of two created edges.  Returns cost delta based on 
-- previous and new node resolution caches
-- delete n1 -> n2, create u -> v, u' -> v'
-- assumes original is edge n1 -> n2, u' -> (n2, X), n1 -> (n2,v), u (n1,Y)
heuristicDeleteDelta :: GlobalSettings -> PhylogeneticGraph -> LG.Edge -> (VertexCost, LG.LNode VertexInfo, LG.LNode VertexInfo)
heuristicDeleteDelta inGS inPhyloGraph (n1, n2) =
  if LG.isEmpty (fst6 inPhyloGraph) then error "Empty graph in heuristicAddDelta"
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
          uLabAfter      = POS.getOutDegree2VertexSoftWired inGS (six6 inPhyloGraph) u (v, vLab) uOtherChild inGraph
          uPrimeLabAfter = POS.getOutDegree2VertexSoftWired inGS (six6 inPhyloGraph) u' (v', vPrimeLab) uPrimeOtherChild inGraph

          -- cost of resolutions
          (_, uCostBefore) = POS.extractDisplayTrees (Just (-1)) False (vertexResolutionData uLab)
          (_, uPrimeCostBefore) = POS.extractDisplayTrees (Just (-1)) False (vertexResolutionData uPrimeLab)
          (_, uCostAfter) = POS.extractDisplayTrees (Just (-1)) False (vertexResolutionData uLabAfter)
          (_, uPrimeCostAfter) = POS.extractDisplayTrees (Just (-1)) False (vertexResolutionData uPrimeLabAfter)

          addNetDelta = uCostAfter - uCostBefore +  uPrimeCostAfter - uPrimeCostBefore 


      in
      -- this should not happen--should try to crete new edges from children of net edges
      if null (LG.parents inGraph n1) || null (filter (/= n1) $ LG.parents inGraph n2) || null (LG.descendants inGraph n2) || null (filter (/= n2) $ LG.descendants inGraph n1) || null (filter ((/= n2) . fst) $ LG.labDescendants inGraph (u', uPrimeLab)) || null (filter ((/= n1) . fst) $ LG.labDescendants inGraph (u, uLab)) then (infinity, dummyNode, dummyNode)
      -- this should not happen--should try to crete new edges from children of net edges
      else if (length (LG.parents inGraph n1) /= 1) || (length (LG.parents inGraph n2) /= 2) || (length (LG.descendants inGraph n2) /= 1) || (length (LG.descendants inGraph n1) /= 2) then error ("Graph malformation in numbersof parents and children in heuristicDeleteDelta")
      else 
         (addNetDelta, (u, uLabAfter), (u', uPrimeLabAfter))
           