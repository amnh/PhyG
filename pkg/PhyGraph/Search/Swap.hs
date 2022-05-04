{- |
Module      :  Swap.hs
Description :  Module specifying graph swapping rearrangement functions
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

module Search.Swap  ( swapSPRTBR
                    , reoptimizeSplitGraphFromVertexTuple
                    , rejoinGraphKeepBestTuple
                    ) where

import           Control.Parallel.Strategies
import qualified Data.List                            as L
import           Data.Maybe
import qualified Data.Vector                          as V
import           Debug.Trace
import           GeneralUtilities
import qualified GraphOptimization.Medians            as M
import qualified GraphOptimization.PostOrderFunctions as POS
import qualified GraphOptimization.PreOrderFunctions  as PRE
import qualified GraphOptimization.Traversals         as T
import qualified Graphs.GraphOperations               as GO
import qualified ParallelUtilities                    as PU
import           Types.Types
import qualified Utilities.LocalGraph                 as LG
import           Utilities.Utilities                  as U


-- | swapSPRTBR perfomrs SPR or TBR branch (edge) swapping on graphs
-- runs both SPR and TBR depending on argument since so much duplicated functionality
-- 'steepest' abandons swap graph and switces to found graph as soon as anyhting 'better'
-- is found. The alternative (all) examines the entire neighborhood and retuns the best result
-- the retuns is a list of better graphs and the number of swapping rounds were required to ge there
swapSPRTBR  :: String
            -> GlobalSettings
            -> ProcessedData
            -> Int
            -> Int
            -> Bool
            -> Bool
            -> Bool
            -> Bool
            -> (Maybe SAParams, PhylogeneticGraph)
            -> ([PhylogeneticGraph], Int)
swapSPRTBR swapType inGS inData numToKeep maxMoveEdgeDist steepest hardwiredSPR doIA returnMutated (inSimAnnealParams, inGraph) =
   -- trace ("In swapSPRTBR:") (
   if LG.isEmpty (fst6 inGraph) then ([], 0)
   else
      let numLeaves = V.length $ fst3 inData
          leafGraph = T.makeSimpleLeafGraph inData
          leafDecGraph = T.makeLeafGraph inData
          leafGraphSoftWired = T.makeLeafGraphSoftWired inData
          charInfoVV = six6 inGraph


          inGraphNetPenalty = if (graphType inGS == Tree) then 0.0
                             else if (graphFactor inGS) == NoNetworkPenalty then 0.0
                             else if (graphFactor inGS) == Wheeler2015Network then T.getW15NetPenalty Nothing inGraph
                             else if (graphType inGS == HardWired) then error ("Graph type not implemented: " ++ (show $ graphType inGS))
                             else error ("Network penalty type " ++ (show $ graphFactor inGS) ++ " is not yet implemented")
          inGraphNetPenaltyFactor = inGraphNetPenalty / (snd6 inGraph)
      in
      -- trace ("SSPRTBR:" ++ (show inGraphNetPenaltyFactor)) (

      if inSimAnnealParams == Nothing then
       -- trace ("Non SA swap") (
      -- steepest takes immediate best--does not keep equall cost
      -- NOthing for SimAnneal Params
          if steepest then
             let (swappedGraphs, counter) = swapSteepest swapType hardwiredSPR inGS inData numToKeep maxMoveEdgeDist True 0 (snd6 inGraph) [] [inGraph] numLeaves leafGraph leafDecGraph leafGraphSoftWired charInfoVV doIA inGraphNetPenaltyFactor Nothing

                 -- swap "all" after steepest descent
                 -- (swappedGraphs', counter') = swapAll swapType hardwiredSPR inGS inData numToKeep maxMoveEdgeDist True counter (snd6 $ head swappedGraphs) [] swappedGraphs numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA inGraphNetPenaltyFactor
             in
             -- trace ("Steepest SSPRTBR: " ++ (show (length swappedGraphs, counter)))
             (swappedGraphs, counter)

          -- All does all swaps before taking best
          else
             -- trace ("Going into SwapAll") (
             let (swappedGraphs, counter) = swapAll swapType hardwiredSPR inGS inData numToKeep maxMoveEdgeDist False 0 (snd6 inGraph) [] [inGraph] numLeaves leafGraph leafDecGraph leafGraphSoftWired charInfoVV doIA inGraphNetPenaltyFactor
             in
             -- trace ("All SSPRTBR: " ++ (show (length swappedGraphs, counter)))
             (swappedGraphs, counter)
             --)
             -- )

      -- simulated annealing/drifting acceptance does a steepest with SA acceptance
      -- then a swap steepest and all on annealed graph
      -- same at this level method (SA, Drift) choice occurs at lower level
      else
         -- annealed should only yield a single graph
         --trace ("\tAnnealing/Drifting Swap") (
         let -- create list of params with unique list of random values for rounds of annealing
             annealingRounds = rounds $ fromJust inSimAnnealParams
             newSimAnnealParamList = U.generateUniqueRandList annealingRounds inSimAnnealParams

             -- this to ensure current step set to 0
             (annealedGraphs', anealedCounter) = unzip $ (fmap (swapSteepest swapType hardwiredSPR inGS inData 1 maxMoveEdgeDist True 0 (snd6 inGraph) [] [inGraph] numLeaves leafGraph leafDecGraph leafGraphSoftWired charInfoVV doIA inGraphNetPenaltyFactor) newSimAnnealParamList `using` PU.myParListChunkRDS)

             annealedGraphs = take numToKeep $ GO.selectPhylogeneticGraph [("unique","")] 0 ["unique"] $ concat annealedGraphs'

             (swappedGraphs, counter) = swapSteepest swapType hardwiredSPR inGS inData numToKeep maxMoveEdgeDist True 0 (min (snd6 inGraph) (snd6 $ head annealedGraphs)) [] annealedGraphs numLeaves leafGraph leafDecGraph leafGraphSoftWired charInfoVV doIA inGraphNetPenaltyFactor Nothing

             -- swap "all" after steepest descent
             (swappedGraphs', counter') = swapAll swapType hardwiredSPR inGS inData numToKeep maxMoveEdgeDist True counter (snd6 $ head swappedGraphs) [] swappedGraphs numLeaves leafGraph leafDecGraph leafGraphSoftWired charInfoVV doIA inGraphNetPenaltyFactor

             uniqueGraphs = take numToKeep $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (inGraph : swappedGraphs')
         in
         -- trace ("Steepest SSPRTBR: " ++ (show (length swappedGraphs, counter)))
         --trace ("AC:" ++ (show $ fmap snd6 $ concat annealedGraphs') ++ " -> " ++ (show $ fmap snd6 $ swappedGraphs')) (

         -- this Bool for Genetic Algorithm mutation step
         if not returnMutated then (uniqueGraphs, counter' + (sum anealedCounter))
         else (annealedGraphs, counter' + (sum anealedCounter))
         -- )
         -- )




-- | swapAll performs branch swapping on all 'break' edges and all readditions
-- edges are unsorted since doing all of them
swapAll  :: String
         -> Bool
         -> GlobalSettings
         -> ProcessedData
         -> Int
         -> Int
         -> Bool
         -> Int
         -> VertexCost
         -> [PhylogeneticGraph]
         -> [PhylogeneticGraph]
         -> Int
         -> SimpleGraph
         -> DecoratedGraph
         -> DecoratedGraph
         -> V.Vector (V.Vector CharInfo)
         -> Bool
         -> VertexCost
         -> ([PhylogeneticGraph], Int)
swapAll swapType hardwiredSPR inGS inData numToKeep maxMoveEdgeDist steepest counter curBestCost curSameBetterList inGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor =
   -- trace ("ALL") (
   if null inGraphList then
      (take numToKeep $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] curSameBetterList, counter)
   else
      let firstGraph = head inGraphList
          firstDecoratedGraph = thd6 firstGraph
          (firstRootIndex, _) = head $ LG.getRoots firstDecoratedGraph

          -- determine edges to break on--'bridge' edges only for network
          -- filter out edges from root since no use--would just rejoin
          breakEdgeList = if (graphType inGS) == Tree then filter ((/= firstRootIndex) . fst3) $ LG.labEdges firstDecoratedGraph
                          else filter ((/= firstRootIndex) . fst3) $ GO.getEdgeSplitList firstDecoratedGraph

          -- create list of breaks
          splitTupleList = fmap (GO.splitGraphOnEdge firstDecoratedGraph) breakEdgeList `using` PU.myParListChunkRDS
          (splitGraphList, graphRootList, prunedGraphRootIndexList,  originalConnectionOfPruned) = L.unzip4 splitTupleList

          reoptimizedSplitGraphList = zipWith3 (reoptimizeSplitGraphFromVertex inGS inData doIA netPenaltyFactor) splitGraphList graphRootList prunedGraphRootIndexList `using` PU.myParListChunkRDS

          -- create rejoins-- adds in break list so don't remake the initial graph
          -- didn't concatMap so can parallelize later
          -- this cost prob doesn't include the root/net penalty--so need to figure out
          swapPairList = concat $ L.zipWith4 (rejoinGraphKeepBest inGS swapType hardwiredSPR curBestCost maxMoveEdgeDist doIA charInfoVV) reoptimizedSplitGraphList prunedGraphRootIndexList originalConnectionOfPruned (fmap head $ fmap ((LG.parents $ thd6 firstGraph).fst3) breakEdgeList)

          -- keeps better heuristic swap costs graphs based on current best as opposed to minimum heuristic costs
          -- minimumCandidateGraphCost = if (null swapPairList) then infinity
          --                             else minimum $ fmap snd swapPairList
          candidateSwapGraphList = if (graphType inGS /= Tree) then
                                        -- checks for cylces--rare but can occur
                                        filter ((== False) . (GO.parentInChain . fst)) $ filter ((== False) . (LG.cyclic . fst)) $ filter ((<= curBestCost). snd) swapPairList
                                   else filter ((<= curBestCost). snd) swapPairList


          -- this should be incremental--full 2-pass for now
          reoptimizedSwapGraphList = if (graphType inGS /= Tree) then
                                        -- checks for anc/desc time violation--can occur with some edge splits
                                        let newGraphs = fmap (T.multiTraverseFullyLabelGraph inGS inData False False Nothing) (fmap fst candidateSwapGraphList) `using` PU.myParListChunkRDS
                                        in filter ((== False) . (GO.hasNetNodeAncestorViolation . thd6)) newGraphs
                                     else fmap (T.multiTraverseFullyLabelGraph inGS inData False False Nothing) (fmap fst candidateSwapGraphList) `using` PU.myParListChunkRDS


          -- selects best graph list based on full optimization
          bestSwapGraphList = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] reoptimizedSwapGraphList

          bestSwapCost = if null breakEdgeList || null swapPairList || null candidateSwapGraphList || null reoptimizedSwapGraphList || null bestSwapGraphList then infinity
                         else snd6 $ head bestSwapGraphList

      in
      -- trace ("Cycles: "  ++ (show $ fmap LG.cyclic $ fmap fst candidateSwapGraphList)) (
      -- trace ("Breakable Edges :" ++ (show $ fmap LG.toEdge breakEdgeList) ++ "\nIn graph:\n" ++ (LG.prettify $ fst6 firstGraph)) (
      -- trace ("(Est, [FP]): " ++ (show minimumCandidateGraphCost) ++ " " ++ (show $ fmap snd6 reoptimizedSwapGraphList)) (
      -- either no better or more of same cost graphs
      -- trace ("BSG: " ++ " simple " ++ (LG.prettify $ fst6 $ head bestSwapGraphList) ++ " Decorated " ++ (LG.prettify $ thd6 $ head bestSwapGraphList) ++ "\nCharinfo\n" ++ (show $ charType $ V.head $ V.head $ six6 $ head bestSwapGraphList)) (
      trace ("All--Choosing what to do: " ++ (show (bestSwapCost, curBestCost, length curSameBetterList, numToKeep)) ++ " " ++ (show (length reoptimizedSplitGraphList, length swapPairList, length candidateSwapGraphList, length reoptimizedSwapGraphList, length bestSwapGraphList))) (
      if bestSwapCost == curBestCost then
         if (length curSameBetterList == numToKeep) then 
            -- trace ("Same cost and maxed out")
            swapAll swapType hardwiredSPR inGS inData numToKeep maxMoveEdgeDist steepest (counter + 1) curBestCost curSameBetterList (tail inGraphList) numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor
         else
             --equality informed by zero-length edges
             let newCurSameBestList = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] (firstGraph : curSameBetterList)
                                      -- if firstGraph `notElem` curSameBetterList then (firstGraph : curSameBetterList)
                                      -- else curSameBetterList
                 graphsToSwap = take numToKeep $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (tail inGraphList)

             in
             -- trace ("Same cost: " ++ (show bestSwapCost) ++ " with " ++ (show $ length $ graphsToSwap) ++ " more to swap and " ++ (show $ length newCurSameBestList)
             --   ++ " graphs in 'best' list " ++ " max size " ++ (show numToKeep) ) -- ++ "\n" ++ (concat prettyBestSwapGraphList))
             swapAll swapType hardwiredSPR inGS inData numToKeep maxMoveEdgeDist steepest (counter + 1) curBestCost newCurSameBestList graphsToSwap numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor

      -- better cost graphs
      else if (bestSwapCost < curBestCost) then
         trace ("\t->" ++ (show bestSwapCost)) --  ++ "\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph $ thd6 $ head bestSwapGraphList))
         swapAll swapType hardwiredSPR inGS inData numToKeep maxMoveEdgeDist steepest (counter + 1) bestSwapCost bestSwapGraphList (bestSwapGraphList ++ (tail inGraphList))  numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor

      -- didn't find equal or better graphs
      else
         -- trace ("Worse cost")(
         let newCurSameBestList = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] (firstGraph : curSameBetterList)
                                  -- if firstGraph `notElem` curSameBetterList then (firstGraph : curSameBetterList)
                                  -- else curSameBetterList
         in
         swapAll swapType hardwiredSPR inGS inData numToKeep maxMoveEdgeDist steepest (counter + 1) curBestCost newCurSameBestList (tail inGraphList) numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor
         -- )
      -- )
      )

-- | swapSteepest performs branch swapping greedily switching to found graph if better
-- infomrs evaluation--less parallelism
swapSteepest   :: String
               -> Bool
               -> GlobalSettings
               -> ProcessedData
               -> Int
               -> Int
               -> Bool
               -> Int
               -> VertexCost
               -> [PhylogeneticGraph]
               -> [PhylogeneticGraph]
               -> Int
               -> SimpleGraph
               -> DecoratedGraph
               -> DecoratedGraph
               -> V.Vector (V.Vector CharInfo)
               -> Bool
               -> VertexCost
               -> Maybe SAParams
               -> ([PhylogeneticGraph], Int)
swapSteepest swapType hardwiredSPR inGS inData numToKeep maxMoveEdgeDist steepest counter curBestCost curSameBetterList inGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor inSimAnnealVals =
   --trace ("steepest") (
   if null inGraphList then
      (take numToKeep $ GO.getBVUniqPhylogeneticGraph True curSameBetterList, counter)
   else
      let firstGraph = head inGraphList
          firstDecoratedGraph = thd6 firstGraph
          (firstRootIndex, _) = head $ LG.getRoots firstDecoratedGraph

          -- filter out edges from root since no use--would just rejoin
          firstEdgeList = filter ((/= firstRootIndex) . fst3) $ LG.labEdges firstDecoratedGraph

          -- determine edges to break on--'bridge' edges only for network
          -- longest edges first first
          breakEdgeList = if (graphType inGS) == Tree then GO.sortEdgeListByLength firstEdgeList
                          else GO.sortEdgeListByLength $ GO.getEdgeSplitList firstDecoratedGraph

          -- create list of breaks
          splitTupleList = fmap (GO.splitGraphOnEdge firstDecoratedGraph) breakEdgeList

          (splitGraphList, graphRootList, prunedGraphRootIndexList,  originalConnectionOfPruned) = L.unzip4 splitTupleList

          reoptimizedSplitGraphList = zipWith3 (reoptimizeSplitGraphFromVertex inGS inData doIA netPenaltyFactor) splitGraphList graphRootList prunedGraphRootIndexList

          -- create rejoins-- reoptimized fully in steepest returns PhylogheneticGraph
          (reoptimizedSwapGraphList, newAnnealVals) = rejoinGraphKeepBestSteepest inGS inData swapType hardwiredSPR curBestCost numToKeep maxMoveEdgeDist True doIA charInfoVV inSimAnnealVals $ L.zip4 reoptimizedSplitGraphList prunedGraphRootIndexList originalConnectionOfPruned (fmap head $ fmap ((LG.parents $ thd6 firstGraph).fst3) breakEdgeList)

          -- this should be incremental--full 2-pass for now
          -- reoptimizedSwapGraph = T.multiTraverseFullyLabelGraph inGS inData False False Nothing $ fst $ head swapPairList


          bestSwapCost = if null breakEdgeList || null reoptimizedSplitGraphList || null reoptimizedSwapGraphList then infinity
                         else snd $ head reoptimizedSwapGraphList

      in
      -- trace ("Breakable Edges :" ++ (show $ fmap LG.toEdge breakEdgeList) ++ "\nIn graph:\n" ++ (LG.prettify $ fst6 firstGraph)) (

      -- either no better or more of same cost graphs
      -- trace ("BSG: " ++ " simple " ++ (LG.prettify $ fst6 $ head bestSwapGraphList) ++ " Decorated " ++ (LG.prettify $ thd6 $ head bestSwapGraphList) ++ "\nCharinfo\n" ++ (show $ charType $ V.head $ V.head $ six6 $ head bestSwapGraphList)) (

      -- trace ("Steepest--Choosing what to do: " ++ (show (bestSwapCost, curBestCost, length curSameBetterList, numToKeep)) ++ " " ++ (show (length reoptimizedSplitGraphList, length reoptimizedSwapGraphList))) (
      -- check that graphs were returned.  If nothing then return in Graph
      if bestSwapCost == infinity then (inGraphList, counter + 1)

      else if inSimAnnealVals == Nothing then
          if (bestSwapCost < curBestCost) then
             --trace ("Steepest better")
             trace ("\t->" ++ (show bestSwapCost))
             swapSteepest swapType hardwiredSPR inGS inData numToKeep maxMoveEdgeDist steepest (counter + 1) bestSwapCost (fmap fst reoptimizedSwapGraphList) (fmap fst reoptimizedSwapGraphList) numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor Nothing

          else if (bestSwapCost == curBestCost) then
            -- buffer full
            if (length curSameBetterList == numToKeep) then
                  swapSteepest swapType hardwiredSPR inGS inData numToKeep maxMoveEdgeDist steepest (counter + 1) curBestCost curSameBetterList (tail inGraphList) numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor inSimAnnealVals
            -- room for another Graph if unnique
            else
               let newCurSameBestList = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] (firstGraph : curSameBetterList)
                                      -- if firstGraph `notElem` curSameBetterList then (firstGraph : curSameBetterList)
                                      -- else curSameBetterList
                   graphsToSwap = take numToKeep $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] ((tail inGraphList) L.\\ newCurSameBestList)
               in
               swapSteepest swapType hardwiredSPR inGS inData numToKeep maxMoveEdgeDist steepest (counter + 1) curBestCost curSameBetterList graphsToSwap numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor inSimAnnealVals

          -- didn't find equal or better graphs
          else (inGraphList, counter + 1)

      -- Simulated annealing
      else
            -- abstract stopping criterion to continue
            let numDone = if (method $ fromJust inSimAnnealVals) == SimAnneal then currentStep $ fromJust newAnnealVals
                          else driftChanges $ fromJust newAnnealVals
                numMax  = if (method $ fromJust inSimAnnealVals) == SimAnneal then numberSteps $ fromJust newAnnealVals
                          else driftMaxChanges $ fromJust newAnnealVals
            in

            --simulated Annealing/Drift check if at end then return
            if (numDone < numMax)  then
                swapSteepest swapType hardwiredSPR inGS inData numToKeep maxMoveEdgeDist steepest (counter + 1) bestSwapCost (fmap fst reoptimizedSwapGraphList) (fmap fst reoptimizedSwapGraphList) numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor newAnnealVals

            else
                (fmap fst reoptimizedSwapGraphList, counter + 1)


      -- )

-- | rejoinGraphKeepBestTuple wrapper for rejoinGraphKeepBest but with last 5 arguments as a tuple
rejoinGraphKeepBestTuple :: GlobalSettings
                            -> String
                            -> Bool
                            -> VertexCost
                            -> Int
                            -> Bool
                            -> V.Vector (V.Vector CharInfo)
                            -> ((DecoratedGraph, VertexCost), LG.Node, LG.Node, LG.Node)
                            -> [(SimpleGraph, VertexCost)]
rejoinGraphKeepBestTuple inGS swapType hardwiredSPR curBestCost  maxMoveEdgeDist  doIA charInfoVV ((splitGraph, splitCost), prunedGraphRootIndex, nakedNode, originalSplitNode) =
   rejoinGraphKeepBest inGS swapType hardwiredSPR curBestCost  maxMoveEdgeDist  doIA charInfoVV (splitGraph, splitCost) prunedGraphRootIndex nakedNode originalSplitNode


-- | rejoinGraphKeepBest rejoins split trees on available edges (non-root, and not original split)
-- if steepest is False does not sort order of edges, other wise sorts in order of closeness to original edge
-- uses delta
-- NNI sorts edges on propinquity taking first 2 edges
-- TBR does the rerooting of pruned subtree
-- originalConnectionOfPruned is the "naked" node that was creted when teh graph was split and will
-- be used for the rejoin node in the middle of th einvaded edge
rejoinGraphKeepBest :: GlobalSettings
                    -> String
                    -> Bool
                    -> VertexCost
                    -> Int
                    -> Bool
                    -> V.Vector (V.Vector CharInfo)
                    -> (DecoratedGraph, VertexCost)
                    -> LG.Node
                    -> LG.Node
                    -> LG.Node
                    -> [(SimpleGraph, VertexCost)]
rejoinGraphKeepBest inGS swapType hardWiredSPR curBestCost  maxMoveEdgeDist  doIA charInfoVV (splitGraph, splitCost) prunedGraphRootIndex nakedNode originalSplitNode =
   -- case where swap split retunred empty because too few nodes in remaining graph to add to
   if LG.isEmpty splitGraph then []
   -- this more for fusing situations
   else if splitCost > curBestCost then []
   else
      let -- outgroupEdges = LG.out splitGraph graphRoot
          (_, prunedSubTreeEdges) = LG.nodesAndEdgesAfter splitGraph [(nakedNode, fromJust $ LG.lab splitGraph nakedNode)]

          --only sort if limited egde rejoins
          splitEdges = if (graphType inGS) == Tree then LG.labEdges splitGraph
                       else GO.getEdgeSplitList splitGraph
          edgesToInvade = if maxMoveEdgeDist == (maxBound :: Int) then splitEdges L.\\ prunedSubTreeEdges -- L.\\ (outgroupEdges ++ prunedSubTreeEdges)
                          else take maxMoveEdgeDist $ L.intersect splitEdges ((GO.sortEdgeListByDistance splitGraph [originalSplitNode] [originalSplitNode]) L.\\ prunedSubTreeEdges)

          prunedGraphRootNode = (prunedGraphRootIndex, fromJust $ LG.lab splitGraph prunedGraphRootIndex)
          (nodesAfterPrunedRoot, edgesInPrunedSubGraph) = LG.nodesAndEdgesAfter splitGraph [prunedGraphRootNode]

          -- pruned subgraph too small to reroot in TBR
          onlySPR = length nodesAfterPrunedRoot < 3


          candidateEditList = concatMap (addSubGraph inGS swapType hardWiredSPR doIA splitGraph prunedGraphRootNode splitCost nakedNode onlySPR edgesInPrunedSubGraph charInfoVV) edgesToInvade `using` PU.myParListChunkRDS


          minCandidateCost = if (not $ null candidateEditList) then minimum $ fmap fst3 candidateEditList
                             else infinity
      in
      -- trace ("RGKB: " ++ (show $ fmap LG.toEdge edgesToInvade) ++ "\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph splitGraph)) (
      if minCandidateCost > curBestCost then []
      else
         let bestEdits = filter ((<= curBestCost). fst3) candidateEditList -- not minimum cancidate cost--better if check all equal or better than curent best
             splitGraphSimple = GO.convertDecoratedToSimpleGraph splitGraph
             swapSimpleGraphList = fmap (applyGraphEdits splitGraphSimple) bestEdits
         in
         zip swapSimpleGraphList (L.replicate (length swapSimpleGraphList) minCandidateCost)
      -- )

-- | rejoinGraphKeepBestSteepest rejoins split trees on available edges (non-root, and not original split)
-- if steepest is False does not sort order of edges, other wise sorts in order of closeness to original edge
-- uses delta
-- NNI sorts edges on propinquity taking first 2 edges
-- TBR does the rerooting of pruned subtree
-- originalConnectionOfPruned is the "naked" node that was creted when teh graph was split and will
-- be used for the rejoin node in the middle of the invaded edge
rejoinGraphKeepBestSteepest :: GlobalSettings
                             -> ProcessedData
                             -> String
                             -> Bool
                             -> VertexCost
                             -> Int
                             -> Int
                             -> Bool
                             -> Bool
                             -> V.Vector (V.Vector CharInfo)
                             -> Maybe SAParams
                             -> [((DecoratedGraph, VertexCost), LG.Node , LG.Node , LG.Node)]
                             -> ([(PhylogeneticGraph, VertexCost)], Maybe SAParams)
rejoinGraphKeepBestSteepest inGS inData swapType hardwiredSPR curBestCost numToKeep maxMoveEdgeDist steepest doIA charInfoVV inSimAnnealVals splitInfoList =
   if null splitInfoList then ([], inSimAnnealVals)
   else
      -- trace ("In rejoin steepes with split node " ++ (show $ fft5 $ head splitInfoList)) (
      let ((splitGraph, splitCost), prunedGraphRootIndex, nakedNode, originalSplitNode) = head splitInfoList
          -- outgroupEdges = LG.out splitGraph graphRoot
          (_, prunedSubTreeEdges) = LG.nodesAndEdgesAfter splitGraph [(nakedNode, fromJust $ LG.lab splitGraph nakedNode)]

          splitEdges = if (graphType inGS) == Tree then LG.labEdges splitGraph
                       else GO.getEdgeSplitList splitGraph

          edgesToInvade = take maxMoveEdgeDist $ L.intersect splitEdges ((GO.sortEdgeListByDistance splitGraph [originalSplitNode] [originalSplitNode]) L.\\ prunedSubTreeEdges) -- L.\\ (outgroupEdges ++ prunedSubTreeEdges)


          prunedGraphRootNode = (prunedGraphRootIndex, fromJust $ LG.lab splitGraph prunedGraphRootIndex)
          (nodesAfterPrunedRoot, edgesInPrunedSubGraph) = LG.nodesAndEdgesAfter splitGraph [prunedGraphRootNode]
          onlySPR = length nodesAfterPrunedRoot < 3

          (candidateGraphList, newAnnealVals) = addSubGraphSteepest inGS inData swapType hardwiredSPR doIA splitGraph prunedGraphRootNode splitCost curBestCost nakedNode onlySPR edgesInPrunedSubGraph charInfoVV inSimAnnealVals edgesToInvade

          -- numLeaves = fromIntegral $ V.length $ fst3 inData
      in
      -- experimental union-type split exclusion
         -- skip rearrangements if split delta too small VAron and Wheeler (2013)
      -- if (curBestCost / numLeaves) > 1.17 * (curBestCost - splitCost) then ([], inSimAnnealVals)

      -- case where swap split returned empty because too few nodes in remaining graph to add to
      -- else
      -- trace ("RGKBS:" ++ (show (length candidateGraphList, fmap snd6 candidateGraphList)) ++ " " ++ (show $ LG.isEmpty splitGraph)) (
      if LG.isEmpty splitGraph || null candidateGraphList then ([], inSimAnnealVals)

      -- normal steepest--only return if better if equal does not return, but will return multiple
      else if inSimAnnealVals == Nothing then
        let candidateGraphList' = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] candidateGraphList
        in
        if (snd6 $ head candidateGraphList') < curBestCost then
            -- ([(head candidateGraphList, snd6 $ head candidateGraphList)], Nothing)
            (zip candidateGraphList' (fmap snd6 candidateGraphList'), Nothing)
            -- ([(head candidateGraphList, snd6 $ head candidateGraphList)], Nothing)

        else rejoinGraphKeepBestSteepest inGS inData swapType hardwiredSPR curBestCost numToKeep maxMoveEdgeDist steepest doIA charInfoVV Nothing (tail splitInfoList)

      -- simulated Annealing/Drifting--recurse if not at max steps
      else
         ([(head candidateGraphList, snd6 $ head candidateGraphList)], newAnnealVals)

      -- )

      -- ))

-- | addSubGraphSteepest "adds" a subtree back into an edge calculating the cost of the graph via the delta of the add and costs of the two components
-- used in "steepest" descendt swapping
addSubGraphSteepest :: GlobalSettings
                     -> ProcessedData
                     -> String
                     -> Bool
                     -> Bool
                     -> DecoratedGraph
                     -> LG.LNode VertexInfo
                     -> VertexCost
                     -> VertexCost
                     -> LG.Node
                     -> Bool
                     -> [LG.LEdge EdgeInfo]
                     -> V.Vector (V.Vector CharInfo)
                     -> Maybe SAParams
                     -> [LG.LEdge EdgeInfo]
                     -> ([PhylogeneticGraph], Maybe SAParams)
addSubGraphSteepest inGS inData swapType hardwiredSPR doIA inGraph prunedGraphRootNode splitCost curBestCost nakedNode onlySPR edgesInPrunedSubGraph charInfoVV inSimAnnealVals targetEdgeList =
   if null targetEdgeList then ([], inSimAnnealVals)
   -- this more for graph fusing checks
   else if (splitCost > curBestCost) && (inSimAnnealVals == Nothing) then ([], Nothing)
   else
      let targetEdge@(eNode, vNode, _) = head targetEdgeList
          -- existingEdgeCost = minLength targetlabel
          edge0 = (nakedNode, vNode, 0.0)
          edge1 = (eNode, nakedNode, 0.0)
          targetEdgeData = M.makeEdgeData doIA inGraph charInfoVV targetEdge
          -- edge2 = (nakedNode, prunedGraphRootIndex, 0.0)
          --newNode = (nakedNode, TL.pack ("HTU" ++ (show nakedNode)))

          -- for SPRT/NNI only need preliminary state of root-pruned node
          -- for TBR ther are multiple verticces created for each edge
          doSPR = if swapType == "spr" || swapType == "nni" || onlySPR then True
                  else False
          subGraphEdgeVertDataTripleList = if doSPR then [((-1, -1), vertData $ snd prunedGraphRootNode, ([],[]))]
                                           else if hardwiredSPR then take 2 $ getPrunedEdgeData (graphType inGS) doIA inGraph prunedGraphRootNode edgesInPrunedSubGraph charInfoVV
                                           else getPrunedEdgeData (graphType inGS) doIA inGraph prunedGraphRootNode edgesInPrunedSubGraph charInfoVV


          -- this is SPR/NNI
          (delta, _, (tbrEdgesAdd, tbrEdgesDelete)) = getSubGraphDelta targetEdgeData doIA charInfoVV (head subGraphEdgeVertDataTripleList)

          -- this if TBR --need to recurse through all reootings of pruned tree
          (newTBRList, newAnnealVals) = getSubGraphDeltaTBR inGS inData targetEdgeData [edge0, edge1] (eNode, vNode) doIA inGraph splitCost curBestCost charInfoVV inSimAnnealVals subGraphEdgeVertDataTripleList

      in
      -- trace ("ASGR: " ++ (show (delta, splitCost, delta + splitCost, curBestCost))) (
      -- do not redo origal edge so retun infinite cost and dummy edits
      if doSPR && (eNode == nakedNode) then addSubGraphSteepest inGS inData swapType hardwiredSPR doIA inGraph prunedGraphRootNode splitCost curBestCost nakedNode onlySPR edgesInPrunedSubGraph charInfoVV inSimAnnealVals (tail targetEdgeList)

      --TBR case
      else if length subGraphEdgeVertDataTripleList > 1 then
         if null newTBRList then addSubGraphSteepest inGS inData swapType hardwiredSPR doIA inGraph prunedGraphRootNode splitCost curBestCost nakedNode onlySPR edgesInPrunedSubGraph charInfoVV newAnnealVals (tail targetEdgeList)

         -- must be better if non-empty list
         else (newTBRList, newAnnealVals)

      -- SPR case
      -- regular non-SA case
      else if inSimAnnealVals == Nothing then
          -- better heursitic cost
          -- reoptimize to check  cost
          if (delta + splitCost <= curBestCost) then
             let splitGraphSimple = GO.convertDecoratedToSimpleGraph inGraph
                 swapSimpleGraph = applyGraphEdits splitGraphSimple (delta + splitCost, [edge0, edge1] ++ tbrEdgesAdd, (eNode, vNode) : tbrEdgesDelete)
                 reoptimizedCandidateGraph = if (graphType inGS == Tree) then T.multiTraverseFullyLabelGraph inGS inData False False Nothing swapSimpleGraph
                                             else
                                                if (not . LG.cyclic) swapSimpleGraph && (not . GO.parentInChain) swapSimpleGraph then  T.multiTraverseFullyLabelGraph inGS inData False False Nothing swapSimpleGraph
                                                else emptyPhylogeneticGraph
             in
             if (snd6 reoptimizedCandidateGraph < curBestCost) then
                {-
                if (graphType inGS /= Tree) then
                    -- check for cycles and anc/desc time violations in Networks
                    if (not . LG.cyclic) swapSimpleGraph && (not . GO.parentInChainthen swapSimpleGraph) then (filter ((== False) . (GO.hasNetNodeAncestorViolation . thd6)) [reoptimizedCandidateGraph], Nothing)
                    else ([], Nothing)
                else ([reoptimizedCandidateGraph], Nothing)
                -}
                ([reoptimizedCandidateGraph], Nothing)
             else addSubGraphSteepest inGS inData swapType hardwiredSPR doIA inGraph prunedGraphRootNode splitCost curBestCost nakedNode onlySPR edgesInPrunedSubGraph charInfoVV Nothing (tail targetEdgeList)

          -- not better heuristic cost
          else addSubGraphSteepest inGS inData swapType hardwiredSPR doIA inGraph prunedGraphRootNode splitCost curBestCost nakedNode onlySPR edgesInPrunedSubGraph charInfoVV Nothing (tail targetEdgeList)


      -- simulated annealing/Drift case
      else
        --simulated Annealing check if at end of steps then return
        -- original sim anneal  values since is SPR case and not updated yet
        let splitGraphSimple = GO.convertDecoratedToSimpleGraph inGraph
            swapSimpleGraph = applyGraphEdits splitGraphSimple (delta + splitCost, [edge0, edge1] ++ tbrEdgesAdd, (eNode, vNode) : tbrEdgesDelete)
            reoptimizedCandidateGraph = if (graphType inGS == Tree) then T.multiTraverseFullyLabelGraph inGS inData False False Nothing swapSimpleGraph
                                        else
                                            if (not . LG.cyclic) swapSimpleGraph  && (not . GO.parentInChain) swapSimpleGraph then T.multiTraverseFullyLabelGraph inGS inData False False Nothing swapSimpleGraph
                                            else emptyPhylogeneticGraph

        in
        ([reoptimizedCandidateGraph], newAnnealVals)


      -- )

-- | getSubGraphDeltaTBR calculated cost of adding a subgraph into and edge
-- for SPR use the preliminary of subGraph to final of e and v nodes
-- can use median fruntions for postorder if set final-> prelim or e and f
getSubGraphDeltaTBR :: GlobalSettings
                    -> ProcessedData
                    -> VertexBlockData
                    -> [LG.LEdge Double]
                    -> LG.Edge
                    -> Bool
                    -> DecoratedGraph
                    -> VertexCost
                    -> VertexCost
                    -> V.Vector (V.Vector CharInfo)
                    -> Maybe SAParams
                    -> [(LG.Edge, VertexBlockData, ([LG.LEdge Double],[LG.Edge]))]
                    -> ([PhylogeneticGraph], Maybe SAParams)
getSubGraphDeltaTBR inGS inData evEdgeData edgeToAddInList edgeToDeleteIn doIA inGraph splitCost curBestCost charInfoVV inSimAnnealVals subGraphEdgeVertDataTripleList =
   -- trace ("SGD") (
   -- found nothing better
   if null subGraphEdgeVertDataTripleList then ([], inSimAnnealVals)
   else
      let (_, subGraphVertData, (tbrEdgesAdd, tbrEdgesDelete)) = head subGraphEdgeVertDataTripleList
         -- Use edge union data for delta to edge data
          -- costMethod = if doIA then ImpliedAlignment
          --             else DirectOptimization

          subGraphEdgeUnionCost = if (not doIA) then V.sum $ fmap V.sum $ fmap (fmap snd) $ POS.createVertexDataOverBlocks subGraphVertData evEdgeData charInfoVV []
                                  else V.sum $ fmap V.sum $ fmap (fmap snd) $ POS.createVertexDataOverBlocksStaticIA subGraphVertData evEdgeData charInfoVV []

      in

      -- regular TBR case
      if inSimAnnealVals == Nothing then
          if subGraphEdgeUnionCost + splitCost < curBestCost then
             -- reoptimize and check
             let splitGraphSimple = GO.convertDecoratedToSimpleGraph inGraph
                 swapSimpleGraph = applyGraphEdits splitGraphSimple (subGraphEdgeUnionCost + splitCost, edgeToAddInList ++ tbrEdgesAdd, edgeToDeleteIn : tbrEdgesDelete)
                 reoptimizedCandidateGraph = T.multiTraverseFullyLabelGraph inGS inData False False Nothing swapSimpleGraph
             in
             if snd6 reoptimizedCandidateGraph < curBestCost then ([reoptimizedCandidateGraph], Nothing)
             else getSubGraphDeltaTBR inGS inData evEdgeData edgeToAddInList edgeToDeleteIn doIA inGraph splitCost curBestCost charInfoVV Nothing (tail subGraphEdgeVertDataTripleList)

          -- move on
         else
             getSubGraphDeltaTBR inGS inData evEdgeData edgeToAddInList edgeToDeleteIn doIA inGraph splitCost curBestCost charInfoVV Nothing (tail subGraphEdgeVertDataTripleList)

      -- Simulated Annealing/Drifting case
      else
        let -- update annealing values for next round
            -- abstract stopping criterion to continue
            numDone = if (method $ fromJust inSimAnnealVals) == SimAnneal then currentStep $ fromJust inSimAnnealVals
                      else driftChanges $ fromJust inSimAnnealVals
            numMax  = if (method $ fromJust inSimAnnealVals) == SimAnneal then numberSteps $ fromJust inSimAnnealVals
                      else driftMaxChanges $ fromJust inSimAnnealVals

            -- get acceptance based on heuristic costs
            (acceptGraph, nextAnnealVals) = simAnnealAccept inSimAnnealVals curBestCost (subGraphEdgeUnionCost + splitCost)

            -- optimize graph
            splitGraphSimple = GO.convertDecoratedToSimpleGraph inGraph
            swapSimpleGraph = applyGraphEdits splitGraphSimple (subGraphEdgeUnionCost + splitCost, edgeToAddInList ++ tbrEdgesAdd, edgeToDeleteIn : tbrEdgesDelete)
            reoptimizedCandidateGraph = T.multiTraverseFullyLabelGraph inGS inData False False Nothing swapSimpleGraph

        in
        if (numDone < numMax)  then
            -- accept based on heurisrtic cost
            if acceptGraph then
                ([reoptimizedCandidateGraph], nextAnnealVals)
            else
                getSubGraphDeltaTBR inGS inData evEdgeData edgeToAddInList edgeToDeleteIn doIA inGraph splitCost curBestCost charInfoVV nextAnnealVals (tail subGraphEdgeVertDataTripleList)
        else
            ([reoptimizedCandidateGraph], nextAnnealVals)



-- | addSubTree "adds" a subtree back into an edge calculating the cost of the graph via the delta of the add and costs of the two components
-- does NOT reoptimize candidate trees--happens after return since tlooking at "all"
addSubGraph :: GlobalSettings
            -> String
            -> Bool
            -> Bool
            -> DecoratedGraph
            -> LG.LNode VertexInfo
            -> VertexCost
            -> LG.Node
            -> Bool
            -> [LG.LEdge EdgeInfo]
            -> V.Vector (V.Vector CharInfo)
            -> LG.LEdge EdgeInfo
            -> [(VertexCost, [LG.LEdge Double], [LG.Edge])]
addSubGraph inGS swapType hardWiredSPR doIA inGraph prunedGraphRootNode splitCost parentPrunedRoot onlySPR edgesInPrunedSubGraph charInfoVV targetEdge@(eNode, vNode, _) =
   let -- existingEdgeCost = minLength targetlabel
       edge0 = (parentPrunedRoot, vNode, 0.0)
       edge1 = (eNode, parentPrunedRoot, 0.0)
       targetEdgeData = M.makeEdgeData doIA inGraph charInfoVV targetEdge
       -- edge2 = (nakedNode, prunedGraphRootIndex, 0.0)
       --newNode = (nakedNode, TL.pack ("HTU" ++ (show nakedNode)))
       -- if subtree fewer than 3 leaves then can only do an SPR rearragement--no rerro0ts
       -- prunedGraphRootNode = (prunedGraphRootIndex, fromJust $ LG.lab inGraph prunedGraphRootIndex)
       -- (nodesAfterPrunedRoot, edgesInPrunedSubGraph) = LG.nodesAndEdgesAfter inGraph [prunedGraphRootNode]
       -- onlySPR = length nodesAfterPrunedRoot < 3
       doSPR = if swapType == "spr" || swapType == "nni" || onlySPR then True
               else False
       subGraphEdgeVertDataTripleList = if doSPR then [((-1, -1), vertData $ snd prunedGraphRootNode, ([],[]))]
                                        else if hardWiredSPR then take 2 $ getPrunedEdgeData (graphType inGS) doIA inGraph prunedGraphRootNode edgesInPrunedSubGraph charInfoVV
                                        else getPrunedEdgeData (graphType inGS) doIA inGraph prunedGraphRootNode edgesInPrunedSubGraph charInfoVV

       -- get deltas and edges for TBR rerooting of pruned subgraph
       deltaEdgeTripleList = fmap (getSubGraphDelta targetEdgeData doIA charInfoVV) subGraphEdgeVertDataTripleList `using` PU.myParListChunkRDS

       delta = minimum $ fmap fst3 deltaEdgeTripleList

       minDeltaEditList = fmap thd3 $ filter ((== delta) . fst3) deltaEdgeTripleList


       -- get TBR edits if rerooting took place
       tbrEdgeEditList = if length subGraphEdgeVertDataTripleList == 1 then []
                         else minDeltaEditList

       (tbrEdgesToAddList, tbrEdgesToDeleteList) = unzip tbrEdgeEditList

       deltaCostList = replicate (length tbrEdgeEditList) (delta + splitCost)
       edgesToAddList = zipWith (++) (replicate (length tbrEdgeEditList) [edge0, edge1]) tbrEdgesToAddList
       edgesToDeleteList = zipWith (:) (replicate (length tbrEdgeEditList) (eNode, vNode)) tbrEdgesToDeleteList

   in

   -- do not redo origal edge so retun infinite cost and dummy edits
   --trace ("ASG " ++ (show (length subGraphEdgeVertDataTripleList, delta, length tbrEdgeEditList, length deltaCostList, length edgesToAddList, length edgesToDeleteList))) (
   {-
   if (eNode == parentPrunedRoot) then
      -- trace ("ASG: break edge")
      [(infinity, [], [])]
   else
    -}
      -- SPR or single reroot TBR case
      if null tbrEdgeEditList then
        -- trace ("all single roots")
        [(delta + splitCost, [edge0, edge1], [(eNode, vNode)])]

      -- TBR reroots
      else
        -- trace ("all multiple roots")
        zip3 deltaCostList edgesToAddList edgesToDeleteList
      -- )


-- | getTBREdgeEdits takes and edge and returns the list of edit to pruned subgraph
-- as a pair of edges to add and those to delete
-- since reroot edge is directed (e,v), edges away from v will have correct
-- orientation. Edges between 'e' and the root will have to be flipped
-- original root edges and reroort edge are deleted and new root and edge spanning orginal root created
-- returns ([add], [delete])
getTBREdgeEdits :: DecoratedGraph -> LG.LNode VertexInfo -> LG.Edge -> ([LG.LEdge Double],[LG.Edge])
getTBREdgeEdits inGraph prunedGraphRootNode  rerootEdge =
   --trace ("Gettiung TBR Edits for " ++ (show rerootEdge)) (
   let prunedGraphRootIndex = fst prunedGraphRootNode
       originalRootEdgeNodes = LG.descendants inGraph prunedGraphRootIndex
       originalRootEdges = LG.out inGraph prunedGraphRootIndex

       -- get path from new root edge fst vertex to orginal root and flip those edges
       closerToPrunedRootEdgeNode = (fst rerootEdge, fromJust $ LG.lab inGraph $ fst rerootEdge)
       (nodesInPath, edgesinPath) = LG.postOrderPathToNode inGraph closerToPrunedRootEdgeNode prunedGraphRootNode

       -- don't want original root edges to be flipped since deleted
       edgesToFlip = edgesinPath L.\\ originalRootEdges
       flippedEdges = fmap GO.convertToSimpleEdge $ fmap LG.flipLEdge edgesToFlip

       -- new edges on new root position and spanning old root
       -- ad in closer vertex to root to make sure direction of edge is correct
       newEdgeOnOldRoot = if (snd3 $ head originalRootEdges) `elem` ((fst rerootEdge) : (fmap fst nodesInPath)) then (snd3 $ head originalRootEdges, snd3 $ last originalRootEdges, 0.0)
                          else (snd3 $ last originalRootEdges, snd3 $ head originalRootEdges, 0.0)
       newRootEdges = [(prunedGraphRootIndex, fst rerootEdge, 0.0 ),(prunedGraphRootIndex, snd rerootEdge, 0.0)]


   in
   -- original root edge so no change
   if (fst rerootEdge) `elem` originalRootEdgeNodes &&  (snd rerootEdge) `elem` originalRootEdgeNodes then ([],[])

   -- rerooted
   else
      -- delete orignal root edges and rerootEdge
      -- add new root edges
      -- and new edge on old root--but need orientation
      -- flip edges from new root to old (delete and add list)
      --trace ("\n\nIn Graph:\n"++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inGraph) ++ "\nTBR Edits: " ++ (show (rerootEdge, prunedGraphRootIndex, fmap LG.toEdge flippedEdges))
      --   ++ "\nEdges to add: " ++ (show $ fmap LG.toEdge $ newEdgeOnOldRoot : (flippedEdges ++ newRootEdges)) ++ "\nEdges to delete: " ++ (show $ rerootEdge : (fmap LG.toEdge (edgesToFlip ++ originalRootEdges))))
      (newEdgeOnOldRoot : (flippedEdges ++ newRootEdges), rerootEdge : (fmap LG.toEdge (edgesToFlip ++ originalRootEdges)))
      -- )


-- | getPrunedEdgeData takes fully optimized pruned data and returns edge list, edge data for TBR additions,
-- and graph edits for that edge reroot
getPrunedEdgeData :: GraphType
                  -> Bool
                  -> DecoratedGraph
                  -> LG.LNode VertexInfo
                  -> [LG.LEdge EdgeInfo]
                  -> V.Vector (V.Vector CharInfo)
                  -> [(LG.Edge, VertexBlockData, ([LG.LEdge Double],[LG.Edge]))]
getPrunedEdgeData inGraphType doIA inGraph prunedGraphRootNode edgesInPrunedSubGraph charInfoVV =
   if LG.isEmpty inGraph then error "Empty graph in getPrunedEdgeData"
   else
      -- trace ("PED") (
      let -- (_, edgeAfterList) = LG.nodesAndEdgesAfter inGraph [prunedGraphRootNode]

          -- this so not rerooted on network edge--screws things up
          nonNetWorkEdgeList = if inGraphType /= Tree then filter ((== False) . (LG.isNetworkLabEdge inGraph)) edgesInPrunedSubGraph
                               else edgesInPrunedSubGraph

          -- oringal pruned root
          prunedRootEdges = LG.out inGraph $ fst prunedGraphRootNode

          -- new virtual edge of oringal root 2 edges
          virtualRootEdge = (snd3 $ head prunedRootEdges, snd3 $ last prunedRootEdges, thd3 $ head prunedRootEdges)

          -- edges availabel for testing
          edgeAfterList' = virtualRootEdge : (nonNetWorkEdgeList L.\\ prunedRootEdges)

          -- could be parallelized
          edgeDataList = fmap (M.makeEdgeData doIA inGraph charInfoVV) edgeAfterList' `using` PU.myParListChunkRDS

          -- get potential TBR edits--here so not recalculated multiple times for each prune
          edgeEdits = fmap (getTBREdgeEdits inGraph prunedGraphRootNode) (fmap LG.toEdge edgeAfterList') `using` PU.myParListChunkRDS
      in
      if length prunedRootEdges /= 2 then error ("Incorrect number of out edges (should be 2) in root of pruned graph: " ++ (show $ length prunedRootEdges))
      else
         -- trace ("PED Lengths: " ++ (show  $ (length edgeAfterList', length edgeDataList, length edgeEdits)))
         zip3 (fmap LG.toEdge edgeAfterList') edgeDataList edgeEdits
    -- )


-- | getSubGraphDelta calculated cost of adding a subgraph into and edge
-- for SPR use the preliminary of subGraph to final of e and v nodes
-- can use median fruntions for postorder if set final-> prelim or e and f
getSubGraphDelta :: VertexBlockData
                 -> Bool
                 -> V.Vector (V.Vector CharInfo)
                 -> (LG.Edge, VertexBlockData, ([LG.LEdge Double],[LG.Edge]))
                 -> (VertexCost, LG.Edge, ([LG.LEdge Double],[LG.Edge]))
getSubGraphDelta evEdgeData doIA charInfoVV (edgeToJoin, subGraphVertData, edgeSubGraphEdits) =
   let -- existingEdgeCost = minLength edgeToJoin
       --eNodeVertData = vertData $ fromJust $ LG.lab inGraph (fst edgeToJoin)
       --vNodeVertData = vertData $ fromJust $ LG.lab inGraph (snd edgeToJoin)
       -- subGraphVertData = snd subGraphEdgeVertDataPair

       -- create edge union 'character' blockData
       -- based on final assignments but set to preliminary
       -- need to filter gaps if DO, not itIA
       --edgeUnionVertData = M.makeEdgeData doIA inGraph charInfoVV evEdge
       -- edgeUnionVertData = M.createEdgeUnionOverBlocks (not doIA) eNodeVertData vNodeVertData charInfoVV []

       {--
       Use edge union data for delta to edge data
       costMethod = if doIA then ImpliedAlignment
                    else DirectOptimization
       -}

       subGraphEdgeUnionCost = if (not doIA) then V.sum $ fmap V.sum $ fmap (fmap snd) $ POS.createVertexDataOverBlocks subGraphVertData evEdgeData charInfoVV []
                              else V.sum $ fmap V.sum $ fmap (fmap snd) $ POS.createVertexDataOverBlocksStaticIA subGraphVertData evEdgeData charInfoVV []

       -- subGraphEdgeUnionCost = sum $ fmap fst $ V.zipWith3 (PRE.getBlockCostPairsFinal costMethod) subGraphVertData edgeUnionVertData charInfoVV


       {-
       This estimate very close to subGraphEdgeUnionCost.  Seems to always underestmate where as subGraphEdgeUnionCost can sometimes overestimate
       but 3x as much work if use costEV, 2x if use existingEdgeCost

       dummyE = M.createEdgeUnionOverBlocks (not doIA) eNodeVertData eNodeVertData charInfoVV []
       dummyV = M.createEdgeUnionOverBlocks (not doIA) vNodeVertData vNodeVertData charInfoVV []
       dummySGV = M.createEdgeUnionOverBlocks (not doIA) (PRE.setFinalToPreliminaryStates subGraphVertData) (PRE.setFinalToPreliminaryStates subGraphVertData) charInfoVV []

       costNewE = V.sum $ fmap V.sum $ fmap (fmap snd) $ POS.createVertexDataOverBlocks dummyE subGraphVertData charInfoVV []
       costNewV = V.sum $ fmap V.sum $ fmap (fmap snd) $ POS.createVertexDataOverBlocks dummyV subGraphVertData charInfoVV []
       costEV = V.sum $ fmap V.sum $ fmap (fmap snd) $ POS.createVertexDataOverBlocks dummyE dummyV charInfoVV []

       subGraphEdgeUnionCost' = (costNewE + costNewV - costEV) / 2.0
       -}



   in
   -- remove this check when things are working
   {-
   if null subGraphVertData || null evEdgeData || (subGraphEdgeUnionCost == 0.0) then
      trace ("SGD null or 0 :" ++ (show (edgeToJoin, length subGraphVertData, length evEdgeData, subGraphEdgeUnionCost)) )
      --  ++ "\nInGraph:\n"
      --  ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inGraph) ++ "\n" ++ (show inGraph))
      (subGraphEdgeUnionCost, edgeToJoin, edgeSubGraphEdits)
    --trace ("GSD:" ++ (show ((costNewE, costNewV, costEV))) ++ " -> " ++ (show subGraphEdgeUnionCost') ++  " v " ++ (show subGraphEdgeUnionCost))
    -- trace ("Delta: " ++ (show subGraphEdgeUnionCost))
    --min raphEdgeUnionCost  subGraphEdgeUnionCost'
   else
   -}
      -- trace ("SGD no 0:" ++ (show (subGraphEdgeUnionCost,edgeToJoin)))
    (subGraphEdgeUnionCost, edgeToJoin, edgeSubGraphEdits)


-- | reoptimizeSplitGraphFromVertex fully labels the component graph that is connected to the specified vertex
-- retuning that graph with 2 optimized components and their cost
-- both components goo through multi-traversal optimizations
-- doIA option to only do IA optimization as opposed to full thing--should be enormously faster--but yet more approximate
-- creates final for both base graph and priunned component due to rerooting non-concordance of preorder and post order assignments
-- terminology bse graph is the component with the original root, pruned that which has been removed form the original
-- graph to be readded to edge set
-- The function
      -- 1) optimizes two components seprately from their "root"
      -- 2) takes nodes and edges for each and cretes new graph
      -- 3) returns graph and summed cost of two components
      -- 4) adds in root and netPenalty factor estimates since net penalty can only be calculated on full graph
            -- part of this is turning off net penalty cost when optimizing base and pruned graph components
-- if doIA is TRUE then call function that onl;y optimizes the IA assignments on the "original graph" after split.
-- this keeps teh IA chracters in sync across the two graphs
reoptimizeSplitGraphFromVertex :: GlobalSettings
                          -> ProcessedData
                          -> Bool
                          -> VertexCost
                          -> DecoratedGraph
                          -> Int
                          -> Int
                          -> (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertex inGS inData doIA netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex =
   if doIA then
      -- only reoptimize the IA states for dynamic characters
      reoptimizeSplitGraphFromVertexIA inGS inData netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex
   else
      -- perform full optimizations of nodes
      -- these required for full optimization
      let nonExactCharacters = U.getNumberSequenceCharacters (thd3 inData)
          origGraph = inSplitGraph -- thd6 origPhyloGraph
          leafGraph = LG.extractLeafGraph origGraph
          calcBranchLengths = False

          -- create simple graph version of split for post order pass
          splitGraphSimple = GO.convertDecoratedToSimpleGraph inSplitGraph


          -- create optimized base graph
          -- False for staticIA
          (postOrderBaseGraph, localRootCost, _) = T.generalizedGraphPostOrderTraversal (inGS {graphFactor = NoNetworkPenalty}) nonExactCharacters inData leafGraph False (Just startVertex) splitGraphSimple


          fullBaseGraph = PRE.preOrderTreeTraversal (inGS {graphFactor = NoNetworkPenalty}) (finalAssignment inGS) False calcBranchLengths (nonExactCharacters > 0) startVertex True postOrderBaseGraph

          -- create fully optimized pruned graph.  Post order tehn preorder

          -- get root node of pruned graph--parent since that is the full pruned piece (keeping that node for addition to base graph and edge creation)
          startPrunedNode = (prunedSubGraphRootVertex, fromJust $ LG.lab origGraph prunedSubGraphRootVertex)
          startPrunedParentNode = head $ LG.labParents origGraph prunedSubGraphRootVertex
          startPrunedParentEdge = (fst startPrunedParentNode, prunedSubGraphRootVertex, dummyEdge)


          -- False for staticIA
          (postOrderPrunedGraph, _, _) = T.generalizedGraphPostOrderTraversal (inGS {graphFactor = NoNetworkPenalty}) nonExactCharacters inData leafGraph False (Just prunedSubGraphRootVertex) splitGraphSimple


          -- False for staticIA
          fullPrunedGraph = PRE.preOrderTreeTraversal (inGS {graphFactor = NoNetworkPenalty}) (finalAssignment inGS) False calcBranchLengths (nonExactCharacters > 0) prunedSubGraphRootVertex True postOrderPrunedGraph

          -- get root node of base graph
          startBaseNode = (startVertex, fromJust $ LG.lab (thd6 fullBaseGraph) startVertex)



          -- get nodes and edges in base and pruned graph (both PhylogeneticGrapgs so thd6)
          (baseGraphNonRootNodes, baseGraphEdges) = LG.nodesAndEdgesAfter (thd6 fullBaseGraph) [startBaseNode]

          (prunedGraphNonRootNodes, prunedGraphEdges) = if LG.isLeaf origGraph prunedSubGraphRootVertex then ([],[])
                                                        else LG.nodesAndEdgesAfter (thd6 fullPrunedGraph) [startPrunedNode]

          -- make fully optimized graph from base and split components
          fullSplitGraph = LG.mkGraph ([startBaseNode, startPrunedNode, startPrunedParentNode] ++  baseGraphNonRootNodes ++ prunedGraphNonRootNodes) (startPrunedParentEdge : (baseGraphEdges ++ prunedGraphEdges))

          -- cost of split graph to be later combined with re-addition delta for heuristic graph cost
          prunedCost = if LG.isLeaf origGraph prunedSubGraphRootVertex then 0
                       else snd6 fullPrunedGraph
          splitGraphCost = ((1.0 + netPenaltyFactor) * ((snd6 fullBaseGraph) + prunedCost)) + localRootCost

      in

      {-
      trace ("Orig graph cost " ++ (show $ subGraphCost $ fromJust $ LG.lab origGraph startVertex) ++ " Base graph cost " ++ (show $ snd6 fullBaseGraph) ++ " pruned subgraph cost " ++ (show prunedCost) ++ " at node " ++ (show prunedSubGraphRootVertex) ++ " parent " ++ (show $ fst startPrunedParentNode)

         ++ "\nBaseGraphNodes\n" ++ (show $ L.sort  $ fmap fst baseGraphNonRootNodes) ++ "\nPruned nodes from root: " ++ "\n" ++ (show $ fmap fst $ startPrunedNode : prunedGraphNonRootNodes)
         ++ "\nSplit Graph\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph fullSplitGraph)
         ++ "\nOrig graph:\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph origGraph))
      -}
      --trace ("reoptimizeSplitGraphFromVertex: " ++ (show splitGraphCost))
      (fullSplitGraph, splitGraphCost)

-- | reoptimizeSplitGraphFromVertexTuple wrapper for reoptimizeSplitGraphFromVertex with last 3 args as tuple
reoptimizeSplitGraphFromVertexTuple :: GlobalSettings
                          -> ProcessedData
                          -> Bool
                          -> VertexCost
                          -> (DecoratedGraph, Int , Int)
                          -> (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertexTuple inGS inData doIA netPenaltyFactor (inSplitGraph, startVertex, prunedSubGraphRootVertex) =
   reoptimizeSplitGraphFromVertex inGS inData doIA netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex


-- | reoptimizeSplitGraphFromVertexIA performs operations of reoptimizeSplitGraphFromVertex for static charcaters
-- but dynamic characters--only update IA assignments and initialized from origPhylo graph (at leaves) to keep IA characters in sync
-- since all "static" only need single traversal post order pass
reoptimizeSplitGraphFromVertexIA :: GlobalSettings
                          -> ProcessedData
                          -> VertexCost
                          -> DecoratedGraph
                          -> Int
                          -> Int
                          -> (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertexIA inGS inData netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex =
   --if graphType inGS /= Tree then error "Networks not yet implemented in reoptimizeSplitGraphFromVertexIA"
   --else
      let   nonExactCharacters = U.getNumberSequenceCharacters (thd3 inData)
            origGraph = inSplitGraph -- thd6 origPhyloGraph

            -- create leaf graphs--but copy IA final to prelim
            leafGraph = GO.copyIAFinalToPrelim $ LG.extractLeafGraph origGraph
            calcBranchLengths = False

            -- create simple graph version of split for post order pass
            splitGraphSimple = GO.convertDecoratedToSimpleGraph inSplitGraph

            --Create base graph
            -- create postorder assignment--but only from single traversal
            -- True flag fior staticIA
            postOrderBaseGraph = T.postOrderTreeTraversal (inGS {graphFactor = NoNetworkPenalty}) inData leafGraph True (Just startVertex) splitGraphSimple
            baseGraphCost = snd6 postOrderBaseGraph

            -- True flag fior staticIA
            fullBaseGraph = PRE.preOrderTreeTraversal (inGS {graphFactor = NoNetworkPenalty}) (finalAssignment inGS) True calcBranchLengths (nonExactCharacters > 0) startVertex True postOrderBaseGraph

            localRootCost = if (rootCost inGS) == NoRootCost then 0.0
                              else if (rootCost inGS) == Wheeler2015Root then T.getW15RootCost inData postOrderBaseGraph
                              else error ("Root cost type " ++ (show $ rootCost inGS) ++ " is not yet implemented")

            -- get root node of base graph
            startBaseNode = (startVertex, fromJust $ LG.lab (thd6 fullBaseGraph) startVertex)

            --Create pruned graph
            -- get root node of pruned graph--parent since that is the full pruned piece (keeping that node for addition to base graph and edge creation)
            startPrunedNode = GO.makeIAPrelimFromFinal (prunedSubGraphRootVertex, fromJust $ LG.lab origGraph prunedSubGraphRootVertex)
            startPrunedParentNode =  head $ LG.labParents origGraph prunedSubGraphRootVertex
            startPrunedParentEdge = (fst startPrunedParentNode, prunedSubGraphRootVertex, dummyEdge)


            -- True flag fior staticIA
            postOrderPrunedGraph =  T.postOrderTreeTraversal (inGS {graphFactor = NoNetworkPenalty}) inData leafGraph True (Just prunedSubGraphRootVertex) splitGraphSimple
            prunedGraphCost = snd6 postOrderPrunedGraph

            -- True flag fior staticIA
            fullPrunedGraph = PRE.preOrderTreeTraversal (inGS {graphFactor = NoNetworkPenalty}) (finalAssignment inGS) True calcBranchLengths (nonExactCharacters > 0) prunedSubGraphRootVertex True postOrderPrunedGraph

            -- get nodes and edges in base and pruned graph (both PhylogeneticGrapgs so thd6)
            (baseGraphNonRootNodes, baseGraphEdges) = LG.nodesAndEdgesAfter (thd6 fullBaseGraph) [startBaseNode]

            (prunedGraphNonRootNodes, prunedGraphEdges) = if LG.isLeaf origGraph prunedSubGraphRootVertex then ([],[])
                                                          else LG.nodesAndEdgesAfter (thd6 fullPrunedGraph) [startPrunedNode]

            -- make fully optimized graph from base and split components
            fullSplitGraph = LG.mkGraph ([startBaseNode, startPrunedNode, startPrunedParentNode] ++  baseGraphNonRootNodes ++ prunedGraphNonRootNodes) (startPrunedParentEdge : (baseGraphEdges ++ prunedGraphEdges))

            splitGraphCost = ((1.0 + netPenaltyFactor) * (baseGraphCost + prunedGraphCost)) + localRootCost

      in
      -- remove when working
      -- trace ("ROGFVIA split costs:" ++ (show (baseGraphCost, prunedGraphCost, localRootCost)) ++ " -> " ++ (show splitGraphCost)) (
      if splitGraphCost == 0 then
         error ("Split costs:" ++ (show (baseGraphCost, prunedGraphCost, localRootCost)) ++ " -> " ++ (show splitGraphCost)
            ++ " Split graph simple:\n" ++ (LG.prettify splitGraphSimple)
            ++ "\nFull:\n" ++ (show inSplitGraph)
            ++ "\nOriginal Graph:\n" ++ (show origGraph))
      else (fullSplitGraph, splitGraphCost)
      -- )

{-
-- | applyGraphEdits' takes a  graphs and list of nodes and edges to add and delete and creates new graph
applyGraphEdits' :: (Show a, Show b) => LG.Gr a b -> (VertexCost, [LG.LEdge b], [LG.Edge]) ->  LG.Gr a b
applyGraphEdits' inGraph (_, edgesToAdd, edgesToDelete) =
   let editedGraph = LG.insEdges edgesToAdd $ LG.delEdges edgesToDelete inGraph
   in
   -- trace ("AGE: " ++ (show editStuff) ++ "\nIn graph:\n" ++ (LG.prettify inGraph) ++ "New Graph:\n" ++ (LG.prettify editedGraph))
   editedGraph
-}

-- | applyGraphEdits takes a  graphs and list of nodes and edges to add and delete and creates new graph
applyGraphEdits :: (Show a, Show b) => LG.Gr a b -> (VertexCost, [LG.LEdge b], [LG.Edge]) ->  LG.Gr a b
applyGraphEdits inGraph (_, edgesToAdd, edgesToDelete) = LG.insertDeleteEdges inGraph (edgesToAdd, edgesToDelete)

