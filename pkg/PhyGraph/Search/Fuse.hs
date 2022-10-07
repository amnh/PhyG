{- |
Module      :  Fuse.hs
Description :  Module specifying graph fusing recombination functions
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

module Search.Fuse  ( fuseAllGraphs
                    ) where

import           Control.Parallel.Strategies
import           Data.Bits
import qualified Data.BitVector.LittleEndian          as BV
import qualified Data.List                            as L
import qualified Data.Map                             as MAP
import           Data.Maybe
import qualified Data.Text.Lazy                       as TL
import qualified Data.Vector                          as V
import           Debug.Trace
import           GeneralUtilities
import qualified GraphOptimization.Traversals         as T
import qualified Graphs.GraphOperations               as GO
import qualified ParallelUtilities                    as PU
import qualified Search.Swap                          as S
import           Types.Types
import qualified Utilities.LocalGraph                 as LG
import qualified Data.InfList                as IL

-- | fuseAllGraphs takes a list of phylogenetic graphs and performs all pairwise fuses
-- later--could limit by options making random choices for fusing
-- keeps results according to options (best, unique, etc)
-- unique is unique of "best" from individual fusings
-- singleRound short circuits recursive continuation on newly found graphs
fuseAllGraphs :: GlobalSettings
              -> ProcessedData
              -> [Int]
              -> Int
              -> Int
              -> Int
              -> Bool
              -> Bool
              -> Bool
              -> Bool
              -> Bool
              -> Bool
              -> Bool
              -> Bool
              -> Maybe Int
              -> Bool
              -> [PhylogeneticGraph]
              -> ([PhylogeneticGraph], Int)
fuseAllGraphs inGS inData rSeedList keepNum maxMoveEdgeDist counter doNNI doSPR doTBR doSteepest doAll returnBest returnUnique singleRound fusePairs randomPairs inGraphList =
   if null inGraphList then ([], 0)
   else if length inGraphList == 1 then (inGraphList, 0)
   else
      let -- getting values to be passed for graph diagnosis later
         numLeaves = V.length $ fst3 inData
         -- leafGraph = T.makeSimpleLeafGraph inData
         -- leafDecGraph = T.makeLeafGraph inData
         -- leafGraphSoftWired = T.makeLeafGraphSoftWired inData
         -- hasNonExactChars = U.getNumberSequenceCharacters (thd3 inData) > 0
         charInfoVV = six6 $ head inGraphList

         curBest = minimum $ fmap snd6 inGraphList

         curBestGraph = head $ filter ((== curBest) . snd6) inGraphList

         -- get net penalty estimate from optimal graph for delta recombine later
         inGraphNetPenalty = if (graphType inGS == Tree) then 0.0
                             else if (graphFactor inGS) == NoNetworkPenalty then 0.0
                             else if (graphFactor inGS) == Wheeler2015Network then 
                                 if (graphType inGS) == HardWired then 0.0
                                 else T.getW15NetPenalty Nothing curBestGraph
                             else if (graphFactor inGS) == Wheeler2023Network then 
                                 if (graphType inGS) == HardWired then 0.0
                                 else T.getW23NetPenalty Nothing curBestGraph
                             else if (graphFactor inGS) == PMDLGraph then 
                                 let (_, _, _, networkNodeList) = LG.splitVertexList (fst6 curBestGraph)
                                 in
                                 if (graphType inGS) == Tree then fst $ IL.head (graphComplexityList inGS)
                                 else if (graphType inGS) == SoftWired then fst $ (graphComplexityList inGS) IL.!!! (length networkNodeList)
                                 else if (graphType inGS) == HardWired then snd $ (graphComplexityList inGS) IL.!!! (length networkNodeList)
                                 else error ("Graph type " ++ (show $ graphType inGS) ++ " is not yet implemented in fuseAllGraphs")
                             else error ("Network penalty type " ++ (show $ graphFactor inGS) ++ " is not yet implemented")
         inGraphNetPenaltyFactor = inGraphNetPenalty / curBest


         -- get fuse pairs
         graphPairList' = getListPairs inGraphList
         (graphPairList, randString) = if isNothing fusePairs then (graphPairList', "")
                         else if randomPairs then (takeRandom (head rSeedList) (fromJust fusePairs) graphPairList', " randomized")
                         else (takeNth (fromJust fusePairs) graphPairList', "")


         newGraphList = concat (PU.seqParMap rdeepseq (fusePair inGS inData numLeaves charInfoVV inGraphNetPenaltyFactor keepNum maxMoveEdgeDist doNNI doSPR doTBR) graphPairList) -- `using` PU.myParListChunkRDS)

         fuseBest = if not (null newGraphList) then  minimum $ fmap snd6 newGraphList
                    else infinity

      in

      trace ("\tFusing " ++ (show $ length graphPairList) ++ randString ++ " graph pairs") (
      if null newGraphList then (inGraphList, counter + 1)
      else if returnUnique then
         let uniqueList = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (inGraphList ++ newGraphList)
         in
         if fuseBest < curBest then
               trace ("\t->" ++ (show fuseBest)) --  ++ "\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph $ thd6 $ head bestSwapGraphList))
               fuseAllGraphs inGS inData (drop 2 rSeedList) keepNum maxMoveEdgeDist (counter + 1) doNNI doSPR doTBR doSteepest doAll returnBest returnUnique singleRound fusePairs randomPairs uniqueList
         else (uniqueList, counter + 1)

      else -- return best
         -- only do one round of fusing
         if singleRound then (take keepNum $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] (inGraphList ++ newGraphList), counter + 1)

         -- recursive rounds
         else
            let allBestList = take keepNum $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] (inGraphList ++ newGraphList)
            in

            -- found worse
            if fuseBest > curBest then (allBestList, counter + 1)

            -- found better
            else if fuseBest < curBest then
               trace ("\t->" ++ (show fuseBest)) --  ++ "\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph $ thd6 $ head bestSwapGraphList))
               fuseAllGraphs inGS inData (drop 2  rSeedList) keepNum maxMoveEdgeDist (counter + 1) doNNI doSPR doTBR doSteepest doAll returnBest returnUnique singleRound fusePairs randomPairs allBestList

            -- equal cost just return--could keep finding equal
            else (allBestList, counter + 1)
      )


-- | fusePair recombines a single pair of graphs
-- this is done by coopting the split and readd functinos from the Swap.Swap functions and exchanging
-- pruned subgraphs with the same leaf complement (as recorded by the subtree root node bit vector field)
-- spr-like and tbr-like readds can be performed as with options
fusePair :: GlobalSettings
         -> ProcessedData
         -> Int
         -> V.Vector (V.Vector CharInfo)
         -> VertexCost
         -> Int
         -> Int
         -> Bool
         -> Bool
         -> Bool
         -> (PhylogeneticGraph, PhylogeneticGraph)
         -> [PhylogeneticGraph]
fusePair inGS inData numLeaves charInfoVV netPenalty keepNum maxMoveEdgeDist doNNI doSPR doTBR (leftGraph, rightGraph) =
   if (LG.isEmpty $ fst6 leftGraph) || (LG.isEmpty $ fst6 rightGraph) then error "Empty graph in fusePair"
   else if (fst6 leftGraph) == (fst6 rightGraph) then []
   else
      -- split graphs at all bridge edges (all edges for Tree)
      let -- left graph splits
          leftDecoratedGraph = thd6 leftGraph
          (leftRootIndex, _) = head $ LG.getRoots leftDecoratedGraph
          leftBreakEdgeList = if (graphType inGS) == Tree then filter ((/= leftRootIndex) . fst3) $ LG.labEdges leftDecoratedGraph
                              else filter ((/= leftRootIndex) . fst3) $ LG.getEdgeSplitList leftDecoratedGraph
          leftSplitTupleList = PU.seqParMap rdeepseq  (LG.splitGraphOnEdge leftDecoratedGraph) leftBreakEdgeList -- `using` PU.myParListChunkRDS
          (_, _, leftPrunedGraphRootIndexList,  leftOriginalConnectionOfPrunedList) = L.unzip4 leftSplitTupleList
          --leftPrunedGraphRootIndexList = fmap thd4 leftSplitTupleList
          leftPrunedGraphBVList = fmap bvLabel $ fmap fromJust $ fmap (LG.lab leftDecoratedGraph) leftPrunedGraphRootIndexList


          -- right graph splits
          rightDecoratedGraph = thd6 rightGraph
          (rightRootIndex, _) = head $ LG.getRoots rightDecoratedGraph
          rightBreakEdgeList = if (graphType inGS) == Tree then filter ((/= rightRootIndex) . fst3) $ LG.labEdges rightDecoratedGraph
                              else filter ((/= rightRootIndex) . fst3) $ LG.getEdgeSplitList rightDecoratedGraph
          rightSplitTupleList = PU.seqParMap rdeepseq (LG.splitGraphOnEdge rightDecoratedGraph) rightBreakEdgeList -- `using` PU.myParListChunkRDS
          (_, _, rightPrunedGraphRootIndexList,  rightOriginalConnectionOfPrunedList) = L.unzip4 rightSplitTupleList
          -- rightPrunedGraphRootIndexList = fmap thd4 rightSplitTupleList
          rightPrunedGraphBVList = fmap bvLabel $ fmap fromJust $ fmap (LG.lab rightDecoratedGraph) rightPrunedGraphRootIndexList


          -- need to get all pairs of split graphs
          (leftSplitTupleList', rightSplitTupleList') =  unzip $ cartProd leftSplitTupleList rightSplitTupleList
          (leftPrunedGraphBVList', rightPrunedGraphBVList') = unzip $ cartProd leftPrunedGraphBVList rightPrunedGraphBVList
          -- (leftBaseBVList, rightBaseBVList) = unzip $ cartProd leftBaseGraphBVList rightBaseGraphBVList


          -- get compatible split pairs via checking bv of root index of pruned subgraphs
          leftRightMatchList = zipWith (==) leftPrunedGraphBVList' rightPrunedGraphBVList'

          -- only take compatible, non-identical pairs with > 2 terminal--otherwise basically SPR move or nmothing (if identical)
            -- oalso checks that prune and splits don't match between the grap[hs to be recombined]
          recombinablePairList = L.zipWith (getCompatibleNonIdenticalSplits numLeaves) leftRightMatchList leftPrunedGraphBVList'
          (leftValidTupleList, rightValidTupleList, _) = L.unzip3 $ filter ((==True) . thd3) $ zip3 leftSplitTupleList' rightSplitTupleList' recombinablePairList


          -- create new "splitgraphs" by replacing nodes and edges of pruned subgraph in reciprocal graphs
          -- retuns reindexed list of base graph root, pruned component root,  parent of pruned component root, original graph break edge
          (leftBaseRightPrunedSplitGraphList, leftRightGraphRootIndexList, leftRightPrunedParentRootIndexList, leftRightPrunedRootIndexList, leftRightOriginalConnectionOfPrunedList) = L.unzip5 (PU.seqParMap rdeepseq (exchangePrunedGraphs numLeaves) (zip3 leftValidTupleList rightValidTupleList leftOriginalConnectionOfPrunedList)) -- `using` PU.myParListChunkRDS)

          (rightBaseLeftPrunedSplitGraphList, rightLeftGraphRootIndexList, rightLeftPrunedParentRootIndexList, rightLeftPrunedRootIndexList, rightLeftOriginalConnectionOfPrunedList) = L.unzip5 (PU.seqParMap rdeepseq (exchangePrunedGraphs numLeaves) (zip3 rightValidTupleList leftValidTupleList rightOriginalConnectionOfPrunedList)) -- `using` PU.myParListChunkRDS)

          -- reoptimize splitGraphs so ready for readdition--using updated base and prune indices
          -- False for doIA
          leftRightOptimizedSplitGraphCostList = PU.seqParMap rdeepseq (S.reoptimizeSplitGraphFromVertexTuple inGS inData False netPenalty) (zip3 leftBaseRightPrunedSplitGraphList leftRightGraphRootIndexList leftRightPrunedRootIndexList) -- `using` PU.myParListChunkRDS

          rightLeftOptimizedSplitGraphCostList = PU.seqParMap rdeepseq (S.reoptimizeSplitGraphFromVertexTuple inGS inData False netPenalty) (zip3 rightBaseLeftPrunedSplitGraphList rightLeftGraphRootIndexList rightLeftPrunedRootIndexList) -- `using` PU.myParListChunkRDS

          -- Check if base graphs are different as well (nneded to be reoptimized to get base root bv)
          -- otherwise no point in recombination
          {-
          leftBaseGraphBVList = fmap bvLabel $ fmap fromJust $ zipWith LG.lab (fmap fst leftRightOptimizedSplitGraphCostList) leftRightGraphRootIndexList
          rightBaseGraphBVList =fmap bvLabel $ fmap fromJust $ zipWith LG.lab (fmap fst rightLeftOptimizedSplitGraphCostList) rightLeftGraphRootIndexList
          baseGraphDifferentList = zipWith (/=) leftBaseBVList rightBaseBVList
          -}
          baseGraphDifferentList = L.replicate (length leftRightOptimizedSplitGraphCostList) True

          (_, leftRightOptimizedSplitGraphCostList', _, leftRightPrunedRootIndexList', leftRightPrunedParentRootIndexList', leftRightOriginalConnectionOfPrunedList') = L.unzip6 $ filter ((== True) . fst6) $ L.zip6 baseGraphDifferentList leftRightOptimizedSplitGraphCostList leftRightGraphRootIndexList leftRightPrunedRootIndexList leftRightPrunedParentRootIndexList leftRightOriginalConnectionOfPrunedList

          (_, rightLeftOptimizedSplitGraphCostList', _, rightLeftPrunedRootIndexList', rightLeftPrunedParentRootIndexList', rightLeftOriginalConnectionOfPrunedList') = L.unzip6 $ filter ((== True) . fst6) $ L.zip6 baseGraphDifferentList rightLeftOptimizedSplitGraphCostList rightLeftGraphRootIndexList rightLeftPrunedRootIndexList rightLeftPrunedParentRootIndexList rightLeftOriginalConnectionOfPrunedList

          -- re-add pruned component to base component left-right and right-left
          -- need cure best cost
          curBetterCost = min (snd6 leftGraph) (snd6 rightGraph)

          -- get network penalty factors to pass on
          networkCostFactor = min (getNetworkPentaltyFactor inGS (snd6 leftGraph) leftGraph) (getNetworkPentaltyFactor inGS (snd6 rightGraph) rightGraph) 


          -- left and right root indices shold be the same
          leftRightFusedGraphList = recombineComponents inGS inData keepNum maxMoveEdgeDist doNNI doSPR doTBR charInfoVV curBetterCost leftRightOptimizedSplitGraphCostList' leftRightPrunedRootIndexList' leftRightPrunedParentRootIndexList' leftRightOriginalConnectionOfPrunedList' leftRootIndex networkCostFactor
          rightLeftFusedGraphList = recombineComponents inGS inData keepNum maxMoveEdgeDist doNNI doSPR doTBR charInfoVV curBetterCost rightLeftOptimizedSplitGraphCostList' rightLeftPrunedRootIndexList' rightLeftPrunedParentRootIndexList' rightLeftOriginalConnectionOfPrunedList' rightRootIndex networkCostFactor


          -- get "best" fused graphs from leftRight and rightLeft
          bestFusedGraphs = take keepNum $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] (leftRightFusedGraphList ++ rightLeftFusedGraphList)

          -- | get fuse graphs via swap function
      in
      if null leftValidTupleList then []
      else
         {-
         trace ("FP: " ++ (show (length leftValidTupleList, length rightValidTupleList)) ++ " num (Left,Right) " ++ (show (length leftSplitTupleList, length rightSplitTupleList))
            ++ "\nLeftRight splitCost " ++ (show $ fmap snd leftRightOptimizedSplitGraphCostList)
            ++ "\nrightLeft splitCost " ++ (show $ fmap snd rightLeftOptimizedSplitGraphCostList))
         -}
         bestFusedGraphs

-- | recombineComponents takes readdition arguments (swap, steepest etc) and wraps the swap-stype rejoining of components
-- ignores doSteepeast for now--doesn't seem to have meaning in rejoining since not then taking that graph for fusion and shortcircuiting
-- not doing original connection first (originalConnectionOfPrunedComponentList)-since so much work might as well do soem SPR at least
recombineComponents :: GlobalSettings
                    -> ProcessedData
                    -> Int
                    -> Int
                    -> Bool
                    -> Bool
                    -> Bool
                    -> V.Vector (V.Vector CharInfo)
                    -> VertexCost
                    -> [(DecoratedGraph, VertexCost)]
                    -> [Int]
                    -> [Int]
                    -> [Int]
                    -> LG.Node
                    -> VertexCost
                    -> [PhylogeneticGraph]
recombineComponents inGS inData numToKeep inMaxMoveEdgeDist doNNI doSPR doTBR charInfoVV curBestCost splitGraphCostPairList prunedRootIndexList prunedParentRootIndexList _ graphRoot networkCostFactor =
   -- check and see if any reconnecting to do
   --trace ("RecombineComponents " ++ (show $ length splitGraphCostPairList)) (
   if null splitGraphCostPairList then []
   else
      -- top line to cover SPR HarWired bug
      let swapType = if doTBR then "tbr"
                     else if doSPR then "spr"
                     else if doNNI then "nni"
                     else "alternate"

          doIA = False --- since splits not created together, IA won't be consistent between components
          steepest = False -- should look at all better

          -- network costs--using an input value that is minimum of inputs
          netPenaltyFactorList = L.replicate (length splitGraphCostPairList) networkCostFactor

          -- no simulated annealling functionality infuse
          inSimAnnealParams = Nothing

          -- get edges in pruned (to be exchanged) graphs
          edgesInPrunedList = fmap LG.getEdgeListAfter $ zip (fmap fst splitGraphCostPairList) prunedParentRootIndexList

          -- get edges in base (not to be exchanged) graphs
          rejoinEdgesList = fmap (getBaseGraphEdges graphRoot) $ zip (fmap fst splitGraphCostPairList) edgesInPrunedList

          --huge zip to fit arguments into revised join function
          graphDataList = zip9 (fmap fst splitGraphCostPairList) 
                               (fmap GO.convertDecoratedToSimpleGraph $ fmap fst splitGraphCostPairList) 
                               (fmap snd splitGraphCostPairList)  
                               (L.replicate (length splitGraphCostPairList) graphRoot)
                               prunedRootIndexList
                               prunedParentRootIndexList
                               rejoinEdgesList
                               edgesInPrunedList
                               netPenaltyFactorList

          -- do "all additions" -
          recombinedGraphList = concat $ PU.seqParMap rdeepseq (S.rejoinGraphTuple swapType inGS inData numToKeep inMaxMoveEdgeDist steepest curBestCost [] doIA charInfoVV inSimAnnealParams) graphDataList
                                                              

          -- this based on heuristic deltas
          bestFuseCost = if null recombinedGraphList then infinity
                         else minimum $ fmap snd6 recombinedGraphList

      in
      --trace ("Checking in fusing") (
      if null recombinedGraphList then []
      else if bestFuseCost <= curBestCost then
         take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] recombinedGraphList
      else []
      -- )
      -- )

-- | getNetworkPentaltyFactor get scale network penalty for graph
getNetworkPentaltyFactor :: GlobalSettings -> VertexCost -> PhylogeneticGraph -> VertexCost 
getNetworkPentaltyFactor inGS graphCost inGraph =
   if LG.isEmpty $ thd6 inGraph then 0.0
   else
        let inGraphNetPenalty = if (graphType inGS == Tree) || (graphType inGS == HardWired) then 0.0
                                else if (graphFactor inGS) == NoNetworkPenalty then 0.0
                                else if (graphFactor inGS) == Wheeler2015Network then T.getW15NetPenalty Nothing inGraph
                                else if (graphFactor inGS) == Wheeler2023Network then T.getW23NetPenalty Nothing inGraph
                                else error ("Network penalty type " ++ (show $ graphFactor inGS) ++ " is not yet implemented")
        in
        inGraphNetPenalty / graphCost




-- | getBaseGraphEdges gets teh edges in the base graph teh trhe exchanged sub graphs can be rejoined
-- basically all edges except at root and those in the subgraph
getBaseGraphEdges :: (Eq b) => LG.Node -> (LG.Gr a b, [LG.LEdge b]) -> [LG.LEdge b]
getBaseGraphEdges graphRoot (inGraph, edgesInSubGraph) =
   if LG.isEmpty inGraph then []
   else 
      filter ((/= graphRoot) . fst3) $ (LG.labEdges inGraph) L.\\ edgesInSubGraph

-- | getCompatibleNonIdenticalSplits takes the number of leaves, splitGraph of the left graph, the splitGraph if the right graph,
-- the bitVector equality list of pruned roots, the bitvector of the root of the pruned graph on left
-- (this could be either since filter for identity--just to check leaf numbers)
-- checks that the leaf sets of the pruned subgraphs are equal, greater than 1 leaf, fewer thanm nleaves - 2, and non-identical
-- removed identity check fo now--so much time to do that (O(n)) may not be worth it
getCompatibleNonIdenticalSplits :: Int
                                -> Bool
                                -> BV.BitVector
                                -> Bool
getCompatibleNonIdenticalSplits numLeaves leftRightMatch leftPrunedGraphBV =

   if not leftRightMatch then False
   else if popCount leftPrunedGraphBV < 3 then False
   else if popCount leftPrunedGraphBV > (numLeaves - 3) then False
   else True
   {-
      -- check for pruned components non-identical

      let -- (leftNodesInPrunedGraph, _) = LG.nodesAndEdgesAfter (fst4 leftSplitTuple) [((thd4 leftSplitTuple), fromJust $ LG.lab (fst4 leftSplitTuple) (thd4 leftSplitTuple))]
          -- leftPrunedBVNodeList = L.sort $ filter ((> 1) . popCount) $ fmap (bvLabel . snd)  leftNodesInPrunedGraph
          -- (rightNodesInPrunedGraph, _) = LG.nodesAndEdgesAfter (fst4 rightSplitTuple) [((thd4 rightSplitTuple), fromJust $ LG.lab (fst4 rightSplitTuple) (thd4 rightSplitTuple))]
          -- rightPrunedBVNodeList = L.sort $ filter ((> 1) . popCount) $ fmap (bvLabel . snd)  rightNodesInPrunedGraph


          -- (leftNodesInBaseGraph, _) = LG.nodesAndEdgesAfter (fst4 leftSplitTuple) [((snd4 leftSplitTuple), fromJust $ LG.lab (fst4 leftSplitTuple) (snd4 leftSplitTuple))]
          -- leftBaseBVNodeList = L.sort $ filter ((> 1) . popCount) $ fmap (bvLabel . snd)  leftNodesInPrunedGraph
          -- (rightNodesInBaseGraph, _) = LG.nodesAndEdgesAfter (fst4 rightSplitTuple) [((snd4 rightSplitTuple), fromJust $ LG.lab (fst4 rightSplitTuple) (snd4 rightSplitTuple))]
          -- rightBaseBVNodeList = L.sort $ filter ((> 1) . popCount) $ fmap (bvLabel . snd)  rightNodesInPrunedGraph

      in
      --if leftPrunedBVNodeList == rightPrunedBVNodeList then False
      -- else if leftBaseBVNodeList == rightBaseBVNodeList then False
      --else True
      True
      -}

-- | exchangePrunedGraphs creates a new "splitGraph" containing both first (base) and second (pruned) graph components
-- both components need to have HTU and edges reindexed to be in sync, oringal edge terminal node is also reindexed and returned for limit readd distance
exchangePrunedGraphs :: Int -> ((DecoratedGraph, LG.Node, LG.Node, LG.Node), (DecoratedGraph, LG.Node, LG.Node, LG.Node), LG.Node) -> (DecoratedGraph, Int , Int, Int, Int)
exchangePrunedGraphs numLeaves (firstGraphTuple, secondGraphTuple, breakEdgeNode) =
   let (firstSplitGraph, firstGraphRootIndex, _, _) = firstGraphTuple
       (secondSplitGraph, _, secondPrunedGraphRootIndex, _) = secondGraphTuple

       -- get nodes and edges of firstBase
       firstGraphRootLabel = fromJust $ LG.lab firstSplitGraph firstGraphRootIndex
       firstGraphRootNode = (firstGraphRootIndex, firstGraphRootLabel)
       (firstBaseGraphNodeList', firstBaseGraphEdgeList) = LG.nodesAndEdgesAfter firstSplitGraph [firstGraphRootNode]

       --add in root nodes of partitions since not included in "nodesAfter" function
       firstBaseGraphNodeList = firstGraphRootNode : firstBaseGraphNodeList'


       -- get nodes and edges of second pruned
       secondPrunedGraphRootLabel = fromJust $ LG.lab secondSplitGraph secondPrunedGraphRootIndex
       secondPrunedGraphRootNode = (secondPrunedGraphRootIndex, secondPrunedGraphRootLabel)
       secondPrunedParentNode = head $ LG.labParents secondSplitGraph secondPrunedGraphRootIndex
       (secondPrunedGraphNodeList', secondPrunedGraphEdgeList') = LG.nodesAndEdgesAfter secondSplitGraph [secondPrunedGraphRootNode]

       -- add root node of second pruned since not included in "nodesAfter" function
       -- add in gandparent nodes of pruned and its edges to pruned graphs
       secondPrunedGraphNodeList = [secondPrunedGraphRootNode, secondPrunedParentNode] ++ secondPrunedGraphNodeList'
       secondPrunedGraphEdgeList = (head $ LG.inn secondSplitGraph secondPrunedGraphRootIndex) : secondPrunedGraphEdgeList'

       -- reindex base and pruned partitions (HTUs and edges) to get in sync and make combinable
       -- 0 is dummy since won't be in base split
       (baseGraphNodes, baseGraphEdges, numBaseHTUs, reindexedBreakEdgeNodeBase) = reindexSubGraph numLeaves 0 firstBaseGraphNodeList firstBaseGraphEdgeList breakEdgeNode
       (prunedGraphNodes, prunedGraphEdges, _, _) = reindexSubGraph numLeaves numBaseHTUs secondPrunedGraphNodeList secondPrunedGraphEdgeList breakEdgeNode

       -- should always be in base graph--should be in first (base) component--if not use original node
       reindexedBreakEdgeNode = if (reindexedBreakEdgeNodeBase /= Nothing) then fromJust reindexedBreakEdgeNodeBase
                                else breakEdgeNode


       -- create and reindex new split graph
       newSplitGraph = LG.mkGraph (baseGraphNodes ++ prunedGraphNodes) (baseGraphEdges ++ prunedGraphEdges)

       -- get graph root Index, pruned root index, pruned root parent index
       -- firstGraphRootIndex should not have changed in reindexing--same as numLeaves
       prunedParentRootIndex = fst $ head $ (LG.getRoots newSplitGraph) L.\\ [firstGraphRootNode]
       prunedRootIndex = head $ LG.descendants newSplitGraph prunedParentRootIndex

   in
   if (length $ LG.getRoots newSplitGraph) /= 2 then error ("Not 2 components in split graph: " ++ "\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph newSplitGraph))
   else if (length $ LG.descendants newSplitGraph prunedParentRootIndex) /= 1 then error ("Too many children of parentPrunedNode: " ++ "\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph newSplitGraph))
   else if (length $ LG.parents secondSplitGraph secondPrunedGraphRootIndex) /= 1 then error ("Parent number not equal to 1 in node "
      ++ (show secondPrunedGraphRootIndex) ++ " of second graph\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph secondSplitGraph))
   else if (length $ LG.inn secondSplitGraph secondPrunedGraphRootIndex) /= 1 then error ("Edge incedent tor pruned graph not equal to 1 in node "
      ++ (show $ fmap LG.toEdge $  LG.inn secondSplitGraph secondPrunedGraphRootIndex) ++ " of second graph\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph secondSplitGraph))
   else
     {-
     } trace ("Nodes: " ++ (show (firstGraphRootIndex, prunedParentRootIndex, prunedRootIndex)) ++ " First Graph\n:" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph firstSplitGraph)
         ++ "\nSecond Graph\n:" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph secondSplitGraph)
         ++ "\nNew split graph\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph newSplitGraph)
         )
      -}
      (newSplitGraph, firstGraphRootIndex, prunedParentRootIndex, prunedRootIndex, reindexedBreakEdgeNode)


-- | reindexSubGraph reindexes the non-leaf nodes and edges of a subgraph to allow topological combination of subgraphs
-- the leaf indices are unchanges but HTUs are changes ot in order enumeration statting with an input offset
-- new BreakEdge is returned as a Maybe becuase may be either in base or pruned subgraphs
reindexSubGraph :: Int -> Int -> [LG.LNode VertexInfo] -> [LG.LEdge b] -> LG.Node -> ([LG.LNode VertexInfo], [LG.LEdge b], Int, Maybe LG.Node)
reindexSubGraph numLeaves offset nodeList edgeList origBreakEdge =
   if null nodeList || null edgeList then ([],[], offset, Nothing)
   else
      -- create map of node indices from list
      let (newNodeList, indexList) = unzip $ getPairList numLeaves offset nodeList
          indexMap = MAP.fromList indexList
          newEdgeList = fmap (reIndexEdge indexMap) edgeList
          newBreakEdge = MAP.lookup origBreakEdge indexMap
      in
      {-
      if newBreakEdge == Nothing then error  ("Map index for break edge node not found: " ++ (show origBreakEdge) ++ " in Map " ++ (show $ MAP.toList indexMap))
      else
      -}
         -- trace ("RISG:" ++ (show (fmap fst nodeList, fmap fst newNodeList, numLeaves)) ++ " map " ++ (show $ MAP.toList indexMap))
         (newNodeList, newEdgeList, 1 + (maximum $ fmap fst newNodeList) - numLeaves, newBreakEdge)

-- | reIndexEdge takes a map and a labelled edge and returns new indices same label edge based on map
reIndexEdge :: MAP.Map Int Int -> LG.LEdge b -> LG.LEdge b
reIndexEdge indexMap (u,v,l) =
   let u' = MAP.lookup u indexMap
       v' = MAP.lookup v indexMap
   in
   if u' == Nothing || v' == Nothing then error ("Error in map lookup in reindexEdge: " ++ show (u,v))
   else (fromJust u', fromJust v', l)


-- | getPairList returns an original index new index lits of pairs
-- assumes leaf nmodes are first numleaves
getPairList :: Int -> Int -> [LG.LNode VertexInfo] -> [(LG.LNode VertexInfo, (Int, Int))]
getPairList numLeaves counter nodeList =
   if null nodeList then []
   else
      let (firstIndex, firstLabel) = head nodeList
          newLabel = firstLabel {vertName = TL.pack ("HTU" ++ (show $ counter + numLeaves))}
      in
      if firstIndex < numLeaves then (head nodeList, (firstIndex, firstIndex)) : getPairList numLeaves counter (tail nodeList)
      else ((counter + numLeaves, newLabel) , (firstIndex, (counter + numLeaves))) : getPairList numLeaves (counter + 1) (tail nodeList)

