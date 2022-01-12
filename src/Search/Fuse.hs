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

import Types.Types
import qualified ParallelUtilities       as PU
import Control.Parallel.Strategies
import GeneralUtilities
import qualified Graphs.GraphOperations  as GO
import qualified Utilities.LocalGraph    as LG
import Utilities.Utilities               as U
import Debug.Trace
import           Data.Char
import           Text.Read
import           Data.Maybe
import qualified GraphOptimization.Traversals as T
import qualified Data.Vector as V
import qualified GraphOptimization.PreOrderFunctions as PRE
import qualified GraphOptimization.PostOrderFunctions as POS
import qualified Data.List as L
import qualified Data.Text.Lazy              as TL
import qualified GraphOptimization.Medians as M
import qualified Search.Swap as S
import Data.Bits
import qualified Data.Map as MAP
import qualified Data.Text.Lazy              as TL
import qualified Data.BitVector.LittleEndian as BV

-- | fuseAllGraphs takes a list of phylogenetic graphs and performs all pairwise fuses
-- later--could limit by options making random choices for fusing
-- keeps results according to options (best, unique, etc)
-- unique is unique of "best" from individual fusings
-- singleRound short circuits recursive continuation on newly found graphs
fuseAllGraphs :: GlobalSettings 
              -> ProcessedData 
              -> Int 
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
              -> [PhylogeneticGraph] 
              -> ([PhylogeneticGraph], Int)
fuseAllGraphs inGS inData rSeed keepNum maxMoveEdgeDist counter doNNI doSPR doTBR doSteepest doAll returnBest returnUnique singleRound inGraphList = 
   if null inGraphList then ([], 0)
   else 
      let -- getting values to be passed for graph diagnosis later
         numLeaves = V.length $ fst3 inData
         leafGraph = T.makeSimpleLeafGraph inData
         leafDecGraph = T.makeLeafGraph inData
         leafGraphSoftWired = T.makeLeafGraphSoftWired inData
         hasNonExactChars = U.getNumberNonExactCharacters (thd3 inData) > 0
         charInfoVV = six6 $ head inGraphList

         curBest = minimum $ fmap snd6 inGraphList

         -- get fuse pairs
         graphPairList = getListPairs inGraphList

         newGraphList = concat (fmap (fusePair inGS inData numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV rSeed keepNum maxMoveEdgeDist doNNI doSPR doTBR doSteepest doAll) graphPairList `using` PU.myParListChunkRDS)

         fuseBest = if not (null newGraphList) then  minimum $ fmap snd6 newGraphList
                    else infinity

      in
      trace ("\tFusing " ++ (show $ length graphPairList) ++ " graph pairs") (
      if null newGraphList then (inGraphList, counter + 1)
      else if returnUnique then (GO.selectPhylogeneticGraph [("unique", (show keepNum))] 0 ["unique"] (inGraphList ++ newGraphList), counter + 1)
      else -- return best
         -- only do one round of fusing 
         if singleRound then (GO.selectPhylogeneticGraph [("best", (show keepNum))] 0 ["best"] (inGraphList ++ newGraphList), counter + 1)

         -- recursive rounds
         else 
            let allBestList = GO.selectPhylogeneticGraph [("best", (show keepNum))] 0 ["best"] (inGraphList ++ newGraphList)
            in
            
            -- found worse
            if fuseBest > curBest then (allBestList, counter + 1)

            -- found better   
            else if fuseBest < curBest then fuseAllGraphs inGS inData rSeed keepNum maxMoveEdgeDist (counter + 1) doNNI doSPR doTBR doSteepest doAll returnBest returnUnique singleRound allBestList
            
            -- equal cost
            else 
               -- found novel graph of equal cost
               if length allBestList /= length inGraphList then fuseAllGraphs inGS inData rSeed keepNum maxMoveEdgeDist (counter + 1) doNNI doSPR doTBR doSteepest doAll returnBest returnUnique singleRound allBestList

               -- nothing novel
               else (allBestList, counter + 1)
      )


-- | fusePair recombines a single pair of graphs 
-- this is done by cooptinmg the split and readd functinos from the Swap.Swap functions and exchanging 
-- pruned subgraphs with the same leaf complement (as recorded by the subbtree roiot node bit vector field)
-- spr-like and tbr-like readds can be perfomred as with options
fusePair :: GlobalSettings 
         -> ProcessedData 
         -> Int 
         -> SimpleGraph
         -> DecoratedGraph
         -> DecoratedGraph
         -> Bool
         -> V.Vector (V.Vector CharInfo) 
         -> Int 
         -> Int 
         -> Int 
         -> Bool 
         -> Bool 
         -> Bool 
         -> Bool 
         -> Bool 
         -> (PhylogeneticGraph, PhylogeneticGraph)
         -> [PhylogeneticGraph]
fusePair inGS inData numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV rSeed keepNum maxMoveEdgeDist doNNI doSPR doTBR doSteepest doAll (leftGraph, rightGraph) =
   if (LG.isEmpty $ fst6 leftGraph) || (LG.isEmpty $ fst6 rightGraph) then error "Empty graph in fusePair"
   else if (fst6 leftGraph) == (fst6 rightGraph) then []
   else
      -- split graphs at all bridge edges (all edges for Tree) 
      let -- left graph splits
          leftDecoratedGraph = thd6 leftGraph
          (leftRootIndex, _) = head $ LG.getRoots leftDecoratedGraph
          leftBreakEdgeList = if (graphType inGS) == Tree then filter ((/= leftRootIndex) . fst3) $ LG.labEdges leftDecoratedGraph
                              else filter ((/= leftRootIndex) . fst3) $ GO.getEdgeSplitList leftDecoratedGraph
          leftSplitTupleList = fmap (GO.splitGraphOnEdge leftDecoratedGraph) leftBreakEdgeList `using` PU.myParListChunkRDS
          (leftSplitGraphList, leftGraphRootIndexList, leftPrunedGraphRootIndexList,  leftOriginalConnectionOfPrunedList) = L.unzip4 leftSplitTupleList
          --leftPrunedGraphRootIndexList = fmap thd4 leftSplitTupleList
          leftPrunedGraphBVList = fmap bvLabel $ fmap fromJust $ fmap (LG.lab leftDecoratedGraph) leftPrunedGraphRootIndexList
          

          -- right graph splits
          rightDecoratedGraph = thd6 rightGraph
          (rightRootIndex, _) = head $ LG.getRoots rightDecoratedGraph
          rightBreakEdgeList = if (graphType inGS) == Tree then filter ((/= rightRootIndex) . fst3) $ LG.labEdges rightDecoratedGraph
                              else filter ((/= rightRootIndex) . fst3) $ GO.getEdgeSplitList rightDecoratedGraph
          rightSplitTupleList = fmap (GO.splitGraphOnEdge rightDecoratedGraph) rightBreakEdgeList `using` PU.myParListChunkRDS
          (rightSplitGraphList, rightGraphRootIndexList, rightPrunedGraphRootIndexList,  rightOriginalConnectionOfPrunedList) = L.unzip4 rightSplitTupleList
          -- rightPrunedGraphRootIndexList = fmap thd4 rightSplitTupleList
          rightPrunedGraphBVList = fmap bvLabel $ fmap fromJust $ fmap (LG.lab rightDecoratedGraph) rightPrunedGraphRootIndexList
          

          -- need to get all pairs of split graphs
          (leftSplitTupleList', rightSplitTupleList') =  unzip $ cartProd leftSplitTupleList rightSplitTupleList
          (leftPrunedGraphBVList', rightPrunedGraphBVList') = unzip $ cartProd leftPrunedGraphBVList rightPrunedGraphBVList
          -- (leftBaseBVList, rightBaseBVList) = unzip $ cartProd leftBaseGraphBVList rightBaseGraphBVList


          -- get compatible split pairs via checking bv of root index of pruned subgraphs
          leftRightMatchList = zipWith (==) leftPrunedGraphBVList' rightPrunedGraphBVList'

          -- only take compatible, non-identical pairs with > 1 terminal--otherwise basically SPR move or nmothing (if identical)
          recombinablePairList = L.zipWith4 (getCompatibleNonIdenticalSplits numLeaves) leftSplitTupleList' rightSplitTupleList' leftRightMatchList leftPrunedGraphBVList'
          (leftValidTupleList, rightValidTupleList, _) = L.unzip3 $ filter ((==True) . thd3) $ zip3 leftSplitTupleList' rightSplitTupleList' recombinablePairList
          

          -- create new "splitgraphs" by replacing nodes and edges of pruned subgraph in reciprocal graphs
          -- retuns reindexed list of base graph root, pruned component root,  parent of pruned component root, original graph break edge
          (leftBaseRightPrunedSplitGraphList, leftRightGraphRootIndexList, leftRightPrunedParentRootIndexList, leftRightPrunedRootIndexList, leftRightOriginalConnectionOfPrunedList) = L.unzip5 (fmap (exchangePrunedGraphs numLeaves) (zip3 leftValidTupleList rightValidTupleList leftOriginalConnectionOfPrunedList) `using` PU.myParListChunkRDS)
          (rightBaseLeftPrunedSplitGraphList, rightLeftGraphRootIndexList, rightLeftPrunedParentRootIndexList, rightLeftPrunedRootIndexList, rightLeftOriginalConnectionOfPrunedList) = L.unzip5 (fmap (exchangePrunedGraphs numLeaves) (zip3 rightValidTupleList leftValidTupleList rightOriginalConnectionOfPrunedList) `using` PU.myParListChunkRDS)

          -- reoptimize splitGraphs so ready for readdition--using updated base and prune indices
          -- False for doIA
          leftRightOptimizedSplitGraphCostList = fmap (S.reoptimizeSplitGraphFromVertexTuple inGS inData False charInfoVV) (zip3 leftBaseRightPrunedSplitGraphList leftRightGraphRootIndexList leftRightPrunedRootIndexList) `using` PU.myParListChunkRDS

          rightLeftOptimizedSplitGraphCostList = fmap (S.reoptimizeSplitGraphFromVertexTuple inGS inData False charInfoVV) (zip3 rightBaseLeftPrunedSplitGraphList rightLeftGraphRootIndexList rightLeftPrunedRootIndexList) `using` PU.myParListChunkRDS

          -- Check if base graphs are different as well (nneded to be reoptimized to get base root bv)
          -- otherwise no point in recombination
          {-
          leftBaseGraphBVList = fmap bvLabel $ fmap fromJust $ zipWith LG.lab (fmap fst leftRightOptimizedSplitGraphCostList) leftRightGraphRootIndexList
          rightBaseGraphBVList =fmap bvLabel $ fmap fromJust $ zipWith LG.lab (fmap fst rightLeftOptimizedSplitGraphCostList) rightLeftGraphRootIndexList
          baseGraphDifferentList = zipWith (/=) leftBaseBVList rightBaseBVList 
          -}
          baseGraphDifferentList = L.replicate (length leftRightOptimizedSplitGraphCostList) True 

          (_, leftRightOptimizedSplitGraphCostList', leftRightGraphRootIndexList', leftRightPrunedRootIndexList', leftRightPrunedParentRootIndexList', leftRightOriginalConnectionOfPrunedList') = L.unzip6 $ filter ((== True) . fst6) $ L.zip6 baseGraphDifferentList leftRightOptimizedSplitGraphCostList leftRightGraphRootIndexList leftRightPrunedRootIndexList leftRightPrunedParentRootIndexList leftRightOriginalConnectionOfPrunedList

          (_, rightLeftOptimizedSplitGraphCostList', rightLeftGraphRootIndexList', rightLeftPrunedRootIndexList', rightLeftPrunedParentRootIndexList', rightLeftOriginalConnectionOfPrunedList') = L.unzip6 $ filter ((== True) . fst6) $ L.zip6 baseGraphDifferentList rightLeftOptimizedSplitGraphCostList rightLeftGraphRootIndexList rightLeftPrunedRootIndexList rightLeftPrunedParentRootIndexList rightLeftOriginalConnectionOfPrunedList

          -- re-add pruned component to base component left-right and right-left
          -- need cure best cost
          curBetterCost = min (snd6 leftGraph) (snd6 rightGraph)
          charInfoVV = six6 leftGraph

          leftRightFusedGraphList = recombineComponents inGS inData keepNum maxMoveEdgeDist doNNI doSPR doTBR doSteepest doAll charInfoVV curBetterCost leftRightOptimizedSplitGraphCostList' leftRightGraphRootIndexList' leftRightPrunedRootIndexList' leftRightPrunedParentRootIndexList' leftRightOriginalConnectionOfPrunedList'
          rightLeftFusedGraphList = recombineComponents inGS inData keepNum maxMoveEdgeDist doNNI doSPR doTBR doSteepest doAll charInfoVV curBetterCost rightLeftOptimizedSplitGraphCostList' rightLeftGraphRootIndexList' rightLeftPrunedRootIndexList' rightLeftPrunedParentRootIndexList' rightLeftOriginalConnectionOfPrunedList'


          -- get "best" fused graphs from leftRight and rightLeft
          bestFusedGraphs = GO.selectPhylogeneticGraph [("best", (show keepNum))] 0 ["best"] (leftRightFusedGraphList ++ rightLeftFusedGraphList)

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
recombineComponents :: GlobalSettings 
                    -> ProcessedData
                    -> Int
                    -> Int
                    -> Bool 
                    -> Bool 
                    -> Bool 
                    -> Bool 
                    -> Bool
                    -> V.Vector (V.Vector CharInfo) 
                    -> VertexCost
                    -> [(DecoratedGraph, VertexCost)]
                    -> [Int]
                    -> [Int]
                    -> [Int]
                    -> [Int]
                    -> [PhylogeneticGraph]
recombineComponents inGS inData numToKeep inMaxMoveEdgeDist doNNI doSPR doTBR doSteepest doAll charInfoVV curBestCost splitGraphCostPairList baseRootIndexList prunedRootIndexList prunedParentRootIndexList originalConnectionOfPrunedComponentList = 
   -- check and see if any reconnecting to do
   -- trace ("RecombineComponents " ++ (show $ length splitGraphCostPairList)) (
   if null splitGraphCostPairList then []
   else
      let swapType = if doTBR then "tbr"
                     else if doSPR then "spr"
                     else if doNNI then "nni"
                     else "spr" -- will be set with 2 as maxMoveEdgeDist
          maxMoveEdgeDist = if not doTBR && not doSPR && not doNNI then 2
                          else inMaxMoveEdgeDist
          doIA = False --- since splits not created together, IA won't be consistent between components

          graphDataList = L.zip5 splitGraphCostPairList baseRootIndexList prunedRootIndexList prunedParentRootIndexList originalConnectionOfPrunedComponentList

          -- do "all additions" --steepest really doens't have meaning here since will not pop out to recombine on new graph
          -- False for doSteepest
          recombinedSimpleGraphCostPairList = concat (fmap (S.rejoinGraphKeepBestTuple inGS swapType curBestCost numToKeep maxMoveEdgeDist False doIA charInfoVV) graphDataList `using` PU.myParListChunkRDS)

          -- this based on heuristic deltas
          bestFuseCost = minimum $ fmap snd recombinedSimpleGraphCostPairList
          bestFuseSimpleGraphs = fmap fst $ filter ((== bestFuseCost) . snd) recombinedSimpleGraphCostPairList
          
      in
      --trace ("Checking in fusing") (
      if null recombinedSimpleGraphCostPairList then []
      else if bestFuseCost <= curBestCost then
         let rediagnodedGraphList = fmap (T.multiTraverseFullyLabelGraph inGS inData False False Nothing) bestFuseSimpleGraphs `using` PU.myParListChunkRDS
             bestRediagnosedGraphList = GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] rediagnodedGraphList
         in
         if (snd6 $ head bestRediagnosedGraphList) <= curBestCost then bestRediagnosedGraphList
         else []
      else []
      -- )
      -- )


   



-- | getCompatibleNonIdenticalSplits takes the number of leaves, splitGraph of the left graph, the splitGraph if the right graph, 
-- the bitVector equality list of pruned roots, the bitvector of the root of the pruned graph on left 
-- (this could be either since filter for identity--just to check leaf numbers)
-- checks that the leaf sets of the pruned suubgraphs are equal, greater than 1 leaf, fewer thanm nuleaves - 2, and non-identical
getCompatibleNonIdenticalSplits :: Int 
                                -> (DecoratedGraph, LG.Node, LG.Node, LG.Node) 
                                -> (DecoratedGraph, LG.Node, LG.Node, LG.Node) 
                                -> Bool 
                                -> BV.BitVector 
                                -> Bool
getCompatibleNonIdenticalSplits numLeaves leftSplitTuple rightSplitTuple leftRightMatch leftPrunedGraphBV = 
   
   if not leftRightMatch then False
   else if popCount leftPrunedGraphBV < 2 then False
   else if popCount leftPrunedGraphBV > (numLeaves - 3) then False 
   else 
      let (leftNodesInPrunedGraph, _) = LG.nodesAndEdgesAfter (fst4 leftSplitTuple) [((thd4 leftSplitTuple), fromJust $ LG.lab (fst4 leftSplitTuple) (thd4 leftSplitTuple))]
          leftBVNodeList = L.sort $ filter ((> 1) . popCount) $ fmap (bvLabel . snd)  leftNodesInPrunedGraph
          (rightNodesInPrunedGraph, _) = LG.nodesAndEdgesAfter (fst4 rightSplitTuple) [((thd4 rightSplitTuple), fromJust $ LG.lab (fst4 rightSplitTuple) (thd4 rightSplitTuple))]
          rightBVNodeList = L.sort $ filter ((> 1) . popCount) $ fmap (bvLabel . snd)  rightNodesInPrunedGraph
      in
      if leftBVNodeList == rightBVNodeList then False
      else True
        
           


-- | exchangePrunedGraphs creates a new "splitGraph" containing both first (base) and second (pruned) graph components
-- both components need to have HTU and edges reindexed to be in sync, oringal edge terminal node is also reindexed and returned for limit readd distance 
exchangePrunedGraphs :: Int -> ((DecoratedGraph, LG.Node, LG.Node, LG.Node), (DecoratedGraph, LG.Node, LG.Node, LG.Node), LG.Node) -> (DecoratedGraph, Int , Int, Int, Int)
exchangePrunedGraphs numLeaves (firstGraphTuple, secondGraphTuple, breakEdgeNode) =
   let (firstSplitGraph, firstGraphRootIndex, firstPrunedGraphRootIndex, firstOriginalConnectionOfPruned) = firstGraphTuple
       (secondSplitGraph, secondGraphRootIndex, secondPrunedGraphRootIndex, secondOriginalConnectionOfPruned) = secondGraphTuple

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
       (prunedGraphNodes, prunedGraphEdges, _, reindexedBreakEdgeNodePruned) = reindexSubGraph numLeaves numBaseHTUs secondPrunedGraphNodeList secondPrunedGraphEdgeList breakEdgeNode

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
   if u' == Nothing || v' == Nothing then error ("Eror in map lookup in reindexEdge: " ++ show (u,v))
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
