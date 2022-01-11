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
          (leftSplitGraphList, leftGraphRootIndexList, leftPrunedGraphRootIndexList,  leftOriginalConnectionOfPruned) = L.unzip4 leftSplitTupleList
          --leftPrunedGraphRootIndexList = fmap thd4 leftSplitTupleList
          leftPrunedGraphBVList = fmap bvLabel $ fmap fromJust $ fmap (LG.lab leftDecoratedGraph) leftPrunedGraphRootIndexList

          -- right graph splits
          rightDecoratedGraph = thd6 rightGraph
          (rightRootIndex, _) = head $ LG.getRoots rightDecoratedGraph
          rightBreakEdgeList = if (graphType inGS) == Tree then filter ((/= rightRootIndex) . fst3) $ LG.labEdges rightDecoratedGraph
                              else filter ((/= rightRootIndex) . fst3) $ GO.getEdgeSplitList rightDecoratedGraph
          rightSplitTupleList = fmap (GO.splitGraphOnEdge rightDecoratedGraph) rightBreakEdgeList `using` PU.myParListChunkRDS
          (rightSplitGraphList, rightGraphRootIndexList, rightPrunedGraphRootIndexList,  rightOriginalConnectionOfPruned) = L.unzip4 rightSplitTupleList
          -- rightPrunedGraphRootIndexList = fmap thd4 rightSplitTupleList
          rightPrunedGraphBVList = fmap bvLabel $ fmap fromJust $ fmap (LG.lab rightDecoratedGraph) rightPrunedGraphRootIndexList

          -- get compatible split pairs via checking bv of root index of pruned subgraphs
          leftRightMatchList = zipWith (==) leftPrunedGraphBVList rightPrunedGraphBVList

          -- only take compatible, non-identical pairs with > 1 terminal--otherwise basically SPR move or nmothing (if identical)
          -- (leftValidTupleList, rightValidTupleList, _, _)  = L.unzip4 $ filter ((> 2) . popCount . fth4) $ filter ((== True) . thd4) $ L.zip4 leftSplitTupleList rightSplitTupleList leftRightMatchList leftRightOrList
          switchableList = L.zipWith4 (getCompatibleNonIdenticalSplits numLeaves) leftSplitTupleList rightSplitTupleList leftRightMatchList leftPrunedGraphBVList
          (leftValidTupleList, rightValidTupleList, _) = L.unzip3 $ filter ((==True) . thd3) $ zip3 leftSplitTupleList rightSplitTupleList switchableList
          -- (leftValidTupleList, rightValidTupleList) = L.unzip $ concat $ L.zipWith4 (getCompatibleNonIdenticalSplits numLeaves) leftSplitTupleList rightSplitTupleList leftRightMatchList leftPrunedGraphBV


          -- create new "splitgraphs" by replacing nodes and edges of pruned subgraph in reciprocal graphs
          leftBaseRightPrunedSplitGraphList = fmap (exchangePrunedGraphs numLeaves) (zip leftValidTupleList rightValidTupleList) `using` PU.myParListChunkRDS
          rightBaseLeftPrunedSplitGraphList = fmap (exchangePrunedGraphs numLeaves) (zip rightValidTupleList leftValidTupleList) `using` PU.myParListChunkRDS

          -- reoptimize splitGraphs so ready for readdition
          -- False for doIA
          leftRightOptimizedSplitGraphCostList = fmap (S.reoptimizeSplitGraphFromVertexTuple inGS inData False charInfoVV) (zip3 leftBaseRightPrunedSplitGraphList leftGraphRootIndexList rightPrunedGraphRootIndexList) `using` PU.myParListChunkRDS

          rightLeftOptimizedSplitGraphCostList = fmap (S.reoptimizeSplitGraphFromVertexTuple inGS inData False charInfoVV) (zip3 rightBaseLeftPrunedSplitGraphList rightGraphRootIndexList leftPrunedGraphRootIndexList) `using` PU.myParListChunkRDS

          -- perform swap-like fuse operations
          -- needs to be abstracted outn to a function given all th ecomplex arguments
          {-
          swapType = if doTBR then "tbr"
                     else if doSPPR then "spr"
                     else if doNNI then "nni"
                     else "nni"
          maxMoveEdgeDist' = if not doTBR && not doSPR && not doNNI then 2
                             else maxMoveEdgeDist
          doIA = False
          
          leftRightFuseList = if doSteepest then rejoinGraphKeepBestSteepest inGS inData swapType (snd6 leftGraph) numToKeep maxMoveEdgeDist' True doIA (six6 leftGraph) $ L.zip5 leftRightOptimizedSplitGraphCostList leftGraphRootIndexList prunedGraphRootIndexList originalConnectionOfPruned (fmap head $ fmap ((LG.parents $ thd6 firstGraph).fst3) breakEdgeList)
                              else error "Blah"
          -}


          -- | get fuse graphs via swap function
      in
      if null leftValidTupleList then []
      else 
         trace ("FP: " ++ (show $ length leftValidTupleList) ++ "num (Left,Right) " ++ (show (length leftBaseRightPrunedSplitGraphList, length rightBaseLeftPrunedSplitGraphList)) 
            ++ "\nLeftRight splitCost " ++ (show $ fmap snd leftRightOptimizedSplitGraphCostList)
            ++ "\nrightLeft splitCost " ++ (show $ fmap snd rightLeftOptimizedSplitGraphCostList)) 
         [leftGraph, rightGraph]  

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
-- both components need to have HTU and edges reindexed to be in sync
exchangePrunedGraphs :: Int -> ((DecoratedGraph, LG.Node, LG.Node, LG.Node), (DecoratedGraph, LG.Node, LG.Node, LG.Node)) -> DecoratedGraph
exchangePrunedGraphs numLeaves (firstGraphTuple, secondGraphTuple) =
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
       (baseGraphNode, baseGraphEdges, numBaseHTUs) = reindexSubGraph numLeaves 0 firstBaseGraphNodeList firstBaseGraphEdgeList
       (prunedGraphNode, prunedGraphEdges, _) = reindexSubGraph numLeaves numBaseHTUs secondPrunedGraphNodeList secondPrunedGraphEdgeList
   
       -- create and reindex new split graph
       newSplitGraph = LG.mkGraph (baseGraphNode ++ prunedGraphNode) (baseGraphEdges ++ prunedGraphEdges) 

   in
   if (length $ LG.parents secondSplitGraph secondPrunedGraphRootIndex) /= 1 then error ("Parent number not equal to 1 in node " 
      ++ (show secondPrunedGraphRootIndex) ++ " of second graph\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph secondSplitGraph))
   else if (length $ LG.inn secondSplitGraph secondPrunedGraphRootIndex) /= 1 then error ("Edge incedent tor p[runed graph not equal to 1 in node " 
      ++ (show $ fmap LG.toEdge $  LG.inn secondSplitGraph secondPrunedGraphRootIndex) ++ " of second graph\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph secondSplitGraph))
   else
      trace ("First Graph\n:" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph firstSplitGraph)
         ++ "\nSecond Graph\n:" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph secondSplitGraph)
         ++ "\nNew split graph\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph newSplitGraph)
         )
      newSplitGraph


-- | reindexSubGraph reindexes the non-leaf nodes and edges of a subgraph to allow topological combination of subgraphs 
-- the leaf indices are unchanges but HTUs are changes ot in order enumeration statting with an input offset
reindexSubGraph :: Int -> Int -> [LG.LNode VertexInfo] -> [LG.LEdge b] -> ([LG.LNode VertexInfo], [LG.LEdge b], Int)
reindexSubGraph numLeaves offset nodeList edgeList =
   if null nodeList || null edgeList then ([],[], offset)
   else 
      -- create map of node indices from list
      let (newNodeList, indexList) = unzip $ getPairList numLeaves offset nodeList
          indexMap = MAP.fromList indexList
          newEdgeList = fmap (reIndexEdge indexMap) edgeList
      in
      trace ("RISG:" ++ (show (fmap fst newNodeList, numLeaves)))
      (newNodeList, newEdgeList, 1 + (maximum $ fmap fst newNodeList) - numLeaves)

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
