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

-- | fuseAllGraphs takes a list of phylogenetic graphs and performs all pairwise fuses
-- later--could limit by options making random choices for fusing
-- keeps results according to options (best, unique, etc)
-- unique is unique of "best" from individual fusings
-- singleRound short circuits recursive continuation on newly found graphs
fuseAllGraphs :: GlobalSettings 
              -> ProcessedData 
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
fuseAllGraphs inGS inData rSeed keepNum doNNI doSPR doTBR doSteepest doAll returnBest returnUnique singleRound inGraphList = 
   if null inGraphList then ([], 0)
   else 
      let -- getting values to be passed for graph diagnosis later
         numLeaves = V.length $ fst3 inData
         leafGraph = T.makeSimpleLeafGraph inData
         leafDecGraph = T.makeLeafGraph inData
         leafGraphSoftWired = T.makeLeafGraphSoftWired inData
         hasNonExactChars = U.getNumberNonExactCharacters (thd3 inData) > 0
         charInfoVV = six6 $ head inGraphList

         -- get fuse pairs
         graphPairList = getListPairs inGraphList

         newGraphList = concat (fmap (fusePair inGS inData numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV rSeed keepNum doNNI doSPR doTBR doSteepest doAll) graphPairList `using` PU.myParListChunkRDS)

      in
      trace ("\tFusing " ++ (show $ length graphPairList) ++ " graph pairs")    
      (newGraphList, 1)


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
         -> Bool 
         -> Bool 
         -> Bool 
         -> Bool 
         -> Bool 
         -> (PhylogeneticGraph, PhylogeneticGraph)
         -> [PhylogeneticGraph]
fusePair inGS inData numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV rSeed keepNum doNNI doSPR doTBR doSteepest doAll (leftGraph, rightGraph) =
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
          (_, leftGraphRootIndexList, leftPrunedGraphRootIndexList,  leftOriginalConnectionOfPruned) = L.unzip4 leftSplitTupleList
          --leftPrunedGraphRootIndexList = fmap thd4 leftSplitTupleList
          leftPrunedGraphBV = fmap bvLabel $ fmap fromJust $ fmap (LG.lab leftDecoratedGraph) leftPrunedGraphRootIndexList

          -- right graph splits
          rightDecoratedGraph = thd6 rightGraph
          (rightRootIndex, _) = head $ LG.getRoots rightDecoratedGraph
          rightBreakEdgeList = if (graphType inGS) == Tree then filter ((/= rightRootIndex) . fst3) $ LG.labEdges rightDecoratedGraph
                              else filter ((/= rightRootIndex) . fst3) $ GO.getEdgeSplitList rightDecoratedGraph
          rightSplitTupleList = fmap (GO.splitGraphOnEdge rightDecoratedGraph) rightBreakEdgeList `using` PU.myParListChunkRDS
          (_, rightGraphRootIndexList, rightPrunedGraphRootIndexList,  rightOriginalConnectionOfPruned) = L.unzip4 rightSplitTupleList
          -- rightPrunedGraphRootIndexList = fmap thd4 rightSplitTupleList
          rightPrunedGraphBV = fmap bvLabel $ fmap fromJust $ fmap (LG.lab rightDecoratedGraph) rightPrunedGraphRootIndexList

          -- get compatible split pairs
          -- check bv of root index of pruned subgraphs
          leftRightMatchList = zipWith (==) leftPrunedGraphBV rightPrunedGraphBV
          leftRightOrList =  zipWith (.|.) leftPrunedGraphBV rightPrunedGraphBV
          (leftValidTupleList, rightValidTupleList, _, _)  = L.unzip4 $ filter ((> 2) . popCount . fth4) $ filter ((== True) . thd4) $ L.zip4 leftSplitTupleList rightSplitTupleList leftRightMatchList leftRightOrList

          -- create new "splitgraphs" by replacing nodes and edges of pruned subgraph in reciprocal graphs
          leftBaseRightPrunedSplitGraphList = fmap exchangePrunedGraphs (zip leftValidTupleList rightValidTupleList) `using` PU.myParListChunkRDS
          rightBaseLeftPrunedSplitGraphList = fmap exchangePrunedGraphs (zip rightValidTupleList leftValidTupleList) `using` PU.myParListChunkRDS

          -- reoptimize splitGraphs so ready for readdition
          -- False for doIA

          leftRightOptimizedSplitGraphCostList = fmap (S.reoptimizeSplitGraphFromVertexTuple inGS inData False charInfoVV) (zip3 leftBaseRightPrunedSplitGraphList leftGraphRootIndexList rightPrunedGraphRootIndexList) `using` PU.myParListChunkRDS

          rightLeftOptimizedSplitGraphCostList = fmap (S.reoptimizeSplitGraphFromVertexTuple inGS inData False charInfoVV) (zip3 rightBaseLeftPrunedSplitGraphList rightGraphRootIndexList leftPrunedGraphRootIndexList) `using` PU.myParListChunkRDS


      in
      trace ("FP: " ++ (show $ length leftValidTupleList) ++ "\nLeftRight splitCost " ++ (show $ fmap snd leftRightOptimizedSplitGraphCostList)
         ++ "\nrightLeft splitCost " ++ (show $ fmap snd rightLeftOptimizedSplitGraphCostList)) 
      [leftGraph, rightGraph]  

-- | exchangePrunedGraphs creates a new "splitGraph" containing both first base and second pruned graph components
exchangePrunedGraphs :: ((DecoratedGraph, LG.Node, LG.Node, LG.Node), (DecoratedGraph, LG.Node, LG.Node, LG.Node)) -> DecoratedGraph
exchangePrunedGraphs (firstGraphTuple, secondGraphTuple) =
   let (firstSplitGraph, firstGraphRootIndex, firstPrunedGraphRootIndex, firstOriginalConnectionOfPruned) = firstGraphTuple
       (secondSplitGraph, secondGraphRootIndex, secondPrunedGraphRootIndex, secondOriginalConnectionOfPruned) = secondGraphTuple

       -- get nodes and edges of firstBase and secondPruned graphs
       -- need to add in root nodes of partitions since not included in "nodesAfter" function
       -- need to add in gandparent nodes of pruned and its edges to pruned graphs

       firstGraphRootLabel = fromJust $ LG.lab firstSplitGraph firstGraphRootIndex
       firstGraphRootNode = (firstGraphRootIndex, firstGraphRootLabel)
       (firstBaseGraphNodeList', firstBaseGraphEdgeList) = LG.nodesAndEdgesAfter firstSplitGraph [firstGraphRootNode]
       firstBaseGraphNodeList = firstGraphRootNode : firstBaseGraphNodeList'
   
       secondPrunedGraphRootLabel = fromJust $ LG.lab secondSplitGraph secondPrunedGraphRootIndex
       secondPrunedGraphRootNode = (secondPrunedGraphRootIndex, secondPrunedGraphRootLabel) 
       secondPrunedParentNode = head $ LG.labParents secondSplitGraph secondPrunedGraphRootIndex
       (secondPrunedGraphNodeList', secondPrunedGraphEdgeList') = LG.nodesAndEdgesAfter secondSplitGraph [secondPrunedGraphRootNode]
       secondPrunedGraphNodeList = [secondPrunedGraphRootNode, secondPrunedParentNode] ++ secondPrunedGraphNodeList'
       secondPrunedGraphEdgeList = (head $ LG.inn secondSplitGraph secondPrunedGraphRootIndex) : secondPrunedGraphEdgeList'
   
       -- create and reindex new split graph
       newSplitGraph = LG.mkGraph (firstBaseGraphNodeList ++ secondPrunedGraphNodeList) (firstBaseGraphEdgeList ++ secondPrunedGraphEdgeList) 

   in
   if (length $ LG.labParents secondSplitGraph secondPrunedGraphRootIndex) /= 1 then error ("Parent number not equal to 1 in node " 
      ++ (show secondPrunedGraphRootIndex) ++ " of second graph\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph secondSplitGraph))
   else if (length $ LG.inn secondSplitGraph secondPrunedGraphRootIndex) /= 1 then error ("Edge incedent tor p[runed graph not equal to 1 in node " 
      ++ (show $ fmap LG.toEdge $  LG.inn secondSplitGraph secondPrunedGraphRootIndex) ++ " of second graph\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph secondSplitGraph))
   else
      trace ("First Graph\n:" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph firstSplitGraph)
         ++ "\nSecond Graph\n:" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph secondSplitGraph)
         ++ "\nNew split graph\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph newSplitGraph)
         )
      newSplitGraph