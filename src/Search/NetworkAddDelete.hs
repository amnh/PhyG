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

-- Can create a delta (outgroup rooted cost) before and after insertion with incremental post-order based on 
-- nodes to optimize being total nodes L.\\ the descendents of the two v and v'
-- would be a "true" value that ciould then be rerooted etc

-- | insertAllNetEdges deletes network edges one each each round until no better or additional 
-- graphs are found
-- call with ([], infinity) [single input graph]
insertAllNetEdges :: GlobalSettings -> ProcessedData -> Int -> ([PhylogeneticGraph], VertexCost) -> [PhylogeneticGraph] -> [PhylogeneticGraph]
insertAllNetEdges inGS inData numToKeep (curBestGraphList, curBestGraphCost) inPhyloGraphList =
   if null inPhyloGraphList then take numToKeep curBestGraphList
   else
      let currentCost = min curBestGraphCost (snd6 $ head inPhyloGraphList)
          (newGraphList, newGraphCost) = insertEachNetEdge inGS inData numToKeep (head inPhyloGraphList)
      in
      -- worse graphs found--go on 
      if  newGraphCost > currentCost then insertAllNetEdges inGS inData numToKeep ((head inPhyloGraphList) : curBestGraphList, currentCost) (tail inPhyloGraphList)

      -- "steepest style descent" abandons existing list if better cost found
      else if newGraphCost < currentCost then insertAllNetEdges inGS inData numToKeep (newGraphList, newGraphCost) newGraphList 

      -- equal cost
      -- not sure if should add new graphs to queue to do edge deletion again
      else 
         -- new grapjh list contains the input graph if equal and filterd unique already in deleteEachNetEdge
         let newCurSameBestList = GO.getBVUniqPhylogeneticGraph True (curBestGraphList ++ newGraphList)
         in
         insertAllNetEdges inGS inData numToKeep (newCurSameBestList, currentCost) (tail inPhyloGraphList)

-- | insertEachNetEdge takes a phylogenetic graph and inserts all permissible network edges one at time
-- and returns best list of new Phylogenetic Graphs and cost
-- even if worse--could be used for simulated annealing later
-- if equal returns unique graph list
insertEachNetEdge :: GlobalSettings -> ProcessedData -> Int -> PhylogeneticGraph -> ([PhylogeneticGraph], VertexCost)
insertEachNetEdge inGS inData numToKeep inPhyloGraph =
   if LG.isEmpty $ thd6 inPhyloGraph then error "Empty input insertEachNetEdge graph in deleteAllNetEdges"
   else
      let currentCost = snd6 inPhyloGraph

          candidateNetworkEdgeList = getPermissibleEdgePairs (thd6 inPhyloGraph)

          newGraphList = fmap (insertNetEdge inGS inData inPhyloGraph) candidateNetworkEdgeList `using`  PU.myParListChunkRDS

          minCostGraphList = GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] newGraphList
          minCost = snd6 $ head minCostGraphList
      in
      -- no network edges to delete
      if null candidateNetworkEdgeList then ([inPhyloGraph], currentCost)
      else if minCost /= currentCost then (minCostGraphList, minCost)
      else (GO.selectPhylogeneticGraph [("unique", (show numToKeep))] 0 ["unique"] $ inPhyloGraph : minCostGraphList, currentCost)

-- | getPermissibleEdgePairs takes a DecoratedGraph and returns the list of all pairs
-- of edges that can be joined by a network edge and meet all necessary conditions
getPermissibleEdgePairs :: DecoratedGraph -> [(LG.LEdge EdgeInfo, LG.LEdge EdgeInfo)]
getPermissibleEdgePairs inGraph = 
   if LG.isEmpty inGraph then error "Empty input graph in isEdgePairPermissible"
   else 
       let edgeList = LG.labEdges inGraph
           edgePairs = cartProd edgeList edgeList
           contraintList = getGraphCoevalConstraints inGraph
           edgeTestList = fmap (isEdgePairPermissible inGraph contraintList) edgePairs
           pairList = filter ((== True) . snd) $ zip edgePairs edgeTestList
       in
       fmap fst pairList

-- | isEdgePairPermissible takes a graph and two edges, coeval contraints, and tests whether a
-- pair of edges can be linked by a new edge and satify three consitions:
--    1) neither edge is a network edge
--    2) one edge cannot be "before" while the other is "after" in any of the constraint pairs
--    3) neither neither edge is an ancestor or descndent edge of the other (tested via bv of nodes)
-- the result should apply to a new edge in either direction
isEdgePairPermissible :: DecoratedGraph -> [([LG.LEdge EdgeInfo],[LG.LEdge EdgeInfo])] -> (LG.LEdge EdgeInfo, LG.LEdge EdgeInfo) -> Bool
isEdgePairPermissible inGraph constraintList (edge1, edge2) =
   if LG.isEmpty inGraph then error "Empty input graph in isEdgePairPermissible"
   else 
       if (LG.isNetworkLabEdge inGraph edge1) || (LG.isNetworkLabEdge inGraph edge2) then False
       else if not (meetsAllCoevalConstraints constraintList edge1 edge2) then False
       else if (isAncDescEdge inGraph edge1 edge2) then False
       else True

-- | getCoevalConstraintEdges takes a graph and a network node and creates two lists: one of edges
-- "before" (ie towards root) and a second "after (ie away from root) 
-- this defines a coeval constraint.  No network edge can be added that would be directed
-- from the before group to the after
getCoevalConstraintEdges :: (Eq a, Eq b, Show a) => LG.Gr a b -> LG.LNode a -> ([LG.LEdge b],[LG.LEdge b])
getCoevalConstraintEdges inGraph inNode =
   if LG.isEmpty inGraph then error "Empty input graph in getCoevalConstraintEdges"
   else 
       let (_, edgeBeforeList) = LG.nodesAndEdgesBefore inGraph [inNode]
           (_, edgeAfterList) = LG.nodesAndEdgesAfter inGraph [inNode]
       in
       (edgeBeforeList, edgeAfterList)


-- | getGraphCoevalConstraints takes a greaph and returns coeval constraints based on network nodes
getGraphCoevalConstraints :: (Eq a, Eq b, Show a) => LG.Gr a b -> [([LG.LEdge b],[LG.LEdge b])]
getGraphCoevalConstraints inGraph =
   if LG.isEmpty inGraph then error "Empty input graph in getGraphCoevalConstraints"
   else 
       let (_, _, _, networkNodeList) = LG.splitVertexList inGraph
       in
       if null networkNodeList then []
       else fmap (getCoevalConstraintEdges inGraph) networkNodeList

-- | meetsAllCoevalConstraints checks constraint pair list and examines
-- whether one edge is fomr before and one after--if so fails False
-- else True if all pass
meetsAllCoevalConstraints :: (Eq b) =>[([LG.LEdge b],[LG.LEdge b])] -> LG.LEdge b -> LG.LEdge b -> Bool
meetsAllCoevalConstraints constraintList edge1 edge2 = 
   if null constraintList then True
   else 
       let (beforeList, afterList) = head constraintList
       in
       if edge1 `elem` beforeList && edge2 `elem` afterList then False
       else if edge2 `elem` beforeList && edge1 `elem` afterList then False
       else meetsAllCoevalConstraints (tail constraintList) edge1 edge2

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

-- | insertNetEdge inserts an edge between two other edges, creating 2 new nodes and rediagnoses graph
-- contacts deletes 2 orginal edges and adds 2 nodes and 5 new edges
-- does not check any egde reasonable-ness properties
-- new edge directed from first to second edge
-- naive for now
insertNetEdge :: GlobalSettings -> ProcessedData -> PhylogeneticGraph -> (LG.LEdge b, LG.LEdge b) -> PhylogeneticGraph
insertNetEdge inGS inData inPhyloGraph ((u,v, _), (u',v', _)) =
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
           newPhyloGraph = if (graphType inGS == SoftWired) then T.multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex newSimple
                           else if (graphType inGS == HardWired) then T.multiTraverseFullyLabelHardWired inGS inData leafGraph startVertex newSimple
                           else error "Unsupported graph type in deleteNetEdge.  Must be soft or hard wired"                   
       in
       newPhyloGraph

-- | deleteAllNetEdges deletes network edges one each each round until no better or additional 
-- graphs are found
-- call with ([], infinity) [single input graph]
deleteAllNetEdges :: GlobalSettings -> ProcessedData -> Int -> ([PhylogeneticGraph], VertexCost) -> [PhylogeneticGraph] -> [PhylogeneticGraph]
deleteAllNetEdges inGS inData numToKeep (curBestGraphList, curBestGraphCost) inPhyloGraphList =
   if null inPhyloGraphList then take numToKeep curBestGraphList
   else
      let currentCost = min curBestGraphCost (snd6 $ head inPhyloGraphList)
          (newGraphList, newGraphCost) = deleteEachNetEdge inGS inData numToKeep (head inPhyloGraphList)
      in
      -- worse graphs found--go on 
      if  newGraphCost > currentCost then deleteAllNetEdges inGS inData numToKeep ((head inPhyloGraphList) : curBestGraphList, currentCost) (tail inPhyloGraphList)

      -- "steepest style descent" abandons existing list if better cost found
      else if newGraphCost < currentCost then deleteAllNetEdges inGS inData numToKeep (newGraphList, newGraphCost) newGraphList 

      -- equal cost
      -- not sure if should add new graphs to queue to do edge deletion again
      else 
         -- new grapjh list contains the input graph if equal and filterd unique already in deleteEachNetEdge
         let newCurSameBestList = GO.getBVUniqPhylogeneticGraph True (curBestGraphList ++ newGraphList)
         in
         deleteAllNetEdges inGS inData numToKeep (newCurSameBestList, currentCost) (tail inPhyloGraphList)

-- | deleteEachNetEdge takes a phylogenetic graph and deletes all network edges one at time
-- and returns best list of new Phylogenetic Graphs and cost
-- even if worse--could be used for simulated annealing later
-- if equal returns unique graph list
deleteEachNetEdge :: GlobalSettings -> ProcessedData -> Int -> PhylogeneticGraph -> ([PhylogeneticGraph], VertexCost)
deleteEachNetEdge inGS inData numToKeep inPhyloGraph =
   if LG.isEmpty $ thd6 inPhyloGraph then error "Empty input phylogenetic graph in deleteAllNetEdges"
   else
      let currentCost = snd6 inPhyloGraph
          networkEdgeList = LG.netEdges $ thd6 inPhyloGraph
          newGraphList = fmap (deleteNetEdge inGS inData inPhyloGraph) networkEdgeList `using`  PU.myParListChunkRDS
          minCostGraphList = GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] newGraphList
          minCost = snd6 $ head minCostGraphList
      in
      -- no network edges to delete
      if null networkEdgeList then ([inPhyloGraph], currentCost)
      else if minCost /= currentCost then (minCostGraphList, minCost)
      else (GO.selectPhylogeneticGraph [("unique", (show numToKeep))] 0 ["unique"] $ inPhyloGraph : minCostGraphList, currentCost)



-- | deleteEdge deletes an edge (checking if network) and rediagnoses graph
-- contacts in=out=1 edgfes and removes node, reindexing nodes and edges
-- naive for now
deleteNetEdge :: GlobalSettings -> ProcessedData -> PhylogeneticGraph -> LG.Edge -> PhylogeneticGraph
deleteNetEdge inGS inData inPhyloGraph edgeToDelete =
   if LG.isEmpty $ thd6 inPhyloGraph then error "Empty input phylogenetic graph in deleteNetEdge"
   else if not (LG.isNetworkEdge (fst6 inPhyloGraph) edgeToDelete) then error ("Edge to delete: " ++ (show edgeToDelete) ++ " not in graph:\n" ++ (LG.prettify $ fst6 inPhyloGraph))
   else 
       let delSimple = LG.contractIn1Out1Edges $ LG.delEdge edgeToDelete $ fst6 inPhyloGraph
           leafGraph = LG.extractLeafGraph $ thd6 inPhyloGraph

           -- prune other edges if now unused
           pruneEdges = True

           -- don't warn that edges are being pruned
           warnPruneEdges = False

           -- graph optimization from root
           startVertex = Nothing

           -- full two-pass optimization
           newPhyloGraph = if (graphType inGS == SoftWired) then T.multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex delSimple
                           else if (graphType inGS == HardWired) then T.multiTraverseFullyLabelHardWired inGS inData leafGraph startVertex delSimple
                           else error "Unsupported graph type in deleteNetEdge.  Must be soft or hard wired"
       in
       newPhyloGraph

