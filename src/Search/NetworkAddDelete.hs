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

-- Need a "steepest" that takes first better netowrk for add and delete.
-- Choose order in max branch length, root-end, leaf-end, and at random


-- | moveAllNetEdges removes each edge and adds an edge to all possible plasses each round 
-- until no better or additional graphs are found
-- call with ([], infinity) [single input graph]
moveAllNetEdges :: GlobalSettings -> ProcessedData -> Int -> Int -> ([PhylogeneticGraph], VertexCost) -> [PhylogeneticGraph] -> ([PhylogeneticGraph], Int)
moveAllNetEdges inGS inData numToKeep counter (curBestGraphList, curBestGraphCost) inPhyloGraphList =
   if null inPhyloGraphList then (take numToKeep curBestGraphList, counter)
   else 
      let firstPhyloGraph = head inPhyloGraphList
          currentCost = min curBestGraphCost (snd6 firstPhyloGraph)
          netEdgeList = LG.labNetEdges (thd6 firstPhyloGraph)
          newGraphList' = concatMap (deleteOneNetAddAll inGS inData numToKeep currentCost firstPhyloGraph) netEdgeList
          newGraphList = GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] newGraphList'
          newGraphCost = snd6 $ head newGraphList
      in
      if null netEdgeList then trace ("\tNo network edges to move") ([], counter)
      else if  newGraphCost > currentCost then moveAllNetEdges inGS inData numToKeep (counter + 1) ((head inPhyloGraphList) : curBestGraphList, currentCost) (tail inPhyloGraphList)
      else if newGraphCost < currentCost then
         -- trace ("Found better move")
         moveAllNetEdges inGS inData numToKeep (counter + 1) (newGraphList, newGraphCost) (newGraphList ++ (tail inPhyloGraphList))

      else 
         -- new grapjh list contains the input graph if equal and filterd unique already in deleteEachNetEdge
         let newCurSameBestList = GO.getBVUniqPhylogeneticGraph True (curBestGraphList ++ newGraphList)
         in
         moveAllNetEdges inGS inData numToKeep (counter + 1) (newCurSameBestList, currentCost) (tail inPhyloGraphList)

      
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
      if minNewCost <= currentCost then GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] insertedGraphList
      else [] 



-- | insertAllNetEdges adds network edges one each each round until no better or additional 
-- graphs are found
-- call with ([], infinity) [single input graph]
insertAllNetEdges :: GlobalSettings -> ProcessedData -> Int -> Int -> ([PhylogeneticGraph], VertexCost) -> [PhylogeneticGraph] -> ([PhylogeneticGraph], Int)
insertAllNetEdges inGS inData numToKeep counter (curBestGraphList, curBestGraphCost) inPhyloGraphList =
   if null inPhyloGraphList then (take numToKeep curBestGraphList, counter)
   else
      let currentCost = min curBestGraphCost (snd6 $ head inPhyloGraphList)
          (newGraphList, newGraphCost) = insertEachNetEdge inGS inData numToKeep Nothing (head inPhyloGraphList)
      in
      -- worse graphs found--go on 
      if  newGraphCost > currentCost then insertAllNetEdges inGS inData numToKeep (counter + 1) ((head inPhyloGraphList) : curBestGraphList, currentCost) (tail inPhyloGraphList) 

      -- "steepest style descent" abandons existing list if better cost found
      else if newGraphCost < currentCost then 
         -- trace ("Found a new edge to insert")
         insertAllNetEdges inGS inData numToKeep (counter + 1) (newGraphList, newGraphCost) (newGraphList ++ (tail inPhyloGraphList))

      -- equal cost
      -- not sure if should add new graphs to queue to do edge deletion again
      else 
         -- new grapjh list contains the input graph if equal and filterd unique already in deleteEachNetEdge
         let newCurSameBestList = GO.getBVUniqPhylogeneticGraph True (curBestGraphList ++ newGraphList)
         in
         insertAllNetEdges inGS inData numToKeep (counter + 1) (newCurSameBestList, currentCost) (tail inPhyloGraphList)

-- | insertEachNetEdge takes a phylogenetic graph and inserts all permissible network edges one at time
-- and returns best list of new Phylogenetic Graphs and cost
-- even if worse--could be used for simulated annealing later
-- if equal returns unique graph list
insertEachNetEdge :: GlobalSettings -> ProcessedData -> Int -> Maybe VertexCost -> PhylogeneticGraph -> ([PhylogeneticGraph], VertexCost)
insertEachNetEdge inGS inData numToKeep preDeleteCost inPhyloGraph =
   if LG.isEmpty $ thd6 inPhyloGraph then error "Empty input insertEachNetEdge graph in deleteAllNetEdges"
   else
      let currentCost = snd6 inPhyloGraph

          candidateNetworkEdgeList = getPermissibleEdgePairs (thd6 inPhyloGraph)

          -- newGraphList = concat (fmap (insertNetEdgeBothDirections inGS inData inPhyloGraph) candidateNetworkEdgeList `using`  PU.myParListChunkRDS)
          newGraphList = fmap (insertNetEdge inGS inData inPhyloGraph preDeleteCost) candidateNetworkEdgeList `using`  PU.myParListChunkRDS

          minCostGraphList = GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] newGraphList
          minCost = if null minCostGraphList then infinity
                    else snd6 $ head minCostGraphList
      in
      trace ("\tExamining " ++ (show $ length candidateNetworkEdgeList) ++ " candidate edge pairs") (
      -- no network edges to insert
      if null candidateNetworkEdgeList then ([inPhyloGraph], currentCost)
      else if minCost /= currentCost then 
         -- trace ("IENE: " ++ (show (minCost, currentCost)) ++ " from " ++ (show $ fmap snd6 newGraphList) ++ " yeilding " ++ (show $ fmap snd6  minCostGraphList))
         (minCostGraphList, minCost)
      else 
         (GO.selectPhylogeneticGraph [("unique", (show numToKeep))] 0 ["unique"] $ inPhyloGraph : minCostGraphList, currentCost)
      )

-- | getPermissibleEdgePairs takes a DecoratedGraph and returns the list of all pairs
-- of edges that can be joined by a network edge and meet all necessary conditions
getPermissibleEdgePairs :: DecoratedGraph -> [(LG.LEdge EdgeInfo, LG.LEdge EdgeInfo)]
getPermissibleEdgePairs inGraph = 
   if LG.isEmpty inGraph then error "Empty input graph in isEdgePairPermissible"
   else 
       let edgeList = LG.labEdges inGraph
           edgePairs = cartProd edgeList edgeList
           contraintList = getGraphCoevalConstraints inGraph
           edgeTestList = fmap (isEdgePairPermissible inGraph contraintList) edgePairs `using`  PU.myParListChunkRDS
           pairList = filter ((== True) . snd) $ zip edgePairs edgeTestList
       in
       -- trace ("CArtProd:" ++ (show $ cartProd (fmap LG.toEdge edgeList) (fmap LG.toEdge edgeList)))
       fmap fst pairList

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
getGraphCoevalConstraints :: (Eq a, Eq b, Show a, NFData b) => LG.Gr a b -> [([LG.LEdge b],[LG.LEdge b])]
getGraphCoevalConstraints inGraph =
   if LG.isEmpty inGraph then error "Empty input graph in getGraphCoevalConstraints"
   else 
       let (_, _, _, networkNodeList) = LG.splitVertexList inGraph
       in
       if null networkNodeList then []
       else fmap (getCoevalConstraintEdges inGraph) networkNodeList `using`  PU.myParListChunkRDS

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

-- | insertNetEdgeBothDirections calls insertNetEdge for both u -> v and v -> u new edge orientations
insertNetEdgeBothDirections :: GlobalSettings -> ProcessedData -> PhylogeneticGraph ->  Maybe VertexCost -> (LG.LEdge b, LG.LEdge b) -> [PhylogeneticGraph]
insertNetEdgeBothDirections inGS inData inPhyloGraph preDeleteCost (u,v) = fmap (insertNetEdge inGS inData inPhyloGraph preDeleteCost) [(u,v), (v,u)]

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
           newPhyloGraph = if (graphType inGS == SoftWired) then 
                                -- trace ("NewSimple\n:" ++ (LG.prettify newSimple)) 
                                T.multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex newSimple
                           else if (graphType inGS == HardWired) then T.multiTraverseFullyLabelHardWired inGS inData leafGraph startVertex newSimple
                           else error "Unsupported graph type in deleteNetEdge.  Must be soft or hard wired" 

            
           -- calculates heursitic graph delta
           (heuristicDelta, _, _, _, _)  = heuristicAddDelta inGS inPhyloGraph edgePair (fst newNodeOne) (fst newNodeTwo) 

           edgeAddDelta = deltaPenaltyAdjustment (graphFactor inGS) (V.length $ fst3 inData) inPhyloGraph


       in
       -- trace ("INE Deltas: " ++ (show (heuristicDelta, edgeAddDelta)) ++ " preDelete " ++ (show preDeleteCost)
       --  ++ "New Nodes " ++ (show [newNodeOne, newNodeTwo]) ++ " delete edges " ++ (show [(u,v), (u',v')]) ++ " New edges " ++ (show newEdgeList)
       --  ++ "\nInGraph:\n" ++ (LG.prettify inSimple) ++ "\nNewGraph:\n" ++ (LG.prettify newSimple) ) (

       -- preDelete cost changes criterion for edge move
       if preDeleteCost == Nothing then 
          if heuristicDelta + edgeAddDelta < 0 then newPhyloGraph
          else emptyPhylogeneticGraph
          
       else 
         -- no net add cost becasue the numbe rof net nodes is unchaned in add/delete when preDelete cost /= Noting
          if heuristicDelta + (snd6 inPhyloGraph) <= fromJust preDeleteCost then newPhyloGraph
          else emptyPhylogeneticGraph
       --   )

-- | heuristicAddDelta takes teh existing graph, edge pair, and new nodes to create and makes
-- the new nodes and reoprtimizes starting nodes of two edges.  Returns cost delta based on 
-- previous and new node resolution caches
-- returns cost delta and the reoptimized nodes for use in incremental optimization
-- creates (v', n2), (n2, X)u' new, (v, n2)n1  (n1,Y) new  
-- results of new edge n1 -> n2, u' -> (n2, X), n1 -> (n2,v), u (n1,Y)
heuristicAddDelta :: GlobalSettings -> PhylogeneticGraph -> (LG.LEdge b, LG.LEdge b) -> LG.Node -> LG.Node -> (VertexCost, LG.LNode VertexInfo, LG.LNode VertexInfo, LG.LNode VertexInfo, LG.LNode VertexInfo)
heuristicAddDelta inGS inPhyloGraph ((u,v, _), (u',v', _)) n1 n2 =
  if LG.isEmpty (fst6 inPhyloGraph) then error "Empty graph in heuristicAddDelta"
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
      -- this should not happen--should try to crete new edges from children of net edges
      if (length $ LG.descendants (thd6 inPhyloGraph) u) < 2 ||  (length $ LG.descendants (thd6 inPhyloGraph) u') < 2 then error ("Outdegree 1 nodes in heuristicAddDelta")
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
   else error ("Netowrk edge cost model not yet implemented: " ++ (show edgeCostModel))


-- | deleteAllNetEdges deletes network edges one each each round until no better or additional 
-- graphs are found
-- call with ([], infinity) [single input graph]
deleteAllNetEdges :: GlobalSettings -> ProcessedData -> Int -> Int -> ([PhylogeneticGraph], VertexCost) -> [PhylogeneticGraph] -> ([PhylogeneticGraph], Int)
deleteAllNetEdges inGS inData numToKeep counter (curBestGraphList, curBestGraphCost) inPhyloGraphList =
   if null inPhyloGraphList then (take numToKeep curBestGraphList, counter)
   else
      let currentCost = min curBestGraphCost (snd6 $ head inPhyloGraphList)
          (newGraphList, newGraphCost) = deleteEachNetEdge inGS inData numToKeep False (head inPhyloGraphList)
      in
      -- worse graphs found--go on 
      if  newGraphCost > currentCost then deleteAllNetEdges inGS inData numToKeep (counter + 1) ((head inPhyloGraphList) : curBestGraphList, currentCost) (tail inPhyloGraphList) 

      -- "steepest style descent" abandons existing list if better cost found
      else if newGraphCost < currentCost then deleteAllNetEdges inGS inData numToKeep (counter + 1) (newGraphList, newGraphCost) (newGraphList ++ (tail inPhyloGraphList))

      -- equal cost
      -- not sure if should add new graphs to queue to do edge deletion again
      else 
         -- new grapjh list contains the input graph if equal and filterd unique already in deleteEachNetEdge
         let newCurSameBestList = GO.getBVUniqPhylogeneticGraph True (curBestGraphList ++ newGraphList)
         in
         deleteAllNetEdges inGS inData numToKeep  (counter + 1)  (newCurSameBestList, currentCost) (tail inPhyloGraphList)

-- | deleteEachNetEdge takes a phylogenetic graph and deletes all network edges one at time
-- and returns best list of new Phylogenetic Graphs and cost
-- even if worse--could be used for simulated annealing later
-- if equal returns unique graph list
deleteEachNetEdge :: GlobalSettings -> ProcessedData -> Int -> Bool ->  PhylogeneticGraph -> ([PhylogeneticGraph], VertexCost)
deleteEachNetEdge inGS inData numToKeep force inPhyloGraph =
   if LG.isEmpty $ thd6 inPhyloGraph then error "Empty input phylogenetic graph in deleteAllNetEdges"
   else
      let currentCost = snd6 inPhyloGraph
          networkEdgeList = LG.netEdges $ thd6 inPhyloGraph
          newGraphList = fmap (deleteNetEdge inGS inData inPhyloGraph force) networkEdgeList `using`  PU.myParListChunkRDS
          minCostGraphList = GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] newGraphList
          minCost = if null minCostGraphList then infinity 
                    else snd6 $ head minCostGraphList
      in
      -- no network edges to delete
      if null networkEdgeList then trace ("\tNo network edges to delete") ([inPhyloGraph], currentCost)
      else if minCost /= currentCost then (minCostGraphList, minCost)
      else (GO.selectPhylogeneticGraph [("unique", (show numToKeep))] 0 ["unique"] $ inPhyloGraph : minCostGraphList, currentCost)



-- | deleteEdge deletes an edge (checking if network) and rediagnoses graph
-- contacts in=out=1 edgfes and removes node, reindexing nodes and edges
-- naive for now
-- force requires reoptimization no matter what--used for net move
deleteNetEdge :: GlobalSettings -> ProcessedData -> PhylogeneticGraph -> Bool -> LG.Edge -> PhylogeneticGraph
deleteNetEdge inGS inData inPhyloGraph force edgeToDelete =
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

           (heuristicDelta, _, _) = heuristicDeleteDelta inGS inPhyloGraph edgeToDelete

           edgeAddDelta = deltaPenaltyAdjustment (graphFactor inGS) (V.length $ fst3 inData) inPhyloGraph
           
           -- full two-pass optimization
           newPhyloGraph = if (graphType inGS == SoftWired) then T.multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex delSimple
                           else if (graphType inGS == HardWired) then T.multiTraverseFullyLabelHardWired inGS inData leafGraph startVertex delSimple
                           else error "Unsupported graph type in deleteNetEdge.  Must be soft or hard wired"
       in
       if force then newPhyloGraph
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
      if (length (LG.parents inGraph n1) /= 1) || (length (LG.parents inGraph n2) /= 2) || (length (LG.descendants inGraph n2) /= 1) || (length (LG.descendants inGraph n1) /= 2) then error ("Graph malformation in numbersof parents and children in heuristicDeleteDelta")
      else 
         (addNetDelta, (u, uLabAfter), (u', uPrimeLabAfter))
           