{- |
Module      :  GraphOperations.hs
Description :  Module specifying general graph functions
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

module Graphs.GraphOperations (  ladderizeGraph
                               , rerootTree
                               , rerootDisplayTree
                               , generateDisplayTrees
                               , getNodeType
                               , convertDecoratedToSimpleGraph
                               , convertToSimpleEdge
                               , graphCostFromNodes
                               , switchRootTree
                               , dichotomizeRoot
                               , showDecGraphs
                               , sortEdgeListByLength
                               , sortEdgeListByDistance
                               , splitGraphOnEdge
                               , joinGraphOnEdge
                               , splitGraphOnEdge'
                               , getEdgeSplitList
                               , selectPhylogeneticGraph
                               , getUniqueGraphs
                               , copyIAFinalToPrelim
                               , copyIAPrelimToFinal
                               , makeIAFinalFromPrelim
                               , makeIAPrelimFromFinal
                               , topologicalEqual
                               , getTopoUniqPhylogeneticGraph
                               , getBVUniqPhylogeneticGraph
                               , makeDummyLabEdge
                               , contractIn1Out1EdgesRename
                               , renameSimpleGraphNodes
                               , renameSimpleGraphNodesString
                               , generateDisplayTreesRandom
                               , hasNetNodeAncestorViolation
                               , convertGeneralGraphToPhylogeneticGraph
                               , parentInChain
                               , selectGraphStochastic
                               , makeNewickList
                               ) where

import           Bio.DynamicCharacter
import           Control.Parallel.Strategies
import           Data.Bits
import qualified Data.BitVector.LittleEndian as BV
import qualified Data.Char                   as C
import qualified Data.List                   as L
import           Data.Maybe
import qualified Data.Text.Lazy              as T
import qualified Data.Vector                 as V
import qualified Data.Vector.Generic         as GV
import           Debug.Trace
import           GeneralUtilities
import qualified GraphFormatUtilities        as GFU
import qualified GraphOptimization.Medians   as M
import qualified ParallelUtilities           as PU
import           Text.Read
import           Types.Types
import qualified Utilities.LocalGraph        as LG
import qualified Utilities.Utilities         as U

-- | makeNewickList takes a list of fgl trees and outputs a single String cointaining the graphs in Newick format
makeNewickList ::  Bool -> Bool -> Int -> [SimpleGraph] -> [VertexCost] -> String
makeNewickList writeEdgeWeight writeNodeLabel' rootIndex graphList costList =
    let allTrees = L.foldl' (&&) True (fmap LG.isTree graphList)

        -- check for network HTU label requirement
        writeNodeLabel = if allTrees then writeNodeLabel'
                         else if writeNodeLabel' then writeNodeLabel'
                         else 
                            trace ("HTU labels are required for ENewick Output")
                            True

        graphString = GFU.fglList2ForestEnhancedNewickString (fmap (rerootTree rootIndex) graphList)  writeEdgeWeight writeNodeLabel
        newickStringList = fmap init $ filter (not . null) $ lines graphString
        costStringList  = fmap (('[' :) . (++ "];\n")) (fmap show costList)
        graphStringCost = concat $ zipWith (++) newickStringList costStringList
    in
    graphStringCost

-- | convertGeneralGraphToPhylogeneticGraph inputs a SimpleGraph and converts it to a Phylogenetic graph by:
--  1) transitive reduction -- removes anc <-> desc netork edges
--  2) ladderizes -- all vertices are (in degree, outdegree) (0,1|2) (1,2) (2,1) (1,0)
--        by adding extra HTIs and edges
--  3) checks time consistency and removes edges stepwise from
--        those that violate the most  before/after splits of network edges
--        arbitrary but deterministic
--  4) contracts out any remaning indegree 1 outdegree 1 nodes and renames HTUs in order
convertGeneralGraphToPhylogeneticGraph :: SimpleGraph -> SimpleGraph
convertGeneralGraphToPhylogeneticGraph inGraph =
  if LG.isEmpty inGraph then LG.empty
  else
    let -- remove single "tail" edge from root with single child, replace child node with root
        noTailGraph = LG.contractRootOut1Edge inGraph

        -- remove indeg 1 out deg 1 edges
        noIn1Out1Graph = LG.contractIn1Out1Edges noTailGraph

        -- transitive reduction
        reducedGraph = LG.transitiveReduceGraph noIn1Out1Graph

        -- laderization of indegree and outdegree edges
        ladderGraph = ladderizeGraph reducedGraph

        -- time consistency (after those removed by transitrive reduction)
        timeConsistentGraph = makeGraphTimeConsistent ladderGraph

        -- removes ancestor descendent edges
        noParentChainGraph = removeParentsInChain timeConsistentGraph

        -- remove sister-sister edge.  where two network nodes have same parents
        noSisterSisterGraph = removeSisterSisterEdges noParentChainGraph

    in
    -- trace ("CGP orig:\n" ++ (LG.prettify inGraph) ++ "\nNew:" ++ (LG.prettify timeConsistentGraph))
    -- cycle check to make sure
    if LG.cyclic noSisterSisterGraph then error ("Cycle in graph : \n" ++ (LG.prettify noSisterSisterGraph))
    else noSisterSisterGraph

-- | parentsInChain checks for parents in chain ie network edges
-- that implies a network event between nodes where one is the ancestor of the other
-- a time violation
parentInChain :: (Show a, Eq a, Eq b) => LG.Gr a b -> Bool
parentInChain inGraph =
  if LG.isEmpty inGraph then error "Null graph in parentInChain"
  else
    let   (_, _, _, netVertexList) = LG.splitVertexList inGraph
          parentNetVertList = fmap (LG.labParents inGraph) $ fmap fst netVertexList

          -- get list of nodes that are transitively equal in age
          concurrentList = mergeConcurrentNodeLists parentNetVertList []
          concurrentPairList = concatMap getListPairs concurrentList

          -- get pairs that violate concurrency
          violatingConcurrentPairs = concatMap (concurrentViolatePair inGraph) concurrentPairList
    in
    if null violatingConcurrentPairs then False
    else True

-- | removeParentsInChain checks the parents of each netowrk node are not anc/desc of each other
removeParentsInChain :: SimpleGraph -> SimpleGraph
removeParentsInChain inGraph =
  if LG.isEmpty inGraph then LG.empty
  else
      let (_, _, _, netVertexList) = LG.splitVertexList inGraph
          parentNetVertList = fmap (LG.labParents inGraph) $ fmap fst netVertexList

          -- get list of nodes that are transitively equal in age
          concurrentList = mergeConcurrentNodeLists parentNetVertList []
          concurrentPairList = concatMap getListPairs concurrentList

          -- get pairs that violate concurrency
          violatingConcurrentPairs = concatMap (concurrentViolatePair inGraph) concurrentPairList

          -- get netowrk nodes with violations
          parentNodeViolateList = concatMap pairToList violatingConcurrentPairs
          childNodeViolateList = concatMap (LG.descendants inGraph) parentNodeViolateList
          netNodeViolateList = filter (LG.isNetworkNode inGraph) childNodeViolateList

          netEdgesThatViolate = fmap LG.toEdge $ LG.inn inGraph $ head netNodeViolateList

      in
      if null violatingConcurrentPairs then inGraph
      else if null netNodeViolateList then error ("Should be neNode that violate")
      else if null netEdgesThatViolate then error "Should be violating in edges"
      else
        let edgeDeletedGraph = LG.delEdge (head netEdgesThatViolate) inGraph
            newGraph = contractIn1Out1EdgesRename edgeDeletedGraph
        in
        -- trace ("PIC")
        removeParentsInChain newGraph
    where pairToList (a,b) = [fst a, fst b]

-- | removeSisterSisterEdges takes a graph and recursively removes a single edge fomr where two network
-- edges have the same two parents
removeSisterSisterEdges :: SimpleGraph -> SimpleGraph
removeSisterSisterEdges inGraph =
  if LG.isEmpty inGraph then LG.empty
  else
    let sisterSisterEdges = getSisterSisterEdgeList inGraph
        newGraph = LG.delEdge (head sisterSisterEdges) inGraph
        newGraph' = contractIn1Out1EdgesRename newGraph
    in
    if null sisterSisterEdges then inGraph
    else
      -- trace ("Sister")
      removeSisterSisterEdges  newGraph'



-- | getSisterSisterEdgeList take a graph and returns list of edges where two network nodes
-- have the same two parents
getSisterSisterEdgeList :: LG.Gr a b -> [LG.Edge]
getSisterSisterEdgeList inGraph =
  if LG.isEmpty inGraph then []
  else
      let (_, _, _, netVertexList) = LG.splitVertexList inGraph
          netVertPairs = getListPairs $ fmap fst netVertexList

          parentsList = fmap (LG.parents inGraph) (fmap fst netVertexList)
          parentPairs = getListPairs $ parentsList

          netAndParentPairs = zip netVertPairs parentPairs
          sisterSisterPairs = filter sameParents netAndParentPairs
          sisterSisterEdges = concatMap makeChildParentEdges sisterSisterPairs
      in
      if null sisterSisterPairs then []
      else sisterSisterEdges
      where sameParents (_, (b, c)) = if b == c then True else False
            makeChildParentEdges ((a1,a2), (b, _)) = [(head b,a1), (last b,a1), (head b,a2), (last b,a2)]


-- | concurrentViolatePair takes a pair of nodes and sees if either is ancetral to the other--if so returns pair
-- as list otherwise null list
concurrentViolatePair :: (Eq a, Show a, Eq b) => LG.Gr a b  -> (LG.LNode a, LG.LNode a) -> [(LG.LNode a, LG.LNode a)]
concurrentViolatePair inGraph (node1, node2) =
  if LG.isEmpty inGraph then error "Empty graph in concurrentViolatePair"
  else
    let (nodesBeforeFirst, _)  = LG.nodesAndEdgesBefore inGraph [node1]
        (nodesBeforeSecond, _) = LG.nodesAndEdgesBefore inGraph [node2]
    in
    if node2 `elem` nodesBeforeFirst then [(node1, node2)]
    else if node1 `elem` nodesBeforeSecond then [(node1, node2)]
    else []


-- | mergeConcurrentNodeLists takes a list os lists  and returns a list os lists of merged lists
-- lists are merged if they share any elements
mergeConcurrentNodeLists :: (Eq a) => [[LG.LNode a]] -> [[LG.LNode a]] -> [[LG.LNode a]]
mergeConcurrentNodeLists inListList currentListList =
  if null inListList then
    -- trace ("MCNL:" ++ (show $ fmap (fmap fst) currentListList))
    currentListList

  -- first case
  else if null currentListList then mergeConcurrentNodeLists (tail inListList) [head inListList]
  else
    let firstList = head inListList
        (intersectList, _) = unzip $ filter ((== True) . snd) $ zip (currentListList) (fmap (not . null) $ fmap (L.intersect firstList) currentListList)

        (noIntersectLists, _) = unzip $ filter ((== False) . snd) $ zip (currentListList) (fmap (not . null) $ fmap (L.intersect firstList) currentListList)

        mergedList = if null intersectList then firstList
                     else L.foldl' L.union firstList intersectList
    in
    -- trace ("MCL-F:" ++ (show $ fmap fst firstList) ++ " inter " ++ (show $ fmap (fmap fst) intersectList) ++
    --   " noInter " ++ (show $ fmap (fmap fst) noIntersectLists) ++ " curList " ++ (show $ fmap (fmap fst) currentListList))
    mergeConcurrentNodeLists (tail inListList) (mergedList : noIntersectLists)

{-
-- | checkParentsChain takes a graph vertexNode and its parents and checks if one parent is descnedent of the other
-- a form of time violation
checkParentsChain :: (Show a, Eq a, Eq b) => LG.Gr a b -> LG.LNode a -> [LG.LNode a] -> [LG.Edge]
checkParentsChain inGraph netNode parentNodeList =
  if LG.isEmpty inGraph then error "Empty graph in checkParentsChain"
  else if length parentNodeList /= 2 then error ("Need to have 2 parents for net node: " ++ (show (fst netNode)) ++ " <- " ++ (show $ fmap fst parentNodeList))
  else
    let firstParent = head parentNodeList
        secondParent = last parentNodeList
        (nodesBeforeFirst, _)  = LG.nodesAndEdgesBefore inGraph [firstParent]
        (nodesBeforeSecond, _) = LG.nodesAndEdgesBefore inGraph [secondParent]
    in
    -- trace ("CPC:" ++ (show $ fst netNode) ++ " <- " ++ (show $ fmap fst parentNodeList) ++ "\nfirstParentBefore: "
    --  ++ (show $ fmap fst nodesBeforeFirst) ++ "\nsecondParentBefore: " ++ (show $ fmap fst nodesBeforeSecond)) (
    if secondParent `elem` nodesBeforeFirst then [(fst firstParent, fst netNode)]
    else if firstParent `elem` nodesBeforeSecond then [(fst secondParent, fst netNode)]
    else []
    -- )
-}

-- | makeGraphTimeConsistent takes laderized, trasitive reduced graph and deletes
-- network edges in an arbitrary but deterministic sequence to produce a phylogentic graphs suitable
-- for swapping etc
-- recursively removes the `most' inconsistent edge (breaks the highest number of time consistent conditions)
-- contracts and renames edges at each stage
-- O (n *m) n nodes, m network nodes.  Probbably could be done more efficiently by retainliung beforeafter list and not
-- remaking graph each time.
makeGraphTimeConsistent :: SimpleGraph -> SimpleGraph
makeGraphTimeConsistent inGraph =
  if LG.isEmpty inGraph then LG.empty
  else
    let coevalNodeConstraintList = LG.getGraphCoevalConstraintsNodes inGraph
        networkEdgeList = concatMap (LG.inn inGraph) $ fmap fst $ fmap fst3 coevalNodeConstraintList
        networkEdgeViolationList = fmap (numberTimeViolations coevalNodeConstraintList 0) networkEdgeList
        maxViolations = maximum networkEdgeViolationList
        edgeMaxViolations = fst $ head $ filter ((== maxViolations) . snd) $ zip networkEdgeList networkEdgeViolationList
    in
    -- is a tree
    if null coevalNodeConstraintList then
      -- trace ("MTC Null:\n" ++ (LG.prettify inGraph))
      inGraph

    -- is time consistent
    else if maxViolations == 0 then
      -- trace ("MTC 0:\n" ++ (LG.prettify inGraph))
      removeParentsInChain inGraph

    -- has time violations
    else
      let edgeDeletedGraph = LG.delLEdge edgeMaxViolations inGraph
          newGraph = contractIn1Out1EdgesRename edgeDeletedGraph
      in
      -- trace ("MTC V:" ++ (show maxViolations))
      makeGraphTimeConsistent newGraph

-- | numberTimeViolations takes a directed edge (u,v) and pairs of before after edge lists
-- if u is in the before list and v in the after list of a pir--then there is a time violation
-- recursively counts and retunns the number of violations
numberTimeViolations :: [(LG.LNode a, [LG.LEdge b], [LG.LEdge b])] -> Int -> LG.LEdge b -> Int
numberTimeViolations inTripleList counter inEdge@(u,v,_) =
  if null inTripleList then counter
  else
    let ((tripleNode, _), beforeEdgeList, afterEdgeList) = head inTripleList
        uInBefore = u `elem`  ((fmap fst3 beforeEdgeList) ++ (fmap snd3 beforeEdgeList))
        vInAfter  = v `elem`  ((fmap fst3 afterEdgeList) ++ (fmap snd3 afterEdgeList))
    in
    --trace ("NTV: " ++ (show (u,v)) ++ " node " ++ (show tripleNode) ++ "\nbefore: " ++ (show $ ((fmap fst3 beforeEdgeList) ++ (fmap snd3 beforeEdgeList))) ++ "\nafter: " ++ (show $ ((fmap fst3 afterEdgeList) ++ (fmap snd3 afterEdgeList)))) (

    -- skipping its own split
    if v == tripleNode || u == tripleNode then numberTimeViolations (tail inTripleList) counter inEdge

    -- violates this time pair
    else if uInBefore && vInAfter then numberTimeViolations (tail inTripleList) (counter + 1) inEdge

    -- does not violate pair
    else numberTimeViolations (tail inTripleList) counter inEdge
    -- )

-- | contractIn1Out1EdgesRename contracts in degree and outdegree edges and renames HTUs in index order
-- does one at a time and makes a graph and recurses
contractIn1Out1EdgesRename :: SimpleGraph -> SimpleGraph
contractIn1Out1EdgesRename inGraph =
  if LG.isEmpty inGraph then LG.empty
  else
    let newGraph = LG.contractIn1Out1Edges inGraph
    in
    renameSimpleGraphNodes newGraph


-- | renameSimpleGraphNodes takes nodes and renames HTU nodes based on index
renameSimpleGraphNodes :: SimpleGraph -> SimpleGraph
renameSimpleGraphNodes inGraph =
  if LG.isEmpty inGraph then LG.empty
    else
      let inNodes = LG.labNodes inGraph
          nodeLabels = fmap (makeSimpleLabel inGraph) inNodes
          newNodes = zip (fmap fst inNodes) nodeLabels
          newEdges = LG.labEdges inGraph
    in
    --newGraph
    -- trace ("C11: " ++ (show $ LG.getIsolatedNodes newGraph) ++ " => " ++ (show newNodes) ++ " " ++ (show $ fmap LG.toEdge newEdges))
    LG.mkGraph newNodes newEdges
    where makeSimpleLabel g (a, b)  = if (not $ LG.isLeaf g a) then T.pack $ "HTU"  ++ show a
                                      else b

-- | renameSimpleGraphNodesString takes nodes and renames HTU nodes based on index
renameSimpleGraphNodesString :: LG.Gr String String -> LG.Gr String String
renameSimpleGraphNodesString inGraph =
  if LG.isEmpty inGraph then LG.empty
    else
      let inNodes = LG.labNodes inGraph
          nodeLabels = fmap (makeSimpleLabel inGraph) inNodes
          newNodes = zip (fmap fst inNodes) nodeLabels
          newEdges = LG.labEdges inGraph
    in
    --newGraph
    -- trace ("C11: " ++ (show $ LG.getIsolatedNodes newGraph) ++ " => " ++ (show newNodes) ++ " " ++ (show $ fmap LG.toEdge newEdges))
    LG.mkGraph newNodes newEdges
    where makeSimpleLabel g (a, b)  = if (not $ LG.isLeaf g a) then "HTU"  ++ show a
                                    else b

-- | getEdgeSplitList takes a graph and returns list of edges
-- that split a graph increasing the number of components by 1
-- this is quadratic
-- should change to Tarjan's algorithm (linear)
-- everyhting else in there is O(n^2-3) so maybe doesn't matter
-- filters out edges with parent nodes that are out degree 1 and root edges
getEdgeSplitList :: (Show a, Show b, Eq b) => LG.Gr a b -> [LG.LEdge b]
getEdgeSplitList inGraph =
  if LG.isEmpty inGraph then error ("Empty graph in getEdgeSplitList")
  else
      let origNumComponents = LG.noComponents inGraph
          origEdgeList = LG.labEdges inGraph
          edgeDeleteComponentNumberList = fmap LG.noComponents $ fmap (flip LG.delEdge inGraph) (fmap LG.toEdge origEdgeList)
          bridgeList = fmap snd $ filter ((> origNumComponents) . fst) $ zip edgeDeleteComponentNumberList origEdgeList

          -- filter out edges starting in an outdegree 1 node (network or in out 1) node
          -- this would promote an HTU to a leaf later.  Its a bridge, but not what need
          bridgeList' = filter ((not . LG.isRoot inGraph) .fst3 ) $ filter ((not. LG.isNetworkNode inGraph) . snd3)  $ filter  ((not . LG.isOutDeg1Node inGraph) . fst3) bridgeList
      in

       -- trace ("BridgeList" ++ (show $ fmap LG.toEdge bridgeList') ++ "\nGraph\n" ++ (LG.prettyIndices inGraph))
       bridgeList'

-- | splitGraphOnEdge takes a graph and an edge and returns a single graph but with two components
-- the roots of each component are retuned with two graphs, with broken edge contraced, and 'naked'
-- node returned.  The naked node is used for rejoining the two components during rearrangement
-- (SplitGraph, root of component that has original root, root of component that was cut off, naked node left over)
-- this function does not check whether edge is a 'bridge'
splitGraphOnEdge :: (Show b) => LG.Gr a b -> LG.LEdge b -> (LG.Gr a b, LG.Node, LG.Node, LG.Node)
splitGraphOnEdge inGraph (e,v,l) =
  if LG.isEmpty inGraph then error "Empty graph in splitGraphOnEdge"
  else if (length $ LG.getRoots inGraph) /= 1 then error ("Incorrect number roots in splitGraphOnEdge--must be 1: " ++ (show $ fmap fst $ LG.getRoots inGraph))
  else
      let childrenENode = (LG.descendants inGraph e) L.\\ [v]
          parentsENode = LG.parents inGraph e
          newEdge = (head parentsENode, head childrenENode, l)
          edgesToDelete = [(head parentsENode, e), (e, head childrenENode)] -- (e,v)

          -- make new graph
          splitGraph = LG.insEdge newEdge $ LG.delEdges edgesToDelete inGraph
      in
      if length childrenENode /= 1 then error ("Incorrect number of children of edge to split--must be 1: " ++ (show ((e,v), childrenENode)))
      else if length parentsENode /= 1 then error ("Incorrect number of parents of edge to split--must be 1: " ++ (show ((e,v), parentsENode)))
      else
          -- trace ("SGE:" ++ (show (childrenENode, parentsENode, newEdge, edgesToDelete)))
          (splitGraph, fst $ head $ LG.getRoots inGraph, v, e)

-- | splitGraphOnEdge' like splitGrpahOnEdge above but returns edges creted and destroyed as well
-- used in Goodman-Bermer and could make swap more efficient as well.
splitGraphOnEdge' :: (Show b) => LG.Gr a b -> LG.LEdge b -> (LG.Gr a b, LG.Node, LG.Node, LG.Node, LG.LEdge b, [LG.Edge])
splitGraphOnEdge' inGraph (e,v,l) =
  if LG.isEmpty inGraph then error "Empty graph in splitGraphOnEdge"
  else if (length $ LG.getRoots inGraph) /= 1 then error ("Incorrect number roots in splitGraphOnEdge--must be 1: " ++ (show $ fmap fst $ LG.getRoots inGraph))
  else
      let childrenENode = (LG.descendants inGraph e) L.\\ [v]
          parentsENode = LG.parents inGraph e
          newEdge = (head parentsENode, head childrenENode, l)
          edgesToDelete = [(head parentsENode, e), (e, head childrenENode)] -- (e,v)

          -- make new graph
          splitGraph = LG.insEdge newEdge $ LG.delEdges edgesToDelete inGraph
      in
      if length childrenENode /= 1 then error ("Incorrect number of children of edge to split--must be 1: " ++ (show ((e,v), childrenENode)))
      else if length parentsENode /= 1 then error ("Incorrect number of parents of edge to split--must be 1: " ++ (show ((e,v), parentsENode)))
      else
          -- trace ("SGE:" ++ (show (childrenENode, parentsENode, newEdge, edgesToDelete)))
          (splitGraph, fst $ head $ LG.getRoots inGraph, v, e, newEdge, edgesToDelete)

-- | joinGraphOnEdge takes a graph and adds an edge reducing the component number
-- expected ot be two components to one in SPR/TBR
-- assumes that first node of edge (e,v,l) is 'naked' ie avaiable to make edges but is in graph
-- created from splitGraphOnEdge
joinGraphOnEdge :: (Show a,Show b) => LG.Gr a b -> LG.LEdge b -> LG.Node ->LG.Gr a b
joinGraphOnEdge inGraph (x,y,l) parentofPrunedSubGraph =
  if LG.isEmpty inGraph then error ("Empty graph in joinGraphOnEdge")
  else
      let edgeToCreate0 = (x, parentofPrunedSubGraph, l)
          edgeToCreate1 = (parentofPrunedSubGraph, y, l)
          -- edgeToCreate2 = (parentofPrunedSubGraph, graphToJoinRoot, l)
      in
      -- make new graph
      -- trace ("JGE:" ++ (show edgeToInvade) ++ " " ++ (show (parentofPrunedSubGraph, graphToJoinRoot))) --  ++ "\n" ++ (LG.prettify inGraph))
      LG.insEdges [edgeToCreate0, edgeToCreate1] $ LG.delEdge (x,y) inGraph

-- | sortEdgeListByLength sorts edge list by length (midRange), highest to lowest
sortEdgeListByLength :: [LG.LEdge EdgeInfo] -> [LG.LEdge EdgeInfo]
sortEdgeListByLength inEdgeList =
  if null inEdgeList then []
  else
    reverse $ L.sortOn (midRangeLength . thd3) inEdgeList

-- | sortEdgeListByDistance sorts edges by distance (in edges) from edge pair of vertices
-- cretes a list of edges into (but traveling away from) an initial eNOde and away from
-- an initial vNode adding new nodes to those lists as encountered by traversing edges.
-- the eidea is theat the nodes from a directed edge (eNode, vNode)
-- the list is creted at each round from the "in" and "out" edge lists
-- so they are in order of 1 edge 2 edges etc.
sortEdgeListByDistance :: LG.Gr a b -> [LG.Node] -> [LG.Node] -> [LG.LEdge b]
sortEdgeListByDistance inGraph eNodeList vNodeList =
  if LG.isEmpty inGraph then error ("Empty graph in edgeListByDistance")
  else if (null eNodeList && null vNodeList) then []
  else
        -- get edges 'in' to eNodeList
    let inEdgeList = concatMap (LG.inn inGraph) eNodeList
        newENodeList = fmap fst3 inEdgeList

        -- get edges 'out' from vNodeList
        outEdgeList = concatMap (LG.out inGraph) vNodeList
        newVNodeList = fmap snd3 outEdgeList
    in
    inEdgeList ++ outEdgeList ++ (sortEdgeListByDistance inGraph newENodeList newVNodeList)

-- | ladderizeGraph is a wrapper around ladderizeGraph' to allow for mapping with
-- local nodelist
ladderizeGraph :: SimpleGraph -> SimpleGraph
ladderizeGraph inGraph = ladderizeGraph' inGraph (LG.nodes inGraph)

-- | ladderize takes an input graph and ensures/creates nodes
-- such that all vertices are (indegree, outdegree) (0,>0), (1,2) (2,1) (1,0)
--ladderizeGraph' :: SimpleGraph -> [LG.Node] -> SimpleGraph
ladderizeGraph' :: SimpleGraph -> [LG.Node] -> SimpleGraph
ladderizeGraph' inGraph nodeList
  | LG.isEmpty inGraph = LG.empty
  | null nodeList = inGraph
  | otherwise =
  let -- these are roots, network, tree, leaf nodes
      okNodeDegrees = [(0,1),(0,2),(1,2),(2,1),(1,0)]
      firstNode = head nodeList
      (inEdgeList, outEdgeList) = LG.getInOutEdges inGraph firstNode
      inOutPairLength = (length inEdgeList, length outEdgeList)
  in
  --trace ("node " ++ (show firstNode) ++ " " ++ (show inOutPair)) (
  -- node ok to keep
  if inOutPairLength `elem` okNodeDegrees then ladderizeGraph' inGraph (tail nodeList)
  -- node edges need modification
  else
    let newGraph = resolveNode inGraph firstNode (inEdgeList, outEdgeList) inOutPairLength
    in
    ladderizeGraph' newGraph (LG.nodes newGraph)


-- | resolveNode takes a graph and node and inbound edgelist and outbound edge list
-- and converts node to one of (indeg, outdeg) (0,1),(0,2),(1,2),(2,1),(1,0)
-- this only resolves a single nodes edges at a time and then returns new graph
-- when more hase to be done--that will occur on lultiple passes through nodes.
-- perhaps not the most efficient, but only done once per input graph
-- contracts indegree 1 outdegree 1 nodes
resolveNode :: SimpleGraph -> LG.Node -> ([LG.LEdge Double], [LG.LEdge Double]) -> (Int, Int) -> SimpleGraph
resolveNode inGraph curNode inOutPair@(inEdgeList, outEdgeList) (inNum, outNum) =
  if LG.isEmpty inGraph then LG.empty
  else
    --trace ("Resolveing " ++ show (inNum, outNum)) (
    let numNodes = length $ LG.nodes inGraph
    in
    -- isolated node -- throw error
    if inNum == 0 && outNum == 0 then error ("ResolveNode error: Isolated vertex " ++ show curNode ++ " in graph\n" ++ LG.prettify inGraph )

    -- indegree 1 outdegree 1 node to contract
    else if inNum == 1 && outNum == 1 then
      let newEdge = (fst3 $ head inEdgeList, snd3 $ head outEdgeList, 0.0 :: Double)
          newGraph = LG.insEdge newEdge $ LG.delNode curNode $ LG.delLEdges (inEdgeList ++ outEdgeList) inGraph
      in
      newGraph

    -- leaf leaf with too many parents
    else if (inNum > 1) && (outNum == 0) || (inNum > 2) && (outNum == 1) || (inNum > 1) && (outNum == 2) then
      let first2Edges = take 2 inEdgeList
          newNode = (numNodes , T.pack $ ("HTU" ++ (show numNodes)))
          newEdge1 = (fst3 $ head first2Edges, numNodes, 0.0 :: Double)
          newEdge2 = (fst3 $ last first2Edges, numNodes, 0.0 :: Double)
          newEdge3 = (numNodes, curNode, 0.0 :: Double)
          newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
      in
      newGraph

    else if (inNum < 2 || outNum > 2) then
      let  first2Edges = take 2 outEdgeList
           newNode = (numNodes , T.pack $ ("HTU" ++ (show numNodes)))
           newEdge1 = (numNodes, snd3 $ head first2Edges, 0.0 :: Double)
           newEdge2 = (numNodes, snd3 $ last first2Edges, 0.0 :: Double)
           newEdge3 = (curNode, numNodes, 0.0 :: Double)
           newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
      in
      newGraph

      -- node with too parents and too many children
      -- converts to tree node--biased in that direction
      else if (inNum > 2) && (outNum > 2) then
        let first2Edges = take 2 inEdgeList
            newNode = (numNodes , T.pack $ ("HTU" ++ (show numNodes)))
            newEdge1 = (fst3 $ head first2Edges, numNodes, 0.0 :: Double)
            newEdge2 = (fst3 $ last first2Edges, numNodes, 0.0 :: Double)
            newEdge3 = (numNodes, curNode, 0.0 :: Double)
            newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
        in
        newGraph


        -- root or simple network indegree node
      else if (inNum == 0 || outNum > 2) then
          let first2Edges = take 2 outEdgeList
              newNode = (numNodes , T.pack $ ("HTU" ++ (show numNodes)))
              newEdge1 = (numNodes, snd3 $ head first2Edges, 0.0 :: Double)
              newEdge2 = (numNodes, snd3 $ last first2Edges, 0.0 :: Double)
              newEdge3 = (curNode, numNodes, 0.0 :: Double)
              newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
            in
            newGraph


      else error ("This can't happen in resolveNode in/out edge lists don't need to be resolved " ++ show inOutPair ++ "\n" ++ LG.prettify inGraph)
    --)

-- | rerootTree takes a graph and reroots based on a vertex index (usually leaf outgroup)
--   if input is a forest then only roots the component that contains the vertex wil be rerooted
--   unclear how will effect network edges--will need to verify that does not create cycles
--   multi-rooted components (as opposed to forests) are unaffected with trace warning thrown
--   after checking for existing root and multiroots, should be O(n) where 'n is the length
--   of the path between the old and new root
rerootTree :: (Show a, Show b, Eq b) => Int -> LG.Gr a b -> LG.Gr a b
rerootTree rerootIndex inGraph =
  --trace ("In reroot Graph: " ++ show rerootIndex) (
  if LG.isEmpty inGraph then inGraph
  else
    let componentList = LG.components inGraph
        parentNewRootList = LG.pre inGraph rerootIndex
        newRootOrigEdge = head $ LG.inn inGraph rerootIndex
        parentRootList = fmap (LG.isRoot inGraph) parentNewRootList
        outgroupInComponent = fmap (rerootIndex `elem`) componentList
        componentWithOutgroup = filter ((== True).fst) $ zip outgroupInComponent componentList
        -- (_, inNewRoot, outNewRoot) = LG.getInOutDeg inGraph (LG.labelNode inGraph rerootIndex)
    in

    -- rerooting on root so no indegree edges
    -- this for wagner build reroots where can try to reroot on leaf not yet added
    if null $ LG.inn inGraph rerootIndex then inGraph -- error ("Rerooting on indegree 0 node " ++ (show rerootIndex) ++ "\n" ++ LG.prettyIndices inGraph) -- LG.empty


    else if null componentWithOutgroup then inGraph  -- error ("Error rooting wierdness in rerootTree " ++ (show rerootIndex) ++ "\n" ++ LG.prettyIndices inGraph) -- LG.empty

    -- check if new outtaxon has a parent--shouldn't happen-but could if its an internal node reroot
    else if null parentNewRootList || (True `elem` parentRootList) then inGraph
                                                              else (if null componentWithOutgroup then error ("Outgroup index " ++ show rerootIndex ++ " not found in graph")
    else
        --trace ("RRT: " ++ (show (rerootIndex, inNewRoot, outNewRoot))) ( 
        -- reroot component with new outtaxon
        let componentWithNewOutgroup = snd $ head componentWithOutgroup
            (_, originalRootList) =  unzip $ filter ((==True).fst) $ zip (fmap (LG.isRoot inGraph) componentWithNewOutgroup) componentWithNewOutgroup
            numRoots = length originalRootList
            orginalRoot = head originalRootList
            originalRootEdges = LG.out inGraph orginalRoot

        in

        if numRoots == 0 then error ("No root in rerootTree: Attempting to reroot on edge to node " ++ (show rerootIndex) ++ "\n" ++ LG.prettyIndices inGraph) --LG.empty

        -- check if outgroup in a multirooted component
        -- if wagner build this is ok
        else if numRoots > 1 then inGraph -- error ("Error: Attempting to reroot multi-rooted component") -- inGraph
        else
          --reroot graph safely automatically will only affect the component with the outgroup
          -- delete old root edge and create two new edges from oringal root node.
          -- keep orignl root node and delte/crete new edges when they are encounterd
          --trace ("Moving root from " ++ (show orginalRoot) ++ " to " ++  (show rerootIndex)) (
          let leftChildEdge = (orginalRoot, rerootIndex, LG.edgeLabel $ head originalRootEdges)
              rightChildEdge = (orginalRoot, fst3 newRootOrigEdge, LG.edgeLabel $ last originalRootEdges)

              --  this assumes 2 children of old root -- shouled be correct as Phylogenetic Graph
              newEdgeOnOldRoot = if (length originalRootEdges) /= 2 then error ("Number of root out edges /= 2 in rerootGraph: " ++ (show $ length originalRootEdges)
                ++ " root index: " ++ (show (orginalRoot, rerootIndex)) ++ "\nGraph:\n" ++ (LG.prettyIndices inGraph))
                                 else (snd3 $ head originalRootEdges, snd3 $ last originalRootEdges, thd3 $ head originalRootEdges)

              newRootEdges = [leftChildEdge, rightChildEdge, newEdgeOnOldRoot]
              newGraph = LG.insEdges newRootEdges $ LG.delLEdges (newRootOrigEdge : originalRootEdges) inGraph

              -- get edges that need reversing
              newGraph' = preTraverseAndFlipEdges [leftChildEdge,rightChildEdge] newGraph

          in
          --trace ("=")
          --trace ("Deleting " ++ (show (newRootOrigEdge : originalRootEdges)) ++ "\nInserting " ++ (show newRootEdges))
          --trace ("In " ++ (GFU.showGraph inGraph) ++ "\nNew " ++  (GFU.showGraph newGraph) ++ "\nNewNew "  ++  (GFU.showGraph newGraph'))
          newGraph')
        -- ) -- )

-- | rerootDisplayTree like reroot but inputs original root position instead of figuring it out.
-- assumes graph is tree--not useable fo Wagner builds since they have multiple components while building
rerootDisplayTree :: (Show a, Show b, Eq a, Eq b) => LG.Node -> LG.Node -> LG.Gr a b -> LG.Gr a b
rerootDisplayTree orginalRootIndex rerootIndex inGraph =
  --trace ("In reroot Graph: " ++ show rerootIndex) (
  if LG.isEmpty inGraph then inGraph
  else
    let -- componentList = LG.components inGraph'
        
        -- hack---remove when figured out
        {-
        inGraph = if inGraphType == SoftWired then inGraph' -- LG.removeDuplicateEdges inGraph'
                  else inGraph'
        -}

        parentNewRootList = LG.pre inGraph rerootIndex
        newRootOrigEdge = head $ LG.inn inGraph rerootIndex
        parentRootList = fmap (LG.isRoot inGraph) parentNewRootList
        -- outgroupInComponent = fmap (rerootIndex `elem`) componentList
        -- componentWithOutgroup = filter ((== True).fst) $ zip outgroupInComponent componentList
        (_, inNewRoot, outNewRoot) = LG.getInOutDeg inGraph (LG.labelNode inGraph rerootIndex)

        {- Checks for valid trees
        cyclicString = if LG.cyclic inGraph then " Is cyclic "
                          else  " Not cyclic "

        parentInCharinString = if parentInChain inGraph then " Is Parent in Chain "
                                  else " Not Parent in Chain "

        duplicatedsEdgeString = if not (LG.hasDuplicateEdge inGraph)  then " No duplicate edges "
                                else 
                                  let dupEdgeList' = LG.getDuplicateEdges inGraph
                                  in (" Has duplicate edges: " ++ (show dupEdgeList')) 
        -}
    in

    -- trace ("RRT In: " ++ cyclicString ++ parentInCharinString ++ duplicatedsEdgeString) (

    -- check if cycle and exit if so
    -- don't reroot on in=out=1 since same as it descendent edge 
    if (inNewRoot == 1) && (outNewRoot == 1) then 
      -- trace ("RRT: in 1 out 1") 
      inGraph

    -- else if LG.cyclic inGraph then LG.empty -- inGraph 

    
    -- rerooting on root so no indegree edges
    -- this for wagner build reroots where can try to reroot on leaf not yet added
    else if null $ LG.inn inGraph rerootIndex then inGraph -- error ("Rerooting on indegree 0 node " ++ (show rerootIndex) ++ "\n" ++ LG.prettyIndices inGraph) -- LG.empty


    -- else if null componentWithOutgroup then inGraph  -- error ("Error rooting wierdness in rerootTree " ++ (show rerootIndex) ++ "\n" ++ LG.prettyIndices inGraph) -- LG.empty

    -- check if new outtaxon has a parent--shouldn't happen-but could if its an internal node reroot
    else if null parentNewRootList || (True `elem` parentRootList) then inGraph
                                                              -- else if null componentWithOutgroup then error ("Outgroup index " ++ show rerootIndex ++ " not found in graph")
    else
        -- trace ("RRT: " ++ (show (rerootIndex, inNewRoot, outNewRoot))) ( 
        --reroot component with new outtaxon
        let -- componentWithNewOutgroup = snd $ head componentWithOutgroup
            -- (_, originalRootList) =  unzip $ filter ((==True).fst) $ zip (fmap (LG.isRoot inGraph) componentWithNewOutgroup) componentWithNewOutgroup
            -- numRoots = 1 -- length originalRootList
            orginalRoot = orginalRootIndex -- head originalRootList
            originalRootEdges = (LG.out inGraph orginalRoot)

        in

        {-
        if numRoots == 0 then error ("No root in rerootDisplayTree: Attempting to reroot on edge to node " ++ (show (orginalRoot,rerootIndex)) ++ LG.prettyIndices inGraph) --LG.empty

        -- check if outgroup in a multirooted component
        -- if wagner build this is ok
        -- else if numRoots > 1 then inGraph -- error ("Error: Attempting to reroot multi-rooted component") -- inGraph
        else
        -}
          --reroot graph safely automatically will only affect the component with the outgroup
          -- delete old root edge and create two new edges from oringal root node.
          -- keep orignl root node and delte/crete new edges when they are encounterd
          --trace ("Moving root from " ++ (show orginalRoot) ++ " to " ++  (show rerootIndex)) (
          let leftChildEdge = (orginalRoot, rerootIndex, LG.edgeLabel $ head originalRootEdges)
              rightChildEdge = (orginalRoot, fst3 newRootOrigEdge, LG.edgeLabel $ last originalRootEdges)

              --  this assumes 2 children of old root -- shouled be correct as Phylogenetic Graph
              newEdgeOnOldRoot = if (length originalRootEdges) /= 2 then error ("Number of root out edges /= 2 in rerootGraph: " ++ (show $ length originalRootEdges)
                ++ " root index: " ++ (show (orginalRoot, rerootIndex)) ++ "\nGraph:\n" ++ (LG.prettyIndices inGraph))
                                 else (snd3 $ head originalRootEdges, snd3 $ last originalRootEdges, thd3 $ head originalRootEdges)

              newRootEdges = [leftChildEdge, rightChildEdge, newEdgeOnOldRoot]
              newGraph = LG.insEdges newRootEdges $ LG.delLEdges (newRootOrigEdge : originalRootEdges) inGraph

              -- get edges that need reversing
              newGraph' = preTraverseAndFlipEdgesTree orginalRootIndex [leftChildEdge,rightChildEdge] newGraph


              {- Check for valid tree
              cyclicString' = if LG.cyclic newGraph' then " Is cyclic "
                          else  " Not cyclic "

              parentInCharinString'= if parentInChain newGraph' then " Is Parent in Chain "
                                  else " Not Parent in Chain "

              duplicatedsEdgeString' = if not (LG.hasDuplicateEdge newGraph') then " No duplicate edges "
                                       else let dupEdgeList' = LG.getDuplicateEdges newGraph'
                                            in
                                            (" Has duplicate edges: " ++ (show dupEdgeList') ++ "\nDeleting " ++ (show $ fmap LG.toEdge $ (newRootOrigEdge : originalRootEdges)) ++ "\nInserting " ++ (show $ fmap LG.toEdge $ newRootEdges)) 
              -}
          in
          --trace ("=")
          --trace ("In " ++ (GFU.showGraph inGraph) ++ "\nNew " ++  (GFU.showGraph newGraph) ++ "\nNewNew "  ++  (GFU.showGraph newGraph'))

          -- trace ("Deleting " ++ (show $ fmap LG.toEdge (newRootOrigEdge : originalRootEdges)) ++ "\nInserting " ++ (show $ fmap LG.toEdge newRootEdges) 
          --   ++ "\nRRT Out: " ++ cyclicString' ++ parentInCharinString' ++ duplicatedsEdgeString') (

          -- if cyclicString' == " Is cyclic " then inGraph 
          -- else if LG.hasDuplicateEdge newGraph' then LG.removeDuplicateEdges newGraph'
          -- else 
          {-Cycle check
          if LG.cyclic newGraph' then 
            trace ("Orignal root: " ++ (show orginalRootIndex) ++ "New root: " ++ (show rerootIndex) ++ " Deleting " ++ (show $ fmap LG.toEdge (newRootOrigEdge : originalRootEdges)) ++ "\nInserting " ++ (show $ fmap LG.toEdge newRootEdges) 
              ++ "\nOrigGraph: " ++ (LG.prettyIndices inGraph) ++ "\nNewGraph: " ++ (LG.prettyIndices newGraph)++ "\nNewGraph': " ++ (LG.prettyIndices newGraph'))
            LG.empty
          else newGraph'-}
          newGraph'
          -- )
          -- ) -- )


-- | preTraverseAndFlipEdgesTree traverses a tree from starting edge flipping in-edges since they should
-- be out-edges 
-- when recursion its edges that don't need to be fliped then stops
-- assumes input edge is directed correctly
-- follows  traversal out "pre" order updating graph as edges flipped
preTraverseAndFlipEdgesTree :: (Eq b) => LG.Node -> [LG.LEdge b] ->  LG.Gr a b -> LG.Gr a b
preTraverseAndFlipEdgesTree rootIndex inEdgeList inGraph  =
  if null inEdgeList then inGraph
  else
    let -- first edge directled correctly
        inEdge@(_,v,_) = head inEdgeList

        -- edges "in" to child node of first edge--these should be out and need to be flipped
        childEdges = filter ((/= rootIndex) . fst3) $ filter (/= inEdge) $ LG.inn inGraph v
        
        -- flip to "in" to "out" edges
        flippedEdges = fmap LG.flipLEdge childEdges

        -- -- modify graph accordingly
        newGraph = LG.insEdges flippedEdges $ LG.delLEdges childEdges inGraph
    in
    --trace ("PTFE: flipped " ++ (show $ fmap LG.toEdge flippedEdges)) (
    -- edge terminates in leaf or edges in correct orientation
    if null childEdges then preTraverseAndFlipEdgesTree rootIndex (tail inEdgeList) inGraph

    -- edge needs to be reversed to follow through its children from a new graph
    else preTraverseAndFlipEdgesTree rootIndex (flippedEdges ++ (tail inEdgeList)) newGraph
    -- )


-- | preTraverseAndFlipEdges traverses graph from starting edge flipping edges as needed
-- when recursion its edges that don't need to be fliped then stops
-- assumes input edge is directed correctly
-- follows  traversal out "pre" order ish from edges
-- have to check edge orientatins--make sure thay haven't changed as graph was updated earlier
preTraverseAndFlipEdges :: (Eq b) => [LG.LEdge b] ->  LG.Gr a b -> LG.Gr a b
preTraverseAndFlipEdges inEdgelist inGraph  =
  if null inEdgelist then inGraph
  else
    let inEdge@(_,v,_) = head inEdgelist
        childEdges = (LG.out inGraph v) ++ (filter (/= inEdge) $ LG.inn inGraph v)

        -- returns list of edges that had to be flipped
        edgesToFlip = getToFlipEdges v childEdges 
        flippedEdges = fmap LG.flipLEdge edgesToFlip
        newGraph = LG.insEdges flippedEdges $ LG.delLEdges edgesToFlip inGraph
    in
    --trace ("PTFE: flipped " ++ (show $ fmap LG.toEdge flippedEdges)) (
    -- edge terminates in leaf or edges in correct orientation
    if null childEdges  || null edgesToFlip then preTraverseAndFlipEdges (tail inEdgelist) inGraph
    -- edge needs to be reversed to follow through its children from a new graph
    else preTraverseAndFlipEdges (flippedEdges ++ (tail inEdgelist)) newGraph
    -- )

-- | getToFlipEdges takes an index and check edge list
-- and creates new list of edges that need to be flipped
getToFlipEdges ::  LG.Node -> [LG.LEdge b] -> [LG.LEdge b]
getToFlipEdges parentNodeIndex inEdgeList =
  if null inEdgeList then []
  else
    let firstEdge@(u,_,_) = head inEdgeList
    in
    if parentNodeIndex /= u then firstEdge : getToFlipEdges parentNodeIndex (tail inEdgeList)
    else getToFlipEdges parentNodeIndex (tail inEdgeList)

-- | Random generates display trees up to input number by choosing
-- to keep indegree nodes > 1 unifomaly at random
generateDisplayTreesRandom :: (Show a, Show b, Eq a, Eq b, NFData a, NFData b) => Int -> Int -> LG.Gr a b -> [LG.Gr a b]
generateDisplayTreesRandom rSeed numDisplayTrees inGraph =
  if LG.isEmpty inGraph then error "Empty graph in generateDisplayTreesRandom"
  else
    let atRandomList = take numDisplayTrees $ randomIntList rSeed
        randDisplayTreeList = fmap (randomlyResolveGraphToTree inGraph) atRandomList `using` PU.myParListChunkRDS
    in
    randDisplayTreeList

-- | randomlyResolveGraphToTree resolves a single graph to a tree by choosing single indegree edges
-- uniformly at random and deleting all others from graph
-- in=out=1 nodes are contracted, HTU's withn outdegree 0 removed, graph reindexed
-- but not renamed--edges from root are left alone.
randomlyResolveGraphToTree :: (Show a, Show b, Eq a, Eq b) => LG.Gr a b -> Int -> LG.Gr a b
randomlyResolveGraphToTree inGraph randVal =
  if LG.isEmpty inGraph then error "Empty graph in randomlyResolveGraphToTree"
  else
    let (_, leafList, _, _) = LG.splitVertexList inGraph
        -- rootEdgeList = fmap (LG.out inGraph) $ fmap fst rootList
        inEdgeListByVertex = (fmap (LG.inn inGraph) (LG.nodes inGraph)) -- L.\\ rootEdgeList
        randList = fmap abs $ randomIntList randVal
        edgesToDelete = concat $ zipWith  chooseOneDumpRest randList  (fmap (fmap LG.toEdge) inEdgeListByVertex)
        newTree = LG.delEdges edgesToDelete inGraph
        newTree' = LG.removeNonLeafOut0Nodes leafList newTree
        newTree'' = LG.contractIn1Out1Edges newTree'
        reindexTree = LG.reindexGraph newTree''
    in
    -- trace ("RRGT\n" ++ (LG.prettify inGraph) ++ "\n to delete " ++ (show edgesToDelete) ++ "\nNew graph:\n" ++ (LG.prettify newTree)
    --   ++ "\nnewTree'\n" ++ (LG.prettify newTree') ++ "\nnewTree''\n" ++ (LG.prettify newTree'') ++ "\reindex\n" ++ (LG.prettify reindexTree))
    reindexTree

-- | chooseOneDumpRest takes random val and chooses to keep tht edge in list returniong list of edges to delete
chooseOneDumpRest :: Int -> [LG.Edge] -> [LG.Edge]
chooseOneDumpRest randVal inEdgeList =
  if null inEdgeList then []
  else if length inEdgeList == 1 then []
  else
    let numEdgesIn = length inEdgeList
        (_, indexToKeep) = divMod randVal numEdgesIn
    in
    -- trace ("CODR: " ++ (show inEdgeList) ++ " keeping " ++ (show $ inEdgeList !! indexToKeep) ++ " deleting " ++ (show $ inEdgeList L.\\ [inEdgeList !! indexToKeep]) ++ "based on " ++ (show (randVal, indexToKeep)))
    inEdgeList L.\\ [inEdgeList !! indexToKeep]



-- | generateDisplayTrees nice wrapper around generateDisplayTrees' with clean interface
generateDisplayTrees :: (Eq a) => LG.Gr a b -> [LG.Gr a b]
generateDisplayTrees inGraph =
    let (_, leafList, _, _) = LG.splitVertexList inGraph
    in
    generateDisplayTrees' leafList [inGraph] []

-- | generateDisplayTrees' takes a graph list and recursively generates
-- a list of trees created by progresively resolving each network vertex into a tree vertex
-- in each input graph
-- creating up to 2**m (m network vertices) trees.
-- call -> generateDisplayTrees'  [startGraph] []
-- the second and third args contain graphs that need more work and graphs that are done (ie trees)
generateDisplayTrees' :: (Eq a) => [LG.LNode a] -> [LG.Gr a b] -> [LG.Gr a b] -> [LG.Gr a b]
generateDisplayTrees' leafList curGraphList treeList  =
  if null curGraphList then
      let treeList' = fmap (LG.removeNonLeafOut0Nodes leafList) treeList
          treeList'' = fmap LG.contractIn1Out1Edges treeList'
          reindexedTreeList = fmap LG.reindexGraph treeList''
      in
      reindexedTreeList

  else
    let firstGraph = head curGraphList
    in
      if LG.isEmpty firstGraph then []
      else
        let nodeList = LG.labNodes firstGraph
            inNetEdgeList = filter ((>1).length) $ fmap (LG.inn firstGraph) $ fmap fst nodeList
        in
        if null inNetEdgeList then generateDisplayTrees' leafList (tail curGraphList) (firstGraph : treeList)
        else
          let newGraphList = splitGraphListFromNode inNetEdgeList [firstGraph]
          in
          generateDisplayTrees' leafList (newGraphList ++ (tail curGraphList)) treeList

-- | splitGraphListFromNode take a graph and a list of edges for indegree > 1 node
-- removes each in edge in turn to create a new graph and maintains any in 1 out 1 nodes
-- if these edges are to be contracted out  use 'contractOneOneEdges'
-- should be a single traversal 'splitting' graph each time. removing edges anmd recursing
-- to do it again untill all edges are indegree 1.
-- the edges in list for recurvise deleteion should always be in all graphs uuntil
-- the end
splitGraphListFromNode :: [[LG.LEdge b]] -> [LG.Gr a b] -> [LG.Gr a b]
splitGraphListFromNode inEdgeListList inGraphList =
  if null inEdgeListList then inGraphList
  else if null inGraphList then error "Empty graph lits in splitGraphListFromNode"
  else
    let firstNetEdgeList = head inEdgeListList
        indexList = [0.. (length firstNetEdgeList - 1)]
        repeatedEdgeList = replicate (length firstNetEdgeList) firstNetEdgeList
        netEdgeIndexPairList = zip repeatedEdgeList indexList
        newGraphList = concat $ fmap (deleteEdgesCreateGraphs netEdgeIndexPairList 0) inGraphList
    in
    splitGraphListFromNode (tail inEdgeListList) newGraphList

-- | deleteEdgesCreateGraphs takes a list of edges and an index list and a graph,
-- recursively keeps the index-th edge and deltes the others creating a new graph list
deleteEdgesCreateGraphs :: [([LG.LEdge b], Int)] -> Int -> LG.Gr a b -> [LG.Gr a b]
deleteEdgesCreateGraphs netEdgeIndexPairList counter inGraph =
  if LG.isEmpty inGraph then error "Empty graph in  "deleteEdgesCreateGraphs
  else if null netEdgeIndexPairList then []
  else
    let (edgeList, lIndex) = head netEdgeIndexPairList
        --edgeToKeep = edgeList !! index
        edgesToDelete = (take lIndex edgeList) ++ (drop (lIndex + 1) edgeList)
        newGraph = LG.delLEdges edgesToDelete inGraph
    in
    newGraph : deleteEdgesCreateGraphs (tail netEdgeIndexPairList) (counter + 1) inGraph

-- | getNodeType returns node type for Node
getNodeType :: (Show a, Show b) => LG.Gr a b -> LG.Node -> NodeType
getNodeType inGraph inNode =
    if not $ LG.gelem inNode inGraph then error ("Node " ++ (show inNode) ++ " not in graph\n" ++ (GFU.showGraph inGraph))
    else if LG.isLeaf inGraph inNode then LeafNode
    else if LG.isTreeNode inGraph inNode then TreeNode
    else if LG.isNetworkNode inGraph inNode then NetworkNode
    else if LG.isRoot inGraph inNode then RootNode
    else error ("Node type " ++ (show inNode) ++ " not Leaf, Tree, Network, or Root in graph\n" ++ (GFU.showGraph inGraph))

-- | convertDecoratedToSimpleGraph
convertDecoratedToSimpleGraph :: DecoratedGraph -> SimpleGraph
convertDecoratedToSimpleGraph inDec =
  if LG.isEmpty inDec then LG.empty
  else
    let decNodeList = LG.labNodes inDec
        newNodeLabels = fmap vertName $ fmap snd decNodeList
        simpleNodes = zip (fmap fst decNodeList) newNodeLabels

        {-
        decEdgeList = LG.labEdges inDec
        sourceList = fmap fst3 decEdgeList
        sinkList = fmap snd3 decEdgeList
        newEdgeLables = replicate (length sourceList) 0.0  -- fmap midRangeLength $ fmap thd3 decEdgeList
        simpleEdgeList = zip3 sourceList sinkList newEdgeLables
        -}
        simpleEdgeList = fmap convertToSimpleEdge $ LG.labEdges inDec
    in
    LG.mkGraph simpleNodes simpleEdgeList


-- | convertToSimpleEdge takes a lables edge and relabels with 0.0
convertToSimpleEdge :: LG.LEdge EdgeInfo -> LG.LEdge Double
convertToSimpleEdge (a, b, c) = (a, b, minLength c)

-- | graphCostFromNodes takes a Decorated graph and returns its cost by summing up the local costs
--  of its nodes
graphCostFromNodes :: DecoratedGraph -> Double
graphCostFromNodes inGraph =
  if LG.isEmpty inGraph then 0.0
  else
    sum $ fmap vertexCost $ fmap snd $ LG.labNodes inGraph

-- | switchRootTree takes a new root vertex index of a tree and switches the existing root (and all relevent edges)
-- to new index
switchRootTree :: (Show a) => Int -> LG.Gr a b -> LG.Gr a b
switchRootTree newRootIndex inGraph =
  if LG.isEmpty inGraph then LG.empty
  else
    let rootList = LG.getRoots inGraph
        (newRootCurInEdges, newRootCurOutEdges) = LG.getInOutEdges inGraph newRootIndex
        oldRootEdges = LG.out inGraph $ fst $ head rootList

    in
    -- not a directed tree
    if length rootList /= 1 then error ("Graph input to switchRootTree is not a tree--not single root:" ++ (show rootList))

    -- same root
    else if newRootIndex == (fst $ head rootList) then inGraph
    else
        -- create new edges and delete the old ones
        let newEdgesToAdd = fmap (flipVertices (fst $ head rootList) newRootIndex) (newRootCurInEdges ++ newRootCurOutEdges ++ oldRootEdges)
        in
        LG.insEdges newEdgesToAdd $ LG.delLEdges (newRootCurInEdges ++ newRootCurOutEdges ++ oldRootEdges) inGraph

-- | flipVertices takes an old vertex index and a new vertex index and inserts one for the other
-- in a labelled edge
flipVertices :: Int -> Int ->  LG.LEdge b -> LG.LEdge b
flipVertices a b (u,v,l) =
  let newU = if u == a then b
             else if u == b then a
             else u
      newV = if v == a then b
            else if v == b then a
            else v
  in
  -- trace (show (u,v,l) ++ "->" ++ show (newU, newV, l))
  (newU, newV, l)

-- | dichotomizeRoot takes greaph and dichotimizes not dichotomous roots in graph
dichotomizeRoot :: Int -> SimpleGraph -> SimpleGraph
dichotomizeRoot lOutgroupIndex inGraph =
  if LG.isEmpty inGraph then LG.empty
  else
    let rootList = LG.getRoots inGraph
        currentRoot = fst $ head rootList
        rootEdgeList = LG.out inGraph $ currentRoot
    in
    -- not a tree error
    if (length rootList /= 1) then error ("Graph input to dichotomizeRoot is not a tree--not single root:" ++ (show rootList))

    -- nothing to do
    else if (length rootEdgeList < 3) then inGraph
    else
      let numVertices = length $ LG.nodes inGraph
          newNode = (numVertices, T.pack $ show numVertices)
          edgesToDelete = filter ((/=lOutgroupIndex) . snd3) rootEdgeList
          newEdgeDestinations = fmap snd3 edgesToDelete
          newEdgeStarts = replicate (length newEdgeDestinations) numVertices
          newEdgeLabels = replicate (length newEdgeDestinations) 0.0
          -- nub for case where root edge in "wrong" direction
          -- doesn't filter edges to delete properly
          newEdgesNewNode = L.nub $ zip3 newEdgeStarts newEdgeDestinations newEdgeLabels
          newRootEdge = (currentRoot, numVertices, 0.0)
      in
      LG.delLEdges edgesToDelete $ LG.insEdges (newRootEdge : newEdgesNewNode) $ LG.insNode newNode inGraph


-- | showBlockGraphs takes a vector of vector of DecoratedGraphs and converte and prettifies outputting a String
showDecGraphs :: V.Vector (V.Vector DecoratedGraph) -> String
showDecGraphs inDecVV =
    if V.null inDecVV then []
    else
        concat $ fmap concat $ V.toList $ fmap V.toList $ fmap (fmap LG.prettify) $ fmap (fmap convertDecoratedToSimpleGraph) inDecVV

-- | selectPhylogeneticGraph takes  a series OF arguments and an input list ot PhylogeneticGraphs
-- and returns or filters that list based on options.
-- uses selectListCostPairs in GeneralUtilities
selectPhylogeneticGraph :: [Argument] -> Int -> [String] -> [PhylogeneticGraph] -> [PhylogeneticGraph]
selectPhylogeneticGraph inArgs rSeed selectArgList curGraphs =
    if null curGraphs then []
    else
        let fstArgList = fmap (fmap C.toLower . fst) inArgs
            sndArgList = fmap (fmap C.toLower . snd) inArgs
            lcArgList = zip fstArgList sndArgList
            checkCommandList = checkCommandArgs "select" fstArgList selectArgList
        in
           -- check for valid command options
           if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'select': " ++ show inArgs)
           else if length inArgs > 1 then errorWithoutStackTrace ("Can only have a single select type per command: "  ++ show inArgs)
           else
                let doBest    = not $ not (any ((=="best").fst) lcArgList)
                    doAll     = not $ not (any ((=="all").fst) lcArgList)
                    doRandom  = not $ not (any ((=="atrandom").fst) lcArgList)
                    doUnique  = not $ not (any ((=="unique").fst) lcArgList)
                    numberToKeep
                      | null lcArgList = Just (maxBound :: Int)
                      | null $ snd $ head lcArgList = Just (maxBound :: Int)
                      | otherwise = readMaybe (snd $ head lcArgList) :: Maybe Int
                in
                if doAll then curGraphs
                else if isNothing numberToKeep then errorWithoutStackTrace ("Number to keep specification not an integer: "  ++ show (snd $ head lcArgList))
                else
                    let -- minimum graph cost
                        minGraphCost = minimum $ fmap snd6 curGraphs

                        -- collapse zero-length branchs for unique
                        curGraphsCollapsed = fmap U.collapseGraph curGraphs

                        -- keep only unique graphs based on non-zero edges--in sorted by cost
                        uniqueGraphList = L.sortOn snd6 $ getUniqueGraphs'' (zip curGraphs curGraphsCollapsed)-- curGraphs --  True curGraphs -- getBVUniqPhylogeneticGraph True curGraphs -- getTopoUniqPhylogeneticGraph True curGraphs
                        
                        -- this to avaoid alot of unncesesary graph comparisons for 'best' graphs
                        bestCostGraphs = filter ((== minGraphCost).snd6) curGraphs
                        uniqueBestGraphs = getUniqueGraphs'' (zip bestCostGraphs (fmap U.collapseGraph bestCostGraphs))

                      in
                    if doUnique then take (fromJust numberToKeep) uniqueGraphList
                    else if doBest then
                     -- trace ("SPG: " ++ (show (minGraphCost, length uniqueGraphList, fmap snd6 uniqueGraphList)))
                      take (fromJust numberToKeep) uniqueBestGraphs
                    else if doRandom then
                         let randList = head $ shuffleInt rSeed 1 [0..(length curGraphs - 1)]
                             (_, shuffledGraphs) = unzip $ L.sortOn fst $ zip randList curGraphs
                         in
                         take (fromJust numberToKeep) shuffledGraphs
                    -- default is all best and unique
                    else
                        uniqueBestGraphs

-- | getUniqueGraphs takes each pair of non-zero edges and conpares them--if equal not added to list
-- maybe chnge to nub LG.pretify graphList?
getUniqueGraphs :: Bool -> [PhylogeneticGraph] -> [PhylogeneticGraph]
getUniqueGraphs removeZeroEdges inGraphList =
  if null inGraphList then []
  else
    let inGraphEdgeList = if removeZeroEdges then fmap (filter ((> 0.0) . minLength . thd3)) $ fmap LG.labEdges $ fmap thd6 inGraphList
                          else fmap LG.labEdges $ fmap thd6 inGraphList
    in
    getUniqueGraphs' (zip inGraphEdgeList inGraphList) []


-- | getUniqueGraphs Using fgl ==
-- basically a nub
-- need to add a collapse function for compare as well
-- takes pairs of (noCollapsed, collapsed) phylogenetic graphs,
-- mke strings based on collapsed and returns not collpased
getUniqueGraphs'' :: [(PhylogeneticGraph, PhylogeneticGraph)] -> [PhylogeneticGraph] 
getUniqueGraphs'' inList = nubGraph [] inList

-- | keeps and returns unique graphs based on Eq of Topological Simple Graph
-- String newick w/0 HTU names and branch lengths
-- arbitrarily rooted on 0 for oonsistency
--reversed to keep original order in case sorted on length
nubGraph :: [(PhylogeneticGraph, PhylogeneticGraph, String)] -> [(PhylogeneticGraph, PhylogeneticGraph)] -> [PhylogeneticGraph]
nubGraph curList inList =
  if null inList then reverse $ fmap fst3 curList
  else 
    let (firstGraphNC, firstGraphC) = head inList
        firstString = makeNewickList False False 0 [fst6 firstGraphC] [snd6 firstGraphNC] 
        isMatch = filter (== firstString) (fmap thd3 curList)
    in
    -- trace ("NG: " ++ (show $ null isMatch) ++ " " ++ firstString) (
    if null curList then nubGraph [(firstGraphNC, firstGraphC, firstString)] (tail inList)
    else if null isMatch then nubGraph ((firstGraphNC, firstGraphC, firstString) : curList) (tail inList)
    else nubGraph curList (tail inList)
    -- )

-- | getUniqueGraphs takes each pair of non-zero edges and compares them--if equal not added to list
getUniqueGraphs' :: [([LG.LEdge EdgeInfo], PhylogeneticGraph)] -> [([LG.LEdge EdgeInfo], PhylogeneticGraph)]  -> [PhylogeneticGraph]
getUniqueGraphs' inGraphPairList currentUniquePairs =
    if null inGraphPairList then fmap snd currentUniquePairs
    else
        let firstPair@(firstEdges, _) = head inGraphPairList
        in
        if null currentUniquePairs then getUniqueGraphs' (tail inGraphPairList) [firstPair]
        else
            let equalList = filter (== True) $ fmap ((== firstEdges) . fst) currentUniquePairs
            in
            if null equalList then getUniqueGraphs' (tail inGraphPairList) (firstPair : currentUniquePairs)
            else getUniqueGraphs' (tail inGraphPairList) currentUniquePairs



{-
-- getNonZeroEdges takes a DecortatedGraph and returns the sorted list of non-zero length (< epsilon) edges
getNonZeroEdges :: PhylogeneticGraph -> ([LG.LEdge EdgeInfo], PhylogeneticGraph)
getNonZeroEdges inGraph =
    if LG.isEmpty $ thd6 inGraph then ([], (LG.empty,0.0, LG.empty, V.empty, V.empty, V.empty))
    else
        let edgeList = LG.labEdges (thd6 inGraph)
            minCostEdgeList = fmap (minLength . thd3) edgeList
            (_, nonZeroEdgeList) = unzip $ filter ((>epsilon) . fst) $ zip  minCostEdgeList edgeList
        in
        (L.sortOn fst3 nonZeroEdgeList, inGraph)
-}

-- | copyIAFinalToPrelim takes a Decorated graph and copies
-- the IA final fields to preliminary IA states--this for IA only optimization
-- inswapping and other operations.  Thi sis done becasue the "preliminary" IA states
-- are only known after full post/pre traversals
copyIAFinalToPrelim :: DecoratedGraph -> DecoratedGraph
copyIAFinalToPrelim inGraph =
  if LG.isEmpty inGraph then error "Empty input graph in copyIAFinalToPrelim"
  else
    let nodes = LG.labNodes inGraph
        edges = LG.labEdges inGraph
        newNodes = fmap makeIAPrelimFromFinal nodes
    in
    LG.mkGraph newNodes edges

-- | makeIAPrelimFromFinal updates the label of a node for IA states
-- setting preliminary to final
makeIAPrelimFromFinal :: LG.LNode VertexInfo -> LG.LNode VertexInfo
makeIAPrelimFromFinal (inIndex, label) =
  let labData = vertData label
      newLabData = fmap (fmap f) labData
  in
  (inIndex, label {vertData = newLabData})
  where f c = if (GV.null $ slimIAFinal c) && (GV.null  $ wideIAFinal c) && (GV.null  $ hugeIAFinal c) then c
              else if (not $ GV.null $ slimIAFinal c) then c {slimIAPrelim = M.makeDynamicCharacterFromSingleVector $ slimIAFinal c}
              else if (not $ GV.null $ wideIAFinal c) then c {wideIAPrelim = M.makeDynamicCharacterFromSingleVector $ wideIAFinal c}
              else c {hugeIAPrelim = M.makeDynamicCharacterFromSingleVector $ hugeIAFinal c}


-- | copyIAPrelimToFinal takes a Decorated graph and copies
-- the IA prelim fields to final IA states--this for IA only optimization
-- inswapping and other operations.  THis is fdone for root and leaf vertices
copyIAPrelimToFinal :: DecoratedGraph -> DecoratedGraph
copyIAPrelimToFinal inGraph =
  if LG.isEmpty inGraph then error "Empty input graph in copyIAFinalToPrelim"
  else
    let nodes = LG.labNodes inGraph
        edges = LG.labEdges inGraph
        newNodes = fmap makeIAFinalFromPrelim nodes
    in
    LG.mkGraph newNodes edges

-- | makeIAFinalFomPrelim updates the label of a node for IA states
-- setting final to preliminary
makeIAFinalFromPrelim:: LG.LNode VertexInfo -> LG.LNode VertexInfo
makeIAFinalFromPrelim (inIndex, label) =
  let labData = vertData label
      newLabData = fmap (fmap f) labData
  in
  (inIndex, label {vertData = newLabData})
  where f c = let newSlimIAFinal = extractMediansGapped $ slimIAPrelim c
                  newWideIAFinal = extractMediansGapped $ wideIAPrelim c
                  newHugeIAFinal = extractMediansGapped $ hugeIAPrelim c
              in
              if (GV.null $ snd3 $ slimIAPrelim c) && (GV.null $ snd3 $ wideIAPrelim c) && (GV.null $ snd3 $ hugeIAPrelim c) then c
              else if (not $ GV.null $ snd3 $ slimIAPrelim c) then c {slimIAFinal = newSlimIAFinal}
              else if (not $ GV.null $ snd3 $ wideIAPrelim c) then c {wideIAFinal = newWideIAFinal}
              else c {hugeIAFinal = newHugeIAFinal}


-- | getTopoUniqPhylogeneticGraph takes a list of phylogenetic graphs and returns
-- list of topologically unique graphs--operatres on simple graph field
-- noZeroEdges flag passed to remove zero weight edges
getTopoUniqPhylogeneticGraph :: Bool -> [PhylogeneticGraph] -> [PhylogeneticGraph]
getTopoUniqPhylogeneticGraph nonZeroEdges inPhyloGraphList =
  if null inPhyloGraphList then []
  else
      let uniqueBoolList = createUniqueBoolList nonZeroEdges (fmap fst6 inPhyloGraphList) []
          boolPair = zip inPhyloGraphList uniqueBoolList
      in
      fmap fst $ filter ((== True) . snd) boolPair

-- | createUniqueBoolList creates a list of Bool if graphs are unique--first occurrence is True, others False
createUniqueBoolList :: Bool -> [SimpleGraph] -> [(SimpleGraph,Bool)] -> [Bool]
createUniqueBoolList nonZeroEdges inGraphList boolAccum =
  if null inGraphList then reverse $ fmap snd boolAccum
  else
    let firstGraph = head inGraphList
    in
    if null boolAccum then createUniqueBoolList nonZeroEdges (tail inGraphList) ((firstGraph,True) : boolAccum)
    else
        let checkList = filter (== True) $ fmap (topologicalEqual nonZeroEdges firstGraph) (fmap fst boolAccum)
        in
        if null checkList then createUniqueBoolList nonZeroEdges (tail inGraphList) ((firstGraph,True) : boolAccum)
        else createUniqueBoolList nonZeroEdges (tail inGraphList) ((firstGraph, False) : boolAccum)



-- | topologicalEqual takes two simple graphs and returns True if graphs have same nodes and edges
-- option to exclude zero weight edges
topologicalEqual :: Bool -> SimpleGraph -> SimpleGraph -> Bool
topologicalEqual nonZeroEdges g1 g2 =
  if LG.isEmpty g1 && LG.isEmpty g2 then True
  else if  LG.isEmpty g1 || LG.isEmpty g2 then False
  else
      let nodesG1 = LG.labNodes g1
          nodesG2 = LG.labNodes g2
          edgesG1 = if nonZeroEdges then fmap LG.toEdge $ filter ((> 0). thd3) $ LG.labEdges g1
                    else LG.edges g1
          edgesG2 = if nonZeroEdges then fmap LG.toEdge $ filter ((> 0). thd3) $ LG.labEdges g2
                    else LG.edges g2
      in
      if nodesG1 == nodesG2 && edgesG1 == edgesG2 then True
      else False

{-
-- | topologicalBVEqual takes two Decorated graphs and returns True if graphs have same nodes
-- by bitvector (sorted)
-- option to exclude nodes that are end of zero weight edges
topologicalBVEqual :: Bool -> DecoratedGraph -> DecoratedGraph -> Bool
topologicalBVEqual nonZeroEdges g1 g2 =
  if LG.isEmpty g1 && LG.isEmpty g2 then True
  else if  LG.isEmpty g1 || LG.isEmpty g2 then False
  else
      let bvNodesG1 = L.sort $ fmap bvLabel $ fmap snd $ LG.labNodes g1
          bvNodesG2 = L.sort $ fmap bvLabel $ fmap snd $ LG.labNodes g2
      in
      if bvNodesG1 == bvNodesG2 then True
      else False
-}

-- | getEdgeMinLengthToNode takes a labelled node and returns the min length of
-- the edge leading to the node
getEdgeMinLengthToNode ::[LG.LEdge EdgeInfo] ->  LG.LNode a -> Double
getEdgeMinLengthToNode  edgeList (node, _)=
  let foundEdge = L.find ((== node) . snd3) edgeList
  in
  -- root node will be nor be in in edge set and need so set > 0
  if foundEdge == Nothing then 1.0 --  error ("Edge not found in getEdgeMinLengthToNode: node " ++ (show node) ++ " edge list " ++ (show edgeList))
  else minLength $ thd3 $ fromJust foundEdge

-- | getBVUniqPhylogeneticGraph takes a list of phylogenetic graphs and returns
-- list of topologically unique graphs based on their node bitvector assignments
-- operatres on Decorated graph field
-- noZeroEdges flag passed to remove zero weight edges
getBVUniqPhylogeneticGraph :: Bool -> [PhylogeneticGraph] -> [PhylogeneticGraph]
getBVUniqPhylogeneticGraph nonZeroEdges inPhyloGraphList =
  if null inPhyloGraphList then []
  else
      let bvGraphList = fmap (getBVNodeList nonZeroEdges) $ fmap thd6 inPhyloGraphList
          uniqueBoolList = createBVUniqueBoolList bvGraphList []
          boolPair = zip inPhyloGraphList uniqueBoolList
      in
      fmap fst $ filter ((== True) . snd) boolPair


-- | getBVNodeList takes a DecoratedGraph and returns sorted list (by BV) of nodes
-- removes node with zero edge weight to them if specified
getBVNodeList :: Bool -> DecoratedGraph -> [BV.BitVector]
getBVNodeList nonZeroEdges inGraph =
  if LG.isEmpty inGraph then []
  else
      let nodeList =  LG.labNodes inGraph
          edgeList = LG.labEdges inGraph
          minLengthList = fmap (getEdgeMinLengthToNode edgeList) nodeList
          nodePairList = filter ((> 0) . snd) $ zip nodeList minLengthList
          bvNodeList  =  if nonZeroEdges then L.sort $ fmap bvLabel $ fmap snd $ fmap fst nodePairList
                         else L.sort $ fmap bvLabel $ fmap snd nodeList
          in
          bvNodeList

-- | createBVUniqueBoolList creates a list of Bool if graphs are unique by bitvecector node list
-- first occurrence is True, others False
-- assumes edges filterd b=y lenght already
createBVUniqueBoolList :: [[BV.BitVector]] -> [([BV.BitVector],Bool)] -> [Bool]
createBVUniqueBoolList inBVGraphListList boolAccum =
  if null inBVGraphListList then reverse $ fmap snd boolAccum
  else
    let firstGraphList = head inBVGraphListList
    in
    if null boolAccum then createBVUniqueBoolList  (tail inBVGraphListList) ((firstGraphList,True) : boolAccum)
    else
        let checkList = filter (== True) $ fmap (== firstGraphList) (fmap fst boolAccum)
        in
        if null checkList then createBVUniqueBoolList  (tail inBVGraphListList) ((firstGraphList,True) : boolAccum)
        else createBVUniqueBoolList  (tail inBVGraphListList) ((firstGraphList, False) : boolAccum)

-- | makeDummyLabEdge takes an unlabelled edge and adds a dummy label
makeDummyLabEdge :: EdgeInfo -> LG.Edge -> LG.LEdge EdgeInfo
makeDummyLabEdge edgeLab (u,v) = (u,v,edgeLab)

-- | netNodeAncestorViolation checks whether one of the edge into a netowrk node (in 2)
-- is cinnected to an ancestor (via the other parent) of the node
-- this is a form of time violation since the parents of a network node must be
-- at least possibly coeval
-- this uses the bit vector label of nodes.  If the other child of either parent node
-- of a network node has non-zero intersection between the BV label of the network node
-- and that other child of parent then they conecting edge is from and ancestral node hence a time violation
-- O(n) n netork nodes in Graph, but checks all nodes to see if network
hasNetNodeAncestorViolation :: LG.Gr VertexInfo b -> Bool
hasNetNodeAncestorViolation inGraph =
  if LG.isEmpty inGraph then error "Empty graph in hasNetNodeAncestorViolation"
  else
    let (_, _, _, netWorkNodeList) =  LG.splitVertexList inGraph
        hasAncViolationList = filter (== True) $ fmap (nodeAncViolation inGraph) netWorkNodeList
    in
    -- trace ("HNV: " ++ (show $ (not . null) hasAncViolationList))
    (not . null) hasAncViolationList

-- | nodeAncViolation checks a single node fo ancestrpo connection--he ceviolation
-- should be O(1).  Return True if violation
nodeAncViolation :: LG.Gr VertexInfo b -> LG.LNode VertexInfo -> Bool
nodeAncViolation inGraph inNode =
  let parentList = LG.labParents inGraph (fst inNode)
  in
  if length parentList /= 2 then error ("Parent number should be 2: " ++ (show $ fst inNode) ++ " <- " ++ (show $ fmap fst parentList))
  else
    let sisterNodes = concatMap (LG.sisterLabNodes inGraph) parentList
        sisterBVData = fmap (bvLabel . snd) sisterNodes
        inNodeBVData = bvLabel $ snd inNode
        sisterBVIntersections = fmap (.&. inNodeBVData) sisterBVData
        isAncInNode = filter (== inNodeBVData) sisterBVIntersections
    in
    (not . null) isAncInNode

-- | selectGraphStochastic takes a list of graphs and retuns a list of graphs chosen at Random
-- using an exponential distribution based on graph cost difference divided by an input factor
-- if factor is 0 then stringth graphs cost
-- mprob acceptance = -exp [(cost - minCost)/ factor]
-- returns n graphs by random criterion without replacment
selectGraphStochastic :: Int -> Int -> Double -> [PhylogeneticGraph] -> [PhylogeneticGraph]
selectGraphStochastic rSeed number factor inGraphList =
  if null inGraphList then inGraphList
  else if number >= length inGraphList then inGraphList
  else
    let randList' = randomIntList rSeed
        randList = fmap abs (tail randList')
        newSeed = head randList'
        minCost = minimum $ fmap snd6 inGraphList
        deltaList = fmap ((-) minCost) $ fmap snd6 inGraphList
        probAcceptList = fmap (getProb factor) deltaList

        -- multiplier for resolution 1000, 100 prob be ok
        randMultiplier = 1000
        randMultiplier' = fromIntegral randMultiplier
        intAcceptList = fmap floor $ fmap (* randMultiplier')  probAcceptList
        (_, intRandValList) = unzip $ zipWith divMod randList (replicate (length inGraphList) randMultiplier)
        acceptList = zipWith (<) intRandValList intAcceptList

        -- zip graphs with Bools
        (returnGraphList, _) = unzip $ filter ((== True) .snd) $ zip inGraphList acceptList

        -- takes some random remainder to fiill out length of list
        numLucky = number - (length returnGraphList)
        luckyList = if numLucky > 0 then takeRandom newSeed numLucky (fmap fst $ filter ((== False) .snd) $ zip inGraphList acceptList)
                    else []


    in
    trace ("SGS " ++ (show intAcceptList) ++ " " ++ (show intRandValList) ++ " -> " ++ (show acceptList))
    -- so no more than specified
    take number $ returnGraphList ++ luckyList

    where getProb a b = exp ((-1) * b / a)
