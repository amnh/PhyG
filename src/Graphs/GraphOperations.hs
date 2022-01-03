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

module Graphs.GraphOperations ( ladderizeGraph
                               , verifyTimeConsistency
                               , rerootTree
                               , rerootTree'
                               , generateDisplayTrees
                               , contractOneOneEdges
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
                               ) where

import qualified Data.List                 as L
import qualified Data.BitVector.LittleEndian as BV
import qualified Data.Text.Lazy            as T
import qualified Data.Vector               as V
import           Debug.Trace
import           GeneralUtilities
import qualified GraphFormatUtilities      as GFU
import           Types.Types
import qualified Utilities.LocalGraph      as LG
import qualified Data.Char              as C
import qualified Utilities.Utilities    as U
import           Data.Maybe
import qualified Data.Vector.Generic                                         as GV
import qualified Data.Vector.Storable         as SV
import           Text.Read
import qualified GraphOptimization.Medians as M
import           Bio.DynamicCharacter

-- import Debug.Debug

-- | getEdgeSplitList takes a graph and returns list of edges
-- that split a graph increasing the number of components by 1
-- this is quadratic 
-- should change to Tarjan's algorithm (linear)
--everyhting else in there is O(n^2-3) so maybe doesn't matter 
getEdgeSplitList :: (Eq b) => LG.Gr a b -> [LG.LEdge b]
getEdgeSplitList inGraph = 
  if LG.isEmpty inGraph then error ("Empty graph in getEdgeSplitList")
  else 
      let origNumComponents = LG.noComponents inGraph
          origEdgeList = LG.labEdges inGraph
          edgeDeleteComponentNumberList = fmap LG.noComponents $ fmap (flip LG.delEdge inGraph) (fmap LG.toEdge origEdgeList)
          bridgeList =  fmap snd $ filter ((> origNumComponents) . fst) $ zip edgeDeleteComponentNumberList origEdgeList

          -- filkter out edges starting in an outdegree 1 node (network or in out 1) node
          -- this would promote an HTU to a leaf
          bridgeList' = filter  ((not . LG.isOutDeg1Node inGraph) . fst3) bridgeList
      in
      -- trace ("AP: " ++ (show $ LG.ap $ LG.undir inGraph) ++ "GESL: Components: " ++ (show edgeDeleteComponentNumberList))
      bridgeList'


-- | splitGraphOnEdge takes a graph and an edge and returns a single graph but with two components
-- the roots of each component are retuned with two graphs, with broken edge contraced, and 'naked'
-- node returned.  The naked node is used for rejoining the two components during rearrangement
-- (SplitGraph, root of component that has original root, root of component that was cut off, naked node left over)
-- this function does not check whether edge is a 'bridge'
splitGraphOnEdge :: LG.Gr a b -> LG.LEdge b -> (LG.Gr a b, LG.Node, LG.Node, LG.Node)
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
          (splitGraph, fst $ head $ LG.getRoots inGraph, v, e)


-- | joinGraphOnEdge takes a graph with and adds an edge reducing the component number
-- expected ot be two components to one in SPR/TBR
-- assumes that first node of edge (e,v,l) is 'naked' ie avaiable to make edges but is in graph
-- created from splitGraphOnEdge
joinGraphOnEdge :: LG.Gr a b -> LG.LEdge b -> LG.Node -> LG.Node -> LG.Gr a b
joinGraphOnEdge inGraph edgeToInvade@(x,y,l) nakedNode graphToJoinRoot =
  if LG.isEmpty inGraph then error ("Empty graph in joinGraphOnEdge")
  else 
      let edgeToCreate0 = (x, nakedNode, l)
          edgeToCreate1 = (nakedNode, y, l)
          edgeToCreate2 = (nakedNode, graphToJoinRoot, l)     
      in
      -- make new graph
      LG.insEdges [edgeToCreate0, edgeToCreate1, edgeToCreate2] $ LG.delEdge (x,y) inGraph
 
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

    -- leaf leaf with too manuy parents
    else if ((inNum > 1) && (outNum == 0)) || ((inNum > 2) && (outNum == 1)) then (
      let first2Edges = take 2 inEdgeList
          newNode = (numNodes , T.pack $ ("HTU" ++ (show numNodes)))
          newEdge1 = (fst3 $ head first2Edges, numNodes, 0.0 :: Double)
          newEdge2 = (fst3 $ last first2Edges, numNodes, 0.0 :: Double)
          newEdge3 = (numNodes, curNode, 0.0 :: Double)
          newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
      in
      newGraph) else (if (inNum < 2 || outNum > 2) then
                   let first2Edges = take 2 outEdgeList
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


                 else error ("This can't happen in resolveNode in/out edge lists don't need to be resolved " ++ show inOutPair ++ "\n" ++ LG.prettify inGraph))
    --)

-- | rerootTree' flipped version of rerootGraph
rerootTree' :: (Eq b) => LG.Gr a b -> Int -> LG.Gr a b
rerootTree' inGraph rerootIndex = rerootTree rerootIndex inGraph

{-
--NB--DOES NOT WORK FOR GRAPOHS at least
-- | rerootGraph reroots a single rooted graph
-- for a graph this may result in multiple graphs if rerooting on indegree > 1 (netowrk) node
-- rootIndex is an invariant always the root--no nodes are creted or destroyed

-- 1) removes exisiting root edges (must be only 1) (rootIndex -> a, rootINdex -> b)
-- 2) creates new edge where root was (a->b)
-- 3) checks if rerootindex node has indegree > 1, if so splits call to subreroot function
-- 4) sub reroot removes edge terminating in rerootIndex (c -> rerootIndex)
-- 5) adds two new reroot edges (rootIndex -> c, rootIndex -> rerootIndex)
-- 6) flips edges as needed via preorder traversal till edge flipping no longer needed (incremetnal optimization) 
-- each time remaking graph
rerootGraph :: (Eq b) => LG.Gr a b -> Int -> Int -> [LG.Gr a b]
rerootGraph inGraph inRootIndex rerootIndex = 
  if LG.isEmpty inGraph then [inGraph]
  else 
    let rootList = LG.getRoots inGraph
        origRootIndex = if inRootIndex /= (-1) then inRootIndex
                        else if length rootList /= 1 then error ("Root number <> 1: " ++ (show $ fmap fst rootList))
                        else fst $ head rootList
        origRootEdgeList = if (length $ LG.out inGraph origRootIndex) /= 2 then error ("Require two edges at root:" ++ (show $ fmap LG.toEdge $ LG.out inGraph origRootIndex))
                           else LG.out inGraph origRootIndex

        -- this if to check for leaf as end of edge
        origRootNewEdge = (snd3 $ head origRootEdgeList, snd3 $ last origRootEdgeList, thd3 $ head origRootEdgeList)
        rerootEdgeList = LG.inn inGraph rerootIndex
        newGraphList = fmap (rerootGraph' inGraph origRootIndex rerootIndex origRootEdgeList origRootNewEdge) rerootEdgeList
    in
    trace ("NRE:" ++ (show $ LG.toEdge origRootNewEdge))
    newGraphList

-- | rerootGraph' takes single edge and roots on that edge, creting new edges and directing remaining edges
rerootGraph' :: (Eq b) => LG.Gr a b -> Int -> Int -> [LG.LEdge b] -> LG.LEdge b -> LG.LEdge b-> LG.Gr a b 
rerootGraph' inGraph origRootIndex rerootIndex origRootEdgeList origRootNewEdge rerootEdge =
  let -- this due to edge not getting redirected if connected to edge involved in new root or is leaf
      origRootNewEdge' = if (LG.isLeaf inGraph $ fst3 origRootNewEdge) then LG.flipLEdge origRootNewEdge
                         else if ((snd3 origRootNewEdge) `elem` [fst3 rerootEdge, snd3 rerootEdge]) then LG.flipLEdge origRootNewEdge 
                         else origRootNewEdge

      newRerootEdges = [(origRootIndex, fst3 rerootEdge, thd3 rerootEdge) , (origRootIndex, snd3 rerootEdge, thd3 rerootEdge)]

      newGraph = LG.insEdges (origRootNewEdge' : newRerootEdges) $ LG.delLEdges (rerootEdge : origRootEdgeList) inGraph

      -- need to check for edges indegree to root children and flip
      inEdgesToRootChldren = concatMap (LG.inn newGraph) [fst3 rerootEdge,  snd3 rerootEdge]
      notFromRootInEdges = inEdgesToRootChldren L.\\ newRerootEdges

      -- update graph
      newGraph' = LG.insEdges (fmap LG.flipLEdge notFromRootInEdges) $ LG.delLEdges notFromRootInEdges newGraph
  in
  trace ("Inserting :" ++ (show $ fmap LG.toEdge  (origRootNewEdge' : newRerootEdges)) ++ " root flip: " ++ (show $ fmap LG.toEdge notFromRootInEdges))
  redirectEdgesDamped newGraph' newRerootEdges

-- | redirectEdgesDamped takes a graph and a list of vertices
-- the oout edges from the vertices are checked for correct direction
-- if edges already correct--then returns, else recurses to descendant vertices of edge
redirectEdgesDamped :: LG.Gr a b -> [LG.LEdge b] -> LG.Gr a b
redirectEdgesDamped inGraph edgeList =
  trace ("RDI:" ++ (show $ fmap LG.toEdge edgeList)) (
  if null edgeList then inGraph
  else 
     let nodeList = fmap snd3 edgeList
         outEdgeList = concatMap (LG.out inGraph) nodeList
         (flippedEdges, needToBeFlippedEdges) = unzip $ flipIfNeeded  nodeList outEdgeList
         newGraph = LG.insEdges flippedEdges $ LG.delLEdges needToBeFlippedEdges inGraph
      in
      redirectEdgesDamped newGraph flippedEdges 
  )

-- | flipIfNeeded takes a list of nodes and an edge and if the first field of edge is in node list
-- flips the vertex.  Returns flipeped and original edge
flipIfNeeded ::[Int] -> [LG.LEdge b] -> [(LG.LEdge b, LG.LEdge b)]
flipIfNeeded nodeList inEdgeList =
  trace ("Flip:" ++ (show nodeList) ++ " edge " ++ (show $ fmap LG.toEdge inEdgeList)) (
  if null inEdgeList then []
  else 
    let firstEdge@(a,b,c) = head inEdgeList
    in
    if b `elem` nodeList then ((b,a,c),firstEdge) : flipIfNeeded nodeList (tail inEdgeList)
    else flipIfNeeded nodeList (tail inEdgeList)
  )
-}

-- | rerootTree takes a graph and reroots based on a vertex index (usually leaf outgroup)
--   if input is a forest then only roots the component that contains the vertex wil be rerooted
--   unclear how will effect network edges--will need to verify that does not create cycles
--   multi-rooted components (as opposed to forests) are unaffected with trace warning thrown
--   after checking for existing root and multiroots, should be O(n) where 'n is the length
--   of the path between the old and new root
rerootTree :: (Eq b) => Int -> LG.Gr a b -> LG.Gr a b
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
    in

    -- check if new outtaxon has a parent--shouldn't happen-but could if its an internal node reroot
    if null parentNewRootList || (True `elem` parentRootList) then inGraph 
                                                              else (if null componentWithOutgroup then error ("Outgroup index " ++ show rerootIndex ++ " not found in graph")
    else
        --reroot component with new outtaxon
        let componentWithNewOutgroup = snd $ head componentWithOutgroup
            (_, originalRootList) =  unzip $ filter ((==True).fst) $ zip (fmap (LG.isRoot inGraph) componentWithNewOutgroup) componentWithNewOutgroup
            numRoots = length originalRootList
            orginalRoot = head originalRootList
            originalRootEdges = LG.out inGraph orginalRoot

        in
        
        -- check if outgroup in a multirooted component
        if numRoots > 1 then trace ("Warning: Ignoring reroot of multi-rooted component") inGraph
        else
          --reroot graph safely automatically will only affect the component with the outgroup
          -- delete old root edge and create two new edges from oringal root node.
          -- keep orignl root node and delte/crete new edges when they are encounterd
          --trace ("Moving root from " ++ (show orginalRoot) ++ " to " ++  (show rerootIndex)) (
          let leftChildEdge = (orginalRoot, rerootIndex, LG.edgeLabel $ head originalRootEdges)
              rightChildEdge = (orginalRoot, fst3 newRootOrigEdge, LG.edgeLabel $ last originalRootEdges)

              --  this assumes 2 children of old root -- shouled be correct as Phylogenetic Graph
              newEdgeOnOldRoot = if (length originalRootEdges) /= 2 then error ("Number of root out edges /= 1 in rerootGraph")
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
        -- retursn list of edges that had to be flipped
        edgesToFlip = getToFlipEdges v childEdges
        flippedEdges = fmap LG.flipLEdge edgesToFlip
        newGraph = LG.insEdges flippedEdges $ LG.delLEdges edgesToFlip inGraph
    in
    --trace ("+") (
    -- edge terminates in leaf or edges in correct orientation
    if null childEdges  || null edgesToFlip then preTraverseAndFlipEdges (tail inEdgelist) inGraph
    -- edge needs to be reversed to follow through its children from a new graph
    else preTraverseAndFlipEdges (flippedEdges ++ (tail inEdgelist)) newGraph
    -- )

-- | getToFlipEdges takes an index and ceck edge list
-- and cretes new list of edges that need to be flipped
getToFlipEdges ::  LG.Node -> [LG.LEdge b] -> [LG.LEdge b]
getToFlipEdges parentNodeIndex inEdgeList =
  if null inEdgeList then []
  else
    let firstEdge@(u,_,_) = head inEdgeList
    in
    if parentNodeIndex /= u then firstEdge : getToFlipEdges parentNodeIndex (tail inEdgeList)
    else getToFlipEdges parentNodeIndex (tail inEdgeList)


-- | contractOneOneEdges removes indegree 1 outdegree 1 nodes and its edges creating a new edge
-- connecting the nodes on either side
contractOneOneEdges :: LG.Gr a b -> LG.Gr a b
contractOneOneEdges inGraph =
  if LG.isEmpty inGraph then LG.empty
  else
    let nodeList = LG.nodes inGraph
        inOutEdgeList = zip (fmap (LG.getInOutEdges inGraph) nodeList) nodeList
        (edgesToDelete, edgesToAdd, nodesToDelete) = getContractGraphEdits inOutEdgeList ([],[],[])
        newGraph = LG.insEdges edgesToAdd $ LG.delNodes nodesToDelete $ LG.delLEdges edgesToDelete inGraph
    in
    newGraph

-- | getContractGraphEdits takes a series of pair of indegree outdegree and nodes and
-- returns list of graph edits to contract edges
getContractGraphEdits :: [(([LG.LEdge b], [LG.LEdge b]), LG.Node)] -> ([LG.LEdge b],[LG.LEdge b],[LG.Node]) -> ([LG.LEdge b],[LG.LEdge b],[LG.Node])
getContractGraphEdits inEdgeNodeList curEdits@(edgesToDelete, edgesToAdd, nodesToDelete) =
  if null inEdgeNodeList then curEdits
  else
    let ((firstInEdges, firstOutEdges), firstNode) = head inEdgeNodeList
    in
    if  (length firstInEdges, length firstOutEdges) /= (1,1) then getContractGraphEdits (tail inEdgeNodeList) curEdits
    else
      let newEdge = (fst3 $ head firstInEdges, snd3 $ head firstOutEdges, thd3 $ head firstInEdges)
      in
      getContractGraphEdits (tail inEdgeNodeList) (firstInEdges ++ firstOutEdges ++ edgesToDelete, newEdge : edgesToAdd, firstNode : nodesToDelete)


-- | generateDisplayTrees takes a graph list and recursively generates
-- a list of trees created by progresively resolving each network vertex into a tree vertex
-- in each input graph
-- creating up to 2**m (m network vertices) trees.  This only resolves the indegree
-- edges.  Indegree 1 outdegree 1 edges ARE NOT contracted when created or encountered.
-- call -> generateDisplayTrees  [startGraph] []
-- the second and third args contain graphs that need more work and graphs that are done (ie trees)
generateDisplayTrees :: [LG.Gr a b] -> [LG.Gr a b] -> [LG.Gr a b]
generateDisplayTrees curGraphList treeList  =
  if null curGraphList then treeList
  else
    let firstGraph = head curGraphList
    in
      if LG.isEmpty firstGraph then []
      else
        let nodeList = LG.labNodes firstGraph
            inNetEdgeList = filter ((>1).length) $ fmap (LG.inn firstGraph) $ fmap fst nodeList
        in
        if null inNetEdgeList then generateDisplayTrees (tail curGraphList) (firstGraph : treeList)
        else
          let newGraphList = splitGraphListFromNode inNetEdgeList [firstGraph]
          in
          generateDisplayTrees (newGraphList ++ (tail curGraphList)) treeList

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


-- | verifyTimeConsistency take a SimpleGraph and checks for time consistency
-- of network nodes to verify network nodes are not definately not coeval
-- basically, each set of network edges (usually 2) create a split in the set of nodes (and edges)
-- to those 'before' (ie leading to) and after (ie leadgin from) the heads of the network edges
-- thse sets can be generated for each vertex (cost n for traversal) and these sets of nodes
-- must be compatible A intersect B = A|B|Empty for teh graph to be time consistant.
-- would be edge based to check before and a network edge were to be added

--Logic wrong--may have to look at each pair of in-nodes to network edge

verifyTimeConsistency :: SimpleGraph -> SimpleGraph
verifyTimeConsistency inGraph =
   if LG.isEmpty inGraph then error ("Input Graph is empty in verifyTimeConsistency")
   else inGraph

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
selectPhylogeneticGraph inArgs seed selectArgList curGraphs =
    if null curGraphs then []
    else
        let fstArgList = fmap (fmap C.toLower . fst) inArgs
            sndArgList = fmap (fmap C.toLower . snd) inArgs
            lcArgList = zip fstArgList sndArgList
            checkCommandList = U.checkCommandArgs "select" fstArgList selectArgList
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

                        -- nonZeroEdgeLists for graphs
                        nonZeroEdgeListGraphPairList = fmap getNonZeroEdges curGraphs

                        -- keep only unique graphs based on non-zero edges
                        uniqueGraphList = L.sortOn snd6 $ getBVUniqPhylogeneticGraph True curGraphs -- getTopoUniqPhylogeneticGraph True curGraphs -- ngetUniqueGraphs nonZeroEdgeListGraphPairList []
                    in
                    if doUnique then take (fromJust numberToKeep) uniqueGraphList
                    else if doBest then take (fromJust numberToKeep) $ filter ((== minGraphCost).snd6) uniqueGraphList
                    else if doRandom then
                         let randList = head $ shuffleInt seed 1 [0..(length curGraphs - 1)]
                             (_, shuffledGraphs) = unzip $ L.sortOn fst $ zip randList curGraphs
                         in
                         take (fromJust numberToKeep) shuffledGraphs
                    -- default is best and unique
                    else
                        filter ((== minGraphCost).snd6) uniqueGraphList


-- | could use FGL '==' ?
-- | getUniqueGraphs takes each pair of non-zero edges and conpares them--if equal not added to list
getUniqueGraphs :: [([LG.LEdge EdgeInfo], PhylogeneticGraph)] -> [([LG.LEdge EdgeInfo], PhylogeneticGraph)]  -> [PhylogeneticGraph]
getUniqueGraphs inGraphPairList currentUniquePairs =
    if null inGraphPairList then fmap snd currentUniquePairs
    else
        let firstPair@(firstEdges, _) = head inGraphPairList
        in
        if null currentUniquePairs then getUniqueGraphs (tail inGraphPairList) [firstPair]
        else
            let equalList = filter (== True) $ fmap ((== firstEdges) . fst) currentUniquePairs
            in
            if null equalList then getUniqueGraphs (tail inGraphPairList) (firstPair : currentUniquePairs)
            else getUniqueGraphs (tail inGraphPairList) currentUniquePairs


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
makeIAPrelimFromFinal (index, label) =
  let labData = vertData label
      newLabData = fmap (fmap f) labData
  in
  (index, label {vertData = newLabData})
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
makeIAFinalFromPrelim (index, label) =
  let labData = vertData label
      newLabData = fmap (fmap f) labData
  in
  (index, label {vertData = newLabData})
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
          edgesG1 = LG.labEdges g1
          edgesG2 = LG.labEdges g2
      in
      if bvNodesG1 == bvNodesG2 then True
      else False

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



