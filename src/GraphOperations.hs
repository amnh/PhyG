{- |
Module      :  GraphOperations.hs
Description :  Module specifying data type medians
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

{--
TODO:

  Parallelize  median2Vect
--}

module GraphOperations ( ladderizeGraph
                       , fullyLabelGraph
                       , verifyTimeConsistency
                       , rerootGraph
                       ) where

import           Debug.Trace
import qualified Data.Vector as V
import qualified DOWrapper as DOW
import Types
import qualified LocalGraph as LG
import qualified Data.Text.Lazy as T
import GeneralUtilities
import qualified LocalSequence as LS
import qualified GraphFormatUtilities as GFU

-- | ladderizeGraph is a wrapper around ladderizeGraph' to allow for mapping with 
-- local nodelist
ladderizeGraph :: SimpleGraph -> SimpleGraph
ladderizeGraph inGraph = ladderizeGraph' inGraph (LG.nodes inGraph)



-- | ladderize takes an input graph and ensures/creates nodes
-- such that all vertices are (indegree, outdegree) (0,>0), (1,2) (2,1) (1,0)
ladderizeGraph' :: SimpleGraph -> [LG.Node] -> SimpleGraph
ladderizeGraph' inGraph nodeList = 
    if LG.isEmpty inGraph then LG.empty
    else if null nodeList then inGraph
    else 
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
        --)

-- | resolveNode takes a graph and node and inbound edgelist and outbound edge list
-- and converts node to one off (indeg, outdeg) (0,1),(0,2),(1,2),(2,1),(1,0)
-- this only resolves a single nodes edges at a time and then returns new graph
-- when more hase to be done--that will occur on lultiple passes through nodes.
-- perhaps not the most efficient, but only done once per input graph
resolveNode :: SimpleGraph -> LG.Node -> ([LG.LEdge Double], [LG.LEdge Double]) -> (Int, Int) -> SimpleGraph
resolveNode inGraph curNode inOutPair@(inEdgeList, outEdgeList) (inNum, outNum) =
  if LG.isEmpty inGraph then LG.empty
  else 
    --trace ("Resolveing " ++ show curNode) (
    let numNodes = length $ LG.nodes inGraph
    in
    -- leaf or simple tree node in outdegree 
    if outNum == 0 || outNum == 1 then
      let first2Edges = take 2 inEdgeList
          newNode = (numNodes , T.pack $ ("HTU" ++ (show numNodes)))
          newEdge1 = (fst3 $ head first2Edges, numNodes, 0.0 :: Double)
          newEdge2 = (fst3 $ last first2Edges, numNodes, 0.0 :: Double)
          newEdge3 = (numNodes, curNode, 0.0 :: Double)
          newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
      in
      newGraph 
    -- root or simple network indegree node
    else if inNum == 0 || inNum == 2 then 
      let first2Edges = take 2 outEdgeList
          newNode = (numNodes , T.pack $ ("HTU" ++ (show numNodes)))
          newEdge1 = (numNodes, snd3 $ head first2Edges, 0.0 :: Double)
          newEdge2 = (numNodes, snd3 $ last first2Edges, 0.0 :: Double)
          newEdge3 = (curNode, numNodes, 0.0 :: Double)
          newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
      in 
      newGraph
    else error ("This can't happen in resolveNode in/out edge lists don't need to be resolved " ++ show inOutPair)
    --)

-- | verifyTimeConsistency take a SimpleGraph and checks for time consistency
-- of network nodes to verify network nodes are not definately not coeval
verifyTimeConsistency :: SimpleGraph -> SimpleGraph
verifyTimeConsistency inGraph =
   if LG.isEmpty inGraph then error ("Input Graph is empty in verifyTimeConsistency")
   else inGraph
   
   {-
   errorWithoutStackTrace ("Graph violates time consistency")
   -}

-- | fullyLabelGraph takes an unlabelled "simple' graph, performs post and preorder passes to 
-- fully label the graph and return a PhylogeenticGraph
fullyLabelGraph :: GlobalSettings -> ProcessedData -> SimpleGraph -> PhylogeneticGraph
fullyLabelGraph inGS inData inGraph = 
    if LG.isEmpty inGraph then (LG.empty, 0.0, V.empty, V.empty, inData)
    else 
        (inGraph, 0.0, V.empty, V.empty, inData)


-- | rerootGraph takes a graph and reroots based on a vertex index (usually leaf outgroup)
--   if input is a forest then only roots the component that contains the vertex wil be rerooted
--   unclear how will effect network edges--will need to verify that does not create cycles
--   multi-rooted components (as opposed to forests) are unaffected with trace warning thrown
--   after checking for existing root and multiroots, should be O(n) where 'n is the length
--   of the path between the old and new root 
rerootGraph :: Int -> SimpleGraph -> SimpleGraph
rerootGraph rerootIndex inGraph = 
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
    if null parentNewRootList then inGraph
    -- check if outgroup doesn't change rooting--ie its parent is a root somewhere
    else if True `elem` parentRootList then inGraph
    -- this can't happen but whatever....
    else if null componentWithOutgroup then error ("Outgroup index " ++ show rerootIndex ++ " not found in graph")
    else 
      --rooroot component with new outtaxon 
      let componentWithNewOutgroup = snd $ head componentWithOutgroup
          (_, originalRootList) =  unzip $ filter ((==True).fst) $ zip (fmap (LG.isRoot inGraph) componentWithNewOutgroup) componentWithNewOutgroup
          numRoots = length originalRootList
          orginalRoot = head originalRootList
          originalRootEdges = LG.out inGraph orginalRoot

      in
      {-These to prevent heads of empty lists 
      if null originalRootList then error "No orginal roots"
      else if null (LG.inn inGraph rerootIndex) then error ("Can't find parent of " ++ show rerootIndex)
      else if null componentWithOutgroup then error ("Can't findcomponent with " ++ show rerootIndex)
      else if null originalRootEdges then error ("Null original root edges ")
      else 
        -}
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
            newEdgeOnOldRoot = (snd3 $ head originalRootEdges, snd3 $ last originalRootEdges, thd3 $ head originalRootEdges)

            newRootEdges = [leftChildEdge, rightChildEdge, newEdgeOnOldRoot]
            newGraph = LG.insEdges newRootEdges $ LG.delLEdges (newRootOrigEdge : originalRootEdges) inGraph

            -- get edges that need reversing
            newGraph' = preTraverseAndFlipEdges newGraph [leftChildEdge,rightChildEdge]
        in
        --trace ("Deleting " ++ (show (newRootOrigEdge : originalRootEdges)) ++ "\nInserting " ++ (show newRootEdges)) (
        --trace ("In " ++ (GFU.showGraph inGraph) ++ "\nNew " ++  (GFU.showGraph newGraph) ++ "\nNewNew "  ++  (GFU.showGraph newGraph'))
        newGraph'
        --))

-- | preTraverseAndFlipEdges traverses graph from starting edge flipping edges as needed
-- when recursion its edges that don't need to be fliped then stops
-- assumes input edge is directed correctly
-- follows  traversal out "pre" order ish from edges
-- have to check edge orientatins--make sure thay haven't changed as graph was updated earlier
preTraverseAndFlipEdges :: (Eq b) => LG.Gr a b -> [LG.LEdge b] ->  LG.Gr a b
preTraverseAndFlipEdges inGraph inEdgelist  = 
  if null inEdgelist then inGraph
  else 
    let inEdge@(_,v,_) = head inEdgelist
        childEdges = (LG.out inGraph v) ++ (filter (/= inEdge) $ LG.inn inGraph v)
        -- retursn list of edges that had to be flipped
        edgesToFlip = getToFlipEdges v childEdges
        flippedEdges = fmap LG.flipLEdge edgesToFlip
        newGraph = LG.insEdges flippedEdges $ LG.delLEdges edgesToFlip inGraph
    in
    -- edge terminates in leaf or edges in correct orientation
    if null childEdges  || null edgesToFlip then preTraverseAndFlipEdges inGraph (tail inEdgelist) 
    -- edge needs to be reversed to follow through its children from a new graph
    else preTraverseAndFlipEdges newGraph (flippedEdges ++ (tail inEdgelist))

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
      
