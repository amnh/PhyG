{- |
Module      :  Utilities.LocalGraph.hs
Description :  Module specifying graph types and functionality
                This is for indirection so can change underlying graph library
                without  polutting the rest of the code
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

{-# LANGUAGE ScopedTypeVariables #-}

module Utilities.LocalGraph  where

import qualified Data.Graph.Inductive.PatriciaTree as P
import qualified Data.Graph.Inductive.Query.DFS as DFS
import qualified GraphFormatUtilities              as GFU
--import qualified Data.Text as T
import qualified Data.Text.Lazy as T
import           Data.GraphViz                     as GV
--import           Data.GraphViz.Attributes.Complete (Attribute (Label),
                                                    --Label (..))
import           Data.GraphViz.Commands.IO         as GVIO
import qualified Data.Graph.Inductive.Graph        as G
import           System.IO
import GeneralUtilities
import Data.Maybe
import Debug.Trace
import qualified Data.List as L



-- | Gr local graph definition using FGL
type Gr a b = P.Gr a b
type Node = G.Node
type LNode a = G.LNode a
type DotGraph = GV.DotGraph
type Edge = G.Edge
type LEdge b = G.LEdge b


-- | getFENLocal maps to forestEnhancedNewickStringList2FGLList in GraphFormatUtilities
-- to allow for potnetial swapping FGL graph backend
-- requires leading and trailing space and newlines to be removed
getFENLocal :: T.Text -> [Utilities.LocalGraph.Gr T.Text Double]
getFENLocal = GFU.forestEnhancedNewickStringList2FGLList


-- | readDotLocal calls GrapvViz function to allow for substitution later
readDotLocal :: String -> IO (Utilities.LocalGraph.DotGraph Utilities.LocalGraph.Node)
readDotLocal = GVIO.readDotFile

-- | dotToGraph local mapo dor 
dotToGraph ::  Utilities.LocalGraph.DotGraph Utilities.LocalGraph.Node -> Utilities.LocalGraph.Gr Attributes Attributes
dotToGraph = GV.dotToGraph

-- | hGetDotLocal calls hGetDot from GraphVoz
hGetDotLocal :: Handle -> IO (Utilities.LocalGraph.DotGraph Utilities.LocalGraph.Node)
hGetDotLocal = GVIO.hGetDot

-- | fglToPrettyString calls prettify from FGL
fglToPrettyString :: (Show a, Show b) => P.Gr a b -> String
fglToPrettyString inGraph =
    if G.isEmpty inGraph then "Empty Graph"
    else G.prettify inGraph

-- Wrapper functions for fgl so could swap out later if want to

-- | maps to isEmpty
isEmpty :: Gr a b -> Bool
isEmpty = G.isEmpty

-- | maps to empty
empty :: Gr a b
empty = G.empty

-- | gelem is a node in a graph
gelem :: Node -> Gr a b -> Bool
gelem = G.gelem

-- | maps to labNodes
labNodes :: Gr a b -> [LNode a]
labNodes = G.labNodes

-- | maps to labEdges
labEdges :: Gr a b -> [LEdge b]
labEdges = G.labEdges

-- | toEdge removes label from edge
toEdge :: LEdge b -> Edge
toEdge = G.toEdge

-- | toLEdge adds a label to an edge
toLEdge :: Edge -> b -> LEdge b
toLEdge = G.toLEdge

-- | toLEdge' flipped version of toLEdge
toLEdge' :: b -> Edge -> LEdge b
toLEdge' inLabel inEdge = G.toLEdge inEdge inLabel

-- | maps to indeg
indeg :: Gr a b -> LNode a -> Int
indeg inGraph inLNode = G.indeg inGraph $ fst inLNode

-- | maps to outdeg
outdeg :: Gr a b -> LNode a -> Int
outdeg inGraph inLNode = G.outdeg inGraph $ fst inLNode

-- | getInOut takes a node and a graph and returns 
-- a triple (LNode a, indegree, outDegree)
getInOutDeg :: Gr a b -> LNode a -> (LNode a, Int, Int)
getInOutDeg inGraph inLNode =
    if isEmpty inGraph then error "Empty graph in getInOut"
    else
        let inDeg = indeg inGraph inLNode
            outDeg = outdeg inGraph inLNode
        in
        (inLNode, inDeg, outDeg)

-- | in-bound edge list to node, maps to inn
inn :: Gr a b -> Node -> [LEdge b]
inn = G.inn

-- | lab returns label of node as Maybe
lab :: Gr a b -> Node -> Maybe a
lab = G.lab

-- | out-bound edge list from node, maps to out
out :: Gr a b -> Node -> [LEdge b]
out = G.out

-- | parents of unlabelled node
parents :: Gr a b -> Node -> [Node]
parents inGraph inNode = fst3 <$> G.inn inGraph inNode

-- | descendants of unlabelled node
descendants :: Gr a b -> Node -> [Node]
descendants inGraph inNode = snd3 <$> G.out inGraph inNode

-- | labDescendants labelled descendents of labelled node
labDescendants :: (Eq a) => Gr a b -> LNode a -> [LNode a]
labDescendants inGraph inNode =
    let nodeList = snd3 <$> G.out inGraph (fst inNode)
        maybeLabelList = fmap (lab inGraph) nodeList
        nothingList = filter (== Nothing) maybeLabelList
        labelList = fmap fromJust maybeLabelList
        labNodeList = zip nodeList labelList
    in
    if not $ null nothingList then error "Unlabeled nodes in labDescendants"
    else labNodeList

-- | takes a graph and node and returns pair of inbound and noutbound labelled edges 
getInOutEdges :: Gr a b -> Node -> ([LEdge b], [LEdge b])
getInOutEdges inGraph inNode = (inn inGraph inNode, out inGraph inNode )

-- | nodes returns list of unlabbeled nodes, maps to nodes
nodes ::  Gr a b -> [Node]
nodes = G.nodes

-- | edges returns list of unlabeled nodes, maps to nodes
edges ::  Gr a b -> [Edge]
edges = G.edges

-- | insEdges inserts a list of labelled edges into a graph
insEdges :: [LEdge b] -> Gr a b -> Gr a b
insEdges = G.insEdges

-- | insEdge inserts a labelled edge into a graph
insEdge :: LEdge b -> Gr a b -> Gr a b
insEdge = G.insEdge

-- | delLEdges delete a labelled edge from a graph
-- wrapps around delEdge
delLEdge :: LEdge b -> Gr a b -> Gr a b
delLEdge inEdge = G.delEdge (G.toEdge inEdge)

-- | delLEdge deletes a list of labelled edges from a graph
-- wrapps around delEdges
delLEdges :: [LEdge b] -> Gr a b -> Gr a b
delLEdges inEdgeList = G.delEdges (fmap G.toEdge inEdgeList)

-- | insNode inserts a labelled  node into a graph
insNode :: LNode a -> Gr a b -> Gr a b
insNode = G.insNode

-- | insNodes inserts multiple labelled nodes into a graph
insNodes :: [LNode a] -> Gr a b -> Gr a b
insNodes = G.insNodes

-- | delLNode deletes a labelled node from a graph
-- NB  Removes any edges involving this node
delLNode :: LNode a -> Gr a b -> Gr a b
delLNode inNode = G.delNode (fst inNode)

-- | delNode deletes an unlabelled node from a graph
-- NB  Removes any edges involving this node
delNode :: Node -> Gr a b -> Gr a b
delNode = G.delNode

-- | delNodes deletes a list of unlabelled nodes from a graph
-- NB  I beleive removes any edges involving these nodes
delNodes :: [Node] -> Gr a b -> Gr a b
delNodes = G.delNodes

-- | mkGraph creates a graph from list of nodes and list of edges
mkGraph :: [LNode a] -> [LEdge b] -> Gr a b
mkGraph = G.mkGraph


-- | mkGraphPair creates a graph from pair of list of nodes and list of edges
mkGraphPair ::  ([LNode a], [LEdge b]) -> Gr a b
mkGraphPair (nodeList, edgeList) = G.mkGraph nodeList edgeList

-- | components list of list of nodes (graphalyze can return graph list)
components :: Gr a b -> [[Node]]
components = DFS.components

-- | componentGraphs takes a graph and returns its compnent gaphs
componentGraphs :: Gr a b -> [Gr a b]
componentGraphs inGraph =
    if isEmpty inGraph then []
    else
        let componentNodeListList = components inGraph
            labComponentNodeListList = fmap (fmap (labelNode inGraph)) componentNodeListList
            edgeListListList = fmap (fmap (inn inGraph)) componentNodeListList
            componentEdgeList = fmap concat edgeListListList
            componentGraphList = zipWith mkGraph labComponentNodeListList componentEdgeList
        in
        if length componentNodeListList == 1 then [inGraph]
        else componentGraphList


-- | labelNode uses lab but checks for Nothing and returns labelled node
labelNode :: Gr a b -> Node -> LNode a
labelNode inGraph inNode =
    if isEmpty inGraph then error "Empty graph for label source"
    else
        let label = lab inGraph inNode
        in
        if isNothing label then error ("No label for node " ++ show inNode)
        else
            (inNode, fromJust label)


-- | noComponents returns number of components
noComponents :: Gr a b -> Int
noComponents = DFS.noComponents

-- | isLeaf checks if node is root 
isLeaf  :: Gr a b -> Node -> Bool
isLeaf  inGraph inNode = G.outdeg inGraph inNode == 0

-- | isNetworkNode checks if node is network node 
isNetworkNode  :: Gr a b -> Node -> Bool
isNetworkNode  inGraph inNode = (G.indeg inGraph inNode > 1) && (G.outdeg inGraph inNode > 0)

-- | isTreeNode checks if node is network node 
isTreeNode  :: Gr a b -> Node -> Bool
isTreeNode  inGraph inNode = (G.indeg inGraph inNode == 1) && (G.outdeg inGraph inNode > 0)

-- getRoots returns list of graph roots (labelled)
getRoots :: Gr a b -> [LNode a]
getRoots inGraph =
    if isEmpty inGraph then []
    else
        let nodeList =  labNodes inGraph
            rootBoolList = fmap (isRoot inGraph . fst) nodeList
            pairList = zip rootBoolList nodeList
            rootPairList =  filter ((==True).fst) pairList
            rootList = fmap snd rootPairList
        in
        rootList

-- | isRoot checks if node is root 
isRoot :: Gr a b -> Node-> Bool
isRoot inGraph inNode =
    G.gelem inNode inGraph && (G.indeg inGraph inNode == 0)

-- | pre returns list of nodes linking to a node 
pre :: Gr a b -> Node -> [Node]
pre = G.pre

-- | suc returns list of nodes linking from a node 
suc :: Gr a b -> Node -> [Node]
suc = G.suc

-- | edgeLabel returns label of edge
edgeLabel :: LEdge b -> b
edgeLabel = G.edgeLabel

-- | getOtherVertex retuns the edge vertex /= index
getOtherVertex :: LEdge b -> Int -> Int
getOtherVertex (u,v,_) index = if u == index then v else u

-- | flipEdge flips orientation of unlabelled edge
flipEdge :: Edge -> Edge
flipEdge (u,v) = (v,u)

-- flipLEdge flips orientation of labelled edge
flipLEdge :: LEdge b ->LEdge b
flipLEdge (u,v,w) = (v,u,w)


 -- | splitVertexList splits the vertices of a graph into ([root], [leaf], [tree], [network])
splitVertexList ::  Gr a b -> ([LNode a], [LNode a], [LNode a], [LNode a])
splitVertexList inGraph =
  if G.isEmpty inGraph then ([],[],[],[])
  else
    let -- leaves
        degOutList = outdeg inGraph <$> labNodes inGraph
        newNodePair = zip degOutList (labNodes inGraph)
        leafPairList = filter ((==0).fst ) newNodePair
        (_, leafList) = unzip leafPairList

        -- roots
        degInList = indeg inGraph <$> labNodes inGraph
        newRootPair = zip degInList (labNodes inGraph)
        rootPairList = filter ((==0).fst ) newRootPair
        (_, rootList) = unzip rootPairList

        -- tree nodes
        nodeTripleList = zip3 degInList degOutList (labNodes inGraph)
        treeTripleList = filter ((==1).fst3 ) $ filter ((>0).snd3 ) nodeTripleList
        (_, _, treeVertexList) = unzip3 treeTripleList

         -- network nodes
        networkTripleList = filter ((>1).fst3 ) $ filter ((>0).snd3 ) nodeTripleList
        (_, _, networkVertexList) = unzip3 networkTripleList
    in
    (rootList, leafList, treeVertexList, networkVertexList)

-- | pretty prints graph to String
prettify :: (Show a, Show b) => Gr a b -> String
prettify inGraph =
    if G.isEmpty inGraph then "Empty Graph"
    else G.prettify inGraph

-- | pathToRoot takes a greaph and a vertex and reurns a pair of lists 
-- of vertices and edges to root(s) in order of encountering them to root
-- if a tree--not necessrily if network
pathToRoot :: (Eq a, Eq b) => Gr a b -> LNode a -> ([LNode a], [LEdge b])
pathToRoot inGraph inNode =
    if G.isEmpty inGraph then error "Empty graph in pathToRoot"
    else pathToRoot' inGraph [inNode] [] []

-- | pathToRoot with accumulators
pathToRoot' :: (Eq a, Eq b) => Gr a b -> [LNode a] -> [LNode a] -> [LEdge b] -> ([LNode a], [LEdge b])
pathToRoot' inGraph inNodeList curNodeList curEdgeList =
    if null inNodeList then (L.nub $ reverse curNodeList, L.nub $ reverse curEdgeList)
    else
        let inNode = head inNodeList
        in
        -- root would already be inlist of nodes visited
        if isRoot inGraph (fst inNode) then pathToRoot' inGraph (tail inNodeList) curNodeList curEdgeList
        else
            let inLEdges = inn inGraph (fst inNode)
                inNodes = fmap fst3 inLEdges
                inLabNodes = zip inNodes (fmap (fromJust . lab inGraph) inNodes)
            in
            pathToRoot' inGraph (L.nub $ inLabNodes ++ tail inNodeList) (inLabNodes ++ curNodeList) (inLEdges ++ curEdgeList)




-- | nodesAndEdgesBefore takes a graph and list of nodes to get list of nodes
-- and edges 'before' in the sense of leading to--ie between root and
-- (not including)) that node
-- call with ([], [])
nodesAndEdgesBefore :: (Eq a, Show a) => Gr a b -> ([LNode a], [LEdge b]) -> [LNode a] -> ([LNode a], [LEdge b])
nodesAndEdgesBefore inGraph curResults@(curNodes, curEdges) inNodeList
  | G.isEmpty inGraph = error "Input Graph is empty in nodesAndEdgesBefore"
  | null inNodeList = curResults
  | otherwise =
    let intoEdgeList = inn inGraph (fst $ head inNodeList)
        intoNodeList = fmap fst3 intoEdgeList
        labelMaybeList = fmap (lab inGraph) intoNodeList
        labelList = fmap fromJust labelMaybeList
        intoLabNodeList = zip intoNodeList labelList
    in
    if Nothing `elem` labelMaybeList then error ("Empty node label in nodesAndEdgesBefore" ++ show intoLabNodeList)
    else nodesAndEdgesBefore inGraph (intoLabNodeList ++ curNodes, intoEdgeList ++ curEdges) (intoLabNodeList ++ tail inNodeList)

-- | nodesAndEdgesAfter takes a graph and list of nodes to get list of nodes
-- and edges 'after' in the sense of leading from-ie between (not including)) that node
-- and all the way to any leaves is connects to.
-- call with ([], [])
nodesAndEdgesAfter :: (Eq a, Show a) => Gr a b -> ([LNode a], [LEdge b]) -> [LNode a] -> ([LNode a], [LEdge b])
nodesAndEdgesAfter inGraph curResults@(curNodes, curEdges) inNodeList
  | G.isEmpty inGraph = error "Input Graph is empty in nodesAndEdgesAfter"
  | null inNodeList = curResults
  | otherwise =
    let fromEdgeList = out inGraph (fst $ head inNodeList)
        fromNodeList = fmap snd3 fromEdgeList
        labelMaybeList = fmap (lab inGraph) fromNodeList
        labelList = fmap fromJust labelMaybeList
        fromLabNodeList = zip fromNodeList labelList
    in
    if Nothing `elem` labelMaybeList then error ("Empty node label in nodesAndEdgesAfter" ++ show fromLabNodeList)
    else nodesAndEdgesAfter inGraph (fromLabNodeList ++ curNodes, fromEdgeList ++ curEdges) (fromLabNodeList ++ tail inNodeList)

-- | contractIn1Out1Edges contracts indegree 1, outdegree 1, edges and removes the node in the middle
-- does one at a time and makes a graph and recurses
contractIn1Out1Edges :: (Show a, Show b) => Gr a b -> Gr a b
contractIn1Out1Edges inGraph =
    if G.isEmpty inGraph then G.empty
    else
        let inOutDeg = getInOutDeg inGraph <$> labNodes inGraph
            degree11VertexList = filter ((==1) . thd3) $ filter ((==1) . snd3) inOutDeg
        in
        trace ("11:" ++ show degree11VertexList) (
        if null degree11VertexList then inGraph
        else
                let nodeToDelete = fst3 $ head degree11VertexList
                    inEdgeToDelete  = head $ inn inGraph $ fst nodeToDelete
                    outEdgeToDelete = head $ out inGraph $ fst nodeToDelete
                    newEdgeToAdd    = (fst3 inEdgeToDelete, snd3 outEdgeToDelete, thd3 inEdgeToDelete)
                    newGraph = insEdge newEdgeToAdd $ delLNode nodeToDelete $ delLEdges [inEdgeToDelete, outEdgeToDelete] inGraph
                in
                trace ("Deleting Node " ++ show (fst nodeToDelete) ++ " " ++ show (inEdgeToDelete, outEdgeToDelete) ++ " inserting " ++ show  newEdgeToAdd)
                contractIn1Out1Edges newGraph
                )
