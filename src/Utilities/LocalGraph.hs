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
getFENLocal inText = GFU.forestEnhancedNewickStringList2FGLList inText


-- | readDotLocal calls GrapvViz function to allow for substitution later
readDotLocal :: String -> IO (Utilities.LocalGraph.DotGraph Utilities.LocalGraph.Node)
readDotLocal fileName = GVIO.readDotFile fileName

-- | dotToGraph local mapo dor 
dotToGraph ::  (Utilities.LocalGraph.DotGraph Utilities.LocalGraph.Node) -> (Utilities.LocalGraph.Gr Attributes Attributes)
dotToGraph dotGraphList = GV.dotToGraph dotGraphList

-- | hGetDotLocal calls hGetDot from GraphVoz
hGetDotLocal :: Handle -> IO (Utilities.LocalGraph.DotGraph Utilities.LocalGraph.Node)
hGetDotLocal inFileHandle = GVIO.hGetDot inFileHandle

-- | fglToPrettyString calls prettify from FGL
fglToPrettyString :: (Show a, Show b) => P.Gr a b -> String
fglToPrettyString inGraph = G.prettify inGraph

-- Wrapper functions for fgl so could swap out later if want to

-- | maps to isEmpty
isEmpty :: Gr a b -> Bool
isEmpty inGraph = G.isEmpty inGraph

-- | maps to empty
empty :: Gr a b
empty = G.empty

-- | gelem is a node in a graph
gelem :: Node -> Gr a b -> Bool
gelem inNode inGraph = G.gelem inNode inGraph 

-- | maps to labNodes
labNodes :: Gr a b -> [LNode a] 
labNodes inGraph = G.labNodes inGraph

-- | maps to labEdges
labEdges :: Gr a b -> [LEdge b] 
labEdges inGraph = G.labEdges inGraph

-- | toEdge removes label from edge
toEdge :: LEdge b -> Edge
toEdge inEdge = G.toEdge inEdge 

-- | toLEdge adds a label to an edge
toLEdge :: Edge -> b -> LEdge b 
toLEdge inEdge inLabel = G.toLEdge inEdge inLabel

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
inn inGraph inNode = G.inn inGraph inNode

-- | lab returns label of node as Maybe
lab :: Gr a b -> Node -> Maybe a 
lab inGraph inNode = G.lab inGraph inNode

-- | out-bound edge list from node, maps to out
out :: Gr a b -> Node -> [LEdge b]
out inGraph inNode = G.out inGraph inNode

-- | parents of unlabelled node
parents :: Gr a b -> Node -> [Node]
parents inGraph inNode = fmap fst3 $ G.inn inGraph inNode

-- | descendants of unlabelled node
descendants :: Gr a b -> Node -> [Node]
descendants inGraph inNode = fmap snd3 $ G.out inGraph inNode

-- | labDescendants labelled descendents of labelled node
labDescendants :: (Eq a) => Gr a b -> LNode a -> [LNode a]
labDescendants inGraph inNode = 
    let nodeList = fmap snd3 $ G.out inGraph (fst inNode)
        maybeLabelList = fmap (lab inGraph) nodeList
        nothingList = filter (== Nothing) maybeLabelList
        labelList = fmap fromJust maybeLabelList
        labNodeList = zip nodeList labelList
    in
    if (not $ null nothingList) then error "UNlabeled nodes in labDescendants" 
    else labNodeList

-- | takes a graph and node and returns pair of inbound and noutbound labelled edges 
getInOutEdges :: Gr a b -> Node -> ([LEdge b], [LEdge b])
getInOutEdges inGraph inNode = (inn inGraph inNode, out inGraph inNode )

-- | nodes returns list of unlabbeled nodes, maps to nodes
nodes ::  Gr a b -> [Node]
nodes inGraph = G.nodes inGraph

-- | edges returns list of unlabbeled nodes, maps to nodes
edges ::  Gr a b -> [Edge]
edges inGraph = G.edges inGraph

-- | insEdges inserts a list of labelled edges into a graph
insEdges :: [LEdge b] -> Gr a b -> Gr a b
insEdges edgeList inGraph = G.insEdges edgeList inGraph

-- | insEdge inserts a labelled edge into a graph
insEdge :: LEdge b -> Gr a b -> Gr a b
insEdge inEdge inGraph = G.insEdge inEdge inGraph

-- | delLEdges delete a labelled edge from a graph
-- wrapps around delEdge
delLEdge :: LEdge b -> Gr a b -> Gr a b
delLEdge inEdge inGraph = G.delEdge (G.toEdge inEdge) inGraph

-- | delLEdge deletes a list of labelled edges from a graph
-- wrapps around delEdges
delLEdges :: [LEdge b] -> Gr a b -> Gr a b
delLEdges inEdgeList inGraph = G.delEdges (fmap G.toEdge inEdgeList) inGraph

-- | insNode inserts a labelled  node into a graph
insNode :: LNode a -> Gr a b -> Gr a b
insNode inNode inGraph = G.insNode inNode inGraph

-- | insNodes inserts multiple labelled nodes into a graph
insNodes :: [LNode a] -> Gr a b -> Gr a b
insNodes inNodeList inGraph = G.insNodes inNodeList inGraph

-- | delLNode deletes a labelled node from a graph
-- NB  I beleive removes any edges involving this node
delLNode :: LNode a -> Gr a b -> Gr a b
delLNode inNode inGraph = G.delNode (fst inNode) inGraph

-- | delNode deletes an unlabelled node from a graph
-- NB  I beleive removes any edges involving this node
delNode :: Node -> Gr a b -> Gr a b
delNode inNode inGraph = G.delNode inNode inGraph

-- | delNodes deletes a list of unlabelled nodes from a graph
-- NB  I beleive removes any edges involving these nodes
delNodes :: [Node] -> Gr a b -> Gr a b
delNodes inNodeList inGraph = G.delNodes inNodeList inGraph

-- | mkGraph creates a greaph from list of nodes and list of edges
mkGraph :: [LNode a] -> [LEdge b] -> Gr a b
mkGraph nodeList edgeList = G.mkGraph nodeList edgeList

-- | components list of list of nodes (graphalyze can return graph list)
components :: Gr a b -> [[Node]]
components inGraph = DFS.components inGraph

-- | noComponents returns number of components
noComponents :: Gr a b -> Int
noComponents inGraph = DFS.noComponents inGraph

-- | isLeaf checks if node is root 
isLeaf  :: Gr a b -> Node -> Bool
isLeaf  inGraph inNode = (G.outdeg inGraph inNode) == 0

-- | isNetworkNode checks if node is network node 
isNetworkNode  :: Gr a b -> Node -> Bool
isNetworkNode  inGraph inNode = (((G.indeg inGraph inNode) > 1) && ((G.outdeg inGraph inNode) > 0))

-- | isTreeNode checks if node is network node 
isTreeNode  :: Gr a b -> Node -> Bool
isTreeNode  inGraph inNode = (((G.indeg inGraph inNode) == 1) && ((G.outdeg inGraph inNode) > 0))

-- getRoots returns list of graph roots (labelled)
getRoots :: Gr a b -> [LNode a]
getRoots inGraph = 
    if isEmpty inGraph then []
    else 
        let nodeList =  labNodes inGraph
            rootBoolList = fmap (isRoot inGraph) (fmap fst nodeList)
            pairList = zip rootBoolList nodeList
            rootPairList =  filter ((==True).fst) pairList
            rootList = fmap snd rootPairList
        in
        rootList 

-- | isRoot checks if node is root 
isRoot :: Gr a b -> Node-> Bool
isRoot inGraph inNode = 
    if not $ G.gelem inNode inGraph then False
    else (G.indeg inGraph inNode) == 0

-- | pre returns list of nodes linking to a node 
pre :: Gr a b -> Node -> [Node]
pre inGraph inNode = G.pre inGraph inNode 

-- | suc returns list of nodes linking from a node 
suc :: Gr a b -> Node -> [Node]
suc inGraph inNode = G.suc inGraph inNode 

-- | edgeLabel returns label of edge
edgeLabel :: LEdge b -> b
edgeLabel inEdge = G.edgeLabel inEdge

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
        nodeTripleList = zip3 degOutList degInList (labNodes inGraph)
        treeTripleList = filter ((==1).fst3 ) $ filter ((>0).snd3 ) nodeTripleList
        (_, _, treeVertexList) = unzip3 treeTripleList

         -- network nodes
        networkTripleList = filter ((>1).fst3 ) $ filter ((>0).snd3 ) nodeTripleList
        (_, _, networkVertexList) = unzip3 networkTripleList
    in
    (rootList, leafList, treeVertexList, networkVertexList)

-- | pretty prints graph to String
prettify :: (Show a, Show b) => Gr a b -> String
prettify inGraph = G.prettify inGraph

-- | pathToRoot takes a greaph and a vertex and reurns a pair of lists 
-- of vertices and edges to root(s)
pathToRoot :: (Eq a, Show a) => Gr a b -> LNode a -> ([LNode a], [LEdge b])
pathToRoot inGraph inNode =
    if G.isEmpty inGraph then error "Empty graph in pathToRoot"
    else nodesAndEdgesBefore inGraph ([],[]) [inNode]

-- | nodesAndEdgesBefore takes a graph and list of nodes to get list of nodes
-- and edges 'before' in the sense of leading to--ie between root and
-- (not including)) that node
-- call with ([], [])
nodesAndEdgesBefore :: (Eq a, Show a) => Gr a b -> ([LNode a], [LEdge b]) -> [LNode a] -> ([LNode a], [LEdge b])
nodesAndEdgesBefore inGraph curResults@(curNodes, curEdges) inNodeList =
  if G.isEmpty inGraph then error ("Input Graph is empty in nodesAndEdgesBefore")
  else if null inNodeList then curResults
  else 
    let intoEdgeList = inn inGraph (fst $ head inNodeList)
        intoNodeList = fmap fst3 intoEdgeList
        labelMaybeList = fmap (lab inGraph) intoNodeList
        labelList = fmap fromJust labelMaybeList
        intoLabNodeList = zip intoNodeList labelList
    in
    if Nothing `elem` labelMaybeList then error ("Empty node label in nodesAndEdgesBefore" ++ show intoLabNodeList)
    else nodesAndEdgesBefore inGraph (intoLabNodeList ++ curNodes, intoEdgeList ++ curEdges) (intoLabNodeList ++ (tail inNodeList)) 

-- | nodesAndEdgesAfter takes a graph and list of nodes to get list of nodes
-- and edges 'after' in the sense of leading from-ie between (not including)) that node
-- and all the way to any leaves is connects to.
-- call with ([], [])
nodesAndEdgesAfter :: (Eq a, Show a) => Gr a b -> ([LNode a], [LEdge b]) -> [LNode a] -> ([LNode a], [LEdge b])
nodesAndEdgesAfter inGraph curResults@(curNodes, curEdges) inNodeList =
  if G.isEmpty inGraph then error ("Input Graph is empty in nodesAndEdgesAfter")
  else if null inNodeList then curResults
  else 
    let fromEdgeList = out inGraph (fst $ head inNodeList)
        fromNodeList = fmap snd3 fromEdgeList
        labelMaybeList = fmap (lab inGraph) fromNodeList
        labelList = fmap fromJust labelMaybeList
        fromLabNodeList = zip fromNodeList labelList
    in
    if Nothing `elem` labelMaybeList then error ("Empty node label in nodesAndEdgesAfter" ++ show fromLabNodeList)
    else nodesAndEdgesAfter inGraph (fromLabNodeList ++ curNodes, fromEdgeList ++ curEdges) (fromLabNodeList ++ (tail inNodeList)) 

