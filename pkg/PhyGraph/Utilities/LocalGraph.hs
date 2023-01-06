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

import qualified Data.Graph.Inductive.Basic          as B
import qualified Data.Graph.Inductive.PatriciaTree   as P
import qualified Data.Graph.Inductive.Query.ArtPoint as AP
import qualified Data.Graph.Inductive.Query.BCC      as BCC
import qualified Data.Graph.Inductive.Query.BFS      as BFS
import qualified Data.Graph.Inductive.Query.DFS      as DFS
import qualified GraphFormatUtilities                as GFU
--import qualified Data.Text as T
import           Data.GraphViz                       as GV
import qualified Data.Text.Lazy                      as T
--import           Data.GraphViz.Attributes.Complete (Attribute (Label),
                                                    --Label (..))
import qualified Data.Graph.Inductive.Graph          as G
import           Data.GraphViz.Commands.IO           as GVIO
import           Control.Parallel.Strategies
import qualified Cyclic                              as C
import qualified Data.List                           as L
import qualified Data.Map                            as MAP
import           Data.Maybe
import qualified Data.Vector                         as V
import           GeneralUtilities
import qualified ParallelUtilities                   as PU
import           System.IO
import           Debug.Trace




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

-- | pretty prints graph to String
prettify :: (Show a, Show b) => Gr a b -> String
prettify inGraph =
    if G.isEmpty inGraph then "Empty Graph"
    else G.prettify inGraph

-- | prettyIndices prints graph to String only using indices
prettyIndices :: Gr a b -> String
prettyIndices inGraph =
    if G.isEmpty inGraph then "Empty Graph"
    else
        let nodeList = concat $ fmap (++ ",") $ fmap show $ nodes inGraph
            edgeList = concat $ fmap (++ ",") $ fmap show $ edges inGraph
        in
        nodeList ++ "\n" ++ edgeList

-- | prettyDot prints generic graph to generic dot format
prettyDot :: Gr a b -> String
prettyDot inGraph =
    if G.isEmpty inGraph then error ("Empty Graph")
    else
        let topPart = "digraph G {\n\trankdir = LR; node [ shape = rect];\n"
            nodeList = concatMap ("\t" ++ ) $ fmap (++ ";\n") $ fmap show $ nodes inGraph
            edgeList = concatMap ("\t" ++ ) $ fmap makeEdgeString $ edges inGraph
            endPart = "}\n"
        in
        topPart ++ nodeList ++ edgeList ++ endPart
    where makeEdgeString (a, b) = (show a) ++ " -> " ++ (show b) ++ ";\n"

-- these duplicate edge functions should be O(n log n) based on sort--rest linear


-- hasDuplicateEdgesNub hasDuplicateEdge but based on nub
hasDuplicateEdgesNub :: Gr a b -> Bool
hasDuplicateEdgesNub inGraph =
    if isEmpty inGraph then False
    else 
        let edgeList = fmap toEdge $ labEdges inGraph
            uniqueEdgeList = L.nub edgeList
        in
        not (length edgeList == length uniqueEdgeList)

-- hasDuplicateEdge checked for duplicate edges based on indices (not label)
hasDuplicateEdge :: Gr a b -> Bool
hasDuplicateEdge inGraph =
    if isEmpty inGraph then False
    else 
        let sortedEdges = L.sort $ fmap toEdge $ labEdges inGraph
            groupedEdges = L.group sortedEdges
            dupEdges = filter ((>1) . length ) groupedEdges
            --dupEdges = (fmap toEdge $ labEdges inGraph) L.\\ (L.nub $ fmap toEdge $ labEdges  inGraph)
        in 
        (not . null) dupEdges

-- | getDuplicateEdges retuns a list of edges that are duplicated
-- by indeices--no label comparison
-- can be used to delete the extra
getDuplicateEdges :: Gr a b -> [Edge]
getDuplicateEdges inGraph = 
    if isEmpty inGraph then []
    else 
        let sortedEdges = L.sort $ fmap toEdge $ labEdges inGraph
            groupedEdges = L.group sortedEdges
            dupEdges = filter ((>1) . length ) groupedEdges
            --dupEdges = (fmap toEdge $ labEdges inGraph) L.\\ (L.nub $ fmap toEdge $ labEdges  inGraph)
        in 
        fmap head dupEdges
        -- (fmap toEdge $ labEdges inGraph) L.\\ (L.nub $ fmap toEdge $ labEdges  inGraph)

-- | removeDuplicateEdges removes duplicate edges from graph
removeDuplicateEdges :: Gr a b -> Gr a b 
removeDuplicateEdges inGraph =
    if isEmpty inGraph then inGraph
    else 
        let dupEdges = getDuplicateEdges inGraph
        in
        if null dupEdges then inGraph
        else delEdges dupEdges inGraph

-- | hasTreeNodeWithAllNetworkChildren checks treenodes for all (should be 2) children that
-- are netork nodes
hasTreeNodeWithAllNetworkChildren :: Gr a b -> (Bool, [Node])
hasTreeNodeWithAllNetworkChildren inGraph =
    if isEmpty inGraph then (False, [])
    else 
        let (_, _, treeNodeList, _) = splitVertexList inGraph
            hasAllNetChildrenList = fmap (hasAllNetChildren inGraph) (fmap fst treeNodeList)
            nodesWithAllNetChildren = fmap (fst . fst) $ filter ((== True) .snd) $ zip treeNodeList hasAllNetChildrenList
        in
        ((not . null) nodesWithAllNetChildren, nodesWithAllNetChildren)

-- | hasAllNetChildren checks whether all (usually 2) childrenb of a vertex are network nodes
hasAllNetChildren :: Gr a b -> Node -> Bool
hasAllNetChildren inGraph inNode =
    let children = descendants inGraph inNode
        childVertNodes = filter (== True) $ fmap (isNetworkNode inGraph) children
    in
    length children == length childVertNodes

-- | removeTreeEdgeFromTreeNodeWithAllNetworkChildren takes a graph and removes the first edge (head) 
-- from each tree node with all network children, then contracts those edges and nodes, 
-- then reindexes -- but does not rename graph nodes
removeTreeEdgeFromTreeNodeWithAllNetworkChildren :: Gr a b -> Gr a b
removeTreeEdgeFromTreeNodeWithAllNetworkChildren inGraph = 
    let (toDo, nodesWithEdgesToDelete) = hasTreeNodeWithAllNetworkChildren inGraph
        outEdgesToDeleteList = fmap toEdge $ fmap head $ fmap (out inGraph) nodesWithEdgesToDelete
        newGraph = delEdges outEdgesToDeleteList inGraph
        newGraph' = reindexGraph $ contractIn1Out1Edges newGraph
    in
    if not toDo then inGraph
    else newGraph'

-- | hasChainedNetworkNodes checks if a graph has network nodes with at least one parent that is also a network node
hasChainedNetworkNodes :: Gr a b -> Bool
hasChainedNetworkNodes inGraph = 
    if isEmpty inGraph then False
    else
        let (_, _, _, netVertexList) = splitVertexList inGraph
            chainedNodeList = filter (== True) $ fmap (hasNetParent inGraph) $ fmap fst netVertexList
            netWithChildNetList = filter (== True) $ fmap (hasAllNetChildren inGraph) $ fmap fst netVertexList
        in
        if null netVertexList then False
        else (not . null) (chainedNodeList ++ netWithChildNetList)

-- | hasNetParent checks parent of node and retuens True if one or both are network nodes
hasNetParent :: Gr a b -> Node -> Bool
hasNetParent inGraph inNode =
    let parentList = parents inGraph inNode
        parentNetList = filter (== True) $ fmap (isNetworkNode inGraph) parentList
    in
    (not . null) parentNetList

-- | removeChainedNetworkNodes detectes and fixes (if possible) chained network edges
-- if 1 parent of network edge is tree node can be fixed by delete and contracting that node/edge 
-- else if both parent are netowrks--cannot be fixed and errors out
-- does NOT rename nodes since need vertex info on that--but are reindexed
removeChainedNetworkNodes :: (Show a, Show b) => Bool -> Gr a b -> Maybe (Gr a b)
removeChainedNetworkNodes showWarning inGraph = 
    if isEmpty inGraph then Just inGraph
    else
        let (_, _, _, netVertexList) = splitVertexList inGraph
            parentNetNodeList = fmap (hasNetParent inGraph) $ fmap fst netVertexList
            chainedNodeList =  fmap fst $ filter ((== True) . snd) $ zip netVertexList parentNetNodeList
            fixableChainedEdgeList = concatMap (getTreeEdgeParent inGraph) (fmap fst chainedNodeList)
            newGraph = delEdges fixableChainedEdgeList inGraph
            newGraph' = reindexGraph $ contractIn1Out1Edges newGraph
        in
        if null netVertexList then Just inGraph
        else if null chainedNodeList then Just inGraph
        else if null fixableChainedEdgeList then 
            trace ("Warning: Unfixable chained network nodes (both parent and child nodes are indegree > 1). Skipping graph")
            Nothing
        else 
            let warningString = if showWarning then "Warning: Chained network nodes (both parent and child nodes are indegree > 1), removing edges to tree node parents (this may affect graph cost)\n"
                            else ""
            in
            traceNoLF warningString
             -- : " ++ (show fixableChainedEdgeList))--  ++ "\n" ++ (prettyIndices inGraph))
            Just newGraph'

-- | getTreeEdgeParent gets the tree edge (as list) into a network node as opposed to the edge from a network parent
-- if both parents are netowrk nodes then returns []
getTreeEdgeParent :: Gr a b -> Node -> [Edge]
getTreeEdgeParent inGraph inNode =
    let parentList = parents inGraph inNode
        parentTreeList = fmap fst $ filter ((== False) . snd) $ zip parentList (fmap (isNetworkNode inGraph) parentList)
    in
    if null parentTreeList then []
    else [(head parentTreeList, inNode)]

-- Wrapper functions for fgl so could swap out later if want to

-- | maps to isEmpty
isEmpty :: Gr a b -> Bool
isEmpty = G.isEmpty

-- | maps to empty
empty :: Gr a b
empty = G.empty

-- | maps to equal
equal :: (Eq a, Eq b) => Gr a b -> Gr a b -> Bool
equal x y = G.equal x y

-- | gelem is a node in a graph
gelem :: Node -> Gr a b -> Bool
gelem = G.gelem

-- | maps to labNodes-- I beleive in vertex index order
labNodes :: Gr a b -> [LNode a]
labNodes = G.labNodes

-- | maps to labEdges
labEdges :: Gr a b -> [LEdge b]
labEdges = G.labEdges

-- | toEdge removes label from edge
toEdge :: LEdge b -> Edge
toEdge = G.toEdge

-- | toNode returns fst of LNode a
toNode :: LNode a -> Node
toNode inNode = fst inNode

-- | toLEdge adds a label to an edge
toLEdge :: Edge -> b -> LEdge b
toLEdge = G.toLEdge

-- | toLEdge' flipped version of toLEdge
toLEdge' :: b -> Edge -> LEdge b
toLEdge' inLabel inEdge = G.toLEdge inEdge inLabel

-- | deg mapes to fgl deg
deg :: Gr a b -> Node -> Int
deg inGraph inNode = G.deg inGraph inNode

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

-- | getInOutDegNoLabel takes a node and a graph and returns
-- a pair (indegree, outDegree)
getInOutDegNoLabel :: Gr a b -> Node -> (Int, Int)
getInOutDegNoLabel inGraph inNode =
    if isEmpty inGraph then error "Empty graph in getInOut"
    else
        let inDeg = length $ inn inGraph inNode
            outDeg = length $ out inGraph inNode
        in
        (inDeg, outDeg)

-- | in-bound edge list to node, maps to inn
inn :: Gr a b -> Node -> [LEdge b]
inn = G.inn

-- | lab returns label of node as Maybe
lab :: Gr a b -> Node -> Maybe a
lab = G.lab

-- | out-bound edge list from node, maps to out
out :: Gr a b -> Node -> [LEdge b]
out = G.out

-- | hasEdge  maps to fgl function returns True if graphs has directed edge between nodes
hasEdge :: Gr a b -> Edge -> Bool
hasEdge = G.hasEdge

-- | updateNodeLabel updated teh label infomatino on a  node
-- this is done by deleting that nmode and addingl back in to graph
-- when a node is deleted all edges inceident on it are also deleted
-- so they must be saved and added back
updateNodeLabel :: Gr a b -> Node  -> a -> Gr a b
updateNodeLabel inGraph inNodeIndex newLabel =
    if isEmpty inGraph then error "Empty graph in sisterLabNodes"
    else
        let incidentEdges = (inn inGraph inNodeIndex) ++ (out inGraph inNodeIndex)
        in
        insEdges incidentEdges $ insNode (inNodeIndex, newLabel) $ delNode inNodeIndex inGraph

-- | sisterLabNodes returns list of nodes that are "sister" ie share same parent
-- as input node
sisterLabNodes :: (Eq a) => Gr a b -> LNode a -> [LNode a]
sisterLabNodes inGraph inNode =
    if isEmpty inGraph then error "Empty graph in sisterLabNodes"
    else
        let parentNodeList = labParents inGraph (fst inNode)
            otherChildrenOfParentsList = (concatMap (labDescendants inGraph) parentNodeList) L.\\ [inNode]
        in
        -- shoudl not need nub for phylogenetic graphs but whatever
        L.nub otherChildrenOfParentsList

-- | parents of unlabelled node
parents :: Gr a b -> Node -> [Node]
parents inGraph inNode = fst3 <$> G.inn inGraph inNode

-- | grandParents of unlabelled node
grandParents :: Gr a b -> Node -> [Node]
grandParents inGraph inNode =
    let nodeParents = parents inGraph inNode
    in
    concatMap (parents inGraph) nodeParents

-- | sharedGrandParents sees if a node has common grandparents
sharedGrandParents :: Gr a b -> Node -> Bool
sharedGrandParents inGraph inNode =
    if isEmpty inGraph then error "Empty gaph in sharedGrandParents"
    else if isRoot inGraph inNode then False
    else
        let parentList = parents inGraph inNode
            grandParentLists = fmap (parents inGraph) parentList
            intersection = L.foldl1' L.intersect grandParentLists
        in
        --True
        if null intersection then False else True

-- | labParents returns the labelled parents of a node
labParents :: (Eq a) => Gr a b -> Node -> [LNode a]
labParents inGraph inNode =
    let parentNodeList = parents inGraph inNode
        parentLabelList = fmap (lab inGraph) parentNodeList
        hasNothing = Nothing `elem` parentLabelList
        parentLabelList' = fmap fromJust parentLabelList
    in
    if hasNothing then error "Unlabeled nodes in labParents"
    else zip parentNodeList parentLabelList'

-- | isPhylogeneticGraph checks various issues to see if 
-- there is wierdness in graph
isPhylogeneticGraph :: (Show a, Eq a, NFData a, Show b, Eq b) => Gr a b -> Bool
isPhylogeneticGraph inGraph =
    if isEmpty inGraph then False
    else 
        let nodeList = fmap fst $ labNodes inGraph
            indegreeList = fmap (inn inGraph) nodeList
            outdegreeList = fmap (out inGraph) nodeList
        in    
        if hasDuplicateEdgesNub inGraph then False
        else if length (getRoots inGraph) /= 1 then False
        else if (not . null) (getIsolatedNodes inGraph) then False
        else if (not . null) (filter ((> 2) . length) indegreeList) then False
        else if (not . null) (filter ((> 2) . length) outdegreeList) then False
        else if not (isGraphTimeConsistent inGraph) then False
        else True


-- | descendants of unlabelled node
descendants :: Gr a b -> Node -> [Node]
descendants inGraph inNode = snd3 <$> G.out inGraph inNode

-- | labDescendants labelled descendents of labelled node
labDescendants :: (Eq a) => Gr a b -> LNode a -> [LNode a]
labDescendants inGraph inNode =
    let nodeList = snd3 <$> G.out inGraph (fst inNode)
        maybeLabelList = fmap (lab inGraph) nodeList
        hasNothing = Nothing `elem` maybeLabelList
        labelList = fmap fromJust maybeLabelList
        labNodeList = zip nodeList labelList
    in
    if hasNothing then error "Unlabeled nodes in labDescendants"
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

-- | delEdges delete an unlabelled edge from a graph
-- wrapps around delEdge
delEdge :: Edge -> Gr a b -> Gr a b
delEdge inEdge = G.delEdge inEdge

-- | delLEdge deletes a list of unlabelled edges from a graph
-- wrapps around delEdges
delEdges :: [Edge] -> Gr a b -> Gr a b
delEdges inEdgeList = G.delEdges inEdgeList

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

-- | components list of list of nodes (G.Graphalyze can return graph list)
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

-- | labelNodeFlip uses lab but checks for Nothing and returns labelled node
labelNodeFlip :: Node -> Gr a b -> LNode a
labelNodeFlip = flip labelNode
    

-- | noComponents returns number of components
noComponents :: Gr a b -> Int
noComponents = DFS.noComponents

-- | isLeaf checks if node is root
isLeaf  :: Gr a b -> Node -> Bool
isLeaf  inGraph inNode = G.outdeg inGraph inNode == 0

-- | isOutDeg1Node checks if node has a single child
isOutDeg1Node :: Gr a b -> Node -> Bool
isOutDeg1Node inGraph inNode =  G.outdeg inGraph inNode == 1

-- | isNetworkNode checks if node is network node
isNetworkNode  :: Gr a b -> Node -> Bool
isNetworkNode  inGraph inNode = (G.indeg inGraph inNode > 1) && (G.outdeg inGraph inNode > 0)

-- | isNetworkLeaf checks if node is network node and a leaf--usually an error condition in phylogenetic networks
isNetworkLeaf  :: Gr a b -> Node -> Bool
isNetworkLeaf  inGraph inNode = (G.indeg inGraph inNode > 1) && (G.outdeg inGraph inNode == 0)


-- | - | isNetworkEdge checks if edge is network edge
isNetworkEdge  :: Gr a b -> Edge -> Bool
isNetworkEdge  inGraph inEdge = (G.indeg inGraph (snd inEdge) > 1) && (G.outdeg inGraph (snd inEdge) > 0)

-- | - | isNetworkLabEdge checks if edge is network edge
isNetworkLabEdge  :: Gr a b -> LEdge b -> Bool
isNetworkLabEdge  inGraph inEdge = (G.indeg inGraph (snd3 inEdge) > 1) && (G.outdeg inGraph (snd3 inEdge) > 0)

-- | labNetEdges takes a graph and returns list of network labelled edges
labNetEdges :: Gr a b -> [LEdge b]
labNetEdges inGraph =
    if isEmpty inGraph then error "Empty graph in labNetEdges"
    else filter (isNetworkLabEdge inGraph) $ labEdges inGraph

-- | netEdges takes a graph and returns list of network edges
netEdges :: Gr a b -> [Edge]
netEdges inGraph =
    if isEmpty inGraph then error "Empty graph in labNetEdges"
    else filter (isNetworkEdge inGraph) $ edges inGraph

-- | isTreeNode checks if node is network node
isTreeNode  :: Gr a b -> Node -> Bool
isTreeNode  inGraph inNode = (G.indeg inGraph inNode == 1) && (G.outdeg inGraph inNode > 0)

-- getRoots returns list of graph roots (labelled)
getRoots :: Gr a b -> [LNode a]
getRoots inGraph =
    if isEmpty inGraph then []
    else
        let nodeList =  labNodes inGraph
            -- rootBoolList = fmap (isRoot inGraph . fst) nodeList
            rootBoolList = fmap ((== 0).length) $ fmap (inn inGraph) $ fmap fst nodeList
            pairList = zip rootBoolList nodeList
            rootPairList =  filter ((==True).fst) pairList
            rootList = fmap snd rootPairList
        in
        rootList

-- | getIsolatedNodes returns list of labelled nodes with indegree=outdegree=0
-- should change to use G.deg == 0
getIsolatedNodes :: Gr a b -> [LNode a]
getIsolatedNodes inGraph =
    if isEmpty inGraph then []
    else
        let nodeList =  labNodes inGraph
            -- rootBoolList = fmap (isRoot inGraph . fst) nodeList
            in0BoolList = fmap ((== 0).length) $ fmap (inn inGraph) $ fmap fst nodeList
            out0BoolList = fmap ((== 0).length) $ fmap (out inGraph) $ fmap fst nodeList
            isolateBoolList = zipWith (&&) in0BoolList out0BoolList

            pairList = zip isolateBoolList nodeList
            isolatePairList =  filter ((==True).fst) pairList
            isolateList = fmap snd isolatePairList
        in
        isolateList

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

-- | flipLEdge flips orientation of labelled edge
flipLEdge :: LEdge b -> LEdge b
flipLEdge (u,v,w) = (v,u,w)


-- | isTree takes a graph and checks if there are anmy network nodes--if not returns True
isTree :: Gr a b -> Bool
isTree inGraph =
    if G.isEmpty inGraph then error "Empty graph in isTree"
    else
        let (_, _, _, netNodes) = splitVertexList inGraph
        in
        null netNodes

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

-- | pathToRoot takes a graph and a vertex and returns a pair of lists
-- of vertices and edges to root(s) in order of encountering them to root
-- if a tree--not necessarily if network--should work
pathToRoot :: (Eq a, Eq b) => Gr a b -> LNode a -> ([LNode a], [LEdge b])
pathToRoot inGraph inNode =
    if G.isEmpty inGraph then error "Empty graph in pathToRoot"
    else pathToRoot' inGraph [inNode] [] []

-- | pathToRoot' with accumulators
-- filter operators basically for networks so not retrace paths
-- includes roots as nodes
pathToRoot' :: (Eq a, Eq b) => Gr a b -> [LNode a] -> [LNode a] -> [LEdge b] -> ([LNode a], [LEdge b])
pathToRoot' inGraph inNodeList curNodeList curEdgeList =
    if null inNodeList then (reverse curNodeList, reverse curEdgeList)
    else
        let inNode = head inNodeList
        in
        -- root would already be inlist of nodes visited
        if isRoot inGraph (fst inNode) then pathToRoot' inGraph (tail inNodeList) curNodeList curEdgeList
        else
            let inLEdges = filter (`notElem` curEdgeList) $ inn inGraph (fst inNode)
                inNodes = filter (`notElem` (fmap fst curNodeList)) $ fmap fst3 inLEdges
                --inLabNodes = concatMap (labParents inGraph) (fmap fst3 inLEdges)
                inLabNodes = zip inNodes (fmap (fromJust . lab inGraph) inNodes)
            in
            pathToRoot' inGraph (inLabNodes ++ tail inNodeList) (inLabNodes ++ curNodeList) (inLEdges ++ curEdgeList)

-- | postOrderPathToNode takes a graph and two vertices nd returns a pair of lists
-- of vertices and edges to beteween them in order of encountering them from first to second
-- the path is post order to root so if second vertex is leaf-side of first node will hit root and fail
postOrderPathToNode :: (Eq a, Eq b) => Gr a b -> LNode a -> LNode a -> ([LNode a], [LEdge b])
postOrderPathToNode inGraph startNode endNode =
    if G.isEmpty inGraph then error "Empty graph in pathToRoot"
    else postOrderPathToNode' inGraph endNode [startNode] [] []

-- | postOrderPathToNode' with accumulators
-- filter operators basically for networks so not retrace paths
postOrderPathToNode' :: (Eq a, Eq b) => Gr a b -> LNode a -> [LNode a] -> [LNode a] -> [LEdge b] -> ([LNode a], [LEdge b])
postOrderPathToNode' inGraph endNode inNodeList curNodeList curEdgeList =
    if null inNodeList then (reverse curNodeList, reverse curEdgeList)
    else
        let inNode = head inNodeList
        in
        -- root would already be inlist of nodes visited
        if (fst inNode) == (fst endNode) then postOrderPathToNode' inGraph endNode (tail inNodeList) curNodeList curEdgeList
        else if isRoot inGraph (fst inNode) then error ("postOrderPathToNode hit root before end node.  Root index  " ++ (show $ fst inNode)
            ++  " edges " ++ (show $ fmap toEdge curEdgeList))
        else
            let inLEdges = filter (`notElem` curEdgeList) $ inn inGraph (fst inNode)
                inNodes = filter (`notElem` (fmap fst curNodeList)) $ fmap fst3 inLEdges
                inLabNodes = zip inNodes (fmap (fromJust . lab inGraph) inNodes)
            in
            postOrderPathToNode' inGraph endNode (inLabNodes ++ tail inNodeList) (inLabNodes ++ curNodeList) (inLEdges ++ curEdgeList)

-- | nodesAndEdgesBefore takes a graph and list of nodes to get list of nodes
-- and edges 'before' in the sense of leading to--ie between root and
-- (not including)) that node
-- call with ([], [])
-- filter operators basically for networks so not retrace paths
nodesAndEdgesBefore' :: (Eq a, Eq b, Show a) => Gr a b -> ([LNode a], [LEdge b]) -> [LNode a] -> ([LNode a], [LEdge b])
nodesAndEdgesBefore' inGraph curResults@(curNodes, curEdges) inNodeList
  | G.isEmpty inGraph = error "Input Graph is empty in nodesAndEdgesBefore"
  | null inNodeList = curResults
  | otherwise =
    let intoEdgeList = filter (`notElem` curEdges) $ inn inGraph (fst $ head inNodeList)
        intoNodeList = filter (`notElem` (fmap fst curNodes)) $ fmap fst3 intoEdgeList
        labelMaybeList = fmap (lab inGraph) intoNodeList
        labelList = fmap fromJust labelMaybeList
        intoLabNodeList = zip intoNodeList labelList
    in
    if Nothing `elem` labelMaybeList then error ("Empty node label in nodesAndEdgesBefore" ++ show intoLabNodeList)
    else nodesAndEdgesBefore' inGraph (intoLabNodeList ++ curNodes, intoEdgeList ++ curEdges) (intoLabNodeList ++ tail inNodeList)

-- | nodesAndEdgesBefore takes a graph and list of nodes to get list of nodes
-- and edges 'before' in the sense of leading to--ie between root and
-- (not including)) that node
-- filter operators basically for networks so not retrace paths
-- wrapper without accuulator
-- Does NOT Contain starting nodes
nodesAndEdgesBefore :: (Eq a, Eq b, Show a) => Gr a b -> [LNode a] -> ([LNode a], [LEdge b])
nodesAndEdgesBefore inGraph inNodeList = nodesAndEdgesBefore' inGraph ([],[]) inNodeList

-- | getEdgeListAfter wrapper around nodesAndEdgesAfter to return only edges
getEdgeListAfter :: (Eq a, Eq b,Show a) => (Gr a b, Node) -> [LEdge b]
getEdgeListAfter (inGraph, inNode) =
    snd $ nodesAndEdgesAfter inGraph [(inNode, fromJust $ lab inGraph inNode)]


-- | nodesAndEdgesAfter' takes a graph and list of nodes to get list of nodes
-- and edges 'after' in the sense of leading from-ie between (not including)) that node
-- and all the way to any leaves is connects to.
-- Does NOT Contain starting nodes
-- call with ([], [])
-- filter operators basically for networks so not retrace paths
nodesAndEdgesAfter' :: (Eq a, Eq b, Show a) => Gr a b -> ([LNode a], [LEdge b]) -> [LNode a] -> ([LNode a], [LEdge b])
nodesAndEdgesAfter' inGraph curResults@(curNodes, curEdges) inNodeList
  | G.isEmpty inGraph = error "Input Graph is empty in nodesAndEdgesAfter"
  | null inNodeList = curResults
  | otherwise =
    let fromEdgeList = filter (`notElem` curEdges) $ out inGraph (fst $ head inNodeList)
        fromNodeList = filter (`notElem` (fmap fst curNodes)) $ fmap snd3 fromEdgeList
        labelMaybeList = fmap (lab inGraph) fromNodeList
        
    in
    if Nothing `elem` labelMaybeList then error ("Empty node label in nodesAndEdgesAfter" ++ show (filter ((== Nothing) .snd) $ zip fromNodeList labelMaybeList) ++ " " ++ (show fromNodeList) ++ "\n" ++ (prettyIndices inGraph))
    else 
        let labelList = fmap fromJust labelMaybeList
            fromLabNodeList = zip fromNodeList labelList
        in
        nodesAndEdgesAfter' inGraph (fromLabNodeList ++ curNodes, fromEdgeList ++ curEdges) (fromLabNodeList ++ tail inNodeList)

-- | nodesAndEdgesAfter takes a graph and list of nodes to get list of nodes
-- and edges 'after' in the sense of leading from-ie between (not including)) that node
-- and all the way to any leaves is connects to.
-- wrapper wihtout accumulator
-- Does NOT Contain starting nodes
nodesAndEdgesAfter :: (Eq a, Eq b,Show a) => Gr a b -> [LNode a] -> ([LNode a], [LEdge b])
nodesAndEdgesAfter inGraph inNodeList = nodesAndEdgesAfter' inGraph ([],[]) inNodeList

-- | indexMatchNode returns True if two labelled nodes have same index
indexMatchNode :: LNode a -> LNode a -> Bool
indexMatchNode (a, _) (b, _) = if a == b then True else False


-- | coevalNodePairs generatres a list of pairs of nodes that must be potentially equal in
-- age (ie parents of networkNode)
coevalNodePairs :: (Eq a) => Gr a b -> [(LNode a, LNode a)]
coevalNodePairs inGraph =
    if G.isEmpty inGraph then []
    else
        let (_, _, _, netVertexList) = splitVertexList inGraph
            pairListList = fmap (labParents inGraph) $ fmap fst netVertexList
        in
        if null netVertexList then []
        else fmap makePairs pairListList
    where makePairs a = if length a /= 2 then error ("Not two parents for coevalNodePairs")
                        else (head a, a !! 1)
        

-- | indexMatchEdge returns True if two labelled edges have same node indices
indexMatchEdge :: LEdge b -> LEdge b -> Bool
indexMatchEdge (a,b,_) (c,d,_) = if a == c && b == d then True else False

-- | contractRootOut1Edge contracts indegree 0, outdegree 1, edges and removes the node in the middle
-- does one at a time and makes a graph and recurses
-- removes "tail" edges (single) from root to single child
contractRootOut1Edge :: (Show a, Show b) => Gr a b -> Gr a b
contractRootOut1Edge inGraph =
    if G.isEmpty inGraph then G.empty
    else
        let inOutDeg = getInOutDeg inGraph <$> labNodes inGraph
            out1RootList = filter ((==1) . thd3) $ filter ((==0) . snd3) inOutDeg
            out2RootList = filter ((>1) . thd3) $ filter ((==0) . snd3) inOutDeg
        in
        -- trace ("CRO1E :" ++ (show (inOutDeg, out1RootList, out2RootList))) (
        if null out1RootList then inGraph
        else if length out2RootList > 1 then error ("Multiple roots in graph in contractRootOut1Edge: " ++ (show $ length out2RootList))
        else if null out2RootList then
            let -- get root with out = 1 and its child, that childs' children, and edges
                rootVertex = head out1RootList
                childOfRoot = snd3 $ head $ out inGraph ((fst . fst3) rootVertex)
                childOfRootEdges = out inGraph childOfRoot

                newEdgeToAdd0 = ((fst . fst3) rootVertex, snd3 $ head childOfRootEdges, thd3 $ head childOfRootEdges)
                newEdgeToAdd1 = ((fst . fst3) rootVertex, snd3 $ last childOfRootEdges, thd3 $ last childOfRootEdges)

                -- create new Graph, deleting child node deltes three edges around it
                newGraph = insEdges [newEdgeToAdd0, newEdgeToAdd1] $ delNode childOfRoot inGraph
                

            in
            -- trace ("Removing tail edge root :" ++ (show $ snd3 $ head $ out inGraph ((fst . fst3) rootVertex)))
            contractRootOut1Edge $ reindexGraph newGraph
        else  -- case where mupltiple roots--combine by deleting in0out1 node and creting edge to its child from regular root.
            let rootVertex = head out1RootList
                in0out1RootIndex = (fst . fst3) rootVertex
                
                out2RootVertex = (fst . fst3 . head) out2RootList

                -- removes one of root out edges and inserts teh in0 out1 node adding two oew edges 
                root02EdgeToDelete = last $ out inGraph out2RootVertex
                newEdgeFrom02Root = (out2RootVertex, in0out1RootIndex, thd3 root02EdgeToDelete)
                newEdgeFrom01Root = (in0out1RootIndex, snd3 root02EdgeToDelete, thd3 root02EdgeToDelete)

                -- this relies on (a) root in first HTU
                nonOTUOut0Nodes = fmap (fst . fst3) $ filter ((>= out2RootVertex) . (fst . fst3)) $ filter ((==0) . thd3) $ filter ((>0) . snd3) inOutDeg

                newGraph = insEdges [newEdgeFrom02Root, newEdgeFrom01Root] $ delEdge (toEdge root02EdgeToDelete) $ delNodes nonOTUOut0Nodes inGraph

            in
            -- trace ("Removing extra edge root :" ++ (show $ (root02EdgeToDelete, newEdgeFrom02Root, newEdgeFrom01Root)))
            contractRootOut1Edge $ reindexGraph $ contractIn1Out1Edges $ newGraph

        -- )



-- | contractIn1Out1Edges contracts indegree 1, outdegree 1, edges and removes the node in the middle
-- does one at a time and makes a graph and recurses
contractIn1Out1Edges :: Gr a b -> Gr a b
contractIn1Out1Edges inGraph =
    if G.isEmpty inGraph then G.empty
    else
        let inOutDeg = getInOutDeg inGraph <$> labNodes inGraph
            degree11VertexList = filter ((==1) . thd3) $ filter ((==1) . snd3) inOutDeg
        in
        -- trace ("vertex 11:" ++ show degree11VertexList) (
        if null degree11VertexList then inGraph
        else
                let nodeToDelete = fst3 $ head degree11VertexList
                    inEdgeToDelete  = head $ inn inGraph $ fst nodeToDelete
                    outEdgeToDelete = head $ out inGraph $ fst nodeToDelete
                    newEdgeToAdd    = (fst3 inEdgeToDelete, snd3 outEdgeToDelete, thd3 inEdgeToDelete)
                    reindexedNodes = reindexNodes (fst nodeToDelete) [] $ labNodes inGraph
                    reindexedEdges = reindexEdges (fst nodeToDelete) [] (newEdgeToAdd : (labEdges inGraph))

                    newGraph = mkGraph reindexedNodes reindexedEdges
                    --newGraph = insEdge newEdgeToAdd $ delLNode nodeToDelete inGraph
                in
                -- trace ("Deleting Node " ++ show (fst nodeToDelete) ++ " " ++ show (inEdgeToDelete, outEdgeToDelete) ++ " inserting " ++ show  newEdgeToAdd)
                contractIn1Out1Edges newGraph
                -- )

-- | reindexNodes takes a node (assumes index and fst of node are the same) and a list of
--   nodes deleting the input node and reindexing all the other nodes with indices > than the input are reduced by 1
reindexNodes :: Int -> [LNode a] ->  [LNode a] -> [LNode a]
reindexNodes inNodeIndex curList nodeList =
    if null nodeList then reverse curList
    else
        let firstNode@(index, label) = head nodeList
        in
        if index < inNodeIndex then reindexNodes inNodeIndex (firstNode : curList) (tail nodeList)
        else if index == inNodeIndex then reindexNodes inNodeIndex curList (tail nodeList)
        else reindexNodes inNodeIndex ((index - 1, label) : curList) (tail nodeList)


-- | reindexEdges takes the index of a node that has/is being delted and reindexes indices
-- that are not incident on the node and deleted if incedent on node index
reindexEdges :: Int -> [LEdge b] -> [LEdge b] -> [LEdge b]
reindexEdges inNodeIndex curList edgeList =
    if null edgeList then curList
    else
        let (a,b,c) = head edgeList
            a' = if a < inNodeIndex then a
                 else a - 1
            b' = if b < inNodeIndex then b
                 else b - 1
        in
        -- incident on node to be deleted
        if a == inNodeIndex || b == inNodeIndex then reindexEdges inNodeIndex curList (tail edgeList)

        -- reindexed edge added in
        else reindexEdges inNodeIndex ((a', b', c) : curList) (tail edgeList)


-- | artPoint calls ap to get articulation points of graph
artPoint :: (Eq b) => Gr a b -> [Node]
artPoint inGraph = AP.ap $ B.undir inGraph


-- | makeNodeEdgePair takes a graph and extracts list of nodes and edges
-- returning a pair
makeNodeEdgePair :: Gr a b -> ([LNode a], [LEdge b])
makeNodeEdgePair inGraph =
    if isEmpty inGraph then ([],[])
    else (labNodes inGraph, labEdges inGraph)


-- | makeNodeEdgePairVect takes a graph and extracts list of nodes and edges
-- returning a pair of a vector of nodes and Vector of edges
makeNodeEdgePairVect :: Gr a b -> (V.Vector (LNode a), V.Vector (LEdge b))
makeNodeEdgePairVect inGraph =
    if isEmpty inGraph then (V.empty, V.empty)
    else (V.fromList $ labNodes inGraph, V.fromList $ labEdges inGraph)

-- | extractLeafGraph takes a Decorated graphs and cretes a graph from the leaves and no edges
extractLeafGraph :: Gr a b -> Gr a b
extractLeafGraph inGraph =
    if isEmpty inGraph then G.empty
    else
        let (_, leafNodes, _, _) = splitVertexList inGraph
        in
        mkGraph leafNodes []

-- | makes graph undirected
undir :: (Eq b) => Gr a b -> Gr a b
undir inGraph = B.undir inGraph

-- | finds bi connected components of a graph
bcc ::  Gr a b -> [Gr a b]
bcc inGraph = BCC.bcc inGraph

-- | removeNonLeafOut0NodesAfterRoot removed nodes (and edges attached) that are ourtdegree = zero
-- but have index > root
removeNonLeafOut0NodesAfterRoot :: (Eq a) => Gr a b -> Gr a b
removeNonLeafOut0NodesAfterRoot  inGraph =
    if isEmpty inGraph then empty
    else
        let (rootNodeList, putativeLeafNodeList, _, _) = splitVertexList inGraph
            rootIndex = (fst . head)  rootNodeList
            
            -- this should catch anything out = 0
            zeroOutNodeList  = filter ((> rootIndex) . fst) putativeLeafNodeList 
        in
        trace ("RNLOoN: " ++ (show (rootIndex, fmap fst zeroOutNodeList))) (
        if null zeroOutNodeList then inGraph
        else
            let newGraph = delNodes (fmap fst zeroOutNodeList) inGraph
            in
            removeNonLeafOut0NodesAfterRoot $ reindexGraph newGraph
        )
        
-- | removeNonLeafOut0Nodes removed nodes (and edges attached) that are ourtdegree = zero
-- but not in the leaf node list
-- does not reindex graph
removeNonLeafOut0Nodes :: (Eq a) => [LNode a] -> Gr a b -> Gr a b
removeNonLeafOut0Nodes leafList inGraph =
    if null leafList then error "Null leaf list in removeNonLeafOut0Nodes"
    else if isEmpty inGraph then empty
    else
        let nonLeafList = (labNodes inGraph) L.\\ leafList
            outdegreePairList = zip nonLeafList (fmap (outdeg inGraph) nonLeafList)
            (zeroOutNodeList, _) = unzip $ filter ((==0) .snd) outdegreePairList
        in
        if null zeroOutNodeList then inGraph
        else
            let newGraph = delNodes (fmap fst zeroOutNodeList) inGraph
            in
            removeNonLeafOut0Nodes leafList newGraph
        

-- | reindexGraph takes a graph and reindexes nodes and edges such that nodes
-- are sequential and the firt field matches their node index
reindexGraph :: Gr a b -> Gr a b
reindexGraph inGraph =
    if isEmpty inGraph then empty
    else
        let nodeList = labNodes inGraph
            newIndexList = [0..(length nodeList - 1)]
            nodeIndexPair = zip (fmap fst nodeList) newIndexList
            nodeIndexMap = MAP.fromList nodeIndexPair
            newNodeList = fmap (makeNewNode nodeIndexMap) nodeList
            newEdgeList = fmap (makeNewEdge nodeIndexMap) (labEdges inGraph)
        in
        mkGraph newNodeList newEdgeList

    where makeNewNode indexMap (a,b) = (fromJust $ MAP.lookup a indexMap, b)
          makeNewEdge indexMap (a,b,c) = (fromJust $ MAP.lookup a indexMap, fromJust $ MAP.lookup b indexMap, c)

-- | isBridge uses naive (component number) procedure to determine if edge is a bridge O(n)
isBridge :: Gr a b -> Edge -> Bool
isBridge inGraph inNode =
  if isEmpty inGraph then error ("Empty graph in isBridge")
  else
     let numComponents =  noComponents inGraph
         numComponents' = noComponents $ delEdge inNode inGraph
     in
     numComponents' > numComponents



-- FGL articulation point code--could be modified to get brisge edges in linear time}
------------------------------------------------------------------------------
-- Tree for storing the DFS numbers and back edges for each node in the graph.
-- Each node in this tree is of the form (v,n,b) where v is the vertex number,
-- n is its DFS number and b is the list of nodes (and their DFS numbers) that
-- lead to back back edges for that vertex v.
------------------------------------------------------------------------------
data DFSTree a = B (a,a,[(a,a)]) [DFSTree a]
     deriving (Eq, Show, Read)

------------------------------------------------------------------------------
-- Tree for storing the DFS and low numbers for each node in the graph.
-- Each node in this tree is of the form (v,n,l) where v is the vertex number,
-- n is its DFS number and l is its low number.
------------------------------------------------------------------------------
data LOWTree a = Brc (a,a,a) [LOWTree a]
     deriving (Eq, Show, Read)

------------------------------------------------------------------------------
-- Finds the back edges for a given node.
------------------------------------------------------------------------------
getBackEdges :: Node -> [[(Node,Int)]] -> [(Node,Int)]
getBackEdges _ [] = []
getBackEdges v ls = map head (filter (elem (v,0)) (tail ls))

------------------------------------------------------------------------------
-- Builds a DFS tree for a given graph. Each element (v,n,b) in the tree
-- contains: the node number v, the DFS number n, and a list of backedges b.
------------------------------------------------------------------------------
dfsTree :: Int -> Node -> [Node] -> [[(Node,Int)]] -> Gr a b -> ([DFSTree Int],Gr a b,Int)
dfsTree n _ []      _ g             = ([],g,n)
dfsTree n _ _       _ g | isEmpty g = ([],g,n)
dfsTree n u (v:vs) ls g = case G.match v g of
                            (Nothing, g1) -> dfsTree n u vs ls g1
                            (Just c , g1) -> (B (v,n+1,bck) ts:ts', g3, k)
                             where  bck        = getBackEdges v ls
                                    (ts, g2,m) = dfsTree (n+1) v sc ls' g1
                                    (ts',g3,k) = dfsTree m v vs ls g2
                                    ls'        = ((v,n+1):sc'):ls
                                    sc'        = map (\x->(x,0)) sc
                                    sc         = G.suc' c

------------------------------------------------------------------------------
-- Finds the minimum between a dfs number and a list of back edges' dfs
-- numbers.
------------------------------------------------------------------------------
minbckEdge :: Int -> [(Node,Int)] -> Int
minbckEdge n [] = n
minbckEdge n bs = min n (minimum (map snd bs))

------------------------------------------------------------------------------
-- Returns the low number for a node in a subtree.
------------------------------------------------------------------------------
getLow :: LOWTree Int -> Int
getLow (Brc (_,_,l) _) = l

------------------------------------------------------------------------------
-- Builds a low tree from a DFS tree. Each element (v,n,low) in the tree
-- contains: the node number v, the DFS number n, and the low number low.
------------------------------------------------------------------------------
lowTree :: DFSTree Int -> LOWTree Int
lowTree (B (v,n,[]  ) [] ) = Brc (v,n,n) []
lowTree (B (v,n,bcks) [] ) = Brc (v,n,minbckEdge n bcks) []
lowTree (B (v,n,bcks) trs) = Brc (v,n,lowv) ts
                             where lowv     = min (minbckEdge n bcks) lowChild
                                   lowChild = minimum (map getLow ts)
                                   ts       = map lowTree trs

------------------------------------------------------------------------------
-- Builds a low tree for a given graph. Each element (v,n,low) in the tree
-- contains: the node number v, the DFS number n, and the low number low.
------------------------------------------------------------------------------
getLowTree :: Gr a b -> Node -> LOWTree Int
getLowTree g v = lowTree (head dfsf)
                  where (dfsf, _, _) = dfsTree 0 0 [v] [] g

------------------------------------------------------------------------------
-- Tests if a node in a subtree is an articulation point. An non-root node v
-- is an articulation point iff there exists at least one child w of v such
-- that lowNumber(w) >= dfsNumber(v). The root node is an articulation point
-- iff it has two or more children.
------------------------------------------------------------------------------
isap :: LOWTree Int -> Bool
isap (Brc (_,_,_) []) = False
isap (Brc (_,1,_) ts) = length ts > 1
isap (Brc (_,n,_) ts) = not (null ch)
                        -- modify for bridges
                        -- where ch = filter ( >= n) (map getLow ts)
                        where ch = filter ( >= n) (map getLow ts)
------------------------------------------------------------------------------
-- Finds the articulation points by traversing the low tree.
------------------------------------------------------------------------------
arp :: LOWTree Int -> [Node]
arp (Brc (v,1,_) ts) | length ts > 1         = v:concatMap arp ts
                     | otherwise             =   concatMap arp ts
arp (Brc (v,n,l) ts) | isap (Brc (v,n,l) ts) = v:concatMap arp ts
                     | otherwise             =   concatMap arp ts

------------------------------------------------------------------------------
-- Finds the articulation points of a graph starting at a given node.
------------------------------------------------------------------------------
artpoints :: Gr a b -> Node -> [Node]
artpoints g v = arp (getLowTree g v)

{-|
   Finds the articulation points for a connected undirected graph,
   by using the low numbers criteria:

   a) The root node is an articulation point iff it has two or more children.

   b) An non-root node v is an articulation point iff there exists at least
      one child w of v such that lowNumber(w) >= dfsNumber(v).
-}
ap' :: Gr a b -> [Node]
ap' g = artpoints g v where ((_,v,_,_),_) = G.matchAny g

-- | cyclic maps to cyclic function in module Cyclic.hs
cyclic :: Gr a b -> Bool
cyclic inGraph =
    -- trace ("Cyclic:" ++ (show $ C.cyclic inGraph))
    C.cyclic inGraph

-- | testEdge nodeList fullEdgeList) counter
-- chnage to input graph and delete edge from graph as opposed to making new graphs each time.
-- should be much faster using P.delLEdge (since only one edge to delete)
testEdge :: (Eq b) => P.Gr a b -> G.LEdge b -> [G.LEdge b]
testEdge fullGraph candidateEdge@(e,u,_) =
  let newGraph = G.delLEdge candidateEdge fullGraph
      bfsNodes = BFS.bfs e newGraph
      foundU = L.find (== u) bfsNodes
  in
  [candidateEdge | isNothing foundU]

-- | transitiveReduceGraph take list of nodes and edges, deletes each edge (e,u) in turn makes graph,
-- checks for path between nodes e and u, if there is delete edge otherwise keep edge in list for new graph
-- transitive reduction  Aho et al. 1972
-- this not iterative with new graphs--shold it be?
transitiveReduceGraph ::  (Eq b) => Gr a b -> Gr a b
transitiveReduceGraph fullGraph =
  let requiredEdges = fmap (testEdge fullGraph) (labEdges fullGraph)
      newGraph = G.mkGraph (labNodes fullGraph) (concat requiredEdges)
  in
  newGraph

-- | getCoevalConstraintEdges takes a graph and a network node and creates two lists: one of edges
-- "before" (ie towards root) and a second "after (ie away from root)
-- this defines a coeval constraint.  No network edge can be added that would be directed
-- from the before group to the after
getCoevalConstraintEdges :: (Eq a, Eq b, Show a) => Gr a b -> LNode a -> ([LEdge b],[LEdge b])
getCoevalConstraintEdges inGraph inNode =
   if isEmpty inGraph then error "Empty input graph in getCoevalConstraintEdges"
   else
       let (_, edgeBeforeList) = nodesAndEdgesBefore inGraph [inNode]
           (_, edgeAfterList) = nodesAndEdgesAfter inGraph [inNode]
       in
       (edgeBeforeList, edgeAfterList)


-- | getGraphCoevalConstraints takes a greaph and returns coeval constraints based on network nodes
getGraphCoevalConstraints :: (Eq a, Eq b, Show a, NFData b) => Gr a b -> [([LEdge b],[LEdge b])]
getGraphCoevalConstraints inGraph =
   if isEmpty inGraph then error "Empty input graph in getGraphCoevalConstraints"
   else
       let (_, _, _, networkNodeList) = splitVertexList inGraph
       in
       if null networkNodeList then []
       else PU.seqParMap rdeepseq (getCoevalConstraintEdges inGraph) networkNodeList -- `using`  PU.myParListChunkRDS

-- | getGraphCoevalConstraintsNodes takes a graph and returns coeval constraints based on network nodes
-- and nodes as a triple
getGraphCoevalConstraintsNodes :: (Eq a, Eq b, Show a, NFData b) => Gr a b -> [(LNode a, [LEdge b],[LEdge b])]
getGraphCoevalConstraintsNodes inGraph =
   if isEmpty inGraph then error "Empty input graph in getGraphCoevalConstraints"
   else
       let (_, _, _, networkNodeList) = splitVertexList inGraph
       in
       if null networkNodeList then []
       else
            let (edgeBeforeList, edgeAfterList) = unzip (PU.seqParMap rdeepseq  (getCoevalConstraintEdges inGraph) networkNodeList) --  `using`  PU.myParListChunkRDS)
            in zip3 networkNodeList edgeBeforeList edgeAfterList

-- | meetsAllCoevalConstraintsNodes checks constraint pair list and examines
-- whether one edge is from before and one after--if so fails False
-- else True if all pass
-- new edghe that woudl be creatated in edge1 -> edge2
-- checks if starting and ending vertices of edges to be linked are on same side of each 
-- coeval contraint
meetsAllCoevalConstraintsNodes :: (Eq b) => [(Node, Node, [Node], [Node], [Node], [Node])] -> LEdge b -> LEdge b -> Bool
meetsAllCoevalConstraintsNodes constraintList edge1@(u,v,_) edge2@(u',v',_) =
   if null constraintList then True
   else
       let (_, _, aNodesBefore, bNodesBefore, aNodesAfter, bNodesAfter) = head constraintList
       in
       --checks if edge nodes would stradle existing coeveal nodes
       if (u `elem` aNodesAfter) && (u' `elem` bNodesBefore) then False
       else if (u' `elem` aNodesAfter) && (u `elem` bNodesBefore) then False
       else if (v `elem` aNodesAfter) && (v' `elem` bNodesBefore) then False
       else if (v' `elem` aNodesAfter) && (v `elem` bNodesBefore) then False

       else if (u `elem` aNodesBefore) && (u' `elem` bNodesAfter) then False
       else if (u' `elem` aNodesBefore) && (u `elem` bNodesAfter) then False
       else if (v `elem` aNodesBefore) && (v' `elem` bNodesAfter) then False
       else if (v' `elem` aNodesBefore) && (v `elem` bNodesAfter) then False


       else meetsAllCoevalConstraintsNodes (tail constraintList) edge1 edge2


-- | meetsAllCoevalConstraintsEdges checks constraint pair list and examines
-- whether one edge is from before and one after--if so fails False
-- else True if all pass
    -- not correct I don't htink
meetsAllCoevalConstraintsEdges :: (Eq b) =>[([LEdge b],[LEdge b])] -> LEdge b -> LEdge b -> Bool
meetsAllCoevalConstraintsEdges constraintList edge1 edge2 =
   if null constraintList then True
   else
       let (beforeList, afterList) = head constraintList
       in
       if edge1 `elem` beforeList && edge2 `elem` afterList then False
       else if edge2 `elem` beforeList && edge1 `elem` afterList then False
       else meetsAllCoevalConstraintsEdges (tail constraintList) edge1 edge2

-- | insertDeleteEdges takes a  graphs and list of nodes and edges to add and delete and creates new graph
insertDeleteEdges :: (Show a, Show b) => Gr a b -> ([LEdge b], [Edge]) ->  Gr a b
insertDeleteEdges inGraph (edgesToAdd, edgesToDelete) =
   let editedGraph = insEdges edgesToAdd $ delEdges edgesToDelete inGraph
   in
   -- trace ("AGE: " ++ (show editStuff) ++ "\nIn graph:\n" ++ (prettify inGraph) ++ "New Graph:\n" ++ (prettify editedGraph))
   editedGraph


-- | notMatchEdgeIndices reutnr True if not in edge list but doesn't compare label only indices
notMatchEdgeIndices :: [Edge] -> LEdge b -> Bool
notMatchEdgeIndices unlabeledEdegList labelledEdge =
    if toEdge labelledEdge `elem` unlabeledEdegList then False
    else True

-- | isGraphTimeConsistent retuns False if graph fails time consistency
isGraphTimeConsistent ::  (Show a,Eq a,Eq b, NFData a) => Gr a b -> Bool
isGraphTimeConsistent inGraph =
  if isEmpty inGraph then True
  else 
    if isTree inGraph then True
    else 
      let coevalNodeConstraintList = coevalNodePairs inGraph
          coevalNodeConstraintList' = PU.seqParMap rdeepseq (addBeforeAfterToPair inGraph) coevalNodeConstraintList -- `using`  PU.myParListChunkRDS
          coevalPairsToCompareList = getListPairs coevalNodeConstraintList'
          timeOffendingEdgeList = getEdgesToRemoveForTime inGraph coevalPairsToCompareList
      in
      null timeOffendingEdgeList

-- | addBeforeAfterToPair adds before and after node list to pari of nodes for later use
-- in time contraint edge removal
addBeforeAfterToPair :: (Show a,Eq a,Eq b) 
                            => Gr a b 
                            -> (LNode a, LNode a) 
                            -> (LNode a, LNode a, [LNode a], [LNode a], [LNode a], [LNode a])
addBeforeAfterToPair inGraph (a,b) = 
  if isEmpty inGraph then error "Empty graph in addBeforeAfterToPair"
  else
      let  (aNodesBefore, _) = nodesAndEdgesBefore inGraph [a]
           (bNodesBefore, _) = nodesAndEdgesBefore inGraph [b]
           (aNodesAfter, _) = nodesAndEdgesAfter inGraph [a]
           (bNodesAfter, _) = nodesAndEdgesAfter inGraph [b]
      in
      (a,b, aNodesBefore, bNodesBefore, aNodesAfter, bNodesAfter)

-- | getEdgesToRemoveForTime recursive looks at each pair of nodes that should
-- be potentially coeval based on network vertices. And see if it conflicts with
-- other pairs that should also be coeval.  
-- two pairs of nodes (a,b) and (a',b'), one for each net node tio be compared
-- if a "before" a'  then b must be before b'
-- if a "after" a' then b must be after b'
-- otherwise its a time violation
-- returns net edge to delte in second pair
getEdgesToRemoveForTime :: (Show a,Eq a,Eq b) 
                        => Gr a b 
                        -> [((LNode a, LNode a, [LNode a], [LNode a], [LNode a], [LNode a]),(LNode a, LNode a, [LNode a], [LNode a], [LNode a], [LNode a]))] 
                        -> [Edge]
getEdgesToRemoveForTime inGraph inNodePairList =
  if isEmpty inGraph then []
  else if null inNodePairList then []
  else
    let ((_,_, aNodesBefore, bNodesBefore, aNodesAfter, bNodesAfter),(a',b', _, _,_,_)) = head inNodePairList
    in
    -- trace ("GETRT: " ++ (show inNodePairList)) (
    -- condition if a before a' then b before b' to be ok
    if (a' `elem` aNodesAfter) && (b' `elem` bNodesBefore) then 
      let edgeToRemove = (head $ filter (isNetworkEdge inGraph) $ fmap toEdge $ out inGraph $ fst b') 
      in
      -- trace ("\tRemoving network edge due to time consistancy: " ++ (show edgeToRemove))
      -- trace ("GERT Edges0:" ++ (show edgeToRemove) ++ " " ++ (show $ fmap toEdge $ out inGraph $ fst b') ++ " Net: " ++ (show $ fmap (isNetworkEdge inGraph) $ fmap toEdge $ out inGraph $ fst b'))
      edgeToRemove : getEdgesToRemoveForTime inGraph (tail inNodePairList)
    else if (a' `elem` aNodesBefore) && (b' `elem` bNodesAfter) then 
      let edgeToRemove = (head $ filter (isNetworkEdge inGraph) $ fmap toEdge $ out inGraph $ fst b') 
      in 
      -- trace ("GERT Edges1:" ++ (show edgeToRemove) ++ " " ++ (show $ fmap toEdge $ out inGraph $ fst b'))
      -- trace ("\tRemoving network edge due to time consistancy : "++ (show edgeToRemove))
      edgeToRemove : getEdgesToRemoveForTime inGraph (tail inNodePairList)
    else if (b' `elem` aNodesAfter) && (a' `elem` bNodesBefore) then 
      let edgeToRemove = (head $ filter (isNetworkEdge inGraph) $ fmap toEdge $ out inGraph $ fst b') 
      in
      -- trace ("\tRemoving network edge due to time consistancy: " ++ (show edgeToRemove))
      -- trace ("GERT Edges0:" ++ (show edgeToRemove) ++ " " ++ (show $ fmap toEdge $ out inGraph $ fst b') ++ " Net: " ++ (show $ fmap (isNetworkEdge inGraph) $ fmap toEdge $ out inGraph $ fst b'))
      edgeToRemove : getEdgesToRemoveForTime inGraph (tail inNodePairList)
    else if (b' `elem` aNodesBefore) && (a' `elem` bNodesAfter) then 
      let edgeToRemove = (head $ filter (isNetworkEdge inGraph) $ fmap toEdge $ out inGraph $ fst b') 
      in 
      -- trace ("GERT Edges1:" ++ (show edgeToRemove) ++ " " ++ (show $ fmap toEdge $ out inGraph $ fst b'))
      -- trace ("\tRemoving network edge due to time consistancy : "++ (show edgeToRemove))
      edgeToRemove : getEdgesToRemoveForTime inGraph (tail inNodePairList)
    else getEdgesToRemoveForTime inGraph (tail inNodePairList)
    -- )
                                      
-- | getEdgeSplitList takes a graph and returns list of edges
-- that split a graph increasing the number of components by 1
-- this is quadratic
-- should change to Tarjan's algorithm (linear)
-- everything else in there is O(n^2-3) so maybe doesn't matter
-- filters out edges with parent nodes that are out degree 1 and root edges
getEdgeSplitList :: (Show a, Show b, Eq b) => Gr a b -> [LEdge b]
getEdgeSplitList inGraph =
  if isEmpty inGraph then error ("Empty graph in getEdgeSplitList")
  else
      let origNumComponents = noComponents inGraph
          origEdgeList = labEdges inGraph
          edgeDeleteComponentNumberList = fmap noComponents $ fmap (flip delEdge inGraph) (fmap toEdge origEdgeList)
          bridgeList = fmap snd $ filter ((> origNumComponents) . fst) $ zip edgeDeleteComponentNumberList origEdgeList

          -- filter out edges starting in an outdegree 1 node (network or in out 1) node
          -- this would promote an HTU to a leaf later.  Its a bridge, but not what need
          -- bridgeList' = filter ((not . isRoot inGraph) .fst3 ) $ filter ((not. isNetworkNode inGraph) . snd3)  $ filter  ((not . isOutDeg1Node inGraph) . fst3) bridgeList
          bridgeList' = filter ((not . isRoot inGraph) .fst3 ) $ filter ((not . isNetworkEdge inGraph) . toEdge) $ filter  ((not . isOutDeg1Node inGraph) . fst3) bridgeList
      in

       -- trace ("BridgeList" ++ (show $ fmap toEdge bridgeList') ++ "\nGraph\n" ++ (prettyIndices inGraph))
       bridgeList'

-- | splitGraphOnEdge takes a graph and an edge and returns a single graph but with two components
-- the roots of each component are returned with two graphs, with broken edge contracted, and 'naked'
-- node returned (this is teh original connection vertex on base graph connected to pruned graph).  
-- The naked node is used for rejoining the two components during rearrangement
-- (SplitGraph, root of component that has original root, root of component that was cut off, naked node left over)
-- this function does not check whether edge is a 'bridge'
splitGraphOnEdge :: (Show b) => Gr a b -> LEdge b -> (Gr a b, Node, Node, Node)
splitGraphOnEdge inGraph (e,v,l) =
  if isEmpty inGraph then error "Empty graph in splitGraphOnEdge"
  else if (length $ getRoots inGraph) /= 1 then error ("Incorrect number roots in splitGraphOnEdge--must be 1: " ++ (show $ fmap fst $ getRoots inGraph))
  else
      let childrenENode = (descendants inGraph e) L.\\ [v]
          parentsENode = parents inGraph e
          newEdge = (head parentsENode, head childrenENode, l)
          edgesToDelete = [(head parentsENode, e), (e, head childrenENode)] -- (e,v)

          -- make new graph
          splitGraph = insEdge newEdge $ delEdges edgesToDelete inGraph
      in
      if length childrenENode /= 1 then error ("Incorrect number of children of edge to split--must be 1: " ++ (show ((e,v), childrenENode)))
      else if length parentsENode /= 1 then error ("Incorrect number of parents of edge to split--must be 1: " ++ (show ((e,v), parentsENode)))
      else
          -- trace ("SGE:" ++ (show (childrenENode, parentsENode, newEdge, edgesToDelete)))
          (splitGraph, fst $ head $ getRoots inGraph, v, e)

-- | splitGraphOnEdge' like splitGraphOnEdge above but returns edges created and destroyed as well
-- used in Goodman-Bermer and could make swap more efficient as well.
splitGraphOnEdge' :: (Show b) => Gr a b -> LEdge b -> (Gr a b, Node, Node, Node, LEdge b, [Edge])
splitGraphOnEdge' inGraph (e,v,l) =
  if isEmpty inGraph then error "Empty graph in splitGraphOnEdge"
  else if (length $ getRoots inGraph) /= 1 then error ("Incorrect number roots in splitGraphOnEdge--must be 1: " ++ (show $ fmap fst $ getRoots inGraph))
  else
      let childrenENode = (descendants inGraph e) L.\\ [v]
          parentsENode = parents inGraph e
          newEdge = (head parentsENode, head childrenENode, l)
          edgesToDelete = [(head parentsENode, e), (e, head childrenENode)] -- (e,v)

          -- make new graph
          splitGraph = insEdge newEdge $ delEdges edgesToDelete inGraph
      in
      if length childrenENode /= 1 then error ("Incorrect number of children of edge to split--must be 1: " ++ (show ((e,v), childrenENode)))
      else if length parentsENode /= 1 then error ("Incorrect number of parents of edge to split--must be 1: " ++ (show ((e,v), parentsENode)))
      else
          -- trace ("SGE:" ++ (show (childrenENode, parentsENode, newEdge, edgesToDelete)))
          (splitGraph, fst $ head $ getRoots inGraph, v, e, newEdge, edgesToDelete)

-- | joinGraphOnEdge takes a graph and adds an edge reducing the component number
-- expected ot be two components to one in SPR/TBR
-- assumes that first node of edge (e,v,l) is 'naked' ie avaiable to make edges but is in graph
-- created from splitGraphOnEdge
joinGraphOnEdge :: (Show a,Show b) => Gr a b -> LEdge b -> Node ->Gr a b
joinGraphOnEdge inGraph (x,y,l) parentofPrunedSubGraph =
  if isEmpty inGraph then error ("Empty graph in joinGraphOnEdge")
  else
      let edgeToCreate0 = (x, parentofPrunedSubGraph, l)
          edgeToCreate1 = (parentofPrunedSubGraph, y, l)
          -- edgeToCreate2 = (parentofPrunedSubGraph, graphToJoinRoot, l)
      in
      -- make new graph
      -- trace ("JGE:" ++ (show edgeToInvade) ++ " " ++ (show (parentofPrunedSubGraph, graphToJoinRoot))) --  ++ "\n" ++ (prettify inGraph))
      insEdges [edgeToCreate0, edgeToCreate1] $ delEdge (x,y) inGraph

-- | parentsInChain checks for parents in chain ie network edges
-- that implies a network event between nodes where one is the ancestor of the other
-- a time violation
parentInChain :: (Show a, Eq a, Eq b) => Gr a b -> Bool
parentInChain inGraph =
  if isEmpty inGraph then error "Null graph in parentInChain"
  else
    let   (_, _, _, netVertexList) = splitVertexList inGraph
          parentNetVertList = fmap (labParents inGraph) $ fmap fst netVertexList

          -- get list of nodes that are transitively equal in age
          concurrentList = mergeConcurrentNodeLists parentNetVertList []
          concurrentPairList = concatMap getListPairs concurrentList

          -- get pairs that violate concurrency
          violatingConcurrentPairs = concatMap (concurrentViolatePair inGraph) concurrentPairList
    in
    if null violatingConcurrentPairs then False
    else True

-- | getSisterSisterEdgeList take a graph and returns list of edges where two network nodes
-- have the same two parents
getSisterSisterEdgeList :: Gr a b -> [Edge]
getSisterSisterEdgeList inGraph =
  if isEmpty inGraph then []
  else
      let (_, _, _, netVertexList) = splitVertexList inGraph
      in
      if null netVertexList then []
      else getSisterSisterEdgeByNetVertex inGraph netVertexList 

-- | getSisterSisterEdgeByNetVertex takes a list of vertices and recursively generate a list of edges to delete
getSisterSisterEdgeByNetVertex ::  Gr a b -> [LNode a] -> [Edge]
getSisterSisterEdgeByNetVertex inGraph netNodeList =
  if null netNodeList then []
  else 
    let firstNode = fst $ head netNodeList
        parentsList = parents inGraph firstNode
        grandParentsList = fmap (parents inGraph) parentsList
        sameGrandParentList = L.foldl1' L.intersect grandParentsList
    in
    if null sameGrandParentList then getSisterSisterEdgeByNetVertex inGraph (tail netNodeList)
    else
      -- trace ("Found sister-sister")
      (toEdge $ head $ inn inGraph firstNode) : getSisterSisterEdgeByNetVertex inGraph (tail netNodeList)

-- | concurrentViolatePair takes a pair of nodes and sees if either is ancestral to the other--if so returns pair
-- as list otherwise null list
concurrentViolatePair :: (Eq a, Show a, Eq b) => Gr a b  -> (LNode a, LNode a) -> [(LNode a, LNode a)]
concurrentViolatePair inGraph (node1, node2) =
  if isEmpty inGraph then error "Empty graph in concurrentViolatePair"
  else
    let (nodesBeforeFirst, _)  = nodesAndEdgesBefore inGraph [node1]
        (nodesBeforeSecond, _) = nodesAndEdgesBefore inGraph [node2]
    in
    if node2 `elem` nodesBeforeFirst then [(node1, node2)]
    else if node1 `elem` nodesBeforeSecond then [(node1, node2)]
    else []


-- | mergeConcurrentNodeLists takes a list os lists  and returns a list os lists of merged lists
-- lists are merged if they share any elements
mergeConcurrentNodeLists :: (Eq a) => [[LNode a]] -> [[LNode a]] -> [[LNode a]]
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

-- | sortEdgeListByDistance sorts edges by distance (in edges) from edge pair of vertices
-- cretes a list of edges into (but traveling away from) an initial eNOde and away from
-- an initial vNode adding new nodes to those lists as encountered by traversing edges.
-- the eidea is theat the nodes from a directed edge (eNode, vNode)
-- the list is creted at each round from the "in" and "out" edge lists
-- so they are in order of 1 edge 2 edges etc.
sortEdgeListByDistance :: Gr a b -> [Node] -> [Node] -> [LEdge b]
sortEdgeListByDistance inGraph eNodeList vNodeList =
  if isEmpty inGraph then error ("Empty graph in edgeListByDistance")
  else if (null eNodeList && null vNodeList) then []
  else
        -- get edges 'in' to eNodeList
    let inEdgeList = concatMap (inn inGraph) eNodeList
        newENodeList = fmap fst3 inEdgeList

        -- get edges 'out' from vNodeList
        outEdgeList = concatMap (out inGraph) vNodeList
        newVNodeList = fmap snd3 outEdgeList
    in
    inEdgeList ++ outEdgeList ++ (sortEdgeListByDistance inGraph newENodeList newVNodeList)

-- | switchRootTree takes a new root vertex index of a tree and switches the existing root (and all relevent edges)
-- to new index
switchRootTree :: (Show a) => Int -> Gr a b -> Gr a b
switchRootTree newRootIndex inGraph =
  if isEmpty inGraph then empty
  else
    let rootList = getRoots inGraph
        (newRootCurInEdges, newRootCurOutEdges) = getInOutEdges inGraph newRootIndex
        oldRootEdges = out inGraph $ fst $ head rootList

    in
    -- not a directed tree
    if length rootList /= 1 then error ("Graph input to switchRootTree is not a tree--not single root:" ++ (show rootList))

    -- same root
    else if newRootIndex == (fst $ head rootList) then inGraph
    else
        -- create new edges and delete the old ones
        let newEdgesToAdd = fmap (flipVertices (fst $ head rootList) newRootIndex) (newRootCurInEdges ++ newRootCurOutEdges ++ oldRootEdges)
        in
        insEdges newEdgesToAdd $ delLEdges (newRootCurInEdges ++ newRootCurOutEdges ++ oldRootEdges) inGraph

-- | flipVertices takes an old vertex index and a new vertex index and inserts one for the other
-- in a labelled edge
flipVertices :: Node -> Node ->  LEdge b -> LEdge b
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

-- | rerootTree takes a graph and reroots based on a vertex index (usually leaf outgroup)
--   if input is a forest then only roots the component that contains the vertex wil be rerooted
--   unclear how will effect network edges--will need to verify that does not create cycles
--   multi-rooted components (as opposed to forests) are unaffected with trace warning thrown
--   after checking for existing root and multiroots, should be O(n) where 'n is the length
--   of the path between the old and new root
rerootTree :: (Show a, Show b, Eq b) => Node -> Gr a b -> Gr a b
rerootTree rerootIndex inGraph =
  --trace ("In reroot Graph: " ++ show rerootIndex) (
  if isEmpty inGraph then inGraph
  else
    let componentList = components inGraph
        parentNewRootList = pre inGraph rerootIndex
        newRootOrigEdge = head $ inn inGraph rerootIndex
        parentRootList = fmap (isRoot inGraph) parentNewRootList
        outgroupInComponent = fmap (rerootIndex `elem`) componentList
        componentWithOutgroup = filter ((== True).fst) $ zip outgroupInComponent componentList
        -- (_, inNewRoot, outNewRoot) = getInOutDeg inGraph (labelNode inGraph rerootIndex)
    in

    -- rerooting on root so no indegree edges
    -- this for wagner build reroots where can try to reroot on leaf not yet added
    if null $ inn inGraph rerootIndex then inGraph -- error ("Rerooting on indegree 0 node " ++ (show rerootIndex) ++ "\n" ++ prettyIndices inGraph) -- empty


    else if null componentWithOutgroup then inGraph  -- error ("Error rooting wierdness in rerootTree " ++ (show rerootIndex) ++ "\n" ++ prettyIndices inGraph) -- empty

    -- check if new outtaxon has a parent--shouldn't happen-but could if its an internal node reroot
    else if null parentNewRootList || (True `elem` parentRootList) then inGraph
                                                              else (if null componentWithOutgroup then error ("Outgroup index " ++ show rerootIndex ++ " not found in graph")
    else
        --trace ("RRT: " ++ (show (rerootIndex, inNewRoot, outNewRoot))) ( 
        -- reroot component with new outtaxon
        let componentWithNewOutgroup = snd $ head componentWithOutgroup
            (_, originalRootList) =  unzip $ filter ((==True).fst) $ zip (fmap (isRoot inGraph) componentWithNewOutgroup) componentWithNewOutgroup
            numRoots = length originalRootList
            orginalRoot = head originalRootList
            originalRootEdges = out inGraph orginalRoot

        in

        if numRoots == 0 then error ("No root in rerootTree: Attempting to reroot on edge to node " ++ (show rerootIndex) ++ "\n" ++ prettyIndices inGraph) --empty

        -- check if outgroup in a multirooted component
        -- if wagner build this is ok
        else if numRoots > 1 then inGraph -- error ("Error: Attempting to reroot multi-rooted component") -- inGraph
        else
          --reroot graph safely automatically will only affect the component with the outgroup
          -- delete old root edge and create two new edges from oringal root node.
          -- keep orignl root node and delte/crete new edges when they are encounterd
          --trace ("Moving root from " ++ (show orginalRoot) ++ " to " ++  (show rerootIndex)) (
          let leftChildEdge = (orginalRoot, rerootIndex, edgeLabel $ head originalRootEdges)
              rightChildEdge = (orginalRoot, fst3 newRootOrigEdge, edgeLabel $ last originalRootEdges)

              --  this assumes 2 children of old root -- shouled be correct as Phylogenetic Graph
              newEdgeOnOldRoot = if (length originalRootEdges) /= 2 then error ("Number of root out edges /= 2 in rerootGraph: " ++ (show $ length originalRootEdges)
                ++ " root index: " ++ (show (orginalRoot, rerootIndex)) ++ "\nGraph:\n" ++ (prettyIndices inGraph))
                                 else (snd3 $ head originalRootEdges, snd3 $ last originalRootEdges, thd3 $ head originalRootEdges)

              newRootEdges = [leftChildEdge, rightChildEdge, newEdgeOnOldRoot]
              newGraph = insEdges newRootEdges $ delLEdges (newRootOrigEdge : originalRootEdges) inGraph

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
rerootDisplayTree :: (Show a, Show b, Eq a, Eq b) => Node -> Node -> Gr a b -> Gr a b
rerootDisplayTree orginalRootIndex rerootIndex inGraph =
  --trace ("In reroot Graph: " ++ show rerootIndex) (
  if isEmpty inGraph then inGraph
  else
    let -- componentList = components inGraph'
        
        -- hack---remove when figured out
        {-
        inGraph = if inGraphType == SoftWired then inGraph' -- removeDuplicateEdges inGraph'
                  else inGraph'
        -}

        parentNewRootList = pre inGraph rerootIndex
        newRootOrigEdge = head $ inn inGraph rerootIndex
        parentRootList = fmap (isRoot inGraph) parentNewRootList
        -- outgroupInComponent = fmap (rerootIndex `elem`) componentList
        -- componentWithOutgroup = filter ((== True).fst) $ zip outgroupInComponent componentList
        (_, inNewRoot, outNewRoot) = getInOutDeg inGraph (labelNode inGraph rerootIndex)

        {- Checks for valid trees
        cyclicString = if cyclic inGraph then " Is cyclic "
                          else  " Not cyclic "

        parentInCharinString = if parentInChain inGraph then " Is Parent in Chain "
                                  else " Not Parent in Chain "

        duplicatedsEdgeString = if not (hasDuplicateEdge inGraph)  then " No duplicate edges "
                                else 
                                  let dupEdgeList' = getDuplicateEdges inGraph
                                  in (" Has duplicate edges: " ++ (show dupEdgeList')) 
        -}
    in

    -- trace ("RRT In: " ++ cyclicString ++ parentInCharinString ++ duplicatedsEdgeString) (

    -- check if cycle and exit if so
    -- don't reroot on in=out=1 since same as it descendent edge 
    if (inNewRoot == 1) && (outNewRoot == 1) then 
      -- trace ("RRT: in 1 out 1") 
      inGraph

    -- else if cyclic inGraph then empty -- inGraph 

    
    -- rerooting on root so no indegree edges
    -- this for wagner build reroots where can try to reroot on leaf not yet added
    else if null $ inn inGraph rerootIndex then inGraph -- error ("Rerooting on indegree 0 node " ++ (show rerootIndex) ++ "\n" ++ prettyIndices inGraph) -- empty


    -- else if null componentWithOutgroup then inGraph  -- error ("Error rooting wierdness in rerootTree " ++ (show rerootIndex) ++ "\n" ++ prettyIndices inGraph) -- empty

    -- check if new outtaxon has a parent--shouldn't happen-but could if its an internal node reroot
    else if null parentNewRootList || (True `elem` parentRootList) then inGraph
                                                              -- else if null componentWithOutgroup then error ("Outgroup index " ++ show rerootIndex ++ " not found in graph")
    else
        -- trace ("RRT: " ++ (show (rerootIndex, inNewRoot, outNewRoot))) ( 
        --reroot component with new outtaxon
        let -- componentWithNewOutgroup = snd $ head componentWithOutgroup
            -- (_, originalRootList) =  unzip $ filter ((==True).fst) $ zip (fmap (isRoot inGraph) componentWithNewOutgroup) componentWithNewOutgroup
            -- numRoots = 1 -- length originalRootList
            orginalRoot = orginalRootIndex -- head originalRootList
            originalRootEdges = (out inGraph orginalRoot)

        in

        {-
        if numRoots == 0 then error ("No root in rerootDisplayTree: Attempting to reroot on edge to node " ++ (show (orginalRoot,rerootIndex)) ++ prettyIndices inGraph) --empty

        -- check if outgroup in a multirooted component
        -- if wagner build this is ok
        -- else if numRoots > 1 then inGraph -- error ("Error: Attempting to reroot multi-rooted component") -- inGraph
        else
        -}
          --reroot graph safely automatically will only affect the component with the outgroup
          -- delete old root edge and create two new edges from oringal root node.
          -- keep orignl root node and delte/crete new edges when they are encounterd
          --trace ("Moving root from " ++ (show orginalRoot) ++ " to " ++  (show rerootIndex)) (
          let leftChildEdge = (orginalRoot, rerootIndex, edgeLabel $ head originalRootEdges)
              rightChildEdge = (orginalRoot, fst3 newRootOrigEdge, edgeLabel $ last originalRootEdges)

              --  this assumes 2 children of old root -- shouled be correct as Phylogenetic Graph
              newEdgeOnOldRoot = if (length originalRootEdges) /= 2 then error ("Number of root out edges /= 2 in rerootGraph: " ++ (show $ length originalRootEdges)
                ++ " root index: " ++ (show (orginalRoot, rerootIndex)) ++ "\nGraph:\n" ++ (prettyIndices inGraph))
                                 else (snd3 $ head originalRootEdges, snd3 $ last originalRootEdges, thd3 $ head originalRootEdges)

              newRootEdges = [leftChildEdge, rightChildEdge, newEdgeOnOldRoot]
              newGraph = insEdges newRootEdges $ delLEdges (newRootOrigEdge : originalRootEdges) inGraph

              -- get edges that need reversing
              newGraph' = preTraverseAndFlipEdgesTree orginalRootIndex [leftChildEdge,rightChildEdge] newGraph


              {- Check for valid tree
              cyclicString' = if cyclic newGraph' then " Is cyclic "
                          else  " Not cyclic "

              parentInCharinString'= if parentInChain newGraph' then " Is Parent in Chain "
                                  else " Not Parent in Chain "

              duplicatedsEdgeString' = if not (hasDuplicateEdge newGraph') then " No duplicate edges "
                                       else let dupEdgeList' = getDuplicateEdges newGraph'
                                            in
                                            (" Has duplicate edges: " ++ (show dupEdgeList') ++ "\nDeleting " ++ (show $ fmap toEdge $ (newRootOrigEdge : originalRootEdges)) ++ "\nInserting " ++ (show $ fmap toEdge $ newRootEdges)) 
              -}
          in
          --trace ("=")
          --trace ("In " ++ (GFU.showGraph inGraph) ++ "\nNew " ++  (GFU.showGraph newGraph) ++ "\nNewNew "  ++  (GFU.showGraph newGraph'))

          -- trace ("Deleting " ++ (show $ fmap toEdge (newRootOrigEdge : originalRootEdges)) ++ "\nInserting " ++ (show $ fmap toEdge newRootEdges) 
          --   ++ "\nRRT Out: " ++ cyclicString' ++ parentInCharinString' ++ duplicatedsEdgeString') (

          -- if cyclicString' == " Is cyclic " then inGraph 
          -- else if hasDuplicateEdge newGraph' then removeDuplicateEdges newGraph'
          -- else 
          {-Cycle check
          if cyclic newGraph' then 
            trace ("Orignal root: " ++ (show orginalRootIndex) ++ "New root: " ++ (show rerootIndex) ++ " Deleting " ++ (show $ fmap toEdge (newRootOrigEdge : originalRootEdges)) ++ "\nInserting " ++ (show $ fmap toEdge newRootEdges) 
              ++ "\nOrigGraph: " ++ (prettyIndices inGraph) ++ "\nNewGraph: " ++ (prettyIndices newGraph)++ "\nNewGraph': " ++ (prettyIndices newGraph'))
            empty
          else newGraph'-}
          newGraph'
          -- )
          -- ) -- )


-- | preTraverseAndFlipEdgesTree traverses a tree from starting edge flipping in-edges since they should
-- be out-edges 
-- when recursion its edges that don't need to be fliped then stops
-- assumes input edge is directed correctly
-- follows  traversal out "pre" order updating graph as edges flipped
preTraverseAndFlipEdgesTree :: (Eq b) => Node -> [LEdge b] ->  Gr a b -> Gr a b
preTraverseAndFlipEdgesTree rootIndex inEdgeList inGraph  =
  if null inEdgeList then inGraph
  else
    let -- first edge directled correctly
        inEdge@(_,v,_) = head inEdgeList

        -- edges "in" to child node of first edge--these should be out and need to be flipped
        childEdges = filter ((/= rootIndex) . fst3) $ filter (/= inEdge) $ inn inGraph v
        
        -- flip to "in" to "out" edges
        flippedEdges = fmap flipLEdge childEdges

        -- -- modify graph accordingly
        newGraph = insEdges flippedEdges $ delLEdges childEdges inGraph
    in
    --trace ("PTFE: flipped " ++ (show $ fmap toEdge flippedEdges)) (
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
preTraverseAndFlipEdges :: (Eq b) => [LEdge b] ->  Gr a b -> Gr a b
preTraverseAndFlipEdges inEdgelist inGraph  =
  if null inEdgelist then inGraph
  else
    let inEdge@(_,v,_) = head inEdgelist
        childEdges = (out inGraph v) ++ (filter (/= inEdge) $ inn inGraph v)

        -- returns list of edges that had to be flipped
        edgesToFlip = getToFlipEdges v childEdges 
        flippedEdges = fmap flipLEdge edgesToFlip
        newGraph = insEdges flippedEdges $ delLEdges edgesToFlip inGraph
    in
    --trace ("PTFE: flipped " ++ (show $ fmap toEdge flippedEdges)) (
    -- edge terminates in leaf or edges in correct orientation
    if null childEdges  || null edgesToFlip then preTraverseAndFlipEdges (tail inEdgelist) inGraph
    -- edge needs to be reversed to follow through its children from a new graph
    else preTraverseAndFlipEdges (flippedEdges ++ (tail inEdgelist)) newGraph
    -- )

-- | getToFlipEdges takes an index and check edge list
-- and creates new list of edges that need to be flipped
getToFlipEdges ::  Node -> [LEdge b] -> [LEdge b]
getToFlipEdges parentNodeIndex inEdgeList =
  if null inEdgeList then []
  else
    let firstEdge@(u,_,_) = head inEdgeList
    in
    if parentNodeIndex /= u then firstEdge : getToFlipEdges parentNodeIndex (tail inEdgeList)
    else getToFlipEdges parentNodeIndex (tail inEdgeList)

-- | Random generates display trees up to input number by choosing
-- to keep indegree nodes > 1 unifomaly at random
generateDisplayTreesRandom :: (Show a, Show b, Eq a, Eq b, NFData a, NFData b) => Int -> Int -> Gr a b -> [Gr a b]
generateDisplayTreesRandom rSeed numDisplayTrees inGraph =
  if isEmpty inGraph then error "Empty graph in generateDisplayTreesRandom"
  else if isTree inGraph then [inGraph]
  else
    let atRandomList = take numDisplayTrees $ randomIntList rSeed
        randDisplayTreeList = PU.seqParMap rdeepseq (randomlyResolveGraphToTree inGraph) atRandomList -- `using` PU.myParListChunkRDS
    in
    randDisplayTreeList

-- | randomlyResolveGraphToTree resolves a single graph to a tree by choosing single indegree edges
-- uniformly at random and deleting all others from graph
-- in=out=1 nodes are contracted, HTU's withn outdegree 0 removed, graph reindexed
-- but not renamed--edges from root are left alone.
randomlyResolveGraphToTree :: (Show a, Show b, Eq a, Eq b) => Gr a b -> Int -> Gr a b
randomlyResolveGraphToTree inGraph randVal =
  if isEmpty inGraph then error "Empty graph in randomlyResolveGraphToTree"
  else
    let (_, leafList, _, _) = splitVertexList inGraph
        -- rootEdgeList = fmap (out inGraph) $ fmap fst rootList
        inEdgeListByVertex = (fmap (inn inGraph) (nodes inGraph)) -- L.\\ rootEdgeList
        randList = fmap abs $ randomIntList randVal
        edgesToDelete = concat $ zipWith  chooseOneDumpRest randList  (fmap (fmap toEdge) inEdgeListByVertex)
        newTree = delEdges edgesToDelete inGraph
        newTree' = removeNonLeafOut0Nodes leafList newTree
        newTree'' = contractIn1Out1Edges newTree'
        reindexTree = reindexGraph newTree''
    in
    -- trace ("RRGT\n" ++ (prettify inGraph) ++ "\n to delete " ++ (show edgesToDelete) ++ "\nNew graph:\n" ++ (prettify newTree)
    --   ++ "\nnewTree'\n" ++ (prettify newTree') ++ "\nnewTree''\n" ++ (prettify newTree'') ++ "\reindex\n" ++ (prettify reindexTree))
    reindexTree

-- | chooseOneDumpRest takes random val and chooses to keep the edge in list returning list of edges to delete
chooseOneDumpRest :: Int -> [Edge] -> [Edge]
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
generateDisplayTrees :: (Eq a) => Bool -> Gr a b -> [Gr a b]
generateDisplayTrees contractEdges inGraph =
    if isTree inGraph then [inGraph]
    else 
        let (_, leafList, _, _) = splitVertexList inGraph
        in
        generateDisplayTrees' contractEdges leafList [inGraph] []

-- | generateDisplayTrees' takes a graph list and recursively generates
-- a list of trees created by progresively resolving each network vertex into a tree vertex
-- in each input graph
-- creating up to 2**m (m network vertices) trees.
-- call -> generateDisplayTrees'  [startGraph] []
-- the second and third args contain graphs that need more work and graphs that are done (ie trees)
generateDisplayTrees' :: (Eq a) => Bool -> [LNode a] -> [Gr a b] -> [Gr a b] -> [Gr a b]
generateDisplayTrees' contractEdges leafList curGraphList treeList  =
  if null curGraphList then
      let treeList' = fmap (removeNonLeafOut0Nodes leafList) treeList
          treeList'' = if contractEdges then fmap contractIn1Out1Edges treeList'
                       else treeList'
          reindexedTreeList = fmap reindexGraph treeList''
      in
      reindexedTreeList

  else
    let firstGraph = head curGraphList
    in
      if isEmpty firstGraph then []
      else
        let nodeList = labNodes firstGraph
            inNetEdgeList = filter ((>1).length) $ fmap (inn firstGraph) $ fmap fst nodeList
        in
        if null inNetEdgeList then generateDisplayTrees' contractEdges leafList (tail curGraphList) (firstGraph : treeList)
        else
          let newGraphList = splitGraphListFromNode inNetEdgeList [firstGraph]
          in
          generateDisplayTrees' contractEdges leafList (newGraphList ++ (tail curGraphList)) treeList

-- | splitGraphListFromNode take a graph and a list of edges for indegree > 1 node
-- removes each in edge in turn to create a new graph and maintains any in 1 out 1 nodes
-- if these edges are to be contracted out  use 'contractOneOneEdges'
-- should be a single traversal 'splitting' graph each time. removing edges anmd recursing
-- to do it again untill all edges are indegree 1.
-- the edges in list for recurvise deleteion should always be in all graphs uuntil
-- the end
splitGraphListFromNode :: [[LEdge b]] -> [Gr a b] -> [Gr a b]
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
-- recursively keeps the index-th edge and deletes the others creating a new graph list
deleteEdgesCreateGraphs :: [([LEdge b], Int)] -> Int -> Gr a b -> [Gr a b]
deleteEdgesCreateGraphs netEdgeIndexPairList counter inGraph =
  if isEmpty inGraph then error "Empty graph in  "deleteEdgesCreateGraphs
  else if null netEdgeIndexPairList then []
  else
    let (edgeList, lIndex) = head netEdgeIndexPairList
        --edgeToKeep = edgeList !! index
        edgesToDelete = (take lIndex edgeList) ++ (drop (lIndex + 1) edgeList)
        newGraph = delLEdges edgesToDelete inGraph
    in
    newGraph : deleteEdgesCreateGraphs (tail netEdgeIndexPairList) (counter + 1) inGraph

-- | undirectedEdgeEquality checks edgse for equality irrespective of direction
undirectedEdgeEquality :: Edge -> Edge -> Bool
undirectedEdgeEquality (a,b) (c,d) = if a == c && b == d then True
                                               else if a == d && b == c then True
                                               else False

-- | undirectedEdgeMinus subtracts edges in the second list from those in the first using
-- undirected matching
undirectedEdgeMinus :: [Edge] -> [Edge] -> [Edge]
undirectedEdgeMinus firstList secondList =
    if null firstList then []
    else
        let firstEdge@(a,b) = head firstList
        in
        if firstEdge `L.elem` secondList then undirectedEdgeMinus (tail firstList) secondList
        else if (b,a) `L.elem` secondList then undirectedEdgeMinus (tail firstList) secondList
        else firstEdge : undirectedEdgeMinus (tail firstList) secondList
