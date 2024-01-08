 {- |
Module      :  Eun.hs
Description :  Module to calculate various graph reconciliation methods Wheeler (2021)
               input graphviz dot files and newick
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

module Reconciliation.Eun ( reconcile
                          , makeProcessedGraph
                          , addGraphLabels)
                          where

import PHANE.Evaluation
import PHANE.Evaluation.Verbosity (Verbosity (..))
import Control.Monad (when)
import Control.Monad.IO.Class (MonadIO (..))
import Control.Parallel.Strategies
import Data.BitVector    qualified                as BV
import Data.Bits   qualified                      as B
import Data.Graph.Inductive.Graph   qualified     as G
import Data.Graph.Inductive.PatriciaTree qualified as P
import Data.Graph.Inductive.Query.BFS  qualified  as BFS
import Data.GraphViz                     as GV
import Data.GraphViz.Printing
import Data.List     qualified                    as L
import Data.Map.Strict   qualified                as Map
import Data.Maybe
import Data.Set       qualified                   as S
import Data.Text.Lazy  qualified                  as T
import Data.Vector    qualified                   as V
import GraphFormatUtilities   qualified           as PhyP
import Graphs.GraphOperations qualified           as GO
import Reconciliation.Adams   qualified           as A
import Types.Types
import Utilities.LocalGraph  qualified            as LG
--import           Debug.Trace
--import           ParallelUtilities                 as PU


{-
-- | turnOnOutZeroBit turns on the bit 'nleaves" signifying that
-- the node is outdegree 1
-- this so outdegree one nodes and their child have differnet bit sets
turnOnOutZeroBit :: BV.BV -> Int -> BV.BV
turnOnOutZeroBit inBitVect nLeaves = BV.or [B.bit nLeaves, inBitVect]
-}

{-
-- | turnOffOutZeroBit turns off the bit 'nleaves" signifying that
-- the node is outdegree /= 1
-- this so outdegree one nodes and their child have differnet bit sets
turnOffOutZeroBit :: BV.BV -> Int -> BV.BV
turnOffOutZeroBit inBitVect nLeaves = BV.extract (nLeaves - 1) 0 inBitVect
-}

-- | setOutDegreeOneBits assumes outdegree of vertex = 1, takes the number of leaves, bitvector
-- repersentation (rawBV) of vertex, a L.sorted list of outdegree=1 vertices and the vertex index
-- and creates a bitvector to prepend of length the number of outdgree=1 vertices where
-- the correpsonding vertex index position in ilist is 'On' and the remainder are 'Off'
-- this insures a unique lablelling for all outdegree=1 vertices
setOutDegreeOneBits :: BV.BV -> [Int] -> Int -> BV.BV
setOutDegreeOneBits inBitVect out1VertexList vertexIndex =
  if null out1VertexList then error "Empty outdegree=1 vertex list in setOutDegreeOneBits"
  else
    let vertListIndex = L.elemIndex vertexIndex out1VertexList
        vertexPosition = fromJust vertListIndex
        boolList = replicate vertexPosition False <> [True] <> replicate (length out1VertexList - vertexPosition - 1) False
        prependBitVect = BV.fromBits boolList
    in
    if isNothing vertListIndex then error ("Vertex " <> show vertexIndex <> " not found in list " <> show out1VertexList)
    else BV.append prependBitVect inBitVect

-- | getRoot takes a graph and list of nodes and returns vertex with indegree 0
-- so assumes a connected graph--with a single root--not a forest
-- this does noit include unconnected leaves
getRoots :: P.Gr a b -> [G.Node] -> [Int]
getRoots inGraph nodeList =
  if null nodeList then []
  else
    let firstNode = head nodeList
    in
    if (G.indeg inGraph firstNode == 0) && (G.outdeg inGraph firstNode > 0) then firstNode : getRoots inGraph (tail nodeList)
    else getRoots inGraph (tail nodeList)


-- | getUnConnectedLeaves takes a graph and list of nodes and returns vertex with indegree 0
-- and outdegree == 0
getUnConnectedLeaves :: P.Gr a b -> [G.Node] -> [Int]
getUnConnectedLeaves inGraph nodeList =
  if null nodeList then []
  else
    let firstNode = head nodeList
    in
    if (G.indeg inGraph firstNode == 0) && (G.outdeg inGraph firstNode == 0) then firstNode : getUnConnectedLeaves inGraph (tail nodeList)
    else getUnConnectedLeaves inGraph (tail nodeList)



-- | getUnconmnectedNOdes takes a graph and list of nodes and returns vertex with indegree 0
-- and outdegeee 0
getUnConnectedNodes :: P.Gr String String -> Int -> [G.Node] -> [G.LNode BV.BV]
getUnConnectedNodes inGraph nLeaves nodeList =
  if null nodeList then []
  else
    let firstNode = head nodeList
        newNode =  (firstNode, B.bit firstNode)
    in
    if G.deg inGraph firstNode == 0 then
      newNode : getUnConnectedNodes inGraph nLeaves (tail nodeList)
    else getUnConnectedNodes inGraph nLeaves (tail nodeList)


-- | makeNodeFromChildren gets bit vectors as union of children in a post order traversal from leaves
-- the prepending of a single 'On' bit if there is only once child (setOutDegreeOneBit)
-- is modified to allow for multiple outdegree 1 vertices as parent of single vertex
makeNodeFromChildren :: P.Gr String String -> Int -> V.Vector (G.LNode BV.BV) -> [Int] -> Int -> PhyG [G.LNode BV.BV]
makeNodeFromChildren inGraph nLeaves leafNodes out1VertexList myVertex =
  if myVertex < nLeaves then pure [leafNodes V.! myVertex]
  else
      let myChildren = G.suc inGraph myVertex
          -- parallel
          action :: Int -> PhyG [G.LNode BV.BV]
          action = makeNodeFromChildren inGraph nLeaves leafNodes out1VertexList

      in do
          myChildrenNodes <- getParallelChunkTraverse >>= \pTraverse ->
               action `pTraverse` myChildren
            -- PU.seqParMap PU.myStrategyRDS (makeNodeFromChildren inGraph nLeaves leafNodes out1VertexList) myChildren -- `using` PU.myParListChunkRDS

          let rawBV = BV.or $ fmap (snd . head) myChildrenNodes
          let myBV = if length myChildren /= 1 then rawBV
                     else setOutDegreeOneBits rawBV out1VertexList myVertex
      
          pure $ (myVertex, myBV) : concat myChildrenNodes

-- | getNodesFromARoot follows nodes connected to a root.
-- can be fmapped over roots to hit all--should be ok if multiple hits on nodes
-- since all labeled by BV.BVs  need to fuse them if multiple roots to make sure nodes are consistent
-- and only one per root--should be ok for multikple jhists of nodes since BVs are from childre
-- just wasted work.  Should L.nub after to maeksure only unique (by BV) nodes in list at end
getNodesFromARoot :: P.Gr String String -> Int -> [G.LNode BV.BV] -> Int -> PhyG [G.LNode BV.BV]
getNodesFromARoot inGraph nLeaves leafNodes rootVertex =
  if  G.isEmpty inGraph then error "Input graph is empty in getLabelledNodes"
  else
    let rootChildVerts = G.suc inGraph rootVertex

        -- get outdree = 1 node list for creting prepended bit vectors
        out1VertexList = L.sort $ filter ((==1).G.outdeg inGraph) $ G.nodes inGraph

         -- parallel
        action :: Int -> PhyG [G.LNode BV.BV]
        action = makeNodeFromChildren inGraph nLeaves (V.fromList leafNodes) out1VertexList
    in do
        -- recurse to children since assume only leaves can be labbeled with BV.BVs
        -- fmap becasue could be > 2 (as in at root)
        rootChildNewNodes <- getParallelChunkTraverse >>= \pTraverse ->
            action `pTraverse` rootChildVerts
            -- rootChildNewNodes = PU.seqParMap PU.myStrategyRDS (makeNodeFromChildren inGraph nLeaves (V.fromList leafNodes) out1VertexList) rootChildVerts -- `using` PU.myParListChunkRDS

        -- check if outdegree = 1
        let rawBV = BV.or $ fmap (snd . head) rootChildNewNodes
        let rootBV = if length rootChildVerts /= 1 then rawBV
                     else setOutDegreeOneBits rawBV out1VertexList rootVertex
    
        pure $ (rootVertex, rootBV) : concat rootChildNewNodes

-- | getLabelledNodes labels nodes with bit vectors union of subtree leaves via post order traversal
-- adds nodes to reDoneNodes as they are preocessed
-- reorder NOdes is n^2 should be figured out how to keep them in order more efficeintly
getLabelledNodes :: P.Gr String String -> Int -> [G.LNode BV.BV] -> PhyG [G.LNode BV.BV]
getLabelledNodes inGraph nLeaves leafNodes  =
  -- trace ("getLabbeled graph with " <> (show $ G.noNodes inGraph) <> " nodes in " <> (showGraph inGraph)) (
  if  G.isEmpty inGraph then error "Input graph is empty in getLabelledNodes"
  else
    let rootVertexList = getRoots inGraph (G.nodes inGraph)

    in do
        htuList' <- mapM (getNodesFromARoot inGraph nLeaves leafNodes) rootVertexList
        let htuList = L.nub $ concat htuList'
    
         -- this for adding in missing data
        let unConnectedNodeList = getUnConnectedNodes inGraph nLeaves (G.nodes inGraph)
        
        pure $ reorderLNodes (htuList <> unConnectedNodeList)  0


-- | findLNode takes an index and looks for node with that as vertex and retuirns that node
findLNode :: Int -> [G.LNode BV.BV] -> G.LNode BV.BV
findLNode vertex lNodeList =
  if null lNodeList then error ("Node " <> show vertex <> " not found")
  else
      let (a,b) = head lNodeList
      in
      if a == vertex then (a,b)
      else findLNode vertex (tail lNodeList)

-- | reorderLNodes takes a list of nodes and reorders and order based on node vertex number
-- n^2 ugh
reorderLNodes :: [G.LNode BV.BV]  -> Int -> [G.LNode BV.BV]
reorderLNodes inNodeList inIndex
  | null inNodeList = []
  | inIndex == length inNodeList = []
  | otherwise =
    let newNode =  findLNode inIndex inNodeList
    in
    newNode : reorderLNodes inNodeList (inIndex + 1)
   -- )

-- | relabelEdgs creates (BV.BV, BV.BV) labnels for an edges
relabelEdge :: V.Vector (G.LNode BV.BV) -> G.LEdge String -> G.LEdge (BV.BV, BV.BV)
relabelEdge allNodesVect inLEdge =
  let (e,u,_) = inLEdge
      eNodeBV = snd (allNodesVect V.! e)
      uNodeBV = snd (allNodesVect V.! u)
  in
  (e,u,(eNodeBV,uNodeBV))

-- | changeLabelEdge labels edges by descendent vertex label
-- assumes leaves first then vertices labeled in order
-- offset for numLeaves assumes only labeling non-leaves so smaller set
-- assumes that if an "urroot" has been added via `union` then it is last vertex
-- this for condition where root has been added to majority consensus tree
changeLabelEdge :: Int -> V.Vector Double -> [G.LEdge b] -> [G.LEdge Double]
changeLabelEdge numLeaves freqVect edgeList =
  if null edgeList then []
  else
    let (e,u,_) =head edgeList
        newLabel
          | u < numLeaves = 1
          | (u - numLeaves) >= V.length freqVect = 1
          | otherwise = freqVect V.! (u - numLeaves)
    in
    --trace (show (e,u) <> " " <> (show (u - numLeaves)) <> " " <> show freqVect) -- <> " " <> show newLabel)
    (e,u, newLabel) : changeLabelEdge numLeaves freqVect (tail edgeList)

-- | addEdgeFrequenciesToGraph takes a greaph and edge frequencies and relables edges
-- with node frequencies of discedendet node
addEdgeFrequenciesToGraph :: P.Gr a b -> Int -> [Double] -> P.Gr a Double
addEdgeFrequenciesToGraph inGraph numLeaves freqList =
  let inNodes = G.labNodes inGraph
      inEdges = G.labEdges inGraph
      newEdges = changeLabelEdge numLeaves (V.fromList freqList) inEdges
  in
  --trace (show inEdges)
  G.mkGraph inNodes newEdges

-- | getLeafNumber take Graph and gets nu,ber of leaves (outdegree = 0)
getLeafNumber :: P.Gr BV.BV (BV.BV, BV.BV) -> Int
getLeafNumber inGraph =
  let degOutList = G.outdeg inGraph <$> G.nodes inGraph
  in length $ filter (==0) degOutList

{-
-- | findStrLabel checks Attributes (list f Attribute) from Graphvz to extract the String label of node
-- returns Maybe Text
findStrLabel :: Attributes -> Maybe T.Text
findStrLabel = getFirst . foldMap getStrLabel


-- | getStrLabel takes an Attribute and reurns Text if StrLabel found, mempty otherwise
getStrLabel :: Attribute -> First T.Text
getStrLabel (Label (StrLabel txt)) = First . Just $ txt
getStrLabel _                      = mempty


-- | getLeafString takes a pairs (node vertex number, graphViz Attributes)
-- and returns String name of leaf of Stringified nude number if unlabbeled
getLeafString :: (Int, Attributes) -> String
getLeafString (nodeIndex, nodeLabel) =
  let maybeTextLabel = findStrLabel nodeLabel
  in
  maybe (show nodeIndex) T.unpack maybeTextLabel

-- | getLeafList returns leaf complement of graph from DOT file
getLeafList ::  P.Gr Attributes Attributes -> [G.LNode String]
getLeafList inGraph =
  if G.isEmpty inGraph then []
  else
    let degOutList = G.outdeg inGraph <$> G.nodes inGraph
        newNodePair = zip degOutList (G.labNodes inGraph)
        leafPairList = filter ((==0).fst ) newNodePair
        (_, leafList) = unzip leafPairList
        (nodeVerts, _) = unzip leafList
        newLabels = fmap getLeafString leafList
        leafList' = zip nodeVerts newLabels
    in
    leafList'
-}

-- | getLeafListNewick returns leaf complement of graph from newick file
-- difference from above is in the leaf label type
getLeafListNewick ::  P.Gr a b -> [G.LNode a]
getLeafListNewick inGraph =
  if G.isEmpty inGraph then []
  else
    let degOutList = G.outdeg inGraph <$> G.nodes inGraph
        newNodePair = zip degOutList (G.labNodes inGraph)
        leafPairList = filter ((==0).fst ) newNodePair
        (_, leafList) = unzip leafPairList
        (nodeVerts, _) = unzip leafList
        -- only different line
        newLabels = fmap snd leafList
        leafList' = zip nodeVerts newLabels
    in
    leafList'

{-
-- | checkNodesSequential takes a list of nodes and returns booolean
-- True if nodes are input with sequential numerical indices
-- False if not--screws up reindexing later which assumes they are successive
checkNodesSequential :: G.Node -> [G.Node] -> Bool
checkNodesSequential prevNode inNodeList
  | null inNodeList = True
  | (head inNodeList - prevNode) /= 1 = trace ("Index or indices missing between " <> (show $ head inNodeList) <> " and " <> (show prevNode))  False
  | otherwise = checkNodesSequential (head inNodeList) (tail inNodeList)
-}

-- | reAnnotateGraphs takes parsed graph input and reformats for EUN
reAnnotateGraphs :: P.Gr String String -> PhyG (P.Gr BV.BV (BV.BV, BV.BV))
reAnnotateGraphs inGraph =
  -- trace ("Reannotating " <> (showGraph inGraph)) (
  if G.isEmpty inGraph then error "Input graph is empty in reAnnotateGraphs"
  else
    let degOutList = G.outdeg inGraph <$> G.nodes inGraph
        nLeaves = length $ filter (==0) degOutList
        leafVerts = [0..(nLeaves - 1)]
        leafIntegers = fmap B.bit leafVerts
        leafBitVects =  leafIntegers  -- fmap (BV.bitVec nLeaves) leafIntegers
        leafNodes = Prelude.zip leafVerts leafBitVects
    in do

        allNodes <- getLabelledNodes inGraph nLeaves leafNodes
        let allEdges = fmap (relabelEdge (V.fromList allNodes)) (G.labEdges inGraph)
        -- assign HTU BV via postorder pass.
        pure $ G.mkGraph allNodes allEdges

-- | checkBVs looks at BV.BV of node and retuns FALSE if found True if not
checkBVs :: BV.BV -> [G.LNode BV.BV] -> Bool
checkBVs inBV nodeList =
  null nodeList || (
  let (_, bv) = head nodeList
  in
  inBV /= bv && checkBVs inBV (tail nodeList))

-- | checkBVs looks at BV.BV of node and retuns FALSE if found True if not
checkEdgeBVs :: (BV.BV, BV.BV) -> [G.LEdge (BV.BV, BV.BV)] -> Bool
checkEdgeBVs (inABV, inBBV) edgeList =
  null edgeList || (
  let (_, _, (aBV,bBV)) = head edgeList
  in
  not ((inABV == aBV) && (inBBV == bBV)) && checkEdgeBVs (inABV, inBBV) (tail edgeList))

-- | addAndReIndexUniqueNodes takes an inital list of nodes and adds new nodes reindexed
-- check identity by BV.BV
-- index bigger than size becasue starting after number of leaves
addAndReIndexUniqueNodes :: Int -> [G.LNode BV.BV] -> [G.LNode BV.BV] -> [G.LNode BV.BV]
addAndReIndexUniqueNodes newIndex nodesToExamine uniqueReIndexedNodes =
  if null nodesToExamine then uniqueReIndexedNodes
  else
    let (_, inBV) = head nodesToExamine
        isUnique = checkBVs inBV uniqueReIndexedNodes
    in
    if isUnique then
      let newNode = (newIndex, inBV)
      in
      addAndReIndexUniqueNodes (newIndex + 1) (tail nodesToExamine) (newNode : uniqueReIndexedNodes)
    else addAndReIndexUniqueNodes newIndex (tail nodesToExamine) uniqueReIndexedNodes


-- | getNodeIndex takes a BV.BV and returns a node with the same BV.BV
getNodeIndex :: BV.BV -> [G.LNode BV.BV] -> Int
getNodeIndex inBV nodeList =
  if null nodeList then error ("Node  with BV " <> show inBV <> " not found in getNodeIndex")
  else
    let (inIndex, bv) = head nodeList
    in
    if bv == inBV then inIndex
    else getNodeIndex inBV (tail nodeList)

-- | addAndReIndexEdges  takes list of indexed nodes and BV, a list of edges to examine and a list of edges to keep
-- checks for uniqueness of edges by BV.BVs on (e,u) and reindexes the edge nodes based on the node set with bit vectors
-- keep method either 'unique" or "all" to keep lists of unique or all edges
addAndReIndexEdges :: String -> [G.LNode BV.BV] -> [G.LEdge (BV.BV,BV.BV)] -> [G.LEdge (BV.BV,BV.BV)] -> [G.LEdge (BV.BV,BV.BV)]
addAndReIndexEdges keepMethod indexedNodes edgesToExamine uniqueReIndexedEdges =
  if null edgesToExamine then uniqueReIndexedEdges
  else
      let (_, _, (eBV, uBV)) = head edgesToExamine
          isUnique = checkEdgeBVs (eBV, uBV) uniqueReIndexedEdges
      in
      if (keepMethod == "all") || isUnique then
          -- Find nodes with BVs of edge
          let eNode = getNodeIndex eBV indexedNodes
              uNode = getNodeIndex uBV indexedNodes
              newEdge = (eNode, uNode, (eBV, uBV))
          in
          addAndReIndexEdges keepMethod indexedNodes (tail edgesToExamine) (newEdge : uniqueReIndexedEdges)
      else addAndReIndexEdges keepMethod indexedNodes (tail edgesToExamine) uniqueReIndexedEdges

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

-- | makeEUN take list of nodes and edges, deletes each edge (e,u) in turn makes graph,
-- checks for path between nodes e and u, if there is delete edge otherwise keep edge in list for new graph
makeEUN ::  (Eq b, NFData b) => [G.LNode a] -> [G.LEdge b] -> P.Gr a b -> PhyG (P.Gr a b)
makeEUN nodeList fullEdgeList fullGraph =
  let -- counterList = [0..(length fullEdgeList - 1)]
      -- requiredEdges = concat $ fmap (testEdge nodeList fullEdgeList) counterList
      -- parallel 
      -- action :: G.LEdge b -> [G.LEdge b]
      action = testEdge fullGraph
  in do
      testPar <- getParallelChunkMap
      let requiredEdges = testPar action fullEdgeList -- PU.seqParMap PU.myStrategyRDS (testEdge fullGraph) fullEdgeList -- `using` PU.myParListChunkRDS
      let newGraph = G.mkGraph nodeList (concat requiredEdges)
      pure newGraph

-- | getLeafLabelMatches tyakes the total list and looks for elements in the smaller local leaf set
-- retuns int index of the match or (-1) if not found so that leaf can be added in orginal order
getLeafLabelMatches ::[G.LNode String] -> G.LNode String -> (Int, Int)
getLeafLabelMatches localLeafList totNode =
  if null localLeafList then (-1, fst totNode)
  else
    let (inIndex, leafString) = head localLeafList
    in
    if snd totNode == leafString then (inIndex, fst totNode)
    else getLeafLabelMatches (tail localLeafList) totNode

-- | reIndexEdge takes an (Int, Int) map, labelled edge, and returns a new labelled edge with new e,u vertices
reIndexLEdge ::  Map.Map Int Int -> G.LEdge a -> G.LEdge String
reIndexLEdge vertexMap inEdge =
  if Map.null vertexMap then error "Null vertex map"
  else
    let (e,u,_) = inEdge
        newE = Map.lookup e vertexMap
        newU = Map.lookup u vertexMap
    in
    if isNothing newE then error ("Error looking up vertex " <> show e <> " in " <> show (e,u))
    else if isNothing newU then error ("Error looking up vertex " <> show u <> " in " <> show (e,u))
    else (fromJust newE, fromJust newU, "")

-- | reIndexAndAddLeaves takes rawGraphs and total input leaf sets and reindexes node, and edges, and adds
-- in leaves (with out edges) so that later processing can get bit vectors correct and match from
-- graph to graph.
-- new node set in teh total leaf set form all graphs plus teh local HTUs renumbered up based on added leaves
-- the map contains leaf mappings based on label of leaf, the HTUs extend that map with stright integers.
-- edges are re-indexed based on that map
reIndexAndAddLeavesEdges :: [G.LNode String] -> ([G.LNode String], P.Gr a b) -> PhyG (P.Gr String String)
reIndexAndAddLeavesEdges totallLeafSet (inputLeafList, inGraph) =
  if G.isEmpty inGraph then pure G.empty
  else
      -- reindex nodes and edges and add in new nodes (total leaf set + local HTUs)
      -- create a map between inputLeafSet and totalLeafSet which is the canonical enumeration
      -- then add in local HTU nodes and for map as well
      -- trace ("Original graph: " <> (showGraph inGraph)) (
      let --parallel
          action :: G.LNode String -> (Int, Int)
          action = getLeafLabelMatches inputLeafList
      in do
          labelPar <- getParallelChunkMap
          let correspondanceList = labelPar action totallLeafSet
            -- PU.seqParMap PU.myStrategyRDS (getLeafLabelMatches inputLeafList) totallLeafSet -- `using` PU.myParListChunkRDS
          let matchList = filter ((/=(-1)).fst) correspondanceList
          --remove order dependancey
          -- htuList = [(length inputLeafList)..(length inputLeafList + htuNumber - 1)]
          let htuList = fmap fst (G.labNodes inGraph) L.\\ fmap fst inputLeafList
          let htuNumber =  length (G.labNodes inGraph) - length inputLeafList
          let newHTUNumbers = [(length totallLeafSet)..(length totallLeafSet + htuNumber - 1)]
          let htuMatchList = zip htuList newHTUNumbers
          let vertexMap = Map.fromList (matchList <> htuMatchList)
          let reIndexedEdgeList = fmap (reIndexLEdge vertexMap) (G.labEdges inGraph)

          let newNodeNumbers = [0..(length totallLeafSet + htuNumber - 1)]
          let attributeList = replicate (length totallLeafSet + htuNumber) "" -- origAttribute
          let newNodeList = zip newNodeNumbers attributeList
      
          pure $ G.mkGraph newNodeList reIndexedEdgeList

-- | relabelNode takes nofde list and labels leaves with label and HTUs with String of HexCode of BV label
relabelNodes :: [G.LNode BV.BV] -> [G.LNode String] -> [G.LNode String]
relabelNodes inNodes leafLabelledNodes
  | null inNodes = []
  | not $ null leafLabelledNodes = head leafLabelledNodes : relabelNodes (tail inNodes) (tail leafLabelledNodes)
  | otherwise =
  let (vertex, _) = head inNodes
  in
  (vertex, "HTU" <> show vertex) : relabelNodes (tail inNodes) []

-- | addGraphLabels take Graph and changes to add nodes labelled wiyth String, edges as well
addGraphLabels :: P.Gr BV.BV (BV.BV, BV.BV) -> [G.LNode String] -> P.Gr String String
addGraphLabels inGraph totallLeafSet
  | G.isEmpty inGraph = error "Empty graph in addGraphLabels"
  | null totallLeafSet = error "Empty leaf set in addGraphLabels"
  | otherwise =
  let newNodes = relabelNodes (G.labNodes inGraph) totallLeafSet
    -- newNodes = totallLeafSet <> newHTUList
      (eList, uList) = unzip (G.edges inGraph)

      newEdges = zip3 eList uList (replicate (length eList) "")
  in
  -- trace ("Relabelled EUN : " <> (showGraph $ G.mkGraph newNodes newEdges) <> " from " <> (show totallLeafSet))
  G.mkGraph newNodes newEdges

-- | intermediateNodeExists takes two node bitvetors and the full bitvector list
-- and checks to see if there is an intermediate node  between first two that would
-- remove need for the edge between the first two.
-- This should reduce time complexity of vertex-based reconciliation to O(n^3) from O(n^4)
intermediateNodeExists :: BV.BV -> BV.BV -> [BV.BV] -> Bool
intermediateNodeExists aBV cBV fullNodeBVList =
  not (null fullNodeBVList) &&
  let bBV = head fullNodeBVList
      leftIntersection = BV.and [aBV, bBV]
      rightIntersection = BV.and [bBV, cBV]
  in
  if (bBV == aBV) || (bBV == cBV) then intermediateNodeExists aBV cBV (tail fullNodeBVList)
  else ((leftIntersection == bBV) && (rightIntersection == cBV)) || intermediateNodeExists aBV cBV (tail fullNodeBVList)

-- | getIntersectionEdges takes a node A and cretes directed edges to each other edge in [B]
-- with rulkesLEdge
--  if A intesect B = empty then no edge
--  else if A intesect B = B then create edge A->B
--  else if A intesect B = A then create edge B->A
--  else --can't happen
-- Added in check for intermediate (by bitvector) node that shold obviate need for
-- breadth first search for vertex-based reconciliation
--  if A > B and B > C and A intersect B = B and B intersect C = C
--    then edge  A->C is redundant and is not added to edge set
getIntersectionEdges ::[BV.BV] -> [G.LNode BV.BV] -> G.LNode BV.BV -> [G.LEdge (BV.BV,BV.BV)]
getIntersectionEdges fullNodeBVList bNodeList aNode =
  if null bNodeList then []
  else
      let (aIndex, aBV) = aNode
          (bIndex, bBV) = head bNodeList
          intersection = BV.and [aBV, bBV]
      in
      -- only do the directed 1/2 so no L.nub issues later
      if (bBV >= aBV) || (intersection == 0) then getIntersectionEdges fullNodeBVList (tail bNodeList) aNode
      else if intersection == bBV then
        if intermediateNodeExists aBV bBV fullNodeBVList then getIntersectionEdges fullNodeBVList (tail bNodeList) aNode
        else (aIndex, bIndex, (aBV, bBV)) : getIntersectionEdges fullNodeBVList (tail bNodeList) aNode
      else  getIntersectionEdges fullNodeBVList (tail bNodeList) aNode

-- | combinable tales a list of bitvecotrs and a single bitvector
-- and checks each of the first to see if combinable
-- if A and B == A,B, or 0 then True else False
-- if True return [bitvector] else []  if not
combinable :: String -> [BV.BV] -> BV.BV -> [BV.BV]
combinable comparison bvList bvIn
  | comparison == "identity" =
    if null bvList then []
    else [bvIn | bvIn `elem` bvList]
  | comparison == "combinable" = -- combinable sensu Nelson 1979
    if null bvList then [bvIn]
    else
      let intersectList = fmap (checkBitVectors bvIn) bvList -- took out paralleism here
          isCombinable = L.foldl' (&&) True intersectList
      in
      [bvIn | isCombinable]
  | otherwise = errorWithoutStackTrace ("Comparison method " <> comparison <> " unrecongnized (combinable/identity)")
      where checkBitVectors a b
              = let c = BV.and [a, b] in
                  c == a || c == b || c == 0

-- | getGraphCompatibleList takes a list of graphs (list of node Bitvectors)
-- and retuns a list of each graph a bitvector node is compatible with
-- this isued later for majority rule consensus
-- each bit vector node will have a list of length 1..number of graphs
getGraphCompatibleList :: String -> [[BV.BV]] -> BV.BV-> [BV.BV]
getGraphCompatibleList comparison inBVListList bvToCheck =
  if null inBVListList then error "Null list of list of bitvectors in getGraphCompatibleList"
  else
    let compatibleList = concatMap (flip (combinable comparison) bvToCheck) inBVListList
    in
    -- trace (show $ length compatibleList)
    compatibleList

-- | getCompatibleList takes a list of graph node bitvectors as lists
-- retuns a list of lists of bitvectors where the length of the list of the individual bitvectors
-- is the number of graphs it is compatible with
getCompatibleList :: String -> [[BV.BV]] -> [[BV.BV]]
getCompatibleList comparison inBVListList =
  if null inBVListList then error "Null list of list of bitvectors in getCompatibleList"
  else
    let uniqueBVList = L.nub $ concat inBVListList
        bvCompatibleListList = fmap (getGraphCompatibleList comparison inBVListList) uniqueBVList
    in
    filter (not . null) bvCompatibleListList

-- | getThresholdNodes takes a threshold and keeps those unique objects present in the threshold percent or
-- higher.  L.sorted by frequency (low to high)
-- urRoot added to make sure there will be a single connected graph
getThresholdNodes :: String -> Int -> Int -> [[G.LNode BV.BV]] -> ([G.LNode BV.BV], [Double])
getThresholdNodes comparison thresholdInt numLeaves objectListList
  | thresholdInt < 0 || thresholdInt > 100 = errorWithoutStackTrace "Threshold must be in range [0,100]"
  | null objectListList = error "Empty list of object lists in getThresholdObjects"
  | otherwise =
    let numGraphs = fromIntegral $ length objectListList
        indexList = [numLeaves..(numLeaves + length objectGroupList - 1)]
        objectGroupList
          | comparison == "combinable" = getCompatibleList comparison (fmap (fmap snd) objectListList)
          | comparison == "identity" = L.group $ L.sort (snd <$> concat objectListList)
          | otherwise = errorWithoutStackTrace ("Comparison method " <> comparison <> " unrecognized (combinable/identity)")
        uniqueList = zip indexList (fmap head objectGroupList)
        frequencyList =fmap (((/ numGraphs) . fromIntegral) . length) objectGroupList  -- removed parallel
        fullPairList = zip uniqueList frequencyList
        threshold = (fromIntegral thresholdInt / 100.0) :: Double
    in
    --trace ("There are " <> (show $ length objectListList) <> " to filter: " <> (show uniqueList) <> "\n" <> (show objectGroupList) <> " " <> (show frequencyList))
    (fst <$> filter ((>= threshold). snd) fullPairList, snd <$> fullPairList)

-- |  getThresholdEdges takes a threshold and number of graphs and keeps those unique edges present in the threshold percent or
-- higher.  L.sorted by frequency (low to high)
-- modified from getThresholdNodes due to type change in edges
-- used and number from numleaves so can use BV
getThresholdEdges :: (Show a, Ord a) => Int -> Int -> [a] -> ([a], [Double])
getThresholdEdges thresholdInt numGraphsInput objectList
  | thresholdInt < 0 || thresholdInt > 100 = errorWithoutStackTrace "Threshold must be in range [0,100]"
  | null objectList = error "Empty list of object lists in getThresholdEdges"
  | otherwise =
  let threshold = (fromIntegral thresholdInt / 100.0) :: Double
      numGraphs = fromIntegral numGraphsInput
      objectGroupList = L.group $ L.sort objectList
      uniqueList = fmap head objectGroupList
      frequencyList = fmap (((/ numGraphs) . fromIntegral) . length) objectGroupList -- removed parallel
      fullPairList = zip uniqueList frequencyList
  in
  --trace ("There are " <> (show numGraphsIn) <> " to filter: " <> (show uniqueList) <> "\n" <> (show $ fmap length objectGroupList) <> " " <> (show frequencyList))
  (fst <$> filter ((>= threshold). snd) fullPairList, snd <$> fullPairList)


-- | getPostOrderVerts takes a vertex and traverses postorder to root places all visirted nodes in a set of found
-- vertices. Keeps placing new nodes in recursion list until a root is hit.  If a node is already in found set
-- it is not added to list of nodes to recurse
-- returns set of visited nodes
getPostOrderVerts :: P.Gr BV.BV (BV.BV, BV.BV) -> S.Set G.Node -> [G.Node] -> S.Set G.Node
getPostOrderVerts inGraph foundVertSet inVertexList =
  if null inVertexList then foundVertSet
  else
    let firstVertex = head inVertexList
    in
    if S.member firstVertex foundVertSet then getPostOrderVerts inGraph foundVertSet (tail inVertexList)
    else
      let newFoundSet = S.insert firstVertex foundVertSet
          parentVerts = G.pre inGraph firstVertex
      in
      getPostOrderVerts inGraph newFoundSet (inVertexList <> parentVerts)

-- | verticesByPostorder takes a graph and a leaf set and an initially empty found vertex set
-- as the postorder pass takes place form each leaf, each visited vertex is placed in foundVertSet
-- when roots are hit, it recurses back untill all paths are traced to a root.
-- final final rgaph is created and retuyrned from foundVertSet and input list
-- could have edges unconnected to leaves if consistent edge leading to a subtree with inconsistent configuration
-- so are filtered out by making sure each vertex in an edge is in the vertex list
verticesByPostorder :: P.Gr BV.BV (BV.BV, BV.BV) -> [G.LNode BV.BV] ->  S.Set G.Node -> P.Gr BV.BV (BV.BV, BV.BV)
verticesByPostorder inGraph leafNodes foundVertSet
  | G.isEmpty inGraph = error "Empty graph in verticesByPostorder"
  | null leafNodes =
    let vertexIndexList = S.toList foundVertSet
        vertexLabelList = fmap (fromJust . G.lab inGraph) vertexIndexList
        vertexList = zip vertexIndexList vertexLabelList
        edgeList = fmap (verifyEdge vertexIndexList) (G.labEdges inGraph) -- removed parallel
    in G.mkGraph vertexList (concat edgeList)
      | otherwise =
    let firstLeaf = fst $ head leafNodes
        firstVertices = getPostOrderVerts inGraph foundVertSet [firstLeaf]
    in
    verticesByPostorder inGraph (tail leafNodes) (S.union foundVertSet firstVertices)

-- | verifyEdge takes a vertex index list and an edge and checks to see if
-- the subtyending vertices are in the vertex list nad returns teh edge as asingleton list
-- if yes--else empty list (for mapping purposes)
verifyEdge :: [G.Node] -> G.LEdge (BV.BV, BV.BV) -> [G.LEdge (BV.BV, BV.BV)]
verifyEdge vertIndexList inEdge@(e,u,_)
  | e `notElem` vertIndexList = []
  | u `notElem` vertIndexList = []
  | otherwise = [inEdge]

{-
-- | sortInputArgs takes a list of arguments (Strings) nd retuns a pair of lists
-- of strings that are newick or graphviz dotFile filenames for later parsing
sortInputArgs :: [String] -> [String] -> ([T.Text],[T.Text],[String],[String],[String]) -> ([T.Text],[T.Text],[String],[String],[String])
sortInputArgs inContents inArgs (curFEN, curNewick, curDot, curNewFiles, curFENFILES) =
  if null inArgs then (curFEN, curNewick, curDot, curNewFiles, curFENFILES)
  else
    let firstFileName = head inArgs
        firstContents = filter (not . isSpace) $ head inContents
    in
    if head firstContents == '(' then -- Newick/EnhancedNewick
      sortInputArgs (tail inContents) (tail inArgs) (curFEN, T.pack firstContents : curNewick, curDot, firstFileName : curNewFiles, curFENFILES)
    else if head firstContents == '<' then -- ForestEnhancedNewick
      sortInputArgs (tail inContents) (tail inArgs) (T.pack firstContents : curFEN, curNewick, curDot, curNewFiles, firstFileName : curFENFILES)
    else if (head firstContents == 's') || (head firstContents == 'g') || (head firstContents == 'd') then --Dot
      sortInputArgs (tail inContents) (tail inArgs) (curFEN, curNewick, firstFileName : curDot, curNewFiles, curFENFILES)
    else errorWithoutStackTrace("Input file " <> firstFileName <> " does not appear to be Newick, Enhanced Newick, Forest Enhanced Newick or dot format ")

-- | nodeText2String takes a node with text label and returns a node with String label
nodeText2String :: G.LNode T.Text -> G.LNode String
nodeText2String (inIndex, label) = (inIndex, T.unpack label)

-- | fglTextA2TextString converts the graph types from Text A to Text String
fglTextB2Text :: P.Gr b Double -> P.Gr b T.Text
fglTextB2Text inGraph =
  if G.isEmpty inGraph then G.empty
  else
    let labNodes = G.labNodes inGraph
        labEdges = G.labEdges inGraph
        (eList, uList, labelList) = unzip3 labEdges
        --- newLabels = fmap toShortest labelList
        newLabels = fmap (T.pack . show) labelList
        newEdges = zip3 eList uList newLabels
    in
    G.mkGraph labNodes newEdges
-}

-- | addUrRootAndEdges creates a single root and adds edges to existing roots
-- and unconnected leaves
addUrRootAndEdges :: P.Gr String Double -> P.Gr String Double
addUrRootAndEdges inGraph =
  let origLabVerts = G.labNodes inGraph
      origLabEdges = G.labEdges inGraph
      origRootList = getRoots inGraph (fst <$> origLabVerts)
      unconnectedLeafList = getUnConnectedLeaves inGraph (fst <$> origLabVerts)
  in
  -- all ok--no unconnected vertices
  if (length origRootList == 1) && null unconnectedLeafList then inGraph

  -- add edges to unconencted leaves
  else if length origRootList == 1 then
    let newEdgeList = zip3 (replicate (length unconnectedLeafList) (head origRootList)) unconnectedLeafList (replicate (length unconnectedLeafList) 0.0)
    in
    G.mkGraph origLabVerts (origLabEdges <> newEdgeList)

  -- add UR root, edges to existing roots, and edges to unconnected leaves
  else
    let unRootedVertices = origRootList <> unconnectedLeafList
        numOrigVerts = length origLabVerts
        newRoot = (numOrigVerts, "HTU" <> show numOrigVerts)
        newEdgeList = zip3 (replicate (length unRootedVertices) numOrigVerts) unRootedVertices (replicate (length unRootedVertices) 0.0)
    in
    G.mkGraph (origLabVerts <> [newRoot]) (origLabEdges <> newEdgeList)

-- | changeVertexEdgeLabels keeps or removes vertex and edge labels
changeVertexEdgeLabels :: (Show b) => Bool -> Bool -> P.Gr String b -> P.Gr String String
changeVertexEdgeLabels keepVertexLabel keepEdgeLabel inGraph =
  --trace ("CVL: " <> (show (keepVertexLabel, keepEdgeLabel))) $
  let inLabNodes = G.labNodes inGraph
      degOutList = G.outdeg inGraph <$> G.nodes inGraph
      nodeOutList = zip  degOutList inLabNodes
      leafNodeList = snd <$> filter ((==0).fst) nodeOutList
      nonLeafNodeList = snd <$> filter ((>0).fst) nodeOutList
      newNonLeafNodes = if not keepVertexLabel then 
                          zip (fmap fst nonLeafNodeList) (replicate (length nonLeafNodeList) "")
                        else fmap checkMakeLabel nonLeafNodeList
      inLabEdges = G.labEdges inGraph
      inEdges = fmap G.toEdge inLabEdges
      newEdges = if keepEdgeLabel then fmap showLabel inLabEdges
                 else fmap (`G.toLEdge` "") inEdges
  in
  -- trace ("CVEL " <> (show (keepVertexLabel, keepEdgeLabel )))
  G.mkGraph (leafNodeList <> newNonLeafNodes) newEdges
    where showLabel (e,u,l) = (e,u,show l)
          checkMakeLabel (a,b) = if head b /= 'H' then (a, "HTU" <> show a)
                                 else (a,b)

-- | reconcile is the overall function to drive all methods
reconcile :: (String, String, Int, Bool, Bool, Bool, String, [P.Gr String String]) -> PhyG (String, P.Gr String String)
reconcile (localMethod, compareMethod, threshold, connectComponents, edgeLabel, vertexLabel, outputFormat, inputGraphList) = 
  let  --parallel 
        reAnnotate :: P.Gr String String -> PhyG (P.Gr BV.BV (BV.BV, BV.BV))
        reAnnotate = reAnnotateGraphs

        -- intersectionAction :: G.LNode BV.BV -> [G.LEdge (BV.BV,BV.BV)]
        -- intersectionAction = getIntersectionEdges (fmap snd thresholdNodes) thresholdNodes
  in do
        -- Reformat graphs with appropriate annotations, BV.BVs, etc
        processedGraphs <- getParallelChunkTraverse >>= \pTraverse ->
            reAnnotate `pTraverse` inputGraphList
          -- PU.seqParMap PU.myStrategyRDS reAnnotateGraphs inputGraphList -- `using` PU.myParListChunkRDS

        -- Create lists of reindexed unique nodes and edges, identity by BV.BVs
        -- The drops to not reexamine leaves repeatedly
        -- Assumes leaves are first in list
        let numLeaves = getLeafNumber (head processedGraphs)
        let leafNodes = take numLeaves (G.labNodes $ head processedGraphs)
        let firstNodes = G.labNodes $ head processedGraphs
        let numFirstNodes = length firstNodes
        let unionNodes = L.sort $ leafNodes <> addAndReIndexUniqueNodes numFirstNodes (concatMap (drop numLeaves) (G.labNodes <$> tail processedGraphs)) (drop numLeaves firstNodes)
        -- unionEdges = addAndReIndexEdges "unique" unionNodes (concatMap G.labEdges (tail processedGraphs)) (G.labEdges $ head processedGraphs)

        let totallLeafString = L.foldl' L.union [] (fmap (fmap snd . getLeafListNewick) inputGraphList)
        let totallLeafSet = zip [0..(length totallLeafString - 1)] totallLeafString
        
        -- Create Adams II consensus
        --
        adamsII <- A.makeAdamsII totallLeafSet (fmap PhyP.relabelFGLEdgesDouble inputGraphList)
        -- adamsIIInfo = "There are " <> show (length $ G.nodes adamsII) <> " nodes present in Adams II consensus"
        let adamsII' = changeVertexEdgeLabels vertexLabel edgeLabel adamsII
        let adamsIIOutDotString = T.unpack $ renderDot $ toDot $ GV.graphToDot GV.quickParams adamsII
        let adamsIIOutFENString = PhyP.fglList2ForestEnhancedNewickString [PhyP.stringGraph2TextGraph $ PhyP.relabelFGLEdgesDouble adamsII'] True True

        --
        -- Create thresholdMajority rule Consensus and dot string
        -- vertex-based CUN-> Majority rule ->Strict
        --
        let (thresholdNodes', nodeFreqs) = getThresholdNodes compareMethod threshold numLeaves (fmap (drop numLeaves . G.labNodes) processedGraphs)
        let thresholdNodes = leafNodes <> thresholdNodes'

        intersectionPar <- getParallelChunkMap
        let intersectionAction = getIntersectionEdges (fmap snd thresholdNodes) thresholdNodes
        let thresholdEdgesList = intersectionPar intersectionAction thresholdNodes
          -- PU.seqParMap PU.myStrategyRDS (getIntersectionEdges (fmap snd thresholdNodes) thresholdNodes) thresholdNodes  -- `using` PU.myParListChunkRDS
        let thresholdEdges = L.nub $ concat thresholdEdgesList
        -- numPossibleEdges =  ((length thresholdNodes * length thresholdNodes) - length thresholdNodes) `div` 2
        let thresholdConsensusGraph = G.mkGraph thresholdNodes thresholdEdges -- O(n^3)

        -- thresholdConInfo =  "There are " <> show (length thresholdNodes) <> " nodes present in >= " <> (show threshold <> "%") <> " of input graphs and " <> show numPossibleEdges <> " candidate edges"
        --                  <> " yielding a final graph with " <> show (length (G.labNodes thresholdConsensusGraph)) <> " nodes and " <> show (length (G.labEdges thresholdConsensusGraph)) <> " edges"

        -- add back labels for vertices and "GV.quickParams" for G.Gr String Double or whatever
        let labelledTresholdConsensusGraph' = addGraphLabels thresholdConsensusGraph totallLeafSet
        let labelledTresholdConsensusGraph'' = addEdgeFrequenciesToGraph labelledTresholdConsensusGraph' (length leafNodes) nodeFreqs

        -- Add urRoot and edges to existing roots if there are unconnected components and connnectComponets is True
        let labelledTresholdConsensusGraph = if not connectComponents then labelledTresholdConsensusGraph''
                                         else addUrRootAndEdges labelledTresholdConsensusGraph''
        let gvRelabelledConsensusGraph = GO.renameSimpleGraphNodesString $ LG.reindexGraph $ changeVertexEdgeLabels vertexLabel edgeLabel labelledTresholdConsensusGraph
        let thresholdConsensusOutDotString = T.unpack $ renderDot $ toDot $ GV.graphToDot GV.quickParams gvRelabelledConsensusGraph
        let thresholdConsensusOutFENString = PhyP.fglList2ForestEnhancedNewickString [PhyP.stringGraph2TextGraph labelledTresholdConsensusGraph] edgeLabel True

        --
        -- Create threshold EUN and dot string, orignial EUN is threshold = 0
        --
        let allEdges = addAndReIndexEdges "all" unionNodes (concatMap G.labEdges (tail processedGraphs)) (G.labEdges $ head processedGraphs)
        let (thresholdEUNEdges, edgeFreqs) = getThresholdEdges threshold (length processedGraphs) allEdges
        thresholdEUNGraph' <- makeEUN unionNodes thresholdEUNEdges (G.mkGraph unionNodes thresholdEUNEdges)

        -- Remove unnconnected HTU nodes via postorder pass from leaves
        let thresholdEUNGraph = verticesByPostorder thresholdEUNGraph' leafNodes S.empty
        -- thresholdEUNInfo =  "\nThreshold EUN deleted " <> show (length unionEdges - length (G.labEdges thresholdEUNGraph) ) <> " of " <> show (length unionEdges) <> " total edges"
        --                    <> " for a final graph with " <> show (length (G.labNodes thresholdEUNGraph)) <> " nodes and " <> show (length (G.labEdges thresholdEUNGraph)) <> " edges"

        -- add back labels for vertices and "GV.quickParams" for G.Gr String Double or whatever
        let thresholdLabelledEUNGraph' = addGraphLabels thresholdEUNGraph totallLeafSet
        let thresholdLabelledEUNGraph'' = addEdgeFrequenciesToGraph thresholdLabelledEUNGraph' (length leafNodes) edgeFreqs

        -- Add urRoot and edges to existing roots if there are unconnected components and connnectComponets is True
        let thresholdLabelledEUNGraph = if not connectComponents then thresholdLabelledEUNGraph''
                                        else addUrRootAndEdges thresholdLabelledEUNGraph''

        -- Create EUN Dot String
        let gvRelabelledEUNGraph = GO.renameSimpleGraphNodesString $ LG.reindexGraph $ changeVertexEdgeLabels vertexLabel edgeLabel thresholdLabelledEUNGraph
        let thresholdEUNOutDotString = T.unpack $ renderDot $ toDot $ GV.graphToDot GV.quickParams gvRelabelledEUNGraph -- eunGraph
        let thresholdEUNOutFENString = PhyP.fglList2ForestEnhancedNewickString [PhyP.stringGraph2TextGraph thresholdLabelledEUNGraph] edgeLabel True

  
        -- Create Adams II consensus
        --
        adamsII <- A.makeAdamsII totallLeafSet (fmap PhyP.relabelFGLEdgesDouble inputGraphList)
        -- adamsIIInfo = "There are " <> show (length $ G.nodes adamsII) <> " nodes present in Adams II consensus"
        let adamsII' = changeVertexEdgeLabels vertexLabel False adamsII
        let adamsIIOutDotString = T.unpack $ renderDot $ toDot $ GV.graphToDot GV.quickParams adamsII'
        let adamsIIOutFENString = PhyP.fglList2ForestEnhancedNewickString [PhyP.stringGraph2TextGraph $ PhyP.relabelFGLEdgesDouble adamsII'] False False

        if localMethod == "eun" then
          if outputFormat == "dot" then pure (thresholdEUNOutDotString,  gvRelabelledEUNGraph)
          else if outputFormat == "fenewick" then pure (thresholdEUNOutFENString, gvRelabelledEUNGraph)
          else errorWithoutStackTrace ("Output graph format " <> outputFormat <> " is not implemented")

        else if localMethod == "adams" then
          if outputFormat == "dot" then pure (adamsIIOutDotString,  adamsII')
          else if outputFormat == "fenewick" then pure (adamsIIOutFENString, adamsII')
          else errorWithoutStackTrace ("Output graph format " <> outputFormat <> " is not implemented")

        else if (localMethod == "majority") || (localMethod == "cun") || (localMethod == "strict") then
            if outputFormat == "dot" then pure (thresholdConsensusOutDotString, gvRelabelledConsensusGraph)
            else if outputFormat == "fenewick" then pure (thresholdConsensusOutFENString, gvRelabelledConsensusGraph)
            else errorWithoutStackTrace ("Output graph format " <> outputFormat <> " is not implemented")

        else errorWithoutStackTrace ("Graph combination method " <> localMethod <> " is not implemented")


-- | makeProcessedGraph takes a set of graphs and a leaf set and adds teh missing leafws to teh graphs and reindexes
-- the nodes and edges of the input graphs consistenly
-- String as oposed to Text due tyo reuse of code in Eun.c
makeProcessedGraph :: [LG.LNode T.Text] -> SimpleGraph -> PhyG SimpleGraph
makeProcessedGraph leafTextList inGraph
  | null leafTextList = error "Null leaf list in makeFullLeafSetGraph"
  | LG.isEmpty inGraph = error "Empty graph in makeFullLeafSetGraph"
  | otherwise = let (_, graphleafTextList, _, _) = LG.splitVertexList inGraph
                    leafStringList = fmap nodeToString leafTextList
                    graphLeafStringList = fmap nodeToString graphleafTextList
                in do
                      reIndexedGraph <- reIndexAndAddLeavesEdges leafStringList (graphLeafStringList, inGraph)
                      let textNodes = (nodeToText <$> LG.labNodes reIndexedGraph)
                      let doubleEdges = (edgeToDouble <$> LG.labEdges reIndexedGraph)
                      pure $ LG.mkGraph textNodes doubleEdges
  where
      nodeToString (a, b) = (a, T.unpack b)
      nodeToText (a, b) = (a, T.pack b)
      edgeToDouble (a, b, c) = (a, b, read c :: Double)
