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

module Graphs.GraphOperations ( ladderizeGraph
                       , verifyTimeConsistency
                       , rerootGraph
                       , rerootGraph'
                       , rerootPhylogeneticGraph
                       , rerootPhylogeneticGraph'
                       , generateDisplayTrees
                       , contractOneOneEdges
                       , getNodeType
                       , convertDecoratedToSimpleGraph
                       ) where

import           Debug.Trace
import Data.List
import Data.Maybe
import Data.Bits ((.&.), (.|.))
import qualified Data.Vector as V
import qualified DirectOptimization.DOWrapper as DOW
import Types.Types
import qualified Utilities.LocalGraph as LG
import qualified Data.Text.Lazy as T
import GeneralUtilities
import qualified Utilities.LocalSequence as LS
import qualified GraphFormatUtilities as GFU
import qualified GraphOptimization.Medians as M

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
    if inNum == 0 || outNum == 0 then error ("ResolveNode error: Isolated vertex " ++ (show curNode ) ++ " in graph\n" ++ (LG.prettify inGraph) )

    -- indegree 1 outdegree 1 node to contract
    else if inNum == 1 && outNum == 1 then
      let newEdge = (fst3 $ head inEdgeList, snd3 $ head outEdgeList, 0.0 :: Double)
          newGraph = LG.insEdge newEdge $ LG.delNode curNode $ LG.delLEdges (inEdgeList ++ outEdgeList) inGraph
      in
      newGraph

    -- leaf leaf with too manuy parents
    else if (inNum > 1) && (outNum == 0) then
      let first2Edges = take 2 inEdgeList
          newNode = (numNodes , T.pack $ ("HTU" ++ (show numNodes)))
          newEdge1 = (fst3 $ head first2Edges, numNodes, 0.0 :: Double)
          newEdge2 = (fst3 $ last first2Edges, numNodes, 0.0 :: Double)
          newEdge3 = (numNodes, curNode, 0.0 :: Double)
          newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
      in
      newGraph 

    -- network node with too many parents
    else if (inNum > 2) && (outNum == 1) then
      let first2Edges = take 2 inEdgeList
          newNode = (numNodes , T.pack $ ("HTU" ++ (show numNodes)))
          newEdge1 = (fst3 $ head first2Edges, numNodes, 0.0 :: Double)
          newEdge2 = (fst3 $ last first2Edges, numNodes, 0.0 :: Double)
          newEdge3 = (numNodes, curNode, 0.0 :: Double)
          newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
      in
      newGraph 

    -- tree node or root node with too many children
    else if (inNum < 2 || outNum > 2) then 
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


    else error ("This can't happen in resolveNode in/out edge lists don't need to be resolved " ++ show inOutPair ++ "\n" ++ LG.prettify inGraph)
    --)

-- | rerootGraph' flipped version of rerootGraph 
rerootGraph' :: (Eq b) => LG.Gr a b -> Int -> LG.Gr a b 
rerootGraph' inGraph rerootIndex = rerootGraph rerootIndex inGraph

-- | rerootGraph takes a graph and reroots based on a vertex index (usually leaf outgroup)
--   if input is a forest then only roots the component that contains the vertex wil be rerooted
--   unclear how will effect network edges--will need to verify that does not create cycles
--   multi-rooted components (as opposed to forests) are unaffected with trace warning thrown
--   after checking for existing root and multiroots, should be O(n) where 'n is the length
--   of the path between the old and new root 
rerootGraph :: (Eq b) => Int -> LG.Gr a b -> LG.Gr a b
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
            newEdgeOnOldRoot = if (length originalRootEdges) /= 2 then error ("Number of root out edges /= 1 in rerootGraph")
                               else (snd3 $ head originalRootEdges, snd3 $ last originalRootEdges, thd3 $ head originalRootEdges)

            newRootEdges = [leftChildEdge, rightChildEdge, newEdgeOnOldRoot]
            newGraph = LG.insEdges newRootEdges $ LG.delLEdges (newRootOrigEdge : originalRootEdges) inGraph

            -- get edges that need reversing
            newGraph' = preTraverseAndFlipEdges [leftChildEdge,rightChildEdge] newGraph 

        in
        --trace ("=") 
        --trace ("Deleting " ++ (show (newRootOrigEdge : originalRootEdges)) ++ "\nInserting " ++ (show newRootEdges)) (
        --trace ("In " ++ (GFU.showGraph inGraph) ++ "\nNew " ++  (GFU.showGraph newGraph) ++ "\nNewNew "  ++  (GFU.showGraph newGraph'))
        newGraph'
        --))

-- | rerootPhylogeneticGraph' flipped version of rerootPhylogeneticGraph 
rerootPhylogeneticGraph' :: PhylogeneticGraph -> Int -> PhylogeneticGraph
rerootPhylogeneticGraph' inGraph rerootIndex = rerootPhylogeneticGraph rerootIndex inGraph

-- | rerootGraph takes a pphylogenetic graph and reroots based on a vertex index (usually leaf outgroup)
--   if input is a forest then only roots the component that contains the vertex wil be rerooted
--   unclear how will effect network edges--will need to verify that does not create cycles
--   multi-rooted components (as opposed to forests) are unaffected with trace warning thrown
--   after checking for existing root and multiroots, should be O(n) where 'n is the length
--   of the path between the old and new root 
--   the PhyloGenetic graph is (minimally) reoptimized along the spine of edegs that are redirected
--   in order that the root costr be correct.  The block display forest (components are always trees--for soft-wired
--   graphs only) is also rerooted, and
--   Character foci set to the new root edge
--   NB--Uses calls to rerootGraph since traversals are for different graphs so wouldn't save
--   much time by consolidating--also since labels are all different--can'treuse alot of info
--   from graph to graph.
rerootPhylogeneticGraph :: Int -> PhylogeneticGraph -> PhylogeneticGraph
rerootPhylogeneticGraph rerootIndex inPhyGraph@(inSimple, inCost, inDecGraph, blockDisplayForestVect, inFociVect, charInfoVectVect) = 
  if LG.isEmpty inSimple || LG.isEmpty inDecGraph then error "Empty graph in rerootPhylogeneticGraph"
  else if inCost == 0 then error "Input graph with cost zero--likely non decorated input garph in rerootPhylogeneticGraph"
  else 
    let -- simple graph rerooted Boolean to specify that non-exact characters need NOT be reoptimized if affected
        newSimpleGraph = rerootGraph rerootIndex inSimple

        -- decorated graph Boolean to specify that non-exact characters need to be reoptimized if affected
        newDecGraph = rerootGraph rerootIndex inDecGraph

        -- reoptimize nodes here
        -- nodes on spine from new root to old root that needs to be reoptimized
        (nodesToOptimize, _) = LG.pathToRoot inDecGraph (rerootIndex, fromJust $ LG.lab inDecGraph rerootIndex)
        -- reversed because ll these node edges are reversed so preorder would be in reverse orientation
        reoptimizedNodes = reOptimizeNodes charInfoVectVect newDecGraph (reverse nodesToOptimize)

        newDecGraph' =  LG.insNodes reoptimizedNodes $ LG.delNodes (fmap fst nodesToOptimize) newDecGraph

        -- sum of root costs on Decorated graph
        newGraphCost = sum $ fmap subGraphCost $ fmap snd $ LG.getRoots newDecGraph'

        -- rerooted diplay forests--don't care about costs--I hope (hence Bool False)
        newBlockDisplayForestVect = if V.null blockDisplayForestVect then V.empty
                                    else V.map (rerootGraph rerootIndex) blockDisplayForestVect
        -- the edge the rerooting was switched to (vect vect vect)
        newRootOrigEdge = head $ LG.inn inSimple rerootIndex
        newCharacterFoci = makeCharFociVVV (LG.toEdge newRootOrigEdge) (V.map V.length charInfoVectVect)
        in
        --trace ("=") 
        --trace ("Deleting " ++ (show (newRootOrigEdge : originalRootEdges)) ++ "\nInserting " ++ (show newRootEdges)) (
        --trace ("In " ++ (GFU.showGraph inGraph) ++ "\nNew " ++  (GFU.showGraph newGraph) ++ "\nNewNew "  ++  (GFU.showGraph newGraph'))
        (newSimpleGraph, newGraphCost, newDecGraph', newBlockDisplayForestVect, newCharacterFoci, charInfoVectVect)


-- | reOptimizeNodes takes a decorated graph and a list of nodes and reoptimizes (relabels)
-- them based on children in input graph
-- simple recursive since each node depends on children
-- remove check for debbubg after it works
-- check for out-degree 1, doens't matter for trees however.
reOptimizeNodes :: V.Vector (V.Vector CharInfo) -> DecoratedGraph -> [LG.LNode VertexInfo] -> [LG.LNode VertexInfo]
reOptimizeNodes charInfoVectVect inGraph oldNodeList =
  if null oldNodeList then []
  else 
    let curNode = head oldNodeList
        -- debugging check of unoptimized nodes  if childern of cur node in list of ones yet to be done
        curNodeChildren = V.toList $ children $ snd curNode 
        foundCurChildern = filter (`elem` curNodeChildren) $ fmap fst (tail oldNodeList)
    in
    if (not $ null foundCurChildern) then error ("Current node has childer in optimize list" ++ show oldNodeList)
    -- code form postDecorateTree
    else 
        let nodeChildren = LG.descendants inGraph (fst curNode)  -- should be 1 or 2, not zero since all leaves already in graph
            leftChild = (head nodeChildren)
            rightChild = (last nodeChildren)
            leftChildLabel = fromJust $ LG.lab inGraph leftChild
            rightChildLabel = fromJust $ LG.lab inGraph rightChild
            curnodeLabel = snd curNode
            newVertexData = createVertexDataOverBlocksNonExact  (vertData leftChildLabel) (vertData  rightChildLabel) (vertData curnodeLabel) charInfoVectVect []
        in
        --debug remove when not needed--checking to see if node should not be re optimized
        if (sort nodeChildren) == (sort $ V.toList $ children curnodeLabel) then 
          trace ("Children for vertex unchanged " ++ (show $ fst curNode))
          curNode : reOptimizeNodes charInfoVectVect inGraph (tail oldNodeList)
        else 
           let newCost =  V.sum $ V.map (V.sum) $ V.map (V.map snd) newVertexData
               newVertexLabel = VertexInfo {  index = fst curNode
                                            , bvLabel = (bvLabel leftChildLabel) .&. (bvLabel rightChildLabel)
                                            , parents = V.fromList $ LG.parents inGraph (fst curNode)
                                            , children = V.fromList nodeChildren
                                            , nodeType = nodeType curnodeLabel
                                            , vertName = vertName curnodeLabel
                                            , vertData = V.map (V.map fst) newVertexData
                                            , vertexCost = newCost
                                            , subGraphCost = (subGraphCost leftChildLabel) + (subGraphCost rightChildLabel) + newCost
                                            }   
            in
            (fst curNode, newVertexLabel) : reOptimizeNodes charInfoVectVect inGraph (tail oldNodeList)

-- | createVertexDataOverBlocksNonExact takes data in blocks and block vector of char info and 
-- extracts the quartet for each block and creates new block data for parent node (usually)
-- not checking if vectors are equal in length
-- this only reoptimized the nonexact characters (sequence characters for now, perhpas otehrs later)
-- and takes the existing optimization for exact (Add, NonAdd, Matrix) for the others. 
createVertexDataOverBlocksNonExact :: VertexBlockData -> VertexBlockData -> VertexBlockData -> V.Vector (V.Vector CharInfo) -> [V.Vector (CharacterData, VertexCost)] -> V.Vector (V.Vector (CharacterData, VertexCost))
createVertexDataOverBlocksNonExact leftBlockData rightBlockData existingBlockData blockCharInfoVect curBlockData =
    if V.null leftBlockData then 
        --trace ("Blocks: " ++ (show $ length curBlockData) ++ " Chars  B0: " ++ (show $ V.map snd $ head curBlockData))
        V.fromList $ reverse curBlockData
    else
        let firstBlock = V.zip4 (V.head leftBlockData) (V.head rightBlockData) (V.head existingBlockData) (V.head blockCharInfoVect)
            firstBlockMedian = M.median2NonExact firstBlock
        in
        createVertexDataOverBlocksNonExact (V.tail leftBlockData) (V.tail rightBlockData) (V.tail existingBlockData) (V.tail blockCharInfoVect) (firstBlockMedian : curBlockData)

--M.median2 $ V.zip3 (vertData leftChildLabel) (vertData  rightChildLabel) blockCharInfo

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
      
-- | makeCharFociVVV takes an edge and creates the Vector Vector Vector structure for that edge based
-- on charInfo 
makeCharFociVVV :: LG.Edge -> V.Vector Int -> V.Vector (V.Vector (V.Vector LG.Edge))
makeCharFociVVV inEdge lengthList =
  if V.null lengthList then V.empty
  else 
    let base = V.singleton inEdge
        result = V.map (localReplicate base) lengthList
    in
    result
      where localReplicate a n =  V.replicate n a


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
    let firstGroup@((firstInEdges, firstOutEdges), firstNode) = head inEdgeNodeList
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
    let (edgeList, index) = head netEdgeIndexPairList
        --edgeToKeep = edgeList !! index
        edgesToDelete = (take index edgeList) ++ (drop (index + 1) edgeList)
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
    {-

      let (_, _, _, netVertList) = LG.splitVertexList inGraph
          beforeLists = [] -- fmap fst $ fmap (nodesAndEdgesBefore inGraph ([],[])) (fmap (:[]) $ (fmap fst netVertList))
          afterLists = fmap fst $ fmap (nodesAndEdgesAfter inGraph ([],[])) (fmap (:[]) $ (fmap fst netVertList))
          allLists = beforeLists ++ afterLists
          areCompatible = checkCompatible allLists 
      in
      if areCompatible then inGraph
      else errorWithoutStackTrace ("Graph violates time consistency" ++ GFU.showGraph inGraph)
    -}

-- | checkCompaptible takes a list of a list of nodes and checks the node lists
-- to see if compatible-- ie A intersect B = A|B|[]
checkCompatible :: [[LG.Node]] -> Bool
checkCompatible inNodeListList = 
  if null inNodeListList then True
  else 
    let firstList = head inNodeListList
        compatibleList = fmap (interCheck firstList) (tail inNodeListList)
    in
    if (not $ foldl' (&&) True compatibleList) then False
    else checkCompatible (tail inNodeListList)

  where 
    interCheck a b = 
      let c = a `intersect` b 
      in
      if null c then True
      else if (c == a) || (c == b) then True
      else False

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

        decEdgeList = LG.labEdges inDec
        sourceList = fmap fst3 decEdgeList
        sinkList = fmap snd3 decEdgeList
        newEdgeLables = fmap minLength $ fmap thd3 decEdgeList
        simpleEdgeList = zip3 sourceList sinkList newEdgeLables
    in
    LG.mkGraph simpleNodes simpleEdgeList
