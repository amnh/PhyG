{- |
Module      :  PostOrderSoftWiredFunctions.hs
Description :  Module specifying post-order softwiired graph functions
Copyright   :  (c) 2022 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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

{-
ToDo:
   Add parallel optimization overblocks and characters?
-}


module GraphOptimization.PostOrderSoftWiredFunctions  ( updateAndFinalizePostOrderSoftWired
                                                      , postOrderSoftWiredTraversal
                                                      , postDecorateSoftWired'
                                                      , postDecorateSoftWired
                                                      , assignPostOrderToDisplayTree
                                                      , softWiredPostOrderTraceBack
                                                      , makeLeafGraphSoftWired
                                                      , getOutDegree1VertexAndGraph
                                                      , getOutDegree1VertexSoftWired
                                                      , getOutDegree2VertexSoftWired
                                                      , createBlockResolutions
                                                      , extractDisplayTrees
                                                      , getOutDegree1VertexAndGraph
                                                      , createBlockResolutions
                                                      , makeCharacterGraph
                                                      , getAllResolutionList
                                                      ) where

import           Data.Bits
import qualified Data.BitVector.LittleEndian as BV
import qualified Data.List                   as L
import           Data.Maybe
import qualified Data.Text.Lazy              as T
import qualified Data.Vector                 as V
import           GeneralUtilities
import qualified GraphOptimization.Medians   as M
import qualified Graphs.GraphOperations      as GO
import           Types.Types
import qualified Utilities.LocalGraph        as LG
import qualified Utilities.Utilities         as U
import           Control.Parallel.Strategies
import qualified ParallelUtilities           as PU
-- import Debug.Debug
import           Debug.Trace


-- | updateAndFinalizePostOrderSoftWired performs the pre-order traceback on the resolutions to create the correct vertex states,
-- ports the post order assignments to the canonical tree and display trees, and creates the character trees from the block trees
updateAndFinalizePostOrderSoftWired :: Maybe Int -> Int -> PhylogeneticGraph -> PhylogeneticGraph
updateAndFinalizePostOrderSoftWired startVertex rootIndex inGraph =
    if LG.isEmpty $ thd6 inGraph then inGraph
    else
        let -- is this repetetive--aftet post order?  Done already?
            -- outgroupRootLabel =  fromJust $ LG.lab (thd6 inGraph) rootIndex
            -- (displayGraphVL, lDisplayCost) = PO.extractDisplayTrees startVertex True (vertexResolutionData outgroupRootLabel)
            displayGraphVL = fth6 inGraph
            lDisplayCost = snd6 inGraph

            -- traceback on resolutions
            newGraph = softWiredPostOrderTraceBack rootIndex (thd6 inGraph)

            -- propagate updated post-order assignments to display trees, which updates dispay tree and cost ealier
            displayGraphVL' = V.zipWith (assignPostOrderToDisplayTree (fmap (vertData . snd) (LG.labNodes newGraph) )) displayGraphVL (V.fromList [0..(V.length displayGraphVL - 1)])

            -- create new, fully  updated post-order graph
            finalPreOrderGraph = (fst6 inGraph, lDisplayCost, newGraph, displayGraphVL', divideDecoratedGraphByBlockAndCharacterSoftWired displayGraphVL', six6 inGraph)
        in
        -- trace ("UFPOSW: " ++ (show $ fmap length displayGraphVL) ++ " " ++ (show $ fmap length displayGraphVL') ++ " " ++ (show $ fmap V.length $ fft6 finalPreOrderGraph))
        finalPreOrderGraph

-- | divideDecoratedGraphByBlockAndCharacterSoftWired takes a Vector of a list of DecoratedGraph
-- continaing a list of decorated trees that are the display trees for that block
-- with (potentially) multiple blocks
-- and (potentially) multiple character per block and creates a Vector of Vector of Decorated Graphs
-- over blocks and characters with the block diplay graph, but only a single block and character for each graph
-- this to be used to create the "best" cost over alternate graph traversals
-- vertexCost and subGraphCost will be taken from characterData localcost/localcostVect and globalCost
-- for this assignment purpose for pre-order a single (head) member of list is used to create the
-- character graphs
divideDecoratedGraphByBlockAndCharacterSoftWired :: V.Vector [DecoratedGraph] -> V.Vector (V.Vector DecoratedGraph)
divideDecoratedGraphByBlockAndCharacterSoftWired inGraphVL =
  if V.null inGraphVL then mempty
  else
    let blockGraphList = fmap head inGraphVL
        characterGraphList = fmap makeCharacterGraph blockGraphList
    in
    characterGraphList

-- | postOrderSoftWiredTraversal performs postorder traversal on Soft-wired graph
-- staticIA is ignored--but kept for functional polymorphism
-- ur-root = ntaxa is an invariant
postOrderSoftWiredTraversal :: GlobalSettings -> ProcessedData -> DecoratedGraph -> Bool -> Maybe Int -> SimpleGraph -> PhylogeneticGraph
postOrderSoftWiredTraversal inGS inData@(_, _, blockDataVect) leafGraph _ startVertex inSimpleGraph =
    if LG.isEmpty inSimpleGraph then emptyPhylogeneticGraph
    else
         -- Assumes root is Number of Leaves
        let rootIndex = if startVertex == Nothing then V.length $ fst3 inData
                        else fromJust startVertex
            blockCharInfo = V.map thd3 blockDataVect
            newSoftWired = postDecorateSoftWired inGS inSimpleGraph leafGraph blockCharInfo rootIndex rootIndex
        in
        --trace ("It Begins at " ++ show rootIndex) (
        -- trace ("POSWT:\n" ++ (LG.prettify inSimpleGraph) ++ "\nVertices:\n" ++ (show $ LG.labNodes $ thd6 newSoftWired)) (
        if (startVertex == Nothing) && (not $ LG.isRoot inSimpleGraph rootIndex) then
            let localRootList = fst <$> LG.getRoots inSimpleGraph
                localRootEdges = concatMap (LG.out inSimpleGraph) localRootList
                currentRootEdges = LG.out inSimpleGraph rootIndex
            in
            error ("Index "  ++ show rootIndex ++ " with edges " ++ show currentRootEdges ++ " not root in graph:" ++ show localRootList ++ " edges:" ++ show localRootEdges ++ "\n" ++ LG.prettify inSimpleGraph)
        else
            -- trace ("POSW:" ++ (show $ fmap V.length $ fft6 newSoftWired))
            newSoftWired
        -- )

-- | postDecorateSoftWired' wrapper for postDecorateSoftWired with args in differnt order for mapping
postDecorateSoftWired' :: GlobalSettings -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.Node -> LG.Node -> SimpleGraph -> PhylogeneticGraph
postDecorateSoftWired' inGS curDecGraph blockCharInfo rootIndex curNode simpleGraph = postDecorateSoftWired inGS simpleGraph curDecGraph blockCharInfo rootIndex curNode

-- | postDecorateSoftWired begins at start index (usually root, but could be a subtree) and moves preorder till children are labelled
-- and then recurses to root postorder labelling vertices and edges as it goes
-- this for a single root
postDecorateSoftWired :: GlobalSettings -> SimpleGraph -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.Node -> LG.Node -> PhylogeneticGraph
postDecorateSoftWired inGS simpleGraph curDecGraph blockCharInfo rootIndex curNode =
    -- if node in current decortated graph then nothing to do and return
    if LG.gelem curNode curDecGraph then
        let nodeLabel = LG.lab curDecGraph curNode
        in
        if isNothing nodeLabel then error ("Null label for node " ++ show curNode)
        else (simpleGraph, subGraphCost (fromJust nodeLabel), curDecGraph, mempty, mempty, blockCharInfo)

    else
        -- get postodre assignmens of children
        -- checks for single child of node
        -- result is single graph afer left and right child traversals
        -- trace ("PDSW making node " ++ show curNode ++ " in\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph curDecGraph)) (
        let nodeChildren = LG.descendants simpleGraph curNode  -- should be 1 or 2, not zero since all leaves already in graph
            leftChild = head nodeChildren
            rightChild = last nodeChildren
            leftChildTree = postDecorateSoftWired inGS simpleGraph curDecGraph blockCharInfo rootIndex leftChild
            rightLeftChildTree = if length nodeChildren == 2 then postDecorateSoftWired inGS simpleGraph (thd6 leftChildTree) blockCharInfo rootIndex rightChild
                                 else leftChildTree
        in
        -- Checks on children
        if length nodeChildren > 2 then error ("Graph not dichotomous in postDecorateSoftWired node " ++ show curNode ++ "\n" ++ LG.prettify simpleGraph)
        else if null nodeChildren then error ("Leaf not in graph in postDecorateSoftWired node " ++ show curNode ++ "\n" ++ LG.prettify simpleGraph)

        else
            -- make node from child block resolutions
            -- child resolutin made ealeri in post roder pass
            let newSubTree = thd6 rightLeftChildTree
            in

            -- single child of node (can certinly happen with soft-wired networks
            if length nodeChildren == 1 then
                -- trace ("Outdegree 1: " ++ (show curNode) ++ " " ++ (show $ GO.getNodeType simpleGraph curNode) ++ " Child: " ++ (show nodeChildren)) (
                let (newGraph, _, _, _, _) = getOutDegree1VertexAndGraph curNode (fromJust $ LG.lab newSubTree leftChild) simpleGraph nodeChildren newSubTree
                in
                (simpleGraph, 0, newGraph, mempty, mempty, blockCharInfo)


            -- 2 children
            else
                -- trace ("Outdegree 2: " ++ (show curNode) ++ " " ++ (show $ GO.getNodeType simpleGraph curNode) ++ " Children: " ++ (show nodeChildren)) (
                -- need to create new resolutions and add to existing sets
                let -- this ensures that left/right choices are based on leaf BV for consistency and label invariance
                    -- larger bitvector is Right, smaller or equal Left
                    ((leftChild', leftChildLabel), (rightChild', rightChildLabel)) = U.leftRightChildLabelBVNode ((leftChild, fromJust $ LG.lab newSubTree leftChild), (rightChild, fromJust $ LG.lab newSubTree rightChild))

                    -- create resolution caches for blocks
                    leftChildNodeType  = nodeType leftChildLabel
                    rightChildNodeType = nodeType rightChildLabel
                    resolutionBlockVL = V.zipWith3 (createBlockResolutions (compressResolutions inGS) curNode leftChild' rightChild' leftChildNodeType rightChildNodeType (GO.getNodeType simpleGraph curNode)) (vertexResolutionData leftChildLabel) (vertexResolutionData rightChildLabel) blockCharInfo

                    -- create canonical Decorated Graph vertex
                    -- 0 cost becasue can't know cosrt until hit root and get best valid resolutions
                    newVertexLabel = VertexInfo { index = curNode
                                                , bvLabel = bvLabel leftChildLabel .|. bvLabel rightChildLabel
                                                , parents = V.fromList $ LG.parents simpleGraph curNode
                                                , children = V.fromList nodeChildren
                                                , nodeType = GO.getNodeType simpleGraph curNode
                                                , vertName = T.pack $ "HTU" ++ show curNode
                                                , vertData = mempty --empty because of resolution data
                                                , vertexResolutionData = resolutionBlockVL
                                                , vertexCost = 0.0 --newCost
                                                , subGraphCost = 0.0 -- (subGraphCost leftChildLabel) + (subGraphCost rightChildLabel) + newCost
                                                }

                    leftEdgeType
                      | leftChildNodeType == NetworkNode = NetworkEdge
                      | leftChildNodeType == LeafNode = PendantEdge
                      | otherwise = TreeEdge
                    rightEdgeType
                      | rightChildNodeType == NetworkNode = NetworkEdge
                      | rightChildNodeType == LeafNode = PendantEdge
                      | otherwise = TreeEdge

                    edgeLable = EdgeInfo { minLength = 0.0
                                         , maxLength = 0.0
                                         , midRangeLength = 0.0
                                         , edgeType = TreeEdge
                                         }

                    leftEdge =  (curNode, leftChild', edgeLable {edgeType = leftEdgeType})
                    rightEdge = (curNode, rightChild', edgeLable {edgeType = rightEdgeType})
                    newGraph =  LG.insEdges [leftEdge, rightEdge] $ LG.insNode (curNode, newVertexLabel) newSubTree

                    (displayGraphVL, lDisplayCost) = if curNode == rootIndex then extractDisplayTrees (Just curNode) True resolutionBlockVL
                                                     else (mempty, 0.0)

                in
                (simpleGraph, lDisplayCost, newGraph, displayGraphVL, mempty, blockCharInfo)

-- | assignPostOrderToDisplayTree takes the post-order (preliminary) block data from canonical Decorated graph
-- to Block display graphs (through list of more than one for a block)
-- input is canonical Decorated Graph and a pair containing the display tree and its Block index
-- this could be integrated into preorder traceback to remove the extra pass
-- do this here is a bit simpler
assignPostOrderToDisplayTree :: [VertexBlockData] -> [DecoratedGraph] -> Int -> [DecoratedGraph]
assignPostOrderToDisplayTree canonicalVertexBlockData displayTreeList displayTreeIndex =
    if null canonicalVertexBlockData then []
    else
        let
            updatedDispayTreeList = fmap (assignVertexBlockData canonicalVertexBlockData displayTreeIndex) displayTreeList
        in
        -- trace ("Index: " ++ (show displayTreeIndex) ++ " Blocks : " ++ (show $ V.length $ head canonicalVertexBlockData) ++ " Nodes: "  ++ (show $ length canonicalVertexBlockData))
        updatedDispayTreeList

-- | assignVertexBlockData assigns the block data of input index and assigns to all nodes of a tree
assignVertexBlockData :: [VertexBlockData] -> Int -> DecoratedGraph -> DecoratedGraph
assignVertexBlockData nodeDataList blockIndex inGraph =
    if null nodeDataList || LG.isEmpty inGraph then LG.empty
    else
        -- trace ("In Index: " ++ (show blockIndex) ++ " BL: " ++ (show $ V.length $ head nodeDataList) ++ " lengths " ++ (show $ fmap V.length nodeDataList)) (
        let blockIndexDataList = fmap (V.! blockIndex) nodeDataList
            (displayNodeIndexList, displayNodeLabelList) = L.unzip $ LG.labNodes inGraph
            updatedNodeLabelList = zipWith updateVertData displayNodeLabelList blockIndexDataList
        in
        LG.mkGraph (zip displayNodeIndexList updatedNodeLabelList) (LG.labEdges inGraph)
        -- )
        where updateVertData inVertexInfo newVertData = inVertexInfo {vertData = V.singleton newVertData}

-- | softWiredPostOrderTraceBack takes resolution data and assigns correct resolution median to preliminary
-- data ssignments.  Proceeds via typical pre-order pass over display tree
softWiredPostOrderTraceBack :: Int -> DecoratedGraph -> DecoratedGraph
softWiredPostOrderTraceBack rootIndex inGraph  =
    if LG.isEmpty inGraph then LG.empty
    else
            -- get edges to remake graph after nodes are updated with preliminary states
        let rootLabel = fromJust $ LG.lab inGraph rootIndex
            inEdgeList = LG.labEdges inGraph

            -- root first--choose best resolutions--returns a list and takes head of that list of equal cost/traceback preliminary
            -- assignments.  Later could look at multiple.
            (rootVertData, rootSubGraphCost, rootResolutionCost, leftRightIndexVect) = getResolutionDataAndIndices rootLabel (V.singleton (Just (-1)))
            newRootLabel = rootLabel { vertData = rootVertData
                                     , vertexCost = rootResolutionCost
                                     , subGraphCost = rootSubGraphCost
                                     }
            newRootNode = (rootIndex, newRootLabel)

            -- get child/children of root
            rootChildren = LG.labDescendants inGraph newRootNode

            -- left / right to match post-order
            rootChildrenBV = fmap (bvLabel . snd) rootChildren
            rootChildrenIsLeft
              | length rootChildrenBV == 1 = [True]
              | head rootChildrenBV > (rootChildrenBV !! 1) = [False, True]
              | otherwise = [True, False]

            rootChildrenTuples = zip3 rootChildren (replicate (length rootChildren) leftRightIndexVect) rootChildrenIsLeft

            -- recurse to children with resolution index from parent
            softWiredUpdatedNodes = softWiredPrelimTraceback inGraph rootChildrenTuples [newRootNode]

        in
        LG.mkGraph softWiredUpdatedNodes inEdgeList

-- | softWiredPrelimTraceback takes a list of nodes to update (with left right index info) based
-- on resolution data, recurses with left right indices pre-order to leaves, keeping list of updated nodes
softWiredPrelimTraceback :: DecoratedGraph
                        -> [(LG.LNode VertexInfo, V.Vector (Maybe Int, Maybe Int), Bool)]
                        -> [LG.LNode VertexInfo]
                        -> [LG.LNode VertexInfo]
softWiredPrelimTraceback inGraph nodesToUpdate updatedNodes =
    if null nodesToUpdate then updatedNodes
    else
        let (firstNode, firstLeftRight, isLeft) = head nodesToUpdate

            -- ensure consistent left/right from post-order
            resolutionIndexVect = if isLeft then fmap fst firstLeftRight
                                  else fmap snd firstLeftRight

            -- get resolution info
            (newBlockData, newSubGraphCost, newVertexCost, childLeftRightIndexVect) = getResolutionDataAndIndices (snd firstNode) resolutionIndexVect

            -- make new node
            newNodeLabel = (snd firstNode) { vertData = newBlockData
                                           , vertexCost = newVertexCost
                                           , subGraphCost = newSubGraphCost
                                           }

            newFirstNode = (fst firstNode, newNodeLabel)
        in

        -- not really necessary, but avoids the graph query operations
        if nodeType (snd firstNode) == LeafNode then
            softWiredPrelimTraceback inGraph (tail nodesToUpdate) (newFirstNode : updatedNodes)

        -- checks if network node (and its children) has been visited already
        else if (nodeType (snd firstNode) == NetworkNode) && isJust (L.find ((== fst firstNode). fst) updatedNodes) then
            softWiredPrelimTraceback inGraph (tail nodesToUpdate) updatedNodes

        else
            let -- get children
                firstChildren = LG.labDescendants inGraph firstNode
                firstChildrenBV = fmap (bvLabel . snd) firstChildren
                firstChildrenIsLeft
                  | length firstChildrenBV == 1 = [True] -- seems dumb but may assume lenght == 2 later
                  | head firstChildrenBV > (firstChildrenBV !! 1) = [False, True]
                  | otherwise = [True, False]

                childrenTuples = zip3 firstChildren (replicate (length firstChildren) childLeftRightIndexVect) firstChildrenIsLeft


        in
        {-
        if V.head resolutionIndexVect == Nothing then error ("'Nothing' index in softWiredPrelimTraceback " ++ (show resolutionIndexVect)
            ++ " Node " ++ (show $ fst firstNode) ++ " isLeft?: " ++ (show isLeft)  ++ " " ++ (show firstLeftRight))
        else
        -}
            softWiredPrelimTraceback inGraph (childrenTuples ++ tail nodesToUpdate) (newFirstNode : updatedNodes)


-- | getResolutionDataAndIndices takes a vertex label (VertexInfo) and returns the resolution data corresponding to
-- the index taken from its child resolution (data, subgraph cost, local resolutoin cost, left/right pairs).
-- Index = (-1) denotes that it is a root label and in that case
-- the best (lowest cost) resolutions are returned
getResolutionDataAndIndices :: VertexInfo -> V.Vector (Maybe Int) -> (VertexBlockData, VertexCost, VertexCost, V.Vector (Maybe Int, Maybe Int))
getResolutionDataAndIndices nodeLabel parentResolutionIndexVect =

    -- should not happen
    --mtrace ("NL " ++ (show $ index nodeLabel) ++ " PRIL: " ++ " length " ++ (show $ V.length parentResolutionIndexVect) ++ " " ++ show parentResolutionIndexVect) (
    if nodeType nodeLabel == LeafNode then
        let leafVertData = fmap (displayData . V.head) (vertexResolutionData nodeLabel)
        in
        (leafVertData, 0, 0, V.singleton (Just 0, Just 0))

    -- root node--take lowest cost
    else if V.head parentResolutionIndexVect == Just (-1) then
        let rootBlockResolutionPair = getBestBlockResolution <$> vertexResolutionData nodeLabel
            (charDataVV, subGraphCostV, resCostV, leftRightIndexVect) = V.unzip4 rootBlockResolutionPair
        in
        (charDataVV, V.sum subGraphCostV, V.sum resCostV, leftRightIndexVect)

    -- non-root node--return the index resolution information
    else
        -- trace ("GRD Length:" ++ (show $ V.length $ vertexResolutionData nodeLabel) ++ " " ++ (show $ fmap fromJust parentResolutionIndexVect) ++ "\n" ++ (show $ vertexResolutionData nodeLabel)) (
        let parentIndexVect = fmap fromJust parentResolutionIndexVect

            -- get resolution data from node label
            resolutionData = vertexResolutionData nodeLabel

            -- get the correct (via index) resolution data for each block
            -- complex for network node since keeps left right sort of array, but only first element maters--this hack keepsm thingfs ok for
            -- tree-like traceback assignment
            resolutionsByBlockV = if nodeType nodeLabel == NetworkNode then
                                        V.zipWith (V.!) resolutionData (V.replicate (V.length parentIndexVect) (V.head parentIndexVect))
                                  else V.zipWith (V.!) resolutionData parentIndexVect

            -- get other resolution info
            charDataVV = fmap displayData resolutionsByBlockV
            lSubGraphCost = V.sum $ fmap displayCost resolutionsByBlockV
            localResolutionCost = V.sum $ fmap resolutionCost resolutionsByBlockV

            -- only takes first left right pair--although others in there
            -- uses first for preliminary asiignment--but could be more
            -- if euqla cost display trees, may have multiple possible preliminary states
            leftRightIndexVect = fmap (head . childResolutions) resolutionsByBlockV
        in
        {-
        if V.null resolutionData then
           error ("Null resolution data in getResolutionDataAndIndices at node with label:" ++ (show nodeLabel))
        else
        -}
        -- trace ("GRDI: " ++ (show $ nodeType nodeLabel) ++ " " ++ (show $ fmap V.length resolutionData) ++ " " ++ (show parentIndexVect))
        (charDataVV, lSubGraphCost, localResolutionCost, leftRightIndexVect)
        -- )

-- | getBestBlockResolution takes vertexResolutionData and returns the best (lowest cost) resolution and associated data
-- for a single block of characters.  A single left right index pair is returns for child resolution sources. Could be multiple
-- from resolution compression (equal leaf set, equal median) and avaialble for potenitally implementd in future, multiple
-- preliminary assignments.  Only valid for root node with all leaves in graph
getBestBlockResolution :: ResolutionBlockData -> (V.Vector CharacterData, VertexCost, VertexCost, (Maybe Int, Maybe Int))
getBestBlockResolution inResBlockData =
    if V.null inResBlockData then (mempty, 0.0, 0.0, (Nothing, Nothing))
    else
        let -- makes sure all leaves in resolution
            -- displayPopList = fmap (complement . displayBVLabel) inResBlockData

            -- this for non full graph root--uses highest number of bits on--make sure all taxa in
            -- should be fine for a subgraph that has a single edge a base--may not be correct
            -- for a sunbgaph that has conections outside of its graph root.
            displayPopList' = fmap (popCount . displayBVLabel) inResBlockData
            maxPop = maximum displayPopList'



            -- subgraph cost
            displayCostList = fmap displayCost inResBlockData

            -- resolution local cost
            resolutionCostList = fmap resolutionCost inResBlockData

            -- takes only first resolution index pair
            childResolutionList = fmap (head . childResolutions) inResBlockData

            -- resolution medians
            displayDataList = fmap displayData inResBlockData

            {-
            -- take only those will all leaves in, then minimum cost
            quintVect = V.zip5 displayPopList displayCostList resolutionCostList childResolutionList displayDataList
            validVect = V.filter (BV.isZeroVector . fst5) quintVect
            validMinCost = V.minimum $ fmap snd5 validVect
            -}

            -- these for "best" this will largest leaf set
            quintVect' =  V.zip5 displayPopList' displayCostList resolutionCostList childResolutionList displayDataList
            validVect' = V.filter ((== maxPop) . fst5) quintVect'
            validMinCost' = V.minimum $ fmap snd5 validVect'

            -- ONly takes first of potentially multiple soliutions to begin traceback
            (_, displayCostV, resCostV, childIndexPairV, displayMedianV) = V.unzip5 $ V.filter  ((== validMinCost') . snd5) validVect'
        in
        if null validVect' then error "Null valid quad in getBestBlockResolution--perhaps not root node or forest component"
        else (V.head displayMedianV, V.head displayCostV, V.head resCostV, V.head childIndexPairV)

-- | makeLeafGraphSoftWired takes input data and creates a 'graph' of leaves with Vertex information
-- but with zero edges.  This 'graph' can be reused as a starting structure for graph construction
-- to avoid remaking of leaf vertices
-- includes leave resolution data
makeLeafGraphSoftWired :: ProcessedData -> DecoratedGraph
makeLeafGraphSoftWired (nameVect, bvNameVect, blocDataVect) =
    if V.null nameVect then error "Empty ProcessedData in makeLeafGraph"
    else
        let leafVertexList = V.toList $ V.map (makeLeafVertexSoftWired nameVect bvNameVect blocDataVect) (V.fromList [0.. V.length nameVect - 1])
        in
        LG.mkGraph leafVertexList []

-- | makeLeafVertexSoftWired makes a single unconnected vertex for a leaf in a Soft-wired graph
makeLeafVertexSoftWired :: V.Vector NameText -> V.Vector NameBV -> V.Vector BlockData -> Int -> LG.LNode VertexInfo
makeLeafVertexSoftWired nameVect bvNameVect inData localIndex =
    --trace ("Making leaf " ++ (show localIndex) ++ " Data " ++ (show $ length inData) ++ " " ++ (show $ fmap length $ fmap snd3 inData)) (
    let centralData = V.map snd3 inData
        thisData = V.map (V.! localIndex) centralData
        thisBVLabel = bvNameVect V.! localIndex
        thisResolutionData = makeLeafResolutionBlockData thisBVLabel ([(localIndex, minimalVertex)],[]) thisData
        minimalVertex = VertexInfo  { index = localIndex
                                    , bvLabel = thisBVLabel
                                    , parents = V.empty
                                    , children = V.empty
                                    , nodeType = LeafNode
                                    , vertName = nameVect V.! localIndex
                                    , vertData = mempty
                                    , vertexResolutionData = mempty
                                    , vertexCost = 0.0
                                    , subGraphCost = 0.0
                                    }
        newVertex = VertexInfo  { index = localIndex
                                , bvLabel = thisBVLabel
                                , parents = V.empty
                                , children = V.empty
                                , nodeType = LeafNode
                                , vertName =  nameVect V.! localIndex
                                , vertData = mempty
                                , vertexResolutionData = thisResolutionData
                                , vertexCost = 0.0
                                , subGraphCost = 0.0
                                }
        in
        --trace ("RD" ++ show $ thisResolutionData)
        (localIndex, newVertex)
        -- )


-- | makeLeafResolutionBlockData creates leaf resolution data from leav BVLabel, leave node, and data.
-- The return type is a vertor over character blocks each containing a list of potential resolutions (display trees) for that block
-- the resolutoins include subtree (display) for that resolution the bv l;abel for the node given that resolution, and the character data
-- in the block (Vector CharacterData) also given that resolution
-- thiis is repeaatted for each bloick in VertexBlockData
makeLeafResolutionBlockData :: NameBV -> ([LG.LNode VertexInfo], [LG.LEdge EdgeInfo]) -> VertexBlockData -> V.Vector ResolutionBlockData
makeLeafResolutionBlockData inBV inSubGraph inVertData =
    let defaultResolutionData = ResolutionData  { displaySubGraph = inSubGraph
                                                , displayBVLabel = inBV
                                                , displayData = mempty
                                                , childResolutions = [(Just 0, Just 0)]
                                                , resolutionCost = 0.0
                                                , displayCost = 0.0
                                                }

        blockIndexList = [0..(V.length inVertData - 1)]
        blockDataList = fmap (inVertData V.!) blockIndexList
        resolutionDataList = modifyDisplayData defaultResolutionData blockDataList []
        resolutionData =   V.fromList $ fmap (V.fromList . (:[])) resolutionDataList
    in
    resolutionData

-- | modifyDisplayData modifies displatData filed in ResolutionData
-- stas list doesn't change number of V.fromList calls
modifyDisplayData :: ResolutionData -> [V.Vector CharacterData] -> [ResolutionData] -> [ResolutionData]
modifyDisplayData resolutionTemplate characterDataVList curResolutionList =
    if null characterDataVList then reverse curResolutionList
    else
        let curBlockData = head characterDataVList
        in
        modifyDisplayData resolutionTemplate (tail characterDataVList) ((resolutionTemplate {displayData = curBlockData}) : curResolutionList)

-- | getOutDegree1VertexAndGraph makes parent node fomr single child for soft-wired resolutions
getOutDegree1VertexAndGraph :: (Show a, Show b)
                            => LG.Node
                            -> VertexInfo
                            -> LG.Gr a b
                            -> [LG.Node]
                            -> DecoratedGraph
                            -> (DecoratedGraph, Bool, VertexInfo, VertexCost, V.Vector [DecoratedGraph])
getOutDegree1VertexAndGraph curNode childLabel simpleGraph nodeChildren subTree =

    -- trace ("In out=1: " ++ (show curNode)) (
    let childResolutionData = vertexResolutionData childLabel

        curNodeResolutionData = addNodeAndEdgeToResolutionData newDisplayNode newLEdge childResolutionData

        newEdgeLabel = EdgeInfo { minLength = 0.0
                                 , maxLength = 0.0
                                 , midRangeLength = 0.0
                                 , edgeType = TreeEdge
                                 }
        newMinVertex = VertexInfo  { index = curNode
                                    , bvLabel = bvLabel childLabel
                                    , parents = V.fromList $ LG.parents simpleGraph curNode
                                    , children = V.fromList nodeChildren
                                    , nodeType = GO.getNodeType simpleGraph curNode
                                    , vertName = T.pack $ "HTU" ++ show curNode
                                    , vertData = mempty
                                    , vertexResolutionData = mempty
                                    , vertexCost = 0.0
                                    , subGraphCost = 0.0
                                    }

        newVertex  = VertexInfo { index = curNode
                                , bvLabel = bvLabel childLabel
                                , parents = V.fromList $ LG.parents simpleGraph curNode
                                , children = V.fromList nodeChildren
                                , nodeType = GO.getNodeType simpleGraph curNode
                                , vertName = T.pack $ "HTU" ++ show curNode
                                , vertData = mempty
                                , vertexResolutionData = curNodeResolutionData
                                , vertexCost = 0.0
                                , subGraphCost = subGraphCost childLabel
                                }

        newLEdge = (curNode, index childLabel, newEdgeLabel)
        newLNode = (curNode, newVertex)
        newDisplayNode = (curNode, newMinVertex)
        newGraph =  LG.insEdge newLEdge $ LG.insNode newLNode subTree

        (displayGraphVL, lDisplayCost) = if nodeType newVertex == RootNode then extractDisplayTrees Nothing True (vertexResolutionData childLabel)
                                        else (mempty, 0.0)


    in
    --trace ("NV1: " ++ show newVertex)
    (newGraph, nodeType newVertex == RootNode, newVertex, lDisplayCost, displayGraphVL)
    -- (newGraph, False, newVertex, 0.0, mempty)
    -- )

-- | getOutDegree1VertexSoftWired returns new vertex only from single child for soft-wired resolutions
getOutDegree1VertexSoftWired :: (Show a, Show b)
                    => LG.Node
                    -> VertexInfo
                    -> LG.Gr a b
                    -> [LG.Node]
                    -> VertexInfo
getOutDegree1VertexSoftWired curNode childLabel simpleGraph nodeChildren =

    -- trace ("In out=1: " ++ (show curNode)) (
    let childResolutionData = vertexResolutionData childLabel


        newEdgeLabel = EdgeInfo { minLength = 0.0
                                 , maxLength = 0.0
                                 , midRangeLength = 0.0
                                 , edgeType = TreeEdge
                                 }
        newMinVertex = VertexInfo  { index = curNode
                                    , bvLabel = bvLabel childLabel
                                    , parents = V.fromList $ LG.parents simpleGraph curNode
                                    , children = V.fromList nodeChildren
                                    , nodeType = NetworkNode
                                    , vertName = T.pack $ "HTU" ++ show curNode
                                    , vertData = mempty
                                    , vertexResolutionData = mempty
                                    , vertexCost = 0.0
                                    , subGraphCost = 0.0
                                    }

        newDisplayNode = (curNode, newMinVertex)
        newLEdge = (curNode, index childLabel, newEdgeLabel)

        curNodeResolutionData = addNodeAndEdgeToResolutionData newDisplayNode newLEdge childResolutionData


        newVertexLabel  = VertexInfo { index = curNode
                                , bvLabel = bvLabel childLabel
                                , parents = V.fromList $ LG.parents simpleGraph curNode
                                , children = V.fromList nodeChildren
                                , nodeType = NetworkNode
                                , vertName = T.pack $ "HTU" ++ show curNode
                                , vertData = mempty
                                , vertexResolutionData = curNodeResolutionData
                                , vertexCost = 0.0
                                , subGraphCost = subGraphCost childLabel
                                }


    in
    newVertexLabel

-- | getOutDegree2VertexSoftWired returns new vertex only from two child nodes for soft-wired resolutions
getOutDegree2VertexSoftWired :: GlobalSettings
                             -> V.Vector (V.Vector CharInfo)
                             -> LG.Node
                             -> LG.LNode VertexInfo
                             -> LG.LNode VertexInfo
                             -> DecoratedGraph
                             -> VertexInfo
getOutDegree2VertexSoftWired inGS charInfoVectVect curNodeIndex leftChild@(leftChildIndex, _) rightChild@(rightChildIndex, _) inGraph =

    let -- this ensures that left/right choices are based on leaf BV for consistency and label invariance
        -- larger bitvector is Right, smaller or equal Left
        ((leftChild', leftChildLabel'), (rightChild', rightChildLabel')) = U.leftRightChildLabelBVNode (leftChild, rightChild)

        -- create resolution caches for blocks
        leftChildNodeType  = nodeType leftChildLabel'
        rightChildNodeType = nodeType rightChildLabel'
        resolutionBlockVL = V.zipWith3 (createBlockResolutions (compressResolutions inGS) curNodeIndex leftChild' rightChild' leftChildNodeType rightChildNodeType TreeNode) (vertexResolutionData leftChildLabel') (vertexResolutionData rightChildLabel') charInfoVectVect

        -- create canonical Decorated Graph vertex
        -- 0 cost becasue can't know cosrt until hit root and get best valid resolutions
        newVertexLabel = VertexInfo {  index = curNodeIndex
                                    , bvLabel = bvLabel leftChildLabel' .|. bvLabel rightChildLabel'
                                    , parents = V.fromList $ LG.parents inGraph curNodeIndex
                                    , children = V.fromList [leftChildIndex, rightChildIndex]
                                    , nodeType = TreeNode
                                    , vertName = T.pack $ "HTU" ++ show curNodeIndex
                                    , vertData = mempty --empty because of resolution data
                                    , vertexResolutionData = resolutionBlockVL
                                    , vertexCost = 0.0 --newCost
                                    , subGraphCost = 0.0 -- (subGraphCost leftChildLabel) + (subGraphCost rightChildLabel) + newCost
                                    }
        in
        newVertexLabel

-- | extractDisplayTrees takes resolutions and pulls out best cost (head for now) need to change type for multiple best
-- option for filter based on pop-count for root cost and complete display tree check
extractDisplayTrees :: Maybe Int -> Bool -> V.Vector ResolutionBlockData -> (V.Vector [DecoratedGraph], VertexCost)
extractDisplayTrees startVertex checkPopCount inRBDV =
    if V.null inRBDV then (V.empty, 0.0)
    else
        let (bestBlockDisplayResolutionList, costVect) = V.unzip $ fmap (getBestResolutionList startVertex checkPopCount) inRBDV
        in
        (bestBlockDisplayResolutionList, V.sum costVect)

-- | createBlockResolutions takes left and right child resolution data for a block (same display tree)
-- and generates node resolution data
createBlockResolutions :: Bool
                       -> LG.Node
                       -> Int
                       -> Int
                       -> NodeType
                       -> NodeType
                       -> NodeType
                       -> ResolutionBlockData
                       -> ResolutionBlockData
                       -> V.Vector CharInfo
                       -> ResolutionBlockData
createBlockResolutions
  compress
  curNode
  leftIndex
  rightIndex
  leftChildNodeType
  rightChildNodeType
  curNodeNodeType
  leftChild
  rightChild
  charInfoV
  | null leftChild && null rightChild = mempty
  | null leftChild = rightChild
  | null rightChild = leftChild
  | otherwise =
    -- trace ("CBR:" ++ (show (leftIndex, leftChildNodeType, rightIndex, rightChildNodeType)) ++ (show $fmap BV.toBits $ fmap displayBVLabel leftChild) ++ " and " ++  (show $fmap BV.toBits $ fmap displayBVLabel rightChild)) (
    let childResolutionPairs = cartProd (V.toList leftChild) (V.toList rightChild)
        -- need to keep these indices correct (hence reverse in checkLeafOverlap) for traceback and compress
        childResolutionIndices = cartProd [0.. (length leftChild - 1)] [0.. (length rightChild - 1)]
        validPairs = checkLeafOverlap (zip childResolutionPairs childResolutionIndices) []
        newResolutionList = fmap (createNewResolution curNode leftIndex rightIndex leftChildNodeType rightChildNodeType charInfoV) validPairs

        --need to add in node and edge to left and right
        edgeLable = EdgeInfo { minLength = 0.0
                             , maxLength = 0.0
                             , midRangeLength = 0.0
                             , edgeType = TreeEdge
                             }
        newMinVertex = VertexInfo { index = curNode
                                  , bvLabel = BV.fromBits [False]
                                  , parents = mempty
                                  , children = mempty
                                  , nodeType = curNodeNodeType
                                  , vertName = T.pack $ "HTU" ++ show curNode
                                  , vertData = mempty
                                  , vertexResolutionData = mempty
                                  , vertexCost = 0.0
                                  , subGraphCost = 0.0
                                  }

        newNode = (curNode, newMinVertex)

        addLeft =  if leftChildNodeType == NetworkNode then
                        let newEdge = (curNode, rightIndex, edgeLable)
                            newRightChildBlockResolutionData = addNodeEdgeToResolutionList newNode newEdge 0 [] rightChild
                        in
                        -- trace ("ANEL:" ++ (show $ (curNode, rightIndex)))
                        newRightChildBlockResolutionData
                   else -- trace ("ANEL-Nothing")
                        mempty

        addRight = if rightChildNodeType == NetworkNode then
                        let newEdge = (curNode, leftIndex, edgeLable)
                            newLeftChildBlockResolutionData = addNodeEdgeToResolutionList newNode newEdge 0 [] leftChild
                        in
                        -- trace ("ANER:" ++ (show $ (curNode, leftIndex)))
                        newLeftChildBlockResolutionData
                   else -- trace ("ANER-Nothing")
                        mempty
    in
    -- trace ("=> " ++ (show $fmap BV.toBits $ fmap displayBVLabel totalResolutions) )
    -- compress to unique resolutions--can loose display trees, but speed up post-order a great deal
    -- trace ("CBR Num res left: " ++ (show $ V.length leftChild) ++ " Num res right: " ++ (show $ V.length rightChild) ++ " =>NRL " ++ (show $ length newResolutionList) ++ " addleft " ++ (show $ length addLeft)  ++ " addright " ++ (show $ length addRight)) (
    if compress then  V.fromList (nubResolutions newResolutionList []) V.++ (addLeft V.++ addRight)
    else V.fromList newResolutionList V.++ (addLeft V.++ addRight)
    -- )

-- | createNewResolution takes a pair of resolutions and creates the median resolution
-- need to watch let/right (based on BV) for preorder stuff
createNewResolution :: LG.Node -> Int -> Int -> NodeType -> NodeType -> V.Vector CharInfo -> ((ResolutionData, ResolutionData),(Int, Int)) -> ResolutionData
createNewResolution curNode leftIndex rightIndex leftChildNodeType rightChildNodeType charInfoV ((leftRes, rightRes), (leftResIndex, rightResIndex)) =
    let -- make  bvLabel for resolution
        resBV = displayBVLabel leftRes .|. displayBVLabel rightRes

        -- Make resolution Display tree infomation
        leftEdgeType
          | leftChildNodeType == NetworkNode = NetworkEdge
          | leftChildNodeType == LeafNode = PendantEdge
          | otherwise = TreeEdge
        rightEdgeType
          | rightChildNodeType == NetworkNode = NetworkEdge
          | rightChildNodeType == LeafNode = PendantEdge
          | otherwise = TreeEdge

        edgeLable = EdgeInfo { minLength = 0.0
                             , maxLength = 0.0
                             , midRangeLength = 0.0
                             , edgeType = TreeEdge
                             }

        leftEdge =  (curNode, leftIndex, edgeLable {edgeType = leftEdgeType})
        rightEdge = (curNode, rightIndex, edgeLable {edgeType = rightEdgeType})
        leftChildTree = displaySubGraph leftRes
        rightChildTree = displaySubGraph rightRes

        -- Data fields empty for display tree data--not needed and muptiple copies of everything
        newNodeLabel = VertexInfo { index = curNode
                                  , bvLabel = resBV
                                  , parents = V.empty
                                  , children = V.fromList [leftIndex, rightIndex]
                                  , nodeType = TreeNode
                                  , vertName = T.pack $ "HTU" ++ show curNode
                                  , vertData = mempty
                                  , vertexResolutionData = mempty
                                  , vertexCost = 0.0
                                  , subGraphCost = 0.0
                                  }

        newNode = (curNode, newNodeLabel)

        resolutionEdgeList = leftEdge : (rightEdge: (snd leftChildTree ++ snd rightChildTree))
        resolutionNodeList = newNode : (fst leftChildTree ++ fst rightChildTree)

        -- Make the data and cost for the resolution
        leftBlockLength = V.length $ displayData leftRes
        rightBlockLength = V.length $ displayData rightRes
        resolutionMedianCostV
          | (leftBlockLength == 0) = V.zip (displayData rightRes) (V.replicate rightBlockLength 0)
          | (rightBlockLength == 0) = V.zip (displayData leftRes) (V.replicate leftBlockLength 0)
          | otherwise = M.median2 (displayData leftRes) (displayData rightRes) charInfoV
        (resolutionMedianV, resolutionCostV) = V.unzip resolutionMedianCostV
        thisResolutionCost = V.sum resolutionCostV
        displaySubTreeCost = displayCost leftRes + displayCost rightRes + thisResolutionCost

    in
    ResolutionData { displaySubGraph = (resolutionNodeList, resolutionEdgeList)
                   , displayBVLabel = resBV
                   , displayData = resolutionMedianV
                   , childResolutions = [(Just leftResIndex, Just rightResIndex)]
                   , resolutionCost = thisResolutionCost
                   , displayCost = displaySubTreeCost
                   }

-- | nubResolutions takes a list of resolulutions and returns 'nub' based on leaf set (bvLabel) and Charinfo vector
nubResolutions :: [ResolutionData] -> [ResolutionData] -> [ResolutionData]
nubResolutions inData curData
  | null inData = reverse curData
  | null curData = nubResolutions (tail inData) [head inData]
  | otherwise =
    let firstData = head inData
        (isUnique, matchIndex) = hasResolutionMatch (displayBVLabel firstData) (displayData firstData) curData 0
    in
    if isUnique then nubResolutions (tail inData) (firstData : curData)
    else
        -- need to update the childResolutoin field of the match
        let firstPart = L.take matchIndex curData
            lastPart = L.drop (matchIndex + 1) curData
            matchChildResolutionData =  childResolutions (curData !! matchIndex)
            firstDataChildResolutionData = childResolutions firstData
            newMatchResolution = (curData !! matchIndex) {childResolutions = firstDataChildResolutionData ++ matchChildResolutionData}
        in
        nubResolutions (tail inData) (firstPart ++ (newMatchResolution : lastPart))

-- | hasResolutionMatch checks for match of bvlabel and Vect charinfo with list
hasResolutionMatch :: NameBV -> V.Vector CharacterData -> [ResolutionData] -> Int -> (Bool, Int)
hasResolutionMatch inBV inCD rDList curIndex =
    if null rDList then (True, curIndex)
    else
        let existingBV =  displayBVLabel $ head rDList
            existingCharData = displayData $ head rDList
        in
        if (existingBV /= inBV) || (existingCharData /= inCD) then hasResolutionMatch inBV inCD (tail rDList) (curIndex + 1) else (False, curIndex)

-- | checkLeafOverlap takes a left right resolution pair list and checks if
-- there is leaf overlap via comparing displayBVLabel if & = 0 then no
-- overlap, and adds to resulting list--reverses order--sholdn't matter
checkLeafOverlap :: [((ResolutionData, ResolutionData), (Int, Int))] -> [((ResolutionData, ResolutionData), (Int, Int))] -> [((ResolutionData, ResolutionData), (Int, Int))]
checkLeafOverlap inPairList curPairList =
    if null inPairList then reverse curPairList
    else
        let inPair@((leftRes, rightRes), (_, _)) = head inPairList
            leftBV = displayBVLabel leftRes
            rightBV = displayBVLabel rightRes
        in
        if BV.isZeroVector (leftBV .&. rightBV) then checkLeafOverlap (tail inPairList) (inPair : curPairList)
        else checkLeafOverlap (tail inPairList) curPairList




-- | addNodeAndEdgeToResolutionData adds new node and edge to resolution data in outdegree = 1 nodes
-- striaght copy would not add this node or edge to subtree in resolutions
addNodeAndEdgeToResolutionData :: LG.LNode VertexInfo -> LG.LEdge EdgeInfo -> V.Vector ResolutionBlockData -> V.Vector ResolutionBlockData
addNodeAndEdgeToResolutionData newNode newEdge = fmap (addNodeEdgeToResolutionList newNode newEdge 0 [])

-- | addNodeEdgeToResolutionList adds new node and edge to single subGraph in ResolutionData
-- adds resolutoin pairs to be equal to the child straight one-for-one correpondance
addNodeEdgeToResolutionList :: LG.LNode VertexInfo -> LG.LEdge EdgeInfo -> Int -> [ResolutionData] -> V.Vector ResolutionData -> V.Vector ResolutionData
addNodeEdgeToResolutionList newNode newEdge resolutionIndex curData inData =
    if null inData then V.fromList $ reverse  curData
    else
        let firstInData = V.head inData
            (inNodeList, inEdgeList) = displaySubGraph firstInData
            -- childResolutionIndexPairList = childResolutions firstInData
            newNodeList = newNode : inNodeList
            newEdgeList = newEdge : inEdgeList
            newFirstData = firstInData { displaySubGraph  = (newNodeList, newEdgeList)
                                        -- this apir in case of left/right issues later
                                        -- not sure this is correct--LWys first for children of out = 1 node
                                       , childResolutions = [(Just 0, Just 0)] -- [(Just resolutionIndex, Just resolutionIndex)]
                                       }
    in
    addNodeEdgeToResolutionList newNode newEdge (resolutionIndex + 1) (newFirstData : curData) (V.tail inData)


-- | getAllResolutionList takes ResolutionBlockData and retuns a list of the all valid (ie all leaves in subtree) display trees
-- for that block- and costs
getAllResolutionList :: ResolutionBlockData -> [(DecoratedGraph, VertexCost)]
getAllResolutionList  inRDList =
    --trace ("GBRL: " ++ (show inRDList)) (
    if null inRDList then error "Null resolution list"
    else
        let displayTreeList = fmap displaySubGraph inRDList
            displayCostList = fmap displayCost inRDList
            displayPopList = fmap (complement . displayBVLabel) inRDList
        in
            let displayBVList = V.zip3 displayTreeList displayCostList displayPopList
                validDisplayList = V.filter (BV.isZeroVector . thd3) displayBVList
                (displayList, costList, _) = V.unzip3 validDisplayList
            in
            --trace ("Valid display list number:" ++ (show $ length validDisplayList)) (
            if V.null validDisplayList then error ("Null validDisplayList in getAllResolutionList" ++ show inRDList)
            else
                let lDisplayTreeList = fmap LG.mkGraphPair (V.toList displayList)
                    -- displayTreeList' = fmap (updateRootCost validMinCost) displayTreeList
                in
                zip lDisplayTreeList (V.toList costList)


-- | getBestResolutionList takes ResolutionBlockData and retuns a list of the best valid (ie all leaves in subtree) display trees
-- for that block-- if checkPopCount is True--otherwise all display trees of any cost and contitution
-- startVertex for a component-- to allow for not every leaf being in componnet but still getting softwired cost
getBestResolutionList :: Maybe Int -> Bool -> ResolutionBlockData -> ([DecoratedGraph], VertexCost)
getBestResolutionList startVertex checkPopCount inRDList =
    --trace ("GBRL: " ++ (show inRDList)) (
    if null inRDList then error "Null resolution list"
    else
        let displayTreeList = fmap displaySubGraph inRDList
            displayCostList = fmap displayCost inRDList
            displayPopList = fmap (complement . displayBVLabel) inRDList
        in
        if not checkPopCount then
            let minCost = minimum displayCostList
                displayCostPairList = V.zip displayTreeList displayCostList
                (bestDisplayList, _) = V.unzip $ V.filter ((== minCost) . snd) displayCostPairList
            in
            (fmap LG.mkGraphPair (V.toList bestDisplayList), minCost)
        else
            let minPopCount = minimum $ fmap popCount displayPopList  --max since complemented above
                displayBVList = V.zip3 displayTreeList displayCostList displayPopList

                -- must have all leaves if startvzertex == Nothing, component maximum otherwise
                -- this for getting cost of component of a softwired network
                validDisplayList = if startVertex == Nothing then V.filter (BV.isZeroVector . thd3) displayBVList
                                   else V.filter ((== minPopCount) . (popCount . thd3)) displayBVList
                validMinCost = V.minimum $ fmap snd3 validDisplayList
                (bestDisplayList, _, _) = V.unzip3 $ V.filter ((== validMinCost) . snd3) validDisplayList
            in
            --trace ("Valid display list number:" ++ (show $ length validDisplayList)) (
            if (startVertex == Nothing) && (V.null validDisplayList) then error ("Null root validDisplayList in getBestResolutionList" ++ (show (startVertex,inRDList)) ++ " This can be caused if the graphType not set correctly.")
            else
                let lDisplayTreeList = fmap LG.mkGraphPair (V.toList bestDisplayList)

                    -- update root cost of display trees for use later (e.g. net penalties, outputting display forrests)
                    lDisplayTreeList' = fmap (updateRootCost validMinCost) lDisplayTreeList
                in
                (lDisplayTreeList', validMinCost)
            -- )

            -- )

-- | getBestResolutionListPair takes ResolutionBlockData and retuns a list of the best valid (ie all leaves in subtree) display trees
-- for that block-- if checkPopCount is True--otherwise all display trees of any cost and contitution
-- startVertex for a component-- to allow for not every leaf being in componnet but still getting softwired cost
-- returns list of pairs
getBestResolutionListPair :: Maybe Int -> Bool -> ResolutionBlockData -> [(DecoratedGraph, VertexCost)]
getBestResolutionListPair startVertex checkPopCount inRDList =
    --trace ("GBRL: " ++ (show inRDList)) (
    if null inRDList then error "Null resolution list"
    else
        let displayTreeList = fmap displaySubGraph inRDList
            displayCostList = fmap displayCost inRDList
            displayPopList = fmap (complement . displayBVLabel) inRDList
        in
        if not checkPopCount then
            let minCost = minimum displayCostList
                displayCostPairList = V.zip displayTreeList displayCostList
                (bestDisplayList, minCostList) = V.unzip $ V.filter ((== minCost) . snd) displayCostPairList
            in
            zip (fmap LG.mkGraphPair (V.toList bestDisplayList)) (V.toList minCostList)
        else
            let minPopCount = minimum $ fmap popCount displayPopList  --max since complemented above
                displayBVList = V.zip3 displayTreeList displayCostList displayPopList

                -- must have all leaves if startvzertex == Nothing, component maximum otherwise
                -- this for getting cost of component of a softwired network
                validDisplayList = if startVertex == Nothing then V.filter (BV.isZeroVector . thd3) displayBVList
                                   else V.filter ((== minPopCount) . (popCount . thd3)) displayBVList
                validMinCost = V.minimum $ fmap snd3 validDisplayList
                (bestDisplayList, minCostList, _) = V.unzip3 $ V.filter ((== validMinCost) . snd3) validDisplayList
            in
            --trace ("Valid display list number:" ++ (show $ length validDisplayList)) (
            if (startVertex == Nothing) && (V.null validDisplayList) then error ("Null root validDisplayList in getBestResolutionListPair" ++ (show (startVertex,inRDList)) ++ " This can be caused if the graphType not set correctly.")
            else
                let lDisplayTreeList = fmap LG.mkGraphPair (V.toList bestDisplayList)

                    -- update root cost of display trees for use later (e.g. net penalties, outputting display forrests)
                    lDisplayTreeList' = fmap (updateRootCost validMinCost) lDisplayTreeList
                in
                zip lDisplayTreeList' (V.toList minCostList)
            -- )

            -- )

-- | updateRootCost updates the subGraphCost of the root node(s) with input value
-- new node is created, so original is deleted, new added, and original edges added back
-- since deleted when node is
-- assumes its a tree wiht a single root
updateRootCost :: VertexCost -> DecoratedGraph -> DecoratedGraph
updateRootCost newRootCost inGraph =
    let (rootIndex, rootLabel) = head $ LG.getRoots inGraph
        rootEdges = LG.out inGraph rootIndex
        newRootLabel = rootLabel {subGraphCost = newRootCost}
    in
    -- trace ("DCC: " ++ (show newRootCost))
    LG.insEdges rootEdges $ LG.insNode (rootIndex, newRootLabel) $ LG.delNode rootIndex inGraph

-- | makeCharacterGraph takes a blockGraph and creates a vector of character graphs
-- each with a single block and single character
-- updating costs
makeCharacterGraph :: DecoratedGraph -> V.Vector DecoratedGraph
makeCharacterGraph inBlockGraph =
  if LG.isEmpty inBlockGraph then V.empty
  else
    let numCharacters =  V.length $ V.head $ vertData $ snd $ head $ LG.labNodes inBlockGraph
        characterGraphList = if numCharacters > 0 then fmap (pullCharacter False inBlockGraph) [0.. (numCharacters - 1)]
                             -- missing data
                             else [pullCharacter True inBlockGraph 0]
    in
    if V.length (vertData $ snd $ head $ LG.labNodes inBlockGraph) /= 1 then error ("Number of blocks /= 1 in makeCharacterGraph :" ++ (show $ V.length (vertData $ snd $ head $ LG.labNodes inBlockGraph)))
    else
      -- trace ("Chars: " ++ show numCharacters)
      V.fromList characterGraphList

-- | pullCharacter takes a DecoratedGraph with a single block and
-- creates a new DecoratedGraph with a single character from the input index
pullCharacter :: Bool -> DecoratedGraph -> Int -> DecoratedGraph
pullCharacter isMissing inBlockGraph characterIndex =
  if LG.isEmpty inBlockGraph then LG.empty
  else
    let (inNodeIndexList, inNodeLabelList) = unzip $ LG.labNodes inBlockGraph
        characterLabelList = fmap (makeCharacterLabels isMissing characterIndex) inNodeLabelList
    in
    LG.mkGraph (zip inNodeIndexList characterLabelList) (LG.labEdges inBlockGraph)

-- | makeCharacterLabels pulls the index character label form the singleton block (via head)
-- and creates a singleton character label, updating costs to that of the character
-- NB the case of missing data is answered here by an "empty charcter"
-- could be better to have V.empty
-- isMIssingChar seems to be extraneous--not sure whey it was there.
makeCharacterLabels :: Bool -> Int -> VertexInfo -> VertexInfo
makeCharacterLabels isMissing characterIndex inVertexInfo =
  -- trace ("MCl in:" ++ (show inVertexInfo) ++ " " ++ (show characterIndex)) (
  let -- isMissingChar = (V.length $ (vertData inVertexInfo) V.! characterIndex) == 0
      newVertexData = V.head (vertData inVertexInfo) V.! characterIndex
      (newVertexCost, newSubGraphCost) = if isMissing then (0, 0)
                                         --if isMissing || isMissingChar then (0, 0)
                                         else (localCost newVertexData, globalCost newVertexData)
      -- newVertexCost = localCost newVertexData
      -- newSubGraphCost = globalCost newVertexData
  in
  -- trace ("MCL " ++ (show $ V.length $ vertData inVertexInfo) ++ " " ++ (show $ fmap  V.length $ vertData inVertexInfo) ) (
  -- trace ("MCL: " ++ (show isMissing) ++ " CI: " ++ (show characterIndex) ++ " " ++ (show $ V.length $ (vertData inVertexInfo) V.! characterIndex))
  inVertexInfo { vertData     = if not isMissing then V.singleton $ V.singleton newVertexData
                                -- if not isMissing && not isMissingChar then V.singleton $ V.singleton newVertexData
                                else V.singleton $ V.singleton emptyCharacter --V.empty
               , vertexCost   = newVertexCost
               , subGraphCost = newSubGraphCost
               }


    -- ) )
