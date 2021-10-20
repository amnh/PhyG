{- |
Module      :  Traversals.hs
Description :  Module specifying graph traversal functions for PhyGraph
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

module GraphOptimization.Traversals ( postOrderTreeTraversal
                                    , preOrderTreeTraversal
                                    , makeLeafGraph
                                    , multiTraverseFullyLabelTree
                                    , multiTraverseFullyLabelGraph
                                    ) where

import           Types.Types
import qualified Data.List as L
import qualified Data.Vector as V
import qualified Utilities.LocalGraph as LG
import GeneralUtilities
import qualified GraphFormatUtilities as GFU
-- import qualified GraphOptimization.Medians as M
import qualified Graphs.GraphOperations  as GO
import qualified GraphOptimization.PostOrderFunctions as PO
import qualified GraphOptimization.PreOrderFunctions as PRE
-- import qualified SymMatrix as SM
-- import qualified Data.BitVector.LittleEndian as BV
import qualified Data.Text.Lazy  as T
import Data.Bits 
import Data.Maybe
import Debug.Debug as D
import Utilities.Utilities as U
import qualified GraphOptimization.Medians as M
import qualified SymMatrix                   as S
import qualified Data.Vector.Generic         as GV
import Debug.Trace
import qualified Data.BitVector.LittleEndian as BV
import           Data.Word



-- | multiTraverseFullyLabelGraph is a wrapper around multi-traversal functions for Tree, 
-- Soft-wired network graph, and Hard-wired network graph
multiTraverseFullyLabelGraph :: GlobalSettings -> ProcessedData -> SimpleGraph -> PhylogeneticGraph
multiTraverseFullyLabelGraph inGS inData inGraph =
    if LG.isEmpty inGraph then emptyPhylogeneticGraph
    else if (graphType inGS == Tree) then 
        -- test for Tree
        let (_, _, _, networkVertexList) = LG.splitVertexList inGraph
        in
        if (null networkVertexList) then multiTraverseFullyLabelTree inGS inData inGraph
        else errorWithoutStackTrace ("Input graph is not a tree/forest, but graph type has been specified (perhaps by default) as Tree. Modify input graph or use 'set()' command to specify network type")
    else if (graphType inGS == SoftWired) then multiTraverseFullyLabelSoftWired  inGS inData inGraph
    else if (graphType inGS == HardWired) then errorWithoutStackTrace "Hard-wired graph optimization not yet supported"
    else errorWithoutStackTrace ("Unknown graph type specified: " ++ show (graphType inGS))

-- | multiTraverseFullyLabelSoftWired fully labels a softwired network component forest
-- including traversal rootings 
-- allows indegree=outdegree=1 vertices
multiTraverseFullyLabelSoftWired :: GlobalSettings -> ProcessedData -> SimpleGraph -> PhylogeneticGraph
multiTraverseFullyLabelSoftWired inGS inData inSimpleGraph = 
    if LG.isEmpty inSimpleGraph then emptyPhylogeneticGraph
    else 
        let -- starting leaves
            leafGraph = makeLeafGraphSoftWired inData

            -- for special casing of nonexact and single exact characters
            nonExactChars = U.getNumberNonExactCharacters (thd3 inData)
            exactCharacters = U.getNumberExactCharacters (thd3 inData)

            -- post order pass
            outgroupRootedSoftWiredPostOrder = postOrderSoftWiredTraversal inGS inData leafGraph inSimpleGraph -- $ GO.rerootGraph' inSimpleGraph (outgroupIndex inGS)
            childrenOfRoot = concat $ fmap (LG.descendants (thd6 outgroupRootedSoftWiredPostOrder)) (fmap fst $ LG.getRoots $ thd6 outgroupRootedSoftWiredPostOrder)
            grandChildrenOfRoot = concat $ fmap (LG.descendants (thd6 outgroupRootedSoftWiredPostOrder)) childrenOfRoot

            -- create list of multi-traversals with original rooting first
            -- subsequent rerooting do not reoptimize exact characters (add nonadd etc) 
            -- they are taken from the fully labelled first input decorated graph later when output graph created
            -- it is important that the first graph be the ourgroup rooted graph (outgroupRootedPhyloGraph) so this 
            -- will have the preoder assignmentsd for th eoutgroup rooted graph as 3rd field.  This can be used for incremental
            -- optimization to get O(log n) initial postorder assingment when mutsating graph.
            recursiveRerootList = outgroupRootedSoftWiredPostOrder : minimalReRootPhyloGraph SoftWired outgroupRootedSoftWiredPostOrder grandChildrenOfRoot

            minCostRecursive = minimum $ fmap snd6 recursiveRerootList
            minCostGraphListRecursive = filter ((== minCostRecursive).snd6) recursiveRerootList

            -- create optimal final graph with best costs and best traversal (rerooting) forest for each character
            -- traversal for exact characters (and costs) are the first of each least since exact only optimizaed for that 
            -- traversal graph.  The result has approprotate post-order assignments for traversals, preorder "final" assignments
            -- are propagated to the Decorated graph field after the preorder pass.
            graphWithBestAssignments' = L.foldl1' setBetterGraphAssignment recursiveRerootList 

            
        in
        --preOrderTreeTraversal (finalAssignment inGS) (nonExactChars > 0) outgroupRootedSoftWiredPostOrder
        --preOrderTreeTraversal (finalAssignment inGS) True  $ head minCostGraphListRecursive
        
        if nonExactChars == 0 then preOrderTreeTraversal (finalAssignment inGS) False outgroupRootedSoftWiredPostOrder
        else if (nonExactChars == 1) then preOrderTreeTraversal (finalAssignment inGS) True  $ head minCostGraphListRecursive
        else preOrderTreeTraversal (finalAssignment inGS) True  graphWithBestAssignments'
        

-- | postOrderSoftWiredTraversal performs postorder traversal on Soft-wired graph
postOrderSoftWiredTraversal :: GlobalSettings -> ProcessedData -> DecoratedGraph -> SimpleGraph -> PhylogeneticGraph
postOrderSoftWiredTraversal inGS inData@(_, _, blockDataVect) leafGraph inSimpleGraph = 
    if LG.isEmpty inSimpleGraph then emptyPhylogeneticGraph  
    else 
         -- Assumes root is Number of Leaves  
        let rootIndex = V.length $ fst3 inData
            blockCharInfo = V.map thd3 blockDataVect
            newSoftWired = postDecorateSoftWired inGS inData inSimpleGraph leafGraph blockCharInfo rootIndex
        in
        --trace ("It Begins at " ++ show rootIndex) (
        if not $ LG.isRoot inSimpleGraph rootIndex then 
            let localRootList = fmap fst $ LG.getRoots inSimpleGraph
                localRootEdges = concatMap (LG.out inSimpleGraph) localRootList
                currentRootEdges = LG.out inSimpleGraph rootIndex
            in
            error ("Index "  ++ (show rootIndex) ++ " with edges " ++ (show currentRootEdges) ++ " not root in graph:" ++ (show localRootList) ++ " edges:" ++ (show localRootEdges) ++ "\n" ++ (GFU.showGraph inSimpleGraph))
        else newSoftWired 
        -- )

-- | postDecorateSoftWired begins at start index (usually root, but could be a subtree) and moves preorder till children are labelled 
-- and then recurses to root postorder labelling vertices and edges as it goes
-- this for a single root
postDecorateSoftWired :: GlobalSettings -> ProcessedData -> SimpleGraph -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.Node -> PhylogeneticGraph
postDecorateSoftWired inGS inData simpleGraph curDecGraph blockCharInfo curNode = 
    -- if node in there nothing to do and return
    if LG.gelem curNode curDecGraph then 
        let nodeLabel = LG.lab curDecGraph curNode
        in
        if nodeLabel == Nothing then error ("Null label for node " ++ show curNode)
        else (simpleGraph, subGraphCost (fromJust nodeLabel), curDecGraph, mempty, mempty, blockCharInfo)

    else 
        -- trace ("PDSW making node " ++ show curNode ++ " in\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph curDecGraph)) (
        let nodeChildren = LG.descendants simpleGraph curNode  -- should be 1 or 2, not zero since all leaves already in graph
            leftChild = (head nodeChildren)
            rightChild = (last nodeChildren)
            leftChildTree = postDecorateSoftWired inGS inData simpleGraph curDecGraph blockCharInfo leftChild
            rightLeftChildTree = if (length nodeChildren == 2) then postDecorateSoftWired inGS inData simpleGraph (thd6 $ leftChildTree) blockCharInfo rightChild
                                 else leftChildTree
        in 
        -- Checks on children
        if length nodeChildren > 2 then error ("Graph not dichotomous in postDecorateSoftWired node " ++ (show curNode) ++ "\n" ++ LG.prettify simpleGraph)
        else if length nodeChildren == 0 then error ("Leaf not in graph in postDecorateSoftWired node " ++ (show curNode) ++ "\n" ++ LG.prettify simpleGraph)

        else
            -- make node from child block resolutions
            let newSubTree = thd6 $ rightLeftChildTree
            in

            -- single child of node (can certinly happen with soft-wired networks
            if length nodeChildren == 1 then 
                -- trace ("Outdegree 1: " ++ (show curNode) ++ " " ++ (show $ GO.getNodeType simpleGraph curNode) ++ " Child: " ++ (show nodeChildren)) (
                let (newGraph, isRoot, newVertexLabel, llocalCost, displayGraphVL) = getOutDegree1VertexAndGraph curNode (fromJust $ LG.lab newSubTree leftChild) simpleGraph nodeChildren newSubTree
                in
                if isRoot then 
                     let newGraph' = softWiredPostOrderTraceBack newGraph (curNode, newVertexLabel)

                         -- propagate post-order assignment to display trees for pre-order pass 
                         -- displayGraphVL' = fmap (assignPostOrderToDisplayTree (LG.labNodes newGraph')) $ V.zip displayGraphVL (V.fromList [0..(V.length displayGraphVL - 1)])
                         displayGraphVL' = V.zipWith (assignPostOrderToDisplayTree (fmap vertData $ fmap snd $ LG.labNodes newGraph')) displayGraphVL (V.fromList [0.. (V.length displayGraphVL - 1)])
                    in
                    (simpleGraph, llocalCost, newGraph', displayGraphVL', PO.divideDecoratedGraphByBlockAndCharacterSoftWired displayGraphVL', blockCharInfo)
                else (simpleGraph, 0, newGraph, mempty, mempty, blockCharInfo)
                -- )

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
                    resolutionBlockVL = fmap (createBlockResolutions (compressResolutions inGS) curNode leftChild' rightChild' leftChildNodeType rightChildNodeType) $  V.zip3 (vertexResolutionData leftChildLabel) (vertexResolutionData rightChildLabel) blockCharInfo
                
                    -- create canonical Decorated Graph vertex
                    -- 0 cost becasue can't know cosrt until hit root and get best valid resolutions
                    newVertexLabel = VertexInfo {  index = curNode
                                            , bvLabel = (bvLabel leftChildLabel) .|. (bvLabel rightChildLabel)
                                            , parents = V.fromList $ LG.parents simpleGraph curNode
                                            , children = V.fromList nodeChildren
                                            , nodeType = GO.getNodeType simpleGraph curNode
                                            , vertName = T.pack $ "HTU" ++ (show curNode)
                                            , vertData = mempty --empty because of resolution data
                                            , vertexResolutionData = resolutionBlockVL
                                            , vertexCost = 0.0 --newCost
                                            , subGraphCost = 0.0 -- (subGraphCost leftChildLabel) + (subGraphCost rightChildLabel) + newCost
                                            }   
                
                    leftEdgeType  = if leftChildNodeType == NetworkNode then NetworkEdge
                                    else if leftChildNodeType == LeafNode then PendantEdge
                                    else TreeEdge
                    rightEdgeType = if rightChildNodeType == NetworkNode then NetworkEdge
                                    else if rightChildNodeType == LeafNode then PendantEdge
                                    else TreeEdge

                    edgeLable = EdgeInfo { minLength = 0.0
                                         , maxLength = 0.0
                                         , midRangeLength = 0.0
                                         , edgeType = TreeEdge
                                         }

                    leftEdge =  (curNode, leftChild', edgeLable {edgeType = leftEdgeType})
                    rightEdge = (curNode, rightChild', edgeLable {edgeType = rightEdgeType})
                    newGraph =  LG.insEdges [leftEdge, rightEdge] $ LG.insNode (curNode, newVertexLabel) newSubTree 

                    (displayGraphVL, lDisplayCost) = if (nodeType newVertexLabel) == RootNode then extractDisplayTrees True resolutionBlockVL
                                                    else (mempty, 0.0)

                in
                
                if (nodeType newVertexLabel) == RootNode then 
                    -- perform "traceback" pre-order pass to choose best resolutions and assign preliminary data
                    -- to complete post-order decoration
                    let newGraph' = softWiredPostOrderTraceBack newGraph (curNode, newVertexLabel)

                        -- propagate post-order assignment to display trees for pre-order pass 
                        -- displayGraphVL' = fmap (assignPostOrderToDisplayTree (LG.labNodes newGraph')) $ V.zip displayGraphVL (V.fromList [0..(V.length displayGraphVL - 1)])
                        displayGraphVL' = V.zipWith (assignPostOrderToDisplayTree (fmap vertData $ fmap snd $ LG.labNodes newGraph')) displayGraphVL (V.fromList [0..(V.length displayGraphVL - 1)])
                    in
                    (simpleGraph, lDisplayCost, newGraph', displayGraphVL', PO.divideDecoratedGraphByBlockAndCharacterSoftWired displayGraphVL', blockCharInfo)

                else (simpleGraph, lDisplayCost, newGraph, mempty, mempty, blockCharInfo)

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
-- data ssignments.  Proceeds via typical pre-order pass over tree
softWiredPostOrderTraceBack :: DecoratedGraph -> LG.LNode VertexInfo -> DecoratedGraph
softWiredPostOrderTraceBack inGraph (rootIndex, rootLabel)  = 
    if LG.isEmpty inGraph then LG.empty
    else 
            -- get edges to remake graph after nodes are updated with preliminary states
        let inEdgeList = LG.labEdges inGraph
            
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
            rootChildrenBV = fmap bvLabel $ fmap snd rootChildren
            rootChildrenIsLeft = if length rootChildrenBV == 1 then [True]
                                 else if (rootChildrenBV !! 0) > (rootChildrenBV !! 1) then [False, True]
                                 else [True, False]
            
            rootChildrenTuples = zip3 rootChildren (replicate (length rootChildren) leftRightIndexVect) rootChildrenIsLeft
            
            -- recurse to children with resolution index from parent
            softWiredUpdatedNodes = softWiredPrelimTraceback inGraph rootChildrenTuples [newRootNode]

        in
        LG.mkGraph softWiredUpdatedNodes inEdgeList

-- | softWiredPrelimTraceback takes a list of nodes to update (with left right index info) based
-- on resolution data, recurses with left right indices pre-order to leaves, keeping list og updated nodes
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
        else if (nodeType (snd firstNode) == NetworkNode) && ((L.find ((== (fst firstNode)). fst) updatedNodes) /= Nothing) then
            softWiredPrelimTraceback inGraph (tail nodesToUpdate) updatedNodes

        else 
            let -- get children 
                firstChildren = LG.labDescendants inGraph firstNode
                firstChildrenBV = fmap bvLabel $ fmap snd firstChildren
                firstChildrenIsLeft = if length firstChildrenBV == 1 then [True]
                                      else if (firstChildrenBV !! 0) > (firstChildrenBV !! 1) then [False, True]
                                      else [True, False]

                childrenTuples = zip3 firstChildren (replicate (length firstChildren) childLeftRightIndexVect) firstChildrenIsLeft
            

        in
        {-
        if V.head resolutionIndexVect == Nothing then error ("'Nothing' index in softWiredPrelimTraceback " ++ (show resolutionIndexVect) 
            ++ " Node " ++ (show $ fst firstNode) ++ " isLeft?: " ++ (show isLeft)  ++ " " ++ (show firstLeftRight))
        else
        -} 
            softWiredPrelimTraceback inGraph (childrenTuples ++ (tail nodesToUpdate)) (newFirstNode : updatedNodes)

           
-- | getResolutionDataAndIndices takes a vertex label (VertexInfo) and returns the resolution data corresponding to
-- the index taken from its child resolution (data, subgraph cost, local resolutoin cost, left/right pairs).  
-- Index = (-1) denotes that it is a root label and in that case
-- the best (lowest cost) resolutions are returned 
getResolutionDataAndIndices :: VertexInfo -> V.Vector (Maybe Int) -> (VertexBlockData, VertexCost, VertexCost, V.Vector (Maybe Int, Maybe Int))
getResolutionDataAndIndices nodeLabel parentResolutionIndexVect = 

    -- should not happen
    --mtrace ("NL " ++ (show $ index nodeLabel) ++ " PRIL: " ++ " length " ++ (show $ V.length parentResolutionIndexVect) ++ " " ++ show parentResolutionIndexVect) (
    if nodeType nodeLabel == LeafNode then 
        let leafVertData = fmap displayData $ fmap (V.head) (vertexResolutionData nodeLabel)
        in
        (leafVertData, 0, 0, V.singleton (Just 0, Just 0))

    -- root node--take lowest cost 
    else if V.head parentResolutionIndexVect == Just (-1) then 
        let rootBlockResolutionPair = fmap getBestBlockResolution $ vertexResolutionData nodeLabel
            (charDataVV, subGraphCostV, resCostV, leftRightIndexVect) = V.unzip4 rootBlockResolutionPair
        in
        (charDataVV, V.sum subGraphCostV, V.sum resCostV, leftRightIndexVect)
        
    -- non-root node--return the index resolution information
    else 
        let parentIndexVect = fmap fromJust parentResolutionIndexVect

            -- get resolution data from node label
            resolutionData = vertexResolutionData nodeLabel

            -- get the correct (via index) resolution data for each block
            resolutionsByBlockV = V.zipWith (V.!) resolutionData parentIndexVect

            -- get other resolution info
            charDataVV = fmap displayData resolutionsByBlockV
            lSubGraphCost = V.sum $ fmap displayCost resolutionsByBlockV
            localResolutionCost = V.sum $ fmap resolutionCost resolutionsByBlockV

            -- only takes first left right pair--although others in there
            -- uses first for preliminary asiignment--but could be more
            -- if euqla cost display trees, may have multiple possible preliminary states
            leftRightIndexVect = fmap head $ fmap childResolutions resolutionsByBlockV
        in
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
            displayPopList = fmap (complement) $ fmap displayBVLabel inResBlockData

            -- subgraph cost
            displayCostList = fmap displayCost inResBlockData

            -- resolution local cost
            resolutionCostList = fmap resolutionCost inResBlockData
            
            -- takes only first resolution index pair
            childResolutionList = fmap head $ fmap childResolutions inResBlockData
            
            -- resolution medians
            displayDataList = fmap displayData inResBlockData

            -- take only those will all leaves in, then minimum cost
            quintVect = V.zip5 displayPopList displayCostList resolutionCostList childResolutionList displayDataList
            validVect = V.filter (BV.isZeroVector . fst5) quintVect
            validMinCost = V.minimum $ fmap snd5 validVect

            -- ONly takes first of potentially multiple soliutions to begin traceback
            (_, displayCostV, resCostV, childIndexPairV, displayMedianV) = V.unzip5 $ V.filter  ((== validMinCost) . snd5) validVect
        in
        if null validVect then error "Null valid quad in getBestBlockResolution--perhaps not root node or forest component"
        else (V.head displayMedianV, V.head displayCostV, V.head resCostV, V.head childIndexPairV)

-- | createBlockResolutions takes left and right child resolution data for a block (same display tree)
-- and generates node resolution data
createBlockResolutions :: Bool -> LG.Node -> Int -> Int -> NodeType -> NodeType -> (ResolutionBlockData, ResolutionBlockData, V.Vector CharInfo) -> ResolutionBlockData
createBlockResolutions compress curNode leftIndex rightIndex leftChildNodeType rightChildNodeType (leftChild, rightChild, charInfoV)  =
    if (null leftChild) && (null rightChild) then mempty
    else if null leftChild then rightChild
    else if null rightChild then leftChild
    else
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
                                      , nodeType = TreeNode
                                      , vertName = T.pack $ "HTU" ++ (show curNode)
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
        if compress then  (V.fromList $ nubResolutions newResolutionList []) V.++ (addLeft V.++ addRight)
        else (V.fromList newResolutionList) V.++ (addLeft V.++ addRight)
        -- )

-- | createNewResolution takes a pair of resolutions and creates the median resolution
-- need to watch let/right (based on BV) for preorder stuff
createNewResolution :: LG.Node -> Int -> Int -> NodeType -> NodeType -> V.Vector CharInfo -> ((ResolutionData, ResolutionData),(Int, Int)) -> ResolutionData
createNewResolution curNode leftIndex rightIndex leftChildNodeType rightChildNodeType charInfoV ((leftRes, rightRes), (leftResIndex, rightResIndex)) = 
    let -- make  bvLabel for resolution 
        resBV = (displayBVLabel leftRes) .|. (displayBVLabel rightRes)

        -- Make resolution Display tree infomation 
        leftEdgeType  = if leftChildNodeType == NetworkNode then NetworkEdge
                        else if leftChildNodeType == LeafNode then PendantEdge
                        else TreeEdge
        rightEdgeType = if rightChildNodeType == NetworkNode then NetworkEdge
                        else if rightChildNodeType == LeafNode then PendantEdge
                        else TreeEdge

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
                                  , vertName = T.pack $ "HTU" ++ (show curNode)
                                  , vertData = mempty
                                  , vertexResolutionData = mempty
                                  , vertexCost = 0.0
                                  , subGraphCost = 0.0
                                  }

        newNode = (curNode, newNodeLabel)

        resolutionEdgeList = leftEdge : (rightEdge: ((snd leftChildTree) ++ (snd rightChildTree)))
        resolutionNodeList = newNode : ((fst leftChildTree) ++ (fst rightChildTree))

        -- Make the data and cost for the resolution
        leftBlockLength = V.length $ displayData leftRes
        rightBlockLength = V.length $ displayData rightRes
        resolutionMedianCostV = if (leftBlockLength == 0) then V.zip (displayData rightRes) (V.replicate rightBlockLength 0)
                                else if (rightBlockLength == 0) then V.zip (displayData leftRes) (V.replicate leftBlockLength 0)
                                else M.median2 (displayData leftRes) (displayData rightRes) charInfoV
        (resolutionMedianV, resolutionCostV) = V.unzip resolutionMedianCostV
        thisResolutionCost = V.sum resolutionCostV
        displaySubTreeCost = (displayCost leftRes) + (displayCost rightRes) + thisResolutionCost

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
nubResolutions inData curData =
    if null inData then reverse curData
    else if null curData then nubResolutions (tail inData) [head inData]
    else
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
        if existingBV /= inBV then hasResolutionMatch inBV inCD (tail rDList) (curIndex + 1)
        else if existingCharData /= inCD then hasResolutionMatch inBV inCD (tail rDList)  (curIndex + 1)
        else (False, curIndex)

-- | checkLeafOverlap takes a left right resolution pair list and checks if 
-- there is leaf overlap via comparing displayBVLabel if & = 0 then no
-- overlap, and adds to resulting list--reverses order--sholdn't matter
checkLeafOverlap :: [((ResolutionData, ResolutionData), (Int, Int))] -> [((ResolutionData, ResolutionData), (Int, Int))] -> [((ResolutionData, ResolutionData), (Int, Int))] 
checkLeafOverlap inPairList curPairList = 
    if null inPairList then reverse curPairList 
    else 
        let inPair@((leftRes, rightRes), (leftResIndex, rightResIndex)) = head inPairList
            leftBV = displayBVLabel leftRes
            rightBV = displayBVLabel rightRes
        in
        if BV.isZeroVector (leftBV .&. rightBV) then checkLeafOverlap (tail inPairList) (inPair : curPairList) 
        else checkLeafOverlap (tail inPairList) curPairList

-- | getOutDegree1VertexAndGraph makes              
getOutDegree1VertexAndGraph :: LG.Node 
                            -> VertexInfo 
                            -> SimpleGraph 
                            -> [LG.Node] 
                            -> DecoratedGraph
                            -> (DecoratedGraph, Bool, VertexInfo, VertexCost, V.Vector [BlockDisplayForest])
getOutDegree1VertexAndGraph curNode childLabel simpleGraph nodeChildren subTree =
    
    --trace ("In out=1") (
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
                                    , vertName = T.pack $ "HTU" ++ (show curNode)
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
                                , vertName = T.pack $ "HTU" ++ (show curNode)
                                , vertData = mempty
                                , vertexResolutionData = curNodeResolutionData
                                , vertexCost = 0.0
                                , subGraphCost = subGraphCost childLabel
                                }   
        
        newLEdge = (curNode, index childLabel, newEdgeLabel) 
        newLNode = (curNode, newVertex)
        newDisplayNode = (curNode, newMinVertex)
        newGraph =  LG.insEdge newLEdge $ LG.insNode newLNode subTree

        (displayGraphVL, lDisplayCost) = if (nodeType newVertex) == RootNode then extractDisplayTrees True (vertexResolutionData childLabel)
                                        else (mempty, 0.0)


    in
    --trace ("NV1: " ++ show newVertex)
    (newGraph, (nodeType newVertex) == RootNode, newVertex, lDisplayCost, displayGraphVL)
    --)

-- | addNodeAndEdgeToResolutionData adds new node and edge to resolution data in outdegree = 1 nodes
-- staright copy would not add this node or edge to subtree in resolutions
addNodeAndEdgeToResolutionData :: LG.LNode VertexInfo -> LG.LEdge EdgeInfo -> V.Vector ResolutionBlockData -> V.Vector ResolutionBlockData
addNodeAndEdgeToResolutionData newNode newEdge inDataLV = 
    fmap (addNodeEdgeToResolutionList newNode newEdge 0 []) inDataLV 

-- | addNodeEdgeToResolutionList adds new node and edge to single subGraph in ResolutionData
-- adds resolutoin pairs to be equal to the child straight one-for-one correpondance
addNodeEdgeToResolutionList :: LG.LNode VertexInfo -> LG.LEdge EdgeInfo -> Int -> [ResolutionData] -> V.Vector ResolutionData -> V.Vector ResolutionData
addNodeEdgeToResolutionList newNode newEdge resolutionIndex curData inData =
    if null inData then V.fromList $ reverse  curData
    else 
        let firstInData = V.head inData
            (inNodeList, inEdgeList) = displaySubGraph firstInData
            childResolutionIndexPairList = childResolutions firstInData
            newNodeList = newNode : inNodeList
            newEdgeList = newEdge : inEdgeList
            newFirstData = firstInData { displaySubGraph  = (newNodeList, newEdgeList)
                                        -- this apir in case of left/right issues later
                                       , childResolutions = [(Just resolutionIndex, Just resolutionIndex)]
                                       }
    in
    addNodeEdgeToResolutionList newNode newEdge (resolutionIndex + 1) (newFirstData : curData) (V.tail inData) 


-- | makeNothingZeroForLeafChild 
-- this to propagate resolution indices for outdegree 1 nodes--but since leaves have (Nothing, Nothing)
-- parents of leaves must be the 0th resoluition of their child
makeNothingZeroForLeafChild :: (Maybe Int, Maybe Int) ->  (Maybe Int, Maybe Int)
makeNothingZeroForLeafChild (a, b) = 
    let a' = if a == Nothing then Just 0 else a
        b' = if b == Nothing then Just 0 else b
        in (a', b') 

-- | extractDisplayTrees takes resolutions and pulls out best cost (head for now) need to change type for multiple best
-- option for filter based on pop-count for root cost and complete display tree check
extractDisplayTrees :: Bool -> V.Vector ResolutionBlockData -> (V.Vector [BlockDisplayForest], VertexCost)
extractDisplayTrees checkPopCount inRBDV = 
    if V.null inRBDV then (V.empty, 0.0) 
    else 
        let (bestBlockDisplayResolutionList, costVect) = V.unzip $ fmap (getBestResolutionList checkPopCount) inRBDV
        in
        (bestBlockDisplayResolutionList, V.sum costVect)

-- | getBestResolutionList takes ResolutionBlockData and retuns a list of the best valid (ie all leaves in subtree) display trees 
-- for that block-- if checkPopCount is True--otherwise all display trees of any cost and contitution
getBestResolutionList :: Bool -> ResolutionBlockData -> ([BlockDisplayForest], VertexCost)
getBestResolutionList checkPopCount inRDList =
    if null inRDList then ([], 0.0)
    else 
        let displayTreeList = fmap displaySubGraph inRDList
            displayCostList = fmap displayCost inRDList
            displayPopList = fmap (complement) $ fmap displayBVLabel inRDList
        in
        if not checkPopCount then 
            let minCost = minimum displayCostList
                displayCostPairList = V.zip displayTreeList displayCostList
                (bestDisplayList, _) = V.unzip $ V.filter ((== minCost) . snd) displayCostPairList
            in
            (fmap LG.mkGraphPair (V.toList bestDisplayList), minCost)
        else
            let displayBVList = V.zip3 displayTreeList displayCostList displayPopList
                validDisplayList = V.filter (BV.isZeroVector . thd3) displayBVList
                validMinCost = V.minimum $ fmap snd3 validDisplayList
                (bestDisplayList, _, _) = V.unzip3 $ V.filter ((== validMinCost) . snd3) validDisplayList
            in
            -- trace ("GBR:" ++ (show $ length displayTreeList) ++ " " ++ (show $ length displayCostList) ++ " " ++ (show $ fmap BV.toBits displayPopList)) (
            if V.null validDisplayList then error "Null validDisplayList in getBestResolutionList"
            else (fmap LG.mkGraphPair (V.toList bestDisplayList), validMinCost)
            -- )

-- | makeLeafGraphSoftWired takes input data and creates a 'graph' of leaves with Vertex information
-- but with zero edges.  This 'graph' can be reused as a starting structure for graph construction
-- to avoid remaking of leaf vertices
-- includes leave resolution data
makeLeafGraphSoftWired :: ProcessedData -> DecoratedGraph
makeLeafGraphSoftWired (nameVect, bvNameVect, blocDataVect) =
    if V.null nameVect then error "Empty ProcessedData in makeLeafGraph"
    else 
        let leafVertexList = V.toList $ V.map (makeLeafVertexSoftWired nameVect bvNameVect blocDataVect) (V.fromList [0.. (V.length nameVect) - 1])
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
        --)


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

        blockIndexList = [0..((V.length inVertData) - 1)]
        blockDataList = fmap (inVertData V.!) blockIndexList
        resolutionDataList = modifyDisplayData defaultResolutionData blockDataList []
        resolutionData =   V.fromList $ fmap V.fromList $ fmap (:[]) resolutionDataList
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


-- | multiTraverseFullyLabelTree performs potorder on default root and other traversal foci, taking the minimum 
-- traversal cost for all nonexact charcters--the initial rooting is used for exact characters 
-- operates with Tree functions
-- need to add forest funcgtinoality--in principle just split into components and optimize them independently
-- but get into root index issues theway htis is written now. 
multiTraverseFullyLabelTree :: GlobalSettings -> ProcessedData -> SimpleGraph -> PhylogeneticGraph
multiTraverseFullyLabelTree inGS inData inSimpleGraph =
    if LG.isEmpty inSimpleGraph then emptyPhylogeneticGraph
    else 
        let leafGraph = makeLeafGraph inData

            -- for special casing of nonexact and single exact characters
            nonExactChars = U.getNumberNonExactCharacters (thd3 inData)
            exactCharacters = U.getNumberExactCharacters (thd3 inData)

            -- initial traversal based on global outgroup and the "next" traversal points as children of existing traversal
            -- here initial root.
            outgroupRootedPhyloGraph = postOrderTreeTraversal inGS inData leafGraph inSimpleGraph -- $ GO.rerootGraph' inSimpleGraph (outgroupIndex inGS)
            childrenOfRoot = concat $ fmap (LG.descendants (thd6 outgroupRootedPhyloGraph)) (fmap fst $ LG.getRoots $ thd6 outgroupRootedPhyloGraph)
            grandChildrenOfRoot = concat $ fmap (LG.descendants (thd6 outgroupRootedPhyloGraph)) childrenOfRoot

            -- create list of multi-traversals with original rooting first
            -- subsequent rerooting do not reoptimize exact characters (add nonadd etc) 
            -- they are taken from the fully labelled first input decorated graph later when output graph created
            -- it is important that the first graph be the ourgroup rooted graph (outgroupRootedPhyloGraph) so this 
            -- will have the preoder assignmentsd for th eoutgroup rooted graph as 3rd field.  This can be used for incremental
            -- optimization to get O(log n) initial postorder assingment when mutsating graph.
            recursiveRerootList = outgroupRootedPhyloGraph : minimalReRootPhyloGraph Tree outgroupRootedPhyloGraph grandChildrenOfRoot

            minCostRecursive = minimum $ fmap snd6 recursiveRerootList
            minCostGraphListRecursive = filter ((== minCostRecursive).snd6) recursiveRerootList

            -- create optimal final graph with best costs and best traversal (rerooting) forest for each character
            -- traversal for exact characters (and costs) are the first of each least since exact only optimizaed for that 
            -- traversal graph.  The result has approprotate post-order assignments for traversals, preorder "final" assignments
            -- are propagated to the Decorated graph field after the preorder pass.
            graphWithBestAssignments' = L.foldl1' setBetterGraphAssignment recursiveRerootList -- (recursiveRerootList !! 0) (recursiveRerootList !! 1) 

            -- this for debuggin purposes
            --allPreorderList = fmap (preOrderTreeTraversal (finalAssignment inGS)) recursiveRerootList

        in
        -- Uncomment this to (and comment the following three cases) avoid traversal rerooting stuff for debugging
        --preOrderTreeTraversal (finalAssignment inGS) outgroupRootedPhyloGraph
        -- trace ("Nums:" ++ (show $ length minCostGraphListRecursive) ++ " " ++ (show $ fmap fft6 allPreorderList))
        --preOrderTreeTraversal  (finalAssignment inGS) $ head minCostGraphListRecursive 
        
        -- special cases that don't require all the work
        
        if nonExactChars == 0 then preOrderTreeTraversal (finalAssignment inGS) False outgroupRootedPhyloGraph
        else if (nonExactChars == 1) then preOrderTreeTraversal (finalAssignment inGS) True  $ head minCostGraphListRecursive
        else preOrderTreeTraversal (finalAssignment inGS) True  graphWithBestAssignments'
        
        

-- | setBetterGraphAssignment takes two phylogenetic graphs and returns the lower cost optimization of each character,
-- with traversal focus etc to get best overall graph
-- since this is meant to work with graphs that have or do not have reoptimized exact (=static-Add/NonAdd/MAtrix) characters 
-- the criterion is lower cost character is taken, unless the cost is zero, then non-zero is taken
-- this function is expected to be used in a fold over a list of graphs
-- the basic comparison is over the costs of the root(s) cost for each  of the character decorated (traversal) graphs

-- May change
-- assumes that a single decorated graph comes in for each Phylogenetic graph from the fully and reroot optimize (V.singleton (V.singleton DecGraph))
-- and goes through the block-character-cost data and reassigns based on that creating a unique (although there could be more than one) decorated 
-- graph for each character in each block.
-- postorder assignments in traversal set of block character trees are NOT propagated back to first decorated graph. 
-- the third field of phylogenetic Graph is set to the 3rd fieled of the first of two inputs--so if startiong fold with outgroup
-- rooted graph--that is what stays which can be used as a preorder graph for incremental optimization
-- when perfoming that sort of operation
-- The traversal graphs
-- are used for the pre-order final assignments which will be propagated back to set those of the 3rd field decorated graph

-- this will have to be modified for solf-wired since incoming blocks will not all be the same underlying gaph
-- unclear how hardwired will be affected
setBetterGraphAssignment :: PhylogeneticGraph -> PhylogeneticGraph -> PhylogeneticGraph
setBetterGraphAssignment firstGraph@(fSimple, _, fDecGraph, fBlockDisplay, fTraversal, fCharInfo) secondGraph@(_, _, sDecGraph, _, sTraversal, _) =
    if LG.isEmpty fDecGraph then secondGraph
    else if LG.isEmpty sDecGraph then firstGraph
    else 
        --trace ("setBetter (" ++ (show fCost) ++ "," ++ (show sCost) ++ ")" ) ( -- ++ " CharInfo blocks:" ++ (show $ length fCharInfo) ++ " characters: " ++ (show $ fmap length fCharInfo) ++ " " ++ (show $ fmap (fmap name) fCharInfo)) (
        let (mergedBlockVect, costVector) = V.unzip $ V.zipWith makeBetterBlock fTraversal sTraversal
        in
         --trace ("setBetter (" ++ (show fCost) ++ "," ++ (show sCost) ++ ") ->" ++ (show $ V.sum costVector) ++ " nt:" ++ (show $ length fTraversal)
         --   ++ "length blocks " ++ (show $ fmap length fTraversal))
        (fSimple, V.sum costVector, fDecGraph, fBlockDisplay, mergedBlockVect, fCharInfo)
        --)
        
-- | makeBetterBlocktakes two verctors of traversals. Each vector contains a decorated graph (=traversla graph) for each 
-- character.  This can be a single sequence or series of exact characters
-- the function applies a character cost comparison to get the better 
makeBetterBlock :: V.Vector DecoratedGraph -> V.Vector DecoratedGraph -> ((V.Vector DecoratedGraph), VertexCost)
makeBetterBlock firstBlockGraphVect secondBlockGraphVect = 
    let (mergedCharacterVect, costVector) = V.unzip $ V.zipWith chooseBetterCharacter firstBlockGraphVect secondBlockGraphVect
    in
    (mergedCharacterVect, V.sum costVector)

-- | chooseBetterCharacter takes a pair of character decorated graphs and chooses teh "better" one as in lower cost, or non-zero cost
--  if one is zer (for exact characters) and returns the better character and cost
--  graph can have multiple roots
chooseBetterCharacter :: DecoratedGraph -> DecoratedGraph -> (DecoratedGraph, VertexCost)
chooseBetterCharacter firstGraph secondGraph =
    if LG.isEmpty firstGraph then error "Empty first graph in chooseBetterCharacter"
    else if LG.isEmpty secondGraph then error "Empty second graph in chooseBetterCharacter"
    else    
        let firstGraphCost = sum $ fmap subGraphCost $ fmap snd $ LG.getRoots firstGraph
            secondGraphCost = sum $ fmap subGraphCost $ fmap snd $ LG.getRoots secondGraph
        in
        -- trace ("Costs " ++ show (firstGraphCost, secondGraphCost)) (
        if firstGraphCost == 0 then (secondGraph, secondGraphCost) 
        else if secondGraphCost == 0 then (firstGraph, firstGraphCost)
        else if secondGraphCost < firstGraphCost then (secondGraph, secondGraphCost)
        else (firstGraph, firstGraphCost)
        -- )

-- | minimalReRootPhyloGraph takes an inialtial post-order labelled phylogenetic graph
-- and "intelligently" reroots by traversing through adjacent edges, hopefully
-- reoptimizing the minimum number of vertices each time (2) but could be more depending
-- on graph topology
-- NB--only deals with post-order assignments
minimalReRootPhyloGraph :: GraphType -> PhylogeneticGraph -> [LG.Node] -> [PhylogeneticGraph]
minimalReRootPhyloGraph localGraphType inGraph nodesToRoot =
    if null nodesToRoot then []
    else 
        let firstRerootIndex = head nodesToRoot
            nextReroots = (LG.descendants (thd6 inGraph) firstRerootIndex) ++ (tail nodesToRoot)
            newGraph = if localGraphType == Tree then PO.rerootPhylogeneticTree' False (-1) False (-1) inGraph firstRerootIndex
                       else if localGraphType == SoftWired then PO.rerootPhylogeneticNetwork' inGraph firstRerootIndex
                       else errorWithoutStackTrace ("Grpah type not implemented/recognized: " ++ show localGraphType)
        in
        --trace ("New cost:" ++ show (snd6 newGraph) ++ " vs " ++ (show $ GO.graphCostFromNodes $ thd6 newGraph))
        if fst6 newGraph == LG.empty then minimalReRootPhyloGraph localGraphType newGraph nextReroots
        else newGraph : minimalReRootPhyloGraph localGraphType newGraph nextReroots


-- | makeLeafGraph takes input data and creates a 'graph' of leaves with Vertex informnation
-- but with zero edges.  This 'graph' can be reused as a starting structure for graph construction
-- to avoid remaking of leaf vertices
makeLeafGraph :: ProcessedData -> DecoratedGraph
makeLeafGraph (nameVect, bvNameVect, blocDataVect) =
    if V.null nameVect then error "Empty ProcessedData in makeLeafGraph"
    else 
        let leafVertexList = V.toList $ V.map (makeLeafVertex nameVect bvNameVect blocDataVect) (V.fromList [0.. (V.length nameVect) - 1])
        in
        LG.mkGraph leafVertexList []

-- | makeLeafVertex makes a single unconnected vertex for a leaf
makeLeafVertex :: V.Vector NameText -> V.Vector NameBV -> V.Vector BlockData -> Int -> LG.LNode VertexInfo
makeLeafVertex nameVect bvNameVect inData localIndex =
    --trace ("Making leaf " ++ (show localIndex) ++ " Data " ++ (show $ length inData) ++ " " ++ (show $ fmap length $ fmap snd3 inData)) (
    let centralData = V.map snd3 inData 
        thisData = V.map (V.! localIndex) centralData
        newVertex = VertexInfo  { index = localIndex
                                , bvLabel = bvNameVect V.! localIndex
                                , parents = V.empty
                                , children = V.empty
                                , nodeType = LeafNode
                                , vertName =  nameVect V.! localIndex
                                , vertData = thisData
                                , vertexResolutionData = mempty
                                , vertexCost = 0.0
                                , subGraphCost = 0.0
                                }   
        in
        -- trace (show (length thisData) ++ (show $ fmap length thisData))
        (localIndex, newVertex)
        --)

-- | postOrderTreeTraversal takes a 'simple' graph and generates 'preliminary' assignments
-- vi post-order traversal, yields cost as well
-- for a binary tree only
-- depending on optimality criterion--will calculate root cost
postOrderTreeTraversal :: GlobalSettings -> ProcessedData -> DecoratedGraph -> SimpleGraph -> PhylogeneticGraph
postOrderTreeTraversal inGS inData@(_, _, blockDataVect) leafGraph inGraph = 
    if LG.isEmpty inGraph then emptyPhylogeneticGraph
    else
        -- Assumes root is Number of Leaves  
        let (rootIndex, _) = head $ LG.getRoots inGraph
            blockCharInfo = V.map thd3 blockDataVect
            newTree = postDecorateTree inGS inData inGraph leafGraph blockCharInfo rootIndex
        in
        -- trace ("It Begins at " ++ (show $ fmap fst $ LG.getRoots inGraph) ++ "\n" ++ show inGraph) (
        if not $ LG.isRoot inGraph rootIndex then 
            let localRootList = fmap fst $ LG.getRoots inGraph
                localRootEdges = concatMap (LG.out inGraph) localRootList
                currentRootEdges = LG.out inGraph rootIndex
            in
            error ("Index "  ++ (show rootIndex) ++ " with edges " ++ (show currentRootEdges) ++ " not root in graph:" ++ (show localRootList) ++ " edges:" ++ (show localRootEdges) ++ "\n" ++ (GFU.showGraph inGraph))
        else newTree 
        --)

-- | postDecorateTree begins at start index (usually root, but could be a subtree) and moves preorder till children are labelled and then reurns postorder
-- labelling vertices and edges as it goes back to root
-- this for a tree so single root
postDecorateTree :: GlobalSettings -> ProcessedData -> SimpleGraph -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.Node -> PhylogeneticGraph
postDecorateTree inGS inData simpleGraph curDecGraph blockCharInfo curNode = 
    -- if node in there (leaf) nothing to do and return
    if LG.gelem curNode curDecGraph then 
        let nodeLabel = LG.lab curDecGraph curNode
            -- identify/create the virtual edge this node would have been created from
            -- this for traversal focus use for character 
            -- impliedRootEdge = getVirtualRootEdge curDecGraph curNode
        in
        if nodeLabel == Nothing then error ("Null label for node " ++ show curNode)
        else
            --  replicating curaDecGraph with number opf blocks--but all the same for tree 
            (simpleGraph, subGraphCost (fromJust nodeLabel), curDecGraph, mempty, mempty, blockCharInfo)
            
    -- Need to make node
    else 
        
        -- check if children in graph
        let nodeChildren = LG.descendants simpleGraph curNode  -- should be 1 or 2, not zero since all leaves already in graph
            leftChild = (head nodeChildren)
            rightChild = (last nodeChildren)
            leftChildTree = postDecorateTree inGS inData simpleGraph curDecGraph blockCharInfo leftChild
            rightLeftChildTree = if (length nodeChildren == 2) then postDecorateTree inGS inData simpleGraph (thd6 $ leftChildTree) blockCharInfo rightChild
                                 else leftChildTree
            newSubTree = thd6 $ rightLeftChildTree
            (leftChildLabel, rightChildLabel) = U.leftRightChildLabelBV (fromJust $ LG.lab newSubTree leftChild, fromJust $ LG.lab newSubTree rightChild)
                   
        in 

        if length nodeChildren > 2 then error ("Graph not dichotomous in postDecorateTree node " ++ (show curNode) ++ "\n" ++ LG.prettify simpleGraph)
        else if length nodeChildren == 0 then error ("Leaf not in graph in postDecorateTree node " ++ (show curNode) ++ "\n" ++ LG.prettify simpleGraph)

        -- make node from childern
        else    
            -- trace ("Making " ++ (show curNode) ++ " from " ++ (show nodeChildren) ++ "Labels " ++ "\n" ++ (show leftChildLabel) ++ "\n" ++ (show rightChildLabel)) (
            -- make node from children and new edges to children
            -- takes characters in blocks--but for tree really all same block
            let -- leftChildLabel = fromJust $ LG.lab newSubTree leftChild
                -- rightChildLabel = fromJust $ LG.lab newSubTree rightChild

                -- this ensures that left/right choices are based on leaf BV for consistency and label invariance
                -- larger bitvector is Right, smaller or equal Left 
                
                newCharData = PO.createVertexDataOverBlocks  (vertData leftChildLabel) (vertData  rightChildLabel) blockCharInfo []
                newCost =  V.sum $ V.map (V.sum) $ V.map (V.map snd) newCharData
                newVertex = VertexInfo {  index = curNode
                                        , bvLabel = (bvLabel leftChildLabel) .|. (bvLabel rightChildLabel)
                                        , parents = V.fromList $ LG.parents simpleGraph curNode
                                        , children = V.fromList nodeChildren
                                        , nodeType = GO.getNodeType simpleGraph curNode
                                        , vertName = T.pack $ "HTU" ++ (show curNode)
                                        , vertData = V.map (V.map fst) newCharData
                                        , vertexResolutionData = mempty
                                        , vertexCost = newCost
                                        , subGraphCost = (subGraphCost leftChildLabel) + (subGraphCost rightChildLabel) + newCost
                                        }   
                newEdgesLabel = EdgeInfo {    minLength = newCost / 2.0
                                            , maxLength = newCost / 2.0
                                            , midRangeLength = newCost / 2.0
                                            , edgeType = TreeEdge
                                         }
                newEdges = fmap LG.toEdge $ LG.out simpleGraph curNode 
                newLEdges =  fmap (LG.toLEdge' newEdgesLabel) newEdges
                newGraph =  LG.insEdges newLEdges $ LG.insNode (curNode, newVertex) newSubTree 
              
            in
            -- trace ("New vertex:" ++ (show newVertex)) (
            -- Do we need to PO.divideDecoratedGraphByBlockAndCharacterTree if not root?  probbaly not
            if (nodeType newVertex) == RootNode then (simpleGraph, (subGraphCost newVertex), newGraph, mempty, PO.divideDecoratedGraphByBlockAndCharacterTree newGraph, blockCharInfo)
            else (simpleGraph, (subGraphCost newVertex), newGraph, mempty, mempty, blockCharInfo)
             
            -- ))
            

-- | preOrderTreeTraversal takes a preliminarily labelled PhylogeneticGraph
-- and returns a full labels with 'final' assignments based on character decorated graphs
-- created postorder (5th of 6 fields).  
-- the preorder states are creted by traversing the traversal DecoratedGraphs in the 5th filed of PhylogeneticGraphs
-- these are by block and character,  Exact charcaters are vectors of standard characters and each seqeunce (non-exact) 
-- has its own traversal graph. These should be treees here (could be forests later) and should all have same root (number taxa)
-- but worth checking to make sure.
-- these were creted by "splitting" them after preorder 

-- For sequence charcaters (slim/wide/huge) final states are created either by DirectOptimization or ImpliedAlignment
-- if DO--then does a  median between parent andchgaps wherre needed--then  doing a 3way state assignmeent filteringgaps out
-- if IA--then a separate post and preorder pass are donne on the slim/wide/huge/AI fields to crete full IA assignments 
-- that are then filtered of gaps and assigned to th efinal fields

-- The final states are propagated back to the second field
-- DecoratedGraph of the full Phylogenetic graph--which does NOT have the character based preliminary assignments
-- ie postorder--since those are traversal specific
-- the character specific decorated graphs have appropriate post and pre-order assignments
-- the traversal begins at the root (for a tree) and proceeds to leaves.
preOrderTreeTraversal :: AssignmentMethod -> Bool -> PhylogeneticGraph -> PhylogeneticGraph
preOrderTreeTraversal finalMethod hasNonExact inPGraph@(inSimple, inCost, inDecorated, blockDisplayV, blockCharacterDecoratedVV, inCharInfoVV) = 
    --trace ("PreO: " ++ (show finalMethod) ++ " " ++ (show $ fmap (fmap charType) inCharInfoVV)) (
    if LG.isEmpty (thd6 inPGraph) then emptyPhylogeneticGraph
    else 
        -- trace ("In PreOrder\n" ++ "Simple:\n" ++ (LG.prettify inSimple) ++ "Decorated:\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inDecorated) ++ "\n" ++ (GFU.showGraph inDecorated)) (
        -- mapped recursive call over blkocks, later characters
        let -- preOrderBlockVect = fmap doBlockTraversal $ Debug.debugVectorZip inCharInfoVV blockCharacterDecoratedVV
            preOrderBlockVect = V.zipWith (doBlockTraversal finalMethod) inCharInfoVV blockCharacterDecoratedVV

            -- if final non-exact states determined by IA then perform passes and assignments of final and final IA fields
            preOrderBlockVect' = if (finalMethod == ImpliedAlignment) && hasNonExact then V.zipWith makeIAAssignments preOrderBlockVect inCharInfoVV
                                else preOrderBlockVect    

            fullyDecoratedGraph = assignPreorderStatesAndEdges finalMethod preOrderBlockVect' inCharInfoVV inDecorated 
        in
        {-
        let blockPost = GO.showDecGraphs blockCharacterDecoratedVV
            blockPre = GO.showDecGraphs preOrderBlockVect
        in
        trace ("BlockPost:\n" ++ blockPost ++ "BlockPre:\n" ++ blockPre ++ "After Preorder\n" ++  (LG.prettify $ GO.convertDecoratedToSimpleGraph fullyDecoratedGraph))
        -}
        (inSimple, inCost, fullyDecoratedGraph, blockDisplayV, preOrderBlockVect, inCharInfoVV)
        -- )
        
-- | makeIAAssignments takes the vector of vector of character trees and (if) slim/wide/huge
-- does an additional pot and pre order pass to assign IA fileds and final fields in slim/wide/huge
makeIAAssignments :: V.Vector DecoratedGraph -> V.Vector CharInfo -> V.Vector DecoratedGraph
makeIAAssignments preOrderBlock inCharInfoV = V.zipWith makeCharacterIA preOrderBlock inCharInfoV

-- | makeCharacterIA takes an individual character postorder tree and if non-exact perform post and preorder IA passes
-- and assignment to final field in slim/wide/huge
makeCharacterIA :: DecoratedGraph -> CharInfo -> DecoratedGraph
makeCharacterIA inGraph charInfo = 
    if (charType charInfo) `notElem` nonExactCharacterTypes then inGraph
    else 
        let postOrderIATree = postOrderIA inGraph charInfo (LG.getRoots inGraph)
            preOrderIATree = preOrderIA postOrderIATree charInfo $ zip (LG.getRoots postOrderIATree) (LG.getRoots postOrderIATree)
        in
        -- trace ("MCIA roots:" ++ (show $ LG.getRoots postOrderIATree) ++ "\n" ++ (show $ LG.edges postOrderIATree))
        preOrderIATree

-- | postOrderIA performs a post-order IA pass assigning leaf preliminary states
-- from the "alignment" fields and setting HTU preliminary by calling the apropriate 2-way
-- matrix
postOrderIA :: DecoratedGraph -> CharInfo -> [LG.LNode VertexInfo] -> DecoratedGraph
postOrderIA inGraph charInfo inNodeList  =
    if null inNodeList then inGraph
    else 
        let inNode@(nodeIndex, nodeLabel) = head inNodeList
            (inNodeEdges, outNodeEdges) = LG.getInOutEdges inGraph nodeIndex
            characterType = charType charInfo
            symbols = length $ costMatrix charInfo
            inCharacter = V.head $ V.head $ vertData nodeLabel
            inCharacter' = inCharacter
        in

        -- trace ("POIA Node: " ++ (show nodeIndex) ++ " " ++ (show $ nodeType nodeLabel) ++ " " ++ (show  $ fmap fst inNodeList)) (
        -- checking sanity of data
        if V.null $ vertData nodeLabel then error "Null vertData in postOrderIA"
        else if V.null $ V.head $ vertData nodeLabel then error "Null vertData data in postOrderIA"
        
        -- leaf take assignment from alignment field
        else if (nodeType nodeLabel)  == LeafNode then
            let newCharacter = if characterType `elem` [SlimSeq, NucSeq] then 
                                    inCharacter { slimIAPrelim = slimAlignment inCharacter'
                                                , slimIAFinal = fst3 $ slimAlignment inCharacter'
                                                , slimFinal = M.createUngappedMedianSequence symbols $ slimAlignment inCharacter'
                                                }
                               else if characterType `elem` [WideSeq, AminoSeq] then 
                                    inCharacter { wideIAPrelim = wideAlignment inCharacter'
                                                , wideIAFinal = fst3 $ wideAlignment inCharacter'
                                                , wideFinal = M.createUngappedMedianSequence symbols $ wideAlignment inCharacter'
                                                }
                               else if characterType == HugeSeq then 
                                    inCharacter { hugeIAPrelim = hugeAlignment inCharacter'
                                                , hugeIAFinal = fst3 $ hugeAlignment inCharacter'
                                                , hugeFinal = M.createUngappedMedianSequence symbols $ hugeAlignment inCharacter'
                                                }
                               else error ("Unrecognized character type " ++ (show characterType))
                newLabel = nodeLabel  {vertData = V.singleton (V.singleton newCharacter)}
                newGraph = LG.insEdges (inNodeEdges ++ outNodeEdges) $ LG.insNode (nodeIndex, newLabel) $ LG.delNode nodeIndex inGraph
            in
            postOrderIA newGraph charInfo (tail inNodeList)

        -- HTU take create assignment from children
        else 
            let childNodes = LG.labDescendants inGraph inNode
                childTree = postOrderIA inGraph charInfo childNodes
            in
            --trace ("Children: " ++ (show  $ fmap fst childNodes)) (
            
            if length childNodes > 2 then error ("Too many children in postOrderIA: " ++ (show $ length childNodes))

            -- in 1 out 1 vertex
            else if length childNodes == 1 then
                let childIndex = fst $ head childNodes
                    childLabel = fromJust $ LG.lab childTree childIndex
                    childCharacter = V.head $ V.head $ vertData childLabel 
                in
                -- sanity checks
                if LG.lab childTree (fst $ head childNodes) == Nothing then error ("No label for node: " ++ (show $ fst $ head childNodes))
                else if V.null $ vertData childLabel then error "Null vertData in postOrderIA"
                else if V.null $ V.head $ vertData childLabel then error "Null vertData data in postOrderIA"
                else 
                    let newLabel = nodeLabel  {vertData = V.singleton (V.singleton childCharacter)}
                        newGraph = LG.insEdges (inNodeEdges ++ outNodeEdges) $ LG.insNode (nodeIndex, newLabel) $ LG.delNode nodeIndex childTree
                    in
                    postOrderIA newGraph charInfo (tail inNodeList)

            -- two children 
            else 
                let childIndices = fmap fst childNodes
                    childlabels = fmap fromJust $ fmap (LG.lab childTree) childIndices
                    childCharacters = fmap vertData childlabels
                    leftChar = V.head $ V.head $ head childCharacters
                    rightChar = V.head $ V.head $ last childCharacters
                    newCharacter = makeIAPrelimCharacter charInfo inCharacter leftChar rightChar
                    newLabel = nodeLabel  {vertData = V.singleton (V.singleton newCharacter)}
                    newGraph = LG.insEdges (inNodeEdges ++ outNodeEdges) $ LG.insNode (nodeIndex, newLabel) $ LG.delNode nodeIndex childTree
                in
                postOrderIA newGraph charInfo (tail inNodeList)
            -- )
    -- )

-- | makeIAPrelimCharacter takes two characters and performs 2-way assignment 
-- based on character type and nodeChar--only IA fields are modified
makeIAPrelimCharacter :: CharInfo -> CharacterData -> CharacterData -> CharacterData -> CharacterData
makeIAPrelimCharacter charInfo nodeChar leftChar rightChar =
     let characterType = charType charInfo
         symbols = (length $ costMatrix charInfo)
     in
     if characterType `elem` [SlimSeq, NucSeq] then 
        let prelimChar = PRE.get2WaySlim (slimTCM charInfo) symbols (fst3 $ slimIAPrelim leftChar) (fst3 $ slimIAPrelim rightChar)
        in
        nodeChar {slimIAPrelim = (prelimChar,  fst3 $ slimIAPrelim leftChar, fst3 $ slimIAPrelim rightChar)}
     else if characterType `elem` [WideSeq, AminoSeq] then 
        let prelimChar = PRE.get2WayWideHuge (wideTCM charInfo) symbols (fst3 $ wideIAPrelim leftChar) (fst3 $ wideIAPrelim rightChar)
        in
        nodeChar {wideIAPrelim = (prelimChar, fst3 $ wideIAPrelim leftChar, fst3 $ wideIAPrelim rightChar)}
     else if characterType == HugeSeq then 
        let prelimChar = PRE.get2WayWideHuge (hugeTCM charInfo) symbols (fst3 $ hugeIAPrelim leftChar) (fst3 $ hugeIAPrelim rightChar)
        in
        nodeChar {hugeIAPrelim = (prelimChar, fst3 $ hugeIAPrelim leftChar, fst3 $ hugeIAPrelim rightChar)}
     else error ("Unrecognized character type " ++ (show characterType))

-- | makeIAFinalharacter takes two characters and performs 2-way assignment 
-- based on character type and nodeChar--only IA fields are modified
makeIAFinalCharacter :: CharInfo -> CharacterData -> CharacterData -> CharacterData -> CharacterData -> CharacterData
makeIAFinalCharacter charInfo nodeChar parentChar leftChar rightChar =
     let characterType = charType charInfo
         symbols = (length $ costMatrix charInfo)
     in
     if characterType `elem` [SlimSeq, NucSeq] then 
        let finalIAChar = PRE.getFinal3WaySlim (slimTCM charInfo) symbols (slimIAFinal parentChar) (fst3 $ slimIAPrelim leftChar) (fst3 $ slimIAPrelim rightChar)
            finalChar =  M.createUngappedMedianSequence symbols $ (finalIAChar, finalIAChar, finalIAChar)
        in
        nodeChar { slimIAFinal = finalIAChar
                 , slimFinal = finalChar 
                 }
     else if characterType `elem` [WideSeq, AminoSeq] then 
        let finalIAChar = PRE.getFinal3WayWideHuge (wideTCM charInfo) symbols  (wideIAFinal parentChar) (fst3 $ wideIAPrelim leftChar) (fst3 $ wideIAPrelim rightChar)
            finalChar = M.createUngappedMedianSequence symbols $ (finalIAChar, finalIAChar, finalIAChar)
        in
        nodeChar { wideIAFinal = finalIAChar
                 , wideFinal = finalChar
                 }
     else if characterType == HugeSeq then 
        let finalIAChar = PRE.getFinal3WayWideHuge (hugeTCM charInfo) symbols  (hugeIAFinal parentChar) (fst3 $ hugeIAPrelim leftChar) (fst3 $ hugeIAPrelim rightChar)
            finalChar = M.createUngappedMedianSequence symbols $ (finalIAChar, finalIAChar, finalIAChar)
        in
        nodeChar { hugeIAFinal = finalIAChar
                 , hugeFinal = finalChar
                 }
     else error ("Unrecognized character type " ++ (show characterType))



-- | preOrderIA performs a pre-order IA pass assigning via the apropriate 3-way matrix
-- the "final" fields are also set by filtering out gaps and 0.
preOrderIA :: DecoratedGraph -> CharInfo -> [(LG.LNode VertexInfo, LG.LNode VertexInfo)] -> DecoratedGraph
preOrderIA inGraph charInfo inNodePairList = 
    if null inNodePairList then inGraph
    else 
        let (inNode@(nodeIndex, nodeLabel), (parentNodeIndex, parentNodeLabel)) = head inNodePairList
            (inNodeEdges, outNodeEdges) = LG.getInOutEdges inGraph nodeIndex
            characterType = charType charInfo
            symbols = (length $ costMatrix charInfo)
            inCharacter = V.head $ V.head $ vertData nodeLabel
            inCharacter' = inCharacter
            parentCharacter = V.head $ V.head $ vertData parentNodeLabel
            childNodes = LG.labDescendants inGraph inNode
        in
        --trace ("PreIA Node:" ++ (show nodeIndex) ++ " " ++ (show $ nodeType nodeLabel) ++ " " ++ (show (fmap fst $ fmap fst inNodePairList,fmap fst $ fmap snd inNodePairList))) (
        -- checking sanity of data
        if V.null $ vertData nodeLabel then error "Null vertData in preOrderIA"
        else if V.null $ V.head $ vertData nodeLabel then error "Null vertData data in preOrderIA"
        else if length childNodes > 2 then error ("Too many children in preOrderIA: " ++ (show $ length childNodes))

        -- leaf done in post-order
        else if (nodeType nodeLabel) == LeafNode then preOrderIA inGraph charInfo (tail inNodePairList)

        else if (nodeType nodeLabel) == RootNode then 
            let newCharacter = if characterType `elem` [SlimSeq, NucSeq] then 
                                    inCharacter { slimFinal = M.createUngappedMedianSequence symbols $ slimAlignment inCharacter'
                                                , slimIAFinal = fst3 $ slimAlignment inCharacter'
                                                }
                               else if characterType `elem` [WideSeq, AminoSeq] then 
                                    inCharacter { wideFinal = M.createUngappedMedianSequence symbols $ wideAlignment inCharacter'
                                                , wideIAFinal = fst3 $ wideAlignment inCharacter'
                                                }
                               else if characterType == HugeSeq then 
                                    inCharacter { hugeFinal = M.createUngappedMedianSequence symbols $ hugeAlignment inCharacter'
                                                , hugeIAFinal = fst3 $ hugeAlignment inCharacter'
                                                }
                               else error ("Unrecognized character type " ++ (show characterType))

                newLabel = nodeLabel  {vertData = V.singleton (V.singleton newCharacter)}
                newGraph = LG.insEdges (inNodeEdges ++ outNodeEdges) $ LG.insNode (nodeIndex, newLabel) $ LG.delNode nodeIndex inGraph
                parentNodeList = replicate (length childNodes) (nodeIndex, newLabel)
            in
            preOrderIA newGraph charInfo ((tail inNodePairList) ++ (zip childNodes parentNodeList))

        -- single child, take parent final assignments, but keep postorder assignments    
        else if length childNodes == 1 then 
                let newCharacter = if characterType `elem` [SlimSeq, NucSeq] then 
                                      inCharacter { slimFinal = slimFinal parentCharacter
                                                  , slimIAFinal = slimIAFinal parentCharacter
                                                  }
                                   else if characterType `elem` [WideSeq, AminoSeq] then 
                                      inCharacter { wideFinal = wideFinal parentCharacter
                                                  , wideIAFinal = wideIAFinal parentCharacter
                                                  }
                                   else if characterType == HugeSeq then 
                                      inCharacter { hugeFinal = hugeFinal parentCharacter
                                                  , hugeIAFinal = hugeIAFinal parentCharacter
                                                  }
                                   else error ("Unrecognized character type " ++ (show characterType))

                    newLabel = nodeLabel  {vertData = V.singleton (V.singleton newCharacter)}
                    newGraph = LG.insEdges (inNodeEdges ++ outNodeEdges) $ LG.insNode (nodeIndex, newLabel) $ LG.delNode nodeIndex inGraph
                    parenNodeList = replicate (length childNodes) (nodeIndex, newLabel)
                in 
                preOrderIA newGraph charInfo ((tail inNodePairList) ++ (zip childNodes parenNodeList))

        -- 2 children, make 3-way 
        else 
            let childLabels = fmap snd childNodes
                leftChar = V.head $ V.head $ vertData $ head childLabels
                rightChar = V.head $ V.head $ vertData $ last childLabels
                finalCharacter = makeIAFinalCharacter charInfo inCharacter parentCharacter leftChar rightChar 

                newLabel = nodeLabel  {vertData = V.singleton (V.singleton finalCharacter)}
                newGraph = LG.insEdges (inNodeEdges ++ outNodeEdges) $ LG.insNode (nodeIndex, newLabel) $ LG.delNode nodeIndex inGraph
                parenNodeList = replicate (length childNodes) (nodeIndex, newLabel)
            in 
            preOrderIA newGraph charInfo ((tail inNodePairList) ++ (zip childNodes parenNodeList)) 
            
        -- )

-- | doBlockTraversal takes a block of postorder decorated character trees character info  
-- could be moved up preOrderTreeTraversal, but like this for legibility
doBlockTraversal :: AssignmentMethod -> V.Vector CharInfo -> V.Vector DecoratedGraph -> V.Vector DecoratedGraph
doBlockTraversal finalMethod inCharInfoV traversalDecoratedVect =
    --trace ("BlockT:" ++ (show $ fmap charType inCharInfoV)) 
    V.zipWith (doCharacterTraversal finalMethod) inCharInfoV traversalDecoratedVect

-- | doCharacterTraversal perfoms preorder traversal on single character tree
-- with single charInfo
-- this so each character can be independently "rooted" for optimal traversals.
doCharacterTraversal :: AssignmentMethod -> CharInfo -> DecoratedGraph -> DecoratedGraph 
doCharacterTraversal finalMethod inCharInfo inGraph =
    -- find root--index should = number of leaves 
    --trace ("charT:" ++ (show $ charType inCharInfo)) (
    let rootVertexList = LG.getRoots inGraph
        (_, leafVertexList, _, _)  = LG.splitVertexList inGraph
        rootIndex = fst $ head rootVertexList
        inEdgeList = LG.labEdges inGraph
    in
    -- remove these two lines if working
    if length rootVertexList /= 1 then error ("Root number not = 1 in doCharacterTraversal" ++ show (rootVertexList))
    else if rootIndex /=  length leafVertexList then error ("Root index not =  number leaves in doCharacterTraversal" ++ show (rootIndex, length rootVertexList))
    else 
        -- root vertex, repeat of label info to avoid problem with zero length zip later, second info ignored for root
        let rootLabel = snd $ head rootVertexList
            rootFinalVertData = PRE.createFinalAssignmentOverBlocks finalMethod RootNode (vertData rootLabel) (vertData rootLabel) inCharInfo True False
            rootChildren =LG.labDescendants inGraph (head rootVertexList)

            -- left / right to match post-order
            rootChildrenBV = fmap bvLabel $ fmap snd rootChildren
            rootChildrenIsLeft = if length rootChildrenBV == 1 then [True]
                                 else if (rootChildrenBV !! 0) > (rootChildrenBV !! 1) then [False, True]
                                 else [True, False]
            newRootNode = (rootIndex, rootLabel {vertData = rootFinalVertData})
            rootChildrenPairs = zip3 rootChildren (replicate (length rootChildren) newRootNode) rootChildrenIsLeft
            upDatedNodes = makeFinalAndChildren finalMethod inGraph rootChildrenPairs [newRootNode] inCharInfo
        in
        -- hope this is the most efficient way since all nodes have been remade
        -- trace (U.prettyPrintVertexInfo $ snd newRootNode)
        LG.mkGraph upDatedNodes inEdgeList
        --)

-- | makeFinalAndChildren takes a graph, list of pairs of (labelled nodes,parent node) to make final assignment and a liss of updated nodes
-- the input nodes are relabelled by preorder functions and added to the list of processed nodes and recursed to their children
-- nodes are retuned in reverse order at they are made--need to check if this will affect graph identity or indexing in fgl
makeFinalAndChildren :: AssignmentMethod 
                     -> DecoratedGraph 
                     -> [(LG.LNode VertexInfo, LG.LNode VertexInfo, Bool)] 
                     -> [LG.LNode VertexInfo] 
                     -> CharInfo 
                     -> [LG.LNode VertexInfo]
makeFinalAndChildren finalMethod inGraph nodesToUpdate updatedNodes inCharInfo =
    --trace ("mFAC:" ++ (show $ charType inCharInfo)) (
    if null nodesToUpdate then updatedNodes
    else 
        let (firstNode, firstParent, isLeft) = head nodesToUpdate
            firstLabel = snd firstNode
            firstNodeType = nodeType firstLabel
            firstVertData = vertData firstLabel
            firstParentVertData = vertData $ snd firstParent
            firstChildren = LG.labDescendants inGraph firstNode
            
            -- this OK with one or two children
            firstChildrenBV = fmap bvLabel $ fmap snd firstChildren
            firstChildrenIsLeft = if length firstChildrenBV == 1 then [True]
                                  else if (firstChildrenBV !! 0) > (firstChildrenBV !! 1) then [False, True]
                                  else [True, False]
            firstFinalVertData = PRE.createFinalAssignmentOverBlocks finalMethod firstNodeType firstVertData firstParentVertData inCharInfo isLeft (length firstChildren == 1)
            newFirstNode = (fst firstNode, firstLabel {vertData = firstFinalVertData})
            childrenPairs = zip3 firstChildren (replicate (length firstChildren) newFirstNode) firstChildrenIsLeft
        in
        -- trace (U.prettyPrintVertexInfo $ snd newFirstNode)
        makeFinalAndChildren finalMethod inGraph (childrenPairs ++ (tail nodesToUpdate)) (newFirstNode : updatedNodes) inCharInfo
        --)

-- | assignPreorderStatesAndEdges takes a postorder decorated graph (should be but not required) and propagates 
-- preorder character states from individual character trees.  Exact characters (Add, nonAdd, matrix) postorder
-- states should be based on the outgroup rooted tree.  
-- root should be median of finals of two descendets--for non-exact based on final 'alignments' field with gaps filtered
-- postorder assignment and preorder will be out of whack--could change to update with correponding postorder
-- but that would not allow use of base decorated graph for incremental optimization (which relies on postorder assignments) in other areas
-- optyion code ikn there to set root final to outgropu final--but makes thigs scewey in matrix character and some pre-order assumptions
assignPreorderStatesAndEdges :: AssignmentMethod -> V.Vector (V.Vector DecoratedGraph) -> V.Vector (V.Vector CharInfo) -> DecoratedGraph  -> DecoratedGraph
assignPreorderStatesAndEdges finalMethd preOrderBlockTreeVV inCharInfoVV inGraph =
    --trace ("aPSAE:" ++ (show $ fmap (fmap charType) inCharInfoVV)) (
    if LG.isEmpty inGraph then error "Empty graph in assignPreorderStatesAndEdges"
    else 
        -- trace ("In assign") (
        let postOrderNodes = LG.labNodes inGraph
            postOrderEdgeList = LG.labEdges inGraph

            -- update node labels
            newNodeList = fmap (updateNodeWithPreorder preOrderBlockTreeVV inCharInfoVV) postOrderNodes
            
            -- update edge labels
            newEdgeList = fmap (updateEdgeInfo finalMethd inCharInfoVV (V.fromList $ L.sortOn fst newNodeList)) postOrderEdgeList
        in
        -- make new graph
        -- LG.mkGraph newNodeList' newEdgeList
        LG.mkGraph newNodeList newEdgeList
        --)

-- | updateNodeWithPreorder takes the preorder decorated graphs (by block and character) and updates the
-- the preorder fields only using character info.  This leaves post and preorder assignment out of sync.
-- but that so can use incremental optimizaytion on base decorated graph in other areas.
updateNodeWithPreorder :: V.Vector (V.Vector DecoratedGraph) -> V.Vector (V.Vector CharInfo) -> LG.LNode VertexInfo -> LG.LNode VertexInfo
updateNodeWithPreorder preOrderBlockTreeVV inCharInfoVV postOrderNode =
    let nodeLabel = snd postOrderNode
        nodeVertData = vertData nodeLabel
        newNodeVertData = V.zipWith3 (updateVertexBlock (fst postOrderNode)) preOrderBlockTreeVV nodeVertData inCharInfoVV 
    in
    (fst postOrderNode, nodeLabel {vertData = newNodeVertData})

-- | updateVertexBlock takes a block of vertex data and updates preorder states of charactes via fmap
updateVertexBlock :: Int -> V.Vector DecoratedGraph -> V.Vector CharacterData -> V.Vector CharInfo -> V.Vector CharacterData 
updateVertexBlock nodeIndex blockTraversalTreeV nodeCharacterDataV charInfoV =
    V.zipWith3 (updatePreorderCharacter nodeIndex) blockTraversalTreeV nodeCharacterDataV charInfoV

-- | updatePreorderCharacter updates the pre-order fields of character data for a vertex from a traversal
-- since there is single character optimized for each character decorated graph-- it is always teh 0th 0th character
-- exact are vectors so take care of multiple there.
-- need to care for issues of missing data
updatePreorderCharacter :: Int -> DecoratedGraph -> CharacterData -> CharInfo -> CharacterData 
updatePreorderCharacter nodeIndex preOrderTree postOrderCharacter charInfo =
    --trace ("N:" ++ (show nodeIndex) ++ " B:" ++ (show blockIndex) ++ " C:" ++ (show characterIndex) ++ "\n" ++ (show $ vertData $ fromJust $ LG.lab preOrderTree nodeIndex)) (
    let maybePreOrderNodeLabel = LG.lab preOrderTree nodeIndex
        preOrderVertData = vertData $ fromJust maybePreOrderNodeLabel
        preOrderCharacterData = if V.null preOrderVertData then emptyCharacter
                                else if V.null $ V.head preOrderVertData then emptyCharacter
                                else V.head $ V.head preOrderVertData -- (preOrderVertData V.! 0) V.! 0
                                
    in
    if maybePreOrderNodeLabel == Nothing then error ("Nothing node label in updatePreorderCharacter node: " ++ show nodeIndex)
    else
        updateCharacter postOrderCharacter preOrderCharacterData (charType charInfo) 
    --)

-- | updateCharacter takes a postorder character and updates the preorder (final) fields with preorder data and character type
-- only updating preorder assignment--except for root, that is needed to draw final state for brnach lengths
updateCharacter :: CharacterData -> CharacterData -> CharType  -> CharacterData
updateCharacter postOrderCharacter preOrderCharacter localCharType  =
    if localCharType == Add then
        postOrderCharacter { rangeFinal = rangeFinal preOrderCharacter }
    
    else if localCharType == NonAdd then
        postOrderCharacter { stateBVFinal = stateBVFinal preOrderCharacter }
    
    else if localCharType == Matrix then
        postOrderCharacter { matrixStatesFinal = matrixStatesFinal preOrderCharacter }
    
    else if (localCharType == SlimSeq || localCharType == NucSeq) then
        postOrderCharacter { slimAlignment = slimAlignment preOrderCharacter
                           , slimFinal = slimFinal preOrderCharacter
                           , slimIAFinal = slimIAFinal preOrderCharacter
                       }
    
    else if (localCharType == WideSeq || localCharType == AminoSeq) then
        postOrderCharacter { wideAlignment = wideAlignment preOrderCharacter
                           , wideFinal = wideFinal preOrderCharacter
                           , wideIAFinal = wideIAFinal preOrderCharacter
                           }

    else if localCharType == HugeSeq then
        postOrderCharacter { hugeAlignment = hugeAlignment preOrderCharacter
                           , hugeFinal = hugeFinal preOrderCharacter
                           , hugeIAFinal = hugeIAFinal preOrderCharacter
                           }

    else error ("Character type unimplemented : " ++ show localCharType)


-- | updateEdgeInfo takes a Decorated graph--fully labelled post and preorder and and edge and 
-- gets edge info--basically lengths
updateEdgeInfo :: AssignmentMethod -> V.Vector (V.Vector CharInfo) -> V.Vector (LG.LNode VertexInfo) -> LG.LEdge EdgeInfo -> LG.LEdge EdgeInfo 
updateEdgeInfo finalMethod inCharInfoVV nodeVector (uNode, vNode, edgeLabel) =
    if V.null nodeVector then error "Empty node list in updateEdgeInfo"
    else 
        let (minW, maxW) = getEdgeWeight finalMethod inCharInfoVV nodeVector (uNode, vNode)
            midW = (minW + maxW) / 2.0
            localEdgeType = edgeType edgeLabel
            newEdgeLabel = EdgeInfo { minLength = minW
                                    , maxLength = maxW
                                    , midRangeLength = midW
                                    , edgeType  = localEdgeType
                                    }
        in
        (uNode, vNode, newEdgeLabel)
        
-- | getEdgeWeight takes a preorder decorated decorated graph and an edge and gets the weight information for that edge
-- basically a min/max distance between the two
getEdgeWeight :: AssignmentMethod -> V.Vector (V.Vector CharInfo) -> V.Vector (LG.LNode VertexInfo) -> (Int, Int) -> (VertexCost, VertexCost)
getEdgeWeight finalMethod inCharInfoVV nodeVector (uNode, vNode) = 
    if V.null nodeVector then error "Empty node list in getEdgeWeight"
    else 
        let uNodeInfo = vertData $ snd $ nodeVector V.! uNode
            vNodeInfo = vertData $ snd $ nodeVector V.! vNode
            blockCostPairs = V.zipWith3 (getBlockCostPairs finalMethod) uNodeInfo vNodeInfo inCharInfoVV
            minCost = sum $ fmap fst blockCostPairs
            maxCost = sum $ fmap snd blockCostPairs
        in
        (minCost, maxCost)
        
-- | getBlockCostPairs takes a block of two nodes and character infomation and returns the min and max block branch costs
getBlockCostPairs :: AssignmentMethod -> V.Vector CharacterData -> V.Vector CharacterData -> V.Vector CharInfo -> (VertexCost, VertexCost)
getBlockCostPairs finalMethod uNodeCharDataV vNodeCharDataV charInfoV = 
    let characterCostPairs = V.zipWith3 (getCharacterDist finalMethod) uNodeCharDataV vNodeCharDataV charInfoV
        minCost = sum $ fmap fst characterCostPairs
        maxCost = sum $ fmap snd characterCostPairs
    in
    (minCost, maxCost)   

-- | getCharacterDist takes a pair of characters and character type, retunring teh minimum and maximum character distances
-- for sequence charcaters this is based on slim/wide/hugeAlignment field, hence all should be n in num characters/seqeunce length
getCharacterDist :: AssignmentMethod -> CharacterData -> CharacterData -> CharInfo -> (VertexCost, VertexCost)
getCharacterDist finalMethod uCharacter vCharacter charInfo =
    let thisWeight = weight charInfo
        thisMatrix = costMatrix charInfo
        thisCharType = charType charInfo
        gapChar = bit $ (length thisMatrix) - 1
        gapCharWide = (bit $ (length thisMatrix) - 1) :: Word64
        gapCharBV = (bit $ (length thisMatrix) - 1) :: BV.BitVector
    in
    if thisCharType == Add then
        let minCost = localCost (M.intervalAdd thisWeight uCharacter vCharacter)
            maxDiff = V.sum $ V.zipWith maxIntervalDiff  (rangeFinal uCharacter) (rangeFinal vCharacter)
            maxCost = thisWeight * (fromIntegral maxDiff)
        in
        (minCost, maxCost)

        
    else if thisCharType == NonAdd then
        let minCost = localCost (M.interUnion thisWeight uCharacter vCharacter)
            maxDiff = length $ V.filter (==False) $ V.zipWith (==) (stateBVFinal uCharacter) (stateBVFinal vCharacter)
            maxCost = thisWeight * (fromIntegral maxDiff)
        in
        (minCost, maxCost)
    
    else if thisCharType == Matrix then
        let minMaxListList= V.zipWith (minMaxMatrixDiff thisMatrix)  (fmap (fmap fst3) $ matrixStatesFinal uCharacter) (fmap (fmap fst3) $ matrixStatesFinal vCharacter)
            minDiff = V.sum $ fmap fst minMaxListList
            maxDiff = V.sum $ fmap snd minMaxListList
            minCost = thisWeight * (fromIntegral minDiff)
            maxCost = thisWeight * (fromIntegral maxDiff)
        in
        (minCost, maxCost)


    else if (thisCharType == SlimSeq || thisCharType == NucSeq) then
        let minMaxDiffList = if finalMethod == DirectOptimization then 
                                let uFinal = slimFinal uCharacter
                                    vFinal = slimFinal vCharacter
                                    newEdgeCharacter = M.getDOMedianCharInfo charInfo (uCharacter {slimPrelim = uFinal}) (vCharacter {slimPrelim = vFinal}) 
                                    (_, newU, newV) = slimGapped newEdgeCharacter
                                in
                                --trace ("GCD:\n" ++ (show m) ++ "\n" ++ (show (uFinal, newU)) ++ "\n" ++ (show (vFinal, newV)))
                                zipWith  (generalSequenceDiff thisMatrix (length thisMatrix)) (GV.toList $ GV.map (zero2Gap gapChar) newU) (GV.toList $ GV.map (zero2Gap gapChar) newV)
                             else zipWith  (generalSequenceDiff thisMatrix (length thisMatrix)) (GV.toList $ slimIAFinal uCharacter) (GV.toList $ slimIAFinal vCharacter)
            (minDiff, maxDiff) = unzip minMaxDiffList
            minCost = thisWeight * (fromIntegral $ sum minDiff)
            maxCost = thisWeight * (fromIntegral $ sum maxDiff)
        in
        --trace ("MMDL: " ++ (show $ (GV.toList $ slimFinal uCharacter)) ++ " " ++ (show $ (GV.toList $ slimFinal vCharacter)) ++ "\n" ++ (show minCost) ++ " " ++ (show maxCost)) 
        (minCost, maxCost)
        
    
    else if (thisCharType == WideSeq || thisCharType == AminoSeq) then
        let minMaxDiffList = if finalMethod == DirectOptimization then 
                                let uFinal = wideFinal uCharacter
                                    vFinal = wideFinal vCharacter
                                    newEdgeCharacter = M.getDOMedianCharInfo charInfo (uCharacter {widePrelim = uFinal}) (vCharacter {widePrelim = vFinal}) 
                                    (_, newU, newV) = wideGapped newEdgeCharacter
                                in
                                --trace ("GCD:\n" ++ (show m) ++ "\n" ++ (show (uFinal, newU)) ++ "\n" ++ (show (vFinal, newV)))
                                zipWith  (generalSequenceDiff thisMatrix (length thisMatrix)) (GV.toList $ GV.map (zero2GapWide gapCharWide) newU) (GV.toList $ GV.map (zero2GapWide gapCharWide) newV)
                             else GV.toList $ GV.zipWith (generalSequenceDiff thisMatrix (length thisMatrix))  (wideIAFinal uCharacter) (wideIAFinal vCharacter)
            (minDiff, maxDiff) = unzip minMaxDiffList
            minCost = thisWeight * (fromIntegral $ sum minDiff)
            maxCost = thisWeight * (fromIntegral $ sum maxDiff)
        in
        (minCost, maxCost)

    else if thisCharType == HugeSeq then
        let minMaxDiffList = if finalMethod == DirectOptimization then 
                                let uFinal = hugeFinal uCharacter
                                    vFinal = hugeFinal vCharacter
                                    newEdgeCharacter = M.getDOMedianCharInfo charInfo (uCharacter {hugePrelim = uFinal}) (vCharacter {hugePrelim = vFinal}) 
                                    (_, newU, newV) = hugeGapped newEdgeCharacter
                                in
                                --trace ("GCD:\n" ++ (show m) ++ "\n" ++ (show (uFinal, newU)) ++ "\n" ++ (show (vFinal, newV)))
                                zipWith  (generalSequenceDiff thisMatrix (length thisMatrix)) (GV.toList $ GV.map (zero2GapBV gapCharBV) newU) (GV.toList $ GV.map (zero2GapBV gapCharBV) newV)
                             else GV.toList $ GV.zipWith (generalSequenceDiff thisMatrix (length thisMatrix)) (hugeIAFinal uCharacter) (hugeIAFinal vCharacter)
            (minDiff, maxDiff) = unzip minMaxDiffList
            minCost = thisWeight * (fromIntegral $ sum minDiff)
            maxCost = thisWeight * (fromIntegral $ sum maxDiff)
        in
        (minCost, maxCost)

    else error ("Character type unimplemented : " ++ show thisCharType)

-- | zero2Gap converts a '0' or no bits set to gap (indel) value 
zero2Gap :: (FiniteBits a) => a -> a -> a
zero2Gap gapChar inVal = if popCount inVal == 0 then gapChar
                         else inVal

-- | zero2GapWide converts a '0' or no bits set to gap (indel) value 
zero2GapWide :: Word64 -> Word64 -> Word64
zero2GapWide gapChar inVal = if popCount inVal == 0 then gapChar
                         else inVal

-- | zero2GapBV converts a '0' or no bits set to gap (indel) value 
zero2GapBV :: BV.BitVector -> BV.BitVector -> BV.BitVector
zero2GapBV gapChar inVal = if popCount inVal == 0 then gapChar
                         else inVal

-- | maxIntervalDiff takes two ranges and gets the maximum difference between the two based on differences
-- in upp and lower ranges.
maxIntervalDiff :: (Int, Int)-> (Int, Int) -> Int 
maxIntervalDiff (a,b) (x,y) =
    let upper = (max b y) - (min b y)
        lower = (max a x) - (min a x)
    in
    max upper lower 

-- | minMaxMatrixDiff takes twovetors of states and calculates the minimum and maximum state differnce cost 
-- between the two
minMaxMatrixDiff :: S.Matrix Int -> V.Vector Int -> V.Vector Int -> (Int, Int)
minMaxMatrixDiff localCostMatrix uStatesV vStatesV =
    let statePairs = (V.toList uStatesV, V.toList vStatesV)
        cartesianPairs = cartProdPair statePairs
        costList = fmap (localCostMatrix S.!) cartesianPairs
    in 
    --trace (show cartesianPairs  ++ " " ++ show costList) 
    (minimum costList, maximum costList)
   


-- | generalSequenceDiff  takes two sequnce elemental bit types and retuns min and max integer 
-- cost differences using matrix values
generalSequenceDiff :: (FiniteBits a) => S.Matrix Int -> Int -> a -> a -> (Int, Int)
generalSequenceDiff thisMatrix numStates uState vState = 
    let uStateList = fmap snd $ filter ((== True).fst) $ zip (fmap (testBit uState) [0.. numStates - 1]) [0.. numStates - 1]
        vStateList = fmap snd $ filter ((== True).fst) $ zip (fmap (testBit vState) [0.. numStates - 1]) [0.. numStates - 1]    
        uvCombinations = cartProd uStateList vStateList
        costOfPairs = fmap (thisMatrix S.!) uvCombinations
    in
    -- trace ("GSD: " ++ (show uStateList) ++ " " ++ (show vStateList) ++ " min " ++ (show $ minimum costOfPairs) ++ " max " ++ (show $  maximum costOfPairs))
    (minimum costOfPairs, maximum costOfPairs)

