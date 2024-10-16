{- |
Module      :  PostOrderSoftWiredFunctionsNew.hs
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
module GraphOptimization.PostOrderSoftWiredFunctionsNew (
    postDecorateSoftWired,
    softWiredPostOrderTraceBack,
    createBlockResolutions,
    addNodeAndEdgeToResolutionData,
    updateRootCost,
    getOutDegree1VertexAndGraph,
    getOutDegree1VertexSoftWired,
    getOutDegree2VertexSoftWired,
    extractDisplayTrees,
    backPortBlockTreeNodesToCanonicalGraph,
) where

import Control.DeepSeq
import Control.Parallel.Strategies
import Data.BitVector.LittleEndian qualified as BV
import Data.Bits
import Data.List qualified as L
import Data.Maybe
import Data.Text.Lazy qualified as T
import Data.Vector qualified as V
import GeneralUtilities
import GraphOptimization.Medians qualified as M
import Graphs.GraphOperations qualified as GO
import PHANE.Evaluation
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import PHANE.Evaluation.Verbosity (Verbosity (..))
import Types.Types
import Utilities.LocalGraph qualified as LG
import Utilities.Utilities qualified as U


-- import Debug.Trace

{-Intial Postorder softwired pass.  All functions with 'New" appended-}

{- | postDecorateSoftWired begins at start index (usually root, but could be a subtree) and moves preorder till children are labelled
and then recurses to root postorder labelling vertices and edges as it goes
this for a single root
-}
postDecorateSoftWired
    ∷ GlobalSettings → Maybe (DecoratedGraph, LG.Node) → SimpleGraph → DecoratedGraph → V.Vector (V.Vector CharInfo) → LG.Node → LG.Node → PhyG PhylogeneticGraph
postDecorateSoftWired inGS incrementalGraph simpleGraph curDecGraph blockCharInfo rootIndex curNode =
    -- if node in current decorated graph then nothing to do and return it
    --   this because will hit node twice if network node
    if LG.gelem curNode curDecGraph
        then
            let nodeLabel = LG.lab curDecGraph curNode
            in  if isNothing nodeLabel
                    then error ("Null label for node " <> show curNode)
                    else pure (simpleGraph, subGraphCost (fromJust nodeLabel), curDecGraph, mempty, mempty, blockCharInfo)
        else -- node is not in decorated graph--ie has not been creted/optimized

        -- get postorder assignment of children
        -- checks for single child of node
        -- result is single graph after left and right child traversals
        -- trace ("PDSW making node " <> show curNode <> " in\n" <> (LG.prettify $ GO.convertDecoratedToSimpleGraph curDecGraph)) (

            let nodeChildren = LG.descendants simpleGraph curNode -- should be 1 or 2, not zero since all leaves already in graph
                leftChild = head nodeChildren
                rightChild = last nodeChildren -- will be same is first for out 1 (network) node
            in  do
                    leftChildTree ← postDecorateSoftWired inGS incrementalGraph simpleGraph curDecGraph blockCharInfo rootIndex leftChild
                    rightLeftChildTree ←
                        if length nodeChildren == 2
                            then postDecorateSoftWired inGS incrementalGraph simpleGraph (thd6 leftChildTree) blockCharInfo rootIndex rightChild
                            else pure leftChildTree

                    -- Checks on children
                    if length nodeChildren > 2
                        then error ("Graph not dichotomous in postDecorateSoftWired node " <> show curNode <> "\n" <> LG.prettyIndices simpleGraph)
                        else
                            if null nodeChildren
                                then error ("Leaf not in graph in postDecorateSoftWired node " <> show curNode <> "\n" <> LG.prettyIndices simpleGraph)
                                else -- after recursing to any children can optimize current node

                                -- make node from child block resolutions
                                -- child resolution made earler in post order pass
                                --  sub tree is updated graph from children--ie not iuncluding current node

                                    let newSubTree = thd6 rightLeftChildTree
                                    in  -- single child of node--network node
                                        if length nodeChildren == 1
                                            then -- use left child for out degree = 1 nmodes, right should be "Nothing"

                                                let (newGraph, _, _, _, _) = getOutDegree1VertexAndGraph curNode (fromJust $ LG.lab newSubTree leftChild) simpleGraph nodeChildren newSubTree
                                                in  -- display graphs and character block are not done yet since will need traceback to get preliminary states
                                                    --  graph ahas all soft wired decorations
                                                    pure (simpleGraph, 0, newGraph, mempty, mempty, blockCharInfo)
                                            else -- 2 children--tree node

                                            -- trace ("Outdegree 2: " <> (show curNode) <> " " <> (show $ GO.getNodeType simpleGraph curNode) <> " Children: " <> (show nodeChildren)) (
                                            -- need to create  resolutions and add to existing sets

                                                let -- this ensures that left/right choices are based on leaf BV for consistency and label invariance
                                                    -- larger bitvector is Right, smaller or equal Left
                                                    ((leftChild', leftChildLabel), (rightChild', rightChildLabel)) =
                                                        U.leftRightChildLabelBVNode
                                                            ((leftChild, fromJust $ LG.lab newSubTree leftChild), (rightChild, fromJust $ LG.lab newSubTree rightChild))

                                                    -- create resolution caches for blocks
                                                    leftChildNodeType = GO.getNodeType simpleGraph leftChild' -- nodeType leftChildLabel
                                                    rightChildNodeType = GO.getNodeType simpleGraph rightChild' -- nodeType rightChildLabel
                                                    leftEdgeType
                                                        | leftChildNodeType == NetworkNode = NetworkEdge
                                                        | leftChildNodeType == LeafNode = PendantEdge
                                                        | otherwise = TreeEdge
                                                    rightEdgeType
                                                        | rightChildNodeType == NetworkNode = NetworkEdge
                                                        | rightChildNodeType == LeafNode = PendantEdge
                                                        | otherwise = TreeEdge

                                                    edgeLable =
                                                        EdgeInfo
                                                            { minLength = 0.0
                                                            , maxLength = 0.0
                                                            , midRangeLength = 0.0
                                                            , edgeType = TreeEdge
                                                            }
                                                in  do
                                                        resolutionBlockVL ←
                                                            mapM
                                                                ( createBlockResolutions'
                                                                    inGS
                                                                    (compressResolutions inGS)
                                                                    curNode
                                                                    leftChild'
                                                                    rightChild'
                                                                    leftChildNodeType
                                                                    rightChildNodeType
                                                                    (GO.getNodeType simpleGraph curNode)
                                                                )
                                                                (V.zip3 (vertexResolutionData leftChildLabel) (vertexResolutionData rightChildLabel) blockCharInfo)

                                                        -- create canonical Decorated Graph vertex
                                                        -- 0 cost becasue can't know cosrt until hit root and get best valid resolutions
                                                        let newVertexLabel =
                                                                VertexInfo
                                                                    { index = curNode
                                                                    , bvLabel = bvLabel leftChildLabel .|. bvLabel rightChildLabel
                                                                    , parents = V.fromList $ LG.parents simpleGraph curNode
                                                                    , children = V.fromList nodeChildren
                                                                    , nodeType = GO.getNodeType simpleGraph curNode -- TreeNode
                                                                    , vertName = T.pack $ "HTU" <> show curNode
                                                                    , vertData = mempty -- empty because of resolution data
                                                                    , vertexResolutionData = resolutionBlockVL
                                                                    , vertexCost = 0.0 -- newCost
                                                                    , subGraphCost = 0.0 -- (subGraphCost leftChildLabel) + (subGraphCost rightChildLabel) + newCost
                                                                    }
                                                        let leftEdge = (curNode, leftChild', edgeLable{edgeType = leftEdgeType})
                                                        let rightEdge = (curNode, rightChild', edgeLable{edgeType = rightEdgeType})
                                                        let newGraph = LG.insEdges [leftEdge, rightEdge] $ LG.insNode (curNode, newVertexLabel) newSubTree

                                                        let (displayGraphVL, lDisplayCost) =
                                                                if curNode == rootIndex
                                                                    then
                                                                        let (displayG, displayCost') = extractDisplayTrees (Just curNode) True resolutionBlockVL
                                                                        in  (displayG, displayCost')
                                                                    else -- (fmap (fmap LG.removeDuplicateEdges) displayG, displayCost')
                                                                        (mempty, 0.0)

                                                        pure (simpleGraph, lDisplayCost, newGraph, displayGraphVL, mempty, blockCharInfo)


-- | getOutDegree1VertexAndGraph makes parent node from single child for soft-wired resolutions
getOutDegree1VertexAndGraph
    ∷ (Show a, Show b)
    ⇒ LG.Node
    → VertexInfo
    → LG.Gr a b
    → [LG.Node]
    → DecoratedGraph
    → (DecoratedGraph, Bool, VertexInfo, VertexCost, V.Vector [DecoratedGraph])
getOutDegree1VertexAndGraph curNode childLabel simpleGraph nodeChildren subTree =
    -- trace ("In out=1: " <> (show curNode)) (
    let childResolutionData = vertexResolutionData childLabel

        curNodeResolutionData = addNodeAndEdgeToResolutionData newDisplayNode newLEdge childResolutionData

        newEdgeLabel =
            EdgeInfo
                { minLength = 0.0
                , maxLength = 0.0
                , midRangeLength = 0.0
                , edgeType = TreeEdge
                }
        newMinVertex =
            VertexInfo
                { index = curNode
                , bvLabel = bvLabel childLabel
                , parents = V.fromList $ LG.parents simpleGraph curNode
                , children = V.fromList nodeChildren
                , nodeType = GO.getNodeType simpleGraph curNode
                , vertName = T.pack $ "HTU" <> show curNode
                , vertData = mempty
                , vertexResolutionData = mempty
                , vertexCost = 0.0
                , subGraphCost = 0.0
                }

        newVertex =
            VertexInfo
                { index = curNode
                , bvLabel = bvLabel childLabel
                , parents = V.fromList $ LG.parents simpleGraph curNode
                , children = V.fromList nodeChildren
                , nodeType = GO.getNodeType simpleGraph curNode -- NetworkNode
                , vertName = T.pack $ "HTU" <> show curNode
                , vertData = mempty
                , vertexResolutionData = curNodeResolutionData
                , vertexCost = 0.0
                , subGraphCost = subGraphCost childLabel
                }

        newLEdge = (curNode, index childLabel, newEdgeLabel)
        newLNode = (curNode, newVertex)
        newDisplayNode = (curNode, newMinVertex)
        newGraph = LG.insEdge newLEdge $ LG.insNode newLNode subTree

        -- Root node should be out degree 2 so this should not happen in general--but could during some
        -- graph rearangemnts in fusing and swapping
        (displayGraphVL, lDisplayCost) =
            if nodeType newVertex == RootNode
                then extractDisplayTrees (Just curNode) True (vertexResolutionData childLabel)
                else (mempty, 0.0)
    in  -- trace ("NV1: " <> show newVertex)
        -- trace ("GOD1VG: " <> (show $ LG.toEdge newLEdge) <> " has edges " <> (show $ LG.hasEdge subTree $  LG.toEdge newLEdge) <> "Resolutions " <> (show $ fmap (fmap U.hasResolutionDuplicateEdges) curNodeResolutionData))
        -- trace ("PDSW-1:" <> (show $ bvLabel newVertex))
        (newGraph, nodeType newVertex == RootNode, newVertex, lDisplayCost, displayGraphVL)


-- (newGraph, False, newVertex, 0.0, mempty)
-- )

-- | getOutDegree1VertexSoftWired returns new vertex only from single child for soft-wired resolutions
getOutDegree1VertexSoftWired
    ∷ (Show a, Show b)
    ⇒ LG.Node
    → VertexInfo
    → LG.Gr a b
    → [LG.Node]
    → VertexInfo
getOutDegree1VertexSoftWired curNode childLabel simpleGraph nodeChildren =
    -- trace ("In out=1: " <> (show curNode)) (
    let childResolutionData = vertexResolutionData childLabel

        newEdgeLabel =
            EdgeInfo
                { minLength = 0.0
                , maxLength = 0.0
                , midRangeLength = 0.0
                , edgeType = TreeEdge
                }
        newMinVertex =
            VertexInfo
                { index = curNode
                , bvLabel = bvLabel childLabel
                , parents = V.fromList $ LG.parents simpleGraph curNode
                , children = V.fromList nodeChildren
                , nodeType = NetworkNode
                , vertName = T.pack $ "HTU" <> show curNode
                , vertData = mempty
                , vertexResolutionData = mempty
                , vertexCost = 0.0
                , subGraphCost = 0.0
                }

        newDisplayNode = (curNode, newMinVertex)
        newLEdge = (curNode, index childLabel, newEdgeLabel)

        curNodeResolutionData = addNodeAndEdgeToResolutionData newDisplayNode newLEdge childResolutionData

        newVertexLabel =
            VertexInfo
                { index = curNode
                , bvLabel = bvLabel childLabel
                , parents = V.fromList $ LG.parents simpleGraph curNode
                , children = V.fromList nodeChildren
                , nodeType = NetworkNode
                , vertName = T.pack $ "HTU" <> show curNode
                , vertData = mempty
                , vertexResolutionData = curNodeResolutionData
                , vertexCost = 0.0
                , subGraphCost = subGraphCost childLabel
                }
    in  newVertexLabel


{- | getOutDegree2VertexSoftWired returns new vertex only from two child nodes for soft-wired resolutions
used in Net Add Delete heuristics
-}
getOutDegree2VertexSoftWired
    ∷ GlobalSettings
    → V.Vector (V.Vector CharInfo)
    → LG.Node
    → LG.LNode VertexInfo
    → LG.LNode VertexInfo
    → DecoratedGraph
    → PhyG VertexInfo
getOutDegree2VertexSoftWired inGS charInfoVectVect curNodeIndex leftChild@(leftChildIndex, _) rightChild@(rightChildIndex, _) inGraph =
    let -- this ensures that left/right choices are based on leaf BV for consistency and label invariance
        -- larger bitvector is Right, smaller or equal Left
        ((leftChild', leftChildLabel'), (rightChild', rightChildLabel')) = U.leftRightChildLabelBVNode (leftChild, rightChild)

        -- create resolution caches for blocks
        leftChildNodeType = nodeType leftChildLabel'
        rightChildNodeType = nodeType rightChildLabel'
    in  do
            -- TODO PArallelize? its parallel in lower call
            resolutionBlockVL ←
                mapM
                    ( createBlockResolutions'
                        inGS
                        (compressResolutions inGS)
                        curNodeIndex
                        leftChild'
                        rightChild'
                        leftChildNodeType
                        rightChildNodeType
                        TreeNode
                    )
                    (V.zip3 (vertexResolutionData leftChildLabel') (vertexResolutionData rightChildLabel') charInfoVectVect)

            -- create canonical Decorated Graph vertex
            -- 0 cost becasue can't know cosrt until hit root and get best valid resolutions
            let newVertexLabel =
                    VertexInfo
                        { index = curNodeIndex
                        , bvLabel = bvLabel leftChildLabel' .|. bvLabel rightChildLabel'
                        , parents = V.fromList $ LG.parents inGraph curNodeIndex
                        , children = V.fromList [leftChildIndex, rightChildIndex]
                        , nodeType = TreeNode
                        , vertName = T.pack $ "HTU" <> show curNodeIndex
                        , vertData = mempty -- empty because of resolution data
                        , vertexResolutionData = resolutionBlockVL
                        , vertexCost = 0.0 -- newCost
                        , subGraphCost = 0.0 -- (subGraphCost leftChildLabel) + (subGraphCost rightChildLabel) + newCost
                        }

            pure newVertexLabel


{- | addNodeAndEdgeToResolutionData adds  node and edge to resolution data in outdegree = 1 nodes
straight copy would not add this node or edge to subtree in resolutions
-}
addNodeAndEdgeToResolutionData
    ∷ LG.LNode VertexInfo → LG.LEdge EdgeInfo → V.Vector ResolutionBlockData → V.Vector ResolutionBlockData
addNodeAndEdgeToResolutionData newNode newEdge = fmap (addNodeEdgeToResolutionBlock newNode newEdge True)


-- | addNodeEdgeToResolutionBlock adds  node and edge to resolutoin block data
addNodeEdgeToResolutionBlock ∷ LG.LNode VertexInfo → LG.LEdge EdgeInfo → Bool → ResolutionBlockData → ResolutionBlockData
addNodeEdgeToResolutionBlock newNode newEdge isIn1Out1Node inResBlockData =
    V.zipWith
        (addNodeEdgeToResolutionList newNode newEdge isIn1Out1Node)
        inResBlockData
        (V.fromList [0 .. V.length inResBlockData - 1])


{- | addNodeEdgeToResolutionList adds  node and edge to single subGraph in ResolutionData
adds resolution pairs to be equal to the child straight one-for-one correpondance
although only a single child-both indieces are set to resolution index since singels can be added to
paired resolitions if one or other is a network node
-}
addNodeEdgeToResolutionList ∷ LG.LNode VertexInfo → LG.LEdge EdgeInfo → Bool → ResolutionData → Int → ResolutionData
addNodeEdgeToResolutionList newNode newEdge _ inResData resolutionIndex =
    let (inNodeList, inEdgeList) = displaySubGraph inResData

        -- childResolutionIndexPairList = childResolutions inResData
        newNodeList = newNode : inNodeList

        -- this check for redundant edges in resolution cash from combinations
        newEdgeList =
            if newEdge `notElem` inEdgeList
                then newEdge : inEdgeList
                else -- trace "Should not happen: Extra edge in addNodeEdgeToResolutionListNew"
                    inEdgeList
        newFirstData =
            inResData
                { displaySubGraph = (newNodeList, newEdgeList)
                , -- both set because can be a display node left right added to 2 child resolutoins
                  childResolutionIndices = (Just resolutionIndex, Just resolutionIndex)
                }
    in  -- trace ("ANETRL:" <> (show $ Just resolutionIndex))
        newFirstData


-- | createBlockResolutions' is a wrapper around createBlockResolutions
createBlockResolutions'
    ∷ GlobalSettings
    → Bool
    → LG.Node
    → Int
    → Int
    → NodeType
    → NodeType
    → NodeType
    → (ResolutionBlockData, ResolutionBlockData, V.Vector CharInfo)
    → PhyG ResolutionBlockData
createBlockResolutions' inGS compressResolutions curNode leftIndex rightIndex leftChildNodeType rightChildNodeType curNodeNodeType (leftChild, rightChild, charInfoV) =
    createBlockResolutions
        inGS
        compressResolutions
        curNode
        leftIndex
        rightIndex
        leftChildNodeType
        rightChildNodeType
        curNodeNodeType
        leftChild
        rightChild
        charInfoV


{- | createBlockResolutions takes left and right child resolution data for a block (same display tree)
and generates node resolution data
-}
createBlockResolutions
    ∷ GlobalSettings
    → Bool
    → LG.Node
    → Int
    → Int
    → NodeType
    → NodeType
    → NodeType
    → ResolutionBlockData
    → ResolutionBlockData
    → V.Vector CharInfo
    → PhyG ResolutionBlockData
createBlockResolutions
    inGS
    compressResolutions
    curNode
    leftIndex
    rightIndex
    leftChildNodeType
    rightChildNodeType
    curNodeNodeType
    leftChild
    rightChild
    charInfoV
        | null leftChild && null rightChild = pure mempty
        | null leftChild = pure rightChild
        | null rightChild = pure leftChild
        | otherwise =
            -- trace ("CBR:" <> (show (leftIndex, leftChildNodeType, rightIndex, rightChildNodeType)) <> (show $fmap BV.toBits $ fmap displayBVLabel leftChild) <> " and " <>  (show $fmap BV.toBits $ fmap displayBVLabel rightChild)) (
            -- trace ("CNR: " <> (show (length leftChild, length rightChild))) (
            let childResolutionPairs = cartProd (V.toList leftChild) (V.toList rightChild)
                -- need to keep these indices correct (hence reverse in checkLeafOverlap ) for traceback and compress
                childResolutionIndices = cartProd [0 .. (length leftChild - 1)] [0 .. (length rightChild - 1)]
                validPairs = concatMap checkLeafOverlap (zip childResolutionPairs childResolutionIndices)

                -- need to add in node and edge to left and right
                edgeLable =
                    EdgeInfo
                        { minLength = 0.0
                        , maxLength = 0.0
                        , midRangeLength = 0.0
                        , edgeType = TreeEdge
                        }
                newMinVertex =
                    VertexInfo
                        { index = curNode
                        , bvLabel = BV.fromBits [False]
                        , parents = mempty
                        , children = mempty
                        , nodeType = curNodeNodeType
                        , vertName = T.pack $ "HTU" <> show curNode
                        , vertData = mempty
                        , vertexResolutionData = mempty
                        , vertexCost = 0.0
                        , subGraphCost = 0.0
                        }

                newNode = (curNode, newMinVertex)

                addLeft =
                    if leftChildNodeType == NetworkNode
                        then
                            let newEdge = (curNode, rightIndex, edgeLable)
                                newRightChildBlockResolutionData = addNodeEdgeToResolutionBlock newNode newEdge False rightChild
                            in  -- trace ("ANEL:" <> (show $ (curNode, rightIndex)))
                                newRightChildBlockResolutionData
                        else -- trace ("ANEL-Nothing")
                            mempty

                addRight =
                    if rightChildNodeType == NetworkNode
                        then
                            let newEdge = (curNode, leftIndex, edgeLable)
                                newLeftChildBlockResolutionData = addNodeEdgeToResolutionBlock newNode newEdge False leftChild
                            in  -- trace ("ANER:" <> (show $ (curNode, leftIndex)))
                                newLeftChildBlockResolutionData
                        else -- trace ("ANER-Nothing")
                            mempty

                resolutionAction ∷ ((ResolutionData, ResolutionData), (Int, Int)) → PhyG ResolutionData
                resolutionAction = createNewResolution inGS curNode leftIndex rightIndex leftChildNodeType rightChildNodeType charInfoV
            in  do
                    resolutionPar ← getParallelChunkTraverse
                    newResolutionList <- resolutionPar resolutionAction validPairs
                    -- newResolutionList = PU.seqParMap PU.myStrategy (createNewResolution curNode leftIndex rightIndex leftChildNodeType rightChildNodeType charInfoV) validPairs

                    -- trace ("CNR: " <> (show (length leftChild, length rightChild))) ( --  <> "\n" <> (show childResolutionIndices) <> "\n" <> (show $ fmap snd validPairs)) (
                    if compressResolutions
                        then pure $ compressBlockResolution (newResolutionList <> V.toList addLeft <> V.toList addRight)
                        else pure $ V.fromList newResolutionList <> (addLeft <> addRight)


{- | compressBlockResolution 'compresses' resolutions of a block by taking only the first of resolutions with the
same set of leaves (via bitvector) and lowest cost
can speed up graph diagnosis, but at the cost of potentially loosing resolutions whihc would be better later
(ie closer to root)
-}
compressBlockResolution ∷ [ResolutionData] → V.Vector ResolutionData
compressBlockResolution inResList =
    if null inResList
        then V.empty
        else
            let -- group by bitvectors of subtree
                resLL = L.groupBy compareBVLabel inResList
                minCostVV = fmap getMinCostList resLL
            in  V.fromList minCostVV
    where
        compareBVLabel a b = displayBVLabel a == displayBVLabel b


-- | getMinCostList takes a list resolutions and returns first lowest cost resolution
getMinCostList ∷ [ResolutionData] → ResolutionData
getMinCostList inList =
    if null inList
        then error "Empty resolution list in getMinCostList"
        else
            let minResCost = minimum $ fmap displayCost inList
            in  head $ filter ((== minResCost) . displayCost) inList


{- | createNewResolution takes a pair of resolutions and creates the median resolution
need to watch let/right (based on BV) for preorder stuff
-}
createNewResolution
    ∷ GlobalSettings
    → LG.Node
    → Int
    → Int
    → NodeType
    → NodeType
    → V.Vector CharInfo
    → ((ResolutionData, ResolutionData), (Int, Int))
    → PhyG ResolutionData
createNewResolution inGS curNode leftIndex rightIndex leftChildNodeType rightChildNodeType charInfoV ((leftRes, rightRes), (leftResIndex, rightResIndex)) =
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

        edgeLable =
            EdgeInfo
                { minLength = 0.0
                , maxLength = 0.0
                , midRangeLength = 0.0
                , edgeType = TreeEdge
                }

        leftEdge = (curNode, leftIndex, edgeLable{edgeType = leftEdgeType})
        rightEdge = (curNode, rightIndex, edgeLable{edgeType = rightEdgeType})
        leftChildTree = displaySubGraph leftRes
        rightChildTree = displaySubGraph rightRes

        -- Data fields empty for display tree data--not needed and multiple copies of everything
        newNodeLabel =
            VertexInfo
                { index = curNode
                , bvLabel = resBV
                , parents = V.empty
                , children = V.fromList [leftIndex, rightIndex]
                , nodeType = TreeNode
                , vertName = T.pack $ "HTU" <> show curNode
                , vertData = mempty
                , vertexResolutionData = mempty
                , vertexCost = 0.0
                , subGraphCost = 0.0
                }

        newNode = (curNode, newNodeLabel)

        -- this check for redundant edges in resoluton cash from combinations
        -- resolutionEdgeList = leftEdge : (rightEdge: (snd leftChildTree <> snd rightChildTree))
        -- is this required?
        existingEdges = snd leftChildTree <> snd rightChildTree
        resolutionEdgeList
            | (leftEdge `notElem` existingEdges) && (rightEdge `notElem` existingEdges) = leftEdge : (rightEdge : existingEdges)
            | (leftEdge `elem` existingEdges) && (rightEdge `elem` existingEdges) = existingEdges
            | (leftEdge `notElem` existingEdges) = leftEdge : existingEdges
            | otherwise = rightEdge : existingEdges

        resolutionNodeList = newNode : (fst leftChildTree <> fst rightChildTree)

        -- Make the data and cost for the resolution
        -- No chnage cost adjustment is True here if PMDL/SI
        leftBlockLength = V.length $ displayData leftRes
        rightBlockLength = V.length $ displayData rightRes
    in do
        resolutionMedianCostV <-
            if (leftBlockLength == 0) then pure $ V.zip (displayData rightRes) (V.replicate rightBlockLength 0)
            else if (rightBlockLength == 0) then pure $ V.zip (displayData leftRes) (V.replicate leftBlockLength 0)
            else M.median2M (U.needTwoEdgeNoCostAdjust inGS True) (displayData leftRes) (displayData rightRes) charInfoV
        let (resolutionMedianV, resolutionCostV) = V.unzip resolutionMedianCostV
        let thisResolutionCost = V.sum resolutionCostV
        let displaySubTreeCost = displayCost leftRes + displayCost rightRes + thisResolutionCost
   
        pure ResolutionData
            { displaySubGraph = (resolutionNodeList, resolutionEdgeList)
            , displayBVLabel = resBV
            , displayData = resolutionMedianV
            , childResolutionIndices = (Just leftResIndex, Just rightResIndex)
            , resolutionCost = thisResolutionCost
            , displayCost = displaySubTreeCost
            }


{- | extractDisplayTrees  takes resolutions and pulls out best cost (head for now) need to change type for multiple best
option for filter based on pop-count for root cost and complete display tree check
-}
extractDisplayTrees ∷ Maybe Int → Bool → V.Vector ResolutionBlockData → (V.Vector [DecoratedGraph], VertexCost)
extractDisplayTrees startVertex checkPopCount inRBDV =
    if V.null inRBDV
        then (V.empty, 0.0)
        else
            let (bestBlockDisplayResolutionList, costVect, _) = V.unzip3 $ fmap (getBestResolutionList startVertex checkPopCount) inRBDV
            in  -- trace ("EDT: " <> (show (V.length bestBlockDisplayResolutionList, fmap length bestBlockDisplayResolutionList, V.sum costVect)))
                (bestBlockDisplayResolutionList, V.sum costVect)


{- | checkLeafOverlap  takes a left right resolution pair list and checks if
there is leaf overlap via comparing displayBVLabel if & = 0 then no
overlap, and adds to resulting list--reverses order--sholdn't matter
-}
checkLeafOverlap ∷ ((ResolutionData, ResolutionData), (Int, Int)) → [((ResolutionData, ResolutionData), (Int, Int))]
checkLeafOverlap inPair@((leftRes, rightRes), (_, _)) =
    let leftBV = displayBVLabel leftRes
        rightBV = displayBVLabel rightRes
    in  ([inPair | BV.isZeroVector (leftBV .&. rightBV)])


{- | getBestResolutionList  takes ResolutionBlockData and retuns a list of the best valid (ie all leaves in subtree) display trees
for that block with cost and resolujtion data in triple-- if checkPopCount is True--otherwise all display trees of any cost and contitution
startVertex for a component-- to allow for not every leaf being in componnet but still getting softwired cost
-}
getBestResolutionList ∷ Maybe Int → Bool → ResolutionBlockData → ([DecoratedGraph], VertexCost, ResolutionBlockData)
getBestResolutionList startVertex checkPopCount inRDList =
    -- trace ("GBRL: " <> (show $ V.length inRDList)) (
    if null inRDList
        then error "Null resolution list"
        else
            let displayTreeList = fmap displaySubGraph inRDList
                displayCostList = fmap displayCost inRDList
                displayPopList = fmap (complement . displayBVLabel) inRDList
            in  if not checkPopCount
                    then
                        let minCost = minimum displayCostList
                            displayCostTripleList = V.zip3 displayTreeList displayCostList inRDList
                            (bestDisplayList, _, bestResList) = V.unzip3 $ V.filter ((== minCost) . snd3) displayCostTripleList
                        in  (fmap LG.mkGraphPair (V.toList bestDisplayList), minCost, bestResList)
                    else
                        let minPopCount = minimum $ fmap popCount displayPopList -- max since complemented above
                            displayBVList = V.zip4 displayTreeList displayCostList displayPopList inRDList

                            -- must have all leaves if startvzertex == Nothing, component maximum otherwise
                            -- this for getting cost of component of a softwired network
                            validDisplayList =
                                if isNothing startVertex
                                    then V.filter (BV.isZeroVector . thd4) displayBVList
                                    else V.filter ((== minPopCount) . (popCount . thd4)) displayBVList
                            validMinCost = V.minimum $ fmap snd4 validDisplayList
                            (bestDisplayList, _, _, bestResList) = V.unzip4 $ V.filter ((== validMinCost) . snd4) validDisplayList
                        in  -- trace ("Valid display list number:" <> (show $ length validDisplayList)) (
                            if V.null validDisplayList
                                then
                                    error
                                        ( "Null root validDisplayList in getBestResolutionList "
                                            <> show (startVertex, inRDList)
                                            <> " This can be caused if the graphType not set correctly."
                                        )
                                else
                                    let lDisplayTreeList = fmap LG.mkGraphPair (V.toList bestDisplayList)

                                        -- update root cost of display trees for use later (e.g. net penalties, outputting display forrests)
                                        lDisplayTreeList' = fmap (updateRootCost validMinCost) lDisplayTreeList
                                    in  -- trace ("GBRL: " <> (show (length lDisplayTreeList', validMinCost)))
                                        (lDisplayTreeList', validMinCost, bestResList)


-- )

{- | updateRootCost  updates the subGraphCost of the root node(s) with input value
 node is created, so original is deleted,  added, and original edges added back
since deleted when node is
assumes its a tree wiht a single root
-}
updateRootCost ∷ VertexCost → DecoratedGraph → DecoratedGraph
updateRootCost newRootCost inGraph =
    let (rootIndex, rootLabel) = head $ LG.getRoots inGraph
        rootEdges = LG.out inGraph rootIndex
        newRootLabel = rootLabel{subGraphCost = newRootCost}
    in  -- trace ("DCC: " <> (show newRootCost))
        LG.insEdges rootEdges $ LG.insNode (rootIndex, newRootLabel) $ LG.delNode rootIndex inGraph


{-
   Traceback code for after intial Postorder softwired pass.  All functions with 'New" appended
-}

{- | softWiredPostOrderTraceBack  takes resolution data and assigns correct resolution median
from vertexResolutionData to preliminary data assignments in vertData.
Proceeds via typical pre-order pass over display tree for each block
using the indices of left and right (if present) of resolutions
first gets root assignment from resolution data and then each block is traversed given its block display tree
-}
softWiredPostOrderTraceBack ∷ Int → PhylogeneticGraph → PhyG PhylogeneticGraph
softWiredPostOrderTraceBack rootIndex inGraph@(inSimpleGraph, b, canonicalGraph, _, _, f) =
    if LG.isEmpty canonicalGraph
        then pure emptyPhylogeneticGraph
        else -- this condition can arise due to strictness in graph evaluation in parallel

            if length (LG.descendants canonicalGraph rootIndex) /= 2
                then pure emptyPhylogeneticGraph
                else -- error ("Root node has improper number of children: " <> show (LG.descendants canonicalGraph rootIndex) <>"\n" <> (LG.prettyIndices canonicalGraph))

                    let -- extract display trees and bloxck char trees from PhylogeneticGraph
                        -- block character trees do not exist yet
                        displayTreeV = (head <$> fth6 inGraph)

                        -- get root node resolution data from canonical Graph created in postorder
                        rootLabel = fromJust $ LG.lab canonicalGraph rootIndex
                        rootResData = vertexResolutionData rootLabel

                        -- traceback for each block based on its display tree, updating trees as it goes, left descendent then right
                        -- at this stage all character trees will have same root descendents sionce all rooted from outgropu postorder traversal
                        -- later (after reroot pass) this will not be the case since each charcater may have a unique traversal root/edge
                        leftChild = head $ LG.descendants canonicalGraph rootIndex
                        rightChild = last $ LG.descendants canonicalGraph rootIndex

                        -- get left right from BV as in postorder
                        ((leftChild', _), (rightChild', _)) =
                            U.leftRightChildLabelBVNode
                                ((leftChild, fromJust $ LG.lab canonicalGraph leftChild), (rightChild, fromJust $ LG.lab canonicalGraph rightChild))

                        getBestAction ∷ ResolutionBlockData → ([DecoratedGraph], VertexCost, ResolutionBlockData)
                        getBestAction = getBestResolutionList (Just rootIndex) True

                        updateAction ∷ (ResolutionData, DecoratedGraph) → (DecoratedGraph, V.Vector DecoratedGraph)
                        updateAction = updateRootBlockTrees rootIndex

                        traceBackAction
                            ∷ LG.Node → (DecoratedGraph, V.Vector DecoratedGraph, Maybe Int, Int) → (DecoratedGraph, V.Vector DecoratedGraph)
                        traceBackAction = traceBackBlock canonicalGraph
                    in  do
                            -- let (rootNodes, leafNode, treeNodes,networkNodes) = LG.splitVertexList inSimpleGraph
                            -- logWith LogInfo ("SWPOT: " <> (show (length rootNodes, length leafNode, length treeNodes, length networkNodes)))

                            -- extract (first) best resolution for each block--there can be more than one for each, but only use the first for
                            -- traceback, preliminary and final assignment etc--part of the heuristic
                            resolutionPar ← getParallelChunkMap
                            let resolutionResult = resolutionPar getBestAction $ V.toList rootResData
                            let (_, _, rootDisplayBlockCharResolutionV) = unzip3 resolutionResult
                            -- \$ PU.seqParMap PU.myStrategy (getBestResolutionList (Just rootIndex) True) rootResData
                            let firstOfEachRootRes = fmap V.head rootDisplayBlockCharResolutionV

                            -- get preliminary character data for blocks
                            -- these should be ok wihtout left right check since were creted with that check on post order
                            let (leftIndexList, rightIndexList) = V.unzip $ fmap childResolutionIndices $ V.fromList firstOfEachRootRes

                            -- update root vertex info for display and character trees for each block
                            -- this includes preliminary data and other fields
                            updatePar ← getParallelChunkMap
                            let updateResult = updatePar updateAction (zip firstOfEachRootRes (V.toList displayTreeV))

                            let (rootUpdatedDisplayTreeV, rootUpdatedCharTreeVV) = unzip updateResult
                            -- \$ PU.seqParMap PU.myStrategy (updateRootBlockTrees rootIndex) (V.zip  (V.fromList firstOfEachRootRes) displayTreeV)

                            traceLeftPar ← getParallelChunkMap
                            let leftResult =
                                    traceLeftPar
                                        (traceBackAction leftChild')
                                        (L.zip4 rootUpdatedDisplayTreeV rootUpdatedCharTreeVV (V.toList leftIndexList) ([0 .. (length rootUpdatedDisplayTreeV - 1)]))
                            let (traceBackDisplayTreeVLeft, traceBackCharTreeVVLeft) = unzip leftResult
                            -- \$ PU.seqParMap PU.myStrategy (traceBackBlock canonicalGraph leftChild') (V.zip4 rootUpdatedDisplayTreeV rootUpdatedCharTreeVV leftIndexList (V.fromList [0..(V.length rootUpdatedDisplayTreeV - 1)]))

                            traceRightPar ← getParallelChunkMap
                            let rightResult =
                                    traceRightPar
                                        (traceBackAction rightChild')
                                        ( L.zip4 traceBackDisplayTreeVLeft traceBackCharTreeVVLeft (V.toList rightIndexList) ([0 .. (length rootUpdatedDisplayTreeV - 1)])
                                        )
                            let (traceBackDisplayTreeV, traceBackCharTreeVV) = unzip rightResult
                            -- \$ PU.seqParMap PU.myStrategy (traceBackBlock canonicalGraph rightChild') (V.zip4 traceBackDisplayTreeVLeft traceBackCharTreeVVLeft rightIndexList (V.fromList [0..(V.length rootUpdatedDisplayTreeV - 1)]))

                            let newCanonicalGraph = backPortBlockTreeNodesToCanonicalGraph canonicalGraph (V.fromList traceBackDisplayTreeV)

                            -- this is a hack due to missing nodes in some character trees--perhpas issue with postorder resolutions?
                            if LG.isEmpty newCanonicalGraph
                                then pure emptyPhylogeneticGraph
                                else
                                    pure $
                                        (inSimpleGraph, b, newCanonicalGraph, fmap (: []) (V.fromList traceBackDisplayTreeV), (V.fromList traceBackCharTreeVV), f)


{- | traceBackBlock performs softwired traceback on block data returns updated display and character trees
the block index specifies which resolution list from the canonical tree at node nodeIndex
the resoliution index is the resolution element that was used to create the parent's state
need to account for out degree 2 and 1 in recursive calls
-}
traceBackBlock
    ∷ DecoratedGraph
    → LG.Node
    → (DecoratedGraph, V.Vector DecoratedGraph, Maybe Int, Int)
    → (DecoratedGraph, V.Vector DecoratedGraph)
traceBackBlock canonicalGraph nodeIndex (displayTree, charTreeV, resolutionIndex, blockIndex) =
    -- trace ("TBB: " <> "Node " <> (show nodeIndex) <> " Block " <> (show blockIndex) <> " Resolution " <> (show (resolutionIndex, LG.descendants displayTree nodeIndex))
    -- <> " nrd: " <> (show (length (vertexResolutionData (fromJust $ LG.lab canonicalGraph nodeIndex)), fmap length (vertexResolutionData (fromJust $ LG.lab canonicalGraph nodeIndex)))) ) (
    -- <> "\n" <> (show (vertexResolutionData (fromJust $ LG.lab canonicalGraph nodeIndex)))) (
    if LG.isEmpty displayTree || V.null charTreeV
        then error "Null data in traceBackBlock"
        else
            let -- get block resolution data from canonical graph
                nodeCanonicalLabel = fromJust $ LG.lab canonicalGraph nodeIndex
                nodeResolutionData = vertexResolutionData nodeCanonicalLabel

                -- hack checking for too high res index here
                newIndex =
                    if (fromJust resolutionIndex) < length (nodeResolutionData V.! blockIndex)
                        then fromJust resolutionIndex
                        else -- traceNoLF ("(" <> (show (fromJust resolutionIndex, length (nodeResolutionData V.! blockIndex))) <> ")")
                            (length (nodeResolutionData V.! blockIndex)) - 1

                -- blockResolutionData = (nodeResolutionData V.! blockIndex) V.! fromJust resolutionIndex
                blockResolutionData = (nodeResolutionData V.! blockIndex) V.! newIndex

                -- update display tree and character tree nodes
                (newDisplayTree, newCharTreeV) = updateNodeBlockTrees nodeIndex (blockResolutionData, displayTree, charTreeV)

                -- get left and right display resolution indices
                (leftResIndex, rightResIndex) = childResolutionIndices blockResolutionData

                -- get left and right node children, won't recurse if leaf so no issue there
                childList = LG.descendants displayTree nodeIndex
            in  -- trace ("TBB: " <> (show (length childList, length nodeResolutionData, blockIndex, length (nodeResolutionData V.! blockIndex), fromJust resolutionIndex))) $
                if isNothing resolutionIndex
                    then
                        error
                            ( "Nothing resolution in traceBackBlock of node "
                                <> show nodeIndex
                                <> " with children "
                                <> show childList
                                <> LG.prettyIndices displayTree
                            )
                    else
                        if length childList > 2
                            then error ("Node " <> show nodeIndex <> " with > 2 children: " <> show childList)
                            else
                                if null childList
                                    then -- its a leaf
                                        (newDisplayTree, newCharTreeV)
                                    else
                                        if length childList == 1
                                            then -- recurse to left of in 1 out 1 node
                                                traceBackBlock canonicalGraph (head childList) (newDisplayTree, newCharTreeV, leftResIndex, blockIndex)
                                            else -- two children to recurse to

                                                let leftChild = head childList
                                                    rightChild = last childList
                                                    -- get left right from BV as in postorder
                                                    ((leftChild', _), (rightChild', _)) =
                                                        U.leftRightChildLabelBVNode
                                                            ((leftChild, fromJust $ LG.lab canonicalGraph leftChild), (rightChild, fromJust $ LG.lab canonicalGraph rightChild))

                                                    (leftDisplayTree, leftCharTreeV) = traceBackBlock canonicalGraph leftChild' (newDisplayTree, newCharTreeV, leftResIndex, blockIndex)
                                                    (rightDisplayTree, rightCharTreeV) = traceBackBlock canonicalGraph rightChild' (leftDisplayTree, leftCharTreeV, rightResIndex, blockIndex)
                                                in  (rightDisplayTree, rightCharTreeV)


-- )

{- | updateNodeBlockTrees takes root resolution data and sets various fields in block display
and creates character trees from block display tree
-}
updateNodeBlockTrees
    ∷ LG.Node → (ResolutionData, DecoratedGraph, V.Vector DecoratedGraph) → (DecoratedGraph, V.Vector DecoratedGraph)
updateNodeBlockTrees nodeIndex (nodeRes, displayTree, charTreeV) =
    if LG.isEmpty displayTree
        then error "Null data in updateNodeBlockTrees"
        else -- update Display tree vertex info
        -- data are a singleton vector since only one "block" (with characters) per block
        -- same, but double for characters V.singleton of V.singleton of CharacterData

            let origNodeLabel = fromJust $ LG.lab displayTree nodeIndex
                newNodeLabel =
                    origNodeLabel
                        { vertData = V.singleton $ displayData nodeRes
                        , vertexResolutionData = mempty
                        , vertexCost = resolutionCost nodeRes
                        , subGraphCost = displayCost nodeRes
                        }
                newDisplayTree = LG.updateNodeLabel displayTree nodeIndex newNodeLabel
                newCharTreeV =
                    V.zipWith
                        (updateNodeCharacterTree nodeIndex (displayData nodeRes))
                        charTreeV
                        (V.fromList [0 .. (V.length (displayData nodeRes) - 1)])
            in  -- trace ("URBT: " <> (show $ LG.prettyIndices displayTree) <> "\n" <> (show $ LG.prettyIndices newDisplayTree))
                (newDisplayTree, newCharTreeV)


-- | updateNodeCharacterTree updates an individual character tree with node info from a display tree
updateNodeCharacterTree ∷ LG.Node → V.Vector CharacterData → DecoratedGraph → Int → DecoratedGraph
updateNodeCharacterTree nodeIndex nodeBlockCharData charTree charIndex =
    if LG.isEmpty charTree
        then error "Empty tree in updateNodeCharacterTree"
        else
            let charData = nodeBlockCharData V.! charIndex
                origNodeLabel = fromJust $ LG.lab charTree nodeIndex
                newNodeLabel =
                    origNodeLabel
                        { vertData = V.singleton (V.singleton charData)
                        , vertexResolutionData = mempty
                        , vertexCost = localCost charData
                        , subGraphCost = globalCost charData
                        }
                newCharTree = LG.updateNodeLabel charTree nodeIndex newNodeLabel
            in  -- trace ("URCT: " <> (show charData))
                newCharTree


{- | updateRootBlockTrees takes root resolution data and sets various fields in block display
and creates character trees from block display tree
-}
updateRootBlockTrees ∷ LG.Node → (ResolutionData, DecoratedGraph) → (DecoratedGraph, V.Vector DecoratedGraph)
updateRootBlockTrees nodeIndex (nodeRes, displayTree) =
    if LG.isEmpty displayTree
        then error "Null data in updateNodeBlockTrees"
        else -- update Display tree vertex info
        -- data are a singleton vector since only one "block" (with characters) per block
        -- same, but double for characters V.singleton of V.singleton of CharacterData

            let origNodeLabel = fromJust $ LG.lab displayTree nodeIndex
                newNodeLabel =
                    origNodeLabel
                        { vertData = V.singleton $ displayData nodeRes
                        , vertexResolutionData = mempty
                        , vertexCost = resolutionCost nodeRes
                        , subGraphCost = displayCost nodeRes
                        }
                newDisplayTree = LG.updateNodeLabel displayTree nodeIndex newNodeLabel
                newCharTreeV =
                    fmap
                        (createNodeCharacterTree nodeIndex (displayData nodeRes) newDisplayTree)
                        (V.fromList [0 .. (V.length (displayData nodeRes) - 1)])
            in  -- trace ("URBT: " <> (show $ LG.prettyIndices displayTree) <> "\n" <> (show $ LG.prettyIndices newDisplayTree))
                (newDisplayTree, newCharTreeV)


-- | createNodeCharacterTree creates an individual character tree with node info from a display tree
createNodeCharacterTree ∷ LG.Node → V.Vector CharacterData → DecoratedGraph → Int → DecoratedGraph
createNodeCharacterTree nodeIndex nodeBlockCharData displayTree charIndex =
    if LG.isEmpty displayTree
        then error "Empty tree in updateNodeCharacterTree"
        else
            let charData = nodeBlockCharData V.! charIndex
                origNodeLabel = fromJust $ LG.lab displayTree nodeIndex
                newNodeLabel =
                    origNodeLabel
                        { vertData = V.singleton (V.singleton charData)
                        , vertexResolutionData = mempty
                        , vertexCost = localCost charData
                        , subGraphCost = globalCost charData
                        }
                newCharTree = LG.updateNodeLabel displayTree nodeIndex newNodeLabel
            in  -- trace ("URCT: " <> (show charData))
                newCharTree


{- | backPortBlockTreeNodesToCanonicalGraph takes block display trees (updated presumably) and ports the block tree node
labels to the cononical Graph--very similar to backPortCharTreeNodesToBlockTree but character vector is not singleton
back port is based on indices of block characters that may have fewer nodes than canonical in a
split graph/ sub graph situation like in swap and fuse
hence not all canonical nodes may be updated--left unchanged, so need to add back those unchanged
-}
backPortBlockTreeNodesToCanonicalGraph ∷ DecoratedGraph → V.Vector DecoratedGraph → DecoratedGraph
backPortBlockTreeNodesToCanonicalGraph inCanonicalGraph blockTreeVect =
    let canonicalNodes = LG.labNodes inCanonicalGraph
        canonicalEdges = LG.labEdges inCanonicalGraph

        -- vector (characters) of vector (nodes) of labels
        blockTreeNodeLabelsVV = fmap ((V.fromList . fmap snd) . LG.labNodes) blockTreeVect
        -- blockTreeNodeIndicesV = V.fromList $ V.head $ fmap (fmap fst) $ fmap LG.labNodes blockTreeVect
        blockTreeNodeIndicesV = V.fromList $ fmap fst $ LG.labNodes $ V.head blockTreeVect

        -- setting max index if things aren't all same length--that's the prob
        updatedCanonicalNodes =
            V.toList $
                fmap
                    (updateCanonicalNodes inCanonicalGraph blockTreeNodeLabelsVV)
                    (V.zip blockTreeNodeIndicesV (V.fromList [0 .. (length blockTreeNodeIndicesV - 1)]))
        -- maxIndex = minimum $ fmap length blockTreeNodeLabelsVV
        -- updatedCanonicalNodes = V.toList $ fmap (updateCanonicalNodes inCanonicalGraph blockTreeNodeLabelsVV) (V.zip blockTreeNodeIndicesV (V.fromList [0..(maxIndex- 1)]))

        -- add in nodes not in characters trees--can happen during swap split trees
        -- based on vertex index first field of node
        unModifiedNodes = orderedNodeMinus canonicalNodes updatedCanonicalNodes
    in  -- update BV vector of unModified nodes?

        -- trace ("BPTCG: " <> (show (fmap fst unModifiedNodes)))

        -- this is a hack if graph is imporper--not sure why htis happens yet
        if (V.null blockTreeVect) || (not $ allSameLength blockTreeNodeLabelsVV)
            then -- trace ("BPTCG: " <> (show (fmap length blockTreeNodeLabelsVV)))
                LG.empty
            else LG.mkGraph (updatedCanonicalNodes <> unModifiedNodes) canonicalEdges
    where
        allSameLength a =
            if V.null a
                then True
                else
                    if length a == 1
                        then True
                        else
                            if (length $ V.head a) /= (length $ a V.! 1)
                                then False
                                else allSameLength (V.tail a)


{- | updateCanonicalNodes takes a pair of block node index and vector of labels and
assigns data to canonical node of same index
-}
updateCanonicalNodes ∷ DecoratedGraph → V.Vector (V.Vector VertexInfo) → (LG.Node, Int) → LG.LNode VertexInfo
updateCanonicalNodes canonicalGraph blockNodeLabelVV (blockNodeIndex, vectIndex) =
    let canonicalLabel = fromJust $ LG.lab canonicalGraph blockNodeIndex
        blockNodeLabelV = fmap (V.! vectIndex) blockNodeLabelVV
        vertDataV = V.concatMap vertData blockNodeLabelV
        vertCostV = fmap vertexCost blockNodeLabelV
        subGraphCostV = fmap subGraphCost blockNodeLabelV

        -- this for Naive softwired where there is no canonical bitvecotr labelling for nodes
        -- and set to default [False]
        canonicalBV =
            if bvLabel canonicalLabel /= BV.fromBits [False]
                then bvLabel canonicalLabel
                else L.foldl1 (.|.) $ fmap bvLabel blockNodeLabelV

        -- update Info
        newVertCost = V.sum vertCostV
        newSubGraphCost = V.sum subGraphCostV
        newLabel = canonicalLabel{bvLabel = canonicalBV, vertData = vertDataV, vertexCost = newVertCost, subGraphCost = newSubGraphCost}
    in  -- trace ("UCN:" <> (show (blockNodeIndex, vectIndex, fmap length blockNodeLabelVV)))
        (blockNodeIndex, newLabel)


{- | orderedNodeMinus takes two lists of pairs where first pair is Int and
pairs are orderd by first element and returns a list of nodes in first list
not in second
assumes lists are orderd (haven't thought if non-unique elements)
should be O(n)
-}
orderedNodeMinus ∷ [(Int, a)] → [(Int, b)] → [(Int, a)]
orderedNodeMinus firstList secondList
    | null firstList = []
    | null secondList = firstList
    | otherwise =
        let firstFirst@(af, _) = head firstList
            (as, _) = head secondList
        in  if af < as
                then firstFirst : orderedNodeMinus (tail firstList) secondList
                else
                    if af == as
                        then orderedNodeMinus (tail firstList) (tail secondList)
                        else -- asf > as
                            orderedNodeMinus firstList (tail secondList)
