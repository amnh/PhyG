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

module GraphOptimization.PostOrderSoftWiredFunctionsNew  ( postDecorateSoftWiredNew
                                                         , softWiredPostOrderTraceBackNew
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
import qualified GraphFormatUtilities        as GFU
-- import qualified GraphOptimization.PreOrderFunctions as PRE
import           Debug.Trace

{-Intial Postorder softwired pass.  All functions with 'New" appended-}


-- | postDecorateSoftWiredNew begins at start index (usually root, but could be a subtree) and moves preorder till children are labelled
-- and then recurses to root postorder labelling vertices and edges as it goes
-- this for a single root
postDecorateSoftWiredNew:: GlobalSettings -> SimpleGraph -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.Node -> LG.Node -> PhylogeneticGraph
postDecorateSoftWiredNew inGS simpleGraph curDecGraph blockCharInfo rootIndex curNode =
    -- if node in current decorated graph then nothing to do and return it
    --   this because will hit node twice if network node
    if LG.gelem curNode curDecGraph then
        let nodeLabel = LG.lab curDecGraph curNode
        in
        if isNothing nodeLabel then error ("Null label for node " ++ show curNode)
        else (simpleGraph, subGraphCost (fromJust nodeLabel), curDecGraph, mempty, mempty, blockCharInfo)

    -- node is not in decorated graph--ie has not been creted/optimized    
    else
        -- get postorder assignment of children
        -- checks for single child of node
        -- result is single graph after left and right child traversals
        -- trace ("PDSW making node " ++ show curNode ++ " in\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph curDecGraph)) (
        let nodeChildren = LG.descendants simpleGraph curNode  -- should be 1 or 2, not zero since all leaves already in graph
            leftChild = head nodeChildren
            rightChild = last nodeChildren -- will be same is first for out 1 (network) node
            leftChildTree = postDecorateSoftWiredNew inGS simpleGraph curDecGraph blockCharInfo rootIndex leftChild
            rightLeftChildTree = if length nodeChildren == 2 then postDecorateSoftWiredNew inGS simpleGraph (thd6 leftChildTree) blockCharInfo rootIndex rightChild
                                 else leftChildTree
        in
        -- Checks on children
        if length nodeChildren > 2 then error ("Graph not dichotomous in postDecorateSoftWiredNew node " ++ show curNode ++ "\n" ++ LG.prettyIndices simpleGraph)
        else if null nodeChildren then error ("Leaf not in graph in postDecorateSoftWiredNew node " ++ show curNode ++ "\n" ++ LG.prettyIndices simpleGraph)

        -- after recursing to any children can optimize current node
        else
            -- make node from child block resolutions
            -- child resolution made earler in post order pass
            -- new sub tree is updated graph from children--ie not iuncluding current node
            let newSubTree = thd6 rightLeftChildTree
            in

            -- single child of node--network node
            if length nodeChildren == 1 then
                -- use left child for out degree = 1 nmodes, right should be "Nothing"
                let (newGraph, _, _, _, _) = getOutDegree1VertexAndGraphNew curNode (fromJust $ LG.lab newSubTree leftChild) simpleGraph nodeChildren newSubTree
                in
                -- display graphs and character block are not done yet since will need traceback to get preliminary states
                -- new graph ahas all soft wired decorations
                (simpleGraph, 0, newGraph, mempty, mempty, blockCharInfo)

            -- 2 children--tree node
            else
                -- trace ("Outdegree 2: " ++ (show curNode) ++ " " ++ (show $ GO.getNodeType simpleGraph curNode) ++ " Children: " ++ (show nodeChildren)) (
                -- need to create new resolutions and add to existing sets
                let -- this ensures that left/right choices are based on leaf BV for consistency and label invariance
                    -- larger bitvector is Right, smaller or equal Left
                    ((leftChild', leftChildLabel), (rightChild', rightChildLabel)) = U.leftRightChildLabelBVNode ((leftChild, fromJust $ LG.lab newSubTree leftChild), (rightChild, fromJust $ LG.lab newSubTree rightChild))

                    -- create resolution caches for blocks
                    leftChildNodeType  = GO.getNodeType simpleGraph leftChild' -- nodeType leftChildLabel
                    rightChildNodeType = GO.getNodeType simpleGraph rightChild' -- nodeType rightChildLabel
                    resolutionBlockVL = V.zipWith3 (createBlockResolutionsNew (compressResolutions inGS) curNode leftChild' rightChild' leftChildNodeType rightChildNodeType (GO.getNodeType simpleGraph curNode)) (vertexResolutionData leftChildLabel) (vertexResolutionData rightChildLabel) blockCharInfo

                    -- create canonical Decorated Graph vertex
                    -- 0 cost becasue can't know cosrt until hit root and get best valid resolutions
                    newVertexLabel = VertexInfo { index = curNode
                                                , bvLabel = bvLabel leftChildLabel .|. bvLabel rightChildLabel
                                                , parents = V.fromList $ LG.parents simpleGraph curNode
                                                , children = V.fromList nodeChildren
                                                , nodeType = GO.getNodeType simpleGraph curNode --TreeNode
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

                    (displayGraphVL, lDisplayCost) = if curNode == rootIndex then 
                                                        let (displayG, displayCost') = extractDisplayTreesNew  (Just curNode) True resolutionBlockVL
                                                        in
                                                        (displayG, displayCost')
                                                        -- (fmap (fmap LG.removeDuplicateEdges) displayG, displayCost')
                                                     else (mempty, 0.0)

                in
                -- trace ("PDWS: " ++ (show $ fmap LG.toEdge [leftEdge, rightEdge]) ++ " has edges " ++ (show $ fmap (LG.hasEdge newSubTree) $ fmap LG.toEdge [leftEdge, rightEdge]) ++ " dupes: " ++ (show $ fmap LG.getDuplicateEdges $ V.head $ fst (extractDisplayTreesNew  (Just curNode) True resolutionBlockVL)) ++ " Resolutions " ++ (show $ fmap (fmap U.hasResolutionDuplicateEdges) resolutionBlockVL))
    
                (simpleGraph, lDisplayCost, newGraph, displayGraphVL, mempty, blockCharInfo)

-- | getOutDegree1VertexAndGraphNew makes parent node from single child for soft-wired resolutions
getOutDegree1VertexAndGraphNew :: (Show a, Show b)
                            => LG.Node
                            -> VertexInfo
                            -> LG.Gr a b
                            -> [LG.Node]
                            -> DecoratedGraph
                            -> (DecoratedGraph, Bool, VertexInfo, VertexCost, V.Vector [DecoratedGraph])
getOutDegree1VertexAndGraphNew curNode childLabel simpleGraph nodeChildren subTree =

    -- trace ("In out=1: " ++ (show curNode)) (
    let childResolutionData = vertexResolutionData childLabel

        curNodeResolutionData = addNodeAndEdgeToResolutionDataNew newDisplayNode newLEdge childResolutionData

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
                                , nodeType = GO.getNodeType simpleGraph curNode -- NetworkNode
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

        -- Root node should be out degree 2 so this should not happen in general--but could during some 
        -- graph rearangemnts in fusing and swapping
        (displayGraphVL, lDisplayCost) = if nodeType newVertex == RootNode then 
                                            extractDisplayTreesNew  Nothing True (vertexResolutionData childLabel)
                                         else (mempty, 0.0)


    in
    --trace ("NV1: " ++ show newVertex)
    --trace ("GOD1VG: " ++ (show $ LG.toEdge newLEdge) ++ " has edges " ++ (show $ LG.hasEdge subTree $  LG.toEdge newLEdge) ++ "Resolutions " ++ (show $ fmap (fmap U.hasResolutionDuplicateEdges) curNodeResolutionData))
    (newGraph, nodeType newVertex == RootNode, newVertex, lDisplayCost, displayGraphVL)
    -- (newGraph, False, newVertex, 0.0, mempty)
    -- )

-- | addNodeAndEdgeToResolutionDataNew adds new node and edge to resolution data in outdegree = 1 nodes
-- straight copy would not add this node or edge to subtree in resolutions
addNodeAndEdgeToResolutionDataNew :: LG.LNode VertexInfo -> LG.LEdge EdgeInfo -> V.Vector ResolutionBlockData -> V.Vector ResolutionBlockData
addNodeAndEdgeToResolutionDataNew newNode newEdge inResBlockDataV = fmap (addNodeEdgeToResolutionBlockNew newNode newEdge True) inResBlockDataV

-- | addNodeEdgeToResolutionBlockNew adds new node and edge to resolutoin block data
addNodeEdgeToResolutionBlockNew :: LG.LNode VertexInfo -> LG.LEdge EdgeInfo -> Bool -> ResolutionBlockData -> ResolutionBlockData
addNodeEdgeToResolutionBlockNew newNode newEdge isIn1Out1Node inResBlockData  =
   V.zipWith (addNodeEdgeToResolutionListNew newNode newEdge isIn1Out1Node) inResBlockData (V.fromList [0..(V.length inResBlockData) - 1])

-- | addNodeEdgeToResolutionListNew adds new node and edge to single subGraph in ResolutionData
-- adds resolution pairs to be equal to the child straight one-for-one correpondance
-- although only a single child-both indieces are set to resolution index since singels can be added to
--paired resolitions if one or other is a network node
addNodeEdgeToResolutionListNew :: LG.LNode VertexInfo -> LG.LEdge EdgeInfo -> Bool -> ResolutionData -> Int -> ResolutionData
addNodeEdgeToResolutionListNew newNode newEdge isIn1Out1Node inResData resolutionIndex =
    let (inNodeList, inEdgeList) = displaySubGraph inResData
            
        -- childResolutionIndexPairList = childResolutions inResData
        newNodeList = newNode : inNodeList
            
        -- this check for redundant edges in resolution cash from combinations
        newEdgeList = if (newEdge `notElem` inEdgeList) then newEdge : inEdgeList
                      else trace ("Should not happen: Extra edge in addNodeEdgeToResolutionListNew") inEdgeList
        newFirstData = inResData { displaySubGraph  = (newNodeList, newEdgeList)
                                   -- both set because can be a display node left right added to 2 child resolutoins
                                   , childResolutions = [(Just resolutionIndex, Just resolutionIndex)]
                                   , childResolutionIndices = (Just resolutionIndex, Just resolutionIndex)
                                   }
    in
    newFirstData 

-- | createBlockResolutionsNew takes left and right child resolution data for a block (same display tree)
-- and generates node resolution data
createBlockResolutionsNew :: Bool
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
createBlockResolutionsNew
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
        -- need to keep these indices correct (hence reverse in checkLeafOverlapNew ) for traceback and compress
        childResolutionIndices = cartProd [0.. (length leftChild - 1)] [0.. (length rightChild - 1)]
        validPairs = concatMap checkLeafOverlapNew  (zip childResolutionPairs childResolutionIndices)

        -- either parallel seems about the same
        -- newResolutionList = fmap (createNewResolutionNew curNode leftIndex rightIndex leftChildNodeType rightChildNodeType charInfoV) validPairs `using` PU.myParListChunkRDS
        newResolutionList = PU.seqParMap rdeepseq  (createNewResolutionNew curNode leftIndex rightIndex leftChildNodeType rightChildNodeType charInfoV) validPairs 

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
                            newRightChildBlockResolutionData = addNodeEdgeToResolutionBlockNew newNode newEdge False rightChild 
                        in
                        -- trace ("ANEL:" ++ (show $ (curNode, rightIndex)))
                        newRightChildBlockResolutionData
                   else -- trace ("ANEL-Nothing")
                        mempty

        addRight = if rightChildNodeType == NetworkNode then
                        let newEdge = (curNode, leftIndex, edgeLable)
                            newLeftChildBlockResolutionData = addNodeEdgeToResolutionBlockNew newNode newEdge False leftChild 
                        in
                        -- trace ("ANER:" ++ (show $ (curNode, leftIndex)))
                        newLeftChildBlockResolutionData
                   else -- trace ("ANER-Nothing")
                        mempty


    in
   --  trace ("CBR:" ++ (show $ leftIndex == rightIndex)) (
    -- trace ("=> " ++ (show $fmap BV.toBits $ fmap displayBVLabel totalResolutions) )
    -- compress to unique resolutions--can loose display trees, but speed up post-order a great deal
    -- trace ("CBR Num res left: " ++ (show $ V.length leftChild) ++ " Num res right: " ++ (show $ V.length rightChild) ++ " =>NRL " ++ (show $ length newResolutionList) ++ " addleft " ++ (show $ length addLeft)  ++ " addright " ++ (show $ length addRight)) (
    -- trace ("CBR: " ++ (show (length childResolutionPairs, length childResolutionIndices, length validPairs)))
    V.fromList newResolutionList V.++ (addLeft V.++ addRight)
    -- )
    -- )

-- | createNewResolutionNew takes a pair of resolutions and creates the median resolution
-- need to watch let/right (based on BV) for preorder stuff
createNewResolutionNew :: LG.Node 
                    -> Int 
                    -> Int 
                    -> NodeType 
                    -> NodeType 
                    -> V.Vector CharInfo 
                    -> ((ResolutionData, ResolutionData),(Int, Int)) 
                    -> ResolutionData
createNewResolutionNew curNode leftIndex rightIndex leftChildNodeType rightChildNodeType charInfoV ((leftRes, rightRes), (leftResIndex, rightResIndex)) =
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

        -- Data fields empty for display tree data--not needed and multiple copies of everything
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

        -- this check for redundant edges in resoluton cash from combinations
        -- resolutionEdgeList = leftEdge : (rightEdge: (snd leftChildTree ++ snd rightChildTree))
        -- is this required?
        existingEdges = snd leftChildTree ++ snd rightChildTree
        resolutionEdgeList = if (leftEdge `notElem` existingEdges) && (rightEdge `notElem` existingEdges) then leftEdge : (rightEdge : existingEdges)
                             else if (leftEdge `elem` existingEdges) && (rightEdge `elem` existingEdges) then existingEdges
                             else if (leftEdge `notElem` existingEdges) then leftEdge : existingEdges
                             else rightEdge : existingEdges

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
    -- trace ("CNR: " ++ (show (leftResIndex, rightResIndex)))
    ResolutionData { displaySubGraph = (resolutionNodeList, resolutionEdgeList)
                   , displayBVLabel = resBV
                   , displayData = resolutionMedianV
                   , childResolutions = [(Just leftResIndex, Just rightResIndex)]
                   , childResolutionIndices = (Just leftResIndex, Just rightResIndex)
                   , resolutionCost = thisResolutionCost
                   , displayCost = displaySubTreeCost
                   }

-- | extractDisplayTreesNew  takes resolutions and pulls out best cost (head for now) need to change type for multiple best
-- option for filter based on pop-count for root cost and complete display tree check
extractDisplayTreesNew  :: Maybe Int -> Bool -> V.Vector ResolutionBlockData -> (V.Vector [DecoratedGraph], VertexCost)
extractDisplayTreesNew  startVertex checkPopCount inRBDV =
    if V.null inRBDV then (V.empty, 0.0)
    else
        let (bestBlockDisplayResolutionList, costVect, _) = V.unzip3 $ fmap (getBestResolutionListNew  startVertex checkPopCount) inRBDV
        in
        -- trace ("EDT: " ++ (show (V.length bestBlockDisplayResolutionList, fmap length bestBlockDisplayResolutionList, V.sum costVect)))
        (bestBlockDisplayResolutionList, V.sum costVect)

-- | checkLeafOverlapNew  takes a left right resolution pair list and checks if
-- there is leaf overlap via comparing displayBVLabel if & = 0 then no
-- overlap, and adds to resulting list--reverses order--sholdn't matter
checkLeafOverlapNew  :: ((ResolutionData, ResolutionData), (Int, Int))-> [((ResolutionData, ResolutionData), (Int, Int))]
checkLeafOverlapNew  inPair@((leftRes, rightRes), (_, _))  =
    let leftBV = displayBVLabel leftRes
        rightBV = displayBVLabel rightRes
    in
    if BV.isZeroVector (leftBV .&. rightBV) then [inPair]
    else []

-- | getBestResolutionListNew  takes ResolutionBlockData and retuns a list of the best valid (ie all leaves in subtree) display trees
-- for that block with cost and resolujtion data in triple-- if checkPopCount is True--otherwise all display trees of any cost and contitution
-- startVertex for a component-- to allow for not every leaf being in componnet but still getting softwired cost
getBestResolutionListNew  :: Maybe Int -> Bool -> ResolutionBlockData -> ([DecoratedGraph], VertexCost, ResolutionBlockData)
getBestResolutionListNew  startVertex checkPopCount inRDList =
    -- trace ("GBRL: " ++ (show $ V.length inRDList)) (
    if null inRDList then error "Null resolution list"
    else
        let displayTreeList = fmap displaySubGraph inRDList
            displayCostList = fmap displayCost inRDList
            displayPopList = fmap (complement . displayBVLabel) inRDList
        in
        if not checkPopCount then
            let minCost = minimum displayCostList
                displayCostTripleList = V.zip3 displayTreeList displayCostList inRDList
                (bestDisplayList, _, bestResList) = V.unzip3 $ V.filter ((== minCost) . snd3) displayCostTripleList
            in
            (fmap LG.mkGraphPair (V.toList bestDisplayList), minCost, bestResList)
        else
            let minPopCount = minimum $ fmap popCount displayPopList  --max since complemented above
                displayBVList = V.zip4 displayTreeList displayCostList displayPopList inRDList

                -- must have all leaves if startvzertex == Nothing, component maximum otherwise
                -- this for getting cost of component of a softwired network
                validDisplayList = if startVertex == Nothing then V.filter (BV.isZeroVector . thd4) displayBVList
                                   else V.filter ((== minPopCount) . (popCount . thd4)) displayBVList
                validMinCost = V.minimum $ fmap snd4 validDisplayList
                (bestDisplayList, _, _, bestResList) = V.unzip4 $ V.filter ((== validMinCost) . snd4) validDisplayList
            in
            --trace ("Valid display list number:" ++ (show $ length validDisplayList)) (
            if (startVertex == Nothing) && (V.null validDisplayList) then error ("Null root validDisplayList in getBestResolutionListNew " ++ (show (startVertex,inRDList)) ++ " This can be caused if the graphType not set correctly.")
            else
                let lDisplayTreeList = fmap LG.mkGraphPair (V.toList bestDisplayList)

                    -- update root cost of display trees for use later (e.g. net penalties, outputting display forrests)
                    lDisplayTreeList' = fmap (updateRootCostNew  validMinCost) lDisplayTreeList
                in
                -- trace ("GBRL: " ++ (show (length lDisplayTreeList', validMinCost)))
                (lDisplayTreeList', validMinCost, bestResList)
      -- )

-- | updateRootCostNew  updates the subGraphCost of the root node(s) with input value
-- new node is created, so original is deleted, new added, and original edges added back
-- since deleted when node is
-- assumes its a tree wiht a single root
updateRootCostNew  :: VertexCost -> DecoratedGraph -> DecoratedGraph
updateRootCostNew  newRootCost inGraph =
    let (rootIndex, rootLabel) = head $ LG.getRoots inGraph
        rootEdges = LG.out inGraph rootIndex
        newRootLabel = rootLabel {subGraphCost = newRootCost}
    in
    -- trace ("DCC: " ++ (show newRootCost))
    LG.insEdges rootEdges $ LG.insNode (rootIndex, newRootLabel) $ LG.delNode rootIndex inGraph


{-
   Traceback code for after intial Postorder softwired pass.  All functions with 'New" appended
-}

-- | softWiredPostOrderTraceBackNew  takes resolution data and assigns correct resolution median 
-- from vertexResolutionData to preliminary data assignments in vertData.  
-- Proceeds via typical pre-order pass over display tree for each block
-- using the indices of left and right (if present) of resolutions
-- first gets root assignment form resolution data and then each block is traversed given its block display tree
softWiredPostOrderTraceBackNew  :: Int -> PhylogeneticGraph -> PhylogeneticGraph
softWiredPostOrderTraceBackNew  rootIndex inGraph@(a, b, canonicalGraph, _, _, f)  =
    if LG.isEmpty canonicalGraph then emptyPhylogeneticGraph
    else
      let -- extract display trees and bloxck char trees from PhylogeneticGraph
          -- block character trees do not exist yet
          displayTreeV = fmap head $ fth6 inGraph
          
          -- get root node resolution data from canonical Graph created in postorder
          rootLabel = fromJust $ LG.lab canonicalGraph rootIndex
          rootResData = vertexResolutionData rootLabel
          
          -- extract (first) best resolution for each block--there can be more than one for each, but only use the first for
          -- traceback, preliminary and final assignment etc--part of the heuristic
          (_, _, rootDisplayBlockCharResolutionV) = V.unzip3 $ PU.seqParMap rdeepseq (getBestResolutionListNew Nothing True) rootResData
          firstOfEachRootRes = fmap V.head rootDisplayBlockCharResolutionV

          -- get preliminary character data for blocks
          -- rootPreliminaryDataVV = fmap displayData firstOfEachRootRes
          (leftIndexList, rightIndexList) = V.unzip $ fmap childResolutionIndices firstOfEachRootRes

          -- update root vertex info for display and character trees for each block
          -- this includes preliminary data and other fields
          (rootUpdatedDisplayTreeV, rootUpdatedCharTreeVV) = V.unzip $ PU.seqParMap rdeepseq (updateRootBlockTrees rootIndex) (V.zip  firstOfEachRootRes displayTreeV)

          -- traceback for each block based on its display tree, updating trees as it goes, left descendent then right 
          -- at this stage all character trees will have same root descendents sionce all rooted from outgropu postorder traversal
          -- later (after reroot pass) this will not be the case since each charcater may have a unique traversal root/edge
          [leftChild, rightChild] = take 2 $ LG.descendants canonicalGraph rootIndex
          (traceBackDisplayTreeVLeft, traceBackCharTreeVVLeft) = V.unzip $ PU.seqParMap rdeepseq (traceBackBlock canonicalGraph leftChild) (V.zip4 rootUpdatedDisplayTreeV rootUpdatedCharTreeVV leftIndexList (V.fromList [0..(V.length rootUpdatedDisplayTreeV - 1)]))
          (traceBackDisplayTreeV, traceBackCharTreeVV) = V.unzip $ PU.seqParMap rdeepseq (traceBackBlock canonicalGraph rightChild) (V.zip4 traceBackDisplayTreeVLeft traceBackCharTreeVVLeft rightIndexList (V.fromList [0..(V.length rootUpdatedDisplayTreeV - 1)]))

      in
      if length (LG.descendants canonicalGraph rootIndex) /= 2 then error ("Root has improper number of children: " ++ (show $ LG.descendants canonicalGraph rootIndex))
      else 
         -- trace ("SWTN: " ++ (show (V.length rootDisplayBlockCharResolutionV,  V.length rootPreliminaryDataVV, fmap V.length rootPreliminaryDataVV )) 
         -- ++ "\n" ++ (show rootBlockChildIndicesV))
         --trace ("SWTN: " ++ (show (leftChild, rightChild)) ++ "\nLeft: " ++ (show leftIndexList) ++ "\nRight: " ++ (show rightIndexList)) 
         -- trace ("SWTN: " ++ (show $ V.length traceBackDisplayTreeVLeft) ++ " " ++ (show $ V.length traceBackCharTreeVVLeft) ++ " " ++ (show $ V.length traceBackDisplayTreeV)
         --    ++ " " ++ (show $ fmap length traceBackCharTreeVV))
         (a, b, canonicalGraph, fmap (:[]) traceBackDisplayTreeV, traceBackCharTreeVV, f)

-- | traceBackBlock performs softwired traceback on block data returns updated display and character trees
-- the block index specifies which resolution listy from the caninical tree at node nodeIndex
-- the resoliution index is teh resolution element that was used to create the parent's state
-- need to account for out degree 2 and 1 in recursive calls
traceBackBlock :: DecoratedGraph -> LG.Node -> (DecoratedGraph, V.Vector DecoratedGraph, Maybe Int, Int) -> (DecoratedGraph, V.Vector DecoratedGraph)
traceBackBlock canonicalGraph nodeIndex (displayTree, charTreeV, resolutionIndex, blockIndex) =
   -- trace ("TBB: " ++ "Node " ++ (show nodeIndex) ++ " Block " ++ (show blockIndex) ++ " Resolution " ++ (show (resolutionIndex, LG.descendants displayTree nodeIndex))) ( 
   if (LG.isEmpty displayTree) || V.null charTreeV then error "Null data in traceBackBlock"
   else 
      let -- get block resolution data from canonical graph
          nodeCanonicalLabel = fromJust $ LG.lab canonicalGraph nodeIndex
          nodeResolutionData = vertexResolutionData nodeCanonicalLabel
          blockResolutionData = (nodeResolutionData V.! blockIndex) V.! (fromJust resolutionIndex)

          -- update display tree and character tree nodes
          (newDisplayTree, newCharTreeV) = updateNodeBlockTrees nodeIndex (blockResolutionData, displayTree, charTreeV)

          -- get left and right display resolution indices
          (leftResIndex, rightResIndex) = childResolutionIndices blockResolutionData

          -- get left and right node children, won't recurse if leaf so no issue there
          childList = LG.descendants displayTree nodeIndex

      in
      if isNothing resolutionIndex then error ("Nothing resolution in traceBackBlock of node " ++ (show nodeIndex) ++ " with children " ++ (show childList)
         ++ (LG.prettyIndices displayTree))
      else if length childList > 2 then error ("Node " ++ (show nodeIndex) ++ " with > 2 children: " ++ (show childList))
      else 
         if null childList then 
            -- its a leaf 
            (newDisplayTree, newCharTreeV)
         else if length childList == 1 then
            -- recurse to left of in 1 out 1 node
            traceBackBlock canonicalGraph (head childList) (newDisplayTree, newCharTreeV, leftResIndex, blockIndex)
         else 
            -- two children to recurse to
            let [leftChild, rightChild] = take 2 childList
                -- get left right from BV as in postorder
                ((leftChild', _), (rightChild', _)) = U.leftRightChildLabelBVNode ((leftChild, fromJust $ LG.lab canonicalGraph leftChild), (rightChild, fromJust $ LG.lab canonicalGraph rightChild))

                (leftDisplayTree, leftCharTreeV) = traceBackBlock canonicalGraph leftChild' (newDisplayTree, newCharTreeV, leftResIndex, blockIndex)
                (rightDisplayTree, rightCharTreeV) = traceBackBlock canonicalGraph rightChild' (leftDisplayTree, leftCharTreeV, rightResIndex, blockIndex)
            in
            (rightDisplayTree, rightCharTreeV)

        -- )

-- | updateNodeBlockTrees takes root resolution data and sets various fields in block display 
-- and creates character trees from block display tree
updateNodeBlockTrees :: LG.Node -> (ResolutionData, DecoratedGraph, V.Vector DecoratedGraph) -> (DecoratedGraph, V.Vector DecoratedGraph)
updateNodeBlockTrees nodeIndex (nodeRes, displayTree, charTreeV) =
   if LG.isEmpty displayTree then error "Null data in updateNodeBlockTrees"
   else 
      -- update Display tree vertex info
      -- data are a singleton vector since only one "block" (with characters) per block
      -- same, but double for characters V.singleton of V.singleton of CharacterData
      let origNodeLabel = fromJust $ LG.lab displayTree nodeIndex
          newNodeLabel = origNodeLabel { vertData     = V.singleton $ displayData nodeRes
                                       , vertexResolutionData = mempty
                                       , vertexCost   = resolutionCost nodeRes
                                       , subGraphCost = displayCost nodeRes
                                       }
          newDisplayTree = LG.updateNodeLabel displayTree nodeIndex newNodeLabel
          newCharTreeV = V.zipWith (updateNodeCharacterTree nodeIndex (displayData nodeRes)) charTreeV (V.fromList [0..((V.length (displayData nodeRes)) - 1)])
      in
      -- trace ("URBT: " ++ (show $ LG.prettyIndices displayTree) ++ "\n" ++ (show $ LG.prettyIndices newDisplayTree))
      (newDisplayTree, newCharTreeV)

-- | updateNodeCharacterTree updates an individual character tree with node info from a display tree
updateNodeCharacterTree :: LG.Node -> V.Vector CharacterData -> DecoratedGraph -> Int -> DecoratedGraph
updateNodeCharacterTree nodeIndex nodeBlockCharData charTree charIndex =
   if LG.isEmpty charTree then error "Empty tree in updateNodeCharacterTree"
   else
      let charData = nodeBlockCharData V.! charIndex
          origNodeLabel = fromJust $ LG.lab charTree nodeIndex
          newNodeLabel = origNodeLabel { vertData     = V.singleton (V.singleton charData)
                                       , vertexResolutionData = mempty
                                       , vertexCost   = localCost charData
                                       , subGraphCost = globalCost charData
                                       }
          newCharTree = LG.updateNodeLabel charTree nodeIndex newNodeLabel
      in
      -- trace ("URCT: " ++ (show charData))
      newCharTree

-- | updateRootBlockTrees takes root resolution data and sets various fields in block display 
-- and creates character trees from block display tree
updateRootBlockTrees :: LG.Node -> (ResolutionData, DecoratedGraph) -> (DecoratedGraph, V.Vector DecoratedGraph)
updateRootBlockTrees nodeIndex (nodeRes, displayTree) =
   if LG.isEmpty displayTree then error "Null data in updateNodeBlockTrees"
   else 
      -- update Display tree vertex info
      -- data are a singleton vector since only one "block" (with characters) per block
      -- same, but double for characters V.singleton of V.singleton of CharacterData
      let origNodeLabel = fromJust $ LG.lab displayTree nodeIndex
          newNodeLabel = origNodeLabel { vertData     = V.singleton $ displayData nodeRes
                                       , vertexResolutionData = mempty
                                       , vertexCost   = resolutionCost nodeRes
                                       , subGraphCost = displayCost nodeRes
                                       }
          newDisplayTree = LG.updateNodeLabel displayTree nodeIndex newNodeLabel
          newCharTreeV = fmap (createNodeCharacterTree nodeIndex (displayData nodeRes) newDisplayTree) (V.fromList [0..((V.length (displayData nodeRes)) - 1)])
      in
      -- trace ("URBT: " ++ (show $ LG.prettyIndices displayTree) ++ "\n" ++ (show $ LG.prettyIndices newDisplayTree))
      (newDisplayTree, newCharTreeV)

-- | createNodeCharacterTree creates an individual character tree with node info from a display tree
createNodeCharacterTree :: LG.Node -> V.Vector CharacterData -> DecoratedGraph -> Int -> DecoratedGraph
createNodeCharacterTree nodeIndex nodeBlockCharData displayTree charIndex =
   if LG.isEmpty displayTree then error "Empty tree in updateNodeCharacterTree"
   else
      let charData = nodeBlockCharData V.! charIndex
          origNodeLabel = fromJust $ LG.lab displayTree nodeIndex
          newNodeLabel = origNodeLabel { vertData     = V.singleton (V.singleton charData)
                                       , vertexResolutionData = mempty
                                       , vertexCost   = localCost charData
                                       , subGraphCost = globalCost charData
                                       }
          newCharTree = LG.updateNodeLabel displayTree nodeIndex newNodeLabel
      in
      -- trace ("URCT: " ++ (show charData))
      newCharTree


-- | softWiredPostOrderTraceBackNew'  takes resolution data and assigns correct resolution median to preliminary
-- data ssignments.  Proceeds via typical pre-order pass over display tree for each block
softWiredPostOrderTraceBackNew'  :: Int -> DecoratedGraph -> DecoratedGraph
softWiredPostOrderTraceBackNew'  rootIndex inGraph  =
    if LG.isEmpty inGraph then LG.empty
    else
            -- get edges to remake graph after nodes are updated with preliminary states
        let rootLabel = fromJust $ LG.lab inGraph rootIndex
            inEdgeList = LG.labEdges inGraph

            -- root first--choose best resolutions--returns a list and takes head of that list of equal cost/traceback preliminary
            -- assignments.  Later could look at multiple.
            (rootVertData, rootSubGraphCost, rootResolutionCost, leftRightIndexVect) = getResolutionDataAndIndicesNew  rootLabel (V.singleton (Just (-1)))
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
            softWiredUpdatedNodes = softWiredPrelimTracebackNew  inGraph rootChildrenTuples [newRootNode]

        in
        LG.mkGraph softWiredUpdatedNodes inEdgeList

-- | softWiredPrelimTracebackNew  takes a list of nodes to update (with left right index info) based
-- on resolution data, recurses with left right indices pre-order to leaves, keeping list of updated nodes
softWiredPrelimTracebackNew  :: DecoratedGraph
                        -> [(LG.LNode VertexInfo, V.Vector (Maybe Int, Maybe Int), Bool)]
                        -> [LG.LNode VertexInfo]
                        -> [LG.LNode VertexInfo]
softWiredPrelimTracebackNew  inGraph nodesToUpdate updatedNodes =
    if null nodesToUpdate then updatedNodes
    else
        let (firstNode, firstLeftRight, isLeft) = head nodesToUpdate

            -- ensure consistent left/right from post-order
            resolutionIndexVect = if isLeft then fmap fst firstLeftRight
                                  else fmap snd firstLeftRight

            -- get resolution info
            (newBlockData, newSubGraphCost, newVertexCost, childLeftRightIndexVect) = getResolutionDataAndIndicesNew  (snd firstNode) resolutionIndexVect

            -- make new node
            newNodeLabel = (snd firstNode) { vertData = newBlockData
                                           , vertexCost = newVertexCost
                                           , subGraphCost = newSubGraphCost
                                           }

            newFirstNode = (fst firstNode, newNodeLabel)
        in

        -- not really necessary, but avoids the graph query operations
        if nodeType (snd firstNode) == LeafNode then
            softWiredPrelimTracebackNew  inGraph (tail nodesToUpdate) (newFirstNode : updatedNodes)

        -- checks if network node (and its children) has been visited already
        else if (nodeType (snd firstNode) == NetworkNode) && isJust (L.find ((== fst firstNode). fst) updatedNodes) then
            softWiredPrelimTracebackNew  inGraph (tail nodesToUpdate) updatedNodes

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
        if V.head resolutionIndexVect == Nothing then error ("'Nothing' index in softWiredPrelimTracebackNew  " ++ (show resolutionIndexVect)
            ++ " Node " ++ (show $ fst firstNode) ++ " isLeft?: " ++ (show isLeft)  ++ " " ++ (show firstLeftRight))
        else
        -}
            softWiredPrelimTracebackNew  inGraph (childrenTuples ++ tail nodesToUpdate) (newFirstNode : updatedNodes)


-- | getResolutionDataAndIndicesNew  takes a vertex label (VertexInfo) and returns the resolution data corresponding to
-- the index taken from its child resolution (data, subgraph cost, local resolution cost, left/right pairs).
-- Index = (-1) denotes that it is a root label and in that case
-- the best (lowest cost) resolutions are returned
getResolutionDataAndIndicesNew  :: VertexInfo -> V.Vector (Maybe Int) -> (VertexBlockData, VertexCost, VertexCost, V.Vector (Maybe Int, Maybe Int))
getResolutionDataAndIndicesNew  nodeLabel parentResolutionIndexVect =

    -- should not happen
    --mtrace ("NL " ++ (show $ index nodeLabel) ++ " PRIL: " ++ " length " ++ (show $ V.length parentResolutionIndexVect) ++ " " ++ show parentResolutionIndexVect) (
    if nodeType nodeLabel == LeafNode then
        let leafVertData = fmap (displayData . V.head) (vertexResolutionData nodeLabel)
        in
        -- trace ("GRDI: Leaf " ++ (show $ fmap V.length (vertexResolutionData nodeLabel)))
        (leafVertData, 0, 0, V.singleton (Just 0, Just 0))

    -- root node--take lowest cost
    else if V.head parentResolutionIndexVect == Just (-1) then
        let rootBlockResolutionPair = getBestBlockResolutionNew  <$> vertexResolutionData nodeLabel
            (charDataVV, subGraphCostV, resCostV, leftRightIndexVect) = V.unzip4 rootBlockResolutionPair
        in
        -- trace ("GRDI: Root " ++ (show $ fmap show parentResolutionIndexVect) ++ " " ++ (show $ fmap V.length (vertexResolutionData nodeLabel)))
        (charDataVV, V.sum subGraphCostV, V.sum resCostV, leftRightIndexVect)

    -- non-root node--return the index resolution information
    else
        -- trace ("GRD Length:" ++ (show $ V.length $ vertexResolutionData nodeLabel) ++ " " ++ (show $ fmap fromJust parentResolutionIndexVect) ++ "\n" ++ (show $ vertexResolutionData nodeLabel)) (
        let parentIndexVect = fmap fromJust parentResolutionIndexVect

            -- get resolution data from node label
            resolutionData = vertexResolutionData nodeLabel

            -- get the correct (via index) resolution data for each block
            -- complex for network node since keeps left right sort of array, but only first element maters--this hack keeps things ok for
            -- tree-like traceback assignment
            -- not sure if first or last would be better--hopefully not matter, or arbitrary equal solutions
            resolutionsByBlockV = if nodeType nodeLabel == NetworkNode then
                                        -- trace ("-" ++ (show (V.length resolutionData, V.length parentIndexVect, V.head parentIndexVect)) ++ " " ++ (show $ parentIndexVect))
                                        -- V.zipWith (V.!) resolutionData (V.replicate (V.length parentIndexVect) (V.head parentIndexVect))
                                        V.zipWith (V.!) resolutionData (V.replicate (V.length parentIndexVect) 0)
                                  else 
                                    -- trace ("+" ++ (show (V.length resolutionData, V.length parentIndexVect)) ++ " " ++ (show $ parentIndexVect))
                                    V.zipWith (V.!) resolutionData parentIndexVect

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
           error ("Null resolution data in getResolutionDataAndIndicesNew  at node with label:" ++ (show nodeLabel))
        else
        -}
        -- trace ("GRDI: " ++ (show $ nodeType nodeLabel) ++ " " ++ (show $ fmap V.length resolutionData) ++ " " ++ (show parentIndexVect))
        (charDataVV, lSubGraphCost, localResolutionCost, leftRightIndexVect)
        -- )

-- | getBestBlockResolutionNew  takes vertexResolutionData and returns the best (lowest cost) resolution and associated data
-- for a single block of characters.  A single left right index pair is returns for child resolution sources. Could be multiple
-- from resolution compression (equal leaf set, equal median) and avaialble for potenitally implementd in future, multiple
-- preliminary assignments.  Only valid for root node with all leaves in graph
getBestBlockResolutionNew  :: ResolutionBlockData -> (V.Vector CharacterData, VertexCost, VertexCost, (Maybe Int, Maybe Int))
getBestBlockResolutionNew  inResBlockData =
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
        if null validVect' then error "Null valid quad in getBestBlockResolutionNew --perhaps not root node or forest component"
        else (V.head displayMedianV, V.head displayCostV, V.head resCostV, V.head childIndexPairV)
