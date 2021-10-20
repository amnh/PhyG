{- |
Module      :  PostOrderFunctions.hs
Description :  Module specifying post-order graph functions
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

{-
ToDo:
   Add parallel optimization overblocks and characters?
-}


module GraphOptimization.PostOrderFunctions  ( rerootPhylogeneticTree
                                             , rerootPhylogeneticTree'
                                             , rerootPhylogeneticNetwork
                                             , rerootPhylogeneticNetwork'
                                             , createVertexDataOverBlocks
                                             , divideDecoratedGraphByBlockAndCharacterTree
                                             , divideDecoratedGraphByBlockAndCharacterSoftWired
                                             , getOutDegree1VertexAndGraph
                                             , addNodeEdgeToResolutionList
                                             , extractDisplayTrees
                                             , createBlockResolutions
                                             ) where

import           Data.Bits                 
import           Data.Maybe
import qualified Data.Vector               as V
import qualified GraphOptimization.Medians as M
import           Types.Types
import qualified Utilities.LocalGraph      as LG
import qualified Utilities.Utilities       as U
import qualified Graphs.GraphOperations    as GO
import           GeneralUtilities
import qualified Data.List as L
import qualified Data.Text.Lazy  as T
import qualified Data.BitVector.LittleEndian as BV


-- import Debug.Debug
import Debug.Trace

-- | reOptimizeNodes takes a decorated graph and a list of nodes and reoptimizes (relabels)
-- them based on children in input graph
-- simple recursive since each node depends on children
-- remove check for debugging after it works
-- check for out-degree 1, doesn't matter for trees however.
reOptimizeNodes :: GlobalSettings -> GraphType -> V.Vector (V.Vector CharInfo) -> DecoratedGraph -> [LG.LNode VertexInfo] -> DecoratedGraph
reOptimizeNodes inGS localGraphType charInfoVectVect inGraph oldNodeList =
    trace ("RON:" ++ (show $ fmap fst oldNodeList)) (
  if null oldNodeList then inGraph
  else
    -- make sure that nodes are optimized in correct order so that nodes are only reoptimized using updated children
    -- this should really not have to happen--order should be determined a priori
    let curNode@(curNodeIndex, curNodeLabel) = head oldNodeList
        nodeChildren = LG.descendants inGraph curNodeIndex  -- should be 1 or 2, not zero since all leaves already in graph
        foundCurChildern = filter (`elem` nodeChildren) $ fmap fst (tail oldNodeList)
    in
    if (not $ null foundCurChildern) then
      -- trace ("Current node " ++ (show curNodeIndex) ++ " has children " ++ (show nodeChildren) ++ " in optimize list (optimization order error)" ++ (show $ fmap fst $ tail oldNodeList))
      reOptimizeNodes inGS localGraphType charInfoVectVect inGraph ((tail oldNodeList) ++ [curNode])
    
    else
        let leftChild = (head nodeChildren)
            rightChild = (last nodeChildren)
            -- leftChildLabel = fromJust $ LG.lab inGraph leftChild
            -- rightChildLabel = fromJust $ LG.lab inGraph rightChild

            -- this ensures that left/right choices are based on leaf BV for consistency and label invariance
            (leftChildLabel, rightChildLabel) = U.leftRightChildLabelBV (fromJust $ LG.lab inGraph leftChild, fromJust $ LG.lab inGraph rightChild)
            newVertexData = createVertexDataOverBlocksNonExact (vertData leftChildLabel) (vertData  rightChildLabel) charInfoVectVect []
            --newVertexData = createVertexDataOverBlocks  (vertData leftChildLabel) (vertData  rightChildLabel) charInfoVectVect []
        in
        {-
        --debug remove when not needed--checking to see if node should not be re optimized
        if (sort nodeChildren) == (sort $ V.toList $ children curnodeLabel) then
          trace ("Children for vertex unchanged " ++ (show curNodeIndex)
          reOptimizeNodes localGraphType charInfoVectVect inGraph (tail oldNodeList)
        else
        -}
        if localGraphType == Tree then 
           let newCost =  if (length nodeChildren < 2) then 0
                          else V.sum $ V.map (V.sum) $ V.map (V.map snd) newVertexData
               newVertexLabel = VertexInfo {  index = curNodeIndex
                                            -- this bit labelling incorect for outdegree = 1, need to prepend bits
                                            , bvLabel = (bvLabel leftChildLabel) .|. (bvLabel rightChildLabel)
                                            , parents = V.fromList $ LG.parents inGraph curNodeIndex
                                            , children = V.fromList nodeChildren
                                            , nodeType = nodeType curNodeLabel
                                            , vertName = vertName curNodeLabel
                                            , vertexResolutionData = mempty
                                            , vertData = if (length nodeChildren < 2) then vertData leftChildLabel
                                                         else V.map (V.map fst) newVertexData
                                            , vertexCost = newCost
                                            , subGraphCost = if (length nodeChildren < 2) then subGraphCost leftChildLabel
                                                             else (subGraphCost leftChildLabel) + (subGraphCost rightChildLabel) + newCost
                                            }
               -- this to add back edges deleted with nodes (undocumented but sensible in fgl)
               replacementEdges = (LG.inn inGraph curNodeIndex) ++ (LG.out inGraph curNodeIndex)
               newGraph = LG.insEdges replacementEdges $ LG.insNode (curNodeIndex, newVertexLabel) $ LG.delNode curNodeIndex inGraph
            in
            --trace ("New vertexCost " ++ show newCost) --  ++ " lcn " ++ (show (vertData leftChildLabel, vertData rightChildLabel, vertData curnodeLabel)))
            reOptimizeNodes inGS localGraphType charInfoVectVect newGraph (tail oldNodeList)
        
        else if  localGraphType == SoftWired then 

            -- single child of node (can certinly happen with soft-wired networks
            if length nodeChildren == 1 then 
                let (_,_, newVertexLabel, _, _) = getOutDegree1VertexAndGraph curNodeIndex leftChildLabel inGraph nodeChildren inGraph

                    -- this to add back edges deleted with nodes (undocumented but sensible in fgl)
                    replacementEdges = (LG.inn inGraph curNodeIndex) ++ (LG.out inGraph curNodeIndex)
                    newGraph = LG.insEdges replacementEdges $ LG.insNode (curNodeIndex, newVertexLabel) $ LG.delNode curNodeIndex inGraph
                in
                reOptimizeNodes inGS localGraphType charInfoVectVect newGraph (tail oldNodeList)
                
            -- two children
            else
                let -- this ensures that left/right choices are based on leaf BV for consistency and label invariance
                    -- larger bitvector is Right, smaller or equal Left 
                    ((leftChild', leftChildLabel), (rightChild', rightChildLabel)) = U.leftRightChildLabelBVNode ((leftChild, fromJust $ LG.lab inGraph leftChild), (rightChild, fromJust $ LG.lab inGraph rightChild))
                    
                    -- create resolution caches for blocks
                    leftChildNodeType  = nodeType leftChildLabel
                    rightChildNodeType = nodeType rightChildLabel
                    resolutionBlockVL = fmap (createBlockResolutions (compressResolutions inGS) curNodeIndex leftChild' rightChild' leftChildNodeType rightChildNodeType) $  V.zip3 (vertexResolutionData leftChildLabel) (vertexResolutionData rightChildLabel) charInfoVectVect
                
                    -- create canonical Decorated Graph vertex
                    -- 0 cost becasue can't know cosrt until hit root and get best valid resolutions
                    newVertexLabel = VertexInfo {  index = curNodeIndex
                                            , bvLabel = (bvLabel leftChildLabel) .|. (bvLabel rightChildLabel)
                                            , parents = V.fromList $ LG.parents inGraph curNodeIndex
                                            , children = V.fromList nodeChildren
                                            , nodeType = GO.getNodeType inGraph curNodeIndex
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

                    leftEdge =  (curNodeIndex, leftChild', edgeLable {edgeType = leftEdgeType})
                    rightEdge = (curNodeIndex, rightChild', edgeLable {edgeType = rightEdgeType})
                    newGraph =  LG.insEdges [leftEdge, rightEdge] $ LG.insNode (curNodeIndex, newVertexLabel) inGraph 
                in
                reOptimizeNodes inGS localGraphType charInfoVectVect newGraph (tail oldNodeList)    
            

        else  errorWithoutStackTrace ("Graph type unrecognized/not yet implemented: " ++ show localGraphType)
        )

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





-- | getOutDegree1VertexAndGraph makes parent node fomr single child for soft-wired resolutions            
getOutDegree1VertexAndGraph :: (Show a, Show b) 
                            => LG.Node 
                            -> VertexInfo 
                            -> LG.Gr a b 
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


-- | createVertexDataOverBlocks is a partial application of generalCreateVertexDataOverBlocks with full (all charcater) median calculation
createVertexDataOverBlocks :: VertexBlockData
                           -> VertexBlockData
                           -> V.Vector (V.Vector CharInfo)
                           -> [V.Vector (CharacterData, VertexCost)]
                           -> V.Vector (V.Vector (CharacterData, VertexCost))
createVertexDataOverBlocks = generalCreateVertexDataOverBlocks M.median2

-- | createVertexDataOverBlocksNonExact is a partial application of generalCreateVertexDataOverBlocks with partial (non-exact charcater) median calculation
createVertexDataOverBlocksNonExact :: VertexBlockData
                                   -> VertexBlockData
                                   -> V.Vector (V.Vector CharInfo)
                                   -> [V.Vector (CharacterData, VertexCost)]
                                   -> V.Vector (V.Vector (CharacterData, VertexCost))
createVertexDataOverBlocksNonExact = generalCreateVertexDataOverBlocks M.median2NonExact

-- | generalCreateVertexDataOverBlocks is a genreal version for optimizing all (Add, NonAdd, Matrix)  
-- and only non-exact (basically sequence) characters based on the median function passed
-- The function takes data in blocks and block vector of char info and
-- extracts the triple for each block and creates new block data for parent node (usually)
-- not checking if vectors are equal in length
generalCreateVertexDataOverBlocks :: (V.Vector CharacterData -> V.Vector CharacterData -> V.Vector CharInfo -> V.Vector (CharacterData, VertexCost)) 
                                  -> VertexBlockData 
                                  -> VertexBlockData 
                                  -> V.Vector (V.Vector CharInfo) 
                                  -> [V.Vector (CharacterData, VertexCost)] 
                                  -> V.Vector (V.Vector (CharacterData, VertexCost))
generalCreateVertexDataOverBlocks medianFunction leftBlockData rightBlockData blockCharInfoVect curBlockData =
    if V.null leftBlockData then
        --trace ("Blocks: " ++ (show $ length curBlockData) ++ " Chars  B0: " ++ (show $ V.map snd $ head curBlockData))
        V.fromList $ reverse curBlockData
    else
        let leftBlockLength = length $ V.head leftBlockData
            rightBlockLength =  length $ V.head rightBlockData
            -- firstBlock = V.zip3 (V.head leftBlockData) (V.head rightBlockData) (V.head blockCharInfoVect)

            -- missing data cases first or zip defaults to zero length
            firstBlockMedian = if (leftBlockLength == 0) then V.zip (V.head rightBlockData) (V.replicate rightBlockLength 0)
                               else if (rightBlockLength == 0) then V.zip (V.head leftBlockData) (V.replicate leftBlockLength 0)
                               else medianFunction (V.head leftBlockData) (V.head rightBlockData) (V.head blockCharInfoVect)
        in
        generalCreateVertexDataOverBlocks medianFunction (V.tail leftBlockData) (V.tail rightBlockData) (V.tail blockCharInfoVect) (firstBlockMedian : curBlockData)



-- | rerootPhylogeneticNetwork take a vertex index and reroots phylogenetic network
-- roots as for tree, but when a netwrok edge is chosen as root edge, its 'sister' network edge
-- is deleted (leaving two in=out=1 nodes), th ereooting takes place,
-- and the edge is re-added in both directions.  Nodes may change between
-- tree and network.  This is updated 
rerootPhylogeneticNetwork :: GlobalSettings -> Int -> PhylogeneticGraph -> PhylogeneticGraph
rerootPhylogeneticNetwork inGS rerootIndex inGraph@(inSimple, _, inDecGraph, blockDisplayForestVV, _, charInfoVectVect) = 
    if LG.isEmpty inSimple then inGraph
    else
        let originalRoots = LG.getRoots inSimple
            originalRootIndex = fst $ head originalRoots
            originalRootChildren = LG.descendants inSimple originalRootIndex
            newRootChildren = LG.descendants inSimple rerootIndex
            newRootParents = LG.parents inSimple rerootIndex
            newRootLabel = LG.lab inDecGraph rerootIndex
            parentNodeLabel = LG.lab inDecGraph $ head newRootParents
        in

        -- sanity check on graph
        if length originalRoots /= 1 then error ("Input graph has <>1 roots: " ++ (show originalRoots) ++ "\n" ++ (LG.prettify inSimple))

        -- would produce same root position
        else if rerootIndex `elem` originalRootChildren then emptyPhylogeneticGraph

        else if newRootLabel == Nothing then error ("New root has no label: " ++ show rerootIndex) 

        else if parentNodeLabel == Nothing then error ("Parent of new root has no label: " ++ show rerootIndex) 

        -- OK to reroot
        else 
            --if length newRootParents > 1 then trace ("Root on network edge") inGraph
            --else 
                -- currently just bags out on network branch
                let isNetworkNode = (nodeType $ fromJust newRootLabel) == NetworkNode 
                    parentIsNetworkNode = (nodeType $ fromJust parentNodeLabel) == NetworkNode
                in
                if isNetworkNode || parentIsNetworkNode then 
                    --trace ("RRN No roots:" ++ (show originalRootIndex) ++ " children " ++ (show originalRootChildren) ++ "-> " ++ (show rerootIndex) ++ " parents " 
                    --    ++ (show newRootParents) ++ " children " ++ (show newRootChildren)) 
                    rerootPhylogeneticTree inGS SoftWired isNetworkNode (fst $ head originalRoots) parentIsNetworkNode (head newRootParents) rerootIndex inGraph
                else rerootPhylogeneticTree inGS SoftWired False (-1) False (-1) rerootIndex inGraph

              

-- | rerootPhylogeneticNetwork' flipped version of rerootPhylogeneticNetwork
rerootPhylogeneticNetwork' :: GlobalSettings -> PhylogeneticGraph -> Int -> PhylogeneticGraph
rerootPhylogeneticNetwork' inGS inGraph rerootIndex = rerootPhylogeneticNetwork inGS rerootIndex inGraph

-- | rerootPhylogeneticTree' flipped version of rerootPhylogeneticTree
rerootPhylogeneticTree' :: GlobalSettings -> GraphType -> Bool -> Int ->  Bool -> Int -> PhylogeneticGraph -> Int -> PhylogeneticGraph
rerootPhylogeneticTree' inGS inGraphType isNetworkNode originalRootIndex parentIsNetworkNode parentIndex inGraph rerootIndex = rerootPhylogeneticTree inGS inGraphType isNetworkNode originalRootIndex parentIsNetworkNode parentIndex rerootIndex inGraph

-- | rerootGraph takes a phylogenetic graph and reroots based on a vertex index (usually leaf outgroup)
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
--   much time by consolidating--also since labels are all different--can't re-use alot of info
--   from graph to graph.
--   NNB only deals with post-order states
rerootPhylogeneticTree ::  GlobalSettings -> GraphType -> Bool -> Int ->  Bool -> Int -> Int -> PhylogeneticGraph -> PhylogeneticGraph
rerootPhylogeneticTree  inGS inGraphType isNetworkNode originalRootIndex parentIsNetworkNode parentIndex rerootIndex inPhyGraph@(inSimple, _, inDecGraph, blockDisplayForestVV, _, charInfoVectVect) =
  if LG.isEmpty inSimple || LG.isEmpty inDecGraph then error "Empty graph in rerootPhylogeneticGraph"
  --else if inCost == 0 then error ("Input graph with cost zero--likely non decorated input graph in rerootPhylogeneticGraph\n" ++ (LG.prettify $ convertDecoratedToSimpleGraph inDecGraph))
  else 
    let -- Check for netowrk node reroot which can cause cycle at root
        -- Does not triger relabelling of node-type since the branch flip keeps it a network node
        newSimpleGraph = rectifyGraph isNetworkNode originalRootIndex parentIsNetworkNode parentIndex rerootIndex inSimple

        -- decorated graph Boolean to specify that non-exact characters need to be reoptimized if affected
        -- could just update with needges? from simple graph rerooting
        newDecGraph = rectifyGraphDecorated isNetworkNode originalRootIndex parentIsNetworkNode parentIndex rerootIndex inDecGraph


        -- reoptimize nodes here
        -- nodes on spine from new root to old root that needs to be reoptimized
        (nodesToOptimize, _) = LG.pathToRoot inDecGraph (rerootIndex, fromJust $ LG.lab inDecGraph rerootIndex)

        -- reversed because ll these node edges are reversed so preorder would be in reverse orientation
        -- this only reoptimizes non-exact characters since rerooting doesn't affect 'exact" character optimization'
        newDecGraph' = reOptimizeNodes inGS inGraphType charInfoVectVect newDecGraph (reverse $ L.nub nodesToOptimize)

        -- sum of root costs on Decorated graph
        newGraphCost = sum $ fmap subGraphCost $ fmap snd $ LG.getRoots newDecGraph'

        -- rerooted diplay forests--don't care about costs--I hope (hence Bool False)
        newBlockDisplayForestVV = if V.null blockDisplayForestVV then mempty
                                  --else fmap (fmap (GO.rerootTree rerootIndex)) blockDisplayForestVV
                                  else fmap (fmap (rectifyGraph isNetworkNode originalRootIndex parentIsNetworkNode parentIndex rerootIndex)) blockDisplayForestVV

        in
        -- trace ("RPT\n" ++ (LG.prettify newSimpleGraph)) (
        if newSimpleGraph == LG.empty then emptyPhylogeneticGraph
        
        -- Same root, so no need to redo
        else if (length nodesToOptimize == 1) then emptyPhylogeneticGraph
        else
          {-
          trace ("To optimize:" ++ (show nodesToOptimize) ++ "\nOG " ++ (show inCost) ++ " :"
          ++ (LG.prettify inDecGraph) ++ "\nRRG:" ++ ((LG.prettify newDecGraph)) ++ "\nNG " ++ (show newGraphCost) ++ " :" ++ (LG.prettify newDecGraph')
          ++ "\nSG:" ++ (LG.prettify newSimpleGraph))
          -}
          -- (newSimpleGraph, newGraphCost, newDecGraph', newBlockDisplayForestVect, V.replicate (length charInfoVectVect) (V.singleton newDecGraph'), charInfoVectVect)
          if inGraphType == Tree then (newSimpleGraph, newGraphCost, newDecGraph', newBlockDisplayForestVV, divideDecoratedGraphByBlockAndCharacterTree newDecGraph', charInfoVectVect)
          else
            error "Need to update final stuff for soft-wired"
            --(newSimpleGraph, newGraphCost, newDecGraph', newBlockDisplayForestVV, divideDecoratedGraphByBlockAndCharacterSoftWired displayGraphVL', charInfoVectVect)
          

          --)

-- | rectifyGraph 'fixes' (flips) edges where a network edge has be chosen as a reroot edge
-- bnasicall at root and network edge originally 'to' network edge
rectifyGraph :: (Show a, Show b, Eq b) => Bool -> Int ->  Bool -> Int -> Int -> LG.Gr a b -> LG.Gr a b 
rectifyGraph isNetworkNode originalRootIndex parentIsNetworkNode parentIndex rerootIndex inGraph =
    if LG.isEmpty inGraph then LG.empty
    else 
        -- sanity check of phylogenetic graph 
        if not isNetworkNode && not parentIsNetworkNode then GO.rerootTree rerootIndex inGraph
        
        -- simply reroot if not network edge reroot
        else if isNetworkNode && parentIsNetworkNode then error ("Graph with parent and child nodes network vertices: " ++ (show (originalRootIndex, parentIsNetworkNode)))

        -- root point is network node
        else if isNetworkNode then
            let reRootGraph = GO.rerootTree rerootIndex inGraph
                inRootEdges = LG.inn reRootGraph originalRootIndex
                newRootEdges = fmap LG.flipLEdge inRootEdges
            in
            LG.insEdges newRootEdges $ LG.delLEdges inRootEdges reRootGraph

        -- parent of root point is network--flip and retype nodes
        else 
            let reRootGraph = GO.rerootTree rerootIndex inGraph

                -- reroot edgeFlip as above 
                inRootEdges = LG.inn reRootGraph originalRootIndex
                newRootEdges = fmap LG.flipLEdge inRootEdges

                -- new flip edge and node label change
                (inParentEdges, outParentEdges) = LG.getInOutEdges reRootGraph parentIndex
                edgeToFlip = head $ filter ((/= originalRootIndex).fst3) inParentEdges
                edgeToNetParent = head $ filter ((== originalRootIndex).fst3) inParentEdges
                nodesToRemake = [fst3 edgeToFlip, parentIndex]
                newEdge = LG.flipLEdge edgeToFlip
                newFstNodeLabel = fromJust $ LG.lab inGraph (fst3 edgeToFlip)
                newSndNodeLable = fromJust $ LG.lab inGraph parentIndex
                (fstNodeInEdges, fstNodeOutEdges) = LG.getInOutEdges reRootGraph (fst3 edgeToFlip)
                fstNodeEdgesToKeep = fstNodeInEdges ++ (fstNodeOutEdges L.\\ [edgeToFlip]) 
                newGraph = LG.insEdges (newEdge : newRootEdges) $ LG.delLEdges (edgeToFlip : inRootEdges) reRootGraph
                           
            in
            -- trace ("Edge to flip " ++ (show $ LG.toEdge edgeToFlip)) (
            if parentIndex == fst3 edgeToFlip then error ("Nodeindex confusion " ++ show parentIndex)
            else newGraph
            -- )

-- | rectifyGraph 'fixes' (flips) edges where a network edge has be chosen as a reroot edge For Decorated Graph--thee should be able to be combined
rectifyGraphDecorated :: Bool -> Int ->  Bool -> Int -> Int -> DecoratedGraph -> DecoratedGraph 
rectifyGraphDecorated isNetworkNode originalRootIndex parentIsNetworkNode parentIndex rerootIndex inGraph =
    if LG.isEmpty inGraph then LG.empty
    else 
        -- sanity check of phylogenetic graph 
        if not isNetworkNode && not parentIsNetworkNode then GO.rerootTree rerootIndex inGraph
        
        -- simply reroot if not network edge reroot
        else if isNetworkNode && parentIsNetworkNode then error ("Graph with parent and child nodes network vertices: " ++ (show (originalRootIndex, parentIsNetworkNode)))

        -- root point is network node
        else if isNetworkNode then
            let reRootGraph = GO.rerootTree rerootIndex inGraph
                inRootEdges = LG.inn reRootGraph originalRootIndex
                newRootEdges = fmap LG.flipLEdge inRootEdges
            in
            LG.insEdges newRootEdges $ LG.delLEdges inRootEdges reRootGraph

        -- parent of root point is network--flip and retype nodes
        else 
            let reRootGraph = GO.rerootTree rerootIndex inGraph

                -- reroot edgeFlip as above 
                inRootEdges = LG.inn reRootGraph originalRootIndex
                newRootEdges = fmap LG.flipLEdge inRootEdges

                -- new flip edge and node label change
                (inParentEdges, outParentEdges) = LG.getInOutEdges reRootGraph parentIndex
                edgeToFlip = head $ filter ((/= originalRootIndex).fst3) inParentEdges
                edgeToNetParent = head $ filter ((== originalRootIndex).fst3) inParentEdges
                nodesToRemake = [fst3 edgeToFlip, parentIndex]
                newEdge = LG.flipLEdge edgeToFlip
                newFstNodeLabel = fromJust $ LG.lab inGraph (fst3 edgeToFlip)
                newSndNodeLable = fromJust $ LG.lab inGraph parentIndex
                (fstNodeInEdges, fstNodeOutEdges) = LG.getInOutEdges reRootGraph (fst3 edgeToFlip)
                fstNodeEdgesToKeep = fstNodeInEdges ++ (fstNodeOutEdges L.\\ (edgeToFlip : inRootEdges))
                newFstNode = (fst3 edgeToFlip, newFstNodeLabel {nodeType = NetworkNode})
                newSndNode = (parentIndex, newSndNodeLable {nodeType = TreeNode})
                newGraph = LG.insEdges (newEdge : (fstNodeEdgesToKeep ++ newRootEdges)) $ LG.insNodes [newFstNode, newSndNode] $ LG.delNodes nodesToRemake reRootGraph
            in
            -- trace ("Edge to flip " ++ (show $ LG.toEdge edgeToFlip) ++ " Nodes to remake " ++ (show nodesToRemake) ++ " New edges " ++ (show $ fmap LG.toEdge (newEdge : fstNodeEdgesToKeep))) (
            if parentIndex == fst3 edgeToFlip then error ("Nodeindex confusion " ++ show parentIndex)
            else newGraph
            -- )
            

-- | divideDecoratedGraphByBlockAndCharacterSoftWired takes a Vector of a list of DecoratedGraph
-- continaing a list of decorated tryees that are the display trees for that block
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

-- | divideDecoratedGraphByBlockAndCharacterTree takes a DecoratedGraph with (potentially) multiple blocks
-- and (potentially) multiple character per block and creates a Vector of Vector of Decorated Graphs
-- over blocks and characyets with the same graph, but only a single block and character for each graph
-- this to be used to create the "best" cost over alternate graph traversals
-- vertexCost and subGraphCost will be taken from characterData localcost/localcostVect and globalCost
divideDecoratedGraphByBlockAndCharacterTree :: DecoratedGraph -> V.Vector (V.Vector DecoratedGraph)
divideDecoratedGraphByBlockAndCharacterTree inGraph = 
  if LG.isEmpty inGraph then V.empty
  else 
    let numBlocks = V.length $ vertData $ snd $ head $ LG.labNodes inGraph
        blockGraphList = fmap (pullBlock inGraph) [0.. (numBlocks - 1)] 
        characterGraphList = fmap makeCharacterGraph blockGraphList
    in
    -- trace ("Blocks " ++ (show numBlocks) ++ " Characters " ++ (show $ fmap length $ vertData $ snd $ head $ LG.labNodes inGraph)) 
    V.fromList characterGraphList

-- | pullBlocks take a DecoratedGraph and creates a newDecorated graph with
-- only data from the input block index
pullBlock :: DecoratedGraph -> Int -> DecoratedGraph
pullBlock inGraph blockIndex =
  if LG.isEmpty inGraph then LG.empty
  else 
    let (inNodeIndexList, inNodeLabelList) = unzip $ LG.labNodes inGraph
        blockNodeLabelList = fmap (makeBlockNodeLabels blockIndex) inNodeLabelList
    in
    LG.mkGraph (zip inNodeIndexList blockNodeLabelList) (LG.labEdges inGraph)

-- | makeBlockNodeLabels takes a block index and an orginal nodel label
-- and cretes a new list of a singleton block from the input block index
makeBlockNodeLabels :: Int -> VertexInfo -> VertexInfo
makeBlockNodeLabels blockIndex inVertexInfo =
  let newVertexData = (vertData inVertexInfo) V.! blockIndex
      newVertexCost = V.sum $ fmap localCost newVertexData
      newsubGraphCost = V.sum $ fmap globalCost newVertexData
  in
  -- trace ("MBD " ++ (show $ length newVertexData) ++ " from " ++ (show $ length (vertData inVertexInfo)))
  inVertexInfo { vertData     = V.singleton newVertexData
               , vertexCost   = newVertexCost
               , subGraphCost = newsubGraphCost
               }

-- | makeCharacterGraph takes a blockGraph and creates a vector of character graphs
-- each with a single block and single character 
-- updating costs
makeCharacterGraph :: DecoratedGraph -> V.Vector DecoratedGraph
makeCharacterGraph inBlockGraph =
  if LG.isEmpty inBlockGraph then V.empty
  else 
    let numCharacters =  V.length $ V.head $ vertData $ snd $ head $ LG.labNodes inBlockGraph
        characterGraphList = if (numCharacters > 0) then fmap (pullCharacter False inBlockGraph) [0.. (numCharacters - 1)]
                             -- missing data
                             else [pullCharacter True inBlockGraph 0]
    in
    if (V.length $ vertData $ snd $ head $ LG.labNodes inBlockGraph) /= 1 then error "Number of blocks /= 1 in makeCharacterGraph"
    else 
      -- trace ("Chars: " ++ show numCharacters)
      V.fromList characterGraphList

-- | pullCharacter takes a DecoratedGraph with a single block and
-- creates a new DecoratedGraph with a single character form the input index
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
makeCharacterLabels :: Bool -> Int -> VertexInfo -> VertexInfo
makeCharacterLabels isMissing characterIndex inVertexInfo =
  let newVertexData = (V.head $ vertData inVertexInfo) V.! characterIndex
      newVertexCost = localCost newVertexData
      newSubGraphCost = globalCost newVertexData
  in
  -- trace ("MCL " ++ (show $ V.length $ vertData inVertexInfo) ++ " " ++ (show $ fmap  V.length $ vertData inVertexInfo) )
  {-
  if V.null $ V.head $ vertData inVertexInfo then 
      inVertexInfo { vertData    = V.singleton $ V.singleton emptyCharacter
                  , vertexCost   = 0
                  , subGraphCost = 0
                  }
  else
  -} 
     inVertexInfo { vertData     = if (not isMissing) then V.singleton $ V.singleton newVertexData
                                   else V.singleton V.empty
                  , vertexCost   = newVertexCost
                  , subGraphCost = newSubGraphCost
                  }
