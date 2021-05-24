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

module GraphOptimization.Traversals  (  fullyLabelGraph
                    , postOrderTreeTraversal
                    , preOrderTreeTraversal
                    , makeLeafGraph
                   ) where

import           Types.Types
import qualified Data.Vector as V
import qualified Utilities.LocalGraph as LG
import GeneralUtilities
import Debug.Trace
import qualified GraphFormatUtilities as GFU
import qualified GraphOptimization.Medians as M
import Data.Maybe


-- | fullyLabelGraph takes an unlabelled "simple' graph, performs post and preorder passes to 
-- fully label the graph and return a PhylogeenticGraph
fullyLabelGraph :: GlobalSettings -> ProcessedData -> SimpleGraph -> PhylogeneticGraph
fullyLabelGraph inGS inData inGraph = 
    if LG.isEmpty inGraph then (LG.empty, 0.0, LG.empty, V.empty, V.empty, V.empty)
    else 
        let leafGraph = makeLeafGraph inData
            postOrderTree = postOrderTreeTraversal inGS inData leafGraph inGraph
            preOrderTree = preOrderTreeTraversal inGS inData postOrderTree
        in
        preOrderTree

-- | makeLeafGraph takes input dsata and creates a 'graph' of leaves wiyth Vertex informnation
-- biut with zero edges.  This 'graph' can be reused as a starting structure for graph construction
-- to avoid remaking of leave vertices
makeLeafGraph :: ProcessedData -> DecoratedGraph
makeLeafGraph inData@(nameVect, bvNameVect, blocDataVect) =
    if V.null nameVect then error "Empty ProcessedData in makeLeafGraph"
    else 
        let leafVertexList = V.toList $ V.map (makeLeafVertex nameVect bvNameVect blocDataVect) (V.fromList [0.. (V.length nameVect) - 1])
        in
        LG.mkGraph leafVertexList []

-- | makeLeafVertex
makeLeafVertex :: V.Vector NameText -> V.Vector NameBV -> V.Vector BlockData -> Int -> LG.LNode VertexInfo
makeLeafVertex nameVect bvNameVect inData index =
    let centralData = V.map snd3 inData 
        thisData = V.map (V.! index) centralData
        newVertex = VertexInfo  { index = index
                                , bvLabel = bvNameVect V.! index
                                , parents = V.empty
                                , children = V.empty
                                , nodeType = LeafNode
                                , vertName =  nameVect V.! index
                                , vertData = thisData
                                , vertexCost = 0.0
                                , subGraphCost = 0.0
                                }   
        in
        (index, newVertex)



-- | postOrderTreeTraversal takes a 'simple' graph and generates 'preliminary' assignments
-- vi post-order traversal, yields cost as well
-- for a binary tree only
-- depending on optimality criterion--will calculate root cost
postOrderTreeTraversal :: GlobalSettings -> ProcessedData -> DecoratedGraph -> SimpleGraph -> PhylogeneticGraph
postOrderTreeTraversal inGS inData@(nameVect, bvNameVect, blocDataVect) leafGraph inGraph = 
    if LG.isEmpty inGraph then (LG.empty, 0.0, LG.empty, V.empty, V.empty, V.empty)
    else
        -- Assumes root is Number of Leaves  
        let rootIndex = V.length $ fst3 inData
            blockCharInfo = V.map thd3 blocDataVect
            newTree = postDecorateTree inGS inData inGraph leafGraph blockCharInfo rootIndex
        in
        trace ("It Begins at " ++ show rootIndex) (
        if not $ LG.isRoot inGraph rootIndex then error ("Index "  ++ (show rootIndex) ++ " not root in graph:\n" ++ (GFU.showGraph inGraph))
        else newTree 
        )


-- | postDecorateTree begins at start index (usually root, but could be a subtree) and moves preorder till childrend ane labelled and then reurns postorder
-- labelling vertices and edges as it goes back to root
-- this for a tree so single rtoot
postDecorateTree :: GlobalSettings -> ProcessedData -> SimpleGraph -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.Node -> PhylogeneticGraph
postDecorateTree inGS inData simpleGraph curDecGraph blockCharInfo curNode = 
    -- if node in there nothing to do and return
    if LG.gelem curNode curDecGraph then 
        let nodeLabel = LG.lab curDecGraph curNode
        in
        if nodeLabel == Nothing then error ("Null label for node " ++ show curNode)
        else (simpleGraph, subGraphCost (fromJust nodeLabel), curDecGraph, V.empty, V.empty, blockCharInfo)

    -- Need to make node
    else 
        -- check if children in graph
        let nodeChildren = LG.descendants simpleGraph curNode  -- should be 1 or 2, not zero since all leaves already in graph
            leftChild = (head nodeChildren)
            rightChild = (last nodeChildren)
            leftChildTree = postDecorateTree inGS inData simpleGraph curDecGraph blockCharInfo leftChild
            rightLeftChildTree = if (length nodeChildren == 2) then postDecorateTree inGS inData simpleGraph leftChildTree blockCharInfo rightChild
                                 else leftChildTree
        in 

        if length nodeChildren > 2 then error ("Graph not dichotomous in postDecorateTree node " ++ (show curNode) ++ "\n" ++ GFU.showGraph simpleGraph)
        else if length nodeChildren == 0 then error ("Leaf not in graph in postDecorateTree node " ++ (show curNode) ++ "\n" ++ GFU.showGraph simpleGraph)

        -- make node from childern
        else    
            -- make node from children and new edges to children
            -- characters in blocks--but for tree really all same block
            let newCharData = M.median2 $ V.zip3 
                newVertex = VertexInfo {  index = curNode
                                        , bvLabel = (bvLabel leftChild) .&. (bvLabel rightChild)
                                        , parents = V.fromList $ LG.parents simpleGraph curNode
                                        , children = V.fromList nodeChildren
                                        , nodeType = LG.getNodeType simpleGraph curNode
                                        , vertName = "HTU" ++ (show curNode)
                                        , vertData = newCharData
                                        , vertexCost = 
                                        , subGraphCost = 
                                        }   
                newEdges = makeNewEdges simpleGraph nodeChildren
                newGraph =  LG.insEdges newEdges $ LG.insNode (curNode, newVertex) (thd6 curDecGraph)                         
            in
            -- return new graph
            (simpleGraph, subGraphCost, newGraph, V.empty, V.empty, blockCharInfo)

-- | preOrderTreeTraversal takes a preliminarily labelled PhylogeneticGraph
-- and returns a full labbels with 'final' assignments
-- invafiant that root is HTU !! nLeaves?
preOrderTreeTraversal :: GlobalSettings -> ProcessedData -> PhylogeneticGraph -> PhylogeneticGraph
preOrderTreeTraversal inGS inData inPGraph = 
    if LG.isEmpty (thd6 inPGraph) then (fst6 inPGraph, 0.0, LG.empty, V.empty, V.empty, V.empty)
    else 
        inPGraph
        



