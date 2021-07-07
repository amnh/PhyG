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

module GraphOptimization.Traversals ( fullyLabelGraph
                                    , postOrderTreeTraversal
                                    , preOrderTreeTraversal
                                    , makeLeafGraph
                                    , multiTraverseFullyLabelGraph
                                    ) where

import           Types.Types
import Data.List
import qualified Data.Vector as V
import qualified Utilities.LocalGraph as LG
import GeneralUtilities
import Debug.Trace
import qualified GraphFormatUtilities as GFU
import qualified GraphOptimization.Medians as M
import qualified Graphs.GraphOperations  as GO

import qualified Data.BitVector.LittleEndian as BV
import qualified Data.Text.Lazy  as T
import Data.Bits ((.&.), (.|.))
import Data.Maybe

-- | multiTraverseFullyLabelGraph naively reroot (no progessive reroot) graph on all vertices and chooses lowest cost
-- does not yet minimize over characters to get multi-min
multiTraverseFullyLabelGraph :: GlobalSettings -> ProcessedData -> SimpleGraph -> PhylogeneticGraph
multiTraverseFullyLabelGraph inGS inData inGraph =
    if LG.isEmpty inGraph then (LG.empty, 0.0, LG.empty, V.empty, V.empty, V.empty)
    else 
        let leafGraph = makeLeafGraph inData

            {-
            -- brute force
            rootList = [0.. (( 2 * (V.length $ fst3 inData)) - 1)] -- need a smarter list going to adjecent edges
            rerootSimpleList = fmap (GO.rerootGraph' inGraph) rootList
            rerootedPhyloGraphList = fmap (fullyLabelGraph inGS inData leafGraph) rerootSimpleList --  (take 1 rerootSimpleList)
            minCost = minimum $ fmap snd6 rerootedPhyloGraphList
            minCostGraphList = filter ((== minCost).snd6) rerootedPhyloGraphList
            
            -- minimal reoptimize--but order naive
            rerootPhyloGraphListDirect = fmap (GO.rerootPhylogeneticGraph' (head rerootedPhyloGraphList)) rootList
            minCostDirect = minimum $ fmap snd6 rerootPhyloGraphListDirect
            minCostGraphListDirect = filter ((== minCost).snd6) rerootPhyloGraphListDirect
            -}

            -- minimal reoptimize with smart order to minimize node reoptimization

            -- initial travesal based on global outgroup and the "next" traversal points as children of existing traversal
            -- here initial root.
            outgroupRootedPhyloGraph = fullyLabelGraph inGS inData leafGraph $ GO.rerootGraph' inGraph (outgroupIndex inGS)
            childrenOfRoot = concat $ fmap (LG.descendants (thd6 outgroupRootedPhyloGraph)) (fmap fst $ LG.getRoots $ thd6 outgroupRootedPhyloGraph)

            -- create list of multi-traversals with original rooting first
            -- subsequent rerooting do not reoptimize exact characters (add nonadd etc) 
            -- they are taken from the fully labelled first input decorated graph later when output graph created
            recursiveRerootList = outgroupRootedPhyloGraph : minimalReRootPhyloGraph outgroupRootedPhyloGraph childrenOfRoot

            minCostRecursive = minimum $ fmap snd6 recursiveRerootList
            minCostGraphListRecursive = filter ((== minCostRecursive).snd6) recursiveRerootList

            -- create optimal final graph with best costs and best traversal (rerooting) forest for each character
            --  travesal for exact charcaters (and costs) are the first of each least since exact onluy optimizaed for that 
            --  traversal graph.  The result has approprotate post-order assignments for traversals, preorder "final" assignments
            -- are propagated to the Decorated graph field after the preorder pass.
            graphWithBestAssignments = setBestGraphAssignments recursiveRerootList (fst6 outgroupRootedPhyloGraph, -1.0, thd6 outgroupRootedPhyloGraph, [], [], six6 outgroupRootedPhyloGraph) 0 0
        in
        --trace ("Outgroup cost:" ++ show (snd6 outgroupRootedPhyloGraph))
        --trace ("Initial Children: " ++ show childrenOfRoot)
        --trace (show minCost ++ ":" ++ (show $ sort $ fmap snd6 rerootedPhyloGraphList) ++ "\n" ++
        --    show minCostDirect ++ ":" ++ (show $ sort $ fmap snd6 rerootPhyloGraphListDirect)
        --    ++ "\n" ++ show minCostRecursive ++ ":" ++ (show $ sort $ fmap snd6 recursiveRerootList))
        graphWithBestAssignments

-- | setBestGraphAssignments take a list of Phylogenetics Graphs and an initial graph with counter as 7 fields
--   and returns teh Phylogenetic graph with optimal exact an dnon-exact costs and travesal greaphs
--   exact character and traversal graphs are taken as those of the head of the list
setBestGraphAssignments :: [PhylogeneticGraph] -> (SimpleGraph, VertexCost, DecoratedGraph, [BlockDisplayForest], [[DecoratedGraph]], V.Vector (V.Vector CharInfo)) -> Int -> Int -> PhylogeneticGraph
setBestGraphAssignments inGraphList (curSimpleGraph, curCost, curDecGraph, curBlockDisplayForestList, curTraversalGraphList, curCharInfo) blockIndex charIndex =
    if null inGraphList then (LG.empty, 0.0, LG.empty, V.empty, V.empty, V.empty)
    else 
        -- Get best costs,  assignments, traversal tree(s) for each non-exact character, head for exact
        -- Get total best cost
        head inGraphList   


-- | minimalReRootPhyloGraph takes an inialtial fully labelled phylogenetic graph
-- and "intelligently" reroots by traversing through adjacent edges, hopefully
-- reoptimizing the minimum number of vertices each time (2) but could be more depending
-- on graph topology
minimalReRootPhyloGraph :: PhylogeneticGraph -> [LG.Node] -> [PhylogeneticGraph]
minimalReRootPhyloGraph inGraph nodesToRoot =
    if null nodesToRoot then []
    else 
        let firstRerootIndex = head nodesToRoot
            nextReroots = (LG.descendants (thd6 inGraph) firstRerootIndex) ++ (tail nodesToRoot)
            newGraph = GO.rerootPhylogeneticGraph' inGraph firstRerootIndex
        in
        --trace ("New cost:" ++ show (snd6 newGraph) ++ " vs " ++ (show $ GO.graphCostFromNodes $ thd6 newGraph))
        newGraph : minimalReRootPhyloGraph newGraph nextReroots

-- | fullyLabelGraph takes an unlabelled "simple' graph, performs post and preorder passes to 
-- fully label the graph and return a PhylogeenticGraph
fullyLabelGraph :: GlobalSettings -> ProcessedData -> DecoratedGraph -> SimpleGraph -> PhylogeneticGraph
fullyLabelGraph inGS inData leafGraph inGraph = 
    if LG.isEmpty inGraph then (LG.empty, 0.0, LG.empty, V.empty, V.empty, V.empty)
    else 
        let postOrderTree = postOrderTreeTraversal inGS inData leafGraph inGraph
            preOrderTree = preOrderTreeTraversal inGS inData postOrderTree
        in
        preOrderTree

-- | makeLeafGraph takes input data and creates a 'graph' of leaves with Vertex informnation
-- but with zero edges.  This 'graph' can be reused as a starting structure for graph construction
-- to avoid remaking of leaf vertices
makeLeafGraph :: ProcessedData -> DecoratedGraph
makeLeafGraph inData@(nameVect, bvNameVect, blocDataVect) =
    if V.null nameVect then error "Empty ProcessedData in makeLeafGraph"
    else 
        let leafVertexList = V.toList $ V.map (makeLeafVertex nameVect bvNameVect blocDataVect) (V.fromList [0.. (V.length nameVect) - 1])
        in
        LG.mkGraph leafVertexList []

-- | makeLeafVertex makes a single unconnected vertex for a leaf
makeLeafVertex :: V.Vector NameText -> V.Vector NameBV -> V.Vector BlockData -> Int -> LG.LNode VertexInfo
makeLeafVertex nameVect bvNameVect inData localIndex =
    let centralData = V.map snd3 inData 
        thisData = V.map (V.! localIndex) centralData
        newVertex = VertexInfo  { index = localIndex
                                , bvLabel = bvNameVect V.! localIndex
                                , parents = V.empty
                                , children = V.empty
                                , nodeType = LeafNode
                                , vertName =  nameVect V.! localIndex
                                , vertData = thisData
                                , vertexCost = 0.0
                                , subGraphCost = 0.0
                                }   
        in
        -- trace (show (length thisData) ++ (show $ fmap length thisData))
        (localIndex, newVertex)



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
        --trace ("It Begins at " ++ show rootIndex) (
        if not $ LG.isRoot inGraph rootIndex then error ("Index "  ++ (show rootIndex) ++ " not root in graph:\n" ++ (GFU.showGraph inGraph))
        else newTree 
        -- )

-- | getVirtualRootEdge cretes tehe virtual edge that would have existed to crete that node rott
-- the ide is that if the tree had been rooted somewhere else -- what edge would have existied and was
-- "divided" to crete this node.  This is used for individual charcater graph traversal foci
-- edge is basically undirected since orientation is unknown
getVirtualRootEdge :: (Show a, Show b) => LG.Gr a b -> LG.Node -> LG.Edge
getVirtualRootEdge inGraph inNode = 
    if LG.isEmpty inGraph then error "Empty graph in getVirtualRootEdge"
    else 
        let childList = LG.descendants inGraph inNode
            parentList = LG.parents inGraph inNode
        in
        -- this really shouldn't be used but in recursion for traversal may be needed
        if LG.isLeaf inGraph inNode then (head parentList, inNode)
        else 
            if length childList /= 2 then error ("Root node with /= 2 childern in getVirtualRootEdge\n" ++ (GFU.showGraph inGraph)) 
            else (head childList, last childList)


-- | postDecorateTree begins at start index (usually root, but could be a subtree) and moves preorder till childrend ane labelled and then reurns postorder
-- labelling vertices and edges as it goes back to root
-- this for a tree so single rtoot
postDecorateTree :: GlobalSettings -> ProcessedData -> SimpleGraph -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.Node -> PhylogeneticGraph
postDecorateTree inGS inData simpleGraph curDecGraph blockCharInfo curNode = 
    -- if node in there nothing to do and return
    if LG.gelem curNode curDecGraph then 
        let nodeLabel = LG.lab curDecGraph curNode
            -- identify/create the virtual edge this node would have been created from
            -- this for traversal focus use for character 
            impliedRootEdge = getVirtualRootEdge curDecGraph curNode
        in
        if nodeLabel == Nothing then error ("Null label for node " ++ show curNode)
        else (simpleGraph, subGraphCost (fromJust nodeLabel), curDecGraph, V.empty, V.replicate (length blockCharInfo) (V.singleton curDecGraph), blockCharInfo)

    -- Need to make node
    else 
        
        -- check if children in graph
        let nodeChildren = LG.descendants simpleGraph curNode  -- should be 1 or 2, not zero since all leaves already in graph
            leftChild = (head nodeChildren)
            rightChild = (last nodeChildren)
            leftChildTree = postDecorateTree inGS inData simpleGraph curDecGraph blockCharInfo leftChild
            rightLeftChildTree = if (length nodeChildren == 2) then postDecorateTree inGS inData simpleGraph (thd6 $ leftChildTree) blockCharInfo rightChild
                                 else leftChildTree
        in 

        if length nodeChildren > 2 then error ("Graph not dichotomous in postDecorateTree node " ++ (show curNode) ++ "\n" ++ LG.prettify simpleGraph)
        else if length nodeChildren == 0 then error ("Leaf not in graph in postDecorateTree node " ++ (show curNode) ++ "\n" ++ LG.prettify simpleGraph)

        -- make node from childern
        else    
            -- make node from children and new edges to children
            -- median2 takes characters in blocks--but for tree really all same block
            let newSubTree = thd6 $ rightLeftChildTree
                leftChildLabel = fromJust $ LG.lab newSubTree leftChild
                rightChildLabel = fromJust $ LG.lab newSubTree rightChild
                
                newCharData = GO.createVertexDataOverBlocks  (vertData leftChildLabel) (vertData  rightChildLabel) blockCharInfo []
                newCost =  V.sum $ V.map (V.sum) $ V.map (V.map snd) newCharData
                newVertex = VertexInfo {  index = curNode
                                        , bvLabel = (bvLabel leftChildLabel) .|. (bvLabel rightChildLabel)
                                        , parents = V.fromList $ LG.parents simpleGraph curNode
                                        , children = V.fromList nodeChildren
                                        , nodeType = GO.getNodeType simpleGraph curNode
                                        , vertName = T.pack $ "HTU" ++ (show curNode)
                                        , vertData = V.map (V.map fst) newCharData
                                        , vertexCost = newCost
                                        , subGraphCost = (subGraphCost leftChildLabel) + (subGraphCost rightChildLabel) + newCost
                                        }   
                newEdgesLabel = EdgeInfo {    minLength = newCost / 2.0
                                            , maxLength = newCost / 2.0
                                            , edgeType = TreeEdge
                                         }
                newEdges = fmap LG.toEdge $ LG.out simpleGraph curNode 
                newLEdges =  fmap (LG.toLEdge' newEdgesLabel) newEdges
                newGraph =  LG.insEdges newLEdges $ LG.insNode (curNode, newVertex) newSubTree                    
            in
            -- return new graph
            --trace ("Node stuff: " ++ show (curNode, index leftChildLabel, index rightChildLabel)) -- , vertData leftChildLabel, vertData  rightChildLabel, V.map (V.map fst) newCharData))
            --trace ("Chars: " ++ (show $ V.length  $ vertData $ fromJust $ LG.lab newSubTree leftChild) ++ " " ++ (show $ fmap V.length $ vertData $ fromJust $ LG.lab newSubTree leftChild)
               -- ++ (show $ V.map (V.map snd) newCharData)) (
            -- trace ("Node " ++ (show curNode) ++ " Cost: " ++ (show $ subGraphCost newVertex) ++ " " 
              --  ++ show ((subGraphCost $ fromJust $ LG.lab newSubTree leftChild), (subGraphCost $ fromJust $ LG.lab newSubTree rightChild), newCost))
            (simpleGraph, (subGraphCost newVertex), newGraph, V.empty, V.empty, blockCharInfo)
            --)
            

-- | preOrderTreeTraversal takes a preliminarily labelled PhylogeneticGraph
-- and returns a full labbels with 'final' assignments
-- invafiant that root is HTU !! nLeaves?
preOrderTreeTraversal :: GlobalSettings -> ProcessedData -> PhylogeneticGraph -> PhylogeneticGraph
preOrderTreeTraversal inGS inData inPGraph = 
    if LG.isEmpty (thd6 inPGraph) then (fst6 inPGraph, 0.0, LG.empty, V.empty, V.empty, V.empty)
    else 
        inPGraph
        
