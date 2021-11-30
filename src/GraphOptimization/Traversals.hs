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
                                    , postOrderSoftWiredTraversal
                                    , multiTraverseFullyLabelTree
                                    , multiTraverseFullyLabelGraph
                                    , multiTraverseFullyLabelSoftWired
                                    , checkUnusedEdgesPruneInfty
                                    , makeLeafGraph
                                    , makeSimpleLeafGraph
                                    , makeLeafGraphSoftWired
                                    , postDecorateTree
                                    , postDecorateTree'
                                    , postDecorateSoftWired
                                    , postDecorateSoftWired'
                                    , updateAndFinalizePostOrderSoftWired
                                    , updatePhylogeneticGraphCost
                                    , getW15NetPenalty
                                    , getW15RootCost
                                    , generalizedGraphPostOrderTraversal
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
-- import Debug.Debug as D
import Utilities.Utilities as U
import Debug.Trace
import qualified Data.BitVector.LittleEndian as BV



-- | multiTraverseFullyLabelGraph is a wrapper around multi-traversal functions for Tree, 
-- Soft-wired network graph, and Hard-wired network graph
-- can either find root, be given root, or start somwhere else (startVertex) do optimize only a component of a forest
multiTraverseFullyLabelGraph :: GlobalSettings -> ProcessedData -> Bool -> Bool -> Maybe Int -> SimpleGraph -> PhylogeneticGraph
multiTraverseFullyLabelGraph inGS inData pruneEdges warnPruneEdges startVertex inGraph
  | LG.isEmpty inGraph = emptyPhylogeneticGraph
  | graphType inGS == Tree =
    -- test for Tree
    let (_, _, _, networkVertexList) = LG.splitVertexList inGraph
    in
    if null networkVertexList then 
        let leafGraph = makeLeafGraph inData
        in multiTraverseFullyLabelTree inGS inData leafGraph startVertex inGraph
    else errorWithoutStackTrace "Input graph is not a tree/forest, but graph type has been specified (perhaps by default) as Tree. Modify input graph or use 'set()' command to specify network type"
      | graphType inGS == SoftWired = 
        let leafGraph = makeLeafGraphSoftWired inData
        in multiTraverseFullyLabelSoftWired  inGS inData pruneEdges warnPruneEdges leafGraph startVertex inGraph
      | graphType inGS == HardWired = errorWithoutStackTrace "Hard-wired graph optimization not yet supported"
      | otherwise = errorWithoutStackTrace ("Unknown graph type specified: " ++ show (graphType inGS))

-- | multiTraverseFullyLabelSoftWired fully labels a softwired network component forest
-- including traversal rootings-- does not reroot on network edges
-- allows indegree=outdegree=1 vertices
-- pruneEdges and warnPruneEdges specify if unused edges (ie not in diuaplytrees) are pruned from
-- canonical tree or if an infinity cost is returned and if a trace warning is thrown if so.
-- in general--input trees should use "pruneEdges" during search--not
-- can either find root, be given root, or start somwhere else (startVertex) do optimize only a component of a forest
multiTraverseFullyLabelSoftWired :: GlobalSettings -> ProcessedData -> Bool -> Bool -> DecoratedGraph -> Maybe Int -> SimpleGraph -> PhylogeneticGraph
multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex inSimpleGraph =
    if LG.isEmpty inSimpleGraph then emptyPhylogeneticGraph
    else 
        let nonExactChars = U.getNumberNonExactCharacters (thd3 inData)
            (postOrderGraph, localRootCost, localStartVertex) = generalizedGraphPostOrderTraversal inGS nonExactChars inData leafGraph startVertex inSimpleGraph
            fullyOptimizedGraph = PRE.preOrderTreeTraversal inGS (finalAssignment inGS) (nonExactChars > 0) localStartVertex postOrderGraph
            in
            checkUnusedEdgesPruneInfty inGS inData pruneEdges warnPruneEdges leafGraph $ updatePhylogeneticGraphCost fullyOptimizedGraph (localRootCost + (snd6 fullyOptimizedGraph))

-- | multiTraverseFullyLabelTree performs potorder on default root and other traversal foci, taking the minimum 
-- traversal cost for all nonexact charcters--the initial rooting is used for exact characters 
-- operates with Tree functions
-- need to add forest functionality--in principle just split into components and optimize them independently
-- but get into root index issues the way this is written now. 
-- can either find root, be given root, or start somwhere else (startVertex) do optimize only a component of a forest
multiTraverseFullyLabelTree :: GlobalSettings -> ProcessedData -> DecoratedGraph -> Maybe Int -> SimpleGraph -> PhylogeneticGraph
multiTraverseFullyLabelTree inGS inData leafGraph startVertex inSimpleGraph =
    if LG.isEmpty inSimpleGraph then emptyPhylogeneticGraph
    else
        let nonExactChars = U.getNumberNonExactCharacters (thd3 inData)
            (postOrderGraph, _, localStartVertex) = generalizedGraphPostOrderTraversal inGS nonExactChars inData leafGraph startVertex inSimpleGraph
        in
        PRE.preOrderTreeTraversal inGS (finalAssignment inGS) (nonExactChars > 0) localStartVertex postOrderGraph
            
    
-- | generalizedGraphPostOrderTraversal performs the postorder pass 
-- on a graph (tree, softWired, or hardWired) to determine the "preliminary" character states
-- include penalty factor cost but not root cost which may or may not be wanted depending on context
-- if full graph--yes, if a component yes or no.
-- hence returnde das pair
generalizedGraphPostOrderTraversal :: GlobalSettings -> Int -> ProcessedData -> DecoratedGraph -> Maybe Int -> SimpleGraph -> (PhylogeneticGraph, VertexCost, Int)
generalizedGraphPostOrderTraversal inGS nonExactChars inData leafGraph startVertex inSimpleGraph =

    -- select postOrder function based on graph type 
    let postOrderFunction = if (graphType inGS) == Tree then postOrderTreeTraversal
                         else if (graphType inGS) == SoftWired then postOrderSoftWiredTraversal
                         else error ("Graph type not implemented: " ++ (show $ graphType inGS))

        -- first traversal on outgroup roo
        outgroupRooted = if startVertex == Nothing then postOrderFunction inGS inData leafGraph startVertex $ GO.rerootTree' inSimpleGraph (outgroupIndex inGS)
                        else postOrderFunction inGS inData leafGraph startVertex inSimpleGraph

        -- start at start vertex--for components or ur-root for full graph
        startVertexList = if startVertex == Nothing then fmap fst $ LG.getRoots $ thd6 outgroupRooted
                              else [fromJust startVertex]

        -- next edges (to vertex in list) to perform rerroting
        -- progresses recursivey over adjacent edges to minimize node reoptimization
        childrenOfRoot = concatMap (LG.descendants (thd6 outgroupRooted)) startVertexList
        grandChildrenOfRoot = concatMap (LG.descendants (thd6 outgroupRooted)) childrenOfRoot

        -- create list of multi-traversals with original rooting first
        -- subsequent rerooting do not reoptimize exact characters (add nonadd etc) 
        -- they are taken from the fully labelled first input decorated graph later when output graph created
        -- it is important that the first graph be the ourgroup rooted graph (outgroupRootedPhyloGraph) so this 
        -- will have the preoder assignmentsd for th eoutgroup rooted graph as 3rd field.  This can be used for incremental
        -- optimization to get O(log n) initial postorder assingment when mutsating graph.
        recursiveRerootList = outgroupRooted : minimalReRootPhyloGraph inGS (graphType inGS) outgroupRooted (head startVertexList) grandChildrenOfRoot

        -- perform traceback on resolution caches is graphtype = softWired
        recursiveRerootList' = if (graphType inGS) == Tree then recursiveRerootList
                               else if (graphType inGS) == SoftWired then fmap (updateAndFinalizePostOrderSoftWired startVertex (head startVertexList)) recursiveRerootList
                               else error ("Graph type not implemented: " ++ (show $ graphType inGS))


        finalizedPostOrderGraphList = L.sortOn snd6 recursiveRerootList'

        -- create optimal final graph with best costs and best traversal (rerooting) forest for each character
        -- traversal for exact characters (and costs) are the first of each least since exact only optimizaed for that 
        -- traversal graph.  The result has approprotate post-order assignments for traversals, preorder "final" assignments
        -- are propagated to the Decorated graph field after the preorder pass.
        -- doesn't have to be sorted, but should minimize assignments
        graphWithBestAssignments = L.foldl1' setBetterGraphAssignment finalizedPostOrderGraphList

        -- same root cost if same data and number of roots
        localRootCost = if (rootCost inGS) == NoRootCost then 0.0
                        else if (rootCost inGS) == Wheeler2015Root then getW15RootCost inData outgroupRooted
                        else error ("Root cost type " ++ (show $ rootCost inGS) ++ " is not yet implemented")

    in
    -- only static characters
    if nonExactChars == 0 then 
        let penaltyFactor  = if (graphType inGS == Tree) then 0.0
                             else if (graphType inGS == HardWired) then error ("Graph type not implemented: " ++ (show $ graphType inGS))
                             else if (graphFactor inGS) == NoNetworkPenalty then 0.0
                             else if (graphFactor inGS) == Wheeler2015Network then getW15NetPenalty outgroupRooted
                             else error ("Network penalty type " ++ (show $ graphFactor inGS) ++ " is not yet implemented")

            outgroupRooted' = updatePhylogeneticGraphCost outgroupRooted (penaltyFactor + (snd6 outgroupRooted))
        in
        (outgroupRooted', localRootCost, head startVertexList)
    else if nonExactChars == 1 then
        let penaltyFactorList  = if (graphType inGS == Tree) then replicate (length finalizedPostOrderGraphList) 0.0
                                 else if (graphType inGS == HardWired) then error ("Graph type not implemented: " ++ (show $ graphType inGS))
                                 else if (graphFactor inGS) == NoNetworkPenalty then replicate (length finalizedPostOrderGraphList) 0.0
                                 else if (graphFactor inGS) == Wheeler2015Network then fmap getW15NetPenalty finalizedPostOrderGraphList
                                 else error ("Network penalty type " ++ (show $ graphFactor inGS) ++ " is not yet implemented")
            newCostList = zipWith (+) penaltyFactorList (fmap snd6 finalizedPostOrderGraphList)

            finalizedPostOrderGraph = head $ L.sortOn snd6 $ zipWith updatePhylogeneticGraphCost finalizedPostOrderGraphList newCostList
        in
        (finalizedPostOrderGraph, localRootCost, head startVertexList)
    -- multiple dynamic characters--cjecks for best root for each character
    else 
        let penaltyFactor  = if (graphType inGS == Tree) then 0.0
                             else if (graphType inGS == HardWired) then error ("Graph type not implemented: " ++ (show $ graphType inGS))
                             else if (graphFactor inGS) == NoNetworkPenalty then 0.0
                             else if (graphFactor inGS) == Wheeler2015Network then getW15NetPenalty graphWithBestAssignments
                             else error ("Network penalty type " ++ (show $ graphFactor inGS) ++ " is not yet implemented")

            graphWithBestAssignments' = updatePhylogeneticGraphCost graphWithBestAssignments (penaltyFactor + (snd6 graphWithBestAssignments))
        in
        (graphWithBestAssignments', localRootCost, head startVertexList)

               


-- | updatePhylogeneticGraphCost takes a PhylgeneticGrtaph and Double and replaces the cost (snd of 6 fields)
-- and returns Phylogenetic graph
updatePhylogeneticGraphCost :: PhylogeneticGraph -> VertexCost -> PhylogeneticGraph
updatePhylogeneticGraphCost (a, _, b, c, d, e) newCost = (a, newCost, b, c, d, e)

-- | getW15RootCost creates a root cost as the 'insertion' of character data.  For sequence data averaged over
-- leaf taxa
getW15RootCost :: ProcessedData -> PhylogeneticGraph -> VertexCost
getW15RootCost (_, _, blockDataV) inGraph =
    if LG.isEmpty $ thd6 inGraph then 0.0
    else
        let (rootList, leafList, _, _) = LG.splitVertexList $ fst6 inGraph
            numRoots = length rootList
            numLeaves = length leafList
            insertDataCost = V.sum $ fmap U.getblockInsertDataCost blockDataV
        in
        (fromIntegral numRoots) * insertDataCost /  (fromIntegral numLeaves)

-- | getW15NetPenalty takes a Phylogenetic tree and returns the network penalty of Wheeler (2015)
-- modified to take theg union of all edges of trees of minimal length
getW15NetPenalty :: PhylogeneticGraph -> VertexCost
getW15NetPenalty inGraph =
    if LG.isEmpty $ thd6 inGraph then 0.0
    else
        let (bestTreeList, _) = extractLowestCostDisplayTree inGraph
            bestTreesEdgeList = L.nubBy undirectedEdgeEquality $ concat $ fmap LG.edges bestTreeList
            rootIndex = fst $ head $ LG.getRoots (fst6 inGraph)
            blockPenaltyList = fmap (getBlockW2015 bestTreesEdgeList rootIndex) (fth6 inGraph)
            (_, leafList, _, _) = LG.splitVertexList (fst6 inGraph)
            numLeaves = length leafList
            divisor = 4.0 * (fromIntegral numLeaves) - 4.0
        in
        (sum $ blockPenaltyList) / divisor


-- | getBlockW2015 takes teh list of trees for a block, gets teh root cost and determines the individual
-- penlaty cost of that block
getBlockW2015 :: [LG.Edge] -> Int -> [BlockDisplayForest] -> VertexCost
getBlockW2015 treeEdgeList rootIndex blockTreeList =
    if null treeEdgeList || null blockTreeList then 0.0
    else 
        let blockTreeEdgeList = L.nubBy undirectedEdgeEquality $ concatMap LG.edges blockTreeList
            numExtraEdges = length $ undirectedEdgeMinus blockTreeEdgeList treeEdgeList
            blockCost = vertexCost $ fromJust $ LG.lab (head blockTreeList) rootIndex
        in
        blockCost * (fromIntegral numExtraEdges) 

-- | checkUnusedEdgesPruneInfty checks if a softwired phylogenetic graph has 
-- "unused" edges sensu Wheeler 2015--that an edge in the canonical graph is 
-- not present in any of the block display trees (that are heursticlly optimal)
-- the optrions specify if the cost returned is Infinity (really max boun Double)
-- with no pruning of edges or th ecost is left unchanged and unused edges are
-- pruned from the canonical graph
-- this is unDirected due to rerooting heuristic in post/preorder optimization
-- inifinity defined in Types.hs
checkUnusedEdgesPruneInfty :: GlobalSettings -> ProcessedData -> Bool -> Bool -> DecoratedGraph-> PhylogeneticGraph -> PhylogeneticGraph
checkUnusedEdgesPruneInfty inGS inData pruneEdges warnPruneEdges leafGraph inGraph@(inSimple, _, inCanonical, blockTreeV, charTreeVV, charInfoVV) =
    let simpleEdgeList = LG.edges inSimple
        displayEdgeSet =  L.nubBy undirectedEdgeEquality $ concat $ concat $ fmap (fmap LG.edges) blockTreeV
        unusedEdges = undirectedEdgeMinus simpleEdgeList displayEdgeSet
    in
    -- no unused edges all OK
    if null unusedEdges then inGraph

    -- unused edges--do not prune return "infinite cost"
    else if not pruneEdges then (inSimple, infinity, inCanonical, blockTreeV, charTreeVV, charInfoVV) 

    -- unused but pruned--need to prune nodes and reoptimize to get final assignments correct
    else 
        let newSimpleGraph = LG.delEdges unusedEdges inSimple
            contractedSimple = LG.contractIn1Out1Edges newSimpleGraph
        in
        if warnPruneEdges then 
            trace ("Pruning " ++ (show $ length unusedEdges) ++ " unused edges and reoptimizing graph") 
            multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph Nothing contractedSimple 
            
        else multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph Nothing contractedSimple 

-- | undirectedEdgeEquality checks edgse for equality irrespective of direction
undirectedEdgeEquality :: LG.Edge -> LG.Edge -> Bool
undirectedEdgeEquality (a,b) (c,d) = if a == c && b == d then True
                                               else if a == d && b == c then True
                                               else False 

-- | undirectedEdgeMinus subtracts edges in the second list from those in the first using 
-- undirected matching
undirectedEdgeMinus :: [LG.Edge] -> [LG.Edge] -> [LG.Edge]
undirectedEdgeMinus firstList secondList =
    if null firstList then []
    else 
        let firstEdge@(a,b) = head firstList
        in
        if firstEdge `L.elem` secondList then undirectedEdgeMinus (tail firstList) secondList
        else if (b,a) `L.elem` secondList then undirectedEdgeMinus (tail firstList) secondList
        else firstEdge : undirectedEdgeMinus (tail firstList) secondList

-- | extractLowestCostDisplayTree takes a phylogenetic graph and takes all valid (complete) resolutions 
-- (display trees) and their costs
-- and determines the total cost (over all blocks) of each display tree 
-- the lowest cost display tree(s) as list are returned with cost
-- this is used in Wheeler (2015) network penalty
extractLowestCostDisplayTree :: PhylogeneticGraph -> ([BlockDisplayForest], VertexCost) 
extractLowestCostDisplayTree inGraph =
 if LG.isEmpty $ thd6 inGraph then error "Empty graph in extractLowestCostDisplayTree" 
 else 
    let (_, outgroupRootLabel) =  head $ LG.getRoots (thd6 inGraph)
        blockResolutionLL = V.toList $ fmap PO.getAllResolutionList (vertexResolutionData outgroupRootLabel)
        displayTreeBlockList = L.transpose blockResolutionLL
        displayTreePairList = L.foldl1' sumTreeCostLists displayTreeBlockList
        minimumCost = minimum $ fmap snd displayTreePairList
        (bestDisplayTreeList, _) = unzip $ filter ((== minimumCost) . snd) displayTreePairList
    in
    -- trace ("FC: " ++ (show $ fmap snd displayTreePairList)) 
    (bestDisplayTreeList, minimumCost)

-- | sumTreeCostLists takes two lists of (Graph, Cost) pairs and sums the costs and keeps the trees the same
-- does not check that graphs are the same after debug
sumTreeCostLists :: (Eq a, Eq b) => [(LG.Gr a b, VertexCost)] ->  [(LG.Gr a b, VertexCost)] ->  [(LG.Gr a b, VertexCost)] 
sumTreeCostLists firstList secondList = 
    if null firstList || null secondList then error "Empty list in sumTreeCostLists"
    else 
        let (firstGraphList, firstCostList) = unzip firstList
            (secondGraphList, secondCostList) =  unzip secondList
            newCostList = zipWith (+)  firstCostList secondCostList 

            -- remove once working
            checkList = filter (== False) $ zipWith LG.equal firstGraphList secondGraphList
        in
        if null checkList then error ("Graph lists not same : " ++ (show checkList))
        else trace ("Graphs match ") zip firstGraphList newCostList


-- | updateAndFinalizePostOrderSoftWired performs the pre-order traceback on the resolutions to create the correct vertex states,
-- ports the post order assignments to the canonical tree and display trees, and creates the character trees from the block trees
updateAndFinalizePostOrderSoftWired :: Maybe Int -> Int -> PhylogeneticGraph -> PhylogeneticGraph
updateAndFinalizePostOrderSoftWired startVertex rootIndex inGraph =
    if LG.isEmpty $ thd6 inGraph then inGraph
    else
        let outgroupRootLabel =  fromJust $ LG.lab (thd6 inGraph) rootIndex
            (displayGraphVL, lDisplayCost) = PO.extractDisplayTrees startVertex True (vertexResolutionData outgroupRootLabel)

            -- traceback on resolutions
            newGraph = softWiredPostOrderTraceBack rootIndex (thd6 inGraph)

            -- propagate updated post-order assignments to display trees, which updates dispay tree and cost ealier
            displayGraphVL' = V.zipWith (assignPostOrderToDisplayTree (fmap (vertData . snd) (LG.labNodes newGraph) )) displayGraphVL (V.fromList [0..(V.length displayGraphVL - 1)])

            -- create new, fully  updated post-order graph
            finalPreOrderGraph = (fst6 inGraph, lDisplayCost, newGraph, displayGraphVL', PO.divideDecoratedGraphByBlockAndCharacterSoftWired displayGraphVL', six6 inGraph)
        in
        finalPreOrderGraph


-- | postOrderSoftWiredTraversal performs postorder traversal on Soft-wired graph
postOrderSoftWiredTraversal :: GlobalSettings -> ProcessedData -> DecoratedGraph -> Maybe Int -> SimpleGraph -> PhylogeneticGraph
postOrderSoftWiredTraversal inGS inData@(_, _, blockDataVect) leafGraph startVertex inSimpleGraph =
    if LG.isEmpty inSimpleGraph then emptyPhylogeneticGraph
    else
         -- Assumes root is Number of Leaves  
        let rootIndex = if startVertex == Nothing then V.length $ fst3 inData
                        else fromJust startVertex
            blockCharInfo = V.map thd3 blockDataVect
            newSoftWired = postDecorateSoftWired inGS inSimpleGraph leafGraph blockCharInfo rootIndex
        in
        --trace ("It Begins at " ++ show rootIndex) (
        if not $ LG.isRoot inSimpleGraph rootIndex then
            let localRootList = fst <$> LG.getRoots inSimpleGraph
                localRootEdges = concatMap (LG.out inSimpleGraph) localRootList
                currentRootEdges = LG.out inSimpleGraph rootIndex
            in
            error ("Index "  ++ show rootIndex ++ " with edges " ++ show currentRootEdges ++ " not root in graph:" ++ show localRootList ++ " edges:" ++ show localRootEdges ++ "\n" ++ GFU.showGraph inSimpleGraph)
        else newSoftWired
        -- )

-- | postDecorateSoftWired' wrapper for postDecorateSoftWired with args in differnt order for mapping
postDecorateSoftWired' :: GlobalSettings -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.Node -> SimpleGraph -> PhylogeneticGraph
postDecorateSoftWired' inGS curDecGraph blockCharInfo curNode simpleGraph = postDecorateSoftWired inGS simpleGraph curDecGraph blockCharInfo curNode

-- | postDecorateSoftWired begins at start index (usually root, but could be a subtree) and moves preorder till children are labelled 
-- and then recurses to root postorder labelling vertices and edges as it goes
-- this for a single root
postDecorateSoftWired :: GlobalSettings -> SimpleGraph -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.Node -> PhylogeneticGraph
postDecorateSoftWired inGS simpleGraph curDecGraph blockCharInfo curNode =
    -- if node in there nothing to do and return
    if LG.gelem curNode curDecGraph then
        let nodeLabel = LG.lab curDecGraph curNode
        in
        if isNothing nodeLabel then error ("Null label for node " ++ show curNode)
        else (simpleGraph, subGraphCost (fromJust nodeLabel), curDecGraph, mempty, mempty, blockCharInfo)

    else
        -- trace ("PDSW making node " ++ show curNode ++ " in\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph curDecGraph)) (
        let nodeChildren = LG.descendants simpleGraph curNode  -- should be 1 or 2, not zero since all leaves already in graph
            leftChild = head nodeChildren
            rightChild = last nodeChildren
            leftChildTree = postDecorateSoftWired inGS simpleGraph curDecGraph blockCharInfo leftChild
            rightLeftChildTree = if length nodeChildren == 2 then postDecorateSoftWired inGS simpleGraph (thd6 leftChildTree) blockCharInfo rightChild
                                 else leftChildTree
        in
        -- Checks on children
        if length nodeChildren > 2 then error ("Graph not dichotomous in postDecorateSoftWired node " ++ show curNode ++ "\n" ++ LG.prettify simpleGraph)
        else if null nodeChildren then error ("Leaf not in graph in postDecorateSoftWired node " ++ show curNode ++ "\n" ++ LG.prettify simpleGraph)

        else
            -- make node from child block resolutions
            let newSubTree = thd6 rightLeftChildTree
            in

            -- single child of node (can certinly happen with soft-wired networks
            if length nodeChildren == 1 then
                -- trace ("Outdegree 1: " ++ (show curNode) ++ " " ++ (show $ GO.getNodeType simpleGraph curNode) ++ " Child: " ++ (show nodeChildren)) (
                let (newGraph, _, _, _, _) = PO.getOutDegree1VertexAndGraph curNode (fromJust $ LG.lab newSubTree leftChild) simpleGraph nodeChildren newSubTree
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
                    resolutionBlockVL = V.zipWith3 (PO.createBlockResolutions (compressResolutions inGS) curNode leftChild' rightChild' leftChildNodeType rightChildNodeType (GO.getNodeType simpleGraph curNode)) (vertexResolutionData leftChildLabel) (vertexResolutionData rightChildLabel) blockCharInfo

                    -- create canonical Decorated Graph vertex
                    -- 0 cost becasue can't know cosrt until hit root and get best valid resolutions
                    newVertexLabel = VertexInfo {  index = curNode
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

                    -- (displayGraphVL, lDisplayCost) = if (nodeType newVertexLabel) == RootNode then PO.extractDisplayTrees True resolutionBlockVL
                                                     -- else (mempty, 0.0)

                in
                (simpleGraph, 0.0, newGraph, mempty, mempty, blockCharInfo)

                {-
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
                -}

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
        else if (nodeType (snd firstNode) == NetworkNode) && isJust (L.find ((== fst firstNode). fst) updatedNodes) then
            softWiredPrelimTraceback inGraph (tail nodesToUpdate) updatedNodes

        else
            let -- get children 
                firstChildren = LG.labDescendants inGraph firstNode
                firstChildrenBV = fmap (bvLabel . snd) firstChildren
                firstChildrenIsLeft
                  | length firstChildrenBV == 1 = [True]
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
            leftRightIndexVect = fmap (head . childResolutions) resolutionsByBlockV
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
            displayPopList = fmap (complement . displayBVLabel) inResBlockData



            -- subgraph cost
            displayCostList = fmap displayCost inResBlockData

            -- resolution local cost
            resolutionCostList = fmap resolutionCost inResBlockData

            -- takes only first resolution index pair
            childResolutionList = fmap (head . childResolutions) inResBlockData

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
    -- trace ("SBGA:" ++  (show $ (length  fTraversal, length sTraversal))) (
    if LG.isEmpty fDecGraph then secondGraph
    else if LG.isEmpty sDecGraph then firstGraph
    else
        --trace ("setBetter (" ++ (show fCost) ++ "," ++ (show sCost) ++ ")" ) ( -- ++ " CharInfo blocks:" ++ (show $ length fCharInfo) ++ " characters: " ++ (show $ fmap length fCharInfo) ++ " " ++ (show $ fmap (fmap name) fCharInfo)) (
        let (mergedBlockVect, costVector) = V.unzip $ V.zipWith makeBetterBlock fTraversal sTraversal
        in
         --trace ("setBetter (" ++ (show fCost) ++ "," ++ (show sCost) ++ ") ->" ++ (show $ V.sum costVector) ++ " nt:" ++ (show $ length fTraversal)
         --   ++ "length blocks " ++ (show $ fmap length fTraversal))
        (fSimple, V.sum costVector, fDecGraph, fBlockDisplay, mergedBlockVect, fCharInfo)
        -- )

-- | makeBetterBlocktakes two verctors of traversals. Each vector contains a decorated graph (=traversla graph) for each 
-- character.  This can be a single sequence or series of exact characters
-- the function applies a character cost comparison to get the better 
makeBetterBlock :: V.Vector DecoratedGraph -> V.Vector DecoratedGraph -> (V.Vector DecoratedGraph, VertexCost)
makeBetterBlock firstBlockGraphVect secondBlockGraphVect =
    let (mergedCharacterVect, costVector) = V.unzip $ V.zipWith chooseBetterCharacter firstBlockGraphVect secondBlockGraphVect
    in
    -- trace ("MBB: " ++ (show $ (length  firstBlockGraphVect, length firstBlockGraphVect)))
    (mergedCharacterVect, V.sum costVector)

-- | chooseBetterCharacter takes a pair of character decorated graphs and chooses teh "better" one as in lower cost, or non-zero cost
--  if one is zer (for exact characters) and returns the better character and cost
--  graph can have multiple roots
chooseBetterCharacter :: DecoratedGraph -> DecoratedGraph -> (DecoratedGraph, VertexCost)
chooseBetterCharacter firstGraph secondGraph
  | LG.isEmpty firstGraph = error "Empty first graph in chooseBetterCharacter"
  | LG.isEmpty secondGraph = error "Empty second graph in chooseBetterCharacter"
  | otherwise =
    let firstGraphCost = sum $ fmap (subGraphCost . snd) (LG.getRoots firstGraph)
        secondGraphCost = sum $ fmap (subGraphCost . snd) (LG.getRoots secondGraph)
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
minimalReRootPhyloGraph :: GlobalSettings -> GraphType -> PhylogeneticGraph -> Int -> [LG.Node] -> [PhylogeneticGraph]
minimalReRootPhyloGraph inGS localGraphType inGraph originalRoot nodesToRoot =
    -- trace ("MRR: " ++ (show nodesToRoot) ++ " " ++ (show $ fmap (LG.descendants (thd6 inGraph)) nodesToRoot)) (
    if null nodesToRoot then []
    else
        let firstRerootIndex = head nodesToRoot
            nextReroots = LG.descendants (thd6 inGraph) firstRerootIndex ++ tail nodesToRoot
            newGraph
              | localGraphType == Tree = PO.rerootPhylogeneticGraph' inGS Tree False False inGraph originalRoot firstRerootIndex
              | localGraphType == SoftWired = PO.rerootPhylogeneticNetwork' inGS inGraph originalRoot firstRerootIndex
              | otherwise = errorWithoutStackTrace ("Graph type not implemented/recognized: " ++ show localGraphType)
        in
        -- trace ("NRR: " ++ " " ++ (show (LG.descendants (thd6 inGraph) firstRerootIndex)) ) ( -- ++ " -> " ++ (show nextReroots) ++ "\n" ++ (LG.prettify $ fst6 inGraph) ++ "\n" ++ (LG.prettify $ fst6 newGraph)) (
        if fst6 newGraph == LG.empty then minimalReRootPhyloGraph inGS localGraphType inGraph originalRoot nextReroots
        else newGraph : minimalReRootPhyloGraph inGS localGraphType newGraph originalRoot nextReroots
        -- ) -- )


-- | makeLeafGraph takes input data and creates a 'graph' of leaves with Vertex informnation
-- but with zero edges.  This 'graph' can be reused as a starting structure for graph construction
-- to avoid remaking of leaf vertices
makeLeafGraph :: ProcessedData -> DecoratedGraph
makeLeafGraph (nameVect, bvNameVect, blocDataVect) =
    if V.null nameVect then error "Empty ProcessedData in makeLeafGraph"
    else
        let leafVertexList = V.toList $ V.map (makeLeafVertex nameVect bvNameVect blocDataVect) (V.fromList [0.. V.length nameVect - 1])
        in
        LG.mkGraph leafVertexList []


-- | makeSimpleLeafGraph takes input data and creates a 'graph' of leaves with Vertex informnation
-- but with zero edges.  This 'graph' can be reused as a starting structure for graph construction
-- to avoid remaking of leaf vertices
makeSimpleLeafGraph :: ProcessedData -> SimpleGraph
makeSimpleLeafGraph (nameVect, _, _) =
    if V.null nameVect then error "Empty ProcessedData in makeSimpleLeafGraph"
    else
        let leafVertexList = V.toList $ V.map (makeSimpleLeafVertex nameVect) (V.fromList [0.. V.length nameVect - 1])
        in
        LG.mkGraph leafVertexList []
        where makeSimpleLeafVertex a b = (b, a V.! b)


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
postOrderTreeTraversal :: GlobalSettings ->  ProcessedData -> DecoratedGraph -> Maybe Int -> SimpleGraph -> PhylogeneticGraph
postOrderTreeTraversal inGS (_, _, blockDataVect) leafGraph startVertex inGraph  =
    if LG.isEmpty inGraph then emptyPhylogeneticGraph
    else
        -- Assumes root is Number of Leaves  
        let rootIndex = if startVertex == Nothing then  fst $ head $ LG.getRoots inGraph
                        else fromJust startVertex
            blockCharInfo = V.map thd3 blockDataVect
            newTree = postDecorateTree inGraph leafGraph blockCharInfo rootIndex
        in
        -- trace ("It Begins at " ++ (show $ fmap fst $ LG.getRoots inGraph) ++ "\n" ++ show inGraph) (
        if not $ LG.isRoot inGraph rootIndex then
            let localRootList = fst <$> LG.getRoots inGraph
                localRootEdges = concatMap (LG.out inGraph) localRootList
                currentRootEdges = LG.out inGraph rootIndex
            in
            error ("Index "  ++ show rootIndex ++ " with edges " ++ show currentRootEdges ++ " not root in graph:" ++ show localRootList ++ " edges:" ++ show localRootEdges ++ "\n" ++ GFU.showGraph inGraph)
        else newTree
        --)

-- | postDecorateTree' is wrapper for postDecorateTree to alow for mapping
postDecorateTree' :: DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.Node -> SimpleGraph -> PhylogeneticGraph
postDecorateTree' curDecGraph blockCharInfo curNode simpleGraph = postDecorateTree simpleGraph curDecGraph blockCharInfo curNode 

-- | postDecorateTree begins at start index (usually root, but could be a subtree) and moves preorder till children are labelled and then returns postorder
-- labelling vertices and edges as it goes back to root
-- this for a tree so single root
postDecorateTree :: SimpleGraph -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.Node -> PhylogeneticGraph
postDecorateTree simpleGraph curDecGraph blockCharInfo curNode =
    -- if node in there (leaf) nothing to do and return
    if LG.gelem curNode curDecGraph then
        let nodeLabel = LG.lab curDecGraph curNode
            -- identify/create the virtual edge this node would have been created from
            -- this for traversal focus use for character 
            -- impliedRootEdge = getVirtualRootEdge curDecGraph curNode
        in
        if isNothing nodeLabel then error ("Null label for node " ++ show curNode)
        else
            -- trace ("In graph :" ++ (show curNode) ++ " " ++ (show nodeLabel))
            --  replicating curaDecGraph with number opf blocks--but all the same for tree 
            (simpleGraph, subGraphCost (fromJust nodeLabel), curDecGraph, mempty, mempty, blockCharInfo)

    -- Need to make node
    else

        -- check if children in graph
        let nodeChildren = LG.descendants simpleGraph curNode  -- should be 1 or 2, not zero since all leaves already in graph
            leftChild = head nodeChildren
            rightChild = last nodeChildren
            leftChildTree = postDecorateTree simpleGraph curDecGraph blockCharInfo leftChild
            rightLeftChildTree = if length nodeChildren == 2 then postDecorateTree simpleGraph (thd6 leftChildTree) blockCharInfo rightChild
                                 else leftChildTree
            newSubTree = thd6 rightLeftChildTree
            (leftChildLabel, rightChildLabel) = U.leftRightChildLabelBV (fromJust $ LG.lab newSubTree leftChild, fromJust $ LG.lab newSubTree rightChild)

        in

        if length nodeChildren > 2 then error ("Graph not dichotomous in postDecorateTree node " ++ show curNode ++ "\n" ++ LG.prettify simpleGraph)
        else if null nodeChildren then error ("Leaf not in graph in postDecorateTree node " ++ show curNode ++ "\n" ++ LG.prettify simpleGraph)

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
                newCost =  V.sum $ V.map V.sum $ V.map (V.map snd) newCharData
                newVertex = VertexInfo {  index = curNode
                                        , bvLabel = bvLabel leftChildLabel .|. bvLabel rightChildLabel
                                        , parents = V.fromList $ LG.parents simpleGraph curNode
                                        , children = V.fromList nodeChildren
                                        , nodeType = GO.getNodeType simpleGraph curNode
                                        , vertName = T.pack $ "HTU" ++ show curNode
                                        , vertData = V.map (V.map fst) newCharData
                                        , vertexResolutionData = mempty
                                        , vertexCost = newCost
                                        , subGraphCost = subGraphCost leftChildLabel + subGraphCost rightChildLabel + newCost
                                        }
                newEdgesLabel = EdgeInfo {    minLength = newCost / 2.0
                                            , maxLength = newCost / 2.0
                                            , midRangeLength = newCost / 2.0
                                            , edgeType = TreeEdge
                                         }
                newEdges = LG.toEdge <$> LG.out simpleGraph curNode
                newLEdges =  fmap (LG.toLEdge' newEdgesLabel) newEdges
                newGraph =  LG.insEdges newLEdges $ LG.insNode (curNode, newVertex) newSubTree

            in
            -- trace ("New vertex:" ++ (show newVertex) ++ " at cost " ++ (show newCost)) (
            -- Do we need to PO.divideDecoratedGraphByBlockAndCharacterTree if not root?  probbaly not
            if nodeType newVertex == RootNode then (simpleGraph, subGraphCost newVertex, newGraph, mempty, PO.divideDecoratedGraphByBlockAndCharacterTree newGraph, blockCharInfo)
            else (simpleGraph, subGraphCost newVertex, newGraph, mempty, mempty, blockCharInfo)

            -- ) -- )

{-  This refactored out above
multiTraverseFullyLabelSoftWired' :: GlobalSettings -> ProcessedData -> Bool -> Bool -> DecoratedGraph -> Maybe Int -> SimpleGraph -> PhylogeneticGraph
multiTraverseFullyLabelSoftWired' inGS inData pruneEdges warnPruneEdges leafGraph startVertex inSimpleGraph =
    if LG.isEmpty inSimpleGraph then emptyPhylogeneticGraph
    else
        let -- for special casing of nonexact and single exact characters
            nonExactChars = U.getNumberNonExactCharacters (thd3 inData)
            -- exactCharacters = U.getNumberExactCharacters (thd3 inData)


            -- post order pass-- rerooted because to keep track of rerooting in Graphs
            outgroupRootedSoftWiredPostOrder = if startVertex == Nothing then postOrderSoftWiredTraversal inGS inData leafGraph startVertex $ GO.rerootTree' inSimpleGraph (outgroupIndex inGS)
                                               else postOrderSoftWiredTraversal inGS inData leafGraph startVertex inSimpleGraph

            -- root node--invariant even when rerooted
            -- (outgroupRootIndex, _) = head $ LG.getRoots $ thd6 outgroupRootedSoftWiredPostOrder

            startVertexList = if startVertex == Nothing then fmap fst $ LG.getRoots $ thd6 outgroupRootedSoftWiredPostOrder
                              else [fromJust startVertex]

            childrenOfRoot = concatMap (LG.descendants (thd6 outgroupRootedSoftWiredPostOrder)) startVertexList
            grandChildrenOfRoot = concatMap (LG.descendants (thd6 outgroupRootedSoftWiredPostOrder)) childrenOfRoot

            -- create list of multi-traversals with original rooting first
            -- subsequent rerooting do not reoptimize exact characters (add nonadd etc) 
            -- they are taken from the fully labelled first input decorated graph later when output graph created
            -- it is important that the first graph be the ourgroup rooted graph (outgroupRootedPhyloGraph) so this 
            -- will have the preoder assignmentsd for th eoutgroup rooted graph as 3rd field.  This can be used for incremental
            -- optimization to get O(log n) initial postorder assingment when mutsating graph.
            recursiveRerootList = outgroupRootedSoftWiredPostOrder : minimalReRootPhyloGraph inGS SoftWired outgroupRootedSoftWiredPostOrder (head startVertexList) grandChildrenOfRoot

            finalizedPostOrderGraphList = L.sortOn snd6 $ fmap (updateAndFinalizePostOrderSoftWired startVertex (head startVertexList)) recursiveRerootList

            -- create optimal final graph with best costs and best traversal (rerooting) forest for each character
            -- traversal for exact characters (and costs) are the first of each least since exact only optimizaed for that 
            -- traversal graph.  The result has approprotate post-order assignments for traversals, preorder "final" assignments
            -- are propagated to the Decorated graph field after the preorder pass.
            -- doesn't have to be sorted, but should minimize assignments
            graphWithBestAssignments = L.foldl1' setBetterGraphAssignment finalizedPostOrderGraphList

            -- same root cost if same data and number of roots
            localRootCost = if (rootCost inGS) == NoRootCost then 0.0
                       else if (rootCost inGS) == Wheeler2015Root then getW15RootCost inData outgroupRootedSoftWiredPostOrder
                       else error ("Root cost type " ++ (show $ rootCost inGS) ++ " is not yet implemented")

        in

        -- Update post-order:
        --  1) traceback to get best resolutions
        --  2) back proppagate resolutions to vertex info
        --  3) update display/character trees
        -- Preorder on updated post-order graph 
        -- trace ("Pre-penalty: " ++ show (fmap snd6 finalizedPostOrderGraphList)) (

        -- only static characters
        if nonExactChars == 0 then 
            let penaltyFactor  = if (graphFactor inGS) == NoNetworkPenalty then 0.0
                                 else if (graphFactor inGS) == Wheeler2015Network then getW15NetPenalty outgroupRootedSoftWiredPostOrder
                                 else error ("Network penalty type " ++ (show $ graphFactor inGS) ++ " is not yet implemented")

                outgroupRootedSoftWiredPostOrder' = updatePhylogeneticGraphCost outgroupRootedSoftWiredPostOrder (penaltyFactor + localRootCost + (snd6 outgroupRootedSoftWiredPostOrder))

                fullyOptimizedGraph = PRE.preOrderTreeTraversal inGS (finalAssignment inGS) False (head startVertexList) outgroupRootedSoftWiredPostOrder'
            in
            checkUnusedEdgesPruneInfty inGS inData pruneEdges warnPruneEdges leafGraph fullyOptimizedGraph

        -- single dynamic character--checks rerooted values
        else if nonExactChars == 1 then
            let penaltyFactorList  = if (graphFactor inGS) == NoNetworkPenalty then replicate (length finalizedPostOrderGraphList) 0.0
                                     else if (graphFactor inGS) == Wheeler2015Network then fmap getW15NetPenalty finalizedPostOrderGraphList
                                     else error ("Network penalty type " ++ (show $ graphFactor inGS) ++ " is not yet implemented")
                newCostList = zipWith3 U.add3 penaltyFactorList (replicate (length finalizedPostOrderGraphList) localRootCost) (fmap snd6 finalizedPostOrderGraphList)

                finalizedPostOrderGraph = head $ L.sortOn snd6 $ zipWith updatePhylogeneticGraphCost finalizedPostOrderGraphList newCostList

                fullyOptimizedGraph = PRE.preOrderTreeTraversal inGS (finalAssignment inGS) True (head startVertexList) finalizedPostOrderGraph
            in
            checkUnusedEdgesPruneInfty inGS inData pruneEdges warnPruneEdges leafGraph fullyOptimizedGraph

        -- multiple dynamic characters--cjecks for best root for each character
        else 
            let penaltyFactor  = if (graphFactor inGS) == NoNetworkPenalty then 0.0
                                 else if (graphFactor inGS) == Wheeler2015Network then getW15NetPenalty graphWithBestAssignments
                                 else error ("Network penalty type " ++ (show $ graphFactor inGS) ++ " is not yet implemented")

                graphWithBestAssignments' = updatePhylogeneticGraphCost graphWithBestAssignments (penaltyFactor +  localRootCost + (snd6 graphWithBestAssignments))

                fullyOptimizedGraph = PRE.preOrderTreeTraversal inGS (finalAssignment inGS) True (head startVertexList) graphWithBestAssignments'
            in
            checkUnusedEdgesPruneInfty inGS inData pruneEdges warnPruneEdges leafGraph fullyOptimizedGraph
        -- )

-}
{-Refactored above
multiTraverseFullyLabelTree' :: GlobalSettings -> ProcessedData -> DecoratedGraph -> Maybe Int -> SimpleGraph -> PhylogeneticGraph
multiTraverseFullyLabelTree' inGS inData leafGraph startVertex inSimpleGraph =
    if LG.isEmpty inSimpleGraph then emptyPhylogeneticGraph
    else
        let nonExactChars = U.getNumberNonExactCharacters (thd3 inData)
            -- exactCharacters = U.getNumberExactCharacters (thd3 inData)

            -- initial traversal based on global outgroup and the "next" traversal points as children of existing traversal
            -- here initial root.
            outgroupRootedPhyloGraph = if startVertex == Nothing then postOrderTreeTraversal inGS inData leafGraph startVertex $ GO.rerootTree' inSimpleGraph (outgroupIndex inGS)
                                       else postOrderTreeTraversal inGS inData leafGraph startVertex inSimpleGraph

            startVertexList = if startVertex == Nothing then fmap fst $ LG.getRoots $ thd6 outgroupRootedPhyloGraph
                              else [fromJust startVertex]

                                      
            childrenOfRoot = concatMap (LG.descendants (thd6 outgroupRootedPhyloGraph)) startVertexList
            grandChildrenOfRoot = concatMap (LG.descendants (thd6 outgroupRootedPhyloGraph)) childrenOfRoot

            -- create list of multi-traversals with original rooting first
            -- subsequent rerooting do not reoptimize exact characters (add nonadd etc) 
            -- they are taken from the fully labelled first input decorated graph later when output graph created
            -- it is important that the first graph be the ourgroup rooted graph (outgroupRootedPhyloGraph) so this 
            -- will have the preoder assignmentsd for th eoutgroup rooted graph as 3rd field.  This can be used for incremental
            -- optimization to get O(log n) initial postorder assingment when mutsating graph.
            recursiveRerootList = outgroupRootedPhyloGraph : minimalReRootPhyloGraph inGS (graphType inGS) outgroupRootedPhyloGraph (head startVertexList) grandChildrenOfRoot

            finalizedPostOrderGraphList = L.sortOn snd6 recursiveRerootList

            -- create optimal final graph with best costs and best traversal (rerooting) forest for each character
            -- traversal for exact characters (and costs) are the first of each least since exact only optimizaed for that 
            -- traversal graph.  The result has approprotate post-order assignments for traversals, preorder "final" assignments
            -- are propagated to the Decorated graph field after the preorder pass.
            -- sorting may help with assignments (fewer)
            graphWithBestAssignments' = L.foldl1' setBetterGraphAssignment finalizedPostOrderGraphList -- (recursiveRerootList !! 0) (recursiveRerootList !! 1) 

            -- this for debuggin purposes
            --allPreorderList = fmap (preOrderTreeTraversal inGS (finalAssignment inGS)) recursiveRerootList

        in
        -- Uncomment this to (and comment the following three cases) avoid traversal rerooting stuff for debugging
        --preOrderTreeTraversal inGS (finalAssignment inGS) outgroupRootedPhyloGraph
        --preOrderTreeTraversal inGS  (finalAssignment inGS) $ head minCostGraphListRecursive 

        -- special cases that don't require all the work

        -- trace ("Nums:" ++ show (length minCostGraphListRecursive) ++ " " ++ show (fmap snd6 minCostGraphListRecursive)) (
        if nonExactChars == 0 then PRE.preOrderTreeTraversal inGS (finalAssignment inGS) False (head startVertexList) outgroupRootedPhyloGraph
        else if nonExactChars == 1 then PRE.preOrderTreeTraversal inGS (finalAssignment inGS) True  (head startVertexList) $ head finalizedPostOrderGraphList
        else PRE.preOrderTreeTraversal inGS (finalAssignment inGS) True  (head startVertexList) graphWithBestAssignments'
        -- )

-}

