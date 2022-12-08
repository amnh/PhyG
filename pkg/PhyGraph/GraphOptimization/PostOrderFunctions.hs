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

module GraphOptimization.PostOrderFunctions  ( rerootPhylogeneticGraph
                                             , rerootPhylogeneticGraph'
                                             , rerootPhylogeneticNetwork
                                             , rerootPhylogeneticNetwork'
                                             , updateDisplayTreesAndCost
                                             ) where

import           Data.Bits
import qualified Data.List                   as L
import           Data.Maybe
import qualified Data.Text.Lazy              as T
import qualified Data.Vector                 as V
import           GeneralUtilities
import qualified Graphs.GraphOperations      as GO
import           Types.Types
import qualified Utilities.LocalGraph        as LG
import qualified Utilities.Utilities         as U
import qualified GraphOptimization.PostOrderSoftWiredFunctions as POSW
import qualified GraphOptimization.PostOrderSoftWiredFunctionsNew as NEW
import           Debug.Trace


-- | updateDisplayTreesAndCost takes a softwired graph and updates
-- display trees and graph cost based on resolutions at root
updateDisplayTreesAndCost :: PhylogeneticGraph -> PhylogeneticGraph
updateDisplayTreesAndCost inGraph =
    if LG.isEmpty $ thd6 inGraph then emptyPhylogeneticGraph
    else
        -- True for check popCount at root fort valid resolution (all leaves in graph)
        let (_, outgroupRootLabel) =  head $ LG.getRoots (thd6 inGraph)
            (displayGraphVL, lDisplayCost) = NEW.extractDisplayTrees Nothing True (vertexResolutionData outgroupRootLabel)
        in
        --trace ("UDTC: " ++ show lDisplayCost)
        (fst6 inGraph, lDisplayCost, thd6 inGraph, displayGraphVL, fft6 inGraph, six6 inGraph)

-- | reOptimizeNodes takes a decorated graph and a list of nodes and reoptimizes (relabels)
-- them based on children in input graph
-- simple recursive since each node depends on children
-- remove check for debugging after it works
-- check for out-degree 1, doesn't matter for trees however.
reOptimizeNodes :: GlobalSettings -> V.Vector (V.Vector CharInfo) -> DecoratedGraph -> [LG.LNode VertexInfo] -> DecoratedGraph
reOptimizeNodes inGS charInfoVectVect inGraph oldNodeList =
  -- trace ("RON:" ++ (show $ fmap fst oldNodeList)) (
  if null oldNodeList then inGraph
  else
    -- make sure that nodes are optimized in correct order so that nodes are only reoptimized using updated children
    -- this should really not have to happen--order should be determined a priori
    let curNode@(curNodeIndex, curNodeLabel) = head oldNodeList
        nodeChildren = LG.descendants inGraph curNodeIndex  -- should be 1 or 2, not zero since all leaves already in graph
        foundCurChildern = filter (`elem` nodeChildren) $ fmap fst (tail oldNodeList)
    in
    if LG.isLeaf inGraph curNodeIndex then trace ("Should not be a leaf in reoptimize nodes: " ++ (show curNodeIndex) ++ " children " ++ (show nodeChildren) ++ "\nGraph:\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inGraph)) inGraph

    -- if node in multiple times due to network--put off optimizatin till last time
    else if curNodeIndex `elem` (fmap fst $ tail oldNodeList) then reOptimizeNodes inGS charInfoVectVect inGraph (tail oldNodeList)

    else if not $ null foundCurChildern then
      -- trace ("Current node " ++ (show curNodeIndex) ++ " has children " ++ (show nodeChildren) ++ " in optimize list (optimization order error)" ++ (show $ fmap fst $ tail oldNodeList))
      reOptimizeNodes inGS charInfoVectVect inGraph (tail oldNodeList ++ [curNode])


    -- somehow root before others -- remove if not needed after debug
    else if LG.isRoot inGraph curNodeIndex && length oldNodeList > 1 then
        error ("Root first :" ++ (show $ fmap fst oldNodeList) ++ "RC " ++ show (LG.descendants inGraph curNodeIndex)) -- ++ "\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inGraph))
        --reOptimizeNodes inGS localGraphType charInfoVectVect inGraph ((tail oldNodeList) ++ [curNode])

    else
        -- trace ("RON: " ++ (show curNodeIndex) ++ " children " ++ (show nodeChildren)) (
        let leftChild = head nodeChildren
            rightChild = last nodeChildren
            -- leftChildLabel = fromJust $ LG.lab inGraph leftChild
            -- rightChildLabel = fromJust $ LG.lab inGraph rightChild

            -- this ensures that left/right choices are based on leaf BV for consistency and label invariance
            (leftChildLabel, rightChildLabel) = U.leftRightChildLabelBV (fromJust $ LG.lab inGraph leftChild, fromJust $ LG.lab inGraph rightChild)
            newVertexData = POSW.createVertexDataOverBlocksNonExact (vertData leftChildLabel) (vertData  rightChildLabel) charInfoVectVect []
            -- newVertexData = createVertexDataOverBlocks  (vertData leftChildLabel) (vertData  rightChildLabel) charInfoVectVect []
        in
        {-
        --debug remove when not needed--checking to see if node should not be re optimized
        if (sort nodeChildren) == (sort $ V.toList $ children curnodeLabel) then
          trace ("Children for vertex unchanged " ++ (show curNodeIndex)
          reOptimizeNodes localGraphType charInfoVectVect inGraph (tail oldNodeList)
        else
        -}
        if (graphType inGS) == Tree || (graphType inGS) == HardWired then
           let newCost =  if length nodeChildren < 2 then 0
                          else V.sum $ V.map V.sum $ V.map (V.map snd) newVertexData
               newVertexLabel = VertexInfo {  index = curNodeIndex
                                            -- this bit labelling incorect for outdegree = 1, need to prepend bits
                                            , bvLabel = if length nodeChildren < 2 then bvLabel leftChildLabel
                                                        else bvLabel leftChildLabel .|. bvLabel rightChildLabel
                                            , parents = V.fromList $ LG.parents inGraph curNodeIndex
                                            , children = V.fromList nodeChildren
                                            , nodeType = GO.getNodeType inGraph curNodeIndex -- nodeType curNodeLabel
                                            , vertName = vertName curNodeLabel
                                            , vertexResolutionData = mempty
                                            , vertData = if length nodeChildren < 2 then vertData leftChildLabel
                                                         else V.map (V.map fst) newVertexData
                                            , vertexCost = newCost
                                            , subGraphCost = if length nodeChildren < 2 then subGraphCost leftChildLabel
                                                             else subGraphCost leftChildLabel + subGraphCost rightChildLabel + newCost
                                            }
               -- this to add back edges deleted with nodes (undocumented but sensible in fgl)
               replacementEdges = LG.inn inGraph curNodeIndex ++ LG.out inGraph curNodeIndex
               newGraph = LG.insEdges replacementEdges $ LG.insNode (curNodeIndex, newVertexLabel) $ LG.delNode curNodeIndex inGraph
            in
            --trace ("New vertexCost " ++ show newCost) --  ++ " lcn " ++ (show (vertData leftChildLabel, vertData rightChildLabel, vertData curnodeLabel)))
            reOptimizeNodes inGS charInfoVectVect newGraph (tail oldNodeList)

        else if  (graphType inGS) == SoftWired then
            -- trace ("Reoptimizing " ++ (show curNodeIndex)) (
            -- single child of node (can certinly happen with soft-wired networks
            if length nodeChildren == 1 then
                --trace ("Out=1\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inGraph)) (
                let (_,_, newVertexLabel, _, _) = NEW.getOutDegree1VertexAndGraph curNodeIndex leftChildLabel inGraph nodeChildren inGraph

                    -- this to add back edges deleted with nodes (undocumented but sensible in fgl)
                    replacementEdges = LG.inn inGraph curNodeIndex ++ LG.out inGraph curNodeIndex
                    newGraph = LG.insEdges replacementEdges $ LG.insNode (curNodeIndex, newVertexLabel) $ LG.delNode curNodeIndex inGraph
                in
                reOptimizeNodes inGS charInfoVectVect newGraph (tail oldNodeList)
                -- )

            -- two children
            else
                -- trace ("Out=2 " ++ (show $ length nodeChildren)) (
                let -- this ensures that left/right choices are based on leaf BV for consistency and label invariance
                    -- larger bitvector is Right, smaller or equal Left
                    ((leftChild', leftChildLabel'), (rightChild', rightChildLabel')) = U.leftRightChildLabelBVNode ((leftChild, fromJust $ LG.lab inGraph leftChild), (rightChild, fromJust $ LG.lab inGraph rightChild))

                    -- create resolution caches for blocks
                    leftChildNodeType  = GO.getNodeType inGraph leftChild' -- nodeType leftChildLabel'
                    rightChildNodeType = GO.getNodeType inGraph rightChild' -- nodeType rightChildLabel'
                    resolutionBlockVL = V.zipWith3 (NEW.createBlockResolutions (compressResolutions inGS) curNodeIndex leftChild' rightChild' leftChildNodeType rightChildNodeType (nodeType curNodeLabel)) (vertexResolutionData leftChildLabel') (vertexResolutionData rightChildLabel') charInfoVectVect

                    -- create canonical Decorated Graph vertex
                    -- 0 cost becasue can't know cosrt until hit root and get best valid resolutions
                    newVertexLabel = VertexInfo {  index = curNodeIndex
                                                , bvLabel = bvLabel leftChildLabel' .|. bvLabel rightChildLabel'
                                                , parents = V.fromList $ LG.parents inGraph curNodeIndex
                                                , children = V.fromList nodeChildren
                                                , nodeType = GO.getNodeType inGraph curNodeIndex -- nodeType curNodeLabel
                                                , vertName = T.pack $ "HTU" ++ show curNodeIndex
                                                , vertData = mempty --empty because of resolution data
                                                , vertexResolutionData = resolutionBlockVL
                                                , vertexCost = 0.0 -- newCost
                                                , subGraphCost = 0.0 -- (subGraphCost leftChildLabel) + (subGraphCost rightChildLabel) + newCost
                                                }

                    {-
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
                    -}

                    replacementEdges = LG.inn inGraph curNodeIndex ++ LG.out inGraph curNodeIndex
                    newGraph =  LG.insEdges replacementEdges $ LG.insNode (curNodeIndex, newVertexLabel) $ LG.delNode curNodeIndex inGraph
                in
                -- trace ("Resolution Data: \n" ++ "left child\n" ++ (show $ vertexResolutionData leftChildLabel) ++ "\nright child\n" ++ (show $ vertexResolutionData rightChildLabel)
                --    ++ "\nCur Node\n" ++ (show $ vertexResolutionData newVertexLabel))
                reOptimizeNodes inGS charInfoVectVect newGraph (tail oldNodeList)

                -- ) -- )

        else  errorWithoutStackTrace ("Graph type unrecognized/not yet implemented: " ++ show (graphType inGS))
        -- )


-- | rerootPhylogeneticNetwork take a vertex index and reroots phylogenetic network
-- roots as for tree, but when a netwrok edge is chosen as root edge, its 'sister' network edge
-- is deleted (leaving two in=out=1 nodes), th ereooting takes place,
-- and the edge is re-added in both directions.  Nodes may change between
-- tree and network.  This is updated
rerootPhylogeneticNetwork :: GlobalSettings -> Int -> Int -> PhylogeneticGraph -> PhylogeneticGraph
rerootPhylogeneticNetwork inGS originalRootIndex rerootIndex inGraph@(inSimple, _, inDecGraph, _, _, _) =
    if LG.isEmpty inSimple then inGraph
    else
        let newRootParents = LG.parents inSimple rerootIndex
            newRootLabel = LG.lab inDecGraph rerootIndex
            parentNodeLabel = LG.lab inDecGraph $ head newRootParents
        in

        if isNothing newRootLabel then error ("New root has no label: " ++ show rerootIndex)

        else if isNothing parentNodeLabel then error ("Parent of new root has no label: " ++ show rerootIndex)

        -- OK to reroot
        else
            --if length newRootParents > 1 then trace ("Root on network edge") inGraph
            --else
                let isNetworkNode = nodeType (fromJust newRootLabel) == NetworkNode
                    parentIsNetworkNode = nodeType (fromJust parentNodeLabel) == NetworkNode
                in
                -- if (not isNetworkNode) && (not parentIsNetworkNode) then
                rerootPhylogeneticGraph inGS isNetworkNode originalRootIndex parentIsNetworkNode rerootIndex inGraph
                -- else inGraph



-- | rerootPhylogeneticNetwork' flipped version of rerootPhylogeneticNetwork
rerootPhylogeneticNetwork' :: GlobalSettings -> PhylogeneticGraph -> Int -> Int -> PhylogeneticGraph
rerootPhylogeneticNetwork' inGS inGraph originalRootIndex rerootIndex = rerootPhylogeneticNetwork inGS originalRootIndex rerootIndex inGraph

-- | rerootPhylogeneticGraph' flipped version of rerootPhylogeneticGraph
rerootPhylogeneticGraph' :: GlobalSettings -> Bool -> Bool -> PhylogeneticGraph -> Int -> Int -> PhylogeneticGraph
rerootPhylogeneticGraph' inGS isNetworkNode parentIsNetworkNode inGraph originalRootIndex rerootIndex = rerootPhylogeneticGraph inGS isNetworkNode originalRootIndex parentIsNetworkNode rerootIndex inGraph

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
rerootPhylogeneticGraph ::  GlobalSettings -> Bool -> Int ->  Bool -> Int -> PhylogeneticGraph -> PhylogeneticGraph
rerootPhylogeneticGraph  inGS isNetworkNode originalRootIndex parentIsNetworkNode rerootIndex inPhyGraph@(inSimple, _, inDecGraph, blockGraphVV, _, charInfoVectVect) =
  if LG.isEmpty inSimple then inPhyGraph
  -- else if inCost == 0 then error ("Input graph with cost zero--likely non decorated input graph in rerootPhylogeneticGraph\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inDecGraph))
  else
    let -- decorated graph Boolean to specify that non-exact characters need to be reoptimized if affected
        -- could just update with needges? from simple graph rerooting
        (newDecGraph, touchedNodes)  = if (graphType inGS) == Tree then (LG.rerootTree rerootIndex inDecGraph, [])
                                       else if (graphType inGS) == SoftWired then rectifyGraphDecorated isNetworkNode originalRootIndex parentIsNetworkNode rerootIndex inDecGraph
                                       else if (graphType inGS) == HardWired then
                                         if (not . LG.cyclic) inSimple then
                                            let (newDecGraph'', touchedNodes') = rectifyGraphDecorated isNetworkNode originalRootIndex parentIsNetworkNode rerootIndex inDecGraph
                                            in
                                            if (not . LG.cyclic) newDecGraph'' then (newDecGraph'', touchedNodes')
                                            else (LG.empty, [])
                                         else (LG.empty, [])
                                       else error ("Error--Graph type unimplemented: " ++ (show (graphType inGS)))

        newSimpleGraph = GO.convertDecoratedToSimpleGraph newDecGraph


        -- reoptimize nodes here
        -- nodes on spine from new root to old root that needs to be reoptimized
        -- THIS IS WRONG FOR SOFTWIRED--extra nodes for  rectifying graph
        (nodesToOptimize, _) = if (graphType inGS) == Tree then LG.pathToRoot inDecGraph (rerootIndex, fromJust $ LG.lab inDecGraph rerootIndex)
                               else if (graphType inGS) == SoftWired then (touchedNodes, [])
                               else if (graphType inGS) == HardWired then (touchedNodes, [])
                               else error ("Error--Graph type unimplemented: " ++ (show (graphType inGS)))

        -- this only reoptimizes non-exact characters since rerooting doesn't affect 'exact" character optimization'
        newDecGraph' = reOptimizeNodes inGS charInfoVectVect newDecGraph nodesToOptimize -- (L.nub $ nodesToOptimize ++ touchedNodes)
        -- newDecGraph' = reOptimizeNodes inGS charInfoVectVect newDecGraph nodesToOptimize (L.nub $ nodesToOptimize ++ touchedNodes)

        -- sum of root costs on Decorated graph
        newGraphCost = sum $ fmap subGraphCost $ fmap snd $ LG.getRoots newDecGraph'

        -- rerooted display forests--don't care about costs--I hope (hence Bool False)
        newblockGraphVV = if V.null blockGraphVV then mempty
                                  --else fmap (fmap (GO.rerootTree rerootIndex)) blockGraphVV
                                  else fmap (fmap (rectifyGraphDecorated' isNetworkNode originalRootIndex parentIsNetworkNode rerootIndex)) blockGraphVV

        in
        --trace ("rerootPhylogeneticGraph:\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inDecGraph) ++ "\nNew\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph newDecGraph) (

        -- this for forbiden condition where rerooting a graph creates parent and child nodes boyth network
        if (null touchedNodes) && ((graphType inGS) /= Tree) then emptyPhylogeneticGraph
        else if newSimpleGraph == LG.empty then emptyPhylogeneticGraph

        -- Same root, so no need to redo
        -- else if (length nodesToOptimize == 1) then inPhyGraph
        else

          {-
          trace ("New RootIndex :" ++ (show (rerootIndex, isNetworkNode, parentIsNetworkNode)) ++ " To optimize:" ++ (show $ fmap fst nodesToOptimize) ++ "\nOG:\n"
          ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inDecGraph) ++ "\nRRG:" ++ ((LG.prettify $ GO.convertDecoratedToSimpleGraph newDecGraph))) (
          -}
          -- (newSimpleGraph, newGraphCost, newDecGraph', newDecoratedGraphVect, V.replicate (length charInfoVectVect) (V.singleton newDecGraph'), charInfoVectVect)
          if ((graphType inGS) == Tree || (graphType inGS) == HardWired) then (newSimpleGraph, newGraphCost, newDecGraph', newblockGraphVV, (snd $ POSW.divideDecoratedGraphByBlockAndCharacterTree newDecGraph'), charInfoVectVect)
          else
            -- get root resolutions and cost
            let (displayGraphVL, lDisplayCost) = NEW.extractDisplayTrees (Just originalRootIndex) True (vertexResolutionData $ fromJust $ LG.lab newDecGraph' originalRootIndex)
            in
            (newSimpleGraph, lDisplayCost, newDecGraph', displayGraphVL, mempty, charInfoVectVect)
            -- )


-- | rectifyGraphDecorated' wrapper around rectifyGraphDecorated for graph only
rectifyGraphDecorated' :: Bool -> Int ->  Bool -> Int -> DecoratedGraph -> DecoratedGraph
rectifyGraphDecorated' isNetworkNode originalRootIndex parentIsNetworkNode rerootIndex inGraph = fst $ rectifyGraphDecorated isNetworkNode originalRootIndex parentIsNetworkNode rerootIndex inGraph


-- | rectifyGraphDecorated 'fixes' (flips) edges where a network edge has be chosen as a reroot edge For Decorated Graph--thee should be able to be combined
-- this can be abstracted if dummy graph set to some input edge label
rectifyGraphDecorated :: Bool -> Int ->  Bool -> Int -> DecoratedGraph -> (DecoratedGraph, [LG.LNode VertexInfo])
rectifyGraphDecorated isNetworkNode originalRootIndex parentIsNetworkNode rerootIndex inGraph =
    if LG.isEmpty inGraph then (LG.empty, [])
    else
        -- sanity check of phylogenetic graph
        if isNetworkNode || parentIsNetworkNode then
            -- trace ("Graph with parent and child nodes network vertices--skipping reroot")
            (LG.empty, [])
        -- can't reroot on network node--can cause alot of problems
        else
            -- get nodes and edged on path from old to new root
            let (nodePathToRoot, edgePathToRoot') = LG.pathToRoot inGraph (rerootIndex, fromJust $ LG.lab inGraph rerootIndex)
                edgePathToRoot =  fmap LG.toEdge edgePathToRoot'
                origRootEdges = fmap LG.toEdge $ LG.out inGraph originalRootIndex
                origVirtualRootEdge = if ((snd $ head origRootEdges) `elem` (fmap fst nodePathToRoot)) then GO.makeDummyLabEdge dummyEdge (snd $ head origRootEdges, snd $ last origRootEdges)
                                      else GO.makeDummyLabEdge dummyEdge (snd $ last origRootEdges, snd $ head origRootEdges)

                -- arbitrarily chooses one of multiple parent vertices is network edge
                -- dummy third field for new root edges
                parentNewRoot = LG.parents inGraph rerootIndex
                newRootEdge =  (head parentNewRoot ,rerootIndex)
                otherEdgesFromParentNewRoot = fmap LG.toEdge $ filter ((/= rerootIndex) . snd3) $ LG.out inGraph (head parentNewRoot)
                newRootChildEdges = fmap (GO.makeDummyLabEdge dummyEdge) [(originalRootIndex, head parentNewRoot), (originalRootIndex, rerootIndex)]

                    --make new graph--init here to not flip original edge to root
                flippedEdgesToOldRoot = fmap (GO.makeDummyLabEdge dummyEdge) $ fmap LG.flipEdge $ edgePathToRoot L.\\ (newRootEdge : origRootEdges)

                --nodesTouchedFlippedEdgesIndices = (fmap fst3 flippedEdgesToOldRoot) `L.union` (fmap snd3 flippedEdgesToOldRoot)
                    --nodesTouchedFlippedLabels = fmap fromJust $ fmap (LG.lab inGraph) nodesTouchedFlippedEdgesIndices

                edgesToDelete = (newRootEdge :  origRootEdges) ++ edgePathToRoot

                -- check if network edghe  to be creted by rerooting is deleted
                edgesToAddBack = fmap (GO.makeDummyLabEdge dummyEdge) $ filter (`elem` edgesToDelete) otherEdgesFromParentNewRoot

                edgesToInsert = ((origVirtualRootEdge : newRootChildEdges) ++ flippedEdgesToOldRoot ++ edgesToAddBack) L.\\ (fmap LG.flipLEdge edgesToAddBack)

                newGraph = LG.insEdges edgesToInsert $ LG.delEdges edgesToDelete inGraph

                -- check for HTU with outdegree 0 due to rerooting issues--could have nested networks
                hasNetLeaf = True `elem` (fmap (LG.isNetworkLeaf newGraph) (LG.nodes newGraph))

                -- get touched nodes
                newRootNodeIndexList = [head parentNewRoot, rerootIndex]
                newRootNodeLabelIndex = fmap (fromJust . LG.lab inGraph) newRootNodeIndexList
                newRootNodeList = filter ((not . LG.isLeaf newGraph) . fst) $ zip newRootNodeIndexList newRootNodeLabelIndex

                nodesToReoptimize = nodePathToRoot `L.union` newRootNodeList


            in
            if length origRootEdges /= 2 then error ("Root does not have two children in rectifyGraphDecorated: " ++ (show origRootEdges))
            else if hasNetLeaf then
                --trace ("Graph with HTU network vertex--skipping reroot")
                (LG.empty, [])

            {-
            else if LG.cyclic newGraph then
                trace ("Cycle")
                (LG.empty, [])
            -}
            else
                {-
                trace ("Original root edges:" ++ (show origRootEdges)
                    ++ " Insertions:" ++ (show (LG.toEdge origVirtualRootEdge, fmap LG.toEdge newRootChildEdges, fmap LG.toEdge flippedEdgesToOldRoot, fmap LG.toEdge edgesToAddBack))
                    ++ "\nDeletions:" ++ (show ((newRootEdge, origRootEdges,edgePathToRoot))))
                -}
                {-
                if (length $ LG.getIsolatedNodes newGraph) > 0 || (length $ LG.getRoots newGraph) > 1 then
                    trace ("Isolated nodes: " ++ (show $ fmap fst $ LG.getIsolatedNodes newGraph) ++ " roots " ++ (show $ fmap fst $ LG.getRoots newGraph)) (newGraph, nodesToReoptimize)
                else
                -}
                (newGraph, nodesToReoptimize)

