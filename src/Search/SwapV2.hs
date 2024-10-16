{- |
Module specifying graph swapping rearrangement functions
-}
module Search.SwapV2 (
    reoptimizeSplitGraphFromVertexIANew,
    reoptimizeSplitGraphFromVertexNew,
    reoptimizeSplitGraphFromVertexTupleNew,
    swapV2,
) where

import Control.Monad (filterM)
import Control.Monad.IO.Class
import Control.Monad.Random.Class
import Data.Foldable (fold, toList)
import Data.Foldable1 (Foldable1)
import Data.Foldable1 qualified as F1
import Data.Functor (($>), (<&>))
import Data.List qualified as L
import Data.Maybe
import Data.Ord (comparing)
import Data.Vector qualified as V
import GHC.Real qualified as Real (infinity)
import GeneralUtilities
import GraphOptimization.Medians qualified as M
import GraphOptimization.PostOrderSoftWiredFunctions qualified as POSW
import GraphOptimization.PreOrderFunctions qualified as PRE
import GraphOptimization.Traversals qualified as T
import Graphs.GraphOperations qualified as GO
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import PHANE.Evaluation.Verbosity (Verbosity (..))
import Types.Types
import Utilities.LocalGraph qualified as LG
import Utilities.Utilities as U


{- | New Swap functions that are based on PHANE paralleization routines.
    1) Naive version
    2) Incremental for split and move evaluation SPR/TBR where possible
    3) steepest descent
        Order edges by distace to original break point on rejoin
        Continue when swtich to new (order of edges)
        Keep single graph only until no better, then multiple
    4) Randomized trajectories
    5) heuristic cost calculations
    6) unions
    7) SA/Drift    
-}
swapV2
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Int
    → [ReducedPhylogeneticGraph]
    → [(Maybe SAParams, ReducedPhylogeneticGraph)]
    → PhyG ([ReducedPhylogeneticGraph], Int)
swapV2 swapParams inGS inData inCounter curBestGraphList inSimAnnealParams =
    if null curBestGraphList
        then pure ([], inCounter)
        else swapNaive swapParams inGS inData inCounter curBestGraphList inSimAnnealParams


{- | Naive swap functions to create reference for later algorithmic improvements 
    1) Take first graph

    2) Create list of splits of graph 
        LG.getEdgeSplitList O(n)

    3) Rejoin all places for all splits
        This is "all around" in that not switching to lower cost graphs
        at first opportunity (ie. 'steepest").
        LG.splitGraphOnEdge (check base graph leaves >=3)
        Curretn function to split graph and reoptimize components: reoptimizeSplitGraphFromVertex
        LG.joinGraphOnEdge

    4) Full evaluation of graphs 
        Time complexities should be with O(n) for graph traversals
            NNI O(n^2)
            SPR O(n^3)
            TBR O(n^4)
-}
swapNaive 
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Int
    → [ReducedPhylogeneticGraph]
    → [(Maybe SAParams, ReducedPhylogeneticGraph)]
    → PhyG ([ReducedPhylogeneticGraph], Int)
swapNaive swapParams inGS inData inCounter curBestGraphList inSimAnnealParams =
    let inGraphPair = L.uncons curBestGraphList
    in
    if isNothing inGraphPair
        then pure ([], inCounter)
        else do
            -- if graph list not empty proceed 
            let (firstGraph, graphsRemaining) = fromJust inGraphPair

            -- get list of splittable edges (bridging, not adjascent to root)
            let edgeList' = L.uncons $ LG.getEdgeSplitList $ thd5 firstGraph
            if isNothing edgeList' then do
                                        logWith LogInfo $ "\tNo Breakable Edges"
                                        pure (graphsRemaining, inCounter + 1)
            else do
                let curBestCost = minimum $ fmap snd5 curBestGraphList
                let fullFirstGraph = GO.convertReduced2PhylogeneticGraph (head curBestGraphList)
                inGraphNetPenalty ← T.getPenaltyFactor inGS inData Nothing fullFirstGraph
                let inGraphNetPenaltyFactor = inGraphNetPenalty / curBestCost

                let edgeList@(firstEdge, restEdges) = fromJust edgeList'
                logWith LogInfo $ "\tBreakable Edges: " <> (show $ 1 + (length restEdges))

                -- split graph on the first edge
                let (splitGraph, graphRoot, prunedGraphRootIndex, originalConnectionOfPruned) = LG.splitGraphOnEdge (thd5 firstGraph) firstEdge

                -- split and optimize graph components (original for time complexity check)
                {-
                (reoptimizedSplitGraph, splitCost) ←
                        reoptimizeSplitGraphFromVertexOrig inGS inData (doIA swapParams) inGraphNetPenaltyFactor splitGraph graphRoot prunedGraphRootIndex
                logWith LogInfo $ "\tSplit Cost: " <> (show splitCost)
                -}
                (reoptimizedSplitGraph', splitCost') ←
                        reoptimizeSplitGraphFromVertexNew swapParams inGS inData (doIA swapParams) inGraphNetPenaltyFactor fullFirstGraph splitGraph graphRoot prunedGraphRootIndex 
                logWith LogInfo $ "\tSplit Cost New: " <> (show splitCost')
                pure (graphsRemaining, inCounter + 1)


{- | reoptimizeSplitGraphFromVertex 
    Original version of reoptimizeSplitGraphFromVertex from Swap
    copied here for testing convenience.
    roughly O(n)
-}
reoptimizeSplitGraphFromVertexNew
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Bool
    → VertexCost
    → PhylogeneticGraph
    → DecoratedGraph
    → LG.Node
    → LG.Node
    → PhyG (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertexNew swapParams inGS inData doIA netPenaltyFactor curGraph inSplitGraph startVertex prunedSubGraphRootVertex =
    -- trace ("RSGFV: " <> (show startVertex)) (
    if doIA
        then -- only reoptimize the IA states for dynamic characters
            reoptimizeSplitGraphFromVertexIANew swapParams inGS inData netPenaltyFactor curGraph inSplitGraph startVertex prunedSubGraphRootVertex 
        else -- perform full optimizations of nodes
            
            -- determine position to start incremental optimization for base graphs
            -- grandparent of pruned root due to recontruction of base graph--so use orig graph
            let parentPoint = L.uncons $ concat $ fmap (LG.parents (thd6 curGraph)) $ LG.parents (thd6 curGraph) prunedSubGraphRootVertex
                
                parentOriginalConnection = if isNothing parentPoint then error ("Split graph involving near-root edges: " <> (show prunedSubGraphRootVertex))
                                           else if (not . null . snd) $ fromJust parentPoint then error ("Split graph involving network edge: " <> (show (prunedSubGraphRootVertex, snd $ fromJust parentPoint)))
                                           else fst $ fromJust parentPoint

                -- these required for full optimization

                nonExactCharacters = U.getNumberSequenceCharacters (thd3 inData)
                origGraph = inSplitGraph -- thd5 origPhyloGraph
                leafGraph =
                    if graphType inGS == SoftWired
                        then LG.extractLeafGraph origGraph -- POSW.makeLeafGraphSoftWired inGS inData -- LG.extractLeafGraph origGraph
                        else LG.extractLeafGraph origGraph
                calcBranchLengths = False

                -- this for multitravers in swap for softwired to turn off
                multiTraverse =
                    if graphType inGS /= HardWired
                        then multiTraverseCharacters inGS
                        else False

                -- create simple graph version of split for post order pass
                splitGraphSimple = GO.convertDecoratedToSimpleGraph inSplitGraph
            in  do
                    -- create optimized base graph
                    -- False for staticIA
                    (postOrderBaseGraph, _) ←
                        T.generalizedGraphPostOrderTraversal
                            (inGS{graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse})
                            nonExactCharacters
                            inData
                            (Just (thd6 curGraph, parentOriginalConnection))
                            leafGraph
                            False
                            (Just startVertex)
                            splitGraphSimple

                    fullBaseGraph ←
                        PRE.preOrderTreeTraversal
                            (inGS{graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse})
                            (finalAssignment inGS)
                            False
                            calcBranchLengths
                            (nonExactCharacters > 0)
                            startVertex
                            True
                            postOrderBaseGraph

                    -- create fully optimized pruned graph.  Post order then preorder

                    -- get root node of pruned graph--parent since that is the full pruned piece (keeping that node for addition to base graph and edge creation)
                    let startPrunedNode = (prunedSubGraphRootVertex, fromJust $ LG.lab origGraph prunedSubGraphRootVertex)
                    let startPrunedParentNode = head $ LG.labParents origGraph prunedSubGraphRootVertex
                    let startPrunedParentEdge = (fst startPrunedParentNode, prunedSubGraphRootVertex, dummyEdge)

                    -- False for staticIA
                    -- Can use existing postOrder assignments for pruned
                    -- updating pruned root final to preliminary for preorder pass
                    let postOrderPrunedGraph = if LG.isLeaf origGraph prunedSubGraphRootVertex then curGraph
                                               else 
                                                    -- updating pruned root final to preliminary for preorder pass
                                                    let origLabelPrunedRoot = fromJust $ (LG.lab (thd6 curGraph) prunedSubGraphRootVertex)
                                                        originalVertInfoPrunedRoot = vertData origLabelPrunedRoot
                                                        prunedGraphCost = subGraphCost origLabelPrunedRoot
                                                        newVertInfoPrunedRoot = PRE.setFinalToPreliminaryStates originalVertInfoPrunedRoot
                                                        newFinalPrunedRoot = (prunedSubGraphRootVertex, origLabelPrunedRoot {vertData = newVertInfoPrunedRoot})
                                                        newPrunedRootOutEdges = LG.out (thd6 curGraph) prunedSubGraphRootVertex
                                                        newPrunedGraph =  LG.insEdges newPrunedRootOutEdges $ LG.insNode newFinalPrunedRoot $ LG.delNode prunedSubGraphRootVertex (thd6 curGraph)
                                                    in (fst6 curGraph, prunedGraphCost, newPrunedGraph, fth6 curGraph, fft6 curGraph, six6 curGraph)
                    
                    -- if not TBR then don't need preorder assignments
                    fullPrunedGraph ← if swapType swapParams `elem` [NNI, SPR] then 
                                            pure postOrderPrunedGraph
                                      else          
                                            PRE.preOrderTreeTraversal
                                                (inGS{graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse})
                                                (finalAssignment inGS)
                                                False
                                                calcBranchLengths
                                                (nonExactCharacters > 0)
                                                prunedSubGraphRootVertex
                                                True
                                                postOrderPrunedGraph
                                        

                    --logWith LogInfo $ "\nPruned Cost: " <> (show $ subGraphCost $ fromJust $ (LG.lab (thd6 curGraph) prunedSubGraphRootVertex))
                    --logWith LogInfo $ "\nPruned Cost: " <> (show $ subGraphCost $ fromJust $ (LG.lab (thd6 fullPrunedGraph) prunedSubGraphRootVertex))
                    

                    -- get root node of base graph
                    let startBaseNode = (startVertex, fromJust $ LG.lab (thd6 fullBaseGraph) startVertex)

                    -- get nodes and edges in base and pruned graph (both PhylogeneticGrapgs so thd5)
                    let (baseGraphNonRootNodes, baseGraphEdges) = LG.nodesAndEdgesAfter (thd6 fullBaseGraph) [startBaseNode]

                    let (prunedGraphNonRootNodes, prunedGraphEdges) =
                            if LG.isLeaf origGraph prunedSubGraphRootVertex
                                then ([], [])
                                -- else LG.nodesAndEdgesAfter (thd6 fullPrunedGraph) [startPrunedNode]
                                else LG.nodesAndEdgesAfter (thd6 fullPrunedGraph) [startPrunedNode]

                    -- make fully optimized graph from base and split components
                    let fullSplitGraph =
                            LG.mkGraph
                                ([startBaseNode, startPrunedNode, startPrunedParentNode] <> baseGraphNonRootNodes <> prunedGraphNonRootNodes)
                                (startPrunedParentEdge : (baseGraphEdges <> prunedGraphEdges))

                    -- cost of split graph to be later combined with re-addition delta for heuristic graph cost
                    let prunedCost = 
                            if LG.isLeaf origGraph prunedSubGraphRootVertex
                                then 0
                                else snd6 postOrderPrunedGraph
                    let splitGraphCost = ((1.0 + netPenaltyFactor) * ((snd6 fullBaseGraph) + prunedCost))

                    {-
                    -- check fo unlabbeld nodes
                    coninicalNodes =  LG.labNodes fullSplitGraph
                    nodeLabels = fmap (LG.lab fullSplitGraph) (fmap fst coninicalNodes)
                    unlabelledNodes = filter ((== Nothing) .snd) $ (zip (fmap fst coninicalNodes) nodeLabels)
                    -}

                    if prunedCost == infinity || (snd6 fullBaseGraph) == infinity
                        then pure (LG.empty, infinity)
                        else pure (fullSplitGraph, splitGraphCost)


{- | reoptimizeSplitGraphFromVertexIAOrig performs operations of reoptimizeSplitGraphFromVertex for static charcaters
but dynamic characters--only update IA assignments and initialized from origPhylo graph (at leaves) to keep IA characters in sync
since all "static" only need single traversal post order pass

uses PhylogenetiGraph internally
-}
reoptimizeSplitGraphFromVertexIANew
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → VertexCost
    → PhylogeneticGraph
    → DecoratedGraph
    → LG.Node
    → LG.Node
    → PhyG (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertexIANew swapParams inGS inData netPenaltyFactor curGraph inSplitGraph startVertex prunedSubGraphRootVertex =
    -- if graphType inGS /= Tree then error "Networks not yet implemented in reoptimizeSplitGraphFromVertexIA"
    -- else
    let nonExactCharacters = U.getNumberSequenceCharacters (thd3 inData)
        origGraph = inSplitGraph -- thd5 origPhyloGraph

        -- create leaf graphs--but copy IA final to prelim
        leafGraph =
            if graphType inGS == SoftWired
                then GO.copyIAFinalToPrelim $ LG.extractLeafGraph origGraph -- POSW.makeLeafGraphSoftWired inGS inData -- LG.extractLeafGraph origGraph
                else GO.copyIAFinalToPrelim $ LG.extractLeafGraph origGraph
        calcBranchLengths = False

        -- this for multitravers in swap for softwired to turn off
        multiTraverse =
            if graphType inGS == Tree
                then multiTraverseCharacters inGS
                else False

        -- create simple graph version of split for post order pass
        splitGraphSimple = GO.convertDecoratedToSimpleGraph inSplitGraph

        parentPoint = L.uncons $ concat $ fmap (LG.parents (thd6 curGraph)) $ LG.parents (thd6 curGraph) prunedSubGraphRootVertex
                
        parentOriginalConnection = if isNothing parentPoint then error ("Split graph involving near-root edges: " <> (show prunedSubGraphRootVertex))
                                           else if (not . null . snd) $ fromJust parentPoint then error ("Split graph involving network edge: " <> (show (prunedSubGraphRootVertex, snd $ fromJust parentPoint)))
                                           else fst $ fromJust parentPoint
    in  do
            -- Create base graph
            -- create postorder assignment--but only from single traversal
            -- True flag fior staticIA
            postOrderBaseGraph ←
                POSW.postOrderTreeTraversal
                    (inGS{graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse})
                    inData
                    (Just (thd6 curGraph, parentOriginalConnection))
                    leafGraph
                    True
                    (Just startVertex)
                    splitGraphSimple
            let baseGraphCost = snd6 postOrderBaseGraph

            -- True flag fior staticIA
            fullBaseGraph ←
                PRE.preOrderTreeTraversal
                    (inGS{graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse})
                    (finalAssignment inGS)
                    True
                    calcBranchLengths
                    (nonExactCharacters > 0)
                    startVertex
                    True
                    postOrderBaseGraph

            {-
            localRootCost = if (rootCost inGS) == NoRootCost then 0.0
                              else error ("Root cost type " <> (show $ rootCost inGS) <> " is not yet implemented")
            -}

            -- get root node of base graph
            let startBaseNode = (startVertex, fromJust $ LG.lab (thd6 fullBaseGraph) startVertex)

            -- Create pruned graph
            -- get root node of pruned graph--parent since that is the full pruned piece (keeping that node for addition to base graph and edge creation)
            let startPrunedNode = GO.makeIAPrelimFromFinal (prunedSubGraphRootVertex, fromJust $ LG.lab origGraph prunedSubGraphRootVertex)
            let startPrunedParentNode = head $ LG.labParents origGraph prunedSubGraphRootVertex
            let startPrunedParentEdge = (fst startPrunedParentNode, prunedSubGraphRootVertex, dummyEdge)

            -- True flag fior staticIA
            let postOrderPrunedGraph = if LG.isLeaf origGraph prunedSubGraphRootVertex then curGraph
                                               else 
                                                    -- updating pruned root final to preliminary for preorder pass
                                                    let origLabelPrunedRoot = fromJust $ (LG.lab (thd6 curGraph) prunedSubGraphRootVertex)
                                                        originalVertInfoPrunedRoot = vertData origLabelPrunedRoot
                                                        prunedGraphCost = subGraphCost origLabelPrunedRoot
                                                        newVertInfoPrunedRoot = PRE.setFinalToPreliminaryStates originalVertInfoPrunedRoot
                                                        newFinalPrunedRoot = (prunedSubGraphRootVertex, origLabelPrunedRoot {vertData = newVertInfoPrunedRoot})
                                                        newPrunedRootOutEdges = LG.out (thd6 curGraph) prunedSubGraphRootVertex
                                                        newPrunedGraph =  LG.insEdges newPrunedRootOutEdges $ LG.insNode newFinalPrunedRoot $ LG.delNode prunedSubGraphRootVertex (thd6 curGraph)
                                                    in (fst6 curGraph, prunedGraphCost, newPrunedGraph, fth6 curGraph, fft6 curGraph, six6 curGraph)
            let prunedGraphCost = snd6 postOrderPrunedGraph
            
            -- True flag fior staticIA
            fullPrunedGraph ← if swapType swapParams `elem` [NNI, SPR] then 
                                            pure postOrderPrunedGraph
                                      else          
                                            PRE.preOrderTreeTraversal
                                                (inGS{graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse})
                                                (finalAssignment inGS)
                                                True
                                                calcBranchLengths
                                                (nonExactCharacters > 0)
                                                prunedSubGraphRootVertex
                                                True
                                                postOrderPrunedGraph

            -- get nodes and edges in base and pruned graph (both PhylogeneticGrapgs so thd5)
            let (baseGraphNonRootNodes, baseGraphEdges) = LG.nodesAndEdgesAfter (thd6 fullBaseGraph) [startBaseNode]

            let (prunedGraphNonRootNodes, prunedGraphEdges) =
                    if LG.isLeaf origGraph prunedSubGraphRootVertex
                        then ([], [])
                        else LG.nodesAndEdgesAfter (thd6 fullPrunedGraph) [startPrunedNode]

            -- make fully optimized graph from base and split components
            let fullSplitGraph =
                    LG.mkGraph
                        ([startBaseNode, startPrunedNode, startPrunedParentNode] <> baseGraphNonRootNodes <> prunedGraphNonRootNodes)
                        (startPrunedParentEdge : (baseGraphEdges <> prunedGraphEdges))

            let splitGraphCost = ((1.0 + netPenaltyFactor) * (baseGraphCost + prunedGraphCost))

            -- remove when working
            -- trace ("ROGFVIA split costs:" <> (show (baseGraphCost, prunedGraphCost, localRootCost)) <> " -> " <> (show splitGraphCost)) (
            if prunedGraphCost == infinity || baseGraphCost == infinity
                then pure (LG.empty, infinity)
                else
                    if splitGraphCost == 0
                        then
                            error
                                ( "Split costs:"
                                    <> (show (baseGraphCost, prunedGraphCost))
                                    <> " -> "
                                    <> (show splitGraphCost)
                                    <> " Split graph simple:\n"
                                    <> (LG.prettify splitGraphSimple)
                                    <> "\nFull:\n"
                                    <> (show inSplitGraph)
                                    <> "\nOriginal Graph:\n"
                                    <> (show origGraph)
                                )
                        else pure (fullSplitGraph, splitGraphCost)

-- | reoptimizeSplitGraphFromVertexTupleNew wrapper for reoptimizeSplitGraphFromVertex with last 3 args as tuple
reoptimizeSplitGraphFromVertexTupleNew
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Bool
    → VertexCost
    → (PhylogeneticGraph, DecoratedGraph, LG.Node, LG.Node)
    → PhyG (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertexTupleNew swapParams inGS inData doIA netPenaltyFactor (curGraph, inSplitGraph, startVertex, prunedSubGraphRootVertex) =
    reoptimizeSplitGraphFromVertexNew swapParams inGS inData doIA netPenaltyFactor curGraph inSplitGraph startVertex prunedSubGraphRootVertex 



-- | reoptimizeSplitGraphFromVertexTupleOrig wrapper for reoptimizeSplitGraphFromVertex with last 3 args as tuple
reoptimizeSplitGraphFromVertexTupleOrig
    ∷ GlobalSettings
    → ProcessedData
    → Bool
    → VertexCost
    → (DecoratedGraph, Int, Int)
    → PhyG (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertexTupleOrig inGS inData doIA netPenaltyFactor (inSplitGraph, startVertex, prunedSubGraphRootVertex) =
    reoptimizeSplitGraphFromVertexOrig inGS inData doIA netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex


{- | reoptimizeSplitGraphFromVertex 
    Original version of reoptimizeSplitGraphFromVertex from Swap
    copied here for testing convenience.
    roughly O(n)
-}
reoptimizeSplitGraphFromVertexOrig
    ∷ GlobalSettings
    → ProcessedData
    → Bool
    → VertexCost
    → DecoratedGraph
    → Int
    → Int
    → PhyG (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertexOrig inGS inData doIA netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex =
    -- trace ("RSGFV: " <> (show startVertex)) (
    if doIA
        then -- only reoptimize the IA states for dynamic characters
            reoptimizeSplitGraphFromVertexIAOrig inGS inData netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex
        else -- perform full optimizations of nodes
        -- these required for full optimization

            let nonExactCharacters = U.getNumberSequenceCharacters (thd3 inData)
                origGraph = inSplitGraph -- thd5 origPhyloGraph
                leafGraph =
                    if graphType inGS == SoftWired
                        then LG.extractLeafGraph origGraph -- POSW.makeLeafGraphSoftWired inGS inData -- LG.extractLeafGraph origGraph
                        else LG.extractLeafGraph origGraph
                calcBranchLengths = False

                -- this for multitravers in swap for softwired to turn off
                multiTraverse =
                    if graphType inGS /= HardWired
                        then multiTraverseCharacters inGS
                        else False

                -- create simple graph version of split for post order pass
                splitGraphSimple = GO.convertDecoratedToSimpleGraph inSplitGraph
            in  do
                    -- create optimized base graph
                    -- False for staticIA
                    (postOrderBaseGraph, _) ←
                        T.generalizedGraphPostOrderTraversal
                            (inGS{graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse})
                            nonExactCharacters
                            inData
                            Nothing
                            leafGraph
                            False
                            (Just startVertex)
                            splitGraphSimple

                    fullBaseGraph ←
                        PRE.preOrderTreeTraversal
                            (inGS{graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse})
                            (finalAssignment inGS)
                            False
                            calcBranchLengths
                            (nonExactCharacters > 0)
                            startVertex
                            True
                            postOrderBaseGraph

                    -- create fully optimized pruned graph.  Post order then preorder

                    -- get root node of pruned graph--parent since that is the full pruned piece (keeping that node for addition to base graph and edge creation)
                    let startPrunedNode = (prunedSubGraphRootVertex, fromJust $ LG.lab origGraph prunedSubGraphRootVertex)
                    let startPrunedParentNode = head $ LG.labParents origGraph prunedSubGraphRootVertex
                    let startPrunedParentEdge = (fst startPrunedParentNode, prunedSubGraphRootVertex, dummyEdge)

                    -- False for staticIA
                    (postOrderPrunedGraph, _) ←
                        T.generalizedGraphPostOrderTraversal
                            (inGS{graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse})
                            nonExactCharacters
                            inData
                            Nothing
                            leafGraph
                            False
                            (Just prunedSubGraphRootVertex)
                            splitGraphSimple

                    -- False for staticIA
                    fullPrunedGraph ←
                        PRE.preOrderTreeTraversal
                            (inGS{graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse})
                            (finalAssignment inGS)
                            False
                            calcBranchLengths
                            (nonExactCharacters > 0)
                            prunedSubGraphRootVertex
                            True
                            postOrderPrunedGraph

                    -- get root node of base graph
                    let startBaseNode = (startVertex, fromJust $ LG.lab (thd6 fullBaseGraph) startVertex)

                    -- get nodes and edges in base and pruned graph (both PhylogeneticGrapgs so thd5)
                    let (baseGraphNonRootNodes, baseGraphEdges) = LG.nodesAndEdgesAfter (thd6 fullBaseGraph) [startBaseNode]

                    let (prunedGraphNonRootNodes, prunedGraphEdges) =
                            if LG.isLeaf origGraph prunedSubGraphRootVertex
                                then ([], [])
                                else LG.nodesAndEdgesAfter (thd6 fullPrunedGraph) [startPrunedNode]

                    -- make fully optimized graph from base and split components
                    let fullSplitGraph =
                            LG.mkGraph
                                ([startBaseNode, startPrunedNode, startPrunedParentNode] <> baseGraphNonRootNodes <> prunedGraphNonRootNodes)
                                (startPrunedParentEdge : (baseGraphEdges <> prunedGraphEdges))

                    -- cost of split graph to be later combined with re-addition delta for heuristic graph cost
                    let prunedCost =
                            if LG.isLeaf origGraph prunedSubGraphRootVertex
                                then 0
                                else snd6 fullPrunedGraph
                    let splitGraphCost = ((1.0 + netPenaltyFactor) * ((snd6 fullBaseGraph) + prunedCost))

                    {-
                    -- check fo unlabbeld nodes
                    coninicalNodes =  LG.labNodes fullSplitGraph
                    nodeLabels = fmap (LG.lab fullSplitGraph) (fmap fst coninicalNodes)
                    unlabelledNodes = filter ((== Nothing) .snd) $ (zip (fmap fst coninicalNodes) nodeLabels)
                    -}

                    if prunedCost == infinity || (snd6 fullBaseGraph) == infinity
                        then pure (LG.empty, infinity)
                        else pure (fullSplitGraph, splitGraphCost)

{- | reoptimizeSplitGraphFromVertexIAOrig performs operations of reoptimizeSplitGraphFromVertex for static charcaters
but dynamic characters--only update IA assignments and initialized from origPhylo graph (at leaves) to keep IA characters in sync
since all "static" only need single traversal post order pass

uses PhylogenetiGraph internally
-}
reoptimizeSplitGraphFromVertexIAOrig
    ∷ GlobalSettings
    → ProcessedData
    → VertexCost
    → DecoratedGraph
    → Int
    → Int
    → PhyG (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertexIAOrig inGS inData netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex =
    -- if graphType inGS /= Tree then error "Networks not yet implemented in reoptimizeSplitGraphFromVertexIA"
    -- else
    let nonExactCharacters = U.getNumberSequenceCharacters (thd3 inData)
        origGraph = inSplitGraph -- thd5 origPhyloGraph

        -- create leaf graphs--but copy IA final to prelim
        leafGraph =
            if graphType inGS == SoftWired
                then GO.copyIAFinalToPrelim $ LG.extractLeafGraph origGraph -- POSW.makeLeafGraphSoftWired inGS inData -- LG.extractLeafGraph origGraph
                else GO.copyIAFinalToPrelim $ LG.extractLeafGraph origGraph
        calcBranchLengths = False

        -- this for multitravers in swap for softwired to turn off
        multiTraverse =
            if graphType inGS == Tree
                then multiTraverseCharacters inGS
                else False

        -- create simple graph version of split for post order pass
        splitGraphSimple = GO.convertDecoratedToSimpleGraph inSplitGraph
    in  do
            -- Create base graph
            -- create postorder assignment--but only from single traversal
            -- True flag fior staticIA
            postOrderBaseGraph ←
                POSW.postOrderTreeTraversal
                    (inGS{graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse})
                    inData
                    Nothing
                    leafGraph
                    True
                    (Just startVertex)
                    splitGraphSimple
            let baseGraphCost = snd6 postOrderBaseGraph

            -- True flag fior staticIA
            fullBaseGraph ←
                PRE.preOrderTreeTraversal
                    (inGS{graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse})
                    (finalAssignment inGS)
                    True
                    calcBranchLengths
                    (nonExactCharacters > 0)
                    startVertex
                    True
                    postOrderBaseGraph

            {-
            localRootCost = if (rootCost inGS) == NoRootCost then 0.0
                              else error ("Root cost type " <> (show $ rootCost inGS) <> " is not yet implemented")
            -}

            -- get root node of base graph
            let startBaseNode = (startVertex, fromJust $ LG.lab (thd6 fullBaseGraph) startVertex)

            -- Create pruned graph
            -- get root node of pruned graph--parent since that is the full pruned piece (keeping that node for addition to base graph and edge creation)
            let startPrunedNode = GO.makeIAPrelimFromFinal (prunedSubGraphRootVertex, fromJust $ LG.lab origGraph prunedSubGraphRootVertex)
            let startPrunedParentNode = head $ LG.labParents origGraph prunedSubGraphRootVertex
            let startPrunedParentEdge = (fst startPrunedParentNode, prunedSubGraphRootVertex, dummyEdge)

            -- True flag fior staticIA
            postOrderPrunedGraph ←
                POSW.postOrderTreeTraversal
                    (inGS{graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse})
                    inData
                    Nothing
                    leafGraph
                    True
                    (Just prunedSubGraphRootVertex)
                    splitGraphSimple
            let prunedGraphCost = snd6 postOrderPrunedGraph

            -- True flag fior staticIA
            fullPrunedGraph ←
                PRE.preOrderTreeTraversal
                    (inGS{graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse})
                    (finalAssignment inGS)
                    True
                    calcBranchLengths
                    (nonExactCharacters > 0)
                    prunedSubGraphRootVertex
                    True
                    postOrderPrunedGraph

            -- get nodes and edges in base and pruned graph (both PhylogeneticGrapgs so thd5)
            let (baseGraphNonRootNodes, baseGraphEdges) = LG.nodesAndEdgesAfter (thd6 fullBaseGraph) [startBaseNode]

            let (prunedGraphNonRootNodes, prunedGraphEdges) =
                    if LG.isLeaf origGraph prunedSubGraphRootVertex
                        then ([], [])
                        else LG.nodesAndEdgesAfter (thd6 fullPrunedGraph) [startPrunedNode]

            -- make fully optimized graph from base and split components
            let fullSplitGraph =
                    LG.mkGraph
                        ([startBaseNode, startPrunedNode, startPrunedParentNode] <> baseGraphNonRootNodes <> prunedGraphNonRootNodes)
                        (startPrunedParentEdge : (baseGraphEdges <> prunedGraphEdges))

            let splitGraphCost = ((1.0 + netPenaltyFactor) * (baseGraphCost + prunedGraphCost))

            -- remove when working
            -- trace ("ROGFVIA split costs:" <> (show (baseGraphCost, prunedGraphCost, localRootCost)) <> " -> " <> (show splitGraphCost)) (
            if prunedGraphCost == infinity || baseGraphCost == infinity
                then pure (LG.empty, infinity)
                else
                    if splitGraphCost == 0
                        then
                            error
                                ( "Split costs:"
                                    <> (show (baseGraphCost, prunedGraphCost))
                                    <> " -> "
                                    <> (show splitGraphCost)
                                    <> " Split graph simple:\n"
                                    <> (LG.prettify splitGraphSimple)
                                    <> "\nFull:\n"
                                    <> (show inSplitGraph)
                                    <> "\nOriginal Graph:\n"
                                    <> (show origGraph)
                                )
                        else pure (fullSplitGraph, splitGraphCost)
