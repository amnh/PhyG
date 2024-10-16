{- |
Module specifying graph traversal functions for PhyGraph
-}
module GraphOptimization.Traversals (
    multiTraverseFullyLabelTree,
    multiTraverseFullyLabelGraph,
    multiTraverseFullyLabelGraph',
    multiTraverseFullyLabelSoftWired,
    multiTraverseFullyLabelSoftWiredReduced,
    multiTraverseFullyLabelHardWired,
    multiTraverseFullyLabelHardWiredReduced,
    checkUnusedEdgesPruneInfty,
    generalizedGraphPostOrderTraversal,
    getPenaltyFactor,
    getPenaltyFactorList,
    updateGraphCostsComplexities,
    updatePhylogeneticGraphCost,
    updatePhylogeneticGraphCostReduced,
    multiTraverseFullyLabelGraphReduced,
    multiTraverseFullyLabelGraphPair,
) where

import Control.Monad (when)
import Data.InfList qualified as IL
import Data.List qualified as L
import Data.Maybe
import GeneralUtilities
import GraphOptimization.PostOrderSoftWiredFunctions qualified as POSW
import GraphOptimization.PreOrderFunctions qualified as PRE
import Graphs.GraphOperations qualified as GO
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging
-- import PHANE.Evaluation.Verbosity (Verbosity (..))
import Types.Types
import Utilities.LocalGraph qualified as LG
import Utilities.Utilities as U
--import Debug.Trace


{- |
'multiTraverseFullyLabelGraphReduced' wrapper to return 'ReducedPhylogeneticGraph'.
-}
multiTraverseFullyLabelGraphReduced
    ∷ GlobalSettings → ProcessedData → Bool → Bool → Maybe Int → SimpleGraph → PhyG ReducedPhylogeneticGraph
multiTraverseFullyLabelGraphReduced inGS inData pruneEdges warnPruneEdges startVertex inGraph =
    GO.convertPhylogeneticGraph2Reduced <$> multiTraverseFullyLabelGraph inGS inData pruneEdges warnPruneEdges startVertex inGraph



-- \$ multiTraverseFullyLabelGraph inGS inData pruneEdges warnPruneEdges startVertex inGraph

{- | multiTraverseFullyLabelGraph is a wrapper around multi-traversal functions for Tree,
Soft-wired network graph, and Hard-wired network graph
can either find root, be given root, or start somwhere else (startVertex) do optimize only a component of a forest
-}
multiTraverseFullyLabelGraph
    ∷ GlobalSettings → ProcessedData → Bool → Bool → Maybe Int → SimpleGraph → PhyG PhylogeneticGraph
multiTraverseFullyLabelGraph inGS inData pruneEdges warnPruneEdges startVertex inGraph
    | LG.isEmpty inGraph = pure emptyPhylogeneticGraph
    | otherwise =
        {-# SCC multiTraverseFullyLabelGraph_TOP_DEF #-}
        case graphType inGS of
            SoftWired →
                let leafGraph = POSW.makeLeafGraphSoftWired inGS inData
                in  multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex inGraph
            HardWired →
                let leafGraph = GO.makeLeafGraph inData
                in  multiTraverseFullyLabelHardWired inGS inData leafGraph startVertex inGraph
            Tree →
                {-# SCC multiTraverseFullyLabelGraph_CASE_OF_Tree #-}
                -- test for Tree
                let (_, _, _, networkVertexList) = LG.splitVertexList inGraph
                in  --            in  do  when (not $ null networkVertexList) . failWithPhase Computing $ unlines
                    do
                        when (not $ null networkVertexList) $ do
                            logWith LogFail $
                                unlines
                                    [ "Input graph is not a tree/forest, but graph type has been specified (perhaps by default) as Tree."
                                    , "Modify input graph or use 'set()' command to specify network type."
                                    , "\tNetwork vertices: " <> show (fst <$> networkVertexList)
                                    , LG.prettify inGraph
                                    ]
                            error "Exceptional state reached: MALFORMED GRAPH"

                        let leafGraph = GO.makeLeafGraph inData
                        multiTraverseFullyLabelTree inGS inData leafGraph startVertex inGraph
            val → failWithPhase Computing $ "Unknown graph type specified: " <> show val


-- | multiTraverseFullyLabelGraphPair maps to multiTraverseFullyLabelGraph with paired  arguments used by report IA and tnt output
multiTraverseFullyLabelGraphPair
    ∷ GlobalSettings → Bool → Bool → Maybe Int → (ProcessedData, SimpleGraph) → PhyG PhylogeneticGraph
multiTraverseFullyLabelGraphPair inGS pruneEdges warnPruneEdges startVertex (inData, inGraph) = multiTraverseFullyLabelGraph inGS inData pruneEdges warnPruneEdges startVertex inGraph


-- | multiTraverseFullyLabelGraph' maps to multiTraverseFullyLabelGraph with differnet order of arguments used by report IA and tnt output
multiTraverseFullyLabelGraph'
    ∷ GlobalSettings → Bool → Bool → Maybe Int → ProcessedData → SimpleGraph → PhyG PhylogeneticGraph
multiTraverseFullyLabelGraph' inGS pruneEdges warnPruneEdges startVertex inData inGraph = multiTraverseFullyLabelGraph inGS inData pruneEdges warnPruneEdges startVertex inGraph


-- | multiTraverseFullyLabelHardWiredReduced is a wrapper for ReducedPhylogeneticTrees
multiTraverseFullyLabelHardWiredReduced
    ∷ GlobalSettings → ProcessedData → DecoratedGraph → Maybe Int → SimpleGraph → PhyG ReducedPhylogeneticGraph
multiTraverseFullyLabelHardWiredReduced inGS inData leafGraph startVertex inSimpleGraph = do
    result ← multiTraverseFullyLabelTree inGS inData leafGraph startVertex inSimpleGraph
    pure $ GO.convertPhylogeneticGraph2Reduced result


multiTraverseFullyLabelHardWired
    ∷ GlobalSettings → ProcessedData → DecoratedGraph → Maybe Int → SimpleGraph → PhyG PhylogeneticGraph
multiTraverseFullyLabelHardWired inGS inData leafGraph startVertex inSimpleGraph = multiTraverseFullyLabelTree inGS inData leafGraph startVertex inSimpleGraph


-- | multiTraverseFullyLabelSoftWiredReduced is a wrapper for ReducedPhylogeneticTrees
multiTraverseFullyLabelSoftWiredReduced
    ∷ GlobalSettings → ProcessedData → Bool → Bool → DecoratedGraph → Maybe Int → SimpleGraph → PhyG ReducedPhylogeneticGraph
multiTraverseFullyLabelSoftWiredReduced inGS inData pruneEdges warnPruneEdges leafGraph startVertex inSimpleGraph = do
    result ← multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex inSimpleGraph
    pure $ GO.convertPhylogeneticGraph2Reduced result


{- | multiTraverseFullyLabelSoftWired fully labels a softwired network component forest
including traversal rootings-- does not reroot on network edges
allows indegree=outdegree=1 vertices
pruneEdges and warnPruneEdges specify if unused edges (ie not in diuaplytrees) are pruned from
canonical tree or if an infinity cost is returned and if a trace warning is thrown if so.
in general--input trees should use "pruneEdges" during search--not
can either find root, be given root, or start somwhere else (startVertex) do optimize only a component of a forest
first Bool for calcualting breanch edger weights
-}
multiTraverseFullyLabelSoftWired
    ∷ GlobalSettings → ProcessedData → Bool → Bool → DecoratedGraph → Maybe Int → SimpleGraph → PhyG PhylogeneticGraph
multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex inSimpleGraph =
    if LG.isEmpty inSimpleGraph
        then pure emptyPhylogeneticGraph
        else do
            let sequenceChars = U.getNumberSequenceCharacters (thd3 inData)
            (postOrderGraph, localStartVertex) ←
                generalizedGraphPostOrderTraversal inGS sequenceChars inData Nothing leafGraph False startVertex inSimpleGraph
            fullyOptimizedGraph ←
                PRE.preOrderTreeTraversal inGS (finalAssignment inGS) False True (sequenceChars > 0) localStartVertex False postOrderGraph
            checkUnusedEdgesPruneInfty inGS inData pruneEdges warnPruneEdges leafGraph $
                updatePhylogeneticGraphCost fullyOptimizedGraph (snd6 fullyOptimizedGraph)


{- | multiTraverseFullyLabelTree performs potorder on default root and other traversal foci, taking the minimum
traversal cost for all nonexact charcters--the initial rooting is used for exact characters
operates with Tree functions
need to add forest functionality--in principle just split into components and optimize them independently
but get into root index issues the way this is written now.
can either find root, be given root, or start somwhere else (startVertex) do optimize only a component of a forest
-}
multiTraverseFullyLabelTree
    ∷ GlobalSettings → ProcessedData → DecoratedGraph → Maybe Int → SimpleGraph → PhyG PhylogeneticGraph
multiTraverseFullyLabelTree inGS inData leafGraph startVertex inSimpleGraph =
    if LG.isEmpty inSimpleGraph
        then pure emptyPhylogeneticGraph
        else
            let sequenceChars = U.getNumberSequenceCharacters (thd3 inData)
                staticIA = False
            in  do
                    (postOrderGraph, localStartVertex) ←
                        generalizedGraphPostOrderTraversal inGS sequenceChars inData Nothing leafGraph staticIA startVertex inSimpleGraph
                    PRE.preOrderTreeTraversal inGS (finalAssignment inGS) False True (sequenceChars > 0) localStartVertex False postOrderGraph


{- | generalizedGraphPostOrderTraversal performs the postorder pass
on a graph (tree, softWired, or hardWired) to determine the "preliminary" character states
include penalty factor cost but not root cost which may or may not be wanted depending on context
if full graph--yes, if a component yes or no.
hence returns the pair
Adds root complexity cost if the start vertex is Nothing (e.g. the graph root), so should be corrext when graphs are split in swap.

Extra input MaybeGraph (incrementalGraph) for incremental optimization for initial post-order 
(no need for reroots--already contant time for each edge)
-}
generalizedGraphPostOrderTraversal
    ∷ GlobalSettings → Int → ProcessedData → Maybe (DecoratedGraph, LG.Node) → DecoratedGraph → Bool → Maybe Int → SimpleGraph → PhyG (PhylogeneticGraph, Int)
generalizedGraphPostOrderTraversal inGS sequenceChars inData incrementalGraph leafGraph staticIA startVertex inSimpleGraph = do
    -- next edges (to vertex in list) to perform rerooting
    -- progresses recursivey over adjacent edges to minimize node reoptimization
    -- childrenOfRoot = concatMap (LG.descendants (thd6 outgroupRooted)) startVertexList
    -- grandChildrenOfRoot = concatMap (LG.descendants (thd6 outgroupRooted)) childrenOfRoot

    -- create list of multi-traversals with original rooting first
    -- subsequent rerooting do not reoptimize exact characters (add nonadd etc)
    -- they are taken from the fully labelled first input decorated graph later when output graph created
    -- it is important that the first graph be the outgroup rooted graph (outgroupRootedPhyloGraph) so this
    -- will have the preorder assignments for the outgroup rooted graph as 3rd field.  This can be used for incremental
    -- optimization to get O(log n) initial postorder assingment when mutating graph.
    -- hardwired reroot cause much pain
    -- the head startvertex list for reoptimizing spit trees in swapping

    -- create optimal final graph with best costs and best traversal (rerooting) forest for each character
    -- traversal for exact characters (and costs) are the first of each least since exact only optimizaed for that
    -- traversal graph.  The result has approprotate post-order assignments for traversals, preorder "final" assignments
    -- are propagated to the Decorated graph field after the preorder pass.
    -- doesn't have to be sorted, but should minimize assignments
    -- graphWithBestAssignments = head recursiveRerootList -- L.foldl1' setBetterGraphAssignment recursiveRerootList'

    {-  root and model complexities moved to output
    -- same root cost if same data and number of roots
    localRootCost = rootComplexity inGS
                    {-if (rootCost inGS) == NoRootCost then 0.0
                    else error ("Root cost type " <> (show $ rootCost inGS) <> " is not yet implemented")
                    -}
    -}
    -- let (rootNodes, leafNode, treeNodes,networkNodes) = LG.splitVertexList inSimpleGraph
    -- logWith LogInfo ("GPOT: " <> (show (length rootNodes, length leafNode, length treeNodes, length networkNodes)))

    -- first traversal on outgroup root
    outgroupRooted ←
        if (graphType inGS) `elem` [Tree, HardWired]
            then POSW.postOrderTreeTraversal inGS inData incrementalGraph leafGraph staticIA startVertex inSimpleGraph
            else
                if (graphType inGS) == SoftWired
                    then POSW.postOrderSoftWiredTraversal inGS inData incrementalGraph leafGraph staticIA startVertex inSimpleGraph
                    else error ("Graph type not implemented: " <> (show $ graphType inGS))

    -- start at start vertex--for components or ur-root for full graph
    -- root cost for overall cost is only added for thopse greaophs that include overall root
    -- model complexity for PMDL is also added here
    let (startVertexList, rootAndModelCost) =
            if isJust startVertex
                then ([fromJust startVertex], 0)
                else
                    if optimalityCriterion inGS == PMDL
                        then (fmap fst $ LG.getRoots $ thd6 outgroupRooted, rootComplexity inGS + modelComplexity inGS)
                        else (fmap fst $ LG.getRoots $ thd6 outgroupRooted, rootComplexity inGS)

    -- only static characters
    if sequenceChars == 0
        then do
            penaltyFactor ← getPenaltyFactor inGS inData startVertex outgroupRooted
            
            staticOnlyGraph ←
                if (graphType inGS) == SoftWired
                    then POSW.updateAndFinalizePostOrderSoftWired startVertex (head startVertexList) outgroupRooted
                    else pure outgroupRooted
            
            let staticOnlyGraph' =
                    if startVertex == Nothing
                        then updatePhylogeneticGraphCost staticOnlyGraph (rootAndModelCost + penaltyFactor + (snd6 staticOnlyGraph))
                        else updatePhylogeneticGraphCost staticOnlyGraph (rootAndModelCost + penaltyFactor + (snd6 staticOnlyGraph))
            pure (staticOnlyGraph', head startVertexList)
        else do
            recursiveRerootList ←
                if (graphType inGS == HardWired)
                    then pure [outgroupRooted]
                    else
                        if (graphType inGS == SoftWired)
                            then do
                                displayResult ← POSW.getDisplayBasedRerootSoftWired inGS SoftWired (head startVertexList) outgroupRooted
                                pure [displayResult]
                            else
                                if (graphType inGS == Tree)
                                    then do
                                        displayResult ← POSW.getDisplayBasedRerootSoftWired inGS Tree (head startVertexList) outgroupRooted
                                        -- logWith LogInfo ("RRL: " <> (show (snd6 displayResult, snd6 outgroupRooted)) <> "\n")
                                        pure [displayResult]
                                    else error ("Graph type not implemented: " <> (show $ graphType inGS))

            -- single sequence (prealigned, dynamic) only (ie no static)

            if sequenceChars == 1 && (U.getNumberExactCharacters (thd3 inData) == 0)
                then do
                    let finalizedPostOrderGraphList = L.sortOn snd6 recursiveRerootList
                    penaltyFactorList ← getPenaltyFactorList inGS inData startVertex finalizedPostOrderGraphList
                    
                    let newCostList =
                            zipWith3
                                sum3
                                penaltyFactorList
                                (fmap snd6 finalizedPostOrderGraphList)
                                (replicate (length finalizedPostOrderGraphList) rootAndModelCost)

                    let finalizedPostOrderGraph = head $ L.sortOn snd6 $ zipWith updatePhylogeneticGraphCost finalizedPostOrderGraphList newCostList

                    pure (finalizedPostOrderGraph, head startVertexList)
                else do
                    -- multiple dynamic characters--checks for best root for each character
                    -- important to have outgroup rooted graph first for fold so don't use sorted recursive list
                    let graphWithBestAssignments = head recursiveRerootList
                    penaltyFactor ← getPenaltyFactor inGS inData startVertex graphWithBestAssignments
                    
                    let graphWithBestAssignments' = updatePhylogeneticGraphCost graphWithBestAssignments (rootAndModelCost + penaltyFactor + (snd6 graphWithBestAssignments))
                    -- logWith LogInfo ("GGPOT: " <> (show (snd6 graphWithBestAssignments', snd6 graphWithBestAssignments)) <> "\n")
                    pure (graphWithBestAssignments', head startVertexList)
    where
        sum3 ∷ VertexCost → VertexCost → VertexCost → VertexCost
        sum3 a b c = a + b + c


{- | getPenaltyFactorList takes graph type and other options and a list of graphs to return a
    list of penalty factors
-}
getPenaltyFactorList ∷ GlobalSettings → ProcessedData → Maybe Int → [PhylogeneticGraph] → PhyG [VertexCost]
getPenaltyFactorList inGS inData startVertex inGraphList =
    if null inGraphList
        then pure []
        else mapM (getPenaltyFactor inGS inData startVertex) inGraphList


{- | getPenaltyFactor take graph information and optimality criterion and returns penalty
or graph complexity for single graph
-}
getPenaltyFactor ∷ GlobalSettings → ProcessedData → Maybe Int → PhylogeneticGraph → PhyG VertexCost
getPenaltyFactor inGS inData startVertex inGraph =
    if LG.isEmpty $ fst6 inGraph
        then pure 0.0
        else
            if (optimalityCriterion inGS) `elem` [PMDL, SI]
                then
                    let (_, _, _, networkNodeList) = LG.splitVertexList (fst6 inGraph)
                    in  if graphType inGS == Tree
                            then pure $ fst $ IL.head (graphComplexityList inGS)
                            else
                                if graphType inGS == SoftWired
                                    then pure $ fst $ (graphComplexityList inGS) IL.!!! (length networkNodeList)
                                    else
                                        if graphType inGS == HardWired
                                            then pure $ snd $ (graphComplexityList inGS) IL.!!! (length networkNodeList)
                                            else error ("Graph type " <> (show (graphType inGS)) <> " is not yet implemented")
                else
                    if graphType inGS == Tree
                        then pure 0.0
                        else
                            if (graphFactor inGS) == NoNetworkPenalty
                                then pure 0.0
                                else
                                    if (graphFactor inGS) == Wheeler2015Network
                                        then POSW.getW15NetPenaltyFull Nothing inGS inData startVertex inGraph
                                        else
                                            if (graphFactor inGS) == Wheeler2023Network
                                                then pure $ POSW.getW23NetPenalty inGraph
                                                else error ("Network penalty type " <> (show (graphFactor inGS)) <> " is not yet implemented")


{- | checkUnusedEdgesPruneInfty checks if a softwired phylogenetic graph has
"unused" edges sensu Wheeler 2015--that an edge in the canonical graph is
not present in any of the block display trees (that are heurstically optimal)
the options specify if the cost returned is Infinity (really max bound Double)
with no pruning of edges or the cost is left unchanged and unused edges are
pruned from the canonical graph
this is unDirected due to rerooting heuristic in post/preorder optimization
inifinity defined in Types.hs
-}
checkUnusedEdgesPruneInfty
    ∷ GlobalSettings → ProcessedData → Bool → Bool → DecoratedGraph → PhylogeneticGraph → PhyG PhylogeneticGraph
checkUnusedEdgesPruneInfty inGS inData pruneEdges warnPruneEdges leafGraph inGraph@(inSimple, _, inCanonical, blockTreeV, charTreeVV, charInfoVV) =
    let simpleEdgeList = LG.edges inSimple
        displayEdgeSet = L.nubBy LG.undirectedEdgeEquality $ concat $ concat $ fmap (fmap LG.edges) blockTreeV
        unusedEdges = LG.undirectedEdgeMinus simpleEdgeList displayEdgeSet
    in  -- no unused edges all OK
        if null unusedEdges
            then pure inGraph
            else -- unused edges--do not prune return "infinite cost"

                if not pruneEdges
                    then -- trace ("Unused edge->Infinity")
                        pure (inSimple, infinity, inCanonical, blockTreeV, charTreeVV, charInfoVV)
                    else -- unused but pruned--need to prune nodes and reoptimize to get final assignments correct

                        let newSimpleGraph = LG.delEdges unusedEdges inSimple
                            contractedSimple = GO.contractIn1Out1EdgesRename newSimpleGraph
                        in  if warnPruneEdges
                                then do
                                    -- too lazy to thread PhyG logging throuhg everything
                                    logWith LogWarn ("Pruning " <> (show $ length unusedEdges) <> " unused edges and reoptimizing graph")
                                    multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph Nothing contractedSimple
                                else multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph Nothing contractedSimple


{- | updateGraphCostsComplexities adds root and model complexities if appropriate to graphs
updates NCM with original data due to weights of bitpacking
-}
updateGraphCostsComplexities
    ∷ GlobalSettings → ProcessedData → ProcessedData → Bool → [ReducedPhylogeneticGraph] → PhyG [ReducedPhylogeneticGraph]
updateGraphCostsComplexities inGS reportingData processedData rediagnoseWithReportingData inGraphList =
    -- parallel setup
    -- trace ("UGCC: " <> (show (optimalityCriterion inGS, rootComplexity inGS))) $
    let traverseAction ∷ SimpleGraph → PhyG ReducedPhylogeneticGraph
        traverseAction = multiTraverseFullyLabelGraphReduced inGS reportingData False False Nothing
    in  if optimalityCriterion inGS /= NCM
            then do
                pure inGraphList
            else -- NCM re-check graph cost (not root or model) due to bit packing
            do
                updatedGraphList ←
                    if (reportingData == emptyProcessedData) || (not rediagnoseWithReportingData) || (not $ U.has4864PackedChars (thd3 processedData))
                        then -- trace ("\t\tCannot update cost with original data--skipping")
                            pure inGraphList
                        else
                            getParallelChunkTraverse >>= \pTraverse →
                                (traverseAction . fst5) `pTraverse` inGraphList

                logWith LogInfo ("\tFinalizing graph cost (updating NCM)" <> "\n")
                pure updatedGraphList


{- |
'updatePhylogeneticGraphCostList' is a list wrapper for 'updatePhylogeneticGraphCost'.

updatePhylogeneticGraphCostList ∷ VertexCost → [ReducedPhylogeneticGraph] → [ReducedPhylogeneticGraph]
updatePhylogeneticGraphCostList rootCost inGraphList =
    fmap (updateCost rootCost) inGraphList
    where
        updateCost ∷ ∀ {b} {a} {c} {d} {e}. (Num b) ⇒ b → (a, b, c, d, e) → (a, b, c, d, e)
        updateCost z (a, oldCost, b, c, e) = (a, oldCost + z, b, c, e)
-}

{- |
'updatePhylogeneticGraphCost' takes a 'PhylgeneticGraph' and 'Double' and replaces the cost (snd of 6 fields)
and returns Phylogenetic graph.
-}
updatePhylogeneticGraphCost ∷ PhylogeneticGraph → VertexCost → PhylogeneticGraph
updatePhylogeneticGraphCost (a, _, b, c, d, e) newCost = (a, newCost, b, c, d, e)


{- |
'updatePhylogeneticGraphCost' takes a 'ReducedPhylogeneticGraph' and 'Double' and replaces the cost (snd of 6 fields)
and returns Phylogenetic graph.
-}
updatePhylogeneticGraphCostReduced ∷ ReducedPhylogeneticGraph → VertexCost → ReducedPhylogeneticGraph
updatePhylogeneticGraphCostReduced (a, _, b, c, e) newCost = (a, newCost, b, c, e)
