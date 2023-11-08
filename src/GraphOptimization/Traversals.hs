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

module GraphOptimization.Traversals ( multiTraverseFullyLabelTree
                                    , multiTraverseFullyLabelGraph
                                    , multiTraverseFullyLabelGraph'
                                    , multiTraverseFullyLabelSoftWired
                                    , multiTraverseFullyLabelSoftWiredReduced
                                    , multiTraverseFullyLabelHardWired
                                    , multiTraverseFullyLabelHardWiredReduced
                                    , checkUnusedEdgesPruneInfty
                                    , generalizedGraphPostOrderTraversal
                                    , updateGraphCostsComplexities
                                    , updatePhylogeneticGraphCost
                                    , updatePhylogeneticGraphCostReduced
                                    , multiTraverseFullyLabelGraphReduced
                                    , multiTraverseFullyLabelGraphPair
                                    ) where


import Control.Evaluation
import Control.Monad (when)
import Control.Monad.IO.Class (MonadIO (..))
import Control.Monad.Logger (LogLevel (..), Logger (..), Verbosity (..))
import Data.List qualified as L
import Data.Maybe
import GeneralUtilities
import GraphOptimization.PostOrderSoftWiredFunctions qualified as POSW
import GraphOptimization.PreOrderFunctions qualified as PRE
import Graphs.GraphOperations qualified as GO
import System.ErrorPhase (ErrorPhase (..))
import Types.Types
import Utilities.LocalGraph qualified as LG
import Utilities.Utilities as U

-- import ParallelUtilities qualified as PU
-- import Debug.Trace


-- | multiTraverseFullyLabelGraphReduced wrapper to return ReducedPhylogeneticGraph
multiTraverseFullyLabelGraphReduced :: GlobalSettings -> ProcessedData -> Bool -> Bool -> Maybe Int -> SimpleGraph -> PhyG ReducedPhylogeneticGraph
multiTraverseFullyLabelGraphReduced inGS inData pruneEdges warnPruneEdges startVertex inGraph = do
    result <- multiTraverseFullyLabelGraph inGS inData pruneEdges warnPruneEdges startVertex inGraph
    pure $ GO.convertPhylogeneticGraph2Reduced result -- $ multiTraverseFullyLabelGraph inGS inData pruneEdges warnPruneEdges startVertex inGraph

-- | multiTraverseFullyLabelGraph is a wrapper around multi-traversal functions for Tree,
-- Soft-wired network graph, and Hard-wired network graph
-- can either find root, be given root, or start somwhere else (startVertex) do optimize only a component of a forest
multiTraverseFullyLabelGraph :: GlobalSettings -> ProcessedData -> Bool -> Bool -> Maybe Int -> SimpleGraph -> PhyG PhylogeneticGraph
multiTraverseFullyLabelGraph inGS inData pruneEdges warnPruneEdges startVertex inGraph
  | LG.isEmpty inGraph = pure emptyPhylogeneticGraph
  | graphType inGS == Tree =
    -- test for Tree
    let (_, _, _, networkVertexList) = LG.splitVertexList inGraph
    in
    if null networkVertexList then
        let leafGraph = GO.makeLeafGraph inData
        in multiTraverseFullyLabelTree inGS inData leafGraph startVertex inGraph
    else errorWithoutStackTrace ("Input graph is not a tree/forest, but graph type has been specified (perhaps by default) as Tree. Modify input graph or use 'set()' command to specify network type\n\tNetwork vertices: " <> (show $ fmap fst networkVertexList) <> "\n" <> (LG.prettify inGraph))
  | graphType inGS == SoftWired =
    let leafGraph = POSW.makeLeafGraphSoftWired inGS inData
    in multiTraverseFullyLabelSoftWired  inGS inData pruneEdges warnPruneEdges leafGraph startVertex inGraph
  | graphType inGS == HardWired =
    let leafGraph = GO.makeLeafGraph inData
    in multiTraverseFullyLabelHardWired inGS inData leafGraph startVertex inGraph
  | otherwise = errorWithoutStackTrace ("Unknown graph type specified: " <> show (graphType inGS))


-- | multiTraverseFullyLabelGraphPair maps to multiTraverseFullyLabelGraph with paired  arguments used by report IA and tnt output
multiTraverseFullyLabelGraphPair :: GlobalSettings -> Bool -> Bool -> Maybe Int -> (ProcessedData, SimpleGraph) -> PhyG PhylogeneticGraph
multiTraverseFullyLabelGraphPair inGS pruneEdges warnPruneEdges startVertex (inData, inGraph) = multiTraverseFullyLabelGraph inGS inData pruneEdges warnPruneEdges startVertex inGraph

-- | multiTraverseFullyLabelGraph' maps to multiTraverseFullyLabelGraph with differnet order of arguments used by report IA and tnt output
multiTraverseFullyLabelGraph' :: GlobalSettings -> Bool -> Bool -> Maybe Int -> ProcessedData -> SimpleGraph -> PhyG PhylogeneticGraph
multiTraverseFullyLabelGraph' inGS pruneEdges warnPruneEdges startVertex inData inGraph = multiTraverseFullyLabelGraph inGS inData pruneEdges warnPruneEdges startVertex inGraph


-- | multiTraverseFullyLabelHardWiredReduced is a wrapper for ReducedPhylogeneticTrees
multiTraverseFullyLabelHardWiredReduced :: GlobalSettings -> ProcessedData -> DecoratedGraph -> Maybe Int -> SimpleGraph -> PhyG ReducedPhylogeneticGraph
multiTraverseFullyLabelHardWiredReduced inGS inData leafGraph startVertex inSimpleGraph = do
    result <- multiTraverseFullyLabelTree inGS inData leafGraph startVertex inSimpleGraph
    pure $ GO.convertPhylogeneticGraph2Reduced result 


multiTraverseFullyLabelHardWired :: GlobalSettings -> ProcessedData -> DecoratedGraph -> Maybe Int -> SimpleGraph -> PhyG PhylogeneticGraph
multiTraverseFullyLabelHardWired inGS inData leafGraph startVertex inSimpleGraph = multiTraverseFullyLabelTree inGS inData leafGraph startVertex inSimpleGraph


-- | multiTraverseFullyLabelSoftWiredReduced is a wrapper for ReducedPhylogeneticTrees
multiTraverseFullyLabelSoftWiredReduced :: GlobalSettings -> ProcessedData -> Bool -> Bool -> DecoratedGraph -> Maybe Int -> SimpleGraph -> PhyG ReducedPhylogeneticGraph
multiTraverseFullyLabelSoftWiredReduced inGS inData pruneEdges warnPruneEdges leafGraph startVertex inSimpleGraph = do
    result <-  multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex inSimpleGraph
    pure $ GO.convertPhylogeneticGraph2Reduced result

-- | multiTraverseFullyLabelSoftWired fully labels a softwired network component forest
-- including traversal rootings-- does not reroot on network edges
-- allows indegree=outdegree=1 vertices
-- pruneEdges and warnPruneEdges specify if unused edges (ie not in diuaplytrees) are pruned from
-- canonical tree or if an infinity cost is returned and if a trace warning is thrown if so.
-- in general--input trees should use "pruneEdges" during search--not
-- can either find root, be given root, or start somwhere else (startVertex) do optimize only a component of a forest
-- first Bool for calcualting breanch edger weights
multiTraverseFullyLabelSoftWired :: GlobalSettings -> ProcessedData -> Bool -> Bool -> DecoratedGraph -> Maybe Int -> SimpleGraph -> PhyG PhylogeneticGraph
multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex inSimpleGraph =
    if LG.isEmpty inSimpleGraph then pure emptyPhylogeneticGraph
    else do
        let sequenceChars = U.getNumberSequenceCharacters (thd3 inData)
        let (postOrderGraph, localStartVertex) = generalizedGraphPostOrderTraversal inGS sequenceChars inData leafGraph False startVertex inSimpleGraph
        fullyOptimizedGraph <- PRE.preOrderTreeTraversal inGS (finalAssignment inGS) False True (sequenceChars > 0) localStartVertex False postOrderGraph
        checkUnusedEdgesPruneInfty inGS inData pruneEdges warnPruneEdges leafGraph $ updatePhylogeneticGraphCost fullyOptimizedGraph (snd6 fullyOptimizedGraph)

-- | multiTraverseFullyLabelTree performs potorder on default root and other traversal foci, taking the minimum
-- traversal cost for all nonexact charcters--the initial rooting is used for exact characters
-- operates with Tree functions
-- need to add forest functionality--in principle just split into components and optimize them independently
-- but get into root index issues the way this is written now.
-- can either find root, be given root, or start somwhere else (startVertex) do optimize only a component of a forest
multiTraverseFullyLabelTree :: GlobalSettings -> ProcessedData -> DecoratedGraph -> Maybe Int -> SimpleGraph -> PhyG PhylogeneticGraph
multiTraverseFullyLabelTree inGS inData leafGraph startVertex inSimpleGraph =
    if LG.isEmpty inSimpleGraph then pure emptyPhylogeneticGraph
    else
        let sequenceChars = U.getNumberSequenceCharacters (thd3 inData)
            -- False for staticIA
            (postOrderGraph, localStartVertex) = generalizedGraphPostOrderTraversal inGS sequenceChars inData leafGraph False startVertex inSimpleGraph
        in do
            PRE.preOrderTreeTraversal inGS (finalAssignment inGS) False True (sequenceChars > 0) localStartVertex False postOrderGraph
            


-- | generalizedGraphPostOrderTraversal performs the postorder pass
-- on a graph (tree, softWired, or hardWired) to determine the "preliminary" character states
-- include penalty factor cost but not root cost which may or may not be wanted depending on context
-- if full graph--yes, if a component yes or no.
-- hence returnde das pair
generalizedGraphPostOrderTraversal :: GlobalSettings -> Int -> ProcessedData -> DecoratedGraph -> Bool -> Maybe Int -> SimpleGraph -> (PhylogeneticGraph, Int)
generalizedGraphPostOrderTraversal inGS sequenceChars inData leafGraph staticIA startVertex inSimpleGraph =

    -- select postOrder function based on graph type
    let postOrderFunction = if (graphType inGS) == Tree then POSW.postOrderTreeTraversal
                         else if (graphType inGS) == SoftWired then POSW.postOrderSoftWiredTraversal
                         else if (graphType inGS) == HardWired then POSW.postOrderTreeTraversal
                         else error ("Graph type not implemented: " <> (show $ graphType inGS))

        -- first traversal on outgroup root
        outgroupRooted = postOrderFunction inGS inData leafGraph staticIA startVertex inSimpleGraph

        -- start at start vertex--for components or ur-root for full graph
        startVertexList = if startVertex == Nothing then fmap fst $ LG.getRoots $ thd6 outgroupRooted
                          else [fromJust startVertex]

        -- next edges (to vertex in list) to perform rerroting
        -- progresses recursivey over adjacent edges to minimize node reoptimization
        -- childrenOfRoot = concatMap (LG.descendants (thd6 outgroupRooted)) startVertexList
        -- grandChildrenOfRoot = concatMap (LG.descendants (thd6 outgroupRooted)) childrenOfRoot

        -- create list of multi-traversals with original rooting first
        -- subsequent rerooting do not reoptimize exact characters (add nonadd etc)
        -- they are taken from the fully labelled first input decorated graph later when output graph created
        -- it is important that the first graph be the ourgroup rooted graph (outgroupRootedPhyloGraph) so this
        -- will have the preorder assignments for the outgroup rooted graph as 3rd field.  This can be used for incremental
        -- optimization to get O(log n) initial postorder assingment when mutsating graph.
        -- hardwired reroot cause much pain
        -- the head startvertex list for reoptimizing spit trees ni swapping
        recursiveRerootList = if (graphType inGS == HardWired) then [outgroupRooted]
                              else if (graphType inGS == SoftWired) then [POSW.getDisplayBasedRerootSoftWired inGS SoftWired (head startVertexList) outgroupRooted]
                              else if (graphType inGS == Tree) then [POSW.getDisplayBasedRerootSoftWired inGS Tree (head startVertexList) outgroupRooted]
                              else error ("Graph type not implemented: " <> (show $ graphType inGS))


        -- remove if tree reroot code pans out
        finalizedPostOrderGraphList = L.sortOn snd6 recursiveRerootList

        -- create optimal final graph with best costs and best traversal (rerooting) forest for each character
        -- traversal for exact characters (and costs) are the first of each least since exact only optimizaed for that
        -- traversal graph.  The result has approprotate post-order assignments for traversals, preorder "final" assignments
        -- are propagated to the Decorated graph field after the preorder pass.
        -- doesn't have to be sorted, but should minimize assignments
        graphWithBestAssignments = head recursiveRerootList -- L.foldl1' setBetterGraphAssignment recursiveRerootList'

        {-  root and model complexities moved to output
        -- same root cost if same data and number of roots
        localRootCost = rootComplexity inGS
                        {-if (rootCost inGS) == NoRootCost then 0.0
                        else if (rootCost inGS) == Wheeler2015Root then getW15RootCost inGS outgroupRooted
                        else error ("Root cost type " <> (show $ rootCost inGS) <> " is not yet implemented")
                        -}
        -}

    in
    -- trace ("GPOT length: " <> (show $ fmap snd6 recursiveRerootList) <> " " <> (show $ graphType inGS)) (
    -- trace ("TRAV:" <> (show startVertex) <> " " <> (show sequenceChars) <> " " <> (show (snd6 outgroupRooted, fmap snd6 finalizedPostOrderGraphList, snd6 graphWithBestAssignments))
    --    <> "\nTraversal root costs: " <> (show (getTraversalCosts outgroupRooted, fmap getTraversalCosts recursiveRerootList', getTraversalCosts graphWithBestAssignments))) (

    -- only static characters
    if sequenceChars == 0 then
        let penaltyFactor  = if (graphType inGS == Tree) then 0.0

                             --it is its own penalty due to counting all changes in in2 out 1 nodes
                             -- else if (graphType inGS == HardWired) then 0.0

                             -- softwired versions
                             else if (graphFactor inGS) == NoNetworkPenalty then 0.0
                             else if (graphFactor inGS) == Wheeler2015Network then POSW.getW15NetPenaltyFull Nothing inGS inData startVertex outgroupRooted
                             else if (graphFactor inGS) == Wheeler2023Network then POSW.getW23NetPenalty outgroupRooted

                             else error ("Network penalty type " <> (show $ graphFactor inGS) <> " is not yet implemented")

            staticOnlyGraph = if (graphType inGS) == SoftWired then POSW.updateAndFinalizePostOrderSoftWired startVertex (head startVertexList) outgroupRooted
                              else outgroupRooted
            -- staticOnlyGraph = head recursiveRerootList'
            staticOnlyGraph' = if startVertex == Nothing then updatePhylogeneticGraphCost  staticOnlyGraph (penaltyFactor + (snd6 staticOnlyGraph))
                               else updatePhylogeneticGraphCost staticOnlyGraph (penaltyFactor + (snd6 staticOnlyGraph))
        in
        -- trace ("Only static: " <> (snd6 staticOnlyGraph'))
        (staticOnlyGraph', head startVertexList)

    -- single sequence (prealigned, dynamic) only (ie no static)
    else if sequenceChars == 1 && (U.getNumberExactCharacters (thd3 inData) == 0) then
        let penaltyFactorList  = if (graphType inGS == Tree) then replicate (length finalizedPostOrderGraphList) 0.0
                                 -- else if (graphType inGS == HardWired) then replicate (length finalizedPostOrderGraphList) 0.0
                                 else if (graphFactor inGS) == NoNetworkPenalty then replicate (length finalizedPostOrderGraphList) 0.0
                                 else if (graphFactor inGS) == Wheeler2015Network then fmap (POSW.getW15NetPenaltyFull Nothing inGS inData  startVertex) finalizedPostOrderGraphList
                                 else if (graphFactor inGS) == Wheeler2023Network then fmap POSW.getW23NetPenalty finalizedPostOrderGraphList
                                 else error ("Network penalty type " <> (show $ graphFactor inGS) <> " is not yet implemented")
            newCostList = zipWith (+) penaltyFactorList (fmap snd6 finalizedPostOrderGraphList)


            finalizedPostOrderGraph = head $ L.sortOn snd6 $ zipWith updatePhylogeneticGraphCost finalizedPostOrderGraphList newCostList
        in
        -- trace ("GPOT-1: " <> (show (snd6 finalizedPostOrderGraph)))
        (finalizedPostOrderGraph, head startVertexList)

    -- multiple dynamic characters--checks for best root for each character
    -- important to have outgroup rooted graph first for fold so don't use sorted recursive list
    else
        let penaltyFactor  = if (graphType inGS == Tree) then 0.0
                             -- else if (graphType inGS == HardWired) then 0.0
                             else if (graphFactor inGS) == NoNetworkPenalty then 0.0
                             else if (graphFactor inGS) == Wheeler2015Network then POSW.getW15NetPenaltyFull Nothing inGS inData  startVertex graphWithBestAssignments
                             else if (graphFactor inGS) == Wheeler2023Network then POSW.getW23NetPenalty graphWithBestAssignments
                             else error ("Network penalty type " <> (show $ graphFactor inGS) <> " is not yet implemented")

            graphWithBestAssignments' = updatePhylogeneticGraphCost graphWithBestAssignments (penaltyFactor + (snd6 graphWithBestAssignments))
        in
        -- trace ("GPOT-2: " <> (show (penaltyFactor + (snd6 graphWithBestAssignments))))
        (graphWithBestAssignments', head startVertexList)

    -- )


-- | checkUnusedEdgesPruneInfty checks if a softwired phylogenetic graph has
-- "unused" edges sensu Wheeler 2015--that an edge in the canonical graph is
-- not present in any of the block display trees (that are heurstically optimal)
-- the options specify if the cost returned is Infinity (really max bound Double)
-- with no pruning of edges or the cost is left unchanged and unused edges are
-- pruned from the canonical graph
-- this is unDirected due to rerooting heuristic in post/preorder optimization
-- inifinity defined in Types.hs
checkUnusedEdgesPruneInfty :: GlobalSettings -> ProcessedData -> Bool -> Bool -> DecoratedGraph-> PhylogeneticGraph -> PhyG PhylogeneticGraph
checkUnusedEdgesPruneInfty inGS inData pruneEdges warnPruneEdges leafGraph inGraph@(inSimple, _, inCanonical, blockTreeV, charTreeVV, charInfoVV) =
    let simpleEdgeList = LG.edges inSimple
        displayEdgeSet =  L.nubBy LG.undirectedEdgeEquality $ concat $ concat $ fmap (fmap LG.edges) blockTreeV
        unusedEdges = LG.undirectedEdgeMinus simpleEdgeList displayEdgeSet
    in
    -- no unused edges all OK
    if null unusedEdges then pure inGraph

    -- unused edges--do not prune return "infinite cost"
    else if not pruneEdges then
        -- trace ("Unused edge->Infinity")
        pure (inSimple, infinity, inCanonical, blockTreeV, charTreeVV, charInfoVV)

    -- unused but pruned--need to prune nodes and reoptimize to get final assignments correct
    else
        let newSimpleGraph = LG.delEdges unusedEdges inSimple
            contractedSimple = GO.contractIn1Out1EdgesRename newSimpleGraph
        in
        if warnPruneEdges then do
            -- too lazy to thread PhyG logging throuhg everything
            logWith LogWarn ("Pruning " <> (show $ length unusedEdges) <> " unused edges and reoptimizing graph")
            multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph Nothing contractedSimple

        else multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph Nothing contractedSimple

-- | updateGraphCostsComplexities adds root and model complexities if appropriate to graphs
-- updates NCM with roig data due to weights of bitpacking
updateGraphCostsComplexities :: GlobalSettings -> ProcessedData -> ProcessedData -> Bool -> [ReducedPhylogeneticGraph] -> PhyG [ReducedPhylogeneticGraph]
updateGraphCostsComplexities inGS reportingData processedData rediagnoseWithReportingData inGraphList =
    --parallel setup
    let traverseAction :: SimpleGraph -> PhyG ReducedPhylogeneticGraph
        traverseAction = multiTraverseFullyLabelGraphReduced inGS reportingData False False Nothing
    in

    if optimalityCriterion inGS == Parsimony then do
        pure inGraphList

    else if optimalityCriterion inGS `elem` [SI, MAPA] then do
        logWith LogInfo ("\tFinalizing graph cost with root priors" <> "\n")
        pure $ updatePhylogeneticGraphCostList (rootComplexity inGS) inGraphList

    else if optimalityCriterion inGS `elem` [NCM] then do
        updatedGraphList <- if (reportingData == emptyProcessedData) || (not rediagnoseWithReportingData) || (not $ U.has4864PackedChars (thd3 processedData)) then
                                 -- trace ("\t\tCannot update cost with original data--skipping")
                                 pure $ updatePhylogeneticGraphCostList (rootComplexity inGS) inGraphList
                            else do
                                --let newGraphList = PU.seqParMap PU.myStrategy  (multiTraverseFullyLabelGraphReduced inGS reportingData False False Nothing) (fmap fst5 inGraphList)
                                --in
                                traversePar <- getParallelChunkTraverse
                                newGraphList <- traversePar traverseAction (fmap fst5 inGraphList)
                                pure $ updatePhylogeneticGraphCostList (rootComplexity inGS) newGraphList
        

        logWith LogInfo  ("\tFinalizing graph cost (updating NCM) with root priors" <> "\n") 
        pure updatedGraphList

    else if optimalityCriterion inGS == PMDL then
        -- trace ("\tFinalizing graph cost with model and root complexities")
        pure $ updatePhylogeneticGraphCostList ((rootComplexity inGS) + (modelComplexity inGS)) inGraphList

    else error ("Optimality criterion not recognized/implemented: " <> (show $ optimalityCriterion inGS))


-- | updatePhylogeneticGraphCostList is a list wrapper for updatePhylogeneticGraphCost
updatePhylogeneticGraphCostList :: VertexCost -> [ReducedPhylogeneticGraph] -> [ReducedPhylogeneticGraph]
updatePhylogeneticGraphCostList rootCost inGraphList =
    fmap (updateCost rootCost) inGraphList
    where updateCost :: forall {b} {a} {c} {d} {e}. Num b => b -> (a, b, c, d, e) -> (a, b, c, d, e)
          updateCost z (a, oldCost, b, c, e) = (a, oldCost + z, b, c, e)

-- | updatePhylogeneticGraphCost takes a PhylgeneticGrtaph and Double and replaces the cost (snd of 6 fields)
-- and returns Phylogenetic graph
updatePhylogeneticGraphCost :: PhylogeneticGraph -> VertexCost -> PhylogeneticGraph
updatePhylogeneticGraphCost (a, _, b, c, d, e) newCost = (a, newCost, b, c, d, e)

-- | updatePhylogeneticGraphCost takes a ReducedPhylogeneticGraph and Double and replaces the cost (snd of 6 fields)
-- and returns Phylogenetic graph
updatePhylogeneticGraphCostReduced :: ReducedPhylogeneticGraph -> VertexCost -> ReducedPhylogeneticGraph
updatePhylogeneticGraphCostReduced (a, _, b, c, e) newCost = (a, newCost, b, c, e)

