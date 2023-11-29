{- |
Module specifying charcter-based Wagner tree building functions.
-}
module Search.WagnerBuild (
    wagnerTreeBuild,
    wagnerTreeBuild',
    rasWagnerBuild,
) where

import Control.Monad (when)
import Control.Monad.IO.Class (MonadIO (..))
import Data.Maybe
import Data.Text.Lazy qualified as TL
import Data.Vector qualified as V
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
import Utilities.Utilities qualified as U


-- import Debug.Trace

-- Haven't added in unions--but prob should

{- | rasWagnerBuild generates a series of random addition sequences and then calls wagnerTreeBuild to construct them.
Does not filter by best, unique etc.  That happens with the select() command specified separately.
-}
rasWagnerBuild ∷ GlobalSettings → ProcessedData → Int → Int → PhyG [ReducedPhylogeneticGraph]
rasWagnerBuild inGS inData rSeed numReplicates =
    if numReplicates == 0
        then do
            pure []
        else
            let numLeaves = V.length $ fst3 inData
                randomizedAdditionSequences = V.fromList <$> shuffleInt rSeed numReplicates [0 .. numLeaves - 1]

                -- "graph" of leaf nodes without any edges
                leafGraph = GO.makeSimpleLeafGraph inData
                leafDecGraph = GO.makeLeafGraph inData

                hasNonExactChars = U.getNumberSequenceCharacters (thd3 inData) > 0

                wagnerTreeAction ∷ (V.Vector Int, Int) → PhyG ReducedPhylogeneticGraph
                wagnerTreeAction = wagnerTreeBuild' inGS inData leafGraph leafDecGraph numLeaves hasNonExactChars
            in  do
                    logWith LogInfo ("\t\tBuilding " <> show numReplicates <> " character Wagner replicates" <> "\n")
                    -- seqParMap better for high level parallel stuff
                    -- PU.seqParMap PU.myStrategy (wagnerTreeBuild inGS inData) randomizedAdditionSequences
                    -- zipWith (wagnerTreeBuild inGS inData leafGraph leafDecGraph numLeaves hasNonExactChars) randomizedAdditionSequences [0..numReplicates - 1] `using` PU.myParListChunkRDS

                    -- TODO
                    wagnerTreeActionTraverse ← getParallelChunkTraverse
                    rasResult ← wagnerTreeActionTraverse wagnerTreeAction (zip randomizedAdditionSequences [0 .. numReplicates - 1])

                    -- rasResult <- mapM (wagnerTreeBuild' inGS inData leafGraph leafDecGraph numLeaves hasNonExactChars) (zip randomizedAdditionSequences [0..numReplicates - 1])
                    -- PU.seqParMap (parStrategy $ strictParStrat inGS) (wagnerTreeBuild' inGS inData leafGraph leafDecGraph numLeaves hasNonExactChars) (zip randomizedAdditionSequences [0..numReplicates - 1])
                    -- fmap (wagnerTreeBuild' inGS inData leafGraph leafDecGraph numLeaves hasNonExactChars) (zip randomizedAdditionSequences [0..numReplicates - 1]) `using` PU.myParListChunkRDS
                    pure rasResult


-- | wagnerTreeBuild' is a wrapper around wagnerTreeBuild to allow for better parallelization--(zipWith not doing so well?)
wagnerTreeBuild'
    ∷ GlobalSettings
    → ProcessedData
    → SimpleGraph
    → DecoratedGraph
    → Int
    → Bool
    → (V.Vector Int, Int)
    → PhyG ReducedPhylogeneticGraph
wagnerTreeBuild' inGS inData leafSimpleGraph leafDecGraph numLeaves hasNonExactChars (additionSequence, replicateIndex) = do
    wagResult ←
        wagnerTreeBuild inGS inData leafSimpleGraph leafDecGraph numLeaves hasNonExactChars additionSequence replicateIndex
    pure $ GO.convertPhylogeneticGraph2Reduced $ wagResult


{- | wagnerTreeBuild builds a wagner tree (Farris 1970--but using random addition seqeuces--not "best" addition)
from a leaf addition sequence. Always produces a tree that can be converted to a soft/hard wired network
afterwards
basic procs is to add edges to unresolved tree
currently naive wrt candidate tree costs
-}
wagnerTreeBuild
    ∷ GlobalSettings → ProcessedData → SimpleGraph → DecoratedGraph → Int → Bool → V.Vector Int → Int → PhyG PhylogeneticGraph
wagnerTreeBuild inGS inData leafSimpleGraph leafDecGraph numLeaves hasNonExactChars additionSequence replicateIndex =
    let rootHTU = (numLeaves, TL.pack $ "HTU" <> show numLeaves)
        nextHTU = (numLeaves + 1, TL.pack $ "HTU" <> show (numLeaves + 1))

        edge0 = (numLeaves, additionSequence V.! 0, 0.0)
        edge1 = (numLeaves, numLeaves + 1, 0.0)
        edge2 = (numLeaves + 1, additionSequence V.! 1, 0.0)
        edge3 = (numLeaves + 1, additionSequence V.! 2, 0.0)

        initialTree = LG.insEdges [edge0, edge1, edge2, edge3] $ LG.insNodes [rootHTU, nextHTU] leafSimpleGraph

        blockCharInfo = V.map thd3 $ thd3 inData

        -- initialFullyDecoratedTree = T.multiTraverseFullyLabelTree inGS inData initialTree
        -- False flag for staticIA--can't be done in build
        calculateBranchLengths = False -- must be True for delata using existing edge
        initialPostOrderTree = POSW.postDecorateTree inGS False initialTree leafDecGraph blockCharInfo numLeaves numLeaves
    in  do
            initialFullyDecoratedTree ←
                PRE.preOrderTreeTraversal
                    inGS
                    (finalAssignment inGS)
                    False
                    calculateBranchLengths
                    hasNonExactChars
                    numLeaves
                    False
                    initialPostOrderTree
            -- this used for missing data adjustments during build
            maxDistance ← U.getMaxNumberObservations (thd3 inData)

            logWith LogInfo ("\tBuilding Wagner replicate " <> show replicateIndex <> "\n")
            wagnerTree ←
                recursiveAddEdgesWagner
                    maxDistance
                    (useIA inGS)
                    (V.drop 3 additionSequence)
                    numLeaves
                    (numLeaves + 2)
                    inGS
                    inData
                    hasNonExactChars
                    leafDecGraph
                    initialFullyDecoratedTree
            pure wagnerTree


{- | recursiveAddEdgesWagner adds edges until 2n -1 (n leaves) vertices in graph
this tested by null additin sequence list
interface will change with correct final states--using post-order pass for now
-}
recursiveAddEdgesWagner
    ∷ VertexCost
    → Bool
    → V.Vector Int
    → Int
    → Int
    → GlobalSettings
    → ProcessedData
    → Bool
    → DecoratedGraph
    → PhylogeneticGraph
    → PhyG PhylogeneticGraph
recursiveAddEdgesWagner maxDistance useIA additionSequence numLeaves numVerts inGS inData hasNonExactChars leafDecGraph inGraph@(inSimple, _, inDecGraph, _, _, charInfoVV) =
    -- all edges/ taxa in graph
    -- trace ("To go " <> (show additionSequence) <> " verts " <> (show numVerts)) (
    if null additionSequence
        then pure inGraph
        else -- trace ("RAEW-In: " <> (show $ length additionSequence)) (
        -- edges/taxa to add, but not the edges that leads to outgroup--redundant with its sister edge

            let -- outgroupEdges = filter ((< numLeaves) . snd3) $ LG.out inDecGraph numLeaves
                edgesToInvade = LG.labEdges inDecGraph -- L.\\ outgroupEdges
                leafToAdd = V.head additionSequence
                leafToAddVertData = vertData $ fromJust $ LG.lab inDecGraph leafToAdd

                -- since this is apporximate--can get a bit off
                {-
                --add unions here
                -- not clear what the delta to compare is--have an existing tree cost and the leaf to add would be 0.
                -- use a single tree after first addition?
                unionEdgeList = S.getUnionRejoinEdgeList inGS inDecGraph charInfoVV [numLeaves] splitDeltaValue (unionThreshold inGS) leafToAddVertData []
                -}

                addTaxonAction ∷ LG.LEdge EdgeInfo → (VertexCost, LG.LNode TL.Text, [LG.LEdge Double], LG.Edge)
                addTaxonAction = addTaxonWagner maxDistance useIA numVerts inGraph leafToAddVertData leafToAdd
            in  do
                    --- TODO
                    addTaxonWagnerPar ← getParallelChunkMap
                    let candidateEditList = addTaxonWagnerPar addTaxonAction edgesToInvade
                    -- let candidateEditList = PU.seqParMap  (parStrategy $ lazyParStrat inGS)  (addTaxonWagner maxDistance useIA numVerts inGraph leafToAddVertData leafToAdd) edgesToInvade

                    let minDelta = minimum $ fmap fst4 candidateEditList
                    let (_, nodeToAdd, edgesToAdd, edgeToDelete) = head $ filter ((== minDelta) . fst4) candidateEditList

                    -- create new tree
                    let newSimple = LG.insEdges edgesToAdd $ LG.insNode nodeToAdd $ LG.delEdge edgeToDelete inSimple

                    -- this reroot since could add taxon sister to outgroup
                    let newSimple' = LG.rerootTree (outgroupIndex inGS) newSimple

                    -- create fully labelled tree, if all taxa in do full multi-labelled for correct graph type
                    -- False flag for static IA--can't do when adding in new leaves
                    let calculateBranchLengths = False -- must be True for delata using existing edge
                    newPhyloGraph ← -- T.multiTraverseFullyLabelTree inGS inData leafDecGraph (Just numLeaves) newSimple'
                        if V.length additionSequence > 1
                            then
                                PRE.preOrderTreeTraversal inGS (finalAssignment inGS) False calculateBranchLengths hasNonExactChars numLeaves False $
                                    POSW.postDecorateTree inGS False newSimple' leafDecGraph charInfoVV numLeaves numLeaves
                            else T.multiTraverseFullyLabelTree inGS inData leafDecGraph (Just numLeaves) newSimple'

                    if isNothing (LG.lab inDecGraph leafToAdd)
                        then error "Missing label data for vertices"
                        else
                            recursiveAddEdgesWagner
                                maxDistance
                                useIA
                                (V.tail additionSequence)
                                numLeaves
                                (numVerts + 1)
                                inGS
                                inData
                                hasNonExactChars
                                leafDecGraph
                                newPhyloGraph


-- )

{- | addTaxonWagner adds a taxon (really edges) by 'invading' and edge, deleting that adege and creteing 3 more
to existing tree and gets cost (for now by postorder traversal--so wasteful but will be by final states later)
returns a tuple of the cost, node to add, edges to add, edge to delete
-}
addTaxonWagner
    ∷ VertexCost
    → Bool
    → Int
    → PhylogeneticGraph
    → VertexBlockData
    → Int
    → LG.LEdge EdgeInfo
    → (VertexCost, LG.LNode TL.Text, [LG.LEdge Double], LG.Edge)
addTaxonWagner maxDistance useIA numVerts (_, _, inDecGraph, _, _, charInfoVV) leafToAddVertData leafToAdd targetEdge =
    let edge0 = (numVerts, leafToAdd, 0.0)
        edge1 = (fst3 targetEdge, numVerts, 0.0)
        edge2 = (numVerts, snd3 targetEdge, 0.0)
        newNode = (numVerts, TL.pack ("HTU" <> show numVerts))

        -- full post order
        -- newSimpleGraph =  LG.insEdges [edge0, edge1, edge2] $ LG.insNode newNode $ LG.delEdge (LG.toEdge targetEdge) inSimple
        -- newCost = snd6 $ T.postDecorateTree newSimpleGraph leafDecGraph charInfoVV numLeaves

        -- heuristic delta
        (delta, edgeUnionVertData) = getDelta useIA leafToAddVertData targetEdge inDecGraph charInfoVV

        -- modification for missing data
        nonMissingDistance = U.getPairwiseObservationsGraph leafToAddVertData edgeUnionVertData

        deltaNormalized = delta * maxDistance / (max 1.0 nonMissingDistance)
    in  -- trace ("ATW :" <> (show $ maxDistance / (max 1.0 nonMissingDistance))) $
        (deltaNormalized, newNode, [edge0, edge1, edge2], LG.toEdge targetEdge)


-- (newCost, newNode, [edge0, edge1, edge2], LG.toEdge targetEdge)

{- | getDelta estimates the delta in tree cost by adding a leaf taxon in Wagner build
must be DO for this--isolated leaves won't have IA
-}
getDelta
    ∷ Bool → VertexBlockData → LG.LEdge EdgeInfo → DecoratedGraph → V.Vector (V.Vector CharInfo) → (VertexCost, VertexBlockData)
getDelta useIA leafToAddVertData (eNode, vNode, _) inDecGraph charInfoVV =
    let eNodeVertData = vertData $ fromJust $ LG.lab inDecGraph eNode
        vNodeVertData = vertData $ fromJust $ LG.lab inDecGraph vNode

        -- create edge union 'character' blockData
        -- filters gaps (True argument) because using DOm (as must) to add taxa not in IA framework
        -- edge union based on final IA assignments filtering gaps (True True)
        edgeUnionVertData = M.createEdgeUnionOverBlocks useIA True eNodeVertData vNodeVertData charInfoVV []
    in  -- trace ("GD: " <> (show edgeUnionVertData)) (
        if isNothing (LG.lab inDecGraph eNode) || isNothing (LG.lab inDecGraph vNode)
            then error "Missing label data for vertices"
            else
                let -- Use edge union data for delta to edge data
                    dLeafEdgeUnionCost = sum (fst <$> V.zipWith3 (PRE.getBlockCostPairsFinal DirectOptimization) leafToAddVertData edgeUnionVertData charInfoVV)
                in  -- should be able to use existing information--but for now using this
                    -- existingEdgeCost' = sum $ fmap fst $ V.zipWith3 (PRE.getBlockCostPairsFinal DirectOptimization) eNodeVertData vNodeVertData charInfoVV

                    -- trace ("Delta: " <> (show (dLeafENode, dLeafVNode, existingEdgeCost)))
                    -- dLeafENode + dLeafVNode - existingEdgeCost
                    -- trace ("Delta: " <> (show dLeafEdgeUnionCost) <> " vs " <> (show dLeafEVAddCost))

                    -- min dLeafEdgeUnionCost dLeafEVAddCost
                    -- dLeafEVAddCost
                    (dLeafEdgeUnionCost, edgeUnionVertData)

-- )
