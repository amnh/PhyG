{- |
Module specifying charcter-based Wagner tree building functions.
-}
module Search.WagnerBuild (
    wagnerTreeBuild,
    wagnerTreeBuild',
    rasWagnerBuild,
) where

import Control.Monad (replicateM)
import Data.List qualified as L
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
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import Types.Types
import Utilities.LocalGraph qualified as LG
import Utilities.Utilities qualified as U


-- import Debug.Trace

-- Haven't added in unions--but prob should

{- | rasWagnerBuild generates a series of random addition sequences and then calls wagnerTreeBuild to construct them.
Does not filter by best, unique etc.  That happens with the select() command specified separately.
-}
rasWagnerBuild ∷ GlobalSettings → ProcessedData → Int → PhyG [ReducedPhylogeneticGraph]
rasWagnerBuild inGS inData numReplicates =
    if numReplicates == 0
        then do
            pure []
        else
            let numLeaves = V.length $ fst3 inData

                -- "graph" of leaf nodes without any edges
                leafGraph = GO.makeSimpleLeafGraph inData
                leafDecGraph = GO.makeLeafGraph inData
                leafIndexVec = V.generate numLeaves id

                hasNonExactChars = U.getNumberSequenceCharacters (thd3 inData) > 0

                wagnerTreeAction ∷ (V.Vector Int, Int) → PhyG ReducedPhylogeneticGraph
                wagnerTreeAction = wagnerTreeBuild' inGS inData leafGraph leafDecGraph numLeaves hasNonExactChars
            in  do
                    randomizedAdditionSequences ← replicateM numReplicates $ shuffleList leafIndexVec
                    logWith LogInfo ("\t\tBuilding " <> show numReplicates <> " character Wagner replicates" <> "\n")
                    getParallelChunkTraverseBy U.strict2of5 >>= \pTraverse →
                        pTraverse wagnerTreeAction $ zip randomizedAdditionSequences [0 .. numReplicates - 1]


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


{- | wagnerTreeBuild builds a wagner tree (Farris 1970--but using random addition sequeces--not "best" addition)
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
    in  do
            initialPostOrderTree ← POSW.postDecorateTree inGS False initialTree leafDecGraph blockCharInfo numLeaves numLeaves
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

                addTaxonAction ∷ LG.LEdge EdgeInfo → PhyG (VertexCost, LG.LNode TL.Text, [LG.LEdge Double], LG.Edge)
                addTaxonAction = addTaxonWagner maxDistance useIA numVerts inGraph leafToAddVertData leafToAdd

                -- False flag for static IA--can't do when adding in new leaves
                postOrderAction :: SimpleGraph → PhyG PhylogeneticGraph
                postOrderAction = POSW.postDecorateTreeForList inGS False leafDecGraph charInfoVV numLeaves numLeaves
            in  do
                    --- TODO
                    addTaxonWagnerPar ← getParallelChunkTraverseBy U.strict1of4
                    candidateEditList ← addTaxonWagnerPar addTaxonAction edgesToInvade
                    -- let candidateEditList = PU.seqParMap  (parStrategy $ lazyParStrat inGS)  (addTaxonWagner maxDistance useIA numVerts inGraph leafToAddVertData leafToAdd) edgesToInvade

                    let bestNCandidateEdgesList = take (graphsSteepest inGS) $ L.sortOn fst4 candidateEditList
                   
                    let newSimpleRerootedList = fmap (createNewSimpleGraph (outgroupIndex inGS) inSimple) bestNCandidateEdgesList
                    
                    -- create fully labelled tree, if all taxa in do full multi-labelled for correct graph type
                    -- False flag for static IA--can't do when adding in new leaves
                    let calculateBranchLengths = False -- must be True for delaa using existing edge
                    
                    -- parallel postorder check on cost
                    postOrderPar <- getParallelChunkTraverseBy U.strict2of6
                    postOrderCandidateList <- postOrderPar postOrderAction newSimpleRerootedList
                    let bestCandCost = minimum $ fmap snd6 $ postOrderCandidateList
                    let postOrderStuff =  head $ filter ((== bestCandCost) . snd6) postOrderCandidateList

                    -- postOrderStuff ← POSW.postDecorateTree inGS False newSimple' leafDecGraph charInfoVV numLeaves numLeaves


                    newPhyloGraph ← -- T.multiTraverseFullyLabelTree inGS inData leafDecGraph (Just numLeaves) newSimple'
                        if V.length additionSequence > 1
                            then
                                PRE.preOrderTreeTraversal
                                    inGS
                                    (finalAssignment inGS)
                                    False
                                    calculateBranchLengths
                                    hasNonExactChars
                                    numLeaves
                                    False
                                    postOrderStuff
                            else T.multiTraverseFullyLabelTree inGS inData leafDecGraph (Just numLeaves) (fst6 postOrderStuff)

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


{- | createNewSimpleGraph perfomrs the addition, new edges, rerotting tasks for a candidate tree
-}
createNewSimpleGraph :: LG.Node -> SimpleGraph -> (VertexCost, LG.LNode TL.Text, [LG.LEdge Double], LG.Edge) -> SimpleGraph
createNewSimpleGraph outgroupIndex inSimple (_, nodeToAdd, edgesToAdd, edgeToDelete) =
    -- create new tree
    let newSimple = LG.insEdges edgesToAdd $ LG.insNode nodeToAdd $ LG.delEdge edgeToDelete inSimple

    -- this reroot since could add taxon sister to outgroup
    in
    LG.rerootTree outgroupIndex newSimple

{- | addTaxonWagner adds a taxon (really edges) by 'invading' and edge, deleting that adege and creating 3 more
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
    → PhyG (VertexCost, LG.LNode TL.Text, [LG.LEdge Double], LG.Edge)
addTaxonWagner maxDistance useIA numVerts (_, _, inDecGraph, _, _, charInfoVV) leafToAddVertData leafToAdd targetEdge =
    let edge0 = (numVerts, leafToAdd, 0.0)
        edge1 = (fst3 targetEdge, numVerts, 0.0)
        edge2 = (numVerts, snd3 targetEdge, 0.0)
        newNode = (numVerts, TL.pack ("HTU" <> show numVerts))
    in  -- full post order
        do
            -- heuristic delta
            (delta, edgeUnionVertData) ← getDelta useIA leafToAddVertData targetEdge inDecGraph charInfoVV

            -- modification for missing data
            let nonMissingDistance = U.getPairwiseObservationsGraph leafToAddVertData edgeUnionVertData

            let deltaNormalized = delta * maxDistance / (max 1.0 nonMissingDistance)

            pure (deltaNormalized, newNode, [edge0, edge1, edge2], LG.toEdge targetEdge)


-- (newCost, newNode, [edge0, edge1, edge2], LG.toEdge targetEdge)

{- | getDelta estimates the delta in tree cost by adding a leaf taxon in Wagner build
must be DO for this--isolated leaves won't have IA
-}
getDelta
    ∷ Bool → VertexBlockData → LG.LEdge EdgeInfo → DecoratedGraph → V.Vector (V.Vector CharInfo) → PhyG (VertexCost, VertexBlockData)
getDelta useIA leafToAddVertData (eNode, vNode, _) inDecGraph charInfoVV =
    let eNodeVertData = vertData $ fromJust $ LG.lab inDecGraph eNode
        vNodeVertData = vertData $ fromJust $ LG.lab inDecGraph vNode
    in  -- create edge union 'character' blockData
        -- filters gaps (True argument) because using DOm (as must) to add taxa not in IA framework
        -- edge union based on final IA assignments filtering gaps (True True)

        if isNothing (LG.lab inDecGraph eNode) || isNothing (LG.lab inDecGraph vNode)
            then error "Missing label data for vertices"
            else do
                edgeUnionVertData ← M.createEdgeUnionOverBlocksM useIA True eNodeVertData vNodeVertData charInfoVV
                -- Use edge union data for delta to edge data

                let dLeafEdgeUnionCost = sum (fst <$> V.zipWith3 (PRE.getBlockCostPairsFinal DirectOptimization) leafToAddVertData edgeUnionVertData charInfoVV)

                pure (dLeafEdgeUnionCost, edgeUnionVertData)
