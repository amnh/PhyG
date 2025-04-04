{- |
Module specifying graph egde adding and deleting functions
-}
module Search.NetworkAddDeleteV2 (
    --deleteAllNetEdges,
    insertAllNetEdges,
    --moveAllNetEdges,
    deltaPenaltyAdjustment,
    --deleteNetEdge,
    --deleteOneNetAddAll,
    --addDeleteNetEdges,
    getCharacterDelta,
    getBlockDelta,
    -- these are not used but to quiet warnings
    --heuristicDeleteDelta,
    heuristicAddDelta,
    heuristicAddDelta',
) where

import Control.Arrow ((&&&))
import Control.Monad (when)
import Control.Monad.IO.Class (MonadIO (..))
import Control.Monad.Random.Class
import Data.Bits
import Data.Foldable (fold)
import Data.Functor ((<&>))
import Data.InfList qualified as IL
import Data.List qualified as L
import Data.Maybe
import Data.Text.Lazy qualified as TL
import Data.Vector qualified as V
import GeneralUtilities
import GraphOptimization.Medians qualified as M
import GraphOptimization.PostOrderSoftWiredFunctions qualified as POSW
import GraphOptimization.PostOrderSoftWiredFunctionsNew qualified as NEW
import GraphOptimization.PreOrderFunctions qualified as PRE
import GraphOptimization.Traversals qualified as T
import Graphs.GraphOperations qualified as GO
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging
import PHANE.Evaluation.Verbosity (Verbosity (..))
import Types.Types
import Utilities.LocalGraph qualified as LG
import Utilities.Utilities qualified as U



-- | (curBestGraphList, annealBestCost) is a wrapper for moveAllNetEdges' allowing for multiple simulated annealing rounds
insertAllNetEdges
    ∷ GlobalSettings
    → ProcessedData
    → NetParams
    → Int
    → Int
    → ([ReducedPhylogeneticGraph], VertexCost)
    → (Maybe SAParams, [ReducedPhylogeneticGraph])
    → PhyG ([ReducedPhylogeneticGraph], Int)
insertAllNetEdges inGS inData netParams maxRounds counter (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
    
    insertAllNetEdges'
        inGS
        inData
        netParams
        counter
        (curBestGraphList, curBestGraphCost)
        inSimAnnealParams
        inPhyloGraphList
            


{- | insertAllNetEdges' adds network edges one each each round until no better or additional
graphs are found
call with ([], infinity) [single input graph]
-}
insertAllNetEdges'
    ∷ GlobalSettings
    → ProcessedData
    → NetParams
    → Int
    → ([ReducedPhylogeneticGraph], VertexCost)
    → Maybe SAParams
    → [ReducedPhylogeneticGraph]
    → PhyG ([ReducedPhylogeneticGraph], Int)
insertAllNetEdges' inGS inData netParams counter (curBestGraphList, curBestGraphCost) inSimAnnealParams inPhyloGraphList = 
    if null inPhyloGraphList then 
        pure (take (netKeepNum netParams) curBestGraphList, counter)

    else 
        let firstPhyloGraph : otherPhyloGraphs = inPhyloGraphList

            currentCost = min curBestGraphCost $ snd5 firstPhyloGraph

            -- check for max net edges
            netNodes = fth4 $ LG.splitVertexList $ thd5 firstPhyloGraph
        in  do
                logWith LogInfo ("\n\tNumber of network edges: " <> (show $ length netNodes) <> " Number of graphs: " <> (show $ length inPhyloGraphList) <> "\n")

                if length netNodes >= (netMaxEdges netParams) then do
                    logWith LogInfo $ unwords ["Maximum number of network edges reached:", show $ length netNodes, "\n"]
                    pure (take (netKeepNum netParams) curBestGraphList, counter)

                else do
                    -- this examines all possible net adds via heursitc first then verifies
                    (newGraphList, _, newSAParams) ←
                        insertEachNetEdgeHeuristicGather 
                            inGS
                            inData
                            netParams
                            Nothing
                            inSimAnnealParams
                            firstPhyloGraph

                    bestNewGraphList ← GO.selectGraphs Best (outgroupIndex inGS) (netKeepNum netParams) 0 newGraphList
                    let newGraphCost = case bestNewGraphList of
                            [] → infinity
                            (_, c, _, _, _) : _ → c

                    
                    postProcessNetworkAdd
                            inGS
                            inData
                            netParams
                            counter
                            (curBestGraphList, curBestGraphCost)
                            (bestNewGraphList, newGraphCost)
                            inSimAnnealParams
                            currentCost
                            otherPhyloGraphs
                    



-- | postProcessNetworkAdd prcesses non-simanneal/drift--so no updating of SAParams
postProcessNetworkAdd
    ∷ GlobalSettings
    → ProcessedData
    -> NetParams
    → Int
    → ([ReducedPhylogeneticGraph], VertexCost)
    → ([ReducedPhylogeneticGraph], VertexCost)
    → Maybe SAParams
    → VertexCost
    → [ReducedPhylogeneticGraph]
    → PhyG ([ReducedPhylogeneticGraph], Int)
postProcessNetworkAdd inGS inData netParams counter (curBestGraphList, _) (newGraphList, newGraphCost) inSimAnnealParams currentCost inPhyloGraphList = case newGraphCost `compare` currentCost of
    -- "steepest style descent" abandons existing list if better cost found

    -- check if graph OK--done in insert function
    LT →
        let graphsToInsert
                | (netSteepest netParams) = newGraphList
                | otherwise = take (netKeepNum netParams) $ newGraphList <> inPhyloGraphList
        in  do
                -- prob do this if steepest, "all" does in middle for more frequent output
                -- logWith LogInfo ("\t-> " <> show newGraphCost)
                insertAllNetEdges'
                    inGS
                    inData
                    netParams
                    (counter + 1)
                    (newGraphList, newGraphCost)
                    inSimAnnealParams
                    graphsToInsert

    -- worse graphs found--go on
    GT →
        insertAllNetEdges'
            inGS
            inData
            netParams
            (counter + 1)
            (curBestGraphList, currentCost)
            inSimAnnealParams
            inPhyloGraphList
    -- equal cost
    EQ → do
        -- new graph list contains the input graph if equal and filterd unique already in insertAllNetEdges
        newCurSameBestList ← GO.selectGraphs Unique (outgroupIndex inGS) (netKeepNum netParams) 0 $ curBestGraphList <> newGraphList
        insertAllNetEdges'
            inGS
            inData
            netParams
            (counter + 1)
            (newCurSameBestList, currentCost)
            inSimAnnealParams
            inPhyloGraphList

{- | insertEachNetEdgeHeuristicGather takes a phylogenetic graph and inserts all permissible network edges all
in parallel and returns unique list of new Phylogenetic Graphs and cost
if equal returns unique graph list

Operates by first perfomring heuriastic cost estimate for all possibilities,
then fully evaluates the best (ie. lowest delta number) graphs based on
HeursticCheck option.
-}
insertEachNetEdgeHeuristicGather
    ∷ GlobalSettings
    → ProcessedData
    -> NetParams
    → Maybe VertexCost
    → Maybe SAParams
    → ReducedPhylogeneticGraph
    → PhyG ([ReducedPhylogeneticGraph], VertexCost, Maybe SAParams)
insertEachNetEdgeHeuristicGather inGS inData netParams preDeleteCost inSimAnnealParams inPhyloGraph =
    if LG.isEmpty $ fst5 inPhyloGraph
        then error "Empty input insertEachNetEdge graph in deleteAllNetEdges"
        else
            let currentCost =
                    if isNothing preDeleteCost
                        then snd5 inPhyloGraph
                        else fromJust preDeleteCost

                (_, _, _, netNodes) = LG.splitVertexList (thd5 inPhyloGraph)

                -- parallel stuff
                heuristicAction ∷ (LG.LEdge b, LG.LEdge b) → PhyG (VertexCost, SimpleGraph)
                heuristicAction = insertNetEdgeHeuristic inGS inData inPhyloGraph

                -- this needs to start with root hence Nothing
                diagnoseAction :: SimpleGraph → PhyG ReducedPhylogeneticGraph
                diagnoseAction = T.multiTraverseFullyLabelGraphReduced inGS inData False False Nothing

            in  do
                    if (length netNodes >= (netMaxEdges netParams)) then do
                            logWith LogInfo ("Maximum number of network edges reached: " <> (show $ length netNodes) <> "\n")
                            pure ([inPhyloGraph], currentCost, inSimAnnealParams)
                    else do
                        candidateNetworkEdgeList' ← getPermissibleEdgePairs inGS (thd5 inPhyloGraph)

                        -- radomize pair list--not needed if doing all (as in here).
                        candidateNetworkEdgeList ←
                            if (netRandom netParams)
                                then shuffleList candidateNetworkEdgeList'
                                else pure candidateNetworkEdgeList'

                        logWith LogInfo ("\t\tExamining at most " <> (show $ length candidateNetworkEdgeList) <> " candidate edge pairs" <> "\n")

                        -- get heuristic costs and simple graphs
                        heurCostSimpleGraphPairList <- 
                            getParallelChunkTraverse >>= \pTraverse →
                                heuristicAction `pTraverse` candidateNetworkEdgeList

                        -- filter out non-phylo graphs and sort on heuristic cost
                        let candidatePairList = L.sortOn fst $ filter ((/= infinity) .fst) heurCostSimpleGraphPairList

                        let graphsToBeEvaluated = 
                                    if (netCheckHeuristic netParams) == BestOnly then
                                        fmap snd $ take 1 candidatePairList

                                    else if (netCheckHeuristic netParams) == Better then
                                        fmap snd $ filter ((< 0) . fst) candidatePairList

                                    else if (netCheckHeuristic netParams) == BetterN then 
                                        fmap snd $ take (graphsSteepest inGS) candidatePairList

                                    else --BestAll
                                        fmap snd candidatePairList

                        if null graphsToBeEvaluated then 
                            pure ([inPhyloGraph], currentCost, inSimAnnealParams)

                        else do
                                -- logWith LogInfo ("\nEvaluating " <> (show $ length graphsToBeEvaluated) <> " candidate graphs" <> " of " <> (show $ length candidatePairList) <> "\n")
                                -- rediagnose some fraction of returned simple graphs--lazy in cost so return only thos need nlater
                                diagnoseActionPar <- (getParallelChunkTraverseBy snd5)
                                checkedGraphCosts <- diagnoseActionPar diagnoseAction graphsToBeEvaluated


                                --pure ([inPhyloGraph], currentCost, inSimAnnealParams)


                                
                                (newGraphList, newSAParams) <-
                                        -- do all in prallel
                                       let genNewSimAnnealParams =
                                            if isNothing inSimAnnealParams
                                                then Nothing
                                                else U.incrementSimAnnealParams inSimAnnealParams
                                       in do
                                            pure (checkedGraphCosts, genNewSimAnnealParams)

                                       
                                let minCost =
                                        if null graphsToBeEvaluated || null newGraphList
                                            then infinity
                                            else minimum $ fmap snd5 newGraphList

                                
                                if minCost < (snd5 inPhyloGraph) then 
                                    logWith LogInfo ("\t\t -> " <> (show minCost) <> "\n")
                                else logWith LogInfo ("")
                                
                                
                                -- hit max network edges to insert
                                if (length netNodes >= (netMaxEdges netParams))
                                    then do
                                        logWith LogInfo ("Maximum number of network edges reached: " <> (show $ length netNodes) <> "\n")
                                        if minCost < (snd5 inPhyloGraph) then do
                                            newBestGraphList <- GO.selectGraphs Best (outgroupIndex inGS) (netKeepNum netParams) 0 newGraphList
                                            pure (newBestGraphList, currentCost, inSimAnnealParams)
                                        else do
                                            newBestGraphList <- GO.selectGraphs Best (outgroupIndex inGS) (netKeepNum netParams) 0 (inPhyloGraph : newGraphList)
                                            pure (newBestGraphList, min minCost (snd5 inPhyloGraph), inSimAnnealParams)

                                -- not reached max netedges if found better recurse, if not then return input graph
                                else 
                                    -- retunr if better or equal empty if not--sorted upon return
                                    if minCost <= (snd5 inPhyloGraph) then
                                        let genNewSimAnnealParams =
                                                if isNothing inSimAnnealParams
                                                    then Nothing
                                                    else U.incrementSimAnnealParams inSimAnnealParams
                                        in do
                                        newBestGraphList <- GO.selectGraphs Best (outgroupIndex inGS) (netKeepNum netParams) 0 newGraphList
                                        pure (newBestGraphList, snd5 inPhyloGraph, genNewSimAnnealParams)
                                        

                                    else do
                                        pure ([], minCost, inSimAnnealParams)
                                

                                

{- | insertNetEdgeHeuristic inserts an edge between two other edges, creating 2 new nodes 
contacts deletes 2 orginal edges and adds 2 nodes and 5 new edges
does not check any edge reasonable-ness properties
new edge directed from first to second edge
naive for now
predeletecost of edge move

returns heuristic cost and simple graph. No rediagnosis
If simpleGraph is not PhylogenticGraph then infinity and emptyGraph
-}
insertNetEdgeHeuristic
    ∷ GlobalSettings
    → ProcessedData
    → ReducedPhylogeneticGraph
    → (LG.LEdge b, LG.LEdge b)
    → PhyG (VertexCost, SimpleGraph)
insertNetEdgeHeuristic inGS inData inPhyloGraph edgePair@((u, v, _), (u', v', _)) =
    if LG.isEmpty $ thd5 inPhyloGraph
        then error "Empty input phylogenetic graph in insertNetEdgeHeuristic"
        else
            let inSimple = fst5 inPhyloGraph

                -- get children of u' to make sure no net children--moved to permissiable edges
                -- u'ChildrenNetNodes = filter (== True) $ fmap (LG.isNetworkNode inSimple) $ LG.descendants inSimple u'

                numNodes = length $ LG.nodes inSimple
                newNodeOne = (numNodes, TL.pack ("HTU" <> (show numNodes)))
                newNodeTwo = (numNodes + 1, TL.pack ("HTU" <> (show $ numNodes + 1)))
                newEdgeList =
                    [ (u, fst newNodeOne, 0.0)
                    , (fst newNodeOne, v, 0.0)
                    , (u', fst newNodeTwo, 0.0)
                    , (fst newNodeTwo, v', 0.0)
                    , (fst newNodeOne, fst newNodeTwo, 0.0)
                    ]
                edgesToDelete = [(u, v), (u', v')]
                newSimple = LG.delEdges edgesToDelete $ LG.insEdges newEdgeList $ LG.insNodes [newNodeOne, newNodeTwo] inSimple

            in  do
                    -- remove these checks when working-
                    isPhyloGraph ← LG.isPhylogeneticGraph newSimple
                    
                    -- unfortunetely need this when there are network edges present before add
                    if not isPhyloGraph
                        then do
                            pure (infinity, LG.empty) -- emptyReducedPhylogeneticGraph
                        else do
                            -- calculates heursitic graph delta
                            --(heuristicDelta', _, _, _, _) <-  heuristicAddDelta' inGS inPhyloGraph edgePair (fst newNodeOne) (fst newNodeTwo)

                            heuristicDelta <- heuristicAddDelta inGS inPhyloGraph edgePair
                                              
                            let edgeAddDelta = deltaPenaltyAdjustment inGS inPhyloGraph "add"

                            -- return heuristic cost and simple graph
                            pure (heuristicDelta + edgeAddDelta,  newSimple) 



{- | getPermissibleEdgePairs takes a DecoratedGraph and returns the list of all pairs
of edges that can be joined by a network edge and meet all necessary conditions
-}
-- add in other conditions
--   reproducable--ie not tree node with two net node children--other stuff
getPermissibleEdgePairs ∷ GlobalSettings → DecoratedGraph → PhyG [(LG.LEdge EdgeInfo, LG.LEdge EdgeInfo)]
getPermissibleEdgePairs inGS inGraph =
    if LG.isEmpty inGraph
        then error "Empty input graph in isEdgePairPermissible"
        else
            let edgeList = LG.labEdges inGraph

                -- edges to potentially conenct
                edgePairs = cartProd edgeList edgeList

                -- get coeval node pairs in existing grap
                coevalNodeConstraintList = LG.coevalNodePairs inGraph

                -- parallel
                -- action :: (LNode a, LNode a) -> (LNode a, LNode a, [LNode a], [LNode a], [LNode a], [LNode a])
                action = LG.addBeforeAfterToPair inGraph
            in  do
                    actionPar ← getParallelChunkMap
                    let coevalNodeConstraintList' = actionPar action coevalNodeConstraintList
                    let edgeAction ∷ (LG.LEdge EdgeInfo, LG.LEdge EdgeInfo) → Bool
                        edgeAction = isEdgePairPermissible inGraph coevalNodeConstraintList'

                    edgePar ← getParallelChunkMap
                    let edgeTestList = edgePar edgeAction edgePairs
                    
                    let pairList = fmap fst $ filter ((== True) . snd) $ zip edgePairs edgeTestList

                    pure pairList


-- where f (a, b) = (LG.toEdge a, LG.toEdge b)

{- | isEdgePairPermissible takes a graph and two edges, coeval contraints, and tests whether a
pair of edges can be linked by a new edge and satify three consitions:
   1) neither edge is a network edge
   2) one edge cannot be "before" while the other is "after" in any of the constraint pairs
   3) neither edge is an ancestor or descendent edge of the other (tested via bv of nodes)
the result should apply to a new edge in either direction
new edge to be creted is edge1 -> ege2
Could change to LG.isPhylogeneticGraph
-}
isEdgePairPermissible
    ∷ DecoratedGraph
    → [(LG.LNode a, LG.LNode a, [LG.LNode a], [LG.LNode a], [LG.LNode a], [LG.LNode a])]
    → (LG.LEdge EdgeInfo, LG.LEdge EdgeInfo)
    → Bool
isEdgePairPermissible inGraph constraintList (edge1@(u, v, _), edge2@(u', v', _)) =
    if LG.isEmpty inGraph
        then error "Empty input graph in isEdgePairPermissible"
        else
            if u == u'
                then False
                else
                    if v == v'
                        then False
                        else -- equality implied in above two
                        -- else if LG.toEdge edge1 == LG.toEdge edge2 then False

                            if (LG.isNetworkNode inGraph u) || (LG.isNetworkNode inGraph u')
                                then False
                                else
                                    if (LG.isNetworkLabEdge inGraph edge1) || (LG.isNetworkLabEdge inGraph edge2)
                                        then False
                                        else
                                            if not (LG.meetsAllCoevalConstraintsNodes (fmap removeNodeLabels constraintList) edge1 edge2)
                                                then False
                                                else
                                                    if (isAncDescEdge inGraph edge1 edge2)
                                                        then False
                                                        else -- get children of u' to make sure no net children

                                                            if (not . null) $ filter (== True) $ fmap (LG.isNetworkNode inGraph) $ LG.descendants inGraph u'
                                                                then False
                                                                else True
    where
        removeNodeLabels (a, b, c, d, e, f) = (LG.toNode a, LG.toNode b, fmap LG.toNode c, fmap LG.toNode d, fmap LG.toNode e, fmap LG.toNode f)


{- | isAncDescEdge takes a graph and two edges and examines whethe either edge is the ancestor or descendent of the other
this is done via examination of teh bitvector fields of the node
-}
isAncDescEdge ∷ DecoratedGraph → LG.LEdge EdgeInfo → LG.LEdge EdgeInfo → Bool
isAncDescEdge inGraph (a, _, _) (b, _, _) =
    if LG.isEmpty inGraph
        then error "Empty input graph in isAncDescEdge"
        else
            let aBV = bvLabel $ fromJust $ LG.lab inGraph a
                bBV = bvLabel $ fromJust $ LG.lab inGraph b
            in  -- trace ("IADE: " <> (show (a, aBV, b, bBV, aBV .&. bBV))) (
                if aBV .&. bBV == aBV
                    then True
                    else
                        if aBV .&. bBV == bBV
                            then True
                            else False

{- These heuristics do not seem to work well at all-}

{- | heuristic add delta based on new display tree and delta from existing costs by block--assumming < 0
original edges subtree1 ((u,l),(u,v)) and subtree2 ((u',v'),(u',l')) create a directed edge from
subtree 1 to subtree 2 via
1) Add node x and y, delete edges (u,v) and (u'v') and create edges (u,x), (x,v), (u',y), and (y,v')
2) real cost is the sum of block costs that are lower for new graph versus older
3) heuristic is when new subtree is lower than existing block by block
   so calculate d(u,v) + d(u',v') [existing display tree cost estimate] compared to
   d((union u,v), v') - d(u'.v') [New display tree cost estimate] over blocks
   so blockDelta = if d((union u,v), v') - d(u'.v') < d(u,v) + d(u',v') then d((union u,v), v') - d(u'.v')
                    else 0 [existing better]
   graphDelta = egdeAddDelta (separately calculated) + sum [blockDelta]
   Compare to real delta to check behavior
original subtrees u -> (a,v) and u' -> (v',b)
-}
heuristicAddDelta ∷ GlobalSettings → ReducedPhylogeneticGraph → (LG.LEdge b, LG.LEdge b) → PhyG VertexCost
heuristicAddDelta _ inPhyloGraph ((u, v, _), (u', v', _)) =
    if LG.isEmpty (fst5 inPhyloGraph)
        then error "Empty graph in heuristicAddDelta"
        else
            let a = head $ filter (/= v) $ LG.descendants (fst5 inPhyloGraph) u
                b = head $ filter (/= v') $ LG.descendants (fst5 inPhyloGraph) u'
                blockTrees =
                    fmap V.fromList $
                        fmap (GO.getDecoratedDisplayTreeList (thd5 inPhyloGraph)) $
                            V.zip (fth5 inPhyloGraph) $
                                V.fromList [0 .. (V.length (fft5 inPhyloGraph) - 1)]
                action :: (V.Vector DecoratedGraph, V.Vector CharInfo) → PhyG VertexCost
                action = getBlockDelta (u, v, u', v', a, b)
            in do
                actionPar <- getParallelChunkTraverse 
                -- blockDeltaV <- V.zipWith (getBlockDelta (u, v, u', v', a, b)) (zip blockTrees (fft5 inPhyloGraph))
                blockDeltaL <- actionPar action (V.toList $ V.zip blockTrees (fft5 inPhyloGraph))
                pure $ sum blockDeltaL


{- | getBlockDelta determines the network add delta for each block (vector of characters)
if existing is lower then zero, else (existing - new)
-}
getBlockDelta
    ∷ (LG.Node, LG.Node, LG.Node, LG.Node, LG.Node, LG.Node) → (V.Vector DecoratedGraph, V.Vector CharInfo) → PhyG VertexCost
getBlockDelta (u, v, u', v', a, b) (inCharV, charInfoV) =
    if V.null inCharV
        then error "Empty charcter tree vector in getBlockDelta"
        else
            let action :: (DecoratedGraph, CharInfo) → (VertexCost, VertexCost)
                action = getCharacterDelta (u, v, u', v', a, b)
            in do
                 actionPar ← getParallelChunkMap
                 let result = actionPar action (V.toList $ V.zip inCharV charInfoV)
                 let (charNewV, charExistingV) = unzip result
                 let newCost = sum charNewV
                 let existingCost = sum charExistingV
                 if (newCost < existingCost)
                    then pure $ newCost - existingCost
                    else pure $ 0.0

{- | getCharacterDelta determines the network add delta for each block (vector of characters)
if existing is lower then zero, else (existing - new)
 calculate d(u,v) + d(u',v') [existing display tree cost estimate] compared to
 d((union u,v), v') - d(u'.v')
need to use final assignemnts--so set prelim to final first
Since a distance--no need for No chanage cost adjustment
-}
getCharacterDelta
    ∷ (LG.Node, LG.Node, LG.Node, LG.Node, LG.Node, LG.Node) → (DecoratedGraph, CharInfo) → (VertexCost, VertexCost)
getCharacterDelta (_, v, _, v', a, b) (inCharTree, charInfo) =
    -- getCharacterDelta (u,v,u',v',a,b) inCharTree charInfo =
    let doIA = False
        noChangeCostAdjust = False --this since want a distance not a median
        -- filterGaps = True
        -- uData = V.head $ V.head $ vertData $ fromJust $ LG.lab inCharTree u
        vData = V.head $ V.head $ vertData $ fromJust $ LG.lab inCharTree v
        vFinalData = V.head $ V.head $ PRE.setPreliminaryToFinalStates $ vertData $ fromJust $ LG.lab inCharTree v
        -- u'Data = V.head $ V.head $ vertData $ fromJust $ LG.lab inCharTree u'
        v'Data = V.head $ V.head $ vertData $ fromJust $ LG.lab inCharTree v'
        v'FinalData = V.head $ V.head $ PRE.setPreliminaryToFinalStates $ vertData $ fromJust $ LG.lab inCharTree v'
        aData = V.head $ V.head $ vertData $ fromJust $ LG.lab inCharTree a
        aFinalData = V.head $ V.head $ PRE.setPreliminaryToFinalStates $ vertData $ fromJust $ LG.lab inCharTree a
        bData = V.head $ V.head $ vertData $ fromJust $ LG.lab inCharTree b

        -- unionUV = M.union2Single doIA filterGaps uData vData charInfo
        -- (_,dUV) =  M.median2Single doIA uData vData charInfo
        -- dUV = vertexCost $ fromJust $ LG.lab inCharTree u
        -- dU'V' = vertexCost $ fromJust $ LG.lab inCharTree u'
        -- (_, dUnionUVV') = M.median2Single doIA unionUV v'Data charInfo

        (newX, dVV') = M.median2Single noChangeCostAdjust doIA vFinalData v'FinalData charInfo
        (_, dAX) = M.median2Single noChangeCostAdjust doIA aFinalData newX charInfo
        (_, dAV) = M.median2Single noChangeCostAdjust doIA aData vData charInfo
        (_, dV'B) = M.median2Single noChangeCostAdjust doIA v'Data bData charInfo
    in  -- trace ("GCD: " <> (show (dVV' + dAX, dAV + dV'B))) (
        (dVV' + dAX, dAV + dV'B)


{- This seems worse than heuristicAddDelta above-}

{- | heuristicAddDelta' takes the existing graph, edge pair, and new nodes to create and makes
the new nodes and reoptimizes starting nodes of two edges.  Returns cost delta based on
previous and new node resolution caches
returns cost delta and the reoptimized nodes for use in incremental optimization
original edges (to be deleted) (u,v) and (u',v'), n1 inserted in (u,v) and n2 inserted into (u',v')
creates (n1, n2), (u,n1), (n1,v), (u',n2), (n2, v')
-}
heuristicAddDelta'
    ∷ GlobalSettings
    → ReducedPhylogeneticGraph
    → (LG.LEdge b, LG.LEdge b)
    → LG.Node
    → LG.Node
    → PhyG (VertexCost, LG.LNode VertexInfo, LG.LNode VertexInfo, LG.LNode VertexInfo, LG.LNode VertexInfo)
heuristicAddDelta' inGS inPhyloGraph ((u, v, _), (u', v', _)) n1 n2 =
    if LG.isEmpty (fst5 inPhyloGraph)
        then error "Empty graph in heuristicAddDelta"
        else
            if graphType inGS == HardWired
                then do
                    uvVertData <- M.makeEdgeDataM False True (thd5 inPhyloGraph) (fft5 inPhyloGraph) (u, v, dummyEdge)
                    uvPrimeData <- M.makeEdgeDataM False True (thd5 inPhyloGraph) (fft5 inPhyloGraph) (u', v', dummyEdge)
                    postOrderVertex <- POSW.createVertexDataOverBlocks inGS uvVertData uvPrimeData (fft5 inPhyloGraph) []
                    let hardDelta = V.sum $ fmap V.sum $ fmap (fmap snd) postOrderVertex
                        
                    pure (hardDelta, dummyNode, dummyNode, dummyNode, dummyNode)
                else -- softwired

                    let uLab = fromJust $ LG.lab (thd5 inPhyloGraph) u
                        uPrimeLab = fromJust $ LG.lab (thd5 inPhyloGraph) u'
                        vLab = fromJust $ LG.lab (thd5 inPhyloGraph) v
                        vPrimeLab = fromJust $ LG.lab (thd5 inPhyloGraph) v'
                        uPrimeOtherChild = head $ filter ((/= v') . fst) $ LG.labDescendants (thd5 inPhyloGraph) (u', uPrimeLab)
                        uOtherChild = head $ filter ((/= v) . fst) $ LG.labDescendants (thd5 inPhyloGraph) (u, uLab)
                    in  do
                            -- direction first edge to second so n2 is outdegree 1 to v'
                            let n2Lab = NEW.getOutDegree1VertexSoftWired n2 vPrimeLab (thd5 inPhyloGraph) [n2]
                            uPrimeLabAfter ← NEW.getOutDegree2VertexSoftWired inGS (fft5 inPhyloGraph) u' (n2, n2Lab) uPrimeOtherChild (thd5 inPhyloGraph)
                            n1Lab ← NEW.getOutDegree2VertexSoftWired inGS (fft5 inPhyloGraph) n1 (v, vLab) (n2, n2Lab) (thd5 inPhyloGraph)
                            uLabAfter ← NEW.getOutDegree2VertexSoftWired inGS (fft5 inPhyloGraph) u uOtherChild (n1, n1Lab) (thd5 inPhyloGraph)

                            -- cost of resolutions
                            let (_, uCostBefore) = NEW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uLab)
                            let (_, uPrimeCostBefore) = NEW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uPrimeLab)
                            let (_, uCostAfter) = NEW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uLabAfter)
                            let (_, uPrimeCostAfter) = NEW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uPrimeLabAfter)

                            let addNetDelta = (uCostAfter - uCostBefore) + (uPrimeCostAfter - uPrimeCostBefore)
                            -- trace ("HAD: " <> (show (uCostAfter, uCostBefore, uPrimeCostAfter, uPrimeCostBefore)) <> " -> " <> (show addNetDelta)) $
                            if null (filter ((/= v') . fst) $ LG.labDescendants (thd5 inPhyloGraph) (u', uPrimeLab))
                                || null (filter ((/= v) . fst) $ LG.labDescendants (thd5 inPhyloGraph) (u, uLab))
                                then pure (infinity, dummyNode, dummyNode, dummyNode, dummyNode)
                                else -- this should not happen--should try to create new edges from children of net edges

                                    if (length $ LG.descendants (thd5 inPhyloGraph) u) < 2 || (length $ LG.descendants (thd5 inPhyloGraph) u') < 2
                                        then error ("Outdegree 1 nodes in heuristicAddDelta")
                                        else pure (addNetDelta, (u, uLabAfter), (u', uPrimeLabAfter), (n1, n1Lab), (n2, n2Lab))


{- | deltaPenaltyAdjustment takes number of leaves and Phylogenetic graph and returns a heuristic graph penalty for adding a single network edge
if Wheeler2015Network, this is based on all changes affecting a single block (most permissive) and Wheeler 2015 calculation of penalty
if PMDLGraph -- graph complexity
if NoNetworkPenalty then 0
modification "add" or subtrct to calculate delta
always delta is positive--whether neg or pos is deltermined when used
-}
deltaPenaltyAdjustment
    ∷ GlobalSettings
    → ReducedPhylogeneticGraph
    → String
    → VertexCost
deltaPenaltyAdjustment inGS inGraph modification =
    -- trace ("DPA: entering: " <> (show $ graphFactor inGS)) (
    let numLeaves = numDataLeaves inGS
        edgeCostModel = graphFactor inGS
        (_, _, _, networkNodeList) = LG.splitVertexList (fst5 inGraph)
    in  if edgeCostModel == NoNetworkPenalty
            then -- trace ("DPA: No penalty")
                0.0
            else -- else if length networkNodeList == 0 then
            -- trace ("DPA: No cost")
            --   0.0

                if edgeCostModel == Wheeler2015Network
                    then (snd5 inGraph) / (fromIntegral $ 2 * ((2 * numLeaves) - 2) + (2 * (length networkNodeList)))
                    else
                        if (optimalityCriterion inGS) `elem` [PMDL, SI]
                            then -- trace  ("DPW: In PMDLGraph") (

                                if graphType inGS == Tree
                                    then fst $ IL.head (graphComplexityList inGS)
                                    else
                                        if graphType inGS == SoftWired
                                            then
                                                let currentComplexity = fst $ (graphComplexityList inGS) IL.!!! (length networkNodeList)
                                                    nextComplexity =
                                                        if modification == "add"
                                                            then fst $ (graphComplexityList inGS) IL.!!! ((length networkNodeList) + 1)
                                                            else
                                                                if modification == "delete"
                                                                    then fst $ (graphComplexityList inGS) IL.!!! ((length networkNodeList) - 1)
                                                                    else error ("SoftWired deltaPenaltyAdjustment modification not recognized: " <> modification)
                                                in  abs (currentComplexity - nextComplexity)
                                            else
                                                if graphType inGS == HardWired
                                                    then
                                                        let currentComplexity = snd $ (graphComplexityList inGS) IL.!!! (length networkNodeList)
                                                            nextComplexity =
                                                                if modification == "add"
                                                                    then snd $ (graphComplexityList inGS) IL.!!! ((length networkNodeList) + 1)
                                                                    else
                                                                        if modification == "delete"
                                                                            then snd $ (graphComplexityList inGS) IL.!!! ((length networkNodeList) - 1)
                                                                            else error ("HardWired deltaPenaltyAdjustment modification not recognized: " <> modification)
                                                        in  abs (currentComplexity - nextComplexity)
                                                    else error ("Graph type not yet implemented: " <> (show $ graphType inGS))
                            else -- )

                                if edgeCostModel == Wheeler2023Network
                                    then -- same as W15 for heuristic penalty for single edge
                                        (snd5 inGraph) / (fromIntegral $ 2 * ((2 * numLeaves) - 2) + (2 * (length networkNodeList)))
                                    else error ("Network edge cost model not yet implemented: " <> (show edgeCostModel))


-- )

