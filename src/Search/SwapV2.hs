{- |
Module specifying graph swapping rearrangement functions
-}
module Search.SwapV2 (
    reoptimizeSplitGraphFromVertexIANew,
    reoptimizeSplitGraphFromVertexNew,
    reoptimizeSplitGraphFromVertexTupleNew,
    swapV2,
) where

import Control.Monad (filterM, liftM)
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


{- | New Swap functions that are based on strict PHANE parallelization routines.
    7) unions--only for TBR moves?
        
-}
swapV2
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Int
    → [ReducedPhylogeneticGraph]
    → Maybe SAParams
    → PhyG ([ReducedPhylogeneticGraph], Int)
swapV2 swapParams inGS inData inCounter curBestGraphList saParams =
    if null curBestGraphList
        then pure ([], inCounter)
        else 
            if isNothing saParams then 
                swapByType swapParams inGS inData inCounter curBestGraphList Nothing
            else
                let numberInstances = rounds $ fromJust saParams
                    -- parallel setup
                    action ∷ ReducedPhylogeneticGraph → PhyG ([ReducedPhylogeneticGraph], Int)
                    action = swapSimulatedAnnealing swapParams inGS inData inCounter saParams

                in do
                -- simulated annealing/drift
                -- spawn rounds 
                actionPar <- getParallelChunkTraverse 
                annealGraphPairList <- actionPar action (L.replicate numberInstances (head curBestGraphList))

                -- collect results
                let (annealGraphList, counterList) = unzip annealGraphPairList

                bestGraphs <- GO.selectGraphs Best (keepNum swapParams) 0.0 $ curBestGraphList <> (concat annealGraphList)

                pure (bestGraphs, sum counterList)
                

{- | Controls the simulated annealing/drift swapping phases
    1) swap potentially accepting worse cost graphs according to SA?Drif
    2) swapping back (normally) if graph cost > input best

    This should really be operating only on a single graph
-}
swapSimulatedAnnealing
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Int
    → Maybe SAParams
    → ReducedPhylogeneticGraph
    → PhyG ([ReducedPhylogeneticGraph], Int)
swapSimulatedAnnealing swapParams inGS inData inCounter saParams inBestGraph = do
    --if null inBestGraphList
    --    then pure ([], inCounter)
    --    else do
            let inBestCost = snd5 inBestGraph

            -- this to make SA/Drift  simpler
            let driftSwapType = if swapType swapParams `elem` [NNI, SPR, TBR] then swapType swapParams
                                else TBR

            --logWith LogInfo $ "\tStarting with " <> (show inBestCost)

            -- generate SA/Drift graphs, use "Best" since want the "worse" results to some extent, and minimizes work
            (saGraphs, saCounter) <- swapByType (swapParams {checkHeuristic = BestOnly, swapType = driftSwapType}) inGS inData inCounter [inBestGraph] saParams

            --logWith LogInfo $ "\tGoing back with " <> (show $ fmap snd5 saGraphs)

            -- swap "back" to optimal costs
            (backGraphs, backCounter) <- swapByType swapParams inGS inData saCounter saGraphs Nothing

            -- take SA/Drif or Input solutions based on cost 
            let finalGraphs =   if (minimum $ fmap snd5 backGraphs) < inBestCost then 
                                    backGraphs
                                else [inBestGraph]

            pure (finalGraphs, backCounter)

{- | Controls the use of alternate versus other types of swap
    makes things cleaner for SA/Drift
    No SA/Drift for Alternate--stick with NNI/SPR/TBR celaner
-}
swapByType
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Int
    → [ReducedPhylogeneticGraph]
    → Maybe SAParams
    → PhyG ([ReducedPhylogeneticGraph], Int)
swapByType swapParams inGS inData inCounter curBestGraphList saParams =
    if null curBestGraphList
        then pure ([], inCounter)
        else 
            if swapType swapParams `elem` [NNI, SPR, TBR] then 
                swapNaive swapParams inGS inData inCounter curBestGraphList curBestGraphList saParams
            else if swapType swapParams == Alternate then
                -- alternate between SPR and TBROnly
                swapAlternate swapParams inGS inData inCounter curBestGraphList curBestGraphList
            else error ("Swap type not recognized/implemented: " <> (show $ swapType swapParams))

{- | swapAlternate uses SPR and TBROnly to sear where if TBROnly finds a new solution, it 
    returns to SPR on new graph and does this untill TBROnly finds no new solutions.

    No SAPArams since not used in drift/SA
-}
swapAlternate 
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Int
    → [ReducedPhylogeneticGraph]
    → [ReducedPhylogeneticGraph]
    → PhyG ([ReducedPhylogeneticGraph], Int)
swapAlternate swapParams inGS inData inCounter graphsToSwap curBestGraphList = 
    if null graphsToSwap then pure ([], inCounter)
    else 
        let curBestCost = minimum $ fmap snd5 graphsToSwap
        in do
            (sprGraphResult, sprCount) <- swapNaive (swapParams {swapType = SPR}) inGS inData inCounter curBestGraphList curBestGraphList Nothing
            (tbrOnlyResult, _) <- swapNaive (swapParams {swapType = TBR}) inGS inData sprCount sprGraphResult sprGraphResult Nothing
            
            let sprBestCost = minimum $ fmap snd5 sprGraphResult
            let tbrBestCost = minimum $ fmap snd5 tbrOnlyResult

            -- The to see if TBR rearrangements did anything--if not end
            if tbrBestCost < sprBestCost then
                swapAlternate swapParams inGS inData sprCount tbrOnlyResult tbrOnlyResult 

            else 
                pure (tbrOnlyResult, sprCount)
        
{- | Naive swap functions to create reference for later algorithmic improvements 
    1) Take first graph

    2) Create list of splits of graph 
        LG.getEdgeSplitList O(n)

    3) Rejoin all places for all splits
        This is "all around" in that not switching to lower cost graphs
    
    4) Take new graph at first opportunity (ie. 'steepest").
        
        LG.splitGraphOnEdge (check base graph leaves >=3)
        Current function to split graph and reoptimize components: reoptimizeSplitGraphFromVertex
        LG.joinGraphOnEdge

    4) Full evaluation of graphs 
        Time complexities should be with O(n) for graph traversals
            NNI O(n)
            SPR O(n^2)
            TBR O(n^3)

    SA/Drift uses the recursive call not full parallel on splits
-}
swapNaive 
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Int
    → [ReducedPhylogeneticGraph]
    → [ReducedPhylogeneticGraph]
    → Maybe SAParams
    → PhyG ([ReducedPhylogeneticGraph], Int)
swapNaive swapParams inGS inData inCounter graphsToSwap curBestGraphList saParams =
    let inGraphPair = L.uncons graphsToSwap
    in
    if isNothing inGraphPair
        then pure (curBestGraphList, inCounter)
        else do
            -- if graph list not empty proceed 
            let (firstGraph, graphsRemaining) = fromJust inGraphPair

            -- get list of splittable edges (bridging, not adjascent to root)
            let edgeList' = L.uncons $ LG.getEdgeSplitList $ thd5 firstGraph
            if isNothing edgeList' then do
                                        logWith LogInfo $ "\tNo Breakable Edges"
                                        pure (graphsRemaining, inCounter + 1)
            else do
                let curBestCost = minimum $ fmap snd5 graphsToSwap

                -- do not shortcicuit here based on cost -- need for SA/Drift

                let fullFirstGraph = GO.convertReduced2PhylogeneticGraph firstGraph
                inGraphNetPenalty ← T.getPenaltyFactor inGS inData Nothing fullFirstGraph
                let inGraphNetPenaltyFactor = inGraphNetPenalty / curBestCost

                -- sort overrides random
                -- SA/Drif uses randomized order
                edgeList <- if (atRandom swapParams) && (not $ sortEdgesSplitCost swapParams) || (isJust saParams) then 
                                shuffleList (LG.getEdgeSplitList $ thd5 firstGraph)
                            else pure (LG.getEdgeSplitList $ thd5 firstGraph)

                let nonExactCharacters = U.getNumberSequenceCharacters (thd3 inData)
                
                -- this is doing all splits in parallel so can sort on those with greatest cost difference
                    -- probably should be an option since does alot of work that would discarded in a "steepeast descent"
                    -- though my still be worth it
                    -- could do in batches also
                -- if random order perhaps do parallel in the join phase in tranches
                -- SA/Drift uses the recursive to save work if terminates early
                rejoinResult <- if (splitParallel swapParams) && (isNothing saParams) then do
                                -- splitAction ::  LG.LEdge EdgeInfo → PhyG (DecoratedGraph, VertexCost, LG.Node, LG.Node, LG.Node)
                                    let splitAction = doASplit swapParams inGS inData (doIA swapParams) nonExactCharacters inGraphNetPenaltyFactor fullFirstGraph
                                    spltActionPar <- (getParallelChunkTraverseBy snd5)
                                    resultListP <- spltActionPar splitAction edgeList

                                    -- order of split cost
                                    if sortEdgesSplitCost swapParams then
                                        rejoinFromOptSplitList swapParams inGS inData (doIA swapParams) inGraphNetPenaltyFactor graphsToSwap curBestCost edgeList (L.sortOn snd5 resultListP) -- ([head splitList])

                                    -- could have been randomized order
                                    else rejoinFromOptSplitList swapParams inGS inData (doIA swapParams) inGraphNetPenaltyFactor graphsToSwap curBestCost edgeList resultListP

                                -- this is recursive on splitting as opposed to parallelized
                                -- used by SA/Drift (with randomization) as well to save work
                                else 
                                    doAllSplitsAndRejoin swapParams inGS inData (doIA swapParams) nonExactCharacters inGraphNetPenaltyFactor graphsToSwap curBestCost fullFirstGraph saParams edgeList

                -- Always return the whole lot for SA, since BEstOnly has been specified should be single graph
                (newValList, returnCost) <- if isJust saParams then 
                                                pure (fmap fst rejoinResult, minimum $ fmap snd rejoinResult)
                                            else if (minimum $ fmap snd rejoinResult) < curBestCost then 
                                                pure (fmap fst rejoinResult, minimum $ fmap snd rejoinResult)
                                            else if (minimum $ fmap snd rejoinResult) == curBestCost then 
                                                pure ((fmap fst rejoinResult) <> curBestGraphList, curBestCost)
                                            else pure (curBestGraphList, curBestCost)

                -- Conditions to keep going or return
                -- SA?Drift hit max numbers
                if isJust saParams then
                    pure (newValList, inCounter + 1)
                -- Better graphs than input
                else if returnCost < curBestCost then 
                    if (TBROnly == swapType swapParams) then
                        pure (newValList, inCounter + 1)
                    else
                        swapNaive swapParams inGS inData (inCounter + 1) newValList newValList saParams
                -- equal graph costs to input
                else if returnCost == curBestCost then do
                    uniqueBestGraphs <- GO.selectGraphs Unique (keepNum swapParams) 0.0 $ newValList <> curBestGraphList
                    swapNaive swapParams inGS inData (inCounter + 1) (graphsRemaining L.\\ uniqueBestGraphs) uniqueBestGraphs saParams

                -- worse graphs than input
                else swapNaive swapParams inGS inData (inCounter + 1) graphsRemaining curBestGraphList saParams

{- | rejoinFromOptSplitList\plitList takes a list of optimized split graphs and
    calls rejoin function
    Good for parallel exectution

    Does not have SA/Drift functionality
-}
rejoinFromOptSplitList 
    :: SwapParams 
    -> GlobalSettings 
    → ProcessedData 
    -> Bool 
    -> VertexCost 
    -> [ReducedPhylogeneticGraph]
    -> VertexCost 
    -> [LG.LEdge EdgeInfo] 
    -> [(DecoratedGraph, VertexCost, LG.Node, LG.Node, LG.Node)] 
    -> PhyG [(ReducedPhylogeneticGraph, VertexCost)]
rejoinFromOptSplitList swapParams inGS inData doIA inGraphNetPenaltyFactor curBestGraphList curBestCost splitEdgeList splitInfoList'' =
    let splitInfoList' = L.uncons splitInfoList''
    in
    if isNothing splitInfoList' then pure $ zip curBestGraphList (fmap snd5 curBestGraphList)
    else 
        let (firstSplit, restSplits) = fromJust splitInfoList'
            (splitGraphOptimized, splitCost, graphRoot, prunedGraphRootIndex, originalConnectionOfPruned) = firstSplit

            -- get root in base (for readdition) and edges in pruned section for rerooting during readdition
            (_, edgesInPrunedGraph) = LG.nodesAndEdgesAfter splitGraphOptimized [(originalConnectionOfPruned, fromJust $ LG.lab splitGraphOptimized originalConnectionOfPruned)]

            -- Filter for bridge edges and desendents of root of pruned graph -- for TBR when needed
            -- remeber to add in root of pruned graph data for TBR edge data so SPR is subset TBR
            tbrRerootEdges
                | swapType swapParams `elem` [NNI, SPR] = []
                | (graphType inGS == Tree) || LG.isTree splitGraphOptimized = edgesInPrunedGraph L.\\ ((LG.out splitGraphOptimized originalConnectionOfPruned) <> (LG.out splitGraphOptimized prunedGraphRootIndex))
                | otherwise = (fmap fst . filter snd . zip edgesInPrunedGraph $ LG.isBridge splitGraphOptimized . LG.toEdge <$> edgesInPrunedGraph) L.\\ ((LG.out splitGraphOptimized originalConnectionOfPruned) <> (LG.out splitGraphOptimized prunedGraphRootIndex))

            charInfoVV = fmap thd3 $ thd3 inData

            (_, edgesInBaseGraph) = LG.nodesAndEdgesAfter splitGraphOptimized [(graphRoot, fromJust $ LG.lab splitGraphOptimized graphRoot)]

            -- Functions for edge data and rejoin function
            (makeEdgeDataFunction, edgeJoinFunction) =
                if graphType inGS == HardWired then 
                (M.makeEdgeDataM False True, edgeJoinDelta inGS False)
                else if not doIA then 
                    (M.makeEdgeDataM False True, edgeJoinDelta inGS False)
                else 
                    (M.makeEdgeDataM True True, edgeJoinDelta inGS True)
            
           -- this needs to sytart with prunedGraphRoot
            diagnoseAction :: SimpleGraph → PhyG ReducedPhylogeneticGraph
            diagnoseAction = T.multiTraverseFullyLabelGraphReduced inGS inData False False Nothing

            -- Remake Graphs from joins
            makeSPRGraphAction :: LG.LEdge EdgeInfo -> PhyG DecoratedGraph
            makeSPRGraphAction = makeSprNewGraph inGS splitGraphOptimized prunedGraphRootIndex originalConnectionOfPruned
            
            makeTBRGraphAction :: (LG.LEdge EdgeInfo, LG.LEdge EdgeInfo) -> PhyG DecoratedGraph
            makeTBRGraphAction = makeTBRNewGraph inGS splitGraphOptimized prunedGraphRootIndex originalConnectionOfPruned 

            -- generate TBR reroot info for the pruned graph so can apply to all rejoins.
            -- and pass to the singleJoin for TBR
            -- for heuristic returns edge state and reroot edge for later graph reconstruction if needed
            makeTBREdgeData :: LG.LEdge EdgeInfo → PhyG VertexBlockData
            makeTBREdgeData =  makeEdgeDataFunction splitGraphOptimized charInfoVV

        in do
            if (null tbrRerootEdges) && (swapType swapParams == TBROnly) then
                    rejoinFromOptSplitList swapParams inGS inData doIA inGraphNetPenaltyFactor curBestGraphList curBestCost splitEdgeList restSplits
            else do
                rejoinEdges <- if atRandom swapParams then 
                                 shuffleList edgesInBaseGraph
                               -- should re-add close to original placement first
                               else if swapType swapParams == NNI then
                                    pure $ take 3 $ LG.sortEdgesByIndexDistance originalConnectionOfPruned edgesInBaseGraph
                                                        
                               else 
                                    pure $ LG.sortEdgesByIndexDistance originalConnectionOfPruned edgesInBaseGraph

                
                {-Make TBR EdgeData-}
                makeTBREdgeDataPar <- getParallelChunkTraverse 
                tbrEdgeDataList <- makeTBREdgeDataPar makeTBREdgeData tbrRerootEdges

                -- this will be [] if SPR
                let tbrEdgeDataPairList = zip tbrRerootEdges tbrEdgeDataList

                -- Parallel set up
                --heuristicAction :: LG.LEdge EdgeInfo → PhyG ([(LG.LEdge EdgeInfo, VertexCost)], [(LG.LEdge EdgeInfo, LG.LEdge EdgeInfo, VertexCost)])
                let heuristicAction = singleJoinHeuristic makeEdgeDataFunction edgeJoinFunction swapParams inGS inData splitGraphOptimized splitCost prunedGraphRootIndex tbrEdgeDataPairList

                {-Heuristic cost calculations-}
                heuristicActionPar <- getParallelChunkTraverse 
                heuristicResultList <- heuristicActionPar heuristicAction rejoinEdges

                -- process returns
                -- changed this to keep list of all not better for BestOnly and SA/Drift
                -- BestOnly now gets best evein if not better than current cost
                let bettterHeuristicEdgesSPR = filter ((< curBestCost) . snd) $ concat $ fmap fst heuristicResultList
                let minHeuristicCostSPR =   if (null $ fmap fst heuristicResultList) then infinity
                                            else minimum $ fmap snd $ concat $ fmap fst heuristicResultList -- bettterHeuristicEdgesSPR

                let tbrPart = concat $ fmap snd heuristicResultList
                let bettterHeuristicEdgesTBR = filter ((< curBestCost) . snd) tbrPart
                let minHeuristicCostTBR =   if (null tbrPart) then infinity
                                            else minimum $ fmap snd tbrPart -- $ fmap snd bettterHeuristicEdgesTBR

                let minHeuristicCost =  min minHeuristicCostSPR minHeuristicCostTBR

                let toReoptimizeAndCheckCostSPR =   if (checkHeuristic swapParams == Better) then bettterHeuristicEdgesSPR
                                                    else if (checkHeuristic swapParams == BestOnly) then filter ((== minHeuristicCost) . snd) $ concat $ fmap fst heuristicResultList -- bettterHeuristicEdgesSPR
                                                    else if (checkHeuristic swapParams == BetterN) then take (graphsSteepest inGS) $ L.sortOn snd $ concat $ fmap fst heuristicResultList 
                                                    -- BestAll
                                                    else concat $ fmap fst heuristicResultList

                --- Make full graphs that are needed for reoptimizaiton
                makeSPRGraphPar <- getParallelChunkTraverse
                toReoptimizeGraphsSPR <- makeSPRGraphPar makeSPRGraphAction (fmap fst toReoptimizeAndCheckCostSPR)

                -- process TBR returns (if done)
                let toReoptimizeAndCheckCostTBR =   if (checkHeuristic swapParams == Better) then bettterHeuristicEdgesTBR
                                                    else if (checkHeuristic swapParams == BestOnly) then filter ((== minHeuristicCost) . snd) tbrPart -- bettterHeuristicEdgesTBR
                                                    else if (checkHeuristic swapParams == BetterN) then take (graphsSteepest inGS) $ L.sortOn snd $ tbrPart
                                                    -- BestAll
                                                    else tbrPart

                makeTBRGraphPar <- getParallelChunkTraverse
                toReoptimizeGraphsTBR <- makeTBRGraphPar makeTBRGraphAction (fmap fst toReoptimizeAndCheckCostTBR)

                -- combine graaphs for reoptimization
                let toReoptimizeGraphs = toReoptimizeGraphsSPR <> toReoptimizeGraphsTBR
                --logWith LogInfo $ "\n\tRediagnosing: " <> (show $ length toReoptimizeGraphs) <> " of " <> (show (length heuristicResultList, length tbrPart))
                
                {-Diagnose and check costs-}
                diagnoseActionPar <- (getParallelChunkTraverseBy snd5)
                checkedGraphCosts <- diagnoseActionPar diagnoseAction (fmap GO.convertDecoratedToSimpleGraph toReoptimizeGraphs)

                let minimumCheckedCost = if (null checkedGraphCosts) then infinity
                                         else minimum $ fmap snd5 checkedGraphCosts

                (newBestGraphs, newBestCost) <- if minimumCheckedCost < curBestCost then do
                                                    logWith LogInfo $ "\t-> " <> (show minimumCheckedCost)
                                                    pure (filter ((== minimumCheckedCost) .snd5) checkedGraphCosts, minimumCheckedCost)
                                                else 
                                                    pure (curBestGraphList, curBestCost)
                                            
                if ((swapType swapParams == TBROnly) || steepest swapParams) && minimumCheckedCost < curBestCost then
                    pure $  zip newBestGraphs (replicate (length newBestGraphs) newBestCost)         
                else rejoinFromOptSplitList swapParams inGS inData doIA inGraphNetPenaltyFactor newBestGraphs newBestCost splitEdgeList restSplits


{- doAllSplitsAndRejoin generates all split reoptimizations
    Calls the rejoin function after creating the optimizaed split graph
    recurses to next edge to split
    Good for space conservation
-}
doAllSplitsAndRejoin 
    :: SwapParams 
    -> GlobalSettings 
    → ProcessedData 
    -> Bool 
    -> Int 
    -> VertexCost 
    -> [ReducedPhylogeneticGraph]
    -> VertexCost 
    -> PhylogeneticGraph 
    -> Maybe SAParams
    -> [LG.LEdge EdgeInfo] 
    -> PhyG [(ReducedPhylogeneticGraph, VertexCost)]
doAllSplitsAndRejoin swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor curBestGraphList curBestCost firstFullGraph saParams edgeList'' = 
                let edgeList' = L.uncons edgeList''
                in
                if isNothing edgeList' then pure $ zip curBestGraphList (fmap snd5 curBestGraphList)
                else 
                    let edgeList@(firstEdge, restEdges) = fromJust edgeList'
                    in do
                        
                        -- split graph on the first edge
                        let (splitGraph, graphRoot, prunedGraphRootIndex, originalConnectionOfPruned) = LG.splitGraphOnEdge (thd6 firstFullGraph) firstEdge
                        
                        (splitGraphOptimized, splitCost) ←
                                reoptimizeSplitGraphFromVertexNew swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor firstFullGraph splitGraph graphRoot prunedGraphRootIndex 
                        --
                        --logWith LogInfo $ "\tSplit Cost New: " <> (show splitCost) -- <> "\n" <> LG.prettyDot reoptimizedSplitGraph'
                        -- avoid shorcircuiting on SA/Drif
                        if (splitCost >= curBestCost) && (isNothing saParams) then 
                            doAllSplitsAndRejoin swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor curBestGraphList curBestCost firstFullGraph saParams restEdges
                        else 
                            -- rejoin function

                            -- get root in base (for readdition) and edges in pruned section for rerooting during readdition
                            let (_, edgesInPrunedGraph) = LG.nodesAndEdgesAfter splitGraphOptimized [(originalConnectionOfPruned, fromJust $ LG.lab splitGraphOptimized originalConnectionOfPruned)]

                                -- Filter for bridge edges and desendents of root of pruned graph -- for TBR when needed
                                -- remeber to add in root of pruned graph data for TBR edge data so SPR is subset TBR
                                tbrRerootEdges
                                    | swapType swapParams /= TBR = []
                                    | (graphType inGS == Tree) || LG.isTree splitGraphOptimized = edgesInPrunedGraph L.\\ ((LG.out splitGraphOptimized originalConnectionOfPruned) <> (LG.out splitGraphOptimized prunedGraphRootIndex))
                                    | otherwise = (fmap fst . filter snd . zip edgesInPrunedGraph $ LG.isBridge splitGraphOptimized . LG.toEdge <$> edgesInPrunedGraph) L.\\ ((LG.out splitGraphOptimized originalConnectionOfPruned) <> (LG.out splitGraphOptimized prunedGraphRootIndex))

                                (_, edgesInBaseGraph) = LG.nodesAndEdgesAfter splitGraphOptimized [(graphRoot, fromJust $ LG.lab splitGraphOptimized graphRoot)]

                                charInfoVV = fmap thd3 $ thd3 inData

                                -- Functions for edge data and rejoin function
                                (makeEdgeDataFunction, edgeJoinFunction) =
                                    if graphType inGS == HardWired then 
                                        (M.makeEdgeDataM False True, edgeJoinDelta inGS False)
                                    else if not doIA then 
                                        (M.makeEdgeDataM False True, edgeJoinDelta inGS False)
                                    else 
                                        (M.makeEdgeDataM True True, edgeJoinDelta inGS True)

                                -- this needs to sytart with prunedGraphRoot
                                diagnoseAction :: SimpleGraph → PhyG ReducedPhylogeneticGraph
                                diagnoseAction = T.multiTraverseFullyLabelGraphReduced inGS inData False False Nothing

                                -- Remake Graphs from joins
                                makeSPRGraphAction :: LG.LEdge EdgeInfo -> PhyG DecoratedGraph
                                makeSPRGraphAction = makeSprNewGraph inGS splitGraphOptimized prunedGraphRootIndex originalConnectionOfPruned

                                makeTBRGraphAction :: (LG.LEdge EdgeInfo, LG.LEdge EdgeInfo) -> PhyG DecoratedGraph
                                makeTBRGraphAction = makeTBRNewGraph inGS splitGraphOptimized prunedGraphRootIndex originalConnectionOfPruned 
                                -- generate TBR reroot info for the pruned graph so can apply to all rejoins.
                                -- and pass to the singleJoin for TBR
                                -- for heuristic returns edge state and reroot edge for later graph reconstruction if needed
                                makeTBREdgeData :: LG.LEdge EdgeInfo → PhyG VertexBlockData
                                makeTBREdgeData =  makeEdgeDataFunction splitGraphOptimized charInfoVV


                            in do
                                rejoinEdges <- if atRandom swapParams then 
                                                 shuffleList edgesInBaseGraph
                                               -- should re-add close to original placement first
                                               else if swapType swapParams == NNI then
                                                    pure $ take 3 $ LG.sortEdgesByIndexDistance originalConnectionOfPruned edgesInBaseGraph
                                                                        
                                               else 
                                                    pure $ LG.sortEdgesByIndexDistance originalConnectionOfPruned edgesInBaseGraph

                                -- shot circuit for TBTROnly in Alternate
                                if (null tbrRerootEdges) && (swapType swapParams == TBROnly) then
                                    doAllSplitsAndRejoin swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor curBestGraphList curBestCost firstFullGraph saParams restEdges

                                --Make TBR EdgeData
                                else do 
                                    makeTBREdgeDataPar <- getParallelChunkTraverse 
                                    tbrEdgeDataList <- makeTBREdgeDataPar makeTBREdgeData tbrRerootEdges

                                    let tbrEdgeDataPairList = zip tbrRerootEdges tbrEdgeDataList

                                    -- Parallel set up
                                    --heuristicAction :: LG.LEdge EdgeInfo → PhyG ([(LG.LEdge EdgeInfo, VertexCost)], [(LG.LEdge EdgeInfo, LG.LEdge EdgeInfo, VertexCost)])
                                    let heuristicAction = singleJoinHeuristic makeEdgeDataFunction edgeJoinFunction swapParams inGS inData splitGraphOptimized splitCost prunedGraphRootIndex tbrEdgeDataPairList


                                    {-Heuristic cost calculations-}
                                    heuristicActionPar <- getParallelChunkTraverse 
                                    heuristicResultList <- heuristicActionPar heuristicAction rejoinEdges

                                    -- process returns
                                    -- changed this to keep list of all not better for BestOnly and SA/Drift
                                    -- BestOnly now gets best evein if not better than current cost
                                    let bettterHeuristicEdgesSPR = filter ((< curBestCost) . snd) $ concat $ fmap fst heuristicResultList
                                    let minHeuristicCostSPR =   if (null $ fmap fst heuristicResultList) then infinity
                                                                else minimum $ fmap snd $ concat $ fmap fst heuristicResultList -- bettterHeuristicEdgesSPR

                                    let tbrPart = concat $ fmap snd heuristicResultList
                                    let bettterHeuristicEdgesTBR = filter ((< curBestCost) . snd) tbrPart
                                    let minHeuristicCostTBR =   if (null tbrPart) then infinity
                                                                else minimum $ fmap snd tbrPart -- $ fmap snd bettterHeuristicEdgesTBR

                                    let minHeuristicCost =  min minHeuristicCostSPR minHeuristicCostTBR

                                    let toReoptimizeAndCheckCostSPR =   if (checkHeuristic swapParams == Better) then bettterHeuristicEdgesSPR
                                                                        else if (checkHeuristic swapParams == BestOnly) then filter ((== minHeuristicCost) . snd) $ concat $ fmap fst heuristicResultList -- bettterHeuristicEdgesSPR
                                                                        else if (checkHeuristic swapParams == BetterN) then take (graphsSteepest inGS) $ L.sortOn snd $ concat $ fmap fst heuristicResultList 
                                                                        -- BestAll
                                                                        else concat $ fmap fst heuristicResultList

                                    --- Make full graphs that are needed for reoptimizaiton
                                    makeSPRGraphPar <- getParallelChunkTraverse
                                    toReoptimizeGraphsSPR <- makeSPRGraphPar makeSPRGraphAction (fmap fst toReoptimizeAndCheckCostSPR)

                                    -- process TBR returns (if done)
                                    let toReoptimizeAndCheckCostTBR =   if (checkHeuristic swapParams == Better) then bettterHeuristicEdgesTBR
                                                                        else if (checkHeuristic swapParams == BestOnly) then filter ((== minHeuristicCost) . snd) tbrPart -- bettterHeuristicEdgesTBR
                                                                        else if (checkHeuristic swapParams == BetterN) then take (graphsSteepest inGS) $ L.sortOn snd $ tbrPart
                                                                        -- BestAll
                                                                        else tbrPart

                                    makeTBRGraphPar <- getParallelChunkTraverse
                                    toReoptimizeGraphsTBR <- makeTBRGraphPar makeTBRGraphAction (fmap fst toReoptimizeAndCheckCostTBR)

                                    -- combine graaphs for reoptimization
                                    let toReoptimizeGraphs = toReoptimizeGraphsSPR <> toReoptimizeGraphsTBR

                                    {-Diagnose and check costs-}
                                    diagnoseActionPar <- (getParallelChunkTraverseBy snd5)
                                    checkedGraphCosts <- diagnoseActionPar diagnoseAction (fmap GO.convertDecoratedToSimpleGraph toReoptimizeGraphs)

                                    let minimumCheckedCost = if (null checkedGraphCosts) then infinity
                                                             else minimum $ fmap snd5 checkedGraphCosts

                                    -- SA/Drift process graphs
                                    if isJust saParams then do
                                        (acceptGraph, newSAParams) <- U.simAnnealAccept saParams curBestCost minimumCheckedCost
                                        (newBestGraphs, newBestCost) <- if acceptGraph then do
                                                                            logWith LogInfo $ "\t-> " <> (show minimumCheckedCost)
                                                                            pure (filter ((== minimumCheckedCost) .snd5) checkedGraphCosts, minimumCheckedCost)
                                                                            
                                                                        else 
                                                                            pure (curBestGraphList, curBestCost)

                                        if U.isSimAnnealTerminated saParams then
                                            pure $ zip newBestGraphs (replicate (length newBestGraphs) newBestCost)
                                        else
                                            doAllSplitsAndRejoin swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor newBestGraphs curBestCost firstFullGraph newSAParams restEdges


                                    -- regular swap
                                    else do
                                        (newBestGraphs, newBestCost) <- if minimumCheckedCost < curBestCost then do
                                                                            logWith LogInfo $ "\t-> " <> (show minimumCheckedCost)
                                                                            pure (filter ((== minimumCheckedCost) .snd5) checkedGraphCosts, minimumCheckedCost)
                                                                        else 
                                                                            pure (curBestGraphList, curBestCost)

                                        if ((swapType swapParams == TBROnly) || steepest swapParams) && minimumCheckedCost < curBestCost then
                                            pure $  zip newBestGraphs (replicate (length newBestGraphs) newBestCost)         
                                        
                                        else doAllSplitsAndRejoin swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor curBestGraphList curBestCost firstFullGraph saParams restEdges
                        

{- | makeSprNewGraph takes split graph and readded edge making new complete graph for rediagnosis etc
-}
makeSprNewGraph :: GlobalSettings -> DecoratedGraph -> LG.Node -> LG.Node -> LG.LEdge EdgeInfo -> PhyG DecoratedGraph
makeSprNewGraph inGS splitGraph prunedGraphRootIndex originalConnectionOfPruned targetEdge@(u, v, _) =
    -- create new edges
    let newEdgeList =
                    [ (u, originalConnectionOfPruned, dummyEdge)
                    , (originalConnectionOfPruned, v, dummyEdge)
                    --, (originalConnectionOfPruned, prunedGraphRootIndex, dummyEdge)
                    ]
        sprNewGraph = LG.insEdges newEdgeList $ LG.delEdges [(u, v)] splitGraph
    in do
        getCheckedGraphNewSPR <-
            if graphType inGS == Tree then pure sprNewGraph
            else do
                isPhyloGraph ← LG.isPhylogeneticGraph sprNewGraph
                if isPhyloGraph then pure sprNewGraph
                else pure LG.empty    
    
        pure getCheckedGraphNewSPR


{- | makeTBRNewGraph takes split graph, rerooted edge and readded edge making new complete graph for rediagnosis etc
-}
makeTBRNewGraph :: GlobalSettings -> DecoratedGraph -> LG.Node -> LG.Node -> (LG.LEdge EdgeInfo, LG.LEdge EdgeInfo) -> PhyG DecoratedGraph
makeTBRNewGraph inGS splitGraph prunedGraphRootIndex originalConnectionOfPruned (targetEdge@(u, v, _), rerootEdge) =

    -- pruned graph rerooted edges
    let (prunedEdgesToAdd, prunedEdgesToDelete) = getTBREdgeEditsDec splitGraph prunedGraphRootIndex rerootEdge

    -- create new edges readdition edges
        newEdgeList =
                    [ (u, originalConnectionOfPruned, dummyEdge)
                    , (originalConnectionOfPruned, v, dummyEdge)
                    --, (originalConnectionOfPruned, prunedGraphRootIndex, dummyEdge)
                    ]
        tbrNewGraph =
            LG.insEdges (newEdgeList <> prunedEdgesToAdd) $
                LG.delEdges ((u, v) : prunedEdgesToDelete) splitGraph
    in do
        getCheckedGraphNewTBR <-
            if graphType inGS == Tree then pure tbrNewGraph
            else do
                isPhyloGraph ← LG.isPhylogeneticGraph tbrNewGraph
                if isPhyloGraph then pure tbrNewGraph
                else pure LG.empty    
    
        pure getCheckedGraphNewTBR


{- | getTBREdgeEditsDec takes and edge and returns the list of edit to pruned subgraph
as a pair of edges to add and those to delete
since reroot edge is directed (e,v), edges away from v will have correct
orientation. Edges between 'e' and the root will have to be flipped
original root edges and reroort edge are deleted and new root and edge spanning orginal root created
delete original connection edge and creates a new one--like SPR
returns ([add], [delete])
-}
getTBREdgeEditsDec ∷ DecoratedGraph → LG.Node → LG.LEdge EdgeInfo → ([LG.LEdge EdgeInfo], [LG.Edge])
getTBREdgeEditsDec inGraph prunedGraphRootIndex rerootEdge =
    -- trace ("Getting TBR Edits for " <> (show rerootEdge)) (
    let -- originalRootEdgeNodes = LG.descendants inGraph prunedGraphRootIndex
        originalRootEdges = LG.out inGraph prunedGraphRootIndex

        -- get path from new root edge fst vertex to orginal root and flip those edges
        -- since (u,v) is u -> v u "closer" to root
        closerToPrunedRootEdgeNode = (fst3 rerootEdge, fromJust $ LG.lab inGraph $ fst3 rerootEdge)
        (nodesInPath, edgesinPath) =
            LG.postOrderPathToNode inGraph closerToPrunedRootEdgeNode (prunedGraphRootIndex, fromJust $ LG.lab inGraph prunedGraphRootIndex)

        -- don't want original root edges to be flipped since later deleted
        edgesToFlip = edgesinPath L.\\ originalRootEdges
        flippedEdges = fmap LG.flipLEdge edgesToFlip

        -- new edges on new root position and spanning old root
        -- add in closer vertex to root to make sure direction of edge is correct
        newEdgeOnOldRoot =
            if (snd3 $ head originalRootEdges) `elem` ((fst3 rerootEdge) : (fmap fst nodesInPath))
                then (snd3 $ head originalRootEdges, snd3 $ last originalRootEdges, dummyEdge)
                else (snd3 $ last originalRootEdges, snd3 $ head originalRootEdges, dummyEdge)
        newRootEdges = [(prunedGraphRootIndex, fst3 rerootEdge, dummyEdge), (prunedGraphRootIndex, snd3 rerootEdge, dummyEdge)]
    in  -- assumes we are not checking original root
        -- rerooted
        -- delete orignal root edges and rerootEdge
        -- add new root edges
        -- and new edge on old root--but need orientation
        -- flip edges from new root to old (delete and add list)
       
        (newEdgeOnOldRoot : (flippedEdges <> newRootEdges), LG.toEdge rerootEdge : (fmap LG.toEdge (edgesToFlip <> originalRootEdges)))


{- | singleJoinHeuristic "rejoins" a pruned graph to base graph based on an edge target in the base graph
    returns the joined graph and heuristic graph cost
    DOES NOT fully optimize graph

    return is pair of SPR result and list of tbr results (if any)

    if swaptype is SPR then the tbrReRootEdgePairList is [] but can be in 
        a TBR search for many splits

    if swapType is TBROnly--does not perform the spr rejoin, for use 
        in Alternate type swap procedure

-}
singleJoinHeuristic 
    :: (DecoratedGraph → V.Vector (V.Vector CharInfo) → LG.LEdge EdgeInfo → PhyG VertexBlockData)
    -> (V.Vector (V.Vector CharInfo) → VertexBlockData → VertexBlockData → PhyG VertexCost)
    -> SwapParams
    -> GlobalSettings
    → ProcessedData
    → DecoratedGraph
    → VertexCost
    → LG.Node
    -> [(LG.LEdge EdgeInfo, VertexBlockData)]
    → LG.LEdge EdgeInfo
    → PhyG ([(LG.LEdge EdgeInfo, VertexCost)], [((LG.LEdge EdgeInfo, LG.LEdge EdgeInfo), VertexCost)])
singleJoinHeuristic makeEdgeDataFunction edgeJoinFunction swapParams inGS inData splitGraph splitCost prunedGraphRootIndex tbrReRootEdgePairList targetEdge =
    let charInfoVV = fmap thd3 $ thd3 inData
    in do
            if swapType swapParams == TBROnly && null tbrReRootEdgePairList then 
                pure ([],[])
            else do
                -- Create union-thype data for target edge
                targetEdgeData ← makeEdgeDataFunction splitGraph charInfoVV targetEdge

                -- get the data from the root of the pruned graph (or leaf if only that)
                let prunedRootVertexData = vertData $ fromJust $ LG.lab splitGraph prunedGraphRootIndex

                -- calculate heuristics join cost
                sprReJoinCost ← edgeJoinFunction charInfoVV prunedRootVertexData targetEdgeData

                
                -- SPR return new graph with heuristic cost
                if swapType swapParams `elem` [NNI, SPR] then 
                    pure ([(targetEdge, splitCost + sprReJoinCost)], []) 

                -- tbr (includes SPR result) and TBROnly which does not
                else if swapType swapParams `elem` [TBR, TBROnly] then 
                    let -- Parallel set up
                        joinAction :: VertexBlockData → PhyG VertexCost
                        joinAction = edgeJoinFunction charInfoVV targetEdgeData 

                    in do
                        {-Heuristic cost calculations-}
                        joinActionPar <- getParallelChunkTraverse 
                        joinActionResultList <- joinActionPar joinAction (fmap snd tbrReRootEdgePairList)

                        let tbrHeuristicCosts = fmap (+ splitCost) joinActionResultList

                        let tbrResult = zip (zip (replicate (length tbrReRootEdgePairList) targetEdge) (fmap fst tbrReRootEdgePairList)) tbrHeuristicCosts

                        -- logWith LogInfo $ " TBR reroots: " <> (show $ length tbrReRootEdgePairList) <> " " <> (show tbrHeuristicCosts)
                        if swapType swapParams == TBR then
                            pure ([(targetEdge, splitCost + sprReJoinCost)], tbrResult)   
                        else 
                            pure ([], tbrResult)
                else error ("Unimplemented swap type: " <> (show $ swapType swapParams)) 

{- | edgeJoinDelta calculates heuristic cost for joining pair edges
    IA field is faster--but has to be there and not for harwired
-}
edgeJoinDelta ∷ GlobalSettings → Bool → V.Vector (V.Vector CharInfo) → VertexBlockData → VertexBlockData → PhyG VertexCost
edgeJoinDelta inGS doIA charInfoVV edgeA edgeB =
    if (not doIA)
        then do
            vertexStuff ← POSW.createVertexDataOverBlocks inGS edgeA edgeB charInfoVV []
            pure $ V.sum $ fmap V.sum $ fmap (fmap snd) vertexStuff
        else do
            vertexStuff ← POSW.createVertexDataOverBlocksStaticIA inGS edgeA edgeB charInfoVV []
            pure $ V.sum $ fmap V.sum $ fmap (fmap snd) vertexStuff

{- doASplit a test function to check a single split reoptimizations
-}
doASplit 
    :: SwapParams 
    -> GlobalSettings 
    → ProcessedData 
    -> Bool 
    -> Int
    -> VertexCost 
    -> PhylogeneticGraph 
    -> (LG.LEdge EdgeInfo) 
    -> PhyG (DecoratedGraph, VertexCost, LG.Node, LG.Node, LG.Node)
doASplit swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor firstFullGraph edgeToSplit = 
    let (splitGraph, graphRoot, prunedGraphRootIndex, originalConnectionOfPruned) = LG.splitGraphOnEdge (thd6 firstFullGraph) edgeToSplit
    in do
        (splitGraphOptimized, splitCost) <- reoptimizeSplitGraphFromVertexNew swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor firstFullGraph splitGraph graphRoot prunedGraphRootIndex 
        pure (splitGraphOptimized, splitCost, graphRoot, prunedGraphRootIndex, originalConnectionOfPruned)

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
    -> Int
    → VertexCost
    → PhylogeneticGraph
    → DecoratedGraph
    → LG.Node
    → LG.Node
    → PhyG (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertexNew swapParams inGS inData doIA nonExactCharacters netPenaltyFactor curGraph inSplitGraph startVertex prunedSubGraphRootVertex =
    -- trace ("RSGFV: " <> (show startVertex)) (
    if doIA
        then -- only reoptimize the IA states for dynamic characters
            reoptimizeSplitGraphFromVertexIANew swapParams inGS inData nonExactCharacters netPenaltyFactor curGraph inSplitGraph startVertex prunedSubGraphRootVertex 
        else -- perform full optimizations of nodes
            
            -- determine position to start incremental optimization for base graphs
            -- grandparent of pruned root due to recontruction of base graph--so use orig graph
            let parentPoint = L.uncons $ concat $ fmap (LG.parents (thd6 curGraph)) $ LG.parents (thd6 curGraph) prunedSubGraphRootVertex
                
                parentOriginalConnection = if isNothing parentPoint then error ("Split graph involving near-root edges: " <> (show prunedSubGraphRootVertex))
                                           else if (not . null . snd) $ fromJust parentPoint then error ("Split graph involving network edge: " <> (show (prunedSubGraphRootVertex, snd $ fromJust parentPoint)))
                                           else fst $ fromJust parentPoint

                -- these required for full optimization

                --nonExactCharacters = U.getNumberSequenceCharacters (thd3 inData)
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
                            (Just (inSplitGraph, parentOriginalConnection))
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
                    -- logWith LogInfo $ "\n\tBase/pruned costs new: " <> (show (snd6 fullBaseGraph, prunedCost))
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
    -> Int
    → VertexCost
    → PhylogeneticGraph
    → DecoratedGraph
    → LG.Node
    → LG.Node
    → PhyG (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertexIANew swapParams inGS inData nonExactCharacters netPenaltyFactor curGraph inSplitGraph startVertex prunedSubGraphRootVertex =
    -- if graphType inGS /= Tree then error "Networks not yet implemented in reoptimizeSplitGraphFromVertexIA"
    -- else
    let origGraph = inSplitGraph -- thd5 origPhyloGraph

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
                    (Just (inSplitGraph, parentOriginalConnection))
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
    -> Int
    → VertexCost
    → (PhylogeneticGraph, DecoratedGraph, LG.Node, LG.Node)
    → PhyG (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertexTupleNew swapParams inGS inData doIA nonExactCharacters netPenaltyFactor (curGraph, inSplitGraph, startVertex, prunedSubGraphRootVertex) =
    reoptimizeSplitGraphFromVertexNew swapParams inGS inData doIA nonExactCharacters netPenaltyFactor curGraph inSplitGraph startVertex prunedSubGraphRootVertex 



