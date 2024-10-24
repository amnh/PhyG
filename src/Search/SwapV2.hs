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
        else swapNaive swapParams inGS inData inCounter curBestGraphList curBestGraphList inSimAnnealParams


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
    → [ReducedPhylogeneticGraph]
    → [(Maybe SAParams, ReducedPhylogeneticGraph)]
    → PhyG ([ReducedPhylogeneticGraph], Int)
swapNaive swapParams inGS inData inCounter graphsToSwap curBestGraphList inSimAnnealParams =
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

                -- skip if already have better graphs
                if curBestCost < snd5 firstGraph then
                    swapNaive swapParams inGS inData (inCounter + 1) graphsRemaining curBestGraphList inSimAnnealParams

                else do

                    let fullFirstGraph = GO.convertReduced2PhylogeneticGraph firstGraph
                    inGraphNetPenalty ← T.getPenaltyFactor inGS inData Nothing fullFirstGraph
                    let inGraphNetPenaltyFactor = inGraphNetPenalty / curBestCost

                    -- sort overrides random
                    edgeList <- if (atRandom swapParams) && (not $ sortEdgesSplitCost swapParams) then 
                                    shuffleList (LG.getEdgeSplitList $ thd5 firstGraph)
                                else pure (LG.getEdgeSplitList $ thd5 firstGraph)

                    let nonExactCharacters = U.getNumberSequenceCharacters (thd3 inData)
                
                    -- this is doing all splits in parallel so can sort on those with greatest cost difference
                        -- probably should be an option since does alot of work that would discarded in a "steepeast descent"
                        -- though my still be worth it
                        -- could do in batches also
                    -- if random order perhaps do parallel in the join phase in tranches
                    rejoinResult <- if splitParallel swapParams then do
                                        -- splitAction ::  LG.LEdge EdgeInfo → PhyG (DecoratedGraph, VertexCost, LG.Node, LG.Node, LG.Node)
                                        let splitAction = doASplit swapParams inGS inData (doIA swapParams) nonExactCharacters inGraphNetPenaltyFactor fullFirstGraph
                                        spltActionPar <- (getParallelChunkTraverseBy snd5)
                                        resultListP <- spltActionPar splitAction edgeList

                                        -- order of split cost
                                        if sortEdgesSplitCost swapParams then
                                            rejoinFromOptSplitList swapParams inGS inData (doIA swapParams) inGraphNetPenaltyFactor graphsToSwap curBestCost edgeList (L.sortOn snd5 resultListP) -- ([head splitList])

                                        -- could have been randomized order
                                        else rejoinFromOptSplitList swapParams inGS inData (doIA swapParams) inGraphNetPenaltyFactor graphsToSwap curBestCost edgeList resultListP

                                    else -- this is recursive on splitting as opposed to paeallelized
                                        doAllSplitsAndRejoin swapParams inGS inData (doIA swapParams) nonExactCharacters inGraphNetPenaltyFactor graphsToSwap curBestCost fullFirstGraph edgeList
                    
                    {- This is recursive with joins inside splits
                    result <- doAllSplits swapParams inGS inData (doIA swapParams) inGraphNetPenaltyFactor fullFirstGraph edgeList
                    -}
                    (newValList, returnCost) <-  if (minimum $ fmap snd rejoinResult) < curBestCost then 
                                                    pure (fmap fst rejoinResult, minimum $ fmap snd rejoinResult )
                                                else if (minimum $ fmap snd rejoinResult) == curBestCost then 
                                                    pure ((fmap fst rejoinResult) <> curBestGraphList, curBestCost)
                                                else pure (curBestGraphList, curBestCost)

                    if returnCost < curBestCost then 
                        swapNaive swapParams inGS inData (inCounter + 1) newValList newValList inSimAnnealParams
                    else if returnCost == curBestCost then do
                        uniqueBestGraphs <- GO.selectGraphs Unique (keepNum swapParams) 0.0 $ newValList <> curBestGraphList
                        swapNaive swapParams inGS inData (inCounter + 1) (graphsRemaining L.\\ uniqueBestGraphs) uniqueBestGraphs inSimAnnealParams
                    else swapNaive swapParams inGS inData (inCounter + 1) graphsRemaining curBestGraphList inSimAnnealParams

{- | rejoinFromOptSplitList\plitList takes a list of optimized split graphs and
    calls rejoin function

    Good for parallel exectution
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
                | swapType swapParams /= TBR = []
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
            rejoinEdges <- if atRandom swapParams then 
                             shuffleList edgesInBaseGraph
                           else 
                             -- should re-add close to original placement first 
                             pure $ LG.sortEdgesByIndexDistance originalConnectionOfPruned edgesInBaseGraph


            {-Make TBR EdgeData-}
            makeTBREdgeDataPar <- getParallelChunkTraverse 
            tbrEdgeDataList <- makeTBREdgeDataPar makeTBREdgeData tbrRerootEdges

            -- this will be [] if SPR
            let tbrEdgeDataPairList = zip tbrRerootEdges tbrEdgeDataList

            -- Parallel set up
            --heuristicAction :: LG.LEdge EdgeInfo → PhyG [(LG.LEdge EdgeInfo, VertexCost)]
            let heuristicAction = singleJoinHeuristic makeEdgeDataFunction edgeJoinFunction swapParams inGS inData splitGraphOptimized splitCost prunedGraphRootIndex tbrEdgeDataPairList

            {-Heuristic cost calculations-}
            heuristicActionPar <- getParallelChunkTraverse 
            heuristicResultList <- heuristicActionPar heuristicAction rejoinEdges

            -- process SPR returns
            let bettterHeuristicEdgesSPR = filter ((< curBestCost) . snd) $ fmap fst heuristicResultList
            let minHeuristicCostSPR =   if (null bettterHeuristicEdgesSPR) then infinity
                                        else minimum $ fmap snd bettterHeuristicEdgesSPR

            let tbrPart = concat $ fmap snd heuristicResultList
            let bettterHeuristicEdgesTBR = filter ((< curBestCost) . snd) tbrPart
            let minHeuristicCostTBR =   if (null bettterHeuristicEdgesTBR) then infinity
                                        else minimum $ fmap snd bettterHeuristicEdgesTBR

            let minHeuristicCost =  min minHeuristicCostSPR minHeuristicCostTBR

            let toReoptimizeAndCheckCostSPR =   if (checkHeuristic swapParams == Better) then 
                                                    bettterHeuristicEdgesSPR
                                                else if (checkHeuristic swapParams == BestOnly) then 
                                                    filter ((== minHeuristicCost) . snd) bettterHeuristicEdgesSPR
                                                else if (checkHeuristic swapParams == BetterN) then 
                                                    take (graphsSteepest inGS) $ L.sortOn snd $ fmap fst heuristicResultList 
                                                -- BestAll
                                                else fmap fst heuristicResultList

            --logWith LogInfo $ "\n\tTaking:" <> (show  ((length $ take (graphsSteepest inGS) $ L.sortOn snd $ fmap fst heuristicResultList), (length $ take (graphsSteepest inGS) $ L.sortOn snd $ tbrPart)))
                                                    
            --- Make full graphs that are needed for reoptimizaiton
            makeSPRGraphPar <- getParallelChunkTraverse
            toReoptimizeGraphsSPR <- makeSPRGraphPar makeSPRGraphAction (fmap fst toReoptimizeAndCheckCostSPR)

            -- process TBR returns (if done)
            let toReoptimizeAndCheckCostTBR =   if (checkHeuristic swapParams == Better) then bettterHeuristicEdgesTBR
                                                else if (checkHeuristic swapParams == BestOnly) then filter ((== minHeuristicCost) . snd) bettterHeuristicEdgesTBR
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

            let minimimCheckedCost = if (null checkedGraphCosts) then infinity
                                     else minimum $ fmap snd5 checkedGraphCosts

            (newBestGraphs, newBestCost) <- if minimimCheckedCost < curBestCost then do
                                                logWith LogInfo $ "\t-> " <> (show minimimCheckedCost)
                                                pure (filter ((== minimimCheckedCost) .snd5) checkedGraphCosts, minimimCheckedCost)
                                            else 
                                                pure (curBestGraphList, curBestCost)
                                        
            --logWith LogInfo $ "\n\tcurBestCost: " <> (show (curBestCost, splitCost)) <> " " <> (show $ zip (fmap snd $ concat heuristicResultList) (fmap snd5 checkedGraphCosts))
            if steepest swapParams && minimimCheckedCost < curBestCost then
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
    -> [LG.LEdge EdgeInfo] 
    -> PhyG [(ReducedPhylogeneticGraph, VertexCost)]
doAllSplitsAndRejoin swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor curBestGraphList curBestCost firstFullGraph edgeList'' = 
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

                        if (splitCost >= curBestCost) then 
                            doAllSplitsAndRejoin swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor curBestGraphList curBestCost firstFullGraph restEdges
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

                                {-
                                -- Parallel set up
                                heuristicAction :: LG.LEdge EdgeInfo → PhyG [(LG.LEdge EdgeInfo, VertexCost)]
                                heuristicAction = singleJoinHeuristic makeEdgeDataFunction edgeJoinFunction swapParams inGS inData splitGraphOptimized splitCost prunedGraphRootIndex originalConnectionOfPruned
                                -}

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
                                               else 
                                                -- should re-add close to original placement first 
                                                pure $ LG.sortEdgesByIndexDistance originalConnectionOfPruned edgesInBaseGraph


                                {-Make TBR EdgeData-}
                                makeTBREdgeDataPar <- getParallelChunkTraverse 
                                tbrEdgeDataList <- makeTBREdgeDataPar makeTBREdgeData tbrRerootEdges

                                let tbrEdgeDataPairList = zip tbrRerootEdges tbrEdgeDataList

                                -- Parallel set up
                                --heuristicAction :: LG.LEdge EdgeInfo → PhyG [(LG.LEdge EdgeInfo, VertexCost)]
                                let heuristicAction = singleJoinHeuristic makeEdgeDataFunction edgeJoinFunction swapParams inGS inData splitGraphOptimized splitCost prunedGraphRootIndex tbrEdgeDataPairList


                                {-Heuristic cost calculations-}
                                heuristicActionPar <- getParallelChunkTraverse 
                                heuristicResultList <- heuristicActionPar heuristicAction rejoinEdges

                                -- process SPR returns
                                let bettterHeuristicEdgesSPR = filter ((< curBestCost) . snd) $ fmap fst heuristicResultList
                                let minHeuristicCostSPR =   if (null bettterHeuristicEdgesSPR) then infinity
                                                            else minimum $ fmap snd bettterHeuristicEdgesSPR

                                let tbrPart = concat $ fmap snd heuristicResultList
                                let bettterHeuristicEdgesTBR = filter ((< curBestCost) . snd) tbrPart
                                let minHeuristicCostTBR =   if (null bettterHeuristicEdgesTBR) then infinity
                                                            else minimum $ fmap snd bettterHeuristicEdgesTBR

                                let minHeuristicCost =  min minHeuristicCostSPR minHeuristicCostTBR

                                let toReoptimizeAndCheckCostSPR =   if (checkHeuristic swapParams == Better) then bettterHeuristicEdgesSPR
                                                                    else if (checkHeuristic swapParams == BestOnly) then filter ((== minHeuristicCost) . snd) bettterHeuristicEdgesSPR
                                                                    else if (checkHeuristic swapParams == BetterN) then take (graphsSteepest inGS) $ L.sortOn snd $ fmap fst heuristicResultList 
                                                                    -- BestAll
                                                                    else fmap fst heuristicResultList

                                --- Make full graphs that are needed for reoptimizaiton
                                makeSPRGraphPar <- getParallelChunkTraverse
                                toReoptimizeGraphsSPR <- makeSPRGraphPar makeSPRGraphAction (fmap fst toReoptimizeAndCheckCostSPR)

                                -- process TBR returns (if done)
                                let toReoptimizeAndCheckCostTBR =   if (checkHeuristic swapParams == Better) then bettterHeuristicEdgesTBR
                                                                    else if (checkHeuristic swapParams == BestOnly) then filter ((== minHeuristicCost) . snd) bettterHeuristicEdgesTBR
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

                                let minimimCheckedCost = if (null checkedGraphCosts) then infinity
                                                         else minimum $ fmap snd5 checkedGraphCosts

                                (newBestGraphs, newBestCost) <- if minimimCheckedCost < curBestCost then do
                                                                    logWith LogInfo $ "\t-> " <> (show minimimCheckedCost)
                                                                    pure (filter ((== minimimCheckedCost) .snd5) checkedGraphCosts, minimimCheckedCost)
                                                                else 
                                                                    pure (curBestGraphList, curBestCost)

                                -- recurse or shortcut if better and steepest
                                if steepest swapParams && minimimCheckedCost < curBestCost then
                                    pure $  zip newBestGraphs (replicate (length newBestGraphs) newBestCost)         
                                else doAllSplitsAndRejoin swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor curBestGraphList curBestCost firstFullGraph restEdges
                        

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
    → PhyG ((LG.LEdge EdgeInfo, VertexCost), [((LG.LEdge EdgeInfo, LG.LEdge EdgeInfo), VertexCost)])
singleJoinHeuristic makeEdgeDataFunction edgeJoinFunction swapParams inGS inData splitGraph splitCost prunedGraphRootIndex tbrReRootEdgePairList targetEdge =
    let charInfoVV = fmap thd3 $ thd3 inData
    in do
            -- Create union-thype data for target edge
            targetEdgeData ← makeEdgeDataFunction splitGraph charInfoVV targetEdge

            -- get the data from the root of the pruned graph (or leaf if only that)
            let prunedRootVertexData = vertData $ fromJust $ LG.lab splitGraph prunedGraphRootIndex

            -- calculate heuristics join cost
            sprReJoinCost ← edgeJoinFunction charInfoVV prunedRootVertexData targetEdgeData

            
            -- SPR return new graph with heuristic cost
            if null tbrReRootEdgePairList then 
                pure ((targetEdge, splitCost + sprReJoinCost), []) 

            -- logWith LogInfo $ " TBR reroots: " <> (show $ length tbrReRootEdgePairList)

            -- tbr (includes SPR result)
            else 
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
                
                    pure ((targetEdge, splitCost + sprReJoinCost), tbrResult)    

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



