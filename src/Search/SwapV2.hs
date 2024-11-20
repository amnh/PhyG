{- |
Module specifying graph swapping rearrangement functions
-}
module Search.SwapV2 (
    rejoinGraphTuple,
    reoptimizeSplitGraphFromVertexIANew,
    reoptimizeSplitGraphFromVertexNew,
    reoptimizeSplitGraphFromVertexTupleFuse,
    reoptimizeSplitGraphFromVertexTupleNew,
    swapDriver,
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

{- | SwapDriver
    Top levl formtesting and swithcing between functions

    SA stuff should be changed in swapMaster to reflect the single SAParams
-}
swapDriver
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Int
    → [ReducedPhylogeneticGraph]
    → [(Maybe SAParams, ReducedPhylogeneticGraph)]
    → PhyG ([ReducedPhylogeneticGraph], Int)
swapDriver swapParams inGS inData inCounter curBestGraphList inSimAnnealParams = 
    let saList = L.uncons inSimAnnealParams
    in
    if isNothing saList then
    -- swapDriver' swapParams inGS inData inCounter curBestGraphList inSimAnnealParams 
        swapV2  swapParams inGS inData inCounter curBestGraphList Nothing 
    else 
        swapV2  swapParams inGS inData inCounter curBestGraphList ((fst . fst . fromJust) saList)

{- | New Swap functions that are based on strict PHANE parallelization routines.
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

                -- this for genetic algorithm to generate mutated graphs
                if returnMutated swapParams then
                    pure (concat $ annealGraphList,sum counterList)

                -- for regular swapping
                else do 
                    bestGraphs <- GO.selectGraphs Best (outgroupIndex inGS) (keepNum swapParams) 0.0 $ curBestGraphList <> (concat annealGraphList)

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
            -- this for genetic algorithm to generate mutated graphs
            if returnMutated swapParams then 
                pure (saGraphs, saCounter) 

            -- swap "back" to optimal costs
            else do 
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
                swapNaive swapParams inGS inData inCounter 0 curBestGraphList curBestGraphList saParams
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
            (sprGraphResult, sprCount) <- swapNaive (swapParams {swapType = SPR}) inGS inData inCounter 0 curBestGraphList curBestGraphList Nothing
            (tbrOnlyResult, _) <- swapNaive (swapParams {swapType = TBR}) inGS inData sprCount 0 sprGraphResult sprGraphResult Nothing
            
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

    SplitCounter to continue splits from point where previous ended in steepest (Farris suggestion)
-}
swapNaive 
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Int
    -> Int 
    → [ReducedPhylogeneticGraph]
    → [ReducedPhylogeneticGraph]
    → Maybe SAParams
    → PhyG ([ReducedPhylogeneticGraph], Int)
swapNaive swapParams inGS inData inCounter splitCounter graphsToSwap curBestGraphList saParams =
    let inGraphPair = L.uncons graphsToSwap
    in
    if isNothing inGraphPair
        then pure (curBestGraphList, inCounter)
        else do
            --logWith LogInfo $ " Split Counter " <> (show splitCounter)
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
                rejoinResult' <- if (splitParallel swapParams) && (isNothing saParams) then do
                                -- splitAction ::  LG.LEdge EdgeInfo → PhyG (DecoratedGraph, VertexCost, LG.Node, LG.Node, LG.Node)
                                    let splitAction = doASplit swapParams inGS inData (doIA swapParams) nonExactCharacters inGraphNetPenaltyFactor fullFirstGraph
                                    spltActionPar <- (getParallelChunkTraverseBy snd5)
                                    resultListP <- spltActionPar splitAction edgeList

                                    -- this filter for malformed graphs from split graph
                                    -- can happen with networks when there are lots of network edges
                                    let resultListP' = 
                                            let fiveList = filter (not . (LG.isEmpty . fst5)) resultListP
                                                (a,b,c,d,e) = L.unzip5 fiveList
                                                f = L.replicate (length fiveList) Nothing
                                            in L.zip6 a b c d e f 

                                    -- add Nothing for fuse edges


                                    -- Reorder list based on previous splits to keep going where previous stopped (Farris suggestion)
                                    let resultListP'' = if atRandom swapParams then resultListP'
                                                        else if isJust saParams then resultListP'
                                                        else 
                                                            -- order of split cost or not
                                                            let newEdgeOrder =  if sortEdgesSplitCost swapParams then 
                                                                                        L.sortOn snd6 resultListP'
                                                                                else resultListP'
                                                                splitNumber = rem splitCounter (length resultListP')
                                                            in
                                                            (drop splitNumber newEdgeOrder) <> (take splitNumber newEdgeOrder)

                                    rejoinFromOptSplitList swapParams inGS inData (doIA swapParams) inGraphNetPenaltyFactor graphsToSwap curBestCost splitCounter resultListP''

                                -- this is recursive on splitting as opposed to parallelized
                                -- used by SA/Drift (with randomization) as well to save work
                                else do
                                    -- Reorder list based on previous splits to keep going where previous stopped (Farris suggestion)
                                    let edgeList' = if atRandom swapParams then edgeList
                                                    else if isJust saParams then edgeList
                                                    else
                                                        let splitNumber = rem splitCounter (length edgeList)
                                                        in
                                                        (drop splitNumber edgeList) <> (take splitNumber edgeList)
                                    doAllSplitsAndRejoin swapParams inGS inData (doIA swapParams) nonExactCharacters inGraphNetPenaltyFactor graphsToSwap curBestCost splitCounter fullFirstGraph saParams edgeList'

                let (rejoinResult, splitCounter) = rejoinResult'

                -- Always return the whole lot for SA, since BestOnly has been specified should be single graph
                let (newValList, returnCost) = if isJust saParams then 
                                                 (fmap fst rejoinResult, minimum $ fmap snd rejoinResult)

                                                -- found better
                                                else if (minimum $ fmap snd rejoinResult) < curBestCost then 
                                                     (fmap fst rejoinResult, minimum $ fmap snd rejoinResult)

                                                -- found more of same cost
                                                else if (minimum $ fmap snd rejoinResult) == curBestCost then 
                                                     ((fmap fst rejoinResult) <> curBestGraphList, curBestCost)

                                                -- did not find any batter or equal
                                                else  (curBestGraphList, curBestCost)

                -- Conditions to keep going or return
                -- SA/Drift hit max numbers
                if isJust saParams then
                    pure (newValList, inCounter + 1)

                -- Better graphs than input
                else if returnCost < curBestCost then 
                    if (TBROnly == swapType swapParams) then
                        pure (newValList, inCounter + 1)
                    else
                        swapNaive swapParams inGS inData (inCounter + 1) splitCounter newValList newValList saParams

                -- equal graph costs to input
                else if returnCost == curBestCost then do
                    uniqueBestGraphs <- GO.selectGraphs Unique (outgroupIndex inGS) (keepNum swapParams) 0.0 $ newValList <> curBestGraphList
                    swapNaive swapParams inGS inData (inCounter + 1) splitCounter (graphsRemaining L.\\ uniqueBestGraphs) uniqueBestGraphs saParams

                -- worse graphs than input
                else swapNaive swapParams inGS inData (inCounter + 1) splitCounter graphsRemaining curBestGraphList saParams

{- | rejoinFromOptSplitList takes a list of optimized split graphs and
    calls rejoin function
    Good for parallel execution

    Does not have SA/Drift functionality

    last arg in split info is edges for use by fuse reconnection--first of which is original fuse 
    connection
-}
rejoinFromOptSplitList 
    :: SwapParams 
    -> GlobalSettings 
    → ProcessedData 
    -> Bool 
    -> VertexCost 
    -> [ReducedPhylogeneticGraph]
    -> VertexCost 
    -> Int
    -- -> [LG.LEdge EdgeInfo] 
    -> [(DecoratedGraph, VertexCost, LG.Node, LG.Node, LG.Node, Maybe [LG.LEdge EdgeInfo])] 
    -> PhyG ([(ReducedPhylogeneticGraph, VertexCost)], Int)
rejoinFromOptSplitList swapParams inGS inData doIA inGraphNetPenaltyFactor curBestGraphList curBestCost splitCounter splitInfoList'' =
    let splitInfoList' = L.uncons splitInfoList''
    in
    -- split counter back to 0 when end of list of splits
    if isNothing splitInfoList' then pure $ (zip curBestGraphList (fmap snd5 curBestGraphList), 0)
    else 
        let (firstSplit, restSplits) = fromJust splitInfoList'
            (splitGraphOptimized, splitCost, graphRoot, prunedGraphRootIndex, originalConnectionOfPruned, fuseEdgesToJoin) = firstSplit

            -- get root in base (for readdition) and edges in pruned section for rerooting during readdition
            (_, edgesInPrunedGraph) = LG.nodesAndEdgesAfter splitGraphOptimized [(originalConnectionOfPruned, fromJust $ LG.lab splitGraphOptimized originalConnectionOfPruned)]

            -- Filter for bridge edges and desendents of root of pruned graph -- for TBR when needed
            -- remeber to add in root of pruned graph data for TBR edge data so SPR is subset TBR
            tbrRerootEdges
                | swapType swapParams `elem` [NNI, SPR] = []
                | (graphType inGS == Tree) || LG.isTree splitGraphOptimized = edgesInPrunedGraph L.\\ ((LG.out splitGraphOptimized originalConnectionOfPruned) <> (LG.out splitGraphOptimized prunedGraphRootIndex))
                | otherwise = (fmap fst . filter snd . zip edgesInPrunedGraph $ LG.isBridge splitGraphOptimized . LG.toEdge <$> edgesInPrunedGraph) L.\\ ((LG.out splitGraphOptimized originalConnectionOfPruned) <> (LG.out splitGraphOptimized prunedGraphRootIndex))

            charInfoVV = fmap thd3 $ thd3 inData

            -- check for fuse edges input.  Only use if NoSwap which contains the initial fuse edge
            -- take 3 for NoSwap to do a litle more work
            (_, edgesInBaseGraph) = if isNothing fuseEdgesToJoin then                                
                                        LG.nodesAndEdgesAfter splitGraphOptimized [(graphRoot, fromJust $ LG.lab splitGraphOptimized graphRoot)]
                                    else ([], fromJust fuseEdgesToJoin)

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
            -- logWith LogInfo $ " Split Counter (Opt) " <> (show splitCounter)
            
            if (null tbrRerootEdges) && (swapType swapParams == TBROnly) then
                    rejoinFromOptSplitList swapParams inGS inData doIA inGraphNetPenaltyFactor curBestGraphList curBestCost (splitCounter + 1)  restSplits

            -- empty edges can happen in networks due to non-briding edge deficit
            else if null edgesInBaseGraph || null edgesInPrunedGraph then 
                    rejoinFromOptSplitList swapParams inGS inData doIA inGraphNetPenaltyFactor curBestGraphList curBestCost (splitCounter + 1)  restSplits

            else do
                {-  pruning of edge via unions
                    Does not seem to help as much as I would have thought--perhaps the union process is not coreect
                    Leaving in to perhaps use/fix later or perhaps really only significant for large taxon sets

                    May not work properly for networks yielding Nothing on edge label call--hence the check
                -}
                edgesInBaseGraph' <- -- edges fuse 
                                     if (isJust fuseEdgesToJoin) && (swapType swapParams == NoSwap) then 
                                        -- for fuse--first is fuse connection w/o swap
                                        pure $ take 3 edgesInBaseGraph

                                     -- network
                                     else if graphType inGS /= Tree then 
                                        pure edgesInBaseGraph
                                     
                                     else if joinType swapParams == JoinPruned && isJust (LG.lab splitGraphOptimized prunedGraphRootIndex) then 
                                        do 
                                            let prunedToRejoinUnionData = vertData $ fromJust $ LG.lab splitGraphOptimized prunedGraphRootIndex
                                            getUnionRejoinEdgeListNew
                                                inGS
                                                splitGraphOptimized
                                                charInfoVV
                                                [graphRoot]
                                                (curBestCost - splitCost)
                                                prunedToRejoinUnionData
                                                []
                                                             
                                     else pure edgesInBaseGraph
                

                let maxMoveEdgeDistance = min (maxMoveEdgeDist swapParams) (maxBound ∷ Int)

                -- reorder/shuffle edge list if desired
                rejoinEdges <-  if isJust fuseEdgesToJoin then 
                                    pure edgesInBaseGraph'

                                else if atRandom swapParams then 
                                    shuffleList edgesInBaseGraph'
                                    
                                -- should re-add close to original placement first
                                else if swapType swapParams == NNI then
                                    pure $ take 3 $ LG.sortEdgesByIndexDistance originalConnectionOfPruned edgesInBaseGraph'
                                                                                    
                                else 
                                    pure $ take maxMoveEdgeDistance $ LG.sortEdgesByIndexDistance originalConnectionOfPruned edgesInBaseGraph'


                
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
                    pure $ (zip newBestGraphs (replicate (length newBestGraphs) newBestCost) , splitCounter + 1)

                else rejoinFromOptSplitList swapParams inGS inData doIA inGraphNetPenaltyFactor newBestGraphs newBestCost (splitCounter + 1)  restSplits


{- | rejoinGraphTuple is a wrapper around rejoinFromOptSplitList for use
in fusing
-}
rejoinGraphTuple 
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → VertexCost
    → [ReducedPhylogeneticGraph]
    → Maybe SAParams
    → (DecoratedGraph, SimpleGraph, VertexCost, LG.Node, LG.Node, LG.Node, [LG.LEdge EdgeInfo], [LG.LEdge EdgeInfo], VertexCost)
    → PhyG [ReducedPhylogeneticGraph]
rejoinGraphTuple swapParams inGS inData curBestCost curBestGraphs inSimAnnealParams
    ( reoptimizedSplitGraph
        , splitGraphSimple
        , splitGraphCost
        , graphRoot
        , prunedGraphRootIndex
        , originalConnectionOfPruned
        , rejoinEdges -- first element is oringal connection--edge that can be rejoined to
        , edgesInPrunedGraph
        , netPenaltyFactor
        ) = 
    let splitInfo = (reoptimizedSplitGraph, splitGraphCost,  graphRoot, prunedGraphRootIndex, originalConnectionOfPruned, Just rejoinEdges)
    in do
        rejoinedReturnList <- rejoinFromOptSplitList swapParams {joinType = JoinAll, steepest = False} inGS inData False netPenaltyFactor curBestGraphs curBestCost 0 [splitInfo]
        let rejoinGraphsList = fmap fst $ fst rejoinedReturnList
        graphsToReturn <- GO.selectGraphs Unique (outgroupIndex inGS) (keepNum swapParams) 0.0 $ curBestGraphs <> rejoinGraphsList
        pure graphsToReturn

{- | doAllSplitsAndRejoin generates all split reoptimizations
    Calls the rejoin function after creating the optimizaed split graph
    recurses to next edge to split
    Good for space conservation

    Used by SA/Drift
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
    -> Int
    -> PhylogeneticGraph 
    -> Maybe SAParams
    -> [LG.LEdge EdgeInfo] 
    -> PhyG ([(ReducedPhylogeneticGraph, VertexCost)], Int)
doAllSplitsAndRejoin swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor curBestGraphList curBestCost splitCounter firstFullGraph saParams edgeList'' = 
                let edgeList' = L.uncons edgeList''
                in
                -- split counter back to 0 when complete splits list
                if isNothing edgeList' then pure $ (zip curBestGraphList (fmap snd5 curBestGraphList), 0)
                else 
                    let edgeList@(firstEdge, restEdges) = fromJust edgeList'
                    in do
                        -- logWith LogInfo $ " Split Counter (Edge) " <> (show splitCounter)
                        -- split graph on the first edge
                        let (splitGraph, graphRoot, prunedGraphRootIndex, originalConnectionOfPruned) = LG.splitGraphOnEdge (thd6 firstFullGraph) firstEdge
                        
                        (splitGraphOptimized, splitCost) ←
                                reoptimizeSplitGraphFromVertexNew swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor firstFullGraph splitGraph graphRoot prunedGraphRootIndex 
                        
                        -- this if no split graph--can happen if lots of network edges
                        if LG.isEmpty splitGraphOptimized then
                            doAllSplitsAndRejoin swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor curBestGraphList curBestCost (splitCounter + 1) firstFullGraph saParams restEdges

                        -- avoid shorcircuiting on SA/Drift
                        else if (splitCost >= curBestCost) && (isNothing saParams) then 
                            doAllSplitsAndRejoin swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor curBestGraphList curBestCost (splitCounter + 1) firstFullGraph saParams restEdges
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


                            in 
                            -- empty edges can happen in networks due to non-briding edge deficit
                            if null edgesInBaseGraph || null edgesInPrunedGraph then 
                                doAllSplitsAndRejoin swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor curBestGraphList curBestCost (splitCounter + 1) firstFullGraph saParams restEdges

                            else do

                                {-pruning of edge via unions
                                    Does not seem to help as much as I would have thought--perhaps the union process is not coreect
                                    Leaving in to perhaps use/fix later or perhaps really only significant for large taxon sets

                                    May not work properly for networks yielding Nothing on edge label call--hence the check
                                -}
                                edgesInBaseGraph' <- if graphType inGS /= Tree then pure edgesInBaseGraph
                                                     else if joinType swapParams == JoinPruned && isJust (LG.lab splitGraphOptimized prunedGraphRootIndex) then 
                                                        do 
                                                            let prunedToRejoinUnionData = vertData $ fromJust $ LG.lab splitGraphOptimized prunedGraphRootIndex
                                                            getUnionRejoinEdgeListNew
                                                                inGS
                                                                splitGraphOptimized
                                                                charInfoVV
                                                                [graphRoot]
                                                                ((snd6 firstFullGraph) - splitCost)
                                                                prunedToRejoinUnionData
                                                                []
                                                             
                                                     else pure edgesInBaseGraph

                                
                                let maxMoveEdgeDistance = min (maxMoveEdgeDist swapParams) (maxBound ∷ Int)

                                rejoinEdges <- if atRandom swapParams then 
                                                 shuffleList edgesInBaseGraph'
                                               -- should re-add close to original placement first
                                               else if swapType swapParams == NNI then
                                                    pure $ take 3 $ LG.sortEdgesByIndexDistance originalConnectionOfPruned edgesInBaseGraph'
                                                                        
                                               else 
                                                    pure $ take maxMoveEdgeDistance $ LG.sortEdgesByIndexDistance originalConnectionOfPruned edgesInBaseGraph'
                                
                                -- logWith LogInfo $ "\nEdge to rejoin: " <> (show $ length rejoinEdges) <> " Union: " <> (show $ length unionEdgeList)

                                -- shot circuit for TBTROnly in Alternate
                                if (null tbrRerootEdges) && (swapType swapParams == TBROnly) then
                                    doAllSplitsAndRejoin swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor curBestGraphList curBestCost (splitCounter + 1) firstFullGraph saParams restEdges

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
                                            pure $ (zip newBestGraphs (replicate (length newBestGraphs) newBestCost), splitCounter + 1)
                                        else
                                            doAllSplitsAndRejoin swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor newBestGraphs curBestCost (splitCounter + 1) firstFullGraph newSAParams restEdges


                                    -- regular swap
                                    else do
                                        (newBestGraphs, newBestCost) <- if minimumCheckedCost < curBestCost then do
                                                                            logWith LogInfo $ "\t-> " <> (show minimumCheckedCost)
                                                                            pure (filter ((== minimumCheckedCost) .snd5) checkedGraphCosts, minimumCheckedCost)
                                                                        else 
                                                                            pure (curBestGraphList, curBestCost)

                                        if ((swapType swapParams == TBROnly) || steepest swapParams) && minimumCheckedCost < curBestCost then
                                            pure $  ((zip newBestGraphs (replicate (length newBestGraphs) newBestCost)), splitCounter + 1)
                                        
                                        else doAllSplitsAndRejoin swapParams inGS inData doIA nonExactCharacters inGraphNetPenaltyFactor curBestGraphList curBestCost (splitCounter + 1) firstFullGraph saParams restEdges


{- | getUnionRejoinEdgeListNew takes a graph (split and reoptimized usually), the overall root index (of split),
split cost, and union threshold value and returns list of edges that have union distance <= threshold factor
assumes that root edges are not network edges (an invariant)
checks node then recurses to children
-}
getUnionRejoinEdgeListNew
    ∷ GlobalSettings
    → DecoratedGraph
    → V.Vector (V.Vector CharInfo)
    → [LG.Node]
    → Double
    → VertexBlockData
    → [LG.LEdge EdgeInfo]
    → PhyG [LG.LEdge EdgeInfo]
getUnionRejoinEdgeListNew inGS inGraph charInfoVV nodeIndexList splitDiffCost nodeToJoinUnionData curEdgeList =
    -- returns edges in post order since prepending--could reverse if want pre-order edges
    -- might be better post order, unclear
    if null nodeIndexList
        then pure curEdgeList
        else
            let nodeIndex = head nodeIndexList
                childEdges = LG.out inGraph nodeIndex
                childNodeIndexList = fmap snd3 childEdges

                -- node data and union distance
                nodeData = vertData $ fromJust $ LG.lab inGraph nodeIndex
            in  -- nodeDataString = U.getUnionFieldsNode nodeData
                -- toJoinString = U.getUnionFieldsNode nodeToJoinUnionData

                -- traceNoLF ("GURE: " <> (show nodeIndex)) (
                -- trace ("GUREL: " <> (show (unionDistance, splitDiffCost, unionThreshold * splitDiffCost))) (
                if (length childEdges) `notElem` [0, 1, 2]
                    then error ("Node has improper number of children : " <> (show $ length childEdges))
                    else do
                        -- add non-outgroup root edge to list for rejoin after checking for acceptable union distance
                        -- if edge is not within union distance factor then stop--no recursion
                        -- this since distance cannot get lower further upt the graph given union creation

                        unionDistance ← getUnionDistanceM nodeToJoinUnionData nodeData charInfoVV
                        let metThreshold = unionDistance < splitDiffCost * (unionThreshold inGS)

                        if LG.isRoot inGraph nodeIndex
                            then
                                if metThreshold
                                    then
                                        if (null $ LG.out inGraph (head childNodeIndexList))
                                            then
                                                getUnionRejoinEdgeListNew
                                                    inGS
                                                    inGraph
                                                    charInfoVV
                                                    [(snd3 $ last childEdges)]
                                                    splitDiffCost
                                                    nodeToJoinUnionData
                                                    ((last childEdges) : curEdgeList)
                                            else
                                                getUnionRejoinEdgeListNew
                                                    inGS
                                                    inGraph
                                                    charInfoVV
                                                    [(snd3 $ head childEdges)]
                                                    splitDiffCost
                                                    nodeToJoinUnionData
                                                    ((head childEdges) : curEdgeList)
                                    else pure curEdgeList
                            else do
                                -- non-root node--process childre 1/2
                                -- recurses to their children if union condition met--but doesn't add network edges
                                -- check current node--then recurse to children

                                -- let -- first and second child data child
                                newCurEdgeListChild ←
                                    if metThreshold
                                        then getUnionRejoinEdgeListNew inGS inGraph charInfoVV childNodeIndexList splitDiffCost nodeToJoinUnionData (childEdges <> curEdgeList)
                                        else pure curEdgeList
                                -- in  -- recurse remaining nodes
                                getUnionRejoinEdgeListNew inGS inGraph charInfoVV (tail nodeIndexList) splitDiffCost nodeToJoinUnionData newCurEdgeListChild


{- | getUnionDistanceM gets distance between two the union fields of two characters
since its a distance no need for no change cost adjustment
-}
getUnionDistanceM ∷ VertexBlockData → VertexBlockData → V.Vector (V.Vector CharInfo) → PhyG Double
getUnionDistanceM union1 union2 charInfoVV =
    let noChangeCostAdjut = False
    in  M.distance2UnionsM noChangeCostAdjut union1 union2 charInfoVV


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
                -- NoSwap for fusing
                if swapType swapParams `elem` [NoSwap, NNI, SPR] then 
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

-- | reoptimizeSplitGraphFromVertexTupleNew wrapper for reoptimizeSplitGraphFromVertex with last 4 args as tuple
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




-- | reoptimizeSplitGraphFromVertexTupleFuse wrapper for reoptimizeSplitGraphFromVertexFuse with last 3 args as tuple
reoptimizeSplitGraphFromVertexTupleFuse
    ∷ GlobalSettings
    → ProcessedData
    → Bool
    → VertexCost
    → (DecoratedGraph, Int, Int)
    → PhyG (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertexTupleFuse inGS inData doIA netPenaltyFactor (inSplitGraph, startVertex, prunedSubGraphRootVertex) =
    reoptimizeSplitGraphFromVertexFuse inGS inData doIA netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex

{- The functions below do more work (hence less efficent--no incremetal e.g.) but are used by fuse when swapping is added.
    They are used (and here) because the fusing of two graphs has unexpected vertex indices which are difficult to track
    (or redo without O(n) cost).

-}


{- | reoptimizeSplitGraphFromVertex fully labels the component graph that is connected to the specified vertex
retuning that graph with 2 optimized components and their cost
both components goo through multi-traversal optimizations
doIA option to only do IA optimization as opposed to full thing--should be enormously faster--but yet more approximate
creates final for both base graph and priunned component due to rerooting non-concordance of preorder and post order assignments
terminology bse graph is the component with the original root, pruned that which has been removed form the original
graph to be readded to edge set
The function
1) optimizes two components seprately from their "root"
2) takes nodes and edges for each and cretes new graph
3) returns graph and summed cost of two components
4) adds in root and netPenalty factor estimates since net penalty can only be calculated on full graph
part of this is turning off net penalty cost when optimizing base and pruned graph components
if doIA is TRUE then call function that onl;y optimizes the IA assignments on the "original graph" after split.
this keeps teh IA chracters in sync across the two graphs
NB uses PhylogeneticGraph internally
This should return infinity for split graph cost if either component is emptyGraph

Tested 3 doublings of metazoa and roughtly O(n)
-}
reoptimizeSplitGraphFromVertexFuse
    ∷ GlobalSettings
    → ProcessedData
    → Bool
    → VertexCost
    → DecoratedGraph
    → Int
    → Int
    → PhyG (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertexFuse inGS inData doIA netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex =
    -- trace ("RSGFV: " <> (show startVertex)) (
    if doIA
        then -- only reoptimize the IA states for dynamic characters
            reoptimizeSplitGraphFromVertexIA inGS inData netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex
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




{- | reoptimizeSplitGraphFromVertexIA performs operations of reoptimizeSplitGraphFromVertex for static charcaters
but dynamic characters--only update IA assignments and initialized from origPhylo graph (at leaves) to keep IA characters in sync
since all "static" only need single traversal post order pass

uses PhylogenetiGraph internally
-}
reoptimizeSplitGraphFromVertexIA
    ∷ GlobalSettings
    → ProcessedData
    → VertexCost
    → DecoratedGraph
    → Int
    → Int
    → PhyG (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertexIA inGS inData netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex =
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