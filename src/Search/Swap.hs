{- |
Module specifying graph swapping rearrangement functions
-}
module Search.Swap (
    --getUnionRejoinEdgeList,
    rejoinGraphTuple,
    reoptimizeSplitGraphFromVertexTuple,
    swapDriver,
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
import Search.SwapV2 qualified as SV2
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
        SV2.swapV2  swapParams inGS inData inCounter curBestGraphList Nothing 
    else 
        SV2.swapV2  swapParams inGS inData inCounter curBestGraphList ((fst . fst . fromJust) saList)

{- | SwapDriver
    Uses compnent functions but with alternate high-level logic
    Generates new simmanneal params for recursive rounds
-}
swapDriver'
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Int
    → [ReducedPhylogeneticGraph]
    → [(Maybe SAParams, ReducedPhylogeneticGraph)]
    → PhyG ([ReducedPhylogeneticGraph], Int)
swapDriver' swapParams inGS inData inCounter curBestGraphList inSimAnnealParams = 
    if null curBestGraphList
        then pure ([], inCounter)
        else do
            let inBestCost = minimum $ fmap snd5 curBestGraphList
            -- liftIO $ putStr ("SD: " <> (show $ swapType swapParams)) --  <> " " <> (show $ fmap snd5 curBestGraphList)) -- <> "\n" <> (LG.prettyIndices $ head $ fmap fst5 curBestGraphList))

            (newGraphList', newCounter, newSAParams) ←
                swapMaster swapParams inGS inData inCounter curBestGraphList inSimAnnealParams

            newGraphList ← GO.selectGraphs Best (outgroupIndex inGS) (keepNum swapParams) 0.0 newGraphList'
            -- liftIO $ putStr ("SD-After: " <> (show $ fmap snd5 newGraphList))
            -- found no better

            if null newGraphList
                then pure (curBestGraphList, newCounter)
                else
                    let newCost = minimum $ fmap snd5 newGraphList
                    in  -- found worse for some reason
                        if newCost > inBestCost
                            then pure (curBestGraphList, newCounter)
                            else -- found same (ie additional)

                                if newCost == inBestCost
                                    then do
                                        allGraphs <- GO.selectGraphs Best (outgroupIndex inGS) (keepNum swapParams) 0.0 (newGraphList <> curBestGraphList)
                                        pure (allGraphs, newCounter)
                                    else -- found better-- go around again
                                    do
                                        let newSimAnnealParamList = zip (U.generateUniqueRandList (length newGraphList) newSAParams) newGraphList
                                        swapDriver' swapParams inGS inData newCounter newGraphList newSimAnnealParamList


{- swapMAster
    Working through the complex logic of swapSPRTBR
-}
swapMaster
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Int
    → [ReducedPhylogeneticGraph]
    → [(Maybe SAParams, ReducedPhylogeneticGraph)]
    → PhyG ([ReducedPhylogeneticGraph], Int, Maybe SAParams)
swapMaster swapParams inGS inData@(leafNames, _, _) inCounter curBestGraphList inSimAnnealParams =
    -- swapSPRTBR swapParams inGS inData inCounter curBestGraphList inSimAnnealParams
    let curBestCost = minimum $ fmap snd5 curBestGraphList
        numLeaves = V.length leafNames
    in  do
            inGraphNetPenalty ← T.getPenaltyFactor inGS inData Nothing $ GO.convertReduced2PhylogeneticGraph (head curBestGraphList)
            let inGraphNetPenaltyFactor = inGraphNetPenalty / curBestCost

            -- (a,b,c) <- swapAll swapParams inGS inData inCounter curBestCost curBestGraphList curBestGraphList numLeaves netPenaltyFactor (fst $ head inSimAnnealParams)
            -- a' <- GO.selectGraphs Best (keepNum swapParams) 0.0 a
            -- pure (a,b,c)

            -- THis seems to screw up
            case (fst $ head inSimAnnealParams) of
                Nothing → do
                    -- steepest takes immediate best--does not keep equall cost-- for now--disabled not working correctly so goes to "all"
                    (swappedGraphs, counter, swapSAPArams) ←
                        swapAll
                            swapParams
                            inGS
                            inData
                            inCounter
                            curBestCost
                            curBestGraphList
                            curBestGraphList
                            numLeaves
                            inGraphNetPenaltyFactor
                            Nothing
                    pure $ case swappedGraphs of
                        [] → (curBestGraphList, counter, Nothing)
                        gs → (gs, counter, swapSAPArams)
                -- simulated annealing/drifting acceptance does a steepest with SA acceptance
                Just simAnneal →
                    -- then a swap steepest and all on annealed graph
                    -- same at this level method (SA, Drift) choice occurs at lower level
                    -- annealed should only yield a single graph
                    let -- create list of params with unique list of random values for rounds of annealing
                        annealDriftRounds = rounds simAnneal
                        newSimAnnealParamList = U.generateUniqueRandList annealDriftRounds (fst $ head inSimAnnealParams)
                        -- parallel setup
                        action ∷ Maybe SAParams → PhyG ([ReducedPhylogeneticGraph], Int, Maybe SAParams)
                        action = swapAll swapParams inGS inData 0 curBestCost curBestGraphList curBestGraphList numLeaves inGraphNetPenaltyFactor

                        extractGraphTopoCost
                            ∷ ([ReducedPhylogeneticGraph], Int, Maybe SAParams)
                            → ([ReducedPhylogeneticGraph], Int, Maybe SAParams)
                        extractGraphTopoCost = applyOver1of3 (listApplying strict1and2of5)
                    in  -- this to ensure current step set to 0
                        do
                            swapPar ← getParallelChunkTraverseBy extractGraphTopoCost
                            (annealDriftGraphs', anealDriftCounterList, annealDriftParams) ← unzip3 <$> swapPar action newSimAnnealParamList

                            -- annealed/Drifted 'mutated' graphs
                            annealDriftGraphs ← GO.selectGraphs Unique (outgroupIndex inGS) (keepNum swapParams) 0.0 $ concat annealDriftGraphs'

                            -- swap back "normally" if desired for full drifting/annealing
                            (swappedGraphs, counter, _) ←
                                swapAll
                                    swapParams
                                    inGS
                                    inData
                                    (sum anealDriftCounterList)
                                    (min curBestCost (minimum $ snd5 <$> annealDriftGraphs))
                                    curBestGraphList
                                    annealDriftGraphs
                                    numLeaves
                                    inGraphNetPenaltyFactor
                                    Nothing

                            bestGraphs ← GO.selectGraphs Best (outgroupIndex inGS) (keepNum swapParams) 0.0 $ curBestGraphList <> swappedGraphs
                            -- this Bool for Genetic Algorithm mutation step
                            pure $
                                if not $ returnMutated swapParams
                                    then (bestGraphs, counter, head annealDriftParams)
                                    else (annealDriftGraphs, sum anealDriftCounterList, head annealDriftParams)


{- | swapAll is a high level function that basically deals with portioning out swap-type swaps
and performs the high level options for Alternate where SPR is perfomred first, then TBR,
but whenever a better (or additional) graph is found during TBR, an SPR swap of that graph
is performed before returning to TBR again.  THis contibues until no new graphs are found in the
SPR + TBR swap.
each call to swapAll' sets break edge number to 0
this swaps untill none found better
-}
swapAll
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Int
    → VertexCost
    → [ReducedPhylogeneticGraph]
    → [ReducedPhylogeneticGraph]
    → Int
    → VertexCost
    → Maybe SAParams
    → PhyG ([ReducedPhylogeneticGraph], Int, Maybe SAParams)
swapAll swapParams inGS inData counter curBestCost curSameBetterList [] numLeaves netPenaltyFactor inSimAnnealParams = pure (curSameBetterList, counter, inSimAnnealParams)
swapAll swapParams inGS inData counter curBestCost curSameBetterList inGraphList@(firstGraph : otherGraphs) numLeaves netPenaltyFactor inSimAnnealParams =
    {-# SCC swapAll_TOP_DEF #-}
    -- nni, spr, tbr
    case swapType swapParams of
        Alternate → do
            -- logWith LogInfo "In swapAll-Alt"
            (sprGraphs, sprCounter, sprSAPArams) ←
                swapAll'
                    (swapParams{swapType = SPR})
                    inGS
                    inData
                    counter
                    curBestCost
                    curSameBetterList
                    inGraphList
                    numLeaves
                    netPenaltyFactor
                    0
                    inSimAnnealParams
            graphsToTBR ← GO.selectGraphs Best (outgroupIndex inGS) (keepNum swapParams) 0.0 $ sprGraphs <> inGraphList
            let sprBestCost = getGraphCost graphsToTBR

            -- tbr until find better or novel equal
            (tbrGraphs, tbrCounter, tbrSAPArams) ←
                swapAll'
                    (swapParams{swapType = TBRAlternate})
                    inGS
                    inData
                    sprCounter
                    sprBestCost
                    graphsToTBR
                    graphsToTBR
                    numLeaves
                    netPenaltyFactor
                    0
                    sprSAPArams
            let tbrBestCost = case tbrGraphs of
                    [] → infinity
                    (_, c, _, _, _) : _ → c

            {- This isn't improving performance so turned off in SwapSPRTBR-}
            -- if found better and alternating union pruning then return so can go back to start union pruning again
            case joinType swapParams of
                JoinPruned | sprBestCost < curBestCost → pure (sprGraphs, sprCounter, sprSAPArams)
                JoinPruned | tbrBestCost < curBestCost → pure (tbrGraphs, tbrCounter, tbrSAPArams)
                -- JoinAlternate | sprBestCost < curBestCost → pure (sprGraphs, sprCounter, sprSAPArams)
                -- JoinAlternate | tbrBestCost < curBestCost → pure (tbrGraphs, tbrCounter, tbrSAPArams)
                _ → case tbrBestCost `compare` sprBestCost of
                    -- found nothing better or equal
                    GT → pure (graphsToTBR, tbrCounter, tbrSAPArams)
                    -- if TBR found better go around again with SPR first--since returned if found better during TBR rejoin
                    LT →
                        swapAll
                            swapParams
                            inGS
                            inData
                            tbrCounter
                            tbrBestCost
                            tbrGraphs
                            tbrGraphs
                            numLeaves
                            netPenaltyFactor
                            tbrSAPArams
                    -- check if found additional
                    EQ → do
                        bestTBRGraphs ← GO.selectGraphs Best (outgroupIndex inGS) (keepNum swapParams) 0.0 tbrGraphs
                        let newTBRGraphs = GO.reducedphylogeneticGraphListMinus bestTBRGraphs $ sprGraphs <> curSameBetterList
                        case newTBRGraphs of
                            g : gs
                                | length bestTBRGraphs < keepNum swapParams →
                                    swapAll
                                        swapParams
                                        inGS
                                        inData
                                        tbrCounter
                                        tbrBestCost
                                        tbrGraphs
                                        newTBRGraphs
                                        numLeaves
                                        netPenaltyFactor
                                        tbrSAPArams
                            -- found nothing new
                            _ → pure (tbrGraphs, tbrCounter, tbrSAPArams)

        -- 0 is for list of edges so move through list past stable edges
        _ →
            swapAll'
                swapParams
                inGS
                inData
                counter
                curBestCost
                curSameBetterList
                inGraphList
                numLeaves
                netPenaltyFactor
                0
                inSimAnnealParams


{- | swapAll' performs branch swapping on all 'break' edges and all readditions
this not a "map" version to reduce memory footprint to a more mangeable level
"steepest" a passable option to short circuit readdition action to return immediately
if better graph found
1) takes first graph
2) if steepeast checks to make sure <= current best cost
3) gets list of "break-able" edges
   all non-root edges if tree
   all non-root bridge edges if network
4) send list (and other info) to split-join function
   goes on if empty list returned or > current best
   add graphs todo list if == current best cost
5) returns all of minimum cost found
if Alternate then when found better do SPR first then TBR
assumes SPR done before  Alternate entering so can star with TBR and iff get better
go back to SPR. NBest for "steepest" descent
For drift and anneal need to randomize order of splits and rejoins
-}
swapAll'
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Int
    → VertexCost
    → [ReducedPhylogeneticGraph]
    → [ReducedPhylogeneticGraph]
    → Int
    → VertexCost
    → Int
    → Maybe SAParams
    → PhyG ([ReducedPhylogeneticGraph], Int, Maybe SAParams)
swapAll' swapParams inGS inData counter curBestCost curSameBetterList inGraphList numLeaves netPenaltyFactor breakEdgeNumber inSimAnnealParams =
    let selectionOf x = GO.selectGraphs x (outgroupIndex inGS) (keepNum swapParams) 0.0
    in  -- don't need to check for mutated here since checked above
        case inGraphList of
            [] →
                let tag x = (x, counter, inSimAnnealParams)
                in  tag <$> Unique `selectionOf` curSameBetterList
            firstGraph : tailGraphs
                | LG.isEmpty $ thd5 firstGraph →
                    swapAll'
                        swapParams
                        inGS
                        inData
                        counter
                        curBestCost
                        curSameBetterList
                        tailGraphs
                        numLeaves
                        netPenaltyFactor
                        breakEdgeNumber
                        inSimAnnealParams
            firstGraph : tailGraphs →
                let firstDecoratedGraph = thd5 firstGraph
                    (firstRootIndex, _) = head $ LG.getRoots firstDecoratedGraph

                    -- determine edges to break on--'bridge' edges only for network
                    -- filter out edges from root since no use--would just rejoin
                    -- sort longest edge to shortest--option to speeed up steepest and conditions for all as well
                    -- this edge sort from Varon and Wheeler 2013
                    breakEdgeList' =
                        let filtration = filter ((/= firstRootIndex) . fst3)
                            extractor
                                | graphType inGS == Tree || LG.isTree firstDecoratedGraph = LG.labEdges
                                | otherwise = LG.getEdgeSplitList
                            sorting
                                | not $ atRandom swapParams = GO.sortEdgeListByLength
                                | otherwise = id
                        in  sorting . filtration $ extractor firstDecoratedGraph
                in  do
                        -- logWith LogInfo "In swapAll'"
                        -- randomize edges list order for anneal and drift
                        breakEdgeList'' ←
                            if atRandom swapParams
                                then shuffleList breakEdgeList'
                                else pure breakEdgeList'

                        -- move first "breakEdgeFactor" edges in split list to end
                        -- since breakEdgeFactor can get incremented past number of edges the integer remainder is determined
                        -- to move to end
                        -- this to reduces the revisiting of stable edges (by moving them to the end of the list)
                        -- yet still insures that all edges will be visited in final (or ay time needed) split.
                        -- used in POY v 1-3, Came from Steve Farris pers. com.
                        let breakEdgeFactor = snd $ divMod breakEdgeNumber (length breakEdgeList'')
                        let breakEdgeList =
                                let (prefix, suffix) = splitAt breakEdgeFactor breakEdgeList''
                                in  suffix <> prefix

                        -- perform intial split and rejoin on each edge in first graph
                        splitJoinResult ←
                            {-# SCC splitJoinResult #-}
                                splitJoinGraph
                                    swapParams
                                    inGS
                                    inData
                                    curBestCost
                                    curSameBetterList
                                    numLeaves
                                    netPenaltyFactor
                                    inSimAnnealParams
                                    firstGraph
                                    breakEdgeNumber
                                    breakEdgeList
                                    breakEdgeList
                        let (newGraphList', newSAParams, newBreakEdgeNumber) = splitJoinResult

                        -- get best return graph list-can be empty if nothing better ort smame cost
                        newGraphList ← Best `selectionOf` newGraphList'

                        -- get unique return graph list-can be empty if nothing better ort same cost
                        let newMinCost = getGraphCost `orInfinity` newGraphList

                        -- logic for returning normal swap operations (better only)
                        -- versus simulated annealin/Drifing returning potentially sub-optimal
                        case inSimAnnealParams of
                            Nothing → case newMinCost `compare` curBestCost of
                                -- found better cost graph
                                LT → do
                                    logWith LogInfo $ "\t->" <> (show newMinCost)
                                    -- for alternate do SPR first then TBR
                                    -- for alternate in TBR or prune union alternate if found better return immediately
                                    if swapType swapParams == TBRAlternate || joinType swapParams == JoinPruned || steepest swapParams
                                        then pure (newGraphList, counter, newSAParams)
                                        else -- regular swap--keep going with better graphs

                                            swapAll'
                                                swapParams
                                                inGS
                                                inData
                                                (counter + 1)
                                                newMinCost
                                                newGraphList
                                                newGraphList
                                                numLeaves
                                                netPenaltyFactor
                                                newBreakEdgeNumber
                                                newSAParams

                                -- found only worse graphs--never happens due to the way splitjoin returns only better or equal
                                -- but could change
                                GT →
                                    swapAll'
                                        swapParams
                                        inGS
                                        inData
                                        (counter + 1)
                                        curBestCost
                                        curSameBetterList
                                        (tail inGraphList)
                                        numLeaves
                                        netPenaltyFactor
                                        0 -- breakEdgeNumber set to zero for new graph to look at
                                        newSAParams
                                -- found same cost graphs
                                EQ →
                                    -- Important to not limit curSameBest since may rediscover graphs via swapping on equal when limiting the number to keep
                                    -- can be a cause of infinite running issues.
                                    -- let newCurSameBetterList = GO.selectGraphs Best (keepNum swapParams) 0.0 (-1) (curSameBetterList <> newGraphList)
                                    let graphChoices = (tailGraphs <> newGraphList) `GO.reducedphylogeneticGraphListMinus` curSameBetterList
                                    in  do
                                            newCurSameBetterList ← (Best `selectionOf`) $ curSameBetterList <> newGraphList
                                            graphsToDo ← Best `selectionOf` graphChoices

                                            -- these conditions help to prevent recswapping endlessly on new graphs thatare not in buffers,
                                            -- but have same cost
                                            let graphsToDo'
                                                    | keepNum swapParams - 1 <= length graphsToDo = tailGraphs
                                                    -- found nothing better that is new
                                                    | keepNum swapParams - 1 <= length newCurSameBetterList = tailGraphs
                                                    | otherwise = graphsToDo

                                            swapAll'
                                                swapParams
                                                inGS
                                                inData
                                                (counter + 1)
                                                curBestCost
                                                newCurSameBetterList
                                                graphsToDo'
                                                numLeaves
                                                netPenaltyFactor
                                                newBreakEdgeNumber
                                                newSAParams

                            -- simulated annealing/Drift post processing
                            Just simAnneal → case newMinCost `compare` curBestCost of
                                -- found better cost graph
                                LT → logWith LogInfo ("\t->" <> show newMinCost) $> (newGraphList, counter, newSAParams)
                                -- not better so check for drift changes or annealing steps and return if reached maximum number
                                _
                                    | currentStep simAnneal >= numberSteps simAnneal || driftChanges simAnneal >= driftMaxChanges simAnneal → do
                                        g ← selectionOf Unique $ newGraphList <> curSameBetterList
                                        pure (g, counter, newSAParams)

                                -- didn't hit stopping numbers so continuing--but based on current best cost not whatever was found
                                _ → do
                                    newBestGraph ← selectionOf Unique $ newGraphList <> curSameBetterList
                                    graphsToDo ← selectionOf Unique $ (newGraphList <> tailGraphs) `GO.reducedphylogeneticGraphListMinus` curSameBetterList
                                    swapAll'
                                        swapParams
                                        inGS
                                        inData
                                        (counter + 1)
                                        curBestCost
                                        newBestGraph
                                        graphsToDo
                                        numLeaves
                                        netPenaltyFactor
                                        newBreakEdgeNumber
                                        newSAParams

                            -- found only worse graphs--never happens due to the way splitjoin returns only better or equal
                            -- but could change
                            Nothing
                                | newMinCost > curBestCost →
                                    -- breakEdgeNumber set to zero for new graph to look at
                                    swapAll'
                                        swapParams
                                        inGS
                                        inData
                                        (counter + 1)
                                        curBestCost
                                        curSameBetterList
                                        tailGraphs
                                        numLeaves
                                        netPenaltyFactor
                                        0
                                        newSAParams
                            -- found same cost graphs
                            -- Important to not limit curSameBest since may rediscover graphs via swapping on equal when limiting the number to keep
                            -- can be a cause of infinite running issues.
                            Nothing →
                                let graphChoices = (tailGraphs <> newGraphList) `GO.reducedphylogeneticGraphListMinus` curSameBetterList
                                in  do
                                        newCurSameBetterList ← (Best `selectionOf`) $ curSameBetterList <> newGraphList
                                        graphsToDo ← Best `selectionOf` graphChoices

                                        -- these conditions help to prevent recswapping endlessly on new graphs thatare not in buffers,
                                        -- but have same cost
                                        let graphsToDo'
                                                | keepNum swapParams - 1 <= length graphsToDo = tailGraphs
                                                -- found nothing better that is new
                                                | keepNum swapParams - 1 <= length newCurSameBetterList = tailGraphs
                                                | otherwise = graphsToDo

                                        swapAll'
                                            swapParams
                                            inGS
                                            inData
                                            (counter + 1)
                                            curBestCost
                                            newCurSameBetterList
                                            graphsToDo'
                                            numLeaves
                                            netPenaltyFactor
                                            newBreakEdgeNumber
                                            newSAParams

                            -- simulated annealing/Drift post processing

                            -- found better cost graph
                            Just simAnneal | newMinCost < curBestCost → do
                                logWith LogInfo $ "\t->" <> show newMinCost
                                pure (newGraphList, counter, newSAParams)

                            -- not better so check for drift changes or annealing steps and return if reached maximum number
                            Just simAnneal | currentStep simAnneal >= numberSteps simAnneal || driftChanges simAnneal >= driftMaxChanges simAnneal → do
                                g ← selectionOf Unique $ newGraphList <> curSameBetterList
                                pure (g, counter, newSAParams)

                            -- didn't hit stopping numbers so continuing--but based on current best cost not whatever was found
                            Just simAnneal → do
                                newBestGraph ← selectionOf Unique $ newGraphList <> curSameBetterList
                                graphsToDo ← selectionOf Unique $ (newGraphList <> tailGraphs) `GO.reducedphylogeneticGraphListMinus` curSameBetterList
                                swapAll'
                                    swapParams
                                    inGS
                                    inData
                                    (counter + 1)
                                    curBestCost
                                    newBestGraph
                                    graphsToDo
                                    numLeaves
                                    netPenaltyFactor
                                    breakEdgeNumber
                                    newSAParams


{- | splitJoinGraph splits a graph on a single input edge (recursively though edge list) and rejoins to all possible other edges
if steepest == True then returns on finding a better graph (lower cost)
this will traverse entire SPR neighbohood if nothing better found (or steepest == False)
different from swapALL (original) in that it doesn't build up split list so lower memory footprint
breakEdgeList Complete keeps original edge list so can create readdition edge lists more easily
parallel map on rejoin if not steepest, if steepest do number of parallel threads so can reurn if any one is better
NB -- need to verify NNI/SPR/TBR rearrangement numbers
assumes break edges are bridge edges
graph split into two peices "base" graph with original root and "pruned" graph that was split off.
the edge connecting the two is (originalConnectionOfPruned -> prunedGraphRootIndex)
the edges in pruned graph do not contaiun that edge since are enumerated via preorder pass from prunedGraphRootIndex
this is the edge that is reconnected when graphs are joined, it is often delted and rejoined to update info and to deal with
conditions where the pruned graph is a single terminal
returns teh "breakEdgeNumber" so that in steepest, the edge breaking can continue in where it left off so to speak.
this can speed up SPR/TBR by a contant factor by not revisiting stable edges (not used in SA/drifting)
used in POY v1-3
-}
splitJoinGraph
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → VertexCost
    → [ReducedPhylogeneticGraph]
    → Int
    → VertexCost
    → Maybe SAParams
    → ReducedPhylogeneticGraph
    → Int
    → [LG.LEdge EdgeInfo]
    → [LG.LEdge EdgeInfo]
    → PhyG ([ReducedPhylogeneticGraph], Maybe SAParams, Int)
splitJoinGraph swapParams inGS inData curBestCost curSameBetterList numLeaves netPenaltyFactor inSimAnnealParams firstGraph breakEdgeNumber' breakEdgeListComplete = \case
    [] → pure (curSameBetterList, inSimAnnealParams, 0)
    -- split on first input edge
    edgeToBreakOn : otherEdges →
        let -- this so breaking edges can contnue where current left off
            -- since "rotates" edges all will be done.
            breakEdgeNumber = breakEdgeNumber' + 1

            -- split input graph into a part with the original root ("base") and the "pruned" graph -- the piece split off w/o original root
            (splitGraph, graphRoot, prunedGraphRootIndex, originalConnectionOfPruned) = LG.splitGraphOnEdge (thd5 firstGraph) edgeToBreakOn
        in  do
                -- reoptimize split graph for re-addition heuristics
                (reoptimizedSplitGraph, splitCost) ←
                    reoptimizeSplitGraphFromVertex inGS inData (doIA swapParams) netPenaltyFactor splitGraph graphRoot prunedGraphRootIndex

                -- check for malformed network split--do nothing if malformed
                if splitCost == infinity
                    then pure ([], inSimAnnealParams, breakEdgeNumber)
                    else do
                        -- regular swap
                        -- get root in base (for readdition) and edges in pruned section for rerooting during TBR readdition
                        let (_, edgesInPrunedGraph) = LG.nodesAndEdgesAfter splitGraph [(originalConnectionOfPruned, fromJust $ LG.lab splitGraph originalConnectionOfPruned)]

                        let edgesInBaseGraph = breakEdgeListComplete L.\\ (edgeToBreakOn : edgesInPrunedGraph)

                        -- insert here union calcuations based on Varon and Wheeler (2013)
                        -- basically--rebuild edge to rejoin list based on critical value, totalCost - splitCost, and
                        -- edge union distance
                        -- build edges pre-order and add to rejoin list if
                        -- 1) not network (but still recurse to children)
                        -- 2) union delta below threshold and recurse to children
                        -- if > threshold then stop, no add, no recurse since children can only get hihger ubnion distance
                        -- use split graph (with reoptimized nodes) and overall graph root to get avialbel edges in base graph for rejoin

                        let prunedToRejoinUnionData = vertData $ fromJust $ LG.lab reoptimizedSplitGraph prunedGraphRootIndex
                        -- prunedToRejoinUnionData = vertData $ fromJust $ LG.lab (thd5 firstGraph) prunedGraphRootIndex
                        let charInfoVV = fft5 firstGraph
                        unionEdgeList ←
                            getUnionRejoinEdgeList
                                inGS
                                reoptimizedSplitGraph
                                charInfoVV
                                [graphRoot]
                                ((snd5 firstGraph) - splitCost)
                                prunedToRejoinUnionData
                                []

                        -- builds graph edge list with unions--need to be able to turn off and just used edges in base graph for some sort
                        -- of "no-union" swap
                        -- determine those edges within distance of original if limited (ie NNI etc)
                        let veryBigDist = (maxMoveEdgeDist swapParams) >= ((maxBound ∷ Int) `div` 3)
                        -- NOTE: There might be a strictness issue here
                        let ~candidateEdges = take (maxMoveEdgeDist swapParams) $ LG.sortEdgeListByDistance splitGraph [graphRoot] [graphRoot]
                        let rejoinEdges' = case (veryBigDist, joinType swapParams) of
                                (True, JoinAll) → edgesInBaseGraph
                                (True, _) → unionEdgeList
                                (False, JoinAll) → candidateEdges
                                (False, _) → L.intersect candidateEdges unionEdgeList

                        -- randomize edges list order for anneal and drift
                        rejoinEdges ←
                            if atRandom swapParams
                                then shuffleList rejoinEdges'
                                else pure rejoinEdges'

                        -- rejoin graph to all possible edges in base graph
                        rejoinResult ←
                            {-# SCC rejoinResult #-}
                                rejoinGraph
                                    swapParams
                                    inGS
                                    inData
                                    curBestCost
                                    []
                                    netPenaltyFactor
                                    reoptimizedSplitGraph
                                    (GO.convertDecoratedToSimpleGraph splitGraph)
                                    splitCost
                                    graphRoot
                                    prunedGraphRootIndex
                                    originalConnectionOfPruned
                                    rejoinEdges
                                    edgesInPrunedGraph
                                    inSimAnnealParams
                        let (newGraphList, newSAParams) = rejoinResult
                        {-
                        trace ("Edge to break on:" <> (show $ LG.toEdge edgeToBreakOn)
                        <> "\nBase graph edges: " <> (show $ fmap LG.toEdge edgesInBaseGraph)
                        <> "\nPruned graph edges: " <> (show $ fmap LG.toEdge edgesInPrunedGraph)
                        <> "\nTarget edges to rejoin: " <> (show $ fmap LG.toEdge rejoinEdges)
                        <> "\nFull edgelist: " <> (show $ fmap LG.toEdge breakEdgeListComplete))
                        -}

                        newGraphList' ← GO.selectGraphs Best (outgroupIndex inGS) (keepNum swapParams) 0.0 newGraphList

                        case inSimAnnealParams of
                            -- only returns graphs if same of better else empty
                            -- adds null o\r better graphs to reurn list
                            Nothing
                                | (not . null) newGraphList && (steepest swapParams) →
                                    pure (newGraphList', inSimAnnealParams, breakEdgeNumber)
                            Nothing → do
                                splitJoinResult ←
                                    splitJoinGraph
                                        swapParams
                                        inGS
                                        inData
                                        curBestCost
                                        curSameBetterList
                                        numLeaves
                                        netPenaltyFactor
                                        inSimAnnealParams
                                        firstGraph
                                        breakEdgeNumber
                                        breakEdgeListComplete
                                        otherEdges
                                let (recurseGraphList, _, newEdgeBreakNumber) = splitJoinResult
                                pure (newGraphList' <> recurseGraphList, inSimAnnealParams, newEdgeBreakNumber)
                            -- Annealing/Drift swap
                            -- return if better
                            -- return if other chosen probabalistically
                            -- recurse if nothing returned

                            -- if better than current graph return
                            Just simAnnealParams →
                                let newMinCost
                                        | null newGraphList = infinity
                                        | otherwise = (snd5 . head) newGraphList'

                                    result
                                        | newMinCost < curBestCost = pure (newGraphList', newSAParams, breakEdgeNumber)
                                        -- if SA returned graphs--return them
                                        | (not . null) newGraphList = pure (newGraphList, newSAParams, breakEdgeNumber)
                                        -- keep going if nothing
                                        | otherwise =
                                            splitJoinGraph
                                                swapParams
                                                inGS
                                                inData
                                                -- (tail randomIntListSwap)
                                                curBestCost
                                                curSameBetterList
                                                numLeaves
                                                netPenaltyFactor
                                                newSAParams
                                                firstGraph
                                                breakEdgeNumber
                                                breakEdgeListComplete
                                                otherEdges
                                in  result


{- | getUnionRejoinEdgeList takes a graph (split and reoptimized usually), the overall root index (of split),
split cost, and union threshold value and returns list of edges that have union distance <= threshold factor
assumes that root edges are not network edges (an invariant)
checks node then recurses to children
-}
getUnionRejoinEdgeList
    ∷ GlobalSettings
    → DecoratedGraph
    → V.Vector (V.Vector CharInfo)
    → [LG.Node]
    → Double
    → VertexBlockData
    → [LG.LEdge EdgeInfo]
    → PhyG [LG.LEdge EdgeInfo]
getUnionRejoinEdgeList inGS inGraph charInfoVV nodeIndexList splitDiffCost nodeToJoinUnionData curEdgeList =
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
                                                getUnionRejoinEdgeList
                                                    inGS
                                                    inGraph
                                                    charInfoVV
                                                    [(snd3 $ last childEdges)]
                                                    splitDiffCost
                                                    nodeToJoinUnionData
                                                    ((last childEdges) : curEdgeList)
                                            else
                                                getUnionRejoinEdgeList
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
                                        then getUnionRejoinEdgeList inGS inGraph charInfoVV childNodeIndexList splitDiffCost nodeToJoinUnionData (childEdges <> curEdgeList)
                                        else pure curEdgeList
                                -- in  -- recurse remaining nodes
                                getUnionRejoinEdgeList inGS inGraph charInfoVV (tail nodeIndexList) splitDiffCost nodeToJoinUnionData newCurEdgeListChild


{- | getUnionDistanceM gets distance between two the union fields of two characters
since its a distance no need for no change cost adjustment
-}
getUnionDistanceM ∷ VertexBlockData → VertexBlockData → V.Vector (V.Vector CharInfo) → PhyG Double
getUnionDistanceM union1 union2 charInfoVV =
    let noChangeCostAdjut = False
    in  M.distance2UnionsM noChangeCostAdjut union1 union2 charInfoVV


-- | rejoinGraphTuple is a wrapper around rejoinGraph for fmapping--only returns graph list not simulated annealing params
rejoinGraphTuple
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → VertexCost
    → [ReducedPhylogeneticGraph]
    → Maybe SAParams
    → (DecoratedGraph, SimpleGraph, VertexCost, LG.Node, LG.Node, LG.Node, [LG.LEdge EdgeInfo], [LG.LEdge EdgeInfo], VertexCost)
    → PhyG [ReducedPhylogeneticGraph]
rejoinGraphTuple
    swapParams
    inGS
    inData
    curBestCost
    curBestGraphs
    inSimAnnealParams
    ( reoptimizedSplitGraph
        , splitGraphSimple
        , splitGraphCost
        , graphRoot
        , prunedGraphRootIndex
        , originalConnectionOfPruned
        , rejoinEdges
        , edgesInPrunedGraph
        , netPenaltyFactor
        ) = do
        result ←
            rejoinGraph
                swapParams
                inGS
                inData
                curBestCost
                curBestGraphs
                netPenaltyFactor
                reoptimizedSplitGraph
                splitGraphSimple
                splitGraphCost
                graphRoot
                prunedGraphRootIndex
                originalConnectionOfPruned
                rejoinEdges
                edgesInPrunedGraph
                inSimAnnealParams
        pure $ fst result


{- | rejoinGraph rejoins a split graph at all edges (if not steepest and found better)
in "base" graph.
if not steepest then do all as map, else recursive on base graph edge list
nni doesn't apper to be correct here--maybe loose it--doing nothing
-}
rejoinGraph
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → VertexCost
    → [ReducedPhylogeneticGraph]
    → VertexCost
    → DecoratedGraph
    → SimpleGraph
    → VertexCost
    → LG.Node
    → LG.Node
    → LG.Node
    → [LG.LEdge EdgeInfo]
    → [LG.LEdge EdgeInfo]
    → Maybe SAParams
    → PhyG ([ReducedPhylogeneticGraph], Maybe SAParams)
rejoinGraph swapParams inGS inData curBestCost curBestGraphs netPenaltyFactor reoptimizedSplitGraph splitGraphSimple splitGraphCost graphRoot prunedGraphRootIndex originalConnectionOfPruned rejoinEdges' edgesInPrunedGraph inSimAnnealParams =
    {-# SCC rejoinGraph_TOP_DEF #-}
    -- found no better--but return equal cost graphs
    -- trace ("In rejoinGraph with num rejoining edges: " <> (show $ length rejoinEdges')) (
    if null rejoinEdges'
        then pure (curBestGraphs, inSimAnnealParams)
        else -- this is for no  swapping option in fuse and genetic algorithm-fuse

            let rejoinEdges =
                    if (swapType swapParams) == NoSwap
                        then take 6 rejoinEdges'
                        else rejoinEdges'
            in  -- regular swapping
                if isNothing inSimAnnealParams
                    then -- check if split graph cost same as graph then return since graph can only get longer on readdition

                        if splitGraphCost >= curBestCost
                            then pure ([], inSimAnnealParams)
                            else -- fmap over all edges in base graph

                                if not (steepest swapParams)
                                    then
                                        let {-
                                             TODO: Add back safe parallelism here
                                                 * Old implementation (unsafe paralel):
                                                   rejoinGraphList = concat $ fmap fst $ PU.seqParMap (parStrategy $ lazyParStrat inGS) (singleJoin swapParams inGS inData reoptimizedSplitGraph splitGraphSimple splitGraphCost prunedGraphRootIndex originalConnectionOfPruned curBestCost edgesInPrunedGraph inSimAnnealParams) rejoinEdges
                                                 * New implementation (safe sequential):
                                             -}
                                            -- parallel stuff
                                            action ∷ LG.LEdge EdgeInfo → PhyG [ReducedPhylogeneticGraph]
                                            action =
                                                {-# SCC rejoinGraph_action_of_singleJoin_1 #-}
                                                fmap fst
                                                    . singleJoin
                                                        swapParams
                                                        inGS
                                                        inData
                                                        reoptimizedSplitGraph
                                                        splitGraphSimple
                                                        splitGraphCost
                                                        prunedGraphRootIndex
                                                        originalConnectionOfPruned
                                                        curBestCost
                                                        edgesInPrunedGraph
                                                        inSimAnnealParams
                                        in  do
                                                -- logWith LogInfo "In rejoinGraph-not-steepest"
                                                rejoinGraphList ←
                                                    getParallelChunkTraverseBy (listApplying strict2of5) >>= \pTraverse →
                                                        fold <$> pTraverse action rejoinEdges
                                                {-
                                                rejoinOperation = fst . singleJoin
                                                    swapParams
                                                    inGS
                                                    inData
                                                    reoptimizedSplitGraph
                                                    splitGraphSimple
                                                    splitGraphCost
                                                    prunedGraphRootIndex
                                                    originalConnectionOfPruned
                                                    curBestCost
                                                    edgesInPrunedGraph
                                                    inSimAnnealParams
                                                 rejoinGraphList = foldMap rejoinOperation rejoinEdges
                                                 -}

                                                {-Checking only min but seems to make slower
                                                newMinCost = if null rejoinGraphList then infinity
                                                             else minimum $ fmap snd rejoinGraphList
                                                (minEstCostNewGraphList, _) = unzip $ filter ((== newMinCost) . snd) rejoinGraphList
                                                -}

                                                -- newGraphList = fmap (T.multiTraverseFullyLabelGraphReduced inGS inData False False Nothing) (fmap fst rejoinGraphList) `using` PU.myParListChunkRDS
                                                newGraphList ← GO.selectGraphs Best (outgroupIndex inGS) (keepNum swapParams) 0 rejoinGraphList

                                                -- will only return graph if <= curBest cost
                                                case rejoinGraphList of
                                                    [] → pure ([], inSimAnnealParams)
                                                    xs → case getGraphCost newGraphList `compare` curBestCost of
                                                        LT → pure (newGraphList, inSimAnnealParams)
                                                        _ → pure ([], inSimAnnealParams)
                                    else -- famp over number of threads edges in base graph
                                    -- then recurse

                                        let -- this could be made a little parallel--but if lots of threads basically can do all
                                            -- to not overload paralle threads
                                            {-  This not so efficient is swapping in single graphs so leaving it be
                                            saRounds = if isNothing inSimAnnealParams then 1
                                                       else rounds $ fromJust inSimAnnealParams

                                            (numGraphsToExamine, _) = divMod PU.getNumThreads saRounds -- this may not "drift" if finds alot better, but that's how its supposed to work
                                            -}
                                            numGraphsToExamine = (graphsSteepest inGS) -- min (graphsSteepest inGS) PU.getNumThreads
                                            rejoinEdgeList = take numGraphsToExamine rejoinEdges

                                            -- parallel stuff
                                            action ∷ LG.LEdge EdgeInfo → PhyG [ReducedPhylogeneticGraph]
                                            action =
                                                {-# SCC rejoinGraph_action_of_singleJoin_2 #-}
                                                fmap fst
                                                    . singleJoin
                                                        swapParams
                                                        inGS
                                                        inData
                                                        reoptimizedSplitGraph
                                                        splitGraphSimple
                                                        splitGraphCost
                                                        prunedGraphRootIndex
                                                        originalConnectionOfPruned
                                                        curBestCost
                                                        edgesInPrunedGraph
                                                        inSimAnnealParams
                                        in  do
                                                -- logWith LogInfo $ "In rejoinGraph-steepest " <> (show $ length rejoinEdgeList)
                                                rejoinGraphList ←
                                                    getParallelChunkTraverseBy (listApplying strict1and2of5) >>= \pTraverse →
                                                        fold <$> pTraverse action rejoinEdgeList

                                                newGraphList' ← GO.selectGraphs Best (outgroupIndex inGS) (keepNum swapParams) 0.0 rejoinGraphList

                                                -- found nothing better or equal
                                                if null rejoinGraphList
                                                    then -- trace ("In steepest worse: " <> (show $ length (drop PU.getNumThreads rejoinEdges)))

                                                        rejoinGraph
                                                            swapParams
                                                            inGS
                                                            inData
                                                            curBestCost
                                                            curBestGraphs
                                                            netPenaltyFactor
                                                            reoptimizedSplitGraph
                                                            splitGraphSimple
                                                            splitGraphCost
                                                            graphRoot
                                                            prunedGraphRootIndex
                                                            originalConnectionOfPruned
                                                            (drop numGraphsToExamine rejoinEdges)
                                                            edgesInPrunedGraph
                                                            inSimAnnealParams
                                                    else -- found better graph

                                                        if (snd5 . head) newGraphList' < curBestCost
                                                            then -- trace ("Steepest better")
                                                                pure (newGraphList', inSimAnnealParams)
                                                            else -- found equal cost graph

                                                                if (snd5 . head) newGraphList' == curBestCost
                                                                    then do
                                                                        newBestList ← GO.selectGraphs Best (outgroupIndex inGS) (keepNum swapParams) 0 $ curBestGraphs <> newGraphList'
                                                                        rejoinGraph
                                                                            swapParams
                                                                            inGS
                                                                            inData
                                                                            curBestCost
                                                                            newBestList
                                                                            netPenaltyFactor
                                                                            reoptimizedSplitGraph
                                                                            splitGraphSimple
                                                                            splitGraphCost
                                                                            graphRoot
                                                                            prunedGraphRootIndex
                                                                            originalConnectionOfPruned
                                                                            (drop numGraphsToExamine rejoinEdges)
                                                                            edgesInPrunedGraph
                                                                            inSimAnnealParams
                                                                    else -- found worse graphs only

                                                                    -- trace ("In steepest worse (after recalculation): " <> (show $ length (drop PU.getNumThreads rejoinEdges)))

                                                                        rejoinGraph
                                                                            swapParams
                                                                            inGS
                                                                            inData
                                                                            curBestCost
                                                                            curBestGraphs
                                                                            netPenaltyFactor
                                                                            reoptimizedSplitGraph
                                                                            splitGraphSimple
                                                                            splitGraphCost
                                                                            graphRoot
                                                                            prunedGraphRootIndex
                                                                            originalConnectionOfPruned
                                                                            (drop numGraphsToExamine rejoinEdges)
                                                                            edgesInPrunedGraph
                                                                            inSimAnnealParams
                    else -- Drifting/Simulated annealing
                    -- basically if it accepted a graph (better or probabalistically) then pass up
                    -- otherwise move to next rejoin, changes in graphs counted at higher level

                        let -- based on "steepest"
                            -- to not overload paralle threads
                            {-  This not so efficient is swapping in single graphs so leaving it be
                            saRounds = if isNothing inSimAnnealParams then 1
                                       else rounds $ fromJust inSimAnnealParams

                            (numGraphsToExamine, _) = divMod PU.getNumThreads saRounds -- this may not "drift" if finds alot better, but that's how its supposed to work
                            -}
                            numGraphsToExamine = graphsSteepest inGS -- min (graphsSteepest inGS) PU.getNumThreads
                            rejoinEdgeList = take numGraphsToExamine rejoinEdges
                            simAnnealParamList = U.generateUniqueRandList numGraphsToExamine inSimAnnealParams

                            -- parallel stuff
                            action ∷ (Maybe SAParams, LG.LEdge EdgeInfo) → PhyG ([ReducedPhylogeneticGraph], Maybe SAParams)
                            action =
                                {-# SCC rejoinGraph_action_of_singleJoinPrime #-}
                                singleJoin'
                                    swapParams
                                    inGS
                                    inData
                                    reoptimizedSplitGraph
                                    splitGraphSimple
                                    splitGraphCost
                                    prunedGraphRootIndex
                                    originalConnectionOfPruned
                                    curBestCost
                                    edgesInPrunedGraph
                        in  do
                                rejoinGraphPairList ←
                                    getParallelChunkTraverseBy (applyOver1of2 (listApplying strict2of5)) >>= \pTraverse →
                                        pTraverse action $ zip simAnnealParamList rejoinEdgeList

                                -- mechanics to see if trhere is a better graph in return set
                                -- only taking first of each list--so can keep sa params with them--really all should have length == 1 anyway
                                -- making sure remove all null lists that nothing was found
                                let nonEmptyPairList = filter (not . null . fst) rejoinGraphPairList

                                let rejoinGraphList =
                                        if (not . null) nonEmptyPairList
                                            then fmap (head . fst) nonEmptyPairList
                                            else []

                                let newMinCost =
                                        if (not . null) rejoinGraphList
                                            then minimum $ fmap snd5 rejoinGraphList
                                            else infinity

                                -- head should only be called when non-empty--so should never get runtime head error
                                let (newMinGraph, newMinGraphSAParams) = head $ filter ((== newMinCost) . snd5 . fst) (zip rejoinGraphList (fmap snd nonEmptyPairList))

                                -- if better than current--pass up and on
                                if newMinCost < curBestCost
                                    then pure ([newMinGraph], newMinGraphSAParams)
                                    else -- check if hit step limit--more for SA than drift

                                        if ((currentStep $ fromJust inSimAnnealParams) >= (numberSteps $ fromJust inSimAnnealParams))
                                            || ((driftChanges $ fromJust inSimAnnealParams) >= (driftMaxChanges $ fromJust inSimAnnealParams))
                                            then pure (curBestGraphs, inSimAnnealParams)
                                            else -- not better so go to SA results

                                            -- return first non-empty result

                                                if (not . null) nonEmptyPairList
                                                    then pure $ head nonEmptyPairList
                                                    else -- if nothing returned (no better or probabalistically chosen) go on with updated SA params

                                                        let newSAParams = (snd . head) rejoinGraphPairList
                                                        in  rejoinGraph
                                                                swapParams
                                                                inGS
                                                                inData
                                                                curBestCost
                                                                curBestGraphs
                                                                netPenaltyFactor
                                                                reoptimizedSplitGraph
                                                                splitGraphSimple
                                                                splitGraphCost
                                                                graphRoot
                                                                prunedGraphRootIndex
                                                                originalConnectionOfPruned
                                                                (drop numGraphsToExamine rejoinEdges)
                                                                edgesInPrunedGraph
                                                                newSAParams


-- | singleJoin' is a wrapper arounds singleJoin to allow parMap with individual SAParams
singleJoin'
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → DecoratedGraph
    → SimpleGraph
    → VertexCost
    → LG.Node
    → LG.Node
    → VertexCost
    → [LG.LEdge EdgeInfo]
    → (Maybe SAParams, LG.LEdge EdgeInfo)
    → PhyG ([ReducedPhylogeneticGraph], Maybe SAParams)
singleJoin' swapParams inGS inData splitGraph splitGraphSimple splitCost prunedGraphRootIndex originalConnectionOfPruned curBestCost edgesInPrunedGraph (inSimAnnealParams, targetEdge) =
    singleJoin
        swapParams
        inGS
        inData
        splitGraph
        splitGraphSimple
        splitCost
        prunedGraphRootIndex
        originalConnectionOfPruned
        curBestCost
        edgesInPrunedGraph
        inSimAnnealParams
        targetEdge


{- | singleJoin takes optimized split graph, split cost, target edge, swap type (ie TBR/SPR/NNI)
and "rejoins" the split graph to a single graph--creates joined graph and calculates a heuristic graph cost
based on the union assignment of the edge and its distance to the root vertex of the pruned graph
if TBR checks all edges in pruned graph with readdition edge (shortcircuits if steepest  == True)
always deletes connecting edge to pruned part and readds--this because sometimes it is there and sometimes not (depending on
if SPR for terminal etc) and can create parallel edges with different weights (0.0 or not) so just remove to be sure.
TBR uses dynamic epsilon even in SPR moves--SPR does not
-}
singleJoin
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → DecoratedGraph
    → SimpleGraph
    → VertexCost
    → LG.Node
    → LG.Node
    → VertexCost
    → [LG.LEdge EdgeInfo]
    → Maybe SAParams
    → LG.LEdge EdgeInfo
    → PhyG ([ReducedPhylogeneticGraph], Maybe SAParams)
singleJoin swapParams inGS inData splitGraph splitGraphSimple splitCost prunedGraphRootIndex originalConnectionOfPruned curBestCost edgesInPrunedGraph inSimAnnealParams targetEdge@(u, v, _)
    | LG.isEmpty splitGraphSimple = pure ([], inSimAnnealParams)
    -- do redo orginal graph join
    | originalConnectionOfPruned `elem` [u, v] = pure ([], inSimAnnealParams)
    | otherwise = case LG.lab splitGraph prunedGraphRootIndex of
        Nothing → pure ([], inSimAnnealParams)
        -- regular swap
        Just label →
            let defaultResult = (mempty, inSimAnnealParams)

                newEdgeList =
                    [ (u, originalConnectionOfPruned, 0.0)
                    , (originalConnectionOfPruned, v, 0.0)
                    , (originalConnectionOfPruned, prunedGraphRootIndex, 0.0)
                    ]

                charInfoVV = fmap thd3 $ thd3 inData

                -- Filter for bridge edges for TBR when needed
                edgesInPrunedGraph'
                    | (graphType inGS == Tree) || LG.isTree splitGraphSimple = edgesInPrunedGraph
                    | otherwise = fmap fst . filter snd . zip edgesInPrunedGraph $ LG.isBridge splitGraphSimple . LG.toEdge <$> edgesInPrunedGraph

                sprNewGraph = LG.insEdges newEdgeList $ LG.delEdges [(u, v), (originalConnectionOfPruned, prunedGraphRootIndex)] splitGraphSimple

                -- here when needed--correct graph is issue in network
                -- swap can screw up time consistency and other issues
                getCheckedGraphNewSPR ∷ PhyG (LG.Gr NameText VertexCost)
                getCheckedGraphNewSPR = do
                    isPhyloGraph ← LG.isPhylogeneticGraph sprNewGraph
                    let result
                            | graphType inGS == Tree = sprNewGraph
                            | isPhyloGraph = sprNewGraph
                            | otherwise = LG.empty
                    pure result

                decide ∷ ReducedPhylogeneticGraph → ([ReducedPhylogeneticGraph], Maybe SAParams)
                decide input@(_, newCost, _, _, _)
                    | newCost <= curBestCost = ([input], inSimAnnealParams)
                    | otherwise = defaultResult

                action ∷ PhyG ([ReducedPhylogeneticGraph], Maybe SAParams)
                action = do
                    sprNewGraphChecked ← getCheckedGraphNewSPR
                    if LG.isEmpty sprNewGraphChecked
                        then pure defaultResult
                        else decide <$> T.multiTraverseFullyLabelGraphReduced inGS inData False False Nothing sprNewGraphChecked
                -- graphType with IA field
                -- only uswe wqhere they exist

                (makeEdgeDataFunction, edgeJoinFunction) =
                    if graphType inGS == HardWired
                        then (M.makeEdgeDataM False True, edgeJoinDelta inGS False)
                        else
                            if not (useIA inGS)
                                then (M.makeEdgeDataM False True, edgeJoinDelta inGS False)
                                else (M.makeEdgeDataM True True, edgeJoinDelta inGS True)
            in  do
                    -- set edge union creation type to IA-based, filtering gaps (should be linear)
                    -- hence True True
                    targetEdgeData ← makeEdgeDataFunction splitGraph charInfoVV targetEdge
                    -- this for using DO for edge O(n^2)
                    -- targetEdgeData = M.makeEdgeData doIA (not doIA) splitGraph charInfoVV targetEdge

                    -- this for SPR/NNI only
                    let prunedRootVertexData = vertData $ fromJust $ LG.lab splitGraph prunedGraphRootIndex

                    -- rejoin should always be DO based on edge and pruned root but can be different lengths (unless Static Approx)
                    sprReJoinCost ← edgeJoinFunction charInfoVV prunedRootVertexData targetEdgeData

                    case inSimAnnealParams of
                        -- SPR or no TBR rearrangements
                        -- wierdness here with first case should have been taking longer--but is quicker for prealigned
                        Nothing → case swapType swapParams of
                            -- _ | length edgesInPrunedGraph < 4 → action
                            SPR | sprReJoinCost + splitCost <= curBestCost → action
                            SPR → pure defaultResult
                            -- do TBR stuff returning SPR results if heuristic better
                            TBR | (length edgesInPrunedGraph < 4) && (sprReJoinCost + splitCost <= curBestCost) → action
                            TBR | (length edgesInPrunedGraph < 4) → pure defaultResult
                            TBR → do
                                -- check if spr better always return if so
                                sprResult ←
                                    if sprReJoinCost + splitCost <= curBestCost + (sprReJoinCost * (dynamicEpsilon inGS))
                                        then fst <$> action
                                        else pure mempty
                                case toList sprResult of
                                    _ : _ → pure (sprResult, inSimAnnealParams)
                                    [] → do
                                        tbrResult' ←
                                            tbrJoin
                                                swapParams
                                                inGS
                                                inData
                                                splitGraph
                                                splitGraphSimple
                                                splitCost
                                                prunedGraphRootIndex
                                                originalConnectionOfPruned
                                                curBestCost
                                                edgesInPrunedGraph'
                                                inSimAnnealParams
                                                targetEdge
                                        let (tbrResult, _) = tbrResult'
                                        pure (tbrResult, inSimAnnealParams)

                            -- TBRAlternate can skip SPR moves since done already in alternate scenario
                            _ → do
                                tbrResult' ←
                                    tbrJoin
                                        swapParams
                                        inGS
                                        inData
                                        splitGraph
                                        splitGraphSimple
                                        splitCost
                                        prunedGraphRootIndex
                                        originalConnectionOfPruned
                                        curBestCost
                                        edgesInPrunedGraph'
                                        inSimAnnealParams
                                        targetEdge
                                let (tbrResult, _) = tbrResult'
                                pure (tbrResult, inSimAnnealParams)

                        -- simulated annealing/Drift swap
                        Just _simAnnealParams → do
                            -- check if spr better always return if so
                            sprNewGraphChecked ← getCheckedGraphNewSPR
                            rediagnosedSPRGraph@(_, newCost, _, _, _) ←
                                T.multiTraverseFullyLabelGraphReduced inGS inData False False Nothing sprNewGraphChecked

                            let sprResult
                                    | sprReJoinCost + splitCost <= curBestCost + (sprReJoinCost * (dynamicEpsilon inGS)) && newCost < curBestCost =
                                        [rediagnosedSPRGraph]
                                    | otherwise = mempty
                            case toList sprResult of
                                -- spr better than current
                                (_, cost, _, _, _) : _ → do
                                    (_, newSAParams) ← U.simAnnealAccept inSimAnnealParams curBestCost cost
                                    pure (sprResult, newSAParams)
                                -- do simAnneal/Drift for SPR and on to tbr if not accept
                                [] → do
                                    acceptance ← U.simAnnealAccept inSimAnnealParams curBestCost (sprReJoinCost + splitCost)
                                    case acceptance of
                                        -- if accepted (better or random) then return with updated annealing/Drift parameters
                                        (True, newSAParams) → pure ([rediagnosedSPRGraph], newSAParams)
                                        -- rejected-recurse with updated SA params
                                        (False, newSAParams) → case swapType swapParams of
                                            SPR → pure (mempty, newSAParams)
                                            _ | length edgesInPrunedGraph < 4 → pure (mempty, newSAParams)
                                            _ →
                                                tbrJoin
                                                    swapParams
                                                    inGS
                                                    inData
                                                    splitGraph
                                                    splitGraphSimple
                                                    splitCost
                                                    prunedGraphRootIndex
                                                    originalConnectionOfPruned
                                                    curBestCost
                                                    edgesInPrunedGraph'
                                                    newSAParams
                                                    targetEdge


{- | edgeJoinDelta calculates heuristic cost for joining pair edges
    IA field is faster--but has to be there and not for harwired
-}
edgeJoinDelta ∷ GlobalSettings → Bool → V.Vector (V.Vector CharInfo) → VertexBlockData → VertexBlockData → PhyG VertexCost
edgeJoinDelta inGS useIA charInfoVV edgeA edgeB =
    if (not useIA)
        then do
            vertexStuff ← POSW.createVertexDataOverBlocks inGS edgeA edgeB charInfoVV []
            pure $ V.sum $ fmap V.sum $ fmap (fmap snd) vertexStuff
        else do
            vertexStuff ← POSW.createVertexDataOverBlocksStaticIA inGS edgeA edgeB charInfoVV []
            pure $ V.sum $ fmap V.sum $ fmap (fmap snd) vertexStuff


{- | tbrJoin performs TBR rearrangements on pruned graph component
"reroots" pruned graph on each bridge edge and tries join to
target edge as in SPR
each edge is tried in turn (except for original root edge covered by singleJoin SPR function)
if heuristic edge join cost is below current best cost then the component is rerooted, joined to
target edge and graph fully diagnosed to verify cost
"steepest" short circuits checking edges to return better verified cost graph immediately
otherwise ("all") returns all graphs better than or equal to current better cost
if nothing equal or better found returns empty list
tests if reroot edges are bridge edges
uses dynamic epsilon--seems the delta estimate is high
-}
tbrJoin
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → DecoratedGraph
    → SimpleGraph
    → VertexCost
    → LG.Node
    → LG.Node
    → VertexCost
    → [LG.LEdge EdgeInfo]
    → Maybe SAParams
    → LG.LEdge EdgeInfo
    → PhyG ([ReducedPhylogeneticGraph], Maybe SAParams)
tbrJoin swapParams inGS inData splitGraph splitGraphSimple splitCost prunedGraphRootIndex originalConnectionOfPruned curBestCost edgesInPrunedGraph inSimAnnealParams targetEdge =
    -- this is for networks stopping TBR rearrangements with there are network edegs involved in pruned part of graph
    let hasNetEdges
            | graphType inGS == Tree = False
            | otherwise = null $ filter ((== True) . LG.isNetworkLabEdge splitGraph) edgesInPrunedGraph
    in  case edgesInPrunedGraph of
            [] → pure ([], inSimAnnealParams)
            _ | hasNetEdges → pure ([], inSimAnnealParams)
            _ →
                -- get target edge data\
                -- always using IA for union, but filtering out gaps (doIA (not doIA))
                let charInfoVV = fmap thd3 $ thd3 inData

                    -- graphType with IA field
                    -- only uswe wqhere they exist\
                    (makeEdgeDataFunction, edgeJoinFunction)
                        | graphType inGS == HardWired = (M.makeEdgeDataM False True, edgeJoinDelta inGS False)
                        | not (useIA inGS) = (M.makeEdgeDataM False True, edgeJoinDelta inGS False)
                        | otherwise = (M.makeEdgeDataM True True, edgeJoinDelta inGS True)
                in  do
                        targetEdgeData ← makeEdgeDataFunction splitGraph charInfoVV targetEdge

                        -- parallell stuff
                        let makeEdgeAction ∷ LG.LEdge b → PhyG VertexBlockData
                            makeEdgeAction = makeEdgeDataFunction splitGraph charInfoVV

                        let joinAction ∷ VertexBlockData → PhyG VertexCost
                            joinAction = edgeJoinFunction charInfoVV targetEdgeData

                        let rerootAction ∷ LG.LEdge EdgeInfo → SimpleGraph
                            rerootAction = rerootPrunedAndMakeGraph splitGraphSimple prunedGraphRootIndex originalConnectionOfPruned targetEdge

                        -- Debugging info
                        -- debugger ∷ (Logger m, Show a, Show b, Show c) ⇒ String -> a → LG.Gr b c → m ()
                        -- debugger str n g = logWith LogTech $ fold
                        --    [ "In: 'tbrJoin' with '", str, "'\n  Graph [ G_", show n, " ]:\n", LG.prettify g ]

                        -- reoptimizeAction ∷ SimpleGraph → PhyG ReducedPhylogeneticGraph
                        let reoptimizeAction g0 = do
                                -- debugger "reoptimizeAction" 0 g0
                                result@(g1, _, _, _, _) ← T.multiTraverseFullyLabelGraphReduced inGS inData False False Nothing g0
                                -- debugger "reoptimizeAction" 1 g1
                                pure result
                        -- logic for annealing/Drift  regular swap first
                        case inSimAnnealParams of
                            -- get heuristic delta joins for edges in pruned graph
                            Nothing
                                | not (steepest swapParams) →
                                    let rerootEdgeList = filter ((/= prunedGraphRootIndex) . fst3) $ filter ((/= originalConnectionOfPruned) . fst3) edgesInPrunedGraph
                                    in  do
                                            -- debugger "CASE OF -> Nothing( 1 )" 0 splitGraphSimple
                                            -- True True to use IA fields and filter gaps
                                            rerootEdgeDataList ←
                                                getParallelChunkTraverse >>= \pTraverse →
                                                    makeEdgeAction `pTraverse` rerootEdgeList

                                            rerootEdgeDeltaList' ←
                                                getParallelChunkTraverse >>= \pTraverse →
                                                    joinAction `pTraverse` rerootEdgeDataList

                                            let rerootEdgeDeltaList = fmap (+ splitCost) rerootEdgeDeltaList'

                                            let joinBestCost = minimum $ curBestCost : rerootEdgeDeltaList
                                            -- logWith LogInfo $ "TBRJoin JBC CBC->" <> show (curBestCost, joinBestCost)

                                            -- check for possible better/equal graphs and verify
                                            let deltaAdjustmentJoinCost = (curBestCost - splitCost) * (dynamicEpsilon inGS)
                                            let candidateEdgeList = fmap fst $ filter ((<= (curBestCost + deltaAdjustmentJoinCost)) . snd) (zip rerootEdgeList rerootEdgeDeltaList)

                                            -- /NOTE:/ All returned fields are required by 'reoptimizeAction', so evaluate strictly
                                            candidateJoinedGraphList' ←
                                                getParallelChunkMap <&> \pMap →
                                                    rerootAction `pMap` candidateEdgeList

                                            -- check for graph wierdness if network
                                            candidateJoinedGraphList ←
                                                if graphType inGS == Tree
                                                    then pure candidateJoinedGraphList'
                                                    else filterM LG.isPhylogeneticGraph candidateJoinedGraphList'

                                            -- check for graph wierdness
                                            rediagnosedGraphList' ←
                                                getParallelChunkTraverseBy strict2of5 >>= \pTraverse →
                                                    reoptimizeAction `pTraverse` candidateJoinedGraphList
                                            let rediagnosedGraphList = filter ((<= curBestCost) . snd5) rediagnosedGraphList'

                                            let result
                                                    | null candidateEdgeList = []
                                                    | null rediagnosedGraphList = []
                                                    | otherwise = rediagnosedGraphList
                                            pure (result, Nothing)

                            -- get steepest edges
                            Nothing →
                                let -- to not overload paralle threads
                                    {-  This not so efficient is swapping in single graphs so leaving it be
                                    saRounds = if isNothing inSimAnnealParams then 1
                                               else rounds $ fromJust inSimAnnealParams

                                    (numGraphsToExamine, _) = divMod PU.getNumThreads saRounds -- this may not "drift" if finds alot better, but that's how its supposed to work
                                    -}
                                    numEdgesToExamine = graphsSteepest inGS -- min (graphsSteepest inGS) PU.getNumThreads
                                    firstSetEdges = take numEdgesToExamine edgesInPrunedGraph

                                    -- get heuristic delta joins for steepest edge set
                                    rerootEdgeList = filter ((/= prunedGraphRootIndex) . fst3) $ filter ((/= originalConnectionOfPruned) . fst3) firstSetEdges
                                in  do
                                        -- debugger "CASE OF -> Nothing( 2 )" 0 splitGraphSimple
                                        -- True True to use IA fields and filter gaps

                                        rerootEdgeDataList ←
                                            getParallelChunkTraverse >>= \pTraverse →
                                                makeEdgeAction `pTraverse` rerootEdgeList

                                        rerootEdgeDeltaList' ←
                                            getParallelChunkTraverse >>= \pTraverse →
                                                joinAction `pTraverse` rerootEdgeDataList

                                        let rerootEdgeDeltaList = fmap (+ splitCost) rerootEdgeDeltaList'

                                        -- check for possible better/equal graphs and verify
                                        let deltaAdjustmentJoinCost = (curBestCost - splitCost) * (dynamicEpsilon inGS)
                                        let candidateEdgeList = fmap fst $ filter ((<= (curBestCost + deltaAdjustmentJoinCost)) . snd) (zip rerootEdgeList rerootEdgeDeltaList)

                                        -- /NOTE:/ All returned fields are required by 'reoptimizeAction', so evaluate strictly
                                        candidateJoinedGraphList ←
                                            getParallelChunkMap <&> \pMap →
                                                rerootAction `pMap` candidateEdgeList

                                        rediagnosedGraphList' ←
                                            getParallelChunkTraverseBy strict2of5 >>= \pTraverse →
                                                reoptimizeAction `pTraverse` candidateJoinedGraphList
                                        let rediagnosedGraphList = filter ((<= curBestCost) . snd5) rediagnosedGraphList'

                                        -- trace ("TBR steepest: " <> (show $ length rerootEdgeList) <> " edges to go " <> (show $ length $ (drop numEdgesToExamine edgesInPrunedGraph))) (
                                        if null candidateEdgeList
                                            then
                                                tbrJoin
                                                    swapParams
                                                    inGS
                                                    inData
                                                    splitGraph
                                                    splitGraphSimple
                                                    splitCost
                                                    prunedGraphRootIndex
                                                    originalConnectionOfPruned
                                                    curBestCost
                                                    (drop numEdgesToExamine edgesInPrunedGraph)
                                                    inSimAnnealParams
                                                    targetEdge
                                            else
                                                if null rediagnosedGraphList
                                                    then
                                                        tbrJoin
                                                            swapParams
                                                            inGS
                                                            inData
                                                            splitGraph
                                                            splitGraphSimple
                                                            splitCost
                                                            prunedGraphRootIndex
                                                            originalConnectionOfPruned
                                                            curBestCost
                                                            (drop numEdgesToExamine edgesInPrunedGraph)
                                                            inSimAnnealParams
                                                            targetEdge
                                                    else -- trace ("TBR: " <> (show $ minimum $ fmap snd5 rediagnosedGraphList))
                                                        pure (rediagnosedGraphList, Nothing)

                            -- simulated annealing/Drift stuff
                            -- based on steepest type swapping
                            Just simAnnealParams →
                                let -- to not overload paralle threads
                                    {-  This not so efficient is swapping in single graphs so leaving it be
                                    saRounds = if isNothing inSimAnnealParams then 1
                                               else rounds $ fromJust inSimAnnealParams

                                    (numGraphsToExamine, _) = divMod PU.getNumThreads saRounds -- this may not "drift" if finds alot better, but that's how its supposed to work
                                    -}
                                    numEdgesToExamine = graphsSteepest inGS -- min (graphsSteepest inGS) PU.getNumThreads
                                    firstSetEdges = take numEdgesToExamine edgesInPrunedGraph

                                    -- get heuristic delta joins for steepest edge set
                                    rerootEdgeList = filter ((/= prunedGraphRootIndex) . fst3) $ filter ((/= originalConnectionOfPruned) . fst3) firstSetEdges
                                in  do
                                        -- debugger "CASE OF -> Just" 0 splitGraphSimple
                                        -- True True to use IA fields and filter gaps
                                        rerootEdgeDataList ←
                                            getParallelChunkTraverse >>= \pTraverse →
                                                makeEdgeAction `pTraverse` rerootEdgeList

                                        rerootEdgeDeltaList' ←
                                            getParallelChunkTraverse >>= \pTraverse →
                                                joinAction `pTraverse` rerootEdgeDataList

                                        let rerootEdgeDeltaList = fmap (+ splitCost) rerootEdgeDeltaList'
                                        -- PU.seqParMap (parStrategy $ lazyParStrat inGS) (edgeJoinFunction charInfoVV targetEdgeData) rerootEdgeDataList

                                        let minDelta =
                                                if (not . null) rerootEdgeDeltaList
                                                    then minimum $ rerootEdgeDeltaList
                                                    else infinity
                                        let minEdgeList =
                                                if (not . null) rerootEdgeDeltaList
                                                    then fmap fst $ filter ((== minDelta) . snd) (zip rerootEdgeList rerootEdgeDeltaList)
                                                    else []

                                        -- check for possible better/equal graphs and verify
                                        rediagnosedGraphList ← case minEdgeList of
                                            [] → pure []
                                            xs →
                                                getParallelChunkTraverseBy strict2of5 >>= \pTraverse →
                                                    (reoptimizeAction . rerootAction) `pTraverse` xs

                                        let newMinCost =
                                                if (not . null) minEdgeList
                                                    then minimum $ fmap snd5 rediagnosedGraphList
                                                    else infinity

                                        -- only taking one for SA/Drift check
                                        let newMinGraph =
                                                if newMinCost /= infinity
                                                    then head $ filter ((== newMinCost) . snd5) rediagnosedGraphList
                                                    else emptyReducedPhylogeneticGraph

                                        -- if better always return it--hope this conditions short circuits so don't fully diagnose graph all the time
                                        if minDelta < curBestCost && newMinCost < curBestCost
                                            then do
                                                (_, newSAParams) ← U.simAnnealAccept inSimAnnealParams curBestCost newMinCost
                                                pure ([newMinGraph], newSAParams)
                                            else -- check if hit step limit--more for SA than drift

                                                if ((currentStep simAnnealParams) >= (numberSteps simAnnealParams))
                                                    || ((driftChanges simAnnealParams) >= (driftMaxChanges simAnnealParams))
                                                    then pure ([], inSimAnnealParams)
                                                    else do
                                                        (acceptGraph, newSAParams) ← U.simAnnealAccept inSimAnnealParams curBestCost minDelta
                                                        -- banner = if newMinCost <  curBestCost then "TBR heur better"
                                                        --         else "TBR Accepted not better"
                                                        -- if accepted (better or random) then return with updated annealing/Drift parameters
                                                        if acceptGraph
                                                            then -- trace (banner <> (show $ driftChanges $ fromJust newSAParams))
                                                                pure ([newMinGraph], newSAParams)
                                                            else -- rejected--recurse wirth updated params

                                                            -- trace ("TBR SA Drift Reject: " <> (show $ driftChanges $ fromJust newSAParams))

                                                                tbrJoin
                                                                    swapParams
                                                                    inGS
                                                                    inData
                                                                    splitGraph
                                                                    splitGraphSimple
                                                                    splitCost
                                                                    prunedGraphRootIndex
                                                                    originalConnectionOfPruned
                                                                    curBestCost
                                                                    (drop numEdgesToExamine edgesInPrunedGraph)
                                                                    newSAParams
                                                                    targetEdge


-- | rerootPrunedAndMakeGraph reroots the pruned graph component on the rerootEdge and joins to base gaph at target edge
rerootPrunedAndMakeGraph ∷ SimpleGraph → LG.Node → LG.Node → LG.LEdge EdgeInfo → LG.LEdge EdgeInfo → SimpleGraph
rerootPrunedAndMakeGraph splitGraphSimple prunedGraphRootIndex originalConnectionOfPruned (u, v, _) rerootEdge =
    -- get edges to delete and edges to add
    let (prunedEdgesToAdd, prunedEdgesToDelete) = getTBREdgeEditsSimple splitGraphSimple prunedGraphRootIndex rerootEdge

        -- edges to connect rerooted pruned component and base graph
        connectingEdges =
            [ (u, originalConnectionOfPruned, 0.0)
            , (originalConnectionOfPruned, v, 0.0)
            , (originalConnectionOfPruned, prunedGraphRootIndex, 0.0)
            ]

        tbrNewGraph =
            LG.insEdges (connectingEdges <> prunedEdgesToAdd) $
                LG.delEdges ([(u, v), (originalConnectionOfPruned, prunedGraphRootIndex)] <> prunedEdgesToDelete) splitGraphSimple
    in  tbrNewGraph


{- | getTBREdgeEditsSimple takes and edge and returns the list of edit to pruned subgraph
as a pair of edges to add and those to delete
since reroot edge is directed (e,v), edges away from v will have correct
orientation. Edges between 'e' and the root will have to be flipped
original root edges and reroort edge are deleted and new root and edge spanning orginal root created
delete original connection edge and creates a new one--like SPR
returns ([add], [delete])
-}
getTBREdgeEditsSimple ∷ SimpleGraph → LG.Node → LG.LEdge a → ([LG.LEdge Double], [LG.Edge])
getTBREdgeEditsSimple inGraph prunedGraphRootIndex rerootEdge =
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
                then (snd3 $ head originalRootEdges, snd3 $ last originalRootEdges, 0.0)
                else (snd3 $ last originalRootEdges, snd3 $ head originalRootEdges, 0.0)
        newRootEdges = [(prunedGraphRootIndex, fst3 rerootEdge, 0.0), (prunedGraphRootIndex, snd3 rerootEdge, 0.0)]
    in  -- assumes we are not checking original root
        -- rerooted
        -- delete orignal root edges and rerootEdge
        -- add new root edges
        -- and new edge on old root--but need orientation
        -- flip edges from new root to old (delete and add list)
       
        (newEdgeOnOldRoot : (flippedEdges <> newRootEdges), LG.toEdge rerootEdge : (fmap LG.toEdge (edgesToFlip <> originalRootEdges)))


-- | reoptimizeSplitGraphFromVertexTuple wrapper for reoptimizeSplitGraphFromVertex with last 3 args as tuple
reoptimizeSplitGraphFromVertexTuple
    ∷ GlobalSettings
    → ProcessedData
    → Bool
    → VertexCost
    → (DecoratedGraph, Int, Int)
    → PhyG (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertexTuple inGS inData doIA netPenaltyFactor (inSplitGraph, startVertex, prunedSubGraphRootVertex) =
    reoptimizeSplitGraphFromVertex inGS inData doIA netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex


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
reoptimizeSplitGraphFromVertex
    ∷ GlobalSettings
    → ProcessedData
    → Bool
    → VertexCost
    → DecoratedGraph
    → Int
    → Int
    → PhyG (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertex inGS inData doIA netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex =
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


getGraphCost ∷ [ReducedPhylogeneticGraph] → VertexCost
getGraphCost = snd5 . head


orInfinity ∷ ∀ a r t. (Foldable t, Fractional r) ⇒ (t a → r) → t a → r
orInfinity f xs
    | null xs = fromRational Real.infinity
    | otherwise = f xs

{- Old swap driving code that got overly complex and
    was not working properly--very slow
-}

{-
{- | swapSPRTBR performs SPR or TBR branch (edge) swapping on graphs
runs both SPR and TBR depending on argument since so much duplicated functionality
'steepest' abandons swap graph and switches to found graph as soon as anyhting 'better'
is found. The alternative (all) examines the entire neighborhood and retuns the best result
the return is a list of better graphs and the number of swapping rounds were required to ge there
if joinType ==  JoinAll is specified a single round is performed--otherwise a union rounds
alternate between joinPruned and joinAll.  This to be rapid but complete.
joinType = JoinAll for annealing/drifting
-}
swapSPRTBR
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Int
    → [ReducedPhylogeneticGraph]
    → [(Maybe SAParams, ReducedPhylogeneticGraph)]
    → PhyG ([ReducedPhylogeneticGraph], Int)
swapSPRTBR swapParams inGS inData inCounter currBestGraphs = \case
    [] → pure (currBestGraphs, inCounter)
    firstPair@(inSimAnnealParams, _) : _
        | joinType swapParams == JoinAll || isJust inSimAnnealParams →
            swapSPRTBR' (swapParams{joinType = JoinAll}) inGS inData inCounter firstPair

    firstPair@(inSimAnnealParams, inGraph) : morePairs → do
        -- join with union pruing first then followed by joinAll, but joinAlternate will return on better gaphs to return to join prune
        (firstList, firstCounter) ←
            swapSPRTBR' (swapParams{joinType = JoinPruned}) inGS inData inCounter firstPair
        -- the + 5 is to allow for extra buffer room with input graph and multiple equally costly solutions, can help
        bestFirstList ← GO.selectGraphs Best (keepNum swapParams) 0.0 $ inGraph : firstList
        (alternateList, alternateCounter') ←
            swapSPRTBRList (swapParams{joinType = JoinAlternate}) inGS inData firstCounter bestFirstList $
                zip
                    (U.generateUniqueRandList (length bestFirstList) inSimAnnealParams)
                    bestFirstList
        (bestAlternateList, alternateCounter) ← case joinType swapParams of
            JoinAlternate →
                let graphs = inGraph : (alternateList <> bestFirstList)
                    tag x = (x, alternateCounter')
                in  tag <$> GO.selectGraphs Best (keepNum swapParams) 0.0 graphs
            -- JoinPruned-don't alternate during search
            _ → pure (bestFirstList, firstCounter)

        -- recursive list version as opposed ot parMap version
        -- should reduce memory footprint at cost of less parallelism--but random replicates etc should take care of that
        (afterSecondList, afterSecondCounter) ←
            swapSPRTBRList (swapParams{joinType = JoinAll}) inGS inData alternateCounter bestAlternateList $
                zip
                    (U.generateUniqueRandList (length bestAlternateList) inSimAnnealParams)
                    bestAlternateList

        bestSecondList ← GO.selectGraphs Best (keepNum swapParams) 0.0 afterSecondList

        -- add final joinall if buffer full?
        let nextBestGraphs = case comparing getGraphCost bestSecondList currBestGraphs of
                LT → bestSecondList
                _ → currBestGraphs
        --liftIO $ putStr ("SSPR: " <> (show $ fmap snd5 nextBestGraphs) ) -- <> " " <> (show $ fmap (snd5 .snd) morePairs))
        swapSPRTBR swapParams inGS inData (afterSecondCounter + inCounter) nextBestGraphs morePairs

{- | swapSPRTBRList is a wrapper around swapSPRTBR' allowing for a list of graphs and a current best cost
reduce time of swap
-}
swapSPRTBRList
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Int
    → [ReducedPhylogeneticGraph]
    → [(Maybe SAParams, ReducedPhylogeneticGraph)]
    → PhyG ([ReducedPhylogeneticGraph], Int)
swapSPRTBRList swapParams inGS inData inCounter curBestGraphs = \case
    [] → pure (curBestGraphs, inCounter)
    -- currrent graph worse than best saved
    (_, inGraph) : otherDoubles
        | snd5 inGraph > getGraphCost curBestGraphs →
            swapSPRTBRList swapParams inGS inData inCounter curBestGraphs otherDoubles

    firstPair@(inSimAnnealParams, inGraph) : otherDoubles → do
        --logWith LogInfo "In SwapSPRTBRList"

        (graphList, swapCounter) ← swapSPRTBR' swapParams inGS inData inCounter firstPair
        bestNewGraphList ← GO.selectGraphs Best (keepNum swapParams) 0.0 graphList
        --liftIO $ putStr (" SSPRTBRList: " <> (show $ fmap snd5 (inGraph : bestNewGraphList)))
        let recurse = swapSPRTBRList swapParams inGS inData swapCounter
        let recurse' = flip recurse otherDoubles
        case comparing getGraphCost bestNewGraphList curBestGraphs of
            -- found equal
            EQ → recurse' =<< GO.selectGraphs Unique (keepNum swapParams) 0.0 (curBestGraphs <> bestNewGraphList)
            -- found worse
            GT → recurse' curBestGraphs
            -- found better
            LT →
                pure (bestNewGraphList, swapCounter)
                {-
                let betterResult =
                        let doubleToSwap = zip (U.generateUniqueRandList (length bestNewGraphList) inSimAnnealParams) bestNewGraphList
                        in  recurse bestNewGraphList $ doubleToSwap <> otherDoubles
                in  case joinType swapParams of
                        JoinAlternate → betterResult
                        JoinPruned → betterResult
                        _ → recurse' bestNewGraphList
                -}

{- | swapSPRTBR' is the central functionality of swapping allowing for repeated calls with alternate
options such as joinType to ensure complete swap but with an edge unions pass to
reduce time of swap
also manages SA/Drifting versus greedy swap
-}
swapSPRTBR'
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Int
    → (Maybe SAParams, ReducedPhylogeneticGraph)
    → PhyG ([ReducedPhylogeneticGraph], Int)
swapSPRTBR' swapParams inGS inData@(leafNames, _, _) inCounter (inSimAnnealParams, inGraph@(simpleG, costG, _, _, _))
    | LG.isEmpty simpleG = pure ([], 0)
    | otherwise = do
        --liftIO $ putStr (" SSPR': " <> (show costG))
        --logWith LogInfo "In SwapSPRTBR'"
        -- inGraphNetPenalty <- POSW.getNetPenaltyReduced inGS inData inGraph
        inGraphNetPenalty ← T.getPenaltyFactor inGS inData Nothing $ GO.convertReduced2PhylogeneticGraph inGraph
        let inGraphNetPenaltyFactor = inGraphNetPenalty / costG
        let numLeaves = V.length leafNames
        case inSimAnnealParams of
            Nothing → do
                -- steepest takes immediate best--does not keep equall cost-- for now--disabled not working correctly so goes to "all"
                (swappedGraphs, counter, _) ←
                    swapAll
                        swapParams
                        inGS
                        inData
                        inCounter
                        costG
                        []
                        [inGraph]
                        numLeaves
                        inGraphNetPenaltyFactor
                        Nothing
                pure $ case swappedGraphs of
                    [] → ([inGraph], counter)
                    gs → (gs, counter)
            -- simulated annealing/drifting acceptance does a steepest with SA acceptance
            Just simAnneal →
                -- then a swap steepest and all on annealed graph
                -- same at this level method (SA, Drift) choice occurs at lower level
                -- annealed should only yield a single graph
                let -- create list of params with unique list of random values for rounds of annealing
                    annealDriftRounds = rounds simAnneal
                    newSimAnnealParamList = U.generateUniqueRandList annealDriftRounds inSimAnnealParams
                    -- parallel setup
                    action ∷ Maybe SAParams → PhyG ([ReducedPhylogeneticGraph], Int, Maybe SAParams)
                    action = swapAll swapParams inGS inData 0 costG [] [inGraph] numLeaves inGraphNetPenaltyFactor
                in  -- this to ensure current step set to 0
                    do
                        swapPar ← getParallelChunkTraverse
                        (annealDriftGraphs', anealDriftCounterList, _) ← unzip3 <$> swapPar action newSimAnnealParamList

                        -- annealed/Drifted 'mutated' graphs
                        annealDriftGraphs ← GO.selectGraphs Unique (keepNum swapParams) 0.0 $ concat annealDriftGraphs'

                        -- swap back "normally" if desired for full drifting/annealing
                        (swappedGraphs, counter, _) ←
                            swapAll
                                swapParams
                                inGS
                                inData
                                (sum anealDriftCounterList)
                                (min costG (minimum $ snd5 <$> annealDriftGraphs))
                                []
                                annealDriftGraphs
                                numLeaves
                                inGraphNetPenaltyFactor
                                Nothing

                        bestGraphs ← GO.selectGraphs Best (keepNum swapParams) 0.0 $ inGraph : swappedGraphs
                        -- this Bool for Genetic Algorithm mutation step
                        pure $
                            if not $ returnMutated swapParams
                                then (bestGraphs, counter)
                                else (annealDriftGraphs, sum anealDriftCounterList)
-}
