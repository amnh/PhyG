{- |
Module specifying graph egde adding and deleting functions
-}
module Search.NetworkAddDeleteV2 (
    deleteAllNetEdges,
    insertAllNetEdges,
    moveAllNetEdges,
    deltaPenaltyAdjustment,
    --deleteNetEdge,
    --deleteOneNetAddAll,
    addDeleteNetEdges,
    getCharacterAddDelta,
    getBlockDeltaAdd,
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

{-----------------------------------
    Move network edge functions
------------------------------------}

-- | moveAllNetEdges is a wrapper for moveAllNetEdges' allowing for multiple simulated annealing rounds
moveAllNetEdges
    ∷ GlobalSettings
    → ProcessedData
    -> NetParams
    → Int
    → ([ReducedPhylogeneticGraph], VertexCost)
    → (Maybe SAParams, [ReducedPhylogeneticGraph])
    → PhyG ([ReducedPhylogeneticGraph], Int)
moveAllNetEdges inGS inData netParams counter (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) = case inSimAnnealParams of
    Nothing →
        moveAllNetEdges'
            inGS
            inData
            netParams 
            counter
            (curBestGraphList, curBestGraphCost)
            Nothing
            inPhyloGraphList
    Just simAnneal →
        let -- create list of params with unique list of random values for rounds of annealing
            annealingRounds = rounds simAnneal
            saParamList = replicate annealingRounds inSimAnnealParams

            -- parallel setup
            action ∷ (Maybe SAParams, [ReducedPhylogeneticGraph]) → PhyG ([ReducedPhylogeneticGraph], Int)
            action =
                moveAllNetEdges''
                    inGS
                    inData
                    netParams
                    counter
                    (curBestGraphList, curBestGraphCost)
        in  do
                moveResult ←
                    getParallelChunkTraverse >>= \pTraverse →
                        pTraverse action . zip saParamList $ replicate annealingRounds inPhyloGraphList
                let (annealRoundsList, counterList) = unzip moveResult
                GO.selectGraphs Best (outgroupIndex inGS) (netKeepNum netParams) 0 (fold annealRoundsList) <&> \x → (x, sum counterList)

-- | moveAllNetEdges'' is wrapper around moveAllNetEdges' to use parmap
moveAllNetEdges''
    ∷ GlobalSettings
    → ProcessedData
    -> NetParams
    → Int
    → ([ReducedPhylogeneticGraph], VertexCost)
    → (Maybe SAParams, [ReducedPhylogeneticGraph])
    → PhyG ([ReducedPhylogeneticGraph], Int)
moveAllNetEdges'' inGS inData netParams counter (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
    moveAllNetEdges'
        inGS
        inData
        netParams
        counter
        (curBestGraphList, curBestGraphCost)
        inSimAnnealParams
        inPhyloGraphList


{- | moveAllNetEdges' removes each edge and adds an edge to all possible places (or steepest) each round
until no better or additional graphs are found
call with ([], infinity) [single input graph]
-}
moveAllNetEdges'
    ∷ GlobalSettings
    → ProcessedData
    -> NetParams
    → Int
    → ([ReducedPhylogeneticGraph], VertexCost)
    → Maybe SAParams
    → [ReducedPhylogeneticGraph]
    → PhyG ([ReducedPhylogeneticGraph], Int)
moveAllNetEdges' inGS inData netParams counter (curBestGraphList, curBestGraphCost) inSimAnnealParams = \case
    [] → pure (take (netKeepNum netParams) curBestGraphList, counter)
    firstPhyloGraph : otherPhyloGraphs
        | LG.isEmpty $ fst5 firstPhyloGraph →
            moveAllNetEdges'
                inGS
                inData
                netParams
                counter
                (curBestGraphList, curBestGraphCost)
                inSimAnnealParams
                otherPhyloGraphs
    firstPhyloGraph : otherPhyloGraphs →
        let currentCost = min curBestGraphCost $ snd5 firstPhyloGraph
            -- parallel setup
            action ∷ LG.Edge → PhyG [ReducedPhylogeneticGraph]
            action = deleteOneNetAddAll' inGS inData netParams firstPhyloGraph inSimAnnealParams
        in  do
                -- randomize order of edges to try moving
                netEdgeList ←
                    let edges = LG.labNetEdges $ thd5 firstPhyloGraph
                        permutationOf
                            | (netRandom netParams) = shuffleList
                            | otherwise = pure
                     in permutationOf edges

                deleteResult ←
                    getParallelChunkTraverse >>= \pTraverse →
                        pTraverse (action . LG.toEdge) netEdgeList

                let newGraphList' = fold deleteResult
                newGraphList ← GO.selectGraphs Best (outgroupIndex inGS) (netKeepNum netParams) 0 newGraphList'
                let newGraphCost = case newGraphList' of
                        [] → infinity
                        (_, c, _, _, _) : _ → c

                --logWith LogInfo (" MoveAll': " <> (show newGraphCost) <> " from " <> (show $ fmap snd5 $ concat deleteResult) <>"\n")
                -- if graph is a tree no edges to delete
                case netEdgeList of
                    [] → (firstPhyloGraph : otherPhyloGraphs, counter) <$ logWith LogInfo "\t\tGraph in move has no network edges to move--skipping\n"
                    e : es → case inSimAnnealParams of
                        Nothing → case newGraphCost `compare` currentCost of
                            GT →
                                moveAllNetEdges'
                                    inGS
                                    inData
                                    netParams
                                    (counter + 1)
                                    (firstPhyloGraph : curBestGraphList, currentCost)
                                    inSimAnnealParams
                                    otherPhyloGraphs
                            LT
                                | (netSteepest netParams) →
                                    moveAllNetEdges'
                                        inGS
                                        inData
                                        netParams
                                        (counter + 1)
                                        (newGraphList, newGraphCost)
                                        inSimAnnealParams
                                        newGraphList
                            LT →
                                moveAllNetEdges'
                                    inGS
                                    inData
                                    netParams
                                    (counter + 1)
                                    (newGraphList, newGraphCost)
                                    inSimAnnealParams
                                    (newGraphList <> otherPhyloGraphs)
                            EQ → do
                                -- new graph list contains the input graph if equal and filterd unique already in moveAllNetEdges
                                newCurSameBestList ← GO.selectGraphs Unique (outgroupIndex inGS) (netKeepNum netParams) 0 $ curBestGraphList <> newGraphList
                                moveAllNetEdges'
                                    inGS
                                    inData
                                    netParams
                                    (counter + 1)
                                    (newCurSameBestList, currentCost)
                                    inSimAnnealParams
                                    otherPhyloGraphs

                        -- sim anneal choice
                        Just simAnneal →
                            -- not sure why this was here--seems to work
                            let -- abstract stopping criterion to continue
                                (numDone, numMax) = case method simAnneal of
                                    SimAnneal → currentStep &&& numberSteps $ simAnneal
                                    _ → driftChanges &&& driftMaxChanges $ simAnneal
                            in  do
                                    uniqueGraphList ← GO.selectGraphs Unique (outgroupIndex inGS) (netKeepNum netParams) 0 newGraphList'

                                    let (annealBestCost, nextUniqueList) = case uniqueGraphList of
                                            [] → (curBestGraphCost, [])
                                            (_, cost, _, _, _) : more → (min curBestGraphCost cost, more)

                                    (acceptFirstGraph, newSAParams) ← case uniqueGraphList of
                                        [] → pure (False, U.incrementSimAnnealParams inSimAnnealParams)
                                        (_, c, _, _, _) : _ → U.simAnnealAccept inSimAnnealParams annealBestCost c

                                    case numDone `compare` numMax of
                                        LT | acceptFirstGraph → do
                                            moveAllNetEdges'
                                                inGS
                                                inData
                                                netParams
                                                (counter + 1)
                                                ((head uniqueGraphList) : curBestGraphList, annealBestCost)
                                                newSAParams
                                                (nextUniqueList <> otherPhyloGraphs)
                                        LT → do
                                            moveAllNetEdges'
                                                inGS
                                                inData
                                                netParams
                                                (counter + 1)
                                                (curBestGraphList, annealBestCost)
                                                newSAParams
                                                (nextUniqueList <> otherPhyloGraphs)

                                        -- if want non-optimized list for GA or whatever
                                        _ | (netReturnMutated netParams)  → pure (take (netKeepNum netParams) curBestGraphList, counter)
                                        -- optimize list and return
                                        _ → do
                                            (bestMoveList', counter') ←
                                                moveAllNetEdges'
                                                    inGS
                                                    inData
                                                    netParams {netReturnMutated = False}
                                                    (counter + 1)
                                                    ([], annealBestCost)
                                                    Nothing
                                                    $ take (netKeepNum netParams) curBestGraphList
                                            bestMoveList ← GO.selectGraphs Best (outgroupIndex inGS) (netKeepNum netParams) 0 bestMoveList'
                                            pure (take (netKeepNum netParams) bestMoveList, counter')

-- | deleteOneNetAddAll' wrapper on deleteOneNetAddAll to allow for parmap
deleteOneNetAddAll'
    ∷ GlobalSettings
    → ProcessedData
    -> NetParams
    → ReducedPhylogeneticGraph
    → Maybe SAParams
    → LG.Edge
    → PhyG [ReducedPhylogeneticGraph]
deleteOneNetAddAll' inGS inData netParams inPhyloGraph inSimAnnealParams edgeToDelete =
    deleteOneNetAddAll
        inGS
        inData
        netParams
        inPhyloGraph
        [edgeToDelete]
        inSimAnnealParams

{- | deleteOneNetAddAll version deletes net edges in turn and reAdds-based on original cost
but this cost in graph (really not correct) but allows logic of insert edge to function better
unlike deleteOneNetAddAll' only deals with single edge deletion at a time
-}
deleteOneNetAddAll
    ∷ GlobalSettings
    → ProcessedData
    -> NetParams
    → ReducedPhylogeneticGraph
    → [LG.Edge]
    → Maybe SAParams
    → PhyG [ReducedPhylogeneticGraph]
deleteOneNetAddAll inGS inData netParams inPhyloGraph edgeToDeleteList inSimAnnealParams =
    if null edgeToDeleteList
        then do
            --logWith LogInfo ("\tGraph has no edges to move---skipping")
            pure [inPhyloGraph]
        else
            if LG.isEmpty $ thd5 inPhyloGraph
                then error "Empty graph in deleteOneNetAddAll"
                else do
                    -- trace ("DONAA-New: " <> (show $ snd5 inPhyloGraph) <> " Steepest:" <> (show doSteepest)) (
                    logWith
                        LogInfo
                        ("Moving " <> (show $ length edgeToDeleteList) <> " network edges, current best cost: " <> (show $ snd5 inPhyloGraph) <> "\n")
                    -- start with initial graph cost
                    let inGraphCost = snd5 inPhyloGraph

                    -- get deleted simple graphs and bool for changed
                    delGraphBoolPair ← deleteNetworkEdge (fst5 inPhyloGraph) (head edgeToDeleteList)

                    -- no change in network structure
                    if snd delGraphBoolPair == False
                        then do
                            deleteOneNetAddAll
                                inGS
                                inData
                                netParams
                                inPhyloGraph
                                (tail edgeToDeleteList)
                                inSimAnnealParams
                        else
                            let simpleGraphToInsert = fst delGraphBoolPair

                                (_, _, _, curNetNodes) = LG.splitVertexList simpleGraphToInsert
                                curNumNetNodes = length curNetNodes

                                -- optimize deleted graph and update cost with input cost
                                leafGraph = LG.extractLeafGraph $ thd5 inPhyloGraph
                            in  do
                                    graphToInsert ← T.multiTraverseFullyLabelSoftWiredReduced inGS inData False False leafGraph Nothing simpleGraphToInsert -- `using` PU.myParListChunkRDS

                                    -- keep same cost and just keep better--check if better than original later
                                    let graphToInsert' = T.updatePhylogeneticGraphCostReduced graphToInsert inGraphCost

                                    if isNothing inSimAnnealParams then 
                                        do

                                        insertedGraphTripleList ←
                                            insertEachNetEdgeHeuristicGather --insertEachNetEdge
                                                inGS
                                                inData
                                                netParams {netMaxEdges = curNumNetNodes + 1}
                                                Nothing
                                                inSimAnnealParams
                                                graphToInsert'

                                        let newMinimumCost = snd3 insertedGraphTripleList
                                        let newBestGraphs = filter ((== newMinimumCost) . snd5) $ fst3 insertedGraphTripleList

                                        --logWith LogInfo  ("DONAA-New: " <> (show (inGraphCost, snd5 graphToInsert, fmap snd5 $ fst3 insertedGraphTripleList, newMinimumCost)))
                                        if newMinimumCost < inGraphCost
                                            then do
                                                --logWith LogInfo ("DONA=> " <> (show newMinimumCost))
                                                pure newBestGraphs
                                            else do
                                                deleteOneNetAddAll
                                                    inGS
                                                    inData
                                                    netParams
                                                    inPhyloGraph
                                                    (tail edgeToDeleteList)
                                                    inSimAnnealParams

                                    -- SA/Drift stuff--single call recurses within insertRoundsSA
                                    else do
                                        (saGraphList, _) <- insertRoundsSA inGS inData (netParams {netEditType = NetMove}) 0 ([inPhyloGraph], snd5 inPhyloGraph) inSimAnnealParams  [graphToInsert']
                                        pure saGraphList
                                         



{---------------------------------------
    Add-Delete network edge functions
----------------------------------------}

{- |
'addDeleteNetEdges is a wrapper for 'addDeleteNetEdges'' allowing for multiple simulated annealing rounds.
-}
addDeleteNetEdges
    ∷ GlobalSettings
    → ProcessedData
    -> NetParams
    → Int
    → Int
    → ([ReducedPhylogeneticGraph], VertexCost)
    → (Maybe SAParams, [ReducedPhylogeneticGraph])
    → PhyG ([ReducedPhylogeneticGraph], Int)
addDeleteNetEdges inGS inData netParams maxRounds counter (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) = case inSimAnnealParams of
    Nothing →
        addDeleteNetEdges'
            inGS
            inData
            netParams
            maxRounds
            counter
            (curBestGraphList, curBestGraphCost)
            Nothing
            inPhyloGraphList
    Just simAnneal →
        let -- create list of params with unique list of random values for rounds of annealing
            annealingRounds = rounds simAnneal
            saParamList = replicate annealingRounds inSimAnnealParams

            -- parallel setup
            action ∷ (Maybe SAParams, [ReducedPhylogeneticGraph]) → PhyG ([ReducedPhylogeneticGraph], Int)
            action =
                addDeleteNetEdges''
                    inGS
                    inData
                    netParams 
                    maxRounds
                    counter
                    (curBestGraphList, curBestGraphCost)
        in  do
                addDeleteResult ←
                    getParallelChunkTraverse >>= \pTraverse →
                        pTraverse action . zip saParamList $ replicate annealingRounds inPhyloGraphList
                let (annealRoundsList, counterList) = unzip addDeleteResult
                GO.selectGraphs Best (outgroupIndex inGS) (netKeepNum netParams) 0 (fold annealRoundsList) <&> \x → (x, sum counterList)


-- | addDeleteNetEdges'' is wrapper around addDeleteNetEdges' to use parallel above
addDeleteNetEdges''
    ∷ GlobalSettings
    → ProcessedData
    -> NetParams
    → Int
    → Int
    → ([ReducedPhylogeneticGraph], VertexCost)
    → (Maybe SAParams, [ReducedPhylogeneticGraph])
    → PhyG ([ReducedPhylogeneticGraph], Int)
addDeleteNetEdges'' inGS inData netParams maxRounds counter (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
    addDeleteNetEdges'
        inGS
        inData
        netParams
        maxRounds
        counter
        (curBestGraphList, curBestGraphCost)
        inSimAnnealParams
        inPhyloGraphList


{- | addDeleteNetEdges' alternates addition and deletion of network edges
until no better or additional graphs are found (or max rounds met)
call with ([], infinity) [single input graph]
-}
addDeleteNetEdges'
    ∷ GlobalSettings
    → ProcessedData
    -> NetParams
    → Int
    → Int
    → ([ReducedPhylogeneticGraph], VertexCost)
    → Maybe SAParams
    → [ReducedPhylogeneticGraph]
    → PhyG ([ReducedPhylogeneticGraph], Int)
addDeleteNetEdges' inGS inData netParams maxRounds counter (curBestGraphList, curBestGraphCost) inSimAnnealParams = \case
    [] → pure (take (netKeepNum netParams) curBestGraphList, counter)
    -- if hit maxmimum rounds then return
    inPhyloGraphList → case counter `compare` maxRounds of
        EQ → pure (take (netKeepNum netParams) curBestGraphList, counter)
        -- other wise add/delete
        _ → do
            -- insert edges first
            (insertGraphList, _) ←
                insertAllNetEdges'
                    inGS
                    inData
                    netParams
                    counter
                    (curBestGraphList, curBestGraphCost)
                    inSimAnnealParams
                    inPhyloGraphList

            -- this to update randlists in SAPArams for subsequent calls
            let updatedSAParamList = case inSimAnnealParams of
                    Nothing → [Nothing, Nothing]
                    _ → replicate 2 inSimAnnealParams

            -- if no better--take input for delte phase
            (insertGraphList', insertGraphCost, toDeleteList) ← case insertGraphList of
                [] → pure (curBestGraphList, curBestGraphCost, inPhyloGraphList)
                gs → do
                    newList ← GO.selectGraphs Best (outgroupIndex inGS) (maxBound ∷ Int) 0 gs
                    pure (newList, snd5 $ head newList, newList)

            -- delete edges
            (deleteGraphList, _) ←
                deleteAllNetEdges'
                    inGS
                    inData
                    netParams
                    counter
                    (insertGraphList', insertGraphCost)
                    (head updatedSAParamList)
                    toDeleteList

            -- gather beter if any
            (newBestGraphList, newBestGraphCost, graphsToDoNext) ← case deleteGraphList of
                [] → pure (curBestGraphList, curBestGraphCost, inPhyloGraphList)
                gs → do
                    newDeleteGraphs ← GO.selectGraphs Best (outgroupIndex inGS) (maxBound ∷ Int) 0 gs
                    pure (newDeleteGraphs, snd5 $ head newDeleteGraphs, newDeleteGraphs)

            -- check is same then return
            case newBestGraphCost `compare` curBestGraphCost of
                EQ → pure (take (netKeepNum netParams) curBestGraphList, counter)
                -- if better (or nothing) keep going
                _ →
                    addDeleteNetEdges'
                        inGS
                        inData
                        netParams
                        maxRounds
                        (counter + 1)
                        (newBestGraphList, newBestGraphCost)
                        (last updatedSAParamList)
                        graphsToDoNext


{-----------------------------------
    Delete network edge functions
------------------------------------}

-- | deleteAllNetEdges is a wrapper for moveAllNetEdges' allowing for multiple simulated annealing rounds
deleteAllNetEdges
    ∷ GlobalSettings
    → ProcessedData
    -> NetParams
    → Int
    → ([ReducedPhylogeneticGraph], VertexCost)
    → (Maybe SAParams, [ReducedPhylogeneticGraph])
    → PhyG ([ReducedPhylogeneticGraph], Int)
deleteAllNetEdges inGS inData netParams counter (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
    
    deleteAllNetEdges'
                inGS
                inData
                netParams
                counter
                (curBestGraphList, curBestGraphCost)
                inSimAnnealParams
                inPhyloGraphList


{- | deleteAllNetEdges deletes network edges one each each round until no better or additional
graphs are found
call with ([], infinity) [single input graph]
-}
deleteAllNetEdges'
    ∷ GlobalSettings
    → ProcessedData
    -> NetParams
    → Int
    → ([ReducedPhylogeneticGraph], VertexCost)
    → Maybe SAParams
    → [ReducedPhylogeneticGraph]
    → PhyG ([ReducedPhylogeneticGraph], Int)
deleteAllNetEdges' inGS inData netParams counter (curBestGraphList, curBestGraphCost) inSimAnnealParams inPhyloGraphList = 
    if null inPhyloGraphList then 
        pure (take (netKeepNum netParams) curBestGraphList, counter)

    else 
        let firstPhyloGraph : otherPhyloGraphs = inPhyloGraphList
            currentCost = min curBestGraphCost $ snd5 firstPhyloGraph
            netNodes = fth4 $ LG.splitVertexList $ thd5 firstPhyloGraph
        in
        -- if a tree then no edges to delete and recurse
        if length netNodes == 0 then do
            logWith LogInfo "\tGraph in delete network edges is tree--skipping\n"
            deleteAllNetEdges'
                inGS
                inData
                netParams
                (counter + 1)
                (firstPhyloGraph : curBestGraphList, currentCost)
                inSimAnnealParams
                otherPhyloGraphs

        else if isNothing inSimAnnealParams then 
            do 
                logWith LogInfo $ "\tNumber of network edges: " <> (show $ length netNodes) <> " Number of graphs: " <> (show $ length inPhyloGraphList) <> "\n"
                (newGraphList', _, newSAParams) ←
                    deleteEachNetEdge
                        inGS
                        inData
                        netParams 
                        False
                        inSimAnnealParams
                        firstPhyloGraph

                newGraphList ← GO.selectGraphs Best (outgroupIndex inGS) (netKeepNum netParams) 0 newGraphList'
                let newGraphCost = case newGraphList of
                        [] → infinity
                        (_, c, _, _, _) : _ → c

                if null newGraphList then
                    pure (take (netKeepNum netParams) curBestGraphList, counter + 1)

                else 
                    postProcessNetworkDelete
                            inGS
                            inData
                            netParams
                            counter
                            (curBestGraphList, curBestGraphCost)
                            inSimAnnealParams
                            inPhyloGraphList
                            newGraphList
                            newGraphCost
                            currentCost

        -- SA/Drift
        else 
            deleteRoundsSA inGS inData netParams counter (curBestGraphList, curBestGraphCost) inSimAnnealParams inPhyloGraphList

{- | deleteRoundsSA deletes network edges one each each round using Simulated Annealing/Drift

    Randomizes graphs to examine in deleteEdge so can do multiple rounds there
    After SA delete can return mutated or do a addEdges to get back to "optimal" edges
-}
deleteRoundsSA
    ∷ GlobalSettings
    → ProcessedData
    → NetParams
    → Int
    → ([ReducedPhylogeneticGraph], VertexCost)
    → Maybe SAParams
    → [ReducedPhylogeneticGraph]
    → PhyG ([ReducedPhylogeneticGraph], Int)
deleteRoundsSA inGS inData netParams counter (curBestGraphList, curBestGraphCost) inSimAnnealParams inPhyloGraphList = 
    let -- create list of params with unique list of random values for rounds of annealing
        annealingRounds = rounds $ fromJust inSimAnnealParams
        annealParamList = replicate annealingRounds inSimAnnealParams

        -- set up parallel
        deleteAction :: (Maybe SAParams, ReducedPhylogeneticGraph) -> PhyG ([ReducedPhylogeneticGraph], VertexCost, Maybe SAParams)
        deleteAction = deleteEachNetEdge' inGS inData netParams False
        
    in do
        graphsToAnneal <- if length inPhyloGraphList == 1 then 
                                pure $ replicate annealingRounds $ head  inPhyloGraphList
                          else do
                                shuffledLists <- mapM shuffleList $ replicate annealingRounds inPhyloGraphList
                                pure $ fmap head shuffledLists

        deleteActionPar <- getParallelChunkTraverse 
        newGraphTripleList <- deleteActionPar deleteAction (zip annealParamList graphsToAnneal)

        uniqueList ← GO.selectGraphs Unique (outgroupIndex inGS) (netKeepNum netParams) 0.0 (concat $ fmap fst3 newGraphTripleList)

        -- return mutated for GA
        if netReturnMutated netParams then do
           pure (uniqueList, counter + annealingRounds)

        else do
            logWith LogInfo "\t\tAdding edges after network deletion Simulated Annealing/Drifing\n"
            -- delete edges to get back to "best" edge lists
            --delActionPar <- getParallelChunkTraverse
            --deletedTripleList <- delActionPar deleteAction uniqueList

            -- set steapest to false so tries all resulting graphs
            let uniqueCost = if (not .null) uniqueList then
                                minimum $ fmap snd5 uniqueList
                          else infinity

            -- if found nothing
            if null uniqueList then 
                pure (inPhyloGraphList, counter)

            else do
                (insertedGraphList, _)  <- insertAllNetEdges' inGS inData (netParams {netSteepest = False}) counter (uniqueList, uniqueCost) Nothing uniqueList 

                -- return better and equal including inputs
                finalList ← GO.selectGraphs Best (outgroupIndex inGS) (netKeepNum netParams) 0.0 $ inPhyloGraphList <> insertedGraphList -- (concat $ fmap fst3 deletedTripleList)

                let netNodesList = fmap length $ fmap (fth4 . LG.splitVertexList . thd5) finalList

                let finalCost = if (not .null) finalList then
                                    minimum $ fmap snd5 finalList
                              else infinity            

                if finalCost < uniqueCost then
                    logWith LogInfo ("\t\t-> " <> (show $ min finalCost uniqueCost) <> " with " <> (show netNodesList) <> " netWork nodes\n")
                else 
                    logWith LogInfo ""
                
                pure (finalList, counter + annealingRounds)

{- | deleteEachNetEdge' is a wrapper around deleteEachNetEdge to allow for zipping new simAnneal params for each
replicate
-}
deleteEachNetEdge'
    ∷ GlobalSettings
    → ProcessedData
    -> NetParams
    → Bool
    → (Maybe SAParams, ReducedPhylogeneticGraph)
    → PhyG ([ReducedPhylogeneticGraph], VertexCost, Maybe SAParams)
deleteEachNetEdge' inGS inData netParams force (inSimAnnealParams, inPhyloGraph) =
    deleteEachNetEdge inGS inData netParams force inSimAnnealParams inPhyloGraph

{- | deleteEachNetEdge takes a phylogenetic graph and deletes all network edges one at time
and returns best list of new Phylogenetic Graphs and cost
even if worse--used for simulated annealing later
if equal returns unique graph list
-}
deleteEachNetEdge
    ∷ GlobalSettings
    → ProcessedData
    -> NetParams
    → Bool
    → Maybe SAParams
    → ReducedPhylogeneticGraph
    → PhyG ([ReducedPhylogeneticGraph], VertexCost, Maybe SAParams)
deleteEachNetEdge inGS inData netParams force inSimAnnealParams inPhyloGraph =
    let currentCost = snd5 inPhyloGraph
    in

    -- empty graph
    if LG.isEmpty $ thd5 inPhyloGraph then
        pure ([], infinity, inSimAnnealParams) -- error "Empty input phylogenetic graph in deleteAllNetEdges"

    -- SA/Drif and done all steps
    else if isJust inSimAnnealParams && U.isSimAnnealTerminated inSimAnnealParams then 
        pure ([inPhyloGraph], currentCost, inSimAnnealParams)

    -- stuff to do
    else
            let currentCost = snd5 inPhyloGraph

                -- potentially randomize order of list
                networkEdgeList = LG.netEdges $ thd5 inPhyloGraph

                --- parallel
                action ∷ LG.Edge → PhyG ReducedPhylogeneticGraph
                action = deleteNetEdge inGS inData inPhyloGraph force inSimAnnealParams
            in  do
                if null networkEdgeList
                    then do
                        -- for SA/Drift if deletes all
                        if (isJust inSimAnnealParams) then
                            pure ([inPhyloGraph], snd5 inPhyloGraph, inSimAnnealParams)
                        else do
                            logWith LogInfo ("\tNo network edges to delete" <> "\n")
                            pure ([], infinity, inSimAnnealParams)
                else do
                    logWith LogInfo ("\tNetwork edges to delete: " <> (show $ length networkEdgeList) <> "\n")
                    --could shuffle edge list if not doain all at once--but are now
                    delNetEdgeList ←
                        getParallelChunkTraverse >>= \pTraverse →
                            action `pTraverse` networkEdgeList

                    -- list shoul;d always have graphs since deleting always yields a valid graph
                    let minCost =
                                if null delNetEdgeList
                                    then infinity
                                    else minimum $ fmap snd5 delNetEdgeList

                    let (newGraphList, newSAParams) =  (delNetEdgeList, U.incrementSimAnnealParams inSimAnnealParams)
                        
                    -- Always return best if better
                    if isNothing inSimAnnealParams || minCost < (snd5 inPhyloGraph) then do

                        bestCostGraphList ← filter ((/= infinity) . snd5) <$> GO.selectGraphs Best (outgroupIndex inGS) (netKeepNum netParams) 0 (inPhyloGraph : newGraphList)
                        

                        --logWith LogInfo ("\tNew costs:" <> (show $ fmap snd5 newGraphList) <>"\n")
                        pure (bestCostGraphList, minCost, newSAParams)

                    -- SA/Drift and Worse or equal
                    else 
                        do
                            newGraphList ← filter ((/= infinity) . snd5) <$> GO.selectGraphs AtRandom (outgroupIndex inGS) (netKeepNum netParams) 0 newGraphList

                            let bestNewGraph = head newGraphList

                            (acceptGraph, newSAParams) <- U.simAnnealAccept inSimAnnealParams currentCost minCost

                            if acceptGraph then
                                deleteEachNetEdge inGS inData netParams force newSAParams bestNewGraph

                            else  
                                deleteEachNetEdge inGS inData netParams force newSAParams inPhyloGraph
                        
                    
                    
{- | deleteEdge deletes an edge (checking if network) and rediagnoses graph
contacts in=out=1 edgfes and removes node, reindexing nodes and edges
naive for now
force requires reoptimization no matter what--used for net move
skipping heuristics for now--awful
calls deleteNetworkEdge that has various graph checks
-}
deleteNetEdge
    ∷ GlobalSettings
    → ProcessedData
    → ReducedPhylogeneticGraph
    → Bool
    -> Maybe SAParams
    → LG.Edge
    → PhyG ReducedPhylogeneticGraph
deleteNetEdge inGS inData inPhyloGraph force inSimAnnealParams edgeToDelete =
    if LG.isEmpty $ thd5 inPhyloGraph
        then error "Empty input phylogenetic graph in deleteNetEdge"
        else
            if not (LG.isNetworkEdge (fst5 inPhyloGraph) edgeToDelete)
                then error ("Edge to delete: " <> (show edgeToDelete) <> " not in graph:\n" <> (LG.prettify $ fst5 inPhyloGraph))
                else do
                    -- trace ("DNE: " <> (show edgeToDelete)) (
                    (delSimple, wasModified) ← deleteNetworkEdge (fst5 inPhyloGraph) edgeToDelete

                    -- delSimple = GO.contractIn1Out1EdgesRename $ LG.delEdge edgeToDelete $ fst5 inPhyloGraph

                    -- prune other edges if now unused
                    let pruneEdges = False

                    -- don't warn that edges are being pruned
                    let warnPruneEdges = False

                    -- graph optimization from root
                    let startVertex = Nothing

                    heuristicDelta <- heuristicDeleteDelta inGS inPhyloGraph edgeToDelete

                    let edgeDeleteDelta = deltaPenaltyAdjustment inGS inPhyloGraph "delete"

                    let totalHeuristicCost = (snd5 inPhyloGraph) + heuristicDelta - edgeDeleteDelta

                    -- full two-pass optimization--cycles checked in edge deletion function
                    let leafGraph = LG.extractLeafGraph $ thd5 inPhyloGraph

                    {-
                    newPhyloGraph ←
                        -- check if deletion modified graph
                        if not wasModified then 
                            pure emptyReducedPhylogeneticGraph

                        else if (graphType inGS == SoftWired)
                            then T.multiTraverseFullyLabelSoftWiredReduced inGS inData pruneEdges warnPruneEdges leafGraph startVertex delSimple
                            else
                                if (graphType inGS == HardWired)
                                    then T.multiTraverseFullyLabelHardWiredReduced inGS inData leafGraph startVertex delSimple
                                    else error "Unsupported graph type in deleteNetEdge.  Must be soft or hard wired"
                    -}

                    -- check heuristic
                    --logWith LogInfo $ "DNE: " <> -- (show (heuristicDelta, edgeDeleteDelta, heuristicDelta - edgeDeleteDelta)) <> " -> " <>
                    --    (show totalHeuristicCost) <> " vs " <> (show $ (snd5 newPhyloGraph))

                    if force || totalHeuristicCost < (snd5 inPhyloGraph) || isJust inSimAnnealParams then 
                        do
                            newPhyloGraph ←
                                -- check if deletion modified graph
                                if not wasModified then 
                                    pure emptyReducedPhylogeneticGraph

                                else if (graphType inGS == SoftWired)
                                    then T.multiTraverseFullyLabelSoftWiredReduced inGS inData pruneEdges warnPruneEdges leafGraph startVertex delSimple
                                    else
                                        if (graphType inGS == HardWired)
                                            then T.multiTraverseFullyLabelHardWiredReduced inGS inData leafGraph startVertex delSimple
                                            else error "Unsupported graph type in deleteNetEdge.  Must be soft or hard wired"
                                        -- trace ("DNE forced")
                            if (snd5 newPhyloGraph) <= (snd5 inPhyloGraph) then
                                pure newPhyloGraph

                            -- mutation or moving
                            else if force || isJust inSimAnnealParams then 
                                pure newPhyloGraph

                            else 
                                pure emptyReducedPhylogeneticGraph

                    else 
                        pure emptyReducedPhylogeneticGraph


{- | deleteNetworkEdge deletes a network edges from a simple graph
retuns newGraph if can be modified or input graph with Boolean to tell if modified
and contracts, reindexes/names internaledges/veritices around deletion
can't raise to general graph level due to vertex info
in edges (b,a) (c,a) (a,d), deleting (a,b) deletes node a, inserts edge (b,d)
contacts node c since  now in1out1 vertex
checks for chained network edges--can be created by progressive deletion
checks for cycles now
shouldn't need for check for creating a node with children that are both network nodes
since that would require that condition coming in and shodl be there--ie checked earlier in addition and input
-}
deleteNetworkEdge ∷ SimpleGraph → LG.Edge → PhyG (SimpleGraph, Bool)
deleteNetworkEdge inGraph inEdge@(p1, nodeToDelete) =
    if LG.isEmpty inGraph
        then error ("Cannot delete edge from empty graph")
        else
            let childrenNodeToDelete = LG.descendants inGraph nodeToDelete
                parentsNodeToDelete = LG.parents inGraph nodeToDelete
                newGraph = LG.delEdge inEdge inGraph

                -- conversion as if input--see if affects length
                newGraph'' = GO.contractIn1Out1EdgesRename newGraph

            in  -- error conditions and creation of chained network edges (forbidden in phylogenetic graph--causes resolution cache issues)
                if length childrenNodeToDelete /= 1
                    then error ("Cannot delete non-network edge in deleteNetworkEdge: (1)" <> (show inEdge) <> "\n" <> (LG.prettyIndices inGraph))
                    else
                        if length parentsNodeToDelete /= 2
                            then error ("Cannot delete non-network edge in deleteNetworkEdge (2): " <> (show inEdge) <> "\n" <> (LG.prettyIndices inGraph))
                            else -- warning if chained on input, skip if chained net edges in output

                                if (LG.isNetworkNode inGraph p1)
                                    then do
                                        logWith LogWarn ("\tWarning: Chained network nodes in deleteNetworkEdge skipping deletion" <> "\n")
                                        pure (LG.empty, False)
                                    else
                                        if LG.hasChainedNetworkNodes newGraph''
                                            then do
                                                logWith LogWarn ("\tWarning: Chained network nodes in deleteNetworkEdge skipping deletion (2)" <> "\n")
                                                pure (LG.empty, False)
                                            else
                                                if LG.isEmpty newGraph''
                                                    then do
                                                        pure (LG.empty, False)
                                                    else do
                                                        pure (newGraph'', True)



-- | postProcessNetworkDelete postprocesses results from delete actions for "regular" ie non-annealing/Drift network delete operations
postProcessNetworkDelete
    ∷ GlobalSettings
    → ProcessedData
    -> NetParams
    → Int
    → ([ReducedPhylogeneticGraph], VertexCost)
    → Maybe SAParams
    → [ReducedPhylogeneticGraph]
    → [ReducedPhylogeneticGraph]
    → VertexCost
    → VertexCost
    → PhyG ([ReducedPhylogeneticGraph], Int)
postProcessNetworkDelete inGS inData netParams counter (curBestGraphList, _) inSimAnnealParams inPhyloGraphList newGraphList newGraphCost currentCost =
    -- worse graphs found--go on
    if newGraphCost > currentCost
        then do
            deleteAllNetEdges'
                inGS
                inData
                netParams
                (counter + 1)
                ((head inPhyloGraphList) : curBestGraphList, currentCost)
                inSimAnnealParams
                (tail inPhyloGraphList)
        else -- "steepest style descent" abandons existing list if better cost found

            if newGraphCost < currentCost
                then do
                    logWith LogInfo ("\t-> " <> (show newGraphCost))
                    if netSteepest netParams
                        then do
                            deleteAllNetEdges'
                                inGS
                                inData
                                netParams
                                (counter + 1)
                                (newGraphList, newGraphCost)
                                inSimAnnealParams
                                newGraphList
                        else do
                            deleteAllNetEdges'
                                inGS
                                inData
                                netParams
                                (counter + 1)
                                (newGraphList, newGraphCost)
                                inSimAnnealParams
                                (newGraphList <> (tail inPhyloGraphList))
                else -- equal cost
                -- new graph list contains the input graph if equal and filtered unique already in deleteEachNetEdge
                do
                    newCurSameBestList ← GO.selectGraphs Unique (outgroupIndex inGS) (netKeepNum netParams) 0.0 $ curBestGraphList <> newGraphList
                    deleteAllNetEdges'
                        inGS
                        inData
                        netParams
                        (counter + 1)
                        (newCurSameBestList, currentCost)
                        inSimAnnealParams
                        (tail inPhyloGraphList)


{-----------------------------------
    Insert network edge functions
------------------------------------}

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
                    -- this recurses until max net edges or no new edges reached
                    if isNothing inSimAnnealParams then do

                        -- this examines all possible net adds via heursitc first then verifies
                        -- could get cost here in second tuple place but need to get best anyway later
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


                        --logWith LogInfo ("IANE: " <> (show $ fmap snd5 newGraphList) <> " " <> (show newGraphCost))
                        -- this will recurse back if required
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

                    -- SA/Drift stuff--single call recurses within insertRoundsSA
                    else 
                        insertRoundsSA inGS inData (netParams {netEditType = NetAdd}) counter (curBestGraphList, curBestGraphCost) inSimAnnealParams inPhyloGraphList

{- | insertRoundsSA adds network edges one each each round using Simulated Annealing/Drift

    Randomizes graphs to examine in insertEachNetEdgeHeuristicGather so can do multiple rounds there
    After SA insert can return mutated or do a delete to get back to "optimal" edges
-}
insertRoundsSA
    ∷ GlobalSettings
    → ProcessedData
    → NetParams
    → Int
    → ([ReducedPhylogeneticGraph], VertexCost)
    → Maybe SAParams
    → [ReducedPhylogeneticGraph]
    → PhyG ([ReducedPhylogeneticGraph], Int)
insertRoundsSA inGS inData netParams counter (curBestGraphList, curBestGraphCost) inSimAnnealParams inPhyloGraphList = 
    let -- create list of params with unique list of random values for rounds of annealing
        annealingRounds = rounds $ fromJust inSimAnnealParams
        annealParamList = replicate annealingRounds inSimAnnealParams

        -- set up parallel
        addAction :: (Maybe SAParams, ReducedPhylogeneticGraph)→ PhyG ([ReducedPhylogeneticGraph], Int)
        addAction = insertAllNetEdgesSA inGS inData netParams counter curBestGraphCost 

        deleteAction :: ReducedPhylogeneticGraph -> PhyG ([ReducedPhylogeneticGraph], VertexCost, Maybe SAParams)
        deleteAction = deleteEachNetEdge inGS inData netParams False Nothing
        
    in do
        graphsToAnneal <- if length inPhyloGraphList == 1 then 
                                pure $ replicate annealingRounds $ head  inPhyloGraphList
                          else do
                                shuffledLists <- mapM shuffleList $ replicate annealingRounds inPhyloGraphList
                                pure $ fmap head shuffledLists

        addActionPar <- getParallelChunkTraverse 
        newGraphTripleList <- addActionPar addAction (zip annealParamList graphsToAnneal)

        uniqueList ← GO.selectGraphs Unique (outgroupIndex inGS) (netKeepNum netParams) 0.0 (concat $ fmap fst newGraphTripleList)

        let uniqueCost = if (not .null) uniqueList then
                                minimum $ fmap snd5 uniqueList
                         else infinity

        -- return mutated for GA
        if netReturnMutated netParams then do
           pure (uniqueList, counter + annealingRounds)

        
        else do
            -- netADD SA
            if netEditType netParams /= NetMove then do
                logWith LogInfo "\t\tDeleting extraneous edges after network addition Simulated Annealing/Drifing\n"
                -- delete edges to get back to "best" edge lists
                
                (deletedGraphList, _) <- deleteAllNetEdges' inGS inData (netParams {netSteepest = False}) counter (uniqueList, uniqueCost) Nothing uniqueList 

                -- return better and equal including inputs
                finalList ← GO.selectGraphs Best (outgroupIndex inGS) (netKeepNum netParams) 0.0 $ inPhyloGraphList <> deletedGraphList -- (concat $ fmap fst3 deletedTripleList)

                let netNodesList = fmap length $ fmap (fth4 . LG.splitVertexList . thd5) finalList

                {-
                let finalCost = if (not .null) finalList then
                                    minimum $ fmap snd5 finalList
                                else infinity     
                if (minimum $ fmap snd5 finalList) < (minimum $ fmap snd5 uniqueList) then
                    logWith LogInfo ("\t\t-> " <> (show $ min (minimum $ fmap snd5 finalList) (minimum $ fmap snd5 uniqueList)) <> " with " <> (show netNodesList) <> " netWork nodes\n")
                else 
                    logWith LogInfo ""
                -}
                pure (finalList, counter + annealingRounds)

            -- netMOve so regular Move back asfter SA
            else do
                logWith LogInfo "\t\tMoving edges after network addition Simulated Annealing/Drifing\n"
                -- delete edges to get back to "best" edge lists
                
                (movedGraphList, _) <- moveAllNetEdges' inGS inData (netParams {netSteepest = False}) counter (uniqueList, uniqueCost) Nothing uniqueList 

                -- return better and equal including inputs
                finalList ← GO.selectGraphs Best (outgroupIndex inGS) (netKeepNum netParams) 0.0 $ inPhyloGraphList <> movedGraphList -- (concat $ fmap fst3 deletedTripleList)

                {-
                let netNodesList = fmap length $ fmap (fth4 . LG.splitVertexList . thd5) finalList

                
                let finalCost = if (not .null) finalList then
                                    minimum $ fmap snd5 finalList
                                else infinity     
                if (minimum $ fmap snd5 finalList) < (minimum $ fmap snd5 uniqueList) then
                    logWith LogInfo ("\t\t-> " <> (show $ min (minimum $ fmap snd5 finalList) (minimum $ fmap snd5 uniqueList)) <> " with " <> (show netNodesList) <> " netWork nodes\n")
                else 
                    logWith LogInfo ""
                -}
                pure (finalList, counter + 1)


{- insertaAllNetEdgesSA performes a single SA/Drift Add trajectory 
    Will continue adding edges if selected by SA/Drift until
    edges hit maximum number, or SA/Drift 
    steps end
-}
insertAllNetEdgesSA
    ∷ GlobalSettings
    → ProcessedData
    → NetParams
    → Int
    → VertexCost
    → (Maybe SAParams, ReducedPhylogeneticGraph)
    → PhyG ([ReducedPhylogeneticGraph], Int)
insertAllNetEdgesSA inGS inData netParams counter curBestCost (inSimAnnealParams, inPhyloGraph) =
    if U.isSimAnnealTerminated inSimAnnealParams then
        pure ([inPhyloGraph], counter)  
    else 
        let netNodes = fth4 $ LG.splitVertexList $ thd5 inPhyloGraph
        in
        if length netNodes >= (netMaxEdges netParams) then do
                logWith LogInfo $ unwords ["Maximum number of network edges reached:", show $ length netNodes, "\n"]
                pure ([inPhyloGraph], counter)  
        else 
            do
                -- generate graphs with randomization on heuristic costs
                (newGraphList, newCost, newSAParams) ←
                    insertEachNetEdgeHeuristicGather 
                        inGS
                        inData
                        netParams
                        Nothing
                        inSimAnnealParams
                        inPhyloGraph

                --logWith LogInfo ("SA/Drift graph return: " <> (show $ fmap snd5 newGraphList) <> " " <> (show newCost))
                bestNewGraphList ← GO.selectGraphs Best (outgroupIndex inGS) (netKeepNum netParams) 0 newGraphList

                -- terminate if not getting any new graphs (could happen if lots of edges already)
                if null bestNewGraphList then
                        pure ([inPhyloGraph], counter + 1)

                -- check if accept based on SA/Drift
                 else do
                        let bestNewGraph = head bestNewGraphList
                        (acceptGraph, newSAParams) <- U.simAnnealAccept inSimAnnealParams curBestCost newCost

                        let newCurBestCost = min curBestCost newCost

                        if acceptGraph then 
                            insertAllNetEdgesSA inGS inData netParams (counter + 1) newCurBestCost (newSAParams, bestNewGraph) 

                        else 
                            insertAllNetEdgesSA inGS inData netParams (counter + 1) newCurBestCost (newSAParams, inPhyloGraph)

                                    

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

Operates by first performing heuristic cost estimate for all possibilities,
then fully evaluates the best (ie. lowest delta number) graphs based on
HeursticCheck option.

For SA/Drift just returns all best afte verify whether better or not 
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
                            candidateNetworkEdgeList ← getPermissibleEdgePairs inGS (thd5 inPhyloGraph)

                            -- only take some smaller sample for SA/Drift
                            candidateNetworkEdgeList' <- if isNothing inSimAnnealParams then 
                                                            pure candidateNetworkEdgeList
                                                         else do
                                                            shuffledList <- shuffleList candidateNetworkEdgeList
                                                            pure $ take ((graphsSteepest inGS) * (graphsSteepest inGS)) shuffledList
                            logWith LogInfo ("\t\tExamining at most " <> (show $ length candidateNetworkEdgeList') <> " candidate edge pairs" <> "\n")


                            -- get heuristic costs and simple graphs
                            heurCostSimpleGraphPairList <- 
                                getParallelChunkTraverse >>= \pTraverse →
                                    heuristicAction `pTraverse` candidateNetworkEdgeList'

                            -- filter out non-phylo graphs and sort on heuristic cost
                            let nonInfinitePairList = filter ((/= infinity) .fst) heurCostSimpleGraphPairList
                            let candidatePairList = L.sortOn fst nonInfinitePairList

                            let graphsToBeEvaluated = 
                                        -- make sure take some number of SA/Drif
                                        if isJust inSimAnnealParams then
                                             fmap snd $ take (graphsSteepest inGS) candidatePairList

                                        else if (netCheckHeuristic netParams) == BestOnly then
                                            fmap snd $ take 1 candidatePairList

                                        else if (netCheckHeuristic netParams) == Better then
                                            fmap snd $ filter ((< 0) . fst) candidatePairList

                                        else if (netCheckHeuristic netParams) == BetterN then 
                                            fmap snd $ take (graphsSteepest inGS) candidatePairList

                                        else --BestAll
                                            fmap snd candidatePairList

                            
                            --logWith LogInfo ("GTBE: " <> (show $ fmap fst candidatePairList))

                            if null graphsToBeEvaluated then
                                    pure ([inPhyloGraph], currentCost, inSimAnnealParams)

                            else do
                                    -- rediagnose some fraction of returned simple graphs--lazy in cost so return only thos need nlater
                                    diagnoseActionPar <- (getParallelChunkTraverseBy snd5)
                                    checkedGraphCosts <- diagnoseActionPar diagnoseAction graphsToBeEvaluated

                                    (newGraphList, newSAParams) <-
                                                let genNewSimAnnealParams =
                                                        if isNothing inSimAnnealParams
                                                            then Nothing
                                                        else U.incrementSimAnnealParams inSimAnnealParams
                                                in
                                                pure (checkedGraphCosts, genNewSimAnnealParams)
                                           
                                    let minCost =
                                            if null graphsToBeEvaluated || null newGraphList
                                                then infinity
                                                else minimum $ fmap snd5 newGraphList
                                    
                                    
                                    if minCost < (snd5 inPhyloGraph) then 
                                        logWith LogInfo ("\t\t-> " <> (show minCost) <> "\n")
                                    else if isJust inSimAnnealParams 
                                        then logWith LogInfo ("\t\t-> " <> (show minCost) <> "\n")
                                    else 
                                        logWith LogInfo ("")

                                    let genNewSimAnnealParams =
                                            if isNothing inSimAnnealParams then 
                                                Nothing
                                            else 
                                                U.incrementSimAnnealParams inSimAnnealParams
                                        
                                    
                                    -- always return if better in any conditions
                                    if minCost <= (snd5 inPhyloGraph) then
                                        --logWith LogInfo ("IENEHG: " <> (show (minCost, (snd5 inPhyloGraph) )))
                                        do
                                            newBestGraphList <- GO.selectGraphs Best (outgroupIndex inGS) (netKeepNum netParams) 0 newGraphList
                                            pure (newBestGraphList, minCost, genNewSimAnnealParams)
                                    
                                    -- return if worse for SA/Drift
                                    else do
                                        if isNothing inSimAnnealParams then 
                                            pure ([], minCost, inSimAnnealParams)
                                        else do
                                            newRandGraphList <- GO.selectGraphs AtRandom (outgroupIndex inGS) (netKeepNum netParams) 0 newGraphList
                                            pure (newRandGraphList, minCost, genNewSimAnnealParams)
                                
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
        {-
        removeNodeLabels :: forall {f1 :: * -> *} {f2 :: * -> *}
                                   {f3 :: * -> *} {f4 :: * -> *} {a1} {a2} {a3} {a4} {a5} {a6}.
                            (Functor f1, Functor f2, Functor f3, Functor f4) =>
                            (LG.LNode a1, LG.LNode a2, f1 (LG.LNode a3), f2 (LG.LNode a4),
                             f3 (LG.LNode a5), f4 (LG.LNode a6))
                            -> (LG.Node, LG.Node, f1 LG.Node, f2 LG.Node, f3 LG.Node,
                                f4 LG.Node)
        -}
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

{- These heuristics do not seem to work very well--high variance with reoptimized costs-}

{- Reverse of heuristicAdd Delta

Start with edges (u',b), (u', y), (y, v'), (u,a), (u, x), (x, y), (x,v)
(x,y) is the edge to delete resulting in subtrees U' -> b and v'; u -> v and a
existing cost for graph fragment would be preliminary cost (full subtree) 
would be cost u'  + cost u.  Cost of new fragment is cost (full subtree)
of new u and u' wiich would be cost of (v" + b) and (v + a)

This estimate should be better (less variance) than add.
Does include the benefit of deleting a network edge
-}

heuristicDeleteDelta ∷ GlobalSettings → ReducedPhylogeneticGraph → LG.Edge → PhyG VertexCost
heuristicDeleteDelta inGS inPhyloGraph (x,y) =
    if LG.isEmpty (fst5 inPhyloGraph)
        then error "Empty graph in heuristicDeleteDelta"
        else
            let -- base vertices of graph fragment
                u' = head $ LG.parents (fst5 inPhyloGraph) y
                u  = head $ LG.parents (fst5 inPhyloGraph) x
                
                -- single child
                v' = head $ LG.descendants (fst5 inPhyloGraph) y
                v  = head $ filter (/= y) $ LG.descendants (fst5 inPhyloGraph) x
                a  = head $ filter (/= x) $ LG.descendants (fst5 inPhyloGraph) u
                b  = head $ filter (/= y) $ LG.descendants (fst5 inPhyloGraph) u'

                -- use minimia of resolution data for each block/character

                --uResolutionData = vertexResolutionData $ fromJust $ LG.lab (thd5 inPhyloGraph) u
                --u'ResolutionData = vertexResolutionData $ fromJust $ LG.lab (thd5 inPhyloGraph) u'

                aVertexInfo = fromJust $ LG.lab (thd5 inPhyloGraph) a
                bVertexInfo = fromJust $ LG.lab (thd5 inPhyloGraph) b

                --aResolutionData = vertexResolutionData aVertexInfo 
                --bResolutionData = vertexResolutionData bVertexInfo

                vVertexInfo  = fromJust $ LG.lab (thd5 inPhyloGraph) v
                v'VertexInfo = fromJust $ LG.lab (thd5 inPhyloGraph) v'
                
                --vResolutionData  = vertexResolutionData  vVertexInfo 
                --v'ResolutionData = vertexResolutionData  v'VertexInfo 


                -- get best resolution cost dor each character and block

                costSubGraphU  = GO.getSoftWiredNodeSubGraphCost inPhyloGraph u
                costSubGraphU' = GO.getSoftWiredNodeSubGraphCost inPhyloGraph u'

            in do
                newDataU  <- NEW.getOutDegree2VertexSoftWired inGS (fft5 inPhyloGraph) (-1) (v, vVertexInfo) (a, aVertexInfo) (thd5 inPhyloGraph)
                newDataU' <- NEW.getOutDegree2VertexSoftWired inGS (fft5 inPhyloGraph) (-1) (v', v'VertexInfo) (b, bVertexInfo) (thd5 inPhyloGraph)

                let costSubGraphNewU  = V.sum $ fmap GO.minBlockResolutionCost (vertexResolutionData newDataU)
                let costSubGraphNewU' = V.sum $ fmap GO.minBlockResolutionCost (vertexResolutionData newDataU')
                
                logWith LogInfo $ "HDD: " <> (show (costSubGraphU,costSubGraphU')) <> " to " <> (show (costSubGraphNewU, costSubGraphNewU')) <> " net " <> (show $ (costSubGraphNewU +  costSubGraphNewU') - (costSubGraphU + costSubGraphU')) <> "\n"

                -- new should be higher or equal cost than existing
                if (costSubGraphNewU +  costSubGraphNewU') - (costSubGraphU + costSubGraphU') < 0.0 then 
                    pure 0.0
                else 
                    pure $ (costSubGraphNewU +  costSubGraphNewU') - (costSubGraphU + costSubGraphU')


{- | heuristic add delta based on new display tree and delta from existing costs by block--assumming < 0
original edges subtree1 ((u,l),(u,v)) and subtree2 ((u',v'),(u',l')) create a directed edge from
subtree 1 to subtree 2 via (x -> y)
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
                action = getBlockDeltaAdd (u, v, u', v', a, b)
            in do
                actionPar <- getParallelChunkTraverse 
                -- blockDeltaV <- V.zipWith (getBlockDeltaAdd (u, v, u', v', a, b)) (zip blockTrees (fft5 inPhyloGraph))
                blockDeltaL <- actionPar action (V.toList $ V.zip blockTrees (fft5 inPhyloGraph))
                pure $ sum blockDeltaL


{- | getBlockDeltaAdd determines the network add delta for each block (vector of characters)
if existing is lower then zero, else (existing - new)
-}
getBlockDeltaAdd
    ∷ (LG.Node, LG.Node, LG.Node, LG.Node, LG.Node, LG.Node) → (V.Vector DecoratedGraph, V.Vector CharInfo) → PhyG VertexCost
getBlockDeltaAdd (u, v, u', v', a, b) (inCharV, charInfoV) =
    if V.null inCharV
        then error "Empty charcter tree vector in getBlockDeltaAdd"
        else
            let action :: (DecoratedGraph, CharInfo) → (VertexCost, VertexCost)
                action = getCharacterAddDelta (u, v, u', v', a, b)
            in do
                 actionPar ← getParallelChunkMap
                 let result = actionPar action (V.toList $ V.zip inCharV charInfoV)
                 let (charNewV, charExistingV) = unzip result
                 let newCost = sum charNewV
                 let existingCost = sum charExistingV
                 if (newCost < existingCost)
                    then pure $ newCost - existingCost
                    else pure $ 0.0

{- | getCharacterAddDelta determines the network add delta for each block (vector of characters)
if existing is lower then zero, else (existing - new)
 calculate d(u,v) + d(u',v') [existing display tree cost estimate] compared to
 d((union u,v), v') - d(u'.v')
need to use final assignemnts--so set prelim to final first
Since a distance--no need for No change cost adjustment
-}
getCharacterAddDelta
    ∷ (LG.Node, LG.Node, LG.Node, LG.Node, LG.Node, LG.Node) → (DecoratedGraph, CharInfo) → (VertexCost, VertexCost)
getCharacterAddDelta (_, v, _, v', a, b) (inCharTree, charInfo) =
    -- getCharacterAddDelta (u,v,u',v',a,b) inCharTree charInfo =
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


{-
-- | insertNetEdgeBothDirections calls insertNetEdge for both u -> v and v -> u new edge orientations
insertNetEdgeBothDirections :: GlobalSettings -> ProcessedData -> ReducedPhylogeneticGraph ->  Maybe VertexCost -> (LG.LEdge b, LG.LEdge b) -> [ReducedPhylogeneticGraph]
insertNetEdgeBothDirections inGS inData inPhyloGraph preDeleteCost (u,v) = fmap (insertNetEdge inGS inData inPhyloGraph preDeleteCost) [(u,v), (v,u)]
-}



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


