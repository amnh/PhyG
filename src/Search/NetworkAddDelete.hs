{- |
Module specifying graph egde adding and deleting functions
-}
module Search.NetworkAddDelete (
    deleteAllNetEdges,
    insertAllNetEdges,
    moveAllNetEdges,
    deltaPenaltyAdjustment,
    deleteNetEdge,
    deleteOneNetAddAll,
    addDeleteNetEdges,
    getCharacterDelta,
    getBlockDelta,
    -- these are not used but to quiet warnings
    heuristicDeleteDelta,
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


{- |
'addDeleteNetEdges' is a wrapper for 'addDeleteNetEdges'' allowing for multiple simulated annealing rounds.
-}
addDeleteNetEdges
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Int
    → Int
    → Int
    → Bool
    → Bool
    → Bool
    → ([ReducedPhylogeneticGraph], VertexCost)
    → (Maybe SAParams, [ReducedPhylogeneticGraph])
    → PhyG ([ReducedPhylogeneticGraph], Int)
addDeleteNetEdges inGS inData maxNetEdges numToKeep maxRounds counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) = case inSimAnnealParams of
    Nothing →
        addDeleteNetEdges'
            inGS
            inData
            maxNetEdges
            numToKeep
            maxRounds
            counter
            returnMutated
            doSteepest
            doRandomOrder
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
                    maxNetEdges
                    numToKeep
                    maxRounds
                    counter
                    returnMutated
                    doSteepest
                    doRandomOrder
                    (curBestGraphList, curBestGraphCost)
        in  do
                addDeleteResult ←
                    getParallelChunkTraverse >>= \pTraverse →
                        pTraverse action . zip saParamList $ replicate annealingRounds inPhyloGraphList
                let (annealRoundsList, counterList) = unzip addDeleteResult
                GO.selectGraphs Best numToKeep 0 (fold annealRoundsList) <&> \x → (x, sum counterList)


-- | addDeleteNetEdges'' is wrapper around addDeleteNetEdges' to use parmap
addDeleteNetEdges''
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Int
    → Int
    → Int
    → Bool
    → Bool
    → Bool
    → ([ReducedPhylogeneticGraph], VertexCost)
    → (Maybe SAParams, [ReducedPhylogeneticGraph])
    → PhyG ([ReducedPhylogeneticGraph], Int)
addDeleteNetEdges'' inGS inData maxNetEdges numToKeep maxRounds counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
    addDeleteNetEdges'
        inGS
        inData
        maxNetEdges
        numToKeep
        maxRounds
        counter
        returnMutated
        doSteepest
        doRandomOrder
        (curBestGraphList, curBestGraphCost)
        inSimAnnealParams
        inPhyloGraphList


{- | addDeleteNetEdges' removes each edge and adds an edge to all possible places (or steepest) each round
until no better or additional graphs are found (or max rounds met)
call with ([], infinity) [single input graph]
doesn't have to be random, but likely to converge quickly if not
-}
addDeleteNetEdges'
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Int
    → Int
    → Int
    → Bool
    → Bool
    → Bool
    → ([ReducedPhylogeneticGraph], VertexCost)
    → Maybe SAParams
    → [ReducedPhylogeneticGraph]
    → PhyG ([ReducedPhylogeneticGraph], Int)
addDeleteNetEdges' inGS inData maxNetEdges numToKeep maxRounds counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) inSimAnnealParams = \case
    [] → pure (take numToKeep curBestGraphList, counter)
    -- if hit maxmimum rounds then return
    inPhyloGraphList → case counter `compare` maxRounds of
        EQ → pure (take numToKeep curBestGraphList, counter)
        -- other wise add/delete
        _ → do
            -- insert edges first
            (insertGraphList, _) ←
                insertAllNetEdges'
                    inGS
                    inData
                    maxNetEdges
                    numToKeep
                    counter
                    returnMutated
                    doSteepest
                    doRandomOrder
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
                    newList ← GO.selectGraphs Best (maxBound ∷ Int) 0 gs
                    pure (newList, snd5 $ head newList, newList)

            -- delete edges
            (deleteGraphList, _) ←
                deleteAllNetEdges'
                    inGS
                    inData
                    maxNetEdges
                    numToKeep
                    counter
                    returnMutated
                    doSteepest
                    doRandomOrder
                    (insertGraphList', insertGraphCost)
                    (head updatedSAParamList)
                    toDeleteList

            -- gather beter if any
            (newBestGraphList, newBestGraphCost, graphsToDoNext) ← case deleteGraphList of
                [] → pure (curBestGraphList, curBestGraphCost, inPhyloGraphList)
                gs → do
                    newDeleteGraphs ← GO.selectGraphs Best (maxBound ∷ Int) 0 gs
                    pure (newDeleteGraphs, snd5 $ head newDeleteGraphs, newDeleteGraphs)

            -- check is same then return
            case newBestGraphCost `compare` curBestGraphCost of
                EQ → pure (take numToKeep curBestGraphList, counter)
                -- if better (or nothing) keep going
                _ →
                    addDeleteNetEdges'
                        inGS
                        inData
                        maxNetEdges
                        numToKeep
                        maxRounds
                        (counter + 1)
                        returnMutated
                        doSteepest
                        doRandomOrder
                        (newBestGraphList, newBestGraphCost)
                        (last updatedSAParamList)
                        graphsToDoNext


-- | moveAllNetEdges is a wrapper for moveAllNetEdges' allowing for multiple simulated annealing rounds
moveAllNetEdges
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Int
    → Int
    → Bool
    → Bool
    → Bool
    → ([ReducedPhylogeneticGraph], VertexCost)
    → (Maybe SAParams, [ReducedPhylogeneticGraph])
    → PhyG ([ReducedPhylogeneticGraph], Int)
moveAllNetEdges inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) = case inSimAnnealParams of
    Nothing →
        moveAllNetEdges'
            inGS
            inData
            maxNetEdges
            numToKeep
            counter
            returnMutated
            doSteepest
            doRandomOrder
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
                    maxNetEdges
                    numToKeep
                    counter
                    returnMutated
                    doSteepest
                    doRandomOrder
                    (curBestGraphList, curBestGraphCost)
        in  do
                moveResult ←
                    getParallelChunkTraverse >>= \pTraverse →
                        pTraverse action . zip saParamList $ replicate annealingRounds inPhyloGraphList
                let (annealRoundsList, counterList) = unzip moveResult
                GO.selectGraphs Best numToKeep 0 (fold annealRoundsList) <&> \x → (x, sum counterList)


-- | moveAllNetEdges'' is wrapper around moveAllNetEdges' to use parmap
moveAllNetEdges''
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Int
    → Int
    → Bool
    → Bool
    → Bool
    → ([ReducedPhylogeneticGraph], VertexCost)
    → (Maybe SAParams, [ReducedPhylogeneticGraph])
    → PhyG ([ReducedPhylogeneticGraph], Int)
moveAllNetEdges'' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
    moveAllNetEdges'
        inGS
        inData
        maxNetEdges
        numToKeep
        counter
        returnMutated
        doSteepest
        doRandomOrder
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
    → Int
    → Int
    → Int
    → Bool
    → Bool
    → Bool
    → ([ReducedPhylogeneticGraph], VertexCost)
    → Maybe SAParams
    → [ReducedPhylogeneticGraph]
    → PhyG ([ReducedPhylogeneticGraph], Int)
moveAllNetEdges' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) inSimAnnealParams = \case
    [] → pure (take numToKeep curBestGraphList, counter)
    firstPhyloGraph : otherPhyloGraphs
        | LG.isEmpty $ fst5 firstPhyloGraph →
            moveAllNetEdges'
                inGS
                inData
                maxNetEdges
                numToKeep
                counter
                returnMutated
                doSteepest
                doRandomOrder
                (curBestGraphList, curBestGraphCost)
                inSimAnnealParams
                otherPhyloGraphs
    firstPhyloGraph : otherPhyloGraphs →
        let currentCost = min curBestGraphCost $ snd5 firstPhyloGraph
            -- parallel setup
            action ∷ LG.Edge → PhyG [ReducedPhylogeneticGraph]
            action = deleteOneNetAddAll' inGS inData maxNetEdges numToKeep doSteepest doRandomOrder firstPhyloGraph inSimAnnealParams
        in  do
                -- randomize order of edges to try moving
                netEdgeList ←
                    let edges = LG.labNetEdges $ thd5 firstPhyloGraph
                        permutationOf
                            | doRandomOrder = shuffleList
                            | otherwise = pure
                     in permutationOf edges

                deleteResult ←
                    getParallelChunkTraverse >>= \pTraverse →
                        pTraverse (action . LG.toEdge) netEdgeList

                let newGraphList' = fold deleteResult
                newGraphList ← GO.selectGraphs Best numToKeep 0 newGraphList'
                let newGraphCost = case newGraphList' of
                        [] → infinity
                        (_, c, _, _, _) : _ → c

                -- if graph is a tree no edges to delete
                case netEdgeList of
                    [] → (firstPhyloGraph : otherPhyloGraphs, counter) <$ logWith LogInfo "\t\tGraph in move has no network edges to move--skipping\n"
                    e : es → case inSimAnnealParams of
                        Nothing → case newGraphCost `compare` currentCost of
                            GT →
                                moveAllNetEdges'
                                    inGS
                                    inData
                                    maxNetEdges
                                    numToKeep
                                    (counter + 1)
                                    returnMutated
                                    doSteepest
                                    doRandomOrder
                                    (firstPhyloGraph : curBestGraphList, currentCost)
                                    inSimAnnealParams
                                    otherPhyloGraphs
                            LT
                                | doSteepest →
                                    moveAllNetEdges'
                                        inGS
                                        inData
                                        maxNetEdges
                                        numToKeep
                                        (counter + 1)
                                        returnMutated
                                        doSteepest
                                        doRandomOrder
                                        (newGraphList, newGraphCost)
                                        inSimAnnealParams
                                        newGraphList
                            LT →
                                moveAllNetEdges'
                                    inGS
                                    inData
                                    maxNetEdges
                                    numToKeep
                                    (counter + 1)
                                    returnMutated
                                    doSteepest
                                    doRandomOrder
                                    (newGraphList, newGraphCost)
                                    inSimAnnealParams
                                    (newGraphList <> otherPhyloGraphs)
                            EQ → do
                                -- new graph list contains the input graph if equal and filterd unique already in moveAllNetEdges
                                newCurSameBestList ← GO.selectGraphs Unique numToKeep 0 $ curBestGraphList <> newGraphList
                                moveAllNetEdges'
                                    inGS
                                    inData
                                    maxNetEdges
                                    numToKeep
                                    (counter + 1)
                                    returnMutated
                                    doSteepest
                                    doRandomOrder
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
                                    uniqueGraphList ← GO.selectGraphs Unique numToKeep 0 newGraphList'

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
                                                maxNetEdges
                                                numToKeep
                                                (counter + 1)
                                                returnMutated
                                                doSteepest
                                                doRandomOrder
                                                ((head uniqueGraphList) : curBestGraphList, annealBestCost)
                                                newSAParams
                                                (nextUniqueList <> otherPhyloGraphs)
                                        LT → do
                                            moveAllNetEdges'
                                                inGS
                                                inData
                                                maxNetEdges
                                                numToKeep
                                                (counter + 1)
                                                returnMutated
                                                doSteepest
                                                doRandomOrder
                                                (curBestGraphList, annealBestCost)
                                                newSAParams
                                                (nextUniqueList <> otherPhyloGraphs)

                                        -- if want non-optimized list for GA or whatever
                                        _ | returnMutated → pure (take numToKeep curBestGraphList, counter)
                                        -- optimize list and return
                                        _ → do
                                            (bestMoveList', counter') ←
                                                moveAllNetEdges'
                                                    inGS
                                                    inData
                                                    maxNetEdges
                                                    numToKeep
                                                    (counter + 1)
                                                    False
                                                    doSteepest
                                                    doRandomOrder
                                                    ([], annealBestCost)
                                                    Nothing
                                                    $ take numToKeep curBestGraphList
                                            bestMoveList ← GO.selectGraphs Best numToKeep 0 bestMoveList'
                                            pure (take numToKeep bestMoveList, counter')


-- | (curBestGraphList, annealBestCost) is a wrapper for moveAllNetEdges' allowing for multiple simulated annealing rounds
insertAllNetEdges
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Int
    → Int
    → Int
    → Bool
    → Bool
    → Bool
    → ([ReducedPhylogeneticGraph], VertexCost)
    → (Maybe SAParams, [ReducedPhylogeneticGraph])
    → PhyG ([ReducedPhylogeneticGraph], Int)
insertAllNetEdges inGS inData maxNetEdges numToKeep maxRounds counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
    let -- parallel setup
        randAction ∷ [Int] → PhyG ([ReducedPhylogeneticGraph], Int)
        randAction =
            insertAllNetEdgesRand
                inGS
                inData
                maxNetEdges
                numToKeep
                counter
                returnMutated
                doSteepest
                (curBestGraphList, curBestGraphCost)
                Nothing
                inPhyloGraphList

        action ∷ (Maybe SAParams, [ReducedPhylogeneticGraph]) → PhyG ([ReducedPhylogeneticGraph], Int)
        action =
            insertAllNetEdges''
                inGS
                inData
                maxNetEdges
                numToKeep
                counter
                returnMutated
                doSteepest
                doRandomOrder
                (curBestGraphList, curBestGraphCost)
    in  if isNothing inSimAnnealParams
            then -- check for multiple rounds of addition--if > 1 then need to randomize order

                if maxRounds == 1
                    then
                        insertAllNetEdges'
                            inGS
                            inData
                            maxNetEdges
                            numToKeep
                            counter
                            returnMutated
                            doSteepest
                            doRandomOrder
                            (curBestGraphList, curBestGraphCost)
                            Nothing
                            inPhyloGraphList
                    else do
                        -- this odd contruction is to ensure that the correct number of annealing rounds are
                        -- done even though the values in the lists are ignored.  Should be refactored.
                        -- this happened during migration to monadic getRandom
                        let intList = replicate maxRounds (0 ∷ Int)
                        let intListList = replicate maxRounds intList
                        insertGraphResult ←
                            getParallelChunkTraverse >>= \pTraverse →
                                randAction `pTraverse` intListList

                        let (insertGraphList, counterList) = unzip insertGraphResult
                        -- insert functions take care of returning "better" or empty
                        -- should be empty if nothing better
                        GO.selectGraphs Best numToKeep 0 (fold insertGraphList) <&> \x → (x, sum counterList)
            else
                let -- create list of params with unique list of random values for rounds of annealing
                    annealingRounds = rounds $ fromJust inSimAnnealParams
                    annealParamGraphList = replicate annealingRounds inSimAnnealParams
                in  do
                        insertGraphResult ←
                            getParallelChunkTraverse >>= \pTraverse →
                                pTraverse action . zip annealParamGraphList $ replicate annealingRounds inPhyloGraphList
                        let (annealRoundsList, counterList) = unzip insertGraphResult
                        let selectionType
                                | not returnMutated || isNothing inSimAnnealParams = Best
                                | otherwise = Unique

                        GO.selectGraphs selectionType numToKeep 0 (fold annealRoundsList) <&> \x → (x, sum counterList)


-- | insertAllNetEdgesRand is a wrapper around insertAllNetEdges'' to pass unique randomLists to insertAllNetEdges'
insertAllNetEdgesRand
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Int
    → Int
    → Bool
    → Bool
    → ([ReducedPhylogeneticGraph], VertexCost)
    → Maybe SAParams
    → [ReducedPhylogeneticGraph]
    → [Int]
    → PhyG ([ReducedPhylogeneticGraph], Int)
insertAllNetEdgesRand inGS inData maxNetEdges numToKeep counter returnMutated doSteepest (curBestGraphList, curBestGraphCost) inSimAnnealParams inPhyloGraphList _ =
    insertAllNetEdges'
        inGS
        inData
        maxNetEdges
        numToKeep
        counter
        returnMutated
        doSteepest
        True
        (curBestGraphList, curBestGraphCost)
        inSimAnnealParams
        inPhyloGraphList


-- | insertAllNetEdges'' is a wrapper around insertAllNetEdges' to allow for seqParMap
insertAllNetEdges''
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Int
    → Int
    → Bool
    → Bool
    → Bool
    → ([ReducedPhylogeneticGraph], VertexCost)
    → (Maybe SAParams, [ReducedPhylogeneticGraph])
    → PhyG ([ReducedPhylogeneticGraph], Int)
insertAllNetEdges'' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
    insertAllNetEdges'
        inGS
        inData
        maxNetEdges
        numToKeep
        counter
        returnMutated
        doSteepest
        doRandomOrder
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
    → Int
    → Int
    → Int
    → Bool
    → Bool
    → Bool
    → ([ReducedPhylogeneticGraph], VertexCost)
    → Maybe SAParams
    → [ReducedPhylogeneticGraph]
    → PhyG ([ReducedPhylogeneticGraph], Int)
insertAllNetEdges' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) inSimAnnealParams = \case
    -- this logic so don't return mutated if finish insertion before hitting other stopping points
    -- and don't want mutated--straight deletion on current best graphs
    [] → case inSimAnnealParams of
        Just _
            | not returnMutated →
                deleteAllNetEdges'
                    inGS
                    inData
                    maxNetEdges
                    numToKeep
                    counter
                    False
                    doSteepest
                    doRandomOrder
                    ([], curBestGraphCost)
                    Nothing
                    (take numToKeep curBestGraphList)
        _ → pure (take numToKeep curBestGraphList, counter)
    firstPhyloGraph : otherPhyloGraphs →
        let currentCost = min curBestGraphCost $ snd5 firstPhyloGraph

            -- check for max net edges
            (_, _, _, netNodes) = LG.splitVertexList $ thd5 firstPhyloGraph
        in  do
                logWith LogInfo ("\t\tNumber of network edges: " <> (show $ length netNodes) <> "\n")

                (newGraphList, _, newSAParams) ←
                    insertEachNetEdge
                        inGS
                        inData
                        maxNetEdges
                        numToKeep
                        doSteepest
                        doRandomOrder
                        Nothing
                        inSimAnnealParams
                        firstPhyloGraph

                bestNewGraphList ← GO.selectGraphs Best numToKeep 0 newGraphList
                let newGraphCost = case bestNewGraphList of
                        [] → infinity
                        (_, c, _, _, _) : _ → c

                -- logWith LogInfo ("\t\tNumber of network edges: " <> (show $ length netNodes) <> "\n")

                case length netNodes `compare` maxNetEdges of
                    LT → case newGraphList of
                        [] → do
                            logWith LogInfo $ unwords ["\t\tNumber of network edges:", show $ length netNodes, "\n"]
                            pure (take numToKeep curBestGraphList, counter)
                        -- regular insert keeping best
                        g : gs → case inSimAnnealParams of
                            Nothing →
                                postProcessNetworkAdd
                                    inGS
                                    inData
                                    maxNetEdges
                                    numToKeep
                                    counter
                                    returnMutated
                                    doSteepest
                                    doRandomOrder
                                    (curBestGraphList, curBestGraphCost)
                                    (newGraphList, newGraphCost)
                                    inSimAnnealParams
                                    netNodes
                                    currentCost
                                    otherPhyloGraphs
                            -- simulated annealing--needs new SAParams
                            Just _ → do
                                (saGraphs, saCounter) ←
                                    postProcessNetworkAddSA
                                        inGS
                                        inData
                                        maxNetEdges
                                        numToKeep
                                        counter
                                        returnMutated
                                        doSteepest
                                        doRandomOrder
                                        (curBestGraphList, curBestGraphCost)
                                        (newGraphList, newGraphCost)
                                        newSAParams
                                        netNodes
                                        currentCost
                                        otherPhyloGraphs

                                -- if want mutated then return
                                case returnMutated of
                                    True → pure (saGraphs, saCounter)
                                    -- delete non-minimal edges if any
                                    -- not sim anneal/drift regular optimal searching
                                    _ → do
                                        let annealBestCost = minimum $ fmap snd5 saGraphs
                                        (bestList', counter') ←
                                            deleteAllNetEdges'
                                                inGS
                                                inData
                                                maxNetEdges
                                                numToKeep
                                                saCounter
                                                False
                                                doSteepest
                                                doRandomOrder
                                                (saGraphs, annealBestCost)
                                                Nothing
                                                saGraphs
                                        bestList ← GO.selectGraphs Best numToKeep 0 bestList'
                                        pure (bestList, counter')
                    _ → do
                        logWith LogInfo $ unwords ["Maximum number of network edges reached:", show $ length netNodes, "\n"]
                        pure (take numToKeep curBestGraphList, counter)


{- | postProcessNetworkAddSA processes simaneal/drift
assumes SAParams are updated during return of graph list above
-}
postProcessNetworkAddSA
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Int
    → Int
    → Bool
    → Bool
    → Bool
    → ([ReducedPhylogeneticGraph], VertexCost)
    → ([ReducedPhylogeneticGraph], VertexCost)
    → Maybe SAParams
    → [LG.LNode VertexInfo]
    → VertexCost
    → [ReducedPhylogeneticGraph]
    → PhyG ([ReducedPhylogeneticGraph], Int)
postProcessNetworkAddSA inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (newGraphList, newGraphCost) inSimAnnealParams _ _ inPhyloGraphList =
    -- trace ("\t\tNumber of network edges: " <> (show $ length netNodes)) (
    -- this to deal with empty list issues if nothing found
    let (nextNewGraphList, firstNewGraphList) =
            if (not . null) newGraphList
                then (tail newGraphList, [head newGraphList])
                else ([], [])
        graphsToInsert =
            if doSteepest
                then newGraphList
                else take numToKeep $ newGraphList <> inPhyloGraphList
    in  -- always accept if found better
        if newGraphCost < curBestGraphCost
            then do
                logWith LogInfo ("\t-> " <> (show newGraphCost))
                insertAllNetEdges'
                    inGS
                    inData
                    maxNetEdges
                    numToKeep
                    (counter + 1)
                    returnMutated
                    doSteepest
                    doRandomOrder
                    (newGraphList, newGraphCost)
                    inSimAnnealParams
                    graphsToInsert
            else -- check if hit max change/ cooling steps

                if ((currentStep $ fromJust inSimAnnealParams) >= (numberSteps $ fromJust inSimAnnealParams))
                    || ((driftChanges $ fromJust inSimAnnealParams) >= (driftMaxChanges $ fromJust inSimAnnealParams))
                    then GO.selectGraphs Unique numToKeep 0 (newGraphList <> curBestGraphList) <&> \x → (x, counter)
                    else -- more to do

                        let annealBestCost = min curBestGraphCost newGraphCost
                        in  do
                                insertAllNetEdges'
                                    inGS
                                    inData
                                    maxNetEdges
                                    numToKeep
                                    (counter + 1)
                                    returnMutated
                                    doSteepest
                                    doRandomOrder
                                    (firstNewGraphList <> curBestGraphList, annealBestCost)
                                    inSimAnnealParams
                                    (nextNewGraphList <> inPhyloGraphList)


-- )

-- | postProcessNetworkAdd prcesses non-simaneal/drift--so no updating of SAParams
postProcessNetworkAdd
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Int
    → Int
    → Bool
    → Bool
    → Bool
    → ([ReducedPhylogeneticGraph], VertexCost)
    → ([ReducedPhylogeneticGraph], VertexCost)
    → Maybe SAParams
    → [LG.LNode VertexInfo]
    → VertexCost
    → [ReducedPhylogeneticGraph]
    → PhyG ([ReducedPhylogeneticGraph], Int)
postProcessNetworkAdd inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, _) (newGraphList, newGraphCost) inSimAnnealParams _ currentCost inPhyloGraphList = case newGraphCost `compare` currentCost of
    -- "steepest style descent" abandons existing list if better cost found

    -- check if graph OK--done in insert function
    LT →
        let graphsToInsert
                | doSteepest = newGraphList
                | otherwise = take numToKeep $ newGraphList <> inPhyloGraphList
        in  do
                logWith LogInfo ("\t-> " <> show newGraphCost)
                insertAllNetEdges'
                    inGS
                    inData
                    maxNetEdges
                    numToKeep
                    (counter + 1)
                    returnMutated
                    doSteepest
                    doRandomOrder
                    (newGraphList, newGraphCost)
                    inSimAnnealParams
                    graphsToInsert

    -- worse graphs found--go on
    GT →
        insertAllNetEdges'
            inGS
            inData
            maxNetEdges
            numToKeep
            (counter + 1)
            returnMutated
            doSteepest
            doRandomOrder
            (curBestGraphList, currentCost)
            inSimAnnealParams
            inPhyloGraphList
    -- equal cost
    -- not sure if should add new graphs to queue to do edge deletion again
    EQ → do
        -- new graph list contains the input graph if equal and filterd unique already in insertAllNetEdges
        newCurSameBestList ← GO.selectGraphs Unique numToKeep 0 $ curBestGraphList <> newGraphList
        insertAllNetEdges'
            inGS
            inData
            maxNetEdges
            numToKeep
            (counter + 1)
            returnMutated
            doSteepest
            doRandomOrder
            (newCurSameBestList, currentCost)
            inSimAnnealParams
            inPhyloGraphList


{- | insertEachNetEdge takes a phylogenetic graph and inserts all permissible network edges one at time
and returns unique list of new Phylogenetic Graphs and cost
even if worse--could be used for simulated annealing later
if equal returns unique graph list
-}
insertEachNetEdge
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Int
    → Bool
    → Bool
    → Maybe VertexCost
    → Maybe SAParams
    → ReducedPhylogeneticGraph
    → PhyG ([ReducedPhylogeneticGraph], VertexCost, Maybe SAParams)
insertEachNetEdge inGS inData maxNetEdges numToKeep doSteepest doRandomOrder preDeleteCost inSimAnnealParams inPhyloGraph =
    if LG.isEmpty $ fst5 inPhyloGraph
        then error "Empty input insertEachNetEdge graph in deleteAllNetEdges"
        else
            let currentCost =
                    if isNothing preDeleteCost
                        then snd5 inPhyloGraph
                        else fromJust preDeleteCost

                (_, _, _, netNodes) = LG.splitVertexList (thd5 inPhyloGraph)

                -- parallel stuff
                action ∷ (LG.LEdge b, LG.LEdge b) → PhyG ReducedPhylogeneticGraph
                action = insertNetEdge inGS inData inPhyloGraph preDeleteCost
            in  do
                    candidateNetworkEdgeList' ← getPermissibleEdgePairs inGS (thd5 inPhyloGraph)

                    -- radomize pair list
                    candidateNetworkEdgeList ←
                        if doRandomOrder
                            then shuffleList candidateNetworkEdgeList'
                            else pure candidateNetworkEdgeList'

                    inNetEdRList ←
                        insertNetEdgeRecursive
                            inGS
                            inData
                            maxNetEdges
                            doSteepest
                            doRandomOrder
                            inPhyloGraph
                            preDeleteCost
                            inSimAnnealParams
                            candidateNetworkEdgeList

                    inNetEdRListMAP ←
                        getParallelChunkTraverse >>= \pTraverse →
                            action `pTraverse` candidateNetworkEdgeList

                    let (newGraphList, newSAParams) =
                            if not doSteepest
                                then
                                    let genNewSimAnnealParams =
                                            if isNothing inSimAnnealParams
                                                then Nothing
                                                else U.incrementSimAnnealParams inSimAnnealParams
                                    in  -- TODO
                                        (filter (/= emptyReducedPhylogeneticGraph) inNetEdRListMAP, genNewSimAnnealParams)
                                else inNetEdRList
                    let minCost =
                            if null candidateNetworkEdgeList || null newGraphList
                                then infinity
                                else minimum $ fmap snd5 newGraphList

                    logWith LogInfo ("\tExamining at most " <> (show $ length candidateNetworkEdgeList) <> " candidate edge pairs" <> "\n")

                    -- no network edges to insert
                    -- trace ("IENE: " <> (show minCost)) (
                    if (length netNodes >= maxNetEdges)
                        then do
                            logWith LogInfo ("Maximum number of network edges reached: " <> (show $ length netNodes) <> "\n")
                            pure ([inPhyloGraph], snd5 inPhyloGraph, inSimAnnealParams)
                        else -- no edges to add

                            if null candidateNetworkEdgeList
                                then do
                                    -- trace ("IENE num cand edges:" <> (show $ length candidateNetworkEdgeList))
                                    pure ([inPhyloGraph], currentCost, newSAParams)
                                else -- single if steepest so no need to unique

                                    if doSteepest
                                        then GO.selectGraphs Best numToKeep 0 newGraphList <&> \x → (x, minCost, newSAParams)
                                        else -- "all" option needs to recurse since does all available edges at each step
                                        -- logic is here since not in the deleteNetEdge function

                                            if isNothing inSimAnnealParams
                                                then
                                                    let -- parallel stuff
                                                        insertAction
                                                            ∷ (Maybe SAParams, ReducedPhylogeneticGraph) → PhyG ([ReducedPhylogeneticGraph], VertexCost, Maybe SAParams)
                                                        insertAction = insertEachNetEdge' inGS inData maxNetEdges numToKeep doSteepest doRandomOrder preDeleteCost
                                                    in  if minCost < currentCost
                                                            then do
                                                                let annealParamList = replicate (length newGraphList) newSAParams
                                                                --
                                                                insertResult ←
                                                                    getParallelChunkTraverse >>= \pTraverse →
                                                                        pTraverse insertAction $ zip annealParamList newGraphList
                                                                let (allGraphListList, costList, allSAParamList) = unzip3 insertResult
                                                                (allMinCost, allMinCostGraphs) ← case fold allGraphListList of
                                                                    [] → pure (infinity, [])
                                                                    gs → GO.selectGraphs Unique numToKeep 0.0 gs <&> \x → (minimum costList, x)

                                                                pure (allMinCostGraphs, allMinCost, U.incrementSimAnnealParams $ head allSAParamList)
                                                            else GO.selectGraphs Unique numToKeep 0 newGraphList <&> \x → (x, minCost, newSAParams)
                                                else -- SA anneal/Drift

                                                -- always take better

                                                    if minCost < currentCost
                                                        then do
                                                            pure (newGraphList, minCost, newSAParams)
                                                        else -- check if hit step limit--more for SA than drift

                                                            if ((currentStep $ fromJust inSimAnnealParams) >= (numberSteps $ fromJust inSimAnnealParams))
                                                                || ((driftChanges $ fromJust inSimAnnealParams) >= (driftMaxChanges $ fromJust inSimAnnealParams))
                                                                then do
                                                                    pure ([inPhyloGraph], snd5 inPhyloGraph, inSimAnnealParams)
                                                                else -- otherwise do the anneal/Drift accept, or keep going on input graph

                                                                do
                                                                    (acceptGraph, nextSAParams) ← U.simAnnealAccept inSimAnnealParams currentCost minCost
                                                                    case acceptGraph of
                                                                        True → pure (newGraphList, minCost, newSAParams)
                                                                        _ →
                                                                            insertEachNetEdge
                                                                                inGS
                                                                                inData
                                                                                maxNetEdges
                                                                                numToKeep
                                                                                doSteepest
                                                                                doRandomOrder
                                                                                preDeleteCost
                                                                                nextSAParams
                                                                                inPhyloGraph


-- | insertEachNetEdge' is a wrapper around insertEachNetEdge to allow for parmapping with multiple parameters
insertEachNetEdge'
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Int
    → Bool
    → Bool
    → Maybe VertexCost
    → (Maybe SAParams, ReducedPhylogeneticGraph)
    → PhyG ([ReducedPhylogeneticGraph], VertexCost, Maybe SAParams)
insertEachNetEdge' inGS inData maxNetEdges numToKeep doSteepest doRandomOrder preDeleteCost (inSimAnnealParams, inPhyloGraph) =
    insertEachNetEdge inGS inData maxNetEdges numToKeep doSteepest doRandomOrder preDeleteCost inSimAnnealParams inPhyloGraph


{- | insertNetEdgeRecursive recursively inserts edges and returns new graph only if better
if parallel evaluated numThreads each time (steepest scenario)
-}
insertNetEdgeRecursive
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Bool
    → Bool
    → ReducedPhylogeneticGraph
    → Maybe VertexCost
    → Maybe SAParams
    → [(LG.LEdge EdgeInfo, LG.LEdge EdgeInfo)]
    → PhyG ([ReducedPhylogeneticGraph], Maybe SAParams)
insertNetEdgeRecursive inGS inData maxNetEdges doSteepest doRandomOrder inPhyloGraph preDeleteCost inSimAnnealParams inEdgePairList =
    -- trace ("Edges pairs to go : " <> (show $ length edgePairList)) (
    if null inEdgePairList
        then do
            pure ([inPhyloGraph], inSimAnnealParams)
        else -- don't want to over saturate the parallel thread system

            let {-saRounds = if isNothing inSimAnnealParams then 1
                           else rounds $ fromJust inSimAnnealParams
                (numGraphsToExamine, _) = divMod PU.getNumThreads saRounds -- this may not "drift" if finds alot better, but that's how its supposed to work
                -}

                numGraphsToExamine = graphsSteepest inGS -- min (graphsSteepest inGS) PU.getNumThreads
                -- firstEdgePair = head edgePairList
                edgePairList = take numGraphsToExamine inEdgePairList

                -- check for max net edges
                (_, _, _, netNodes) = LG.splitVertexList (thd5 inPhyloGraph)

                -- parallel seup
                action ∷ (LG.LEdge b, LG.LEdge b) → PhyG ReducedPhylogeneticGraph
                action = insertNetEdge inGS inData inPhyloGraph preDeleteCost
            in  do
                    -- need to check display/character trees not conical graph
                    -- these graph costs are "exact" or at least non-heuristic--needs to be updated when get a good heuristic
                    newGraphList'' ←
                        getParallelChunkTraverse >>= \pTraverse →
                            action `pTraverse` edgePairList

                    let newGraphList' = filter (/= emptyReducedPhylogeneticGraph) newGraphList''
                    newGraphList ← GO.selectGraphs Best (maxBound ∷ Int) 0 newGraphList'
                    let newGraphCost = snd5 $ head newGraphList

                    if length netNodes >= maxNetEdges
                        then do
                            logWith LogInfo $ unwords ["Maximum number of network edges reached:", show $ length netNodes, "\n"]
                            pure ([inPhyloGraph], inSimAnnealParams)
                        else -- malformed graph--returns nothing for either regular or simAnneal/drift

                            if null newGraphList'
                                then do
                                    -- trace ("INER: Empty more to go : " <> (show $ length $ tail edgePairList))
                                    insertNetEdgeRecursive
                                        inGS
                                        inData
                                        maxNetEdges
                                        doSteepest
                                        doRandomOrder
                                        inPhyloGraph
                                        preDeleteCost
                                        inSimAnnealParams
                                        (drop numGraphsToExamine inEdgePairList)
                                else -- "regular" insert, within steepest

                                    if isNothing inSimAnnealParams
                                        then -- better cost

                                            if newGraphCost < snd5 inPhyloGraph
                                                then do
                                                    -- cyclic check in insert edge function
                                                    -- trace ("INER: Better -> " <> (show $ snd5 newGraph))
                                                    pure (newGraphList, inSimAnnealParams)
                                                else -- not better
                                                do
                                                    -- trace ("INER: Really Not Better")
                                                    insertNetEdgeRecursive
                                                        inGS
                                                        inData
                                                        maxNetEdges
                                                        doSteepest
                                                        doRandomOrder
                                                        inPhyloGraph
                                                        preDeleteCost
                                                        inSimAnnealParams
                                                        (drop numGraphsToExamine inEdgePairList)
                                        else -- sim annealing/drift

                                        -- trace ("IENR:" <> (show (newGraphCost, snd5 inPhyloGraph)) <> " params: " <> (show (currentStep $ fromJust inSimAnnealParams, numberSteps $ fromJust inSimAnnealParams, driftChanges $ fromJust inSimAnnealParams, driftMaxChanges $ fromJust inSimAnnealParams))) (
                                        -- if better always accept

                                            if newGraphCost < snd5 inPhyloGraph
                                                then do
                                                    -- cyclic check in insert edge function
                                                    -- these graph costs are "exact" or at least non-heuristic--needs to be updated when get a good heuristic
                                                    (_, nextSAParams) ← U.simAnnealAccept inSimAnnealParams (snd5 inPhyloGraph) newGraphCost
                                                    pure (newGraphList, nextSAParams)
                                                else -- check if hit step limit--more for SA than drift

                                                    if ((currentStep $ fromJust inSimAnnealParams) >= (numberSteps $ fromJust inSimAnnealParams))
                                                        || ((driftChanges $ fromJust inSimAnnealParams) >= (driftMaxChanges $ fromJust inSimAnnealParams))
                                                        then do
                                                            pure ([inPhyloGraph], inSimAnnealParams)
                                                        else do
                                                            -- otherwise do the anneal/Drift accept

                                                            (acceptGraph, nextSAParams) ← U.simAnnealAccept inSimAnnealParams (snd5 inPhyloGraph) newGraphCost
                                                            case acceptGraph of
                                                                True → pure (newGraphList, nextSAParams)
                                                                _ →
                                                                    insertNetEdgeRecursive
                                                                        inGS
                                                                        inData
                                                                        maxNetEdges
                                                                        doSteepest
                                                                        doRandomOrder
                                                                        inPhyloGraph
                                                                        preDeleteCost
                                                                        nextSAParams
                                                                        (drop numGraphsToExamine inEdgePairList)


{- | insertNetEdge inserts an edge between two other edges, creating 2 new nodes and rediagnoses graph
contacts deletes 2 orginal edges and adds 2 nodes and 5 new edges
does not check any edge reasonable-ness properties
new edge directed from first to second edge
naive for now
predeletecost of edge move
no choice of graph--just makes and returns
-}
insertNetEdge
    ∷ GlobalSettings
    → ProcessedData
    → ReducedPhylogeneticGraph
    → Maybe VertexCost
    → (LG.LEdge b, LG.LEdge b)
    → PhyG ReducedPhylogeneticGraph
insertNetEdge inGS inData inPhyloGraph _ edgePair@((u, v, _), (u', v', _)) =
    if LG.isEmpty $ thd5 inPhyloGraph
        then error "Empty input phylogenetic graph in insNetEdge"
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

                -- do not prune other edges if now unused
                pruneEdges = False

                -- don't warn that edges are being pruned
                warnPruneEdges = False

                -- graph optimization from root
                startVertex = Nothing

                -- full two-pass optimization
                leafGraph = LG.extractLeafGraph $ thd5 inPhyloGraph
            in  do
                    newPhyloGraph ← T.multiTraverseFullyLabelSoftWiredReduced inGS inData pruneEdges warnPruneEdges leafGraph startVertex newSimple

                    -- calculates heursitic graph delta
                    -- (heuristicDelta, _, _, _, _)  = heuristicAddDelta inGS inPhyloGraph edgePair (fst newNodeOne) (fst newNodeTwo)
                    heuristicDelta' <- heuristicAddDelta' inGS inPhyloGraph edgePair

                    let edgeAddDelta = deltaPenaltyAdjustment inGS inPhyloGraph "add"

                    let heuristicFactor = (heuristicDelta' + edgeAddDelta) / edgeAddDelta

                    -- use or not Net add heuristics
                    let metHeuristicThreshold = not (useNetAddHeuristic inGS) || heuristicFactor < (2 / 3)

                    -- remove these checks when working
                    isPhyloGraph ← LG.isPhylogeneticGraph newSimple
                    if not isPhyloGraph
                        then do
                            pure emptyReducedPhylogeneticGraph
                        else
                            if metHeuristicThreshold
                                then -- if (GO.parentsInChainGraph . thd5) newPhyloGraph then emptyPhylogeneticGraph
                                -- else

                                    if (snd5 newPhyloGraph <= snd5 inPhyloGraph)
                                        then do
                                            pure newPhyloGraph
                                        else do
                                            pure emptyReducedPhylogeneticGraph
                                else do
                                    pure emptyReducedPhylogeneticGraph


-- | (curBestGraphList, annealBestCost) is a wrapper for moveAllNetEdges' allowing for multiple simulated annealing rounds
deleteAllNetEdges
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Int
    → Int
    → Bool
    → Bool
    → Bool
    → ([ReducedPhylogeneticGraph], VertexCost)
    → (Maybe SAParams, [ReducedPhylogeneticGraph])
    → PhyG ([ReducedPhylogeneticGraph], Int)
deleteAllNetEdges inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
    if isNothing inSimAnnealParams
        then
            deleteAllNetEdges'
                inGS
                inData
                maxNetEdges
                numToKeep
                counter
                returnMutated
                doSteepest
                doRandomOrder
                (curBestGraphList, curBestGraphCost)
                inSimAnnealParams
                inPhyloGraphList
        else
            let -- create list of params with unique list of random values for rounds of annealing
                annealingRounds = rounds $ fromJust inSimAnnealParams
                annealParamGraphList = replicate annealingRounds inSimAnnealParams

                -- parallel
                action ∷ (Maybe SAParams, [ReducedPhylogeneticGraph]) → PhyG ([ReducedPhylogeneticGraph], Int)
                action =
                    deleteAllNetEdges''
                        inGS
                        inData
                        maxNetEdges
                        numToKeep
                        counter
                        returnMutated
                        doSteepest
                        doRandomOrder
                        (curBestGraphList, curBestGraphCost)
            in  do
                    deleteResult ←
                        getParallelChunkTraverse >>= \pTraverse →
                            pTraverse action . zip annealParamGraphList $ replicate annealingRounds inPhyloGraphList

                    let (annealRoundsList, counterList) = unzip deleteResult
                    GO.selectGraphs Best numToKeep 0 (fold annealRoundsList) <&> \x → (x, sum counterList)


-- | deleteAllNetEdges'' is a wrapper around deleteAllNetEdges' to allow use of seqParMap
deleteAllNetEdges''
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Int
    → Int
    → Bool
    → Bool
    → Bool
    → ([ReducedPhylogeneticGraph], VertexCost)
    → (Maybe SAParams, [ReducedPhylogeneticGraph])
    → PhyG ([ReducedPhylogeneticGraph], Int)
deleteAllNetEdges'' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) (inSimAnnealParams, inPhyloGraphList) =
    deleteAllNetEdges'
        inGS
        inData
        maxNetEdges
        numToKeep
        counter
        returnMutated
        doSteepest
        doRandomOrder
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
    → Int
    → Int
    → Int
    → Bool
    → Bool
    → Bool
    → ([ReducedPhylogeneticGraph], VertexCost)
    → Maybe SAParams
    → [ReducedPhylogeneticGraph]
    → PhyG ([ReducedPhylogeneticGraph], Int)
deleteAllNetEdges' inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) inSimAnnealParams = \case
    [] → pure (take numToKeep curBestGraphList, counter)
    firstPhyloGraph : otherPhyloGraphs
        | LG.isEmpty $ fst5 firstPhyloGraph →
            deleteAllNetEdges'
                inGS
                inData
                maxNetEdges
                numToKeep
                counter
                returnMutated
                doSteepest
                doRandomOrder
                (curBestGraphList, curBestGraphCost)
                inSimAnnealParams
                otherPhyloGraphs
    inputPhyloGraphs@(firstPhyloGraph : otherPhyloGraphs) → do
        let currentCost = min curBestGraphCost $ snd5 firstPhyloGraph

        (newGraphList', _, newSAParams) ←
            deleteEachNetEdge
                inGS
                inData
                numToKeep
                doSteepest
                doRandomOrder
                False
                inSimAnnealParams
                firstPhyloGraph

        newGraphList ← GO.selectGraphs Best numToKeep 0 newGraphList'
        let newGraphCost = case newGraphList of
                [] → infinity
                (_, c, _, _, _) : _ → c

        -- if graph is a tree no edges to delete
        case LG.isTree $ fst5 firstPhyloGraph of
            True → do
                logWith LogInfo "\tGraph in delete network edges is tree--skipping\n"
                deleteAllNetEdges'
                    inGS
                    inData
                    maxNetEdges
                    numToKeep
                    (counter + 1)
                    returnMutated
                    doSteepest
                    doRandomOrder
                    (firstPhyloGraph : curBestGraphList, currentCost)
                    inSimAnnealParams
                    otherPhyloGraphs

            -- is this an issue for SA?
            _ → case newGraphList of
                [] → pure (take numToKeep curBestGraphList, counter + 1)
                g : gs → case inSimAnnealParams of
                    -- regular delete wihtout simulated annealing
                    Nothing →
                        postProcessNetworkDelete
                            inGS
                            inData
                            maxNetEdges
                            numToKeep
                            counter
                            returnMutated
                            doSteepest
                            doRandomOrder
                            (curBestGraphList, curBestGraphCost)
                            inSimAnnealParams
                            inputPhyloGraphs
                            newGraphList
                            newGraphCost
                            currentCost
                    -- simulated annealing
                    Just simAnneal → do
                        (saGraphs, saCounter) ←
                            postProcessNetworkDeleteSA
                                inGS
                                inData
                                maxNetEdges
                                numToKeep
                                counter
                                returnMutated
                                doSteepest
                                doRandomOrder
                                (curBestGraphList, curBestGraphCost)
                                newSAParams
                                inputPhyloGraphs
                                newGraphList
                                newGraphCost
                                currentCost

                        case returnMutated of
                            -- if want mutated then return
                            True → pure (saGraphs, saCounter)
                            -- insert non-minimal edges if any
                            -- not sim anneal/drift regular optimal searching
                            False → do
                                let annealBestCost = minimum $ fmap snd5 saGraphs
                                insertedGraphs ←
                                    insertAllNetEdges'
                                        inGS
                                        inData
                                        maxNetEdges
                                        numToKeep
                                        saCounter
                                        False
                                        doSteepest
                                        doRandomOrder
                                        (saGraphs, annealBestCost)
                                        Nothing
                                        saGraphs

                                let (bestList', counter') = insertedGraphs
                                bestList ← GO.selectGraphs Best numToKeep 0 bestList'
                                pure (bestList, counter')


{- | postProcessNetworkDeleteSA postprocesses results from delete actions for non-annealing/Drift network delete operations
assumes SAParams are updated during return of graph list above
-}
postProcessNetworkDeleteSA
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Int
    → Int
    → Bool
    → Bool
    → Bool
    → ([ReducedPhylogeneticGraph], VertexCost)
    → Maybe SAParams
    → [ReducedPhylogeneticGraph]
    → [ReducedPhylogeneticGraph]
    → VertexCost
    → VertexCost
    → PhyG ([ReducedPhylogeneticGraph], Int)
postProcessNetworkDeleteSA inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, curBestGraphCost) inSimAnnealParams inPhyloGraphList newGraphList newGraphCost currentCost =
    -- this to deal with empty list issues if nothing found
    let (nextNewGraphList, firstNewGraphList) = case newGraphList of
            [] → ([], [])
            g : gs → (gs, [g])

        graphsToDelete
            | doSteepest = newGraphList
            | otherwise = take numToKeep $ newGraphList <> inPhyloGraphList
    in  -- always accept if found better
        if newGraphCost < currentCost
            then do
                logWith LogInfo ("\t-> " <> (show newGraphCost))
                deleteAllNetEdges'
                    inGS
                    inData
                    maxNetEdges
                    numToKeep
                    (counter + 1)
                    returnMutated
                    doSteepest
                    doRandomOrder
                    (newGraphList, newGraphCost)
                    inSimAnnealParams
                    graphsToDelete
            else -- check if hit max change/ cooling steps

                if ((currentStep $ fromJust inSimAnnealParams) >= (numberSteps $ fromJust inSimAnnealParams))
                    || ((driftChanges $ fromJust inSimAnnealParams) >= (driftMaxChanges $ fromJust inSimAnnealParams))
                    then GO.selectGraphs Unique numToKeep 0 (newGraphList <> curBestGraphList) <&> \x → (x, counter)
                    else -- more to do

                        let annealBestCost = min curBestGraphCost newGraphCost
                        in  do
                                deleteAllNetEdges'
                                    inGS
                                    inData
                                    maxNetEdges
                                    numToKeep
                                    (counter + 1)
                                    returnMutated
                                    doSteepest
                                    doRandomOrder
                                    (firstNewGraphList <> curBestGraphList, annealBestCost)
                                    inSimAnnealParams
                                    (nextNewGraphList <> inPhyloGraphList)


-- | postProcessNetworkDelete postprocesses results from delete actions for "regular" ie non-annealing/Drift network delete operations
postProcessNetworkDelete
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Int
    → Int
    → Bool
    → Bool
    → Bool
    → ([ReducedPhylogeneticGraph], VertexCost)
    → Maybe SAParams
    → [ReducedPhylogeneticGraph]
    → [ReducedPhylogeneticGraph]
    → VertexCost
    → VertexCost
    → PhyG ([ReducedPhylogeneticGraph], Int)
postProcessNetworkDelete inGS inData maxNetEdges numToKeep counter returnMutated doSteepest doRandomOrder (curBestGraphList, _) inSimAnnealParams inPhyloGraphList newGraphList newGraphCost currentCost =
    -- worse graphs found--go on
    if newGraphCost > currentCost
        then do
            deleteAllNetEdges'
                inGS
                inData
                maxNetEdges
                numToKeep
                (counter + 1)
                returnMutated
                doSteepest
                doRandomOrder
                ((head inPhyloGraphList) : curBestGraphList, currentCost)
                inSimAnnealParams
                (tail inPhyloGraphList)
        else -- "steepest style descent" abandons existing list if better cost found

            if newGraphCost < currentCost
                then do
                    logWith LogInfo ("\t-> " <> (show newGraphCost))
                    if doSteepest
                        then do
                            deleteAllNetEdges'
                                inGS
                                inData
                                maxNetEdges
                                numToKeep
                                (counter + 1)
                                returnMutated
                                doSteepest
                                doRandomOrder
                                (newGraphList, newGraphCost)
                                inSimAnnealParams
                                newGraphList
                        else do
                            deleteAllNetEdges'
                                inGS
                                inData
                                maxNetEdges
                                numToKeep
                                (counter + 1)
                                returnMutated
                                doSteepest
                                doRandomOrder
                                (newGraphList, newGraphCost)
                                inSimAnnealParams
                                (newGraphList <> (tail inPhyloGraphList))
                else -- equal cost
                -- not sure if should add new graphs to queue to do edge deletion again

                -- new graph list contains the input graph if equal and filterd unique already in deleteEachNetEdge
                do
                    newCurSameBestList ← GO.selectGraphs Unique numToKeep 0.0 $ curBestGraphList <> newGraphList
                    deleteAllNetEdges'
                        inGS
                        inData
                        maxNetEdges
                        numToKeep
                        (counter + 1)
                        returnMutated
                        doSteepest
                        doRandomOrder
                        (newCurSameBestList, currentCost)
                        inSimAnnealParams
                        (tail inPhyloGraphList)


-- | deleteOneNetAddAll' wrapper on deleteOneNetAddAll to allow for parmap
deleteOneNetAddAll'
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Int
    → Bool
    → Bool
    → ReducedPhylogeneticGraph
    → Maybe SAParams
    → LG.Edge
    → PhyG [ReducedPhylogeneticGraph]
deleteOneNetAddAll' inGS inData maxNetEdges numToKeep doSteepest doRandomOrder inPhyloGraph inSimAnnealParams edgeToDelete =
    deleteOneNetAddAll
        inGS
        inData
        maxNetEdges
        numToKeep
        doSteepest
        doRandomOrder
        inPhyloGraph
        [edgeToDelete]
        inSimAnnealParams


{- | deleteOneNetAddAll version deletes net edges in turn and readds-based on original cost
but this cost in graph (really not correct) but allows logic of insert edge to function better
unlike deleteOneNetAddAll' only deals with single edge deletion at a time
-}
deleteOneNetAddAll
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Int
    → Bool
    → Bool
    → ReducedPhylogeneticGraph
    → [LG.Edge]
    → Maybe SAParams
    → PhyG [ReducedPhylogeneticGraph]
deleteOneNetAddAll inGS inData maxNetEdges numToKeep doSteepest doRandomOrder inPhyloGraph edgeToDeleteList inSimAnnealParams =
    if null edgeToDeleteList
        then do
            -- trace ("\tGraph has no edges to move---skipping")
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
                                maxNetEdges
                                numToKeep
                                doSteepest
                                doRandomOrder
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

                                    insertedGraphTripleList ←
                                        insertEachNetEdge
                                            inGS
                                            inData
                                            (curNumNetNodes + 1)
                                            numToKeep
                                            doSteepest
                                            doRandomOrder
                                            Nothing
                                            inSimAnnealParams
                                            graphToInsert'

                                    let newMinimumCost = snd3 insertedGraphTripleList

                                    let newBestGraphs = filter ((== newMinimumCost) . snd5) $ fst3 insertedGraphTripleList

                                    -- trace ("DONAA-New: " <> (show (inGraphCost, fmap snd5 graphsToInsert, fmap snd5 graphsToInsert', newMinimumCost))) (
                                    if newMinimumCost < inGraphCost
                                        then do
                                            -- trace ("DONA-> ")
                                            pure newBestGraphs
                                        else do
                                            deleteOneNetAddAll
                                                inGS
                                                inData
                                                maxNetEdges
                                                numToKeep
                                                doSteepest
                                                doRandomOrder
                                                inPhyloGraph
                                                (tail edgeToDeleteList)
                                                inSimAnnealParams


{- | getPermissibleEdgePairs takes a DecoratedGraph and returns the list of all pairs
of edges that can be joined by a network edge and meet all necessary conditions
-}

-- add in other conditions
--   reproducable--ie not tree noide with two net node children--other stuff
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
                    -- PU.seqParMap (parStrategy $ lazyParStrat inGS) (isEdgePairPermissible inGraph coevalNodeConstraintList') edgePairs -- `using`  PU.myParListChunkRDS

                    let pairList = fmap fst $ filter ((== True) . snd) $ zip edgePairs edgeTestList

                    -- trace ("Edge Pair list :" <> (show $ fmap f pairList) <> "\n"
                    --  <> "GPEP\n" <> (LG.prettify $ GO.convertDecoratedToSimpleGraph inGraph))
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

{- These heuristics do not seem tom work well at all-}

{- | heuristic add delta' based on new display tree and delta from existing costs by block--assumming < 0
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
heuristicAddDelta' ∷ GlobalSettings → ReducedPhylogeneticGraph → (LG.LEdge b, LG.LEdge b) → PhyG VertexCost
heuristicAddDelta' _ inPhyloGraph ((u, v, _), (u', v', _)) =
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


-- if dUnionUVV' - dU'V' < dU'V' then dUnionUVV' - dU'V'
-- else 0.0
-- )

{- | heuristicAddDelta takes the existing graph, edge pair, and new nodes to create and makes
the new nodes and reoptimizes starting nodes of two edges.  Returns cost delta based on
previous and new node resolution caches
returns cost delta and the reoptimized nodes for use in incremental optimization
original edges (to be deleted) (u,v) and (u',v'), n1 inserted in (u,v) and n2 inserted into (u',v')
creates (n1, n2), (u,n1), (n1,v), (u',n2), (n2, v')
-}
heuristicAddDelta
    ∷ GlobalSettings
    → ReducedPhylogeneticGraph
    → (LG.LEdge b, LG.LEdge b)
    → LG.Node
    → LG.Node
    → PhyG (VertexCost, LG.LNode VertexInfo, LG.LNode VertexInfo, LG.LNode VertexInfo, LG.LNode VertexInfo)
heuristicAddDelta inGS inPhyloGraph ((u, v, _), (u', v', _)) n1 n2 =
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

{- | deleteEachNetEdge takes a phylogenetic graph and deletes all network edges one at time
and returns best list of new Phylogenetic Graphs and cost
even if worse--could be used for simulated annealing later
if equal returns unique graph list
-}
deleteEachNetEdge
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Bool
    → Bool
    → Bool
    → Maybe SAParams
    → ReducedPhylogeneticGraph
    → PhyG ([ReducedPhylogeneticGraph], VertexCost, Maybe SAParams)
deleteEachNetEdge inGS inData numToKeep doSteepest doRandomOrder force inSimAnnealParams inPhyloGraph =
    -- trace ("DENE start") (
    if LG.isEmpty $ thd5 inPhyloGraph
        then do
            pure ([], infinity, inSimAnnealParams) -- error "Empty input phylogenetic graph in deleteAllNetEdges"
        else
            let currentCost = snd5 inPhyloGraph

                -- potentially randomize order of list
                networkEdgeList' = LG.netEdges $ thd5 inPhyloGraph

                --- parallel
                action ∷ LG.Edge → PhyG ReducedPhylogeneticGraph
                action = deleteNetEdge inGS inData inPhyloGraph force
            in  do
                    networkEdgeList ←
                        if not doRandomOrder
                            then pure networkEdgeList'
                            else shuffleList networkEdgeList'

                    delNetEdgeList ←
                        getParallelChunkTraverse >>= \pTraverse →
                            action `pTraverse` networkEdgeList

                    deleteNetEdgeRecursiveList ← deleteNetEdgeRecursive inGS inData inPhyloGraph force inSimAnnealParams networkEdgeList
                    let (newGraphList, newSAParams) =
                            if not doSteepest
                                then (delNetEdgeList, U.incrementSimAnnealParams inSimAnnealParams)
                                else deleteNetEdgeRecursiveList

                    bestCostGraphList ← filter ((/= infinity) . snd5) <$> GO.selectGraphs Best numToKeep 0 newGraphList
                    let minCost =
                            if null bestCostGraphList
                                then infinity
                                else minimum $ fmap snd5 bestCostGraphList

                    -- no network edges to delete
                    if null networkEdgeList
                        then do
                            logWith LogInfo ("\tNo network edges to delete" <> "\n")
                            pure ([inPhyloGraph], currentCost, inSimAnnealParams)
                        else -- single if steepest so no neeed to unique--and have run through all options (including SA stuff) via recursive call

                            if doSteepest
                                then do
                                    pure (newGraphList, minCost, newSAParams)
                                else -- "all" option needs to recurse since does all available edges at each step
                                -- logic is here since not in the deleteNetEdge function

                                    if isNothing inSimAnnealParams
                                        then
                                            if minCost < currentCost
                                                then -- trace ("DENE--Delete net edge return:" <> (show (minCost,length uniqueCostGraphList))) (

                                                    let annealParamList = replicate (length bestCostGraphList) newSAParams

                                                        -- parallel
                                                        deleteAction ∷ (Maybe SAParams, ReducedPhylogeneticGraph) → PhyG ([ReducedPhylogeneticGraph], VertexCost, Maybe SAParams)
                                                        deleteAction = deleteEachNetEdge' inGS inData numToKeep doSteepest doRandomOrder force
                                                    in  do
                                                            -- TODO
                                                            nextGraphDoubleList ←
                                                                getParallelChunkTraverse >>= \pTraverse →
                                                                    deleteAction `pTraverse` zip annealParamList bestCostGraphList

                                                            let newMinCost = minimum $ fmap snd3 nextGraphDoubleList
                                                            let newGraphListBetter = filter ((== newMinCost) . snd5) $ concatMap fst3 nextGraphDoubleList

                                                            GO.selectGraphs Unique numToKeep 0 newGraphListBetter <&> \x → (x, newMinCost, newSAParams)
                                                else do
                                                    pure (bestCostGraphList, currentCost, newSAParams)
                                        else -- SA anneal/Drift

                                        -- always take better

                                            if minCost < currentCost
                                                then do
                                                    pure (bestCostGraphList, minCost, newSAParams)
                                                else -- check if hit step limit--more for SA than drift

                                                    if ((currentStep $ fromJust inSimAnnealParams) >= (numberSteps $ fromJust inSimAnnealParams))
                                                        || ((driftChanges $ fromJust inSimAnnealParams) >= (driftMaxChanges $ fromJust inSimAnnealParams))
                                                        then do
                                                            pure ([inPhyloGraph], snd5 inPhyloGraph, inSimAnnealParams)
                                                        else do
                                                            -- otherwise do the anneal/Drift accept, or keep going on input graph
                                                            (acceptGraph, nextSAParams) ← U.simAnnealAccept inSimAnnealParams currentCost minCost
                                                            case acceptGraph of
                                                                True → pure (bestCostGraphList, minCost, newSAParams)
                                                                _ → deleteEachNetEdge inGS inData numToKeep doSteepest doRandomOrder force nextSAParams inPhyloGraph


{- | deleteEachNetEdge' is a wrapper around deleteEachNetEdge to allow for zipping new random seeds for each
replicate
-}
deleteEachNetEdge'
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Bool
    → Bool
    → Bool
    → (Maybe SAParams, ReducedPhylogeneticGraph)
    → PhyG ([ReducedPhylogeneticGraph], VertexCost, Maybe SAParams)
deleteEachNetEdge' inGS inData numToKeep doSteepest doRandomOrder force (inSimAnnealParams, inPhyloGraph) =
    deleteEachNetEdge inGS inData numToKeep doSteepest doRandomOrder force inSimAnnealParams inPhyloGraph


{- | deleteNetEdgeRecursive like deleteEdge, deletes an edge (checking if network) and rediagnoses graph
contacts in=out=1 edges and removes node, reindexing nodes and edges
except returns on first better (as opposed to do all deletes first)
or sim annleal/drift
-}
deleteNetEdgeRecursive
    ∷ GlobalSettings
    → ProcessedData
    → ReducedPhylogeneticGraph
    → Bool
    → Maybe SAParams
    → [LG.Edge]
    → PhyG ([ReducedPhylogeneticGraph], Maybe SAParams)
deleteNetEdgeRecursive inGS inData inPhyloGraph force inSimAnnealParams inEdgeToDeleteList =
    if null inEdgeToDeleteList
        then do
            pure ([], inSimAnnealParams)
        else
            let {- Unclear if should adjust to number of rounds if already limiting to graphsSteepest value
                 saRounds = if isNothing inSimAnnealParams then 1
                           else rounds $ fromJust inSimAnnealParams

                 (numGraphsToExamine, _) = divMod PU.getNumThreads saRounds -- this may not "drift" if finds alot better, but that's how its supposed to work
                -}
                numGraphsToExamine = graphsSteepest inGS -- min (graphsSteepest inGS) PU.getNumThreads
                -- edgeToDelete = head inEdgeToDeleteList
                edgeToDeleteList = take numGraphsToExamine inEdgeToDeleteList

                leafGraph = LG.extractLeafGraph $ thd5 inPhyloGraph

                -- prune other edges if now unused
                pruneEdges = False

                -- don't warn that edges are being pruned
                warnPruneEdges = False

                -- graph optimization from root
                startVertex = Nothing

                -- parallel
                deleteAction ∷ LG.Edge → PhyG (SimpleGraph, Bool)
                deleteAction = deleteNetworkEdge (fst5 inPhyloGraph)

                softTraverse ∷ SimpleGraph → PhyG ReducedPhylogeneticGraph
                softTraverse = T.multiTraverseFullyLabelSoftWiredReduced inGS inData pruneEdges warnPruneEdges leafGraph startVertex

                hardTraverse ∷ SimpleGraph → PhyG ReducedPhylogeneticGraph
                hardTraverse = T.multiTraverseFullyLabelHardWiredReduced inGS inData leafGraph startVertex
            in  do
                    -- calls general funtion to remove network graph edge
                    simpleGraphList' ←
                        getParallelChunkTraverse >>= \pTraverse →
                            deleteAction `pTraverse` edgeToDeleteList

                    let simpleGraphList = fmap fst $ filter ((== True) . snd) simpleGraphList'
                    -- \$ PU.seqParMap (parStrategy $ lazyParStrat inGS) (deleteNetworkEdge (fst5 inPhyloGraph)) edgeToDeleteList

                    -- delSimple = GO.contractIn1Out1EdgesRename $ LG.delEdge edgeToDelete $ fst5 inPhyloGraph
                    -- full two-pass optimization
                    let leafGraph = LG.extractLeafGraph $ thd5 inPhyloGraph

                    -- (heuristicDelta, _, _) = heuristicDeleteDelta inGS inPhyloGraph edgeToDelete
                    -- heuristicDelta = 0.0

                    -- can treat as negative for delete
                    -- edgeAddDelta = deltaPenaltyAdjustment inGS inPhyloGraph "delete"

                    newPhyloGraphList' ← case graphType inGS of
                        SoftWired →
                            getParallelChunkTraverse >>= \pTraverse →
                                softTraverse `pTraverse` simpleGraphList
                        HardWired →
                            getParallelChunkTraverse >>= \pTraverse →
                                hardTraverse `pTraverse` simpleGraphList
                        val →
                            failWithPhase Computing $
                                fold
                                    ["Unsupported graph type '", show val, "' in deleteNetEdge. Must be soft- or hard-wired"]

                    newPhyloGraphList ← GO.selectGraphs Best (maxBound ∷ Int) 0 newPhyloGraphList'

                    -- if not modified return original graph
                    -- This check seems to be issue with delete not functioning properly
                    if null simpleGraphList
                        then do
                            pure ([inPhyloGraph], inSimAnnealParams)
                        else -- forcing delete for move

                            if force
                                then do
                                    -- trace ("DNERec forced")
                                    pure (newPhyloGraphList, inSimAnnealParams)
                                else -- regular search not sim anneal/drift

                                    if (isNothing inSimAnnealParams)
                                        then -- return if better

                                            if (snd5 $ head newPhyloGraphList) < (snd5 inPhyloGraph)
                                                then do
                                                    -- trace  ("DNERec Better -> " <> (show $ snd5 newPhyloGraph))
                                                    pure (newPhyloGraphList, inSimAnnealParams)
                                                else do
                                                    -- need to update edge list for new graph
                                                    -- potentially randomize order of list
                                                    deleteNetEdgeRecursive inGS inData inPhyloGraph force inSimAnnealParams (drop numGraphsToExamine inEdgeToDeleteList)
                                        else -- sim anneal/drift

                                        -- if better always accept

                                            if (snd5 $ head newPhyloGraphList) < (snd5 inPhyloGraph)
                                                then do
                                                    -- these graph costs are "exact" or at least non-heuristic--needs to be updated when get a good heuristic
                                                    (_, nextSAParams) ← U.simAnnealAccept inSimAnnealParams (snd5 inPhyloGraph) . snd5 $ head newPhyloGraphList
                                                    pure (newPhyloGraphList, nextSAParams)
                                                else -- check if hit step limit--more for SA than drift

                                                    if ((currentStep $ fromJust inSimAnnealParams) >= (numberSteps $ fromJust inSimAnnealParams))
                                                        || ((driftChanges $ fromJust inSimAnnealParams) >= (driftMaxChanges $ fromJust inSimAnnealParams))
                                                        then do
                                                            pure ([inPhyloGraph], inSimAnnealParams)
                                                        else do
                                                            -- otherwise do the anneal/Drift accept

                                                            (acceptGraph, nextSAParams) ← U.simAnnealAccept inSimAnnealParams (snd5 inPhyloGraph) . snd5 $ head newPhyloGraphList
                                                            case acceptGraph of
                                                                True → pure (newPhyloGraphList, nextSAParams)
                                                                _ → deleteNetEdgeRecursive inGS inData inPhyloGraph force nextSAParams (drop numGraphsToExamine inEdgeToDeleteList)


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
    → LG.Edge
    → PhyG ReducedPhylogeneticGraph
deleteNetEdge inGS inData inPhyloGraph force edgeToDelete =
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

                    -- (heuristicDelta, _, _) = heuristicDeleteDelta inGS inPhyloGraph edgeToDelete

                    -- edgeAddDelta = deltaPenaltyAdjustment inGS inPhyloGraph "delete"

                    -- full two-pass optimization--cycles checked in edge deletion function
                    let leafGraph = LG.extractLeafGraph $ thd5 inPhyloGraph

                    newPhyloGraph ←
                        if (graphType inGS == SoftWired)
                            then T.multiTraverseFullyLabelSoftWiredReduced inGS inData pruneEdges warnPruneEdges leafGraph startVertex delSimple
                            else
                                if (graphType inGS == HardWired)
                                    then T.multiTraverseFullyLabelHardWiredReduced inGS inData leafGraph startVertex delSimple
                                    else error "Unsupported graph type in deleteNetEdge.  Must be soft or hard wired"
                    -- check if deletino modified graph
                    if not wasModified
                        then do
                            pure inPhyloGraph
                        else -- else if force || (graphType inGS) == HardWired then

                            if force
                                then do
                                    -- trace ("DNE forced")
                                    pure newPhyloGraph
                                else -- if (heuristicDelta / (dynamicEpsilon inGS)) - edgeAddDelta < 0 then newPhyloGraph

                                    if (snd5 newPhyloGraph) < (snd5 inPhyloGraph)
                                        then do
                                            -- trace ("DNE Better: " <> (show $ snd5 newPhyloGraph))
                                            pure newPhyloGraph
                                        else do
                                            -- trace ("DNE Not Better: " <> (show $ snd5 newPhyloGraph))
                                            pure inPhyloGraph


-- )

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
                -- parentNodeToKeep = head $ filter (/= p1) parentsNodeToDelete
                -- newEdge = (parentNodeToKeep, head childrenNodeToDelete, 0.0)
                -- newGraph = LG.insEdge newEdge $ LG.delNode nodeToDelete inGraph
                newGraph = LG.delEdge inEdge inGraph
                -- newGraph' = GO.contractIn1Out1EdgesRename newGraph

                -- conversion as if input--see if affects length
                -- newGraph'' = GO.convertGeneralGraphToPhylogeneticGraph False newGraph
                newGraph'' = GO.contractIn1Out1EdgesRename newGraph
            in  -- error conditions and creation of chained network edges (forbidden in phylogenetic graph--causes resolutoin cache issues)
                if length childrenNodeToDelete /= 1
                    then error ("Cannot delete non-network edge in deleteNetworkEdge: (1)" <> (show inEdge) <> "\n" <> (LG.prettyIndices inGraph))
                    else
                        if length parentsNodeToDelete /= 2
                            then error ("Cannot delete non-network edge in deleteNetworkEdge (2): " <> (show inEdge) <> "\n" <> (LG.prettyIndices inGraph))
                            else -- warning if chained on input, skip if chained net edges in output

                                if (LG.isNetworkNode inGraph p1)
                                    then do
                                        -- error ("Error: Chained network nodes in deleteNetworkEdge : " <> (show inEdge) <> "\n" <> (LG.prettyIndices inGraph) <> " skipping")
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
                                                        {-trace ("DNE: Edge to delete " <> (show inEdge) <> " cnd " <> (show childrenNodeToDelete) <> " pnd " <> (show parentsNodeToDelete) <> " pntk " <> (show parentNodeToKeep)
                                                           <> " ne " <> (show newEdge) <> "\nInGraph: " <> (LG.prettyIndices inGraph) <> "\nNewGraph: " <> (LG.prettyIndices newGraph) <> "\nNewNewGraph: "
                                                           <> (LG.prettyIndices newGraph')) -}
                                                        pure (newGraph'', True)


{- | heuristicDeleteDelta takes the existing graph, edge to delete,
reoptimizes starting nodes of two created edges.  Returns cost delta based on
previous and new node resolution caches
delete n1 -> n2, create u -> v, u' -> v'
assumes original is edge n1 -> n2, u' -> (n2, X), n1 -> (n2,v), u (n1,Y)
-}
heuristicDeleteDelta
    ∷ GlobalSettings
    → ReducedPhylogeneticGraph
    → LG.Edge
    → PhyG (VertexCost, LG.LNode VertexInfo, LG.LNode VertexInfo)
heuristicDeleteDelta inGS inPhyloGraph (n1, n2) =
    if LG.isEmpty (fst5 inPhyloGraph)
        then error "Empty graph in heuristicDeleteDelta"
        else
            if graphType inGS == HardWired
                then -- ensures delete--will always be lower or equakl cost if delete edge from HardWired
                    pure (-1, dummyNode, dummyNode)
                else
                    let inGraph = thd5 inPhyloGraph
                        u = head $ LG.parents inGraph n1
                        u' = head $ filter (/= n1) $ LG.parents inGraph n2
                        v' = head $ LG.descendants inGraph n2
                        v = head $ filter (/= n2) $ LG.descendants inGraph n1

                        uLab = fromJust $ LG.lab inGraph u
                        uPrimeLab = fromJust $ LG.lab inGraph u'
                        vLab = fromJust $ LG.lab inGraph v
                        vPrimeLab = fromJust $ LG.lab inGraph v'

                        uOtherChild = head $ filter ((/= n1) . fst) $ LG.labDescendants inGraph (u, uLab)
                        uPrimeOtherChild = head $ filter ((/= n2) . fst) $ LG.labDescendants inGraph (u', uPrimeLab)
                    in  do
                            -- skip over netnodes
                            uLabAfter ← NEW.getOutDegree2VertexSoftWired inGS (fft5 inPhyloGraph) u (v, vLab) uOtherChild inGraph
                            uPrimeLabAfter ← NEW.getOutDegree2VertexSoftWired inGS (fft5 inPhyloGraph) u' (v', vPrimeLab) uPrimeOtherChild inGraph

                            -- cost of resolutions
                            let (_, uCostBefore) = NEW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uLab)
                            let (_, uPrimeCostBefore) = NEW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uPrimeLab)
                            let (_, uCostAfter) = NEW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uLabAfter)
                            let (_, uPrimeCostAfter) = NEW.extractDisplayTrees (Just (-1)) False (vertexResolutionData uPrimeLabAfter)

                            let addNetDelta = uCostAfter - uCostBefore + uPrimeCostAfter - uPrimeCostBefore

                            -- this should not happen--should try to crete new edges from children of net edges
                            if null (LG.parents inGraph n1)
                                || null (filter (/= n1) $ LG.parents inGraph n2)
                                || null (LG.descendants inGraph n2)
                                || null (filter (/= n2) $ LG.descendants inGraph n1)
                                || null (filter ((/= n2) . fst) $ LG.labDescendants inGraph (u', uPrimeLab))
                                || null (filter ((/= n1) . fst) $ LG.labDescendants inGraph (u, uLab))
                                then pure (infinity, dummyNode, dummyNode)
                                else -- this should not happen--should try to crete new edges from children of net edges

                                    if (length (LG.parents inGraph n1) /= 1)
                                        || (length (LG.parents inGraph n2) /= 2)
                                        || (length (LG.descendants inGraph n2) /= 1)
                                        || (length (LG.descendants inGraph n1) /= 2)
                                        then error ("Graph malformation in numbers of parents and children in heuristicDeleteDelta")
                                        else pure (addNetDelta, (u, uLabAfter), (u', uPrimeLabAfter))

{-
-- | insertNetEdgeBothDirections calls insertNetEdge for both u -> v and v -> u new edge orientations
insertNetEdgeBothDirections :: GlobalSettings -> ProcessedData -> ReducedPhylogeneticGraph ->  Maybe VertexCost -> (LG.LEdge b, LG.LEdge b) -> [ReducedPhylogeneticGraph]
insertNetEdgeBothDirections inGS inData inPhyloGraph preDeleteCost (u,v) = fmap (insertNetEdge inGS inData inPhyloGraph preDeleteCost) [(u,v), (v,u)]
-}
