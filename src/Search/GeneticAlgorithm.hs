{- |
Module specifying graph sGeneticAlgorithm functions
-}
module Search.GeneticAlgorithm (
    geneticAlgorithm,
) where

import Control.Monad (when)
import Control.Monad.IO.Class (MonadIO (..))
import Control.Monad.Random.Class
import Data.Foldable (fold)
import GeneralUtilities
import Graphs.GraphOperations qualified as GO
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import PHANE.Evaluation.Verbosity (Verbosity (..))
import Search.Fuse qualified as F
import Search.NetworkAddDelete qualified as N
import Search.Swap qualified as S
import Types.Types
import Utilities.LocalGraph qualified as LG


{- | geneticAlgorithm takes arguments and performs genetic algorithm on input graphs
the process follows several steps
1) input graphs are mutated
      this step is uncharacteristically first so that is can operate on
      graphs that have been "fused" (recombined) already
      mutated graphs are added up to popsize
      if input graphs are already at the population size, an equal number of mutants are added (exceeding input popsize)
2) graph are recombined using fusing operations
3) population undergoes selection to population size (unique graphs)
      selection based on delta with best graph and severity factor on (0,Inf) 1 pure cost delta < 1 more severe, > 1 less severe
      if "elitist" (default) 'best' graphs are always selected to ensure no worse.
4) operation repearts for number of generations
-}
geneticAlgorithm
    ∷ GlobalSettings
    → ProcessedData
    → Bool
    → Int
    → Int
    → Int
    → Int
    → Int
    → Double
    → Int
    → Int
    → Int
    → [ReducedPhylogeneticGraph]
    → PhyG ([ReducedPhylogeneticGraph], Int)
geneticAlgorithm inGS inData doElitist maxNetEdges keepNum popSize generations generationCounter severity recombinations stopCount stopNum inGraphList =
    if null inGraphList
        then return ([], 0)
        else
            if generationCounter == generations
                then return (inGraphList, generationCounter)
                else
                    if stopCount >= stopNum
                        then return (inGraphList, generationCounter)
                        else
                            let -- get elite list of best solutions
                                initialEliteList = GO.selectGraphs Best (maxBound ∷ Int) 0.0 (-1) inGraphList
                            in  do
                                    logWith LogInfo ("Genetic algorithm generation: " <> (show generationCounter) <> "\n")

                                    seedList ← getRandoms
                                    -- mutate input graphs, produces number input, limited to popsize
                                    mutatedGraphList' ←
                                        traverse (mutateGraph inGS inData maxNetEdges) $ (takeRandom (seedList !! 1) popSize inGraphList)

                                    let numShort = popSize - (length mutatedGraphList')
                                    let randList = (randomIntList $ seedList !! 2)
                                    let graphList = takeRandom (seedList !! 3) numShort inGraphList

                                    additionalMutated ← traverse (mutateGraph inGS inData maxNetEdges) graphList
                                    -- adjust to correct populationsize if input number < popSize
                                    let mutatedGraphList =
                                            if length mutatedGraphList' >= popSize
                                                then mutatedGraphList'
                                                else mutatedGraphList' <> additionalMutated

                                    -- get unique graphs, no point in recombining repetitions
                                    let uniqueMutatedGraphList = GO.selectGraphs Unique (maxBound ∷ Int) 0.0 (-1) (mutatedGraphList <> inGraphList)

                                    -- recombine elite with mutated and mutated with mutated
                                    let recombineSwap = getRandomElement (seedList !! 4) [NoSwap, NNI, SPR] --  these take too long, "tbr", "alternate"]

                                    -- options to join via union choices or all in fuse
                                    -- this is ignored for now in fuse--JoinAll is what it does
                                    let joinType = getRandomElement (seedList !! 6) [JoinAlternate, JoinAll]

                                    let doSteepest = True
                                    let returnBest = False
                                    let returnUnique = True
                                    let singleRound = False
                                    let fusePairs = Just recombinations
                                    let randomPairs = True
                                    let reciprocal = False

                                    -- populate SwapParams structure
                                    let swapParams =
                                            SwapParams
                                                { swapType = recombineSwap
                                                , joinType = joinType
                                                , atRandom = True -- randomize swap order
                                                , keepNum = (2 * popSize)
                                                , maxMoveEdgeDist = (maxBound ∷ Int)
                                                , steepest = doSteepest
                                                , joinAlternate = False -- not working now
                                                , doIA = False
                                                , returnMutated = False
                                                }

                                    (recombinedGraphList, _) ←
                                        F.fuseAllGraphs
                                            swapParams
                                            inGS
                                            inData
                                            0
                                            returnBest
                                            returnUnique
                                            singleRound
                                            fusePairs
                                            randomPairs
                                            reciprocal
                                            uniqueMutatedGraphList

                                    -- selection of graphs population
                                    -- unique sorted on cost so getting unique with lowest cost
                                    let selectedGraphs = GO.selectGraphs Unique popSize 0.0 (-1) recombinedGraphList
                                    let newCost = snd5 $ head selectedGraphs

                                    -- if new graphs better cost then take those
                                    if newCost < (snd5 $ head initialEliteList)
                                        then
                                            geneticAlgorithm
                                                inGS
                                                inData
                                                doElitist
                                                maxNetEdges
                                                keepNum
                                                popSize
                                                generations
                                                (generationCounter + 1)
                                                severity
                                                recombinations
                                                0
                                                stopNum
                                                selectedGraphs
                                        else -- if new graphs not better then add in elites to ensure monotonic decrease in cost

                                            let newGraphList = GO.selectGraphs Unique keepNum 0.0 (-1) (initialEliteList <> selectedGraphs)
                                            in  geneticAlgorithm
                                                    inGS
                                                    inData
                                                    doElitist
                                                    maxNetEdges
                                                    keepNum
                                                    popSize
                                                    generations
                                                    (generationCounter + 1)
                                                    severity
                                                    recombinations
                                                    (stopCount + 1)
                                                    stopNum
                                                    newGraphList


-- | mutateGraph mutates a graph using drift functionality
mutateGraph ∷ GlobalSettings → ProcessedData → Int → ReducedPhylogeneticGraph → PhyG ReducedPhylogeneticGraph
mutateGraph inGS inData maxNetEdges inGraph
    | LG.isEmpty (fst5 inGraph) = failWithPhase Computing "Empty graph in mutateGraph"
    | otherwise =
        let getRandomFrom es = (`getRandomElement` es) <$> getRandom

            -- static parameter values
            valAlternate = False
            valDoIA = False
            valDoRandomOrder = True
            valJoinType = JoinAll -- keep selection of rejoins based on all possibilities
            valMaxMoveEdgeDist = 10000
            valNumToKeep = 5
            valReturnMutated = True
            valSteepest = True

            -- randomize simulated annealing parameters
            getRandomSAParams = do
                rDrift ← getRandomFrom [5, 10, 20]
                rSteps ← getRandomFrom [5, 10, 20]
                rMethod ← getRandomFrom [Drift, SimAnneal]
                rStream ← getRandoms
                pure . Just $
                    SAParams
                        { method = rMethod
                        , numberSteps = rSteps
                        , currentStep = 0
                        , randomIntegerList = rStream
                        , rounds = 1
                        , driftAcceptEqual = 0.5
                        , driftAcceptWorse = 2.0
                        , -- this could be an important factor don't want too severe, but significant
                          driftMaxChanges = rDrift
                        , driftChanges = 0
                        }

            -- randomize 'swap' parameters
            getRandomSwapParams = do
                swapType ← getRandomFrom [SPR, Alternate]
                -- randomize network edit parameters
                let doRandomOrder = True
                -- populate SwapParams structure
                pure $
                    SwapParams
                        { swapType = swapType
                        , joinType = valJoinType
                        , atRandom = True -- randomize split and rejoin edge orders
                        , keepNum = valNumToKeep
                        , maxMoveEdgeDist = valMaxMoveEdgeDist
                        , steepest = valSteepest
                        , joinAlternate = False
                        , doIA = valDoIA
                        , returnMutated = valReturnMutated
                        }

            firstOrOldIfNoneExists =
                pure . \case
                    ([], _) → inGraph
                    (x : _, _) → x

            mutateOption1 = do
                rSAParams ← getRandomSAParams
                rStream ← getRandoms
                rSwapParams ← getRandomSwapParams
                firstOrOldIfNoneExists =<< S.swapSPRTBR rSwapParams inGS inData 0 [inGraph] [(rStream, rSAParams, inGraph)]

            mutateOption2 = do
                rSAParams ← getRandomSAParams
                rStream ← getRandoms
                rSwapParams ← getRandomSwapParams
                firstOrOldIfNoneExists =<< S.swapSPRTBR rSwapParams inGS inData 0 [inGraph] [(rStream, rSAParams, inGraph)]

            mutateOption3 = do
                rSAParams ← getRandomSAParams
                rSwapParams ← getRandomSwapParams
                rVal ← getRandom
                firstOrOldIfNoneExists
                    =<< N.moveAllNetEdges
                        inGS
                        inData
                        rVal
                        maxNetEdges
                        valNumToKeep
                        0
                        valReturnMutated
                        valSteepest
                        valDoRandomOrder
                        ([], infinity)
                        (rSAParams, [inGraph])

            mutateOption4 = do
                rSAParams ← getRandomSAParams
                rSwapParams ← getRandomSwapParams
                rVal ← getRandom
                firstOrOldIfNoneExists
                    =<< N.moveAllNetEdges
                        inGS
                        inData
                        rVal
                        maxNetEdges
                        valNumToKeep
                        0
                        valReturnMutated
                        valSteepest
                        valDoRandomOrder
                        ([], infinity)
                        (rSAParams, [inGraph])

            mutateOption5 = do
                rMaxRounds ← getRandomFrom [1 .. 5]
                rSAParams ← getRandomSAParams
                rSwapParams ← getRandomSwapParams
                rVal ← getRandom
                firstOrOldIfNoneExists
                    =<< N.insertAllNetEdges
                        inGS
                        inData
                        rVal
                        maxNetEdges
                        valNumToKeep
                        rMaxRounds
                        0
                        valReturnMutated
                        valSteepest
                        valDoRandomOrder
                        ([], infinity)
                        (rSAParams, [inGraph])

            mutateOption6 = do
                rMaxRounds ← getRandomFrom [1 .. 5]
                rSAParams ← getRandomSAParams
                rSwapParams ← getRandomSwapParams
                rVal ← getRandom
                firstOrOldIfNoneExists
                    =<< N.addDeleteNetEdges
                        inGS
                        inData
                        rVal
                        maxNetEdges
                        valNumToKeep
                        rMaxRounds
                        0
                        valReturnMutated
                        valSteepest
                        valDoRandomOrder
                        ([], infinity)
                        (rSAParams, [inGraph])

            mutateOption7 = do
                rSAParams ← getRandomSAParams
                rSwapParams ← getRandomSwapParams
                rVal ← getRandom
                firstOrOldIfNoneExists
                    =<< N.deleteAllNetEdges
                        inGS
                        inData
                        rVal
                        maxNetEdges
                        valNumToKeep
                        0
                        valReturnMutated
                        valSteepest
                        valDoRandomOrder
                        ([], infinity)
                        (rSAParams, [inGraph])
        in  do
                -- randomize edit type
                editType ← getRandomFrom ["swap", "netEdge"]
                netEditType ← getRandomFrom ["netAdd", "netDelete", "netAddDelete"] -- , "netMove"]
                case (graphType inGS, editType, netEditType) of
                    -- only swap mutation stuff for tree
                    (Tree, _, nEdit) | nEdit `notElem` ["netAdd", "netDelete", "netAddDelete", "netMove"] → mutateOption1
                    -- graphs choose what type of mutation at random
                    (_, "swap", _) → mutateOption2
                    -- move only for Hardwired
                    (HardWired, _, _) → mutateOption3
                    -- SoftWired
                    (SoftWired, _, "netMove") → mutateOption4
                    (SoftWired, _, "netAdd") → mutateOption5
                    (SoftWired, _, "netAddDelete") → mutateOption6
                    (SoftWired, _, "netDelete") → mutateOption7
                    (SoftWired, _, val) →
                        failWithPhase Parsing $
                            fold
                                ["Unrecognized edit type '", val, "' for sofwired network"]
