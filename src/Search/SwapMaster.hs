{- |
Module exposing the functionality to swap sub-graphs of a phylogenetic graph.
-}
module Search.SwapMaster (
    swapMaster,
) where

import Commands.Verify qualified as VER
import Control.Monad (when)
import Control.Monad.IO.Class (MonadIO (..))
import Control.Monad.Random.Class
import Data.Bifunctor (first)
import Data.Char
import Data.Foldable (fold)
import Data.Functor ((<&>))
import Data.Maybe
import GeneralUtilities
import Graphs.GraphOperations qualified as GO
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import PHANE.Evaluation.Verbosity (Verbosity (..))
import Search.Swap qualified as S
import Text.Read
import Types.Types
import Utilities.Utilities as U


{- | swapMaster processes and spawns the swap functions
the 2 x maxMoveDist since distance either side to list 2* dist on sorted edges
-}
swapMaster
    ∷ [Argument]
    → GlobalSettings
    → ProcessedData
    → [ReducedPhylogeneticGraph]
    → PhyG [ReducedPhylogeneticGraph]
swapMaster inArgs inGS inData inGraphListInput =
    {-# SCC swapMaster_TOP_DEF #-}
    if null inGraphListInput
        then do
            logWith LogInfo "No graphs to swap\n"
            pure []
        else -- else if graphType inGS == HardWired then trace ("Swapping hardwired graphs is currenty not implemented") inGraphList

            let -- process args for swap
                ( keepNum
                    , maxMoveEdgeDist'
                    , steps'
                    , annealingRounds'
                    , doDrift
                    , driftRounds'
                    , acceptEqualProb
                    , acceptWorseFactor
                    , maxChanges
                    , replicateNumber
                    , lcArgList
                    ) = getSwapParams inArgs

                swapType
                    | any ((== "nni") . fst) lcArgList = NNI
                    | any ((== "spr") . fst) lcArgList = SPR
                    | any ((== "tbr") . fst) lcArgList = TBR
                    | any ((== "alternate") . fst) lcArgList = Alternate
                    | otherwise = Alternate

                maxMoveEdgeDist =
                    if swapType == NNI
                        then 2
                        else fromJust maxMoveEdgeDist'

                -- randomized orders of split and join-- not implemented
                -- doRandomized = any ((=="randomized").fst) lcArgList

                -- set implied alignment swapping
                doIA' = any ((== "ia") . fst) lcArgList
                doIA'' = doIA'

                --- steepest/all options
                doSteepest' = any ((== "steepest") . fst) lcArgList
                doAll = any ((== "all") . fst) lcArgList

                -- steepest default
                doSteepest = ((not doSteepest' && not doAll) || doSteepest')

                -- simulated annealing parameters
                -- returnMutated to return annealed Graphs before swapping fir use in Genetic Algorithm
                doAnnealing = any ((== "annealing") . fst) lcArgList

                returnMutated = any ((== "returnmutated") . fst) lcArgList

                -- checking of heuristic graph costs
                heuristicCheck 
                    | any ((== "bestonly") . fst) lcArgList = BestOnly
                    | any ((== "better") . fst) lcArgList = Better
                    | any ((== "bettern") . fst) lcArgList = BetterN
                    | any ((== "bestall") . fst) lcArgList = BestAll
                    | otherwise = BetterN

                -- turn off union selection of rejoin--default to do both, union first
                joinType
                    | graphType inGS == HardWired = JoinAll
                    | any ((== "joinall") . fst) lcArgList = JoinAll
                    | any ((== "joinpruned") . fst) lcArgList = JoinPruned
                    | any ((== "joinalternate") . fst) lcArgList = JoinAlternate
                    | otherwise = JoinAll

                -- randomize split graph and rejoin edges, defualt to randomize
                atRandom
                    | any ((== "atrandom") . fst) lcArgList = True
                    | any ((== "inOrder") . fst) lcArgList = False
                    | swapType == NNI = False
                    | otherwise = True

                -- split edge order based on greartest diffenrece in costr when graph is split
                    -- does all of them before sorting
                    -- since comes after testing for random will override
                sortEdgesSplitCost
                    | any ((== "splitsequential") . fst) lcArgList = False
                    | any ((== "sortsplit") . fst) lcArgList = True
                    | atRandom = False
                    | otherwise = True

                -- when plitting base graph--do in parallel or via recursive sequential
                -- might save on memeory, coulod be a bit more efficient time-wise
                -- definately affects trajectory--small examples had worse optimality outcomes
                parallelSplit
                    | sortEdgesSplitCost = True
                    | any ((== "splitparallel") . fst) lcArgList = True
                    | any ((== "splitsequential") . fst) lcArgList = False
                    | otherwise = True


                -- populate SwapParams structure
                localSwapParams = 
                    SwapParams
                        { swapType = swapType
                        , joinType = joinType
                        , atRandom = atRandom
                        , checkHeuristic = heuristicCheck
                        , sortEdgesSplitCost = sortEdgesSplitCost
                        , keepNum = (fromJust keepNum)
                        , maxMoveEdgeDist = maxMoveEdgeDist
                        , splitParallel = parallelSplit
                        , steepest = doSteepest
                        , joinAlternate = False -- join prune alternates--turned off for now
                        , doIA = doIA''
                        , returnMutated = returnMutated
                        }

                -- swap replicates is meant to allow multiple randomized swap trajectories
                -- set to 1 if not randomized swap or SA/Drifting (set by their own options)
                replicates
                    | not atRandom = 1
                    | doAnnealing = 1
                    | otherwise = fromJust replicateNumber

                -- replicate inGraphList based on 'replicates' for randomized trajectories
                inGraphList = concat $ replicate replicates inGraphListInput
                numGraphs = length inGraphList

                -- parallel setup
                action ∷ [(Maybe SAParams, ReducedPhylogeneticGraph)] → PhyG ([ReducedPhylogeneticGraph], Int)
                -- action = {-# SCC swapMaster_action_swapSPRTBR #-} S.swapSPRTBR localSwapParams inGS inData 0 inGraphList
                action = {-# SCC swapMaster_action_swapSPRTBR #-} S.swapDriver localSwapParams inGS inData 0 inGraphList
            in  do
                    simAnnealParams ←
                        getSimAnnealParams doAnnealing doDrift steps' annealingRounds' driftRounds' acceptEqualProb acceptWorseFactor maxChanges

                    -- create simulated annealing random lists uniquely for each fmap
                    let newSimAnnealParamList = replicate numGraphs simAnnealParams

                    let progressString
                            | (not doAnnealing && not doDrift) =
                                ( "Swapping "
                                    <> show (length inGraphListInput)
                                    <> " input graph(s) with "
                                    <> show replicates
                                    <> " trajectories at minimum cost "
                                    <> show (minimum $ fmap snd5 inGraphList)
                                    <> " keeping maximum of "
                                    <> show (fromJust keepNum)
                                    <> " graphs per input graph"
                                    <> "\n"
                                )
                            | method (fromJust simAnnealParams) == SimAnneal =
                                ( "Simulated Annealing (Swapping) "
                                    <> show (rounds $ fromJust simAnnealParams)
                                    <> " rounds with "
                                    <> show (numberSteps $ fromJust simAnnealParams)
                                    <> " cooling steps "
                                    <> show (length inGraphList)
                                    <> " input graph(s) at minimum cost "
                                    <> show (minimum $ fmap snd5 inGraphList)
                                    <> " keeping maximum of "
                                    <> show (fromJust keepNum)
                                    <> " graphs"
                                    <> "\n"
                                )
                            | otherwise =
                                "Drifting (Swapping) "
                                    <> show (rounds $ fromJust simAnnealParams)
                                    <> " rounds with "
                                    <> show (driftMaxChanges $ fromJust simAnnealParams)
                                    <> " maximum changes per round on "
                                    <> show (length inGraphList)
                                    <> " input graph(s) at minimum cost "
                                    <> show (minimum $ fmap snd5 inGraphList)
                                    <> " keeping maximum of "
                                    <> show (fromJust keepNum)
                                    <> " graphs"
                                    <> "\n"

                    logWith LogInfo progressString

                    let simAnnealList = (: []) <$> zip newSimAnnealParamList inGraphList
                    graphPairList ←
                        getParallelChunkTraverse >>= \pTraverse →
                            action `pTraverse` simAnnealList

                    let (graphListList, counterList) = first fold $ unzip graphPairList
                    (newGraphList, counter) ← GO.selectGraphs Best (fromJust keepNum) 0 graphListList <&> \x → (x, sum counterList)

                    let finalGraphList = case newGraphList of
                            [] → inGraphList
                            _ → newGraphList

                    let fullBuffWarning =
                            if length newGraphList >= (fromJust keepNum)
                                then
                                    "\n\tWarning--Swap returned as many minimum cost graphs as the 'keep' number.  \n\tThis may have limited the effectiveness of the swap. \n\tConsider increasing the 'keep' value or adding an additional swap."
                                else ""

                    let endString
                            | (not doAnnealing && not doDrift) =
                                ( "\n\tAfter swap: "
                                    <> show (length finalGraphList)
                                    <> " resulting graphs with minimum cost "
                                    <> show (minimum $ fmap snd5 finalGraphList)
                                    <> " with swap rounds (total): "
                                    <> show counter
                                    <> " "
                                    <> show swapType
                                )
                            | method (fromJust simAnnealParams) == SimAnneal =
                                ( "\n\tAfter Simulated Annealing: "
                                    <> show (length finalGraphList)
                                    <> " resulting graphs with minimum cost "
                                    <> show (minimum $ fmap snd5 finalGraphList)
                                    <> " with swap rounds (total): "
                                    <> show counter
                                    <> " "
                                    <> show swapType
                                )
                            | otherwise =
                                "\n\tAfter Drifting: "
                                    <> show (length finalGraphList)
                                    <> " resulting graphs with minimum cost "
                                    <> show (minimum $ fmap snd5 finalGraphList)
                                    <> " with swap rounds (total): "
                                    <> show counter
                                    <> " "
                                    <> show swapType

                    logWith LogInfo (endString <> fullBuffWarning <> "\n")
                    pure finalGraphList


-- | getSimumlatedAnnealingParams returns SA parameters
-- set SA?Drif max changes / number steps > 0 so can always check if one otr other matches 
-- to terminate in swapping etc
getSimAnnealParams
    ∷ Bool
    → Bool
    → Maybe Int
    → Maybe Int
    → Maybe Int
    → Maybe Double
    → Maybe Double
    → Maybe Int
    → PhyG (Maybe SAParams)
getSimAnnealParams doAnnealing doDrift steps' annealingRounds' driftRounds' acceptEqualProb acceptWorseFactor maxChanges
    | not doAnnealing && not doDrift = pure Nothing
    | otherwise =
        let steps = max 3 (fromJust steps')

            annealingRounds = case annealingRounds' of
                Just v | 1 <= v → v
                _ → 1

            driftRounds = case driftRounds' of
                Just v | 1 <= v → v
                _ → 1

            saMethod
                | doDrift && doAnnealing = Drift
                | doDrift = Drift
                | otherwise = SimAnneal

            equalProb = case acceptEqualProb of
                Nothing → 0
                Just v | v > 1 → 1
                Just v | v < 0 → 0
                Just v → v

            worseFactor = max (fromJust acceptWorseFactor) 0.0

            changes = case maxChanges of
                Just num | num > 0 → num
                _ → 15

            getResult =
                Just $
                    SAParams
                        { method = saMethod
                        , numberSteps = steps
                        , currentStep = 0
                        , rounds = max annealingRounds driftRounds
                        , driftAcceptEqual = equalProb
                        , driftAcceptWorse = worseFactor
                        , driftMaxChanges = changes
                        , driftChanges = 0
                        }
        in  do
                when (doDrift && doAnnealing) $
                    logWith
                        LogWarn
                        "\tSpecified both Simulated Annealing (with temperature steps) and Drifting (without)--defaulting to drifting.\n"
                pure $ getResult


-- | getSwapParams takes areg list and preocesses returning parameter values
getSwapParams
    ∷ [Argument]
    → ( Maybe Int
      , Maybe Int
      , Maybe Int
      , Maybe Int
      , Bool
      , Maybe Int
      , Maybe Double
      , Maybe Double
      , Maybe Int
      , Maybe Int
      , [(String, String)]
      )
getSwapParams inArgs =
    let fstArgList = fmap (fmap toLower . fst) inArgs
        sndArgList = fmap (fmap toLower . snd) inArgs
        lcArgList = zip fstArgList sndArgList
        checkCommandList = checkCommandArgs "swap" fstArgList VER.swapArgList
    in  -- check for valid command options
        if not checkCommandList
            then errorWithoutStackTrace ("Unrecognized command in 'swap': " <> show inArgs)
            else
                let keepList = filter ((== "keep") . fst) lcArgList
                    keepNum
                        | length keepList > 1 =
                            errorWithoutStackTrace ("Multiple 'keep' number specifications in swap command--can have only one: " <> show inArgs)
                        | null keepList = Just 10
                        | otherwise = readMaybe (snd $ head keepList) ∷ Maybe Int

                    moveLimitList = filter (not . null) (snd <$> filter ((`elem` ["alternate", "spr", "tbr", "nni"]) . fst) lcArgList)
                    maxMoveEdgeDist'
                        | length moveLimitList > 1 =
                            errorWithoutStackTrace
                                ("Multiple maximum edge distance number specifications in swap command--can have only one (e.g. spr:2): " <> show inArgs)
                        | null moveLimitList = Just ((maxBound ∷ Int) `div` 3)
                        | otherwise = readMaybe (head moveLimitList) ∷ Maybe Int

                    -- simulated anealing options
                    stepsList = filter ((== "steps") . fst) lcArgList
                    steps'
                        | length stepsList > 1 =
                            errorWithoutStackTrace
                                ("Multiple annealing steps value specifications in swap command--can have only one (e.g. steps:10): " <> show inArgs)
                        | null stepsList = Just 10
                        | otherwise = readMaybe (snd $ head stepsList) ∷ Maybe Int

                    annealingList = filter ((== "annealing") . fst) lcArgList
                    annealingRounds'
                        | length annealingList > 1 =
                            errorWithoutStackTrace ("Multiple 'annealing' rounds number specifications in swap command--can have only one: " <> show inArgs)
                        | null annealingList = Just 1
                        | otherwise = readMaybe (snd $ head annealingList) ∷ Maybe Int

                    -- drift options
                    doDrift = any ((== "drift") . fst) lcArgList

                    driftList = filter ((== "drift") . fst) lcArgList
                    driftRounds'
                        | length driftList > 1 =
                            errorWithoutStackTrace ("Multiple 'drift' rounds number specifications in swap command--can have only one: " <> show inArgs)
                        | null driftList = Just 1
                        | otherwise = readMaybe (snd $ head driftList) ∷ Maybe Int

                    acceptEqualList = filter ((== "acceptequal") . fst) lcArgList
                    acceptEqualProb
                        | length acceptEqualList > 1 =
                            errorWithoutStackTrace ("Multiple 'drift' acceptEqual specifications in swap command--can have only one: " <> show inArgs)
                        | null acceptEqualList = Just 0.5
                        | otherwise = readMaybe (snd $ head acceptEqualList) ∷ Maybe Double

                    acceptWorseList = filter ((== "acceptworse") . fst) lcArgList
                    acceptWorseFactor
                        | length acceptWorseList > 1 =
                            errorWithoutStackTrace ("Multiple 'drift' acceptWorse specifications in swap command--can have only one: " <> show inArgs)
                        | null acceptWorseList = Just 20.0
                        | otherwise = readMaybe (snd $ head acceptWorseList) ∷ Maybe Double

                    maxChangesList = filter ((== "maxchanges") . fst) lcArgList
                    maxChanges
                        | length maxChangesList > 1 =
                            errorWithoutStackTrace ("Multiple 'drift' maxChanges number specifications in swap command--can have only one: " <> show inArgs)
                        | null maxChangesList = Just 15
                        | otherwise = readMaybe (snd $ head maxChangesList) ∷ Maybe Int

                    replicatesList = filter ((== "replicates") . fst) lcArgList
                    replicates
                        | length replicatesList > 1 =
                            errorWithoutStackTrace ("Multiple 'swap' replicates number specifications in swap command--can have only one: " <> show inArgs)
                        | null replicatesList = Just 1
                        | otherwise = readMaybe (snd $ head replicatesList) ∷ Maybe Int
                in  -- check inputs
                    if isNothing keepNum
                        then errorWithoutStackTrace ("Keep specification not an integer in swap: " <> show (head keepList))
                        else
                            if isNothing maxMoveEdgeDist'
                                then errorWithoutStackTrace ("Maximum edge move distance specification not an integer (e.g. spr:2): " <> show (head moveLimitList))
                                else
                                    if isNothing steps'
                                        then errorWithoutStackTrace ("Annealing steps specification not an integer (e.g. steps:10): " <> show (snd $ head stepsList))
                                        else
                                            if isNothing acceptEqualProb
                                                then
                                                    errorWithoutStackTrace
                                                        ("Drift 'acceptEqual' specification not a float (e.g. acceptEqual:0.75): " <> show (snd $ head acceptEqualList))
                                                else
                                                    if isNothing acceptWorseFactor
                                                        then
                                                            errorWithoutStackTrace
                                                                ("Drift 'acceptWorse' specification not a float (e.g. acceptWorse:1.0): " <> show (snd $ head acceptWorseList))
                                                        else
                                                            if isNothing maxChanges
                                                                then
                                                                    errorWithoutStackTrace
                                                                        ("Drift 'maxChanges' specification not an integer (e.g. maxChanges:10): " <> show (snd $ head maxChangesList))
                                                                else
                                                                    if isNothing replicates
                                                                        then
                                                                            errorWithoutStackTrace
                                                                                ("Swap 'replicates' specification not an integer (e.g. replicates:5): " <> show (snd $ head replicatesList))
                                                                        else -- trace ("GSP: " <> (show inArgs) <> " " <> (show )(keepNum, maxMoveEdgeDist', steps', annealingRounds', doDrift, driftRounds', acceptEqualProb, acceptWorseFactor, maxChanges, lcArgList))

                                                                            ( keepNum
                                                                            , maxMoveEdgeDist'
                                                                            , steps'
                                                                            , annealingRounds'
                                                                            , doDrift
                                                                            , driftRounds'
                                                                            , acceptEqualProb
                                                                            , acceptWorseFactor
                                                                            , maxChanges
                                                                            , replicates
                                                                            , lcArgList
                                                                            )
