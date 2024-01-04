{- |
Module exposing the functionality to swap sub-graphs of a phylogenetic graph.
-}
module Search.SwapMaster (
    swapMaster,
) where

import Commands.Verify qualified as VER
import Control.Monad (when)
import Data.Char
import Data.Foldable (fold)
import Data.Functor (($>))
import Data.Maybe
import GeneralUtilities
import Graphs.GraphOperations qualified as GO
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import Search.Swap qualified as S
import Text.Read
import Types.Types


{- | swapMaster processes and spawns the swap functions
the 2 x maxMoveDist since distance either side to list 2* dist on sorted edges
-}
swapMaster
    ∷ [Argument]
    → GlobalSettings
    → ProcessedData
    → [ReducedPhylogeneticGraph]
    → PhyG [ReducedPhylogeneticGraph]
swapMaster inArgs inGS inData [] = logWith LogInfo "No graphs to swap\n" $> []
swapMaster inArgs inGS inData inGraphListInput = do
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
        ) ←
        getSwapParams inArgs

    let swapType
            | any ((== "nni") . fst) lcArgList = NNI
            | any ((== "spr") . fst) lcArgList = SPR
            | any ((== "tbr") . fst) lcArgList = TBR
            | any ((== "alternate") . fst) lcArgList = Alternate
            | otherwise = Alternate

    let maxMoveEdgeDist = case swapType of
            NNI → 2
            _ → maxMoveEdgeDist'

    -- randomized orders of split and join-- not implemented
    -- doRandomized = any ((=="randomized").fst) lcArgList

    -- set implied alignment swapping
    let doIA'' = any ((== "ia") . fst) lcArgList

    -- steepest/all options
    let doSteepest' = any ((== "steepest") . fst) lcArgList
    let doAll = any ((== "all") . fst) lcArgList

    -- steepest default
    let doSteepest = ((not doSteepest' && not doAll) || doSteepest')

    -- simulated annealing parameters
    -- returnMutated to return annealed Graphs before swapping fir use in Genetic Algorithm
    let doAnnealing = any ((== "annealing") . fst) lcArgList

    let returnMutated = any ((== "returnmutated") . fst) lcArgList

    -- turn off union selection of rejoin--default to do both, union first
    let joinType
            | graphType inGS == HardWired = JoinAll
            | any ((== "joinall") . fst) lcArgList = JoinAll
            | any ((== "joinpruned") . fst) lcArgList = JoinPruned
            | any ((== "joinalternate") . fst) lcArgList = JoinAlternate
            | otherwise = JoinAlternate

    -- randomize split graph and rejoin edges, defualt to randomize
    let atRandom
            | any ((== "atrandom") . fst) lcArgList = True
            | any ((== "inOrder") . fst) lcArgList = False
            | otherwise = True

    -- populate SwapParams structure
    let localSwapParams =
            SwapParams
                { swapType = swapType
                , joinType = joinType
                , atRandom = atRandom
                , keepNum = keepNum
                , maxMoveEdgeDist = maxMoveEdgeDist
                , steepest = doSteepest
                , joinAlternate = False -- join prune alternates--turned off for now
                , doIA = doIA''
                , returnMutated = returnMutated
                }

    -- swap replicates is meant to allow multiple randomized swap trajectories
    -- set to 1 if not randomized swap or SA/Drifting (set by their own options)
    let replicates
            | not atRandom = 1
            | doAnnealing = 1
            | otherwise = replicateNumber

    -- replicate inGraphList based on 'replicates' for randomized trajectories
    let inGraphList = concat $ replicate replicates inGraphListInput
    let numGraphs = length inGraphList

    -- parallel setup
    let action ∷ [(Maybe SAParams, ReducedPhylogeneticGraph)] → PhyG ([ReducedPhylogeneticGraph], Int)
        action = {-# SCC swapMaster_action_swapSPRTBR #-} S.swapSPRTBR localSwapParams inGS inData 0 inGraphList
    simAnnealParams ←
        getSimAnnealParams doAnnealing doDrift steps' annealingRounds' driftRounds' acceptEqualProb acceptWorseFactor maxChanges

    -- create simulated annealing random lists uniquely for each fmap
    let newSimAnnealParamList = replicate numGraphs simAnnealParams

    let progressString
            | (not doAnnealing && not doDrift) =
                unwords
                    [ "Swapping"
                    , show $ length inGraphListInput
                    , "input graph(s) with"
                    , show replicates
                    , "trajectories at minimum cost"
                    , show . minimum $ snd5 <$> inGraphList
                    , "keeping maximum of"
                    , show keepNum
                    , "graphs per input graph\n"
                    ]
            | otherwise =
                case simAnnealParams of
                    Nothing → "Unknown procesing"
                    Just simAnneal → case method simAnneal of
                        SimAnneal →
                            unwords
                                [ "Simulated Annealing (Swapping)"
                                , show $ rounds simAnneal
                                , "rounds with"
                                , show $ numberSteps simAnneal
                                , "cooling steps"
                                , show $ length inGraphList
                                , "input graph(s) at minimum cost"
                                , show . minimum $ snd5 <$> inGraphList
                                , "keeping maximum of"
                                , show keepNum
                                , "graphs\n"
                                ]
                        _ →
                            unwords
                                [ "Drifting (Swapping)"
                                , show $ rounds simAnneal
                                , "rounds with"
                                , show $ driftMaxChanges simAnneal
                                , "maximum changes per round on"
                                , show $ length inGraphList
                                , "input graph(s) at minimum cost"
                                , show . minimum $ snd5 <$> inGraphList
                                , "keeping maximum of"
                                , show keepNum
                                , "graphs\n"
                                ]

    logWith LogInfo progressString
    -- TODO
    -- let graphPairList = PU.seqParMap (parStrategy $ strictParStrat inGS) (S.swapSPRTBR localSwapParams inGS inData 0 inGraphList) ((:[]) <$> zip3 (U.generateRandIntLists (head randomIntListSwap) numGraphs) newSimAnnealParamList inGraphList)

    let simAnnealList ∷ [[(Maybe SAParams, ReducedPhylogeneticGraph)]] -- [(Maybe SAParams,
        simAnnealList = fmap pure $ zip newSimAnnealParamList inGraphList
    pTraverse ← getParallelChunkTraverse
    graphPairList ← pTraverse action simAnnealList
    -- mapM (S.swapSPRTBR localSwapParams inGS inData 0 inGraphList) simAnnealList

    let (graphListList, counterList) = unzip graphPairList
    newGraphList ← GO.selectGraphs Best keepNum 0.0 $ concat graphListList
    let counter = sum counterList

    let finalGraphList
            | null newGraphList = inGraphList
            | otherwise = newGraphList

    let fullBuffWarning
            | length newGraphList < keepNum = ""
            | otherwise =
                "\n\tWarning--Swap returned as many minimum cost graphs as the 'keep' number.  \n\tThis may have limited the effectiveness of the swap. \n\tConsider increasing the 'keep' value or adding an additional swap."

    let endString
            | (not doAnnealing && not doDrift) =
                unwords
                    [ "\n\tAfter swap:"
                    , show $ length finalGraphList
                    , "resulting graphs with minimum cost"
                    , show . minimum $ snd5 <$> finalGraphList
                    , "with swap rounds (total):"
                    , show counter
                    , show swapType
                    ]
            | otherwise =
                case simAnnealParams of
                    Nothing → "Unknown procesing"
                    Just simAnneal → case method simAnneal of
                        SimAnneal →
                            unwords
                                [ "\n\tAfter Simulated Annealing:"
                                , show $ length finalGraphList
                                , "resulting graphs with minimum cost"
                                , show . minimum $ snd5 <$> finalGraphList
                                , "with swap rounds (total):"
                                , show counter
                                , show swapType
                                ]
                        _ →
                            unwords
                                [ "\n\tAfter Drifting:"
                                , show $ length finalGraphList
                                , "resulting graphs with minimum cost"
                                , show . minimum $ snd5 <$> finalGraphList
                                , "with swap rounds (total):"
                                , show counter
                                , show swapType
                                ]

    logWith LogInfo $ endString <> fullBuffWarning <> "\n"
    pure finalGraphList


-- | getSimumlatedAnnealingParams returns SA parameters
getSimAnnealParams
    ∷ Bool
    → Bool
    → Int
    → Int
    → Int
    → Double
    → Double
    → Int
    → PhyG (Maybe SAParams)
getSimAnnealParams doAnnealing doDrift steps annealingRounds' driftRounds' acceptEqualProb acceptWorseFactor maxChanges
    | not doAnnealing && not doDrift = pure Nothing
    | otherwise =
        let annealingRounds
                | annealingRounds' < 1 = 1
                | otherwise = annealingRounds'

            driftRounds
                | driftRounds' < 1 = 1
                | otherwise = driftRounds'

            saMethod
                | doDrift && doAnnealing = Drift
                | doDrift = Drift
                | otherwise = SimAnneal

            equalProb
                | acceptEqualProb < 0.0 = 0.0
                | acceptEqualProb > 1.0 = 1.0
                | otherwise = acceptEqualProb

            worseFactor = max 0 acceptWorseFactor

            changes
                | maxChanges >= 0 = maxChanges
                | otherwise = 15

            saValues =
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
                pure saValues


-- | getSwapParams takes arg list and preocesses returning parameter values
getSwapParams
    ∷ [Argument]
    → PhyG
        ( Int
        , Int
        , Int
        , Int
        , Bool
        , Int
        , Double
        , Double
        , Int
        , Int
        , [(String, String)]
        )
getSwapParams inArgs =
    let fstArgList = fmap (fmap toLower . fst) inArgs
        sndArgList = fmap (fmap toLower . snd) inArgs
        lcArgList = zip fstArgList sndArgList
        checkCommandList = checkCommandArgs "swap" fstArgList VER.swapArgList

        gatherKeys ∷ [String] → [String]
        gatherKeys keys = snd <$> filter ((`elem` keys) . fst) lcArgList

        tryExtractValueKeyedAs ∷ (Read a) ⇒ Bool → a → String → [String] → PhyG a
        tryExtractValueKeyedAs isReal def label keys = case gatherKeys keys of
            [] → pure def
            str : [] → case readMaybe str of
                Nothing →
                    failWithPhase Parsing $
                        fold
                            [ "'"
                            , label
                            , "' specification not "
                            , if isReal then "a real value" else "an integer value"
                            , " in swap:\n\t"
                            , show str
                            ]
                Just val → pure val
            _ : _ →
                failWithPhase Parsing $
                    fold
                        [ "Multiple '"
                        , label
                        , "'"
                        , "specifications in 'swap' command -- only one is permitted:\n\t"
                        , show inArgs
                        ]

        limit ∷ Int
        limit = maxBound `div` 3
    in  do
            -- check for valid command options
            when (not checkCommandList) $ failWithPhase Parsing $ "Unrecognized command in 'swap': " <> show inArgs

            (,,,,,,,,,,)
                <$> tryExtractValueKeyedAs False 10 "keep" ["keep"]
                <*> tryExtractValueKeyedAs False limit "maximum edge distances" ["alternate", "spr", "tbr", "nni"]
                <*> tryExtractValueKeyedAs False 10 "annealing steps" ["steps"]
                <*> tryExtractValueKeyedAs False 10 "annealing" ["annealing"]
                <*> pure (not . null $ gatherKeys ["drift"])
                <*> tryExtractValueKeyedAs False 1 "drift rounds" ["drift"]
                <*> tryExtractValueKeyedAs True 0.5 "drift acceptEqual" ["acceptequal"]
                <*> tryExtractValueKeyedAs True 20 "drift acceptWorse" ["acceptWorse"]
                <*> tryExtractValueKeyedAs False 15 "drift maxChanges" ["maxchanges"]
                <*> tryExtractValueKeyedAs False 1 "swap replicates" ["replicates"]
                <*> pure lcArgList
