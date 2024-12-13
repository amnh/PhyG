{- |
Module exposing the functionality to swap sub-graphs of a phylogenetic graph.
-}
module Search.SwapMaster (
    swapMaster,
) where

import Commands.Verify qualified as VER
import Control.Monad (when)
import Data.Bifunctor (first)
import Data.Char
import Data.Foldable (fold)
import Data.Functor ((<&>))
import Data.Maybe
import GeneralUtilities
import Graphs.GraphOperations qualified as GO
import GraphOptimization.Traversals qualified as T
import PHANE.Evaluation
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import Search.SwapV2 qualified as SV2
import Text.Read
import Types.Types
-- import Debug.Trace

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
                inputBestCost = minimum $ fmap snd5 inGraphList
                ( keepNum
                    , maxMoveEdgeDist'
                    , steps'
                    , doAnnealing
                    , annealingRounds'
                    , doDrift
                    , driftRounds'
                    , acceptEqualProb
                    , acceptWorseFactor
                    , maxChanges
                    , replicateNumber
                    , levelNumber
                    , lcArgList
                    ) = getSwapParams inArgs

                -- local multiTraverse control option
                -- Default MultiTraverse global setting--need to rediagnose if set differnet from swap or global option
                multiTraverseValue = filter ((== "multitraverse") . fst) lcArgList
                doMultiTraverse'  
                    | length multiTraverseValue > 1 =
                                errorWithoutStackTrace ("Multiple multiTraverse specifications in swap--can have only one: " <> show inArgs)
                    | null multiTraverseValue = Just $ fmap toLower $ show (multiTraverseCharacters inGS)
                    | null (snd $ head multiTraverseValue) = errorWithoutStackTrace ("1-MultiTraverse swap option must be 'True' or 'False'" <> show inArgs)
                    | otherwise = readMaybe (show $ snd $ head multiTraverseValue) ∷ Maybe String
                
                doMultiTraverse = if isNothing doMultiTraverse' then errorWithoutStackTrace ("2-MultiTraverse swap option must be 'True' or 'False'" <> show inArgs)
                                  else if fromJust doMultiTraverse'  == "true" then True
                                  else if fromJust doMultiTraverse'  == "false" then False
                                  else errorWithoutStackTrace ("3-MultiTraverse swap option must be 'True' or 'False'" <> show inArgs)

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

                inSupport = any ((== "support") . fst) lcArgList

                -- set implied alignment swapping
                doIA' = any ((== "ia") . fst) lcArgList
                doIA'' = doIA'

                --- steepest/all options
                doSteepest' = any ((== "steepest") . fst) lcArgList
                doAll = any ((== "all") . fst) lcArgList

                -- steepest default
                doSteepest = ((not doSteepest' && not doAll) || doSteepest')

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
                    | doAnnealing || doDrift = JoinAll
                    | any ((== "joinpruned") . fst) lcArgList = JoinPruned
                    -- | any ((== "joinalternate") . fst) lcArgList = JoinPruned
                    | otherwise = JoinAll

                -- randomize split graph and rejoin edges, default to inOrder (due to sotSplit below)
                atRandom
                    | any ((== "atrandom") . fst) lcArgList = True
                    | doAnnealing || doDrift = True
                    | any ((== "inOrder") . fst) lcArgList = False
                    | any ((== "sortsplit") . fst) lcArgList = False
                    | swapType == NNI = False
                    | otherwise = False

                -- split edge order based on greartest diffenrece in costr when graph is split
                    -- does all of them before sorting
                    -- since comes after testing for random will override
                sortEdgesSplitCost
                    | doAnnealing || doDrift = False
                    | any ((== "splitsequential") . fst) lcArgList = False
                    | any ((== "sortsplit") . fst) lcArgList = True
                    | atRandom = False
                    | otherwise = True

                -- when plitting base graph--do in parallel or via recursive sequential
                -- might save on memeory, coulod be a bit more efficient time-wise
                -- definately affects trajectory--small examples had worse optimality outcomes
                parallelSplit
                    | doAnnealing || doDrift = False
                    | sortEdgesSplitCost = True
                    | any ((== "splitparallel") . fst) lcArgList = True
                    | any ((== "splitsequential") . fst) lcArgList = False
                    | otherwise = True

                -- set level of swap heristric intensity
                swapLevel
                    | all ((/= "level") . fst) lcArgList = (-1)
                    | fromJust levelNumber < 0 = 0
                    | fromJust levelNumber > 3 = 3
                    | otherwise = fromJust levelNumber

                -- populate SwapParams structure
                -- levels may have > 1 swap pass
                standardSwap = SwapParams
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

                -- levels may have > 1 swap pass
                -- and different values 
                (localSwapParams, localMultiTraverse) = if swapLevel == (-1) then
                                                            (standardSwap, doMultiTraverse) 

                                                        -- Basicall zero heuristics (adds factor of n)
                                                        else if swapLevel == 0 then
                                                            (standardSwap 
                                                                { joinType = JoinAll
                                                                , checkHeuristic = BestAll
                                                                }, True)

                                                        -- Other heuristics from 1 (slowest) -> 3 (fastest)
                                                        -- will make a call to reoptimize using MultiTraverse after swap
                                                        else if swapLevel == 1 then
                                                            (standardSwap 
                                                                { joinType = JoinAll
                                                                , checkHeuristic = BetterN
                                                                }, False)

                                                        -- This is default behavior
                                                        else if swapLevel == 2 then
                                                            (standardSwap 
                                                                { joinType = JoinPruned
                                                                , checkHeuristic = BetterN
                                                                }, False)

                                                        -- Fastest least effective
                                                        else if swapLevel == 3 then
                                                            (standardSwap 
                                                                { joinType = JoinPruned
                                                                , checkHeuristic = BestOnly
                                                                }, False)

                                                        else error ("Unimplemented swap level: " <> (show swapLevel))
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
                --action ∷ [(Maybe SAParams, ReducedPhylogeneticGraph)] → PhyG ([ReducedPhylogeneticGraph], Int)
                --action = {-# SCC swapMaster_action_swapSPRTBR #-} S.swapDriver localSwapParams inGS inData 0 inGraphList

                reoptimizeAction ∷ GlobalSettings → ProcessedData → Bool → Bool → Maybe Int → SimpleGraph → PhyG ReducedPhylogeneticGraph
                reoptimizeAction = T.multiTraverseFullyLabelGraphReduced
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

                    -- Rediagnose with MultiTraverse on if that is the setting
                        -- don't bother to if in support--only care about topology for jacknife and bootstrap
                    inGraphList' <- if inSupport then
                                        pure inGraphList
                                    else if swapLevel == (-1) then
                                            if (localMultiTraverse == multiTraverseCharacters inGS) then 
                                                pure inGraphList
                                            else do
                                                logWith LogInfo $ "\tMultiTraverse to " <> (show localMultiTraverse)  <> "\n"
                                                getParallelChunkTraverse >>= \pTraverse →
                                                    pTraverse
                                                        (reoptimizeAction (inGS{multiTraverseCharacters = localMultiTraverse}) inData False False Nothing . fst5) inGraphList

                                    -- swap level 0 uses MultiTraverse
                                    else if swapLevel == 0  && (not $ multiTraverseCharacters inGS) then do
                                                logWith LogInfo $ "\tMultiTraverse to True for swap level 0 " <> "\n"
                                                getParallelChunkTraverse >>= \pTraverse →
                                                    pTraverse
                                                        (reoptimizeAction (inGS{multiTraverseCharacters = True}) inData False False Nothing . fst5) inGraphList

                                    else if swapLevel == 0 then 
                                                pure inGraphList

                                    -- already not MultiTraverse and swap levels (1-3) that do not use it
                                    else if (not $ multiTraverseCharacters inGS) then do
                                                pure inGraphList

                                    -- is MultiTraverse and swap levels 1-3, need to reoptimize to false
                                    else do
                                                logWith LogInfo $ "\tMultiTraverse to False for swap level " <> (show swapLevel) <> "\n"
                                                getParallelChunkTraverse >>= \pTraverse →
                                                    pTraverse
                                                        (reoptimizeAction (inGS{multiTraverseCharacters = False}) inData False False Nothing . fst5) inGraphList
                                        

                    -- parallel setup
                    --action ∷ [(Maybe SAParams, ReducedPhylogeneticGraph)] → PhyG ([ReducedPhylogeneticGraph], Int)
                    let action = {-# SCC swapMaster_action_swapSPRTBR #-} SV2.swapDriver localSwapParams (inGS {multiTraverseCharacters = localMultiTraverse}) inData 0 inGraphList'

                    let simAnnealList = (: []) <$> zip newSimAnnealParamList inGraphList'
                    graphPairList ←
                        getParallelChunkTraverse >>= \pTraverse →
                            action `pTraverse` simAnnealList

                    let (graphListList, counterList) = first fold $ unzip graphPairList
                    (newGraphList, counter) ← GO.selectGraphs Best (outgroupIndex inGS) (fromJust keepNum) 0 graphListList <&> \x → (x, sum counterList)

                    let finalGraphList = case newGraphList of
                            [] → inGraphList'
                            _ → newGraphList

                    let fullBuffWarning =
                            if (length newGraphList >= (fromJust keepNum)) && (not inSupport)
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

                    -- add in second round for higher swap levels, but not if inSupport (jackknife and bootstrap)
                    if inSupport then 
                         pure finalGraphList
                    else if swapLevel == (-1) then 
                        if (localMultiTraverse == multiTraverseCharacters inGS) then 
                            pure finalGraphList
                        else 
                            do
                                logWith LogInfo $ "\tMultiTraverse to " <> (show $ multiTraverseCharacters inGS)  <> "\n"
                                getParallelChunkTraverse >>= \pTraverse →
                                            pTraverse
                                                    (reoptimizeAction inGS inData False False Nothing . fst5) finalGraphList
                    else if swapLevel == 0 then
                        if multiTraverseCharacters inGS then 
                            pure finalGraphList

                        else do
                                logWith LogInfo  "\tMultiTraverse to False\n"
                                getParallelChunkTraverse >>= \pTraverse →
                                                pTraverse
                                                    (reoptimizeAction (inGS{multiTraverseCharacters = False}) inData False False Nothing . fst5) finalGraphList
                    else if swapLevel == 1 then
                        if not $ multiTraverseCharacters inGS then 
                            pure finalGraphList
                        else do
                            logWith LogInfo  "\tMultiTraverse to True\n"
                            reDiagGraphs <- getParallelChunkTraverse >>= \pTraverse →
                                                pTraverse
                                                    (reoptimizeAction (inGS{multiTraverseCharacters = True}) inData False False Nothing . fst5) finalGraphList
                            if inputBestCost < (minimum $ fmap snd5 reDiagGraphs) then 
                                pure inGraphListInput
                            else 
                                pure reDiagGraphs

                    -- swap levels 2 and 3 are followed by a 1
                    else do
                        logWith LogInfo $ "\tSecond round level swap " <> "\n"
                        
                        let swapParamsLevel = standardSwap { joinType = JoinAll
                                                            , checkHeuristic = BetterN
                                                            }
                        let actionLevel = {-# SCC swapMaster_action_swapSPRTBR #-} SV2.swapDriver swapParamsLevel (inGS {multiTraverseCharacters = localMultiTraverse}) inData 0 finalGraphList

                        let simAnnealListLevel = (: []) <$> zip newSimAnnealParamList finalGraphList
                        graphPairListLevel ←
                            getParallelChunkTraverse >>= \pTraverse →
                                actionLevel `pTraverse` simAnnealListLevel

                        let (graphListListLevel, counterListLevel) = first fold $ unzip graphPairListLevel

                        -- Rediagnose with MultiTraverse on if that is the setting
                        reoptimizedGraphList <- if not $ multiTraverseCharacters inGS then do
                                                    logWith LogInfo "\n"
                                                    pure graphListListLevel

                                                else do
                                                    logWith LogInfo $ "\n\tMultiTraverse  to True\n"
                                                    reDiagGraphs <- getParallelChunkTraverse >>= \pTraverse →
                                                            pTraverse
                                                                (reoptimizeAction (inGS{multiTraverseCharacters = True}) inData False False Nothing . fst5) graphListListLevel
                                                    if inputBestCost < (minimum $ fmap snd5 reDiagGraphs) then 
                                                        pure inGraphListInput
                                                    else 
                                                        pure reDiagGraphs


                        (newGraphListLevel, _) ← GO.selectGraphs Best (outgroupIndex inGS) (fromJust keepNum) 0 reoptimizedGraphList <&> \x → (x, sum counterListLevel)

                        let finalGraphListLevel = case newGraphListLevel of
                                [] → finalGraphList
                                _ → newGraphListLevel

                        pure finalGraphListLevel

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
      , Bool
      , Maybe Int
      , Bool
      , Maybe Int
      , Maybe Double
      , Maybe Double
      , Maybe Int
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

                    -- simulated annealing parameters
                    doAnnealing = any ((== "annealing") . fst) lcArgList

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

                    levelList = filter ((== "level") . fst) lcArgList
                    levelNumber
                        | length levelList > 1 =
                            errorWithoutStackTrace ("Multiple 'level' number specifications in swap command--can have only one: " <> show inArgs)
                        | null levelList = Just 2
                        | otherwise = readMaybe (snd $ head levelList) ∷ Maybe Int


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
                                                                    else if isNothing levelNumber
                                                                        then
                                                                            errorWithoutStackTrace
                                                                                ("Swap 'level' specification not an integer (e.g. level:2): " <> show (snd $ head replicatesList))
                                                                        else -- trace ("GSP: " <> (show inArgs) <> " " <> (show )(keepNum, maxMoveEdgeDist', steps', annealingRounds', doDrift, driftRounds', acceptEqualProb, acceptWorseFactor, maxChanges, lcArgList))

                                                                            ( keepNum
                                                                            , maxMoveEdgeDist'
                                                                            , steps'
                                                                            , doAnnealing
                                                                            , annealingRounds'
                                                                            , doDrift
                                                                            , driftRounds'
                                                                            , acceptEqualProb
                                                                            , acceptWorseFactor
                                                                            , maxChanges
                                                                            , replicates
                                                                            , levelNumber
                                                                            , lcArgList
                                                                            )
