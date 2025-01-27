{-# OPTIONS_GHC -fmax-pmcheck-models=50 #-}

{- |
Module controlling timed randomized search functions
-}
module Search.Search (
    search,
) where

import Commands.Transform qualified as TRANS
import Commands.Verify qualified as VER
import Control.Arrow ((&&&))
import Control.Concurrent.Timeout (timeout)
import Control.DeepSeq
import Control.Monad (when)
import Control.Monad.IO.Unlift
import Control.Monad.Random.Class qualified as Sample
import Data.Bifunctor (bimap)
import Data.Char
import Data.Foldable
--import Data.Foldable1 qualified as F1
import Data.Functor (($>), (<&>))
import Data.List qualified as L
import Data.List.NonEmpty (NonEmpty (..))
import Data.List.Split qualified as LS
import Data.Maybe
import Data.Vector qualified as V
import GeneralUtilities
import Graphs.GraphOperations qualified as GO
import PHANE.Evaluation
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import Search.Build qualified as B
import Search.Refinement qualified as R
import System.Timing
import Text.Read
import Types.Types
import Utilities.LocalGraph qualified as LG
import Utilities.Utilities qualified as U


-- Bandit lists are concateenated for each of use andd update--first 3 are graph evaluation
-- and remainder are sesrch bandits.
-- they are updated separately

-- | treeBanditList is list of search types to be chosen from if graphType is tree
treeBanditList ∷ [String]
treeBanditList =
    [ "buildCharacter"
    , "buildSPR"        -- uses quickest options so other swaps should be more exhaustinve
    , "buildTBR"        -- uses quickest options so other swaps should be more exhaustinve
    , "swapSPR"
    , "swapTBR"
    , "swapAlternate"
    , "fuse"
    , "fuseSPR"
    , "fuseTBR"
    , "driftSPR"
    , "driftTBR"
    , "annealSPR"
    , "annealTBR"
    , "geneticAlgorithm"
    -- , "buildDistance" -- distance only up front to reduce memory footprint
    ]


-- | netWorkBanditList is list of search types unique to graphType network
netWorkBanditList ∷ [String]
netWorkBanditList = ["networkAdd", "networkDelete", "networkAddDelete", "driftNetwork", "annealNetwork", "networkMove"]


-- | fullBanditList is list of search types to be chosen from if Network
fullBanditList ∷ [String]
fullBanditList = treeBanditList <> netWorkBanditList


-- | graphBanditList list for more rapid, or more thorough methods
graphBanditList ∷ [String]
graphBanditList = fmap show [MultiTraverse, SingleTraverse, StaticApproximation]


-- | A strict, three-way version of 'uncurry'.
uncurry' ∷ (Functor f, NFData d) ⇒ (a → b → f d) → (a, b) → f d
uncurry' f (a, b) = force <$> f a b


{-
-- | A strict, three-way version of 'uncurry'.
uncurry3' ∷ (Functor f, NFData d) ⇒ (a → b → c → f d) → (a, b, c) → f d
uncurry3' f (a, b, c) = force <$> f a b c
-}

-- | search timed randomized search returns graph list and comment list with info String for each search instance
search
    ∷ [Argument]
    → GlobalSettings
    → ProcessedData
    → [ReducedPhylogeneticGraph]
    → PhyG ([ReducedPhylogeneticGraph], [[String]])
search inArgs inGS inData inGraphList' =
    -- flatThetaList is the initial prior list (flat) of search (bandit) choices
    -- can also be used in search for non-Thomspon search
    let flattenList :: Fractional b => [a] -> [(a, b)]
        flattenList xs =
            let count = length xs
                limit = 1 / fromIntegral count
            in  zip xs $ L.replicate count limit

        flatGraphBanditList = flattenList graphBanditList

        flatThetaList = flattenList $ case graphType inGS of
            Tree → treeBanditList
            _ → fullBanditList

        totalFlatTheta = flatGraphBanditList <> flatThetaList
    in  do
            (searchTime, keepNum, instances, thompsonSample, mFactor, mFunction, maxNetEdges, stopNum) ← getSearchParams inGS inArgs

            let threshold = fromSeconds . fromIntegral $ (100 * searchTime) `div` 100
            let initialSeconds = fromSeconds . fromIntegral $ (0 ∷ Int)
            let searchTimed ∷ (Int, ([ReducedPhylogeneticGraph], [String])) → Evaluation () ([ReducedPhylogeneticGraph], [String])
                searchTimed =
                    uncurry' $
                        searchForDuration
                            inGS
                            inData
                            [[]]
                            keepNum
                            thompsonSample
                            mFactor
                            mFunction
                            totalFlatTheta
                            1
                            maxNetEdges
                            initialSeconds
                            threshold
                            0
                            stopNum

            logWith
                LogInfo
                ( "Randomized search for "
                    <> (show searchTime)
                    <> " seconds with "
                    <> (show instances)
                    <> " instances keeping at most "
                    <> (show keepNum)
                    <> " graphs"
                    <> "\n"
                )
            -- if initial graph list is empty make some

            inGraphList ← case length inGraphList' `compare` keepNum of
                LT → do
                    dWagGraphList ←
                        B.buildGraph
                            [("distance", ""), ("replicates", show (1000)), ("rdwag", ""), ("best", show keepNum), ("return", show keepNum)]
                            inGS
                            inData

                    fmap (take keepNum) . GO.selectGraphs Unique (outgroupIndex inGS) (maxBound ∷ Int) 0 $ dWagGraphList <> inGraphList'
                _ → pure inGraphList'

            let threadCount = instances -- <- (max 1) <$> getNumCapabilities
            let startGraphs = replicate threadCount (inGraphList, mempty)
            let threadInits ∷ [(Int, ([ReducedPhylogeneticGraph], [String]))]
                threadInits = zip [1 ..] startGraphs
            -- If there are no input graphs--make some via distance
            -- resultList ← pooledMapConcurrently searchTimed threadInits
            resultList ←
                getParallelChunkTraverse >>= \pTraverse →
                    searchTimed `pTraverse` threadInits
            let (newGraphList, commentList) = unzip resultList
            let newCostList = L.group $ L.sort $ fmap getMinGraphListCost newGraphList

            let iterationHitString =
                    ( "Hit minimum cost "
                        <> (show $ minimum $ fmap snd5 $ concat newGraphList)
                        <> " in "
                        <> (show $ length $ head newCostList)
                        <> " of "
                        <> (show $ length newGraphList)
                        <> " iterations"
                        <> "\n"
                    )
            let completeGraphList = inGraphList <> fold newGraphList
            filteredGraphList ← GO.selectGraphs Unique (outgroupIndex inGS) (maxBound ∷ Int) 0 completeGraphList
            let selectedGraphList = take keepNum filteredGraphList

            logWith LogInfo iterationHitString

            pure (selectedGraphList, commentList <> [[iterationHitString]])
    where
        getMinGraphListCost ∷ (Foldable t, Functor t) ⇒ t (a, Double, c, d, e) → Double
        getMinGraphListCost a
            | not $ null a = minimum $ fmap snd5 a
            | otherwise = infinity


-- unneeded
-- instance NFData  (IO (Maybe ([ReducedPhylogeneticGraph], [String]))) where rnf x = seq x ()

{- $ toMicroseconds allotedSeconds
this CPUtime is total over all threads--not wall clock
so changed to crappier getCurrentTime in System.Timing to
get wall clock-like ellapsed time
addied time out to terminate when exceeeded time remaining.
never teminates due to time
-}


searchForDuration
    ∷ ∀ r
     . (Fractional r, Real r)
    ⇒ GlobalSettings
    → ProcessedData
    → [[VertexCost]]
    → Int
    → Bool
    → Int
    → String
    → [(String, r)]
    → Int
    → Int
    → CPUTime
    → CPUTime
    → Int
    → Int
    → Int
    → ([ReducedPhylogeneticGraph], [String])
    → PhyG ([ReducedPhylogeneticGraph], [String])
searchForDuration inGS inData pairwiseDistances keepNum thompsonSample mFactor mFunction totalThetaList counter maxNetEdges inTotalSeconds allotedSeconds stopCount stopNum refIndex (inGraphList, infoStringList) =
    let timeLimit = fromIntegral $ toMicroseconds allotedSeconds

        logWarning ∷ b → [String] → PhyG b
        logWarning val tokens = logWith LogWarn ((unwords $ "Thread" : show refIndex : tokens) <> "\n") $> val

        runForDuration ∷ PhyG a → PhyG (Maybe a)
        runForDuration = liftIOOp (timeout timeLimit)
    in  do
            -- this line to keep control of graph number
            inGraphList' ← take keepNum <$> GO.selectGraphs Unique (outgroupIndex inGS) (maxBound ∷ Int) 0.0 inGraphList

            -- searchingInnerOp ∷ PhyG ([ReducedPhylogeneticGraph], [String])
            let searchingInnerOp =
                    force $
                        performSearch
                            inGS
                            inData
                            pairwiseDistances
                            keepNum
                            totalThetaList
                            maxNetEdges
                            inTotalSeconds
                            (inGraphList', infoStringList)

            -- searchingForDuration ∷ PhyG ([ReducedPhylogeneticGraph], [String])
            let searchingForDuration = do
                    -- result = force $ performSearch inGS inData pairwiseDistances keepNum thompsonSample thetaList maxNetEdges (head seedList) inTotalSeconds (inGraphList', infoStringList)
                    result ← runForDuration searchingInnerOp
                    case result of
                        Nothing → logWarning (inGraphList, []) ["terminated due to time"]
                        Just gs → logWarning gs ["is OK", show allotedSeconds, "->", show . fromIntegral $ toMicroseconds allotedSeconds]

            -- (elapsedSeconds, output) <- timeOpUT $
            (elapsedSeconds, elapsedSecondsCPU, output) ← timeOpCPUWall searchingForDuration

            -- update theta list based on performance
            let outTotalSeconds = timeSum inTotalSeconds elapsedSecondsCPU
            let bestCost =
                    if (not $ null $ fst output)
                        then minimum $ fmap snd5 $ fst output
                        else infinity
            let finalTimeString = ",Final Values,," <> (show bestCost) <> "," <> (show $ toSeconds outTotalSeconds)

            case snd output of
                -- [] → pure output
                [] ->   do
                        let thetaString = L.intercalate "," . fmap (showRealValue . snd) $ totalThetaList
                        pure (fst output, infoStringList <> (snd output) <> [finalTimeString <> "," <> thetaString <> "," <> "*"]) 
                x : xs → do
                    -- passing time as CPU time not wall clock so parallel timings change to elapsedSeconds for wall clock
                    (updatedThetaList, newStopCount) ←
                        updateTheta
                            SearchBandit
                            thompsonSample
                            mFactor
                            mFunction
                            counter
                            (x :| xs)
                            (drop 3 totalThetaList)
                            elapsedSecondsCPU
                            outTotalSeconds
                            stopCount
                            stopNum
                    (updatedGraphTheta, _) ←
                        updateTheta
                            GraphBandit
                            thompsonSample
                            mFactor
                            mFunction
                            counter
                            (x :| xs)
                            (take 3 totalThetaList)
                            elapsedSecondsCPU
                            outTotalSeconds
                            stopCount
                            stopNum
                    -- add lists together properly so first three are the graph bandits
                    let combinedThetaList = updatedGraphTheta <> updatedThetaList
                    let thetaString = L.intercalate "," . fmap (showRealValue . snd) $ case snd output of
                            [] → totalThetaList
                            _ → combinedThetaList

                    let remainingTime = allotedSeconds `timeLeft` elapsedSeconds
                    logWith LogMore $
                        unlines
                            [ "Thread   \t" <> show refIndex
                            , "Alloted  \t" <> show allotedSeconds
                            , "Ellapsed \t" <> show elapsedSeconds
                            , "Remaining\t" <> show remainingTime
                            , "\n"
                            ]

                    if (toPicoseconds remainingTime) == 0 || newStopCount >= stopNum || (null $ snd output)
                        -- output with strings correctly added together
                        then pure (fst output, infoStringList <> (snd output) <> [finalTimeString <> "," <> thetaString <> "," <> "*"]) 
                        else
                            searchForDuration
                                inGS
                                inData
                                pairwiseDistances
                                keepNum
                                thompsonSample
                                mFactor
                                mFunction
                                combinedThetaList
                                (counter + 1)
                                maxNetEdges
                                outTotalSeconds
                                remainingTime
                                newStopCount
                                stopNum
                                refIndex
                                $ bimap (inGraphList <>) (infoStringList <>) output


-- | updateTheta updates the expected success parameters for the bandit search list
updateTheta
    ∷ ∀ r
     . (Fractional r, Real r)
    ⇒ BanditType
    → Bool
    → Int
    → String
    → Int
    → NonEmpty String
    → [(String, r)]
    → CPUTime
    → CPUTime
    → Int
    → Int
    → PhyG ([(String, r)], Int)
updateTheta thisBandit thompsonSample mFactor mFunction counter (infoString :| infoStringList) inPairList elapsedSeconds totalSeconds stopCount stopNum = case inPairList of
    [] → pure ([], stopCount)
    _ →
        let searchBandit = case thisBandit of
                SearchBandit → takeWhile (/= ',') $ tail infoString
                -- GraphBandit
                _ → case LS.splitOn "," infoString of
                    toks | "StaticApprox" `elem` toks → "StaticApproximation"
                    toks | "MultiTraverse:False" `elem` toks → "SingleTraverse"
                    _ → "MultiTraverse"

            searchDeltaString = takeWhile (/= ',') $ tail $ dropWhile (/= ',') (tail infoString)
            searchDelta ∷ r
            searchDelta = fromRational . toRational $ (read searchDeltaString ∷ Double)
        in  if not thompsonSample
                then
                    let newStopCount
                            | searchDelta <= 0.0 = stopCount + 1
                            | otherwise = 0

                        stopString
                            | newStopCount >= stopNum =
                                "\n\tSearch iterations have not improved for " <> (show newStopCount) <> " iterations--terminating this search command" <> "\n"
                            | otherwise = ""
                    in  do
                            when (thisBandit == SearchBandit) $ logWith LogInfo stopString
                            pure (inPairList, newStopCount)
                else -- update via results, previous history, memory \factor and type of memory "loss"

                    let -- get timing and benefit accounting for 0's
                        -- average time ratio in time factor for benefit adjustment
                        totalTime = (fromIntegral $ toSeconds totalSeconds) / (fromIntegral counter) ∷ r

                        timeFactor ∷ r
                        timeFactor
                            | counter == 1 = 1.0
                            | toSeconds elapsedSeconds == 0 = 0.1
                            | otherwise = (fromIntegral $ toSeconds elapsedSeconds) / totalTime

                        durationTime =
                            if toSeconds elapsedSeconds <= 1
                                then 1
                                else toSeconds elapsedSeconds

                        -- get bandit index and orignial theta value from pair list
                        indexBandit' = L.elemIndex searchBandit $ fmap fst inPairList
                        indexBandit = fromJust indexBandit'

                        inThetaBandit = snd $ inPairList !! indexBandit

                        newStopCount =
                            if searchDelta <= 0.0
                                then stopCount + 1
                                else 0

                        stopString =
                            if newStopCount >= stopNum
                                then
                                    ( "\n\tSearch iterations have not improved for " <> (show newStopCount) <> " iterations--terminating this search command" <> "\n"
                                    )
                                else ""
                    in  -- check error
                        if isNothing indexBandit'
                            then
                                error
                                    ("Bandit index not found: " <> searchBandit <> " in " <> showRealValuedPairs inPairList <> (show $ tail $ head infoStringList))
                            else -- "simple" for testing, sucess=1, no time factor, full average overall interations

                                if mFunction == "simple"
                                    then -- first Simple update based on 0 = no better. 1 == any better irrespective of magnitude or time, full memory
                                    -- all thetas (second field) sum to one.
                                    -- if didn't do anything but increment everyone else with 1/(num bandits - 1), then renormalize
                                    -- dowenweights unsiucessful bandit
                                    -- need to add more downweight if search long

                                        let previousSuccessList = fmap (* (fromIntegral counter)) $ fmap snd inPairList

                                            benefit =
                                                if searchDelta <= 0.0
                                                    then 0.0
                                                    else searchDelta / (fromIntegral durationTime)

                                            -- see if bandit was sucessful and if so set increment
                                            incrementBandit = if benefit == 0.0 then 0.0 else 1.0
                                            newBanditVal ∷ r
                                            newBanditVal = incrementBandit + ((fromIntegral counter) * inThetaBandit)

                                            -- for nonBandit if not successful increment them (basically to downweight bandit), other wise not
                                            incrementNonBanditVals =
                                                if benefit == 0.0
                                                    then 1.0 / ((fromIntegral $ length inPairList) - 1.0)
                                                    else 0.0
                                            updatedSuccessList = fmap (+ incrementNonBanditVals) previousSuccessList

                                            -- uopdate the bandit from list by splitting and rejoining
                                            firstBanditPart = take indexBandit updatedSuccessList
                                            thirdBanditPart = drop (indexBandit + 1) updatedSuccessList
                                            newSuccessList = firstBanditPart <> (newBanditVal : thirdBanditPart)
                                            totalTheta = sum newSuccessList

                                            newThetaList = (/ totalTheta) <$> newSuccessList
                                        in  do
                                                logWith LogInfo stopString
                                                pure (zip (fmap fst inPairList) newThetaList, newStopCount)
                                    else -- more complex 'recency' options

                                        if mFunction `elem` ["linear", "exponential"]
                                            then
                                                let -- weight factors for previous theta (wN-1)* \theta
                                                    -- maxed ot counter so averaging not wierd for early iterations
                                                    mFactor' ∷ r
                                                    mFactor' = min (fromIntegral counter) $ fromIntegral mFactor
                                                    (wN_1, wN) = case mFunction of
                                                        "linear" → (mFactor' / (mFactor' + 1), 1.0 / (mFactor' + 1))
                                                        "exponential" →
                                                            let f ∷ r → Double
                                                                f = fromRational . toRational
                                                                g ∷ Double → r
                                                                g = fromRational . toRational
                                                                e ∷ Double
                                                                e = f mFactor'
                                                                v ∷ r
                                                                v = 1.0 / g (2.0 ** e)
                                                            in  (1.0 - v, v)
                                                        _ →
                                                            error $
                                                                unwords
                                                                    [ "Thompson search option"
                                                                    , mFunction
                                                                    , "not recognized:"
                                                                    , show ["simple", "linear", "exponential"]
                                                                    ]

                                                    -- simple "success-based" benefit, scaled to average time of search iteration
                                                    searchBenefit ∷ r
                                                    searchBenefit
                                                        | searchDelta <= 0.0 = 0.0
                                                        | otherwise = searchDelta / timeFactor

                                                    previousSuccessList = fmap (* (fromIntegral counter)) $ fmap snd inPairList

                                                    -- average of new am]nd previous if m=1, no memory if m=0, longer memory with larger m
                                                    -- m should be limited to counter
                                                    -- linear ((m * previous value) + new value) / m+1
                                                    -- exponential ((2^(-m)  * previous value) + new value) / (2^(-m) + 1)
                                                    newBanditVal ∷ r
                                                    newBanditVal
                                                        | searchBenefit > inThetaBandit = (wN * searchBenefit) + (wN_1 * inThetaBandit)
                                                        | otherwise = (wN_1 * inThetaBandit) + (wN * (inThetaBandit + searchBenefit))

                                                    -- for nonBandit if not successful increment them (basically to downweight bandit), other wise not
                                                    incrementNonBanditVals
                                                        | searchDelta <= 0.0 = 1.0 / ((fromIntegral $ length inPairList) - 1.0)
                                                        | otherwise = 0.0

                                                    updatedSuccessList = fmap (+ incrementNonBanditVals) previousSuccessList

                                                    -- uopdate the bandit from list by splitting and rejoining, then normalizing to 1.0
                                                    firstBanditPart ∷ [r]
                                                    firstBanditPart = take indexBandit updatedSuccessList
                                                    thirdBanditPart ∷ [r]
                                                    thirdBanditPart = drop (indexBandit + 1) updatedSuccessList
                                                    newSuccessList ∷ [r]
                                                    newSuccessList = firstBanditPart <> (newBanditVal : thirdBanditPart)
                                                    totalTheta = sum newSuccessList

                                                    newThetaList = (/ totalTheta) <$> newSuccessList
                                                in  do
                                                        -- trace ("Update : \n" <> (show $ fmap snd inPairList) <> "\n" <> (show previousSuccessList) <> "\n" <> (show updatedSuccessList) <> "\n" <> (show newThetaList) <> "\n") $
                                                        -- trace ("Not simple: " <> mFunction <> " search benefit " <> (show searchBenefit) <> " " <> searchBandit <> " index " <> (show indexBandit) <> " total time: " <> (show $ toSeconds totalSeconds) <> " elapsed time: " <> (show $ toSeconds elapsedSeconds) <> " -> " <> (show (searchDelta, toSeconds elapsedSeconds)) <> "\n" <> (show $ fmap snd inPairList) <> "\n" <> (show newThetaList) <> "\n" <> (head infoStringList)) (
                                                        logWith LogInfo stopString
                                                        pure (zip (fmap fst inPairList) newThetaList, newStopCount)
                                            else -- )

                                                errorWithoutStackTrace
                                                    ("Thompson search option " <> mFunction <> " not recognized " <> (show ["simple", "linear", "exponential"]))


-- | This exponentiation functionn from http://www.haskell.org/haskellwiki/Generic_number_type#squareRoot
(^!) ∷ (Num a) ⇒ a → Int → a
(^!) x n = x ^ n


-- | squareRoot integer square root from http://www.haskell.org/haskellwiki/Generic_number_type#squareRoot
squareRoot ∷ Integer → Integer
squareRoot 0 = 0
squareRoot 1 = 1
squareRoot n =
    let twopows = iterate (^! 2) 2
        (lowerRoot, lowerN) =
            last $ takeWhile ((n >=) . snd) $ zip (1 : twopows) twopows
        newtonStep x = div (x + div n x) 2
        iters = iterate newtonStep (squareRoot (div n lowerN) * lowerRoot)
        isRoot r = r ^! 2 <= n && n < (r + 1) ^! 2
    in  head $ dropWhile (not . isRoot) iters


{- | performSearch takes input graphs and performs randomized build and search with time limit
Thompson sampling and mFactor to pick strategy from updated theta success values
the random calls return the tail of the input list to avoid long list access--can do vector since infinite
if no input graphs then do a unitary distance build to get a quick start
-}
performSearch
    ∷ ∀ r
     . (Real r)
    ⇒ GlobalSettings
    → ProcessedData
    → [[VertexCost]]
    → Int
    → [(String, r)]
    → Int
    → CPUTime
    → ([ReducedPhylogeneticGraph], [String])
    → PhyG ([ReducedPhylogeneticGraph], [String])
performSearch inGS' inData' _pairwiseDistances keepNum totalThetaList maxNetEdges inTime (inGraphList', _) =
    -- set up basic parameters for search/refine methods
    let thetaList = drop 3 totalThetaList
        numLeaves = V.length $ fst3 inData'
        -- set up log for sample
        thompsonString = "," <> showRealValuedPairs totalThetaList

        getRandomTheta = sampleRandomChoices thetaList

        -- choose search type from list with frequencies as input from searchForDuration
        -- adjust if all trees and networks chose to ensure some net stuff tried
        getSearchBandit
            | null inGraphList' = pure "buildDistance"
            | otherwise = do
                rTheta ← getRandomTheta
                case graphType inGS' of
                    Tree → pure rTheta
                    _ →
                        let thoseValues =
                                ["networkMove", "networkDelete", "driftNetwork", "annealNetwork"]
                        in  case filter (== False) $ LG.isTree . fst5 <$> inGraphList' of
                                _ : _
                                    | rTheta `notElem` thoseValues →
                                        sampleRandomChoices
                                            [("networkAdd", 0.5), ("networkAddDelete", 0.5)]
                                _ → pure rTheta

        numToCharBuild = fromInteger $ squareRoot $ toInteger numLeaves
        numToDistBuild = min 1000 (numLeaves * numLeaves)
        numDistToKeep = keepNum

        showArg a = fst a <> ":" <> snd a

        -- common swap arguments
        getSwapKeep = min keepNum <$> sampleRandomChoices [(1, 0.50), (2, 0.33), (4, 0.17)]
        genSwapOpts = do
            -- this commented because of multitraverse randomization below
            -- getSwapMultiTraverse <- sampleRandomChoices [("true", 0.33), ("false", 0.67)]
            getSwapCheck <- sampleRandomChoices [("bestonly", 0.2),("better", 0.2),("bettern", 0.60)]
            getSwapJoin <- sampleRandomChoices [("joinall", 0.50), ("joinpruned", 0.50)]
            pure [(getSwapCheck, ""), (getSwapJoin, "")]

        swapGeneralOpt = pure [("steepest", ""),("atrandom", ""), ("splitsequential", "")]
        fastSwapOpts = swapGeneralOpt <> pure [("bestonly",""), ("multitraverse","false")] -- , ("joinpruned","")]
        regSwapOpts = swapGeneralOpt <> genSwapOpts

        -- common drift arguments
        getDriftArgs = do
            maxChanges ← sampleRandomChoices [("5", 0.33), ("10", 0.34), ("20", 0.33)]
            acceptEqual ← sampleRandomChoices [("0.1", 0.5), ("0.5", 0.5)]
            acceptWorse ← sampleRandomChoices [("10.0", 0.33), ("20.0", 0.34), ("40", 0.33)]
            pure
                [ ("drift", "")
                , ("maxChanges", maxChanges)
                , ("acceptEqual", acceptEqual)
                , ("acceptWorse", acceptWorse)
                ]

        -- common annealing arguments
        getAnnealArgs = do
            tempSteps ← sampleRandomChoices [("5", 0.33), ("10", 0.34), ("20", 0.33)]
            pure [("annealing", ""), ("steps", tempSteps)]

        -- common fuse options
        getFusePairs = sampleRandomChoices [("5", 0.45), ("10", 0.45), ("20", 0.1)]
        fuseKeep = 2 * keepNum

        -- network edit options
        netGeneralArgs ∷ [(String, String)]
        netGeneralArgs = [("keep", show keepNum), ("steepest", ""), ("atRandom", ""), ("maxnetedges", show maxNetEdges)]
        netMoveArgs = ("netMove", "") : netGeneralArgs
        netAddArgs = ("netAdd", "") : netGeneralArgs
        netDelArgs = ("netDel", "") : netGeneralArgs
        netAddDelArgs = ("netAddDel", "") : netGeneralArgs
        getNetDriftAnnealMethod ∷ PhyG [(String, String)]
        getNetDriftAnnealMethod = pure . (id &&& const "") <$> sampleRandomChoices [("netAdd", 0.5), ("netDel", 0.5)] -- ,("netMove", 0.25),  ("netAddDel", 0.25),]

        -- Genetic Algorithm Arguments
        -- stops after 2 rounds with no improvement (if generations > 2)
        getGeneticAlgArgs = do
            popSize ← sampleRandomChoices [("10", 0.50), ("20", 0.25), ("40", 0.25)]
            generations ← sampleRandomChoices [("1", 1.0)] -- , "2" , "4"]
            severity ← sampleRandomChoices [("0.0", 0.33), ("1.0", 0.34), ("2.0", 0.33)]
            recombinations ← sampleRandomChoices [("10", 0.45), ("20", 0.45), ("40", 0.1)]
            pure
                [ ("popsize", popSize)
                , ("generations", generations)
                , ("severity", severity)
                , ("recombinations", recombinations)
                , ("stop", "2")
                , ("maxnetedges", show maxNetEdges)
                ]
    in  do
            searchBandit ← getSearchBandit
            -- unless fuse or genetic algorithm, only operate on "best" input graphs
            -- this to reduce memory footprint when have multiple iterations
            inGraphList'' ←
                if searchBandit `elem` ["fuse", "fuseSPR", "fuseTBR", "geneticAlgorithm"]
                    then pure inGraphList'
                    else GO.selectGraphs Best (outgroupIndex inGS') keepNum 0 inGraphList'

            -- Can't do both static approx and multitraverse:False
            let transformBy xs = TRANS.transform xs inGS' inData' inData' inGraphList''
            newDataMTF ← transformBy [("multitraverse", "false")]
            newDataSA ← transformBy [("staticapprox", [])]

            -- set graph valuation bandit
            graphEvaluationBandit ← sampleRandomChoices $ take 3 totalThetaList
            -- apply graph evaluation bandit
            let transformToStaticApproximation =
                    U.getNumberNonExactCharacters (thd3 inData') /= 0
                        && graphEvaluationBandit == "StaticApproximation"

            let transformMultiTraverse = graphEvaluationBandit == "SingleTraverse"

            let ((inGS, origData, inData, inGraphList), transformString)
                    | transformToStaticApproximation && (useIA inGS') = (newDataSA, ",StaticApprox")
                    | transformMultiTraverse = (newDataMTF, ",MultiTraverse:False")
                    | otherwise = ((inGS', inData', inData', inGraphList''), "")

            -- only doing build character now in search due to memory footprint of distnce matrix
            buildType ← case searchBandit of
                "buildCharacter" → pure "character"
                "buildDistance" → pure "distance"
                _ → sampleRandomChoices [("distance", 0.5), ("character", 0.5)]

            -- common build arguments including block and distance
            --    Tree does not use block--doesn't work very well for tree building
            buildMethod ← case inGraphList' of
                [] → pure "unitary"
                _ : _ → case graphType inGS' of
                    Tree → pure "unitary"
                    _ → sampleRandomChoices [("unitary", 0.8), ("block", 0.2)]

            let wagnerOptions = case (buildType, buildMethod) of
                    ("distance", "block") →
                        [ ("replicates", show numToCharBuild)
                        , ("rdwag", "")
                        , ("best", show (1 ∷ Int))
                        ]
                    ("distance", _) →
                        [ ("replicates", show numToDistBuild)
                        , ("rdwag", "")
                        , ("best", show numDistToKeep)
                        , ("return", show numToCharBuild)
                        ]
                    ("character", "block") →
                        [("replicates", show (1 ∷ Int))]
                    ("character", _) → []

            blockOptions ← case buildMethod of
                -- to resolve block build graphs
                "block" →
                    sampleRandomChoices [("eun", 0.5), ("cun", 0.5)] >>= \reconciliationMethod →
                        pure
                            [ ("block", "")
                            , ("atRandom", "")
                            , ("displaytrees", show numToCharBuild)
                            , (reconciliationMethod, "")
                            ]
                _ → pure []

            let builder bArgs = B.buildGraph bArgs inGS' inData'
            let attach :: b -> a -> (a, b)
                attach = flip (,)
            let selectUniqueGraphs = GO.selectGraphs Unique (outgroupIndex inGS) (maxBound ∷ Int) 0.0

            --logWith LogInfo $ "\nSearch setting (multitraverse, SttaicApprox) to " <> (show $ (multiTraverseCharacters inGS, transformToStaticApproximation)) <> "\n"

            -- bandit list with search arguments set
            -- primes (') for build to start with untransformed data
            (searchGraphs, searchArgs) ← case searchBandit of
                "buildCharacter" →
                    let buildArgs = [(buildType, "")] <> wagnerOptions <> blockOptions
                    in  attach buildArgs <$> builder buildArgs

                {- this not used in searh due to memory footprint of ditance matrix
                "buildDistance" →
                    -- search for dist builds 1000, keeps 10 best distance then selects 10 best after rediagnosis
                    -- this line in here to allow for returning lots of rediagnosed distance trees, then
                    -- reducing to unique best cost trees--but is a memory pig
                    let buildArgs = [(buildType, "")] <> wagnerOptions <> blockOptions
                    in  attach buildArgs <$> builder buildArgs
                -}
                "buildSPR" →
                    let -- build part
                        buildArgs = [(buildType, "")] <> wagnerOptions <> blockOptions
                        -- swap options
                        swapType = "spr"
                    in  -- search
                        do
                            buildGraphs ← builder buildArgs
                            buildGraphs' ← selectUniqueGraphs buildGraphs
                            swapKeep ← getSwapKeep
                            swapArgs <- fastSwapOpts
                            swapList ← R.swapMaster (swapArgs <> [(swapType, ""),("keep", show swapKeep)]) inGS inData buildGraphs'
                            pure (swapList, buildArgs <> swapArgs)
                "buildTBR" →
                    let -- build part
                        buildArgs = [(buildType, "")] <> wagnerOptions <> blockOptions
                        -- swap options
                        swapType = "tbr"
                    in  -- search
                        do
                            buildGraphs ← builder buildArgs
                            buildGraphs' ← selectUniqueGraphs buildGraphs
                            swapKeep ← getSwapKeep
                            swapArgs <- fastSwapOpts
                            swapList ← R.swapMaster (swapArgs <> [(swapType, ""),("keep", show swapKeep)]) inGS inData buildGraphs'
                            pure (swapList, buildArgs <> swapArgs)
                "buildAlternate" →
                    let -- build part
                        buildArgs = [(buildType, "")] <> wagnerOptions <> blockOptions
                        -- swap options
                        swapType = "alternate" -- default anyway
                    in  -- search
                        do
                            buildGraphs ← builder buildArgs
                            buildGraphs' ← selectUniqueGraphs buildGraphs
                            swapKeep ← getSwapKeep
                            swapArgs <- fastSwapOpts
                            swapList ← R.swapMaster (swapArgs <> [(swapType, ""),("keep", show swapKeep)]) inGS inData buildGraphs'
                            pure (swapList, buildArgs <> swapArgs)
                "swapSPR" →
                    let -- swap options
                        swapType = "spr"
                    in  -- search
                        do
                            swapKeep ← getSwapKeep
                            swapArgs <- regSwapOpts
                            swapList ← R.swapMaster (swapArgs <> [(swapType, ""),("keep", show swapKeep)]) inGS inData inGraphList
                            pure (swapList, swapArgs)
                "swapTBR" →
                    let -- swap options
                        swapType = "tbr"
                    in  -- search
                        do
                            swapKeep ← getSwapKeep
                            swapArgs <- regSwapOpts
                            swapList ← R.swapMaster (swapArgs <> [(swapType, ""),("keep", show swapKeep)]) inGS inData inGraphList
                            pure (swapList, swapArgs)
                "swapAlternate" →
                    let -- swap options
                        swapType = "alternate" -- default anyway
                    in  -- search
                        do
                            swapKeep ← getSwapKeep
                            swapArgs <- regSwapOpts
                            swapList ← R.swapMaster (swapArgs <> [(swapType, ""),("keep", show swapKeep)]) inGS inData inGraphList
                            pure (swapList, swapArgs)
                -- drift only best graphs
                "driftSPR" →
                    let -- swap args
                        swapType = "spr"
                    in  -- perform search
                        do
                            driftArgs ← getDriftArgs
                            swapKeep ← getSwapKeep
                            swapArgs <- fastSwapOpts 
                            
                            -- swap with drift (common) arguments
                            let swapDriftArgs = swapArgs <> [(swapType, ""),("keep", show swapKeep)] <> driftArgs
                            swapList ← R.swapMaster swapDriftArgs inGS inData inGraphList
                            pure (swapList, swapArgs)
                -- drift only best graphs
                "driftTBR" →
                    let -- swap args
                        swapType = "tbr"
                    in  -- perform search
                        do
                            driftArgs ← getDriftArgs
                            swapKeep ← getSwapKeep
                            swapArgs <- fastSwapOpts

                            -- swap with drift (common) arguments
                            let swapDriftArgs = swapArgs <> [(swapType, ""),("keep", show swapKeep)] <> driftArgs
                            swapList ← R.swapMaster swapDriftArgs inGS inData inGraphList
                            pure (swapList, swapDriftArgs)
                -- anneal only best graphs
                "annealSPR" →
                    let -- swap args
                        swapType = "spr"
                    in  -- perform search
                        do
                            annealArgs ← getAnnealArgs
                            swapKeep ← getSwapKeep
                            swapArgs <- fastSwapOpts

                            -- swap with anneal (common) arguments
                            let swapAnnealArgs = swapArgs <> [(swapType, ""),("keep", show swapKeep)] <> annealArgs
                            swapList ← R.swapMaster swapAnnealArgs inGS inData inGraphList
                            pure (swapList, swapAnnealArgs)
                -- anneal only best graphs
                "annealTBR" →
                    let -- swap args
                        swapType = "tbr"
                    in  -- perform search
                        do
                            annealArgs ← getAnnealArgs
                            swapKeep ← getSwapKeep
                            swapArgs <- fastSwapOpts

                            -- swap with anneal (common) arguments
                            let swapAnnealArgs = swapArgs <> [(swapType, ""),("keep", show swapKeep)] <> annealArgs
                            swapList ← R.swapMaster swapAnnealArgs inGS inData inGraphList
                            pure (swapList, swapAnnealArgs)
                "geneticAlgorithm" →
                    do
                        -- args from above
                        -- perform search
                        gaArgs ← getGeneticAlgArgs
                        gaReturn ← R.geneticAlgorithmMaster gaArgs inGS inData inGraphList
                        pure (gaReturn, gaArgs)
                "fuse" → do
                    -- should more graphs be added if only one?  Would downweight fuse perhpas too much
                    fusePairs ← getFusePairs
                    -- fuse arguments
                    -- this to limit memory footprint of fuse during search
                    -- gsNum = min (graphsSteepest inGS) 5
                    -- inGSgs1 = inGS{graphsSteepest = gsNum}
                    let fuseArgs =
                            [ ("none", "")
                            , ("unique", "")
                            , ("atrandom", "")
                            , ("pairs", fusePairs)
                            , ("keep", show fuseKeep)
                            , ("noreciprocal", "")
                            , ("multitraverse", show $ multiTraverseCharacters inGS)
                            , ("maxparallel", "false")
                            ]
                    -- perform search
                    R.fuseGraphs fuseArgs inGS inData inGraphList <&> (\x → (x, fuseArgs))
                "fuseSPR" → do
                    -- fuse arguments
                    -- inGSgs1 = inGS{graphsSteepest = 1}
                    fusePairs ← getFusePairs
                    let fuseArgs =
                            [ ("spr", "")
                            , ("unique", "")
                            , ("atrandom", "")
                            , ("pairs", fusePairs)
                            , ("keep", show fuseKeep)
                            , ("noreciprocal", "")
                            , ("multitraverse", show $ multiTraverseCharacters inGS)
                            , ("maxparallel", "false")
                            ]
                    -- perform search
                    R.fuseGraphs fuseArgs inGS inData inGraphList <&> (\x → (x, fuseArgs))
                "fuseTBR" → do
                    -- fuse arguments
                    -- inGSgs1 = inGS{graphsSteepest = 1}
                    fusePairs ← getFusePairs
                    let fuseArgs =
                            [ ("tbr", "")
                            , ("unique", "")
                            , ("atrandom", "")
                            , ("pairs", fusePairs)
                            , ("keep", show fuseKeep)
                            , ("noreciprocal", "")
                            , ("multitraverse", show $ multiTraverseCharacters inGS)
                            , ("maxparallel", "false")
                            ]
                    -- perform search
                    R.fuseGraphs fuseArgs inGS inData inGraphList <&> (\x → (x, fuseArgs))
                "networkAdd" →
                    let -- network add args
                        netEditArgs = netAddArgs
                    in  -- perform search
                        R.netEdgeMaster netEditArgs inGS inData inGraphList <&> (\x → (x, netEditArgs))
                "networkDelete" →
                    let -- network delete args
                        netEditArgs = netDelArgs
                    in  -- perform search
                        R.netEdgeMaster netEditArgs inGS inData inGraphList <&> (\x → (x, netEditArgs))
                "networkAddDelete" →
                    let -- network add/delete args
                        netEditArgs = netAddDelArgs
                    in  -- perform search
                        R.netEdgeMaster netEditArgs inGS inData inGraphList <&> (\x → (x, netEditArgs))
                "networkMove" →
                    let -- network move args
                        netEditArgs = netMoveArgs
                    in  -- perform search
                        R.netEdgeMaster netEditArgs inGS inData inGraphList <&> (\x → (x, netEditArgs))
                "driftNetwork" → do
                    driftArgs ← getDriftArgs
                    netDriftAnnealArgs ← getNetDriftAnnealMethod
                    -- network add/delete  + drift args
                    let netEditArgs = fold [netDriftAnnealArgs, netGeneralArgs, driftArgs]
                    -- perform search
                    R.netEdgeMaster netEditArgs inGS inData inGraphList <&> (\x → (x, netEditArgs))
                "annealNetwork" → do
                    annealArgs ← getAnnealArgs
                    netDriftAnnealArgs ← getNetDriftAnnealMethod
                    -- network add/delete  + annealing  args
                    let netEditArgs = fold [netDriftAnnealArgs, netGeneralArgs, annealArgs]
                    -- perform search
                    R.netEdgeMaster netEditArgs inGS inData inGraphList <&> (\x → (x, netEditArgs))
                _ → error ("Unknown/unimplemented method in search: " <> searchBandit)

            -- process
            uniqueGraphs' ← fmap (take keepNum) . selectUniqueGraphs $ searchGraphs <> inGraphList
            let transformBy' xs = TRANS.transform xs inGS origData inData uniqueGraphs'
            newDataMT ← transformBy' [("multiTraverse", "true")]
            newDataT ← transformBy' [("dynamic", [])]
            let (uniqueGraphs, transString)
                    | (not transformToStaticApproximation && not transformMultiTraverse) = (uniqueGraphs', "")
                    | transformToStaticApproximation = (fth4 newDataT, ",Dynamic")
                    | otherwise = (fth4 newDataMT, ",MultiTraverse:True")

            -- string of delta and cost of graphs
            let extract :: (Foldable f, Functor f, Ord b)  => f (a, b, c, d, e) -> b
                extract = minimum . fmap snd5

            let deltaString
                    | null inGraphList' = "10.0,"
                    | otherwise = show $ extract inGraphList' - extract uniqueGraphs

            let currentBestString
                    | not $ null uniqueGraphs = show $ extract uniqueGraphs
                    | otherwise = show infinity

            -- create string for search stats
            let searchString =
                    ","
                        <> searchBandit
                        <> ","
                        <> deltaString
                        <> ","
                        <> currentBestString
                        <> ","
                        <> (show $ toSeconds inTime)
                        <> ","
                        <> (L.intercalate "," $ fmap showArg searchArgs)
                        <> transformString
                        <> transString

            pure (uniqueGraphs, [searchString <> thompsonString])


{- |
'getSearchParams' takes arguments and returns search params.
-}
getSearchParams ∷ GlobalSettings -> [Argument] → PhyG (Int, Int, Int, Bool, Int, String, Int, Int)
getSearchParams inGS inArgs =
    let fstArgList = fmap (fmap toLower . fst) inArgs
        sndArgList = fmap (fmap toLower . snd) inArgs
        lcArgList = zip fstArgList sndArgList
        checkCommandList = checkCommandArgs "search" fstArgList VER.searchArgList
    in  -- check for valid command options
        if not checkCommandList
            then errorWithoutStackTrace ("Unrecognized command in 'search': " <> show inArgs)
            else
                let instancesList = filter ((== "instances") . fst) lcArgList
                    instances
                        | length instancesList > 1 =
                            errorWithoutStackTrace ("Multiple 'keep' number specifications in search command--can have only one: " <> show inArgs)
                        | null instancesList = Just 1
                        | otherwise = readMaybe (snd $ head instancesList) ∷ Maybe Int

                    keepList = filter ((== "keep") . fst) lcArgList
                    keepNum
                        | length keepList > 1 =
                            errorWithoutStackTrace ("Multiple 'keep' number specifications in search command--can have only one: " <> show inArgs)
                        | null keepList = Just $ keepGraphs inGS
                        | otherwise = readMaybe (snd $ head keepList) ∷ Maybe Int

                    daysList = filter ((== "days") . fst) lcArgList
                    days
                        | length daysList > 1 =
                            errorWithoutStackTrace ("Multiple 'days' number specifications in search command--can have only one: " <> show inArgs)
                        | null daysList = Just 0
                        | otherwise = readMaybe (snd $ head daysList) ∷ Maybe Int

                    hoursList = filter ((== "hours") . fst) lcArgList
                    hours
                        | length hoursList > 1 =
                            errorWithoutStackTrace ("Multiple 'hours' number specifications in search command--can have only one: " <> show inArgs)
                        | null hoursList = Just 0
                        | otherwise = readMaybe (snd $ head hoursList) ∷ Maybe Int

                    minutesList = filter ((== "minutes") . fst) lcArgList
                    minutes
                        | length minutesList > 1 =
                            errorWithoutStackTrace ("Multiple 'minutes' number specifications in search command--can have only one: " <> show inArgs)
                        | null minutesList = Just 1
                        | otherwise = readMaybe (snd $ head minutesList) ∷ Maybe Int

                    secondsList = filter ((== "seconds") . fst) lcArgList
                    seconds
                        | length secondsList > 1 =
                            errorWithoutStackTrace ("Multiple 'seconds' number specifications in search command--can have only one: " <> show inArgs)
                        | null secondsList = Just 0
                        | otherwise = readMaybe (snd $ head secondsList) ∷ Maybe Int

                    maxNetEdgesList = filter ((== "maxnetedges") . fst) lcArgList
                    maxNetEdges
                        | length maxNetEdgesList > 1 =
                            errorWithoutStackTrace ("Multiple 'maxNetEdges' number specifications in netEdge command--can have only one: " <> show inArgs)
                        | null maxNetEdgesList = Just 10
                        | otherwise = readMaybe (snd $ head maxNetEdgesList) ∷ Maybe Int

                    thompsonList = filter ((== "thompson") . fst) lcArgList
                    mFactor
                        | length thompsonList > 1 =
                            errorWithoutStackTrace ("Multiple 'Thompson' number specifications in search command--can have only one: " <> show inArgs)
                        | null thompsonList = Just 1
                        | otherwise = readMaybe (snd $ head thompsonList) ∷ Maybe Int

                    stopList = filter ((== "stop") . fst) lcArgList
                    stopNum
                        | length stopList > 1 =
                            errorWithoutStackTrace ("Multiple 'stop' number specifications in search command--can have only one: " <> show inArgs)
                        | null stopList = Just (maxBound ∷ Int)
                        | otherwise = readMaybe (snd $ head stopList) ∷ Maybe Int

                    thompson = any ((== "thompson") . fst) lcArgList
                    mLinear = any ((== "linear") . fst) lcArgList
                    mExponential = any ((== "exponential") . fst) lcArgList
                    mSimple = any ((== "simple") . fst) lcArgList

                    mFunction =
                        if mLinear && mExponential
                            then "linear"
                            else
                                if mLinear
                                    then "linear"
                                    else
                                        if mExponential
                                            then "exponential"
                                            else
                                                if mSimple
                                                    then "simple"
                                                    else "linear"
                in  if isNothing keepNum
                        then errorWithoutStackTrace ("Keep specification not an integer in search: " <> show (head keepList))
                        else
                            if isNothing instances
                                then errorWithoutStackTrace ("Instances specification not an integer in search: " <> show (head instancesList))
                                else
                                    if isNothing days
                                        then errorWithoutStackTrace ("Days specification not an integer in search: " <> show (head daysList))
                                        else
                                            if isNothing hours
                                                then errorWithoutStackTrace ("Hours specification not an integer in search: " <> show (head hoursList))
                                                else
                                                    if isNothing minutes
                                                        then errorWithoutStackTrace ("Minutes specification not an integer in search: " <> show (head minutesList))
                                                        else
                                                            if isNothing seconds
                                                                then errorWithoutStackTrace ("Seconds factor specification not an integer in search: " <> show (head secondsList))
                                                                else
                                                                    if isNothing mFactor
                                                                        then
                                                                            errorWithoutStackTrace
                                                                                ("Thompson mFactor specification not an integer or not found in search (e.g. Thompson:1) " <> show (head thompsonList))
                                                                        else
                                                                            if isNothing maxNetEdges
                                                                                then
                                                                                    errorWithoutStackTrace
                                                                                        ("Search 'maxNetEdges' specification not an integer or not found (e.g. maxNetEdges:8): " <> show (snd $ head maxNetEdgesList))
                                                                                else
                                                                                    if isNothing stopNum
                                                                                        then
                                                                                            errorWithoutStackTrace
                                                                                                ("Search stop specification not an integer or not found in search (e.g. stop:10) " <> show (head stopList))
                                                                                        else
                                                                                            let seconds' =
                                                                                                    if ((fromJust minutes > 0) || (fromJust hours > 0) || (fromJust days > 0)) && (null secondsList)
                                                                                                        then Just 0
                                                                                                        else seconds
                                                                                                searchTime = (fromJust seconds') + (60 * (fromJust minutes)) + (3600 * (fromJust hours))
                                                                                            in  do
                                                                                                    when (mLinear && mExponential) $
                                                                                                        logWith LogWarn ("Thompson recency function specification has both 'linear' and 'exponential', defaulting to 'linear'\n")
                                                                                                    pure $
                                                                                                        ( searchTime
                                                                                                        , fromJust keepNum
                                                                                                        , fromJust instances
                                                                                                        , thompson
                                                                                                        , fromJust mFactor
                                                                                                        , mFunction
                                                                                                        , fromJust maxNetEdges
                                                                                                        , fromJust stopNum
                                                                                                        )


showRealValue ∷ ∀ r. (Real r) ⇒ r → String
showRealValue =
    let convert ∷ r → Double
        convert = fromRational . toRational
    in  show . convert


showRealValuedPairs ∷ ∀ r. (Real r) ⇒ [(String, r)] → String
showRealValuedPairs =
    let asTuple ∷ (String, r) → String
        asTuple (x, y) = fold ["( ", show x, ", ", showRealValue y, " )"]

        enclose ∷ String → String
        enclose = ("[ " <>) . (<> " ]")
    in  enclose . L.intercalate ", " . fmap asTuple


sampleRandomChoices ∷ ∀ a r. (Real r) ⇒ [(a, r)] → PhyG a
sampleRandomChoices = Sample.fromList . fmap (fmap toRational)
