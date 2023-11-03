{- |
Module controlling timed randomized search functions
-}
module Search.Search (
    search,
) where

import Commands.Transform qualified as TRANS
import Commands.Verify qualified as VER
import Control.Concurrent.Async
import Control.DeepSeq
import Control.Evaluation
import Control.Exception
import Control.Monad (join, when)
import Control.Monad.IO.Class (MonadIO (..))
import Control.Monad.IO.Unlift
import Control.Monad.Logger (LogLevel (..), Logger (..), Verbosity (..))
import Data.Bifunctor (bimap)
import Data.Char
import Data.Functor (($>), (<&>))
import Data.Foldable
import Data.List qualified as L
import Data.List.Split qualified as LS
import Data.Maybe
-- import qualified GraphOptimization.Traversals as T
import Data.Vector qualified as V
import GeneralUtilities
import Graphs.GraphOperations qualified as GO
import Search.Build qualified as B
import Search.Refinement qualified as R
import System.ErrorPhase (ErrorPhase (..))
import System.IO
import System.Random
import System.Timeout
import System.Timing
import Text.Read
import Types.Types
import UnliftIO.Async (pooledMapConcurrently)
import Utilities.LocalGraph qualified as LG
import Utilities.Utilities qualified as U
-- import Debug.Trace


-- Bandit lists are concateenated for each of use andd update--first 3 are graph evaluation
-- and remainder are sesrch bandits.
-- they are updated separately

-- | treeBanditList is list of search types to be chosen from if graphType is tree
treeBanditList ∷ [String]
treeBanditList =
    [ "buildCharacter"
    , "buildDistance" -- "buildSPR", "buildAlternate",
    , "swapSPR"
    , "swapAlternate"
    , "fuse"
    , "fuseSPR"
    , "fuseTBR"
    , "driftSPR"
    , "driftAlternate"
    , "annealSPR"
    , "annealAlternate"
    , "geneticAlgorithm"
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
uncurry3' ∷ (Functor f, NFData d) ⇒ (a → b → c → f d) → (a, b, c) → f d
uncurry3' f (a, b, c) = force <$> f a b c


-- | search timed randomized search returns graph list and comment list with info String for each search instance
search
    ∷ [Argument]
    → GlobalSettings
    → ProcessedData
    → [[VertexCost]]
    → Int
    → [ReducedPhylogeneticGraph]
    → PhyG ([ReducedPhylogeneticGraph], [[String]])
search inArgs inGS inData pairwiseDistances rSeed inGraphList' = 

    -- flatThetaList is the initial prior list (flat) of search (bandit) choices
    -- can also be used in search for non-Thomspon search
    let flatThetaList =
            if graphType inGS == Tree
                then zip treeBanditList (L.replicate (length treeBanditList) (1.0 / (fromIntegral $ length treeBanditList)))
            else zip fullBanditList (L.replicate (length fullBanditList) (1.0 / (fromIntegral $ length fullBanditList)))

        flatGraphBanditList = zip graphBanditList (L.replicate (length graphBanditList) (1.0 / (fromIntegral $ length graphBanditList)))

        totalFlatTheta = flatGraphBanditList <> flatThetaList
    in do

        (searchTime, keepNum, instances, thompsonSample, mFactor, mFunction, maxNetEdges, stopNum) <- getSearchParams inArgs


        let threshold = fromSeconds . fromIntegral $ (100 * searchTime) `div` 100
        let initialSeconds = fromSeconds . fromIntegral $ (0 ∷ Int)
        let searchTimed = uncurry3' $
                searchForDuration
                    inGS
                    inData
                    pairwiseDistances
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
        let infoIndices = [1 ..]
        let seadStreams = randomIntList <$> randomIntList rSeed
    
        logWith LogInfo ("Randomized seach for " <> (show searchTime) <> " seconds with " <> (show instances) <> " instances keeping at most " <> (show keepNum) <> " graphs" <> "\n")
        -- if initial graph list is empty make some
        dWagGraphList ←
            B.buildGraph
                        [("distance", ""), ("replicates", show (1000)), ("rdwag", ""), ("best", show keepNum), ("return", show keepNum)]
                        inGS
                        inData
                        pairwiseDistances
                        rSeed
        let inGraphList =
                if (not . null) inGraphList'
                    then inGraphList'
                else take keepNum $ GO.selectGraphs Unique (maxBound ∷ Int) 0.0 (-1) dWagGraphList

        --  threadCount <- (max 1) <$> getNumCapabilities
        let threadCount = instances -- <- (max 1) <$> getNumCapabilities
        let startGraphs = replicate threadCount (inGraphList, mempty)
        let threadInits = zip3 infoIndices seadStreams startGraphs
        resultList ← pooledMapConcurrently searchTimed threadInits -- If there are no input graphs--make some via distance
        let (newGraphList, commentList) = unzip resultList
        let newCostList = L.group $ L.sort $ fmap getMinGraphListCost newGraphList

        let iterationHitString = ("Hit minimum cost " <> (show $ minimum $ fmap snd5 $ concat newGraphList) <> " in " <> (show $ length $ head newCostList) <> " of " <> (show $ length newGraphList) <> " iterations" <> "\n")
        let completeGraphList = inGraphList <> fold newGraphList
        let filteredGraphList = GO.selectGraphs Unique (maxBound ∷ Int) 0.0 (-1) completeGraphList
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
    ∷ GlobalSettings
    → ProcessedData
    → [[VertexCost]]
    → Int
    → Bool
    → Int
    → String
    → [(String, Double)]
    → Int
    → Int
    → CPUTime
    → CPUTime
    → Int
    → Int
    → Int
    → [Int]
    → ([ReducedPhylogeneticGraph], [String])
    → PhyG ([ReducedPhylogeneticGraph], [String])
searchForDuration inGS inData pairwiseDistances keepNum thompsonSample mFactor mFunction totalThetaList counter maxNetEdges inTotalSeconds allotedSeconds stopCount stopNum refIndex seedList (inGraphList, infoStringList) =
    let timeLimit = fromIntegral $ toMicroseconds allotedSeconds

        -- this line to keep control of graph number
        inGraphList' = take keepNum $ GO.selectGraphs Unique (maxBound ∷ Int) 0.0 (-1) inGraphList

        logWarning :: b -> [String] -> PhyG b
        logWarning val tokens = logWith LogInfo ((unwords $ "Thread" : show refIndex : tokens)  <> "\n") $> val

        runForDuration :: PhyG a -> PhyG (Maybe a)
        runForDuration = liftIOOp (timeout timeLimit)

        searchingInnerOp :: PhyG ([ReducedPhylogeneticGraph], [String])
        searchingInnerOp = force $ performSearch
            inGS
            inData
            pairwiseDistances
            keepNum
            thompsonSample
            totalThetaList
            maxNetEdges
            (head seedList)
            inTotalSeconds
            (inGraphList', infoStringList)

        searchingForDuration :: PhyG ([ReducedPhylogeneticGraph], [String])
        searchingForDuration = do
            -- result = force $ performSearch inGS inData pairwiseDistances keepNum thompsonSample thetaList maxNetEdges (head seedList) inTotalSeconds (inGraphList', infoStringList)
            result <- runForDuration searchingInnerOp 
            case result of
                Nothing -> logWarning (inGraphList, []) ["terminated due to time" ]
                Just gs -> logWarning gs [ "is OK", show allotedSeconds, "->", show . fromIntegral $ toMicroseconds allotedSeconds]
    in  do
            -- (elapsedSeconds, output) <- timeOpUT $
            (elapsedSeconds, elapsedSecondsCPU, output) ← timeOpCPUWall searchingForDuration

            -- update theta list based on performance
            let outTotalSeconds = timeSum inTotalSeconds elapsedSecondsCPU
            let bestCost =
                    if (not $ null $ fst output)
                        then minimum $ fmap snd5 $ fst output
                        else infinity
            let finalTimeString = ",Final Values,," <> (show bestCost) <> "," <> (show $ toSeconds outTotalSeconds)
            -- passing time as CPU time not wall clock so parallel timings change to elapsedSeconds for wall clock
            let (updatedThetaList, newStopCount) =
                    updateTheta
                        SearchBandit
                        thompsonSample
                        mFactor
                        mFunction
                        counter
                        (snd output)
                        (drop 3 totalThetaList)
                        elapsedSecondsCPU
                        outTotalSeconds
                        stopCount
                        stopNum
            let (updatedGraphTheta, _) =
                    updateTheta
                        GraphBandit
                        thompsonSample
                        mFactor
                        mFunction
                        counter
                        (snd output)
                        (take 3 totalThetaList)
                        elapsedSecondsCPU
                        outTotalSeconds
                        stopCount
                        stopNum
            -- add lists together properly so first three are the graph bandits
            let combinedThetaList = updatedGraphTheta <> updatedThetaList
            let thetaString =
                    if (null $ snd output)
                        then L.intercalate "," $ fmap (show . snd) totalThetaList
                        else L.intercalate "," $ fmap (show . snd) combinedThetaList

            let remainingTime = allotedSeconds `timeLeft` elapsedSeconds
            logWith LogInfo $ unlines
                [ "Thread   \t" <> show refIndex
                , "Alloted  \t" <> show allotedSeconds
                , "Ellapsed \t" <> show elapsedSeconds
                , "Remaining\t" <> show remainingTime
                , "\n"
                ]

            if (toPicoseconds remainingTime) == 0 || newStopCount >= stopNum || (null $ snd output)
                then pure (fst output, infoStringList <> (snd output) <> [finalTimeString <> "," <> thetaString <> "," <> "*"]) -- output with strings correctly added together
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
                        (tail seedList) $
                        bimap (inGraphList <>) (infoStringList <>) output

{-Momadic  version has a bug with logging search-}
{-
-- | updateTheta updates the expected success parameters for the bandit search list
updateTheta
    ∷ BanditType
    → Bool
    → Int
    → String
    → Int
    → [String]
    → [(String, Double)]
    → CPUTime
    → CPUTime
    → Int
    → Int
    → PhyG ([(String, Double)], Int)
updateTheta thisBandit thompsonSample mFactor mFunction counter infoStringList inPairList elapsedSeconds totalSeconds stopCount stopNum =
    if null inPairList
        then do
            pure ([], stopCount)
        else
            let searchBandit =
                    if thisBandit == SearchBandit
                        then takeWhile (/= ',') (tail $ head infoStringList)
                    else -- GraphBandit

                            if "StaticApprox" `elem` (LS.splitOn "," $ head infoStringList)
                                then 
                                    --trace
                                    --    ("GraphBandit is StaticApprox")
                                        "StaticApproximation"
                                else
                                    if "MultiTraverse:False" `elem` (LS.splitOn "," $ head infoStringList)
                                        then
                                            --trace
                                            --    ("GraphBandit is SingleTraverse")
                                                "SingleTraverse"
                                        else
                                            --trace
                                            --    ("GraphBandit is MultiTraverse ")
                                                "MultiTraverse"
                searchDeltaString = takeWhile (/= ',') $ tail $ dropWhile (/= ',') (tail $ head infoStringList)
                searchDelta = read searchDeltaString ∷ Double
            in  if not thompsonSample
                    then
                        let newStopCount =
                                if searchDelta <= 0.0
                                    then stopCount + 1
                                    else 0
                            stopString =
                                 if newStopCount >= stopNum
                                    then ("\n\tSearch iterations have not improved for " <> (show newStopCount) <> " iterations--terminating this search command" <> "\n")
                                 else ""
                        in  if thisBandit == SearchBandit
                                then do
                                    logWith LogInfo stopString
                                    pure (inPairList, newStopCount)
                            else do
                                    pure (inPairList, newStopCount)
                    else -- trace ("UT1: " <> (show infoStringList)) (
                    -- update via results, previous history, memory \factor and type of memory "loss"

                        let -- get timing and benefit accounting for 0's
                            -- average time ratio in time factor for benefit adjustment
                            totalTime = (fromIntegral $ toSeconds totalSeconds) / (fromIntegral counter) ∷ Double
                            timeFactor =
                                if counter == 1
                                    then 1.0
                                    else
                                        if toSeconds elapsedSeconds == 0
                                            then 0.1
                                            else (fromIntegral $ toSeconds elapsedSeconds) / totalTime

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
                                    then ("\n\tSearch iterations have not improved for " <> (show newStopCount) <> " iterations--terminating this search command" <> "\n")
                                    else ""
                        in  -- check error
                            if isNothing indexBandit'
                                then error ("Bandit index not found: " <> searchBandit <> " in " <> (show inPairList) <> (show $ tail $ head infoStringList))
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

                                                newThetaList = fmap (/ totalTheta) newSuccessList
                                            in  do -- trace ("UT2: " <> searchBandit <> " index " <> (show indexBandit) <> " total time: " <> (show $ toSeconds totalSeconds) <> " elapsed time: " <> (show $ toSeconds elapsedSeconds) <> " -> " <> (show (searchDelta, toSeconds elapsedSeconds, benefit)) <> "\n" <> (show $ fmap snd inPairList) <> "\n" <> (show newThetaList) <> "\n" <> (head infoStringList)) (
                                                logWith LogInfo stopString
                                                pure (zip (fmap fst inPairList) newThetaList, newStopCount)
                                        else -- )

                                        -- more complex 'recency' options

                                            if mFunction `elem` ["linear", "exponential"]
                                                then
                                                    let -- weight factors for previous theta (wN-1)* \theta
                                                        -- maxed ot counter so averaging not wierd for early iterations
                                                        mFactor' = min (fromIntegral counter ∷ Double) (fromIntegral mFactor ∷ Double)
                                                        (wN_1, wN) =
                                                            if mFunction == "linear"
                                                                then (mFactor' / (mFactor' + 1), 1.0 / (mFactor' + 1))
                                                                else
                                                                    if mFunction == "exponential"
                                                                        then (1.0 - (1.0 / (2.0 ** mFactor')), (1.0 / (2.0 ** mFactor')))
                                                                        else error ("Thompson search option " <> mFunction <> " not recognized " <> (show ["simple", "linear", "exponential"]))

                                                        -- simple "success-based" benefit, scaled to average time of search iteration
                                                        searchBenefit =
                                                            if searchDelta <= 0.0
                                                                then 0.0
                                                                else searchDelta / timeFactor

                                                        previousSuccessList = fmap (* (fromIntegral counter)) $ fmap snd inPairList

                                                        -- average of new am]nd previous if m=1, no memory if m=0, longer memory with larger m
                                                        -- m should be limited to counter
                                                        -- linear ((m * previous value) + new value) / m+1
                                                        -- exponential ((2^(-m)  * previous value) + new value) / (2^(-m) + 1)
                                                        newBanditVal =
                                                            if searchBenefit > inThetaBandit
                                                                then (wN * searchBenefit) + (wN_1 * inThetaBandit)
                                                                else (wN_1 * inThetaBandit) + (wN * (inThetaBandit + searchBenefit))

                                                        -- for nonBandit if not successful increment them (basically to downweight bandit), other wise not
                                                        incrementNonBanditVals =
                                                            if searchDelta <= 0.0
                                                                then 1.0 / ((fromIntegral $ length inPairList) - 1.0)
                                                                else 0.0
                                                        updatedSuccessList = fmap (+ incrementNonBanditVals) previousSuccessList

                                                        -- uopdate the bandit from list by splitting and rejoining, then normalizing to 1.0
                                                        firstBanditPart = take indexBandit updatedSuccessList
                                                        thirdBanditPart = drop (indexBandit + 1) updatedSuccessList
                                                        newSuccessList = firstBanditPart <> (newBanditVal : thirdBanditPart)
                                                        totalTheta = sum newSuccessList

                                                        newThetaList = fmap (/ totalTheta) newSuccessList
                                                    in do -- trace ("Update : \n" <> (show $ fmap snd inPairList) <> "\n" <> (show previousSuccessList) <> "\n" <> (show updatedSuccessList) <> "\n" <> (show newThetaList) <> "\n") $
                                                        -- trace ("Not simple: " <> mFunction <> " search benefit " <> (show searchBenefit) <> " " <> searchBandit <> " index " <> (show indexBandit) <> " total time: " <> (show $ toSeconds totalSeconds) <> " elapsed time: " <> (show $ toSeconds elapsedSeconds) <> " -> " <> (show (searchDelta, toSeconds elapsedSeconds)) <> "\n" <> (show $ fmap snd inPairList) <> "\n" <> (show newThetaList) <> "\n" <> (head infoStringList)) (
                                                        logWith LogInfo stopString
                                                        pure (zip (fmap fst inPairList) newThetaList, newStopCount)
                                                else -- )

                                                    errorWithoutStackTrace
                                                        ("Thompson search option " <> mFunction <> " not recognized " <> (show ["simple", "linear", "exponential"]))

-}

{- non monadic update theta -}
-- | updateTheta updates the expected success parameters for the bandit search list
updateTheta
    ∷ BanditType
    → Bool
    → Int
    → String
    → Int
    → [String]
    → [(String, Double)]
    → CPUTime
    → CPUTime
    → Int
    → Int
    → ([(String, Double)], Int)
updateTheta thisBandit thompsonSample mFactor mFunction counter infoStringList inPairList elapsedSeconds totalSeconds stopCount stopNum =
    if null inPairList
        then ([], stopCount)
        else
            let searchBandit =
                    if thisBandit == SearchBandit
                        then takeWhile (/= ',') (tail $ head infoStringList)
                        else -- GraphBandit

                            if "StaticApprox" `elem` (LS.splitOn "," $ head infoStringList)
                                then
                                    -- trace
                                    --    ("GraphBandit is StaticApprox")
                                        "StaticApproximation"
                                else
                                    if "MultiTraverse:False" `elem` (LS.splitOn "," $ head infoStringList)
                                        then
                                            --trace
                                            --    ("GraphBandit is SingleTraverse")
                                                "SingleTraverse"
                                        else
                                            -- trace
                                            --    ("GraphBandit is MultiTraverse ")
                                                "MultiTraverse"
                searchDeltaString = takeWhile (/= ',') $ tail $ dropWhile (/= ',') (tail $ head infoStringList)
                searchDelta = read searchDeltaString ∷ Double
            in  if not thompsonSample
                    then
                        let newStopCount =
                                if searchDelta <= 0.0
                                    then stopCount + 1
                                    else 0
                            stopString =
                                if newStopCount >= stopNum
                                    then ("\n\tSearch iterations have not improved for " <> (show newStopCount) <> " iterations--terminating this search command")
                                    else ""
                        in  if thisBandit == SearchBandit
                                then
                                    -- trace
                                    --  stopString
                                        (inPairList, newStopCount)
                                else (inPairList, newStopCount)
                    else -- trace ("UT1: " <> (show infoStringList)) (
                    -- update via results, previous history, memory \factor and type of memory "loss"

                        let -- get timing and benefit accounting for 0's
                            -- average time ratio in time factor for benefit adjustment
                            totalTime = (fromIntegral $ toSeconds totalSeconds) / (fromIntegral counter) ∷ Double
                            timeFactor =
                                if counter == 1
                                    then 1.0
                                    else
                                        if toSeconds elapsedSeconds == 0
                                            then 0.1
                                            else (fromIntegral $ toSeconds elapsedSeconds) / totalTime

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
                                    then ("\n\tSearch iterations have not improved for " <> (show newStopCount) <> " iterations--terminating this search command")
                                    else ""
                        in  -- check error
                            if isNothing indexBandit'
                                then error ("Bandit index not found: " <> searchBandit <> " in " <> (show inPairList) <> (show $ tail $ head infoStringList))
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

                                                newThetaList = fmap (/ totalTheta) newSuccessList
                                            in  -- trace ("UT2: " <> searchBandit <> " index " <> (show indexBandit) <> " total time: " <> (show $ toSeconds totalSeconds) <> " elapsed time: " <> (show $ toSeconds elapsedSeconds) <> " -> " <> (show (searchDelta, toSeconds elapsedSeconds, benefit)) <> "\n" <> (show $ fmap snd inPairList) <> "\n" <> (show newThetaList) <> "\n" <> (head infoStringList)) (
                                                --trace
                                                --    stopString
                                                    (zip (fmap fst inPairList) newThetaList, newStopCount)
                                        else -- )

                                        -- more complex 'recency' options

                                            if mFunction `elem` ["linear", "exponential"]
                                                then
                                                    let -- weight factors for previous theta (wN-1)* \theta
                                                        -- maxed ot counter so averaging not wierd for early iterations
                                                        mFactor' = min (fromIntegral counter ∷ Double) (fromIntegral mFactor ∷ Double)
                                                        (wN_1, wN) =
                                                            if mFunction == "linear"
                                                                then (mFactor' / (mFactor' + 1), 1.0 / (mFactor' + 1))
                                                                else
                                                                    if mFunction == "exponential"
                                                                        then (1.0 - (1.0 / (2.0 ** mFactor')), (1.0 / (2.0 ** mFactor')))
                                                                        else error ("Thompson search option " <> mFunction <> " not recognized " <> (show ["simple", "linear", "exponential"]))

                                                        -- simple "success-based" benefit, scaled to average time of search iteration
                                                        searchBenefit =
                                                            if searchDelta <= 0.0
                                                                then 0.0
                                                                else searchDelta / timeFactor

                                                        previousSuccessList = fmap (* (fromIntegral counter)) $ fmap snd inPairList

                                                        -- average of new am]nd previous if m=1, no memory if m=0, longer memory with larger m
                                                        -- m should be limited to counter
                                                        -- linear ((m * previous value) + new value) / m+1
                                                        -- exponential ((2^(-m)  * previous value) + new value) / (2^(-m) + 1)
                                                        newBanditVal =
                                                            if searchBenefit > inThetaBandit
                                                                then (wN * searchBenefit) + (wN_1 * inThetaBandit)
                                                                else (wN_1 * inThetaBandit) + (wN * (inThetaBandit + searchBenefit))

                                                        -- for nonBandit if not successful increment them (basically to downweight bandit), other wise not
                                                        incrementNonBanditVals =
                                                            if searchDelta <= 0.0
                                                                then 1.0 / ((fromIntegral $ length inPairList) - 1.0)
                                                                else 0.0
                                                        updatedSuccessList = fmap (+ incrementNonBanditVals) previousSuccessList

                                                        -- uopdate the bandit from list by splitting and rejoining, then normalizing to 1.0
                                                        firstBanditPart = take indexBandit updatedSuccessList
                                                        thirdBanditPart = drop (indexBandit + 1) updatedSuccessList
                                                        newSuccessList = firstBanditPart <> (newBanditVal : thirdBanditPart)
                                                        totalTheta = sum newSuccessList

                                                        newThetaList = fmap (/ totalTheta) newSuccessList
                                                    in  -- trace ("Update : \n" <> (show $ fmap snd inPairList) <> "\n" <> (show previousSuccessList) <> "\n" <> (show updatedSuccessList) <> "\n" <> (show newThetaList) <> "\n") $
                                                        -- trace ("Not simple: " <> mFunction <> " search benefit " <> (show searchBenefit) <> " " <> searchBandit <> " index " <> (show indexBandit) <> " total time: " <> (show $ toSeconds totalSeconds) <> " elapsed time: " <> (show $ toSeconds elapsedSeconds) <> " -> " <> (show (searchDelta, toSeconds elapsedSeconds)) <> "\n" <> (show $ fmap snd inPairList) <> "\n" <> (show newThetaList) <> "\n" <> (head infoStringList)) (
                                                        --trace
                                                        --    stopString
                                                            (zip (fmap fst inPairList) newThetaList, newStopCount)
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
    ∷ GlobalSettings
    → ProcessedData
    → [[VertexCost]]
    → Int
    → Bool
    → [(String, Double)]
    → Int
    → Int
    → CPUTime
    → ([ReducedPhylogeneticGraph], [String])
    → PhyG ([ReducedPhylogeneticGraph], [String])
performSearch inGS' inData' pairwiseDistances keepNum _ totalThetaList maxNetEdges rSeed inTime (inGraphList', _) =
    -- set up basic parameters for search/refine methods
    let thetaList = drop 3 totalThetaList
        numLeaves = V.length $ fst3 inData'
        -- set up log for sample
        thompsonString = "," <> (show totalThetaList)

        -- get infinite lists if integers and doubles
        -- TODO: This is problemetic since the source of randomness (rSeed) is used twice).
        randIntList = randomIntList rSeed
        randDoubleList = randoms (mkStdGen rSeed) ∷ [Double]
        randSeed0 = randIntList !! 0
        randSeed1 = randIntList !! 1

        -- this for constant access to random doubles need take for infinite list
        -- need to update as more random doubles are needed
        randDoubleVect = V.fromList $ take 20 randDoubleList

        -- choose search type from list with frequencies as input from searchForDuration
        -- adjust if all trees and networks chose to ensure some net stuff tried
        searchBandit =
            if null inGraphList'
                then "buildDistance"
                else
                    if (graphType inGS' == Tree)
                        then chooseElementAtRandomPair (randDoubleVect V.! 0) thetaList
                        else
                            let someGraphsNetwork = filter (== False) $ fmap LG.isTree $ fmap fst5 inGraphList'
                                tempBandit = chooseElementAtRandomPair (randDoubleVect V.! 0) thetaList
                            in  if (not . null) someGraphsNetwork
                                    then tempBandit
                                    else
                                        if tempBandit `notElem` ["networkMove", "networkDelete", "driftNetwork", "annealNetwork"]
                                            then tempBandit
                                            else chooseElementAtRandomPair (randDoubleVect V.! 0) [("networkAdd", 0.5), ("networkAddDelete", 0.5)]

        -- set graph valuation bandit
        graphEvaluationBandit = chooseElementAtRandomPair (randDoubleVect V.! 11) (take 3 totalThetaList)

        -- common build arguments including block and distance
        --    Tree does not use block--doesn't work very well for tree building
        buildMethod =
            if null inGraphList'
                then "unitary"
                else
                    if (graphType inGS' /= Tree)
                        then chooseElementAtRandomPair (randDoubleVect V.! 10) [("unitary", 0.8), ("block", 0.2)]
                        else "unitary"

        buildType =
            if searchBandit == "buildCharacter"
                then "character"
                else
                    if searchBandit == "buildDistance"
                        then "distance"
                        else chooseElementAtRandomPair (randDoubleVect V.! 12) [("distance", 0.5), ("character", 0.5)]

        numToCharBuild = fromInteger $ squareRoot $ toInteger numLeaves
        numToDistBuild = min 1000 (numLeaves * numLeaves)
        numDistToKeep = 50

        -- to resolve block build graphs
        reconciliationMethod = chooseElementAtRandomPair (randDoubleVect V.! 13) [("eun", 0.5), ("cun", 0.5)]

        wagnerOptions =
            if buildType == "distance"
                then
                    if buildMethod == "block"
                        then [("replicates", show numToCharBuild), ("rdwag", ""), ("best", show (1 ∷ Int))]
                        else [("replicates", show numToDistBuild), ("rdwag", ""), ("best", show numDistToKeep), ("return", show numToCharBuild)]
                else
                    if buildType == "character"
                        then
                            if buildMethod == "block"
                                then [("replicates", show (1 ∷ Int))]
                                else [("replicates", show numToCharBuild)]
                        else []

        showArg a = (fst a) <> ":" <> (snd a)

        blockOptions =
            if buildMethod == "block"
                then [("block", ""), ("atRandom", ""), ("displaytrees", show numToCharBuild), (reconciliationMethod, "")]
                else []

        -- common swap arguments
        swapKeep = min keepNum (chooseElementAtRandomPair (randDoubleVect V.! 14) [(1, 0.50), (2, 0.33), (4, 0.17)])

        -- common drift arguments
        maxChanges = chooseElementAtRandomPair (randDoubleVect V.! 1) [("5", 0.33), ("10", 0.34), ("20", 0.33)]
        acceptEqual = chooseElementAtRandomPair (randDoubleVect V.! 2) [("0.1", 0.5), ("0.5", 0.5)]
        acceptWorse = chooseElementAtRandomPair (randDoubleVect V.! 3) [("10.0", 0.33), ("20.0", 0.34), ("40", 0.33)]
        driftArgs = [("drift", ""), ("maxChanges", maxChanges), ("acceptEqual", acceptEqual), ("acceptWorse", acceptWorse)]

        -- common annealing arguments
        tempSteps = chooseElementAtRandomPair (randDoubleVect V.! 4) [("5", 0.33), ("10", 0.34), ("20", 0.33)]
        annealArgs = [("annealing", ""), ("steps", tempSteps)]

        -- common fuse options
        fusePairs = chooseElementAtRandomPair (randDoubleVect V.! 5) [("5", 0.45), ("10", 0.45), ("20", 0.1)]
        fuseKeep = 2 * keepNum

        -- network edit options
        netGeneralArgs = [("keep", show keepNum), ("steepest", ""), ("atRandom", ""), ("maxnetedges", show maxNetEdges)]
        netMoveArgs = ("netMove", "") : netGeneralArgs
        netAddArgs = ("netAdd", "") : netGeneralArgs
        netDelArgs = ("netDel", "") : netGeneralArgs
        netAddDelArgs = ("netAddDel", "") : netGeneralArgs
        netDriftAnnealMethod = chooseElementAtRandomPair (randDoubleVect V.! 17) [("netAdd", 0.5), ("netDel", 0.5)] -- ,("netMove", 0.25),  ("netAddDel", 0.25),]

        -- Genetic Algorithm Arguments
        -- stops after 2 rounds with no improvement (if generations > 2)
        popSize = chooseElementAtRandomPair (randDoubleVect V.! 6) [("10", 0.50), ("20", 0.25), ("40", 0.25)]
        generations = chooseElementAtRandomPair (randDoubleVect V.! 7) [("1", 1.0)] -- , "2" , "4"]
        severity = chooseElementAtRandomPair (randDoubleVect V.! 8) [("0.0", 0.33), ("1.0", 0.34), ("2.0", 0.33)]
        recombinations = chooseElementAtRandomPair (randDoubleVect V.! 9) [("10", 0.45), ("20", 0.45), ("40", 0.1)]

        gaArgs =
            [("popsize", popSize), ("generations", generations), ("severity", severity), ("recombinations", recombinations), ("stop", "2")]
                <> [("maxnetedges", show maxNetEdges)]

        -- unless fuse or genetic algorithm, only operate on "best" input graphs
        -- this to reduce memory footrpint when have multiple iterations
        inGraphList'' =
            if searchBandit `notElem` ["fuse", "fuseSPR", "fuseTBR", "geneticAlgorithm"]
                then GO.selectGraphs Best keepNum 0.0 (-1) inGraphList'
                else inGraphList'

        -- apply graph evaluation bandit
        transformToStaticApproximation =
            if U.getNumberNonExactCharacters (thd3 inData') == 0
                then False
                else
                    if graphEvaluationBandit == "StaticApproximation"
                        then True
                        else False

        transformMultiTraverse =
            if graphEvaluationBandit == "SingleTraverse"
                then True
                else False

    in do
            -- Can't do both static approx and multitraverse:False
            newDataMTF <- TRANS.transform [("multitraverse", "false")] inGS' inData' inData' 0 inGraphList''
            newDataSA  <- TRANS.transform [("staticapprox", [])] inGS' inData' inData' 0 inGraphList''
            let ((inGS, origData, inData, inGraphList), transformString) =
                    if transformToStaticApproximation && (useIA inGS') then (newDataSA, ",StaticApprox")
                    else
                        if transformMultiTraverse
                            then (newDataMTF, ",MultiTraverse:False")
                            else ((inGS', inData', inData', inGraphList''), "")
            -- bandit list with search arguments set
            -- primes (') for build to start with untransformed data
            (searchGraphs, searchArgs) <- case searchBandit of
                "buildCharacter" →
                    let buildArgs = [(buildType, "")] <> wagnerOptions <> blockOptions
                    in  --    graphList = B.buildGraph buildArgs inGS' inData' pairwiseDistances randSeed0
                        do
                            -- search
                            graphList ← B.buildGraph buildArgs inGS' inData' pairwiseDistances randSeed0
                            pure (graphList, buildArgs)
                "buildDistance" →
                    let -- build options
                        buildArgs = [(buildType, "")] <> wagnerOptions <> blockOptions
                        -- search for dist builds 1000, keeps 10 best distance then selects 10 best after rediagnosis
                        -- this line in here to allow for returning lots of rediagnosed distance trees, then
                        -- reducing to unique best cost trees--but is a memory pig
                    in  (\gList -> (gList, buildArgs)) <$> B.buildGraph buildArgs inGS' inData' pairwiseDistances randSeed0
                "buildSPR" →
                    let -- build part
                        buildArgs = [(buildType, "")] <> wagnerOptions <> blockOptions
                        buildGraphs = B.buildGraph buildArgs inGS' inData' pairwiseDistances randSeed0
                        -- swap options
                        swapType = "spr"
                        swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep), ("atrandom", "")]
                    in  -- search
                        do  buildGraphs' <- GO.selectGraphs Unique (maxBound ∷ Int) 0.0 (-1) <$> buildGraphs
                            swapList <- R.swapMaster swapArgs inGS inData randSeed1 buildGraphs'
                            pure (swapList, buildArgs <> swapArgs)
                "buildAlternate" →
                    let -- build part
                        buildArgs = [(buildType, "")] <> wagnerOptions <> blockOptions
                        buildGraphs = B.buildGraph buildArgs inGS' inData' pairwiseDistances randSeed0
                        -- swap options
                        swapType = "alternate" -- default anyway
                        swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep), ("atrandom", "")]
                    in  -- search
                        do  buildGraphs' <- GO.selectGraphs Unique (maxBound ∷ Int) 0.0 (-1) <$> buildGraphs
                            swapList <- R.swapMaster swapArgs inGS inData randSeed1 buildGraphs'
                            pure (swapList, buildArgs <> swapArgs)
                "swapSPR" →
                    let -- swap options
                        swapType = "spr"
                        swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep), ("atrandom", "")]
                    in  -- search
                        do swapList <- R.swapMaster swapArgs inGS inData randSeed1 inGraphList
                           pure (swapList, swapArgs)
                "swapAlternate" →
                    let -- swap options
                        swapType = "alternate" -- default anyway
                        swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep), ("atrandom", "")]
                    in  -- search
                        do swapList <- R.swapMaster swapArgs inGS inData randSeed1 inGraphList
                           pure (swapList, swapArgs)
                -- drift only best graphs
                "driftSPR" →
                    let -- swap args
                        swapType = "spr"
                        swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep), ("atrandom", "")]
                        -- swap with drift (common) arguments
                        swapDriftArgs = swapArgs <> driftArgs
                    in  -- perform search
                        do swapList <- R.swapMaster swapDriftArgs inGS inData randSeed1 inGraphList
                           pure (swapList, swapArgs)
                -- drift only best graphs
                "driftAlternate" →
                    let -- swap args
                        swapType = "alternate"
                        swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep), ("atrandom", "")]
                        -- swap with drift (common) arguments
                        swapDriftArgs = swapArgs <> driftArgs
                    in  -- perform search
                        do swapList <- R.swapMaster swapDriftArgs inGS inData randSeed1 inGraphList
                           pure (swapList, swapDriftArgs)
                -- anneal only best graphs
                "annealSPR" →
                    let -- swap args
                        swapType = "spr"
                        swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep), ("atrandom", "")]
                        -- swap with anneal (common) arguments
                        swapAnnealArgs = swapArgs <> annealArgs
                    in  -- perform search
                        do swapList <- R.swapMaster swapAnnealArgs inGS inData randSeed1 inGraphList
                           pure (swapList, swapAnnealArgs)
                -- anneal only best graphs
                "annealAlternate" →
                    let -- swap args
                        swapType = "alternate"
                        swapArgs = [(swapType, ""), ("steepest", ""), ("keep", show swapKeep), ("atrandom", "")]
                        -- swap with anneal (common) arguments
                        swapAnnealArgs = swapArgs <> annealArgs
                    in  -- perform search
                        do swapList <- R.swapMaster swapAnnealArgs inGS inData randSeed1 inGraphList
                           pure (swapList, swapAnnealArgs)
                "geneticAlgorithm" →
                    do
                    -- args from above
                    -- perform search
                        gaReturn <- R.geneticAlgorithmMaster gaArgs inGS inData randSeed1 inGraphList
                        pure (gaReturn, gaArgs)
                "fuse" →
                    -- should more graphs be added if only one?  Would downweight fuse perhpas too much
                    let -- fuse arguments
                        -- this to limit memory footprint of fuse during search
                        --gsNum = min (graphsSteepest inGS) 5
                        --inGSgs1 = inGS{graphsSteepest = gsNum}
                        fuseArgs =
                            [("none", ""), ("all", ""), ("unique", ""), ("atrandom", ""), ("pairs", fusePairs), ("keep", show fuseKeep), ("noreciprocal", "")]
                    in  -- perform search
                        R.fuseGraphs fuseArgs inGS inData randSeed1 inGraphList <&> (\x -> (x, fuseArgs))
                "fuseSPR" →
                    let -- fuse arguments
                        --inGSgs1 = inGS{graphsSteepest = 1}
                        fuseArgs =
                            [("spr", ""), ("all", ""), ("unique", ""), ("atrandom", ""), ("pairs", fusePairs), ("keep", show fuseKeep), ("noreciprocal", "")]
                    in  -- perform search
                        R.fuseGraphs fuseArgs inGS inData randSeed1 inGraphList <&> (\x -> (x, fuseArgs))
                "fuseTBR" →
                    let -- fuse arguments
                        --inGSgs1 = inGS{graphsSteepest = 1}
                        fuseArgs =
                            [("tbr", ""), ("all", ""), ("unique", ""), ("atrandom", ""), ("pairs", fusePairs), ("keep", show fuseKeep), ("noreciprocal", "")]
                    in  -- perform search
                        R.fuseGraphs fuseArgs inGS inData randSeed1 inGraphList <&> (\x -> (x, fuseArgs))
                "networkAdd" →
                    let -- network add args
                        netEditArgs = netAddArgs
                    in  -- perform search
                        R.netEdgeMaster netEditArgs inGS inData randSeed1 inGraphList <&> (\x -> (x, netEditArgs))
                "networkDelete" →
                    let -- network delete args
                        netEditArgs = netDelArgs
                    in  -- perform search
                        R.netEdgeMaster netEditArgs inGS inData randSeed1 inGraphList <&> (\x -> (x, netEditArgs))
                "networkAddDelete" →
                    let -- network add/delete args
                        netEditArgs = netAddDelArgs
                    in  -- perform search
                        R.netEdgeMaster netEditArgs inGS inData randSeed1 inGraphList <&> (\x -> (x, netEditArgs))
                "networkMove" →
                    let -- network move args
                        netEditArgs = netMoveArgs
                    in  -- perform search
                        R.netEdgeMaster netEditArgs inGS inData randSeed1 inGraphList <&> (\x -> (x, netEditArgs))
                "driftNetwork" →
                    let -- network add/delete  + drift args
                        netEditArgs = [(netDriftAnnealMethod, "")] <> netGeneralArgs <> driftArgs
                    in  -- perform search
                        R.netEdgeMaster netEditArgs inGS inData randSeed1 inGraphList <&> (\x -> (x, netEditArgs))
                "annealNetwork" →
                    let -- network add/delete  + annealing  args
                        netEditArgs = [(netDriftAnnealMethod, "")] <> netGeneralArgs <> annealArgs
                    in  -- perform search
                        R.netEdgeMaster netEditArgs inGS inData randSeed1 inGraphList <&> (\x -> (x, netEditArgs))
                _ → error ("Unknown/unimplemented method in search: " <> searchBandit)

            -- process
            let uniqueGraphs' = take keepNum $ GO.selectGraphs Unique (maxBound ∷ Int) 0.0 (-1) (searchGraphs <> inGraphList)
            newDataMT <- TRANS.transform [("multiTraverse", "true")] inGS origData inData 0 uniqueGraphs'
            newDataT  <- TRANS.transform [("dynamic", [])] inGS' origData inData 0 uniqueGraphs'
            let (uniqueGraphs, transString)
                    | (not transformToStaticApproximation && not transformMultiTraverse) = (uniqueGraphs', "")
                    | transformToStaticApproximation = (fth4 newDataT, ",Dynamic")
                    | otherwise = (fth4  newDataMT, ",MultiTraverse:True")

            -- string of delta and cost of graphs
            let deltaString
                    | null inGraphList' = "10.0,"
                    | otherwise =  show ((minimum $ fmap snd5 inGraphList') - (minimum $ fmap snd5 uniqueGraphs)) -- <> ","
                    
            let currentBestString
                    | not $ null uniqueGraphs = show $ minimum $ fmap snd5 uniqueGraphs
                    | otherwise =  show infinity

            -- create string for search stats
            let searchString = ","
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



-- | getSearchParams takes arguments and returns search params
getSearchParams ∷ [Argument] → PhyG (Int, Int, Int, Bool, Int, String, Int, Int)
getSearchParams inArgs =
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
                        | null keepList = Just 10
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
                                                                                                when (mLinear && mExponential) $ logWith LogWarn ("Thompson recency function specification has both 'linear' and 'exponential', defaulting to 'linear'\n")
                                                                                                pure $ ( searchTime
                                                                                                        , fromJust keepNum
                                                                                                        , fromJust instances
                                                                                                        , thompson
                                                                                                        , fromJust mFactor
                                                                                                        , mFunction
                                                                                                        , fromJust maxNetEdges
                                                                                                        , fromJust stopNum
                                                                                                        )
