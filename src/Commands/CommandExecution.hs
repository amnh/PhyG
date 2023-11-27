{-# LANGUAGE CPP #-}

{- |
Module to coordinate command execution.
-}
module Commands.CommandExecution (
    executeCommands,
    executeRenameReblockCommands,
    getDataListList,
) where

import Commands.CommandUtilities
import Commands.Transform qualified as TRANS
import Commands.Verify qualified as VER
import Control.Arrow ((&&&))
import Control.Monad (when)
import Control.Monad.IO.Class (MonadIO (..))
import Data.CSV qualified as CSV
import Data.Char
import Data.Char qualified as C
import Data.InfList qualified as IL
import Data.List qualified as L
import Data.List.Split qualified as SL
import Data.Maybe
import Data.Ord
import Data.Text.Lazy qualified as T
import Data.Vector qualified as V
import Data.Version qualified as DV
import GeneralUtilities
import GraphOptimization.Traversals qualified as TRAV
import Graphs.GraphOperations qualified as GO
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import PHANE.Evaluation.Verbosity (Verbosity (..))
import Reconciliation.ReconcileGraphs qualified as R
import Search.Build qualified as B
import Search.Refinement qualified as REF
import Search.Search qualified as S
import Support.Support qualified as SUP
import System.CPU qualified as SC
import System.IO
import System.IO.Unsafe qualified as SIOU
import System.Info qualified as SI
import System.Timing
import Text.Read
import Types.Types
import Utilities.Distances qualified as D
import Utilities.Utilities qualified as U


{- | executeCommands reads input files and returns raw data
need to close files after read
-}
executeCommands
    ∷ GlobalSettings
    → ([NameText], [(NameText, NameText)])
    → Int
    → String
    → ProcessedData
    → ProcessedData
    → ProcessedData
    → [ReducedPhylogeneticGraph]
    → [[VertexCost]]
    → [Int]
    → [ReducedPhylogeneticGraph]
    → [Command]
    → PhyG ([ReducedPhylogeneticGraph], GlobalSettings, [Int], [ReducedPhylogeneticGraph])
executeCommands globalSettings excludeRename numInputFiles crossReferenceString origProcessedData processedData reportingData curGraphs pairwiseDist seedList supportGraphList commandList = do
    if null commandList
        then pure (curGraphs, globalSettings, seedList, supportGraphList)
        else do
            let (firstOption, firstArgs) = head commandList

            -- skip "Read" and "Rename "commands already processed
            if firstOption == Read
                then error ("Read command should already have been processed: " <> show (firstOption, firstArgs))
                else
                    if firstOption == Rename
                        then error ("Rename command should already have been processed: " <> show (firstOption, firstArgs))
                        else
                            if firstOption == Reblock
                                then error ("Reblock command should already have been processed: " <> show (firstOption, firstArgs))
                                else
                                    if firstOption == Run
                                        then error ("Run command should already have been processed: " <> show (firstOption, firstArgs))
                                        else -- other commands

                                            if firstOption == Build
                                                then do
                                                    (elapsedSeconds, newGraphList') ←
                                                        timeOp $ pure $ B.buildGraph firstArgs globalSettings processedData pairwiseDist (head seedList)

                                                    newGraphList ← newGraphList'
                                                    let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
                                                    let newSearchData = searchInfo : searchData globalSettings

                                                    executeCommands
                                                        (globalSettings{searchData = newSearchData})
                                                        excludeRename
                                                        numInputFiles
                                                        crossReferenceString
                                                        origProcessedData
                                                        processedData
                                                        reportingData
                                                        (curGraphs <> newGraphList)
                                                        pairwiseDist
                                                        (tail seedList)
                                                        supportGraphList
                                                        (tail commandList)
                                                else
                                                    if firstOption == Refine
                                                        then do
                                                            (elapsedSeconds, newGraphList') ←
                                                                timeOp $ pure $ REF.refineGraph firstArgs globalSettings processedData (head seedList) curGraphs

                                                            newGraphList ← newGraphList'
                                                            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
                                                            let newSearchData = searchInfo : searchData globalSettings

                                                            executeCommands
                                                                (globalSettings{searchData = newSearchData})
                                                                excludeRename
                                                                numInputFiles
                                                                crossReferenceString
                                                                origProcessedData
                                                                processedData
                                                                reportingData
                                                                newGraphList
                                                                pairwiseDist
                                                                (tail seedList)
                                                                supportGraphList
                                                                (tail commandList)
                                                        else
                                                            if firstOption == Fuse
                                                                then do
                                                                    (elapsedSeconds, newGraphList) ←
                                                                        timeOp $ REF.fuseGraphs firstArgs globalSettings processedData (head seedList) curGraphs

                                                                    let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
                                                                    let newSearchData = searchInfo : searchData globalSettings

                                                                    executeCommands
                                                                        (globalSettings{searchData = newSearchData})
                                                                        excludeRename
                                                                        numInputFiles
                                                                        crossReferenceString
                                                                        origProcessedData
                                                                        processedData
                                                                        reportingData
                                                                        newGraphList
                                                                        pairwiseDist
                                                                        (tail seedList)
                                                                        supportGraphList
                                                                        (tail commandList)
                                                                else
                                                                    if firstOption == Report
                                                                        then do
                                                                            let doDotPDF = elem "dotpdf" $ fmap (fmap toLower . fst) firstArgs
                                                                            let collapse' = elem "collapse" $ fmap (fmap toLower . fst) firstArgs
                                                                            let noCollapse' = elem "nocollapse" $ fmap (fmap toLower . fst) firstArgs
                                                                            let reconcile = any ((== "reconcile") . fst) firstArgs

                                                                            -- set default collapse for dotPDF to True, False otherwise
                                                                            let collapse -- this will casue problems with reconcile
                                                                                    | reconcile = False
                                                                                    | collapse' = True
                                                                                    | noCollapse' = False
                                                                                    -- \| doDotPDF = True
                                                                                    | otherwise = False

                                                                            let curGraphs' =
                                                                                    if not collapse
                                                                                        then curGraphs
                                                                                        else fmap U.collapseReducedGraph curGraphs

                                                                            -- use 'temp' updated graphs s don't repeatedly add model and root complexityies
                                                                            -- reporting collapsed
                                                                            -- reverse sorting graphs by cost
                                                                            let rediagnoseWithReportingData = optimalityCriterion globalSettings == NCM && U.has4864PackedChars (thd3 processedData)
                                                                            updatedCostGraphs ←
                                                                                TRAV.updateGraphCostsComplexities globalSettings reportingData processedData rediagnoseWithReportingData curGraphs'
                                                                            let graphsWithUpdatedCosts =
                                                                                    L.sortOn
                                                                                        (Data.Ord.Down . snd5)
                                                                                        updatedCostGraphs
                                                                            -- (TRAV.updateGraphCostsComplexities globalSettings reportingData processedData rediagnoseWithReportingData curGraphs')
                                                                            reportStuff@(reportString, outFile, writeMode) ←
                                                                                reportCommand
                                                                                    globalSettings
                                                                                    firstArgs
                                                                                    excludeRename
                                                                                    numInputFiles
                                                                                    crossReferenceString
                                                                                    reportingData
                                                                                    graphsWithUpdatedCosts
                                                                                    supportGraphList
                                                                                    pairwiseDist

                                                                            if null reportString
                                                                                then do
                                                                                    executeCommands
                                                                                        globalSettings
                                                                                        excludeRename
                                                                                        numInputFiles
                                                                                        crossReferenceString
                                                                                        origProcessedData
                                                                                        processedData
                                                                                        reportingData
                                                                                        curGraphs
                                                                                        pairwiseDist
                                                                                        seedList
                                                                                        supportGraphList
                                                                                        (tail commandList)
                                                                                else do
                                                                                    logWith LogInfo ("Report writing to " <> outFile <> "\n")

                                                                                    if doDotPDF
                                                                                        then do
                                                                                            let reportString' = changeDotPreamble "digraph {" "digraph G {\n\trankdir = LR;\tnode [ shape = none];\n" reportString
                                                                                            printGraphVizDot reportString' outFile
                                                                                            executeCommands
                                                                                                globalSettings
                                                                                                excludeRename
                                                                                                numInputFiles
                                                                                                crossReferenceString
                                                                                                origProcessedData
                                                                                                processedData
                                                                                                reportingData
                                                                                                curGraphs
                                                                                                pairwiseDist
                                                                                                seedList
                                                                                                supportGraphList
                                                                                                (tail commandList)
                                                                                        else do
                                                                                            if outFile == "stderr"
                                                                                                then liftIO $ hPutStr stderr reportString
                                                                                                else
                                                                                                    if outFile == "stdout"
                                                                                                        then liftIO $ putStr reportString
                                                                                                        else
                                                                                                            if writeMode == "overwrite"
                                                                                                                then liftIO $ writeFile outFile reportString
                                                                                                                else
                                                                                                                    if writeMode == "append"
                                                                                                                        then liftIO $ appendFile outFile reportString
                                                                                                                        else failWithPhase Parsing ("Error 'read' command not properly formatted" <> show reportStuff)
                                                                                            executeCommands
                                                                                                globalSettings
                                                                                                excludeRename
                                                                                                numInputFiles
                                                                                                crossReferenceString
                                                                                                origProcessedData
                                                                                                processedData
                                                                                                reportingData
                                                                                                curGraphs
                                                                                                pairwiseDist
                                                                                                seedList
                                                                                                supportGraphList
                                                                                                (tail commandList)
                                                                        else
                                                                            if firstOption == Search
                                                                                then do
                                                                                    (elapsedSeconds, output) ←
                                                                                        timeOp $ S.search firstArgs globalSettings processedData pairwiseDist (head seedList) curGraphs
                                                                                    -- in pure result
                                                                                    -- (newGraphList, serchInfoList) <- S.search firstArgs globalSettings origProcessedData processedData reportingDatapairwiseDist (head seedList) curGraphs
                                                                                    let searchInfo =
                                                                                            makeSearchRecord
                                                                                                firstOption
                                                                                                firstArgs
                                                                                                curGraphs
                                                                                                (fst output)
                                                                                                (fromIntegral $ toMilliseconds elapsedSeconds)
                                                                                                (concatMap (L.intercalate "\n") (snd output))
                                                                                    let newSearchData = searchInfo : searchData globalSettings
                                                                                    executeCommands
                                                                                        (globalSettings{searchData = newSearchData})
                                                                                        excludeRename
                                                                                        numInputFiles
                                                                                        crossReferenceString
                                                                                        origProcessedData
                                                                                        processedData
                                                                                        reportingData
                                                                                        (fst output)
                                                                                        pairwiseDist
                                                                                        (tail seedList)
                                                                                        supportGraphList
                                                                                        (tail commandList)
                                                                                else
                                                                                    if firstOption == Select
                                                                                        then do
                                                                                            (elapsedSeconds, newGraphList) ← timeOp $ pure $ GO.selectPhylogeneticGraphReduced firstArgs (head seedList) curGraphs

                                                                                            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
                                                                                            let newSearchData = searchInfo : searchData globalSettings

                                                                                            let typeSelected =
                                                                                                    if null firstArgs
                                                                                                        then "best"
                                                                                                        else fmap C.toLower $ fst $ head firstArgs

                                                                                            logWith LogInfo ("Selecting " <> typeSelected <> " graphs" <> "\n")
                                                                                            executeCommands
                                                                                                (globalSettings{searchData = newSearchData})
                                                                                                excludeRename
                                                                                                numInputFiles
                                                                                                crossReferenceString
                                                                                                origProcessedData
                                                                                                processedData
                                                                                                reportingData
                                                                                                newGraphList
                                                                                                pairwiseDist
                                                                                                (tail seedList)
                                                                                                supportGraphList
                                                                                                (tail commandList)
                                                                                        else
                                                                                            if firstOption == Set
                                                                                                then do
                                                                                                    -- if set changes graph aspects--may nned to reoptimize
                                                                                                    (newGlobalSettings, newProcessedData, seedList') ← setCommand firstArgs globalSettings origProcessedData processedData seedList
                                                                                                    newGraphList ←
                                                                                                        if not (requireReoptimization globalSettings newGlobalSettings)
                                                                                                            then pure curGraphs
                                                                                                            else -- TODO should be parallel
                                                                                                                mapM (TRAV.multiTraverseFullyLabelGraphReduced newGlobalSettings newProcessedData True True Nothing) (fmap fst5 curGraphs)

                                                                                                    let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList 0 "No Comment"
                                                                                                    let newSearchData = searchInfo : searchData newGlobalSettings

                                                                                                    if not (requireReoptimization globalSettings newGlobalSettings)
                                                                                                        then do logWith LogInfo "No need to reoptimize graphs\n"
                                                                                                        else do logWith LogInfo "Reoptimizing gaphs\n"

                                                                                                    executeCommands
                                                                                                        (newGlobalSettings{searchData = newSearchData})
                                                                                                        excludeRename
                                                                                                        numInputFiles
                                                                                                        crossReferenceString
                                                                                                        origProcessedData
                                                                                                        processedData
                                                                                                        reportingData
                                                                                                        newGraphList
                                                                                                        pairwiseDist
                                                                                                        seedList'
                                                                                                        supportGraphList
                                                                                                        (tail commandList)
                                                                                                else
                                                                                                    if firstOption == Swap
                                                                                                        then do
                                                                                                            (elapsedSeconds, newGraphList) ←
                                                                                                                timeOp $ REF.swapMaster firstArgs globalSettings processedData (head seedList) curGraphs

                                                                                                            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
                                                                                                            let newSearchData = searchInfo : searchData globalSettings

                                                                                                            executeCommands
                                                                                                                (globalSettings{searchData = newSearchData})
                                                                                                                excludeRename
                                                                                                                numInputFiles
                                                                                                                crossReferenceString
                                                                                                                origProcessedData
                                                                                                                processedData
                                                                                                                reportingData
                                                                                                                newGraphList
                                                                                                                pairwiseDist
                                                                                                                (tail seedList)
                                                                                                                supportGraphList
                                                                                                                (tail commandList)
                                                                                                        else
                                                                                                            if firstOption == Support
                                                                                                                then do
                                                                                                                    (elapsedSeconds, newSupportGraphList') ←
                                                                                                                        timeOp $ pure $ SUP.supportGraph firstArgs globalSettings processedData (head seedList) curGraphs

                                                                                                                    newSupportGraphList ← newSupportGraphList'
                                                                                                                    let searchInfo =
                                                                                                                            makeSearchRecord firstOption firstArgs curGraphs newSupportGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
                                                                                                                    let newSearchData = searchInfo : searchData globalSettings

                                                                                                                    executeCommands
                                                                                                                        (globalSettings{searchData = newSearchData})
                                                                                                                        excludeRename
                                                                                                                        numInputFiles
                                                                                                                        crossReferenceString
                                                                                                                        origProcessedData
                                                                                                                        processedData
                                                                                                                        reportingData
                                                                                                                        curGraphs
                                                                                                                        pairwiseDist
                                                                                                                        (tail seedList)
                                                                                                                        (supportGraphList <> newSupportGraphList)
                                                                                                                        (tail commandList)
                                                                                                                else
                                                                                                                    if firstOption == Transform
                                                                                                                        then do
                                                                                                                            (elapsedSeconds, (newGS, newOrigData, newProcessedData, newGraphs)) ←
                                                                                                                                timeOp $ TRANS.transform firstArgs globalSettings origProcessedData processedData (head seedList) curGraphs

                                                                                                                            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphs (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
                                                                                                                            let newSearchData = searchInfo : searchData globalSettings

                                                                                                                            executeCommands
                                                                                                                                (newGS{searchData = newSearchData})
                                                                                                                                excludeRename
                                                                                                                                numInputFiles
                                                                                                                                crossReferenceString
                                                                                                                                newOrigData
                                                                                                                                newProcessedData
                                                                                                                                reportingData
                                                                                                                                newGraphs
                                                                                                                                pairwiseDist
                                                                                                                                (tail seedList)
                                                                                                                                supportGraphList
                                                                                                                                (tail commandList)
                                                                                                                        else error ("Command " <> show firstOption <> " not recognized/implemented")


-- | makeSearchRecord take sbefore and after data of a commend and returns SearchData record
makeSearchRecord
    ∷ Instruction → [Argument] → [ReducedPhylogeneticGraph] → [ReducedPhylogeneticGraph] → Int → String → SearchData
makeSearchRecord firstOption firstArgs curGraphs newGraphList elapsedTime comment =
    SearchData
        { instruction = firstOption
        , arguments = firstArgs
        , minGraphCostIn =
            if null curGraphs
                then infinity
                else minimum $ fmap snd5 curGraphs
        , maxGraphCostIn =
            if null curGraphs
                then infinity
                else maximum $ fmap snd5 curGraphs
        , numGraphsIn = length curGraphs
        , minGraphCostOut =
            if null newGraphList
                then infinity
                else minimum $ fmap snd5 newGraphList
        , maxGraphCostOut =
            if null newGraphList
                then infinity
                else maximum $ fmap snd5 newGraphList
        , numGraphsOut = length newGraphList
        , commentString = comment
        , duration = elapsedTime
        }


{- | setCommand takes arguments to change globalSettings and multiple data aspects (e.g. 'blocks')
needs to be abtracted--too long
if seed list is empty [] then processes first set--confusing--shold be refactored
-}
setCommand
    ∷ [Argument] → GlobalSettings → ProcessedData → ProcessedData → [Int] → PhyG (GlobalSettings, ProcessedData, [Int])
setCommand argList globalSettings origProcessedData processedData inSeedList =
    let processor f = fmap (fmap C.toLower) . filter (/= "") . fmap f
        commandList = processor fst argList
        optionList = processor snd argList
        checkCommandList = checkCommandArgs "set" commandList VER.setArgList
        leafNameVect = fst3 processedData
    in  do
            when (not checkCommandList) . failWithPhase Parsing $
                "Unrecognized command in 'set': " <> show argList

            when (length commandList > 1 || length optionList > 1) . failWithPhase Parsing $
                "Set option error: can only have one set argument for each command: " <> show (commandList, optionList)

            -- early extraction of partition character and bc2-gt64 follows from null inputs
            -- this due to not having all info required for all global settings,
            -- so options restricted and repeated
            -- needs to be fixed to be more clear and clean
            case inSeedList of
                [] → case head commandList of
                    "partitioncharacter" → case head optionList of
                        localPartitionChar@[x] →
                            do
                                logWith LogInfo $ "PartitionCharacter set to '" <> localPartitionChar <> "'\n"
                                let x = globalSettings{partitionCharacter = localPartitionChar}
                                pure (x, processedData, inSeedList)
                        val →
                            failWithPhase Parsing $
                                "Error in 'set' command. Partitioncharacter '" <> val <> "' must be a single character\n"
                    {--}

                    "missingthreshold" → case readMaybe (head optionList) ∷ Maybe Int of
                        Nothing →
                            failWithPhase Parsing $
                                "Set option 'missingThreshold' must be set to an integer value (e.g. missingThreshold:50): " <> head optionList
                        Just val
                            | val < 0 || 100 < val →
                                failWithPhase Parsing $
                                    "Set option 'missingThreshold' must be set to an integer value between 0 and 100: " <> head optionList
                        Just val → do
                            logWith LogInfo $ "MissingThreshold set to " <> head optionList
                            pure (globalSettings{missingThreshold = val}, processedData, inSeedList)

                    -- sets root cost as well-- need in both places--one to process data and one to
                    -- keep in current global
                    "criterion" → do
                        localCriterion ← case head optionList of
                            "parsimony" → pure Parsimony
                            "pmdl" → pure PMDL
                            "si" → pure SI
                            "mapa" → pure MAPA
                            "ncm" → pure NCM
                            val →
                                failWithPhase Parsing $
                                    "Error in 'set' command. Criterion '" <> val <> "' is not 'parsimony', 'ml', or 'pmdl'"

                        -- create lazy list of graph complexity indexed by number of network nodes--need leaf number for base tree complexity
                        (lGraphComplexityList, lRootComplexity) ← case localCriterion of
                            NCM
                                | origProcessedData /= emptyProcessedData →
                                    pure
                                        (IL.repeat (0.0, 0.0), U.calculateNCMRootCost origProcessedData)
                            NCM →
                                pure
                                    (IL.repeat (0.0, 0.0), U.calculateNCMRootCost processedData)
                            Parsimony → pure $ (IL.repeat (0.0, 0.0), 0.0)
                            val
                                | val `elem` [PMDL, SI, MAPA] →
                                    pure $
                                        (U.calculateGraphComplexity &&& U.calculateW15RootCost) processedData
                            val → failWithPhase Parsing $ "Optimality criterion not recognized: " <> show val

                        let lGraphFactor
                                | localCriterion `elem` [PMDL, SI, MAPA] = PMDLGraph
                                | otherwise = graphFactor globalSettings

                        logWith LogInfo $ case localCriterion of
                            NCM → unwords ["Optimality criterion set to", show NCM, "in -log (base 10) likelihood units"]
                            val
                                | val `elem` [PMDL, SI] →
                                    unwords ["Optimality criterion set to", show val, "Tree Complexity =", show . fst $ IL.head lGraphComplexityList, "bits"]
                            val → "Optimality criterion set to " <> show val

                        pure $
                            ( globalSettings
                                { optimalityCriterion = localCriterion
                                , graphComplexityList = lGraphComplexityList
                                , rootComplexity = lRootComplexity
                                , graphFactor = lGraphFactor
                                }
                            , processedData
                            , inSeedList
                            )
                    "bc2" →
                        let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                            changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            changeMaybe = readMaybe changeString ∷ Maybe Double
                        in  do
                                when (length commandList /= length optionList) . failWithPhase Parsing $
                                    "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                when (null $ head optionList) . failWithPhase Parsing $
                                    "Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): no values found "
                                when (',' `notElem` head optionList) . failWithPhase Parsing $
                                    "Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): no comma found "
                                case liftA2 (,) noChangeMaybe changeMaybe of
                                    Nothing →
                                        failWithPhase Parsing $
                                            "Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): "
                                                <> head optionList
                                    Just val@(noChangeValue, changeValue) → do
                                        when (bc2 globalSettings /= val) . logWith LogInfo $
                                            "bit cost 2 state set to " <> show val
                                        pure (globalSettings{bc2 = val}, processedData, inSeedList)
                    "bc4" →
                        let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                            changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            changeMaybe = readMaybe changeString ∷ Maybe Double
                        in  do
                                when (length commandList /= length optionList) . failWithPhase Parsing $
                                    "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                when (null $ head optionList) . failWithPhase Parsing $
                                    "Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): no values found "
                                when (',' `notElem` head optionList) . failWithPhase Parsing $
                                    "Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): no comma found "
                                case liftA2 (,) noChangeMaybe changeMaybe of
                                    Nothing →
                                        failWithPhase Parsing $
                                            "Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): "
                                                <> head optionList
                                    Just val@(noChangeValue, changeValue) → do
                                        when (bc4 globalSettings /= val) . logWith LogInfo $
                                            "bit cost 4 state set to " <> show val
                                        pure (globalSettings{bc4 = val}, processedData, inSeedList)
                    "bc5" →
                        let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                            changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            changeMaybe = readMaybe changeString ∷ Maybe Double
                        in  do
                                when (length commandList /= length optionList) . failWithPhase Parsing $
                                    "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                when (null $ head optionList) . failWithPhase Parsing $
                                    "Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): no values found "
                                when (',' `notElem` head optionList) . failWithPhase Parsing $
                                    "Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): no comma found "
                                case liftA2 (,) noChangeMaybe changeMaybe of
                                    Nothing →
                                        failWithPhase Parsing $
                                            "Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): "
                                                <> head optionList
                                    Just val@(noChangeValue, changeValue) → do
                                        when (bc5 globalSettings /= val) . logWith LogInfo $
                                            "bit cost 5 state set to " <> show val
                                        pure (globalSettings{bc5 = val}, processedData, inSeedList)
                    "bc8" →
                        let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                            changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            changeMaybe = readMaybe changeString ∷ Maybe Double
                        in  do
                                when (length commandList /= length optionList) . failWithPhase Parsing $
                                    "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                when (null $ head optionList) . failWithPhase Parsing $
                                    "Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): no values found "
                                when (',' `notElem` head optionList) . failWithPhase Parsing $
                                    "Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): no comma found "
                                case liftA2 (,) noChangeMaybe changeMaybe of
                                    Nothing →
                                        failWithPhase Parsing $
                                            "Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): "
                                                <> head optionList
                                    Just val@(noChangeValue, changeValue) → do
                                        when (bc8 globalSettings /= val) . logWith LogInfo $
                                            "bit cost 8 state set to " <> show val
                                        pure (globalSettings{bc8 = val}, processedData, inSeedList)
                    "bc64" →
                        let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                            changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            changeMaybe = readMaybe changeString ∷ Maybe Double
                        in  do
                                when (length commandList /= length optionList) . failWithPhase Parsing $
                                    "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                when (null $ head optionList) . failWithPhase Parsing $
                                    "Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): no values found "
                                when (',' `notElem` head optionList) . failWithPhase Parsing $
                                    "Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): no comma found "
                                case liftA2 (,) noChangeMaybe changeMaybe of
                                    Nothing →
                                        failWithPhase Parsing $
                                            "Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): "
                                                <> head optionList
                                    Just val@(noChangeValue, changeValue) → do
                                        when (bc64 globalSettings /= val) . logWith LogInfo $
                                            "bit cost 64 state set to " <> show val
                                        pure (globalSettings{bc64 = val}, processedData, inSeedList)
                    "bcgt64" →
                        let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                            changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            changeMaybe = readMaybe changeString ∷ Maybe Double
                        in  do
                                when (length commandList /= length optionList) . failWithPhase Parsing $
                                    "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                when (null $ head optionList) . failWithPhase Parsing $
                                    "Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): no values found "
                                when (',' `notElem` head optionList) . failWithPhase Parsing $
                                    "Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): no comma found "
                                case liftA2 (,) noChangeMaybe changeMaybe of
                                    Nothing →
                                        failWithPhase Parsing $
                                            "Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): "
                                                <> head optionList
                                    Just val@(noChangeValue, changeValue) → do
                                        when (bcgt64 globalSettings /= val) . logWith LogInfo $
                                            "bit cost > 64 state set to " <> show val
                                        pure (globalSettings{bcgt64 = val}, processedData, inSeedList)

                    -- partition character to reset
                    "bcgt64" → pure (globalSettings, processedData, inSeedList)
                -- =-=-=-=-=-=-=-=-=-=-=-=-=
                -- =                       =
                -- = regular command stuff =
                -- = not initial at start  =
                -- =                       =
                -- =-=-=-=-=-=-=-=-=-=-=-=-=
                seed : otherSeeds → case head commandList of
                    "bc2" →
                        let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                            changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            changeMaybe = readMaybe changeString ∷ Maybe Double
                        in  do
                                when (length commandList /= length optionList) . failWithPhase Parsing $
                                    "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                when (null $ head optionList) . failWithPhase Parsing $
                                    "Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): no values found "
                                when (',' `notElem` head optionList) . failWithPhase Parsing $
                                    "Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): no comma found "
                                case liftA2 (,) noChangeMaybe changeMaybe of
                                    Nothing →
                                        failWithPhase Parsing $
                                            "Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): "
                                                <> head optionList
                                    Just val@(noChangeValue, changeValue) → do
                                        when (bc2 globalSettings /= val) . logWith LogInfo $
                                            "bit cost 2 state set to " <> show val
                                        pure (globalSettings{bc2 = val}, processedData, inSeedList)
                    "bc4" →
                        let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                            changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            changeMaybe = readMaybe changeString ∷ Maybe Double
                        in  do
                                when (length commandList /= length optionList) . failWithPhase Parsing $
                                    "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                when (null $ head optionList) . failWithPhase Parsing $
                                    "Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): no values found "
                                when (',' `notElem` head optionList) . failWithPhase Parsing $
                                    "Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): no comma found "
                                case liftA2 (,) noChangeMaybe changeMaybe of
                                    Nothing →
                                        failWithPhase Parsing $
                                            "Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): "
                                                <> head optionList
                                    Just val@(noChangeValue, changeValue) → do
                                        when (bc4 globalSettings /= val) . logWith LogInfo $
                                            "bit cost 4 state set to " <> show val
                                        pure (globalSettings{bc4 = val}, processedData, inSeedList)
                    "bc5" →
                        let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                            changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            changeMaybe = readMaybe changeString ∷ Maybe Double
                        in  do
                                when (length commandList /= length optionList) . failWithPhase Parsing $
                                    "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                when (null $ head optionList) . failWithPhase Parsing $
                                    "Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): no values found "
                                when (',' `notElem` head optionList) . failWithPhase Parsing $
                                    "Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): no comma found "
                                case liftA2 (,) noChangeMaybe changeMaybe of
                                    Nothing →
                                        failWithPhase Parsing $
                                            "Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): "
                                                <> head optionList
                                    Just val@(noChangeValue, changeValue) → do
                                        when (bc5 globalSettings /= val) . logWith LogInfo $
                                            "bit cost 5 state set to " <> show val
                                        pure (globalSettings{bc5 = val}, processedData, inSeedList)
                    "bc8" →
                        let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                            changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            changeMaybe = readMaybe changeString ∷ Maybe Double
                        in  do
                                when (length commandList /= length optionList) . failWithPhase Parsing $
                                    "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                when (null $ head optionList) . failWithPhase Parsing $
                                    "Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): no values found "
                                when (',' `notElem` head optionList) . failWithPhase Parsing $
                                    "Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): no comma found "
                                case liftA2 (,) noChangeMaybe changeMaybe of
                                    Nothing →
                                        failWithPhase Parsing $
                                            "Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): "
                                                <> head optionList
                                    Just val@(noChangeValue, changeValue) → do
                                        when (bc8 globalSettings /= val) . logWith LogInfo $
                                            "bit cost 8 state set to " <> show val
                                        pure (globalSettings{bc8 = val}, processedData, inSeedList)
                    "bc64" →
                        let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                            changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            changeMaybe = readMaybe changeString ∷ Maybe Double
                        in  do
                                when (length commandList /= length optionList) . failWithPhase Parsing $
                                    "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                when (null $ head optionList) . failWithPhase Parsing $
                                    "Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): no values found "
                                when (',' `notElem` head optionList) . failWithPhase Parsing $
                                    "Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): no comma found "
                                case liftA2 (,) noChangeMaybe changeMaybe of
                                    Nothing →
                                        failWithPhase Parsing $
                                            "Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): "
                                                <> head optionList
                                    Just val@(noChangeValue, changeValue) → do
                                        when (bc64 globalSettings /= val) . logWith LogInfo $
                                            "bit cost 64 state set to " <> show val
                                        pure (globalSettings{bc64 = val}, processedData, inSeedList)
                    "bcgt64" →
                        let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                            changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                            changeMaybe = readMaybe changeString ∷ Maybe Double
                        in  do
                                when (length commandList /= length optionList) . failWithPhase Parsing $
                                    "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                when (null $ head optionList) . failWithPhase Parsing $
                                    "Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): no values found "
                                when (',' `notElem` head optionList) . failWithPhase Parsing $
                                    "Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): no comma found "
                                case liftA2 (,) noChangeMaybe changeMaybe of
                                    Nothing →
                                        failWithPhase Parsing $
                                            "Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): "
                                                <> head optionList
                                    Just val@(noChangeValue, changeValue) → do
                                        when (bcgt64 globalSettings /= val) . logWith LogInfo $
                                            "bit cost > 64 state set to " <> show val
                                        pure (globalSettings{bcgt64 = val}, processedData, inSeedList)

                    -- processed above, but need here since put in different value
                    "criterion" → do
                        localCriterion ← case head optionList of
                            "parsimony" → pure Parsimony
                            "pmdl" → pure PMDL
                            "si" → pure SI
                            "mapa" → pure MAPA
                            "ncm" → pure NCM
                            val →
                                failWithPhase Parsing $
                                    "Error in 'set' command. Criterion '" <> val <> "' is not 'parsimony', 'ml', or 'pmdl'"

                        -- create lazy list of graph complexity indexed by number of network nodes--need leaf number for base tree complexity
                        (lGraphComplexityList, lRootComplexity) ← case localCriterion of
                            NCM
                                | origProcessedData /= emptyProcessedData →
                                    pure
                                        (IL.repeat (0.0, 0.0), U.calculateNCMRootCost origProcessedData)
                            NCM →
                                pure
                                    (IL.repeat (0.0, 0.0), U.calculateNCMRootCost processedData)
                            Parsimony → pure $ (IL.repeat (0.0, 0.0), 0.0)
                            val
                                | val `elem` [PMDL, SI, MAPA] →
                                    pure $
                                        (U.calculateGraphComplexity &&& U.calculateW15RootCost) processedData
                            val → failWithPhase Parsing $ "Optimality criterion not recognized: " <> show val

                        let lGraphFactor
                                | localCriterion `elem` [PMDL, SI, MAPA] = PMDLGraph
                                | otherwise = graphFactor globalSettings

                        logWith LogInfo $ case localCriterion of
                            NCM → unwords ["Optimality criterion set to", show NCM, "in -log (base 10) likelihood units"]
                            val
                                | val `elem` [PMDL, SI] →
                                    unwords ["Optimality criterion set to", show val, "Tree Complexity =", show . fst $ IL.head lGraphComplexityList, "bits"]
                            val → "Optimality criterion set to " <> show val

                        pure $
                            ( globalSettings
                                { optimalityCriterion = localCriterion
                                , graphComplexityList = lGraphComplexityList
                                , rootComplexity = lRootComplexity
                                , graphFactor = lGraphFactor
                                }
                            , processedData
                            , inSeedList
                            )

                    -- modify the behavior of resolutionCache softwired optimization
                    "compressresolutions" → do
                        localCriterion ← case toLower <$> head optionList of
                            "true" → pure True
                            "false" → pure False
                            val →
                                failWithPhase Parsing $
                                    "Error in 'set' command. CompressResolutions '" <> val <> "' is not 'true' or 'false'"
                        logWith LogInfo $ "CompressResolutions set to " <> show localCriterion
                        pure (globalSettings{compressResolutions = localCriterion}, processedData, inSeedList)

                    -- this not intended to be for users
                    "dynamicepsilon" → case readMaybe (head optionList) ∷ Maybe Double of
                        Nothing →
                            failWithPhase Parsing $
                                "Set option 'dynamicEpsilon' must be set to a double value >= 0.0 (e.g. dynamicepsilon:0.02): " <> head optionList
                        Just val
                            | val < 0.0 →
                                failWithPhase Parsing $
                                    "Set option 'dynamicEpsilon' must be set to a double value >= 0.0 (e.g. dynamicepsilon:0.02): " <> show val
                        Just localValue → do
                            logWith LogInfo $ "Dynamic Epsilon factor set to " <> head optionList
                            pure (globalSettings{dynamicEpsilon = 1.0 + (localValue * fractionDynamic globalSettings)}, processedData, inSeedList)
                    "finalassignment" → do
                        localMethod ← case head optionList of
                            "do" → pure DirectOptimization
                            "directoptimization" → pure DirectOptimization
                            "ia" → pure ImpliedAlignment
                            "impliedalignment" → pure ImpliedAlignment
                            _ →
                                failWithPhase Parsing $
                                    "Error in 'set' command. FinalAssignment  '"
                                        <> head optionList
                                        <> "' is not 'DirectOptimization (DO)' or 'ImpliedAlignment (IA)'"

                        when (graphType globalSettings == Tree) . logWith LogInfo $
                            "FinalAssignment set to " <> show localMethod

                        if (graphType globalSettings == Tree || localMethod == DirectOptimization)
                            then pure (globalSettings{finalAssignment = localMethod}, processedData, inSeedList)
                            else do
                                logWith LogInfo "FinalAssignment set to DO (ignoring IA option) for non-Tree graphs"
                                pure (globalSettings{finalAssignment = DirectOptimization}, processedData, inSeedList)
                    "graphfactor" → do
                        localMethod ← case toLower <$> head optionList of
                            "nopenalty" → pure NoNetworkPenalty
                            "w15" → pure Wheeler2015Network
                            "w23" → pure Wheeler2023Network
                            "pmdl" → pure PMDLGraph
                            val →
                                failWithPhase Parsing $
                                    "Error in 'set' command. GraphFactor  '" <> val <> "' is not 'NoPenalty', 'W15', 'W23', or 'PMDL'"
                        logWith LogInfo $ "GraphFactor set to " <> show localMethod
                        pure (globalSettings{graphFactor = localMethod}, processedData, inSeedList)
                    "graphssteepest" → case readMaybe (head optionList) ∷ Maybe Int of
                        Nothing →
                            failWithPhase Parsing $
                                "Set option 'graphsSteepest' must be set to an integer value (e.g. graphsSteepest:5): " <> head optionList
                        Just localValue → do
                            logWith LogInfo $ "GraphsStreepest set to " <> show localValue
                            pure (globalSettings{graphsSteepest = localValue}, processedData, inSeedList)
                    "graphtype" → do
                        localGraphType ← case head optionList of
                            "tree" → pure Tree
                            "softwired" → pure SoftWired
                            "hardwired" → pure HardWired
                            val →
                                failWithPhase Parsing $
                                    "Error in 'set' command. Graphtype '" <> val <> "' is not 'tree', 'hardwired', or 'softwired'"

                        let netPenalty = case localGraphType of
                                HardWired → NoNetworkPenalty
                                _ → graphFactor globalSettings

                        let settingResult = case localGraphType of
                                Tree → globalSettings{graphType = localGraphType}
                                _ →
                                    globalSettings
                                        { graphType = localGraphType
                                        , finalAssignment = DirectOptimization
                                        , graphFactor = netPenalty
                                        }
                        when (localGraphType /= Tree) $
                            logWith LogInfo $
                                unwords
                                    ["Graphtype set to", show localGraphType, "with graph factor NoPenalty and final assignment to DO"]

                        pure (settingResult, processedData, inSeedList)

                    -- In first to do stuff above also
                    "missingthreshold" → case readMaybe (head optionList) ∷ Maybe Int of
                        Nothing →
                            failWithPhase Parsing $
                                "Set option 'missingThreshold' must be set to an integer value (e.g. missingThreshold:50): " <> head optionList
                        Just localValue | localValue == missingThreshold globalSettings → pure (globalSettings, processedData, inSeedList)
                        Just localValue → do
                            logWith LogWarn $ "MissingThreshold set to " <> show localValue
                            pure (globalSettings{missingThreshold = localValue}, processedData, inSeedList)
                    "modelcomplexity" → case readMaybe (head optionList) ∷ Maybe Double of
                        Nothing →
                            failWithPhase Parsing $
                                "Set option 'modelComplexity' must be set to a double value (e.g. modelComplexity:123.456): " <> head optionList
                        Just localValue → do
                            logWith LogInfo $ "Model Complexity set to " <> head optionList
                            pure (globalSettings{modelComplexity = localValue}, processedData, inSeedList)

                    -- modify the behavior of rerooting character trees for all graph types
                    "multitraverse" → do
                        localCriterion ← case toLower <$> head optionList of
                            "true" → pure True
                            "false" → pure False
                            _ →
                                failWithPhase Parsing $
                                    "Error in 'set' command. MultiTraverse '" <> head optionList <> "' is not 'true' or 'false'"
                        logWith LogInfo $ "MultiTraverse set to " <> show localCriterion
                        pure (globalSettings{multiTraverseCharacters = localCriterion}, processedData, inSeedList)
                    "outgroup" →
                        let outTaxonName = T.pack $ filter (/= '"') $ head $ filter (/= "") $ fmap snd argList
                        in  case V.elemIndex outTaxonName leafNameVect of
                                Nothing →
                                    failWithPhase Parsing $
                                        unwords
                                            ["Error in 'set' command. Out-taxon", T.unpack outTaxonName, "not found in input leaf list", show $ T.unpack <$> leafNameVect]
                                Just outTaxonIndex → do
                                    logWith LogInfo $ "Outgroup set to " <> T.unpack outTaxonName
                                    pure (globalSettings{outgroupIndex = outTaxonIndex, outGroupName = outTaxonName}, processedData, inSeedList)
                    "partitioncharacter" → case head optionList of
                        localPartitionChar@[_] → do
                            when (localPartitionChar /= partitionCharacter globalSettings) . logWith LogInfo $
                                "PartitionCharacter set to '" <> head optionList <> "'"
                            pure (globalSettings{partitionCharacter = localPartitionChar}, processedData, inSeedList)
                        val →
                            failWithPhase Parsing $
                                "Error in 'set' command. Partitioncharacter '" <> val <> "' must be a single character"
                    "reportnaivedata" → do
                        localMethod ← case toLower <$> head optionList of
                            "true" → pure True
                            "false" → pure False
                            val → failWithPhase Parsing $ "Error in 'set' command. NeportNaive  '" <> val <> "' is not 'True' or 'False'"
                        logWith LogInfo $ "ReportNaiveData set to " <> show localMethod
                        pure (globalSettings{reportNaiveData = localMethod}, processedData, inSeedList)
                    "rootcost" → do
                        localMethod ← case toLower <$> head optionList of
                            "norootcost" → pure NoRootCost
                            "w15" → pure Wheeler2015Root
                            "pmdl" → pure PMDLRoot
                            "ml" → pure MLRoot
                            val → failWithPhase Parsing $ "Error in 'set' command. RootCost '" <> val <> "' is not 'NoRootCost', 'W15', or 'PMDL'"
                        lRootComplexity ← case localMethod of
                            NoRootCost → pure 0.0
                            val | val `elem` [Wheeler2015Root, PMDLRoot, MLRoot] → pure $ U.calculateW15RootCost processedData
                            val → failWithPhase Parsing $ "Error in 'set' command. No determined root complexity of '" <> show val <> "'"

                        logWith LogInfo $
                            unwords
                                ["RootCost set to", show localMethod, show lRootComplexity, "bits"]

                        pure (globalSettings{rootCost = localMethod, rootComplexity = lRootComplexity}, processedData, inSeedList)
                    "seed" → case readMaybe (head optionList) ∷ Maybe Int of
                        Nothing →
                            failWithPhase Parsing $
                                "Set option 'seed' must be set to an integer value (e.g. seed:123): " <> head optionList
                        Just localValue → do
                            logWith LogInfo $ "Random Seed set to " <> head optionList
                            pure (globalSettings{seed = localValue}, processedData, randomIntList localValue)
                    "softwiredmethod" → do
                        localMethod ← case toLower <$> head optionList of
                            "naive" → pure Naive
                            "exhaustive" → pure Naive
                            "resolutioncache" → pure ResolutionCache
                            val →
                                failWithPhase Parsing $
                                    "Error in 'set' command. SoftwiredMethod  '" <> val <> "' is not 'Exhaustive' or 'ResolutionCache'"
                        logWith LogInfo $ "SoftwiredMethod " <> show localMethod
                        pure (globalSettings{softWiredMethod = localMethod}, processedData, inSeedList)

                    -- modify the use of Network Add heurisitcs in network optimization
                    "usenetaddheuristic" → do
                        localCriterion ← case toLower <$> head optionList of
                            "true" → pure True
                            "false" → pure False
                            val →
                                failWithPhase Parsing $
                                    "Error in 'set' command. UseNetAddHeuristic '" <> val <> "' is not 'true' or 'false'"
                        logWith LogInfo $ "UseNetAddHeuristic set to " <> show localCriterion
                        pure (globalSettings{useNetAddHeuristic = localCriterion}, processedData, inSeedList)

                    -- these not intended for users
                    "jointhreshold" → case readMaybe (head optionList) ∷ Maybe Double of
                        Nothing →
                            failWithPhase Parsing $
                                "Set option 'joinThreshold' must be set to an double value >= 1.0 (e.g. joinThreshold:1.17): " <> head optionList
                        Just localValue
                            | localValue < 1.0 →
                                failWithPhase Parsing $
                                    "Set option 'joinThreshold' must be set to a double value >= 1.0 (e.g. joinThreshold:1.17): " <> show localValue
                        Just localValue → do
                            logWith LogInfo $ "JoinThreshold set to " <> show localValue
                            pure (globalSettings{unionThreshold = localValue}, processedData, inSeedList)

                    -- parallel strategy settings options
                    "defparstrat" → do
                        localMethod ← case head optionList of
                            "r0" → pure R0
                            "rpar" → pure RPar
                            "rseq" → pure RSeq
                            "rdeepseq" → pure RDeepSeq
                            val →
                                failWithPhase Parsing $
                                    "Error in 'set' command. DefParStrat  '" <> val <> "' is not 'r0', 'WrPar', 'rSeq', or 'rDeepSeq'"
                        logWith LogInfo $ "DefParStrat set to " <> show localMethod
                        pure (globalSettings{defaultParStrat = localMethod}, processedData, inSeedList)
                    "lazyparstrat" → do
                        localMethod ← case head optionList of
                            "r0" → pure R0
                            "rpar" → pure RPar
                            "rseq" → pure RSeq
                            "rdeepseq" → pure RDeepSeq
                            val →
                                failWithPhase Parsing $
                                    "Error in 'set' command. DefParStrat  '" <> val <> "' is not 'r0', 'WrPar', 'rSeq', or 'rDeepSeq'"
                        logWith LogInfo $ "LazyParStrat set to " <> show localMethod
                        pure (globalSettings{lazyParStrat = localMethod}, processedData, inSeedList)
                    "strictparstrat" → do
                        localMethod ← case head optionList of
                            "r0" → pure R0
                            "rpar" → pure RPar
                            "rseq" → pure RSeq
                            "rdeepseq" → pure RDeepSeq
                            val →
                                failWithPhase Parsing $
                                    "Error in 'set' command. DefParStrat  '" <> val <> "' is not 'r0', 'WrPar', 'rSeq', or 'rDeepSeq'"
                        logWith LogInfo $ "StrictParStrat set to " <> show localMethod
                        pure (globalSettings{strictParStrat = localMethod}, processedData, inSeedList)

                    -- modify the use of implied alkignemnt in heuristics
                    "useia" → do
                        localCriterion ← case toLower <$> head optionList of
                            "true" → pure True
                            "false" → pure False
                            val →
                                failWithPhase Parsing $
                                    "Error in 'set' command. UseIA '" <> head optionList <> "' is not 'true' or 'false'"
                        logWith LogInfo $ "UseIA set to " <> show localCriterion
                        pure (globalSettings{useIA = localCriterion}, processedData, inSeedList)
                    val → do
                        logWith LogWarn $ "Warning: Unrecognized/missing 'set' option in " <> show argList
                        pure (globalSettings, processedData, inSeedList)


{- |
'reportCommand' takes report options, current data and graphs and returns a
(potentially large) String to print and the channel to print it to and write mode
overwrite/append if global settings reportNaiveData is True then need to rediagnose
graph with processed data since naiveData was sent to command and will not match
what is in the optimized graphs.
-}
reportCommand
    ∷ GlobalSettings
    → [Argument]
    → ([NameText], [(NameText, NameText)])
    → Int
    → String
    → ProcessedData
    → [ReducedPhylogeneticGraph]
    → [ReducedPhylogeneticGraph]
    → [[VertexCost]]
    → PhyG (String, String, String)
reportCommand globalSettings argList excludeRename numInputFiles crossReferenceString processedData curGraphs supportGraphs pairwiseDistanceMatrix =
    let argListWithoutReconcileCommands = filter ((`notElem` VER.reconcileArgList) . fst) argList
        -- check for balances double quotes and only one pair
        outFileNameList = filter (/= "") $ fmap snd argListWithoutReconcileCommands -- argList
        commandList = fmap (fmap C.toLower) $ filter (/= "") $ fmap fst argListWithoutReconcileCommands
    in  -- reconcileList = filter (/= "") $ fmap fst argList

        if length outFileNameList > 1
            then do failWithPhase Outputting ("Report can only have one file name: " <> show outFileNameList <> " " <> show argList)
            else
                let checkCommandList = checkCommandArgs "report" commandList VER.reportArgList
                    outfileName =
                        if null outFileNameList
                            then "stderr"
                            else tail $ L.init $ head outFileNameList
                    writeMode =
                        if "overwrite" `elem` commandList
                            then "overwrite"
                            else "append"
                in  -- error too harsh, lose everything else
                    -- if (null $ filter (/= "overwrite") $ filter (/= "append") commandList) then errorWithoutStackTrace ("Error: Missing 'report' option in " <> show commandList)
                    -- else
                    if not checkCommandList
                        then errorWithoutStackTrace ("Unrecognized command in report: " <> show argList)
                        else -- This for reconciled data

                            if "crossrefs" `elem` commandList
                                then
                                    let dataString = crossReferenceString
                                    in  pure (dataString, outfileName, writeMode)
                                else
                                    if "data" `elem` commandList
                                        then
                                            let blocks = thd3 processedData
                                                numChars = V.sum $ fmap (V.length . thd3) blocks
                                                dataString = phyloDataToString 0 $ thd3 processedData
                                                baseData =
                                                    [ ["Input data contained:"]
                                                    , ["", show (length $ fst3 processedData) <> " terminal taxa"]
                                                    , ["", show numInputFiles <> " input data files"]
                                                    , ["", show (length blocks) <> " character blocks"]
                                                    , ["", show numChars <> " total characters"]
                                                    ]
                                                leafNames = V.toList (T.unpack <$> fst3 processedData)
                                                leafField = ["Terminal taxa:"] : fmap (" " :) (SL.chunksOf 10 leafNames)
                                                excludedTaxa =
                                                    if (not . null . fst) excludeRename
                                                        then T.unpack <$> fst excludeRename
                                                        else ["None"]
                                                excludedField = ["Excluded taxa:"] : fmap (" " :) (SL.chunksOf 10 excludedTaxa)
                                                renameFirstList = fmap (((: []) . T.unpack) . fst) (snd excludeRename)
                                                renameSecondList = fmap (((: []) . T.unpack) . snd) (snd excludeRename)
                                                renamePairList =
                                                    if (not . null . snd) excludeRename
                                                        then (" " :) <$> zipWith (<>) renameFirstList renameSecondList
                                                        else [[" ", "None", "None"]]
                                                renameField = ["Renamed taxa:", "New Name", "Original Name"] : renamePairList
                                                charInfoFields = ["Index", "Block", "Name", "Type", "Activity", "Weight", "Prealigned", "Alphabet", "TCM"]
                                            in  pure
                                                    ( CSV.genCsvFile
                                                        (baseData <> [[""]] <> leafField <> [[""]] <> excludedField <> [[""]] <> renameField <> [[""]] <> (charInfoFields : dataString))
                                                    , outfileName
                                                    , writeMode
                                                    )
                                        else
                                            if "diagnosis" `elem` commandList
                                                then do
                                                    -- action :: SimpleGraph -> ReducedPhylogeneticGraph
                                                    let action = TRAV.multiTraverseFullyLabelGraphReduced globalSettings processedData False False Nothing
                                                    pTraverse ← getParallelChunkTraverse
                                                    rediagnosedGraph ← pTraverse action (fmap fst5 curGraphs)

                                                    let curGraphs' =
                                                            if not (reportNaiveData globalSettings)
                                                                then curGraphs
                                                                else rediagnosedGraph
                                                    -- else PU.seqParMap (parStrategy $ strictParStrat globalSettings) (TRAV.multiTraverseFullyLabelGraphReduced globalSettings processedData False False Nothing) (fmap fst5 curGraphs)

                                                    let dataString = CSV.genCsvFile $ concatMap (getGraphDiagnosis globalSettings processedData) (zip curGraphs' [0 .. (length curGraphs' - 1)])
                                                    if null curGraphs
                                                        then do
                                                            logWith LogInfo "No graphs to diagnose\n"
                                                            pure ("No graphs to diagnose", outfileName, writeMode)
                                                        else do
                                                            logWith
                                                                LogInfo
                                                                ("Diagnosing " <> show (length curGraphs) <> " graphs at minimum cost " <> show (minimum $ fmap snd5 curGraphs) <> "\n")
                                                            pure (dataString, outfileName, writeMode)
                                                else
                                                    if "displaytrees" `elem` commandList
                                                        then do
                                                            -- need to specify -O option for multiple graphs
                                                            -- TODO parallelize
                                                            rediagnodesGraphs ←
                                                                mapM (TRAV.multiTraverseFullyLabelGraph globalSettings processedData False False Nothing) (fmap fst5 curGraphs)
                                                            let inputDisplayVVList = fmap fth6 rediagnodesGraphs
                                                            let costList = fmap snd5 curGraphs
                                                            let displayCostListList = fmap GO.getDisplayTreeCostList rediagnodesGraphs
                                                            let displayInfoString =
                                                                    if ("dot" `elem` commandList) || ("dotpdf" `elem` commandList)
                                                                        then ("//DisplayTree costs : " <> show (fmap (sum . fst) displayCostListList, displayCostListList))
                                                                        else -- newick

                                                                            let middle = fmap bracketToCurly $ show (fmap (sum . fst) displayCostListList, displayCostListList)
                                                                            in  ("[DisplayTree costs : " <> middle <> "]")

                                                            let treeIndexStringList =
                                                                    if ("dot" `elem` commandList) || ("dotpdf" `elem` commandList)
                                                                        then fmap (((<> "\n") . ("//Canonical Tree " <>)) . show) [0 .. (length inputDisplayVVList - 1)]
                                                                        else -- newick
                                                                            fmap (((<> "]\n") . ("[Canonical Tree " <>)) . show) [0 .. (length inputDisplayVVList - 1)]
                                                            let canonicalGraphPairList = zip treeIndexStringList inputDisplayVVList
                                                            let blockStringList = unlines (fmap (outputBlockTrees commandList costList (outgroupIndex globalSettings)) canonicalGraphPairList)
                                                            -- graphString = outputGraphString commandList (outgroupIndex globalSettings) (fmap thd6 curGraphs) (fmap snd6 curGraphs)

                                                            if null curGraphs || graphType globalSettings /= SoftWired
                                                                then do
                                                                    logWith LogInfo "No soft-wired graphs to report display trees\n"
                                                                    pure ("No soft-wired graphs to report display trees", outfileName, writeMode)
                                                                else pure (displayInfoString <> "\n" <> blockStringList, outfileName, writeMode)
                                                        else
                                                            if ("graphs" `elem` commandList) && ("reconcile" `notElem` commandList)
                                                                then -- else if (not .null) (L.intersect ["graphs", "newick", "dot", "dotpdf"] commandList) then

                                                                    let graphString = outputGraphString commandList (outgroupIndex globalSettings) (fmap thd5 curGraphs) (fmap snd5 curGraphs)
                                                                    in  if null curGraphs
                                                                            then do
                                                                                logWith LogInfo "No graphs to report\n"
                                                                                pure ("No graphs to report", outfileName, writeMode)
                                                                            else do
                                                                                logWith
                                                                                    LogInfo
                                                                                    ("Reporting " <> show (length curGraphs) <> " graph(s) at minimum cost " <> show (minimum $ fmap snd5 curGraphs) <> "\n")
                                                                                pure (graphString, outfileName, writeMode)
                                                                else
                                                                    if "ia" `elem` commandList || "impliedalignment" `elem` commandList
                                                                        then
                                                                            if null curGraphs
                                                                                then do
                                                                                    logWith LogInfo "No graphs to create implied alignments\n"
                                                                                    pure ("No impliedAlgnments to report", outfileName, writeMode)
                                                                                else
                                                                                    let includeMissing = elem "includemissing" commandList
                                                                                        concatSeqs = elem "concatenate" commandList
                                                                                    in  do
                                                                                            iaContentList ←
                                                                                                mapM
                                                                                                    (getImpliedAlignmentString globalSettings (includeMissing || concatSeqs) concatSeqs processedData)
                                                                                                    (zip curGraphs [0 .. (length curGraphs - 1)])
                                                                                            logWith
                                                                                                LogInfo
                                                                                                "\tWarning: Prealigned sequence data with non-additive type costs (all change values equal) have been recoded to non-additive characters and will not appear in implied alignment output.\n"
                                                                                            pure (concat iaContentList, outfileName, writeMode)
                                                                        else
                                                                            if "pairdist" `elem` commandList
                                                                                then
                                                                                    let nameData = L.intercalate "," (V.toList (T.unpack <$> fst3 processedData)) <> "\n"
                                                                                    in  do
                                                                                            pairwiseDistanceMatrix' ← D.getPairwiseDistances processedData
                                                                                            let dataString = CSV.genCsvFile $ fmap (fmap show) pairwiseDistanceMatrix'
                                                                                            pure (nameData <> dataString, outfileName, writeMode)
                                                                                else
                                                                                    if "reconcile" `elem` commandList
                                                                                        then do
                                                                                            recResult ← R.makeReconcileGraph VER.reconcileArgList argList (fmap fst5 curGraphs)
                                                                                            -- let (reconcileString, ) = recResult
                                                                                            let (_, reconcileGraph) = recResult
                                                                                            let reconcileString = outputGraphString commandList (outgroupIndex globalSettings) [GO.convertSimpleToDecoratedGraph reconcileGraph] [0]
                                                                                            if null curGraphs
                                                                                                then do
                                                                                                    logWith LogInfo "No graphs to reconcile\n"
                                                                                                    pure ([], outfileName, writeMode)
                                                                                                else pure (reconcileString, outfileName, writeMode)
                                                                                        else
                                                                                            if "search" `elem` commandList
                                                                                                then
                                                                                                    let dataString' = fmap showSearchFields $ reverse $ searchData globalSettings
                                                                                                        -- reformat the "search" command fields a bit
                                                                                                        dataString = processSearchFields dataString'
                                                                                                        sysInfoData =
                                                                                                            "System Info, OS: "
                                                                                                                <> SI.os
                                                                                                                <> ", Chip Arch: "
                                                                                                                <> SI.arch
                                                                                                                <> ", Compiler: "
                                                                                                                <> SI.compilerName
                                                                                                                <> " "
                                                                                                                <> DV.showVersion SI.compilerVersion
                                                                                                                <> ", Compile Date: "
                                                                                                                <> (__DATE__ <> " " <> __TIME__)
                                                                                                        cpuInfoString =
                                                                                                            if SI.os /= "linux"
                                                                                                                then "CPU Info, No /proc/cpuinfo on darwin"
                                                                                                                else
                                                                                                                    let cpuInfoM = SIOU.unsafePerformIO SC.tryGetCPUs
                                                                                                                    in  if isNothing cpuInfoM
                                                                                                                            then "CPU Info, Couldn't parse CPU Info"
                                                                                                                            else
                                                                                                                                "CPU Info, Physical Processors: "
                                                                                                                                    <> show (SC.physicalProcessors (fromJust cpuInfoM))
                                                                                                                                    <> ", Physical Cores: "
                                                                                                                                    <> show (SC.physicalCores (fromJust cpuInfoM))
                                                                                                                                    <> ", Logical Cores: "
                                                                                                                                    <> show (SC.logicalCores (fromJust cpuInfoM))
                                                                                                        baseData = sysInfoData <> "\n" <> cpuInfoString <> "\nSearchData\nRandom seed, " <> show (seed globalSettings) <> "\n"
                                                                                                        charInfoFields =
                                                                                                            [ "Command"
                                                                                                            , "Arguments"
                                                                                                            , "Min cost in"
                                                                                                            , "Max cost in"
                                                                                                            , "Num graphs in"
                                                                                                            , "Min cost out"
                                                                                                            , "Max cost out"
                                                                                                            , "Num graphs out"
                                                                                                            , "CPU time (secs)"
                                                                                                            , "Comment"
                                                                                                            ]
                                                                                                    in  pure (baseData <> CSV.genCsvFile (charInfoFields : dataString), outfileName, writeMode)
                                                                                                else
                                                                                                    if "support" `elem` commandList
                                                                                                        then
                                                                                                            let graphString = outputGraphStringSimple commandList (outgroupIndex globalSettings) (fmap fst5 supportGraphs) (fmap snd5 supportGraphs)
                                                                                                            in  -- trace ("Rep Sup: " <> (LG.prettify $ fst5 $ head supportGraphs)) (
                                                                                                                if null supportGraphs
                                                                                                                    then do
                                                                                                                        logWith LogInfo "\tNo support graphs to report\n"
                                                                                                                        pure ([], outfileName, writeMode)
                                                                                                                    else do
                                                                                                                        logWith LogInfo ("Reporting " <> show (length curGraphs) <> " support graph(s)" <> "\n")
                                                                                                                        pure (graphString, outfileName, writeMode)
                                                                                                        else -- )

                                                                                                            if "tnt" `elem` commandList
                                                                                                                then
                                                                                                                    if null curGraphs
                                                                                                                        then do
                                                                                                                            logWith LogInfo "No graphs to create implied alignments for TNT output\n"
                                                                                                                            pure ("No impliedAlgnments for TNT to report", outfileName, writeMode)
                                                                                                                        else do
                                                                                                                            -- action :: SinmpleGraph -> ReducedPhylogeneticGraph
                                                                                                                            let action = TRAV.multiTraverseFullyLabelGraph globalSettings processedData False False Nothing
                                                                                                                            pTraverse ← getParallelChunkTraverse
                                                                                                                            reoptimizedGraphs ← pTraverse action (fmap fst5 curGraphs)
                                                                                                                            let curGraphs' =
                                                                                                                                    if not (reportNaiveData globalSettings)
                                                                                                                                        then (fmap GO.convertReduced2PhylogeneticGraph curGraphs)
                                                                                                                                        else reoptimizedGraphs
                                                                                                                            -- PU.seqParMap
                                                                                                                            --     (parStrategy $ strictParStrat globalSettings)
                                                                                                                            --     (TRAV.multiTraverseFullyLabelGraph globalSettings processedData False False Nothing)
                                                                                                                            --     (fmap fst5 curGraphs)

                                                                                                                            tntContentList' ← mapM (getTNTString globalSettings processedData) (zip curGraphs' [0 .. (length curGraphs' - 1)])
                                                                                                                            let tntContentList = concat tntContentList'
                                                                                                                            pure (tntContentList, outfileName, writeMode)
                                                                                                                else do
                                                                                                                    logWith LogWarn ("\nUnrecognized/missing report option in " <> show commandList <> " defaulting to 'graphs'" <> "\n")
                                                                                                                    let graphString = outputGraphString commandList (outgroupIndex globalSettings) (fmap thd5 curGraphs) (fmap snd5 curGraphs)
                                                                                                                    if null curGraphs
                                                                                                                        then do
                                                                                                                            logWith LogInfo "No graphs to report\n"
                                                                                                                            pure ("No graphs to report", outfileName, writeMode)
                                                                                                                        else do
                                                                                                                            logWith
                                                                                                                                LogInfo
                                                                                                                                ("Reporting " <> show (length curGraphs) <> " graph(s) at minimum cost " <> show (minimum $ fmap snd5 curGraphs) <> "\n")
                                                                                                                            pure (graphString, outfileName, writeMode)
    where
        bracketToCurly a =
            if a == '('
                then '{'
                else
                    if a == ')'
                        then '}'
                        else a
