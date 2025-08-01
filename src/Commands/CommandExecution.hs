{-# LANGUAGE CPP #-}

{- |
Module to coordinate command execution.
-}
module Commands.CommandExecution (
    executeCommands,
    executeRenameReblockCommands,
    getDataListList,
    getDataListList',
) where

import Commands.CommandUtilities
import Commands.Transform qualified as TRANS
import Commands.Verify qualified as VER
import Control.Arrow ((&&&))
import Control.Monad (unless, when)
import Control.Monad.IO.Class (MonadIO (..))
import Data.Bifunctor (bimap)
import Data.CSV qualified as CSV
import Data.Char
import Data.Char qualified as C
import Data.Foldable (fold)
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
    → [ReducedPhylogeneticGraph]
    → [Command]
    → Bool
    → PhyG ([ReducedPhylogeneticGraph], GlobalSettings, [ReducedPhylogeneticGraph])
executeCommands globalSettings excludeRename numInputFiles crossReferenceString origProcessedData processedData reportingData curGraphs supportGraphList commandList isFirst = case commandList of
    [] → pure (curGraphs, globalSettings, supportGraphList)
    (firstOption, firstArgs) : otherCommands → case firstOption of
        -- skip "Read" and "Rename "commands already processed
        Read → error ("Read command should already have been processed: " <> show (firstOption, firstArgs))
        Rename → error ("Rename command should already have been processed: " <> show (firstOption, firstArgs))
        Reblock → error ("Reblock command should already have been processed: " <> show (firstOption, firstArgs))
        Run → error ("Run command should already have been processed: " <> show (firstOption, firstArgs))
        -- other commands
        Build → do
            (elapsedSeconds, newGraphList') ←
                timeOp . pure $
                    B.buildGraph firstArgs globalSettings processedData
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
                supportGraphList
                otherCommands
                isFirst
        Refine → do
            (elapsedSeconds, newGraphList') ←
                timeOp . pure $
                    REF.refineGraph firstArgs globalSettings processedData curGraphs

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
                supportGraphList
                otherCommands
                isFirst
        Fuse → do
            (elapsedSeconds, newGraphList) ←
                timeOp $
                    REF.fuseGraphs firstArgs globalSettings processedData curGraphs

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
                supportGraphList
                otherCommands
                isFirst
        Report → do
            let doDotPDF = elem "dotpdf" $ fmap (fmap toLower . fst) firstArgs
            let collapse' = elem "collapse" $ fmap (fmap toLower . fst) firstArgs
            let noCollapse' = elem "nocollapse" $ fmap (fmap toLower . fst) firstArgs
            let reconcile = any ((== "reconcile") . fst) firstArgs

            -- set default collapse for dotPDF to True, False otherwise
            let collapse -- this will casue problems with reconcile
                    | reconcile = False
                    | collapse' = True
                    | noCollapse' = False
                    --  | doDotPDF = True
                    | otherwise = False

            let curGraphs'
                    | not collapse = curGraphs
                    | otherwise = U.collapseReducedGraph <$> curGraphs

            -- use 'temp' updated graphs s don't repeatedly add model and root complexities
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

            unless (null reportString) . logWith LogInfo $ "Report writing to \"" <> outFile <> "\"\n"

            case reportString of
                "" →
                    executeCommands
                        globalSettings
                        excludeRename
                        numInputFiles
                        crossReferenceString
                        origProcessedData
                        processedData
                        reportingData
                        curGraphs
                        supportGraphList
                        otherCommands
                        isFirst
                _ | doDotPDF → do
                    let reportString' = changeDotPreamble "digraph {" "digraph G {\n\trankdir = LR;\tedge [colorscheme=spectral11];\tnode [shape = none];\n" reportString
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
                        supportGraphList
                        otherCommands
                        isFirst
                _ → do
                    case outFile of
                        "stderr" → liftIO $ hPutStr stderr reportString
                        "stdout" → liftIO $ putStr reportString
                        _ → case writeMode of
                            "overwrite" → liftIO $ writeFile outFile reportString
                            "append" → liftIO $ appendFile outFile reportString
                            _ →
                                failWithPhase Parsing $
                                    "Error 'report' command not properly formatted" <> show reportStuff

                    executeCommands
                        globalSettings
                        excludeRename
                        numInputFiles
                        crossReferenceString
                        origProcessedData
                        processedData
                        reportingData
                        curGraphs
                        supportGraphList
                        otherCommands
                        isFirst
        Search → do
            (elapsedSeconds, output) ←
                timeOp $
                    S.search firstArgs globalSettings processedData curGraphs
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
                supportGraphList
                otherCommands
                isFirst
        Select → do
            (elapsedSeconds, newGraphList) ←
                timeOp $
                    GO.selectPhylogeneticGraphReduced firstArgs (outgroupIndex globalSettings) curGraphs
            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
            let newSearchData = searchInfo : searchData globalSettings
            let typeSelected = case firstArgs of
                    [] → "best"
                    (a, _) : _ → C.toLower <$> a
            logWith LogInfo $ unwords ["Selecting", typeSelected, "graphs", "\n"]
            executeCommands
                (globalSettings{searchData = newSearchData})
                excludeRename
                numInputFiles
                crossReferenceString
                origProcessedData
                processedData
                reportingData
                newGraphList
                supportGraphList
                otherCommands
                isFirst
        Set → do
            -- if set changes graph aspects--may need to reoptimize
            (newGlobalSettings, newProcessedData) ← setCommand firstArgs globalSettings reportingData processedData isFirst
            let needReoptimize = requireReoptimization globalSettings newGlobalSettings
            newGraphList ←
                if not needReoptimize
                    then pure curGraphs
                    else -- then logWith LogInfo "No need to reoptimize graphs\n" $> curGraphs
                    do
                        logWith LogInfo "Reoptimizing graphs\n"
                        -- TODO should be parallel
                        mapM (TRAV.multiTraverseFullyLabelGraphReduced newGlobalSettings newProcessedData True True Nothing) $ fst5 <$> curGraphs

            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList 0 "No Comment"
            let newSearchData = searchInfo : searchData newGlobalSettings

            executeCommands
                (newGlobalSettings{searchData = newSearchData})
                excludeRename
                numInputFiles
                crossReferenceString
                origProcessedData
                processedData
                reportingData
                newGraphList
                supportGraphList
                otherCommands
                False
        Swap → do
            (elapsedSeconds, newGraphList) ←
                timeOp $
                    REF.swapMaster firstArgs globalSettings processedData curGraphs
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
                supportGraphList
                otherCommands
                isFirst
        Support → do
            (elapsedSeconds, newSupportGraphList') ←
                timeOp . pure $
                    SUP.supportGraph firstArgs globalSettings processedData curGraphs

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
                (supportGraphList <> newSupportGraphList)
                otherCommands
                isFirst
        Transform → do
            (elapsedSeconds, (newGS, newOrigData, newProcessedData, newGraphs)) ←
                timeOp $
                    TRANS.transform firstArgs globalSettings origProcessedData processedData curGraphs

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
                supportGraphList
                otherCommands
                isFirst
        val → error $ "Command " <> show val <> " not recognized/implemented"


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
if seed list is empty [] then processes first set--confusing--should be refactored
-}
setCommand ∷ [Argument] → GlobalSettings → ProcessedData → ProcessedData → Bool → PhyG (GlobalSettings, ProcessedData)
setCommand argList globalSettings origProcessedData processedData isFirst =
    let material ∷ Argument → Bool
        material (k, v) = not $ null k || null v
        normalize = fmap C.toLower
        niceArgs = bimap normalize normalize <$> filter material argList
        (commandList, optionList) = unzip niceArgs
        checkCommandList = checkCommandArgs "set" commandList VER.setArgList
        leafNameVect = fst3 processedData
    in  case niceArgs of
            [] → if (not $ null (fst $ head argList)) then failWithPhase Computing $ "Error processing 'set' command, requires second argument: " <> (fst $ head argList)
                 else failWithPhase Computing $ "Attempting to process a non-existant 'set' command " <> (fst $ head argList)
            (firstCommand, firstOption) : _otherArgs → do
                when (not checkCommandList) . failWithPhase Parsing $
                    "Unrecognized command in 'set': " <> show argList

                when (length commandList > 1 || length optionList > 1) . failWithPhase Parsing $
                    "Set option error: can only have one set argument for each command: " <> show (commandList, optionList)

                case isFirst of
                    True → case firstCommand of
                        "partitioncharacter" → case firstOption of
                            localPartitionChar@[_] →
                                do
                                    logWith LogInfo $ "PartitionCharacter set to '" <> localPartitionChar <> "'\n"
                                    let x = globalSettings{partitionCharacter = localPartitionChar}
                                    pure (x, processedData)
                            val →
                                failWithPhase Parsing $
                                    "Error in 'set' command. Partitioncharacter '" <> val <> "' must be a single character\n"
                        "missingthreshold" → case readMaybe (firstOption) ∷ Maybe Int of
                            Nothing →
                                failWithPhase Parsing $
                                    "Set option 'missingThreshold' must be set to an integer value (e.g. missingThreshold:50): " <> firstOption
                            Just val
                                | val < 0 || 100 < val →
                                    failWithPhase Parsing $
                                        "Set option 'missingThreshold' must be set to an integer value between 0 and 100: " <> firstOption
                            Just val → do
                                logWith LogInfo $ "MissingThreshold set to " <> firstOption <> "\n"
                                pure (globalSettings{missingThreshold = val}, processedData)

                        -- sets root cost as well-- need in both places--one to process data and one to
                        -- keep in current global
                        -- MUST be set aheard of data packing so correct--otherwise alhobaet etc modified
                        "criterion" → do
                            localCriterion ← case firstOption of
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
                                -- NCM | origProcessedData /= emptyProcessedData →
                                --    pure (IL.repeat (0.0, 0.0), U.calculateNCMRootCost origProcessedData)

                                NCM → pure (IL.repeat (0.0, 0.0), U.calculateNCMRootCost origProcessedData)
                                Parsimony → pure $ (IL.repeat (0.0, 0.0), 0.0)
                                MAPA → pure $ (IL.repeat (0.0, 0.0), U.calculateMAPARootCost origProcessedData)
                                val
                                    | val `elem` [PMDL, SI] →
                                        pure $ (U.calculateGraphComplexity &&& (U.calculatePMDLRootCost True)) origProcessedData
                                val → failWithPhase Parsing $ "Optimality criterion not recognized: " <> show val

                            let lGraphFactor
                                    | localCriterion `elem` [PMDL, SI, MAPA] = PMDLGraph
                                    | otherwise = graphFactor globalSettings

                            logWith LogInfo $ case localCriterion of
                                NCM → unwords ["Optimality criterion set to", show NCM, "in -log (base 10) likelihood units\n"]
                                val
                                    | val `elem` [PMDL, SI] →
                                        unwords ["Optimality criterion set to", show val, "Tree Complexity =", show . fst $ IL.head lGraphComplexityList, "bits\n"]
                                val → "Optimality criterion set to " <> show val <> "\n"

                            pure $
                                ( globalSettings
                                    { optimalityCriterion = localCriterion
                                    , graphComplexityList = lGraphComplexityList
                                    , rootComplexity = lRootComplexity
                                    , graphFactor = lGraphFactor
                                    }
                                , processedData
                                )
                        "bc2" →
                            let (noChangeString, changeString) = changingStrings firstOption
                                noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                                changeMaybe = readMaybe changeString ∷ Maybe Double
                            in  do
                                    when (length commandList /= length optionList) . failWithPhase Parsing $
                                        "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                    when (null $ firstOption) . failWithPhase Parsing $
                                        "Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): no values found "
                                    when (',' `notElem` firstOption) . failWithPhase Parsing $
                                        "Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): no comma found "
                                    case liftA2 (,) noChangeMaybe changeMaybe of
                                        Nothing →
                                            failWithPhase Parsing $
                                                "Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): "
                                                    <> firstOption
                                        Just val@(_, _) → do
                                            when (bc2 globalSettings /= val) . logWith LogInfo $
                                                "bit cost 2 state set to " <> show val <> "\n"
                                            pure (globalSettings{bc2 = val}, processedData)
                        "bc4" →
                            let (noChangeString, changeString) = changingStrings firstOption
                                noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                                changeMaybe = readMaybe changeString ∷ Maybe Double
                            in  do
                                    when (length commandList /= length optionList) . failWithPhase Parsing $
                                        "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                    when (null $ firstOption) . failWithPhase Parsing $
                                        "Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): no values found "
                                    when (',' `notElem` firstOption) . failWithPhase Parsing $
                                        "Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): no comma found "
                                    case liftA2 (,) noChangeMaybe changeMaybe of
                                        Nothing →
                                            failWithPhase Parsing $
                                                "Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): "
                                                    <> firstOption
                                        Just val@(_, _) → do
                                            when (bc4 globalSettings /= val) . logWith LogInfo $
                                                "bit cost 4 state set to " <> show val <> "\n"
                                            pure (globalSettings{bc4 = val}, processedData)
                        "bc5" →
                            let (noChangeString, changeString) = changingStrings firstOption
                                noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                                changeMaybe = readMaybe changeString ∷ Maybe Double
                            in  do
                                    when (length commandList /= length optionList) . failWithPhase Parsing $
                                        "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                    when (null $ firstOption) . failWithPhase Parsing $
                                        "Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): no values found "
                                    when (',' `notElem` firstOption) . failWithPhase Parsing $
                                        "Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): no comma found "
                                    case liftA2 (,) noChangeMaybe changeMaybe of
                                        Nothing →
                                            failWithPhase Parsing $
                                                "Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): "
                                                    <> firstOption
                                        Just val@(_, _) → do
                                            when (bc5 globalSettings /= val) . logWith LogInfo $
                                                "bit cost 5 state set to " <> show val <> "\n"
                                            pure (globalSettings{bc5 = val}, processedData)
                        "bc8" →
                            let (noChangeString, changeString) = changingStrings firstOption
                                noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                                changeMaybe = readMaybe changeString ∷ Maybe Double
                            in  do
                                    when (length commandList /= length optionList) . failWithPhase Parsing $
                                        "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                    when (null $ firstOption) . failWithPhase Parsing $
                                        "Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): no values found "
                                    when (',' `notElem` firstOption) . failWithPhase Parsing $
                                        "Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): no comma found "
                                    case liftA2 (,) noChangeMaybe changeMaybe of
                                        Nothing →
                                            failWithPhase Parsing $
                                                "Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): "
                                                    <> firstOption
                                        Just val@(_, _) → do
                                            when (bc8 globalSettings /= val) . logWith LogInfo $
                                                "bit cost 8 state set to " <> show val <> "\n"
                                            pure (globalSettings{bc8 = val}, processedData)
                        "bc64" →
                            let (noChangeString, changeString) = changingStrings firstOption
                                noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                                changeMaybe = readMaybe changeString ∷ Maybe Double
                            in  do
                                    when (length commandList /= length optionList) . failWithPhase Parsing $
                                        "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                    when (null $ firstOption) . failWithPhase Parsing $
                                        "Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): no values found "
                                    when (',' `notElem` firstOption) . failWithPhase Parsing $
                                        "Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): no comma found "
                                    case liftA2 (,) noChangeMaybe changeMaybe of
                                        Nothing →
                                            failWithPhase Parsing $
                                                "Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): "
                                                    <> firstOption
                                        Just val@(_, _) → do
                                            when (bc64 globalSettings /= val) . logWith LogInfo $
                                                "bit cost 64 state set to " <> show val <> "\n"
                                            pure (globalSettings{bc64 = val}, processedData)
                        "bcgt64" →
                            let (noChangeString, changeString) = changingStrings firstOption
                                noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                                changeMaybe = readMaybe changeString ∷ Maybe Double
                            in  do
                                    when (length commandList /= length optionList) . failWithPhase Parsing $
                                        "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                    when (null $ firstOption) . failWithPhase Parsing $
                                        "Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): no values found "
                                    when (',' `notElem` firstOption) . failWithPhase Parsing $
                                        "Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): no comma found "
                                    case liftA2 (,) noChangeMaybe changeMaybe of
                                        Nothing →
                                            failWithPhase Parsing $
                                                "Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): "
                                                    <> firstOption
                                        Just val@(_, _) → do
                                            when (bcgt64 globalSettings /= val) . logWith LogInfo $
                                                "bit cost > 64 state set to " <> show val <> "\n"
                                            pure (globalSettings{bcgt64 = val}, processedData)

        
                        val → do
                            when (val `notElem` VER.setArgList) . logWith LogWarn $
                                fold ["Warning: Unrecognized/missing 'set' option '", val, "' in ", show argList, "\n"]
                            pure (globalSettings, processedData)

                    -- =-=-=-=-=-=-=-=-=-=-=-=-=
                    -- =                       =
                    -- = regular command stuff =
                    -- = not initial at start  =
                    -- =                       =
                    -- =-=-=-=-=-=-=-=-=-=-=-=-=
                    _ → case firstCommand of
                        "bc2" →
                            let (noChangeString, changeString) = changingStrings firstOption
                                noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                                changeMaybe = readMaybe changeString ∷ Maybe Double
                            in  do
                                    when (length commandList /= length optionList) . failWithPhase Parsing $
                                        "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                    when (null $ firstOption) . failWithPhase Parsing $
                                        "Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): no values found "
                                    when (',' `notElem` firstOption) . failWithPhase Parsing $
                                        "Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): no comma found "
                                    case liftA2 (,) noChangeMaybe changeMaybe of
                                        Nothing →
                                            failWithPhase Parsing $
                                                "Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): "
                                                    <> firstOption
                                        Just val@(_, _) → do
                                            when (bc2 globalSettings /= val) . logWith LogInfo $
                                                "bit cost 2 state set to " <> show val <> "\n"
                                            pure (globalSettings{bc2 = val}, processedData)
                        "bc4" →
                            let (noChangeString, changeString) = changingStrings firstOption
                                noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                                changeMaybe = readMaybe changeString ∷ Maybe Double
                            in  do
                                    when (length commandList /= length optionList) . failWithPhase Parsing $
                                        "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                    when (null $ firstOption) . failWithPhase Parsing $
                                        "Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): no values found "
                                    when (',' `notElem` firstOption) . failWithPhase Parsing $
                                        "Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): no comma found "
                                    case liftA2 (,) noChangeMaybe changeMaybe of
                                        Nothing →
                                            failWithPhase Parsing $
                                                "Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): "
                                                    <> firstOption
                                        Just val@(_, _) → do
                                            when (bc4 globalSettings /= val) . logWith LogInfo $
                                                "bit cost 4 state set to " <> show val <> "\n"
                                            pure (globalSettings{bc4 = val}, processedData)
                        "bc5" →
                            let (noChangeString, changeString) = changingStrings firstOption
                                noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                                changeMaybe = readMaybe changeString ∷ Maybe Double
                            in  do
                                    when (length commandList /= length optionList) . failWithPhase Parsing $
                                        "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                    when (null $ firstOption) . failWithPhase Parsing $
                                        "Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): no values found "
                                    when (',' `notElem` firstOption) . failWithPhase Parsing $
                                        "Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): no comma found "
                                    case liftA2 (,) noChangeMaybe changeMaybe of
                                        Nothing →
                                            failWithPhase Parsing $
                                                "Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): "
                                                    <> firstOption
                                        Just val@(_, _) → do
                                            when (bc5 globalSettings /= val) . logWith LogInfo $
                                                "bit cost 5 state set to " <> show val <> "\n"
                                            pure (globalSettings{bc5 = val}, processedData)
                        "bc8" →
                            let (noChangeString, changeString) = changingStrings firstOption
                                noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                                changeMaybe = readMaybe changeString ∷ Maybe Double
                            in  do
                                    when (length commandList /= length optionList) . failWithPhase Parsing $
                                        "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                    when (null $ firstOption) . failWithPhase Parsing $
                                        "Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): no values found "
                                    when (',' `notElem` firstOption) . failWithPhase Parsing $
                                        "Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): no comma found "
                                    case liftA2 (,) noChangeMaybe changeMaybe of
                                        Nothing →
                                            failWithPhase Parsing $
                                                "Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): "
                                                    <> firstOption
                                        Just val@(_, _) → do
                                            when (bc8 globalSettings /= val) . logWith LogInfo $
                                                "bit cost 8 state set to " <> show val <> "\n"
                                            pure (globalSettings{bc8 = val}, processedData)
                        "bc64" →
                            let (noChangeString, changeString) = changingStrings firstOption
                                noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                                changeMaybe = readMaybe changeString ∷ Maybe Double
                            in  do
                                    when (length commandList /= length optionList) . failWithPhase Parsing $
                                        "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                    when (null $ firstOption) . failWithPhase Parsing $
                                        "Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): no values found "
                                    when (',' `notElem` firstOption) . failWithPhase Parsing $
                                        "Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): no comma found "
                                    case liftA2 (,) noChangeMaybe changeMaybe of
                                        Nothing →
                                            failWithPhase Parsing $
                                                "Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): "
                                                    <> firstOption
                                        Just val@(_, _) → do
                                            when (bc64 globalSettings /= val) . logWith LogInfo $
                                                "bit cost 64 state set to " <> show val <> "\n"
                                            pure (globalSettings{bc64 = val}, processedData)
                        "bcgt64" →
                            let (noChangeString, changeString) = changingStrings firstOption
                                noChangeMaybe = readMaybe noChangeString ∷ Maybe Double
                                changeMaybe = readMaybe changeString ∷ Maybe Double
                            in  do
                                    when (length commandList /= length optionList) . failWithPhase Parsing $
                                        "Set option error: number of values and options do not match:" <> show (commandList, optionList)
                                    when (null $ firstOption) . failWithPhase Parsing $
                                        "Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): no values found "
                                    when (',' `notElem` firstOption) . failWithPhase Parsing $
                                        "Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): no comma found "
                                    case liftA2 (,) noChangeMaybe changeMaybe of
                                        Nothing →
                                            failWithPhase Parsing $
                                                "Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): "
                                                    <> firstOption
                                        Just val@(_, _) → do
                                            when (bcgt64 globalSettings /= val) . logWith LogInfo $
                                                "bit cost > 64 state set to " <> show val <> "\n"
                                            pure (globalSettings{bcgt64 = val}, processedData)

                        -- processed above, but need here since put in different value
                        "criterion" → do
                            localCriterion ← case firstOption of
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
                                -- NCM | origProcessedData /= emptyProcessedData →
                                --    pure (IL.repeat (0.0, 0.0), U.calculateNCMRootCost origProcessedData)

                                NCM → pure (IL.repeat (0.0, 0.0), U.calculateNCMRootCost origProcessedData)
                                Parsimony → pure $ (IL.repeat (0.0, 0.0), 0.0)
                                MAPA → pure $ (IL.repeat (0.0, 0.0), U.calculateMAPARootCost origProcessedData)
                                val
                                    | val `elem` [PMDL, SI] →
                                        pure $ (U.calculateGraphComplexity &&& (U.calculatePMDLRootCost True)) origProcessedData
                                val → failWithPhase Parsing $ "Optimality criterion not recognized: " <> show val

                            let lGraphFactor
                                    | localCriterion `elem` [PMDL, SI, MAPA] = PMDLGraph
                                    | otherwise = graphFactor globalSettings

                            logWith LogInfo $ case localCriterion of
                                NCM → unwords ["Optimality criterion set to", show NCM, "in -log (base 10) likelihood units\n"]
                                val
                                    | val `elem` [PMDL, SI] →
                                        unwords ["Optimality criterion set to", show val, "Tree Complexity =", show . fst $ IL.head lGraphComplexityList, "bits\n"]
                                val → "Optimality criterion set to " <> show val <> "\n"

                            pure $
                                ( globalSettings
                                    { optimalityCriterion = localCriterion
                                    , graphComplexityList = lGraphComplexityList
                                    , rootComplexity = lRootComplexity
                                    , graphFactor = lGraphFactor
                                    }
                                , processedData
                                )

                        -- modify the behavior of resolutionCache softwired optimization
                        "compressresolutions" → do
                            localCriterion ← case toLower <$> firstOption of
                                "true" → pure True
                                "false" → pure False
                                val →
                                    failWithPhase Parsing $
                                        "Error in 'set' command. CompressResolutions '" <> val <> "' is not 'true' or 'false'"
                            logWith LogInfo $ "CompressResolutions set to " <> show localCriterion <> "\n"
                            pure (globalSettings{compressResolutions = localCriterion}, processedData)

                        -- this not intended to be for users
                        "dynamicepsilon" → case readMaybe (firstOption) ∷ Maybe Double of
                            Nothing →
                                failWithPhase Parsing $
                                    "Set option 'dynamicEpsilon' must be set to a double value >= 0.0 (e.g. dynamicepsilon:0.02): " <> firstOption
                            Just val
                                | val < 0.0 →
                                    failWithPhase Parsing $
                                        "Set option 'dynamicEpsilon' must be set to a double value >= 0.0 (e.g. dynamicepsilon:0.02): " <> show val
                            Just localValue → do
                                logWith LogInfo $ "Dynamic Epsilon factor set to " <> firstOption <> "\n"
                                pure (globalSettings{dynamicEpsilon = 1.0 + (localValue * fractionDynamic globalSettings)}, processedData)
                        "finalassignment" → do
                            localMethod ← case firstOption of
                                "do" → pure DirectOptimization
                                "directoptimization" → pure DirectOptimization
                                "ia" → pure ImpliedAlignment
                                "impliedalignment" → pure ImpliedAlignment
                                val →
                                    failWithPhase Parsing $
                                        fold
                                            [ "Error in 'set' command. FinalAssignment  '"
                                            , val
                                            , "' is not 'DirectOptimization (DO)' or 'ImpliedAlignment (IA)'"
                                            ]

                            case graphType globalSettings of
                                Tree → do
                                    logWith LogInfo $ "FinalAssignment set to " <> show localMethod <> "\n"
                                    pure (globalSettings{finalAssignment = localMethod}, processedData)
                                _ → do
                                    unless (localMethod == DirectOptimization) $
                                        logWith LogInfo "FinalAssignment set to DO (ignoring IA option) for non-Tree graphs\n"
                                    pure (globalSettings{finalAssignment = DirectOptimization}, processedData)
                        "graphfactor" → do
                            localMethod ← case toLower <$> firstOption of
                                "nopenalty" → pure NoNetworkPenalty
                                "w15" → pure Wheeler2015Network
                                "w23" → pure Wheeler2023Network
                                "pmdl" → pure PMDLGraph
                                val →
                                    failWithPhase Parsing $
                                        "Error in 'set' command. GraphFactor  '" <> val <> "' is not 'NoPenalty', 'W15', 'W23', or 'PMDL'"
                            logWith LogInfo $ "GraphFactor set to " <> show localMethod <> "\n"
                            pure (globalSettings{graphFactor = localMethod}, processedData)
                        "graphssteepest" → case readMaybe (firstOption) ∷ Maybe Int of
                            Nothing →
                                failWithPhase Parsing $
                                    "Set option 'graphsSteepest' must be set to an integer value (e.g. graphsSteepest:5): " <> firstOption
                            Just localValue → do
                                logWith LogInfo $ "GraphsStreepest set to " <> show localValue <> "\n"
                                pure (globalSettings{graphsSteepest = localValue}, processedData)
                        "graphtype" → do
                            localGraphType ← case firstOption of
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
                                        ["Graphtype set to", show localGraphType, "with graph factor NoPenalty and final assignment to DO\n"]

                            pure (settingResult, processedData)

                        "keep" -> case readMaybe (firstOption) ∷ Maybe Int of
                            Nothing →
                                failWithPhase Parsing $
                                    "Set option 'keep' must be set to a integer value (e.g. keep:20): " <> firstOption
                            Just localValue → do
                                logWith LogInfo $ "Maximum number of graphs (keep) set to " <> firstOption <> "\n"
                                pure (globalSettings{keepGraphs = localValue}, processedData)


                        -- In first to do stuff above also
                        "missingthreshold" → case readMaybe (firstOption) ∷ Maybe Int of
                            Nothing →
                                failWithPhase Parsing $
                                    "Set option 'missingThreshold' must be set to an integer value (e.g. missingThreshold:50): " <> firstOption
                            Just localValue | localValue == missingThreshold globalSettings → pure (globalSettings, processedData)
                            Just localValue → do
                                logWith LogWarn $ "MissingThreshold set to " <> show localValue <> "\n"
                                pure (globalSettings{missingThreshold = localValue}, processedData)
                        "modelcomplexity" → case readMaybe (firstOption) ∷ Maybe Double of
                            Nothing →
                                failWithPhase Parsing $
                                    "Set option 'modelComplexity' must be set to a double value (e.g. modelComplexity:123.456): " <> firstOption
                            Just localValue → do
                                logWith LogInfo $ "Model Complexity set to " <> firstOption <> "\n"
                                pure (globalSettings{modelComplexity = localValue}, processedData)

                        -- modify the behavior of rerooting character trees for all graph types
                        "multitraverse" → do
                            localCriterion ← case toLower <$> firstOption of
                                "true" → pure True
                                "false" → pure False
                                _ →
                                    failWithPhase Parsing $
                                        "Error in 'set' command. MultiTraverse '" <> firstOption <> "' is not 'true' or 'false'"
                            logWith LogInfo $ "MultiTraverse set to " <> show localCriterion <> "\n"
                            pure (globalSettings{multiTraverseCharacters = localCriterion}, processedData)
                        "outgroup" →
                            let outTaxonName = T.pack $ filter (/= '"') $ head $ filter (/= "") $ fmap snd argList
                            in  case V.elemIndex outTaxonName leafNameVect of
                                    Nothing →
                                        failWithPhase Parsing $
                                            unwords
                                                ["Error in 'set' command. Out-taxon", T.unpack outTaxonName, "not found in input leaf list", show $ T.unpack <$> leafNameVect]
                                    Just outTaxonIndex → do
                                        logWith LogInfo $ "Outgroup set to " <> T.unpack outTaxonName <> "\n"
                                        pure (globalSettings{outgroupIndex = outTaxonIndex, outGroupName = outTaxonName}, processedData)
                        "partitioncharacter" → case firstOption of
                            localPartitionChar@[_] → do
                                when (localPartitionChar /= partitionCharacter globalSettings) . logWith LogInfo $
                                    "PartitionCharacter set to '" <> firstOption <> "'"
                                pure (globalSettings{partitionCharacter = localPartitionChar}, processedData)
                            val →
                                failWithPhase Parsing $
                                    "Error in 'set' command. Partitioncharacter '" <> val <> "' must be a single character"
                        "reportheuristics" → do
                            localMethod ← case toLower <$> firstOption of
                                "true" → pure True
                                "false" → pure False
                                val → failWithPhase Parsing $ "Error in 'set' command. ReportHeursitics  '" <> val <> "' is not 'True' or 'False'"
                            logWith LogInfo $ "ReportHeursitics set to " <> show localMethod <> "\n"
                            pure (globalSettings{reportHeuristics = localMethod}, processedData)
                        "reportnaivedata" → do
                            localMethod ← case toLower <$> firstOption of
                                "true" → pure True
                                "false" → pure False
                                val → failWithPhase Parsing $ "Error in 'set' command. ReportNaive  '" <> val <> "' is not 'True' or 'False'"
                            logWith LogInfo $ "ReportNaiveData set to " <> show localMethod <> "\n"
                            pure (globalSettings{reportNaiveData = localMethod}, processedData)
                        "rootcost" → do
                            localMethod ← case toLower <$> firstOption of
                                "mapa" → pure MAPARoot
                                "ncm" → pure NCMRoot
                                "norootcost" → pure NoRootCost
                                "pmdl" → pure PMDLRoot
                                "si" → pure SIRoot
                                "w15" → pure Wheeler2015Root
                                val → failWithPhase Parsing $ "Error in 'set' command. RootCost '" <> val <> "' is not 'NoRootCost', 'W15', or 'PMDL'"

                            lRootComplexity ← case localMethod of
                                NoRootCost → pure 0.0
                                val | val `elem` [Wheeler2015Root, PMDLRoot] → pure $ U.calculatePMDLRootCost True origProcessedData
                                val → failWithPhase Parsing $ "Error in 'set' command. No determined root complexity of '" <> show val <> "'"

                            logWith LogInfo $ unwords ["RootCost set to", show localMethod, show lRootComplexity, "bits\n"]
                            pure (globalSettings{rootCost = localMethod, rootComplexity = lRootComplexity}, processedData)
                        "seed" → case readMaybe firstOption ∷ Maybe Int of
                            Nothing → failWithPhase Parsing $ "Set option 'seed' must be set to an integer value (e.g. seed:123): " <> firstOption
                            Just localValue → do
                                logWith LogInfo $ "Random Seed set to " <> firstOption <> "\n"
                                setRandomSeed localValue
                                pure (globalSettings, processedData)
                        "softwiredmethod" → do
                            localMethod ← case firstOption of
                                "naive" → pure Naive
                                "exhaustive" → pure Naive
                                "resolutioncache" → pure ResolutionCache
                                _ →
                                    failWithPhase Parsing $
                                        fold ["Error in 'set' command. SoftwiredMethod '", firstOption, "' is not 'Exhaustive' or 'ResolutionCache'"]

                            logWith LogInfo $ "SoftwiredMethod " <> show localMethod <> "\n"
                            pure (globalSettings{softWiredMethod = localMethod}, processedData)

                        -- modify the use of Network Add heurisitcs in network optimization
                        "usenetaddheuristic" → do
                            localCriterion ← case firstOption of
                                "true" → pure True
                                "false" → pure False
                                val →
                                    failWithPhase Parsing $
                                        "Error in 'set' command. UseNetAddHeuristic '" <> val <> "' is not 'true' or 'false'"
                            logWith LogInfo $ "UseNetAddHeuristic set to " <> show localCriterion <> "\n"
                            pure (globalSettings{useNetAddHeuristic = localCriterion}, processedData)

                        -- these not intended for users
                        "jointhreshold" → case readMaybe (firstOption) ∷ Maybe Double of
                            Nothing →
                                failWithPhase Parsing $
                                    "Set option 'joinThreshold' must be set to an double value >= 1.0 (e.g. joinThreshold:1.17): " <> firstOption
                            Just localValue
                                | localValue < 1.0 →
                                    failWithPhase Parsing $
                                        "Set option 'joinThreshold' must be set to a double value >= 1.0 (e.g. joinThreshold:1.17): " <> show localValue
                            Just localValue → do
                                logWith LogInfo $ "JoinThreshold set to " <> show localValue <> "\n"
                                pure (globalSettings{unionThreshold = localValue}, processedData)

                        -- parallel strategy settings options
                        "defparstrat" → do
                            localMethod ← case firstOption of
                                "r0" → pure R0
                                "rpar" → pure RPar
                                "rseq" → pure RSeq
                                "rdeepseq" → pure RDeepSeq
                                val →
                                    failWithPhase Parsing $
                                        "Error in 'set' command. DefParStrat  '" <> val <> "' is not 'r0', 'WrPar', 'rSeq', or 'rDeepSeq'"
                            logWith LogInfo $ "DefParStrat set to " <> show localMethod <> "\n"
                            pure (globalSettings{defaultParStrat = localMethod}, processedData)
                        "lazyparstrat" → do
                            localMethod ← case firstOption of
                                "r0" → pure R0
                                "rpar" → pure RPar
                                "rseq" → pure RSeq
                                "rdeepseq" → pure RDeepSeq
                                val →
                                    failWithPhase Parsing $
                                        "Error in 'set' command. DefParStrat  '" <> val <> "' is not 'r0', 'WrPar', 'rSeq', or 'rDeepSeq'"
                            logWith LogInfo $ "LazyParStrat set to " <> show localMethod <> "\n"
                            pure (globalSettings{lazyParStrat = localMethod}, processedData)
                        "strictparstrat" → do
                            localMethod ← case firstOption of
                                "r0" → pure R0
                                "rpar" → pure RPar
                                "rseq" → pure RSeq
                                "rdeepseq" → pure RDeepSeq
                                val →
                                    failWithPhase Parsing $
                                        "Error in 'set' command. DefParStrat  '" <> val <> "' is not 'r0', 'WrPar', 'rSeq', or 'rDeepSeq'"
                            logWith LogInfo $ "StrictParStrat set to " <> show localMethod <> "\n"
                            pure (globalSettings{strictParStrat = localMethod}, processedData)

                        -- modify the use of implied alkignemnt in heuristics
                        "useia" → do
                            localCriterion ← case firstOption of
                                "true" → pure True
                                "false" → pure False
                                val →
                                    failWithPhase Parsing $
                                        "Error in 'set' command. UseIA '" <> val <> "' is not 'true' or 'false'"
                            logWith LogInfo $ "UseIA set to " <> show localCriterion <> "\n"
                            pure (globalSettings{useIA = localCriterion}, processedData)
                        val → do
                            logWith LogWarn $ fold ["Warning: Unrecognized/missing 'set' option '", val, "' in ", show argList, "\n"]
                            pure (globalSettings, processedData)


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
    → PhyG (String, String, String)
reportCommand globalSettings argList excludeRename numInputFiles crossReferenceString processedData curGraphs supportGraphs =
    let argListWithoutReconcileCommands = filter ((`notElem` VER.reconcileArgList) . fst) argList
        -- check for balances double quotes and only one pair
        outFileNameList = filter (`notElem` ["", "min", "max", "mid"]) $ fmap snd argListWithoutReconcileCommands -- argList
        edgeWeigthVal = 
            let valList = filter (`elem` ["min", "max", "mid"]) $ fmap snd argListWithoutReconcileCommands 
            in
            if null valList then "min"
            else head valList
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
                        else -- This for edge complexities 
                            if ("complexity" `elem` commandList) && (optimalityCriterion globalSettings `notElem` [PMDL, SI])  
                                then do
                                    logWith LogInfo "Cannot report edge complexities unless optimality criterion is either PMDL or SI\n"
                                    pure ("Cannot report edge complexities unless optimality criterion is either PMDL or SI", outfileName, writeMode)
                            else if ("complexity" `elem` commandList) && (optimalityCriterion globalSettings `elem` [PMDL, SI])
                                then -- else if (not .null) (L.intersect ["graphs", "newick", "dot", "dotpdf"] commandList) then
                                    let relabelledDecoratedGraph = fmap (relabelEdgeComplexity globalSettings processedData ) (fmap thd5 curGraphs)
                                        graphString = outputGraphString edgeWeigthVal commandList (outgroupIndex globalSettings) relabelledDecoratedGraph (fmap snd5 curGraphs)
                                    in  if null curGraphs
                                        then do
                                            logWith LogInfo "No graphs to report edge compleities \n"
                                            pure ("No graphs to edge complexities report", outfileName, writeMode)
                                        else do
                                            logWith LogInfo ("Reporting " <> show (length curGraphs) <> " edge complexities of graphs at minimum cost " <> show (minimum $ fmap snd5 curGraphs) <> "\n")
                                            pure (graphString, outfileName, writeMode)
                            else if "crossrefs" `elem` commandList
                                then
                                    let dataString = crossReferenceString
                                    in  pure (dataString, outfileName, writeMode)
                                else
                                    if "data" `elem` commandList
                                        then
                                            let blocks = thd3 processedData
                                                numChars = V.sum $ fmap (V.length . thd3) blocks
                                                dataString = phyloDataToString 0 blocks
                                                -- \$ thd3 processedData
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
                                                    curGraphs' ←
                                                        if False -- not (reportNaiveData globalSettings)
                                                            then pure curGraphs
                                                            else
                                                                let action ∷ SimpleGraph → PhyG ReducedPhylogeneticGraph
                                                                    action = TRAV.multiTraverseFullyLabelGraphReduced globalSettings processedData False False Nothing
                                                                in  getParallelChunkTraverse >>= \pTraverse →
                                                                        (action . fst5) `pTraverse` curGraphs

                                                    dataStringList <- 
                                                        let action :: (ReducedPhylogeneticGraph, Int) → PhyG [[String]]
                                                            action = getGraphDiagnosis globalSettings processedData 
                                                        in do
                                                            diagPar ← getParallelChunkTraverse
                                                            diagPar action (zip curGraphs' [0 .. (length curGraphs' - 1)])

                                                    let dataString = CSV.genCsvFile $ concat dataStringList
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

                                                                    let graphString = outputGraphString edgeWeigthVal commandList (outgroupIndex globalSettings) (fmap thd5 curGraphs) (fmap snd5 curGraphs)
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
                                                                            if "metadata" `elem` commandList
                                                                                then do
                                                                                    curGraphs' ←
                                                                                        if False -- not (reportNaiveData globalSettings)
                                                                                            then pure curGraphs
                                                                                            else
                                                                                                let action ∷ SimpleGraph → PhyG ReducedPhylogeneticGraph
                                                                                                    action = TRAV.multiTraverseFullyLabelGraphReduced globalSettings processedData False False Nothing
                                                                                                in  getParallelChunkTraverse >>= \pTraverse →
                                                                                                        (action . fst5) `pTraverse` curGraphs

                                                                                    dataStringList <- 
                                                                                        let action :: (ReducedPhylogeneticGraph, Int) → PhyG [[String]]
                                                                                            action = getGraphMetaData globalSettings processedData 
                                                                                        in do
                                                                                            diagPar ← getParallelChunkTraverse
                                                                                            diagPar action (zip curGraphs' [0 .. (length curGraphs' - 1)])

                                                                                    let dataString = CSV.genCsvFile $ concat dataStringList
                                                                                    if null curGraphs
                                                                                        then do
                                                                                            logWith LogInfo "No graphs to get metaData\n"
                                                                                            pure ("No graphs to get metaData", outfileName, writeMode)
                                                                                        else do
                                                                                            logWith
                                                                                                LogInfo
                                                                                                ("Getting metaData from " <> show (length curGraphs) <> " graphs at minimum cost " <> show (minimum $ fmap snd5 curGraphs) <> "\n")
                                                                                            pure (dataString, outfileName, writeMode)
                                                                            else if "pairdist" `elem` commandList
                                                                                then
                                                                                    let nameData = L.intercalate "," (V.toList (T.unpack <$> fst3 processedData)) <> "\n"
                                                                                    in  do
                                                                                            pairwiseDistanceMatrix' ← D.getPairwiseDistances processedData
                                                                                            let dataString = CSV.genCsvFile $ fmap (fmap show) pairwiseDistanceMatrix'
                                                                                            pure (nameData <> dataString, outfileName, writeMode)
                                                                                else
                                                                                    if "parameterestimation"  `elem` commandList
                                                                                        then do
                                                                                            curGraphs' ←
                                                                                                if False -- not (reportNaiveData globalSettings)
                                                                                                    then pure curGraphs
                                                                                                    else
                                                                                                        let action ∷ SimpleGraph → PhyG ReducedPhylogeneticGraph
                                                                                                            action = TRAV.multiTraverseFullyLabelGraphReduced globalSettings processedData False False Nothing
                                                                                                        in  getParallelChunkTraverse >>= \pTraverse →
                                                                                                                (action . fst5) `pTraverse` curGraphs

                                                                                            dataStringList <- 
                                                                                                let action :: (ReducedPhylogeneticGraph, Int) → PhyG [[String]]
                                                                                                    action = getGraphParameters globalSettings processedData 
                                                                                                in do
                                                                                                    diagPar ← getParallelChunkTraverse
                                                                                                    diagPar action (zip curGraphs' [0 .. (length curGraphs' - 1)])

                                                                                            let dataString = CSV.genCsvFile $ concat dataStringList
                                                                                            
                                                                                            if null curGraphs
                                                                                                then do
                                                                                                    logWith LogInfo "No graphs for metaData\n"
                                                                                                    pure ("No graphs for metaData", outfileName, writeMode)
                                                                                                else do
                                                                                                    logWith
                                                                                                        LogInfo
                                                                                                        ("Getting metaData for " <> show (length curGraphs) <> " graphs at minimum cost " <> show (minimum $ fmap snd5 curGraphs) <> "\n")
                                                                                                    pure (dataString, outfileName, writeMode)

                                                                                    else if "reconcile" `elem` commandList
                                                                                        then do
                                                                                            recResult ← R.makeReconcileGraph VER.reconcileArgList argList (fmap fst5 curGraphs)
                                                                                            -- let (reconcileString, ) = recResult
                                                                                            let (_, reconcileGraph) = recResult
                                                                                            let reconcileString = outputGraphString edgeWeigthVal commandList (outgroupIndex globalSettings) [GO.convertSimpleToDecoratedGraph reconcileGraph] [0]
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
                                                                                                                            curGraphs' ←
                                                                                                                                if not (reportNaiveData globalSettings)
                                                                                                                                    then pure $ GO.convertReduced2PhylogeneticGraph <$> curGraphs
                                                                                                                                    else
                                                                                                                                        let action ∷ SimpleGraph → PhyG PhylogeneticGraph
                                                                                                                                            action = TRAV.multiTraverseFullyLabelGraph globalSettings processedData False False Nothing
                                                                                                                                        in  getParallelChunkTraverse >>= \pTraverse →
                                                                                                                                                pTraverse (action . fst5) curGraphs

                                                                                                                            tntContentList' ← traverse (getTNTString globalSettings processedData) $ zip curGraphs' [0 .. length curGraphs' - 1]
                                                                                                                            let tntContentList = concat tntContentList'
                                                                                                                            pure (tntContentList, outfileName, writeMode)
                                                                                                                else do
                                                                                                                    logWith LogWarn ("\nUnrecognized/missing report option in " <> show commandList <> " defaulting to 'graphs'" <> "\n")
                                                                                                                    let graphString = outputGraphString edgeWeigthVal commandList (outgroupIndex globalSettings) (fmap thd5 curGraphs) (fmap snd5 curGraphs)
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
        bracketToCurly = \case
            '(' → '{'
            ')' → '}'
            val → val


changingStrings ∷ String → (String, String)
changingStrings str = case span (/= ',') $ filter (`notElem` ['(', ')']) str of
    (prefix, ',' : suffix) → (prefix, suffix)
    (prefix, _) → (prefix, mempty)
