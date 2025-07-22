{-# LANGUAGE ImportQualifiedPost #-}

{- |
Program to perform phylogenetic searches on general graphs with diverse data types.
-}
module Main (main) where

import CommandLineOptions
import Commands.CommandExecution qualified as CE
import Commands.ProcessCommands qualified as PC
import Commands.Verify qualified as V
import Control.Monad (when)
import Control.Monad.IO.Class (MonadIO (..))
import Data.CSV qualified as CSV
import Data.Foldable (fold)
import Data.Foldable1 (head)
import Data.InfList qualified as IL
import Data.List qualified as L
import Data.List.NonEmpty (NonEmpty (..))
import Data.String (fromString)
import Data.Text.Builder.Linear (runBuilder)
import Data.Text.Lazy qualified as Text
import Data.Text.Short qualified as ST
import Data.Time.Clock
import GeneralUtilities
import GraphFormatUtilities qualified as GFU
import GraphOptimization.Traversals qualified as T
import Graphs.GraphOperations qualified as GO
import Input.BitPack qualified as BP
import Input.DataTransformation qualified as DT
import Input.ReadInputFiles qualified as RIF
import Input.Reorganize qualified as R
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import PHANE.Evaluation.Verbosity (Verbosity (..))
import Software.Preamble
import System.CPUTime
import System.IO
import Types.Types
import Utilities.LocalGraph qualified as LG
import Utilities.Utilities qualified as U
import Prelude hiding (head)


{- |
Main entry point.
-}
main ∷ IO ()
main = do
    hSetEncoding stdout utf8
    hSetEncoding stderr utf8
    opts ← getArgsCLI <$> parseCommandLineOptions
    either printInformationDisplay evaluatePhyG opts


{- |
Fully evaluate the 'PhyG' computation 'performSearch' by encapsulating effects
within the 'Evaluation' monad in order to correctly control logging, parallelism,
and randomness.
-}
evaluatePhyG ∷ FilePath → IO ()
evaluatePhyG path = do
    logConfig ← initializeLogging Dump Warn Nothing
    firstSeed ← initializeRandomSeed
    runEvaluation logConfig firstSeed () $ performSearch firstSeed path


{- |
Perform phylogenetic search using the supplied input file.
-}
performSearch ∷ RandomSeed → FilePath → PhyG ()
performSearch initialSeed inputFilePath = do
    printProgramPreamble
    logWith LogInfo $ "\nCommand script file: '" <> inputFilePath <> "'\n"
    logWith LogInfo $ "Initial random seed set to " <> show initialSeed <> "\n"
    timeCDBegin ← liftIO getCurrentTime

    -- Process commands to get list of actions
    commandContents ← liftIO $ readFile inputFilePath

    -- Process run commands to create one list of things to do
    commandContents' ← PC.expandRunCommands [] (lines commandContents)
    thingsToDo'' ← PC.getCommandList commandContents'
    -- let thingsToDo'' = PC.getCommandList commandContents'
    -- mapM_ (logWith LogTech . show) thingsToDo'

    -- preprocess commands for non-parsimony optimality criteria
    let thingsToDo' = PC.preprocessOptimalityCriteriaScripts thingsToDo''

    -- Process Read commands (with prealigned and tcm flags first)
    -- expand read commands for wildcards
    expandedReadCommands ← mapM (RIF.expandReadCommands []) $ filter ((== Read) . fst) thingsToDo'

    -- sort added to sort input read commands for left right consistency
    let thingsToDo = PC.moveSetOutgroupFirst $ L.sort (fold expandedReadCommands) <> filter ((/= Read) . fst) thingsToDo'
    -- logWith LogDump . show $ fold expandedReadCommands

    -- check commands and options for basic correctness
    logWith LogMore "\nChecking command file syntax\n"
    let !commandsOK = V.verifyCommands thingsToDo [] []

    if commandsOK
        then logWith LogMore "\tCommands appear to be properly specified--file availability and contents not checked.\n"
        else failWithPhase Parsing "Commands not properly specified\n"

    movedPrealignedList ← mapM (PC.movePrealignedTCM . snd) (filter ((== Read) . fst) thingsToDo)
    dataGraphList ← mapM RIF.executeReadCommands movedPrealignedList
    let (rawData, rawGraphs, terminalsToInclude, terminalsToExclude, renameFilePairs, reBlockPairs) = RIF.extractInputTuple dataGraphList

    if null rawData && null rawGraphs
        then failWithPhase Unifying "\n\nNeither data nor graphs entered.  Nothing can be done.\n"
        else logWith LogInfo $ unwords ["Entered", show (length rawData), "data file(s) and", show (length rawGraphs), "input graphs\n"]

    -- get set partitions character from Set commands early, the isFirst==True puts in first section--only processing a few fields
    -- confusing and should be changed
    let setCommands = filter ((== Set) . fst) thingsToDo
    (_, partitionCharOptimalityGlobalSettings, _) ←
        CE.executeCommands emptyGlobalSettings mempty 0 [] mempty mempty mempty mempty mempty setCommands True

    -- Split fasta/fastc sequences into corresponding pieces based on '#' partition character
    rawDataSplit ← DT.partitionSequences (ST.fromString (partitionCharacter partitionCharOptimalityGlobalSettings)) rawData

    -- Process Rename Commands
    newNamePairList ← liftIO $ CE.executeRenameReblockCommands Rename renameFilePairs thingsToDo
    if thereExistsSome newNamePairList
        then logWith LogInfo $ unwords ["Renaming", show $ length newNamePairList, "terminals\n"]
        else logWith LogInfo "No terminals to be renamed\n"

    let renamedData = fmap (DT.renameData newNamePairList) rawDataSplit
    let renamedGraphs = fmap (GFU.relabelGraphLeaves newNamePairList) rawGraphs

    let numInputFiles = length renamedData

    let thingsToDoAfterReadRename = filter ((/= Read) . fst) $ filter ((/= Rename) . fst) thingsToDo

    -- Reconcile Data and Graphs (if input) including ladderization
    -- could be sorted, but no real need
    -- get taxa to include in analysis
    let renderTerminals pref = (<> "\n") . (("Terminals to " <> pref <> ": ") <>) . unwords . fmap Text.unpack
    when (thereExistsSome terminalsToInclude)
        . logWith LogInfo
        $ renderTerminals "include" terminalsToInclude

    when (thereExistsSome terminalsToExclude)
        . logWith LogInfo
        $ renderTerminals "exclude" terminalsToExclude

    -- Uses names from terminal list if non-null, and remove excluded terminals
    let dataLeafNames' =
            if thereExistsSome terminalsToInclude
                then L.sort $ L.nub terminalsToInclude
                else L.sort $ DT.getDataTerminalNames renamedData
    let dataLeafNames'' = dataLeafNames' L.\\ terminalsToExclude
    logWith LogInfo ("Data were input for " <> show (length dataLeafNames'') <> " terminals\n")

    -- check data for missing data threshold and remove those above
    let missingToExclude = DT.checkLeafMissingData (missingThreshold partitionCharOptimalityGlobalSettings) rawData
    let dataLeafNames =
            if (missingThreshold partitionCharOptimalityGlobalSettings) == 100
                then dataLeafNames''
                else -- get number of occurrences of name in rawData Terminfo lists
                    dataLeafNames'' L.\\ missingToExclude

    when (thereExistsSome missingToExclude)
        . logWith LogInfo
        $ unlines
            [ "Terminals above missing data threshold and excluded: "
            , unwords [show (length dataLeafNames'' - length dataLeafNames), show missingToExclude]
            , unwords [show $ length dataLeafNames, "terminals remain to be analyzed"]
            ]

    (crossReferenceString, defaultGlobalSettings, naiveData, reconciledData, reconciledGraphs) ← case dataLeafNames of
        [] → failWithPhase Unifying "No leaf data to be analyzed--all excluded"
        x : xs →
            {-
            Data processing here-- there are multiple steps not composed so that
            large data files can be precessed and intermediate data goes out
            of scope and can be freed back to system.

            This created here and passed to command execution later to
            remove dependency of renamed data in command execution to
            reduce memory footprint keeoing that stuff around.
            -}
            let dNames = x :| xs
                -- crossReferenceString = CSV.genCsvFile $ CE.getDataListList renamedData dNames
                -- Add in missing terminals to raw data where required
                reconciledData' = DT.addMissingTerminalsToInput dNames [] <$> renamedData
                reconciledGraphs = fmap (GFU.reIndexLeavesEdges dNames . GFU.checkGraphsAndData dNames) renamedGraphs
            in  do
                    -- Check for data file with all missing data--as in had no terminals with data in termainals list
                    reconciledData ← fold <$> traverse DT.removeAllMissingCharacters reconciledData'

                    -- Create unique bitvector names for leaf taxa.
                    let leafBitVectorNames = DT.createBVNames reconciledData

                    -- Create naive data
                    -- basic usable format organized into blocks,
                    -- but not grouped by types, or packed (bit, sankoff, prealigned etc)
                    -- Need to check data for equal in character number
                    naiveData ← DT.createNaiveData partitionCharOptimalityGlobalSettings reconciledData leafBitVectorNames []
                    let crossReferenceString = CSV.genCsvFile $ CE.getDataListList' renamedData naiveData

                    {-  for debugging fiole names after partitions
                    let bD = thd3  naiveData
                    let bN = fmap (Text.unpack . fst3)  bD
                    logWith LogInfo $ "NaivRecoded: " <> (concat bN)
                    -}

                    -- get mix of static/dynamic characters to adjust dynmaicEpsilon
                    -- doing on naive data so no packing etc
                    let fractionDynamicData = U.getFractionDynamic naiveData

                    -- Set global values before search--should be integrated with executing commands
                    -- only stuff that is data dependent here (and seed)
                    let defaultGlobalSettings =
                            emptyGlobalSettings
                                { outgroupIndex = 0
                                , outGroupName = head dNames
                                , seed = fromEnum initialSeed
                                , numDataLeaves = length leafBitVectorNames
                                , fractionDynamic = fractionDynamicData
                                , dynamicEpsilon = 1.0 + ((dynamicEpsilon emptyGlobalSettings - 1.0) * fractionDynamicData)
                                }

                    pure (crossReferenceString, defaultGlobalSettings, naiveData, reconciledData, reconciledGraphs)

    -- logWith LogInfo ("Fraction characters that are dynamic: " <> (show $ (fromIntegral lengthDynamicCharacters) / (fromIntegral $ lengthDynamicCharacters + numStaticCharacters)))

    -- Check to see if there are taxa without any observations. Would become total wildcards
    let taxaDataSizeList =
            filter ((== 0) . snd)
                . zip dataLeafNames
                . foldl1 (zipWith (+))
                $ fmap (fmap (snd3 . U.filledDataFields (0, 0)) . fst) reconciledData

    case taxaDataSizeList of
        [] → logWith LogInfo "All taxa contain data\n"
        xs →
            failWithPhase
                Unifying
                $ fold
                    [ "\nError: There are taxa without any data: "
                    , L.intercalate ", " $ Text.unpack . fst <$> xs
                    , "\n"
                    ]

    -- Ladderizes (resolves) input graphs and ensures that networks are time-consistent
    -- chained network nodes should never be introduced later so only checked no
    -- checks for children of tree node that are all netowork nodee (causes displayu problem)
    -- let noChainNetNodesList = fmap fromJust $ filter (/=Nothing) $ fmap (LG.removeChainedNetworkNodes True) reconciledGraphs
    let noSisterNetworkNodes = fmap LG.removeTreeEdgeFromTreeNodeWithAllNetworkChildren reconciledGraphs -- noChainNetNodesList
    ladderizedGraphList ← mapM (GO.convertGeneralGraphToPhylogeneticGraph True) noSisterNetworkNodes

    {-To do
    -- Remove any not "selected" taxa from both data and graphs (easier to remove from fgl)
    let reconciledData' = removeTaxaFromData includeList reconciledData
    let reconciledGraphs' = removeTaxaFromGraphs includeList reconciledData
    -}

    -- Group Data--all nonadditives to single character, additives with
    -- alphabet < 64 recoded to nonadditive binary, additives with same alphabet
    -- combined,
    naiveDataGrouped ← R.combineDataByType partitionCharOptimalityGlobalSettings naiveData -- R.groupDataByType naiveData

    -- Bit pack non-additive data
    naiveDataPacked ← BP.packNonAdditiveData partitionCharOptimalityGlobalSettings naiveDataGrouped

    -- Optimize Data convert
    -- prealigned to non-additive or matrix
    -- bitPack resulting non-additive
    optimizedPrealignedData ← R.optimizePrealignedData partitionCharOptimalityGlobalSettings naiveDataPacked

    -- Execute any 'Block' change commands--make reBlockedNaiveData
    newBlockPairList ← liftIO $ CE.executeRenameReblockCommands Reblock reBlockPairs thingsToDo

    reBlockedNaiveData ← R.reBlockData newBlockPairList optimizedPrealignedData -- naiveData
    let thingsToDoAfterReblock = filter ((/= Reblock) . fst) $ filter ((/= Rename) . fst) thingsToDoAfterReadRename

    -- Combines data of exact types into single vectors in each block
    -- this is final data processing step
    optDataNBPL ← R.combineDataByType partitionCharOptimalityGlobalSettings reBlockedNaiveData

    when (thereExistsSome newBlockPairList) $ logWith LogInfo "Reorganizing Block data"
    let optimizedData =
            if thereExistsSome newBlockPairList
                then optDataNBPL
                else optimizedPrealignedData

    -- set outgroup needs to be the firt set command (I think due to stricness)
    let initialSetCommands = filter ((== Set) . fst) thingsToDoAfterReblock

    let commandsAfterInitialDiagnose = filter ((/= Set) . fst) thingsToDoAfterReblock

    -- Set reporting data for qualitative characters to Naive data (usually but not if huge data set), empty if packed
    reportingData <-
            if reportNaiveData partitionCharOptimalityGlobalSettings then R.reBlockData newBlockPairList naiveData
            else pure emptyProcessedData

    -- Check if no dynamic charcters--if so--then set multitraverse to False
    {-
    let (multiTraverseSetting, softWiredMethodSetting) = if 0 /= (U.getNumberNonExactCharacters $ thd3 reBlockedNaiveData)
                                                                then (multiTraverseCharacters defaultGlobalSettings, softWiredMethod defaultGlobalSettings)
                                                         else (False, ResolutionCache) 
    -}
    let multiTraverseSetting = if 0 /= (U.getNumberNonExactCharacters $ thd3 reBlockedNaiveData)
                                                                then multiTraverseCharacters defaultGlobalSettings
                               else False

    --logWith LogInfo $ "New Vals: " <> (show $ (False, ResolutionCache))

    -- This rather awkward syntax makes sure global settings (outgroup, criterion etc) are in place for initial input graph diagnosis
    (_, initialGlobalSettings, _) ←
        CE.executeCommands
            (defaultGlobalSettings {multiTraverseCharacters = multiTraverseSetting}) -- , softWiredMethod = softWiredMethodSetting})
            (terminalsToExclude, renameFilePairs)
            numInputFiles
            crossReferenceString
            optimizedData
            optimizedData
            reportingData
            []
            []
            initialSetCommands
            False

    -- Get CPUTime so far ()data input and processing
    dataCPUTime ← liftIO getCPUTime

    -- Diagnose any input graphs
    let action = T.multiTraverseFullyLabelGraphReduced initialGlobalSettings optimizedData True True Nothing
    actionPar ← getParallelChunkTraverse

    inputGraphList ←
        actionPar action (fmap (LG.rerootTree (outgroupIndex initialGlobalSettings)) ladderizedGraphList)
    {-PU.seqParMap
        PU.myStrategy
        (T.multiTraverseFullyLabelGraphReduced initialGlobalSettings optimizedData True True Nothing)
        (fmap (LG.rerootTree (outgroupIndex initialGlobalSettings)) ladderizedGraphList)
    -}

    -- Get CPUTime for input graphs
    afterGraphDiagnoseTCPUTime ← liftIO getCPUTime
    let (inputGraphTime, inGraphNumber, minOutCost, maxOutCost) = case inputGraphList of
            [] → (0, 0, infinity, infinity)
            _ →
                ( fromIntegral afterGraphDiagnoseTCPUTime - fromIntegral dataCPUTime
                , length inputGraphList
                , minimum $ fmap snd5 inputGraphList
                , maximum $ fmap snd5 inputGraphList
                )

    let inputProcessingData = emptySearchData{commentString = "Input and data processing", duration = fromIntegral dataCPUTime}
    let inputGraphProcessing =
            emptySearchData
                { minGraphCostOut = minOutCost
                , maxGraphCostOut = maxOutCost
                , numGraphsOut = inGraphNumber
                , commentString = "Input graph processing"
                , duration = inputGraphTime
                }

    -- Execute Following Commands (searches, reports etc)
    (finalGraphList, _, _) ←
        CE.executeCommands
            (initialGlobalSettings{searchData = [inputGraphProcessing, inputProcessingData]})
            (terminalsToExclude, renameFilePairs)
            numInputFiles
            crossReferenceString
            optimizedData
            optimizedData
            reportingData
            inputGraphList
            []
            commandsAfterInitialDiagnose -- (transformString <> commandsAfterInitialDiagnose)
            False

    -- print global setting just to check
    -- logWith LogInfo (show _finalGlobalSettings)

    -- Add in model and root cost if optimality criterion needs it
    -- if (rootComplexity initialGlobalSettings) /= 0.0 then logWith LogInfo ("\tUpdating final graph with any root priors")
    -- else logWith LogInfo ""

    {-
    This should not be necessary--moved to traversal with root cost adjustment
        rediagnose for NCM and PMDL due to packing, in most cases not required, just being sure etc
    -}
    let rediagnoseWithReportingdata = True
    finalGraphList' ←
        if rediagnoseWithReportingdata
            then -- if optimalityCriterion initialGlobalSettings `elem` [SI, NCM, PMDL] then
                T.updateGraphCostsComplexities initialGlobalSettings reportingData optimizedData rediagnoseWithReportingdata finalGraphList
            else pure finalGraphList

    -- let finalGraphList' = finalGraphList

    let minCost = if null finalGraphList then 0.0 else minimum $ fmap snd5 finalGraphList'
    let maxCost = if null finalGraphList then 0.0 else maximum $ fmap snd5 finalGraphList'

    -- get network numbers for graph complexities (PMDL, SI)
    let grabber = length . fth4 . LG.splitVertexList . fst5 <$> finalGraphList'
    let pairFunction ∷ ∀ {a}. (a, a) → a
        (netWorkVertexList, pairFunction, units)
            | optimalityCriterion initialGlobalSettings == Parsimony = (replicate (length finalGraphList') 0, fst, "")
            -- PMDL and DI in base 2
            | optimalityCriterion initialGlobalSettings `elem` [PMDL, SI] = (grabber, fst, " bits")
            -- \| graphType initialGlobalSettings == SoftWired = (grabber, fst, " bits")
            -- NCM and MAPA in base 10
            | otherwise = (grabber, snd, " dits")

    -- final results reporting to stderr
    logWith LogInfo $
        unwords
            [ "Execution returned"
            , show $ length finalGraphList'
            , "graph(s) at cost range"
            , (show (minCost, maxCost)) <> units
            , "\n"
            ]

    -- insures model complexity 0 if not PMDL  correctly accounted for in traversals
    let adjModelComplexity =
            if optimalityCriterion initialGlobalSettings == PMDL
                then modelComplexity initialGlobalSettings
                else 0.0

    when (optimalityCriterion initialGlobalSettings /= Parsimony) $
        logWith LogInfo $
            unwords
                [ "\tModel complexity " <> (show adjModelComplexity) <> units <> "\n"
                , "\tRoot complexity " <> (show $ rootComplexity initialGlobalSettings) <> units <> "\n"
                , "\tGraph complexities "
                    <> (show $ fmap pairFunction $ fmap ((graphComplexityList initialGlobalSettings) IL.!!!) netWorkVertexList)
                    <> units
                , "\n\n"
                ]

    -- Final Stderr report
    timeCPUEnd ← liftIO getCPUTime
    timeCDEnd ← liftIO getCurrentTime

    -- logWith LogInfo ("CPU Time " <> (show timeCPUEnd))
    let wallClockDuration = floor (1000000000000 * nominalDiffTimeToSeconds (diffUTCTime timeCDEnd timeCDBegin)) ∷ Integer
    -- logWith LogInfo ("Current time " <>  (show wallClockDuration))
    let cpuUsage = fromIntegral timeCPUEnd / fromIntegral wallClockDuration ∷ Double
    -- logWith LogInfo ("CPU % " <> (show cpuUsage))

    logWith LogInfo . unlines $
        ("\t" <>)
            <$> [ unwords ["Wall-Clock time ", show ((fromIntegral wallClockDuration ∷ Double) / 1000000000000.0), "second(s)"]
                , unwords ["CPU time", show ((fromIntegral timeCPUEnd ∷ Double) / 1000000000000.0), "second(s)"]
                , unwords ["CPU usage", show (floor (100.0 * cpuUsage) ∷ Integer) <> "%"]
                ]


thereExistsSome ∷ (Foldable f) ⇒ f a → Bool
thereExistsSome = not . null


printProgramPreamble ∷ PhyG ()
printProgramPreamble = liftIO preambleText >>= logWith LogDone . runBuilder . (fromString "\n\n" <>)
