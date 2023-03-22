{- |
Module      :  phygraph.hs
Description :  Progam to perform phylogenetic searchs on general graphs with diverse data types
Copyright   :  (c) 2021 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
License     :

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

Maintainer  :  Ward Wheeler <wheeler@amnh.org>
Stability   :  unstable
Portability :  portable (I hope)

-}

{-# LANGUAGE CPP #-}

module Main (main) where

import qualified Commands.CommandExecution    as CE
import qualified Commands.ProcessCommands     as PC
import qualified Commands.Verify              as V
import qualified Data.CSV                     as CSV
import qualified Data.List                    as L
import           Data.Maybe
import qualified Data.Text.Lazy               as Text
import qualified Data.Text.Short              as ST
import           Data.Time.Clock
import           Debug.Trace
import           GeneralUtilities
import qualified GraphFormatUtilities         as GFU
import qualified GraphOptimization.Traversals as T
import qualified Graphs.GraphOperations       as GO
import qualified Input.BitPack                as BP
import qualified Input.DataTransformation     as DT
import qualified Input.ReadInputFiles         as RIF
import qualified Input.Reorganize             as R
import qualified ParallelUtilities            as PU
import           System.CPUTime
import           System.Environment
import           System.IO
import           Types.Types
import qualified Utilities.Distances          as D
import qualified Utilities.LocalGraph         as LG
import qualified Utilities.Utilities          as U
-- import           GHC.Stats

-- | main driver
main :: IO ()
main = do
    let compileDate = (__DATE__ ++ " " ++ __TIME__)
    -- let compileDate = (__DATE__ ++ " " ++ __TIME__)
    let splash = "\nPhyG version " ++ pgVersion ++ "\nCopyright(C) 2022-2023 Ward Wheeler and The American Museum of Natural History\n"
    let splash2 = "PhyG comes with ABSOLUTELY NO WARRANTY; This is free software, and may be \nredistributed "
    let splash3 = "\tunder the 3-Clause BSD License.\nCompiled " ++ compileDate
    hPutStrLn stderr (splash ++ splash2 ++ splash3)

    -- Process arguments--a single file containing commands
    args <- getArgs

    if length args /= 1 then errorWithoutStackTrace "\nProgram requires a single argument--the name of command script file.\n\n"
    else hPutStr stderr "\nCommand script file: "
    hPutStrLn stderr $ head args

    -- System time for Random seed
    timeD <- getSystemTimeSeconds
    timeCD <- getCurrentTime
    hPutStrLn stderr ("Initial random seed set to " ++ show timeD)

    -- hPutStrLn stderr ("Current time is " ++ show timeD)
    let seedList = randomIntList timeD

    -- Process commands to get list of actions
    commandContents <- readFile $ head args

    -- Process run commands to create one list of things to do
    commandContents' <- PC.expandRunCommands [] (lines commandContents)
    let thingsToDo'' = PC.getCommandList  commandContents'
    --mapM_ (hPutStrLn stderr) (fmap show thingsToDo')

    -- preprocess commands for non-parsimony optimality criteria
    let thingsToDo' = PC.preprocessOptimalityCriteriaScripts thingsToDo''

    -- Process Read commands (with prealigned and tcm flags first)
        --expand read commands for wildcards
    expandedReadCommands <- mapM (RIF.expandReadCommands []) $ filter ((== Read) . fst) thingsToDo'

    -- sort added to sort input read commands for left right consistancy
    let thingsToDo = L.sort (concat expandedReadCommands) ++ filter ((/= Read) . fst) thingsToDo'
    --hPutStrLn stderr (show $ concat expandedReadCommands)

    -- check commands and options for basic correctness
    hPutStrLn stderr "\tChecking command file syntax"
    let !commandsOK = V.verifyCommands thingsToDo [] []

    if commandsOK then hPutStrLn stderr "Commands appear to be properly specified--file availability and contents not checked.\n"
    else errorWithoutStackTrace "Commands not properly specified"

    dataGraphList <- mapM RIF.executeReadCommands $ fmap (PC.movePrealignedTCM . snd) (filter ((== Read) . fst) thingsToDo)
    let (rawData, rawGraphs, terminalsToInclude, terminalsToExclude, renameFilePairs, reBlockPairs) = RIF.extractInputTuple dataGraphList

    if null rawData && null rawGraphs then errorWithoutStackTrace "\n\nNeither data nor graphs entered.  Nothing can be done."
    else hPutStrLn stderr ("Entered " ++ show (length rawData) ++ " data file(s) and " ++ show (length rawGraphs) ++ " input graphs")

    -- get set partitions character from Set commands early, the enpty seed list puts in first section--only processing a few fields
    -- confusing and should be changed
    let setCommands = filter ((== Set).fst) thingsToDo
    (_, partitionCharOptimalityGlobalSettings, _, _) <- CE.executeCommands emptyGlobalSettings mempty 0 [] mempty mempty mempty mempty mempty mempty mempty setCommands

    -- Split fasta/fastc sequences into corresponding pieces based on '#' partition character
    let rawDataSplit = DT.partitionSequences (ST.fromString (partitionCharacter partitionCharOptimalityGlobalSettings)) rawData

    -- Process Rename Commands
    newNamePairList <- CE.executeRenameReblockCommands Rename renameFilePairs thingsToDo
    if not $ null newNamePairList then hPutStrLn stderr ("Renaming " ++ show (length newNamePairList) ++ " terminals")
    else hPutStrLn stderr "No terminals to be renamed"

    let renamedData   = fmap (DT.renameData newNamePairList) rawDataSplit
    let renamedGraphs = fmap (GFU.relabelGraphLeaves  newNamePairList) rawGraphs

    let numInputFiles = length renamedData

    let thingsToDoAfterReadRename = filter ((/= Read) .fst) $ filter ((/= Rename) .fst) thingsToDo

    -- Reconcile Data and Graphs (if input) including ladderization
        -- could be sorted, but no real need
        -- get taxa to include in analysis
    if not $ null terminalsToInclude then hPutStrLn stderr ("Terminals to include:" ++ concatMap (++ " ") (fmap Text.unpack terminalsToInclude))
    else hPutStrLn stderr ""
    if not $ null terminalsToExclude then hPutStrLn stderr ("Terminals to exclude:" ++ concatMap (++ " ") (fmap Text.unpack terminalsToExclude))
    else hPutStrLn stderr ""

    -- Uses names from terminal list if non-null, and remove excluded terminals
    let dataLeafNames' = if not $ null terminalsToInclude then L.sort $ L.nub terminalsToInclude
                        else L.sort $ DT.getDataTerminalNames renamedData
    let dataLeafNames = dataLeafNames' L.\\ terminalsToExclude
    hPutStrLn stderr ("Data were input for " ++ show (length dataLeafNames) ++ " terminals")

    -- this created here and passed to command execution later to remove dependency of renamed data in command execution to
    -- reduce memory footprint keeoing that stuff around.
    let crossReferenceString = CSV.genCsvFile $ CE.getDataListList renamedData dataLeafNames

    -- add in missing terminals to raw data where required
    let reconciledData' = fmap (DT.addMissingTerminalsToInput dataLeafNames []) renamedData

    -- check for data file with all missing data--as in had no terminals with data in termainals list
    let reconciledData = concatMap DT.removeAllMissingCharacters reconciledData'

    let reconciledGraphs = fmap (GFU.reIndexLeavesEdges dataLeafNames . GFU.checkGraphsAndData dataLeafNames) renamedGraphs

    -- Check to see if there are taxa without any observations. Would become total wildcards
    let taxaDataSizeList = filter ((==0).snd) $ zip dataLeafNames $ foldl1 (zipWith (+)) $ fmap (fmap (snd3 . U.filledDataFields (0,0)) . fst) reconciledData
    if not (null taxaDataSizeList) then errorWithoutStackTrace ("\nError: There are taxa without any data: "
            ++ L.intercalate ", " (fmap (Text.unpack . fst) taxaDataSizeList) ++ "\n")
    else hPutStrLn stderr "All taxa contain data"

    -- Ladderizes (resolves) input graphs and ensures that networks are time-consistent
    -- chained netowrk nodes should never be introduced later so only checked no
    -- checks for children of tree node that are all netowork nodee (causes displayu problem)
    let noChainNetNodesList = fmap fromJust $ filter (/=Nothing) $ fmap (LG.removeChainedNetworkNodes True) reconciledGraphs
    let noSisterNetworkNodes = fmap LG.removeTreeEdgeFromTreeNodeWithAllNetworkChildren noChainNetNodesList
    let ladderizedGraphList = fmap (GO.convertGeneralGraphToPhylogeneticGraph "correct") noSisterNetworkNodes

    {-To do
    -- Remove any not "selected" taxa from both data and graphs (easier to remove from fgl)
    let reconciledData' = removeTaxaFromData includeList reconciledData
    let reconciledGraphs' = removeTaxaFromGraphs includeList reconciledData
    -}

    -- Create unique bitvector names for leaf taxa.
    let leafBitVectorNames = DT.createBVNames reconciledData

    {-
    Data processing here-- there are multiple steps not composed so that
    large data files can be precessed and intermediate data goes out
    of scope and can be freed back to system
    -}

    -- Create Naive data -- basic usable format organized into blocks, but not grouped by types, or packed (bit, sankoff, prealigned etc)
    -- Need to check data for equal in character number
    let naiveData = DT.createNaiveData partitionCharOptimalityGlobalSettings reconciledData leafBitVectorNames []

    -- Set reporting data for qualitative characaters to Naive data (usually but not is huge), empty if packed
    let reportingData = if reportNaiveData partitionCharOptimalityGlobalSettings then
                             naiveData
                        else emptyProcessedData

    -- get mix of static/dynamic characters to adjust dynmaicEpsilon
    -- doing on naive data so no packing etc
    let fractionDynamicData = U.getFractionDynamic naiveData

    -- Group Data--all nonadditives to single character, additives with
    -- alphabet < 64 recoded to nonadditive binary, additives with same alphabet
    -- combined,
    let naiveDataGrouped = R.combineDataByType partitionCharOptimalityGlobalSettings naiveData -- R.groupDataByType naiveData

    -- Bit pack non-additive data
    let naiveDataPacked = BP.packNonAdditiveData partitionCharOptimalityGlobalSettings naiveDataGrouped

    -- Optimize Data convert
        -- prealigned to non-additive or matrix
        -- bitPack resulting non-additive
    let optimizedPrealignedData = R.optimizePrealignedData partitionCharOptimalityGlobalSettings naiveDataPacked


    -- Execute any 'Block' change commands--make reBlockedNaiveData
    newBlockPairList <- CE.executeRenameReblockCommands Reblock reBlockPairs thingsToDo

    let reBlockedNaiveData = R.reBlockData newBlockPairList optimizedPrealignedData -- naiveData
    let thingsToDoAfterReblock = filter ((/= Reblock) .fst) $ filter ((/= Rename) .fst) thingsToDoAfterReadRename

    -- Combines data of exact types into single vectors in each block
    -- this is final data processing step
    let optimizedData = if (not . null) newBlockPairList then
                            trace "Reorganizing Block data"
                            R.combineDataByType partitionCharOptimalityGlobalSettings reBlockedNaiveData
                        else optimizedPrealignedData



    -- Set global values before search--should be integrated with executing commands
    -- only stuff that is data dependent here (and seed)
    let defaultGlobalSettings = emptyGlobalSettings { outgroupIndex = 0
                                                    , outGroupName = head dataLeafNames
                                                    , seed = timeD
                                                    , numDataLeaves = length leafBitVectorNames
                                                    , fractionDynamic = fractionDynamicData
                                                    , dynamicEpsilon = 1.0 + ((dynamicEpsilon emptyGlobalSettings - 1.0) * fractionDynamicData)
                                                    }
    -- hPutStrLn stderr ("Fraction characters that are dynamic: " ++ (show $ (fromIntegral lengthDynamicCharacters) / (fromIntegral $ lengthDynamicCharacters + numStaticCharacters)))

    let initialSetCommands = filter ((== Set).fst) thingsToDoAfterReblock
    let commandsAfterInitialDiagnose = filter ((/= Set).fst) thingsToDoAfterReblock

    -- This rather awkward syntax makes sure global settings (outgroup, criterion etc) are in place for initial input graph diagnosis
    (_, initialGlobalSettings, seedList', _) <- CE.executeCommands defaultGlobalSettings (terminalsToExclude, renameFilePairs) numInputFiles crossReferenceString optimizedData optimizedData reportingData [] [] seedList [] initialSetCommands

    -- Get CPUTime so far ()data input and processing
    dataCPUTime <- getCPUTime

    -- Diagnose any input graphs
    let inputGraphList = PU.seqParMap PU.myStrategy  (T.multiTraverseFullyLabelGraph initialGlobalSettings optimizedData True True Nothing) (fmap (LG.rerootTree (outgroupIndex initialGlobalSettings)) ladderizedGraphList)

    -- Get CPUTime for input graphs
    afterGraphDiagnoseTCPUTime <- getCPUTime
    let (inputGraphTime, inGraphNumber, minOutCost, maxOutCost) = if null inputGraphList then (0, 0, infinity, infinity)
                                                                else (fromIntegral afterGraphDiagnoseTCPUTime - fromIntegral dataCPUTime, length inputGraphList, minimum $ fmap snd6 inputGraphList, maximum $ fmap snd6 inputGraphList)

    let inputProcessingData   = emptySearchData {commentString = "Input and data processing", duration = fromIntegral dataCPUTime}
    let inputGraphProcessing  = emptySearchData {minGraphCostOut = minOutCost, maxGraphCostOut = maxOutCost, numGraphsOut = inGraphNumber, commentString = "Input graph processing", duration = inputGraphTime}


    -- Create lazy pairwise distances if needed later for build or report
    let pairDist = D.getPairwiseDistances optimizedData

    -- Execute Following Commands (searches, reports etc)
    (finalGraphList, _, _, _) <- CE.executeCommands (initialGlobalSettings {searchData = [inputGraphProcessing, inputProcessingData]}) (terminalsToExclude, renameFilePairs) numInputFiles crossReferenceString optimizedData optimizedData reportingData inputGraphList pairDist seedList' [] commandsAfterInitialDiagnose -- (transformString ++ commandsAfterInitialDiagnose)

    -- print global setting just to check
    --hPutStrLn stderr (show _finalGlobalSettings)

    -- Add in model and root cost if optimality criterion needs it
    -- if (rootComplexity initialGlobalSettings) /= 0.0 then hPutStrLn stderr ("\tUpdating final graph with any root priors")
    -- else hPutStrLn stderr ""

    -- rediagnose for NCM due to packing, in most cases not required, just being sure etc
    let rediagnoseWithReportingdata = True
    let finalGraphList' = T.updateGraphCostsComplexities initialGlobalSettings reportingData optimizedData rediagnoseWithReportingdata finalGraphList

    let minCost = if null finalGraphList then 0.0 else minimum $ fmap snd6 finalGraphList'
    let maxCost = if null finalGraphList then 0.0 else maximum $ fmap snd6 finalGraphList'


    -- final results reporting to stderr
    hPutStrLn stderr ("Execution returned " ++ show (length finalGraphList') ++ " graph(s) at cost range " ++ show (minCost, maxCost))

    -- Final Stderr report
    timeCPUEnd <- getCPUTime
    timeCDEnd <- getCurrentTime

    --hPutStrLn stderr ("CPU Time " ++ (show timeCPUEnd))
    let wallClockDuration = floor (1000000000000 * nominalDiffTimeToSeconds (diffUTCTime timeCDEnd timeCD)) :: Integer
    --hPutStrLn stderr ("Current time " ++  (show wallClockDuration))
    let cpuUsage = fromIntegral timeCPUEnd / fromIntegral wallClockDuration :: Double
    --hPutStrLn stderr ("CPU % " ++ (show cpuUsage))

    hPutStrLn stderr ("\n\tWall-Clock time " ++ show ((fromIntegral wallClockDuration :: Double) / 1000000000000.0) ++ " second(s)"
        ++ "\n\tCPU time " ++ show ((fromIntegral timeCPUEnd :: Double) / 1000000000000.0) ++ " second(s)"
        ++ "\n\tCPU usage " ++ show (floor (100.0 * cpuUsage) :: Integer) ++ "%"
        )

