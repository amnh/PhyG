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

module Main (main) where

import qualified Commands.CommandExecution    as CE
import qualified Commands.ProcessCommands     as PC
import qualified Data.List                    as L
import           GeneralUtilities
import qualified GraphFormatUtilities         as GFU
import qualified GraphOptimization.Traversals as T
import qualified Graphs.GraphOperations       as GO
import qualified Input.DataTransformation     as DT
import qualified Input.ReadInputFiles         as RIF
import           System.Environment
import           System.IO
import           Types.Types
import qualified Utilities.Distances          as D
import qualified Data.Text.Lazy               as Text
import qualified Data.Text.Short              as ST
import qualified Utilities.Utilities          as U
import qualified Input.Reorganize             as R
-- import qualified Utilities.LocalGraph         as LG

-- | debugRandomFlag sets flag for debugging randomization
debugRandomFlag :: Bool
debugRandomFlag = False

-- | main driver
main :: IO ()
main = do
    let splash = "\nPhyG version " ++ pgVersion ++ "\nCopyright(C) 2021 Ward Wheeler and The American Museum of Natural History\n"
    let splash2 = "PhyG comes with ABSOLUTELY NO WARRANTY; This is free software, and may be \nredistributed "
    let splash3 = "under the 3-Clause BSD License.\n"
    hPutStrLn stderr (splash ++ splash2 ++ splash3)

    -- Process arguments--a single file containing commands
    args <- getArgs

    if length args /= 1 then errorWithoutStackTrace "\nProgram requires a single argument--the name of command script file.\n\n"
    else hPutStr stderr "\nCommand script file: "
    hPutStrLn stderr $ head args

    -- System time for Random seed
    timeD <- getSystemTimeSeconds
    -- hPutStrLn stderr ("Current time is " ++ show timeD)
    let seedList = if not debugRandomFlag then randomIntList timeD
                   else randomIntList 0


    -- Process commands to get list of actions
    commandContents <- readFile $ head args

      -- Process so one command per line?

    commandContents' <- PC.expandRunCommands [] (lines commandContents)
    let thingsToDo' = PC.getCommandList  commandContents'
    --mapM_ (hPutStrLn stderr) (fmap show thingsToDo')

    -- Process Read commands (with prealigned and tcm flags first)
        --expand read commands for wildcards
    expandedReadCommands <- mapM (RIF.expandReadCommands []) $ filter ((== Read) . fst) thingsToDo'
    let thingsToDo = (concat expandedReadCommands) ++ (filter ((/= Read) . fst) thingsToDo')
    --hPutStrLn stderr (show $ concat expandedReadCommands)

    dataGraphList <- mapM RIF.executeReadCommands $ fmap PC.movePrealignedTCM $ fmap snd $ filter ((== Read) . fst) thingsToDo
    let (rawData, rawGraphs, terminalsToInclude, terminalsToExclude, renameFilePairs, reBlockPairs) = RIF.extractInputTuple dataGraphList

    if null rawData && null rawGraphs then errorWithoutStackTrace "\n\nNeither data nor graphs entered.  Nothing can be done."
    else hPutStrLn stderr ("Entered " ++ (show $ length rawData) ++ " data file(s) and " ++ (show $ length rawGraphs) ++ " input graphs")

    -- Split fasta/fastc sequences into corresponding pieces based on '#' partition character
    let rawDataSplit = DT.partitionSequences (ST.fromString "#") rawData

    -- Process Rename Commands
    newNamePairList <- CE.executeRenameReblockCommands renameFilePairs thingsToDo
    if (not $ null newNamePairList) then hPutStrLn stderr ("Renaming " ++ (show $ length newNamePairList) ++ " terminals") -- ++ (show $ L.sortBy (\(a,_) (b,_) -> compare a b) newNamePairList))
    else hPutStrLn stderr ("No terminals to be renamed")

    let renamedData   = fmap (DT.renameData newNamePairList) rawDataSplit
    let renamedGraphs = fmap (GFU.relabelGraphLeaves  newNamePairList) rawGraphs

    let thingsToDoAfterReadRename = (filter ((/= Read) .fst) $ filter ((/= Rename) .fst) thingsToDo)

    -- Reconcile Data and Graphs (if input) including ladderization
        -- could be sorted, but no real need
        -- get taxa to include in analysis
    if not $ null terminalsToInclude then hPutStrLn stderr ("Terminals to include:" ++show terminalsToInclude)
    else hPutStrLn stderr ("")   
    if not $ null terminalsToExclude then hPutStrLn stderr ("Terminals to exclude:" ++show terminalsToExclude)
    else hPutStrLn stderr ("")

    -- Uses names form terminal list if non-null, and remove exckuded terminals
    let dataLeafNames' = if (not $ null terminalsToInclude) then L.sort $ L.nub terminalsToInclude
                        else L.sort $ DT.getDataTerminalNames renamedData
    let dataLeafNames = dataLeafNames' L.\\ terminalsToExclude
    hPutStrLn stderr ("Data were input for " ++ (show $ length dataLeafNames) ++ " terminals")


    let reconciledData = fmap (DT.addMissingTerminalsToInput dataLeafNames []) renamedData
    let reconciledGraphs = fmap (GFU.reIndexLeavesEdges dataLeafNames) $ fmap (GFU.checkGraphsAndData dataLeafNames) renamedGraphs
    
    -- Check to see if there are taxa without any observations. Would become total wildcards
    let taxaDataSizeList = filter ((==0).snd) $ zip dataLeafNames $ foldl1 (zipWith (+)) $ fmap (fmap snd3) $ fmap (fmap (U.filledDataFields (0,0))) $ fmap fst reconciledData
    if (length taxaDataSizeList /= 0) then hPutStrLn stderr ("\nWarning (but a serious one): There are taxa without any data: " 
            ++ (L.intercalate ", " $ fmap Text.unpack $ fmap fst taxaDataSizeList) ++ "\n")
    else hPutStrLn stderr "All taxa contain data"

    -- Ladderizes (resolves) input graphs and verifies that networks are time-consistent
    let ladderizedGraphList = fmap GO.verifyTimeConsistency $ fmap GO.ladderizeGraph reconciledGraphs

    {-To do
    -- Remove any not "selected" taxa from both data and graphs (easier to remove from fgl)
    let reconciledData' = removeTaxaFromData includeList reconciledData
    let reconciledGraphs' = removeTaxaFromGraphs includeList reconciledData
    -}

    -- Create unique bitvector names for leaf taxa.
    let leafBitVectorNames = DT.createBVNames reconciledData

    -- Create Naive data -- basic usable format organized into blocks, but not grouped by types, or packed (bit, sankoff, prealigned etc)
    -- Need to check data for equal in character number
    let naiveData = DT.createNaiveData reconciledData leafBitVectorNames []

    
    -- Execute any 'Block' change commands--make reBlockedNaiveData
    newBlockPairList <- CE.executeRenameReblockCommands reBlockPairs thingsToDo
    let reBlockedNaiveData = R.reBlockData newBlockPairList naiveData
    let thingsToDoAfterReblock = filter ((/= Reblock) .fst) $ filter ((/= Rename) .fst) thingsToDoAfterReadRename


    -- Group Data--all nonadditives to single character, additives with same alphabet, 
    let groupedData = R.groupDataByType reBlockedNaiveData

    -- Optimize Data convert
        -- Additive characters with alphabets < 64 to multiple binary nonadditive
        -- all binary characters to nonadditive
        -- matrix 2 states to non-additive with weight
        -- prealigned to non-additive or matrix
        -- bitPack non-additive
    let optimizedData = groupedData --  place holder (consolidate all add, tcms, non-add etc chars in blocks)


    -- Set global vaues before search--should be integrated with executing commands
    let defaultGlobalSettings = GlobalSettings { outgroupIndex = 0
                                               , outGroupName = head dataLeafNames
                                               , optimalityCriterion = Parsimony
                                               , graphType = Tree
                                               , compressResolutions = True
                                               , finalAssignment = ImpliedAlignment
                                               , graphFactor = Wheeler2015Network
                                               , rootCost = NoRootCost
                                               }
    --hPutStrLn stderr (show defaultGlobalSettings)

    let initialSetCommands = filter ((== Set).fst) thingsToDoAfterReblock
    let commandsAfterInitialDiagnose = filter ((/= Set).fst) thingsToDoAfterReblock

    -- This rather awkward syntax makes sure global settings (outgroup, criterion etc) are in place for initial input graph diagnosis
    (_, initialGlobalSettings) <- CE.executeCommands defaultGlobalSettings renamedData optimizedData [] [] seedList initialSetCommands
    let inputGraphList = map (T.multiTraverseFullyLabelGraph initialGlobalSettings optimizedData True True Nothing) (fmap (GO.rerootTree (outgroupIndex initialGlobalSettings)) ladderizedGraphList)
    

    -- Create lazy pairwise distances if needed later for build or report
    let pairDist = D.getPairwiseDistances optimizedData
    

    -- Execute Following Commands (searches, reports etc)
    (finalGraphList, finalGlobalSettings) <- CE.executeCommands initialGlobalSettings renamedData optimizedData inputGraphList pairDist seedList commandsAfterInitialDiagnose

    -- print global setting just to check
    --hPutStrLn stderr (show _finalGlobalSettings)

    -- Final Stderr report
    timeDN <- getSystemTimeSeconds
    let minCost = if null finalGraphList then 0.0 else minimum $ fmap snd6 finalGraphList
    let maxCost = if null finalGraphList then 0.0 else maximum $ fmap snd6 finalGraphList
    hPutStrLn stderr ("Execution returned " ++ (show $ length finalGraphList) ++ " graph(s) at cost range " ++ (show (minCost, maxCost)) ++ " in "++ show (timeDN - timeD) ++ " second(s)")

