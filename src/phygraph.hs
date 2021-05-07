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

module Main where

--import           Debug.Trace
import           System.IO
import           System.Environment
import           Data.List
import           Types
import qualified ProcessCommands as PC
import qualified ReadInputFiles as RIF
import           GeneralUtilities
import qualified CommandExecution as CE
import qualified GraphFormatUtilities as GFU
import qualified DataTransformation as DT
import qualified Distances as D
import qualified GraphOperations as GO


-- | main driver
main :: IO ()
main =
  do
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
    hPutStrLn stderr ("Current time is " ++ show timeD)
    

    -- Process commands to get list of actions
    commandContents <- readFile $ head args

      -- Process so one command per line?

    commandContents' <- PC.expandRunCommands [] (lines commandContents) 
    let thingsToDo = PC.getCommandList  commandContents'
    --mapM_ (hPutStrLn stderr) (fmap show thingsToDo)

    -- Process Read commands (with prealigned and tcm flags first)
    dataGraphList <- mapM (RIF.executeReadCommands [] [] False ([],[])) $ fmap PC.movePrealignedTCM $ fmap snd $ filter ((== Read) . fst) thingsToDo
    let (rawData, rawGraphs) = RIF.extractDataGraphPair dataGraphList

    if null rawData && null rawGraphs then errorWithoutStackTrace "\n\nNeither data nor graphs entered.  Nothing can be done."
    else hPutStrLn stderr ("Entered " ++ (show $ length rawData) ++ " data file(s) and " ++ (show $ length rawGraphs) ++ " input graphs")

    -- Process Rename Commands
    newNamePairList <- CE.executeRenameCommands [] thingsToDo
    if (not $ null newNamePairList) then hPutStrLn stderr ("Renaming " ++ (show $ length newNamePairList) ++ " terminals")
    else hPutStrLn stderr ("No terminals to be renamed")

    let renamedData = fmap (DT.renameData newNamePairList) rawData
    let renamedGraphs =  fmap (GFU.relabelGraphLeaves  newNamePairList) rawGraphs

    let thingsToDoAfterReadRename = (filter ((/= Read) .fst) $ filter ((/= Rename) .fst) thingsToDo)

    -- Reconcile Data and Graphs (if input) including ladderization
        -- could be sorted, but no real need
    let dataLeafNames = DT.getDataTerminalNames renamedData
    hPutStrLn stderr ("Data were input for " ++ (show $ length dataLeafNames) ++ " terminals")
    --hPutStrLn stderr (show dataLeafNames)

    
    let reconciledData = fmap (DT.addMissingTerminalsToInput dataLeafNames) renamedData 
    let reconciledGraphs = fmap (GFU.reIndexLeavesEdges dataLeafNames) $ fmap (GFU.checkGraphsAndData dataLeafNames) renamedGraphs  

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
    -- Need to check data for equal in charcater number
    let naiveData = DT.createNaiveData reconciledData leafBitVectorNames []
  
    {-To do
      Execute any 'Block' change commands--make reBlockedNaiveData
    -}

    -- Optimize Data
    let optimizedData = naiveData --  place holder (consolidate all add, non-add etc chars in blocks)


    -- Set global vaues before search--should be integrated with executing commands
    let defaultGlobalSettings = GlobalSettings {outgroupIndex = 0, outGroupName = head dataLeafNames, optimalityCriterion = Parsimony, graphType = Tree}
    --hPutStrLn stderr (show defaultGlobalSettings)

    let initialSetCommands = takeWhile ((== Set).fst) thingsToDoAfterReadRename 
    let commandsAfterInitialDiagnose = dropWhile ((== Set).fst) thingsToDoAfterReadRename 
    
    -- This rather awkward syntax makes sure global settings (outgroup, criterion etc) are in place for initial input graph diagnosis
    (_, initialGlobalSettings) <- CE.executeCommands defaultGlobalSettings rawData optimizedData [] [] initialSetCommands
    let inputGraphList = map (GO.fullyLabelGraph initialGlobalSettings optimizedData) (fmap (GO.rerootGraph (outgroupIndex initialGlobalSettings)) ladderizedGraphList)

    -- Create lazy pairwise distances if needed later for build or report
    let pairDist = D.getPairwiseDistances optimizedData
    --hPutStrLn stderr (show pairDist)

    
    -- Execute Following Commands (searches, reports etc)
    (finalGraphList, finalGlobalSettings) <- CE.executeCommands initialGlobalSettings rawData optimizedData inputGraphList pairDist commandsAfterInitialDiagnose

    -- print global setting just to check
    --hPutStrLn stderr (show finalGlobalSettings)

    -- Final Stderr report
    timeDN <- getSystemTimeSeconds 
    hPutStrLn stderr ("Execution returned " ++ (show $ length finalGraphList) ++ " graph(s) in "++ show (timeDN - timeD) ++ " second(s)")
    

{-
    hPutStrLn stderr ("\tData for " ++ (show $ fmap length $ fst $ head rawData))
    hPutStrLn stderr ("\tAlp[habet] for " ++ (show  $ fmap snd rawData))

-}