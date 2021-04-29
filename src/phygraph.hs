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

import           System.IO
import           System.Environment

import qualified ProcessCommands as PC
import           Types
import qualified ReadInputFiles as RIF
import           GeneralUtilities
import qualified CommandExecution as CE
import qualified GraphFormatUtilities as GFU
import qualified DataTransformation as DT
import qualified Distances as D
import           Data.List
--import           Debug.Trace


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

    hPutStrLn stderr ("Entered " ++ (show $ length rawData) ++ " data file(s) and " ++ (show $ length rawGraphs) ++ " input graphs")

    -- Process Rename Commands
    newNamePairList <- CE.executeRenameCommands [] thingsToDo
    if (not $ null newNamePairList) then hPutStrLn stderr ("Renaming " ++ (show $ length newNamePairList) ++ " terminals")
    else hPutStrLn stderr ("No terminals to be renamed")

    let renamedData = fmap (DT.renameData newNamePairList) rawData
    let renamedGraphs =  fmap (GFU.relabelGraphLeaves  newNamePairList) rawGraphs

    -- Reconcile Data and Graphs (if input)
    let dataLeafNames = sort $ DT.getDataTerminalNames renamedData
    hPutStrLn stderr ("Data were input for " ++ (show $ length dataLeafNames) ++ " terminals")

    let reconciledData = fmap (DT.addMissingTerminalsToInput dataLeafNames) renamedData 
    let reconciledGraphs = fmap (GFU.checkGraphsAndData dataLeafNames) renamedGraphs -- fmap GFU.toPhylogeneticGraph $ 

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

    -- To test data recoding and basic median2
    let pairDist = D.getPairwiseDistances naiveData
    --hPutStrLn stderr (show pairDist)


    -- Optimize Data
    let optimizedData = naiveData --  place holder

    -- Execute Following Commands (searches, reports etc)
    finalGraphList <- CE.executeCommands rawData optimizedData reconciledGraphs pairDist (filter ((/= Read) .fst) $ filter ((/= Rename) .fst) thingsToDo)

    timeDN <- getSystemTimeSeconds 
    hPutStrLn stderr ("Execution time " ++ show (timeDN - timeD))
    
    -- Final Stderr report
    hPutStrLn stderr ("Execution returned " ++ (show $ length finalGraphList) ++ " graphs")
    hPutStrLn stderr "\nDone"


{-
    hPutStrLn stderr ("\tData for " ++ (show $ fmap length $ fst $ head rawData))
    hPutStrLn stderr ("\tAlp[habet] for " ++ (show  $ fmap snd rawData))

-}