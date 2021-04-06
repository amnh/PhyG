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

import           ProcessCommands
import           Types
import qualified ReadInputFiles as RIF
import           GeneralUtilities
import qualified CommandExecution as CE
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
    let thingsToDo = getCommandList  commandContents
    mapM_ (hPutStrLn stderr) (fmap show thingsToDo)

    -- Process Read, rename commands
    (rawData, rawGraphs) <- RIF.executeReadCommands [] [] $ concat $ fmap snd $ filter ((== Read) . fst) thingsToDo
    hPutStrLn stderr ("Entered " ++ (show $ length rawData) ++ " data file(s) and " ++ (show $ length rawGraphs) ++ " input graphs")

    -- Reconcile Data and Graphs (if input)
    let (reconciledData, reconciledGraphs) = (rawData, rawGraphs) -- place holder

    -- Optimize Data
    let optimizedData = reconciledData --  place holder

    -- Execute Following Commands (searches, reports etc)
    finalGraphList <- CE.executeCommands optimizedData reconciledGraphs thingsToDo

    
    -- Final Stderr report
    hPutStrLn stderr ("Search returned " ++ (show $ length finalGraphList) ++ " phylogenetic graphs")
    hPutStrLn stderr "\nDone"


{-
    hPutStrLn stderr ("\tData for " ++ (show $ fmap length $ fst $ head rawData))
    hPutStrLn stderr ("\tAlp[habet] for " ++ (show  $ fmap snd rawData))

-}