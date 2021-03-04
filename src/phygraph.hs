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

--import           Debug.Trace


-- | main driver
main :: IO ()
main =
  do
    let splash = "\nPhyGraph version 0.1\nCopyright(C) 2021 Ward Wheeler and The American Museum of Natural History\n"
    let splash2 = "PhyGraph comes with ABSOLUTELY NO WARRANTY; This is free software, and may be \nredistributed "
    let splash3 = "under the 3-Clause BSD License.\n"
    hPutStrLn stderr (splash ++ splash2 ++ splash3)
    
    -- Process arguments--a single file containing commands
    args <- getArgs

    if length args /= 1 then errorWithoutStackTrace "\nProgram requires a single argument--the name of command script file.\n\n"
    else hPutStr stderr "\nCommand script file: "
    hPutStrLn stderr $ head args

    -- Process commands to get list of actions
    commandContents <- readFile $ head args
    let thingsToDo = commandList commandContents

    -- Prcess Read, rename commands

    -- Reconcile Data

    -- Optimize Data

    -- Searchs Actions (wit intermeduate reports)

    -- Final Report actions

    mapM_ (hPutStrLn stderr) (fmap show thingsToDo)

    hPutStrLn stderr "\nDone"


