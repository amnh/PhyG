{- |
Module      :  ReadInputFiles.hs
Description :  Module to read input files for phylogenetic analysis
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



module ReadInputFiles  where

import           Types
import           Debug.Trace
import           Data.Char


-- | Read arg list allowable modifiers in read
readArgList :: [String]
readArgList = ["tcm", "prealigned", "fasta", "fastc", "tnt", "csv", "dot", "newick" , "enewick", "fenewick", "rename"]

-- | getReadArgs processes arguments ofr the 'read' command
-- should allow mulitple files and gracefully error check
-- also allows tcm file specification (limit 1 tcm per command?)
-- as fasta, fastc, tnt, tcm, prealigned
getReadArgs :: String -> [(String, String)] -> [Argument]
getReadArgs fullCommand argList = 
    if null argList then []
    else 
        let (firstPart, secondPart) = head argList 
        in
        -- plain file name with no modifier
        if (null firstPart) then
            if (head secondPart == '"') || (last secondPart == '"') then (firstPart, (init $ tail secondPart)) : getReadArgs fullCommand (tail argList)
            else errorWithoutStackTrace ("\n\n'Read' command error '" ++ (secondPart) ++"' : Need to specify filename in double quotes") 
        -- Change to allowed modifiers
        else if ((fmap toLower firstPart) `notElem` readArgList)  then errorWithoutStackTrace ("\n\n'Read' command error: " ++ fullCommand ++ " contains unrecognized option '" ++ firstPart ++ "'")
        else if (null secondPart) || (head secondPart /= '"') || (last secondPart /= '"') then errorWithoutStackTrace ("\n\n'Read' command error '" ++ (secondPart) ++"' : Need to specify filename in double quotes") 
        else (firstPart, (init $ tail secondPart)) : getReadArgs fullCommand (tail argList)


-- | executeReadCommands
executeReadCommands :: [RawData] -> [RawGraph] -> [Argument] -> IO ([RawData], [RawGraph])
executeReadCommands curData curGraphs argList = 
    if null argList then return (curData, curGraphs)
    else 
        let (firstOption, firstFile) = head argList
        in 
        -- try to figure out file type
        if null firstOption then executeReadCommands curData curGraphs (tail argList)
        else if firstOption == "fasta" then executeReadCommands curData curGraphs (tail argList)
        else if firstOption == "fastc" then executeReadCommands curData curGraphs (tail argList)
        else if firstOption == "tnt" then executeReadCommands curData curGraphs (tail argList)
        else if firstOption == "dot" then executeReadCommands curData curGraphs (tail argList)
        else if firstOption == "tcm" then executeReadCommands curData curGraphs (tail argList)
        else if firstOption == "prealigned" then executeReadCommands curData curGraphs (tail argList)
        else if firstOption == "rename" then executeReadCommands curData curGraphs (tail argList)
        else if (reverse $ take 3 (reverse firstOption)) == "ick" then executeReadCommands curData curGraphs (tail argList)
        else errorWithoutStackTrace ("\n\n'Read' command error: option " ++ firstOption ++ " not recognized/implemented")
        