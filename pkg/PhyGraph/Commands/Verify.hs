{- |
Module      :  Verify.hs
Description :  Module to verify (more or less) input commands 
Copyright   :  (c) 2022 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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



module Commands.Verify
    ( verifyCommands
    , allowedCommandList
    , buildArgList
    , fuseArgList
    , geneticAlgorithmArgList
    , netEdgeArgList
    , readArgList
    , reconcileArgList
    , refineArgList
    , reportArgList
    , searchArgList
    , setArgList
    , selectArgList
    , supportArgList
    , swapArgList
    , transformArgList
    ) where

import           Types.Types
import           GeneralUtilities
import qualified Data.List as L
import qualified Data.Char as C
import           Text.Read

-- import           Debug.Trace

-- | allowedCommandList is the permitted command string list
allowedCommandList :: [String]
allowedCommandList = ["build", "fuse", "read", "reblock", "refine", "rename", "report", "run", "search", "select", "set", "support", "swap"]

-- list of valid instructions
validInstructionList :: [Instruction]
validInstructionList =  [Build, Fuse, Read, Reblock, Refine, Rename, Report, Run, Select, Set, Swap, Search, Support, Transform]

-- | buildArgList is the list of valid build arguments
buildArgList :: [String]
buildArgList = ["replicates", "nj", "wpgma", "dwag", "rdwag", "distance", "character", "best","none","otu","spr","tbr", "block","cun", "eun", "atrandom", "first", "displaytrees", "graph"]

-- | fuseArgList arguments
fuseArgList :: [String]
fuseArgList = ["spr","tbr", "keep", "steepest", "all", "nni", "best", "unique", "once", "atrandom", "pairs"]

-- | geneticAlgorithm arguments
geneticAlgorithmArgList :: [String]
geneticAlgorithmArgList = ["popsize", "generations", "elitist", "severity", "recombinations","geneticalgorithm", "ga", "maxnetedges"]

-- | netEdgeArgList arguments for network edge add/delete operations
netEdgeArgList :: [String]
netEdgeArgList = ["keep", "steepest", "all", "netadd", "netdel", "netdelete", "netadddel", "netadddelete", "netmove", "annealing", "steps", "returnmutated", "drift", "acceptequal", "acceptworse", "maxchanges","steepest","atrandom", "maxnetedges", "rounds"]

-- | Read arg list allowable modifiers in read
readArgList :: [String]
readArgList = ["tcm", "nucleotide", "aminoacid", "fasta", "fastc", "tnt", "csv",
    "dot", "newick" , "enewick", "fenewick", "include", "exclude", "rename", "block", "prefasta", 
    "prefastc", "preaminoacid", "prenucleotide"] -- "prealigned", 

-- should be moved to a single file for import
-- | reconcileCommandList list of allowable commands
reconcileArgList :: [String]
reconcileArgList = ["method", "compare", "threshold", "outformat", "connect", "edgelabel", "vertexlabel"] -- "outfile" 

-- | reconcileOptionsList list of allowable command options of method, compare, threshhold, and outformat
reconcileOptionsList :: [String]
reconcileOptionsList = ["eun", "cun", "strict", "majority", "adams", "dot" ,"dotpdf", "fen", "newick", "true", "false", "combinable", "identity"]


-- | refinement arguments
refineArgList :: [String]
refineArgList = fuseArgList ++ netEdgeArgList ++ geneticAlgorithmArgList

-- | reportArgList contains valid 'report' arguments
reportArgList :: [String]
reportArgList = ["all", "data", "search", "graphs", "overwrite", "append", "dot", "dotpdf", "newick", "ascii", "crossrefs", "pairdist", "diagnosis","displaytrees", "reconcile", "support", "ia", "impliedalignment", "tnt", "includemissing", "concatenate", "htulabels", "branchlengths", "nohtulabels", "nobranchlengths", "collapse", "nocollapse"] ++ reconcileArgList

-- | search arguments
searchArgList :: [String]
searchArgList = ["days", "hours", "minutes", "seconds", "instances", "thompson", "mfactor", "linear", "exponential", "maxnetadges"]

-- | buildArgList is the list of valid build arguments
selectArgList :: [String]
selectArgList = ["best", "all", "unique", "atrandom"]

-- | setArgLIst contains valid 'set' arguments
setArgList :: [String]
setArgList = ["outgroup", "criterion", "graphtype", "compressresolutions", "finalassignment", "graphfactor", "rootcost", "seed","partitioncharacter", "modelcomplexity", 
    "bc2", "bc4", "bc5", "bc8", "bc64", "bcgt64", "dynamicepsilon"]

-- | refinement arguments
supportArgList :: [String]
supportArgList = ["bootstrap", "jackknife", "goodmanbremer", "gb", "gbsample", "replicates", "buildonly", "atrandom"] -- "bootstrap", 

-- | buildArgList is the list of valid build arguments
swapArgList :: [String]
swapArgList = ["spr","tbr", "alternate", "keep", "steepest", "all", "nni", "ia", "annealing", "maxtemp", "mintemp", "steps", "returnmutated", "drift", "acceptequal", "acceptworse", "maxchanges"]

-- | transform arguments
transformArgList :: [String]
transformArgList = ["totree", "tosoftwired", "tohardwired", "staticapprox", "dynamic", "atrandom", "first", "displaytrees", "weight", "name", "type", "dynamicepsilon", "outgroup"]


-- | verifyCommands takes a command list and tests whether the commands 
-- and arguments are permissible before program execution--prevents late failure 
-- after alot of processing time.
-- bit does not check for files existance, write/read ability, or contents for format or 
-- anyhhting else for that matter
-- does check if files are both read from and written to
verifyCommands :: [Command] -> [String] -> [String] -> Bool
verifyCommands inCommandList inFilesToRead inFilesToWrite = 
    if null inCommandList then True
    else 
        let firstCommand = head inCommandList
            commandInstruction = fst firstCommand
            inArgs = snd firstCommand

            -- check valid commandInstructions
            -- this is done earlier but might get oved so putting here just in case
            checkInstruction = commandInstruction `elem` validInstructionList
     
        in
        if not checkInstruction then errorWithoutStackTrace ("Invalid command was specified : " ++ (show commandInstruction))
        else 
            -- check each command for valid arguments
            -- make lower-case arguments
            let fstArgList = filter (/= []) $ fmap (fmap C.toLower . fst) inArgs
                sndArgList = filter (/= []) $ fmap (fmap C.toLower . snd) inArgs
                fileNameList = fmap (filter (/= '"')) $ filter (/= []) $ fmap snd inArgs

                -- Read
                (checkOptions, filesToReadFrom, filesToWriteTo) = 

                                   -- Build
                                   if commandInstruction == Build then 
                                        (checkCommandArgs "build" fstArgList buildArgList, [""], [""])

                                   -- Fuse 
                                   else if commandInstruction == Fuse then 
                                        (checkCommandArgs "fuse" fstArgList fuseArgList, [""], [""])

                                   else if commandInstruction == Read then 
                                        let fileArgs = concat $ filter (/= []) $ fmap snd inArgs
                                            numDoubleQuotes = length $ filter (== '"') fileArgs
                                            has02DoubleQuotes = (numDoubleQuotes == 2) || (numDoubleQuotes == 0)
                                        in
                                        if not has02DoubleQuotes then errorWithoutStackTrace ("Unbalanced quotation marks in 'read' file argument: " ++ fileArgs)
                                        else (checkCommandArgs "read"  fstArgList readArgList, fileNameList, [""])

                                   -- Reblock -- no arguments--but reads string and blocknames
                                   else if commandInstruction == Reblock then 
                                        let fileArgs = concat $ filter (/= []) $ fmap snd inArgs
                                            numDoubleQuotes = length $ filter (== '"') fileArgs
                                            (numDouble, numUnbalanced) = divMod numDoubleQuotes 2
                                        in
                                        -- trace (show (fstArgList, sndArgList,fileNameList)) (
                                        if numDouble < 2 then errorWithoutStackTrace ("Need at least two fields in 'rebock' command, new block name and old block(s) in double quotes: " ++ fileArgs)
                                        else if numUnbalanced /= 0 then errorWithoutStackTrace ("Unbalanced quotation marks in 'reblock' command: " ++ fileArgs)
                                        else (True,[""], [""])
                                        -- )

                                   -- Reconcile -- part of report

                                   -- Refine
                                   else if commandInstruction == Refine then 
                                        (checkCommandArgs "refine" fstArgList refineArgList,[""], [""])

                                   -- Rename -- -- no arguments--but reads string and taxon names
                                   else if commandInstruction == Rename then 
                                        let fileArgs = concat $ filter (/= []) $ fmap snd inArgs
                                            numDoubleQuotes = length $ filter (== '"') fileArgs
                                            (numDouble, numUnbalanced) = divMod numDoubleQuotes 2
                                        in
                                        -- trace (show (fstArgList, sndArgList,fileNameList)) (
                                        if numDouble < 2 then errorWithoutStackTrace ("Need at least two fields in 'rename' command, new taxon name and old taxon name(s) in double quotes: " ++ fileArgs)
                                        else if numUnbalanced /= 0 then errorWithoutStackTrace ("Unbalanced quotation marks in 'rename' command: " ++ fileArgs)
                                        else (True,[""], [""])
                                        -- )

                                   -- Report
                                   else if commandInstruction == Report then 
                                        let fileArgs = concat $ filter (/= []) $ fmap snd inArgs
                                            numDoubleQuotes = length $ filter (== '"') fileArgs
                                            has02DoubleQuotes = (numDoubleQuotes == 2) || (numDoubleQuotes == 0)
                                        in
                                        if not has02DoubleQuotes then errorWithoutStackTrace ("Unbalanced quotation marks in report file argument: " ++ fileArgs)
                                        else (checkCommandArgs "report" fstArgList reportArgList, [""], fileNameList)

                                   -- Run  -- processed out before this into command list

                                   -- Search
                                   else if commandInstruction == Search then 
                                        if "reconcile" `notElem` fstArgList then (checkCommandArgs "report" fstArgList searchArgList, [""], [""])
                                        else 
                                            let reconcilePairList = zip fstArgList sndArgList
                                                nonThresholdreconcileModPairList = filter ((/= "threshold"). fst) $ reconcilePairList
                                                thresholdreconcileModPairList = filter ((== "threshold"). fst) $ reconcilePairList
                                                checkReconcile1 = checkCommandArgs "reconcile"  (fmap fst nonThresholdreconcileModPairList) reconcileArgList
                                                checkReconcile2 = checkCommandArgs "reconcile"  (fmap fst thresholdreconcileModPairList) reconcileArgList
                                                checkReconcile3 = checkCommandArgs "reconcile modifier (method, compare, outformat, connect, edgelabel, vertexlabel)"  
                                                    (fmap snd nonThresholdreconcileModPairList) reconcileOptionsList
                                                checkReconcile4 = L.foldl1' (&&) $ True : (fmap isInt (filter (/= []) (fmap snd thresholdreconcileModPairList)))
                                                checkReconcile = checkReconcile1 && checkReconcile2 && checkReconcile3 && checkReconcile4
                                            in
                                            if checkReconcile then (checkCommandArgs "report" fstArgList searchArgList,[""], [""])
                                            else (False, [], [])
                                   -- Select
                                   else if commandInstruction == Select then 
                                        (checkCommandArgs "select" fstArgList selectArgList, [""], [""])

                                   -- Set
                                   else if commandInstruction == Set then 
                                        (checkCommandArgs "set" fstArgList setArgList, [""], [""])

                                   -- Support
                                   else if commandInstruction == Support then 
                                        (checkCommandArgs "support" fstArgList supportArgList, [""], [""])

                                   -- Swap
                                   else if commandInstruction == Swap then 
                                        (checkCommandArgs "swap" fstArgList swapArgList, [""], [""])

                                   -- Transform 
                                   else if commandInstruction == Transform then 
                                        (checkCommandArgs "transform" fstArgList transformArgList, [""], [""])
                                         
                                   else errorWithoutStackTrace ("Unrecognized command was specified : " ++ (show commandInstruction))

                                   

            in
            if checkOptions then
                let allFilesToReadFrom = filter (/= "") $ filesToReadFrom ++ inFilesToRead
                    allFilesToWriteTo = filter (/= "") $ filesToWriteTo ++ inFilesToWrite
                    readAndWriteFileList = L.intersect allFilesToReadFrom allFilesToWriteTo
                in
                -- trace (show (allFilesToReadFrom, allFilesToWriteTo)) (
                if (not .null) readAndWriteFileList then 
                    errorWithoutStackTrace ("Error--Both reading from and writing to files (could cause errors and/or loss of data): " ++ (show readAndWriteFileList))
                else verifyCommands (tail inCommandList) allFilesToReadFrom allFilesToWriteTo
                -- )

            else 
                -- Won't get to here--will error at earlier stages
                False
            
        where isInt a = if (readMaybe a :: Maybe Int) /= Nothing then True else False

