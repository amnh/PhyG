{- |
Module      :  Verify.hs
Description :  Module to verify (more or less) input commands
Copyright   :  (c) 2022-2023 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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

import Data.Char qualified as C
import Data.List qualified as L
import GeneralUtilities
import Text.Read
import Types.Types

-- import           Debug.Trace

-- | allowedCommandList is the permitted command string list
allowedCommandList :: [String]
allowedCommandList = ["build", "fuse", "read", "reblock", "refine", "rename", "report", "run", "search", "select", "set", "support", "swap", "transform"]

-- list of valid instructions
validInstructionList :: [Instruction]
validInstructionList =  [Build, Fuse, Read, Reblock, Refine, Rename, Report, Run, Search, Select, Set, Support, Swap, Transform]

-- | buildArgList is the list of valid build arguments
buildArgList :: [String]
buildArgList = ["atrandom", "best", "block", "character", "cun", "displaytrees", "distance", "dwag", "eun", "first", "graph", "none","nj",
                "otu", "rdwag","return", "replicates", "spr","tbr", "wpgma"]

-- | fuseArgList arguments
fuseArgList :: [String]
fuseArgList = ["all", "atrandom", "best", "joinall", "joinsome", "keep", "once", "pairs", "none", "spr",  "steepest", "tbr", "unique", "reciprocal", "noreciprocal"]

-- | geneticAlgorithm arguments
geneticAlgorithmArgList :: [String]
geneticAlgorithmArgList = ["elitist", "ga", "generations", "geneticalgorithm", "maxnetedges", "popsize", "recombinations", "severity", "stop"]

-- | netEdgeArgList arguments for network edge add/delete operations
netEdgeArgList :: [String]
netEdgeArgList = ["acceptequal", "acceptworse", "all",  "annealing", "atrandom", "drift", "inorder", "keep", "maxchanges", "maxnetedges", "netadd",
    "netadddel", "netadddelete", "netdel", "netdelete", "netmove", "returnmutated", "rounds", "steepest", "steps"]

-- | Read arg list allowable modifiers in read
readArgList :: [String]
readArgList = ["aminoacid", "block", "dot", "enewick", "exclude", "fasta", "fastc", "fenewick", "hugeseq", "include", "newick" , "nucleotide",
                "preaminoacid", "prefasta", "prefastc", "prehugeseq", "prenucleotide", "rename", "tcm", "tnt"] -- "prealigned", "csv",

-- should be moved to a single file for import
-- | reconcileCommandList list of allowable commands
reconcileArgList :: [String]
reconcileArgList = ["compare", "connect", "edgelabel", "method", "outformat", "threshold", "vertexlabel"] -- "outfile"

-- | reconcileOptionsList list of allowable command options of method, compare, threshhold, and outformat
reconcileOptionsList :: [String]
reconcileOptionsList = ["adams", "combinable", "cun", "dot" ,"dotpdf", "eun", "false", "fen", "identity", "majority", "newick", "strict", "true"]


-- | refinement arguments
refineArgList :: [String]
refineArgList = fuseArgList <> netEdgeArgList <> geneticAlgorithmArgList

-- | reportArgList contains valid 'report' arguments
reportArgList :: [String]
reportArgList = ["append", "ascii", "branchlengths", "collapse", "concatenate", "crossrefs", "data", "diagnosis", "displaytrees",
    "dot", "dotpdf", "graphs", "htulabels", "ia", "includemissing", "impliedalignment", "newick", "nobranchlengths", "nocollapse",
    "nohtulabels","overwrite", "pairdist", "reconcile", "search", "support", "tnt"] <> reconcileArgList

-- | search arguments
searchArgList :: [String]
searchArgList = ["days", "exponential", "hours", "instances", "linear", "maxnetedges", "minutes", "mfactor", "seconds", "simple",
    "stop", "thompson"]

-- | selectArgList is the list of valid select arguments
selectArgList :: [String]
selectArgList = ["all", "atrandom", "best", "threshold", "unique"]

-- | setArgList contains valid 'set' arguments
    -- joinThreshold and dynamicEpsilon are not intended for users--but could be of course
setArgList :: [String]
setArgList = ["bc2", "bc4", "bc5", "bc8", "bc64", "bcgt64", "compressresolutions", "criterion", "defparstrat", "dynamicepsilon", "finalassignment", "graphfactor", "graphssteepest", "graphtype", "jointhreshold", "lazyparstrat", "missingthreshold", "modelcomplexity", "multitraverse", "outgroup", "partitioncharacter", "reportnaive", "rootcost", "seed", "softwiredmethod", "strictparstrat", "usenetaddheuristic", "useia"]

-- | refinement arguments
supportArgList :: [String]
supportArgList = ["atrandom","bootstrap", "buildonly", "gb", "gbsample", "goodmanbremer", "jackknife", "replicates", "spr", "tbr"]

-- | swapArgList is the list of valid swap arguments
swapArgList :: [String]
swapArgList = ["acceptequal", "acceptworse", "all", "alternate", "annealing", "atrandom", "drift", "ia", "inorder", "joinall", "joinpruned", "joinalternate", "keep", "maxchanges",
    "nni", "replicates", "returnmutated", "spr", "steepest", "steps", "tbr"]

-- | transform arguments
transformArgList :: [String]
transformArgList = ["atrandom", "compressresolutions", "displaytrees", "dynamic", "dynamicepsilon", "first", "graphfactor",
    "graphssteepest", "jointhreshold", "multitraverse", "name", "outgroup", "softwiredmethod", "staticapprox", "tohardwired",
    "tosoftwired", "totree", "type", "usenetaddheuristic", "weight"]


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
        if not checkInstruction then errorWithoutStackTrace ("Invalid command was specified : " <> (show commandInstruction))
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
                                        if not has02DoubleQuotes then errorWithoutStackTrace ("Unbalanced quotation marks in 'read' file argument: " <> fileArgs)
                                        else (checkCommandArgs "read"  fstArgList readArgList, fileNameList, [""])

                                   -- Reblock -- no arguments--but reads string and blocknames
                                   else if commandInstruction == Reblock then
                                        let fileArgs = concat $ filter (/= []) $ fmap snd inArgs
                                            numDoubleQuotes = length $ filter (== '"') fileArgs
                                            (numDouble, numUnbalanced) = divMod numDoubleQuotes 2
                                        in
                                        -- trace (show (fstArgList, sndArgList,fileNameList)) (
                                        if numDouble < 2 then errorWithoutStackTrace ("Need at least two fields in 'rebock' command, new block name and old block(s) in double quotes: " <> fileArgs)
                                        else if numUnbalanced /= 0 then errorWithoutStackTrace ("Unbalanced quotation marks in 'reblock' command: " <> fileArgs)
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
                                        if numDouble < 2 then errorWithoutStackTrace ("Need at least two fields in 'rename' command, new taxon name and old taxon name(s) in double quotes: " <> fileArgs)
                                        else if numUnbalanced /= 0 then errorWithoutStackTrace ("Unbalanced quotation marks in 'rename' command: " <> fileArgs)
                                        else (True,[""], [""])
                                        -- )

                                   -- Report
                                   else if commandInstruction == Report then
                                        let fileArgs = concat $ filter (/= []) $ fmap snd inArgs
                                            numDoubleQuotes = length $ filter (== '"') fileArgs
                                            has02DoubleQuotes = (numDoubleQuotes == 2) || (numDoubleQuotes == 0)
                                        in
                                        if not has02DoubleQuotes then errorWithoutStackTrace ("Unbalanced quotation marks in report file argument: " <> fileArgs)
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

                                   else errorWithoutStackTrace ("Unrecognized command was specified : " <> (show commandInstruction))



            in
            if checkOptions then
                let allFilesToReadFrom = filter (/= "") $ filesToReadFrom <> inFilesToRead
                    allFilesToWriteTo = filter (/= "") $ filesToWriteTo <> inFilesToWrite
                    readAndWriteFileList = L.intersect allFilesToReadFrom allFilesToWriteTo
                in
                -- trace (show (allFilesToReadFrom, allFilesToWriteTo)) (
                if (not .null) readAndWriteFileList then
                    errorWithoutStackTrace ("Error--Both reading from and writing to files (could cause errors and/or loss of data): " <> (show readAndWriteFileList))
                else verifyCommands (tail inCommandList) allFilesToReadFrom allFilesToWriteTo
                -- )

            else
                -- Won't get to here--will error at earlier stages
                False

        where isInt a = if (readMaybe a :: Maybe Int) /= Nothing then True else False

