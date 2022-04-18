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
geneticAlgorithmArgList = ["popsize", "generations", "elitist", "severity", "recombinations","geneticalgorithm", "ga"]

-- | netEdgeArgList arguments for network edge add/delete operations
netEdgeArgList :: [String]
netEdgeArgList = ["keep", "steepest", "all", "netadd", "netdel", "netdelete", "netadddel", "netmove", "annealing", "steps", "returnmutated", "drift", "acceptequal", "acceptworse", "maxchanges","steepest","atrandom"]

-- | Read arg list allowable modifiers in read
readArgList :: [String]
readArgList = ["tcm", "prealigned", "nucleotide", "aminoacid", "custom_alphabet", "fasta", "fastc", "tnt", "csv",
    "dot", "newick" , "enewick", "fenewick", "terminals", "include", "exclude", "rename", "block", "prefasta", 
    "prefastc", "preaminoacid", "prenucleotide", "precustom_alphabet"]

-- should be moved to a single file for import
-- | reconcileCommandList list of allowable commands
reconcileArgList :: [String]
reconcileArgList = ["method", "compare", "threshold", "outformat", "connect", "edgelabel", "vertexlabel"] -- "outfile" 

-- | reconcileOptionsList list of allowable command options of method, compare, threshhold, and outformat
reconcileOptionsList :: [String]
reconcileOptionsList = ["eun", "cun", "strict", "majority", "adams", "dot" ,"dotpdf", "fen", "newick", "true", "false", "combinable", "identity"]


-- | refinement arguments
refineArgList :: [String]
refineArgList = ["netadd", "netdel", "netdelete", "netadddel", "netmove","geneticalgorithm", "ga"] ++ fuseArgList ++ netEdgeArgList ++ geneticAlgorithmArgList

-- | reportArgList contains valid 'report' arguments
reportArgList :: [String]
reportArgList = ["all", "data", "search", "graphs", "overwrite", "append", "dot", "dotpdf", "newick", "ascii", "crossrefs", "pairdist", "diagnosis","displaytrees", "reconcile", "support", "ia", "impliedalignment", "tnt"] ++ reconcileArgList

-- | search arguments
searchArgList :: [String]
searchArgList = ["days", "hours", "minutes", "seconds", "instances"]

-- | buildArgList is the list of valid build arguments
selectArgList :: [String]
selectArgList = ["best", "all", "unique", "atrandom"]

-- | setArgLIst contains valid 'set' arguments
setArgList :: [String]
setArgList = ["outgroup", "criterion", "graphtype", "compressresolutions", "finalassignment", "graphfactor", "rootcost", "seed"]

-- | refinement arguments
supportArgList :: [String]
supportArgList = ["bootstrap", "jackknife", "goodmanbremer", "gb", "gbsample", "replicates", "buildonly", "atrandom"] -- "bootstrap", 

-- | buildArgList is the list of valid build arguments
swapArgList :: [String]
swapArgList = ["spr","tbr", "keep", "steepest", "all", "nni", "ia", "annealing", "maxtemp", "mintemp", "steps", "returnmutated", "drift", "acceptequal", "acceptworse", "maxchanges"]


-- | makeArgsLowerCase reformats a command
makeArgsLowerCase :: Command -> Command
makeArgsLowerCase inCommand =
    let argList = snd inCommand
        commandList = fmap (fmap C.toLower) $ filter (/= "") $ fmap fst argList
        optionList = fmap (fmap C.toLower) $ filter (/= "") $ fmap snd argList
    in
    (fst inCommand, zip commandList optionList)
     
-- | verifyCommands takes a command list and tests whether the commands 
-- and arguments are permissible before program execution--prevents late failure 
-- after alot of processing time.
-- bit does not check for files existance, write/read ability, or contents for format or 
-- anyhhting else for that matter
-- does check if files are both read from and written to
verifyCommands :: [Command] -> Bool
verifyCommands inCommandList = 
    let instructionList = fmap fst inCommandList

        -- check valid instructions
        -- this is done earlier but might get oved so putting here just in case
        checkInstructionList = fmap (`elem` validInstructionList) instructionList
        instructionInvalidList = fmap fst $ filter ((== False) . snd) $ zip instructionList checkInstructionList

    in
    if (not . null) instructionInvalidList then errorWithoutStackTrace ("Invalid commands were specified : " ++ (show instructionInvalidList))
    else 
        -- check each command for valid arguments
        -- make lower-case arguments
        let newCommandList = fmap makeArgsLowerCase inCommandList

            -- Read
            readPairList = fmap snd $ filter ((== Read) .fst) newCommandList
            (readCommandList, readFileList) = unzip $ concat readPairList
            checkRead = checkCommandArgs "read"  (filter (/= []) readCommandList) readArgList


            -- Build
            buildPairList = fmap snd $ filter ((== Build) .fst) newCommandList
            (buildCommandList, _) = unzip $ concat buildPairList
            checkBuild = checkCommandArgs "build"  (filter (/= []) buildCommandList) buildArgList
        
            -- Fuse
            fusePairList = fmap snd $ filter ((== Fuse) .fst) newCommandList
            (fuseCommandList, _) = unzip $ concat fusePairList
            checkFuse = checkCommandArgs "fuse"  (filter (/= []) fuseCommandList) fuseArgList

            -- Reblock -- part of read

            -- Reconcile -- part of report
            
            -- Refine
            refinePairList = fmap snd $ filter ((== Refine) .fst) newCommandList
            (refineCommandList, _) = unzip $ concat refinePairList
            checkRefine = checkCommandArgs "refine"  (filter (/= []) refineCommandList) refineArgList

            -- Rename -- part of read

            -- Report
            reportPairList = fmap snd $ filter ((== Report) .fst) newCommandList
            (reportCommandList, reportFileList) = unzip $ concat reportPairList
            checkReport = checkCommandArgs "report"  (filter (/= []) reportCommandList) reportArgList

            -- check reconcile options 
            reconcilePairList = filter ((`elem` reconcileArgList). fst) $ concat reportPairList
            nonThresholdreconcileModPairList = filter ((/= "threshold"). fst) $ reconcilePairList
            thresholdreconcileModPairList = filter ((== "threshold"). fst) $ reconcilePairList
            checkReconcile1 = checkCommandArgs "reconcile"  (filter (/= []) (fmap fst nonThresholdreconcileModPairList)) reconcileArgList
            checkReconcile2 = checkCommandArgs "reconcile"  (filter (/= []) (fmap fst thresholdreconcileModPairList)) reconcileArgList
            checkReconcile3 = checkCommandArgs "reconcile modifier (method, compare, outformat, connect, edgelabel, vertexlabel)"  (filter (/= []) (fmap snd nonThresholdreconcileModPairList)) reconcileOptionsList
            checkReconcile4 = L.foldl1' (&&) $ fmap isInt (filter (/= []) (fmap snd thresholdreconcileModPairList))
            checkReconcile = checkReconcile1 && checkReconcile2 && checkReconcile3 && checkReconcile4
            -- reconcileOptionsList

            -- Run  -- processed out before this into command list

            -- Search
            searchPairList = fmap snd $ filter ((== Search) .fst) newCommandList
            (searchCommandList, _) = unzip $ concat searchPairList
            checkSearch = checkCommandArgs "search"  (filter (/= []) searchCommandList) searchArgList


            -- Select
            selectPairList = fmap snd $ filter ((== Select) .fst) newCommandList
            (selectCommandList, _) = unzip $ concat selectPairList
            checkSelect = checkCommandArgs "select"  (filter (/= []) selectCommandList) selectArgList

            -- Set
            setPairList = fmap snd $ filter ((== Set) .fst) newCommandList
            (setCommandList, _) = unzip $ concat setPairList
            checkSet = checkCommandArgs "set"  (filter (/= []) setCommandList) setArgList

            -- Support
            supportPairList = fmap snd $ filter ((== Support) .fst) newCommandList
            (supportCommandList, _) = unzip $ concat supportPairList
            checkSupport = checkCommandArgs "support"  (filter (/= []) supportCommandList) supportArgList

            -- Swap
            swapPairList = fmap snd $ filter ((== Swap) .fst) newCommandList
            (swapCommandList, _) = unzip $ concat swapPairList
            checkSwap = checkCommandArgs "swap"  (filter (/= []) swapCommandList) swapArgList

        in
        if not checkReconcile4 then errorWithoutStackTrace ("Reconcile modifier 'threshold' requires and integer option (e.g. 51): " ++ (show (concat $ fmap snd thresholdreconcileModPairList)))
        else if checkRead && checkBuild && checkFuse && checkReconcile && checkReport && checkRefine && checkSearch && checkSelect && checkSet && checkSupport && checkSwap then
            let filesToReadFrom = readFileList
                filesToWriteTo  = reportFileList 
                readAndWriteFileList = L.intersect filesToReadFrom filesToWriteTo
            in
            if (not .null) readAndWriteFileList then 
                errorWithoutStackTrace ("Error--Both reading from and writing to files (could cause errors and/or loss of data): " ++ (show readAndWriteFileList))
            else True

        else 
            -- Won't get to here--will error at earlier stages
            False
        
        where noQuotes a = if '"' `notElem` a then True else False
              isInt a = if (readMaybe a :: Maybe Int) /= Nothing then True else False

        