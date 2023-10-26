{-# LANGUAGE CPP #-}

{- |
Module      :  CommandExecution.hs
Description :  Module to coordinate command execution
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
module Commands.CommandExecution (
    executeCommands,
    executeRenameReblockCommands,
    getDataListList,
) where

import Commands.CommandUtilities
import Commands.Transform qualified as TRANS
import Commands.Verify qualified as VER
import Control.Evaluation
import Control.Monad (when)
import Control.Monad.IO.Class (MonadIO (..))
import Control.Monad.Logger (LogLevel (..), Logger (..), Verbosity (..))
import Data.CSV qualified as CSV
import Data.Char
import Data.Char qualified as C
import Data.InfList qualified as IL
import Data.List qualified as L
import Data.List.Split qualified as SL
import Data.Maybe
import Data.Ord
import Data.Text.Lazy qualified as T
import Data.Vector qualified as V
import Data.Version qualified as DV
import GeneralUtilities
import GraphOptimization.Traversals qualified as TRAV
import Graphs.GraphOperations qualified as GO
import ParallelUtilities qualified as PU
import Reconciliation.ReconcileGraphs qualified as R
import Search.Build qualified as B
import Search.Refinement qualified as REF
import Search.Search qualified as S
import Support.Support qualified as SUP
import System.CPU qualified as SC
import System.ErrorPhase (ErrorPhase (..))
import System.IO
import System.IO.Unsafe qualified as SIOU
import System.Info qualified as SI
import System.Timing
import Text.Read
import Types.Types
import Utilities.Utilities qualified as U


{- | executeCommands reads input files and returns raw data
need to close files after read
-}
executeCommands
    ∷ GlobalSettings
    → ([NameText], [(NameText, NameText)])
    → Int
    → String
    → ProcessedData
    → ProcessedData
    → ProcessedData
    → [ReducedPhylogeneticGraph]
    → [[VertexCost]]
    → [Int]
    → [ReducedPhylogeneticGraph]
    → [Command]
    → PhyG ([ReducedPhylogeneticGraph], GlobalSettings, [Int], [ReducedPhylogeneticGraph])
executeCommands globalSettings excludeRename numInputFiles crossReferenceString origProcessedData processedData reportingData curGraphs pairwiseDist seedList supportGraphList commandList = do
    if null commandList
        then pure (curGraphs, globalSettings, seedList, supportGraphList)
        else do
            let (firstOption, firstArgs) = head commandList

            -- skip "Read" and "Rename "commands already processed
            if firstOption == Read
                then error ("Read command should already have been processed: " <> show (firstOption, firstArgs))
                else
                    if firstOption == Rename
                        then error ("Rename command should already have been processed: " <> show (firstOption, firstArgs))
                        else
                            if firstOption == Reblock
                                then error ("Reblock command should already have been processed: " <> show (firstOption, firstArgs))
                                else
                                    if firstOption == Run
                                        then error ("Run command should already have been processed: " <> show (firstOption, firstArgs))
                                        else -- other commands

                                            if firstOption == Build
                                                then do
                                                    (elapsedSeconds, newGraphList') ←
                                                        timeOp $ pure $ B.buildGraph firstArgs globalSettings processedData pairwiseDist (head seedList)

                                                    newGraphList <- newGraphList'
                                                    let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
                                                    let newSearchData = searchInfo : searchData globalSettings

                                                    executeCommands
                                                        (globalSettings{searchData = newSearchData})
                                                        excludeRename
                                                        numInputFiles
                                                        crossReferenceString
                                                        origProcessedData
                                                        processedData
                                                        reportingData
                                                        (curGraphs <> newGraphList)
                                                        pairwiseDist
                                                        (tail seedList)
                                                        supportGraphList
                                                        (tail commandList)
                                                else
                                                    if firstOption == Refine
                                                        then do
                                                            (elapsedSeconds, newGraphList') ←
                                                                timeOp $ pure $ REF.refineGraph firstArgs globalSettings processedData (head seedList) curGraphs

                                                            newGraphList <- newGraphList'
                                                            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
                                                            let newSearchData = searchInfo : searchData globalSettings

                                                            executeCommands
                                                                (globalSettings{searchData = newSearchData})
                                                                excludeRename
                                                                numInputFiles
                                                                crossReferenceString
                                                                origProcessedData
                                                                processedData
                                                                reportingData
                                                                newGraphList
                                                                pairwiseDist
                                                                (tail seedList)
                                                                supportGraphList
                                                                (tail commandList)
                                                        else
                                                            if firstOption == Fuse
                                                                then do
                                                                    (elapsedSeconds, newGraphList) ←
                                                                        timeOp $ REF.fuseGraphs firstArgs globalSettings processedData (head seedList) curGraphs

                                                                    let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
                                                                    let newSearchData = searchInfo : searchData globalSettings

                                                                    executeCommands
                                                                        (globalSettings{searchData = newSearchData})
                                                                        excludeRename
                                                                        numInputFiles
                                                                        crossReferenceString
                                                                        origProcessedData
                                                                        processedData
                                                                        reportingData
                                                                        newGraphList
                                                                        pairwiseDist
                                                                        (tail seedList)
                                                                        supportGraphList
                                                                        (tail commandList)
                                                                else
                                                                    if firstOption == Report
                                                                        then do
                                                                            let doDotPDF = elem "dotpdf" $ fmap (fmap toLower . fst) firstArgs
                                                                            let collapse' = elem "collapse" $ fmap (fmap toLower . fst) firstArgs
                                                                            let noCollapse' = elem "nocollapse" $ fmap (fmap toLower . fst) firstArgs

                                                                            -- set default collapse for dotPDF to True, False otherwise
                                                                            let collapse
                                                                                    | collapse' = True
                                                                                    | noCollapse' = False
                                                                                    | doDotPDF = True
                                                                                    | otherwise = False

                                                                            let curGraphs' =
                                                                                    if not collapse
                                                                                        then curGraphs
                                                                                        else fmap U.collapseReducedGraph curGraphs

                                                                            -- use 'temp' updated graphs s don't repeatedly add model and root complexityies
                                                                            -- reporting collapsed
                                                                            -- reverse sorting graphs by cost
                                                                            let rediagnoseWithReportingData = optimalityCriterion globalSettings == NCM && U.has4864PackedChars (thd3 processedData)
                                                                            let graphsWithUpdatedCosts =
                                                                                    L.sortOn
                                                                                        (Data.Ord.Down . snd5)
                                                                                        (TRAV.updateGraphCostsComplexities globalSettings reportingData processedData rediagnoseWithReportingData curGraphs')
                                                                            reportStuff@(reportString, outFile, writeMode) ←
                                                                                reportCommand
                                                                                    globalSettings
                                                                                    firstArgs
                                                                                    excludeRename
                                                                                    numInputFiles
                                                                                    crossReferenceString
                                                                                    reportingData
                                                                                    graphsWithUpdatedCosts
                                                                                    supportGraphList
                                                                                    pairwiseDist

                                                                            if null reportString
                                                                                then do
                                                                                    executeCommands
                                                                                        globalSettings
                                                                                        excludeRename
                                                                                        numInputFiles
                                                                                        crossReferenceString
                                                                                        origProcessedData
                                                                                        processedData
                                                                                        reportingData
                                                                                        curGraphs
                                                                                        pairwiseDist
                                                                                        seedList
                                                                                        supportGraphList
                                                                                        (tail commandList)
                                                                                else do
                                                                                    logWith LogInfo ("Report writing to " <> outFile)

                                                                                    if doDotPDF
                                                                                        then do
                                                                                            let reportString' = changeDotPreamble "digraph {" "digraph G {\n\trankdir = LR;\tnode [ shape = none];\n" reportString
                                                                                            printGraphVizDot reportString' outFile
                                                                                            executeCommands
                                                                                                globalSettings
                                                                                                excludeRename
                                                                                                numInputFiles
                                                                                                crossReferenceString
                                                                                                origProcessedData
                                                                                                processedData
                                                                                                reportingData
                                                                                                curGraphs
                                                                                                pairwiseDist
                                                                                                seedList
                                                                                                supportGraphList
                                                                                                (tail commandList)
                                                                                        else do
                                                                                            if outFile == "stderr"
                                                                                                then liftIO $ hPutStr stderr reportString
                                                                                                else
                                                                                                    if outFile == "stdout"
                                                                                                        then liftIO $ putStr reportString
                                                                                                        else
                                                                                                            if writeMode == "overwrite"
                                                                                                                then liftIO $ writeFile outFile reportString
                                                                                                                else
                                                                                                                    if writeMode == "append"
                                                                                                                        then liftIO $ appendFile outFile reportString
                                                                                                                        else failWithPhase Parsing ("Error 'read' command not properly formatted" <> show reportStuff)
                                                                                            executeCommands
                                                                                                globalSettings
                                                                                                excludeRename
                                                                                                numInputFiles
                                                                                                crossReferenceString
                                                                                                origProcessedData
                                                                                                processedData
                                                                                                reportingData
                                                                                                curGraphs
                                                                                                pairwiseDist
                                                                                                seedList
                                                                                                supportGraphList
                                                                                                (tail commandList)
                                                                        else
                                                                            if firstOption == Search
                                                                                then do
                                                                                    (elapsedSeconds, output) ←
                                                                                        timeOp $ S.search firstArgs globalSettings processedData pairwiseDist (head seedList) curGraphs
                                                                                    -- in pure result
                                                                                    -- (newGraphList, serchInfoList) <- S.search firstArgs globalSettings origProcessedData processedData reportingDatapairwiseDist (head seedList) curGraphs
                                                                                    let searchInfo =
                                                                                            makeSearchRecord
                                                                                                firstOption
                                                                                                firstArgs
                                                                                                curGraphs
                                                                                                (fst output)
                                                                                                (fromIntegral $ toMilliseconds elapsedSeconds)
                                                                                                (concatMap (L.intercalate "\n") (snd output))
                                                                                    let newSearchData = searchInfo : searchData globalSettings
                                                                                    executeCommands
                                                                                        (globalSettings{searchData = newSearchData})
                                                                                        excludeRename
                                                                                        numInputFiles
                                                                                        crossReferenceString
                                                                                        origProcessedData
                                                                                        processedData
                                                                                        reportingData
                                                                                        (fst output)
                                                                                        pairwiseDist
                                                                                        (tail seedList)
                                                                                        supportGraphList
                                                                                        (tail commandList)
                                                                                else
                                                                                    if firstOption == Select
                                                                                        then do
                                                                                            (elapsedSeconds, newGraphList) ← timeOp $ pure $ GO.selectPhylogeneticGraphReduced firstArgs (head seedList) curGraphs

                                                                                            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
                                                                                            let newSearchData = searchInfo : searchData globalSettings

                                                                                            let typeSelected =
                                                                                                    if null firstArgs
                                                                                                        then "best"
                                                                                                        else fmap C.toLower $ fst $ head firstArgs

                                                                                            logWith LogInfo ("Selecting " <> typeSelected <> " graphs")
                                                                                            executeCommands
                                                                                                (globalSettings{searchData = newSearchData})
                                                                                                excludeRename
                                                                                                numInputFiles
                                                                                                crossReferenceString
                                                                                                origProcessedData
                                                                                                processedData
                                                                                                reportingData
                                                                                                newGraphList
                                                                                                pairwiseDist
                                                                                                (tail seedList)
                                                                                                supportGraphList
                                                                                                (tail commandList)
                                                                                        else
                                                                                            if firstOption == Set
                                                                                                then do
                                                                                                    -- if set changes graph aspects--may nned to reoptimize
                                                                                                    (newGlobalSettings, newProcessedData, seedList') ← setCommand firstArgs globalSettings origProcessedData processedData seedList
                                                                                                    let newGraphList =
                                                                                                            if not (requireReoptimization globalSettings newGlobalSettings)
                                                                                                                then curGraphs
                                                                                                                else -- logWith LogInfo "Reoptimizing gaphs"
                                                                                                                    fmap (TRAV.multiTraverseFullyLabelGraphReduced newGlobalSettings newProcessedData True True Nothing) (fmap fst5 curGraphs)

                                                                                                    let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList 0 "No Comment"
                                                                                                    let newSearchData = searchInfo : searchData newGlobalSettings

                                                                                                    if not (requireReoptimization globalSettings newGlobalSettings)
                                                                                                        then do logWith LogInfo "No need to reoptimize graphs"
                                                                                                        else do logWith LogInfo "Reoptimizing gaphs"

                                                                                                    executeCommands
                                                                                                        (newGlobalSettings{searchData = newSearchData})
                                                                                                        excludeRename
                                                                                                        numInputFiles
                                                                                                        crossReferenceString
                                                                                                        origProcessedData
                                                                                                        processedData
                                                                                                        reportingData
                                                                                                        newGraphList
                                                                                                        pairwiseDist
                                                                                                        seedList'
                                                                                                        supportGraphList
                                                                                                        (tail commandList)
                                                                                                else
                                                                                                    if firstOption == Swap
                                                                                                        then do
                                                                                                            (elapsedSeconds, newGraphList) ←
                                                                                                                timeOp $ REF.swapMaster firstArgs globalSettings processedData (head seedList) curGraphs

                                                                                                            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
                                                                                                            let newSearchData = searchInfo : searchData globalSettings

                                                                                                            executeCommands
                                                                                                                (globalSettings{searchData = newSearchData})
                                                                                                                excludeRename
                                                                                                                numInputFiles
                                                                                                                crossReferenceString
                                                                                                                origProcessedData
                                                                                                                processedData
                                                                                                                reportingData
                                                                                                                newGraphList
                                                                                                                pairwiseDist
                                                                                                                (tail seedList)
                                                                                                                supportGraphList
                                                                                                                (tail commandList)
                                                                                                        else
                                                                                                            if firstOption == Support
                                                                                                                then do
                                                                                                                    (elapsedSeconds, newSupportGraphList') ←
                                                                                                                        timeOp $ pure $ SUP.supportGraph firstArgs globalSettings processedData (head seedList) curGraphs

                                                                                                                    newSupportGraphList <- newSupportGraphList'
                                                                                                                    let searchInfo =
                                                                                                                            makeSearchRecord firstOption firstArgs curGraphs newSupportGraphList (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
                                                                                                                    let newSearchData = searchInfo : searchData globalSettings

                                                                                                                    executeCommands
                                                                                                                        (globalSettings{searchData = newSearchData})
                                                                                                                        excludeRename
                                                                                                                        numInputFiles
                                                                                                                        crossReferenceString
                                                                                                                        origProcessedData
                                                                                                                        processedData
                                                                                                                        reportingData
                                                                                                                        curGraphs
                                                                                                                        pairwiseDist
                                                                                                                        (tail seedList)
                                                                                                                        (supportGraphList <> newSupportGraphList)
                                                                                                                        (tail commandList)
                                                                                                                else
                                                                                                                    if firstOption == Transform
                                                                                                                        then do
                                                                                                                            (elapsedSeconds, (newGS, newOrigData, newProcessedData, newGraphs)) ←
                                                                                                                                timeOp $ TRANS.transform firstArgs globalSettings origProcessedData processedData (head seedList) curGraphs

                                                                                                                            let searchInfo = makeSearchRecord firstOption firstArgs curGraphs newGraphs (fromIntegral $ toMilliseconds elapsedSeconds) "No Comment"
                                                                                                                            let newSearchData = searchInfo : searchData globalSettings

                                                                                                                            executeCommands
                                                                                                                                (newGS{searchData = newSearchData})
                                                                                                                                excludeRename
                                                                                                                                numInputFiles
                                                                                                                                crossReferenceString
                                                                                                                                newOrigData
                                                                                                                                newProcessedData
                                                                                                                                reportingData
                                                                                                                                newGraphs
                                                                                                                                pairwiseDist
                                                                                                                                (tail seedList)
                                                                                                                                supportGraphList
                                                                                                                                (tail commandList)
                                                                                                                        else error ("Command " <> show firstOption <> " not recognized/implemented")


-- | makeSearchRecord take sbefore and after data of a commend and returns SearchData record
makeSearchRecord
    ∷ Instruction → [Argument] → [ReducedPhylogeneticGraph] → [ReducedPhylogeneticGraph] → Int → String → SearchData
makeSearchRecord firstOption firstArgs curGraphs newGraphList elapsedTime comment =
    SearchData
        { instruction = firstOption
        , arguments = firstArgs
        , minGraphCostIn =
            if null curGraphs
                then infinity
                else minimum $ fmap snd5 curGraphs
        , maxGraphCostIn =
            if null curGraphs
                then infinity
                else maximum $ fmap snd5 curGraphs
        , numGraphsIn = length curGraphs
        , minGraphCostOut =
            if null newGraphList
                then infinity
                else minimum $ fmap snd5 newGraphList
        , maxGraphCostOut =
            if null newGraphList
                then infinity
                else maximum $ fmap snd5 newGraphList
        , numGraphsOut = length newGraphList
        , commentString = comment
        , duration = elapsedTime
        }


{- | setCommand takes arguments to change globalSettings and multiple data aspects (e.g. 'blocks')
needs to be abtracted--too long
if seed list is empty [] then processes first set--confusing--shold be refactored
-}
setCommand
    ∷ [Argument] → GlobalSettings → ProcessedData → ProcessedData → [Int] → PhyG (GlobalSettings, ProcessedData, [Int])
setCommand argList globalSettings origProcessedData processedData inSeedList =
    let commandList = fmap (fmap C.toLower) $ filter (/= "") $ fmap fst argList
        optionList = fmap (fmap C.toLower) $ filter (/= "") $ fmap snd argList
        checkCommandList = checkCommandArgs "set" commandList VER.setArgList
        leafNameVect = fst3 processedData
    in  if not checkCommandList
            then do failWithPhase Parsing ("Unrecognized command in 'set': " <> show argList)
            else -- this could be changed later

                if length commandList > 1 || length optionList > 1
                    then do
                        failWithPhase Parsing ("Set option error: can only have one set argument for each command: " <> show (commandList, optionList))
                    else -- early extraction of partition character and bc2-gt64 follows from null inputs
                    -- this due to not having all info required for all global settings, so options restricted and repeated
                    -- needs to be fixed to be more clear and clean

                        if null inSeedList
                            then
                                if head commandList == "partitioncharacter"
                                    then
                                        let localPartitionChar = head optionList
                                        in  if length localPartitionChar /= 1
                                                then do
                                                    failWithPhase
                                                        Parsing
                                                        ("Error in 'set' command. Partitioncharacter '" <> show localPartitionChar <> "' must be a single character")
                                                else do
                                                    logWith LogInfo ("PartitionCharacter set to '" <> head optionList <> "'")
                                                    pure (globalSettings{partitionCharacter = localPartitionChar}, processedData, inSeedList)
                                    else
                                        if head commandList == "missingthreshold"
                                            then
                                                let localValue = readMaybe (head optionList) ∷ Maybe Int
                                                in  if isNothing localValue
                                                        then do
                                                            failWithPhase
                                                                Parsing
                                                                ("Set option 'missingThreshold' must be set to an integer value (e.g. missingThreshold:50): " <> head optionList)
                                                        else
                                                            if (fromJust localValue < 0) || (fromJust localValue > 100)
                                                                then do
                                                                    failWithPhase Parsing ("Set option 'missingThreshold' must be set to an integer value between 0 and 100: " <> head optionList)
                                                                else do
                                                                    logWith LogInfo ("MissingThreshold set to " <> head optionList)
                                                                    pure (globalSettings{missingThreshold = fromJust localValue}, processedData, inSeedList)
                                            else -- sets root cost as well-- need in both places--one to process data and one to
                                            -- keep in current global

                                                if head commandList == "criterion" -- (globalSettings, processedData, inSeedList)
                                                    then
                                                        let localCriterion
                                                                | (head optionList == "parsimony") = Just Parsimony
                                                                | (head optionList == "pmdl") = Just PMDL
                                                                | (head optionList == "si") = Just SI
                                                                | (head optionList == "mapa") = Just MAPA
                                                                | (head optionList == "ncm") = Just NCM
                                                                | otherwise = Nothing

                                                            -- create lazy list of graph complexity indexed by number of network nodes--need leaf number for base tree complexity
                                                            lGraphComplexityList
                                                                | localCriterion `elem` [Just Parsimony, Just NCM] = Just $ IL.repeat (0.0, 0.0)
                                                                | localCriterion `elem` [Just PMDL, Just SI, Just MAPA] = Just $ U.calculateGraphComplexity processedData
                                                                | otherwise = Nothing

                                                            lRootComplexity
                                                                | localCriterion == Just Parsimony = Just 0.0
                                                                | localCriterion `elem` [Just PMDL, Just SI, Just MAPA] = Just $ U.calculateW15RootCost processedData
                                                                | localCriterion == Just NCM =
                                                                    if origProcessedData /= emptyProcessedData
                                                                        then Just $ U.calculateNCMRootCost origProcessedData
                                                                        else Just $ U.calculateNCMRootCost processedData
                                                                | otherwise = Nothing

                                                            lGraphFactor =
                                                                if localCriterion `elem` [Just PMDL, Just SI, Just MAPA]
                                                                    then PMDLGraph
                                                                    else graphFactor globalSettings
                                                        in  if isNothing localCriterion
                                                                then do
                                                                    failWithPhase Parsing ("Error in 'set' command. Criterion '" <> head optionList <> "' is not 'parsimony', 'ml', or 'pmdl'")
                                                                else
                                                                    if isNothing lGraphComplexityList
                                                                        then do
                                                                            failWithPhase Parsing ("Optimality criterion not recognized: " <> show localCriterion)
                                                                        else
                                                                            if isNothing lRootComplexity
                                                                                then do
                                                                                    failWithPhase Parsing ("Optimality criterion not recognized: " <> show localCriterion)
                                                                                else
                                                                                    pure
                                                                                        ( globalSettings
                                                                                            { optimalityCriterion = fromJust localCriterion
                                                                                            , graphComplexityList = fromJust lGraphComplexityList
                                                                                            , rootComplexity = fromJust lRootComplexity
                                                                                            , graphFactor = lGraphFactor
                                                                                            }
                                                                                        , processedData
                                                                                        , inSeedList
                                                                                        )
                                                    else
                                                        if head commandList == "bc2"
                                                            then
                                                                let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                                    noChangeValue = readMaybe noChangeString ∷ Maybe Double
                                                                    changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                                    changeValue = readMaybe changeString ∷ Maybe Double
                                                                in  if length commandList /= length optionList
                                                                        then do
                                                                            failWithPhase Parsing ("Set option error: number of values and options do not match: " <> show (commandList, optionList))
                                                                        else
                                                                            if (null . head) optionList
                                                                                then
                                                                                    failWithPhase
                                                                                        Parsing
                                                                                        "Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): no values found "
                                                                                else
                                                                                    if ',' `notElem` head optionList
                                                                                        then do
                                                                                            failWithPhase
                                                                                                Parsing
                                                                                                "Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): no comma found "
                                                                                        else
                                                                                            if isNothing noChangeValue || isNothing changeValue
                                                                                                then do
                                                                                                    failWithPhase
                                                                                                        Parsing
                                                                                                        ( "Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): "
                                                                                                            <> head optionList
                                                                                                        )
                                                                                                else
                                                                                                    if bc2 globalSettings /= (fromJust noChangeValue, fromJust changeValue)
                                                                                                        then do
                                                                                                            logWith LogInfo ("bit cost 2 state set to " <> show (fromJust noChangeValue, fromJust changeValue))
                                                                                                            pure (globalSettings{bc2 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                                                        else pure (globalSettings{bc2 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                            else
                                                                if head commandList == "bc4"
                                                                    then
                                                                        let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                                            noChangeValue = readMaybe noChangeString ∷ Maybe Double
                                                                            changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                                            changeValue = readMaybe changeString ∷ Maybe Double
                                                                        in  if (null . head) optionList
                                                                                then do
                                                                                    failWithPhase
                                                                                        Parsing
                                                                                        "Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): no values found "
                                                                                else
                                                                                    if ',' `notElem` head optionList
                                                                                        then do
                                                                                            failWithPhase
                                                                                                Parsing
                                                                                                "Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): no comma found "
                                                                                        else
                                                                                            if isNothing noChangeValue || isNothing changeValue
                                                                                                then do
                                                                                                    failWithPhase
                                                                                                        Parsing
                                                                                                        ( "Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): "
                                                                                                            <> head optionList
                                                                                                        )
                                                                                                else
                                                                                                    if bc4 globalSettings /= (fromJust noChangeValue, fromJust changeValue)
                                                                                                        then do
                                                                                                            logWith LogInfo ("bit cost 4 state set to " <> show (fromJust noChangeValue, fromJust changeValue))
                                                                                                            pure (globalSettings{bc4 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                                                        else pure (globalSettings{bc4 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                    else
                                                                        if head commandList == "bc5"
                                                                            then
                                                                                let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                                                    noChangeValue = readMaybe noChangeString ∷ Maybe Double
                                                                                    changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                                                    changeValue = readMaybe changeString ∷ Maybe Double
                                                                                in  if (null . head) optionList
                                                                                        then do
                                                                                            failWithPhase
                                                                                                Parsing
                                                                                                "Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): no values found "
                                                                                        else
                                                                                            if ',' `notElem` head optionList
                                                                                                then do
                                                                                                    failWithPhase
                                                                                                        Parsing
                                                                                                        "Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): no comma found "
                                                                                                else
                                                                                                    if isNothing noChangeValue || isNothing changeValue
                                                                                                        then do
                                                                                                            failWithPhase
                                                                                                                Parsing
                                                                                                                ( "Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): "
                                                                                                                    <> head optionList
                                                                                                                )
                                                                                                        else
                                                                                                            if bc5 globalSettings /= (fromJust noChangeValue, fromJust changeValue)
                                                                                                                then do
                                                                                                                    logWith LogInfo ("bit cost 5 state set to " <> show (fromJust noChangeValue, fromJust changeValue))
                                                                                                                    pure (globalSettings{bc5 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                                                                else pure (globalSettings{bc5 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                            else
                                                                                if head commandList == "bc8"
                                                                                    then
                                                                                        let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                                                            noChangeValue = readMaybe noChangeString ∷ Maybe Double
                                                                                            changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                                                            changeValue = readMaybe changeString ∷ Maybe Double
                                                                                        in  if (null . head) optionList
                                                                                                then do
                                                                                                    failWithPhase
                                                                                                        Parsing
                                                                                                        "Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): no values found "
                                                                                                else
                                                                                                    if ',' `notElem` head optionList
                                                                                                        then do
                                                                                                            failWithPhase
                                                                                                                Parsing
                                                                                                                "Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): no comma found "
                                                                                                        else
                                                                                                            if isNothing noChangeValue || isNothing changeValue
                                                                                                                then do
                                                                                                                    failWithPhase
                                                                                                                        Parsing
                                                                                                                        ( "Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): "
                                                                                                                            <> head optionList
                                                                                                                        )
                                                                                                                else
                                                                                                                    if bc8 globalSettings /= (fromJust noChangeValue, fromJust changeValue)
                                                                                                                        then do
                                                                                                                            logWith LogInfo ("bit cost 8 state set to " <> show (fromJust noChangeValue, fromJust changeValue))
                                                                                                                            pure (globalSettings{bc8 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                                                                        else pure (globalSettings{bc8 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                                    else
                                                                                        if head commandList == "bc64"
                                                                                            then
                                                                                                let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                                                                    noChangeValue = readMaybe noChangeString ∷ Maybe Double
                                                                                                    changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                                                                    changeValue = readMaybe changeString ∷ Maybe Double
                                                                                                in  if (null . head) optionList
                                                                                                        then do
                                                                                                            failWithPhase
                                                                                                                Parsing
                                                                                                                "Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): no values found "
                                                                                                        else
                                                                                                            if ',' `notElem` head optionList
                                                                                                                then do
                                                                                                                    failWithPhase
                                                                                                                        Parsing
                                                                                                                        "Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): no comma found "
                                                                                                                else
                                                                                                                    if isNothing noChangeValue || isNothing changeValue
                                                                                                                        then do
                                                                                                                            failWithPhase
                                                                                                                                Parsing
                                                                                                                                ( "Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): "
                                                                                                                                    <> head optionList
                                                                                                                                )
                                                                                                                        else
                                                                                                                            if bc64 globalSettings /= (fromJust noChangeValue, fromJust changeValue)
                                                                                                                                then do
                                                                                                                                    logWith LogInfo ("bit cost 64 state set to " <> show (fromJust noChangeValue, fromJust changeValue))
                                                                                                                                    pure (globalSettings{bc64 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                                                                                else pure (globalSettings{bc64 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                                            else
                                                                                                if head commandList == "bcgt64"
                                                                                                    then
                                                                                                        let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                                                                            noChangeValue = readMaybe noChangeString ∷ Maybe Double
                                                                                                            changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                                                                            changeValue = readMaybe changeString ∷ Maybe Double
                                                                                                        in  if (null . head) optionList
                                                                                                                then do
                                                                                                                    failWithPhase
                                                                                                                        Parsing
                                                                                                                        "Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): no values found "
                                                                                                                else
                                                                                                                    if ',' `notElem` head optionList
                                                                                                                        then do
                                                                                                                            failWithPhase
                                                                                                                                Parsing
                                                                                                                                "Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): no comma found "
                                                                                                                        else
                                                                                                                            if isNothing noChangeValue || isNothing changeValue
                                                                                                                                then do
                                                                                                                                    failWithPhase
                                                                                                                                        Parsing
                                                                                                                                        ( "Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): "
                                                                                                                                            <> head optionList
                                                                                                                                        )
                                                                                                                                else
                                                                                                                                    if bcgt64 globalSettings /= (fromJust noChangeValue, fromJust changeValue)
                                                                                                                                        then do
                                                                                                                                            logWith LogInfo ("bit cost > 64 state set to " <> show (fromJust noChangeValue, fromJust changeValue))
                                                                                                                                            pure (globalSettings{bcgt64 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                                                                                        else do pure (globalSettings{bcgt64 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                                                    else -- partition character to reset
                                                                                                    do
                                                                                                        -- trace ("PartitionCharacter set to '" <> (partitionCharacter globalSettings) <> "'")
                                                                                                        pure (globalSettings, processedData, inSeedList)
                            else -- regular command stuff not initial at start

                                if head commandList == "bc2"
                                    then
                                        let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                            noChangeValue = readMaybe noChangeString ∷ Maybe Double
                                            changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                            changeValue = readMaybe changeString ∷ Maybe Double
                                        in  if length commandList /= length optionList
                                                then do
                                                    failWithPhase Parsing ("Set option error: number of values and options do not match: " <> show (commandList, optionList))
                                                else
                                                    if (null . head) optionList
                                                        then
                                                            failWithPhase
                                                                Parsing
                                                                "Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): no values found "
                                                        else
                                                            if ',' `notElem` head optionList
                                                                then do
                                                                    failWithPhase
                                                                        Parsing
                                                                        "Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): no comma found "
                                                                else
                                                                    if isNothing noChangeValue || isNothing changeValue
                                                                        then do
                                                                            failWithPhase
                                                                                Parsing
                                                                                ( "Set option 'bc2' must be set to a pair of double values in parens, separated by a comma (e.g. bc2:(0.1, 1.1): "
                                                                                    <> head optionList
                                                                                )
                                                                        else
                                                                            if bc2 globalSettings /= (fromJust noChangeValue, fromJust changeValue)
                                                                                then do
                                                                                    logWith LogInfo ("bit cost 2 state set to " <> show (fromJust noChangeValue, fromJust changeValue))
                                                                                    pure (globalSettings{bc2 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                                else pure (globalSettings{bc2 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                    else
                                        if head commandList == "bc4"
                                            then
                                                let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                    noChangeValue = readMaybe noChangeString ∷ Maybe Double
                                                    changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                    changeValue = readMaybe changeString ∷ Maybe Double
                                                in  if (null . head) optionList
                                                        then do
                                                            failWithPhase
                                                                Parsing
                                                                "Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): no values found "
                                                        else
                                                            if ',' `notElem` head optionList
                                                                then do
                                                                    failWithPhase
                                                                        Parsing
                                                                        "Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): no comma found "
                                                                else
                                                                    if isNothing noChangeValue || isNothing changeValue
                                                                        then do
                                                                            failWithPhase
                                                                                Parsing
                                                                                ( "Set option 'bc4' must be set to a pair of double values in parens, separated by a comma (e.g. bc4:(0.1, 1.1): "
                                                                                    <> head optionList
                                                                                )
                                                                        else
                                                                            if bc4 globalSettings /= (fromJust noChangeValue, fromJust changeValue)
                                                                                then do
                                                                                    logWith LogInfo ("bit cost 4 state set to " <> show (fromJust noChangeValue, fromJust changeValue))
                                                                                    pure (globalSettings{bc4 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                                else pure (globalSettings{bc4 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                            else
                                                if head commandList == "bc5"
                                                    then
                                                        let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                            noChangeValue = readMaybe noChangeString ∷ Maybe Double
                                                            changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                            changeValue = readMaybe changeString ∷ Maybe Double
                                                        in  if (null . head) optionList
                                                                then do
                                                                    failWithPhase
                                                                        Parsing
                                                                        "Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): no values found "
                                                                else
                                                                    if ',' `notElem` head optionList
                                                                        then do
                                                                            failWithPhase
                                                                                Parsing
                                                                                "Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): no comma found "
                                                                        else
                                                                            if isNothing noChangeValue || isNothing changeValue
                                                                                then do
                                                                                    failWithPhase
                                                                                        Parsing
                                                                                        ( "Set option 'bc5' must be set to a pair of double values in parens, separated by a comma (e.g. bc5:(0.1, 1.1): "
                                                                                            <> head optionList
                                                                                        )
                                                                                else
                                                                                    if bc5 globalSettings /= (fromJust noChangeValue, fromJust changeValue)
                                                                                        then do
                                                                                            logWith LogInfo ("bit cost 5 state set to " <> show (fromJust noChangeValue, fromJust changeValue))
                                                                                            pure (globalSettings{bc5 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                                        else pure (globalSettings{bc5 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                    else
                                                        if head commandList == "bc8"
                                                            then
                                                                let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                                    noChangeValue = readMaybe noChangeString ∷ Maybe Double
                                                                    changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                                    changeValue = readMaybe changeString ∷ Maybe Double
                                                                in  if (null . head) optionList
                                                                        then do
                                                                            failWithPhase
                                                                                Parsing
                                                                                "Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): no values found "
                                                                        else
                                                                            if ',' `notElem` head optionList
                                                                                then do
                                                                                    failWithPhase
                                                                                        Parsing
                                                                                        "Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): no comma found "
                                                                                else
                                                                                    if isNothing noChangeValue || isNothing changeValue
                                                                                        then do
                                                                                            failWithPhase
                                                                                                Parsing
                                                                                                ( "Set option 'bc8' must be set to a pair of double values in parens, separated by a comma (e.g. bc8:(0.1, 1.1): "
                                                                                                    <> head optionList
                                                                                                )
                                                                                        else
                                                                                            if bc8 globalSettings /= (fromJust noChangeValue, fromJust changeValue)
                                                                                                then do
                                                                                                    logWith LogInfo ("bit cost 8 state set to " <> show (fromJust noChangeValue, fromJust changeValue))
                                                                                                    pure (globalSettings{bc8 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                                                else pure (globalSettings{bc8 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                            else
                                                                if head commandList == "bc64"
                                                                    then
                                                                        let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                                            noChangeValue = readMaybe noChangeString ∷ Maybe Double
                                                                            changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                                            changeValue = readMaybe changeString ∷ Maybe Double
                                                                        in  if (null . head) optionList
                                                                                then do
                                                                                    failWithPhase
                                                                                        Parsing
                                                                                        "Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): no values found "
                                                                                else
                                                                                    if ',' `notElem` head optionList
                                                                                        then do
                                                                                            failWithPhase
                                                                                                Parsing
                                                                                                "Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): no comma found "
                                                                                        else
                                                                                            if isNothing noChangeValue || isNothing changeValue
                                                                                                then do
                                                                                                    failWithPhase
                                                                                                        Parsing
                                                                                                        ( "Set option 'bc64' must be set to a pair of double values in parens, separated by a comma (e.g. bc64:(0.1, 1.1): "
                                                                                                            <> head optionList
                                                                                                        )
                                                                                                else
                                                                                                    if bc64 globalSettings /= (fromJust noChangeValue, fromJust changeValue)
                                                                                                        then do
                                                                                                            logWith LogInfo ("bit cost 64 state set to " <> show (fromJust noChangeValue, fromJust changeValue))
                                                                                                            pure (globalSettings{bc64 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                                                        else pure (globalSettings{bc64 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                    else
                                                                        if head commandList == "bcgt64"
                                                                            then
                                                                                let noChangeString = takeWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                                                    noChangeValue = readMaybe noChangeString ∷ Maybe Double
                                                                                    changeString = tail $ dropWhile (/= ',') $ filter (`notElem` ['(', ')']) $ head optionList
                                                                                    changeValue = readMaybe changeString ∷ Maybe Double
                                                                                in  if (null . head) optionList
                                                                                        then do
                                                                                            failWithPhase
                                                                                                Parsing
                                                                                                "Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): no values found "
                                                                                        else
                                                                                            if ',' `notElem` head optionList
                                                                                                then do
                                                                                                    failWithPhase
                                                                                                        Parsing
                                                                                                        "Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): no comma found "
                                                                                                else
                                                                                                    if isNothing noChangeValue || isNothing changeValue
                                                                                                        then do
                                                                                                            failWithPhase
                                                                                                                Parsing
                                                                                                                ( "Set option 'bcgt64' must be set to a pair of double values in parens, separated by a comma (e.g. bcgt64:(0.1, 1.1): "
                                                                                                                    <> head optionList
                                                                                                                )
                                                                                                        else
                                                                                                            if bcgt64 globalSettings /= (fromJust noChangeValue, fromJust changeValue)
                                                                                                                then do
                                                                                                                    logWith LogInfo ("bit cost > 64 state set to " <> show (fromJust noChangeValue, fromJust changeValue))
                                                                                                                    pure (globalSettings{bcgt64 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                                                                else pure (globalSettings{bcgt64 = (fromJust noChangeValue, fromJust changeValue)}, processedData, inSeedList)
                                                                            else -- processed above, but need here since put in different value

                                                                                if head commandList == "criterion" -- (globalSettings, processedData, inSeedList)
                                                                                -- {-
                                                                                -- trace ("In Set: " <> (show optionList)) $
                                                                                    then
                                                                                        let localCriterion
                                                                                                | (head optionList == "parsimony") = Just Parsimony
                                                                                                | (head optionList == "pmdl") = Just PMDL
                                                                                                | (head optionList == "si") = Just SI
                                                                                                | (head optionList == "mapa") = Just MAPA
                                                                                                | (head optionList == "ncm") = Just NCM
                                                                                                | otherwise = Nothing

                                                                                            -- create lazy list of graph complexity indexed by number of network nodes--need leaf number for base tree complexity
                                                                                            lGraphComplexityList
                                                                                                | localCriterion `elem` [Just Parsimony, Just NCM] = Just $ IL.repeat (0.0, 0.0)
                                                                                                | localCriterion `elem` [Just PMDL, Just SI, Just MAPA] = Just $ U.calculateGraphComplexity processedData
                                                                                                | otherwise = Nothing

                                                                                            lRootComplexity
                                                                                                | localCriterion == Just Parsimony = Just 0.0
                                                                                                | localCriterion `elem` [Just PMDL, Just SI, Just MAPA] = Just $ U.calculateW15RootCost processedData
                                                                                                | localCriterion == Just NCM =
                                                                                                    if origProcessedData /= emptyProcessedData
                                                                                                        then Just $ U.calculateNCMRootCost origProcessedData
                                                                                                        else Just $ U.calculateNCMRootCost processedData
                                                                                                | otherwise = Nothing

                                                                                            lGraphFactor =
                                                                                                if localCriterion `elem` [Just PMDL, Just SI, Just MAPA]
                                                                                                    then PMDLGraph
                                                                                                    else graphFactor globalSettings
                                                                                        in  if isNothing localCriterion
                                                                                                then do
                                                                                                    failWithPhase Parsing ("Error in 'set' command. Criterion '" <> head optionList <> "' is not 'parsimony', 'ml', or 'pmdl'")
                                                                                                else
                                                                                                    if isNothing lGraphComplexityList
                                                                                                        then do
                                                                                                            failWithPhase Parsing ("Optimality criterion not recognized: " <> show localCriterion)
                                                                                                        else
                                                                                                            if isNothing lRootComplexity
                                                                                                                then do
                                                                                                                    failWithPhase Parsing ("Optimality criterion not recognized: " <> show localCriterion)
                                                                                                                else
                                                                                                                    pure
                                                                                                                        ( globalSettings
                                                                                                                            { optimalityCriterion = fromJust localCriterion
                                                                                                                            , graphComplexityList = fromJust lGraphComplexityList
                                                                                                                            , rootComplexity = fromJust lRootComplexity
                                                                                                                            , graphFactor = lGraphFactor
                                                                                                                            }
                                                                                                                        , processedData
                                                                                                                        , inSeedList
                                                                                                                        )
                                                                                    else -- }

                                                                                    -- modify the behavior of resolutionCache softwired optimization

                                                                                        if head commandList == "compressresolutions"
                                                                                            then
                                                                                                let localCriterion
                                                                                                        | (head optionList == "true") = Just True
                                                                                                        | (head optionList == "false") = Just False
                                                                                                        | otherwise = Nothing
                                                                                                in  if isNothing localCriterion
                                                                                                        then do
                                                                                                            failWithPhase Parsing ("Error in 'set' command. CompressResolutions '" <> head optionList <> "' is not 'true' or 'false'")
                                                                                                        else do
                                                                                                            logWith LogInfo ("CompressResolutions set to " <> head optionList)
                                                                                                            pure (globalSettings{compressResolutions = fromJust localCriterion}, processedData, inSeedList)
                                                                                            else -- this not intended to be for users

                                                                                                if head commandList == "dynamicepsilon"
                                                                                                    then
                                                                                                        let localValue = readMaybe (head optionList) ∷ Maybe Double
                                                                                                        in  if isNothing localValue
                                                                                                                then error ("Set option 'dynamicEpsilon' must be set to an double value >= 0.0 (e.g. dynamicepsilon:0.02): " <> head optionList)
                                                                                                                else
                                                                                                                    if fromJust localValue < 0.0
                                                                                                                        then
                                                                                                                            errorWithoutStackTrace
                                                                                                                                ("Set option 'dynamicEpsilon' must be set to a double value >= 0.0 (e.g. dynamicepsilon:0.02): " <> head optionList)
                                                                                                                        else do
                                                                                                                            logWith LogInfo ("Dynamic Epsilon factor set to " <> head optionList)
                                                                                                                            pure
                                                                                                                                (globalSettings{dynamicEpsilon = 1.0 + (fromJust localValue * fractionDynamic globalSettings)}, processedData, inSeedList)
                                                                                                    else
                                                                                                        if head commandList == "finalassignment"
                                                                                                            then
                                                                                                                let localMethod
                                                                                                                        | ((head optionList == "do") || (head optionList == "directoptimization")) = Just DirectOptimization
                                                                                                                        | ((head optionList == "ia") || (head optionList == "impliedalignment")) = Just ImpliedAlignment
                                                                                                                        | otherwise = Nothing
                                                                                                                in  if isNothing localMethod
                                                                                                                        then do
                                                                                                                            failWithPhase
                                                                                                                                Parsing
                                                                                                                                ( "Error in 'set' command. FinalAssignment  '"
                                                                                                                                    <> head optionList
                                                                                                                                    <> "' is not 'DirectOptimization (DO)' or 'ImpliedAlignment (IA)'"
                                                                                                                                )
                                                                                                                        else
                                                                                                                            if graphType globalSettings == Tree
                                                                                                                                then do
                                                                                                                                    logWith LogInfo ("FinalAssignment set to " <> head optionList)
                                                                                                                                    pure (globalSettings{finalAssignment = fromJust localMethod}, processedData, inSeedList)
                                                                                                                                else
                                                                                                                                    if localMethod == Just DirectOptimization
                                                                                                                                        then do
                                                                                                                                            pure (globalSettings{finalAssignment = fromJust localMethod}, processedData, inSeedList)
                                                                                                                                        else do
                                                                                                                                            logWith LogInfo "FinalAssignment set to DO (ignoring IA option) for non-Tree graphs"
                                                                                                                                            pure (globalSettings{finalAssignment = DirectOptimization}, processedData, inSeedList)
                                                                                                            else
                                                                                                                if head commandList == "graphfactor"
                                                                                                                    then
                                                                                                                        let localMethod
                                                                                                                                | (head optionList == "nopenalty") = Just NoNetworkPenalty
                                                                                                                                | (head optionList == "w15") = Just Wheeler2015Network
                                                                                                                                | (head optionList == "w23") = Just Wheeler2023Network
                                                                                                                                | (head optionList == "pmdl") = Just PMDLGraph
                                                                                                                                | otherwise = Nothing
                                                                                                                        in  if isNothing localMethod
                                                                                                                                then do
                                                                                                                                    failWithPhase
                                                                                                                                        Parsing
                                                                                                                                        ("Error in 'set' command. GraphFactor  '" <> head optionList <> "' is not 'NoPenalty', 'W15', 'W23', or 'PMDL'")
                                                                                                                                else do
                                                                                                                                    logWith LogInfo ("GraphFactor set to " <> show localMethod)
                                                                                                                                    pure (globalSettings{graphFactor = fromJust localMethod}, processedData, inSeedList)
                                                                                                                    else
                                                                                                                        if head commandList == "graphssteepest"
                                                                                                                            then
                                                                                                                                let localValue = readMaybe (head optionList) ∷ Maybe Int
                                                                                                                                in  if isNothing localValue
                                                                                                                                        then do
                                                                                                                                            failWithPhase
                                                                                                                                                Parsing
                                                                                                                                                ("Set option 'graphsSteepest' must be set to an integer value (e.g. graphsSteepest:5): " <> head optionList)
                                                                                                                                        else do
                                                                                                                                            logWith LogInfo ("GraphsStreepest set to " <> head optionList)
                                                                                                                                            pure (globalSettings{graphsSteepest = fromJust localValue}, processedData, inSeedList)
                                                                                                                            else
                                                                                                                                if head commandList == "graphtype"
                                                                                                                                    then
                                                                                                                                        let localGraphType
                                                                                                                                                | (head optionList == "tree") = Just Tree
                                                                                                                                                | (head optionList == "softwired") = Just SoftWired
                                                                                                                                                | (head optionList == "hardwired") = Just HardWired
                                                                                                                                                | otherwise = Nothing
                                                                                                                                        in  if isNothing localGraphType
                                                                                                                                                then do
                                                                                                                                                    failWithPhase
                                                                                                                                                        Parsing
                                                                                                                                                        ("Error in 'set' command. Graphtype '" <> head optionList <> "' is not 'tree', 'hardwired', or 'softwired'")
                                                                                                                                                else
                                                                                                                                                    if localGraphType /= Just Tree
                                                                                                                                                        then do
                                                                                                                                                            let netPenalty =
                                                                                                                                                                    if localGraphType == Just HardWired
                                                                                                                                                                        then NoNetworkPenalty
                                                                                                                                                                        else graphFactor globalSettings

                                                                                                                                                            logWith LogInfo ("Graphtype set to " <> head optionList <> " with graph factor NoPenalty and final assignment to DO")
                                                                                                                                                            pure
                                                                                                                                                                ( globalSettings{graphType = fromJust localGraphType, finalAssignment = DirectOptimization, graphFactor = netPenalty}
                                                                                                                                                                , processedData
                                                                                                                                                                , inSeedList
                                                                                                                                                                )
                                                                                                                                                        else do
                                                                                                                                                            logWith LogInfo ("Graphtype set to " <> head optionList)
                                                                                                                                                            pure (globalSettings{graphType = fromJust localGraphType}, processedData, inSeedList)
                                                                                                                                    else -- In first to do stuff above also

                                                                                                                                        if head commandList == "missingthreshold"
                                                                                                                                            then
                                                                                                                                                let localValue = readMaybe (head optionList) ∷ Maybe Int
                                                                                                                                                in  if isNothing localValue
                                                                                                                                                        then error ("Set option 'missingThreshold' must be set to an integer value (e.g. missingThreshold:50): " <> head optionList)
                                                                                                                                                        else
                                                                                                                                                            if fromJust localValue == missingThreshold globalSettings
                                                                                                                                                                then pure (globalSettings, processedData, inSeedList)
                                                                                                                                                                else do
                                                                                                                                                                    logWith LogInfo ("MissingThreshold set to " <> head optionList)
                                                                                                                                                                    pure (globalSettings{missingThreshold = fromJust localValue}, processedData, inSeedList)
                                                                                                                                            else
                                                                                                                                                if head commandList == "modelcomplexity"
                                                                                                                                                    then
                                                                                                                                                        let localValue = readMaybe (head optionList) ∷ Maybe Double
                                                                                                                                                        in  if isNothing localValue
                                                                                                                                                                then error ("Set option 'modelComplexity' must be set to a double value (e.g. modelComplexity:123.456): " <> head optionList)
                                                                                                                                                                else do
                                                                                                                                                                    logWith LogInfo ("Model Complexity set to " <> head optionList)
                                                                                                                                                                    pure (globalSettings{modelComplexity = fromJust localValue}, processedData, inSeedList)
                                                                                                                                                    else -- modify the behavior of rerooting character trees for all graph types

                                                                                                                                                        if head commandList == "multitraverse"
                                                                                                                                                            then do
                                                                                                                                                                let localCriterion
                                                                                                                                                                        | (head optionList == "true") = Just True
                                                                                                                                                                        | (head optionList == "false") = Just False
                                                                                                                                                                        | otherwise = Nothing

                                                                                                                                                                if isNothing localCriterion
                                                                                                                                                                    then do failWithPhase Parsing ("Error in 'set' command. MultiTraverse '" <> head optionList <> "' is not 'true' or 'false'")
                                                                                                                                                                    else do
                                                                                                                                                                        logWith LogInfo ("MultiTraverse set to " <> head optionList)
                                                                                                                                                                        pure (globalSettings{multiTraverseCharacters = fromJust localCriterion}, processedData, inSeedList)
                                                                                                                                                            else
                                                                                                                                                                if head commandList == "outgroup"
                                                                                                                                                                    then
                                                                                                                                                                        let outTaxonName = T.pack $ filter (/= '"') $ head $ filter (/= "") $ fmap snd argList
                                                                                                                                                                            outTaxonIndex = V.elemIndex outTaxonName leafNameVect
                                                                                                                                                                        in  if isNothing outTaxonIndex
                                                                                                                                                                                then do
                                                                                                                                                                                    failWithPhase
                                                                                                                                                                                        Parsing
                                                                                                                                                                                        ( "Error in 'set' command. Out-taxon "
                                                                                                                                                                                            <> T.unpack outTaxonName
                                                                                                                                                                                            <> " not found in input leaf list"
                                                                                                                                                                                            <> show (fmap T.unpack leafNameVect)
                                                                                                                                                                                        )
                                                                                                                                                                                else do
                                                                                                                                                                                    logWith LogInfo ("Outgroup set to " <> T.unpack outTaxonName)
                                                                                                                                                                                    pure (globalSettings{outgroupIndex = fromJust outTaxonIndex, outGroupName = outTaxonName}, processedData, inSeedList)
                                                                                                                                                                    else
                                                                                                                                                                        if head commandList == "partitioncharacter"
                                                                                                                                                                            then
                                                                                                                                                                                let localPartitionChar = head optionList
                                                                                                                                                                                in  if length localPartitionChar /= 1
                                                                                                                                                                                        then
                                                                                                                                                                                            errorWithoutStackTrace
                                                                                                                                                                                                ("Error in 'set' command. Partitioncharacter '" <> show localPartitionChar <> "' must be a single character")
                                                                                                                                                                                        else
                                                                                                                                                                                            if localPartitionChar /= partitionCharacter globalSettings
                                                                                                                                                                                                then do
                                                                                                                                                                                                    logWith LogInfo ("PartitionCharacter set to '" <> head optionList <> "'")
                                                                                                                                                                                                    pure (globalSettings{partitionCharacter = localPartitionChar}, processedData, inSeedList)
                                                                                                                                                                                                else pure (globalSettings, processedData, inSeedList)
                                                                                                                                                                            else
                                                                                                                                                                                if head commandList == "reportnaivedata"
                                                                                                                                                                                    then do
                                                                                                                                                                                        let localMethod
                                                                                                                                                                                                | (head optionList == "true") = True
                                                                                                                                                                                                | (head optionList == "false") = False
                                                                                                                                                                                                | otherwise =
                                                                                                                                                                                                    errorWithoutStackTrace ("Error in 'set' command. NeportNaive  '" <> head optionList <> "' is not 'True' or 'False'")

                                                                                                                                                                                        logWith LogInfo ("ReportNaiveData set to " <> show localMethod)
                                                                                                                                                                                        pure (globalSettings{reportNaiveData = localMethod}, processedData, inSeedList)
                                                                                                                                                                                    else
                                                                                                                                                                                        if head commandList == "rootcost"
                                                                                                                                                                                            then do
                                                                                                                                                                                                let localMethod
                                                                                                                                                                                                        | (head optionList == "norootcost") = NoRootCost
                                                                                                                                                                                                        | (head optionList == "w15") = Wheeler2015Root
                                                                                                                                                                                                        | (head optionList == "pmdl") = PMDLRoot
                                                                                                                                                                                                        | (head optionList == "ml") = MLRoot
                                                                                                                                                                                                        | otherwise =
                                                                                                                                                                                                            errorWithoutStackTrace ("Error in 'set' command. RootCost  '" <> head optionList <> "' is not 'NoRootCost', 'W15', or 'PMDL'")

                                                                                                                                                                                                let lRootComplexity
                                                                                                                                                                                                        | localMethod == NoRootCost = 0.0
                                                                                                                                                                                                        | localMethod `elem` [Wheeler2015Root, PMDLRoot, MLRoot] = U.calculateW15RootCost processedData
                                                                                                                                                                                                        | otherwise = error ("Root cost method not recognized: " <> show localMethod)

                                                                                                                                                                                                logWith LogInfo ("RootCost set to " <> show localMethod <> " " <> show lRootComplexity <> " bits")
                                                                                                                                                                                                pure (globalSettings{rootCost = localMethod, rootComplexity = lRootComplexity}, processedData, inSeedList)
                                                                                                                                                                                            else
                                                                                                                                                                                                if head commandList == "seed"
                                                                                                                                                                                                    then
                                                                                                                                                                                                        let localValue = readMaybe (head optionList) ∷ Maybe Int
                                                                                                                                                                                                        in  if isNothing localValue
                                                                                                                                                                                                                then error ("Set option 'seed' must be set to an integer value (e.g. seed:123): " <> head optionList)
                                                                                                                                                                                                                else do
                                                                                                                                                                                                                    logWith LogInfo ("Random Seed set to " <> head optionList)
                                                                                                                                                                                                                    pure (globalSettings{seed = fromJust localValue}, processedData, randomIntList (fromJust localValue))
                                                                                                                                                                                                    else
                                                                                                                                                                                                        if head commandList == "softwiredmethod"
                                                                                                                                                                                                            then do
                                                                                                                                                                                                                let localMethod
                                                                                                                                                                                                                        | (head optionList == "naive") = Naive
                                                                                                                                                                                                                        | (head optionList == "exhaustive") = Naive
                                                                                                                                                                                                                        | (head optionList == "resolutioncache") = ResolutionCache
                                                                                                                                                                                                                        | otherwise =
                                                                                                                                                                                                                            errorWithoutStackTrace
                                                                                                                                                                                                                                ("Error in 'set' command. SoftwiredMethod  '" <> head optionList <> "' is not 'Exhaustive' or 'ResolutionCache'")

                                                                                                                                                                                                                logWith LogInfo ("SoftwiredMethod " <> show localMethod)
                                                                                                                                                                                                                pure (globalSettings{softWiredMethod = localMethod}, processedData, inSeedList)
                                                                                                                                                                                                            else -- modify the use of Network Add heurisitcs in network optimization

                                                                                                                                                                                                                if head commandList == "usenetaddheuristic"
                                                                                                                                                                                                                    then do
                                                                                                                                                                                                                        let localCriterion
                                                                                                                                                                                                                                | (head optionList == "true") = True
                                                                                                                                                                                                                                | (head optionList == "false") = False
                                                                                                                                                                                                                                | otherwise =
                                                                                                                                                                                                                                    errorWithoutStackTrace ("Error in 'set' command. UseNetAddHeuristic '" <> head optionList <> "' is not 'true' or 'false'")

                                                                                                                                                                                                                        logWith LogInfo ("UseNetAddHeuristic set to " <> head optionList)
                                                                                                                                                                                                                        pure (globalSettings{useNetAddHeuristic = localCriterion}, processedData, inSeedList)
                                                                                                                                                                                                                    else -- these not intended for users

                                                                                                                                                                                                                        if head commandList == "jointhreshold"
                                                                                                                                                                                                                            then
                                                                                                                                                                                                                                let localValue = readMaybe (head optionList) ∷ Maybe Double
                                                                                                                                                                                                                                in  if isNothing localValue
                                                                                                                                                                                                                                        then error ("Set option 'joinThreshold' must be set to an double value >= 1.0 (e.g. joinThreshold:1.17): " <> head optionList)
                                                                                                                                                                                                                                        else
                                                                                                                                                                                                                                            if fromJust localValue < 1.0
                                                                                                                                                                                                                                                then
                                                                                                                                                                                                                                                    errorWithoutStackTrace
                                                                                                                                                                                                                                                        ("Set option 'joinThreshold' must be set to a double value >= 1.0 (e.g. joinThreshold:1.17): " <> head optionList)
                                                                                                                                                                                                                                                else do
                                                                                                                                                                                                                                                    logWith LogInfo ("JoinThreshold set to " <> head optionList)
                                                                                                                                                                                                                                                    pure (globalSettings{unionThreshold = fromJust localValue}, processedData, inSeedList)
                                                                                                                                                                                                                            else -- parallel strategy settings options

                                                                                                                                                                                                                                if head commandList == "defparstrat"
                                                                                                                                                                                                                                    then do
                                                                                                                                                                                                                                        let localMethod
                                                                                                                                                                                                                                                | (head optionList == "r0") = R0
                                                                                                                                                                                                                                                | (head optionList == "rpar") = RPar
                                                                                                                                                                                                                                                | (head optionList == "rseq") = RSeq
                                                                                                                                                                                                                                                | (head optionList == "rdeepseq") = RDeepSeq
                                                                                                                                                                                                                                                | otherwise =
                                                                                                                                                                                                                                                    errorWithoutStackTrace
                                                                                                                                                                                                                                                        ("Error in 'set' command. DefParStrat  '" <> head optionList <> "' is not 'r0', 'WrPar', 'rSeq', or 'rDeepSeq'")

                                                                                                                                                                                                                                        logWith LogInfo ("DefParStrat set to " <> show localMethod)
                                                                                                                                                                                                                                        pure (globalSettings{defaultParStrat = localMethod}, processedData, inSeedList)
                                                                                                                                                                                                                                    else
                                                                                                                                                                                                                                        if head commandList == "lazyparstrat"
                                                                                                                                                                                                                                            then do
                                                                                                                                                                                                                                                let localMethod
                                                                                                                                                                                                                                                        | (head optionList == "r0") = R0
                                                                                                                                                                                                                                                        | (head optionList == "rpar") = RPar
                                                                                                                                                                                                                                                        | (head optionList == "rseq") = RSeq
                                                                                                                                                                                                                                                        | (head optionList == "rdeepseq") = RDeepSeq
                                                                                                                                                                                                                                                        | otherwise =
                                                                                                                                                                                                                                                            errorWithoutStackTrace
                                                                                                                                                                                                                                                                ("Error in 'set' command. LazyParStrat  '" <> head optionList <> "' is not 'r0', 'WrPar', 'rSeq', or 'rDeepSeq'")

                                                                                                                                                                                                                                                logWith LogInfo ("LazyParStrat set to " <> show localMethod)
                                                                                                                                                                                                                                                pure (globalSettings{lazyParStrat = localMethod}, processedData, inSeedList)
                                                                                                                                                                                                                                            else
                                                                                                                                                                                                                                                if head commandList == "strictparstrat"
                                                                                                                                                                                                                                                    then do
                                                                                                                                                                                                                                                        let localMethod
                                                                                                                                                                                                                                                                | (head optionList == "r0") = R0
                                                                                                                                                                                                                                                                | (head optionList == "rpar") = RPar
                                                                                                                                                                                                                                                                | (head optionList == "rseq") = RSeq
                                                                                                                                                                                                                                                                | (head optionList == "rdeepseq") = RDeepSeq
                                                                                                                                                                                                                                                                | otherwise =
                                                                                                                                                                                                                                                                    errorWithoutStackTrace
                                                                                                                                                                                                                                                                        ("Error in 'set' command. StrictParStrat  '" <> head optionList <> "' is not 'r0', 'WrPar', 'rSeq', or 'rDeepSeq'")

                                                                                                                                                                                                                                                        logWith LogInfo ("StrictParStrat set to " <> show localMethod)
                                                                                                                                                                                                                                                        pure (globalSettings{strictParStrat = localMethod}, processedData, inSeedList)
                                                                                                                                                                                                                                                    else -- modify the use of implied alkignemnt in heuristics

                                                                                                                                                                                                                                                        if head commandList == "useia"
                                                                                                                                                                                                                                                            then do
                                                                                                                                                                                                                                                                let localCriterion
                                                                                                                                                                                                                                                                        | (head optionList == "true") = True
                                                                                                                                                                                                                                                                        | (head optionList == "false") = False
                                                                                                                                                                                                                                                                        | otherwise = errorWithoutStackTrace ("Error in 'set' command. UseIA '" <> head optionList <> "' is not 'true' or 'false'")

                                                                                                                                                                                                                                                                logWith LogInfo ("UseIA set to " <> head optionList)
                                                                                                                                                                                                                                                                pure (globalSettings{useIA = localCriterion}, processedData, inSeedList)
                                                                                                                                                                                                                                                            else do
                                                                                                                                                                                                                                                                logWith LogInfo ("Warning: Unrecognized/missing 'set' option in " <> show argList)
                                                                                                                                                                                                                                                                pure (globalSettings, processedData, inSeedList)


{- | reportCommand takes report options, current data and graphs and returns a
(potentially large) String to print and the channel to print it to
and write mode overwrite/append
if global settings reportNaiveData is True then need to rediagnose graph with processed data since
naiveData was sent to command and will not match what is in the optimized graphs
-}
reportCommand
    ∷ GlobalSettings
    → [Argument]
    → ([NameText], [(NameText, NameText)])
    → Int
    → String
    → ProcessedData
    → [ReducedPhylogeneticGraph]
    → [ReducedPhylogeneticGraph]
    → [[VertexCost]]
    → PhyG (String, String, String)
reportCommand globalSettings argList excludeRename numInputFiles crossReferenceString processedData curGraphs supportGraphs pairwiseDistanceMatrix =
    let argListWithoutReconcileCommands = filter ((`notElem` VER.reconcileArgList) . fst) argList
        -- check for balances double quotes and only one pair
        outFileNameList = filter (/= "") $ fmap snd argListWithoutReconcileCommands -- argList
        commandList = fmap (fmap C.toLower) $ filter (/= "") $ fmap fst argListWithoutReconcileCommands
    in  -- reconcileList = filter (/= "") $ fmap fst argList

        if length outFileNameList > 1
            then do failWithPhase Outputting ("Report can only have one file name: " <> show outFileNameList <> " " <> show argList)
            else
                let checkCommandList = checkCommandArgs "report" commandList VER.reportArgList
                    outfileName =
                        if null outFileNameList
                            then "stderr"
                            else tail $ L.init $ head outFileNameList
                    writeMode =
                        if "overwrite" `elem` commandList
                            then "overwrite"
                            else "append"
                in  -- error too harsh, lose everything else
                    -- if (null $ filter (/= "overwrite") $ filter (/= "append") commandList) then errorWithoutStackTrace ("Error: Missing 'report' option in " <> show commandList)
                    -- else
                    if not checkCommandList
                        then errorWithoutStackTrace ("Unrecognized command in report: " <> show argList)
                        else -- This for reconciled data

                            if "crossrefs" `elem` commandList
                                then
                                    let dataString = crossReferenceString
                                    in  pure (dataString, outfileName, writeMode)
                                else
                                    if "data" `elem` commandList
                                        then
                                            let blocks = thd3 processedData
                                                numChars = V.sum $ fmap (V.length . thd3) blocks
                                                dataString = phyloDataToString 0 $ thd3 processedData
                                                baseData =
                                                    [ ["Input data contained:"]
                                                    , ["", show (length $ fst3 processedData) <> " terminal taxa"]
                                                    , ["", show numInputFiles <> " input data files"]
                                                    , ["", show (length blocks) <> " character blocks"]
                                                    , ["", show numChars <> " total characters"]
                                                    ]
                                                leafNames = V.toList (T.unpack <$> fst3 processedData)
                                                leafField = ["Terminal taxa:"] : fmap (" " :) (SL.chunksOf 10 leafNames)
                                                excludedTaxa =
                                                    if (not . null . fst) excludeRename
                                                        then T.unpack <$> fst excludeRename
                                                        else ["None"]
                                                excludedField = ["Excluded taxa:"] : fmap (" " :) (SL.chunksOf 10 excludedTaxa)
                                                renameFirstList = fmap (((: []) . T.unpack) . fst) (snd excludeRename)
                                                renameSecondList = fmap (((: []) . T.unpack) . snd) (snd excludeRename)
                                                renamePairList =
                                                    if (not . null . snd) excludeRename
                                                        then (" " :) <$> zipWith (<>) renameFirstList renameSecondList
                                                        else [[" ", "None", "None"]]
                                                renameField = ["Renamed taxa:", "New Name", "Original Name"] : renamePairList
                                                charInfoFields = ["Index", "Block", "Name", "Type", "Activity", "Weight", "Prealigned", "Alphabet", "TCM"]
                                            in  pure
                                                    ( CSV.genCsvFile
                                                        (baseData <> [[""]] <> leafField <> [[""]] <> excludedField <> [[""]] <> renameField <> [[""]] <> (charInfoFields : dataString))
                                                    , outfileName
                                                    , writeMode
                                                    )
                                        else
                                            if "diagnosis" `elem` commandList
                                                then -- need to rediagnose if reportNaiveData

                                                    let curGraphs' =
                                                            if not (reportNaiveData globalSettings)
                                                                then curGraphs
                                                                else
                                                                    PU.seqParMap
                                                                        (parStrategy $ strictParStrat globalSettings)
                                                                        (TRAV.multiTraverseFullyLabelGraphReduced globalSettings processedData False False Nothing)
                                                                        (fmap fst5 curGraphs)
                                                        dataString = CSV.genCsvFile $ concatMap (getGraphDiagnosis globalSettings processedData) (zip curGraphs' [0 .. (length curGraphs' - 1)])
                                                    in  if null curGraphs
                                                            then do
                                                                logWith LogInfo "No graphs to diagnose"
                                                                pure ("No graphs to diagnose", outfileName, writeMode)
                                                            else do
                                                                logWith
                                                                    LogInfo
                                                                    ("Diagnosing " <> show (length curGraphs) <> " graphs at minimum cost " <> show (minimum $ fmap snd5 curGraphs))
                                                                pure (dataString, outfileName, writeMode)
                                                else
                                                    if "displaytrees" `elem` commandList
                                                        then -- need to specify -O option for multiple graphs

                                                            let rediagnodesGraphs = fmap (TRAV.multiTraverseFullyLabelGraph globalSettings processedData False False Nothing) (fmap fst5 curGraphs)
                                                                inputDisplayVVList = fmap fth6 rediagnodesGraphs
                                                                costList = fmap snd5 curGraphs
                                                                displayCostListList = fmap GO.getDisplayTreeCostList rediagnodesGraphs
                                                                displayInfoString =
                                                                    if ("dot" `elem` commandList) || ("dotpdf" `elem` commandList)
                                                                        then ("//DisplayTree costs : " <> show (fmap (sum . fst) displayCostListList, displayCostListList))
                                                                        else -- newick

                                                                            let middle = fmap bracketToCurly $ show (fmap (sum . fst) displayCostListList, displayCostListList)
                                                                            in  ("[DisplayTree costs : " <> middle <> "]")

                                                                treeIndexStringList =
                                                                    if ("dot" `elem` commandList) || ("dotpdf" `elem` commandList)
                                                                        then fmap (((<> "\n") . ("//Canonical Tree " <>)) . show) [0 .. (length inputDisplayVVList - 1)]
                                                                        else -- newick
                                                                            fmap (((<> "]\n") . ("[Canonical Tree " <>)) . show) [0 .. (length inputDisplayVVList - 1)]
                                                                canonicalGraphPairList = zip treeIndexStringList inputDisplayVVList
                                                                blockStringList = unlines (fmap (outputBlockTrees commandList costList (outgroupIndex globalSettings)) canonicalGraphPairList)
                                                            in  -- graphString = outputGraphString commandList (outgroupIndex globalSettings) (fmap thd6 curGraphs) (fmap snd6 curGraphs)

                                                                if null curGraphs || graphType globalSettings /= SoftWired
                                                                    then do
                                                                        logWith LogInfo "No soft-wired graphs to report display trees"
                                                                        pure ("No soft-wired graphs to report display trees", outfileName, writeMode)
                                                                    else pure (displayInfoString <> "\n" <> blockStringList, outfileName, writeMode)
                                                        else
                                                            if "graphs" `elem` commandList
                                                                then -- else if (not .null) (L.intersect ["graphs", "newick", "dot", "dotpdf"] commandList) then

                                                                    let graphString = outputGraphString commandList (outgroupIndex globalSettings) (fmap thd5 curGraphs) (fmap snd5 curGraphs)
                                                                    in  if null curGraphs
                                                                            then do
                                                                                logWith LogInfo "No graphs to report"
                                                                                pure ("No graphs to report", outfileName, writeMode)
                                                                            else do
                                                                                logWith
                                                                                    LogInfo
                                                                                    ("Reporting " <> show (length curGraphs) <> " graph(s) at minimum cost " <> show (minimum $ fmap snd5 curGraphs))
                                                                                pure (graphString, outfileName, writeMode)
                                                                else
                                                                    if "ia" `elem` commandList || "impliedalignment" `elem` commandList
                                                                        then
                                                                            if null curGraphs
                                                                                then do
                                                                                    logWith LogInfo "No graphs to create implied alignments"
                                                                                    pure ("No impliedAlgnments to report", outfileName, writeMode)
                                                                                else
                                                                                    let includeMissing = elem "includemissing" commandList
                                                                                        concatSeqs = elem "concatenate" commandList
                                                                                        
                                                                                    in  do
                                                                                        iaContentList <- mapM (getImpliedAlignmentString globalSettings (includeMissing || concatSeqs) concatSeqs processedData) (zip curGraphs [0 .. (length curGraphs - 1)])
                                                                                        logWith
                                                                                                LogInfo
                                                                                                "\tWarning: Prealigned sequence data with non-additive type costs (all change values equal) have been recoded to non-additive characters and will not appear in implied alignment output."
                                                                                        pure (concat iaContentList, outfileName, writeMode)
                                                                        else
                                                                            if "pairdist" `elem` commandList
                                                                                then
                                                                                    let nameData = L.intercalate "," (V.toList (T.unpack <$> fst3 processedData)) <> "\n"
                                                                                        dataString = CSV.genCsvFile $ fmap (fmap show) pairwiseDistanceMatrix
                                                                                    in  pure (nameData <> dataString, outfileName, writeMode)
                                                                                else
                                                                                    if "reconcile" `elem` commandList
                                                                                        then
                                                                                            let (reconcileString, _) = R.makeReconcileGraph VER.reconcileArgList argList (fmap fst5 curGraphs)
                                                                                            in  if null curGraphs
                                                                                                    then do
                                                                                                        logWith LogInfo "No graphs to reconcile"
                                                                                                        pure ([], outfileName, writeMode)
                                                                                                    else pure (reconcileString, outfileName, writeMode)
                                                                                        else
                                                                                            if "search" `elem` commandList
                                                                                                then
                                                                                                    let dataString' = fmap showSearchFields $ reverse $ searchData globalSettings
                                                                                                        -- reformat the "search" command fields a bit
                                                                                                        dataString = processSearchFields dataString'
                                                                                                        sysInfoData =
                                                                                                            "System Info, OS: "
                                                                                                                <> SI.os
                                                                                                                <> ", Chip Arch: "
                                                                                                                <> SI.arch
                                                                                                                <> ", Compiler: "
                                                                                                                <> SI.compilerName
                                                                                                                <> " "
                                                                                                                <> DV.showVersion SI.compilerVersion
                                                                                                                <> ", Compile Date: "
                                                                                                                <> (__DATE__ <> " " <> __TIME__)
                                                                                                        cpuInfoString =
                                                                                                            if SI.os /= "linux"
                                                                                                                then "CPU Info, No /proc/cpuinfo on darwin"
                                                                                                                else
                                                                                                                    let cpuInfoM = SIOU.unsafePerformIO SC.tryGetCPUs
                                                                                                                    in  if isNothing cpuInfoM
                                                                                                                            then "CPU Info, Couldn't parse CPU Info"
                                                                                                                            else
                                                                                                                                "CPU Info, Physical Processors: "
                                                                                                                                    <> show (SC.physicalProcessors (fromJust cpuInfoM))
                                                                                                                                    <> ", Physical Cores: "
                                                                                                                                    <> show (SC.physicalCores (fromJust cpuInfoM))
                                                                                                                                    <> ", Logical Cores: "
                                                                                                                                    <> show (SC.logicalCores (fromJust cpuInfoM))
                                                                                                        baseData = sysInfoData <> "\n" <> cpuInfoString <> "\nSearchData\nRandom seed, " <> show (seed globalSettings) <> "\n"
                                                                                                        charInfoFields =
                                                                                                            [ "Command"
                                                                                                            , "Arguments"
                                                                                                            , "Min cost in"
                                                                                                            , "Max cost in"
                                                                                                            , "Num graphs in"
                                                                                                            , "Min cost out"
                                                                                                            , "Max cost out"
                                                                                                            , "Num graphs out"
                                                                                                            , "CPU time (secs)"
                                                                                                            , "Comment"
                                                                                                            ]
                                                                                                    in  pure (baseData <> CSV.genCsvFile (charInfoFields : dataString), outfileName, writeMode)
                                                                                                else
                                                                                                    if "support" `elem` commandList
                                                                                                        then
                                                                                                            let graphString = outputGraphStringSimple commandList (outgroupIndex globalSettings) (fmap fst5 supportGraphs) (fmap snd5 supportGraphs)
                                                                                                            in  -- trace ("Rep Sup: " <> (LG.prettify $ fst5 $ head supportGraphs)) (
                                                                                                                if null supportGraphs
                                                                                                                    then do
                                                                                                                        logWith LogInfo "\tNo support graphs to report"
                                                                                                                        pure ([], outfileName, writeMode)
                                                                                                                    else do
                                                                                                                        logWith LogInfo ("Reporting " <> show (length curGraphs) <> " support graph(s)")
                                                                                                                        pure (graphString, outfileName, writeMode)
                                                                                                        else -- )

                                                                                                            if "tnt" `elem` commandList
                                                                                                                then
                                                                                                                    if null curGraphs
                                                                                                                        then do
                                                                                                                            logWith LogInfo "No graphs to create implied alignments for TNT output"
                                                                                                                            pure ("No impliedAlgnments for TNT to report", outfileName, writeMode)
                                                                                                                        else
                                                                                                                            let curGraphs' =
                                                                                                                                    if not (reportNaiveData globalSettings)
                                                                                                                                        then (fmap GO.convertReduced2PhylogeneticGraph curGraphs)
                                                                                                                                        else
                                                                                                                                            PU.seqParMap
                                                                                                                                                (parStrategy $ strictParStrat globalSettings)
                                                                                                                                                (TRAV.multiTraverseFullyLabelGraph globalSettings processedData False False Nothing)
                                                                                                                                                (fmap fst5 curGraphs)
                                                                                                                            in do
                                                                                                                            tntContentList' <- mapM (getTNTString globalSettings processedData) (zip curGraphs' [0 .. (length curGraphs' - 1)])
                                                                                                                            let tntContentList = concat tntContentList'
                                                                                                                            pure (tntContentList, outfileName, writeMode)
                                                                                                                else do
                                                                                                                    logWith LogWarn ("\nUnrecognized/missing report option in " <> show commandList <> " defaulting to 'graphs'")
                                                                                                                    let graphString = outputGraphString commandList (outgroupIndex globalSettings) (fmap thd5 curGraphs) (fmap snd5 curGraphs)
                                                                                                                    if null curGraphs
                                                                                                                        then do
                                                                                                                            logWith LogInfo "No graphs to report"
                                                                                                                            pure ("No graphs to report", outfileName, writeMode)
                                                                                                                        else do
                                                                                                                            logWith
                                                                                                                                LogInfo
                                                                                                                                ("Reporting " <> show (length curGraphs) <> " graph(s) at minimum cost " <> show (minimum $ fmap snd5 curGraphs) <> "\n")
                                                                                                                            pure (graphString, outfileName, writeMode)
    where
        bracketToCurly a =
            if a == '('
                then '{'
                else
                    if a == ')'
                        then '}'
                        else a
