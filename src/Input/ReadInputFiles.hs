{- |
Module exposing functionality for reading input files for phylogenetic analysis.
-}
module Input.ReadInputFiles (
    executeReadCommands,
    getReadArgs,
    extractInputTuple,
    expandReadCommands,
) where

import Commands.Verify qualified as V
import Control.Monad (when)
import Control.Monad.IO.Class (MonadIO (..))
import Control.Monad.Random.Class
import Data.Char
import Data.Char qualified as C
import Data.Foldable
import Data.Graph.Inductive.Basic qualified as B
import Data.List qualified as L
import Data.Maybe
import Data.Text.Lazy qualified as T
import Data.Text.Lazy.IO qualified as TIO
import Data.Text.Short qualified as ST
import Debug.Trace
import GeneralUtilities qualified as GU
import GraphFormatUtilities qualified as GFU
import Input.FastAC qualified as FAC
import Input.TNTUtilities qualified as TNT
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import PHANE.Evaluation.Verbosity (Verbosity (..))
import System.Directory
import System.IO
import System.Path.Glob qualified as SPG
import Text.Read
import Types.Types
import Utilities.LocalGraph qualified as LG
import Utilities.Utilities qualified as U


{- | expandReadCommands expands read commands to multiple satisfying wild cards
read command can have multiple file names
-}
expandReadCommands ∷ [Command] → Command → PhyG [Command]
expandReadCommands _newReadList inCommand@(commandType, argList') =
    let argList = filter ((`notElem` ["tcm"]) . fst) argList'
        tcmArgList = filter ((`elem` ["tcm"]) . fst) argList'
        fileNames = fmap snd $ filter ((/= "tcm") . fst) $ filter ((/= "") . snd) argList'
        modifierList = fmap fst argList
    in  case commandType of
          Read -> do
              globbedFileNames ← liftIO $ mapM SPG.glob fileNames
              when (all null globbedFileNames) . failWithPhase Inputting $ unwords
                  [ "File(s) not found in 'read' command (could be due to incorrect filename or missing closing double quote '\"''):"
                  , show fileNames
                  ]

              newArgPairs ← mapM makeNewArgs (zip modifierList globbedFileNames)

              let commandList = replicate (length newArgPairs) commandType
              let tcmArgListL = replicate (length newArgPairs) tcmArgList

              pure $ zip commandList (zipWith (<>) newArgPairs tcmArgListL)

          _ -> failWithPhase Inputting $ unwords [ "Incorrect command type in expandReadCommands:", show inCommand ]


{- | makeNewArgs takes an argument modifier (first in pair) and replicates is and zips with
globbed file name list to create a list of arguments
-}
makeNewArgs ∷ (String, [String]) → PhyG [(String, String)]
makeNewArgs (modifier, fileNameList) =
    if null fileNameList
        then do failWithPhase Inputting ("Null filename list in makeNewArgs: " <> show (modifier, fileNameList))
        else
            let modList = replicate (length fileNameList) modifier
            in  do return $ zip modList fileNameList


{-
Added strictness on read so could close files after reading to
   allow for large number (1000's) of files to be input.
   this may not be the way to go, especially if files are big and
   closing them is not an issue
-}

{- | extractInputTuple  takes the list of pairs from mapped executeReadCommands
and returns ([RawData], [SimpleGraph])
-}
extractInputTuple
    ∷ [([RawData], [SimpleGraph], [NameText], [NameText], [(NameText, NameText)], [(NameText, NameText)])]
    → ([RawData], [SimpleGraph], [NameText], [NameText], [(NameText, NameText)], [(NameText, NameText)])
extractInputTuple dataGraphList =
    let (inDataList, inGraphList, inTerminalsList, inExcludeList, inRenamePairs, inReBlockPairs) = L.unzip6 dataGraphList
        rawData = L.sort $ concat inDataList
        rawGraphs = concat inGraphList
        rawTerminals = concat inTerminalsList
        excludeTerminals = concat inExcludeList
        renamePairs = concat inRenamePairs
        reBlockPairs = concat inReBlockPairs
    in  (rawData, rawGraphs, rawTerminals, excludeTerminals, renamePairs, reBlockPairs)


-- | executeReadCommands wrapper for executeReadCommands' with out all the empty list imput
executeReadCommands
    ∷ [Argument]
    → PhyG ([RawData], [SimpleGraph], [NameText], [NameText], [(NameText, NameText)], [(NameText, NameText)])
executeReadCommands = executeReadCommands' [] [] [] [] [] [] False ([], [], 1.0)


{- | executeReadCommands' reads input files and returns raw data, input graphs, and terminal taxa to include
assumes that "prealigned" and "tcm:File" are the first arguments if they are specified
so that they can apply to all the files in the command without order depence
-}
executeReadCommands'
    ∷ [RawData]
    → [SimpleGraph]
    → [NameText]
    → [NameText]
    → [(NameText, NameText)]
    → [(NameText, NameText)]
    → Bool
    → ([ST.ShortText], [[Int]], Double)
    → [Argument]
    → PhyG ([RawData], [SimpleGraph], [NameText], [NameText], [(NameText, NameText)], [(NameText, NameText)])
executeReadCommands' curData curGraphs curTerminals curExcludeList curRenamePairs curReBlockPairs _ tcmPair argList = do
    if null argList
        then return (curData, curGraphs, curTerminals, curExcludeList, curRenamePairs, curReBlockPairs)
        else do
            -- logWith LogInfo ("ERC: " <> (show argList) <> " " <> (show tcmPair))
            let isPrealigned' = False -- removed this option and m,asde specific to sequence types
            --      | isPrealigned = True
            --      | "prealigned" `elem`  fmap fst argList = True
            --      | otherwise = False
            let (firstOption', firstFile) = head argList
            let firstOption = fmap C.toLower firstOption'
            -- Check for prealigned
            -- if firstOption == "prealigned" then
            --    executeReadCommands' curData curGraphs curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)
            -- else do
            do
                fileHandle ←
                    liftIO $
                        if ',' `notElem` firstFile
                            then openFile firstFile ReadMode
                            else return (stdin ∷ Handle)

                canBeReadFrom ← liftIO $ hIsReadable fileHandle
                if not canBeReadFrom
                    then do failWithPhase Inputting ("\n\n'Read' error: file " <> firstFile <> " cannot be read")
                    else
                        if not $ null firstOption
                            then do logWith LogInfo ("Reading " <> firstFile <> " with option " <> firstOption <> "\n")
                            else do logWith LogInfo ("Reading " <> firstFile <> " with no options" <> "\n")

                -- this is awkward but need to use dot utilities
                -- added split, write, read hack for multiple input graphs in single file
                if firstOption == "dot"
                    then do
                        fileContents ← liftIO $ TIO.hGetContents fileHandle
                        let splitGraphsText = splitDotGraphs fileContents
                        -- logWith LogInfo $ (concat $ fmap T.unpack splitGraphsText)

                        inputDotList <- mapM (processDotFile firstFile) splitGraphsText

                        executeReadCommands'
                            curData
                            (inputDotList <> curGraphs)
                            curTerminals
                            curExcludeList
                            curRenamePairs
                            curReBlockPairs
                            isPrealigned'
                            tcmPair
                            (tail argList)

                    else -- not "dot" files
                    do
                        -- logWith LogInfo ("FC1: " <> firstFile)
                        fileContents ←
                            liftIO $
                                if ',' `notElem` firstFile
                                    then TIO.hGetContents fileHandle
                                    else return (T.pack firstFile)

                        {--Strict version--don't think necessary--but could getvinot open fikle number issues?
                         -- destroys lazyness but allows closing right away
                         -- this so don't have huge numbers of open files for large data sets
                        fileContents <- hGetContents' fileHandle
                        hClose fileHandle
                        -}

                        if T.null fileContents
                            then do failWithPhase Inputting ("Error: Input file " <> firstFile <> " is empty")
                            else -- try to figure out file type based on first and following characters

                                if firstOption == "tcm"
                                    then do
                                        newTCMPair ← processTCMContents (',' `elem` firstFile) (T.unpack fileContents) firstFile
                                        executeReadCommands'
                                            curData
                                            curGraphs
                                            curTerminals
                                            curExcludeList
                                            curRenamePairs
                                            curReBlockPairs
                                            isPrealigned'
                                            newTCMPair
                                            (tail argList)
                                    else
                                        if null firstOption
                                            then
                                                let firstChar = T.head $ T.dropWhile (== ' ') fileContents
                                                in  if (toLower firstChar == '/') || (toLower firstChar == 'd') || (toLower firstChar == 'g')
                                                        then do
                                                            logWith LogInfo ("\tTrying to parse " <> firstFile <> " as dot" <> "\n")
                                                            -- destroys lazyness but allows closing right away
                                                            -- this so don't have huge numbers of open files for large data sets
                                                            fileHandle2 ← liftIO $ openFile firstFile ReadMode
                                                            dotGraph ← liftIO $ LG.hGetDotLocal fileHandle2
                                                            liftIO $ hClose fileHandle2
                                                            let inputDot = GFU.relabelFGL $ LG.dotToGraph dotGraph
                                                            let hasLoops = B.hasLoop inputDot
                                                            if hasLoops
                                                                then do failWithPhase Parsing ("Input graph in " <> firstFile <> "  has loops/self-edges" <> "\n")
                                                                else do logWith LogInfo ""
                                                            let hasCycles = GFU.cyclic inputDot
                                                            if hasCycles
                                                                then do failWithPhase Parsing ("Input graph in " <> firstFile <> " has at least one cycle" <> "\n")
                                                                else
                                                                    executeReadCommands'
                                                                        curData
                                                                        (inputDot : curGraphs)
                                                                        curTerminals
                                                                        curExcludeList
                                                                        curRenamePairs
                                                                        curReBlockPairs
                                                                        isPrealigned'
                                                                        tcmPair
                                                                        (tail argList)
                                                        else
                                                            if (toLower firstChar == '<') || (toLower firstChar == '(')
                                                                then
                                                                    let thisGraphList = getFENewickGraphText fileContents
                                                                        hasCycles = filter id $ fmap GFU.cyclic thisGraphList
                                                                        hasLoops = filter id $ fmap B.hasLoop thisGraphList
                                                                    in  if not $ null hasLoops
                                                                            then do failWithPhase Parsing ("Input graph in " <> firstFile <> "  has loops/self-edges" <> "\n")
                                                                            else
                                                                                if not $ null hasCycles
                                                                                    then do failWithPhase Parsing ("Input graph in " <> firstFile <> " has at least one cycle" <> "\n")
                                                                                    else
                                                                                        executeReadCommands'
                                                                                            curData
                                                                                            (thisGraphList <> curGraphs)
                                                                                            curTerminals
                                                                                            curExcludeList
                                                                                            curRenamePairs
                                                                                            curReBlockPairs
                                                                                            isPrealigned'
                                                                                            tcmPair
                                                                                            (tail argList)
                                                                else
                                                                    if toLower firstChar == 'x'
                                                                        then do
                                                                            tntData ← TNT.getTNTDataText fileContents firstFile
                                                                            logWith LogInfo ("\tTrying to parse " <> firstFile <> " as TNT" <> "\n")
                                                                            executeReadCommands'
                                                                                (tntData : curData)
                                                                                curGraphs
                                                                                curTerminals
                                                                                curExcludeList
                                                                                curRenamePairs
                                                                                curReBlockPairs
                                                                                isPrealigned'
                                                                                tcmPair
                                                                                (tail argList)
                                                                        else
                                                                            let fileContents' = T.unlines $ filter (not . T.null) $ T.takeWhile (/= ';') <$> T.lines fileContents
                                                                            in  if T.null (T.dropWhile (== ' ') fileContents')
                                                                                    then do failWithPhase Parsing ("Null file '" <> firstFile <> "' input")
                                                                                    else
                                                                                        if T.head (T.dropWhile (== ' ') fileContents') == '>'
                                                                                            then
                                                                                                let secondLine = T.lines fileContents' !! 1
                                                                                                    hasSpaces = T.find (== ' ') secondLine
                                                                                                in  -- numWords = length $ words secondLine

                                                                                                    -- spaces between alphabet elements suggest fastc
                                                                                                    if isJust hasSpaces
                                                                                                        then
                                                                                                            let fastcData = FAC.getFastCText fileContents firstFile isPrealigned'
                                                                                                            in  do
                                                                                                                    fastcCharInfo ← FAC.getFastcCharInfo fastcData firstFile isPrealigned' tcmPair
                                                                                                                    logWith LogInfo ("\tTrying to parse " <> firstFile <> " as fastc--if it should be fasta specify 'fasta:' on input." <> "\n")
                                                                                                                    executeReadCommands'
                                                                                                                        ((fastcData, [fastcCharInfo]) : curData)
                                                                                                                        curGraphs
                                                                                                                        curTerminals
                                                                                                                        curExcludeList
                                                                                                                        curRenamePairs
                                                                                                                        curReBlockPairs
                                                                                                                        isPrealigned'
                                                                                                                        tcmPair
                                                                                                                        (tail argList)
                                                                                                        else
                                                                                                            let fastaData' = FAC.getFastAText fileContents firstFile isPrealigned'
                                                                                                            in  do
                                                                                                                    (fastaCharInfo, fastaData) ←
                                                                                                                        {-# SCC "getFastaCharInfo_1" #-} FAC.getFastaCharInfo fastaData' firstFile firstOption isPrealigned' tcmPair
                                                                                                                    logWith LogInfo ("\tTrying to parse " <> firstFile <> " as fasta" <> "\n")
                                                                                                                    executeReadCommands'
                                                                                                                        ((fastaData, [fastaCharInfo]) : curData)
                                                                                                                        curGraphs
                                                                                                                        curTerminals
                                                                                                                        curExcludeList
                                                                                                                        curRenamePairs
                                                                                                                        curReBlockPairs
                                                                                                                        isPrealigned'
                                                                                                                        tcmPair
                                                                                                                        (tail argList)
                                                                                            else do failWithPhase Parsing ("Cannot determine file type for " <> firstFile <> " need to prepend type")
                                            else -- fasta
                                                let firstChar = T.head $ T.dropWhile (== ' ') fileContents
                                                in
                                                if firstOption `elem` ["fasta", "nucleotide", "aminoacid", "hugeseq"]
                                                    then
                                                        let fastaData' = FAC.getFastAText fileContents firstFile isPrealigned'
                                                        in  do
                                                                (fastaCharInfo, fastaData) ←
                                                                    {-# SCC "getFastaCharInfo_2" #-} FAC.getFastaCharInfo fastaData' firstFile firstOption isPrealigned' tcmPair
                                                                executeReadCommands'
                                                                    ((fastaData, [fastaCharInfo]) : curData)
                                                                    curGraphs
                                                                    curTerminals
                                                                    curExcludeList
                                                                    curRenamePairs
                                                                    curReBlockPairs
                                                                    isPrealigned'
                                                                    tcmPair
                                                                    (tail argList)
                                                    else -- fastc

                                                        if firstOption == "fastc"
                                                            then
                                                                let fastcData = FAC.getFastCText fileContents firstFile isPrealigned'
                                                                in  do
                                                                        fastcCharInfo ← FAC.getFastcCharInfo fastcData firstFile isPrealigned' tcmPair
                                                                        executeReadCommands'
                                                                            ((fastcData, [fastcCharInfo]) : curData)
                                                                            curGraphs
                                                                            curTerminals
                                                                            curExcludeList
                                                                            curRenamePairs
                                                                            curReBlockPairs
                                                                            isPrealigned'
                                                                            tcmPair
                                                                            (tail argList)
                                                            else -- prealigned fasta

                                                                if firstOption `elem` ["prefasta", "prenucleotide", "preaminoacid", "prehugeseq"]
                                                                    then
                                                                        let fastaData' = FAC.getFastAText fileContents firstFile True
                                                                        in  do
                                                                                (fastaCharInfo, fastaData) ←
                                                                                    {-# SCC "getFastaCharInfo_3" #-} FAC.getFastaCharInfo fastaData' firstFile firstOption True tcmPair
                                                                                -- trace ("POSTREAD:" <> (show fastaCharInfo) <> "\n" <> (show fastaData))
                                                                                executeReadCommands'
                                                                                    ((fastaData, [fastaCharInfo]) : curData)
                                                                                    curGraphs
                                                                                    curTerminals
                                                                                    curExcludeList
                                                                                    curRenamePairs
                                                                                    curReBlockPairs
                                                                                    isPrealigned'
                                                                                    tcmPair
                                                                                    (tail argList)
                                                                    else -- prealigned fastc

                                                                        if firstOption == "prefastc"
                                                                            then
                                                                                let fastcData = FAC.getFastCText fileContents firstFile True
                                                                                in  do
                                                                                        fastcCharInfo ← FAC.getFastcCharInfo fastcData firstFile True tcmPair
                                                                                        executeReadCommands'
                                                                                            ((fastcData, [fastcCharInfo]) : curData)
                                                                                            curGraphs
                                                                                            curTerminals
                                                                                            curExcludeList
                                                                                            curRenamePairs
                                                                                            curReBlockPairs
                                                                                            isPrealigned'
                                                                                            tcmPair
                                                                                            (tail argList)
                                                                            else -- tnt
                                                                            -- tnt

                                                                                if firstOption == "tnt"
                                                                                    then do
                                                                                        tntData ← TNT.getTNTDataText fileContents firstFile
                                                                                        executeReadCommands'
                                                                                            (tntData : curData)
                                                                                            curGraphs
                                                                                            curTerminals
                                                                                            curExcludeList
                                                                                            curRenamePairs
                                                                                            curReBlockPairs
                                                                                            isPrealigned'
                                                                                            tcmPair
                                                                                            (tail argList)
                                                                                    else -- else if firstOption == "prealigned" then executeReadCommands' curData curGraphs curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)
                                                                                    -- FENEwick

                                                                                        if firstOption `elem` ["newick", "enewick", "fenewick"]
                                                                                            then
                                                                                                let thisGraphList = getFENewickGraphText fileContents
                                                                                                    hasLoops = filter id $ fmap B.hasLoop thisGraphList
                                                                                                    hasCycles = filter id $ fmap GFU.cyclic thisGraphList
                                                                                                in  if not $ null hasLoops
                                                                                                        then do failWithPhase Parsing ("Input graphin " <> firstFile <> "  has loops/self-edges")
                                                                                                        else
                                                                                                            if not $ null hasCycles
                                                                                                                then do failWithPhase Parsing ("Input graph in " <> firstFile <> " has at least one cycle")
                                                                                                                else
                                                                                                                    executeReadCommands'
                                                                                                                        curData
                                                                                                                        (thisGraphList <> curGraphs)
                                                                                                                        curTerminals
                                                                                                                        curExcludeList
                                                                                                                        curRenamePairs
                                                                                                                        curReBlockPairs
                                                                                                                        isPrealigned'
                                                                                                                        tcmPair
                                                                                                                        (tail argList)
                                                                                            else -- reading terminals list to include--must be "new" names if taxa are renamed

                                                                                                if firstOption == "include" then
                                                                                                    
                                                                                                    if (toLower firstChar) `elem` ['/','d','g','<','(', 'x']
                                                                                                            then do 
                                                                                                                failWithPhase Parsing ("Input 'include' file " <> firstFile <> " does not look like one.")
                                                                                                    else
                                                                                                        let terminalsList = fmap ((T.pack . filter (/= '"')) . filter C.isPrint) (words $ unlines $ U.stripComments (T.unpack <$> T.lines fileContents))
                                                                                                        in  executeReadCommands'
                                                                                                                curData
                                                                                                                curGraphs
                                                                                                                (terminalsList <> curTerminals)
                                                                                                                curExcludeList
                                                                                                                curRenamePairs
                                                                                                                curReBlockPairs
                                                                                                                isPrealigned'
                                                                                                                tcmPair
                                                                                                                (tail argList)
                                                                                                    else
                                                                                                        if firstOption == "exclude" then
                                                                                                                if (toLower firstChar) `elem` ['/','d','g','<','(', 'x']
                                                                                                                    then do 
                                                                                                                        failWithPhase Parsing ("Input 'exclude' file " <> firstFile <> " does not look like one.")
                                                                                                        else
                                                                                                                let excludeList = fmap ((T.pack . filter (/= '"')) . filter C.isPrint) (words $ unlines $ U.stripComments (T.unpack <$> T.lines fileContents))
                                                                                                                in  executeReadCommands'
                                                                                                                        curData
                                                                                                                        curGraphs
                                                                                                                        curTerminals
                                                                                                                        (excludeList <> curExcludeList)
                                                                                                                        curRenamePairs
                                                                                                                        curReBlockPairs
                                                                                                                        isPrealigned'
                                                                                                                        tcmPair
                                                                                                                        (tail argList)
                                                                                                            else
                                                                                                                if firstOption == "rename" then 
                                                                                                                    if (toLower firstChar) `elem` ['/','d','g','<','(', 'x']
                                                                                                                    then do 
                                                                                                                        failWithPhase Parsing ("Input 'rename' file " <> firstFile <> " does not look like one.")
                                                                                                                    else do
                                                                                                                        let renameLines = U.stripComments (T.unpack <$> T.lines fileContents)
                                                                                                                        namePairsLL ← mapM (makeNamePairs firstFile) renameLines
                                                                                                                        let namePairs = concat namePairsLL
                                                                                                                        let changeNamePairs = filter areDiff (namePairs <> curRenamePairs)
                                                                                                                        let newChangeNames = fmap fst changeNamePairs
                                                                                                                        let origChangeNames = fmap snd changeNamePairs
                                                                                                                        let intersectionChangeNames = L.nub $ L.intersect newChangeNames origChangeNames

                                                                                                                        if (not $ null intersectionChangeNames)
                                                                                                                            then errorWithoutStackTrace ("Error: Renaming of " <> (show intersectionChangeNames) <> " as both a new name and to be renamed")
                                                                                                                            else
                                                                                                                                executeReadCommands'
                                                                                                                                    curData
                                                                                                                                    curGraphs
                                                                                                                                    curTerminals
                                                                                                                                    curExcludeList
                                                                                                                                    (namePairs <> curRenamePairs)
                                                                                                                                    curReBlockPairs
                                                                                                                                    isPrealigned'
                                                                                                                                    tcmPair
                                                                                                                                    (tail argList)
                                                                                                                    else
                                                                                                                        if firstOption == "block" then
                                                                                                                            if (toLower firstChar) `elem` ['/','d','g','<','(', 'x']
                                                                                                                            then do 
                                                                                                                                failWithPhase Parsing ("Input 'block' file " <> firstFile <> " does not look like one.")
                                                                                                                            else do
                                                                                                                                let renameLines = U.stripComments (T.unpack <$> T.lines fileContents)
                                                                                                                                blockPairsLL ← mapM (makeNamePairs firstFile) renameLines
                                                                                                                                let blockPairs = concat blockPairsLL

                                                                                                                                executeReadCommands'
                                                                                                                                    curData
                                                                                                                                    curGraphs
                                                                                                                                    curTerminals
                                                                                                                                    curExcludeList
                                                                                                                                    curRenamePairs
                                                                                                                                    (blockPairs <> curReBlockPairs)
                                                                                                                                    isPrealigned'
                                                                                                                                    tcmPair
                                                                                                                                    (tail argList)
                                                                                                                            else -- else if firstOption == "prealigned" then
                                                                                                                            --    executeReadCommands' curData curGraphs curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)
                                                                                                                            do failWithPhase Parsing ("\n\n'Read' command error: option " <> firstOption <> " not recognized/implemented")
    where
        areDiff (a, b) =
            if a /= b
                then True
                else False


{- | makeNamePairs takes lines of rename or reblock file and returns pairs of names
for renaming or reblocking
-}
makeNamePairs ∷ String → String → PhyG [(T.Text, T.Text)]
makeNamePairs inFileName inLine =
    if null inLine
        then do return []
        else
            let lineWords = fmap T.strip $ T.words $ T.pack $ filter (/= '"') $ filter C.isPrint inLine
            in  if length lineWords < 2
                    then
                        failWithPhase
                            Parsing
                            ("Rename file " <> inFileName <> " line needs at least two Strings to rename the second as the first: " <> inLine)
                    else
                        let targetNameList = replicate (length $ tail lineWords) (head lineWords)
                            renamePairList = zip targetNameList (tail lineWords)
                        in  do
                                return renamePairList


{- | getReadArgs processes arguments ofr the 'read' command
should allow mulitple files and gracefully error check
also allows tcm file specification (limit 1 tcm per command?)
as fasta, fastc, tnt, tcm, prealigned
-}
getReadArgs ∷ String → [(String, String)] → PhyG [Argument]
getReadArgs fullCommand argList =
    if null argList
        then return []
        else -- trace ("GRA: " <> fullCommand <> " " <> (show argList)) (

            let (firstPart, secondPart) = head argList
            in  do
                    restPart ← getReadArgs fullCommand (tail argList)

                    -- command in wrong place like prealigned or rename after file name
                    if (not . null) firstPart && null secondPart
                        then
                            failWithPhase
                                Parsing
                                ( "\n\n'Read' command error: possibly incorrect placement of option specification '"
                                    <> firstPart
                                    <> "' should be before filename as in '"
                                    <> firstPart
                                    <> ":filename'\n"
                                )
                        else -- plain file name with no modifier

                            if null firstPart
                                then
                                    if (head secondPart == '"') || (last secondPart == '"')
                                        then do
                                            return $ (firstPart, init $ tail secondPart) : restPart
                                        else failWithPhase Parsing ("\n\n'Read' command error: '" <> secondPart <> "' : Need to specify filename in double quotes\n")
                                else -- Change to allowed modifiers

                                    if fmap toLower firstPart `notElem` V.readArgList
                                        then failWithPhase Parsing ("\n\n'Read' command error: " <> fullCommand <> " contains unrecognized option '" <> firstPart <> "'\n")
                                        else
                                            if null secondPart && (firstPart == "prealigned")
                                                then do
                                                    return $ (firstPart, []) : restPart
                                                else
                                                    if null (tail secondPart)
                                                        then return argList
                                                        else do
                                                            return $ (firstPart, init $ tail secondPart) : restPart


{- | processDotFile Function for preocessing dot files to allow for multiple graphs
-}
processDotFile :: String -> T.Text -> PhyG SimpleGraph
processDotFile fileName inGraphText = 
    if T.null inGraphText then pure LG.empty
    else do
        -- create temp file name
        randomInt <- getRandom
        let tempFileName = fileName <> (show $ abs (randomInt :: Int)) <> ".tmp" 
        tempHandle <- liftIO $ openFile tempFileName WriteMode
        liftIO $ hPutStr tempHandle (T.unpack inGraphText)
        liftIO $ hClose tempHandle
        tempHandle' <- liftIO $ openFile tempFileName ReadMode
        
        -- create dot graph format from temp file                
        dotGraph ← liftIO $ LG.hGetDotLocal tempHandle' -- fileHandle
        liftIO $ hClose tempHandle'
        
        -- delete temp file
        liftIO $ liftIO $  removeFile tempFileName

        let inputDot = GFU.relabelFGL $ LG.dotToGraph dotGraph
        let hasLoops = B.hasLoop inputDot
        if hasLoops
                then do failWithPhase Parsing ("Input graph in " <> fileName <> "  has loops/self-edges" <> "\n")
                else do logWith LogInfo ""
        let hasCycles = GFU.cyclic inputDot
        if hasCycles
                then do failWithPhase Parsing ("Input graph in " <> fileName <> " has at least one cycle" <> "\n")
                else pure inputDot




{- | splitDotGraphs splits dot graph input into a list of dot format graphs
    based around '{' and '}'
-}
splitDotGraphs :: T.Text -> [T.Text]
splitDotGraphs inText =
    if T.null inText then []
    else 
        let linesList = T.lines inText
            graphList = splitTextGraphLines [] $ filter notCommentLine linesList
        in
        graphList 

    where notCommentLine a = if T.null a then False
                             else if (T.head a == '/') && (T.head (T.tail a) == '/') then False
                             else True

{- | splitTextGraphLines takes line of Txt (list) and splits on lines with '}'
-}
splitTextGraphLines :: [T.Text] -> [T.Text] -> [T.Text]
splitTextGraphLines curGraphs inLines =
    if null inLines then curGraphs
    -- Terminal comment
    else if length inLines == 1 then curGraphs
    else 
        let firstGraphLines' = L.takeWhile ('}' `textNotElem`) inLines
            restLines' = L.dropWhile ('}' `textNotElem`) inLines
            (firstGraphLines, restLines) = if (not . null) restLines' then
                                                (firstGraphLines' <> [head restLines'], tail restLines')
                                            else error ("Error parsing multiple graphs in dot file:" <> (concat $ fmap T.unpack inLines))

            firstGraph = T.unlines firstGraphLines
        in
        splitTextGraphLines (firstGraph : curGraphs) restLines

    where textNotElem a b = if a `T.elem` b then False
                            else True

{- | getFENewickGraphText takes graph contents and returns local graph format
could be mulitple input graphs
-}
getFENewickGraphText ∷ T.Text → [LG.Gr T.Text Double]
getFENewickGraphText fileString =
    -- trace (fileString)
    LG.getFENLocal (T.filter (/= '\n') $ T.strip fileString)


{- | processTCMContents
functionality added to integerize tcm and add weight factor to allow for decimal tcm values
-}
processTCMContents ∷ Bool → String → String → PhyG ([ST.ShortText], [[Int]], Double)
processTCMContents indelGap inContents fileName =
    if null inContents
        then do failWithPhase Parsing ("\n\n'Read' 'tcm' command error: empty tcmfile `" <> fileName)
        else
            if indelGap
                then
                    let indelString = L.takeWhile (/= ',') inContents
                        substString = tail $ L.dropWhile (/= ',') inContents
                        indelMaybe = readMaybe indelString ∷ Maybe Int
                        substMaybe = readMaybe substString ∷ Maybe Int
                    in  if isNothing indelMaybe
                            then do failWithPhase Parsing ("Specification of indel cost must be an Integer (Indel cost, Substitution cost): " <> indelString)
                            else
                                if isNothing substMaybe
                                    then do
                                        failWithPhase Parsing ("Specification of substitution cost must be an Integer (Indel cost, Substitution cost): " <> substString)
                                    else do
                                        return ([], [[fromJust indelMaybe, fromJust substMaybe], []], 1.0)
                else -- First line is alphabet
                do
                    let tcmLines = lines inContents
                    let localAlphabet = fmap ST.pack $ words $ head tcmLines
                    -- to account for '-'
                    let numElements = 1 + length localAlphabet
                    let costMatrixStrings = words <$> tail tcmLines
                    -- localCostMatrix = filter (/= []) $ fmap (fmap (GU.stringToInt fileName)) costMatrixStrings
                    (scaleFactor, localCostMatrix) ← getCostMatrixAndScaleFactor fileName costMatrixStrings
                    let numLines = length localCostMatrix
                    let lengthRowList = fmap length localCostMatrix
                    let rowsCorrectLength = foldl' (&&) True $ fmap (== numElements) lengthRowList

                    if ST.singleton '-' `elem` localAlphabet
                        then do
                            failWithPhase
                                Parsing
                                "\n\n'Read' 'tcm' file format error: '-' (InDel/Gap) should not be specified as an alphabet element.  It is added in automatically and assumed to be the last row and column of tcm matrix"
                        else
                            if numLines /= numElements
                                then do
                                    failWithPhase
                                        Parsing
                                        ( "\n\n'Read' 'tcm' file format error: incorrect number of lines in matrix from "
                                            <> fileName
                                            <> " this could be due to a missing element symbol line at the beginning of the file (e.g. A C G T, '-' is assumed) or there is a mismatch in the dimensions of the tcm (including '-') "
                                            <> show numElements
                                            <> " elements are implied and there are "
                                            <> show numLines
                                            <> "\n" <> (concat $ L.intersperse " " $ words $ head tcmLines)
                                            <> "\n " <> (show $ length $ words $ head tcmLines)
                                        )
                                else
                                    if not rowsCorrectLength
                                        then do
                                            failWithPhase
                                                Parsing
                                                ( "\n\n'Read' 'tcm' file format error: incorrect lines length(s) in matrix from "
                                                    <> fileName
                                                    <> " there should be (including '-') "
                                                    <> show numElements
                                                )
                                        else -- trace (show (scaleFactor, localCostMatrix))
                                        -- trace (show localAlphabet <> "\n" <> show localCostMatrix)
                                            return (localAlphabet <> [ST.singleton '-'], localCostMatrix, scaleFactor)


{- | getCostMatrixAndScaleFactor' takes [[String]] and returns cost matrix as
[[Int]] but if there are decimal values, a scalerFactor is determined and the
cost matrix are integerized by multiplication by 1/scaleFactor
Unlike version above is more flexible with Double format

Adds precision to 100--fixes an issue when smallest and other values are similar
Not sure how this related to 2**16 (CShort) for ffi for <= 8 alphabets
Also can be integer overflow if the number are large since treated as integers 
    during DO for DNA--seems like pairwise alignment cost may be limited to 2^31 or 2^32
-}
getCostMatrixAndScaleFactor ∷ String → [[String]] → PhyG (Double, [[Int]])
getCostMatrixAndScaleFactor fileName = \case
    [] -> failWithPhase Parsing "Empty list in inStringListList"
    inStringListList ->
        let maxDecimalPlaces = maximum $ getDecimals <$> concat inStringListList
            doubleMatrix = filter (/= []) $ fmap (fmap (GU.stringToDouble fileName)) inStringListList
            minDouble = minimum $ fmap minimum $ fmap (filter (> 0.0)) doubleMatrix
            maxDouble = maximum $ fmap maximum doubleMatrix

            -- set precision based on range of max and min values
            -- to deal with some pathplogical situatinos with very high cost of no change and lo cost for change
            range = maxDouble / minDouble
            precisionFactor = if range < 10.0 then 10.0
                              else 1.0 
                              {-
                              if range >= 100000.0 then 1.0
                              else if range >= 10000.0 then 10.0
                              else if range >= 1000.0 then 100.0
                              else if range >= 100.0 then 1000.0
                              else if range >= 10.0 then 10000.0
                              else 100000.0
                              -}
            
            rescaledDoubleMatrix = fmap (fmap (* (precisionFactor / minDouble))) doubleMatrix
            integerizedMatrix = fmap (fmap round) rescaledDoubleMatrix

            scaleFactor
                | maxDecimalPlaces == 0 = 1.0
                | otherwise = (minDouble / precisionFactor)

            outputMatrix = case maxDecimalPlaces of
                0 -> filter (/= []) $ fmap (GU.stringToInt fileName) <$> inStringListList
                _ -> integerizedMatrix

        in  --trace ("GCMSC: " <> (show (precisionFactor,minDouble,scaleFactor,integerizedMatrix,rescaledDoubleMatrix) )) $
            pure $ (scaleFactor, outputMatrix)


{-
- | diagonalsEqualZero takes an integer matrix [[Int]] 
-- and checks if any diagonal values are == 0
diagonalsEqualZero :: [[Int]] -> Int -> Bool
diagonalsEqualZero inMatrix index =
    if null inMatrix then False
    else 
        if (head inMatrix) !! index /= 0 then True
        else diagonalsEqualZero (tail inMatrix) (index + 1)
-}

-- | getDecimals returns the number of decimal places in a string rep of a number
getDecimals ∷ String → Int
getDecimals inString =
    if null inString
        then 0
        else
            let decimalPart = dropWhile (/= '.') inString
            in  if null decimalPart
                    then 0
                    else length decimalPart - 1


{- Unused code

{- | getCostMatrixAndScaleFactor takes [[String]] and returns cost matrix as
[[Int]] but if there are decimal values, a scalerFactor is determined and the
cost matrix are integerized by multiplication by 1/scaleFactor
scaleFactor will alway be a factor of 10 to allow for easier
compilation of tcm characters later (weights the same...)
-}
getCostMatrixAndScaleFactor' ∷ String → [[String]] → PhyG (Double, [[Int]])
getCostMatrixAndScaleFactor' fileName inStringListList =
    if null inStringListList
        then failWithPhase Parsing "Empty list in inStringListList"
        else
            let maxDecimalPlaces = maximum $ getDecimals <$> concat inStringListList
                scaleFactor =
                    if maxDecimalPlaces == 0
                        then 1.0
                        else 0.1 ** fromIntegral maxDecimalPlaces
            in  if maxDecimalPlaces == 0
                    then do
                        return $ (scaleFactor, filter (/= []) $ fmap (fmap (GU.stringToInt fileName)) inStringListList)
                    else do
                        let newStringListList = fmap (fmap (integerizeString maxDecimalPlaces)) inStringListList
                        return $ (scaleFactor, filter (/= []) $ fmap (fmap (GU.stringToInt fileName)) newStringListList)

{- | integerizeString integerizes a String by removing the '.' and adding the number of '0's to pad out
adds maxDecimalPlaces - the nymber in the String
-}
integerizeString ∷ Int → String → String
integerizeString maxDecimalPlaces inString =
    if null inString
        then error "Null string in integerizeString"
        else
            let decimalPart = dropWhile (/= '.') inString
            in  if inString == "0"
                    then inString
                    else
                        if null decimalPart
                            then inString <> replicate maxDecimalPlaces '0'
                            else filter (/= '.') inString <> replicate (maxDecimalPlaces - (length decimalPart - 1)) '0'
-}