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

{-
Functions to peform the input file reading for PhyG
-}

module Input.ReadInputFiles
(  executeReadCommands
 , getReadArgs
 , extractInputTuple
 , expandReadCommands
) where

import Commands.Verify qualified as V
import Control.Evaluation
import Control.Monad.IO.Class (MonadIO (..))
import Control.Monad.Logger (LogLevel (..), Logger (..), Verbosity (..))
import Data.Char
import Data.Char qualified as C
import Data.Foldable
import Data.Graph.Inductive.Basic qualified as B
import Data.List qualified as L
import Data.Maybe
import Data.Text.Lazy qualified as T
import Data.Text.Lazy.IO qualified as TIO
import Data.Text.Short qualified as ST
import GeneralUtilities qualified as GU
import GraphFormatUtilities qualified as GFU
import Input.FastAC qualified as FAC
import Input.TNTUtilities qualified as TNT
import System.ErrorPhase (ErrorPhase (..))
import System.IO
import System.Path.Glob qualified as SPG
import Text.Read
import Types.Types
import Utilities.LocalGraph qualified as LG
import Utilities.Utilities qualified      as U
 
-- import Debug.Trace


-- | expandReadCommands expands read commands to multiple satisfying wild cards
-- read command can have multiple file names
expandReadCommands ::  [Command] -> Command -> PhyG [Command]
expandReadCommands _newReadList inCommand@(commandType, argList') =
    let argList = filter ((`notElem` ["tcm"]) . fst) argList'
        tcmArgList = filter ((`elem` ["tcm"]) . fst) argList'
        fileNames = fmap snd $ filter ((/="tcm") . fst) $ filter ((/= "") . snd) argList'
        modifierList = fmap fst argList
    in
    -- trace ("ERC: " <> (show fileNames)) (
    if commandType /= Read then error ("Incorrect command type in expandReadCommands: " <> show inCommand)
    else do
        globbedFileNames <- liftIO $ mapM SPG.glob fileNames
        if all null globbedFileNames then do failWithPhase Inputting ("File(s) not found in 'read' command (could be due to incorrect filename or missing closing double quote '\"''): " <> show fileNames)
        else do
            newArgPairs <- mapM makeNewArgs (zip modifierList globbedFileNames)

            let commandList = replicate (length newArgPairs) commandType
            let tcmArgListL = replicate (length newArgPairs) tcmArgList
            
            return $ zip commandList (zipWith (<>) newArgPairs tcmArgListL)
    -- )

-- | makeNewArgs takes an argument modifier (first in pair) and replicates is and zips with
-- globbed file name list to create a list of arguments
makeNewArgs :: (String, [String]) -> PhyG [(String, String)]
makeNewArgs (modifier, fileNameList) =
    if null fileNameList then do failWithPhase Inputting ("Null filename list in makeNewArgs: " <> show (modifier, fileNameList))
    else let modList = replicate (length fileNameList) modifier
         in do return $ zip modList fileNameList

{-
Added strictness on read so could close files after reading to
   allow for large number (1000's) of files to be input.
   this may not be the way to go, especially if files are big and
   closing them is not an issue
-}

-- | extractInputTuple  takes the list of pairs from mapped executeReadCommands
-- and returns ([RawData], [SimpleGraph])
extractInputTuple :: [([RawData], [SimpleGraph], [NameText], [NameText], [(NameText, NameText)], [(NameText, NameText)])]
                  -> ([RawData], [SimpleGraph], [NameText], [NameText], [(NameText, NameText)], [(NameText, NameText)])
extractInputTuple dataGraphList =
    let (inDataList, inGraphList, inTerminalsList, inExcludeList, inRenamePairs, inReBlockPairs) = L.unzip6 dataGraphList
        rawData   = L.sort $ concat inDataList
        rawGraphs = concat inGraphList
        rawTerminals = concat inTerminalsList
        excludeTerminals = concat inExcludeList
        renamePairs = concat inRenamePairs
        reBlockPairs = concat inReBlockPairs
    in  (rawData, rawGraphs, rawTerminals, excludeTerminals, renamePairs, reBlockPairs)


-- | executeReadCommands wrapper for executeReadCommands' with out all the empty list imput
executeReadCommands :: [Argument]
                    -> PhyG ([RawData], [SimpleGraph], [NameText], [NameText], [(NameText, NameText)], [(NameText, NameText)])
executeReadCommands = executeReadCommands' [] [] [] [] [] [] False ([],[],1.0)

-- | executeReadCommands' reads input files and returns raw data, input graphs, and terminal taxa to include
-- assumes that "prealigned" and "tcm:File" are the first arguments if they are specified
-- so that they can apply to all the files in the command without order depence
executeReadCommands' :: [RawData]
                    -> [SimpleGraph]
                    -> [NameText]
                    -> [NameText]
                    -> [(NameText, NameText)]
                    -> [(NameText, NameText)]
                    -> Bool
                    -> ([ST.ShortText], [[Int]], Double)
                    -> [Argument]
                    -> PhyG ([RawData], [SimpleGraph], [NameText], [NameText], [(NameText, NameText)], [(NameText, NameText)])
executeReadCommands' curData curGraphs curTerminals curExcludeList curRenamePairs curReBlockPairs _ tcmPair argList = do
    if null argList then return (curData, curGraphs, curTerminals, curExcludeList, curRenamePairs, curReBlockPairs)
    else do
        -- logWith LogInfo ("ERC: " <> (show argList) <> " " <> (show tcmPair))
        let isPrealigned'= False  --removed this option and m,asde specific to sequence types
        --      | isPrealigned = True
        --      | "prealigned" `elem`  fmap fst argList = True
        --      | otherwise = False
        let (firstOption', firstFile) = head argList
        let firstOption = fmap C.toLower firstOption'
        --Check for prealigned
        -- if firstOption == "prealigned" then
        --    executeReadCommands' curData curGraphs curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)
        -- else do
        do
            fileHandle <- liftIO $ if ',' `notElem` firstFile then openFile firstFile ReadMode
                                   else return (stdin :: Handle)
            canBeReadFrom <- liftIO $ hIsReadable fileHandle
            if not canBeReadFrom then do failWithPhase Inputting ("\n\n'Read' error: file " <> firstFile <> " cannot be read")
            else
                if not $ null firstOption then do logWith LogInfo ("Reading " <> firstFile <> " with option " <> firstOption)
                else do logWith LogInfo ("Reading " <> firstFile <> " with no options")

            -- this is awkward but need to use dot utilities
            if firstOption == "dot" then do
                dotGraph <- liftIO $ LG.hGetDotLocal fileHandle
                let inputDot = GFU.relabelFGL $ LG.dotToGraph dotGraph
                let hasLoops = B.hasLoop inputDot
                if hasLoops then do failWithPhase Parsing  ("Input graph in " <> firstFile <> "  has loops/self-edges")
                else do logWith LogInfo ""
                let hasCycles = GFU.cyclic inputDot
                if hasCycles then do failWithPhase Parsing("Input graph in " <> firstFile <> " has at least one cycle")
                else executeReadCommands' curData (inputDot : curGraphs) curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)
            -- not "dot" files
            else do
                -- logWith LogInfo ("FC1: " <> firstFile)
                fileContents <- liftIO $ if ',' `notElem` firstFile then TIO.hGetContents fileHandle
                                         else return (T.pack firstFile)

                {--Strict version--don't think necessary--but could getvinot open fikle number issues?
                 -- destroys lazyness but allows closing right away
                 -- this so don't have huge numbers of open files for large data sets
                fileContents <- hGetContents' fileHandle
                hClose fileHandle
                -}

                if T.null fileContents then do failWithPhase Inputting ("Error: Input file " <> firstFile <> " is empty")
                else
                    -- try to figure out file type based on first and following characters
                    if firstOption == "tcm" then do
                        newTCMPair <- processTCMContents (',' `elem` firstFile) (T.unpack fileContents) firstFile
                        executeReadCommands' curData curGraphs curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' newTCMPair (tail argList)
                    else if null firstOption then
                        let firstChar = T.head $ T.dropWhile (== ' ') fileContents
                        in
                        if (toLower firstChar == '/') || (toLower firstChar == 'd') || (toLower firstChar == 'g')then do
                            logWith LogInfo ("\tTrying to parse " <> firstFile <> " as dot")
                            -- destroys lazyness but allows closing right away
                            -- this so don't have huge numbers of open files for large data sets
                            fileHandle2 <- liftIO $ openFile  firstFile ReadMode
                            dotGraph    <- liftIO $ LG.hGetDotLocal fileHandle2
                            liftIO $ hClose fileHandle2
                            let inputDot = GFU.relabelFGL $ LG.dotToGraph dotGraph
                            let hasLoops = B.hasLoop inputDot
                            if hasLoops then do failWithPhase Parsing ("Input graph in " <> firstFile <> "  has loops/self-edges")
                            else do logWith LogInfo ""
                            let hasCycles = GFU.cyclic inputDot
                            if hasCycles then do failWithPhase Parsing ("Input graph in " <> firstFile <> " has at least one cycle")
                            else executeReadCommands' curData (inputDot : curGraphs) curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)
                        else if (toLower firstChar == '<') || (toLower firstChar == '(')  then
                            let thisGraphList = getFENewickGraphText fileContents
                                hasCycles = filter id $ fmap GFU.cyclic thisGraphList
                                hasLoops = filter id $ fmap B.hasLoop thisGraphList
                            in
                            if not $ null hasLoops then do failWithPhase Parsing ("Input graph in " <> firstFile <> "  has loops/self-edges")
                            else if not $ null hasCycles then do failWithPhase Parsing ("Input graph in " <> firstFile <> " has at least one cycle")
                            else executeReadCommands' curData (thisGraphList <> curGraphs) curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)

                        else if toLower firstChar == 'x' then
                            let tntData = TNT.getTNTDataText fileContents firstFile
                            in do
                            logWith LogInfo ("\tTrying to parse " <> firstFile <> " as TNT")
                            executeReadCommands' (tntData : curData) curGraphs curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)
                        else
                            let fileContents' =  T.unlines $ filter (not . T.null) $ T.takeWhile (/= ';') <$> T.lines fileContents
                            in
                              if T.null (T.dropWhile (== ' ') fileContents') then do failWithPhase Parsing ("Null file '" <> firstFile <> "' input")
                              else if T.head (T.dropWhile (== ' ') fileContents') == '>' then
                                let secondLine = T.lines fileContents' !! 1
                                    hasSpaces = T.find (== ' ') secondLine
                                    -- numWords = length $ words secondLine
                                in
                                -- spaces between alphabet elements suggest fastc
                                if isJust hasSpaces then
                                    let fastcData = FAC.getFastCText fileContents firstFile isPrealigned'
                                        fastcCharInfo = FAC.getFastcCharInfo fastcData firstFile isPrealigned' tcmPair
                                    in do
                                    logWith LogInfo ("\tTrying to parse " <> firstFile <> " as fastc--if it should be fasta specify 'fasta:' on input.")
                                    executeReadCommands' ((fastcData, [fastcCharInfo]) : curData) curGraphs curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)
                                else
                                    let fastaData' = FAC.getFastAText fileContents firstFile isPrealigned'
                                        (fastaCharInfo, fastaData)  = FAC.getFastaCharInfo fastaData' firstFile firstOption isPrealigned' tcmPair
                                    in do
                                    logWith LogInfo ("\tTrying to parse " <> firstFile <> " as fasta")
                                    executeReadCommands' ((fastaData, [fastaCharInfo]) : curData) curGraphs curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)

                        else do failWithPhase Parsing ("Cannot determine file type for " <> firstFile <> " need to prepend type")
                    -- fasta
                    else if firstOption `elem` ["fasta", "nucleotide", "aminoacid", "hugeseq"] then
                        let fastaData' = FAC.getFastAText fileContents firstFile isPrealigned'
                            (fastaCharInfo, fastaData) = FAC.getFastaCharInfo fastaData' firstFile firstOption isPrealigned' tcmPair
                        in
                        executeReadCommands' ((fastaData, [fastaCharInfo]) : curData) curGraphs curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)
                    -- fastc
                    else if firstOption == "fastc"  then
                        let fastcData = FAC.getFastCText fileContents firstFile isPrealigned'
                            fastcCharInfo = FAC.getFastcCharInfo fastcData firstFile isPrealigned' tcmPair
                        in
                        executeReadCommands' ((fastcData, [fastcCharInfo]) : curData) curGraphs curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)

                    --prealigned fasta
                    else if firstOption `elem` ["prefasta", "prenucleotide", "preaminoacid", "prehugeseq"] then
                        let fastaData' = FAC.getFastAText fileContents firstFile True
                            (fastaCharInfo, fastaData) = FAC.getFastaCharInfo fastaData' firstFile firstOption True tcmPair
                        in
                        -- trace ("POSTREAD:" <> (show fastaCharInfo) <> "\n" <> (show fastaData))
                        executeReadCommands' ((fastaData, [fastaCharInfo]) : curData) curGraphs curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)

                    -- prealigned fastc
                    else if firstOption == "prefastc"  then
                        let fastcData = FAC.getFastCText fileContents firstFile True
                            fastcCharInfo = FAC.getFastcCharInfo fastcData firstFile True tcmPair
                        in
                        executeReadCommands' ((fastcData, [fastcCharInfo]) : curData) curGraphs curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)
                    -- tnt
                    -- tnt
                    else if firstOption == "tnt" then
                        let tntData = TNT.getTNTDataText fileContents firstFile
                        in
                        executeReadCommands' (tntData : curData) curGraphs curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)
                    -- else if firstOption == "prealigned" then executeReadCommands' curData curGraphs curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)
                    -- FENEwick
                    else if firstOption `elem` ["newick" , "enewick", "fenewick"]  then
                        let thisGraphList = getFENewickGraphText fileContents
                            hasLoops = filter id $ fmap B.hasLoop thisGraphList
                            hasCycles = filter id $ fmap GFU.cyclic thisGraphList
                        in
                        if not $ null hasLoops then do failWithPhase Parsing ("Input graphin " <> firstFile <> "  has loops/self-edges")
                        else if not $ null hasCycles then do failWithPhase Parsing ("Input graph in " <> firstFile <> " has at least one cycle")
                        else executeReadCommands' curData (thisGraphList <> curGraphs) curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)

                    -- reading terminals list to include--must be "new" names if taxa are renamed
                    else if firstOption == "include"  then
                        let terminalsList = fmap ((T.pack . filter (/= '"')) . filter C.isPrint) (words $ unlines $ U.stripComments (T.unpack <$> T.lines fileContents))
                        in
                        executeReadCommands' curData curGraphs (terminalsList <> curTerminals) curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)
                    else if firstOption == "exclude"  then
                        let excludeList = fmap ((T.pack . filter (/= '"')) . filter C.isPrint) (words $ unlines $ U.stripComments (T.unpack <$> T.lines fileContents))
                        in
                        executeReadCommands' curData curGraphs curTerminals (excludeList <> curExcludeList) curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)
                    else if firstOption == "rename"  then do
                        let renameLines = U.stripComments (T.unpack <$> T.lines fileContents)
                        namePairsLL <- mapM (makeNamePairs firstFile) renameLines
                        let namePairs = concat namePairsLL
                        
                        executeReadCommands' curData curGraphs curTerminals curExcludeList (namePairs <> curRenamePairs) curReBlockPairs isPrealigned' tcmPair (tail argList)
                    else if firstOption == "block"  then do
                        let renameLines = U.stripComments (T.unpack <$> T.lines fileContents)
                        blockPairsLL <- mapM (makeNamePairs firstFile) renameLines
                        let blockPairs = concat blockPairsLL
                        
                        executeReadCommands' curData curGraphs curTerminals curExcludeList curRenamePairs (blockPairs <> curReBlockPairs) isPrealigned' tcmPair (tail argList)
                    --else if firstOption == "prealigned" then
                    --    executeReadCommands' curData curGraphs curTerminals curExcludeList curRenamePairs curReBlockPairs isPrealigned' tcmPair (tail argList)
                    else do failWithPhase Parsing ("\n\n'Read' command error: option " <> firstOption <> " not recognized/implemented")

-- | makeNamePairs takes lines of rename or reblock file and returns pairs of names
-- for renaming or reblocking
makeNamePairs :: String -> String -> PhyG [(T.Text, T.Text)]
makeNamePairs inFileName inLine =
    if null inLine then do return []
    else
        let lineWords = fmap T.strip $ T.words $ T.pack $ filter (/= '"') $ filter C.isPrint inLine
        in
        if length lineWords < 2 then failWithPhase Parsing ("Rename file " <> inFileName <> " line needs at least two Strings to rename the second as the first: " <> inLine)
        else
            let targetNameList = replicate (length $ tail lineWords) (head lineWords)
                renamePairList = zip targetNameList (tail lineWords)
            in do
            return renamePairList



-- | getReadArgs processes arguments ofr the 'read' command
-- should allow mulitple files and gracefully error check
-- also allows tcm file specification (limit 1 tcm per command?)
-- as fasta, fastc, tnt, tcm, prealigned
getReadArgs :: String -> [(String, String)] -> PhyG [Argument]
getReadArgs fullCommand argList =
    if null argList then return []
    else
        -- trace ("GRA: " <> fullCommand <> " " <> (show argList)) (
        let (firstPart, secondPart) = head argList
        in do
        restPart <- getReadArgs fullCommand (tail argList)
            
        -- command in wrong place like prealigned or rename after file name
        if (not . null) firstPart && null secondPart then
            failWithPhase Parsing ("\n\n'Read' command error: possibly incorrect placement of option specification '" <> firstPart
                <> "' should be before filename as in '" <> firstPart <> ":filename'\n")

        -- plain file name with no modifier
        else if null firstPart then
            if (head secondPart == '"') || (last secondPart == '"') then do
                return $ (firstPart, init $ tail secondPart) : restPart
            else failWithPhase Parsing  ("\n\n'Read' command error: '" <> secondPart <>"' : Need to specify filename in double quotes\n")

        -- Change to allowed modifiers
        else if fmap toLower firstPart `notElem` V.readArgList then
            failWithPhase Parsing  ("\n\n'Read' command error: " <> fullCommand <> " contains unrecognized option '" <> firstPart <> "'\n")

        else if null secondPart && (firstPart == "prealigned")  then do
            return $ (firstPart, []) : restPart

        else if null (tail secondPart) then return argList

        else do 
            return $ (firstPart, init $ tail secondPart) : restPart
        -- )

{-  Not used when migrated to Text input from String

-- | getFENewickGraph takes graph contents and returns local graph format
-- could be mulitple input graphs
getFENewickGraph :: String -> [LG.Gr T.Text Double]
getFENewickGraph fileString =
    getFENewickGraphText (T.pack fileString)
-}

-- | getFENewickGraphText takes graph contents and returns local graph format
-- could be mulitple input graphs
getFENewickGraphText :: T.Text -> [LG.Gr T.Text Double]
getFENewickGraphText fileString =
    -- trace (fileString)
    LG.getFENLocal (T.filter (/= '\n') $ T.strip fileString)

-- | processTCMContents
-- functionality added to integerize tcm and add weight factor to allow for decimal tcm values
processTCMContents :: Bool -> String -> String -> PhyG ([ST.ShortText], [[Int]], Double)
processTCMContents indelGap inContents fileName =
    if null inContents then do failWithPhase Parsing ("\n\n'Read' 'tcm' command error: empty tcmfile `" <> fileName)
    else
        if indelGap then
            let indelString = L.takeWhile (/= ',') inContents
                substString = tail $ L.dropWhile (/= ',') inContents
                indelMaybe = readMaybe indelString :: Maybe Int
                substMaybe = readMaybe substString :: Maybe Int
            in
            if isNothing indelMaybe then do failWithPhase Parsing ("Specification of indel cost must be an Integer (Indel cost, Substitution cost): " <> indelString)
            else if isNothing substMaybe then do failWithPhase Parsing ("Specification of substitution cost must be an Integer (Indel cost, Substitution cost): " <> substString)
            else do
                return ([], [[fromJust indelMaybe, fromJust substMaybe],[]], 1.0)
        --First line is alphabet
        else do
            let tcmLines = lines inContents
            let localAlphabet =  fmap ST.pack $ words $ head tcmLines
            -- to account for '-'
            let numElements = 1 + length localAlphabet
            let costMatrixStrings = words <$> tail tcmLines
            -- localCostMatrix = filter (/= []) $ fmap (fmap (GU.stringToInt fileName)) costMatrixStrings
            (scaleFactor, localCostMatrix) <- getCostMatrixAndScaleFactor fileName costMatrixStrings
            let numLines = length localCostMatrix
            let lengthRowList = fmap length localCostMatrix
            let rowsCorrectLength = foldl' (&&) True $ fmap (==  numElements) lengthRowList
            
            if ST.singleton '-' `elem` localAlphabet then do failWithPhase Parsing "\n\n'Read' 'tcm' file format error: '-' (InDel/Gap) should not be specifid as an alphabet element.  It is added in automatically and assumed to be the last row and column of tcm matrix"
            else if numLines /= numElements then do 
                failWithPhase Parsing ("\n\n'Read' 'tcm' file format error: incorrect number of lines in matrix from " <> fileName <> " this could be due to a missing element symbol line at the beginning of the file (e.g. A C G T, '-' is assumed) or there is a mismatch in the dimensions of the tcm (including '-') " <> show numElements <> " elements are implied and there are " <> show numLines)
            else if not rowsCorrectLength then do 
                failWithPhase Parsing ("\n\n'Read' 'tcm' file format error: incorrect lines length(s) in matrix from " <> fileName <> " there should be (including '-') " <> show numElements)
            else
                --trace (show (scaleFactor, localCostMatrix))
                -- trace (show localAlphabet <> "\n" <> show localCostMatrix)
                return (localAlphabet <> [ST.singleton '-'], localCostMatrix, scaleFactor)


-- | getCostMatrixAndScaleFactor takes [[String]] and returns cost matrix as
-- [[Int]] but if there are decimal values, a scalerFactor is determined and the
-- cost matrix are integerized by multiplication by 1/scaleFactor
-- scaleFactor will alway be a factor of 10 to allow for easier
-- compilation of tcm charcaters later (weights the same...)
getCostMatrixAndScaleFactor :: String -> [[String]] -> PhyG (Double, [[Int]])
getCostMatrixAndScaleFactor fileName inStringListList =
    if null inStringListList then failWithPhase Parsing "Empty list in inStringListList"
    else
        let maxDecimalPlaces = maximum $ getDecimals <$> concat inStringListList
            scaleFactor = if maxDecimalPlaces == 0 then 1.0
                          else 0.1 ** fromIntegral maxDecimalPlaces
            in
            if maxDecimalPlaces == 0 then do
                return $ (scaleFactor, filter (/= []) $ fmap (fmap (GU.stringToInt fileName)) inStringListList)
            else do
                let newStringListList = fmap (fmap (integerizeString maxDecimalPlaces)) inStringListList
                return $ (scaleFactor, filter (/= []) $ fmap (fmap (GU.stringToInt fileName)) newStringListList)


-- | getDecimals returns the number of decimal places in a string rep of a number
getDecimals :: String -> Int
getDecimals inString =
    if null inString then 0
    else
        let decimalPart = dropWhile (/= '.') inString
        in
        if null decimalPart then 0
        else length decimalPart - 1

-- | integerizeString integerizes a String by removing the '.' and adding the number of '0's to pad out
-- adds maxDecimalPlaces - the nymber in the String
integerizeString :: Int -> String -> String
integerizeString maxDecimalPlaces inString =
    if null inString then error "Null string in integerizeString"
    else
        let decimalPart = dropWhile (/= '.') inString
        in
        if inString == "0" then inString
        else if null decimalPart then inString <> replicate maxDecimalPlaces '0'
        else filter (/= '.') inString <> replicate (maxDecimalPlaces - (length decimalPart - 1)) '0'
