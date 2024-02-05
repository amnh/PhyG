{-  To add commands:
    1) add new command in Types.hs
    2) add string name of command to allowedCommandList
    3) add instruction processing to getInstruction
    4) add argument processing to parseCommandArg
    5) add argument processing function
        with meaningful errors
    6) Add amchinery of command in general code
-}

{- |
Module to process commands.
-}
module Commands.ProcessCommands (
    expandRunCommands,
    getCommandList,
    movePrealignedTCM,
    preprocessOptimalityCriteriaScripts,
) where

import Commands.Verify qualified as V
import Control.Monad.IO.Class (MonadIO (..))
import Data.Char
import Data.Foldable
import Data.List qualified as L
import Data.List.Split
import Data.Maybe
import GeneralUtilities
import Input.ReadInputFiles qualified as RIF
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import Types.Types


{- | preprocessOptimalityCriteriaScripts takes a processed command list and
processes for optimlity criteria that change tcms and such for
Parsimony, PMDL, SI, MAPA, and NCM
-}
preprocessOptimalityCriteriaScripts ∷ [Command] → [Command]
preprocessOptimalityCriteriaScripts inCommandList = inCommandList


{- | expandRunCommands takes raw coomands and if a "run" command is found it reads that file
and adds those commands in place
ensures one command per line
-}
expandRunCommands ∷ [String] → [String] → PhyG [String]
expandRunCommands curLines inLines =
    -- trace ("EXP " <> (show curLines) <> show inLines) (
    if null inLines
        then return $ reverse curLines
        else do
            let firstLineRead = removeComments [filter (/= ' ') $ head inLines]
            (firstLine, restLine) ←
                if null firstLineRead
                    then return ([], [])
                    else splitCommandLine $ head firstLineRead

            let leftParens = length $ filter (== '(') firstLine
            let rightParens = length $ filter (== ')') firstLine

            -- trace ("FL " <> firstLine) (
            -- only deal with run lines
            if leftParens /= rightParens
                then do failWithPhase Parsing ("Command line with unbalances parens '()': " <> firstLine <> "\n")
                else
                    if null firstLine
                        then expandRunCommands curLines (tail inLines)
                        else
                            if take 3 (fmap toLower firstLine) /= "run"
                                then expandRunCommands (firstLine : curLines) (restLine : tail inLines)
                                else do
                                    -- is a "run command"
                                    parsedFirst ← parseCommand firstLine
                                    let (_, runFileList) = head parsedFirst
                                    runFileNames ← mapM (checkFileNames . snd) runFileList
                                    fileListContents ← liftIO $ mapM readFile runFileNames
                                    let newLines = concatMap lines fileListContents
                                    expandRunCommands (newLines <> curLines) (restLine : tail inLines)


-- )
-- )

{- | splitCommandLine takes a line with potentially multiple commands and splits
between the first command and all others.
-}
splitCommandLine ∷ String → PhyG (String, String)
splitCommandLine inLine =
    if null inLine
        then return ([], [])
        else
            let leftParens = length $ filter (== '(') inLine
                rightParens = length $ filter (== ')') inLine
                firstPart = takeWhile (/= '(') inLine
                parenPart = getBalancedParenPart "" (dropWhile (/= '(') inLine) 0 0
                firstCommand = firstPart <> parenPart
                restPart = drop (length firstCommand) inLine
            in  if leftParens /= rightParens
                    then do failWithPhase Parsing ("Command line with unbalances parens '()': " <> inLine <> "\n")
                    else return (firstCommand, restPart)


-- | checkFileNames checks if first and last element of String are double quotes and removes them
checkFileNames ∷ String → PhyG String
checkFileNames inName
    | null inName = do failWithPhase Parsing "Error: Null file name"
    | head inName /= '"' = do failWithPhase Parsing ("Error: File name must be in double quotes (b): " <> inName <> "\n")
    | last inName /= '"' = do failWithPhase Parsing ("Error: File name must be in double quotes (e): " <> inName <> "\n")
    | otherwise = return $ init $ tail inName


{- | getCommandList takes a String from a file and returns a list of commands and their arguments
these are syntactically verified, but any input files are not checked
commands in lines one to a line
-}
getCommandList ∷ [String] → PhyG [Command]
getCommandList rawContents =
    if null rawContents
        then do
            -- errorWithoutStackTrace "Error: Empty command file"
            failWithPhase Parsing "Empty command file"
        else do
            let rawList = removeComments $ fmap (filter (/= ' ')) rawContents
            -- expand for read wildcards here--cretge a new, potentially longer list
            parsedRaw ← mapM parseCommand rawList
            let processedCommands = concat parsedRaw

            let reportCommands = filter ((== Report) . fst) processedCommands

            let reportGraphsArgs = filter (== "graphs") $ fmap (fmap toLower) $ fmap fst $ concatMap snd reportCommands
            let reportNewickArgs = filter (== "newick") $ fmap (fmap toLower) $ fmap fst $ concatMap snd reportCommands
            let reportDotArgs = filter (== "dot") $ fmap (fmap toLower) $ fmap fst $ concatMap snd reportCommands

            if null reportCommands || null (reportGraphsArgs <> reportNewickArgs <> reportDotArgs)
                then -- trace ("Warning: No reporting of resulting graphs is specified.  Adding default report graph file 'defaultGraph.dot'") $

                    let addedReport = (Report, [("graphs", []), ([], "_defaultGraph.dot_"), ("dotpdf", [])])
                    in  do
                            logWith LogWarn "Warning: No reporting of resulting graphs is specified.  Adding default report graph file 'defaultGraph.dot'\n"
                            return (processedCommands <> [addedReport])
                else -- trace (show rawList)
                do
                    return processedCommands


{- | removeComments deletes anyhting on line (including line)
after double dash "--"
-}
removeComments ∷ [String] → [String]
removeComments = \case
    [] → []
    firstLine : otherLines → case firstLine of
        [] → removeComments otherLines
        '-' : '-' : _ → removeComments otherLines
        _ →
            let nonComment = case splitOn "--" firstLine of
                    [] → firstLine
                    prefix : _ → filter isPrint prefix
            in  nonComment : removeComments otherLines


{- | getInstruction returns the command type from an input String
all operations on lower case
-}
getInstruction ∷ String → [String] → PhyG Instruction
getInstruction inString possibleCommands
    | null inString = do failWithPhase Parsing "Empty command String"
    | fmap toLower inString == "build" = return Build
    | fmap toLower inString == "fuse" = return Fuse
    | fmap toLower inString == "read" = return Read
    | fmap toLower inString == "reblock" = return Reblock
    | fmap toLower inString == "refine" = return Refine
    | fmap toLower inString == "rename" = return Rename
    | fmap toLower inString == "report" = return Report
    | fmap toLower inString == "run" = return Run
    | fmap toLower inString == "search" = return Search
    | fmap toLower inString == "select" = return Select
    | fmap toLower inString == "set" = return Set
    | fmap toLower inString == "support" = return Support
    | fmap toLower inString == "swap" = return Swap
    | fmap toLower inString == "transform" = return Transform
    | otherwise =
        let errorMatch = snd $ getBestMatch (maxBound ∷ Int, "no suggestion") possibleCommands inString
        in  do
                failWithPhase Parsing $
                    fold
                        ["\nError: Unrecognized command. By \'", inString, "\' did you mean \'", errorMatch, "\'?\n"]


{- | parseCommand takes a command file line and processes the String into a command and its arguemnts
assumes single command per line
-}
parseCommand ∷ String → PhyG [Command]
parseCommand = \case
    [] → pure []
    inLine →
        let (firstString, restString) = getSubCommand inLine False
            instructionString = takeWhile (/= '(') firstString -- inLine
        in  do
                -- this does not allow recursive multi-option arguments
                -- NEED TO FIX
                -- make in to a more sophisticated split outside of parens
                argList ← argumentSplitter inLine $ init $ tail $ dropWhile (/= '(') $ filter (/= ' ') firstString

                localInstruction ← getInstruction instructionString V.allowedCommandList
                processedArg ← parseCommandArg firstString localInstruction argList
                parsedRest ← parseCommand restString

                pure $ (localInstruction, processedArg) : parsedRest


{- | getSubCommand takes a string and extracts the first occurrence of the
structure bleh(...), and splits the string on that, the sub command can contain
parens and commas
-}
getSubCommand ∷ String → Bool → (String, String)
getSubCommand inString hasComma =
    if null inString
        then ([], [])
        else
            let firstPart = takeWhile (/= '(') inString
                secondPart = dropWhile (/= '(') inString
                parenPart = getBalancedParenPart "" secondPart 0 0
                incrCounter = if hasComma then 1 else 0
                remainderPart = drop (length (firstPart <> parenPart) + incrCounter) inString -- to remove ','
            in  (firstPart <> parenPart, remainderPart)


{- | getBalancedParenPart stakes a string starting with '(' and takes all
characters until (and including) the balancing ')'
call with  getBalancedParenPart "" inString 0 0
-}
getBalancedParenPart ∷ String → String → Int → Int → String
getBalancedParenPart curString inString countLeft countRight =
    if null inString
        then reverse curString
        else
            let firstChar = head inString
            in  if firstChar == '('
                    then getBalancedParenPart (firstChar : curString) (tail inString) (countLeft + 1) countRight
                    else
                        if firstChar == ')'
                            then
                                if countLeft == countRight + 1
                                    then reverse (firstChar : curString)
                                    else getBalancedParenPart (firstChar : curString) (tail inString) countLeft (countRight + 1)
                            else getBalancedParenPart (firstChar : curString) (tail inString) countLeft countRight


{- | argumentSplitter takes argument string and returns individual strings of arguments
which can include null, single, multiple or sub-command arguments
these are each pairs of an option string (could be null) and a subarguments String (also could be null)
-}
argumentSplitter ∷ String → String → PhyG [(String, String)]
argumentSplitter commandLineString inString
    | null inString = return []
    | otherwise =
        let commaIndex = fromMaybe (maxBound ∷ Int) (L.elemIndex ',' inString)
            semiIndex = fromMaybe (maxBound ∷ Int) (L.elemIndex ':' inString)
            leftParenIndex = fromMaybe (maxBound ∷ Int) (L.elemIndex '(' inString)
            firstDivider = minimum [commaIndex, semiIndex, leftParenIndex]
        in  do
                checkErrors ← freeOfSimpleErrors inString
                if not checkErrors
                    then do failWithPhase Parsing "Error in command specification format\n"
                    else -- simple no argument arg

                        if firstDivider == (maxBound ∷ Int)
                            then
                                if head inString == '"'
                                    then return [([], inString)]
                                    else return [(inString, [])]
                            else
                                if commaIndex == firstDivider
                                    then -- no arg

                                        if null (take firstDivider inString)
                                            then do failWithPhase Parsing ("Error in command '" <> commandLineString <> "' perhaps due to extraneous commas (',')\n")
                                            else
                                                if head (take firstDivider inString) == '"'
                                                    then do
                                                        restPart ← argumentSplitter commandLineString (drop (firstDivider + 1) inString)
                                                        return $ ([], take firstDivider inString) : restPart
                                                    else do
                                                        restPart ← argumentSplitter commandLineString (drop (firstDivider + 1) inString)
                                                        return $ (take firstDivider inString, []) : restPart
                                    else
                                        if semiIndex == firstDivider
                                            then -- has arg after ':'

                                                if inString !! (semiIndex + 1) == '('
                                                    then do
                                                        restPart ← argumentSplitter commandLineString (drop 2 $ dropWhile (/= ')') inString)
                                                        return $ (take firstDivider inString, takeWhile (/= ')') (drop (firstDivider + 1) inString) <> ")") : restPart
                                                    else
                                                        let nextStuff = dropWhile (/= ',') inString
                                                            remainder = if null nextStuff then [] else tail nextStuff
                                                        in  do
                                                                restPart ← argumentSplitter commandLineString remainder
                                                                return $ (take firstDivider inString, takeWhile (/= ',') (drop (firstDivider + 1) inString)) : restPart
                                            else -- arg is sub-commnd

                                                let (subCommand, remainderString) = getSubCommand inString True
                                                in  do
                                                        restPart ← argumentSplitter commandLineString remainderString
                                                        return $ (subCommand, []) : restPart


{- | freeOfSimpleErrors take command string and checks for simple for atting errors
lack of separators ',' between args
add as new errors are found
-}
freeOfSimpleErrors ∷ String → PhyG Bool
freeOfSimpleErrors commandString
    | null commandString = errorWithoutStackTrace "\n\nError in command string--empty"
    | isSequentialSubsequence ['"', '"'] commandString = do
        failWithPhase
            Parsing
            ( "\n\nCommand format error: " <> commandString <> "\n\tPossibly missing comma ',' between arguments or extra double quotes'\"'."
            )
    | isSequentialSubsequence [',', ')'] commandString = do
        failWithPhase
            Parsing
            ("\n\nCommand format error: " <> commandString <> "\n\tPossibly terminal comma ',' after or within arguments.")
    | isSequentialSubsequence [',', '('] commandString = do failWithPhase Parsing ("\n\nCommand format error: " <> commandString <> "\n\tPossibly comma ',' before '('.")
    | isSequentialSubsequence ['(', ','] commandString = do failWithPhase Parsing ("\n\nCommand format error: " <> commandString <> "\n\tPossibly starting comma ',' before arguments.")
    | otherwise =
        let beforeDoubleQuotes = dropWhile (/= '"') commandString
        in  if null beforeDoubleQuotes
                then do return True
                else -- likely more conditions to develop
                do return True


{- | parseCommandArg takes an Instruction and arg list of Strings and returns list
of parsed srguments for that instruction
-}
parseCommandArg ∷ String → Instruction → [(String, String)] → PhyG [Argument]
parseCommandArg fullCommand localInstruction argList
    | localInstruction == Read =
        if not $ null argList
            then RIF.getReadArgs fullCommand argList
            else
                failWithPhase
                    Parsing
                    ("\n\n'Read' command error '" <> fullCommand <> "': Need to specify at least one filename in double quotes")
    | otherwise = return argList


-- | movePrealignedTCM move prealigned and tcm commands to front of argument list
movePrealignedTCM ∷ [Argument] → PhyG [Argument]
movePrealignedTCM inArgList =
    if null inArgList
        then return []
        else
            let firstPart = filter ((== "prealigned") . fst) inArgList
                secondPart = filter ((== "tcm") . fst) inArgList
                restPart = filter ((/= "tcm") . fst) $ filter ((/= "prealigned") . fst) inArgList
            in  if length secondPart > 1
                    then failWithPhase Parsing ("\n\n'Read' command error '" <> show inArgList <> "': can only specify a single tcm file")
                    else return $ firstPart <> secondPart <> restPart
