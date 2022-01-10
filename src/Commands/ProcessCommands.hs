{- |
Module      :  ProcessCommands.hs
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


{-  To add commands:
    1) add new command in Types.hs
    2) add string name of command to allowedCommandList
    3) add instruction processing to getInstruction
    4) add argument processing to parseCommandArg
    5) add argument processing function
        with meaningful errors
    6) Add amchinery of command in general code
-}


module Commands.ProcessCommands  where

import           Data.Char
import           Data.Foldable
import qualified Data.List            as L
import           Data.List.Split
import           Data.Maybe
import           GeneralUtilities
import qualified Input.ReadInputFiles as RIF
import           Types.Types
import           Debug.Trace


-- | expandRunCommands takes raw coomands and if a "run" command is found it reads that file
-- and adds those commands in place
-- ensures one command per line
expandRunCommands :: [String] -> [String] -> IO [String]
expandRunCommands curLines inLines =
    --trace ("EXP " ++ (show curLines) ++ show inLines) (
    if null inLines then return $ reverse curLines
    else
        let firstLineRead = removeComments [filter (/= ' ') $ head inLines]
            (firstLine, restLine) =  if null firstLineRead then ([],[])
                                     else splitCommandLine $ head firstLineRead
        in
        --trace ("FL " ++ firstLine) (
        -- only deal with run lines
        if null firstLine then expandRunCommands curLines (tail inLines)
        else if take 3 (fmap toLower firstLine) /= "run" then expandRunCommands (firstLine : curLines) (restLine : tail inLines)
        else do -- is a "run command"
             let (_, runFileList) = head $ parseCommand firstLine
             let runFileNames = fmap (checkFileNames . snd) runFileList
             fileListContents <- mapM readFile runFileNames
             let newLines = concatMap lines fileListContents
             expandRunCommands (newLines ++ curLines)  (restLine : tail inLines)
             --)
        --)


-- | splitCommandLine takes a line with potentially multiple commands and splits
-- between the first command and all others.
splitCommandLine :: String -> (String, String)
splitCommandLine inLine =
    if null inLine then ([],[])
    else
        let firstPart = takeWhile (/= '(') inLine
            parenPart = getBalancedParenPart "" (dropWhile (/= '(') inLine) 0 0
            firstCommand = firstPart ++ parenPart
            restPart = drop (length firstCommand) inLine
        in
        (firstCommand, restPart)

-- | checkFileNames checks if forst and last element of String are double quotes and remomves them
checkFileNames :: String -> String
checkFileNames inName
  | null inName = errorWithoutStackTrace "Error: Null file name"
  | head inName /= '"' = errorWithoutStackTrace ("Error: File name must be in double quotes (b): " ++ inName)
  | last inName /= '"' = errorWithoutStackTrace ("Error: File name must be in double quotes (e): " ++ inName)
  | otherwise = init $ tail inName


-- | getCommandList takes a String from a file and returns a list of commands and their arguments
-- these are syntactically verified, but any input files are not checked
-- commands in lines one to a line
getCommandList  :: [String] -> [Command]
getCommandList  rawContents =
    if null rawContents then errorWithoutStackTrace "Error: Empty command file"
    else
        let rawList = removeComments $ fmap (filter (/= ' ')) rawContents
            -- expand for read wildcards here--cretge a new, potentially longer list
            processedCommands = concatMap parseCommand rawList
        in
        --trace (show rawList)
        processedCommands


-- | removeComments deletes anyhting on line (including line)
-- after double dash "--"
removeComments :: [String] -> [String]
removeComments inLineList =
    if null inLineList then []
    else
        let firstLine = head inLineList
            firstTwo  = take 2 firstLine
        in
        -- Comment line
        if (firstTwo == "--") || null firstLine then removeComments $ tail inLineList else (
       -- Remove commments from line to end
       let nonComment = filter isPrint $ head $ splitOn "--" firstLine
       in
       nonComment : removeComments (tail inLineList))


-- | allowedCommandList is the permitted command string list
allowedCommandList :: [String]
allowedCommandList = ["build", "fuse", "read", "reblock", "refine", "rename", "report", "run", "select", "set", "swap"]


-- | getInstruction returns the command type from an input String
-- all operations on lower case
getInstruction :: String -> [String] -> Instruction
getInstruction inString possibleCommands
    | null inString = error "Empty command String"
    | fmap toLower inString == "build"     = Build
    | fmap toLower inString == "fuse"      = Fuse
    | fmap toLower inString == "read"      = Read
    | fmap toLower inString == "reblock"   = Reblock
    | fmap toLower inString == "refine"    = Rename
    | fmap toLower inString == "rename"    = Rename
    | fmap toLower inString == "report"    = Report
    | fmap toLower inString == "run"       = Run
    | fmap toLower inString == "select"    = Select
    | fmap toLower inString == "set"       = Set
    | fmap toLower inString == "swap"      = Swap
    | otherwise =
        let errorMatch = snd $ getBestMatch (maxBound :: Int ,"no suggestion") possibleCommands inString
        in  errorWithoutStackTrace $ fold
              ["\nError: Unrecognized command. By \'", inString, "\' did you mean \'", errorMatch, "\'?\n"]


-- | parseCommand takes a command file line and processes the String into a command and its arguemnts
-- asumes single command per line
parseCommand :: String -> [Command]
parseCommand inLine =
    if null inLine then []
    else
        let (firstString, restString) =  getSubCommand inLine False
            instructionString = takeWhile (/= '(') firstString --inLine
            -- this doesn not allow recursive multi-option arguments
            -- NEED TO FIX
            -- make in to a more sophisticated split outside of parens
            argList = argumentSplitter  $ init $ tail $ dropWhile (/= '(') $ filter (/= ' ') firstString
            instruction = getInstruction instructionString allowedCommandList
            processedArg = parseCommandArg firstString instruction argList
        in
        --trace (instructionString ++ " " ++  show argList)
        (instruction, processedArg) : parseCommand restString


-- | getSubCommand takes a string ans extracts the first occurrence of the
-- structure bleh(...), and splits the string on that, th esub command can contain
-- parens and commas
getSubCommand :: String -> Bool -> (String, String)
getSubCommand inString hasComma =
    if null inString then ([],[])
    else
        let firstPart = takeWhile (/= '(') inString
            secondPart = dropWhile (/= '(') inString
            parenPart = getBalancedParenPart "" secondPart 0 0
            incrCounter = if hasComma then 1 else 0
            remainderPart = drop (length (firstPart ++ parenPart) + incrCounter) inString -- to remove ','
        in
        (firstPart ++ parenPart, remainderPart)


-- | getBalancedParenPart stakes a string starting with '(' and takes all
-- characters until (and including) the balancing ')'
-- call with  getBalancedParenPart "" inString 0 0
getBalancedParenPart :: String -> String -> Int -> Int -> String
getBalancedParenPart curString inString countLeft countRight =
    if null inString then reverse curString
    else
        let firstChar = head inString
        in
        if firstChar == '(' then getBalancedParenPart (firstChar : curString) (tail inString) (countLeft + 1) countRight
        else if firstChar == ')' then
            if countLeft == countRight + 1 then reverse (firstChar : curString)
            else getBalancedParenPart  (firstChar : curString) (tail inString) countLeft (countRight + 1)
        else getBalancedParenPart (firstChar : curString) (tail inString) countLeft countRight


-- | argumentSplitter takes argument string and returns individual strings of arguments
-- which can include null, single, multiple or sub-command arguments
-- these are eac pairs of an option string (could be null) and a subarguments String (also could be null)
argumentSplitter :: String -> [(String, String)]
argumentSplitter inString
  | null inString = []
  | (freeOfSimpleErrors inString) == False = errorWithoutStackTrace ("Error in command specification format")
  | otherwise =
    let commaIndex = if (L.elemIndex ',' inString) == Nothing then (maxBound :: Int) else fromJust (L.elemIndex ',' inString)
        semiIndex = if (L.elemIndex ':' inString) == Nothing then (maxBound :: Int) else fromJust (L.elemIndex ':' inString)
        leftParenIndex = if (L.elemIndex '(' inString) == Nothing then (maxBound :: Int) else fromJust (L.elemIndex '(' inString)
        firstDivider = minimum [commaIndex, semiIndex, leftParenIndex]
    in
    -- simple no argument arg
    if firstDivider == (maxBound :: Int) then
        if head inString == '"' then [([], inString)]
        else [(inString, [])]
    else if commaIndex == firstDivider then
        -- no arg
        if head (take firstDivider inString) == '"' then ([], (take firstDivider inString)) : argumentSplitter (drop (firstDivider + 1) inString)
        else ((take firstDivider inString), []) : argumentSplitter (drop (firstDivider + 1) inString)
    else if semiIndex == firstDivider then
        -- has arg after ':'
        if inString !! (semiIndex + 1) == '(' then
            (take firstDivider inString,(takeWhile (/= ')') (drop (firstDivider + 1) inString) ++ ")")) : argumentSplitter (drop 2 $ dropWhile (/= ')') inString)
        else
            let nextStuff =  dropWhile (/= ',') inString
                remainder = if null nextStuff then [] else tail nextStuff
            in
            (take firstDivider inString, takeWhile (/= ',') (drop (firstDivider + 1) inString)) : argumentSplitter remainder
    else -- arg is sub-commnd
        let (subCommand, remainderString) = getSubCommand inString True
        in
        (subCommand, []) : argumentSplitter remainderString


-- | freeOfSimpleErrors take command string and checks for simple for atting errors
-- lack of separators ',' between args
-- add as new errors are found
freeOfSimpleErrors :: String -> Bool
freeOfSimpleErrors commandString
  | null commandString = errorWithoutStackTrace ("\n\nError in command string--empty")
  | isSequentialSubsequence  ['"','"'] commandString = errorWithoutStackTrace ("\n\nCommand format error: " ++ commandString ++ "\n\tPossibly missing comma ',' between arguments.")
  | isSequentialSubsequence  [',',')'] commandString = errorWithoutStackTrace ("\n\nCommand format error: " ++ commandString ++ "\n\tPossibly terminal comma ',' after arguments.")
  | isSequentialSubsequence  [',',')'] commandString = errorWithoutStackTrace ("\n\nCommand format error: " ++ commandString ++ "\n\tPossibly terminal comma ',' after arguments.")
  | isSequentialSubsequence  [',','('] commandString = errorWithoutStackTrace ("\n\nCommand format error: " ++ commandString ++ "\n\tPossibly comma ',' before '('.")
  | isSequentialSubsequence  ['(',','] commandString = errorWithoutStackTrace ("\n\nCommand format error: " ++ commandString ++ "\n\tPossibly starting comma ',' before arguments.")
  | otherwise =
    let beforeDoubleQuotes = dropWhile (/= '"') commandString
    in
    if null beforeDoubleQuotes then True
        -- likely more conditions to develop
    else True


-- | parseCommandArg takes an Instruction and arg list of Strings and returns list
-- of parsed srguments for that instruction
parseCommandArg :: String -> Instruction -> [(String, String)] -> [Argument]
parseCommandArg fullCommand instruction argList
    | instruction == Read = if not $ null argList then RIF.getReadArgs fullCommand argList
                            else errorWithoutStackTrace ("\n\n'Read' command error '" ++ fullCommand ++ "': Need to specify at least one filename in double quotes")
    | otherwise = argList


-- | movePrealignedTCM move prealigned and tcm commands to front of argument list
movePrealignedTCM :: [Argument] -> [Argument]
movePrealignedTCM inArgList =
    if null inArgList then []
    else
        let firstPart = filter ((== "prealigned").fst) inArgList
            secondPart = filter ((== "tcm").fst) inArgList
            restPart = filter ((/= "tcm").fst) $ filter ((/= "prealigned").fst) inArgList
        in
        if length secondPart > 1 then errorWithoutStackTrace ("\n\n'Read' command error '" ++ show inArgList ++ "': can only specify a single tcm file")
        else firstPart ++ secondPart ++ restPart


