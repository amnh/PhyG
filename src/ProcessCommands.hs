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


module ProcessCommands  where

import           Data.Char
import           Data.List
import           Debug.Trace
import           Data.List.Split
import           Data.Maybe

import           Types
import           GeneralUtilities
import           ReadInputFiles


-- | getCommandList takes a String from a file and returns a list of commands and their arguments
-- these are syntactically verified, but any input files are not checked
getCommandList  :: String -> [Command]
getCommandList  rawContents =
    if null rawContents then errorWithoutStackTrace ("Error: Empty command file")
    else 
        let rawList = removeComments $ fmap (filter (/= ' ')) $ lines rawContents
            processedCommands = concat $ fmap parseCommand rawList
        in
        trace (show rawList)
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
        if firstTwo == "--" then removeComments $ tail inLineList
        -- Empty line
        else if null firstLine then removeComments $ tail inLineList
        else 
            -- Remove commments from line to end
            let nonComment = head $ splitOn "--" firstLine
            in
            nonComment : (removeComments $ tail inLineList)

-- | allowedCommandList is the permitted command string list
allowedCommandList :: [String]
allowedCommandList = ["read", "build", "swap", "refine", "run", "report", "set", "transform", "support"]

-- | getInstruction returns teh command type forom an input String
-- all operations on lower case
getInstruction :: String -> [String] -> Instruction 
getInstruction inString possibleCommands 
    | null inString = error "Empty command String"
    | (fmap toLower inString) == "read" = Read
    | (fmap toLower inString) == "run" = Run
    | (fmap toLower inString) == "build" = Build
    | (fmap toLower inString) == "swap" = Swap
    | (fmap toLower inString) == "refine" = Refine
    | (fmap toLower inString) == "report" = Report
    | (fmap toLower inString) == "set" = Set
    | (fmap toLower inString) == "transform" = Transform
    | (fmap toLower inString) == "support" = Support
    | otherwise = 
        let errorMatch = snd $ getBestMatch (maxBound :: Int ,"no suggestion") possibleCommands inString
        in
        errorWithoutStackTrace ("\nError: Unrecognized command. By \'" ++ inString ++ "\' did you mean \'" ++ errorMatch ++ "\'?\n") 

   
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
        trace (instructionString ++ " " ++  show argList)
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
            remainderPart = drop ((length (firstPart ++ parenPart)) + incrCounter) inString -- to remove ','
        in
        (firstPart ++ parenPart, remainderPart)

-- | getBalancedParenPart stakes a string starting with '(' and takes all
-- characters until (and including) the balancing ')'
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
argumentSplitter inString =
    if null inString then []
    else if (freeOfSimpleErrors inString) == False then errorWithoutStackTrace ("Error in command specification format")
    else 
        let commaIndex = if (elemIndex ',' inString) == Nothing then (maxBound :: Int) else fromJust (elemIndex ',' inString)
            semiIndex = if (elemIndex ':' inString) == Nothing then (maxBound :: Int) else fromJust (elemIndex ':' inString)
            leftParenIndex = if (elemIndex '(' inString) == Nothing then (maxBound :: Int) else fromJust (elemIndex '(' inString) 
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
freeOfSimpleErrors commandString =
    if null commandString then  errorWithoutStackTrace ("\n\nError in command string--empty")
    else if isSequentialSubsequence  ['"','"'] commandString then  errorWithoutStackTrace ("\n\nCommand format error: " ++ commandString ++ "\n\tPossibly missing comma ',' between arguments.")
    else if isSequentialSubsequence  [',',')'] commandString then  errorWithoutStackTrace ("\n\nCommand format error: " ++ commandString ++ "\n\tPossibly terminal comma ',' after arguments.")
    else if isSequentialSubsequence  [',',')'] commandString then  errorWithoutStackTrace ("\n\nCommand format error: " ++ commandString ++ "\n\tPossibly terminal comma ',' after arguments.")
    else if isSequentialSubsequence  [',','('] commandString then  errorWithoutStackTrace ("\n\nCommand format error: " ++ commandString ++ "\n\tPossibly comma ',' before '('.")
    else if isSequentialSubsequence  ['(',','] commandString then  errorWithoutStackTrace ("\n\nCommand format error: " ++ commandString ++ "\n\tPossibly starting comma ',' before arguments.")
    --else if (last commandString) /= ')' then errorWithoutStackTrace ("\n\nCommand format error: " ++ commandString ++ "\n\tMissing final ')'.")
    else 
        let beforeDoubleQuotes = dropWhile (/= '"') commandString
        in
        if null beforeDoubleQuotes then True
            -- likely more conditions to develop
        else True

-- | parseCommandArg takes an Instruction and arg list of Strings and returns list
-- of parsed srguments for that instruction
parseCommandArg :: String -> Instruction -> [(String, String)] -> [Argument]
parseCommandArg fullCommand instruction argList 
    | instruction == Read = if (not $ null argList) then getReadArgs fullCommand argList 
                            else errorWithoutStackTrace ("\n\n'Read' command error '" ++ fullCommand ++ "': Need to specify at least one filename in double quotes") 
    | otherwise = argList


