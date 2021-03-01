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
    6) Add amchinery of command in general code
-}


module ProcessCommands  where

import           Control.Exception
import           Data.Typeable
import           Data.Char
import           Data.List
import           Debug.Trace
import           Data.List.Split
import           Data.Maybe

import           Types
import           GeneralUtilities


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

-- | commandList takes a String from a file and returns a list of commands and their arguments
-- these are syntactically verified, but any input files are not checked
commandList :: String -> [Command]
commandList rawContents =
    if null rawContents then errorWithoutStackTrace ("Error: Empty command file")
    else 
        let rawList = removeComments $ fmap (filter (/= ' ')) $ lines rawContents
            processedCommands = fmap parseCommand rawList
        in
        trace (show rawList)
        processedCommands

-- | allowedCommandList is the permitted command string list
allowedCommandList :: [String]
allowedCommandList = ["read", "build", "swap", "refine", "run", "report", "set", "transform"]

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
    | otherwise = 
        let errorMatch = snd $ getBestMatch (maxBound :: Int ,"no suggestion") possibleCommands inString
        in
        errorWithoutStackTrace ("\nError: Unrecognized command. By \'" ++ inString ++ "\' did you mean \'" ++ errorMatch ++ "\'?\n") 

   
-- | parseCommand takes a command file line and processes the String into a command and its arguemnts
parseCommand :: String -> Command
parseCommand inLine =
    if null inLine then error "Null command line"
    else 
        let instructionString = takeWhile (/= '(') inLine
            -- this doesn not allow recursive multi-option arguments
            -- NEED TO FIX 
            -- make in to a more sophisticated split outside of parens
            argList = argumentSplitter inLine -- splitOn "," $ init $ tail $ dropWhile (/= '(') inLine
            instruction = getInstruction instructionString allowedCommandList
            processedArg = parseCommandArg instruction argList
        in
        trace (show argList)
        (instruction, processedArg)

-- | getSubCommand takes a string ans extracts the first occurrence of the 
-- structure bleh(...), and splits the string on that, th esub command can contain 
-- parens and commas
getSubCommand :: String -> (String, String)
getSubCommand inString = 
    if null inString then ([],[])
    else 
        let firstPart = takeWhile (/= '(') inString
            secondPart = '(' : dropWhile (/= '(') inString
            parenPart = getBalancedParenPart secondPart "" 0 0
            remainderPart = drop ((length parenPart) + 1) inString -- to remove ','
        in
        (firstPart ++ parenPart, remainderPart)

-- | getBalancedParenPart stakes a string starting with '(' and takes all
-- characters until (and including) the balancing ')'
getBalancedParenPart :: String -> String -> Int -> Int -> String
getBalancedParenPart inString curString countLeft countRight =
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
argumentSplitter :: String -> [String] 
argumentSplitter inString =
    if null inString then []
    else 
        let commaIndex = if (elemIndex ',' inString) == Nothing then (maxBound :: Int) else fromJust (elemIndex ',' inString)
            semiIndex = if (elemIndex ':' inString) == Nothing then (maxBound :: Int) else fromJust (elemIndex ':' inString)
            leftParenIndex = if (elemIndex '(' inString) == Nothing then (maxBound :: Int) else fromJust (elemIndex '(' inString) 
            firstDivider = minimum [commaIndex, semiIndex, leftParenIndex]
        in
        if commaIndex == firstDivider then 
            -- no arg
            (take firstDivider inString) : argumentSplitter (drop (firstDivider + 1) inString)
        else if semiIndex == firstDivider then
            -- has arg after ':'
            if inString !! (semiIndex + 1) == '(' then 
                ((takeWhile (/= ')') inString) ++ ")") : argumentSplitter (drop 2 $ dropWhile (/= ')') inString)
            else (takeWhile (/= ',') inString) : argumentSplitter (tail $ dropWhile (/= ',') inString)
        else -- arg is sub-commnd
            let (subCommand, remainderString) = getSubCommand inString
            in
            subCommand : argumentSplitter remainderString

-- | parseCommandArg takes an Instruction and arg list of Strings and returns list
-- of parsed srguments for that instruction
parseCommandArg :: Instruction -> [String] -> [Argument]
parseCommandArg instruction argList 
    | instruction == Read = getReadArgs argList
    | otherwise = argList

-- | getReadArgs processes arguments ofr the 'read' command
-- should allow mulitp0le files and gracefully error check
getReadArgs :: [String] -> [Argument]
getReadArgs argList = 
    if null argList then []
    else 
        let firstArg = filter (/= ' ') $ head argList 
        in
        if (length firstArg) == 0 then errorWithoutStackTrace ("\n'Read' command error: Need to specify at least one filename in double quotes") 
        else 
            if (head firstArg /= '"') || (last firstArg /= '"') then errorWithoutStackTrace ("\n'Read' command error '" ++ (firstArg) ++"' : Need to specify filename in double quotes") 
            else (init $ tail firstArg) : getReadArgs (tail argList)



