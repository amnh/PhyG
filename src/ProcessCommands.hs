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

module ProcessCommands  where

import           Control.Exception
import           Data.Typeable
import           Control.Monad.Catch (throwM)
import           Data.Char
import           Debug.Trace
import           Data.List.Split
import           Data.Maybe

import           Types
import           GeneralUtilities


-- | Exception machinery for bad intial command line
data BadCommandLine = BadCommandLine
    deriving Typeable
instance Show BadCommandLine where
    show BadCommandLine = "Error: Program requires a single argument--the name of command script file.\n"
instance Exception BadCommandLine

-- | Exception machinery for empty command file
data BadCommandFile = BadCommandFile
    deriving Typeable
instance Show BadCommandFile where
    show BadCommandFile = "Error: Error in processing command file.\n"
instance Exception BadCommandFile

-- | Exception machinery for bad command
data BadCommand = BadCommand
    deriving Typeable
instance Show BadCommand where
    show BadCommand = "Error: Unrecognized command entered.\n"
instance Exception BadCommand

-- | allowedCommandList is the permitted command string list
allowedCommandList :: [String]
allowedCommandList = ["read", "build", "swap", "refine", "run", "report"]

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
commandList :: String -> Maybe [Command]
commandList rawContents =
    if null rawContents then trace ("Error: Empty command file") Nothing
    else 
        let rawList = removeComments $ fmap (filter (/= ' ')) $ lines rawContents
            processedCommands = fmap parseCommand rawList
        in
        trace (show rawList)
        Just processedCommands

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
    | otherwise = 
        let errorMatch = snd $ getBestMatch (maxBound :: Int ,"no suggestion") possibleCommands inString
        in
        trace ("\nError: Unrecognized command. By \'" ++ inString ++ "\' did you mean \'" ++ errorMatch ++ "\'?\n") NotACommand

   
-- | parseCommand takes a command file line and processes the String into a command and its arguemnts
parseCommand :: String -> Command
parseCommand inLine =
    if null inLine then error "Null command line"
    else 
        let instructionString = takeWhile (/= '(') inLine
            argList = splitOn "," $ init $ tail $ dropWhile (/= '(') inLine
            instruction = getInstruction instructionString allowedCommandList
        in
        (instruction, argList)
        


