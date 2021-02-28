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
import           Control.Monad.Catch
import           Data.Char
import           Debug.Trace
import           Data.List.Split
import           Data.Maybe

import           Types


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

-- | allowedCommandList is teh permitted command string list
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
        if firstTwo == "--" then removeComments $ tail inLineList
        else 
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
        Just (fmap fromJust processedCommands)


-- | parseCommand takes a command file line and processes the String into a command and its arguemnts
parseCommand :: String -> Maybe Command
parseCommand inLine =
    if null inLine then error "Null command line"
    else 
        let instructionString = takeWhile (/= '(') inLine
            argList = splitOn "," $ init $ tail $ dropWhile (/= '(') inLine
        in
        Just (Read,[])



