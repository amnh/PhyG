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

module ReadInputFiles  
( executeReadCommands
 ,getReadArgs
) where

import           Types
import           Debug.Trace
import           Data.Char
import           System.IO
import           Data.List
import qualified Data.Text.Lazy  as T
import qualified Data.Text.Short as ST
import qualified LocalGraph as LG
import qualified GraphFormatUtilities as GFU
import qualified TNTUtilities as TNT
import qualified DataTransformation as DT


-- | executeReadCommands reads iput files and returns raw data 
-- need to close files after read
executeReadCommands :: [PhyloData] -> [SimpleGraph] -> [Argument] -> IO ([PhyloData], [SimpleGraph])
executeReadCommands curData curGraphs argList = do
    if null argList then return (curData, curGraphs)
    else do
        let (firstOption, firstFile) = head argList
        fileHandle <- openFile  firstFile ReadMode
        canBeReadFrom <- hIsReadable fileHandle
        if not canBeReadFrom then errorWithoutStackTrace ("\n\n'Read' error: file " ++ firstFile ++ " cannot be read")
        else hPutStrLn stderr ("Reading " ++ firstFile ++ " with option " ++ firstOption)

        -- this is awkward but need to use dot utilities
        if firstOption == "dot" then do
            dotGraph <- LG.hGetDotLocal fileHandle
            let inputDot = GFU.relabelFGL $ LG.dotToGraph dotGraph
            executeReadCommands curData (inputDot : curGraphs) (tail argList)
        -- not "dot" files
        else do
            fileContents <- hGetContents fileHandle
            -- try to figure out file type
            if null firstOption then 
                -- first cahr == 'x' then tnt
                -- first cahr '>' fasta/c if lots of spaces/length then fastc else fasta
                executeReadCommands curData curGraphs (tail argList)
            -- fasta
            else if (firstOption `elem` ["fasta", "nucleotide", "aminoacid"]) then 
                let fastaData = getFastA firstOption fileContents firstFile
                    fastaCharInfo = getFastaCharInfo fastaData firstFile firstOption
                in
                executeReadCommands ((fastaData, [fastaCharInfo]) : curData) curGraphs (tail argList)
            -- fastc
            else if (firstOption `elem` ["fastc", "custom_alphabet"])  then 
                let fastaData = getFastC firstOption fileContents firstFile
                    fastaCharInfo = getFastaCharInfo fastaData firstFile firstOption
                in
                executeReadCommands ((fastaData, [fastaCharInfo]) : curData) curGraphs (tail argList)
            -- tnt
            else if firstOption == "tnt" then
                let tntData = TNT.getTNTData fileContents firstFile
                in
                executeReadCommands (tntData : curData) curGraphs (tail argList)
            else if firstOption == "tcm" then executeReadCommands curData curGraphs (tail argList)
            else if firstOption == "prealigned" then executeReadCommands curData curGraphs (tail argList)
            -- FENEwick
            else if (firstOption `elem` ["newick" , "enewick", "fenewick"])  then 
                let thisGraphList = getFENewickGraph fileContents
                in 
                executeReadCommands curData (thisGraphList ++ curGraphs) (tail argList)
            else errorWithoutStackTrace ("\n\n'Read' command error: option " ++ firstOption ++ " not recognized/implemented")

-- | Read arg list allowable modifiers in read
readArgList :: [String]
readArgList = ["tcm", "prealigned", "nucleotide", "aminoacid", "custom_alphabet", "fasta", "fastc", "tnt", "csv", 
    "dot", "newick" , "enewick", "fenewick"]

-- | getReadArgs processes arguments ofr the 'read' command
-- should allow mulitple files and gracefully error check
-- also allows tcm file specification (limit 1 tcm per command?)
-- as fasta, fastc, tnt, tcm, prealigned
getReadArgs :: String -> [(String, String)] -> [Argument]
getReadArgs fullCommand argList = 
    if null argList then []
    else 
        let (firstPart, secondPart) = head argList 
        in
        -- plain file name with no modifier
        if (null firstPart) then
            if (head secondPart == '"') || (last secondPart == '"') then (firstPart, (init $ tail secondPart)) : getReadArgs fullCommand (tail argList)
            else errorWithoutStackTrace ("\n\n'Read' command error '" ++ (secondPart) ++"' : Need to specify filename in double quotes") 
        -- Change to allowed modifiers
        else if ((fmap toLower firstPart) `notElem` readArgList)  then errorWithoutStackTrace ("\n\n'Read' command error: " ++ fullCommand ++ " contains unrecognized option '" ++ firstPart ++ "'")
        else if (null secondPart) || (head secondPart /= '"') || (last secondPart /= '"') then errorWithoutStackTrace ("\n\n'Read' command error '" ++ (secondPart) ++"' : Need to specify filename in double quotes") 
        else (firstPart, (init $ tail secondPart)) : getReadArgs fullCommand (tail argList)

-- | getAlphabet takse a list of short-text lists and returns alphabet as list of short-text
getAlphabet :: [String] -> [ST.ShortText] -> [ST.ShortText] 
getAlphabet curList inList =
    if null inList then fmap ST.fromString $ (sort curList) `union` ["-"]
    else 
        let firstChars = fmap (:[]) $ nub $ ST.toString $ head inList 
        in
        getAlphabet (firstChars `union` curList) (tail inList)


-- | getFastaCharInfo get alphabet , names etc from processed fasta data
-- this doesn't separate ambiguities from elements--processed later
getFastaCharInfo :: [TermData] -> String -> String -> CharInfo
getFastaCharInfo inData dataName dataType = 
    if null inData then error "Empty inData in getFastaCharInfo"
    else 
        let nucleotideAlphabet = fmap ST.fromString ["A","C","G","T","U","R","Y","S","W","K","M","B","D","H","V","N","?","-"]
            aminoAcidAlphabet  = fmap ST.fromString ["A","B","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","X","Y","Z","Y", "-","?"]
            --onlyInNucleotides = [ST.fromString "U"]
            --onlyInAminoAcids = fmap ST.fromString ["E","F","I","L","P","Q","X","Z"]
            sequenceData = getAlphabet [] $ concat $ fmap snd inData
            seqType = if dataType == "nucleotide" then NucSeq
                      else if dataType == "aminoacid" then AminoSeq
                      else if dataType == "custom_alphabet" then GenSeq
                      else if (sequenceData `intersect` nucleotideAlphabet == sequenceData) then trace ("Assuming file " ++ dataName 
                          ++ " is nucleotide data. Specify `aminoacid' filetype if this is incorrect.") NucSeq
                      else if (sequenceData `intersect` aminoAcidAlphabet == sequenceData) then trace ("Assuming file " ++ dataName 
                          ++ " is amino acid data. Specify `nucleotide' filetype if this is incorrect.") AminoSeq
                      -- can fit in byte with reasonable pre-calculation costs
                      else if (length sequenceData) < 9 then SmallAlphSeq
                      else GenSeq
            defaultGenSeqCharInfo = CharInfo {
                                       charType = seqType
                                     , activity = True
                                     , weight = 1.0
                                     , costMatrix = []
                                     , name = T.pack dataName
                                     , alphabet = sequenceData
                                     }
        in
        trace (show sequenceData)
        defaultGenSeqCharInfo

-- | getFastA processes fasta file 
-- assumes single character alphabet
-- deletes '-' (unless "prealigned"), and spaces
getFastA :: String -> String -> String -> [TermData] 
getFastA modifier fileContents fileName  =
    if null fileContents then errorWithoutStackTrace ("\n\n'Read' command error: empty file")
    else if (head fileContents) /= '>' then errorWithoutStackTrace ("\n\n'Read' command error: fasta file must start with '>'")
    else 
        let terminalSplits = T.split (=='>') $ T.pack fileContents 
            pairData =  getPhyloDataPairsFastA modifier (tail terminalSplits)
            (hasDupTerminals, dupList) = DT.checkDuplicatedTerminals pairData
        in
        -- tail because initial split will an empty text
        if not hasDupTerminals then pairData
        else errorWithoutStackTrace ("\tInput file " ++ fileName ++ " has duplicate terminals: " ++ show dupList)
        
       
        

-- | getPhyloDataPairsFastA takes splits of Text and returns terminalName, Data pairs--minimal error checking
getPhyloDataPairsFastA :: String -> [T.Text] -> [TermData]
getPhyloDataPairsFastA modifier inTextList =
    if null inTextList then []
    else 
        let firstText = head inTextList
            firstName = T.takeWhile (/= '$') $ head $ T.lines firstText
            firstData = T.filter (/= ' ') $ T.toUpper $ T.concat $ tail $ T.lines firstText
            firstDataNoGaps = T.filter (/= '-') firstData
        in
        --trace (T.unpack firstName ++ "\n"  ++ T.unpack firstData) (
        if modifier == "prealigned" then (firstName, [ST.fromText $ T.toStrict firstData]) : getPhyloDataPairsFastA modifier (tail inTextList)
        else (firstName, [ST.fromText  $ T.toStrict firstDataNoGaps]) : getPhyloDataPairsFastA modifier (tail inTextList)
        --)
        
-- | getFastC processes fasta file 
-- assumes spaces between alphabet elements
-- deletes '-' (unless "prealigned")
getFastC :: String -> String -> String -> [TermData] 
getFastC modifier fileContents fileName =
    if null fileContents then errorWithoutStackTrace ("\n\n'Read' command error: empty file")
    else if (head fileContents) /= '>' then errorWithoutStackTrace ("\n\n'Read' command error: fasta file must start with '>'")
    else 
        let terminalSplits = T.split (=='>') $ T.pack fileContents 
            pairData = getPhyloDataPairsFastC modifier (tail terminalSplits)
            (hasDupTerminals, dupList) = DT.checkDuplicatedTerminals pairData
        in
        -- tail because initial split will an empty text
        if not hasDupTerminals then pairData
        else errorWithoutStackTrace ("\tInput file " ++ fileName ++ " has duplicate terminals: " ++ show dupList)

-- | getPhyloDataPairsFastA takes splits of Text and returns terminalName, Data pairs--minimal error checking
-- this splits on spaces in sequences
getPhyloDataPairsFastC :: String -> [T.Text] -> [TermData]
getPhyloDataPairsFastC modifier inTextList =
    if null inTextList then []
    else 
        let firstText = head inTextList
            firstName = T.takeWhile (/= '$') $ head $ T.lines firstText
            firstData = T.split (== ' ') $ T.concat $ tail $ T.lines firstText
            firstDataNoGaps = fmap (T.filter (/= '-')) firstData
        in
        trace (T.unpack firstName ++ "\n"  ++ (T.unpack $ T.intercalate (T.pack " ") firstData)) (
        if modifier == "prealigned" then (firstName, fmap ST.fromText  $ fmap T.toStrict firstData) : getPhyloDataPairsFastC modifier (tail inTextList)
        else (firstName, fmap ST.fromText $ fmap T.toStrict firstDataNoGaps) : getPhyloDataPairsFastC modifier (tail inTextList)
        )
        
-- | getFENewickGraph takes graph contents and returns local graph format
-- could be mulitple input graphs
getFENewickGraph :: String -> [LG.Gr T.Text Double] 
getFENewickGraph fileString = 
    -- trace (fileString)
    LG.getFENLocal (T.filter (/= '\n') $ T.strip $ T.pack fileString) 

