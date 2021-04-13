{- |
Module      :  TNTUtilities.hs
Description :  Module to read tnt input files for phylogenetic analysis
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

module TNTUtilities  (getTNTData
                      ) where

import           Types
import           Debug.Trace
import           Data.Char
import           Data.List
import           Data.Maybe
import qualified Data.Text.Lazy  as T
import qualified Data.Text.Short as ST
import qualified DataTransformation as DT
import qualified Data.Vector as V
import           Text.Read


-- getTNTData take file contents and returns raw data and char info form TNT file
getTNTData :: String -> String -> PhyloData
getTNTData inString fileName =
    if null inString then errorWithoutStackTrace ("\n\nTNT input processing error--empty file")
    else 
        let inText = T.strip $ T.pack inString
        in
        if (toLower $ T.head inText) /= 'x' then errorWithoutStackTrace ("\n\nTNT input processing error--must begin with 'xread'")
        else 
            -- look for quoted message
            let singleQuotes = T.count (T.pack "'") inText
                quotedMessage = if singleQuotes == 0 then T.pack "No TNT title message" 
                              else if singleQuotes > 2 then errorWithoutStackTrace ("\n\nTNT input processing error--too many single quotes in title")
                              else (T.split (== '\'') inText) !! 1
                restFile = tail $ T.lines $ (T.split (== '\'') inText) !! 2
                firstLine = head restFile
                numChar = read (T.unpack $ head $ T.words firstLine) :: Int 
                numTax = read (T.unpack $ last $ T.words firstLine) :: Int
            in
            trace ("\nTNT message : " ++ (T.unpack quotedMessage) ++ " with " ++ (show numTax) ++ " taxa and " ++ (show numChar) ++ " characters") (
            let semiColonLineNumber = findIndex (== T.pack ";") restFile
            in
            if semiColonLineNumber == Nothing then  errorWithoutStackTrace ("\n\nTNT input processing error--can't find ';' to end data block")
            else 
                let dataBlock = filter ((>0).T.length) $ tail $ take (fromJust semiColonLineNumber) restFile
                    charInfoBlock = filter ((>0).T.length) $ tail $ drop (fromJust semiColonLineNumber) restFile
                    numDataLines = length dataBlock
                    (interleaveNumber, interleaveRemainder) = numDataLines `quotRem` numTax
                in
                --trace (show dataBlock ++ "\n" ++ show (interleaveNumber, interleaveRemainder)) (
                if interleaveRemainder /= 0 then errorWithoutStackTrace ("\n\nTNT input processing error--number of taxa mis-specified or interleaved format error")
                else 
                    let sortedData = glueInterleave dataBlock numTax numChar []
                        (hasDupTerminals, dupList) = DT.checkDuplicatedTerminals sortedData
                        renamedDefaultCharInfo = renameTNTChars fileName 0 (replicate numChar defaultTNTCharInfo)
                        charInfoData = getTNTCharInfo fileName 0 renamedDefaultCharInfo charInfoBlock
                        checkInfo = (length charInfoData) == numChar
                in
                if not checkInfo then error "Character information number not equal to input character number"
                else if not hasDupTerminals then (sortedData,charInfoData)
                else errorWithoutStackTrace ("\tInput file " ++ fileName ++ " has duplicate terminals: " ++ show dupList)
                )--)
            
-- | glueInterleave takes interleves lines and puts them together with name error checking based on number of taxa
glueInterleave :: [T.Text] -> Int -> Int -> [(T.Text, String)] -> [TermData]
glueInterleave lineList numTax numChars curData =
    if null lineList then
        -- recheck num taxa
        -- check chars after process due to ambiguities 
        if length curData /= numTax then error ("Error in glueInterleave: final taxon number error: " ++ show numTax ++ " vs. " ++ (show $ length curData))
        else 
            let nameList = fmap fst curData
                charShortText = (fmap ST.fromString $ fmap snd curData)
            in
            trace ((show $ length curData) ++ " " ++ (show $ length $ snd $ head curData))
            zip nameList (fmap (:[]) charShortText)

    else if length lineList < numTax then error ("Error in glueInterleave: line number error")
    else 
        let thisDataBlock = fmap T.words $ take numTax lineList
            blockNames = fmap head thisDataBlock
            blockStrings = fmap T.unpack $ fmap last thisDataBlock
            canonicalNames = if (length curData > 0) then fmap fst curData
                             else blockNames
            canonicalStrings = fmap snd curData
        in
        --trace (show blockNames ++ "\n" ++ show canonicalNames ++ "\n" ++ show blockStrings ++ "\n" ++ show canonicalStrings) (
        --check for raxon order
        if (blockNames /= canonicalNames) then errorWithoutStackTrace ("\n\nTNT input processing error: interleaved taxon order error")
        else 
            let newChars = if (length curData > 0) then zipWith (++) canonicalStrings blockStrings
                           else blockStrings
            in
            glueInterleave (drop numTax lineList) numTax numChars (zip canonicalNames newChars)
        --)

-- | defaultTNTCharInfo default values for TNT characters
defaultTNTCharInfo :: CharInfo 
defaultTNTCharInfo = CharInfo { charType = NonAdd
                                , activity = True
                                , weight = 1.0
                                , costMatrix = []
                                , name = T.empty
                                , alphabet = []
                                , prealigned = True
                                }

-- | renameTNTChars creates a unique name for each character from fileNamer:Number
renameTNTChars :: String -> Int -> [CharInfo] -> [CharInfo]
renameTNTChars fileName charIndex inCharInfo =
    if null inCharInfo then []
    else 
        let newName = T.pack $ fileName ++ show charIndex
            firstCharInfo = head inCharInfo
            newCharInfo = firstCharInfo {name = newName}
        in 
        newCharInfo : renameTNTChars fileName (charIndex + 1) (tail inCharInfo)


-- | getTNTCharInfo numChar charInfoBlock
-- bit of  but needs to update as it goes along
getTNTCharInfo :: String -> Int -> [CharInfo] -> [T.Text] -> [CharInfo]
getTNTCharInfo fileName charNumber curCharInfo inLines =
    if null inLines then []
    else 
        let firstLine = T.strip $ head inLines
        in
        if (T.last firstLine) /= ';' then errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " processing error--ccode/costs lines must end with semicolon ';': " ++ T.unpack firstLine)
        else -- have a valid line
            let wordList = fmap T.toLower $ T.words $ T.init firstLine
                command2 = T.take 2 $ head wordList
                newCharInfo  = if command2 == (T.pack "cc") then getCCodes fileName charNumber (tail wordList) curCharInfo
                               else curCharInfo
                newCharInfo' = if command2 == (T.pack "co") then getCosts fileName charNumber (tail wordList) newCharInfo
                               else newCharInfo

            in
            if  (command2 /= (T.pack "cc")) && (command2 /= (T.pack "co")) then errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " processing error--unrecognized command: " ++ T.unpack firstLine)
            else newCharInfo'

-- | ccodeChars are the TNT ccode control characters
ccodeChars :: [Char]
ccodeChars = ['+', '-', '[', ']', '(', ')', '/']

-- | costsCharacters are teh TNT costs charcaters for sankoff matrix characters
costsCharacters :: [Char]
costsCharacters = ['>','/']

-- | getCCodes takes aline form TNT and modifies charac ters according to cc-code option
-- assumes single command (ccodeChars) per line
-- could sort and num so only hit each char once--but would be n^2 then.
getCCodes :: String -> Int -> [T.Text] -> [CharInfo] -> [CharInfo] 
getCCodes fileName charNumber commandWordList curCharInfo =
    if null curCharInfo then []
    else 
        let charStatus = last commandWordList
            scopeList = init commandWordList
            charIndices = concat $ fmap (scopeToIndex charNumber) scopeList
            updatedCharInfo = newCharInfo curCharInfo charStatus charIndices 0
        in
        updatedCharInfo

-- | getCosts takes aline form TNT and modifies charac ters according to cc-code option
getCosts :: String -> Int -> [T.Text] -> [CharInfo] -> [CharInfo]
getCosts fileName charNumber commandWordList curCharInfo =
    if null curCharInfo then []
    else 
    curCharInfo

-- | scopeToIndex takes the number of characters and converts to a list of indices
scopeToIndex :: Int -> T.Text -> [Int]
scopeToIndex numChars scopeText =
    if T.null scopeText then []
    else    -- stop will include '.' if present`
        let (start, stop) = T.breakOn (T.pack ".") scopeText
        in
        --single integer index
        if start == scopeText then 
            let scopeSingleton = [read (T.unpack scopeText) :: Int]
            in
            if (head scopeSingleton) < numChars then scopeSingleton
            else errorWithoutStackTrace ("\n\nTNT ccode processing error: scope greater than char number " ++ show scopeSingleton ++ " > " ++ (show (numChars - 1)))
        else 
            let startIndex = if T.null start then Just 0
                             else readMaybe (T.unpack start) :: Maybe Int
                stopIndex = if stop == T.pack "." then Just (numChars - 1)
                            else readMaybe (T.unpack $ T.tail stop) :: Maybe Int
            in
            if (startIndex == Nothing) || (stopIndex == Nothing) then errorWithoutStackTrace ("\n\nTNT ccode processing error: ccode " ++ (T.unpack scopeText) ++ " contains non-integer")
            else [(max 0 (fromJust startIndex))..(min (fromJust stopIndex) numChars)]

-- | updateCharInfo ujpdates the a specific vector character element
-- if that char is not in index list it is unaffected and added back to create th new vector
-- in a single pass. 
newCharInfo :: [CharInfo] -> T.Text -> [Int] -> Int -> [CharInfo] 
newCharInfo inCharList newStatus indexList charIndex =
    if null indexList then []
    else 
        let firstIndex = head indexList
            firstCharInfo =  head inCharList
        in
        if charIndex /= firstIndex then firstCharInfo : (newCharInfo (tail inCharList) newStatus indexList (charIndex + 1))
        else if newStatus == T.pack "-" then 
            let updatedCharInfo = firstCharInfo {charType = NonAdd}
            in
            updatedCharInfo :  (newCharInfo (tail inCharList) newStatus (tail indexList) (charIndex + 1))
        else if newStatus == T.pack "+" then 
            let updatedCharInfo = firstCharInfo {charType = Add}
            in
            updatedCharInfo :  (newCharInfo (tail inCharList) newStatus (tail indexList) (charIndex + 1))
        else if newStatus == T.pack "[" then 
            let updatedCharInfo = firstCharInfo {activity = True}
            in
            updatedCharInfo :  (newCharInfo (tail inCharList) newStatus (tail indexList) (charIndex + 1))
        else if newStatus == T.pack "]" then 
            let updatedCharInfo = firstCharInfo {activity = False}
            in
            updatedCharInfo :  (newCharInfo (tail inCharList) newStatus (tail indexList) (charIndex + 1))
        else if newStatus == T.pack "(" then 
            let updatedCharInfo = firstCharInfo {charType = Matrix}
            in
            updatedCharInfo :  (newCharInfo (tail inCharList) newStatus (tail indexList) (charIndex + 1))
        else if newStatus == T.pack ")" then 
            let updatedCharInfo = firstCharInfo {charType = NonAdd}
            in
            updatedCharInfo :  (newCharInfo (tail inCharList) newStatus (tail indexList) (charIndex + 1))
        else if newStatus == T.pack "/" then 
            let updatedCharInfo = firstCharInfo {weight = 1.0}
            in
            updatedCharInfo :  (newCharInfo (tail inCharList) newStatus (tail indexList) (charIndex + 1))
        else if (T.head newStatus) ==  '/' then 
            let newWeight = readMaybe (tail $ T.unpack newStatus) :: Maybe Double
                updatedCharInfo =  firstCharInfo {weight = fromJust newWeight}
            in
            if newWeight == Nothing then errorWithoutStackTrace ("\n\nTNT ccode processing error: weight " ++ (tail $ T.unpack newStatus) ++ " not an integer")
            else updatedCharInfo :  (newCharInfo (tail inCharList) newStatus (tail indexList) (charIndex + 1))
        else 
            trace ("Warning: TNT ccodes command " ++ (T.unpack newStatus) ++ " is unrecognized/not implemented--skipping")
            firstCharInfo : (newCharInfo (tail inCharList) newStatus indexList (charIndex + 1))
        
        
            
