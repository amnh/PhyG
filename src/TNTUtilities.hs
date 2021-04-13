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
import qualified SymMatrix as SM


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
                        renamedDefaultCharInfo = reverse $ renameTNTChars fileName 0 (replicate numChar defaultTNTCharInfo)
                        charInfoData = getTNTCharInfo fileName numChar renamedDefaultCharInfo charInfoBlock
                        checkInfo = (length charInfoData) == numChar
                in
                if not checkInfo then error ("Character information number not equal to input character number: " ++ show numChar ++ " v " ++ (show $ length charInfoData))
                else if not hasDupTerminals then (sortedData,charInfoData)
                else errorWithoutStackTrace ("\tInput file " ++ fileName ++ " has duplicate terminals: " ++ show dupList)
                )
            
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
            --trace ((show $ length curData) ++ " " ++ (show $ length $ snd $ head curData))
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
        let newName = T.pack $ fileName ++ ":" ++ show charIndex
            firstCharInfo = head inCharInfo
            newCharInfo = firstCharInfo {name = newName}
        in 
        newCharInfo : renameTNTChars fileName (charIndex + 1) (tail inCharInfo)


-- | getTNTCharInfo numChar charInfoBlock
-- bit of  but needs to update as it goes along
getTNTCharInfo :: String -> Int -> [CharInfo] -> [T.Text] -> [CharInfo]
getTNTCharInfo fileName charNumber curCharInfo inLines =
    if null inLines then curCharInfo
    else 
        let firstLine = T.strip $ head inLines
        in
        if T.null firstLine then getTNTCharInfo fileName charNumber curCharInfo (tail inLines)
        -- hit 'proc /;' line at end
        else if T.head firstLine == 'p' then curCharInfo
        else if (T.last firstLine) /= ';' then errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " processing error--ccode/costs lines must end with semicolon ';': " ++ T.unpack firstLine)
        else -- have a valid line
            let wordList = fmap T.toLower $ T.words $ T.init firstLine
                command2 = T.take 2 $ head wordList
                newCharInfo  = if command2 == (T.pack "cc") then getCCodes fileName charNumber (tail wordList) curCharInfo
                               else curCharInfo
                newCharInfo' = if command2 == (T.pack "co") then getCosts fileName charNumber (tail wordList) newCharInfo
                               else newCharInfo

            in
            if  (command2 /= (T.pack "cc")) && (command2 /= (T.pack "co")) then 
                 trace ("\n\nWarning: TNT input file " ++ fileName ++ " unrecognized/not implemented command ignored : " ++ T.unpack firstLine)
                 getTNTCharInfo fileName charNumber curCharInfo (tail inLines)
            else getTNTCharInfo fileName charNumber newCharInfo' (tail inLines)

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
        let charStatus =  head commandWordList
            scopeList = tail commandWordList
            charIndices = nub $ sort $ concat $ fmap (scopeToIndex charNumber) scopeList
            updatedCharInfo = newCharInfo curCharInfo charStatus charIndices 0 []
        in
        --trace (show charStatus ++ " " ++ (show scopeList) ++ " " ++ show charIndices)
        updatedCharInfo

-- | getCosts takes aline form TNT and modifies charac ters according to cc-code option
-- command format : costs A.B = X/Y Z U>V Q;
-- assumes X/Y and U>V have no sapces (= one word)
getCosts :: String -> Int -> [T.Text] -> [CharInfo] -> [CharInfo]
getCosts fileName charNumber commandWordList curCharInfo =
    --trace ("Costs " ++ show commandWordList) (
    if null curCharInfo then []
    else 
        let scopeList =  takeWhile (/= (T.pack "=")) commandWordList
            charIndices = nub $ sort $ concat $ fmap (scopeToIndex charNumber) scopeList
            (localAlphabet, localMatrix) = processCostsLine fileName $  tail $ dropWhile (/= (T.pack "=")) commandWordList
            updatedCharInfo = newCharInfoMatrix curCharInfo localAlphabet localMatrix charIndices 0 [] 
        in
        --trace ("Alph " ++ show localAlphabet)
        updatedCharInfo
        --)

-- | processCostsLine takes the transformatino commands of TNT and creates a TCM matrix from that
-- does not check for alphabet size or order so sorts states to get matrix
-- TNT states (alphabet elements) must be single characters
processCostsLine :: String -> [T.Text] -> ([ST.ShortText],[[Int]])
processCostsLine fileName wordList =
    if null wordList then errorWithoutStackTrace ("\n\nTNT ccode processing error:  'costs' command without transfomation costs specified")
    else 
        let localAlphabet = sort $ nub $ concat $ fmap getAlphabetElements wordList
            transCosts = getTransformationCosts localAlphabet wordList
            localMatrix = makeMatrix localAlphabet transCosts
        in
        --trace (show localAlphabet ++ " " ++ show transCosts ++ "\n\t" ++ show localMatrix)
        (localAlphabet, localMatrix)

-- | makeMatrix takes triples and makes a square matrix with 0 diagonals if  not specified
makeMatrix :: [ST.ShortText]->  [(Int, Int, Int)] -> [[Int]]
makeMatrix localAlphabet transCosts =
    if null transCosts then []
    else 
        let initialMatrix = replicate (length localAlphabet) $ replicate (length localAlphabet) 0
            newMatrix = SM.toFullLists $ SM.updateMatrix  (SM.fromLists initialMatrix) transCosts
        in 
        -- check for uninitialized, non-diagnonal cells--maybe metricity as warning
        newMatrix

-- | getTransformationCosts get state pairs and their costs
-- retuns as (i,j,k) = i->j costs k
getTransformationCosts :: [ST.ShortText]-> [T.Text] -> [(Int, Int, Int)]
getTransformationCosts localAlphabet wordList =
    if null wordList then []
    else if length wordList == 1 then errorWithoutStackTrace ("\n\nTNT ccode processing error:  'costs' command imporperly formated (transformation and costs in pairs)")
    else 
        let [transText, costText] = take 2 wordList
            transCost = readMaybe (T.unpack costText) :: Maybe Int
            directedOperator = T.find (== '>') transText
            symetricalOperator = T.find (== '/') transText
        in
        if (directedOperator == Nothing) && (symetricalOperator == Nothing) then errorWithoutStackTrace ("\n\nTNT ccode processing error:  'costs' command requires '/' or '>'")
        else 
            let fromStateText =  if symetricalOperator == Nothing then T.takeWhile (/= '>') transText
                                 else (T.takeWhile (/= '/') transText)
                toStateText   =  if symetricalOperator == Nothing then T.tail $ T.dropWhile (/= '>') transText
                                 else T.tail $ T.dropWhile (/= '/') transText
                fromState = readMaybe (T.unpack fromStateText) :: Maybe Int 
                toState   = readMaybe (T.unpack toStateText)   :: Maybe Int
            in
            if transCost == Nothing then errorWithoutStackTrace ("\n\nTNT ccode processing error:  'costs' command transformation " ++ (T.unpack costText) ++ " does not appear to be an integer.")
            else if fromState == Nothing then errorWithoutStackTrace ("\n\nTNT ccode processing error:  'costs' command transformation " ++ (T.unpack fromStateText) ++ " does not appear to be an integer.")
            else if toState == Nothing then errorWithoutStackTrace ("\n\nTNT ccode processing error:  'costs' command transformation " ++ (T.unpack toStateText) ++ " does not appear to be an integer.")
            else 
                let newTripleList = if symetricalOperator == Nothing then [(fromJust fromState, fromJust toState, fromJust transCost)] 
                                    else  [(fromJust fromState, fromJust toState, fromJust transCost), (fromJust toState, fromJust fromState, fromJust transCost)] 
                in
                newTripleList ++ getTransformationCosts localAlphabet (drop 2 wordList)

-- | getAlphabetElements takes  Text and returns non '/' '>' elements
getAlphabetElements :: T.Text -> [ST.ShortText]
getAlphabetElements inText =
    if T.null inText then []
    else 
        let hasForwardSlash = T.find (== '/') inText
            hasGreaterThan = T.find (== '>') inText
        in
        if (hasForwardSlash == Nothing) && (hasGreaterThan == Nothing) then []
        else 
            let symbolList = T.filter (/= '/') $ T.filter (/= '>') inText
            in
            fmap ST.singleton $ T.unpack symbolList
        


-- | scopeToIndex takes the number of characters and converts to a list of indices
scopeToIndex :: Int -> T.Text -> [Int]
scopeToIndex numChars scopeText =
    if T.null scopeText then []
    else    -- stop will include '.' if present`
        let (start, stop) = T.breakOn (T.pack ".") scopeText
        in
        --trace (show (start, stop)) (
        --single integer index
        if (not $ T.null  start) && (T.head start `elem` ccodeChars) then errorWithoutStackTrace ("\n\nTNT ccode processing error: ccode '" ++ (T.unpack scopeText) ++ "' incorrect format.  Scope must follow operator ")
        else if start == scopeText then 
            let scopeSingleton = readMaybe (T.unpack scopeText) :: Maybe Int
            in
            if scopeSingleton == Nothing then errorWithoutStackTrace ("\n\nTNT ccode processing error: ccode '" ++ (T.unpack scopeText) ++ "' contains non-integer")
            else if (fromJust scopeSingleton) < numChars then [fromJust scopeSingleton] 
            else errorWithoutStackTrace ("\n\nTNT ccode processing error: scope greater than char number " ++ show scopeSingleton ++ " > " ++ (show (numChars - 1)))
        else 
            let startIndex = if T.null start then Just 0
                             else readMaybe (T.unpack start) :: Maybe Int
                stopIndex = if stop == T.pack "." then Just (numChars - 1)
                            else readMaybe (T.unpack $ T.tail stop) :: Maybe Int
            in
            if (startIndex == Nothing) || (stopIndex == Nothing) then errorWithoutStackTrace ("\n\nTNT ccode processing error: ccode '" ++ (T.unpack scopeText) ++ "' contains non-integer")
            else [(max 0 (fromJust startIndex))..(min (fromJust stopIndex) numChars)]
            --)

-- | newCharInfo updates the a specific list character element
-- if that char is not in index list it is unaffected and added back to create the new list
-- in a single pass. 
-- if nothing to do (and nothing done so curCharLIst == []) then return original charInfo
-- othewise retiurn rteh reverse since new values are prepended 
newCharInfo :: [CharInfo] -> T.Text -> [Int] -> Int -> [CharInfo] -> [CharInfo] 
newCharInfo inCharList newStatus indexList charIndex curCharList =
    --trace (show charIndex ++ " " ++ show indexList ++ " " ++ (show $ length inCharList)) (
    if null inCharList then curCharList
    else if null indexList then 
        if null curCharList then inCharList
        else 
            (reverse curCharList) ++ inCharList
    else 
        let firstIndex = head indexList
            firstCharInfo =  head inCharList
        in
        if charIndex /= firstIndex then newCharInfo (tail inCharList) newStatus indexList (charIndex + 1) (firstCharInfo : curCharList) 
        else 
            let updatedCharInfo = if      newStatus == T.pack "-" then firstCharInfo {charType = NonAdd}
                                  else if newStatus == T.pack "+" then firstCharInfo {charType = Add}
                                  else if newStatus == T.pack "[" then firstCharInfo {activity = True}
                                  else if newStatus == T.pack "]" then firstCharInfo {activity = False}
                                  else if newStatus == T.pack "(" then firstCharInfo {charType = Matrix}
                                  else if newStatus == T.pack ")" then firstCharInfo {charType = NonAdd}
                                  else if newStatus == T.pack "/" then firstCharInfo {weight = 1.0}
                                  else if (T.head newStatus) ==  '/' then
                                    let newWeight = readMaybe (tail $ T.unpack newStatus) :: Maybe Double    
                                    in
                                    if newWeight == Nothing then errorWithoutStackTrace ("\n\nTNT ccode processing error: weight " ++ (tail $ T.unpack newStatus) ++ " not an integer")
                                    else firstCharInfo {weight = fromJust newWeight}
                                  else 
                                    trace ("Warning: TNT ccodes command " ++ (T.unpack newStatus) ++ " is unrecognized/not implemented--skipping")
                                    firstCharInfo
            in
            newCharInfo (tail inCharList) newStatus (tail indexList) (charIndex + 1) (updatedCharInfo : curCharList) 
            

-- | newCharInfoMatrix updates alphabet and rcm matrix for characters in indexLisyt
-- if that char is not in index list it is unaffected and added back to create the new list
-- in a single pass. 
-- if nothing to do (and nothing done so curCharLIst == []) then return original charInfo
-- othewise retiurn the reverse since new values are prepended 
newCharInfoMatrix :: [CharInfo] -> [ST.ShortText] -> [[Int]] -> [Int] -> Int -> [CharInfo] -> [CharInfo] 
newCharInfoMatrix inCharList localAlphabet localMatrix indexList charIndex curCharList =
    --trace (show charIndex ++ " " ++ show indexList ++ " " ++ (show $ length inCharList)) (
    if null inCharList then curCharList
    else if null indexList then 
        if null curCharList then inCharList
        else 
            (reverse curCharList) ++ inCharList
    else 
        let firstIndex = head indexList
            firstCharInfo =  head inCharList
        in
        if charIndex /= firstIndex then newCharInfoMatrix (tail inCharList) localAlphabet localMatrix  indexList (charIndex + 1) (firstCharInfo : curCharList) 
        else 
            let updatedCharInfo = firstCharInfo {alphabet = localAlphabet, costMatrix = localMatrix}
            in
            newCharInfoMatrix (tail inCharList) localAlphabet localMatrix  (tail indexList) (charIndex + 1) (updatedCharInfo : curCharList) 
            
