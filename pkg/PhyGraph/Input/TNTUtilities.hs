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
Functions to input TNT file reading for PhyG
    This is far from complete wrt TNT functionality
    Deals with scope and add/nonadd/sankloff
               ccode, ccosts
               Ambiguities in "dense" tnt rows (single character states, no spaces)
               Ambiguiities and multi-character states (e.g. CYTB, 1.23)
               Will limit continuous character reps to 9 sig digits
                    this to allow 2x32 bit representations ina single Word64 later

One Big thing--
    For multicharacter states--
        the parsing assumes there is more than one multicharacter character.
        Can easily add a char from single-state (with spaces)
        But the parsing relies on counting hte number of "words" in a line of data
            if 2--single character states
            if > 2 --multicharacter states.
            A singl multi-character states would be parsed as single charcter states
            NEED TO FIX
-}

module Input.TNTUtilities  (getTNTData
                      ) where

import           Data.Alphabet
import           Data.Char
import qualified Data.Char                 as C
import           Data.Foldable
import qualified Data.List                 as L
import           Data.Maybe
import           Data.MetricRepresentation
import qualified Data.Set                  as Set
import qualified Data.TCM                  as TCM
import qualified Data.Text.Lazy            as T
import qualified Data.Text.Short           as ST
import           Debug.Trace
import qualified GeneralUtilities          as GU
import qualified Input.DataTransformation  as DT
import qualified Input.FastAC              as FAC
import qualified SymMatrix                 as SM
import qualified Data.Vector               as V
import           Text.Read
import           Types.Types




-- getTNTData take file contents and returns raw data and char info form TNT file
getTNTData :: String -> String -> RawData
getTNTData inString fileName =
    if null inString then errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " processing error--empty file")
    else
        let inString' = unlines $ filter (not . null) $ GU.stripString <$> lines inString
            inText = T.strip $ T.pack inString'
        in
        if toLower (T.head inText) /= 'x' then errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " processing error--must begin with 'xread'")
        else
            -- look for quoted message
            let singleQuotes = T.count (T.pack "'") inText
                quotedMessage
                  | singleQuotes == 0 = T.pack "No TNT title message"
                  | singleQuotes > 2 = errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " processing error--too many single quotes in title")
                  | otherwise = T.split (== '\'') inText !! 1
                (firstNum, secondNum, remainderText) = removeNCharNTax $ T.split (== '\'') inText !! 2
                numCharM = readMaybe (T.unpack firstNum) :: Maybe Int
                numTaxM = readMaybe (T.unpack secondNum) :: Maybe Int
                restFile = filter ((>0). T.length) $ T.lines remainderText

                numChar = fromJust numCharM
                numTax = fromJust numTaxM
            in
            -- trace (show quotedMessage ++ " " ++ (show remainderText) ++ "\n" ++ show restFile) (
            if T.null inText then errorWithoutStackTrace ("n\nTNT input file " ++ fileName ++ " processing error--empty TNT contents")
            else if null restFile then errorWithoutStackTrace ("n\nTNT input file " ++ fileName ++ " processing error--empty TNT contents after first line")
            else if isNothing numCharM then errorWithoutStackTrace ("n\nTNT input file " ++ fileName ++ " processing error--number of characters:" ++ show (T.unpack firstNum))
            else if isNothing numTaxM then errorWithoutStackTrace ("n\nTNT input file " ++ fileName ++ " processing error--number of taxa:" ++ show (T.unpack secondNum))
            else
                trace ("\nTNT file file " ++ fileName ++ " message : " ++ T.unpack quotedMessage ++ " with " ++ show numTax ++ " taxa and " ++ show numChar ++ " characters") (
                let semiColonLineNumber = L.findIndex ((== ';').T.head) restFile -- (== T.pack ";") restFile
                in
                if isNothing semiColonLineNumber then  errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " processing error--can't find ';' to end data block" ++ show restFile)
                else
                    let dataBlock = filter ((>0).T.length) $ fmap (T.filter printOrSpace) $ take (fromJust semiColonLineNumber) restFile
                        --dataBlock = filter ((>0).T.length) $ fmap (T.filter C.isPrint) $ take (fromJust semiColonLineNumber) restFile
                        charInfoBlock = filter (/= T.pack ";") $ filter ((>0).T.length) $ tail $ drop (fromJust semiColonLineNumber) restFile
                        numDataLines = length dataBlock
                        -- (interleaveNumber, interleaveRemainder) = numDataLines `quotRem` numTax
                        (_, interleaveRemainder) = numDataLines `quotRem` numTax
                    in
                    -- trace (show dataBlock ++ "\n" ++ show (interleaveNumber, interleaveRemainder, numDataLines, numTax)) (
                    if interleaveRemainder /= 0 then errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " processing error--number of taxa mis-specified or interleaved format error ")
                    else
                        let sortedData = glueInterleave fileName dataBlock numTax numChar []
                            charNumberList = fmap (length . snd) sortedData
                            nameLengthList = zip (fmap fst sortedData) charNumberList
                            incorrectLengthList = filter ((/= numChar).snd) nameLengthList
                            (hasDupTerminals, dupList) = DT.checkDuplicatedTerminals sortedData
                            renamedDefaultCharInfo = renameTNTChars fileName 0 (replicate numChar defaultTNTCharInfo)
                            charInfoData = getTNTCharInfo fileName numChar renamedDefaultCharInfo charInfoBlock
                            checkInfo = length charInfoData == numChar
                    in
                    -- trace ("Sorted data:" ++ show sortedData) (
                    --trace ("Alph2  " ++ (show $ fmap alphabet charInfoData)) (
                    if not checkInfo then error ("Character information number not equal to input character number: " ++ show numChar ++ " v " ++ show (length charInfoData))
                    else if not $ null incorrectLengthList then errorWithoutStackTrace ("\tInput file " ++ fileName ++ " has terminals with incorrect or varying numbers of characters (should be "
                        ++ show numChar ++ "):" ++ show  incorrectLengthList)
                    else if hasDupTerminals then errorWithoutStackTrace ("\tInput file " ++ fileName ++ " has duplicate terminals: " ++ show dupList)
                    else
                        -- check non-Additive alphabet to be numbers
                        -- integerize and reweight additive chars (including in ambiguities)
                        let curNames = fmap ((T.filter (/= '"') . T.filter C.isPrint) . fst) sortedData
                            curData = fmap snd sortedData
                            (curData',charInfoData') = checkAndRecodeCharacterAlphabets fileName curData charInfoData [] []
                        in
                        -- trace (show (curNames, curData'))
                        (zip curNames curData',charInfoData')
                    ) -- ) )
                    where printOrSpace a = (C.isPrint a || C.isSpace a) && (a /= '\r')

-- | removeNCharNTax removes teh first two "words" of nachr and ntx, but leaves text  with line feeds so can use
-- lines later
removeNCharNTax :: T.Text -> (T.Text, T.Text, T.Text)
removeNCharNTax inText =
    let noLeadingSpaces = T.dropWhile (not . C.isDigit) inText
        nCharacters = T.takeWhile C.isDigit noLeadingSpaces
        remainder = T.dropWhile C.isDigit noLeadingSpaces
        noLeadingSpaces' = T.dropWhile (not . C.isDigit) remainder
        nTaxa = T.takeWhile C.isDigit noLeadingSpaces'
        remainder' = T.dropWhile C.isDigit noLeadingSpaces'
    in
    (nCharacters, nTaxa, remainder')

-- | glueInterleave takes interleves lines and puts them together with name error checking based on number of taxa
-- needs to be more robust on detecting multichar blocks--now if only a single multicahr in a block would think
--its a regular block
glueInterleave :: String -> [T.Text] -> Int -> Int -> [(T.Text, [String])] -> [TermData]
glueInterleave fileName lineList numTax numChars curData
  | null lineList =
    -- recheck num taxa
    -- check chars after process due to ambiguities
    if length curData /= numTax then error ("Error in glueInterleave: final taxon number error: " ++ show numTax ++ " vs. " ++ show (length curData))
    else
        let nameList = fmap fst curData
            charShortTextList = fmap (fmap ST.fromString . snd) curData
        in
        --trace ((show $ length curData) ++ " " ++ (show $ length $ snd $ head curData))
        zip nameList charShortTextList
  | length lineList < numTax = error "Error in glueInterleave: line number error"
  | otherwise =
    let thisDataBlock = T.words <$> take numTax lineList
        blockNames = fmap head thisDataBlock
        -- if two words then regular TNT without space, otherwise spaces between states
        blockStrings = if length (head thisDataBlock) == 2 then fmap (((collectAmbiguities fileName . fmap (:[])) . T.unpack) . last) thisDataBlock
                       else fmap ((collectMultiCharAmbiguities fileName . fmap T.unpack) . tail) thisDataBlock
        canonicalNames = if not (null curData) then fmap fst curData
                         else blockNames
        canonicalStrings = fmap snd curData
    in
    --trace (show blockNames ++ "\n" ++ show canonicalNames ++ "\n" ++ show blockStrings ++ "\n" ++ show canonicalStrings) (
    --check for raxon order
    --trace ("Words:" ++ (show $ length $ head thisDataBlock)) (
    if blockNames /= canonicalNames then errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ "processing error--interleaved taxon order error or mispecified number of taxa")
    else
        let newChars = if not (null curData) then zipWith (++) canonicalStrings blockStrings
                       else blockStrings
        in
        glueInterleave fileName (drop numTax lineList) numTax numChars (zip canonicalNames newChars)
            --)

-- | collectMultiCharAmbiguities take a list of Strings and collects TNT ambiguities [X Y]  into single Strings
-- this only for multicharacter TNT characters as opposed to collectAmbiguities
collectMultiCharAmbiguities :: String -> [String] -> [String]
collectMultiCharAmbiguities fileName inStringList =
 if null inStringList then []
    else
        let firstString = head inStringList
        in
        -- might be better to check for head and last [] for better error processing
        if ('[' `elem` firstString) && (']' `elem` firstString) then
            if '-' `elem`  firstString then
                let firstPart = takeWhile (/= '-') firstString
                    secondPart = tail $ dropWhile (/= '-') firstString
                in
                concat [firstPart, " ", secondPart] : collectMultiCharAmbiguities fileName (tail inStringList)
            else errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " processing error: ambiguity format problem no ambiguity or range '-' in : " ++ firstString)
        else if ']' `elem` firstString then errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " processing error: ambiguity format problem naked ']' in : " ++ firstString)
        else if '[' `elem` firstString then
            let firstStringList = takeWhile (']' `notElem`) inStringList
                ambiguityStringList = firstStringList ++ [inStringList !! max 0 (length firstStringList)]
            in
            --trace (show firstStringList ++ show ambiguityStringList ++ show (head $ drop (length firstStringList) inStringList))
            unwords ambiguityStringList : collectMultiCharAmbiguities fileName (drop (length ambiguityStringList) inStringList)
        else firstString : collectMultiCharAmbiguities fileName (tail inStringList)
        --)

-- | collectAmbiguities take a list of Strings and collects TNT ambiguities [XY]
-- into single Strings
-- this only for single 'block' TNT characters where states are a single character
collectAmbiguities :: String -> [String] -> [String]
collectAmbiguities fileName inStringList =
    if null inStringList then []
    else
        let firstString = head inStringList
        in
        --trace ("CA " ++ firstString) (
        if firstString == "]" then errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " processing error: ambiguity format problem naked ']' in : " ++ firstString)
        else if firstString == "[" then
            let ambiguityStringList = takeWhile (/="]") inStringList ++ ["]"]
            in
            -- trace ("CA:" ++ (concat ambiguityStringList)) --  ++ " " ++ concat (drop (length $ concat ambiguityStringList) inStringList))
            concat ambiguityStringList : collectAmbiguities fileName (drop (length $ concat ambiguityStringList) inStringList)
        else firstString : collectAmbiguities fileName (tail inStringList)
        -- )


-- | defaultTNTCharInfo default values for TNT characters
defaultTNTCharInfo :: CharInfo
defaultTNTCharInfo = CharInfo { charType = NonAdd
                                , activity = True
                                , weight = 1.0
                                , costMatrix = SM.empty
                                , name = T.empty
                                , alphabet = fromSymbols []
                                , prealigned = True
                                , slimTCM    = FAC.genDiscreteDenseOfDimension (0 :: Word)
                                , wideTCM    = snd $ metricRepresentation <$> TCM.fromRows [[0::Word]]
                                , hugeTCM    = snd $ metricRepresentation <$> TCM.fromRows [[0::Word]]
                                , origInfo   = V.singleton (T.empty, NonAdd, fromSymbols [])
                                }

-- | renameTNTChars creates a unique name for each character from fileNamer:Number
renameTNTChars :: String -> Int -> [CharInfo] -> [CharInfo]
renameTNTChars fileName charIndex inCharInfo =
    if null inCharInfo then []
    else
        let newName = T.pack $ filter (/= ' ') fileName ++ ":" ++ show charIndex
            firstCharInfo = head inCharInfo
            localCharInfo = firstCharInfo {name = newName}
        in
        localCharInfo : renameTNTChars fileName (charIndex + 1) (tail inCharInfo)


-- | getTNTCharInfo numChar charInfoBlock
-- bit of  but needs to update as it goes along
getTNTCharInfo :: String -> Int -> [CharInfo] -> [T.Text] -> [CharInfo]
getTNTCharInfo fileName charNumber curCharInfo inLines =
    if null inLines then curCharInfo
    else
        let firstLine' = T.strip $ head inLines
            multipleCommandsInLine = fmap ((T.reverse . T.cons ';') . T.reverse) (filter ((>0).T.length) $ T.strip <$> T.splitOn (T.pack ";") firstLine')
            firstLine = head multipleCommandsInLine

        in
        -- trace (show multipleCommandsInLine) (
        if T.null firstLine then getTNTCharInfo fileName charNumber curCharInfo (tail inLines)
        -- hit 'proc /;' line at end
        else if T.head firstLine == 'p' then curCharInfo
        else if T.last firstLine /= ';' then errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " processing error--ccode/costs lines must end with semicolon ';': " ++ T.unpack firstLine)
        else -- have a valid line
            let wordList = T.words $ T.init firstLine
                command2 = T.toLower $ T.take 2 $ head wordList
                localCharInfo  = if command2 == T.pack "cc" then getCCodes fileName charNumber (tail wordList) curCharInfo
                               else curCharInfo
                localCharInfo' = if command2 == T.pack "co" then getCosts fileName charNumber (tail wordList) localCharInfo
                               else localCharInfo

            in
            if  (command2 /= T.pack "cc") && (command2 /= T.pack "co") then
                 trace ("\n\nWarning: TNT input file " ++ fileName ++ " unrecognized/not implemented command ignored : " ++ T.unpack firstLine)
                 getTNTCharInfo fileName charNumber curCharInfo (tail multipleCommandsInLine ++ tail inLines)
            else getTNTCharInfo fileName charNumber localCharInfo' (tail multipleCommandsInLine ++ tail inLines)
            -- )

-- | ccodeChars are the TNT ccode control characters
ccodeChars :: [Char]
ccodeChars = ['+', '-', '[', ']', '(', ')', '/']

-- | getCCodes takes aline form TNT and modifies characters according to cc-code option
-- assumes single command (ccodeChars) per line
-- could sort and num so only hit each char once--but would be n^2 then.
-- teh sinfglton stuff for compount things like "+."
getCCodes :: String -> Int -> [T.Text] -> [CharInfo] -> [CharInfo]
getCCodes fileName charNumber commandWordList curCharInfo =
    if null curCharInfo then []
    else
        let charStatus =  if T.length (head commandWordList) == 1 then head commandWordList
                          else T.singleton $ T.head $ head commandWordList
            scopeList = if T.length (head commandWordList) == 1 then tail commandWordList
                        else T.tail (head commandWordList) : tail commandWordList
            charIndices = L.nub $ L.sort $ concatMap (scopeToIndex fileName charNumber) scopeList
            updatedCharInfo = getNewCharInfo fileName curCharInfo charStatus charIndices 0 []
        in
        --trace (show charStatus ++ " " ++ (show scopeList) ++ " " ++ show charIndices)
        --if T.length charStatus > 1 then errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " character status processing error:  option not recognized/implemented " ++ (T.unpack charStatus))
        --else
        updatedCharInfo


-- | getCosts takes a line from TNT and modifies characters according to cc-code option
-- command format : costs A.B = X/Y Z U>V Q;
-- assumes X/Y and U>V have no sapces (= one word)
getCosts :: String -> Int -> [T.Text] -> [CharInfo] -> [CharInfo]
getCosts fileName charNumber commandWordList curCharInfo =
    --trace ("Costs " ++ show commandWordList) (
    if null curCharInfo then []
    else
        let scopeList =  takeWhile (/= T.pack "=") commandWordList
            charIndices = L.nub $ L.sort $ concatMap (scopeToIndex fileName charNumber) scopeList
            (localAlphabet, localMatrix) = processCostsLine fileName $  tail $ dropWhile (/= T.pack "=") commandWordList
            updatedCharInfo = newCharInfoMatrix curCharInfo localAlphabet localMatrix charIndices 0 []
        in
        --trace ("Alph " ++ (show $ fmap alphabet updatedCharInfo))
        updatedCharInfo
        --)

-- | processCostsLine takes the transformation commands of TNT and creates a TCM matrix from that
-- does not check for alphabet size or order so sorts states to get matrix
-- TNT states (alphabet elements) must be single characters
processCostsLine :: String -> [T.Text] -> ([ST.ShortText],[[Int]])
processCostsLine fileName wordList =
    if null wordList then errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " costs processing error:  'costs' command without transfomation costs specified")
    else
        let localAlphabet = L.sort $ L.nub $ concatMap getAlphabetElements wordList
            transCosts = getTransformationCosts fileName localAlphabet wordList
            localMatrix = makeMatrix fileName localAlphabet transCosts
        in
        --trace ("TNT" ++ show localAlphabet ++ " " ++ show transCosts ++ "\n\t" ++ show localMatrix)
        (localAlphabet, localMatrix)

-- | makeMatrix takes triples and makes a square matrix with 0 diagonals if  not specified
makeMatrix :: String -> [ST.ShortText]->  [(Int, Int, Int)] -> [[Int]]
makeMatrix fileName localAlphabet transCosts
  | null transCosts = []
  | length transCosts < (length localAlphabet * length localAlphabet) - length localAlphabet =
    errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " costs processing error: 'costs' command not all pairwise (non-diagnonal) transformation costs specified")
  | length transCosts > (length localAlphabet * length localAlphabet) =
    errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " costs processing error: 'costs' command too many pairwise transformation costs specified")
  | otherwise =
    let initialMatrix = replicate (length localAlphabet) $ replicate (length localAlphabet) 0
        newMatrix = SM.toFullLists $ SM.updateMatrix  (SM.fromLists initialMatrix) transCosts
    in
    -- check for uninitialized, non-diagnonal cells--maybe metricity as warning
    newMatrix

-- | getTransformationCosts get state pairs and their costs
-- retuns as (i,j,k) = i->j costs k
-- need to make this gerneral to letter states
getTransformationCosts :: String -> [ST.ShortText]-> [T.Text] -> [(Int, Int, Int)]
getTransformationCosts fileName localAlphabet wordList
  | null wordList = []
  | length wordList == 1 = errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " ccode processing error:  'costs' command imporperly formated (transformation and costs in pairs)")
  | otherwise =
    let -- wordlist must be >= 2 due to above tests
        transText = head wordList
        costText = wordList !! 1
        transCost = readMaybe (T.unpack costText) :: Maybe Int
        directedOperator = T.find (== '>') transText
        symetricalOperator = T.find (== '/') transText
    in
    if isNothing directedOperator && isNothing symetricalOperator then errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " ccode processing error:  'costs' command requires '/' or '>': " ++ T.unpack transText)
    else
        let fromStateText =  if isNothing symetricalOperator then T.takeWhile (/= '>') transText
                             else T.takeWhile (/= '/') transText
            toStateText   =  if isNothing symetricalOperator then T.tail $ T.dropWhile (/= '>') transText
                             else T.tail $ T.dropWhile (/= '/') transText
            fromState = readMaybe (T.unpack fromStateText) :: Maybe Int
            toState   = readMaybe (T.unpack toStateText)   :: Maybe Int
        in
        if isNothing transCost then errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " ccode processing error:  'costs' command transformation " ++ T.unpack costText ++ " does not appear to be an integer.")
        else if isJust fromState && isJust toState then
            -- states are numerical
            let newTripleList = if isNothing symetricalOperator then [(fromJust fromState, fromJust toState, fromJust transCost)]
                                else  [(fromJust fromState, fromJust toState, fromJust transCost), (fromJust toState, fromJust fromState, fromJust transCost)]
            in
            newTripleList ++ getTransformationCosts fileName localAlphabet (drop 2 wordList)
        else
            -- states are characters (or multicharacters)
            let fromStateIndex = L.elemIndex (ST.fromText $ T.toStrict fromStateText) localAlphabet
                toStateIndex   = L.elemIndex (ST.fromText $ T.toStrict toStateText) localAlphabet
                newTripleList  = if isNothing symetricalOperator then [(fromJust fromStateIndex, fromJust toStateIndex, fromJust transCost)]
                                else  [(fromJust fromStateIndex, fromJust toStateIndex, fromJust transCost), (fromJust toStateIndex, fromJust fromStateIndex, fromJust transCost)]
            in
            if isNothing fromStateIndex then errorWithoutStackTrace ("\n\nTNT input file " ++ fileName ++ " ccode processing error:  'costs' command " ++ show (T.unwords wordList) ++ " transformation state " ++ T.unpack fromStateText ++ " was not found in charcater alphabet " ++ show localAlphabet)
            else if isNothing toStateIndex then errorWithoutStackTrace  ("\n\nTNT input file " ++ fileName ++ " ccode processing error:  'costs' command " ++ show (T.unwords wordList) ++ " transformation state " ++ T.unpack toStateText ++ " was not found in charcater alphabet " ++ show localAlphabet)
            else newTripleList ++ getTransformationCosts fileName localAlphabet (drop 2 wordList)


-- | getAlphabetElements takes  Text and returns non '/' '>' elements
-- this is for a single "block" as in A1>G2 or C>T
getAlphabetElements :: T.Text -> [ST.ShortText]
getAlphabetElements inText =
    if T.null inText then []
    else
        let hasForwardSlash = T.find (== '/') inText
            hasGreaterThan = T.find (== '>') inText
        in
        if isNothing hasForwardSlash && isNothing hasGreaterThan then []
        else
            let firstSymbol = T.takeWhile (`notElem` ['/','>']) inText
                secondSymbol = T.tail $ T.dropWhile (`notElem` ['/','>']) inText
            in
            [ST.fromText $ T.toStrict firstSymbol, ST.fromText  $ T.toStrict secondSymbol]
            {-
            let symbolList = T.filter (/= '/') $ T.filter (/= '>') inText
            in
            fmap ST.singleton $ T.unpack symbolList
            -}


-- | scopeToIndex takes the number of characters and converts to a list of indices
scopeToIndex :: String -> Int -> T.Text -> [Int]
scopeToIndex fileName numChars scopeText =
    if T.null scopeText then []
    else    -- stop will include '.' if present`
        let (start, stop) = T.breakOn (T.pack ".") scopeText
        in
        --trace (show (start, stop)) (
        --single integer index
        if not (T.null  start) && (T.head start `elem` ccodeChars) then errorWithoutStackTrace ("\n\nTNT file " ++ fileName ++ " ccode processing error: ccode '" ++ T.unpack scopeText ++ "' incorrect format.  Scope must follow operator ")
        else if start == scopeText then
            let scopeSingleton = readMaybe (T.unpack scopeText) :: Maybe Int
            in
            if isNothing scopeSingleton then errorWithoutStackTrace ("\n\nTNT file " ++ fileName ++ " ccode processing error: ccode '" ++ T.unpack scopeText ++ "' contains non-integer")
            else if fromJust scopeSingleton < numChars then [fromJust scopeSingleton]
            else errorWithoutStackTrace ("\n\nTNT file " ++ fileName ++ " ccode processing error: scope greater than char number " ++ show (fromJust scopeSingleton) ++ " > " ++ show (numChars - 1))
        else
            let startIndex = if T.null start then Just 0
                             else readMaybe (T.unpack start) :: Maybe Int
                stopIndex = if stop == T.pack "." then Just (numChars - 1)
                            else readMaybe (T.unpack $ T.tail stop) :: Maybe Int
            in
            if isNothing startIndex || isNothing stopIndex then errorWithoutStackTrace ("\n\nTNT file " ++ fileName ++ " ccode processing error: ccode '" ++ T.unpack scopeText ++ "' contains non-integer")
            else if fromJust startIndex >= numChars then errorWithoutStackTrace ("\n\nTNT file " ++ fileName ++ " ccode processing error: scope greater than char number " ++  show (fromJust startIndex) ++ " > " ++ show (numChars - 1))
            else if fromJust stopIndex >= numChars then errorWithoutStackTrace ("\n\nTNT file " ++ fileName ++ " ccode processing error: scope greater than char number " ++  show (fromJust stopIndex) ++ " > " ++ show (numChars - 1))
            else [(max 0 (fromJust startIndex))..(min (fromJust stopIndex) numChars)]
            --)

-- | getNewCharInfo updates the a specific list character element
-- if that char is not in index list it is unaffected and added back to create the new list
-- in a single pass.
-- if nothing to do (and nothing done so curCharLIst == []) then return original charInfo
-- othewise return the reverse since new values are prepended
getNewCharInfo :: String -> [CharInfo] -> T.Text -> [Int] -> Int -> [CharInfo] -> [CharInfo]
getNewCharInfo fileName inCharList newStatus indexList charIndex curCharList =
    --trace (show charIndex ++ " " ++ show indexList ++ " " ++ (show $ length inCharList)) (
    if null inCharList then reverse curCharList
    else if null indexList then
        if null curCharList then inCharList
        else
            reverse curCharList ++ inCharList
    else
        let firstIndex = head indexList
            firstCharInfo =  head inCharList
        in
        if charIndex /= firstIndex then getNewCharInfo fileName (tail inCharList) newStatus indexList (charIndex + 1) (firstCharInfo : curCharList)
        else
            let updatedCharInfo
                  | newStatus == T.pack "-" = firstCharInfo {charType = NonAdd}
                  | newStatus == T.pack "+" = firstCharInfo {charType = Add}
                  | newStatus == T.pack "[" = firstCharInfo {activity = True}
                  | newStatus == T.pack "]" = firstCharInfo {activity = False}
                  | newStatus == T.pack "(" = firstCharInfo {charType = Matrix}
                  | newStatus == T.pack ")" = firstCharInfo {charType = NonAdd}
                  | newStatus == T.pack "/" = firstCharInfo {weight = 1.0}
                  | T.head newStatus ==  '/' =
                     let newWeight = readMaybe (tail $ T.unpack newStatus) :: Maybe Double
                     in
                     if isNothing newWeight then errorWithoutStackTrace ("\n\nTNT file " ++ fileName ++ " ccode processing error: weight " ++ tail (T.unpack newStatus) ++ " not an integer")
                     else firstCharInfo {weight = fromJust newWeight}
                          | otherwise =
                     trace ("Warning: TNT file " ++ fileName ++ " ccodes command " ++ T.unpack newStatus ++ " is unrecognized/not implemented--skipping")
                     firstCharInfo
                    in
                    getNewCharInfo fileName (tail inCharList) newStatus (tail indexList) (charIndex + 1) (updatedCharInfo : curCharList)


-- | newCharInfoMatrix updates alphabet and tcm matrix for characters in indexList
-- if that character is not in index list it is unaffected and added back to create the new list
-- in a single pass.
-- if nothing to do (and nothing done so curCharLIst == []) then return original charInfo
-- othewise retiurn the reverse since new values are prepended
newCharInfoMatrix :: [CharInfo] -> [ST.ShortText] -> [[Int]] -> [Int] -> Int -> [CharInfo] -> [CharInfo]
newCharInfoMatrix inCharList localAlphabet localMatrix indexList charIndex curCharList =
    --trace (show charIndex ++ " " ++ show indexList ++ " " ++ (show $ length inCharList)) (
    if null inCharList then reverse curCharList
    else if null indexList then
        if null curCharList then inCharList
        else
            reverse curCharList ++ inCharList
    else
        let firstIndex = head indexList
            firstCharInfo =  head inCharList
        in
        if charIndex /= firstIndex then newCharInfoMatrix (tail inCharList) localAlphabet localMatrix indexList (charIndex + 1) (firstCharInfo : curCharList)
        else
            let updatedCharInfo = firstCharInfo {alphabet = fromSymbols localAlphabet, costMatrix = SM.fromLists localMatrix}
            in
            -- trace ("TNT2" ++ (show $ alphabet updatedCharInfo))
            newCharInfoMatrix (tail inCharList) localAlphabet localMatrix  (tail indexList) (charIndex + 1) (updatedCharInfo : curCharList)

-- | reconcileAlphabetAndCostMatrix trakes the original charcater alphabet created from the cost matrix and compares to the
-- observed states.  If observed states are a subset of the inferred, the inferred are used to replace the original
-- this could happen if a matrix is specified for arange of characters, some of which do not exhibit all the states
-- otherwsie an error is thrown since states done't agree with m,artrix specification
-- this could happen for a DNA character (ACGT-) with a martix specified of numerical values (01234)
reconcileAlphabetAndCostMatrix
  :: ( Ord s
     , Show s
     )
  => String
  -> String
  -> Alphabet s -> Alphabet s -> Alphabet s
reconcileAlphabetAndCostMatrix fileName charName observedAlphabet inferredAlphabet
    | observedAlphabet `isAlphabetSubsetOf` inferredAlphabet = inferredAlphabet
    | otherwise = errorWithoutStackTrace $ fold
        [ "Error: TNT file "
        , fileName
        , " character "
        , charName, " Observed alphabet "
        , show observedAlphabet
        , " is incompatible with matrix specification states "
        , show inferredAlphabet
        ]

isAlphabetSubsetOf :: Ord s => Alphabet s -> Alphabet s -> Bool
isAlphabetSubsetOf specialAlphabet queryAlphabet =
    let querySet   = alphabetSymbols   queryAlphabet
        specialSet = alphabetSymbols specialAlphabet
    in  querySet `Set.isSubsetOf` specialSet

-- | checkAndRecodeCharacterAlphabets take RawData and checks the data with char info.
-- verifies that states (including in ambiguity) are Numerical for additive, and checks alphabets and cost matrices
-- and assigns correct alphabet to all characters
checkAndRecodeCharacterAlphabets :: String -> [[ST.ShortText]] -> [CharInfo] -> [[ST.ShortText]] -> [CharInfo] -> ([[ST.ShortText]], [CharInfo])
checkAndRecodeCharacterAlphabets fileName inData inCharInfo newData newCharInfo
  | null inCharInfo && null newCharInfo = error "Empty inCharInfo on input in checkAndRecodeCharacterAlphabets"
  | null inData = error "Empty inData in checkAndRecodeCharacterAlphabets"
  | null inCharInfo =
    --trace (show $ L.transpose newData)
    -- (reverse $ fmap reverse newData, reverse newCharInfo)
    (L.transpose $ reverse newData, reverse newCharInfo)
  | otherwise =
    let firstColumn = fmap head inData
        firstCharInfo =  head inCharInfo
        originalAlphabet = alphabet firstCharInfo
        thisName = name firstCharInfo
        (firstAlphabet, newWeight, newColumn) = getAlphabetFromSTList fileName firstColumn firstCharInfo
        updatedCharInfo = if (Matrix == charType firstCharInfo) && (firstAlphabet /= originalAlphabet) then
                            firstCharInfo {alphabet = reconcileAlphabetAndCostMatrix fileName (T.unpack thisName) firstAlphabet originalAlphabet}
                          else firstCharInfo {alphabet = firstAlphabet, weight = newWeight}
    in
    -- checkAndRecodeCharacterAlphabets fileName (fmap tail inData) (tail inCharInfo) (prependColumn newColumn newData []) (updatedCharInfo : newCharInfo)
    checkAndRecodeCharacterAlphabets fileName (fmap tail inData) (tail inCharInfo) (newColumn : newData) (updatedCharInfo : newCharInfo)


-- | getAlphabetFromSTList take a list of ST.ShortText and returns list of unique alphabet elements,
-- recodes decimat AB.XYZ to ABXYZ and reweights by that factor 1/1000 for .XYZ 1/10 for .X etc
-- checks if char is additive for numerical alphabet
getAlphabetFromSTList :: String -> [ST.ShortText] -> CharInfo -> (Alphabet ST.ShortText, Double, [ST.ShortText])
getAlphabetFromSTList fileName inStates inCharInfo =
  if null inStates then error "Empty column data in getAlphabetFromSTList"
  else
    let thisType = charType inCharInfo
        thisWeight = weight inCharInfo
        mostDecimals = if thisType == Add then maximum $ fmap getDecimals inStates
                       else 0
        (thisAlphabet, newColumn) = getAlphWithAmbiguity fileName inStates thisType mostDecimals [] []
        newWeight = if mostDecimals > 0 then  thisWeight / (10.0 ** fromIntegral mostDecimals)
                    else thisWeight
    in
    --trace (show (thisAlphabet, newWeight, newColumn, mostDecimals))
    (fromSymbols thisAlphabet, newWeight, newColumn)


-- | getDecimals tkase a state ShortText and return number decimals--if ambiguous then the most of range
getDecimals :: ST.ShortText -> Int
getDecimals inChar =
    if ST.null inChar then 0
    else
        let inCharText = ST.toString inChar
        in
        -- trace (inCharText) (
        if '.' `notElem` inCharText then 0
        else if head inCharText /= '[' then
            -- no ambiguity
            length (dropWhile (/= '.') inCharText) - 1
        else
            -- has ambiguity (range)
            let rangeStringList = words $ filter (`notElem` ['[',']']) inCharText
                fstChar = if '.' `elem` head rangeStringList then dropWhile (/= '.') (head rangeStringList)
                          else ['.']
                sndChar = if '.' `elem` last rangeStringList then dropWhile (/= '.') (last rangeStringList)
                          else ['.']
            in
            -- trace (fstChar ++ " " ++ sndChar)
            max (length fstChar) (length sndChar) - 1
            --)


-- | getAlphWithAmbiguity take a list of ShortText with information and accumulatiors
-- For both nonadditive and additve looks for [] to denote ambiguity and splits states
--  if splits on spaces if there are spaces (within []) (ala fastc or multicharacter states)
--  else if no spaces
--    if non-additive then each symbol is split out as an alphabet element -- as in TNT
--    if is additive splits on '-' to denote range
-- rescales (integerizes later) additive characters with decimal places to an integer type rep
-- for additive charcaters if states are not nummerical then throws an error
getAlphWithAmbiguity ::  String -> [ST.ShortText] -> CharType  -> Int -> [ST.ShortText] -> [ST.ShortText] -> ([ST.ShortText], [ST.ShortText])
getAlphWithAmbiguity fileName inStates thisType mostDecimals newAlph newStates =
    if null inStates then (L.sort $ L.nub newAlph, reverse newStates)
    else
        let firstState = ST.toString $ head inStates
        in
            if thisType /= Add then
                if head firstState /= '['  then
                    if firstState == "?" then getAlphWithAmbiguity fileName (tail inStates) thisType  mostDecimals newAlph (head inStates : newStates)
                    else getAlphWithAmbiguity fileName (tail inStates) thisType  mostDecimals (head inStates : newAlph) (head inStates : newStates)
                else -- ambiguity
                    let newAmbigStates  = fmap ST.fromString $ words $ filter (`notElem` ['[',']']) firstState
                    in
                    getAlphWithAmbiguity fileName (tail inStates) thisType  mostDecimals (newAmbigStates ++ newAlph) (head inStates : newStates)

        else -- additive character
            let scaleFactor = 1.0 / (10.0 ** fromIntegral mostDecimals)
            in
            if head firstState /= '['  then
                if mostDecimals == 0 then getAlphWithAmbiguity fileName (tail inStates) thisType  mostDecimals (head inStates : newAlph) (head inStates : newStates)
                else
                    let stateNumber = readMaybe firstState :: Maybe Double
                        newStateNumber = takeWhile (/='.') $ show (fromJust stateNumber / scaleFactor)
                    in
                    if isNothing stateNumber then
                        if firstState == "?" then getAlphWithAmbiguity fileName (tail inStates) thisType  mostDecimals (ST.fromString "-1" : newAlph) (ST.fromString "-1" : newStates)
                        else errorWithoutStackTrace ("\n\nTNT file " ++ fileName ++ " ccode processing error: Additive character not a number (Int/Float) " ++ firstState)
                    else getAlphWithAmbiguity fileName (tail inStates) thisType  mostDecimals (ST.fromString newStateNumber : newAlph) (ST.fromString newStateNumber : newStates)
            else
                let gutsList  = words $ filter (`notElem` ['[',']']) firstState
                    newStateNumberList = fmap readMaybe gutsList :: [Maybe Double]
                    newStateNumberStringList = fmap (((takeWhile (/='.') . show) . (/ scaleFactor)) . fromJust) newStateNumberList
                in
                if Nothing `elem` newStateNumberList then errorWithoutStackTrace ("\n\nTNT file " ++ fileName ++ " ccode processing error: Additive character not a number (Int/Float) " ++ firstState)
                else
                    let newAmbigState =  ST.fromString $ '[' : unwords newStateNumberStringList ++ "]"
                    in
                    getAlphWithAmbiguity fileName (tail inStates) thisType  mostDecimals (fmap ST.fromString newStateNumberStringList ++ newAlph) (newAmbigState : newStates)



