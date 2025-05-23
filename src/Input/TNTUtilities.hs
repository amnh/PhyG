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

{- |
Module to read TNT input files for phylogenetic analysis.
-}
module Input.TNTUtilities (
    getTNTDataText,
) where

import Control.Monad (replicateM, when)
import Control.Monad.IO.Class (MonadIO (..))
import Data.Alphabet
import Data.Char
import Data.Char qualified as C
import Data.Foldable
import Data.List qualified as L
import Data.List.NonEmpty (NonEmpty (..))
import Data.Maybe
import Data.MetricRepresentation
import Data.Set qualified as Set
import Data.TCM qualified as TCM
import Data.Text.Lazy qualified as T
import Data.Text.Short qualified as ST
import Data.Vector qualified as V
import Debug.Trace
import Input.DataTransformation qualified as DT
import Input.FastAC qualified as FAC
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import PHANE.Evaluation.Verbosity (Verbosity (..))
import SymMatrix qualified as SM
import Text.Read
import Types.Types


-- getTNTDataText take file contents and returns raw data and char info form TNT file
getTNTDataText ∷ T.Text → String → PhyG RawData
getTNTDataText inString fileName =
    if T.null inString
        then errorWithoutStackTrace ("\n\nTNT input file " <> fileName <> " processing error--empty file")
        else
            let inString' = T.unlines $ filter ((/= '&') . T.head) $ filter (not . T.null) $ fmap T.strip (T.lines inString)
                inText' = T.strip inString'

                -- this for leading command like "taxname ujsed in report tnt"
                inText =
                    if (T.head inText') /= 'x'
                        then T.unlines $ tail $ T.lines inText'
                        else inText'
            in  if toLower (T.head inText) /= 'x'
                    then errorWithoutStackTrace ("\n\nTNT input file " <> fileName <> " processing error--must begin with 'xread'")
                    else -- look for quoted message

                        let singleQuotes = T.count (T.pack "'") inText
                            quotedMessage
                                | singleQuotes == 0 = T.pack "No TNT title message"
                                | singleQuotes > 2 =
                                    errorWithoutStackTrace ("\n\nTNT input file " <> fileName <> " processing error--too many single quotes in title")
                                | otherwise = T.split (== '\'') inText !! 1
                            (firstNum, secondNum, remainderText) = if singleQuotes /= 0 then removeNCharNTax $ T.split (== '\'') inText !! 2
                                                                   else removeNCharNTax $ inText
                            numCharM = readMaybe (T.unpack firstNum) ∷ Maybe Int
                            numTaxM = readMaybe (T.unpack secondNum) ∷ Maybe Int
                            restFile = filter ((> 0) . T.length) $ T.lines remainderText

                            numChar = fromJust numCharM
                            numTax = fromJust numTaxM
                        in  -- trace (show quotedMessage <> " " <> (show remainderText) <> "\n" <> show restFile) (
                            if T.null inText
                                then errorWithoutStackTrace ("n\nTNT input file " <> fileName <> " processing error--empty TNT contents")
                                else
                                    if null restFile
                                        then errorWithoutStackTrace ("n\nTNT input file " <> fileName <> " processing error--empty TNT contents after first line")
                                        else
                                            if isNothing numCharM
                                                then
                                                    errorWithoutStackTrace
                                                        ("n\nTNT input file " <> fileName <> " processing error--number of characters:" <> show (T.unpack firstNum))
                                                else
                                                    if isNothing numTaxM
                                                        then errorWithoutStackTrace ("n\nTNT input file " <> fileName <> " processing error--number of taxa:" <> show (T.unpack secondNum))
                                                        else
                                                            let semiColonLineNumber = L.findIndex ((== ';') . T.head) restFile -- (== T.pack ";") restFile
                                                            in  if isNothing semiColonLineNumber
                                                                    then
                                                                        errorWithoutStackTrace
                                                                            ("\n\nTNT input file " <> fileName <> " processing error--can't find ';' to end data block" <> show restFile)
                                                                    else
                                                                        let dataBlock = filter ((> 0) . T.length) (T.filter printOrSpace <$> take (fromJust semiColonLineNumber) restFile)
                                                                            -- dataBlock = filter ((>0).T.length) $ fmap (T.filter C.isPrint) $ take (fromJust semiColonLineNumber) restFile
                                                                            charInfoBlock = filter (/= T.pack ";") $ filter ((> 0) . T.length) $ tail $ drop (fromJust semiColonLineNumber) restFile
                                                                            numDataLines = length dataBlock
                                                                            -- (interleaveNumber, interleaveRemainder) = numDataLines `quotRem` numTax
                                                                            (_, interleaveRemainder) = numDataLines `quotRem` numTax
                                                                        in  -- trace (show dataBlock <> "\n" <> show (interleaveNumber, interleaveRemainder, numDataLines, numTax)) (
                                                                            if interleaveRemainder /= 0
                                                                                then
                                                                                    errorWithoutStackTrace
                                                                                        ("\n\nTNT input file " <> fileName <> " processing error--number of taxa mis-specified or interleaved format error ")
                                                                                else
                                                                                    let sortedData = glueInterleave fileName dataBlock numTax numChar []
                                                                                        charNumberList = fmap (length . snd) sortedData
                                                                                        nameLengthList = zip (fmap fst sortedData) charNumberList
                                                                                        incorrectLengthList = filter ((/= numChar) . snd) nameLengthList
                                                                                        (hasDupTerminals, dupList) = DT.checkDuplicatedTerminals sortedData
                                                                                    in  do
                                                                                            renamedDefaultCharInfo ← renameTNTChars fileName 0 <$> (replicateM numChar defaultTNTCharInfo)
                                                                                            charInfoData ← getTNTCharInfo fileName numChar renamedDefaultCharInfo charInfoBlock
                                                                                            let checkInfo = length charInfoData == numChar

                                                                                            -- trace ("Sorted data:" <> show sortedData) (
                                                                                            -- trace ("Alph2  " <> (show $ fmap alphabet charInfoData)) (
                                                                                            if not checkInfo
                                                                                                then
                                                                                                    error
                                                                                                        ("Character information number not equal to input character number: " <> show numChar <> " v " <> show (length charInfoData))
                                                                                                else
                                                                                                    if not $ null incorrectLengthList
                                                                                                        then
                                                                                                            errorWithoutStackTrace
                                                                                                                ( "\tInput file "
                                                                                                                    <> fileName
                                                                                                                    <> " has terminals with incorrect or varying numbers of characters (should be "
                                                                                                                    <> show numChar
                                                                                                                    <> "):"
                                                                                                                    <> show incorrectLengthList
                                                                                                                )
                                                                                                        else
                                                                                                            if hasDupTerminals
                                                                                                                then errorWithoutStackTrace ("\tInput file " <> fileName <> " has duplicate terminals: " <> show dupList)
                                                                                                                else do
                                                                                                                    -- check non-Additive alphabet to be numbers
                                                                                                                    -- integerize and reweight additive chars (including in ambiguities)
                                                                                                                    let curNames = fmap ((T.filter (/= '"') . T.filter C.isPrint) . fst) sortedData
                                                                                                                    let curData = fmap snd sortedData
                                                                                                                    let (curData', charInfoData') = checkAndRecodeCharacterAlphabets fileName curData charInfoData [] []
                                                                                                                    logWith
                                                                                                                        LogInfo
                                                                                                                        ( "\nTNT file file "
                                                                                                                            <> fileName
                                                                                                                            <> " message : "
                                                                                                                            <> T.unpack quotedMessage
                                                                                                                            <> " with "
                                                                                                                            <> show numTax
                                                                                                                            <> " taxa and "
                                                                                                                            <> show numChar
                                                                                                                            <> " characters"
                                                                                                                            <> "\n"
                                                                                                                        )
                                                                                                                    -- logWith LogInfo ("TNTU:" <> (show (fmap length curData', length charInfoData')))
                                                                                                                    let (newData, newCharInfoL) = filterInvariant curData' charInfoData' ([], [])
                                                                                                                    -- logWith LogInfo ("TNTUF:" <> (show (fmap length newData, length newCharInfoL)))
                                                                                                                    if null newCharInfoL
                                                                                                                        then failWithPhase Parsing ("TNT file " <> fileName <> " has no variant data")
                                                                                                                        else pure $ (zip curNames newData, newCharInfoL)
    where
        -- ) )
        printOrSpace a = (C.isPrint a || C.isSpace a) && (a /= '\r')


-- | filterInvariant filters out characters that are identical in all terminals
filterInvariant
    ∷ [[ST.ShortText]] → [charInfoData] → ([[ST.ShortText]], [charInfoData]) → ([[ST.ShortText]], [charInfoData])
filterInvariant inDataLL inCharInfoL newData@(newDataLL, newCharInfoL) =
    if null inCharInfoL
        then (fmap reverse newDataLL, reverse newCharInfoL)
        else
            let firstCharData = fmap head inDataLL
                firstCharData' = filter (`notElem` [ST.pack "-", ST.pack "?"]) firstCharData
                allSameL = fmap ((== (head firstCharData'))) (tail firstCharData')
                allSame = if null firstCharData' then True
                          else foldl' (&&) True allSameL
            in  -- trace ("FI:" <> (show (allSameL, allSame))) $
                if allSame
                    then filterInvariant (fmap tail inDataLL) (tail inCharInfoL) (newDataLL, newCharInfoL)
                    else
                        if null newCharInfoL
                            then filterInvariant (fmap tail inDataLL) (tail inCharInfoL) (fmap (: []) firstCharData, head inCharInfoL : newCharInfoL)
                            else filterInvariant (fmap tail inDataLL) (tail inCharInfoL) (zipWith (:) firstCharData newDataLL, head inCharInfoL : newCharInfoL)


{- | removeNCharNTax removes the first two "words" of nchar and ntax, but leaves text  with line feeds so can use
lines later
-}
removeNCharNTax ∷ T.Text → (T.Text, T.Text, T.Text)
removeNCharNTax inText =
    let noLeadingSpaces = T.dropWhile (not . C.isDigit) inText
        nCharacters = T.takeWhile C.isDigit noLeadingSpaces
        remainder = T.dropWhile C.isDigit noLeadingSpaces
        noLeadingSpaces' = T.dropWhile (not . C.isDigit) remainder
        nTaxa = T.takeWhile C.isDigit noLeadingSpaces'
        remainder' = T.dropWhile C.isDigit noLeadingSpaces'
    in  (nCharacters, nTaxa, remainder')


{- | glueInterleave takes interleves lines and puts them together with name error checking based on number of taxa
needs to be more robust on detecting multichar blocks--now if only a single multichar in a block would think
its a regular block
-}
glueInterleave ∷ String → [T.Text] → Int → Int → [(T.Text, [String])] → [TermData]
glueInterleave fileName lineList numTax numChars curData
    | null lineList =
        -- recheck num taxa
        -- check chars after process due to ambiguities
        if length curData /= numTax
            then error ("Error in glueInterleave: final taxon number error: " <> show numTax <> " vs. " <> show (length curData))
            else
                let nameList = fmap (T.strip . fst) curData
                    charShortTextList = fmap (fmap ST.fromString . snd) curData
                in  -- trace ((show $ length curData) <> " " <> (show $ length $ snd $ head curData))
                    zip nameList charShortTextList
    | length lineList < numTax = error "Error in glueInterleave: line number error"
    | otherwise =
        let thisDataBlock = T.words <$> take numTax lineList
            blockNames = fmap head thisDataBlock
            -- if two words then regular TNT without space, otherwise spaces between states
            blockStrings =
                if length (head thisDataBlock) == 2
                    then fmap (((collectAmbiguities fileName . fmap (: [])) . T.unpack) . last) thisDataBlock
                    else fmap ((collectMultiCharAmbiguities fileName . fmap T.unpack) . tail) thisDataBlock
            canonicalNames =
                if not (null curData)
                    then fmap fst curData
                    else blockNames
            canonicalStrings = fmap snd curData
        in  -- trace ("GIL: " <> show blockNames <> "\n" <> show canonicalNames <> "\n" <> show blockStrings <> "\n" <> show canonicalStrings) (
            -- check for taxon order
            -- trace ("Words:" <> (show $ length $ head thisDataBlock)) (
            if blockNames /= canonicalNames
                then
                    errorWithoutStackTrace
                        ("\n\nTNT input file " <> fileName <> "processing error--interleaved taxon order error or mispecified number of taxa")
                else
                    let newChars =
                            if not (null curData)
                                then zipWith (<>) canonicalStrings blockStrings
                                else blockStrings
                    in  glueInterleave fileName (drop numTax lineList) numTax numChars (zip canonicalNames newChars)


-- )

{- | collectMultiCharAmbiguities take a list of Strings and collects TNT ambiguities [X Y]  into single Strings
this only for multicharacter TNT characters as opposed to collectAmbiguities
-}
collectMultiCharAmbiguities ∷ String → [String] → [String]
collectMultiCharAmbiguities fileName inStringList =
    if null inStringList
        then []
        else
            let firstString = head inStringList
            in  -- might be better to check for head and last [] for better error processing
                if ('[' `elem` firstString) && (']' `elem` firstString)
                    then
                        if '-' `elem` firstString
                            then
                                let firstPart = takeWhile (/= '-') firstString
                                    secondPart = tail $ dropWhile (/= '-') firstString
                                in  concat [firstPart, " ", secondPart] : collectMultiCharAmbiguities fileName (tail inStringList)
                            else
                                errorWithoutStackTrace
                                    ( "\n\nTNT input file "
                                        <> fileName
                                        <> " processing error: ambiguity format problem no ambiguity or range '-' in : "
                                        <> firstString
                                    )
                    else
                        if ']' `elem` firstString
                            then
                                errorWithoutStackTrace
                                    ("\n\nTNT input file " <> fileName <> " processing error: ambiguity format problem naked ']' in : " <> firstString)
                            else
                                if '[' `elem` firstString
                                    then
                                        let firstStringList = takeWhile (']' `notElem`) inStringList
                                            ambiguityStringList = firstStringList <> [inStringList !! max 0 (length firstStringList)]
                                        in  -- trace (show firstStringList <> show ambiguityStringList <> show (head $ drop (length firstStringList) inStringList))
                                            unwords ambiguityStringList : collectMultiCharAmbiguities fileName (drop (length ambiguityStringList) inStringList)
                                    else firstString : collectMultiCharAmbiguities fileName (tail inStringList)


-- )

{- | collectAmbiguities take a list of Strings and collects TNT ambiguities [XY]
into single Strings
this only for single 'block' TNT characters where states are a single character
-}
collectAmbiguities ∷ String → [String] → [String]
collectAmbiguities fileName inStringList =
    if null inStringList
        then []
        else
            let firstString = head inStringList
            in  -- trace ("CA " <> firstString) (
                if firstString == "]"
                    then
                        errorWithoutStackTrace
                            ("\n\nTNT input file " <> fileName <> " processing error: ambiguity format problem naked ']' in : " <> firstString)
                    else
                        if firstString == "["
                            then
                                let ambiguityStringList = takeWhile (/= "]") inStringList <> ["]"]
                                in  -- trace ("CA:" <> (concat ambiguityStringList)) --  <> " " <> concat (drop (length $ concat ambiguityStringList) inStringList))
                                    concat ambiguityStringList : collectAmbiguities fileName (drop (length $ concat ambiguityStringList) inStringList)
                            else firstString : collectAmbiguities fileName (tail inStringList)


-- )

-- | defaultTNTCharInfo default values for TNT characters
defaultTNTCharInfo ∷ (MonadIO m) ⇒ m CharInfo
defaultTNTCharInfo =
    let a = fromSymbols $ ST.fromString "0" :| []
        f info =
            info
                { charType = NonAdd
                , activity = True
                , weight = 1.0
                , costMatrix = SM.empty
                , name = T.empty
                , alphabet = a
                , prealigned = True
                , origInfo = V.singleton (T.empty, NonAdd, a)
                }
    in  f <$> emptyCharInfo


-- | renameTNTChars creates a unique name for each character from fileNamer:Number
renameTNTChars ∷ String → Int → [CharInfo] → [CharInfo]
renameTNTChars fileName charIndex = \case
    [] → []
    firstInfo : otherInfo →
        let newName = T.pack $ filter (/= ' ') fileName <> "#" <> show charIndex
            localInfo = firstInfo{name = newName}
        in  localInfo : renameTNTChars fileName (charIndex + 1) otherInfo


{- | getTNTCharInfo numChar charInfoBlock
bit of  but needs to update as it goes along
-}
getTNTCharInfo ∷ String → Int → [CharInfo] → [T.Text] → PhyG [CharInfo]
getTNTCharInfo fileName charNumber curCharInfo inLines =
    if null inLines
        then do
            pure curCharInfo
        else
            let firstLine' = T.strip $ head inLines
                multipleCommandsInLine = fmap ((T.reverse . T.cons ';') . T.reverse) (filter ((> 0) . T.length) $ T.strip <$> T.splitOn (T.pack ";") firstLine')
                firstLine = head multipleCommandsInLine
            in  -- trace ("GTCI:" <> (show multipleCommandsInLine)) $
                if T.null firstLine
                    then getTNTCharInfo fileName charNumber curCharInfo (tail inLines)
                    else -- hit 'proc /;' line at end

                        if T.head firstLine == 'p'
                            then do
                                pure curCharInfo
                            else
                                if T.last firstLine /= ';'
                                    then
                                        errorWithoutStackTrace
                                            ( "\n\nTNT input file " <> fileName <> " processing error--ccode/costs lines must end with semicolon ';': " <> T.unpack firstLine
                                            )
                                    else do
                                        -- have a valid line
                                        let wordList = T.words $ T.init firstLine
                                        let command2 = T.toLower $ T.take 2 $ head wordList
                                        -- localCharInfoResult ← getCCodes fileName charNumber (tail wordList) curCharInfo
                                        localCharInfo ←
                                            if command2 == T.pack "cc"
                                                then getCCodes fileName charNumber (tail wordList) curCharInfo
                                                else pure curCharInfo
                                        let localCharInfo' =
                                                if command2 == T.pack "co"
                                                    then getCosts fileName charNumber (tail wordList) localCharInfo
                                                    else localCharInfo

                                        if (command2 /= T.pack "cc") && (command2 /= T.pack "co")
                                            then do
                                                logWith
                                                    LogInfo
                                                    ("\n\nWarning: TNT input file " <> fileName <> " unrecognized/not implemented command ignored : " <> T.unpack firstLine <> "\n")
                                                getTNTCharInfo fileName charNumber curCharInfo (tail multipleCommandsInLine <> tail inLines)
                                            else getTNTCharInfo fileName charNumber localCharInfo' (tail multipleCommandsInLine <> tail inLines)


-- )

-- | ccodeChars are the TNT ccode control characters
ccodeChars ∷ [Char]
ccodeChars = ['+', '-', '[', ']', '(', ')', '/']


{- | getCCodes takes a line from TNT and modifies characters according to cc-code option
assumes single command (ccodeChars) per line
could sort and num so only hit each char once--but would be n^2 then.
the singleton stuff for compount things like "+."
-}
getCCodes ∷ String → Int → [T.Text] → [CharInfo] → PhyG [CharInfo]
getCCodes fileName charNumber commandWordList curCharInfo =
    -- trace ("getCCodes " <> show commandWordList) $
    if null curCharInfo
        then do
            pure []
        else
            let charStatus =
                    if T.length (head commandWordList) == 1
                        then head commandWordList
                        else T.singleton $ T.head $ head commandWordList
                scopeList =
                    if T.length (head commandWordList) == 1
                        then tail commandWordList
                        else -- not a weight--weight gets added to scope without special case

                            if T.head (head commandWordList) /= '/'
                                then T.tail (head commandWordList) : tail commandWordList
                                else -- a weight '/'--weight gets added to scope without special case

                                    T.unwords (tail $ T.words $ T.tail (head commandWordList)) : tail commandWordList
                charIndices = L.nub $ L.sort $ concatMap (scopeToIndex fileName charNumber) scopeList
            in  do
                    updatedCharInfo ← getNewCharInfo fileName curCharInfo charStatus (head commandWordList) charIndices 0 []
                    -- trace (show charStatus <> " " <> (show scopeList) <> " " <> show charIndices)
                    -- if T.length charStatus > 1 then errorWithoutStackTrace ("\n\nTNT input file " <> fileName <> " character status processing error:  option not recognized/implemented " <> (T.unpack charStatus))
                    -- else
                    pure updatedCharInfo


{- | getCosts takes a line from TNT and modifies characters according to cc-code option
command format : costs A.B = X/Y Z U>V Q;
assumes X/Y and U>V have no sapces (= one word)
-}
getCosts ∷ String → Int → [T.Text] → [CharInfo] → [CharInfo]
getCosts fileName charNumber commandWordList curCharInfo =
    -- trace ("getCosts " <> show commandWordList) $
    if null curCharInfo
        then []
        else
            let scopeList = takeWhile (/= T.pack "=") commandWordList
                charIndices = L.nub $ L.sort $ concatMap (scopeToIndex fileName charNumber) scopeList
                (localAlphabet, localMatrix) = processCostsLine fileName $ tail $ dropWhile (/= T.pack "=") commandWordList
                updatedCharInfo = newCharInfoMatrix curCharInfo localAlphabet localMatrix charIndices 0 []
            in  --trace
                    --("Alph " <> (show $ fmap alphabet updatedCharInfo))
                    updatedCharInfo


{- | processCostsLine takes the transformation commands of TNT and creates a TCM matrix from that
does not check for alphabet size or order so sorts states to get matrix
TNT states (alphabet elements) must be single characters
-}
processCostsLine ∷ String → [T.Text] → (NonEmpty ST.ShortText, [[Int]])
processCostsLine fileName wordList =
    if null wordList
        then
            errorWithoutStackTrace
                ("\n\nTNT input file " <> fileName <> " costs processing error:  'costs' command without transfomation costs specified")
        else case L.sort . L.nub $ foldMap getAlphabetElements wordList of
            [] → errorWithoutStackTrace ("\n\nTNT input file " <> fileName <> " No symbols found!")
            s : ss →
                let localAlphabet = s :| ss
                    transCosts = getTransformationCosts fileName localAlphabet wordList
                    localMatrix = makeMatrix fileName localAlphabet transCosts
                in  (localAlphabet, localMatrix)


-- trace ("TNT" <> show localAlphabet <> " " <> show transCosts <> "\n\t" <> show localMatrix)

-- | makeMatrix takes triples and makes a square matrix with 0 diagonals if  not specified
makeMatrix ∷ String → NonEmpty ST.ShortText → [(Int, Int, Int)] → [[Int]]
makeMatrix fileName localAlphabet transCosts
    | null transCosts = []
    | length transCosts < (length localAlphabet * length localAlphabet) - length localAlphabet =
        errorWithoutStackTrace
            ( "\n\nTNT input file "
                <> fileName
                <> " costs processing error: 'costs' command not all pairwise (non-diagnonal) transformation costs specified"
            )
    | length transCosts > (length localAlphabet * length localAlphabet) =
        errorWithoutStackTrace
            ( "\n\nTNT input file " <> fileName <> " costs processing error: 'costs' command too many pairwise transformation costs specified"
            )
    | otherwise =
        let initialMatrix = replicate (length localAlphabet) $ replicate (length localAlphabet) 0
            newMatrix = SM.toFullLists $ SM.updateMatrix (SM.fromLists initialMatrix) transCosts
        in  -- check for uninitialized, non-diagnonal cells--maybe metricity as warning
            newMatrix


{- | getTransformationCosts get state pairs and their costs
retuns as (i,j,k) = i->j costs k
need to make this gerneral to letter states
-}
getTransformationCosts ∷ String → NonEmpty ST.ShortText → [T.Text] → [(Int, Int, Int)]
getTransformationCosts fileName localAlphabet wordList
    | null wordList = []
    | length wordList == 1 =
        errorWithoutStackTrace
            ( "\n\nTNT input file "
                <> fileName
                <> " ccode processing error:  'costs' command imporperly formated (transformation and costs in pairs)"
            )
    | otherwise =
        let -- wordlist must be >= 2 due to above tests
            transText = head wordList
            costText = wordList !! 1
            transCost = readMaybe (T.unpack costText) ∷ Maybe Int
            directedOperator = T.find (== '>') transText
            symetricalOperator = T.find (== '/') transText
        in  if isNothing directedOperator && isNothing symetricalOperator
                then
                    errorWithoutStackTrace
                        ("\n\nTNT input file " <> fileName <> " ccode processing error:  'costs' command requires '/' or '>': " <> T.unpack transText)
                else
                    let fromStateText =
                            if isNothing symetricalOperator
                                then T.takeWhile (/= '>') transText
                                else T.takeWhile (/= '/') transText
                        toStateText =
                            if isNothing symetricalOperator
                                then T.tail $ T.dropWhile (/= '>') transText
                                else T.tail $ T.dropWhile (/= '/') transText
                        fromState = readMaybe (T.unpack fromStateText) ∷ Maybe Int
                        toState = readMaybe (T.unpack toStateText) ∷ Maybe Int
                    in  if isNothing transCost
                            then
                                errorWithoutStackTrace
                                    ( "\n\nTNT input file "
                                        <> fileName
                                        <> " ccode processing error:  'costs' command transformation "
                                        <> T.unpack costText
                                        <> " does not appear to be an integer."
                                    )
                            else
                                if isJust fromState && isJust toState
                                    then -- states are numerical

                                        let newTripleList =
                                                if isNothing symetricalOperator
                                                    then [(fromJust fromState, fromJust toState, fromJust transCost)]
                                                    else [(fromJust fromState, fromJust toState, fromJust transCost), (fromJust toState, fromJust fromState, fromJust transCost)]
                                        in  newTripleList <> getTransformationCosts fileName localAlphabet (drop 2 wordList)
                                    else -- states are characters (or multicharacters)

                                        let fromStateIndex = L.elemIndex (ST.fromText $ T.toStrict fromStateText) $ toList localAlphabet
                                            toStateIndex = L.elemIndex (ST.fromText $ T.toStrict toStateText) $ toList localAlphabet
                                            newTripleList =
                                                if isNothing symetricalOperator
                                                    then [(fromJust fromStateIndex, fromJust toStateIndex, fromJust transCost)]
                                                    else
                                                        [ (fromJust fromStateIndex, fromJust toStateIndex, fromJust transCost)
                                                        , (fromJust toStateIndex, fromJust fromStateIndex, fromJust transCost)
                                                        ]
                                        in  if isNothing fromStateIndex
                                                then
                                                    errorWithoutStackTrace
                                                        ( "\n\nTNT input file "
                                                            <> fileName
                                                            <> " ccode processing error:  'costs' command "
                                                            <> show (T.unwords wordList)
                                                            <> " transformation state "
                                                            <> T.unpack fromStateText
                                                            <> " was not found in charcater alphabet "
                                                            <> show localAlphabet
                                                        )
                                                else
                                                    if isNothing toStateIndex
                                                        then
                                                            errorWithoutStackTrace
                                                                ( "\n\nTNT input file "
                                                                    <> fileName
                                                                    <> " ccode processing error:  'costs' command "
                                                                    <> show (T.unwords wordList)
                                                                    <> " transformation state "
                                                                    <> T.unpack toStateText
                                                                    <> " was not found in charcater alphabet "
                                                                    <> show localAlphabet
                                                                )
                                                        else newTripleList <> getTransformationCosts fileName localAlphabet (drop 2 wordList)


{- | getAlphabetElements takes  Text and returns non '/' '>' elements
this is for a single "block" as in A1>G2 or C>T
-}
getAlphabetElements ∷ T.Text → [ST.ShortText]
getAlphabetElements inText =
    if T.null inText
        then []
        else
            let hasForwardSlash = T.find (== '/') inText
                hasGreaterThan = T.find (== '>') inText
            in  if isNothing hasForwardSlash && isNothing hasGreaterThan
                    then []
                    else
                        let firstSymbol = T.takeWhile (`notElem` ['/', '>']) inText
                            secondSymbol = T.tail $ T.dropWhile (`notElem` ['/', '>']) inText
                        in  [ST.fromText $ T.toStrict firstSymbol, ST.fromText $ T.toStrict secondSymbol]


{-
let symbolList = T.filter (/= '/') $ T.filter (/= '>') inText
in
fmap ST.singleton $ T.unpack symbolList
-}

-- | scopeToIndex takes the number of characters and converts to a list of indices
scopeToIndex ∷ String → Int → T.Text → [Int]
scopeToIndex fileName numChars scopeText =
    if T.null scopeText
        then []
        else -- stop will include '.' if present`

            let (start, stop) = T.breakOn (T.pack ".") scopeText
            in  -- trace (show (start, stop)) (
                -- single integer index
                if not (T.null start) && (T.head start `elem` ccodeChars)
                    then
                        errorWithoutStackTrace
                            ( "\n\nTNT file "
                                <> fileName
                                <> " ccode processing error: ccode '"
                                <> T.unpack scopeText
                                <> "' incorrect format.  Scope must follow operator "
                            )
                    else
                        if start == scopeText
                            then
                                let scopeSingleton = readMaybe (T.unpack scopeText) ∷ Maybe Int
                                in  if isNothing scopeSingleton
                                        then
                                            errorWithoutStackTrace
                                                ("\n\nTNT file " <> fileName <> " ccode processing error: ccode '" <> T.unpack scopeText <> "' contains non-integer (0)")
                                        else
                                            if fromJust scopeSingleton < numChars
                                                then [fromJust scopeSingleton]
                                                else
                                                    errorWithoutStackTrace
                                                        ( "\n\nTNT file "
                                                            <> fileName
                                                            <> " ccode processing error: scope greater than char number "
                                                            <> show (fromJust scopeSingleton)
                                                            <> " > "
                                                            <> show (numChars - 1)
                                                        )
                            else
                                let startIndex =
                                        if T.null start
                                            then Just 0
                                            else readMaybe (T.unpack start) ∷ Maybe Int
                                    stopIndex =
                                        if stop == T.pack "."
                                            then Just (numChars - 1)
                                            else readMaybe (T.unpack $ T.tail stop) ∷ Maybe Int
                                in  if isNothing startIndex || isNothing stopIndex
                                        then
                                            errorWithoutStackTrace
                                                ("\n\nTNT file " <> fileName <> " ccode processing error: ccode '" <> T.unpack scopeText <> "' contains non-integer (1)")
                                        else
                                            if fromJust startIndex >= numChars
                                                then
                                                    errorWithoutStackTrace
                                                        ( "\n\nTNT file "
                                                            <> fileName
                                                            <> " ccode processing error: scope greater than char number "
                                                            <> show (fromJust startIndex)
                                                            <> " > "
                                                            <> show (numChars - 1)
                                                        )
                                                else
                                                    if fromJust stopIndex >= numChars
                                                        then
                                                            errorWithoutStackTrace
                                                                ( "\n\nTNT file "
                                                                    <> fileName
                                                                    <> " ccode processing error: scope greater than char number "
                                                                    <> show (fromJust stopIndex)
                                                                    <> " > "
                                                                    <> show (numChars - 1)
                                                                )
                                                        else [(max 0 (fromJust startIndex)) .. (min (fromJust stopIndex) numChars)]


-- )

{- | getNewCharInfo updates the a specific list character element
if that char is not in index list it is unaffected and added back to create the new list
in a single pass.
if nothing to do (and nothing done so curCharLIst == []) then return original charInfo
othewise return the reverse since new values are prepended
-}
getNewCharInfo ∷ String → [CharInfo] → T.Text → T.Text → [Int] → Int → [CharInfo] → PhyG [CharInfo]
getNewCharInfo fileName inCharList newStatus newStatusFull indexList charIndex curCharList =
    -- trace (show charIndex <> " " <> show indexList <> " " <> (show $ length inCharList)) (
    if null inCharList
        then do
            pure $ reverse curCharList
        else
            if null indexList
                then
                    if null curCharList
                        then do
                            pure inCharList
                        else do
                            pure $ reverse curCharList <> inCharList
                else
                    let firstIndex = head indexList
                        firstInfo = head inCharList
                    in  if charIndex /= firstIndex
                            then getNewCharInfo fileName (tail inCharList) newStatus newStatusFull indexList (charIndex + 1) (firstInfo : curCharList)
                            else
                                let updatedCharInfo
                                        | newStatus == T.pack "-" = firstInfo{charType = NonAdd}
                                        | newStatus == T.pack "+" = firstInfo{charType = Add}
                                        | newStatus == T.pack "[" = firstInfo{activity = True}
                                        | newStatus == T.pack "]" = firstInfo{activity = False}
                                        | newStatus == T.pack "(" = firstInfo{charType = Matrix}
                                        | newStatus == T.pack ")" = firstInfo{charType = NonAdd}
                                        -- \| newStatus == T.pack "/" = firstInfo {weight = 1.0}
                                        | T.head newStatus == '/' =
                                            let newWeight = readMaybe (tail $ T.unpack newStatusFull) ∷ Maybe Double
                                            in  if isNothing newWeight
                                                    then
                                                        errorWithoutStackTrace
                                                            ("\n\nTNT file " <> fileName <> " ccode processing error: weight " <> tail (T.unpack newStatusFull) <> " not a double")
                                                    else firstInfo{weight = fromJust newWeight}
                                        | otherwise =
                                            -- trace ("Warning: TNT file " <> fileName <> " ccodes command " <> T.unpack newStatus <> " is unrecognized/not implemented--skipping")
                                            firstInfo
                                in  do
                                        when ((T.unpack newStatus `notElem` ["-", "+", "[", "]", "(", ")"]) && (T.head newStatus /= '/')) $
                                            logWith
                                                LogWarn
                                                ( "Warning: TNT file "
                                                    <> fileName
                                                    <> " ccodes command "
                                                    <> T.unpack newStatus
                                                    <> " is unrecognized/not implemented--skipping"
                                                    <> "\n"
                                                )
                                        getNewCharInfo
                                            fileName
                                            (tail inCharList)
                                            newStatus
                                            newStatusFull
                                            (tail indexList)
                                            (charIndex + 1)
                                            (updatedCharInfo : curCharList)


{- | newCharInfoMatrix updates alphabet and tcm matrix for characters in indexList
if that character is not in index list it is unaffected and added back to create the new list
in a single pass.
if nothing to do (and nothing done so curCharList == []) then return original charInfo
othewise retiurn the reverse since new values are prepended
-}
newCharInfoMatrix ∷ [CharInfo] → NonEmpty ST.ShortText → [[Int]] → [Int] → Int → [CharInfo] → [CharInfo]
newCharInfoMatrix inCharList localAlphabet localMatrix indexList charIndex curCharList =
    -- trace (show charIndex <> " " <> show indexList <> " " <> (show $ length inCharList)) (
    if null inCharList
        then reverse curCharList
        else
            if null indexList
                then
                    if null curCharList
                        then inCharList
                        else reverse curCharList <> inCharList
                else
                    let firstIndex = head indexList
                        firstInfo = head inCharList
                    in  if charIndex /= firstIndex
                            then newCharInfoMatrix (tail inCharList) localAlphabet localMatrix indexList (charIndex + 1) (firstInfo : curCharList)
                            else -- let updatedCharInfo = firstInfo {alphabet = fromSymbols localAlphabet, costMatrix = SM.fromLists localMatrix}

                                let updatedCharInfo = firstInfo{alphabet = fromSymbols localAlphabet, costMatrix = SM.fromLists localMatrix}
                                in  -- trace ("TNT2" <> (show $ alphabet updatedCharInfo))
                                    newCharInfoMatrix (tail inCharList) localAlphabet localMatrix (tail indexList) (charIndex + 1) (updatedCharInfo : curCharList)


{- | reconcileAlphabetAndCostMatrix trakes the original charcater alphabet created from the cost matrix and compares to the
observed states.  If observed states are a subset of the inferred, the inferred are used to replace the original
this could happen if a matrix is specified for arange of characters, some of which do not exhibit all the states
otherwsie an error is thrown since states don't agree with m,artrix specification
this could happen for a DNA character (ACGT-) with a martix specified of numerical values (01234)
-}
reconcileAlphabetAndCostMatrix
    ∷ ( Ord s
      , Show s
      )
    ⇒ String
    → String
    → Alphabet s
    → Alphabet s
    → Alphabet s
reconcileAlphabetAndCostMatrix fileName charName observedAlphabet inferredAlphabet
    | observedAlphabet `isAlphabetSubsetOf` inferredAlphabet = inferredAlphabet
    | otherwise =
        errorWithoutStackTrace $
            fold
                [ "Error: TNT file "
                , fileName
                , " character number "
                , tail $ dropWhile (/= '#') charName
                , " Observed "
                , show observedAlphabet
                , " is incompatible with matrix specification states "
                , show inferredAlphabet
                ]


-- checks is observed alphabet is sub set of inferred (ie that from cost matrix string)
-- this can happen if a general DNA cost is sspecified for many characters, but some
-- characters may have only a few of the states.
isAlphabetSubsetOf ∷ (Show s, Ord s) ⇒ Alphabet s → Alphabet s → Bool
isAlphabetSubsetOf specialAlphabet queryAlphabet =
    let querySet = alphabetSymbols queryAlphabet
        specialSet = alphabetSymbols specialAlphabet
    in  not $ querySet `Set.isSubsetOf` specialSet


{-
-- this seems to produce backwards logic
isAlphabetSubsetOf' :: Ord s => Alphabet s -> Alphabet s -> Bool
isAlphabetSubsetOf' specialAlphabet queryAlphabet =
    let querySet   = alphabetSymbols   queryAlphabet
        specialSet = alphabetSymbols specialAlphabet
    in  trace ("RACM: " <> (show $ querySet `Set.isSubsetOf` specialSet)) querySet `Set.isSubsetOf` specialSet
-}

{- | checkAndRecodeCharacterAlphabets take RawData and checks the data with char info.
verifies that states (including in ambiguity) are Numerical for additive, and checks alphabets and cost matrices
and assigns correct alphabet to all characters
-}
checkAndRecodeCharacterAlphabets
    ∷ String → [[ST.ShortText]] → [CharInfo] → [[ST.ShortText]] → [CharInfo] → ([[ST.ShortText]], [CharInfo])
checkAndRecodeCharacterAlphabets fileName inData inCharInfo newData newCharInfo
    | null inCharInfo && null newCharInfo = error "Empty inCharInfo on input in checkAndRecodeCharacterAlphabets"
    | null inData = error "Empty inData in checkAndRecodeCharacterAlphabets"
    | null inCharInfo =
        -- trace (show $ L.transpose newData)
        -- (reverse $ fmap reverse newData, reverse newCharInfo)
        (L.transpose $ reverse newData, reverse newCharInfo)
    | otherwise =
        let firstColumn = fmap head inData
            firstInfo = head inCharInfo
            originalAlphabet = alphabet firstInfo
            thisName = name firstInfo
            (firstAlphabet, newWeight, newColumn) = getAlphabetFromSTList fileName firstColumn firstInfo
            updatedCharInfo =
                if (Matrix == charType firstInfo) && (firstAlphabet /= originalAlphabet)
                    then firstInfo{alphabet = reconcileAlphabetAndCostMatrix fileName (T.unpack thisName) firstAlphabet originalAlphabet}
                    else firstInfo{alphabet = firstAlphabet, weight = newWeight}
        in  -- checkAndRecodeCharacterAlphabets fileName (fmap tail inData) (tail inCharInfo) (prependColumn newColumn newData []) (updatedCharInfo : newCharInfo)
            checkAndRecodeCharacterAlphabets
                fileName
                (fmap tail inData)
                (tail inCharInfo)
                (newColumn : newData)
                (updatedCharInfo : newCharInfo)


{- | getAlphabetFromSTList take a list of ST.ShortText and returns list of unique alphabet elements,
recodes decimat AB.XYZ to ABXYZ and reweights by that factor 1/1000 for .XYZ 1/10 for .X etc
checks if char is additive for numerical alphabet
-}
getAlphabetFromSTList ∷ String → [ST.ShortText] → CharInfo → (Alphabet ST.ShortText, Double, [ST.ShortText])
getAlphabetFromSTList fileName inStates inCharInfo =
    if null inStates
        then error "Empty column data in getAlphabetFromSTList"
        else
            let thisType = charType inCharInfo
                thisWeight = weight inCharInfo
                mostDecimals =
                    if thisType == Add
                        then maximum $ fmap getDecimals inStates
                        else 0
                (thisAlphabet', newColumn) = getAlphWithAmbiguity fileName inStates thisType mostDecimals [] []

                newWeight
                    | mostDecimals > 0 = thisWeight / (10.0 ** fromIntegral mostDecimals)
                    | otherwise = thisWeight

                thisAlphabet = case thisAlphabet' of
                    [] → ST.fromString "0" :| []
                    s : ss → s :| ss
            in  -- trace (show (thisAlphabet, newWeight, newColumn, mostDecimals))
                -- (fromSymbols thisAlphabet, newWeight, newColumn)
                -- trace ("getAlph weight: " <> (show thisWeight))
                (fromSymbols thisAlphabet, newWeight, newColumn)


-- | getDecimals tkase a state ShortText and return number decimals--if ambiguous then the most of range
getDecimals ∷ ST.ShortText → Int
getDecimals inChar =
    if ST.null inChar
        then 0
        else
            let inCharText = ST.toString inChar
            in  -- trace (inCharText) (
                if '.' `notElem` inCharText
                    then 0
                    else
                        if head inCharText /= '['
                            then -- no ambiguity
                                length (dropWhile (/= '.') inCharText) - 1
                            else -- has ambiguity (range)

                                let rangeStringList = words $ filter (`notElem` ['[', ']']) inCharText
                                    fstChar =
                                        if '.' `elem` head rangeStringList
                                            then dropWhile (/= '.') (head rangeStringList)
                                            else ['.']
                                    sndChar =
                                        if '.' `elem` last rangeStringList
                                            then dropWhile (/= '.') (last rangeStringList)
                                            else ['.']
                                in  -- trace (fstChar <> " " <> sndChar)
                                    max (length fstChar) (length sndChar) - 1


-- )

{- | getAlphWithAmbiguity take a list of ShortText with information and accumulatiors
For both nonadditive and additive. Searches for [] to denote ambiguity and splits states
 if splits on spaces if there are spaces (within []) (ala fastc or multicharacter states)
 else if no spaces
   if non-additive then each symbol is split out as an alphabet element -- as in TNT
   if is additive splits on '-' to denote range
rescales (integerizes later) additive characters with decimal places to an integer type rep
for additive characters if states are not nummerical then throws an error
-}
getAlphWithAmbiguity
    ∷ String → [ST.ShortText] → CharType → Int → [ST.ShortText] → [ST.ShortText] → ([ST.ShortText], [ST.ShortText])
getAlphWithAmbiguity fileName inStates thisType mostDecimals newAlph newStates =
    if null inStates
        then (L.sort $ L.nub newAlph, reverse newStates)
        else
            let firstState = ST.toString $ head inStates
            in  if thisType /= Add
                    then
                        if head firstState /= '['
                            then
                                if firstState `elem` ["?", "-"]
                                    then getAlphWithAmbiguity fileName (tail inStates) thisType mostDecimals newAlph (head inStates : newStates)
                                    else getAlphWithAmbiguity fileName (tail inStates) thisType mostDecimals (head inStates : newAlph) (head inStates : newStates)
                            else -- ambiguity
                            -- let newAmbigStates  = fmap ST.fromString $ words $ filter (`notElem` ['[',']']) firstState

                                let newAmbigStates = fmap ST.fromString $ fmap (: []) $ filter (`notElem` ['[', ']']) firstState
                                in  getAlphWithAmbiguity fileName (tail inStates) thisType mostDecimals (newAmbigStates <> newAlph) (head inStates : newStates)
                    else -- additive character

                        let scaleFactor = 1.0 / (10.0 ** fromIntegral mostDecimals)
                        in  if head firstState /= '['
                                then
                                    if mostDecimals == 0
                                        then getAlphWithAmbiguity fileName (tail inStates) thisType mostDecimals (head inStates : newAlph) (head inStates : newStates)
                                        else
                                            let stateNumber = readMaybe firstState ∷ Maybe Double
                                                newStateNumber = takeWhile (/= '.') $ show (fromJust stateNumber / scaleFactor)
                                            in  if isNothing stateNumber
                                                    then
                                                        if firstState `elem` ["?", "-"]
                                                            then
                                                                getAlphWithAmbiguity
                                                                    fileName
                                                                    (tail inStates)
                                                                    thisType
                                                                    mostDecimals
                                                                    (ST.fromString "-1" : newAlph)
                                                                    (ST.fromString "-1" : newStates)
                                                            else
                                                                errorWithoutStackTrace
                                                                    ("\n\nTNT file " <> fileName <> " ccode processing error: Additive character not a number (Int/Float) " <> firstState)
                                                    else
                                                        getAlphWithAmbiguity
                                                            fileName
                                                            (tail inStates)
                                                            thisType
                                                            mostDecimals
                                                            (ST.fromString newStateNumber : newAlph)
                                                            (ST.fromString newStateNumber : newStates)
                                else -- trace ("GAlphAmb: " <> (show firstState)) $

                                    let hasDecimalorDash = (elem '.' firstState) || (elem '-' firstState)
                                        gutsList =
                                            if hasDecimalorDash
                                                then fmap (: []) $ filter (`notElem` ['[', ']', '.', '-']) firstState
                                                else (: []) <$> filter (`notElem` ['[', ']']) firstState
                                        newStateNumberList = fmap readMaybe gutsList ∷ [Maybe Double]
                                        newStateNumberStringList = fmap (((takeWhile (`notElem` ['.', '-']) . show) . (/ scaleFactor)) . fromJust) newStateNumberList
                                    in  if Nothing `elem` newStateNumberList
                                            then
                                                errorWithoutStackTrace
                                                    ("\n\nTNT file " <> fileName <> " ccode processing error: Additive character not a number (Int/Float) " <> firstState)
                                            else
                                                let newAmbigState =
                                                        if hasDecimalorDash
                                                            then ST.filter (/= ' ') $ ST.fromString $ '[' : unwords newStateNumberStringList <> "]"
                                                            else ST.fromString $ '[' : concat newStateNumberStringList <> "]"
                                                in  getAlphWithAmbiguity
                                                        fileName
                                                        (tail inStates)
                                                        thisType
                                                        mostDecimals
                                                        (fmap ST.fromString newStateNumberStringList <> newAlph)
                                                        (newAmbigState : newStates)

-- )
