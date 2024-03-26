{-# LANGUAGE OverloadedStrings #-}

{- |
Module with functionality to transform phylogenetic data
-}
module Input.DataTransformation (
    renameData,
    getDataTerminalNames,
    addMissingTerminalsToInput,
    checkDuplicatedTerminals,
    createNaiveData,
    createBVNames,
    partitionSequences,
    missingAligned,
    setMissingBits,
    removeAllMissingCharacters,
    checkLeafMissingData,
) where

import Bio.DynamicCharacter
import Bio.DynamicCharacter.Element
import Control.Monad (when)
import Control.Monad.IO.Class (MonadIO (..))
import Data.Coerce (coerce)
import Data.Alphabet
import Data.Alphabet.Codec
import Data.Alphabet.IUPAC
import Data.Alphabet.Special
import Data.Bifunctor
import Data.Bimap (Bimap)
import Data.Bimap qualified as BM
import Bio.DynamicCharacter.Element qualified as Elem
import Data.BitVector.LittleEndian qualified as BV
import Data.BitVector.LittleEndian (BitVector)
import Data.Bits
import Data.Char qualified as C
import Data.Foldable
import Data.Functor (($>))
import Data.List qualified as L
import Data.List.NonEmpty (NonEmpty (..))
import Data.List.NonEmpty qualified as NE
import Data.Maybe
import Data.String
import Data.Text.Lazy qualified as T
import Data.Text.Short qualified as ST
import Data.Vector qualified as V
import Data.Vector.Storable qualified as SV
import Data.Vector.Unboxed qualified as UV
import Data.Vector.Generic qualified as GV
import Foreign.C.Types
import GeneralUtilities
import Numeric.Natural
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import PHANE.Evaluation.Verbosity (Verbosity (..))
import Text.Read
import Types.Types
import Utilities.Utilities qualified as U


-- import Debug.Trace

{- | checkLeafMissingData checks missing data in inputs
missing data is defined as not in an input file
if missing number / numInputFiles > threshold then put in list of
leaves to exclude
-}
checkLeafMissingData ∷ Int → [RawData] → [NameText]
checkLeafMissingData theshold inDataList =
    if theshold == 100
        then []
        else
            let numInputFiles = length inDataList
                criterion = (fromIntegral (100 - theshold)) / 100.0
                minOccurence = (ceiling $ criterion * (fromIntegral numInputFiles)) ∷ Int
                leafList = fmap fst $ concat $ fmap fst inDataList
                groupedLeafList = L.group $ L.sort leafList
                leafOccurence = fmap length groupedLeafList
                leafNameOccurencePair = zip (fmap head groupedLeafList) leafOccurence
                leafToExclude = fmap fst $ filter ((< minOccurence) . snd) leafNameOccurencePair
            in  -- trace ("CLMD :" <> (show (theshold, numInputFiles, criterion, minOccurence, leafOccurence, leafToExclude))) $
                leafToExclude


{- | removeAllMissingCharacters removes characters from list in rawData if all taxa are missing
this can happen when taxa are renamed or added in terminals file
only checks a list length of 1 basically a sequence character
static chars passed on
-}
removeAllMissingCharacters ∷ RawData → PhyG [RawData]
removeAllMissingCharacters inData@(termData, charData) =
    let lengthCheck = exists . filter (> 0) $ length . snd <$> termData
        -- check for non-single sequence character
        multipleSeqChar = length charData /= 1 || ((charType $ head charData) `elem` exactCharacterTypes)
        warnMessage =
            unwords
                [ "Input file"
                , T.unpack . name $ head charData
                , "contains all missing data (perhaps due to renaming or adding/deleting terminals) and has been skipped."
                , "\n"
                ]
        result
            | multipleSeqChar || lengthCheck = pure [inData]
            | otherwise = logWith LogWarn warnMessage $> []
    in  result


{- | partitionSequences takes a character to split sequnces, usually '#'' as in POY, but can be changed
and divides the seqeunces into corresponding partitions.  Replicate character info appending
a number to character name
assumes that input rawdata are a single character (as in form a single file) for sequence data
-}
partitionSequences ∷ ST.ShortText → [RawData] → PhyG [RawData]
partitionSequences partChar inDataList =
    if null inDataList
        then do
            pure []
        else
            let firstRawData@(taxDataList, charInfoList) = head inDataList
            in  -- for raw seqeunce data this will always be a single character
                if (length charInfoList > 1) || (charType (head charInfoList) `notElem` sequenceCharacterTypes)
                    then do
                        restStuff ← partitionSequences partChar (tail inDataList)
                        pure $ firstRawData : restStuff
                    else
                        ( let (leafNameList, leafDataList) = unzip taxDataList
                              partitionCharList = fmap (U.splitSequence partChar) leafDataList
                              partitionCharListByPartition = makePartitionList partitionCharList
                              firstPartNumber = length $ head partitionCharList
                              allSame = filter (== firstPartNumber) $ length <$> tail partitionCharList
                              pairPartitions = zip (fmap T.unpack leafNameList) (fmap length partitionCharList)
                          in  -- check partition numbers consistent + 1 because of tail
                              if (length allSame + 1) /= length partitionCharList
                                then
                                    errorWithoutStackTrace
                                        ( "Number of sequence partitions not consistent in "
                                            <> T.unpack (name $ head charInfoList)
                                            <> " "
                                            <> (show pairPartitions)
                                            <> "\n\tThis may be due to occurence of the partition character '"
                                            <> (show partChar)
                                            <> "' in sequence data.  This value can be changed with the 'set' command."
                                        )
                                else -- if single partition then nothing to do

                                    if firstPartNumber == 1
                                        then do
                                            restStuff ← partitionSequences partChar (tail inDataList)
                                            pure $ firstRawData : restStuff
                                        else -- split data

                                        -- make new structures to create RawData list

                                            let leafNameListList = replicate firstPartNumber leafNameList

                                                -- these filtered from terminal partitions
                                                leafDataListList = fmap (fmap (filter (/= partChar))) partitionCharListByPartition

                                                -- create TermData
                                                newTermDataList = joinLists leafNameListList leafDataListList

                                                -- filter out taxa with empty data-- so can be reconciled proerly later
                                                newTermDataList' = fmap removeTaxaWithNoData newTermDataList

                                                -- create final [RawData]
                                                charInfoListList = replicate firstPartNumber charInfoList
                                                newRawDataList = zip newTermDataList' charInfoListList
                                            in  do
                                                    logWith LogInfo ("\nPartitioning " <> T.unpack (name $ head charInfoList) <> " into " <> show firstPartNumber <> " segments\n")
                                                    restStuff ← partitionSequences partChar (tail inDataList)
                                                    -- trace (" NCI " <> (show $ newTermDataList'))
                                                    pure $ newRawDataList <> restStuff
                        )


{- | removeTaxaWithNoData takes a single TermData list and removes taxa with empty data
these can be created from paritioning sequences where there are no data in a
partitition.  This allows for data reconciliation/renaming later.
-}
removeTaxaWithNoData ∷ [TermData] → [TermData]
removeTaxaWithNoData inTermData =
    if null inTermData
        then []
        else
            let newData = filter (exists . snd) inTermData
            in  -- trace ((show $ length inTermData) <> " -> " <> (show $ length newData))
                newData


{- | joinLists takes two lists of lists (of same length) and zips the
heads of each, then continues till all joined
-}
joinLists ∷ [[a]] → [[b]] → [[(a, b)]]
joinLists listA listB
    | length listA /= length listB = error ("Input lists not equal " <> show (length listA, length listB))
    | null listA = []
    | otherwise =
        let firstList = zip (head listA) (head listB)
        in  firstList : joinLists (tail listA) (tail listB)


-- | makePartitionList take list by taxon and retuns list by partition
makePartitionList ∷ [[[ST.ShortText]]] → [[[ST.ShortText]]]
makePartitionList inListList =
    if null $ head inListList
        then []
        else
            let firstParList = fmap head inListList
            in  firstParList : makePartitionList (fmap tail inListList)


{- | renameData takes a list of rename Text pairs (new name, oldName)
and replaces the old name with the new
-}
renameData ∷ [(T.Text, T.Text)] → RawData → RawData
renameData newNamePairList inData =
    if null newNamePairList
        then inData
        else
            let terminalData = fst inData
            in  if null terminalData
                    then inData
                    else
                        let newTerminalData = fmap (relabelterminalData newNamePairList) terminalData
                        in  (newTerminalData, snd inData)


{- | relabelterminalData takes a list of Text pairs and the terminals with the
second name in the pairs is changed to the first
-}
relabelterminalData ∷ [(T.Text, T.Text)] → TermData → TermData
relabelterminalData namePairList terminalData@(leafName, leafData) =
    if null namePairList
        then terminalData
        else
            let foundName = find ((== leafName) . snd) namePairList
            in  if isNothing foundName
                    then -- trace ("Not renaming " <> (T.unpack leafName)) --  <> " " <> show namePairList)
                        terminalData
                    else -- trace ("Renaming " <> (T.unpack leafName) <> " to " <> (T.unpack $ fst $ fromJust foundName))
                        (fst $ fromJust foundName, leafData)


{- | getDataTerminalNames takes all input data and gets full terminal list
and adds missing data for terminals not in input files
-}
getDataTerminalNames ∷ [RawData] → [T.Text]
getDataTerminalNames inDataList =
    if null inDataList
        then []
        else L.sort $ L.nub $ fst <$> concatMap fst inDataList


-- | addMissingTerminalsToInput dataLeafNames renamedData
addMissingTerminalsToInput ∷ NonEmpty T.Text → [TermData] → RawData → RawData
addMissingTerminalsToInput dataLeafNames@(firstLeafName :| xs) curTermData inData@(termDataList, charInfoList) =
    let rData = case find ((== firstLeafName) . fst) termDataList of
            Just vs → vs : curTermData
            Nothing → (firstLeafName, []) : curTermData
    in  case xs of
            [] → (reverse rData, charInfoList)
            y : ys → addMissingTerminalsToInput (y :| ys) rData inData


-- | checkDuplicatedTerminals takes list TermData and checks for repeated terminal names
checkDuplicatedTerminals ∷ [TermData] → (Bool, [T.Text])
checkDuplicatedTerminals inData =
    if null inData
        then (False, [])
        else
            let nameList = L.group $ L.sort $ fmap fst inData
                dupList = filter ((> 1) . length) nameList
            in  if null dupList
                    then (False, [])
                    else (True, fmap head dupList)


{- | joinSortFileData takes list of list of short text and merges line by line to join leaf states
and sorts the result
-}
joinSortFileData ∷ [[ST.ShortText]] → [String]
joinSortFileData inFileLists =
    if null (head inFileLists)
        then []
        else
            let -- changed sort order (now by data file input, more or less) to reduce time complexity
                -- should still be label invariant
                -- firstLeaf = L.sort $ ST.toString $ ST.concat $ fmap head inFileLists
                firstLeaf = ST.toString $ ST.concat $ fmap head inFileLists
            in  -- firstLeaf = show $ L.sort $ fmap head inFileLists

                firstLeaf : joinSortFileData (fmap tail inFileLists)


{- | createBVNames takes input data, sorts the raw data, hashes, sorts those to create
unique, label invariant (but data related so arbitrary but consistent)
Assumes the rawData come in sorted by the data reconciliation process
These used for vertex labels, caching, left/right DO issues
-}
createBVNames ∷ [RawData] → [(T.Text, BitVector)]
createBVNames inDataList =
    let rawDataList = fmap fst inDataList
        textNameList = fst <$> head rawDataList
        textNameList' = fst <$> last rawDataList

        fileLeafCharList = fmap (fmap snd) rawDataList
        fileLeafList = fmap (fmap ST.concat) fileLeafCharList
        leafList = reverse $ joinSortFileData fileLeafList

        -- hash not guaranteed to be stable over OS or library version
        -- leafHash = fmap H.hash leafList
        leafHash = leafList
        leafHashPair = L.sortOn fst $ zip leafHash [0 .. (length textNameList - 1)] -- textNameList
        (_, leafReoderedList) = unzip leafHashPair
        -- leafOrder = sortOn fst $ zip leafReoderedList [0..((length textNameList) - 1)]
        -- (nameList, intList) = unzip leafOrder

        -- bv1 = BV.bitVec (length textNameList) (1 :: Integer)
        boolList = replicate (length textNameList - 1) False
        bv1 = BV.fromBits $ True : boolList
        bvList = fmap (shiftL bv1) leafReoderedList -- [0..((length textNameList) - 1)]
    in  if textNameList /= textNameList'
            then error "Taxa are not properly ordered in createBVNames"
            else -- trace (show $ fmap BV.toBits bvList)
                zip textNameList bvList


{- | createNaiveData takes input RawData and transforms to "Naive" data.
these data are organized into blocks (set to input filenames initially)
and are bitvector coded, but are not organized by character type, packed ot
optimized in any other way (prealigned-> nonadd, Sankoff.  2 state sankoff to binary,
constant characters skipped etc)
these processes take place later
these data can be input to any data optimization commands and are useful
for data output as they haven't been reordered or transformed in any way.
the RawData is a list since it is organized by input file
the list accumulator is to avoid Vector snoc/cons O(n)
-}
createNaiveData ∷ GlobalSettings → [RawData] → [(T.Text, BitVector)] → [BlockData] → PhyG ProcessedData
createNaiveData inGS inDataList leafBitVectorNames curBlockData =
    -- trace ("CND: " <> (show $ (optimalityCriterion inGS))) $
    if null inDataList
        then -- trace ("Naive data with " <> (show $ length curBlockData) <> " blocks and " <> (show $ fmap length $ fmap V.head $ fmap snd3 curBlockData) <> " characters")

            pure
                ( V.fromList $ fmap fst leafBitVectorNames
                , V.fromList $ fmap snd leafBitVectorNames
                , V.fromList $ reverse curBlockData
                )
        else
            let (firstData, firstCharInfo) = head inDataList
            in  -- empty file should have been caught earlier, but avoids some head/tail errors
                if null firstCharInfo
                    then -- trace "Empty CharInfo"
                        createNaiveData inGS (tail inDataList) leafBitVectorNames curBlockData
                    else -- process data as come in--each of these should be from a single file
                    -- and initially assigned to a single, unique block

                        let thisBlockName = name $ head firstCharInfo
                            -- thisBlockName   = T.append (T.pack "block-")  (name $ head firstCharInfo)
                            thisBlockCharInfo = V.fromList firstCharInfo
                            maxCharacterLength = maximum $ fmap length $ (fmap snd firstData)
                            recodedCharacters = recodeRawData (fmap fst firstData) (fmap snd firstData) firstCharInfo maxCharacterLength []
                            -- thisBlockGuts = V.zip (V.fromList $ fmap snd leafBitVectorNames) recodedCharacters
                            previousBlockName =
                                if not $ null curBlockData
                                    then fst3 $ head curBlockData
                                    else T.empty
                            thisBlockName' =
                                if T.takeWhile (/= '#') previousBlockName /= T.takeWhile (/= '#') thisBlockName
                                    then thisBlockName
                                    else
                                        let oldSuffix = T.dropWhile (/= '#') previousBlockName
                                            indexSuffix =
                                                if T.null oldSuffix
                                                    then T.pack "#0"
                                                    else
                                                        let oldIndex = readMaybe (T.unpack $ T.tail oldSuffix) ∷ Maybe Int
                                                            newIndex = 1 + fromJust oldIndex
                                                        in  if isNothing oldIndex
                                                                then error "Bad suffix in createNaiveData"
                                                                else T.pack ("#" <> show newIndex)
                                        in  T.append (T.takeWhile (/= '#') thisBlockName) indexSuffix

                            thisBlockCharInfo'' = V.zipWith (resetAddNonAddAlphabets recodedCharacters) thisBlockCharInfo (V.fromList [0 .. (V.length thisBlockCharInfo - 1)])

                            -- create "orginal" character info for later use in outputs after character recoding and transformation etc.
                            thisBlockCharInfo' = fmap setOrigCharInfo thisBlockCharInfo''

                            -- recode with appropriate missing data codes
                            recodedCharacters' = fmap (recodeNonAddMissingBlock thisBlockCharInfo') recodedCharacters

                            -- reweight characters by NCM (Tuffley and Steel, 1997) factor (does not include root factor--must be separately calculated)
                            thisBlockCharInfoNCM =
                                if (optimalityCriterion inGS) == NCM
                                    then fmap reweightNCM thisBlockCharInfo'
                                    else thisBlockCharInfo'

                            thisBlockData = (thisBlockName', recodedCharacters', thisBlockCharInfoNCM)

                            (prealignedDataEqualLength, nameMinPairList, nameNonMinPairList) = checkPrealignedEqualLength (fmap fst leafBitVectorNames) thisBlockData
                        in  -- trace ("CND:" <> (show $ fmap length $ (fmap snd firstData))) (
                            if not prealignedDataEqualLength
                                then
                                    errorWithoutStackTrace
                                        ( "Error on input of prealigned sequence characters in file "
                                            <> (takeWhile (/= '#') $ T.unpack thisBlockName')
                                            <> "--not equal length [(Taxon, Length)]: \nMinimum length taxa: "
                                            <> (show nameMinPairList)
                                            <> "\nNon Minimum length taxa: "
                                            <> (show nameNonMinPairList)
                                        )
                                else
                                    logWith LogInfo ("Recoding input block: " <> T.unpack thisBlockName' <> "\n")
                                        *> createNaiveData inGS (tail inDataList) leafBitVectorNames (thisBlockData : curBlockData)


-- | reweightNCM takes character info and reweights via NCM (Tuffley and Steel, 1997) -log_10 1/(alphabet size)
reweightNCM ∷ CharInfo → CharInfo
reweightNCM inCharInfo =
    let originalWeight = weight inCharInfo
        alphabetSize = length $ alphabet inCharInfo
        newWeight = originalWeight * (-1.0) * logBase 10 (1.0 / (fromIntegral alphabetSize))
    in  -- trace ("RWNCM: " <> (show (originalWeight, newWeight)))
        inCharInfo{weight = newWeight}


{- | setOrigCharInfo takes fields from charInfo and sets the initial original charcter infomatin field
as a singleton Vector
-}
setOrigCharInfo ∷ CharInfo → CharInfo
setOrigCharInfo inCharInfo =
    let origData = (name inCharInfo, charType inCharInfo, alphabet inCharInfo)
    in  inCharInfo{origInfo = V.singleton origData}


-- | recodeAddNonAddMissing takes Block data and recodes missing for additive and non-additive characters
recodeNonAddMissingBlock ∷ V.Vector CharInfo → V.Vector CharacterData → V.Vector CharacterData
recodeNonAddMissingBlock blockCharInfo singleTaxonBlockData =
    V.zipWith recodeNonAddMissingCharacter blockCharInfo singleTaxonBlockData


-- | recodeAddNonAddMissingCharacter recodes additive and non-additive missing data
recodeNonAddMissingCharacter ∷ CharInfo → CharacterData → CharacterData
recodeNonAddMissingCharacter charInfo inCharData =
    let inCharType = charType charInfo
    in  if inCharType /= NonAdd
            then inCharData
            else
                let nonAddState = (V.head . snd3 . stateBVPrelim) inCharData
                    newState =
                        if nonAddState == (nonAddState `xor` nonAddState)
                            then V.singleton (complement nonAddState)
                            else V.singleton nonAddState
                in  if nonAddState == (nonAddState `xor` nonAddState)
                        then
                            inCharData
                                { stateBVPrelim = (newState, newState, newState)
                                , stateBVFinal = newState
                                }
                        else inCharData


{- | getAddNonAddAlphabets takes recoded character data and resets the alphabet
field in charInfo to reflect observed states.  This is used to properly set missing and
bit packing values
-}
resetAddNonAddAlphabets ∷ V.Vector (V.Vector CharacterData) → CharInfo → Int → CharInfo
resetAddNonAddAlphabets taxonByCharData charInfo charIndex =
    let inCharType = charType charInfo
    in  if inCharType `notElem` [Add, NonAdd]
            then charInfo
            else
                if inCharType == NonAdd
                    then
                        let -- get actual states
                            inCharV = fmap (V.head . snd3 . stateBVPrelim) $ fmap V.head taxonByCharData
                            missingVal = BV.fromBits $ L.replicate (fromEnum . BV.dimension . V.head . snd3 . stateBVPrelim . V.head $ V.head taxonByCharData) True
                            -- missingVal = V.foldl1' (.|.) inCharV
                            nonMissingBV = V.foldl1' (.|.) $ V.filter (/= missingVal) inCharV

                            -- max in case of all missing character
                            numStates = max 1 $ popCount nonMissingBV

                            -- numBits = UV.length . coerce $ (V.head . snd3 . stateBVPrelim) $ (V.head taxonByCharData) V.! charIndex
                            foundSymbols = ST.fromString . show <$> (0 :| [1 .. pred numStates])
                            stateAlphabet = fromSymbols foundSymbols -- fromSymbolsWOGap foundSymbols
                        in  -- trace ("RNA: " <> (show stateAlphabet))
                            charInfo{alphabet = stateAlphabet}
                    else
                        if inCharType == Add
                            then
                                let (minRangeL, maxRangeL) = V.unzip $ fmap (V.head . snd3 . rangePrelim) $ fmap (V.! charIndex) taxonByCharData

                                    minRange =
                                        if minimum minRangeL < (maxBound ∷ Int)
                                            then minimum minRangeL
                                            else 0

                                    maxRange =
                                        if maximum maxRangeL > (minBound ∷ Int)
                                            then maximum maxRangeL
                                            else 0

                                    foundSymbols = ST.fromString . show <$> (minRange :| [succ minRange .. maxRange])
                                    stateAlphabet = fromSymbols foundSymbols -- fromSymbolsWOGap foundSymbols
                                in  if maxRange < minRange
                                        then error ("Error in processing of additive character states " <> (show (minRange, maxRange)))
                                        else -- trace ("RA: " <> (show (minimum minRangeL, maximum maxRangeL)) <> " -> " <> (show (minRange, maxRange)) <> " " <> (show foundSymbols)) -- <> " -> " <> (show stateAlphabet))
                                            charInfo{alphabet = stateAlphabet}
                            else error ("Unrecognized character type in resetAddNonAddAlphabets: " <> (show inCharType))


{- | checkPrealignedEqualLength checks prealigned type for equal length
at this stage (called before reblocking) there should only be a single charcter per block
but more general--if not great if > 1 character with naming
-}
checkPrealignedEqualLength
    ∷ [T.Text] → (NameText, V.Vector (V.Vector CharacterData), V.Vector CharInfo) → (Bool, [(T.Text, Int)], [(T.Text, Int)])
checkPrealignedEqualLength nameTextList (_, taxByCharacterDataVV, charInfoV) =
    let numCharsIndexList = [0 .. (V.length charInfoV) - 1]
        sameLengthPairList = zipWith (verifyPrealignedCharacterLength nameTextList taxByCharacterDataVV) (V.toList charInfoV) numCharsIndexList
        badOnes = filter ((== False) . fst3) sameLengthPairList
    in  if null badOnes
            then (True, [], [])
            else (False, concat $ fmap snd3 badOnes, concat $ fmap thd3 badOnes)


{- | verifyPrealignedCharacterLength takes an index for character and examines theat character--if prealigned checks for
equal length if prealigned then "True"
-}
verifyPrealignedCharacterLength
    ∷ [T.Text] → V.Vector (V.Vector CharacterData) → CharInfo → Int → (Bool, [(T.Text, Int)], [(T.Text, Int)])
verifyPrealignedCharacterLength nameTextList taxByCharacterDataVV charInfo charIndex =
    let inCharType = charType charInfo
        inCharV = fmap (V.! charIndex) taxByCharacterDataVV
    in  if inCharType `notElem` prealignedCharacterTypes
            then (True, [], [])
            else -- True for IA field use but doesn't matter in this case since looking at prealigned fields

                let prealigedDataLengthList = V.toList $ fmap (U.getCharacterLength' True charInfo) inCharV
                    {-
                    if inCharType == AlignedSlim then V.toList $ fmap SV.length $ fmap (snd3 . alignedSlimPrelim) inCharV
                      else if inCharType == AlignedWide then V.toList $ fmap UV.length $ fmap (snd3 . alignedWidePrelim) inCharV
                      else if inCharType == AlignedHuge then V.toList $ fmap V.length $ fmap (snd3 . alignedHugePrelim) inCharV
                      else  error ("Character type " <> show inCharType <> " unrecongized/not implemented")
                      -}

                    nameLengthPairList = zip nameTextList prealigedDataLengthList
                    lMinLength = minimum prealigedDataLengthList
                    haveMinLength = filter ((== lMinLength) . snd) nameLengthPairList
                    notMinMinLength = filter ((/= lMinLength) . snd) nameLengthPairList
                in  -- all min length then all same length
                    -- trace ("VPCL:" <> (show $ (haveMinLength, notMinMinLength))) (
                    if null notMinMinLength
                        then (True, [], [])
                        else (False, haveMinLength, notMinMinLength)


-- )

{- | recodeRawData takes the ShortText representation of character states/ranges etc
and recodes the appropriate fields in CharacterData (from Types)
the list accumulator is to avoid Vector cons/snoc O(n)
differentiates between seqeunce type and others with char info
-}
recodeRawData ∷ [NameText] → [[ST.ShortText]] → [CharInfo] → Int → [[CharacterData]] → V.Vector (V.Vector CharacterData)
recodeRawData inTaxNames inData inCharInfo maxCharLength curCharData =
    -- trace ("RRD: ") (
    if null inTaxNames
        then V.fromList $ reverse $ fmap V.fromList curCharData
        else
            let firstData = head inData
                firstDataRecoded = createLeafCharacter inCharInfo firstData maxCharLength
            in  -- trace ("RRD:" <> (show firstData))
                -- trace ("Recoding " <> (T.unpack $ head inTaxNames) <> " as " <> (show $ charType $ head inCharInfo) <> "\n\t" <> show firstDataRecoded)
                -- trace ((show $ length inData) <> " " <> (show $ length firstData) <> " " <> (show $ length inCharInfo)
                recodeRawData (tail inTaxNames) (tail inData) inCharInfo maxCharLength (firstDataRecoded : curCharData)


-- )

-- | missingNonAdditive is non-additive missing character value, all 1's based on alphabet size
missingNonAdditive ∷ CharInfo → CharacterData
missingNonAdditive inCharInfo =
    let missingChar = V.singleton (BV.fromBits $ replicate (length $ alphabet inCharInfo) True)
    in  emptyCharacter
            { stateBVPrelim = (missingChar, missingChar, missingChar)
            , stateBVFinal = missingChar
            }


-- | missingAdditive is additive missing character value, all 1's based on alphabet size
missingAdditive ∷ CharInfo → CharacterData
missingAdditive inCharInfo =
    let minState' = readMaybe (ST.toString . head . toList $ alphabet inCharInfo) ∷ Maybe Int
        maxState' = readMaybe (ST.toString . last . toList $ alphabet inCharInfo) ∷ Maybe Int
        minState =
            if isNothing minState'
                then 0
                else fromJust minState'
        maxState =
            if isNothing maxState'
                then 0
                else fromJust maxState'

        missingRange =
            V.zip
                (V.singleton minState)
                (V.singleton maxState)
    in  emptyCharacter
            { rangePrelim = (missingRange, missingRange, missingRange)
            , rangeFinal = missingRange
            }


{- | missingMatrix is additive missing character value, all 1's based on alphabet size
setrting stateBVPrelim/Final for approx DO-like costs (lookup)
-}
missingMatrix ∷ CharInfo → CharacterData
missingMatrix inCharInfo =
    let numStates = length $ alphabet inCharInfo
        missingState = (0 ∷ StateCost, [], [])
    in  emptyCharacter
            { matrixStatesPrelim = V.singleton (V.replicate numStates missingState)
            , matrixStatesFinal = V.singleton (V.replicate numStates missingState)
            }


-- | getMissingValue takes the character type and returns the appropriate missing data value
getMissingValue ∷ [CharInfo] → Int → [CharacterData]
getMissingValue inChar maxCharLength
    | null inChar = []
    | charType (head inChar) `elem` nonExactCharacterTypes = [emptyCharacter] -- []
    | charType (head inChar) `elem` prealignedCharacterTypes = [missingAligned (head inChar) maxCharLength]
    | charType (head inChar) == NonAdd = missingNonAdditive (head inChar) : getMissingValue (tail inChar) maxCharLength
    | charType (head inChar) == Add = missingAdditive (head inChar) : getMissingValue (tail inChar) maxCharLength
    | charType (head inChar) == Matrix = missingMatrix (head inChar) : getMissingValue (tail inChar) maxCharLength
    | otherwise = error ("Datatype " <> show (charType $ head inChar) <> " not recognized")


{- | missingAligned creates missing data (all bits on) for prealigned data
important n ot to sett all bits--then run of in to seg fault land
-}
missingAligned ∷ CharInfo → Int → CharacterData
missingAligned inChar charLength =
    let alphabetLen = length $ alphabet inChar
        missingVals x = GV.replicate charLength $ setMissingBits x 0 alphabetLen
        missingChar x = let v = missingVals x in ( v, v, v )
    in  case charType inChar of
            AlignedSlim -> emptyCharacter{alignedSlimPrelim = missingChar (0 ∷ SlimState) }
            AlignedWide -> emptyCharacter{alignedWidePrelim = missingChar (0 ∷ WideState) }
            AlignedHuge -> emptyCharacter{alignedHugePrelim = missingChar (fromNumber (toEnum alphabetLen) 0 ∷ HugeState) }
            val ->  error $ unwords [ "Datatype", show val, "not recognized" ]


-- | setMissingBits sets the first bits by index to '1' rest left as is (0 on input)
setMissingBits ∷ (Show a, FiniteBits a) ⇒ a → Int → Int → a
setMissingBits inVal curIndex alphSize
    | curIndex == alphSize = inVal
    | otherwise = setMissingBits (setBit inVal curIndex) (curIndex + 1) alphSize


{- | getStateBitVectorList takes the alphabet of a character ([ShorText])
and returns bitvectors (with of size alphabet) for each state in order of states in alphabet
-}
getStateBitVectorList ∷ Alphabet ST.ShortText → V.Vector (ST.ShortText, BitVector)
getStateBitVectorList localStates
    | null localStates = error "Character with empty alphabet in getStateBitVectorList"
    | otherwise =
            let stateCount = toEnum $ length localStates
                stateIndexList = V.fromList [0 .. stateCount - 1]
                genNum = (2 ^) ∷ Word → Natural
                bvList ∷ V.Vector BitVector
                bvList = (BV.fromNumber stateCount . genNum) <$> stateIndexList
            in  V.zip (alphabetStateNames localStates) bvList


iupacToBVPairs
    ∷ (IsString s, Ord s)
    ⇒ Alphabet s
    → Bimap (NonEmpty s) (NonEmpty s)
    → V.Vector (s, BitVector)
iupacToBVPairs inputAlphabet iupac = V.fromList $ bimap NE.head encoder <$> BM.toAscList iupac
    where
        constructor = flip BV.fromNumber (0 ∷ Int)
        encoder = encodeState inputAlphabet constructor


{- | nucleotideBVPairs for recoding DNA sequences
this done to insure not recalculating everything for each base
-}
nucleotideBVPairs ∷ V.Vector (ST.ShortText, BitVector)
nucleotideBVPairs = iupacToBVPairs baseAlphabet iupacToDna
    where
        baseAlphabet = fromSymbols $ ST.fromString <$> "-" :| ["A", "C", "G", "T"]


{- | aminoAcidBVPairs for recoding protein sequences
this done to insure not recalculating everything for each residue
B, Z, X, ? for ambiguities
-}
aminoAcidBVPairs ∷ V.Vector (ST.ShortText, BitVector)
aminoAcidBVPairs = iupacToBVPairs acidAlphabet iupacToAminoAcid
    where
        acidAlphabet =
            fromSymbols $
                fromString
                    <$> "A" :| ["C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-"]


{- | getBVCode take a Vector of (ShortText, BV) and returns bitvector code for
ShortText state
-}
getBVCode ∷ V.Vector (ST.ShortText, BitVector) → ST.ShortText → BitVector
getBVCode bvCodeVect inState =
    let newCode = V.find ((== inState) . fst) bvCodeVect
    in  maybe (error ("State " <> ST.toString inState <> " not found in bitvect code " <> show bvCodeVect)) snd newCode


getNucleotideSequenceChar ∷ Bool → [ST.ShortText] → [CharacterData]
getNucleotideSequenceChar isPrealigned stateList =
    let sequenceVect
            | null stateList = mempty
            | otherwise =
                SV.fromList $
                    BV.toUnsignedNumber . getBVCode nucleotideBVPairs <$> (fmap ST.pack $ fmap (fmap C.toUpper) $ fmap ST.unpack stateList)
        newSequenceChar =
            if not isPrealigned
                then
                    emptyCharacter
                        { slimPrelim = sequenceVect
                        , slimGapped = (sequenceVect, sequenceVect, sequenceVect)
                        }
                else
                    emptyCharacter
                        { alignedSlimPrelim = (sequenceVect, sequenceVect, sequenceVect)
                        }
    in  [newSequenceChar]


getAminoAcidSequenceChar ∷ Bool → [ST.ShortText] → [CharacterData]
getAminoAcidSequenceChar isPrealigned stateList =
    let sequenceVect
            | null stateList = mempty
            | otherwise =
                UV.fromList $
                    BV.toUnsignedNumber . getBVCode aminoAcidBVPairs <$> (fmap ST.pack $ fmap (fmap C.toUpper) $ fmap ST.unpack stateList)
        newSequenceChar =
            if not isPrealigned
                then
                    emptyCharacter
                        { widePrelim = sequenceVect
                        , wideGapped = (sequenceVect, sequenceVect, sequenceVect)
                        }
                else
                    emptyCharacter
                        { alignedWidePrelim = (sequenceVect, sequenceVect, sequenceVect)
                        }
    in  [newSequenceChar]


{- | getGeneralBVCode take a Vector of (ShortText, BV) and returns bitvector code for
ShortText state.  These states can be ambiguous as in general sequences
so states need to be parsed first
the AND to all states makes ambiguityoes only observed states
-}
getGeneralBVCode ∷ V.Vector (ST.ShortText, BitVector) → ST.ShortText → (SlimState, WideState, HugeState)
getGeneralBVCode bvCodeVect inState =
    let inStateString = ST.toString inState
        alphabetLen = toEnum $ length bvCodeVect
        construct x = (BV.toUnsignedNumber x, BV.toUnsignedNumber x, Elem.fromNumber alphabetLen $ BV.toUnsignedNumber x) 
    in  -- if '[' `notElem` inStateString then --single state
        if (head inStateString /= '[') && (last inStateString /= ']') -- single state
            then
                let newCode = V.find ((== inState) . fst) bvCodeVect
                    allBVStates = V.foldl1' (.|.) (fmap snd bvCodeVect)
                    bvDimension = fromEnum . BV.dimension . snd $ V.head bvCodeVect
                in  case newCode of
                        Just (_,x) -> construct x
                        Nothing -> case inState of

                            -- B is Aspartic Acid or Asparagine if '-' =0 then states 3 and 12.
                            "B" ->
                                let x = BV.fromBits $ (replicate 3 False) <> [True] <> (replicate 8 False) <> [True] <> (replicate (bvDimension - 13) False)
                                in  construct x
                            -- any amino acid but not '-'
                            "X" -> 
                                let x = allBVStates .&. BV.fromBits (False : (replicate (bvDimension - 1) True))
                                in  construct x

                            -- any state including '-'
                            "?" -> 
                                let x = allBVStates .&. (BV.fromBits (replicate bvDimension True))
                                in  construct x

                            _ -> error $ unwords [ "State", ST.toString inState, "not found in bitvect code", show bvCodeVect ]

            else
                let statesStringList = words $ tail $ init inStateString
                    stateList = fmap ST.fromString statesStringList
                    maybeBVList = fmap getBV stateList
                    stateBVList = fmap (snd . fromJust) maybeBVList
                    ambiguousBVState = foldr1 (.|.) stateBVList
                in  if Nothing `elem` maybeBVList
                        then error ("Ambiguity group " <> inStateString <> " contained states not found in bitvect code " <> show bvCodeVect)
                        else construct ambiguousBVState
    where
        getBV s = V.find ((== s) . fst) bvCodeVect


{- | getGeneralSequenceChar encode general (ie not nucleotide or amino acid) sequences
as bitvectors.  Main difference with getSequenceChar is in dealing wioth ambiguities
they need to be parsed and "or-ed" differently
need to have all three preliminary fields populatied for some reason--prob shouldn't need that
-}
getGeneralSequenceChar ∷ CharInfo → [ST.ShortText] → [CharacterData]
getGeneralSequenceChar inCharInfo stateList =
    let cType = charType inCharInfo
        -- isAligned = prealigned inCharInfo
        stateBVPairVect :: V.Vector (ST.ShortText, BitVector)
        stateBVPairVect = getStateBitVectorList $ alphabet inCharInfo
        (slimVec, wideVec, hugeVec)
            | not $ null stateList = 
                    (\(x, y, z) → (SV.fromList $ toList x, UV.fromList $ toList y, z)) . V.unzip3 . V.fromList $
                        fmap (getGeneralBVCode stateBVPairVect) stateList
            | otherwise = (mempty, mempty, mempty)
        newSequenceChar =
            emptyCharacter
                { slimPrelim = if cType `elem` [SlimSeq, NucSeq] then slimVec else mempty
                , slimGapped = if cType `elem` [SlimSeq, NucSeq] then (slimVec, slimVec, slimVec) else (mempty, mempty, mempty)
                , slimFinal = if cType `elem` [SlimSeq, NucSeq] then slimVec else mempty
                , widePrelim = if cType `elem` [WideSeq, AminoSeq] then wideVec else mempty
                , wideGapped = if cType `elem` [WideSeq, AminoSeq] then (wideVec, wideVec, wideVec) else (mempty, mempty, mempty)
                , wideFinal = if cType `elem` [WideSeq, AminoSeq] then wideVec else mempty
                , hugePrelim = if cType == HugeSeq then hugeVec else mempty
                , hugeGapped = if cType == HugeSeq then (hugeVec, hugeVec, hugeVec) else (mempty, mempty, mempty)
                , hugeFinal = if cType == HugeSeq then hugeVec else mempty
                , alignedSlimPrelim = if cType `elem` [AlignedSlim] then (slimVec, slimVec, slimVec) else (mempty, mempty, mempty)
                , alignedSlimFinal = if cType `elem` [AlignedSlim] then slimVec else mempty
                , alignedWidePrelim = if cType `elem` [AlignedWide] then (wideVec, wideVec, wideVec) else (mempty, mempty, mempty)
                , alignedWideFinal = if cType `elem` [AlignedWide] then wideVec else mempty
                , alignedHugePrelim = if cType `elem` [AlignedHuge] then (hugeVec, hugeVec, hugeVec) else (mempty, mempty, mempty)
                , alignedHugeFinal = if cType `elem` [AlignedHuge] then hugeVec else mempty
                }
    in  -- trace ("GGSC" <> (show stateList) <> "\n" <> (show newSequenceChar ))
        [newSequenceChar]


{- | getStateBitVector takes teh alphabet of a character ([ShorText])
and returns then bitvectorfor that state in order of states in alphabet
-}
getStateBitVector ∷ Alphabet ST.ShortText → ST.ShortText → BitVector
getStateBitVector localAlphabet = encodeState localAlphabet (const constructor) . (: [])
    where
        constructor :: BitVector
        constructor =
            let x = bit $ length localAlphabet - 1
            in  x `xor` x


-- getMinMaxStates takes  list of strings and determines the minimum and maximum integer values
getMinMaxStates ∷ [String] → (Int, Int) → (Int, Int)
getMinMaxStates inStateStringList (curMin, curMax) =
    if null inStateStringList
        then (curMin, curMax)
        else
            let firstString = head inStateStringList
            in  -- missing data
                if firstString == "-" || firstString == "?"
                    then getMinMaxStates (tail inStateStringList) (curMin, curMax)
                    else -- single state

                        if '[' `notElem` firstString
                            then
                                let onlyInt = readMaybe firstString ∷ Maybe Int
                                in  if isNothing onlyInt
                                        then error ("State not an integer in getIntRange: " <> firstString)
                                        else
                                            let minVal =
                                                    if fromJust onlyInt < curMin
                                                        then fromJust onlyInt
                                                        else curMin
                                                maxVal =
                                                    if fromJust onlyInt > curMax
                                                        then fromJust onlyInt
                                                        else curMax
                                            in  getMinMaxStates (tail inStateStringList) (minVal, maxVal)
                            else -- range of states

                                let statesStringList = words $ tail $ init firstString
                                    stateInts = fmap readMaybe statesStringList ∷ [Maybe Int]
                                in  if Nothing `elem` stateInts
                                        then error ("Non-integer in range " <> firstString)
                                        else
                                            let localMin = minimum $ fmap fromJust stateInts
                                                localMax = maximum $ fmap fromJust stateInts
                                                minVal =
                                                    if localMin < curMin
                                                        then localMin
                                                        else curMin
                                                maxVal =
                                                    if localMax > curMax
                                                        then localMax
                                                        else curMax
                                            in  getMinMaxStates (tail inStateStringList) (minVal, maxVal)


-- getIntRange takes the local states and returns the Integer range of an additive character
-- in principle allows for > 2 states
getIntRange ∷ ST.ShortText → Alphabet ST.ShortText → (Int, Int)
getIntRange localState totalAlphabet =
    let stateString = ST.toString localState
    in  -- single state
        if (stateString == "?") || (stateString == "-")
            then getMinMaxStates (ST.toString <$> toList totalAlphabet) (maxBound ∷ Int, minBound ∷ Int)
            else
                if '[' `notElem` stateString
                    then
                        let onlyInt = readMaybe stateString ∷ Maybe Int
                        in  if isNothing onlyInt
                                then error ("State not an integer in getIntRange: " <> ST.toString localState)
                                else (fromJust onlyInt, fromJust onlyInt)
                    else -- Range of states

                        let hasDash = ST.any (== '-') localState
                            statesStringList =
                                if hasDash
                                    then fmap ST.toString $ fmap (ST.filter (`notElem` ['[', ']'])) $ ST.split (== '-') localState
                                    else fmap (: []) $ ST.toString $ ST.filter (`notElem` ['[', ']']) localState
                            -- words $ tail $ init stateString
                            stateInts = fmap readMaybe statesStringList ∷ [Maybe Int]
                        in  -- trace ("GIR:" <> (show localState) <> " -> " <> (show (minimum $ fmap fromJust stateInts, maximum $ fmap fromJust stateInts))) (
                            if Nothing `elem` stateInts
                                then error ("Non-integer in range " <> ST.toString localState)
                                else (minimum $ fmap fromJust stateInts, maximum $ fmap fromJust stateInts)


-- )

-- | getTripleList
getTripleList ∷ MatrixTriple → MatrixTriple → [ST.ShortText] → [ST.ShortText] → [MatrixTriple]
getTripleList hasState notHasState localAlphabet stateList =
    if null localAlphabet
        then []
        else
            let firstAlphState = head localAlphabet
            in  if firstAlphState `elem` stateList
                    then -- trace ("State " <> show firstAlphState <> " in " <> show localAlphabet)
                        hasState : getTripleList hasState notHasState (tail localAlphabet) stateList
                    else notHasState : getTripleList hasState notHasState (tail localAlphabet) stateList


-- | getInitialMatrixVector gets matrix vector
getInitialMatrixVector ∷ Alphabet ST.ShortText → ST.ShortText → V.Vector MatrixTriple
getInitialMatrixVector alphabet' localState =
    let hasState = (0 ∷ StateCost, [], [])
        notHasState = (maxBound ∷ StateCost, [], [])
        localAlphabet = toList alphabet'
    in  let stateString = ST.toString localState
        in  -- single state
            if '[' `notElem` stateString
                then -- trace ("GIMV: " <> (show $ V.fromList $ getTripleList hasState notHasState localAlphabet [localState]))
                    V.fromList $ getTripleList hasState notHasState localAlphabet [localState]
                else -- polylorphic/ambiguous

                    let statesStringList = words $ tail $ init stateString
                        stateList = fmap ST.fromString statesStringList
                    in  V.fromList $ getTripleList hasState notHasState localAlphabet stateList


{- | getQualitativeCharacters processes non-sequence characters (non-additive, additive, sankoff/matrix)
and recodes returning list of encoded characters
reverses order due to prepending
matrix stateBVPrelim/Final for approx matrix costs
adddded in code for ambiguities for non-additive--somehow got lost--alphabet [robbaly wrong now as well
-}
getQualitativeCharacters ∷ [CharInfo] → [ST.ShortText] → [CharacterData] → [CharacterData]
getQualitativeCharacters inCharInfoList inStateList curCharList =
    if null inCharInfoList
        then reverse curCharList
        else
            let firstCharInfo = head inCharInfoList
                firstState = head inStateList
                firstCharType = charType firstCharInfo
                totalAlphabet = alphabet firstCharInfo
            in  -- single state
                if firstCharType == NonAdd
                    then
                        let stateBV
                                | ST.length firstState == 1 = getStateBitVector (alphabet firstCharInfo) firstState
                                | otherwise =
                                        let ambiguousStateST = ST.filter (`notElem` ['[', ']']) firstState
                                            ambiguousStateString = ST.toString ambiguousStateST
                                            stateSTList = fmap ST.singleton ambiguousStateString
                                            stateBVList = getStateBitVector (alphabet firstCharInfo) <$> stateSTList
                                        in  -- showStuff = show (firstState, ambiguousStateST, ambiguousStateString, stateSTList, stateBVList)

                                            -- trace ("GQC: " <> (show ambiguousStateString) <> " " <> (show stateSTList) <> " " <> (show stateBVList))
                                            L.foldl1' (.|.) stateBVList
                            newCharacter = emptyCharacter{stateBVPrelim = (V.singleton stateBV, V.singleton stateBV, V.singleton stateBV)}
                        in  -- trace (" -> " <> (show stateBV) <> " from " <> (show totalAlphabet))
                            getQualitativeCharacters (tail inCharInfoList) (tail inStateList) (newCharacter : curCharList)
                    else
                        if firstCharType == Add
                            then
                                if firstState == ST.fromString "-1"
                                    then getQualitativeCharacters (tail inCharInfoList) (tail inStateList) (missingAdditive firstCharInfo : curCharList)
                                    else
                                        let (minRange, maxRange) = getIntRange firstState totalAlphabet
                                            newCharacter =
                                                emptyCharacter
                                                    { rangePrelim = (V.singleton (minRange, maxRange), V.singleton (minRange, maxRange), V.singleton (minRange, maxRange))
                                                    }
                                        in  if ST.length firstState > 1
                                                then -- trace ("GQC: " <> show firstState)
                                                    getQualitativeCharacters (tail inCharInfoList) (tail inStateList) (newCharacter : curCharList)
                                                else getQualitativeCharacters (tail inCharInfoList) (tail inStateList) (newCharacter : curCharList)
                            else
                                if firstCharType == Matrix
                                    then
                                        if firstState `elem` fmap ST.fromString ["?", "-"]
                                            then getQualitativeCharacters (tail inCharInfoList) (tail inStateList) (missingMatrix firstCharInfo : curCharList)
                                            else
                                                let initialMatrixVector = getInitialMatrixVector (alphabet firstCharInfo) firstState
                                                    newCharacter = emptyCharacter{matrixStatesPrelim = V.singleton initialMatrixVector}
                                                in  -- trace (show initialMatrixVector) (
                                                    -- trace ((show $ alphabet firstCharInfo) <> " " <> (ST.toString firstState)) (
                                                    -- trace ("GQC " <> (T.unpack $ name firstCharInfo) <> (show $ alphabet firstCharInfo) <> " " <> (show $ costMatrix firstCharInfo)) (
                                                    if null (costMatrix firstCharInfo)
                                                        then
                                                            errorWithoutStackTrace
                                                                ( "\n\nMatrix character input error: No cost matrix has been specified for character "
                                                                    <> T.unpack (name firstCharInfo)
                                                                    <> " perhaps file "
                                                                    <> (takeWhile (/= '#') $ T.unpack (name firstCharInfo))
                                                                    <> " character "
                                                                    <> ((tail $ dropWhile (/= '#') $ T.unpack (name firstCharInfo)))
                                                                )
                                                        else getQualitativeCharacters (tail inCharInfoList) (tail inStateList) (newCharacter : curCharList)
                                    else -- )

                                        error ("Character type " <> show firstCharType <> " not recongnized/implemented")


{- | createLeafCharacter takes rawData and Charinfo and returns CharacterData type
need to add in missing data as well
-}
createLeafCharacter ∷ [CharInfo] → [ST.ShortText] → Int → [CharacterData]
createLeafCharacter inCharInfoList rawDataList maxCharLength
    | null inCharInfoList =
        error "Null data in charInfoList createLeafCharacter"
    | null rawDataList -- missing data
        =
        getMissingValue inCharInfoList maxCharLength
    | otherwise =
        let localCharType = charType $ head inCharInfoList
            localAlphabet = alphabet $ head inCharInfoList
            isNucleotideData = isAlphabetDna localAlphabet
            isAminoAcidData = isAlphabetAminoAcid localAlphabet
        in  -- trace ("CLC: " <> (show localCharType)) (
            if localCharType `elem` sequenceCharacterTypes
                then -- in if length inCharInfoList == 1 then  -- should this be `elem` sequenceCharacterTypes
                case localCharType of
                    NucSeq → getNucleotideSequenceChar False rawDataList
                    AminoSeq → getAminoAcidSequenceChar False rawDataList
                    -- ambiguities different, and alphabet varies with character (potentially)
                    AlignedSlim →
                        if isNucleotideData
                            then getNucleotideSequenceChar True rawDataList
                            else getGeneralSequenceChar (head inCharInfoList) rawDataList
                    AlignedWide →
                        if isAminoAcidData
                            then getAminoAcidSequenceChar True rawDataList
                            else getGeneralSequenceChar (head inCharInfoList) rawDataList
                    AlignedHuge → getGeneralSequenceChar (head inCharInfoList) rawDataList
                    SlimSeq → getGeneralSequenceChar (head inCharInfoList) rawDataList
                    WideSeq → getGeneralSequenceChar (head inCharInfoList) rawDataList
                    HugeSeq → getGeneralSequenceChar (head inCharInfoList) rawDataList
                    _ → getQualitativeCharacters inCharInfoList rawDataList []
                else
                    if length inCharInfoList /= length rawDataList
                        then error "Mismatch in number of characters and character info"
                        else getQualitativeCharacters inCharInfoList rawDataList []


-- )

exists ∷ (Foldable f) ⇒ f a → Bool
exists = not . null
