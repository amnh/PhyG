{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE OverloadedStrings #-}

{- |
Module exposing fasta/c sequence parsing and importing functionality.
-}
module Input.FastAC (
    getFastAText,
    getFastaCharInfo,
    getFastC,
    getFastCText,
    getFastcCharInfo,
    genDiscreteDenseOfDimension,
) where

import Control.DeepSeq
import Control.Monad (when)
import Control.Monad.IO.Class (MonadIO (..))
import Data.Alphabet
import Data.Bits
import Data.Char qualified as C
import Data.Hashable
import Data.Foldable (fold)
import Data.List qualified as L
import Data.List.NonEmpty (NonEmpty (..))
import Data.List.NonEmpty qualified as NE
import Data.MetricRepresentation
import Data.MetricRepresentation qualified as MR
import Data.TCM (TCMDiagnosis (..), TCMStructure (..))
import Data.TCM qualified as TCM
import Data.TCM.Dense qualified as TCMD
import Data.Text.Lazy qualified as T
import Data.Text.Short qualified as ST
import Data.Vector qualified as V
import GeneralUtilities
import Input.DataTransformation qualified as DT
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import PHANE.Evaluation.Verbosity (Verbosity (..))
import SymMatrix qualified as S
import Types.Types

import Debug.Trace

{- | getAlphabet takse a list of short-text lists and returns alphabet as list of short-text
although with multicharacter alphabets that contain '[' or ']' this would be a problem,
its only used for single character alphabets in fasta formats.
'#' for partitions in fasta sequences
-}
getAlphabet ∷ [String] → [ST.ShortText] → [ST.ShortText]
getAlphabet curList inList =
    let notAlphElement = ST.fromString <$> ["?", "[", "]", "#"]
    in  if null inList
            then filter (`notElem` notAlphElement) $ fmap ST.fromString $ L.sort curList `L.union` ["-"]
            else
                let firstChars = fmap (: []) $ L.nub $ ST.toString $ head inList
                in  getAlphabet (firstChars `L.union` curList) (tail inList)


{- | generateDefaultMatrix takes an alphabet and generates cost matrix (assuming '-'
  in already)
-}
generateDefaultMatrix ∷ Alphabet ST.ShortText → Int → Int → Int → [[Int]]
generateDefaultMatrix inAlph rowCount indelCost substitutionCost
    | null inAlph = []
    | rowCount == length inAlph = []
    | otherwise =
        let firstPart =
                if rowCount < (length inAlph - 1)
                    then replicate rowCount substitutionCost
                    else replicate rowCount indelCost
            thirdPart =
                if rowCount < (length inAlph - 1)
                    then replicate (length inAlph - rowCount - 1 - 1) substitutionCost <> [indelCost]
                    else []
        in  (firstPart <> [0] <> thirdPart) : generateDefaultMatrix inAlph (rowCount + 1) indelCost substitutionCost


{- | getFastaCharInfo get alphabet , names etc from processed fasta data
this doesn't separate ambiguities from elements--processed later
need to read in TCM or default
only for single character element sequecnes
-}
getFastaCharInfo ∷ [TermData] → String → String → Bool → ([ST.ShortText], [[Int]], Double) → PhyG (CharInfo, [TermData])
getFastaCharInfo inData dataName dataType isPrealigned localTCM =
    if null inData
        then error "Empty inData in getFastaCharInfo"
        else
            let nucleotideAlphabet = fmap ST.fromString ["A", "C", "G", "T", "U", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N", "?", "-"]
                lAminoAcidAlphabet = fmap ST.fromString ["A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "X", "Y", "Z", "-", "?"]
                -- onlyInNucleotides = [ST.fromString "U"]
                -- onlyInAminoAcids = fmap ST.fromString ["E","F","I","L","P","Q","X","Z"]
                sequenceData = getAlphabet [] $ foldMap snd inData

                seqType
                    | dataType == "nucleotide" {- trace ("File " <> dataName <> " is nucleotide data.") -} = NucSeq
                    | dataType == "aminoacid" {- trace ("File " <> dataName <> " is aminoacid data.") -} = AminoSeq
                    | dataType == "hugeseq" {- trace ("File " <> dataName <> " is large alphabet data.") -} = HugeSeq
                    | dataType == "custom_alphabet" {- trace ("File " <> dataName <> " is large alphabet data.") -} = HugeSeq
                    | (sequenceData `L.intersect` nucleotideAlphabet == sequenceData {- trace ("Assuming file " <> dataName
                                                                                     <> " is nucleotide data. Specify `aminoacid' filetype if this is incorrect.") -}) =
                        NucSeq
                    | (sequenceData `L.intersect` lAminoAcidAlphabet == sequenceData {- trace ("Assuming file " <> dataName
                                                                                     <> " is amino acid data. Specify `nucleotide' filetype if this is incorrect.") -}) =
                        AminoSeq
                    | length sequenceData <= 8 {- trace ("File " <> dataName <> " is small alphabet data.") -} = SlimSeq
                    | length sequenceData <= 64 {- trace ("File " <> dataName <> " is wide alphabet data.") -} = WideSeq
                    | otherwise {- trace ("File " <> dataName <> " is large alphabet data.") -} = HugeSeq

                seqAlphabet = fromSymbols seqSymbols

                seqSymbols =
                    let toSymbols = fmap ST.fromString
                    in  case seqType of
                            NucSeq → toSymbols $ "A" :| ["C", "G", "T", "-"]
                            AminoSeq → toSymbols $ "A" :| ["C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-"]
                            _ → "-" :| sequenceData

                thisAlphabet =
                    case fst3 localTCM of
                        [] → seqAlphabet
                        a : as → fromSymbols $ a :| as

                -- capitalize input data if NucSeq or AminoSeq
                outData =
                    if seqType `notElem` [NucSeq, AminoSeq]
                        then inData
                        else fmap makeUpperCaseTermData inData
            in  do
                    case seqType of
                        NucSeq → logWith LogInfo ("File " <> dataName <> " is nucleotide data." <> "\n")
                        AminoSeq → logWith LogInfo ("File " <> dataName <> " is aminoacid data." <> "\n")
                        HugeSeq → logWith LogInfo ("File " <> dataName <> " is large alphabet data." <> "\n")
                        WideSeq → logWith LogInfo ("File " <> dataName <> " is wide alphabet data." <> "\n")
                        SlimSeq → logWith LogInfo ("File " <> dataName <> " is slim alphabet data." <> "\n")
                        _ → failWithPhase Parsing ("File " <> dataName <> " is of unknown data type." <> "\n")

                    processedCharInfo ← commonFastCharInfo dataName isPrealigned localTCM seqType thisAlphabet
                    pure (processedCharInfo, outData)


-- | makeUpperCaseTermData
makeUpperCaseTermData ∷ TermData → TermData
makeUpperCaseTermData (taxName, dataList) =
    let newData = fmap (fmap C.toUpper) $ fmap ST.unpack dataList
    in  (taxName, fmap ST.pack newData)


-- | commonFastCharInfo breaks out common functions between fasta and fastc parsing
commonFastCharInfo ∷ String → Bool → ([ST.ShortText], [[Int]], Double) → CharType → Alphabet ST.ShortText → PhyG CharInfo
commonFastCharInfo dataName isPrealigned localTCM seqType thisAlphabet =
    let localCostMatrix ∷ S.Matrix Int
        localCostMatrix =
            force $
                if null $ fst3 localTCM
                    then
                        let (indelCost, substitutionCost) =
                                if null $ snd3 localTCM
                                    then (1, 1)
                                    else ((head . head . snd3) localTCM, (last . head . snd3) localTCM)
                        in  S.fromLists $ generateDefaultMatrix thisAlphabet 0 indelCost substitutionCost
                    else S.fromLists $ snd3 localTCM

        localCostMatrixTransformed = transformGapLastToGapFirst' localCostMatrix

        scmForSlim ∷ ∀ {a1} {a2}. (Enum a1, Enum a2) ⇒ a1 → a2 → Word
        scmForSlim i j = localCostMatrixTransformed S.! (fromEnum i, fromEnum j)

        tcmWeightFactor = thd3 localTCM

        nilValuesTCM ∷ [[Word]]
        nilValuesTCM = [[0,0],[0,0]]

        nilTCM = snd $ TCM.fromRows nilValuesTCM

        alignedSeqType
            | not isPrealigned = seqType
            | seqType `elem` [NucSeq, SlimSeq] = AlignedSlim
            | seqType `elem` [WideSeq, AminoSeq] = AlignedWide
            | seqType == HugeSeq = AlignedHuge
            | otherwise = error "Unrecognozed data type in getFastcCharInfo"

        notNullFromTCM f = not . null $ f localTCM
        noNeedToEmitWarning = notNullFromTCM fst3 || notNullFromTCM snd3

        logMessage
            | noNeedToEmitWarning = fold [ "Processing TCM data for file : ", dataName,  "\n" ]
            | otherwise = fold
                [ "Warning: no tcm file specified for use with fasta/c file : "
                , dataName
                , ". Using default, all 1 diagonal 0 cost matrix."
                , "\n"
                ]

        getAllTCMs = 
                let dimension = force . fromIntegral $ V.length localCostMatrix
    
                    -- this 2x2 so if some Show instances are called don't get error
                    slimMetricNil = genDiscreteDenseOfDimension dimension
    
                    slimMetric = TCMD.generateDenseTransitionCostMatrix 0 dimension scmForSlim
    
                in  do  wideMetricNil <- metricRepresentation nilTCM
                        hugeMetricNil <- metricRepresentation nilTCM
                        let resultSlim = pure (1, slimMetric, wideMetricNil, hugeMetricNil)
                        let resultWide = do
                                (wideWeight, wideMetric) <- getTCMMemo (thisAlphabet, localCostMatrixTransformed)
                                pure (wideWeight, slimMetricNil, wideMetric, hugeMetricNil)
                      
                        let resultHuge = do
                                (hugeWeight, hugeMetric) <- getTCMMemo (thisAlphabet, localCostMatrixTransformed)
                                pure (hugeWeight, slimMetricNil, wideMetricNil, hugeMetric)

                        case seqType of
                            NucSeq → resultSlim
                            SlimSeq → resultSlim
                            WideSeq → resultWide
                            AminoSeq → resultWide
                            HugeSeq → resultHuge
                            _ → error $ "getFastaCharInfo: Failure proceesing the CharType: '" <> show seqType <> "'"


    in  do  if noNeedToEmitWarning
            then logWith LogInfo logMessage 
            else logWith LogWarn logMessage
            defInfo <- emptyCharInfo

            (additionalWeight, localDenseCostMatrix, localWideTCM, localHugeTCM) <- getAllTCMs

            pure $ defInfo
                    { charType = alignedSeqType
                    , activity = True
                    , weight = tcmWeightFactor * fromRational additionalWeight
                    , costMatrix = localCostMatrix
                    , slimTCM = localDenseCostMatrix
                    , wideTCM = localWideTCM
                    , hugeTCM = localHugeTCM
                    , name = T.pack (filter (/= ' ') dataName <> "#0")
                    , alphabet = thisAlphabet
                    , prealigned = isPrealigned
                    , origInfo = V.singleton (T.pack (filter (/= ' ') dataName <> "#0"), alignedSeqType, thisAlphabet)
                    }


-- | getTCMMemo creates the memoized tcm for large alphabet sequences
getTCMMemo
    ∷ ( FiniteBits b
      , Hashable b
      , Integral i
      , MonadIO m
      , NFData b
      )
    ⇒ (a, S.Matrix i)
    → m (Rational, MR.MetricRepresentation b)
getTCMMemo (_inAlphabet, inMatrix) =
    let (coefficient, tcm) = force . TCM.fromRows $ S.getFullVects inMatrix
    in  do  metric <- case tcmStructure $ TCM.diagnoseTcm tcm of
                NonAdditive → pure discreteMetric
                Additive → pure . linearNorm . toEnum $ TCM.size tcm
                _ → {-# SCC "tcmStructure_OTHER" #-} metricRepresentation tcm
            pure (coefficient, metric)


{- | getSequenceAphabet take a list of ShortText with information and accumulatiors
For both nonadditive and additve looks for [] to denote ambiguity and splits states
 if splits on spaces if there are spaces (within []) (ala fastc or multicharacter states)
 else if no spaces
   if non-additive then each symbol is split out as an alphabet element -- as in TNT
   if is additive splits on '-' to denote range
rescales (integerizes later) additive characters with decimal places to an integer type rep
for additive charcaters if states are not nummerical then throuws an error
DOES not check for partitin characters
-}
getSequenceAphabet ∷ [ST.ShortText] → [ST.ShortText] → [ST.ShortText]
getSequenceAphabet newAlph inStates =
    if null inStates
        then -- removes indel gap from alphabet if present and then (re) adds at end
        -- (filter (/= (ST.singleton '-')) $ sort $ nub newAlph) <> [ST.singleton '-']
            L.sort (L.nub newAlph) <> [ST.singleton '-']
        else
            let firstState = ST.toString $ head inStates
            in  if head firstState /= '['
                    then
                        if firstState `elem` ["?", "-"]
                            then getSequenceAphabet newAlph (tail inStates)
                            else getSequenceAphabet (head inStates : newAlph) (tail inStates)
                    else -- ambiguity

                        let newAmbigStates = fmap ST.fromString $ words $ filter (`notElem` ['[', ']']) firstState
                        in  getSequenceAphabet (newAmbigStates <> newAlph) (tail inStates)


{- | getFastcCharInfo get alphabet , names etc from processed fasta data
this doesn't separate ambiguities from elements--processed later
need to read in TCM or default
-}

-- Not correct with default alphabet and matrix now after tcm recodeing added to low for decmials I htink.
-- only for multi-character element seqeunces
getFastcCharInfo ∷ [TermData] → String → Bool → ([ST.ShortText], [[Int]], Double) → PhyG CharInfo
getFastcCharInfo inData dataName isPrealigned localTCM =
    if null inData
        then error "Empty inData in getFastcCharInfo"
        else -- if null $ fst localTCM then errorWithoutStackTrace ("Must specify a tcm file with fastc data for fie : " <> dataName)

            let symbolsFound
                    | not $ null $ fst3 localTCM = fst3 localTCM
                    | otherwise = getSequenceAphabet [] $ concatMap snd inData

                thisAlphabet = fromSymbols $ NE.fromList symbolsFound

                seqType =
                    case length thisAlphabet of
                        n | n <= 8 → SlimSeq
                        n | n <= 64 → WideSeq
                        _ → HugeSeq
            in  commonFastCharInfo dataName isPrealigned localTCM seqType thisAlphabet


{- | getFastAText processes fasta file
assumes single character alphabet
deletes '-' (unless "prealigned"), and spaces
-}
getFastAText ∷ T.Text → String → Bool → [TermData]
getFastAText fileContents' fileName isPreligned =
    if T.null fileContents'
        then errorWithoutStackTrace "\n\n'Read' command error: empty file"
        else -- removes ';' comments and spaces

            let fileContents = T.unlines $ filter (not . T.null) $ T.takeWhile (/= ';') <$> T.lines fileContents'
            in  if T.head fileContents /= '>'
                    then errorWithoutStackTrace "\n\n'Read' command error: fasta file must start with '>'"
                    else
                        let terminalSplits = T.split (== '>') fileContents
                            pairData = getRawDataPairsFastA isPreligned (tail terminalSplits)
                            (hasDupTerminals, dupList) = DT.checkDuplicatedTerminals pairData
                        in  -- tail because initial split will an empty text
                            if hasDupTerminals
                                then errorWithoutStackTrace ("\tInput file " <> fileName <> " has duplicate terminals: " <> show dupList)
                                else pairData


{- | getRawDataPairsFastA takes splits of Text and returns terminalName, Data pairs--minimal error checking
taxon nmame unitil finds ' ' , '$' or ';'
-}
getRawDataPairsFastA ∷ Bool → [T.Text] → [TermData]
getRawDataPairsFastA isPreligned inTextList =
    if null inTextList
        then []
        else
            let firstText = head inTextList
                firstName =
                    T.strip $
                        T.filter (/= '"') $
                            T.filter C.isPrint $
                                T.takeWhile (/= ' ') $
                                    T.takeWhile (/= '$') $
                                        T.takeWhile (/= ';') $
                                            head $
                                                T.lines firstText
                firstData = T.strip $ T.filter C.isPrint $ T.filter (/= ' ') $ T.concat $ tail $ T.lines firstText
                firstDataNoGaps = T.filter (/= '-') firstData
                firtDataSTList = fmap (ST.fromText . T.toStrict) (T.chunksOf 1 firstData)
                firstDataNoGapsSTList = fmap (ST.fromText . T.toStrict) (T.chunksOf 1 firstDataNoGaps)
            in  -- trace (T.unpack firstName <> "\n"  <> T.unpack firstData) (
                -- trace ("FA " <> (show firtDataSTList)) (
                if isPreligned
                    then -- trace ("GRDPF: " <> (show isPreligned))
                        (firstName, firtDataSTList) : getRawDataPairsFastA isPreligned (tail inTextList)
                    else (firstName, firstDataNoGapsSTList) : getRawDataPairsFastA isPreligned (tail inTextList)


-- )

{- | getFastC processes fasta file
assumes spaces between alphabet elements
deletes '-' (unless "prealigned")
-}
getFastC ∷ String → String → Bool → [TermData]
getFastC fileContents' fileName isPreligned =
    if null fileContents'
        then errorWithoutStackTrace "\n\n'Read' command error: empty file"
        else
            let fileContentLines = filter (not . null) $ stripString <$> lines fileContents'
            in  if null fileContentLines
                    then
                        errorWithoutStackTrace
                            ( "File "
                                <> show fileName
                                <> " is having problems reading as 'fastc'.  If this is a 'fasta' file, "
                                <> "prepend `fasta:' to the file name as in 'fasta:\"bleh.fas\"'"
                            )
                    else -- ';' comments if in terminal name are removed by getRawDataPairsFastC--otherwise leaves in there--unless its first character of line
                    --  because of latexIPA encodings using ';'(and '$')

                        let fileContents = unlines $ filter ((/= ';') . head) fileContentLines
                        in  if null fileContents
                                then
                                    errorWithoutStackTrace
                                        ( "File "
                                            <> show fileName
                                            <> " is having problems reading as 'fastc'.  If this is a 'fasta' file, "
                                            <> "prepend `fasta:' to the file name as in 'fasta:\"bleh.fas\"'"
                                        )
                                else
                                    if head fileContents /= '>'
                                        then errorWithoutStackTrace "\n\n'Read' command error: fasta file must start with '>'"
                                        else
                                            let terminalSplits = T.split (== '>') $ T.pack fileContents
                                                pairData = recodeFASTCAmbiguities fileName $ getRawDataPairsFastC isPreligned (tail terminalSplits)
                                                (hasDupTerminals, dupList) = DT.checkDuplicatedTerminals pairData
                                            in  -- tail because initial split will an empty text
                                                if hasDupTerminals
                                                    then errorWithoutStackTrace ("\tInput file " <> fileName <> " has duplicate terminals: " <> show dupList)
                                                    else pairData


{- | getFastCText processes fasta file
assumes spaces between alphabet elements
deletes '-' (unless "prealigned")
-}
getFastCText ∷ T.Text → String → Bool → [TermData]
getFastCText fileContents' fileName isPreligned =
    if T.null fileContents'
        then errorWithoutStackTrace "\n\n'Read' command error: empty file"
        else
            let fileContentLines = filter (not . T.null) $ fmap T.strip (T.lines fileContents')
            in  if null fileContentLines
                    then
                        errorWithoutStackTrace
                            ( "File "
                                <> show fileName
                                <> " is having problems reading as 'fastc'.  If this is a 'fasta' file, "
                                <> "prepend `fasta:' to the file name as in 'fasta:\"bleh.fas\"'"
                            )
                    else -- ';' comments if in terminal name are removed by getRawDataPairsFastC--otherwise leaves in there--unless its first character of line
                    --  because of latexIPA encodings using ';'(and '$')

                        let fileContents = T.unlines $ filter ((/= ';') . T.head) fileContentLines
                        in  if T.null fileContents
                                then
                                    errorWithoutStackTrace
                                        ( "File "
                                            <> show fileName
                                            <> " is having problems reading as 'fastc'.  If this is a 'fasta' file, "
                                            <> "prepend `fasta:' to the file name as in 'fasta:\"bleh.fas\"'"
                                        )
                                else
                                    if T.head fileContents /= '>'
                                        then errorWithoutStackTrace "\n\n'Read' command error: fasta file must start with '>'"
                                        else
                                            let terminalSplits = T.split (== '>') fileContents
                                                pairData = recodeFASTCAmbiguities fileName $ getRawDataPairsFastC isPreligned (tail terminalSplits)
                                                (hasDupTerminals, dupList) = DT.checkDuplicatedTerminals pairData
                                            in  -- tail because initial split will an empty text
                                                if hasDupTerminals
                                                    then errorWithoutStackTrace ("\tInput file " <> fileName <> " has duplicate terminals: " <> show dupList)
                                                    else pairData


-- | recodeFASTCAmbiguities take list of TermData and scans for ambiguous groups staring with '['' and ending with ']
recodeFASTCAmbiguities ∷ String → [TermData] → [TermData]
recodeFASTCAmbiguities fileName inData =
    if null inData
        then []
        else
            let (firstName, firstData) = head inData
                newData = concatAmbig fileName firstData
            in  (firstName, newData) : recodeFASTCAmbiguities fileName (tail inData)


{- | concatAmbig takes a list of ShortText and concatanates ambiguyous states '['X Y Z...']' into a
single Short Tex for later processing
-}
concatAmbig ∷ String → [ST.ShortText] → [ST.ShortText]
concatAmbig fileName inList =
    if null inList
        then []
        else
            let firstGroup = ST.toString $ head inList
            in  -- not ambiguity group
                -- trace (firstGroup <> show inList) (
                if null firstGroup
                    then concatAmbig fileName (tail inList)
                    else
                        if head firstGroup /= '['
                            then head inList : concatAmbig fileName (tail inList)
                            else
                                let ambiguityGroup = head inList : getRestAmbiguityGroup fileName (tail inList)
                                in  -- trace (show ambiguityGroup)
                                    ST.concat ambiguityGroup : concatAmbig fileName (drop (length ambiguityGroup) inList)


-- )

-- | getRestAmbiguityGroup takes a list of ShorText and keeps added them until one is found with ']'
getRestAmbiguityGroup ∷ String → [ST.ShortText] → [ST.ShortText]
getRestAmbiguityGroup fileName inList =
    if null inList
        then errorWithoutStackTrace ("\n\n'Read' command error: fastc file " <> fileName <> " with unterminated ambiguity specification ']'")
        else
            let firstGroup = ST.toString $ head inList
            in  if ']' `notElem` firstGroup
                    then ST.cons ' ' (head inList) : getRestAmbiguityGroup fileName (tail inList)
                    else [ST.cons ' ' $ head inList]


{- | getRawDataPairsFastC takes splits of Text and returns terminalName, Data pairs--minimal error checking
this splits on spaces in sequences
takes taxon name until encouters '' '. '$', or ';'
-}
getRawDataPairsFastC ∷ Bool → [T.Text] → [TermData]
getRawDataPairsFastC isPreligned inTextList =
    if null inTextList
        then []
        else
            let firstText = head inTextList
                firstName =
                    T.strip $
                        T.filter (/= '"') $
                            T.filter C.isPrint $
                                T.takeWhile (/= ' ') $
                                    T.takeWhile (/= '$') $
                                        T.takeWhile (/= ';') $
                                            head $
                                                T.lines firstText
                firstData = T.split (== ' ') $ T.concat $ tail $ T.lines firstText
                firstDataNoGaps = filter (/= "-") firstData
            in  -- trace (show firstData) (
                -- trace (T.unpack firstName <> "\n"  <> (T.unpack $ T.intercalate " " firstData)) (
                if isPreligned
                    then (firstName, fmap (ST.fromText . T.toStrict) firstData) : getRawDataPairsFastC isPreligned (tail inTextList)
                    else (firstName, fmap (ST.fromText . T.toStrict) firstDataNoGaps) : getRawDataPairsFastC isPreligned (tail inTextList)


-- | add to tnt
genDiscreteDenseOfDimension
    ∷ (Enum i)
    ⇒ i
    → TCMD.DenseTransitionCostMatrix
genDiscreteDenseOfDimension d =
    let n = toEnum $ fromEnum d
        r = [0 .. n - 1]
        m = [[if i == j then 0 else 1 | j ← r] | i ← r]
    in  TCMD.generateDenseTransitionCostMatrix n n . S.getCost $ V.fromList <$> V.fromList m


transformGapLastToGapFirst' ∷ S.Matrix Int → S.Matrix Word
transformGapLastToGapFirst' mat =
    let (n, _) = S.dim mat
        m = S.getFullVects mat

        f 0 0 = m S.! (n - 1, n - 1)
        f 0 j = m S.! (n - 1, j - 1)
        f i 0 = m S.! (i - 1, n - 1)
        f i j = m S.! (i - 1, j - 1)
    in  S.fromLists [[fromIntegral $ f i j | j ← [0 .. n - 1]] | i ← [0 .. n - 1]]

{-
transformGapLastToGapFirst :: TCM.TCM -> TCM.TCM
transformGapLastToGapFirst tcm =
  let n = TCM.size tcm
  in  TCM.generate n $ \case
        (0,0) -> tcm TCM.! (n - 1, n - 1)
        (0,j) -> tcm TCM.! (n - 1, j - 1)
        (i,0) -> tcm TCM.! (i - 1, n - 1)
        (i,j) -> tcm TCM.! (i - 1, j - 1)
-}
