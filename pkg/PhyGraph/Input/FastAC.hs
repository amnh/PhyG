{- |
Module      :  FastAC.hs
Description :  Module proving fasta/c sequence import functions
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

module Input.FastAC
  ( getFastA
  , getFastaCharInfo
  , getFastC
  , getFastcCharInfo
  , genDiscreteDenseOfDimension
  ) where

import           Control.DeepSeq
import           Data.Alphabet
import           Data.Bits
import qualified Data.Char                 as C
import           Data.Hashable
import qualified Data.List                 as L
import           Data.MetricRepresentation
import qualified Data.MetricRepresentation as MR
import           Data.TCM                  (TCMDiagnosis (..),
                                            TCMStructure (..))
import qualified Data.TCM                  as TCM
import qualified Data.TCM.Dense            as TCMD
import qualified Data.Text.Lazy            as T
import qualified Data.Text.Short           as ST
import qualified Data.Vector               as V
import           Debug.Trace
import           GeneralUtilities
import qualified Input.DataTransformation  as DT
import qualified SymMatrix                 as S
import           Types.Types


-- | getAlphabet takse a list of short-text lists and returns alphabet as list of short-text
-- although with multicharacter alphabets that contain '[' or ']' this would be a problem,
-- its only used for single character alphabets in fasta formats.
-- '#' for partitions in fasta sequences
getAlphabet :: [String] -> [ST.ShortText] -> [ST.ShortText]
getAlphabet curList inList =
    let notAlphElement = fmap ST.fromString ["?", "[", "]", "#"]
    in
    if null inList then filter (`notElem` notAlphElement) $ fmap ST.fromString $ L.sort curList `L.union` ["-"]
    else
        let firstChars = fmap (:[]) $ L.nub $ ST.toString $ head inList
        in  getAlphabet (firstChars `L.union` curList) (tail inList)


-- | generateDefaultMatrix takes an alphabet and generates cost matrix (assuming '-'
--   in already)
generateDefaultMatrix :: Alphabet ST.ShortText -> Int -> Int -> Int -> [[Int]]
generateDefaultMatrix inAlph rowCount indelCost substitutionCost =
    if null inAlph then []
    else if rowCount == length inAlph then []
    else
        let firstPart = if rowCount < ((length inAlph) - 1) then replicate rowCount substitutionCost
                        else replicate rowCount indelCost
            thirdPart = if rowCount < ((length inAlph) - 1) then (replicate ((length inAlph) - rowCount - 1 - 1) substitutionCost) ++ [indelCost]
                        else []
        in
        (firstPart ++ [0] ++ thirdPart) : generateDefaultMatrix inAlph (rowCount + 1) indelCost substitutionCost

-- | getFastaCharInfo get alphabet , names etc from processed fasta data
-- this doesn't separate ambiguities from elements--processed later
-- need to read in TCM or default
getFastaCharInfo :: [TermData] -> String -> String -> Bool -> ([ST.ShortText], [[Int]], Double) -> CharInfo
getFastaCharInfo inData dataName dataType isPrealigned localTCM =
    if null inData then error "Empty inData in getFastaCharInfo"
    else
        let nucleotideAlphabet = fmap ST.fromString ["A","C","G","T","U","R","Y","S","W","K","M","B","D","H","V","N","?","-"]
            lAminoAcidAlphabet  = fmap ST.fromString ["A","B","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","X","Y","Z", "-","?"]
            --onlyInNucleotides = [ST.fromString "U"]
            --onlyInAminoAcids = fmap ST.fromString ["E","F","I","L","P","Q","X","Z"]
            sequenceData = getAlphabet [] $ foldMap snd inData

            seqType
              | dataType == "nucleotide" = trace ("File " ++ dataName ++ " is nucleotide data.")  NucSeq
              | dataType == "aminoacid" = trace ("File " ++ dataName ++ " is aminoacid data.") AminoSeq
              | dataType == "custom_alphabet" = trace ("File " ++ dataName ++ " is large alphabet data.") HugeSeq
              | (sequenceData `L.intersect` nucleotideAlphabet == sequenceData) = trace ("Assuming file " ++ dataName
                ++ " is nucleotide data. Specify `aminoacid' filetype if this is incorrect.") NucSeq
              | (sequenceData `L.intersect` lAminoAcidAlphabet == sequenceData) = trace ("Assuming file " ++ dataName
                ++ " is amino acid data. Specify `nucleotide' filetype if this is incorrect.") AminoSeq
              | length sequenceData <=  8 = trace ("File " ++ dataName ++ " is small alphabet data.") SlimSeq
              | length sequenceData <= 64 = trace ("File " ++ dataName ++ " is wide alphabet data.") WideSeq
              | otherwise = trace ("File " ++ dataName ++ " is large alphabet data.") HugeSeq

            seqAlphabet = fromSymbols seqSymbols

            seqSymbols =
              let toSymbols = fmap ST.fromString
              in  case seqType of
                    NucSeq   -> toSymbols ["A","C","G","T","-"]
                    AminoSeq -> toSymbols ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y", "-"]
                    _        -> sequenceData

            localCostMatrix = if null $ fst3 localTCM then
                                let (indelCost, substitutionCost) = if null $ snd3 localTCM then (1,1)
                                                                    else ((head . head . snd3) localTCM,  (last . head . snd3) localTCM)
                                in
                                S.fromLists $ generateDefaultMatrix seqAlphabet 0 indelCost substitutionCost
                              else S.fromLists $ snd3 localTCM

            tcmDense = TCMD.generateDenseTransitionCostMatrix 0 (fromIntegral $ V.length localCostMatrix) (getCost localCostMatrix)
            -- not sure of this
            tcmNaught = genDiscreteDenseOfDimension (length sequenceData)
            localDenseCostMatrix = if seqType `elem` [NucSeq, SlimSeq] then tcmDense
                                   else tcmNaught

            (wideWeightFactor, localWideTCM)
              | seqType `elem` [WideSeq, AminoSeq] = getTCMMemo (thisAlphabet, localCostMatrix)
              | otherwise                          = metricRepresentation <$> TCM.fromRows [[0::Word, 0::Word],[0::Word, 0::Word]] -- this 2x2 so if some Show instances are called don't get error

            (hugeWeightFactor, localHugeTCM)
              | seqType == HugeSeq           = getTCMMemo (thisAlphabet, localCostMatrix)
              | otherwise                    = metricRepresentation <$> TCM.fromRows [[0::Word, 0::Word],[0::Word, 0::Word]] -- this 2x2 so if some Show instances are called don't get error

            tcmWeightFactor = thd3 localTCM
            thisAlphabet =
              case fst3 localTCM of
                a | null a -> seqAlphabet
                a          -> fromSymbols a

            alignedSeqType = if not isPrealigned then seqType
                             else
                                if seqType `elem` [NucSeq, SlimSeq] then AlignedSlim
                                else if seqType `elem` [WideSeq, AminoSeq] then AlignedWide
                                else if seqType == HugeSeq then AlignedHuge
                                else error "Unrecognozed data type in getFastaCharInfo"

            defaultSeqCharInfo = CharInfo {
                                               charType = alignedSeqType
                                             , activity = True
                                             , weight = tcmWeightFactor *
                                                        if seqType == HugeSeq
                                                        then fromRational hugeWeightFactor
                                                        else if seqType `elem` [WideSeq, AminoSeq]
                                                        then fromRational wideWeightFactor
                                                        else 1
                                             , costMatrix = localCostMatrix
                                             , slimTCM = localDenseCostMatrix
                                             , wideTCM = localWideTCM
                                             , hugeTCM = localHugeTCM
                                             , name = T.pack (filter (/= ' ') dataName <> "#0")
                                             , alphabet = thisAlphabet
                                             , prealigned = isPrealigned
                                             , origInfo = V.singleton (T.pack (filter (/= ' ') dataName <> "#0"), alignedSeqType, thisAlphabet)
                                             }
        in
        -- trace ("FASTCINFO:" ++ (show $ charType defaultSeqCharInfo)) (
        if (null . fst3) localTCM && (null . snd3) localTCM then 
            trace ("Warning: no tcm file specified for use with fasta file : " ++ dataName ++ ". Using default, all 1 diagonal 0 cost matrix.") 
            defaultSeqCharInfo
        else 
            trace ("Processing TCM data for file : "  ++ dataName) 
            defaultSeqCharInfo
        -- )


-- | getCost is helper function for generartion for a dense TCM
getCost :: S.Matrix Int -> Word -> Word -> Word
getCost localCM i j =
    let x = S.getFullVects localCM
    in  toEnum $ (x V.! fromEnum i) V.! fromEnum j


-- | getTCMMemo creates the memoized tcm for large alphabet sequences
getTCMMemo
  :: ( FiniteBits b
     , Hashable b
     , NFData b
     )
  => (a, S.Matrix Int)
  -> (Rational, MR.MetricRepresentation b)
getTCMMemo (_inAlphabet, inMatrix) =
    let (coefficient, tcm) = TCM.fromRows $ S.getFullVects inMatrix
        metric = case tcmStructure $ TCM.diagnoseTcm tcm of
                   NonAdditive -> discreteMetric
                   Additive    -> linearNorm . toEnum $ TCM.size tcm
                   _           -> metricRepresentation tcm
    in (coefficient, metric)


-- | getSequenceAphabet take a list of ShortText with inform ation and accumulatiors
-- For both nonadditive and additve looks for [] to denote ambiguity and splits states
--  if splits on spaces if there are spaces (within []) (ala fastc or multicharacter states)
--  else if no spaces
--    if non-additive then each symbol is split out as an alphabet element -- as in TNT
--    if is additive splits on '-' to denote range
-- rescales (integerizes later) additive characters with decimal places to an integer type rep
-- for additive charcaters if states are not nummerical then throuws an error
getSequenceAphabet :: [ST.ShortText]  -> [ST.ShortText] -> [ST.ShortText]
getSequenceAphabet newAlph inStates =
    if null inStates then
        -- removes indel gap from alphabet if present and then (re) adds at end
        -- (filter (/= (ST.singleton '-')) $ sort $ nub newAlph) ++ [ST.singleton '-']
        L.sort (L.nub newAlph) ++ [ST.singleton '-']
    else
        let firstState = ST.toString $ head inStates
        in
        if head firstState /= '['  then
            if firstState `elem` ["?","-"] then getSequenceAphabet  newAlph (tail inStates)
            else getSequenceAphabet (head inStates : newAlph) (tail inStates)
        else -- ambiguity
            let newAmbigStates  = fmap ST.fromString $ words $ filter (`notElem` ['[',']']) firstState
            in
            getSequenceAphabet (newAmbigStates ++ newAlph) (tail inStates)


-- | getFastcCharInfo get alphabet , names etc from processed fasta data
-- this doesn't separate ambiguities from elements--processed later
-- need to read in TCM or default

--Not correct with default alphabet and matrix now after tcm recodeing added to loow for decmials I htink.
getFastcCharInfo :: [TermData] -> String -> Bool -> ([ST.ShortText], [[Int]], Double) -> CharInfo
getFastcCharInfo inData dataName isPrealigned localTCM =
    if null inData then error "Empty inData in getFastcCharInfo"
    else
        --if null $ fst localTCM then errorWithoutStackTrace ("Must specify a tcm file with fastc data for fie : " ++ dataName)
        let thisAlphabet = fromSymbols symbolsFound

            symbolsFound
              | not $ null $ fst3 localTCM = fst3 localTCM
              | otherwise = getSequenceAphabet [] $ concatMap snd inData

            inMatrix
              | not $ null $ fst3 localTCM = S.fromLists $ snd3 localTCM
              | otherwise = let (indelCost, substitutionCost) = if null $ snd3 localTCM then (1,1)
                                                                else ((head . head . snd3) localTCM,  (last . head . snd3) localTCM)
                            in
                            S.fromLists $ generateDefaultMatrix thisAlphabet 0 indelCost substitutionCost

            tcmWeightFactor = thd3 localTCM
            tcmDense = TCMD.generateDenseTransitionCostMatrix 0 (fromIntegral $ V.length inMatrix) (getCost inMatrix)

            -- not sure of this
            tcmNaught = genDiscreteDenseOfDimension (length thisAlphabet)
            localDenseCostMatrix = if length thisAlphabet < 9  then tcmDense
                                   else tcmNaught
            seqType =
                case length thisAlphabet of
                  n | n <=  8 -> SlimSeq
                  n | n <= 64 -> WideSeq
                  _           -> HugeSeq

            (wideWeightFactor, localWideTCM)
              | seqType `elem` [WideSeq, AminoSeq] = getTCMMemo (thisAlphabet, inMatrix)
              | otherwise                          = metricRepresentation <$> TCM.fromRows [[0::Word]]

            (hugeWeightFactor, localHugeTCM)
              | seqType == HugeSeq           = getTCMMemo (thisAlphabet, inMatrix)
              | otherwise                          = metricRepresentation <$> TCM.fromRows [[0::Word]]

            alignedSeqType = if not isPrealigned then seqType
                             else
                                if seqType `elem` [NucSeq, SlimSeq] then AlignedSlim
                                else if seqType `elem` [WideSeq, AminoSeq] then AlignedWide
                                else if seqType == HugeSeq then AlignedHuge
                                else error "Unrecognozed data type in getFastaCharInfo"

            defaultHugeSeqCharInfo = CharInfo {
                                       charType = alignedSeqType
                                     , activity = True
                                     , weight = tcmWeightFactor *
                                                if   seqType == HugeSeq
                                                then fromRational hugeWeightFactor
                                                else if seqType `elem` [WideSeq, AminoSeq]
                                                then fromRational wideWeightFactor
                                                else 1
                                     , costMatrix = inMatrix
                                     , slimTCM = localDenseCostMatrix
                                     , wideTCM = localWideTCM
                                     , hugeTCM = localHugeTCM
                                     , name = T.pack (filter (/= ' ') dataName ++ "#0")
                                     , alphabet = thisAlphabet
                                     , prealigned = isPrealigned
                                     , origInfo = V.singleton (T.pack (filter (/= ' ') dataName ++ "#0"), alignedSeqType, thisAlphabet)
                                     }
        in
        --trace ("FCI " ++ (show $ length thisAlphabet) ++ " alpha size" ++ show thisAlphabet) (
        if null (fst3 localTCM) then trace ("Warning: no tcm file specified for use with fastc file : " ++ dataName ++ ". Using default, all 1 diagonal 0 cost matrix.") defaultHugeSeqCharInfo
        else defaultHugeSeqCharInfo
        --)

-- | getSequenceAphabet takes a list of ShortText and returns the alp[habet and adds '-' if not present

-- | getFastA processes fasta file
-- assumes single character alphabet
-- deletes '-' (unless "prealigned"), and spaces
getFastA :: String -> String -> Bool-> [TermData]
getFastA fileContents' fileName isPreligned=
    if null fileContents' then errorWithoutStackTrace "\n\n'Read' command error: empty file"
    else
        -- removes ';' comments
        let fileContents =  unlines $ filter (not.null) $ takeWhile (/= ';') <$> lines fileContents'
        in
        if head fileContents /= '>' then errorWithoutStackTrace "\n\n'Read' command error: fasta file must start with '>'"
        else
            let terminalSplits = T.split (=='>') $ T.pack fileContents
                pairData =  getRawDataPairsFastA isPreligned (tail terminalSplits)
                (hasDupTerminals, dupList) = DT.checkDuplicatedTerminals pairData
            in
            -- tail because initial split will an empty text
            if hasDupTerminals then errorWithoutStackTrace ("\tInput file " ++ fileName ++ " has duplicate terminals: " ++ show dupList)
            else pairData

-- | getRawDataPairsFastA takes splits of Text and returns terminalName, Data pairs--minimal error checking
getRawDataPairsFastA :: Bool -> [T.Text] -> [TermData]
getRawDataPairsFastA isPreligned inTextList =
    if null inTextList then []
    else
        let firstText = head inTextList
            firstName = T.strip $ T.filter (/= '"') $ T.filter C.isPrint $ T.takeWhile (/= '$') $ T.takeWhile (/= ';') $ head $ T.lines firstText
            firstData = T.strip $ T.filter C.isPrint $ T.filter (/= ' ') $ T.toUpper $ T.concat $ tail $ T.lines firstText
            firstDataNoGaps = T.filter (/= '-') firstData
            firtDataSTList = fmap (ST.fromText . T.toStrict) (T.chunksOf 1 firstData)
            firstDataNoGapsSTList = fmap (ST.fromText . T.toStrict) (T.chunksOf 1 firstDataNoGaps)
        in
        --trace (T.unpack firstName ++ "\n"  ++ T.unpack firstData) (
        -- trace ("FA " ++ (show firtDataSTList)) (
        if isPreligned then (firstName, firtDataSTList) : getRawDataPairsFastA isPreligned (tail inTextList)
        else (firstName, firstDataNoGapsSTList) : getRawDataPairsFastA isPreligned (tail inTextList)
        -- )

-- | getFastC processes fasta file
-- assumes spaces between alphabet elements
-- deletes '-' (unless "prealigned")
-- NEED TO ADD AMBIGUITY
getFastC :: String -> String -> Bool -> [TermData]
getFastC fileContents' fileName isPreligned =
    if null fileContents' then errorWithoutStackTrace "\n\n'Read' command error: empty file"
    else
        let fileContentLines = filter (not.null) $ stripString <$> lines fileContents'
        in
        if null fileContentLines then errorWithoutStackTrace ("File " ++ show fileName ++ " is having problems reading as 'fastc'.  If this is a 'fasta' file, "
            ++ "prepend `fasta:' to the file name as in 'fasta:\"bleh.fas\"'")
        -- ';' comments if in terminal name are removed by getRawDataPairsFastC--otherwise leaves in there--unless its first character of line
        --  because of latexIPA encodings using ';'(and '$')
        else
            let fileContents = unlines $ filter ((/=';').head) fileContentLines
            in
            if null fileContents then errorWithoutStackTrace ("File " ++ show fileName ++ " is having problems reading as 'fastc'.  If this is a 'fasta' file, "
                ++ "prepend `fasta:' to the file name as in 'fasta:\"bleh.fas\"'")
            else if head fileContents /= '>' then errorWithoutStackTrace "\n\n'Read' command error: fasta file must start with '>'"
            else
                let terminalSplits = T.split (=='>') $ T.pack fileContents
                    pairData = recodeFASTCAmbiguities fileName $ getRawDataPairsFastC isPreligned (tail terminalSplits)
                    (hasDupTerminals, dupList) = DT.checkDuplicatedTerminals pairData
                in
                -- tail because initial split will an empty text
                if hasDupTerminals then errorWithoutStackTrace ("\tInput file " ++ fileName ++ " has duplicate terminals: " ++ show dupList)
                else pairData

-- | recodeFASTCAmbiguities take list of TermData and scans for ambiguous groups staring with '['' and ending with ']
recodeFASTCAmbiguities :: String -> [TermData] -> [TermData]
recodeFASTCAmbiguities fileName inData =
    if null inData then []
    else
        let (firstName, firstData) = head inData
            newData = concatAmbig fileName firstData
        in
        (firstName, newData) : recodeFASTCAmbiguities fileName (tail inData)

-- | concatAmbig takes a list of ShortText and concatanates ambiguyous states '['X Y Z...']' into a
-- single Short Tex for later processing
concatAmbig :: String -> [ST.ShortText] -> [ST.ShortText]
concatAmbig fileName inList =
    if null inList then []
    else
        let firstGroup = ST.toString $ head inList
        in
        -- not ambiguity group
        -- trace (firstGroup ++ show inList) (
        if null firstGroup then concatAmbig fileName (tail inList)
        else if head firstGroup /= '[' then head inList : concatAmbig fileName (tail inList)
        else
            let ambiguityGroup = head inList : getRestAmbiguityGroup fileName (tail inList)
            in
            --trace (show ambiguityGroup)
            ST.concat ambiguityGroup : concatAmbig fileName (drop (length ambiguityGroup) inList)
            --)

-- | getRestAmbiguityGroup takes a list of ShorText and keeps added them until one is found with ']'
getRestAmbiguityGroup :: String -> [ST.ShortText] -> [ST.ShortText]
getRestAmbiguityGroup fileName inList =
    if null inList then errorWithoutStackTrace ("\n\n'Read' command error: fastc file " ++ fileName ++ " with unterminated ambiguity specification ']'")
    else
        let firstGroup = ST.toString $ head inList
        in
        if ']' `notElem` firstGroup then ST.cons ' ' (head inList) : getRestAmbiguityGroup fileName (tail inList)
        else [ST.cons ' ' $ head inList]

-- | getRawDataPairsFastA takes splits of Text and returns terminalName, Data pairs--minimal error checking
-- this splits on spaces in sequences
getRawDataPairsFastC :: Bool -> [T.Text] -> [TermData]
getRawDataPairsFastC isPreligned inTextList =
    if null inTextList then []
    else
        let firstText = head inTextList
            firstName = T.strip $ T.filter (/= '"') $ T.filter C.isPrint $ T.takeWhile (/= '$') $ T.takeWhile (/= ';') $ head $ T.lines firstText
            firstData = T.split (== ' ') $ T.concat $ tail $ T.lines firstText
            firstDataNoGaps = filter (/= T.pack "-") firstData
        in
        --trace (show firstData) (
        -- trace (T.unpack firstName ++ "\n"  ++ (T.unpack $ T.intercalate (T.pack " ") firstData)) (
        if isPreligned then (firstName, fmap (ST.fromText . T.toStrict) firstData) : getRawDataPairsFastC isPreligned (tail inTextList)
        else (firstName, fmap (ST.fromText . T.toStrict) firstDataNoGaps) : getRawDataPairsFastC isPreligned (tail inTextList)

-- | add to tnt
genDiscreteDenseOfDimension
  :: Enum i
  => i
  -> TCMD.DenseTransitionCostMatrix
genDiscreteDenseOfDimension d =
  let n = toEnum $ fromEnum d
      r = [0 .. n - 1]
      m = [ [ if i==j then 0 else 1 | j <- r] | i <- r]
  in  TCMD.generateDenseTransitionCostMatrix n n . getCost $ V.fromList <$> V.fromList m
