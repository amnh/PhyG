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

{--
TODO:

  Parallelize  median2Vect
--}

module FastAC (   getFastA
                , getFastaCharInfo
                , getFastC
                , getFastcCharInfo
                ) where

import           Types
import           Debug.Trace
import qualified SymMatrix as S 
import qualified Data.Text.Lazy  as T
import qualified Data.Text.Short as ST
import qualified DataTransformation as DT
import           Data.List
import qualified Data.TCM as TCM
import qualified Data.TCM.Dense as TCMD
import qualified Data.Vector as V
import qualified SymMatrix as SM
import qualified Data.MetricRepresentation as MR
import qualified Bio.Character.Encodable.Dynamic.AmbiguityGroup as AG
import Data.MetricRepresentation
import qualified Data.TCM.Memoized.FFI as TCMFFI
import qualified Data.TCM.Memoized as TCMM
 




-- | getAlphabet takse a list of short-text lists and returns alphabet as list of short-text
-- filters out '?' '[' and ']' adds in '-' for indel Gap
getAlphabet :: [String] -> [ST.ShortText] -> [ST.ShortText] 
getAlphabet curList inList =
    let notAlphElement = fmap ST.fromString ["?", "[", "]"]
    in
    if null inList then 
        filter (`notElem` notAlphElement) $ fmap ST.fromString $ (sort curList) `union` ["-"]
    else 
        let firstChars = fmap (:[]) $ nub $ ST.toString $ head inList 
        in
        getAlphabet (firstChars `union` curList) (tail inList)


-- | generateDefaultMatrix takes an alphabet and generates cost matrix (assuming '-'
--   in already)
generateDefaultMatrix :: [ST.ShortText] -> Int -> [[Int]]
generateDefaultMatrix inAlph rowCount =
    if null inAlph then []
    else if rowCount == length inAlph then []
    else 
        let firstPart = replicate rowCount 1
            thirdPart = replicate ((length inAlph) - rowCount - 1) 1
        in
        (firstPart ++ [0] ++ thirdPart) : generateDefaultMatrix inAlph (rowCount + 1)

-- | getFastaCharInfo get alphabet , names etc from processed fasta data
-- this doesn't separate ambiguities from elements--processed later
-- need to read in TCM or default
getFastaCharInfo :: [TermData] -> String -> String -> Bool -> ([ST.ShortText], [[Int]]) -> CharInfo
getFastaCharInfo inData dataName dataType isPrealigned localTCM = 
    if null inData then error "Empty inData in getFastaCharInfo"
    else 
        let nucleotideAlphabet = fmap ST.fromString ["A","C","G","T","U","R","Y","S","W","K","M","B","D","H","V","N","?","-"]
            aminoAcidAlphabet  = fmap ST.fromString ["A","B","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","X","Y","Z", "-","?"]
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
            seqAlphabet = if seqType == NucSeq then fmap ST.fromString ["A","C","G","T","-"]
                       else if seqType == AminoSeq then fmap ST.fromString ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y", "-"]
                       else if seqType == Binary then fmap ST.fromString ["0","1"]
                       else fst localTCM
            localCostMatrix = if localTCM == ([],[]) then S.fromLists $ generateDefaultMatrix seqAlphabet 0
                              else S.fromLists $ snd localTCM
            tcmDense = TCMD.generateDenseTransitionCostMatrix 0 (fromIntegral $ V.length localCostMatrix) (getCost localCostMatrix)
            -- not sure of this
            tcmNaught = TCMD.generateDenseTransitionCostMatrix 0 0 (getCost V.empty)
            localDenseCostMatrix = if seqType == NucSeq || seqType == SmallAlphSeq then tcmDense
                                   else tcmNaught

            (weightFactor, localMemoTCM) = if seqType == GenSeq then getTCMMemo localTCM
                           else getTCMMemo ([],[])

            defaultGenSeqCharInfo = CharInfo {
                                       charType = seqType
                                     , activity = True
                                     , weight = if seqType == GenSeq then weightFactor
                                                else 1.0
                                     , costMatrix = localCostMatrix
                                     , denseTCM = localDenseCostMatrix
                                     , memoTCM = localMemoTCM
                                     , name = T.pack ((filter (/= ' ') dataName) ++ ":0")
                                     , alphabet = if localTCM == ([],[]) then seqAlphabet
                                                  else fst localTCM
                                     , prealigned = isPrealigned
                                     }
        in
        trace ("Warning: no tcm file specified for use with fasta file : " ++ dataName ++ ". Using default, all 1 diagonal 0 cost matrix.")
        defaultGenSeqCharInfo
    
-- | getCost is helper function for generartion for a dense TCM
getCost :: S.Matrix Int -> Word -> Word -> Word     
getCost localCM i j = 
    let x = SM.getFullVects localCM
    in  toEnum $ (x V.! fromEnum i) V.! fromEnum j


-- | getTCMMemo creates the memoized tcm for large alphabet sequences
getTCMMemo :: ([ST.ShortText], [[Int]]) -> (Double, MR.MetricRepresentation (AG.AmbiguityGroup -> AG.AmbiguityGroup -> (AG.AmbiguityGroup, Word))) 
getTCMMemo (inAlphabet, inMatrix) =
        let (lweight, tcm) = TCM.fromRows $ SM.getFullVects $ SM.fromLists inMatrix
            sigma i j       = toEnum . fromEnum $ tcm TCM.! (fromEnum i, fromEnum j)
            memoMatrixValue = TCMM.generateMemoizedTransitionCostMatrix (toEnum $ length inAlphabet) sigma
        in  (fromRational lweight, ExplicitLayout tcm (TCMFFI.getMedianAndCost2D memoMatrixValue))
        


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
        (filter (/= (ST.singleton '-')) $ sort $ nub newAlph) ++ [ST.singleton '-']
    else
        let firstState = ST.toString $ head inStates
        in
        if (head firstState) /= '['  then 
            if (firstState `elem` ["?"]) then getSequenceAphabet  newAlph (tail inStates)
            else getSequenceAphabet ((head inStates) : newAlph) (tail inStates)
        else -- ambiguity
            let newAmbigStates  = fmap ST.fromString $ words $ filter (`notElem` ['[',']']) firstState
            in
            getSequenceAphabet (newAmbigStates ++ newAlph) (tail inStates) 




-- | getFastcCharInfo get alphabet , names etc from processed fasta data
-- this doesn't separate ambiguities from elements--processed later
-- need to read in TCM or default
getFastcCharInfo :: [TermData] -> String -> Bool -> ([ST.ShortText], [[Int]]) -> CharInfo
getFastcCharInfo inData dataName isPrealigned localTCM = 
    if null inData then error "Empty inData in getFastcCharInfo"
    else 
        --if null $ fst localTCM then errorWithoutStackTrace ("Must specify a tcm file with fastc data for fie : " ++ dataName)
        let thisAlphabet = if (not $ null $ fst localTCM) then fst localTCM
                           else getSequenceAphabet [] $ concat $ fmap snd inData
            inMatrix = if (not $ null $ fst localTCM) then S.fromLists $ snd localTCM
                       else S.fromLists $ generateDefaultMatrix thisAlphabet 0
            tcmDense = TCMD.generateDenseTransitionCostMatrix 0 (fromIntegral $ V.length inMatrix) (getCost inMatrix)

            -- not sure of this
            tcmNaught = TCMD.generateDenseTransitionCostMatrix 0 0 (getCost V.empty)
            localDenseCostMatrix = if (length $ thisAlphabet) < 9  then tcmDense
                                   else tcmNaught

            (weightFactor, localMemoTCM) = if (length $ thisAlphabet) > 8 then getTCMMemo localTCM
                           else getTCMMemo ([],[])

            defaultGenSeqCharInfo = CharInfo {
                                       charType = if (length $ thisAlphabet) < 9 then SmallAlphSeq
                                                     else GenSeq
                                     , activity = True
                                     , weight = if (length $ thisAlphabet) > 8 then weightFactor
                                                else 1.0
                                     , costMatrix = inMatrix
                                     , denseTCM = localDenseCostMatrix
                                     , memoTCM = localMemoTCM
                                     , name = T.pack ((filter (/= ' ') dataName) ++ ":0")
                                     , alphabet = thisAlphabet
                                     , prealigned = isPrealigned
                                     }
        in
        trace ("Warning: no tcm file specified for use with fastc file : " ++ dataName ++ ". Using default, all 1 diagonal 0 cost matrix.")
        defaultGenSeqCharInfo

-- | getSequenceAphabet takes a list of ShortText and returns the alp[habet and adds '-' if not present  

-- | getFastA processes fasta file 
-- assumes single character alphabet
-- deletes '-' (unless "prealigned"), and spaces
getFastA :: String -> String -> String -> [TermData] 
getFastA modifier fileContents' fileName  =
    if null fileContents' then errorWithoutStackTrace ("\n\n'Read' command error: empty file")
    else 
        -- removes ';' comments   
        let fileContents =  unlines $ filter (not.null) $ fmap (takeWhile (/= ';')) $ lines fileContents'
        in 
        if (head fileContents) /= '>' then errorWithoutStackTrace ("\n\n'Read' command error: fasta file must start with '>'")
        else 
            let terminalSplits = T.split (=='>') $ T.pack fileContents 
                pairData =  getRawDataPairsFastA modifier (tail terminalSplits)
                (hasDupTerminals, dupList) = DT.checkDuplicatedTerminals pairData
            in
            -- tail because initial split will an empty text
            if hasDupTerminals then errorWithoutStackTrace ("\tInput file " ++ fileName ++ " has duplicate terminals: " ++ show dupList)
            else pairData 
        
       
        

-- | getRawDataPairsFastA takes splits of Text and returns terminalName, Data pairs--minimal error checking
getRawDataPairsFastA :: String -> [T.Text] -> [TermData]
getRawDataPairsFastA modifier inTextList =
    if null inTextList then []
    else 
        let firstText = head inTextList
            firstName = T.takeWhile (/= '$') $ T.takeWhile (/= ';') $ head $ T.lines firstText
            firstData = T.filter (/= ' ') $ T.toUpper $ T.concat $ tail $ T.lines firstText
            firstDataNoGaps = T.filter (/= '-') firstData
            firtDataSTList = fmap ST.fromText $ fmap T.toStrict $ T.chunksOf 1 firstData
            firstDataNoGapsSTList = fmap ST.fromText $ fmap T.toStrict $ T.chunksOf 1 firstDataNoGaps
        in
        --trace (T.unpack firstName ++ "\n"  ++ T.unpack firstData) (
        --trace ("FA " ++ (show $ length firstDataNoGapsSTList)) (
        if modifier == "prealigned" then (firstName, firtDataSTList) : getRawDataPairsFastA modifier (tail inTextList)
        else (firstName, firstDataNoGapsSTList) : getRawDataPairsFastA modifier (tail inTextList)
        --)
        
-- | getFastC processes fasta file 
-- assumes spaces between alphabet elements
-- deletes '-' (unless "prealigned")
-- NEED TO ADD AMBIGUITY
getFastC :: String -> String -> String -> [TermData] 
getFastC modifier fileContents' fileName =
    if null fileContents' then errorWithoutStackTrace ("\n\n'Read' command error: empty file")
    else 
        -- removes ';' comments   
        let fileContents =  unlines $ filter (not.null) $ fmap (takeWhile (/= ';')) $ lines fileContents'
        in 
        if (head fileContents) /= '>' then errorWithoutStackTrace ("\n\n'Read' command error: fasta file must start with '>'")
        else 
            let terminalSplits = T.split (=='>') $ T.pack fileContents 
                pairData = recodeFASTCAmbiguities fileName $ getRawDataPairsFastC modifier (tail terminalSplits)
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
        else if head firstGroup /= '[' then (head inList) : concatAmbig fileName (tail inList)
        else 
            let ambiguityGroup = (head inList) : getRestAmbiguityGroup fileName (tail inList)
            in
            --trace (show ambiguityGroup) 
            (ST.concat ambiguityGroup) : concatAmbig fileName (drop (length ambiguityGroup) inList)
            --)

-- | getRestAmbiguityGroup takes a list of ShorText and keeps added them until one is found with ']'
getRestAmbiguityGroup :: String -> [ST.ShortText] -> [ST.ShortText]
getRestAmbiguityGroup fileName inList = 
    if null inList then errorWithoutStackTrace ("\n\n'Read' command error: fastc file " ++ fileName ++ " with unterminated ambiguity specification ']'")
    else 
        let firstGroup = ST.toString $ head inList
        in
        if ']' `notElem` firstGroup then (ST.cons ' ' $ head inList) : getRestAmbiguityGroup fileName (tail inList)
        else [ST.cons ' ' $ head inList]

-- | getRawDataPairsFastA takes splits of Text and returns terminalName, Data pairs--minimal error checking
-- this splits on spaces in sequences
getRawDataPairsFastC :: String -> [T.Text] -> [TermData]
getRawDataPairsFastC modifier inTextList =
    if null inTextList then []
    else 
        let firstText = head inTextList
            firstName = T.takeWhile (/= '$') $ T.takeWhile (/= ';') $ head $ T.lines firstText
            firstData = T.split (== ' ') $ T.concat $ tail $ T.lines firstText
            firstDataNoGaps = fmap (T.filter (/= '-')) firstData
        in
        -- trace (T.unpack firstName ++ "\n"  ++ (T.unpack $ T.intercalate (T.pack " ") firstData)) (
        if modifier == "prealigned" then (firstName, fmap ST.fromText  $ fmap T.toStrict firstData) : getRawDataPairsFastC modifier (tail inTextList)
        else (firstName, fmap ST.fromText $ fmap T.toStrict firstDataNoGaps) : getRawDataPairsFastC modifier (tail inTextList)
        --
        
