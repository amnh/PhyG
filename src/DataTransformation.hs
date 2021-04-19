{- |
Module      :  DataTransformation.hs
Description :  Module with functionality to transform phylogenetic data
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


module DataTransformation ( renameData
                          , getDataTerminalNames
                          , addMissingTerminalsToInput
                          , checkDuplicatedTerminals
                          , createNaiveData
                          , createBVNames
                          ) where


import qualified Data.Text.Lazy as T
import           Types
import           Data.List
import           Data.Maybe
import qualified Data.BitVector as BV
import qualified Data.Vector    as V
import qualified Data.Text.Short as ST
import qualified Data.Hashable as H
import           Debug.Trace

-- | renameData takes a list of rename Text pairs (new name, oldName)
-- and replaces the old name with the new
renameData :: [(T.Text, T.Text)] -> RawData -> RawData
renameData newNamePairList inData =
  if null newNamePairList then inData
  else
      let terminalData =  fst inData
      in
      if null terminalData then inData
      else 
          let newTerminalData = fmap (relabelterminalData newNamePairList) terminalData
          in
          (newTerminalData, snd inData)

-- | relabelterminalData takes a list of Text pairs and the terminals with the
-- second name in the pairs is changed to the first
relabelterminalData :: [(T.Text, T.Text)] -> TermData -> TermData
relabelterminalData namePairList terminalData@(leafName, leafData) = 
     if null namePairList then terminalData
     else 
        let foundName = find ((== leafName) .snd) namePairList 
        in
        if foundName == Nothing then terminalData
        else (fst $ fromJust foundName, leafData)

-- | getDataTerminalNames takes all input data and getss full terminal list
-- and adds missing data for trerminals not in input files 
getDataTerminalNames :: [RawData] -> [T.Text]
getDataTerminalNames inDataList =
    if null inDataList then []
    else 
        sort $ nub $ fmap fst $ concat $ fmap fst inDataList

-- | addMissingTerminalsToInput dataLeafNames renamedData 
addMissingTerminalsToInput :: [T.Text] -> RawData -> RawData
addMissingTerminalsToInput dataLeafNames inData@(termDataList, charInfoList) = 
    if null dataLeafNames then (sortOn fst termDataList, charInfoList)
    else 
        let firstLeafName = head dataLeafNames
            foundLeaf = find ((== firstLeafName) .fst)  termDataList
        in
        if foundLeaf /= Nothing then addMissingTerminalsToInput (tail dataLeafNames) inData
        else addMissingTerminalsToInput (tail dataLeafNames) ((firstLeafName, []) : termDataList, charInfoList)

-- | checkDuplicatedTerminals takes list TermData and checks for repeated terminal names
checkDuplicatedTerminals :: [TermData] -> (Bool, [T.Text]) 
checkDuplicatedTerminals inData =
    if null inData then (False, []) 
    else 
        let nameList = group $ sort $ fmap fst inData
            dupList = filter ((>1).length) nameList
        in
        if null dupList then (False, [])
        else (True, fmap head dupList)

-- | joinSortFileData takes list if list of short text and merges line by line to joing leaf states
-- and sorts the result
joinSortFileData :: [[ST.ShortText]] -> [String]
joinSortFileData inFileLists =
    if ((length $ head inFileLists) == 0) then []
    else     
        let firstLeaf = sort $ ST.toString $ ST.concat $ fmap head inFileLists
        in
        firstLeaf : joinSortFileData (fmap tail inFileLists)


-- | createBVNames takes input data, sorts the raw data, hashes, sorts those to create
-- unique, label invariant (but data related so arbitrary but consistent)
-- Assumes the rawData come in sorted by the data reconciliation process
-- These used for vertex labels, caching, left/right DO issues
createBVNames :: [RawData] -> [(T.Text, BV.BV)]
createBVNames inDataList =
    let rawDataList = fmap fst inDataList
        textNameList = fmap fst $ head rawDataList
        textNameList' = fmap fst $ last rawDataList
        fileLeafCharList = fmap (fmap snd) rawDataList
        fileLeafList =  fmap (fmap ST.concat) fileLeafCharList
        leafList = reverse $ joinSortFileData fileLeafList
        leafHash = fmap H.hash leafList 
        leafHashPair = sortOn fst $ zip leafHash textNameList
        (_, leafReoderedList) = unzip leafHashPair
        leafOrder = sortOn fst $ zip leafReoderedList [0..((length textNameList) - 1)]
        (nameList, intList) = unzip leafOrder
        bv1 = BV.bitVec (length nameList) 1
        bvList = fmap (bv1 BV.<<.) (fmap (BV.bitVec (length nameList)) intList)
    in
    if textNameList /= textNameList' then error "Taxa are not properly ordered in createBVNames"
    else zip nameList bvList

-- | createNaiveData takes input RawData and transforms to "Naive" data.
-- these data are otganized into bloicks (set to input filenames initially)
-- and are bitvector coded, but are not organized by charcter type, packed ot
-- optimized in any other way (prealigned-> nonadd, Sankoff.  2 state sankoff to binary, 
-- constant charcaters skipped etc)
-- these processes take place latet
-- these data can be input to any data optimization commands and are useful
-- for data output as they haven't been reordered or transformed in any way.
-- the RawData is a list since it is organized by input file
-- the list accumulator is to avoid Vector snoc/cons O(n)
createNaiveData :: [RawData] -> [(T.Text, BV.BV)] -> [BlockData] -> ProcessedData
createNaiveData inDataList leafBitVectorNames curBlockData = 
    if null inDataList then (V.fromList $ fmap fst leafBitVectorNames, V.fromList $ reverse curBlockData)
    else 
        let firstInput@(firstData, firstCharInfo) = head inDataList
        in
        -- empty file should have been caught earlier, but avoids some head/tail errors
        if null firstCharInfo then createNaiveData (tail inDataList) leafBitVectorNames  curBlockData
        else 
            -- process data as come in--each of these should be from a single file
            -- and initially assigned to a single, unique block
            let thisBlockName = T.takeWhile (/= ':') $ name $ head firstCharInfo
                thisBlockCharInfo = V.fromList firstCharInfo
                recodedCharacters = recodeRawData (fmap snd firstData) firstCharInfo []
                thisBlockGuts = V.zip (V.fromList $ fmap snd leafBitVectorNames) recodedCharacters
                thisBlockData = (thisBlockName, thisBlockGuts, thisBlockCharInfo)
            in
            trace ("Recoding block: " ++ T.unpack thisBlockName)
            createNaiveData (tail inDataList) leafBitVectorNames  (thisBlockData : curBlockData)

-- | recodeRawData takes the ShortText representation of character states/ranges etc
-- and recodes the apporpriate fields in CharacterData (from Types) 
-- the list accumulator is to avoid Vectotr cons/snoc O(n)
-- differentiates between seqeunce type and others with char info
recodeRawData :: [[ST.ShortText]] -> [CharInfo] -> [[CharacterData]] -> V.Vector (V.Vector CharacterData)
recodeRawData inData inCharInfo curCharData =
    if null inData then V.fromList $ reverse $ fmap V.fromList curCharData
    else 
        let firstData = head inData
            firstDataRecoded = createLeafCharacter inCharInfo firstData
        in
        --trace ((show $ length inData) ++ " " ++ (show $ length firstData) ++ " " ++ (show $ length inCharInfo))
        recodeRawData (tail inData) inCharInfo (firstDataRecoded : curCharData)  

-- | createLeafCharacter takes rawData and Charinfo and returns CharcaterData type
-- need to add in missing data as well
createLeafCharacter :: [CharInfo] -> [ST.ShortText] -> [CharacterData]
createLeafCharacter inCharInfoList rawDataList =
    if null inCharInfoList then error "Null data in charInfoList createLeafCharacter"
    else if null rawDataList then
        -- missing data
        [] --getMissingValue (charType $ head inCharInfoList)
    else 
        if (length inCharInfoList == 1) && (charType (head inCharInfoList) `elem` [SmallAlphSeq, NucSeq, AminoSeq, GenSeq]) then
            trace ("Sequence character")
            []
        else 
            trace ("Non-sequence character")
            []