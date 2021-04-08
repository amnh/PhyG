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


module DataTransformation (renameData
                          ,getDataTerminalNames
                          ,addMissingTerminalsToInput
                          ,checkDuplicatedTerminals
                          ) where


import qualified Data.Text.Lazy as T
import           Types
import           Data.List
import           Data.Maybe

-- | renameData takes a list of rename Text pairs (new name, oldName)
-- and replaces the old name with the new
renameData :: [(T.Text, T.Text)] -> PhyloData -> PhyloData
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
getDataTerminalNames :: [PhyloData] -> [T.Text]
getDataTerminalNames inDataList =
    if null inDataList then []
    else 
        sort $ nub $ fmap fst $ concat $ fmap fst inDataList

-- | addMissingTerminalsToInput dataLeafNames renamedData 
addMissingTerminalsToInput :: [T.Text] -> PhyloData -> PhyloData
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

