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


module Input.Reorganize
  ( combineDataByType
  , reBlockData
  , removeConstantCharactersPrealigned
  , removeConstantCharsPrealigned
  , optimizePrealignedData
  , convert2BV
  , getRecodingType
  ) where

import           Control.Parallel.Strategies
import           Data.Alphabet
import           Data.Bits
import qualified Data.BitVector.LittleEndian as BV
import qualified Data.List                   as L
import           Data.Maybe
import qualified Data.Text.Lazy              as T
import qualified Data.Text.Short             as ST
import qualified Data.Vector                 as V
import qualified Data.Vector.Generic         as GV
import qualified Data.Vector.Storable        as SV
import qualified Data.Vector.Unboxed         as UV
import           Data.Word
import           Debug.Trace
import           Foreign.C.Types             (CUInt)
import           GeneralUtilities
import qualified GraphOptimization.Medians   as M
import qualified Input.BitPack               as BP
import qualified ParallelUtilities           as PU
import qualified SymMatrix                   as S
import           Text.Read
import           Types.Types
import qualified Utilities.Utilities         as U
import qualified Data.MetricRepresentation   as MR


-- | optimizePrealignedData convert
        -- prealigned to non-additive or matrix
            -- here
        -- bitPack new non-additive
            -- packNonAdditive
optimizePrealignedData :: GlobalSettings -> ProcessedData -> ProcessedData
optimizePrealignedData inGS inData@(_, _, blockDataVect) = 
    let -- remove constant characters from prealigned
        inData' = removeConstantCharactersPrealigned inData
        
        -- convert prealigned to nonadditive if all 1 tcms
        inData'' = convertPrealignedToNonAdditive inData'
    
        -- bit packing for non-additivecharacters
        inData''' = BP.packNonAdditiveData inGS inData''

    in
    
    if U.getNumberPrealignedCharacters blockDataVect == 0 then 
        -- trace ("Not Bitpacking...") 
        inData
    else 
        --trace ("Bitpacking...")
        inData'''
    
    

-- | convertPrealignedToNonAdditive converts prealigned data to non-additive
-- if homogeneous TCM (all 1's non-diagnoal)
convertPrealignedToNonAdditive :: ProcessedData -> ProcessedData
convertPrealignedToNonAdditive (nameVect, bvNameVect, blockDataVect) = (nameVect, bvNameVect, fmap convertPrealignedToNonAdditiveBlock blockDataVect)

-- | convertPrealignedToNonAdditiveBlock takes a character block and convertes prealigned to non-add if tcms all 1's
-- this is done taxon by taxon and character by character since can convert with only local infomation
convertPrealignedToNonAdditiveBlock :: BlockData -> BlockData
convertPrealignedToNonAdditiveBlock (nameBlock, charDataVV, charInfoV) =
    let codingTypeV = fmap fst $ fmap getRecodingType (fmap costMatrix charInfoV)
        (newCharDataVV, newCharInfoVV) = V.unzip $ fmap (convertTaxonPrealignedToNonAdd charInfoV codingTypeV) charDataVV
    in
    (nameBlock, newCharDataVV, V.head newCharInfoVV)

-- | convertTaxonPrealignedToNonAdd takes a vector of character info and vector of charcter data for a taxon
-- and recodes prealigned as non-additve if tcm's all 1s
convertTaxonPrealignedToNonAdd :: V.Vector CharInfo -> V.Vector String -> V.Vector CharacterData -> (V.Vector CharacterData, V.Vector CharInfo)
convertTaxonPrealignedToNonAdd charInfoV codingTypeV charDataV =
    V.unzip $ V.zipWith3 convertTaxonPrealignedToNonAddCharacter charInfoV codingTypeV charDataV

-- | convertTaxonPrealignedToNonAddCharacter takes a taxon character and char info and cost matrix type
-- and transforms to non-additive if all tcms are 1's
convertTaxonPrealignedToNonAddCharacter :: CharInfo -> String -> CharacterData -> (CharacterData, CharInfo)
convertTaxonPrealignedToNonAddCharacter charInfo matrixType charData =
    -- trace ("CTP: " ++ (show $ charType charInfo) ++ " " ++ (show $ matrixType)) (
    if charType charInfo `notElem` prealignedCharacterTypes then (charData, charInfo)
    else if matrixType /= "nonAdd" then (charData, charInfo)
    else
            -- this coefficient do to reducion by integer factors in matrix representation
            -- need to be multipled or 2:2 goes to 1:1 and weigh/cost incorrect (non for slim--ffi etc)
        let matrixCoefficient = if (charType charInfo) == AlignedSlim then 1
                                else if (charType charInfo) == AlignedWide then MR.minInDelCost (wideTCM charInfo)
                                else if (charType charInfo) == AlignedHuge then MR.minInDelCost (hugeTCM charInfo)
                                else error ("Unrecognized--unimplemented character type in convertTaxonPrealignedToNonAddCharacter: " ++ (show $ charType charInfo))

            charWeight = weight charInfo

            newStateBV = if charType charInfo == AlignedSlim then
                            convert2BVTriple 32 $ (snd3 . alignedSlimPrelim) charData
                         else if charType charInfo == AlignedWide then
                            convert2BVTriple 64 $ (snd3 . alignedWidePrelim) charData
                         else if charType charInfo == AlignedHuge then
                            (snd3 $ alignedHugePrelim charData, snd3 $ alignedHugePrelim charData, snd3 $ alignedHugePrelim charData)
                         else error ("Unrecognized character type in convertTaxonPrealignedToNonAddCharacter: " ++ (show $ charType charInfo))
        in
        (emptyCharacter {stateBVPrelim = newStateBV}, charInfo {charType = NonAdd, weight = charWeight * (fromIntegral $ fromEnum matrixCoefficient)})
        -- )


-- | convert2BVTriple takes CUInt or Word64 and converts to Triple Vector of bitvectors
convert2BVTriple :: (Integral a, GV.Vector v a) => Word -> v a -> (V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector)
convert2BVTriple size inM =
   let inMList = GV.toList inM
       inMBV = fmap (BV.fromNumber size) inMList

   in
   (V.fromList inMBV, V.fromList inMBV, V.fromList inMBV)

-- | convert2BV takes CUInt or Word64 and converts to Vector of bitvectors
-- this for leaves so assume M only one needed really
convert2BV :: (Integral a, GV.Vector v a) => Word -> (v a, v a, v a) -> (V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector)
convert2BV size (_, inM, _) =
   let inMList = GV.toList inM
       inMBV = fmap (BV.fromNumber size) inMList

   in
   (V.fromList inMBV, V.fromList inMBV, V.fromList inMBV)

-- | getRecodingType takes a cost matrix and detemines if it can be recodes as non-additive,
-- non-additive with gap chars, or matrix
-- assumes indel costs are in last row and column
getRecodingType :: S.Matrix Int -> (String, Int)
getRecodingType inMatrix =
   if S.null inMatrix then error "Null matrix in getRecodingType"
   else
      if (not . S.isSymmetric) inMatrix then ("matrix",  0)
      else
         let matrixLL = S.toFullLists inMatrix
             lastRow = L.last matrixLL
             numUniqueCosts = length $ L.group $ L.sort $ (filter (/= 0) $ concat matrixLL)

         in
         -- trace  ("GRT: " ++ (show numUniqueCosts)) (
         -- all same except for 0
         if numUniqueCosts == 1 then ("nonAdd", 0)

         else ("matrix",  head lastRow)
         {-
         -- all same except for gaps
         else if numUniqueCosts == 2 then
            trace ("NAG: " ++ (show $ length $ L.group $ L.sort $ filter (/= 0) lastRow)) (
            if  (length $ L.group $ filter (/= 0) lastRow) == 1 then ("nonAddGap", head lastRow)

            -- some no gaps different
            else ("matrix",  head lastRow)
            )
         -- to many types for nonadd coding
         else ("matrix",  head lastRow)

         -}
         -- )

-- | reBlockData takes original block assignments--each input file is a block--
-- and combines, creates new, deletes empty blocks from user input
-- reblock pair names may contain wildcards
reBlockData :: [(NameText, NameText)] -> ProcessedData -> ProcessedData
reBlockData reBlockPairs inData@(leafNames, leafBVs, blockDataV) =
    -- trace ("RBD:" ++ (show $ fmap fst3 blockDataV )) (
    if null reBlockPairs then trace "Character Blocks as input files" inData
    else
        let -- those block to be reassigned--nub in case repeated names
            toBeReblockedNames = fmap (T.filter (/= '"')) $ L.nub $ fmap snd reBlockPairs
            unChangedBlocks = V.filter ((`notElemWildcards` toBeReblockedNames).fst3) blockDataV
            blocksToChange = V.filter ((`elemWildcards` toBeReblockedNames).fst3) blockDataV
            newBlocks = makeNewBlocks reBlockPairs blocksToChange []
            reblockedBlocks = unChangedBlocks V.++ V.fromList newBlocks
        in
        trace ("\nReblocking: " ++ show toBeReblockedNames ++ " leaving unchanged: " ++ show (fmap fst3 unChangedBlocks)
            ++ "\n\tNew blocks: " ++ show (fmap fst3 reblockedBlocks) ++ "\n")
        (leafNames, leafBVs, reblockedBlocks)
        -- )

-- | makeNewBlocks takes lists of reblock pairs and existing relevant blocks and creates new blocks returned as a list
makeNewBlocks :: [(NameText, NameText)] -> V.Vector BlockData -> [BlockData] -> [BlockData]
makeNewBlocks reBlockPairs inBlockV curBlockList
  | null reBlockPairs = curBlockList
  | V.null inBlockV && null curBlockList =
    errorWithoutStackTrace ("Reblock pair names do not have a match for any input block--perhaps missing '#0/N'? Blocks: " ++ (show $ fmap snd reBlockPairs))
  | V.null inBlockV = curBlockList
  | otherwise =
    let firstBlock = V.head inBlockV
        firstName = fst3 firstBlock
        newPairList = fst <$> filter (textMatchWildcards firstName.snd) reBlockPairs
    in
    if null newPairList then errorWithoutStackTrace ("Reblock pair names do not have a match for any input block--perhaps missing '#0'? Specified pairs: " ++ show reBlockPairs
        ++ " input block name: " ++ T.unpack firstName)
    else if length newPairList > 1 then errorWithoutStackTrace ("Multiple reblock destinations for single input block" ++ show newPairList)
    else
        let newBlockName = head newPairList
            existingBlock = filter ((==newBlockName).fst3) curBlockList
        in
        -- new block to be created
        if null existingBlock then
            --trace("NBlocks:" ++ (show $ fmap fst3 curBlockList))
            makeNewBlocks reBlockPairs (V.tail inBlockV) ((newBlockName, snd3 firstBlock, thd3 firstBlock) : curBlockList)

        -- existing block to be added to
        else if length existingBlock > 1 then error ("Error: Block to be added to more than one block-should not happen: " ++ show reBlockPairs)
        else
            -- need to add character vectors to vertex vectors and add to CharInfo
            -- could be multiple  'characteres' if non-exact data in inpuit file (Add, NonAdd, MAtrix etc)
            let blockToAddTo = head existingBlock
                newCharData  = V.zipWith (V.++) (snd3 blockToAddTo) (snd3 firstBlock)
                newCharInfo  = thd3 blockToAddTo V.++ thd3 firstBlock
            in
            --trace("EBlocks:" ++ (show $ fmap fst3 curBlockList))
            makeNewBlocks reBlockPairs (V.tail inBlockV) ((newBlockName, newCharData, newCharInfo) : filter ((/=newBlockName).fst3) curBlockList)


-- | combineDataByType combines data of same type for (exact types) into 
-- same vectors in character so have one non-add, one add, one of each packed type, 
-- can have multiple matrix (due to cost matrix differneces)
-- simialr result to groupDataByType, but does not assume single characters.
combineDataByType :: GlobalSettings -> ProcessedData -> ProcessedData
combineDataByType inGS inData@(taxNames, taxBVNames, _) =
    --recode add to non-add before combine-- takes care wor integer weighting 
    let (_, _, blockDataV') = recodeAddToNonAddCharacters inGS maxAddStatesToRecode inData
        recodedData = fmap combineData blockDataV' 
    in
    (taxNames, taxBVNames, recodedData)

-- | combineData creates for a block) lists of each data type and concats then creating new data and new char info
combineData :: BlockData -> BlockData
combineData (blockName, blockDataVV, charInfoV) =
    let (newBlockDataLV, newCharInfoLV) = unzip (PU.seqParMap rdeepseq (combineBlockData charInfoV) (V.toList blockDataVV)) -- `using` PU.myParListChunkRDS)
    in 
    (blockName, V.fromList newBlockDataLV, head newCharInfoLV)

-- | combineBlockData takes a vector of char info and vector or charcater data for a taxon and 
-- combined exact data types into single characters
-- additive characters should have already been converted (if less that maxAddStatesToRecode) to binary
-- non-additive characters with integer weights could be  repeated before combining-- but this has been disabled 
-- due to memory issues in large data sets
-- non-integer additive and non-additive are grouped together by weight so they can be combined and bit packed
-- all other types are grouped by weight for efficiency of optimization and reductionn of multiplies
-- zero weight characters are filtered out
combineBlockData :: V.Vector CharInfo -> V.Vector CharacterData -> (V.Vector CharacterData, V.Vector CharInfo)
combineBlockData inCharInfoV inCharDataV = 
    let pairCharsInfo = V.zip inCharInfoV inCharDataV

        -- characters to not be reorganized-- nbasically the sequence characters
        sequenceCharacters = V.toList $ V.filter ((> 0) . weight . fst) $ V.filter ((`elem` sequenceCharacterTypes) . charType . fst)  pairCharsInfo

        -- matrix characters are more complex--can only join if same matrix
        matrixCharsPair = V.filter ((> 0) . weight . fst) $ V.filter ((== Matrix) . charType . fst) pairCharsInfo
        (newMatrixCharL, newMatrixCharInfoL) = if (not . null) matrixCharsPair then unzip $ organizeMatrixCharsByMatrix (V.toList matrixCharsPair)
                                               else ([],[])

        -- non-additive characters
            -- multiple characters by weight, if only 1 weight then all together
        nonAddChars = V.filter ((> 0) . weight . fst) $ V.filter ((== NonAdd) . charType . fst) pairCharsInfo
        (newNonAddCharInfo, newNonAddChar) = unzip $ V.toList $ groupCharactersByWeight nonAddChars

        -- additive characters
            -- multiple characters by weight, if only 1 weight then all together
            
        addChars = V.filter ((> 0) . weight . fst) $ V.filter ((== Add) . charType . fst) pairCharsInfo
        (newAddCharInfo, newAddChar) = unzip $ V.toList $ groupCharactersByWeight addChars

        -- Packed2 characters
        packed2Chars = V.filter ((> 0) . weight . fst) $ V.filter ((== Packed2) . charType . fst) pairCharsInfo
        (newPacked2CharInfo, newPacked2Char) = unzip $ V.toList $ groupCharactersByWeight packed2Chars

        -- Packed4 characters
        packed4Chars = V.filter ((> 0) . weight . fst) $ V.filter ((== Packed4) . charType . fst) pairCharsInfo
        (newPacked4CharInfo, newPacked4Char) = unzip $ V.toList $ groupCharactersByWeight packed4Chars


        -- Packed5 characters
        packed5Chars = V.filter ((> 0) . weight . fst) $ V.filter ((== Packed5) . charType . fst) pairCharsInfo
        (newPacked5CharInfo, newPacked5Char) = unzip $ V.toList $ groupCharactersByWeight packed5Chars


        -- Packed8 characters
        packed8Chars = V.filter ((> 0) . weight . fst) $ V.filter ((== Packed8) . charType . fst) pairCharsInfo
        (newPacked8CharInfo, newPacked8Char) = unzip $ V.toList $ groupCharactersByWeight packed8Chars

        -- Packed64 characters
        packed64Chars = V.filter ((> 0) . weight . fst) $ V.filter ((== Packed64) . charType . fst) pairCharsInfo
        (newPacked64CharInfo, newPacked64Char) = unzip $ V.toList $ groupCharactersByWeight packed64Chars

        -- Add together all new characters, seqeunce characters and char info
        -- newCharList = newNonAddCharL ++ (V.toList nonAddCharsWeightNotInt) ++ newAddCharL ++ (V.toList addCharsWeightNot1) ++ newPacked2CharL ++ newPacked4CharL ++ newPacked5CharL ++ newPacked8CharL ++ newPacked64CharL ++ newMatrixCharL ++ (fmap snd sequenceCharacters)
        -- newCharInfoList = newNonAddCharInfoL ++ (V.toList nonAddCharsWeightNotIntInfo) ++ newAddCharInfoL ++ (V.toList addCharsWeightNot1Info) ++ newPacked2CharInfoL ++ newPacked4CharInfoL ++ newPacked5CharInfoL ++ newPacked8CharInfoL ++ newPacked64CharInfoL ++ newMatrixCharInfoL ++ (fmap fst sequenceCharacters)
        
        newCharList = newNonAddChar ++ newAddChar ++ newPacked2Char ++ newPacked4Char ++ newPacked5Char ++ newPacked8Char ++ newPacked64Char ++ newMatrixCharL ++ (fmap snd sequenceCharacters)
        newCharInfoList = newNonAddCharInfo ++ newAddCharInfo ++ newPacked2CharInfo ++ newPacked4CharInfo ++ newPacked5CharInfo ++ newPacked8CharInfo ++ newPacked64CharInfo ++ newMatrixCharInfoL ++ (fmap fst sequenceCharacters)

    in
    (V.fromList newCharList, V.fromList newCharInfoList)

-- | groupCharactersByWeight takes a list of characters and returns a list of lists of charcter with same weight
-- checked as Double.  
-- NB--Does not check for character type---assuming this is desired it must be assured before input
groupCharactersByWeight :: V.Vector (CharInfo, CharacterData) -> V.Vector (CharInfo, CharacterData)
groupCharactersByWeight inCharsPairList =
    if V.null inCharsPairList then V.empty
    else 
        let weightList = L.nub $ V.toList $ fmap weight (fmap fst inCharsPairList)
            charListByWeight = fmap (getSameWeightChars inCharsPairList) $ V.fromList weightList
        in
        V.concatMap mergeCharacters charListByWeight

-- | mergeCharacters merges the data fieed of characters based on type
-- NB--Does not check all chars are same type or weight--will silently combine mempty values
-- with whatever weight.  
-- returns list with single member so can concat later
mergeCharacters :: V.Vector (CharInfo, CharacterData) -> V.Vector (CharInfo, CharacterData)
mergeCharacters inCharsPairList =
    if V.null inCharsPairList then V.empty
    else 
        let thisCharType = (charType . fst . V.head) inCharsPairList
            
            -- non-add data
            dataFieldNonAdd = V.concatMap (snd3 . stateBVPrelim . snd) inCharsPairList
            newNonAddChar = ((snd . V.head) inCharsPairList) {stateBVPrelim = (dataFieldNonAdd, dataFieldNonAdd, dataFieldNonAdd), stateBVFinal = dataFieldNonAdd}
            
            -- add data
            dataFieldAdd = V.concatMap (snd3 . rangePrelim . snd) inCharsPairList
            newAddChar = ((snd . V.head) inCharsPairList) {rangePrelim = (dataFieldAdd, dataFieldAdd, dataFieldAdd), rangeFinal = dataFieldAdd}
            
            -- packed data
            dataFieldPacked = UV.concat $ V.toList $ fmap (snd3 . packedNonAddPrelim . snd) inCharsPairList
            newPackedChar = ((snd . V.head) inCharsPairList) {packedNonAddPrelim = (dataFieldPacked, dataFieldPacked, dataFieldPacked), packedNonAddFinal = dataFieldPacked}

            -- add info of merging to character info
            newOrigInfo = V.concatMap (origInfo . fst) inCharsPairList
            newCharInfo = ((fst . V.head) inCharsPairList) {origInfo = newOrigInfo}

        in
        if thisCharType == NonAdd then V.singleton (newCharInfo, newNonAddChar)
        else if thisCharType == Add then V.singleton (newCharInfo, newAddChar)
        else if thisCharType `elem` packedNonAddTypes then V.singleton (newCharInfo, newPackedChar)
        else error ("Error in mergeCharacters: Character type " ++ show thisCharType ++ " unrecognized/not implemented")



-- | getSameWeightChars returns character pairs with same matrix as testMatrix
getSameWeightChars ::  V.Vector (CharInfo, CharacterData) -> Double -> V.Vector (CharInfo, CharacterData)
getSameWeightChars inCharsPairList testWeight =
    if V.null inCharsPairList then V.empty
    else 
        let inWeightList = fmap weight (fmap fst inCharsPairList)
            weightPairPair = V.zip inWeightList inCharsPairList
            matchList = V.filter ((== testWeight) . weight . fst . snd) weightPairPair
        in
        fmap snd matchList

{- Currently unused--employed to reduce chrater umbers by replicating characters with integer weight
    can cause memory issues with large data sets
-- | replicateCharPairByWeight replicates characters by integer weight
replicateCharPairByWeight :: (CharInfo, CharacterData) -> [(CharInfo, CharacterData)] 
replicateCharPairByWeight firstPair = 
    let charIntWeight = doubleAsInt $ (weight . fst) firstPair
    in
    if isNothing charIntWeight then error ("Character weight not an integer in replicateCharPair: " ++ (show $ (weight . fst) firstPair))
    else 
        (replicate (fromJust charIntWeight) firstPair)
-}

-- | organizeMatrixCharsByMatrix combines matrix charcters if they have the same cost matrix
organizeMatrixCharsByMatrix :: [(CharInfo, CharacterData)] -> [(CharacterData, CharInfo)]
organizeMatrixCharsByMatrix inCharsPairList = 
    if null inCharsPairList then []
    else 
        let costMatrixList = L.nub $ fmap costMatrix (fmap fst inCharsPairList)
            charMatrixLL = fmap (getSameMatrixChars inCharsPairList) costMatrixList
            newMatrixPairs = fmap combineMatrixCharsByMatrix charMatrixLL
        in
        newMatrixPairs 

-- | combineMatrixCharsByMatrix combines matrix characters--assumes cost matrices are the same
combineMatrixCharsByMatrix :: [(CharInfo, CharacterData)] -> (CharacterData, CharInfo)
combineMatrixCharsByMatrix inCharList =
    let newMatrixcharData = V.concat $ fmap matrixStatesPrelim $ fmap snd inCharList
        newMatrixChar = ((snd . head) inCharList) { matrixStatesPrelim = newMatrixcharData
                                                  , matrixStatesFinal = newMatrixcharData
                                                  }
        newMatrixCharInfo = ((fst . head) inCharList) {origInfo = V.concat $ fmap (origInfo . fst) inCharList}
    in
    (newMatrixChar, newMatrixCharInfo)

-- | getSameMatrixChars returns character pairs with same matrix as testMatrix
getSameMatrixChars ::  [(CharInfo, CharacterData)] -> S.Matrix Int -> [(CharInfo, CharacterData)] 
getSameMatrixChars inCharsPairList testMatrix =
    let inMatrixList = fmap costMatrix (fmap fst inCharsPairList)
        matrixPairPair = zip inMatrixList inCharsPairList
        matchList = filter ((== testMatrix) . costMatrix . fst . snd) matrixPairPair
    in
    fmap snd matchList

-- | removeConstantCharactersPrealigned takes processed data and removes constant characters
-- from prealignedCharacterTypes
removeConstantCharactersPrealigned :: ProcessedData -> ProcessedData
removeConstantCharactersPrealigned (nameVect, bvNameVect, blockDataVect) =
    let newBlockData = V.fromList (PU.seqParMap rdeepseq  removeConstantBlockPrealigned (V.toList blockDataVect)) -- `using` PU.myParListChunkRDS)
    in
    (nameVect, bvNameVect, newBlockData)

-- | removeConstantBlockPrealigned takes block data and removes constant characters
removeConstantBlockPrealigned :: BlockData -> BlockData
removeConstantBlockPrealigned inBlockData@(blockName, taxVectByCharVect, charInfoV) =
    -- check for null data--really reallyu shouldn't happen
    if V.null taxVectByCharVect then trace ("Warning: Null block data in removeConstantBlockPrealigned") inBlockData
    
    -- check for prealigned data in block
    else if U.getNumberPrealignedCharacters (V.singleton inBlockData) == 0 then inBlockData
    else 
        let numChars = V.length $ V.head taxVectByCharVect

            -- create vector of single characters with vector of taxon data of sngle character each
            -- like a standard matrix with a single character
            singleCharVect = fmap (U.getSingleCharacter taxVectByCharVect) (V.fromList [0.. numChars - 1])

            -- actually remove constants form chaarcter list
            singleCharVect' = V.zipWith removeConstantCharsPrealigned singleCharVect charInfoV

            -- recreate the taxa vext by character vect block data expects
            -- should filter out length zero characters
            newTaxVectByCharVect = U.glueBackTaxChar singleCharVect'
         in
         (blockName, newTaxVectByCharVect, charInfoV)

-- | removeConstantCharsPrealigned takes a single 'character' and if proper type removes if all values are the same
-- could be done if character has max lenght of 0 as well.
-- packed types already filtered when created
removeConstantCharsPrealigned :: V.Vector CharacterData -> CharInfo -> V.Vector CharacterData
removeConstantCharsPrealigned singleChar charInfo =
    let inCharType = charType charInfo
    in

    -- dynamic characters don't do this
    if inCharType `notElem` prealignedCharacterTypes then singleChar
    else
        let variableVect = getVariableChars inCharType singleChar
        in
        variableVect

-- | getVariableChars checks identity of states in a vector positin in all taxa
-- and returns True if variable, False if constant
-- bit packed and non-exact should not get in here
getVariableChars :: CharType -> V.Vector CharacterData -> V.Vector CharacterData
getVariableChars inCharType singleChar =
    let nonAddV = fmap snd3 $ fmap stateBVPrelim singleChar
        addV    = fmap snd3 $ fmap rangePrelim singleChar
        matrixV = fmap matrixStatesPrelim singleChar
        alSlimV = fmap snd3 $ fmap alignedSlimPrelim singleChar
        alWideV = fmap snd3 $ fmap alignedWidePrelim singleChar
        alHugeV = fmap snd3 $ fmap alignedHugePrelim singleChar

        -- get identity vect
        boolVar =   if inCharType == NonAdd then getVarVectBits inCharType nonAddV []
                    else if inCharType == Add then getVarVectAdd addV []
                    else if inCharType == Matrix then getVarVectMatrix matrixV []
                    else if inCharType == AlignedSlim then getVarVectBits inCharType alSlimV []
                    else if inCharType == AlignedWide then getVarVectBits inCharType alWideV []
                    else if inCharType == AlignedHuge then getVarVectBits inCharType alHugeV []
                    else error ("Char type unrecognized in getVariableChars: " ++ show inCharType)

        -- get Variable characters by type
        nonAddVariable = fmap (filterConstantsV (V.fromList boolVar)) nonAddV
        addVariable    = fmap (filterConstantsV (V.fromList boolVar)) addV
        matrixVariable = fmap (filterConstantsV (V.fromList boolVar)) matrixV
        alSlimVariable = fmap (filterConstantsSV (V.fromList boolVar)) alSlimV
        alWideVariable = fmap (filterConstantsUV (V.fromList boolVar)) alWideV
        alHugeVariable = fmap (filterConstantsV (V.fromList boolVar)) alHugeV

        -- assign to propoer character fields
        outCharVect = V.zipWith (assignNewField inCharType) singleChar (V.zip6 nonAddVariable addVariable matrixVariable alSlimVariable alWideVariable alHugeVariable)

    in
    -- trace ("GVC:" ++ (show $ length boolVar) ++ " -> " ++ (show $ length $ filter (== False) boolVar))
    outCharVect

-- | getVarVectAdd takes a vector of a vector additive ranges and returns False if range overlap
-- True if not (short circuits)
-- based on range overlap
getVarVectAdd :: V.Vector (V.Vector (Int, Int)) -> [Bool] -> [Bool]
getVarVectAdd stateVV curBoolList =
    if V.null (V.head stateVV) then L.reverse curBoolList

    else
        let firstChar = fmap V.head stateVV
            isVariable = checkIsVariableAdditive (V.head firstChar) (V.tail firstChar)

        in
        getVarVectAdd (fmap V.tail stateVV) (isVariable : curBoolList)


-- | getVarVectMatrix takes a generic vector and returns False if values are same
-- True if not (short circuits)
-- based on simple identity not max cost zero
getVarVectMatrix :: V.Vector (V.Vector (V.Vector MatrixTriple)) -> [Bool] -> [Bool]
getVarVectMatrix stateVV curBoolList =
    if V.null (V.head stateVV) then L.reverse curBoolList

    else
        let firstChar = fmap V.head stateVV
            isVariable = checkIsVariableMatrix (getMatrixStateList $ V.head firstChar) (V.tail firstChar)

        in
        getVarVectMatrix (fmap V.tail stateVV) (isVariable : curBoolList)


-- | getVarVectBits takes a generic vector and returns False if values are same
-- True if not (short circuits)
-- based on simple identity not max cost zero
getVarVectBits :: (FiniteBits a, Eq a, GV.Vector v a) => CharType -> V.Vector (v a) -> [Bool] -> [Bool]
getVarVectBits inCharType stateVV curBoolList =
    if GV.null (V.head stateVV) then L.reverse curBoolList

    else
        let firstChar = fmap GV.head stateVV
            isVariable = if inCharType      == NonAdd       then checkIsVariableBit (GV.head firstChar) (GV.tail firstChar)
                         else if inCharType == AlignedSlim  then checkIsVariableBit (GV.head firstChar) (GV.tail firstChar)
                         else if inCharType == AlignedWide  then checkIsVariableBit (GV.head firstChar) (GV.tail firstChar)
                         else if inCharType == AlignedHuge  then checkIsVariableBit (GV.head firstChar) (GV.tail firstChar)
                         else error ("Char type unrecognized in getVariableChars: " ++ show inCharType)

        in
        getVarVectBits inCharType (fmap GV.tail stateVV) (isVariable : curBoolList)

{-
-- | checkIsVariableIdentity takes a generic vector and sees if all elements are identical
checkIsVariableIdentity ::  (Eq a, GV.Vector v a) => a -> v a -> Bool
checkIsVariableIdentity firstElement inVect =
    if GV.null inVect then False
    else
        if firstElement /= GV.head inVect then True
        else checkIsVariableIdentity firstElement (GV.tail inVect)
-}

-- | checkIsVariableAdditive checks if additive charcter is variable
-- by taking new ranges and range costs of first element with all others
-- if summed cost > 0 then variable
checkIsVariableAdditive :: (Int, Int) -> V.Vector (Int, Int) -> Bool
checkIsVariableAdditive (ir1, ir2) rangeList =
    if V.null rangeList then False
    else
        let (nr1, nr2) = V.head rangeList
            (newMin, newMax, newCost) = M.getNewRange ir1 ir2 nr1 nr2
        in
        if newCost > 0 then True
        else checkIsVariableAdditive (newMin, newMax) (V.tail rangeList)

-- | getMatrixStateList returns minimum matrix characters states as integers
getMatrixStateList :: V.Vector MatrixTriple -> [Int]
getMatrixStateList inState =
    let statePairList = zip (fmap fst3 $ V.toList inState) [0..(V.length inState - 1)]
        minCost = minimum $ fmap fst3 $ V.toList inState
        stateList = fmap snd $ filter ((== minCost) .fst) statePairList
    in
    stateList

-- | checkIsVariableMatrix checks if matrix charcter is variable
-- by taking min cost states of first element and
-- checks for overlap--if empty list intersect then variable
checkIsVariableMatrix :: [Int] -> V.Vector (V.Vector MatrixTriple) -> Bool
checkIsVariableMatrix inStateList restStatesV =
    if V.null restStatesV then False
    else
        let nextStateList = getMatrixStateList $ V.head restStatesV

            newStateList = L.intersect inStateList nextStateList
        in
        if null newStateList then True
        else checkIsVariableMatrix newStateList (V.tail restStatesV)

-- | checkIsVariableBit takes a generic vector and checks for
-- state overlap via bit AND (.&.)
checkIsVariableBit ::  (FiniteBits a, GV.Vector v a) => a -> v a -> Bool
checkIsVariableBit firstElement restVect =
    if GV.null restVect then False
    else
        let newState = firstElement .&. (GV.head restVect)
        in
        if popCount newState == 0 then True
        else checkIsVariableBit newState (GV.tail restVect)

-- | these need to be abstracted but had problems with the bool list -> Generic vector, and SV pair

-- | filerConstantsV takes the charcter data and filters out teh constants
-- uses filter to keep O(n)
--filterConstantsV :: (GV.Vector v a) => [Bool] -> v a -> v a
filterConstantsV :: V.Vector Bool -> V.Vector a -> V.Vector a
filterConstantsV inVarBoolV charVect =
    let pairVect = V.zip charVect inVarBoolV
        variableCharV = V.map fst $ V.filter ((== True) . snd) pairVect
    in
    variableCharV


-- | filerConstantsSV takes the charcter data and filters out teh constants
-- uses filter to keep O(n)
--filterConstantsV :: (GV.Vector v a) => [Bool] -> v a -> v a
filterConstantsSV ::  (SV.Storable a) => V.Vector Bool -> SV.Vector a -> SV.Vector a
filterConstantsSV inVarBoolV charVect =
    let varVect = filterConstantsV inVarBoolV (V.fromList $ SV.toList charVect)
    in
    SV.fromList $ V.toList varVect

-- | filerConstantsUV takes the charcter data and filters out teh constants
-- uses filter to keep O(n)
--filterConstantsV :: (GV.Vector v a) => [Bool] -> v a -> v a
filterConstantsUV ::  (UV.Unbox a) => V.Vector Bool -> UV.Vector a -> UV.Vector a
filterConstantsUV inVarBoolV charVect =
    let varVect = filterConstantsV inVarBoolV (V.fromList $ UV.toList charVect)
    in
    UV.fromList $ V.toList varVect

-- | assignNewField takes character type and a 6-tuple of charcter fields and assigns the appropriate
-- to the correct field
-- neither bit packed nor nno-exact should het here
assignNewField :: CharType
               -> CharacterData
               -> (V.Vector BV.BitVector, V.Vector (Int, Int), V.Vector (V.Vector MatrixTriple), SV.Vector CUInt, UV.Vector Word64, V.Vector BV.BitVector)
               -> CharacterData
assignNewField inCharType charData (nonAddData, addData, matrixData, alignedSlimData, alignedWideData, alignedHugeData) =
    if inCharType == NonAdd then charData {stateBVPrelim = (nonAddData, nonAddData, nonAddData)}
    else if inCharType == Add then charData {rangePrelim = (addData, addData, addData)}
    else if inCharType == Matrix then charData {matrixStatesPrelim = matrixData}
    else if inCharType == AlignedSlim then charData {alignedSlimPrelim = (alignedSlimData, alignedSlimData, alignedSlimData)}
    else if inCharType == AlignedWide then charData {alignedWidePrelim = (alignedWideData, alignedWideData, alignedWideData)}
    else if inCharType == AlignedHuge then charData {alignedHugePrelim = (alignedHugeData, alignedHugeData, alignedHugeData)}
    else error ("Char type unrecognized in assignNewField: " ++ show inCharType)

-- | recodeAddToNonAddCharacters takes an max states number and processsed data
-- and recodes additive characters with max state < input max (0..input max - 1)
-- as a series of binary non-additive characters
recodeAddToNonAddCharacters :: GlobalSettings -> Int -> ProcessedData -> ProcessedData
recodeAddToNonAddCharacters inGS maxStateToRecode (nameVect, nameBVVect, blockDataVect) =
    let newBlockDataVect = fmap (convertAddToNonAddBlock inGS maxStateToRecode) blockDataVect
    in
    (nameVect, nameBVVect, newBlockDataVect)

-- | convertAddToNonAddBlock converts additive charcters to no-additive in a block
convertAddToNonAddBlock :: GlobalSettings ->  Int -> BlockData -> BlockData
convertAddToNonAddBlock inGS maxStateToRecode (blockName, taxByCharDataVV, charInfoV) =
    let (newTaxByCharDataVV, newCharInfoVV) = V.unzip $ fmap (recodeTaxonData inGS maxStateToRecode charInfoV) taxByCharDataVV
    in
    -- trace ("CNAB: " ++ (show (V.length $ V.head newTaxByCharDataVV, V.length $ V.head newCharInfoVV)))
    (blockName, newTaxByCharDataVV, V.head newCharInfoVV)

-- | recodeTaxonData recodes Add as nonAdd for each taxon in turn
recodeTaxonData :: GlobalSettings -> Int -> V.Vector CharInfo -> V.Vector CharacterData -> (V.Vector CharacterData, V.Vector CharInfo)
recodeTaxonData inGS maxStateToRecode charInfoV taxonCharacterDataV =
    let (newCharDataVV, newCharInfoVV) = unzip $ zipWith (recodeAddToNonAddCharacter inGS maxStateToRecode) (V.toList taxonCharacterDataV) (V.toList charInfoV)
    in
    -- trace ("RTD: " ++ (show (V.length $ V.concat newCharDataVV, V.length $ V.concat newCharInfoVV)))
    (V.concat newCharDataVV, V.concat newCharInfoVV)

-- |recodeAddToNonAddCharacter takes a single character for single taxon and recodes if non-additive with
-- fewer than maxStateToRecode states.
-- assumes states in linear order
-- replicatee charinfo for multiple new characters after recoding
recodeAddToNonAddCharacter :: GlobalSettings -> Int -> CharacterData -> CharInfo -> (V.Vector CharacterData,  V.Vector CharInfo)
recodeAddToNonAddCharacter inGS maxStateToRecode inCharData inCharInfo =
    let inCharType = charType inCharInfo
        numStates = 1 + (L.maximum $ fmap makeInt $ alphabet inCharInfo) -- min 2 (1 + (L.last $ L.sort $ fmap makeInt $ alphabetSymbols $ alphabet inCharInfo))
        -- numStates = 1 + (L.last $ L.sort $ fmap makeInt $ alphabetSymbols $ alphabet inCharInfo)
        origName = name inCharInfo
    in
    -- if a single state recodes to a single uninfomative binary 
        -- removed || ((not . doubleIsInt . weight) inCharInfo)  to allow for recodding (leaving weight) for non-integer weights
    if (inCharType /= Add) then (V.singleton inCharData, V.singleton inCharInfo)
    
        -- the limit on recoded states is removed for PMDL/ML since otherwise bit costs will be incorrect
    else if (numStates > maxStateToRecode) && ((optimalityCriterion inGS)  `notElem` [PMDL, Likelihood]) then (V.singleton inCharData, V.singleton inCharInfo)
    else if numStates < 2 then (V.empty, V.empty)
    else 
        -- create numStates - 1 no-additve chaaracters (V.singleton inCharData, V.singleton inCharInfo)
        -- bits ON-- [0.. snd range]
        let minRangeIndex = fst $ V.head $ snd3 $ rangePrelim inCharData
            maxRangeIndex = snd $ V.head $ snd3 $ rangePrelim inCharData
            inCharOrigData = origInfo inCharInfo
            newCharInfo = inCharInfo { name = (T.pack $ (T.unpack origName) ++ "RecodedToNonAdd")
                                     , charType = NonAdd
                                     , alphabet = fromSymbols $ fmap ST.fromString $ fmap show [(0 :: Int), (1 :: Int)]
                                     , origInfo = inCharOrigData
                                     }

            -- create new characters and new character info 
            newCharList = fmap (makeNewNonAddChar minRangeIndex maxRangeIndex) [0..numStates - 2]

            newCharInfoList =  replicate (numStates - 1) newCharInfo                  
            
                                              
        in
        -- trace ("RTNA: Numstates " ++ (show numStates) ++ " " ++ (show $ (snd3 . rangePrelim) inCharData) ++ " -> " ++ (show $ fmap (snd3 . stateBVPrelim) newCharList))
            -- (show (length newCharList, V.length $ V.replicate (numStates - 1) newCharInfo)) ++ "\n" ++ (show newCharList) ++ "\n" ++ (show $ charType newCharInfo))
        (V.fromList newCharList, V.fromList newCharInfoList)
        where makeInt a = let newA = readMaybe (ST.toString a) :: Maybe Int
                          in
                          if isNothing newA then error ("State '" ++ (ST.toString  a) ++ "' not recoding to Int")
                          else fromJust newA

-- | makeNewNonAddCharacter takes a stateIndex and charcatear number
-- and makes a non-additive character with 0 or 1 coding
-- based on stateIndex versus state number
-- if stateIndex > charNumber then 1 else 0 (coded as bit 0 for 0, bit 1 for 1)
makeNewNonAddChar :: Int -> Int -> Int -> CharacterData
makeNewNonAddChar minStateIndex maxStateIndex charIndex =
    let bvMinState = if minStateIndex <= charIndex then BV.fromBits [True, False]
                     else BV.fromBits [False, True]

        bvMaxState = if maxStateIndex <= charIndex then BV.fromBits [True, False]
                     else BV.fromBits [False, True]

        bvState = bvMinState .|. bvMaxState
    in
    emptyCharacter { stateBVPrelim = (V.singleton bvState, V.singleton bvState, V.singleton bvState)
                   , stateBVFinal = V.singleton bvState
                   }

{-  There is an error in this replaced by combineData
import qualified Data.Bifunctor as BF

-- | groupDataByType takes naive data (ProcessedData) and returns PrcessedData
-- with characters reorganized (within blocks)
    -- all non-additive (with same weight) merged to a single vector character
    -- all additive with same alphabet (ie numberical) recoded to single vector
    -- all matrix characters with same costmatrix recoded to single character
    -- removes innactive characters
groupDataByType :: ProcessedData -> ProcessedData
groupDataByType inData =
    let -- before reorganizing--convert additive with <= maxAddStatesToRecode in TYpes.hs states to non-additive ala Farris (1970)
        (nameVect, nameBVVect, blockDataVect) = recodeAddToNonAddCharacters maxAddStatesToRecode inData

        -- reorganize data into same type
        organizedBlockData =  fmap organizeBlockData' (V.toList blockDataVect) `using` PU.myParListChunkRDS
    in
    --trace ("Before Taxa:" ++ (show $ length nameBVVect) ++ " Blocks:" ++ (show $ length blockDataVect) ++ " Characters:" ++ (show $ fmap length $ fmap thd3 blockDataVect)
    --    ++ "\nAfter Taxa:" ++ (show $ length nameBVVect) ++ " Blocks:" ++ (show $ length organizedBlockData) ++ " Characters:" ++ (show $ fmap length $ fmap thd3 organizedBlockData))
    (nameVect, nameBVVect, V.fromList organizedBlockData)

-- | organizeBlockData' special cases and divides characters so that exact characters
-- are reorgnized into single characters by type and cost matrix, while non-exact sequence
-- characters are unchanged.  Characters are reorganized with exact first in block then non-exact
organizeBlockData' :: BlockData -> BlockData
organizeBlockData' localBlockData =
    let numExactChars = U.getNumberExactCharacters (V.singleton localBlockData)
        numNonExactChars = U.getNumberSequenceCharacters (V.singleton localBlockData)
    in
    -- if no exact--nothing to combine
    if numExactChars == 0 then localBlockData

    -- if only non-exact--split and recombine
    else if numNonExactChars == 0 then organizeBlockData [] [] [] [] localBlockData

    -- if both nonexact and exact--pull out non-exact and recombine exact
    else if (numExactChars > 0) && (numNonExactChars > 0) then
        let (exactCharacters, nonSequenceCharacters) = U.splitBlockCharacters (snd3 localBlockData) (thd3 localBlockData) 0 [] []
            newExactCharacters = organizeBlockData [] [] [] [] exactCharacters
            newCharData = V.zipWith (V.++) (snd3 newExactCharacters) (snd3 nonSequenceCharacters)
            newCharInfo = thd3 newExactCharacters V.++ thd3 nonSequenceCharacters
        in
        (fst3 localBlockData, newCharData, newCharInfo)
    else error "This shouldn't happen in organizeBlockData'"

-- | organizeBlockData takes a BlockData element and organizes its character by character type
-- to single add, non-add, matrix, non-exact characters (and those with non-integer weights) are left as is due to their need for
-- individual traversal graphs
-- second element of tuple is a vector over taxa (leaves on input) with
-- a character vector for each leaf/taxon-- basically a matrix with taxon rows and character columns
-- the character info vector is same size as one for each leaf
-- the first 4 args are accumulators for the character types.  Matrix type is list of list since can have multiple
-- matrices.  All non-Exact are in same pile.
-- characters with weight > 1 are recoded as multiples of same character, if weight non-integer geoes into the "unchanged" pile
-- when bit packed later (if non-additive) will have only log 64 operations impact
-- the pairs for some data types are to keep track of things that vary--like matrices and non-exact character information
-- there should be no bit-packed data created yet
organizeBlockData :: [([CharacterData], CharInfo)]
                  -> [([CharacterData], CharInfo)]
                  -> [([[CharacterData]], CharInfo)]
                  -> [([CharacterData], CharInfo)]
                  -> BlockData
                  -> BlockData
organizeBlockData nonAddCharList addCharList matrixCharListList unchangedCharList (blockName, characterDataVectVect, charInfoVect) =
    -- Bit of a cop out but not managing the missing data in blocks thing for multiple non-exact in block
    -- need to add multiple non-exact in block
    if null charInfoVect then
        -- concatenate all new characters, reverse (for good measure), and convert to vectors
        -- with unrecoded non-Exact characters and new CharInfo vector (reversed)
        -- need to make sure the character info is in the order of return types--nonAdd, Add, Matrix etc
        {-Will need a function to add all this stuff back together
        (blockNamne, newCharacterVector, newCharInfoVect)
        -}
        --trace ("New Char lengths :" ++ (show (length nonAddCharList, length addCharList, length matrixCharListList, length unchangedCharList))) (
        let (newCharDataVectVect, newCharInfoVect) = makeNewCharacterData nonAddCharList addCharList matrixCharListList unchangedCharList
        in
        (blockName, newCharDataVectVect, newCharInfoVect)
        -- )
    else
        -- proceed character by character increasing accumulators and consuming character data vector and character infoVect
        -- maybe only accumulate for matrix and non additives?
        let firstCharacter = V.head charInfoVect
            fCharType = charType firstCharacter
            fCharWeight = weight firstCharacter
            intWeight = doubleAsInt fCharWeight
            fCharMatrix = costMatrix firstCharacter
            fCharActivity = activity firstCharacter
            fAlphabet = alphabet firstCharacter
            firstCharacterTaxa = fmap U.safeVectorHead characterDataVectVect
        in
        -- trace ("FCT: " ++ (show $ length firstCharacterTaxa) ++ " " ++ (show characterDataVectVect)) (
        --trace ("CVDD: " ++ (show (length characterDataVectVect, fmap length characterDataVectVect))) (

        -- remove inactive characters
        if not fCharActivity || (length fAlphabet < 2) || fCharWeight == (0 :: Double) then
            --trace ("Innactive/Weight 0")
            organizeBlockData nonAddCharList addCharList matrixCharListList unchangedCharList (blockName, V.map V.tail characterDataVectVect, V.tail charInfoVect)
        else 
            (-- trace ("OBD: " ++ (show intWeight)) (
            if isNothing intWeight then
               -- add to unchanged pile
               let currentUnchangedCharacter = (V.toList firstCharacterTaxa, firstCharacter)

               in
               -- trace ("Unchanged character:" ++ (show $ length $ fst currentUnchangedCharacter) ++ " Name:" ++ (T.unpack $ name firstCharacter) ++ " " ++ (show (charType firstCharacter))
               --    ++ " " ++ (show $ fst currentUnchangedCharacter))
               -- trace ("Character Weight non-integer:" ++ show fCharWeight)
               organizeBlockData nonAddCharList addCharList matrixCharListList (currentUnchangedCharacter : unchangedCharList)  (blockName, V.map V.tail characterDataVectVect, V.tail charInfoVect)

           -- issue with the line "firstCharacterTaxa = fmap V.head characterDataVectVect" since missing character will be empoty and throw an error on V.head
           else if fCharType `notElem` exactCharacterTypes
               then errorWithoutStackTrace "Blocks with more than one Non-Exact Character not yet implemented"

           -- non-additive characters
           else if fCharType == NonAdd then
               let replicateNumber = fromJust intWeight
                   currentNonAdditiveCharacter = (V.toList $ fmap V.head characterDataVectVect, firstCharacter)
               in
               -- trace ("Non-Additive") (
               if replicateNumber == 1 then organizeBlockData (currentNonAdditiveCharacter : nonAddCharList) addCharList matrixCharListList unchangedCharList  (blockName, V.map V.tail characterDataVectVect, V.tail charInfoVect)
               else organizeBlockData (replicate replicateNumber currentNonAdditiveCharacter ++ nonAddCharList) addCharList matrixCharListList unchangedCharList  (blockName, V.map V.tail characterDataVectVect, V.tail charInfoVect)
               -- )

           -- additive characters
           else if fCharType == Add then
               let replicateNumber = fromJust intWeight
                   currentAdditiveCharacter = (V.toList $ fmap V.head characterDataVectVect, firstCharacter)
               in
               -- trace ("Additive") (
               if replicateNumber == 1 then organizeBlockData nonAddCharList (currentAdditiveCharacter : addCharList) matrixCharListList unchangedCharList  (blockName, V.map V.tail characterDataVectVect, V.tail charInfoVect)
               else organizeBlockData nonAddCharList (replicate replicateNumber currentAdditiveCharacter ++ addCharList) matrixCharListList unchangedCharList  (blockName, V.map V.tail characterDataVectVect, V.tail charInfoVect)
               -- )

           -- matrix characters--more complex since need to check for matrix identity
           else if fCharType == Matrix then
               let replicateNumber = fromJust intWeight
                   currentMatrixCharacter = (V.toList $ fmap V.head characterDataVectVect, firstCharacter)
                   newMatrixCharacterList = addMatrixCharacter matrixCharListList fCharMatrix currentMatrixCharacter replicateNumber
               in
               -- trace ("Matrix") (
               organizeBlockData nonAddCharList addCharList newMatrixCharacterList unchangedCharList (blockName, V.map V.tail characterDataVectVect, V.tail charInfoVect)
               -- )

           -- error in type--there should be no bit-packed data created yet.
           else error ("Unrecognized/not implemented charcter type: " ++ show fCharType)) -- )
        


-- | makeNewCharacterData takes nonAddCharList addCharList matrixCharListList unchangedCharList and synthesises them into new character data
-- with a single character for the exact types (nonAdd, Add, Matrix) and mulitple characters for the "unchanged" which includes
-- non-exact characters and those with non-integer weights
-- and character Information vector
-- these only update preliminary of their type--meant to happen before decoration processes
-- emptyCharacter defined in Types
makeNewCharacterData :: [([CharacterData], CharInfo)]
                     -> [([CharacterData], CharInfo)]
                     -> [([[CharacterData]], CharInfo)]
                     -> [([CharacterData], CharInfo)]
                     -> (V.Vector (V.Vector CharacterData), V.Vector CharInfo)
makeNewCharacterData nonAddCharList addCharList matrixCharListList  unchangedCharList =
    let
        -- Non-Additive Characters
        nonAddCharacter = combineNonAdditveCharacters nonAddCharList emptyCharacter []

        -- keep track or original data in origInfo field
        origNonAddData = V.fromList $ zip3 (fmap name $ fmap snd nonAddCharList) (fmap charType $ fmap snd nonAddCharList) (fmap alphabet $ fmap snd nonAddCharList)

        nonAddCharInfo = V.singleton $ (snd $ head nonAddCharList) {name = T.pack "CombinedNonAdditiveCharacters", origInfo = origNonAddData}

        -- Additive Characters
        addCharacter = combineAdditveCharacters addCharList emptyCharacter []

        -- keep track or original data in origInfo field
        origAddData = V.fromList $ zip3 (fmap name $ fmap snd addCharList) (fmap charType $ fmap snd addCharList) (fmap alphabet $ fmap snd addCharList)

        addCharInfo = V.singleton $ (snd $ head addCharList) {name = T.pack "CombinedAdditiveCharacters", origInfo = origAddData}
        
        -- Matrix Characters
        (matrixCharacters, matrixCharInfoList) = mergeMatrixCharacters matrixCharListList emptyCharacter

        -- Unchanged characters
        (unchangedCharacters, unchangeCharacterInfoList) = unzip unchangedCharList

        -- buildList incrementally
        newCharacterList' = [nonAddCharacter | not (null nonAddCharacter)]
        newCharacterList'' = if null addCharacter then newCharacterList'
                             else addCharacter : newCharacterList'
        newCharacterList''' = newCharacterList'' ++ matrixCharacters ++ unchangedCharacters

        newChararacterInfoList' = [nonAddCharInfo | not (null nonAddCharacter)]
        newChararacterInfoList'' = if null addCharacter then newChararacterInfoList'
                                  else addCharInfo : newChararacterInfoList'
        newChararacterInfoList''' = newChararacterInfoList'' ++ (fmap V.singleton matrixCharInfoList) ++ (fmap V.singleton unchangeCharacterInfoList)

    in
    {-
    trace ("Recoded Non-Additive: " ++ (show $ length nonAddCharList) ++ "->" ++ (show (length nonAddCharacter, fmap length $ fmap stateBVPrelim nonAddCharacter))
        ++ " Additive: " ++ (show $ length addCharList) ++ "->" ++ (show (length addCharacter, fmap length $ fmap rangePrelim addCharacter))
        ++ " Matrix " ++ (show  $length matrixCharListList) ++ "->" ++ (show $ length matrixCharacters)
        ++ " total list: " ++ (show (length newCharacterList''', fmap length newCharacterList''')) ++ " CI " ++ (show $ length newChararacterInfoList'''))
    -}
    (V.fromList $ V.fromList <$> L.transpose newCharacterList''', V.concat newChararacterInfoList''')


-- | combineMatrixCharacters cretes a series of lists of characters each of which has a different cost matrix
-- each character "type" (based on matrix) can have 1 or more characters
mergeMatrixCharacters :: [([[CharacterData]], CharInfo)] -> CharacterData -> ([[CharacterData]], [CharInfo])
mergeMatrixCharacters inMatrixCharListList charTemplate =
    -- should probably reverse the characters to maintian similar ordering to input
    let (charDataList, charInfoList) = unzip inMatrixCharListList
        combinedMatrixCharList = fmap (combineMatrixCharacters charTemplate []) charDataList
    in
    (combinedMatrixCharList, charInfoList)

-- | combineMatrixCharacters takes all matrix characters with same cost matrix and combines into
-- a single character with vector of original characters
combineMatrixCharacters :: CharacterData -> [[V.Vector MatrixTriple]] -> [[CharacterData]] -> [CharacterData]
combineMatrixCharacters charTemplate currentTripleList inMatrixCharDataList =
   if null inMatrixCharDataList then
      -- create character vector for preliminary states concatenating by taxon
      let taxRowCharList = L.transpose currentTripleList
          newCharacterData = fmap (makeMatrixCharacterList charTemplate) taxRowCharList
      in
      newCharacterData
   else
        -- first Character
        let charDataList = head inMatrixCharDataList
            prelimTripleList = fmap (V.head . matrixStatesPrelim) charDataList
        in
        combineMatrixCharacters charTemplate (prelimTripleList : currentTripleList) (tail inMatrixCharDataList)

-- | makeMatrixCharacterList takes a taxon list of matrix characters
-- and converts to single vector and makes new character for the taxon
makeMatrixCharacterList :: CharacterData -> [V.Vector MatrixTriple] -> CharacterData
makeMatrixCharacterList charTemplate tripleList = charTemplate {matrixStatesPrelim = V.fromList tripleList}

-- | combineNonAdditveCharacters takes a list of character data with singleton non-additive characters and puts
-- them together in a single character for each taxon
combineNonAdditveCharacters :: [([CharacterData], CharInfo)] -> CharacterData -> [[BV.BitVector]] -> [CharacterData]
combineNonAdditveCharacters nonAddCharList charTemplate currentBVList =
    if null nonAddCharList then
        -- create character vector for preliminary states concatenating by taxon
        -- single created and redone twice with prepend no need to reverse (that there really is anyway)
        let taxRowCharList = L.transpose currentBVList
            newCharacterData = fmap (makeNonAddCharacterList charTemplate) taxRowCharList
        in
        newCharacterData
    else
        -- first Character
        let (charDataList, _) = head nonAddCharList
            prelimBVList = fmap ((V.head . snd3) . stateBVPrelim) charDataList
        in
        combineNonAdditveCharacters (tail nonAddCharList) charTemplate (prelimBVList : currentBVList)

-- | combineAdditveCharacters takes a list of character data with singleton non-additive characters and puts
-- them together in a single character for each taxon
combineAdditveCharacters :: [([CharacterData], CharInfo)] -> CharacterData -> [[(Int, Int)]] -> [CharacterData]
combineAdditveCharacters addCharList charTemplate currentRangeList =
    if null addCharList then
        -- create character vector for preliminary states concatenating by taxon
        -- single created and redone twice with prepend no need to reverse (that there really is anyway)
        let taxRowCharList = L.transpose currentRangeList
            newCharacterData = fmap (makeAddCharacterList charTemplate) taxRowCharList
        in
        newCharacterData
    else
        -- first Character
        let (charDataList, _) = head addCharList
            prelimRangeList = fmap ((V.head . snd3) . rangePrelim) charDataList
        in
        combineAdditveCharacters (tail addCharList) charTemplate (prelimRangeList : currentRangeList)

-- | makeNonAddCharacterList takes a taxon list of characters
-- convertes chars to single vector and makes new character for the taxon
-- assumes a leaf so all fields same
makeNonAddCharacterList :: CharacterData -> [BV.BitVector] -> CharacterData
makeNonAddCharacterList charTemplate bvList = charTemplate {stateBVPrelim = (V.fromList bvList, V.fromList bvList,V.fromList bvList)}

-- | makeAddCharacterList takes a taxon list of characters
-- to single vector and makes new character for the taxon
-- assums a leaf so so all fields same
makeAddCharacterList :: CharacterData -> [(Int, Int)] -> CharacterData
makeAddCharacterList charTemplate rangeList = charTemplate {rangePrelim = (V.fromList rangeList, V.fromList rangeList, V.fromList rangeList)}

-- | addMatrixCharacter adds a matrix character to the appropriate (by cost matrix) list of matrix characters
-- replicates character by integer weight
addMatrixCharacter :: [([[CharacterData]], CharInfo)] -> S.Matrix Int -> ([CharacterData], CharInfo)-> Int -> [([[CharacterData]], CharInfo)]
addMatrixCharacter inMatrixCharacterList currentCostMatrix currentMatrixCharacter replicateNumber =
    if null inMatrixCharacterList then
        -- didn't find a match --so need to add new type to list of matrix character types
        if replicateNumber == 1 then
                [([fst currentMatrixCharacter], snd currentMatrixCharacter)]

            else
                [BF.first (replicate replicateNumber) currentMatrixCharacter]

    else
        let firstList@(firstMatrixCharList, localCharInfo) = head inMatrixCharacterList
            firstMatrix = costMatrix localCharInfo
        in

        -- matrices match--so correct matrix character type
        if firstMatrix == currentCostMatrix then
            if replicateNumber == 1 then
                (fst currentMatrixCharacter : firstMatrixCharList, localCharInfo) : tail inMatrixCharacterList

            else
                (replicate replicateNumber (fst currentMatrixCharacter) ++ firstMatrixCharList, localCharInfo) : tail inMatrixCharacterList

        -- matrices don't match so recurse to next one
        else firstList : addMatrixCharacter (tail inMatrixCharacterList) currentCostMatrix currentMatrixCharacter replicateNumber
-}
