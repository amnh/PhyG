{-# OPTIONS_GHC -Wno-missed-specialisations #-}

{- |
Module with functionality to transform phylogenetic data
-}
module Input.Reorganize (
    combineDataByType,
    reBlockData,
    removeConstantCharactersPrealigned,
    removeConstantCharsPrealigned,
    optimizePrealignedData,
    convert2BV,
    getRecodingType,
) where

import Bio.DynamicCharacter.Element (SlimState, WideState)
import Data.Alphabet
import Data.BitVector.LittleEndian qualified as BV
import Data.Bits
import Data.Functor (($>))
import Data.List qualified as L
import Data.List.NonEmpty (NonEmpty (..))
import Data.Maybe
import Data.Text.Lazy qualified as T
import Data.Text.Short qualified as ST
import Data.Vector qualified as V
import Data.Vector.Generic qualified as GV
import Data.Vector.Storable qualified as SV
import Data.Vector.Unboxed qualified as UV
import Debug.Trace
import GeneralUtilities
import GraphOptimization.Medians qualified as M
import Input.BitPack qualified as BP
import Measure.Transition qualified as Measure
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import SymMatrix qualified as S
import Text.Read
import TransitionMatrix.Metricity
import Types.Types
import Utilities.Utilities qualified as U


{- | optimizePrealignedData convert
prealigned to non-additive or matrix
here
bitPack new non-additive
packNonAdditive
-}
optimizePrealignedData ∷ GlobalSettings → ProcessedData → PhyG ProcessedData
optimizePrealignedData inGS inData@(_, _, blockDataVect)
    | U.getNumberPrealignedCharacters blockDataVect == 0 = logWith LogMore "Not bitpacking..." $> inData
    | otherwise = do
        logWith LogMore "Bitpacking..."

        -- remove constant characters from prealigned
        inData' ← removeConstantCharactersPrealigned inData

        -- convert prealigned to nonadditive if all 1 tcms
        let inData'' = convertPrealignedToNonAdditive inData'

        -- bit packing for non-additivecharacters
        BP.packNonAdditiveData inGS inData''


{- | convertPrealignedToNonAdditive converts prealigned data to non-additive
if homogeneous TCM (all 1's non-diagnoal)
-}
convertPrealignedToNonAdditive ∷ ProcessedData → ProcessedData
convertPrealignedToNonAdditive (nameVect, bvNameVect, blockDataVect) = (nameVect, bvNameVect, fmap convertPrealignedToNonAdditiveBlock blockDataVect)


{- | convertPrealignedToNonAdditiveBlock takes a character block and convertes prealigned to non-add if tcms all 1's
this is done taxon by taxon and character by character since can convert with only local infomation
-}
convertPrealignedToNonAdditiveBlock ∷ BlockData → BlockData
convertPrealignedToNonAdditiveBlock (nameBlock, charDataVV, charInfoV) =
    let codingTypeV = (metricity <$> charInfoV)
        (newCharDataVV, newCharInfoVV) = V.unzip $ fmap (convertTaxonPrealignedToNonAdd charInfoV codingTypeV) charDataVV
    in  (nameBlock, newCharDataVV, V.head newCharInfoVV)


{- | convertTaxonPrealignedToNonAdd takes a vector of character info and vector of charcter data for a taxon
and recodes prealigned as non-additve if tcm's all 1s
-}
convertTaxonPrealignedToNonAdd
    ∷ V.Vector CharInfo → V.Vector Metricity → V.Vector CharacterData → (V.Vector CharacterData, V.Vector CharInfo)
convertTaxonPrealignedToNonAdd charInfoV codingTypeV charDataV =
    V.unzip $ V.zipWith3 convertTaxonPrealignedToNonAddCharacter charInfoV codingTypeV charDataV


{- | convertTaxonPrealignedToNonAddCharacter takes a taxon character and char info and cost matrix type
and transforms to non-additive if all tcms are 1's
-}
convertTaxonPrealignedToNonAddCharacter ∷ CharInfo → Metricity → CharacterData → (CharacterData, CharInfo)
convertTaxonPrealignedToNonAddCharacter charInfo matrixType charData
    | charType charInfo `notElem` prealignedCharacterTypes = (charData, charInfo)
    | otherwise =
        let charWeight = weight charInfo
            newStateBV = case charType charInfo of
                AlignedSlim → convert2BVTriple 32 $ (snd3 . alignedSlimPrelim) charData
                AlignedWide → convert2BVTriple 64 $ (snd3 . alignedWidePrelim) charData
                AlignedHuge → error "No point in converting huge--can't pack anyway"
                _ →
                    error $
                        "Unrecognized character type in convertTaxonPrealignedToNonAddCharacter: "
                            <> show (charType charInfo)

            modifiedChar = emptyCharacter{stateBVPrelim = newStateBV}
            modifiedType =
                charInfo
                    { charType = NonAdd
                    , weight = charWeight
                    }
        in  case matrixType of
                Special (DiscreteMetric _) → (modifiedChar, modifiedType)
                _ → (charData, charInfo)


-- | convert2BVTriple takes SlimState or WideState and converts to Triple Vector of bitvectors
convert2BVTriple
    ∷ (Integral a, GV.Vector v a) ⇒ Word → v a → (V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector)
convert2BVTriple size inM =
    let inMList = GV.toList inM
        inMBV = fmap (BV.fromNumber size) inMList
    in  (V.fromList inMBV, V.fromList inMBV, V.fromList inMBV)


{- | convert2BV takes SlimState or WideState and converts to Vector of bitvectors
this for leaves so assume M only one needed really
-}
convert2BV
    ∷ (Integral a, GV.Vector v a) ⇒ Word → (v a, v a, v a) → (V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector)
convert2BV size (_, inM, _) =
    let inMList = GV.toList inM
        inMBV = fmap (BV.fromNumber size) inMList
    in  (V.fromList inMBV, V.fromList inMBV, V.fromList inMBV)


{- | getRecodingType takes a cost matrix and detemines if it can be recodes as non-additive,
non-additive with gap chars, or matrix
assumes indel costs are in first row and column
-}
getRecodingType ∷ S.Matrix Int → (String, Int)
getRecodingType inMatrix =
    if S.null inMatrix
        then error "Null matrix in getRecodingType"
        else
            if (not . S.isSymmetric) inMatrix
                then ("matrix", 0)
                else
                    let matrixLL = S.toFullLists inMatrix
                        secondRow = matrixLL L.!! 1
                        numUniqueCosts = length $ L.group $ L.sort $ (filter (/= 0) $ concat matrixLL)
                    in  -- trace  ("GRT: " <> (show numUniqueCosts)) (

                        -- don't recode huge--no point
                        if S.rows inMatrix > 32
                            then ("matrix", head secondRow)
                            else -- all same except for 0

                                if numUniqueCosts == 1
                                    then ("nonAdd", 0)
                                    else ("matrix", head secondRow)


{-
-- all same except for gaps
else if numUniqueCosts == 2 then
   trace ("NAG: " <> (show $ length $ L.group $ L.sort $ filter (/= 0) lastRow)) (
   if  (length $ L.group $ filter (/= 0) lastRow) == 1 then ("nonAddGap", head lastRow)

   -- some no gaps different
   else ("matrix",  head lastRow)
   )
-- to many types for nonadd coding
else ("matrix",  head lastRow)

-}
-- )

{- | reBlockData takes original block assignments--each input file is a block--
and combines, creates new, deletes empty blocks from user input
reblock pair names may contain wildcards
-}
reBlockData ∷ [(NameText, NameText)] → ProcessedData → PhyG ProcessedData
reBlockData reBlockPairs inData@(leafNames, leafBVs, blockDataV)
    | null reBlockPairs = logWith LogInfo "Character Blocks as input files\n" $> inData
    | otherwise -- those block to be reassigned--nub in case repeated names
        =
        let toBeReblockedNames = fmap (T.filter (/= '"')) $ L.nub $ fmap snd reBlockPairs
            unChangedBlocks = V.filter ((`notElemWildcards` toBeReblockedNames) . fst3) blockDataV
            blocksToChange = V.filter ((`elemWildcards` toBeReblockedNames) . fst3) blockDataV
            newBlocks = makeNewBlocks reBlockPairs blocksToChange []
            reblockedBlocks = unChangedBlocks <> V.fromList newBlocks
            message =
                unwords
                    [ "\nReblocking:"
                    , show toBeReblockedNames
                    , "leaving unchanged:"
                    , show (fmap fst3 unChangedBlocks)
                    , "\n\tNew blocks:"
                    , show (fmap fst3 reblockedBlocks)
                    , "\n"
                    ]
        in  logWith LogInfo message $> (leafNames, leafBVs, reblockedBlocks)


-- | makeNewBlocks takes lists of reblock pairs and existing relevant blocks and creates new blocks returned as a list
makeNewBlocks ∷ [(NameText, NameText)] → V.Vector BlockData → [BlockData] → [BlockData]
makeNewBlocks reBlockPairs inBlockV curBlockList
    | null reBlockPairs = curBlockList
    | V.null inBlockV && null curBlockList =
        errorWithoutStackTrace
            ( "Reblock pair names do not have a match for any input block--perhaps missing '#0/N'? Blocks: " <> (show $ fmap snd reBlockPairs)
            )
    | V.null inBlockV = curBlockList
    | otherwise =
        let firstBlock = V.head inBlockV
            firstName = fst3 firstBlock
            newPairList = fst <$> filter (textMatchWildcards firstName . snd) reBlockPairs
        in  if null newPairList
                then
                    errorWithoutStackTrace
                        ( "Reblock pair names do not have a match for any input block--perhaps missing '#0'? Specified pairs: "
                            <> show reBlockPairs
                            <> " input block name: "
                            <> T.unpack firstName
                        )
                else
                    if length newPairList > 1
                        then errorWithoutStackTrace ("Multiple reblock destinations for single input block" <> show newPairList)
                        else
                            let newBlockName = head newPairList
                                existingBlock = filter ((== newBlockName) . fst3) curBlockList
                            in  -- new block to be created
                                if null existingBlock
                                    then -- trace("NBlocks:" <> (show $ fmap fst3 curBlockList))
                                        makeNewBlocks reBlockPairs (V.tail inBlockV) ((newBlockName, snd3 firstBlock, thd3 firstBlock) : curBlockList)
                                    else -- existing block to be added to

                                        if length existingBlock > 1
                                            then error ("Error: Block to be added to more than one block-should not happen: " <> show reBlockPairs)
                                            else -- need to add character vectors to vertex vectors and add to CharInfo
                                            -- could be multiple  'characteres' if non-exact data in inpuit file (Add, NonAdd, MAtrix etc)

                                                let blockToAddTo = head existingBlock
                                                    newCharData = V.zipWith (<>) (snd3 blockToAddTo) (snd3 firstBlock)
                                                    newCharInfo = thd3 blockToAddTo <> thd3 firstBlock
                                                in  -- trace("EBlocks:" <> (show $ fmap fst3 curBlockList))
                                                    makeNewBlocks
                                                        reBlockPairs
                                                        (V.tail inBlockV)
                                                        ((newBlockName, newCharData, newCharInfo) : filter ((/= newBlockName) . fst3) curBlockList)


{- | combineDataByType combines data of same type for (exact types) into
same vectors in character so have one non-add, one add, one of each packed type,
can have multiple matrix (due to cost matrix differneces)
similar result to groupDataByType, but does not assume single characters.
-}
combineDataByType ∷ GlobalSettings → ProcessedData → PhyG ProcessedData
combineDataByType inGS inData@(taxNames, taxBVNames, _) =
    -- recode add to non-add before combine-- takes care wor integer weighting
    let (_, _, blockDataV') = recodeAddToNonAddCharacters inGS maxAddStatesToRecode inData
    in  do
            recodedData ← mapM combineData blockDataV'
            pure (taxNames, taxBVNames, recodedData)


-- | combineData creates for a block) lists of each data type and concats then creating new data and new char info
combineData ∷ BlockData → PhyG BlockData
combineData (blockName, blockDataVV, charInfoV) =
    let -- parallel setup
        action ∷ V.Vector CharacterData → (V.Vector CharacterData, V.Vector CharInfo)
        action = combineBlockData charInfoV
    in  do
            pTraverse ← getParallelChunkMap
            let result = pTraverse action (V.toList blockDataVV)
            let (newBlockDataLV, newCharInfoLV) = unzip result
            -- (newBlockDataLV, newCharInfoLV) = unzip (PU.seqParMap PU.myStrategyRDS (combineBlockData charInfoV) (V.toList blockDataVV)) -- `using` PU.myParListChunkRDS)
            pure (blockName, V.fromList newBlockDataLV, head newCharInfoLV)


{- | combineBlockData takes a vector of char info and vector or charcater data for a taxon and
combined exact data types into single characters
additive characters should have already been converted (if less that maxAddStatesToRecode) to binary
non-additive characters with integer weights could be  repeated before combining-- but this has been disabled
due to memory issues in large data sets
non-integer additive and non-additive are grouped together by weight so they can be combined and bit packed
all other types are grouped by weight for efficiency of optimization and reductionn of multiplies
zero weight characters are filtered out
-}
combineBlockData ∷ V.Vector CharInfo → V.Vector CharacterData → (V.Vector CharacterData, V.Vector CharInfo)
combineBlockData inCharInfoV inCharDataV =
    let pairCharsInfo = V.zip inCharInfoV inCharDataV

        -- characters to not be reorganized-- nbasically the sequence characters
        sequenceCharacters = V.toList $ V.filter ((> 0) . weight . fst) $ V.filter ((`elem` sequenceCharacterTypes) . charType . fst) pairCharsInfo

        -- matrix characters are more complex--can only join if same matrix
        matrixCharsPair = V.filter ((> 0) . weight . fst) $ V.filter ((== Matrix) . charType . fst) pairCharsInfo
        (newMatrixCharL, newMatrixCharInfoL) =
            if (not . null) matrixCharsPair
                then unzip $ organizeMatrixCharsByMatrix (V.toList matrixCharsPair)
                else ([], [])

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
        -- newCharList = newNonAddCharL <> (V.toList nonAddCharsWeightNotInt) <> newAddCharL <> (V.toList addCharsWeightNot1) <> newPacked2CharL <> newPacked4CharL <> newPacked5CharL <> newPacked8CharL <> newPacked64CharL <> newMatrixCharL <> (fmap snd sequenceCharacters)
        -- newCharInfoList = newNonAddCharInfoL <> (V.toList nonAddCharsWeightNotIntInfo) <> newAddCharInfoL <> (V.toList addCharsWeightNot1Info) <> newPacked2CharInfoL <> newPacked4CharInfoL <> newPacked5CharInfoL <> newPacked8CharInfoL <> newPacked64CharInfoL <> newMatrixCharInfoL <> (fmap fst sequenceCharacters)

        newCharList =
            newNonAddChar
                <> newAddChar
                <> newPacked2Char
                <> newPacked4Char
                <> newPacked5Char
                <> newPacked8Char
                <> newPacked64Char
                <> newMatrixCharL
                <> (fmap snd sequenceCharacters)
        newCharInfoList =
            newNonAddCharInfo
                <> newAddCharInfo
                <> newPacked2CharInfo
                <> newPacked4CharInfo
                <> newPacked5CharInfo
                <> newPacked8CharInfo
                <> newPacked64CharInfo
                <> newMatrixCharInfoL
                <> (fmap fst sequenceCharacters)
    in  (V.fromList newCharList, V.fromList newCharInfoList)


{- | groupCharactersByWeight takes a list of characters and returns a list of lists of charcter with same weight
checked as Double.
NB--Does not check for character type---assuming this is desired it must be assured before input
-}
groupCharactersByWeight ∷ V.Vector (CharInfo, CharacterData) → V.Vector (CharInfo, CharacterData)
groupCharactersByWeight inCharsPairList =
    if V.null inCharsPairList
        then V.empty
        else
            let weightList = L.nub $ V.toList $ fmap weight (fmap fst inCharsPairList)
                charListByWeight = fmap (getSameWeightChars inCharsPairList) $ V.fromList weightList
            in  V.concatMap mergeCharacters charListByWeight


{- | mergeCharacters merges the data fieed of characters based on type
NB--Does not check all chars are same type or weight--will silently combine mempty values
with whatever weight.
returns list with single member so can concat later
-}
mergeCharacters ∷ V.Vector (CharInfo, CharacterData) → V.Vector (CharInfo, CharacterData)
mergeCharacters inCharsPairList =
    if V.null inCharsPairList
        then V.empty
        else
            let thisCharType = (charType . fst . V.head) inCharsPairList

                -- non-add data
                dataFieldNonAdd = V.concatMap (snd3 . stateBVPrelim . snd) inCharsPairList
                newNonAddChar =
                    ((snd . V.head) inCharsPairList)
                        { stateBVPrelim = (dataFieldNonAdd, dataFieldNonAdd, dataFieldNonAdd)
                        , stateBVFinal = dataFieldNonAdd
                        }

                -- add data
                dataFieldAdd = V.concatMap (snd3 . rangePrelim . snd) inCharsPairList
                newAddChar = ((snd . V.head) inCharsPairList){rangePrelim = (dataFieldAdd, dataFieldAdd, dataFieldAdd), rangeFinal = dataFieldAdd}

                -- packed data
                dataFieldPacked = UV.concat $ V.toList $ fmap (snd3 . packedNonAddPrelim . snd) inCharsPairList
                newPackedChar =
                    ((snd . V.head) inCharsPairList)
                        { packedNonAddPrelim = (dataFieldPacked, dataFieldPacked, dataFieldPacked)
                        , packedNonAddFinal = dataFieldPacked
                        }

                -- add info of merging to character info
                newOrigInfo = V.concatMap (origInfo . fst) inCharsPairList
                newCharInfo = ((fst . V.head) inCharsPairList){origInfo = newOrigInfo}
            in  if thisCharType == NonAdd
                    then V.singleton (newCharInfo, newNonAddChar)
                    else
                        if thisCharType == Add
                            then V.singleton (newCharInfo, newAddChar)
                            else
                                if thisCharType `elem` packedNonAddTypes
                                    then V.singleton (newCharInfo, newPackedChar)
                                    else error ("Error in mergeCharacters: Character type " <> show thisCharType <> " unrecognized/not implemented")


-- | getSameWeightChars returns character pairs with same matrix as testMatrix
getSameWeightChars ∷ V.Vector (CharInfo, CharacterData) → Double → V.Vector (CharInfo, CharacterData)
getSameWeightChars inCharsPairList testWeight =
    if V.null inCharsPairList
        then V.empty
        else
            let inWeightList = fmap weight (fmap fst inCharsPairList)
                weightPairPair = V.zip inWeightList inCharsPairList
                matchList = V.filter ((== testWeight) . weight . fst . snd) weightPairPair
            in  fmap snd matchList


{- Currently unused--employed to reduce chrater umbers by replicating characters with integer weight
    can cause memory issues with large data sets
-- | replicateCharPairByWeight replicates characters by integer weight
replicateCharPairByWeight :: (CharInfo, CharacterData) -> [(CharInfo, CharacterData)]
replicateCharPairByWeight firstPair =
    let charIntWeight = doubleAsInt $ (weight . fst) firstPair
    in
    if isNothing charIntWeight then error ("Character weight not an integer in replicateCharPair: " <> (show $ (weight . fst) firstPair))
    else
        (replicate (fromJust charIntWeight) firstPair)
-}

-- | organizeMatrixCharsByMatrix combines matrix charcters if they have the same cost matrix
organizeMatrixCharsByMatrix ∷ [(CharInfo, CharacterData)] → [(CharacterData, CharInfo)]
organizeMatrixCharsByMatrix inCharsPairList =
    if null inCharsPairList
        then []
        else
            let costMatrixList = L.nub $ fmap costMatrix (fmap fst inCharsPairList)
                charMatrixLL = fmap (getSameMatrixChars inCharsPairList) costMatrixList
                newMatrixPairs = fmap combineMatrixCharsByMatrix charMatrixLL
            in  newMatrixPairs


-- | combineMatrixCharsByMatrix combines matrix characters--assumes cost matrices are the same
combineMatrixCharsByMatrix ∷ [(CharInfo, CharacterData)] → (CharacterData, CharInfo)
combineMatrixCharsByMatrix inCharList =
    let newMatrixcharData = V.concat $ fmap matrixStatesPrelim $ fmap snd inCharList
        newMatrixChar =
            ((snd . head) inCharList)
                { matrixStatesPrelim = newMatrixcharData
                , matrixStatesFinal = newMatrixcharData
                }
        newMatrixCharInfo = ((fst . head) inCharList){origInfo = V.concat $ fmap (origInfo . fst) inCharList}
    in  (newMatrixChar, newMatrixCharInfo)


-- | getSameMatrixChars returns character pairs with same matrix as testMatrix
getSameMatrixChars ∷ [(CharInfo, CharacterData)] → S.Matrix Int → [(CharInfo, CharacterData)]
getSameMatrixChars inCharsPairList testMatrix =
    let inMatrixList = fmap costMatrix (fmap fst inCharsPairList)
        matrixPairPair = zip inMatrixList inCharsPairList
        matchList = filter ((== testMatrix) . costMatrix . fst . snd) matrixPairPair
    in  fmap snd matchList


{- | removeConstantCharactersPrealigned takes processed data and removes constant characters
from prealignedCharacterTypes
-}
removeConstantCharactersPrealigned ∷ ProcessedData → PhyG ProcessedData
removeConstantCharactersPrealigned (nameVect, bvNameVect, blockDataVect) =
    let -- parallel setup
        action ∷ BlockData → PhyG BlockData
        action = removeConstantBlockPrealigned
    in  do
            pTraverse ← getParallelChunkTraverse
            newBlockData <- pTraverse action blockDataVect
            pure (nameVect, bvNameVect, newBlockData)


-- | removeConstantBlockPrealigned takes block data and removes constant characters
removeConstantBlockPrealigned ∷ BlockData → PhyG BlockData
removeConstantBlockPrealigned inBlockData@(blockName, taxVectByCharVect, charInfoV)
    -- check for null data--really really shouldn't happen
    | V.null taxVectByCharVect = logWith LogWarn "Null block data in removeConstantBlockPrealigned" $> inBlockData
    -- check for prealigned data in block
    | U.getNumberPrealignedCharacters (V.singleton inBlockData) == 0 = pure inBlockData
    -- standard case for removal
    | otherwise =
        let numChars = V.length $ V.head taxVectByCharVect
            -- create vector of single characters with vector of taxon data of sngle character each
            -- like a standard matrix with a single character
            singleCharVect = fmap (U.getSingleCharacter taxVectByCharVect) $ V.fromList [0 .. numChars - 1]
        in  do  -- actually remove constants form chaarcter list
                singleCharVect' <- V.zipWithM removeConstantCharsPrealigned singleCharVect charInfoV
                -- recreate the taxa vext by character vect block data expects
                -- should filter out length zero characters
                let newTaxVectByCharVect = U.glueBackTaxChar singleCharVect'
                pure (blockName, newTaxVectByCharVect, charInfoV)


{- | removeConstantCharsPrealigned takes a single 'character' and if proper type removes if all values are the same
could be done if character has max lenght of 0 as well.
packed types already filtered when created
-}
removeConstantCharsPrealigned ∷ V.Vector CharacterData → CharInfo → PhyG (V.Vector CharacterData)
removeConstantCharsPrealigned singleChar charInfo =
    let inCharType = charType charInfo
    in  -- dynamic characters don't do this
        if inCharType `notElem` prealignedCharacterTypes
            then pure singleChar
            else getVariableChars inCharType singleChar


{- | getVariableChars checks identity of states in a vector positin in all taxa
and returns True if variable, False if constant
bit packed and non-exact should not get in here
-}
getVariableChars ∷ CharType → V.Vector CharacterData → PhyG (V.Vector CharacterData)
getVariableChars inCharType singleChar =
    let nonAddV = snd3 . stateBVPrelim <$> singleChar
        addV = snd3 . rangePrelim <$> singleChar
        alSlimV = snd3 . alignedSlimPrelim <$> singleChar
        alWideV = snd3 . alignedWidePrelim <$> singleChar
        alHugeV = snd3 . alignedHugePrelim <$> singleChar
        matrixV = matrixStatesPrelim <$> singleChar

        -- get identity vect
        boolVar =
            if inCharType == NonAdd
                then getVarVectBits inCharType nonAddV []
                else
                    if inCharType == Add
                        then getVarVectAdd addV []
                        else
                            if inCharType == Matrix
                                then getVarVectMatrix matrixV []
                                else
                                    if inCharType == AlignedSlim
                                        then getVarVectBits inCharType alSlimV []
                                        else
                                            if inCharType == AlignedWide
                                                then getVarVectBits inCharType alWideV []
                                                else
                                                    if inCharType == AlignedHuge
                                                        then getVarVectBits inCharType alHugeV []
                                                        else error ("Char type unrecognized in getVariableChars: " <> show inCharType)

        boolVar' = V.fromList boolVar

        -- get Variable characters by type
        nonAddVariable = filterConstantsV boolVar' <$> nonAddV
        addVariable = filterConstantsV boolVar' <$> addV
        matrixVariable = filterConstantsV boolVar' <$> matrixV
        alSlimVariable = filterConstantsSV boolVar' <$> alSlimV
        alWideVariable = filterConstantsUV boolVar' <$> alWideV

        -- this is a hack-not sure why popCount etc don't work with HugeVector little endian BVs
        -- these should be short seqs in general so no ral impact on effeciancy
        alHugeVariable = alHugeV -- fmap (filterConstantsV (V.fromList boolVar)) alHugeV

        -- assign to proper character fields
    in  V.zipWithM (assignNewField inCharType) singleChar $
            V.zip6 nonAddVariable addVariable matrixVariable alSlimVariable alWideVariable alHugeVariable


{- | getVarVectAdd takes a vector of a vector additive ranges and returns False if range overlap
True if not (short circuits)
based on range overlap
-}
getVarVectAdd ∷ V.Vector (V.Vector (Int, Int)) → [Bool] → [Bool]
getVarVectAdd stateVV curBoolList =
    if V.null (V.head stateVV)
        then L.reverse curBoolList
        else
            let firstChar = fmap V.head stateVV
                isVariable = checkIsVariableAdditive (V.head firstChar) (V.tail firstChar)
            in  getVarVectAdd (fmap V.tail stateVV) (isVariable : curBoolList)


{- | getVarVectMatrix takes a generic vector and returns False if values are same
True if not (short circuits)
based on simple identity not max cost zero
-}
getVarVectMatrix ∷ V.Vector (V.Vector (V.Vector MatrixTriple)) → [Bool] → [Bool]
getVarVectMatrix stateVV curBoolList =
    if V.null (V.head stateVV)
        then L.reverse curBoolList
        else
            let firstChar = fmap V.head stateVV
                isVariable = checkIsVariableMatrix (getMatrixStateList $ V.head firstChar) (V.tail firstChar)
            in  getVarVectMatrix (fmap V.tail stateVV) (isVariable : curBoolList)


{- |
getVarVectBits takes a generic vector and returns False if values are same
True if not (short circuits)
based on simple identity not max cost zero
-}
getVarVectBits ∷ (FiniteBits a, GV.Vector v a) ⇒ CharType → V.Vector (v a) → [Bool] → [Bool]
getVarVectBits inCharType stateVV curBoolList =
    if GV.null (V.head stateVV)
        then L.reverse curBoolList
        else
            let firstChar = fmap GV.head stateVV
                isVariable =
                    if inCharType == NonAdd
                        then checkIsVariableBit (GV.head firstChar) (GV.tail firstChar)
                        else
                            if inCharType == AlignedSlim
                                then checkIsVariableBit (GV.head firstChar) (GV.tail firstChar)
                                else
                                    if inCharType == AlignedWide
                                        then checkIsVariableBit (GV.head firstChar) (GV.tail firstChar)
                                        else
                                            if inCharType == AlignedHuge
                                                then checkIsVariableBit (GV.head firstChar) (GV.tail firstChar)
                                                else error ("Char type unrecognized in getVariableChars: " <> show inCharType)
            in  getVarVectBits inCharType (fmap GV.tail stateVV) (isVariable : curBoolList)


{-
-- | checkIsVariableIdentity takes a generic vector and sees if all elements are identical
checkIsVariableIdentity ::  (Eq a, GV.Vector v a) => a -> v a -> Bool
checkIsVariableIdentity firstElement inVect =
    if GV.null inVect then False
    else
        if firstElement /= GV.head inVect then True
        else checkIsVariableIdentity firstElement (GV.tail inVect)
-}

{- | checkIsVariableAdditive checks if additive charcter is variable
by taking new ranges and range costs of first element with all others
if summed cost > 0 then variable
-}
checkIsVariableAdditive ∷ (Int, Int) → V.Vector (Int, Int) → Bool
checkIsVariableAdditive (ir1, ir2) rangeList =
    if V.null rangeList
        then False
        else
            let (nr1, nr2) = V.head rangeList
                (newMin, newMax, newCost) = M.getNewRange ir1 ir2 nr1 nr2
            in  if newCost > 0
                    then True
                    else checkIsVariableAdditive (newMin, newMax) (V.tail rangeList)


-- | getMatrixStateList returns minimum matrix characters states as integers
getMatrixStateList ∷ V.Vector MatrixTriple → [Int]
getMatrixStateList inState =
    let statePairList = zip (fmap fst3 $ V.toList inState) [0 .. (V.length inState - 1)]
        minCost = minimum $ fmap fst3 $ V.toList inState
        stateList = fmap snd $ filter ((== minCost) . fst) statePairList
    in  stateList


{- | checkIsVariableMatrix checks if matrix charcter is variable
by taking min cost states of first element and
checks for overlap--if empty list intersect then variable
-}
checkIsVariableMatrix ∷ [Int] → V.Vector (V.Vector MatrixTriple) → Bool
checkIsVariableMatrix inStateList restStatesV =
    if V.null restStatesV
        then False
        else
            let nextStateList = getMatrixStateList $ V.head restStatesV

                newStateList = L.intersect inStateList nextStateList
            in  if null newStateList
                    then True
                    else checkIsVariableMatrix newStateList (V.tail restStatesV)


{- | checkIsVariableBit takes a generic vector and checks for
state overlap via bit AND (.&.)
-}
checkIsVariableBit ∷ (FiniteBits a, GV.Vector v a) ⇒ a → v a → Bool
checkIsVariableBit firstElement restVect =
    if GV.null restVect
        then False
        else
            let newState = firstElement .&. (GV.head restVect)
            in  if popCount newState == 0
                    then True
                    else checkIsVariableBit newState (GV.tail restVect)


-- | these need to be abstracted but had problems with the bool list -> Generic vector, and SV pair

{- | filerConstantsV takes the charcter data and filters out teh constants
uses filter to keep O(n)
filterConstantsV :: (GV.Vector v a) => [Bool] -> v a -> v a
-}
filterConstantsV ∷ V.Vector Bool → V.Vector a → V.Vector a
filterConstantsV inVarBoolV charVect =
    let pairVect = V.zip charVect inVarBoolV
        variableCharV = V.map fst $ V.filter ((== True) . snd) pairVect
    in  variableCharV


{- | filerConstantsSV takes the charcter data and filters out teh constants
uses filter to keep O(n)
filterConstantsV :: (GV.Vector v a) => [Bool] -> v a -> v a
-}
filterConstantsSV ∷ (SV.Storable a) ⇒ V.Vector Bool → SV.Vector a → SV.Vector a
filterConstantsSV inVarBoolV charVect =
    let varVect = filterConstantsV inVarBoolV (V.fromList $ SV.toList charVect)
    in  SV.fromList $ V.toList varVect


{- | filerConstantsUV takes the charcter data and filters out teh constants
uses filter to keep O(n)
filterConstantsV :: (GV.Vector v a) => [Bool] -> v a -> v a
-}
filterConstantsUV ∷ (UV.Unbox a) ⇒ V.Vector Bool → UV.Vector a → UV.Vector a
filterConstantsUV inVarBoolV charVect =
    let varVect = filterConstantsV inVarBoolV (V.fromList $ UV.toList charVect)
    in  UV.fromList $ V.toList varVect


{- | assignNewField takes character type and a 6-tuple of charcter fields and assigns the appropriate
to the correct field
neither bit packed nor nno-exact should het here
-}
assignNewField
    ∷ CharType
    → CharacterData
    → ( V.Vector BV.BitVector
      , V.Vector (Int, Int)
      , V.Vector (V.Vector MatrixTriple)
      , SV.Vector SlimState
      , UV.Vector WideState
      , V.Vector BV.BitVector
      )
    → PhyG CharacterData
assignNewField inCharType charData (nonAddData, addData, matrixData, alignedSlimData, alignedWideData, alignedHugeData) = case inCharType of
    NonAdd → pure charData{stateBVPrelim = (nonAddData, nonAddData, nonAddData)}
    Add → pure charData{rangePrelim = (addData, addData, addData)}
    Matrix → pure charData{matrixStatesPrelim = matrixData}
    AlignedSlim → pure charData{alignedSlimPrelim = (alignedSlimData, alignedSlimData, alignedSlimData)}
    AlignedWide → pure charData{alignedWidePrelim = (alignedWideData, alignedWideData, alignedWideData)}
    AlignedHuge → pure charData{alignedHugePrelim = (alignedHugeData, alignedHugeData, alignedHugeData)}
    val → failWithPhase Parsing $ "Char type unrecognized in assignNewField: " <> show val


{- | recodeAddToNonAddCharacters takes an max states number and processsed data
and recodes additive characters with max state < input max (0..input max - 1)
as a series of binary non-additive characters
-}
recodeAddToNonAddCharacters ∷ GlobalSettings → Int → ProcessedData → ProcessedData
recodeAddToNonAddCharacters inGS maxStateToRecode (nameVect, nameBVVect, blockDataVect) =
    let newBlockDataVect = fmap (convertAddToNonAddBlock inGS maxStateToRecode) blockDataVect
    in  (nameVect, nameBVVect, newBlockDataVect)


-- | convertAddToNonAddBlock converts additive characters to non-additive in a block
convertAddToNonAddBlock ∷ GlobalSettings → Int → BlockData → BlockData
convertAddToNonAddBlock inGS maxStateToRecode (blockName, taxByCharDataVV, charInfoV) =
    let (newTaxByCharDataVV, newCharInfoVV) = V.unzip $ fmap (recodeTaxonData inGS maxStateToRecode charInfoV) taxByCharDataVV
    in  -- trace ("CNAB: " <> (show (V.length $ V.head newTaxByCharDataVV, V.length $ V.head newCharInfoVV)))
        (blockName, newTaxByCharDataVV, V.head newCharInfoVV)


-- | recodeTaxonData recodes Add as nonAdd for each taxon in turn
recodeTaxonData
    ∷ GlobalSettings → Int → V.Vector CharInfo → V.Vector CharacterData → (V.Vector CharacterData, V.Vector CharInfo)
recodeTaxonData inGS maxStateToRecode charInfoV taxonCharacterDataV =
    let (newCharDataVV, newCharInfoVV) = unzip $ zipWith (recodeAddToNonAddCharacter inGS maxStateToRecode) (V.toList taxonCharacterDataV) (V.toList charInfoV)
    in  -- trace ("RTD: " <> (show (V.length $ V.concat newCharDataVV, V.length $ V.concat newCharInfoVV)))
        (V.concat newCharDataVV, V.concat newCharInfoVV)


{- | recodeAddToNonAddCharacter takes a single character for single taxon and recodes if non-additive with
 fewer than maxStateToRecode states.
 assumes states in linear order
 replicatee charinfo for multiple new characters after recoding
-}
recodeAddToNonAddCharacter ∷ GlobalSettings → Int → CharacterData → CharInfo → (V.Vector CharacterData, V.Vector CharInfo)
recodeAddToNonAddCharacter inGS maxStateToRecode inCharData inCharInfo =
    let inCharType = charType inCharInfo
        numStates = 1 + (L.maximum $ fmap makeInt $ alphabet inCharInfo) -- min 2 (1 + (L.last $ L.sort $ fmap makeInt $ alphabetSymbols $ alphabet inCharInfo))
        -- numStates = 1 + (L.last $ L.sort $ fmap makeInt $ alphabetSymbols $ alphabet inCharInfo)
        origName = name inCharInfo
    in  -- if a single state recodes to a single uninfomative binary
        -- removed || ((not . doubleIsInt . weight) inCharInfo)  to allow for recodding (leaving weight) for non-integer weights
        if (inCharType /= Add)
            then (V.singleton inCharData, V.singleton inCharInfo)
            else -- the limit on recoded states is removed for PMDL/ML since otherwise bit costs will be incorrect

                if (numStates > maxStateToRecode) && ((optimalityCriterion inGS) `notElem` [PMDL, SI, MAPA])
                    then (V.singleton inCharData, V.singleton inCharInfo)
                    else
                        if numStates < 2
                            then (V.empty, V.empty)
                            else -- create numStates - 1 no-additve chaaracters (V.singleton inCharData, V.singleton inCharInfo)
                            -- bits ON-- [0.. snd range]

                                let minRangeIndex = fst $ V.head $ snd3 $ rangePrelim inCharData
                                    maxRangeIndex = snd $ V.head $ snd3 $ rangePrelim inCharData
                                    inCharOrigData = origInfo inCharInfo
                                    newCharInfo =
                                        inCharInfo
                                            { name = (T.pack $ (T.unpack origName) <> "RecodedToNonAdd")
                                            , charType = NonAdd
                                            , alphabet = fromSymbols . fmap ST.fromString $ "0" :| ["1"]
                                            , origInfo = inCharOrigData
                                            }

                                    -- create new characters and new character info
                                    newCharList = fmap (makeNewNonAddChar minRangeIndex maxRangeIndex) [0 .. numStates - 2]

                                    newCharInfoList = replicate (numStates - 1) newCharInfo
                                in  -- trace ("RTNA: Numstates " <> (show numStates) <> " " <> (show $ (snd3 . rangePrelim) inCharData) <> " -> " <> (show $ fmap (snd3 . stateBVPrelim) newCharList))
                                    -- (show (length newCharList, V.length $ V.replicate (numStates - 1) newCharInfo)) <> "\n" <> (show newCharList) <> "\n" <> (show $ charType newCharInfo))
                                    (V.fromList newCharList, V.fromList newCharInfoList)
    where
        makeInt a =
            let newA = readMaybe (ST.toString a) ∷ Maybe Int
            in  if isNothing newA
                    then error ("State '" <> (ST.toString a) <> "' not recoding to Int")
                    else fromJust newA


{- | makeNewNonAddCharacter takes a stateIndex and charcatear number
and makes a non-additive character with 0 or 1 coding
based on stateIndex versus state number
if stateIndex > charNumber then 1 else 0 (coded as bit 0 for 0, bit 1 for 1)
-}
makeNewNonAddChar ∷ Int → Int → Int → CharacterData
makeNewNonAddChar minStateIndex maxStateIndex charIndex =
    let bvMinState =
            if minStateIndex <= charIndex
                then BV.fromBits [True, False]
                else BV.fromBits [False, True]

        bvMaxState =
            if maxStateIndex <= charIndex
                then BV.fromBits [True, False]
                else BV.fromBits [False, True]

        bvState = bvMinState .|. bvMaxState
    in  emptyCharacter
            { stateBVPrelim = (V.singleton bvState, V.singleton bvState, V.singleton bvState)
            , stateBVFinal = V.singleton bvState
            }
