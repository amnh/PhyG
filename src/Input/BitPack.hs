{- |
Module      :  BitPack.hs
Description :  Module with functionality to transform NonAdditive data to bit packed
               Word64 structures
Copyright   :  (c) 2022 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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


module Input.BitPack
  ( packNonAdditiveData
  , median2Packed
  ) where

import qualified Data.List                   as L
import qualified Data.List.Split             as SL
import           Data.Maybe
import           Types.Types
import qualified Data.BitVector.LittleEndian as BV
import qualified Data.Vector                 as V
import qualified Utilities.Utilities         as U
import           GeneralUtilities
import           Debug.Trace
import           Data.Word
import qualified ParallelUtilities as PU
import Control.Parallel.Strategies
import qualified Data.Text.Lazy            as T
import Data.Bits

{-
This module contains structures and functions for bit-packing operations
the basic idea is to transform the Bit-vector representation of non-additive
characters (whihc have effectively unlimited number of states) to more efficient 
operations based on Word64 encodings. Additionally, this should reduce
memory footprint when large number of non-additive characters are
input (such as genomic SNP data).

Several new types are cretged and manipulated based on the n umber of states
in the character--Packed2, Packed4, Packed5, Packed8, Packed64. These types hold
as subsets of bits (1 per state) multiple (other than for Packed64) original
non-additive characters.

The optimizations are based on using bi-wise operations specific
to each packed type to create preliminary (post-order with cost) and final
(pre-order) states.  Ideally the functions should contain no logic branches
or recursion (= loops) so that operations are solely bit-based.

Methods are similar too those of Lamport 1975, Ronquist 1998, Moelannen 1999, Goloboff 2002, 
and White and Holland 2011, but differ in detail due to
the programming language (Haskell) and need to maintain data structures 
useful for network analysis.

The basic character data structure is a vector of Word64, the original vector of bit-vectors
is split into a series of new characters depending on the number of states (in non-missing cells).

A single Word64 can then hold 32 2-state. 16 4-state, 12 5-state (awkward but useful for DNA),
4 8-state, 1 64-state, and a bit-vector for >64 states. 

Character weights are all = 1 in static charcaters.  This is ensured by organizeBlockData
in Input.Reorganize.hs basically--characters are multiplied by weight (if integer--otherwise not recoded)
So can check and only recode characters with weight of 1.
-}


{-
Functions for median2 calculations of packed types 
These are used in post-order graph traversals and pairwise
distance functions among others.
-}

-- | median2Packed takes two characters of packedNonAddTypes
-- and retuns new character data based on 2-median and cost
median2Packed :: CharType -> CharacterData -> CharacterData -> CharacterData
median2Packed inCharType leftChar rightChar =
    --assumes all weight 1
    let (newStateVect, newCost) = if inCharType == Packed2       then median2Word64 andOR2  (snd3 $ packedNonAddPrelim leftChar) (snd3 $ packedNonAddPrelim rightChar)
                                  else if inCharType == Packed4  then median2Word64 andOR4  (snd3 $ packedNonAddPrelim leftChar) (snd3 $ packedNonAddPrelim rightChar)
                                  --else if inCharType == Packed5  then median2Packed andOR5  (snd3 $ packedNonAddPrelim leftChar) (snd3 $ packedNonAddPrelim rightChar)
                                  --else if inCharType == Packed8  then median2Packed andOR8  (snd3 $ packedNonAddPrelim leftChar) (snd3 $ packedNonAddPrelim rightChar)
                                  else if inCharType == Packed64 then median2Word64 andOR64 (snd3 $ packedNonAddPrelim leftChar) (snd3 $ packedNonAddPrelim rightChar)
                                  else error ("Character type " ++ show inCharType ++ " unrecongized/not implemented")

        newCharacter = emptyCharacter { packedNonAddPrelim = (snd3 $ packedNonAddPrelim leftChar, newStateVect, snd3 $ packedNonAddPrelim rightChar)
                                      , localCost = fromIntegral newCost
                                      , globalCost = (fromIntegral newCost) + globalCost leftChar + globalCost rightChar
                                      }
    in
    newCharacter

{-
andOrN functions derived from White and Holland 2011
-}

-- | mask2A first Mask for 2 state 64 bit
-- 32 x (01)
-- 6148914691236517205
mask2A :: Word64
mask2A = 0x5555555555555555

-- | mask2B second Mask for 2 state 64 bit
-- 32 x (10)
-- 12297829382473034410
mask2B :: Word64
mask2B = 0xAAAAAAAAAAAAAAAA

-- | mask4A first mask for 4 state  64 bits
-- 8608480567731124087 
-- 16 X (0111)
mask4A :: Word64
mask4A = 0x7777777777777777

-- | mask4B second mask for 4 states 64 bits
-- 9838263505978427528 
-- 16 X (1000)
mask4B :: Word64
mask4B = 0x8888888888888888

-- | andOR2 and or function for Packed2 encoding
andOR2 :: Word64 -> Word64 -> (Word64, Int)
andOR2 x y = 
    let u = shiftR ((((x .&. y .&. mask2A) + mask2A) .|. (x .&. y)) .&. mask2B) 1 
        z = (x .&. y) .|. ((x .|. y) .&. ((u + mask2A) `xor` mask2B)) 
    in
    (z, popCount u)

-- | andOR4 and or function for Packed4 encoding
andOR4 :: Word64 -> Word64 -> (Word64, Int)
andOR4 x y = 
    let u = shiftR ((((x .&. y .&. mask4A) + mask4A) .|. (x .&. y)) .&. mask4B) 3 
        z = (x .&. y) .|. ((x .|. y) .&. ((u + mask4A) `xor` mask4B)) 
    in
    (z, popCount u)



-- | andOR64 and or function for Packed64 encoding
andOR64 :: Word64 -> Word64 -> (Word64, Int)
andOR64 x y =
     if  (x .&. y) /= zeroBits then (x .&. y, 0)
     else (x .|. y, popCount $ x .|. y)

-- | median2Word64 driver function for median of two PackedN states
median2Word64 :: (Word64 -> Word64 -> (Word64, Int)) -> V.Vector Word64 -> V.Vector Word64 -> (V.Vector Word64, Int)
median2Word64 andOrFun leftVect rightVect =
    let (stateVect, costVect) = V.unzip $ V.zipWith andOrFun leftVect rightVect
    in
    (stateVect, V.sum costVect)

{-
Functions for median3 calculations of packed types 
These are used in pre-order graph traversals and final state assignment
among others.
-}


{-
Functions to encode ("pack") non-additive characters into new Word64 characters
based on their number of states
-}

-- | packData takes input data and creates a variety of bit-packed data types
-- to increase efficiency and reduce footprint of non-additive characters
-- that are encoded as bitvectors
packNonAdditiveData :: ProcessedData -> ProcessedData
packNonAdditiveData (nameVect, bvNameVect, blockDataVect) = 
    let newBlockDataVect = fmap recodeNonAddCharacters blockDataVect
    in
    (nameVect, bvNameVect, newBlockDataVect)

-- | recodeNonAddCharacters takes block data, goes through characters
-- and recodes NonAdditive. 
-- Concat and list for charInfoV because new charcaters can be created
-- and newCharInfo then as well, could be multiple per input 'charcater'
recodeNonAddCharacters :: BlockData -> BlockData
recodeNonAddCharacters (nameBlock, charDataVV, charInfoV) =
    let numChars = V.length charInfoV
        
        -- create vector of single characters with vector of taxon data of sngle character each
        singleCharVectList = V.toList $ fmap (U.getSingleCharacter charDataVV) (V.fromList [0.. numChars - 1])

        -- bit pack the nonadd
        (recodedSingleVecList, newCharInfoLL) = unzip $ zipWith packNonAdd singleCharVectList (V.toList charInfoV)
        
        -- recreate BlockData, tacxon dominant structure
        newTaxVectByCharVect = U.glueBackTaxChar (V.fromList $ concat recodedSingleVecList)
        
    in
    (nameBlock, newTaxVectByCharVect, V.fromList $ concat newCharInfoLL)

-- | packNonAdd takes taxon by vector character data and list of character information
-- and returns bit packed and recoded non-additive characters and charInfo
-- input int is character index in block
-- the weight is skipping because of the weight replication in reorganize
-- if characters have non integer weight then they were not reorganized and left
-- as single BV--here as well. Should be very few (if any) of them.
packNonAdd :: V.Vector CharacterData -> CharInfo -> ([V.Vector CharacterData], [CharInfo])
packNonAdd inCharDataV charInfo =
    if (charType charInfo /= NonAdd) || (weight charInfo > 1)  then ([inCharDataV],[charInfo])
    else 
        -- recode non-additive characters
        let leafNonAddV =  fmap (snd3 . stateBVPrelim) inCharDataV
            numNonAdd = (V.length . V.head) leafNonAddV

            -- split characters into groups by states number 2,4,5,8,64, >64 (excluding missing)
            stateNumDataPairList = V.toList $ fmap (getStateNumber leafNonAddV) (V.fromList [0.. numNonAdd - 1]) 

            -- sort characters by states number (2, 4, 5, 8, 64, >64 -> 128)
            (state2CharL, state4CharL, state5CharL, state8CharL, state64CharL, state128CharL) = binStateNumber stateNumDataPairList ([],[],[],[],[],[])

            -- make new characters based on state size
            (newStateCharListList, newCharInfoList) = unzip $ (zipWith (makeStateNCharacter charInfo) [2,4,5,8,64,128] [state2CharL, state4CharL, state5CharL, state8CharL, state64CharL, state128CharL] `using` PU.myParListChunkRDS)

        in
        trace ("PNA: " ++ (show $ fmap fst stateNumDataPairList) ) --  ++ "\n" ++ (show (newStateCharListList, newCharInfoList) ))
        (fmap V.fromList newStateCharListList, concat newCharInfoList)

-- | makeStateNCharacter takes a list of characters each of which is a list of taxon character values and
-- creates a new character of all characters for give taxon and packs (64/ state number) characters into a 64 bit Word64
-- via chuncksOf--or if 64, not packing, if 128 stays bitvector
-- check for non-sequential states (A,T) or (0,2) etc
-- return is list of taxa x single new (packed) character
makeStateNCharacter ::  CharInfo -> Int -> [[BV.BitVector]] -> ([CharacterData], [CharInfo])
makeStateNCharacter charInfo stateNumber charDataLL = 
    let (recodeList, newCharInfo) = if stateNumber > 64 then recodeBV2BV charInfo charDataLL
                                    else if stateNumber == 64 then recodeBV2Word64Single charInfo charDataLL
                                    else recodeBV2Word64 charInfo stateNumber charDataLL
    in
    (recodeList, newCharInfo)

-- | recodeBV2BV take a list of BV.bitvector non-add characters and creates a list (taxa)
-- BV non-additive characters of type NonAdd.
-- this results in a single character and charInfo in list so can be concatenated
-- and removed if empty
recodeBV2BV :: CharInfo -> [[BV.BitVector]] -> ([CharacterData], [CharInfo])
recodeBV2BV charInfo charTaxBVLL =
    if null charTaxBVLL then ([],[])
    else 
        let -- convert to taxon by characgter data list
            newStateList = makeNewCharacterData charTaxBVLL

            -- rename with thype
            newCharName = T.append (name charInfo) $ T.pack "LargeState"

            -- create new characters for each taxon
            newCharDataList = fmap (makeNewData emptyCharacter) newStateList
        in
        (newCharDataList, [charInfo {name = newCharName, charType = NonAdd}])
        where makeNewData a b = a {stateBVPrelim = (b,b,b), stateBVFinal = b}


-- | recodeBV2Word64Single take a list of BV.bitvector non-add characters and creates a list (taxa)
-- of Word64 unpacked non-additive characters of type Packed64.
-- this results in a single character and charInfo in list so can be concatenated
-- and removed if empty
recodeBV2Word64Single :: CharInfo -> [[BV.BitVector]] -> ([CharacterData], [CharInfo])
recodeBV2Word64Single charInfo charTaxBVLL =
    if null charTaxBVLL then ([],[])
    else 
        let newCharName = T.append (name charInfo) $ T.pack "64State"
        
            -- convert BV to Word64
            taxWord64BLL = fmap (fmap BV.toUnsignedNumber) charTaxBVLL

            -- convert to taxon by character data lisyt
            newStateList = makeNewCharacterData taxWord64BLL

            -- make new character data
            newCharDataList = fmap (makeNewData emptyCharacter) newStateList
        in
        (newCharDataList, [charInfo {name = newCharName, charType = Packed64}])
        where makeNewData a b = a {packedNonAddPrelim = (b,b,b), packedNonAddFinal = b}

-- | makeNewCharacterData takes a list of characters, each of which is a list of taxon states
-- of type a (bitvector or Word64) and returns a list of taxa each of which is a vector
-- of type a charactyer data
makeNewCharacterData :: [[a]] -> [V.Vector a]
makeNewCharacterData charByTaxSingleCharData  =
    let taxonByCharL = L.transpose charByTaxSingleCharData
        taxonByCharV = fmap V.fromList taxonByCharL
    in
    taxonByCharV


-- | recodeBV2Word64 take a list of BV.bitvector non-add characters and the states number of creates
-- Word64 representaions where subcharcaters are created and shifted to proper positions and ORd
-- to create packed reresentation--new character types Packed2, Packed4, Packed5, and Packed8. 
-- this results in a single character and charInfo in list so can be concatenated
-- and removed if empty
recodeBV2Word64 :: CharInfo -> Int -> [[BV.BitVector]] -> ([CharacterData], [CharInfo])
recodeBV2Word64 charInfo stateNumber charTaxBVLL =
    trace ("Enter RBV2W64: " ++ (show stateNumber) ++ " " ++ (show (length charTaxBVLL, fmap length charTaxBVLL))) (
    if null charTaxBVLL then ([],[])
    else
        let newCharType = if stateNumber == 2 then Packed2
                          else if stateNumber == 4 then Packed4
                          else if stateNumber == 5 then Packed5
                          else if stateNumber == 8 then Packed8
                          else error ("State number " ++ (show stateNumber) ++ " not to be packed in recodeBV2Word64")

            newCharName = T.append (name charInfo) $ T.pack ((show stateNumber) ++ "State")
            
            -- get number of characters that can be packed into Word64 for that state number
            numCanPack = fst $ divMod 64 stateNumber

            -- convert to taxon by character data list
            taxCharBVLL = L.transpose charTaxBVLL

            -- get state index list for all characters (could be non sequential 0|2; A|T etc)
            stateIndexLL = fmap getStateIndexList charTaxBVLL 

            -- convert data each taxon into packedWord64
            packedDataL = fmap (packIntoWord64 stateNumber numCanPack stateIndexLL) taxCharBVLL

        in
        -- trace ("RBV2W64: " ++ (show packedDataL))
        (packedDataL, [charInfo {name = newCharName, charType = newCharType}])
        )


-- | packIntoWord64 takes a list of bitvectors for a taxon, the state number and number that can be packed into 
-- a Word64 and performs appropriate bit settting and shifting to create  Word64
packIntoWord64 :: Int -> Int -> [[Int]] -> [BV.BitVector] -> CharacterData
packIntoWord64 stateNumber numToPack stateCharacterIndexL inBVList =
    -- get packable chunk of bv and correcsponding state indices
    let packBVList = SL.chunksOf numToPack inBVList
        packIndexLL = SL.chunksOf numToPack stateCharacterIndexL

        -- pack each chunk 
        packedWordVect = V.fromList $ zipWith (makeWord64FromChunk stateNumber) packIndexLL packBVList

    in    
    emptyCharacter { packedNonAddPrelim = (packedWordVect, packedWordVect, packedWordVect)
                   , packedNonAddFinal = packedWordVect
                   }


-- | makeWord64FromChunk takes a list (= chunk) of bitvectors and cretes bit subcharacter (Word64)
-- with adjacent bits for each BV in chunk.  It then bit shifts the appropriate number of bits for each member
-- of the chunk and finally ORs all (64/stateNumber) Word64s to  make the final packed representation
makeWord64FromChunk ::  Int -> [[Int]] ->  [BV.BitVector] -> Word64
makeWord64FromChunk stateNumber stateIndexLL bvList =
    if null bvList then (0 :: Word64)
    else
        let subCharacterList = zipWith3 (makeSubCharacter stateNumber) stateIndexLL bvList [0..(length bvList - 1)]
        in
        trace ("MW64FC: " ++ (show subCharacterList) ++ " " ++ (show $ L.foldl1' (.|.) subCharacterList))
        L.foldl1' (.|.) subCharacterList

-- | makeSubCharacter makes sub-character (ie only those bits for states) from single bitvector and shifts appropriate number of bits
-- to make Word64 with sub character bits set and all other bits OFF and in correct bit positions for that sub-character
makeSubCharacter :: Int -> [Int] -> BV.BitVector -> Int -> Word64
makeSubCharacter stateNumber stateIndexList inBV subCharacterIndex =
    trace ("Making sub character:" ++ (show stateNumber ++ " " ++ (show stateIndexList) ++ " " ++ (show subCharacterIndex) ++ (show inBV))) (
    let -- get bit of state indices
        bitStates = fmap (testBit inBV) stateIndexList

        -- get index of states when only minimally bit encoded (0101, 0001 -> 11, 01)
        newBitStates = setOnBits (zeroBits :: Word64) bitStates 0 
        subCharacter = shiftL newBitStates (subCharacterIndex * stateNumber)
    in
    trace ("MSC: " ++ (show subCharacterIndex) ++ " " ++ (show bitStates) ++ " " ++ (show newBitStates) ++ " " ++ (show subCharacter)) (
    -- cna remove this check when working
    if length stateIndexList `notElem` [((fst $ divMod 2 stateNumber) + 1) .. stateNumber] then error ("State number of index list do not match: " ++ (show (stateNumber, length stateIndexList, stateIndexList)))
    else 
        subCharacter
    )
    )
        
-- | setOnBits recursively sets On bits in a list of Bool
setOnBits :: Word64 -> [Bool] -> Int -> Word64
setOnBits baseVal onList bitIndex =
    if null onList then baseVal
    else
        let newVal = if (head onList == True) then setBit baseVal bitIndex
                     else baseVal
        in
        setOnBits newVal (tail onList) (bitIndex + 1)

-- | getStateIndexList takes list of list bit vectors and for each taxon for a given bv character 
-- and returns a list of 
-- bit indices of states in the bv this because states can be non-seqeuntial (0|3)
-- used to have a list of all states used (ON) in a character in all taxa
getStateIndexList :: [BV.BitVector] -> [Int]
getStateIndexList taxBVL = 
    if null taxBVL then []
    else 
        let inBV = L.foldl1' (.|.) taxBVL
            onList = fmap (testBit inBV) [0.. (finiteBitSize inBV) - 1]
            onIndexPair = zip onList [0.. (finiteBitSize inBV) - 1]
            indexList = fmap snd $ filter ((== True) .fst) onIndexPair

        in
        trace ("GSIL: " ++ (show indexList))
        indexList 


-- | binStateNumber takes a list of pairs of char states number and data column as list of bitvectors and 
-- into list for 2,4,5,8,64,>64
binStateNumber :: [(Int, [BV.BitVector])] 
               -> ([[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]]) 
               -> ([[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]])
binStateNumber inPairList (cur2, cur4, cur5, cur8, cur64, cur128) =
    if null inPairList then 
        --dont' really need to reverse here but seems hygenic
        trace ("BSN: " ++ (show (length cur2, length cur4, length cur5, length cur8, length cur64,  length cur128)))
        (L.reverse cur2, L.reverse cur4, L.reverse cur5, L.reverse cur8, L.reverse cur64,  L.reverse cur128)
    else 
        let (stateNum, stateData) = head inPairList
        in
        -- skip--constant
        if stateNum < 2 then binStateNumber (tail inPairList) (cur2, cur4, cur5, cur8, cur64, cur128)
        else if stateNum == 2 then binStateNumber (tail inPairList) (stateData : cur2, cur4, cur5, cur8, cur64, cur128)
        else if stateNum <= 4 then binStateNumber (tail inPairList) (cur2, stateData : cur4, cur5, cur8, cur64, cur128)
        else if stateNum <= 5 then binStateNumber (tail inPairList) (cur2, cur4, stateData : cur5, cur8, cur64, cur128)
        else if stateNum <= 8 then binStateNumber (tail inPairList) (cur2, cur4, cur5, stateData : cur8, cur64, cur128)
        else if stateNum <= 64 then binStateNumber (tail inPairList) (cur2, cur4, cur5, cur8, stateData : cur64, cur128)
        else binStateNumber (tail inPairList) (cur2, cur4, cur5, cur8, cur64, stateData : cur128)

-- | getStateNumber returns the number of uniqe (non missing) states for a 'column' 
-- of nonadd bitvector values
-- the charState values are in ranges for 2,4,5,8,64 and bigger numbers
-- missingVal not takeen from alphabet size since that was not updated in reorganize.
-- So take OR of all on bits--may be non-sequential--ie 0 2 7 so need to watch that.
-- returns pair of stateNUmber class (2,4,5,8,64, >64 as 128) and list of states
-- for efficient glueing back together later

getStateNumber :: V.Vector (V.Vector BV.BitVector) -> Int -> (Int, [BV.BitVector])
getStateNumber  characterDataVV characterIndex =
    let thisCharV = V.map (V.! characterIndex) characterDataVV
        missingVal = V.foldl1' (.|.) thisCharV
        nonMissingBV = V.foldl1' (.|.) $ V.filter (/= missingVal) thisCharV
        numStates = popCount nonMissingBV

        -- this turns off non-missing bits
        thisCharL = V.toList $ (fmap (.&. nonMissingBV) thisCharV)
    in
    if numStates <= 2 then (2, thisCharL)
    else if numStates <= 4 then (4, thisCharL)
    else if numStates <= 5 then (5, thisCharL)
    else if numStates <= 8 then (8, thisCharL)
    else if numStates <= 64 then (64, thisCharL)
    else (128, thisCharL)

