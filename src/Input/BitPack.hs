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
-}

{-
The basic character data structure is a vector of Word64, the original vector of bit-vectors
is split into a series of new characters depending on the number of states (in non-missing cells).

A single Word64 can then hold 32 2-state. 16 4-state, 12 5-state (awkward but useful for DNA),
4 8-state, 1 64-state, and a bit-vector for >64 states.

The AND/OR operations for post-order 2-medians occur by total word AND/OR
then through bit shifting and masking examine each sub-character (bit range) in turn.

in brief, each state in each sub-character is bit shifted right to the 0th position (big-endian)
AND-ed with '1' (ie all 0 nbit except the 0th) these state values are OR-ed over the states.
If there was overlap between the two inputs (ie no change) this will be 1, else 0.
The complement of this value is taken and added to that for all the sub-characters yielding
the unweighted cost of the pair of inputs. This value (call it overlapBit for now 1 for yes, 0 no)
is used for median calculation

The median states are calculated by keeping the AND and OR sub character states (all other masked to zero
via an AND operation) and mulitplying the AND substates by the overlapBit and
the OR states by the complement of the overlap bit and ORing the two results

THis could no doubt be made more efficient, but at least this has no logic
although worth testing the multiplicatoin versus a function (hence logic):
f overlapBit a
where f 0 a = 0 
f _ a = a.

so (F overlapBit subcharater AND) OR (F ~overlapBit subcharactr OR)
-} 

{-
Character weights are all = 1 in static charcaters.  This is ensured by organizeBlockData
in Input.Reorganize.hs basically--charcters are nultiplied by weight (if integer--otherwise not recoded)
So can check and only recode characters with weight of 1.
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
            numNonAdd = V.length $ V.head leafNonAddV

            -- split characters into groups by states number 2,4,5,8,64, >64 (excluding missing)
            stateNumDataPairList = V.toList $ fmap (getStateNumber leafNonAddV) (V.fromList [0.. numNonAdd - 1]) 

            -- sort characters by states number (2, 4, 5, 8, 64, >64 -> 128)
            (state2CharL, state4CharL, state5CharL, state8CharL, state64CharL, state128CharL) = binStateNumber stateNumDataPairList ([],[],[],[],[],[])

            -- make new characters based on state size
            (newStateCharListList, newCharInfoList) = unzip $ (zipWith (makeStateNCharacter charInfo) [2,4,5,8,64,128] [state2CharL, state4CharL, state5CharL, state8CharL, state64CharL, state128CharL] `using` PU.myParListChunkRDS)

        in
        trace ("PNA: " ++ (show $ fmap fst stateNumDataPairList) ++ "\n" ++ (show(newStateCharListList, newCharInfoList) ))
        (fmap V.fromList newStateCharListList, concat newCharInfoList)

-- | makeStateNCharacter takes a list of charcaters each of which is a list of taxon caracter values and
-- creates a new character of all characters for give taxon and packs (64/ state number) characters into a 64 bit Word64
-- via chuncksOf--or if 64, not packing, if 128 stays bitvector
-- check for non-sequential states (A,T) or (0,2) etc
-- return is list of taxa x single new character
makeStateNCharacter ::  CharInfo -> Int -> [[BV.BitVector]] -> ([CharacterData], [CharInfo])
makeStateNCharacter charInfo stateNumber charDataLL = 
    let (recodeList, newCharInfo) = if stateNumber > 64 then recodeBV2BV charInfo charDataLL
                                    else if stateNumber == 64 then recodeBV2Word64Single charInfo charDataLL
                                    else recodeBV2Word64 charInfo stateNumber charDataLL
    in
    (recodeList, newCharInfo)

-- | recodeBV2BV take a list of BV.bitvector non-add characters and creates a list (taxa)
-- BV non-additive characters of type NonAdd.
recodeBV2BV :: CharInfo -> [[BV.BitVector]] -> ([CharacterData], [CharInfo])
recodeBV2BV charInfo taxBVLL =
    if null taxBVLL then ([],[])
    else 
        let newStateList = fmap V.fromList taxBVLL
            newCharName = T.append (name charInfo) $ T.pack "LargeState"
            newCharDataList = fmap (makeNewData emptyCharacter) newStateList
        in
        (newCharDataList, [charInfo {name = newCharName, charType = NonAdd}])
        where makeNewData a b = a {stateBVPrelim = (b,b,b), stateBVFinal = b}


-- | recodeBV2Word64Single take a list of BV.bitvector non-add characters and creates a list (taxa)
-- of Word64 unpacked non-additive characters of type Packed64.
recodeBV2Word64Single :: CharInfo -> [[BV.BitVector]] -> ([CharacterData], [CharInfo])
recodeBV2Word64Single charInfo taxBVLL =
    if null taxBVLL then ([],[])
    else 
        let newCharName = T.append (name charInfo) $ T.pack "64State"
        
            -- convert BV to Word64
            taxWord64BLL = fmap (fmap BV.toUnsignedNumber) taxBVLL

            -- convert to vector Word64
            newStateList = fmap V.fromList taxWord64BLL

            -- make new character data
            newCharDataList = fmap (makeNewData emptyCharacter) newStateList
        in
        (newCharDataList, [charInfo {name = newCharName, charType = Packed64}])
        where makeNewData a b = a {packedNonAddPrelim = (b,b,b), packedNonAddFinal = b}


-- | recodeBV2Word64 take a list of BV.bitvector non-add characters and the states number of creates
-- Word64 representaions where subcharcaters are created and shifted to proper positions and ORd
-- to create packed reresentation--new character types Packed2, Packed4, Packed5, and Packed8. 
recodeBV2Word64 :: CharInfo -> Int -> [[BV.BitVector]] -> ([CharacterData], [CharInfo])
recodeBV2Word64 charInfo stateNumber taxBVLL =
    trace ("Enter RBV2W64: " ++ (show stateNumber) ++ " " ++ (show (length taxBVLL, fmap length taxBVLL))) (
    if null taxBVLL then ([],[])
    else
        let newCharType = if stateNumber == 2 then Packed2
                          else if stateNumber == 4 then Packed4
                          else if stateNumber == 5 then Packed5
                          else if stateNumber == 8 then Packed8
                          else error ("State number " ++ (show stateNumber) ++ " not to be packed in recodeBV2Word64")

            newCharName = T.append (name charInfo) $ T.pack ((show stateNumber) ++ "State")
            
            -- get number of characters that can be packed into Word64 for that state number
            numCanPack = fst $ divMod 64 stateNumber

            -- get state index list for all characters (could be non seqeuenctial 0|2; A|T etc)
            stateIndexLL = getStateIndexList taxBVLL 
            chunkStateIndexLLL = SL.chunksOf numCanPack stateIndexLL

            -- create chunks of charcaters to be put into single element of vector of Word64
            chunkDataListL = SL.chunksOf numCanPack taxBVLL

            -- pack chunks from each taxon into Word64
            packedDataL = fmap (chunksToWord64 stateNumber chunkStateIndexLLL) chunkDataListL

        in
        trace ("RBV2W64: " ++ (show packedDataL))
        (packedDataL, [charInfo {name = newCharName, charType = newCharType}])
        )


-- | chunksToWord64 take states number and chunk (list of BV.bitvector non-additive)
-- and a list of the indices of the character states  (could be non seqeuential 0|2; A|T etc)
-- and converts to a single Word64 by setting bits and shifting
-- yeilds a single taxon's new character
-- since chunks of both char indices and bitvecots--that's why list of list of list
chunksToWord64 :: Int -> [[[Int]]] ->  [[BV.BitVector]] -> CharacterData
chunksToWord64 stateNumber stateCharIndexLLL taxChunk =
    let word64Vect = V.fromList $ zipWith (makeWord64FromChunk stateNumber) stateCharIndexLLL taxChunk 
        word64Prelim = (word64Vect, word64Vect, word64Vect)
    in
    trace ("C2W64: " ++ (show word64Vect))
    emptyCharacter {packedNonAddPrelim = word64Prelim, packedNonAddFinal = word64Vect }

-- | makeWord64FromChunk takes a list (= chunk) of bitvectors and cretes bit subcharacter (Word64)
-- with adjacent bits for each BV in chunk.  It then bit shifts the appropriate number of bits for each member
-- of the chunk and finally ORs all (64/stateNumber) Word64s to  make the final packed representation
makeWord64FromChunk ::  Int -> [[Int]] ->  [BV.BitVector] -> Word64
makeWord64FromChunk stateNumber chunkIndexLL bvChunkList =
    if null bvChunkList then (0 :: Word64)
    else
        let subCharacterList = zipWith3 (makeSubCharacter stateNumber) chunkIndexLL bvChunkList [0..(length bvChunkList - 1)]
        in
        trace ("MW64FC: " ++ (show subCharacterList) ++ " " ++ (show $ L.foldl1' (.|.) subCharacterList))
        L.foldl1' (.|.) subCharacterList

-- | makeSubCharacter makes subcharcter from single bitvector and shifts appropriate number of bits
-- to make Word64 with sub character bits set and all other bits OFF 
makeSubCharacter :: Int -> [Int] -> BV.BitVector -> Int -> Word64
makeSubCharacter stateNumber stateIndexList inBV chunkIndex =
    trace ("Making sub character:" ++ (show stateNumber ++ " " ++ (show stateIndexList) ++ " " ++ (show chunkIndex) ++ (show inBV))) (
    let -- get bit of state indices
        bitStates = fmap (testBit inBV) stateIndexList

        -- get index of states when only minimally bit encoded (0101, 0001 -> 11, 01)
        newBitStates = setOnBits (zeroBits :: Word64) bitStates 0 
        subCharacter = shiftL newBitStates (chunkIndex * stateNumber)
    in
    trace ("MSC: " ++ (show bitStates) ++ " " ++ (show newBitStates) ++ " " ++ (show subCharacter)) (
    -- cna remove this check when working
    if stateNumber /= length stateIndexList then error ("State number of index list do not match: " ++ (show (stateNumber, length stateIndexList, stateIndexList)))
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

-- | getStateIndexList takes lits of list bit vectors and for each returns a list of 
-- bit indices of states in the bv this becasue starets can be non-seqeuntial (0|3)
getStateIndexList :: [[BV.BitVector]] -> [[Int]]
getStateIndexList taxBVLL = 
    if null taxBVLL then []
    else 
        let firstCharBV = head taxBVLL
            stateIndexL = getStateIndices firstCharBV
        in
        trace ("GSIL: " ++ (show firstCharBV) ++ " " ++ (show stateIndexL))
        stateIndexL : getStateIndexList (tail taxBVLL)

-- | getStateIndices takes list of vectors of BVs and index and gets the states for that index bv
getStateIndices :: [BV.BitVector] -> [Int]
getStateIndices bvList =
    let -- get all ON bits
        onBV = L.foldl1' (.|.) bvList

        -- test for ON bits
        onList = fmap (testBit onBV) [0.. (finiteBitSize onBV) - 1]
        onIndexPair = zip onList [0.. (finiteBitSize onBV) - 1]
        indexList = fmap snd $ filter ((== True) .fst) onIndexPair
    in
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