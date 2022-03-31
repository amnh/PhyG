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

-- | packNonAdd takes taxon by vector cahracter data and list of character information
-- and returns bit packed and recoded non-additive characters and charInfo
-- input int is character index in block
-- the weight is skipping becasue of the weight replication in reorganize
-- if chracters have non integer weight then they were not reorganized and left
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
        in
        trace ("PNA: " ++ (show $ fmap fst stateNumDataPairList))
        ([inCharDataV],[charInfo])

-- | makeStateNCharacter takes a list of charcaters each of which is a list of taxon caracter values and
-- creates a new character of all characters for give taxon and packs (64/ state number) characters into a 64 bit Word64
-- via chuncksOf--or if 64, not packing, if 128 stays bitvector
-- check for non-sequential states (A,T) or (0,2) etc
makeStateNCharacter ::  CharInfo -> Int -> [[BV.BitVector]] -> CharacterData
makeStateNCharacter charInfo stateNumber charDataLL = 
    emptyCharacter 

-- | binStateNumber takes a list of pairs ofg char states number and data column as list of bitvectors and 
-- into list for 2,4,5,8,64,>64
binStateNumber :: [(Int, [BV.BitVector])] 
               -> ([[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]]) 
               -> ([[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]])
binStateNumber inPairList (cur2, cur4, cur5, cur8, cur64, cur128) =
    if null inPairList then 
        --d ont' really need to reverse here but seems hygenic
        (L.reverse cur2, L.reverse cur4, L.reverse cur5, L.reverse cur8, L.reverse cur64,  L.reverse cur128)
    else 
        let (stateNum, stateData) = head inPairList
        in
        -- skip--uninformative
        if stateNum == 1 then binStateNumber (tail inPairList) (cur2, cur4, cur5, cur8, cur64, cur128)
        else if stateNum <= 2 then binStateNumber (tail inPairList) (stateData : cur2, cur4, cur5, cur8, cur64, cur128)
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
        numStates = popCount $ V.foldl1' (.|.) $ V.filter (/= missingVal) thisCharV
        thisCharL = V.toList thisCharV
    in
    if numStates <= 2 then (2, thisCharL)
    else if numStates <= 4 then (4, thisCharL)
    else if numStates <= 5 then (5, thisCharL)
    else if numStates <= 8 then (8, thisCharL)
    else if numStates <= 64 then (64, thisCharL)
    else (128, thisCharL)