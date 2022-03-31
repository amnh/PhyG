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
        -- like a standard matrix with a single character
        singleCharVectList = V.toList $ fmap (U.getSingleCharacter charDataVV) (V.fromList [0.. numChars - 1])

        (recodedSingleVecList, newCharInfoLL) = unzip $ zipWith packNonAdd singleCharVectList (V.toList charInfoV)
        

        newTaxVectByCharVect = U.glueBackTaxChar (V.fromList $ concat recodedSingleVecList)

        
    in
    (nameBlock, newTaxVectByCharVect, V.fromList $ concat newCharInfoLL)

-- | packNonAdd takes taxon by vector cahracter data and list of character information
-- and returns bit packed and recoded non-additive characters and charInfo
-- input int is character index in block
packNonAdd :: V.Vector CharacterData -> CharInfo -> ([V.Vector CharacterData], [CharInfo])
packNonAdd inCharDataV charInfo =
    if charType charInfo /= NonAdd then ([inCharDataV],[charInfo])
    else 
        -- recode non-additive characters
        ([inCharDataV],[charInfo])

