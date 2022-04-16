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
  , packedPreorder
  , threeWayPacked
  , unionPacked
  , maxCharDiff
  ) where

import           Control.Parallel.Strategies
import           Data.Bits
import qualified Data.BitVector.LittleEndian as BV
import           Data.Char                   (intToDigit)
import qualified Data.List                   as L
import qualified Data.List.Split             as SL
import qualified Data.Text.Lazy              as T
import qualified Data.Vector                 as V
import           Data.Word
import           Debug.Trace
import           GeneralUtilities
import           Numeric                     (showIntAtBase)
import qualified ParallelUtilities           as PU
import           Types.Types
import qualified Utilities.Utilities         as U


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

-- | showBits cponverts Value to bits as Srting
showBits :: Word64 -> String
showBits inVal = showIntAtBase 2 intToDigit inVal ""

-- | showBitsV shoiw vector of bits
showBitsV :: V.Vector Word64 -> String
showBitsV inValV = concat $ fmap (++ " ") $ V.toList $ fmap showBits inValV


-- | maxCharDiff get the approximate maximum differnet in number of states
-- could do exact with masking, but this likely good enough for general purposes
maxCharDiff :: CharType -> Word64 -> Word64 -> Int
maxCharDiff inCharType a b =
    let numDiffBits = popCount $ xor a b
        numPacked = if inCharType == Packed2       then 2
                    else if inCharType == Packed4  then 4
                    else if inCharType == Packed5  then 5
                    else if inCharType == Packed8  then 8
                    else if inCharType == Packed64 then 64
                    else error ("Character type " ++ show inCharType ++ " unrecognized/not implemented")
    in
    if inCharType == Packed64 then if a == b then 0 else 1
    else
        let (maxNum, _) = divMod numDiffBits numPacked
        in
        maxNum

-- | median2Packed takes two characters of packedNonAddTypes
-- and retuns new character data based on 2-median and cost
median2Packed :: CharType -> CharacterData -> CharacterData -> CharacterData
median2Packed inCharType leftChar rightChar =
    --assumes all weight 1
    let (newStateVect, newCost) = if inCharType == Packed2       then median2Word64 andOR2  (snd3 $ packedNonAddPrelim leftChar) (snd3 $ packedNonAddPrelim rightChar)
                                  else if inCharType == Packed4  then median2Word64 andOR4  (snd3 $ packedNonAddPrelim leftChar) (snd3 $ packedNonAddPrelim rightChar)
                                  else if inCharType == Packed5  then median2Word64 andOR5  (snd3 $ packedNonAddPrelim leftChar) (snd3 $ packedNonAddPrelim rightChar)
                                  else if inCharType == Packed8  then median2Word64 andOR8  (snd3 $ packedNonAddPrelim leftChar) (snd3 $ packedNonAddPrelim rightChar)
                                  else if inCharType == Packed64 then median2Word64 andOR64 (snd3 $ packedNonAddPrelim leftChar) (snd3 $ packedNonAddPrelim rightChar)
                                  else error ("Character type " ++ show inCharType ++ " unrecognized/not implemented")

        newCharacter = emptyCharacter { packedNonAddPrelim = (snd3 $ packedNonAddPrelim leftChar, newStateVect, snd3 $ packedNonAddPrelim rightChar)
                                      , localCost = fromIntegral newCost
                                      , globalCost = (fromIntegral newCost) + globalCost leftChar + globalCost rightChar
                                      }
    in
    -- trace ("M2P: " ++ (showBitsV $ (snd3 . packedNonAddPrelim) leftChar) ++ " " ++ (showBitsV $ (snd3 . packedNonAddPrelim) rightChar) ++ " -> " ++   (showBitsV $ (snd3 . packedNonAddPrelim) newCharacter) ++ " at cost " ++ (show newCost))
    newCharacter

-- | unionPacked returns character that is the union (== OR) for bit packed characters
-- of the final fields as preliminary and final
unionPacked :: CharacterData -> CharacterData -> CharacterData
unionPacked charL charR =
    let newVect = V.zipWith (.|.) (packedNonAddFinal charL) (packedNonAddFinal charR)
    in
    emptyCharacter { packedNonAddPrelim = (newVect, newVect, newVect)
                   , packedNonAddFinal = newVect
                   }

{-
Masks for vaious operations and state numbers
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

-- | mask4C mak for 4 states an 64 bits
-- 4919131752989213764
-- 16 x (0100)
mask4C :: Word64
mask4C = 0x4444444444444444

-- | mask4D mak for 4 states an 64 bits
-- 2459565876494606882
-- 16 x (0010)
mask4D :: Word64
mask4D = 0x2222222222222222

-- | mask4E mak for 4 states an 64 bits
-- 1229782938247303441
-- 16 x (0001)
mask4E :: Word64
mask4E = 0X1111111111111111

{-
5 state masks top 4 bits OFF may require to mask out top 4 bits for states and cost
top 4 OFF rest ON 0xFFFFFFFFFFFFFFF
                  1152921504606846975
top 4 ON rest OFF 0xF000000000000000
                  17293822569102704640
-}

-- | mask5A first mask for 5 states 64 bits -- need to check what state of top 4  bits--should be OFF I think
-- 12 x (01111)
-- 557865244164603375
mask5A :: Word64
mask5A = 0x7BDEF7BDEF7BDEF

-- | mask5B second mask for 5 states 64 bits -- need to check what state of top 4  bits-these are OFF
-- 12 x (10000)
-- 595056260442243600 (top 4 OFF) v 17888878829544948240 (top 4 ON)
-- 0x842108421084210 (top 4 OFF) v  F842108421084210 (top 4 ON)
mask5B :: Word64
mask5B = 0x842108421084210
-- | mask5C mask 5 states 64 bits
-- 12 x (01000)
-- 297528130221121800
mask5C :: Word64
mask5C = 0x421084210842108

-- | mask5D mask 5 states 64 bits
-- 12 x (00100)
-- 148764065110560900
mask5D :: Word64
mask5D = 0x210842108421084

-- | mask5E mask 5 states 64 bits
-- 12 x (00010)
-- 74382032555280450
mask5E :: Word64
mask5E = 0x108421084210842

-- | mask5F mask 5 states 64 bits
-- 12 x (00001)
-- 37191016277640225
mask5F :: Word64
mask5F = 0x84210842108421

-- | mask8A first mask for 8 states 64 bits
-- 8 x (01111111)
-- 9187201950435737471
mask8A :: Word64
mask8A = 0x7F7F7F7F7F7F7F7F

-- | mask8B second mask for 8 states 64 bits
-- 8 x (10000000)
-- 9259542123273814144
mask8B :: Word64
mask8B = 0x8080808080808080

-- | mask8C mask for 8 states 64 bits
-- 8 x (01000000)
-- 4629771061636907072
mask8C :: Word64
mask8C = 0x4040404040404040

-- | mask8D mask for 8 states 64 bits
-- 8 x (00100000)
-- 2314885530818453536
mask8D :: Word64
mask8D = 0x2020202020202020

-- | mask8E mask for 8 states 64 bits
-- 8 x (00010000)
-- 1157442765409226768
mask8E :: Word64
mask8E = 0x1010101010101010

-- | mask8F mask for 8 states 64 bits
-- 8 x (00001000)
-- 578721382704613384
mask8F :: Word64
mask8F = 0x808080808080808

-- | mask8G mask for 8 states 64 bits
-- 8 x (00000100)
-- 289360691352306692
mask8G :: Word64
mask8G = 0x404040404040404

-- | mask8H mask for 8 states 64 bits
-- 8 x (00000010)
-- 144680345676153346
mask8H :: Word64
mask8H = 0x202020202020202

-- | mask8I mask for 8 states 64 bits
-- 8 x (00000001)
-- 72340172838076673
mask8I :: Word64
mask8I = 0x101010101010101

{-
andOrN functions derived from White and Holland 2011
-}

{-
-- | andOR2 Packed2 modified from Goloboff 2002
-- this is incomplete-- Goloboff uses look up table for on bits to length
andOR2 :: Word64 -> Word64 -> (Word64, Int)
andOR2 x y =
   let  x1 = x .&. y
        c1 = xor mask2B (shiftR ((x1 .&. mask2B) .|. (x1 .&. mask2A)) 1)
        c2 = c1 .|. (shiftL c1 1)
        newState = x1 .|. c2
        numChanges = lookUpLegnth((c1 .|. (shiftR c1 31)) .&. 0xFFFFFFFF)
   in
   (newState, numChanges)
   -}


-- For all biut packed charters--post order median 2
-- for now--for somwe reason either mis-diagnosed or mis-coded (from White and Holland 2011)
-- the "on" bit in u is reflective of number of intersections not unions.
-- hence subtracting the number of unions from numbers of characters
-- determined by leading OFF bits since packing will likely have ragged edges
-- no not always the pack-able number

-- | andOR2 and or function for Packed2 encoding
andOR2 :: Word64 -> Word64 -> (Word64, Int)
andOR2 x y =
    let u = shiftR ((((x .&. y .&. mask2A) + mask2A) .|. (x .&. y)) .&. mask2B) 1
        z = (x .&. y) .|. ((x .|. y) .&. ((u + mask2A) `xor` mask2B))

        -- get number of characters by checking states (may not be full)
        numEmptyBits = countLeadingZeros x --- could be y just as well

        -- shift divide by 2 states
        numNonCharacters = shiftR numEmptyBits 1
        numChars =  32 - numNonCharacters
    in
    {-
    trace ("AO2 numChars:" ++ (show numChars) ++ " x & y:" ++ (showBits $ x .&. y) ++ "\nx .&. y .&. mask2A:" ++ (showBits $ (x .&. y .&. mask2A)) ++ "\n((x .&. y .&. mask2A) + mask2A):" ++ (showBits $ ((x .&. y .&. mask2A) + mask2A))
      ++ "\n:(((x .&. y .&. mask2A) + mask2A) .|. (x .&. y)): " ++ (showBits $ (((x .&. y .&. mask2A) + mask2A) .|. (x .&. y))) ++ "\n:((((x .&. y .&. mask2A) + mask2A) .|. (x .&. y)) .&. mask2B):"
      ++ (showBits $ ((((x .&. y .&. mask2A) + mask2A) .|. (x .&. y)) .&. mask2B)) ++ "\nu: " ++ (showBits u)
      ++"\npc: " ++ (show $ popCount u) ++ " x:" ++ (showBits x) ++ " y:" ++ (showBits y) ++ " => u:" ++ (showBits u) ++ " z:" ++ (showBits z)) -- ++ " mask2A:" ++ (showBits mask2A) ++ " mask2B:" ++ (showBits mask2B))
    -}
    (z, numChars - (popCount u))

-- | andOR4 and or function for Packed4 encoding
andOR4 :: Word64 -> Word64 -> (Word64, Int)
andOR4 x y =
    let u = shiftR ((((x .&. y .&. mask4A) + mask4A) .|. (x .&. y)) .&. mask4B) 3
        z = (x .&. y) .|. ((x .|. y) .&. ((u + mask4A) `xor` mask4B))

        -- get number of characters by checking states (may not be full)
        numEmptyBits = countLeadingZeros x --- could be y just as well

        -- shift divide by 4 states
        numNonCharacters = shiftR numEmptyBits 2
        numChars =  16 - numNonCharacters
    in
    (z, numChars - (popCount u))

-- | andOR5 and or function for Packed5 encoding
-- potential issue with top 4 bits--not sure on mask5B whether top 4 should be on or OFF.
-- can always maske top 4 with AND 0000111... (0xFFFFFFFFFFFFFFF or 1152921504606846975)
-- to remove bits for counting
-- and calcualted state
andOR5:: Word64 -> Word64 -> (Word64, Int)
andOR5 x y =
    let u = shiftR ((((x .&. y .&. mask5A) + mask5A) .|. (x .&. y)) .&. mask5B) 4
        z = (x .&. y) .|. ((x .|. y) .&. ((u + mask5A) `xor` mask5B))

        -- get number of characters by checking states (may not be full)
        numEmptyBits = countLeadingZeros x --- could be y just as well

        -- since top 4 bits always off (can't put anyhting in there) need to subtract those zeros
        -- to get character number. Cant shift to get / 5 so integer divide
        (numNonCharacters, _) =  divMod (numEmptyBits - 4) 5
        numChars =  12 - numNonCharacters
    in
    -- trace ("AO5 numChars:" ++ (show numChars) ++ " x & y:" ++ (showBits $ x .&. y) ++ " u:" ++ (showBits u) ++ " z:" ++ (showBits z) ++ " leading 0:" ++ (show numEmptyBits) ++ " non-chars:" ++ (show numNonCharacters) ++ " popCount u:" ++ (show $ popCount u))
    (z, numChars - (popCount u))

-- | andOR8 and or function for Packed8 encoding
andOR8 :: Word64 -> Word64 -> (Word64, Int)
andOR8 x y =
    let u = shiftR ((((x .&. y .&. mask8A) + mask8A) .|. (x .&. y)) .&. mask8B) 7
        z = (x .&. y) .|. ((x .|. y) .&. ((u + mask8A) `xor` mask8B))

        -- get number of characters by checking states (may not be full)
        numEmptyBits = countLeadingZeros x --- could be y just as well

        -- shift divide by 8 states
        numNonCharacters = shiftR numEmptyBits 3
        numChars =  8 - numNonCharacters
    in
    (z, numChars - (popCount u))

-- | andOR64 and or function for Packed64 encoding
andOR64 :: Word64 -> Word64 -> (Word64, Int)
andOR64 x y =
     if  (x .&. y) /= zeroBits then (x .&. y, 0)
     else (x .|. y, 1)

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


-- | packePreorder takes character type, current node (and children preliminary assignments)
-- and parent final assignment and creates final assignment for current node
-- a bit clumsy since uses Goloboff modifications and have to do some of the preOrder pass
-- in Goloboff but not done here
packedPreorder :: CharType -> (V.Vector Word64, V.Vector Word64, V.Vector Word64) -> V.Vector Word64 -> V.Vector Word64
packedPreorder inCharType (leftPrelim, childPrelim, rightPrelim) parentFinal =
    let newStateVect = if inCharType == Packed2       then V.zipWith4 preOrder2 leftPrelim childPrelim rightPrelim parentFinal
                       else if inCharType == Packed4  then V.zipWith4 preOrder4 leftPrelim childPrelim rightPrelim parentFinal
                       else if inCharType == Packed5  then V.zipWith4 preOrder5 leftPrelim childPrelim rightPrelim parentFinal
                       else if inCharType == Packed8  then V.zipWith4 preOrder8 leftPrelim childPrelim rightPrelim parentFinal
                       else if inCharType == Packed64 then V.zipWith4 preOrder64 leftPrelim childPrelim rightPrelim parentFinal
                       else error ("Character type " ++ show inCharType ++ " unrecognized/not implemented")
    in
    newStateVect


-- | preOrder2 performs bitpacked Fitch preorder based on Goloboff 2002
-- less efficient than it could be due to not using Goloboff for post-order
-- assignment so have to calculate some post-order values that would
-- already exist otherwise.  Given that pre-order should be much less frequent than
-- pre-order shouldn't be that bad
preOrder2 :: Word64 -> Word64 -> Word64 -> Word64 -> Word64
preOrder2 leftPrelim childPrelim rightPrelim parentFinal =
    -- post-order stuff to get "temp" state used to calculate final
    let t = leftPrelim .&. rightPrelim

    -- preOrder values
        x2 = parentFinal .&. (complement childPrelim)
        c3 = (mask2B .&. x2) .|. (shiftR (mask2A .&. x2) 1)
        c4 = c3 .|. (shiftL c3 1)

        finalState = (parentFinal .&. (complement c4)) .|. (c4 .&. (childPrelim .|. parentFinal .&. t))
    in
    finalState

-- | preOrder4 from preOrder2 but for 4 states
preOrder4 :: Word64 -> Word64 -> Word64 -> Word64 -> Word64
preOrder4 leftPrelim childPrelim rightPrelim parentFinal =
    -- post-order stuff to get "temp" state used to calculate final
    let x1 = leftPrelim .&. rightPrelim
        y1 = leftPrelim .|. rightPrelim
        c1 = xor mask4B ((mask4B .&. x1) .|. (shiftR (mask4C .&. x1) 1) .|. (shiftR (mask4D .&. x1) 2) .|. (shiftR (mask4E .&. x1) 3))
        c2 = c1 .|. (shiftL c1 1) .|. (shiftL c1 2) .|. (shiftL c1 3)
        t = c2 .|. y1

    -- preOrder values
        x2 = parentFinal .&. (complement childPrelim)
        c3 = (mask4B .&. x2) .|. (shiftR (mask4C .&. x2) 1) .|. (shiftR (mask4D .&. x2) 2) .|. (shiftR (mask4E .&. x2) 3)
        c4 = c3 .|. (shiftL c3 1) .|. (shiftL c3 2) .|. (shiftL c3 3)

        finalState = (parentFinal .&. (complement c4)) .|. (c4 .&. (childPrelim .|. parentFinal .&. t))
    in
    finalState

-- | preOrder5 from preOrder2 but for 5 states
preOrder5 :: Word64 -> Word64 -> Word64 -> Word64 -> Word64
preOrder5 leftPrelim childPrelim rightPrelim parentFinal =
    -- post-order stuff to get "temp" state used to calculate final
    let x1 = leftPrelim .&. rightPrelim
        y1 = leftPrelim .|. rightPrelim
        c1 = xor mask5B ((mask5B .&. x1) .|. (shiftR (mask5C .&. x1) 1) .|. (shiftR (mask5D .&. x1) 2) .|. (shiftR (mask5E .&. x1) 3) .|. (shiftR (mask5F .&. x1) 4))
        c2 = c1 .|. (shiftL c1 1) .|. (shiftL c1 2) .|. (shiftL c1 3) .|. (shiftL c1 4)
        t = c2 .|. y1

    -- preOrder values
        x2 = parentFinal .&. (complement childPrelim)
        c3 = (mask5B .&. x2) .|. (shiftR (mask5C .&. x2) 1) .|. (shiftR (mask5D .&. x2) 2) .|. (shiftR (mask5E .&. x2) 3) .|. (shiftR (mask5F .&. x2) 4)
        c4 = c3 .|. (shiftL c3 1) .|. (shiftL c3 2) .|. (shiftL c3 3) .|. (shiftL c3 4)

        finalState = (parentFinal .&. (complement c4)) .|. (c4 .&. (childPrelim .|. parentFinal .&. t))
    in
    finalState

-- | preOrder8 from preOrder2 but for 8 states
preOrder8 :: Word64 -> Word64 -> Word64 -> Word64 -> Word64
preOrder8 leftPrelim childPrelim rightPrelim parentFinal =
    -- post-order stuff to get "temp" state used to calculate final
    let x1 = leftPrelim .&. rightPrelim
        y1 = leftPrelim .|. rightPrelim
        c1 = xor mask8B ((mask8B .&. x1) .|. (shiftR (mask8C .&. x1) 1) .|. (shiftR (mask8D .&. x1) 2) .|. (shiftR (mask8E .&. x1) 3) .|. (shiftR (mask8F .&. x1) 4) .|. (shiftR (mask8G .&. x1) 5) .|. (shiftR (mask8H .&. x1) 6) .|. (shiftR (mask8I .&. x1) 7))
        c2 = c1 .|. (shiftL c1 1) .|. (shiftL c1 2) .|. (shiftL c1 3) .|. (shiftL c1 4) .|. (shiftL c1 5) .|. (shiftL c1 6) .|. (shiftL c1 7)
        t = c2 .|. y1

    -- preOrder values
        x2 = parentFinal .&. (complement childPrelim)
        c3 = (mask8B .&. x2) .|. (shiftR (mask8C .&. x2) 1) .|. (shiftR (mask8D .&. x2) 2) .|. (shiftR (mask8E .&. x2) 3) .|. (shiftR (mask8F .&. x2) 4) .|. (shiftR (mask8G .&. x2) 5) .|. (shiftR (mask8H .&. x2) 6) .|. (shiftR (mask8I .&. x2) 7)
        c4 = c1 .|. (shiftL c3 1) .|. (shiftL c3 2) .|. (shiftL c3 3) .|. (shiftL c3 4) .|. (shiftL c3 5) .|. (shiftL c3 6) .|. (shiftL c3 7)

        finalState = (parentFinal .&. (complement c4)) .|. (c4 .&. (childPrelim .|. parentFinal .&. t))
    in
    finalState


-- | preOrder64 performs simple Fitch preorder ("up-pass") on Word64
preOrder64 :: Word64 -> Word64 -> Word64 -> Word64 -> Word64
preOrder64 leftPrelim childPrelim rightPrelim parentFinal =
    let a = parentFinal .&. (complement childPrelim)
        b = leftPrelim .&. rightPrelim
        c = parentFinal .|. childPrelim
        d = childPrelim .|. (parentFinal .&. (leftPrelim .|. rightPrelim))
    in
    if a == (zeroBits:: Word64) then parentFinal
    else if b == (zeroBits:: Word64) then c
    else d

{-
Functions for hard-wired 3-way optimization
 basically
            C & P1 & P2 -> if not 0
            else (C & P1) | (C & P2) | (P1 & P2) -> if not 0
            else C | P1 | P2
 bit operations based on Goloboff (2002) for trichotomous trees
-}
 -- | threeWayPacked median 3 for hard-wired networks
 -- this is based on Goloboff (2002) for trichotomous trees
threeWayPacked :: CharType -> V.Vector Word64 -> V.Vector Word64 -> V.Vector Word64 -> V.Vector Word64
threeWayPacked inCharType parent1 parent2 curNode =
    let newStateVect = if inCharType == Packed2       then V.zipWith3 threeWay2 parent1 parent2 curNode
                       else if inCharType == Packed4  then V.zipWith3 threeWay4 parent1 parent2 curNode
                       else if inCharType == Packed5  then V.zipWith3 threeWay5 parent1 parent2 curNode
                       else if inCharType == Packed8  then V.zipWith3 threeWay8 parent1 parent2 curNode
                       else if inCharType == Packed64 then V.zipWith3 threeWay64 parent1 parent2 curNode
                       else error ("Character type " ++ show inCharType ++ " unrecognized/not implemented")
    in
    newStateVect

-- | threeWay2 3-way hardwired optimization for Packed2 Word64
threeWay2 :: Word64 -> Word64 -> Word64 -> Word64
threeWay2 p1 p2 cN =
    let x = p1 .&. p2 .&. cN
        y = (p1 .&. p2) .|. (p1 .&. cN) .|. (p2 .&. cN)
        z = p1 .|. p2 .|. cN
        c1 = xor mask2B ((mask2B .&. x) .|. (shiftR (mask2A .&. x) 1))
        d1 = xor mask2B ((mask2B .&. y) .|. (shiftR (mask2A .&. y) 1))
        c2 = c1 .|. (shiftL c1 1)
        d2 = d1 .|. (shiftL d1 1)
        newState = x .|. (y .&. c2) .|. (z .&. d2)
    in
    newState

-- | threeWay4 3-way hardwired optimization for Packed4 Word64
threeWay4 :: Word64 -> Word64 -> Word64 -> Word64
threeWay4 p1 p2 cN =
    let x = p1 .&. p2 .&. cN
        y = (p1 .&. p2) .|. (p1 .&. cN) .|. (p2 .&. cN)
        z = p1 .|. p2 .|. cN
        c1 = xor mask4B ((mask4B .&. x) .|. (shiftR (mask4C .&. x) 1) .|. (shiftR (mask4D .&. x) 2) .|. (shiftR (mask4E .&. x) 3))
        d1 = xor mask4B ((mask4B .&. y) .|. (shiftR (mask4C .&. y) 1) .|. (shiftR (mask4D .&. y) 2) .|. (shiftR (mask4E .&. y) 3))
        c2 = c1 .|. (shiftL c1 1) .|. (shiftL c1 2) .|. (shiftL c1 3)
        d2 = d1 .|. (shiftL d1 1) .|. (shiftL d1 2) .|. (shiftL d1 3)
        newState = x .|. (y .&. c2) .|. (z .&. d2)
    in
    newState

-- | threeWay5 3-way hardwired optimization for Packed5 Word64
threeWay5 :: Word64 -> Word64 -> Word64 -> Word64
threeWay5 p1 p2 cN =
    let x = p1 .&. p2 .&. cN
        y = (p1 .&. p2) .|. (p1 .&. cN) .|. (p2 .&. cN)
        z = p1 .|. p2 .|. cN
        c1 = xor mask5B ((mask5B .&. x) .|. (shiftR (mask5C .&. x) 1) .|. (shiftR (mask5D .&. x) 2) .|. (shiftR (mask5E .&. x) 3) .|. (shiftR (mask5F .&. x) 4))
        d1 = xor mask5B ((mask5B .&. y) .|. (shiftR (mask5C .&. y) 1) .|. (shiftR (mask5D .&. y) 2) .|. (shiftR (mask5E .&. y) 3) .|. (shiftR (mask5F .&. y) 4))
        c2 = c1 .|. (shiftL c1 1) .|. (shiftL c1 2) .|. (shiftL c1 3) .|. (shiftL c1 4)
        d2 = d1 .|. (shiftL d1 1) .|. (shiftL d1 2) .|. (shiftL d1 3) .|. (shiftL d1 4)
        newState = x .|. (y .&. c2) .|. (z .&. d2)
    in
    newState

-- | threeWay8 3-way hardwired optimization for Packed8 Word64
threeWay8 :: Word64 -> Word64 -> Word64 -> Word64
threeWay8 p1 p2 cN =
    let x = p1 .&. p2 .&. cN
        y = (p1 .&. p2) .|. (p1 .&. cN) .|. (p2 .&. cN)
        z = p1 .|. p2 .|. cN
        c1 = xor mask8B ((mask8B .&. x) .|. (shiftR (mask8C .&. x) 1) .|. (shiftR (mask8D .&. x) 2) .|. (shiftR (mask8E .&. x) 3) .|. (shiftR (mask8F .&. x) 4) .|. (shiftR (mask8G .&. x) 5) .|. (shiftR (mask8H .&. x) 6) .|. (shiftR (mask8I .&. x) 7))
        d1 = xor mask8B ((mask8B .&. y) .|. (shiftR (mask8C .&. y) 1) .|. (shiftR (mask8D .&. y) 2) .|. (shiftR (mask8E .&. y) 3) .|. (shiftR (mask8F .&. y) 4) .|. (shiftR (mask8G .&. y) 5) .|. (shiftR (mask8H .&. y) 6) .|. (shiftR (mask8I .&. y) 7))
        c2 = c1 .|. (shiftL c1 1) .|. (shiftL c1 2) .|. (shiftL c1 3) .|. (shiftL c1 4) .|. (shiftL c1 5) .|. (shiftL c1 6) .|. (shiftL c1 7)
        d2 = d1 .|. (shiftL d1 1) .|. (shiftL d1 2) .|. (shiftL d1 3) .|. (shiftL d1 4) .|. (shiftL d1 5) .|. (shiftL d1 6) .|. (shiftL d1 7)
        newState = x .|. (y .&. c2) .|. (z .&. d2)
    in
    newState

-- | threeWay64 3-way hardwired optimization for straight Word64
threeWay64 :: Word64 -> Word64 -> Word64 -> Word64
threeWay64 p1 p2 cN =
    let x = p1 .&. p2 .&. cN
        y = (p1 .&. p2) .|. (p1 .&. cN) .|. (p2 .&. cN)
        z = p1 .|. p2 .|. cN
    in
    if x /= (zeroBits:: Word64) then x
    else if y /=  (zeroBits :: Word64) then y
    else z


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
        newTaxVectByCharVect = V.fromList $ fmap V.fromList $ L.transpose $ concat recodedSingleVecList

    in
    -- trace ("RNAC: " ++ (show (length recodedSingleVecList, fmap length recodedSingleVecList)) ++ " -> " ++ (show $ fmap length newTaxVectByCharVect) ++ " " ++ (show $ length $ V.fromList $ concat newCharInfoLL))
    (nameBlock, newTaxVectByCharVect, V.fromList $ concat newCharInfoLL)

-- | packNonAdd takes taxon by vector character data and list of character information
-- and returns bit packed and recoded non-additive characters and charInfo
-- input int is character index in block
-- the weight is skipping because of the weight replication in reorganize
-- if characters have non integer weight then they were not reorganized and left
-- as single BV--here as well. Should be very few (if any) of them.
packNonAdd :: V.Vector CharacterData -> CharInfo -> ([[CharacterData]], [CharInfo])
packNonAdd inCharDataV charInfo =
    if (charType charInfo /= NonAdd) || (weight charInfo > 1)  then ([V.toList inCharDataV],[charInfo])
    else
        -- recode non-additive characters
        let leafNonAddV = V.toList $ fmap (snd3 . stateBVPrelim) inCharDataV
            numNonAdd = (length . head) leafNonAddV

            -- split characters into groups by states number 2,4,5,8,64, >64 (excluding missing)
            stateNumDataPairList = fmap (getStateNumber leafNonAddV) [0.. numNonAdd - 1]

            -- sort characters by states number (2, 4, 5, 8, 64, >64 -> 128)
            (state2CharL, state4CharL, state5CharL, state8CharL, state64CharL, state128CharL) = binStateNumber stateNumDataPairList ([],[],[],[],[],[])

            -- make new characters based on state size
            (newStateCharListList, newCharInfoList) = unzip $ (zipWith (makeStateNCharacter charInfo) [2,4,5,8,64,128] [state2CharL, state4CharL, state5CharL, state8CharL, state64CharL, state128CharL] `using` PU.myParListChunkRDS)

        in
        -- trace ("PNA: " ++ (show numNonAdd))  -- (show $ fmap fst stateNumDataPairList) ) --  ++ "\n" ++ (show (newStateCharListList, newCharInfoList) ))
        (newStateCharListList, concat newCharInfoList)

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
    -- trace ("Enter RBV2W64 In: " ++ (show stateNumber) ++ " " ++ (show (length charTaxBVLL, fmap length charTaxBVLL))) (
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
        -- trace ("RBV2W64 Out: " ++ (show $ fmap (snd3 . packedNonAddPrelim) packedDataL))
        (packedDataL, [charInfo {name = newCharName, charType = newCharType}])
        -- )


-- | packIntoWord64 takes a list of bitvectors for a taxon, the state number and number that can be packed into
-- a Word64 and performs appropriate bit settting and shifting to create  Word64
-- paralle looked to bag out here
packIntoWord64 :: Int -> Int -> [[Int]] -> [BV.BitVector] -> CharacterData
packIntoWord64 stateNumber numToPack stateCharacterIndexL inBVList =
    -- get packable chunk of bv and correcsponding state indices
    let packBVList = SL.chunksOf numToPack inBVList
        packIndexLL = SL.chunksOf numToPack stateCharacterIndexL

        -- pack each chunk
        packedWordVect = V.fromList $ zipWith (makeWord64FromChunk stateNumber) packIndexLL packBVList

    in
    -- trace ("PIW64 chunks/values: " ++ (show $ V.length packedWordVect))
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
        -- trace ("MW64FC: " ++ (show subCharacterList) ++ " " ++ (show $ L.foldl1' (.|.) subCharacterList))
        L.foldl1' (.|.) subCharacterList

-- | makeSubCharacter makes sub-character (ie only those bits for states) from single bitvector and shifts appropriate number of bits
-- to make Word64 with sub character bits set and all other bits OFF and in correct bit positions for that sub-character
makeSubCharacter :: Int -> [Int] -> BV.BitVector -> Int -> Word64
makeSubCharacter stateNumber stateIndexList inBV subCharacterIndex =
    -- trace ("Making sub character:" ++ (show stateNumber ++ " " ++ (show stateIndexList) ++ " " ++ (show subCharacterIndex) ++ (show inBV))) (
    let -- get bit of state indices
        bitStates = fmap (testBit inBV) stateIndexList

        -- get index of states when only minimally bit encoded (0101, 0001 -> 11, 01)
        newBitStates = setOnBits (zeroBits :: Word64) bitStates 0
        subCharacter = shiftL newBitStates (subCharacterIndex * stateNumber)
    in
    -- trace ("MSC: " ++ (show subCharacterIndex) ++ " " ++ (show bitStates) ++ " " ++ (show newBitStates) ++ " " ++ (show subCharacter)) (
    -- cna remove this check when working
    if length stateIndexList `notElem` [((fst $ divMod 2 stateNumber) + 1) .. stateNumber] then error ("State number of index list do not match: " ++ (show (stateNumber, length stateIndexList, stateIndexList)))
    else
        subCharacter
    -- )
    -- )

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
        -- trace ("GSIL: " ++ (show indexList))
        indexList


-- | binStateNumber takes a list of pairs of char states number and data column as list of bitvectors and
-- into list for 2,4,5,8,64,>64
binStateNumber :: [(Int, [BV.BitVector])]
               -> ([[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]])
               -> ([[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]])
binStateNumber inPairList (cur2, cur4, cur5, cur8, cur64, cur128) =
    if null inPairList then
        --dont' really need to reverse here but seems hygenic
        trace ("Recoding NonAdditive Characters : " ++ (show (length cur2, length cur4, length cur5, length cur8, length cur64,  length cur128)))
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

getStateNumber :: [V.Vector BV.BitVector] -> Int -> (Int, [BV.BitVector])
getStateNumber  characterDataVV characterIndex =
    -- trace ("GSN:" ++ (show characterIndex) ++ " " ++ (show $ fmap V.length characterDataVV) ++ "\n\t" ++ (show $ fmap (V.! characterIndex) characterDataVV)) (
    if null characterDataVV then (0, [])
    else
        let thisCharV = fmap (V.! characterIndex) characterDataVV
            missingVal = BV.fromBits $ L.replicate (fromEnum $ BV.dimension (head $ thisCharV)) True
            nonMissingStates = filter (/= missingVal) thisCharV
            nonMissingBV = L.foldl1' (.|.) nonMissingStates
            numStates = popCount nonMissingBV

            -- this turns off non-missing bits
            thisCharL = (fmap (.&. nonMissingBV) thisCharV)
        in
        if null nonMissingStates then (1, [])
        else if numStates == 1 then (1, [])
        else if numStates == 2 then (2, thisCharL)
        else if numStates <= 4 then (4, thisCharL)
        else if numStates == 5 then (5, thisCharL)
        else if numStates <= 8 then (8, thisCharL)
        else if numStates <= 64 then (64, thisCharL)
        else (128, thisCharL)

        -- ) -- )
