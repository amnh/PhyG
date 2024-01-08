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
  , median2PackedUnionField
  , packedPreorder
  , threeWayPacked
  , threeWayPacked'
  , unionPacked
  , minMaxCharDiff
  ) where

import Bio.DynamicCharacter 
import Bio.DynamicCharacter.Element (SlimState, WideState)
import PHANE.Evaluation
import PHANE.Evaluation.Verbosity (Verbosity (..))
import Data.BitVector.LittleEndian qualified as BV
import Data.Bits
import Data.List qualified as L
import Data.List.Split qualified as SL
import Data.Text.Lazy qualified as T
import Data.Vector qualified as V
import Data.Vector.Generic qualified as GV
import Data.Vector.Unboxed qualified as UV
import Data.Word
import GeneralUtilities
import Types.Types
import Utilities.Utilities qualified as U
-- import ParallelUtilities qualified as PU
-- import Debug.Trace


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

Character weights are all = 1 in static characters.  This is ensured by organizeBlockData
in Input.Reorganize.hs basically--characters are multiplied by weight (if integer--otherwise not recoded)
So can check and only recode characters with weight of 1.
-}


{-
Functions for median2 calculations of packed types
These are used in post-order graph traversals and pairwise
distance functions among others.
-}


{-
Masks for various operations and state numbers
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

-- | mask2scN 11(32x00) mask to reveal states of Nth subcharacter in 64Bit Word
mask2sc0 :: Word64
mask2sc0 = 0x3

mask2sc1 :: Word64
mask2sc1 = shiftL mask2sc0 2

mask2sc2 :: Word64
mask2sc2 = shiftL mask2sc0 (2 * 2)

mask2sc3 :: Word64
mask2sc3 = shiftL mask2sc0 (2 * 3)

mask2sc4 :: Word64
mask2sc4 = shiftL mask2sc0 (2 * 4)

mask2sc5 :: Word64
mask2sc5 = shiftL mask2sc0 (2 * 5)

mask2sc6 :: Word64
mask2sc6 = shiftL mask2sc0 (2 * 6)

mask2sc7 :: Word64
mask2sc7 = shiftL mask2sc0 (2 * 7)

mask2sc8 :: Word64
mask2sc8 = shiftL mask2sc0 (2 * 8)

mask2sc9 :: Word64
mask2sc9 = shiftL mask2sc0 (2 * 9)

mask2sc10 :: Word64
mask2sc10 = shiftL mask2sc0 (2 * 10)

mask2sc11 :: Word64
mask2sc11 = shiftL mask2sc0 (2 * 11)

mask2sc12 :: Word64
mask2sc12 = shiftL mask2sc0 (2 * 12)

mask2sc13 :: Word64
mask2sc13 = shiftL mask2sc0 (2 * 13)

mask2sc14 :: Word64
mask2sc14 = shiftL mask2sc0 (2 * 14)

mask2sc15 :: Word64
mask2sc15 = shiftL mask2sc0 (2 * 15)

mask2sc16 :: Word64
mask2sc16 = shiftL mask2sc0 (2 * 16)

mask2sc17 :: Word64
mask2sc17 = shiftL mask2sc0 (2 * 17)

mask2sc18 :: Word64
mask2sc18 = shiftL mask2sc0 (2 * 18)

mask2sc19 :: Word64
mask2sc19 = shiftL mask2sc0 (2 * 19)

mask2sc20 :: Word64
mask2sc20 = shiftL mask2sc0 (2 * 20)

mask2sc21 :: Word64
mask2sc21 = shiftL mask2sc0 (2 * 21)

mask2sc22 :: Word64
mask2sc22 = shiftL mask2sc0 (2 * 22)

mask2sc23 :: Word64
mask2sc23 = shiftL mask2sc0 (2 * 23)

mask2sc24 :: Word64
mask2sc24 = shiftL mask2sc0 (2 * 24)

mask2sc25 :: Word64
mask2sc25 = shiftL mask2sc0 (2 * 25)

mask2sc26 :: Word64
mask2sc26 = shiftL mask2sc0 (2 * 26)

mask2sc27 :: Word64
mask2sc27 = shiftL mask2sc0 (2 * 27)

mask2sc28 :: Word64
mask2sc28 = shiftL mask2sc0 (2 * 28)

mask2sc29 :: Word64
mask2sc29 = shiftL mask2sc0 (2 * 29)

mask2sc30 :: Word64
mask2sc30 = shiftL mask2sc0 (2 * 30)

mask2sc31 :: Word64
mask2sc31 = shiftL mask2sc0 (2 * 31)

-- | mask4scN 1111(16x0000) mask to reveal states of Nth subcharacter in 64Bit Word
mask4sc0 :: Word64
mask4sc0 = 0xF

mask4sc1 :: Word64
mask4sc1 = shiftL mask4sc0 4

mask4sc2 :: Word64
mask4sc2 = shiftL mask4sc0 (4 * 2)

mask4sc3 :: Word64
mask4sc3 = shiftL mask4sc0 (4 * 3)

mask4sc4 :: Word64
mask4sc4 = shiftL mask4sc0 (4 * 4)

mask4sc5 :: Word64
mask4sc5 = shiftL mask4sc0 (4 * 5)

mask4sc6 :: Word64
mask4sc6 = shiftL mask4sc0 (4 * 6)

mask4sc7 :: Word64
mask4sc7 = shiftL mask4sc0 (4 * 7)

mask4sc8 :: Word64
mask4sc8 = shiftL mask4sc0 (4 * 8)

mask4sc9 :: Word64
mask4sc9 = shiftL mask4sc0 (4 * 9)

mask4sc10 :: Word64
mask4sc10 = shiftL mask4sc0 (4 * 10)

mask4sc11 :: Word64
mask4sc11 = shiftL mask4sc0 (4 * 11)

mask4sc12 :: Word64
mask4sc12 = shiftL mask4sc0 (4 * 12)

mask4sc13 :: Word64
mask4sc13 = shiftL mask4sc0 (4 * 13)

mask4sc14 :: Word64
mask4sc14 = shiftL mask4sc0 (4 * 14)

mask4sc15 :: Word64
mask4sc15 = shiftL mask4sc0 (4 * 15)

-- | mask5scN 11111(12x00000) mask to reveal states of Nth subcharacter in 64Bit Word
mask5sc0 :: Word64
mask5sc0 = 0x1F

mask5sc1 :: Word64
mask5sc1 = shiftL mask5sc0 5

mask5sc2 :: Word64
mask5sc2 = shiftL mask5sc0 (5 * 2)

mask5sc3 :: Word64
mask5sc3 = shiftL mask5sc0 (5 * 3)

mask5sc4 :: Word64
mask5sc4 = shiftL mask5sc0 (5 * 4)

mask5sc5 :: Word64
mask5sc5 = shiftL mask5sc0 (5 * 5)

mask5sc6 :: Word64
mask5sc6 = shiftL mask5sc0 (5 * 6)

mask5sc7 :: Word64
mask5sc7 = shiftL mask5sc0 (5 * 7)

mask5sc8 :: Word64
mask5sc8 = shiftL mask5sc0 (5 * 8)

mask5sc9 :: Word64
mask5sc9 = shiftL mask5sc0 (5 * 9)

mask5sc10 :: Word64
mask5sc10 = shiftL mask5sc0 (5 * 10)

mask5sc11 :: Word64
mask5sc11 = shiftL mask5sc0 (5 * 11)

-- | mask8scN 11111111(7x00000000) mask to reveal states of Nth subcharacter in 64Bit Word
mask8sc0 :: Word64
mask8sc0 = 0xFF

mask8sc1 :: Word64
mask8sc1 = shiftL mask8sc0 8

mask8sc2 :: Word64
mask8sc2 = shiftL mask8sc0 (8 * 2)

mask8sc3 :: Word64
mask8sc3 = shiftL mask8sc0 (8 * 3)

mask8sc4 :: Word64
mask8sc4 = shiftL mask8sc0 (8 * 4)

mask8sc5 :: Word64
mask8sc5 = shiftL mask8sc0 (8 * 5)

mask8sc6 :: Word64
mask8sc6 = shiftL mask8sc0 (8 * 6)

mask8sc7 :: Word64
mask8sc7 = shiftL mask8sc0 (8 * 7)

{--
Lists of sub-character masks for operations over packed characters
-}

packed2SubCharList :: [Word64]
packed2SubCharList = [mask2sc0, mask2sc1, mask2sc2, mask2sc3, mask2sc4, mask2sc5, mask2sc6, mask2sc7, mask2sc8, mask2sc9,
                            mask2sc10, mask2sc11, mask2sc12, mask2sc13, mask2sc14, mask2sc15, mask2sc16, mask2sc17, mask2sc18, mask2sc19,
                            mask2sc20, mask2sc21, mask2sc22, mask2sc23, mask2sc24, mask2sc25, mask2sc26, mask2sc27, mask2sc28, mask2sc29,
                            mask2sc30, mask2sc31]

packed4SubCharList :: [Word64]
packed4SubCharList = [mask4sc0, mask4sc1, mask4sc2, mask4sc3, mask4sc4, mask4sc5, mask4sc6, mask4sc7, mask4sc8, mask4sc9,
                            mask4sc10, mask4sc11, mask4sc12, mask4sc13, mask4sc14, mask4sc15]

packed5SubCharList :: [Word64]
packed5SubCharList = [mask5sc0, mask5sc1, mask5sc2, mask5sc3, mask5sc4, mask5sc5, mask5sc6, mask5sc7, mask5sc8, mask5sc9,
                            mask5sc10, mask5sc11]

packed8SubCharList :: [Word64]
packed8SubCharList = [mask8sc0, mask8sc1, mask8sc2, mask8sc3, mask8sc4, mask8sc5, mask8sc6, mask8sc7]

{-
Packed character minimum and maximum length functions
-}

-- | mainMxCharDiff get the approximate minimum and maximum difference in number of states
-- uses masking with &/==
minMaxCharDiff :: CharType -> (Double, Double) -> Word64 -> Word64 -> (Double, Double)
minMaxCharDiff inCharType bitCosts a b =
    let (minVal, maxVal) =  if inCharType == Packed2       then minMaxPacked2 bitCosts a b
                            else if inCharType == Packed4  then minMaxPacked4 bitCosts a b
                            else if inCharType == Packed5  then minMaxPacked5 bitCosts a b
                            else if inCharType == Packed8  then minMaxPacked8 bitCosts a b
                            else if inCharType == Packed64 then minMaxPacked64 bitCosts a b
                            else error ("Character type " <> show inCharType <> " unrecognized/not implemented")
    in
    (minVal, maxVal)


-- | minMaxPacked2 minium and maximum cost 32x2 bit nonadditive character
-- the popcount for equality A/C -> A/C is identical but could be A->C so max 1
-- basically unrolled to make faster
minMaxPacked2 :: (Double, Double) -> Word64 -> Word64 -> (Double, Double)
minMaxPacked2 (lNoChangeCost, lChangeCost) a b =
    let a0 = a .&. mask2sc0
        b0 = b .&. mask2sc0
        max0
          | a0 == (0 :: Word64) = 0
          | (a0 .&. b0) == (0 :: Word64) = 1
          | (a0 == b0) = 0
          | otherwise = 1

        a1 = a .&. mask2sc1
        b1 = b .&. mask2sc1
        max1
          | a1 == (0 :: Word64) = 0
          | (a1 .&. b1) == (0 :: Word64) = 1
          | (a1 == b1) = 0
          | otherwise = 1

        a2 = a .&. mask2sc2
        b2 = b .&. mask2sc2
        max2
          | a2 == (0 :: Word64) = 0
          | (a2 .&. b2) == (0 :: Word64) = 1
          | (a2 == b2) = 0
          | otherwise = 1

        a3 = a .&. mask2sc3
        b3 = b .&. mask2sc3
        max3
          | a3 == (0 :: Word64) = 0
          | (a3 .&. b3) == (0 :: Word64) = 1
          | (a3 == b3) = 0
          | otherwise = 1

        a4 = a .&. mask2sc4
        b4 = b .&. mask2sc4
        max4
          | a4 == (0 :: Word64) = 0
          | (a4 .&. b4) == (0 :: Word64) = 1
          | (a4 == b4) = 0
          | otherwise = 1

        a5 = a .&. mask2sc5
        b5 = b .&. mask2sc5
        max5
          | a5 == (0 :: Word64) = 0
          | (a5 .&. b5) == (0 :: Word64) = 1
          | (a5 == b5) = 0
          | otherwise = 1

        a6 = a .&. mask2sc6
        b6 = b .&. mask2sc6
        max6
          | a6 == (0 :: Word64) = 0
          | (a6 .&. b6) == (0 :: Word64) = 1
          | (a6 == b6) = 0
          | otherwise = 1

        a7 = a .&. mask2sc7
        b7 = b .&. mask2sc7
        max7
          | a7 == (0 :: Word64) = 0
          | (a7 .&. b7) == (0 :: Word64) = 1
          | (a7 == b7) = 0
          | otherwise = 1

        a8 = a .&. mask2sc8
        b8 = b .&. mask2sc8
        max8
          | a8 == (0 :: Word64) = 0
          | (a8 .&. b8) == (0 :: Word64) = 1
          | (a8 == b8) = 0
          | otherwise = 1

        a9 = a .&. mask2sc9
        b9 = b .&. mask2sc9
        max9
          | a9 == (0 :: Word64) = 0
          | (a9 .&. b9) == (0 :: Word64) = 1
          | (a9 == b9) = 0
          | otherwise = 1

        a10 = a .&. mask2sc10
        b10 = b .&. mask2sc10
        max10
          | a10 == (0 :: Word64) = 0
          | (a10 .&. b10) == (0 :: Word64) = 1
          | (a10 == b10) = 0
          | otherwise = 1

        a11 = a .&. mask2sc11
        b11 = b .&. mask2sc11
        max11
          | a11 == (0 :: Word64) = 0
          | (a11 .&. b11) == (0 :: Word64) = 1
          | (a11 == b11) = 0
          | otherwise = 1

        a12 = a .&. mask2sc12
        b12 = b .&. mask2sc12
        max12
          | a12 == (0 :: Word64) = 0
          | (a12 .&. b12) == (0 :: Word64) = 1
          | (a12 == b12) = 0
          | otherwise = 1

        a13 = a .&. mask2sc13
        b13 = b .&. mask2sc13
        max13
          | a13 == (0 :: Word64) = 0
          | (a13 .&. b13) == (0 :: Word64) = 1
          | (a13 == b13) = 0
          | otherwise = 1

        a14 = a .&. mask2sc14
        b14 = b .&. mask2sc14
        max14
          | a14 == (0 :: Word64) = 0
          | (a14 .&. b14) == (0 :: Word64) = 1
          | (a14 == b14) = 0
          | otherwise = 1

        a15 = a .&. mask2sc15
        b15 = b .&. mask2sc15
        max15
          | a15 == (0 :: Word64) = 0
          | (a15 .&. b15) == (0 :: Word64) = 1
          | (a15 == b15) = 0
          | otherwise = 1

        a16 = a .&. mask2sc16
        b16 = b .&. mask2sc16
        max16
          | a16 == (0 :: Word64) = 0
          | (a16 .&. b16) == (0 :: Word64) = 1
          | (a16 == b16) = 0
          | otherwise = 1

        a17 = a .&. mask2sc17
        b17 = b .&. mask2sc17
        max17
          | a17 == (0 :: Word64) = 0
          | (a17 .&. b17) == (0 :: Word64) = 1
          | (a17 == b17) = 0
          | otherwise = 1

        a18 = a .&. mask2sc18
        b18 = b .&. mask2sc18
        max18
          | a18 == (0 :: Word64) = 0
          | (a18 .&. b18) == (0 :: Word64) = 1
          | (a1 == b18) = 0
          | otherwise = 1

        a19 = a .&. mask2sc19
        b19 = b .&. mask2sc19
        max19
          | a19 == (0 :: Word64) = 0
          | (a19 .&. b19) == (0 :: Word64) = 1
          | (a19 == b19) = 0
          | otherwise = 1

        a20 = a .&. mask2sc20
        b20 = b .&. mask2sc20
        max20
          | a20 == (0 :: Word64) = 0
          | (a20 .&. b20) == (0 :: Word64) = 1
          | (a20 == b20) = 0
          | otherwise = 1

        a21 = a .&. mask2sc21
        b21 = b .&. mask2sc21
        max21
          | a21 == (0 :: Word64) = 0
          | (a21 .&. b21) == (0 :: Word64) = 1
          | (a21 == b21) = 0
          | otherwise = 1

        a22 = a .&. mask2sc22
        b22 = b .&. mask2sc22
        max22
          | a22 == (0 :: Word64) = 0
          | (a22 .&. b22) == (0 :: Word64) = 1
          | (a22 == b22) = 0
          | otherwise = 1

        a23 = a .&. mask2sc23
        b23 = b .&. mask2sc23
        max23
          | a23 == (0 :: Word64) = 0
          | (a23 .&. b23) == (0 :: Word64) = 1
          | (a23 == b23) = 0
          | otherwise = 1

        a24 = a .&. mask2sc24
        b24 = b .&. mask2sc24
        max24
          | a24 == (0 :: Word64) = 0
          | (a24 .&. b24) == (0 :: Word64) = 1
          | (a24 == b24) = 0
          | otherwise = 1

        a25 = a .&. mask2sc25
        b25 = b .&. mask2sc25
        max25
          | a25 == (0 :: Word64) = 0
          | (a25 .&. b25) == (0 :: Word64) = 1
          | (a25 == b25) = 0
          | otherwise = 1

        a26 = a .&. mask2sc26
        b26 = b .&. mask2sc26
        max26
          | a26 == (0 :: Word64) = 0
          | (a26 .&. b26) == (0 :: Word64) = 1
          | (a26 == b26) = 0
          | otherwise = 1

        a27 = a .&. mask2sc27
        b27 = b .&. mask2sc27
        max27
          | a27 == (0 :: Word64) = 0
          | (a27 .&. b27) == (0 :: Word64) = 1
          | (a27 == b27) = 0
          | otherwise = 1

        a28 = a .&. mask2sc28
        b28 = b .&. mask2sc28
        max28
          | a28 == (0 :: Word64) = 0
          | (a28 .&. b28) == (0 :: Word64) = 1
          | (a28 == b28) = 0
          | otherwise = 1

        a29 = a .&. mask2sc29
        b29 = b .&. mask2sc29
        max29
          | a29 == (0 :: Word64) = 0
          | (a29 .&. b29) == (0 :: Word64) = 1
          | (a29 == b29) = 0
          | otherwise = 1

        a30 = a .&. mask2sc30
        b30 = b .&. mask2sc30
        max30
          | a30 == (0 :: Word64) = 0
          | (a30 .&. b30) == (0 :: Word64) = 1
          | (a30 == b30) = 0
          | otherwise = 1

        a31 = a .&. mask2sc31
        b31 = b .&. mask2sc31
        max31
          | a31 == (0 :: Word64) = 0
          | (a31 .&. b31) == (0 :: Word64) = 1
          | (a31 == b31) = 0
          | otherwise = 1

        -- sum up values
        (_, minNumNoChange, minNumChange) = andOR2 a b
        maxVal = sum [max0, max1, max2, max3, max4, max5, max6, max7, max8, max9
              , max10, max11, max12, max13, max14, max15, max16, max17, max18, max19
              , max20, max21, max22, max23, max24, max25, max26, max27, max28, max29
              , max30, max31]

    in
    -- trace ("MM2:" <> "\t" <> (showBits a0) <> " " <> (showBits b0) <> "->" <> (showBits $ a0 .&. b0) <> "=>" <> (show max0) <> "\n\t" <> (showBits a10) <> " " <> (showBits b10) <> "->" <> (showBits $ a10 .&. b10) <> "=>" <> (show max10))
    if lNoChangeCost == 0.0 then (fromIntegral minNumChange, fromIntegral maxVal)
    else ((lNoChangeCost * fromIntegral minNumNoChange) + (lChangeCost * fromIntegral minNumChange), (lNoChangeCost * fromIntegral ((32 :: Int) - maxVal)) + (lChangeCost * fromIntegral maxVal))

-- | minMaxPacked4 minium and maximum cost 16x4 bit nonadditive character
-- could add popcount == 1 for equality A/C -> A/C is identical but could be A->C so max 1
-- basically unrolled to make faster
-- any diffference between states gets 1 for max
minMaxPacked4 :: (Double, Double) ->  Word64 -> Word64 -> (Double, Double)
minMaxPacked4 (lNoChangeCost, lChangeCost) a b =
    let a0 = a .&. mask4sc0
        b0 = b .&. mask4sc0
        max0
          | a0 == (0 :: Word64) = 0
          | (a0 .&. b0) == (0 :: Word64) = 1
          | (a0 == b0) = 0
          | otherwise = 1

        a1 = a .&. mask4sc1
        b1 = b .&. mask4sc1
        max1
          | a1 == (0 :: Word64) = 0
          | (a1 .&. b1) == (0 :: Word64) = 1
          | (a1 == b1) = 0
          | otherwise = 1

        a2 = a .&. mask4sc2
        b2 = b .&. mask4sc2
        max2
          | a2 == (0 :: Word64) = 0
          | (a2 .&. b2) == (0 :: Word64) = 1
          | (a2 == b2) = 0
          | otherwise = 1

        a3 = a .&. mask4sc3
        b3 = b .&. mask4sc3
        max3
          | a3 == (0 :: Word64) = 0
          | (a3 .&. b3) == (0 :: Word64) = 1
          | (a3 == b3) = 0
          | otherwise = 1

        a4 = a .&. mask4sc4
        b4 = b .&. mask4sc4
        max4
          | a4 == (0 :: Word64) = 0
          | (a4 .&. b4) == (0 :: Word64) = 1
          | (a4 == b4) = 0
          | otherwise = 1

        a5 = a .&. mask4sc5
        b5 = b .&. mask4sc5
        max5
          | a5 == (0 :: Word64) = 0
          | (a5 .&. b5) == (0 :: Word64) = 1
          | (a5 == b5) = 0
          | otherwise = 1

        a6 = a .&. mask4sc6
        b6 = b .&. mask4sc6
        max6
          | a6 == (0 :: Word64) = 0
          | (a6 .&. b6) == (0 :: Word64) = 1
          | (a6 == b6) = 0
          | otherwise = 1

        a7 = a .&. mask4sc7
        b7 = b .&. mask4sc7
        max7
          | a7 == (0 :: Word64) = 0
          | (a7 .&. b7) == (0 :: Word64) = 1
          | (a7 == b7) = 0
          | otherwise = 1

        a8 = a .&. mask4sc8
        b8 = b .&. mask4sc8
        max8
          | a8 == (0 :: Word64) = 0
          | (a8 .&. b8) == (0 :: Word64) = 1
          | (a8 == b8) = 0
          | otherwise = 1

        a9 = a .&. mask4sc9
        b9 = b .&. mask4sc9
        max9
          | a9 == (0 :: Word64) = 0
          | (a9 .&. b9) == (0 :: Word64) = 1
          | (a9 == b9) = 0
          | otherwise = 1

        a10 = a .&. mask4sc10
        b10 = b .&. mask4sc10
        max10
          | a10 == (0 :: Word64) = 0
          | (a10 .&. b10) == (0 :: Word64) = 1
          | (a10 == b10) = 0
          | otherwise = 1

        a11 = a .&. mask4sc11
        b11 = b .&. mask4sc11
        max11
          | a11 == (0 :: Word64) = 0
          | (a11 .&. b11) == (0 :: Word64) = 1
          | (a11 == b11) = 0
          | otherwise = 1

        a12 = a .&. mask4sc12
        b12 = b .&. mask4sc12
        max12
          | a12 == (0 :: Word64) = 0
          | (a12 .&. b12) == (0 :: Word64) = 1
          | (a12 == b12) = 0
          | otherwise = 1

        a13 = a .&. mask4sc13
        b13 = b .&. mask4sc13
        max13
          | a13 == (0 :: Word64) = 0
          | (a13 .&. b13) == (0 :: Word64) = 1
          | (a13 == b13) = 0
          | otherwise = 1

        a14 = a .&. mask4sc14
        b14 = b .&. mask4sc14
        max14
          | a14 == (0 :: Word64) = 0
          | (a14 .&. b14) == (0 :: Word64) = 1
          | (a14 == b14) = 0
          | otherwise = 1

        a15 = a .&. mask4sc15
        b15 = b .&. mask4sc15
        max15
          | a15 == (0 :: Word64) = 0
          | (a15 .&. b15) == (0 :: Word64) = 1
          | (a15 == b15) = 0
          | otherwise = 1

         -- sum up values
        (_, minNumNoChange, minNumChange) = andOR4 a b
        maxVal = sum [max0, max1, max2, max3, max4, max5, max6, max7, max8, max9
              , max10, max11, max12, max13, max14, max15]

    in
    -- trace ("MM2:" <> "\t" <> (showBits a0) <> " " <> (showBits b0) <> "->" <> (showBits $ a0 .&. b0) <> "=>" <> (show max0) <> "\n\t" <> (showBits a10) <> " " <> (showBits b10) <> "->" <> (showBits $ a10 .&. b10) <> "=>" <> (show max10))
    if lNoChangeCost == 0.0 then (fromIntegral minNumChange, fromIntegral maxVal)
    else ((lNoChangeCost * fromIntegral minNumNoChange) + (lChangeCost * fromIntegral minNumChange), (lNoChangeCost * fromIntegral ((16 :: Int)- maxVal)) + (lChangeCost * fromIntegral maxVal))

-- | minMaxPacked5 minium and maximum cost 12x5 bit nonadditive character
-- the popcount for equality A/C -> A/C is identical but could be A->C so max 1
-- basically unrolled to make faster
minMaxPacked5 :: (Double, Double) ->  Word64 -> Word64 -> (Double, Double)
minMaxPacked5 (lNoChangeCost, lChangeCost) a b =
    let a0 = a .&. mask8sc0
        b0 = b .&. mask5sc0
        max0
          | a0 == (0 :: Word64) = 0
          | (a0 .&. b0) == (0 :: Word64) = 1
          | (a0 == b0) = 0
          | otherwise = 1

        a1 = a .&. mask5sc1
        b1 = b .&. mask5sc1
        max1
          | a1 == (0 :: Word64) = 0
          | (a1 .&. b1) == (0 :: Word64) = 1
          | (a1 == b1) = 0
          | otherwise = 1

        a2 = a .&. mask5sc2
        b2 = b .&. mask5sc2
        max2
          | a2 == (0 :: Word64) = 0
          | (a2 .&. b2) == (0 :: Word64) = 1
          | (a2 == b2) = 0
          | otherwise = 1

        a3 = a .&. mask5sc3
        b3 = b .&. mask5sc3
        max3
          | a3 == (0 :: Word64) = 0
          | (a3 .&. b3) == (0 :: Word64) = 1
          | (a3 == b3) = 0
          | otherwise = 1

        a4 = a .&. mask5sc4
        b4 = b .&. mask5sc4
        max4
          | a4 == (0 :: Word64) = 0
          | (a4 .&. b4) == (0 :: Word64) = 1
          | (a4 == b4) = 0
          | otherwise = 1

        a5 = a .&. mask5sc5
        b5 = b .&. mask5sc5
        max5
          | a5 == (0 :: Word64) = 0
          | (a5 .&. b5) == (0 :: Word64) = 1
          | (a5 == b5) = 0
          | otherwise = 1

        a6 = a .&. mask5sc6
        b6 = b .&. mask5sc6
        max6
          | a6 == (0 :: Word64) = 0
          | (a6 .&. b6) == (0 :: Word64) = 1
          | (a6 == b6) = 0
          | otherwise = 1

        a7 = a .&. mask5sc7
        b7 = b .&. mask5sc7
        max7
          | a7 == (0 :: Word64) = 0
          | (a7 .&. b7) == (0 :: Word64) = 1
          | (a7 == b7) = 0
          | otherwise = 1

        a8 = a .&. mask5sc8
        b8 = b .&. mask5sc8
        max8
          | a8 == (0 :: Word64) = 0
          | (a8 .&. b8) == (0 :: Word64) = 1
          | (a8 == b8) = 0
          | otherwise = 1

        a9 = a .&. mask5sc9
        b9 = b .&. mask5sc9
        max9
          | a9 == (0 :: Word64) = 0
          | (a9 .&. b9) == (0 :: Word64) = 1
          | (a9 == b9) = 0
          | otherwise = 1

        a10 = a .&. mask5sc10
        b10 = b .&. mask5sc10
        max10
          | a10 == (0 :: Word64) = 0
          | (a10 .&. b10) == (0 :: Word64) = 1
          | (a10 == b10) = 0
          | otherwise = 1

        a11 = a .&. mask5sc11
        b11 = b .&. mask5sc11
        max11
          | a11 == (0 :: Word64) = 0
          | (a11 .&. b11) == (0 :: Word64) = 1
          | (a11 == b11) = 0
          | otherwise = 1


         -- sum up values
        (_, minNumNoChange, minNumChange) = andOR5 a b
        maxVal = sum [max0, max1, max2, max3, max4, max5, max6, max7, max8, max9
              , max10, max11]

    in
    -- trace ("MM2:" <> "\t" <> (showBits a0) <> " " <> (showBits b0) <> "->" <> (showBits $ a0 .&. b0) <> "=>" <> (show max0) <> "\n\t" <> (showBits a10) <> " " <> (showBits b10) <> "->" <> (showBits $ a10 .&. b10) <> "=>" <> (show max10))
    if lNoChangeCost == 0.0 then (fromIntegral minNumChange, fromIntegral maxVal)
    else ((lNoChangeCost * fromIntegral minNumNoChange) + (lChangeCost * fromIntegral minNumChange), (lNoChangeCost * fromIntegral ((12 :: Int) - maxVal)) + (lChangeCost * fromIntegral maxVal))

-- | minMaxPacked8 minium and maximum cost 12x5 bit nonadditive character
-- the popcount for equality A/C -> A/C is identical but could be A->C so max 1
-- basically unrolled to make faster
minMaxPacked8 :: (Double, Double) -> Word64 -> Word64 -> (Double, Double)
minMaxPacked8 (lNoChangeCost, lChangeCost) a b =
    let a0 = a .&. mask8sc0
        b0 = b .&. mask8sc0
        max0
          | a0 == (0 :: Word64) = 0
          | (a0 .&. b0) == (0 :: Word64) = 1
          | (a0 == b0) = 0
          | otherwise = 1

        a1 = a .&. mask8sc1
        b1 = b .&. mask8sc1
        max1
          | a1 == (0 :: Word64) = 0
          | (a1 .&. b1) == (0 :: Word64) = 1
          | (a1 == b1) = 0
          | otherwise = 1

        a2 = a .&. mask8sc2
        b2 = b .&. mask8sc2
        max2
          | a2 == (0 :: Word64) = 0
          | (a2 .&. b2) == (0 :: Word64) = 1
          | (a2 == b2) = 0
          | otherwise = 1

        a3 = a .&. mask8sc3
        b3 = b .&. mask8sc3
        max3
          | a3 == (0 :: Word64) = 0
          | (a3 .&. b3) == (0 :: Word64) = 1
          | (a3 == b3) = 0
          | otherwise = 1

        a4 = a .&. mask8sc4
        b4 = b .&. mask8sc4
        max4
          | a4 == (0 :: Word64) = 0
          | (a4 .&. b4) == (0 :: Word64) = 1
          | (a4 == b4) = 0
          | otherwise = 1

        a5 = a .&. mask8sc5
        b5 = b .&. mask8sc5
        max5
          | a5 == (0 :: Word64) = 0
          | (a5 .&. b5) == (0 :: Word64) = 1
          | (a5 == b5) = 0
          | otherwise = 1

        a6 = a .&. mask8sc6
        b6 = b .&. mask8sc6
        max6
          | a6 == (0 :: Word64) = 0
          | (a6 .&. b6) == (0 :: Word64) = 1
          | (a6 == b6) = 0
          | otherwise = 1

        a7 = a .&. mask8sc7
        b7 = b .&. mask8sc7
        max7
          | a7 == (0 :: Word64) = 0
          | (a7 .&. b7) == (0 :: Word64) = 1
          | (a7 == b7) = 0
          | otherwise = 1


       -- sum up values
        (_, minNumNoChange, minNumChange) = andOR8 a b
        maxVal = sum [max0, max1, max2, max3, max4, max5, max6, max7]

    in
    -- trace ("MM2:" <> "\t" <> (showBits a0) <> " " <> (showBits b0) <> "->" <> (showBits $ a0 .&. b0) <> "=>" <> (show max0) <> "\n\t" <> (showBits a10) <> " " <> (showBits b10) <> "->" <> (showBits $ a10 .&. b10) <> "=>" <> (show max10))
    if lNoChangeCost == 0.0 then (fromIntegral minNumChange, fromIntegral maxVal)
    else ((lNoChangeCost * fromIntegral minNumNoChange) + (lChangeCost * fromIntegral minNumChange), (lNoChangeCost * fromIntegral ((7 :: Int)- maxVal)) + (lChangeCost * fromIntegral maxVal))



-- | minMaxPacked64 minium and maximum cost 64 bit nonadditive character
-- the popcount for equality A/C -> A/C is identical but could be A->C so max 1
-- operattion over each sub-character
minMaxPacked64 :: (Double, Double) -> Word64 -> Word64 -> (Double, Double)
minMaxPacked64 (lNoChangeCost, lChangeCost) a b =
    let maxVal = if (a == b) && (popCount a == 1) then (0 :: Int) else (1 :: Int)
        minVal = if (a .&. b) == (0 :: Word64) then (1 :: Int) else (0 :: Int)
    in
    if lNoChangeCost == 0.0 then (fromIntegral minVal, fromIntegral maxVal)
    else ((lNoChangeCost * fromIntegral maxVal) + (lChangeCost * fromIntegral minVal), (lNoChangeCost * fromIntegral minVal) + (lChangeCost * fromIntegral maxVal))


-- | median2Packed takes two characters of packedNonAddTypes
-- and retuns new character data based on 2-median and cost
median2Packed :: CharType -> Double -> (Double, Double) -> CharacterData -> CharacterData -> CharacterData
median2Packed inCharType inCharWeight (thisNoChangeCost, thisChangeCost) leftChar rightChar =
    let (newStateVect, numNoChange, numChange) = if inCharType == Packed2  then median2Word64 andOR2  (snd3 $ packedNonAddPrelim leftChar) (snd3 $ packedNonAddPrelim rightChar)
                                  else if inCharType == Packed4  then median2Word64 andOR4  (snd3 $ packedNonAddPrelim leftChar) (snd3 $ packedNonAddPrelim rightChar)
                                  else if inCharType == Packed5  then median2Word64 andOR5  (snd3 $ packedNonAddPrelim leftChar) (snd3 $ packedNonAddPrelim rightChar)
                                  else if inCharType == Packed8  then median2Word64 andOR8  (snd3 $ packedNonAddPrelim leftChar) (snd3 $ packedNonAddPrelim rightChar)
                                  else if inCharType == Packed64 then median2Word64 andOR64 (snd3 $ packedNonAddPrelim leftChar) (snd3 $ packedNonAddPrelim rightChar)
                                  else error ("Character type " <> show inCharType <> " unrecognized/not implemented")

        -- this for PMDL/ML costs
        newCost = if thisNoChangeCost == 0.0 then inCharWeight * (fromIntegral numChange)
                  else inCharWeight * ((thisChangeCost * (fromIntegral numChange)) + (thisNoChangeCost * (fromIntegral numNoChange)))

        newCharacter = emptyCharacter { packedNonAddPrelim = (snd3 $ packedNonAddPrelim leftChar, newStateVect, snd3 $ packedNonAddPrelim rightChar)
                                      , localCost = newCost
                                      , globalCost =newCost + globalCost leftChar + globalCost rightChar
                                      }
    in
    -- trace ("M2P: " <> (showBitsV $ (snd3 . packedNonAddPrelim) leftChar) <> " " <> (showBitsV $ (snd3 . packedNonAddPrelim) rightChar) <> " -> " <>   (showBitsV $ (snd3 . packedNonAddPrelim) newCharacter) <> " at cost " <> (show newCost))
    newCharacter

-- | median2PackedUnionField takes two characters of packedNonAddTypes
-- and retuns new character data based on 2-median and cost
median2PackedUnionField :: CharType -> Double -> (Double, Double) -> CharacterData -> CharacterData -> CharacterData
median2PackedUnionField inCharType inCharWeight (thisNoChangeCost, thisChangeCost) leftChar rightChar =
    let (newStateVect, numNoChange, numChange) = if inCharType == Packed2  then median2Word64 andOR2  (packedNonAddUnion leftChar) (packedNonAddUnion rightChar)
                                  else if inCharType == Packed4  then median2Word64 andOR4  (packedNonAddUnion leftChar) (packedNonAddUnion rightChar)
                                  else if inCharType == Packed5  then median2Word64 andOR5  (packedNonAddUnion leftChar) (packedNonAddUnion rightChar)
                                  else if inCharType == Packed8  then median2Word64 andOR8  (packedNonAddUnion leftChar) (packedNonAddUnion rightChar)
                                  else if inCharType == Packed64 then median2Word64 andOR64 (packedNonAddUnion leftChar) (packedNonAddUnion rightChar)
                                  else error ("Character type " <> show inCharType <> " unrecognized/not implemented")

        -- this for PMDL/ML costs
        newCost = if thisNoChangeCost == 0.0 then inCharWeight * (fromIntegral numChange)
                  else inCharWeight * ((thisChangeCost * (fromIntegral numChange)) + (thisNoChangeCost * (fromIntegral numNoChange)))

        newCharacter = emptyCharacter { packedNonAddUnion = newStateVect
                                      , localCost = newCost
                                      , globalCost =newCost + globalCost leftChar + globalCost rightChar
                                      }
    in
    newCharacter


-- | unionPacked returns character that is the union (== OR) for bit packed characters
-- of the final fields as preliminary and final
unionPacked :: CharacterData -> CharacterData -> CharacterData
unionPacked charL charR =
    let newVect = UV.zipWith (.|.) (packedNonAddFinal charL) (packedNonAddFinal charR)
    in
    emptyCharacter { packedNonAddPrelim = (newVect, newVect, newVect)
                   , packedNonAddFinal = newVect
                   }


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


-- For all biut packed characters--post order median 2
-- for now--for somwe reason either mis-diagnosed or mis-coded (from White and Holland 2011)
-- the "on" bit in u is reflective of number of intersections not unions.
-- hence subtracting the number of unions from numbers of characters
-- determined by leading OFF bits since packing will likely have ragged edges
-- no not always the pack-able number

-- | median2Word64 driver function for median of two PackedN states
median2Word64 :: (Word64 -> Word64 -> (Word64, Int, Int)) -> UV.Vector Word64 -> UV.Vector Word64 -> (UV.Vector Word64, Int, Int)
median2Word64 andOrFun leftVect rightVect =
    let (stateVect, noChangeVect,changeVect) = UV.unzip3 $ UV.zipWith andOrFun leftVect rightVect
    in
    (stateVect, UV.sum noChangeVect, UV.sum changeVect)

-- | andOR2 and or function for Packed2 encoding
andOR2 :: Word64 -> Word64 -> (Word64, Int, Int)
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
    trace ("AO2 numChars:" <> (show numChars) <> " x & y:" <> (showBits $ x .&. y) <> "\nx .&. y .&. mask2A:" <> (showBits $ (x .&. y .&. mask2A)) <> "\n((x .&. y .&. mask2A) + mask2A):" <> (showBits $ ((x .&. y .&. mask2A) + mask2A))
      <> "\n:(((x .&. y .&. mask2A) + mask2A) .|. (x .&. y)): " <> (showBits $ (((x .&. y .&. mask2A) + mask2A) .|. (x .&. y))) <> "\n:((((x .&. y .&. mask2A) + mask2A) .|. (x .&. y)) .&. mask2B):"
      <> (showBits $ ((((x .&. y .&. mask2A) + mask2A) .|. (x .&. y)) .&. mask2B)) <> "\nu: " <> (showBits u)
      <>"\npc: " <> (show $ popCount u) <> " x:" <> (showBits x) <> " y:" <> (showBits y) <> " => u:" <> (showBits u) <> " z:" <> (showBits z)) -- <> " mask2A:" <> (showBits mask2A) <> " mask2B:" <> (showBits mask2B))
    -}
    (z, popCount u, numChars - popCount u)

-- | andOR4 and or function for Packed4 encoding
andOR4 :: Word64 -> Word64 -> (Word64, Int, Int)
andOR4 x y =
    let u = shiftR ((((x .&. y .&. mask4A) + mask4A) .|. (x .&. y)) .&. mask4B) 3
        z = (x .&. y) .|. ((x .|. y) .&. ((u + mask4A) `xor` mask4B))

        -- get number of characters by checking states (may not be full)
        numEmptyBits = countLeadingZeros x --- could be y just as well

        -- shift divide by 4 states
        numNonCharacters = shiftR numEmptyBits 2
        numChars =  16 - numNonCharacters
    in
    (z, popCount u, numChars - popCount u)

-- | andOR5 and or function for Packed5 encoding
-- potential issue with top 4 bits--not sure on mask5B whether top 4 should be on or OFF.
-- can always mask top 4 with AND 0000111... (0xFFFFFFFFFFFFFFF or 1152921504606846975)
-- to remove bits for counting
-- and calcualted state
andOR5:: Word64 -> Word64 -> (Word64, Int, Int)
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
    -- trace ("AO5 numChars:" <> (show numChars) <> " x & y:" <> (showBits $ x .&. y) <> " u:" <> (showBits u) <> " z:" <> (showBits z) <> " leading 0:" <> (show numEmptyBits) <> " non-chars:" <> (show numNonCharacters) <> " popCount u:" <> (show $ popCount u))
    (z, popCount u, numChars - popCount u)

-- | andOR8 and or function for Packed8 encoding
andOR8 :: Word64 -> Word64 -> (Word64, Int, Int)
andOR8 x y =
    let u = shiftR ((((x .&. y .&. mask8A) + mask8A) .|. (x .&. y)) .&. mask8B) 7
        z = (x .&. y) .|. ((x .|. y) .&. ((u + mask8A) `xor` mask8B))

        -- get number of characters by checking states (may not be full)
        numEmptyBits = countLeadingZeros x --- could be y just as well

        -- shift divide by 8 states
        numNonCharacters = shiftR numEmptyBits 3
        numChars =  8 - numNonCharacters
    in
    (z, popCount u, numChars - popCount u)

-- | andOR64 and or function for Packed64 encoding
-- x `xor` x to make sure all 0 bits
andOR64 :: Word64 -> Word64 -> (Word64, Int, Int)
andOR64 x y =
     if  (x .&. y) /= (x `xor` x)  then (x .&. y, 1, 0)
     else (x .|. y, 0, 1)

{-
Functions for median3 calculations of packed types
These are used in pre-order graph traversals and final state assignment
among others.
-}


-- | packePreorder takes character type, current node (and children preliminary assignments)
-- and parent final assignment and creates final assignment for current node
-- a bit clumsy since uses Goloboff modifications and have to do some of the preOrder pass
-- in Goloboff but not done here
packedPreorder :: CharType -> (UV.Vector Word64, UV.Vector Word64, UV.Vector Word64) -> UV.Vector Word64 -> UV.Vector Word64
packedPreorder inCharType (leftPrelim, childPrelim, rightPrelim) parentFinal =
    let newStateVect
          | inCharType == Packed2 = UV.zipWith4 preOrder2 leftPrelim childPrelim rightPrelim parentFinal
          | inCharType == Packed4 = UV.zipWith4 preOrder4 leftPrelim childPrelim rightPrelim parentFinal
          | inCharType == Packed5 = UV.zipWith4 preOrder5 leftPrelim childPrelim rightPrelim parentFinal
          | inCharType == Packed8 = UV.zipWith4 preOrder8 leftPrelim childPrelim rightPrelim parentFinal
          | inCharType == Packed64 = UV.zipWith4 preOrder64 leftPrelim childPrelim rightPrelim parentFinal
          | otherwise = error ("Character type " <> show inCharType <> " unrecognized/not implemented")
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
        x2 = parentFinal .&. complement childPrelim
        c3 = (mask2A .&. x2) .|. shiftR (mask2B .&. x2) 1
        c4 = c3 .|. shiftL c3 1

        finalState = (parentFinal .&. complement c4) .|. (c4 .&. (childPrelim .|. parentFinal .&. t))
    in
    -- trace ("PO2: " <> " in " <> (show (showBits leftPrelim,  showBits childPrelim, showBits rightPrelim, showBits parentFinal)) <> "->" <> (show $ showBits finalState))
    finalState

-- | preOrder4 from preOrder2 but for 4 states
preOrder4 :: Word64 -> Word64 -> Word64 -> Word64 -> Word64
preOrder4 leftPrelim childPrelim rightPrelim parentFinal =
    -- post-order stuff to get "temp" state used to calculate final
    let x1 = leftPrelim .&. rightPrelim
        y1 = leftPrelim .|. rightPrelim
        c1 = xor mask4E ((mask4E .&. x1) .|. shiftR (mask4D .&. x1) 1 .|. shiftR (mask4C .&. x1) 2 .|. shiftR (mask4B .&. x1) 3)
        c2 = c1 .|. shiftL c1 1 .|. shiftL c1 2 .|. shiftL c1 3
        t = c2 .|. y1

    -- preOrder values
        x2 = parentFinal .&. complement childPrelim
        c3 = (mask4E .&. x2) .|. shiftR (mask4D .&. x2) 1 .|. shiftR (mask4C .&. x2) 2 .|. shiftR (mask4B .&. x2) 3
        c4 = c3 .|. shiftL c3 1 .|. shiftL c3 2 .|. shiftL c3 3

        finalState = (parentFinal .&. complement c4) .|. (c4 .&. (childPrelim .|. parentFinal .&. t))
    in
    finalState

-- | preOrder5 from preOrder2 but for 5 states
preOrder5 :: Word64 -> Word64 -> Word64 -> Word64 -> Word64
preOrder5 leftPrelim childPrelim rightPrelim parentFinal =
    -- post-order stuff to get "temp" state used to calculate final
    let x1 = leftPrelim .&. rightPrelim
        y1 = leftPrelim .|. rightPrelim
        c1 = xor mask5F ((mask5F .&. x1) .|. shiftR (mask5E .&. x1) 1 .|. shiftR (mask5D .&. x1) 2 .|. shiftR (mask5C .&. x1) 3 .|. shiftR (mask5B .&. x1) 4)
        c2 = c1 .|. shiftL c1 1 .|. shiftL c1 2 .|. shiftL c1 3 .|. shiftL c1 4
        t = c2 .|. y1

    -- preOrder values
        x2 = parentFinal .&. complement childPrelim
        c3 = (mask5F .&. x2) .|. shiftR (mask5E .&. x2) 1 .|. shiftR (mask5D .&. x2) 2 .|. shiftR (mask5C .&. x2) 3 .|. shiftR (mask5B .&. x2) 4
        c4 = c3 .|. shiftL c3 1 .|. shiftL c3 2 .|. shiftL c3 3 .|. shiftL c3 4

        finalState = (parentFinal .&. complement c4) .|. (c4 .&. (childPrelim .|. parentFinal .&. t))
    in
    finalState

-- | preOrder8 from preOrder2 but for 8 states
preOrder8 :: Word64 -> Word64 -> Word64 -> Word64 -> Word64
preOrder8 leftPrelim childPrelim rightPrelim parentFinal =
    -- post-order stuff to get "temp" state used to calculate final
    let x1 = leftPrelim .&. rightPrelim
        y1 = leftPrelim .|. rightPrelim
        c1 = xor mask8I ((mask8I .&. x1) .|. shiftR (mask8H .&. x1) 1 .|. shiftR (mask8G .&. x1) 2 .|. shiftR (mask8F .&. x1) 3 .|. shiftR (mask8E .&. x1) 4 .|. shiftR (mask8D .&. x1) 5 .|. shiftR (mask8C .&. x1) 6 .|. shiftR (mask8B .&. x1) 7)
        c2 = c1 .|. shiftL c1 1 .|. shiftL c1 2 .|. shiftL c1 3 .|. shiftL c1 4 .|. shiftL c1 5 .|. shiftL c1 6 .|. shiftL c1 7
        t = c2 .|. y1

    -- preOrder values
        x2 = parentFinal .&. complement childPrelim
        c3 = (mask8I .&. x2) .|. shiftR (mask8H .&. x2) 1 .|. shiftR (mask8G .&. x2) 2 .|. shiftR (mask8F .&. x2) 3 .|. shiftR (mask8E .&. x2) 4 .|. shiftR (mask8D .&. x2) 5 .|. shiftR (mask8C .&. x2) 6 .|. shiftR (mask8B .&. x2) 7
        c4 = c1 .|. shiftL c3 1 .|. shiftL c3 2 .|. shiftL c3 3 .|. shiftL c3 4 .|. shiftL c3 5 .|. shiftL c3 6 .|. shiftL c3 7

        finalState = (parentFinal .&. complement c4) .|. (c4 .&. (childPrelim .|. parentFinal .&. t))
    in
    finalState


-- | preOrder64 performs simple Fitch preorder ("up-pass") on Word64
preOrder64 :: Word64 -> Word64 -> Word64 -> Word64 -> Word64
preOrder64 leftPrelim childPrelim rightPrelim parentFinal =
    let a = parentFinal .&. complement childPrelim
        b = leftPrelim .&. rightPrelim
        c = parentFinal .|. childPrelim
        d = childPrelim .|. (parentFinal .&. (leftPrelim .|. rightPrelim))
    in
    if a == (leftPrelim `xor` leftPrelim :: Word64) then parentFinal
    else if b == (leftPrelim `xor` leftPrelim :: Word64) then c
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
threeWayPacked :: CharType -> UV.Vector Word64 -> UV.Vector Word64 -> UV.Vector Word64 -> UV.Vector Word64
threeWayPacked inCharType parent1 parent2 curNode =
    let newStateVect
          | inCharType == Packed2 = UV.zipWith3 threeWay2 parent1 parent2 curNode
          | inCharType == Packed4 = UV.zipWith3 threeWay4 parent1 parent2 curNode
          | inCharType == Packed5 = UV.zipWith3 threeWay5 parent1 parent2 curNode
          | inCharType == Packed8 = UV.zipWith3 threeWay8 parent1 parent2 curNode
          | inCharType == Packed64 = UV.zipWith3 threeWay64 parent1 parent2 curNode
          | otherwise = error ("Character type " <> show inCharType <> " unrecognized/not implemented")
    in
    newStateVect

-- | threeWayPacked' median 3 for hard-wired networks
-- this uses lists of masks so likely slower than Goloboff
-- this approach could also be used for min/max to be simpler but alos likelu slower since previous is
-- manually unrolled
threeWayPacked' :: CharType -> UV.Vector Word64 -> UV.Vector Word64 -> UV.Vector Word64 -> UV.Vector Word64
threeWayPacked' inCharType parent1 parent2 curNode =
    let newStateVect
          | inCharType == Packed2 = UV.zipWith3 (threeWayNWord64 packed2SubCharList) parent1 parent2 curNode
          | inCharType == Packed4 = UV.zipWith3 (threeWayNWord64 packed4SubCharList) parent1 parent2 curNode
          | inCharType == Packed5 = UV.zipWith3 (threeWayNWord64 packed5SubCharList) parent1 parent2 curNode
          | inCharType == Packed8 = UV.zipWith3 (threeWayNWord64 packed8SubCharList) parent1 parent2 curNode
          | inCharType == Packed64 = UV.zipWith3 threeWay64 parent1 parent2 curNode
          | otherwise = error ("Character type " <> show inCharType <> " unrecognized/not implemented")
    in
    newStateVect


-- | threeWayNWord64 3-way hardwired optimization for Packed N Word64
-- non-additive character--maps over sub-characters with appropriate masks
-- lists of subcharacters with all ovther bits OFF are created via masks
-- then zipped over threeway function and ORed to create 32 bit final state
-- this is an alternate approach to the three node optimization of Golobiff below.
-- both should yield same result-- this is polymoprhic and simple--but not parallel
--as in Goloboff so likely slower
threeWayNWord64 :: [Word64] -> Word64 -> Word64 -> Word64 -> Word64
threeWayNWord64 packedSubCharList p1 p2 cN =
       let p1SubCharList = fmap (p1 .&.) packedSubCharList
           p2SubCharList = fmap (p2 .&.) packedSubCharList
           cNSubCharList = fmap (cN .&.) packedSubCharList
           threeWayList  = zipWith3 threeWay64 p1SubCharList p2SubCharList cNSubCharList
       in
       L.foldl1' (.|.) threeWayList


-- | threeWay2 3-way hardwired optimization for Packed2 Word64
-- but used on subCharacters
threeWay2 :: Word64 -> Word64 -> Word64 -> Word64
threeWay2 p1 p2 cN =
    let x = p1 .&. p2 .&. cN
        y = (p1 .&. p2) .|. (p1 .&. cN) .|. (p2 .&. cN)
        z = p1 .|. p2 .|. cN
        c1 = xor mask2B ((mask2B .&. x) .|. shiftR (mask2A .&. x) 1)
        d1 = xor mask2B ((mask2B .&. y) .|. shiftR (mask2A .&. y) 1)
        c2 = c1 .|. shiftL c1 1
        d2 = d1 .|. shiftL d1 1
        newState = x .|. (y .&. c2) .|. (z .&. d2)
    in
    newState

-- | threeWay4 3-way hardwired optimization for Packed4 Word64
threeWay4 :: Word64 -> Word64 -> Word64 -> Word64
threeWay4 p1 p2 cN =
    let x = p1 .&. p2 .&. cN
        y = (p1 .&. p2) .|. (p1 .&. cN) .|. (p2 .&. cN)
        z = p1 .|. p2 .|. cN
        c1 = xor mask4B ((mask4B .&. x) .|. shiftR (mask4C .&. x) 1 .|. shiftR (mask4D .&. x) 2 .|. shiftR (mask4E .&. x) 3)
        d1 = xor mask4B ((mask4B .&. y) .|. shiftR (mask4C .&. y) 1 .|. shiftR (mask4D .&. y) 2 .|. shiftR (mask4E .&. y) 3)
        c2 = c1 .|. shiftL c1 1 .|. shiftL c1 2 .|. shiftL c1 3
        d2 = d1 .|. shiftL d1 1 .|. shiftL d1 2 .|. shiftL d1 3
        newState = x .|. (y .&. c2) .|. (z .&. d2)
    in
    newState

-- | threeWay5 3-way hardwired optimization for Packed5 Word64
threeWay5 :: Word64 -> Word64 -> Word64 -> Word64
threeWay5 p1 p2 cN =
    let x = p1 .&. p2 .&. cN
        y = (p1 .&. p2) .|. (p1 .&. cN) .|. (p2 .&. cN)
        z = p1 .|. p2 .|. cN
        c1 = xor mask5B ((mask5B .&. x) .|. shiftR (mask5C .&. x) 1 .|. shiftR (mask5D .&. x) 2 .|. shiftR (mask5E .&. x) 3 .|. shiftR (mask5F .&. x) 4)
        d1 = xor mask5B ((mask5B .&. y) .|. shiftR (mask5C .&. y) 1 .|. shiftR (mask5D .&. y) 2 .|. shiftR (mask5E .&. y) 3 .|. shiftR (mask5F .&. y) 4)
        c2 = c1 .|. shiftL c1 1 .|. shiftL c1 2 .|. shiftL c1 3 .|. shiftL c1 4
        d2 = d1 .|. shiftL d1 1 .|. shiftL d1 2 .|. shiftL d1 3 .|. shiftL d1 4
        newState = x .|. (y .&. c2) .|. (z .&. d2)
    in
    newState

-- | threeWay8 3-way hardwired optimization for Packed8 Word64
threeWay8 :: Word64 -> Word64 -> Word64 -> Word64
threeWay8 p1 p2 cN =
    let x = p1 .&. p2 .&. cN
        y = (p1 .&. p2) .|. (p1 .&. cN) .|. (p2 .&. cN)
        z = p1 .|. p2 .|. cN
        c1 = xor mask8B ((mask8B .&. x) .|. shiftR (mask8C .&. x) 1 .|. shiftR (mask8D .&. x) 2 .|. shiftR (mask8E .&. x) 3 .|. shiftR (mask8F .&. x) 4 .|. shiftR (mask8G .&. x) 5 .|. shiftR (mask8H .&. x) 6 .|. shiftR (mask8I .&. x) 7)
        d1 = xor mask8B ((mask8B .&. y) .|. shiftR (mask8C .&. y) 1 .|. shiftR (mask8D .&. y) 2 .|. shiftR (mask8E .&. y) 3 .|. shiftR (mask8F .&. y) 4 .|. shiftR (mask8G .&. y) 5 .|. shiftR (mask8H .&. y) 6 .|. shiftR (mask8I .&. y) 7)
        c2 = c1 .|. shiftL c1 1 .|. shiftL c1 2 .|. shiftL c1 3 .|. shiftL c1 4 .|. shiftL c1 5 .|. shiftL c1 6 .|. shiftL c1 7
        d2 = d1 .|. shiftL d1 1 .|. shiftL d1 2 .|. shiftL d1 3 .|. shiftL d1 4 .|. shiftL d1 5 .|. shiftL d1 6 .|. shiftL d1 7
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
    if x /= (p1 `xor` p1 :: Word64) then x
    else if y /=  (p1 `xor` p1  :: Word64) then y
    else z


{-
Functions to encode ("pack") non-additive characters into new Word64 characters
based on their number of states
-}

-- | packData takes input data and creates a variety of bit-packed data types
-- to increase efficiency and reduce footprint of non-additive characters
-- that are encoded as bitvectors
packNonAdditiveData :: GlobalSettings -> ProcessedData -> PhyG ProcessedData
packNonAdditiveData inGS (nameVect, bvNameVect, blockDataVect) = 
    -- need to check if this blowws out memory on big data sets (e.g. genomic)
    let -- parallel setup
        action ::  BlockData -> PhyG BlockData
        action = recodeNonAddCharacters inGS
    in do
        newBlockDataList <- getParallelChunkTraverse >>= \pTraverse ->
            action `pTraverse` V.toList blockDataVect
            --PU.seqParMap (parStrategy $ strictParStrat inGS) (recodeNonAddCharacters inGS) (V.toList blockDataVect) --  could be an option to save memory etc
    
        pure (nameVect, bvNameVect, V.fromList newBlockDataList)

-- | recodeNonAddCharacters takes block data, goes through characters
-- and recodes NonAdditive.
-- Concat and list for charInfoV because new characters can be created
-- and newCharInfo then as well, could be multiple per input 'charcater'
recodeNonAddCharacters :: GlobalSettings -> BlockData -> PhyG BlockData
recodeNonAddCharacters inGS (nameBlock, charDataVV, charInfoV) =
    let numChars = V.length charInfoV

        -- create vector of single characters with vector of taxon data of sngle character each
        singleCharVectList = V.toList $ fmap (U.getSingleCharacter charDataVV) (V.fromList [0.. numChars - 1])

        packAction :: (V.Vector CharacterData, CharInfo) -> PhyG ([[CharacterData]], [CharInfo])
        packAction = packNonAddPair inGS 
    in do
        -- bit pack the nonadd
        result <- getParallelChunkTraverse >>= \pTraverse ->
            packAction `pTraverse` zip singleCharVectList (V.toList charInfoV)
        let (recodedSingleVecList, newCharInfoLL) = unzip result -- $ zipWith (packNonAdd inGS) singleCharVectList (V.toList charInfoV)

        -- recreate BlockData, tacxon dominant structure
        let newTaxVectByCharVect = V.fromList $ fmap V.fromList $ L.transpose $ concat recodedSingleVecList

        -- trace ("RNAC: " <> (show (length recodedSingleVecList, fmap length recodedSingleVecList)) <> " -> " <> (show $ fmap length newTaxVectByCharVect) <> " " <> (show $ length $ V.fromList $ concat newCharInfoLL))
        pure (nameBlock, newTaxVectByCharVect, V.fromList $ concat newCharInfoLL)

-- | packNonAddPair is a wrapper arpounf packNonAdd
packNonAddPair :: GlobalSettings ->  (V.Vector CharacterData, CharInfo) -> PhyG ([[CharacterData]], [CharInfo])
packNonAddPair inGS (inCharDataV, charInfo) = packNonAdd inGS inCharDataV charInfo

-- | packNonAdd takes (vector of taxa) by character data and list of character information
-- and returns bit packed and recoded non-additive characters and charInfo
-- input int is character index in block
-- the weight is skipping because of the weight replication in reorganize
-- if characters have non integer weight then they were not reorganized and left
-- as single BV--here as well. Should be very few (if any) of them.
packNonAdd ::GlobalSettings ->  V.Vector CharacterData -> CharInfo -> PhyG ([[CharacterData]], [CharInfo])
packNonAdd inGS inCharDataV charInfo =
    -- trace ("PNA in weight: " <> (show $ weight charInfo)) (
    if charType charInfo /= NonAdd then pure ([V.toList inCharDataV],[charInfo])
    else
        -- recode non-additive characters
        let leafNonAddV = V.toList $ fmap (snd3 . stateBVPrelim) inCharDataV

              -- there is a problem with this index--they should all be the same but there are two classes in some cases
              -- I believe due to missing data
            -- numNonAdd = (length . head) leafNonAddV
            numNonAdd = minimum $ fmap length leafNonAddV

            -- parallel setup
            stateAction ::  Int -> (Int, [BV.BitVector])
            stateAction = getStateNumber leafNonAddV

            charAction ::  (Int, [[BV.BitVector]]) -> PhyG ([CharacterData], [CharInfo])
            charAction =  makeStateNCharacterTuple inGS charInfo


        in do
            -- split characters into groups by states number 2,4,5,8,64, >64 (excluding missing)
            statePar <- getParallelChunkMap
            let stateNumDataPairList = statePar stateAction [0.. numNonAdd - 1]
                -- PU.seqParMap (parStrategy $ strictParStrat inGS)   (getStateNumber leafNonAddV) [0.. numNonAdd - 1]

            -- sort characters by states number (2, 4, 5, 8, 64, >64 -> 128)
            let (state2CharL, state4CharL, state5CharL, state8CharL, state64CharL, state128CharL) = binStateNumber stateNumDataPairList ([],[],[],[],[],[])

            -- make new characters based on state size
            result <- getParallelChunkTraverse >>= \pTraverse ->
                charAction `pTraverse` zip [2,4,5,8,64,128] [state2CharL, state4CharL, state5CharL, state8CharL, state64CharL, state128CharL]
            let (newStateCharListList, newCharInfoList) = unzip result
              -- (PU.seqParMap (parStrategy $ strictParStrat inGS)  (makeStateNCharacterTuple inGS charInfo) (zip [2,4,5,8,64,128] [state2CharL, state4CharL, state5CharL, state8CharL, state64CharL, state128CharL]))

            -- this here in case recoding and removing constant (with missing) characters yields no data to be bitpacked
            if (L.foldl1' (&&) $ fmap null [state2CharL, state4CharL, state5CharL, state8CharL, state64CharL, state128CharL]) then 
              pure ([V.toList inCharDataV],[charInfo])
            else pure (newStateCharListList, concat newCharInfoList)
            -- )

-- | makeStateNCharacterTuple is a wrapper for makeStateNCharacter to allow for parMap use
makeStateNCharacterTuple ::  GlobalSettings -> CharInfo -> (Int, [[BV.BitVector]]) -> PhyG ([CharacterData], [CharInfo])
makeStateNCharacterTuple inGS charInfo (stateNumber, charDataLL) = makeStateNCharacter inGS charInfo stateNumber charDataLL

-- | makeStateNCharacter takes a list of characters each of which is a list of taxon character values and
-- creates a new character of all characters for give taxon and packs (64/ state number) characters into a 64 bit Word64
-- via chuncksOf--or if 64, not packing, if 128 stays bitvector
-- check for non-sequential states (A,T) or (0,2) etc
-- return is list of taxa x single new (packed) character
makeStateNCharacter :: GlobalSettings -> CharInfo -> Int -> [[BV.BitVector]] -> PhyG ([CharacterData], [CharInfo])
makeStateNCharacter inGS charInfo stateNumber charDataLL = do
    result <-  if stateNumber > 64 then recodeBV2BV inGS charInfo charDataLL
               else if stateNumber == 64 then recodeBV2Word64Single inGS charInfo charDataLL
               else recodeBV2Word64 inGS charInfo stateNumber charDataLL
  
    pure result 

-- | recodeBV2BV take a list of BV.bitvector non-add characters and creates a list (taxa)
-- BV non-additive characters of type NonAdd.
-- this results in a single character and charInfo in list so can be concatenated
-- and removed if empty
recodeBV2BV :: GlobalSettings -> CharInfo -> [[BV.BitVector]] -> PhyG ([CharacterData], [CharInfo])
recodeBV2BV inGS charInfo charTaxBVLL =
    if null charTaxBVLL then pure ([],[])
    else
        let -- convert to taxon by characgter data list
            newStateList = makeNewCharacterData charTaxBVLL

            -- rename with thype
            newCharName = T.append (name charInfo) $ T.pack "LargeState"

            -- create new characters for each taxon
            newCharDataList = fmap (makeNewData emptyCharacter) newStateList
        in do
        pure (newCharDataList, [charInfo {name = newCharName, charType = NonAdd, noChangeCost = (fst . bcgt64) inGS, changeCost = (snd . bcgt64) inGS}])
        where makeNewData a b = a {stateBVPrelim = (b,b,b), stateBVFinal = b}


-- | recodeBV2Word64Single take a list of BV.bitvector non-add characters and creates a list (taxa)
-- of Word64 unpacked non-additive characters of type Packed64.
-- this results in a single character and charInfo in list so can be concatenated
-- and removed if empty
-- Aassumes a leaf only sets snd3
recodeBV2Word64Single :: GlobalSettings -> CharInfo -> [[BV.BitVector]] -> PhyG ([CharacterData], [CharInfo])
recodeBV2Word64Single inGS charInfo charTaxBVLL =
    if null charTaxBVLL then pure ([],[])
    else
        let newCharName = T.append (name charInfo) $ T.pack "64State"

            -- convert BV to Word64
            taxWord64BLL = fmap (fmap BV.toUnsignedNumber) charTaxBVLL

            -- convert to taxon by character data lisyt
            newStateList = makeNewCharacterData taxWord64BLL

            -- make new character data
            newCharDataList = fmap (makeNewData emptyCharacter) newStateList
        in do
        pure (newCharDataList, [charInfo {name = newCharName, charType = Packed64, noChangeCost = (fst . bc64) inGS, changeCost = (snd . bc64) inGS}])
        where makeNewData a b = a {packedNonAddPrelim = (b,b,b), packedNonAddFinal = b}

-- | makeNewCharacterData takes a list of characters, each of which is a list of taxon states
-- of type a (bitvector or Word64) and returns a list of taxa each of which is a vector
-- of type a charactyer data
-- generic verctor so can have Unboxed V and Boxed V
makeNewCharacterData :: (GV.Vector v a) => [[a]] -> [v a]
makeNewCharacterData charByTaxSingleCharData  =
    let taxonByCharL = L.transpose charByTaxSingleCharData
        taxonByCharV = fmap GV.fromList taxonByCharL
    in
    taxonByCharV

-- | recodeBV2Word64 take a list of BV.bitvector non-add characters and the states number of creates
-- Word64 representaions where subcharacters are created and shifted to proper positions and ORd
-- to create packed reresentation--new character types Packed2, Packed4, Packed5, and Packed8.
-- this results in a single character and charInfo in list so can be concatenated
-- and removed if empty
recodeBV2Word64 ::GlobalSettings ->  CharInfo -> Int -> [[BV.BitVector]] -> PhyG ([CharacterData], [CharInfo])
recodeBV2Word64 inGS charInfo stateNumber charTaxBVLL =
    -- trace ("Enter RBV2W64 In: " <> (show stateNumber) <> " " <> (show (length charTaxBVLL, fmap length charTaxBVLL))) (
    if null charTaxBVLL then pure ([],[])
    else
        let newCharType = if stateNumber == 2 then Packed2
                          else if stateNumber == 4 then Packed4
                          else if stateNumber == 5 then Packed5
                          else if stateNumber == 8 then Packed8
                          else error ("State number " <> (show stateNumber) <> " not to be packed in recodeBV2Word64")

            newCharName = T.append (name charInfo) $ T.pack ((show stateNumber) <> "State")

            -- get number of characters that can be packed into Word64 for that state number
            numCanPack = fst $ divMod 64 stateNumber

            -- convert to taxon by character data list
            taxCharBVLL = L.transpose charTaxBVLL

            -- parallel setup
            stateAction :: [BV.BitVector] -> [Int]
            stateAction = getStateIndexList 

            packAction :: [[Int]] -> [BV.BitVector] -> CharacterData
            packAction = packIntoWord64 stateNumber numCanPack 

        in do
            -- get state index list for all characters (could be non sequential 0|2; A|T etc)
            statePar <- getParallelChunkMap
            let stateIndexLL = statePar stateAction charTaxBVLL
                -- PU.seqParMap (parStrategy $ strictParStrat inGS)  getStateIndexList charTaxBVLL

            -- convert data each taxon into packedWord64
            packPar <- getParallelChunkMap
            let packedDataL = packPar (packAction stateIndexLL) taxCharBVLL 
                -- PU.seqParMap (parStrategy $ strictParStrat inGS)  (packIntoWord64 stateNumber numCanPack stateIndexLL) taxCharBVLL

            -- get noChange and Change cost for char type
            let (lNoChangeCost, lChangeCost) = if stateNumber == 2 then bc2 inGS
                                           else if stateNumber == 4 then bc4 inGS
                                           else if stateNumber == 5 then bc5 inGS
                                           else if stateNumber == 8 then bc8 inGS
                                           else error ("Change/NoChange costs for state number " <> (show stateNumber) <> " not to be set in recodeBV2Word64")

        
            -- trace ("RBV2W64 Out: " <> (show $ fmap (snd3 . packedNonAddPrelim) packedDataL))
            -- trace ("RBV2W64: " <> (show (lNoChangeCost, lChangeCost)))
            pure (packedDataL, [charInfo {name = newCharName, charType = newCharType, noChangeCost = lNoChangeCost, changeCost = lChangeCost}])
            -- )


-- | packIntoWideState takes a list of bitvectors for a taxon, the state number and number that can be packed into
-- a WideState and performs appropriate bit settting and shifting to create  WideState
-- paralle looked to bag out here
packIntoWord64 :: Int -> Int -> [[Int]] -> [BV.BitVector] -> CharacterData
packIntoWord64 stateNumber numToPack stateCharacterIndexL inBVList =
    -- get packable chunk of bv and correcsponding state indices
    let packBVList = SL.chunksOf numToPack inBVList
        packIndexLL = SL.chunksOf numToPack stateCharacterIndexL

        -- pack each chunk
        packedWordVect = UV.fromList $ zipWith (makeWord64FromChunk stateNumber) packIndexLL packBVList

    in  emptyCharacter { packedNonAddPrelim = (packedWordVect, packedWordVect, packedWordVect)
                   , packedNonAddFinal = packedWordVect
                   }


-- | makeWord64FromChunk takes a list (= chunk) of bitvectors and creates bit subcharacter (Word64)
-- with adjacent bits for each BV in chunk.  It then bit shifts the appropriate number of bits for each member
-- of the chunk and finally ORs all (64/stateNumber) Word64s to  make the final packed representation
makeWord64FromChunk ::  Int -> [[Int]] ->  [BV.BitVector] -> Word64
makeWord64FromChunk stateNumber stateIndexLL bvList =
    if null bvList then (0 :: Word64)
    else
        let subCharacterList = zipWith3 (makeSubCharacter stateNumber) stateIndexLL bvList [0..(length bvList - 1)]
        in
        -- trace ("MW64FC: " <> (show subCharacterList) <> " " <> (show $ L.foldl1' (.|.) subCharacterList))
        L.foldl1' (.|.) subCharacterList

-- | makeSubCharacter makes sub-character (ie only those bits for states) from single bitvector and shifts appropriate number of bits
-- to make Word64 with sub character bits set and all other bits OFF and in correct bit positions for that sub-character
makeSubCharacter :: Int -> [Int] -> BV.BitVector -> Int -> Word64
makeSubCharacter stateNumber stateIndexList inBV subCharacterIndex =
    -- trace ("Making sub character:" <> (show stateNumber <> " " <> (show stateIndexList) <> " " <> (show subCharacterIndex) <> (show inBV))) (
    let -- get bit of state indices
        bitStates = fmap (testBit inBV) stateIndexList

        -- get index of states when only minimally bit encoded (0101, 0001 -> 11, 01)
        newBitStates = setOnBits (0 `xor` 0 :: Word64) bitStates 0
        subCharacter = shiftL newBitStates (subCharacterIndex * stateNumber)
    in
    -- trace ("MSC: " <> (show subCharacterIndex) <> " " <> (show bitStates) <> " " <> (show newBitStates) <> " " <> (show subCharacter)) $
    -- cna remove this check when working
    
    if length stateIndexList `notElem` [(fst (divMod 2 stateNumber) + 1) .. stateNumber] then error ("State number of index list do not match: " <> show (stateNumber, length stateIndexList, stateIndexList))
    else
           subCharacter
    -- )
    -- )

-- | setOnBits recursively sets On bits in a list of Bool
setOnBits :: Word64 -> [Bool] -> Int -> Word64
setOnBits baseVal onList bitIndex =
    if null onList then baseVal
    else
        let newVal = if head onList then setBit baseVal bitIndex
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
            onList = fmap (testBit inBV) [0.. finiteBitSize inBV - 1]
            onIndexPair = zip onList [0.. finiteBitSize inBV - 1]
            indexList = (snd <$> filter fst onIndexPair)

        in
        -- trace ("GSIL: " <> (show indexList))
        indexList


-- | binStateNumber takes a list of pairs of char states number and data column as list of bitvectors and
-- into list for 2,4,5,8,64,>64
binStateNumber :: [(Int, [BV.BitVector])]
               -> ([[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]])
               -> ([[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]],[[BV.BitVector]])
binStateNumber inPairList (cur2, cur4, cur5, cur8, cur64, cur128) =
    if null inPairList then
        --dont' really need to reverse here but seems hygenic
        -- trace ("Recoding NonAdditive Characters : " <> (show (length cur2, length cur4, length cur5, length cur8, length cur64,  length cur128)))
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
-- checks for min length vchars to make those missing-- can happen with implied alignment recoding
getStateNumber :: [V.Vector BV.BitVector] -> Int -> (Int, [BV.BitVector])
getStateNumber  characterDataVV characterIndex  =
    -- trace ("GSN:" <> (show characterIndex) <> " " <> (show $ fmap V.length characterDataVV) <> "\n\t" <> (show $ fmap (V.! characterIndex) characterDataVV)) (
    if null characterDataVV then (0, [])
    else
        let thisCharV = fmap (V.! characterIndex) characterDataVV
            missingVal = L.foldl1' (.|.) thisCharV
            nonMissingStates = filter (/= missingVal) thisCharV
            nonMissingBV = L.foldl1' (.|.) nonMissingStates
            numStates = popCount nonMissingBV

            -- this turns off non-missing bits
            thisCharL = fmap (.&. nonMissingBV) thisCharV
        in
        -- trace ("GSN:" <> (show nonMissingStates) <> "\nNMBV " <> (show nonMissingBV) <> " MBV " <> (show missingVal) <> " -> " <> (show numStates) ) (
        (if null nonMissingStates || (numStates == 1) then (1, []) else (if numStates == 2 then (2, thisCharL)
        else if numStates <= 4 then (4, thisCharL)
        else if numStates == 5 then (5, thisCharL)
        else if numStates <= 8 then (8, thisCharL)
        else if numStates <= 64 then (64, thisCharL)
        else (128, thisCharL)))

        -- ) -- )
