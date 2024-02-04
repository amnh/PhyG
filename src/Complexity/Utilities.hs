{- |
Module      :  Utilities
Description :  Functions to generate (algorithmic) complexity of objects
Copyright   :  (c) 2018-2020 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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


module Complexity.Utilities
  (  getInformationContent
  ,  split2Matrix
  --,  invertMatrixExt
  ,  ppMatrix
  )  where

import Codec.Compression.GZip qualified as GZ
import Complexity.Huffman
import Data.ByteString.Lazy qualified as B
import Data.String.Encode qualified as E
import Data.List
-- import Complexity.Constants
-- import Numeric.LinearAlgebra

-- | occurencesInList counts occurences elements in list--makes result double so can divide later
occurencesInList :: Eq a => [a] -> a -> Double
occurencesInList elementList element = fromIntegral $ (length . filter (== element)) elementList

-- | getTotalBits takes occurence and bit encoding to return total bits of program
getTotalBits :: Double -> [Double] -> [Double] -> Double
getTotalBits currentSum occurenceList bitList =
  if null occurenceList || null bitList then currentSum
  else
      let occurrences = head occurenceList
          bitCoding = head bitList
      in
      --trace ("n=" ++ show occurrences ++ " b=" ++ show bitCoding)
      getTotalBits (currentSum + (occurrences * bitCoding)) (tail occurenceList) (tail bitList)

-- | getShannon gets Shannon entropy bits by freqency for string (list of symbols)
getShannon :: String -> Int
getShannon inCharList =
  if null inCharList then error "Input list to getShannon is empty"
  else
    let totalSymbols = length inCharList
        symbolList = nub inCharList
        symbolOccurences = fmap (occurencesInList inCharList) symbolList
        symbolFrequency = fmap (/ fromIntegral totalSymbols) symbolOccurences
        symbolBits = fmap (logBase 2.0) symbolFrequency
        totalBits = getTotalBits 0.0 symbolOccurences symbolBits --This and above line could be a single function
    in
    --trace ("There were " ++ show (length symbolList) ++ " unique symbols in program.")
    ceiling $ abs totalBits

-- | getHuffCode takes the code map list and rturns the binary code
getHuffCode :: Code Char -> Char -> [Bit]
getHuffCode huffCode inChar =
  if null huffCode then error "Code empty or character not found"
  else
    let (candChar, bitList) = head huffCode
    in
    if candChar == inChar then bitList
    else getHuffCode (tail huffCode) inChar

bitToChar :: Bit -> Char
bitToChar inBit =
  if inBit == Zero then '0'
  else '1'

-- | getHuffman takes input list of chars and returns Huffman code representation
getHuffman :: String -> (Int, String)
getHuffman inString =
  if null inString then error "Input list to getHuffman is empty"
  else
    let symbolList = nub inString
        symbolOccurences = fmap (occurencesInList inString) symbolList
        symbolOccurencePairs = zip symbolList symbolOccurences
        huffTree = huffman symbolOccurencePairs
        huffCodes = codewords huffTree
        bitList = concatMap (getHuffCode huffCodes) inString
        bitString = fmap bitToChar bitList
    in
    --trace (ppCode huffCodes)
    (length bitList, bitString)

-- | getInformation takes a program string and returns Shannon, Huffman bits, Huffman binary program, and comopressed length bits .
-- for min on compression--short things can increase in bits when compressed. 
-- Using Shannon bits, but could use Huffman since Shannon can't be realized
getInformationContent :: String -> (Double, Int, String, Double)
getInformationContent programString =
  if null programString then error "Empty program in getInformation"
  else
    let (huffBits, huffBinary) =  getHuffman programString
        compressedStream = GZ.compressWith GZ.defaultCompressParams {GZ.compressLevel = GZ.bestCompression} (E.convertString programString)
        shannonBits = getShannon programString
        compressedBits = min shannonBits (8 * (length $ B.unpack compressedStream))
    in
    (fromIntegral shannonBits, huffBits, huffBinary, fromIntegral compressedBits)
    
-- | split2Matrix takes alist and splits after n elements to make it a list of lists
split2Matrix :: Int -> [Double] -> [[Double]]
split2Matrix lengthPieces inList
  | null inList = []
  | lengthPieces > length inList = error ("List too short to split " ++ show lengthPieces ++ " " ++ show inList)
  | otherwise =
      take lengthPieces inList : split2Matrix lengthPieces (drop lengthPieces inList)

{-
-- | exportable versions of functions
--invertMatrixExt uses library to invert matrix but uses simple
invertMatrixExt :: [[Double]] -> [[Double]]
invertMatrixExt a =
  if null a then error "Null matrix to invert"
  else
    let dimension = length a
        aFlat = concat a
        a' = matrix dimension aFlat
        invA' = inv a'
    in
    toLists invA'
-}

-- | ppMatrix makes a nice(er) string for matrix
ppMatrix :: [[Double]] -> String
ppMatrix inMatrix =
  if null inMatrix then "[]"
  else
    let rowS = fmap show inMatrix
        rowStrings = fmap (++ "\n") rowS
        matString = concat rowStrings
    in
    matString

