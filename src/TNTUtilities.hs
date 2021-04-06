{- |
Module      :  TNTUtilities.hs
Description :  Module to read tnt input files for phylogenetic analysis
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

{-
Functions to peform the input file reading for PhyG
-}

module TNTUtilities  (getTNTData
                      ) where

import           Types
import           Debug.Trace
import           Data.Char
import           Data.List
import           Data.Maybe
import qualified Data.Text.Lazy  as T
import qualified Data.Text.Short as ST

-- getTNTData tkae file contents and returns raw data and char info form TNT file
getTNTData :: String -> PhyloData
getTNTData inString =
    if null inString then errorWithoutStackTrace ("\n\nTNT input processing error--empty file")
    else 
        let inText = T.strip $ T.pack inString
        in
        if (toLower $ T.head inText) /= 'x' then errorWithoutStackTrace ("\n\nTNT input processing error--must begin with 'xread'")
        else 
            -- look for quoted message
            let singleQuotes = T.count (T.pack "'") inText
                quotedMessage = if singleQuotes == 0 then T.pack "No TNT title message" 
                              else if singleQuotes > 2 then errorWithoutStackTrace ("\n\nTNT input processing error--too many single quotes in title")
                              else (T.split (== '\'') inText) !! 1
                restFile = tail $ T.lines $ (T.split (== '\'') inText) !! 2
                firstLine = head restFile
                numChar = read (T.unpack $ head $ T.words firstLine) :: Int 
                numTax = read (T.unpack $ last $ T.words firstLine) :: Int
            in
            trace ("TNT message : " ++ (T.unpack quotedMessage) ++ " with " ++ (show numTax) ++ " taxa and " ++ (show numChar) ++ " characters") (
            let semiColonLineNumber = findIndex (== T.pack ";") restFile
            in
            if semiColonLineNumber == Nothing then  errorWithoutStackTrace ("\n\nTNT input processing error--can't find ';' to end data block")
            else 
                let dataBlock = filter ((>0).T.length) $ tail $ take (fromJust semiColonLineNumber) restFile
                    numDataLines = length dataBlock
                    (interleaveNumber, interleaveRemainder) = numDataLines `quotRem` numTax
                in
                trace (show dataBlock ++ "\n" ++ show (interleaveNumber, interleaveRemainder)) (
                if interleaveRemainder /= 0 then errorWithoutStackTrace ("\n\nTNT input processing error--number of taxa mis-specified or interleaved format error")
                else 
                    let sortedData = glueInterleave dataBlock numTax numChar []
                    in
                    (sortedData,[])
            ))
            
-- | glueInterleave takes interleves lines and puts them together with name error checking based on number of taxa
glueInterleave :: [T.Text] -> Int -> Int -> [(T.Text, String)] -> [TermData]
glueInterleave lineList numTax numChars curData =
    if null lineList then
        -- recheck num taxa
        -- check chars after process due to ambiguities 
        if length curData /= numTax then error ("Error in glueInterleave: final taxon number error: " ++ show numTax ++ " vs. " ++ (show $ length curData))
        else 
            let nameList = fmap fst curData
                charShortText = (fmap ST.fromString $ fmap snd curData)
            in
            trace ((show $ length curData) ++ " " ++ (show $ length $ snd $ head curData))
            zip nameList (fmap (:[]) charShortText)

    else if length lineList < numTax then error ("Error in glueInterleave: line number error")
    else 
        let thisDataBlock = fmap T.words $ take numTax lineList
            blockNames = fmap head thisDataBlock
            blockStrings = fmap T.unpack $ fmap last thisDataBlock
            canonicalNames = if (length curData > 0) then fmap fst curData
                             else blockNames
            canonicalStrings = fmap snd curData
        in
        --trace (show blockNames ++ "\n" ++ show canonicalNames ++ "\n" ++ show blockStrings ++ "\n" ++ show canonicalStrings) (
        --check for raxon order
        if (blockNames /= canonicalNames) then errorWithoutStackTrace ("\n\nTNT input processing error: interleaved taxon order error")
        else 
            let newChars = if (length curData > 0) then zipWith (++) canonicalStrings blockStrings
                           else blockStrings
            in
            glueInterleave (drop numTax lineList) numTax numChars (zip canonicalNames newChars)
        --)
