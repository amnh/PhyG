{- |
Module      :  Utilities.hs
Description :  Module specifying utility functions for use with PhyGraph
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

module Utilities.Utilities  where

import qualified Data.Vector as V
import qualified Data.BitVector.LittleEndian as BV
import Data.List.NonEmpty (NonEmpty(..))
import qualified Data.List.NonEmpty as NE
import Bio.Character.Encodable  -- this has DynamicCharacter reference
import Data.Alphabet
import Data.Foldable



-- | dynamicCharacterTo3Vector takes a DYnamicCharacter and returns three Vectors
dynamicCharacterTo3Vector :: DynamicCharacter -> (Word, V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector)
dynamicCharacterTo3Vector (Missing x) = (x, V.empty, V.empty, V.empty)
dynamicCharacterTo3Vector (DC x) = 
    let neVect = V.fromList $ toList x
        (a,b,c) = V.unzip3 neVect
    in
    (0 :: Word, a, b, c)


convertVectorToDynamicCharacter :: V.Vector BV.BitVector -> DynamicCharacter
convertVectorToDynamicCharacter inVector =
    let lenAlph = BV.dimension $ V.head inVector
        arbitraryAlphabet = fromSymbols $ show <$> 0 :| [1 .. lenAlph - 1]
    in
    encodeStream arbitraryAlphabet $ fmap (NE.fromList . f 0 . BV.toBits) . NE.fromList $ toList  inVector
    where 
        f :: Word -> [Bool] -> [String]
        f _ [] = []
        f n (x:xs)
            | x = show n : f (n+1) xs
            | otherwise = f (n+1) xs