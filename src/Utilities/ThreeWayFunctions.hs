{- |
Module      :  ThreeWayFunctions.hs
Description :  Module specifying three way optimization functions for use in pre-order
               optimization of HardWired graphs and iterative pass-type optimization for Trees 
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

{-
ToDo:
   Add parallel optimization overblocks and characters?
-}


module Utilities.ThreeWayFunctions  ( threeMedian
                                    ) where

import           Bio.DynamicCharacter
import           Data.Alphabet
import           Data.Bits
import qualified Data.BitVector.LittleEndian as BV
import qualified Data.List                   as L
import           Data.Maybe
import qualified Data.Vector                 as V
import qualified Data.Vector.Generic         as GV
import           Debug.Trace
import qualified DirectOptimization.PreOrder as DOP
import           GeneralUtilities
import qualified GraphOptimization.Medians   as M
import           Types.Types
import qualified Utilities.LocalGraph        as LG

-- | threeWayCharMedian calculates a 3-median for data types in a single character
-- for dynamic characters this is done by 3 min-trees
-- taking best result.  Later true 3-way via ffi can be incorporated
threeWayCharMedian :: GlobalSettings -> AssignmentMethod -> CharInfo -> CharacterData -> CharacterData -> CharacterData -> CharacterData
threeWayCharMedian inGS finalMethod charInfo inA inB inC = 
   inA