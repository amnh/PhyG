{- |
Module      :  Types.hs
Description :  Module specifying data types
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

module Types where

import qualified Data.Text.Lazy  as T
import qualified Data.Text.Short as ST
import qualified LocalGraph as LG

-- | Program Version
pgVersion :: String
pgVersion = "0.1"

-- | Types for timed searches
type Days = Int
type Hours = Int
type Minutes = Int
type Seconds = Int 
type Time = (Days, Hours, Minutes, Seconds)

-- | Command types
-- data Argument = String | Double | Int | Bool | Time
--    deriving (Show, Eq)
type Argument = (String, String)

--For rename format rename:(a,b,c,...,y,z) => a-y renamed to z 

data Instruction = NotACommand | Read | Report | Build | Swap | Refine | Run | Set | Transform | Support | Rename
    deriving (Show, Eq)

type Command = (Instruction, [Argument])

-- | CharType data type for input characters
data CharType = Binary | Add | NonAdd | Matrix | SmallAlphSeq | NucSeq | AminoSeq | GenSeq 
    deriving (Read, Show, Eq)

-- | CharInfo information about characters
data CharInfo = CharInfo { charType :: CharType
                         , activity :: Bool
                         , weight :: Double
                         , costMatrix :: [[Int]]
                         , name :: T.Text
                         , alphabet :: [ST.ShortText]
                         } deriving (Show, Eq)

-- | RawData type processed from input to be passed to characterData
--to recode into usable form
--the format is tuple of a list of taxon-data list tuples and charinfo list.
--the data list and charinfo list must have the same length
type PhyloData = ([TermData], [CharInfo])

-- | type TermData type contians termnal name and list of characters
type TermData = (T.Text, [ST.ShortText])

-- | type RawGraph is input graphs with leaf and edge labels
type SimpleGraph = LG.Gr T.Text Double

-- | type RawGraph is input graphs with leaf and edge labels
-- need to establish this later
-- type LabelledGraph =  LG.Gr a b


