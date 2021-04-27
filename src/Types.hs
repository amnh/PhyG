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
import qualified Data.BitVector as BV
import qualified Data.Vector    as V
import qualified SymMatrix as S 

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
data Instruction = NotACommand | Read | Report | Build | Swap | Refine | Run | Set | Transform | Support | Rename | Select
    deriving (Show, Eq)

-- | Node variety
data NodeType = Root | Leaf | Tree | Network
    deriving (Show, Eq)

type Command = (Instruction, [Argument])

-- | CharType data type for input characters
data CharType = Binary | Add | NonAdd | Matrix | SmallAlphSeq | NucSeq | AminoSeq | GenSeq | MatrixApproxSmall | MatrixApproxLarge
    deriving (Read, Show, Eq)


-- | CharInfo information about characters
-- will likely add full (for small alphabets) and hashMap (for large alphabets) tcm's here
data CharInfo = CharInfo { name :: T.Text
                         , charType :: CharType
                         , activity :: Bool
                         , weight :: Double
                         , costMatrix :: S.Matrix Int
                         , alphabet :: [ST.ShortText]
                         , prealigned :: Bool
                         } deriving (Show, Eq)

-- | Types for vertex information
type VertexCost = Double
type StateCost = Int

-- | index of child vertices 
type ChildStateIndex = Int

-- | unique bitvector labelling for vertex based on descednent labellings
-- these labels are used for caching, left/right DO optimizaiton
-- thery should be label invariant
-- a hash of sorted input data for leaves
-- will need a map from NameBV to T.Text name (or not)
type NameBV = BV.BV

-- | Human legibale name for vertices, charcaters, and Blocks
type NameText = T.Text

-- | TYpes for Matrix/Sankoff charcaters
type MatrixTriple = (StateCost, [ChildStateIndex], [ChildStateIndex])
type MatrixState = V.Vector MatrixTriple

-- Only date here that varies by vertex, rest inglobal charcater info
-- vectors so all data of single type can be grouped together
-- will need to add masks for bit-packing non-additive chars
-- may have to add single
-- for approximate (DO-like) costs can use stateBVPrelim/stateBVFinal
data CharacterData = CharacterData {   stateBVPrelim :: V.Vector BV.BV  -- for Non-additive, seqeujnce, and Sankoff/Matrix approximate state
                                     , minRangePrelim :: V.Vector Int -- for Additive
                                     , maxRangePrelim :: V.Vector Int -- for Additive
                                     -- triple for Sankoff optimization--cost, left and right descendant states (could be multiple)
                                     , matrixStatesPrelim :: V.Vector MatrixState -- for Sankoff/Matrix
                                     , stateBVFinal :: V.Vector BV.BV  -- for Non-additive ans Sankoff/Matrix approximate state
                                     , minRangeFinal :: V.Vector Int -- for Additive
                                     , maxRangeFinal :: V.Vector Int -- for Additive
                                     , matrixStatesFinal :: V.Vector (V.Vector (StateCost)) -- for Sankoff/Matrix  keeps delta to "best" states 0 or > 0
                                     , localCostVect :: V.Vector StateCost 
                                     , localCost :: VertexCost -- weight * V.sum localCostVect
                                     , globalCost :: VertexCost -- unclear if need vector version
                                     } deriving (Show, Eq)

-- | type TermData type contians termnal name and list of characters
-- characters as ShortText to save space on input
type TermData = (NameText, [ST.ShortText])
type LeafData = (NameText, V.Vector CharacterData)


data VertexType = VertexType { index :: Int  -- For accessing
                             , bvLabel :: BV.BV -- For comparison of vertices subtrees, etc
                             , parents :: V.Vector Int --indegree indices
                             , children :: V.Vector Int -- outdegree indices
                             , variety :: NodeType  
                             } deriving (Show, Eq)


data EdgeType = EdgeType { sourceIndex :: Int  -- For accessing
                         , sinkIndex :: Int  -- For accessing
                         , minLength :: Double
                         , maxLength :: Double
                         } deriving (Show, Eq)

-- | Type BLockDisplayTree is a Forest of tree components (indegree, outdegree) = (0,1|2),(1,2),(1,0)
-- these are "resolved" from more general graphs
-- will have to allow for indegre=outdegree=1 for disply tryee genration and rteconciliation
type BlockDisplayForest = LG.Gr VertexType EdgeType 

-- | Type BlockFoci are a vector for each character (in a block usually) of a vector of edges since there may be more than 1 "best" focus
-- static charcatsr all are fine--so length 1 defualt value
-- dynamic characters its the edge of traversal focus, a psuedo-root 
type BlockFoci = V.Vector (V.Vector (LG.LEdge EdgeType))

-- | type RawGraph is input graphs with leaf and edge labels
type SimpleGraph = LG.Gr NameText Double

-- | Type phylogentic Graph is a graph with
-- cost, optimality value, 
-- block display trees, character traversal "foci" (could have multiple)
-- Data optimizations exist in Processed Data
-- Question of wheterh traversal foci shold be in Graph or Data section
-- for now in Graph but could easiily be moved to Processed data
-- May need "EdgeData" like VertexData for heuristics.  UNclear if local scope during SPR/TBR will do it.
--	Fields:
--		1) "Simple" graph with fileds useful for outputting graphs  
--		2) Graph optimality value or cost
--		3) Vector of display trees for each data Block
--		4) Vector of traversal foci for each character (Blocks -> Vector of Chars -> Vector of traversal edges)
-- 		5) Data associated with that tree, alwasy same for leaves (Could keep them separate), but vary for HTUs
type PhylogeneticGraph = (SimpleGraph, VertexCost, V.Vector BlockDisplayForest, V.Vector BlockFoci, ProcessedData)


-- | type RawGraph is input graphs with leaf and edge labels
-- need to establish this later
-- type LabelledGraph =  LG.Gr a b
-- | RawData type processed from input to be passed to characterData
-- to recode into usable form
-- the format is tuple of a list of taxon-data list tuples and charinfo list.
-- the data list and charinfo list must have the same length
type RawData = ([TermData], [CharInfo])

-- | Processed data is the basic data structire for analysis
-- can be input to functions
-- based on "blocks" that follow same display tree (soft-wired network)
-- each block has a set of characters (for each vertex eventually) and character info
-- the vector of T.Text are teh names--leaves form input, internal HTU ++ (show index)
-- ablock are initialy set bu input file, and can be later changed by "set(block,...)"
-- command
-- "Naive" "Opyimized" and "Transformed" darta are this type after differnet processing steps
type ProcessedData = (V.Vector NameText, V.Vector BlockData)

-- | Block data  is the basic data unit that is optimized on a display tree
-- it is row, ie vertex dominant
-- it has a bitvector name derived from leaf bitvector labels (union of children)
-- the bitvector name can vary among blocks (for non-leaves) due to alternate display trees 
-- a vector of characterer data where assignments and costs reside, and a vector of character info
-- leaves will alwasy be first (indices 0..n-1) for simpler updating of data during graph optimization
-- NameText is the block label used for assignment and reporting output
-- Initially set to input filename of character
--	Fields:
--		1) name of the block--intially taken from input filenames
--		2) vector of vertiex data with (node Bitvector--for dsisplayu tree stuff, vector of character states)
--		3) Vector of charcater information for characters in the block 
type BlockData = (NameText, V.Vector (NameBV, V.Vector CharacterData), V.Vector CharInfo)



