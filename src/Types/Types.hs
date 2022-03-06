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

{-# LANGUAGE DerivingStrategies, DeriveGeneric, DeriveAnyClass #-}

module Types.Types where

import           Data.Alphabet
import qualified Data.BitVector.LittleEndian as BV
import qualified Data.MetricRepresentation   as MR
import qualified Data.TCM.Dense              as TCMD
import qualified Data.Text.Lazy              as T
import qualified Data.Text.Short             as ST
import qualified Data.Vector                 as V
import qualified Data.Vector.Storable        as SV
import qualified Data.Vector.Unboxed         as UV
import qualified Data.Time.Clock             as TC
import           Data.Word
import           Foreign.C.Types             (CUInt)
import qualified SymMatrix                   as S
import qualified Utilities.LocalGraph        as LG
import GHC.Generics 
import Control.DeepSeq

-- | Debug Flag
isDebug :: Bool
isDebug = False


-- | Program Version
pgVersion :: String
pgVersion = "0.1"

-- | used for comparing graph costs and edge lengths that are Double
epsilon :: Double
epsilon = 0.0001


-- | infinity is a large Double for use with Graph costs
infinity :: Double
infinity = (read "Infinity") :: Double

-- | Types for timed searches
type Days = Int
type Hours = Int
type Minutes = Int
type Seconds = Int
type Time = (Days, Hours, Minutes, Seconds)


-- | Command types
-- data Argument = String | Double | Int | Bool | Time
--    deriving stock (Show, Eq)
type Argument = (String, String)


--For rename format rename:(a,b,c,...,y,z) => a-y renamed to z
data Instruction = NotACommand | Build | Fuse | Read | Reblock | Refine | Rename | Report | Run | Select | Set | Swap | Search | Support | Transform
    deriving stock (Show, Eq)

-- | Node variety
data NodeType = RootNode | LeafNode | TreeNode | NetworkNode
    deriving stock (Show, Eq)

-- | Edge types
data EdgeType = NetworkEdge | TreeEdge | PendantEdge
    deriving stock (Show, Eq, Ord)

-- | Command type structure
type Command = (Instruction, [Argument])

-- | CharType data type for input characters
data CharType = Binary | Add | NonAdd | Matrix | SlimSeq | WideSeq | HugeSeq | NucSeq | AminoSeq | MatrixApproxSmall | MatrixApproxLarge
    deriving stock (Read, Show, Eq)

-- | types for character classes
nonExactCharacterTypes :: [CharType]
nonExactCharacterTypes = [SlimSeq, WideSeq, HugeSeq, NucSeq, AminoSeq]

exactCharacterTypes :: [CharType]
exactCharacterTypes = [Binary, Add, NonAdd,Matrix , MatrixApproxSmall, MatrixApproxLarge]

-- | Graph types for searching etc.  Can be modified by 'Set command
-- HardWired and SoftWired are network types
-- 'Tree'  would be a single tree in the sense as produced by typical phylogentic
--seqrch programs--no forests
data GraphType = Tree | HardWired | SoftWired
    deriving stock (Show, Eq)

-- | Optimality criterion sets the cost function for graphs and potentially models
data OptimalityCriterion = Parsimony | PMDL
    deriving stock (Show, Eq)

data GraphFactor = NoNetworkPenalty | Wheeler2015Network | PMDLGraph
    deriving stock (Show, Eq)

data RootCost = NoRootCost | Wheeler2015Root | PMDLRoot
    deriving stock (Show, Eq)

-- | Method for makeing final seqeujnce charcatert states assignment 
-- do an DO-based method--more exact but higher time complexity--single preorder
-- pass but worst cae O(n^2) in seqeunce length
-- or assign based on Implied alignment --requires additional post/pre order 
-- traversal but is linear in sequence length
data AssignmentMethod = DirectOptimization | ImpliedAlignment
    deriving stock (Show, Eq)

data SearchData 
    = SearchData
    { instruction :: Instruction
    , arguments   :: [Argument]
    , minGraphCostIn :: VertexCost
    , maxGraphCostIn :: VertexCost
    , numGraphsIn   :: Int
    , minGraphCostOut :: VertexCost
    , maxGraphCostOut :: VertexCost
    , numGraphsOut   :: Int
    , commentString :: String
    , duration    :: Int
    } deriving stock (Show, Eq)


data  GlobalSettings
    = GlobalSettings
    { outgroupIndex       :: Int -- Outgroup terminal index, default 0 (first input leaf)
    , outGroupName        :: T.Text -- Outgroup name
    , optimalityCriterion :: OptimalityCriterion
    , graphType           :: GraphType
    , compressResolutions :: Bool -- "nub" resolutions in softwired graph
    , finalAssignment     :: AssignmentMethod
    , graphFactor         :: GraphFactor -- net penalty/graph complexity
    , rootCost            :: RootCost
    , seed                :: Int -- random seed
    , searchData         :: [SearchData]
    } deriving stock (Show, Eq)

instance NFData GlobalSettings where rnf x = seq x ()


-- | CharInfo information about characters
-- null values for these are in Input.FastAC.hs
--  TCMD.DenseTransitionCostMatrix          => genDiscreteDenseOfDimension (length alphabet)
--  MR.MetricRepresentation Word64          => metricRepresentation <$> TCM.fromRows [[0::Word]]
--  MR.MetricRepresentation BV.BitVector    => metricRepresentation <$> TCM.fromRows [[0::Word]]
data CharInfo = CharInfo { name       :: NameText
                         , charType   :: CharType
                         , activity   :: Bool
                         , weight     :: Double
                         , costMatrix :: S.Matrix Int
                         , slimTCM    :: TCMD.DenseTransitionCostMatrix
                         , wideTCM    :: MR.MetricRepresentation Word64
                         , hugeTCM    :: MR.MetricRepresentation BV.BitVector
                         , alphabet   :: Alphabet ST.ShortText
                         , prealigned :: Bool
                         } deriving stock (Show, Eq)

instance NFData CharInfo where rnf x = seq x ()

-- | Types for vertex information
type VertexCost = Double
type StateCost = Int
type VertexIndex = Int

-- | index of child vertices
type ChildStateIndex = Int

-- | unique bitvector labelling for vertex based on descednent labellings
-- these labels are used for caching, left/right DO optimizaiton
-- thery should be label invariant
-- a hash of sorted input data for leaves
-- will need a map from NameBV to T.Text name (or not)
type NameBV = BV.BitVector

-- | Human legibale name for vertices, characters, and Blocks
type NameText = T.Text

-- | TYpes for Matrix/Sankoff characters
    -- Triple contains info from left and right child--could be only one
    -- use fst then
    -- finals state also vector of triple--but need to keep min cost
    -- for final assignments
    -- filter by local costVect to get 'best' states an each node
type MatrixTriple = (StateCost, [ChildStateIndex], [ChildStateIndex])

-- Only date here that varies by vertex, rest inglobal character info
-- vectors so all data of single type can be grouped together
-- will need to add masks for bit-packing non-additive chars
-- may have to add single assignment for hardwired and IP optimization
-- for approximate sakoff (DO-like) costs can use stateBVPrelim/stateBVFinal
-- for matrix/Saknoff characters-- Vector of vector of States
    --BUT all with same cost matrix/tcm
-- triples (add, no-add, sequence) are to keep children of vertex states for pre-order pass 
-- order is always (left parent median, median, right parent median)
-- do not need for matrix since up pass is a traceback from parent
-- sequence characters are a vector of bitvectors--so only a single seqeunce character
--  per "charctaer" this is so the multi-traversal can take place independently for each
--  sequence character, creating a properly "rooted" tree/graph for each non-exact seqeunce character
-- prelim is created from gapped, final from (via 3-way minimization) parent final and child alignment (2nd and 3rd fields).
-- th ea'alignment' fields hold the im plied alignment data
data CharacterData = CharacterData {   stateBVPrelim      :: (V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector)  -- preliminary for Non-additive chars, Sankoff Approx
                                     -- for Non-additive ans Sankoff/Matrix approximate state
                                     , stateBVFinal       :: V.Vector BV.BitVector
                                     -- for Additive
                                     , rangePrelim        :: (V.Vector (Int, Int), V.Vector (Int, Int), V.Vector (Int, Int))
                                     , rangeFinal         :: V.Vector (Int, Int)
                                     -- for multiple Sankoff/Matrix with slim tcm
                                     , matrixStatesPrelim :: V.Vector (V.Vector MatrixTriple)
                                     , matrixStatesFinal  :: V.Vector (V.Vector MatrixTriple)
                                     -- preliminary for m,ultiple seqeunce chars with same TCM
                                     , slimPrelim         :: SV.Vector CUInt
                                     -- gapped mediasn of left, right, and preliminary used in preorder pass
                                     , slimGapped         :: (SV.Vector CUInt, SV.Vector CUInt, SV.Vector CUInt)
                                     , slimAlignment      :: (SV.Vector CUInt, SV.Vector CUInt, SV.Vector CUInt)
                                     , slimFinal          :: SV.Vector CUInt
                                     , slimIAPrelim       :: (SV.Vector CUInt, SV.Vector CUInt, SV.Vector CUInt)
                                     , slimIAFinal        :: SV.Vector CUInt
                                     -- vector of individual character costs (Can be used in reweighting-ratchet)
                                     , widePrelim         :: UV.Vector Word64
                                     -- gapped median of left, right, and preliminary used in preorder pass
                                     , wideGapped         :: (UV.Vector Word64, UV.Vector Word64, UV.Vector Word64)
                                     , wideAlignment      :: (UV.Vector Word64, UV.Vector Word64, UV.Vector Word64)
                                     , wideFinal          :: UV.Vector Word64
                                     , wideIAPrelim       :: (UV.Vector Word64, UV.Vector Word64, UV.Vector Word64)
                                     , wideIAFinal        :: UV.Vector Word64
                                     -- vector of individual character costs (Can be used in reweighting-ratchet)
                                     , hugePrelim         :: V.Vector BV.BitVector
                                     -- gapped mediasn of left, right, and preliminary used in preorder pass
                                     , hugeGapped         :: (V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector)
                                     , hugeAlignment      :: (V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector)
                                     , hugeFinal          :: V.Vector BV.BitVector
                                     , hugeIAPrelim       :: (V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector)
                                     , hugeIAFinal        :: V.Vector BV.BitVector
                                     -- vector of individual character costs (Can be used in reweighting-ratchet)
                                     , localCostVect      :: V.Vector StateCost
                                     -- weight * V.sum localCostVect
                                     , localCost          :: VertexCost
                                     -- unclear if need vector version
                                     , globalCost         :: VertexCost
                                     } deriving stock (Show, Eq, Generic)

instance NFData CharacterData where rnf x = seq x ()


-- | emptyCharcater useful for intialization and missing data
emptyCharacter :: CharacterData
emptyCharacter = CharacterData { stateBVPrelim = (mempty, mempty, mempty)  -- preliminary for Non-additive chars, Sankoff Approx
                         , stateBVFinal       = mempty
                         -- for Additive
                         , rangePrelim        = (mempty, mempty, mempty)
                         , rangeFinal         = mempty
                         -- for multiple Sankoff/Matrix with sme tcm
                         , matrixStatesPrelim = mempty
                         , matrixStatesFinal  = mempty
                         -- preliminary for m,ultiple seqeunce cahrs with same TCM
                         , slimPrelim         = mempty
                         -- gapped mediasn of left, right, and preliminary used in preorder pass
                         , slimGapped         = (mempty, mempty, mempty)
                         , slimAlignment      = (mempty, mempty, mempty)
                         , slimFinal          = mempty
                         , slimIAPrelim       = (mempty, mempty, mempty)
                         , slimIAFinal        = mempty
                         -- gapped median of left, right, and preliminary used in preorder pass
                         , widePrelim         = mempty
                         -- gapped median of left, right, and preliminary used in preorder pass
                         , wideGapped         = (mempty, mempty, mempty)
                         , wideAlignment      = (mempty, mempty, mempty)
                         , wideFinal          = mempty
                         , wideIAPrelim       = (mempty, mempty, mempty)
                         , wideIAFinal        = mempty
                         -- vector of individual character costs (Can be used in reweighting-ratchet)
                         , hugePrelim         = mempty
                         -- gapped mediasn of left, right, and preliminary used in preorder pass
                         , hugeGapped         = (mempty, mempty, mempty)
                         , hugeAlignment      = (mempty, mempty, mempty)
                         , hugeFinal          = mempty
                         , hugeIAPrelim       = (mempty, mempty, mempty)
                         , hugeIAFinal        = mempty
                         -- vector of individual character costs (Can be used in reweighting-ratchet)
                         , localCostVect      = V.singleton 0
                         -- weight * V.sum localCostVect
                         , localCost          = 0
                         -- unclear if need vector version
                         , globalCost         = 0
                         }

-- | 

-- | type TermData type contians termnal name and list of characters
-- characters as ShortText to save space on input
type TermData = (NameText, [ST.ShortText])
type LeafData = (NameText, V.Vector CharacterData)

-- | VertexBlockData vector over blocks of character data in block (Vector)
-- blocks of character data for  a given vertex
type VertexBlockData = V.Vector (V.Vector CharacterData)

-- | VertexBlockDataMaybe vector over maybe  blocks of character data in block (Vector)
-- blocks of character data for  a given vertex
type VertexBlockDataMaybe = V.Vector (V.Vector (Maybe CharacterData))


-- | ResolutionData contains vertex information for soft-wired network components
-- these are used in the idenitification of minimal cost display trees for a block of
-- data that follow the same display tree
type ResolutionVertexData = V.Vector ResolutionBlockData

-- | ResolutionBlockData contains a list of ResolutionData
-- this list contains all the potential resolutions of a softwired
-- networtk vertex
type ResolutionBlockData = V.Vector ResolutionData

-- | ResolutionData contains individual block information for a given resoluton of soft-wired network components
-- these are used in the idenitification of minimal cost display trees for a block of
-- data that follow the same display tree
-- nodes are VertexInfo for ease of conversion--but nthe info is largely bogus and not to be trusted, same with EdgeInfo
data ResolutionData = ResolutionData { displaySubGraph  :: ([LG.LNode VertexInfo], [LG.LEdge EdgeInfo]) -- holds the post-order display sub-tree for the block
                                     , displayBVLabel   :: NameBV -- For comparison of vertices subtrees, left/right, anmd root leaf inclusion
                                     , displayData      :: V.Vector CharacterData -- data for characters in block
                                     -- list of left, right resolution indices to create current index, used in traceback to get prelminary states
                                     -- and in compressing reolutions to keep only those that result in differnet preliminary states
                                     -- but allowing traceback of resoliutions to get preliminary states
                                     , childResolutions :: [(Maybe Int, Maybe Int)]
                                     , resolutionCost   :: VertexCost -- cost of creating the resolution
                                     , displayCost      :: VertexCost -- cost of that display subtree
                                     } deriving stock (Show, Eq)

-- | VertexInfo type -- vertex information for Decorated Graph
data VertexInfo = VertexInfo { index        :: Int  -- For accessing
                             , bvLabel      :: NameBV -- For comparison of vertices subtrees, left/right
                             , parents      :: V.Vector Int --indegree indices
                             , children     :: V.Vector Int -- outdegree indices
                             , nodeType     :: NodeType -- root, leaf, network, tree
                             , vertName     :: NameText --Text name of vertex either input or HTU#
                             , vertData     :: VertexBlockData -- data as vector of blocks (each a vector of characters)
                             , vertexResolutionData :: V.Vector ResolutionBlockData -- soft-wired network component resolution information for Blocks
                             , vertexCost   :: VertexCost -- local cost of vertex
                             , subGraphCost :: VertexCost -- cost of graph to leaves from the vertex
                             } deriving  stock (Generic, Show, Eq)

instance NFData VertexInfo where rnf x = seq x ()

-- | emptyVertex useful for graph rearrangements
emptyVertexInfo :: VertexInfo
emptyVertexInfo = VertexInfo { index        = (-1)  
                             , bvLabel      = BV.fromBits [False]
                             , parents      = mempty 
                             , children     = mempty 
                             , nodeType     = TreeNode -- root, leaf, network, tree
                             , vertName     = T.pack "EmptyVertex"
                             , vertData     = mempty
                             , vertexResolutionData = mempty
                             , vertexCost   = 0.0 
                             , subGraphCost = 0.0
                             } 

-- | usefule in some cases
dummyNode :: LG.LNode VertexInfo
dummyNode = (-1, emptyVertexInfo)

-- | type edge data, source and sink node indices are fst3 and snd3 fields.
data  EdgeInfo = EdgeInfo   { minLength :: VertexCost
                            , maxLength :: VertexCost
                            , midRangeLength :: VertexCost
                            , edgeType  :: EdgeType
                            } deriving stock (Show, Eq, Ord)

instance NFData EdgeInfo where rnf x = seq x ()

-- | dummyEdge for convenience
dummyEdge :: EdgeInfo
dummyEdge = EdgeInfo    { minLength = 0
                        , maxLength = 0
                        , midRangeLength = 0
                        , edgeType  = TreeEdge
                        } 

-- | DecortatedGraph is the canonical graph contining all final information
-- from preorder traversal trees
-- and post-order info usually from an initial root-based traversal
type DecoratedGraph = LG.Gr VertexInfo EdgeInfo

-- | Type BLockDisplayTree is a Forest of tree components (indegree, outdegree) = (0,1|2),(1,2),(1,0)
-- these are "resolved" from more general graphs
-- will have to allow for indegre=outdegree=1 for dispaly tree generation and reconciliation
-- the vertData field will always have a single Bloclk--teh vecor of blocks will be a vector of
-- BlockDisplayForests.  These woulod better have a single Vector of cChracter info as
-- opposed to the Decorated Tree type, but the record naming and use gets screwed up.
type BlockDisplayForest = LG.Gr VertexInfo EdgeInfo

-- | CharacterTraversalForest is a forest of tree compnents for a single character
-- this is used for non-exact character traversal trees
-- there will always be only a single block and single character even though
-- expresed as Vector fo Vector of Chaarcters.  Would be better as a single character as
-- opposed to the Decorated Tree type, but the record naming and use gets screwed up.
type CharacterTraversalForest = LG.Gr VertexInfo EdgeInfo


-- | type RawGraph is input graphs with leaf and edge labels
type SimpleGraph = LG.Gr NameText Double

-- | Type phylogentic Graph is a graph with
-- cost, optimality value,
-- block display trees, character traversal "foci" (could have multiple)
-- Data optimizations exist in Processed Data
-- Question of wheterh traversal foci shold be in Graph or Data section
-- for now in Graph but could easiily be moved to Processed data
-- May need "EdgeData" like VertexData for heuristics.  UNclear if local scope during SPR/TBR will do it.
--    Fields:
--        1) "Simple" graph with fileds useful for outputting graphs
--        2) Graph optimality value or cost
--        3) Decorated Graph with optimized vertex/Node data
--        4) Vector of display trees for each data Block
--                  root and vertex costs not updated in rerooting so cannot be trusted
--                  Each block can have multiple disdpaly trees so extra vector there
--        5) Vector of traversal foci for each character (Vector of Blocks -> Vector of Characters, a single tree for each character)
--               vector is over blocks, then characters (could have have multiple for each character, but only single tracked here)
--               only important for dynamic (ie non-exact) characters whose costs depend on traversal focus
--               one graph per character  
--        6) Vector of Block Character Information (whihc is a Vector itself) required to properly optimize characters
type PhylogeneticGraph = (SimpleGraph, VertexCost, DecoratedGraph, V.Vector [BlockDisplayForest], V.Vector (V.Vector CharacterTraversalForest), V.Vector (V.Vector CharInfo))

-- | emptyPhylogeneticGraph specifies and empty phylogenetic graph
emptyPhylogeneticGraph :: PhylogeneticGraph
emptyPhylogeneticGraph = (LG.empty, infinity, LG.empty, V.empty, V.empty, V.empty)

-- | RawData type processed from input to be passed to characterData
-- to recode into usable form
-- the format is tuple of a list of taxon-data list tuples and charinfo list.
-- the data list and charinfo list must have the same length
type RawData = ([TermData], [CharInfo])

-- | Processed data is the basic data structire for analysis
-- can be input to functions
-- based on "blocks" that follow same display tree (soft-wired network)
-- each block has a set of characters (for each vertex eventually) and character info
-- the vector of T.Text are the names--leaves form input, internal HTU ++ (show index)
-- ablock are initialy set bu input file, and can be later changed by "set(block,...)"
-- command
-- "Naive" "Optimized" and "Transformed" darta are this type after different processing steps
-- the first and second vectors are size number of leaves, teh third is number of blocks
type ProcessedData = (V.Vector NameText, V.Vector NameBV, V.Vector BlockData)

-- | Block data  is the basic data unit that is optimized on a display tree
--  Block data contain data fo all leaves and all characters in the block
-- it is row, ie vertex dominant
-- it has a bitvector name derived from leaf bitvector labels (union of children)
-- the bitvector name can vary among blocks (for non-leaves) due to alternate display trees
-- a vector of characterer data where assignments and costs reside, and a vector of character info
-- leaves will alwasy be first (indices 0..n-1) for simpler updating of data during graph optimization
-- NameText is the block label used for assignment and reporting output
-- Initially set to input filename of character
--    Fields:
--        1) name of the block--intially taken from input filenames
--        2) vector of vertex/leaf data with vector of character data for each leaf
--        3) Vector of character information for characters in the block
type BlockData = (NameText, V.Vector (V.Vector CharacterData), V.Vector CharInfo)

-- | type SAParams parameter structure for simulated alnnealing and Drifting
data SimulatedAnnealingMethod = SimAnneal | Drift
    deriving stock (Read, Show, Eq)

data SAParams = SAParams { method            :: SimulatedAnnealingMethod
                         , numberSteps       :: Int
                         , currentStep       :: Int
                         , randomIntegerList :: [Int]
                         , rounds            :: Int
                         , driftAcceptEqual  :: Double
                         , driftAcceptWorse  :: Double
                         , driftMaxChanges   :: Int
                         , driftChanges      :: Int
                         } deriving stock (Show, Eq)
