{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE OverloadedStrings #-}

{- |
Module specifying data types
-}
module Types.Types (
    module Types.Types,
) where

import Bio.DynamicCharacter (HugeDynamicCharacter, OpenDynamicCharacter, SlimDynamicCharacter, WideDynamicCharacter)
import Bio.DynamicCharacter.Element (HugeState, SlimState, WideState)
import Control.DeepSeq
import Control.Monad.IO.Class (MonadIO)
import Control.Parallel.Strategies
import Data.Alphabet
import Data.BitVector.LittleEndian qualified as BV
import Data.InfList qualified as IL
import Data.List.NonEmpty (NonEmpty (..))
import Data.MetricRepresentation as MR
import Data.TCM qualified as TCM
import Data.TCM.Dense qualified as TCMD
import Data.Text.Lazy qualified as T
import Data.Text.Short qualified as ST
import Data.Vector qualified as V
import Data.Vector.Storable qualified as SV
import Data.Vector.Unboxed qualified as UV
import Data.Word (Word64)
import GHC.Generics
import PHANE.Evaluation
import SymMatrix qualified as S
import Utilities.LocalGraph qualified as LG


-- | Debug Flag
isDebug ∷ Bool
isDebug = False


-- | Program Version
pgVersion ∷ String
pgVersion = "0.1"


-- | used for comparing graph costs and edge lengths that are Double
epsilon ∷ Double
epsilon = 0.0001


-- | infinity is a large Double for use with Graph costs
infinity ∷ Double
infinity = read "Infinity" ∷ Double


{- | maxAddStatesToRecode maximum size of addditive character to recode into
non-additive characters 65 can fit in 4 WideState since nstates - 1 binaries
prob could be bigger based on cost of optimizing additive versus but this
seems a reasonable number (prob should be timed to verify)
-}
maxAddStatesToRecode ∷ Int
maxAddStatesToRecode = 129


{- |
Core monad transformer stack for evaluating computations within the application PhyG.
-}
type PhyG = Evaluation ()


-- | Types for timed searches
type Days = Int


type Hours = Int


type Minutes = Int


type Seconds = Int


type Time = (Days, Hours, Minutes, Seconds)


{- | Command types
data Argument = String | Double | Int | Bool | Time
   deriving stock (Show, Eq)
-}
type Argument = (String, String)


-- For rename format rename:(a,b,c,...,y,z) => a-y renamed to z
data Instruction
    = NotACommand
    | Build
    | Fuse
    | Read
    | Reblock
    | Refine
    | Rename
    | Report
    | Run
    | Select
    | Set
    | Swap
    | Search
    | Support
    | Transform
    deriving stock (Show, Eq, Ord)


-- | Node variety
data NodeType = RootNode | LeafNode | TreeNode | NetworkNode | In1Out1
    deriving stock (Show, Eq)


-- | Edge types
data EdgeType = NetworkEdge | TreeEdge | PendantEdge
    deriving stock (Show, Eq, Ord)


-- | Command type structure
type Command = (Instruction, [Argument])


-- | CharType data type for input characters
data CharType
    = Add
    | NonAdd
    | Matrix
    | SlimSeq
    | WideSeq
    | HugeSeq
    | NucSeq
    | AminoSeq
    | AlignedSlim
    | AlignedWide
    | AlignedHuge
    | Packed2
    | Packed4
    | Packed5
    | Packed8
    | Packed64
    deriving stock (Read, Show, Eq)


{- | Evaluation strategries for a graph, full muti-traversal, single-trwraversal based on outpgroup rooting,
static apporximation of non-exact (unaligned sequence) charcaters.
-}
data GraphEvaluation = MultiTraverse | SingleTraverse | StaticApproximation
    deriving stock (Read, Show, Eq)


-- | bandit types
data BanditType = SearchBandit | GraphBandit
    deriving stock (Read, Show, Eq)


-- non additive bit packed types (64 not really 'packed' but treated as if were)
-- these are not entered but are created by transforming existing non-additive characters
packedNonAddTypes ∷ [CharType]
packedNonAddTypes = [Packed2, Packed4, Packed5, Packed8, Packed64]


-- aligned not in here because they are not reorganized, and would screw up reroot optimization
exactCharacterTypes ∷ [CharType]
exactCharacterTypes = [Add, NonAdd, Matrix] <> packedNonAddTypes


-- | types for character classes
nonExactCharacterTypes ∷ [CharType]
nonExactCharacterTypes = [SlimSeq, WideSeq, HugeSeq, NucSeq, AminoSeq] -- , AlignedSlim, AlignedWide, AlignedHuge]


-- prealigned types
prealignedCharacterTypes ∷ [CharType]
prealignedCharacterTypes = [AlignedSlim, AlignedWide, AlignedHuge]


-- sequence types
sequenceCharacterTypes ∷ [CharType]
sequenceCharacterTypes = nonExactCharacterTypes <> prealignedCharacterTypes


{- | Graph types for searching etc.  Can be modified by 'Set command
HardWired and SoftWired are network types
'Tree'  would be a single tree in the sense as produced by typical phylogentic
seqrch programs--no forests
-}
data GraphType = Tree | HardWired | SoftWired
    deriving stock (Show, Eq)


{- | Optimality criterion sets the cost function for graphs and potentially models
likelihood form is the "Self Information" in context of Kolmogorov complexity/MDL/PMDL
MAPA is Integrate Likelihood/dBaysian stuff via TCM modification
-}
data OptimalityCriterion = Parsimony | PMDL | SI | MAPA | NCM
    deriving stock (Show, Eq)


data GraphFactor = NoNetworkPenalty | Wheeler2015Network | Wheeler2023Network | PMDLGraph
    deriving stock (Show, Eq)


data RootCost = NoRootCost | MAPARoot | NCMRoot | PMDLRoot | SIRoot | Wheeler2015Root
    deriving stock (Show, Eq)


data SoftWiredAlgorithm = Naive | ResolutionCache
    deriving stock (Show, Eq)


data ParallelStrategy = R0 | RSeq | RPar | RDeepSeq
    deriving stock (Show, Eq)


{- | Method for makeing final seqeujnce charactert states assignment
do an DO-based method--more exact but higher time complexity--single preorder
pass but worst cae O(n^2) in seqeunce length
or assign based on Implied alignment --requires additional post/pre order
traversal but is linear in sequence length
-}
data AssignmentMethod = DirectOptimization | ImpliedAlignment
    deriving stock (Show, Eq)


data SearchData = SearchData
    { instruction ∷ Instruction
    , arguments ∷ [Argument]
    , minGraphCostIn ∷ VertexCost
    , maxGraphCostIn ∷ VertexCost
    , numGraphsIn ∷ Int
    , minGraphCostOut ∷ VertexCost
    , maxGraphCostOut ∷ VertexCost
    , numGraphsOut ∷ Int
    , commentString ∷ String
    , duration ∷ Int
    }
    deriving stock (Show, Eq)


{- | maxSimultaneousGraphsSteepest is the maximum number of graphs that are evaluated
at a step in "steepest" algorithms of swap and fuse. Set becasue can increase
run time of these procedurs by delaying finding "better" solutins to move to.
-}
maxSimultaneousGraphsSteepest ∷ Int
maxSimultaneousGraphsSteepest = 10


-- | SwapType types for swapping, TBRAlternate for special casing in Swap
-- TBR only does non-Spr moves in a TBR searc process, would be done after SPR in an 
-- Alternate seaarch
data SwapType = NoSwap | NNI | SPR | TBR | Alternate | TBRAlternate |TBROnly
    deriving stock (Show, Eq)


-- | JoinType types for join methods
data JoinType = JoinPruned | JoinAll -- | JoinAlternate
    deriving stock (Show, Eq)


-- | SelectGraphType types to select gaphs
data SelectGraphType = Best | Unique | AtRandom | All
    deriving stock (Show, Eq)


-- | Support method types
data SupportMethod = Jackknife | Bootstrap | GoodmanBremer
    deriving stock (Show, Eq)


-- | return parallel Strategy
parStrategy ∷ (NFData b) ⇒ ParallelStrategy → Strategy b
parStrategy parStrat
    | parStrat == R0 = r0
    | parStrat == RPar = rpar
    | parStrat == RSeq = rseq
    | parStrat == RDeepSeq = rdeepseq
    | otherwise = rdeepseq


data GlobalSettings = GlobalSettings
    { bc2 ∷ (Double, Double) -- PMDL bitCost for 2 states of no-change and change as pair
    , bc4 ∷ (Double, Double) -- PMDL bitCost for 4 states of no-change and change as pair
    , bc5 ∷ (Double, Double) -- PMDL bitCost for 5 states of no-change and change as pair
    , bc8 ∷ (Double, Double) -- PMDL bitCost for 8 states of no-change and change as pair
    , bc64 ∷ (Double, Double) -- PMDL bitCost for 64 states of no-change and change as pair
    , bcgt64 ∷ (Double, Double) -- PMDL bitCost for > 64 states of no-change and change as pair
    , compressResolutions ∷ Bool -- "nub" resolutions in softwired graph
    , defaultParStrat ∷ ParallelStrategy -- default parallel strategy
    , dynamicEpsilon ∷ Double -- factor of dynamic heuristics overestimating graph deltas determind by fraction of data is dynamic and user value
    , finalAssignment ∷ AssignmentMethod
    , fractionDynamic ∷ Double -- estimated fraction of character length that are dynamic (actually seqeunce) for setting dynamicEpsilon
    , graphComplexityList ∷ IL.InfList (VertexCost, VertexCost) -- complexity of graphs in bits, index for number of network nodes (0= tree etc0 lazy so only evaluate each once when needed O(n) but needlazyness and permanence
    , graphsSteepest ∷ Int -- he maximum number of graphs that are evaluated
    -- at a step in "steepest" algorithms of swap and network add/delete. Set because can increase
    -- run time of these procedures by delaying finding "better" solutions to move to.
    -- also increases memory footprint
    , graphType ∷ GraphType
    , graphFactor ∷ GraphFactor -- net penalty/graph complexity
    , lazyParStrat ∷ ParallelStrategy -- default parallel strategy to WHNF
    , missingThreshold ∷ Int -- default threshold of maximum missing data to keep in data set 100 (keep all, 0 would be no missing data)
    , modelComplexity ∷ Double -- model cost for PMDL, 0.0 for other criteria
    , multiTraverseCharacters ∷ Bool -- If true "reroot" charcter trees to get best cost for (only affects) dynamic characters, if False then no
    , numDataLeaves ∷ Int -- number of leaves  set after data processing--for conveniance really
    , optimalityCriterion ∷ OptimalityCriterion
    , outgroupIndex ∷ Int -- Outgroup terminal index, default 0 (first input leaf)
    , outGroupName ∷ T.Text -- Outgroup name
    , partitionCharacter ∷ String -- 'character' for mparitioning seqeunce data into homologous sections'--checks for length == 1 later
    , reportNaiveData ∷ Bool -- reports using Naive data so preserves character order and codings.  This comes at a cost in memory footprint.  If False,
    -- packed characters are reported--and are somewhat inscrutable. But perhaps 3% of data footprint--useful for large
    -- add/non add dat asets liker SNP genomic data
    , rootComplexity ∷ VertexCost -- complexity of root in bits per root for PMDL/ML calculations
    , rootCost ∷ RootCost
    , searchData ∷ [SearchData]
    , seed ∷ Int -- random seed
    , softWiredMethod ∷ SoftWiredAlgorithm -- algorithm to optimize softwired graphs
    , strictParStrat ∷ ParallelStrategy -- default parallel strategy to Fully evaluate
    , unionThreshold ∷ Double -- this is the edge union cost threshold for rejoing edges during SPR and TBR, and (perhas) character Wagner build
    -- as described by Varon and Wheeler (2013) and set to 1.17 experimentally
    , useIA ∷ Bool -- turn on/off IA everywhere can (mainly for internal testing)
    , useNetAddHeuristic ∷ Bool -- Netowrk addition heuristic--very coarse currently
    }
    deriving stock (Show, Eq)


instance NFData GlobalSettings where rnf x = seq x ()


{- | CharInfo information about characters
null values for these are in Input.FastAC.hs
 TCMD.DenseTransitionCostMatrix          => genDiscreteDenseOfDimension (length alphabet)
 MR.MetricRepresentation WideState          => metricRepresentation <$> TCM.fromRows [[0::Word]]
 MR.MetricRepresentation BV.BitVector    => metricRepresentation <$> TCM.fromRows [[0::Word]]
changeCost and noChange costs are for PMDL costs for packed/non-additive character for
other character types the cost matrix holds this information in comvcert with weight since the matrix values are integers
-}
data CharInfo = CharInfo
    { name ∷ NameText
    , charType ∷ CharType
    , activity ∷ Bool
    , weight ∷ Double
    , costMatrix ∷ S.Matrix Int
    , slimTCM ∷ TCMD.DenseTransitionCostMatrix
    , wideTCM ∷ MR.MetricRepresentation WideState
    , hugeTCM ∷ MR.MetricRepresentation HugeState
    , changeCost ∷ Double
    , noChangeCost ∷ Double
    , alphabet ∷ Alphabet ST.ShortText
    , prealigned ∷ Bool
    , origInfo ∷ V.Vector (NameText, CharType, Alphabet ST.ShortText)
    }
    deriving stock (Show, Eq)


instance NFData CharInfo where rnf x = seq x ()


instance Ord CharInfo where x `compare` y = show x `compare` show y


-- | Types for vertex information
type VertexCost = Double


type StateCost = Int


type VertexIndex = Int


-- | index of child vertices
type ChildStateIndex = Int


{- | unique bitvector labelling for vertex based on descednent labellings
these labels are used for caching, left/right DO optimizaiton
thery should be label invariant
a hash of sorted input data for leaves
will need a map from NameBV to T.Text name (or not)
-}
type NameBV = BV.BitVector


-- | Human legibale name for vertices, characters, and Blocks
type NameText = T.Text


{- | TYpes for Matrix/Sankoff characters
Triple contains info from left and right child--could be only one
use fst then
finals state also vector of triple--but need to keep min cost
for final assignments
filter by local costVect to get 'best' states an each node
-}
type MatrixTriple = (StateCost, [ChildStateIndex], [ChildStateIndex])


-- Only date here that varies by vertex, rest inglobal character info
-- vectors so all data of single type can be grouped together
-- will need to add masks for bit-packing non-additive chars
-- may have to add single assignment for hardwired and IP optimization
-- for approximate sakoff (DO-like) costs can use stateBVPrelim/stateBVFinal
-- for matrix/Saknoff characters-- Vector of vector of States
-- BUT all with same cost matrix/tcm
-- triples (add, no-add, sequence) are to keep children of vertex states for pre-order pass
-- order is always (left parent median, median, right parent median)
-- do not need for matrix since up pass is a traceback from parent
-- sequence characters are a vector of bitvectors--so only a single seqeunce character
--  per "charctaer" this is so the multi-traversal can take place independently for each
--  sequence character, creating a properly "rooted" tree/graph for each non-exact seqeunce character
-- prelim is created from gapped, final from (via 3-way minimization) parent final and child alignment (2nd and 3rd fields).
-- the 'alignment' fields hold the implied alignment data
-- the 'union' fields hold post-order unions of subgraph charcaters (ia for sequence) fir use in uinion threshold
-- during branch addintion/readdition (e.g swapping)
data CharacterData = CharacterData
    { -- for Non-additive
      stateBVPrelim ∷ (V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector) -- HugeDynamicCharacter -- preliminary for Non-additive chars, Sankoff Approx
    , stateBVFinal ∷ V.Vector BV.BitVector
    , stateBVUnion ∷ V.Vector BV.BitVector
    , -- for Additive
      rangePrelim ∷ (V.Vector (Int, Int), V.Vector (Int, Int), V.Vector (Int, Int))
    , rangeFinal ∷ V.Vector (Int, Int)
    , rangeUnion ∷ V.Vector (Int, Int)
    , -- for multiple Sankoff/Matrix with slim tcm
      matrixStatesPrelim ∷ V.Vector (V.Vector MatrixTriple)
    , matrixStatesFinal ∷ V.Vector (V.Vector MatrixTriple)
    , matrixStatesUnion ∷ V.Vector (V.Vector MatrixTriple)
    , -- preliminary for m,ultiple seqeunce chars with same TCM
      slimPrelim ∷ SV.Vector SlimState
    , -- gapped medians of left, right, and preliminary used in preorder pass
      slimGapped ∷ SlimDynamicCharacter
    , slimAlignment ∷ SlimDynamicCharacter
    , slimFinal ∷ SV.Vector SlimState
    , slimIAPrelim ∷ SlimDynamicCharacter
    , slimIAFinal ∷ SV.Vector SlimState
    , slimIAUnion ∷ SV.Vector SlimState
    , -- vector of individual character costs (Can be used in reweighting-ratchet)
      widePrelim ∷ UV.Vector WideState
    , -- gapped median of left, right, and preliminary used in preorder pass
      wideGapped ∷ WideDynamicCharacter
    , wideAlignment ∷ WideDynamicCharacter
    , wideFinal ∷ UV.Vector WideState
    , wideIAPrelim ∷ WideDynamicCharacter
    , wideIAFinal ∷ UV.Vector WideState
    , wideIAUnion ∷ UV.Vector WideState
    , -- vector of individual character costs (Can be used in reweighting-ratchet)
      hugePrelim ∷ V.Vector HugeState
    , -- gapped medians of left, right, and preliminary used in preorder pass
      hugeGapped ∷ HugeDynamicCharacter
    , hugeAlignment ∷ HugeDynamicCharacter
    , hugeFinal ∷ V.Vector HugeState
    , hugeIAPrelim ∷ HugeDynamicCharacter
    , hugeIAFinal ∷ V.Vector HugeState
    , hugeIAUnion ∷ V.Vector HugeState
    , -- vectors for pre-aligned sequences also used in static approx
      alignedSlimPrelim ∷ SlimDynamicCharacter
    , alignedSlimFinal ∷ SV.Vector SlimState
    , alignedSlimUnion ∷ SV.Vector SlimState
    , alignedWidePrelim ∷ WideDynamicCharacter
    , alignedWideFinal ∷ UV.Vector WideState
    , alignedWideUnion ∷ UV.Vector WideState
    , alignedHugePrelim ∷ HugeDynamicCharacter
    , alignedHugeFinal ∷ V.Vector HugeState
    , alignedHugeUnion ∷ V.Vector HugeState
    , -- coiuld be made Storable later is using C or GPU/Accelerate
      packedNonAddPrelim ∷ OpenDynamicCharacter UV.Vector Word64
    , packedNonAddFinal ∷ UV.Vector Word64
    , packedNonAddUnion ∷ UV.Vector Word64
    , -- vector of individual character costs (Can be used in reweighting-ratchet)
      localCostVect ∷ V.Vector StateCost
    , -- weight * V.sum localCostVect
      localCost ∷ VertexCost
    , -- unclear if need vector version
      globalCost ∷ VertexCost
    }
    deriving stock (Show, Eq, Generic)


instance NFData CharacterData where rnf x = seq x ()


{- | type TermData type contains termnal name and list of characters
characters as ShortText to save space on input
-}
type TermData = (NameText, [ST.ShortText])


type LeafData = (NameText, V.Vector CharacterData)


{- | VertexBlockData vector over blocks of character data in block (Vector)
blocks of character data for  a given vertex
-}
type VertexBlockData = V.Vector (V.Vector CharacterData)


{- | VertexBlockDataMaybe vector over maybe  blocks of character data in block (Vector)
blocks of character data for  a given vertex
-}
type VertexBlockDataMaybe = V.Vector (V.Vector (Maybe CharacterData))


{- | ResolutionData contains vertex information for soft-wired network components
these are used in the idenitification of minimal cost display trees for a block of
data that follow the same display tree
-}
type ResolutionVertexData = V.Vector ResolutionBlockData


{- | ResolutionBlockData contains a list of ResolutionData
this list contains all the potential resolutions of a softwired
networtk vertex
-}
type ResolutionBlockData = V.Vector ResolutionData


{- | ResolutionData contains individual block information for a given resoluton of soft-wired network components
these are used in the idenitification of minimal cost display trees for a block of
data that follow the same display tree
nodes are VertexInfo for ease of conversion--but nthe info is largely bogus and not to be trusted, same with EdgeInfo
-}
data ResolutionData = ResolutionData
    { displaySubGraph ∷ ([LG.LNode VertexInfo], [LG.LEdge EdgeInfo]) -- holds the post-order display sub-tree for the block
    , displayBVLabel ∷ NameBV -- For comparison of vertices subtrees, left/right, anmd root leaf inclusion
    , displayData ∷ V.Vector CharacterData -- data for characters in block
    -- left and right indices of child resolution Data for traceback and preliminary state assignment
    , childResolutionIndices ∷ (Maybe Int, Maybe Int)
    , resolutionCost ∷ VertexCost -- cost of creating the resolution
    , displayCost ∷ VertexCost -- cost of that display subtree
    }
    deriving stock (Show, Eq)


instance NFData ResolutionData where rnf x = seq x ()


-- | VertexInfo type -- vertex information for Decorated Graph
data VertexInfo = VertexInfo
    { bvLabel ∷ NameBV -- For comparison of vertices subtrees, left/right
    , children ∷ V.Vector Int -- outdegree indices
    , index ∷ Int -- For accessing
    , nodeType ∷ NodeType -- root, leaf, network, tree
    , parents ∷ V.Vector Int -- indegree indices
    , subGraphCost ∷ VertexCost -- cost of graph to leaves from the vertex
    , vertexCost ∷ VertexCost -- local cost of vertex
    , vertData ∷ VertexBlockData -- data as vector of blocks (each a vector of characters)
    , vertName ∷ NameText -- Text name of vertex either input or HTU#
    , vertexResolutionData ∷ V.Vector ResolutionBlockData -- soft-wired network component resolution information for Blocks
    }
    deriving stock (Generic, Show, Eq)


instance NFData VertexInfo where rnf x = seq x ()


-- | type edge data, source and sink node indices are fst3 and snd3 fields.
data EdgeInfo = EdgeInfo
    { edgeType ∷ EdgeType
    , maxLength ∷ VertexCost
    , midRangeLength ∷ VertexCost
    , minLength ∷ VertexCost
    }
    deriving stock (Show, Eq, Ord)


instance NFData EdgeInfo where rnf x = seq x ()


{- | DecortatedGraph is the canonical graph contining all final information
from preorder traversal trees
and post-order info usually from an initial root-based traversal
-}
type DecoratedGraph = LG.Gr VertexInfo EdgeInfo


{- | Type BLockDisplayTree is a Forest of tree components (indegree, outdegree) = (0,1|2),(1,2),(1,0)
these are "resolved" from more general graphs
will have to allow for indegre=outdegree=1 for dispaly tree generation and reconciliation
the vertData field will always have a single Bloclk--teh vecor of blocks will be a vector of
DecoratedGraphs.  These woulod better have a single Vector of cChracter info as
opposed to the Decorated Tree type, but the record naming and use gets screwed up.
type BlockDisplayForest = LG.Gr VertexInfo EdgeInfo
-}

{- | DecoratedGraph is a forest of tree compnents for a single character
this is used for non-exact character traversal trees
there will always be only a single block and single character even though
expresed as Vector fo Vector of Chaarcters.  Would be better as a single character as
opposed to the Decorated Tree type, but the record naming and use gets screwed up.
type DecoratedGraph = LG.Gr VertexInfo EdgeInfo
-}

-- | type RawGraph is input graphs with leaf and edge labels
type SimpleGraph = LG.Gr NameText VertexCost


{- | Type phylogentic Graph is a graph with
cost, optimality value,
block display trees, character traversal "foci" (could have multiple)
Data optimizations exist in Processed Data
Question of wheterh traversal foci shold be in Graph or Data section
for now in Graph but could easiily be moved to Processed data
May need "EdgeData" like VertexData for heuristics.  UNclear if local scope during SPR/TBR will do it.
   Fields:
       1) "Simple" graph with fileds useful for outputting graphs
       2) Graph optimality value or cost
       3) Decorated Graph with optimized vertex/Node data
       4) Vector of display trees for each data Block
                 root and vertex costs not updated in rerooting so cannot be trusted
                 Each block can have multiple display trees so extra list there
       5) Vector of traversal foci for each character (Vector of Blocks -> Vector of Characters, a single tree for each character)
              vector is over blocks, then characters (could have have multiple for each character, but only single tracked here)
              only important for dynamic (ie non-exact) characters whose costs depend on traversal focus
              one graph per character
       6) Vector of Block Character Information (whihc is a Vector itself) required to properly optimize characters
-}
type PhylogeneticGraph =
    ( SimpleGraph
    , VertexCost
    , DecoratedGraph
    , V.Vector [DecoratedGraph]
    , V.Vector (V.Vector DecoratedGraph)
    , V.Vector (V.Vector CharInfo)
    )


-- | Type phylogentic Graph is a graph with general types
type GenPhyloGraph a b =
    (SimpleGraph, VertexCost, DecoratedGraph, V.Vector [LG.Gr a b], V.Vector (V.Vector (LG.Gr a b)), V.Vector (V.Vector CharInfo))


{- | Type ReducedPhylogenticGraph is a graph
that has most of the information of a PhylogeneticGraph but does not have repeated
decorations.  The only lacking information is the traversal topologies of the charcter graphs (5th field of phylogenetic graph)
the display tree field (4th) does not have decorations but only toplogy infomatin in form of SimpleGraph
the purpose of the type is to remove the redundant decorations to 1/3 of what they are in PhylogeneticGraph
   Fields:
       1) "Simple" graph with fileds useful for outputting graphs
       2) Graph optimality value or cost
       3) Decorated Graph with optimized vertex/Node data
       4) Vector of display trees for each data Block as Simple Graphs
                 Each block can have multiple display trees so extra list there
       5) Vector of Block Character Information (whihc is a Vector itself) required to properly optimize characters
-}
type ReducedPhylogeneticGraph = (SimpleGraph, VertexCost, DecoratedGraph, V.Vector [SimpleGraph], V.Vector (V.Vector CharInfo))


{- | RawData type processed from input to be passed to characterData
to recode into usable form
the format is tuple of a list of taxon-data list tuples and charinfo list.
the data list and charinfo list must have the same length
-}
type RawData = ([TermData], [CharInfo])


{- | Processed data is the basic data structire for analysis
can be input to functions
based on "blocks" that follow same display tree (soft-wired network)
each block has a set of characters (for each vertex eventually) and character info
the vector of T.Text are the names--leaves form input, internal HTU <> (show index)
ablock are initialy set bu input file, and can be later changed by "set(block,...)"
command
"Naive" "Optimized" and "Transformed" darta are this type after different processing steps
the first and second vectors are size number of leaves, teh third is number of blocks
-}
type ProcessedData = (V.Vector NameText, V.Vector NameBV, V.Vector BlockData)


{- | Block data  is the basic data unit that is optimized on a display tree
 Block data contain data fo all leaves and all characters in the block
it is row, ie vertex dominant
it has a bitvector name derived from leaf bitvector labels (union of children)
the bitvector name can vary among blocks (for non-leaves) due to alternate display trees
a vector of characterer data where assignments and costs reside, and a vector of character info
leaves will alwasy be first (indices 0..n-1) for simpler updating of data during graph optimization
NameText is the block label used for assignment and reporting output
Initially set to input filename of character
   Fields:
       1) name of the block--intially taken from input filenames
       2) vector of vertex/leaf data with vector of character data for each leaf
       3) Vector of character information for characters in the block
-}
type BlockData = (NameText, V.Vector (V.Vector CharacterData), V.Vector CharInfo)


-- | type SAParams parameter structure for simulated alnnealing and Drifting
data SimulatedAnnealingMethod = SimAnneal | Drift
    deriving stock (Read, Show, Eq)


-- | Simulated Annealing parameters
data SAParams = SAParams
    { method ∷ SimulatedAnnealingMethod
    , numberSteps ∷ Int
    , currentStep ∷ Int
    , rounds ∷ Int
    , driftAcceptEqual ∷ Double
    , driftAcceptWorse ∷ Double
    , driftMaxChanges ∷ Int
    , driftChanges ∷ Int
    }
    deriving stock (Show, Eq)


instance NFData SAParams where rnf x = seq x ()

-- type for reoptimizinmg candiate graphs after heuristic cost
    -- Best = only the best/lowest heuristic costs get rechecked
    -- Better = all thse graphs with better heuristic scores than the curernt best score
    -- BestN = check the best N scores that are better than the curent best score
    -- BestAll checks all graphs

data HeuristicCheck = BestOnly | Better | BetterN | BestAll
    deriving stock (Read, Show, Eq)

-- | SwapParam type for swap parameers
data SwapParams = SwapParams
    { atRandom ∷ Bool -- randomized splitting and rejoining
    , checkHeuristic :: HeuristicCheck -- for reoptimizing graphs after heuristic costs
    , doIA ∷ Bool -- use Implied alignment fields for rearragement costs
    , joinAlternate ∷ Bool -- in alternate swapping for TBR
    , joinType ∷ JoinType -- Union pruning on or off
    , keepNum ∷ Int -- number equally costly solutoins to keep
    , maxMoveEdgeDist ∷ Int -- maximum rejoin distance from initial mplacement
    , returnMutated ∷ Bool -- return changed graphs for simlated annealing, genetic algorithm
    , sortEdgesSplitCost :: Bool -- sort edges based on split cost-- greatest delta first
    , splitParallel :: Bool -- when splittting graph--do spliots in parallel or sequenctial
    , steepest ∷ Bool -- steepest descent versus "all"
    , swapType ∷ SwapType -- NNI/SPR/TBR/Alternate
    }
    deriving stock (Show, Eq)


instance NFData SwapParams where rnf x = seq x ()


-- | empty structures for convenient use

-- | emptyProcessedData empty processsed data dfor memory saving with large qualitative data sets(e.g. SNPS)
emptyProcessedData ∷ ProcessedData
emptyProcessedData = (V.empty, V.empty, V.empty)


-- | emptySearchData for use in getting basic procesin input data
emptySearchData ∷ SearchData
emptySearchData =
    SearchData
        { instruction = NotACommand
        , arguments = []
        , minGraphCostIn = infinity
        , maxGraphCostIn = infinity
        , numGraphsIn = 0
        , minGraphCostOut = infinity
        , maxGraphCostOut = infinity
        , numGraphsOut = 0
        , commentString = []
        , duration = 0
        }


-- | emptyGlobalSettings for early use in parition charcter.  Can't have full due to data dependency of outpogroup name
emptyGlobalSettings ∷ GlobalSettings
emptyGlobalSettings =
    GlobalSettings
        { outgroupIndex = 0
        , outGroupName = "NoOutgroupSet"
        , optimalityCriterion = Parsimony
        , graphType = Tree
        , compressResolutions = False
        , finalAssignment = DirectOptimization
        , graphFactor = Wheeler2015Network
        , rootCost = NoRootCost
        , rootComplexity = 0.0
        , graphComplexityList = IL.repeat (0.0, 0.0)
        , modelComplexity = 0.0
        , seed = 0
        , searchData = []
        , partitionCharacter = "#"
        , numDataLeaves = 0
        , bc2 = (0.0, 1.0)
        , bc4 = (0.0, 1.0)
        , bc5 = (0.0, 1.0)
        , bc8 = (0.0, 1.0)
        , bc64 = (0.0, 1.0)
        , bcgt64 = (0.0, 1.0)
        , fractionDynamic = 1.0
        , dynamicEpsilon = 1.00
        , graphsSteepest = 10
        , softWiredMethod = ResolutionCache
        , multiTraverseCharacters = True
        , reportNaiveData = True
        , unionThreshold = 1.17
        , defaultParStrat = RSeq
        , lazyParStrat = RPar -- default parallel strategy
        , strictParStrat = RDeepSeq -- high level--basically srtict evaluation
        , useNetAddHeuristic = True
        , useIA = True
        , missingThreshold = 100
        }


{- | emptyPhylogeneticGraph specifies and empty phylogenetic graph
important cost is infinity for filtering operations
-}
emptyPhylogeneticGraph ∷ PhylogeneticGraph
emptyPhylogeneticGraph = (LG.empty, infinity, LG.empty, V.empty, V.empty, V.empty)


{- | emptyReducedPhylogeneticGraph specifies and empty phylogenetic graph
important cost is infinity for filtering operations
-}
emptyReducedPhylogeneticGraph ∷ ReducedPhylogeneticGraph
emptyReducedPhylogeneticGraph = (LG.empty, infinity, LG.empty, V.empty, V.empty)


-- | emptycharacter useful for intialization and missing data
emptyCharacter ∷ CharacterData
emptyCharacter =
    CharacterData -- for non-additive
        { stateBVPrelim = (mempty, mempty, mempty) -- preliminary for Non-additive chars, Sankoff Approx
        , stateBVFinal = mempty
        , stateBVUnion = mempty
        , -- for Additive
          rangePrelim = (mempty, mempty, mempty)
        , rangeFinal = mempty
        , rangeUnion = mempty
        , -- for multiple Sankoff/Matrix with sme tcm
          matrixStatesPrelim = mempty
        , matrixStatesFinal = mempty
        , matrixStatesUnion = mempty
        , -- preliminary for m,ultiple seqeunce cahrs with same TCM
          slimPrelim = mempty
        , -- gapped medians of left, right, and preliminary used in preorder pass
          slimGapped = (mempty, mempty, mempty)
        , slimAlignment = (mempty, mempty, mempty)
        , slimFinal = mempty
        , slimIAPrelim = (mempty, mempty, mempty)
        , slimIAFinal = mempty
        , slimIAUnion = mempty
        , -- gapped median of left, right, and preliminary used in preorder pass
          widePrelim = mempty
        , -- gapped median of left, right, and preliminary used in preorder pass
          wideGapped = (mempty, mempty, mempty)
        , wideAlignment = (mempty, mempty, mempty)
        , wideFinal = mempty
        , wideIAPrelim = (mempty, mempty, mempty)
        , wideIAFinal = mempty
        , wideIAUnion = mempty
        , -- vector of individual character costs (Can be used in reweighting-ratchet)
          hugePrelim = mempty
        , -- gapped mediasn of left, right, and preliminary used in preorder pass
          hugeGapped = (mempty, mempty, mempty)
        , hugeAlignment = (mempty, mempty, mempty)
        , hugeFinal = mempty
        , hugeIAPrelim = (mempty, mempty, mempty)
        , hugeIAFinal = mempty
        , hugeIAUnion = mempty
        , -- vectors for pre-aligned sequences also used in static approx
          alignedSlimPrelim = (mempty, mempty, mempty)
        , alignedSlimFinal = mempty
        , alignedSlimUnion = mempty
        , alignedWidePrelim = (mempty, mempty, mempty)
        , alignedWideFinal = mempty
        , alignedWideUnion = mempty
        , alignedHugePrelim = (mempty, mempty, mempty)
        , alignedHugeFinal = mempty
        , alignedHugeUnion = mempty
        , packedNonAddPrelim = (mempty, mempty, mempty)
        , packedNonAddFinal = mempty
        , packedNonAddUnion = mempty
        , -- vector of individual character costs (Can be used in reweighting-ratchet)
          localCostVect = V.singleton 0
        , -- weight * V.sum localCostVect
          localCost = 0
        , -- unclear if need vector version
          globalCost = 0
        }


-- | emptyVertex useful for graph rearrangements
emptyVertexInfo ∷ VertexInfo
emptyVertexInfo =
    VertexInfo
        { index = -1
        , bvLabel = BV.fromBits [False]
        , parents = mempty
        , children = mempty
        , nodeType = TreeNode -- root, leaf, network, tree
        , vertName = "EmptyVertex"
        , vertData = mempty
        , vertexResolutionData = mempty
        , vertexCost = 0.0
        , subGraphCost = 0.0
        }


-- | usefule in some cases
dummyNode ∷ LG.LNode VertexInfo
dummyNode = (-1, emptyVertexInfo)


-- | dummyEdge for convenience
dummyEdge ∷ EdgeInfo
dummyEdge =
    EdgeInfo
        { minLength = 0
        , maxLength = 0
        , midRangeLength = 0
        , edgeType = TreeEdge
        }


-- emptyCharInfo for convenience
emptyCharInfo ∷ (MonadIO m) ⇒ m CharInfo
emptyCharInfo =
    let minimalMatrix ∷ [[Int]]
        minimalMatrix = [[0, 1], [1, 0]]
        (_, tcm) = TCM.fromRows minimalMatrix
        sTCM = TCMD.generateDenseTransitionCostMatrix 2 2 . S.getCost $ V.fromList <$> V.fromList minimalMatrix
    in  do
            wTCM ← MR.metricRepresentation tcm
            hTCM ← MR.metricRepresentation tcm
            pure
                CharInfo
                    { name = "EmptyCharName"
                    , charType = NonAdd
                    , activity = True
                    , weight = 1.0
                    , costMatrix = S.empty
                    , slimTCM = sTCM
                    , wideTCM = wTCM
                    , hugeTCM = hTCM
                    , changeCost = 1.0
                    , noChangeCost = 0.0
                    , alphabet = fromSymbols $ "0" :| ["1"]
                    , prealigned = False
                    , origInfo = V.empty
                    }
