{- |
Module specifying graph swapping rearrangement functions
-}
module Search.SwapV2 (
    swapV2,
) where

import Control.Monad (filterM)
import Control.Monad.IO.Class
import Control.Monad.Random.Class
import Data.Foldable (fold, toList)
import Data.Foldable1 (Foldable1)
import Data.Foldable1 qualified as F1
import Data.Functor (($>), (<&>))
import Data.List qualified as L
import Data.Maybe
import Data.Ord (comparing)
import Data.Vector qualified as V
import GHC.Real qualified as Real (infinity)
import GeneralUtilities
import GraphOptimization.Medians qualified as M
import GraphOptimization.PostOrderSoftWiredFunctions qualified as POSW
import GraphOptimization.PreOrderFunctions qualified as PRE
import GraphOptimization.Traversals qualified as T
import Graphs.GraphOperations qualified as GO
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import PHANE.Evaluation.Verbosity (Verbosity (..))
import Search.Swap qualified as SV1
import Types.Types
import Utilities.LocalGraph qualified as LG
import Utilities.Utilities as U


{- | New Swap functions that are based on PHANE paralleization routines.
    1) Naive version
    2) steepest descent
        Order edges by distace to original break point on rejoin
        Continue when swtich to new (order of edges)
        Keep single graph only until no better, then multiple
    3) Randomized trajectories
    4) heuristic cost calculations
    5) SA/Drift
    6) unions
    
-}
swapV2
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Int
    → [ReducedPhylogeneticGraph]
    → [(Maybe SAParams, ReducedPhylogeneticGraph)]
    → PhyG ([ReducedPhylogeneticGraph], Int)
swapV2 swapParams inGS inData inCounter curBestGraphList inSimAnnealParams =
    if null curBestGraphList
        then pure ([], inCounter)
        else swapNaive swapParams inGS inData inCounter curBestGraphList inSimAnnealParams


{- | Naive swap functions to create reference for later algorithmic improvements 
    1) Take first graph

    2) Create list of splits of graph 
        LG.getEdgeSplitList

    3) Rejoin all places for all splits
        This is "all around" in that not switching to lower cost graphs
        at first opportunity (ie. 'steepest").
        LG.splitGraphOnEdge (check base graph leaves >=3)
        LG.joinGraphOnEdge

    4) Full evaluation of graphs 
        Time complexities should be with O(n) for graph traversals
            NNI O(n^2)
            SPR O(n^3)
            TBR O(n^4)
-}
swapNaive 
    ∷ SwapParams
    → GlobalSettings
    → ProcessedData
    → Int
    → [ReducedPhylogeneticGraph]
    → [(Maybe SAParams, ReducedPhylogeneticGraph)]
    → PhyG ([ReducedPhylogeneticGraph], Int)
swapNaive swapParams inGS inData inCounter curBestGraphList inSimAnnealParams =
    if null curBestGraphList
        then pure ([], inCounter)
        else pure (curBestGraphList, inCounter)


