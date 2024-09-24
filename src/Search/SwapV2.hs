{- |
Module specifying graph swapping rearrangement functions
-}
module Search.SwapV2 (
    swapV2,
    swapNaive,
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