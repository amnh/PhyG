module DirectOptimization.DOHuge where

import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Huge
import Data.BitVector.LittleEndian                                 (BitVector, dimension)
import Data.Bits
import Data.Foldable
import Data.MetricRepresentation
import Data.Vector                                                 (Vector, head)
import Prelude                                                     hiding (head)

wrapperHugeDO
  :: Vector BitVector
  -> Vector BitVector
  -> MetricRepresentation BitVector
  -> (Vector BitVector, Int)
wrapperHugeDO lhs rhs tcmMemo = (medians, fromEnum cost)
  where
    (cost, (medians,_,_)) = hugePairwiseDO (minInDelCost tcmMemo) gapState (retreivePairwiseTCM tcmMemo) (lhs, lhs, lhs) (rhs, rhs, rhs)

    gapState = bit . fromEnum $ n - 1
    n = case length lhs of
          0 -> case length rhs of
                 0 -> 64
                 _ -> dimension $ head rhs
          _ -> dimension $ head lhs
