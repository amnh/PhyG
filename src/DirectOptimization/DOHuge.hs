module DirectOptimization.DOHuge where

import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Huge
import Data.BitVector.LittleEndian                                 (BitVector)
import Data.MetricRepresentation
import Data.Vector                                                 (Vector)


wrapperHugeDO
  :: Vector BitVector
  -> Vector BitVector
  -> MetricRepresentation BitVector
  -> (Vector BitVector, Int)
wrapperHugeDO lhs rhs tcmMemo = (medians, fromEnum cost)
  where
    (cost, (medians,_,_)) = hugePairwiseDO (retreivePairwiseTCM tcmMemo) (lhs, lhs, lhs) (rhs, rhs, rhs)
