module Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Huge
  ( HugeDynamicCharacter
  , hugePairwiseDO
  ) where

import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.DynamicCharacter2
import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Huge.Ukkonen
import Data.BitVector.LittleEndian
import Data.Word


hugePairwiseDO
  :: (BitVector -> BitVector -> (BitVector, Word))
  -> HugeDynamicCharacter
  -> HugeDynamicCharacter
  -> (Word, HugeDynamicCharacter)
hugePairwiseDO = hugeUkkonenDO
