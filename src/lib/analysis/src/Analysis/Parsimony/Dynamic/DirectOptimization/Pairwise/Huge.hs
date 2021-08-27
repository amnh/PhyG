module Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Huge
  ( HugeDynamicCharacter
  , HugeState
  , hugePairwiseDO
  ) where

import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.DynamicCharacter2
import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Ukkonen


hugePairwiseDO
  :: Word
  -> HugeState
  -> (HugeState -> HugeState -> (HugeState, Word))
  -> HugeDynamicCharacter
  -> HugeDynamicCharacter
  -> (Word, HugeDynamicCharacter)
hugePairwiseDO = ukkonenDO
