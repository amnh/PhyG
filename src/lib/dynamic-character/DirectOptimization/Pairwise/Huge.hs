module DirectOptimization.Pairwise.Huge
  ( HugeDynamicCharacter
  , HugeState
  , hugePairwiseDO
  ) where

import Bio.DynamicCharacter
import DirectOptimization.Pairwise.Ukkonen


hugePairwiseDO
  :: Word
  -> HugeState
  -> (HugeState -> HugeState -> (HugeState, Word))
  -> HugeDynamicCharacter
  -> HugeDynamicCharacter
  -> (Word, HugeDynamicCharacter)
hugePairwiseDO = ukkonenDO
