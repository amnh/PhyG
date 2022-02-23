module DirectOptimization.Pairwise.Huge
  ( HugeDynamicCharacter
  , HugeState
  , hugePairwiseDO
  ) where

import Bio.DynamicCharacter
import DirectOptimization.Pairwise.Internal (AlignmentCost, SymbolChangeCost, TCM2Dλ)
import DirectOptimization.Pairwise.Ukkonen


hugePairwiseDO
  :: SymbolChangeCost
  -> TCM2Dλ HugeState
  -> HugeDynamicCharacter
  -> HugeDynamicCharacter
  -> (AlignmentCost, HugeDynamicCharacter)
hugePairwiseDO = ukkonenDO
