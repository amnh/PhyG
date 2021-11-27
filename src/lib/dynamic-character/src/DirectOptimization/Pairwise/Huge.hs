module DirectOptimization.Pairwise.Huge
  ( HugeDynamicCharacter
  , HugeState
  , hugePairwiseDO
  ) where

import Bio.DynamicCharacter
import DirectOptimization.Pairwise.Internal (Distance, TCM2Dλ)
import DirectOptimization.Pairwise.Ukkonen


hugePairwiseDO
  :: Distance
  -> TCM2Dλ HugeState
  -> HugeDynamicCharacter
  -> HugeDynamicCharacter
  -> (Distance, HugeDynamicCharacter)
hugePairwiseDO = ukkonenDO
