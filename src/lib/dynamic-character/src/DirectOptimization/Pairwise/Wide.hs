module DirectOptimization.Pairwise.Wide
  ( WideDynamicCharacter
  , WideState
  , widePairwiseDO
  ) where

import Bio.DynamicCharacter
import DirectOptimization.Pairwise.Internal (AlignmentCost, SymbolChangeCost, TCM2Dλ)
import DirectOptimization.Pairwise.Ukkonen


widePairwiseDO
  :: SymbolChangeCost
  -> TCM2Dλ WideState
  -> WideDynamicCharacter
  -> WideDynamicCharacter
  -> (AlignmentCost, WideDynamicCharacter)
widePairwiseDO = ukkonenDO
