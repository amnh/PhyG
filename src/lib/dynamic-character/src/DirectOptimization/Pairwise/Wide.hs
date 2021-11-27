module DirectOptimization.Pairwise.Wide
  ( WideDynamicCharacter
  , WideState
  , widePairwiseDO
  ) where

import Bio.DynamicCharacter
import DirectOptimization.Pairwise.Internal (Distance, TCM2Dλ)
import DirectOptimization.Pairwise.Ukkonen


widePairwiseDO
  :: Distance
  -> TCM2Dλ WideState
  -> WideDynamicCharacter
  -> WideDynamicCharacter
  -> (Distance, WideDynamicCharacter)
widePairwiseDO = ukkonenDO
