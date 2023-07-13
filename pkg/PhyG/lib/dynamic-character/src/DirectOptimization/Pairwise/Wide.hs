module DirectOptimization.Pairwise.Wide
  ( WideDynamicCharacter
  , WideState
  , widePairwiseDO
  ) where

import Bio.DynamicCharacter
import DirectOptimization.Pairwise.Ukkonen


widePairwiseDO
  :: Word
  -> (WideState -> WideState -> (WideState, Word))
  -> WideDynamicCharacter
  -> WideDynamicCharacter
  -> (Word, WideDynamicCharacter)
widePairwiseDO = ukkonenDO
