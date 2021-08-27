module Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Wide
  ( WideDynamicCharacter
  , WideState
  , widePairwiseDO
  ) where

import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.DynamicCharacter2
import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Ukkonen


widePairwiseDO
  :: Word
  -> WideState
  -> (WideState -> WideState -> (WideState, Word))
  -> WideDynamicCharacter
  -> WideDynamicCharacter
  -> (Word, WideDynamicCharacter)
widePairwiseDO = ukkonenDO
