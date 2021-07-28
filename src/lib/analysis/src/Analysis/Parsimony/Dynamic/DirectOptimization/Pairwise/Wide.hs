module Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Wide
  ( WideDynamicCharacter
  , widePairwiseDO
  ) where

import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.DynamicCharacter2
import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Wide.Ukkonen
import Data.Word


widePairwiseDO
  :: Word8
  -> (Word64 -> Word64 -> (Word64, Word))
  -> WideDynamicCharacter
  -> WideDynamicCharacter
  -> (Word, WideDynamicCharacter)
widePairwiseDO = wideUkkonenDO
