module DirectOptimization.Pairwise.Slim
  ( SlimDynamicCharacter
  , SlimState
  , slimPairwiseDO
  ) where

import Bio.DynamicCharacter
import DirectOptimization.Pairwise.Internal (AlignmentCost)
import DirectOptimization.Pairwise.Slim.FFI


slimPairwiseDO
  :: TCMρ
  -> SlimDynamicCharacter
  -> SlimDynamicCharacter
  -> (AlignmentCost, SlimDynamicCharacter)
slimPairwiseDO = smallAlphabetPairwiseDO
