module DirectOptimization.Pairwise.Slim
  ( SlimDynamicCharacter
  , SlimState
  , slimPairwiseDO
  ) where

import Bio.DynamicCharacter
import DirectOptimization.Pairwise.Internal (Distance)
import DirectOptimization.Pairwise.Slim.FFI


slimPairwiseDO
  :: TCMρ
  -> SlimDynamicCharacter
  -> SlimDynamicCharacter
  -> (Distance, SlimDynamicCharacter)
slimPairwiseDO = smallAlphabetPairwiseDO
