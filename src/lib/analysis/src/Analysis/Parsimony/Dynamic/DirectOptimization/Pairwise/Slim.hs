module Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Slim
  ( SlimDynamicCharacter
  , slimPairwiseDO
  ) where

import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.DynamicCharacter2
import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Slim.FFI

slimPairwiseDO
  :: DenseTransitionCostMatrix
  -> SlimDynamicCharacter
  -> SlimDynamicCharacter
  -> (Word, SlimDynamicCharacter)
slimPairwiseDO = smallAlphabetPairwiseDO
