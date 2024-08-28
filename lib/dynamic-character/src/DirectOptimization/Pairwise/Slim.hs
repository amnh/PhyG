module DirectOptimization.Pairwise.Slim (
    SlimDynamicCharacter,
    SlimState,
    slimPairwiseDO,
) where

import Bio.DynamicCharacter
import Bio.DynamicCharacter.Element (SlimState)
import DirectOptimization.Pairwise.Slim.FFI


slimPairwiseDO
    ∷ DenseTransitionCostMatrix
    → SlimDynamicCharacter
    → SlimDynamicCharacter
    → (Word, SlimDynamicCharacter)
slimPairwiseDO = smallAlphabetPairwiseDO
