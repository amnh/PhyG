{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE Strict #-}

module DirectOptimization.PreOrder (
    preOrderLogic,
) where

import Bio.DynamicCharacter
import Control.Monad
import Data.Bits
import Data.STRef
import Data.Vector.Generic (Vector)


{- |
Faithful translation of Algorithm 8 (Non-root Node Alignment) from the
"Efficient Implied Alignment" paper found at:
<https://doi.org/10.1186/s12859-020-03595-2 10.1186/s12859-020-03595-2 DOI>
-}
{-# INLINEABLE preOrderLogic #-}
{-# SPECIALIZE preOrderLogic ∷
    Bool → SlimDynamicCharacter → SlimDynamicCharacter → SlimDynamicCharacter → SlimDynamicCharacter
    #-}
{-# SPECIALIZE preOrderLogic ∷
    Bool → WideDynamicCharacter → WideDynamicCharacter → WideDynamicCharacter → WideDynamicCharacter
    #-}
{-# SPECIALIZE preOrderLogic ∷
    Bool → HugeDynamicCharacter → HugeDynamicCharacter → HugeDynamicCharacter → HugeDynamicCharacter
    #-}
preOrderLogic
    ∷ ( FiniteBits a
      , Vector v a
      )
    ⇒ Bool
    → (v a, v a, v a)
    -- ^ Parent Final       Alignment
    → (v a, v a, v a)
    -- ^ Parent Preliminary Context
    → (v a, v a, v a)
    -- ^ Child  Preliminary Context
    → (v a, v a, v a)
    -- ^ Child  Final       Alignment
preOrderLogic isLeftChild pAlignment pContext cContext = forceDynamicCharacter $ unsafeCharacterBuiltByST caLen f
    where
        ccLen = fromEnum $ characterLength cContext
        caLen = characterLength pAlignment
        indices = [0 .. fromEnum caLen - 1]

        f char = pokeTempChar char *> fillTempChar char

        -- We set an arbitrary element from the parent alignemnt context to the first element of the temporary character.
        -- This ensures that at least one element can be querried for the determining "nil" and "gap" values.
        pokeTempChar char = setFrom pAlignment char 0 0

        -- Select character building function
        fillTempChar
            | isMissing cContext = missingλ -- Missing case is all gaps
            | otherwise = alignmentλ -- Standard pre-order logic
        missingλ char =
            forM_ indices (char `setGapped`)

        alignmentλ char = do
            j' ← newSTRef 0
            k' ← newSTRef 0

            forM_ indices $ \i → do
                k ← readSTRef k'
                if k > ccLen || pAlignment `isGapped` i
                    then char `setGapped` i
                    else do
                        j ← readSTRef j'
                        modifySTRef j' succ
                        -- Remember that 'Delete' leaves 'voids' in the 'left' character.
                        if pAlignment `isAlign` i
                            || (not isLeftChild && pAlignment `isDelete` i && pContext `isDelete` j)
                            || (isLeftChild && pAlignment `isInsert` i && pContext `isInsert` j)
                            then modifySTRef k' succ *> setFrom cContext char k i
                            else char `setGapped` i
