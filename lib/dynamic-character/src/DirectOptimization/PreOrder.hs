{-# LANGUAGE ApplicativeDo #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# OPTIONS_GHC -foptimal-applicative-do #-}

{- HLINT ignore "Redundant pure" -}

module DirectOptimization.PreOrder (
    preOrderLogic,
) where

import Bio.DynamicCharacter
import Control.Monad
import Control.Monad.ST
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
    ∷ ∀ a v
     . ( FiniteBits a
       , Vector v a
       )
    ⇒ Bool
    → OpenDynamicCharacter v a
    -- ^ Parent Final       Alignment
    → OpenDynamicCharacter v a
    -- ^ Parent Preliminary Context
    → OpenDynamicCharacter v a
    -- ^ Child  Preliminary Context
    → OpenDynamicCharacter v a
    -- ^ Child  Final       Alignment
preOrderLogic isLeftChild pAlignment pContext cContext =
    let pAlignment' = forceDynamicCharacter pAlignment
        pContext'   = forceDynamicCharacter pContext
        cContext'   = forceDynamicCharacter cContext

        caLen = characterLength pAlignment'
        ccLen = fromEnum $ characterLength cContext'
        range = [0 .. fromEnum caLen - 1]

        -- Constructor definition supplied to 'unsafeCharacterBuiltByST'.
        builder ∷ ∀ s. TempOpenDynamicCharacter (ST s) v a → ST s ()
        builder char = pokeTempChar char *> fillTempChar char

        -- We set an arbitrary element from the parent alignment context
        -- to the first element of the temporary character. This ensures
        -- that at least one element can be querried for the determining
        -- "nil" and "gap" values.
        pokeTempChar ∷ ∀ s. TempOpenDynamicCharacter (ST s) v a → ST s ()
        pokeTempChar char = setFrom pAlignment' char 0 0

        -- Select character building function
        fillTempChar ∷ ∀ s. TempOpenDynamicCharacter (ST s) v a → ST s ()
        fillTempChar
            | isMissing cContext' = missingλ -- Missing case is all gaps
            | otherwise = alignmentλ -- Standard pre-order logic

        -- Construct a missing character value of appropriate length.
        missingλ ∷ ∀ s. TempOpenDynamicCharacter (ST s) v a → ST s ()
        missingλ char =
            forM_ range (char `setGapped`)

        -- Construct alignment derived from parent pre-order and self post-order.
        alignmentλ ∷ ∀ s. TempOpenDynamicCharacter (ST s) v a → ST s ()
        alignmentλ char = do
            j' ← newSTRef 0
            k' ← newSTRef 0
            forM_ range $ \i → do
                k ← readSTRef k'
                if k > ccLen || pAlignment' `isGapped` i
                    then char `setGapped` i
                    else do
                        j ← readSTRef j'
                        modifySTRef j' succ
                        -- Remember that 'Delete' leaves 'voids' in the 'left' character.
                        if pAlignment' `isAlign` i
                            || (not isLeftChild && pAlignment' `isDelete` i && pContext' `isDelete` j)
                            || (isLeftChild && pAlignment' `isInsert` i && pContext' `isInsert` j)
                            then modifySTRef k' succ *> setFrom cContext' char k i
                            else char `setGapped` i
                pure () -- For ApplicativeDo
            pure () -- For ApplicativeDo
    in  forceDynamicCharacter $ unsafeCharacterBuiltByST caLen builder
