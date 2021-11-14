{-# LANGUAGE FlexibleContexts #-}

module DirectOptimization.PreOrder
  ( preOrderLogic
  ) where

import Bio.DynamicCharacter
import Control.Monad
import Data.Bits
import Data.STRef
import Data.Vector.Generic (Vector)


-- |
-- Faithful translation of Algorithm 8 (Non-root Node Alignment) from the
-- "Efficient Implied Alignment" paper found at:
-- <https://doi.org/10.1186/s12859-020-03595-2 10.1186/s12859-020-03595-2 DOI>
{-# INLINEABLE preOrderLogic #-}
{-# SPECIALISE preOrderLogic :: Bool -> SlimDynamicCharacter -> SlimDynamicCharacter -> SlimDynamicCharacter -> SlimDynamicCharacter #-}
{-# SPECIALISE preOrderLogic :: Bool -> WideDynamicCharacter -> WideDynamicCharacter -> WideDynamicCharacter -> WideDynamicCharacter #-}
{-# SPECIALISE preOrderLogic :: Bool -> HugeDynamicCharacter -> HugeDynamicCharacter -> HugeDynamicCharacter -> HugeDynamicCharacter #-}
preOrderLogic
  :: ( FiniteBits a
     , Vector v a
     )
  => Bool
  -> (v a, v a, v a) -- ^ Parent Final       Alignment
  -> (v a, v a, v a) -- ^ Parent Preliminary Context
  -> (v a, v a, v a) -- ^ Child  Preliminary Context
  -> (v a, v a, v a) -- ^ Child  Final       Alignment
preOrderLogic isLeftChild pAlignment pContext cContext = unsafeCharacterBuiltByST caLen f
  where
    paLen = fromEnum $ caLen
    ccLen = fromEnum $ characterLength cContext
    caLen = characterLength pAlignment

    -- Select character building function
    f | isMissing cContext = missingλ   -- Missing case is all gaps
      | otherwise          = alignmentλ -- Standard pre-order logic

    missingλ char = 
      forM_ [0 .. paLen - 1] (char `setGapped`)

    alignmentλ char =  do
        j'  <- newSTRef 0
        k'  <- newSTRef 0

        forM_ [0 .. paLen - 1] $ \i -> do
            k <- readSTRef k'
            if   k > ccLen || pAlignment `isGapped` i
            then char `setGapped` i
            else do
                j <- readSTRef j'
                modifySTRef j' succ
                if    pAlignment `isAlign` i
                  || (    isLeftChild && pAlignment `isDelete` i && pContext `isDelete` j)
                  || (not isLeftChild && pAlignment `isInsert` i && pContext `isInsert` j)
                then char `setGapped` i
                else modifySTRef k' succ *> setFrom cContext char k i
