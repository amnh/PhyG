{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE Strict #-}

module DirectOptimization.PreOrder
  ( preOrderLogic
  ) where

import Bio.DynamicCharacter
import Control.Monad
import Data.Bits
import Data.STRef
import Data.Vector.Generic (Vector)

import Debug.Trace (trace)

{-
tr :: Show a => String -> a -> a
tr key val = trace (key ": \t" <> show val) val


trChar
  :: ( FiniteBits a
     , Show a
     , Vector v a
     )
  => String
  -> OpenDynamicCharacter v a
  -> OpenDynamicCharacter v a
trChar label char = trace (key ":\n\t" <> renderDynamicCharacter val <> "\n") val
-}

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
     , Show a
     , Vector v a
     )
  => Bool
  -> OpenDynamicCharacter v a -- ^ Parent Final       Alignment
  -> OpenDynamicCharacter v a -- ^ Parent Preliminary Context
  -> OpenDynamicCharacter v a -- ^ Child  Preliminary Context
  -> OpenDynamicCharacter v a -- ^ Child  Final       Alignment
preOrderLogic isLeftChild pAlignment pContext cContext = renderComputation . forceDynamicCharacter $ unsafeCharacterBuiltByST caLen f
  where
    renderComputation result =
        let tagged key = "\t" <> key <> ":\t"
            taggedShow key val = tagged key <> show val
            taggedChar key val = tagged key <> renderDynamicCharacter val <> "\n"
            rendering = unlines
                [ "INPUT:"
                , "  preOrderLogic with ..."
                , taggedShow "isLeftChild" isLeftChild
                , taggedChar "pAlignment"  pAlignment
                , taggedChar "pContext"    pContext
                , taggedChar "cContext"    cContext
                , "OUTPUT:"
                , taggedChar "return" result
                ]
        in  result -- trace rendering result

    ccLen   = fromEnum $ characterLength cContext
    caLen   = characterLength pAlignment
    indices = [ 0 .. fromEnum caLen - 1 ]

    f char = pokeTempChar char *> fillTempChar char

    -- We set an arbitrary element from the parent alignemnt context to the first element of the temporary character.
    -- This ensures that at least one element can be querried for the determining "nil" and "gap" values.
    pokeTempChar char = setFrom pAlignment char 0 0

    -- Select character building function
    fillTempChar
      | isMissing cContext = missingλ   -- Missing case is all gaps
      | otherwise          = alignmentλ -- Standard pre-order logic

    missingλ char = 
      forM_ indices (char `setGapped`)

    alignmentλ char =  do
        j'  <- newSTRef 0
        k'  <- newSTRef 0

        forM_ indices $ \i -> do
            k <- readSTRef k'
            if   k > ccLen || pAlignment `isGapped` i
            then char `setGapped` i
            else do
                j <- readSTRef j'
                modifySTRef j' succ
                -- Remember that 'Delete' leaves 'voids' in the 'left' character.
                if    pAlignment `isAlign` i
                  || (not isLeftChild && pAlignment `isDelete` i && pContext `isDelete` j)
                  || (    isLeftChild && pAlignment `isInsert` i && pContext `isInsert` j)
                then modifySTRef k' succ *> setFrom cContext char k i
                else char `setGapped` i
