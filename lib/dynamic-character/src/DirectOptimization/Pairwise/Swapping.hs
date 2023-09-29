{-# LANGUAGE ApplicativeDo #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE Strict #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UnboxedTuples #-}

{- |
Module      :  DirectOptimization.Pairwise.Swapping
Copyright   :  (c) 2015-2021 Ward Wheeler
License     :  BSD-style

Maintainer  :  wheeler@amnh.org
Stability   :  provisional
Portability :  portable

Direct optimization pairwise alignment using the Needleman-Wunsch algorithm.
These functions will allocate an M * N matrix.
-}
module DirectOptimization.Pairwise.Swapping (
    Direction (),
    swappingDO,
    buildDirectionMatrix,
    minimumCostDirection,
) where

import Bio.DynamicCharacter
import Bio.DynamicCharacter.Element (WideState)
import Control.Monad.ST
import Data.Bits
import Data.Foldable
import Data.Matrix.Unboxed (Matrix, unsafeFreeze)
import Data.Matrix.Unboxed.Mutable qualified as M
import Data.Vector qualified as V
import Data.Vector.Generic (Vector, (!))
import Data.Vector.Generic qualified as GV
import Data.Vector.Unboxed qualified as UV
import Data.Vector.Unboxed.Mutable qualified as MUV
import DirectOptimization.Pairwise.Direction
import DirectOptimization.Pairwise.Internal


{- |
Performs a naive direct optimization.
Takes in two characters to run DO on and an overlap function
Returns an assignment character, the cost of that assignment, the assignment
character with gaps included, the aligned version of the first input character,
and the aligned version of the second input character. The process for this
algorithm is to generate a traversal matrix, then perform a traceback.
-}
{-# SCC swappingDO #-}
{-# INLINEABLE swappingDO #-}
{-# SPECIALIZE swappingDO ∷
    (WideState → WideState → (WideState, Word)) → WideDynamicCharacter → WideDynamicCharacter → (Word, WideDynamicCharacter)
    #-}
{-# SPECIALIZE swappingDO ∷
    (HugeState → HugeState → (HugeState, Word)) → HugeDynamicCharacter → HugeDynamicCharacter → (Word, HugeDynamicCharacter)
    #-}
swappingDO
    ∷ ( FiniteBits e
      , Ord (v e)
      , Vector v e
      )
    ⇒ TCMλ e
    → OpenDynamicCharacter v e
    → OpenDynamicCharacter v e
    → (Word, OpenDynamicCharacter v e)
swappingDO = directOptimizationFromDirectionMatrix buildDirectionMatrix


{-# SCC buildDirectionMatrix #-}
{-# INLINEABLE buildDirectionMatrix #-}
{-# SPECIALIZE buildDirectionMatrix ∷
    WideState
    → (WideState → WideState → (WideState, Word))
    → UV.Vector WideState
    → UV.Vector WideState
    → (Word, Matrix Direction)
    #-}
{-# SPECIALIZE buildDirectionMatrix ∷
    HugeState
    → (HugeState → HugeState → (HugeState, Word))
    → V.Vector HugeState
    → V.Vector HugeState
    → (Word, Matrix Direction)
    #-}
buildDirectionMatrix
    ∷ (Vector v e)
    ⇒ e
    -- ^ Gap state
    → TCMλ e
    -- ^ Metric between states producing the medoid of states.
    → v e
    -- ^ Shorter dynamic character related to the "left column"
    → v e
    -- ^ Longer  dynamic character related to the "top row"
    → (Word, Matrix Direction)
buildDirectionMatrix gap tcmλ lesserLeft longerTop = fullMatrix
    where
        costλ x y = snd $ tcmλ x y
        rows = GV.length lesserLeft + 1
        cols = GV.length longerTop + 1

        fullMatrix = runST $ do
            mDir ← M.new (rows, cols)
            vOne ← MUV.new cols
            vTwo ← MUV.new cols

            let write v p@(~(_, j)) c d = MUV.unsafeWrite v j c *> M.unsafeWrite mDir p d

            write vOne (0, 0) 0 DiagArrow

            -- Special case the first row
            -- We need to ensure that there are only Left Arrow values in the directional matrix.
            -- We can also reduce the number of comparisons the first row makes from 3 to 1,
            -- since the diagonal and upward values are "out of bounds."
            for_ [1 .. cols - 1] $ \j →
                let topElement = longerTop ! (j - 1)
                    firstCellCost = costλ gap topElement
                in  do
                        firstPrevCost ← MUV.unsafeRead vOne (j - 1)
                        write vOne (0, j) (firstCellCost + firstPrevCost) LeftArrow

            for_ [1 .. rows - 1] $ \i →
                let (prev, curr)
                        | odd i = (vOne, vTwo)
                        | otherwise = (vTwo, vOne)
                    leftElement = lesserLeft ! (i - 1)
                    -- Special case the first cell of each row
                    -- We need to ensure that there are only Up Arrow values in the directional matrix.
                    -- We can also reduce the number of comparisons the first row makes from 3 to 1,
                    -- since the diagonal and leftward values are "out of bounds."
                    firstCellCost = costλ leftElement gap
                in  do
                        firstPrevCost ← MUV.unsafeRead prev 0
                        write curr (i, 0) (firstCellCost + firstPrevCost) UpArrow
                        -- Finish special case for first cell of each row
                        -- Begin processing all other cells in the curr vector
                        for_ [1 .. cols - 1] $ \j →
                            let topElement = longerTop ! (j - 1)
                                deleteCost = costλ gap topElement
                                alignCost = costλ leftElement topElement
                                insertCost = costλ leftElement gap
                            in  do
                                    diagCost ← MUV.unsafeRead prev $ j - 1
                                    topCost ← MUV.unsafeRead prev j
                                    leftCost ← MUV.unsafeRead curr $ j - 1
                                    let (# c, d #) =
                                            minimumCostDirection
                                                (deleteCost + leftCost)
                                                (alignCost + diagCost)
                                                (insertCost + topCost)
                                    write curr (i, j) c d

            let v
                    | odd rows = vOne
                    | otherwise = vTwo

            c ← MUV.unsafeRead v (cols - 1)
            m ← unsafeFreeze mDir
            pure (c, m)
