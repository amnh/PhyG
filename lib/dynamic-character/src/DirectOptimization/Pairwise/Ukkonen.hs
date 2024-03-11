{-# LANGUAGE ApplicativeDo #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UnboxedTuples #-}

{- |
Module      :  DirectOptimization.Pairwise.Ukkonen
Copyright   :  (c) 2015-2021 Ward Wheeler
License     :  BSD-style

Maintainer  :  wheeler@amnh.org
Stability   :  provisional
Portability :  portable

Direct optimization pairwise alignment using the Needleman-Wunsch algorithm.
These functions will allocate an M * N matrix.
-}
module DirectOptimization.Pairwise.Ukkonen (
    Direction (),
    ukkonenDO,
    createUkkonenMethodMatrix,
) where

import Bio.DynamicCharacter
import Bio.DynamicCharacter.Element (WideState, HugeState)
import Bio.DynamicCharacter.Measure
import Control.Monad (unless, when)
import Control.Monad.Loops (iterateUntilM, whileM_)
import Control.Monad.Primitive
import Control.Monad.ST
import Data.Bits
import Data.Foldable
import Data.Matrix.Unboxed (Matrix, unsafeFreeze)
import Data.Matrix.Unboxed.Mutable (MMatrix)
import Data.Matrix.Unboxed.Mutable qualified as M
import Data.STRef
import Data.Vector qualified as V
import Data.Vector.Generic (Vector, (!))
import Data.Vector.Generic qualified as GV
import Data.Vector.Unboxed qualified as UV
import DirectOptimization.Pairwise.Internal
import DirectOptimization.Pairwise.Swapping


{- |
Performs a naive direct optimization.
Takes in two characters to run DO on and an overlap function
Returns an assignment character, the cost of that assignment, the assignment
character with gaps included, the aligned version of the first input character,
and the aligned version of the second input character. The process for this
algorithm is to generate a traversal matrix, then perform a traceback.
-}
{-# SCC ukkonenDO #-}
{-# INLINEABLE ukkonenDO #-}
{-# SPECIALIZE ukkonenDO ∷
    Word
    → (WideState → WideState → (WideState, Word))
    → WideDynamicCharacter
    → WideDynamicCharacter
    → (Word, WideDynamicCharacter)
    #-}
{-# SPECIALIZE ukkonenDO ∷
    Word
    → (HugeState → HugeState → (HugeState, Word))
    → HugeDynamicCharacter
    → HugeDynamicCharacter
    → (Word, HugeDynamicCharacter)
    #-}
ukkonenDO
    ∷ ( FiniteBits e
      , Ord (v e)
      , Vector v e
      )
    ⇒ Word
    -- ^ Coefficient value, representing the /minimum/ transition cost from a state to gap
    → TCMλ e
    -- ^ Metric between states producing the medoid of states.
    → OpenDynamicCharacter v e
    -- ^ /1st/ dynamic character
    → OpenDynamicCharacter v e
    -- ^ /2nd/ dynamic character
    → (Word, OpenDynamicCharacter v e)
ukkonenDO coefficient tcmλ char1 char2
    | noGainFromUkkonenMethod = buildFullMatrix
    | otherwise = buildBandMatrix
    where
        buildFullMatrix = swappingDO tcmλ char1 char2
        buildBandMatrix = directOptimizationFromDirectionMatrix ukkonenBandλ tcmλ char1 char2

        ukkonenBandλ = createUkkonenMethodMatrix coefficient inputGapAmbiguities

        ~(_, _, _, x, y) = measureCharactersWithoutGaps char1 char2

        lesser = extractMediansGapped x
        longer = extractMediansGapped y

        -- /O(1)/
        --
        -- If the longer character is 50% larger than the shorter character, then
        -- there is no point in using the barriers. Rather, we fill the full matrix
        -- immediately.
        --
        -- Additionally, if the shorter sequence is of length 4 or less, then the
        -- initial barrier will be set adjacent to or beyond the lower left and
        -- upper right corners.
        --
        -- Also, a threshold coefficient is computed as the minimal indel cost from
        -- any symbol in the alphabet to gap. However, if the indel cost for any
        -- symbol is zero, the algorithm will hang, and a naive approach must be taken.
        --
        -- Lastly, if the sum of the gaps in both strings is equal to or exceeds the
        -- length of the longer string, then the threshold criteria will never be met
        -- by definition.
        --
        -- Do not perform Ukkonen's algorithm if and only if:
        --
        -- > longerLen >= 1.5 * lesserLen
        --     OR
        -- > lesserLen <= 4
        --     OR
        -- > coefficient == 0
        --
        noGainFromUkkonenMethod =
            lesserLen <= 4
                || 2 * longerLen >= 3 * lesserLen
                || coefficient == 0
                || inputGapAmbiguities >= longerLen
            where
                longerLen = vLength longer
                lesserLen = vLength lesser

        -- /O(n + m)/
        --
        -- NOTE: There will be no *unambiguous* gap elements in the dynamic characters!
        --   However, there may be abiguous elements which contain gap as a possibility.
        --
        -- If one or more of the character elements contained a gap, diagonal
        -- directions in the matrix have an "indel" cost. 'gapsPresentInInputs' is
        -- necessary in order to decrement the threshold value to account for this.
        -- This was not described in Ukkonen's original paper, as the inputs were assumed
        -- not to contain any gaps.
        inputGapAmbiguities = char1Gaps + char2Gaps
            where
                char1Gaps = countGaps lesser
                char2Gaps = countGaps longer
                countGaps = vLength . GV.filter maybeGap
                maybeGap = (`testBit` 0) -- Zero is the gap bit!
        vLength = toEnum . GV.length


{- |
/O( (n - m + 1 ) * log(n - m + 1) )/, /n/ >= /m/

Generates an /optimal/, partially-filled-in matrix using Ukkonen's string
edit distance algorithm.

Note that the threshold value is lowered more than described in Ukkonen's
paper. This is to handle input elements that contain a gap. In Ukkonen's
original description of the algorithm, there was a subtle assumption that
input did not contain any gap symbols.
-}
{-# SCC createUkkonenMethodMatrix #-}
{-# INLINEABLE createUkkonenMethodMatrix #-}
{-# SPECIALIZE createUkkonenMethodMatrix ∷
    Word
    → Word
    → WideState
    → (WideState → WideState → (WideState, Word))
    → UV.Vector WideState
    → UV.Vector WideState
    → (Word, Matrix Direction)
    #-}
{-# SPECIALIZE createUkkonenMethodMatrix ∷
    Word
    → Word
    → HugeState
    → (HugeState → HugeState → (HugeState, Word))
    → V.Vector HugeState
    → V.Vector HugeState
    → (Word, Matrix Direction)
    #-}
createUkkonenMethodMatrix
    ∷ (Vector v e)
    ⇒ Word
    -- ^ Coefficient value, representing the /minimum/ transition cost from a state to gap
    → Word
    -- ^ Number of abiguous elements in both inputs which contained gap as a possible state
    → e
    -- ^ Gap State
    → TCMλ e
    -- ^ Metric between states producing the medoid of states.
    → v e
    -- ^ Shorter dynamic character
    → v e
    -- ^ Longer  dynamic character
    → (Word, Matrix Direction)
createUkkonenMethodMatrix minimumIndelCost inputGapAmbiguities gap tcmλ lesserLeft longerTop = finalMatrix
    where
        -- General values that need to be in scope for the recursive computations.
        lesserLen = GV.length lesserLeft
        longerLen = GV.length longerTop

        -- We start the offset at four rather than at one so that the first doubling
        -- isn't trivially small.
        startOffset = 2

        -- /O(1)/
        --
        -- Necessary to compute the width of a row in the barrier-constrained matrix.
        quasiDiagonalWidth = toEnum $ differenceInLength + 1
            where
                differenceInLength = longerLen - lesserLen

        extra = (inputGapAmbiguities +)

        finalMatrix = runST $ do
            (mCost, mDir) ← buildInitialBandedMatrix gap tcmλ lesserLeft longerTop $ extra startOffset
            let getAlignmentCost = M.unsafeRead mCost (lesserLen, longerLen)
            offsetRef ← newSTRef startOffset

            let needToResizeBand = do
                    offset ← readSTRef offsetRef
                    -- If the filled row width exceeds the actual row length,
                    -- Then clearly we are done as we have filled the entire matrix.
                    if quasiDiagonalWidth + extra offset > toEnum longerLen
                        then pure False
                        else
                            let partialWidth = quasiDiagonalWidth + offset
                                -- Value that the alignment cost must be less than
                                threshold -- The threshold value must be non-negative
                                    | partialWidth <= inputGapAmbiguities = 0
                                    | otherwise = minimumIndelCost * (partialWidth - inputGapAmbiguities)
                            in  (threshold <=) <$> getAlignmentCost

            whileM_ needToResizeBand $ do
                previousOffset ← readSTRef offsetRef
                let currentOffset = previousOffset `shiftL` 1 -- Multiply by 2
                writeSTRef offsetRef currentOffset
                expandBandedMatrix
                    gap
                    tcmλ
                    lesserLeft
                    longerTop
                    mCost
                    mDir
                    (extra previousOffset)
                    (extra currentOffset)

            c ← getAlignmentCost
            m ← unsafeFreeze mDir
            pure (c, m)


{-# SCC buildInitialBandedMatrix #-}
buildInitialBandedMatrix
    ∷ (Vector v e)
    ⇒ e
    -- ^ Gap
    → TCMλ e
    -- ^ Metric between states producing the medoid of states.
    → v e
    -- ^ Shorter dynamic character
    → v e
    -- ^ Longer dynamic character
    → Word
    → ST s (MMatrix s Word, MMatrix s Direction)
buildInitialBandedMatrix gap tcmλ lesserLeft longerTop o = fullMatrix
    where
        (offset, costλ, rows, cols, width, quasiDiagonalWidth) = ukkonenConstants tcmλ lesserLeft longerTop o

        fullMatrix = do
            ---------------------------------------
            -- Allocate required space           --
            ---------------------------------------

            mCost ← M.new (rows, cols)
            mDir ← M.new (rows, cols)

            ---------------------------------------
            -- Define some generalized functions --
            ---------------------------------------
            let ~(readCost, write, internalCell, leftColumn, leftBoundary, rightBoundary, rightColumn) =
                    edgeCellDefinitions gap costλ longerTop mCost mDir

            -- Define how to compute values to an entire row of the Ukkonen matrix.
            let writeRow i =
                    -- Precompute some values that will be used for the whole row
                    let start = max 0 $ i - offset
                        stop = min (cols - 1) $ i - offset + width - 1
                        leftElement = lesserLeft ! (i - 1)
                        insertCost = costλ leftElement gap

                        -- Each row in the matrix with values in the band has 'width' cells.
                        -- However, the band runs off the left end of the matrix for the first
                        -- several rows of the matrix. How many rows, though?
                        -- There are exactly 'offset' number of cells left of first matrix column
                        -- in the first row. The number of cells to the left of the matrix
                        -- decreases by one in each subsequent row. The means that the first
                        -- row, and then the next 'offset' number of rows require special handling
                        -- of the boundary. The last row index requiring special handling is index
                        -- 'offset'. Subsequent rows have the band begin at least one cell away
                        -- from the matrix boundary.
                        firstCell
                            | i <= offset = leftColumn
                            | otherwise = leftBoundary

                        lastCell
                            | i <= cols - quasiDiagonalWidth - offset = rightBoundary
                            | otherwise = rightColumn
                    in  do
                            -- Write to the first cell of the Ukkonen band
                            firstCell leftElement insertCost i start
                            -- Write to the all the intermediary cells of the Ukkonen band
                            for_ [start + 1 .. stop - 1] $ \j →
                                internalCell leftElement insertCost i j
                            -- Write to the last cell of the Ukkonen band
                            lastCell leftElement insertCost i stop

            ---------------------------------------
            -- Compute all values of the matrix  --
            ---------------------------------------

            -- Write to the origin to seed the first row.
            write (0, 0) (# 0, DiagArrow #)

            -- Each row in the matrix with values in the band has 'width' cells.
            -- However, the band runs off the end of the matrix for the first & last
            -- rows of the matrix. We subtract the 'offest' from the 'width' because
            -- there are exactly 'offset' number of cells left of first matrix column.
            -- Hence the top row's width is 'width' minus 'offset'. The last cell index
            -- in the top row of the band is 'width' minus 'offset' minus 1.
            let topRowWidth = width - offset
            let topRowWrite !j !cost = write (0, j) (# cost, LeftArrow #)

            -- Write the first row to seed subsequent rows.
            for_ [1 .. min (cols - 1) (topRowWidth - 1)] $ \j →
                let topElement = longerTop ! (j - 1)
                    firstCellCost = costλ gap topElement
                in  do
                        firstPrevCost ← readCost 0 $ j - 1
                        topRowWrite j $ firstCellCost + firstPrevCost

            -- Loop through the remaining rows.
            for_ [1 .. rows - 1] writeRow

            -- Return the matricies for possible expansion
            pure (mCost, mDir)


{-# SCC expandBandedMatrix #-}


{- |
Given a partially computed alignment matrix,
will expand the computed region to the new specified offset.


Dimensions: 13 ⨉ 17
 ⊗ ┃  ⁎ α1 α2 α3 α4 α5 α6 α7 α8 α9 α0 α1 α2 α3 α4 α5 α6
━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 ⁎ ┃ 0↖ 0← 0← 0← 0← 0← 0← 0← 0← 0← 0← 0← 0← 0← 0← 0← 0←
α1 ┃ 0↑ 0↖ 0↖ 0↖ 0← 0← 0← 0← 0← 0← 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0←
α2 ┃ 0↑ 0↖ 0↖ 0↖ 0← 0← 0← 0← 0← 0← 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0←
α3 ┃ 0↑ 0↑ 0↑ 0↖ 0↖ 0↖ 0← 0← 0← 0← 0← 0← 0← 0← 0← 0← 0↖
α4 ┃ 0↑ 0↑ 0↑ 0↖ 0↖ 0← 0↖ 0← 0↖ 0← 0← 0← 0← 0← 0← 0← 0←
α5 ┃ 0↑ 0↑ 0↑ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0← 0← 0← 0← 0← 0← 0← 0←
α6 ┃ 0↑ 0↖ 0↖ 0↖ 0↑ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0←
α7 ┃ 0↑ 0↑ 0↑ 0↑ 0↖ 0↖ 0↖ 0↖ 0↖ 0← 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖
α8 ┃ 0↑ 0↑ 0↑ 0↑ 0↑ 0↖ 0← 0↖ 0↑ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖
α9 ┃ 0↑ 0↑ 0↑ 0↑ 0↖ 0↑ 0↖ 0← 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖
α0 ┃ 0↑ 0↑ 0↑ 0↑ 0↑ 0↖ 0↑ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖
α1 ┃ 0↑ 0↑ 0↑ 0↑ 0↖ 0↑ 0↖ 0↖ 0↖ 0← 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖
α2 ┃ 0↑ 0↖ 0↖ 0↖ 0↑ 0↑ 0↑ 0↖ 0↑ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖

     ┌───────────────w───────────────┐
     │              ┏━━━━━━━co━━━━━━━┪
     ┢━━━━━qd━━━━━━┓┠─po─┐┌────Δo────┨
 ⊗ ┃ ┃0  1  2  3  4┃┃5  6││7  8  9 10┃11 12 13 14 15 16
━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 0 ┃ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒ ▒▒ ▒▒ ▒▒
 1 ┃ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒ ▒▒ ▒▒ ▒▒
 2 ┃ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒ ▒▒ ▒▒ ▒▒
 3 ┃ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒ ▒▒ ▒▒ ▒▒
 4 ┃ ▒▒ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒ ▒▒ ▒▒ ▒▒
 5 ┃ ▒▒ ▒▒ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒ ▒▒ ▒▒ ▒▒
 6 ┃ ▒▒ ▒▒ ▒▒ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒ ▒▒ ▒▒ ▒▒
 7 ┃    ▒▒ ▒▒ ▒▒ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒ ▒▒ ▒▒
 8 ┃       ▒▒ ▒▒ ▒▒ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒ ▒▒
 9 ┃          ▒▒ ▒▒ ▒▒ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒
 0 ┃             ▒▒ ▒▒ ▒▒ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓
 1 ┃                ▒▒ ▒▒ ▒▒ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓
 2 ┃                   ▒▒ ▒▒ ▒▒ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██


w  : Width
qd : Quasi-diagonal
co : Current Offset
po : Previous Offset
Δo : Difference in Offset

Note:
w  = qd + co
co = po + Δo

And often:
co = 2*po = 2*Δo

██ : The core band
      * Previously computed, sections may need to be recomputed

▓▓ : The previous extension
      * Previously computed, sections may need to be recomputed

▒▒ : The new extension
      * Needs to be computed
-}
expandBandedMatrix
    ∷ (Vector v e)
    ⇒ e
    -- ^ Gap state
    → TCMλ e
    -- ^ Metric between states producing the medoid of states.
    → v e
    -- ^ Shorter dynamic character
    → v e
    -- ^ Longer  dynamic character
    → MMatrix s Word
    → MMatrix s Direction
    → Word
    → Word
    → ST s ()
expandBandedMatrix gap tcmλ lesserLeft longerTop mCost mDir po co = updatedBand
    where
        (offset, costλ, rows, cols, width, qd) = ukkonenConstants tcmλ lesserLeft longerTop co
        prevOffset = fromEnum po

        updatedBand = do
            ---------------------------------------
            -- Allocate mutable state variables  --
            ---------------------------------------

            tailStart ← newSTRef cols

            t0' ← newSTRef (-1)
            t1' ← newSTRef $ qd + fromEnum po

            ---------------------------------------
            -- Define some generalized functions --
            ---------------------------------------
            let ~(readCost, write, internalCell, leftColumn, leftBoundary, rightBoundary, rightColumn) =
                    edgeCellDefinitions gap costλ longerTop mCost mDir

            let computeCell !leftElement !insertCost !i !j =
                    {-# SCC recomputeCell #-}
                    let !topElement = longerTop ! (j - 1)
                        !deleteCost = costλ gap topElement
                        !alignCost = costλ leftElement topElement
                    in  do
                            diagCost ← readCost (i - 1) $ j - 1
                            topCost ← readCost (i - 1) j
                            leftCost ← readCost i $ j - 1
                            oldCost ← readCost i j
                            let !e@(# c, _ #) =
                                    minimumCostDirection
                                        (deleteCost + leftCost)
                                        (alignCost + diagCost)
                                        (insertCost + topCost)
                            write (i, j) e
                            pure (c == oldCost, j + 1)
            --                  pure (c /= oldCost, j+1)

            let recomputeRange leftElement insertCost i x y = do
                    lastDiff ← newSTRef 0
                    for_ [x .. y] $ \j → do
                        (same, _) ← computeCell leftElement insertCost i j
                        unless same $ writeSTRef lastDiff j
                    readSTRef lastDiff

            -- Define how to compute values to an entire row of the Ukkonen matrix.
            let extendRow i =
                    -- Precopmute some values that will be used for the whole row
                    let start0 = max 0 $ i - offset
                        start3 = min cols $ i + width - offset - prevOffset - 1
                        goUpTo = max 0 (i - prevOffset) - 1
                        stop = min (cols - 1) $ i + width - offset - 1
                        leftElement = lesserLeft ! (i - 1)
                        insertCost = costλ leftElement gap
                        firstCell
                            | i <= offset = leftColumn
                            | otherwise = leftBoundary

                        lastCell
                            | i <= cols - qd - offset = rightBoundary
                            | otherwise = rightColumn

                        b0 = start0
                        e0 = goUpTo
                        b1 = start3
                        e1 = stop

                        continueRecomputing (same, j) = same || j >= stop
                        computeCell' ~(_, j) = computeCell leftElement insertCost i j
                        internalCell' j = internalCell leftElement insertCost i j
                        recomputeUntilSame j = snd <$> iterateUntilM continueRecomputing computeCell' (False, j)
                    in  do
                            -- First, we fill in 0 or more cells of the left region of
                            -- the expanded band. This is the region [b0, e0] computed
                            -- above.
                            --  ⊗ ┃  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
                            -- ━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                            --  0 ┃ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  1 ┃ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  2 ┃ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  3 ┃ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  4 ┃ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --      b0    e0
                            --     ┏━━━━━━━━┓
                            --  5 ┃ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --
                            --  6 ┃ ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  7 ┃    ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒
                            --  8 ┃       ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒
                            --  9 ┃          ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒
                            -- 10 ┃             ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██
                            -- 11 ┃                ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██
                            -- 12 ┃                   ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██
                            --

                            -- Conditionally write to the first cell of the Ukkonen band
                            if i > prevOffset
                                then firstCell leftElement insertCost i b0
                                else pure ()

                            for_ [b0 + 1 .. e0] internalCell'

                            -- Next, we assign to s0 the value t0 from the previous row.
                            -- The cell t0 is up to where the values were recomputed in
                            -- the previous row.
                            -- We recompute the cells in the range [e0 + 1, s0].
                            -- We assign to t0 the last cell in the range [s1, s2] which
                            -- was updated for the next row.
                            -- We remember cells t0 for the next row.
                            --  ⊗ ┃  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
                            -- ━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                            --  0 ┃ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  1 ┃ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  2 ┃ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  3 ┃ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  4 ┃ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --            e0    s0
                            --              ┏━━━━━┓
                            --  5 ┃ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --
                            --  6 ┃ ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  7 ┃    ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒
                            --  8 ┃       ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒
                            --  9 ┃          ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒
                            -- 10 ┃             ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██
                            -- 11 ┃                ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██
                            -- 12 ┃                   ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██
                            --
                            --
                            s0 ← (\x → min (x + 1) e1) <$> readSTRef t0'
                            writeSTRef t0' (-1)

                            when (s0 > e0 && toEnum i > po) $
                                recomputeRange leftElement insertCost i (e0 + 1) s0 >>= writeSTRef t0'
                            t0 ← readSTRef t0'

                            -- If s0 = t0, we recompute the cell (s0 + 1).
                            -- If the cost is the same, we stop here and remember the cell
                            -- before we stopped.
                            -- If the cost is not the same, we update cell (s0 + 1) and
                            -- move on to (s0 + 2).
                            -- This procedure continues until (s0 + n) has the same cost
                            -- as before, or *until we reach b1.*
                            -- We remember the cell (s0 + n - 1) as t0 for the next row.
                            --  ⊗ ┃  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
                            -- ━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                            --  0 ┃ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  1 ┃ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  2 ┃ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  3 ┃ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  4 ┃ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --                  s0    t0
                            --                    ╔═════╗
                            --  5 ┃ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --
                            --  6 ┃ ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  7 ┃    ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒
                            --  8 ┃       ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒
                            --  9 ┃          ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒
                            -- 10 ┃             ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██
                            -- 11 ┃                ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██
                            -- 12 ┃                   ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██
                            --
                            if s0 == t0 && s0 > 0
                                then recomputeUntilSame (s0 + 1) >>= writeSTRef t0' . pred
                                else
                                    if s0 <= e0 && e0 > 0
                                        then recomputeUntilSame (e0 + 1) >>= writeSTRef t0' . pred
                                        else pure ()

                            -- Next, we assign to s1 the value t1 from the previous row.
                            -- We also assign s2 the value t2 from the previous row.
                            -- The range [t1, t2] is where the values were recomputed in
                            -- the previous row.
                            -- We recompute the cells in the range [s1, s2].
                            -- We assign to t2 the last cell in the range [s1, s2] which
                            -- was updated for the next row.
                            -- We remember cells s1 as t1 for the next row.
                            --  ⊗ ┃  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
                            -- ━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                            --  0 ┃ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  1 ┃ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  2 ┃ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  3 ┃ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  4 ┃ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --                                 s1 s2
                            --                                ┏━━━━━┓
                            --  5 ┃ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --
                            --  6 ┃ ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  7 ┃    ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒
                            --  8 ┃       ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒
                            --  9 ┃          ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒
                            -- 10 ┃             ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██
                            -- 11 ┃                ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██
                            -- 12 ┃                   ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██
                            --
                            -- NOPE, Try again
                            --
                            -- Next, we assign to s1 the value t1 from the previous row.
                            -- If s1 is less than  t0, we assign to s1 the value t0 + 1.
                            -- This ensures that we do not start "behind" where we have
                            -- previously computed.
                            -- Then if s1 is greater than e1, we assign to s1 the
                            -- value e1. This ensures one cell is always written to.
                            -- We recompute the cells in the range [s1, b1 - 1].
                            -- If any cell in the range was updated, we assign to s1 to t1.
                            -- We remember cell t1 for the next row.
                            --  ⊗ ┃  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
                            -- ━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                            --  0 ┃ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  1 ┃ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  2 ┃ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  3 ┃ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  4 ┃ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --                                 s1       b1
                            --                                ┏━━━━━━━━┓
                            --  5 ┃ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --
                            --  6 ┃ ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  7 ┃    ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒
                            --  8 ┃       ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒
                            --  9 ┃          ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒
                            -- 10 ┃             ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██
                            -- 11 ┃                ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██
                            -- 12 ┃                   ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██
                            --
                            s1 ← do
                                a ← readSTRef t0'
                                b ← readSTRef t1'
                                pure . min e1 $ max a b

                            t1 ← recomputeRange leftElement insertCost i s1 $ b1 - 1

                            -- If no cells were updated, a zero value is returned.
                            -- In this case, the "last" updated cell for the next row is b1.
                            writeSTRef t1' $ if t1 == 0 then b1 else s1

                            -- Lastly, we fill in 0 or more cells of the left region of
                            -- the expanded band. This is the region [b1, e1] computed
                            -- above.
                            --  ⊗ ┃  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
                            -- ━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                            --  0 ┃ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  1 ┃ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  2 ┃ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  3 ┃ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  4 ┃ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --                                          b1       e1
                            --                                         ┏━━━━━━━━━━━┓
                            --  5 ┃ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --
                            --  6 ┃ ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒ ▒▒
                            --  7 ┃    ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒ ▒▒
                            --  8 ┃       ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒ ▒▒
                            --  9 ┃          ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██ ▒▒
                            -- 10 ┃             ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██ ██
                            -- 11 ┃                ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██ ██
                            -- 12 ┃                   ▒▒ ▒▒ ▒▒ ▒▒ ██ ██ ██ ██ ██ ██ ██
                            --
                            for_ [b1 .. e1 - 1] internalCell'

                            -- Conditionally write to the last cell of the Ukkonen band
                            if i < rows - fromEnum po
                                then lastCell leftElement insertCost i stop
                                else pure ()

            ---------------------------------------
            -- Compute all values of the matrix  --
            ---------------------------------------

            -- We start computation in the top row at the index equal to the
            -- quasi-diagonal width plus the previous offset. This is because the last
            -- cell of the top row which we computed was at the previous index since
            -- we previous computed an number of cells from the quasi-diagonal equal to
            -- the previous offset.
            let topRowStart = qd + prevOffset

            -- Each row in the matrix with values in the band has 'width' cells.
            -- However, the band runs off the end of the matrix for the first & last
            -- rows of the matrix. We subtract the 'offset' from the 'width' because
            -- there are exactly 'offset' number of cells left of first matrix column.
            -- Hence the top row's width is 'width' minus 'offset'. The last cell index
            -- in the top row of the band is 'width' minus 'offset' minus 1.
            let topRowWidth = width - offset

            -- Of course, we must be cetrain that we don't extend past the last column
            -- of the matrix. To prevent this, we take the minimum of the top row width
            -- and the number of columns. So the last index we will compute in the top
            -- row is the minimum of the two options minus one due to zero-indexing.
            let topRowCease = pred $ min cols topRowWidth

            -- Write out the left arrow value on the top row.
            let topRowWrite !j !cost = write (0, j) (# cost, LeftArrow #)

            -- Extend the first row to seed subsequent rows.
            for_ [topRowStart .. topRowCease] $ \j →
                let !topElement = longerTop ! (j - 1)
                    !firstCellCost = costλ gap topElement
                in  do
                        !firstPrevCost ← readCost 0 $ j - 1
                        topRowWrite j $ firstCellCost + firstPrevCost

            writeSTRef tailStart topRowStart

            -- Loop through the remaining rows.
            for_ [1 .. rows - 1] extendRow


edgeCellDefinitions
    ∷ ( PrimMonad m
      , Vector v e
      )
    ⇒ e
    -- ^ Gap state
    → (e → e → Word)
    -- ^ Distance between states
    → v e
    -- ^ Longer dynamic character
    → MMatrix (PrimState m) Word
    → MMatrix (PrimState m) Direction
    → ( Int → Int → m Word
      , (Int, Int) → (# Word, Direction #) → m ()
      , e → Word → Int → Int → m ()
      , e → Word → Int → Int → m ()
      , e → Word → Int → Int → m ()
      , e → Word → Int → Int → m ()
      , e → Word → Int → Int → m ()
      )
edgeCellDefinitions gap costλ longerTop mCost mDir =
    (readCost, write, internalCell, leftColumn, leftBoundary, rightBoundary, rightColumn)
    where
        -- Read the cost of a cell
        readCost = curry $ M.unsafeRead mCost

        -- Write to a single cell of the current vector and directional matrix simultaneously
        write !p (# !c, !d #) = M.unsafeWrite mCost p c *> M.unsafeWrite mDir p d

        -- Write to an internal cell (not on a boundary) of the matrix.
        internalCell !leftElement !insertCost !i !j =
            {-# SCC internalCell_expanding #-}
            let !topElement = longerTop ! (j - 1)
                !deleteCost = costλ gap topElement
                !alignCost = costλ leftElement topElement
            in  do
                    diagCost ← readCost (i - 1) $ j - 1
                    topCost ← readCost (i - 1) j
                    leftCost ← readCost i $ j - 1
                    let v =
                            minimumCostDirection
                                (deleteCost + leftCost)
                                (alignCost + diagCost)
                                (insertCost + topCost)
                    write (i, j) v

        -- Define how to compute the first cell of the first "offset" rows.
        -- We need to ensure that there are only Up Arrow values in the directional matrix.
        -- We can also reduce the number of comparisons the first row makes from 3 to 1,
        -- since the diagonal and leftward values are "out of bounds."
        leftColumn _leftElement !insertCost !i !j =
            {-# SCC leftColumn #-}
            do
                firstPrevCost ← readCost (i - 1) j
                write (i, j) (# insertCost + firstPrevCost, UpArrow #)

        -- Define how to compute the first cell of the remaining rows.
        -- We need to ensure that there are no Left Arrow values in the directional matrix.
        -- We can also reduce the number of comparisons the first row makes from 3 to 2,
        -- since the leftward values are "out of bounds."
        -- Define how to compute the first cell of the remaining rows.
        -- We need to ensure that there are no Left Arrow values in the directional matrix.
        -- We can also reduce the number of comparisons the first row makes from 3 to 2,
        -- since the leftward values are "out of bounds."
        leftBoundary !leftElement !insertCost !i !j =
            {-# SCC leftBoundary #-}
            let topElement = longerTop ! (j - 1)
                alignCost = costλ leftElement topElement
            in  do
                    diagCost ← readCost (i - 1) $ j - 1
                    topCost ← readCost (i - 1) j
                    let v =
                            minimumCostDirection
                                maxBound
                                (alignCost + diagCost)
                                (insertCost + topCost)
                    write (i, j) v

        -- Define how to compute the last cell of the first "rows - offset" rows.
        -- We need to ensure that there are only Left Arrow values in the directional matrix.
        -- We can also reduce the number of comparisons the first row makes from 3 to 1,
        -- since the diagonal and upward values are "out of bounds."
        rightBoundary !leftElement _insertCost !i !j =
            {-# SCC rightBoundary #-}
            let topElement = longerTop ! (j - 1)
                deleteCost = costλ gap topElement
                alignCost = costλ leftElement topElement
            in  do
                    diagCost ← readCost (i - 1) $ j - 1
                    leftCost ← readCost i $ j - 1
                    let v =
                            minimumCostDirection
                                (deleteCost + leftCost)
                                (alignCost + diagCost)
                                maxBound
                    write (i, j) v

        rightColumn = {-# SCC rightColumn #-} internalCell


{- |
Produces a set of reusable values and  functions which are "constant" between
different incarnations of the Ukkonen algorithms.
-}
ukkonenConstants
    ∷ (Vector v e)
    ⇒ TCMλ e
    -- ^ Metric between states producing the medoid of states.
    → v e
    -- ^ Shorter dynamic character
    → v e
    -- ^ Longer  dynamic character
    → Word
    -- ^ Current  offset from quasi-diagonal
    → (Int, e → e → Word, Int, Int, Int, Int)
ukkonenConstants tcmλ lesserLeft longerTop co =
    (offset, costλ, rows, cols, width, quasiDiagonalWidth)
    where
        offset = clampOffset co
        costλ x = snd . tcmλ x
        longerLen = GV.length longerTop
        lesserLen = GV.length lesserLeft
        rows = GV.length lesserLeft + 1
        cols = GV.length longerTop + 1
        width = quasiDiagonalWidth + (offset `shiftL` 1) -- Multiply by 2
        quasiDiagonalWidth = differenceInLength + 1
            where
                differenceInLength = longerLen - lesserLen

        -- Note: "offset" cannot cause "width + quasiDiagonalWidth" to exceed "2 * cols"
        clampOffset o =
            let o' = fromEnum o in min o' $ cols - quasiDiagonalWidth
