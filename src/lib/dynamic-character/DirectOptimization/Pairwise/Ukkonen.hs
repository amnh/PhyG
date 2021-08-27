-----------------------------------------------------------------------------
-- |
-- Module      :  DirectOptimization.Pairwise.Ukkonen
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-- Direct optimization pairwise alignment using the Needleman-Wunsch algorithm.
-- These functions will allocate an M * N matrix.
--
-----------------------------------------------------------------------------

{-# LANGUAGE ApplicativeDo       #-}
{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE ConstraintKinds     #-}
{-# LANGUAGE DerivingStrategies  #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies        #-}

module DirectOptimization.Pairwise.Ukkonen
  ( ukkonenDO
  ) where

import           Bio.DynamicCharacter
import           Control.Monad                        (unless, when)
import           Control.Monad.Loops                  (iterateUntilM, whileM_)
import           Control.Monad.Primitive
import           Control.Monad.ST
import           Data.Bits
import           Data.Foldable
import           Data.Matrix.Unboxed                  (Matrix, unsafeFreeze)
import           Data.Matrix.Unboxed.Mutable          (MMatrix)
import qualified Data.Matrix.Unboxed.Mutable          as M
import           Data.STRef
import qualified Data.Vector                          as V
import           Data.Vector.Generic                  (Vector, (!))
import qualified Data.Vector.Generic                  as GV
import qualified Data.Vector.Unboxed                  as UV
import           DirectOptimization.Pairwise.Internal
import           DirectOptimization.Pairwise.Swapping


-- |
-- Performs a naive direct optimization.
-- Takes in two characters to run DO on and an overlap function
-- Returns an assignment character, the cost of that assignment, the assignment
-- character with gaps included, the aligned version of the first input character,
-- and the aligned version of the second input character. The process for this
-- algorithm is to generate a traversal matrix, then perform a traceback.
{-# SCC        ukkonenDO #-}
{-# INLINEABLE ukkonenDO #-}
{-# SPECIALISE ukkonenDO :: Word -> WideState -> (WideState -> WideState -> (WideState, Word)) -> WideDynamicCharacter -> WideDynamicCharacter -> (Word, WideDynamicCharacter) #-}
{-# SPECIALISE ukkonenDO :: Word -> HugeState -> (HugeState -> HugeState -> (HugeState, Word)) -> HugeDynamicCharacter -> HugeDynamicCharacter -> (Word, HugeDynamicCharacter) #-}
ukkonenDO
  :: ( FiniteBits a
     , Ord a
     , Ord (v a)
     , Vector v a
     , Vector v (a, a, a)
     )
  => Word                    -- ^ Coefficent
  -> a                       -- ^ Gap state
  -> (a -> a -> (a, Word))   -- ^ Metric between states
  -> (v a, v a, v a)         -- ^ /1st/ dynamic character
  -> (v a, v a, v a)         -- ^ /2nd/ dynamic character
  -> (Word, (v a, v a, v a))
ukkonenDO coefficient gap overlapλ char1 char2
  | noGainFromUkkonenMethod = buildFullMatrix
  | otherwise               = directOptimization (const buildPartialMatrixMaybe) gap overlapλ char1 char2
  where
    ~(_, _, _, (lesser,_,_), (longer,_,_)) = measureAndUngapCharacters gap char1 char2

    buildFullMatrix = swappingDO gap overlapλ char1 char2

--  buildPartialMatrixMaybe :: Vector v a => (a -> a -> (a, Word)) -> v a -> v a -> (Word, Matrix Direction)
    buildPartialMatrixMaybe = createUkkonenMethodMatrix coefficient gapsPresentInInputs gap

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
    --     OR
    -- > gapsPresentInInputs >= longerLen
    noGainFromUkkonenMethod =     lesserLen <= 4
                           || 2 * longerLen >= 3 * lesserLen
                           || coefficient == 0
                           || gapsPresentInInputs >= longerLen
      where
        longerLen = toEnum $ GV.length longer
        lesserLen = toEnum $ GV.length lesser

    -- /O(n + m)/
    --
    -- If one or more of the aligned character elements contained a gap, diagonal
    -- directions in the matrix have an "indel" cost. 'gapsPresentInInputs' is
    -- necessary in order to decrement the threshold value to account for this.
    -- This was not described in Ukkonen's original paper, as the inputs were assumed
    -- not to contain any gaps.
    gapsPresentInInputs = char1Gaps + char2Gaps
      where
        char1Gaps = toEnum $ countGaps char1
        char2Gaps = toEnum $ countGaps char2
        countGaps (x,_,_) = GV.length $ GV.filter hasGap x
        hasGap b  = popCount (b .&. gap) > 0


-- |
-- /O( (n - m + 1 ) * log(n - m + 1) )/, /n/ >= /m/
--
-- Generates an /optimal/, partially-filled-in matrix using Ukkonen's string
-- edit distance algorithm.
--
-- Note that the threshold value is lowered more than described in Ukkonen's
-- paper. This is to handle input elements that contain a gap. In Ukkonen's
-- original description of the algorithm, there was a subtle assumption that
-- input did not contain any gap symbols.
{-# SCC        createUkkonenMethodMatrix #-}
{-# INLINEABLE createUkkonenMethodMatrix #-}
{-# SPECIALISE createUkkonenMethodMatrix :: Word -> Word -> WideState -> (WideState -> WideState -> (WideState, Word)) -> UV.Vector WideState -> UV.Vector WideState -> (Word, Matrix Direction) #-}
{-# SPECIALISE createUkkonenMethodMatrix :: Word -> Word -> HugeState -> (HugeState -> HugeState -> (HugeState, Word)) ->  V.Vector HugeState ->  V.Vector HugeState -> (Word, Matrix Direction) #-}
createUkkonenMethodMatrix
  :: ( Eq a
     , Vector v a
     )
  => Word   -- ^ Coefficient value, representing the /minimum/ transition cost from a state to gap
  -> Word   -- ^ Gaps present in input
  -> a
  -> (a -> a -> (a, Word))
  -> v a -- ^ Longer dynamic character
  -> v a -- ^ Shorter dynamic character
--  -> (v a, v a, v a) -- ^ Longer dynamic character
--  -> (v a, v a, v a) -- ^ Shorter dynamic character
  -> (Word, Matrix Direction)
createUkkonenMethodMatrix minimumIndelCost gapsPresentInInputs gap overlapλ longerTop lesserLeft = finalMatrix
  where
    -- General values that need to be in scope for the recursive computations.
    longerLen = GV.length longerTop
    lesserLen = GV.length lesserLeft

    -- We start the offset at two rather than at one so that the first doubling
    -- isn't trivially small.
    startOffset = 2 + gapsPresentInInputs

    -- /O(1)/
    --
    -- Necessary to compute the width of a row in the barrier-constrained matrix.
    quasiDiagonalWidth = toEnum $ differenceInLength + 1
      where
        differenceInLength = longerLen - lesserLen

    needToResizeBand :: forall s. MMatrix s Word -> STRef s Word -> ST s Bool
    needToResizeBand mCost offsetRef = do
        offset        <- readSTRef offsetRef
        if   quasiDiagonalWidth + offset > toEnum longerLen
        then pure False
        else do
                alignmentCost <- M.unsafeRead mCost (lesserLen, longerLen)
                let threshold -- The threshold value must be non-negative
                      | quasiDiagonalWidth + offset <= gapsPresentInInputs = 0
                      | otherwise = minimumIndelCost * (quasiDiagonalWidth + offset - gapsPresentInInputs)
                pure $ threshold <= alignmentCost

    finalMatrix = runST $ do
        (mCost, mDir) <- buildInitialBandedMatrix gap overlapλ longerTop lesserLeft startOffset
        offsetRef <- newSTRef startOffset
        whileM_ (needToResizeBand mCost offsetRef) $ do
          previousOffset <- readSTRef offsetRef
          let currentOffset = previousOffset `shiftL` 1 -- Multiply by 2
          writeSTRef offsetRef currentOffset
          expandBandedMatrix gap overlapλ longerTop lesserLeft mCost mDir previousOffset currentOffset

        c <- M.unsafeRead mCost (lesserLen, longerLen)
        m <- unsafeFreeze mDir
        pure (c, m)


{-# SCC buildInitialBandedMatrix #-}
buildInitialBandedMatrix
  :: ( Eq a
     , Vector v a
     )
  => a   -- ^ Gap
  -> (a -> a -> (a, Word))
  -> v a -- ^ Longer dynamic character
  -> v a -- ^ Shorter dynamic character
  -> Word
  -> ST s (MMatrix s Word, MMatrix s Direction)
buildInitialBandedMatrix gap overlapλ longerTop lesserLeft o = fullMatrix
  where
    (offset, cost, rows, cols, width, quasiDiagonalWidth) = ukkonenConstants overlapλ lesserLeft longerTop o

    fullMatrix = do

      ---------------------------------------
      -- Allocate required space           --
      ---------------------------------------

      mCost <- M.new (rows, cols)
      mDir  <- M.new (rows, cols)

      ---------------------------------------
      -- Define some generalized functions --
      ---------------------------------------
      let (write, internalCell, leftColumn, leftBoundary, rightBoundary, rightColumn) = edgeCellDefinitions mCost mDir overlapλ cost longerTop gap

      -- Define how to compute values to an entire row of the Ukkonen matrix.
      let writeRow i =
            -- Precomute some values that will be used for the whole row
            let start = max  0         $ i - offset
                stop  = min (cols - 1) $ i - offset + width - 1
                leftElement = lesserLeft ! (i - 1)
                insertCost  = cost gap leftElement
                firstCell
                  | i <= offset = leftColumn
                  | otherwise   = leftBoundary

                lastCell
                  | i <= cols - quasiDiagonalWidth - offset = rightBoundary
                  | otherwise = rightColumn
            in  do -- Write to the first cell of the Ukkonen band
                   firstCell leftElement insertCost i start >>= write (i, start)
                   -- Write to the all the intermediary cells of the Ukkonen band
                   for_ [start + 1 .. stop - 1] $ \j ->
                       internalCell leftElement insertCost i j >>= write (i, j)
                   -- Write to the last cell of the Ukkonen band
                   lastCell leftElement insertCost i stop >>= write (i, stop)

      ---------------------------------------
      -- Compute all values of the matrix  --
      ---------------------------------------

      -- Write to the origin to seed the first row.
      write (0, 0) (0, DiagArrow)

      -- Write the first row to seed subsequent rows.
      for_ [1 .. min (cols - 1) (width - offset - 1)] $ \j ->
        let topElement    = longerTop ! (j - 1)
            firstCellCost = cost gap topElement
        in  do firstPrevCost <- M.unsafeRead mCost (0, j - 1)
               write (0,j) (firstCellCost + firstPrevCost, LeftArrow)

      -- Loop through the remaining rows.
      for_ [1 .. rows - 1] writeRow

      -- Return the matricies for possible expansion
      pure (mCost, mDir)


{-# SCC expandBandedMatrix #-}
-- |
-- Given a partially computed alignment matrix,
-- will expand the computed region to the new specified offset.
--
--
-- Dimensions: 13 ⨉ 17
--  ⊗ ┃  ⁎ α1 α2 α3 α4 α5 α6 α7 α8 α9 α0 α1 α2 α3 α4 α5 α6
-- ━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
--  ⁎ ┃ 0↖ 0← 0← 0← 0← 0← 0← 0← 0← 0← 0← 0← 0← 0← 0← 0← 0←
-- α1 ┃ 0↑ 0↖ 0↖ 0↖ 0← 0← 0← 0← 0← 0← 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0←
-- α2 ┃ 0↑ 0↖ 0↖ 0↖ 0← 0← 0← 0← 0← 0← 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0←
-- α3 ┃ 0↑ 0↑ 0↑ 0↖ 0↖ 0↖ 0← 0← 0← 0← 0← 0← 0← 0← 0← 0← 0↖
-- α4 ┃ 0↑ 0↑ 0↑ 0↖ 0↖ 0← 0↖ 0← 0↖ 0← 0← 0← 0← 0← 0← 0← 0←
-- α5 ┃ 0↑ 0↑ 0↑ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0← 0← 0← 0← 0← 0← 0← 0←
-- α6 ┃ 0↑ 0↖ 0↖ 0↖ 0↑ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0←
-- α7 ┃ 0↑ 0↑ 0↑ 0↑ 0↖ 0↖ 0↖ 0↖ 0↖ 0← 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖
-- α8 ┃ 0↑ 0↑ 0↑ 0↑ 0↑ 0↖ 0← 0↖ 0↑ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖
-- α9 ┃ 0↑ 0↑ 0↑ 0↑ 0↖ 0↑ 0↖ 0← 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖
-- α0 ┃ 0↑ 0↑ 0↑ 0↑ 0↑ 0↖ 0↑ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖
-- α1 ┃ 0↑ 0↑ 0↑ 0↑ 0↖ 0↑ 0↖ 0↖ 0↖ 0← 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖
-- α2 ┃ 0↑ 0↖ 0↖ 0↖ 0↑ 0↑ 0↑ 0↖ 0↑ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖ 0↖
--
--      ┌───────────────w───────────────┐
--      │              ┏━━━━━━━co━━━━━━━┪
--      ┢━━━━━qd━━━━━━┓┠─po─┐┌────Δo────┨
--  ⊗ ┃ ┃0  1  2  3  4┃┃5  6││7  8  9 10┃11 12 13 14 15 16
-- ━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
--  0 ┃ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒ ▒▒ ▒▒ ▒▒
--  1 ┃ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒ ▒▒ ▒▒ ▒▒
--  2 ┃ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒ ▒▒ ▒▒ ▒▒
--  3 ┃ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒ ▒▒ ▒▒ ▒▒
--  4 ┃ ▒▒ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒ ▒▒ ▒▒ ▒▒
--  5 ┃ ▒▒ ▒▒ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒ ▒▒ ▒▒ ▒▒
--  6 ┃ ▒▒ ▒▒ ▒▒ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒ ▒▒ ▒▒ ▒▒
--  7 ┃    ▒▒ ▒▒ ▒▒ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒ ▒▒ ▒▒
--  8 ┃       ▒▒ ▒▒ ▒▒ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒ ▒▒
--  9 ┃          ▒▒ ▒▒ ▒▒ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓ ▒▒
--  0 ┃             ▒▒ ▒▒ ▒▒ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓ ▓▓
--  1 ┃                ▒▒ ▒▒ ▒▒ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██ ▓▓
--  2 ┃                   ▒▒ ▒▒ ▒▒ ▒▒ ▓▓ ▓▓ ██ ██ ██ ██ ██
--
--
-- w  : Width
-- qd : Quasi-diagonal
-- co : Current Offset
-- po : Previous Offset
-- Δo : Difference in Offset
--
-- Note:
-- w  = qd + co
-- co = po + Δo
--
-- And often:
-- co = 2*po = 2*Δo
--
-- ██ : The core band
--       * Previously computed, sections may need to be recomputed
--
-- ▓▓ : The previous extension
--       * Previously computed, sections may need to be recomputed
--
-- ▒▒ : The new extension
--       * Needs to be computed
--
expandBandedMatrix
  :: ( Eq a
     , Vector v a
     )
  => a
  -> (a -> a -> (a, Word))
  -> v a -- ^ Longer dynamic character
  -> v a -- ^ Shorter dynamic character
  -> MMatrix s Word
  -> MMatrix s Direction
  -> Word
  -> Word
  -> ST s ()
expandBandedMatrix gap overlapλ longerTop lesserLeft mCost mDir po co = updatedBand
  where
    (offset, cost, rows, cols, width, qd) = ukkonenConstants overlapλ lesserLeft longerTop co
    prevOffset  = fromEnum po

    updatedBand = do

      ---------------------------------------
      -- Allocate mutable state variables  --
      ---------------------------------------

      tailStart <- newSTRef cols

      t0' <- newSTRef (-1)
      t1' <- newSTRef $ qd + fromEnum po

      ---------------------------------------
      -- Define some generalized functions --
      ---------------------------------------
      let (write, internalCell, leftColumn, leftBoundary, rightBoundary, rightColumn) = edgeCellDefinitions mCost mDir overlapλ cost longerTop gap

      let computeCell leftElement insertCost i j = {-# SCC recomputeCell #-}
            let topElement = longerTop ! (j - 1)
                deleteCost = cost topElement    gap
                (alignElem, alignCost) = overlapλ topElement leftElement
            in do
                  diagCost <- M.unsafeRead mCost (i - 1, j - 1)
                  topCost  <- M.unsafeRead mCost (i - 1, j    )
                  leftCost <- M.unsafeRead mCost (i    , j - 1)
                  oldCost  <- M.unsafeRead mCost (i    , j    )
                  let e@(c,_) = getMinimalResult gap alignElem
                                  [ ( alignCost + diagCost, DiagArrow)
                                  , (deleteCost + leftCost, LeftArrow)
                                  , (insertCost +  topCost, UpArrow  )
                                  ]
                  write (i,j) e
                  pure (c == oldCost, j+1)
--                  pure (c /= oldCost, j+1)

      let recomputeRange leftElement insertCost i x y = do
            lastDiff <- newSTRef 0
            for_ [x .. y] $ \j -> do
              (same, _) <- computeCell leftElement insertCost i j
              unless same $ writeSTRef lastDiff j
            readSTRef lastDiff

      -- Define how to compute values to an entire row of the Ukkonen matrix.
      let extendRow i =
            -- Precomute some values that will be used for the whole row
            let start0 =  max  0         $ i - offset
                start3 =  min  cols      $ i + width - offset - prevOffset - 1
                goUpTo =  max  0         ( i - prevOffset) - 1
                stop   =  min (cols - 1) $ i + width - offset - 1
                leftElement = lesserLeft ! (i - 1)
                insertCost  = cost gap leftElement
                firstCell
                  | i <= offset = leftColumn
                  | otherwise   = leftBoundary

                lastCell
                  | i <= cols - qd - offset = rightBoundary
                  | otherwise = rightColumn

                b0 = start0
                e0 = goUpTo
                b1 = start3
                e1 = stop


                continueRecomputing (same, j) = same || j >= stop
                computeCell' ~(_,j) = computeCell leftElement insertCost i j
                internalCell' j = internalCell leftElement insertCost i j >>= write (i,j)
                recomputeUntilSame j = snd <$> iterateUntilM continueRecomputing computeCell' (False, j)
            in  do -- First, we fill in 0 or more cells of the left region of
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
                   if   i > prevOffset
                   then firstCell leftElement insertCost i start0 >>= write (i, b0)
                   else pure ()

                   for_ [b0+1 .. e0] internalCell'

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
                   s0 <- (\x -> min (x+1) e1) <$> readSTRef t0'
                   writeSTRef t0' (-1)

                   when (s0 > e0 && toEnum i > po) $
                       recomputeRange leftElement insertCost i (e0+1) s0 >>= writeSTRef t0'
                   t0 <- readSTRef t0'

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
                   if      s0 == t0 && s0 > 0
                   then recomputeUntilSame (s0 + 1) >>= writeSTRef t0' . pred
                   else if s0 <= e0 && e0 > 0
                   then recomputeUntilSame (e0 + 1) >>= writeSTRef t0' . pred
                   else pure ()

{-
                   headStop' <- if   leadStop >= start1
                                then pure leadStop
                                else recomputeUntilSame start1
-}
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
                   s1 <- readSTRef t1'

                   t1 <- recomputeRange leftElement insertCost i s1 $ b1 - 1

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
                   if   i < rows - fromEnum po
                   then lastCell leftElement insertCost i stop >>= write (i, stop)
                   else pure ()

                   -- Update references for the next row
--                   writeSTRef headStop headStop'
--                   writeSTRef tailStart $ if tailStop' /= start2 then tailStop' else start3

      ---------------------------------------
      -- Compute all values of the matrix  --
      ---------------------------------------

      let start = qd + prevOffset

      -- Extend the first row to seed subsequent rows.
      for_ [start .. min (cols - 1) (width - offset - 1)] $ \j ->
        let topElement    = longerTop ! (j - 1)
            firstCellCost = cost gap topElement
        in  do firstPrevCost <- M.unsafeRead mCost (0, j - 1)
               write (0,j) (firstCellCost + firstPrevCost, LeftArrow)

      writeSTRef tailStart start

      -- Loop through the remaining rows.
      for_ [1 .. rows - 1] extendRow


getMinimalResult
  :: ( Eq a
     , Foldable f
     )
  => a
  -> a
  -> f (Word, Direction)
  -> (Word, Direction)
getMinimalResult gap alignElem opts =
    let v@(~(c,d)) = minimum opts
    in  if   d == DiagArrow && alignElem == gap
        then (c, LeftArrow)
        else v


edgeCellDefinitions
  :: ( Eq a
     , PrimMonad m
     , Vector v a
     )
  => MMatrix (PrimState m) Word
  -> MMatrix (PrimState m) Direction
  -> (a -> a -> (a, Word))
  -> (a -> a -> Word)
  -> v a -- ^ Longer dynamic character
  -> a
  -> ( (Int, Int) -> (Word, Direction) -> m ()
     , a -> Word -> Int -> Int -> m (Word, Direction)
     , a -> Word -> Int -> Int -> m (Word, Direction)
     , a -> Word -> Int -> Int -> m (Word, Direction)
     , a -> Word -> Int -> Int -> m (Word, Direction)
     , a -> Word -> Int -> Int -> m (Word, Direction)
     )
edgeCellDefinitions mCost mDir overlapλ cost longerTop gap = (write, internalCell, leftColumn, leftBoundary, rightBoundary, rightColumn)
  where
    -- Write to a single cell of the current vector and directional matrix simultaneously
    write !p ~(!c, !d) = M.unsafeWrite mCost p c *> M.unsafeWrite mDir p d

    -- Write to an internal cell (not on a boundary) of the matrix.
    internalCell leftElement insertCost i j = {-# SCC internalCell_expanding #-}
        let topElement = longerTop ! (j - 1)
            -- Preserve the gap in the top (longer) string
        in  if topElement == gap
            then (\x -> (x, LeftArrow)) <$> M.unsafeRead mCost (i, j - 1)
            -- Normal Needleman-Wunsch Logic
            else let  deleteCost = cost topElement    gap
                      (alignElem, alignCost) = overlapλ topElement leftElement
                 in  do diagCost <- M.unsafeRead mCost (i - 1, j - 1)
                        topCost  <- M.unsafeRead mCost (i - 1, j    )
                        leftCost <- M.unsafeRead mCost (i    , j - 1)
                        pure $ getMinimalResult gap alignElem
                            [ ( alignCost + diagCost, DiagArrow)
                            , (deleteCost + leftCost, LeftArrow)
                            , (insertCost +  topCost, UpArrow  )
                            ]

    -- Define how to compute the first cell of the first "offset" rows.
    -- We need to ensure that there are only Up Arrow values in the directional matrix.
    -- We can also reduce the number of comparisons the first row makes from 3 to 1,
    -- since the diagonal and leftward values are "out of bounds."
    leftColumn _leftElement insertCost i j = {-# SCC leftColumn #-} do
        firstPrevCost <- M.unsafeRead mCost (i - 1, j)
        pure (insertCost + firstPrevCost, UpArrow)

    -- Define how to compute the first cell of the remaining rows.
    -- We need to ensure that there are no Left Arrow values in the directional matrix.
    -- We can also reduce the number of comparisons the first row makes from 3 to 2,
    -- since the leftward values are "out of bounds."
    -- Define how to compute the first cell of the remaining rows.
    -- We need to ensure that there are no Left Arrow values in the directional matrix.
    -- We can also reduce the number of comparisons the first row makes from 3 to 2,
    -- since the leftward values are "out of bounds."
    leftBoundary leftElement insertCost i j = {-# SCC leftBoundary #-}
        let topElement = longerTop ! (j - 1)
            (alignElem, alignCost) = overlapλ topElement leftElement
        in  do diagCost <- M.unsafeRead mCost (i - 1, j - 1)
               topCost  <- M.unsafeRead mCost (i - 1, j    )
               pure $ getMinimalResult gap alignElem
                   [ ( alignCost + diagCost, DiagArrow)
                   , (insertCost +  topCost, UpArrow  )
                   ]

    -- Define how to compute the last cell of the first "rows - offset" rows.
    -- We need to ensure that there are only Left Arrow values in the directional matrix.
    -- We can also reduce the number of comparisons the first row makes from 3 to 1,
    -- since the diagonal and upward values are "out of bounds."
    rightBoundary leftElement _insertCost i j = {-# SCC rightBoundary #-}
        let topElement = longerTop ! (j - 1)
            deleteCost = cost topElement    gap
            (alignElem, alignCost) = overlapλ topElement leftElement
        in  do diagCost <- M.unsafeRead mCost (i - 1, j - 1)
               leftCost <- M.unsafeRead mCost (i    , j - 1)
               pure $ getMinimalResult gap alignElem
                   [ ( alignCost + diagCost, DiagArrow)
                   , (deleteCost + leftCost, LeftArrow)
                   ]

    rightColumn = {-# SCC rightColumn #-} internalCell


-- |
-- Produces a set of reusable values and  functions which are "constant" between
-- different incarnations of the Ukkonen algorithms.
ukkonenConstants
  :: ( Vector v a
     )
  => (a -> a -> (a, Word))
  -> v a -- ^ Longer dynamic character
  -> v a -- ^ Shorter dynamic character
  -> Word
  -> (Int, a -> a -> Word, Int, Int, Int, Int)
ukkonenConstants overlapλ lesserLeft longerTop o =
    (offset, cost, rows, cols, width, quasiDiagonalWidth)
  where
    -- Note: "offset" cannot cause "width + quasiDiagonalWidth" to exceed "2 * cols"
    offset      = let o' = fromEnum o in  min o' $ cols - quasiDiagonalWidth
    cost x y    = snd $ overlapλ x y
    longerLen   = GV.length longerTop
    lesserLen   = GV.length lesserLeft
    rows        = GV.length lesserLeft + 1
    cols        = GV.length longerTop  + 1
    width       = quasiDiagonalWidth + (offset `shiftL` 1) -- Multiply by 2
    quasiDiagonalWidth = differenceInLength + 1
      where
        differenceInLength = longerLen - lesserLen
