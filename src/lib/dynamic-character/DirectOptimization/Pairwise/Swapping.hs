-----------------------------------------------------------------------------
-- |
-- Module      :  DirectOptimization.Pairwise.Swapping
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

{-# LANGUAGE ApplicativeDo    #-}
{-# LANGUAGE ConstraintKinds  #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE Strict           #-}
{-# LANGUAGE TypeFamilies     #-}

module DirectOptimization.Pairwise.Swapping
  ( swappingDO
  ) where

import           Bio.DynamicCharacter
import           Control.Monad.ST
import           Data.Bits
import           Data.Foldable
import           Data.Matrix.Unboxed                  (Matrix, unsafeFreeze)
import qualified Data.Matrix.Unboxed.Mutable          as M
import qualified Data.Vector                          as V
import           Data.Vector.Generic                  (Vector, (!))
import qualified Data.Vector.Generic                  as GV
import qualified Data.Vector.Unboxed                  as UV
import qualified Data.Vector.Unboxed.Mutable          as MUV
import           DirectOptimization.Pairwise.Internal


-- |
-- Performs a naive direct optimization.
-- Takes in two characters to run DO on and an overlap function
-- Returns an assignment character, the cost of that assignment, the assignment
-- character with gaps included, the aligned version of the first input character,
-- and the aligned version of the second input character. The process for this
-- algorithm is to generate a traversal matrix, then perform a traceback.
{-# SCC        swappingDO #-}
{-# INLINEABLE swappingDO #-}
{-# SPECIALISE swappingDO :: WideState -> (WideState -> WideState -> (WideState, Word)) -> WideDynamicCharacter -> WideDynamicCharacter -> (Word, WideDynamicCharacter) #-}
{-# SPECIALISE swappingDO :: HugeState -> (HugeState -> HugeState -> (HugeState, Word)) -> HugeDynamicCharacter -> HugeDynamicCharacter -> (Word, HugeDynamicCharacter) #-}
swappingDO
  :: ( FiniteBits a
     , Ord a
     , Ord (v a)
     , Vector v a
     , Vector v (a, a, a)
     )
  => a
  -> (a -> a -> (a, Word))
  -> (v a, v a, v a)
  -> (v a, v a, v a)
  -> (Word, (v a, v a, v a))
swappingDO = directOptimization buildDirectionMatrix


{-# SCC        buildDirectionMatrix #-}
{-# INLINEABLE buildDirectionMatrix #-}
{-# SPECIALISE buildDirectionMatrix :: WideState -> (WideState -> WideState -> (WideState, Word)) -> UV.Vector WideState -> UV.Vector WideState -> (Word, Matrix Direction) #-}
{-# SPECIALISE buildDirectionMatrix :: HugeState -> (HugeState -> HugeState -> (HugeState, Word)) ->  V.Vector HugeState ->  V.Vector HugeState -> (Word, Matrix Direction) #-}
buildDirectionMatrix
  :: ( Vector v a
     )
  => a
  -> (a -> a -> (a, Word))
  -> v a
  -> v a
  -> (Word, Matrix Direction)
buildDirectionMatrix gap overlapλ topChar leftChar = fullMatrix
  where
    cost x y = snd $ overlapλ x y
    rows     = GV.length leftChar + 1
    cols     = GV.length topChar  + 1

    fullMatrix = runST $ do
      mDir <- M.new (rows, cols)
      vOne <- MUV.new cols
      vTwo <- MUV.new cols

      let write v p@(_,j) c d = MUV.unsafeWrite v j c *> M.unsafeWrite mDir p d

      write vOne (0,0) 0 DiagArrow

      -- Special case the first row
      -- We need to ensure that there are only Left Arrow values in the directional matrix.
      -- We can also reduce the number of comparisons the first row makes from 3 to 1,
      -- since the diagonal and upward values are "out of bounds."
      for_ [1 .. cols - 1] $ \j ->
        let topElement    = topChar ! (j - 1)
            firstCellCost = cost gap topElement
        in  do firstPrevCost <- MUV.unsafeRead vOne (j - 1)
               write vOne (0,j) (firstCellCost + firstPrevCost) LeftArrow

      for_ [1 .. rows - 1] $ \i ->
        let (prev, curr)
              | odd i     = (vOne, vTwo)
              | otherwise = (vTwo, vOne)
            leftElement   = leftChar ! (i - 1)
            -- Special case the first cell of each row
            -- We need to ensure that there are only Up Arrow values in the directional matrix.
            -- We can also reduce the number of comparisons the first row makes from 3 to 1,
            -- since the diagonal and leftward values are "out of bounds."
            firstCellCost = cost leftElement gap
        in  do firstPrevCost <- MUV.unsafeRead prev 0
               write curr (i,0) (firstCellCost + firstPrevCost) UpArrow
               -- Finish special case for first cell of each row
               -- Begin processing all other cells in the curr vector
               for_ [1 .. cols - 1] $ \j ->
                 let topElement  = topChar ! (j - 1)
                     deleteCost  = cost topElement gap
                     alignCost   = cost topElement leftElement
                     insertCost  = cost gap        leftElement
                 in  do diagCost <- MUV.unsafeRead prev $ j - 1
                        topCost  <- MUV.unsafeRead prev   j
                        leftCost <- MUV.unsafeRead curr $ j - 1
                        let xs = [ ( alignCost + diagCost, DiagArrow)
                                 , (deleteCost + leftCost, LeftArrow)
                                 , (insertCost +  topCost, UpArrow  )
                                 ]
                        let (c,d) = minimum xs
                        write curr (i,j) c d

      let v | odd  rows = vOne
            | otherwise = vTwo
      c <- MUV.unsafeRead v (cols - 1)
      m <- unsafeFreeze mDir
      pure (c, m)
