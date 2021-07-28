-----------------------------------------------------------------------------
-- |
-- Module      :  Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Huge.Swapping
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
{-# LANGUAGE UnboxedTuples    #-}

module Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Huge.Swapping
  ( hugeSwappingDO
  , directOptimization
  , traceback
  ) where

import           Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.DynamicCharacter2
import           Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Huge.Internal
import           Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Internal          (Direction(..))
import           Control.Monad.ST
import           Data.BitVector.LittleEndian
import           Data.Bits
import           Data.DList                                                               (snoc)
import           Data.Foldable
import           Data.Matrix.Unboxed                                                      (Matrix, unsafeFreeze, unsafeIndex)
import qualified Data.Matrix.Unboxed.Mutable                                              as M
import           Data.MonoTraversable
import           Data.Vector                                                              (Vector, (!))
import qualified Data.Vector                                                              as V
import qualified Data.Vector.Unboxed.Mutable                                              as MUV


-- |
-- Performs a naive direct optimization.
-- Takes in two characters to run DO on and an overlap function
-- Returns an assignment character, the cost of that assignment, the assignment
-- character with gaps included, the aligned version of the first input character,
-- and the aligned version of the second input character. The process for this
-- algorithm is to generate a traversal matrix, then perform a traceback.
{-# SCC        hugeSwappingDO #-}
{-# INLINEABLE hugeSwappingDO #-}
hugeSwappingDO
  :: (BitVector -> BitVector -> (BitVector, Word))
  -> HugeDynamicCharacter
  -> HugeDynamicCharacter
  -> (Word, HugeDynamicCharacter)
hugeSwappingDO = directOptimization buildDirectionMatrix


buildDirectionMatrix
  :: (BitVector -> BitVector -> (BitVector, Word))
  -> Vector BitVector
  -> Vector BitVector
  -> (Word, Matrix Direction)
buildDirectionMatrix overlapλ topChar leftChar = fullMatrix
  where
    symbolCount = dimension $ V.head topChar
    cost x y    = snd $ overlapλ x y
    gap         = bit . fromEnum $ symbolCount - 1
    rows        = olength leftChar + 1
    cols        = olength topChar  + 1

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


{-# SCC directOptimization #-}
directOptimization
  :: ((BitVector -> BitVector -> (BitVector, Word)) -> Vector BitVector -> Vector BitVector -> (Word, Matrix Direction))
  -> (BitVector -> BitVector -> (BitVector, Word))
  -> HugeDynamicCharacter -- ^ Longer dynamic character
  -> HugeDynamicCharacter -- ^ Shorter dynamic character
  -> (Word, HugeDynamicCharacter)
directOptimization matrixFunction overlapλ lhs rhs =
    let ~(swapped, gap, gapsLesser, gapsLonger, (shorterChar,_,_), (longerChar,_,_)) = measureAndUngapCharacters lhs rhs
        ~(alignmentCost, ungappedAlignment) =
            if      olength shorterChar == 0
            then if olength  longerChar == 0
                 -- Niether character was Missing, but both are empty when gaps are removed
                 then (0, (mempty,mempty,mempty))
                 -- Niether character was Missing, but one of them is empty when gaps are removed
                 else let vec = V.generate (V.length longerChar) $ \i -> fst (overlapλ (longerChar ! i) gap)
                      in  (0, (vec, V.replicate (V.length longerChar) (gap `xor` gap), longerChar))
                 -- Both have some non-gap elements, perform string alignment
            else let (cost, traversalMatrix) = matrixFunction overlapλ longerChar shorterChar
                 in  (cost, traceback gap overlapλ traversalMatrix longerChar shorterChar)
        transformation    = if swapped then \(m,l,r) -> (m,r,l) else id
        regappedAlignment = insertGaps gap gapsLesser gapsLonger ungappedAlignment
        alignmentContext  = transformation regappedAlignment
    in  (alignmentCost, alignmentContext)


{-# SCC traceback #-}
traceback
  :: BitVector
  -> (BitVector -> BitVector -> (BitVector, Word))
  -> Matrix Direction
  -> Vector BitVector     -- ^ Longer dynamic character
  -> Vector BitVector     -- ^ Shorter dynamic character
  -> HugeDynamicCharacter -- ^ Resulting dynamic character
traceback gap overlapλ alignMatrix longerChar lesserChar = alignmentContext
  where
    f x y = fst $ overlapλ x y

    alignmentContext = dlistToDynamic $ go startPoint
    dlistToDynamic   = V.unzip3 . V.fromList . toList

    longerLen = olength longerChar
    lesserLen = olength lesserChar
    rows      = lesserLen + 1
    cols      = longerLen + 1
    zero      = gap `xor` gap

    startPoint = (rows - 1, cols - 1)

    go p@(~(i, j))
      | p == (0,0) = mempty
      | otherwise  =
        let previousSequence = go (row', col')

            directionArrow = unsafeIndex alignMatrix p

            (# row', col', localContext #) =
                case directionArrow of
                  LeftArrow -> let j' = j - 1
                                   y  = longerChar ! j'
                                   e  = (f gap y, zero,    y)
                               in  (# i , j', e #)
                  UpArrow   -> let i' = i - 1
                                   x  = lesserChar ! i'
                                   e  = (f x gap,    x, zero)
                               in  (# i', j , e #)
                  DiagArrow -> let i' = i - 1
                                   j' = j - 1
                                   x  = lesserChar ! i'
                                   y  = longerChar ! j'
                                   e  = (f x y, x, y)
                               in  (# i', j', e #)
        in  previousSequence `snoc` localContext
