-----------------------------------------------------------------------------
-- |
-- Module      :  Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.NeedlemanWunsch
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

module Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Wide.Swapping
  ( wideSwappingDO
  , directOptimization
  , traceback
  ) where

import           Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.DynamicCharacter2
import           Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Internal (Direction(..))
import           Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Wide.Internal
import           Control.Monad.ST
import           Data.Bits
import           Data.DList                  (snoc)
import           Data.Foldable
import           Data.Matrix.Unboxed         (Matrix, unsafeFreeze, unsafeIndex)
import qualified Data.Matrix.Unboxed.Mutable as M
import           Data.MonoTraversable
import           Data.Vector.Unboxed         (Vector,(!))
import qualified Data.Vector.Unboxed         as  UV
import qualified Data.Vector.Unboxed.Mutable as MUV
import           Data.Word


-- |
-- Performs a naive direct optimization.
-- Takes in two characters to run DO on and an overlap function
-- Returns an assignment character, the cost of that assignment, the assignment
-- character with gaps included, the aligned version of the first input character,
-- and the aligned version of the second input character. The process for this
-- algorithm is to generate a traversal matrix, then perform a traceback.
{-# SCC        wideSwappingDO #-}
{-# INLINEABLE wideSwappingDO #-}
wideSwappingDO
  :: Word8 -- ^ Alphbet size / symbol count
  -> (Word64 -> Word64 -> (Word64, Word))
  -> WideDynamicCharacter
  -> WideDynamicCharacter
  -> (Word, WideDynamicCharacter)
wideSwappingDO = directOptimization buildDirectionMatrix


buildDirectionMatrix
  :: Word8 -- ^ Alphbet size / symbol count
  -> (Word64 -> Word64 -> (Word64, Word))
  -> Vector Word64
  -> Vector Word64
  -> (Word, Matrix Direction)
buildDirectionMatrix symbolCount overlapλ topChar leftChar = fullMatrix
  where
    cost x y   = snd $ overlapλ x y
    gap        = bit . fromEnum $ symbolCount - 1
    rows       = olength leftChar + 1
    cols       = olength topChar  + 1

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
  :: (Word8 -> (Word64 -> Word64 -> (Word64, Word)) -> Vector Word64 -> Vector Word64 -> (Word, Matrix Direction))
  -> Word8
  -> (Word64 -> Word64 -> (Word64, Word))
  -> WideDynamicCharacter -- ^ Longer dynamic character
  -> WideDynamicCharacter -- ^ Shorter dynamic character
  -> (Word, WideDynamicCharacter)
directOptimization matrixFunction symbolCount overlapλ lhs rhs =
    let gap = bit $ fromEnum symbolCount :: Word64
        ~(swapped, gapsLesser, gapsLonger, (shorterChar,_,_), (longerChar,_,_)) = measureAndUngapCharacters symbolCount lhs rhs
        ~(alignmentCost, ungappedAlignment) =
            if      olength shorterChar == 0
            then if olength  longerChar == 0
                 -- Niether character was Missing, but both are empty when gaps are removed
                 then (0, (mempty,mempty,mempty))
                 -- Niether character was Missing, but one of them is empty when gaps are removed
                 else let vec = UV.generate (UV.length longerChar) $ \i -> fst (overlapλ (longerChar ! i) gap)
                      in  (0, (vec, UV.replicate (UV.length longerChar) 0, longerChar))
                 -- Both have some non-gap elements, perform string alignment
            else let (cost, traversalMatrix) =  matrixFunction symbolCount overlapλ longerChar shorterChar
                 in  (cost, traceback gap overlapλ traversalMatrix longerChar shorterChar)
        transformation    = if swapped then \(m,l,r) -> (m,r,l) else id
        regappedAlignment = insertGaps symbolCount gapsLesser gapsLonger ungappedAlignment
        alignmentContext  = transformation regappedAlignment
    in  (alignmentCost, alignmentContext)


{-# SCC traceback #-}
traceback
  :: Word64
  -> (Word64 -> Word64 -> (Word64, Word))
  -> Matrix Direction
  -> Vector Word64        -- ^ Longer dynamic character
  -> Vector Word64        -- ^ Shorter dynamic character
  -> WideDynamicCharacter -- ^ Resulting dynamic character
traceback gap overlapλ alignMatrix longerChar lesserChar = alignmentContext
  where
    f x y = fst $ overlapλ x y

    alignmentContext = dlistToDynamic $ go startPoint
    dlistToDynamic = UV.unzip3 . UV.fromList . toList

    longerLen = olength longerChar
    lesserLen = olength lesserChar
    rows      = lesserLen + 1
    cols      = longerLen + 1

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
                                   e  = (f gap y, 0, y)
                               in  (# i , j', e #)
                  UpArrow   -> let i' = i - 1
                                   x  = lesserChar ! i'
                                   e  = (f x gap, x, 0)
                               in  (# i', j , e #)
                  DiagArrow -> let i' = i - 1
                                   j' = j - 1
                                   x  = lesserChar ! i'
                                   y  = longerChar ! j'
                                   e  = (f x y, x, y)
                               in  (# i', j', e #)
        in  previousSequence `snoc` localContext
