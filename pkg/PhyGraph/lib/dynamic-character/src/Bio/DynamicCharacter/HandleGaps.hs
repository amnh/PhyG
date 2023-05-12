-----------------------------------------------------------------------------
-- |
-- Module      :  Bio.DynamicCharacter.HandleGaps
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

{-# Language ApplicativeDo #-}
{-# Language DerivingStrategies #-}
{-# Language GeneralizedNewtypeDeriving #-}
{-# Language ImportQualifiedPost #-}
{-# Language Strict #-}

module Bio.DynamicCharacter.HandleGaps
  ( GapSet()
  , nullGapSet
  , deleteGaps
  , insertGaps
  ) where

import Bio.DynamicCharacter
import Control.Monad (unless, when)
import Control.Monad.Loops (whileM_)
import Control.Monad.ST
import Data.Bits
import Data.Foldable
import Data.IntMap.Strict (IntMap, fromDistinctAscList, lookupMax, toAscList)
import Data.STRef
import Data.Semigroup
import Data.Vector.Generic (Vector, (!))
import Data.Vector.Generic qualified as GV
import Data.Vector.Generic.Mutable qualified as MGV
import Data.Vector.Unboxed.Mutable qualified as MUV


-- |
-- A set of gaps extracted froma sequence to be re-inserted later.
--
-- /NOTE:/ The 'GapSet' type is closed to ensure invaiants hold.
newtype GapSet = GapSet (IntMap Word)
    deriving newtype (Eq, Show)


-- |
-- The 'GapSet' which contains no gaps.
--
-- /NOTE:/ Useful for making comparisons.
nullGapSet :: GapSet
nullGapSet = GapSet mempty


-- |
-- Strips the gap elements from the supplied character.
--
-- Remembers the locations of the gap characters that were deleted
--
-- If the character contains /only/ gaps, a missing character is returned.
{-# INLINEABLE deleteGaps #-}
{-# SPECIALISE deleteGaps :: SlimDynamicCharacter -> (GapSet, SlimDynamicCharacter) #-}
{-# SPECIALISE deleteGaps :: WideDynamicCharacter -> (GapSet, WideDynamicCharacter) #-}
{-# SPECIALISE deleteGaps :: HugeDynamicCharacter -> (GapSet, HugeDynamicCharacter) #-}
deleteGaps
  :: ( FiniteBits e
     , Vector v e
     )
  => OpenDynamicCharacter v e -- ^ Dynamic character
  -> (GapSet, OpenDynamicCharacter v e)
deleteGaps c@(x,y,z)
  | GV.null x   = (noGaps,         c)
  | null gaps   = (noGaps,         c)
  | newLen == 0 = (gapSet,   missing)
  | otherwise   = (gapSet, newVector)
  where
    noGaps    = GapSet mempty
    missing   = (GV.empty, GV.empty, GV.empty)
    newVector = runST $ do
        j <- newSTRef 0
        let isGapAtJ = isGap c <$> readSTRef j

        let g v = do
              whileM_ isGapAtJ $ modifySTRef j succ
              j' <- readSTRef j
              modifySTRef j succ
              pure $ v ! j'

        let genVec v = do
              vec <- GV.generateM newLen (const (g v))
              writeSTRef j 0
              pure vec

        (,,) <$> genVec x <*> genVec y <*> genVec z

    gapCount = fromEnum . getSum $ foldMap Sum gaps
    charLen  = GV.length x
    newLen   = charLen - gapCount

    gapSet = GapSet gaps
    gaps = fromDistinctAscList $ reverse refs
    refs = runST $ do
       nonGaps <- newSTRef 0
       prevGap <- newSTRef False
       gapLen  <- newSTRef 0
       gapRefs <- newSTRef []

       let handleGapBefore op = do
               gapBefore <- readSTRef prevGap
               when gapBefore $ do
                 j <- readSTRef nonGaps
                 g <- readSTRef gapLen
                 modifySTRef gapRefs ( (j,g): )
                 op

       for_ [0 .. charLen - 1] $ \i ->
          if   isGap c i
          then modifySTRef gapLen succ *> writeSTRef prevGap True
          else do handleGapBefore $ do
                    writeSTRef  gapLen 0
                    writeSTRef prevGap False
                  modifySTRef nonGaps succ

       handleGapBefore $ pure ()
       readSTRef gapRefs


-- |
-- Adds gaps elements to the supplied character.
--
-- /NOTE:/ It is important to have the 'gap' state passed in as a parameter!
-- There is the possibility that the alignment context is empty, but one or more
-- gaps need to be added. We cannot construct a gap state from an empty input so
-- the gap state must be externally supplied.
{-# INLINEABLE insertGaps #-}
{-# SPECIALISE insertGaps :: SlimState -> GapSet -> GapSet -> SlimDynamicCharacter -> SlimDynamicCharacter #-}
{-# SPECIALISE insertGaps :: WideState -> GapSet -> GapSet -> WideDynamicCharacter -> WideDynamicCharacter #-}
{-# SPECIALISE insertGaps :: HugeState -> GapSet -> GapSet -> HugeDynamicCharacter -> HugeDynamicCharacter #-}
insertGaps
  :: ( FiniteBits e
     , Vector v e
     )
  => e      -- ^ Gap State
  -> GapSet -- ^ Removed gap references from "left"  (shorter) dynamic character
  -> GapSet -- ^ Removed gap references from "right" (larger)  dynamic character
  -> OpenDynamicCharacter v e -- ^ Alignment context to have gap references inserted
  -> OpenDynamicCharacter v e -- ^ Fully gapped alignment context
insertGaps gap (GapSet lGaps) (GapSet rGaps) input@(x,y,z)
  | null lGaps && null rGaps = input -- No work needed
  | otherwise =
  -- Use a Let/In binding here since the module is STRICT!
  -- Prevents indexing errors on missing character input.
  let nil       = gap `xor` gap
      totalGaps = fromEnum . getSum . foldMap Sum
      gapVecLen = maybe 0 (succ . fst) . lookupMax
      lGapCount = totalGaps lGaps
      rGapCount = totalGaps rGaps
      medLength = GV.length x
      newLength = lGapCount + rGapCount + medLength

      ins = (gap, gap, nil)
      del = (nil, gap, gap)

  in  runST $ do
--        xVec <- MGV.unsafeNew newLength
--        yVec <- MGV.unsafeNew newLength
--        zVec <- MGV.unsafeNew newLength
        xVec <- MGV.replicate newLength zeroBits
        yVec <- MGV.replicate newLength zeroBits
        zVec <- MGV.replicate newLength zeroBits
        lVec <- MUV.replicate (gapVecLen lGaps) 0
        rVec <- MUV.replicate (gapVecLen rGaps) 0
        lGap <- newSTRef 0
        mPtr <- newSTRef 0
        rGap <- newSTRef 0
        -- Write out to the mutable vectors
        for_ (toAscList lGaps) $ uncurry (MUV.unsafeWrite lVec)
        for_ (toAscList rGaps) $ uncurry (MUV.unsafeWrite rVec)

        let inc r = modifySTRef r succ

        let align i = do
              m <- readSTRef mPtr
              MGV.unsafeWrite xVec i $ x ! m
              MGV.unsafeWrite yVec i $ y ! m
              MGV.unsafeWrite zVec i $ z ! m
              inc mPtr
              -- Deletions leave a void in the left (shorter) character.
              --
              -- This means we are "consuming" an element from the right (longer) character
              -- Hence we increment the right gap reference, as we "progress" through the
              -- right character sequence's elements in the alignment context.
              when (isAlign input m || isDelete input m) $
                inc rGap
              -- Similar logic as above, however,
              -- Insertions leave a void in the right (longer) character.
              when (isAlign input m || isInsert input m) $
                inc lGap

        let insertGapWith i (xe,ye,ze) gapRef gapVec = do
              let len = MGV.length gapVec
              rg <- readSTRef gapRef
              v  <- if rg >= len then pure 0 else MGV.unsafeRead gapVec rg
              if   v == 0
              then pure False
              else do MGV.unsafeWrite xVec i xe
                      MGV.unsafeWrite yVec i ye
                      MGV.unsafeWrite zVec i ze
                      MGV.unsafeWrite gapVec rg $ v - 1
                      pure True

        for_ [0 .. newLength - 1] $ \i -> do
             written <- insertGapWith i ins lGap lVec
             unless written $ do
               written' <- insertGapWith i del rGap rVec
               unless written' $ align i

        x' <- GV.unsafeFreeze xVec
        y' <- GV.unsafeFreeze yVec
        z' <- GV.unsafeFreeze zVec
        pure (x', y', z')
