-----------------------------------------------------------------------------
-- |
-- Module      :  DirectOptimization.Pairwise.Internal
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
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE Strict           #-}

module DirectOptimization.Pairwise.Internal
  ( -- * Alignment types
    Direction(..)
  , Distance
  , TCM2Dλ
    -- * Alignment generic functions
  , directOptimization
  , directOptimizationFromDirectionMatrix
    -- * Reused helper functions
  , getCostλ
  ) where

import           Bio.DynamicCharacter
import           Bio.DynamicCharacter.HandleGaps
import           Bio.DynamicCharacter.Measure
import           Control.Applicative
import           Control.Monad                         (join)
import           Control.Monad.Loops                   (whileM_)
import           Data.Bits
import           Data.Matrix.Class                     (Matrix, dim, unsafeIndex)
import qualified Data.Matrix.Unboxed                   as UM
import           Data.STRef
import qualified Data.Vector                           as V
import           Data.Vector.Generic                   (Vector, (!), (!?))
import qualified Data.Vector.Generic                   as GV
import qualified Data.Vector.Storable                  as SV
import qualified Data.Vector.Unboxed                   as UV
import           DirectOptimization.Pairwise.Direction
import           Measure.Matrix
import           Measure.Unit


{-# SCC directOptimization #-}
{-# INLINEABLE directOptimization #-}
{-# SPECIALISE directOptimization :: (SV.Vector SlimState -> SV.Vector SlimState -> (Distance, SlimDynamicCharacter)) -> TCM2Dλ SlimState -> SlimDynamicCharacter -> SlimDynamicCharacter -> (Distance, SlimDynamicCharacter) #-}
{-# SPECIALISE directOptimization :: (UV.Vector WideState -> UV.Vector WideState -> (Distance, WideDynamicCharacter)) -> TCM2Dλ WideState -> WideDynamicCharacter -> WideDynamicCharacter -> (Distance, WideDynamicCharacter) #-}
{-# SPECIALISE directOptimization :: ( V.Vector HugeState ->  V.Vector HugeState -> (Distance, HugeDynamicCharacter)) -> TCM2Dλ HugeState -> HugeDynamicCharacter -> HugeDynamicCharacter -> (Distance, HugeDynamicCharacter) #-}
directOptimization
  :: ( FiniteBits e
     , Ord (v e)
     , Vector v e
     )
  => (v e -> v e -> (Distance, OpenDynamicCharacter v e)) -- ^ Alignment function
  -> TCM2Dλ e -- ^ Metric for computing state distance and median state
  -> OpenDynamicCharacter v e
  -> OpenDynamicCharacter v e
  -> (Distance, OpenDynamicCharacter v e)
directOptimization alignmentFunction overlapλ = handleMissing generateAlignmentResult
  where
    generateAlignmentResult lhs rhs =
      let -- Build a 'gap' state now that we know that we can access a non-empty sequence.
          gap = let tmp = extractMediansGapped rhs ! 0 in (tmp `xor` tmp) `setBit` 0
          -- Remove gaps from the inputs and measure the results to determine
          -- which ungapped character is longer and which is shorter.
          -- Always pass the shorter character into alignment functions first!
          ~(swapped, gapsLesser, gapsLonger, lesser, longer) = measureCharactersWithoutGaps lhs rhs
          lesserMeds = extractMediansGapped lesser
          longerMeds = extractMediansGapped longer
          ~(alignmentCost, ungappedAlignment) =
              if      GV.length lesserMeds == 0
              -- Neither character was Missing, but one or both are empty when gaps are removed
              then alignmentWithAllGaps overlapλ longerMeds
              -- Both have some non-gap elements, perform string alignment
              else alignmentFunction lesserMeds longerMeds
          regappedAlignment = insertGaps gap gapsLesser gapsLonger ungappedAlignment
          transformation    = if swapped then transposeCharacter else id
          alignmentContext  = transformation regappedAlignment
      in  (alignmentCost, alignmentContext)


{-# SCC directOptimizationFromDirectionMatrix #-}
{-# INLINEABLE directOptimizationFromDirectionMatrix #-}
{-# SPECIALISE directOptimizationFromDirectionMatrix :: (WideState -> TCM2Dλ WideState -> UV.Vector WideState -> UV.Vector WideState -> (Distance, UM.Matrix Direction)) -> TCM2Dλ WideState -> WideDynamicCharacter -> WideDynamicCharacter -> (Distance, WideDynamicCharacter) #-}
{-# SPECIALISE directOptimizationFromDirectionMatrix :: (HugeState -> TCM2Dλ HugeState ->  V.Vector HugeState ->  V.Vector HugeState -> (Distance, UM.Matrix Direction)) -> TCM2Dλ HugeState -> HugeDynamicCharacter -> HugeDynamicCharacter -> (Distance, HugeDynamicCharacter) #-}
directOptimizationFromDirectionMatrix
  :: ( FiniteBits e
     , Ord (v e)
     , Vector v e
     , Matrix m t Direction
     )
  => (e -> TCM2Dλ e -> v e -> v e -> (Distance, m t Direction)) -- ^ Alignment matrix generator function
  -> TCM2Dλ e -- ^ Metric for computing state distance and median state
  -> OpenDynamicCharacter v e
  -> OpenDynamicCharacter v e
  -> (Distance, OpenDynamicCharacter v e)
directOptimizationFromDirectionMatrix matrixGenerator overlapλ =
    handleMissing $ directOptimization alignmentFunction overlapλ
  where
    alignmentFunction lhs rhs =
        let gap = let tmp = rhs ! 0 in (tmp `xor` tmp) `setBit` 0
            (cost, traversalMatrix) = matrixGenerator gap overlapλ lhs rhs
        in  (cost, traceback gap overlapλ traversalMatrix lhs rhs)


{-# INLINEABLE getCostλ #-}
getCostλ :: TCM2Dλ e -> e -> e -> Word
getCostλ =
    let normalize :: (Distance, e) -> Word
        normalize = fromIntegral . fst
    in  (.).(.) $ normalize


{-# SCC traceback #-}
{-# INLINEABLE traceback #-}
{-# SPECIALISE traceback :: WideState -> TCM2Dλ WideState -> UM.Matrix Direction -> UV.Vector WideState -> UV.Vector WideState -> WideDynamicCharacter #-}
{-# SPECIALISE traceback :: HugeState -> TCM2Dλ HugeState -> UM.Matrix Direction ->  V.Vector HugeState ->  V.Vector HugeState -> HugeDynamicCharacter #-}
traceback
  :: ( Bits e
     , Matrix m t Direction
     , Vector v e
     )
  => e
  -> TCM2Dλ e
  -> m t Direction
  -> v e -- ^ Shorter dynamic character related to the "left column"
  -> v e -- ^ Longer  dynamic character related to the "top row"
  -> OpenDynamicCharacter v e -- ^ Resulting dynamic character alignment context
traceback gap overlapλ directionMatrix lesser longer = alignment
  where
    f x y = snd $ overlapλ x y
    getDirection = curry $ unsafeIndex directionMatrix
    -- The maximum size the alignment could be
    bufferLength = toEnum $ GV.length lesser + GV.length longer

    -- Construct the aligned dynamic character by using a buffer of it's maximum
    -- possible length. Computet the length will performing the traceback. Then
    -- after the traceback is performed, copy from the buffer to create a dynamic
    -- character of the correct size.
    --
    -- NOTE: The buffer creation, copying, and freeing are all handled by the
    --       'unsafeCharacterBuiltByBufferedST' function.
    alignment = unsafeCharacterBuiltByBufferedST bufferLength $ \a -> do
        let (m,n) = dim directionMatrix
        iR <- newSTRef $ m - 1
        jR <- newSTRef $ n - 1
        kR <- newSTRef $ fromEnum bufferLength

        -- Set up convenience methods for accessing character elements
        let getElement char ref = (char !) <$> (modifySTRef ref pred *> readSTRef ref)
        let getLesserElement = getElement lesser iR
        let getLongerElement = getElement longer jR

        -- Set up convenience methods for setting alignment elements on traceback
        let setElementAt = modifySTRef kR pred *> readSTRef kR
        let delete    re = setElementAt >>= \k -> setDelete a k    (f gap re) re
        let align  le re = setElementAt >>= \k -> setAlign  a k le (f le  re) re
        let insert le    = setElementAt >>= \k -> setInsert a k le (f le gap)

        -- Determine when to break the alignment loop
        let continueAlignment =
               let notAtOrigin i j = i /= 0 || j /= 0
               in  liftA2 notAtOrigin (readSTRef iR) $ readSTRef jR

        -- Perform traceback
        whileM_ continueAlignment $ do
            arrow <- liftA2 getDirection (readSTRef iR) $ readSTRef jR
            case arrow of
                   LeftArrow -> getLongerElement >>= delete
                   DiagArrow -> join $ liftA2 align getLesserElement getLongerElement
                   UpArrow   -> getLesserElement >>= insert

        -- Return the actual alignment length
        k <- readSTRef kR
        pure $ bufferLength - toEnum k


{-# SCC alignmentWithAllGaps #-}
{-# INLINEABLE alignmentWithAllGaps #-}
alignmentWithAllGaps
  :: ( Bits e
     , Vector v e
     )
  => TCM2Dλ e
  -> v e
  -> (Distance, OpenDynamicCharacter v e)
alignmentWithAllGaps overlapλ character =
    case character !? 0 of
      -- Neither character was Missing, but both are empty when gaps are removed
      Nothing -> (0, (GV.empty, GV.empty, GV.empty))
      -- Neither character was Missing, but one of them is empty when gaps are removed
      Just e  ->
        let len = GV.length character
            nil = e `xor` e
            gap = nil`setBit` 0
            zed = GV.replicate len nil
            med = GV.generate  len $ snd . overlapλ gap . (character !)
        in  (0, (zed, med, character))


{-# SCC handleMissing #-}
handleMissing
  :: Vector v e
  => (OpenDynamicCharacter v e -> OpenDynamicCharacter v e -> (Distance, OpenDynamicCharacter v e))
  -> OpenDynamicCharacter v e
  -> OpenDynamicCharacter v e
  -> (Distance, OpenDynamicCharacter v e)
handleMissing f lhs rhs =
  case (isMissing lhs, isMissing rhs) of
    (True , True ) -> (0, lhs)
    (True , False) -> (0, rhs)
    (False, True ) -> (0, lhs)
    (False, False) -> f lhs rhs
