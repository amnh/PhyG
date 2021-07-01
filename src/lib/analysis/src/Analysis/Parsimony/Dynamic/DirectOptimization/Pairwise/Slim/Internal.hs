-----------------------------------------------------------------------------
-- |
--
-----------------------------------------------------------------------------

{-# LANGUAGE Strict #-}

module Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Slim.Internal
  ( deleteGaps
  , insertGaps
  , measureCharacters
  , measureAndUngapCharacters
  ) where

import           Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.DynamicCharacter2
import           Control.Monad                                  (unless, when)
import           Control.Monad.Loops                            (whileM_)
import           Control.Monad.ST
import           Data.Bits
import           Data.Foldable
import           Data.IntMap          (IntMap)
import qualified Data.IntMap          as IM
import           Data.Ord
import           Data.STRef
import           Data.Semigroup
import qualified Data.Vector.Unboxed.Mutable                    as MUV
import qualified Data.Vector.Storable                           as V
import qualified Data.Vector.Storable.Mutable                   as MV
import           Data.Word


-- |
-- Strips the gap elements from the supplied character.
--
-- Remembers the locations of the gap characters that were deleted
--
-- If the character contains /only/ gaps, a missing character is returned.
{-# INLINEABLE deleteGaps #-}
deleteGaps
  :: Word8
  -> SlimDynamicCharacter
  -> (IntMap Word, SlimDynamicCharacter)
deleteGaps symbolCount c@(x,y,z)
  | V.null x    = (mempty,       c)
  | null gaps   = (gaps,         c)
  | newLen == 0 = (gaps,    mempty)
  | otherwise   = (gaps, newVector)
  where
    newVector = runST $ do
        j <- newSTRef 0
        let isGapAtJ = do
              j' <- readSTRef j
              pure $ j' < charLen && x V.! j' == gap
        let g v = do
              whileM_ isGapAtJ (modifySTRef j succ)
              j' <- readSTRef j
              modifySTRef j succ
              pure $ v V.! j'
        x' <- V.generateM newLen $ const (g x)
        y' <- V.generateM newLen $ const (g y)
        z' <- V.generateM newLen $ const (g z)
        pure (x', y', z')

    gapCount = fromEnum . getSum $ foldMap Sum gaps
    charLen  = V.length x
    newLen   = charLen - gapCount
    gap      = bit . fromEnum $ symbolCount - 1

    gaps = IM.fromDistinctAscList $ reverse refs
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
          if x V.! i == gap
          then modifySTRef gapLen succ *> writeSTRef prevGap True
          else do handleGapBefore $ do
                    writeSTRef  gapLen 0
                    writeSTRef prevGap False
                  modifySTRef nonGaps succ
  
       handleGapBefore $ pure ()
       readSTRef gapRefs


isAlign, isDelete, isInsert :: SlimDynamicCharacter -> Int -> Bool
isInsert (_,y,z) i = y V.! i /= 0 && z V.! i == 0
isDelete (_,y,z) i = y V.! i == 0 && z V.! i /= 0
isAlign  (_,y,z) i = y V.! i /= 0 && z V.! i /= 0


-- |
-- Adds gaps elements to the supplied character.
insertGaps
  :: Word8
  -> IntMap Word
  -> IntMap Word
  -> SlimDynamicCharacter
  -> SlimDynamicCharacter
insertGaps symbolCount lGaps rGaps meds@(x,y,z)
  | null lGaps && null rGaps = meds -- No work needed
  | otherwise                = newVector
  where
    gap      = bit . fromEnum $ symbolCount - 1
    totalGaps = fromEnum . getSum . foldMap Sum
    gapVecLen = maybe 0 (succ . fst) . IM.lookupMax
    lGapCount = totalGaps lGaps
    rGapCount = totalGaps rGaps
    newLength = lGapCount + rGapCount + V.length x

    ins = (gap, gap,   0)
    del = (gap,   0, gap)

    newVector = runST $ do
      xVec <- MV.unsafeNew newLength
      yVec <- MV.unsafeNew newLength
      zVec <- MV.unsafeNew newLength
      lVec <- MUV.replicate (gapVecLen lGaps) 0
      rVec <- MUV.replicate (gapVecLen rGaps) 0
      lGap <- newSTRef 0
      mPtr <- newSTRef 0
      rGap <- newSTRef 0
   -- Write out to the mutable vectors
      for_ (IM.toAscList lGaps) $ uncurry (MUV.unsafeWrite lVec)
      for_ (IM.toAscList rGaps) $ uncurry (MUV.unsafeWrite rVec)

      let align i = do
            m <- readSTRef mPtr
            MV.unsafeWrite xVec i $ x V.! m
            MV.unsafeWrite yVec i $ y V.! m
            MV.unsafeWrite zVec i $ z V.! m
            modifySTRef mPtr succ
            when (isAlign meds m || isDelete meds m) $ do
              modifySTRef lGap succ
            when (isAlign meds m || isInsert meds m) $ do
              modifySTRef rGap succ

      let insertGapWith i (xe,ye,ze) gapRef gapVec = do
            rg <- readSTRef gapRef
            v  <- if rg >= MUV.length gapVec then pure 0 else MUV.unsafeRead gapVec rg
            if   v == 0
            then pure False
            else do MV.unsafeWrite xVec i xe
                    MV.unsafeWrite yVec i ye
                    MV.unsafeWrite zVec i ze
                    MUV.unsafeWrite gapVec rg $ v - 1
                    pure True

      for_ [0 .. newLength - 1] $ \i -> do
           written <- insertGapWith i ins lGap lVec
           unless written $ do
             written' <- insertGapWith i del rGap rVec
             unless written' $ align i

      x' <- V.unsafeFreeze xVec
      y' <- V.unsafeFreeze yVec
      z' <- V.unsafeFreeze zVec
      pure (x', y', z')

-- |
-- /O(1)/ for input characters of differing lengths
--
-- /O(k)/ for input characters of equal length, where /k/ is the shared prefix of
-- both characters.
--
-- Returns the dynamic character that is shorter first, longer second, and notes
-- whether or not the inputs were swapped to place the characters in this ordering.
--
-- Handles equal length characters by considering the lexicographically larger
-- character as longer.
--
-- Handles equality of inputs by /not/ swapping.
{-# INLINEABLE measureCharacters #-}
measureCharacters
  :: SlimDynamicCharacter
  -> SlimDynamicCharacter
  -> (Ordering, SlimDynamicCharacter, SlimDynamicCharacter)
measureCharacters lhs@(lhsMedians,_,_) rhs@(rhsMedians,_,_)
  | lhsOrdering == GT = (lhsOrdering, rhs, lhs)
  | otherwise         = (lhsOrdering, lhs, rhs)
  where
    lhsOrdering =
        -- First, compare inputs by length.
        case comparing V.length lhsMedians rhsMedians of
          -- If the inputs are equal length,
          -- Then compare by the (arbitrary) lexicographical ordering of the median states.
          EQ -> case lhsMedians `compare` rhsMedians of
                  -- If the input median states have the same ordering,
                  -- Lastly, we compare by the lexicographic ordering of the "tagged triples."
                  --
                  -- If they are equal after this step,
                  -- Then the inputs are representationally equal.
                  -- Actually, honest to goodness 100% equal!
                  EQ -> lhs `compare` rhs
                  v  -> v
          v  -> v


-- |
-- /O(1)/ for input characters of differing lengths
--
-- /O(k)/ for input characters of equal length, where /k/ is the shared prefix of
-- both characters.
--
-- Considers the median values of the characters, ignores the left/right tagging.
--
-- First remove the gaps from the input characters.
--
-- If both "ungapped" inputs are empty, we measure the original "gapped" inputs to
-- determine if the inputs need to be swapped. This is required to ensure comutativity
-- of subsequent operations which use this method.
--
-- Returns the "ungapped" dynamic character that is "shorter" first, "longer" second,
-- the removed gap mappings (in the same order), and notes whether or not the inputs
-- were swapped to place the characters in this ordering.
--
-- Handles equal length characters by considering the lexicographically larger
-- character as longer.
--
-- Handles equality of inputs by /not/ swapping.
{-# INLINEABLE measureAndUngapCharacters #-}
measureAndUngapCharacters
  :: Word8
  -> SlimDynamicCharacter
  -> SlimDynamicCharacter
  -> (Bool, IntMap Word, IntMap Word, SlimDynamicCharacter, SlimDynamicCharacter)
measureAndUngapCharacters symbolCount char1 char2
  | swapInputs = (True , gapsChar2, gapsChar1, ungappedChar2, ungappedChar1)
  | otherwise  = (False, gapsChar1, gapsChar2, ungappedChar1, ungappedChar2)
  where
    swapInputs = measure == GT
    (gapsChar1, ungappedChar1) = deleteGaps symbolCount char1
    (gapsChar2, ungappedChar2) = deleteGaps symbolCount char2
    (measure, _, _) =
        case measureCharacters ungappedChar1 ungappedChar2 of
          (EQ,_,_) -> measureCharacters char1 char2
          x        -> x
