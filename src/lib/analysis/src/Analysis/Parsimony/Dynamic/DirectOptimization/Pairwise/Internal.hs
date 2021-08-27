-----------------------------------------------------------------------------
-- |
-- Module      :  Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Internal
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
{-# LANGUAGE DerivingStrategies    #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE Strict           #-}
{-# LANGUAGE TypeFamilies     #-}
{-# LANGUAGE UnboxedTuples    #-}


module Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Internal
  ( Cost
  , Direction(..)
  , OverlapFunction
    -- * General pre-processing
  , deleteGaps
  , insertGaps
  , measureCharacters
  , measureAndUngapCharacters
    -- * Alignment
  , directOptimization
  , traceback
  ) where

import           Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.DynamicCharacter2
import           Control.Monad                                                            (unless, when)
import           Control.Monad.Loops                                                      (whileM_)
import           Control.Monad.ST
import           Data.Bits
import           Data.DList                                                               (snoc)
import           Data.Foldable
import           Data.IntMap                                                              (IntMap)
import qualified Data.IntMap                                                              as IM
import           Data.Matrix.Unboxed                                                      (Matrix, unsafeIndex)
import           Data.Ord
import           Data.STRef
import           Data.Semigroup
import qualified Data.Vector                                                              as V
import           Data.Vector.Generic                                                      (Vector, (!))
import qualified Data.Vector.Generic                                                      as GV
import qualified Data.Vector.Generic.Mutable                                              as MGV
import qualified Data.Vector.Primitive                                                    as PV
import qualified Data.Vector.Unboxed                                                      as UV
import qualified Data.Vector.Unboxed.Mutable                                              as MUV
import           Data.Word                                                               (Word8)
import           Numeric.Extended.Natural


-- |
-- Which direction to align the character at a given matrix point.
--
-- It should be noted that the ordering of the three arrow types are important,
-- as it guarantees that the derived 'Ord' instance will have the following
-- property:
--
-- DiagArrow < LeftArrow < UpArrow
--
-- This means:
--
--   - DiagArrow has highest precedence when one or more costs are equal
--
--   - LeftArrow has second highest precedence when one or more costs are equal
--
--   -   UpArrow has lowest precedence when one or more costs are equal
--
-- Using this 'Ord' instance, we can resolve ambiguous transformations in a
-- deterministic way. Without loss of generality in determining the ordering,
-- we choose the same biasing as the C code called from the FFI for consistency.
data Direction = DiagArrow | LeftArrow | UpArrow
    deriving stock (Eq, Ord)


-- |
-- This internal type used for computing the alignment cost. This type has an
-- "infinity" value that is conveniently used for the barrier costs. The cost is
-- strictly non-negative, and possibly infinite.
type Cost = ExtendedNatural


-- |
-- A generalized function representation: the "overlap" between dynamic character
-- elements, supplying the corresponding median and cost to align the two
-- characters.
type OverlapFunction e = e -> e -> (e, Word)


newtype instance UV.MVector s Direction = MV_Direction (PV.MVector s Word8)


newtype instance UV.Vector   Direction  = V_Direction  (PV.Vector    Word8)


instance UV.Unbox Direction


instance MGV.MVector UV.MVector Direction where

    {-# INLINE basicLength #-}
    basicLength (MV_Direction v) = MGV.basicLength v

    {-# INLINE basicUnsafeSlice #-}
    basicUnsafeSlice i n (MV_Direction v) = MV_Direction $ MGV.basicUnsafeSlice i n v

    {-# INLINE basicOverlaps #-}
    basicOverlaps (MV_Direction v1) (MV_Direction v2) = MGV.basicOverlaps v1 v2

    {-# INLINE basicUnsafeNew #-}
    basicUnsafeNew n = MV_Direction <$> MGV.basicUnsafeNew n

    {-# INLINE basicInitialize #-}
    basicInitialize (MV_Direction v) = MGV.basicInitialize v

    {-# INLINE basicUnsafeReplicate #-}
    basicUnsafeReplicate n x = MV_Direction <$> MGV.basicUnsafeReplicate n (fromDirection x)

    {-# INLINE basicUnsafeRead #-}
    basicUnsafeRead (MV_Direction v) i = toDirection <$> MGV.basicUnsafeRead v i

    {-# INLINE basicUnsafeWrite #-}
    basicUnsafeWrite (MV_Direction v) i x = MGV.basicUnsafeWrite v i (fromDirection x)

    {-# INLINE basicClear #-}
    basicClear (MV_Direction v) = MGV.basicClear v

    {-# INLINE basicSet #-}
    basicSet (MV_Direction v) x = MGV.basicSet v (fromDirection x)

    {-# INLINE basicUnsafeCopy #-}
    basicUnsafeCopy (MV_Direction v1) (MV_Direction v2) = MGV.basicUnsafeCopy v1 v2

    basicUnsafeMove (MV_Direction v1) (MV_Direction v2) = MGV.basicUnsafeMove v1 v2

    {-# INLINE basicUnsafeGrow #-}
    basicUnsafeGrow (MV_Direction v) n = MV_Direction <$> MGV.basicUnsafeGrow v n


instance GV.Vector UV.Vector Direction where

    {-# INLINE basicUnsafeFreeze #-}
    basicUnsafeFreeze (MV_Direction v) = V_Direction <$> GV.basicUnsafeFreeze v

    {-# INLINE basicUnsafeThaw #-}
    basicUnsafeThaw (V_Direction v) = MV_Direction <$> GV.basicUnsafeThaw v

    {-# INLINE basicLength #-}
    basicLength (V_Direction v) = GV.basicLength v

    {-# INLINE basicUnsafeSlice #-}
    basicUnsafeSlice i n (V_Direction v) = V_Direction $ GV.basicUnsafeSlice i n v

    {-# INLINE basicUnsafeIndexM #-}
    basicUnsafeIndexM (V_Direction v) i = toDirection <$> GV.basicUnsafeIndexM v i

    basicUnsafeCopy (MV_Direction mv) (V_Direction v) = GV.basicUnsafeCopy mv v

    {-# INLINE elemseq #-}
    elemseq _ = seq


instance Show Direction where

    show DiagArrow = "↖"
    show LeftArrow = "←"
    show UpArrow   = "↑"


-- |
-- Strips the gap elements from the supplied character.
--
-- Remembers the locations of the gap characters that were deleted
--
-- If the character contains /only/ gaps, a missing character is returned.
{-# INLINEABLE deleteGaps #-}
{-# SPECIALISE deleteGaps :: SlimState -> SlimDynamicCharacter -> (IntMap Word, SlimDynamicCharacter) #-}
{-# SPECIALISE deleteGaps :: WideState -> WideDynamicCharacter -> (IntMap Word, WideDynamicCharacter) #-}
{-# SPECIALISE deleteGaps :: HugeState -> HugeDynamicCharacter -> (IntMap Word, HugeDynamicCharacter) #-}
deleteGaps
  :: ( Ord a
     , Vector v a
     )
  => a               -- ^ Gap state
  -> (v a, v a, v a) -- ^ Dynamic character
  -> (IntMap Word, (v a, v a, v a))
deleteGaps gap c@(x,y,z)
  | GV.null x    = (mempty,       c)
  | null gaps   = (gaps,         c)
  | newLen == 0 = (gaps,   missing)
  | otherwise   = (gaps, newVector)
  where
    missing   = (GV.empty, GV.empty, GV.empty)
    newVector = runST $ do
        j <- newSTRef 0
        let isGapAtJ = do
              j' <- readSTRef j
              pure $ j' < charLen && x ! j' == gap
        let g v = do
              whileM_ isGapAtJ (modifySTRef j succ)
              j' <- readSTRef j
              modifySTRef j succ
              pure $ v ! j'
        x' <- GV.generateM newLen $ const (g x)
        writeSTRef j 0
        y' <- GV.generateM newLen $ const (g y)
        writeSTRef j 0
        z' <- GV.generateM newLen $ const (g z)
        pure (x', y', z')

    gapCount = fromEnum . getSum $ foldMap Sum gaps
    charLen  = GV.length x
    newLen   = charLen - gapCount

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
          if x ! i == gap
          then modifySTRef gapLen succ *> writeSTRef prevGap True
          else do handleGapBefore $ do
                    writeSTRef  gapLen 0
                    writeSTRef prevGap False
                  modifySTRef nonGaps succ

       handleGapBefore $ pure ()
       readSTRef gapRefs


-- |
-- Adds gaps elements to the supplied character.
{-# INLINEABLE insertGaps #-}
{-# SPECIALISE insertGaps :: SlimState -> IntMap Word -> IntMap Word -> SlimDynamicCharacter -> SlimDynamicCharacter #-}
{-# SPECIALISE insertGaps :: WideState -> IntMap Word -> IntMap Word -> WideDynamicCharacter -> WideDynamicCharacter #-}
{-# SPECIALISE insertGaps :: HugeState -> IntMap Word -> IntMap Word -> HugeDynamicCharacter -> HugeDynamicCharacter #-}
insertGaps
  :: ( FiniteBits a
     , Vector v a
     )
  => a               -- ^ Gap state
  -> IntMap Word     -- ^ Removed gap references from "left"  (shorter) dynamic character
  -> IntMap Word     -- ^ Removed gap references from "right" (larger)  dynamic character
  -> (v a, v a, v a) -- ^ Alignment context to have gap references inserted
  -> (v a, v a, v a) -- ^ Fully gapped alignment context
insertGaps gap lGaps rGaps meds@(x,y,z)
  | null lGaps && null rGaps = meds -- No work needed
  | otherwise                = newVector
  where
    zero      = gap `xor` gap
    totalGaps = fromEnum . getSum . foldMap Sum
    gapVecLen = maybe 0 (succ . fst) . IM.lookupMax
    lGapCount = totalGaps lGaps
    rGapCount = totalGaps rGaps
    newLength = lGapCount + rGapCount + GV.length x

    ins = (gap,  gap, zero)
    del = (gap, zero,  gap)

    newVector = runST $ do
      xVec <- MGV.unsafeNew newLength
      yVec <- MGV.unsafeNew newLength
      zVec <- MGV.unsafeNew newLength
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
            MGV.unsafeWrite xVec i $ x ! m
            MGV.unsafeWrite yVec i $ y ! m
            MGV.unsafeWrite zVec i $ z ! m
            modifySTRef mPtr succ
            when (isAlign meds m || isDelete meds m) $ do
              modifySTRef rGap succ
            when (isAlign meds m || isInsert meds m) $ do
              modifySTRef lGap succ

      let insertGapWith i (xe,ye,ze) gapRef gapVec = do
            rg <- readSTRef gapRef
            v  <- if rg >= MGV.length gapVec then pure 0 else MGV.unsafeRead gapVec rg
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


{-# SCC directOptimization #-}
{-# INLINEABLE directOptimization #-}
{-# SPECIALISE directOptimization :: (WideState -> (WideState -> WideState -> (WideState, Word)) -> UV.Vector WideState -> UV.Vector WideState -> (Word, Matrix Direction)) -> WideState -> (WideState -> WideState -> (WideState, Word)) -> WideDynamicCharacter -> WideDynamicCharacter -> (Word, WideDynamicCharacter) #-}
{-# SPECIALISE directOptimization :: (HugeState -> (HugeState -> HugeState -> (HugeState, Word)) ->  V.Vector HugeState ->  V.Vector HugeState -> (Word, Matrix Direction)) -> HugeState -> (HugeState -> HugeState -> (HugeState, Word)) -> HugeDynamicCharacter -> HugeDynamicCharacter -> (Word, HugeDynamicCharacter) #-}
directOptimization
  :: ( FiniteBits a
     , Ord a
     , Ord (v a)
     , Vector v a
     , Vector v (a, a, a)
     )
  => (a -> (a -> a -> (a, Word)) -> v a -> v a -> (Word, Matrix Direction))
  -> a                     -- ^ Gap state
  -> (a -> a -> (a, Word)) -- ^ Metric for coputing state distance and median state
  -> (v a, v a, v a)       -- ^ Longer  dynamic character
  -> (v a, v a, v a)       -- ^ Shorter dynamic character
  -> (Word, (v a, v a, v a))
directOptimization matrixFunction gap overlapλ lhs rhs =
  case (isMissing lhs, isMissing rhs) of
    (True , True ) -> (0, lhs)
    (True , False) -> (0, rhs)
    (False, True ) -> (0, lhs)
    (False, False) ->
      let ~(swapped, gapsLesser, gapsLonger, (shorterChar,_,_), (longerChar,_,_)) = measureAndUngapCharacters gap lhs rhs
          ~(alignmentCost, ungappedAlignment) =
              if      GV.length shorterChar == 0
              then if GV.length  longerChar == 0
                   -- Niether character was Missing, but both are empty when gaps are removed
                   then (0, (GV.empty, GV.empty, GV.empty))
                   -- Niether character was Missing, but one of them is empty when gaps are removed
                   else let len = GV.length longerChar
                            vec = GV.generate len $ \i -> fst (overlapλ (longerChar ! i) gap)
                        in  (0, (vec, GV.replicate len (gap `xor` gap), longerChar))
                   -- Both have some non-gap elements, perform string alignment
              else let (cost, traversalMatrix) = matrixFunction gap overlapλ longerChar shorterChar
                   in  (cost, traceback gap overlapλ traversalMatrix longerChar shorterChar)
          transformation    = if swapped then \(m,l,r) -> (m,r,l) else id
          regappedAlignment = insertGaps gap gapsLesser gapsLonger ungappedAlignment
          alignmentContext  = transformation regappedAlignment
      in  (alignmentCost, alignmentContext)


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
{-# SPECIALISE measureCharacters :: SlimDynamicCharacter -> SlimDynamicCharacter -> (Ordering, SlimDynamicCharacter, SlimDynamicCharacter) #-}
{-# SPECIALISE measureCharacters :: WideDynamicCharacter -> WideDynamicCharacter -> (Ordering, WideDynamicCharacter, WideDynamicCharacter) #-}
{-# SPECIALISE measureCharacters :: HugeDynamicCharacter -> HugeDynamicCharacter -> (Ordering, HugeDynamicCharacter, HugeDynamicCharacter) #-}
measureCharacters
  :: ( Ord (v a)
     , Vector v a
     )
  => (v a, v a, v a)
  -> (v a, v a, v a)
  -> (Ordering, (v a, v a, v a), (v a, v a, v a))
measureCharacters lhs@(lhsMedians,_,_) rhs@(rhsMedians,_,_)
  | lhsOrdering == GT = (lhsOrdering, rhs, lhs)
  | otherwise         = (lhsOrdering, lhs, rhs)
  where
    lhsOrdering =
        -- First, compare inputs by length.
        case comparing GV.length lhsMedians rhsMedians of
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
-- /O(n)/
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
{-# SPECIALISE measureAndUngapCharacters :: SlimState -> SlimDynamicCharacter -> SlimDynamicCharacter -> (Bool, IntMap Word, IntMap Word, SlimDynamicCharacter, SlimDynamicCharacter) #-}
{-# SPECIALISE measureAndUngapCharacters :: WideState -> WideDynamicCharacter -> WideDynamicCharacter -> (Bool, IntMap Word, IntMap Word, WideDynamicCharacter, WideDynamicCharacter) #-}
{-# SPECIALISE measureAndUngapCharacters :: HugeState -> HugeDynamicCharacter -> HugeDynamicCharacter -> (Bool, IntMap Word, IntMap Word, HugeDynamicCharacter, HugeDynamicCharacter) #-}
measureAndUngapCharacters
  :: ( Ord a
     , Ord (v a)
     , Vector v a
     )
  => a                -- ^ Gap state
  -> (v a, v a, v a)  -- ^ First  dynamic character
  -> (v a, v a, v a)  -- ^ Second dynamic character
  -> (Bool, IntMap Word, IntMap Word, (v a, v a, v a), (v a, v a, v a))
measureAndUngapCharacters gap char1 char2
  | swapInputs = (True , gapsChar2, gapsChar1, ungappedChar2, ungappedChar1)
  | otherwise  = (False, gapsChar1, gapsChar2, ungappedChar1, ungappedChar2)
  where
    swapInputs = measure == GT
    (gapsChar1, ungappedChar1) = deleteGaps gap char1
    (gapsChar2, ungappedChar2) = deleteGaps gap char2
    (measure, _, _) =
        case measureCharacters ungappedChar1 ungappedChar2 of
          (EQ,_,_) -> measureCharacters char1 char2
          x        -> x


{-# SCC traceback #-}
{-# INLINEABLE traceback #-}
{-# SPECIALISE traceback :: WideState -> (WideState -> WideState -> (WideState, Word)) -> Matrix Direction -> UV.Vector WideState -> UV.Vector WideState -> WideDynamicCharacter #-}
{-# SPECIALISE traceback :: HugeState -> (HugeState -> HugeState -> (HugeState, Word)) -> Matrix Direction ->  V.Vector HugeState ->  V.Vector HugeState -> HugeDynamicCharacter #-}
traceback
  :: ( Bits a
     , Vector v a
     , Vector v (a, a, a)
     )
  => a
  -> (a -> a -> (a, Word))
  -> Matrix Direction
  -> v a -- ^ Longer  dynamic character states
  -> v a -- ^ Shorter dynamic character states
  -> (v a, v a, v a) -- ^ Resulting dynamic character alignment context
traceback gap overlapλ alignMatrix longerChar lesserChar = alignmentContext
  where
    f x y = fst $ overlapλ x y

    alignmentContext = dlistToDynamic $ go startPoint
    dlistToDynamic   = GV.unzip3 . GV.fromList . toList

    longerLen  = GV.length longerChar
    lesserLen  = GV.length lesserChar
    rows       = lesserLen + 1
    cols       = longerLen + 1
    zero       = gap `xor` gap
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


{--
-- |
-- Serializes an alignment matrix to a 'String'. Uses input characters for row
-- and column labelings.
--
-- Useful for debugging purposes.
renderCostMatrix
  :: ( DOCharConstraint s
     , Foldable f
     , Functor f
     , Indexable f
     , Key f ~ (Int, Int)
     , Show a
     , Show b
     )
  => s
  -> s
  -> f (a, b) -- ^ The Needleman-Wunsch alignment matrix
  -> String
renderCostMatrix lhs rhs mtx = unlines
    [ dimensionPrefix
    , headerRow
    , barRow
    , renderedRows
    ]
  where
    (_,longer,lesser) = measureCharacters lhs rhs
    longerTokens      = toShownIntegers longer
    lesserTokens      = toShownIntegers lesser
--    toShownIntegers   = fmap (show . showBitsValue) . otoList
    toShownIntegers   = fmap (const "#") . otoList
    matrixTokens      = showCell <$> mtx
    showCell (c,d)    = show c <> show d
    maxPrefixWidth    = maxLengthOf lesserTokens
    maxColumnWidth    = max (maxLengthOf longerTokens) . maxLengthOf $ toList matrixTokens
    maxLengthOf       = maximum . fmap length

{-
    showBitsValue :: FiniteBits b => b -> Word
    showBitsValue b = go (finiteBitSize b) 0
      where
        go 0 v = v
        go i v = let i' = i-1
                     v' | b `testBit` i' = v + bit i'
                        | otherwise      = v
                 in  go i' v'
-}

    colCount = olength longer + 1
    rowCount = olength lesser + 1

    dimensionPrefix  = " " <> unwords
        [ "Dimensions:"
        , show rowCount
        , "X"
        , show colCount
        ]

    headerRow = mconcat
        [ " "
        , pad maxPrefixWidth "\\"
        , "| "
        , pad maxColumnWidth "*"
        , concatMap (pad maxColumnWidth) longerTokens
        ]

    barRow    = mconcat
        [ " "
        , bar maxPrefixWidth
        , "+"
        , concatMap (const (bar maxColumnWidth)) $ undefined : longerTokens
        ]
      where
        bar n = replicate (n+1) '-'

    renderedRows = unlines . zipWith renderRow ("*":lesserTokens) $ getRows matrixTokens
      where
        renderRow e vs = " " <> pad maxPrefixWidth e <> "| " <> concatMap (pad maxColumnWidth) vs

        getRows m = (`getRow'` m) <$> [0 .. rowCount - 1]
        getRow' i m = g <$> [0 .. colCount - 1]
          where
            g j = fromMaybe "" $ (i,j) `lookup` m


    pad :: Int -> String -> String
    pad n e = replicate (n - len) ' ' <> e <> " "
      where
        len = length e
--}


{-# INLINE fromDirection #-}
fromDirection :: Direction -> Word8
fromDirection DiagArrow = 0
fromDirection LeftArrow = 1
fromDirection UpArrow   = 2


{-# INLINE toDirection #-}
toDirection :: Word8 -> Direction
toDirection 0 = DiagArrow
toDirection 1 = LeftArrow
toDirection _ = UpArrow
