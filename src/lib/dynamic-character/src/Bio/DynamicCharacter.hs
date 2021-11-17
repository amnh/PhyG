{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs            #-}
{-# LANGUAGE RankNTypes       #-}
{-# LANGUAGE Strict           #-}
{-# LANGUAGE UnboxedTuples    #-}

module Bio.DynamicCharacter
  (  -- * Element Varieties of a Dynamic Character
    SlimState
  , WideState
  , HugeState
    -- * Generic Dynamic Character Constructions
  , OpenDynamicCharacter
  , TempOpenDynamicCharacter
    -- * Dynamic Character Immutable Varieties
  , SlimDynamicCharacter
  , WideDynamicCharacter
  , HugeDynamicCharacter
    -- * Dynamic Character Mutable Varieties
  , TempSlimDynamicCharacter
  , TempWideDynamicCharacter
  , TempHugeDynamicCharacter
    -- * Queries
  , isAlign
  , isDelete
  , isInsert
  , isGapped
  , isGap
  , isMissing
  , characterLength
    -- * Mutators
  , setAlign
  , setDelete
  , setInsert
  , setGapped
  , setFrom
  , transposeCharacter
    -- * Extractors
  , extractMedians
  , extractMediansLeft
  , extractMediansRight
  , extractMediansGapped
  , extractMediansLeftGapped
  , extractMediansRightGapped
    -- * Immutable Constructors
  , encodeDynamicCharacter
  , encodeState
    -- * Mutable Constructors
  , newTempCharacter
  , freezeTempCharacter
  , unsafeCharacterBuiltByST
  , unsafeCharacterBuiltByBufferedST
    -- * Rendering
  , decodeState
  , renderDynamicCharacter
  ) where

import           Control.Monad.Primitive
import           Control.Monad.ST
import           Data.Alphabet
import           Data.Alphabet.Codec
import           Data.BitVector.LittleEndian
import           Data.Bits
import           Data.Foldable
import           Data.Ord
import           Data.STRef
import           Data.Set                    (Set)
import qualified Data.Vector                 as V
import           Data.Vector.Generic         (Mutable, Vector, unsafeFreeze, (!))
import qualified Data.Vector.Generic         as GV
import           Data.Vector.Generic.Mutable (unsafeNew, unsafeRead, unsafeWrite)
import qualified Data.Vector.Generic.Mutable as GMV (length)
import qualified Data.Vector.Storable        as SV
import qualified Data.Vector.Unboxed         as UV
import           Data.Word
import           Foreign.C.Types


-- |
-- Encoding for a dynamic character element with an alphabet size of /8/ or less.
--
-- /NOTE:/ This encoding uses more bits than required! This is due to the C FFI
-- implementation details. It would be possible to reduce this to a 'CUChar' if
-- and only iff the C interface and implementation is updated.
type SlimState = CUInt


-- |
-- Encoding for a dynamic character element with an alphabet size in the range [9, 64].
type WideState = Word64


-- |
-- Encoding for a dynamic character element with an alphabet size of /65/ or greater.
type HugeState = BitVector


-- |
-- This triple of vectors is the general structure of all dynamic character
-- representations. Each vector is equal length.
--
-- A triple of empty vectors represents a missing character. This can be queried
-- by 'isMissing'.
--
-- All triples are arranged thusly:
--
--   * 1: ``Left'' character. When normalized, this will be the shorter character.
--        Delete events leave a void in this character.
--
--   * 2: Median character. The alignment of both the Left and Right characters.
--        Align, Delete, & Insert do *not* leave a void. This vector will always
--        contain states with one or more bits set, i.e. never contain voids.
--
--   * 3: ``Right'' character. When normalized, this will be the longer character.
--        Insert events leave a void in this character.
--
type OpenDynamicCharacter v e = ( v e, v e, v e )


-- |
-- Encoding for dynamic characters with an alphabet size of /8/ or less.
type SlimDynamicCharacter = OpenDynamicCharacter SV.Vector SlimState

-- |
-- Encoding for dynamic characters with an alphabet size in the range [9, 64].
type WideDynamicCharacter = OpenDynamicCharacter UV.Vector WideState


-- |
-- Encoding for dynamic characters with an alphabet size of /65/ or greater.
type HugeDynamicCharacter = OpenDynamicCharacter  V.Vector HugeState


-- |
-- Generic representation of a /mutable/ dynamic character.
type TempOpenDynamicCharacter m v e = OpenDynamicCharacter (Mutable v (PrimState m)) e


-- |
-- Mutable encoding of 'SlimDynamicCharacter'.
type TempSlimDynamicCharacter m = TempOpenDynamicCharacter m SV.Vector SlimState


-- |
-- Mutable encoding of 'WideDynamicCharacter'.
type TempWideDynamicCharacter m = TempOpenDynamicCharacter m UV.Vector WideState


-- |
-- Mutable encoding of 'HugeDynamicCharacter'.
type TempHugeDynamicCharacter m = TempOpenDynamicCharacter m  V.Vector HugeState


isAlign, isDelete, isInsert, isGapped :: (FiniteBits e, Vector v e) => OpenDynamicCharacter v e -> Int -> Bool
{-# INLINEABLE isAlign #-}
{-# SPECIALISE isAlign  :: SlimDynamicCharacter -> Int -> Bool #-}
{-# SPECIALISE isAlign  :: WideDynamicCharacter -> Int -> Bool #-}
{-# SPECIALISE isAlign  :: HugeDynamicCharacter -> Int -> Bool #-}
{-# INLINEABLE isDelete #-}
{-# SPECIALISE isDelete :: SlimDynamicCharacter -> Int -> Bool #-}
{-# SPECIALISE isDelete :: WideDynamicCharacter -> Int -> Bool #-}
{-# SPECIALISE isDelete :: HugeDynamicCharacter -> Int -> Bool #-}
{-# INLINEABLE isInsert #-}
{-# SPECIALISE isInsert :: SlimDynamicCharacter -> Int -> Bool #-}
{-# SPECIALISE isInsert :: WideDynamicCharacter -> Int -> Bool #-}
{-# SPECIALISE isInsert :: HugeDynamicCharacter -> Int -> Bool #-}
{-# INLINEABLE isGapped #-}
{-# SPECIALISE isGapped :: SlimDynamicCharacter -> Int -> Bool #-}
{-# SPECIALISE isGapped :: WideDynamicCharacter -> Int -> Bool #-}
{-# SPECIALISE isGapped :: HugeDynamicCharacter -> Int -> Bool #-}
isAlign  (lc,_,rc) i = i < GV.length lc && popCount (lc ! i) /= 0 && popCount (rc ! i) /= 0
isDelete (lc,_,rc) i = i < GV.length lc && popCount (lc ! i) == 0 && popCount (rc ! i) /= 0
isInsert (lc,_,rc) i = i < GV.length lc && popCount (lc ! i) /= 0 && popCount (rc ! i) == 0
isGapped (lc,_,rc) i = i < GV.length lc && popCount (lc ! i) == 0 && popCount (rc ! i) == 0


{-# INLINEABLE isGap #-}
{-# SPECIALISE isGap :: SlimDynamicCharacter -> Int -> Bool #-}
{-# SPECIALISE isGap :: WideDynamicCharacter -> Int -> Bool #-}
{-# SPECIALISE isGap :: HugeDynamicCharacter -> Int -> Bool #-}
isGap :: (Bits e, Vector v e) => OpenDynamicCharacter v e -> Int -> Bool
isGap (_,mc,_) i = i < GV.length mc &&
    let val = mc ! i
        gap = buildGap val
    in  gap == val


{-# INLINEABLE isMissing #-}
{-# SPECIALISE isMissing :: SlimDynamicCharacter -> Bool #-}
{-# SPECIALISE isMissing :: WideDynamicCharacter -> Bool #-}
{-# SPECIALISE isMissing :: HugeDynamicCharacter -> Bool #-}
isMissing :: Vector v e => OpenDynamicCharacter v e -> Bool
isMissing (x,y,z) = GV.length x == 0 && GV.length y == 0 && GV.length z == 0


{-# INLINEABLE setFrom #-}
{-# SPECIALISE setFrom :: SlimDynamicCharacter -> TempSlimDynamicCharacter (ST s) -> Int -> Int -> ST s () #-}
{-# SPECIALISE setFrom :: WideDynamicCharacter -> TempWideDynamicCharacter (ST s) -> Int -> Int -> ST s () #-}
{-# SPECIALISE setFrom :: HugeDynamicCharacter -> TempHugeDynamicCharacter (ST s) -> Int -> Int -> ST s () #-}
setFrom
  :: ( PrimMonad m
     , Vector v e
     )
  => OpenDynamicCharacter v e  -- ^ source
  -> TempOpenDynamicCharacter m v e -- ^ destination
  -> Int -- ^ Index to read from source
  -> Int -- ^ Index to write to destination
  -> m ()
setFrom (slc,smc,src) (dlc,dmc,drc) i j =
    unsafeWrite dlc j (slc ! i) *> unsafeWrite dmc j (smc ! i) *> unsafeWrite drc j (src ! i)


{-# INLINEABLE setAlign #-}
{-# SPECIALISE setAlign :: TempSlimDynamicCharacter (ST s) -> Int -> SlimState -> SlimState -> SlimState -> ST s () #-}
{-# SPECIALISE setAlign :: TempWideDynamicCharacter (ST s) -> Int -> WideState -> WideState -> WideState -> ST s () #-}
{-# SPECIALISE setAlign :: TempHugeDynamicCharacter (ST s) -> Int -> HugeState -> HugeState -> HugeState -> ST s () #-}
setAlign
  :: ( PrimMonad m
     , Vector v e
     )
  => TempOpenDynamicCharacter m v e
  -> Int -- ^ Index to set
  -> e   -- ^ Aligned ``Left'' element
  -> e   -- ^ Median Element
  -> e   -- ^ Aligned ``Right'' Element
  -> m ()
setAlign (lc,mc,rc) i le me re =
    unsafeWrite lc i le *> unsafeWrite mc i me *> unsafeWrite rc i re


{-# INLINEABLE setDelete #-}
{-# SPECIALISE setDelete :: TempSlimDynamicCharacter (ST s) -> Int -> SlimState -> SlimState -> ST s () #-}
{-# SPECIALISE setDelete :: TempWideDynamicCharacter (ST s) -> Int -> WideState -> WideState -> ST s () #-}
{-# SPECIALISE setDelete :: TempHugeDynamicCharacter (ST s) -> Int -> HugeState -> HugeState -> ST s () #-}
setDelete
  :: ( Bits e
     , PrimMonad m
     , Vector v e
     )
  => TempOpenDynamicCharacter m v e -- ^ Modifiable character
  -> Int -- ^ Index to set
  -> e   -- ^ Deleted ``Right'' element
  -> e   -- ^ Median Element
  -> m ()
setDelete (lc,mc,rc) i me re =
    unsafeWrite lc i (me `xor` me) *> unsafeWrite mc i me *> unsafeWrite rc i re


{-# INLINEABLE setInsert #-}
{-# SPECIALISE setInsert :: TempSlimDynamicCharacter (ST s) -> Int -> SlimState -> SlimState -> ST s () #-}
{-# SPECIALISE setInsert :: TempWideDynamicCharacter (ST s) -> Int -> WideState -> WideState -> ST s () #-}
{-# SPECIALISE setInsert :: TempHugeDynamicCharacter (ST s) -> Int -> HugeState -> HugeState -> ST s () #-}
setInsert
  :: ( Bits e
     , PrimMonad m
     , Vector v e
     )
  => TempOpenDynamicCharacter m v e -- ^ Modifiable character
  -> Int -- ^ Index to set
  -> e   -- ^ Median Element
  -> e   -- ^ Inserted ``Left'' element
  -> m ()
setInsert (lc,mc,rc) i le me =
    unsafeWrite lc i le *> unsafeWrite mc i me *> unsafeWrite rc i (me `xor` me)


{-# INLINEABLE setGapped #-}
{-# SPECIALISE setGapped :: TempSlimDynamicCharacter (ST s) -> Int -> ST s () #-}
{-# SPECIALISE setGapped :: TempWideDynamicCharacter (ST s) -> Int -> ST s () #-}
{-# SPECIALISE setGapped :: TempHugeDynamicCharacter (ST s) -> Int -> ST s () #-}
setGapped
  :: ( Bits e
     , PrimMonad m
     , Vector v e
     )
  => TempOpenDynamicCharacter m v e -- ^ Modifiable character
  -> Int -- ^ Index to set
  -> m ()
setGapped (lc,mc,rc) i = do
    tmp <- unsafeRead mc i
    let (# gap, nil #) = buildGapAndNil tmp
    unsafeWrite lc i nil
    unsafeWrite mc i gap
    unsafeWrite rc i nil


{-# INLINEABLE transposeCharacter #-}
transposeCharacter :: OpenDynamicCharacter v e -> OpenDynamicCharacter v e
transposeCharacter (lc,mc,rc) = (rc,mc,lc)


{-# INLINEABLE characterLength #-}
characterLength :: Vector v e  => OpenDynamicCharacter v e -> Word
characterLength = toEnum . GV.length . extractMediansGapped


-- |
-- Extract the /ungapped/ medians of a dynamic character.
{-# INLINEABLE extractMedians #-}
{-# SPECIALISE extractMedians :: SlimDynamicCharacter -> SV.Vector SlimState #-}
{-# SPECIALISE extractMedians :: WideDynamicCharacter -> UV.Vector WideState #-}
{-# SPECIALISE extractMedians :: HugeDynamicCharacter ->  V.Vector HugeState #-}
extractMedians :: (FiniteBits e, Vector v e) => OpenDynamicCharacter v e -> v e
extractMedians (_,me,_)
    | GV.null me = me
    | otherwise  =
        let gap  = buildGap $me ! 0
        in  GV.filter (/=gap) me


-- |
-- Extract the left child's /ungapped/ medians used to construct the dynamic character.
{-# INLINEABLE extractMediansLeft #-}
{-# SPECIALISE extractMediansLeft :: SlimDynamicCharacter -> SV.Vector SlimState #-}
{-# SPECIALISE extractMediansLeft :: WideDynamicCharacter -> UV.Vector WideState #-}
{-# SPECIALISE extractMediansLeft :: HugeDynamicCharacter ->  V.Vector HugeState #-}
extractMediansLeft :: (FiniteBits e, Vector v e) => OpenDynamicCharacter v e -> v e
extractMediansLeft (lc,_,_)
    | GV.null lc = lc
    | otherwise  =
        let nil  = buildNil $ lc ! 0
        in  GV.filter (/=nil) lc


-- |
-- Extract the right child's /ungapped/ medians used to construct the dynamic character.
{-# INLINEABLE extractMediansRight #-}
{-# SPECIALISE extractMediansRight :: SlimDynamicCharacter -> SV.Vector SlimState #-}
{-# SPECIALISE extractMediansRight :: WideDynamicCharacter -> UV.Vector WideState #-}
{-# SPECIALISE extractMediansRight :: HugeDynamicCharacter ->  V.Vector HugeState #-}
extractMediansRight :: (FiniteBits e, Vector v e) => OpenDynamicCharacter v e -> v e
extractMediansRight (_,_,rc)
    | GV.null rc = rc
    | otherwise  =
        let nil  = buildNil $ rc ! 0
        in  GV.filter (/=nil) rc


-- |
-- Extract the /gapped/ medians of a dynamic character.
{-# INLINEABLE extractMediansGapped #-}
extractMediansGapped :: OpenDynamicCharacter v e -> v e
extractMediansGapped (_,me,_) = me


-- |
-- Extract the left child's /gapped/ medians used to construct the dynamic character.
{-# INLINEABLE extractMediansLeftGapped #-}
extractMediansLeftGapped :: OpenDynamicCharacter v e -> v e
extractMediansLeftGapped (lc,_,_) = lc


-- |
-- Extract the right child's /gapped/ medians used to construct the dynamic character.
{-# INLINEABLE extractMediansRightGapped #-}
extractMediansRightGapped :: OpenDynamicCharacter v e -> v e
extractMediansRightGapped (_,_,rc) = rc


{-# INLINEABLE encodeDynamicCharacter #-}
{-# SPECIALISE encodeDynamicCharacter :: (Foldable f, Foldable g, Ord s) => Alphabet s -> (Word -> SlimState) -> f (g s) -> SlimDynamicCharacter #-}
{-# SPECIALISE encodeDynamicCharacter :: (Foldable f, Foldable g, Ord s) => Alphabet s -> (Word -> WideState) -> f (g s) -> WideDynamicCharacter #-}
{-# SPECIALISE encodeDynamicCharacter :: (Foldable f, Foldable g, Ord s) => Alphabet s -> (Word -> HugeState) -> f (g s) -> HugeDynamicCharacter #-}
{-# SPECIALISE encodeDynamicCharacter :: (Foldable f, Ord s) => Alphabet s -> (Word -> SlimState) -> f (Set s) -> SlimDynamicCharacter #-}
{-# SPECIALISE encodeDynamicCharacter :: (Foldable f, Ord s) => Alphabet s -> (Word -> WideState) -> f (Set s) -> WideDynamicCharacter #-}
{-# SPECIALISE encodeDynamicCharacter :: (Foldable f, Ord s) => Alphabet s -> (Word -> HugeState) -> f (Set s) -> HugeDynamicCharacter #-}
{-# SPECIALISE encodeDynamicCharacter :: (Foldable f) => Alphabet String -> (Word -> SlimState) -> f (Set String) -> SlimDynamicCharacter #-}
{-# SPECIALISE encodeDynamicCharacter :: (Foldable f) => Alphabet String -> (Word -> WideState) -> f (Set String) -> WideDynamicCharacter #-}
{-# SPECIALISE encodeDynamicCharacter :: (Foldable f) => Alphabet String -> (Word -> HugeState) -> f (Set String) -> HugeDynamicCharacter #-}
{-# SPECIALISE encodeDynamicCharacter :: Ord s => Alphabet s -> (Word -> SlimState) -> V.Vector (Set s) -> SlimDynamicCharacter #-}
{-# SPECIALISE encodeDynamicCharacter :: Ord s => Alphabet s -> (Word -> WideState) -> V.Vector (Set s) -> WideDynamicCharacter #-}
{-# SPECIALISE encodeDynamicCharacter :: Ord s => Alphabet s -> (Word -> HugeState) -> V.Vector (Set s) -> HugeDynamicCharacter #-}
{-# SPECIALISE encodeDynamicCharacter :: Alphabet String -> (Word -> SlimState) -> V.Vector (Set String) -> SlimDynamicCharacter #-}
{-# SPECIALISE encodeDynamicCharacter :: Alphabet String -> (Word -> WideState) -> V.Vector (Set String) -> WideDynamicCharacter #-}
{-# SPECIALISE encodeDynamicCharacter :: Alphabet String -> (Word -> HugeState) -> V.Vector (Set String) -> HugeDynamicCharacter #-}
encodeDynamicCharacter
  :: ( Bits e
     , Foldable f
     , Foldable g
     , Ord s
     , Vector v e
     )
  => Alphabet s      -- ^ Alphabet of symbols
  -> (Word -> e)     -- ^ Constructor for an empty element, taking the alphabet size
  -> f (g s)         -- ^ Sequence of ambiguity groups of symbols
  -> OpenDynamicCharacter v e -- ^ Encoded dynamic character
encodeDynamicCharacter alphabet f sequenceOfSymbols = dynamicCharater
  where
    len = length sequenceOfSymbols

    dynamicCharater = unsafeCharacterBuiltByST (toEnum len) $ \char -> do
       iRef <- newSTRef 0

       let writeElement symbols =
             let encodedVal = encodeState alphabet f symbols
             in do i <- readSTRef iRef
                   setAlign char i encodedVal encodedVal encodedVal
                   modifySTRef iRef succ

       traverse_ writeElement sequenceOfSymbols


{-# INLINEABLE newTempCharacter #-}
{-# SPECIALISE newTempCharacter :: Word -> ST s (TempSlimDynamicCharacter (ST s)) #-}
{-# SPECIALISE newTempCharacter :: Word -> ST s (TempWideDynamicCharacter (ST s)) #-}
{-# SPECIALISE newTempCharacter :: Word -> ST s (TempHugeDynamicCharacter (ST s)) #-}
newTempCharacter
  :: ( Vector v e
     , PrimMonad m
     )
  => Word
  -> m (TempOpenDynamicCharacter m v e)
newTempCharacter n =
    let i = fromEnum n
    in  (,,) <$> unsafeNew i <*> unsafeNew i <*> unsafeNew i


{-# INLINEABLE freezeTempCharacter #-}
{-# SPECIALISE freezeTempCharacter :: TempSlimDynamicCharacter (ST s) -> ST s SlimDynamicCharacter #-}
{-# SPECIALISE freezeTempCharacter :: TempWideDynamicCharacter (ST s) -> ST s WideDynamicCharacter #-}
{-# SPECIALISE freezeTempCharacter :: TempHugeDynamicCharacter (ST s) -> ST s HugeDynamicCharacter #-}
freezeTempCharacter
  :: ( PrimMonad m
     ,  Vector v e
     )
  => TempOpenDynamicCharacter m v e
  -> m (OpenDynamicCharacter v e)
freezeTempCharacter (lc,mc,rc) =
    (,,) <$> unsafeFreeze lc <*> unsafeFreeze mc <*> unsafeFreeze rc


{-# INLINEABLE unsafeCharacterBuiltByST #-}
{-# SPECIALISE unsafeCharacterBuiltByST :: Word -> (forall s. TempSlimDynamicCharacter (ST s) -> ST s ()) -> SlimDynamicCharacter #-}
{-# SPECIALISE unsafeCharacterBuiltByST :: Word -> (forall s. TempWideDynamicCharacter (ST s) -> ST s ()) -> WideDynamicCharacter #-}
{-# SPECIALISE unsafeCharacterBuiltByST :: Word -> (forall s. TempHugeDynamicCharacter (ST s) -> ST s ()) -> HugeDynamicCharacter #-}
unsafeCharacterBuiltByST
  :: Vector v e
  => Word
  -> (forall s. TempOpenDynamicCharacter (ST s) v e -> ST s ())
  -> OpenDynamicCharacter v e
unsafeCharacterBuiltByST n f = runST $ do
    char <- newTempCharacter n
    f char
    freezeTempCharacter char


-- |
-- Allocated a buffer of specified size.
-- Uses function to generate character and return the final size.
-- Copies character from buffer to character of final size.
{-# INLINEABLE unsafeCharacterBuiltByBufferedST #-}
{-# SPECIALISE unsafeCharacterBuiltByBufferedST :: Word -> (forall s. TempSlimDynamicCharacter (ST s) -> ST s Word) -> SlimDynamicCharacter #-}
{-# SPECIALISE unsafeCharacterBuiltByBufferedST :: Word -> (forall s. TempWideDynamicCharacter (ST s) -> ST s Word) -> WideDynamicCharacter #-}
{-# SPECIALISE unsafeCharacterBuiltByBufferedST :: Word -> (forall s. TempHugeDynamicCharacter (ST s) -> ST s Word) -> HugeDynamicCharacter #-}
unsafeCharacterBuiltByBufferedST
  :: Vector v e
  => Word -- ^ Buffer length
  -> (forall s. TempOpenDynamicCharacter (ST s) v e -> ST s Word)
  -> OpenDynamicCharacter v e
unsafeCharacterBuiltByBufferedST b f = runST $ do
    buff <- newTempCharacter b
    char <- f buff >>= newTempCharacter
    copyCharacter buff char
    freezeTempCharacter char
  where
    copyCharacter src@(x,_,_) des@(y,_,_) =
        let m = GMV.length x
            n = GMV.length y
            o = m - n
        in  forM_ [ 0 .. n - 1 ] $ copyAt src des o

    copyAt (slc,smc,src) (dlc,dmc,drc) o i =
      let i' = o + i
      in  do  unsafeRead slc i' >>= unsafeWrite dlc i
              unsafeRead smc i' >>= unsafeWrite dmc i
              unsafeRead src i' >>= unsafeWrite drc i


renderDynamicCharacter
  :: ( FiniteBits e
     , Show e
     , Vector v e
     )
  => OpenDynamicCharacter v e
  -> String
renderDynamicCharacter (lc,mc,rc) = unlines
    [ "Character Length: " <> show (GV.length mc)
    , printVector lcStr
    , printVector mcStr
    , printVector rcStr
    ]
  where
    show' x | popCount x > 0 = show x
            | otherwise      = [voidC]
    voidC = '█'
    lcStr = show' <$> GV.toList lc
    mcStr = show' <$> GV.toList mc
    rcStr = show' <$> GV.toList rc
    eSize = length . maximumBy (comparing length) $ lcStr <> mcStr <> rcStr <> [[voidC]]
    pad s =
      let c | s == [voidC] = voidC
            | otherwise    = ' '
      in  replicate (eSize - length s) c <> s

    intercalate'     [] = []
    intercalate'    [x] = x
    intercalate' (x:xs) =
      let sep = case x of
                  e:_ | e == voidC -> [voidC]
                  _                -> " "
      in   x <> sep <> intercalate' xs

    printVector vec = "[ " <> intercalate' (pad <$> vec) <> " ]"


buildGap :: Bits e => e -> e
buildGap e = buildNil e `setBit` gapIndex


buildNil :: Bits e => e -> e
buildNil e = e `xor` e


buildGapAndNil :: Bits e => e -> (# e, e #)
buildGapAndNil e =
  let nil = buildNil e
      gap = nil `setBit` gapIndex
  in  (# gap, nil #)
