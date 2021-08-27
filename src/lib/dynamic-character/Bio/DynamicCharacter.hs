{-# LANGUAGE Strict #-}

module Bio.DynamicCharacter where

import           Data.BitVector.LittleEndian
import           Data.Bits
import qualified Data.Vector                 as V
import           Data.Vector.Generic         (Vector, (!))
import qualified Data.Vector.Generic         as GV
import qualified Data.Vector.Storable        as SV
import qualified Data.Vector.Unboxed         as UV
import           Data.Word
import           Foreign.C.Types

type SlimState = CUInt
type WideState = Word64
type HugeState = BitVector

type SlimDynamicCharacter = (SV.Vector SlimState, SV.Vector SlimState, SV.Vector SlimState)
type WideDynamicCharacter = (UV.Vector WideState, UV.Vector WideState, UV.Vector WideState)
type HugeDynamicCharacter = ( V.Vector HugeState,  V.Vector HugeState,  V.Vector HugeState)


isAlign, isDelete, isInsert, isGapped :: (FiniteBits a, Vector v a) => (v a, v a, v a) -> Int -> Bool
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
isAlign  (_,y,z) i = i < GV.length y && popCount (y ! i) /= 0 && popCount (z ! i) /= 0
isDelete (_,y,z) i = i < GV.length y && popCount (y ! i) == 0 && popCount (z ! i) /= 0
isInsert (_,y,z) i = i < GV.length y && popCount (y ! i) /= 0 && popCount (z ! i) == 0
isGapped (_,y,z) i = i < GV.length y && popCount (y ! i) /= 0 && popCount (z ! i) /= 0


{-# INLINEABLE isMissing #-}
{-# SPECIALISE isMissing :: SlimDynamicCharacter -> Bool #-}
{-# SPECIALISE isMissing :: WideDynamicCharacter -> Bool #-}
{-# SPECIALISE isMissing :: HugeDynamicCharacter -> Bool #-}
isMissing :: Vector v a => (v a, v a, v a) -> Bool
isMissing (x,y,z) = GV.length x == 0 && GV.length y == 0 && GV.length z == 0
