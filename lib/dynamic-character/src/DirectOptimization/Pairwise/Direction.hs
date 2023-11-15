{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE Strict #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UnboxedTuples #-}

{- |
Module      :  DirectOptimization.Pairwise.Direction
Copyright   :  (c) 2015-2021 Ward Wheeler
License     :  BSD-style

Maintainer  :  wheeler@amnh.org
Stability   :  provisional
Portability :  portable

Direct optimization pairwise alignment using the Needleman-Wunsch algorithm.
These functions will allocate an M * N matrix.
-}
module DirectOptimization.Pairwise.Direction (
    Direction (..),

    -- * Querying
    minimumCostDirection,

    -- * Rendering
    boldDirection,
) where

import Data.Int
import Data.Vector.Generic qualified as GV
import Data.Vector.Generic.Mutable qualified as MGV
import Data.Vector.Primitive qualified as PV
import Data.Vector.Unboxed qualified as UV
import Data.Word


{- |
Which direction to align the character at a given matrix point.

It should be noted that the ordering of the three arrow types are important,
as it guarantees that the derived 'Ord' instance will have the following
property:

DiagArrow < LeftArrow < UpArrow

This means:

  - DiagArrow has highest precedence when one or more costs are equal

  - LeftArrow has second highest precedence when one or more costs are equal

  -   UpArrow has lowest precedence when one or more costs are equal

Using this 'Ord' instance, we can resolve ambiguous transformations in a
deterministic way. Without loss of generality in determining the ordering,
we choose the same biasing as the C code called from the FFI for consistency.
-}
data Direction = DiagArrow | LeftArrow | UpArrow
--data Direction = UpArrow | LeftArrow | DiagArrow 
    deriving stock (Eq, Ord)


newtype instance UV.MVector s Direction = MV_Direction (PV.MVector s Word8)


newtype instance UV.Vector Direction = V_Direction (PV.Vector Word8)


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
    show UpArrow = "↑"


{-# INLINEABLE boldDirection #-}
boldDirection ∷ Char → Char
boldDirection '↖' = '⇖'
boldDirection '←' = '⇐'
boldDirection '↑' = '⇑'
boldDirection d = d


{- |
Given the cost of deletion, alignment, and insertion (respectively), selects
the least costly direction. In the case of one or more equal costs, the
direction arrows are returned in the following descending order of priority:

  [ DiagArrow, LeftArrow, UpArrow ]
-}
{-# INLINEABLE minimumCostDirection #-}
{-# SPECIALIZE INLINE minimumCostDirection ∷ Int → Int → Int → (# Int, Direction #) #-}
{-# SPECIALIZE INLINE minimumCostDirection ∷ Int8 → Int8 → Int8 → (# Int8, Direction #) #-}
{-# SPECIALIZE INLINE minimumCostDirection ∷ Int16 → Int16 → Int16 → (# Int16, Direction #) #-}
{-# SPECIALIZE INLINE minimumCostDirection ∷ Int32 → Int32 → Int32 → (# Int32, Direction #) #-}
{-# SPECIALIZE INLINE minimumCostDirection ∷ Int64 → Int64 → Int64 → (# Int64, Direction #) #-}
{-# SPECIALIZE INLINE minimumCostDirection ∷ Word → Word → Word → (# Word, Direction #) #-}
{-# SPECIALIZE INLINE minimumCostDirection ∷ Word8 → Word8 → Word8 → (# Word8, Direction #) #-}
{-# SPECIALIZE INLINE minimumCostDirection ∷ Word16 → Word16 → Word16 → (# Word16, Direction #) #-}
{-# SPECIALIZE INLINE minimumCostDirection ∷ Word32 → Word32 → Word32 → (# Word32, Direction #) #-}
{-# SPECIALIZE INLINE minimumCostDirection ∷ Word64 → Word64 → Word64 → (# Word64, Direction #) #-}
minimumCostDirection
    ∷ (Ord e)
    ⇒ e
    → e
    → e
    → (# e, Direction #)
minimumCostDirection delCost alnCost insCost
    | alnCost <= delCost =
        if alnCost <= insCost
            then (# alnCost, DiagArrow #)
            else (# insCost, UpArrow #)
    | delCost <= insCost = (# delCost, LeftArrow #)
    | otherwise = (# insCost, UpArrow #)


{-# INLINE fromDirection #-}
fromDirection ∷ Direction → Word8
fromDirection DiagArrow = 0
fromDirection LeftArrow = 1
fromDirection UpArrow = 2


{-# INLINE toDirection #-}
toDirection ∷ Word8 → Direction
toDirection 0 = DiagArrow
toDirection 1 = LeftArrow
toDirection _ = UpArrow
