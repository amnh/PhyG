{-# LANGUAGE FlexibleContexts #-}

module DirectOptimization.PreOrder where

import           Bio.DynamicCharacter
import           Control.Monad
import           Control.Monad.ST
import           Data.Bits
import           Data.STRef
import           Data.Vector.Generic         (Vector, (!))
import qualified Data.Vector.Generic         as GV
import qualified Data.Vector.Generic.Mutable as MGV


-- |
-- Faithful translation of Algorithm 8 (Non-root Node Alignment) from the
-- "Efficient Implied Alignment" paper found at:
-- <https://doi.org/10.1186/s12859-020-03595-2 10.1186/s12859-020-03595-2 DOI>
{-# INLINEABLE preOrderLogic #-}
{-# SPECIALISE preOrderLogic :: Word -> Bool -> SlimDynamicCharacter -> SlimDynamicCharacter -> SlimDynamicCharacter -> SlimDynamicCharacter #-}
{-# SPECIALISE preOrderLogic :: Word -> Bool -> WideDynamicCharacter -> WideDynamicCharacter -> WideDynamicCharacter -> WideDynamicCharacter #-}
{-# SPECIALISE preOrderLogic :: Word -> Bool -> HugeDynamicCharacter -> HugeDynamicCharacter -> HugeDynamicCharacter -> HugeDynamicCharacter #-}
preOrderLogic
  :: ( FiniteBits a
     , Vector v a
     )
  => Word
  -> Bool
  -> (v a, v a, v a) -- ^ Parent Final       Alignment
  -> (v a, v a, v a) -- ^ Parent Preliminary Context
  -> (v a, v a, v a) -- ^ Child  Preliminary Context
  -> (v a, v a, v a) -- ^ Child  Final       Alignment
preOrderLogic symbolCount isLeftChild pAlignment@(x,_,_) pContext cContext@(xs,ys,zs)
  | isMissing cContext = mAlignment -- Missing case is all gaps
  | otherwise          = cAlignment
  where
    wlog  = x ! 0
    zero  = wlog `xor` wlog
    gap   = bit . fromEnum $ symbolCount - 1
    paLen = lengthSeq pAlignment
    ccLen = lengthSeq cContext

    mAlignment =
      let zeds = GV.replicate paLen zero
      in  (GV.replicate paLen gap, zeds, zeds)

    cAlignment = runST $ do
      j'  <- newSTRef 0
      k'  <- newSTRef 0
      tempAlign@(xs',ys',zs') <- (,,) <$> MGV.unsafeNew paLen
                                      <*> MGV.unsafeNew paLen
                                      <*> MGV.unsafeNew paLen

      let setAt k i = setAlign tempAlign i (xs ! k) (ys ! k) (zs ! k)

      forM_ [0 .. paLen - 1] $ \i -> do
          k <- readSTRef k'
          if   k > ccLen || pAlignment `isGapped` i
          then tempAlign `setGapped` i
          else do
              j <- readSTRef j'
              modifySTRef j' succ
              if    pAlignment `isAlign` i
                || (    isLeftChild && pAlignment `isDelete` i && pContext `isDelete` j)
                || (not isLeftChild && pAlignment `isInsert` i && pContext `isInsert` j)
              then tempAlign `setGapped` i
              else modifySTRef k' succ *> (k `setAt` i)

      (,,) <$> GV.basicUnsafeFreeze xs'
           <*> GV.basicUnsafeFreeze ys'
           <*> GV.basicUnsafeFreeze zs'


indexSeq
  :: Vector v a
  => (v a, v a, v a)
  -> Int
  -> (a, a, a)
indexSeq (x,y,z) i
  | i >= GV.length x = error "Indexing error in Implied Alignment logic"
  | otherwise        = (x ! i, y ! i, z ! i)


lengthSeq :: Vector v a => (v a, v a, v a) -> Int
lengthSeq (x,_,_) = GV.length x
