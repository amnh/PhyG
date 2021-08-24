{-# LANGUAGE FlexibleContexts #-}

module DirectOptimization.PreOrder where

import           Control.Monad.ST
import           Data.Bits
import           Data.Vector.Generic (Vector, (!))
import qualified Data.Vector.Generic as GV
import           Data.Functor
import           Data.STRef


-- |
-- Faithful translation of Algorithm 8 (Non-root Node Alignment) from the
-- "Efficient Implied Alignment" paper found at <https://doi.org/10.1186/s12859-020-03595-2 10.1186/s12859-020-03595-2> 
preOrderLogic
  :: ( FiniteBits a
     , Vector v a
     , Vector v (a, a, a)
     )
  => Word
  -> Bool
  -> (v a, v a, v a) -- ^ Parent Final       Alignment
  -> (v a, v a, v a) -- ^ Parent Preliminary Context
  -> (v a, v a, v a) -- ^ Child  Preliminary Context
  -> (v a, v a, v a) -- ^ Child  Final       Alignment
preOrderLogic symbolCount isLeftChild pAlignment@(x,_,_) pContext cContext = GV.unzip3 cAlignment
  where
    wlog  = x ! 0
    zero  = wlog `xor` wlog
    gap   = (bit . fromEnum $ symbolCount - 1, zero, zero)
    paLen = lengthSeq pAlignment
    ccLen = lengthSeq cContext
    cAlignment = runST $ do
      j' <- newSTRef 0
      k' <- newSTRef 0
      let go  i = do
              k <- readSTRef k'
              if   k > ccLen || pAlignment `isGappedAt` i
              then pure gap
              else do
                  j <- readSTRef j'
                  modifySTRef j' succ
                  if    pAlignment `isAlignedAt` i
                    || (    isLeftChild && pAlignment `isDeletedAt`  i && pContext `isDeletedAt`  j)
                    || (not isLeftChild && pAlignment `isInsertedAt` i && pContext `isInsertedAt` j)
                  then modifySTRef k' succ $> cContext `indexSeq` k
                  else pure gap

      GV.generateM paLen go 


isAlignedAt, isInsertedAt, isDeletedAt, isGappedAt
  :: ( FiniteBits a
     , Vector v a
     )
  => (v a, v a, v a)
  -> Int
  -> Bool


isAlignedAt (_,y,z) i =
  i < GV.length y && popCount (y ! i) /= 0 && popCount (z ! i) /= 0


isDeletedAt (_,y,z) i =
  i < GV.length y && popCount (y ! i) == 0 && popCount (z ! i) == 0


isInsertedAt (_,y,z) i =
  i < GV.length y && popCount (y ! i) == 0 && popCount (z ! i) == 0


isGappedAt (_,y,z) i =
  i < GV.length y && popCount (y ! i) == 0 && popCount (z ! i) == 0


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
