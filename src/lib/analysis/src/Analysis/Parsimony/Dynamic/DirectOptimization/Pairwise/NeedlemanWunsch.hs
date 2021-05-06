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

{-# LANGUAGE ConstraintKinds  #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies     #-}

module Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.NeedlemanWunsch
  ( naiveDO
  , naiveDOMemo
  , overlap2
  ) where

import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Internal
import Bio.Character.Encodable
import Data.Bits
import Data.Matrix.NotStupid                                           (matrix)
import Data.MonoTraversable
import Data.Foldable.Custom
import Data.List.NonEmpty      (NonEmpty(..))
import Data.Semigroup
import Data.Semigroup.Foldable

-- |
-- Performs a naive direct optimization.
-- Takes in two characters to run DO on and a metadata object
-- Returns an assignment character, the cost of that assignment, the assignment
-- character with gaps included, the aligned version of the first input character,
-- and the aligned version of the second input character. The process for this
-- algorithm is to generate a traversal matrix, then perform a traceback.
{-# INLINE naiveDO #-}
{-# SPECIALISE naiveDO :: (Word -> Word -> Word) -> DynamicCharacter -> DynamicCharacter -> (Word, DynamicCharacter) #-}
naiveDO
  :: ( DOCharConstraint s
     , FiniteBits (Subcomponent (Element s))
     )
  => (Word -> Word -> Word)  -- ^ Structure defining the transition costs between character states
  -> s                       -- ^ First  dynamic character
  -> s                       -- ^ Second dynamic character
  -> (Word, s)               -- ^ The cost and resulting the alignment
naiveDO costStruct char1 char2 = directOptimization (overlap2 costStruct) char1 char2 createNeedlemanWunchMatrix


-- |
-- The same as 'naiveDO' except that the "cost structure" parameter is assumed to
-- be a memoized overlap function.
{-# INLINE naiveDOMemo #-}
{-# SPECIALISE naiveDOMemo :: OverlapFunction AmbiguityGroup -> DynamicCharacter -> DynamicCharacter -> (Word, DynamicCharacter) #-}
naiveDOMemo :: DOCharConstraint s
            => OverlapFunction (Subcomponent (Element s))
            -> s
            -> s
            -> (Word, s)
naiveDOMemo tcm char1 char2 = directOptimization tcm char1 char2 createNeedlemanWunchMatrix


-- |
-- Main function to generate a 'NeedlemanWunchMatrix'. Works as in Needleman-Wunsch,
-- but allows for multiple indel/replacement costs, depending on the symbol change
-- cost function. Also, returns the aligned parent characters, with appropriate
-- ambiguities, as the third of each tuple in the matrix.
--
-- Takes in two 'EncodableDynamicCharacter's and a 'OverlapFunction'. The first
-- character must be the longer of the two and is the top labeling of the matrix.
-- Returns a 'NeedlemanWunchMatrix'.
{-# INLINE createNeedlemanWunchMatrix #-}
{-# SPECIALISE createNeedlemanWunchMatrix ::OverlapFunction AmbiguityGroup ->  DynamicCharacter -> DynamicCharacter -> NeedlemanWunchMatrix #-}
createNeedlemanWunchMatrix
  :: DOCharConstraint s
  => OverlapFunction (Subcomponent (Element s))
  -> s
  -> s
  -> NeedlemanWunchMatrix
createNeedlemanWunchMatrix overlapFunction topChar leftChar = result
  where
    result             = matrix rows cols generatingFunction
    rows               = olength leftChar + 1
    cols               = olength topChar  + 1
    generatingFunction = needlemanWunschDefinition overlapFunction topChar leftChar result



-- |
-- Takes one or more elements of 'FiniteBits' and a symbol change cost function
-- and returns a tuple of a new character, along with the cost of obtaining that
-- character. The return character may be (or is even likely to be) ambiguous.
-- Will attempt to intersect the two characters, but will union them if that is
-- not possible, based on the symbol change cost function.
--
-- To clarify, the return character is an intersection of all possible least-cost
-- combinations, so for instance, if @ char1 == A,T @ and @ char2 == G,C @, and
-- the two (non-overlapping) least cost pairs are A,C and T,G, then the return
-- value is A,C,G,T.
{-# INLINE overlap #-}
{-# SPECIALISE overlap :: FiniteBits e => (Word -> Word -> Word) -> NonEmpty e -> (e, Word) #-}
{-# SPECIALISE overlap :: (Word -> Word -> Word) -> NonEmpty AmbiguityGroup -> (AmbiguityGroup, Word) #-}
overlap
  :: ( FiniteBits e
     , Foldable1 f
     , Functor f
     )
  => (Word -> Word -> Word) -- ^ Symbol change matrix (SCM) to determine cost
  -> f e                    -- ^ List of elements for of which to find the k-median and cost
  -> (e, Word)              -- ^ K-median and cost
overlap sigma xs = go size maxBound zero
  where
    (size, zero) = let wlog = getFirst $ foldMap1 First xs
                   in  (finiteBitSize wlog, wlog `xor` wlog)

    go 0 theCost bits = (bits, theCost)
    go i oldCost bits =
        let i' = i - 1
            newCost = sum' $ getDistance (toEnum i') <$> xs
            (minCost, bits') = case oldCost `compare` newCost of
                                 EQ -> (oldCost, bits `setBit` i')
                                 LT -> (oldCost, bits            )
                                 GT -> (newCost, zero `setBit` i')
        in go i' minCost bits'

    getDistance i b = go' size (maxBound :: Word)
      where
        go' :: Int -> Word -> Word
        go' 0 a = a
        go' j a =
          let j' = j - 1
              a' = if b `testBit` j' then min a $ sigma i (toEnum j') else a
          in  go' j' a'


-- |
-- Calculate the median between /two/ states.
{-# INLINE overlap2 #-}
{-# SPECIALISE overlap2 :: (Word -> Word -> Word) -> AmbiguityGroup -> AmbiguityGroup -> (AmbiguityGroup, Word) #-}
overlap2
  :: FiniteBits e -- EncodableStreamElement e {- , Show e -})
  => (Word -> Word -> Word)
  -> e
  -> e
  -> (e, Word)
overlap2 sigma char1 char2 = overlap sigma $ char1 :| [char2]

