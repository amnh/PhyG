{- |
Module      :  LocalSequence
Description :  Functions that map Data.Sequence to list-like functtions (head, tail etc)
Copyright   :  (c) 2014 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
License     :  

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.

Maintainer  :  Ward Wheeler <wheeler@amnh.org>
Stability   :  unstable
Portability :  portable (I hope)

-}

{-# Language ImportQualifiedPost #-} 

module LocalSequence
    ( module LocalSequence
    ) where

import Data.Sequence ((<|), (|>), (><))
import Data.Sequence qualified as S
import Data.Vector qualified as V
import Data.Foldable qualified as F

-- | sequjence type for exporting 
type Seq = S.Seq

-- | head maps to Take 
head :: Seq a -> a
head inSeq = S.index inSeq 0

-- | tail maps to drop Seq
tail :: Seq a -> Seq a
tail inSeq = S.drop 1 inSeq

-- | (!!) index
(!) :: Seq a -> Int -> a
(!) inSeq index = S.index inSeq index

-- | cons maps to (<|)
cons :: a -> Seq a -> Seq a
cons newElem inSeq = newElem <| inSeq

-- | snoc maps to (|>)
snoc ::Seq a ->  a -> Seq a
snoc inSeq newElem = inSeq |> newElem

-- | snoc with args reversed
snocFlip :: a -> Seq a -> Seq a
snocFlip newElem inSeq = inSeq |> newElem

-- | empty to empty
empty :: Seq a
empty = S.empty

-- | null equal to empty
null :: (Eq a) => Seq a -> Bool
null inSeq = if inSeq == S.empty then True
             else False

-- | singleton to singleton
singleton :: a -> Seq a
singleton newElem = S.singleton newElem

-- | ++ maps to ><
(++) :: Seq a -> Seq a -> Seq a
(++) inSeqA inSeqB = inSeqA >< inSeqB

-- | concat fold over ><
concat :: (Eq a) => Seq (Seq a) -> Seq a
concat inSeqSeq = concatInternal inSeqSeq LocalSequence.empty

-- | concatInternal internal concat function with accumulator
concatInternal :: (Eq a) => Seq (Seq a) -> Seq a -> Seq a
concatInternal inSeqSeq newSeq
    | LocalSequence.null inSeqSeq = newSeq
    | otherwise =
        let firstSeq = LocalSequence.head inSeqSeq
        in  concatInternal (LocalSequence.tail inSeqSeq) (firstSeq >< newSeq)

-- | zip maps to zip
zip :: Seq a -> Seq b -> Seq (a,b)
zip inSeqA inSeqB = S.zip inSeqA inSeqB

-- | length maps to length
length :: Seq a -> Int
length inSeq = S.length inSeq

-- | toList from Foldable
toList :: Seq a -> [a]
toList inSeq = F.toList inSeq

-- | fromList from fromList
fromList :: [a] -> Seq a 
fromList aList  = S.fromList aList

-- | toVector via intemediate List (alas) 
toVector :: Seq a -> V.Vector a
toVector inSeq = V.fromList $ toList inSeq

-- | toVector via intemediate List (alas) 
fromVector :: V.Vector a -> Seq a
fromVector inVect = S.fromList $ V.toList inVect 

-- | reverse to reverse
reverse :: Seq a -> Seq a
reverse inSeq = S.reverse inSeq

-- | last should be connstant time
last :: Seq a -> a
last inSeq = S.index inSeq $ (S.length inSeq) - 1

-- | map to fmap for ease of migrating libraries
map :: Traversable t => (a->b) -> t a -> t b
map f = fmap f

-- | drop maps to drop
drop :: Int -> Seq a -> Seq a
drop number inSeq = S.drop number inSeq

-- | take maps to take
take :: Int -> Seq a -> Seq a
take number inSeq = S.take number inSeq

-- | unsafeTake maps to take
unsafeTake ::  Int -> Seq a -> Seq a
unsafeTake number inSeq = S.take number inSeq

-- | unsafeDrop maps to drop
unsafeDrop :: Int -> Seq a -> Seq a
unsafeDrop number inSeq = S.drop number inSeq

-- | replicate maps to replicate
replicate ::  Int -> a -> Seq a 
replicate number value = S.replicate number value
