{- |
Module      :  ParallelUtilities.hs
Description :  Utilities for parallel traversals, and other related functions
Copyright   :  (c) 2021 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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

{-# Options_GHC -fno-warn-orphans #-}

module ParallelUtilities
    ( parmap
    , seqParMap
    , getNumThreads
    , rnf
    , myStrategy
    , myStrategyRS
    , myStrategyRPAR
    , myParListChunk
    , myParListChunkRDS
    , myChunkParMapRDS
    ) where

import Control.Concurrent
import Control.DeepSeq
import Control.Parallel.Strategies
import Data.BitVector qualified as BV
import System.IO.Unsafe


-- Map a function over a traversable structure in parallel
-- Preferred over parMap which is limited to lists
-- Add chunking (with arguement) (via chunkList) "fmap blah blah `using` parListChunk chunkSize rseq/rpar"
-- but would have to do one for lists (with Chunk) and one for vectors  (splitAt recusively)
parmap :: Traversable t => Strategy b -> (a -> b) -> t a -> t b
parmap strat f = withStrategy (parTraversable strat).fmap f

-- | seqParMap takes strategy,  if numThread == 1 retuns fmap otherwise parmap and
seqParMap :: Traversable t => Strategy b -> (a -> b) -> t a -> t b
seqParMap strat f =
  if getNumThreads > 1 then parmap strat f
  else fmap f

myParListChunk :: Strategy a -> Strategy [a] 
myParListChunk localStrategy = parListChunk getNumThreads localStrategy

myParListChunkRDS :: (NFData a) => Strategy [a] 
myParListChunkRDS = parListChunk getNumThreads myStrategy

myStrategy :: (NFData b) => Strategy b
myStrategy = rdeepseq

myStrategyRS :: Strategy b
myStrategyRS = rseq

myStrategyRPAR :: Strategy b
myStrategyRPAR = rpar

-- | getNumThreads gets number of COncurrent  threads
{-# NOINLINE getNumThreads #-}
getNumThreads :: Int
getNumThreads = unsafePerformIO getNumCapabilities


-- NFData instance for parmap/rdeepseq Bit Vectory types
instance NFData BV.BV where
  
    rnf bv = BV.size bv `seq` BV.nat bv `seq` ()


-- | myChunkParMapRDS chuncked parmap that defaults to fmap if not paralell
myChunkParMapRDS :: NFData c => (b -> c) -> [b] -> [c]
myChunkParMapRDS f inList = 
  if getNumThreads == 1 then fmap f inList
  else fmap f inList `using` myParListChunkRDS
