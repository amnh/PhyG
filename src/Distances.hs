{- |
Module      :  Distances.hs
Description :  Module specifying data types
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

module Distances (  distance2
                , getPairwiseDistances
                , getBlockDistance
                , getPairwiseBlocDistance
                ) where

import           Types
import           Debug.Trace
import qualified Medians as M

-- | distance2 gets distance cost between two vertices
-- leaf or otherwise in Processsed data (data in blocks and bit coded)
distance2 :: ProcessedData -> Int -> Int -> VertexCost
distance2 inData firstIndex secondIndex = 0.0



-- | getPairwiseDistances takes Processed data
-- and retuns a matrix (list of lists of Double) of pairwise
-- distances among vertices in data set over blocks ans all cahracter types
getPairwiseDistances :: ProcessedData ->  [[VertexCost]]
getPairwiseDistances inData = [[0.0]]


-- | getBlockDistance takes Block data and returns distance between
-- vertices based on block data
getBlockDistance :: BlockData -> Int -> Int -> VertexCost
getBlockDistance inData firstIndex secondIndex = 0.0

-- | getPairwiseBlocDistance returns pairwsie ditances among vertices for 
-- a block of data
getPairwiseBlocDistance :: BlockData -> [[VertexCost]]
getPairwiseBlocDistance inData = [[0.0]]