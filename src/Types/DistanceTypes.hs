{- |
Module      :  Types.hs
Description :  Module with Types distance tree construction methods dWag, Neightbor-Joining, UPGMA, and WPGMA
              -- but with added refinement based on 4-point metric
Copyright   :  (c) 2020 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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
module Types.DistanceTypes where

import Data.Vector qualified as V
import SymMatrix qualified as M


type Vertex = Int


type Weight = Double


type Edge = (Vertex, Vertex, Weight)


type Tree = (V.Vector Vertex, V.Vector Edge)


type TreeWithData = (String, Tree, Double, M.Matrix Double)


type SplitTreeData = (V.Vector Edge, V.Vector Edge, Double, V.Vector Edge, M.Matrix Double)


-- | emptyTree
emptyTree ∷ Tree
emptyTree = (V.empty, V.empty)


-- | emptyTreeWithData
emptyTreeWithData ∷ TreeWithData
emptyTreeWithData = ("()[];", (V.empty, V.empty), 0.0, M.empty)


-- | used for comparing tree costs that are Double
epsilon ∷ Double
epsilon = 0.000000000000001


-- | precision for branch lengths display
precision ∷ Int
precision = 8
