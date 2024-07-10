{- |
Types for encoding distances within a graph.
-}
module Types.DistanceTypes
  ( module Types.DistanceTypes,
  ) where

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
