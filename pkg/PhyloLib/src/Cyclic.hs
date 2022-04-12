{- |
Module      :  Cyclic.hs
Description :  Topological sorting involves removing nodes from the graph into a list to
    determine the order they can appear (let nodes be tasks and edges be 
    constraints). This algorithm can only work on Directed Acyclic Graphs. In
    this variation we do not save the nodes in an order, but if we cannot 
    remove all nodes from the graph then a topological sort isn't possible 
    implying the graph has a cycle. In a topological sort, leaf nodes (nodes 
    without a successor) would be last in the ordering. If we remove a leaf
    and the edges that connet to it from a graph, then another leaf must 
    remain. If no other leaves remain, then the graph cannot be topologically
    sorted which indicates the existence of a cycle.

    Source : https://gist.github.com/msanatan/7933189#file-cyclic-hs
-}

module Cyclic 
(leafNode,
hasLeaf,
delLeaf,
cyclic
) where

import Data.Graph.Inductive

{-  Topological sorting involves removing nodes from the graph into a list to
    determine the order they can appear (let nodes be tasks and edges be 
    constraints). This algorithm can only work on Directed Acyclic Graphs. In
    this variation we do not save the nodes in an order, but if we cannot 
    remove all nodes from the graph then a topological sort isn't possible 
    implying the graph has a cycle. In a topological sort, leaf nodes (nodes 
    without a successor) would be last in the ordering. If we remove a leaf
    and the edges that connet to it from a graph, then another leaf must 
    remain. If no other leaves remain, then the graph cannot be topologically
    sorted which indicates the existence of a cycle.
-}


{-  This method determines whether a nodes is a leaf. It receives a graph and a
    node and returns a boolean with True indicating it's a leaf node
-}
leafNode :: (DynGraph g) => g a b -> Node -> Bool
leafNode gr node = suc gr node == []


{-  This method determines whether the graph has a leaf. This is done by 
    testing the leaf_node condition on every node of the graph It receives a 
    graph and returns a boolean with True indicating there is a leaf node
-}
hasLeaf :: (DynGraph g) => g a b -> Bool
hasLeaf gr = checkNodes . nodes $ gr
    where
    checkNodes :: [Node] -> Bool
    checkNodes [] = False
    checkNodes (x:xs)  | leafNode gr x = True
                       | otherwise = checkNodes(xs)
                        

{-  This method deletes a leaf node from a graph. The edges connecting to the
    leaf node are deleted as well. 
-}
delLeaf :: (DynGraph g) => g a b -> g a b
delLeaf gr = delNode leaf gr'
    where
    leaf = head [x | x <- nodes gr, leafNode gr x]
    gr' = delEdges newEdges gr
    newLedges = inn gr leaf
    newEdges = [(x,y) | (x,y,_) <- newLedges]


{-  This method indicates whether a given graph has a cycle or not. If it does
    True is returned, False otherwise. If the graph has no nodes it's acyclic,
    if it has no leaf nodes then cyclic. If it does have leaf nodes, remove it
    and check again
-}
cyclic :: (DynGraph g) => g a b -> Bool
cyclic gr   | isEmpty gr = False
            | not $ hasLeaf gr = True
            | otherwise = cyclic $ delLeaf gr