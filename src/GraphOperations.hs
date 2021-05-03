{- |
Module      :  GraphOperations.hs
Description :  Module specifying data type medians
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

{--
TODO:

  Parallelize  median2Vect
--}

module GraphOperations ( ladderizeGraph
                       , fullyLabelGraph
                       , verifyTimeConsistency
                       , rerootGraph
                       ) where

import           Debug.Trace
import qualified Data.Vector as V
import qualified DOWrapper as DOW
import Types
import qualified LocalGraph as LG
import qualified Data.Text.Lazy as T
import GeneralUtilities
import qualified LocalSequence as LS
import qualified GraphFormatUtilities as GFU

-- | ladderizeGraph is a wrapper around ladderizeGraph' to allow for mapping with 
-- local nodelist
ladderizeGraph :: SimpleGraph -> SimpleGraph
ladderizeGraph inGraph = ladderizeGraph' inGraph (LG.nodes inGraph)



-- | ladderize takes an input graph and ensures/creates nodes
-- such that all vertices are (indegree, outdegree) (0,>0), (1,2) (2,1) (1,0)
ladderizeGraph' :: SimpleGraph -> [LG.Node] -> SimpleGraph
ladderizeGraph' inGraph nodeList = 
    if LG.isEmpty inGraph then LG.empty
    else if null nodeList then inGraph
    else 
        let -- these are roots, network, tree, leaf nodes
            okNodeDegrees = [(0,1),(0,2),(1,2),(2,1),(1,0)]
            firstNode = head nodeList
            (inEdgeList, outEdgeList) = LG.getInOutEdges inGraph firstNode
            inOutPairLength = (length inEdgeList, length outEdgeList)
        in
        --trace ("node " ++ (show firstNode) ++ " " ++ (show inOutPair)) (
        -- node ok to keep
        if inOutPairLength `elem` okNodeDegrees then ladderizeGraph' inGraph (tail nodeList)
        -- node edges need modification
        else 
          let newGraph = resolveNode inGraph firstNode (inEdgeList, outEdgeList) inOutPairLength
          in 
          ladderizeGraph' newGraph (LG.nodes newGraph)
        --)

-- | resolveNode takes a graph and node and inbound edgelist and outbound edge list
-- and converts node to one off (indeg, outdeg) (0,1),(0,2),(1,2),(2,1),(1,0)
-- this only resolves a single nodes edges at a time and then returns new graph
-- when more hase to be done--that will occur on lultiple passes through nodes.
-- perhaps not the most efficient, but only done once per input graph
resolveNode :: SimpleGraph -> LG.Node -> ([LG.LEdge Double], [LG.LEdge Double]) -> (Int, Int) -> SimpleGraph
resolveNode inGraph curNode inOutPair@(inEdgeList, outEdgeList) (inNum, outNum) =
  if LG.isEmpty inGraph then LG.empty
  else 
    --trace ("Resolveing " ++ show curNode) (
    let numNodes = length $ LG.nodes inGraph
    in
    -- leaf or simple tree node in outdegree 
    if outNum == 0 || outNum == 1 then
      let first2Edges = take 2 inEdgeList
          newNode = (numNodes , T.pack $ ("HTU" ++ (show numNodes)))
          newEdge1 = (fst3 $ head first2Edges, numNodes, 0.0 :: Double)
          newEdge2 = (fst3 $ last first2Edges, numNodes, 0.0 :: Double)
          newEdge3 = (numNodes, curNode, 0.0 :: Double)
          newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
      in
      newGraph 
    -- root or simple network indegree node
    else if inNum == 0 || inNum == 2 then 
      let first2Edges = take 2 outEdgeList
          newNode = (numNodes , T.pack $ ("HTU" ++ (show numNodes)))
          newEdge1 = (numNodes, snd3 $ head first2Edges, 0.0 :: Double)
          newEdge2 = (numNodes, snd3 $ last first2Edges, 0.0 :: Double)
          newEdge3 = (curNode, numNodes, 0.0 :: Double)
          newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
      in 
      newGraph
    else error ("This can't happen in resolveNode in/out edge lists don't need to be resolved " ++ show inOutPair)
    --)

-- | verifyTimeConsistency take a SimpleGraph and checks for time consistency
-- of network nodes to verify network nodes are not definately not coeval
verifyTimeConsistency :: SimpleGraph -> SimpleGraph
verifyTimeConsistency inGraph =
   if LG.isEmpty inGraph then error ("Input Graph is empty in verifyTimeConsistency")
   else inGraph
   
   {-
   errorWithoutStackTrace ("Graph violates time consistency")
   -}

-- | fullyLabelGraph takes an unlabelled "simple' graph, performs post and preorder passes to 
-- fully label the graph and return a PhylogeenticGraph
fullyLabelGraph :: GlobalSettings -> ProcessedData -> SimpleGraph -> PhylogeneticGraph
fullyLabelGraph inGS inData inGraph = 
    if LG.isEmpty inGraph then (LG.empty, 0.0, V.empty, V.empty, inData)
    else 
        (inGraph, 0.0, V.empty, V.empty, inData)


-- | rerootGraph takes a graph and reroots based on a vertex index (usually leaf outgroup)
-- if input is a forest then only roots the component that contains the vertex
-- makes use of existing code based on Sequence
rerootGraph :: Int -> LG.Gr a Double -> LG.Gr a Double
rerootGraph rerootIndex inGraph = inGraph
  {-
  if LG.isEmpty inGraph then inGraph
  else 
    let inNodes = LG.labNodes inGraph
        inEdges = LG.labEdges inGraph
        nOTUs = length $ fst3 $ GFU.splitVertexList inGraph
        redirectedEdgeList = LS.toList $ directEdges rerootIndex nOTUs True (LS.fromList inEdges)
    in
    LG.mkGraph inNodes redirectedEdgeList
    -}

  -- | directexEdges takes a Sequence of edges and outgrop index and directs the 
  -- edges (parent -> child vertices) based on that
directEdges ::Int -> Int -> Bool -> LS.Seq (LG.LEdge Double) -> LS.Seq (LG.LEdge Double)
directEdges vertIndex nLeaves isFirst inEdgeSeq
  | LS.null inEdgeSeq = LS.empty
  | isFirst = --to find out group edge order larger to smaller will have outgroup index second
    let outgroupEdge = getEdgeRoot vertIndex inEdgeSeq
        remainingEdgeVect = subtractSequence (LS.singleton outgroupEdge) inEdgeSeq
        (a,b,w) = orderEdge outgroupEdge
    in
    LS.cons (a,b,w) (directEdges a nLeaves False remainingEdgeVect)
  | vertIndex < nLeaves = LS.empty
  | otherwise = -- not outgroup but regular node, get two child edges
    let descdendantEdges = getChildEdges vertIndex nLeaves inEdgeSeq
        remainingEdgeVect = subtractSequence descdendantEdges inEdgeSeq
        newDescEdges = LS.map (directSingleEdge vertIndex) descdendantEdges
    in
    if LS.length newDescEdges /= 2 then error ("There should be 2 child edges for index " ++ show vertIndex ++ " and there are(is) " ++ show (LS.length newDescEdges) ++ " " ++ show newDescEdges)
    else
        let (_, bf, _) = LS.head newDescEdges
            (_, bs, _) = LS.last newDescEdges
            firstSubEdges = directEdges bf nLeaves False remainingEdgeVect
            remainingEdgeVect' = subtractSequence firstSubEdges remainingEdgeVect
            secondSubEdges = directEdges bs nLeaves False remainingEdgeVect'
        in
        (newDescEdges LS.++ (firstSubEdges LS.++ secondSubEdges))

-- | orderEdge takes an Edge and puts high index first then lower
orderEdge :: LG.LEdge b -> LG.LEdge b
orderEdge (a,b,w) =
  if a > b then (a,b,w)
  else (b,a,w)

-- | subtractSequence subtracts elements of sequence a from sequence b
-- is thins n^2 ?
-- edges are directed
subtractSequence :: (Eq a) => LS.Seq a -> LS.Seq a -> LS.Seq a
subtractSequence a b
  | LS.null a = b
  | LS.null b = LS.empty
  | otherwise =
    let firstB = LS.head b
        notFound = LS.notElem firstB a
    in
    if notFound then LS.cons firstB (subtractSequence a (LS.tail b))
    else subtractSequence a (LS.tail b)

-- | getChildEdges returns the two edges that are childre of a vertex
getChildEdges :: (Eq b) => Int -> Int -> LS.Seq (LG.LEdge b) -> LS.Seq (LG.LEdge b)
getChildEdges vertIndex nLeaves inEdgeVect
  | LS.null inEdgeVect = LS.empty
  | vertIndex < nLeaves = error ("Looking for child of leaf " ++ show (vertIndex, nLeaves))
  | otherwise =
    let (a,b,w) = LS.head inEdgeVect
    in
    if (a == vertIndex) || (b == vertIndex) then LS.cons (a,b,w) (getChildEdges vertIndex nLeaves (LS.tail inEdgeVect)) else getChildEdges vertIndex nLeaves (LS.tail inEdgeVect)

-- | directSingleEdge takes an Int and makes that 'e' and otehr vertex as 'u' in edge (e->u)
directSingleEdge :: (Show b) => Int -> LG.LEdge b -> LG.LEdge b
directSingleEdge nodeIndex (a,b,w)
  | a == nodeIndex = (a,b,w)
  | b == nodeIndex = (b,a,w)
  | otherwise = error ("Index " ++ show nodeIndex ++ " doesn't match edge " ++ show (a,b,w) ++ " in directSingleEdge")

-- | getEdges takes root Index and determines edges from root
getEdgeRoot :: Int -> LS.Seq (LG.LEdge Double) -> LG.LEdge Double
getEdgeRoot edgeIndex edgeSeq =
  if LS.null edgeSeq then ( -1 :: Int, -1 :: Int, (-1.0) :: Double)
  else
   let (eVect, uVect, _) = LS.head edgeSeq
   in
   if (eVect == edgeIndex) || (uVect == edgeIndex) then LS.head edgeSeq else getEdgeRoot edgeIndex (LS.tail edgeSeq)


