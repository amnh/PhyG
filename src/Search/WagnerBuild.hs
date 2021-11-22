{- |
Module      :  WagnerBuild.hs
Description :  Module specifying charcter-based Wagner tree building functions
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

module Search.WagnerBuild  ( wagnerTreeBuild
                           , rasWagnerBuild
                     ) where

import Types.Types
import qualified Data.Vector as V
import GeneralUtilities
import Debug.Trace
import qualified GraphOptimization.Traversals as T
import qualified Utilities.LocalGraph as LG
import qualified Data.Text.Lazy              as TL
import qualified Graphs.GraphOperations  as GO
import qualified GraphOptimization.PreOrderFunctions as PRE
import Utilities.Utilities as U
import qualified Data.List as L
import Data.Maybe
import qualified GraphOptimization.Medians as M
import qualified ParallelUtilities as PU
import Control.Parallel.Strategies


-- import qualified ParallelUtilities as PU --need instance for VerexInfo


-- | rasWagnerBuild generates a series of random addition sequences and then calls wagnerTreeBuild to construct them.
-- Does not filter by best, unique etc.  That happens with the select() command specified separately.
rasWagnerBuild :: GlobalSettings -> ProcessedData -> Int -> Int -> [PhylogeneticGraph]
rasWagnerBuild inGS inData seed numReplicates =
   if numReplicates == 0 then []
   else 
      let numLeaves = V.length $ fst3 inData 
          randomizedAdditionSequences = V.fromList <$> shuffleInt seed numReplicates [0..numLeaves - 1]

          -- "graph" of leaf nodes without any edges
          leafGraph = T.makeSimpleLeafGraph inData
          leafDecGraph = T.makeLeafGraph inData

          hasNonExactChars = U.getNumberNonExactCharacters (thd3 inData) > 0
      in
      trace ("Building " ++ (show numReplicates) ++ " character Wagner replicates")
      -- PU.seqParMap PU.myStrategy (wagnerTreeBuild inGS inData) randomizedAdditionSequences
      fmap (wagnerTreeBuild inGS inData leafGraph leafDecGraph numLeaves hasNonExactChars) randomizedAdditionSequences `using` PU.myParListChunkRDS 

-- | wagnerTreeBuild builds a wagner tree (Farris 1970--but using random addition seqeuces--not "best" addition) 
-- from a leaf addition sequence. Always produces a tree that can be converted to a soft/hard wired network
-- afterwards
-- basic procs is to add edges to unresolved tree
-- currently naive wrt candidate tree costs
wagnerTreeBuild :: GlobalSettings -> ProcessedData -> SimpleGraph -> DecoratedGraph -> Int -> Bool -> V.Vector Int -> PhylogeneticGraph
wagnerTreeBuild inGS inData leafSimpleGraph leafDecGraph  numLeaves hasNonExactChars additionSequence = 
   let rootHTU = (numLeaves, TL.pack $ "HTU" ++ (show numLeaves))
       nextHTU = (numLeaves + 1, TL.pack $ "HTU" ++ (show $ numLeaves + 1))

       edge0 = (numLeaves, (additionSequence V.! 0), 0.0)
       edge1 = (numLeaves, numLeaves + 1, 0.0)
       edge2 = (numLeaves + 1, (additionSequence V.! 1), 0.0)
       edge3 = (numLeaves + 1, (additionSequence V.! 2), 0.0)

       initialTree = LG.insEdges [edge0, edge1, edge2, edge3] $ LG.insNodes [rootHTU, nextHTU] leafSimpleGraph

       blockCharInfo = V.map thd3 $ thd3 inData

       -- initialFullyDecoratedTree = T.multiTraverseFullyLabelTree inGS inData initialTree 
       initialFullyDecoratedTree = PRE.preOrderTreeTraversal inGS (finalAssignment inGS) hasNonExactChars $ T.postDecorateTree initialTree leafDecGraph blockCharInfo numLeaves

       wagnerTree = recursiveAddEdgesWagner (V.drop 3 $ additionSequence) numLeaves (numLeaves + 2) inGS inData hasNonExactChars leafDecGraph initialFullyDecoratedTree 
   in
   -- trace ("Initial Tree:\n" ++ (LG.prettify initialTree) ++ "FDT at cost "++ (show $ snd6 initialFullyDecoratedTree) ++":\n" 
   --   ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph $ thd6 initialFullyDecoratedTree))
   wagnerTree


-- | recursiveAddEdgesWagner adds edges until 2n -1 (n leaves) vertices in graph
-- this tested by null additin sequence list
-- interface will change with correct final states--using post-order pass for now
recursiveAddEdgesWagner ::V.Vector Int -> Int -> Int -> GlobalSettings -> ProcessedData -> Bool -> DecoratedGraph -> PhylogeneticGraph -> PhylogeneticGraph 
recursiveAddEdgesWagner additionSequence numLeaves numVerts inGS inData hasNonExactChars leafDecGraph inGraph@(inSimple, inCost, inDecGraph, _, _, charInfoVV) =
   -- all edges/ taxa in graph
   -- trace ("To go " ++ (show additionSequence) ++ " verts " ++ (show numVerts)) (
   if null additionSequence then inGraph
   else
      -- edges/taxa to add, but not the edges that leads to outgroup--redundant with its sister edge
      let outgroupEdges = filter ((< numLeaves) . snd3) $ LG.out inDecGraph numLeaves
          edgesToInvade = (LG.labEdges inDecGraph) L.\\ outgroupEdges
          leafToAdd = V.head additionSequence
          candidateEditList = fmap (addTaxonWagner inGS numLeaves numVerts inGraph leafToAdd leafDecGraph) edgesToInvade
          minDelta = minimum $ fmap fst4 candidateEditList
          (_, nodeToAdd, edgesToAdd, edgeToDelete) = head $ filter  ((== minDelta). fst4) candidateEditList

          -- create new tree
          newSimple = LG.insEdges edgesToAdd $ LG.insNode nodeToAdd $ LG.delEdge edgeToDelete inSimple

          newSimple' = if V.length additionSequence == 1 then GO.rerootTree (outgroupIndex inGS) newSimple
                       else newSimple

          -- create fully labelled tree, if all taxa in do full multi-labelled for correct graph type
          newPhyloGraph = if (V.length additionSequence > 1) then PRE.preOrderTreeTraversal inGS (finalAssignment inGS) hasNonExactChars $ T.postDecorateTree newSimple' leafDecGraph charInfoVV numLeaves
                          else T.multiTraverseFullyLabelGraph inGS inData False False Nothing newSimple'   
      in
      recursiveAddEdgesWagner (V.tail additionSequence)  numLeaves (numVerts + 1) inGS inData hasNonExactChars leafDecGraph newPhyloGraph
      -- )

-- | addTaxonWagner adds a taxon (really edges) by 'invading' and edge, deleting that adege and creteing 3 more
-- to existing tree and gets cost (for now by postorder traversal--so wasteful but will be by final states later)
-- returns a tuple of the cost, node to add, edges to add, edge to delete
addTaxonWagner ::  GlobalSettings -> Int -> Int -> PhylogeneticGraph -> Int -> DecoratedGraph -> LG.LEdge EdgeInfo -> (VertexCost, LG.LNode TL.Text, [LG.LEdge Double], LG.Edge) 
addTaxonWagner inGS numLeaves numVerts inGraph@(inSimple, inCost, inDecGraph, _, _, charInfoVV) leafToAdd leafDecGraph targetEdge =
   let edge0 = (numVerts, leafToAdd, 0.0)
       edge1 = (fst3 targetEdge, numVerts, 0.0)
       edge2 = (numVerts, snd3 targetEdge, 0.0)
       newNode = (numVerts, TL.pack ("HTU" ++ (show numVerts)))

       -- full post order
       --newSimpleGraph =  LG.insEdges [edge0, edge1, edge2] $ LG.insNode newNode $ LG.delEdge (LG.toEdge targetEdge) inSimple
       --newCost = snd6 $ T.postDecorateTree newSimpleGraph leafDecGraph charInfoVV numLeaves

       -- heuristic delta
       delta = getDelta leafToAdd targetEdge inDecGraph charInfoVV
      
   in
   (delta, newNode, [edge0, edge1, edge2], LG.toEdge targetEdge)
   -- (newCost, newNode, [edge0, edge1, edge2], LG.toEdge targetEdge)


-- | getDelta estimates the delta in tree cost by adding a leaf taxon in Wagner build
-- must be DO for this--isolated leaves won't have IA
getDelta :: Int -> LG.LEdge EdgeInfo -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> VertexCost
getDelta leafToAdd (eNode, vNode, targetlabel) inDecGraph charInfoVV =
   let existingEdgeCost = minLength targetlabel
       leafToAddVertData = vertData $ fromJust $ LG.lab inDecGraph leafToAdd
       eNodeVertData = vertData $ fromJust $ LG.lab inDecGraph eNode
       vNodeVertData = vertData $ fromJust $ LG.lab inDecGraph vNode

       -- create edge union 'character' blockData
       -- filters gaps (True argument) because using DOm (as must) to add taxa not in IA framework
       -- based on final assignments
       edgeUnionVertData = M.createEdgeUnionOverBlocks True eNodeVertData vNodeVertData charInfoVV []

   in
   -- trace ("GD: " ++ (show edgeUnionVertData)) (
   if (LG.lab inDecGraph leafToAdd == Nothing) || (LG.lab inDecGraph eNode == Nothing) || (LG.lab inDecGraph vNode == Nothing) then error ("Missing label data for vertices")
   else 
      let dLeafENode = sum $ fmap fst $ V.zipWith3 (PRE.getBlockCostPairsFinal DirectOptimization) leafToAddVertData eNodeVertData charInfoVV
          dLeafVNode = sum $ fmap fst $ V.zipWith3 (PRE.getBlockCostPairsFinal DirectOptimization) leafToAddVertData vNodeVertData charInfoVV

          -- Use edge union data for delta to edge data
          dLeafEdgeUnionCost = sum $ fmap fst $ V.zipWith3 (PRE.getBlockCostPairsFinal DirectOptimization) leafToAddVertData edgeUnionVertData charInfoVV

          -- should be able to use existing information--but for now using this
          -- existingEdgeCost' = sum $ fmap fst $ V.zipWith3 (PRE.getBlockCostPairsFinal DirectOptimization) eNodeVertData vNodeVertData charInfoVV
      in
      -- trace ("Delta: " ++ (show (dLeafENode, dLeafVNode, existingEdgeCost)))
      -- dLeafENode + dLeafVNode - existingEdgeCost
      -- trace ("Delta: " ++ (show dLeafEdgeUnionCost))
      

      dLeafEdgeUnionCost
      -- )


