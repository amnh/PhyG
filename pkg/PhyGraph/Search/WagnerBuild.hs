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
                           , wagnerTreeBuild'
                           , rasWagnerBuild
                     ) where

import           Control.Parallel.Strategies
-- import qualified Data.List                           as L
import           Data.Maybe
import qualified Data.Text.Lazy                      as TL
import qualified Data.Vector                         as V
import           Debug.Trace
import           GeneralUtilities
import qualified GraphOptimization.Medians           as M
import qualified GraphOptimization.PreOrderFunctions as PRE
import qualified GraphOptimization.Traversals        as T
-- import qualified Graphs.GraphOperations              as GO
import qualified ParallelUtilities                   as PU
import           Types.Types
import qualified Utilities.LocalGraph                as LG
import           Utilities.Utilities                 as U
import qualified GraphOptimization.PostOrderSoftWiredFunctions as POSW


-- import qualified ParallelUtilities as PU --need instance for VerexInfo


-- | rasWagnerBuild generates a series of random addition sequences and then calls wagnerTreeBuild to construct them.
-- Does not filter by best, unique etc.  That happens with the select() command specified separately.
rasWagnerBuild :: GlobalSettings -> ProcessedData -> Int -> Int -> [PhylogeneticGraph]
rasWagnerBuild inGS inData rSeed numReplicates =
   if numReplicates == 0 then []
   else
      let numLeaves = V.length $ fst3 inData
          randomizedAdditionSequences = V.fromList <$> shuffleInt rSeed numReplicates [0..numLeaves - 1]

          -- "graph" of leaf nodes without any edges
          leafGraph = T.makeSimpleLeafGraph inData
          leafDecGraph = T.makeLeafGraph inData

          hasNonExactChars = U.getNumberSequenceCharacters (thd3 inData) > 0
      in
      trace ("\t\tBuilding " ++ (show numReplicates) ++ " character Wagner replicates")
      -- seqParMap better for high level parallel stuff
      -- PU.seqParMap PU.myStrategy (wagnerTreeBuild inGS inData) randomizedAdditionSequences
      -- zipWith (wagnerTreeBuild inGS inData leafGraph leafDecGraph numLeaves hasNonExactChars) randomizedAdditionSequences [0..numReplicates - 1] `using` PU.myParListChunkRDS
      PU.seqParMap rdeepseq (wagnerTreeBuild' inGS inData leafGraph leafDecGraph numLeaves hasNonExactChars) (zip randomizedAdditionSequences [0..numReplicates - 1])
      -- fmap (wagnerTreeBuild' inGS inData leafGraph leafDecGraph numLeaves hasNonExactChars) (zip randomizedAdditionSequences [0..numReplicates - 1]) `using` PU.myParListChunkRDS


-- | wagnerTreeBuild' is a wrapper around wagnerTreeBuild to allow for better parallation--(zipWith not doing so well?)
wagnerTreeBuild' :: GlobalSettings -> ProcessedData -> SimpleGraph -> DecoratedGraph -> Int -> Bool -> (V.Vector Int, Int) -> PhylogeneticGraph
wagnerTreeBuild' inGS inData leafSimpleGraph leafDecGraph  numLeaves hasNonExactChars (additionSequence, replicateIndex) =
   wagnerTreeBuild inGS inData leafSimpleGraph leafDecGraph  numLeaves hasNonExactChars additionSequence replicateIndex 

-- | wagnerTreeBuild builds a wagner tree (Farris 1970--but using random addition seqeuces--not "best" addition)
-- from a leaf addition sequence. Always produces a tree that can be converted to a soft/hard wired network
-- afterwards
-- basic procs is to add edges to unresolved tree
-- currently naive wrt candidate tree costs
wagnerTreeBuild :: GlobalSettings -> ProcessedData -> SimpleGraph -> DecoratedGraph -> Int -> Bool -> V.Vector Int -> Int -> PhylogeneticGraph
wagnerTreeBuild inGS inData leafSimpleGraph leafDecGraph  numLeaves hasNonExactChars additionSequence replicateIndex =
   trace ("\tBuilding Wagner replicate " ++ (show replicateIndex)) (
   let rootHTU = (numLeaves, TL.pack $ "HTU" ++ (show numLeaves))
       nextHTU = (numLeaves + 1, TL.pack $ "HTU" ++ (show $ numLeaves + 1))

       edge0 = (numLeaves, (additionSequence V.! 0), 0.0)
       edge1 = (numLeaves, numLeaves + 1, 0.0)
       edge2 = (numLeaves + 1, (additionSequence V.! 1), 0.0)
       edge3 = (numLeaves + 1, (additionSequence V.! 2), 0.0)

       initialTree = LG.insEdges [edge0, edge1, edge2, edge3] $ LG.insNodes [rootHTU, nextHTU] leafSimpleGraph

       blockCharInfo = V.map thd3 $ thd3 inData

       -- initialFullyDecoratedTree = T.multiTraverseFullyLabelTree inGS inData initialTree
       -- False flag for staticIA--can't be done in build
       calculateBranchLengths = False -- must be True for delata using existing edge
       initialFullyDecoratedTree = PRE.preOrderTreeTraversal inGS (finalAssignment inGS) False calculateBranchLengths hasNonExactChars numLeaves False $ POSW.postDecorateTree False initialTree leafDecGraph blockCharInfo numLeaves numLeaves

       wagnerTree = recursiveAddEdgesWagner (V.drop 3 $ additionSequence) numLeaves (numLeaves + 2) inGS inData hasNonExactChars leafDecGraph initialFullyDecoratedTree
   in
   -- trace ("Initial Tree:\n" ++ (LG.prettify initialTree) ++ "FDT at cost "++ (show $ snd6 initialFullyDecoratedTree) ++":\n"
   --    ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph $ thd6 initialFullyDecoratedTree))
   wagnerTree
   )


-- | recursiveAddEdgesWagner adds edges until 2n -1 (n leaves) vertices in graph
-- this tested by null additin sequence list
-- interface will change with correct final states--using post-order pass for now
recursiveAddEdgesWagner ::V.Vector Int -> Int -> Int -> GlobalSettings -> ProcessedData -> Bool -> DecoratedGraph -> PhylogeneticGraph -> PhylogeneticGraph
recursiveAddEdgesWagner additionSequence numLeaves numVerts inGS inData hasNonExactChars leafDecGraph inGraph@(inSimple, _, inDecGraph, _, _, charInfoVV) =
   -- all edges/ taxa in graph
   -- trace ("To go " ++ (show additionSequence) ++ " verts " ++ (show numVerts)) (
   if null additionSequence then inGraph
   else
      -- trace ("RAEW-In: " ++ (show $ length additionSequence)) (
      -- edges/taxa to add, but not the edges that leads to outgroup--redundant with its sister edge
      let -- outgroupEdges = filter ((< numLeaves) . snd3) $ LG.out inDecGraph numLeaves
          edgesToInvade = (LG.labEdges inDecGraph) -- L.\\ outgroupEdges
          leafToAdd = V.head additionSequence

          -- since this is apporximate--can get a bit off
          -- candidateEditList = fmap (addTaxonWagner numVerts inGraph leafToAdd) edgesToInvade
          candidateEditList = PU.seqParMap rdeepseq (addTaxonWagner numVerts inGraph leafToAdd) edgesToInvade
          minDelta = minimum $ fmap fst4 candidateEditList
          (_, nodeToAdd, edgesToAdd, edgeToDelete) = head $ filter  ((== minDelta). fst4) candidateEditList

          -- create new tree
          newSimple = LG.insEdges edgesToAdd $ LG.insNode nodeToAdd $ LG.delEdge edgeToDelete inSimple

          -- this reroot since could add taxon sister to outgroup
          newSimple' = LG.rerootTree (outgroupIndex inGS) newSimple

          -- create fully labelled tree, if all taxa in do full multi-labelled for correct graph type
          -- False flag for static IA--can't do when adding in new leaves
          calculateBranchLengths = False -- must be True for delata using existing edge
          newPhyloGraph = -- T.multiTraverseFullyLabelTree inGS inData leafDecGraph (Just numLeaves) newSimple'
                          if (V.length additionSequence > 1) then PRE.preOrderTreeTraversal inGS (finalAssignment inGS) False calculateBranchLengths hasNonExactChars numLeaves False $ POSW.postDecorateTree False newSimple' leafDecGraph charInfoVV numLeaves numLeaves
                          else T.multiTraverseFullyLabelTree inGS inData leafDecGraph (Just numLeaves) newSimple'

      in
      {-
      let -- progress = takeWhile (/='.') $ show  ((fromIntegral (100 * (newVertexIndex - nOTUs))/fromIntegral (nOTUs - 2)) :: Double)
          (percentAdded, _) = divMod (100 * ((numLeaves - 2) - (V.length additionSequence))) (numLeaves - 2)
          (decileNumber, decileRemainder) = divMod percentAdded 10
          (_, oddRemainder) = divMod ((numLeaves - 2) - (V.length additionSequence)) 2
      in
      --trace (show (percentAdded, decileNumber, decileRemainder)) (
      
      if decileRemainder == 0 && oddRemainder == 0 then
          trace ("\t\t"++ (show $ 10 * decileNumber) ++ "%")
          recursiveAddEdgesWagner (V.tail additionSequence)  numLeaves (numVerts + 1) inGS inData hasNonExactChars leafDecGraph newPhyloGraph
      else 
      -}
      -- trace ("RAEW-Out: " ++ (show $ length additionSequence))
      recursiveAddEdgesWagner (V.tail additionSequence)  numLeaves (numVerts + 1) inGS inData hasNonExactChars leafDecGraph newPhyloGraph
      -- )

-- | addTaxonWagner adds a taxon (really edges) by 'invading' and edge, deleting that adege and creteing 3 more
-- to existing tree and gets cost (for now by postorder traversal--so wasteful but will be by final states later)
-- returns a tuple of the cost, node to add, edges to add, edge to delete
addTaxonWagner ::  Int
               -> PhylogeneticGraph
               -> Int
               -> LG.LEdge EdgeInfo
               -> (VertexCost, LG.LNode TL.Text, [LG.LEdge Double], LG.Edge)
addTaxonWagner numVerts (_, _, inDecGraph, _, _, charInfoVV) leafToAdd  targetEdge =
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
getDelta leafToAdd (eNode, vNode, _) inDecGraph charInfoVV =
   let leafToAddVertData = vertData $ fromJust $ LG.lab inDecGraph leafToAdd
       eNodeVertData = vertData $ fromJust $ LG.lab inDecGraph eNode
       vNodeVertData = vertData $ fromJust $ LG.lab inDecGraph vNode
       
       -- create edge union 'character' blockData
       -- filters gaps (True argument) because using DOm (as must) to add taxa not in IA framework
       -- based on final assignments
       -- not useing IA--False argument since no IA field for all leaves
       edgeUnionVertData = M.createEdgeUnionOverBlocks False True eNodeVertData vNodeVertData charInfoVV []

   in
   -- trace ("GD: " ++ (show edgeUnionVertData)) (
   if (LG.lab inDecGraph leafToAdd == Nothing) || (LG.lab inDecGraph eNode == Nothing) || (LG.lab inDecGraph vNode == Nothing) then error ("Missing label data for vertices")
   else
      let -- Use edge union data for delta to edge data
          dLeafEdgeUnionCost = sum $ fmap fst $ V.zipWith3 (PRE.getBlockCostPairsFinal DirectOptimization) leafToAddVertData edgeUnionVertData charInfoVV

          
          -- should be able to use existing information--but for now using this
          -- existingEdgeCost' = sum $ fmap fst $ V.zipWith3 (PRE.getBlockCostPairsFinal DirectOptimization) eNodeVertData vNodeVertData charInfoVV
      in
      -- trace ("Delta: " ++ (show (dLeafENode, dLeafVNode, existingEdgeCost)))
      -- dLeafENode + dLeafVNode - existingEdgeCost
      -- trace ("Delta: " ++ (show dLeafEdgeUnionCost) ++ " vs " ++ (show dLeafEVAddCost))

      -- min dLeafEdgeUnionCost dLeafEVAddCost
      -- dLeafEVAddCost
      dLeafEdgeUnionCost
      -- )


