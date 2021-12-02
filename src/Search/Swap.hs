{- |
Module      :  Swap.hs
Description :  Module specifying graph swapping rearrangement functions
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

module Search.Swap  ( swapMaster
                    ) where

import Types.Types
import qualified ParallelUtilities       as PU
import Control.Parallel.Strategies
import GeneralUtilities
import qualified Graphs.GraphOperations  as GO
import qualified Utilities.LocalGraph    as LG
import Utilities.Utilities               as U
import Debug.Trace
import           Data.Char
import           Text.Read
import           Data.Maybe
import qualified GraphOptimization.Traversals as T
import qualified Data.Vector as V
import qualified GraphOptimization.PreOrderFunctions as PRE
import qualified Data.List as L
import qualified Data.Text.Lazy              as TL
import qualified GraphOptimization.Medians as M
import qualified Data.List as L

-- | buildArgList is the list of valid build arguments
swapArgList :: [String]
swapArgList = ["spr","tbr", "keep", "steepest", "all", "nni", "ia"]


-- | swapMaster processes and spawns the swap functions
swapMaster ::  [Argument] -> GlobalSettings -> ProcessedData -> Int -> [PhylogeneticGraph] -> [PhylogeneticGraph]
swapMaster inArgs inGS inData rSeed inGraphList =
   trace ("Swapping " ++ (show $ length inGraphList) ++ " graph(s)") (
   if null inGraphList then []
   else 
      let fstArgList = fmap (fmap toLower . fst) inArgs
          sndArgList = fmap (fmap toLower . snd) inArgs
          lcArgList = zip fstArgList sndArgList
          checkCommandList = U.checkCommandArgs "swap" fstArgList swapArgList
   in
   -- check for valid command options
   if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'swap': " ++ show inArgs)
   else 
       let keepList = filter ((=="keep").fst) lcArgList
           keepNum
            | length keepList > 1 =
              errorWithoutStackTrace ("Multiple 'keep' number specifications in swap command--can have only one: " ++ show inArgs)
            | null keepList = Just 1
            | otherwise = readMaybe (snd $ head keepList) :: Maybe Int
      in
      if isNothing keepNum then errorWithoutStackTrace ("Keep specification not an integer: "  ++ show (snd $ head keepList))
      else 
         let -- getting values to be passed for graph diagnosis later
             numLeaves = V.length $ fst3 inData
             leafGraph = T.makeSimpleLeafGraph inData
             leafDecGraph = T.makeLeafGraph inData
             leafGraphSoftWired = T.makeLeafGraphSoftWired inData
             hasNonExactChars = U.getNumberNonExactCharacters (thd3 inData) > 0
             charInfoVV = V.map thd3 $ thd3 inData

             -- process args for swap
             doNNI = any ((=="nni").fst) lcArgList
             doSPR' = any ((=="spr").fst) lcArgList
             doTBR = any ((=="tbr").fst) lcArgList
             doIA = any ((=="ia").fst) lcArgList
             doSteepest' = any ((=="steepest").fst) lcArgList
             doAll = any ((=="all").fst) lcArgList
             doSPR = if (not doNNI && not doSPR' && not doTBR) then True
                     else doSPR'
             doSteepest = if (not doSteepest' && not doAll) then True
                          else doSteepest'
             (newGraphList, counterNNI)  = if doNNI then 
                                             let graphPairList = fmap (swapSPRTBR "nni" inGS inData (fromJust keepNum) doSteepest numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA) inGraphList `using` PU.myParListChunkRDS
                                                 (graphListList, counterList) = unzip graphPairList
                                             in (concat graphListList, sum counterList)
                                           else (inGraphList, 0)
             (newGraphList', counterSPR)  = if doSPR then 
                                             let graphPairList = fmap (swapSPRTBR "spr" inGS inData (fromJust keepNum) doSteepest numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA) newGraphList `using` PU.myParListChunkRDS
                                                 (graphListList, counterList) = unzip graphPairList
                                             in (concat graphListList, sum counterList)
                                           else (newGraphList, 0)

             (newGraphList'', counterTBR) = if doTBR then 
                                             let graphPairList =  fmap (swapSPRTBR "tbr" inGS inData (fromJust keepNum) doSteepest numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA) newGraphList' `using` PU.myParListChunkRDS
                                                 (graphListList, counterList) = unzip graphPairList
                                             in (concat graphListList, sum counterList)
                                           else (newGraphList', 0)
            in
            trace ("Number after swap Graphs: " ++ (show $ length newGraphList'') ++ " Swap rounds: " ++ (show counterNNI) ++ " NNI and " ++ (show counterSPR) ++ " SPR and " ++ (show counterTBR) ++ " TBR")
            newGraphList''
   )


-- | swapSPRTBR perfomrs SPR or TBR branch (edge) swapping on graphs
-- runs both SPR and TBR depending on argument since so much duplicated functionality
-- 'steepest' abandons swap graph and switces to found graph as soon as anyhting 'better'
-- is found. The alternative (all) examines the entire neighborhood and retuns the best result
-- the retuns is a list of better graphs and the number of swapping rounds were required to ge there
swapSPRTBR  :: String 
            -> GlobalSettings 
            -> ProcessedData
            -> Int 
            -> Bool 
            -> Int
            -> SimpleGraph 
            -> DecoratedGraph 
            -> DecoratedGraph 
            -> Bool
            -> V.Vector (V.Vector CharInfo) 
            -> Bool
            -> PhylogeneticGraph 
            -> ([PhylogeneticGraph], Int)
swapSPRTBR swapType inGS inData numToKeep steepest numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA inGraph = 
   trace ("In swapSPRTBR:") (
   if LG.isEmpty (fst6 inGraph) then ([], 0)
   else 
      -- steepest takes immediate best
      if steepest then 
         let (swappedGraphs, counter) = swapSteepest swapType inGS inData numToKeep steepest 0 (snd6 inGraph) [] [inGraph] numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA
         in 
         (swappedGraphs, counter)

      -- All does all swaps before taking best
      else  
         trace ("Going into SwapAll") (
         let (swappedGraphs, counter) = swapAll swapType inGS inData numToKeep steepest 0 (snd6 inGraph) [] [inGraph] numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA
         in 
         trace ("SSPRTBR: " ++ (show counter))
         (swappedGraphs, counter)
         )
         )
      
-- | swapAll performs branch swapping on all 'break' edges and all readditions
-- edges are unsorted since doing all of them
swapAll  :: String 
         -> GlobalSettings 
         -> ProcessedData 
         -> Int 
         -> Bool 
         -> Int 
         -> VertexCost 
         -> [PhylogeneticGraph] 
         -> [PhylogeneticGraph] 
         -> Int 
         -> SimpleGraph 
         -> DecoratedGraph 
         -> DecoratedGraph 
         -> Bool
         -> V.Vector (V.Vector CharInfo) 
         -> Bool
         -> ([PhylogeneticGraph], Int)
swapAll swapType inGS inData numToKeep steepest counter curBestCost curSameBetterList inGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA =
   trace ("SWAPALLL Num: " ++ (show $ length inGraphList) ++ " returning: " ++ (show $ length curSameBetterList) ++ " counter: " ++ (show counter))  (
   if null inGraphList then trace ("Swap All Done") (curSameBetterList, counter)
   else 
      let firstGraph = head inGraphList
          firstDecoratedGraph = thd6 firstGraph
          (firstRootIndex, _) = head $ LG.getRoots firstDecoratedGraph

          -- filter out edges from root since no use--would just rejoin
          firstEdgeList = filter ((/= firstRootIndex) . fst3) $ LG.labEdges firstDecoratedGraph

          -- determine edges to break on--'bridge' edges only for network
          breakEdgeList = if (graphType inGS) == Tree then firstEdgeList
                          else GO.getEdgeSplitList firstDecoratedGraph

          -- create list of breaks
          (splitGraphList, graphRootList, prunedGraphRootIndexList,  originalConnectionOfPruned) = L.unzip4 $ fmap (GO.splitGraphOnEdge firstDecoratedGraph) breakEdgeList

          reoptimizedSplitGraphList = zipWith3 (reoptimizeGraphFromVertex inGS inData swapType doIA charInfoVV firstDecoratedGraph) splitGraphList graphRootList prunedGraphRootIndexList 

          -- create rejoins-- adds in break list so don't remake the initial graph
          -- didn't concatMap so can parallelize later
          -- this cost prob doesn't include the root/net penalty--so need to figure out
          swapPairList = concat $ L.zipWith5 (rejoinGraphKeepBest inGS swapType curBestCost numToKeep steepest doIA charInfoVV) reoptimizedSplitGraphList graphRootList prunedGraphRootIndexList originalConnectionOfPruned breakEdgeList

          -- keeps only "best" heuristic swap costs graphs
          minimumCandidateGraphCost = if (null swapPairList) then infinity
                                      else minimum $ fmap snd swapPairList
          candidateSwapGraphList = filter ((== minimumCandidateGraphCost). snd) swapPairList

          
          -- this should be incremental--full 2-pass for now]
          reoptimizedSwapGraphList = if (graphType inGS == Tree) then fmap (T.multiTraverseFullyLabelTree inGS inData leafDecGraph Nothing) (fmap fst candidateSwapGraphList)
                                     else if (graphType inGS == SoftWired) then fmap (T.multiTraverseFullyLabelSoftWired inGS inData False False leafGraphSoftWired Nothing) (fmap fst candidateSwapGraphList)
                                     else errorWithoutStackTrace "Hard-wired graph optimization not yet supported"

          -- selects best graph list based on full optimization
          bestSwapGraphList = GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] reoptimizedSwapGraphList

          bestSwapCost = if null swapPairList then infinity
                         else snd6 $ head bestSwapGraphList

      in
      trace ("Breakable Edges :" ++ (show $ fmap LG.toEdge breakEdgeList) ++ "\nIn graph:\n" ++ (LG.prettify $ fst6 firstGraph)) (
      
      -- either no better or more of same cost graphs
      trace ("BSG: " ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph $ thd6 $ head bestSwapGraphList)) (
      if bestSwapCost == curBestCost then 
         let newCurSameBestList = if firstGraph `notElem` curSameBetterList then (firstGraph : curSameBetterList)
                                  else curSameBetterList
             graphsToSwap = bestSwapGraphList L.\\ newCurSameBestList               
         in
         trace ("Same cost")
         swapAll swapType inGS inData numToKeep steepest (counter + 1) curBestCost newCurSameBestList ((tail inGraphList) ++ graphsToSwap) numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA

      -- better cost graphs
      else if (bestSwapCost < curBestCost) then 
         trace ("Better cost")
         swapAll swapType inGS inData numToKeep steepest (counter + 1) bestSwapCost [] bestSwapGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA

      -- didn't find equal or better graphs
      else 
         trace ("Worse cost")
         swapAll swapType inGS inData numToKeep steepest (counter + 1) curBestCost [firstGraph] (tail inGraphList) numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA
      )
      )
      )



-- | rejoinGraphKeepBest rejoins split trees on available edges (non-root, and not original split)
-- if steepest is False does not sort order of edges, other wise sorts in order of closeness to original edge
-- uses delta
-- NNI sorts edges on propinquity taking first 2 edges
-- TBR does the rerooting of pruned subtree
-- originalConnectionOfPruned is the "naked" node that was creted when teh graph was split and will 
-- be used for the rejoin node in the middle of th einvaded edge
rejoinGraphKeepBest :: GlobalSettings 
                    -> String 
                    -> VertexCost 
                    -> Int 
                    -> Bool 
                    -> Bool 
                    -> V.Vector (V.Vector CharInfo) 
                    -> (DecoratedGraph, VertexCost) 
                    -> LG.Node 
                    -> LG.Node 
                    -> LG.Node 
                    -> LG.LEdge EdgeInfo 
                    -> [(SimpleGraph, VertexCost)]
rejoinGraphKeepBest inGS swapType curBestCost numToKeep steepest doIA charInfoVV (splitGraph, splitCost) graphRoot prunedGraphRootIndex nakedNode originalSplitEdge = 
   -- case where swap split retunred empty ebcasue too few nodes in remaining graph to add to
   if LG.isEmpty splitGraph then []
   else
      let outgroupEdges = LG.out splitGraph graphRoot
          (_, prunedSubTreeEdges) = LG.nodesAndEdgesAfter splitGraph ([],[]) [(nakedNode, fromJust $ LG.lab splitGraph nakedNode)]
          edgesToInvade = (LG.labEdges splitGraph) L.\\ (outgroupEdges ++ prunedSubTreeEdges)
          candidateEditList = if (not steepest) then fmap (addSubGraph inGS doIA splitGraph prunedGraphRootIndex splitCost nakedNode charInfoVV) edgesToInvade 
                              else error "Steepest not yet implemented"

          minCandidateCost = if (not $ null candidateEditList) then minimum $ fmap fst4 candidateEditList   
                             else infinity
      in
      trace ("RGKB: " ++ (show $ fmap LG.toEdge edgesToInvade) ++ " " ++ (show curBestCost) ++ " v " ++ (show minCandidateCost)) (
      if minCandidateCost > curBestCost then []
      else 
         let bestEdits = filter ((== minCandidateCost). fst4) candidateEditList
             splitGraphSimple = GO.convertDecoratedToSimpleGraph splitGraph
             swapSimpleGraphList = fmap (applyGraphEdits splitGraphSimple) bestEdits
         in
         zip swapSimpleGraphList (L.replicate (length swapSimpleGraphList) minCandidateCost)
      )

-- | applyGraphEdits takes a  graphs and list of nodes and edges to add and delete and creates new graph
applyGraphEdits :: (Show a, Show b) => LG.Gr a b -> (VertexCost, LG.LNode a, [LG.LEdge b], LG.Edge) ->  LG.Gr a b
applyGraphEdits inGraph editStuff@(_, nodeToAdd, edgesToAdd, edgeToDelete) = 
   let editedGraph = LG.insEdges edgesToAdd $ LG.delEdge edgeToDelete inGraph
   in
   trace ("AGE: " ++ (show editStuff) ++ "\nIn graph:\n" ++ (LG.prettify inGraph) ++ "New Graph:\n" ++ (LG.prettify editedGraph)) 
   editedGraph
   


-- | addSubTree "adds" a subtree back into an edge calculating the cost of the graph via the delta of the add and costs of the two components
addSubGraph :: GlobalSettings 
            -> Bool 
            -> DecoratedGraph 
            -> LG.Node 
            -> VertexCost 
            -> LG.Node 
            -> V.Vector (V.Vector CharInfo) 
            -> LG.LEdge EdgeInfo 
            -> (VertexCost, LG.LNode TL.Text, [LG.LEdge Double], LG.Edge) 
addSubGraph inGS doIA inGraph prunedGraphRootIndex splitCost nakedNode charInfoVV targetEdge@(eNode, vNode, targetlabel) =  
   let existingEdgeCost = minLength targetlabel
       edge0 = (nakedNode, vNode, 0.0)
       edge1 = (eNode, nakedNode, 0.0)
       -- edge2 = (nakedNode, prunedGraphRootIndex, 0.0)
       newNode = (nakedNode, TL.pack ("HTU" ++ (show nakedNode)))
       delta = getSubGraphDelta targetEdge doIA inGraph prunedGraphRootIndex charInfoVV
   in
   
   -- do not redo origal edge so retun infinite cost and dummy edits
   if (eNode == nakedNode) then  trace ("ASG: break edge") (infinity, (-1, TL.empty), [], (-1,-1))
   else trace ("ASG: " ++ (show (eNode, vNode, nakedNode)))  (delta + splitCost, newNode, [edge0, edge1], (eNode, vNode)) -- edge 2
   


-- | getSubGraphDelta calcualted cost of adding a subgraph into and edge
getSubGraphDelta :: LG.LEdge EdgeInfo -> Bool -> DecoratedGraph -> LG.Node -> V.Vector (V.Vector CharInfo) -> VertexCost
getSubGraphDelta (eNode, vNode, targetlabel) doIA inGraph prunedGraphRootIndex charInfoVV = 
   let existingEdgeCost = minLength targetlabel
       subGraphVertData = vertData $ fromJust $ LG.lab inGraph prunedGraphRootIndex
       eNodeVertData = vertData $ fromJust $ LG.lab inGraph eNode
       vNodeVertData = vertData $ fromJust $ LG.lab inGraph vNode

       -- create edge union 'character' blockData
       -- based on final assignments--need to filter gaps if DO, not itIA
       edgeUnionVertData = M.createEdgeUnionOverBlocks (not doIA) eNodeVertData vNodeVertData charInfoVV []

       -- Use edge union data for delta to edge data
       costMethod = if doIA then ImpliedAlignment
                    else DirectOptimization

       subGraphEdgeUnionCost = sum $ fmap fst $ V.zipWith3 (PRE.getBlockCostPairsFinal costMethod) subGraphVertData edgeUnionVertData charInfoVV
   in
   subGraphEdgeUnionCost

-- | reoptimizeGraphFromVertex fully labels the component graph that is connected to the specified vertex
-- for softwired--need to deal with popocount at root
-- reooting issues for single component
-- need the cost to calculate the deltas later during rejoin--summed costs of the two comp[onets after splitting
-- doIA option to only do IA optimization as opposed to full thing--should be enormously faster--but yet more approximate
-- cretes finel for base graph but only does preorder pass fo component if TBR swap
reoptimizeGraphFromVertex :: GlobalSettings 
                          -> ProcessedData 
                          -> String 
                          -> Bool 
                          -> V.Vector (V.Vector CharInfo) 
                          -> DecoratedGraph 
                          -> DecoratedGraph 
                          -> Int 
                          -> Int 
                          -> (DecoratedGraph, VertexCost)
reoptimizeGraphFromVertex inGS inData swapType doIA charInfoVV origGraph inGraph startVertex prunedSubGraphRootVertex =

   trace ("RGFV: startVertex " ++ (show startVertex) ++ " prunedVertex " ++ (show prunedSubGraphRootVertex) ++ "\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inGraph)) (

   -- create graph of base (with ur-root) and pruned (non-ur-root) components
   let nonExactCharacters = U.getNumberNonExactCharacters (thd3 inData)
       leafGraph = LG.extractLeafGraph inGraph

       -- DO or IA for reoptimization for use of final sytartes later IA faster but more approximate
       (postOrderBaseGraph, localRootCost, localStartVertex) = if not doIA then T.generalizedGraphPostOrderTraversal inGS nonExactCharacters inData leafGraph (Just startVertex) (GO.convertDecoratedToSimpleGraph inGraph)
                            else error "IA reoptimizeGraphFromVertex not yet implemented"

       -- pruned component cost
       prunedSubGraphCost = vertexCost $ fromJust $ LG.lab inGraph prunedSubGraphRootVertex

       -- get pruned component nodes and edges
       parentPrunedNodeIndex = head $ LG.parents origGraph prunedSubGraphRootVertex

       -- set same label as for pruned node so that edge length is zero
       parentPrunedNode = (parentPrunedNodeIndex, fromJust $ LG.lab origGraph parentPrunedNodeIndex)
       (prunedNodes, prunedEdges) = LG.nodesAndEdgesAfter inGraph ([],[]) [parentPrunedNode]

       -- add back pruned component nodes and edges to post-order base component
       fullPostOrderGraph = LG.mkGraph ((LG.labNodes $ thd6 postOrderBaseGraph) ++ (parentPrunedNode : prunedNodes)) ((LG.labEdges $ thd6 postOrderBaseGraph) ++ prunedEdges) 

       -- this has block and character trees from postOrder of base graph and simple and cononical tree from fullPostOrderGraph
       fullPostOrderPhylogeneticGraph = (GO.convertDecoratedToSimpleGraph fullPostOrderGraph, prunedSubGraphCost + (snd6 postOrderBaseGraph) + localRootCost, fullPostOrderGraph, fth6 postOrderBaseGraph, fft6 postOrderBaseGraph, charInfoVV) 

       -- perform pre-order on base component 
       completeSplitGraph = if (swapType /= "tbr") then PRE.preOrderTreeTraversal inGS (finalAssignment inGS) (nonExactCharacters > 0) localStartVertex fullPostOrderPhylogeneticGraph
                            else -- TBR requires preorder for pruned component
                                 error "TBR not yet implemented"

       -- update and add back label for parentPrunedNode which is removed in the partial preorder pass
       canonicalSplitGraph = thd6 completeSplitGraph
       edgesFromParentPrunedNode = LG.out canonicalSplitGraph parentPrunedNodeIndex

       -- crete new cnonical graph deleting the unlabelled parentPrunedNode and adding labelled version and teh edge from it (should be one)
       canonicalSplitGraph' = LG.insEdges edgesFromParentPrunedNode $ LG.insNode parentPrunedNode $ LG.delNode parentPrunedNodeIndex canonicalSplitGraph


   in
   trace ("RGFV-After: \n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph (thd6 completeSplitGraph)) ++ " Pruned: " ++ (show prunedSubGraphRootVertex) 
      ++ " From: "  ++ (show $ fst parentPrunedNode) ++ "\n" 
      ++ (show $ fmap fst prunedNodes) ++ " " ++ (show $ fmap LG.toEdge prunedEdges) ) ( -- ++ "\n" 
   --    ++ (show $ (LG.labNodes canonicalSplitGraph') !! (fst parentPrunedNode)))

   -- check if base graph has fewer than 3 leaves (5 nodes) -- then nowhere to readd and screwes things up later
   if (length $ LG.nodes $ fst6 postOrderBaseGraph)  - (length prunedNodes ) < 5 then  trace ("Too few nodes") (LG.empty, infinity)
   else (canonicalSplitGraph', prunedSubGraphCost + (snd6 postOrderBaseGraph) + localRootCost)
   )
      )


-- | swapSteepest performs branch swapping greedily swithcing to found graph if better
swapSteepest   :: String 
               -> GlobalSettings 
               -> ProcessedData 
               -> Int 
               -> Bool 
               -> Int 
               -> VertexCost 
               -> [PhylogeneticGraph] 
               -> [PhylogeneticGraph] 
               -> Int 
               -> SimpleGraph 
               -> DecoratedGraph 
               -> DecoratedGraph 
               -> Bool
               -> V.Vector (V.Vector CharInfo) 
               -> Bool
               -> ([PhylogeneticGraph], Int)
swapSteepest swapType inGS inData numToKeep steepest counter curBestCost curBetterList inGraphList  numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA =
   if null inGraphList then (curBetterList, counter)
   else 
      (inGraphList, counter)

