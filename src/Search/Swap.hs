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
import qualified GraphOptimization.PostOrderFunctions as POS
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
   if null inGraphList then []
   else 
      trace ("Swapping " ++ (show $ length inGraphList) ++ " input graph(s) with minimum cost "++ (show $ minimum $ fmap snd6 inGraphList)) (
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
               charInfoVV = six6 $ head inGraphList

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
              trace ("After swap: " ++ (show $ length newGraphList'') ++ " resulting graphs with swap rounds (total): " ++ (show counterNNI) ++ " NNI, " ++ (show counterSPR) ++ " SPR, " ++ (show counterTBR) ++ " TBR")
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
   -- trace ("In swapSPRTBR:") (
   if LG.isEmpty (fst6 inGraph) then ([], 0)
   else 
      -- steepest takes immediate best
      if steepest then 
         let (swappedGraphs, counter) = swapSteepest swapType inGS inData numToKeep True 0 (snd6 inGraph) [] [inGraph] numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA

             -- swap "all" after steepest descent
             (swappedGraphs', counter') = swapAll swapType inGS inData numToKeep True counter (snd6 inGraph) [] swappedGraphs numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA
         in
         (swappedGraphs', counter')

      -- All does all swaps before taking best
      else  
         -- trace ("Going into SwapAll") (
         let (swappedGraphs, counter) = swapAll swapType inGS inData numToKeep False 0 (snd6 inGraph) [] [inGraph] numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA
         in 
         -- trace ("SSPRTBR: " ++ (show (length swappedGraphs, counter)))
         (swappedGraphs, counter)
         -- )
         -- )
      
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
   if null inGraphList then 
      (curSameBetterList, counter)
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
          splitTupleList = fmap (GO.splitGraphOnEdge firstDecoratedGraph) breakEdgeList `using` PU.myParListChunkRDS
          (splitGraphList, graphRootList, prunedGraphRootIndexList,  originalConnectionOfPruned) = L.unzip4 splitTupleList

          reoptimizedSplitGraphList = zipWith3 (reoptimizeGraphFromVertex inGS inData swapType doIA charInfoVV firstDecoratedGraph) splitGraphList graphRootList prunedGraphRootIndexList `using` PU.myParListChunkRDS

          -- create rejoins-- adds in break list so don't remake the initial graph
          -- didn't concatMap so can parallelize later
          -- this cost prob doesn't include the root/net penalty--so need to figure out
          swapPairList = concat $ L.zipWith5 (rejoinGraphKeepBest inGS swapType curBestCost numToKeep steepest doIA charInfoVV) reoptimizedSplitGraphList graphRootList prunedGraphRootIndexList originalConnectionOfPruned (fmap head $ fmap ((LG.parents $ thd6 firstGraph).fst3) breakEdgeList)

          -- keeps only "best" heuristic swap costs graphs
          minimumCandidateGraphCost = if (null swapPairList) then infinity
                                      else minimum $ fmap snd swapPairList
          candidateSwapGraphList = filter ((== minimumCandidateGraphCost). snd) swapPairList

          
          -- this should be incremental--full 2-pass for now
          reoptimizedSwapGraphList = fmap (T.multiTraverseFullyLabelGraph inGS inData False False Nothing) (fmap fst candidateSwapGraphList) `using` PU.myParListChunkRDS
                                     

          -- selects best graph list based on full optimization
          bestSwapGraphList = GO.selectPhylogeneticGraph [("best", (show numToKeep))] 0 ["best"] reoptimizedSwapGraphList

          bestSwapCost = if null swapPairList then infinity
                         else snd6 $ head bestSwapGraphList

      in
      -- trace ("Breakable Edges :" ++ (show $ fmap LG.toEdge breakEdgeList) ++ "\nIn graph:\n" ++ (LG.prettify $ fst6 firstGraph)) (
      -- trace ("(Est, [FP]): " ++ (show minimumCandidateGraphCost) ++ " " ++ (show $ fmap snd6 reoptimizedSwapGraphList)) (
      -- either no better or more of same cost graphs
      -- trace ("BSG: " ++ " simple " ++ (LG.prettify $ fst6 $ head bestSwapGraphList) ++ " Decorated " ++ (LG.prettify $ thd6 $ head bestSwapGraphList) ++ "\nCharinfo\n" ++ (show $ charType $ V.head $ V.head $ six6 $ head bestSwapGraphList)) (
      if bestSwapCost == curBestCost then 
         --this needs to be better--informed by zero-length edges
         let newCurSameBestList = if firstGraph `notElem` curSameBetterList then (firstGraph : curSameBetterList)
                                  else curSameBetterList
             graphsToSwap = ((tail inGraphList) ++ bestSwapGraphList) L.\\ newCurSameBestList               
         in
         --trace ("Same cost: " ++ (show bestSwapCost) ++ " with " ++ (show $ length $ (tail inGraphList) ++ graphsToSwap) ++ " more to swap and " ++ (show $ length newCurSameBestList) 
         --    ++ " graphs in 'best' list")
         swapAll swapType inGS inData numToKeep steepest (counter + 1) curBestCost newCurSameBestList graphsToSwap numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA

      -- better cost graphs
      else if (bestSwapCost < curBestCost) then 
         -- trace ("Better cost: " ++ (show bestSwapCost))
         swapAll swapType inGS inData numToKeep steepest (counter + 1) bestSwapCost [] ((tail inGraphList) ++ bestSwapGraphList)  numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA

      -- didn't find equal or better graphs
      else 
         -- trace ("Worse cost")
         let newCurSameBestList = if firstGraph `notElem` curSameBetterList then (firstGraph : curSameBetterList)
                                  else curSameBetterList
         in
         swapAll swapType inGS inData numToKeep steepest (counter + 1) curBestCost newCurSameBestList (tail inGraphList) numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA
      -- )
      -- )
      -- )

-- | swapSteepest performs branch swapping greedily switching to found graph if better
   -- infomrs evaluation--less parallelism
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
swapSteepest swapType inGS inData numToKeep steepest counter curBestCost curSameBetterList inGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA =
   if null inGraphList then 
      (curSameBetterList, counter)
   else 
      let firstGraph = head inGraphList
          firstDecoratedGraph = thd6 firstGraph
          (firstRootIndex, _) = head $ LG.getRoots firstDecoratedGraph

          -- filter out edges from root since no use--would just rejoin
          firstEdgeList = filter ((/= firstRootIndex) . fst3) $ LG.labEdges firstDecoratedGraph

          -- determine edges to break on--'bridge' edges only for network
          -- longest edges first first
          breakEdgeList = if (graphType inGS) == Tree then GO.sortEdgeListByLength firstEdgeList
                          else GO.sortEdgeListByLength $ GO.getEdgeSplitList firstDecoratedGraph

          -- create list of breaks
          splitTupleList = fmap (GO.splitGraphOnEdge firstDecoratedGraph) breakEdgeList  
          
          (splitGraphList, graphRootList, prunedGraphRootIndexList,  originalConnectionOfPruned) = L.unzip4 splitTupleList

          reoptimizedSplitGraphList = zipWith3 (reoptimizeGraphFromVertex inGS inData swapType doIA charInfoVV firstDecoratedGraph) splitGraphList graphRootList prunedGraphRootIndexList 

          -- create rejoins-- reoptimized fully in steepest returns PhylogheneticGraph 
          reoptimizedSwapGraphList = rejoinGraphKeepBestSteepest inGS inData swapType curBestCost numToKeep True doIA charInfoVV $ L.zip5 reoptimizedSplitGraphList graphRootList prunedGraphRootIndexList originalConnectionOfPruned (fmap head $ fmap ((LG.parents $ thd6 firstGraph).fst3) breakEdgeList)

          -- this should be incremental--full 2-pass for now
          -- reoptimizedSwapGraph = T.multiTraverseFullyLabelGraph inGS inData False False Nothing $ fst $ head swapPairList 
                                     

          bestSwapCost = if null reoptimizedSwapGraphList then infinity
                         else snd $ head reoptimizedSwapGraphList

      in
      -- trace ("Breakable Edges :" ++ (show $ fmap LG.toEdge breakEdgeList) ++ "\nIn graph:\n" ++ (LG.prettify $ fst6 firstGraph)) (
      
      -- either no better or more of same cost graphs
      -- trace ("BSG: " ++ " simple " ++ (LG.prettify $ fst6 $ head bestSwapGraphList) ++ " Decorated " ++ (LG.prettify $ thd6 $ head bestSwapGraphList) ++ "\nCharinfo\n" ++ (show $ charType $ V.head $ V.head $ six6 $ head bestSwapGraphList)) (
      if (bestSwapCost < curBestCost) then 
         --trace ("Steepest better")
         swapSteepest swapType inGS inData numToKeep steepest (counter + 1) bestSwapCost [] (fmap fst reoptimizedSwapGraphList) numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA

      -- didn't find equal or better graphs
      else (inGraphList, counter + 1)
      
      


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
                    -> LG.Node 
                    -> [(SimpleGraph, VertexCost)]
rejoinGraphKeepBest inGS swapType curBestCost numToKeep steepest doIA charInfoVV (splitGraph, splitCost) graphRoot prunedGraphRootIndex nakedNode originalSplitNode = 
   -- case where swap split retunred empty because too few nodes in remaining graph to add to
   if LG.isEmpty splitGraph then []
   else
      let outgroupEdges = LG.out splitGraph graphRoot
          (_, prunedSubTreeEdges) = LG.nodesAndEdgesAfter splitGraph [(nakedNode, fromJust $ LG.lab splitGraph nakedNode)]
          edgesToInvade = if steepest then (GO.sortEdgeListByDistance splitGraph [originalSplitNode] [originalSplitNode])  L.\\ (outgroupEdges ++ prunedSubTreeEdges)
                          else (LG.labEdges splitGraph) L.\\ (outgroupEdges ++ prunedSubTreeEdges)
          candidateEditList = concatMap (addSubGraph inGS swapType doIA splitGraph prunedGraphRootIndex splitCost nakedNode charInfoVV) edgesToInvade `using` PU.myParListChunkRDS


          minCandidateCost = if (not $ null candidateEditList) then minimum $ fmap fst3 candidateEditList   
                             else infinity
      in
      -- trace ("RGKB: " ++ (show $ fmap LG.toEdge edgesToInvade) ++ " " ++ (show curBestCost) ++ " v " ++ (show minCandidateCost)) (
      if minCandidateCost > curBestCost then []
      else 
         let bestEdits = filter ((== minCandidateCost). fst3) candidateEditList
             splitGraphSimple = GO.convertDecoratedToSimpleGraph splitGraph
             swapSimpleGraphList = fmap (applyGraphEdits splitGraphSimple) bestEdits
         in
         zip swapSimpleGraphList (L.replicate (length swapSimpleGraphList) minCandidateCost)
      -- )

-- | rejoinGraphKeepBestSteepest rejoins split trees on available edges (non-root, and not original split)
-- if steepest is False does not sort order of edges, other wise sorts in order of closeness to original edge
-- uses delta
-- NNI sorts edges on propinquity taking first 2 edges
-- TBR does the rerooting of pruned subtree
-- originalConnectionOfPruned is the "naked" node that was creted when teh graph was split and will 
-- be used for the rejoin node in the middle of th einvaded edge
rejoinGraphKeepBestSteepest :: GlobalSettings 
                             -> ProcessedData
                             -> String 
                             -> VertexCost 
                             -> Int 
                             -> Bool 
                             -> Bool 
                             -> V.Vector (V.Vector CharInfo) 
                             -> [((DecoratedGraph, VertexCost) , LG.Node , LG.Node , LG.Node , LG.Node)]
                             -> [(PhylogeneticGraph, VertexCost)]
rejoinGraphKeepBestSteepest inGS inData swapType curBestCost numToKeep steepest doIA charInfoVV splitInfoList = 
   if null splitInfoList then []
   else
      -- trace ("In rejoin steepes with split node " ++ (show $ fft5 $ head splitInfoList)) (
      let ((splitGraph, splitCost), graphRoot, prunedGraphRootIndex, nakedNode, originalSplitNode) = head splitInfoList
          outgroupEdges = LG.out splitGraph graphRoot
          (_, prunedSubTreeEdges) = LG.nodesAndEdgesAfter splitGraph [(nakedNode, fromJust $ LG.lab splitGraph nakedNode)]
          edgesToInvade = (GO.sortEdgeListByDistance splitGraph [originalSplitNode] [originalSplitNode]) L.\\ (outgroupEdges ++ prunedSubTreeEdges)
          -- edgesToInvade =  (LG.labEdges splitGraph) L.\\ (outgroupEdges ++ prunedSubTreeEdges)
          candidateEditList = addSubGraphSteepest inGS swapType doIA splitGraph prunedGraphRootIndex splitCost curBestCost nakedNode charInfoVV edgesToInvade 

          minCandidateCost = if (not $ null candidateEditList) then fst3 $ head candidateEditList   
                             else infinity
      in
      -- trace ("RGKB: " ++ (show $ fmap LG.toEdge edgesToInvade) ++ " " ++ (show curBestCost) ++ " v " ++ (show minCandidateCost)) (
      
      -- case where swap split retunred empty because too few nodes in remaining graph to add to
      if LG.isEmpty splitGraph then []
      else if minCandidateCost < curBestCost then 
         let splitGraphSimple = GO.convertDecoratedToSimpleGraph splitGraph
             swapSimpleGraph = applyGraphEdits splitGraphSimple (head candidateEditList)

             -- reoptimize to check here--should definately be incremental
             reoptimizedCandidateGraph = T.multiTraverseFullyLabelGraph inGS inData False False Nothing swapSimpleGraph
         in
         -- trace ("(Est, FP): " ++ (show (minCandidateCost, snd6 reoptimizedCandidateGraph))) (
         if (snd6 reoptimizedCandidateGraph < curBestCost) then [(reoptimizedCandidateGraph, snd6 reoptimizedCandidateGraph)]
         else rejoinGraphKeepBestSteepest inGS inData swapType curBestCost numToKeep steepest doIA charInfoVV (tail splitInfoList)
         -- )
      else rejoinGraphKeepBestSteepest inGS inData swapType curBestCost numToKeep steepest doIA charInfoVV (tail splitInfoList) 
      -- ))

-- | addSubGraphSteepest "adds" a subtree back into an edge calculating the cost of the graph via the delta of the add and costs of the two components
-- used in "steepest" descendt swapping
addSubGraphSteepest :: GlobalSettings 
                     -> String
                     -> Bool 
                     -> DecoratedGraph 
                     -> LG.Node 
                     -> VertexCost 
                     -> VertexCost 
                     -> LG.Node 
                     -> V.Vector (V.Vector CharInfo) 
                     -> [LG.LEdge EdgeInfo] 
                     -> [(VertexCost, [LG.LEdge Double], [LG.Edge])] 
addSubGraphSteepest inGS swapType doIA inGraph prunedGraphRootIndex splitCost curBestCost nakedNode charInfoVV targetEdgeList =  
   if null targetEdgeList then []
   else 
      let targetEdge@(eNode, vNode, targetlabel) = head targetEdgeList
          existingEdgeCost = minLength targetlabel
          edge0 = (nakedNode, vNode, 0.0)
          edge1 = (eNode, nakedNode, 0.0)
          -- edge2 = (nakedNode, prunedGraphRootIndex, 0.0)
          --newNode = (nakedNode, TL.pack ("HTU" ++ (show nakedNode)))

          -- for SPRT/NNI only need preliminary state of root-pruned node
          -- for TBR ther are multiple verticces created for each edge
          subGraphEdgeVertDataPairList = if swapType == "spr" ||swapType == "nni"  then [((-1, -1), vertData $ fromJust $ LG.lab inGraph prunedGraphRootIndex)]
                                         else error "Steepest TBR not yet implemented"

          (delta, rerootEdge) = getSubGraphDelta targetEdge doIA inGraph charInfoVV (head subGraphEdgeVertDataPairList)
          tbrEdgesAdd = []
          tbrEdgesDelete =[]
      in
      -- trace ("ASGR: " ++ (show (delta, splitCost, delta + splitCost))) (
      -- do not redo origal edge so retun infinite cost and dummy edits
      if (eNode == nakedNode) then addSubGraphSteepest inGS swapType doIA inGraph prunedGraphRootIndex splitCost curBestCost nakedNode charInfoVV (tail targetEdgeList)

      -- better heursitic cost
      else if (delta + splitCost < curBestCost) then [(delta + splitCost, [edge0, edge1] ++ tbrEdgesAdd, (eNode, vNode) : tbrEdgesDelete)] 

      -- not better heuristic cost
      else addSubGraphSteepest inGS swapType doIA inGraph prunedGraphRootIndex splitCost curBestCost nakedNode charInfoVV (tail targetEdgeList)
      -- )

-- | addSubTree "adds" a subtree back into an edge calculating the cost of the graph via the delta of the add and costs of the two components
addSubGraph :: GlobalSettings 
            -> String
            -> Bool 
            -> DecoratedGraph 
            -> LG.Node 
            -> VertexCost 
            -> LG.Node 
            -> V.Vector (V.Vector CharInfo) 
            -> LG.LEdge EdgeInfo 
            -> [(VertexCost, [LG.LEdge Double], [LG.Edge])] 
addSubGraph inGS swapType doIA inGraph prunedGraphRootIndex splitCost nakedNode charInfoVV targetEdge@(eNode, vNode, targetlabel) =  
   let existingEdgeCost = minLength targetlabel
       edge0 = (nakedNode, vNode, 0.0)
       edge1 = (eNode, nakedNode, 0.0)
       -- edge2 = (nakedNode, prunedGraphRootIndex, 0.0)
       --newNode = (nakedNode, TL.pack ("HTU" ++ (show nakedNode)))
       -- if subtree fewer than 3 leaves then can only do an SPR rearragement--no rerro0ts
       (nodesAfterPrunedRoot, _) = LG.nodesAndEdgesAfter inGraph [(prunedGraphRootIndex, fromJust $ LG.lab inGraph prunedGraphRootIndex)]
       onlySPR = length nodesAfterPrunedRoot < 3
       subGraphEdgeVertDataPairList = if swapType == "spr" || swapType == "nni" || onlySPR then [((-1, -1), vertData $ fromJust $ LG.lab inGraph prunedGraphRootIndex)]
                                      else getPrunedEdgeData (graphType inGS) doIA inGraph prunedGraphRootIndex charInfoVV

       -- get deltas and edges for TBR rerooting of pruned subgraph
       deltaEdgePairList = fmap (getSubGraphDelta targetEdge doIA inGraph charInfoVV) subGraphEdgeVertDataPairList

       delta = minimum $ fmap fst deltaEdgePairList

       (_, minDeltaEdgePairList) = unzip $ filter ((== delta) . fst) deltaEdgePairList


       -- get TBR edits if rerooting took place
       tbrEdgeEditList = if length subGraphEdgeVertDataPairList == 1 then []
                         else fmap (getTBREdgeEdits inGraph prunedGraphRootIndex) minDeltaEdgePairList

       (tbrEdgesToAddList, tbeEdgesToDeleteList) = unzip tbrEdgeEditList

       deltaCostList = replicate (length tbrEdgeEditList) (delta + splitCost)
       edgesToAddList = zipWith (++) (replicate (length tbrEdgeEditList) [edge0, edge1]) tbrEdgesToAddList
       edgesToDeleteList = zipWith (:) (replicate (length tbrEdgeEditList) (eNode, vNode)) tbeEdgesToDeleteList
 
   in
   
   -- do not redo origal edge so retun infinite cost and dummy edits
   
   if (eNode == nakedNode) then  
      -- trace ("ASG: break edge") 
      [(infinity, [], [])]
   else
      -- SPR or single reroot TBR case
      if length subGraphEdgeVertDataPairList == 1 then [(delta + splitCost, [edge0, edge1], [(eNode, vNode)])]
      
      -- TBR reroots 
      else zip3 deltaCostList edgesToAddList edgesToDeleteList
   

-- | getTBREdgeEdits takes and edge and returns the list of edita to pruned subgraph 
-- as a pair of edges to add and those to delete
getTBREdgeEdits :: DecoratedGraph -> Int -> LG.Edge -> ([LG.LEdge Double],[LG.Edge])
getTBREdgeEdits inGraph prunedGraphRootIndex rerootEdge =
   ([],[])


-- | getPrunedEdgeData takes fully optimized pruned data and returns edge list
-- and edge data for TBR additions
getPrunedEdgeData ::  GraphType -> Bool -> DecoratedGraph -> Int -> V.Vector (V.Vector CharInfo) -> [(LG.Edge, VertexBlockData)]
getPrunedEdgeData graphType doIA inGraph prunedGraphRootIndex charInfoVV =
   if LG.isEmpty inGraph then error "Empty graph in getPrunedEdgeData"
   else 
      let (_, edgeAfterList) = LG.nodesAndEdgesAfter inGraph [(prunedGraphRootIndex, fromJust $ LG.lab inGraph prunedGraphRootIndex)]

          -- this so not rerooted on network edge--screws things up
          nonNetWorkEdgeList = if graphType /= Tree then filter ((== False) . (LG.isNetworkLabEdge inGraph)) edgeAfterList
                               else edgeAfterList
         
          -- oringal pruned root
          prunedRootEdges = LG.out inGraph prunedGraphRootIndex

          -- new virtual edge of oringal root 2 edges
          virtualRootEdge = (snd3 $ head prunedRootEdges, snd3 $ last prunedRootEdges, thd3 $ head prunedRootEdges)

          -- edges availabel for testing
          edgeAfterList' = virtualRootEdge : (nonNetWorkEdgeList L.\\ prunedRootEdges)
          edgeDataList = fmap (makeEdgeData doIA inGraph charInfoVV) edgeAfterList' 
      in
      if length prunedRootEdges /= 2 then error ("Incorrect number of out edges (should be 2) in root of pruned graph: " ++ (show $ length prunedRootEdges))
      else
         trace ("Pruned reroots: " ++ (show $ fmap LG.toEdge edgeAfterList' ))  
         zip (fmap LG.toEdge edgeAfterList') edgeDataList

-- | makeEdgeData takes and edge and makes the VertData for the edge from the union of the two vertices
makeEdgeData :: Bool -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.LEdge b -> VertexBlockData
makeEdgeData doIA inGraph charInfoVV inEdge=
   let eNode = fst3 inEdge
       vNode = snd3 inEdge
       eNodeVertData = vertData $ fromJust $ LG.lab inGraph eNode
       vNodeVertData = vertData $ fromJust $ LG.lab inGraph eNode
   in
   M.createEdgeUnionOverBlocks (not doIA) eNodeVertData vNodeVertData charInfoVV []



-- | getSubGraphDelta calculated cost of adding a subgraph into and edge
-- for SPR use the preliminary of subGraph to final of e and v nodes
-- can use median fruntions for postorder if set final-> prelim or e and f
getSubGraphDelta :: LG.LEdge EdgeInfo -> Bool -> DecoratedGraph -> V.Vector (V.Vector CharInfo) ->  (LG.Edge, VertexBlockData) -> (VertexCost, LG.Edge)
getSubGraphDelta (eNode, vNode, targetlabel) doIA inGraph charInfoVV subGraphEdgeVertDataPair = 
   let existingEdgeCost = minLength targetlabel
       eNodeVertData = vertData $ fromJust $ LG.lab inGraph eNode
       vNodeVertData = vertData $ fromJust $ LG.lab inGraph vNode
       subGraphVertData = snd subGraphEdgeVertDataPair

       -- create edge union 'character' blockData
       -- based on final assignments but set to preliminary
       -- need to filter gaps if DO, not itIA
       edgeUnionVertData = M.createEdgeUnionOverBlocks (not doIA) eNodeVertData vNodeVertData charInfoVV []

       -- Use edge union data for delta to edge data
       costMethod = if doIA then ImpliedAlignment
                    else DirectOptimization

       subGraphEdgeUnionCost = if (not doIA) then V.sum $ fmap V.sum $ fmap (fmap snd) $ POS.createVertexDataOverBlocks subGraphVertData edgeUnionVertData charInfoVV []
                              else error "IA not yet implemented in getSubGraphDelta"

       -- subGraphEdgeUnionCost = sum $ fmap fst $ V.zipWith3 (PRE.getBlockCostPairsFinal costMethod) subGraphVertData edgeUnionVertData charInfoVV

       
       {-
       This estimate very close to subGraphEdgeUnionCost.  Seems to always underestmate where as subGraphEdgeUnionCost can sometimes overestimate
       but 3x as much work if use costEV, 2x if use existingEdgeCost
       
       dummyE = M.createEdgeUnionOverBlocks (not doIA) eNodeVertData eNodeVertData charInfoVV []
       dummyV = M.createEdgeUnionOverBlocks (not doIA) vNodeVertData vNodeVertData charInfoVV []
       dummySGV = M.createEdgeUnionOverBlocks (not doIA) (PRE.setFinalToPreliminaryStates subGraphVertData) (PRE.setFinalToPreliminaryStates subGraphVertData) charInfoVV []

       costNewE = V.sum $ fmap V.sum $ fmap (fmap snd) $ POS.createVertexDataOverBlocks dummyE subGraphVertData charInfoVV []
       costNewV = V.sum $ fmap V.sum $ fmap (fmap snd) $ POS.createVertexDataOverBlocks dummyV subGraphVertData charInfoVV []
       costEV = V.sum $ fmap V.sum $ fmap (fmap snd) $ POS.createVertexDataOverBlocks dummyE dummyV charInfoVV []

       subGraphEdgeUnionCost' = (costNewE + costNewV - existingEdgeCost) / 2.0
       -}
       

   in
    --trace ("GSD:" ++ (show ((costNewE, costNewV, costEV))) ++ " -> " ++ (show subGraphEdgeUnionCost') ++  " v " ++ (show subGraphEdgeUnionCost))
    -- trace ("Delta: " ++ (show subGraphEdgeUnionCost))
    --min subGraphEdgeUnionCost  subGraphEdgeUnionCost'

   (subGraphEdgeUnionCost, fst subGraphEdgeVertDataPair)


-- | reoptimizeGraphFromVertex fully labels the component graph that is connected to the specified vertex
-- retuning that graph with 2 optimized components and their cost
-- both components goo through multi-traversal optimizations
-- doIA option to only do IA optimization as opposed to full thing--should be enormously faster--but yet more approximate
-- creates final for both base graph and priunned component due to rerooting non-concordance of preorder and post order assignments
-- terminology bse graph is the component with the original root, pruned that which has been removed form the original
-- graph to be readded to edge set
-- The function
      -- 1) optimizes two components seprately fomr their "root"
      -- 2) takes nodes and edges for each and cretes new graph
      -- 3) returns graph and summed cost of two components
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
reoptimizeGraphFromVertex inGS inData swapType doIA charInfoVV origGraph inSplitGraph startVertex prunedSubGraphRootVertex =

       -- these required for full optimization
   let nonExactCharacters = U.getNumberNonExactCharacters (thd3 inData)
       leafGraph = LG.extractLeafGraph origGraph

       -- crete simple graph version of split for post order pass
       splitGraphSimple = GO.convertDecoratedToSimpleGraph inSplitGraph


       -- create optimized base graph
       (postOrderBaseGraph, localRootCost, _) = if not doIA then T.generalizedGraphPostOrderTraversal inGS nonExactCharacters inData leafGraph (Just startVertex) splitGraphSimple
                                                else 
                                                   -- Use IA assingment but ONLY reoptimize the IA states 
                                                   error "IA reoptimizeGraphFromVertex not yet implemented"
       
       fullBaseGraph = PRE.preOrderTreeTraversal inGS (finalAssignment inGS) (nonExactCharacters > 0) startVertex True postOrderBaseGraph

       -- create fully optimized pruned graph
       (postOrderPrunedGraph, _, _) = if not doIA then T.generalizedGraphPostOrderTraversal inGS nonExactCharacters inData leafGraph (Just prunedSubGraphRootVertex) splitGraphSimple
                                        else 
                                          -- Use IA assingment but ONLY reoptimize the IA states 
                                          error "IA reoptimizeGraphFromVertex not yet implemented"
       fullPrunedGraph = PRE.preOrderTreeTraversal inGS (finalAssignment inGS) (nonExactCharacters > 0) prunedSubGraphRootVertex True postOrderPrunedGraph


       -- get root node of base graph
       startBaseNode = (startVertex, fromJust $ LG.lab (thd6 fullBaseGraph) startVertex)

       -- get root node of pruned graph--parent since that is the full pruned piece (keeping that node for addition to base graph and edge creation)
       startPrunedNode = (prunedSubGraphRootVertex, fromJust $ LG.lab origGraph prunedSubGraphRootVertex)
       startPrunedParentNode = head $ LG.labParents origGraph prunedSubGraphRootVertex
       startPrunedParentEdge = (fst startPrunedParentNode, prunedSubGraphRootVertex, dummyEdge)

       
       -- get nodes and edges in base and pruned graph (both PhylogeneticGrapgs so thd6)
       (baseGraphNonRootNodes, baseGraphEdges) = LG.nodesAndEdgesAfter (thd6 fullBaseGraph) [startBaseNode]

       (prunedGraphNonRootNodes, prunedGraphEdges) = LG.nodesAndEdgesAfter (thd6 fullPrunedGraph) [startPrunedNode]

       -- make fully optimized graph from base and split components
       fullSplitGraph = LG.mkGraph ([startBaseNode, startPrunedNode, startPrunedParentNode] ++  baseGraphNonRootNodes ++ prunedGraphNonRootNodes) (startPrunedParentEdge : (baseGraphEdges ++ prunedGraphEdges))

       -- cost of split graph to be later combined with re-addition delta for heuristic graph cost
       splitGraphCost = (snd6 fullBaseGraph) + (snd6 fullPrunedGraph) + localRootCost

   in
   --trace ("Orig graph cost " ++ (show $ subGraphCost $ fromJust $ LG.lab origGraph startVertex) ++ " Base graph cost " ++ (show $ snd6 fullBaseGraph) ++ " pruned subgraph cost " ++ (show $ snd6 fullPrunedGraph) ++ " at node " ++ (show prunedSubGraphRootVertex) ++ " parent " ++ (show $ fst startPrunedParentNode)) -- ++ "\nSplit Graph\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph fullSplitGraph)) (
   --trace ("SG:" ++ (show fullSplitGraph))
   (fullSplitGraph, splitGraphCost)
   
-- | applyGraphEdits takes a  graphs and list of nodes and edges to add and delete and creates new graph
applyGraphEdits :: (Show a, Show b) => LG.Gr a b -> (VertexCost, [LG.LEdge b], [LG.Edge]) ->  LG.Gr a b
applyGraphEdits inGraph editStuff@(_, edgesToAdd, edgesToDelete) = 
   let editedGraph = LG.insEdges edgesToAdd $ LG.delEdges edgesToDelete inGraph
   in
   -- trace ("AGE: " ++ (show editStuff) ++ "\nIn graph:\n" ++ (LG.prettify inGraph) ++ "New Graph:\n" ++ (LG.prettify editedGraph)) 
   editedGraph
   



