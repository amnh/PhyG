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

-- | buildArgList is the list of valid build arguments
swapArgList :: [String]
swapArgList = ["spr","tbr", "keep", "steepest", "all", "nni", "ia"]


-- | swapMaster processes and spawns the swap functions
-- the 2 x maxMoveDist since distance either side to list 2* dist on sorted edges
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
             moveLimitList = filter (not . null) $ fmap snd $ filter ((/="keep").fst) lcArgList
             maxMoveEdgeDist  
              | length moveLimitList > 1 =
                errorWithoutStackTrace ("Multiple maximum edge distance number specifications in swap command--can have only one (e.g. spr:2): " ++ show inArgs)
              | null moveLimitList = Just ((maxBound :: Int) `div` 2) 
              | otherwise = readMaybe (head moveLimitList) :: Maybe Int
        in
        if isNothing keepNum then errorWithoutStackTrace ("Keep specification not an integer in swap: "  ++ show (head keepList))
        else if isNothing maxMoveEdgeDist then errorWithoutStackTrace ("Maximum edge move distance specification not an integer (e.g. spr:2): "  ++ show (snd $ head keepList))
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
                                               let graphPairList = fmap (swapSPRTBR "nni" inGS inData (fromJust keepNum) 2 doSteepest numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA) inGraphList `using` PU.myParListChunkRDS
                                                   (graphListList, counterList) = unzip graphPairList
                                               in (concat graphListList, sum counterList)
                                             else (inGraphList, 0)
               (newGraphList', counterSPR)  = if doSPR then 
                                               let graphPairList = fmap (swapSPRTBR "spr" inGS inData (fromJust keepNum) (2 * (fromJust maxMoveEdgeDist)) doSteepest numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA) newGraphList `using` PU.myParListChunkRDS
                                                   (graphListList, counterList) = unzip graphPairList
                                               in (concat graphListList, sum counterList)
                                             else (newGraphList, 0)

               (newGraphList'', counterTBR) = if doTBR then 
                                               let graphPairList =  fmap (swapSPRTBR "tbr" inGS inData (fromJust keepNum) (2 * (fromJust maxMoveEdgeDist)) doSteepest numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA) newGraphList' `using` PU.myParListChunkRDS
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
swapSPRTBR swapType inGS inData numToKeep maxMoveEdgeDist steepest numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA inGraph = 
   -- trace ("In swapSPRTBR:") (
   if LG.isEmpty (fst6 inGraph) then ([], 0)
   else 
      -- steepest takes immediate best--does not keep equall cost
      if steepest then 
         let (swappedGraphs, counter) = swapSteepest swapType inGS inData numToKeep maxMoveEdgeDist True 0 (snd6 inGraph) [] [inGraph] numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA

             -- swap "all" after steepest descent
             (swappedGraphs', counter') = swapAll swapType inGS inData numToKeep maxMoveEdgeDist True counter (snd6 inGraph) [] swappedGraphs numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA
         in
         (swappedGraphs', counter')

      -- All does all swaps before taking best
      else  
         -- trace ("Going into SwapAll") (
         let (swappedGraphs, counter) = swapAll swapType inGS inData numToKeep maxMoveEdgeDist False 0 (snd6 inGraph) [] [inGraph] numLeaves leafGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA
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
swapAll swapType inGS inData numToKeep maxMoveEdgeDist steepest counter curBestCost curSameBetterList inGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA =
   --trace ("ALL") (
   if null inGraphList then 
      (take numToKeep $ GO.getBVUniqPhylogeneticGraph True curSameBetterList, counter)
   else 
      let firstGraph = head inGraphList
          firstDecoratedGraph = thd6 firstGraph
          (firstRootIndex, _) = head $ LG.getRoots firstDecoratedGraph

          -- determine edges to break on--'bridge' edges only for network
          -- filter out edges from root since no use--would just rejoin
          breakEdgeList = if (graphType inGS) == Tree then filter ((/= firstRootIndex) . fst3) $ LG.labEdges firstDecoratedGraph
                          else filter ((/= firstRootIndex) . fst3) $ GO.getEdgeSplitList firstDecoratedGraph

          -- create list of breaks
          splitTupleList = fmap (GO.splitGraphOnEdge firstDecoratedGraph) breakEdgeList `using` PU.myParListChunkRDS
          (splitGraphList, graphRootList, prunedGraphRootIndexList,  originalConnectionOfPruned) = L.unzip4 splitTupleList

          reoptimizedSplitGraphList = zipWith3 (reoptimizeGraphFromVertex inGS inData swapType doIA charInfoVV firstGraph) splitGraphList graphRootList prunedGraphRootIndexList `using` PU.myParListChunkRDS

          -- create rejoins-- adds in break list so don't remake the initial graph
          -- didn't concatMap so can parallelize later
          -- this cost prob doesn't include the root/net penalty--so need to figure out
          swapPairList = concat $ L.zipWith5 (rejoinGraphKeepBest inGS swapType curBestCost numToKeep maxMoveEdgeDist steepest doIA charInfoVV) reoptimizedSplitGraphList graphRootList prunedGraphRootIndexList originalConnectionOfPruned (fmap head $ fmap ((LG.parents $ thd6 firstGraph).fst3) breakEdgeList)

          -- keeps better heuristic swap costs graphs based on current best as opposed to minimum heuristic costs
          -- minimumCandidateGraphCost = if (null swapPairList) then infinity
          --                             else minimum $ fmap snd swapPairList
          candidateSwapGraphList = filter ((<= curBestCost). snd) swapPairList

          
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
         --equality informed by zero-length edges
         let newCurSameBestList = GO.getBVUniqPhylogeneticGraph True (firstGraph : curSameBetterList)
                                  -- if firstGraph `notElem` curSameBetterList then (firstGraph : curSameBetterList)
                                  -- else curSameBetterList
             graphsToSwap = ((tail inGraphList) ++ bestSwapGraphList) L.\\ newCurSameBestList               
         in
         --trace ("Same cost: " ++ (show bestSwapCost) ++ " with " ++ (show $ length $ (tail inGraphList) ++ graphsToSwap) ++ " more to swap and " ++ (show $ length newCurSameBestList) 
         --    ++ " graphs in 'best' list")
         swapAll swapType inGS inData numToKeep maxMoveEdgeDist steepest (counter + 1) curBestCost newCurSameBestList graphsToSwap numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA

      -- better cost graphs
      else if (bestSwapCost < curBestCost) then 
         -- trace ("Better cost: " ++ (show bestSwapCost))
         swapAll swapType inGS inData numToKeep maxMoveEdgeDist steepest (counter + 1) bestSwapCost bestSwapGraphList (bestSwapGraphList ++ (tail inGraphList))  numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA

      -- didn't find equal or better graphs
      else 
         -- trace ("Worse cost")
         let newCurSameBestList = GO.getBVUniqPhylogeneticGraph True (firstGraph : curSameBetterList)
                                  -- if firstGraph `notElem` curSameBetterList then (firstGraph : curSameBetterList)
                                  -- else curSameBetterList
         in
         swapAll swapType inGS inData numToKeep maxMoveEdgeDist steepest (counter + 1) curBestCost newCurSameBestList (tail inGraphList) numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA
      -- )
      -- )
      -- )

-- | swapSteepest performs branch swapping greedily switching to found graph if better
   -- infomrs evaluation--less parallelism
swapSteepest   :: String 
               -> GlobalSettings 
               -> ProcessedData 
               -> Int 
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
swapSteepest swapType inGS inData numToKeep maxMoveEdgeDist steepest counter curBestCost curSameBetterList inGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA =
   --trace ("steepest") (
   if null inGraphList then 
      (take numToKeep $ GO.getBVUniqPhylogeneticGraph True curSameBetterList, counter)
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

          reoptimizedSplitGraphList = zipWith3 (reoptimizeGraphFromVertex inGS inData swapType doIA charInfoVV firstGraph) splitGraphList graphRootList prunedGraphRootIndexList 

          -- create rejoins-- reoptimized fully in steepest returns PhylogheneticGraph 
          reoptimizedSwapGraphList = rejoinGraphKeepBestSteepest inGS inData swapType curBestCost numToKeep maxMoveEdgeDist True doIA charInfoVV $ L.zip5 reoptimizedSplitGraphList graphRootList prunedGraphRootIndexList originalConnectionOfPruned (fmap head $ fmap ((LG.parents $ thd6 firstGraph).fst3) breakEdgeList)

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
         swapSteepest swapType inGS inData numToKeep maxMoveEdgeDist steepest (counter + 1) bestSwapCost (fmap fst reoptimizedSwapGraphList) (fmap fst reoptimizedSwapGraphList) numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired hasNonExactChars charInfoVV doIA

      -- didn't find equal or better graphs
      else (inGraphList, counter + 1)
      
      --)


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
rejoinGraphKeepBest inGS swapType curBestCost numToKeep maxMoveEdgeDist steepest doIA charInfoVV (splitGraph, splitCost) graphRoot prunedGraphRootIndex nakedNode originalSplitNode = 
   -- case where swap split retunred empty because too few nodes in remaining graph to add to
   if LG.isEmpty splitGraph then []
   else
      let outgroupEdges = LG.out splitGraph graphRoot
          (_, prunedSubTreeEdges) = LG.nodesAndEdgesAfter splitGraph [(nakedNode, fromJust $ LG.lab splitGraph nakedNode)]

          --only sort if limited egde rejoins
          edgesToInvade = if maxMoveEdgeDist == (maxBound :: Int) then (LG.labEdges splitGraph) L.\\ prunedSubTreeEdges -- L.\\ (outgroupEdges ++ prunedSubTreeEdges)
                          else take maxMoveEdgeDist $ (GO.sortEdgeListByDistance splitGraph [originalSplitNode] [originalSplitNode]) L.\\ prunedSubTreeEdges

          prunedGraphRootNode = (prunedGraphRootIndex, fromJust $ LG.lab splitGraph prunedGraphRootIndex)
          (nodesAfterPrunedRoot, edgesInPrunedSubGraph) = LG.nodesAndEdgesAfter splitGraph [prunedGraphRootNode]
          onlySPR = length nodesAfterPrunedRoot < 3
       

          candidateEditList = concatMap (addSubGraph inGS swapType doIA splitGraph prunedGraphRootNode splitCost nakedNode onlySPR edgesInPrunedSubGraph charInfoVV) edgesToInvade `using` PU.myParListChunkRDS


          minCandidateCost = if (not $ null candidateEditList) then minimum $ fmap fst3 candidateEditList   
                             else infinity
      in
      -- trace ("RGKB: " ++ (show $ fmap LG.toEdge edgesToInvade) ++ " " ++ (show curBestCost) ++ " v " ++ (show minCandidateCost)) (
      if minCandidateCost > curBestCost then []
      else 
         let bestEdits = filter ((<= curBestCost). fst3) candidateEditList -- not minimum cancidate cost--better if checkk all equal or better than curent best
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
                             -> Int 
                             -> Bool 
                             -> Bool 
                             -> V.Vector (V.Vector CharInfo) 
                             -> [((DecoratedGraph, VertexCost) , LG.Node , LG.Node , LG.Node , LG.Node)]
                             -> [(PhylogeneticGraph, VertexCost)]
rejoinGraphKeepBestSteepest inGS inData swapType curBestCost numToKeep maxMoveEdgeDist steepest doIA charInfoVV splitInfoList = 
   if null splitInfoList then []
   else
      -- trace ("In rejoin steepes with split node " ++ (show $ fft5 $ head splitInfoList)) (
      let ((splitGraph, splitCost), graphRoot, prunedGraphRootIndex, nakedNode, originalSplitNode) = head splitInfoList
          outgroupEdges = LG.out splitGraph graphRoot
          (_, prunedSubTreeEdges) = LG.nodesAndEdgesAfter splitGraph [(nakedNode, fromJust $ LG.lab splitGraph nakedNode)]
          edgesToInvade = take maxMoveEdgeDist $ (GO.sortEdgeListByDistance splitGraph [originalSplitNode] [originalSplitNode]) L.\\ prunedSubTreeEdges -- L.\\ (outgroupEdges ++ prunedSubTreeEdges)
          

          prunedGraphRootNode = (prunedGraphRootIndex, fromJust $ LG.lab splitGraph prunedGraphRootIndex)
          (nodesAfterPrunedRoot, edgesInPrunedSubGraph) = LG.nodesAndEdgesAfter splitGraph [prunedGraphRootNode]
          onlySPR = length nodesAfterPrunedRoot < 3
          
          candidateGraphList = addSubGraphSteepest inGS inData swapType doIA splitGraph prunedGraphRootNode splitCost curBestCost nakedNode onlySPR edgesInPrunedSubGraph charInfoVV edgesToInvade

      in
      -- trace ("RGKB: " ++ (show $ fmap LG.toEdge edgesToInvade) ++ " " ++ (show curBestCost) ++ " v " ++ (show minCandidateCost)) (
      
      -- case where swap split retunred empty because too few nodes in remaining graph to add to
      if LG.isEmpty splitGraph || null candidateGraphList then []
      else if (snd6 $ head candidateGraphList) < curBestCost then [(head candidateGraphList, snd6 $ head candidateGraphList)]
      else rejoinGraphKeepBestSteepest inGS inData swapType curBestCost numToKeep maxMoveEdgeDist steepest doIA charInfoVV (tail splitInfoList) 
      -- ))

-- | addSubGraphSteepest "adds" a subtree back into an edge calculating the cost of the graph via the delta of the add and costs of the two components
-- used in "steepest" descendt swapping
addSubGraphSteepest :: GlobalSettings 
                     -> ProcessedData
                     -> String
                     -> Bool 
                     -> DecoratedGraph 
                     -> LG.LNode VertexInfo
                     -> VertexCost 
                     -> VertexCost 
                     -> LG.Node 
                     -> Bool
                     -> [LG.LEdge EdgeInfo]
                     -> V.Vector (V.Vector CharInfo) 
                     -> [LG.LEdge EdgeInfo] 
                     -> [PhylogeneticGraph]
addSubGraphSteepest inGS inData swapType doIA inGraph prunedGraphRootNode splitCost curBestCost nakedNode onlySPR edgesInPrunedSubGraph charInfoVV targetEdgeList =  
   if null targetEdgeList then []
   else 
      let targetEdge@(eNode, vNode, targetlabel) = head targetEdgeList
          existingEdgeCost = minLength targetlabel
          edge0 = (nakedNode, vNode, 0.0)
          edge1 = (eNode, nakedNode, 0.0)
          targetEdgeData = makeEdgeData doIA inGraph charInfoVV targetEdge
          -- edge2 = (nakedNode, prunedGraphRootIndex, 0.0)
          --newNode = (nakedNode, TL.pack ("HTU" ++ (show nakedNode)))

          -- for SPRT/NNI only need preliminary state of root-pruned node
          -- for TBR ther are multiple verticces created for each edge
          subGraphEdgeVertDataTripleList = if swapType == "spr" || swapType == "nni" || onlySPR then [((-1, -1), vertData $ snd prunedGraphRootNode, ([],[]))]
                                           else getPrunedEdgeData (graphType inGS) doIA inGraph prunedGraphRootNode edgesInPrunedSubGraph charInfoVV

          -- this is SPR/NNI
          (delta, rerootEdge, (tbrEdgesAdd, tbrEdgesDelete)) = getSubGraphDelta targetEdgeData doIA inGraph charInfoVV (head subGraphEdgeVertDataTripleList)

          -- this if TBR --need to recurse through all reootings of pruned tree
          newTBRList = getSubGraphDeltaTBR inGS inData targetEdgeData [edge0, edge1] (eNode, vNode) doIA inGraph splitCost curBestCost charInfoVV subGraphEdgeVertDataTripleList
          
      in
      -- trace ("ASGR: " ++ (show (delta, splitCost, delta + splitCost))) (
      -- do not redo origal edge so retun infinite cost and dummy edits
      if (eNode == nakedNode) then addSubGraphSteepest inGS inData swapType doIA inGraph prunedGraphRootNode splitCost curBestCost nakedNode onlySPR edgesInPrunedSubGraph charInfoVV (tail targetEdgeList)

      --TBR case
      else if length subGraphEdgeVertDataTripleList > 1 then 
         if null newTBRList then addSubGraphSteepest inGS inData swapType doIA inGraph prunedGraphRootNode splitCost curBestCost nakedNode onlySPR edgesInPrunedSubGraph charInfoVV (tail targetEdgeList)

         -- must be better if non-empty list
         else newTBRList

      -- better heursitic cost
      -- reoptimize to check  cost
      else if (delta + splitCost <= curBestCost) then 
         let splitGraphSimple = GO.convertDecoratedToSimpleGraph inGraph
             swapSimpleGraph = applyGraphEdits splitGraphSimple (delta + splitCost, [edge0, edge1] ++ tbrEdgesAdd, (eNode, vNode) : tbrEdgesDelete)
             reoptimizedCandidateGraph = T.multiTraverseFullyLabelGraph inGS inData False False Nothing swapSimpleGraph
         in
         if (snd6 reoptimizedCandidateGraph < curBestCost) then [reoptimizedCandidateGraph]
         else addSubGraphSteepest inGS inData swapType doIA inGraph prunedGraphRootNode splitCost curBestCost nakedNode onlySPR edgesInPrunedSubGraph charInfoVV (tail targetEdgeList)

      -- not better heuristic cost
      else addSubGraphSteepest inGS inData swapType doIA inGraph prunedGraphRootNode splitCost curBestCost nakedNode onlySPR edgesInPrunedSubGraph charInfoVV (tail targetEdgeList)
      -- )

-- | getSubGraphDeltaTBR calculated cost of adding a subgraph into and edge
-- for SPR use the preliminary of subGraph to final of e and v nodes
-- can use median fruntions for postorder if set final-> prelim or e and f
getSubGraphDeltaTBR :: GlobalSettings 
                    -> ProcessedData
                    -> VertexBlockData 
                    -> [LG.LEdge Double]
                    -> LG.Edge
                    -> Bool 
                    -> DecoratedGraph 
                    -> VertexCost
                    -> VertexCost
                    -> V.Vector (V.Vector CharInfo) 
                    -> [(LG.Edge, VertexBlockData, ([LG.LEdge Double],[LG.Edge]))]
                    -> [PhylogeneticGraph]
getSubGraphDeltaTBR inGS inData evEdgeData edgeToAddInList edgeToDeleteIn doIA inGraph splitCost curBestCost charInfoVV subGraphEdgeVertDataTripleList = 
   -- found nothing better
   if null subGraphEdgeVertDataTripleList then []
   else 
      let (edgeToJoin, subGraphVertData, (tbrEdgesAdd, tbrEdgesDelete)) = head subGraphEdgeVertDataTripleList
         -- Use edge union data for delta to edge data
          costMethod = if doIA then ImpliedAlignment
                       else DirectOptimization

          subGraphEdgeUnionCost = if (not doIA) then V.sum $ fmap V.sum $ fmap (fmap snd) $ POS.createVertexDataOverBlocks subGraphVertData evEdgeData charInfoVV []
                                  else V.sum $ fmap V.sum $ fmap (fmap snd) $ POS.createVertexDataOverBlocksStaticIA subGraphVertData evEdgeData charInfoVV []


      in
      if subGraphEdgeUnionCost + splitCost < curBestCost then 
         -- reoptimize and check
         let splitGraphSimple = GO.convertDecoratedToSimpleGraph inGraph
             swapSimpleGraph = applyGraphEdits splitGraphSimple (subGraphEdgeUnionCost + splitCost, edgeToAddInList ++ tbrEdgesAdd, edgeToDeleteIn : tbrEdgesDelete)
             reoptimizedCandidateGraph = T.multiTraverseFullyLabelGraph inGS inData False False Nothing swapSimpleGraph
         in
         if snd6 reoptimizedCandidateGraph < curBestCost then [reoptimizedCandidateGraph]
         else getSubGraphDeltaTBR inGS inData evEdgeData edgeToAddInList edgeToDeleteIn doIA inGraph splitCost curBestCost charInfoVV (tail subGraphEdgeVertDataTripleList)

      else
         -- move on 
         getSubGraphDeltaTBR inGS inData evEdgeData edgeToAddInList edgeToDeleteIn doIA inGraph splitCost curBestCost charInfoVV (tail subGraphEdgeVertDataTripleList)


-- | addSubTree "adds" a subtree back into an edge calculating the cost of the graph via the delta of the add and costs of the two components
-- does NOT reoptimize candidate trees--happens after return since tlooking at "all"
addSubGraph :: GlobalSettings 
            -> String
            -> Bool 
            -> DecoratedGraph 
            -> LG.LNode VertexInfo
            -> VertexCost 
            -> LG.Node 
            -> Bool
            -> [LG.LEdge EdgeInfo]
            -> V.Vector (V.Vector CharInfo) 
            -> LG.LEdge EdgeInfo 
            -> [(VertexCost, [LG.LEdge Double], [LG.Edge])] 
addSubGraph inGS swapType doIA inGraph prunedGraphRootNode splitCost parentPrunedRoot onlySPR edgesInPrunedSubGraph charInfoVV targetEdge@(eNode, vNode, targetlabel) =  
   let existingEdgeCost = minLength targetlabel
       edge0 = (parentPrunedRoot, vNode, 0.0)
       edge1 = (eNode, parentPrunedRoot, 0.0)
       targetEdgeData = makeEdgeData doIA inGraph charInfoVV targetEdge 
       -- edge2 = (nakedNode, prunedGraphRootIndex, 0.0)
       --newNode = (nakedNode, TL.pack ("HTU" ++ (show nakedNode)))
       -- if subtree fewer than 3 leaves then can only do an SPR rearragement--no rerro0ts
       -- prunedGraphRootNode = (prunedGraphRootIndex, fromJust $ LG.lab inGraph prunedGraphRootIndex)
       -- (nodesAfterPrunedRoot, edgesInPrunedSubGraph) = LG.nodesAndEdgesAfter inGraph [prunedGraphRootNode]
       -- onlySPR = length nodesAfterPrunedRoot < 3
       subGraphEdgeVertDataTripleList = if swapType == "spr" || swapType == "nni" || onlySPR then [((-1, -1), vertData $ snd prunedGraphRootNode, ([],[]))]
                                      else getPrunedEdgeData (graphType inGS) doIA inGraph prunedGraphRootNode edgesInPrunedSubGraph charInfoVV

       -- get deltas and edges for TBR rerooting of pruned subgraph
       deltaEdgeTripleList = fmap (getSubGraphDelta targetEdgeData doIA inGraph charInfoVV) subGraphEdgeVertDataTripleList `using` PU.myParListChunkRDS

       delta = minimum $ fmap fst3 deltaEdgeTripleList

       minDeltaEditList = fmap thd3 $ filter ((== delta) . fst3) deltaEdgeTripleList


       -- get TBR edits if rerooting took place
       tbrEdgeEditList = if length subGraphEdgeVertDataTripleList == 1 then []
                         else minDeltaEditList

       (tbrEdgesToAddList, tbeEdgesToDeleteList) = unzip tbrEdgeEditList

       deltaCostList = replicate (length tbrEdgeEditList) (delta + splitCost)
       edgesToAddList = zipWith (++) (replicate (length tbrEdgeEditList) [edge0, edge1]) tbrEdgesToAddList
       edgesToDeleteList = zipWith (:) (replicate (length tbrEdgeEditList) (eNode, vNode)) tbeEdgesToDeleteList
 
   in
   
   -- do not redo origal edge so retun infinite cost and dummy edits
   
   if (eNode == parentPrunedRoot) then  
      -- trace ("ASG: break edge") 
      [(infinity, [], [])]
   else
      -- SPR or single reroot TBR case
      if length subGraphEdgeVertDataTripleList == 1 then [(delta + splitCost, [edge0, edge1], [(eNode, vNode)])]
      
      -- TBR reroots 
      else zip3 deltaCostList edgesToAddList edgesToDeleteList
   

-- | getTBREdgeEdits takes and edge and returns the list of edita to pruned subgraph 
-- as a pair of edges to add and those to delete
-- since reroot edge is directed (e,v), edges away from v will have correct
-- orientation. Edges between 'e' and the root will have to be flipped
-- original root edges and rerort edge are deleted and new root and edge spanning orginal root created
-- returns ([add], [delete])
getTBREdgeEdits :: DecoratedGraph -> LG.LNode VertexInfo -> [LG.LEdge EdgeInfo] -> LG.Edge -> ([LG.LEdge Double],[LG.Edge])
getTBREdgeEdits inGraph prunedGraphRootNode edgesInPrunedSubGraph rerootEdge =
   --trace ("Gettiung TBR Edits for " ++ (show rerootEdge)) (
   let prunedGraphRootIndex = fst prunedGraphRootNode
       originalRootEdgeNodes = LG.descendants inGraph prunedGraphRootIndex
       originalRootEdges = LG.out inGraph prunedGraphRootIndex
       
       -- get path from new root edge fst vertex to orginal root and flip those edges
       closerToPrunedRootEdgeNode = (fst rerootEdge, fromJust $ LG.lab inGraph $ fst rerootEdge)
       (nodesInPath, edgesinPath) = LG.postOrderPathToNode inGraph closerToPrunedRootEdgeNode prunedGraphRootNode 
       
       -- don't want original root edges to be flipped since deleted
       edgesToFlip = edgesinPath L.\\ originalRootEdges
       flippedEdges = fmap GO.convertToSimpleEdge $ fmap LG.flipLEdge edgesToFlip

       -- new edges on new root position and spanning old root
       -- ad in closer vertex to root to make sure direction of edge is correct
       newEdgeOnOldRoot = if (snd3 $ head originalRootEdges) `elem` ((fst rerootEdge) : (fmap fst nodesInPath)) then (snd3 $ head originalRootEdges, snd3 $ last originalRootEdges, 0.0)
                          else (snd3 $ last originalRootEdges, snd3 $ head originalRootEdges, 0.0)
       newRootEdges = [(prunedGraphRootIndex, fst rerootEdge, 0.0 ),(prunedGraphRootIndex, snd rerootEdge, 0.0)]


   in
   -- original root edge so no change
   if (fst rerootEdge) `elem` originalRootEdgeNodes &&  (snd rerootEdge) `elem` originalRootEdgeNodes then ([],[])

   -- rerooted
   else 
      -- delete orignal root edges and rerootEdge
      -- add new root edges 
      -- and new edge on old root--but need orientation
      -- flip edges from new root to old (delete and add list) 
      --trace ("\n\nIn Graph:\n"++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inGraph) ++ "\nTBR Edits: " ++ (show (rerootEdge, prunedGraphRootIndex, fmap LG.toEdge flippedEdges))
      --   ++ "\nEdges to add: " ++ (show $ fmap LG.toEdge $ newEdgeOnOldRoot : (flippedEdges ++ newRootEdges)) ++ "\nEdges to delete: " ++ (show $ rerootEdge : (fmap LG.toEdge (edgesToFlip ++ originalRootEdges))))
      (newEdgeOnOldRoot : (flippedEdges ++ newRootEdges), rerootEdge : (fmap LG.toEdge (edgesToFlip ++ originalRootEdges)))
      --)


-- | getPrunedEdgeData takes fully optimized pruned data and returns edge list, edge data for TBR additions,
-- and graph edits for that edge reroot
getPrunedEdgeData :: GraphType 
                  -> Bool 
                  -> DecoratedGraph 
                  -> LG.LNode VertexInfo 
                  -> [LG.LEdge EdgeInfo] 
                  -> V.Vector (V.Vector CharInfo) 
                  -> [(LG.Edge, VertexBlockData, ([LG.LEdge Double],[LG.Edge]))]
getPrunedEdgeData graphType doIA inGraph prunedGraphRootNode edgesInPrunedSubGraph charInfoVV =
   if LG.isEmpty inGraph then error "Empty graph in getPrunedEdgeData"
   else 
      let -- (_, edgeAfterList) = LG.nodesAndEdgesAfter inGraph [prunedGraphRootNode]

          -- this so not rerooted on network edge--screws things up
          nonNetWorkEdgeList = if graphType /= Tree then filter ((== False) . (LG.isNetworkLabEdge inGraph)) edgesInPrunedSubGraph
                               else edgesInPrunedSubGraph
         
          -- oringal pruned root
          prunedRootEdges = LG.out inGraph $ fst prunedGraphRootNode

          -- new virtual edge of oringal root 2 edges
          virtualRootEdge = (snd3 $ head prunedRootEdges, snd3 $ last prunedRootEdges, thd3 $ head prunedRootEdges)

          -- edges availabel for testing
          edgeAfterList' = virtualRootEdge : (nonNetWorkEdgeList L.\\ prunedRootEdges)

          -- could be parallelized
          edgeDataList = fmap (makeEdgeData doIA inGraph charInfoVV) edgeAfterList' `using` PU.myParListChunkRDS

          -- get potential TBR edits--here so not recalculated multiple times for each prune
          edgeEdits = fmap (getTBREdgeEdits inGraph prunedGraphRootNode edgesInPrunedSubGraph) (fmap LG.toEdge edgeAfterList') `using` PU.myParListChunkRDS
      in
      if length prunedRootEdges /= 2 then error ("Incorrect number of out edges (should be 2) in root of pruned graph: " ++ (show $ length prunedRootEdges))
      else
         --trace ("PED: " ++ (show  $ zip (fmap length edgeDataList) (fmap LG.toEdge edgeAfterList')) ++ " SG edges: " ++ (show $ fmap LG.toEdge edgesInPrunedSubGraph))  
         zip3 (fmap LG.toEdge edgeAfterList') edgeDataList edgeEdits

-- | makeEdgeData takes and edge and makes the VertData for the edge from the union of the two vertices
makeEdgeData :: Bool -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.LEdge b -> VertexBlockData
makeEdgeData doIA inGraph charInfoVV inEdge =
   let eNode = fst3 inEdge
       vNode = snd3 inEdge
       eNodeVertData = vertData $ fromJust $ LG.lab inGraph eNode
       vNodeVertData = vertData $ fromJust $ LG.lab inGraph eNode
   in
   M.createEdgeUnionOverBlocks doIA (not doIA) eNodeVertData vNodeVertData charInfoVV []



-- | getSubGraphDelta calculated cost of adding a subgraph into and edge
-- for SPR use the preliminary of subGraph to final of e and v nodes
-- can use median fruntions for postorder if set final-> prelim or e and f
getSubGraphDelta :: VertexBlockData 
                 -> Bool 
                 -> DecoratedGraph 
                 -> V.Vector (V.Vector CharInfo) 
                 -> (LG.Edge, VertexBlockData, ([LG.LEdge Double],[LG.Edge])) 
                 -> (VertexCost, LG.Edge, ([LG.LEdge Double],[LG.Edge]))
getSubGraphDelta evEdgeData doIA inGraph charInfoVV (edgeToJoin, subGraphVertData, edgeSubGraphEdits) = 
   let --existingEdgeCost = minLength targetlabel
       --eNodeVertData = vertData $ fromJust $ LG.lab inGraph eNode
       --NodeVertData = vertData $ fromJust $ LG.lab inGraph vNode
       -- subGraphVertData = snd subGraphEdgeVertDataPair

       -- create edge union 'character' blockData
       -- based on final assignments but set to preliminary
       -- need to filter gaps if DO, not itIA
       --edgeUnionVertData = makeEdgeData doIA inGraph charInfoVV evEdge
       -- edgeUnionVertData = M.createEdgeUnionOverBlocks (not doIA) eNodeVertData vNodeVertData charInfoVV []

       -- Use edge union data for delta to edge data
       costMethod = if doIA then ImpliedAlignment
                    else DirectOptimization

       subGraphEdgeUnionCost = if (not doIA) then V.sum $ fmap V.sum $ fmap (fmap snd) $ POS.createVertexDataOverBlocks subGraphVertData evEdgeData charInfoVV []
                              else V.sum $ fmap V.sum $ fmap (fmap snd) $ POS.createVertexDataOverBlocksStaticIA subGraphVertData evEdgeData charInfoVV []

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
   -- remove this check when things are working
   if null subGraphVertData || null evEdgeData || (subGraphEdgeUnionCost == 0.0) then 
      trace ("SGD null or 0 :" ++ (show (edgeToJoin, length subGraphVertData, length evEdgeData, subGraphEdgeUnionCost)) )
      --  ++ "\nInGraph:\n"  
      --  ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inGraph) ++ "\n" ++ (show inGraph))
      (subGraphEdgeUnionCost, edgeToJoin, edgeSubGraphEdits) 
    --trace ("GSD:" ++ (show ((costNewE, costNewV, costEV))) ++ " -> " ++ (show subGraphEdgeUnionCost') ++  " v " ++ (show subGraphEdgeUnionCost))
    -- trace ("Delta: " ++ (show subGraphEdgeUnionCost))
    --min subGraphEdgeUnionCost  subGraphEdgeUnionCost'
   else 
      -- trace ("SGD no 0:" ++ (show (subGraphEdgeUnionCost,edgeToJoin)))
      (subGraphEdgeUnionCost, edgeToJoin, edgeSubGraphEdits)


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
-- if doIA is TRUE then call function that onl;y optimizes the IA assignments on the "original graph" after split.
-- this keeps teh IA chracters in sunc across the two graphs
reoptimizeGraphFromVertex :: GlobalSettings 
                          -> ProcessedData 
                          -> String 
                          -> Bool 
                          -> V.Vector (V.Vector CharInfo) 
                          -> PhylogeneticGraph 
                          -> DecoratedGraph 
                          -> Int 
                          -> Int 
                          -> (DecoratedGraph, VertexCost)
reoptimizeGraphFromVertex inGS inData swapType doIA charInfoVV origPhyloGraph inSplitGraph startVertex prunedSubGraphRootVertex =
   if doIA then 
      -- only reoptimize the IA states for dynamic characters
      reoptimizeGraphFromVertexIA inGS inData swapType charInfoVV origPhyloGraph inSplitGraph startVertex prunedSubGraphRootVertex 
   else
      -- perform full optimizations of nodes
      -- these required for full optimization
      let nonExactCharacters = U.getNumberNonExactCharacters (thd3 inData)
          origGraph = thd6 origPhyloGraph
          leafGraph = LG.extractLeafGraph origGraph
          calcBranchLengths = False

          -- create simple graph version of split for post order pass
          splitGraphSimple = GO.convertDecoratedToSimpleGraph inSplitGraph


          -- create optimized base graph
          -- False for staticIA
          (postOrderBaseGraph, localRootCost, _) = T.generalizedGraphPostOrderTraversal inGS nonExactCharacters inData leafGraph False (Just startVertex) splitGraphSimple
                                                   
          
          fullBaseGraph = PRE.preOrderTreeTraversal inGS (finalAssignment inGS) False calcBranchLengths (nonExactCharacters > 0) startVertex True postOrderBaseGraph

          -- create fully optimized pruned graph.  Post order tehn preorder

          -- get root node of pruned graph--parent since that is the full pruned piece (keeping that node for addition to base graph and edge creation)
          startPrunedNode = (prunedSubGraphRootVertex, fromJust $ LG.lab origGraph prunedSubGraphRootVertex)
          startPrunedParentNode = head $ LG.labParents origGraph prunedSubGraphRootVertex
          startPrunedParentEdge = (fst startPrunedParentNode, prunedSubGraphRootVertex, dummyEdge)


          -- False for staticIA
          (postOrderPrunedGraph, _, _) = T.generalizedGraphPostOrderTraversal inGS nonExactCharacters inData leafGraph False (Just prunedSubGraphRootVertex) splitGraphSimple
                                                

          -- False for staticIA
          fullPrunedGraph = PRE.preOrderTreeTraversal inGS (finalAssignment inGS) False calcBranchLengths (nonExactCharacters > 0) prunedSubGraphRootVertex True postOrderPrunedGraph
         
          -- get root node of base graph
          startBaseNode = (startVertex, fromJust $ LG.lab (thd6 fullBaseGraph) startVertex)

          
          
          -- get nodes and edges in base and pruned graph (both PhylogeneticGrapgs so thd6)
          (baseGraphNonRootNodes, baseGraphEdges) = LG.nodesAndEdgesAfter (thd6 fullBaseGraph) [startBaseNode]

          (prunedGraphNonRootNodes, prunedGraphEdges) = if LG.isLeaf (thd6 origPhyloGraph) prunedSubGraphRootVertex then ([],[])
                                                        else LG.nodesAndEdgesAfter (thd6 fullPrunedGraph) [startPrunedNode]

          -- make fully optimized graph from base and split components
          fullSplitGraph = LG.mkGraph ([startBaseNode, startPrunedNode, startPrunedParentNode] ++  baseGraphNonRootNodes ++ prunedGraphNonRootNodes) (startPrunedParentEdge : (baseGraphEdges ++ prunedGraphEdges))

          -- cost of split graph to be later combined with re-addition delta for heuristic graph cost
          prunedCost = if LG.isLeaf (thd6 origPhyloGraph) prunedSubGraphRootVertex then 0
                       else snd6 fullPrunedGraph
          splitGraphCost = (snd6 fullBaseGraph) + prunedCost + localRootCost

      in
      
      -- trace ("Orig graph cost " ++ (show $ subGraphCost $ fromJust $ LG.lab origGraph startVertex) ++ " Base graph cost " ++ (show $ snd6 fullBaseGraph) ++ " pruned subgraph cost " ++ (show prunedCost) ++ " at node " ++ (show prunedSubGraphRootVertex) ++ " parent " ++ (show $ fst startPrunedParentNode))
         {-
         ++ "\nBaseGraphNodes\n" ++ (show $ L.sort  $ fmap fst baseGraphNonRootNodes) ++ "\nPruned nodes from root: " ++ "\n" ++ (show $ fmap fst $ startPrunedNode : prunedGraphNonRootNodes) 
         ++ "\nSplit Graph\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph fullSplitGraph)
         ++ "\nOrig graph:\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph origGraph))
         -}
      (fullSplitGraph, splitGraphCost)

      
-- | reoptimizeGraphFromVertexIA performs operations of reoptimizeGraphFromVertex for static charcaters
-- but dynamic characters--only update IA assignments and initialized from origPhylo graph (at leaves) to keep IA characters in sync
-- since all "static" only need single traversal post order pass
reoptimizeGraphFromVertexIA :: GlobalSettings 
                          -> ProcessedData 
                          -> String 
                          -> V.Vector (V.Vector CharInfo) 
                          -> PhylogeneticGraph 
                          -> DecoratedGraph 
                          -> Int 
                          -> Int 
                          -> (DecoratedGraph, VertexCost)
reoptimizeGraphFromVertexIA inGS inData swapType charInfoVV origPhyloGraph inSplitGraph startVertex prunedSubGraphRootVertex =
   if graphType inGS /= Tree then error "Networks not yet implemented in reoptimizeGraphFromVertexIA"
   else 
      let   nonExactCharacters = U.getNumberNonExactCharacters (thd3 inData)
            origGraph = thd6 origPhyloGraph

            -- create leaf graphs--but copy IA final to prelim
            leafGraph = GO.copyIAFinalToPrelim $ LG.extractLeafGraph origGraph
            calcBranchLengths = False

            -- create simple graph version of split for post order pass
            splitGraphSimple = GO.convertDecoratedToSimpleGraph inSplitGraph

            --Create base graph
            -- create postorder assignment--but only from single traversal
            -- True flag fior staticIA
            postOrderBaseGraph = T.postOrderTreeTraversal inGS inData leafGraph True (Just startVertex) splitGraphSimple
            baseGraphCost = snd6 postOrderBaseGraph
            
            -- True flag fior staticIA
            fullBaseGraph = PRE.preOrderTreeTraversal inGS (finalAssignment inGS) True calcBranchLengths (nonExactCharacters > 0) startVertex True postOrderBaseGraph

            localRootCost = if (rootCost inGS) == NoRootCost then 0.0
                              else if (rootCost inGS) == Wheeler2015Root then T.getW15RootCost inData postOrderBaseGraph
                              else error ("Root cost type " ++ (show $ rootCost inGS) ++ " is not yet implemented")

            -- get root node of base graph
            startBaseNode = (startVertex, fromJust $ LG.lab (thd6 fullBaseGraph) startVertex)    

            --Create pruned graph
            -- get root node of pruned graph--parent since that is the full pruned piece (keeping that node for addition to base graph and edge creation)
            startPrunedNode = GO.makeIAPrelimFromFinal (prunedSubGraphRootVertex, fromJust $ LG.lab origGraph prunedSubGraphRootVertex)
            startPrunedParentNode =  head $ LG.labParents origGraph prunedSubGraphRootVertex
            startPrunedParentEdge = (fst startPrunedParentNode, prunedSubGraphRootVertex, dummyEdge)


            -- True flag fior staticIA
            postOrderPrunedGraph =  T.postOrderTreeTraversal inGS inData leafGraph True (Just prunedSubGraphRootVertex) splitGraphSimple
            prunedGraphCost = snd6 postOrderPrunedGraph                             

            -- True flag fior staticIA
            fullPrunedGraph = PRE.preOrderTreeTraversal inGS (finalAssignment inGS) True calcBranchLengths (nonExactCharacters > 0) prunedSubGraphRootVertex True postOrderPrunedGraph

            -- get nodes and edges in base and pruned graph (both PhylogeneticGrapgs so thd6)
            (baseGraphNonRootNodes, baseGraphEdges) = LG.nodesAndEdgesAfter (thd6 fullBaseGraph) [startBaseNode]

            (prunedGraphNonRootNodes, prunedGraphEdges) = if LG.isLeaf (thd6 origPhyloGraph) prunedSubGraphRootVertex then ([],[])
                                                          else LG.nodesAndEdgesAfter (thd6 fullPrunedGraph) [startPrunedNode]

            -- make fully optimized graph from base and split components
            fullSplitGraph = LG.mkGraph ([startBaseNode, startPrunedNode, startPrunedParentNode] ++  baseGraphNonRootNodes ++ prunedGraphNonRootNodes) (startPrunedParentEdge : (baseGraphEdges ++ prunedGraphEdges))

            splitGraphCost = baseGraphCost + prunedGraphCost + localRootCost

      in
      -- remove when working
      -- trace ("ROGFVIA split costs:" ++ (show (baseGraphCost, prunedGraphCost, localRootCost)) ++ " -> " ++ (show splitGraphCost)) (
      if splitGraphCost == 0 then
         error ("Split costs:" ++ (show (baseGraphCost, prunedGraphCost, localRootCost)) ++ " -> " ++ (show splitGraphCost)
            ++ " Split graph simple:\n" ++ (LG.prettify splitGraphSimple) 
            ++ "\nFull:\n" ++ (show inSplitGraph)
            ++ "\nOriginal Graph:\n" ++ (show origGraph))
      else (fullSplitGraph, splitGraphCost)
      -- )

-- | applyGraphEdits takes a  graphs and list of nodes and edges to add and delete and creates new graph
applyGraphEdits :: (Show a, Show b) => LG.Gr a b -> (VertexCost, [LG.LEdge b], [LG.Edge]) ->  LG.Gr a b
applyGraphEdits inGraph editStuff@(_, edgesToAdd, edgesToDelete) = 
   let editedGraph = LG.insEdges edgesToAdd $ LG.delEdges edgesToDelete inGraph
   in
   -- trace ("AGE: " ++ (show editStuff) ++ "\nIn graph:\n" ++ (LG.prettify inGraph) ++ "New Graph:\n" ++ (LG.prettify editedGraph)) 
   editedGraph
   

