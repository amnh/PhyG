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

module Search.Swap  ( swapSPRTBR
                    , reoptimizeSplitGraphFromVertexTuple
                    , rejoinGraphTuple
                    ) where

import           Control.Parallel.Strategies
import qualified Data.List                            as L
import           Data.Maybe
import qualified Data.Vector                          as V
import           GeneralUtilities
import qualified GraphOptimization.Medians            as M
import qualified GraphOptimization.PostOrderFunctions as POS
import qualified GraphOptimization.PreOrderFunctions  as PRE
import qualified GraphOptimization.Traversals         as T
import qualified Graphs.GraphOperations               as GO
import qualified ParallelUtilities                    as PU
import           Types.Types
import qualified Utilities.LocalGraph                 as LG
import           Utilities.Utilities                  as U
import qualified GraphOptimization.PostOrderSoftWiredFunctions as POSW
-- import           Debug.Trace


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
            -> Bool
            -> Bool
            -> Bool
            -> (Maybe SAParams, PhylogeneticGraph)
            -> ([PhylogeneticGraph], Int)
swapSPRTBR swapType inGS inData numToKeep maxMoveEdgeDist steepest alternate doIA returnMutated (inSimAnnealParams, inGraph) =
   -- trace ("In swapSPRTBR:") (
   if LG.isEmpty (fst6 inGraph) then ([], 0)
   else
      let numLeaves = V.length $ fst3 inData
          leafGraph = T.makeSimpleLeafGraph inData
          leafDecGraph = T.makeLeafGraph inData
          leafGraphSoftWired = POSW.makeLeafGraphSoftWired inData
          charInfoVV = six6 inGraph


          inGraphNetPenalty = if (graphType inGS == Tree) || (graphType inGS == HardWired) then 0.0
                             else if (graphFactor inGS) == NoNetworkPenalty then 0.0
                             else if (graphFactor inGS) == Wheeler2015Network then T.getW15NetPenalty Nothing inGraph
                             else if (graphFactor inGS) == Wheeler2023Network then T.getW23NetPenalty Nothing inGraph
                             else error ("Network penalty type " ++ (show $ graphFactor inGS) ++ " is not yet implemented")
          inGraphNetPenaltyFactor = inGraphNetPenalty / (snd6 inGraph)
      in
      -- trace ("SSPRTBR:" ++ (show inGraphNetPenaltyFactor)) (

      if inSimAnnealParams == Nothing then
        -- trace ("Non SA swap") (
        -- steepest takes immediate best--does not keep equall cost-- for now--disabled not working correctly so goes to "all"
        -- Nothing for SimAnneal Params
         let (swappedGraphs, counter, _) = swapAll swapType inGS inData numToKeep maxMoveEdgeDist steepest alternate 0 (snd6 inGraph) [] [inGraph] numLeaves leafGraph leafDecGraph leafGraphSoftWired charInfoVV doIA inGraphNetPenaltyFactor inSimAnnealParams
         in
         -- trace ("SWAPSPRTBR: " ++ (show $ fmap snd6 swappedGraphs)) (
         if null swappedGraphs then ([inGraph], counter)
         else (swappedGraphs, counter)
         -- )

         
      -- simulated annealing/drifting acceptance does a steepest with SA acceptance
      -- then a swap steepest and all on annealed graph
      -- same at this level method (SA, Drift) choice occurs at lower level
      else
         -- annealed should only yield a single graph
         -- trace ("\tAnnealing/Drifting Swap: " ++ swapType ++ (show (method $ fromJust inSimAnnealParams, numberSteps $ fromJust inSimAnnealParams, currentStep $ fromJust inSimAnnealParams, head $ randomIntegerList $ fromJust inSimAnnealParams, rounds $ fromJust inSimAnnealParams, driftAcceptEqual $ fromJust inSimAnnealParams, driftAcceptWorse $ fromJust inSimAnnealParams, driftMaxChanges $ fromJust inSimAnnealParams, driftChanges $ fromJust inSimAnnealParams))) (
         let -- create list of params with unique list of random values for rounds of annealing
             annealDriftRounds = rounds $ fromJust inSimAnnealParams
             newSimAnnealParamList = U.generateUniqueRandList annealDriftRounds inSimAnnealParams

             -- this to ensure current step set to 0
             (annealDriftGraphs', anealDriftCounter, _) = unzip3 $ (PU.seqParMap rdeepseq (swapAll swapType inGS inData 1 maxMoveEdgeDist True alternate 0 (snd6 inGraph) [] [inGraph] numLeaves leafGraph leafDecGraph leafGraphSoftWired charInfoVV doIA inGraphNetPenaltyFactor) newSimAnnealParamList) -- `using` PU.myParListChunkRDS)

             -- annealed/Drifted 'mutated' graphs
             annealDriftGraphs = take numToKeep $ GO.selectPhylogeneticGraph [("unique","")] 0 ["unique"] $ concat annealDriftGraphs'

             -- swap back "normally" if desired for full drifting/annealing 
             (swappedGraphs, counter, _) = swapAll swapType inGS inData numToKeep maxMoveEdgeDist True alternate 0 (min (snd6 inGraph) (minimum $ fmap snd6 annealDriftGraphs)) [] annealDriftGraphs numLeaves leafGraph leafDecGraph leafGraphSoftWired charInfoVV doIA inGraphNetPenaltyFactor Nothing

             bestGraphs = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] (inGraph : swappedGraphs)
         in
         -- trace ("Steepest SSPRTBR: " ++ (show (length swappedGraphs, counter)))
         --trace ("AC:" ++ (show $ fmap snd6 $ concat annealedGraphs') ++ " -> " ++ (show $ fmap snd6 $ swappedGraphs')) (

         -- this Bool for Genetic Algorithm mutation step
         if not returnMutated then (bestGraphs, counter + (sum anealDriftCounter))
         else (annealDriftGraphs, sum anealDriftCounter)
         -- )
         -- )

-- | swapAll is a high level function that basically deals with portioning out swap-type swaps 
-- and performs the high level options for "alternate" where SPR is perfomred first, then TBR,
-- but whenever a better (or additional) graph is found during TBR, an SPR swap of that graph
-- is performed before returning to TBR again.  THis contibues untill no new graphs are found in the 
-- SPR + TBR swap.
-- each call to swapAll' sets break edge number to 0
swapAll  :: String
         -> GlobalSettings
         -> ProcessedData
         -> Int
         -> Int
         -> Bool
         -> Bool
         -> Int
         -> VertexCost
         -> [PhylogeneticGraph]
         -> [PhylogeneticGraph]
         -> Int
         -> SimpleGraph
         -> DecoratedGraph
         -> DecoratedGraph
         -> V.Vector (V.Vector CharInfo)
         -> Bool
         -> VertexCost
         -> Maybe SAParams
         -> ([PhylogeneticGraph], Int, Maybe SAParams)
swapAll swapType inGS inData numToKeep maxMoveEdgeDist steepest alternate counter curBestCost curSameBetterList inGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor inSimAnnealParams =
   if null inGraphList then (curSameBetterList, counter, inSimAnnealParams)
   else 
      -- nni, spr, tbr
      if swapType /= "alternate" then 
         swapAll' swapType inGS inData numToKeep maxMoveEdgeDist steepest False counter curBestCost curSameBetterList inGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor 0 inSimAnnealParams

      -- alternate
      else 
         let -- spr first
             (sprGraphs, sprCounter, sprSAPArams) = swapAll' "spr" inGS inData numToKeep maxMoveEdgeDist steepest False counter curBestCost curSameBetterList inGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor 0 inSimAnnealParams
             graphsToTBR = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] (sprGraphs ++ inGraphList)
             sprBestCost = (snd6 . head) graphsToTBR
            
             -- tbr until find better or novel equal
             (tbrGraphs, tbrCounter, tbrSAPArams) = swapAll' "tbr" inGS inData numToKeep maxMoveEdgeDist steepest True sprCounter sprBestCost graphsToTBR graphsToTBR numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor 0 sprSAPArams
             tbrBestCost = if (not . null) tbrGraphs then
                              (snd6 . head) tbrGraphs
                           else infinity

         in
         -- if TBR found better go around again with SPR first
         if tbrBestCost < sprBestCost then
            swapAll swapType inGS inData numToKeep maxMoveEdgeDist steepest alternate tbrCounter tbrBestCost tbrGraphs tbrGraphs numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor tbrSAPArams

         -- check if found additional
         else if tbrBestCost == sprBestCost then
            let newTBRGraphs = tbrGraphs `GO.phylogeneticGraphListMinus` sprGraphs
            in
            if (not . null) newTBRGraphs then
               swapAll swapType inGS inData numToKeep maxMoveEdgeDist steepest alternate tbrCounter tbrBestCost tbrGraphs newTBRGraphs numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor tbrSAPArams

            -- found nothing new
            else (tbrGraphs, tbrCounter, tbrSAPArams)

         -- found nothing better or equal
         else (graphsToTBR, tbrCounter, tbrSAPArams)



-- | swapAll' performs branch swapping on all 'break' edges and all readditions
-- this not a "map" version to reduce memory footprint to a more mangeable level
-- "steepest" a passable option to short circuit readdition action to return immediately 
-- if better graph found
-- 1) takes first graph
-- 2) if steepeast checks to make sure <= current best cost
-- 3) gets list of "break-able" edges
--    all non-root edges if tree
--    all non-root bridge edges if network
-- 4) send list (and other info) to split-join function
--    goes on if empty list returned or > current best
--    add graphs todo list if == current best cost
-- 5) returns all of minimum cost found
-- if Alternate then when found better do SPR first then TBR
   -- assumes SPR done before  "alternate" entering so can star with TBR and iff get better
   -- go back to SPR. NBest for "steepest" descent
-- For drift and anneal need to randomize order of splits and rejoins
swapAll' :: String
         -> GlobalSettings
         -> ProcessedData
         -> Int
         -> Int
         -> Bool
         -> Bool
         -> Int
         -> VertexCost
         -> [PhylogeneticGraph]
         -> [PhylogeneticGraph]
         -> Int
         -> SimpleGraph
         -> DecoratedGraph
         -> DecoratedGraph
         -> V.Vector (V.Vector CharInfo)
         -> Bool
         -> VertexCost
         -> Int
         -> Maybe SAParams
         -> ([PhylogeneticGraph], Int, Maybe SAParams)
swapAll' swapType inGS inData numToKeep maxMoveEdgeDist steepest alternate counter curBestCost curSameBetterList inGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor breakEdgeNumber inSimAnnealParams =
   -- trace (" In cost " ++ (show curBestCost) ++ (" " ++ swapType)) (
   -- don't beed to check for mutated here since checked above
   if null inGraphList then
      -- trace (" Out cost " ++ (show curBestCost) ++ (" " ++ swapType))
      (take numToKeep $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] curSameBetterList, counter, inSimAnnealParams)
   else
      let firstGraph = head inGraphList
          firstDecoratedGraph = thd6 firstGraph
          (firstRootIndex, _) = head $ LG.getRoots firstDecoratedGraph

          -- determine edges to break on--'bridge' edges only for network
          -- filter out edges from root since no use--would just rejoin
          -- sort longest edge to shortest--option to speeed up steepest and conditions for all as well
          -- this edge sort from Varon and Wheeler 2013
          breakEdgeList' = if (graphType inGS) == Tree || LG.isTree firstDecoratedGraph then 
                              GO.sortEdgeListByLength $ filter ((/= firstRootIndex) . fst3) $ LG.labEdges firstDecoratedGraph
                          else 
                              GO.sortEdgeListByLength $ filter ((/= firstRootIndex) . fst3) $ LG.getEdgeSplitList firstDecoratedGraph

          -- randomize edges list order for anneal and drift
          breakEdgeList'' = if isJust inSimAnnealParams then 
                              permuteList (head $ randomIntegerList $ fromJust inSimAnnealParams) breakEdgeList' 
                            else breakEdgeList' 

          -- move first "breakEdgeFactor" edges in split list to end
          -- since breakEdgeFactor can get incremented past number of edges the integer remainder is determined
          -- to move to end
          -- this to reduces the revisiting of stable edges (by moving them to the end of the list)
          -- yet still insures that all edges will be visited in final (or ay time needed) split.
          -- used in POY v 1-3
          breakEdgeFactor = snd $ divMod breakEdgeNumber (length breakEdgeList'')
          breakEdgeList =  (drop breakEdgeFactor breakEdgeList'') ++ (take breakEdgeFactor breakEdgeList'')


          -- perform intial split and rejoin on each edge in first graph
          (newGraphList', newSAParams, newBreakEdgeNumber) = splitJoinGraph swapType inGS inData numToKeep maxMoveEdgeDist steepest curBestCost curSameBetterList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor inSimAnnealParams firstGraph breakEdgeNumber breakEdgeList breakEdgeList

          -- get best return graph list-can be empty if nothing better ort smame cost
          newGraphList = GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] newGraphList'

          -- get unique return graph list-can be empty if nothing better ort same cost
          newGraphListUnique = GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] newGraphList'

          newMinCost = if (not . null) newGraphList' then minimum $ fmap snd6 newGraphList'
                       else infinity 

      in  
      -- trace ("Current min cost: "  ++ (show (newMinCost, curBestCost))) (
      -- traceNoLF ("\tBEF:" ++ (show breakEdgeFactor)) (
      -- logic for returning normal swap operations (better only)
      -- versus simulated annealin/Drifing returning potentially sub-optimal
      if isNothing inSimAnnealParams then
         postProcessSwap swapType inGS inData numToKeep maxMoveEdgeDist steepest alternate counter curBestCost curSameBetterList inGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor Nothing newMinCost newBreakEdgeNumber newGraphList

      -- simulated annealing/Drift post processing      
      else 
         postProcessAnnealDrift swapType inGS inData numToKeep maxMoveEdgeDist steepest alternate counter curBestCost curSameBetterList inGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor newSAParams newMinCost newBreakEdgeNumber newGraphListUnique
      -- )

-- | postProcessSwap factors out the post processing of swap results to allow for clearer code 
-- with "regular" optimal swapping
postProcessSwap :: String
               -> GlobalSettings
               -> ProcessedData
               -> Int
               -> Int
               -> Bool
               -> Bool
               -> Int
               -> VertexCost
               -> [PhylogeneticGraph]
               -> [PhylogeneticGraph]
               -> Int
               -> SimpleGraph
               -> DecoratedGraph
               -> DecoratedGraph
               -> V.Vector (V.Vector CharInfo)
               -> Bool
               -> VertexCost
               -> Maybe SAParams
               -> VertexCost
               -> Int 
               -> [PhylogeneticGraph]
               -> ([PhylogeneticGraph], Int, Maybe SAParams)
postProcessSwap swapType inGS inData numToKeep maxMoveEdgeDist steepest alternate counter curBestCost curSameBetterList inGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor inSimAnnealParams newMinCost breakEdgeNumber newGraphList =
   -- found better cost graph
      if newMinCost < curBestCost then
         traceNoLF ("\t->" ++ (show newMinCost))( -- ++ swapType) (
         -- for alternarte do SPR first then TBR
         let graphsToSwap = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] newGraphList -- (newGraphList ++ (tail inGraphList))
         in
         {- This moved to swapAll functionality
         if alternate || swapType == "alternate" then 
            let (sprList, _, _) = swapAll' "spr" inGS inData numToKeep maxMoveEdgeDist steepest alternate (counter + 1) newMinCost newGraphList graphsToSwap numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor inSimAnnealParams
                sprMinCost = min curBestCost ((snd6 . head) sprList)
                graphsToTBRSwap = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] sprList -- (sprList ++ (tail inGraphList))

            in
            swapAll' "tbr" inGS inData numToKeep maxMoveEdgeDist steepest alternate (counter + 2) sprMinCost sprList graphsToTBRSwap numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor inSimAnnealParams
         -}

         -- for alternate if found better return immediately
         if alternate then (newGraphList, counter, inSimAnnealParams)

         -- regular swap--keep going with better graphs
         else swapAll' swapType inGS inData numToKeep maxMoveEdgeDist steepest alternate (counter + 1) newMinCost newGraphList graphsToSwap numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor breakEdgeNumber inSimAnnealParams
         )

      -- found only worse graphs--never happens due to the way splitjoin returns only better or equal
      else if newMinCost > curBestCost then
         -- trace ("Worse " ++ (show newMinCost)) (
         let newCurSameBetterList = GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] (curSameBetterList ++ newGraphList)
         in
         -- traceNoLF ("\tHolding " ++ (show $ length newCurSameBetterList) ++ " at cost "  ++ (show curBestCost) ++ " with " ++ (show $ tail inGraphList) ++ " remaining to " ++ swapType ++ " swap") 

         -- breakEdgeNUmber set to zero for new graph to look at
         swapAll' swapType inGS inData numToKeep maxMoveEdgeDist steepest alternate (counter + 1) curBestCost newCurSameBetterList (tail inGraphList) numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor 0 inSimAnnealParams
         -- )

      -- found same cost graphs
      else
         -- Important to not limit curSameBest since may rediscover graphs via swapping on equal when limiting the number to keep
         -- can be a cause of infinite running issues.
         let newCurSameBetterList = GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] (curSameBetterList ++ newGraphList)
             -- ( _, newNovelGraphList) = unzip $ filter ((== True) .fst) $ zip (fmap (GO.isNovelGraph (curSameBetterList ++ (tail inGraphList))) newGraphList) newGraphList
                -- newNovelGraphList = newGraphList GO.phylogeneticGraphListMinus curSameBetterList
             graphsToDo  = (GO.selectPhylogeneticGraph [("best", "")] 0 ["best"]  $ (tail inGraphList) ++ newGraphList) `GO.phylogeneticGraphListMinus` curSameBetterList
             -- graphsToDo  = (GO.selectPhylogeneticGraph [("best", "")] 0 ["best"]  $ (tail inGraphList) ++ newNovelGraphList) GO.phylogeneticGraphListMinus curSameBetterList
             graphsToDo' = if length graphsToDo >= (numToKeep - 1) then (tail inGraphList)
                           else graphsToDo
               --graphsToDo' = (tail inGraphList)
         in
         -- trace ("Num in best: " ++ (show $ length curSameBetterList) ++ " Num to do: " ++ (show $ length graphsToDo) ++ " from: " ++ (show (length newNovelGraphList, length newGraphList)))
         -- traceNoLF ("\tRemaining to " ++ swapType ++ " swap " ++ (show $ length graphsToDo') ++ " at cost "  ++ (show curBestCost)) (
         
         -- for alternate if equal return immediately
         if alternate then (newCurSameBetterList, counter, inSimAnnealParams)

         -- regular swap--keep going with novel equal cost graphs
         else swapAll' swapType inGS inData numToKeep maxMoveEdgeDist steepest alternate (counter + 1) curBestCost newCurSameBetterList graphsToDo' numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor breakEdgeNumber inSimAnnealParams
         -- )
         
-- | postProcessAnnealDrift factors out the post processing of swap results to allow for clearer code 
-- with simulated annealing/Drifting returns of potentially
-- non optimal graphs
-- if better gaph cost found then update graph and better cost
-- else if graphs probablitstically found--return them if hit max drfit changes or annealing steps
-- else go on with updated number of changes/steps
-- removes alternate and just does straight NNI/SPR/TBR in order to properly do teh returns of altered and better graphs
-- not sure about curBestGraph and return--new Graphs could be empty 
postProcessAnnealDrift :: String
               -> GlobalSettings
               -> ProcessedData
               -> Int
               -> Int
               -> Bool
               -> Bool
               -> Int
               -> VertexCost
               -> [PhylogeneticGraph]
               -> [PhylogeneticGraph]
               -> Int
               -> SimpleGraph
               -> DecoratedGraph
               -> DecoratedGraph
               -> V.Vector (V.Vector CharInfo)
               -> Bool
               -> VertexCost
               -> Maybe SAParams
               -> VertexCost
               -> Int 
               -> [PhylogeneticGraph]
               -> ([PhylogeneticGraph], Int, Maybe SAParams)
postProcessAnnealDrift swapType inGS inData numToKeep maxMoveEdgeDist steepest alternate counter curBestCost curSameBetterList inGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor inSimAnnealParams newMinCost breakEdgeNumber newGraphList =
   -- trace ("PPA: " ++ (show (newMinCost, curBestCost)) ++ " counts " ++ (show ((currentStep $ fromJust inSimAnnealParams), (numberSteps $ fromJust inSimAnnealParams), (driftChanges $ fromJust inSimAnnealParams), (driftMaxChanges $ fromJust inSimAnnealParams)))) (
   -- found better cost graph
      if newMinCost < curBestCost then
         traceNoLF ("\t->" ++ (show newMinCost)) 
         (newGraphList, counter, inSimAnnealParams)
         {-
         -- for alternate do SPR first then TBR
         if not alternate then 
            swapAll' swapType inGS inData numToKeep maxMoveEdgeDist steepest alternate (counter + 1) newMinCost newGraphList (newGraphList ++ (tail inGraphList)) numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor inSimAnnealParams

         else swapAll' "tbr" inGS inData numToKeep maxMoveEdgeDist steepest alternate (counter + 1) newMinCost newGraphList (newGraphList ++ (tail inGraphList)) numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor inSimAnnealParams
         -}

      -- not better so check for drift changes or annealing steps and return if reached maximum number
      else if ((currentStep $ fromJust inSimAnnealParams) >= (numberSteps $ fromJust inSimAnnealParams)) || ((driftChanges $ fromJust inSimAnnealParams) >= (driftMaxChanges $ fromJust inSimAnnealParams)) then 
         --trace ("PPA return: " ++ (show (newMinCost, curBestCost))) 
         (take numToKeep $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (newGraphList ++ curSameBetterList), counter, inSimAnnealParams)

      -- didn't hit stopping numbers so continuing--but based on current best cost not whatever was found
      else 
         swapAll' swapType inGS inData numToKeep maxMoveEdgeDist steepest alternate (counter + 1) curBestCost (newGraphList ++ curSameBetterList) (newGraphList ++ (tail inGraphList))  numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor breakEdgeNumber inSimAnnealParams

      -- )

-- | splitJoinGraph splits a graph on a single input edge (recursively though edge list) and rejoins to all possible other edges
-- if steepest == True then returns on finding a better graph (lower cost)
-- this will traverse entire SPR neighbohood if nothing better found (or steepest == False)
-- different from swapALL (original) in that it doesn't build up split list so lower memory footprint
-- breakEdgeList Complete keeps original edge list so can create readdition edge lists more easily
-- parallel map on rejoin if not steepest, if steepest do number of parallel threads so can reurn if any one is better
-- NB -- need to verify NNI/SPR/TBR rearrangement numbers
-- assumes break edges are bridge edges
-- graph split into two peices "base" graph with original root and "pruned" graph that was split off.
-- the edge connecting the two is (originalConnectionOfPruned -> prunedGraphRootIndex)
-- the edges in pruned graph do not contaiun that edge since are enumerated via preorder pass from prunedGraphRootIndex
-- this is the edge that is reconnected when graphs are joined, it is often delted and rejoined to update info and to deal with
-- conditions where the pruned graph is a single terminal
-- returns teh "breakEdgeNumber" so that in steepest, the edge breaking can continue in where it left off so to speak.
-- this can speed up SPR/TBR by a contant factor by not revisiting stable edges (not used in SA/drifting)
-- used in POY v1-3
splitJoinGraph :: String
               -> GlobalSettings
               -> ProcessedData
               -> Int
               -> Int
               -> Bool
               -> VertexCost
               -> [PhylogeneticGraph]
               -> Int
               -> SimpleGraph
               -> DecoratedGraph
               -> DecoratedGraph
               -> V.Vector (V.Vector CharInfo)
               -> Bool
               -> VertexCost
               -> Maybe SAParams
               -> PhylogeneticGraph
               -> Int 
               -> [LG.LEdge EdgeInfo]
               -> [LG.LEdge EdgeInfo]
               -> ([PhylogeneticGraph], Maybe SAParams, Int)
splitJoinGraph swapType inGS inData numToKeep maxMoveEdgeDist steepest curBestCost curSameBetterList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor inSimAnnealParams firstGraph breakEdgeNumber' breakEdgeListComplete breakEdgeList =
   if null breakEdgeList then (curSameBetterList, inSimAnnealParams, 0)
   else 
      -- split on first input edge
      let edgeToBreakOn = head breakEdgeList 

          -- this so breaking edges can contnue where current left off 
          -- since "rotates" edges all will be done.
          breakEdgeNumber = breakEdgeNumber' + 1

          -- split input graph into a part with the original root ("base") and the "pruned" graph -- the piece split off w/o original root
          (splitGraph, graphRoot, prunedGraphRootIndex,  originalConnectionOfPruned) = LG.splitGraphOnEdge (thd6 firstGraph) edgeToBreakOn

          -- reoptimize split graph for re-addition heuristics
          (reoptimizedSplitGraph, splitCost) = reoptimizeSplitGraphFromVertex inGS inData doIA netPenaltyFactor splitGraph graphRoot prunedGraphRootIndex

          -- get root in base (for readdition) and edges in pruned section for rerooting during TBR readdition
          (_, edgesInPrunedGraph) = LG.nodesAndEdgesAfter splitGraph [(originalConnectionOfPruned, fromJust $ LG.lab splitGraph originalConnectionOfPruned)]
          
          edgesInBaseGraph = breakEdgeListComplete L.\\ (edgeToBreakOn : edgesInPrunedGraph)
          
          -- determine those edges within distance of original if limited (ie NNI etc)
          rejoinEdges' = if maxMoveEdgeDist >= ((maxBound :: Int) `div` 3) then edgesInBaseGraph
                        else take maxMoveEdgeDist $ (LG.sortEdgeListByDistance splitGraph [graphRoot] [graphRoot]) 

          -- randomize edges list order for anneal and drift
          rejoinEdges = if isJust inSimAnnealParams then 
                           permuteList ((randomIntegerList $ fromJust inSimAnnealParams) !! 1) rejoinEdges' 
                          else rejoinEdges' 

          -- rejoin graph to all possible edges in base graph
          (newGraphList, newSAParams) = 
            {-
            trace ("Edge to break on:" ++ (show $ LG.toEdge edgeToBreakOn) 
            ++ "\nBase graph edges: " ++ (show $ fmap LG.toEdge edgesInBaseGraph) 
            ++ "\nPruned graph edges: " ++ (show $ fmap LG.toEdge edgesInPrunedGraph) 
            ++ "\nTarget edges to rejoin: " ++ (show $ fmap LG.toEdge rejoinEdges) 
            ++ "\nFull edgelist: " ++ (show $ fmap LG.toEdge breakEdgeListComplete)) 
            -}
            rejoinGraph swapType inGS inData numToKeep maxMoveEdgeDist steepest curBestCost [] doIA netPenaltyFactor reoptimizedSplitGraph (GO.convertDecoratedToSimpleGraph splitGraph) splitCost graphRoot prunedGraphRootIndex originalConnectionOfPruned rejoinEdges edgesInPrunedGraph charInfoVV inSimAnnealParams

          newGraphList' = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] newGraphList
      in
      
      -- regular swap
      if isNothing inSimAnnealParams then
         -- only returns graphs if same of better else empty
         -- adds null o\r better graphs to reurn list 
         if (not . null) newGraphList && steepest then (newGraphList', inSimAnnealParams, breakEdgeNumber)
         else 
            let (recurseGraphList, _, newEdgeBreakNumber) = splitJoinGraph swapType inGS inData numToKeep maxMoveEdgeDist steepest curBestCost curSameBetterList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor inSimAnnealParams firstGraph breakEdgeNumber breakEdgeListComplete (tail breakEdgeList)
            in
            (newGraphList' ++ recurseGraphList, inSimAnnealParams, newEdgeBreakNumber) 

      -- Annealing/Drift swap
         -- return if better
         -- return if other chosen probabalistically
         -- recurse if nothing returned
      else
         -- if better than current graph return
         let newMinCost = if (not . null) newGraphList then (snd6 . head) newGraphList'
                          else infinity
         in
         if newMinCost < curBestCost then (newGraphList', newSAParams, breakEdgeNumber)

         -- if SA returned graphs--return them
         else if (not . null) newGraphList then (newGraphList, newSAParams, breakEdgeNumber)

         -- keep going if nothing
         else
            splitJoinGraph swapType inGS inData numToKeep maxMoveEdgeDist steepest curBestCost curSameBetterList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired charInfoVV doIA netPenaltyFactor newSAParams firstGraph breakEdgeNumber breakEdgeListComplete (tail breakEdgeList)



-- | rejoinGraphTuple is a wrapper around rejoinGraph for fmapping--only returns graph list not simulated annealing params
rejoinGraphTuple :: String
                 -> GlobalSettings
                 -> ProcessedData
                 -> Int
                 -> Int
                 -> Bool
                 -> VertexCost
                 -> [PhylogeneticGraph]
                 -> Bool
                 -> V.Vector (V.Vector CharInfo)
                 -> Maybe SAParams
                 -> (DecoratedGraph, SimpleGraph, VertexCost, LG.Node,LG.Node, LG.Node, [LG.LEdge EdgeInfo], [LG.LEdge EdgeInfo], VertexCost)
                 -> [PhylogeneticGraph]
rejoinGraphTuple swapType inGS inData numToKeep maxMoveEdgeDist steepest curBestCost curBestGraphs doIA charInfoVV inSimAnnealParams (reoptimizedSplitGraph, splitGraphSimple, splitGraphCost, graphRoot, prunedGraphRootIndex, originalConnectionOfPruned, rejoinEdges, edgesInPrunedGraph, netPenaltyFactor) = 
   fst $ rejoinGraph swapType inGS inData numToKeep maxMoveEdgeDist steepest curBestCost curBestGraphs doIA netPenaltyFactor reoptimizedSplitGraph splitGraphSimple splitGraphCost graphRoot prunedGraphRootIndex originalConnectionOfPruned rejoinEdges edgesInPrunedGraph charInfoVV inSimAnnealParams 

-- | rejoinGraph rejoins a split graph at all edges (if not steepest and found better)
-- in "base" graph.
-- if not steepest then do all as map, else recursive on base graph edge list 
-- nni doesn't apper to be correct here--maybe loose it--doing nothing
rejoinGraph :: String
            -> GlobalSettings
            -> ProcessedData
            -> Int
            -> Int
            -> Bool
            -> VertexCost
            -> [PhylogeneticGraph]
            -> Bool
            -> VertexCost
            -> DecoratedGraph
            -> SimpleGraph
            -> VertexCost
            -> LG.Node
            -> LG.Node
            -> LG.Node
            -> [LG.LEdge EdgeInfo]
            -> [LG.LEdge EdgeInfo]
            -> V.Vector (V.Vector CharInfo)
            -> Maybe SAParams
            -> ([PhylogeneticGraph], Maybe SAParams)
rejoinGraph swapType inGS inData numToKeep maxMoveEdgeDist steepest curBestCost curBestGraphs doIA netPenaltyFactor reoptimizedSplitGraph splitGraphSimple splitGraphCost graphRoot prunedGraphRootIndex originalConnectionOfPruned rejoinEdges' edgesInPrunedGraph charInfoVV inSimAnnealParams =

   -- found no better--but return equal cost graphs
   -- trace ("In rejoinGraph with num rejoining edges: " ++ (show $ length rejoinEdges)) (
   if null rejoinEdges' then (curBestGraphs, inSimAnnealParams)

   else
      -- this is for no  swapping option in fuse and genetic algorithm-fuse
      let rejoinEdges = if swapType == "none" then take 6 rejoinEdges'
                        else rejoinEdges'
      in
      -- regular swapping
      if isNothing inSimAnnealParams then
      
         -- check if split graph cost same as graph then return since graph can only get longer on readdition 
         if splitGraphCost >= curBestCost then ([], inSimAnnealParams)

         else 
            -- fmap over all edges in base graph
            if not steepest then 
               let -- rejoinGraphList = concatMap (singleJoin swapType steepest inGS inData reoptimizedSplitGraph splitGraphSimple splitGraphCost doIA prunedGraphRootIndex originalConnectionOfPruned charInfoVV curBestCost edgesInPrunedGraph) rejoinEdges `using` PU.myParListChunkRDS
                   rejoinGraphList = concat $ fmap fst $ PU.seqParMap rdeepseq (singleJoin swapType steepest inGS inData reoptimizedSplitGraph splitGraphSimple splitGraphCost doIA prunedGraphRootIndex originalConnectionOfPruned charInfoVV curBestCost edgesInPrunedGraph inSimAnnealParams) rejoinEdges 
                   
                   {-Checking only min but seems to make slower
                   newMinCost = if null rejoinGraphList then infinity
                                else minimum $ fmap snd rejoinGraphList
                   (minEstCostNewGraphList, _) = unzip $ filter ((== newMinCost) . snd) rejoinGraphList
                   -}
                   
                   -- newGraphList = fmap (T.multiTraverseFullyLabelGraph inGS inData False False Nothing) (fmap fst rejoinGraphList) `using` PU.myParListChunkRDS
                   newGraphList' = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] rejoinGraphList -- newGraphList
               in
               -- will only return graph if <= curBest cost
               if null rejoinGraphList then ([], inSimAnnealParams)
               else if (snd6 . head) newGraphList' <= curBestCost then (newGraphList', inSimAnnealParams)
               else ([], inSimAnnealParams)
               
            -- famp over number of threads edges in base graph
            -- then recurse
            else 
               -- trace ("In steepest: " ++ (show PU.getNumThreads) ++ " " ++ (show $ length $ take PU.getNumThreads rejoinEdges)) (
               let -- this could be made a little paralle--but if lots of threads basically can do all 
                   -- to not overload paralle threads
                   {-  This not so efficient is swapping in single graphs so leaving it be
                   saRounds = if isNothing inSimAnnealParams then 1
                              else rounds $ fromJust inSimAnnealParams
                        
                   (numGraphsToExamine, _) = divMod PU.getNumThreads saRounds -- this may not "drift" if finds alot better, but that's how its supposed to work
                   -}
                   numGraphsToExamine = min (graphsSteepest inGS) PU.getNumThreads
                   rejoinEdgeList = take numGraphsToExamine rejoinEdges
                   --rejoinGraphList = concatMap (singleJoin swapType steepest inGS inData reoptimizedSplitGraph splitGraphSimple splitGraphCost doIA prunedGraphRootIndex originalConnectionOfPruned charInfoVV curBestCost edgesInPrunedGraph) rejoinEdgeList `using` PU.myParListChunkRDS
                   rejoinGraphList = concat $ fmap fst $ PU.seqParMap rdeepseq (singleJoin swapType steepest inGS inData reoptimizedSplitGraph splitGraphSimple splitGraphCost doIA prunedGraphRootIndex originalConnectionOfPruned charInfoVV curBestCost edgesInPrunedGraph inSimAnnealParams) rejoinEdgeList 
                  
                   {--Checking only min but seems to make slower
                   newMinCost = if null rejoinGraphList then infinity
                                else minimum $ fmap snd rejoinGraphList
                   (minEstCostNewGraphList, _) = unzip $ filter ((== newMinCost) . snd) rejoinGraphList
                   -}
                  
                   -- newGraphList = fmap (T.multiTraverseFullyLabelGraph inGS inData False False Nothing) (fmap fst rejoinGraphList) `using` PU.myParListChunkRDS
                   newGraphList' = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] rejoinGraphList -- newGraphList
               in
               -- found nothing better or equal 
               if null rejoinGraphList then 
                  -- trace ("In steepest worse: " ++ (show $ length (drop PU.getNumThreads rejoinEdges))) 
                   rejoinGraph swapType inGS inData numToKeep maxMoveEdgeDist steepest curBestCost curBestGraphs doIA netPenaltyFactor reoptimizedSplitGraph splitGraphSimple splitGraphCost graphRoot prunedGraphRootIndex originalConnectionOfPruned (drop numGraphsToExamine rejoinEdges) edgesInPrunedGraph charInfoVV inSimAnnealParams

               -- found better graph
               else if (snd6 . head) newGraphList'  < curBestCost then 
                  -- trace ("Steepest better") 
                  (newGraphList', inSimAnnealParams)

               -- found equal cost graph 
               else if (snd6 . head) newGraphList' == curBestCost then 
                  let newBestList = take numToKeep $ GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] (curBestGraphs ++ newGraphList') 
                  in
                  rejoinGraph swapType inGS inData numToKeep maxMoveEdgeDist steepest curBestCost newBestList doIA netPenaltyFactor reoptimizedSplitGraph splitGraphSimple splitGraphCost graphRoot prunedGraphRootIndex originalConnectionOfPruned (drop numGraphsToExamine rejoinEdges) edgesInPrunedGraph charInfoVV inSimAnnealParams 
               -- found worse graphs only
               else 
                  -- trace ("In steepest worse (after recalculation): " ++ (show $ length (drop PU.getNumThreads rejoinEdges))) 
                  rejoinGraph swapType inGS inData numToKeep maxMoveEdgeDist steepest curBestCost curBestGraphs doIA netPenaltyFactor reoptimizedSplitGraph splitGraphSimple splitGraphCost graphRoot prunedGraphRootIndex originalConnectionOfPruned (drop numGraphsToExamine rejoinEdges) edgesInPrunedGraph charInfoVV inSimAnnealParams
               -- ))

      -- Drifting/Simulated annealing 
      -- basically if it accepted a graph (better or probabalistically) then pass up
      -- otherwise move to next rejoin, changes in graphs counted at higher level
      else 
         let --based on "steepest"
             -- to not overload paralle threads
                   {-  This not so efficient is swapping in single graphs so leaving it be
                   saRounds = if isNothing inSimAnnealParams then 1
                              else rounds $ fromJust inSimAnnealParams
                        
                   (numGraphsToExamine, _) = divMod PU.getNumThreads saRounds -- this may not "drift" if finds alot better, but that's how its supposed to work
                   -}
             numGraphsToExamine = min (graphsSteepest inGS) PU.getNumThreads
             rejoinEdgeList = take numGraphsToExamine rejoinEdges
             simAnnealParamList = U.generateUniqueRandList numGraphsToExamine inSimAnnealParams
             rejoinGraphPairList = PU.seqParMap rdeepseq (singleJoin' swapType steepest inGS inData reoptimizedSplitGraph splitGraphSimple splitGraphCost doIA prunedGraphRootIndex originalConnectionOfPruned charInfoVV curBestCost edgesInPrunedGraph) (zip simAnnealParamList rejoinEdgeList)

             -- mechanics to see if trhere is a better graph in return set
             -- only taking first of each list--so can keep sa params with them--really all should have length == 1 anyway
             -- making sure remove all null lists that nothing was found
             nonEmptyPairList = filter (not . null . fst) rejoinGraphPairList

             rejoinGraphList = if (not . null) nonEmptyPairList then fmap (head . fst) nonEmptyPairList
                               else []

             newMinCost = if (not . null) rejoinGraphList then minimum $ fmap snd6 rejoinGraphList
                          else infinity

             -- head should only be called when non-empty--so should never get runtime head error
             (newMinGraph, newMinGraphSAParams) = head $ filter ((== newMinCost) . snd6 . fst) (zip rejoinGraphList (fmap snd nonEmptyPairList))
         in

         -- if better than current--pass up and on
         if newMinCost < curBestCost then ([newMinGraph], newMinGraphSAParams)

         -- check if hit step limit--more for SA than drift
         else if ((currentStep $ fromJust inSimAnnealParams) >= (numberSteps $ fromJust inSimAnnealParams)) || ((driftChanges $ fromJust inSimAnnealParams) >= (driftMaxChanges $ fromJust inSimAnnealParams)) then 
            (curBestGraphs, inSimAnnealParams)

         -- not better so go to SA results
         else 
            -- return first non-empty result
            if (not . null) nonEmptyPairList then head nonEmptyPairList

            -- if nothing returned (no better or probabalistically chosen) go on with updated SA params
            else
               let newSAParams = (snd . head) rejoinGraphPairList
               in
               rejoinGraph swapType inGS inData numToKeep maxMoveEdgeDist steepest curBestCost curBestGraphs doIA netPenaltyFactor reoptimizedSplitGraph splitGraphSimple splitGraphCost graphRoot prunedGraphRootIndex originalConnectionOfPruned (drop numGraphsToExamine rejoinEdges) edgesInPrunedGraph charInfoVV newSAParams

         
                
-- | singleJoin' is a wrapper arounds singleJoin to allow parMap with individual SAParams
singleJoin' :: String 
           -> Bool
           -> GlobalSettings
           -> ProcessedData
           -> DecoratedGraph
           -> SimpleGraph
           -> VertexCost 
           -> Bool
           -> LG.Node
           -> LG.Node
           -> V.Vector (V.Vector CharInfo)
           -> VertexCost
           -> [LG.LEdge EdgeInfo]
           -> (Maybe SAParams, LG.LEdge EdgeInfo)
           -> ([PhylogeneticGraph], Maybe SAParams)
singleJoin' swapType steepest inGS inData splitGraph splitGraphSimple splitCost doIA prunedGraphRootIndex originalConnectionOfPruned charInfoVV curBestCost edgesInPrunedGraph (inSimAnnealParams, targetEdge) =
   singleJoin swapType steepest inGS inData splitGraph splitGraphSimple splitCost doIA prunedGraphRootIndex originalConnectionOfPruned charInfoVV curBestCost edgesInPrunedGraph inSimAnnealParams targetEdge

-- | singleJoin takes optimized split graph, split cost, target edge, swap type (ie TBR/SPR/NNI)
-- and "rejoins" the split graph to a single graph--creates joined graph and calculates a heuristic graph cost 
-- based on the union assignment of the edge and its distance to the root vertex of the pruned graph
-- if TBR checks all edges in pruned graph with readdition edge (shortcircuits if steepest  == True)
-- always deletes connecting edge to pruned part and readds--this because sometimes it is there and sometimes not (depending on 
-- if SPR for terminal etc) and can create parallel edges with different weights (0.0 or not) so just remove to be sure.
-- TBR uses dynamic epsilon even in SPR moves--SPR does not
singleJoin :: String 
           -> Bool
           -> GlobalSettings
           -> ProcessedData
           -> DecoratedGraph
           -> SimpleGraph
           -> VertexCost 
           -> Bool
           -> LG.Node
           -> LG.Node
           -> V.Vector (V.Vector CharInfo)
           -> VertexCost
           -> [LG.LEdge EdgeInfo]
           -> Maybe SAParams
           -> LG.LEdge EdgeInfo 
           -> ([PhylogeneticGraph], Maybe SAParams)
singleJoin swapType steepest inGS inData splitGraph splitGraphSimple splitCost doIA prunedGraphRootIndex originalConnectionOfPruned charInfoVV curBestCost edgesInPrunedGraph inSimAnnealParams targetEdge@(u,v, _) = 
   -- trace ("Rejoinging: " ++ (show $ LG.toEdge targetEdge)) (
   let newEdgeList = [(u, originalConnectionOfPruned, 0.0),(originalConnectionOfPruned, v, 0.0),(originalConnectionOfPruned, prunedGraphRootIndex, 0.0)]
       targetEdgeData = M.makeEdgeData doIA splitGraph charInfoVV targetEdge

       --this for SPR/NNI only
       prunedRootVertexData = vertData $ fromJust $ LG.lab splitGraph prunedGraphRootIndex

       sprReJoinCost = edgeJoinDelta doIA charInfoVV prunedRootVertexData targetEdgeData

       sprNewGraph = LG.insEdges newEdgeList $ LG.delEdges [(u,v),(originalConnectionOfPruned, prunedGraphRootIndex)] splitGraphSimple

       -- here when needed
       rediagnosedSPRGraph = T.multiTraverseFullyLabelGraph inGS inData False False Nothing sprNewGraph
      
       -- Filter for bridge edges fir TBR when needed
       edgesInPrunedGraph' = if (graphType inGS == Tree) || LG.isTree splitGraphSimple then edgesInPrunedGraph
                             else fmap fst $ filter ((== True) . snd) $ zip edgesInPrunedGraph (fmap (LG.isBridge splitGraphSimple) (fmap LG.toEdge edgesInPrunedGraph))
   in

   -- do redo orginal graph join
   if originalConnectionOfPruned `elem` [u,v] then ([], inSimAnnealParams)

   -- regular swap
   else if isNothing inSimAnnealParams then

      -- SPR or no TBR rearrangements
      if (swapType == "spr") || ((length edgesInPrunedGraph) < 4) then 
         if (sprReJoinCost + splitCost) <= curBestCost then 
            if (graphType inGS /= Tree) && ((not . LG.isGraphTimeConsistent) sprNewGraph)  then ([], inSimAnnealParams)
            else if snd6 rediagnosedSPRGraph <= curBestCost then ([rediagnosedSPRGraph], inSimAnnealParams)
            else ([], inSimAnnealParams)
         else ([], inSimAnnealParams)

      else -- TBR 
         
         -- do TBR stuff returning SPR results if heuristic better
         let sprResult = if (sprReJoinCost + splitCost) <= (curBestCost * (dynamicEpsilon inGS)) then 
                           if (graphType inGS /= Tree) && ((not . LG.isGraphTimeConsistent) sprNewGraph)  then []
                           else if snd6 rediagnosedSPRGraph <= curBestCost then [rediagnosedSPRGraph]
                           else []
                         else []
             (tbrResult, _) = tbrJoin steepest inGS inData splitGraph splitGraphSimple splitCost doIA prunedGraphRootIndex originalConnectionOfPruned charInfoVV curBestCost edgesInPrunedGraph' inSimAnnealParams targetEdge
         in
         if (not . null) sprResult then (sprResult, inSimAnnealParams)

         -- unions implementation here
         -- 1.17 from Varon and Wheeler 2013
         else if sprReJoinCost > 1.17 * (curBestCost - splitCost) then ([], inSimAnnealParams)

         -- else if ((snd6 rediagnosedSPRGraph) - curBestCost) > 1.17 * (curBestCost - splitCost) then ([], inSimAnnealParams)
         
         else (tbrResult, inSimAnnealParams)

   -- simulated annealing/Drift swap
   else 
      -- check if spr better always return if so
      let sprResult = if (sprReJoinCost + splitCost) <= curBestCost then 
                        if (graphType inGS /= Tree) && ((not . LG.isGraphTimeConsistent) sprNewGraph)  then []
                        else if snd6 rediagnosedSPRGraph < curBestCost then [rediagnosedSPRGraph]
                        else []
                      else []
      in
      -- spr better than current
      if (not . null) sprResult then
         let (_, newSAParams) = U.simAnnealAccept inSimAnnealParams curBestCost (snd6 rediagnosedSPRGraph)
         in
         -- trace ("SPR found better " ++ (show $ snd6 rediagnosedSPRGraph))
         (sprResult, newSAParams)

      -- do simAnneal/Drift for SPR and on to tbr if not accept
      else
         let (acceptGraph, newSAParams) = U.simAnnealAccept inSimAnnealParams curBestCost (sprReJoinCost + splitCost)
         in

         -- if accepted (better or random) then return with updated annealing/Drift parameters
         if acceptGraph then 
            {-
            let banner = if snd6 rediagnosedSPRGraph <  curBestCost then "SPR heur better"
                         else "Accepted not better"
            in 
            trace banner
            -}
            ([rediagnosedSPRGraph], newSAParams)

         -- rejected--recurse with updated SA params
         -- SPR or small prune
         else if (swapType == "spr") || ((length edgesInPrunedGraph) < 4) then 
            ([], newSAParams)

         -- tbr
         else 
            tbrJoin steepest inGS inData splitGraph splitGraphSimple splitCost doIA prunedGraphRootIndex originalConnectionOfPruned charInfoVV curBestCost edgesInPrunedGraph' newSAParams targetEdge
   -- )

-- | edgeJoinDelta calculates heuristic cost for jopineing pair edges
edgeJoinDelta :: Bool -> V.Vector (V.Vector CharInfo) -> VertexBlockData -> VertexBlockData -> VertexCost
edgeJoinDelta doIA charInfoVV edgeA edgeB =
   if (not doIA) then V.sum $ fmap V.sum $ fmap (fmap snd) $ POS.createVertexDataOverBlocks edgeA edgeB charInfoVV []
   else V.sum $ fmap V.sum $ fmap (fmap snd) $ POS.createVertexDataOverBlocksStaticIA edgeA edgeB charInfoVV []


-- | tbrJoin performs TBR rearrangements on pruned graph component
-- "reroots" pruned graph on each bridge edge and tries join to 
-- target edge as in SPR
-- each edge is tried in turn (except for original root edge covered by singleJoin SPR function)
-- if heuristic edge join cost is below current best cost then the component is rerooted, joined to 
-- target edge and graph fully diagnosed to verify cost
-- "steepest" short circuits checking edges to return better verified cost graph immediately
-- otherwise ("all") returns all graphs better than or equal to current better cost
-- if nothing equal or better found returns empty list 
-- tests if reroot edges are bridge edges
-- uses dynamic epsilon--seems the delta estimate is high
tbrJoin :: Bool
        -> GlobalSettings
        -> ProcessedData
        -> DecoratedGraph
        -> SimpleGraph
        -> VertexCost 
        -> Bool
        -> LG.Node
        -> LG.Node
        -> V.Vector (V.Vector CharInfo)
        -> VertexCost
        -> [LG.LEdge EdgeInfo]
        -> Maybe SAParams
        -> LG.LEdge EdgeInfo 
        -> ([PhylogeneticGraph], Maybe SAParams)
tbrJoin steepest inGS inData splitGraph splitGraphSimple splitCost doIA prunedGraphRootIndex originalConnectionOfPruned charInfoVV curBestCost edgesInPrunedGraph inSimAnnealParams targetEdge = 
   -- trace ("In tbrJoin: to " ++ (show $ LG.toEdge targetEdge)) (
   if null edgesInPrunedGraph then ([], inSimAnnealParams)
   else 
      -- get target edge data
      let targetEdgeData = M.makeEdgeData doIA splitGraph charInfoVV targetEdge
      in  
      
      -- logic for annealing/Drift  regular swap first
      if isNothing inSimAnnealParams then 
         if not steepest then       
            -- get heuristic delta joins for edges in pruned graph
            let rerootEdgeList = filter ((/= prunedGraphRootIndex) . fst3) $ filter ((/= originalConnectionOfPruned) . fst3) edgesInPrunedGraph
                rerootEdgeDataList = PU.seqParMap rdeepseq (M.makeEdgeData doIA splitGraph charInfoVV) rerootEdgeList
                rerootEdgeDeltaList = fmap (+ splitCost) $ PU.seqParMap rdeepseq (edgeJoinDelta doIA charInfoVV targetEdgeData) rerootEdgeDataList
                   
                -- check for possible better/equal graphs and verify
                candidateEdgeList = fmap fst $ filter ((<= (curBestCost * (dynamicEpsilon inGS))) . snd) (zip rerootEdgeList rerootEdgeDeltaList)
                candidateJoinedGraphList = PU.seqParMap rdeepseq (rerootPrunedAndMakeGraph  splitGraphSimple prunedGraphRootIndex originalConnectionOfPruned targetEdge) candidateEdgeList
                rediagnosedGraphList = filter ((<= curBestCost) . snd6) $ PU.seqParMap rdeepseq (T.multiTraverseFullyLabelGraph inGS inData False False Nothing) candidateJoinedGraphList

                -- for debugging
                -- allRediagnosedList = PU.seqParMap rdeepseq (T.multiTraverseFullyLabelGraph inGS inData False False Nothing) (PU.seqParMap rdeepseq (rerootPrunedAndMakeGraph  splitGraphSimple  prunedGraphRootIndex originalConnectionOfPruned targetEdge) rerootEdgeList)
                   
            in
            {-
            trace ("TBR All equal/better: " ++ (show curBestCost) ++ " " ++ (show rerootEdgeDeltaList) ++ " -> " 
               ++ (show $ fmap snd6 allRediagnosedList) 
               ++ " " ++ (show (length rerootEdgeList, length candidateEdgeList, length rediagnosedGraphList))) (
            -}
            if null candidateEdgeList then ([], Nothing)
            else if null rediagnosedGraphList then ([], Nothing)
            else (rediagnosedGraphList, Nothing)
            -- )
         else 
            -- trace ("TBR steepest") (
            -- get steepest edges
            let -- to not overload paralle threads
                {-  This not so efficient is swapping in single graphs so leaving it be
                saRounds = if isNothing inSimAnnealParams then 1
                           else rounds $ fromJust inSimAnnealParams
                     
                (numGraphsToExamine, _) = divMod PU.getNumThreads saRounds -- this may not "drift" if finds alot better, but that's how its supposed to work
                -}
                numEdgesToExamine = min (graphsSteepest inGS) PU.getNumThreads
                firstSetEdges = take numEdgesToExamine edgesInPrunedGraph

                -- get heuristic delta joins for steepest edge set
                rerootEdgeList = filter ((/= prunedGraphRootIndex) . fst3) $ filter ((/= originalConnectionOfPruned) . fst3) firstSetEdges
                rerootEdgeDataList = PU.seqParMap rdeepseq (M.makeEdgeData doIA splitGraph charInfoVV) rerootEdgeList
                rerootEdgeDeltaList = fmap (+ splitCost) $ PU.seqParMap rdeepseq (edgeJoinDelta doIA charInfoVV targetEdgeData) rerootEdgeDataList
                
                -- check for possible better/equal graphs and verify
                candidateEdgeList = fmap fst $ filter ((<= (curBestCost * (dynamicEpsilon inGS))) . snd) (zip rerootEdgeList rerootEdgeDeltaList)
                candidateJoinedGraphList = PU.seqParMap rdeepseq (rerootPrunedAndMakeGraph splitGraphSimple prunedGraphRootIndex originalConnectionOfPruned targetEdge) candidateEdgeList
                rediagnosedGraphList = filter ((<= curBestCost) . snd6) $ PU.seqParMap rdeepseq (T.multiTraverseFullyLabelGraph inGS inData False False Nothing) candidateJoinedGraphList-- get 

            in
            -- trace ("TBR steepest: " ++ (show $ length rerootEdgeList) ++ " edges to go " ++ (show $ length $ (drop numEdgesToExamine edgesInPrunedGraph))) (
            if null candidateEdgeList then 
               tbrJoin steepest inGS inData splitGraph splitGraphSimple splitCost doIA prunedGraphRootIndex originalConnectionOfPruned charInfoVV curBestCost (drop numEdgesToExamine edgesInPrunedGraph) inSimAnnealParams targetEdge
            else if null rediagnosedGraphList then 
               tbrJoin steepest inGS inData splitGraph splitGraphSimple splitCost doIA prunedGraphRootIndex originalConnectionOfPruned charInfoVV curBestCost (drop numEdgesToExamine edgesInPrunedGraph) inSimAnnealParams targetEdge
            else 
               -- trace ("TBR: " ++ (show $ minimum $ fmap snd6 rediagnosedGraphList))
               (rediagnosedGraphList, Nothing)    
            -- )

         -- simulated annealing/Drift stuff
         -- based on steepest type swapping
         else
            -- trace ("TBR SA/Drift") (
            let -- to not overload paralle threads
                {-  This not so efficient is swapping in single graphs so leaving it be
                saRounds = if isNothing inSimAnnealParams then 1
                           else rounds $ fromJust inSimAnnealParams
                     
                (numGraphsToExamine, _) = divMod PU.getNumThreads saRounds -- this may not "drift" if finds alot better, but that's how its supposed to work
                -}
                numEdgesToExamine = min (graphsSteepest inGS) PU.getNumThreads
                firstSetEdges = take numEdgesToExamine edgesInPrunedGraph

                -- get heuristic delta joins for steepest edge set
                rerootEdgeList = filter ((/= prunedGraphRootIndex) . fst3) $ filter ((/= originalConnectionOfPruned) . fst3) firstSetEdges
                rerootEdgeDataList = PU.seqParMap rdeepseq (M.makeEdgeData doIA splitGraph charInfoVV) rerootEdgeList
                rerootEdgeDeltaList = fmap (+ splitCost) $ PU.seqParMap rdeepseq (edgeJoinDelta doIA charInfoVV targetEdgeData) rerootEdgeDataList

                minDelta = if (not . null) rerootEdgeDeltaList then minimum $ rerootEdgeDeltaList
                           else infinity
                minEdgeList = if (not . null) rerootEdgeDeltaList then fmap fst $ filter ((== minDelta) . snd)  (zip rerootEdgeList rerootEdgeDeltaList)
                              else []
                
                -- check for possible better/equal graphs and verify
                candidateJoinedGraphList = PU.seqParMap rdeepseq (rerootPrunedAndMakeGraph splitGraphSimple prunedGraphRootIndex originalConnectionOfPruned targetEdge) minEdgeList
                rediagnosedGraphList = PU.seqParMap rdeepseq (T.multiTraverseFullyLabelGraph inGS inData False False Nothing) candidateJoinedGraphList

                newMinCost = if (not . null) minEdgeList then minimum $ fmap snd6 rediagnosedGraphList
                             else infinity

                -- only taking one for SA/Drift check
                newMinGraph = if newMinCost /= infinity then head $ filter ((== newMinCost) . snd6) rediagnosedGraphList
                              else emptyPhylogeneticGraph
                
            in
            -- if better always return it--hope this conditions short circuits so don't fully diagnose graph all the time
            if minDelta < curBestCost && newMinCost < curBestCost then
               let (_, newSAParams) = U.simAnnealAccept inSimAnnealParams curBestCost newMinCost
               in
               -- trace ("TBR SA Drift better: " ++ (show $ driftChanges $ fromJust newSAParams))
               ([newMinGraph], newSAParams)

            -- check if hit step limit--more for SA than drift
            else if ((currentStep $ fromJust inSimAnnealParams) >= (numberSteps $ fromJust inSimAnnealParams)) || ((driftChanges $ fromJust inSimAnnealParams) >= (driftMaxChanges $ fromJust inSimAnnealParams)) then 
                 ([], inSimAnnealParams)

            else
               let (acceptGraph, newSAParams) = U.simAnnealAccept inSimAnnealParams curBestCost minDelta
                   -- banner = if newMinCost <  curBestCost then "TBR heur better"
                   --         else "TBR Accepted not better"
               in

               -- if accepted (better or random) then return with updated annealing/Drift parameters
               if acceptGraph then 
                  -- trace (banner ++ (show $ driftChanges $ fromJust newSAParams))
                  ([newMinGraph], newSAParams)

               -- rejected--recurse wirth updated params
               else 
                  -- trace ("TBR SA Drift Reject: " ++ (show $ driftChanges $ fromJust newSAParams))
                  tbrJoin steepest inGS inData splitGraph splitGraphSimple splitCost doIA prunedGraphRootIndex originalConnectionOfPruned charInfoVV curBestCost (drop numEdgesToExamine edgesInPrunedGraph) newSAParams targetEdge
            -- )

-- | rerootPrunedAndMakeGraph reroots the pruned graph component on the rerootEdge and joins to base gaph at target edge
rerootPrunedAndMakeGraph :: SimpleGraph -> LG.Node -> LG.Node -> LG.LEdge EdgeInfo -> LG.LEdge EdgeInfo -> SimpleGraph
rerootPrunedAndMakeGraph splitGraphSimple prunedGraphRootIndex originalConnectionOfPruned (u,v, _) rerootEdge =
   -- get edges to delete and edges to add 
   let (prunedEdgesToAdd, prunedEdgesToDelete) = getTBREdgeEditsSimple splitGraphSimple prunedGraphRootIndex rerootEdge
       
       -- edges to connect rerooted pruned component and base graph
       connectingEdges = [(u, originalConnectionOfPruned, 0.0),(originalConnectionOfPruned, v, 0.0),(originalConnectionOfPruned, prunedGraphRootIndex, 0.0)]
       
       tbrNewGraph = LG.insEdges (connectingEdges ++ prunedEdgesToAdd) $ LG.delEdges ([(u,v),(originalConnectionOfPruned, prunedGraphRootIndex)] ++ prunedEdgesToDelete) splitGraphSimple
   in
   tbrNewGraph

-- | getTBREdgeEditsSimple takes and edge and returns the list of edit to pruned subgraph
-- as a pair of edges to add and those to delete
-- since reroot edge is directed (e,v), edges away from v will have correct
-- orientation. Edges between 'e' and the root will have to be flipped
-- original root edges and reroort edge are deleted and new root and edge spanning orginal root created
-- delete original connection edge and creates a new one--like SPR 
-- returns ([add], [delete])
getTBREdgeEditsSimple :: (Show a) => SimpleGraph -> LG.Node -> LG.LEdge a -> ([LG.LEdge Double],[LG.Edge])
getTBREdgeEditsSimple inGraph prunedGraphRootIndex rerootEdge =
   --trace ("Getting TBR Edits for " ++ (show rerootEdge)) (
   let -- originalRootEdgeNodes = LG.descendants inGraph prunedGraphRootIndex
       originalRootEdges = LG.out inGraph prunedGraphRootIndex

       -- get path from new root edge fst vertex to orginal root and flip those edges
       -- since (u,v) is u -> v u "closer" to root
       closerToPrunedRootEdgeNode = (fst3 rerootEdge, fromJust $ LG.lab inGraph $ fst3 rerootEdge)
       (nodesInPath, edgesinPath) = LG.postOrderPathToNode inGraph closerToPrunedRootEdgeNode (prunedGraphRootIndex, fromJust $ LG.lab inGraph prunedGraphRootIndex)

       -- don't want original root edges to be flipped since later deleted
       edgesToFlip = edgesinPath L.\\ originalRootEdges
       flippedEdges = fmap LG.flipLEdge edgesToFlip

       -- new edges on new root position and spanning old root
       -- add in closer vertex to root to make sure direction of edge is correct
       newEdgeOnOldRoot = if (snd3 $ head originalRootEdges) `elem` ((fst3 rerootEdge) : (fmap fst nodesInPath)) then (snd3 $ head originalRootEdges, snd3 $ last originalRootEdges, 0.0)
                          else (snd3 $ last originalRootEdges, snd3 $ head originalRootEdges, 0.0)
       newRootEdges = [(prunedGraphRootIndex, fst3 rerootEdge, 0.0 ),(prunedGraphRootIndex, snd3 rerootEdge, 0.0)]
       
   in
   -- assumes we are not checking original root
   -- rerooted
   -- delete orignal root edges and rerootEdge
   -- add new root edges
   -- and new edge on old root--but need orientation
   -- flip edges from new root to old (delete and add list)
   {-
   trace ("\n\nIn Graph:\n"++ (LG.prettyIndices inGraph) ++ "\nTBR Edits: " ++ (show (LG.toEdge rerootEdge, prunedGraphRootIndex))
        ++ " NewEdgeOldRoot: " ++ (show $ LG.toEdge newEdgeOnOldRoot)
        ++ " New rootEdges: " ++ (show $ fmap LG.toEdge newRootEdges)
        )
   -}
   --   ++ "\nEdges to add: " ++ (show $ fmap LG.toEdge $ newEdgeOnOldRoot : (flippedEdges ++ newRootEdges)) ++ "\nEdges to delete: " ++ (show $ rerootEdge : (fmap LG.toEdge (edgesToFlip ++ originalRootEdges))))
   (newEdgeOnOldRoot : (flippedEdges ++ newRootEdges), LG.toEdge rerootEdge : (fmap LG.toEdge (edgesToFlip ++ originalRootEdges)))
   -- )

-- | reoptimizeSplitGraphFromVertex fully labels the component graph that is connected to the specified vertex
-- retuning that graph with 2 optimized components and their cost
-- both components goo through multi-traversal optimizations
-- doIA option to only do IA optimization as opposed to full thing--should be enormously faster--but yet more approximate
-- creates final for both base graph and priunned component due to rerooting non-concordance of preorder and post order assignments
-- terminology bse graph is the component with the original root, pruned that which has been removed form the original
-- graph to be readded to edge set
-- The function
      -- 1) optimizes two components seprately from their "root"
      -- 2) takes nodes and edges for each and cretes new graph
      -- 3) returns graph and summed cost of two components
      -- 4) adds in root and netPenalty factor estimates since net penalty can only be calculated on full graph
            -- part of this is turning off net penalty cost when optimizing base and pruned graph components
-- if doIA is TRUE then call function that onl;y optimizes the IA assignments on the "original graph" after split.
-- this keeps teh IA chracters in sync across the two graphs
reoptimizeSplitGraphFromVertex :: GlobalSettings
                          -> ProcessedData
                          -> Bool
                          -> VertexCost
                          -> DecoratedGraph
                          -> Int
                          -> Int
                          -> (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertex inGS inData doIA netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex =
   if doIA then
      -- only reoptimize the IA states for dynamic characters
      reoptimizeSplitGraphFromVertexIA inGS inData netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex
   else
      -- perform full optimizations of nodes
      -- these required for full optimization
      let nonExactCharacters = U.getNumberSequenceCharacters (thd3 inData)
          origGraph = inSplitGraph -- thd6 origPhyloGraph
          leafGraph = LG.extractLeafGraph origGraph
          calcBranchLengths = False

          -- create simple graph version of split for post order pass
          splitGraphSimple = GO.convertDecoratedToSimpleGraph inSplitGraph


          -- create optimized base graph
          -- False for staticIA
          (postOrderBaseGraph, _) = T.generalizedGraphPostOrderTraversal (inGS {graphFactor = NoNetworkPenalty}) nonExactCharacters inData leafGraph False (Just startVertex) splitGraphSimple


          fullBaseGraph = PRE.preOrderTreeTraversal (inGS {graphFactor = NoNetworkPenalty}) (finalAssignment inGS) False calcBranchLengths (nonExactCharacters > 0) startVertex True postOrderBaseGraph

          -- create fully optimized pruned graph.  Post order tehn preorder

          -- get root node of pruned graph--parent since that is the full pruned piece (keeping that node for addition to base graph and edge creation)
          startPrunedNode = (prunedSubGraphRootVertex, fromJust $ LG.lab origGraph prunedSubGraphRootVertex)
          startPrunedParentNode = head $ LG.labParents origGraph prunedSubGraphRootVertex
          startPrunedParentEdge = (fst startPrunedParentNode, prunedSubGraphRootVertex, dummyEdge)


          -- False for staticIA
          (postOrderPrunedGraph, _) = T.generalizedGraphPostOrderTraversal (inGS {graphFactor = NoNetworkPenalty}) nonExactCharacters inData leafGraph False (Just prunedSubGraphRootVertex) splitGraphSimple


          -- False for staticIA
          fullPrunedGraph = PRE.preOrderTreeTraversal (inGS {graphFactor = NoNetworkPenalty}) (finalAssignment inGS) False calcBranchLengths (nonExactCharacters > 0) prunedSubGraphRootVertex True postOrderPrunedGraph

          -- get root node of base graph
          startBaseNode = (startVertex, fromJust $ LG.lab (thd6 fullBaseGraph) startVertex)



          -- get nodes and edges in base and pruned graph (both PhylogeneticGrapgs so thd6)
          (baseGraphNonRootNodes, baseGraphEdges) = LG.nodesAndEdgesAfter (thd6 fullBaseGraph) [startBaseNode]

          (prunedGraphNonRootNodes, prunedGraphEdges) = if LG.isLeaf origGraph prunedSubGraphRootVertex then ([],[])
                                                        else LG.nodesAndEdgesAfter (thd6 fullPrunedGraph) [startPrunedNode]

          -- make fully optimized graph from base and split components
          fullSplitGraph = LG.mkGraph ([startBaseNode, startPrunedNode, startPrunedParentNode] ++  baseGraphNonRootNodes ++ prunedGraphNonRootNodes) (startPrunedParentEdge : (baseGraphEdges ++ prunedGraphEdges))

          -- cost of split graph to be later combined with re-addition delta for heuristic graph cost
          prunedCost = if LG.isLeaf origGraph prunedSubGraphRootVertex then 0
                       else snd6 fullPrunedGraph
          splitGraphCost = ((1.0 + netPenaltyFactor) * ((snd6 fullBaseGraph) + prunedCost))

      in

      {-
      trace ("Orig graph cost " ++ (show $ subGraphCost $ fromJust $ LG.lab origGraph startVertex) ++ " Base graph cost " ++ (show $ snd6 fullBaseGraph) ++ " pruned subgraph cost " ++ (show prunedCost) ++ " at node " ++ (show prunedSubGraphRootVertex) ++ " parent " ++ (show $ fst startPrunedParentNode)

         ++ "\nBaseGraphNodes\n" ++ (show $ L.sort  $ fmap fst baseGraphNonRootNodes) ++ "\nPruned nodes from root: " ++ "\n" ++ (show $ fmap fst $ startPrunedNode : prunedGraphNonRootNodes)
         ++ "\nSplit Graph\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph fullSplitGraph)
         ++ "\nOrig graph:\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph origGraph))
      -}
      --trace ("reoptimizeSplitGraphFromVertex: " ++ (show splitGraphCost))
      (fullSplitGraph, splitGraphCost)

-- | reoptimizeSplitGraphFromVertexTuple wrapper for reoptimizeSplitGraphFromVertex with last 3 args as tuple
reoptimizeSplitGraphFromVertexTuple :: GlobalSettings
                          -> ProcessedData
                          -> Bool
                          -> VertexCost
                          -> (DecoratedGraph, Int , Int)
                          -> (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertexTuple inGS inData doIA netPenaltyFactor (inSplitGraph, startVertex, prunedSubGraphRootVertex) =
   reoptimizeSplitGraphFromVertex inGS inData doIA netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex


-- | reoptimizeSplitGraphFromVertexIA performs operations of reoptimizeSplitGraphFromVertex for static charcaters
-- but dynamic characters--only update IA assignments and initialized from origPhylo graph (at leaves) to keep IA characters in sync
-- since all "static" only need single traversal post order pass
reoptimizeSplitGraphFromVertexIA :: GlobalSettings
                          -> ProcessedData
                          -> VertexCost
                          -> DecoratedGraph
                          -> Int
                          -> Int
                          -> (DecoratedGraph, VertexCost)
reoptimizeSplitGraphFromVertexIA inGS inData netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex =
   --if graphType inGS /= Tree then error "Networks not yet implemented in reoptimizeSplitGraphFromVertexIA"
   --else
      let   nonExactCharacters = U.getNumberSequenceCharacters (thd3 inData)
            origGraph = inSplitGraph -- thd6 origPhyloGraph

            -- create leaf graphs--but copy IA final to prelim
            leafGraph = GO.copyIAFinalToPrelim $ LG.extractLeafGraph origGraph
            calcBranchLengths = False

            -- create simple graph version of split for post order pass
            splitGraphSimple = GO.convertDecoratedToSimpleGraph inSplitGraph

            --Create base graph
            -- create postorder assignment--but only from single traversal
            -- True flag fior staticIA
            postOrderBaseGraph = T.postOrderTreeTraversal (inGS {graphFactor = NoNetworkPenalty}) inData leafGraph True (Just startVertex) splitGraphSimple
            baseGraphCost = snd6 postOrderBaseGraph

            -- True flag fior staticIA
            fullBaseGraph = PRE.preOrderTreeTraversal (inGS {graphFactor = NoNetworkPenalty}) (finalAssignment inGS) True calcBranchLengths (nonExactCharacters > 0) startVertex True postOrderBaseGraph

            {-
            localRootCost = if (rootCost inGS) == NoRootCost then 0.0
                              else if (rootCost inGS) == Wheeler2015Root then T.getW15RootCost inGS postOrderBaseGraph
                              else error ("Root cost type " ++ (show $ rootCost inGS) ++ " is not yet implemented")
            -}

            -- get root node of base graph
            startBaseNode = (startVertex, fromJust $ LG.lab (thd6 fullBaseGraph) startVertex)

            --Create pruned graph
            -- get root node of pruned graph--parent since that is the full pruned piece (keeping that node for addition to base graph and edge creation)
            startPrunedNode = GO.makeIAPrelimFromFinal (prunedSubGraphRootVertex, fromJust $ LG.lab origGraph prunedSubGraphRootVertex)
            startPrunedParentNode =  head $ LG.labParents origGraph prunedSubGraphRootVertex
            startPrunedParentEdge = (fst startPrunedParentNode, prunedSubGraphRootVertex, dummyEdge)


            -- True flag fior staticIA
            postOrderPrunedGraph =  T.postOrderTreeTraversal (inGS {graphFactor = NoNetworkPenalty}) inData leafGraph True (Just prunedSubGraphRootVertex) splitGraphSimple
            prunedGraphCost = snd6 postOrderPrunedGraph

            -- True flag fior staticIA
            fullPrunedGraph = PRE.preOrderTreeTraversal (inGS {graphFactor = NoNetworkPenalty}) (finalAssignment inGS) True calcBranchLengths (nonExactCharacters > 0) prunedSubGraphRootVertex True postOrderPrunedGraph

            -- get nodes and edges in base and pruned graph (both PhylogeneticGrapgs so thd6)
            (baseGraphNonRootNodes, baseGraphEdges) = LG.nodesAndEdgesAfter (thd6 fullBaseGraph) [startBaseNode]

            (prunedGraphNonRootNodes, prunedGraphEdges) = if LG.isLeaf origGraph prunedSubGraphRootVertex then ([],[])
                                                          else LG.nodesAndEdgesAfter (thd6 fullPrunedGraph) [startPrunedNode]

            -- make fully optimized graph from base and split components
            fullSplitGraph = LG.mkGraph ([startBaseNode, startPrunedNode, startPrunedParentNode] ++  baseGraphNonRootNodes ++ prunedGraphNonRootNodes) (startPrunedParentEdge : (baseGraphEdges ++ prunedGraphEdges))

            splitGraphCost = ((1.0 + netPenaltyFactor) * (baseGraphCost + prunedGraphCost))

      in
      -- remove when working
      -- trace ("ROGFVIA split costs:" ++ (show (baseGraphCost, prunedGraphCost, localRootCost)) ++ " -> " ++ (show splitGraphCost)) (
      if splitGraphCost == 0 then
         error ("Split costs:" ++ (show (baseGraphCost, prunedGraphCost)) ++ " -> " ++ (show splitGraphCost)
            ++ " Split graph simple:\n" ++ (LG.prettify splitGraphSimple)
            ++ "\nFull:\n" ++ (show inSplitGraph)
            ++ "\nOriginal Graph:\n" ++ (show origGraph))
      else (fullSplitGraph, splitGraphCost)
      -- )


