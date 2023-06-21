{- |
Module      :  Swap.hs
Description :  Module specifying graph swapping rearrangement functions
Copyright   :  (c) 2021-2023 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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
                    , getUnionRejoinEdgeList
                    ) where

import qualified Data.List                                     as L
import           Data.Maybe
import qualified Data.Vector                                   as V
import           GeneralUtilities
import qualified GraphOptimization.Medians                     as M
import qualified GraphOptimization.PostOrderSoftWiredFunctions as POSW
import qualified GraphOptimization.PreOrderFunctions           as PRE
import qualified GraphOptimization.Traversals                  as T
import qualified Graphs.GraphOperations                        as GO
import qualified ParallelUtilities                             as PU
import           Types.Types
import qualified Utilities.LocalGraph                          as LG
import           Utilities.Utilities                           as U
-- import           Debug.Trace


-- | swapSPRTBR performs SPR or TBR branch (edge) swapping on graphs
-- runs both SPR and TBR depending on argument since so much duplicated functionality
-- 'steepest' abandons swap graph and switces to found graph as soon as anyhting 'better'
-- is found. The alternative (all) examines the entire neighborhood and retuns the best result
-- the retuns is a list of better graphs and the number of swapping rounds were required to ge there
-- if joinType ==  JoinAll is specified a single round is performed--otherwise a union rounds
-- alternate between joinPruned and joinAll.  This to be rapid but complete.
-- joinType = JoinAll for annealing/drifting
swapSPRTBR  :: SwapParams
            -> GlobalSettings
            -> ProcessedData
            -> Int
            -> [PhylogeneticGraph]
            -> [([Int], Maybe SAParams, PhylogeneticGraph)]
            -> ([PhylogeneticGraph], Int)
swapSPRTBR swapParams inGS inData inCounter curBestGraphs inTripleList =
   --trace ("SWAPSPRTBR: " ++ swapType ++ " " ++ joinType ++ " " ++ (show $ snd6 $ thd3 $ head inTripleList)) $
   if null inTripleList then (curBestGraphs, inCounter)
   else
      let (randomIntListSwap, inSimAnnealParams, inGraph) = head inTripleList
      in
      -- dont swap if worse than current best
      -- if snd6 inGraph > (snd6 . head) curBestGraphs then
      --   swapSPRTBR swapType joinType atRandom inGS inData numToKeep maxMoveEdgeDist steepest alternate doIA returnMutated curBestGraphs inCounter (tail inTripleList)
      --else
        if (joinType swapParams) == JoinAll || isJust inSimAnnealParams then
            swapSPRTBR' (swapParams {joinType = JoinAll}) inGS inData inCounter (randomIntListSwap, inSimAnnealParams, inGraph)
        else
            -- join with union pruing first then followed by joinAll, but joinAlternate will return on better gaphs to return to join prune
            let (firstList, firstCounter) = swapSPRTBR' (swapParams {joinType = JoinPruned}) inGS inData inCounter (randomIntListSwap, inSimAnnealParams, inGraph)

                -- the + 5 is to allow for extra buffer room with input graph and multiple equally costly solutions, can help
                bestFirstList = GO.selectGraphs Best (keepNum swapParams) 0.0 (-1) (inGraph : firstList)
                
                ---change JoinAlternate for return to pruned union
                (bestAlternateList, alternateCounter) = 
                  if (joinType swapParams) == JoinAlternate then 
                     let (alternateList, alternateCounter') = swapSPRTBRList (swapParams {joinType = JoinAlternate}) inGS inData firstCounter bestFirstList $ zip3 (U.generateRandIntLists (length bestFirstList) (head $ drop 2 randomIntListSwap)) (U.generateUniqueRandList (length bestFirstList) inSimAnnealParams) bestFirstList
                     in
                     (GO.selectGraphs Best  (keepNum swapParams) 0.0 (-1) $ inGraph : (alternateList ++ bestFirstList), alternateCounter') 

               -- JoinPruned-don't alternate during search
                  else (bestFirstList, firstCounter)             
                
                -- recursive list version as opposed ot parMap version
                -- should reduce memory footprint at cost of less parallelism--but random replicates etc should take care of that
                (afterSecondList, afterSecondCounter) = swapSPRTBRList (swapParams {joinType = JoinAll}) inGS inData alternateCounter bestAlternateList $ zip3 (U.generateRandIntLists (length bestAlternateList) (head $ drop 2 randomIntListSwap)) (U.generateUniqueRandList (length bestAlternateList) inSimAnnealParams) bestAlternateList

                bestSecondList = GO.selectGraphs Best (keepNum swapParams) 0.0 (-1) afterSecondList

                --add final joinall if buffer full?

            in
            if (snd6 $ head bestSecondList) < (snd6 $ head curBestGraphs) then
               swapSPRTBR swapParams inGS inData (afterSecondCounter + inCounter) bestSecondList (tail inTripleList)
            else 
               swapSPRTBR swapParams inGS inData (afterSecondCounter + inCounter) curBestGraphs (tail inTripleList)
            -- (bestSecondList, afterSecondCounter + inCounter)
            

-- | swapSPRTBRList is a wrapper around swapSPRTBR' allowing for a list of graphs and a current best cost
-- reduce time of swap
swapSPRTBRList :: SwapParams
               -> GlobalSettings
               -> ProcessedData
               -> Int
               -> [PhylogeneticGraph]
               -> [([Int], Maybe SAParams, PhylogeneticGraph)]
               -> ([PhylogeneticGraph], Int)
swapSPRTBRList swapParams inGS inData inCounter curBestGraphs tripleList =
   if null tripleList then (curBestGraphs, inCounter)
   else
      let (randomIntListSwap, inSimAnnealParams, inGraph) = head tripleList
      in

      -- currrent graph worse than best saved
      if snd6 inGraph > (snd6 . head) curBestGraphs then
         swapSPRTBRList swapParams inGS inData inCounter curBestGraphs (tail tripleList)
      else
         let (graphList, swapCounter) = swapSPRTBR' swapParams inGS inData inCounter (randomIntListSwap, inSimAnnealParams, inGraph)
             bestNewGraphList = GO.selectGraphs Best (keepNum swapParams) 0.0 (-1) graphList
         in

         -- found better
         if (snd6 $ head bestNewGraphList) < (snd6 $ head curBestGraphs )then
               if (joinType swapParams == JoinAlternate || joinType swapParams == JoinPruned) then
                  let tripleToSwap = zip3 (U.generateRandIntLists (head $ drop (inCounter + 1) $ randomIntListSwap) (length bestNewGraphList)) (U.generateUniqueRandList (length bestNewGraphList) inSimAnnealParams) bestNewGraphList
                  in
                  swapSPRTBRList swapParams inGS inData swapCounter bestNewGraphList (tripleToSwap ++ tail tripleList)
               else 
                  swapSPRTBRList swapParams inGS inData swapCounter bestNewGraphList (tail tripleList)

         -- found worse
         else if (snd6 $ head bestNewGraphList) > (snd6 $ head curBestGraphs) then
               swapSPRTBRList swapParams inGS inData swapCounter curBestGraphs (tail tripleList)

         -- found equal
         else 
            let newUniqueBestGraphList = GO.selectGraphs Unique (keepNum swapParams) 0.0 (-1) (curBestGraphs ++ bestNewGraphList)
            in
            swapSPRTBRList swapParams inGS inData swapCounter newUniqueBestGraphList (tail tripleList)


-- | swapSPRTBR' is the central functionality of swapping allowing for repeated calls with alternate
-- options such as joinType to ensure complete swap but with an edge unions pass to
-- reduce time of swap
-- also manages SA/Drifting versus greedy swap
swapSPRTBR' :: SwapParams
            -> GlobalSettings
            -> ProcessedData
            -> Int
            -> ([Int], Maybe SAParams, PhylogeneticGraph)
            -> ([PhylogeneticGraph], Int)
swapSPRTBR' swapParams inGS inData inCounter (randomIntListSwap, inSimAnnealParams, inGraph) =
   -- trace ("SPRTBRList:" ++ (show $ joinType swapParams)) $
   if LG.isEmpty (fst6 inGraph) then ([], 0)
   else
      let numLeaves = V.length $ fst3 inData
          {-
          leafGraph = GO.makeSimpleLeafGraph inData
          leafDecGraph = GO.makeLeafGraph inData
          leafGraphSoftWired = POSW.makeLeafGraphSoftWired inGS inData
          -}

          inGraphNetPenalty = POSW.getNetPenalty inGS inData inGraph
          inGraphNetPenaltyFactor = inGraphNetPenalty / (snd6 inGraph)
      in
      -- trace ("SSPRTBR:" ++ (show inGraphNetPenaltyFactor)) (

      if inSimAnnealParams == Nothing then
        -- trace ("Non SA swap") (
        -- steepest takes immediate best--does not keep equall cost-- for now--disabled not working correctly so goes to "all"
        -- Nothing for SimAnneal Params
         let (swappedGraphs, counter, _) = swapAll swapParams inGS inData randomIntListSwap inCounter (snd6 inGraph) [] [inGraph] numLeaves inGraphNetPenaltyFactor inSimAnnealParams
         in
         
         if null swappedGraphs then ([inGraph], counter)
         else (swappedGraphs, counter)
         


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
             (annealDriftGraphs', anealDriftCounterList, _) = unzip3 $ (PU.seqParMap (parStrategy $ lazyParStrat inGS) (swapAll swapParams inGS inData randomIntListSwap 0 (snd6 inGraph) [] [inGraph] numLeaves inGraphNetPenaltyFactor) newSimAnnealParamList) -- `using` PU.myParListChunkRDS)

             -- annealed/Drifted 'mutated' graphs
             annealDriftGraphs = GO.selectGraphs Unique (keepNum swapParams) 0.0 (-1) $ concat annealDriftGraphs'

             -- swap back "normally" if desired for full drifting/annealing
             (swappedGraphs, counter, _) = swapAll swapParams inGS inData randomIntListSwap (sum anealDriftCounterList) (min (snd6 inGraph) (minimum $ fmap snd6 annealDriftGraphs)) [] annealDriftGraphs numLeaves inGraphNetPenaltyFactor Nothing

             bestGraphs = GO.selectGraphs Best (keepNum swapParams) 0.0 (-1) (inGraph : swappedGraphs)
         in
         -- trace ("Steepest SSPRTBR: " ++ (show (length swappedGraphs, counter)))
         --trace ("AC:" ++ (show $ fmap snd6 $ concat annealedGraphs') ++ " -> " ++ (show $ fmap snd6 $ swappedGraphs')) (

         -- this Bool for Genetic Algorithm mutation step
         if not (returnMutated swapParams) then (bestGraphs, counter)
         else (annealDriftGraphs, sum anealDriftCounterList)
         -- )
         -- )

-- | swapAll is a high level function that basically deals with portioning out swap-type swaps
-- and performs the high level options for Alternate where SPR is perfomred first, then TBR,
-- but whenever a better (or additional) graph is found during TBR, an SPR swap of that graph
-- is performed before returning to TBR again.  THis contibues untill no new graphs are found in the
-- SPR + TBR swap.
-- each call to swapAll' sets break edge number to 0
-- this swaps untill none found better
swapAll  :: SwapParams
         -> GlobalSettings
         -> ProcessedData
         -> [Int]
         -> Int
         -> VertexCost
         -> [PhylogeneticGraph]
         -> [PhylogeneticGraph]
         -> Int
         -> VertexCost
         -> Maybe SAParams
         -> ([PhylogeneticGraph], Int, Maybe SAParams)
swapAll swapParams inGS inData randomIntListSwap counter curBestCost curSameBetterList inGraphList numLeaves netPenaltyFactor inSimAnnealParams =
   --trace ("SWAPALL: " ++ swapType ++ " " ++ joinType ++ " " ++ (show curBestCost)) $
   if null inGraphList then (curSameBetterList, counter, inSimAnnealParams)
   else
      -- nni, spr, tbr
      if (swapType swapParams) /= Alternate then
         -- 0 is for list of edges so move through list past stable edges 
         swapAll' swapParams inGS inData randomIntListSwap counter curBestCost curSameBetterList inGraphList numLeaves netPenaltyFactor 0 inSimAnnealParams

      -- alternate
      else
         let -- spr first
             -- 0 is for list of edges so move through list past stable edges 
             (sprGraphs, sprCounter, sprSAPArams) = swapAll' (swapParams {swapType = SPR}) inGS inData randomIntListSwap counter curBestCost curSameBetterList inGraphList numLeaves netPenaltyFactor 0 inSimAnnealParams
             graphsToTBR = GO.selectGraphs Best (keepNum swapParams) 0.0 (-1) (sprGraphs ++ inGraphList)
             sprBestCost = (snd6 . head) graphsToTBR

             -- tbr until find better or novel equal
             (tbrGraphs, tbrCounter, tbrSAPArams) = swapAll' (swapParams {swapType = TBRAlternate})  inGS inData (tail randomIntListSwap) sprCounter sprBestCost graphsToTBR graphsToTBR numLeaves netPenaltyFactor 0 sprSAPArams
             tbrBestCost = if (not . null) tbrGraphs then
                              (snd6 . head) tbrGraphs
                           else infinity

         in
         {- This isn't improving performance so turned off in SwapSPRTBR-}
         -- if found better and alternating union pruning then return so can go back to start union pruning again
         if (joinType swapParams) == JoinAlternate && sprBestCost < curBestCost then (sprGraphs, sprCounter, sprSAPArams)

         else if (joinType swapParams) == JoinAlternate && tbrBestCost < curBestCost then (tbrGraphs, tbrCounter, tbrSAPArams)

         -- if TBR found better go around again with SPR first--since returned if found better during TBR rejoin
         else if tbrBestCost < sprBestCost then
            swapAll swapParams inGS inData (drop 2 randomIntListSwap) tbrCounter tbrBestCost tbrGraphs tbrGraphs numLeaves netPenaltyFactor tbrSAPArams

         -- check if found additional
         else if tbrBestCost == sprBestCost then
            let bestTBRGraphs = GO.selectGraphs Best (keepNum swapParams) 0.0 (-1) tbrGraphs
                newTBRGraphs = bestTBRGraphs `GO.phylogeneticGraphListMinus` sprGraphs 
            in
            if (not . null) newTBRGraphs && (length bestTBRGraphs < keepNum swapParams) then
               swapAll swapParams inGS inData (drop 3 randomIntListSwap) tbrCounter tbrBestCost tbrGraphs newTBRGraphs numLeaves netPenaltyFactor tbrSAPArams

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
   -- assumes SPR done before  Alternate entering so can star with TBR and iff get better
   -- go back to SPR. NBest for "steepest" descent
-- For drift and anneal need to randomize order of splits and rejoins
swapAll' :: SwapParams
         -> GlobalSettings
         -> ProcessedData
         -> [Int]
         -> Int
         -> VertexCost
         -> [PhylogeneticGraph]
         -> [PhylogeneticGraph]
         -> Int
         -> VertexCost
         -> Int
         -> Maybe SAParams
         -> ([PhylogeneticGraph], Int, Maybe SAParams)
swapAll' swapParams inGS inData randomIntListSwap counter curBestCost curSameBetterList inGraphList numLeaves netPenaltyFactor breakEdgeNumber inSimAnnealParams =
   -- trace (" In cost " ++ (show curBestCost) ++ (" " ++ swapType)) (
   -- don't beed to check for mutated here since checked above
   if null inGraphList then
      -- trace (" Out cost " ++ (show curBestCost) ++ (" " ++ swapType))
      (GO.selectGraphs Unique (keepNum swapParams) 0.0 (-1) curSameBetterList, counter, inSimAnnealParams)
   else
      let firstGraph = head inGraphList
          firstDecoratedGraph = thd6 firstGraph
          (firstRootIndex, _) = head $ LG.getRoots firstDecoratedGraph

          -- determine edges to break on--'bridge' edges only for network
          -- filter out edges from root since no use--would just rejoin
          -- sort longest edge to shortest--option to speeed up steepest and conditions for all as well
          -- this edge sort from Varon and Wheeler 2013
          breakEdgeList' = if (graphType inGS) == Tree || LG.isTree firstDecoratedGraph then
                              if not (atRandom swapParams) then GO.sortEdgeListByLength $ filter ((/= firstRootIndex) . fst3) $ LG.labEdges firstDecoratedGraph
                              else filter ((/= firstRootIndex) . fst3) $ LG.labEdges firstDecoratedGraph
                          else
                              if not (atRandom swapParams) then GO.sortEdgeListByLength $ filter ((/= firstRootIndex) . fst3) $ LG.getEdgeSplitList firstDecoratedGraph
                              else filter ((/= firstRootIndex) . fst3) $ LG.getEdgeSplitList firstDecoratedGraph

          -- randomize edges list order for anneal and drift
          breakEdgeList'' = if isJust inSimAnnealParams then
                              permuteList (head $ randomIntegerList $ fromJust inSimAnnealParams) breakEdgeList'
                            else if (atRandom swapParams) then
                              permuteList (head randomIntListSwap) breakEdgeList'
                            else breakEdgeList'

          -- move first "breakEdgeFactor" edges in split list to end
          -- since breakEdgeFactor can get incremented past number of edges the integer remainder is determined
          -- to move to end
          -- this to reduces the revisiting of stable edges (by moving them to the end of the list)
          -- yet still insures that all edges will be visited in final (or ay time needed) split.
          -- used in POY v 1-3, Came from Steve Farris pers. com.
          breakEdgeFactor = snd $ divMod breakEdgeNumber (length breakEdgeList'')
          breakEdgeList =  (drop breakEdgeFactor breakEdgeList'') ++ (take breakEdgeFactor breakEdgeList'')


          -- perform intial split and rejoin on each edge in first graph
          (newGraphList', newSAParams, newBreakEdgeNumber) = splitJoinGraph swapParams inGS inData (tail randomIntListSwap) curBestCost curSameBetterList numLeaves netPenaltyFactor inSimAnnealParams firstGraph breakEdgeNumber breakEdgeList breakEdgeList

          -- get best return graph list-can be empty if nothing better ort smame cost
          newGraphList = GO.selectGraphs Best (keepNum swapParams) 0.0 (-1) newGraphList'

          -- get unique return graph list-can be empty if nothing better ort same cost
          -- newGraphListUnique = GO.selectGraphs Unique (maxBound::Int) 0.0 (-1) newGraphList'

          newMinCost = if (not . null) newGraphList' then snd6 $ head newGraphList
                       else infinity

      in
      -- logic for returning normal swap operations (better only)
      -- versus simulated annealin/Drifing returning potentially sub-optimal
      if isNothing inSimAnnealParams then
         -- postProcessSwap swapParams inGS inData (drop 2 randomIntListSwap) counter curBestCost curSameBetterList inGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired netPenaltyFactor Nothing newMinCost newBreakEdgeNumber newGraphList

         -- found better cost graph
         if newMinCost < curBestCost then
            traceNoLF ("\t->" ++ (show newMinCost)) $ --  ++ " " ++ (show curBestCost)) $
            -- for alternarte do SPR first then TBR
            
            -- for alternate in TBR or prune union alternate if found better return immediately
            if (swapType swapParams == TBRAlternate) || (joinType swapParams == JoinAlternate) then (newGraphList, counter, newSAParams)

            -- regular swap--keep going with better graphs
            else swapAll' swapParams inGS inData randomIntListSwap (counter + 1) newMinCost newGraphList newGraphList numLeaves netPenaltyFactor newBreakEdgeNumber newSAParams
            

         -- found only worse graphs--never happens due to the way splitjoin returns only better or equal
         -- but could change
         else if newMinCost > curBestCost then
            -- trace ("Worse " ++ (show newMinCost)) (
            -- breakEdgeNUmber set to zero for new graph to look at
            swapAll' swapParams inGS inData randomIntListSwap (counter + 1) curBestCost curSameBetterList (tail inGraphList) numLeaves netPenaltyFactor 0 newSAParams
            -- )

         -- found same cost graphs
         else
            -- Important to not limit curSameBest since may rediscover graphs via swapping on equal when limiting the number to keep
            -- can be a cause of infinite running issues.
            let newCurSameBetterList = GO.selectGraphs Best (keepNum swapParams) 0.0 (-1) (curSameBetterList ++ newGraphList)

                graphsToDo  = GO.selectGraphs Best (keepNum swapParams) 0.0 (-1)  $ ((tail inGraphList) ++ newGraphList) `GO.phylogeneticGraphListMinus` curSameBetterList

                -- these conditions help to prevent recswapping endlessly on new graphs thatare not in buffers,
                -- but have same cost
                graphsToDo' = if length graphsToDo >= ((keepNum swapParams) - 1) then (tail inGraphList)

                              -- found nothing better that is new
                              else if length newCurSameBetterList >= ((keepNum swapParams) - 1) then (tail inGraphList)

                              else graphsToDo

            in
            swapAll' swapParams inGS inData (tail randomIntListSwap) (counter + 1) curBestCost newCurSameBetterList graphsToDo' numLeaves netPenaltyFactor newBreakEdgeNumber newSAParams

      -- simulated annealing/Drift post processing
      else
         -- postProcessAnnealDrift swapParams inGS inData (drop 2 randomIntListSwap) counter curBestCost curSameBetterList inGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired netPenaltyFactor newSAParams newMinCost newBreakEdgeNumber newGraphListUnique

         -- found better cost graph
         if newMinCost < curBestCost then
            traceNoLF ("\t->" ++ (show newMinCost))
            (newGraphList, counter, newSAParams)

         -- not better so check for drift changes or annealing steps and return if reached maximum number
         else if ((currentStep $ fromJust inSimAnnealParams) >= (numberSteps $ fromJust inSimAnnealParams)) || ((driftChanges $ fromJust inSimAnnealParams) >= (driftMaxChanges $ fromJust inSimAnnealParams)) then
            --trace ("PPA return: " ++ (show (newMinCost, curBestCost)))
            (GO.selectGraphs Unique (keepNum swapParams) 0.0 (-1) (newGraphList ++ curSameBetterList), counter, newSAParams)

         -- didn't hit stopping numbers so continuing--but based on current best cost not whatever was found
         else
            -- traceNoLF ("[" ++ (show $ length (newGraphList ++ (tail inGraphList))) ++ "]")
            swapAll' swapParams inGS inData randomIntListSwap (counter + 1) curBestCost (newGraphList ++ curSameBetterList) (newGraphList ++ (tail inGraphList))  numLeaves netPenaltyFactor breakEdgeNumber newSAParams

      --)

{-
-- | postProcessSwap factors out the post processing of swap results to allow for clearer code
-- with "regular" optimal swapping
postProcessSwap   :: SwapParams
                  -> GlobalSettings
                  -> ProcessedData
                  -> [Int]
                  -> Int
                  -> VertexCost
                  -> [PhylogeneticGraph]
                  -> [PhylogeneticGraph]
                  -> Int
                  -> SimpleGraph
                  -> DecoratedGraph
                  -> DecoratedGraph
                  -> VertexCost
                  -> Maybe SAParams
                  -> VertexCost
                  -> Int
                  -> [PhylogeneticGraph]
                  -> ([PhylogeneticGraph], Int, Maybe SAParams)
postProcessSwap swapParams inGS inData randomIntListSwap counter curBestCost curSameBetterList inGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired netPenaltyFactor inSimAnnealParams newMinCost breakEdgeNumber newGraphList =

      -- found better cost graph
      if newMinCost < curBestCost then
         traceNoLF ("\t->" ++ (show newMinCost)) $ --  ++ " " ++ (show curBestCost)) $
         -- for alternarte do SPR first then TBR
         let graphsToSwap = GO.selectGraphs Best (keepNum swapParams) 0.0 (-1) newGraphList -- (newGraphList ++ (tail inGraphList))
         in

         -- for alternate in TBR or prune union alternate if found better return immediately
         if (swapType swapParams == TBRAlternate) || (joinType swapParams == JoinAlternate) then (newGraphList, counter, inSimAnnealParams)

         -- regular swap--keep going with better graphs
         else swapAll' swapParams inGS inData randomIntListSwap (counter + 1) newMinCost newGraphList graphsToSwap numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired netPenaltyFactor breakEdgeNumber inSimAnnealParams
         

      -- found only worse graphs--never happens due to the way splitjoin returns only better or equal
      -- but could change
      else if newMinCost > curBestCost then
         -- trace ("Worse " ++ (show newMinCost)) (
         -- breakEdgeNUmber set to zero for new graph to look at
         swapAll' swapParams inGS inData randomIntListSwap (counter + 1) curBestCost curSameBetterList (tail inGraphList) numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired netPenaltyFactor 0 inSimAnnealParams
         -- )

      -- found same cost graphs
      else
         -- Important to not limit curSameBest since may rediscover graphs via swapping on equal when limiting the number to keep
         -- can be a cause of infinite running issues.
         let newCurSameBetterList = GO.selectGraphs Best (keepNum swapParams) 0.0 (-1) (curSameBetterList ++ newGraphList)

             graphsToDo  = GO.selectGraphs Best (keepNum swapParams) 0.0 (-1)  $ ((tail inGraphList) ++ newGraphList) `GO.phylogeneticGraphListMinus` curSameBetterList

             -- these conditions help to prevent recswapping endlessly on new graphs thatare not in buffers,
             -- but have same cost
             graphsToDo' = if length graphsToDo >= ((keepNum swapParams) - 1) then (tail inGraphList)
                           else if length newCurSameBetterList == length curSameBetterList then (tail inGraphList)
                           else if length newCurSameBetterList >= ((keepNum swapParams) - 1) then (tail inGraphList)
                           --else if (fmap fst6 newCurSameBetterList) ==  (fmap fst6 curSameBetterList) then (tail inGraphList)
                           else graphsToDo

         in
         swapAll' swapParams inGS inData (tail randomIntListSwap) (counter + 1) curBestCost newCurSameBetterList graphsToDo' numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired netPenaltyFactor breakEdgeNumber inSimAnnealParams
      

-- | postProcessAnnealDrift factors out the post processing of swap results to allow for clearer code
-- with simulated annealing/Drifting returns of potentially
-- non optimal graphs
-- if better gaph cost found then update graph and better cost
-- else if graphs probablitstically found--return them if hit max drfit changes or annealing steps
-- else go on with updated number of changes/steps
-- removes alternate and just does straight NNI/SPR/TBR in order to properly do teh returns of altered and better graphs
-- not sure about curBestGraph and return--new Graphs could be empty
postProcessAnnealDrift  :: SwapParams
                        -> GlobalSettings
                        -> ProcessedData
                        -> [Int]
                        -> Int
                        -> VertexCost
                        -> [PhylogeneticGraph]
                        -> [PhylogeneticGraph]
                        -> Int
                        -> SimpleGraph
                        -> DecoratedGraph
                        -> DecoratedGraph
                        -> VertexCost
                        -> Maybe SAParams
                        -> VertexCost
                        -> Int
                        -> [PhylogeneticGraph]
                        -> ([PhylogeneticGraph], Int, Maybe SAParams)
postProcessAnnealDrift swapParams inGS inData randomIntListSwap counter curBestCost curSameBetterList inGraphList numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired netPenaltyFactor inSimAnnealParams newMinCost breakEdgeNumber newGraphList =
   -- trace ("PPA: " ++ (show (newMinCost, curBestCost)) ++ " counts " ++ (show ((currentStep $ fromJust inSimAnnealParams), (numberSteps $ fromJust inSimAnnealParams), (driftChanges $ fromJust inSimAnnealParams), (driftMaxChanges $ fromJust inSimAnnealParams)))) (
      -- found better cost graph
      if newMinCost < curBestCost then
         traceNoLF ("\t->" ++ (show newMinCost))
         (newGraphList, counter, inSimAnnealParams)

      -- not better so check for drift changes or annealing steps and return if reached maximum number
      else if ((currentStep $ fromJust inSimAnnealParams) >= (numberSteps $ fromJust inSimAnnealParams)) || ((driftChanges $ fromJust inSimAnnealParams) >= (driftMaxChanges $ fromJust inSimAnnealParams)) then
         --trace ("PPA return: " ++ (show (newMinCost, curBestCost)))
         (GO.selectGraphs Unique (keepNum swapParams) 0.0 (-1) (newGraphList ++ curSameBetterList), counter, inSimAnnealParams)

      -- didn't hit stopping numbers so continuing--but based on current best cost not whatever was found
      else
         -- traceNoLF ("[" ++ (show $ length (newGraphList ++ (tail inGraphList))) ++ "]")
         swapAll' swapParams inGS inData randomIntListSwap (counter + 1) curBestCost (newGraphList ++ curSameBetterList) (newGraphList ++ (tail inGraphList))  numLeaves leafSimpleGraph leafDecGraph leafGraphSoftWired netPenaltyFactor breakEdgeNumber inSimAnnealParams

      -- )
-}

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
splitJoinGraph :: SwapParams
               -> GlobalSettings
               -> ProcessedData
               -> [Int]
               -> VertexCost
               -> [PhylogeneticGraph]
               -> Int
               -> VertexCost
               -> Maybe SAParams
               -> PhylogeneticGraph
               -> Int
               -> [LG.LEdge EdgeInfo]
               -> [LG.LEdge EdgeInfo]
               -> ([PhylogeneticGraph], Maybe SAParams, Int)
splitJoinGraph swapParams inGS inData randomIntListSwap curBestCost curSameBetterList numLeaves netPenaltyFactor inSimAnnealParams firstGraph breakEdgeNumber' breakEdgeListComplete breakEdgeList =
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
          (reoptimizedSplitGraph, splitCost) = reoptimizeSplitGraphFromVertex inGS inData (doIA swapParams) netPenaltyFactor splitGraph graphRoot prunedGraphRootIndex

          -- get root in base (for readdition) and edges in pruned section for rerooting during TBR readdition
          (_, edgesInPrunedGraph) = LG.nodesAndEdgesAfter splitGraph [(originalConnectionOfPruned, fromJust $ LG.lab splitGraph originalConnectionOfPruned)]

          edgesInBaseGraph = breakEdgeListComplete L.\\ (edgeToBreakOn : edgesInPrunedGraph)

          -- insert here union calcuations based on Varon and Wheeler (2013)
          -- basically--rebuild edge to rejoin list based on critical value, totalCost - splitCost, and
          -- edge union distance
          -- build edges pre-order and add to rejoin list if
               -- 1) not network (but still recurse to children)
               -- 2) union delta below threshold and recurse to children
             -- if > threshold then stop, no add, no recurse since children can only get hihger ubnion distance
          -- use split graph (with reoptimized nodes) and overall graph root to get avialbel edges in base graph for rejoin

          prunedToRejoinUnionData = vertData $ fromJust $ LG.lab reoptimizedSplitGraph prunedGraphRootIndex
          --prunedToRejoinUnionData = vertData $ fromJust $ LG.lab (thd6 firstGraph) prunedGraphRootIndex
          charInfoVV = six6 firstGraph
          unionEdgeList = getUnionRejoinEdgeList inGS reoptimizedSplitGraph charInfoVV [graphRoot] ((snd6 firstGraph) - splitCost) prunedToRejoinUnionData [] 

          -- builds graph edge list with unions--need to be able to turn off and just used edges in base graph for some sort
          -- of "no-union" swap
          -- determine those edges within distance of original if limited (ie NNI etc)
          rejoinEdges' = if (maxMoveEdgeDist swapParams) >= ((maxBound :: Int) `div` 3) then
                           if (joinType swapParams) == JoinAll then edgesInBaseGraph
                           else unionEdgeList
                         else
                           let candidateEdges = take (maxMoveEdgeDist swapParams) $ (LG.sortEdgeListByDistance splitGraph [graphRoot] [graphRoot])
                           in
                           if (joinType swapParams) == JoinAll then candidateEdges
                           else L.intersect candidateEdges unionEdgeList


          -- randomize edges list order for anneal and drift
          rejoinEdges = if isJust inSimAnnealParams then
                           permuteList ((randomIntegerList $ fromJust inSimAnnealParams) !! 1) rejoinEdges'
                        else if (atRandom swapParams) then
                              permuteList (head randomIntListSwap) rejoinEdges'
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
            rejoinGraph swapParams inGS inData curBestCost [] netPenaltyFactor reoptimizedSplitGraph (GO.convertDecoratedToSimpleGraph splitGraph) splitCost graphRoot prunedGraphRootIndex originalConnectionOfPruned rejoinEdges edgesInPrunedGraph inSimAnnealParams

          newGraphList' = GO.selectGraphs Best (keepNum swapParams) 0.0 (-1) newGraphList
      in
      --trace ("SJG:" ++ (show (length edgesInBaseGraph, length unionEdgeList))) $
      -- regular swap
      if isNothing inSimAnnealParams then
         -- only returns graphs if same of better else empty
         -- adds null o\r better graphs to reurn list
         if (not . null) newGraphList && (steepest swapParams) then (newGraphList', inSimAnnealParams, breakEdgeNumber)
         else
            let (recurseGraphList, _, newEdgeBreakNumber) = splitJoinGraph swapParams inGS inData (tail randomIntListSwap) curBestCost curSameBetterList numLeaves netPenaltyFactor inSimAnnealParams firstGraph breakEdgeNumber breakEdgeListComplete (tail breakEdgeList)
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
            splitJoinGraph swapParams inGS inData (tail randomIntListSwap) curBestCost curSameBetterList numLeaves netPenaltyFactor newSAParams firstGraph breakEdgeNumber breakEdgeListComplete (tail breakEdgeList)
      --)

-- | getUnionRejoinEdgeList takes a graph (split and reoptimized usually), the overall root index (of split),
-- split cost, and union threshold value and returns list of edges that have union distance <= threshold factor
-- assumes that root edges are not network edges (an invariant)
-- checks node then recurses to children
getUnionRejoinEdgeList :: GlobalSettings 
                       -> DecoratedGraph 
                       -> V.Vector (V.Vector CharInfo)
                       -> [LG.Node] 
                       -> Double 
                       -> VertexBlockData 
                       -> [LG.LEdge EdgeInfo] 
                       -> [LG.LEdge EdgeInfo]
getUnionRejoinEdgeList inGS inGraph charInfoVV nodeIndexList splitDiffCost nodeToJoinUnionData curEdgeList =
   -- returns edges in post order since prepending--could reverse if want pre-order edges
   -- might be better post order, unclear
   if null nodeIndexList then curEdgeList
   else
      let nodeIndex = head nodeIndexList
          childEdges = LG.out inGraph nodeIndex
          childNodeIndexList = fmap snd3 childEdges
       
          -- node data and union distance
          nodeData = vertData $ fromJust $ LG.lab inGraph nodeIndex
          unionDistance = getUnionDistance nodeToJoinUnionData nodeData charInfoVV
          metThreshold = unionDistance < splitDiffCost * (unionThreshold inGS)

          -- nodeDataString = U.getUnionFieldsNode nodeData
          -- toJoinString = U.getUnionFieldsNode nodeToJoinUnionData
      in
      -- traceNoLF ("GURE: " ++ (show nodeIndex)) (
      -- trace ("GUREL: " ++ (show (unionDistance, splitDiffCost, unionThreshold * splitDiffCost))) (
      if (length childEdges) `notElem` [0,1,2] then error ("Node has improper number of children : " ++ (show $ length childEdges))
      -- add non-outgroup root edge to list for rejoin after checking for acceptable union distance
      -- if edge is not within union distance factor then stop--no recursion
      -- this since distance cannot get lower further upt the graph given union creation
      else if LG.isRoot inGraph nodeIndex then
               if metThreshold then
                  if (null $ LG.out inGraph (head childNodeIndexList)) then
                     getUnionRejoinEdgeList inGS inGraph charInfoVV [(snd3 $ last childEdges)] splitDiffCost nodeToJoinUnionData ((last childEdges) : curEdgeList)
                  else
                     getUnionRejoinEdgeList inGS inGraph charInfoVV [(snd3 $ head childEdges)] splitDiffCost nodeToJoinUnionData ((head childEdges) : curEdgeList)
               else curEdgeList

      -- non-root node--process childre 1/2
      -- recurses to their children if union condition met--but doesn't add network edges
      -- check current node--then recurse to children
      else
         let -- first and second child data child
             newCurEdgeListChild = if metThreshold then
                                       getUnionRejoinEdgeList inGS inGraph charInfoVV childNodeIndexList splitDiffCost nodeToJoinUnionData (childEdges ++  curEdgeList) 

                                    else curEdgeList


         in
         -- recurse remaining nodes
         getUnionRejoinEdgeList inGS inGraph charInfoVV(tail nodeIndexList) splitDiffCost nodeToJoinUnionData newCurEdgeListChild
         -- )

-- | getUnionDistance gets distance between two the union fields of two characters
getUnionDistance :: VertexBlockData -> VertexBlockData -> V.Vector (V.Vector CharInfo) -> Double
getUnionDistance union1 union2 charInfoVV =
   let (_, newUnionCost) = M.distance2Unions union1 union2 charInfoVV

   in
   newUnionCost

-- | rejoinGraphTuple is a wrapper around rejoinGraph for fmapping--only returns graph list not simulated annealing params
rejoinGraphTuple :: SwapParams 
                 -> GlobalSettings
                 -> ProcessedData
                 -> VertexCost
                 -> [PhylogeneticGraph]
                 -> Maybe SAParams
                 -> (DecoratedGraph, SimpleGraph, VertexCost, LG.Node,LG.Node, LG.Node, [LG.LEdge EdgeInfo], [LG.LEdge EdgeInfo], VertexCost)
                 -> [PhylogeneticGraph]
rejoinGraphTuple swapParams inGS inData curBestCost curBestGraphs inSimAnnealParams (reoptimizedSplitGraph, splitGraphSimple, splitGraphCost, graphRoot, prunedGraphRootIndex, originalConnectionOfPruned, rejoinEdges, edgesInPrunedGraph, netPenaltyFactor) =
   fst $ rejoinGraph swapParams inGS inData curBestCost curBestGraphs netPenaltyFactor reoptimizedSplitGraph splitGraphSimple splitGraphCost graphRoot prunedGraphRootIndex originalConnectionOfPruned rejoinEdges edgesInPrunedGraph inSimAnnealParams

-- | rejoinGraph rejoins a split graph at all edges (if not steepest and found better)
-- in "base" graph.
-- if not steepest then do all as map, else recursive on base graph edge list
-- nni doesn't apper to be correct here--maybe loose it--doing nothing
rejoinGraph :: SwapParams 
            -> GlobalSettings
            -> ProcessedData
            -> VertexCost
            -> [PhylogeneticGraph]
            -> VertexCost
            -> DecoratedGraph
            -> SimpleGraph
            -> VertexCost
            -> LG.Node
            -> LG.Node
            -> LG.Node
            -> [LG.LEdge EdgeInfo]
            -> [LG.LEdge EdgeInfo]
            -> Maybe SAParams
            -> ([PhylogeneticGraph], Maybe SAParams)
rejoinGraph swapParams inGS inData curBestCost curBestGraphs netPenaltyFactor reoptimizedSplitGraph splitGraphSimple splitGraphCost graphRoot prunedGraphRootIndex originalConnectionOfPruned rejoinEdges' edgesInPrunedGraph inSimAnnealParams =

   -- found no better--but return equal cost graphs
   -- trace ("In rejoinGraph with num rejoining edges: " ++ (show $ length rejoinEdges')) (
   if null rejoinEdges' then (curBestGraphs, inSimAnnealParams)

   else
      -- this is for no  swapping option in fuse and genetic algorithm-fuse
      let rejoinEdges = if (swapType swapParams) == None then take 6 rejoinEdges'
                        else rejoinEdges'
      in
      -- regular swapping
      if isNothing inSimAnnealParams then

         -- check if split graph cost same as graph then return since graph can only get longer on readdition
         if splitGraphCost >= curBestCost then ([], inSimAnnealParams)

         else
            -- fmap over all edges in base graph
            if not (steepest swapParams) then
               let -- rejoinGraphList = concatMap (singleJoin swapType steepest inGS inData reoptimizedSplitGraph splitGraphSimple splitGraphCost doIA prunedGraphRootIndex originalConnectionOfPruned charInfoVV curBestCost edgesInPrunedGraph) rejoinEdges `using` PU.myParListChunkRDS
                   rejoinGraphList = concat $ fmap fst $ PU.seqParMap (parStrategy $ lazyParStrat inGS) (singleJoin swapParams inGS inData reoptimizedSplitGraph splitGraphSimple splitGraphCost prunedGraphRootIndex originalConnectionOfPruned curBestCost edgesInPrunedGraph inSimAnnealParams) rejoinEdges

                   {-Checking only min but seems to make slower
                   newMinCost = if null rejoinGraphList then infinity
                                else minimum $ fmap snd rejoinGraphList
                   (minEstCostNewGraphList, _) = unzip $ filter ((== newMinCost) . snd) rejoinGraphList
                   -}

                   -- newGraphList = fmap (T.multiTraverseFullyLabelGraph inGS inData False False Nothing) (fmap fst rejoinGraphList) `using` PU.myParListChunkRDS
                   newGraphList' = GO.selectGraphs Best (keepNum swapParams) 0.0 (-1) rejoinGraphList -- newGraphList
               in
               -- will only return graph if <= curBest cost
               if null rejoinGraphList then ([], inSimAnnealParams)
               else if (snd6 . head) newGraphList' <= curBestCost then (newGraphList', inSimAnnealParams)
               else ([], inSimAnnealParams)

            -- famp over number of threads edges in base graph
            -- then recurse
            else
               -- trace ("In steepest: " ++ (show PU.getNumThreads) ++ " " ++ (show $ length $ take PU.getNumThreads rejoinEdges)) (
               let -- this could be made a little parallel--but if lots of threads basically can do all
                   -- to not overload paralle threads
                   {-  This not so efficient is swapping in single graphs so leaving it be
                   saRounds = if isNothing inSimAnnealParams then 1
                              else rounds $ fromJust inSimAnnealParams

                   (numGraphsToExamine, _) = divMod PU.getNumThreads saRounds -- this may not "drift" if finds alot better, but that's how its supposed to work
                   -}
                   numGraphsToExamine = min (graphsSteepest inGS) PU.getNumThreads
                   rejoinEdgeList = take numGraphsToExamine rejoinEdges
                   --rejoinGraphList = concatMap (singleJoin swapType steepest inGS inData reoptimizedSplitGraph splitGraphSimple splitGraphCost doIA prunedGraphRootIndex originalConnectionOfPruned charInfoVV curBestCost edgesInPrunedGraph) rejoinEdgeList `using` PU.myParListChunkRDS
                   rejoinGraphList = concat $ fmap fst $ PU.seqParMap (parStrategy $ lazyParStrat inGS) (singleJoin swapParams inGS inData reoptimizedSplitGraph splitGraphSimple splitGraphCost prunedGraphRootIndex originalConnectionOfPruned curBestCost edgesInPrunedGraph inSimAnnealParams) rejoinEdgeList

                   {--Checking only min but seems to make slower
                   newMinCost = if null rejoinGraphList then infinity
                                else minimum $ fmap snd rejoinGraphList
                   (minEstCostNewGraphList, _) = unzip $ filter ((== newMinCost) . snd) rejoinGraphList
                   -}

                   -- newGraphList = fmap (T.multiTraverseFullyLabelGraph inGS inData False False Nothing) (fmap fst rejoinGraphList) `using` PU.myParListChunkRDS
                   newGraphList' = GO.selectGraphs Best (keepNum swapParams) 0.0 (-1) rejoinGraphList -- newGraphList
               in
               -- found nothing better or equal
               if null rejoinGraphList then
                  -- trace ("In steepest worse: " ++ (show $ length (drop PU.getNumThreads rejoinEdges)))
                   rejoinGraph swapParams inGS inData curBestCost curBestGraphs netPenaltyFactor reoptimizedSplitGraph splitGraphSimple splitGraphCost graphRoot prunedGraphRootIndex originalConnectionOfPruned (drop numGraphsToExamine rejoinEdges) edgesInPrunedGraph inSimAnnealParams

               -- found better graph
               else if (snd6 . head) newGraphList'  < curBestCost then
                  -- trace ("Steepest better")
                  (newGraphList', inSimAnnealParams)

               -- found equal cost graph
               else if (snd6 . head) newGraphList' == curBestCost then
                  let newBestList = GO.selectGraphs Best (keepNum swapParams) 0.0 (-1) (curBestGraphs ++ newGraphList')
                  in
                  rejoinGraph swapParams inGS inData curBestCost newBestList netPenaltyFactor reoptimizedSplitGraph splitGraphSimple splitGraphCost graphRoot prunedGraphRootIndex originalConnectionOfPruned (drop numGraphsToExamine rejoinEdges) edgesInPrunedGraph inSimAnnealParams
               -- found worse graphs only
               else
                  -- trace ("In steepest worse (after recalculation): " ++ (show $ length (drop PU.getNumThreads rejoinEdges)))
                  rejoinGraph swapParams inGS inData curBestCost curBestGraphs netPenaltyFactor reoptimizedSplitGraph splitGraphSimple splitGraphCost graphRoot prunedGraphRootIndex originalConnectionOfPruned (drop numGraphsToExamine rejoinEdges) edgesInPrunedGraph inSimAnnealParams
               -- ) -- )

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
             rejoinGraphPairList = PU.seqParMap (parStrategy $ lazyParStrat inGS) (singleJoin' swapParams inGS inData reoptimizedSplitGraph splitGraphSimple splitGraphCost prunedGraphRootIndex originalConnectionOfPruned curBestCost edgesInPrunedGraph) (zip simAnnealParamList rejoinEdgeList)

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
               rejoinGraph swapParams inGS inData curBestCost curBestGraphs netPenaltyFactor reoptimizedSplitGraph splitGraphSimple splitGraphCost graphRoot prunedGraphRootIndex originalConnectionOfPruned (drop numGraphsToExamine rejoinEdges) edgesInPrunedGraph newSAParams

         -- )

-- | singleJoin' is a wrapper arounds singleJoin to allow parMap with individual SAParams
singleJoin'   :: SwapParams 
              -> GlobalSettings
              -> ProcessedData
              -> DecoratedGraph
              -> SimpleGraph
              -> VertexCost
              -> LG.Node
              -> LG.Node
              -> VertexCost
              -> [LG.LEdge EdgeInfo]
              -> (Maybe SAParams, LG.LEdge EdgeInfo)
              -> ([PhylogeneticGraph], Maybe SAParams)
singleJoin' swapParams inGS inData splitGraph splitGraphSimple splitCost prunedGraphRootIndex originalConnectionOfPruned curBestCost edgesInPrunedGraph (inSimAnnealParams, targetEdge) =
   singleJoin swapParams inGS inData splitGraph splitGraphSimple splitCost prunedGraphRootIndex originalConnectionOfPruned curBestCost edgesInPrunedGraph inSimAnnealParams targetEdge

-- | singleJoin takes optimized split graph, split cost, target edge, swap type (ie TBR/SPR/NNI)
-- and "rejoins" the split graph to a single graph--creates joined graph and calculates a heuristic graph cost
-- based on the union assignment of the edge and its distance to the root vertex of the pruned graph
-- if TBR checks all edges in pruned graph with readdition edge (shortcircuits if steepest  == True)
-- always deletes connecting edge to pruned part and readds--this because sometimes it is there and sometimes not (depending on
-- if SPR for terminal etc) and can create parallel edges with different weights (0.0 or not) so just remove to be sure.
-- TBR uses dynamic epsilon even in SPR moves--SPR does not
singleJoin :: SwapParams
           -> GlobalSettings
           -> ProcessedData
           -> DecoratedGraph
           -> SimpleGraph
           -> VertexCost
           -> LG.Node
           -> LG.Node
           -> VertexCost
           -> [LG.LEdge EdgeInfo]
           -> Maybe SAParams
           -> LG.LEdge EdgeInfo
           -> ([PhylogeneticGraph], Maybe SAParams)
singleJoin swapParams inGS inData splitGraph splitGraphSimple splitCost prunedGraphRootIndex originalConnectionOfPruned curBestCost edgesInPrunedGraph inSimAnnealParams targetEdge@(u,v, _) =
   -- trace ("Rejoinging: " ++ (show $ LG.toEdge targetEdge)) (
   let newEdgeList = [(u, originalConnectionOfPruned, 0.0),(originalConnectionOfPruned, v, 0.0),(originalConnectionOfPruned, prunedGraphRootIndex, 0.0)]

       charInfoVV = fmap thd3 $ thd3 inData

       -- graphTYpoe with IA field
       -- only uswe wqhere they exist
       (makeEdgeDataFunction, edgeJoinFunction) = if graphType inGS == HardWired then (M.makeEdgeData False True, edgeJoinDelta False) 
                                                  else if not (useIA inGS) then (M.makeEdgeData False True, edgeJoinDelta False)
                                                  else (M.makeEdgeData True True, edgeJoinDelta True)

       
       -- set edge union creation type to IA-based, filtering gaps (should be linear)
       -- hence True True
       targetEdgeData = makeEdgeDataFunction splitGraph charInfoVV targetEdge
       -- this for using DO for edge O(n^2)
       --targetEdgeData = M.makeEdgeData doIA (not doIA) splitGraph charInfoVV targetEdge

       --this for SPR/NNI only
       prunedRootVertexData = vertData $ fromJust $ LG.lab splitGraph prunedGraphRootIndex

       -- rejoin should always be DO based on edge and pruned root but can be different lengths (unless Static Approx)
       sprReJoinCost = edgeJoinFunction charInfoVV prunedRootVertexData targetEdgeData

       sprNewGraph = LG.insEdges newEdgeList $ LG.delEdges [(u,v),(originalConnectionOfPruned, prunedGraphRootIndex)] splitGraphSimple

       -- here when needed--correct graph is issue in network 
       -- swap can screw up time consistency and other issues
       sprNewGraphChecked = if graphType inGS == Tree then sprNewGraph
                            else if not (LG.isPhylogeneticGraph sprNewGraph) then LG.empty
                            else sprNewGraph
                            
       rediagnosedSPRGraph = T.multiTraverseFullyLabelGraph inGS inData False False Nothing sprNewGraphChecked


       -- Filter for bridge edges for TBR when needed
       edgesInPrunedGraph' = if (graphType inGS == Tree) || LG.isTree splitGraphSimple then edgesInPrunedGraph
                             else fmap fst $ filter ((== True) . snd) $ zip edgesInPrunedGraph (fmap (LG.isBridge splitGraphSimple) (fmap LG.toEdge edgesInPrunedGraph))
   in

   -- do redo orginal graph join
   if originalConnectionOfPruned `elem` [u,v] then ([], inSimAnnealParams)

   -- regular swap
   else if isNothing inSimAnnealParams then

      -- SPR or no TBR rearrangements
      if ((swapType swapParams) == SPR) || ((length edgesInPrunedGraph) < 4) then
         if (sprReJoinCost + splitCost) <= curBestCost then
            if LG.isEmpty sprNewGraphChecked then ([], inSimAnnealParams)
            -- else if (GO.parentsInChainGraph . thd6) rediagnosedSPRGraph then ([], inSimAnnealParams)
            else if snd6 rediagnosedSPRGraph <= curBestCost then ([rediagnosedSPRGraph], inSimAnnealParams)
            else ([], inSimAnnealParams)
         else ([], inSimAnnealParams)

      -- Full TBR
      else if ((swapType swapParams) == TBR) then

         -- do TBR stuff returning SPR results if heuristic better
         let sprResult = if (sprReJoinCost + splitCost) <= curBestCost + (sprReJoinCost * (dynamicEpsilon inGS)) then
                           if LG.isEmpty sprNewGraphChecked then []
                           -- else if (GO.parentsInChainGraph . thd6) rediagnosedSPRGraph then []
                           else if snd6 rediagnosedSPRGraph <= curBestCost then [rediagnosedSPRGraph]
                           else []
                         else []
             (tbrResult, _) = tbrJoin swapParams inGS inData splitGraph splitGraphSimple splitCost prunedGraphRootIndex originalConnectionOfPruned curBestCost edgesInPrunedGraph' inSimAnnealParams targetEdge
         in
         if (not . null) sprResult then (sprResult, inSimAnnealParams)

         else (tbrResult, inSimAnnealParams)

      -- TBRAlternate can skip SPR moves since done already in alternate scenario
      else
         let (tbrResult, _) = tbrJoin swapParams inGS inData splitGraph splitGraphSimple splitCost prunedGraphRootIndex originalConnectionOfPruned curBestCost edgesInPrunedGraph' inSimAnnealParams targetEdge
         in
         (tbrResult, inSimAnnealParams)


   -- simulated annealing/Drift swap
   else
      -- check if spr better always return if so
      let sprResult = if (sprReJoinCost + splitCost) <= curBestCost then
                        if LG.isEmpty sprNewGraphChecked then []
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
         else if ((swapType swapParams) == SPR) || ((length edgesInPrunedGraph) < 4) then
            ([], newSAParams)

         -- tbr
         else
            tbrJoin swapParams inGS inData splitGraph splitGraphSimple splitCost prunedGraphRootIndex originalConnectionOfPruned curBestCost edgesInPrunedGraph' newSAParams targetEdge
   -- )

-- | edgeJoinDelta calculates heuristic cost for joining pair edges
-- I a filed is faster--but has to be there and not for harwired
edgeJoinDelta :: Bool -> V.Vector (V.Vector CharInfo) -> VertexBlockData -> VertexBlockData -> VertexCost
edgeJoinDelta useIA charInfoVV edgeA edgeB =
   if (not useIA) then V.sum $ fmap V.sum $ fmap (fmap snd) $ POSW.createVertexDataOverBlocks edgeA edgeB charInfoVV []
   else V.sum $ fmap V.sum $ fmap (fmap snd) $ POSW.createVertexDataOverBlocksStaticIA edgeA edgeB charInfoVV []


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
tbrJoin :: SwapParams
        -> GlobalSettings
        -> ProcessedData
        -> DecoratedGraph
        -> SimpleGraph
        -> VertexCost
        -> LG.Node
        -> LG.Node
        -> VertexCost
        -> [LG.LEdge EdgeInfo]
        -> Maybe SAParams
        -> LG.LEdge EdgeInfo
        -> ([PhylogeneticGraph], Maybe SAParams)
tbrJoin swapParams inGS inData splitGraph splitGraphSimple splitCost prunedGraphRootIndex originalConnectionOfPruned curBestCost edgesInPrunedGraph inSimAnnealParams targetEdge =
   
   -- this is for networks stopping TBR rearrangements with there are nwetowrk edegs incolved in pruned part of graph
   let hasNetEdges = if graphType inGS == Tree then False 
                     else null $ filter ((== True) . LG.isNetworkLabEdge splitGraph) edgesInPrunedGraph

   in

   if null edgesInPrunedGraph then ([], inSimAnnealParams)

   else if hasNetEdges then ([], inSimAnnealParams)

   else
      -- get target edge data\
      -- always using IA for union, but filtering out gaps (doIA (not doIA))
      let charInfoVV = fmap thd3 $ thd3 inData
          
          -- graphTYpoe with IA field
          -- only uswe wqhere they exist
          (makeEdgeDataFunction, edgeJoinFunction) = if graphType inGS == HardWired then (M.makeEdgeData False True, edgeJoinDelta False) 
                                                     else if not (useIA inGS) then (M.makeEdgeData False True, edgeJoinDelta False)
                                                     else (M.makeEdgeData True True, edgeJoinDelta True)

          targetEdgeData = makeEdgeDataFunction splitGraph charInfoVV targetEdge

      in

      -- logic for annealing/Drift  regular swap first
      if isNothing inSimAnnealParams then
         if not (steepest swapParams) then
            -- get heuristic delta joins for edges in pruned graph
            let rerootEdgeList = filter ((/= prunedGraphRootIndex) . fst3) $ filter ((/= originalConnectionOfPruned) . fst3) edgesInPrunedGraph
                
                -- True True to use IA fields and filter gaps
                rerootEdgeDataList = PU.seqParMap (parStrategy $ lazyParStrat inGS) (makeEdgeDataFunction splitGraph charInfoVV) rerootEdgeList

                rerootEdgeDeltaList = fmap (+ splitCost) $ PU.seqParMap (parStrategy $ lazyParStrat inGS) (edgeJoinFunction charInfoVV targetEdgeData) rerootEdgeDataList

                -- check for possible better/equal graphs and verify
                deltaAdjustmentJoinCost = (curBestCost - splitCost) * (dynamicEpsilon inGS)
                candidateEdgeList = fmap fst $ filter ((<= (curBestCost + deltaAdjustmentJoinCost)) . snd) (zip rerootEdgeList rerootEdgeDeltaList)
                candidateJoinedGraphList' = PU.seqParMap (parStrategy $ lazyParStrat inGS) (rerootPrunedAndMakeGraph  splitGraphSimple prunedGraphRootIndex originalConnectionOfPruned targetEdge) candidateEdgeList

                -- check for grap wierdness if network
                candidateJoinedGraphList = if graphType inGS == Tree then candidateJoinedGraphList'
                                           else filter LG.isPhylogeneticGraph candidateJoinedGraphList'

                -- check for graph wierdness
                -- rediagnosedGraphList = filter (not . GO.parentsInChainGraph . thd6) $ filter ((<= curBestCost) . snd6) $ PU.seqParMap (parStrategy $ lazyParStrat inGS) (T.multiTraverseFullyLabelGraph inGS inData False False Nothing) candidateJoinedGraphList

                rediagnosedGraphList = filter ((<= curBestCost) . snd6) $ PU.seqParMap (parStrategy $ lazyParStrat inGS) (T.multiTraverseFullyLabelGraph inGS inData False False Nothing) candidateJoinedGraphList


                -- for debugging
                -- allRediagnosedList = PU.seqParMap PU.myStrategy (T.multiTraverseFullyLabelGraph inGS inData False False Nothing) (PU.seqParMap PU.myStrategy (rerootPrunedAndMakeGraph  splitGraphSimple  prunedGraphRootIndex originalConnectionOfPruned targetEdge) rerootEdgeList)

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
                
                -- True True to use IA fields and filter gaps
                rerootEdgeDataList = PU.seqParMap (parStrategy $ lazyParStrat inGS) (makeEdgeDataFunction splitGraph charInfoVV) rerootEdgeList
                rerootEdgeDeltaList = fmap (+ splitCost) $ PU.seqParMap (parStrategy $ lazyParStrat inGS) (edgeJoinFunction charInfoVV targetEdgeData) rerootEdgeDataList

                -- check for possible better/equal graphs and verify
                deltaAdjustmentJoinCost = (curBestCost - splitCost) * (dynamicEpsilon inGS)
                candidateEdgeList = fmap fst $ filter ((<= (curBestCost + deltaAdjustmentJoinCost)) . snd) (zip rerootEdgeList rerootEdgeDeltaList)
                candidateJoinedGraphList = PU.seqParMap (parStrategy $ lazyParStrat inGS) (rerootPrunedAndMakeGraph splitGraphSimple prunedGraphRootIndex originalConnectionOfPruned targetEdge) candidateEdgeList
                rediagnosedGraphList = filter ((<= curBestCost) . snd6) $ PU.seqParMap (parStrategy $ lazyParStrat inGS) (T.multiTraverseFullyLabelGraph inGS inData False False Nothing) candidateJoinedGraphList-- get

            in
            -- trace ("TBR steepest: " ++ (show $ length rerootEdgeList) ++ " edges to go " ++ (show $ length $ (drop numEdgesToExamine edgesInPrunedGraph))) (
            if null candidateEdgeList then
               tbrJoin swapParams inGS inData splitGraph splitGraphSimple splitCost prunedGraphRootIndex originalConnectionOfPruned curBestCost (drop numEdgesToExamine edgesInPrunedGraph) inSimAnnealParams targetEdge
            else if null rediagnosedGraphList then
               tbrJoin swapParams inGS inData splitGraph splitGraphSimple splitCost prunedGraphRootIndex originalConnectionOfPruned curBestCost (drop numEdgesToExamine edgesInPrunedGraph) inSimAnnealParams targetEdge
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
                
                -- True True to use IA fields and filter gaps
                rerootEdgeDataList = PU.seqParMap (parStrategy $ lazyParStrat inGS) (makeEdgeDataFunction splitGraph charInfoVV) rerootEdgeList
                rerootEdgeDeltaList = fmap (+ splitCost) $ PU.seqParMap (parStrategy $ lazyParStrat inGS) (edgeJoinFunction charInfoVV targetEdgeData) rerootEdgeDataList

                minDelta = if (not . null) rerootEdgeDeltaList then minimum $ rerootEdgeDeltaList
                           else infinity
                minEdgeList = if (not . null) rerootEdgeDeltaList then fmap fst $ filter ((== minDelta) . snd)  (zip rerootEdgeList rerootEdgeDeltaList)
                              else []

                -- check for possible better/equal graphs and verify
                candidateJoinedGraphList = PU.seqParMap (parStrategy $ lazyParStrat inGS) (rerootPrunedAndMakeGraph splitGraphSimple prunedGraphRootIndex originalConnectionOfPruned targetEdge) minEdgeList
                rediagnosedGraphList = PU.seqParMap (parStrategy $ lazyParStrat inGS) (T.multiTraverseFullyLabelGraph inGS inData False False Nothing) candidateJoinedGraphList

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
                  tbrJoin swapParams inGS inData splitGraph splitGraphSimple splitCost prunedGraphRootIndex originalConnectionOfPruned curBestCost (drop numEdgesToExamine edgesInPrunedGraph) newSAParams targetEdge
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
   -- trace ("RSGFV: " ++ (show startVertex)) (
   if doIA then
      -- only reoptimize the IA states for dynamic characters
      reoptimizeSplitGraphFromVertexIA inGS inData netPenaltyFactor inSplitGraph startVertex prunedSubGraphRootVertex
   else
      -- perform full optimizations of nodes
      -- these required for full optimization
      let nonExactCharacters = U.getNumberSequenceCharacters (thd3 inData)
          origGraph = inSplitGraph -- thd6 origPhyloGraph
          leafGraph = if graphType inGS == SoftWired then LG.extractLeafGraph origGraph -- POSW.makeLeafGraphSoftWired inGS inData -- LG.extractLeafGraph origGraph
                      else LG.extractLeafGraph origGraph
          calcBranchLengths = False

          -- this for multitravers in swap for softwired to turn off
          multiTraverse = if graphType inGS == Tree then multiTraverseCharacters inGS
                         else False

          -- create simple graph version of split for post order pass
          splitGraphSimple = GO.convertDecoratedToSimpleGraph inSplitGraph


          -- create optimized base graph
          -- False for staticIA
          (postOrderBaseGraph, _) = T.generalizedGraphPostOrderTraversal (inGS {graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse}) nonExactCharacters inData leafGraph False (Just startVertex) splitGraphSimple


          fullBaseGraph = PRE.preOrderTreeTraversal (inGS {graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse}) (finalAssignment inGS) False calcBranchLengths (nonExactCharacters > 0) startVertex True postOrderBaseGraph

          -- create fully optimized pruned graph.  Post order tehn preorder

          -- get root node of pruned graph--parent since that is the full pruned piece (keeping that node for addition to base graph and edge creation)
          startPrunedNode = (prunedSubGraphRootVertex, fromJust $ LG.lab origGraph prunedSubGraphRootVertex)
          startPrunedParentNode = head $ LG.labParents origGraph prunedSubGraphRootVertex
          startPrunedParentEdge = (fst startPrunedParentNode, prunedSubGraphRootVertex, dummyEdge)


          -- False for staticIA
          (postOrderPrunedGraph, _) = T.generalizedGraphPostOrderTraversal (inGS {graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse}) nonExactCharacters inData leafGraph False (Just prunedSubGraphRootVertex) splitGraphSimple


          -- False for staticIA
          fullPrunedGraph = PRE.preOrderTreeTraversal (inGS {graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse}) (finalAssignment inGS) False calcBranchLengths (nonExactCharacters > 0) prunedSubGraphRootVertex True postOrderPrunedGraph

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

          {-
          -- check fo unlabbeld nodes
          coninicalNodes =  LG.labNodes fullSplitGraph
          nodeLabels = fmap (LG.lab fullSplitGraph) (fmap fst coninicalNodes)
          unlabelledNodes = filter ((== Nothing) .snd) $ (zip (fmap fst coninicalNodes) nodeLabels)
          -}

      in
      -- trace ("RSGV: " ++ (show $ fmap fst unlabelledNodes))
      {-
      trace ("Orig graph cost " ++ (show $ subGraphCost $ fromJust $ LG.lab origGraph startVertex) ++ " Base graph cost " ++ (show $ snd6 fullBaseGraph) ++ " pruned subgraph cost " ++ (show prunedCost) ++ " at node " ++ (show prunedSubGraphRootVertex) ++ " parent " ++ (show $ fst startPrunedParentNode)

         ++ "\nBaseGraphNodes\n" ++ (show $ L.sort  $ fmap fst baseGraphNonRootNodes) ++ "\nPruned nodes from root: " ++ "\n" ++ (show $ fmap fst $ startPrunedNode : prunedGraphNonRootNodes)
         ++ "\nSplit Graph\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph fullSplitGraph)
         ++ "\nOrig graph:\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph origGraph))
      -}
      -- trace ("reoptimizeSplitGraphFromVertex: " ++ (show splitGraphCost))
      (fullSplitGraph, splitGraphCost)

      -- )

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
            leafGraph = if graphType inGS == SoftWired then
                              GO.copyIAFinalToPrelim  $ LG.extractLeafGraph origGraph --POSW.makeLeafGraphSoftWired inGS inData -- LG.extractLeafGraph origGraph
                         else GO.copyIAFinalToPrelim  $ LG.extractLeafGraph origGraph
            calcBranchLengths = False
            
            -- this for multitravers in swap for softwired to turn off
            multiTraverse = if graphType inGS == Tree then multiTraverseCharacters inGS
                         else False

            -- create simple graph version of split for post order pass
            splitGraphSimple = GO.convertDecoratedToSimpleGraph inSplitGraph

            --Create base graph
            -- create postorder assignment--but only from single traversal
            -- True flag fior staticIA
            postOrderBaseGraph = POSW.postOrderTreeTraversal (inGS {graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse}) inData leafGraph True (Just startVertex) splitGraphSimple
            baseGraphCost = snd6 postOrderBaseGraph

            -- True flag fior staticIA
            fullBaseGraph = PRE.preOrderTreeTraversal (inGS {graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse}) (finalAssignment inGS) True calcBranchLengths (nonExactCharacters > 0) startVertex True postOrderBaseGraph

            {-
            localRootCost = if (rootCost inGS) == NoRootCost then 0.0
                              else if (rootCost inGS) == Wheeler2015Root then POSW.getW15RootCost inGS postOrderBaseGraph
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
            postOrderPrunedGraph =  POSW.postOrderTreeTraversal (inGS {graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse}) inData leafGraph True (Just prunedSubGraphRootVertex) splitGraphSimple
            prunedGraphCost = snd6 postOrderPrunedGraph

            -- True flag fior staticIA
            fullPrunedGraph = PRE.preOrderTreeTraversal (inGS {graphFactor = NoNetworkPenalty, multiTraverseCharacters = multiTraverse}) (finalAssignment inGS) True calcBranchLengths (nonExactCharacters > 0) prunedSubGraphRootVertex True postOrderPrunedGraph

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


