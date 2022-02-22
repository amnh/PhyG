{- |
Module      :  GeneticAlgorithm.hs
Description :  Module specifying graph sGeneticAlgorithm functions
Copyright   :  (c) 2022 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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

module Search.GeneticAlgorithm ( geneticAlgorithm
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
import qualified Search.Swap as S
import qualified Search.NetworkAddDelete as N
import qualified Search.Fuse as F


-- | geneticAlgorithm takes arguments and performs genetic algorithm on input graphs
-- the process follows several steps
-- 1) input graphs are mutated 
--       this step is uncharacteristically first so that is can operate on 
--       graphs that have been "fused" (recombined) already
--       mutated graphs are added up to popsize
--       if input graphs are already at the population size, an equal number of mutants are added (exceeding input popsize)
-- 2) graph are recombined using fusing operations
-- 3) population undergoes selection to population size (unique graphs)
--       selection based on delta with best graph and severity factor on (0,Inf) 1 pure cost delta < 1 more severe, > 1 less severe
--       if "elitist" (default) 'best' graphs are always selected to ensure no worse.
-- 4) operation repearts for number of generations
geneticAlgorithm :: GlobalSettings -> ProcessedData -> Int -> Bool -> Int -> Int -> Int -> Int -> Double -> Int -> [PhylogeneticGraph] -> ([PhylogeneticGraph], Int)
geneticAlgorithm inGS inData rSeed doElitist keepNum popSize generations generationCounter severity recombinations inGraphList =
    if null inGraphList then ([], 0)
    else if generationCounter == generations then  (inGraphList, generationCounter)
    else 
        let seedList = randomIntList rSeed

            -- get elite list of best solutions
            initialEliteList = GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] inGraphList

            -- mutate input graphs, produces number input, limited to popsize 
            mutatedGraphList' = zipWith (mutateGraph inGS inData) (randomIntList $ head seedList) $ takeRandom (seedList !! 1) popSize inGraphList

            -- adjust to correct populationsize if input number < popSize
            mutatedGraphList = if length mutatedGraphList' >= popSize then mutatedGraphList'
                               else 
                                    let numShort = popSize - (length mutatedGraphList')
                                        additionalMutated = zipWith (mutateGraph inGS inData) (randomIntList $ seedList !! 2) $ takeRandom (seedList !! 3) numShort inGraphList
                                    in
                                    mutatedGraphList' ++ additionalMutated

            -- get unique graphs, no point in recombining repetitions 
            uniqueMutatedGraphList = GO.selectPhylogeneticGraph [("unique","")] 0 ["unique"] (mutatedGraphList ++ inGraphList)

            -- recombine elite with mutated and mutated with mutated
            recombineSwap = getRandomElement (seedList !! 4) ["nni", "spr", "tbr"]
            (doNNI, doSPR, doTBR) = if recombineSwap == "nni" then (True, False, False)
                                   else if recombineSwap == "spr" then (False, True, False)
                                   else (False, False, True)
            doSteepest = True
            doAll = False
            returnBest = False
            returnUnique = True
            singleRound = False
            fusePairs = Just recombinations
            randomPairs = True

            (recombinedGraphList, _) = F.fuseAllGraphs inGS inData (drop 6 seedList) (2 * popSize) (maxBound :: Int) 0  doNNI doSPR doTBR doSteepest doAll returnBest returnUnique singleRound fusePairs randomPairs uniqueMutatedGraphList

            -- selection of graphs population
            -- unique sorted on cost so getting unique with lowest cost
            selectedGraphs = take popSize $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] recombinedGraphList
            newCost = snd6 $ head selectedGraphs
            
        in
        trace ("\tGA " ++ (show $ snd6 $ head initialEliteList) ++ " -> " ++ (show newCost) ++ "\nInGraphs " ++ (show $ L.sort $ fmap snd6 inGraphList)
            ++ "\nMutated " ++ (show $  L.sort $ fmap snd6 mutatedGraphList)
            ++ "\nRecombined " ++ recombineSwap ++ " " ++ (show $  L.sort $ fmap snd6 recombinedGraphList)) (
        
        -- if new graphs better cost then take those
        if newCost < (snd6 $ head initialEliteList) then 
            geneticAlgorithm inGS inData (seedList !! 5) doElitist keepNum popSize generations (generationCounter + 1) severity recombinations selectedGraphs

        -- if noew graphs not better then add in elites toi ensure monotonic decrease in cost
        else 
            let newGraphList = take keepNum $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] (initialEliteList ++ selectedGraphs)
            in
            geneticAlgorithm inGS inData (seedList !! 5) doElitist keepNum popSize generations (generationCounter + 1) severity recombinations newGraphList
        )

-- | mutateGraph mutates a graph using drift functionality
mutateGraph :: GlobalSettings -> ProcessedData -> Int -> PhylogeneticGraph -> PhylogeneticGraph
mutateGraph inGS inData rSeed inGraph =
    if LG.isEmpty (fst6 inGraph) then error "Empty graph in mutateGraph"
    else 
        let randList = randomIntList rSeed
            saValues = Just $ SAParams { method = Drift
                                , numberSteps = 0
                                , currentStep = 0
                                , randomIntegerList = randomIntList rSeed
                                , rounds      = 1
                                , driftAcceptEqual  = 0.67
                                , driftAcceptWorse  = 0.0
                                -- this could be an important factor don't want too severe, but significant
                                , driftMaxChanges   = getRandomElement (randList !! 1) [1,2,4] -- or something
                                , driftChanges      = 0 
                                } 
        in
        let --randomize edit type
            editType = getRandomElement (randList !! 0) ["swap", "netEdge"]

            -- randomize Swap parameters 
            numToKeep = 1
            maxMoveEdgeDist = 10
            steepest = True
            hardwiredSPR = False 
            doIA = False
            returnMutated = True
            inSimAnnealParams = saValues
            swapType = getRandomElement (randList !! 2) ["spr","tbr"]

            --radnomize network edit parameters
            netEditType = getRandomElement (randList !! 3) ["netAdd", "netDelete", "netMove"]
            doRandomOrder = True

        in
        
        -- only swap stuff for tree
        if graphType inGS == Tree || (LG.isTree (fst6 inGraph) && netEditType /= "netadd") then 
            head $ fst $ S.swapSPRTBR swapType inGS inData numToKeep maxMoveEdgeDist steepest hardwiredSPR doIA returnMutated (inSimAnnealParams, inGraph)

        -- graphs choose what type of mutation at random
        else 
            if editType == "swap" then
                head $ fst $ S.swapSPRTBR swapType inGS inData numToKeep maxMoveEdgeDist steepest hardwiredSPR doIA returnMutated (inSimAnnealParams, inGraph)
            else 
                -- move only for Hardwired
                if (graphType inGS) == HardWired then 
                    head $ fst $ N.moveAllNetEdges inGS inData (randList !! 4) numToKeep 0 returnMutated steepest doRandomOrder ([], infinity) (inSimAnnealParams, [inGraph]) 
                    

                -- SoftWired
                else
                    if netEditType == "netMove" then
                        head $ fst $ N.moveAllNetEdges inGS inData (randList !! 4) numToKeep 0 returnMutated steepest doRandomOrder ([], infinity) (inSimAnnealParams, [inGraph])
                    
                    else if netEditType == "netadd" then
                        head $ fst $ N.insertAllNetEdges inGS inData (randList !! 4) numToKeep 0 returnMutated steepest doRandomOrder ([], infinity) (inSimAnnealParams, [inGraph])
                        
                    -- net delete
                    else 
                        head $ fst $ N.deleteAllNetEdges inGS inData (randList !! 4) numToKeep 0 returnMutated steepest doRandomOrder ([], infinity) (inSimAnnealParams, [inGraph])
