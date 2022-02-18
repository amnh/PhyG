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
geneticAlgorithm :: GlobalSettings -> ProcessedData -> Int -> Bool -> Int -> Int -> Int -> Double -> Int -> [PhylogeneticGraph] -> ([PhylogeneticGraph], Int)
geneticAlgorithm inGS inData rSeed doElitist keepNum popSize generations severity recombinations inGraphList =
    if null inGraphList then ([], 0)
    else 
        let seedList = randomIntList rSeed

            -- get elite list of best solutions
            initialEliteList = GO.selectPhylogeneticGraph [("best", (show keepNum))] 0 ["best"] inGraphList

            -- mutate input graphs
            mutatedGraphList = zipWith (mutateGraph inGS) (randomIntList $ head seedList) $ takeRandom (seedList !! 1) popSize inGraphList

            -- recom,bine elite with mutated and mutated with mutated
            recombinedGraphList = mutatedGraphList

            -- selection of graphs population
            selectedGraphs = recombinedGraphList
            newCost = minimum $ fmap snd6 selectedGraphs

            -- final list with elits added back
            newGraphs = take popSize $ GO.selectPhylogeneticGraph [("unique", (show keepNum))] 0 ["unique"] inGraphList
        in
        trace ("GA " ++ (show $ snd6 $ head initialEliteList) ++ " -> " ++ (show newCost)) (
        if newCost < (snd6 $ head initialEliteList) then 
            (newGraphs, 1)
        else 
            (take keepNum $ GO.selectPhylogeneticGraph [("unique", (show popSize))] 0 ["unique"] (initialEliteList ++ newGraphs), 1)
        )

-- | mutateGraph mutates a graph using drift functionality
mutateGraph :: GlobalSettings -> Int -> PhylogeneticGraph -> PhylogeneticGraph
mutateGraph inGS rSeed inGraph =
    if LG.isEmpty (fst6 inGraph) then error "Empty graph in mutateGraph"
    else 
        let saValues = SAParams { method = Drift
                                , numberSteps = 0
                                , currentStep = 0
                                , randomIntegerList = randomIntList rSeed
                                , rounds      = 1
                                , driftAcceptEqual  = 0.67
                                , driftAcceptWorse  = 0.0
                                -- this could be an important factor don't want too severe, but significant
                                , driftMaxChanges   = 2
                                , driftChanges      = 0 
                                } 
        in
        -- only swap stuff for tree
        if graphType inGS == Tree then inGraph


        -- graphs choose what type of mutation at random
        else inGraph
