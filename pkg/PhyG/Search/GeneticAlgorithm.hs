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

import           Debug.Trace
import           GeneralUtilities
import qualified Graphs.GraphOperations  as GO
import qualified Search.Fuse             as F
import qualified Search.NetworkAddDelete as N
import qualified Search.Swap             as S
import           Types.Types
import qualified Utilities.LocalGraph    as LG


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
geneticAlgorithm :: GlobalSettings
                 -> ProcessedData
                 -> Int
                 -> Bool
                 -> Int
                 -> Int
                 -> Int
                 -> Int
                 -> Int
                 -> Double
                 -> Int
                 -> Int
                 -> Int
                 -> [ReducedPhylogeneticGraph]
                 -> ([ReducedPhylogeneticGraph], Int)
geneticAlgorithm inGS inData rSeed doElitist maxNetEdges keepNum popSize generations generationCounter severity recombinations stopCount stopNum inGraphList =
    if null inGraphList then ([], 0)
    else if generationCounter == generations then  (inGraphList, generationCounter)
    else if stopCount >= stopNum then (inGraphList, generationCounter)
    else
        trace ("Genetic algorithm generation: " ++ (show generationCounter)) (
        let seedList = randomIntList rSeed

            -- get elite list of best solutions
            initialEliteList = GO.selectGraphs Best (maxBound::Int) 0.0 (-1) inGraphList

            -- mutate input graphs, produces number input, limited to popsize
            mutatedGraphList' = zipWith (mutateGraph inGS inData maxNetEdges) (randomIntList $ head seedList) $ takeRandom (seedList !! 1) popSize inGraphList

            -- adjust to correct populationsize if input number < popSize
            mutatedGraphList = if length mutatedGraphList' >= popSize then mutatedGraphList'
                               else
                                    let numShort = popSize - (length mutatedGraphList')
                                        additionalMutated = zipWith (mutateGraph inGS inData maxNetEdges) (randomIntList $ seedList !! 2) $ takeRandom (seedList !! 3) numShort inGraphList
                                    in
                                    mutatedGraphList' ++ additionalMutated

            -- get unique graphs, no point in recombining repetitions
            uniqueMutatedGraphList = GO.selectGraphs Unique (maxBound::Int) 0.0 (-1) (mutatedGraphList ++ inGraphList)

            -- recombine elite with mutated and mutated with mutated
            recombineSwap = getRandomElement (seedList !! 4) [None, NNI, SPR] --  these take too long, "tbr", "alternate"]

            -- options to join via union choices or all in fuse
            -- this is ignored ofr now in fuse--JoinAll is what it does
            joinType =  getRandomElement (seedList !! 6) [JoinAlternate, JoinAll ]

            doSteepest = True
            returnBest = False
            returnUnique = True
            singleRound = False
            fusePairs = Just recombinations
            randomPairs = True
            reciprocal = False

            -- populate SwapParams structure
            swapParams = SwapParams { swapType = recombineSwap
                                    , joinType = joinType 
                                    , atRandom = True -- randomize swap order
                                    , keepNum  = (2 * popSize)
                                    , maxMoveEdgeDist = (maxBound :: Int)
                                    , steepest = doSteepest
                                    , joinAlternate = False -- not working now 
                                    , doIA = False
                                    , returnMutated = False 
                                    }


            (recombinedGraphList, _) = F.fuseAllGraphs swapParams inGS inData (drop 6 seedList) 0 returnBest returnUnique singleRound fusePairs randomPairs reciprocal uniqueMutatedGraphList

            -- selection of graphs population
            -- unique sorted on cost so getting unique with lowest cost
            selectedGraphs = GO.selectGraphs Unique popSize 0.0 (-1) recombinedGraphList
            newCost = snd5 $ head selectedGraphs

        in
        {-
        trace ("\tGA " ++ (show $ snd5 $ head initialEliteList) ++ " -> " ++ (show newCost) ++ "\nInGraphs " ++ (show $ L.sort $ fmap snd5 inGraphList)
            ++ "\nMutated " ++ (show $  L.sort $ fmap snd5 mutatedGraphList)
            ++ "\nRecombined " ++ recombineSwap ++ " " ++ (show $  L.sort $ fmap snd5 recombinedGraphList)) (
        -}
        -- if new graphs better cost then take those
        if newCost < (snd5 $ head initialEliteList) then
            geneticAlgorithm inGS inData (seedList !! 5) doElitist maxNetEdges keepNum popSize generations (generationCounter + 1) severity recombinations 0 stopNum selectedGraphs

        -- if new graphs not better then add in elites to ensure monotonic decrease in cost
        else
            let newGraphList = GO.selectGraphs Unique keepNum 0.0 (-1) (initialEliteList ++ selectedGraphs)
            in
            geneticAlgorithm inGS inData (seedList !! 5) doElitist maxNetEdges keepNum popSize generations (generationCounter + 1) severity recombinations (stopCount + 1) stopNum newGraphList
        )

-- | mutateGraph mutates a graph using drift functionality
mutateGraph :: GlobalSettings -> ProcessedData -> Int -> Int -> ReducedPhylogeneticGraph -> ReducedPhylogeneticGraph
mutateGraph inGS inData maxNetEdges rSeed inGraph =
    if LG.isEmpty (fst5 inGraph) then error "Empty graph in mutateGraph"
    else
        let joinType = JoinAll -- keep selection of rejoins based on all possibilities
            atRandom = True -- randomize split and rejoin edge orders
            randList = randomIntList rSeed
            saValues = Just $ SAParams  { method = getRandomElement (randList !! 0) [Drift, SimAnneal]
                                        , numberSteps = getRandomElement (randList !! 1) [5, 10, 20]
                                        , currentStep = 0
                                        , randomIntegerList = randomIntList rSeed
                                        , rounds      = 1
                                        , driftAcceptEqual  = 0.5
                                        , driftAcceptWorse  = 2.0
                                        -- this could be an important factor don't want too severe, but significant
                                        , driftMaxChanges   = getRandomElement (randList !! 2) [5, 10, 20] -- or something
                                        , driftChanges      = 0
                                        }
        in
        let --randomize edit type
            editType = getRandomElement (randList !! 3) ["swap", "netEdge"]

            -- randomize Swap parameters
            alternate = False
            numToKeep = 5
            maxMoveEdgeDist = 10000
            steepest = True
            doIA = False
            returnMutated = True
            inSimAnnealParams = saValues
            swapType = getRandomElement (randList !! 4) [SPR,Alternate]

            --randomize network edit parameters
            netEditType = getRandomElement (randList !! 5) ["netAdd", "netDelete", "netAddDelete"] -- , "netMove"]
            doRandomOrder = True
            maxRounds = getRandomElement (randList !! 6) [1..5]

            -- populate SwapParams structure
            swapParams = SwapParams { swapType = swapType
                                         , joinType = joinType 
                                         , atRandom = atRandom
                                         , keepNum  = numToKeep
                                         , maxMoveEdgeDist = maxMoveEdgeDist
                                         , steepest = steepest
                                         , joinAlternate = alternate 
                                         , doIA = doIA
                                         , returnMutated = returnMutated 
                                         }

        in

        -- only swap mutation stuff for tree
        if graphType inGS == Tree || netEditType `notElem`  ["netAdd", "netDelete", "netAddDelete", "netMove"] then
            -- trace ("1")
            let (newGraphList, _) =  S.swapSPRTBR swapParams inGS inData 0 [inGraph][(randList, inSimAnnealParams, inGraph)]
            in
            if (not . null) newGraphList then head newGraphList
            else inGraph
            -- head $ fst $ S.swapSPRTBR swapType inGS inData numToKeep maxMoveEdgeDist steepest alternate doIA returnMutated (inSimAnnealParams, inGraph)

        -- graphs choose what type of mutation at random
        else
            if editType == "swap" then
                -- trace ("2")
                let (newGraphList, _) =  S.swapSPRTBR swapParams inGS inData 0 [inGraph] [(randList,inSimAnnealParams, inGraph)]
                in
                if (not . null) newGraphList then head newGraphList
                else inGraph
                -- head $ fst $ S.swapSPRTBR swapType inGS inData numToKeep maxMoveEdgeDist steepest alternate doIA returnMutated (inSimAnnealParams, inGraph)

            else
                -- move only for Hardwired
                if (graphType inGS) == HardWired then
                    -- trace ("3")
                    let (newGraphList, _) =  N.moveAllNetEdges inGS inData (randList !! 7) maxNetEdges numToKeep 0 returnMutated steepest doRandomOrder ([], infinity) (inSimAnnealParams, [inGraph])
                    in
                    if (not . null) newGraphList then head newGraphList
                    else inGraph
                    -- head $ fst $ N.moveAllNetEdges inGS inData (randList !! 7) maxNetEdges numToKeep 0 returnMutated steepest doRandomOrder ([], infinity) (inSimAnnealParams, [inGraph])


                -- SoftWired
                else
                    if netEditType == "netMove" then
                        -- trace ("4")
                        let (newGraphList, _) =  N.moveAllNetEdges inGS inData (randList !! 7) maxNetEdges numToKeep 0 returnMutated steepest doRandomOrder ([], infinity) (inSimAnnealParams, [inGraph])
                        in
                        if (not . null) newGraphList then head newGraphList
                        else inGraph
                        -- head $ fst $ N.moveAllNetEdges inGS inData (randList !! 7) maxNetEdges numToKeep 0 returnMutated steepest doRandomOrder ([], infinity) (inSimAnnealParams, [inGraph])

                    else if netEditType == "netAdd" then
                        -- trace ("5")
                        let (newGraphList, _) =  N.insertAllNetEdges inGS inData (randList !! 7) maxNetEdges numToKeep maxRounds 0 returnMutated steepest doRandomOrder ([], infinity) (inSimAnnealParams, [inGraph])
                        in
                        if (not . null) newGraphList then head newGraphList
                        else inGraph
                        -- head $ fst $ N.insertAllNetEdges inGS inData (randList !! 7) maxNetEdges numToKeep maxRounds 0 returnMutated steepest doRandomOrder ([], infinity) (inSimAnnealParams, [inGraph])

                    else if netEditType == "netAddDelete" then
                        -- trace ("6")
                        let (newGraphList, _) =  N.addDeleteNetEdges inGS inData (randList !! 7) maxNetEdges numToKeep maxRounds 0 returnMutated steepest doRandomOrder ([], infinity) (inSimAnnealParams, [inGraph])
                        in
                        if (not . null) newGraphList then head newGraphList
                        else inGraph
                        -- head $ fst $ N.addDeleteNetEdges inGS inData (randList !! 7) maxNetEdges numToKeep maxRounds 0 returnMutated steepest doRandomOrder ([], infinity) (inSimAnnealParams, [inGraph])

                    -- net delete
                    else
                        -- trace ("7")
                        let (newGraphList, _) =  N.deleteAllNetEdges inGS inData (randList !! 7) maxNetEdges numToKeep 0 returnMutated steepest doRandomOrder ([], infinity) (inSimAnnealParams, [inGraph])
                        in
                        if (not . null) newGraphList then head newGraphList
                        else inGraph
                        -- head $ fst $ N.deleteAllNetEdges inGS inData (randList !! 7) maxNetEdges numToKeep 0 returnMutated steepest doRandomOrder ([], infinity) (inSimAnnealParams, [inGraph])
