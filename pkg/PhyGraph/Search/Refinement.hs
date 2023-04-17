{- |
Module      :  Refinement.hs
Description :  Module controlling graph refinement functions
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

module Search.Refinement  ( refineGraph
                          , netEdgeMaster
                          , fuseGraphs
                          , swapMaster
                          , geneticAlgorithmMaster
                          ) where

import           Debug.Trace
import           GeneralUtilities
import qualified Graphs.GraphOperations      as GO
import qualified ParallelUtilities           as PU
import           Types.Types
-- import qualified Search.Swap as S
import qualified Commands.Verify             as VER
import           Data.Char
import           Data.Maybe
import qualified Search.Fuse                 as F
import qualified Search.GeneticAlgorithm     as GA
import qualified Search.NetworkAddDelete     as N
import qualified Search.SwapMaster           as SM
import           Text.Read
import           Utilities.Utilities         as U

-- | swapMaster moved to Search.SwapMaster due to very long (>20') compile times
-- with --enalble-profinling
swapMaster :: [Argument]
            -> GlobalSettings
            -> ProcessedData
            -> Int
            -> [PhylogeneticGraph]
            -> [PhylogeneticGraph]
swapMaster = SM.swapMaster


-- | driver for overall refinement
refineGraph :: [Argument] 
            -> GlobalSettings 
            -> ProcessedData 
            -> Int 
            -> [PhylogeneticGraph] 
            -> [PhylogeneticGraph]
refineGraph inArgs inGS inData rSeed inGraphList =
   if null inGraphList then errorWithoutStackTrace "No graphs input to refine"
   else
      let fstArgList = fmap (fmap toLower . fst) inArgs
          sndArgList = fmap (fmap toLower . snd) inArgs
          lcArgList = zip fstArgList sndArgList
          checkCommandList = checkCommandArgs "refineGraph" fstArgList VER.refineArgList
     in
     -- check for valid command options
     if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'GeneticAlgorithm': " ++ show inArgs)
     else
      let doNetAdd = any ((=="netadd").fst) lcArgList
          doNetDel = any ((=="netdel").fst) lcArgList || any ((=="netdelete").fst) lcArgList
          doNetAddDel = any ((=="netadddel").fst) lcArgList ||  any ((=="netadddelete").fst) lcArgList
          doNetMov = any ((=="netmove").fst) lcArgList
          doGenAlg = any ((=="ga").fst) lcArgList || any ((=="geneticalgorithm").fst) lcArgList
      in

      -- network edge edits
      if doNetAdd || doNetDel || doNetAddDel || doNetMov then
         netEdgeMaster inArgs inGS inData rSeed inGraphList

      -- genetic algorithm
      else if doGenAlg then
         geneticAlgorithmMaster inArgs inGS inData rSeed inGraphList

      else error "No refinement operation specified"

-- | geneticAlgorithmMaster takes arguments and performs genetic algorithm on input graphs
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
geneticAlgorithmMaster :: [Argument] 
                       -> GlobalSettings 
                       -> ProcessedData 
                       -> Int 
                       -> [PhylogeneticGraph] 
                       -> [PhylogeneticGraph]
geneticAlgorithmMaster inArgs inGS inData rSeed inGraphList =
   if null inGraphList then trace "No graphs to undergo Genetic Algorithm" []
   else
      trace ("Genetic Algorithm operating on population of " ++ show (length inGraphList) ++ " input graph(s) with cost range ("++ show (minimum $ fmap snd6 inGraphList) ++ "," ++ show (maximum $ fmap snd6 inGraphList) ++ ")") (

      -- process args
      let (doElitist, keepNum, popSize, generations, severity, recombinations, maxNetEdges, stopNum) = getGeneticAlgParams inArgs
          (newGraphList, generationCounter) = GA.geneticAlgorithm inGS inData rSeed doElitist (fromJust maxNetEdges) (fromJust keepNum) (fromJust popSize) (fromJust generations) 0 (fromJust severity) (fromJust recombinations) 0 stopNum inGraphList
      in
      trace ("\tGenetic Algorithm: " ++ show (length newGraphList) ++ " resulting graphs with cost range (" ++ show (minimum $ fmap snd6 newGraphList) ++ "," ++ show (maximum $ fmap snd6 newGraphList) ++ ")" ++ " after " ++ show generationCounter ++ " generation(s)")
      newGraphList
      )

-- | getGeneticAlgParams returns paramlist from arglist
getGeneticAlgParams :: [Argument] 
                    -> (Bool, Maybe Int, Maybe Int, Maybe Int, Maybe Double, Maybe Int, Maybe Int, Int)
getGeneticAlgParams inArgs =
      let fstArgList = fmap (fmap toLower . fst) inArgs
          sndArgList = fmap (fmap toLower . snd) inArgs
          lcArgList = zip fstArgList sndArgList
          checkCommandList = checkCommandArgs "geneticalgorithm" fstArgList VER.geneticAlgorithmArgList
      in
      -- check for valid command options
      if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'GeneticAlgorithm': " ++ show inArgs)
      else
         let keepList = filter ((=="keep").fst) lcArgList
             keepNum
              | length keepList > 1 =
                errorWithoutStackTrace ("Multiple 'keep' number specifications in fuse command--can have only one: " ++ show inArgs)
              | null keepList = Just 10
              | otherwise = readMaybe (snd $ head keepList) :: Maybe Int

             popSizeList = filter ((=="popsize").fst) lcArgList
             popSize
              | length popSizeList > 1 =
                errorWithoutStackTrace ("Multiple 'popsize' number specifications in genetic algorithm command--can have only one: " ++ show inArgs)
              | null popSizeList = Just 20
              | otherwise = readMaybe (snd $ head popSizeList) :: Maybe Int

             generationsList = filter ((=="generations").fst) lcArgList
             generations
              | length generationsList > 1 =
                errorWithoutStackTrace ("Multiple 'generations' number specifications in genetic algorithm command--can have only one: " ++ show inArgs)
              | null generationsList = Just 5
              | otherwise = readMaybe (snd $ head generationsList) :: Maybe Int

             severityList = filter ((=="severity").fst) lcArgList
             severity
              | length severityList > 1 =
                errorWithoutStackTrace ("Multiple 'severity' number specifications in genetic algorithm command--can have only one: " ++ show inArgs)
              | null severityList = Just 1.0
              | otherwise = readMaybe (snd $ head severityList) :: Maybe Double

             recombinationsList = filter ((=="recombinations").fst) lcArgList
             recombinations
              | length recombinationsList > 1 =
                errorWithoutStackTrace ("Multiple 'recombinations' number specifications in genetic algorithm command--can have only one: " ++ show inArgs)
              | null recombinationsList = Just 100
              | otherwise = readMaybe (snd $ head recombinationsList) :: Maybe Int

             maxNetEdgesList = filter ((=="maxnetedges").fst) lcArgList
             maxNetEdges
              | length maxNetEdgesList > 1 =
                errorWithoutStackTrace ("Multiple 'maxNetEdges' number specifications in genetic algorithm command--can have only one: " ++ show inArgs)
              | null maxNetEdgesList = Just 10
              | otherwise = readMaybe (snd $ head maxNetEdgesList) :: Maybe Int

             stopList = filter ((=="stop").fst) lcArgList
             stopNum
                | length stopList > 1 =
                  errorWithoutStackTrace ("Multiple 'stop' number specifications in search command--can have only one: " ++ show inArgs)
                | null stopList = Just (maxBound :: Int)
                | otherwise = readMaybe (snd $ head stopList) :: Maybe Int

             -- in case want to make it an option
             -- doElitist' = any ((=="nni").fst) lcArgList
             doElitist = True
         in

         --check arguments
         if isNothing keepNum then errorWithoutStackTrace ("Keep specification not an integer in Genetic Algorithm: "  ++ show (head keepList))
         else if isNothing popSize then errorWithoutStackTrace ("PopSize specification not an integer in Genetic Algorithm: "  ++ show (head popSizeList))
         else if isNothing generations then errorWithoutStackTrace ("Generations specification not an integer in Genetic Algorithm: "  ++ show (head generationsList))
         else if isNothing severity then errorWithoutStackTrace ("Severity factor specification not an integer in Genetic Algorithm: "  ++ show (head severityList))
         else if isNothing recombinations then errorWithoutStackTrace ("Severity factor specification not an integer in Genetic Algorithm: "  ++ show (head recombinationsList))
         else if isNothing stopNum then errorWithoutStackTrace ("Stop specification not an integer or not found in Genetic Algorithm (e.g. stop:10) "  ++ show (head stopList))
         else
            (doElitist, keepNum, popSize, generations, severity, recombinations, maxNetEdges, fromJust stopNum)

-- | fuseGraphs is a wrapper for graph recombination
-- the functions make heavy use of branch swapping functions in Search.Swap
fuseGraphs :: [Argument] 
           -> GlobalSettings 
           -> ProcessedData 
           -> Int 
           -> [PhylogeneticGraph] 
           -> [PhylogeneticGraph]
fuseGraphs inArgs inGS inData rSeed inGraphList
  | null inGraphList = trace "Fusing--skipped: No graphs to fuse" []
  | length inGraphList == 1 = trace "Fusing--skipped: Need > 1 graphs to fuse" inGraphList
  | graphType inGS == HardWired = trace "Fusing hardwired graphs is currenty not implemented" inGraphList
  | otherwise = trace ("Fusing " ++ show (length inGraphList) ++ " input graph(s) with minimum cost "++ show (minimum $ fmap snd6 inGraphList)) (

     -- process args for fuse placement
           let (keepNum, maxMoveEdgeDist, fusePairs, lcArgList) = getFuseGraphParams inArgs

               -- steepest off by default due to wanteing to check all addition points
               doSteepest' = any ((=="steepest").fst) lcArgList
               doAll = any ((=="all").fst) lcArgList

               doSteepest
                 | (not doSteepest' && not doAll) = True
                 | (doSteepest' && doAll) = True
                 | doAll = False
                 | otherwise = doSteepest'

               -- readdition options, specified as swap types
               -- no alternate or nni for fuse--not really meaningful

               swapType
                 | any ((=="tbr").fst) lcArgList = TBR
                 | any ((=="spr").fst) lcArgList = SPR
                 | otherwise = None

               -- turn off union selection of rejoin--default to do both, union first
               joinType
                 | any ((=="joinall").fst) lcArgList = JoinAll
                 | any ((=="joinpruned").fst) lcArgList = JoinPruned
                 | otherwise = JoinAlternate

               -- set implied alignment swapping
               doIA' = any ((=="ia").fst) lcArgList
               doIA = if (graphType inGS /= Tree) && doIA' then trace "\tIgnoring 'IA' swap option for non-Tree" False
                      else doIA'

               returnBest = any ((=="best").fst) lcArgList
               returnUnique = (not returnBest) || (any ((=="unique").fst) lcArgList)
               doSingleRound = any ((=="once").fst) lcArgList
               randomPairs = any ((=="atrandom").fst) lcArgList
               fusePairs' = if fusePairs == Just (maxBound :: Int) then Nothing
                            else fusePairs

               -- this for exchange or one dirction transfer of sub-graph--one half time for noreciprocal
               reciprocal' = any ((=="reciprocal").fst) lcArgList
               notReciprocal = any ((=="notreciprocal").fst) lcArgList
               reciprocal
                 | not reciprocal' = False
                 | notReciprocal = False
                 | otherwise = True

               seedList = randomIntList rSeed

               -- populate SwapParams structure
               swapParams = SwapParams {  swapType = swapType
                                             , joinType = joinType 
                                             , atRandom = randomPairs -- really same as swapping at random not so important here
                                             , keepNum  = (fromJust keepNum)
                                             , maxMoveEdgeDist = (2 * fromJust maxMoveEdgeDist)
                                             , steepest = doSteepest
                                             , joinAlternate = False -- join prune alternates--turned off for now
                                             , doIA = doIA
                                             , returnMutated = False -- no SA/Drift swapping in Fuse
                                             }
           in
           -- perform graph fuse operations
           let (newGraphList, counterFuse) = F.fuseAllGraphs swapParams inGS inData seedList 0 returnBest returnUnique doSingleRound fusePairs' randomPairs reciprocal inGraphList

           in
           trace ("\tAfter fusing: " ++ show (length newGraphList) ++ " resulting graphs with minimum cost " ++ show (minimum $ fmap snd6 newGraphList) ++ " after fuse rounds (total): " ++ show counterFuse)
           newGraphList
      )

-- | getFuseGraphParams returns fuse parameters from arglist
getFuseGraphParams :: [Argument] 
                   -> (Maybe Int, Maybe Int, Maybe Int, [(String, String)])
getFuseGraphParams inArgs =
    let fstArgList = fmap (fmap toLower . fst) inArgs
        sndArgList = fmap (fmap toLower . snd) inArgs
        lcArgList = zip fstArgList sndArgList
        checkCommandList = checkCommandArgs "fuse" fstArgList VER.fuseArgList
     in
     -- check for valid command options
     if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'fuse': " ++ show inArgs)
     else
         let keepList = filter ((=="keep").fst) lcArgList
             keepNum
              | length keepList > 1 =
                errorWithoutStackTrace ("Multiple 'keep' number specifications in fuse command--can have only one: " ++ show inArgs)
              | null keepList = Just 10
              | otherwise = readMaybe (snd $ head keepList) :: Maybe Int

             -- this is awkward syntax but works to chejkc fpor multiple swapping limit commands
             moveLimitList = filter (not . null) (snd <$> filter ((`notElem` ["keep", "pairs"]).fst) lcArgList)
             maxMoveEdgeDist
              | length moveLimitList > 1 =
                errorWithoutStackTrace ("Multiple maximum edge distance number specifications in fuse command--can have only one (e.g. spr:2): " ++ show inArgs)
              | null moveLimitList = Just ((maxBound :: Int) `div` 3)
              | otherwise = readMaybe (head moveLimitList) :: Maybe Int

             pairList = filter ((=="pairs").fst) lcArgList
             fusePairs
              | length pairList > 1 =
                errorWithoutStackTrace ("Multiple 'pair' number specifications in fuse command--can have only one: " ++ show inArgs)
              | null pairList = Just (maxBound :: Int)
              | otherwise = readMaybe (snd $ head pairList) :: Maybe Int

        in

        --check arguments
        if isNothing keepNum then errorWithoutStackTrace ("Keep specification not an integer in swap: "  ++ show (head keepList))
        else if isNothing maxMoveEdgeDist then errorWithoutStackTrace ("Maximum edge move distance specification in fuse command not an integer (e.g. spr:2): "  ++ show (head moveLimitList))
        else if isNothing fusePairs then errorWithoutStackTrace ("fusePairs specification not an integer in fuse: "  ++ show (head pairList))

        else
            (keepNum, maxMoveEdgeDist, fusePairs, lcArgList)

-- | netEdgeMaster overall master for add/delete net edges
netEdgeMaster :: [Argument] 
              -> GlobalSettings 
              -> ProcessedData 
              -> Int 
              -> [PhylogeneticGraph] 
              -> [PhylogeneticGraph]
netEdgeMaster inArgs inGS inData rSeed inGraphList =
   if null inGraphList then trace "No graphs to edit network edges" []
   else if graphType inGS == Tree then trace "\tCannot perform network edge operations on graphtype tree--set graphtype to SoftWired or HardWired" inGraphList

   -- process args for netEdgeMaster
   else
           let (keepNum, steps', annealingRounds', driftRounds', acceptEqualProb, acceptWorseFactor, maxChanges, maxNetEdges, lcArgList, maxRounds) = getNetEdgeParams inArgs
               doNetAdd = any ((=="netadd").fst) lcArgList
               doNetDelete = any ((=="netdel").fst) lcArgList || any ((=="netdelete").fst) lcArgList
               doAddDelete = any ((=="netadddel").fst) lcArgList || any ((=="netadddelete").fst) lcArgList
               doMove = any ((=="netmove").fst) lcArgList
               doSteepest' = any ((=="steepest").fst) lcArgList
               doAll = any ((=="all").fst) lcArgList

               -- do steepest default
               doSteepest
                 | (not doSteepest' && not doAll) = True
                 | doSteepest' && doAll = True
                 | otherwise = doSteepest'

               -- randomized order default
               doRandomOrder'
                 | any ((=="atrandom").fst) lcArgList = True
                 | any ((=="inorder").fst) lcArgList = False
                 | otherwise = True

               -- simulated annealing parameters
               -- returnMutated to return annealed Graphs before swapping fir use in Genetic Algorithm
               doAnnealing = any ((=="annealing").fst) lcArgList

               doDrift     = any ((=="drift").fst) lcArgList

               --ensures random edge order for drift/annealing
               doRandomOrder = doRandomOrder' || doDrift || doAnnealing

               returnMutated = any ((=="returnmutated").fst) lcArgList

               simAnnealParams = if not doAnnealing && not doDrift then Nothing
                                 else
                                    let steps = max 3 (fromJust steps')
                                        annealingRounds
                                          | isNothing annealingRounds' = 1
                                          | fromJust annealingRounds' < 1 = 1
                                          | otherwise = fromJust annealingRounds'

                                        driftRounds
                                          | isNothing driftRounds' = 1
                                          | fromJust driftRounds' < 1 = 1
                                          | otherwise = fromJust driftRounds'

                                        saMethod
                                          | doDrift && doAnnealing = trace "\tSpecified both Simulated Annealing (with temperature steps) and Drifting (without)--defaulting to drifting."
                                                    Drift
                                          | doDrift = Drift
                                          | otherwise = SimAnneal

                                        equalProb
                                          | fromJust acceptEqualProb < 0.0 = 0.0
                                          | fromJust acceptEqualProb > 1.0 = 1.0
                                          | otherwise = fromJust acceptEqualProb

                                        worseFactor = max (fromJust acceptWorseFactor) 0.0

                                        changes = if fromJust maxChanges < 0 then 15
                                                     else fromJust maxChanges

                                        saValues = SAParams { method = saMethod
                                                            , numberSteps = steps
                                                            , currentStep = 0
                                                            , randomIntegerList = randomIntList rSeed
                                                            , rounds      = max annealingRounds driftRounds
                                                            , driftAcceptEqual  = equalProb
                                                            , driftAcceptWorse  = worseFactor
                                                            , driftMaxChanges   = changes
                                                            , driftChanges      = 0
                                                            }
                                    in
                                    Just saValues


               -- create simulated annealing random lists uniquely for each fmap
               newSimAnnealParamList = U.generateUniqueRandList (length inGraphList) simAnnealParams

               -- perform add/delete/move operations
               bannerText
                 | isJust simAnnealParams = let editString
                                                  | doNetAdd = " add) "
                                                  | doNetDelete = " delete) "
                                                  | doAddDelete = " add/delete) "
                                                  | otherwise = " move "
                                            in
                                            if method (fromJust simAnnealParams) == SimAnneal then
                                               "Simulated Annealing (Network edge" ++ editString ++ show (rounds $ fromJust simAnnealParams) ++ " rounds " ++ show (length inGraphList) ++ " with " ++ show (numberSteps $ fromJust simAnnealParams) ++ " cooling steps " ++ show (length inGraphList) ++ " input graph(s) at minimum cost "++ show (minimum $ fmap snd6 inGraphList) ++ " keeping maximum of " ++ show (fromJust keepNum) ++ " graphs"
                                            else
                                               "Drifting (Network edge" ++ editString ++ show (rounds $ fromJust simAnnealParams) ++ " rounds " ++ show (length inGraphList) ++ " with " ++ show (numberSteps $ fromJust simAnnealParams) ++ " cooling steps " ++ show (length inGraphList) ++ " input graph(s) at minimum cost "++ show (minimum $ fmap snd6 inGraphList) ++ " keeping maximum of " ++ show (fromJust keepNum) ++ " graphs"
                 | doNetDelete = ("Network edge delete on " ++ show (length inGraphList) ++ " input graph(s) with minimum cost "++ show (minimum $ fmap snd6 inGraphList))
                 | doNetAdd = ("Network edge add on " ++ show (length inGraphList) ++ " input graph(s) with minimum cost "++ show (minimum $ fmap snd6 inGraphList) ++ " and maximum " ++ show (fromJust maxRounds) ++ " rounds")
                 | doAddDelete = ("Network edge add/delete on " ++ show (length inGraphList) ++ " input graph(s) with minimum cost "++ show (minimum $ fmap snd6 inGraphList)  ++ " and maximum " ++ show (fromJust maxRounds) ++ " rounds")
                 | doMove = ("Network edge move on " ++ show (length inGraphList) ++ " input graph(s) with minimum cost "++ show (minimum $ fmap snd6 inGraphList))
                 | otherwise = ""

            in
            trace bannerText (

            let (newGraphList, counterAdd) = if doNetAdd then
                                                if graphType inGS == HardWired then
                                                    trace "Adding edges to hardwired graphs will always increase cost, skipping"
                                                    (inGraphList, 0)
                                                else
                                                    -- trace ("REFINE Add") (
                                                    let graphPairList = PU.seqParMap (parStrategy $ strictParStrat inGS)  (N.insertAllNetEdges inGS inData rSeed (fromJust maxNetEdges) (fromJust keepNum) (fromJust maxRounds) 0 returnMutated doSteepest doRandomOrder ([], infinity)) (zip newSimAnnealParamList (fmap (: []) inGraphList)) -- `using` PU.myParListChunkRDS
                                                        (graphListList, counterList) = unzip graphPairList
                                                    in
                                                    (GO.selectGraphs Unique (fromJust keepNum) 0.0 (-1) $ concat graphListList, sum counterList)
                                                    -- )
                                             else (inGraphList, 0)


                (newGraphList', counterDelete) = if doNetDelete then
                                                    {-
                                                    if graphType inGS == HardWired then
                                                        trace ("Deleting edges from hardwired graphs will trivially remove all network edges to a tree, skipping")
                                                        (newGraphList, 0)
                                                    else
                                                    -}
                                                  -- trace ("REFINE Delete") (
                                                     let graphPairList = PU.seqParMap (parStrategy $ strictParStrat inGS)  (N.deleteAllNetEdges inGS inData rSeed (fromJust maxNetEdges) (fromJust keepNum) 0 returnMutated doSteepest doRandomOrder ([], infinity)) (zip newSimAnnealParamList (fmap (: []) newGraphList)) -- `using` PU.myParListChunkRDS
                                                         (graphListList, counterList) = unzip graphPairList
                                                     in
                                                     (GO.selectGraphs Unique (fromJust keepNum) 0.0 (-1) $ concat graphListList, sum counterList)
                                                  -- )
                                                else (newGraphList, 0)


                (newGraphList'', counterMove) = if doMove then
                                                    -- trace ("Network move option currently disabled--skipping.")
                                                    -- (newGraphList', 0 :: Int)

                                                    let graphPairList = PU.seqParMap (parStrategy $ strictParStrat inGS)  (N.moveAllNetEdges inGS inData rSeed (fromJust maxNetEdges) (fromJust keepNum) 0 returnMutated doSteepest doRandomOrder ([], infinity)) (zip newSimAnnealParamList (fmap (: []) newGraphList')) -- `using` PU.myParListChunkRDS
                                                        (graphListList, counterList) = unzip graphPairList
                                                    in
                                                    (GO.selectGraphs Unique (fromJust keepNum) 0.0 (-1) $ concat graphListList, sum counterList)

                                                else (newGraphList', 0)

                (newGraphList''', counterAddDelete) = if doAddDelete then
                                                        if graphType inGS == HardWired then
                                                            trace "Adding and Deleting edges to/from hardwired graphs will trivially remove all network edges to a tree, skipping"
                                                            (newGraphList'', 0)
                                                        else
                                                            let graphPairList = PU.seqParMap (parStrategy $ strictParStrat inGS)  (N.addDeleteNetEdges inGS inData rSeed (fromJust maxNetEdges) (fromJust keepNum) (fromJust maxRounds) 0 returnMutated doSteepest doRandomOrder ([], infinity)) (zip newSimAnnealParamList (fmap (: []) newGraphList'')) -- `using` PU.myParListChunkRDS
                                                                (graphListList, counterList) = unzip graphPairList
                                                            in
                                                            (GO.selectGraphs Unique (fromJust keepNum) 0.0 (-1) $ concat graphListList, sum counterList)

                                                else (newGraphList'', 0)

            in
            let resultGraphList = if null newGraphList''' then inGraphList
                                  else GO.selectGraphs Unique (fromJust keepNum) 0.0 (-1) newGraphList'''
            in
            trace ("\tAfter network edge add/delete/move: " ++ show (length resultGraphList) ++ " resulting graphs at cost " ++ show (minimum $ fmap snd6 resultGraphList) ++ " with add/delete/move rounds (total): " ++ show counterAdd ++ " Add, "
            ++ show counterDelete ++ " Delete, " ++ show counterMove ++ " Move, " ++ show counterAddDelete ++ " AddDelete")
            resultGraphList
            )

-- | getNetEdgeParams returns net edge cparameters from argument list
getNetEdgeParams :: [Argument]
                 -> (Maybe Int, Maybe Int, Maybe Int, Maybe Int, Maybe Double, Maybe Double, Maybe Int, Maybe Int, [(String, String)], Maybe Int)
getNetEdgeParams inArgs =
     let fstArgList = fmap (fmap toLower . fst) inArgs
         sndArgList = fmap (fmap toLower . snd) inArgs
         lcArgList = zip fstArgList sndArgList
         checkCommandList = checkCommandArgs "netEdgeMaster" fstArgList VER.netEdgeArgList
     in
     -- check for valid command options
     if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'netEdge': " ++ show inArgs)
     else
         let keepList = filter ((=="keep").fst) lcArgList
             keepNum
              | length keepList > 1 =
                errorWithoutStackTrace ("Multiple 'keep' number specifications in netEdge command--can have only one: " ++ show inArgs)
              | null keepList = Just 10
              | otherwise = readMaybe (snd $ head keepList) :: Maybe Int

             -- simulated anealing options
             stepsList   = filter ((=="steps").fst) lcArgList
             steps'
              | length stepsList > 1 =
                errorWithoutStackTrace ("Multiple annealing steps value specifications in netEdge command--can have only one (e.g. steps:10): " ++ show inArgs)
              | null stepsList = Just 10
              | otherwise = readMaybe (snd $ head stepsList) :: Maybe Int

             annealingList = filter ((=="annealing").fst) lcArgList
             annealingRounds'
              | length annealingList > 1 =
                errorWithoutStackTrace ("Multiple 'annealing' rounds number specifications in netEdge command--can have only one: " ++ show inArgs)
              | null annealingList = Just 1
              | otherwise = readMaybe (snd $ head annealingList) :: Maybe Int

             -- drift options
             driftList = filter ((=="drift").fst) lcArgList
             driftRounds'
              | length driftList > 1 =
                errorWithoutStackTrace ("Multiple 'drift' rounds number specifications in swap command--can have only one: " ++ show inArgs)
              | null driftList = Just 1
              | otherwise = readMaybe (snd $ head driftList) :: Maybe Int

             acceptEqualList = filter ((=="acceptequal").fst) lcArgList
             acceptEqualProb
              | length acceptEqualList > 1 =
                errorWithoutStackTrace ("Multiple 'drift' acceptEqual specifications in swap command--can have only one: " ++ show inArgs)
              | null acceptEqualList = Just 0.5
              | otherwise = readMaybe (snd $ head acceptEqualList) :: Maybe Double

             acceptWorseList = filter ((=="acceptworse").fst) lcArgList
             acceptWorseFactor
              | length acceptWorseList > 1 =
                errorWithoutStackTrace ("Multiple 'drift' acceptWorse specifications in swap command--can have only one: " ++ show inArgs)
              | null acceptWorseList = Just 20.0
              | otherwise = readMaybe (snd $ head acceptWorseList) :: Maybe Double

             maxChangesList = filter ((=="maxchanges").fst) lcArgList
             maxChanges
              | length maxChangesList > 1 =
                errorWithoutStackTrace ("Multiple 'drift' maxChanges number specifications in swap command--can have only one: " ++ show inArgs)
              | null maxChangesList = Just 15
              | otherwise = readMaybe (snd $ head maxChangesList) :: Maybe Int

             maxNetEdgesList = filter ((=="maxnetedges").fst) lcArgList
             maxNetEdges
              | length maxNetEdgesList > 1 =
                errorWithoutStackTrace ("Multiple 'maxNetEdges' number specifications in netEdge command--can have only one: " ++ show inArgs)
              | null maxNetEdgesList = Just 10
              | otherwise = readMaybe (snd $ head maxNetEdgesList) :: Maybe Int

             maxRoundsList = filter ((=="rounds").fst) lcArgList
             maxRounds
              | length maxRoundsList > 1 =
                errorWithoutStackTrace ("Multiple 'rounds' number specifications in netEdge command--can have only one: " ++ show inArgs)
              | null maxRoundsList = Just 1
              | otherwise = readMaybe (snd $ head maxRoundsList) :: Maybe Int

         in

         -- check inputs
         if isNothing keepNum then errorWithoutStackTrace ("Keep specification not an integer in netEdge: "  ++ show (head keepList))
         else if isNothing steps'           then errorWithoutStackTrace ("Annealing steps specification not an integer (e.g. steps:10): "  ++ show (snd $ head stepsList))
         else if isNothing acceptEqualProb  then errorWithoutStackTrace ("Drift 'acceptEqual' specification not a float (e.g. acceptEqual:0.75): "  ++ show (snd $ head acceptEqualList))
         else if isNothing acceptWorseFactor then errorWithoutStackTrace ("Drift 'acceptWorse' specification not a float (e.g. acceptWorse:1.0): "  ++ show (snd $ head acceptWorseList))
         else if isNothing maxChanges       then errorWithoutStackTrace ("Drift 'maxChanges' specification not an integer (e.g. maxChanges:10): "  ++ show (snd $ head maxChangesList))
         else if isNothing maxNetEdges       then errorWithoutStackTrace ("Drift 'maxChanges' specification not an integer (e.g. maxChanges:10): "  ++ show (snd $ head maxNetEdgesList))
         else if isNothing maxRounds       then errorWithoutStackTrace ("Network edit 'rounds' specification not an integer (e.g. rounds:10): "  ++ show (snd $ head maxRoundsList))

         else
            (keepNum, steps', annealingRounds', driftRounds', acceptEqualProb, acceptWorseFactor, maxChanges, maxNetEdges, lcArgList, maxRounds)
