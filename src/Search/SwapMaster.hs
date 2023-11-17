{- |
Module      :  SwapMaster.hs
Description :  Module controlling grapjh swap functions
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

module Search.SwapMaster  ( swapMaster
                          ) where

import qualified Commands.Verify             as VER
import           Data.Char
import           Data.Maybe
import           Debug.Trace
import           GeneralUtilities
import qualified Graphs.GraphOperations      as GO
import qualified ParallelUtilities           as PU
import qualified Search.Swap                 as S
import           Text.Read
import           Types.Types
import           Utilities.Utilities         as U

-- | swapMaster processes and spawns the swap functions
-- the 2 x maxMoveDist since distance either side to list 2* dist on sorted edges
swapMaster ::  [Argument] 
           -> GlobalSettings
           -> ProcessedData 
           -> Int 
           -> [ReducedPhylogeneticGraph] 
           -> [ReducedPhylogeneticGraph]
swapMaster inArgs inGS inData rSeed inGraphListInput =
   if null inGraphListInput then trace "No graphs to swap" []
   -- else if graphType inGS == HardWired then trace ("Swapping hardwired graphs is currenty not implemented") inGraphList
   else

           let -- process args for swap
               (keepNum, maxMoveEdgeDist', steps', annealingRounds', doDrift, driftRounds', acceptEqualProb, acceptWorseFactor, maxChanges, replicateNumber, lcArgList) = getSwapParams inArgs

               swapType
                 | any ((=="nni").fst) lcArgList = NNI
                 | any ((=="spr").fst) lcArgList = SPR
                 | any ((=="tbr").fst) lcArgList = TBR
                 | any ((=="alternate").fst) lcArgList = Alternate
                 | otherwise = Alternate

               maxMoveEdgeDist = if swapType == NNI then 2
                                  else fromJust maxMoveEdgeDist'

               -- randomized orders of split and join-- not implemented
               -- doRandomized = any ((=="randomized").fst) lcArgList

               -- set implied alignment swapping
               doIA' = any ((=="ia").fst) lcArgList
               doIA'' = if (graphType inGS /= Tree) && doIA' then trace "\tIgnoring 'IA' swap option for non-Tree" False
                      else doIA'

               --- steepest/all options
               doSteepest' = any ((=="steepest").fst) lcArgList
               doAll = any ((=="all").fst) lcArgList

               -- steepest default
               doSteepest = ((not doSteepest' && not doAll) || doSteepest')

               -- simulated annealing parameters
               -- returnMutated to return annealed Graphs before swapping fir use in Genetic Algorithm
               doAnnealing = any ((=="annealing").fst) lcArgList

               returnMutated = any ((=="returnmutated").fst) lcArgList

               -- turn off union selection of rejoin--default to do both, union first
               joinType
                 | graphType inGS == HardWired = JoinAll
                 | any ((=="joinall").fst) lcArgList = JoinAll
                 | any ((=="joinpruned").fst) lcArgList = JoinPruned
                 | any ((=="joinalternate").fst) lcArgList = JoinAlternate
                 | otherwise = JoinAlternate 


               -- randomize split graph and rejoin edges, defualt to randomize
               atRandom
                 | any ((=="atrandom").fst) lcArgList = True
                 | any ((=="inOrder").fst) lcArgList = False
                 | otherwise = True

               randomIntListSwap = randomIntList rSeed

               simAnnealParams = getSimAnnealParams doAnnealing doDrift steps' annealingRounds' driftRounds' acceptEqualProb acceptWorseFactor maxChanges rSeed

               -- swap replicates is meant to allow multiple randomized swap trajectories
               -- set to 1 if not randomized swap or SA/Drifting (set by their own options)
               replicates
                 | not atRandom = 1
                 | doAnnealing = 1
                 | otherwise = fromJust replicateNumber

               -- replicate inGraphList based on 'replicates' for randomized trajectories
               inGraphList  = concat $ replicate replicates inGraphListInput
               numGraphs = length inGraphList

               -- create simulated annealing random lists uniquely for each fmap
               newSimAnnealParamList = U.generateUniqueRandList numGraphs simAnnealParams

               progressString
                 | (not doAnnealing && not doDrift) = ("Swapping " <> show (length inGraphListInput) <> " input graph(s) with " <> show replicates <> " trajectories at minimum cost " <> show (minimum $ fmap snd5 inGraphList) <> " keeping maximum of " <> show (fromJust keepNum) <> " graphs per input graph")
                 | method (fromJust simAnnealParams) == SimAnneal = ("Simulated Annealing (Swapping) " <> show (rounds $ fromJust simAnnealParams) <> " rounds with " <> show (numberSteps $ fromJust simAnnealParams) <> " cooling steps " <> show (length inGraphList) <> " input graph(s) at minimum cost " <> show (minimum $ fmap snd5 inGraphList) <> " keeping maximum of " <> show (fromJust keepNum) <> " graphs")
                 | otherwise = "Drifting (Swapping) " <> show (rounds $ fromJust simAnnealParams) <> " rounds with " <> show (driftMaxChanges $ fromJust simAnnealParams) <> " maximum changes per round on " <> show (length inGraphList) <> " input graph(s) at minimum cost " <> show (minimum $ fmap snd5 inGraphList) <> " keeping maximum of " <> show (fromJust keepNum) <> " graphs"

               -- populate SwapParams structure
               localSwapParams = SwapParams {  swapType = swapType
                                             , joinType = joinType 
                                             , atRandom = atRandom
                                             , keepNum  = (fromJust keepNum)
                                             , maxMoveEdgeDist = maxMoveEdgeDist
                                             , steepest = doSteepest
                                             , joinAlternate = False -- join prune alternates--turned off for now
                                             , doIA = doIA''
                                             , returnMutated = returnMutated 
                                             }
           in

           trace progressString (
           let (newGraphList, counter) = let graphPairList = PU.seqParMap (parStrategy $ defaultParStrat inGS) (S.swapSPRTBR localSwapParams inGS inData 0 inGraphList) ((:[]) <$> zip3 (U.generateRandIntLists (head randomIntListSwap) numGraphs) newSimAnnealParamList inGraphList)

                                             -- graphPairList = PU.seqParMap PU.myStrategyHighLevel  (S.swapSPRTBR localSwapParams inGS inData 0 inGraphList) ((:[]) <$> zip3 (U.generateRandIntLists (head randomIntListSwap) numGraphs) newSimAnnealParamList inGraphList) -- `using` PU.myParListChunkRDS
                                             (graphListList, counterList) = unzip graphPairList
                                         in (GO.selectGraphs Best (fromJust keepNum) 0.0 (-1) $ concat graphListList, sum counterList)
              in
              let finalGraphList = if null newGraphList then inGraphList
                                   else newGraphList

                  fullBuffWarning = if length newGraphList >= (fromJust keepNum) then "\n\tWarning--Swap returned as many minimum cost graphs as the 'keep' number.  \n\tThis may have limited the effectiveness of the swap. \n\tConsider increasing the 'keep' value or adding an additional swap."
                                    else ""

                  endString
                    | (not doAnnealing && not doDrift) = ("\n\tAfter swap: " <> show (length finalGraphList) <> " resulting graphs with minimum cost " <> show (minimum $ fmap snd5 finalGraphList) <> " with swap rounds (total): " <> show counter <> " " <> show swapType)
                    | method (fromJust simAnnealParams) == SimAnneal = ("\n\tAfter Simulated Annealing: " <> show (length finalGraphList) <> " resulting graphs with minimum cost " <> show (minimum $ fmap snd5 finalGraphList) <> " with swap rounds (total): " <> show counter <> " " <> show swapType)
                    | otherwise = "\n\tAfter Drifting: " <> show (length finalGraphList) <> " resulting graphs with minimum cost " <> show (minimum $ fmap snd5 finalGraphList) <> " with swap rounds (total): " <> show counter <> " " <> show swapType

              in
              trace (endString <> fullBuffWarning)
              finalGraphList
              )


-- | getSimumlatedAnnealingParams returns SA parameters
getSimAnnealParams :: Bool 
                   -> Bool 
                   -> Maybe Int 
                   -> Maybe Int 
                   -> Maybe Int
                   -> Maybe Double 
                   -> Maybe Double 
                   -> Maybe Int 
                   -> Int 
                   -> Maybe SAParams
getSimAnnealParams doAnnealing doDrift steps' annealingRounds' driftRounds' acceptEqualProb acceptWorseFactor maxChanges rSeed =
    if not doAnnealing && not doDrift then Nothing
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

-- | getSwapParams takes areg list and preocesses returning parameter values
getSwapParams :: [Argument] 
              -> (Maybe Int, Maybe Int, Maybe Int, Maybe Int, Bool, Maybe Int, Maybe Double, Maybe Double, Maybe Int, Maybe Int, [(String, String)])
getSwapParams inArgs =
    let fstArgList = fmap (fmap toLower . fst) inArgs
        sndArgList = fmap (fmap toLower . snd) inArgs
        lcArgList = zip fstArgList sndArgList
        checkCommandList = checkCommandArgs "swap" fstArgList VER.swapArgList
    in
     -- check for valid command options
     if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'swap': " <> show inArgs)
     else
         let keepList = filter ((=="keep").fst) lcArgList
             keepNum
              | length keepList > 1 =
                errorWithoutStackTrace ("Multiple 'keep' number specifications in swap command--can have only one: " <> show inArgs)
              | null keepList = Just 10
              | otherwise = readMaybe (snd $ head keepList) :: Maybe Int

             moveLimitList = filter (not . null) (snd <$> filter ((`elem` ["alternate", "spr", "tbr", "nni"]).fst) lcArgList)
             maxMoveEdgeDist'
              | length moveLimitList > 1 =
                errorWithoutStackTrace ("Multiple maximum edge distance number specifications in swap command--can have only one (e.g. spr:2): " <> show inArgs)
              | null moveLimitList = Just ((maxBound :: Int) `div` 3)
              | otherwise = readMaybe (head moveLimitList) :: Maybe Int

             -- simulated anealing options
             stepsList   = filter ((=="steps").fst) lcArgList
             steps'
              | length stepsList > 1 =
                errorWithoutStackTrace ("Multiple annealing steps value specifications in swap command--can have only one (e.g. steps:10): " <> show inArgs)
              | null stepsList = Just 10
              | otherwise = readMaybe (snd $ head stepsList) :: Maybe Int

             annealingList = filter ((=="annealing").fst) lcArgList
             annealingRounds'
              | length annealingList > 1 =
                errorWithoutStackTrace ("Multiple 'annealing' rounds number specifications in swap command--can have only one: " <> show inArgs)
              | null annealingList = Just 1
              | otherwise = readMaybe (snd $ head annealingList) :: Maybe Int

             -- drift options
             doDrift   = any ((=="drift").fst) lcArgList

             driftList = filter ((=="drift").fst) lcArgList
             driftRounds'
              | length driftList > 1 =
                errorWithoutStackTrace ("Multiple 'drift' rounds number specifications in swap command--can have only one: " <> show inArgs)
              | null driftList = Just 1
              | otherwise = readMaybe (snd $ head driftList) :: Maybe Int

             acceptEqualList = filter ((=="acceptequal").fst) lcArgList
             acceptEqualProb
              | length acceptEqualList > 1 =
                errorWithoutStackTrace ("Multiple 'drift' acceptEqual specifications in swap command--can have only one: " <> show inArgs)
              | null acceptEqualList = Just 0.5
              | otherwise = readMaybe (snd $ head acceptEqualList) :: Maybe Double

             acceptWorseList = filter ((=="acceptworse").fst) lcArgList
             acceptWorseFactor
              | length acceptWorseList > 1 =
                errorWithoutStackTrace ("Multiple 'drift' acceptWorse specifications in swap command--can have only one: " <> show inArgs)
              | null acceptWorseList = Just 20.0
              | otherwise = readMaybe (snd $ head acceptWorseList) :: Maybe Double

             maxChangesList = filter ((=="maxchanges").fst) lcArgList
             maxChanges
              | length maxChangesList > 1 =
                errorWithoutStackTrace ("Multiple 'drift' maxChanges number specifications in swap command--can have only one: " <> show inArgs)
              | null maxChangesList = Just 15
              | otherwise = readMaybe (snd $ head maxChangesList) :: Maybe Int

             replicatesList = filter ((=="replicates").fst) lcArgList
             replicates
              | length replicatesList > 1 =
                errorWithoutStackTrace ("Multiple 'swap' replicates number specifications in swap command--can have only one: " <> show inArgs)
              | null replicatesList = Just 1
              | otherwise = readMaybe (snd $ head replicatesList) :: Maybe Int

        in
        -- check inputs
        if isNothing keepNum               then errorWithoutStackTrace ("Keep specification not an integer in swap: "  <> show (head keepList))
        else if isNothing maxMoveEdgeDist' then errorWithoutStackTrace ("Maximum edge move distance specification not an integer (e.g. spr:2): "  <> show (head moveLimitList))
        else if isNothing steps'           then errorWithoutStackTrace ("Annealing steps specification not an integer (e.g. steps:10): "  <> show (snd $ head stepsList))
        else if isNothing acceptEqualProb  then errorWithoutStackTrace ("Drift 'acceptEqual' specification not a float (e.g. acceptEqual:0.75): "  <> show (snd $ head acceptEqualList))
        else if isNothing acceptWorseFactor then errorWithoutStackTrace ("Drift 'acceptWorse' specification not a float (e.g. acceptWorse:1.0): "  <> show (snd $ head acceptWorseList))
        else if isNothing maxChanges       then errorWithoutStackTrace ("Drift 'maxChanges' specification not an integer (e.g. maxChanges:10): "  <> show (snd $ head maxChangesList))
        else if isNothing replicates       then errorWithoutStackTrace ("Swap 'replicates' specification not an integer (e.g. replicates:5): "  <> show (snd $ head replicatesList))

        else
            -- trace ("GSP: " <> (show inArgs) <> " " <> (show )(keepNum, maxMoveEdgeDist', steps', annealingRounds', doDrift, driftRounds', acceptEqualProb, acceptWorseFactor, maxChanges, lcArgList))
            (keepNum, maxMoveEdgeDist', steps', annealingRounds', doDrift, driftRounds', acceptEqualProb, acceptWorseFactor, maxChanges, replicates, lcArgList)
