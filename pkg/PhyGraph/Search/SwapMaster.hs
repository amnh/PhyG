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
import           Control.Parallel.Strategies
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
swapMaster ::  [Argument] -> GlobalSettings -> ProcessedData -> Int -> [PhylogeneticGraph] -> [PhylogeneticGraph]
swapMaster inArgs inGS inData rSeed inGraphList =
   if null inGraphList then trace ("No graphs to swap") []
   else

           let -- process args for swap
               (keepNum, maxMoveEdgeDist', steps', annealingRounds', doDrift, driftRounds', acceptEqualProb, acceptWorseFactor, maxChanges, lcArgList) = getSwapParams inArgs
               doNNI' = any ((=="nni").fst) lcArgList
               doSPR' = any ((=="spr").fst) lcArgList
               doTBR' = any ((=="tbr").fst) lcArgList

               -- randomized orders of split and join-- not implemented
               -- doRandomized = any ((=="randomized").fst) lcArgList

               doIA' = any ((=="ia").fst) lcArgList
               doIA = if (graphType inGS /= Tree) && doIA' then trace ("\tIgnoring 'IA' swap option for non-Tree") False
                      else doIA'
               doSteepest' = any ((=="steepest").fst) lcArgList
               doAll = any ((=="all").fst) lcArgList
               
               -- doSPR'' = if (not doNNI' && not doSPR' && not doTBR') then True
               --          else doSPR'

               -- steepest default
               doSteepest = if (not doSteepest' && not doAll) then True
                            else doSteepest'

               -- Workaround for Hardwired SPR issue
               -- remove if possible--very confusing
               (doTBR, doSPR, doNNI, maxMoveEdgeDist, hardWiredSPR) = if (graphType inGS /= HardWired) then (doTBR', doSPR', doNNI', maxMoveEdgeDist', False)
                                                                      else
                                                                        if doNNI' then (True, False, False, Just 1, True)
                                                                        else if doSPR then (True, False, False, maxMoveEdgeDist', True)
                                                                        else (doTBR', doSPR', doNNI', maxMoveEdgeDist', False)

               -- alternating rounds of SPR ande TBR is default unless NNI, SPR, or TBR specified 
               -- swap type if alternate will be "TBR"
               doAlternate = if (doTBR || doSPR || doNNI) then False
                            else True

               

               -- simulated annealing parameters
               -- returnMutated to return annealed Graphs before swapping fir use in Genetic Algorithm
               doAnnealing = any ((=="annealing").fst) lcArgList

               returnMutated = any ((=="returnmutated").fst) lcArgList

               simAnnealParams = getSimAnnealParams doAnnealing doDrift steps' annealingRounds' driftRounds' acceptEqualProb acceptWorseFactor maxChanges rSeed

               -- create simulated annealing random lists uniquely for each fmap
               newSimAnnealParamList = U.generateUniqueRandList (length inGraphList) simAnnealParams

               progressString = if not doAnnealing then ("Swapping " ++ (show $ length inGraphList) ++ " input graph(s) with minimum cost "++ (show $ minimum $ fmap snd6 inGraphList) ++ " keeping maximum of " ++ (show $ fromJust keepNum) ++ " graphs")
                             else
                                if (method $ fromJust simAnnealParams) == SimAnneal then
                                    ("Simulated Annealing (Swapping) " ++ (show $ rounds $ fromJust simAnnealParams) ++ " rounds " ++ (show $ length inGraphList) ++ " with " ++ (show $ numberSteps $ fromJust simAnnealParams) ++ " cooling steps " ++ (show $ length inGraphList) ++ " input graph(s) at minimum cost "++ (show $ minimum $ fmap snd6 inGraphList) ++ " keeping maximum of " ++ (show $ fromJust keepNum) ++ " graphs")
                                else
                                    ("Drifting (Swapping) " ++ (show $ rounds $ fromJust simAnnealParams) ++ " rounds " ++ (show $ length inGraphList) ++ " with " ++ (show $ numberSteps $ fromJust simAnnealParams) ++ " cooling steps " ++ (show $ length inGraphList) ++ " input graph(s) at minimum cost "++ (show $ minimum $ fmap snd6 inGraphList) ++ " keeping maximum of " ++ (show $ fromJust keepNum) ++ " graphs")

           in

           trace (progressString) (

           let (newGraphList, counterNNI)  = if doNNI then
                                               let graphPairList = PU.seqParMap rdeepseq (S.swapSPRTBR "nni" inGS inData (fromJust keepNum) 2 doSteepest False hardWiredSPR doIA returnMutated) (zip newSimAnnealParamList inGraphList) -- `using` PU.myParListChunkRDS
                                                   (graphListList, counterList) = unzip graphPairList
                                               in (take (fromJust keepNum) $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] $ concat graphListList, sum counterList)
                                             else (inGraphList, 0)
               (newGraphList', counterSPR')  = if (doSPR || doAlternate) then
                                               let graphPairList = PU.seqParMap rdeepseq  (S.swapSPRTBR "spr" inGS inData (fromJust keepNum) (2 * (fromJust maxMoveEdgeDist)) doSteepest False hardWiredSPR doIA returnMutated) (zip newSimAnnealParamList newGraphList) -- `using` PU.myParListChunkRDS
                                                   (graphListList, counterList) = unzip graphPairList
                                               in
                                               (take (fromJust keepNum) $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] $ concat graphListList, sum counterList)
                                             else (newGraphList, 0)

               (newGraphList'', counterTBR') = if (doTBR || doAlternate) then
                                               let graphPairList =  PU.seqParMap rdeepseq  (S.swapSPRTBR "tbr" inGS inData (fromJust keepNum) (2 * (fromJust maxMoveEdgeDist)) doSteepest doAlternate hardWiredSPR doIA returnMutated) (zip newSimAnnealParamList newGraphList') -- `using` PU.myParListChunkRDS
                                                   (graphListList, counterList) = unzip graphPairList
                                               in
                                               (take (fromJust keepNum) $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] $ concat graphListList, sum counterList)
                                             else (newGraphList', 0)
              in
              --if nothing better or equal
              let (counterSPR, counterTBR, counterAlternate) = if doAlternate then (0,0, counterSPR' + counterTBR')
                                                               else (counterSPR', counterTBR', 0)
                  finalGraphList = if null newGraphList'' then inGraphList
                                   else newGraphList''
                  endString = if not doAnnealing then ("\n\tAfter swap: " ++ (show $ length finalGraphList) ++ " resulting graphs with minimum cost " ++ (show $ minimum $ fmap snd6 finalGraphList) ++ " with swap rounds (total): " ++ (show counterNNI) ++ " NNI, " ++ (show counterSPR) ++ " SPR, " ++ (show counterTBR) ++ " TBR, " ++ (show counterAlternate) ++ " Alternating SPR/TBR")
                              else if (method $ fromJust simAnnealParams) == SimAnneal then
                                ("\n\tAfter Simulated Annealing: " ++ (show $ length finalGraphList) ++ " resulting graphs")
                              else
                                ("\n\tAfter Drifting: " ++ (show $ length finalGraphList) ++ " resulting graphs")
              in
              trace (endString)
              finalGraphList
            )

-- | getSimumlatedAnnealingParams returns SA parameters
getSimAnnealParams :: Bool -> Bool -> Maybe Int -> Maybe Int -> Maybe Int -> Maybe Double -> Maybe Double -> Maybe Int -> Int -> Maybe SAParams
getSimAnnealParams doAnnealing doDrift steps' annealingRounds' driftRounds' acceptEqualProb acceptWorseFactor maxChanges rSeed =
    if (not doAnnealing && not doDrift) then Nothing
                                 else
                                    let steps = max 3 (fromJust steps')
                                        annealingRounds = if annealingRounds' == Nothing then 1
                                                          else if fromJust annealingRounds' < 1 then 1
                                                          else fromJust annealingRounds'

                                        driftRounds = if driftRounds' == Nothing then 1
                                                          else if fromJust driftRounds' < 1 then 1
                                                          else fromJust driftRounds'

                                        saMethod = if doDrift && doAnnealing then
                                                    trace ("\tSpecified both Simulated Annealing (with temperature steps) and Drifting (without)--defaulting to drifting.")
                                                    Drift
                                                 else if doDrift then Drift
                                                 else SimAnneal

                                        equalProb = if fromJust acceptEqualProb < 0.0 then 0.0
                                                    else if fromJust acceptEqualProb > 1.0 then 1.0
                                                    else fromJust acceptEqualProb


                                        worseFactor = if fromJust acceptWorseFactor < 0.0 then 0.0
                                                      else fromJust acceptWorseFactor

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
getSwapParams :: [Argument] -> (Maybe Int, Maybe Int, Maybe Int, Maybe Int, Bool, Maybe Int, Maybe Double, Maybe Double, Maybe Int, [(String, String)])
getSwapParams inArgs =
    let fstArgList = fmap (fmap toLower . fst) inArgs
        sndArgList = fmap (fmap toLower . snd) inArgs
        lcArgList = zip fstArgList sndArgList
        checkCommandList = checkCommandArgs "swap" fstArgList VER.swapArgList
    in
     -- check for valid command options
     if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'swap': " ++ show inArgs)
     else
         let keepList = filter ((=="keep").fst) lcArgList
             keepNum
              | length keepList > 1 =
                errorWithoutStackTrace ("Multiple 'keep' number specifications in swap command--can have only one: " ++ show inArgs)
              | null keepList = Just 10
              | otherwise = readMaybe (snd $ head keepList) :: Maybe Int

             moveLimitList = filter (not . null) $ fmap snd $ filter ((`elem` ["spr", "tbr", "nni"]).fst) lcArgList
             maxMoveEdgeDist'
              | length moveLimitList > 1 =
                errorWithoutStackTrace ("Multiple maximum edge distance number specifications in swap command--can have only one (e.g. spr:2): " ++ show inArgs)
              | null moveLimitList = Just ((maxBound :: Int) `div` 3)
              | otherwise = readMaybe (head moveLimitList) :: Maybe Int

             -- simulated anealing options
             stepsList   = filter ((=="steps").fst) lcArgList
             steps'
              | length stepsList > 1 =
                errorWithoutStackTrace ("Multiple annealing steps value specifications in swap command--can have only one (e.g. steps:10): " ++ show inArgs)
              | null stepsList = Just 10
              | otherwise = readMaybe (snd $ head stepsList) :: Maybe Int

             annealingList = filter ((=="annealing").fst) lcArgList
             annealingRounds'
              | length annealingList > 1 =
                errorWithoutStackTrace ("Multiple 'annealing' rounds number specifications in swap command--can have only one: " ++ show inArgs)
              | null annealingList = Just 1
              | otherwise = readMaybe (snd $ head annealingList) :: Maybe Int

             -- drift options
             doDrift     = any ((=="drift").fst) lcArgList

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
              | null acceptWorseList = Just 1.0
              | otherwise = readMaybe (snd $ head acceptWorseList) :: Maybe Double

             maxChangesList = filter ((=="maxchanges").fst) lcArgList
             maxChanges
              | length maxChangesList > 1 =
                errorWithoutStackTrace ("Multiple 'drift' maxChanges number specifications in swap command--can have only one: " ++ show inArgs)
              | null maxChangesList = Just 15
              | otherwise = readMaybe (snd $ head maxChangesList) :: Maybe Int

        in
        -- check inputs
        if isNothing keepNum               then errorWithoutStackTrace ("Keep specification not an integer in swap: "  ++ show (head keepList))
        else if isNothing maxMoveEdgeDist' then errorWithoutStackTrace ("Maximum edge move distance specification not an integer (e.g. spr:2): "  ++ show (head moveLimitList))
        else if isNothing steps'           then errorWithoutStackTrace ("Annealing steps specification not an integer (e.g. steps:10): "  ++ show (snd $ head stepsList))
        else if isNothing acceptEqualProb  then errorWithoutStackTrace ("Drift 'acceptEqual' specification not a float (e.g. acceptEqual:0.75): "  ++ show (snd $ head acceptEqualList))
        else if isNothing acceptWorseFactor then errorWithoutStackTrace ("Drift 'acceptWorse' specification not a float (e.g. acceptWorse:1.0): "  ++ show (snd $ head acceptWorseList))
        else if isNothing maxChanges       then errorWithoutStackTrace ("Drift 'maxChanges' specification not an integer (e.g. maxChanges:10): "  ++ show (snd $ head maxChangesList))

        else
            (keepNum, maxMoveEdgeDist', steps', annealingRounds', doDrift, driftRounds', acceptEqualProb, acceptWorseFactor, maxChanges, lcArgList)
