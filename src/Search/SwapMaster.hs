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

import Types.Types
import qualified ParallelUtilities       as PU
import Control.Parallel.Strategies
import Debug.Trace
import GeneralUtilities
import qualified Graphs.GraphOperations  as GO
import qualified Search.Swap as S
import qualified Search.NetworkAddDelete as N
import Utilities.Utilities               as U
import Data.Maybe
import           Data.Char
import           Text.Read
import qualified Search.Fuse as F
import qualified Search.GeneticAlgorithm as GA
-- | buildArgList is the list of valid build arguments
swapArgList :: [String]
swapArgList = ["spr","tbr", "keep", "steepest", "all", "nni", "ia", "annealing", "maxtemp", "mintemp", "steps", "returnmutated", "drift", "acceptequal", "acceptworse", "maxchanges"]


-- | swapMaster processes and spawns the swap functions
-- the 2 x maxMoveDist since distance either side to list 2* dist on sorted edges
swapMaster ::  [Argument] -> GlobalSettings -> ProcessedData -> Int -> [PhylogeneticGraph] -> [PhylogeneticGraph]
swapMaster inArgs inGS inData rSeed inGraphList = 
   if null inGraphList then trace ("No graphs to swap") []
   else 
      let fstArgList = fmap (fmap toLower . fst) inArgs
          sndArgList = fmap (fmap toLower . snd) inArgs
          lcArgList = zip fstArgList sndArgList
          checkCommandList = checkCommandArgs "swap" fstArgList swapArgList
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
           let -- process args for swap
               doNNI' = any ((=="nni").fst) lcArgList
               doSPR' = any ((=="spr").fst) lcArgList
               doTBR' = any ((=="tbr").fst) lcArgList
               doIA' = any ((=="ia").fst) lcArgList
               doIA = if (graphType inGS /= Tree) then trace ("\tIgnoring 'IA' swap option for non-Tree") False
                      else doIA'
               doSteepest' = any ((=="steepest").fst) lcArgList
               doAll = any ((=="all").fst) lcArgList
               doSPR'' = if (not doNNI' && not doSPR' && not doTBR') then True
                         else doSPR'
               
               doSteepest = if (not doSteepest' && not doAll) then True
                            else doSteepest'
               
               -- Workaround for Hardwired SPR issue
               (doTBR, doSPR, doNNI, maxMoveEdgeDist, hardWiredSPR) = if (graphType inGS /= HardWired) then (doTBR', doSPR'', doNNI', maxMoveEdgeDist', False)
                                                                      else
                                                                        if doNNI' then (True, False, False, Just 1, True)
                                                                        else if doSPR'' then (True, False, False, maxMoveEdgeDist', True)
                                                                        else (doTBR', doSPR'', doNNI', maxMoveEdgeDist', False)

               -- simulated annealing parameters
               -- returnMutated to return annealed Graphs before swapping fir use in Genetic Algorithm
               doAnnealing = any ((=="annealing").fst) lcArgList
             
               returnMutated = any ((=="returnmutated").fst) lcArgList

               simAnnealParams = if (not doAnnealing && not doDrift) then Nothing
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
                                               let graphPairList1 = fmap (S.swapSPRTBR "nni" inGS inData (fromJust keepNum) 2 doSteepest hardWiredSPR doIA returnMutated) (zip newSimAnnealParamList inGraphList) `using` PU.myParListChunkRDS
                                                   (graphListList, counterList) = unzip graphPairList1
                                               in (take (fromJust keepNum) $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] $ concat graphListList, sum counterList)
                                             else (inGraphList, 0)
               (newGraphList', counterSPR)  = if doSPR then 
                                               let graphPairList2 = fmap (S.swapSPRTBR "spr" inGS inData (fromJust keepNum) (2 * (fromJust maxMoveEdgeDist)) doSteepest hardWiredSPR doIA returnMutated) (zip newSimAnnealParamList newGraphList) `using` PU.myParListChunkRDS
                                                   (graphListList, counterList) = unzip graphPairList2
                                               in 
                                              (take (fromJust keepNum) $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] $ concat graphListList, sum counterList)
                                             else (newGraphList, 0)

               (newGraphList'', counterTBR) = if doTBR then 
                                               let graphPairList3 =  fmap (S.swapSPRTBR "tbr" inGS inData (fromJust keepNum) (2 * (fromJust maxMoveEdgeDist)) doSteepest hardWiredSPR doIA returnMutated) (zip newSimAnnealParamList newGraphList') `using` PU.myParListChunkRDS
                                                   (graphListList, counterList) = unzip graphPairList3
                                               in 
                                               (take (fromJust keepNum) $ GO.selectPhylogeneticGraph [("unique", "")] 0 ["unique"] $ concat graphListList, sum counterList)
                                             else (newGraphList', 0)
              in
              let endString = if not doAnnealing then ("\tAfter swap: " ++ (show $ length newGraphList'') ++ " resulting graphs with swap rounds (total): " ++ (show counterNNI) ++ " NNI, " ++ (show counterSPR) ++ " SPR, " ++ (show counterTBR) ++ " TBR")
                              else if (method $ fromJust simAnnealParams) == SimAnneal then
                                ("\tAfter Simulated Annealing: " ++ (show $ length newGraphList'') ++ " resulting graphs") 
                              else 
                                ("\tAfter Drifting: " ++ (show $ length newGraphList'') ++ " resulting graphs")
              in 
              trace (endString)
              newGraphList''     
            )