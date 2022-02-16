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
                          , refineArgList
                          , netEdgeMaster
                          , fuseGraphs
                          ) where

import Types.Types
import qualified ParallelUtilities       as PU
import Control.Parallel.Strategies
import qualified Utilities.LocalGraph    as LG
import Debug.Trace
import qualified GraphOptimization.Traversals as T
import qualified GraphOptimization.PreOrderFunctions as PRE
import qualified GraphOptimization.PostOrderFunctions as POS
import qualified Data.List as L
import qualified Data.Text.Lazy              as TL
import GeneralUtilities
import qualified Graphs.GraphOperations  as GO
import Data.Maybe
import qualified Search.Swap as S
import qualified Search.NetworkAddDelete as N
import Utilities.Utilities               as U
import qualified Data.Vector as V
import Data.Maybe
import           Data.Char
import           Text.Read
import qualified Search.Fuse as F

-- | fuseArgList arguments
fuseArgList :: [String]
fuseArgList = ["spr","tbr", "keep", "steepest", "all", "nni", "best", "unique", "once"]

-- | fuseGraphs is a wrapper for graph recombination
-- the functions make heavy use of branch swapping functions in Search.Swap 
fuseGraphs :: [Argument] -> GlobalSettings -> ProcessedData -> Int -> [PhylogeneticGraph] -> [PhylogeneticGraph]
fuseGraphs inArgs inGS inData seed inGraphList =
   if null inGraphList then []
   else 
      trace ("Fusing " ++ (show $ length inGraphList) ++ " input graph(s) with minimum cost "++ (show $ minimum $ fmap snd6 inGraphList)) (
      let fstArgList = fmap (fmap toLower . fst) inArgs
          sndArgList = fmap (fmap toLower . snd) inArgs
          lcArgList = zip fstArgList sndArgList
          checkCommandList = U.checkCommandArgs "fuse" fstArgList fuseArgList
     in
     -- check for valid command options
     if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'swap': " ++ show inArgs)
     else 
         let keepList = filter ((=="keep").fst) lcArgList
             keepNum
              | length keepList > 1 =
                errorWithoutStackTrace ("Multiple 'keep' number specifications in fuse command--can have only one: " ++ show inArgs)
              | null keepList = Just 10
              | otherwise = readMaybe (snd $ head keepList) :: Maybe Int
             moveLimitList = filter (not . null) $ fmap snd $ filter ((/="keep").fst) lcArgList
             maxMoveEdgeDist  
              | length moveLimitList > 1 =
                errorWithoutStackTrace ("Multiple maximum edge distance number specifications in fuse command--can have only one (e.g. spr:2): " ++ show inArgs)
              | null moveLimitList = Just ((maxBound :: Int) `div` 3) 
              | otherwise = readMaybe (head moveLimitList) :: Maybe Int
        in
        if isNothing keepNum then errorWithoutStackTrace ("Keep specification not an integer in swap: "  ++ show (head keepList))
        else if isNothing maxMoveEdgeDist then errorWithoutStackTrace ("Maximum edge move distance specification in fuse command not an integer (e.g. spr:2): "  ++ show (snd $ head keepList))
        else 
           let -- process args for fuse placement
               doNNI' = any ((=="nni").fst) lcArgList
               doSPR' = any ((=="spr").fst) lcArgList
               doTBR = any ((=="tbr").fst) lcArgList
               doSteepest' = any ((=="steepest").fst) lcArgList
               doAll = any ((=="all").fst) lcArgList
               doSteepest = if (not doSteepest' && not doAll) then True
                            else doSteepest'
               doSPR = if doTBR then False
                       else doSPR'
               doNNI = if doSPR || doTBR then False
                       else doNNI'
               returnBest = any ((=="best").fst) lcArgList
               returnUnique = any ((=="unique").fst) lcArgList
               doSingleRound = any ((=="once").fst) lcArgList
           in
           -- perform graph fuse operations 
           let (newGraphList, counterFuse) = F.fuseAllGraphs inGS inData seed (fromJust keepNum) (2 * (fromJust maxMoveEdgeDist)) 0 doNNI doSPR doTBR doSteepest doAll returnBest returnUnique doSingleRound inGraphList
           in                             
           trace ("\tAfter fusing: " ++ (show $ length newGraphList) ++ " resulting graphs with minimum cost " ++ (show $ minimum $ fmap snd6 newGraphList) ++ " after fuse rounds (total): " ++ (show counterFuse))
           newGraphList
      )

-- | driver for overall refinement
refineGraph :: [Argument] -> GlobalSettings -> ProcessedData -> Int -> [PhylogeneticGraph] -> [PhylogeneticGraph]
refineGraph inArgs inGS inData seed inGraphList = 
   if null inGraphList then error ("No graphs input to refine") 
   else netEdgeMaster inArgs inGS inData seed inGraphList

-- | refinement arguments
refineArgList :: [String]
refineArgList = ["spr","tbr", "keep", "steepest", "all", "nni", "ia" ,"netadd", "netdel", "netdelete", "nad","netmove"]

-- | netEdgeArgList argiments for network edge add/delete operations
netEdgeArgList :: [String]
netEdgeArgList = ["keep", "steepest", "all", "netadd", "netdel", "netdelete", "nad", "netmove", "annealing", "steps", "returnmutated", "drift", "acceptequal", "acceptworse", "maxchanges","steepest","atrandom"]

-- | netEdgeMaster overall master for add/delete net edges
netEdgeMaster :: [Argument] -> GlobalSettings -> ProcessedData -> Int -> [PhylogeneticGraph] -> [PhylogeneticGraph]
netEdgeMaster inArgs inGS inData rSeed inGraphList =
   if null inGraphList then error "No graphs input to netEdgeMaster"
   else if graphType inGS == Tree then trace ("\tCannot perform network edge operations on graphtype tree--set graphtype to SoftWired or HardWired") inGraphList
   else 
      let fstArgList = fmap (fmap toLower . fst) inArgs
          sndArgList = fmap (fmap toLower . snd) inArgs
          lcArgList = zip fstArgList sndArgList
          checkCommandList = U.checkCommandArgs "netEdgeMaster" fstArgList netEdgeArgList
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
             doAnnealing = any ((=="annealing").fst) lcArgList

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
         if isNothing keepNum then errorWithoutStackTrace ("Keep specification not an integer in netEdge: "  ++ show (head keepList))
         else if isNothing steps'           then errorWithoutStackTrace ("Annealing steps specification not an integer (e.g. steps:10): "  ++ show (snd $ head stepsList))
         else if isNothing acceptEqualProb  then errorWithoutStackTrace ("Drift 'acceptEqual' specification not a float (e.g. acceptEqual:0.75): "  ++ show (snd $ head acceptEqualList))
         else if isNothing acceptWorseFactor then errorWithoutStackTrace ("Drift 'acceptWorse' specification not a float (e.g. acceptWorse:1.0): "  ++ show (snd $ head acceptWorseList))
         else if isNothing maxChanges       then errorWithoutStackTrace ("Drift 'maxChanges' specification not an integer (e.g. maxChanges:10): "  ++ show (snd $ head maxChangesList))
         else 
           let -- getting values to be passed for graph diagnosis later
               numLeaves = V.length $ fst3 inData
               leafGraph = T.makeSimpleLeafGraph inData
               leafDecGraph = T.makeLeafGraph inData
               leafGraphSoftWired = T.makeLeafGraphSoftWired inData
               hasNonExactChars = U.getNumberNonExactCharacters (thd3 inData) > 0
               charInfoVV = fmap thd3 $ thd3 inData -- six6 $ head inGraphList

               -- process args for netEdgeMaster
               doNetAdd = any ((=="netadd").fst) lcArgList
               doNetDelete = (any ((=="netdel").fst) lcArgList) || (any ((=="netdelete").fst) lcArgList)
               doAddDelete = any ((=="nad").fst) lcArgList
               doMove = any ((=="netmove").fst) lcArgList
               doSteepest' = any ((=="steepest").fst) lcArgList
               doAll = any ((=="all").fst) lcArgList

               -- do steepest default
               doSteepest = if (not doSteepest' && not doAll) then True
                            else if doSteepest' && doAll then True
                            else doSteepest'
               doRandomOrder = any ((=="atRandom").fst) lcArgList

               -- simulated annealing parameters
               -- returnMutated to return annealed Graphs before swapping fir use in Genetic Algorithm
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

               -- perform add/delete/move operations
               bannerText = if simAnnealParams /= Nothing then 
                              if (method $ fromJust simAnnealParams) == SimAnneal then
                                 ("Simulated Annealing (Network edge moves) " ++ (show $ rounds $ fromJust simAnnealParams) ++ " rounds " ++ (show $ length inGraphList) ++ " with " ++ (show $ numberSteps $ fromJust simAnnealParams) ++ " cooling steps " ++ (show $ length inGraphList) ++ " input graph(s) at minimum cost "++ (show $ minimum $ fmap snd6 inGraphList) ++ " keeping maximum of " ++ (show $ fromJust keepNum) ++ " graphs")
                              else 
                                 ("Drifting (Network edge moves) " ++ (show $ rounds $ fromJust simAnnealParams) ++ " rounds " ++ (show $ length inGraphList) ++ " with " ++ (show $ numberSteps $ fromJust simAnnealParams) ++ " cooling steps " ++ (show $ length inGraphList) ++ " input graph(s) at minimum cost "++ (show $ minimum $ fmap snd6 inGraphList) ++ " keeping maximum of " ++ (show $ fromJust keepNum) ++ " graphs")
                            else if doNetDelete || doAddDelete then 
                              ("Network edge delete on " ++ (show $ length inGraphList) ++ " input graph(s) with minimum cost "++ (show $ minimum $ fmap snd6 inGraphList))
                            else if doNetAdd || doAddDelete then 
                              ("Network edge add on " ++ (show $ length inGraphList) ++ " input graph(s) with minimum cost "++ (show $ minimum $ fmap snd6 inGraphList))
                            else if doMove then
                              ("Network edge move on " ++ (show $ length inGraphList) ++ " input graph(s) with minimum cost "++ (show $ minimum $ fmap snd6 inGraphList))
                            else ""

            in
            trace (bannerText) (
               
            let (newGraphList, counterDelete) = if doNetDelete || doAddDelete then 
                                                let graphPairList = fmap (N.deleteAllNetEdges inGS inData rSeed (fromJust keepNum) 0 returnMutated doSteepest doRandomOrder ([], infinity)) (zip newSimAnnealParamList (fmap (: []) inGraphList)) `using` PU.myParListChunkRDS
                                                    (graphListList, counterList) = unzip graphPairList
                                                in (GO.selectPhylogeneticGraph [("unique", (show $ fromJust keepNum))] 0 ["unique"] $ concat graphListList, sum counterList)
                                                
                                               else (inGraphList, 0)   


                (newGraphList', counterAdd) = if doNetAdd || doAddDelete then 
                                               let graphPairList = fmap (N.insertAllNetEdges inGS inData rSeed (fromJust keepNum) 0 returnMutated doSteepest doRandomOrder ([], infinity)) (zip newSimAnnealParamList (fmap (: []) inGraphList)) `using` PU.myParListChunkRDS
                                                   (graphListList, counterList) = unzip graphPairList
                                                in (GO.selectPhylogeneticGraph [("unique", (show $ fromJust keepNum))] 0 ["unique"] $ concat graphListList, sum counterList)
                                                
                                             else (newGraphList, 0)   

                (newGraphList'', counterMove) = if doMove then 
                                                let graphPairList = fmap (N.moveAllNetEdges inGS inData rSeed (fromJust keepNum) 0 returnMutated doSteepest doRandomOrder ([], infinity)) (zip newSimAnnealParamList (fmap (: []) inGraphList)) `using` PU.myParListChunkRDS
                                                    (graphListList, counterList) = unzip graphPairList
                                                in (GO.selectPhylogeneticGraph [("unique", (show $ fromJust keepNum))] 0 ["unique"] $ concat graphListList, sum counterList)
                                                
                                             else (newGraphList', 0)   

            in
            let resultGraphList = if null newGraphList'' then inGraphList
                                  else GO.selectPhylogeneticGraph [("unique", (show $ fromJust keepNum))] 0 ["unique"] $ inGraphList ++ newGraphList''
            in
            trace ("\tAfter network edge add/delete/move: " ++ (show $ length resultGraphList) ++ " resulting graphs at cost " ++ (show $ minimum $ fmap snd6 resultGraphList) ++ " with add/delete/move rounds (total): " ++ (show counterAdd) ++ " Add, " 
            ++ (show counterDelete) ++ " Delete, " ++ (show counterMove) ++ " Move")
            resultGraphList
            )
     
