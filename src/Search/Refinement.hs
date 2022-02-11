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
   if null inGraphList then trace ("No graphs input to refine") []
   else netEdgeMaster inArgs inGS inData seed inGraphList

-- | refinement arguments
refineArgList :: [String]
refineArgList = ["spr","tbr", "keep", "steepest", "all", "nni", "ia" ,"netadd", "netdel", "netdelete", "nad","netmove"]

-- | netEdgeArgList argiments for network edge add/delete operations
netEdgeArgList :: [String]
netEdgeArgList = ["keep", "steepest", "all", "netadd", "netdel", "netdelete", "nad", "netmove", "annealing", "maxtemp", "mintemp", "steps", "returnmutated"]

-- | netEdgeMaster overall master for add/delete net edges
netEdgeMaster :: [Argument] -> GlobalSettings -> ProcessedData -> Int -> [PhylogeneticGraph] -> [PhylogeneticGraph]
netEdgeMaster inArgs inGS inData rSeed inGraphList =
   if graphType inGS == Tree then trace ("\tCannot perform network edge operations on graphtype tree--set graphtype to SoftWired or HardWired") inGraphList
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
             maxTempList = filter ((=="maxtemp").fst) lcArgList
             maxTemp'     
              | length maxTempList > 1 =
                errorWithoutStackTrace ("Multiple maximum annealing temperature value specifications in netEdge command--can have only one (e.g. maxTemp:100): " ++ show inArgs)
              | null maxTempList = Just 11.0 
              | otherwise = readMaybe (snd $ head maxTempList) :: Maybe Double

             minTempList = filter ((=="mintemp").fst) lcArgList  
             minTemp'   
              | length minTempList > 1 =
                errorWithoutStackTrace ("Multiple minimum annealing temperature value specifications in netEdge command--can have only one (e.g. minTemp:1.0): " ++ show inArgs)
              | null minTempList = Just 1.0 
              | otherwise = readMaybe (snd $ head minTempList) :: Maybe Double

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

             doAnnealing = any ((=="annealing").fst) lcArgList
         in
         if isNothing keepNum then errorWithoutStackTrace ("Keep specification not an integer in netEdge: "  ++ show (head keepList))
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
               doSteepest = if (not doSteepest' && not doAll) then True
                            else doSteepest'

               -- simulated annealing parameters
               -- returnMutated to return annealed Graphs before swapping fir use in Genetic Algorithm
               returnMutated = any ((=="returnmutated").fst) lcArgList

               simAnnealParams = if not doAnnealing then Nothing
                                 else 
                                    let steps = max 3 (fromJust steps')
                                        maxTemp = max 20.1 (fromJust maxTemp')
                                        minTemp'' = max 0.1 (fromJust minTemp')
                                        minTemp = if minTemp'' > maxTemp then maxTemp / 100.0
                                                  else minTemp''
                                        annealingRounds = if annealingRounds' == Nothing then 1
                                                          else if fromJust annealingRounds' < 1 then 1
                                                          else fromJust annealingRounds'
                                    in
                                    Just (maxTemp, minTemp, steps, 0, randomIntList rSeed, annealingRounds) 

                                 
               -- create simulated annealing random lists uniquely for each fmap
               newSimAnnealParamList = U.generateUniqueRandList (length inGraphList) simAnnealParams

               -- perform add/delete/move operations
               bannerText = if simAnnealParams /= Nothing then 
                     ("Simulated Annealing (Network edge moves) " ++ (show $ six6 $ fromJust simAnnealParams) ++ " rounds " ++ (show $ length inGraphList) ++ " with " ++ (show $ thd6 $ fromJust simAnnealParams) ++ " cooling steps " ++ (show $ length inGraphList) ++ " input graph(s) at minimum cost "++ (show $ minimum $ fmap snd6 inGraphList) ++ " keeping maximum of " ++ (show $ fromJust keepNum) ++ " graphs")
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
                                                let graphPairList = fmap (N.deleteAllNetEdges inGS inData (fromJust keepNum) 0 returnMutated ([], infinity)) (zip newSimAnnealParamList (fmap (: []) inGraphList)) `using` PU.myParListChunkRDS
                                                    (graphListList, counterList) = unzip graphPairList
                                                in (GO.selectPhylogeneticGraph [("unique", (show $ fromJust keepNum))] 0 ["unique"] $ concat graphListList, sum counterList)
                                                
                                               else (inGraphList, 0)   


                (newGraphList', counterAdd) = if doNetAdd || doAddDelete then 
                                               let graphPairList = fmap (N.insertAllNetEdges inGS inData (fromJust keepNum) 0 returnMutated ([], infinity)) (zip newSimAnnealParamList (fmap (: []) inGraphList)) `using` PU.myParListChunkRDS
                                                   (graphListList, counterList) = unzip graphPairList
                                                in (GO.selectPhylogeneticGraph [("unique", (show $ fromJust keepNum))] 0 ["unique"] $ concat graphListList, sum counterList)
                                                
                                             else (newGraphList, 0)   

                (newGraphList'', counterMove) = if doMove then 
                                                let graphPairList = fmap (N.moveAllNetEdges inGS inData (fromJust keepNum) 0 returnMutated ([], infinity)) (zip newSimAnnealParamList (fmap (: []) inGraphList)) `using` PU.myParListChunkRDS
                                                    (graphListList, counterList) = unzip graphPairList
                                                in (GO.selectPhylogeneticGraph [("unique", (show $ fromJust keepNum))] 0 ["unique"] $ concat graphListList, sum counterList)
                                                
                                             else (newGraphList', 0)   

            in
            trace ("\tAfter network edge add/delete/move: " ++ (show $ length newGraphList'') ++ " resulting graphs at cost " ++ (show $ minimum $ fmap snd6 newGraphList'') ++ " with add/delete/move rounds (total): " ++ (show counterAdd) ++ " Add, " 
            ++ (show counterDelete) ++ " Delete, " ++ (show counterMove) ++ " Move")
            newGraphList''
            )
     
