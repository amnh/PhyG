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
              | null keepList = Just 1
              | otherwise = readMaybe (snd $ head keepList) :: Maybe Int
             moveLimitList = filter (not . null) $ fmap snd $ filter ((/="keep").fst) lcArgList
             maxMoveEdgeDist  
              | length moveLimitList > 1 =
                errorWithoutStackTrace ("Multiple maximum edge distance number specifications in fuse command--can have only one (e.g. spr:2): " ++ show inArgs)
              | null moveLimitList = Just ((maxBound :: Int) `div` 2) 
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
           let (newGraphList, counterFuse) = F.fuseAllGraphs inGS inData seed (fromJust keepNum) doNNI doSPR doTBR doSteepest doAll returnBest returnUnique doSingleRound inGraphList
           in                             
           trace ("After graph fusing : " ++ (show $ length newGraphList) ++ " resulting graphs with fuse rounds (total): " ++ (show counterFuse))
           newGraphList
      )

-- | driver for overall refinement
refineGraph :: [Argument] -> GlobalSettings -> ProcessedData -> Int -> [PhylogeneticGraph] -> [PhylogeneticGraph]
refineGraph inArgs inGS inData seed inGraphList = 
   netEdgeMaster inArgs inGS inData seed inGraphList

-- | refinement arguments
refineArgList :: [String]
refineArgList = ["spr","tbr", "keep", "steepest", "all", "nni", "ia" ,"netadd", "netdel", "netdelete", "nad","netmove"]

-- | netEdgeArgList argiments for network edge add/delete operations
netEdgeArgList :: [String]
netEdgeArgList = ["keep", "steepest", "all", "netadd", "netdel", "netdelete", "nad", "netmove"]

-- | netEdgeMaster overall master for add/delete net edges
netEdgeMaster :: [Argument] -> GlobalSettings -> ProcessedData -> Int -> [PhylogeneticGraph] -> [PhylogeneticGraph]
netEdgeMaster inArgs inGS inData rSeed inGraphList =
   if graphType inGS == Tree then trace ("Cannot perform network edge operations on graphtype tree--set graphtype to SoftWired or HardWired") inGraphList
   else 
      trace ("Network edge add/delete on " ++ (show $ length inGraphList) ++ " input graph(s) with minimum cost "++ (show $ minimum $ fmap snd6 inGraphList)) (
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
              | null keepList = Just 1
              | otherwise = readMaybe (snd $ head keepList) :: Maybe Int
         in
         if isNothing keepNum then errorWithoutStackTrace ("Keep specification not an integer in netEdge: "  ++ show (head keepList))
         else 
           let -- getting values to be passed for graph diagnosis later
               numLeaves = V.length $ fst3 inData
               leafGraph = T.makeSimpleLeafGraph inData
               leafDecGraph = T.makeLeafGraph inData
               leafGraphSoftWired = T.makeLeafGraphSoftWired inData
               hasNonExactChars = U.getNumberNonExactCharacters (thd3 inData) > 0
               charInfoVV = six6 $ head inGraphList

               -- process args for netEdgeMaster
               doNetAdd = any ((=="netadd").fst) lcArgList
               doNetDelete = (any ((=="netdel").fst) lcArgList) || (any ((=="netdelete").fst) lcArgList)
               doAddDelete = any ((=="nad").fst) lcArgList
               doMove = any ((=="netmove").fst) lcArgList
               doSteepest' = any ((=="steepest").fst) lcArgList
               doAll = any ((=="all").fst) lcArgList
               doSteepest = if (not doSteepest' && not doAll) then True
                            else doSteepest'

               -- performo add/delete operations 
               (newGraphList, counterDelete) = if doNetDelete || doAddDelete then 
                                                let graphPairList = fmap (N.deleteAllNetEdges inGS inData (fromJust keepNum) 0 ([], infinity)) (fmap (: []) inGraphList) `using` PU.myParListChunkRDS
                                                    (graphListList, counterList) = unzip graphPairList
                                                in (concat graphListList, sum counterList)
                                            else (inGraphList, 0)   


               (newGraphList', counterAdd) = if doNetAdd || doAddDelete then 
                                                let graphPairList = fmap (N.insertAllNetEdges inGS inData (fromJust keepNum) 0 ([], infinity)) (fmap (: []) newGraphList) `using` PU.myParListChunkRDS
                                                    (graphListList, counterList) = unzip graphPairList
                                                in (concat graphListList, sum counterList)
                                            else (newGraphList, 0)   

               (newGraphList'', counterMove) = if doMove then 
                                                let graphPairList = fmap (N.moveAllNetEdges inGS inData (fromJust keepNum) 0 ([], infinity)) (fmap (: []) newGraphList) `using` PU.myParListChunkRDS
                                                    (graphListList, counterList) = unzip graphPairList
                                                in (concat graphListList, sum counterList)
                                            else (newGraphList', 0)   

           in
           trace ("After network edge add/delete/move: " ++ (show $ length newGraphList'') ++ " resulting graphs with add/delete/move rounds (total): " ++ (show counterAdd) ++ " Add, " 
            ++ (show counterDelete) ++ " Delete, " ++ (show counterMove) ++ " Move")
           newGraphList'
     )
