{- |
Module      :  Build.hs
Description :  Module specifying graph building functions
Copyright   :  (c) 2021 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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

{-
Distance builds (Wagner, NJ, WPGMA are imported from Wag2020.
https://www.github.com/wardwheeler/wag2020
-}

module Search.Build  ( buildGraph
                     ) where

import qualified GraphOptimization.Traversals as T
import qualified Data.Text.Lazy              as TL
import Types.Types
import qualified ParallelUtilities            as PU
import Utilities.Utilities as U
import qualified Utilities.DistanceUtilities as DU
import qualified SymMatrix               as M
import           Data.Maybe
import           Text.Read
import           Data.Char
import qualified Data.List               as L
import qualified Data.Vector             as V
import qualified Search.DistanceMethods  as DM
import qualified Search.WagnerBuild      as WB
import GeneralUtilities
import qualified Graphs.GraphOperations  as GO
import qualified Search.DistanceWagner   as DW
import Debug.Trace

-- | buildArgList is the list of valid build arguments
buildArgList :: [String]
buildArgList = ["replicates", "nj", "wpgma", "dwag", "rdwag", "distance", "character", "tree", "graph", "best","none","otu","spr","tbr"]

-- | buildGraph takes build options and returns contructed graphList 
buildGraph :: [Argument] -> GlobalSettings -> ProcessedData ->  [[VertexCost]] -> Int-> [PhylogeneticGraph]
buildGraph inArgs inGS inData@(nameTextVect, _, _) pairwiseDistances seed =
   let fstArgList = fmap (fmap toLower . fst) inArgs
       sndArgList = fmap (fmap toLower . snd) inArgs
       lcArgList = zip fstArgList sndArgList
       checkCommandList = U.checkCommandArgs "build" fstArgList buildArgList
   in
   -- check for valid command options
   if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'build': " ++ show inArgs)
   else

      let buildDistance = filter ((=="distance").fst) lcArgList
          buildCharacter = filter ((=="character").fst) lcArgList
          repPairList = filter ((=="replicates").fst) lcArgList
          numReplicates
            | length repPairList > 1 =
              errorWithoutStackTrace ("Multiple replicate number specifications in command--can have only one: " ++ show inArgs)
            | null repPairList = Just 10
            | otherwise = readMaybe (snd $ head repPairList) :: Maybe Int
          keepPairList = filter ((=="best").fst) lcArgList
          numToSave
            | length keepPairList > 1 =
                  errorWithoutStackTrace ("Multiple best number specifications in command--can have only one: " ++ show inArgs)
            | null keepPairList = numReplicates
            | otherwise = readMaybe (snd $ head keepPairList) :: Maybe Int

      in
      if not (null buildDistance) && not (null buildCharacter) then
         errorWithoutStackTrace ("Cannot specify both 'character' and 'distance' builds in same build command" ++ show inArgs)
      else if isNothing numReplicates then errorWithoutStackTrace ("Replicates specification not an integer: "  ++ show (snd $ head repPairList))
      else if isNothing numToSave then errorWithoutStackTrace ("Best specification not an integer: "  ++ show (snd $ head keepPairList))
      else if not $ null buildDistance then
         -- distance build
         -- do all options in line and add together for return tree list
         let doDWag  = any ((=="dwag").fst) lcArgList
             doRDwag = any ((=="rdwag").fst) lcArgList
             doNJ    = any ((=="nj").fst) lcArgList
             doWPGMA = any ((=="wpgma").fst) lcArgList
             doOTU   = any ((=="otu").fst) lcArgList
             doSPR   = any ((=="spr").fst) lcArgList
             doTBR   = any ((=="tbr").fst) lcArgList
             outgroupElem = outgroupIndex inGS
             nameStringVect = fmap TL.unpack nameTextVect
             distMatrix = M.fromLists pairwiseDistances
         in
         let refinement
               | doTBR = "tbr"
               | doSPR = "spr"
               | doOTU = "otu"
               | otherwise = "none"
             treeList    = [distanceWagner inGS inData nameStringVect distMatrix outgroupElem refinement | doDWag]
             treeList'   = if doRDwag then treeList ++ randomizedDistanceWagner inGS inData nameStringVect distMatrix outgroupElem (fromJust numReplicates) seed (fromJust numToSave) refinement
                           else treeList
             treeList''  = if doNJ then treeList' ++ [neighborJoin inGS inData nameStringVect distMatrix outgroupElem refinement]
                           else treeList'
             treeList''' = if doWPGMA then treeList'' ++ [wPGMA inGS inData nameStringVect distMatrix outgroupElem refinement]
                           else treeList''
         in
         if null treeList''' then errorWithoutStackTrace ("Distance build is specified, but without any method: " ++ show inArgs)
         else
            -- trace (show inArgs ++ " Yielded " ++ (show $ length treeList''') ++ " trees")
            treeList'''

      else 
         -- character build 
         trace ("Building Character Wager") (
         let treeList = WB.rasWagnerBuild inGS inData seed (fromJust numReplicates) 
             graphList = fmap (T.multiTraverseFullyLabelGraph inGS inData False False)  (fmap fst6 treeList)
         in
         if (graphType inGS == SoftWired) then 
            -- multi-traverse reroots graphs
            graphList
         else if (graphType inGS == Tree) then 
            treeList
         else errorWithoutStackTrace ("Graph type " ++ (show $ graphType inGS) ++ " not implemented")
         )

-- | distanceWagner takes Processed data and pairwise distance matrix and returns
-- 'best' addition sequence Wagner (defined in Farris, 1972) as fully decorated tree (as Graph)
distanceWagner :: GlobalSettings -> ProcessedData -> V.Vector String -> M.Matrix Double -> Int -> String -> PhylogeneticGraph
distanceWagner inGS inData leafNames distMatrix outgroupValue refinement =
   let distWagTree = head $ DM.doWagnerS leafNames distMatrix "closest" outgroupValue "best" []
       distWagTree' = head $ DW.performRefinement refinement "best:1" "first" leafNames outgroupValue distWagTree
       distWagTreeSimpleGraph = DU.convertToDirectedGraphText leafNames outgroupValue (snd4 distWagTree')
   in
   T.multiTraverseFullyLabelTree inGS inData (GO.dichotomizeRoot outgroupValue $ GO.switchRootTree (length leafNames) distWagTreeSimpleGraph)

-- | randomizedDistanceWagner takes Processed data and pairwise distance matrix and returns
-- random addition sequence Wagner trees fully decorated tree (as Graph)
randomizedDistanceWagner :: GlobalSettings -> ProcessedData ->  V.Vector String -> M.Matrix Double -> Int -> Int -> Int -> Int -> String -> [PhylogeneticGraph]
randomizedDistanceWagner inGS inData leafNames distMatrix outgroupValue numReplicates seed numToKeep refinement =
   let randomizedAdditionSequences = V.fromList <$> shuffleInt seed numReplicates [0..(length leafNames - 1)]
       randomizedAdditionWagnerTreeList = DM.doWagnerS leafNames distMatrix "random" outgroupValue "random" randomizedAdditionSequences
       randomizedAdditionWagnerTreeList' = take numToKeep $ L.sortOn thd4 randomizedAdditionWagnerTreeList
       randomizedAdditionWagnerTreeList'' = head <$> PU.seqParMap PU.myStrategy (DW.performRefinement refinement "best:1"  "first" leafNames outgroupValue) randomizedAdditionWagnerTreeList'
       randomizedAdditionWagnerSimpleGraphList = fmap (DU.convertToDirectedGraphText leafNames outgroupValue . snd4) randomizedAdditionWagnerTreeList''
   in
   fmap ((T.multiTraverseFullyLabelTree inGS inData . GO.dichotomizeRoot outgroupValue) . GO.switchRootTree (length leafNames)) randomizedAdditionWagnerSimpleGraphList

-- | neighborJoin takes Processed data and pairwise distance matrix and returns
-- Neighbor-Joining tree as fully decorated tree (as Graph)
neighborJoin :: GlobalSettings -> ProcessedData -> V.Vector String -> M.Matrix Double -> Int -> String -> PhylogeneticGraph
neighborJoin inGS inData leafNames distMatrix outgroupValue refinement =
   let njTree = DM.neighborJoining leafNames distMatrix outgroupValue
       njTree' = head $ DW.performRefinement refinement "best:1" "first" leafNames outgroupValue njTree
       njSimpleGraph = DU.convertToDirectedGraphText leafNames outgroupValue (snd4 njTree')
   in
   T.multiTraverseFullyLabelTree inGS inData (GO.dichotomizeRoot outgroupValue $ GO.switchRootTree (length leafNames) njSimpleGraph)

-- | wPGMA takes Processed data and pairwise distance matrix and returns
-- WPGMA tree as fully decorated tree (as Graph)
-- since root index not nOTUs as with other tres--chanegd as with dWag and NJ to make consistent.
wPGMA :: GlobalSettings -> ProcessedData -> V.Vector String -> M.Matrix Double -> Int -> String -> PhylogeneticGraph
wPGMA inGS inData leafNames distMatrix outgroupValue refinement =
   let wpgmaTree = DM.wPGMA leafNames distMatrix outgroupValue
       wpgmaTree' = head $ DW.performRefinement refinement "best:1" "first" leafNames outgroupValue wpgmaTree
       wpgmaSimpleGraph = DU.convertToDirectedGraphText leafNames outgroupValue (snd4 wpgmaTree')
   in
   T.multiTraverseFullyLabelTree inGS inData (GO.dichotomizeRoot outgroupValue $ GO.switchRootTree (length leafNames) wpgmaSimpleGraph)
