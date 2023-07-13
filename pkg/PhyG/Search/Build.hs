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

import qualified Commands.Verify                as VER
import           Data.Char
import qualified Data.List                      as L
import qualified Data.Text.Lazy                 as TL
import qualified Data.Vector                    as V
import           Debug.Trace
import           GeneralUtilities
import qualified GraphOptimization.Traversals   as T
import qualified Graphs.GraphOperations         as GO
import qualified ParallelUtilities              as PU
import qualified Reconciliation.ReconcileGraphs as R
import qualified Search.DistanceMethods         as DM
import qualified Search.DistanceWagner          as DW
import qualified Search.WagnerBuild             as WB
import qualified SymMatrix                      as M
import           Text.Read
import           Types.Types
import qualified Utilities.DistanceUtilities    as DU
import qualified Utilities.Distances            as DD
import qualified Utilities.LocalGraph           as LG
import qualified Utilities.Utilities            as U
import           Data.Maybe

-- | buildGraph wraps around build tree--build trees and adds network edges after build if network
-- with appropriate options
-- transforms graph type to Tree for builds then back to initial graph type
buildGraph :: [Argument] -> GlobalSettings -> ProcessedData ->  [[VertexCost]] -> Int-> [ReducedPhylogeneticGraph]
buildGraph inArgs inGS inData pairwiseDistances rSeed =
   let fstArgList = fmap (fmap toLower . fst) inArgs
       sndArgList = fmap (fmap toLower . snd) inArgs
       lcArgList = zip fstArgList sndArgList
       checkCommandList = checkCommandArgs "build" fstArgList VER.buildArgList
   in
       -- check for valid command options
   if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'build': " ++ show inArgs)
   else
       let -- block build options including number of display trees to return
           buildBlock = filter ((=="block").fst) lcArgList
           displayBlock = filter ((=="displaytrees").fst) lcArgList
           numDisplayTrees
                | length displayBlock > 1 =
                  errorWithoutStackTrace ("Multiple displayTree number specifications in command--can have only one: " ++ show inArgs)
                | null displayBlock = Just 10
                | null (snd $ head displayBlock) = Just 10
                | otherwise = readMaybe (snd $ head displayBlock) :: Maybe Int


           returnList = filter ((=="return").fst) lcArgList
           numReturnTrees
                | length returnList > 1 =
                  errorWithoutStackTrace ("Multiple 'return' number specifications in command--can have only one: " ++ show inArgs)
                | null returnList = Just (maxBound :: Int)
                | null (snd $ head returnList) = Just (maxBound :: Int)
                | otherwise = readMaybe (snd $ head returnList) :: Maybe Int

           doEUN' = any ((=="eun").fst) lcArgList
           doCUN' = any ((=="cun").fst) lcArgList
           doEUN = if not doEUN' && not doCUN' then True
                   else doEUN'
           returnTrees' = any ((=="displaytrees").fst) lcArgList
           returnGraph' = any ((=="graph").fst) lcArgList
           returnRandomDisplayTrees' = any ((=="atrandom").fst) lcArgList
           returnFirst' = any ((=="first").fst) lcArgList
           buildDistance = any ((=="distance").fst) lcArgList

           -- temprary change (if needed) to buyild tree structures
           inputGraphType = graphType inGS
           treeGS = inGS {graphType = Tree}

           -- really only trees now--but maybe later if can ensure phylogenetic graph from recocnile
           (returnGraph, returnTrees)  = if (graphType inGS) == Tree then (False, True)
                                         else
                                           if returnGraph' || returnTrees' then (returnGraph', returnTrees')
                                           else (False, True)

           -- default to return reandom and overrides if both specified
           (returnRandomDisplayTrees, _) = if returnRandomDisplayTrees' || returnFirst' then (returnRandomDisplayTrees', returnFirst')
                                                     else (True, False)

           -- initial build of trees from combined data--or by blocks
           firstGraphs' = if null buildBlock then
                            let simpleTreeOnly = False
                            in
                            buildTree simpleTreeOnly inArgs treeGS inData pairwiseDistances rSeed
                         else -- removing taxa with missing data for block
                            trace ("Block building initial graph(s)") (
                            let simpleTreeOnly = True
                                processedDataList = U.getProcessDataByBlock True inData
                                distanceMatrixList = if buildDistance then PU.seqParMap PU.myStrategyHighLevel DD.getPairwiseDistances processedDataList -- `using` PU.myParListChunkRDS
                                                     else replicate (length processedDataList) []

                                blockTrees = concat (PU.seqParMap PU.myStrategyHighLevel (buildTree' simpleTreeOnly inArgs treeGS rSeed) (zip distanceMatrixList processedDataList)) --  `using` PU.myParListChunkRDS)
                                -- blockTrees = concat (PU.myChunkParMapRDS (buildTree' simpleTreeOnly inArgs treeGS inputGraphType seed) (zip distanceMatrixList processedDataList))

                                -- reconcile trees and return graph and/or display trees (limited by numDisplayTrees) already re-optimized with full data set
                                returnGraphs = reconcileBlockTrees rSeed blockTrees (fromJust numDisplayTrees) returnTrees returnGraph returnRandomDisplayTrees doEUN
                            in
                            -- trace (concatMap LG.prettify returnGraphs)
                            -- trace ("BG: " ++ (concatMap LG.prettyIndices returnGraphs))
                            PU.seqParMap PU.myStrategyHighLevel (T.multiTraverseFullyLabelGraphReduced inGS inData True True Nothing) returnGraphs -- `using` PU.myParListChunkRDS
                            )

           -- this to allow 'best' to return more trees then later 'returned' and contains memory by letting other graphs go out of scope
           firstGraphs = if null buildBlock then
                            GO.selectGraphs Unique (fromJust numReturnTrees) 0.0 (-1) firstGraphs'
                         else firstGraphs'

           -- reporting info
           returnString = if (not . null) firstGraphs then
                            ("\tReturning " ++ (show $ length firstGraphs) ++ " graphs at cost range " ++ (show (minimum $ fmap snd5 firstGraphs, maximum $ fmap snd5 firstGraphs)))
                          else "\t\tReturning 0 graphs"

           costString = if (not . null) firstGraphs then
                            ("\tBlock build yielded " ++ (show $ length firstGraphs) ++ " graphs at cost range " ++ (show (minimum $ fmap snd5 firstGraphs, maximum $ fmap snd5 firstGraphs)))
                        else "\t\tBlock build returned 0 graphs"

       in
       if isNothing numDisplayTrees then
                                errorWithoutStackTrace ("DisplayTrees specification in build not an integer: "  ++ show (snd $ head displayBlock))
       else if isNothing numReturnTrees then
                                errorWithoutStackTrace ("Return number specifications in build not an integer: " ++ show (snd $ head returnList))

       else
            trace returnString (
            if inputGraphType == Tree || (not . null) buildBlock then
              -- trace ("BB: " ++ (concat $ fmap  LG.prettify $ fmap fst6 firstGraphs)) (
              if null buildBlock then firstGraphs
              else trace (costString) firstGraphs
              -- )
            else
              trace ("\tRediagnosing as " ++ (show (graphType inGS)))
              PU.seqParMap PU.myStrategyHighLevel (T.multiTraverseFullyLabelGraphReduced inGS inData False False Nothing) (fmap fst5 firstGraphs) -- `using` PU.myParListChunkRDS
            )


-- | reconcileBlockTrees takes a lists of trees (with potentially varying leave complement) and reconciled them
-- as per the arguments producing a set of displayTrees (ordered or resolved random), and/or the reconciled graph
-- all outputs are re-optimzed and ready to go
reconcileBlockTrees ::  Int -> [ReducedPhylogeneticGraph] -> Int -> Bool -> Bool -> Bool ->  Bool -> [SimpleGraph]
reconcileBlockTrees rSeed blockTrees numDisplayTrees returnTrees returnGraph returnRandomDisplayTrees doEUN =
      -- trace ("Reconcile producing " ++ (show numDisplayTrees)) (
      let -- numLeaves = V.length $ fst3 inData
          -- fullLeafSet = zip [0..(numLeaves - 1)] (V.toList $ fst3 inData)
          simpleGraphList = fmap fst5 blockTrees
          -- fullLeafGraphList = fmap (E.makeProcessedGraph fullLeafSet) simpleGraphList
          reconcileArgList = if doEUN then [("eun", []), ("vertexLabel:true", []), ("connect:True", [])]
                             else [("cun", []), ("vertexLabel:true", []), ("connect:True", [])]

          -- create reconciled graph--NB may NOT be phylogenetic graph--time violations etc.
          reconciledGraphInitial = snd $ R.makeReconcileGraph VER.reconcileArgList reconcileArgList simpleGraphList

          -- ladderize, time consistent-ized, removed chained network edges, removed treenodes with all network edge children
          {-
          reconciledGraph' = GO.convertGeneralGraphToPhylogeneticGraph "correct" reconciledGraphInitial
          noChainedGraph = LG.removeChainedNetworkNodes False reconciledGraph'
          noTreeNodesWithAllNetChildren = LG.removeTreeEdgeFromTreeNodeWithAllNetworkChildren $ fromJust noChainedGraph
          contractedGraph = GO.contractIn1Out1EdgesRename noTreeNodesWithAllNetChildren
          reconciledGraph = GO.convertGeneralGraphToPhylogeneticGraph "correct" contractedGraph -}

          -- chained was separate in past, now in convertGeneralGraphToPhylogeneticGraph
          reconciledGraph = GO.convertGeneralGraphToPhylogeneticGraph True reconciledGraphInitial


          -- this for non-convertable graphs
          reconciledGraph' = if LG.isEmpty reconciledGraph then reconciledGraphInitial
                             else if reconciledGraph /= LG.empty then reconciledGraph
                             else reconciledGraphInitial


          contractIn1Out1Nodes = True
          displayGraphs' = if not returnRandomDisplayTrees then take numDisplayTrees $ LG.generateDisplayTrees contractIn1Out1Nodes reconciledGraph'
                           else LG.generateDisplayTreesRandom rSeed numDisplayTrees reconciledGraph'

          -- need this to fix up some graphs after other stuff chnaged
          displayGraphs = fmap (GO.convertGeneralGraphToPhylogeneticGraph True) displayGraphs'

          -- displayGraphs = fmap GO.ladderizeGraph $ fmap GO.renameSimpleGraphNodes displayGraphs'

          numNetNodes = length $ fth4 (LG.splitVertexList reconciledGraph)
      in
      if LG.isEmpty reconciledGraph && not returnTrees then
        errorWithoutStackTrace "\n\n\tError--reconciled graph could not be converted to phylogenetic graph.  Consider modifying block tree search options or returning display trees."
      -- trace ("RBT: " ++ (LG.prettyIndices reconciledGraph') ++ "\n" ++ (LG.prettyIndices reconciledGraph)) (
      else if (not $ LG.isEmpty reconciledGraph) && not returnTrees then
        trace ("Reconciled graph has " ++ show numNetNodes ++ " network nodes hence up to 2^" ++ show numNetNodes ++ " display trees for softwired network")
        [reconciledGraph]
      else if not returnGraph && returnTrees then
         displayGraphs
      else
         if (not $ LG.isEmpty reconciledGraph) then
            trace ("\n\tReconciled graph has " ++ show numNetNodes ++ " network nodes hence up to 2^" ++ show numNetNodes ++ " display trees")
                -- ++ "\n" ++ (LG.prettyIndices reconciledGraph') ++ "\n" ++ (LG.prettyIndices reconciledGraph))
            reconciledGraph : displayGraphs
         else
            trace "\n\tWarning--reconciled graph could not be converted to phylogenetic graph.  Consider modifying block tree search options or performing standard builds."
            displayGraphs
     -- ))


-- | buildTree' wraps build tree and changes order of arguments for mapping
buildTree' :: Bool-> [Argument] -> GlobalSettings -> Int -> ([[VertexCost]], ProcessedData) -> [ReducedPhylogeneticGraph]
buildTree' simpleTreeOnly inArgs inGS rSeed (pairwiseDistances, inData) =
   buildTree simpleTreeOnly inArgs inGS inData pairwiseDistances rSeed

-- | buildTree takes build options and returns contructed graphList
-- simpleTreeOnly (for block build) returns a single best tree to reduce edges in
-- reconcile step
buildTree :: Bool -> [Argument] -> GlobalSettings ->ProcessedData ->  [[VertexCost]] -> Int-> [ReducedPhylogeneticGraph]
buildTree simpleTreeOnly inArgs inGS inData@(nameTextVect, _, _) pairwiseDistances rSeed =
   let fstArgList = fmap (fmap toLower . fst) inArgs
       sndArgList = fmap (fmap toLower . snd) inArgs
       lcArgList = zip fstArgList sndArgList
       checkCommandList = checkCommandArgs "build" fstArgList VER.buildArgList
   in
   -- check for valid command options
   if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'build': " ++ show inArgs)
   else

      let buildDistance = any ((=="distance").fst) lcArgList
          buildCharacter = any ((=="character").fst) lcArgList
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
      if buildDistance && buildCharacter then
         errorWithoutStackTrace ("Cannot specify both 'character' and 'distance' builds in same build command" ++ show inArgs)
      else if isNothing numReplicates then errorWithoutStackTrace ("Replicates specification not an integer: "  ++ show (snd $ head repPairList))
      else if isNothing numToSave then errorWithoutStackTrace ("Best specification not an integer: "  ++ show (snd $ head keepPairList))
      else if buildDistance then
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
         trace "\tBuilding Distance Wagner" (
         let refinement
               | doTBR = "tbr"
               | doSPR = "spr"
               | doOTU = "otu"
               | otherwise = "none"
             treeList    = [distanceWagner simpleTreeOnly inGS inData nameStringVect distMatrix outgroupElem refinement | doDWag]
             treeList'   = if doRDwag then treeList ++ randomizedDistanceWagner simpleTreeOnly inGS inData nameStringVect distMatrix outgroupElem (fromJust numReplicates) rSeed (fromJust numToSave) refinement
                           else treeList
             treeList''  = if doNJ then treeList' ++ [neighborJoin simpleTreeOnly inGS inData nameStringVect distMatrix outgroupElem refinement]
                           else treeList'
             treeList''' = if doWPGMA then treeList'' ++ [wPGMA simpleTreeOnly inGS inData nameStringVect distMatrix outgroupElem refinement]
                           else treeList''
         in
         if null treeList''' then errorWithoutStackTrace ("Distance build is specified, but without any method: " ++ show inArgs)
         else
            let costRangeString = if not simpleTreeOnly then " at cost range " ++ show (minimum $ fmap snd5 treeList''', maximum $ fmap snd5 treeList''')
                                  else ""
            in
            trace ("\tDistance build yielded " ++ show (length treeList''') ++ " trees" ++ costRangeString) (
            if not simpleTreeOnly then treeList'''
            else GO.selectGraphs Best 1 0.0 (-1) treeList'''
            )
         )

      else
         -- character build
         -- final diagnosis in input graph type
         trace "\tBuilding Character Wagner" (
         let treeList = WB.rasWagnerBuild inGS inData rSeed (fromJust numReplicates)
             costRangeString = if not simpleTreeOnly then " at cost range " ++ show (minimum $ fmap snd5 treeList, maximum $ fmap snd5 treeList)
                               else ""
         in
         if not simpleTreeOnly then
            trace ("\tCharacter build yielded " ++ show (length treeList) ++ " tree(s)" ++ costRangeString)
            treeList
         else
            trace ("\tCharacter build yielded " ++ show (length $ GO.selectGraphs Best 1 0.0 (-1) treeList) ++ " tree(s)" ++ costRangeString)
            GO.selectGraphs Best 1 0.0 (-1) treeList
         )

-- | distanceWagner takes Processed data and pairwise distance matrix and returns
-- 'best' addition sequence Wagner (defined in Farris, 1972) as fully decorated tree (as Graph)
distanceWagner :: Bool -> GlobalSettings -> ProcessedData -> V.Vector String -> M.Matrix Double -> Int -> String -> ReducedPhylogeneticGraph
distanceWagner simpleTreeOnly inGS inData leafNames distMatrix outgroupValue refinement =
   let distWagTree = head $ DM.doWagnerS leafNames distMatrix "closest" outgroupValue "best" []
       distWagTree' = head $ DW.performRefinement refinement "best:1" "first" leafNames outgroupValue distWagTree
       distWagTreeSimpleGraph = DU.convertToDirectedGraphText leafNames outgroupValue (snd4 distWagTree')
       charInfoVV = V.map thd3 $ thd3 inData
   in
   if not simpleTreeOnly then T.multiTraverseFullyLabelGraphReduced inGS inData False False Nothing (GO.renameSimpleGraphNodes $ GO.dichotomizeRoot outgroupValue $ LG.switchRootTree (length leafNames) distWagTreeSimpleGraph)
   else
      let simpleWag = GO.renameSimpleGraphNodes $ GO.dichotomizeRoot outgroupValue $ LG.switchRootTree (length leafNames) distWagTreeSimpleGraph
      in
      (simpleWag, 0.0, LG.empty, V.empty, charInfoVV)

-- | randomizedDistanceWagner takes Processed data and pairwise distance matrix and returns
-- random addition sequence Wagner trees fully decorated tree (as Graph)
randomizedDistanceWagner :: Bool -> GlobalSettings -> ProcessedData ->  V.Vector String -> M.Matrix Double -> Int -> Int -> Int -> Int -> String -> [ReducedPhylogeneticGraph]
randomizedDistanceWagner simpleTreeOnly inGS inData leafNames distMatrix outgroupValue numReplicates rSeed numToKeep refinement =
   let randomizedAdditionSequences = V.fromList <$> shuffleInt rSeed numReplicates [0..(length leafNames - 1)]
       randomizedAdditionWagnerTreeList = DM.doWagnerS leafNames distMatrix "random" outgroupValue "random" randomizedAdditionSequences
       randomizedAdditionWagnerTreeList' = take numToKeep $ L.sortOn thd4 randomizedAdditionWagnerTreeList
       randomizedAdditionWagnerTreeList'' = head <$> PU.seqParMap PU.myStrategyHighLevel (DW.performRefinement refinement "best:1"  "first" leafNames outgroupValue) randomizedAdditionWagnerTreeList'
       randomizedAdditionWagnerSimpleGraphList = fmap (DU.convertToDirectedGraphText leafNames outgroupValue . snd4) randomizedAdditionWagnerTreeList''
       charInfoVV = V.map thd3 $ thd3 inData
   in
   if not simpleTreeOnly then fmap ((T.multiTraverseFullyLabelGraphReduced inGS inData False False Nothing . GO.renameSimpleGraphNodes . GO.dichotomizeRoot outgroupValue) . LG.switchRootTree (length leafNames)) randomizedAdditionWagnerSimpleGraphList
   else
      let numTrees = length randomizedAdditionWagnerSimpleGraphList
          simpleRDWagList = fmap (GO.dichotomizeRoot outgroupValue . LG.switchRootTree (length leafNames)) randomizedAdditionWagnerSimpleGraphList
      in
      L.zip5 simpleRDWagList (replicate numTrees 0.0) (replicate numTrees LG.empty) (replicate numTrees V.empty) (replicate numTrees charInfoVV)

-- | neighborJoin takes Processed data and pairwise distance matrix and returns
-- Neighbor-Joining tree as fully decorated tree (as Graph)
neighborJoin :: Bool -> GlobalSettings -> ProcessedData -> V.Vector String -> M.Matrix Double -> Int -> String -> ReducedPhylogeneticGraph
neighborJoin simpleTreeOnly inGS inData leafNames distMatrix outgroupValue refinement =
   let njTree = DM.neighborJoining leafNames distMatrix outgroupValue
       njTree' = head $ DW.performRefinement refinement "best:1" "first" leafNames outgroupValue njTree
       njSimpleGraph = DU.convertToDirectedGraphText leafNames outgroupValue (snd4 njTree')
       charInfoVV = V.map thd3 $ thd3 inData
   in
   if not simpleTreeOnly then T.multiTraverseFullyLabelGraphReduced inGS inData False False Nothing (GO.renameSimpleGraphNodes $ GO.dichotomizeRoot outgroupValue $ LG.switchRootTree (length leafNames) njSimpleGraph)
   else
      let simpleNJ = GO.dichotomizeRoot outgroupValue $ LG.switchRootTree (length leafNames) njSimpleGraph
      in
      (simpleNJ, 0.0, LG.empty, V.empty, charInfoVV)


-- | wPGMA takes Processed data and pairwise distance matrix and returns
-- WPGMA tree as fully decorated tree (as Graph)
-- since root index not nOTUs as with other tres--chanegd as with dWag and NJ to make consistent.
wPGMA :: Bool -> GlobalSettings -> ProcessedData -> V.Vector String -> M.Matrix Double -> Int -> String -> ReducedPhylogeneticGraph
wPGMA simpleTreeOnly inGS inData leafNames distMatrix outgroupValue refinement =
   let wpgmaTree = DM.wPGMA leafNames distMatrix outgroupValue
       wpgmaTree' = head $ DW.performRefinement refinement "best:1" "first" leafNames outgroupValue wpgmaTree
       wpgmaSimpleGraph = DU.convertToDirectedGraphText leafNames outgroupValue (snd4 wpgmaTree')
       charInfoVV = V.map thd3 $ thd3 inData
   in
   if not simpleTreeOnly then T.multiTraverseFullyLabelGraphReduced inGS inData False False Nothing (GO.renameSimpleGraphNodes $ GO.dichotomizeRoot outgroupValue $ LG.switchRootTree (length leafNames) wpgmaSimpleGraph)
   else
      let simpleWPGMA = GO.dichotomizeRoot outgroupValue $ LG.switchRootTree (length leafNames) wpgmaSimpleGraph
      in
      (simpleWPGMA, 0.0, LG.empty, V.empty, charInfoVV)
