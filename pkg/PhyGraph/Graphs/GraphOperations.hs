{- |
Module      :  GraphOperations.hs
Description :  Module specifying general graph functions--with types specific to Types.hs
               graph functions that a re general are in LocalGraph.hs
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

module Graphs.GraphOperations (  ladderizeGraph
                               , convertDecoratedToSimpleGraph
                               , convertToSimpleEdge
                               , graphCostFromNodes
                               , dichotomizeRoot
                               , showDecGraphs
                               , sortEdgeListByLength
                               , selectPhylogeneticGraph
                               , getUniqueGraphs
                               , copyIAFinalToPrelim
                               , copyIAPrelimToFinal
                               , makeIAFinalFromPrelim
                               , makeIAPrelimFromFinal
                               , topologicalEqual
                               , getTopoUniqPhylogeneticGraph
                               , getBVUniqPhylogeneticGraph
                               , makeDummyLabEdge
                               , contractIn1Out1EdgesRename
                               , renameSimpleGraphNodes
                               , renameSimpleGraphNodesString
                               , hasNetNodeAncestorViolation
                               , convertGeneralGraphToPhylogeneticGraph
                               , selectGraphStochastic
                               , makeNewickList
                               , makeGraphTimeConsistent
                               , isNovelGraph
                               , getNodeType
                               , getDisplayTreeCostList
                               ) where

import           Bio.DynamicCharacter
import           Control.Parallel.Strategies
import           Data.Bits
import qualified Data.BitVector.LittleEndian as BV
import qualified Data.Char                   as C
import qualified Data.List                   as L
import           Data.Maybe
import qualified Data.Text.Lazy              as T
import qualified Data.Vector                 as V
import qualified Data.Vector.Generic         as GV
import           Debug.Trace
import           GeneralUtilities
import qualified GraphFormatUtilities        as GFU
import qualified GraphOptimization.Medians   as M
import qualified ParallelUtilities           as PU
import           Text.Read
import           Types.Types
import qualified Utilities.LocalGraph        as LG
import qualified Utilities.Utilities         as U

-- | makeNewickList takes a list of fgl trees and outputs a single String cointaining the graphs in Newick format
makeNewickList ::  Bool -> Bool -> Int -> [SimpleGraph] -> [VertexCost] -> String
makeNewickList writeEdgeWeight writeNodeLabel' rootIndex graphList costList =
    let allTrees = L.foldl' (&&) True (fmap LG.isTree graphList)

        -- check for network HTU label requirement
        writeNodeLabel = if allTrees then writeNodeLabel'
                         else if writeNodeLabel' then writeNodeLabel'
                         else 
                            trace ("HTU labels are required for ENewick Output")
                            True

        graphString = GFU.fglList2ForestEnhancedNewickString (fmap (LG.rerootTree rootIndex) graphList)  writeEdgeWeight writeNodeLabel
        newickStringList = fmap init $ filter (not . null) $ lines graphString
        costStringList  = fmap (('[' :) . (++ "];\n")) (fmap show costList)
        graphStringCost = concat $ zipWith (++) newickStringList costStringList
    in
    graphStringCost

-- | convertGeneralGraphToPhylogeneticGraph inputs a SimpleGraph and converts it to a Phylogenetic graph by:
--  1) transitive reduction -- removes anc <-> desc netork edges
--  2) ladderizes -- all vertices are (in degree, outdegree) (0,1|2) (1,2) (2,1) (1,0)
--        by adding extra HTIs and edges
--  3) checks time consistency and removes edges stepwise from
--        those that violate the most  before/after splits of network edges
--        arbitrary but deterministic
--  4) contracts out any remaning indegree 1 outdegree 1 nodes and renames HTUs in order
-- these tests can be screwed up by imporperly formated graphs comming in (self edges, chained network edge etc)
convertGeneralGraphToPhylogeneticGraph :: String -> SimpleGraph -> SimpleGraph
convertGeneralGraphToPhylogeneticGraph failCorrect inGraph =
  if LG.isEmpty inGraph then LG.empty
  else
    let -- remove single "tail" edge from root with single child, replace child node with root
        noTailGraph = LG.contractRootOut1Edge inGraph

        -- remove indeg 1 out deg 1 edges
        noIn1Out1Graph = contractIn1Out1EdgesRename noTailGraph

        -- transitive reduction
        -- only wanted to EUN and CUN--but they do it
        -- reducedGraph = LG.transitiveReduceGraph noIn1Out1Graph

        -- laderization of indegree and outdegree edges
        ladderGraph = ladderizeGraph noIn1Out1Graph -- reducedGraph

        -- time consistency (after those removed by transitrive reduction)
        timeConsistentGraph = makeGraphTimeConsistent failCorrect ladderGraph

        -- removes ancestor descendent edges transitiveReduceGraph should do this
        -- but that looks at all nodes not just vertex
        noParentChainGraph = removeParentsInChain failCorrect timeConsistentGraph

        -- remove sister-sister edge.  where two network nodes have same parents
        noSisterSisterGraph = removeSisterSisterEdges failCorrect noParentChainGraph

    in
    if LG.isEmpty timeConsistentGraph then LG.empty

    else if LG.isEmpty noParentChainGraph then LG.empty

    else if LG.isEmpty noSisterSisterGraph then LG.empty

    -- trace ("CGP orig:\n" ++ (LG.prettify inGraph) ++ "\nNew:" ++ (LG.prettify timeConsistentGraph))
    -- cycle check to make sure--can be removed when things working
    -- else if LG.cyclic noSisterSisterGraph then error ("Cycle in graph : \n" ++ (LG.prettify noSisterSisterGraph))

    -- this final need to ladderize or recontract?
    else noSisterSisterGraph
    
-- | removeParentsInChain checks the parents of each netowrk node are not anc/desc of each other
removeParentsInChain :: String -> SimpleGraph -> SimpleGraph
removeParentsInChain failCorrect inGraph =
  if LG.isEmpty inGraph then LG.empty
  else
      let (_, _, _, netVertexList) = LG.splitVertexList inGraph
          parentNetVertList = fmap (LG.labParents inGraph) $ fmap fst netVertexList

          -- get list of nodes that are transitively equal in age
          concurrentList = LG.mergeConcurrentNodeLists parentNetVertList []
          concurrentPairList = concatMap getListPairs concurrentList

          -- get pairs that violate concurrency
          violatingConcurrentPairs = concatMap (LG.concurrentViolatePair inGraph) concurrentPairList

          -- get netowrk nodes with violations
          parentNodeViolateList = concatMap pairToList violatingConcurrentPairs
          childNodeViolateList = concatMap (LG.descendants inGraph) parentNodeViolateList
          netNodeViolateList = filter (LG.isNetworkNode inGraph) childNodeViolateList

          netEdgesThatViolate = fmap LG.toEdge $ LG.inn inGraph $ head netNodeViolateList

      in
      if null violatingConcurrentPairs then inGraph
      else if null netNodeViolateList then error ("Should be neNode that violate")
      else if null netEdgesThatViolate then error "Should be violating in edges"
      else if failCorrect == "fail" then LG.empty
      else
        let edgeDeletedGraph = LG.delEdge (head netEdgesThatViolate) inGraph
            newGraph = contractIn1Out1EdgesRename edgeDeletedGraph
        in
        -- trace ("PIC")
        removeParentsInChain failCorrect newGraph
    where pairToList (a,b) = [fst a, fst b]

-- | removeSisterSisterEdges takes a graph and recursively removes a single edge fomr where two network
-- edges have the same two parents
removeSisterSisterEdges :: String -> SimpleGraph -> SimpleGraph
removeSisterSisterEdges failCorrect inGraph =
  if LG.isEmpty inGraph then LG.empty
  else
    let sisterSisterEdges = LG.getSisterSisterEdgeList inGraph
        -- newGraph = LG.delEdge (head sisterSisterEdges) inGraph
        newGraph = LG.delEdges sisterSisterEdges inGraph
        newGraph' = contractIn1Out1EdgesRename newGraph
    in
    if null sisterSisterEdges then inGraph
    else if failCorrect == "fail" then LG.empty
    else 
      -- trace ("Sister")
      -- removeSisterSisterEdges  
      newGraph'

-- | makeGraphTimeConsistent takes laderized, transitive reduced graph and deletes
-- network edges in an arbitrary but deterministic sequence to produce a phylogentic graphs suitable
-- for swapping etc
-- looks for violation of time between netork edges based on "before" and "after"
-- tests of nodes that should be potentially same age
-- removes second edge of second pair of two network edges in each case adn remakes graph
makeGraphTimeConsistent :: String -> SimpleGraph -> SimpleGraph
makeGraphTimeConsistent failOut inGraph =
  if LG.isEmpty inGraph then LG.empty
  else if LG.isTree inGraph then inGraph
  else
    let coevalNodeConstraintList = LG.coevalNodePairs inGraph
        coevalNodeConstraintList' = PU.seqParMap rdeepseq  (LG.addBeforeAfterToPair inGraph) coevalNodeConstraintList -- `using`  PU.myParListChunkRDS
        coevalPairsToCompareList = getListPairs coevalNodeConstraintList'
        timeOffendingEdgeList = LG.getEdgesToRemoveForTime inGraph coevalPairsToCompareList
        newGraph = LG.delEdges timeOffendingEdgeList inGraph
    in
    -- trace ("MGTC:" ++ (show timeOffendingEdgeList))
    if (failOut == "fail") && ((not . null) timeOffendingEdgeList) then LG.empty
    else contractIn1Out1EdgesRename newGraph

-- | contractIn1Out1EdgesRename contracts in degree and outdegree edges and renames HTUs in index order
-- does one at a time and makes a graph and recurses
contractIn1Out1EdgesRename :: SimpleGraph -> SimpleGraph
contractIn1Out1EdgesRename inGraph =
  if LG.isEmpty inGraph then LG.empty
  else
    let newGraph = LG.contractIn1Out1Edges inGraph
    in
    renameSimpleGraphNodes newGraph


-- | renameSimpleGraphNodes takes nodes and renames HTU nodes based on index
renameSimpleGraphNodes :: SimpleGraph -> SimpleGraph
renameSimpleGraphNodes inGraph =
  if LG.isEmpty inGraph then LG.empty
    else
      let inNodes = LG.labNodes inGraph
          nodeLabels = fmap (makeSimpleLabel inGraph) inNodes
          newNodes = zip (fmap fst inNodes) nodeLabels
          newEdges = LG.labEdges inGraph
    in
    --newGraph
    -- trace ("C11: " ++ (show $ LG.getIsolatedNodes newGraph) ++ " => " ++ (show newNodes) ++ " " ++ (show $ fmap LG.toEdge newEdges))
    LG.mkGraph newNodes newEdges
    where makeSimpleLabel g (a, b)  = if (not $ LG.isLeaf g a) then T.pack $ "HTU"  ++ show a
                                      else b

-- | renameSimpleGraphNodesString takes nodes and renames HTU nodes based on index
renameSimpleGraphNodesString :: LG.Gr String String -> LG.Gr String String
renameSimpleGraphNodesString inGraph =
  if LG.isEmpty inGraph then LG.empty
    else
      let inNodes = LG.labNodes inGraph
          nodeLabels = fmap (makeSimpleLabel inGraph) inNodes
          newNodes = zip (fmap fst inNodes) nodeLabels
          newEdges = LG.labEdges inGraph
    in
    --newGraph
    -- trace ("C11: " ++ (show $ LG.getIsolatedNodes newGraph) ++ " => " ++ (show newNodes) ++ " " ++ (show $ fmap LG.toEdge newEdges))
    LG.mkGraph newNodes newEdges
    where makeSimpleLabel g (a, b)  = if (not $ LG.isLeaf g a) then "HTU"  ++ show a
                                      else b
                                      
-- | sortEdgeListByLength sorts edge list by length (midRange), highest to lowest
sortEdgeListByLength :: [LG.LEdge EdgeInfo] -> [LG.LEdge EdgeInfo]
sortEdgeListByLength inEdgeList =
  if null inEdgeList then []
  else
    reverse $ L.sortOn (midRangeLength . thd3) inEdgeList

-- | ladderizeGraph is a wrapper around ladderizeGraph' to allow for mapping with
-- local nodelist
ladderizeGraph :: SimpleGraph -> SimpleGraph
ladderizeGraph inGraph = ladderizeGraph' inGraph (LG.nodes inGraph)

-- | ladderize takes an input graph and ensures/creates nodes
-- such that all vertices are (indegree, outdegree) (0,>0), (1,2) (2,1) (1,0)
--ladderizeGraph' :: SimpleGraph -> [LG.Node] -> SimpleGraph
ladderizeGraph' :: SimpleGraph -> [LG.Node] -> SimpleGraph
ladderizeGraph' inGraph nodeList
  | LG.isEmpty inGraph = LG.empty
  | null nodeList = inGraph
  | otherwise =
  let -- these are roots, network, tree, leaf nodes
      okNodeDegrees = [(0,2),(1,2),(2,1),(1,0)]
      firstNode = head nodeList
      (inEdgeList, outEdgeList) = LG.getInOutEdges inGraph firstNode
      inOutPairLength = (length inEdgeList, length outEdgeList)
  in
  -- trace ("node " ++ (show firstNode) ++ " " ++ (show inOutPairLength)) (
  -- node ok to keep
  if inOutPairLength `elem` okNodeDegrees then ladderizeGraph' inGraph (tail nodeList)
  -- node edges need modification
  else
    let newGraph = resolveNode inGraph firstNode (inEdgeList, outEdgeList) inOutPairLength
    in
    -- trace ("resolving " ++ "node " ++ (show firstNode) ++ " " ++ (show inOutPairLength) )
    ladderizeGraph' newGraph (LG.nodes newGraph)
    -- )

-- | resolveNode takes a graph and node and inbound edgelist and outbound edge list
-- and converts node to one of (indeg, outdeg) (0,1),(0,2),(1,2),(2,1),(1,0)
-- this only resolves a single nodes edges at a time and then returns new graph
-- when more hase to be done--that will occur on lultiple passes through nodes.
-- perhaps not the most efficient, but only done once per input graph
-- contracts indegree 1 outdegree 1 nodes
resolveNode :: SimpleGraph -> LG.Node -> ([LG.LEdge Double], [LG.LEdge Double]) -> (Int, Int) -> SimpleGraph
resolveNode inGraph curNode inOutPair@(inEdgeList, outEdgeList) (inNum, outNum) =
  if LG.isEmpty inGraph then LG.empty
  else
    -- trace ("Resolving " ++ show (inNum, outNum)) (
    let numNodes = length $ LG.nodes inGraph
    in
    -- isolated node -- throw warning and delete
    if inNum == 0 && outNum == 0 then 
      trace ("Warning: ResolveNode deleting isolated vertex " ++ show curNode) ( --  ++ " in graph\n" ++ LG.prettify inGraph )
      let newGraph = LG.delNode curNode inGraph
      in
      newGraph
      )

     -- node with too many parents and too many children
    -- converts to tree node--biased in that direction
    else if (inNum > 2) && (outNum > 2) then
      let first2Edges = take 2 inEdgeList
          newNode = (numNodes , T.pack $ ("HTU" ++ (show numNodes)))
          newEdge1 = (fst3 $ head first2Edges, numNodes, 0.0 :: Double)
          newEdge2 = (fst3 $ last first2Edges, numNodes, 0.0 :: Double)
          newEdge3 = (numNodes, curNode, 0.0 :: Double)
          newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
      in
      newGraph

    -- leaf leaf with too many parents
    else if (inNum > 1) && (outNum == 0) || (inNum > 2) && (outNum == 1) || (inNum > 1) && (outNum == 2) then
      let first2Edges = take 2 inEdgeList
          newNode = (numNodes , T.pack $ ("HTU" ++ (show numNodes)))
          newEdge1 = (fst3 $ head first2Edges, numNodes, 0.0 :: Double)
          newEdge2 = (fst3 $ last first2Edges, numNodes, 0.0 :: Double)
          newEdge3 = (numNodes, curNode, 0.0 :: Double)
          newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
      in
      newGraph

    -- indegree 1 outdegree 1 node to contract
    else if inNum == 1 && outNum == 1 then
      let newEdge = (fst3 $ head inEdgeList, snd3 $ head outEdgeList, 0.0 :: Double)
          newGraph = LG.insEdge newEdge $ LG.delNode curNode $ LG.delLEdges (inEdgeList ++ outEdgeList) inGraph
      in
      newGraph

    else if (inNum < 2 || outNum > 2) then
      let  first2Edges = take 2 outEdgeList
           newNode = (numNodes , T.pack $ ("HTU" ++ (show numNodes)))
           newEdge1 = (numNodes, snd3 $ head first2Edges, 0.0 :: Double)
           newEdge2 = (numNodes, snd3 $ last first2Edges, 0.0 :: Double)
           newEdge3 = (curNode, numNodes, 0.0 :: Double)
           newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
      in
      newGraph

        -- root or simple network indegree node
    else if (inNum == 0 || outNum > 2) then
      let first2Edges = take 2 outEdgeList
          newNode = (numNodes , T.pack $ ("HTU" ++ (show numNodes)))
          newEdge1 = (numNodes, snd3 $ head first2Edges, 0.0 :: Double)
          newEdge2 = (numNodes, snd3 $ last first2Edges, 0.0 :: Double)
          newEdge3 = (curNode, numNodes, 0.0 :: Double)
          newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
      in
      newGraph

    -- check if indegree 0 is a leaf (ie index < root)
    else if outNum == 0 then
      -- get root index 
      let rootIndex = fst $ head $ LG.getRoots inGraph
      in
      if curNode < rootIndex then inGraph
      else LG.delNode curNode inGraph


    else error ("This can't happen in resolveNode in/out edge lists don't need to be resolved " ++ show inOutPair ++ "\n" ++ LG.prettify inGraph)
    -- )

-- | convertDecoratedToSimpleGraph
convertDecoratedToSimpleGraph :: DecoratedGraph -> SimpleGraph
convertDecoratedToSimpleGraph inDec =
  if LG.isEmpty inDec then LG.empty
  else
    let decNodeList = LG.labNodes inDec
        newNodeLabels = fmap vertName $ fmap snd decNodeList
        simpleNodes = zip (fmap fst decNodeList) newNodeLabels

        {-
        decEdgeList = LG.labEdges inDec
        sourceList = fmap fst3 decEdgeList
        sinkList = fmap snd3 decEdgeList
        newEdgeLables = replicate (length sourceList) 0.0  -- fmap midRangeLength $ fmap thd3 decEdgeList
        simpleEdgeList = zip3 sourceList sinkList newEdgeLables
        -}
        simpleEdgeList = fmap convertToSimpleEdge $ LG.labEdges inDec
    in
    LG.mkGraph simpleNodes simpleEdgeList


-- | convertToSimpleEdge takes a lables edge and relabels with 0.0
convertToSimpleEdge :: LG.LEdge EdgeInfo -> LG.LEdge Double
convertToSimpleEdge (a, b, c) = (a, b, minLength c)

-- | graphCostFromNodes takes a Decorated graph and returns its cost by summing up the local costs
--  of its nodes
graphCostFromNodes :: DecoratedGraph -> Double
graphCostFromNodes inGraph =
  if LG.isEmpty inGraph then 0.0
  else
    sum $ fmap vertexCost $ fmap snd $ LG.labNodes inGraph

-- | dichotomizeRoot takes greaph and dichotimizes not dichotomous roots in graph
dichotomizeRoot :: Int -> SimpleGraph -> SimpleGraph
dichotomizeRoot lOutgroupIndex inGraph =
  if LG.isEmpty inGraph then LG.empty
  else
    let rootList = LG.getRoots inGraph
        currentRoot = fst $ head rootList
        rootEdgeList = LG.out inGraph $ currentRoot
    in
    -- not a tree error
    if (length rootList /= 1) then error ("Graph input to dichotomizeRoot is not a tree--not single root:" ++ (show rootList))

    -- nothing to do
    else if (length rootEdgeList < 3) then inGraph
    else
      let numVertices = length $ LG.nodes inGraph
          newNode = (numVertices, T.pack $ show numVertices)
          edgesToDelete = filter ((/=lOutgroupIndex) . snd3) rootEdgeList
          newEdgeDestinations = fmap snd3 edgesToDelete
          newEdgeStarts = replicate (length newEdgeDestinations) numVertices
          newEdgeLabels = replicate (length newEdgeDestinations) 0.0
          -- nub for case where root edge in "wrong" direction
          -- doesn't filter edges to delete properly
          newEdgesNewNode = L.nub $ zip3 newEdgeStarts newEdgeDestinations newEdgeLabels
          newRootEdge = (currentRoot, numVertices, 0.0)
      in
      LG.delLEdges edgesToDelete $ LG.insEdges (newRootEdge : newEdgesNewNode) $ LG.insNode newNode inGraph


-- | showBlockGraphs takes a vector of vector of DecoratedGraphs and converte and prettifies outputting a String
showDecGraphs :: V.Vector (V.Vector DecoratedGraph) -> String
showDecGraphs inDecVV =
    if V.null inDecVV then []
    else
        concat $ fmap concat $ V.toList $ fmap V.toList $ fmap (fmap LG.prettify) $ fmap (fmap convertDecoratedToSimpleGraph) inDecVV

-- | selectPhylogeneticGraph takes  a series OF arguments and an input list ot PhylogeneticGraphs
-- and returns or filters that list based on options.
-- uses selectListCostPairs in GeneralUtilities
selectPhylogeneticGraph :: [Argument] -> Int -> [String] -> [PhylogeneticGraph] -> [PhylogeneticGraph]
selectPhylogeneticGraph inArgs rSeed selectArgList curGraphs =
    if null curGraphs then []
    else
        let fstArgList = fmap (fmap C.toLower . fst) inArgs
            sndArgList = fmap (fmap C.toLower . snd) inArgs
            lcArgList = zip fstArgList sndArgList
            checkCommandList = checkCommandArgs "select" fstArgList selectArgList
        in
           -- check for valid command options
           if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'select': " ++ show inArgs)
           else if length inArgs > 1 then errorWithoutStackTrace ("Can only have a single select type per command: "  ++ show inArgs)
           else
                let doBest    = not $ not (any ((=="best").fst) lcArgList)
                    doAll     = not $ not (any ((=="all").fst) lcArgList)
                    doRandom  = not $ not (any ((=="atrandom").fst) lcArgList)
                    doUnique  = not $ not (any ((=="unique").fst) lcArgList)
                    numberToKeep
                      | null lcArgList = Just (maxBound :: Int)
                      | null $ snd $ head lcArgList = Just (maxBound :: Int)
                      | otherwise = readMaybe (snd $ head lcArgList) :: Maybe Int
                in
                if doAll then curGraphs
                else if isNothing numberToKeep then errorWithoutStackTrace ("Number to keep specification not an integer: "  ++ show (snd $ head lcArgList))
                else
                    let -- minimum graph cost
                        minGraphCost = minimum $ fmap snd6 curGraphs

                        -- collapse zero-length branchs for unique
                        curGraphsCollapsed = fmap U.collapseGraph curGraphs

                        -- keep only unique graphs based on non-zero edges--in sorted by cost
                        uniqueGraphList = L.sortOn snd6 $ getUniqueGraphs'' (zip curGraphs curGraphsCollapsed)-- curGraphs --  True curGraphs -- getBVUniqPhylogeneticGraph True curGraphs -- getTopoUniqPhylogeneticGraph True curGraphs
                        
                        -- this to avaoid alot of unncesesary graph comparisons for 'best' graphs
                        bestCostGraphs = filter ((== minGraphCost).snd6) curGraphs
                        uniqueBestGraphs = getUniqueGraphs'' (zip bestCostGraphs (fmap U.collapseGraph bestCostGraphs))

                      in
                    if doUnique then take (fromJust numberToKeep) uniqueGraphList
                    else if doBest then
                     -- trace ("SPG: " ++ (show (minGraphCost, length uniqueGraphList, fmap snd6 uniqueGraphList)))
                      take (fromJust numberToKeep) uniqueBestGraphs
                    else if doRandom then
                         let randList = head $ shuffleInt rSeed 1 [0..(length curGraphs - 1)]
                             (_, shuffledGraphs) = unzip $ L.sortOn fst $ zip randList curGraphs
                         in
                         take (fromJust numberToKeep) shuffledGraphs
                    -- default is all best and unique
                    else
                        uniqueBestGraphs

-- | getUniqueGraphs takes each pair of non-zero edges and conpares them--if equal not added to list
-- maybe chnge to nub LG.pretify graphList?
getUniqueGraphs :: Bool -> [PhylogeneticGraph] -> [PhylogeneticGraph]
getUniqueGraphs removeZeroEdges inGraphList =
  if null inGraphList then []
  else
    let inGraphEdgeList = if removeZeroEdges then fmap (filter ((> 0.0) . minLength . thd3)) $ fmap LG.labEdges $ fmap thd6 inGraphList
                          else fmap LG.labEdges $ fmap thd6 inGraphList
    in
    getUniqueGraphs' (zip inGraphEdgeList inGraphList) []


-- | getUniqueGraphs Using fgl ==
-- basically a nub
-- need to add a collapse function for compare as well
-- takes pairs of (noCollapsed, collapsed) phylogenetic graphs,
-- make strings based on collapsed and returns not collpased
getUniqueGraphs'' :: [(PhylogeneticGraph, PhylogeneticGraph)] -> [PhylogeneticGraph] 
getUniqueGraphs'' inList = nubGraph [] inList

-- | isNovelGraph  checks if a graph is in list of existing graphs
-- uses colapsed representation
isNovelGraph :: [PhylogeneticGraph] -> PhylogeneticGraph -> Bool
isNovelGraph graphList testGraph = 
  if null graphList then True
  else 
      let collapsedInGraph = (LG.prettyIndices . fst6 . U.collapseGraph) testGraph
          collapseGraphList = fmap (LG.prettyIndices . fst6 . U.collapseGraph) graphList
          matchList = filter (== collapsedInGraph) collapseGraphList
      in
      -- trace ("IsNovel: " ++ (show $ null matchList))
      null matchList

-- | keeps and returns unique graphs based on Eq of Topological Simple Graph
-- String prettyIndices w/0 HTU names and branch lengths
-- arbitrarily rooted on 0 for oonsistency
--reversed to keep original order in case sorted on length
nubGraph :: [(PhylogeneticGraph, PhylogeneticGraph, String)] -> [(PhylogeneticGraph, PhylogeneticGraph)] -> [PhylogeneticGraph]
nubGraph curList inList =
  if null inList then reverse $ fmap fst3 curList
  else 
    let (firstGraphNC, firstGraphC) = head inList
        firstString = LG.prettyIndices $ thd6 firstGraphNC
        isMatch = filter (== firstString) (fmap thd3 curList)
    in
    -- trace ("NG: " ++ (show $ null isMatch) ++ " " ++ firstString) (
    if null curList then nubGraph [(firstGraphNC, firstGraphC, firstString)] (tail inList)
    else if null isMatch then nubGraph ((firstGraphNC, firstGraphC, firstString) : curList) (tail inList)
    else nubGraph curList (tail inList)
    -- )

-- | getUniqueGraphs takes each pair of non-zero edges and compares them--if equal not added to list
getUniqueGraphs' :: [([LG.LEdge EdgeInfo], PhylogeneticGraph)] -> [([LG.LEdge EdgeInfo], PhylogeneticGraph)]  -> [PhylogeneticGraph]
getUniqueGraphs' inGraphPairList currentUniquePairs =
    if null inGraphPairList then fmap snd currentUniquePairs
    else
        let firstPair@(firstEdges, _) = head inGraphPairList
        in
        if null currentUniquePairs then getUniqueGraphs' (tail inGraphPairList) [firstPair]
        else
            let equalList = filter (== True) $ fmap ((== firstEdges) . fst) currentUniquePairs
            in
            if null equalList then getUniqueGraphs' (tail inGraphPairList) (firstPair : currentUniquePairs)
            else getUniqueGraphs' (tail inGraphPairList) currentUniquePairs

-- | getNodeType returns node type for Node
getNodeType :: (Show a, Show b) => LG.Gr a b -> LG.Node -> NodeType
getNodeType inGraph inNode =
    if not $ LG.gelem inNode inGraph then error ("Node " ++ (show inNode) ++ " not in graph\n" ++ (GFU.showGraph inGraph))
    else if LG.isLeaf inGraph inNode then LeafNode
    else if LG.isTreeNode inGraph inNode then TreeNode
    else if LG.isNetworkNode inGraph inNode then NetworkNode
    else if LG.isRoot inGraph inNode then RootNode
    else error ("Node type " ++ (show inNode) ++ " not Leaf, Tree, Network, or Root in graph\n" ++ (GFU.showGraph inGraph))

-- | copyIAFinalToPrelim takes a Decorated graph and copies
-- the IA final fields to preliminary IA states--this for IA only optimization
-- inswapping and other operations.  Thi sis done becasue the "preliminary" IA states
-- are only known after full post/pre traversals
copyIAFinalToPrelim :: DecoratedGraph -> DecoratedGraph
copyIAFinalToPrelim inGraph =
  if LG.isEmpty inGraph then error "Empty input graph in copyIAFinalToPrelim"
  else
    let nodes = LG.labNodes inGraph
        edges = LG.labEdges inGraph
        newNodes = fmap makeIAPrelimFromFinal nodes
    in
    LG.mkGraph newNodes edges

-- | makeIAPrelimFromFinal updates the label of a node for IA states
-- setting preliminary to final
makeIAPrelimFromFinal :: LG.LNode VertexInfo -> LG.LNode VertexInfo
makeIAPrelimFromFinal (inIndex, label) =
  let labData = vertData label
      newLabData = fmap (fmap f) labData
  in
  (inIndex, label {vertData = newLabData})
  where f c = if (GV.null $ slimIAFinal c) && (GV.null  $ wideIAFinal c) && (GV.null  $ hugeIAFinal c) then c
              else if (not $ GV.null $ slimIAFinal c) then c {slimIAPrelim = M.makeDynamicCharacterFromSingleVector $ slimIAFinal c}
              else if (not $ GV.null $ wideIAFinal c) then c {wideIAPrelim = M.makeDynamicCharacterFromSingleVector $ wideIAFinal c}
              else c {hugeIAPrelim = M.makeDynamicCharacterFromSingleVector $ hugeIAFinal c}


-- | copyIAPrelimToFinal takes a Decorated graph and copies
-- the IA prelim fields to final IA states--this for IA only optimization
-- inswapping and other operations.  THis is fdone for root and leaf vertices
copyIAPrelimToFinal :: DecoratedGraph -> DecoratedGraph
copyIAPrelimToFinal inGraph =
  if LG.isEmpty inGraph then error "Empty input graph in copyIAFinalToPrelim"
  else
    let nodes = LG.labNodes inGraph
        edges = LG.labEdges inGraph
        newNodes = fmap makeIAFinalFromPrelim nodes
    in
    LG.mkGraph newNodes edges

-- | makeIAFinalFomPrelim updates the label of a node for IA states
-- setting final to preliminary
makeIAFinalFromPrelim:: LG.LNode VertexInfo -> LG.LNode VertexInfo
makeIAFinalFromPrelim (inIndex, label) =
  let labData = vertData label
      newLabData = fmap (fmap f) labData
  in
  (inIndex, label {vertData = newLabData})
  where f c = let newSlimIAFinal = extractMediansGapped $ slimIAPrelim c
                  newWideIAFinal = extractMediansGapped $ wideIAPrelim c
                  newHugeIAFinal = extractMediansGapped $ hugeIAPrelim c
              in
              if (GV.null $ snd3 $ slimIAPrelim c) && (GV.null $ snd3 $ wideIAPrelim c) && (GV.null $ snd3 $ hugeIAPrelim c) then c
              else if (not $ GV.null $ snd3 $ slimIAPrelim c) then c {slimIAFinal = newSlimIAFinal}
              else if (not $ GV.null $ snd3 $ wideIAPrelim c) then c {wideIAFinal = newWideIAFinal}
              else c {hugeIAFinal = newHugeIAFinal}


-- | getTopoUniqPhylogeneticGraph takes a list of phylogenetic graphs and returns
-- list of topologically unique graphs--operatres on simple graph field
-- noZeroEdges flag passed to remove zero weight edges
getTopoUniqPhylogeneticGraph :: Bool -> [PhylogeneticGraph] -> [PhylogeneticGraph]
getTopoUniqPhylogeneticGraph nonZeroEdges inPhyloGraphList =
  if null inPhyloGraphList then []
  else
      let uniqueBoolList = createUniqueBoolList nonZeroEdges (fmap fst6 inPhyloGraphList) []
          boolPair = zip inPhyloGraphList uniqueBoolList
      in
      fmap fst $ filter ((== True) . snd) boolPair

-- | createUniqueBoolList creates a list of Bool if graphs are unique--first occurrence is True, others False
createUniqueBoolList :: Bool -> [SimpleGraph] -> [(SimpleGraph,Bool)] -> [Bool]
createUniqueBoolList nonZeroEdges inGraphList boolAccum =
  if null inGraphList then reverse $ fmap snd boolAccum
  else
    let firstGraph = head inGraphList
    in
    if null boolAccum then createUniqueBoolList nonZeroEdges (tail inGraphList) ((firstGraph,True) : boolAccum)
    else
        let checkList = filter (== True) $ fmap (topologicalEqual nonZeroEdges firstGraph) (fmap fst boolAccum)
        in
        if null checkList then createUniqueBoolList nonZeroEdges (tail inGraphList) ((firstGraph,True) : boolAccum)
        else createUniqueBoolList nonZeroEdges (tail inGraphList) ((firstGraph, False) : boolAccum)



-- | topologicalEqual takes two simple graphs and returns True if graphs have same nodes and edges
-- option to exclude zero weight edges
topologicalEqual :: Bool -> SimpleGraph -> SimpleGraph -> Bool
topologicalEqual nonZeroEdges g1 g2 =
  if LG.isEmpty g1 && LG.isEmpty g2 then True
  else if  LG.isEmpty g1 || LG.isEmpty g2 then False
  else
      let nodesG1 = LG.labNodes g1
          nodesG2 = LG.labNodes g2
          edgesG1 = if nonZeroEdges then fmap LG.toEdge $ filter ((> 0). thd3) $ LG.labEdges g1
                    else LG.edges g1
          edgesG2 = if nonZeroEdges then fmap LG.toEdge $ filter ((> 0). thd3) $ LG.labEdges g2
                    else LG.edges g2
      in
      if nodesG1 == nodesG2 && edgesG1 == edgesG2 then True
      else False

-- | getEdgeMinLengthToNode takes a labelled node and returns the min length of
-- the edge leading to the node
getEdgeMinLengthToNode ::[LG.LEdge EdgeInfo] ->  LG.LNode a -> Double
getEdgeMinLengthToNode  edgeList (node, _)=
  let foundEdge = L.find ((== node) . snd3) edgeList
  in
  -- root node will be nor be in in edge set and need so set > 0
  if foundEdge == Nothing then 1.0 --  error ("Edge not found in getEdgeMinLengthToNode: node " ++ (show node) ++ " edge list " ++ (show edgeList))
  else minLength $ thd3 $ fromJust foundEdge

-- | getBVUniqPhylogeneticGraph takes a list of phylogenetic graphs and returns
-- list of topologically unique graphs based on their node bitvector assignments
-- operatres on Decorated graph field
-- noZeroEdges flag passed to remove zero weight edges
getBVUniqPhylogeneticGraph :: Bool -> [PhylogeneticGraph] -> [PhylogeneticGraph]
getBVUniqPhylogeneticGraph nonZeroEdges inPhyloGraphList =
  if null inPhyloGraphList then []
  else
      let bvGraphList = fmap (getBVNodeList nonZeroEdges) $ fmap thd6 inPhyloGraphList
          uniqueBoolList = createBVUniqueBoolList bvGraphList []
          boolPair = zip inPhyloGraphList uniqueBoolList
      in
      fmap fst $ filter ((== True) . snd) boolPair


-- | getBVNodeList takes a DecoratedGraph and returns sorted list (by BV) of nodes
-- removes node with zero edge weight to them if specified
getBVNodeList :: Bool -> DecoratedGraph -> [BV.BitVector]
getBVNodeList nonZeroEdges inGraph =
  if LG.isEmpty inGraph then []
  else
      let nodeList =  LG.labNodes inGraph
          edgeList = LG.labEdges inGraph
          minLengthList = fmap (getEdgeMinLengthToNode edgeList) nodeList
          nodePairList = filter ((> 0) . snd) $ zip nodeList minLengthList
          bvNodeList  =  if nonZeroEdges then L.sort $ fmap bvLabel $ fmap snd $ fmap fst nodePairList
                         else L.sort $ fmap bvLabel $ fmap snd nodeList
          in
          bvNodeList

-- | createBVUniqueBoolList creates a list of Bool if graphs are unique by bitvecector node list
-- first occurrence is True, others False
-- assumes edges filterd b=y lenght already
createBVUniqueBoolList :: [[BV.BitVector]] -> [([BV.BitVector],Bool)] -> [Bool]
createBVUniqueBoolList inBVGraphListList boolAccum =
  if null inBVGraphListList then reverse $ fmap snd boolAccum
  else
    let firstGraphList = head inBVGraphListList
    in
    if null boolAccum then createBVUniqueBoolList  (tail inBVGraphListList) ((firstGraphList,True) : boolAccum)
    else
        let checkList = filter (== True) $ fmap (== firstGraphList) (fmap fst boolAccum)
        in
        if null checkList then createBVUniqueBoolList  (tail inBVGraphListList) ((firstGraphList,True) : boolAccum)
        else createBVUniqueBoolList  (tail inBVGraphListList) ((firstGraphList, False) : boolAccum)

-- | makeDummyLabEdge takes an unlabelled edge and adds a dummy label
makeDummyLabEdge :: EdgeInfo -> LG.Edge -> LG.LEdge EdgeInfo
makeDummyLabEdge edgeLab (u,v) = (u,v,edgeLab)

-- | netNodeAncestorViolation checks whether one of the edge into a netowrk node (in 2)
-- is cinnected to an ancestor (via the other parent) of the node
-- this is a form of time violation since the parents of a network node must be
-- at least possibly coeval
-- this uses the bit vector label of nodes.  If the other child of either parent node
-- of a network node has non-zero intersection between the BV label of the network node
-- and that other child of parent then they conecting edge is from and ancestral node hence a time violation
-- O(n) n netork nodes in Graph, but checks all nodes to see if network
hasNetNodeAncestorViolation :: LG.Gr VertexInfo b -> Bool
hasNetNodeAncestorViolation inGraph =
  if LG.isEmpty inGraph then error "Empty graph in hasNetNodeAncestorViolation"
  else
    let (_, _, _, netWorkNodeList) =  LG.splitVertexList inGraph
        hasAncViolationList = filter (== True) $ fmap (nodeAncViolation inGraph) netWorkNodeList
    in
    -- trace ("HNV: " ++ (show $ (not . null) hasAncViolationList))
    (not . null) hasAncViolationList

-- | nodeAncViolation checks a single node fo ancestrpo connection--he ceviolation
-- should be O(1).  Return True if violation
nodeAncViolation :: LG.Gr VertexInfo b -> LG.LNode VertexInfo -> Bool
nodeAncViolation inGraph inNode =
  let parentList = LG.labParents inGraph (fst inNode)
  in
  if length parentList /= 2 then error ("Parent number should be 2: " ++ (show $ fst inNode) ++ " <- " ++ (show $ fmap fst parentList))
  else
    let sisterNodes = concatMap (LG.sisterLabNodes inGraph) parentList
        sisterBVData = fmap (bvLabel . snd) sisterNodes
        inNodeBVData = bvLabel $ snd inNode
        sisterBVIntersections = fmap (.&. inNodeBVData) sisterBVData
        isAncInNode = filter (== inNodeBVData) sisterBVIntersections
    in
    (not . null) isAncInNode

-- | selectGraphStochastic takes a list of graphs and retuns a list of graphs chosen at Random
-- using an exponential distribution based on graph cost difference divided by an input factor
-- if factor is 0 then stringth graphs cost
-- mprob acceptance = -exp [(cost - minCost)/ factor]
-- returns n graphs by random criterion without replacment
selectGraphStochastic :: Int -> Int -> Double -> [PhylogeneticGraph] -> [PhylogeneticGraph]
selectGraphStochastic rSeed number factor inGraphList =
  if null inGraphList then inGraphList
  else if number >= length inGraphList then inGraphList
  else
    let randList' = randomIntList rSeed
        randList = fmap abs (tail randList')
        newSeed = head randList'
        minCost = minimum $ fmap snd6 inGraphList
        deltaList = fmap ((-) minCost) $ fmap snd6 inGraphList
        probAcceptList = fmap (getProb factor) deltaList

        -- multiplier for resolution 1000, 100 prob be ok
        randMultiplier = 1000
        randMultiplier' = fromIntegral randMultiplier
        intAcceptList = fmap floor $ fmap (* randMultiplier')  probAcceptList
        (_, intRandValList) = unzip $ zipWith divMod randList (replicate (length inGraphList) randMultiplier)
        acceptList = zipWith (<) intRandValList intAcceptList

        -- zip graphs with Bools
        (returnGraphList, _) = unzip $ filter ((== True) .snd) $ zip inGraphList acceptList

        -- takes some random remainder to fiill out length of list
        numLucky = number - (length returnGraphList)
        luckyList = if numLucky > 0 then takeRandom newSeed numLucky (fmap fst $ filter ((== False) .snd) $ zip inGraphList acceptList)
                    else []


    in
    trace ("SGS " ++ (show intAcceptList) ++ " " ++ (show intRandValList) ++ " -> " ++ (show acceptList))
    -- so no more than specified
    take number $ returnGraphList ++ luckyList

    where getProb a b = exp ((-1) * b / a)

-- | getDisplayTreeCostList returns a list of teh "block" costs of display trees
-- in a piar with any graph 'penalty' cost
getDisplayTreeCostList :: PhylogeneticGraph -> ([VertexCost], VertexCost)
getDisplayTreeCostList inGraph =
  if LG.isEmpty $ thd6 inGraph then ([], 0.0)
  else
    let rootIndex = fst $ head $ LG.getRoots $ fst6 inGraph
        displayTreeCharVect = fft6 inGraph
        displayTreeCostVect = fmap (getBlockCost rootIndex) displayTreeCharVect
        nonGraphCost = V.sum  displayTreeCostVect
    in
    (V.toList displayTreeCostVect, (snd6 inGraph) - nonGraphCost)

-- | getBlockCost returns the cost, summed over characters, of a character block
getBlockCost :: LG.Node -> V.Vector DecoratedGraph -> VertexCost
getBlockCost rootIndex charGraphVect = 
  if V.null charGraphVect then 0.0
  else 
    V.sum $ fmap (getCharacterCost rootIndex) charGraphVect

-- | getCharacterCost returns charcter cost as root of character tree
getCharacterCost :: LG.Node -> DecoratedGraph -> VertexCost
getCharacterCost rootIndex inGraph =
  if LG.isEmpty inGraph then 0.0
  else 
    let rootLabel = LG.lab inGraph rootIndex
    in
    if isNothing rootLabel then error ("Root index without label: " ++ (show rootIndex))
    else subGraphCost $ fromJust rootLabel

