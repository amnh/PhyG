{-# OPTIONS_GHC -Wno-missed-specialisations #-}

{- |
Module specifying general graph functions--with types specific to Types.hs
graph functions that a re general are in LocalGraph.hs
-}
module Graphs.GraphOperations (
    contractIn1Out1EdgesRename,
    convertDecoratedToSimpleGraph,
    convertDecoratedToSimpleGraphBranchLength,
    convertGeneralGraphToPhylogeneticGraph,
    convertPhylogeneticGraph2Reduced,
    convertReduced2PhylogeneticGraph,
    convertReduced2PhylogeneticGraphSimple,
    convertSimpleToDecoratedGraph,
    convertToSimpleEdge,
    copyIAFinalToPrelim,
    copyIAPrelimToFinal,
    dichotomizeRoot,
    getBVUniqPhylogeneticGraph,
    getDecoratedDisplayTreeList,
    getDisplayTreeCostList,
    getNodeType,
    getTopoUniqPhylogeneticGraph,
    getUniqueGraphs,
    getSoftWiredNodeSubGraphCost,
    getSoftWiredNodeSubGraphCostWithVertex,
    graphCostFromNodes,
    hasNetNodeAncestorViolation,
    isNovelGraph,
    isPhylogeneticDecoratedGraph,
    ladderizeGraph,
    makeDummyLabEdge,
    makeGraphTimeConsistent,
    makeIAFinalFromPrelim,
    makeIAPrelimFromFinal,
    makeLeafGraph,
    makeNewickList,
    makeSimpleLeafGraph,
    minBlockResolutionCost,
    parentsInChainGraph,
    phylogeneticGraphListMinus,
    reducedphylogeneticGraphListMinus,
    remakePhylogeneticGraph,
    removeParentsInChain,
    renameSimpleGraphNodes,
    renameSimpleGraphNodesString,
    selectGraphStochastic,
    selectGraphs,
    selectGraphsFull,
    selectPhylogeneticGraph,
    selectPhylogeneticGraphReduced,
    showDecGraphs,
    sortEdgeListByLength,
    topologicalEqual,
) where

import Bio.DynamicCharacter
import Commands.Verify qualified as V
import Control.Monad.Random (filterM)
import Control.Monad.Random.Class
import Data.BitVector.LittleEndian qualified as BV
import Data.Bits
import Data.Char qualified as C
import Data.Functor ((<&>))
import Data.List qualified as L
import Data.Maybe
import Data.Ord
import Data.Text.Lazy qualified as T
import Data.Vector qualified as V
import Data.Vector.Generic qualified as GV
import GeneralUtilities
import GraphFormatUtilities qualified as GFU
import GraphOptimization.Medians qualified as M
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (Computing))
import PHANE.Evaluation.Verbosity (Verbosity (..))
import Text.Read
import Types.Types
import Utilities.LocalGraph qualified as LG
import Utilities.Utilities qualified as U
import Debug.Trace

-- | getSoftWiredNodeSubGraphCost  takes a softwired graph and returns the sub graph cost 
-- for that node based on minimium of resolutoncache costs--so underestimate likely

-- Blocks
-- vertexResolutionData ∷ V.Vector ResolutionBlockData 

-- Characters within Blocks
-- ResolutionBlockData = V.Vector ResolutionData

getSoftWiredNodeSubGraphCost:: ReducedPhylogeneticGraph -> Int -> VertexCost
getSoftWiredNodeSubGraphCost inPhyloGraph vertex =
    if LG.isEmpty $ thd5 inPhyloGraph then infinity
    else
        let labelDataMaybe = LG.lab (thd5 inPhyloGraph) vertex
        in
        if isNothing labelDataMaybe then error ("Vertex not in graph: " <> (show vertex))
        else 
            let resData = vertexResolutionData $ fromJust labelDataMaybe
            in
            V.sum $ fmap minBlockResolutionCost resData 

-- | minBlockkResolutionCost returns the average cost subgraph resolutoin cost of a character
minBlockResolutionCost :: ResolutionBlockData -> VertexCost
minBlockResolutionCost inBlockResData = 
    -- 0.5 *  (V.minimum $ fmap minCharResolutionCost inBlockResData) + (V.maximum $ fmap minCharResolutionCost inBlockResData) 
    let sumAll = V.sum $ fmap displayCost inBlockResData
        {-
        trimmedList = if V.length inBlockResData >3 then
                        let first = L.take (V.length inBlockResData -1) $ L.sort $ V.toList $ fmap displayCost inBlockResData
                        in
                        L.drop 1 first
                      else if V.length inBlockResData == 3 then
                         L.drop 1 $ L.sort $ V.toList $ fmap displayCost inBlockResData
                      else V.toList $ fmap displayCost inBlockResData
        sumAll' = L.sum trimmedList 
        -}
    in
    --sumAll' / (fromIntegral $ L.length trimmedList)
    sumAll / (fromIntegral $ V.length inBlockResData)


{-
These functions try to look at complementary resolution costs
to get heuristic but was a large overestimate of delta
-}

-- | getSoftWiredNodeSubGraphCostWithVertex  takes a softwired graph and returns the sub graph cost 
-- for that node based on minimium of resolutoncache costs--so underestimate likely

-- Blocks
-- vertexResolutionData ∷ V.Vector ResolutionBlockData 

-- Characters within Blocks
-- ResolutionBlockData = V.Vector ResolutionData

getSoftWiredNodeSubGraphCostWithVertex:: ReducedPhylogeneticGraph -> Int -> Int -> (VertexCost, VertexCost)
getSoftWiredNodeSubGraphCostWithVertex inPhyloGraph netVertex vertex =
    if LG.isEmpty $ thd5 inPhyloGraph then (infinity, infinity)
    else
        let labelDataMaybe = LG.lab (thd5 inPhyloGraph) vertex
        in
        if isNothing labelDataMaybe then error ("Vertex not in graph: " <> (show vertex))
        else 
            let resData = vertexResolutionData $ fromJust labelDataMaybe
                (hasNetVertexList, noNetVertexList) =  V.unzip $ fmap (minBlockResolutionCostWithVertex netVertex) resData 
            in
            (V.sum hasNetVertexList, V.sum noNetVertexList)

-- | minBlockResolutionCostWithVertex returns the average cost subgraph resolutoin cost of a character
minBlockResolutionCostWithVertex :: Int -> ResolutionBlockData -> (VertexCost, VertexCost)
minBlockResolutionCostWithVertex netVertex inBlockResData  = 
    -- 0.5 *  (V.minimum $ fmap minCharResolutionCost inBlockResData) + (V.maximum $ fmap minCharResolutionCost inBlockResData) 
    let (withVertexList, noVertexList) = V.unzip $ fmap (minCharResolutionCostWithVertex netVertex) inBlockResData
        sumHasVertex = V.sum withVertexList
        sumNoVertex  = V.sum noVertexList
    in
    (V.minimum withVertexList, V.minimum noVertexList)
    --(sumHasVertex / (fromIntegral $ V.length withVertexList), sumNoVertex / (fromIntegral $ V.length noVertexList))

-- | minCharResolutionCostWithVertex returns the minimum cost subgraph resolution cost of a character
minCharResolutionCostWithVertex :: Int -> ResolutionData -> (VertexCost, VertexCost)
minCharResolutionCostWithVertex netVertex inResData = 
    let subgraphNodes = getNodeIndices $ displaySubGraph inResData
        hasNetVertex = netVertex `elem` subgraphNodes
    in
    if hasNetVertex then
        (displayCost inResData, 0.0)
    else 
        (0.0, displayCost inResData)
    where getNodeIndices (a,_) = fmap fst a

 
-- | convertPhylogeneticGraph2Reduced takes a Phylogenetic graph and returns a reduced phylogenetiv graph
convertPhylogeneticGraph2Reduced ∷ PhylogeneticGraph → ReducedPhylogeneticGraph
convertPhylogeneticGraph2Reduced inPhyloGraph@(a, b, c, displayTreeV, _, f) =
    if inPhyloGraph == emptyPhylogeneticGraph
        then emptyReducedPhylogeneticGraph
        else
            let newDisplayTreeV = fmap (fmap convertDecoratedToSimpleGraph) displayTreeV
            in  (a, b, c, newDisplayTreeV, f)


{- | convertReduced2PhylogeneticGraphSimple just adds in a correctly typed fifth field
bit decorartions are not extractged
-}
convertReduced2PhylogeneticGraphSimple ∷ ReducedPhylogeneticGraph → PhylogeneticGraph
convertReduced2PhylogeneticGraphSimple (a, b, c, d, f) =
    let newDisplayTreeV = fmap (fmap convertSimpleToDecoratedGraph) d
    in  (a, b, c, newDisplayTreeV, mempty, f)


-- | convertReduced2PhylogeneticGraph takes a reduced phylogenetic graph and returns a phylogenetiv graph
convertReduced2PhylogeneticGraph ∷ ReducedPhylogeneticGraph → PhylogeneticGraph
convertReduced2PhylogeneticGraph inReducedPhyloGraph@(a, b, canonicalGraph, displayTreeV, charInfoVV) =
    if inReducedPhyloGraph == emptyReducedPhylogeneticGraph
        then emptyPhylogeneticGraph
        else
            let newDisplayTreeVL = fmap (getDecoratedDisplayTreeList canonicalGraph) $ V.zip displayTreeV $ V.fromList [0 .. (V.length charInfoVV - 1)]
                newCharacterTreeVV = fmap getCharacterTree $ V.zip (fmap head newDisplayTreeVL) (fmap V.length charInfoVV)
            in  (a, b, canonicalGraph, newDisplayTreeVL, newCharacterTreeVV, charInfoVV)


{- | getCharacterTree takes an input Decotate disply tree list and returns a vectopr of character trees
based on the first display tree in list
-}
getCharacterTree ∷ (DecoratedGraph, Int) → V.Vector DecoratedGraph
getCharacterTree (inDisplayTree, numCharTrees) =
    if LG.isEmpty inDisplayTree
        then error "Empty display tree getCharaterTree"
        else
            let displayNodeList = LG.labNodes inDisplayTree
                displayEdgeList = LG.labEdges inDisplayTree
                charNodeListList = fmap (makeCharacterNodes displayNodeList) [0 .. numCharTrees - 1]
            in  V.fromList $ zipWith LG.mkGraph charNodeListList (L.replicate numCharTrees displayEdgeList)


-- | makeCharacterNodes makes character nodes from dispaly tree nodes
makeCharacterNodes ∷ [LG.LNode VertexInfo] → Int → [LG.LNode VertexInfo]
makeCharacterNodes displayTreeNodeList charIndex =
    if null displayTreeNodeList
        then error "Null node list in extractCharacterNodeInfo"
        else fmap (extractCharNodeInfo charIndex) displayTreeNodeList


-- | extractCharNodeInfo takes a charcter indfex and produces a node with newinfo of single block and single character
extractCharNodeInfo ∷ Int → LG.LNode VertexInfo → LG.LNode VertexInfo
extractCharNodeInfo charIndex (a, cLabel) =
    let charVertData = V.singleton $ V.singleton $ (V.head (vertData cLabel)) V.! charIndex
        charVertexResolutionData = V.singleton $ V.singleton $ (V.head (vertexResolutionData cLabel)) V.! charIndex
        newLabel =
            cLabel
                { vertData = charVertData
                , vertexResolutionData = charVertexResolutionData
                }
    in  (a, newLabel)


{- | getDecoratedDisplayTreeList takes a cononical graph and pair of list of display trees as ASimpleGraph and an index
and creates a new display tree list from the canonical infomation
all display treelist nodes are assigned same decorations
edge weights are not recalculated
-}
getDecoratedDisplayTreeList ∷ DecoratedGraph → ([SimpleGraph], Int) → [DecoratedGraph]
getDecoratedDisplayTreeList inCanonicalGraph (inDisplayL, blockIndex)
    | LG.isEmpty inCanonicalGraph = error "Empty canonical graph in getDecoratedDisplayTree"
    | otherwise =
        -- get edges and nodes for canonical graph
        let canNodList = LG.labNodes inCanonicalGraph
            canEdgedList = LG.labEdges inCanonicalGraph
            newNodeListList = fmap (makeDisplayNodes canNodList blockIndex) (fmap LG.nodes inDisplayL)
            newDisplayEdgeList = fmap (makeDisplayEdges canEdgedList) inDisplayL
        in  zipWith LG.mkGraph newNodeListList newDisplayEdgeList


{- | makeDisplayNodes takes canincal node label data and pulls out specific block
infomartion and cretes new nodes with single block data only
assumes nodes have same indices
-}
makeDisplayNodes ∷ [LG.LNode VertexInfo] → Int → [LG.Node] → [LG.LNode VertexInfo]
makeDisplayNodes canonicalNodeList blockIndex inDisplayNodeList
    | length canonicalNodeList /= length inDisplayNodeList = error "Canonical and display node lists are unequal in length"
    | otherwise = extractBlockNodeInfo blockIndex <$> canonicalNodeList


-- | extractBlockNodeInfo takes a single node and block inde and extracts block data to make a new, single block
extractBlockNodeInfo ∷ Int → LG.LNode VertexInfo → LG.LNode VertexInfo
extractBlockNodeInfo blockIndex (a, cLabel) =
    let blockVertData = V.singleton $ (vertData cLabel) V.! blockIndex
        blockVertexResolutionData = V.singleton $ (vertexResolutionData cLabel) V.! blockIndex
        newLabel =
            cLabel
                { vertData = blockVertData
                , vertexResolutionData = blockVertexResolutionData
                }
    in  (a, newLabel)


{- | makeDisplayEdges tcretes edges in display tree
by taking corresponding edge in canonical tree
edge weight info is not recalculated and is likely incorrect
-}
makeDisplayEdges ∷ [LG.LEdge EdgeInfo] → SimpleGraph → [LG.LEdge EdgeInfo]
makeDisplayEdges canonicalEdgeList displayTree =
    let displayEdgeList = LG.labEdges displayTree
    in  fmap (copyLabelInfo canonicalEdgeList) displayEdgeList


-- | copyLabelInfo finds corresponding edge (undirected) and takes relevant data from it to make new labelled edge
copyLabelInfo ∷ [LG.LEdge EdgeInfo] → LG.LEdge Double → LG.LEdge EdgeInfo
copyLabelInfo canonicalEdgeList displayEdge@(u, v, dl) =
    if null canonicalEdgeList
        then error "Empty edge list in copyLabelInfo"
        else
            let canonicalEdge = L.find (matchIndices displayEdge) canonicalEdgeList
                newLabel =
                    (thd3 $ fromJust canonicalEdge)
                        { minLength = dl
                        , maxLength = dl
                        , midRangeLength = dl
                        }
            in  if isNothing canonicalEdge
                    then error "Cannot find canonical edge in copyLabelInfo"
                    else (u, v, newLabel)
    where
        matchIndices (a, b, _) (a', b', _) =
            if a == a' && b == b'
                then True
                else
                    if a == b' && b == a'
                        then True
                        else False


{- | isPhylogeneticDecoratedGraph checks various issues to see if
there is wierdness in graph
-}
isPhylogeneticDecoratedGraph ∷ DecoratedGraph → PhyG Bool
isPhylogeneticDecoratedGraph inGraph =
    if LG.isEmpty inGraph
        then pure False
        else
            let nodeList = fmap fst $ LG.labNodes inGraph
                indegreeList = fmap (LG.inn inGraph) nodeList
                outdegreeList = fmap (LG.out inGraph) nodeList
            in  if LG.hasDuplicateEdgesNub inGraph
                    then pure False
                    else
                        if length (LG.getRoots inGraph) /= 1
                            then pure False
                            else
                                if LG.outdeg inGraph (head $ LG.getRoots inGraph) /= 2
                                    then pure False
                                    else
                                        if (not . null) (LG.getIsolatedNodes inGraph)
                                            then pure False
                                            else
                                                if (not . null) (filter ((> 2) . length) indegreeList)
                                                    then pure False
                                                    else
                                                        if (not . null) (filter ((> 2) . length) outdegreeList)
                                                            then pure False
                                                            else
                                                                if parentsInChainGraph inGraph
                                                                    then pure False
                                                                    else do
                                                                        consistent ← LG.isGraphTimeConsistent inGraph
                                                                        if not consistent
                                                                            then pure False
                                                                            else pure True


{- | parentsInChainGraph checks all network vertices in graph to see if
parents in any network vertex is ancestor of other
-}
parentsInChainGraph ∷ DecoratedGraph → Bool
parentsInChainGraph inGraph =
    if LG.isEmpty inGraph
        then False
        else
            let (_, _, _, netVertexList) = LG.splitVertexList inGraph
                parentChainList = fmap (parentsInChainVertex inGraph) $ fmap fst netVertexList
            in  -- if ((not . null) $ filter (== True) parentChainList) then traceNoLF ("NPDG+ ") $ ((not . null) $ filter (== True) parentChainList)
                -- else  traceNoLF ("NPDG- ") $ (not . null) $ filter (== True) parentChainList
                (not . null) $ filter (== True) parentChainList


{- | parentsInChainVertex checks for a network vertex if
one parent is ancestor of other
retuens false if not net vertex
-}
parentsInChainVertex ∷ DecoratedGraph → LG.Node → Bool
parentsInChainVertex inGraph inNode =
    if LG.isEmpty inGraph
        then False
        else
            let parentNodes = LG.labParents inGraph inNode
                firstParent = head parentNodes
                secondParent = last parentNodes
                firstBV = bvLabel $ snd firstParent
                secondBV = bvLabel $ snd secondParent
                oneAncOther = (firstBV .&. secondBV) `elem` [firstBV, secondBV]
            in  if length parentNodes < 2
                    then False
                    else
                        if oneAncOther
                            then True
                            else -- traceNoLF ("PCV:" <> (show (firstBV, secondBV, firstBV .&. secondBV)))
                                False


{- | ReducedphylogeneticGraphListMinus subtracts the secoind argiument list from first
if an element is multiple times in firt list each will be removed
equality comparison is based on String rep of graphs vertes and edges (prettyVertices)
does not take cost into account--or edge weight--only topology.
result like (minuendList - subtrahendList)
-}
reducedphylogeneticGraphListMinus ∷ [ReducedPhylogeneticGraph] → [ReducedPhylogeneticGraph] → [ReducedPhylogeneticGraph]
reducedphylogeneticGraphListMinus minuendList subtrahendList
    | null minuendList = []
    | null subtrahendList = minuendList
    | otherwise =
        let minuendSimpleStringList = fmap (LG.prettyIndices . fst5) minuendList
            subtrahendSinpleStringList = fmap (LG.prettyIndices . fst5) subtrahendList
            inSubtrahendList = fmap (`elem` subtrahendSinpleStringList) minuendSimpleStringList
            differenceList = fmap fst $ filter (not . snd) $ zip minuendList inSubtrahendList
        in  differenceList


{- | phylogeneticGraphListMinus subtracts the secoind argiument list from first
if an element is multiple times in firt list each will be removed
equality comparison is based on String rep of graphs vertes and edges (prettyVertices)
does not take cost into account--or edge weight--only topology.
result like (minuendList - subtrahendList)
-}
phylogeneticGraphListMinus ∷ [PhylogeneticGraph] → [PhylogeneticGraph] → [PhylogeneticGraph]
phylogeneticGraphListMinus minuendList subtrahendList
    | null minuendList = []
    | null subtrahendList = minuendList
    | otherwise =
        let minuendSimpleStringList = fmap (LG.prettyIndices . fst6) minuendList
            subtrahendSinpleStringList = fmap (LG.prettyIndices . fst6) subtrahendList
            inSubtrahendList = fmap (`elem` subtrahendSinpleStringList) minuendSimpleStringList
            differenceList = fmap fst $ filter (not . snd) $ zip minuendList inSubtrahendList
        in  differenceList


-- | makeNewickList takes a list of fgl trees and outputs a single String cointaining the graphs in Newick format
makeNewickList ∷ Bool → Bool → Bool → Int → [SimpleGraph] → [VertexCost] → String
makeNewickList isTNT writeEdgeWeight writeNodeLabel' rootIndex graphList costList =
    let allTrees = L.foldl' (&&) True (fmap LG.isTree graphList)

        -- check for network HTU label requirement
        writeNodeLabel
            | allTrees = writeNodeLabel'
            | writeNodeLabel' = writeNodeLabel'
            | otherwise -- trace "HTU labels are required for ENewick Output"
                =
                True

        graphString = GFU.fglList2ForestEnhancedNewickString (fmap (LG.rerootTree rootIndex) graphList) writeEdgeWeight writeNodeLabel
        newickStringList = fmap init $ filter (not . null) $ lines graphString
        costStringList =
            if not isTNT
                then fmap ((('[' :) . (<> "];\n")) . show) costList
                else L.replicate (length costList) "*\n"

        graphStringCost = concat $ zipWith (<>) newickStringList costStringList
        graphStringCostTNT = fmap commaToSpace graphStringCost

        graphStringCost' =
            if not isTNT
                then graphStringCost
                else (take (length graphStringCostTNT - 2) graphStringCostTNT) <> ";\n"
    in  --trace ("MNL: " <> (show (writeNodeLabel, isTNT, writeEdgeWeight, writeNodeLabel',graphStringCost)))
        graphStringCost'
    where
        commaToSpace a = if a /= ',' then a else ' '


{- | convertGeneralGraphToPhylogeneticGraph inputs a SimpleGraph and converts it to a Phylogenetic graph by:
 1) transitive reduction -- removes anc <-> desc netork edges
 2) ladderizes -- all vertices are (in degree, outdegree) (0,1|2) (1,2) (2,1) (1,0)
       by adding extra HTIs and edges
 3) checks time consistency and removes edges stepwise from
       those that violate the most  before/after splits of network edges
       arbitrary but deterministic
 4) contracts out any remaning indegree 1 outdegree 1 nodes and renames HTUs in order
these tests can be screwed up by imporperly formated graphs comming in (self edges, chained network edge etc)
-}
convertGeneralGraphToPhylogeneticGraph ∷ Bool → SimpleGraph → PhyG SimpleGraph
convertGeneralGraphToPhylogeneticGraph correct inGraph =
    if LG.isEmpty inGraph
        then pure LG.empty
        else
            let -- remove single "tail" edge from root with single child, replace child node with root
                noTailGraph = LG.contractRootOut1Edge inGraph

                -- remove non-leaf nodes (index > root) with outdegree 0
                nonNonLeafOut0 = LG.removeNonLeafOut0NodesAfterRoot noTailGraph

                -- remove indeg 1 out deg 1 edges
                noIn1Out1Graph = contractIn1Out1EdgesRename nonNonLeafOut0 -- noTailGraph

                -- transitive reduction
                -- only wanted to EUN and CUN--but they do it
                -- reducedGraph = LG.transitiveReduceGraph noIn1Out1Graph

                -- caused problems at one point
                -- laderization of indegree and outdegree edges
                ladderGraph = ladderizeGraph noIn1Out1Graph -- reducedGraph
            in  do
                    -- time consistency (after those removed by transitrive reduction)
                    timeConsistentGraph ← makeGraphTimeConsistent correct ladderGraph

                    -- removes parent child network edges
                    let noChainedGraph = LG.removeChainedNetworkNodes False timeConsistentGraph

                    -- removes ancestor descendent edges transitiveReduceGraph should do this
                    -- but that looks at all nodes not just vertex
                    let noParentChainGraph = removeParentsInChain correct $ fromJust noChainedGraph -- timeConsistentGraph --

                    -- deals the nodes with all network children
                    let noTreeEdgeGraph = LG.removeTreeEdgeFromTreeNodeWithAllNetworkChildren noParentChainGraph

                    -- remove sister-sister edge.  where two network nodes have same parents
                    let noSisterSisterGraph = removeSisterSisterEdges correct noTreeEdgeGraph

                    -- remove and new zero nodes
                    let finalGraph = LG.removeNonLeafOut0NodesAfterRoot noSisterSisterGraph

                    if LG.isEmpty timeConsistentGraph
                        then pure LG.empty
                        else
                            if isNothing noChainedGraph
                                then pure LG.empty
                                else
                                    if LG.isEmpty noParentChainGraph
                                        then pure LG.empty
                                        else
                                            if LG.isEmpty noSisterSisterGraph
                                                then pure LG.empty
                                                else -- trace ("CGP orig:\n" <> (LG.prettify inGraph) <> "\nNew:" <> (LG.prettify timeConsistentGraph))
                                                -- cycle check to make sure--can be removed when things working
                                                -- else if LG.cyclic noSisterSisterGraph then error ("Cycle in graph : \n" <> (LG.prettify noSisterSisterGraph))

                                                -- this final need to ladderize or recontract?

                                                    if finalGraph == inGraph
                                                        then pure finalGraph
                                                        else convertGeneralGraphToPhylogeneticGraph correct finalGraph


-- | removeParentsInChain checks the parents of each netowrk node are not anc/desc of each other
removeParentsInChain ∷ Bool → SimpleGraph → SimpleGraph
removeParentsInChain correct inGraph =
    if LG.isEmpty inGraph
        then LG.empty
        else
            let (_, _, _, netVertexList) = LG.splitVertexList inGraph
                parentNetVertList = fmap (LG.labParents inGraph . fst) netVertexList

                -- get list of nodes that are transitively equal in age
                concurrentList = LG.mergeConcurrentNodeLists parentNetVertList []
                concurrentPairList = concatMap getListPairs concurrentList

                -- get pairs that violate concurrency
                violatingConcurrentPairs = concatMap (LG.concurrentViolatePair inGraph) concurrentPairList

                -- get network nodes with violations
                parentNodeViolateList = concatMap pairToList violatingConcurrentPairs
                childNodeViolateList = concatMap (LG.descendants inGraph) parentNodeViolateList
                netNodeViolateList = filter (LG.isNetworkNode inGraph) childNodeViolateList

                netEdgesThatViolate = fmap LG.toEdge $ LG.inn inGraph $ head netNodeViolateList
            in  if null violatingConcurrentPairs
                    then inGraph
                    else
                        if null netNodeViolateList
                            then error "Should be neNode that violate"
                            else
                                if null netEdgesThatViolate
                                    then error "Should be violating in edges"
                                    else
                                        if not correct
                                            then LG.empty
                                            else
                                                let edgeDeletedGraph = LG.delEdge (head netEdgesThatViolate) inGraph
                                                    newGraph = contractIn1Out1EdgesRename edgeDeletedGraph
                                                in  -- trace ("PIC")
                                                    removeParentsInChain correct newGraph
    where
        pairToList (a, b) = [fst a, fst b]


{- | removeSisterSisterEdges takes a graph and recursively removes a single edge fomr where two network
edges have the same two parents
-}
removeSisterSisterEdges ∷ Bool → SimpleGraph → SimpleGraph
removeSisterSisterEdges correct inGraph =
    if LG.isEmpty inGraph
        then LG.empty
        else
            let sisterSisterEdges = LG.getSisterSisterEdgeList inGraph
                -- newGraph = LG.delEdge (head sisterSisterEdges) inGraph
                newGraph = LG.delEdges sisterSisterEdges inGraph
                newGraph' = contractIn1Out1EdgesRename newGraph
            in  if null sisterSisterEdges
                    then inGraph
                    else
                        if not correct
                            then LG.empty
                            else -- trace ("Sister")
                            -- removeSisterSisterEdges
                                newGraph'


{- | makeGraphTimeConsistent takes laderized, transitive reduced graph and deletes
network edges in an arbitrary but deterministic sequence to produce a phylogentic graphs suitable
for swapping etc
looks for violation of time between netork edges based on "before" and "after"
tests of nodes that should be potentially same age
removes second edge of second pair of two network edges in each case adn remakes graph
strict paralle since will need each recursive run
-}
makeGraphTimeConsistent ∷ Bool → SimpleGraph → PhyG SimpleGraph
makeGraphTimeConsistent correct inGraph
    | LG.isEmpty inGraph = pure LG.empty
    | LG.isTree inGraph = pure inGraph
    | otherwise =
        let coevalNodeConstraintList = LG.coevalNodePairs inGraph
            -- parallel setup
            -- action :: (Show a,Eq a,Eq b) => (LG.LNode a, LG.LNode a) -> (LG.LNode a, LG.LNode a, [LG.LNode a], [LG.LNode a], [LG.LNode a], [LG.LNode a])
            action = LG.addBeforeAfterToPair inGraph
        in  do
                pTraverse ← getParallelChunkMap
                let coevalNodeConstraintList' = pTraverse action coevalNodeConstraintList
                -- PU.seqParMap PU.myStrategyRDS (LG.addBeforeAfterToPair inGraph) coevalNodeConstraintList -- `using`  PU.myParListChunkRDS
                let coevalPairsToCompareList = getListPairs coevalNodeConstraintList'
                let timeOffendingEdgeList = LG.getEdgesToRemoveForTime inGraph coevalPairsToCompareList
                let newGraph = LG.delEdges timeOffendingEdgeList inGraph
                -- trace ("MGTC:" <> (show timeOffendingEdgeList))
                if (not correct) && (not . null) timeOffendingEdgeList
                    then pure LG.empty
                    else pure $ contractIn1Out1EdgesRename newGraph


{- | contractIn1Out1EdgesRename contracts in degree and outdegree edges and renames HTUs in index order
does one at a time and makes a graph and recurses
-}
contractIn1Out1EdgesRename ∷ SimpleGraph → SimpleGraph
contractIn1Out1EdgesRename inGraph =
    if LG.isEmpty inGraph
        then LG.empty
        else
            let newGraph = LG.contractIn1Out1Edges inGraph
            in  renameSimpleGraphNodes newGraph


-- | renameSimpleGraphNodes takes nodes and renames HTU nodes based on index
renameSimpleGraphNodes ∷ SimpleGraph → SimpleGraph
renameSimpleGraphNodes inGraph
    | LG.isEmpty inGraph = LG.empty
    | otherwise =
        let inNodes = LG.labNodes inGraph
            nodeLabels = fmap (makeSimpleLabel inGraph) inNodes
            newNodes = zip (fmap fst inNodes) nodeLabels
            newEdges = LG.labEdges inGraph
            makeSimpleLabel g (a, b)
                | not $ LG.isLeaf g a = T.pack $ "HTU" <> show a
                | otherwise = b
        in  LG.mkGraph newNodes newEdges


-- | renameSimpleGraphNodesString takes nodes and renames HTU nodes based on index
renameSimpleGraphNodesString ∷ LG.Gr String String → LG.Gr String String
renameSimpleGraphNodesString inGraph =
    if LG.isEmpty inGraph
        then LG.empty
        else
            let inNodes = LG.labNodes inGraph
                nodeLabels = fmap (makeSimpleLabel inGraph) inNodes
                newNodes = zip (fmap fst inNodes) nodeLabels
                newEdges = LG.labEdges inGraph
            in  -- newGraph
                -- trace ("C11: " <> (show $ LG.getIsolatedNodes newGraph) <> " => " <> (show newNodes) <> " " <> (show $ fmap LG.toEdge newEdges))
                LG.mkGraph newNodes newEdges
    where
        makeSimpleLabel g (a, b) =
            if not $ LG.isLeaf g a
                then "HTU" <> show a
                else b


-- | sortEdgeListByLength sorts edge list by length (midRange), highest to lowest
sortEdgeListByLength ∷ [LG.LEdge EdgeInfo] → [LG.LEdge EdgeInfo]
sortEdgeListByLength inEdgeList =
    if null inEdgeList
        then []
        else L.sortOn (Data.Ord.Down . midRangeLength . thd3) inEdgeList


{- | ladderizeGraph is a wrapper around ladderizeGraph' to allow for mapping with
local nodelist
-}
ladderizeGraph ∷ SimpleGraph → SimpleGraph
ladderizeGraph inGraph = ladderizeGraph' inGraph (LG.nodes inGraph)


{- | ladderize takes an input graph and ensures/creates nodes
such that all vertices are (indegree, outdegree) (0,>0), (1,2) (2,1) (1,0)
ladderizeGraph' :: SimpleGraph -> [LG.Node] -> SimpleGraph
-}
ladderizeGraph' ∷ SimpleGraph → [LG.Node] → SimpleGraph
ladderizeGraph' inGraph nodeList
    | LG.isEmpty inGraph = LG.empty
    | null nodeList = inGraph
    | otherwise =
        let -- these are roots, network, tree, leaf nodes
            okNodeDegrees = [(0, 2), (1, 2), (2, 1), (1, 0)]
            firstNode = head nodeList
            (inEdgeList, outEdgeList) = LG.getInOutEdges inGraph firstNode
            inOutPairLength = (length inEdgeList, length outEdgeList)
        in  -- trace ("node " <> (show firstNode) <> " " <> (show inOutPairLength)) (
            -- node ok to keep
            if inOutPairLength `elem` okNodeDegrees
                then ladderizeGraph' inGraph (tail nodeList)
                else -- node edges need modification

                    let newGraph = resolveNode inGraph firstNode (inEdgeList, outEdgeList) inOutPairLength
                    in  -- trace ("resolving " <> "node " <> (show firstNode) <> " " <> (show inOutPairLength) )
                        ladderizeGraph' newGraph (LG.nodes newGraph)


-- )

{- | resolveNode takes a graph and node and inbound edgelist and outbound edge list
and converts node to one of (indeg, outdeg) (0,1),(0,2),(1,2),(2,1),(1,0)
this only resolves a single nodes edges at a time and then returns new graph
when more hase to be done--that will occur on lultiple passes through nodes.
perhaps not the most efficient, but only done once per input graph
contracts indegree 1 outdegree 1 nodes
-}
resolveNode ∷ SimpleGraph → LG.Node → ([LG.LEdge Double], [LG.LEdge Double]) → (Int, Int) → SimpleGraph
resolveNode inGraph curNode inOutPair@(inEdgeList, outEdgeList) (inNum, outNum) =
    if LG.isEmpty inGraph
        then LG.empty
        else -- trace ("Resolving " <> show (inNum, outNum)) (

            let numNodes = length $ LG.nodes inGraph
            in  -- isolated node -- throw warning and delete
                if inNum == 0 && outNum == 0
                    then -- this is  debug info
                    -- trace ("Warning: ResolveNode deleting isolated vertex " <> show curNode) ( --  <> " in graph\n" <> LG.prettify inGraph )

                        let newGraph = LG.delNode curNode inGraph
                        in  newGraph
                    else -- node with too many parents and too many children
                    -- converts to tree node--biased in that direction

                        if (inNum > 2) && (outNum > 2)
                            then
                                let first2Edges = take 2 inEdgeList
                                    newNode = (numNodes, T.pack ("HTU" <> show numNodes))
                                    newEdge1 = (fst3 $ head first2Edges, numNodes, 0.0 ∷ Double)
                                    newEdge2 = (fst3 $ last first2Edges, numNodes, 0.0 ∷ Double)
                                    newEdge3 = (numNodes, curNode, 0.0 ∷ Double)
                                    newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
                                in  newGraph
                            else -- leaf leaf with too many parents

                                if (inNum > 1) && (outNum == 0) || (inNum > 2) && (outNum == 1) || (inNum > 1) && (outNum == 2)
                                    then
                                        let first2Edges = take 2 inEdgeList
                                            newNode = (numNodes, T.pack ("HTU" <> show numNodes))
                                            newEdge1 = (fst3 $ head first2Edges, numNodes, 0.0 ∷ Double)
                                            newEdge2 = (fst3 $ last first2Edges, numNodes, 0.0 ∷ Double)
                                            newEdge3 = (numNodes, curNode, 0.0 ∷ Double)
                                            newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
                                        in  newGraph
                                    else -- indegree 1 outdegree 1 node to contract

                                        if inNum == 1 && outNum == 1
                                            then
                                                let newEdge = (fst3 $ head inEdgeList, snd3 $ head outEdgeList, 0.0 ∷ Double)
                                                    newGraph = LG.insEdge newEdge $ LG.delNode curNode $ LG.delLEdges (inEdgeList <> outEdgeList) inGraph
                                                in  newGraph
                                            else
                                                if inNum < 2 || outNum > 2
                                                    then
                                                        let first2Edges = take 2 outEdgeList
                                                            newNode = (numNodes, T.pack ("HTU" <> show numNodes))
                                                            newEdge1 = (numNodes, snd3 $ head first2Edges, 0.0 ∷ Double)
                                                            newEdge2 = (numNodes, snd3 $ last first2Edges, 0.0 ∷ Double)
                                                            newEdge3 = (curNode, numNodes, 0.0 ∷ Double)
                                                            newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
                                                        in  newGraph
                                                    else -- root or simple network indegree node

                                                        if inNum == 0 || outNum > 2
                                                            then
                                                                let first2Edges = take 2 outEdgeList
                                                                    newNode = (numNodes, T.pack ("HTU" <> show numNodes))
                                                                    newEdge1 = (numNodes, snd3 $ head first2Edges, 0.0 ∷ Double)
                                                                    newEdge2 = (numNodes, snd3 $ last first2Edges, 0.0 ∷ Double)
                                                                    newEdge3 = (curNode, numNodes, 0.0 ∷ Double)
                                                                    newGraph = LG.insEdges [newEdge1, newEdge2, newEdge3] $ LG.delLEdges first2Edges $ LG.insNode newNode inGraph
                                                                in  newGraph
                                                            else -- check if indegree 0 is a leaf (ie index < root)

                                                                if outNum == 0
                                                                    then -- get root index

                                                                        let rootIndex = fst $ head $ LG.getRoots inGraph
                                                                        in  if curNode < rootIndex
                                                                                then inGraph
                                                                                else LG.delNode curNode inGraph
                                                                    else
                                                                        error
                                                                            ( "This can't happen in resolveNode in/out edge lists don't need to be resolved " <> show inOutPair <> "\n" <> LG.prettify inGraph
                                                                            )


-- )

{- | convertSimpleToDecoratedGraph takes a sinple graph and creates a Decorated graph
but with dummy info--basically just with the correct type and structures
-}
convertSimpleToDecoratedGraph ∷ SimpleGraph → DecoratedGraph
convertSimpleToDecoratedGraph inSimple =
    if LG.isEmpty inSimple
        then LG.empty
        else
            let simpleNodeList = LG.labNodes inSimple
                simpleEdgeList = LG.labEdges inSimple
                decNodeList = fmap simpleNodeToDecorated simpleNodeList
                decEdgeList = fmap simpleEdgeToDecorated simpleEdgeList
            in  LG.mkGraph decNodeList decEdgeList


-- | simpleNodeToDecorated takes a simple node and cretes a decoraetd node with info available
simpleNodeToDecorated ∷ LG.LNode T.Text → LG.LNode VertexInfo
simpleNodeToDecorated (indexNode, nameNode) =
    -- probbaly need to add other info--but cn;t get all of it (e..g. BV and costs)
    ( indexNode
    , emptyVertexInfo
        { index = indexNode
        , vertName = nameNode
        }
    )


-- | simpleEdgeToDecorated takes a Double edge label and returns EdgInfo
simpleEdgeToDecorated ∷ LG.LEdge Double → LG.LEdge EdgeInfo
simpleEdgeToDecorated (a, b, weightDouble) = (a, b, dummyEdge{minLength = weightDouble, maxLength = weightDouble, midRangeLength = weightDouble})

{-convertDecoratedToSimpleGraphBranchLength convertgs decorated to simple graph but allows
    specification of edge/branch weight/length ("min", "max", "mid") -}
convertDecoratedToSimpleGraphBranchLength :: String -> DecoratedGraph → SimpleGraph
convertDecoratedToSimpleGraphBranchLength branchWeight inDec =
    if LG.isEmpty inDec then LG.empty
    else if branchWeight `notElem` ["min", "max", "mid"] then 
        errorWithoutStackTrace ("Edge/Branch weight/length not recognized ('min', 'max', 'mid') : " <> branchWeight)
    else 
        let edgeWeight = if branchWeight == "min" then minLength
                         else if branchWeight == "max" then maxLength
                         else midRangeLength

            decNodeList = LG.labNodes inDec
            newNodeLabels = fmap (vertName . snd) decNodeList
            simpleNodes = zip (fmap fst decNodeList) newNodeLabels
            labEdgeList = LG.labEdges inDec
            edgeWeightList = fmap (edgeWeight . thd3) labEdgeList
            simpleEdgeList = fmap convertToSimpleEdge (zip edgeWeightList labEdgeList)
            in  LG.mkGraph simpleNodes simpleEdgeList


-- | convertDecoratedToSimpleGraph takes a decorated graph and returns the simple graph equivalent
convertDecoratedToSimpleGraph ∷ DecoratedGraph → SimpleGraph
convertDecoratedToSimpleGraph inDec = convertDecoratedToSimpleGraphBranchLength "min" inDec
    
-- | convertToSimpleEdge takes a lables edge and relabels with input
convertToSimpleEdge ∷ (VertexCost, LG.LEdge EdgeInfo) → LG.LEdge Double
convertToSimpleEdge (edgeWeight, (a, b, _)) = (a, b, edgeWeight)


{- | graphCostFromNodes takes a Decorated graph and returns its cost by summing up the local costs
 of its nodes
-}
graphCostFromNodes ∷ DecoratedGraph → Double
graphCostFromNodes inGraph =
    if LG.isEmpty inGraph
        then 0.0
        else sum $ fmap (vertexCost . snd) (LG.labNodes inGraph)

-- | dichotomizeRoot takes greaph and dichotimizes not dichotomous roots in graph
dichotomizeRoot ∷ Int → SimpleGraph → SimpleGraph
dichotomizeRoot lOutgroupIndex inGraph =
    if LG.isEmpty inGraph
        then LG.empty
        else
            let rootList = LG.getRoots inGraph
                currentRoot = fst $ head rootList
                rootEdgeList = LG.out inGraph currentRoot
            in  -- not a tree error
                if length rootList /= 1
                    then error ("Graph input to dichotomizeRoot is not a tree--not single root:" <> show rootList)
                    else -- nothing to do

                        if length rootEdgeList < 3
                            then inGraph
                            else
                                let numVertices = length $ LG.nodes inGraph
                                    newNode = (numVertices, T.pack $ show numVertices)
                                    edgesToDelete = filter ((/= lOutgroupIndex) . snd3) rootEdgeList
                                    newEdgeDestinations = fmap snd3 edgesToDelete
                                    newEdgeStarts = replicate (length newEdgeDestinations) numVertices
                                    newEdgeLabels = replicate (length newEdgeDestinations) 0.0
                                    -- nub for case where root edge in "wrong" direction
                                    -- doesn't filter edges to delete properly
                                    newEdgesNewNode = L.nub $ zip3 newEdgeStarts newEdgeDestinations newEdgeLabels
                                    newRootEdge = (currentRoot, numVertices, 0.0)
                                in  LG.delLEdges edgesToDelete $ LG.insEdges (newRootEdge : newEdgesNewNode) $ LG.insNode newNode inGraph


-- | showBlockGraphs takes a vector of vector of DecoratedGraphs and converte and prettifies outputting a String
showDecGraphs ∷ V.Vector (V.Vector DecoratedGraph) → String
showDecGraphs inDecVV =
    if V.null inDecVV
        then []
        else concatMap concat (V.toList $ fmap ((V.toList . fmap LG.prettify) . fmap convertDecoratedToSimpleGraph) inDecVV)


{- |
A wrapper around selectGraphsFull for ReducedPhylogeneticGraph.
-}
selectGraphs ∷ SelectGraphType → LG.Node -> Int → Double → [ReducedPhylogeneticGraph] → PhyG [ReducedPhylogeneticGraph]
selectGraphs selectType outgroupIndex numberToKeep threshold inGraphList =
    let convertReduced2GenPhyloGraph (a, b, c, d, f) = (a, b, c, d, mempty, f)

        convertGenPhyloGraph2Reduced (a, b, c, d, _, f) = (a, b, c, d, f)

        fullPhyloGraphList = convertReduced2GenPhyloGraph <$> inGraphList

        newFullGraphs = selectGraphsFull selectType outgroupIndex numberToKeep threshold fullPhyloGraphList
    in  fmap convertGenPhyloGraph2Reduced <$> newFullGraphs


-- Basically a Phylogenetic Graph with abstract graph types--can't seem to get a type with this to compile
-- (SimpleGraph, VertexCost, LG.Gr a b, V.Vector [LG.Gr a b], V.Vector (V.Vector (LG.Gr a b)), V.Vector (V.Vector CharInfo))

{- |
A wrapper around selectPhylogeneticGraph with a better interface.
-}
selectGraphsFull
    ∷ SelectGraphType
    -> LG.Node
    → Int
    → Double
    → [GenPhyloGraph a b]
    → PhyG [GenPhyloGraph a b]
selectGraphsFull selectType outgroupIndex numberToKeep threshold = \case
    [] → pure []
    inGraphList →
        let stringArgs
                | threshold > 0.0 = ("threshold", show threshold)
                | otherwise = case selectType of
                    Best → ("best", "")
                    Unique → ("unique", "")
                    AtRandom → ("atrandom", "")
                    All → ("all", "")
                    _ → ("best", "")
        in  take numberToKeep <$> selectPhylogeneticGraph [stringArgs] outgroupIndex [] inGraphList


{- |
A wrapper for ReducedPhylogeneticGraph.

/Note: Inefficient in conversions./
-}
selectPhylogeneticGraphReduced ∷ [Argument] -> LG.Node → [ReducedPhylogeneticGraph] → PhyG [ReducedPhylogeneticGraph]
selectPhylogeneticGraphReduced inArgs outgroupIndex curGraphs =
    let phylographList = fmap convertReduced2PhylogeneticGraph curGraphs
    in  do
            selectedPhylographs ← selectPhylogeneticGraph inArgs outgroupIndex [] phylographList
            pure $ convertPhylogeneticGraph2Reduced <$> selectedPhylographs


{- | selectPhylogeneticGraph takes  a series of arguments and an input list ot PhylogeneticGraphs
and returns or filters that list based on options.
uses selectListCostPairs in GeneralUtilities

Operates on collapsed graphs
-}
selectPhylogeneticGraph
    ∷ [Argument]
    -> LG.Node
    → [String]
    → [GenPhyloGraph a b]
    → PhyG [GenPhyloGraph a b]
selectPhylogeneticGraph inArgs outgroupIndex _ [] = pure []
selectPhylogeneticGraph inArgs outgroupIndex _ curGraphs =
    let fstArgList = fmap (fmap C.toLower . fst) inArgs
        sndArgList = fmap (fmap C.toLower . snd) inArgs
        lcArgList = zip fstArgList sndArgList
        checkCommandList = checkCommandArgs "select" fstArgList V.selectArgList
    in  -- check for valid command options
        if not checkCommandList
            then errorWithoutStackTrace ("Unrecognized command in 'select': " <> show inArgs)
            else
                if length inArgs > 1
                    then errorWithoutStackTrace ("Can only have a single select type per command: " <> show inArgs)
                    else
                        let doBest = keyExists "best" lcArgList
                            doAll = keyExists "all" lcArgList
                            doRandom = keyExists "atrandom" lcArgList
                            doUnique = keyExists "unique" lcArgList
                            doThreshold = keyExists "threshold" lcArgList

                            thresholdList = filter ((== "threshold") . fst) lcArgList
                            thresholdParam
                                | length thresholdList > 1 =
                                    errorWithoutStackTrace ("Multiple 'threshold' number specifications in select command--can have only one: " <> show inArgs)
                                | null thresholdList = Just 0.1
                                | otherwise = readMaybe (snd $ head thresholdList) ∷ Maybe Double

                            nonThresholdList = filter ((/= "threshold") . fst) lcArgList
                            numberToKeep
                                | length nonThresholdList > 1 =
                                    errorWithoutStackTrace
                                        ("Multiple 'best/unique/atRandom' number specifications in select command--can have only one: " <> show inArgs)
                                | null nonThresholdList = Just (maxBound ∷ Int)
                                | null (snd $ head nonThresholdList) = Just (maxBound ∷ Int)
                                | otherwise = readMaybe (snd $ head nonThresholdList) ∷ Maybe Int
                        in  case numberToKeep of
                                Nothing →
                                    failWithPhase Computing $
                                        "Keep specification not an integer in select: " <> show (snd $ head nonThresholdList) <> show lcArgList
                                Just keep → case thresholdParam of
                                    Nothing → failWithPhase Computing $ "Threshold specification not a float in select: " <> show (snd $ head thresholdList)
                                    Just _ | doAll → pure curGraphs
                                    Just threshold →
                                        let -- minimum graph cost
                                            minGraphCost = minimum $ fmap snd6 curGraphs

                                            -- collapse zero-length branchs for unique
                                            curGraphsCollapsed = fmap U.collapseGraph curGraphs

                                            -- keep only unique graphs based on non-zero edges--in sorted by cost
                                            uniqueGraphList = L.sortOn snd6 $ nubGraph outgroupIndex [] (zip curGraphs curGraphsCollapsed) -- curGraphs --  True curGraphs -- getBVUniqPhylogeneticGraph True curGraphs -- getTopoUniqPhylogeneticGraph True curGraphs

                                            -- this to avaoid alot of unncesesary graph comparisons for 'best' graphs
                                            bestCostGraphs = filter ((== minGraphCost) . snd6) curGraphs
                                            uniqueBestGraphs = nubGraph outgroupIndex [] (zip bestCostGraphs (fmap U.collapseGraph bestCostGraphs))

                                            result
                                                | doUnique = pure $ take keep uniqueGraphList
                                                | doThreshold =
                                                    let baseGraphValue = snd6 $ head uniqueGraphList
                                                        thresholdLimit = thresholdValue * baseGraphValue
                                                        thresholdValue
                                                            | threshold < 0.0 = 1
                                                            | threshold > 1.0 = threshold
                                                            | otherwise = 1.0 + threshold

                                                        thresholdGraphList = filter ((<= thresholdLimit) . snd6) uniqueGraphList
                                                    in  pure thresholdGraphList
                                                | doBest = pure $ take keep uniqueBestGraphs
                                                | doRandom = shuffleList $ take keep curGraphs
                                                | otherwise = pure uniqueBestGraphs -- default is all best and unique
                                        in  result


{- | getUniqueGraphs takes each pair of non-zero edges and conpares them--if equal not added to list
maybe chnge to nub LG.pretify graphList?
-}
getUniqueGraphs
    ∷ Bool
    → [GenPhyloGraph a b]
    → [GenPhyloGraph a b]
getUniqueGraphs removeZeroEdges inGraphList =
    if null inGraphList
        then []
        else
            let inGraphEdgeList =
                    if removeZeroEdges
                        then fmap ((filter ((> 0.0) . minLength . thd3) . LG.labEdges) . thd6) inGraphList
                        else fmap (LG.labEdges . thd6) inGraphList
            in  getUniqueGraphs' (zip inGraphEdgeList inGraphList) []


{- | getUniqueGraphs Using fgl ==
basically a nub
need to add a collapse function for compare as well
takes pairs of (noCollapsed, collapsed) phylogenetic graphs,
make strings based on collapsed and returns not collpased

getUniqueGraphs''
    ∷ LG.Node
    -> [(GenPhyloGraph a b, GenPhyloGraph a b)]
    → [GenPhyloGraph a b]
getUniqueGraphs'' = nubGraph outgroupIndex []
-}


{- | isNovelGraph  checks if a graph is in list of existing graphs
uses colapsed representation
-}
isNovelGraph
    ∷ [GenPhyloGraph a b]
    → GenPhyloGraph a b
    → Bool
isNovelGraph graphList testGraph =
    null graphList
        || ( let collapsedInGraph = (LG.prettyIndices . fst6 . U.collapseGraph) testGraph
                 collapseGraphList = fmap (LG.prettyIndices . fst6 . U.collapseGraph) graphList
                 matchList = filter (== collapsedInGraph) collapseGraphList
             in  -- trace ("IsNovel: " <> (show $ null matchList))
                 null matchList
           )


{- | keeps and returns unique graphs based on Eq of Topological Simple Graph
String prettyIndices w/0 HTU names and branch lengths
arbitrarily rooted on 0 for oonsistency
reversed to keep original order in case sorted on length

reroots graph for comparison if needed--but not for networks
-}
nubGraph
    ∷ LG.Node
    -> [(GenPhyloGraph a b, GenPhyloGraph a b, String)]
    → [(GenPhyloGraph a b, GenPhyloGraph a b)]
    → [GenPhyloGraph a b]
nubGraph outgroupIndex curList inList =
    if null inList
        then reverse $ fmap fst3 curList
        else
            let (firstGraphNC, firstGraphC) = head inList

                -- nub on newick topology only--should be collapsed already
                -- reroot for comparison since can move root in search

                -- check curent graph rpoot to see if needs to be rerooted
                inRoot = (fst $ head $ LG.getRoots $ fst6 firstGraphNC)
                inRootChildren = LG.descendants (fst6 firstGraphNC) inRoot
                isOutGroupRooted = outgroupIndex `elem` inRootChildren

                costList = [snd6 firstGraphNC]
                rerootGraph = LG.rerootTree outgroupIndex $ fst6 firstGraphNC
                graphRoot = head $ LG.parents rerootGraph outgroupIndex

                              -- is not a tree so dangerous to reroot
                firstString = --if (not $ LG.isTree (fst6 firstGraphNC)) then 
                              if isOutGroupRooted || (not $ LG.isTree (fst6 firstGraphNC)) then 
                                     makeNewickList False False False (fst $ head $ LG.getRoots $ fst6 firstGraphNC) [fst6 firstGraphNC] [snd6 firstGraphNC]

                              -- tried to reroot but not properly rooted--eg tried to reroot on a network edge
                              else if (null (LG.parents rerootGraph outgroupIndex)) then 
                                    makeNewickList False False False (fst $ head $ LG.getRoots $ fst6 firstGraphNC) [fst6 firstGraphNC] [snd6 firstGraphNC]

                              -- not in correct root but ok to be rerooted--ie tree
                              else
                                    makeNewickList False False False graphRoot [rerootGraph] costList
                
                -- nub on prettty string
                -- firstString = LG.prettyIndices $ thd6 firstGraphNC
                isMatch = filter (== firstString) (fmap thd3 curList)

            in  
                --trace (firstString <> "\n" <> (concatMap thd3 curList)) $ 
                if null curList
                    then nubGraph outgroupIndex [(firstGraphNC, firstGraphC, firstString)] (tail inList)
                    else
                        if null isMatch
                            then nubGraph outgroupIndex ((firstGraphNC, firstGraphC, firstString) : curList) (tail inList)
                            else nubGraph outgroupIndex curList (tail inList) 


-- )

-- | getUniqueGraphs takes each pair of non-zero edges and compares them--if equal not added to list
getUniqueGraphs'
    ∷ [([LG.LEdge EdgeInfo], GenPhyloGraph a b)]
    → [([LG.LEdge EdgeInfo], GenPhyloGraph a b)]
    → [GenPhyloGraph a b]
getUniqueGraphs' inGraphPairList currentUniquePairs =
    if null inGraphPairList
        then fmap snd currentUniquePairs
        else
            let firstPair@(firstEdges, _) = head inGraphPairList
            in  if null currentUniquePairs
                    then getUniqueGraphs' (tail inGraphPairList) [firstPair]
                    else
                        let equalList = filter id $ fmap ((== firstEdges) . fst) currentUniquePairs
                        in  if null equalList
                                then getUniqueGraphs' (tail inGraphPairList) (firstPair : currentUniquePairs)
                                else getUniqueGraphs' (tail inGraphPairList) currentUniquePairs


-- | getNodeType returns node type for Node
getNodeType ∷ (Show a, Show b) ⇒ LG.Gr a b → LG.Node → NodeType
getNodeType inGraph inNode
    | not $ LG.gelem inNode inGraph = error ("Node " <> show inNode <> " not in graph\n" <> GFU.showGraph inGraph)
    | LG.isLeaf inGraph inNode = LeafNode
    | LG.isTreeNode inGraph inNode = TreeNode
    | LG.isNetworkNode inGraph inNode = NetworkNode
    | LG.isRoot inGraph inNode = RootNode
    | LG.isIn1Out1 inGraph inNode = In1Out1
    | otherwise = error ("Node type " <> show inNode <> " not Leaf, Tree, Network, or Root in graph\n" <> GFU.showGraph inGraph)


{- | copyIAFinalToPrelim takes a Decorated graph and copies
the IA final fields to preliminary IA states--this for IA only optimization
inswapping and other operations.  This is done because the "preliminary" IA states
are only known after full post/pre traversals
-}
copyIAFinalToPrelim ∷ DecoratedGraph → DecoratedGraph
copyIAFinalToPrelim inGraph =
    if LG.isEmpty inGraph
        then error "Empty input graph in copyIAFinalToPrelim"
        else
            let nodes = LG.labNodes inGraph
                edges = LG.labEdges inGraph
                newNodes = fmap makeIAPrelimFromFinal nodes
            in  LG.mkGraph newNodes edges


{- | makeIAPrelimFromFinal updates the label of a node for IA states
setting preliminary to final
-}
makeIAPrelimFromFinal ∷ LG.LNode VertexInfo → LG.LNode VertexInfo
makeIAPrelimFromFinal (inIndex, label) =
    let labData = vertData label
        newLabData = fmap (fmap f) labData
    in  (inIndex, label{vertData = newLabData})
    where
        f c
            | GV.null (slimIAFinal c) && GV.null (wideIAFinal c) && GV.null (hugeIAFinal c) = c
            | not $ GV.null $ slimIAFinal c = c{slimIAPrelim = M.makeDynamicCharacterFromSingleVector $ slimIAFinal c}
            | not $ GV.null $ wideIAFinal c = c{wideIAPrelim = M.makeDynamicCharacterFromSingleVector $ wideIAFinal c}
            | otherwise = c{hugeIAPrelim = M.makeDynamicCharacterFromSingleVector $ hugeIAFinal c}


{- | copyIAPrelimToFinal takes a Decorated graph and copies
the IA prelim fields to final IA states--this for IA only optimization
inswapping and other operations.  THis is fdone for root and leaf vertices
-}
copyIAPrelimToFinal ∷ DecoratedGraph → DecoratedGraph
copyIAPrelimToFinal inGraph =
    if LG.isEmpty inGraph
        then error "Empty input graph in copyIAFinalToPrelim"
        else
            let nodes = LG.labNodes inGraph
                edges = LG.labEdges inGraph
                newNodes = fmap makeIAFinalFromPrelim nodes
            in  LG.mkGraph newNodes edges


{- | makeIAFinalFomPrelim updates the label of a node for IA states
setting final to preliminary
-}
makeIAFinalFromPrelim ∷ LG.LNode VertexInfo → LG.LNode VertexInfo
makeIAFinalFromPrelim (inIndex, label) =
    let labData = vertData label
        newLabData = fmap (fmap f) labData
    in  (inIndex, label{vertData = newLabData})
    where
        f c =
            let newSlimIAFinal = extractMediansGapped $ slimIAPrelim c
                newWideIAFinal = extractMediansGapped $ wideIAPrelim c
                newHugeIAFinal = extractMediansGapped $ hugeIAPrelim c
            in  if GV.null (snd3 $ slimIAPrelim c) && GV.null (snd3 $ wideIAPrelim c) && GV.null (snd3 $ hugeIAPrelim c)
                    then c
                    else
                        if not $ GV.null $ snd3 $ slimIAPrelim c
                            then c{slimIAFinal = newSlimIAFinal}
                            else
                                if not $ GV.null $ snd3 $ wideIAPrelim c
                                    then c{wideIAFinal = newWideIAFinal}
                                    else c{hugeIAFinal = newHugeIAFinal}


{- | getTopoUniqPhylogeneticGraph takes a list of phylogenetic graphs and returns
list of topologically unique graphs--operatres on simple graph field
noZeroEdges flag passed to remove zero weight edges
-}
getTopoUniqPhylogeneticGraph ∷ Bool → [PhylogeneticGraph] → [PhylogeneticGraph]
getTopoUniqPhylogeneticGraph nonZeroEdges inPhyloGraphList =
    if null inPhyloGraphList
        then []
        else
            let uniqueBoolList = createUniqueBoolList nonZeroEdges (fmap fst6 inPhyloGraphList) []
                boolPair = zip inPhyloGraphList uniqueBoolList
            in  (fst <$> filter snd boolPair)


-- | createUniqueBoolList creates a list of Bool if graphs are unique--first occurrence is True, others False
createUniqueBoolList ∷ Bool → [SimpleGraph] → [(SimpleGraph, Bool)] → [Bool]
createUniqueBoolList nonZeroEdges inGraphList boolAccum =
    if null inGraphList
        then reverse $ fmap snd boolAccum
        else
            let firstGraph = head inGraphList
            in  if null boolAccum
                    then createUniqueBoolList nonZeroEdges (tail inGraphList) ((firstGraph, True) : boolAccum)
                    else
                        let checkList = filter id $ fmap (topologicalEqual nonZeroEdges firstGraph . fst) boolAccum
                        in  if null checkList
                                then createUniqueBoolList nonZeroEdges (tail inGraphList) ((firstGraph, True) : boolAccum)
                                else createUniqueBoolList nonZeroEdges (tail inGraphList) ((firstGraph, False) : boolAccum)


{- | topologicalEqual takes two simple graphs and returns True if graphs have same nodes and edges
option to exclude zero weight edges
-}
topologicalEqual ∷ Bool → SimpleGraph → SimpleGraph → Bool
topologicalEqual nonZeroEdges g1 g2
    | LG.isEmpty g1 && LG.isEmpty g2 = True
    | LG.isEmpty g1 || LG.isEmpty g2 = False
    | otherwise =
        let nodesG1 = LG.labNodes g1
            nodesG2 = LG.labNodes g2
            edgesG1 =
                if nonZeroEdges
                    then fmap LG.toEdge $ filter ((> 0) . thd3) $ LG.labEdges g1
                    else LG.edges g1
            edgesG2 =
                if nonZeroEdges
                    then fmap LG.toEdge $ filter ((> 0) . thd3) $ LG.labEdges g2
                    else LG.edges g2
        in  nodesG1 == nodesG2 && edgesG1 == edgesG2


{- | getEdgeMinLengthToNode takes a labelled node and returns the min length of
the edge leading to the node
-}
getEdgeMinLengthToNode ∷ [LG.LEdge EdgeInfo] → LG.LNode a → Double
getEdgeMinLengthToNode edgeList (node, _) =
    let foundEdge = L.find ((== node) . snd3) edgeList
    in  -- root node will be nor be in in edge set and need so set > 0
        if isNothing foundEdge
            then 1.0 --  error ("Edge not found in getEdgeMinLengthToNode: node " <> (show node) <> " edge list " <> (show edgeList))
            else minLength $ thd3 $ fromJust foundEdge


{- | getBVUniqPhylogeneticGraph takes a list of phylogenetic graphs and returns
list of topologically unique graphs based on their node bitvector assignments
operatres on Decorated graph field
noZeroEdges flag passed to remove zero weight edges
-}
getBVUniqPhylogeneticGraph ∷ Bool → [PhylogeneticGraph] → [PhylogeneticGraph]
getBVUniqPhylogeneticGraph nonZeroEdges inPhyloGraphList =
    if null inPhyloGraphList
        then []
        else
            let bvGraphList = fmap (getBVNodeList nonZeroEdges . thd6) inPhyloGraphList
                uniqueBoolList = createBVUniqueBoolList bvGraphList []
                boolPair = zip inPhyloGraphList uniqueBoolList
            in  (fst <$> filter snd boolPair)


{- | getBVNodeList takes a DecoratedGraph and returns sorted list (by BV) of nodes
removes node with zero edge weight to them if specified
-}
getBVNodeList ∷ Bool → DecoratedGraph → [BV.BitVector]
getBVNodeList nonZeroEdges inGraph =
    if LG.isEmpty inGraph
        then []
        else
            let nodeList = LG.labNodes inGraph
                edgeList = LG.labEdges inGraph
                minLengthList = fmap (getEdgeMinLengthToNode edgeList) nodeList
                nodePairList = filter ((> 0) . snd) $ zip nodeList minLengthList
                bvNodeList =
                    if nonZeroEdges
                        then L.sort $ fmap ((bvLabel . snd) . fst) nodePairList
                        else L.sort $ fmap (bvLabel . snd) nodeList
            in  bvNodeList


{- | createBVUniqueBoolList creates a list of Bool if graphs are unique by bitvecector node list
first occurrence is True, others False
assumes edges filterd b=y lenght already
-}
createBVUniqueBoolList ∷ [[BV.BitVector]] → [([BV.BitVector], Bool)] → [Bool]
createBVUniqueBoolList inBVGraphListList boolAccum =
    if null inBVGraphListList
        then reverse $ fmap snd boolAccum
        else
            let firstGraphList = head inBVGraphListList
            in  if null boolAccum
                    then createBVUniqueBoolList (tail inBVGraphListList) ((firstGraphList, True) : boolAccum)
                    else
                        let checkList = filter id $ fmap ((== firstGraphList) . fst) boolAccum
                        in  if null checkList
                                then createBVUniqueBoolList (tail inBVGraphListList) ((firstGraphList, True) : boolAccum)
                                else createBVUniqueBoolList (tail inBVGraphListList) ((firstGraphList, False) : boolAccum)


-- | makeDummyLabEdge takes an unlabelled edge and adds a dummy label
makeDummyLabEdge ∷ EdgeInfo → LG.Edge → LG.LEdge EdgeInfo
makeDummyLabEdge edgeLab (u, v) = (u, v, edgeLab)


{- | netNodeAncestorViolation checks whether one of the edge into a netowrk node (in 2)
is cinnected to an ancestor (via the other parent) of the node
this is a form of time violation since the parents of a network node must be
at least possibly coeval
this uses the bit vector label of nodes.  If the other child of either parent node
of a network node has non-zero intersection between the BV label of the network node
and that other child of parent then they conecting edge is from and ancestral node hence a time violation
O(n) n netork nodes in Graph, but checks all nodes to see if network
-}
hasNetNodeAncestorViolation ∷ LG.Gr VertexInfo b → Bool
hasNetNodeAncestorViolation inGraph =
    if LG.isEmpty inGraph
        then error "Empty graph in hasNetNodeAncestorViolation"
        else
            let (_, _, _, netWorkNodeList) = LG.splitVertexList inGraph
                hasAncViolationList = filter id $ fmap (nodeAncViolation inGraph) netWorkNodeList
            in  -- trace ("HNV: " <> (show $ (not . null) hasAncViolationList))
                (not . null) hasAncViolationList


{- | nodeAncViolation checks a single node fo ancestrpo connection--he ceviolation
should be O(1).  Return True if violation
-}
nodeAncViolation ∷ LG.Gr VertexInfo b → LG.LNode VertexInfo → Bool
nodeAncViolation inGraph inNode =
    let parentList = LG.labParents inGraph (fst inNode)
    in  if length parentList /= 2
            then error ("Parent number should be 2: " <> show (fst inNode) <> " <- " <> show (fmap fst parentList))
            else
                let sisterNodes = concatMap (LG.sisterLabNodes inGraph) parentList
                    sisterBVData = fmap (bvLabel . snd) sisterNodes
                    inNodeBVData = bvLabel $ snd inNode
                    sisterBVIntersections = fmap (.&. inNodeBVData) sisterBVData
                    isAncInNode = filter (== inNodeBVData) sisterBVIntersections
                in  (not . null) isAncInNode


{- | selectGraphStochastic takes a list of graphs and returns a list of graphs chosen at Random
using an exponential distribution based on graph cost difference divided by an input factor
if factor is 0 then stringth graphs cost
mprob acceptance = -exp [(cost - minCost)/ factor]
returns n graphs by random criterion without replacment
-}
selectGraphStochastic ∷ Int → Double → [PhylogeneticGraph] → PhyG [PhylogeneticGraph]
selectGraphStochastic number factor inGraphList
    | null inGraphList = pure inGraphList
    | number >= length inGraphList = pure inGraphList
    | otherwise =
        let getProb ∷ Double → Double
            getProb a = exp ((-1) * factor / a)

            minCost = minimum $ snd6 <$> inGraphList

            selectGraphByWeighting ∷ (Double → Double → Bool) → PhylogeneticGraph → PhyG Bool
            selectGraphByWeighting cmp v =
                let extractKey = getProb . (-) minCost . snd6
                in  getRandomR (0, 1) <&> cmp (extractKey v)
        in  do
                -- so no more than specified
                selectedGraphs ← take number <$> filterM (selectGraphByWeighting (<)) inGraphList

                -- takes some random remainder to fill out length of list
                suffixPad ← case number - length selectedGraphs of
                    0 → pure []
                    numLucky → filterM (selectGraphByWeighting (>=)) inGraphList >>= (fmap (take numLucky) . shuffleList)

                pure $ selectedGraphs <> suffixPad

{- | getDisplayTreeCostList returns a list of the "block" costs of display trees
in a pair with any graph 'penalty' cost
-}
getDisplayTreeCostList ∷ PhylogeneticGraph → ([VertexCost], VertexCost)
getDisplayTreeCostList inGraph =
    if LG.isEmpty $ thd6 inGraph
        then ([], 0.0)
        else
            let rootIndex = fst $ head $ LG.getRoots $ fst6 inGraph
                displayTreeCharVect = fft6 inGraph
                displayTreeCostVect = fmap (getBlockCost rootIndex) displayTreeCharVect
                nonGraphCost = V.sum displayTreeCostVect
            in  (V.toList displayTreeCostVect, snd6 inGraph - nonGraphCost)


-- | getBlockCost returns the cost, summed over characters, of a character block
getBlockCost ∷ LG.Node → V.Vector DecoratedGraph → VertexCost
getBlockCost rootIndex charGraphVect =
    if V.null charGraphVect
        then 0.0
        else V.sum $ fmap (getCharacterCost rootIndex) charGraphVect


-- | getCharacterCost returns charcter cost as root of character tree
getCharacterCost ∷ LG.Node → DecoratedGraph → VertexCost
getCharacterCost rootIndex inGraph =
    if LG.isEmpty inGraph
        then 0.0
        else
            let rootLabel = LG.lab inGraph rootIndex
            in  maybe (error ("Root index without label: " <> show rootIndex)) subGraphCost rootLabel


{- | makeLeafGraph takes input data and creates a 'graph' of leaves with Vertex informnation
but with zero edges.  This 'graph' can be reused as a starting structure for graph construction
to avoid remaking of leaf vertices
-}
makeLeafGraph ∷ ProcessedData → DecoratedGraph
makeLeafGraph (nameVect, bvNameVect, blocDataVect) =
    if V.null nameVect
        then error "Empty ProcessedData in makeLeafGraph"
        else
            let leafVertexList = V.toList $ V.map (makeLeafVertex nameVect bvNameVect blocDataVect) (V.fromList [0 .. V.length nameVect - 1])
            in  LG.mkGraph leafVertexList []


{- | makeSimpleLeafGraph takes input data and creates a 'graph' of leaves with Vertex informnation
but with zero edges.  This 'graph' can be reused as a starting structure for graph construction
to avoid remaking of leaf vertices
-}
makeSimpleLeafGraph ∷ ProcessedData → SimpleGraph
makeSimpleLeafGraph (nameVect, _, _) =
    if V.null nameVect
        then error "Empty ProcessedData in makeSimpleLeafGraph"
        else
            let leafVertexList = V.toList $ V.map (makeSimpleLeafVertex nameVect) (V.fromList [0 .. V.length nameVect - 1])
            in  LG.mkGraph leafVertexList []
    where
        makeSimpleLeafVertex ∷ ∀ {b}. V.Vector b → Int → (Int, b)
        makeSimpleLeafVertex a b = (b, a V.! b)


-- | makeLeafVertex makes a single unconnected vertex for a leaf
makeLeafVertex ∷ V.Vector NameText → V.Vector NameBV → V.Vector BlockData → Int → LG.LNode VertexInfo
makeLeafVertex nameVect bvNameVect inData localIndex =
    -- trace ("Making leaf " <> (show localIndex) <> " Data " <> (show $ length inData) <> " " <> (show $ fmap length $ fmap snd3 inData)) (
    let centralData = V.map snd3 inData
        thisData = V.map (V.! localIndex) centralData
        newVertex =
            VertexInfo
                { index = localIndex
                , bvLabel = bvNameVect V.! localIndex
                , parents = V.empty
                , children = V.empty
                , nodeType = LeafNode
                , vertName = nameVect V.! localIndex
                , vertData = thisData
                , vertexResolutionData = mempty
                , vertexCost = 0.0
                , subGraphCost = 0.0
                }
    in  -- trace (show (length thisData) <> (show $ fmap length thisData))
        (localIndex, newVertex)


-- )

{- | remakePhylogeneticGraph remakes (rebuilds from scratch) phylogenetic graph
fst, thd. 4th and 5th fields
-}
remakePhylogeneticGraph ∷ PhylogeneticGraph → PhylogeneticGraph
remakePhylogeneticGraph inGraph@(a, b, c, d, e, f) =
    if inGraph == emptyPhylogeneticGraph
        then inGraph
        else
            let a' = LG.remakeGraph a
                c' = LG.remakeGraph c
                d' = fmap (fmap LG.remakeGraph) d
                e' = fmap (fmap LG.remakeGraph) e
            in  (a', b, c', d', e', f)


keyExists ∷ ∀ f k v. (Eq k, Foldable f) ⇒ k → f (k, v) → Bool
keyExists key = any ((== key) . fst)
