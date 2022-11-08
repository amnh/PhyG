{- |
Module      :  Traversals.hs
Description :  Module specifying graph traversal functions for PhyGraph
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

module GraphOptimization.Traversals ( postOrderTreeTraversal
                                    , multiTraverseFullyLabelTree
                                    , multiTraverseFullyLabelGraph
                                    , multiTraverseFullyLabelGraph'
                                    , multiTraverseFullyLabelSoftWired
                                    , multiTraverseFullyLabelHardWired
                                    , checkUnusedEdgesPruneInfty
                                    , makeLeafGraph
                                    , makeSimpleLeafGraph
                                    , postDecorateTree
                                    , postDecorateTree'
                                    , updatePhylogeneticGraphCost
                                    , updateGraphCostsComplexities
                                    , getW15NetPenalty
                                    , getW23NetPenalty
                                    , getW15RootCost
                                    , generalizedGraphPostOrderTraversal
                                    -- added to quiet warnings
                                    , setBetterGraphAssignment
                                    , makeBetterBlock
                                    , chooseBetterCharacter
                                    , minimalReRootPhyloGraph
                                    ) where


import qualified Data.List                            as L
import qualified Data.Vector                          as V
import           GeneralUtilities
import qualified GraphFormatUtilities                 as GFU
import           Types.Types
import qualified Utilities.LocalGraph                 as LG
import qualified GraphOptimization.PostOrderFunctions as PO
import qualified GraphOptimization.PreOrderFunctions  as PRE
import qualified Graphs.GraphOperations               as GO
import           Data.Bits
import           Data.Maybe
import qualified Data.Text.Lazy                       as T
import           Debug.Trace
import           Utilities.Utilities                  as U
import qualified GraphOptimization.PostOrderSoftWiredFunctions as POSW
import           Control.Parallel.Strategies
import qualified ParallelUtilities                    as PU

-- | multiTraverseFullyLabelGraph is a wrapper around multi-traversal functions for Tree,
-- Soft-wired network graph, and Hard-wired network graph
-- can either find root, be given root, or start somwhere else (startVertex) do optimize only a component of a forest
multiTraverseFullyLabelGraph :: GlobalSettings -> ProcessedData -> Bool -> Bool -> Maybe Int -> SimpleGraph -> PhylogeneticGraph
multiTraverseFullyLabelGraph inGS inData pruneEdges warnPruneEdges startVertex inGraph
  | LG.isEmpty inGraph = emptyPhylogeneticGraph
  | graphType inGS == Tree =
    -- test for Tree
    let (_, _, _, networkVertexList) = LG.splitVertexList inGraph
    in
    if null networkVertexList then
        let leafGraph = makeLeafGraph inData
        in multiTraverseFullyLabelTree inGS inData leafGraph startVertex inGraph
    else errorWithoutStackTrace ("Input graph is not a tree/forest, but graph type has been specified (perhaps by default) as Tree. Modify input graph or use 'set()' command to specify network type\n\tNetwork vertices: " ++ (show $ fmap fst networkVertexList) ++ "\n" ++ (LG.prettify inGraph))
  | graphType inGS == SoftWired =
    let leafGraph = POSW.makeLeafGraphSoftWired inData
    in multiTraverseFullyLabelSoftWired  inGS inData pruneEdges warnPruneEdges leafGraph startVertex inGraph
  | graphType inGS == HardWired =
    let leafGraph = makeLeafGraph inData
    in multiTraverseFullyLabelHardWired inGS inData leafGraph startVertex inGraph
  | otherwise = errorWithoutStackTrace ("Unknown graph type specified: " ++ show (graphType inGS))


-- |multiTraverseFullyLabelGraph' maps to multiTraverseFullyLabelGraph with differnet order of arguments used by report IA and tnt output
multiTraverseFullyLabelGraph' :: GlobalSettings -> Bool -> Bool -> Maybe Int -> ProcessedData -> SimpleGraph -> PhylogeneticGraph
multiTraverseFullyLabelGraph' inGS pruneEdges warnPruneEdges startVertex inData inGraph = multiTraverseFullyLabelGraph inGS inData pruneEdges warnPruneEdges startVertex inGraph


multiTraverseFullyLabelHardWired :: GlobalSettings -> ProcessedData -> DecoratedGraph -> Maybe Int -> SimpleGraph -> PhylogeneticGraph
multiTraverseFullyLabelHardWired inGS inData leafGraph startVertex inSimpleGraph = multiTraverseFullyLabelTree inGS inData leafGraph startVertex inSimpleGraph

-- | multiTraverseFullyLabelSoftWired fully labels a softwired network component forest
-- including traversal rootings-- does not reroot on network edges
-- allows indegree=outdegree=1 vertices
-- pruneEdges and warnPruneEdges specify if unused edges (ie not in diuaplytrees) are pruned from
-- canonical tree or if an infinity cost is returned and if a trace warning is thrown if so.
-- in general--input trees should use "pruneEdges" during search--not
-- can either find root, be given root, or start somwhere else (startVertex) do optimize only a component of a forest
-- first Bool for calcualting breanch edger weights
multiTraverseFullyLabelSoftWired :: GlobalSettings -> ProcessedData -> Bool -> Bool -> DecoratedGraph -> Maybe Int -> SimpleGraph -> PhylogeneticGraph
multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph startVertex inSimpleGraph =
    if LG.isEmpty inSimpleGraph then emptyPhylogeneticGraph
    else
        let sequenceChars = U.getNumberSequenceCharacters (thd3 inData)
            (postOrderGraph, localStartVertex) = generalizedGraphPostOrderTraversal inGS sequenceChars inData leafGraph False startVertex inSimpleGraph
            fullyOptimizedGraph = PRE.preOrderTreeTraversal inGS (finalAssignment inGS) False True (sequenceChars > 0) localStartVertex False postOrderGraph
        in
        --trace ("MTFLS:\n" ++ (show $ thd6 postOrderGraph))
        checkUnusedEdgesPruneInfty inGS inData pruneEdges warnPruneEdges leafGraph $ updatePhylogeneticGraphCost fullyOptimizedGraph (snd6 fullyOptimizedGraph)

-- | multiTraverseFullyLabelTree performs potorder on default root and other traversal foci, taking the minimum
-- traversal cost for all nonexact charcters--the initial rooting is used for exact characters
-- operates with Tree functions
-- need to add forest functionality--in principle just split into components and optimize them independently
-- but get into root index issues the way this is written now.
-- can either find root, be given root, or start somwhere else (startVertex) do optimize only a component of a forest
multiTraverseFullyLabelTree :: GlobalSettings -> ProcessedData -> DecoratedGraph -> Maybe Int -> SimpleGraph -> PhylogeneticGraph
multiTraverseFullyLabelTree inGS inData leafGraph startVertex inSimpleGraph =
    if LG.isEmpty inSimpleGraph then emptyPhylogeneticGraph
    else
        let sequenceChars = U.getNumberSequenceCharacters (thd3 inData)
            -- False for staticIA
            (postOrderGraph, localStartVertex) = generalizedGraphPostOrderTraversal inGS sequenceChars inData leafGraph False startVertex inSimpleGraph
        in
        PRE.preOrderTreeTraversal inGS (finalAssignment inGS) False True (sequenceChars > 0) localStartVertex False postOrderGraph


-- | generalizedGraphPostOrderTraversal performs the postorder pass
-- on a graph (tree, softWired, or hardWired) to determine the "preliminary" character states
-- include penalty factor cost but not root cost which may or may not be wanted depending on context
-- if full graph--yes, if a component yes or no.
-- hence returnde das pair
generalizedGraphPostOrderTraversal :: GlobalSettings -> Int -> ProcessedData -> DecoratedGraph -> Bool -> Maybe Int -> SimpleGraph -> (PhylogeneticGraph, Int)
generalizedGraphPostOrderTraversal inGS sequenceChars inData leafGraph staticIA startVertex inSimpleGraph =

    -- select postOrder function based on graph type
    let postOrderFunction = if (graphType inGS) == Tree then postOrderTreeTraversal
                         else if (graphType inGS) == SoftWired then POSW.postOrderSoftWiredTraversal
                         else if (graphType inGS) == HardWired then postOrderTreeTraversal
                         else error ("Graph type not implemented: " ++ (show $ graphType inGS))

        -- first traversal on outgroup root
        outgroupRooted = postOrderFunction inGS inData leafGraph staticIA startVertex inSimpleGraph

        -- start at start vertex--for components or ur-root for full graph
        startVertexList = if startVertex == Nothing then fmap fst $ LG.getRoots $ thd6 outgroupRooted
                          else [fromJust startVertex]

        -- next edges (to vertex in list) to perform rerroting
        -- progresses recursivey over adjacent edges to minimize node reoptimization
        -- childrenOfRoot = concatMap (LG.descendants (thd6 outgroupRooted)) startVertexList
        -- grandChildrenOfRoot = concatMap (LG.descendants (thd6 outgroupRooted)) childrenOfRoot

        -- create list of multi-traversals with original rooting first
        -- subsequent rerooting do not reoptimize exact characters (add nonadd etc)
        -- they are taken from the fully labelled first input decorated graph later when output graph created
        -- it is important that the first graph be the ourgroup rooted graph (outgroupRootedPhyloGraph) so this
        -- will have the preorder assignments for the outgroup rooted graph as 3rd field.  This can be used for incremental
        -- optimization to get O(log n) initial postorder assingment when mutsating graph.
        -- hardwired reroot cause much pain
        -- the head startvertex list for reoptimizing spit trees ni swapping
        recursiveRerootList = if (graphType inGS == HardWired) then [outgroupRooted]
                              else if (graphType inGS == SoftWired) then [POSW.getDisplayBasedRerootSoftWired inGS SoftWired (head startVertexList) outgroupRooted]
                             -- need to test which is better
                              -- else if (graphType inGS == SoftWired) then outgroupRooted : minimalReRootPhyloGraph inGS outgroupRooted (head startVertexList) grandChildrenOfRoot
                              else if (graphType inGS == Tree) then [POSW.getDisplayBasedRerootSoftWired inGS Tree (head startVertexList) outgroupRooted]
                              --else if (graphType inGS == Tree) then outgroupRooted : minimalReRootPhyloGraph inGS outgroupRooted (head startVertexList) grandChildrenOfRoot
                              else error ("Graph type not implemented: " ++ (show $ graphType inGS))


        -- remove if tree reroot code pans out
        finalizedPostOrderGraphList = L.sortOn snd6 recursiveRerootList

        -- create optimal final graph with best costs and best traversal (rerooting) forest for each character
        -- traversal for exact characters (and costs) are the first of each least since exact only optimizaed for that
        -- traversal graph.  The result has approprotate post-order assignments for traversals, preorder "final" assignments
        -- are propagated to the Decorated graph field after the preorder pass.
        -- doesn't have to be sorted, but should minimize assignments
        graphWithBestAssignments = head recursiveRerootList -- L.foldl1' setBetterGraphAssignment recursiveRerootList'

        {-  root and model complexities moved to output
        -- same root cost if same data and number of roots
        localRootCost = rootComplexity inGS 
                        {-if (rootCost inGS) == NoRootCost then 0.0
                        else if (rootCost inGS) == Wheeler2015Root then getW15RootCost inGS outgroupRooted
                        else error ("Root cost type " ++ (show $ rootCost inGS) ++ " is not yet implemented")
                        -}
        -}

    in
    -- trace ("GPOT length: " ++ (show $ fmap snd6 recursiveRerootList) ++ " " ++ (show $ graphType inGS)) (
    -- trace ("TRAV:" ++ (show startVertex) ++ " " ++ (show sequenceChars) ++ " " ++ (show (snd6 outgroupRooted, fmap snd6 finalizedPostOrderGraphList, snd6 graphWithBestAssignments))
    --    ++ "\nTraversal root costs: " ++ (show (getTraversalCosts outgroupRooted, fmap getTraversalCosts recursiveRerootList', getTraversalCosts graphWithBestAssignments))) (

    -- only static characters
    if sequenceChars == 0 then
        let penaltyFactor  = if (graphType inGS == Tree) then 0.0

                             --it is its own penalty due to counting all changes in in2 out 1 nodes
                             else if (graphType inGS == HardWired) then 0.0

                             -- softwired versions
                             else if (graphFactor inGS) == NoNetworkPenalty then 0.0
                             else if (graphFactor inGS) == Wheeler2015Network then getW15NetPenalty startVertex outgroupRooted
                             else if (graphFactor inGS) == Wheeler2023Network then getW23NetPenalty startVertex outgroupRooted

                             else error ("Network penalty type " ++ (show $ graphFactor inGS) ++ " is not yet implemented")

            staticOnlyGraph = if (graphType inGS) == SoftWired then POSW.updateAndFinalizePostOrderSoftWired startVertex (head startVertexList) outgroupRooted
                              else outgroupRooted
            -- staticOnlyGraph = head recursiveRerootList'
            staticOnlyGraph' = if startVertex == Nothing then updatePhylogeneticGraphCost staticOnlyGraph (penaltyFactor + (snd6 staticOnlyGraph))
                               else updatePhylogeneticGraphCost staticOnlyGraph (penaltyFactor + (snd6 staticOnlyGraph))
        in
        -- trace ("Only static: " ++ (snd6 staticOnlyGraph'))
        (staticOnlyGraph', head startVertexList)

    -- single seuquence (prealigned, dynamic) only (ie no static)
    else if sequenceChars == 1 && (U.getNumberExactCharacters (thd3 inData) == 0) then
        let penaltyFactorList  = if (graphType inGS == Tree) then replicate (length finalizedPostOrderGraphList) 0.0
                                 else if (graphType inGS == HardWired) then replicate (length finalizedPostOrderGraphList) 0.0
                                 else if (graphFactor inGS) == NoNetworkPenalty then replicate (length finalizedPostOrderGraphList) 0.0
                                 else if (graphFactor inGS) == Wheeler2015Network then fmap (getW15NetPenalty startVertex) finalizedPostOrderGraphList
                                 else if (graphFactor inGS) == Wheeler2023Network then fmap (getW23NetPenalty startVertex) finalizedPostOrderGraphList
                                 else error ("Network penalty type " ++ (show $ graphFactor inGS) ++ " is not yet implemented")
            newCostList = zipWith (+) penaltyFactorList (fmap snd6 finalizedPostOrderGraphList)


            finalizedPostOrderGraph = head $ L.sortOn snd6 $ zipWith updatePhylogeneticGraphCost finalizedPostOrderGraphList newCostList
        in
        -- trace ("GPOT-1: " ++ (show (snd6 finalizedPostOrderGraph)))
        (finalizedPostOrderGraph, head startVertexList)

    -- multiple dynamic characters--checks for best root for each character
    -- important to have outgroup rooted graph first for fold so don't use sorted recursive list
    else
        let penaltyFactor  = if (graphType inGS == Tree) then 0.0
                             else if (graphType inGS == HardWired) then 0.0
                             else if (graphFactor inGS) == NoNetworkPenalty then 0.0
                             else if (graphFactor inGS) == Wheeler2015Network then getW15NetPenalty startVertex graphWithBestAssignments
                             else if (graphFactor inGS) == Wheeler2023Network then getW23NetPenalty startVertex graphWithBestAssignments
                             else error ("Network penalty type " ++ (show $ graphFactor inGS) ++ " is not yet implemented")

            graphWithBestAssignments' = updatePhylogeneticGraphCost graphWithBestAssignments (penaltyFactor + (snd6 graphWithBestAssignments))
        in
        -- trace ("GPOT-2: " ++ (show (penaltyFactor + (snd6 graphWithBestAssignments))))
        (graphWithBestAssignments', head startVertexList)

    -- )

-- | updateGraphCostsComplexities adds root and model complexities if appropriate to graphs
updateGraphCostsComplexities :: GlobalSettings -> PhylogeneticGraph -> PhylogeneticGraph
updateGraphCostsComplexities inGS inGraph = 
    if optimalityCriterion inGS == Parsimony then inGraph
    else if optimalityCriterion inGS == Likelihood then
        -- trace ("\tFinalizing graph cost with root priors")
        updatePhylogeneticGraphCost inGraph ((rootComplexity inGS) +  (snd6 inGraph))
    else if optimalityCriterion inGS == PMDL then
        -- trace ("\tFinalizing graph cost with model and root complexities")
        updatePhylogeneticGraphCost inGraph ((rootComplexity inGS) + (modelComplexity inGS) + (snd6 inGraph))
    else error ("Optimality criterion not recognized/implemented: " ++ (show $ optimalityCriterion inGS))

-- | updatePhylogeneticGraphCost takes a PhylgeneticGrtaph and Double and replaces the cost (snd of 6 fields)
-- and returns Phylogenetic graph
updatePhylogeneticGraphCost :: PhylogeneticGraph -> VertexCost -> PhylogeneticGraph
updatePhylogeneticGraphCost (a, _, b, c, d, e) newCost = (a, newCost, b, c, d, e)

-- | getW15RootCost creates a root cost as the 'insertion' of character data.  For sequence data averaged over
-- leaf taxa
getW15RootCost :: GlobalSettings -> PhylogeneticGraph -> VertexCost
getW15RootCost inGS inGraph =
    if LG.isEmpty $ thd6 inGraph then 0.0
    else
        let (rootList, _, _, _) = LG.splitVertexList $ fst6 inGraph
            numRoots = length rootList
            
        in
        (fromIntegral numRoots) * (rootComplexity inGS)

-- | getW15NetPenalty takes a Phylogenetic tree and returns the network penalty of Wheeler (2015)
-- modified to take the union of all edges of trees of minimal length
-- currently modified -- not exactlty W15
getW15NetPenalty :: Maybe Int -> PhylogeneticGraph -> VertexCost
getW15NetPenalty startVertex inGraph =
    if LG.isEmpty $ thd6 inGraph then 0.0
    else
        let (bestTreeList, _) = extractLowestCostDisplayTree startVertex inGraph
            bestTreesEdgeList = L.nubBy undirectedEdgeEquality $ concat $ fmap LG.edges bestTreeList
            rootIndex = if startVertex == Nothing then fst $ head $ LG.getRoots (fst6 inGraph)
                        else fromJust startVertex
            blockPenaltyList = PU.seqParMap rdeepseq  (getBlockW2015 bestTreesEdgeList rootIndex) (fth6 inGraph)

            -- leaf list for normalization
            (_, leafList, _, _) = LG.splitVertexList (fst6 inGraph)
            numLeaves = length leafList
            numTreeEdges = 4.0 * (fromIntegral numLeaves) - 4.0
            divisor = numTreeEdges
        in
        -- trace ("W15:" ++ (show ((sum $ blockPenaltyList) / divisor )) ++ " from " ++ (show (numTreeEdges, numExtraEdges, divisor, sum blockPenaltyList))) (
        (sum $ blockPenaltyList) / divisor
        -- )

-- | getW23NetPenalty takes a Phylogenetic tree and returns the network penalty of Wheeler (2023)
-- basic idea is new edge improvement must be better than average existing edge cost
-- penalty for each added edge (unlike W15 which was on a block by block basis)
-- num extra edges/2 since actually add 2 new edges when one network edge
getW23NetPenalty :: Maybe Int -> PhylogeneticGraph -> VertexCost
getW23NetPenalty startVertex inGraph =
    if LG.isEmpty $ thd6 inGraph then 0.0
    else
        let (bestTreeList, _) = extractLowestCostDisplayTree startVertex inGraph
            bestTreesEdgeList = L.nubBy undirectedEdgeEquality $ concat $ fmap LG.edges bestTreeList
            
            -- rootIndex = if startVertex == Nothing then fst $ head $ LG.getRoots (fst6 inGraph)
            --            else fromJust startVertex
            
            -- blockPenaltyList = PU.seqParMap rdeepseq (getBlockW2015 bestTreesEdgeList rootIndex) (fth6 inGraph)

            -- leaf list for normalization
            (_, leafList, _, _) = LG.splitVertexList (fst6 inGraph)
            numLeaves = length leafList
            numTreeEdges = 2.0 * (fromIntegral numLeaves) - 2.0
            numExtraEdges = ((fromIntegral $ length bestTreesEdgeList) - numTreeEdges) / 2.0
            divisor = numTreeEdges - numExtraEdges
        in
        -- trace ("W23:" ++ (show ((numExtraEdges * (snd6 inGraph)) / (2.0 * numTreeEdges))) ++ " from " ++ (show (numTreeEdges, numExtraEdges))) (
        if divisor == 0.0 then infinity
        -- else (sum blockPenaltyList) / divisor
        -- else (numExtraEdges * (sum blockPenaltyList)) / divisor
        else (numExtraEdges * (snd6 inGraph)) / (2.0 * numTreeEdges)
        -- )


-- | getBlockW2015 takes the list of trees for a block, gets the root cost and determines the individual
-- penalty cost of that block
getBlockW2015 :: [LG.Edge] -> Int -> [DecoratedGraph] -> VertexCost
getBlockW2015 treeEdgeList rootIndex blockTreeList =
    if null treeEdgeList || null blockTreeList then 0.0
    else
        let blockTreeEdgeList = L.nubBy undirectedEdgeEquality $ concatMap LG.edges blockTreeList
            numExtraEdges = length $ undirectedEdgeMinus blockTreeEdgeList treeEdgeList
            blockCost = subGraphCost $ fromJust $ LG.lab (head blockTreeList) rootIndex
        in
        -- trace ("GBW: " ++ (show (numExtraEdges, blockCost, blockTreeEdgeList)) ++ "\n" ++ (show $ fmap (subGraphCost . snd) $ LG.labNodes (head blockTreeList)))
        blockCost * (fromIntegral numExtraEdges)

-- | checkUnusedEdgesPruneInfty checks if a softwired phylogenetic graph has
-- "unused" edges sensu Wheeler 2015--that an edge in the canonical graph is
-- not present in any of the block display trees (that are heurstically optimal)
-- the options specify if the cost returned is Infinity (really max bound Double)
-- with no pruning of edges or the cost is left unchanged and unused edges are
-- pruned from the canonical graph
-- this is unDirected due to rerooting heuristic in post/preorder optimization
-- inifinity defined in Types.hs
checkUnusedEdgesPruneInfty :: GlobalSettings -> ProcessedData -> Bool -> Bool -> DecoratedGraph-> PhylogeneticGraph -> PhylogeneticGraph
checkUnusedEdgesPruneInfty inGS inData pruneEdges warnPruneEdges leafGraph inGraph@(inSimple, _, inCanonical, blockTreeV, charTreeVV, charInfoVV) =
    let simpleEdgeList = LG.edges inSimple
        displayEdgeSet =  L.nubBy undirectedEdgeEquality $ concat $ concat $ fmap (fmap LG.edges) blockTreeV
        unusedEdges = undirectedEdgeMinus simpleEdgeList displayEdgeSet
    in
    -- no unused edges all OK
    if null unusedEdges then inGraph

    -- unused edges--do not prune return "infinite cost"
    else if not pruneEdges then 
        -- trace ("Unused edge->Infinity") 
        (inSimple, infinity, inCanonical, blockTreeV, charTreeVV, charInfoVV)

    -- unused but pruned--need to prune nodes and reoptimize to get final assignments correct
    else
        let newSimpleGraph = LG.delEdges unusedEdges inSimple
            contractedSimple = GO.contractIn1Out1EdgesRename newSimpleGraph
        in
        if warnPruneEdges then
            trace ("Pruning " ++ (show $ length unusedEdges) ++ " unused edges and reoptimizing graph")
            multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph Nothing contractedSimple

        else multiTraverseFullyLabelSoftWired inGS inData pruneEdges warnPruneEdges leafGraph Nothing contractedSimple

-- | undirectedEdgeEquality checks edgse for equality irrespective of direction
undirectedEdgeEquality :: LG.Edge -> LG.Edge -> Bool
undirectedEdgeEquality (a,b) (c,d) = if a == c && b == d then True
                                               else if a == d && b == c then True
                                               else False

-- | undirectedEdgeMinus subtracts edges in the second list from those in the first using
-- undirected matching
undirectedEdgeMinus :: [LG.Edge] -> [LG.Edge] -> [LG.Edge]
undirectedEdgeMinus firstList secondList =
    if null firstList then []
    else
        let firstEdge@(a,b) = head firstList
        in
        if firstEdge `L.elem` secondList then undirectedEdgeMinus (tail firstList) secondList
        else if (b,a) `L.elem` secondList then undirectedEdgeMinus (tail firstList) secondList
        else firstEdge : undirectedEdgeMinus (tail firstList) secondList

-- | extractLowestCostDisplayTree takes a phylogenetic graph and takes all valid (complete) resolutions
-- (display trees) and their costs
-- and determines the total cost (over all blocks) of each display tree
-- the lowest cost display tree(s) as list are returned with cost
-- this is used in Wheeler (2015) network penalty
extractLowestCostDisplayTree :: Maybe Int -> PhylogeneticGraph -> ([DecoratedGraph], VertexCost)
extractLowestCostDisplayTree startVertex inGraph =
 if LG.isEmpty $ thd6 inGraph then error "Empty graph in extractLowestCostDisplayTree"
 else
    let -- get componen t or global root label
        rootLabel = if startVertex == Nothing then snd $ head $ LG.getRoots (thd6 inGraph)
                    else fromJust $ LG.lab (thd6 inGraph) (fromJust startVertex)

        -- get resolution data for start/rpoot vertex
        blockResolutionLL = V.toList $ fmap POSW.getAllResolutionList (vertexResolutionData rootLabel)
        --blockResolutionLL = V.toList $ fmap (PO.getBestResolutionListPair startVertex False) (vertexResolutionData rootLabel)
        displayTreeBlockList = L.transpose blockResolutionLL
        displayTreePairList = L.foldl1' sumTreeCostLists displayTreeBlockList
        minimumCost = minimum $ fmap snd displayTreePairList
        (bestDisplayTreeList, _) = unzip $ filter ((== minimumCost) . snd) displayTreePairList
    in
    -- trace ("FC: " ++ (show $ fmap snd displayTreePairList))
    (bestDisplayTreeList, minimumCost)

-- | sumTreeCostLists takes two lists of (Graph, Cost) pairs and sums the costs and keeps the trees the same
-- does not check that graphs are the same after debug
sumTreeCostLists :: (Eq a, Eq b) => [(LG.Gr a b, VertexCost)] ->  [(LG.Gr a b, VertexCost)] ->  [(LG.Gr a b, VertexCost)]
sumTreeCostLists firstList secondList =
    if null firstList || null secondList then error "Empty list in sumTreeCostLists"
    else
        let (firstGraphList, firstCostList) = unzip firstList
            (secondGraphList, secondCostList) =  unzip secondList
            newCostList = zipWith (+)  firstCostList secondCostList

            -- remove once working
            checkList = filter (== False) $ zipWith LG.equal firstGraphList secondGraphList
        in
        if null checkList then error ("Graph lists not same : " ++ (show checkList))
        else
            -- trace ("Graphs match ")
            zip firstGraphList newCostList

-- | setBetterGraphAssignment takes two phylogenetic graphs and returns the lower cost optimization of each character,
-- with traversal focus etc to get best overall graph
-- since this is meant to work with graphs that have or do not have reoptimized exact (=static-Add/NonAdd/MAtrix) characters
-- the criterion is lower cost character is taken, unless the cost is zero, then non-zero is taken
-- this function is expected to be used in a fold over a list of graphs
-- the basic comparison is over the costs of the root(s) cost for each  of the character decorated (traversal) graphs

-- May change
-- assumes that a single decorated graph comes in for each Phylogenetic graph from the fully and reroot optimize (V.singleton (V.singleton DecGraph))
-- and goes through the block-character-cost data and reassigns based on that creating a unique (although there could be more than one) decorated
-- graph for each character in each block.
-- postorder assignments in traversal set of block character trees are NOT propagated back to first decorated graph.
-- the third field of phylogenetic Graph is set to the 3rd fieled of the first of two inputs--so if startiong fold with outgroup
-- rooted graph--that is what stays which can be used as a preorder graph for incremental optimization
-- when perfoming that sort of operation
-- The traversal graphs
-- are used for the pre-order final assignments which will be propagated back to set those of the 3rd field decorated graph

-- this will have to be modified for solf-wired since incoming blocks will not all be the same underlying gaph
-- unclear how hardwired will be affected
setBetterGraphAssignment :: PhylogeneticGraph -> PhylogeneticGraph -> PhylogeneticGraph
setBetterGraphAssignment firstGraph@(fSimple, _, fDecGraph, fBlockDisplay, fTraversal, fCharInfo) secondGraph@(_, _, sDecGraph, _, sTraversal, _) =
    -- trace ("SBGA:" ++  (show $ (length  fTraversal, length sTraversal))) (
    if LG.isEmpty fDecGraph then secondGraph
    else if LG.isEmpty sDecGraph then firstGraph
    else
        -- trace ("setBetter (" ++ (show fCost) ++ "," ++ (show sCost) ++ ")"  ++ " CharInfo blocks:" ++ (show $ length fCharInfo) ++ " characters: " ++ (show $ fmap length fCharInfo) ++ " "
        --     ++ (show $ fmap (fmap name) fCharInfo)) (
        let (mergedBlockVect, costVector) = V.unzip $ V.zipWith makeBetterBlock fTraversal sTraversal
        in
         --trace ("setBetter (" ++ (show fCost) ++ "," ++ (show sCost) ++ ") ->" ++ (show $ V.sum costVector) ++ " nt:" ++ (show $ length fTraversal)
         --   ++ "length blocks " ++ (show $ fmap length fTraversal))
        (fSimple, V.sum costVector, fDecGraph, fBlockDisplay, mergedBlockVect, fCharInfo)
        -- )

-- | makeBetterBlocktakes two verctors of traversals. Each vector contains a decorated graph (=traversla graph) for each
-- character.  This can be a single sequence or series of exact characters
-- the function applies a character cost comparison to get the better
makeBetterBlock :: V.Vector DecoratedGraph -> V.Vector DecoratedGraph -> (V.Vector DecoratedGraph, VertexCost)
makeBetterBlock firstBlockGraphVect secondBlockGraphVect =
    let (mergedCharacterVect, costVector) = V.unzip $ V.zipWith chooseBetterCharacter firstBlockGraphVect secondBlockGraphVect
    in
    -- trace ("MBB: " ++ (show $ (length  firstBlockGraphVect, length firstBlockGraphVect)))
    (mergedCharacterVect, V.sum costVector)

-- | chooseBetterCharacter takes a pair of character decorated graphs and chooses teh "better" one as in lower cost, or non-zero cost
--  if one is zer (for exact characters) and returns the better character and cost
--  graph can have multiple roots
chooseBetterCharacter :: DecoratedGraph -> DecoratedGraph -> (DecoratedGraph, VertexCost)
chooseBetterCharacter firstGraph secondGraph
  | LG.isEmpty firstGraph = error "Empty first graph in chooseBetterCharacter"
  | LG.isEmpty secondGraph = error "Empty second graph in chooseBetterCharacter"
  | otherwise =
    let firstGraphCost = sum $ fmap (subGraphCost . snd) (LG.getRoots firstGraph)
        secondGraphCost = sum $ fmap (subGraphCost . snd) (LG.getRoots secondGraph)
    in
    -- trace ("Costs " ++ show (firstGraphCost, secondGraphCost)) (
    if firstGraphCost == 0 then (secondGraph, secondGraphCost)
    else if secondGraphCost == 0 then (firstGraph, firstGraphCost)
    else if secondGraphCost < firstGraphCost then (secondGraph, secondGraphCost)
    else (firstGraph, firstGraphCost)
    -- )

-- | minimalReRootPhyloGraph takes an inialtial post-order labelled phylogenetic graph
-- and "intelligently" reroots by traversing through adjacent edges, hopefully
-- reoptimizing the minimum number of vertices each time (2) but could be more depending
-- on graph topology
-- NB--only deals with post-order assignments
minimalReRootPhyloGraph :: GlobalSettings -> PhylogeneticGraph -> Int -> [LG.Node] -> [PhylogeneticGraph]
minimalReRootPhyloGraph inGS inGraph originalRoot nodesToRoot =
    -- trace ("MRR: " ++ (show nodesToRoot) ++ " " ++ (show $ fmap (LG.descendants (thd6 inGraph)) nodesToRoot)) (
    if null nodesToRoot then []
    else
        let firstRerootIndex = head nodesToRoot
            nextReroots = LG.descendants (thd6 inGraph) firstRerootIndex ++ tail nodesToRoot
            newGraph
              | (graphType inGS)  == Tree = PO.rerootPhylogeneticGraph' inGS False False inGraph originalRoot firstRerootIndex
              | (graphType inGS)  == SoftWired = PO.rerootPhylogeneticNetwork' inGS inGraph originalRoot firstRerootIndex
              | (graphType inGS)  == HardWired = PO.rerootPhylogeneticNetwork' inGS inGraph originalRoot firstRerootIndex
              | otherwise = errorWithoutStackTrace ("Graph type not implemented/recognized: " ++ show (graphType inGS))
        in
        -- trace ("NRR: " ++ " " ++ (show (LG.descendants (thd6 inGraph) firstRerootIndex)) ) ( -- ++ " -> " ++ (show nextReroots) ++ "\n" ++ (LG.prettify $ fst6 inGraph) ++ "\n" ++ (LG.prettify $ fst6 newGraph)) (
        -- trace ("MRR: " ++ (show $ snd6 inGraph) ++ " -> " ++ (show $ snd6 newGraph)) (
        if fst6 newGraph == LG.empty then minimalReRootPhyloGraph inGS inGraph originalRoot nextReroots
        else newGraph : minimalReRootPhyloGraph inGS newGraph originalRoot nextReroots
        -- ) )


-- | makeLeafGraph takes input data and creates a 'graph' of leaves with Vertex informnation
-- but with zero edges.  This 'graph' can be reused as a starting structure for graph construction
-- to avoid remaking of leaf vertices
makeLeafGraph :: ProcessedData -> DecoratedGraph
makeLeafGraph (nameVect, bvNameVect, blocDataVect) =
    if V.null nameVect then error "Empty ProcessedData in makeLeafGraph"
    else
        let leafVertexList = V.toList $ V.map (makeLeafVertex nameVect bvNameVect blocDataVect) (V.fromList [0.. V.length nameVect - 1])
        in
        LG.mkGraph leafVertexList []


-- | makeSimpleLeafGraph takes input data and creates a 'graph' of leaves with Vertex informnation
-- but with zero edges.  This 'graph' can be reused as a starting structure for graph construction
-- to avoid remaking of leaf vertices
makeSimpleLeafGraph :: ProcessedData -> SimpleGraph
makeSimpleLeafGraph (nameVect, _, _) =
    if V.null nameVect then error "Empty ProcessedData in makeSimpleLeafGraph"
    else
        let leafVertexList = V.toList $ V.map (makeSimpleLeafVertex nameVect) (V.fromList [0.. V.length nameVect - 1])
        in
        LG.mkGraph leafVertexList []
        where makeSimpleLeafVertex a b = (b, a V.! b)


-- | makeLeafVertex makes a single unconnected vertex for a leaf
makeLeafVertex :: V.Vector NameText -> V.Vector NameBV -> V.Vector BlockData -> Int -> LG.LNode VertexInfo
makeLeafVertex nameVect bvNameVect inData localIndex =
    -- trace ("Making leaf " ++ (show localIndex) ++ " Data " ++ (show $ length inData) ++ " " ++ (show $ fmap length $ fmap snd3 inData)) (
    let centralData = V.map snd3 inData
        thisData = V.map (V.! localIndex) centralData
        newVertex = VertexInfo  { index = localIndex
                                , bvLabel = bvNameVect V.! localIndex
                                , parents = V.empty
                                , children = V.empty
                                , nodeType = LeafNode
                                , vertName =  nameVect V.! localIndex
                                , vertData = thisData
                                , vertexResolutionData = mempty
                                , vertexCost = 0.0
                                , subGraphCost = 0.0
                                }
        in
        -- trace (show (length thisData) ++ (show $ fmap length thisData))
        (localIndex, newVertex)
        -- )

-- | postOrderTreeTraversal takes a 'simple' graph and generates 'preliminary' assignments
-- vi post-order traversal, yields cost as well
-- for a binary tree only
-- depending on optimality criterion--will calculate root cost
postOrderTreeTraversal :: GlobalSettings ->  ProcessedData -> DecoratedGraph -> Bool -> Maybe Int -> SimpleGraph -> PhylogeneticGraph
postOrderTreeTraversal _ (_, _, blockDataVect) leafGraph staticIA startVertex inGraph  =
    if LG.isEmpty inGraph then emptyPhylogeneticGraph
    else
        -- Assumes root is Number of Leaves
        let rootIndex = if startVertex == Nothing then  fst $ head $ LG.getRoots inGraph
                        else fromJust startVertex
            blockCharInfo = V.map thd3 blockDataVect
            newTree = postDecorateTree staticIA inGraph leafGraph blockCharInfo rootIndex rootIndex
        in
        -- trace ("It Begins at " ++ (show $ fmap fst $ LG.getRoots inGraph) ++ "\n" ++ show inGraph) (
        if (startVertex == Nothing) && (not $ LG.isRoot inGraph rootIndex) then
            let localRootList = fst <$> LG.getRoots inGraph
                localRootEdges = concatMap (LG.out inGraph) localRootList
                currentRootEdges = LG.out inGraph rootIndex
            in
            error ("Index "  ++ show rootIndex ++ " with edges " ++ show currentRootEdges ++ " not root in graph:" ++ show localRootList ++ " edges:" ++ show localRootEdges ++ "\n" ++ GFU.showGraph inGraph)
        else newTree
        --)


-- | postDecorateTree' is wrapper for postDecorateTree to alow for mapping
postDecorateTree' :: Bool -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.Node -> LG.Node -> SimpleGraph -> PhylogeneticGraph
postDecorateTree' staticIA curDecGraph blockCharInfo rootIndex curNode simpleGraph = postDecorateTree staticIA simpleGraph curDecGraph blockCharInfo rootIndex curNode

-- | postDecorateTree begins at start index (usually root, but could be a subtree) and moves preorder till children are labelled and then returns postorder
-- labelling vertices and edges as it goes back to root
-- this for a tree so single root
postDecorateTree :: Bool -> SimpleGraph -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.Node -> LG.Node -> PhylogeneticGraph
postDecorateTree staticIA simpleGraph curDecGraph blockCharInfo rootIndex curNode =
    -- if node in there (leaf) nothing to do and return
    if LG.gelem curNode curDecGraph then
        let nodeLabel = LG.lab curDecGraph curNode
        in
        if isNothing nodeLabel then error ("Null label for node " ++ show curNode)
        else
            -- checks for node already in graph--either leaf or pre-optimized node in Hardwired
            -- trace ("In graph :" ++ (show curNode) ++ " " ++ (show nodeLabel))
            (simpleGraph, subGraphCost (fromJust nodeLabel), curDecGraph, mempty, mempty, blockCharInfo)

    -- Need to make node
    else

        -- check if children in graph
        let nodeChildren = LG.descendants simpleGraph curNode  -- should be 1 or 2, not zero since all leaves already in graph
            leftChild = head nodeChildren
            rightChild = last nodeChildren
            leftChildTree = postDecorateTree staticIA simpleGraph curDecGraph blockCharInfo rootIndex leftChild
            rightLeftChildTree = if length nodeChildren == 2 then postDecorateTree staticIA simpleGraph (thd6 leftChildTree) blockCharInfo rootIndex rightChild
                                 else leftChildTree
            newSubTree = thd6 rightLeftChildTree
            (leftChildLabel, rightChildLabel) = U.leftRightChildLabelBV (fromJust $ LG.lab newSubTree leftChild, fromJust $ LG.lab newSubTree rightChild)

        in

        if length nodeChildren > 2 then error ("Graph not dichotomous in postDecorateTree node " ++ show curNode ++ "\n" ++ LG.prettify simpleGraph)
        else if null nodeChildren then error ("Leaf not in graph in postDecorateTree node " ++ show curNode ++ "\n" ++ LG.prettify simpleGraph)

        -- out-degree 1 should not happen with Tree but will with HardWired graph
        else if length nodeChildren == 1 then
            -- make node from single child and single new edge to child
            -- takes characters in blocks--but for tree really all same block
            let childVertexData = vertData leftChildLabel
                newVertex = VertexInfo {  index = curNode
                                        -- same as child--could and perhaps should prepend 1 to make distinct
                                        , bvLabel = bvLabel leftChildLabel
                                        , parents = V.fromList $ LG.parents simpleGraph curNode
                                        , children = V.fromList nodeChildren
                                        , nodeType = GO.getNodeType simpleGraph curNode
                                        , vertName = T.pack $ "HTU" ++ show curNode
                                        , vertData = childVertexData
                                        -- this not used for Hardwired or Tree
                                        , vertexResolutionData = mempty
                                        , vertexCost = 0.0
                                        , subGraphCost = subGraphCost leftChildLabel
                                        }
                newEdgesLabel = EdgeInfo {    minLength = 0.0
                                            , maxLength = 0.0
                                            , midRangeLength = 0.0
                                            , edgeType = TreeEdge
                                         }
                newEdges = LG.toEdge <$> LG.out simpleGraph curNode
                newLEdges =  fmap (LG.toLEdge' newEdgesLabel) newEdges
                newGraph =  LG.insEdges newLEdges $ LG.insNode (curNode, newVertex) newSubTree

                (newDisplayVect, newCharTreeVV) = POSW.divideDecoratedGraphByBlockAndCharacterTree newGraph

            in
            -- th curnode == roiot index for pruned subtrees
            -- trace ("New vertex:" ++ (show newVertex) ++ " at cost " ++ (show newCost)) (
            -- Do we need to PO.divideDecoratedGraphByBlockAndCharacterTree if not root?  probbaly not

            --if nodeType newVertex == RootNode then (simpleGraph, subGraphCost newVertex, newGraph, mempty, PO.divideDecoratedGraphByBlockAndCharacterTree newGraph, blockCharInfo)
            if nodeType newVertex == RootNode || curNode == rootIndex then (simpleGraph, subGraphCost newVertex, newGraph, newDisplayVect, newCharTreeVV, blockCharInfo)
            else (simpleGraph, subGraphCost newVertex, newGraph, mempty, mempty, blockCharInfo)

        -- make node from 2 children
        else
            -- make node from children and new edges to children
            -- takes characters in blocks--but for tree really all same block
            let -- this ensures that left/right choices are based on leaf BV for consistency and label invariance
                -- larger bitvector is Right, smaller or equal Left

                newCharData = if staticIA then PO.createVertexDataOverBlocksStaticIA  (vertData leftChildLabel) (vertData  rightChildLabel) blockCharInfo []
                              else PO.createVertexDataOverBlocks  (vertData leftChildLabel) (vertData  rightChildLabel) blockCharInfo []
                newCost =  V.sum $ V.map V.sum $ V.map (V.map snd) newCharData
                newVertex = VertexInfo {  index = curNode
                                        , bvLabel = bvLabel leftChildLabel .|. bvLabel rightChildLabel
                                        , parents = V.fromList $ LG.parents simpleGraph curNode
                                        , children = V.fromList nodeChildren
                                        , nodeType = GO.getNodeType simpleGraph curNode
                                        , vertName = T.pack $ "HTU" ++ show curNode
                                        , vertData = V.map (V.map fst) newCharData
                                        , vertexResolutionData = mempty
                                        , vertexCost = newCost
                                        , subGraphCost = subGraphCost leftChildLabel + subGraphCost rightChildLabel + newCost
                                        }
                newEdgesLabel = EdgeInfo {    minLength = newCost / 2.0
                                            , maxLength = newCost / 2.0
                                            , midRangeLength = newCost / 2.0
                                            , edgeType = TreeEdge
                                         }
                newEdges = LG.toEdge <$> LG.out simpleGraph curNode
                newLEdges =  fmap (LG.toLEdge' newEdgesLabel) newEdges
                newGraph =  LG.insEdges newLEdges $ LG.insNode (curNode, newVertex) newSubTree

                (newDisplayVect, newCharTreeVV) = POSW.divideDecoratedGraphByBlockAndCharacterTree newGraph

            in
            -- th curnode == roiot index for pruned subtrees
            -- trace ("New vertex:" ++ (show newVertex) ++ " at cost " ++ (show newCost)) (
            -- Do we need to PO.divideDecoratedGraphByBlockAndCharacterTree if not root?  probbaly not

            --if nodeType newVertex == RootNode then (simpleGraph, subGraphCost newVertex, newGraph, mempty, PO.divideDecoratedGraphByBlockAndCharacterTree newGraph, blockCharInfo)
            if nodeType newVertex == RootNode || curNode == rootIndex then (simpleGraph, subGraphCost newVertex, newGraph, newDisplayVect, newCharTreeVV, blockCharInfo)
            else (simpleGraph, subGraphCost newVertex, newGraph, mempty, mempty, blockCharInfo)

            -- ) -- )


