{- |
Module      :  PostOrderSoftWiredFunctions.hs
Description :  Module specifying post-order softwiired graph functions
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

module GraphOptimization.PostOrderSoftWiredFunctions  ( updateAndFinalizePostOrderSoftWired
                                                      , postOrderSoftWiredTraversal
                                                      , postDecorateSoftWired'
                                                      , postDecorateSoftWired
                                                      , assignPostOrderToDisplayTree
                                                      , softWiredPostOrderTraceBack
                                                      , makeLeafGraphSoftWired
                                                      , getOutDegree1VertexAndGraph
                                                      , getOutDegree1VertexSoftWired
                                                      , getOutDegree2VertexSoftWired
                                                      , extractDisplayTrees
                                                      , createBlockResolutions
                                                      , makeCharacterGraph
                                                      , getAllResolutionList
                                                      , getBestResolutionListPair
                                                      , getDisplayBasedRerootSoftWired
                                                      , divideDecoratedGraphByBlockAndCharacterTree
                                                      , postOrderTreeTraversal
                                                      , postDecorateTree
                                                      , postDecorateTree'
                                                      , createVertexDataOverBlocks
                                                      , createVertexDataOverBlocksStaticIA
                                                      , createVertexDataOverBlocksNonExact
                                                      , getBlockCost
                                                      , updatePhylogeneticGraphCost
                                                      , updateGraphCostsComplexities
                                                      , getW15NetPenalty
                                                      , getW23NetPenalty
                                                      , getW15RootCost
                                                      -- included to quiet warnings
                                                      , rerootBlockCharTrees'
                                                      , getCharTreeBestRoot'
                                                      ) where

import           Data.Bits
import qualified Data.BitVector.LittleEndian as BV
import qualified Data.List                   as L
import           Data.Maybe
import qualified Data.Text.Lazy              as T
import qualified Data.Vector                 as V
import           GeneralUtilities
import qualified GraphOptimization.Medians   as M
import qualified Graphs.GraphOperations      as GO
import           Types.Types
import qualified Utilities.LocalGraph        as LG
import qualified Utilities.Utilities         as U
import           Control.Parallel.Strategies
import qualified ParallelUtilities           as PU
import qualified GraphFormatUtilities        as GFU
import qualified GraphOptimization.PreOrderFunctions as PRE
-- import Debug.Debug
import           Debug.Trace


-- | naivePostOrderSoftWiredTraversal produces the post-order result for a softwired graph using
-- a naive algorithm where all display trees are generated and diagnosed, keeping the best results
-- to return a phylogenetic graph
-- any network penalty is not applied as in postOrderSoftWiredTraversal
-- not contracting in=1 out =1 nodes so that indexing will be consistent
-- does not create resolution cache data structures
-- really only for diagnosis and time complexity comparison with resolution cache algorithm
naivePostOrderSoftWiredTraversal :: GlobalSettings -> ProcessedData -> DecoratedGraph -> Maybe Int -> SimpleGraph -> PhylogeneticGraph
naivePostOrderSoftWiredTraversal inGS inData@(_, _, blockDataVect) leafGraph startVertex inSimpleGraph =
        -- this is a lazy list so can be consumed and not an issue with exponential number of Trees
    let contractIn1Out1Nodes = False
        displayTreeList = LG.generateDisplayTrees contractIn1Out1Nodes inSimpleGraph
        
        -- get the best traversals of the best display trees for each block
        (blockCostList, bestDisplayTreeList, charTreeVectList) = unzip3 $ getBestDisplayCharBlockList inGS inData leafGraph startVertex [] displayTreeList
        
        -- extract specific information to create the phylogenetic graph
        graphCost = sum blockCostList
        displayTreeVect = V.fromList bestDisplayTreeList
        charTreeVectVect = V.fromList charTreeVectList

        -- propagate node character assignments back to canoncail graph
        -- (decoratedCanonicalGraph, decoratedDisplayTreeVect) = assignCanonicalNodes inSimpleGraph displayTreeVect charTreeVectVect

        -- create postorder Phylgenetic graph
        postOrderPhyloGraph = (inSimpleGraph, graphCost, GO.convertSimpleToDecoratedGraph inSimpleGraph, fmap (:[]) $ fmap GO.convertSimpleToDecoratedGraph displayTreeVect, charTreeVectVect, (fmap thd3 blockDataVect))

        -- perform preorder pass
        staticIA = False
        calculateBranchLengths = True
        useMap = False
        rootIndex = if isJust startVertex  then fromJust startVertex
                    else fst $ head $ LG.getRoots inSimpleGraph
        hasNonExact = (U.getNumberSequenceCharacters blockDataVect > 0)
        preOrderPhyloGraph = PRE.preOrderTreeTraversal inGS (finalAssignment inGS) staticIA calculateBranchLengths hasNonExact rootIndex useMap postOrderPhyloGraph


        -- add in netowrk penalty if any
        penaltyFactor  = 0.0
                         {- Currently these calcualtions rely on resolution cache structures so skipping
                         
                         if (graphType inGS == Tree) then 0.0
                         else if (graphType inGS == HardWired) then 0.0
                         else if (graphFactor inGS) == NoNetworkPenalty then 0.0
                         else if (graphFactor inGS) == Wheeler2015Network then getW15NetPenalty startVertex preOrderPhyloGraph
                         else if (graphFactor inGS) == Wheeler2023Network then getW23NetPenalty startVertex preOrderPhyloGraph
                         else error ("Network penalty type " ++ (show $ graphFactor inGS) ++ " is not yet implemented")
                         
                         -}

        preOrderPhyloGraph' = updatePhylogeneticGraphCost preOrderPhyloGraph (penaltyFactor + (snd6 preOrderPhyloGraph))
    in
    -- trace ("NPOST: " ++ (show (graphCost, V.length displayTreeVect)) ++ " " ++ (show $ V.length charTreeVectVect) ++ " -> " ++ (show $ fmap V.length charTreeVectVect))
    trace ("Warning: Net penalty setting ignored--no penalty added (currently)")
    preOrderPhyloGraph' 
    
    -- (inSimpleGraph, graphCost, decoratedCanonicalGraph, decoratedDisplayTreeVect, charTreeVectVect, (fmap thd3 blockDataVect))

-- | assignCanonicalNodes assigns charVect node assignment to canonical and display block graphs
-- in essence, the funciton converts a SimpleGraph into  DecoratedGraph by assignment of node information
-- assignCanonicalNodes ::SimpleGraph -> V.Vector SimpleGraph -> V.Vector (V.Vector DecoratedGraph) -> (DecoratedGraph, V.Vector [DecoratedGraph])
-- assignCanonicalNodes inSimple displayTreeVect charVectVect = 
--    (LG.empty, V.empty)


-- | getBestDisplayCharBlockList takes a Tree gets best rootings, compares to input list if there is one and takes better
-- returning triple of block cost, display tree, char vect from better tree
getBestDisplayCharBlockList :: GlobalSettings 
                            -> ProcessedData 
                            -> DecoratedGraph 
                            -> Maybe Int 
                            -> [(VertexCost, SimpleGraph, V.Vector DecoratedGraph)] 
                            -> [SimpleGraph] 
                            -> [(VertexCost, SimpleGraph, V.Vector DecoratedGraph)] 
getBestDisplayCharBlockList inGS inData leafGraph startVertex currentBest displayTreeList =
    if null displayTreeList then currentBest
    else 
        trace ("GBDCBL Trees: " ++ (show $ length displayTreeList)) (
        -- take first graph
        let firstGraph = head displayTreeList
            
            -- diagnose post order as Tree
            staticIA = False
            outgrouDiagnosedTree = postOrderTreeTraversal inGS inData leafGraph staticIA startVertex firstGraph

            -- reroot as Tree
            rootIndex = if isJust startVertex  then fromJust startVertex
                        else fst $ head $ LG.getRoots firstGraph

            multiTraverseTree = getDisplayBasedRerootSoftWired' Tree rootIndex outgrouDiagnosedTree

            -- choose better vs currentBest
            -- this can be folded for a list > 2
            currentBetter = chooseBetterTriple rootIndex currentBest multiTraverseTree
        in
        trace ("GBDCBL: " ++ (show (snd6 outgrouDiagnosedTree, snd6 multiTraverseTree)) ++ "\n" ++ (LG.prettyIndices firstGraph))
        getBestDisplayCharBlockList inGS inData leafGraph startVertex currentBetter (tail displayTreeList)
        )

-- | chooseBetterTriple takes the current best triplet of graph data and compares to Phylogenetic graph
-- and creates a new triple of better block cost, displayGraph for blocks, and character graphs
chooseBetterTriple :: LG.Node
                   -> [(VertexCost, SimpleGraph, V.Vector DecoratedGraph)] 
                   -> PhylogeneticGraph 
                   -> [(VertexCost, SimpleGraph, V.Vector DecoratedGraph)] 
chooseBetterTriple rootIndex inTripleList (inSimpleTree, _, _, _, charTreeVV, _) =

        -- get block costs from input graph
    let blockCostList = V.toList $ fmap (getBlockCost rootIndex) charTreeVV
        newTripleList = zip3 blockCostList (L.replicate (length blockCostList) inSimpleTree)  (V.toList charTreeVV)
    in

    -- if no previous set with current 
    if null inTripleList then newTripleList
       
    -- compare current to previous, take better of two
    else 
        zipWith chooseBetterBlock inTripleList newTripleList
    where chooseBetterBlock (a,b,c) (a',b',c') = if a' > a then (a,b,c) else (a', b',c')

-- | getBlockCost takes a block and returns sum of root character tree costs
getBlockCost :: LG.Node -> V.Vector DecoratedGraph -> VertexCost
getBlockCost rootIndex characterTreeVect = 
    if V.null characterTreeVect then error "Empty character block vector in getBlockCost"
    else V.sum $ fmap (getCharacterTreeCost rootIndex) characterTreeVect

-- | getCharacterTreeCost takes a character Tree and returns costs
getCharacterTreeCost :: LG.Node -> DecoratedGraph -> VertexCost
getCharacterTreeCost rootIndex characterTree =
    if LG.isEmpty characterTree then error "Empty charcter tree in getCharacterTreeCost"
    else (subGraphCost . snd) $ LG.labelNode characterTree rootIndex

-- | naiveGetDisplayBasedRerootSoftWired is the naive (based on all resolution display trees)
-- the work of getDisplayBasedRerootSoftWired' is already done (rerooting and all) in naivePostOrderSoftWiredTraversal
naiveGetDisplayBasedRerootSoftWired ::  GlobalSettings -> GraphType -> LG.Node -> PhylogeneticGraph -> PhylogeneticGraph
naiveGetDisplayBasedRerootSoftWired _ _ _ inPhyloGraph  =
    inPhyloGraph

-- | postOrderSoftWiredTraversal is a wrapper to allow correct function choice for alternate softwired algorithms
postOrderSoftWiredTraversal :: GlobalSettings -> ProcessedData -> DecoratedGraph -> Bool -> Maybe Int -> SimpleGraph -> PhylogeneticGraph
postOrderSoftWiredTraversal inGS inData leafGraph _ startVertex inSimpleGraph =
    if softWiredMethod inGS == ResolutionCache then postOrderSoftWiredTraversal' inGS inData leafGraph startVertex inSimpleGraph 
    else naivePostOrderSoftWiredTraversal inGS inData leafGraph startVertex inSimpleGraph

-- | postOrderSoftWiredTraversal' performs postorder traversal on Soft-wired graph
-- staticIA is ignored--but kept for functional polymorphism
-- ur-root = ntaxa is an invariant
postOrderSoftWiredTraversal' :: GlobalSettings -> ProcessedData -> DecoratedGraph -> Maybe Int -> SimpleGraph -> PhylogeneticGraph
postOrderSoftWiredTraversal' inGS inData@(_, _, blockDataVect) leafGraph startVertex inSimpleGraph =
    if LG.isEmpty inSimpleGraph then emptyPhylogeneticGraph
    else
         -- Assumes root is Number of Leaves
        let rootIndex = if startVertex == Nothing then V.length $ fst3 inData
                        else fromJust startVertex
            blockCharInfo = V.map thd3 blockDataVect
            newSoftWired = postDecorateSoftWired inGS inSimpleGraph leafGraph blockCharInfo rootIndex rootIndex
        in
        --trace ("It Begins at " ++ show rootIndex) (
        -- trace ("POSWT:\n" ++ (LG.prettify inSimpleGraph) ++ "\nVertices:\n" ++ (show $ LG.labNodes $ thd6 newSoftWired)) (
        if (startVertex == Nothing) && (not $ LG.isRoot inSimpleGraph rootIndex) then
            let localRootList = fst <$> LG.getRoots inSimpleGraph
                localRootEdges = concatMap (LG.out inSimpleGraph) localRootList
                currentRootEdges = LG.out inSimpleGraph rootIndex
            in
            error ("Index "  ++ show rootIndex ++ " with edges " ++ show currentRootEdges ++ " not root in graph:" ++ show localRootList ++ " edges:" ++ show localRootEdges ++ "\n" ++ LG.prettify inSimpleGraph)
        else
            -- trace ("POSW:" ++ (show $ fmap V.length $ fft6 newSoftWired))
            newSoftWired
        -- )

-- | getDisplayBasedRerootSoftWired is a wrapper to allow correct function choice for alternate softwired algorithms
getDisplayBasedRerootSoftWired :: GlobalSettings -> GraphType -> LG.Node -> PhylogeneticGraph -> PhylogeneticGraph
getDisplayBasedRerootSoftWired inGS inGraphType rootIndex inPhyloGraph = 
    if softWiredMethod inGS == ResolutionCache then getDisplayBasedRerootSoftWired' inGraphType rootIndex inPhyloGraph
    else naiveGetDisplayBasedRerootSoftWired inGS inGraphType rootIndex inPhyloGraph

-- | getDisplayBasedRerootSoftWired' takes a graph and generates reroot costs for each character of each block
-- based on rerooting the display tree for that block.
-- Written for soft-wired, but could be modified for tree (data split on vertdata not resolution data)
-- this is a differnt approach from that of "tree" where the decorated, canonical tree is rerooted and each character and block 
-- cost determined from that single rerooting
-- this should help avoid the issue of rerooting complex, reticulate graphs and maintaining
-- all the condition (cycles, time consistency etc) that occur.
-- done correcly this should be able to be used for trees (all display trees same = cononical graph) as
-- well as softwired, but not for hardwired where reticulations are maintained.

-- this can be modified for Tree data structres--basically by starting with vertdata initialiiy without 
-- resolutoin data trace back--shouls be more efficeinet in many was than existing code

-- Input display trees are for reporting only and do not contain actual character data so must be "pulled"
-- from cononical Decorated graph (thd field)
-- the list :[] stuff due to potential list of diplay trees not employed here
getDisplayBasedRerootSoftWired' :: GraphType -> LG.Node -> PhylogeneticGraph -> PhylogeneticGraph
getDisplayBasedRerootSoftWired' inGraphType rootIndex inPhyloGraph@(a,b,decGraph,_,_,f)  = 
    if LG.isEmpty (fst6 inPhyloGraph) then inPhyloGraph
    else 
        let -- update with pass to retrieve vert data from resolution data
            -- Trfee allready has data in vertData field
            (inSimpleGraph, _, inDecGraph, inBlockGraphV', inBlockCharGraphVV', charInfoVV) = if inGraphType == Tree then 
                                                                                                let (displayTrees, charTrees) = divideDecoratedGraphByBlockAndCharacterTree decGraph
                                                                                                in
                                                                                                (a, b, decGraph, displayTrees, charTrees, f)
                                                                                            else updateAndFinalizePostOrderSoftWired (Just rootIndex) rootIndex inPhyloGraph
            
            -- purge double edges from display and character graphs
            -- this should not be happening--issue with postorder network resolutions data
            (inBlockGraphV, inBlockCharGraphVV) = if inGraphType == Tree then (inBlockGraphV', inBlockCharGraphVV')
                                                  else (fmap (fmap LG.removeDuplicateEdges) inBlockGraphV', fmap (fmap LG.removeDuplicateEdges) inBlockCharGraphVV')

            -- reroot block character trees
            -- not sure if should be parallelized `using` PU.myParListChunkRDS
            (newBlockDisplayTreeVect, newBlockCharGraphVV, blockCostV) = unzip3 (zipWith3  (rerootBlockCharTrees rootIndex) (V.toList $ fmap head inBlockGraphV) (V.toList inBlockCharGraphVV) (V.toList charInfoVV) `using` PU.myParListChunkRDS) 
            -- This is slower than myParListChunkRDS
            -- (newBlockDisplayTreeVect, newBlockCharGraphVV, blockCostV) = unzip3 (PU.seqParMap rdeepseq (rerootBlockCharTrees' rootIndex) $ zip3 (V.toList $ fmap head inBlockGraphV) (V.toList inBlockCharGraphVV) (V.toList charInfoVV))
            
            newCononicalGraph = backPortBlockTreeNodesToCanonicalGraph inDecGraph (V.fromList newBlockDisplayTreeVect)
        in
        -- trace ("GDBRS:" ++ "Dec graph:" ++ (LG.prettyIndices decGraph) ++ "\n" ++ (concatMap (++ "\nNew: ") $ fmap show $ fmap LG.getDuplicateEdges $ V.toList newBlockDisplayTreeVect))
        (inSimpleGraph, sum blockCostV, newCononicalGraph, V.fromList $ fmap (:[]) newBlockDisplayTreeVect, V.fromList newBlockCharGraphVV, charInfoVV)
   

-- | rerootBlockCharTrees' wrapper for fmap/seqparmap
rerootBlockCharTrees' :: LG.Node -> (DecoratedGraph, V.Vector DecoratedGraph, V.Vector CharInfo) -> (DecoratedGraph, V.Vector DecoratedGraph, VertexCost)
rerootBlockCharTrees' rootIndex (blockDisplayTree, charTreeVect, charInfoVect) = rerootBlockCharTrees rootIndex blockDisplayTree charTreeVect charInfoVect

-- | rerootBlockCharTrees reroots all character trees (via fmap) in block returns best block char trees and costs
-- with best character tree node assignment back ported to display tree
rerootBlockCharTrees ::LG.Node -> DecoratedGraph -> V.Vector DecoratedGraph -> V.Vector CharInfo -> (DecoratedGraph, V.Vector DecoratedGraph, VertexCost)
rerootBlockCharTrees rootIndex blockDisplayTree charTreeVect charInfoVect =
    if V.null charTreeVect then error "Empty tree vector in rerootBlockCharTrees"
    else 
        let -- next edges (to vertex in list) to perform rerroting
            -- progresses recursivey over adjacent edges to minimize node reoptimization
            -- since initially all same graph can get initial reroot nodes from display tree
            childrenOfRoot = LG.descendants blockDisplayTree rootIndex
            grandChildrenOfRoot = concatMap (LG.descendants blockDisplayTree) childrenOfRoot

            (rerootedCharTreeVect, rerootedCostVect) = unzip (zipWith (getCharTreeBestRoot rootIndex grandChildrenOfRoot) (V.toList charTreeVect) (V.toList charInfoVect) `using` PU.myParListChunkRDS)

            -- this is a little slower than myParListChunkRDS
            -- (rerootedCharTreeVect, rerootedCostVect) = unzip (PU.seqParMap rdeepseq (getCharTreeBestRoot' rootIndex grandChildrenOfRoot) (zip (V.toList charTreeVect) (V.toList charInfoVect)))

            updateBlockDisplayTree = backPortCharTreeNodesToBlockTree blockDisplayTree (V.fromList rerootedCharTreeVect)
        in
        (updateBlockDisplayTree, V.fromList rerootedCharTreeVect, sum rerootedCostVect)

-- | getCharTreeBestRoot' wrapper for use with fmap/seqParMap 
getCharTreeBestRoot' :: LG.Node -> [LG.Node] -> (DecoratedGraph, CharInfo) -> (DecoratedGraph, VertexCost)
getCharTreeBestRoot' rootIndex nodesToRoot (inCharacterGraph, charInfo) = getCharTreeBestRoot rootIndex nodesToRoot inCharacterGraph charInfo

-- | getCharTreeBestRoot takes the root index, a character tree (from a block) and its character info
--- and prerforms the rerootings of that character tree to get the best reroot cost and preliminary assignments
getCharTreeBestRoot :: LG.Node -> [LG.Node] -> DecoratedGraph -> CharInfo -> (DecoratedGraph, VertexCost)
getCharTreeBestRoot rootIndex nodesToRoot inCharacterGraph charInfo =
    -- if prealigned should be rerooted?
    let (bestRootCharGraph, bestRootCost) = if (charType charInfo `notElem` sequenceCharacterTypes) then (inCharacterGraph, (subGraphCost . snd) $ LG.labelNode inCharacterGraph rootIndex)
                                            else rerootCharacterTree rootIndex nodesToRoot charInfo inCharacterGraph 
    in
    (bestRootCharGraph, bestRootCost)

-- | rerootCharacterTree wrapper around rerootCharacterTree' with cleaner interface for "best" results
rerootCharacterTree :: LG.Node -> [LG.Node] ->  CharInfo -> DecoratedGraph -> (DecoratedGraph, VertexCost)
rerootCharacterTree  rootIndex nodesToRoot charInfo inCharacterGraph = 
    rerootCharacterTree' rootIndex nodesToRoot charInfo ((subGraphCost . snd) $ LG.labelNode inCharacterGraph rootIndex) inCharacterGraph inCharacterGraph

-- | rerootCharacterTree' takes a character tree and root index and returns best rooted character tree and cost
-- this is recursive taking best cost to save on memory over an fmap and minimum
-- since does reroot stuff over character trees--that component is less efficient
-- root index always same--just edges conenct to change with rerooting
-- graph is prgressively rerooted to be efficient
rerootCharacterTree' ::LG.Node -> [LG.Node] -> CharInfo -> VertexCost -> DecoratedGraph -> DecoratedGraph -> (DecoratedGraph, VertexCost)
rerootCharacterTree' rootIndex nodesToRoot charInfo bestCost bestGraph inGraph =
    if null nodesToRoot then (bestGraph, bestCost)
    else
        let firstRerootIndex = head nodesToRoot
            nextReroots = (LG.descendants inGraph firstRerootIndex) ++ tail nodesToRoot
            newGraph = rerootAndDiagnoseTree rootIndex firstRerootIndex charInfo inGraph
            newGraphCost = ((subGraphCost . snd) $ LG.labelNode newGraph rootIndex)    
            (bestGraph', bestCost') = if newGraphCost < bestCost then (newGraph, newGraphCost)
                                      else (bestGraph, bestCost)
        in
        -- if LG.isEmpty newGraph then rerootCharacterTree' rootIndex (tail nodesToRoot) charInfo bestCost bestGraph inGraph 
        -- trace ("RRCT:" ++ (show (rootIndex, firstRerootIndex, bestCost, newGraphCost))) 
        --else 
        rerootCharacterTree' rootIndex nextReroots charInfo bestCost' bestGraph' newGraph   
    
-- | rerootAndDiagnoseTree takes tree and reroots and reoptimizes nodes
rerootAndDiagnoseTree :: LG.Node -> LG.Node ->  CharInfo -> DecoratedGraph -> DecoratedGraph
rerootAndDiagnoseTree rootIndex newRerootIndex charInfo inGraph =
    let reRootGraph = LG.rerootDisplayTree rootIndex newRerootIndex inGraph
        (nodesToOptimize, _) = LG.pathToRoot inGraph (LG.labelNode inGraph newRerootIndex)
        reOptimizedGraph = reOptimizeCharacterNodes charInfo reRootGraph nodesToOptimize
    in
    if LG.isEmpty reRootGraph then inGraph
    else reOptimizedGraph


-- | reOptimizeCharacterNodes takes a decorated graph and a list of nodes and reoptimizes (relabels)
-- them based on children in input graph
-- simple recursive since each node depends on children
-- check for out-degree 1 since can be resolved form diplay trees 
reOptimizeCharacterNodes :: CharInfo -> DecoratedGraph -> [LG.LNode VertexInfo] -> DecoratedGraph
reOptimizeCharacterNodes charInfo inGraph oldNodeList =
  -- trace ("RON:" ++ (show $ fmap fst oldNodeList)) (
  if null oldNodeList then inGraph
  else
    let curNode@(curNodeIndex, curNodeLabel) = head oldNodeList
        nodeChildren = LG.descendants inGraph curNodeIndex  -- should be 1 or 2, not zero since all leaves already in graph
        foundCurChildern = filter (`elem` nodeChildren) $ fmap fst (tail oldNodeList)
    in
    {-These are checks that were in for network code--should be unncesesary for charactaer trees
    -- make sure that nodes are optimized in correct order so that nodes are only reoptimized using updated children
    -- this should really not have to happen--order should be determined a priori
    -}
    --if LG.isLeaf inGraph curNodeIndex then trace ("Should not be a leaf in reoptimize nodes: " ++ (show curNodeIndex) ++ " children " ++ (show nodeChildren) ++ "\nGraph:\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inGraph)) inGraph
    --else 

    if not $ null foundCurChildern then
      -- trace ("Current node " ++ (show curNodeIndex) ++ " has children " ++ (show nodeChildren) ++ " in optimize list (optimization order error)" ++ (show $ fmap fst $ tail oldNodeList))
      reOptimizeCharacterNodes charInfo inGraph (tail oldNodeList ++ [curNode])

    -- somehow root before others -- remove if not needed after debug
    --else if LG.isRoot inGraph curNodeIndex && length oldNodeList > 1 then
    --    error ("Root first :" ++ (show $ fmap fst oldNodeList) ++ "RC " ++ show (LG.descendants inGraph curNodeIndex)) -- ++ "\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inGraph))
        --reOptimizeNodes localGraphType charInfoVectVect inGraph ((tail oldNodeList) ++ [curNode])
    else if length nodeChildren > 2 then error ("Node has >2 children: " ++ (show nodeChildren))
    else
    
        -- trace ("RON: " ++ (show curNodeIndex) ++ " children " ++ (show nodeChildren)) (
        let leftChild = head nodeChildren
            rightChild = last nodeChildren
            
            -- this ensures that left/right choices are based on leaf BV for consistency and label invariance
            (leftChildLabel, rightChildLabel) = U.leftRightChildLabelBV (fromJust $ LG.lab inGraph leftChild, fromJust $ LG.lab inGraph rightChild)
            (newVertexData, newVertexCost) = M.median2Single False ((V.head . V.head . vertData) leftChildLabel) ((V.head . V.head . vertData)  rightChildLabel) charInfo

        in
       
        let (newCost, newBVLabel, newVertData, newSubGraphCost) = if length nodeChildren < 2 then (0, bvLabel leftChildLabel, vertData leftChildLabel, subGraphCost leftChildLabel)
                                                                   else (newVertexCost, bvLabel leftChildLabel .|. bvLabel rightChildLabel, V.singleton (V.singleton newVertexData), subGraphCost leftChildLabel + subGraphCost rightChildLabel + newCost)
            newVertexLabel = VertexInfo {  index = curNodeIndex
                                         -- this bit labelling incorect for outdegree = 1, need to prepend bits
                                         , bvLabel = newBVLabel
                                         , parents = V.fromList $ LG.parents inGraph curNodeIndex
                                         , children = V.fromList nodeChildren
                                         , nodeType = nodeType curNodeLabel
                                         , vertName = vertName curNodeLabel
                                         , vertexResolutionData = mempty
                                         , vertData = newVertData
                                         , vertexCost = newCost
                                         , subGraphCost = newSubGraphCost
                                         }

            -- this to add back edges deleted with nodes (undocumented but sensible in fgl)
            replacementEdges = LG.inn inGraph curNodeIndex ++ LG.out inGraph curNodeIndex
            newGraph = LG.insEdges replacementEdges $ LG.insNode (curNodeIndex, newVertexLabel) $ LG.delNode curNodeIndex inGraph
         in
         --trace ("New vertexCost " ++ show newCost) --  ++ " lcn " ++ (show (vertData leftChildLabel, vertData rightChildLabel, vertData curnodeLabel)))
         reOptimizeCharacterNodes charInfo newGraph (tail oldNodeList)


-- | backPortCharTreeNodesToBlockTree assigned nodes states (labels) of character trees to block doisplay Tree
-- updates vertData, vertexCost, and subGraphCost for each .  Subgraph cost queationable since relieds on rooting
backPortCharTreeNodesToBlockTree :: DecoratedGraph -> V.Vector DecoratedGraph -> DecoratedGraph
backPortCharTreeNodesToBlockTree blockDisplayTree rerootedCharTreeVect =
    let blockDisplayNodes = LG.labNodes blockDisplayTree
        blockDisplayEdges = LG.labEdges blockDisplayTree
        
        -- vector (characters) of vector (nodes) of labels
        charTreeLabelsVV = fmap V.fromList $ fmap (fmap snd) $ fmap LG.labNodes rerootedCharTreeVect

        -- for each node index extract (head head) vertdata, vertexCost and subgraphcost
        (vertDataVV, vertCostVV, subGraphCostVV) = V.unzip3 $ fmap (extractTripleVect charTreeLabelsVV) (V.fromList [0..(length blockDisplayNodes - 1)])

        -- update labels for block data nodes 
        updatedDisplayNodes = V.zipWith4 updateNodes (V.fromList blockDisplayNodes) vertDataVV vertCostVV subGraphCostVV

    in
    LG.mkGraph (V.toList updatedDisplayNodes) blockDisplayEdges

-- | extractTripleVect takes a vector of vector character tree labels and a node index and
-- retuns a triple of data (vertData, VertCost, and subgraphCost) from a given node index in all labels
extractTripleVect :: V.Vector (V.Vector VertexInfo) -> Int -> (V.Vector CharacterData, V.Vector VertexCost, V.Vector VertexCost)
extractTripleVect inLabelVV charIndex =
    let nodeLabelV = fmap (V.! charIndex) inLabelVV
        vertDataV = fmap vertData nodeLabelV
        vertCostV = fmap vertexCost nodeLabelV
        subGraphCostV = fmap subGraphCost nodeLabelV
    in
    (fmap (V.head . V.head) vertDataV, vertCostV, subGraphCostV)


-- | updateNodes takes vectors of labelled nodes and updates vertData, VerTCost, and subgraphCost fields
updateNodes ::LG.LNode VertexInfo -> V.Vector CharacterData -> V.Vector VertexCost -> V.Vector VertexCost -> LG.LNode VertexInfo
updateNodes (inIndex, inLabel) charDataV vertexCostV subGraphCostV =
    let newVertCost = V.sum vertexCostV
        newSubGraphCost = V.sum subGraphCostV
        newLabel = inLabel {vertData = V.singleton charDataV, vertexCost = newVertCost, subGraphCost = newSubGraphCost}
    in
    (inIndex, newLabel)

-- | backPortBlockTreeNodesToCanonicalGraph takes block display trees (updated presumably) and ports the block tree node 
-- labels to the cononical Graph
-- very similar to backPortCharTreeNodesToBlockTree but character vector is no singleton
backPortBlockTreeNodesToCanonicalGraph :: DecoratedGraph -> V.Vector DecoratedGraph -> DecoratedGraph
backPortBlockTreeNodesToCanonicalGraph inCanonicalGraph blockTreeVect = 
    let canonicalDisplayNodes = LG.labNodes inCanonicalGraph
        canonicalDisplayEdges = LG.labEdges inCanonicalGraph

        -- vector (characters) of vector (nodes) of labels
        blockTreeLabelsVV = fmap V.fromList $ fmap (fmap snd) $ fmap LG.labNodes blockTreeVect

        -- for each node index extract (head head) vertdata, vertexCost and subgraphcost
        (vertDataVV, vertCostVV, subGraphCostVV) = V.unzip3 $ fmap (extractTripleVectBlock blockTreeLabelsVV) (V.fromList [0..(length canonicalDisplayNodes - 1)])

        -- update labels for canonical data nodes 
        updatedCanonicalNodes = V.zipWith4 updateNodesBlock (V.fromList canonicalDisplayNodes) vertDataVV vertCostVV subGraphCostVV

    in
    LG.mkGraph (V.toList updatedCanonicalNodes) canonicalDisplayEdges

-- | updateNodesBlock takes vectors of labelled nodes and updates vertData, VerTCost, and subgraphCost fields
updateNodesBlock ::LG.LNode VertexInfo -> V.Vector (V.Vector CharacterData) -> V.Vector VertexCost -> V.Vector VertexCost -> LG.LNode VertexInfo
updateNodesBlock (inIndex, inLabel) charDataVV vertexCostV subGraphCostV =
    let newVertCost = V.sum vertexCostV
        newSubGraphCost = V.sum subGraphCostV
        newLabel = inLabel {vertData = charDataVV, vertexCost = newVertCost, subGraphCost = newSubGraphCost}
    in
    (inIndex, newLabel)

-- | extractTripleVectBlock takes a vector of vector block tree labels and a node index and
-- retuns a triple of data (vertData, VertCost, and subgraphCost) from a given node index in all labels
extractTripleVectBlock :: V.Vector (V.Vector VertexInfo) -> Int -> (V.Vector (V.Vector CharacterData), V.Vector VertexCost, V.Vector VertexCost)
extractTripleVectBlock inLabelVV charIndex =
    let nodeLabelV = fmap (V.! charIndex) inLabelVV
        vertDataV = fmap vertData nodeLabelV
        vertCostV = fmap vertexCost nodeLabelV
        subGraphCostV = fmap subGraphCost nodeLabelV
    in
    (fmap V.head vertDataV, vertCostV, subGraphCostV)

-- | divideDecoratedGraphByBlockAndCharacterTree takes a DecoratedGraph with (potentially) multiple blocks
-- and (potentially) multiple character per block and creates a Vector of Vector of Decorated Graphs
-- over blocks and characters with the same graph, but only a single block and character for each graph
-- this to be used to create the "best" cost over alternate graph traversals
-- vertexCost and subGraphCost will be taken from characterData localcost/localcostVect and globalCost
divideDecoratedGraphByBlockAndCharacterTree :: DecoratedGraph -> (V.Vector [DecoratedGraph], V.Vector (V.Vector DecoratedGraph))
divideDecoratedGraphByBlockAndCharacterTree inGraph =
  if LG.isEmpty inGraph then (V.empty, V.empty)
  else
    let numBlocks = V.length $ vertData $ snd $ head $ LG.labNodes inGraph
        blockGraphList = fmap (pullBlock inGraph) [0.. (numBlocks - 1)]
        characterGraphList = fmap makeCharacterGraph blockGraphList
    in
    -- trace ("DDGBCT: Blocks " ++ (show numBlocks) ++ " Characters " ++ (show $ fmap length $ vertData $ snd $ head $ LG.labNodes inGraph) ++ "\n" ++ (show characterGraphList))
    (V.fromList (fmap (:[]) blockGraphList), V.fromList characterGraphList)

-- | pullBlocks take a DecoratedGraph and creates a newDecorated graph with
-- only data from the input block index
pullBlock :: DecoratedGraph -> Int -> DecoratedGraph
pullBlock inGraph blockIndex =
  if LG.isEmpty inGraph then LG.empty
  else
    let (inNodeIndexList, inNodeLabelList) = unzip $ LG.labNodes inGraph
        blockNodeLabelList = fmap (makeBlockNodeLabels blockIndex) inNodeLabelList
    in
    LG.mkGraph (zip inNodeIndexList blockNodeLabelList) (LG.labEdges inGraph)

-- | makeBlockNodeLabels takes a block index and an orginal nodel label
-- and creates a new list of a singleton block from the input block index
makeBlockNodeLabels :: Int -> VertexInfo -> VertexInfo
makeBlockNodeLabels blockIndex inVertexInfo =
  let newVertexData = vertData inVertexInfo V.! blockIndex
      newVertexCost = V.sum $ fmap localCost newVertexData
      newsubGraphCost = V.sum $ fmap globalCost newVertexData
  in
  -- trace ("MBD " ++ (show $ length newVertexData) ++ " from " ++ (show $ length (vertData inVertexInfo)))
  inVertexInfo { vertData     = V.singleton newVertexData
               , vertexCost   = newVertexCost
               , subGraphCost = newsubGraphCost
               }

-- | updateAndFinalizePostOrderSoftWired performs the pre-order traceback on the resolutions of a softwired graph to create the correct vertex states,
-- ports the post order assignments to the display trees, and creates the character trees from the block trees
updateAndFinalizePostOrderSoftWired :: Maybe Int -> Int -> PhylogeneticGraph -> PhylogeneticGraph
updateAndFinalizePostOrderSoftWired startVertex rootIndex inGraph =
    if LG.isEmpty $ thd6 inGraph then inGraph
    else
        let -- this step is important--not repetetive since pull out teh resolution data for each block
            outgroupRootLabel =  fromJust $ LG.lab (thd6 inGraph) rootIndex
            (displayGraphVL, lDisplayCost) = extractDisplayTrees startVertex True (vertexResolutionData outgroupRootLabel)
            
            -- traceback on resolutions
            newGraph = softWiredPostOrderTraceBack rootIndex (thd6 inGraph)

            -- propagate updated post-order assignments to display trees, which updates display tree and cost ealier
            displayGraphVL' = V.zipWith (assignPostOrderToDisplayTree (fmap (vertData . snd) (LG.labNodes newGraph) )) displayGraphVL (V.fromList [0..(V.length displayGraphVL - 1)])

            -- create new, fully  updated post-order graph
            finalPreOrderGraph = (fst6 inGraph, lDisplayCost, newGraph, displayGraphVL', divideDecoratedGraphByBlockAndCharacterSoftWired displayGraphVL', six6 inGraph)
        in
        -- trace ("UAFPS: after extraction: " ++ (show $ fmap LG.getDuplicateEdges $ fmap head $ V.toList displayGraphVL) ++ " after assignment: " ++ (show $ fmap LG.getDuplicateEdges $ fmap head $ V.toList displayGraphVL'))
        -- trace ("UFPOSW: " ++ (show $ fmap length displayGraphVL) ++ " " ++ (show $ fmap length displayGraphVL') ++ " " ++ (show $ fmap V.length $ fft6 finalPreOrderGraph))
        finalPreOrderGraph

-- | divideDecoratedGraphByBlockAndCharacterSoftWired takes a Vector of a list of DecoratedGraph
-- containing a list of decorated trees that are the display trees for that block
-- with (potentially) multiple blocks
-- and (potentially) multiple character per block and creates a Vector of Vector of Decorated Graphs
-- over blocks and characters with the block diplay graph, but only a single block and character for each graph
-- this to be used to create the "best" cost over alternate graph traversals
-- vertexCost and subGraphCost will be taken from characterData localcost/localcostVect and globalCost
-- for this assignment purpose for pre-order a single (head) member of list is used to create the
-- character graphs
divideDecoratedGraphByBlockAndCharacterSoftWired :: V.Vector [DecoratedGraph] -> V.Vector (V.Vector DecoratedGraph)
divideDecoratedGraphByBlockAndCharacterSoftWired inGraphVL =
  if V.null inGraphVL then mempty
  else
    let blockGraphList = fmap head inGraphVL
        characterGraphList = fmap makeCharacterGraph blockGraphList
    in
    characterGraphList

-- | postDecorateSoftWired' wrapper for postDecorateSoftWired with args in differnt order for mapping
postDecorateSoftWired' :: GlobalSettings -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.Node -> LG.Node -> SimpleGraph -> PhylogeneticGraph
postDecorateSoftWired' inGS curDecGraph blockCharInfo rootIndex curNode simpleGraph = postDecorateSoftWired inGS simpleGraph curDecGraph blockCharInfo rootIndex curNode

-- | postDecorateSoftWired begins at start index (usually root, but could be a subtree) and moves preorder till children are labelled
-- and then recurses to root postorder labelling vertices and edges as it goes
-- this for a single root
postDecorateSoftWired :: GlobalSettings -> SimpleGraph -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.Node -> LG.Node -> PhylogeneticGraph
postDecorateSoftWired inGS simpleGraph curDecGraph blockCharInfo rootIndex curNode =
    -- if node in current decortated graph then nothing to do and return
    if LG.gelem curNode curDecGraph then
        let nodeLabel = LG.lab curDecGraph curNode
        in
        if isNothing nodeLabel then error ("Null label for node " ++ show curNode)
        else (simpleGraph, subGraphCost (fromJust nodeLabel), curDecGraph, mempty, mempty, blockCharInfo)

    else
        -- get postodre assignmens of children
        -- checks for single child of node
        -- result is single graph afer left and right child traversals
        -- trace ("PDSW making node " ++ show curNode ++ " in\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph curDecGraph)) (
        let nodeChildren = LG.descendants simpleGraph curNode  -- should be 1 or 2, not zero since all leaves already in graph
            leftChild = head nodeChildren
            rightChild = last nodeChildren
            leftChildTree = postDecorateSoftWired inGS simpleGraph curDecGraph blockCharInfo rootIndex leftChild
            rightLeftChildTree = if length nodeChildren == 2 then postDecorateSoftWired inGS simpleGraph (thd6 leftChildTree) blockCharInfo rootIndex rightChild
                                 else leftChildTree
        in
        -- Checks on children
        if length nodeChildren > 2 then error ("Graph not dichotomous in postDecorateSoftWired node " ++ show curNode ++ "\n" ++ LG.prettify simpleGraph)
        else if null nodeChildren then error ("Leaf not in graph in postDecorateSoftWired node " ++ show curNode ++ "\n" ++ LG.prettify simpleGraph)

        else
            -- make node from child block resolutions
            -- child resolutin made ealeri in post roder pass
            let newSubTree = thd6 rightLeftChildTree
            in

            -- single child of node (can certinly happen with soft-wired networks
            if length nodeChildren == 1 then
                -- trace ("Outdegree 1: " ++ (show curNode) ++ " " ++ (show $ GO.getNodeType simpleGraph curNode) ++ " Child: " ++ (show nodeChildren)) (
                let (newGraph, _, _, _, _) = getOutDegree1VertexAndGraph curNode (fromJust $ LG.lab newSubTree leftChild) simpleGraph nodeChildren newSubTree
                in
                (simpleGraph, 0, newGraph, mempty, mempty, blockCharInfo)

            -- check for error condition
            else if length nodeChildren > 2 then error ("Node has >2 children: " ++ (show nodeChildren))
    
            -- 2 children
            else
                -- trace ("Outdegree 2: " ++ (show curNode) ++ " " ++ (show $ GO.getNodeType simpleGraph curNode) ++ " Children: " ++ (show nodeChildren)) (
                -- need to create new resolutions and add to existing sets
                let -- this ensures that left/right choices are based on leaf BV for consistency and label invariance
                    -- larger bitvector is Right, smaller or equal Left
                    ((leftChild', leftChildLabel), (rightChild', rightChildLabel)) = U.leftRightChildLabelBVNode ((leftChild, fromJust $ LG.lab newSubTree leftChild), (rightChild, fromJust $ LG.lab newSubTree rightChild))

                    -- create resolution caches for blocks
                    leftChildNodeType  = nodeType leftChildLabel
                    rightChildNodeType = nodeType rightChildLabel
                    resolutionBlockVL = V.zipWith3 (createBlockResolutions (compressResolutions inGS) curNode leftChild' rightChild' leftChildNodeType rightChildNodeType (GO.getNodeType simpleGraph curNode)) (vertexResolutionData leftChildLabel) (vertexResolutionData rightChildLabel) blockCharInfo

                    -- create canonical Decorated Graph vertex
                    -- 0 cost becasue can't know cosrt until hit root and get best valid resolutions
                    newVertexLabel = VertexInfo { index = curNode
                                                , bvLabel = bvLabel leftChildLabel .|. bvLabel rightChildLabel
                                                , parents = V.fromList $ LG.parents simpleGraph curNode
                                                , children = V.fromList nodeChildren
                                                , nodeType = GO.getNodeType simpleGraph curNode
                                                , vertName = T.pack $ "HTU" ++ show curNode
                                                , vertData = mempty --empty because of resolution data
                                                , vertexResolutionData = resolutionBlockVL
                                                , vertexCost = 0.0 --newCost
                                                , subGraphCost = 0.0 -- (subGraphCost leftChildLabel) + (subGraphCost rightChildLabel) + newCost
                                                }

                    leftEdgeType
                      | leftChildNodeType == NetworkNode = NetworkEdge
                      | leftChildNodeType == LeafNode = PendantEdge
                      | otherwise = TreeEdge
                    rightEdgeType
                      | rightChildNodeType == NetworkNode = NetworkEdge
                      | rightChildNodeType == LeafNode = PendantEdge
                      | otherwise = TreeEdge

                    edgeLable = EdgeInfo { minLength = 0.0
                                         , maxLength = 0.0
                                         , midRangeLength = 0.0
                                         , edgeType = TreeEdge
                                         }

                    leftEdge =  (curNode, leftChild', edgeLable {edgeType = leftEdgeType})
                    rightEdge = (curNode, rightChild', edgeLable {edgeType = rightEdgeType})
                    newGraph =  LG.insEdges [leftEdge, rightEdge] $ LG.insNode (curNode, newVertexLabel) newSubTree

                    (displayGraphVL, lDisplayCost) = if curNode == rootIndex then 
                                                        let (displayG, displayCost') = extractDisplayTrees (Just curNode) True resolutionBlockVL
                                                        in
                                                        (displayG, displayCost')
                                                        -- (fmap (fmap LG.removeDuplicateEdges) displayG, displayCost')
                                                     else (mempty, 0.0)

                in
                -- trace ("PDWS: " ++ (show $ fmap LG.toEdge [leftEdge, rightEdge]) ++ " has edges " ++ (show $ fmap (LG.hasEdge newSubTree) $ fmap LG.toEdge [leftEdge, rightEdge]) ++ " dupes: " ++ (show $ fmap LG.getDuplicateEdges $ V.head $ fst (extractDisplayTrees (Just curNode) True resolutionBlockVL)) ++ " Resolutions " ++ (show $ fmap (fmap U.hasResolutionDuplicateEdges) resolutionBlockVL))
    
                (simpleGraph, lDisplayCost, newGraph, displayGraphVL, mempty, blockCharInfo)

-- | assignPostOrderToDisplayTree takes the post-order (preliminary) block data from canonical Decorated graph
-- to Block display graphs (through list of more than one for a block)
-- input is canonical Decorated Graph and a pair containing the display tree and its Block index
-- this could be integrated into preorder traceback to remove the extra pass
-- do this here is a bit simpler
assignPostOrderToDisplayTree :: [VertexBlockData] -> [DecoratedGraph] -> Int -> [DecoratedGraph]
assignPostOrderToDisplayTree canonicalVertexBlockData displayTreeList displayTreeIndex =
    if null canonicalVertexBlockData then []
    else
        let
            updatedDispayTreeList = fmap (assignVertexBlockData canonicalVertexBlockData displayTreeIndex) displayTreeList
        in
        -- trace ("Index: " ++ (show displayTreeIndex) ++ " Blocks : " ++ (show $ V.length $ head canonicalVertexBlockData) ++ " Nodes: "  ++ (show $ length canonicalVertexBlockData))
        updatedDispayTreeList

-- | assignVertexBlockData assigns the block data of input index and assigns to all nodes of a tree
assignVertexBlockData :: [VertexBlockData] -> Int -> DecoratedGraph -> DecoratedGraph
assignVertexBlockData nodeDataList blockIndex inGraph =
    if null nodeDataList || LG.isEmpty inGraph then LG.empty
    else
        -- trace ("In Index: " ++ (show blockIndex) ++ " BL: " ++ (show $ V.length $ head nodeDataList) ++ " lengths " ++ (show $ fmap V.length nodeDataList)) (
        let blockIndexDataList = fmap (V.! blockIndex) nodeDataList
            (displayNodeIndexList, displayNodeLabelList) = L.unzip $ LG.labNodes inGraph
            updatedNodeLabelList = zipWith updateVertData displayNodeLabelList blockIndexDataList
        in
        LG.mkGraph (zip displayNodeIndexList updatedNodeLabelList) (LG.labEdges inGraph)
        -- )
        where updateVertData inVertexInfo newVertData = inVertexInfo {vertData = V.singleton newVertData}

-- | softWiredPostOrderTraceBack takes resolution data and assigns correct resolution median to preliminary
-- data ssignments.  Proceeds via typical pre-order pass over display tree
softWiredPostOrderTraceBack :: Int -> DecoratedGraph -> DecoratedGraph
softWiredPostOrderTraceBack rootIndex inGraph  =
    if LG.isEmpty inGraph then LG.empty
    else
            -- get edges to remake graph after nodes are updated with preliminary states
        let rootLabel = fromJust $ LG.lab inGraph rootIndex
            inEdgeList = LG.labEdges inGraph

            -- root first--choose best resolutions--returns a list and takes head of that list of equal cost/traceback preliminary
            -- assignments.  Later could look at multiple.
            (rootVertData, rootSubGraphCost, rootResolutionCost, leftRightIndexVect) = getResolutionDataAndIndices rootLabel (V.singleton (Just (-1)))
            newRootLabel = rootLabel { vertData = rootVertData
                                     , vertexCost = rootResolutionCost
                                     , subGraphCost = rootSubGraphCost
                                     }
            newRootNode = (rootIndex, newRootLabel)

            -- get child/children of root
            rootChildren = LG.labDescendants inGraph newRootNode

            -- left / right to match post-order
            rootChildrenBV = fmap (bvLabel . snd) rootChildren
            rootChildrenIsLeft
              | length rootChildrenBV == 1 = [True]
              | head rootChildrenBV > (rootChildrenBV !! 1) = [False, True]
              | otherwise = [True, False]

            rootChildrenTuples = zip3 rootChildren (replicate (length rootChildren) leftRightIndexVect) rootChildrenIsLeft

            -- recurse to children with resolution index from parent
            softWiredUpdatedNodes = softWiredPrelimTraceback inGraph rootChildrenTuples [newRootNode]

        in
        LG.mkGraph softWiredUpdatedNodes inEdgeList

-- | softWiredPrelimTraceback takes a list of nodes to update (with left right index info) based
-- on resolution data, recurses with left right indices pre-order to leaves, keeping list of updated nodes
softWiredPrelimTraceback :: DecoratedGraph
                        -> [(LG.LNode VertexInfo, V.Vector (Maybe Int, Maybe Int), Bool)]
                        -> [LG.LNode VertexInfo]
                        -> [LG.LNode VertexInfo]
softWiredPrelimTraceback inGraph nodesToUpdate updatedNodes =
    if null nodesToUpdate then updatedNodes
    else
        let (firstNode, firstLeftRight, isLeft) = head nodesToUpdate

            -- ensure consistent left/right from post-order
            resolutionIndexVect = if isLeft then fmap fst firstLeftRight
                                  else fmap snd firstLeftRight

            -- get resolution info
            (newBlockData, newSubGraphCost, newVertexCost, childLeftRightIndexVect) = getResolutionDataAndIndices (snd firstNode) resolutionIndexVect

            -- make new node
            newNodeLabel = (snd firstNode) { vertData = newBlockData
                                           , vertexCost = newVertexCost
                                           , subGraphCost = newSubGraphCost
                                           }

            newFirstNode = (fst firstNode, newNodeLabel)
        in

        -- not really necessary, but avoids the graph query operations
        if nodeType (snd firstNode) == LeafNode then
            softWiredPrelimTraceback inGraph (tail nodesToUpdate) (newFirstNode : updatedNodes)

        -- checks if network node (and its children) has been visited already
        else if (nodeType (snd firstNode) == NetworkNode) && isJust (L.find ((== fst firstNode). fst) updatedNodes) then
            softWiredPrelimTraceback inGraph (tail nodesToUpdate) updatedNodes

        else
            let -- get children
                firstChildren = LG.labDescendants inGraph firstNode
                firstChildrenBV = fmap (bvLabel . snd) firstChildren
                firstChildrenIsLeft
                  | length firstChildrenBV == 1 = [True] -- seems dumb but may assume lenght == 2 later
                  | head firstChildrenBV > (firstChildrenBV !! 1) = [False, True]
                  | otherwise = [True, False]

                childrenTuples = zip3 firstChildren (replicate (length firstChildren) childLeftRightIndexVect) firstChildrenIsLeft


        in
        {-
        if V.head resolutionIndexVect == Nothing then error ("'Nothing' index in softWiredPrelimTraceback " ++ (show resolutionIndexVect)
            ++ " Node " ++ (show $ fst firstNode) ++ " isLeft?: " ++ (show isLeft)  ++ " " ++ (show firstLeftRight))
        else
        -}
            softWiredPrelimTraceback inGraph (childrenTuples ++ tail nodesToUpdate) (newFirstNode : updatedNodes)


-- | getResolutionDataAndIndices takes a vertex label (VertexInfo) and returns the resolution data corresponding to
-- the index taken from its child resolution (data, subgraph cost, local resolution cost, left/right pairs).
-- Index = (-1) denotes that it is a root label and in that case
-- the best (lowest cost) resolutions are returned
getResolutionDataAndIndices :: VertexInfo -> V.Vector (Maybe Int) -> (VertexBlockData, VertexCost, VertexCost, V.Vector (Maybe Int, Maybe Int))
getResolutionDataAndIndices nodeLabel parentResolutionIndexVect =

    -- should not happen
    --mtrace ("NL " ++ (show $ index nodeLabel) ++ " PRIL: " ++ " length " ++ (show $ V.length parentResolutionIndexVect) ++ " " ++ show parentResolutionIndexVect) (
    if nodeType nodeLabel == LeafNode then
        let leafVertData = fmap (displayData . V.head) (vertexResolutionData nodeLabel)
        in
        -- trace ("GRDI: Leaf " ++ (show $ fmap V.length (vertexResolutionData nodeLabel)))
        (leafVertData, 0, 0, V.singleton (Just 0, Just 0))

    -- root node--take lowest cost
    else if V.head parentResolutionIndexVect == Just (-1) then
        let rootBlockResolutionPair = getBestBlockResolution <$> vertexResolutionData nodeLabel
            (charDataVV, subGraphCostV, resCostV, leftRightIndexVect) = V.unzip4 rootBlockResolutionPair
        in
        -- trace ("GRDI: Root " ++ (show $ fmap show parentResolutionIndexVect) ++ " " ++ (show $ fmap V.length (vertexResolutionData nodeLabel)))
        (charDataVV, V.sum subGraphCostV, V.sum resCostV, leftRightIndexVect)

    -- non-root node--return the index resolution information
    else
        -- trace ("GRD Length:" ++ (show $ V.length $ vertexResolutionData nodeLabel) ++ " " ++ (show $ fmap fromJust parentResolutionIndexVect) ++ "\n" ++ (show $ vertexResolutionData nodeLabel)) (
        let parentIndexVect = fmap fromJust parentResolutionIndexVect

            -- get resolution data from node label
            resolutionData = vertexResolutionData nodeLabel

            -- get the correct (via index) resolution data for each block
            -- complex for network node since keeps left right sort of array, but only first element maters--this hack keeps things ok for
            -- tree-like traceback assignment
            -- not sure if first or last would be better--hopefully not matter, or arbitrary equal solutions
            resolutionsByBlockV = if nodeType nodeLabel == NetworkNode then
                                        -- trace ("-" ++ (show (V.length resolutionData, V.length parentIndexVect, V.head parentIndexVect)) ++ " " ++ (show $ parentIndexVect))
                                        -- V.zipWith (V.!) resolutionData (V.replicate (V.length parentIndexVect) (V.head parentIndexVect))
                                        V.zipWith (V.!) resolutionData (V.replicate (V.length parentIndexVect) 0)
                                  else 
                                    -- trace ("+" ++ (show (V.length resolutionData, V.length parentIndexVect)) ++ " " ++ (show $ parentIndexVect))
                                    V.zipWith (V.!) resolutionData parentIndexVect

            -- get other resolution info
            charDataVV = fmap displayData resolutionsByBlockV
            lSubGraphCost = V.sum $ fmap displayCost resolutionsByBlockV
            localResolutionCost = V.sum $ fmap resolutionCost resolutionsByBlockV

            -- only takes first left right pair--although others in there
            -- uses first for preliminary asiignment--but could be more
            -- if euqla cost display trees, may have multiple possible preliminary states
            leftRightIndexVect = fmap (head . childResolutions) resolutionsByBlockV
        in
        {-
        if V.null resolutionData then
           error ("Null resolution data in getResolutionDataAndIndices at node with label:" ++ (show nodeLabel))
        else
        -}
        -- trace ("GRDI: " ++ (show $ nodeType nodeLabel) ++ " " ++ (show $ fmap V.length resolutionData) ++ " " ++ (show parentIndexVect))
        (charDataVV, lSubGraphCost, localResolutionCost, leftRightIndexVect)
        -- )

-- | getBestBlockResolution takes vertexResolutionData and returns the best (lowest cost) resolution and associated data
-- for a single block of characters.  A single left right index pair is returns for child resolution sources. Could be multiple
-- from resolution compression (equal leaf set, equal median) and avaialble for potenitally implementd in future, multiple
-- preliminary assignments.  Only valid for root node with all leaves in graph
getBestBlockResolution :: ResolutionBlockData -> (V.Vector CharacterData, VertexCost, VertexCost, (Maybe Int, Maybe Int))
getBestBlockResolution inResBlockData =
    if V.null inResBlockData then (mempty, 0.0, 0.0, (Nothing, Nothing))
    else
        let -- makes sure all leaves in resolution
            -- displayPopList = fmap (complement . displayBVLabel) inResBlockData

            -- this for non full graph root--uses highest number of bits on--make sure all taxa in
            -- should be fine for a subgraph that has a single edge a base--may not be correct
            -- for a sunbgaph that has conections outside of its graph root.
            displayPopList' = fmap (popCount . displayBVLabel) inResBlockData
            maxPop = maximum displayPopList'



            -- subgraph cost
            displayCostList = fmap displayCost inResBlockData

            -- resolution local cost
            resolutionCostList = fmap resolutionCost inResBlockData

            -- takes only first resolution index pair
            childResolutionList = fmap (head . childResolutions) inResBlockData

            -- resolution medians
            displayDataList = fmap displayData inResBlockData

            {-
            -- take only those will all leaves in, then minimum cost
            quintVect = V.zip5 displayPopList displayCostList resolutionCostList childResolutionList displayDataList
            validVect = V.filter (BV.isZeroVector . fst5) quintVect
            validMinCost = V.minimum $ fmap snd5 validVect
            -}

            -- these for "best" this will largest leaf set
            quintVect' =  V.zip5 displayPopList' displayCostList resolutionCostList childResolutionList displayDataList
            validVect' = V.filter ((== maxPop) . fst5) quintVect'
            validMinCost' = V.minimum $ fmap snd5 validVect'

            -- ONly takes first of potentially multiple soliutions to begin traceback
            (_, displayCostV, resCostV, childIndexPairV, displayMedianV) = V.unzip5 $ V.filter  ((== validMinCost') . snd5) validVect'
        in
        if null validVect' then error "Null valid quad in getBestBlockResolution--perhaps not root node or forest component"
        else (V.head displayMedianV, V.head displayCostV, V.head resCostV, V.head childIndexPairV)

-- | makeLeafGraphSoftWired takes input data and creates a 'graph' of leaves with Vertex information
-- but with zero edges.  This 'graph' can be reused as a starting structure for graph construction
-- to avoid remaking of leaf vertices
-- includes leave resolution data
makeLeafGraphSoftWired :: ProcessedData -> DecoratedGraph
makeLeafGraphSoftWired (nameVect, bvNameVect, blocDataVect) =
    if V.null nameVect then error "Empty ProcessedData in makeLeafGraph"
    else
        let leafVertexList = V.toList $ V.map (makeLeafVertexSoftWired nameVect bvNameVect blocDataVect) (V.fromList [0.. V.length nameVect - 1])
        in
        LG.mkGraph leafVertexList []

-- | makeLeafVertexSoftWired makes a single unconnected vertex for a leaf in a Soft-wired graph
makeLeafVertexSoftWired :: V.Vector NameText -> V.Vector NameBV -> V.Vector BlockData -> Int -> LG.LNode VertexInfo
makeLeafVertexSoftWired nameVect bvNameVect inData localIndex =
    --trace ("Making leaf " ++ (show localIndex) ++ " Data " ++ (show $ length inData) ++ " " ++ (show $ fmap length $ fmap snd3 inData)) (
    let centralData = V.map snd3 inData
        thisData = V.map (V.! localIndex) centralData
        thisBVLabel = bvNameVect V.! localIndex
        thisResolutionData = makeLeafResolutionBlockData thisBVLabel ([(localIndex, minimalVertex)],[]) thisData
        minimalVertex = VertexInfo  { index = localIndex
                                    , bvLabel = thisBVLabel
                                    , parents = V.empty
                                    , children = V.empty
                                    , nodeType = LeafNode
                                    , vertName = nameVect V.! localIndex
                                    , vertData = mempty
                                    , vertexResolutionData = mempty
                                    , vertexCost = 0.0
                                    , subGraphCost = 0.0
                                    }
        newVertex = VertexInfo  { index = localIndex
                                , bvLabel = thisBVLabel
                                , parents = V.empty
                                , children = V.empty
                                , nodeType = LeafNode
                                , vertName =  nameVect V.! localIndex
                                , vertData = mempty
                                , vertexResolutionData = thisResolutionData
                                , vertexCost = 0.0
                                , subGraphCost = 0.0
                                }
        in
        --trace ("RD" ++ show $ thisResolutionData)
        (localIndex, newVertex)
        -- )


-- | makeLeafResolutionBlockData creates leaf resolution data from leav BVLabel, leave node, and data.
-- The return type is a vertor over character blocks each containing a list of potential resolutions (display trees) for that block
-- the resolutoins include subtree (display) for that resolution the bv l;abel for the node given that resolution, and the character data
-- in the block (Vector CharacterData) also given that resolution
-- thiis is repeaatted for each bloick in VertexBlockData
makeLeafResolutionBlockData :: NameBV -> ([LG.LNode VertexInfo], [LG.LEdge EdgeInfo]) -> VertexBlockData -> V.Vector ResolutionBlockData
makeLeafResolutionBlockData inBV inSubGraph inVertData =
    let defaultResolutionData = ResolutionData  { displaySubGraph = inSubGraph
                                                , displayBVLabel = inBV
                                                , displayData = mempty
                                                , childResolutions = [(Just 0, Just 0)]
                                                , resolutionCost = 0.0
                                                , displayCost = 0.0
                                                }

        blockIndexList = [0..(V.length inVertData - 1)]
        blockDataList = fmap (inVertData V.!) blockIndexList
        resolutionDataList = modifyDisplayData defaultResolutionData blockDataList []
        resolutionData =   V.fromList $ fmap (V.fromList . (:[])) resolutionDataList
    in
    resolutionData

-- | modifyDisplayData modifies displatData filed in ResolutionData
-- stas list doesn't change number of V.fromList calls
modifyDisplayData :: ResolutionData -> [V.Vector CharacterData] -> [ResolutionData] -> [ResolutionData]
modifyDisplayData resolutionTemplate characterDataVList curResolutionList =
    if null characterDataVList then reverse curResolutionList
    else
        let curBlockData = head characterDataVList
        in
        modifyDisplayData resolutionTemplate (tail characterDataVList) ((resolutionTemplate {displayData = curBlockData}) : curResolutionList)

-- | getOutDegree1VertexAndGraph makes parent node fomr single child for soft-wired resolutions
getOutDegree1VertexAndGraph :: (Show a, Show b)
                            => LG.Node
                            -> VertexInfo
                            -> LG.Gr a b
                            -> [LG.Node]
                            -> DecoratedGraph
                            -> (DecoratedGraph, Bool, VertexInfo, VertexCost, V.Vector [DecoratedGraph])
getOutDegree1VertexAndGraph curNode childLabel simpleGraph nodeChildren subTree =

    -- trace ("In out=1: " ++ (show curNode)) (
    let childResolutionData = vertexResolutionData childLabel

        curNodeResolutionData = addNodeAndEdgeToResolutionData newDisplayNode newLEdge childResolutionData

        newEdgeLabel = EdgeInfo { minLength = 0.0
                                 , maxLength = 0.0
                                 , midRangeLength = 0.0
                                 , edgeType = TreeEdge
                                 }
        newMinVertex = VertexInfo  { index = curNode
                                    , bvLabel = bvLabel childLabel
                                    , parents = V.fromList $ LG.parents simpleGraph curNode
                                    , children = V.fromList nodeChildren
                                    , nodeType = GO.getNodeType simpleGraph curNode
                                    , vertName = T.pack $ "HTU" ++ show curNode
                                    , vertData = mempty
                                    , vertexResolutionData = mempty
                                    , vertexCost = 0.0
                                    , subGraphCost = 0.0
                                    }

        newVertex  = VertexInfo { index = curNode
                                , bvLabel = bvLabel childLabel
                                , parents = V.fromList $ LG.parents simpleGraph curNode
                                , children = V.fromList nodeChildren
                                , nodeType = GO.getNodeType simpleGraph curNode
                                , vertName = T.pack $ "HTU" ++ show curNode
                                , vertData = mempty
                                , vertexResolutionData = curNodeResolutionData
                                , vertexCost = 0.0
                                , subGraphCost = subGraphCost childLabel
                                }

        newLEdge = (curNode, index childLabel, newEdgeLabel)
        newLNode = (curNode, newVertex)
        newDisplayNode = (curNode, newMinVertex)
        newGraph =  LG.insEdge newLEdge $ LG.insNode newLNode subTree

        (displayGraphVL, lDisplayCost) = if nodeType newVertex == RootNode then extractDisplayTrees Nothing True (vertexResolutionData childLabel)
                                        else (mempty, 0.0)


    in
    --trace ("NV1: " ++ show newVertex)
    --trace ("GOD1VG: " ++ (show $ LG.toEdge newLEdge) ++ " has edges " ++ (show $ LG.hasEdge subTree $  LG.toEdge newLEdge) ++ "Resolutions " ++ (show $ fmap (fmap U.hasResolutionDuplicateEdges) curNodeResolutionData))
    (newGraph, nodeType newVertex == RootNode, newVertex, lDisplayCost, displayGraphVL)
    -- (newGraph, False, newVertex, 0.0, mempty)
    -- )

-- | getOutDegree1VertexSoftWired returns new vertex only from single child for soft-wired resolutions
getOutDegree1VertexSoftWired :: (Show a, Show b)
                    => LG.Node
                    -> VertexInfo
                    -> LG.Gr a b
                    -> [LG.Node]
                    -> VertexInfo
getOutDegree1VertexSoftWired curNode childLabel simpleGraph nodeChildren =

    -- trace ("In out=1: " ++ (show curNode)) (
    let childResolutionData = vertexResolutionData childLabel


        newEdgeLabel = EdgeInfo { minLength = 0.0
                                 , maxLength = 0.0
                                 , midRangeLength = 0.0
                                 , edgeType = TreeEdge
                                 }
        newMinVertex = VertexInfo  { index = curNode
                                    , bvLabel = bvLabel childLabel
                                    , parents = V.fromList $ LG.parents simpleGraph curNode
                                    , children = V.fromList nodeChildren
                                    , nodeType = NetworkNode
                                    , vertName = T.pack $ "HTU" ++ show curNode
                                    , vertData = mempty
                                    , vertexResolutionData = mempty
                                    , vertexCost = 0.0
                                    , subGraphCost = 0.0
                                    }

        newDisplayNode = (curNode, newMinVertex)
        newLEdge = (curNode, index childLabel, newEdgeLabel)

        curNodeResolutionData = addNodeAndEdgeToResolutionData newDisplayNode newLEdge childResolutionData


        newVertexLabel  = VertexInfo { index = curNode
                                , bvLabel = bvLabel childLabel
                                , parents = V.fromList $ LG.parents simpleGraph curNode
                                , children = V.fromList nodeChildren
                                , nodeType = NetworkNode
                                , vertName = T.pack $ "HTU" ++ show curNode
                                , vertData = mempty
                                , vertexResolutionData = curNodeResolutionData
                                , vertexCost = 0.0
                                , subGraphCost = subGraphCost childLabel
                                }


    in
    newVertexLabel

-- | getOutDegree2VertexSoftWired returns new vertex only from two child nodes for soft-wired resolutions
getOutDegree2VertexSoftWired :: GlobalSettings
                             -> V.Vector (V.Vector CharInfo)
                             -> LG.Node
                             -> LG.LNode VertexInfo
                             -> LG.LNode VertexInfo
                             -> DecoratedGraph
                             -> VertexInfo
getOutDegree2VertexSoftWired inGS charInfoVectVect curNodeIndex leftChild@(leftChildIndex, _) rightChild@(rightChildIndex, _) inGraph =

    let -- this ensures that left/right choices are based on leaf BV for consistency and label invariance
        -- larger bitvector is Right, smaller or equal Left
        ((leftChild', leftChildLabel'), (rightChild', rightChildLabel')) = U.leftRightChildLabelBVNode (leftChild, rightChild)

        -- create resolution caches for blocks
        leftChildNodeType  = nodeType leftChildLabel'
        rightChildNodeType = nodeType rightChildLabel'
        resolutionBlockVL = V.zipWith3 (createBlockResolutions (compressResolutions inGS) curNodeIndex leftChild' rightChild' leftChildNodeType rightChildNodeType TreeNode) (vertexResolutionData leftChildLabel') (vertexResolutionData rightChildLabel') charInfoVectVect

        -- create canonical Decorated Graph vertex
        -- 0 cost becasue can't know cosrt until hit root and get best valid resolutions
        newVertexLabel = VertexInfo {  index = curNodeIndex
                                    , bvLabel = bvLabel leftChildLabel' .|. bvLabel rightChildLabel'
                                    , parents = V.fromList $ LG.parents inGraph curNodeIndex
                                    , children = V.fromList [leftChildIndex, rightChildIndex]
                                    , nodeType = TreeNode
                                    , vertName = T.pack $ "HTU" ++ show curNodeIndex
                                    , vertData = mempty --empty because of resolution data
                                    , vertexResolutionData = resolutionBlockVL
                                    , vertexCost = 0.0 --newCost
                                    , subGraphCost = 0.0 -- (subGraphCost leftChildLabel) + (subGraphCost rightChildLabel) + newCost
                                    }
        in
        newVertexLabel

-- | extractDisplayTrees takes resolutions and pulls out best cost (head for now) need to change type for multiple best
-- option for filter based on pop-count for root cost and complete display tree check
extractDisplayTrees :: Maybe Int -> Bool -> V.Vector ResolutionBlockData -> (V.Vector [DecoratedGraph], VertexCost)
extractDisplayTrees startVertex checkPopCount inRBDV =
    if V.null inRBDV then (V.empty, 0.0)
    else
        let (bestBlockDisplayResolutionList, costVect) = V.unzip $ fmap (getBestResolutionList startVertex checkPopCount) inRBDV
        in
        (bestBlockDisplayResolutionList, V.sum costVect)

-- | createBlockResolutions takes left and right child resolution data for a block (same display tree)
-- and generates node resolution data
createBlockResolutions :: Bool
                       -> LG.Node
                       -> Int
                       -> Int
                       -> NodeType
                       -> NodeType
                       -> NodeType
                       -> ResolutionBlockData
                       -> ResolutionBlockData
                       -> V.Vector CharInfo
                       -> ResolutionBlockData
createBlockResolutions
  compress
  curNode
  leftIndex
  rightIndex
  leftChildNodeType
  rightChildNodeType
  curNodeNodeType
  leftChild
  rightChild
  charInfoV
  | null leftChild && null rightChild = mempty
  | null leftChild = rightChild
  | null rightChild = leftChild
  | otherwise =
    -- trace ("CBR:" ++ (show (leftIndex, leftChildNodeType, rightIndex, rightChildNodeType)) ++ (show $fmap BV.toBits $ fmap displayBVLabel leftChild) ++ " and " ++  (show $fmap BV.toBits $ fmap displayBVLabel rightChild)) (
    let childResolutionPairs = cartProd (V.toList leftChild) (V.toList rightChild)
        -- need to keep these indices correct (hence reverse in checkLeafOverlap) for traceback and compress
        childResolutionIndices = cartProd [0.. (length leftChild - 1)] [0.. (length rightChild - 1)]
        validPairs = checkLeafOverlap (zip childResolutionPairs childResolutionIndices) []

        -- either parallel seems about the same
        -- newResolutionList = fmap (createNewResolution curNode leftIndex rightIndex leftChildNodeType rightChildNodeType charInfoV) validPairs `using` PU.myParListChunkRDS
        newResolutionList = PU.seqParMap rdeepseq  (createNewResolution curNode leftIndex rightIndex leftChildNodeType rightChildNodeType charInfoV) validPairs 

        --need to add in node and edge to left and right
        edgeLable = EdgeInfo { minLength = 0.0
                             , maxLength = 0.0
                             , midRangeLength = 0.0
                             , edgeType = TreeEdge
                             }
        newMinVertex = VertexInfo { index = curNode
                                  , bvLabel = BV.fromBits [False]
                                  , parents = mempty
                                  , children = mempty
                                  , nodeType = curNodeNodeType
                                  , vertName = T.pack $ "HTU" ++ show curNode
                                  , vertData = mempty
                                  , vertexResolutionData = mempty
                                  , vertexCost = 0.0
                                  , subGraphCost = 0.0
                                  }

        newNode = (curNode, newMinVertex)

        addLeft =  if leftChildNodeType == NetworkNode then
                        let newEdge = (curNode, rightIndex, edgeLable)
                            newRightChildBlockResolutionData = addNodeEdgeToResolutionList newNode newEdge 0 [] rightChild
                        in
                        -- trace ("ANEL:" ++ (show $ (curNode, rightIndex)))
                        newRightChildBlockResolutionData
                   else -- trace ("ANEL-Nothing")
                        mempty

        addRight = if rightChildNodeType == NetworkNode then
                        let newEdge = (curNode, leftIndex, edgeLable)
                            newLeftChildBlockResolutionData = addNodeEdgeToResolutionList newNode newEdge 0 [] leftChild
                        in
                        -- trace ("ANER:" ++ (show $ (curNode, leftIndex)))
                        newLeftChildBlockResolutionData
                   else -- trace ("ANER-Nothing")
                        mempty


    in
   --  trace ("CBR:" ++ (show $ leftIndex == rightIndex)) (
    -- trace ("=> " ++ (show $fmap BV.toBits $ fmap displayBVLabel totalResolutions) )
    -- compress to unique resolutions--can loose display trees, but speed up post-order a great deal
    -- trace ("CBR Num res left: " ++ (show $ V.length leftChild) ++ " Num res right: " ++ (show $ V.length rightChild) ++ " =>NRL " ++ (show $ length newResolutionList) ++ " addleft " ++ (show $ length addLeft)  ++ " addright " ++ (show $ length addRight)) (
    if compress then  V.fromList (nubResolutions newResolutionList []) V.++ (addLeft V.++ addRight)
    else V.fromList newResolutionList V.++ (addLeft V.++ addRight)
    -- )
    -- )

-- | createNewResolution takes a pair of resolutions and creates the median resolution
-- need to watch let/right (based on BV) for preorder stuff
createNewResolution :: LG.Node -> Int -> Int -> NodeType -> NodeType -> V.Vector CharInfo -> ((ResolutionData, ResolutionData),(Int, Int)) -> ResolutionData
createNewResolution curNode leftIndex rightIndex leftChildNodeType rightChildNodeType charInfoV ((leftRes, rightRes), (leftResIndex, rightResIndex)) =
    let -- make  bvLabel for resolution
        resBV = displayBVLabel leftRes .|. displayBVLabel rightRes

        -- Make resolution Display tree infomation
        leftEdgeType
          | leftChildNodeType == NetworkNode = NetworkEdge
          | leftChildNodeType == LeafNode = PendantEdge
          | otherwise = TreeEdge
        rightEdgeType
          | rightChildNodeType == NetworkNode = NetworkEdge
          | rightChildNodeType == LeafNode = PendantEdge
          | otherwise = TreeEdge

        edgeLable = EdgeInfo { minLength = 0.0
                             , maxLength = 0.0
                             , midRangeLength = 0.0
                             , edgeType = TreeEdge
                             }

        leftEdge =  (curNode, leftIndex, edgeLable {edgeType = leftEdgeType})
        rightEdge = (curNode, rightIndex, edgeLable {edgeType = rightEdgeType})
        leftChildTree = displaySubGraph leftRes
        rightChildTree = displaySubGraph rightRes

        -- Data fields empty for display tree data--not needed and muptiple copies of everything
        newNodeLabel = VertexInfo { index = curNode
                                  , bvLabel = resBV
                                  , parents = V.empty
                                  , children = V.fromList [leftIndex, rightIndex]
                                  , nodeType = TreeNode
                                  , vertName = T.pack $ "HTU" ++ show curNode
                                  , vertData = mempty
                                  , vertexResolutionData = mempty
                                  , vertexCost = 0.0
                                  , subGraphCost = 0.0
                                  }

        newNode = (curNode, newNodeLabel)

        -- this check for redundant edges in resoluton cash from combinations
        -- resolutionEdgeList = leftEdge : (rightEdge: (snd leftChildTree ++ snd rightChildTree))
        existingEdges = snd leftChildTree ++ snd rightChildTree
        resolutionEdgeList = if (leftEdge `notElem` existingEdges) && (rightEdge `notElem` existingEdges) then leftEdge : (rightEdge : existingEdges)
                             else if (leftEdge `elem` existingEdges) && (rightEdge `elem` existingEdges) then existingEdges
                             else if (leftEdge `notElem` existingEdges) then leftEdge : existingEdges
                             else rightEdge : existingEdges

        resolutionNodeList = newNode : (fst leftChildTree ++ fst rightChildTree)

        -- Make the data and cost for the resolution
        leftBlockLength = V.length $ displayData leftRes
        rightBlockLength = V.length $ displayData rightRes
        resolutionMedianCostV
          | (leftBlockLength == 0) = V.zip (displayData rightRes) (V.replicate rightBlockLength 0)
          | (rightBlockLength == 0) = V.zip (displayData leftRes) (V.replicate leftBlockLength 0)
          | otherwise = M.median2 (displayData leftRes) (displayData rightRes) charInfoV
        (resolutionMedianV, resolutionCostV) = V.unzip resolutionMedianCostV
        thisResolutionCost = V.sum resolutionCostV
        displaySubTreeCost = displayCost leftRes + displayCost rightRes + thisResolutionCost

    in
    ResolutionData { displaySubGraph = (resolutionNodeList, resolutionEdgeList)
                   , displayBVLabel = resBV
                   , displayData = resolutionMedianV
                   , childResolutions = [(Just leftResIndex, Just rightResIndex)]
                   , resolutionCost = thisResolutionCost
                   , displayCost = displaySubTreeCost
                   }

-- | nubResolutions takes a list of resolulutions and returns 'nub' based on leaf set (bvLabel) and Charinfo vector
nubResolutions :: [ResolutionData] -> [ResolutionData] -> [ResolutionData]
nubResolutions inData curData
  | null inData = reverse curData
  | null curData = nubResolutions (tail inData) [head inData]
  | otherwise =
    let firstData = head inData
        (isUnique, matchIndex) = hasResolutionMatch (displayBVLabel firstData) (displayData firstData) curData 0
    in
    if isUnique then nubResolutions (tail inData) (firstData : curData)
    else
        -- need to update the childResolutoin field of the match
        let firstPart = L.take matchIndex curData
            lastPart = L.drop (matchIndex + 1) curData
            matchChildResolutionData =  childResolutions (curData !! matchIndex)
            firstDataChildResolutionData = childResolutions firstData
            newMatchResolution = (curData !! matchIndex) {childResolutions = firstDataChildResolutionData ++ matchChildResolutionData}
        in
        nubResolutions (tail inData) (firstPart ++ (newMatchResolution : lastPart))

-- | hasResolutionMatch checks for match of bvlabel and Vect charinfo with list
hasResolutionMatch :: NameBV -> V.Vector CharacterData -> [ResolutionData] -> Int -> (Bool, Int)
hasResolutionMatch inBV inCD rDList curIndex =
    if null rDList then (True, curIndex)
    else
        let existingBV =  displayBVLabel $ head rDList
            existingCharData = displayData $ head rDList
        in
        if (existingBV /= inBV) || (existingCharData /= inCD) then hasResolutionMatch inBV inCD (tail rDList) (curIndex + 1) else (False, curIndex)

-- | checkLeafOverlap takes a left right resolution pair list and checks if
-- there is leaf overlap via comparing displayBVLabel if & = 0 then no
-- overlap, and adds to resulting list--reverses order--sholdn't matter
checkLeafOverlap :: [((ResolutionData, ResolutionData), (Int, Int))] -> [((ResolutionData, ResolutionData), (Int, Int))] -> [((ResolutionData, ResolutionData), (Int, Int))]
checkLeafOverlap inPairList curPairList =
    if null inPairList then reverse curPairList
    else
        let inPair@((leftRes, rightRes), (_, _)) = head inPairList
            leftBV = displayBVLabel leftRes
            rightBV = displayBVLabel rightRes
        in
        if BV.isZeroVector (leftBV .&. rightBV) then checkLeafOverlap (tail inPairList) (inPair : curPairList)
        else checkLeafOverlap (tail inPairList) curPairList




-- | addNodeAndEdgeToResolutionData adds new node and edge to resolution data in outdegree = 1 nodes
-- striaght copy would not add this node or edge to subtree in resolutions
addNodeAndEdgeToResolutionData :: LG.LNode VertexInfo -> LG.LEdge EdgeInfo -> V.Vector ResolutionBlockData -> V.Vector ResolutionBlockData
addNodeAndEdgeToResolutionData newNode newEdge = fmap (addNodeEdgeToResolutionList newNode newEdge 0 [])

-- | addNodeEdgeToResolutionList adds new node and edge to single subGraph in ResolutionData
-- adds resolutoin pairs to be equal to the child straight one-for-one correpondance
addNodeEdgeToResolutionList :: LG.LNode VertexInfo -> LG.LEdge EdgeInfo -> Int -> [ResolutionData] -> V.Vector ResolutionData -> V.Vector ResolutionData
addNodeEdgeToResolutionList newNode newEdge resolutionIndex curData inData =
    if null inData then V.fromList $ reverse  curData
    else
        let firstInData = V.head inData
            (inNodeList, inEdgeList) = displaySubGraph firstInData
            
            -- childResolutionIndexPairList = childResolutions firstInData
            newNodeList = newNode : inNodeList
            
            -- this check for redundant edges in resoluton cash from combinations
            newEdgeList = if (newEdge `notElem` inEdgeList) then newEdge : inEdgeList
                          else inEdgeList
            newFirstData = firstInData { displaySubGraph  = (newNodeList, newEdgeList)
                                        -- this apir in case of left/right issues later
                                        -- not sure this is correct--LWys first for children of out = 1 node
                                       , childResolutions = [(Just 0, Just 0)] -- [(Just resolutionIndex, Just resolutionIndex)]
                                       }
    in
    addNodeEdgeToResolutionList newNode newEdge (resolutionIndex + 1) (newFirstData : curData) (V.tail inData)


-- | getAllResolutionList takes ResolutionBlockData and retuns a list of the all valid (ie all leaves in subtree) display trees
-- for that block- and costs
getAllResolutionList :: ResolutionBlockData -> [(DecoratedGraph, VertexCost)]
getAllResolutionList  inRDList =
    --trace ("GBRL: " ++ (show inRDList)) (
    if null inRDList then error "Null resolution list"
    else
        let displayTreeList = fmap displaySubGraph inRDList
            displayCostList = fmap displayCost inRDList
            displayPopList = fmap (complement . displayBVLabel) inRDList
        in
            let displayBVList = V.zip3 displayTreeList displayCostList displayPopList
                validDisplayList = V.filter (BV.isZeroVector . thd3) displayBVList
                (displayList, costList, _) = V.unzip3 validDisplayList
            in
            --trace ("Valid display list number:" ++ (show $ length validDisplayList)) (
            if V.null validDisplayList then error ("Null validDisplayList in getAllResolutionList" ++ show inRDList)
            else
                let lDisplayTreeList = fmap LG.mkGraphPair (V.toList displayList)
                    -- displayTreeList' = fmap (updateRootCost validMinCost) displayTreeList
                in
                zip lDisplayTreeList (V.toList costList)


-- | getBestResolutionList takes ResolutionBlockData and retuns a list of the best valid (ie all leaves in subtree) display trees
-- for that block-- if checkPopCount is True--otherwise all display trees of any cost and contitution
-- startVertex for a component-- to allow for not every leaf being in componnet but still getting softwired cost
getBestResolutionList :: Maybe Int -> Bool -> ResolutionBlockData -> ([DecoratedGraph], VertexCost)
getBestResolutionList startVertex checkPopCount inRDList =
    --trace ("GBRL: " ++ (show inRDList)) (
    if null inRDList then error "Null resolution list"
    else
        let displayTreeList = fmap displaySubGraph inRDList
            displayCostList = fmap displayCost inRDList
            displayPopList = fmap (complement . displayBVLabel) inRDList
        in
        if not checkPopCount then
            let minCost = minimum displayCostList
                displayCostPairList = V.zip displayTreeList displayCostList
                (bestDisplayList, _) = V.unzip $ V.filter ((== minCost) . snd) displayCostPairList
            in
            (fmap LG.mkGraphPair (V.toList bestDisplayList), minCost)
        else
            let minPopCount = minimum $ fmap popCount displayPopList  --max since complemented above
                displayBVList = V.zip3 displayTreeList displayCostList displayPopList

                -- must have all leaves if startvzertex == Nothing, component maximum otherwise
                -- this for getting cost of component of a softwired network
                validDisplayList = if startVertex == Nothing then V.filter (BV.isZeroVector . thd3) displayBVList
                                   else V.filter ((== minPopCount) . (popCount . thd3)) displayBVList
                validMinCost = V.minimum $ fmap snd3 validDisplayList
                (bestDisplayList, _, _) = V.unzip3 $ V.filter ((== validMinCost) . snd3) validDisplayList
            in
            --trace ("Valid display list number:" ++ (show $ length validDisplayList)) (
            if (startVertex == Nothing) && (V.null validDisplayList) then error ("Null root validDisplayList in getBestResolutionList" ++ (show (startVertex,inRDList)) ++ " This can be caused if the graphType not set correctly.")
            else
                let lDisplayTreeList = fmap LG.mkGraphPair (V.toList bestDisplayList)

                    -- update root cost of display trees for use later (e.g. net penalties, outputting display forrests)
                    lDisplayTreeList' = fmap (updateRootCost validMinCost) lDisplayTreeList
                in
                (lDisplayTreeList', validMinCost)
            -- )

            -- )

-- | getBestResolutionListPair takes ResolutionBlockData and retuns a list of the best valid (ie all leaves in subtree) display trees
-- for that block-- if checkPopCount is True--otherwise all display trees of any cost and contitution
-- startVertex for a component-- to allow for not every leaf being in componnet but still getting softwired cost
-- returns list of pairs
getBestResolutionListPair :: Maybe Int -> Bool -> ResolutionBlockData -> [(DecoratedGraph, VertexCost)]
getBestResolutionListPair startVertex checkPopCount inRDList =
    --trace ("GBRL: " ++ (show inRDList)) (
    if null inRDList then error "Null resolution list"
    else
        let displayTreeList = fmap displaySubGraph inRDList
            displayCostList = fmap displayCost inRDList
            displayPopList = fmap (complement . displayBVLabel) inRDList
        in
        if not checkPopCount then
            let minCost = minimum displayCostList
                displayCostPairList = V.zip displayTreeList displayCostList
                (bestDisplayList, minCostList) = V.unzip $ V.filter ((== minCost) . snd) displayCostPairList
            in
            zip (fmap LG.mkGraphPair (V.toList bestDisplayList)) (V.toList minCostList)
        else
            let minPopCount = minimum $ fmap popCount displayPopList  --max since complemented above
                displayBVList = V.zip3 displayTreeList displayCostList displayPopList

                -- must have all leaves if startvzertex == Nothing, component maximum otherwise
                -- this for getting cost of component of a softwired network
                validDisplayList = if startVertex == Nothing then V.filter (BV.isZeroVector . thd3) displayBVList
                                   else V.filter ((== minPopCount) . (popCount . thd3)) displayBVList
                validMinCost = V.minimum $ fmap snd3 validDisplayList
                (bestDisplayList, minCostList, _) = V.unzip3 $ V.filter ((== validMinCost) . snd3) validDisplayList
            in
            --trace ("Valid display list number:" ++ (show $ length validDisplayList)) (
            if (startVertex == Nothing) && (V.null validDisplayList) then error ("Null root validDisplayList in getBestResolutionListPair" ++ (show (startVertex,inRDList)) ++ " This can be caused if the graphType not set correctly.")
            else
                let lDisplayTreeList = fmap LG.mkGraphPair (V.toList bestDisplayList)

                    -- update root cost of display trees for use later (e.g. net penalties, outputting display forrests)
                    lDisplayTreeList' = fmap (updateRootCost validMinCost) lDisplayTreeList
                in
                zip lDisplayTreeList' (V.toList minCostList)
            -- )

            -- )

-- | updateRootCost updates the subGraphCost of the root node(s) with input value
-- new node is created, so original is deleted, new added, and original edges added back
-- since deleted when node is
-- assumes its a tree wiht a single root
updateRootCost :: VertexCost -> DecoratedGraph -> DecoratedGraph
updateRootCost newRootCost inGraph =
    let (rootIndex, rootLabel) = head $ LG.getRoots inGraph
        rootEdges = LG.out inGraph rootIndex
        newRootLabel = rootLabel {subGraphCost = newRootCost}
    in
    -- trace ("DCC: " ++ (show newRootCost))
    LG.insEdges rootEdges $ LG.insNode (rootIndex, newRootLabel) $ LG.delNode rootIndex inGraph

-- | makeCharacterGraph takes a blockGraph and creates a vector of character graphs
-- each with a single block and single character
-- updating costs
makeCharacterGraph :: DecoratedGraph -> V.Vector DecoratedGraph
makeCharacterGraph inBlockGraph =
  if LG.isEmpty inBlockGraph then V.empty
  else
    let numCharacters =  V.length $ V.head $ vertData $ snd $ head $ LG.labNodes inBlockGraph
        characterGraphList = if numCharacters > 0 then fmap (pullCharacter False inBlockGraph) [0.. (numCharacters - 1)]
                             -- missing data
                             else [pullCharacter True inBlockGraph 0]
    in
    if V.length (vertData $ snd $ head $ LG.labNodes inBlockGraph) /= 1 then error ("Number of blocks /= 1 in makeCharacterGraph :" ++ (show $ V.length (vertData $ snd $ head $ LG.labNodes inBlockGraph)))
    else
      -- trace ("Chars: " ++ show numCharacters)
      V.fromList characterGraphList

-- | pullCharacter takes a DecoratedGraph with a single block and
-- creates a new DecoratedGraph with a single character from the input index
pullCharacter :: Bool -> DecoratedGraph -> Int -> DecoratedGraph
pullCharacter isMissing inBlockGraph characterIndex =
  if LG.isEmpty inBlockGraph then LG.empty
  else
    let (inNodeIndexList, inNodeLabelList) = unzip $ LG.labNodes inBlockGraph
        characterLabelList = fmap (makeCharacterLabels isMissing characterIndex) inNodeLabelList
    in
    LG.mkGraph (zip inNodeIndexList characterLabelList) (LG.labEdges inBlockGraph)

-- | makeCharacterLabels pulls the index character label form the singleton block (via head)
-- and creates a singleton character label, updating costs to that of the character
-- NB the case of missing data is answered here by an "empty charcter"
-- could be better to have V.empty
-- isMIssingChar seems to be extraneous--not sure whey it was there.
makeCharacterLabels :: Bool -> Int -> VertexInfo -> VertexInfo
makeCharacterLabels isMissing characterIndex inVertexInfo =
  -- trace ("MCl in:" ++ (show inVertexInfo) ++ " " ++ (show characterIndex)) (
  let -- isMissingChar = (V.length $ (vertData inVertexInfo) V.! characterIndex) == 0
      newVertexData = V.head (vertData inVertexInfo) V.! characterIndex
      (newVertexCost, newSubGraphCost) = if isMissing then (0, 0)
                                         --if isMissing || isMissingChar then (0, 0)
                                         else (localCost newVertexData, globalCost newVertexData)
      -- newVertexCost = localCost newVertexData
      -- newSubGraphCost = globalCost newVertexData
  in
  -- trace ("MCL " ++ (show $ V.length $ vertData inVertexInfo) ++ " " ++ (show $ fmap  V.length $ vertData inVertexInfo) ) (
  -- trace ("MCL: " ++ (show isMissing) ++ " CI: " ++ (show characterIndex) ++ " " ++ (show $ V.length $ (vertData inVertexInfo) V.! characterIndex))
  inVertexInfo { vertData     = if not isMissing then V.singleton $ V.singleton newVertexData
                                -- if not isMissing && not isMissingChar then V.singleton $ V.singleton newVertexData
                                else V.singleton $ V.singleton emptyCharacter --V.empty
               , vertexCost   = newVertexCost
               , subGraphCost = newSubGraphCost
               }
    -- ) )

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
        -- )

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

                (newDisplayVect, newCharTreeVV) = divideDecoratedGraphByBlockAndCharacterTree newGraph

            in
            -- th curnode == root index for pruned subtrees
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

                newCharData = if staticIA then createVertexDataOverBlocksStaticIA  (vertData leftChildLabel) (vertData  rightChildLabel) blockCharInfo []
                              else createVertexDataOverBlocks  (vertData leftChildLabel) (vertData  rightChildLabel) blockCharInfo []
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

                (newDisplayVect, newCharTreeVV) = divideDecoratedGraphByBlockAndCharacterTree newGraph

            in
            -- th curnode == roiot index for pruned subtrees
            -- trace ("New vertex:" ++ (show newVertex) ++ " at cost " ++ (show newCost)) (
            -- Do we need to PO.divideDecoratedGraphByBlockAndCharacterTree if not root?  probbaly not

            --if nodeType newVertex == RootNode then (simpleGraph, subGraphCost newVertex, newGraph, mempty, PO.divideDecoratedGraphByBlockAndCharacterTree newGraph, blockCharInfo)
            if nodeType newVertex == RootNode || curNode == rootIndex then (simpleGraph, subGraphCost newVertex, newGraph, newDisplayVect, newCharTreeVV, blockCharInfo)
            else (simpleGraph, subGraphCost newVertex, newGraph, mempty, mempty, blockCharInfo)

            -- ) -- )

-- | createVertexDataOverBlocks is a partial application of generalCreateVertexDataOverBlocks with full (all charcater) median calculation
createVertexDataOverBlocks :: VertexBlockData
                           -> VertexBlockData
                           -> V.Vector (V.Vector CharInfo)
                           -> [V.Vector (CharacterData, VertexCost)]
                           -> V.Vector (V.Vector (CharacterData, VertexCost))
createVertexDataOverBlocks = generalCreateVertexDataOverBlocks M.median2

-- | createVertexDataOverBlocksNonExact is a partial application of generalCreateVertexDataOverBlocks with partial (non-exact charcater) median calculation
createVertexDataOverBlocksNonExact :: VertexBlockData
                                   -> VertexBlockData
                                   -> V.Vector (V.Vector CharInfo)
                                   -> [V.Vector (CharacterData, VertexCost)]
                                   -> V.Vector (V.Vector (CharacterData, VertexCost))
createVertexDataOverBlocksNonExact = generalCreateVertexDataOverBlocks M.median2NonExact

-- | createVertexDataOverBlocksStaticIA is an  application of generalCreateVertexDataOverBlocks with exact charcater median calculation
-- and IA claculation for dynmaic characters--not full optimizations
createVertexDataOverBlocksStaticIA :: VertexBlockData
                                   -> VertexBlockData
                                   -> V.Vector (V.Vector CharInfo)
                                   -> [V.Vector (CharacterData, VertexCost)]
                                   -> V.Vector (V.Vector (CharacterData, VertexCost))
createVertexDataOverBlocksStaticIA = generalCreateVertexDataOverBlocks M.median2StaticIA


-- | generalCreateVertexDataOverBlocks is a genreal version for optimizing all (Add, NonAdd, Matrix)
-- and only non-exact (basically sequence) characters based on the median function passed
-- The function takes data in blocks and block vector of char info and
-- extracts the triple for each block and creates new block data for parent node (usually)
-- not checking if vectors are equal in length
generalCreateVertexDataOverBlocks :: (V.Vector CharacterData -> V.Vector CharacterData -> V.Vector CharInfo -> V.Vector (CharacterData, VertexCost))
                                  -> VertexBlockData
                                  -> VertexBlockData
                                  -> V.Vector (V.Vector CharInfo)
                                  -> [V.Vector (CharacterData, VertexCost)]
                                  -> V.Vector (V.Vector (CharacterData, VertexCost))
generalCreateVertexDataOverBlocks medianFunction leftBlockData rightBlockData blockCharInfoVect curBlockData =
    if V.null leftBlockData then
        --trace ("Blocks: " ++ (show $ length curBlockData) ++ " Chars  B0: " ++ (show $ V.map snd $ head curBlockData))
        V.fromList $ reverse curBlockData
    else
        let leftBlockLength = length $ V.head leftBlockData
            rightBlockLength =  length $ V.head rightBlockData
            -- firstBlock = V.zip3 (V.head leftBlockData) (V.head rightBlockData) (V.head blockCharInfoVect)

            -- missing data cases first or zip defaults to zero length
            firstBlockMedian
              | (leftBlockLength == 0) = V.zip (V.head rightBlockData) (V.replicate rightBlockLength 0)
              | (rightBlockLength == 0) = V.zip (V.head leftBlockData) (V.replicate leftBlockLength 0)
              | otherwise = medianFunction (V.head leftBlockData) (V.head rightBlockData) (V.head blockCharInfoVect)
        in
        generalCreateVertexDataOverBlocks medianFunction (V.tail leftBlockData) (V.tail rightBlockData) (V.tail blockCharInfoVect) (firstBlockMedian : curBlockData)

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
            bestTreesEdgeList = L.nubBy LG.undirectedEdgeEquality $ concat $ fmap LG.edges bestTreeList
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
            bestTreesEdgeList = L.nubBy LG.undirectedEdgeEquality $ concat $ fmap LG.edges bestTreeList
            
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
        let blockTreeEdgeList = L.nubBy LG.undirectedEdgeEquality $ concatMap LG.edges blockTreeList
            numExtraEdges = length $ LG.undirectedEdgeMinus blockTreeEdgeList treeEdgeList
            blockCost = subGraphCost $ fromJust $ LG.lab (head blockTreeList) rootIndex
        in
        -- trace ("GBW: " ++ (show (numExtraEdges, blockCost, blockTreeEdgeList)) ++ "\n" ++ (show $ fmap (subGraphCost . snd) $ LG.labNodes (head blockTreeList)))
        blockCost * (fromIntegral numExtraEdges)

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
        blockResolutionLL = V.toList $ fmap getAllResolutionList (vertexResolutionData rootLabel)
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

{- Unused but could be with naive comparisons

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

-}