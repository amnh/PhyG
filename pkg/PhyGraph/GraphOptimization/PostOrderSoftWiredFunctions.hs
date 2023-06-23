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
                                                      , makeLeafGraphSoftWired
                                                      , getDisplayBasedRerootSoftWired
                                                      , divideDecoratedGraphByBlockAndCharacterTree
                                                      , postOrderTreeTraversal
                                                      , postDecorateTree
                                                      , createVertexDataOverBlocks
                                                      , createVertexDataOverBlocksStaticIA
                                                      , createVertexDataOverBlocksNonExact
                                                      , getBlockCost
                                                      , getW15NetPenalty
                                                      , getW15NetPenaltyFull
                                                      , getW23NetPenalty
                                                      , getW23NetPenaltyReduced
                                                      , getW15RootCost
                                                      , getNetPenalty
                                                      , getNetPenaltyReduced
                                                      ) where

import           Data.Bits
import qualified Data.List                                        as L
import           Data.Maybe
import qualified Data.Text.Lazy                                   as T
import qualified Data.Vector                                      as V
import           GeneralUtilities
import qualified GraphFormatUtilities                             as GFU
import qualified GraphOptimization.Medians                        as M
import qualified GraphOptimization.PostOrderSoftWiredFunctionsNew as NEW
import qualified Graphs.GraphOperations                           as GO
import qualified ParallelUtilities                                as PU
import           Types.Types
import qualified Utilities.LocalGraph                             as LG
import qualified Utilities.Utilities                              as U
-- import           Debug.Trace



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
        -- (_, _, _, netVertexList) = LG.splitVertexList inSimpleGraph
        displayTreeList = LG.generateDisplayTrees contractIn1Out1Nodes inSimpleGraph

        -- get root index
        rootIndex = if isJust startVertex  then fromJust startVertex
                    else fst $ head $ LG.getRoots inSimpleGraph

        -- get the best traversals of the best display trees for each block
        (bestTripleInfo, _) = getBestDisplayCharBlockList inGS inData leafGraph rootIndex 0 [] [] displayTreeList
        (blockCostList, bestDisplayTreeList, charTreeVectList) = unzip3 bestTripleInfo

        -- extract specific information to create the phylogenetic graph
        graphCost = sum blockCostList
        displayTreeVect = V.fromList bestDisplayTreeList
        charTreeVectVect = V.fromList charTreeVectList

        -- propagate display node assignment to canonical graph
        -- does not have correct VertInfo--just character assignments
        -- to fox would need to propagate (and update other vertinfo like BV) via postorder pass
        newCononicalGraph = NEW.backPortBlockTreeNodesToCanonicalGraph (GO.convertSimpleToDecoratedGraph inSimpleGraph) displayTreeVect


        -- create postorder Phylgenetic graph
        postOrderPhyloGraph = (inSimpleGraph, graphCost, newCononicalGraph, fmap (:[]) displayTreeVect, charTreeVectVect, (fmap thd3 blockDataVect))

    in
    -- trace ("NPOSW: " ++ (show $ fmap bvLabel $ fmap snd $  LG.labNodes newCononicalGraph) ++ "\nDisplay :" ++ (show $ fmap bvLabel $ fmap snd $  LG.labNodes $ V.head displayTreeVect))
    postOrderPhyloGraph

-- | getBestDisplayCharBlockList takes a Tree gets best rootings, compares to input list if there is one and takes better
-- returning triple of block cost, display tree, char vect from better tree
getBestDisplayCharBlockList :: GlobalSettings
                            -> ProcessedData
                            -> DecoratedGraph
                            -> Int
                            -> Int
                            -> [(VertexCost, DecoratedGraph, V.Vector DecoratedGraph)]
                            -> [PhylogeneticGraph]
                            -> [SimpleGraph]
                            -> ([(VertexCost, DecoratedGraph, V.Vector DecoratedGraph)], [PhylogeneticGraph])
getBestDisplayCharBlockList inGS inData leafGraph rootIndex treeCounter currentBestTriple currentBestTreeList displayTreeList =
    if null displayTreeList then
        -- trace ("\tExamined " ++ (show treeCounter) ++ " display trees")
        (currentBestTriple, currentBestTreeList)
    else
        -- trace ("GBDCBL Trees: " ++ (show $ length displayTreeList)) (
        -- take first graph
        let -- get number of threads for parallel evaluation of display trees
            -- can set +RTS -N1 if CPUTime is off
            numDisplayTreesToEvaluate = PU.getNumThreads

            firstGraphList = take numDisplayTreesToEvaluate displayTreeList

            -- diagnose post order as Tree
            staticIA = False
            outgroupDiagnosedTreeList = PU.seqParMap (parStrategy $ lazyParStrat inGS) (postOrderTreeTraversal inGS inData leafGraph staticIA (Just rootIndex)) firstGraphList

            -- do rerooting of character trees
            multiTraverseTreeList = PU.seqParMap (parStrategy $ lazyParStrat inGS) (getDisplayBasedRerootSoftWired' inGS Tree rootIndex) outgroupDiagnosedTreeList

            -- extract triple (relevent info)--sets if multitraverse (reroot characters) or not
            multiTraverseTripleList = if (multiTraverseCharacters inGS == True)  then PU.seqParMap (parStrategy $ lazyParStrat inGS) (getTreeTriple rootIndex) multiTraverseTreeList
                                      else PU.seqParMap (parStrategy $ lazyParStrat inGS) (getTreeTriple rootIndex) outgroupDiagnosedTreeList

            -- choose better vs currentBestTriple
            -- this can be folded for a list > 2
            newBestTriple = L.foldl' chooseBetterTriple currentBestTriple multiTraverseTripleList -- multiTraverseTree

            -- save best overall dysplay trees for later use in penalty phase
            newBestTreeList = GO.selectGraphsFull Best (maxBound::Int) 0.0 (-1) (multiTraverseTreeList ++ currentBestTreeList)
        in
        -- trace ("GBDCBL: " ++ (show (fmap snd6 currentBestTreeList, fmap snd6 newBestTreeList, fmap snd6 multiTraverseTreeList)))
        getBestDisplayCharBlockList inGS inData leafGraph rootIndex (treeCounter + (length firstGraphList)) newBestTriple newBestTreeList (drop numDisplayTreesToEvaluate displayTreeList)
        -- )

-- | getTreeTriple takes a phylogenetic gaph and returns the triple list of block cost, display tree, and character graphs
getTreeTriple :: LG.Node -> PhylogeneticGraph -> [(VertexCost, DecoratedGraph, V.Vector DecoratedGraph)]
getTreeTriple rootIndex inGraph =
    if LG.isEmpty (fst6 inGraph) then []
    else
        let blockCostList = V.toList $ fmap (getBlockCost rootIndex) (fft6 inGraph)
            graphTriple = zip3 blockCostList (L.replicate (length blockCostList) (thd6 inGraph))  ((V.toList . fft6) inGraph)
        in
        graphTriple

-- | chooseBetterTriple takes the current best triplet of graph data and compares to Phylogenetic graph
-- and creates a new triple of better block cost, displayGraph for blocks, and character graphs
chooseBetterTriple :: [(VertexCost, DecoratedGraph, V.Vector DecoratedGraph)]
                   -> [(VertexCost, DecoratedGraph, V.Vector DecoratedGraph)]
                   -> [(VertexCost, DecoratedGraph, V.Vector DecoratedGraph)]
chooseBetterTriple inTripleList newTripleList =

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

-- | postOrderSoftWiredTraversal is a wrapper to allow correct function choice for alternate softwired algorithms
postOrderSoftWiredTraversal :: GlobalSettings -> ProcessedData -> DecoratedGraph -> Bool -> Maybe Int -> SimpleGraph -> PhylogeneticGraph
postOrderSoftWiredTraversal inGS inData leafGraph _ startVertex inSimpleGraph =
    -- firt case shouldn't happen--just checking if naive is chosen
    if graphType inGS == Tree then postOrderSoftWiredTraversal' inGS inData leafGraph startVertex inSimpleGraph
    else if softWiredMethod inGS == ResolutionCache then postOrderSoftWiredTraversal' inGS inData leafGraph startVertex inSimpleGraph
    else naivePostOrderSoftWiredTraversal inGS inData leafGraph startVertex inSimpleGraph

-- | postOrderSoftWiredTraversal' sets up and calls postorder traversal on Soft-wired graph
-- at root
-- staticIA is ignored--but kept for functional polymorphism
-- ur-root = ntaxa is an invariant
postOrderSoftWiredTraversal' :: GlobalSettings -> ProcessedData -> DecoratedGraph -> Maybe Int -> SimpleGraph -> PhylogeneticGraph
postOrderSoftWiredTraversal' inGS inData@(_, _, blockDataVect) leafGraph startVertex inSimpleGraph =
    if LG.isEmpty inSimpleGraph then emptyPhylogeneticGraph
    else
         -- Assumes root is Number of Leaves--should be invariant everywhere
        let rootIndex = if startVertex == Nothing then V.length $ fst3 inData
                        else fromJust startVertex
            blockCharInfo = V.map thd3 blockDataVect
            -- newSoftWired = postDecorateSoftWired inGS inSimpleGraph leafGraph blockCharInfo rootIndex rootIndex
            newSoftWired = NEW.postDecorateSoftWired inGS inSimpleGraph leafGraph blockCharInfo rootIndex rootIndex
        in
        if (startVertex == Nothing) && (not $ LG.isRoot inSimpleGraph rootIndex) then
            let localRootList = fst <$> LG.getRoots inSimpleGraph
                localRootEdges = concatMap (LG.out inSimpleGraph) localRootList
                currentRootEdges = LG.out inSimpleGraph rootIndex
            in
            error ("Index "  ++ show rootIndex ++ " with edges " ++ show currentRootEdges ++ " not root in graph:" ++ show localRootList ++ " edges:" ++ show localRootEdges ++ "\n" ++ LG.prettify inSimpleGraph)
        else
            newSoftWired

-- | getDisplayBasedRerootSoftWired is a wrapper to allow correct function choice for alternate softwired algorithms
getDisplayBasedRerootSoftWired :: GlobalSettings -> GraphType -> LG.Node -> PhylogeneticGraph -> PhylogeneticGraph
getDisplayBasedRerootSoftWired inGS inGraphType rootIndex inPhyloGraph =
    -- check if doing rerooting--if not then return existing graph
    if inGraphType == Tree then getDisplayBasedRerootSoftWired' inGS inGraphType rootIndex inPhyloGraph
    else if softWiredMethod inGS == ResolutionCache then getDisplayBasedRerootSoftWired' inGS inGraphType rootIndex inPhyloGraph
    else naiveGetDisplayBasedRerootSoftWired inGS inGraphType rootIndex inPhyloGraph

-- | naiveGetDisplayBasedRerootSoftWired is the naive (based on all resolution display trees)
-- the work of getDisplayBasedRerootSoftWired' is already done (rerooting and all) in naivePostOrderSoftWiredTraversal
naiveGetDisplayBasedRerootSoftWired ::  GlobalSettings -> GraphType -> LG.Node -> PhylogeneticGraph -> PhylogeneticGraph
naiveGetDisplayBasedRerootSoftWired _ _ _ inPhyloGraph  =
    inPhyloGraph

-- | getDisplayBasedRerootSoftWired' takes a graph and generates reroot costs for each character of each block
-- based on rerooting the display tree for that block.
-- Written for soft-wired, but could be modified for tree (data split on vertdata not resolution data)
-- this is a differnt approach from that of "tree" where the decorated, canonical tree is rerooted and each character and block
-- cost determined from that single rerooting
-- this should help avoid the issue of rerooting complex, reticulate graphs and maintaining
-- all the condition (cycles, time consistency etc) that occur.
-- done correcly this should be able to be used for trees (all display trees same = cononical graph) as
-- well as softwired, but not for hardwired where reticulations are maintained.

-- this can be modified for Tree data structres--basically by starting with vertdata initially without
-- resolutoin data trace back--should be more efficient in many was than existing code

-- Input display trees are for reporting only and do not contain actual character data so must be "pulled"
-- from cononical Decorated graph (thd field)
-- the list :[] stuff due to potential list of diplay trees not employed here
getDisplayBasedRerootSoftWired' :: GlobalSettings -> GraphType -> LG.Node -> PhylogeneticGraph -> PhylogeneticGraph
getDisplayBasedRerootSoftWired' inGS inGraphType rootIndex inPhyloGraph@(a,b, decGraph, _,_,f)  =
    if LG.isEmpty (fst6 inPhyloGraph) then inPhyloGraph
    else
        let -- update with pass to retrieve vert data from resolution data
            -- Trfee allready has data in vertData field
            (inSimpleGraph, _, inDecGraph, inBlockGraphV', inBlockCharGraphVV', charInfoVV) = 
                if inGraphType == Tree then
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
            -- (newBlockDisplayTreeVect, newBlockCharGraphVV, blockCostV) = unzip3 (zipWith3  (rerootBlockCharTrees inGS rootIndex) (V.toList $ fmap head inBlockGraphV) (V.toList inBlockCharGraphVV) (V.toList charInfoVV) `using` PU.myParListChunkRDS)
            -- This is slower than myParListChunkRDS
            (newBlockDisplayTreeVect, newBlockCharGraphVV, blockCostV) = unzip3 (PU.seqParMap (parStrategy $ lazyParStrat inGS) (rerootBlockCharTrees' inGS rootIndex) $ zip3 (V.toList $ fmap head inBlockGraphV) (V.toList inBlockCharGraphVV) (V.toList charInfoVV))

            newCononicalGraph = NEW.backPortBlockTreeNodesToCanonicalGraph inDecGraph (V.fromList newBlockDisplayTreeVect)
        in
        -- trace ("GDBRS:" ++ (show (b, sum blockCostV)))
        (inSimpleGraph, sum blockCostV, newCononicalGraph, V.fromList $ fmap (:[]) newBlockDisplayTreeVect, V.fromList newBlockCharGraphVV, charInfoVV)

-- | rerootBlockCharTrees' wrapper around rerootBlockCharTrees to allow for parMap
rerootBlockCharTrees' ::GlobalSettings -> LG.Node -> (DecoratedGraph, V.Vector DecoratedGraph, V.Vector CharInfo) -> (DecoratedGraph, V.Vector DecoratedGraph, VertexCost)
rerootBlockCharTrees' inGS rootIndex (blockDisplayTree, charTreeVect, charInfoVect) = rerootBlockCharTrees inGS rootIndex blockDisplayTree charTreeVect charInfoVect 

-- | rerootBlockCharTrees reroots all character trees (via fmap) in block returns best block char trees and costs
-- with best character tree node assignment back ported to display tree
rerootBlockCharTrees ::GlobalSettings -> LG.Node -> DecoratedGraph -> V.Vector DecoratedGraph -> V.Vector CharInfo -> (DecoratedGraph, V.Vector DecoratedGraph, VertexCost)
rerootBlockCharTrees inGS rootIndex blockDisplayTree charTreeVect charInfoVect =
    if V.null charTreeVect then error "Empty tree vector in rerootBlockCharTrees"
    else
        let -- next edges (to vertex in list) to perform rerooting
            -- progresses recursively over adjacent edges to minimize node reoptimization
            -- since initially all same graph can get initial reroot nodes from display tree
            childrenOfRoot = LG.descendants blockDisplayTree rootIndex
            grandChildrenOfRoot = concatMap (LG.descendants blockDisplayTree) childrenOfRoot

            -- leaving  parallel since can be few blocks
            -- (rerootedCharTreeVect, rerootedCostVect) = unzip (zipWith (getCharTreeBestRoot rootIndex grandChildrenOfRoot) (V.toList charTreeVect) (V.toList charInfoVect)) `using` PU.myParListChunkRDS)

            -- unclear if faster than than myParListChunkRDS
            (rerootedCharTreeVect, rerootedCostVect) = unzip (PU.seqParMap (parStrategy $ lazyParStrat inGS) (getCharTreeBestRoot' rootIndex grandChildrenOfRoot) (zip (V.toList charTreeVect) (V.toList charInfoVect)))

            (updateBlockDisplayTree, updatedDisplayVect, blockCost) = if multiTraverseCharacters inGS == True then
                                                                        (backPortCharTreeNodesToBlockTree blockDisplayTree (V.fromList rerootedCharTreeVect), V.fromList rerootedCharTreeVect, sum rerootedCostVect)
                                                                      else
                                                                         let rootCharLabelNodes = fmap (LG.labelNodeFlip rootIndex) charTreeVect
                                                                             existingCost = sum $ fmap (subGraphCost . snd) rootCharLabelNodes
                                                                         in
                                                                         (backPortCharTreeNodesToBlockTree blockDisplayTree charTreeVect, charTreeVect, existingCost)
        in
        (updateBlockDisplayTree, updatedDisplayVect, blockCost)

-- | getCharTreeBestRoot' is awrapper around getCharTreeBestRoot to use parMap
getCharTreeBestRoot' :: LG.Node -> [LG.Node] -> (DecoratedGraph, CharInfo) -> (DecoratedGraph, VertexCost)
getCharTreeBestRoot' rootIndex nodesToRoot (inCharacterGraph, charInfo) =
    getCharTreeBestRoot rootIndex nodesToRoot inCharacterGraph charInfo

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
        -- trace ("RRCT: " ++ (show (newGraphCost, bestCost)))
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
                                         , nodeType = GO.getNodeType inGraph curNodeIndex -- nodeType curNodeLabel
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


-- | backPortCharTreeNodesToBlockTree assigned nodes states (labels) of character trees to block display Tree
-- updates vertData, vertexCost, and subGraphCost for each.  Subgraph cost questionable since relies on rooting
backPortCharTreeNodesToBlockTree :: DecoratedGraph -> V.Vector DecoratedGraph -> DecoratedGraph
backPortCharTreeNodesToBlockTree blockDisplayTree rerootedCharTreeVect =
    let blockDisplayNodes = LG.labNodes blockDisplayTree
        blockDisplayEdges = LG.labEdges blockDisplayTree

        -- vector (characters) of vector (nodes) of labels
        charTreeLabelsVV = fmap V.fromList $ fmap (fmap snd) $ fmap LG.labNodes rerootedCharTreeVect

        -- for each node index extract (head head) vertdata, vertexCost and subgraphcost
        -- (vertDataVV, vertCostVV, subGraphCostVV) = V.unzip3 $ fmap (extractTripleVect charTreeLabelsVV) (V.fromList [0..(length blockDisplayNodes - 1)])
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
updateAndFinalizePostOrderSoftWired startVertexMaybe startVertex inGraph =
    if isNothing startVertexMaybe then NEW.softWiredPostOrderTraceBack startVertex inGraph
    else NEW.softWiredPostOrderTraceBack (fromJust startVertexMaybe) inGraph

-- | makeLeafGraphSoftWired takes input data and creates a 'graph' of leaves with Vertex information
-- but with zero edges.  This 'graph' can be reused as a starting structure for graph construction
-- to avoid remaking of leaf vertices
-- includes leave resolution data
makeLeafGraphSoftWired :: GlobalSettings -> ProcessedData -> DecoratedGraph
makeLeafGraphSoftWired inGS inData@(nameVect, bvNameVect, blocDataVect) =
    if V.null nameVect then error "Empty ProcessedData in makeLeafGraph"
    else if softWiredMethod inGS == Naive then GO.makeLeafGraph inData
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
        -- trace ("MVSW" ++ (show (localIndex, vertexResolutionData newVertex)))
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
                                                , childResolutionIndices = (Just 0, Just 0)
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
postOrderTreeTraversal inGS (_, _, blockDataVect) leafGraph staticIA startVertex inGraph  =
    if LG.isEmpty inGraph then emptyPhylogeneticGraph
    else
        -- Assumes root is Number of Leaves
        let rootIndex = if startVertex == Nothing then  fst $ head $ LG.getRoots inGraph
                        else fromJust startVertex
            blockCharInfo = V.map thd3 blockDataVect

            {-
            -- Hardwired--Not sure whey these edges can occur--somethng about adding edges in after not deleting them when assuming so
            inGraph' = if graphType inGS == Tree then inGraph
                       else (LG.removeNonLeafOut0NodesAfterRoot . LG.removeDuplicateEdges) inGraph
            -}
            newTree = postDecorateTree inGS staticIA inGraph leafGraph blockCharInfo rootIndex rootIndex
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

-- | postDecorateTree begins at start index (usually root, but could be a subtree) and moves preorder till children are labelled and then returns postorder
-- labelling vertices and edges as it goes back to root
-- this for a tree so single root
postDecorateTree :: GlobalSettings ->  Bool -> SimpleGraph -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.Node -> LG.Node -> PhylogeneticGraph
postDecorateTree inGS staticIA simpleGraph curDecGraph blockCharInfo rootIndex curNode =
    -- if node in there (leaf) or Hardwired network nothing to do and return
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
            leftChildTree = postDecorateTree inGS staticIA simpleGraph curDecGraph blockCharInfo rootIndex leftChild
            rightLeftChildTree = if length nodeChildren == 2 then postDecorateTree inGS staticIA simpleGraph (thd6 leftChildTree) blockCharInfo rootIndex rightChild
                                 else leftChildTree
            newSubTree = thd6 rightLeftChildTree
            ((_, leftChildLabel), (_, rightChildLabel)) = U.leftRightChildLabelBVNode (LG.labelNode newSubTree leftChild, LG.labelNode newSubTree rightChild)

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
                              else createVertexDataOverBlocks  (vertData leftChildLabel) (vertData rightChildLabel) blockCharInfo []

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
                                        -- this cost is incorrrect for Harwired netwqork fix at root
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

            -- Graph cost is calculated differently for Tree and Hardwired.  Sub trees can be counted multiple times
            -- in hardwired for outdegree two nodes with one or more network nodes as descendents
            -- this cannot be dealt with at the local node since there can be network all over the graph
            -- so simple add up the local costs of all nodes

            if nodeType newVertex == RootNode || curNode == rootIndex then
                -- Need full info for building trees
                let localCostSum = sum $ fmap vertexCost $ fmap snd $ LG.labNodes newGraph
                    -- updatedDisplayVect = V.zipWith NEW.backPortBlockTreeNodesToCanonicalGraph (fmap head newDisplayVect) newCharTreeVV
                    -- updatedCanonicalGraph = NEW.backPortBlockTreeNodesToCanonicalGraph newGraph updatedDisplayVect
                in
                -- trace ("PDT End: " ++ (show (subGraphCost newVertex, localCostSum)))
                -- (LG.removeDuplicateEdges simpleGraph, localCostSum, LG.removeDuplicateEdges newGraph, fmap (fmap LG.removeDuplicateEdges) newDisplayVect, fmap (fmap LG.removeDuplicateEdges) newCharTreeVV, blockCharInfo)
                -- (simpleGraph, localCostSum, updatedCanonicalGraph, fmap (:[]) updatedDisplayVect, newCharTreeVV, blockCharInfo)
                (simpleGraph, localCostSum, newGraph, newDisplayVect, newCharTreeVV, blockCharInfo)
            else (simpleGraph, subGraphCost newVertex, newGraph, mempty, mempty, blockCharInfo)

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

-- | getNetPenaltyReduced returns appropriate network penalty for a reduced graph
getNetPenaltyReduced :: GlobalSettings -> ProcessedData -> ReducedPhylogeneticGraph -> VertexCost
getNetPenaltyReduced  inGS inData inGraph =
    getNetPenalty inGS inData (GO.convertReduced2PhylogeneticGraph inGraph)

-- | getNetPenalty returns appropriate network penalty
getNetPenalty :: GlobalSettings -> ProcessedData -> PhylogeneticGraph -> VertexCost
getNetPenalty inGS inData inGraph =
    if (graphType inGS == Tree) then 0.0
    else if (graphFactor inGS) == NoNetworkPenalty then 0.0
    else if (graphFactor inGS) == Wheeler2015Network then getW15NetPenaltyFull Nothing inGS inData  Nothing inGraph
    else if (graphFactor inGS) == Wheeler2023Network then getW23NetPenalty inGraph
    else error ("Network penalty type " ++ (show $ graphFactor inGS) ++ " is not yet implemented")

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

-- | getW15NetPenaltyFull takes a Phylogenetic tree and returns the network penalty of Wheeler (2015)
-- does not use resolution cache's so can be used with Naive or Resolution cache SoftWired
-- has to be a single display tree--or could have no penalty for network since edges would be in one or other
-- display tree
getW15NetPenaltyFull :: Maybe ([VertexCost], [SimpleGraph], PhylogeneticGraph, Int) -> GlobalSettings -> ProcessedData -> Maybe Int -> PhylogeneticGraph -> VertexCost
getW15NetPenaltyFull blockInfo inGS inData@(nameVect, _, _) startVertex inGraph =
    if LG.isEmpty $ fst6 inGraph then 0.0
    else if LG.isTree $ fst6 inGraph then 0.0
    else
        -- have to do full data pasess on display trees
        if isNothing blockInfo then
            let rootIndex = if startVertex == Nothing then fst $ head $ LG.getRoots (fst6 inGraph)
                            else fromJust startVertex
                numLeaves' = V.length nameVect
                blockTreeList =  V.toList $ fmap GO.convertDecoratedToSimpleGraph $ fmap head $ fth6 inGraph
                blockCostList = V.toList $ fmap (getBlockCost rootIndex) (fft6 inGraph)

                -- get lowest cost display tree
                staticIA = False
                outgroupRootedList =  PU.seqParMap (parStrategy $ lazyParStrat inGS) (postOrderTreeTraversal inGS inData (GO.makeLeafGraph inData) staticIA (Just rootIndex)) blockTreeList
                multiTraverseTreeList = PU.seqParMap (parStrategy $ lazyParStrat inGS) (getDisplayBasedRerootSoftWired' inGS Tree rootIndex) outgroupRootedList
                lowestCostDisplayTree = head $ GO.selectGraphsFull Best 1 0.0 (-1) multiTraverseTreeList

                -- now can do as input (below)
                lowestCostEdgeList = (LG.edges . fst6) lowestCostDisplayTree
                lowestCostEdgesFlipped = fmap LG.flipEdge lowestCostEdgeList
                blockEdgeList = fmap LG.edges blockTreeList
                numBlockExtraEdgesList = fmap length $ fmap (L.\\ (lowestCostEdgeList ++ lowestCostEdgesFlipped)) blockEdgeList
                blockPenalty = sum $ zipWith (*) blockCostList (fmap fromIntegral numBlockExtraEdgesList)
            in
            -- trace ("GW15N: " ++ (show (blockEdgeList, numBlockExtraEdgesList, blockPenalty, blockPenalty / (4.0 * (fromIntegral numLeaves') - 4.0))))
            blockPenalty / (4.0 * (fromIntegral numLeaves') - 4.0)

        else
            let (blockCostList, blockTreeList, lowestCostDisplayTree, numLeaves) = fromJust blockInfo
                lowestCostEdgeList = (LG.edges . fst6) lowestCostDisplayTree
                lowestCostEdgesFlipped = fmap LG.flipEdge lowestCostEdgeList
                blockEdgeList = fmap LG.edges blockTreeList
                numBlockExtraEdgesList = fmap length $ fmap (L.\\ (lowestCostEdgeList ++ lowestCostEdgesFlipped)) blockEdgeList
                blockPenalty = sum $ zipWith (*) blockCostList (fmap fromIntegral numBlockExtraEdgesList)
            in
            -- trace ("GW15N: " ++ (show (blockEdgeList, numBlockExtraEdgesList, blockPenalty, blockPenalty / (4.0 * (fromIntegral numLeaves) - 4.0))))
            blockPenalty / (4.0 * (fromIntegral numLeaves) - 4.0)

-- | getW15NetPenalty takes a Phylogenetic tree and returns the network penalty of Wheeler (2015)
-- modified to take the union of all edges of trees of minimal length
-- currently modified -- not exactlty W15
getW15NetPenalty :: Maybe Int -> PhylogeneticGraph -> VertexCost
getW15NetPenalty startVertex inGraph =
    if LG.isEmpty $ thd6 inGraph then 0.0
    else if LG.isTree $ fst6 inGraph then 0.0
    else
        let -- (bestTreeList, _) = extractLowestCostDisplayTree startVertex inGraph
            bestTreeList =  V.toList $ fmap head $ fth6 inGraph
            bestTreesEdgeList = L.nubBy LG.undirectedEdgeEquality $ concat $ fmap LG.edges bestTreeList
            rootIndex = if startVertex == Nothing then fst $ head $ LG.getRoots (fst6 inGraph)
                        else fromJust startVertex
            blockPenaltyList = PU.seqParMap PU.myStrategy (getBlockW2015 bestTreesEdgeList rootIndex) (fth6 inGraph)

            -- leaf list for normalization
            (_, leafList, _, _) = LG.splitVertexList (fst6 inGraph)
            numLeaves = length leafList
            numTreeEdges = 4.0 * (fromIntegral numLeaves) - 4.0
            divisor = numTreeEdges
        in
        -- trace ("W15:" ++ (show ((sum $ blockPenaltyList) / divisor )))
        (sum $ blockPenaltyList) / divisor

-- | getW23NetPenaltyReduced takes a ReducedPhylogeneticGraph tree and returns the network penalty of Wheeler and Washburn (2023)
-- basic idea is new edge improvement must be better than average existing edge cost
-- penalty for each added edge (unlike W15 which was on a block by block basis--and requires additional tree diagnoses)
-- num extra edges/2 since actually add 2 new edges when one network edge
-- requires resolution cache data structures
getW23NetPenaltyReduced :: ReducedPhylogeneticGraph -> VertexCost
getW23NetPenaltyReduced inGraph =
    if LG.isEmpty $ thd5 inGraph then 0.0
    else if LG.isTree $ fst5 inGraph then 0.0
    else
        let -- (bestTreeList, _) = extractLowestCostDisplayTree startVertex inGraph
            bestTreeList =  V.toList $ fmap head $ fth5 inGraph
            bestTreesEdgeList = L.nubBy LG.undirectedEdgeEquality $ concat $ fmap LG.edges bestTreeList

            -- leaf list for normalization
            (_, leafList, _, _) = LG.splitVertexList (fst5 inGraph)
            numLeaves = length leafList
            numTreeEdges = 2.0 * (fromIntegral numLeaves) - 2.0
            numExtraEdges = ((fromIntegral $ length bestTreesEdgeList) - numTreeEdges) / 2.0
            -- divisor = numTreeEdges - numExtraEdges
        in
       --  trace ("W23:" ++ (show ((numExtraEdges * (snd6 inGraph)) / (2.0 * numTreeEdges))) ++ " from " ++ (show (numTreeEdges, numExtraEdges))) (
        -- if divisor == 0.0 then infinity
        -- else (sum blockPenaltyList) / divisor
        -- else (numExtraEdges * (sum blockPenaltyList)) / divisor
        --else
        (numExtraEdges * (snd5 inGraph)) / (2.0 * numTreeEdges)
        -- )

-- | getW23NetPenalty takes a Phylogenetic tree and returns the network penalty of Wheeler and Washburn (2023)
-- basic idea is new edge improvement must be better than average existing edge cost
-- penalty for each added edge (unlike W15 which was on a block by block basis--and requires additional tree diagnoses)
-- num extra edges/2 since actually add 2 new edges when one network edge
-- requires resolution cache data structures
getW23NetPenalty :: PhylogeneticGraph -> VertexCost
getW23NetPenalty inGraph =
    if LG.isEmpty $ thd6 inGraph then 0.0
    else if LG.isTree $ fst6 inGraph then 0.0
    else
        let -- (bestTreeList, _) = extractLowestCostDisplayTree startVertex inGraph
            bestTreeList =  V.toList $ fmap head $ fth6 inGraph
            bestTreesEdgeList = L.nubBy LG.undirectedEdgeEquality $ concat $ fmap LG.edges bestTreeList

            -- leaf list for normalization
            (_, leafList, _, _) = LG.splitVertexList (fst6 inGraph)
            numLeaves = length leafList
            numTreeEdges = 2.0 * (fromIntegral numLeaves) - 2.0
            numExtraEdges = ((fromIntegral $ length bestTreesEdgeList) - numTreeEdges) / 2.0
            -- divisor = numTreeEdges - numExtraEdges
        in
       --  trace ("W23:" ++ (show ((numExtraEdges * (snd6 inGraph)) / (2.0 * numTreeEdges))) ++ " from " ++ (show (numTreeEdges, numExtraEdges))) (
        -- if divisor == 0.0 then infinity
        -- else (sum blockPenaltyList) / divisor
        -- else (numExtraEdges * (sum blockPenaltyList)) / divisor
        --else
        (numExtraEdges * (snd6 inGraph)) / (2.0 * numTreeEdges)
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

