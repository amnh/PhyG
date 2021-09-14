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
                                    , preOrderTreeTraversal
                                    , makeLeafGraph
                                    , multiTraverseFullyLabelGraph
                                    , fullyLabelGraphWithoutMultiTraversal
                                    ) where

import           Types.Types
import qualified Data.List as L
import qualified Data.Vector as V
import qualified Utilities.LocalGraph as LG
import GeneralUtilities
import qualified GraphFormatUtilities as GFU
-- import qualified GraphOptimization.Medians as M
import qualified Graphs.GraphOperations  as GO
import qualified GraphOptimization.PostOrderFunctions as PO
import qualified GraphOptimization.PreOrderFunctions as PRE
-- import qualified SymMatrix as SM
-- import qualified Data.BitVector.LittleEndian as BV
import qualified Data.Text.Lazy  as T
import Data.Bits 
import Data.Maybe
import Debug.Debug as D
import Utilities.Utilities as U
import qualified GraphOptimization.Medians as M
import qualified SymMatrix                   as S
import qualified Data.Vector.Generic         as GV
import Debug.Trace
import qualified Debug.Debug as Debug

-- | multiTraverseFullyLabelGraph naively reroot (no progessive reroot) graph on all vertices and chooses lowest cost
-- does not yet minimize over characters to get multi-min
multiTraverseFullyLabelGraph :: GlobalSettings -> ProcessedData -> SimpleGraph -> PhylogeneticGraph
multiTraverseFullyLabelGraph inGS inData inGraph =
    if LG.isEmpty inGraph then emptyPhylogeneticGraph
    else 
        let leafGraph = makeLeafGraph inData

            -- for special casing of nonexact and single exact characters
            nonExactChars = U.getNumberNonExactCharacters (thd3 inData)

            -- minimal reoptimize with smart order to minimize node reoptimization

            -- initial traversal based on global outgroup and the "next" traversal points as children of existing traversal
            -- here initial root.
            outgroupRootedPhyloGraph = postOrderTreeTraversal inGS inData leafGraph $ GO.rerootGraph' inGraph (outgroupIndex inGS)
            childrenOfRoot = concat $ fmap (LG.descendants (thd6 outgroupRootedPhyloGraph)) (fmap fst $ LG.getRoots $ thd6 outgroupRootedPhyloGraph)

            -- create list of multi-traversals with original rooting first
            -- subsequent rerooting do not reoptimize exact characters (add nonadd etc) 
            -- they are taken from the fully labelled first input decorated graph later when output graph created
            -- it is important that teh first graph be the ourgroup rooted graph (outgroupRootedPhyloGraph) so this 
            -- will have the preoder assignmentsd for th eoutgroup rooted graph as 3rd filed.  This can be used for incremental
            -- optimization to get log n initial postorder assingment when mutsating graph.
            recursiveRerootList = outgroupRootedPhyloGraph : minimalReRootPhyloGraph outgroupRootedPhyloGraph childrenOfRoot

            minCostRecursive = minimum $ fmap snd6 recursiveRerootList
            minCostGraphListRecursive = filter ((== minCostRecursive).snd6) recursiveRerootList

            -- create optimal final graph with best costs and best traversal (rerooting) forest for each character
            --  travesal for exact characters (and costs) are the first of each least since exact only optimizaed for that 
            --  traversal graph.  The result has approprotate post-order assignments for traversals, preorder "final" assignments
            -- are propagated to the Decorated graph field after the preorder pass.
            --graphWithBestAssignments  = setBestGraphAssignments recursiveRerootList (fst6 outgroupRootedPhyloGraph)  (thd6 outgroupRootedPhyloGraph) (six6 outgroupRootedPhyloGraph) 
            graphWithBestAssignments' = L.foldl1' setBetterGraphAssignment recursiveRerootList -- (recursiveRerootList !! 0) (recursiveRerootList !! 1) 

        in
        --trace ("Outgroup cost:" ++ show (snd6 outgroupRootedPhyloGraph))
        --trace ("Initial Children: " ++ show childrenOfRoot)
        --trace ("Min graph cost :" ++ show minCostRecursive ++ " Merged :" ++ (show $ snd6 graphWithBestAssignments'))
        --    show minCostDirect ++ ":" ++ (show $ sort $ fmap snd6 rerootPhyloGraphListDirect)
        --    ++ "\n" ++ show minCostRecursive ++ ":" ++ (show $ sort $ fmap snd6 recursiveRerootList))

        -- for debugging
        --trace ("Top " ++ (show $ fmap (fmap charType) $ six6 outgroupRootedPhyloGraph))
        -- graphWithBestAssignments'

        -- Uncomment this to (and comment the following three cases) avoid traversal rerooting stuff for debugging
        -- preOrderTreeTraversal outgroupRootedPhyloGraph
        
        -- special cases that don't require all the work
        if nonExactChars == 0 then preOrderTreeTraversal outgroupRootedPhyloGraph
            
        else if nonExactChars == 1 then preOrderTreeTraversal $ head minCostGraphListRecursive
        
        else preOrderTreeTraversal graphWithBestAssignments'
        
        

-- | setBetterGraphAssignment takes two phylogenetic graphs and returns the lower cost optimization of each character,
-- with traversal focus etc to get best overall graph
-- since this is meant to work with graphs that have or do not have reoptimized exact (=static-Add/NonAdd/MAtrix) characters 
-- the criterion is lower cost character is taken, unless the cost is zero, then non-zero is taken
-- this function is expected to be used in a fold over a list of graphs
-- the basic comparison is over the costs of the root(s) cost for each  of the character decorated (traversal) graphs

--May change
-- assumes that a single decorated graph comes in for each Phylogenetic graph from the fully and reroot optimize (V.singleton (V.singleton DecGraph))
-- and goes through the block-character-cost data and reassigns based on that creating a unique (although there could be more than one) decorated 
-- graph for each character in each block.
-- postorder assignments in traversal set of block character trees are NOT propagated back to first decorated graph. 
-- the thrid field of phylogenetic Graph is set to the 3rd fiuled of the first of two inputs--so if startiong fold with outgroup
-- rooted graph--that is waht stays which can be used as a preorder graph for incremental optimization
-- when perfoming that sort of operation
-- The traversal graphs
-- are used for the pre-order final assignments whic will be propagated back to set those of the 3rd filed decorated graph

-- this will have to be modified for solf-wired since incoming blocks will not all be the same underlying gaph
-- unclear how hardwired will be affected
setBetterGraphAssignment :: PhylogeneticGraph -> PhylogeneticGraph -> PhylogeneticGraph
setBetterGraphAssignment firstGraph@(fSimple, _, fDecGraph, fBlockDisplay, fTraversal, fCharInfo) secondGraph@(_, _, sDecGraph, _, sTraversal, _) =
    if LG.isEmpty fDecGraph then secondGraph
    else if LG.isEmpty sDecGraph then firstGraph
    else 
        --trace ("setBetter (" ++ (show fCost) ++ "," ++ (show sCost) ++ ")" ) ( -- ++ " CharInfo blocks:" ++ (show $ length fCharInfo) ++ " characters: " ++ (show $ fmap length fCharInfo) ++ " " ++ (show $ fmap (fmap name) fCharInfo)) (
        let (mergedBlockVect, costVector) = V.unzip $ fmap makeBetterBlock (D.debugVectorZip fTraversal sTraversal)
        in
         --trace ("setBetter (" ++ (show fCost) ++ "," ++ (show sCost) ++ ") ->" ++ (show $ V.sum costVector) ++ " nt:" ++ (show $ length fTraversal)
         --   ++ "length blocks " ++ (show $ fmap length fTraversal))
        (fSimple, V.sum costVector, fDecGraph, fBlockDisplay, mergedBlockVect, fCharInfo)
        --)
        
-- | makeBetterBlocktakes two verctors of traversals. Each vector contains a decorated graph (=traversla graph) for each 
-- character.  This can be a single sequence or series of exact characters
-- the function applies a character cost comparison to get the better 
makeBetterBlock :: ((V.Vector DecoratedGraph), (V.Vector DecoratedGraph)) -> ((V.Vector DecoratedGraph), VertexCost)
makeBetterBlock (firstBlockGraphVect, secondBlockGraphVect) = 
    let (mergedCharacterVect, costVector) = V.unzip $ fmap chooseBetterCharacter (D.debugVectorZip firstBlockGraphVect secondBlockGraphVect)
    in
    (mergedCharacterVect, V.sum costVector)

-- | chooseBetterCharacter takes a pair of character decorated graphs and chooses teh "better" one as in lower cost, or non-zero cost
--  if one is zer (for exact characters) and returns the better character and cost
--  graph can have multiple roots
chooseBetterCharacter :: (DecoratedGraph, DecoratedGraph) -> (DecoratedGraph, VertexCost)
chooseBetterCharacter (firstGraph, secondGraph) =
    if LG.isEmpty firstGraph then error "Empty first graph in chooseBetterCharacter"
    else if LG.isEmpty secondGraph then error "Empty second graph in chooseBetterCharacter"
    else    
        let firstGraphCost = sum $ fmap subGraphCost $ fmap snd $ LG.getRoots firstGraph
            secondGraphCost = sum $ fmap subGraphCost $ fmap snd $ LG.getRoots secondGraph
        in
        -- trace ("Costs " ++ show (firstGraphCost, secondGraphCost)) (
        if firstGraphCost == 0 then (secondGraph, secondGraphCost) 
        else if secondGraphCost == 0 then (firstGraph, firstGraphCost)
        else if firstGraphCost < secondGraphCost then (firstGraph, firstGraphCost)
        else (secondGraph, secondGraphCost) 
        --)

-- | minimalReRootPhyloGraph takes an inialtial post-order labelled phylogenetic graph
-- and "intelligently" reroots by traversing through adjacent edges, hopefully
-- reoptimizing the minimum number of vertices each time (2) but could be more depending
-- on graph topology
-- NB--only deals with post-order assignments
minimalReRootPhyloGraph :: PhylogeneticGraph -> [LG.Node] -> [PhylogeneticGraph]
minimalReRootPhyloGraph inGraph nodesToRoot =
    if null nodesToRoot then []
    else 
        let firstRerootIndex = head nodesToRoot
            nextReroots = (LG.descendants (thd6 inGraph) firstRerootIndex) ++ (tail nodesToRoot)
            newGraph = PO.rerootPhylogeneticGraph' inGraph firstRerootIndex
        in
        --trace ("New cost:" ++ show (snd6 newGraph) ++ " vs " ++ (show $ GO.graphCostFromNodes $ thd6 newGraph))
        newGraph : minimalReRootPhyloGraph newGraph nextReroots

-- | fullyLabelGraphWithoutMultiTraversal takes an unlabelled "simple' graph, performs post and preorder passes to 
-- fully label the graph and return a PhylogenticGraph
-- this does no rerooting to imporve score
fullyLabelGraphWithoutMultiTraversal :: GlobalSettings -> ProcessedData -> DecoratedGraph -> SimpleGraph -> PhylogeneticGraph
fullyLabelGraphWithoutMultiTraversal inGS inData leafGraph inGraph = 
    if LG.isEmpty inGraph then emptyPhylogeneticGraph
    else 
        let postOrderTree = postOrderTreeTraversal inGS inData leafGraph inGraph
            preOrderTree = preOrderTreeTraversal postOrderTree
        in
        preOrderTree

-- | makeLeafGraph takes input data and creates a 'graph' of leaves with Vertex informnation
-- but with zero edges.  This 'graph' can be reused as a starting structure for graph construction
-- to avoid remaking of leaf vertices
makeLeafGraph :: ProcessedData -> DecoratedGraph
makeLeafGraph (nameVect, bvNameVect, blocDataVect) =
    if V.null nameVect then error "Empty ProcessedData in makeLeafGraph"
    else 
        let leafVertexList = V.toList $ V.map (makeLeafVertex nameVect bvNameVect blocDataVect) (V.fromList [0.. (V.length nameVect) - 1])
        in
        LG.mkGraph leafVertexList []

-- | makeLeafVertex makes a single unconnected vertex for a leaf
makeLeafVertex :: V.Vector NameText -> V.Vector NameBV -> V.Vector BlockData -> Int -> LG.LNode VertexInfo
makeLeafVertex nameVect bvNameVect inData localIndex =
    --trace ("Making leaf " ++ (show localIndex) ++ " Data " ++ (show $ length inData) ++ " " ++ (show $ fmap length $ fmap snd3 inData)) (
    let centralData = V.map snd3 inData 
        thisData = V.map (V.! localIndex) centralData
        newVertex = VertexInfo  { index = localIndex
                                , bvLabel = bvNameVect V.! localIndex
                                , parents = V.empty
                                , children = V.empty
                                , nodeType = LeafNode
                                , vertName =  nameVect V.! localIndex
                                , vertData = thisData
                                , vertexCost = 0.0
                                , subGraphCost = 0.0
                                }   
        in
        -- trace (show (length thisData) ++ (show $ fmap length thisData))
        (localIndex, newVertex)
        --)

-- | postOrderTreeTraversal takes a 'simple' graph and generates 'preliminary' assignments
-- vi post-order traversal, yields cost as well
-- for a binary tree only
-- depending on optimality criterion--will calculate root cost
postOrderTreeTraversal :: GlobalSettings -> ProcessedData -> DecoratedGraph -> SimpleGraph -> PhylogeneticGraph
postOrderTreeTraversal inGS inData@(_, _, blockDataVect) leafGraph inGraph = 
    if LG.isEmpty inGraph then emptyPhylogeneticGraph
    else
        -- Assumes root is Number of Leaves  
        let rootIndex = V.length $ fst3 inData
            blockCharInfo = V.map thd3 blockDataVect
            newTree = postDecorateTree inGS inData inGraph leafGraph blockCharInfo rootIndex
        in
        --trace ("It Begins at " ++ show rootIndex) (
        if not $ LG.isRoot inGraph rootIndex then 
            let localRootList = fmap fst $ LG.getRoots inGraph
                localRootEdges = concatMap (LG.out inGraph) localRootList
                currentRootEdges = LG.out inGraph rootIndex
            in
            error ("Index "  ++ (show rootIndex) ++ " with edges " ++ (show currentRootEdges) ++ " not root in graph:" ++ (show localRootList) ++ " edges:" ++ (show localRootEdges) ++ "\n" ++ (GFU.showGraph inGraph))
        else newTree 
        -- )
{-
Was used in postDecorateTree, but removed 

-- | getVirtualRootEdge cretes tehe virtual edge that would have existed to crete that node rott
-- the ide is that if the tree had been rooted somewhere else -- what edge would have existied and was
-- "divided" to crete this node.  This is used for individual character graph traversal foci
-- edge is basically undirected since orientation is unknown
getVirtualRootEdge :: (Show a, Show b) => LG.Gr a b -> LG.Node -> LG.Edge
getVirtualRootEdge inGraph inNode = 
    if LG.isEmpty inGraph then error "Empty graph in getVirtualRootEdge"
    else 
        let childList = LG.descendants inGraph inNode
            parentList = LG.parents inGraph inNode
        in
        -- this really shouldn't be used but in recursion for traversal may be needed
        if LG.isLeaf inGraph inNode then (head parentList, inNode)
        else 
            if length childList /= 2 then error ("Root node with /= 2 children in getVirtualRootEdge\n" ++ (GFU.showGraph inGraph)) 
            else (head childList, last childList)
-}

-- | postDecorateTree begins at start index (usually root, but could be a subtree) and moves preorder till childrend ane labelled and then reurns postorder
-- labelling vertices and edges as it goes back to root
-- this for a tree so single rtoot
postDecorateTree :: GlobalSettings -> ProcessedData -> SimpleGraph -> DecoratedGraph -> V.Vector (V.Vector CharInfo) -> LG.Node -> PhylogeneticGraph
postDecorateTree inGS inData simpleGraph curDecGraph blockCharInfo curNode = 
    -- if node in there nothing to do and return
    if LG.gelem curNode curDecGraph then 
        let nodeLabel = LG.lab curDecGraph curNode
            -- identify/create the virtual edge this node would have been created from
            -- this for traversal focus use for character 
            -- impliedRootEdge = getVirtualRootEdge curDecGraph curNode
        in
        if nodeLabel == Nothing then error ("Null label for node " ++ show curNode)
        else
            --  replicating curaDecGraph with number opf blocks--but all the same for tree 
            (simpleGraph, subGraphCost (fromJust nodeLabel), curDecGraph, V.replicate (length blockCharInfo) curDecGraph, V.singleton (V.singleton curDecGraph), blockCharInfo)

    -- Need to make node
    else 
        
        -- check if children in graph
        let nodeChildren = LG.descendants simpleGraph curNode  -- should be 1 or 2, not zero since all leaves already in graph
            leftChild = (head nodeChildren)
            rightChild = (last nodeChildren)
            leftChildTree = postDecorateTree inGS inData simpleGraph curDecGraph blockCharInfo leftChild
            rightLeftChildTree = if (length nodeChildren == 2) then postDecorateTree inGS inData simpleGraph (thd6 $ leftChildTree) blockCharInfo rightChild
                                 else leftChildTree
        in 

        if length nodeChildren > 2 then error ("Graph not dichotomous in postDecorateTree node " ++ (show curNode) ++ "\n" ++ LG.prettify simpleGraph)
        else if length nodeChildren == 0 then error ("Leaf not in graph in postDecorateTree node " ++ (show curNode) ++ "\n" ++ LG.prettify simpleGraph)

        -- make node from childern
        else    
            -- make node from children and new edges to children
            -- median2 takes characters in blocks--but for tree really all same block
            let newSubTree = thd6 $ rightLeftChildTree
                -- leftChildLabel = fromJust $ LG.lab newSubTree leftChild
                -- rightChildLabel = fromJust $ LG.lab newSubTree rightChild

                -- this ensures that left/right choices are based on leaf BV for consistency and label invariance
                -- larger bitvector is Right, smaller or equal Left 
                (leftChildLabel, rightChildLabel) = U.leftRightChildLabelBV (fromJust $ LG.lab newSubTree leftChild, fromJust $ LG.lab newSubTree rightChild)
                
                newCharData = PO.createVertexDataOverBlocks  (vertData leftChildLabel) (vertData  rightChildLabel) blockCharInfo []
                newCost =  V.sum $ V.map (V.sum) $ V.map (V.map snd) newCharData
                newVertex = VertexInfo {  index = curNode
                                        , bvLabel = (bvLabel leftChildLabel) .|. (bvLabel rightChildLabel)
                                        , parents = V.fromList $ LG.parents simpleGraph curNode
                                        , children = V.fromList nodeChildren
                                        , nodeType = GO.getNodeType simpleGraph curNode
                                        , vertName = T.pack $ "HTU" ++ (show curNode)
                                        , vertData = V.map (V.map fst) newCharData
                                        , vertexCost = newCost
                                        , subGraphCost = (subGraphCost leftChildLabel) + (subGraphCost rightChildLabel) + newCost
                                        }   
                newEdgesLabel = EdgeInfo {    minLength = newCost / 2.0
                                            , maxLength = newCost / 2.0
                                            , midRangeLength = newCost / 2.0
                                            , edgeType = TreeEdge
                                         }
                newEdges = fmap LG.toEdge $ LG.out simpleGraph curNode 
                newLEdges =  fmap (LG.toLEdge' newEdgesLabel) newEdges
                newGraph =  LG.insEdges newLEdges $ LG.insNode (curNode, newVertex) newSubTree                    
            in
            -- trace ("Orig graph: " ++ (show $ U.getDecoratedGraphBlockCharInformation newGraph))

            (simpleGraph, (subGraphCost newVertex), newGraph,V.replicate (length blockCharInfo) newGraph, PO.divideDecoratedGraphByBlockAndCharacter newGraph, blockCharInfo)
            --)
            

-- | preOrderTreeTraversal takes a preliminarily labelled PhylogeneticGraph
-- and returns a full labels with 'final' assignments based on character decorated graphs
-- created postorder (5th of 6 fields).  
-- the preorder states are creted by traversing the traversal DecoratedGraphs in the 5th filed of PhylogeneticGraphs
-- these are by block and character,  Exact charcaters are vectors of standard characters and each seqeunce (non-exact) 
-- has its own traversal graph. These should be treees here (could be forests later) and should all have same root (number taxa)
-- but worth checking to make sure.
-- these were creted by "splitting" them after preorder 
-- The final states are propagated back to the second field
-- DecoratedGraph of the full Phylogenetic graph--which does NOT have the character based preliminary assignments
-- ie postorder--since those are traversal specific
-- the character speciofic decorated graphs have appropriate post and pre-order assignments
-- the traversal begins at the root (for a tree) and proceeds to leaves.
preOrderTreeTraversal :: PhylogeneticGraph -> PhylogeneticGraph
preOrderTreeTraversal inPGraph@(inSimple, inCost, inDecorated, blockDisplayV, blockCharacterDecoratedVV, inCharInfoVV) = 
    --trace ("PreO:" ++ (show $ fmap (fmap charType) inCharInfoVV)) (
    if LG.isEmpty (thd6 inPGraph) then emptyPhylogeneticGraph
    else 
        -- trace ("In PreOrder\n" ++ "Simple:\n" ++ (LG.prettify inSimple) ++ "Decorated:\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inDecorated) ++ "\n" ++ (GFU.showGraph inDecorated)) (
        -- mapped recursive call over blkocks, later characters
        let preOrderBlockVect = fmap doBlockTraversal $ Debug.debugVectorZip inCharInfoVV blockCharacterDecoratedVV
            fullyDecoratedGraph = assignPreorderStatesAndEdges preOrderBlockVect inCharInfoVV inDecorated 
        in
        {-
        let blockPost = GO.showDecGraphs blockCharacterDecoratedVV
            blockPre = GO.showDecGraphs preOrderBlockVect
        in
        trace ("BlockPost:\n" ++ blockPost ++ "BlockPre:\n" ++ blockPre ++ "After Preorder\n" ++  (LG.prettify $ GO.convertDecoratedToSimpleGraph fullyDecoratedGraph))
        -}
        (inSimple, inCost, fullyDecoratedGraph, blockDisplayV, preOrderBlockVect, inCharInfoVV)
        --)
        
-- | doBlockTraversal takes a block of postorder decorated character trees character info  
-- could be moved up preOrderTreeTraversal, but like this for legibility
doBlockTraversal :: (V.Vector CharInfo, V.Vector DecoratedGraph) -> V.Vector DecoratedGraph
doBlockTraversal (inCharInfoV, traversalDecoratedVect) =
    --trace ("BlockT:" ++ (show $ fmap charType inCharInfoV)) 
    fmap doCharacterTraversal $ Debug.debugVectorZip inCharInfoV traversalDecoratedVect

-- | doCharacterTraversal perfoms preorder traversal on single character tree
-- with single charInfo
-- this so each character can be independently "rooted" for optimal traversals.
doCharacterTraversal ::  (CharInfo, DecoratedGraph) -> DecoratedGraph 
doCharacterTraversal (inCharInfo, inGraph) =
    -- find root--index should = number of leaves 
    --trace ("charT:" ++ (show $ charType inCharInfo)) (
    let rootVertexList = LG.getRoots inGraph
        (_, leafVertexList, _, _)  = LG.splitVertexList inGraph
        rootIndex = fst $ head rootVertexList
        inEdgeList = LG.labEdges inGraph
    in
    -- remove these two lines if working
    if length rootVertexList /= 1 then error ("Root number not = 1 in doCharacterTraversal" ++ show (rootVertexList))
    else if rootIndex /=  length leafVertexList then error ("Root index not =  number leaves in doCharacterTraversal" ++ show (rootIndex, length rootVertexList))
    else 
        -- root vertex, repeat of label info to avoid problem with zero length zip later, second info ignored for root
        let rootLabel = snd $ head rootVertexList
            rootFinalVertData = PRE.createFinalAssignmentOverBlocks RootNode (vertData rootLabel) (vertData rootLabel) inCharInfo True
            rootChildren =LG.labDescendants inGraph (head rootVertexList)

            -- this assumes two children
            rootChildrenBV = fmap bvLabel $ fmap snd rootChildren
            rootChildrenIsLeft = if (rootChildrenBV !! 0) > (rootChildrenBV !! 1) then [False, True]
                                 else [True, False]
            newRootNode = (rootIndex, rootLabel {vertData = rootFinalVertData})
            rootChildrenPairs = zip3 rootChildren (replicate (length rootChildren) newRootNode) rootChildrenIsLeft
            upDatedNodes = makeFinalAndChildren inGraph rootChildrenPairs [newRootNode] inCharInfo
        in
        -- hope this is the most efficient way since all nodes have been remade
        --trace (U.prettyPrintVertexInfo $ rootLabel {vertData = rootFinalVertData})
        LG.mkGraph upDatedNodes inEdgeList
        --)

-- | makeFinalAndChildren takes a graph, list of pairs of (labelled nodes,parent node) to make final assignment and a liss of updated nodes
-- the input nodes are relabelled by preorder functions and added to the list of processed nodes and recursed to their children
-- nodes are retuned in reverse order at they are made--need to check if this will affect graph identity or indexing in fgl
makeFinalAndChildren :: DecoratedGraph 
                     -> [(LG.LNode VertexInfo, LG.LNode VertexInfo, Bool)] 
                     -> [LG.LNode VertexInfo] 
                     -> CharInfo 
                     -> [LG.LNode VertexInfo]
makeFinalAndChildren inGraph nodesToUpdate updatedNodes inCharInfo =
    --trace ("mFAC:" ++ (show $ charType inCharInfo)) (
    if null nodesToUpdate then updatedNodes
    else 
        let (firstNode, firstParent, isLeft) = head nodesToUpdate
            firstLabel = snd firstNode
            firstNodeType = nodeType firstLabel
            firstVertData = vertData firstLabel
            firstParentVertData = vertData $ snd firstParent
            firstChildren = LG.labDescendants inGraph firstNode
            -- this assumes two children
            firstChildrenBV = fmap bvLabel $ fmap snd firstChildren
            firstChildrenIsLeft = if (firstChildrenBV !! 0) > (firstChildrenBV !! 1) then [False, True]
                                 else [True, False]
            firstFinalVertData = PRE.createFinalAssignmentOverBlocks firstNodeType firstVertData firstParentVertData inCharInfo isLeft
            newFirstNode = (fst firstNode, firstLabel {vertData = firstFinalVertData})
            childrenPairs = zip3 firstChildren (replicate (length firstChildren) newFirstNode) firstChildrenIsLeft
        in
        --trace (U.prettyPrintVertexInfo $ firstLabel {vertData = firstFinalVertData})
        makeFinalAndChildren inGraph (childrenPairs ++ (tail nodesToUpdate)) (newFirstNode : updatedNodes) inCharInfo
        --)

-- | assignPreorderStatesAndEdges takes a postorder decorated graph (should be but not required) and propagates 
-- preorder character states from individual character trees.  Exact characters (Add, nonAdd, matrix) postorder
-- states should be based on the outgroup rooted tree.  
-- root should be median of finals of two descendets--for non-exact based on final 'alignments' field with gaps filtered
-- postorder assignment and preorder will be out of whack--could change to update with correponding postorder
-- but that would not allow use of base decorated graph for incremental optimization (which relies on postorder assignments) in other areas
-- optyion code ikn there to set root final to outgropu final--but makes thigs scewey in matrix character and some pre-order assumptions
assignPreorderStatesAndEdges :: V.Vector (V.Vector DecoratedGraph) -> V.Vector (V.Vector CharInfo) -> DecoratedGraph  -> DecoratedGraph
assignPreorderStatesAndEdges preOrderBlockTreeVV inCharInfoVV inGraph =
    --trace ("aPSAE:" ++ (show $ fmap (fmap charType) inCharInfoVV)) (
    if LG.isEmpty inGraph then error "Empty graph in assignPreorderStatesAndEdges"
    else 
        -- trace ("In assign") (
        let postOrderNodes = LG.labNodes inGraph
            postOrderEdgeList = LG.labEdges inGraph

            -- update node labels
            newNodeList = fmap (updateNodeWithPreorder preOrderBlockTreeVV inCharInfoVV) postOrderNodes
            
            {-This code will set root final (and preliminary as well) to outgroup final
            --  remake root from two descendants-- assign one of two descendants
            -- if one is a terminal that then left
            rootNode = head $ LG.getRoots inGraph
            rootChildrenIndices = LG.descendants inGraph (fst rootNode)
            rootChildren = filter ((`elem` rootChildrenIndices).fst) newNodeList
            newRootLabel = makeNewRootLabel rootChildren graphCost (snd rootNode) 

            -- remove old root node and add new
            newNodeList' = (fst rootNode, newRootLabel) : (filter ((/= (fst rootNode)).fst) newNodeList)
            
            -- update edge labels
            newEdgeList = fmap (updateEdgeInfo inCharInfoVV (V.fromList $ L.sortOn fst newNodeList')) postOrderEdgeList
            -}
            newEdgeList = fmap (updateEdgeInfo inCharInfoVV (V.fromList $ L.sortOn fst newNodeList)) postOrderEdgeList
        in
        -- make new graph
        -- LG.mkGraph newNodeList' newEdgeList
        LG.mkGraph newNodeList newEdgeList
        --)

{-
-- | makeNewRootLabel takes list of labeled nodes (children of root) and sets te hroot data to the leaf data if one is a leaf
-- otherwise the vertData of the first child.  This to avoid root branch length problems created by preorder pass over travesals.
-- This will result in one edge being zero length--more properly reflecting the edge-as-root idea.
-- need to update only final states.  
makeNewRootLabel :: [LG.LNode VertexInfo] -> VertexCost -> VertexInfo -> VertexInfo
makeNewRootLabel rootChildren graphCost inRootInfo =
    if null rootChildren then error "Null children of root in chooseNewRootLabel"
    else 
        let (_, childLeafList) = unzip $ filter ((==LeafNode).fst) $ zip (fmap nodeType $ fmap snd rootChildren) rootChildren
        in
        if not $ null childLeafList then 
            let childLeafData = snd $ head childLeafList
            in
            -- trace (show $ vertData childLeafData)
            inRootInfo { vertData = vertData childLeafData
                       , vertexCost = 0.0
                       , subGraphCost = graphCost
                       }
        else 
            let childLeafData = snd $ head rootChildren
            in
            inRootInfo { vertData = vertData childLeafData
                       , vertexCost = 0.0
                       , subGraphCost = graphCost
                       }
-}

-- | updateNodeWithPreorder takes the preorder decorated graphs (by block and character) and updates the
-- the preorder fields only using character info.  This leaves post and preorder assignment out of sync.
-- but that so can use incremental optimizaytion on base decorated graph in other areas.
updateNodeWithPreorder :: V.Vector (V.Vector DecoratedGraph) -> V.Vector (V.Vector CharInfo) -> LG.LNode VertexInfo -> LG.LNode VertexInfo
updateNodeWithPreorder preOrderBlockTreeVV inCharInfoVV postOrderNode =
    let nodeLabel = snd postOrderNode
        nodeVertData = vertData nodeLabel
        newNodeVertData = fmap (updateVertexBlock (fst postOrderNode))$ D.debugVectorZip3 preOrderBlockTreeVV nodeVertData inCharInfoVV 
    in
    (fst postOrderNode, nodeLabel {vertData = newNodeVertData})

-- | updateVertexBlock takes a block of vertex data and updates preorder states of charactes via fmap
updateVertexBlock :: Int -> (V.Vector DecoratedGraph, V.Vector CharacterData, V.Vector CharInfo) -> V.Vector CharacterData 
updateVertexBlock nodeIndex (blockTraversalTreeV, nodeCharacterDataV, charInfoV) =
    fmap (updatePreorderCharacter nodeIndex) $ D.debugVectorZip3 blockTraversalTreeV nodeCharacterDataV charInfoV

-- | updatePreorderCharacter updates the pre-order fields of character data for a vertex from a traversal
-- since there is single character optimized for each character decorated graph-- it is always teh 0th 0th character
-- exact are vectors so take care of multipl there.
updatePreorderCharacter :: Int -> (DecoratedGraph, CharacterData, CharInfo) -> CharacterData 
updatePreorderCharacter nodeIndex (preOrderTree, postOrderCharacter, charInfo) =
    --trace ("N:" ++ (show nodeIndex) ++ " B:" ++ (show blockIndex) ++ " C:" ++ (show characterIndex) ++ "\n" ++ (show $ vertData $ fromJust $ LG.lab preOrderTree nodeIndex)) (
    let maybePreOrderNodeLabel = LG.lab preOrderTree nodeIndex
        preOrderVertData = vertData $ fromJust maybePreOrderNodeLabel
        preOrderCharacterData = if V.null preOrderVertData then emptyCharacter
                                else if V.null $ V.head preOrderVertData then emptyCharacter
                                else V.head $ V.head preOrderVertData -- (preOrderVertData V.! 0) V.! 0
    in
    if maybePreOrderNodeLabel == Nothing then error ("Nothing node label in updatePreorderCharacter node: " ++ show nodeIndex)
    else
        updateCharacter postOrderCharacter preOrderCharacterData (charType charInfo)
    --)

-- | updateCharacter takes a postorder character and updates the preorder (final) fields with preorder data and character type
updateCharacter :: CharacterData -> CharacterData -> CharType -> CharacterData
updateCharacter postOrderCharacter preOrderCharacter localCharType =
    if localCharType == Add then
        postOrderCharacter { rangeFinal = rangeFinal preOrderCharacter }
    
    else if localCharType == NonAdd then
        postOrderCharacter { stateBVFinal = stateBVFinal preOrderCharacter }
    
    else if localCharType == Matrix then
        postOrderCharacter { matrixStatesFinal = matrixStatesFinal preOrderCharacter }
    
    else if (localCharType == SlimSeq || localCharType == NucSeq) then
        postOrderCharacter { slimAlignment = slimAlignment preOrderCharacter
                           , slimFinal = slimFinal preOrderCharacter
                           }
    
    else if (localCharType == WideSeq || localCharType == AminoSeq) then
        postOrderCharacter { wideAlignment = wideAlignment preOrderCharacter
                           , wideFinal = wideFinal preOrderCharacter
                           }

    else if localCharType == HugeSeq then
        postOrderCharacter { hugeAlignment = hugeAlignment preOrderCharacter
                           , hugeFinal = hugeFinal preOrderCharacter
                           }

    else error ("Character type unimplemented : " ++ show localCharType)


-- | updateEdgeInfo takes a Decorated graph--fully labelled post and preorder and and edge and 
-- gets edge info--basically lengths
updateEdgeInfo :: V.Vector (V.Vector CharInfo) -> V.Vector (LG.LNode VertexInfo) -> LG.LEdge EdgeInfo -> LG.LEdge EdgeInfo 
updateEdgeInfo  inCharInfoVV nodeVector (uNode, vNode, edgeLabel) =
    if V.null nodeVector then error "Empty node list in updateEdgeInfo"
    else 
        let (minW, maxW) = getEdgeWeight inCharInfoVV nodeVector (uNode, vNode)
            midW = (minW + maxW) / 2.0
            localEdgeType = edgeType edgeLabel
            newEdgeLabel = EdgeInfo { minLength = minW
                                    , maxLength = maxW
                                    , midRangeLength = midW
                                    , edgeType  = localEdgeType
                                    }
        in
        (uNode, vNode, newEdgeLabel)
        
-- | getEdgeWeight takes a preorder decorated decorated graph and an edge and gets the weight information for that edge
-- basically a min/max distance between the two
getEdgeWeight ::   V.Vector (V.Vector CharInfo) -> V.Vector (LG.LNode VertexInfo) -> (Int, Int) -> (VertexCost, VertexCost)
getEdgeWeight inCharInfoVV nodeVector (uNode, vNode) = 
    if V.null nodeVector then error "Empty node list in getEdgeWeight"
    else 
        let uNodeInfo = vertData $ snd $ nodeVector V.! uNode
            vNodeInfo = vertData $ snd $ nodeVector V.! vNode
            blockCostPairs = fmap getBlockCostPairs $ D.debugVectorZip3 uNodeInfo vNodeInfo inCharInfoVV
            minCost = sum $ fmap fst blockCostPairs
            maxCost = sum $ fmap snd blockCostPairs
        in
        (minCost, maxCost)
        
-- | getBlockCostPairs takes a block of two nodes and character infomation and returns the min and max block branch costs
getBlockCostPairs :: (V.Vector CharacterData, V.Vector CharacterData, V.Vector CharInfo) -> (VertexCost, VertexCost)
getBlockCostPairs (uNodeCharDataV, vNodeCharDataV, charInfoV) = 
    let characterCostPairs = fmap getCharacterDist $ D.debugVectorZip3 uNodeCharDataV vNodeCharDataV charInfoV
        minCost = sum $ fmap fst characterCostPairs
        maxCost = sum $ fmap snd characterCostPairs
    in
    (minCost, maxCost)   

-- | getCharacterDist takes a pair of characters and character type, retunring teh minimum and maximum character distances
-- for sequence charcaters this is based on slim/wide/hugeAlignment field, hence all should be n in num characters/seqeunce length
getCharacterDist :: (CharacterData, CharacterData, CharInfo) -> (VertexCost, VertexCost)
getCharacterDist (uCharacter, vCharacter, charInfo) =
    let thisWeight = weight charInfo
        thisMatrix = costMatrix charInfo
        thisCharType = charType charInfo
    in
    if thisCharType == Add then
        let minCost = localCost (M.intervalAdd thisWeight uCharacter vCharacter)
            maxDiff = V.sum $ fmap maxIntervalDiff $ D.debugVectorZip (rangeFinal uCharacter) (rangeFinal vCharacter)
            maxCost = thisWeight * (fromIntegral maxDiff)
        in
        (minCost, maxCost)

        
    else if thisCharType == NonAdd then
        let minCost = localCost (M.interUnion thisWeight uCharacter vCharacter)
            maxDiff = length $ V.filter (==False) $ V.zipWith (==) (stateBVFinal uCharacter) (stateBVFinal vCharacter)
            maxCost = thisWeight * (fromIntegral maxDiff)
        in
        (minCost, maxCost)
    
    else if thisCharType == Matrix then
        let minMaxListList= fmap (minMaxMatrixDiff thisMatrix) $ D.debugVectorZip (fmap (fmap fst3) $ matrixStatesFinal uCharacter) (fmap (fmap fst3) $ matrixStatesFinal vCharacter)
            minDiff = V.sum $ fmap fst minMaxListList
            maxDiff = V.sum $ fmap snd minMaxListList
            minCost = thisWeight * (fromIntegral minDiff)
            maxCost = thisWeight * (fromIntegral maxDiff)
        in
        (minCost, maxCost)


    else if (thisCharType == SlimSeq || thisCharType == NucSeq) then
        let minMaxDiffList = fmap (generalSequenceDiff thisMatrix (length thisMatrix)) $ zip (GV.toList $ fst3 $ slimAlignment uCharacter) (GV.toList $ fst3 $ slimAlignment vCharacter)
            (minDiff, maxDiff) = unzip minMaxDiffList
            minCost = thisWeight * (fromIntegral $ sum minDiff)
            maxCost = thisWeight * (fromIntegral $ sum maxDiff)
        in
        (minCost, maxCost)
        
    
    else if (thisCharType == WideSeq || thisCharType == AminoSeq) then
        let minMaxDiffList = GV.toList $ GV.map (generalSequenceDiff thisMatrix (length thisMatrix)) $ GV.zip (fst3 $ wideAlignment uCharacter) (fst3 $ wideAlignment vCharacter)
            (minDiff, maxDiff) = unzip minMaxDiffList
            minCost = thisWeight * (fromIntegral $ sum minDiff)
            maxCost = thisWeight * (fromIntegral $ sum maxDiff)
        in
        (minCost, maxCost)

    else if thisCharType == HugeSeq then
        let minMaxDiffList = GV.toList $ GV.map (generalSequenceDiff thisMatrix (length thisMatrix)) $ GV.zip (fst3 $ hugeAlignment uCharacter) (fst3 $ hugeAlignment vCharacter)
            (minDiff, maxDiff) = unzip minMaxDiffList
            minCost = thisWeight * (fromIntegral $ sum minDiff)
            maxCost = thisWeight * (fromIntegral $ sum maxDiff)
        in
        (minCost, maxCost)

    else error ("Character type unimplemented : " ++ show thisCharType)

-- | maxIntervalDiff takes two ranges and gets the maximum difference between the two based on differences
-- in upp and lower ranges.
maxIntervalDiff :: ((Int, Int), (Int, Int)) -> Int 
maxIntervalDiff ((a,b), (x,y)) =
    let upper = (max b y) - (min b y)
        lower = (max a x) - (min a x)
    in
    max upper lower 

-- | minMaxMatrixDiff takes twovetors of states and calculates the minimum and maximum state differnce cost 
-- between the two
minMaxMatrixDiff :: S.Matrix Int -> (V.Vector Int, V.Vector Int) -> (Int, Int)
minMaxMatrixDiff localCostMatrix (uStatesV, vStatesV) =
    let statePairs = (V.toList uStatesV, V.toList vStatesV)
        cartesianPairs = cartProdPair statePairs
        costList = fmap (localCostMatrix S.!) cartesianPairs
    in 
    --trace (show cartesianPairs  ++ " " ++ show costList) 
    (minimum costList, maximum costList)
   


-- | generalSequenceDiff  takes two sequnce elemental bit types and retuns min and max integer 
-- cost differences using matrix values
generalSequenceDiff :: (FiniteBits a) => S.Matrix Int -> Int -> (a, a) -> (Int, Int)
generalSequenceDiff thisMatrix numStates (uState, vState) = 
    let uStateList = fmap snd $ filter ((== True).fst) $ zip (fmap (testBit uState) [0.. numStates - 1]) [0.. numStates - 1]
        vStateList = fmap snd $ filter ((== True).fst) $ zip (fmap (testBit vState) [0.. numStates - 1]) [0.. numStates - 1]    
        uvCombinations = cartProd uStateList vStateList
        costOfPairs = fmap (thisMatrix S.!) uvCombinations
    in
    (minimum costOfPairs, maximum costOfPairs)

