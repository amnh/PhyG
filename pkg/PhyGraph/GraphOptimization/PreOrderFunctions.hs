{- |
Module      :  PostOrderFunctions.hs
Description :  Module specifying pre-order graph functions
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

module GraphOptimization.PreOrderFunctions  ( createFinalAssignmentOverBlocks
                                            , preOrderTreeTraversal
                                            , getBlockCostPairsFinal
                                            , setFinalToPreliminaryStates
                                            , setPreliminaryToFinalStates
                                            , zero2Gap
                                            ) where

import           Bio.DynamicCharacter
import           Data.Bits
import qualified Data.BitVector.LittleEndian as BV
import qualified Data.List                   as L
import qualified Data.Map                    as MAP
import           Data.Maybe
import qualified Data.Vector                 as V
import qualified Data.Vector.Generic         as GV
-- import qualified Data.Vector.Storable         as SV
import qualified Data.Vector.Unboxed         as UV
import qualified DirectOptimization.PreOrder as DOP
import           GeneralUtilities
import qualified GraphOptimization.Medians   as M
import qualified Graphs.GraphOperations      as GO
import qualified Input.BitPack               as BP
import qualified SymMatrix                   as S
import           Types.Types
import qualified Utilities.LocalGraph        as LG
import qualified Utilities.ThreeWayFunctions as TW
import qualified Utilities.Utilities         as U
import           Data.Alphabet
import           Debug.Trace
import           Control.Parallel.Strategies
import qualified ParallelUtilities                    as PU


-- | preOrderTreeTraversal takes a preliminarily labelled PhylogeneticGraph
-- and returns a full labels with 'final' assignments based on character decorated graphs
-- created postorder (5th of 6 fields).
-- the preorder states are created by traversing the traversal DecoratedGraphs in the 5th filed of PhylogeneticGraphs
-- these are by block and character,  Exact characters are vectors of standard characters and each sequence (non-exact)
-- has its own traversal graph. These should be trees here (could be forests later) and should all have same root (number taxa)
-- but worth checking to make sure.
-- these were created by "splitting" them after preorder, or by separate passes in softwired graphs

-- For sequence characters (slim/wide/huge) final states are created either by DirectOptimization or ImpliedAlignment
-- if DO--then does a  median between parent andchgaps wherre needed--then  doing a 3way state assignmeent filteringgaps out
-- if IA--then a separate post and preorder pass are donne on the slim/wide/huge/AI fields to crete full IA assignments
-- that are then filtered of gaps and assigned to th efinal fields

-- The final states are propagated back to the second field
-- DecoratedGraph of the full Phylogenetic graph--which does NOT have the character based preliminary assignments
-- ie postorder--since those are traversal specific
-- the character specific decorated graphs have appropriate post and pre-order assignments
-- the traversal begins at the root (for a tree) and proceeds to leaves.
preOrderTreeTraversal :: GlobalSettings -> AssignmentMethod -> Bool -> Bool -> Bool -> Int -> Bool -> PhylogeneticGraph -> PhylogeneticGraph
preOrderTreeTraversal inGS finalMethod staticIA calculateBranchLengths _ rootIndex useMap (inSimple, inCost, inDecorated, blockDisplayV, blockCharacterDecoratedVV, inCharInfoVV) =
    --trace ("PreO: " ++ (show finalMethod) ++ " " ++ (show $ fmap (fmap charType) inCharInfoVV)) (
    -- trace ("PR-OT pre: " ++ (show $ fmap V.length blockCharacterDecoratedVV)) (
    if LG.isEmpty inDecorated then error "Empty tree in preOrderTreeTraversal"
    else
        -- trace ("In PreOrder\n" ++ "Simple:\n" ++ (LG.prettify inSimple) ++ "Decorated:\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inDecorated) ++ "\n" ++ (GFU.showGraph inDecorated)) (
        -- mapped recursive call over blkocks, later characters
        let preOrderBlockVect = V.fromList (PU.seqParMap rdeepseq  (doBlockTraversal' inGS finalMethod staticIA rootIndex) (zip (V.toList inCharInfoVV) (V.toList blockCharacterDecoratedVV)))  -- `using` PU.myParListChunkRDS)

            -- if final non-exact states determined by IA then perform passes and assignments of final and final IA fields
            -- always do IA pass if Tree--but only assign to final if finalMethod == ImpliedAlignment
            -- also assignes unions for use in rearrangemenrts
            -- should be ok for softwired--if done by display trees, and soft and hardwired with outdegree=1
            preOrderBlockVect' = V.zipWith (makeIAUnionAssignments finalMethod rootIndex) preOrderBlockVect inCharInfoVV
                                 -- if hasNonExact && (graphType inGS) == Tree then V.zipWith (makeIAAssignments finalMethod rootIndex) preOrderBlockVect inCharInfoVV
                                 -- else preOrderBlockVect

            fullyDecoratedGraph = assignPreorderStatesAndEdges inGS finalMethod calculateBranchLengths rootIndex preOrderBlockVect' useMap inCharInfoVV inDecorated
        in
        if null blockCharacterDecoratedVV then error ("Empty preOrderBlockVect in preOrderTreeTraversal at root index rootIndex: " ++ (show rootIndex) ++ " This can be caused if the graphType not set correctly: " ++ (show $ graphType inGS))
        else
            {-
            let blockPost = GO.showDecGraphs blockCharacterDecoratedVV
                blockPre = GO.showDecGraphs preOrderBlockVect
            in
            trace ("BlockPost:\n" ++ blockPost ++ "BlockPre:\n" ++ blockPre ++ "After Preorder\n" ++  (LG.prettify $ GO.convertDecoratedToSimpleGraph fullyDecoratedGraph))
            -}
            (inSimple, inCost, fullyDecoratedGraph, blockDisplayV, preOrderBlockVect, inCharInfoVV)
    -- )

-- | makeIAUnionAssignments takes the vector of vector of character trees and (if) slim/wide/huge
-- does an additional post and pre order pass to assign IAand final fields in all sequece types slim/wide/huge
-- assigns union for all types
makeIAUnionAssignments :: AssignmentMethod -> Int -> V.Vector DecoratedGraph -> V.Vector CharInfo -> V.Vector DecoratedGraph
makeIAUnionAssignments finalMethod rootIndex = V.zipWith (makeCharacterIAUnion finalMethod rootIndex)

-- | makeCharacterIAUnion takes an individual character postorder tree and if perform post and preorder IA passes
-- and assignment to final field in slim/wide/huge
-- also assignes unions for all types
makeCharacterIAUnion :: AssignmentMethod -> Int -> DecoratedGraph -> CharInfo -> DecoratedGraph
makeCharacterIAUnion finalMethod rootIndex inGraph charInfo =
    -- if charType charInfo `notElem` nonExactCharacterTypes then inGraph
    if False then inGraph
    else
        let postOrderIATree = postOrderIAUnion inGraph charInfo [(rootIndex, fromJust $ LG.lab inGraph rootIndex)]
            preOrderIATree = preOrderIA postOrderIATree rootIndex finalMethod charInfo $ zip [(rootIndex, fromJust $ LG.lab postOrderIATree rootIndex)] [(rootIndex, fromJust $ LG.lab postOrderIATree rootIndex)]
        in
        -- trace ("MCIAU:" ++ (U.getUnionFieldsNode $ vertData $ fromJust $ LG.lab postOrderIATree rootIndex) ++ "\n" ++ (U.getUnionFieldsNode $ vertData $ fromJust $ LG.lab postOrderIATree 0) 
        --    ++ "\nAfter preorder:\t" ++ (U.getUnionFieldsNode $ vertData $ fromJust $ LG.lab preOrderIATree rootIndex) ++ "\n" ++ (U.getUnionFieldsNode $ vertData $ fromJust $ LG.lab preOrderIATree 0))
        preOrderIATree

-- | postOrderIAUnion performs a post-order IA pass assigning leaf preliminary states
-- from the "alignment" fields and setting HTU preliminary by calling the apropriate 2-way
-- matrix
-- should eb OK for any root--or partial graph as in for branch swapping
-- also sets unions for all character types, IA based for sequence
postOrderIAUnion :: DecoratedGraph -> CharInfo -> [LG.LNode VertexInfo] -> DecoratedGraph
postOrderIAUnion inGraph charInfo inNodeList =
    if null inNodeList then inGraph
    else
        let inNode@(nodeIndex, nodeLabel) = head inNodeList
            (inNodeEdges, outNodeEdges) = LG.getInOutEdges inGraph nodeIndex
            inCharacter = V.head $ V.head $ vertData nodeLabel
            nodeType' = GO.getNodeType inGraph nodeIndex
        in

        -- checking sanity of data
        if V.null $ vertData nodeLabel then error "Null vertData in postOrderIA"
        else if V.null $ V.head $ vertData nodeLabel then
            -- missing data for taxon
            error "Null vertData data in postOrderIA"

        -- leaf take assignment from alignment field
        else if nodeType' == LeafNode then 
            -- set leaf union fields to preliminary or IA fields
            let newCharacter = M.makeIAUnionPrelimLeaf charInfo inCharacter 
                newLabel = nodeLabel  {vertData = V.singleton (V.singleton newCharacter), nodeType = nodeType'}
                newGraph = LG.insEdges (inNodeEdges ++ outNodeEdges) $ LG.insNode (nodeIndex, newLabel) $ LG.delNode nodeIndex inGraph
            in
            postOrderIAUnion newGraph charInfo (tail inNodeList)

        -- HTU take create assignment from children
        else
            let childNodes = LG.labDescendants inGraph inNode
                childTree = postOrderIAUnion inGraph charInfo childNodes
            in
            --trace ("Children: " ++ (show  $ fmap fst childNodes)) (

            if length childNodes > 2 then error ("Too many children in postOrderIA: " ++ show (length childNodes))

            -- in 1 out 1 vertex
            else if length childNodes == 1 then
                let childIndex = fst $ head childNodes
                    childLabel = fromJust $ LG.lab childTree childIndex
                    childCharacter = V.head $ V.head $ vertData childLabel
                in
                -- sanity checks
                if isNothing (LG.lab childTree (fst $ head childNodes)) then error ("No label for node: " ++ show (fst $ head childNodes))
                else if V.null $ vertData childLabel then error "Null vertData in postOrderIA"
                else if V.null $ V.head $ vertData childLabel then error "Null head vertData data in postOrderIA"
                else
                    let newLabel = nodeLabel  {vertData = V.singleton (V.singleton childCharacter), nodeType = nodeType'}
                        newGraph = LG.insEdges (inNodeEdges ++ outNodeEdges) $ LG.insNode (nodeIndex, newLabel) $ LG.delNode nodeIndex childTree
                    in
                    -- trace ("PostO1Child: " ++ (show nodeIndex) ++ " " ++ (show $ slimFinal childCharacter))
                    postOrderIAUnion newGraph charInfo (tail inNodeList)

            -- two children
            else
                let childIndices = fmap fst childNodes
                    childlabels = fmap (fromJust . LG.lab childTree) childIndices
                    childCharacters = fmap vertData childlabels
                    leftChar = V.head $ V.head $ head childCharacters
                    rightChar = V.head $ V.head $ last childCharacters
                    newCharacter = M.makeIAPrelimCharacter charInfo inCharacter leftChar rightChar
                    newLabel = nodeLabel  {vertData = V.singleton (V.singleton newCharacter), nodeType = nodeType'}
                    newGraph = LG.insEdges (inNodeEdges ++ outNodeEdges) $ LG.insNode (nodeIndex, newLabel) $ LG.delNode nodeIndex childTree
                in
                -- trace ("PostO2hildren: " ++ (show nodeIndex) ++ " " ++ (show $ slimFinal newCharacter) ++ " " ++ (show $ nodeType newLabel)) -- ++ " From: " ++ (show childlabels))
                postOrderIAUnion newGraph charInfo (tail inNodeList)
            -- )
    -- )




-- | preOrderIA performs a pre-order IA pass assigning via the apropriate 3-way matrix
-- the "final" fields are also set by filtering out gaps and 0.
-- skips non-unaligned-sequence types
preOrderIA :: DecoratedGraph -> Int -> AssignmentMethod -> CharInfo -> [(LG.LNode VertexInfo, LG.LNode VertexInfo)] -> DecoratedGraph
preOrderIA inGraph rootIndex finalMethod charInfo inNodePairList =
    if null inNodePairList then inGraph
    else
        let (inNode@(nodeIndex, nodeLabel), (_, parentNodeLabel)) = head inNodePairList
            (inNodeEdges, outNodeEdges) = LG.getInOutEdges inGraph nodeIndex
            characterType = charType charInfo
            inCharacter = V.head $ V.head $ vertData nodeLabel
            inCharacter' = inCharacter
            parentCharacter = V.head $ V.head $ vertData parentNodeLabel
            childNodes = LG.labDescendants inGraph inNode
        in
        --trace ("PreIA Node:" ++ (show nodeIndex) ++ " " ++ (show $ nodeType nodeLabel) ++ " " ++ (show (fmap fst $ fmap fst inNodePairList,fmap fst $ fmap snd inNodePairList))) (
        -- checking sanity of data
        if V.null $ vertData nodeLabel then error "Null vertData in preOrderIA"
        else if V.null $ V.head $ vertData nodeLabel then error "Null vertData data in preOrderIA"
        else if length childNodes > 2 then error ("Too many children in preOrderIA: " ++ show (length childNodes))

        -- leaf done in post-order
        else if nodeType nodeLabel == LeafNode then preOrderIA inGraph rootIndex finalMethod charInfo (tail inNodePairList)

        else if nodeType nodeLabel == RootNode || nodeIndex == rootIndex then
            let newCharacter
                  | characterType `elem` [SlimSeq, NucSeq] =
                     inCharacter { slimIAFinal = extractMediansGapped $ slimIAPrelim inCharacter'
                                 , slimFinal = extractMedians $ slimGapped  inCharacter'
                                 }
                  | characterType `elem` [WideSeq, AminoSeq] =
                     inCharacter { wideIAFinal = extractMediansGapped $ wideIAPrelim inCharacter'
                                 , wideFinal = extractMedians $ wideGapped inCharacter'
                                 }
                  | characterType == HugeSeq =
                     inCharacter { hugeIAFinal = extractMediansGapped $ hugeIAPrelim inCharacter'
                                 , hugeFinal = extractMedians $ hugeGapped inCharacter'
                                 }
                  | otherwise = inCharacter -- error ("Unrecognized character type " ++ show characterType)

                newLabel = nodeLabel  {vertData = V.singleton (V.singleton newCharacter)}
                newGraph = LG.insEdges (inNodeEdges ++ outNodeEdges) $ LG.insNode (nodeIndex, newLabel) $ LG.delNode nodeIndex inGraph
                parentNodeList = replicate (length childNodes) (nodeIndex, newLabel)
            in
            -- trace ("PreIARoot: " ++ (show nodeIndex) ++ " IAFinal: " ++ (show $ slimIAFinal newCharacter) ++ " Final: " ++ (show $ slimFinal newCharacter))
            preOrderIA newGraph rootIndex finalMethod charInfo (tail inNodePairList ++ zip childNodes parentNodeList)

        -- single child, take parent final assignments, but keep postorder assignments
        else if length childNodes == 1 then
                let newCharacter
                      | characterType `elem` [SlimSeq, NucSeq] =
                       inCharacter { slimFinal = slimFinal parentCharacter
                                   , slimIAFinal = slimIAFinal parentCharacter
                                   }
                      | characterType `elem` [WideSeq, AminoSeq] =
                        inCharacter { wideFinal = wideFinal parentCharacter
                                    , wideIAFinal = wideIAFinal parentCharacter
                                    }
                      | characterType == HugeSeq =
                        inCharacter { hugeFinal = hugeFinal parentCharacter
                                    , hugeIAFinal = hugeIAFinal parentCharacter
                                    }
                      | otherwise = inCharacter -- error ("Unrecognized character type " ++ show characterType)

                    newLabel = nodeLabel  {vertData = V.singleton (V.singleton newCharacter)}
                    newGraph = LG.insEdges (inNodeEdges ++ outNodeEdges) $ LG.insNode (nodeIndex, newLabel) $ LG.delNode nodeIndex inGraph
                    parentNodeList = replicate (length childNodes) (nodeIndex, newLabel)
                in
                -- trace ("PreIANet: " ++ (show nodeIndex) ++ " IAFinal: " ++ (show $ slimIAFinal newCharacter) ++ " Final: " ++ (show $ slimFinal newCharacter))
                preOrderIA newGraph rootIndex finalMethod charInfo (tail inNodePairList ++ zip childNodes parentNodeList)

        -- 2 children, make 3-way
        else
            let finalCharacter = M.makeIAFinalCharacter finalMethod charInfo inCharacter parentCharacter -- leftChar rightChar

                newLabel = nodeLabel  {vertData = V.singleton (V.singleton finalCharacter)}
                newGraph = LG.insEdges (inNodeEdges ++ outNodeEdges) $ LG.insNode (nodeIndex, newLabel) $ LG.delNode nodeIndex inGraph
                parentNodeList = replicate (length childNodes) (nodeIndex, newLabel)
            in
            -- trace ("PreIATree: " ++ (show nodeIndex) ++ " IAFinal: " ++ (show $ slimIAFinal finalCharacter) ++ " Final: " ++ (show $ slimFinal finalCharacter))
            preOrderIA newGraph rootIndex finalMethod charInfo (tail inNodePairList ++ zip childNodes parentNodeList)

        -- )

-- | doBlockTraversal' is a wrapper around doBlockTraversal fro seqParMap
doBlockTraversal' :: GlobalSettings -> AssignmentMethod -> Bool -> Int -> (V.Vector CharInfo, V.Vector DecoratedGraph) -> V.Vector DecoratedGraph
doBlockTraversal' inGS finalMethod staticIA rootIndex (inCharInfoV, traversalDecoratedVect) =
    doBlockTraversal inGS finalMethod staticIA rootIndex inCharInfoV traversalDecoratedVect

-- | doBlockTraversal takes a block of postorder decorated character trees character info
-- could be moved up preOrderTreeTraversal, but like this for legibility
doBlockTraversal :: GlobalSettings -> AssignmentMethod -> Bool -> Int -> V.Vector CharInfo -> V.Vector DecoratedGraph -> V.Vector DecoratedGraph
doBlockTraversal inGS finalMethod staticIA rootIndex inCharInfoV traversalDecoratedVect =
    --trace ("BlockT:" ++ (show $ fmap charType inCharInfoV))
    V.zipWith (doCharacterTraversal inGS finalMethod staticIA rootIndex) inCharInfoV traversalDecoratedVect

-- | doCharacterTraversal performs preorder traversal on single character tree
-- with single charInfo
-- this so each character can be independently "rooted" for optimal traversals.
doCharacterTraversal :: GlobalSettings -> AssignmentMethod -> Bool -> Int -> CharInfo -> DecoratedGraph -> DecoratedGraph
doCharacterTraversal inGS finalMethod staticIA rootIndex inCharInfo inGraph =
    -- find root--index should = number of leaves
    --trace ("charT:" ++ (show $ charType inCharInfo)) (
    let -- this is a hack--remve after fixed
        --inGraph = LG.removeDuplicateEdges inGraph'

        isolateNodeList = LG.getIsolatedNodes inGraph
        -- (_, leafVertexList, _, _)  = LG.splitVertexList inGraph
        inEdgeList = LG.labEdges inGraph
    in
    -- remove these two lines if working
    -- if rootIndex /=  length leafVertexList then error ("Root index not =  number leaves in doCharacterTraversal" ++ show (rootIndex, length leafVertexList))
    -- else
        -- root vertex, repeat of label info to avoid problem with zero length zip later, second info ignored for root
        -- since root cannot have 2nd parent
        let rootLabel = fromJust $ LG.lab inGraph rootIndex
            nothingVertData = U.copyToNothing (vertData rootLabel)
            rootFinalVertData = createFinalAssignmentOverBlocks inGS finalMethod staticIA RootNode (vertData rootLabel) (vertData rootLabel) nothingVertData inCharInfo True False False
            rootChildren =LG.labDescendants inGraph (rootIndex, rootLabel)

            -- left / right to match post-order
            rootChildrenBV = fmap (bvLabel . snd) rootChildren
            rootChildrenIsLeft
              | length rootChildrenBV == 1 = [True]
              | head rootChildrenBV > (rootChildrenBV !! 1) = [False, True]
              | otherwise = [True, False]
            newRootNode = (rootIndex, rootLabel {vertData = rootFinalVertData})
            rootChildrenPairs = zip3 rootChildren (replicate (length rootChildren) newRootNode) rootChildrenIsLeft
            upDatedNodes = makeFinalAndChildren inGS finalMethod staticIA inGraph rootChildrenPairs [newRootNode] inCharInfo

            -- update isolated nodes with final == preliminary as with root nodes (and leaves, but without postorder logic)
            updatedIsolateNodes = fmap (updateIsolatedNode inGS finalMethod staticIA inCharInfo) isolateNodeList
        in
        -- hope this is the most efficient way since all nodes have been remade
        -- trace (U.prettyPrintVertexInfo $ snd newRootNode)
        LG.mkGraph (upDatedNodes ++ updatedIsolateNodes) inEdgeList
        --)

-- | updateIsolatedNode updates the final states of an isolated node as if it were a root with final=preliminary
-- states without preorder logic as in regular leaves
-- NB IA length won't match if compared since not in graph
updateIsolatedNode :: GlobalSettings -> AssignmentMethod -> Bool -> CharInfo -> LG.LNode VertexInfo -> LG.LNode VertexInfo
updateIsolatedNode inGS finalMethod staticIA inCharInfo (inNodeIndex, inNodeLabel) =
    -- root so final = preliminary
    let nothingVertData = U.copyToNothing (vertData inNodeLabel)
        newVertData = createFinalAssignmentOverBlocks inGS finalMethod staticIA RootNode (vertData inNodeLabel) (vertData inNodeLabel) nothingVertData inCharInfo True False False
    in
    (inNodeIndex, inNodeLabel {vertData = newVertData})

-- | makeFinalAndChildren takes a graph, list of pairs of (labelled nodes,parent node) to make final assignment and a list of updated nodes
-- the input nodes are relabelled by preorder functions and added to the list of processed nodes and recursed to other nodes first, then
-- their children -- thi sis important for preorder of hardwired graphs since can have 2 parents a single child.
-- nodes are retuned in reverse order at they are made--need to check if this will affect graph identity or indexing in fgl
makeFinalAndChildren :: GlobalSettings
                     -> AssignmentMethod
                     -> Bool
                     -> DecoratedGraph
                     -> [(LG.LNode VertexInfo, LG.LNode VertexInfo, Bool)]
                     -> [LG.LNode VertexInfo]
                     -> CharInfo
                     -> [LG.LNode VertexInfo]
makeFinalAndChildren inGS finalMethod staticIA inGraph nodesToUpdate updatedNodes inCharInfo =
    --trace ("mFAC:" ++ (show $ charType inCharInfo)) (
    if null nodesToUpdate then updatedNodes
    else
        let (firstNode, firstParent, isLeft) = head nodesToUpdate

            -- get current node data
            firstLabel = snd firstNode
            firstNodeType' = GO.getNodeType inGraph $ fst firstNode -- nodeType firstLabel
            firstNodeType = if firstNodeType' /= NetworkNode then firstNodeType'
                            else
                                -- not issue if hardwired I don't think
                                if graphType inGS /= HardWired then trace ("NetNode:" ++ (show $ LG.getInOutDeg inGraph firstNode) ++ " DuplicateEdges (?): " ++ (show $ LG.getDuplicateEdges inGraph)) NetworkNode
                                else NetworkNode
            firstVertData = vertData firstLabel

            -- get node parent data--check if more than one
            firstParents = LG.labParents inGraph $ fst firstNode

            -- if single parent then as usual--else take head of two so no confusion as to whichn is which
            -- this is holdover from no indegree 2 nodes--could be simplified and return structures changed
            firstParentVertData = if (length firstParents == 1) then vertData $ snd firstParent
                                  else vertData $ snd $ head firstParents

            secondParentData = if (length firstParents == 1) then U.copyToNothing firstParentVertData
                               else U.copyToJust $ vertData $ snd $ last firstParents

            -- child data
            firstChildren = LG.labDescendants inGraph firstNode

            -- booleans for further pass
            isIn1Out1 = (length firstChildren == 1) && (length firstParents == 1) -- softwired can happen, need to pass "grandparent" node to skip in 1 out 1
            isIn2Out1 = (length firstChildren == 1) && (length firstParents == 2) -- hardwired can happen, need to pass both parents

            
            -- this OK with one or two children
            firstChildrenBV = fmap (bvLabel . snd) firstChildren
            firstChildrenIsLeft
              | length firstChildrenBV == 1 = [True]
              | head firstChildrenBV > (firstChildrenBV !! 1) = [False, True]
              | otherwise = [True, False]
            firstFinalVertData = createFinalAssignmentOverBlocks inGS finalMethod staticIA firstNodeType firstVertData firstParentVertData secondParentData inCharInfo isLeft isIn1Out1 isIn2Out1
            newFirstNode = (fst firstNode, firstLabel {vertData = firstFinalVertData})

            -- check children if indegree == 2 then don't add to nodes to do if in there already

            childrenTriple  = zip3 firstChildren (replicate (length firstChildren) newFirstNode) firstChildrenIsLeft
            childrenTriple' = if (graphType inGS) == HardWired then filter (indeg2NotInNodeList inGraph (tail nodesToUpdate)) childrenTriple
                              else childrenTriple

        in 
        -- trace (U.prettyPrintVertexInfo $ snd newFirstNode)
        -- makeFinalAndChildren inGS finalMethod staticIA inGraph (childrenPairs ++ tail nodesToUpdate) (newFirstNode : updatedNodes) inCharInfo
        -- childrenPair after nodess to do for hardWired to ensure both parent done before child
        makeFinalAndChildren inGS finalMethod staticIA inGraph ((tail nodesToUpdate) ++ childrenTriple') (newFirstNode : updatedNodes) inCharInfo
        --)

-- | indeg2NotInNodeList checcks a node agains a list by index (fst) if node is indegree 2 and
-- already in the list of n odes "todo" filter out as will already be optimized in appropriate pre-order
indeg2NotInNodeList :: LG.Gr a b -> [(LG.LNode a, LG.LNode a, Bool)] -> (LG.LNode a, LG.LNode a, Bool) -> Bool
indeg2NotInNodeList inGraph checkNodeList (childNode@(childIndex, _), _, _) =
 if LG.isEmpty inGraph then error "Empty graph in indeg2NotInNodeList"
 else
    if LG.indeg inGraph childNode < 2 then True
    else if childIndex `elem` (fmap (fst . fst3) checkNodeList) then False
    else True



-- | assignPreorderStatesAndEdges takes a postorder decorated graph (should be but not required) and propagates
-- preorder character states from individual character trees.  Exact characters (Add, nonAdd, matrix) postorder
-- states should be based on the outgroup rooted tree.
-- root should be median of finals of two descendets--for non-exact based on final 'alignments' field with gaps filtered
-- postorder assignment and preorder will be out of whack--could change to update with correponding postorder
-- but that would not allow use of base decorated graph for incremental optimization (which relies on postorder assignments) in other areas
-- optyion code ikn there to set root final to outgropu final--but makes thigs scewey in matrix character and some pre-order assumptions
assignPreorderStatesAndEdges :: GlobalSettings -> AssignmentMethod -> Bool -> Int -> V.Vector (V.Vector DecoratedGraph) -> Bool -> V.Vector (V.Vector CharInfo) -> DecoratedGraph  -> DecoratedGraph
assignPreorderStatesAndEdges inGS finalMethd calculateBranchEdges rootIndex preOrderBlockTreeVV useMap inCharInfoVV inGraph =
    --trace ("aPSAE:" ++ (show $ fmap (fmap charType) inCharInfoVV)) (
    if LG.isEmpty inGraph then error "Empty graph in assignPreorderStatesAndEdges"
    else
        -- trace ("In assign") (
        let postOrderNodes = LG.labNodes inGraph
            postOrderEdgeList = LG.labEdges inGraph

            -- update node labels
            newNodeList = fmap (updateNodeWithPreorder preOrderBlockTreeVV inCharInfoVV) postOrderNodes

            -- create a vector of vector of pair of nodes and edges for display x charcater trees
            blockTreePairVV =  fmap (fmap LG.makeNodeEdgePairVect) preOrderBlockTreeVV

            -- update edge labels--for softwired need to account for not all edges in all block/display trees
            -- map for case where tree does not contain all leaves as in swap procedures
            -- needs to be updated for softwired as well
            nodeMap = MAP.fromList $ zip (fmap fst newNodeList) newNodeList
            newEdgeList = if (graphType inGS == Tree || graphType inGS == HardWired) then
                                if useMap then fmap (updateEdgeInfoTreeMap finalMethd inCharInfoVV nodeMap) postOrderEdgeList
                                else fmap (updateEdgeInfoTree finalMethd inCharInfoVV (V.fromList newNodeList)) postOrderEdgeList
                          else
                                fmap (updateEdgeInfoSoftWired finalMethd inCharInfoVV blockTreePairVV rootIndex) postOrderEdgeList
        in
        -- make new graph
        -- LG.mkGraph newNodeList' newEdgeList
        if calculateBranchEdges then LG.mkGraph newNodeList newEdgeList
        else LG.mkGraph newNodeList postOrderEdgeList
        --)

-- | updateNodeWithPreorder takes the preorder decorated graphs (by block and character) and updates the
-- the preorder fields only using character info.  This leaves post and preorder assignment out of sync.
-- but that so can use incremental optimization on base decorated graph in other areas.
updateNodeWithPreorder :: V.Vector (V.Vector DecoratedGraph) -> V.Vector (V.Vector CharInfo) -> LG.LNode VertexInfo -> LG.LNode VertexInfo
updateNodeWithPreorder preOrderBlockTreeVV inCharInfoVV postOrderNode =
    let nodeLabel = snd postOrderNode
        nodeVertData = vertData nodeLabel
        newNodeVertData = V.zipWith3 (updateVertexBlock (fst postOrderNode)) preOrderBlockTreeVV nodeVertData inCharInfoVV
    in
    (fst postOrderNode, nodeLabel {vertData = newNodeVertData})

-- | updateVertexBlock takes a block of vertex data and updates preorder states of charactes via fmap
updateVertexBlock :: Int -> V.Vector DecoratedGraph -> V.Vector CharacterData -> V.Vector CharInfo -> V.Vector CharacterData
updateVertexBlock nodeIndex = V.zipWith3 (updatePreorderCharacter nodeIndex)

-- | updatePreorderCharacter updates the pre-order fields of character data for a vertex from a traversal
-- since there is single character optimized for each character decorated graph-- it is always teh 0th 0th character
-- exact are vectors so take care of multiple there.
-- need to care for issues of missing data
updatePreorderCharacter :: Int -> DecoratedGraph -> CharacterData -> CharInfo -> CharacterData
updatePreorderCharacter nodeIndex preOrderTree postOrderCharacter charInfo =
    --trace ("N:" ++ (show nodeIndex) ++ " B:" ++ (show blockIndex) ++ " C:" ++ (show characterIndex) ++ "\n" ++ (show $ vertData $ fromJust $ LG.lab preOrderTree nodeIndex)) (
    let maybePreOrderNodeLabel = LG.lab preOrderTree nodeIndex
        preOrderVertData = vertData $ fromJust maybePreOrderNodeLabel
        preOrderCharacterData
          | V.null preOrderVertData = emptyCharacter
          | V.null $ V.head preOrderVertData = emptyCharacter
          | otherwise = V.head $ V.head preOrderVertData -- (preOrderVertData V.! 0) V.! 0

    in
    -- this can heppen in naked parent node of prunned subGraph in branch swapping
    if isNothing maybePreOrderNodeLabel then emptyCharacter
            -- error ("Nothing node label in updatePreorderCharacter node: " ++ show nodeIndex)
    else
        updateCharacter postOrderCharacter preOrderCharacterData (charType charInfo)
    --)

-- | updateCharacter takes a postorder character and updates the preorder (final) fields with preorder data and character type
-- only updating preorder assignment--except for root, that is needed to draw final state for branch lengths
updateCharacter :: CharacterData -> CharacterData -> CharType  -> CharacterData
updateCharacter postOrderCharacter preOrderCharacter localCharType
  | localCharType == Add =
    postOrderCharacter { rangeFinal = rangeFinal preOrderCharacter
                       , rangeUnion = rangeUnion preOrderCharacter }

  | localCharType == NonAdd =
    postOrderCharacter { stateBVFinal = stateBVFinal preOrderCharacter
                       , stateBVUnion = stateBVUnion preOrderCharacter}

  | localCharType `elem` packedNonAddTypes =
    postOrderCharacter { packedNonAddFinal = packedNonAddFinal preOrderCharacter
                       , packedNonAddUnion = packedNonAddUnion preOrderCharacter }

  | localCharType == Matrix =
    postOrderCharacter { matrixStatesFinal = matrixStatesFinal preOrderCharacter
                       , matrixStatesUnion = matrixStatesUnion preOrderCharacter }

  | localCharType == AlignedSlim =
    postOrderCharacter { alignedSlimPrelim = alignedSlimPrelim preOrderCharacter
                       , alignedSlimFinal  = alignedSlimFinal preOrderCharacter
                       , alignedSlimUnion  = alignedSlimUnion preOrderCharacter
                      }

  | localCharType == AlignedWide =
    postOrderCharacter { alignedWidePrelim = alignedWidePrelim preOrderCharacter
                       , alignedWideFinal  = alignedWideFinal preOrderCharacter
                       , alignedWideUnion  = alignedWideUnion preOrderCharacter
                       }

  | localCharType == AlignedHuge =
    postOrderCharacter { alignedHugePrelim = alignedHugePrelim preOrderCharacter
                       , alignedHugeFinal  = alignedHugeFinal preOrderCharacter
                       , alignedHugeUnion  = alignedHugeUnion preOrderCharacter
                       }

  | localCharType == SlimSeq || localCharType == NucSeq =
    postOrderCharacter { slimAlignment = slimAlignment preOrderCharacter
                      , slimFinal = slimFinal preOrderCharacter
                      , slimIAFinal = slimIAFinal preOrderCharacter
                      , slimIAUnion = slimIAUnion preOrderCharacter
                      }

  | localCharType == WideSeq || localCharType == AminoSeq =
    postOrderCharacter { wideAlignment = wideAlignment preOrderCharacter
                      , wideFinal = wideFinal preOrderCharacter
                      , wideIAFinal = wideIAFinal preOrderCharacter
                      , wideIAUnion = wideIAUnion preOrderCharacter
                      }

  | localCharType == HugeSeq =
    postOrderCharacter { hugeAlignment = hugeAlignment preOrderCharacter
                      , hugeFinal = hugeFinal preOrderCharacter
                      , hugeIAFinal = hugeIAFinal preOrderCharacter
                      , hugeIAUnion = hugeIAUnion preOrderCharacter
                      }

  | otherwise = error ("Character type unimplemented : " ++ show localCharType)

-- | updateEdgeInfoSoftWired gets edge weights via block trees as opposed to canonical graph
-- this because not all edges present in all block/display trees
updateEdgeInfoSoftWired :: AssignmentMethod
                        -> V.Vector (V.Vector CharInfo)
                        -> V.Vector (V.Vector (V.Vector (LG.LNode VertexInfo), V.Vector (LG.LEdge EdgeInfo)))
                        -> Int
                        -> LG.LEdge EdgeInfo
                        -> LG.LEdge EdgeInfo
updateEdgeInfoSoftWired finalMethod inCharInfoVV blockTreePairVV rootIndex (uNode, vNode, edgeLabel) =
    if V.null blockTreePairVV then error "Empty node-edge pair vector in updateEdgeInfoSoftWired"
    else
        let (minWList, maxWList) = V.unzip $ V.zipWith (getEdgeBlockWeightSoftWired finalMethod uNode vNode rootIndex) inCharInfoVV blockTreePairVV
            localEdgeType = edgeType edgeLabel
            newEdgeLabel = EdgeInfo { minLength = V.sum minWList
                                    , maxLength = V.sum maxWList
                                    , midRangeLength = (V.sum minWList + V.sum maxWList) / 2.0
                                    , edgeType  = localEdgeType
                                    }
        in
        (uNode, vNode, newEdgeLabel)


-- | getEdgeBlockWeightSoftWired takes a block of character trees and maps character distances if edge exists in block tree
getEdgeBlockWeightSoftWired :: AssignmentMethod
                            -> Int
                            -> Int
                            -> Int
                            -> V.Vector CharInfo
                            -> V.Vector (V.Vector (LG.LNode VertexInfo), V.Vector (LG.LEdge EdgeInfo))
                            -> (VertexCost, VertexCost)
getEdgeBlockWeightSoftWired finalMethod uNode vNode rootIndex inCharInfoV blockTreePairV =
    let (minWList, maxWList) = V.unzip $ V.zipWith (getEdgeCharacterWeightSoftWired finalMethod uNode vNode rootIndex) inCharInfoV blockTreePairV
    in (V.sum minWList, V.sum maxWList)

-- | getEdgeCharacterWeightSoftWired gets the edge weight for an individual character
-- matches edge in either direction
-- need examine edge root as two edges from rootIndex
getEdgeCharacterWeightSoftWired :: AssignmentMethod
                                -> Int
                                -> Int
                                -> Int
                                -> CharInfo
                                -> (V.Vector (LG.LNode VertexInfo), V.Vector (LG.LEdge EdgeInfo))
                                -> (VertexCost, VertexCost)
getEdgeCharacterWeightSoftWired finalMethod uNode vNode rootIndex inCharInfo (nodeVect, edgeVect) =
    let foundVertexPair = getEdgeVerts uNode vNode rootIndex nodeVect edgeVect
        (uLabel, vLabel) = fromJust foundVertexPair
        uCharacter = V.head $ V.head $ vertData uLabel
        vCharacter = V.head $ V.head $ vertData vLabel
    in
    -- if edge not present and not around root then return  no costs
    if foundVertexPair == Nothing then (0,0)
    else getCharacterDistFinal finalMethod uCharacter vCharacter inCharInfo


-- | getEdgeVerts returns vertex labels if edge in vect or if a virtual edge including root
getEdgeVerts :: Int -> Int -> Int -> V.Vector (LG.LNode VertexInfo) -> V.Vector (LG.LEdge EdgeInfo) -> Maybe (VertexInfo, VertexInfo)
getEdgeVerts uNode vNode rootIndex nodeVect edgeVect =
    -- trace ("GEV:" ++ (show (uNode, vNode, rootIndex) ++ " nodes " ++ (show $ fmap fst nodeVect) ++ " edges " ++ (show $ fmap LG.toEdge edgeVect))) (

    --hack or display edge check I'm not sure--not all edges are in all display trees
    if (uNode >= V.length nodeVect) || (vNode >= V.length nodeVect) then Nothing

    else if edgeInVect (uNode, vNode) edgeVect then Just (snd $ nodeVect V.! uNode, snd $ nodeVect V.! vNode)
    else if (edgeInVect (rootIndex, uNode) edgeVect) && (edgeInVect (rootIndex, vNode) edgeVect) then Just (snd $ nodeVect V.! uNode, snd $ nodeVect V.! vNode)
    else Nothing
    -- )

-- | edgeInVect takes an edges and returns True if in Vector, False otherwise
edgeInVect :: (Int , Int) -> V.Vector (LG.LEdge EdgeInfo) -> Bool
edgeInVect (u, v) edgeVect =
    if V.null edgeVect then False
    else
        let (a, b, _) = V.head edgeVect
        in
        if (u, v) == (a, b) then True
        else if (v, u) == (a, b) then True
        else edgeInVect (u, v) (V.tail edgeVect)

-- | updateEdgeInfoTree takes a Decorated graph--fully labelled post and preorder and and edge and
-- gets edge info--basically lengths
-- this for a tree in that all edges are present in all character/block trees
updateEdgeInfoTree :: AssignmentMethod -> V.Vector (V.Vector CharInfo) -> V.Vector (LG.LNode VertexInfo) -> LG.LEdge EdgeInfo -> LG.LEdge EdgeInfo
updateEdgeInfoTree finalMethod inCharInfoVV nodeVector (uNode, vNode, edgeLabel) =
    if V.null nodeVector then error "Empty node list in updateEdgeInfo"
    else
        let (minW, maxW) = getEdgeWeight finalMethod inCharInfoVV nodeVector (uNode, vNode)
            midW = (minW + maxW) / 2.0
            localEdgeType = edgeType edgeLabel
            newEdgeLabel = EdgeInfo { minLength = minW
                                    , maxLength = maxW
                                    , midRangeLength = midW
                                    , edgeType  = localEdgeType
                                    }
        in
        (uNode, vNode, newEdgeLabel)


-- | updateEdgeInfoTreeMap takes a Decorated graph--fully labelled post and preorder and and edge and
-- gets edge info--basically lengths
-- this for a tree in that all edges are present in all character/block trees
-- uses MAP as opposed to index
updateEdgeInfoTreeMap :: AssignmentMethod -> V.Vector (V.Vector CharInfo) -> MAP.Map Int (LG.LNode VertexInfo) -> LG.LEdge EdgeInfo -> LG.LEdge EdgeInfo
updateEdgeInfoTreeMap finalMethod inCharInfoVV nodeMap (uNode, vNode, edgeLabel) =
    if MAP.null nodeMap then error "Empty node MAP in updateEdgeInfo"
    else
        let (minW, maxW) = getEdgeWeightMap finalMethod inCharInfoVV nodeMap (uNode, vNode)
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
-- the indexing depends on the graph having all leaves in the graph which may not happen
-- during graph swapping
getEdgeWeight :: AssignmentMethod -> V.Vector (V.Vector CharInfo) -> V.Vector (LG.LNode VertexInfo) -> (Int, Int) -> (VertexCost, VertexCost)
getEdgeWeight finalMethod inCharInfoVV nodeVector (uNode, vNode) =
    if V.null nodeVector then error "Empty node list in getEdgeWeight"
    else
        -- trace ("GEW: " ++ (show $ fmap fst nodeVector) ++ " " ++ (show (uNode, vNode))) (
        let uNodeInfo = vertData $ snd $ nodeVector V.! uNode
            vNodeInfo = vertData $ snd $ nodeVector V.! vNode
            blockCostPairs = V.zipWith3 (getBlockCostPairsFinal finalMethod) uNodeInfo vNodeInfo inCharInfoVV
            minCost = sum $ fmap fst blockCostPairs
            maxCost = sum $ fmap snd blockCostPairs
        in

        (minCost, maxCost)
        -- )

-- | getEdgeWeightMap takes a preorder decorated decorated graph and an edge and gets the weight information for that edge
-- basically a min/max distance between the two
-- in this case based on map or vertices rather than direct indexing.
-- the indexing depends on the graph having all leaves in the graph which may not happen
-- during graph swapping
getEdgeWeightMap :: AssignmentMethod -> V.Vector (V.Vector CharInfo) -> MAP.Map Int (LG.LNode VertexInfo) -> (Int, Int) -> (VertexCost, VertexCost)
getEdgeWeightMap finalMethod inCharInfoVV nodeMap (uNode, vNode) =
    if MAP.null nodeMap then error "Empty node map in getEdgeWeight"
    else
        -- trace ("GEWM: " ++ (show $ MAP.toList nodeMap) ++ " " ++ (show (uNode, vNode))) (
        let uNodeInfo = vertData $ snd $ nodeMap MAP.! uNode
            vNodeInfo = vertData $ snd $ nodeMap MAP.! vNode
            blockCostPairs = V.zipWith3 (getBlockCostPairsFinal finalMethod) uNodeInfo vNodeInfo inCharInfoVV
            minCost = sum $ fmap fst blockCostPairs
            maxCost = sum $ fmap snd blockCostPairs
        in

        (minCost, maxCost)
        -- )

-- | getBlockCostPairsFinal takes a block of two nodes and character infomation and returns the min and max block branch costs
getBlockCostPairsFinal :: AssignmentMethod -> V.Vector CharacterData -> V.Vector CharacterData -> V.Vector CharInfo -> (VertexCost, VertexCost)
getBlockCostPairsFinal finalMethod uNodeCharDataV vNodeCharDataV charInfoV =
    let characterCostPairs = V.zipWith3 (getCharacterDistFinal finalMethod) uNodeCharDataV vNodeCharDataV charInfoV
        minCost = sum $ fmap fst characterCostPairs
        maxCost = sum $ fmap snd characterCostPairs
    in
    (minCost, maxCost)

-- | getCharacterDistFinal takes a pair of characters and character type, returning the minimum and maximum character distances
-- for sequence characters this is based on slim/wide/hugeAlignment field, hence all should be O(n) in num characters/sequence length
getCharacterDistFinal :: AssignmentMethod -> CharacterData -> CharacterData -> CharInfo -> (VertexCost, VertexCost)
getCharacterDistFinal finalMethod uCharacter vCharacter charInfo =
    let thisWeight = weight charInfo
        thisMatrix = costMatrix charInfo
        thisCharType = charType charInfo
        lChangeCost = changeCost charInfo
        lNoChangeCost = noChangeCost charInfo
    in
    -- no nded to do nochange/change--all recoded in that case
    if thisCharType == Add then
        let --minCost = localCost (M.intervalAdd thisWeight uCharacter vCharacter)
            (minDiffV, maxDiffV) = V.unzip $ V.zipWith maxMinIntervalDiff  (rangeFinal uCharacter) (rangeFinal vCharacter)

            minCost = thisWeight * (fromIntegral $ V.sum minDiffV)
            maxCost = thisWeight * (fromIntegral $ V.sum maxDiffV)
        in
        (minCost, maxCost)

    -- assumes noChangeCost < changeCost for PMDL/ML
    else if thisCharType == NonAdd then
        let -- minCost = localCost (M.interUnion thisWeight uCharacter vCharacter)
            minDiff = length $ V.filter (==False) $ V.zipWith hasBVIntersection (stateBVFinal uCharacter) (stateBVFinal vCharacter)
            maxDiff = length $ V.filter (==False) $ V.zipWith equalAndSingleState (stateBVFinal uCharacter) (stateBVFinal vCharacter)
            maxCost = thisWeight * fromIntegral maxDiff
            minCost = thisWeight * fromIntegral minDiff
            minNoChange = (length (stateBVFinal uCharacter)) - minDiff
            maxNoChange = (length (stateBVFinal uCharacter)) - maxDiff
            minCost' = thisWeight * ((lNoChangeCost * (fromIntegral minNoChange)) + (lChangeCost * (fromIntegral minDiff)))
            maxCost' = thisWeight * ((lNoChangeCost * (fromIntegral maxNoChange)) + (lChangeCost * (fromIntegral maxDiff)))
            
        in
        if lNoChangeCost == 0.0 then (minCost, maxCost)
        else (minCost', maxCost')

    else if thisCharType `elem` packedNonAddTypes then
        let -- minCost = localCost (BP.median2Packed thisCharType uCharacter vCharacter)
            (minDiffV, maxDiffV) = UV.unzip $ UV.zipWith (BP.minMaxCharDiff thisCharType (lNoChangeCost, lChangeCost)) (packedNonAddFinal uCharacter) (packedNonAddFinal vCharacter)
            maxCost = thisWeight * (UV.sum maxDiffV)
            minCost = thisWeight * (UV.sum minDiffV)
        in
        (minCost, maxCost)

    else if thisCharType == Matrix then
        let minMaxListList= V.zipWith (minMaxMatrixDiff thisMatrix)  (fmap getLowestCostMatrixStates (matrixStatesFinal uCharacter)) (fmap getLowestCostMatrixStates (matrixStatesFinal vCharacter))
            minDiff = V.sum $ fmap fst minMaxListList
            maxDiff = V.sum $ fmap snd minMaxListList
            minCost = thisWeight * fromIntegral minDiff
            maxCost = thisWeight * fromIntegral maxDiff
        in
        (minCost, maxCost)

    else if thisCharType `elem` prealignedCharacterTypes then
        let
            (minDiff, maxDiff) = unzip $ zipWith (M.generalSequenceDiff thisMatrix (length thisMatrix)) (GV.toList $ alignedSlimFinal uCharacter) (GV.toList $ alignedSlimFinal vCharacter)
            minCost = thisWeight * fromIntegral (sum minDiff)
            maxCost = thisWeight * fromIntegral (sum maxDiff)
        in
        (minCost, maxCost)

    else if thisCharType == SlimSeq || thisCharType == NucSeq then
        let minMaxDiffList = if finalMethod == DirectOptimization then
                                let uFinal = M.makeDynamicCharacterFromSingleVector (slimFinal uCharacter)
                                    vFinal = M.makeDynamicCharacterFromSingleVector (slimFinal vCharacter)
                                    newEdgeCharacter = M.getDOMedianCharInfo charInfo (uCharacter {slimGapped = uFinal}) (vCharacter {slimGapped = vFinal})
                                    (newU, _, newV) = slimGapped newEdgeCharacter
                                in
                                --trace ("GCD:\n" ++ (show (slimFinal uCharacter, newU)) ++ "\n" ++ (show (slimFinal vCharacter, newV)) ++ "\nDO Cost:" ++ (show doCOST))
                                zipWith  (M.generalSequenceDiff thisMatrix (length thisMatrix)) (GV.toList newU) (GV.toList newV)
                             else zipWith  (M.generalSequenceDiff thisMatrix (length thisMatrix)) (GV.toList $ slimIAFinal uCharacter) (GV.toList $ slimIAFinal vCharacter)
            (minDiff, maxDiff) = unzip minMaxDiffList
            minCost = thisWeight * fromIntegral (sum minDiff)
            maxCost = thisWeight * fromIntegral (sum maxDiff)
        in
        --trace ("MMDL: " ++ (show minCost) ++ " " ++ (show maxCost))
        (minCost, maxCost)


    else if thisCharType == WideSeq || thisCharType == AminoSeq then
        let minMaxDiffList = if finalMethod == DirectOptimization then
                                let uFinal = M.makeDynamicCharacterFromSingleVector (wideFinal uCharacter)
                                    vFinal = M.makeDynamicCharacterFromSingleVector (wideFinal vCharacter)
                                    newEdgeCharacter = M.getDOMedianCharInfo charInfo (uCharacter {wideGapped = uFinal}) (vCharacter {wideGapped = vFinal})
                                    (newU, _, newV) = wideGapped newEdgeCharacter
                                in
                                --trace ("GCD:\n" ++ (show m) ++ "\n" ++ (show (uFinal, newU)) ++ "\n" ++ (show (vFinal, newV)))
                                zipWith  (M.generalSequenceDiff thisMatrix (length thisMatrix)) (GV.toList newU) (GV.toList newV)
                             else GV.toList $ GV.zipWith (M.generalSequenceDiff thisMatrix (length thisMatrix))  (wideIAFinal uCharacter) (wideIAFinal vCharacter)
            (minDiff, maxDiff) = unzip minMaxDiffList
            minCost = thisWeight * fromIntegral (sum minDiff)
            maxCost = thisWeight * fromIntegral (sum maxDiff)
        in
        (minCost, maxCost)

    else if thisCharType == HugeSeq then
        let minMaxDiffList = if finalMethod == DirectOptimization then
                                let uFinal = M.makeDynamicCharacterFromSingleVector (hugeFinal uCharacter)
                                    vFinal = M.makeDynamicCharacterFromSingleVector (hugeFinal vCharacter)
                                    newEdgeCharacter = M.getDOMedianCharInfo charInfo (uCharacter {hugeGapped = uFinal}) (vCharacter {hugeGapped = vFinal})
                                    (newU, _, newV) = hugeGapped newEdgeCharacter
                                in
                                -- trace ("GCD:\n" ++ (show (uFinal, newU)) ++ "\n" ++ (show (vFinal, newV)))
                                zipWith  (M.generalSequenceDiff thisMatrix (length thisMatrix)) (GV.toList newU) (GV.toList newV)
                             else GV.toList $ GV.zipWith (M.generalSequenceDiff thisMatrix (length thisMatrix)) (hugeIAFinal uCharacter) (hugeIAFinal vCharacter)
            (minDiff, maxDiff) = unzip minMaxDiffList
            minCost = thisWeight * fromIntegral (sum minDiff)
            maxCost = thisWeight * fromIntegral (sum maxDiff)
        in
        (minCost, maxCost)

    else error ("Character type not recognized/unimplemented : " ++ show thisCharType)
    where hasBVIntersection a b = (not . BV.isZeroVector) (a .&. b) 
          equalAndSingleState a b = if (a == b) && (popCount a == 1) then True else False

-- | zero2Gap converts a '0' or no bits set to gap (indel) value
zero2Gap :: (FiniteBits a) => a -> a
zero2Gap inVal = if popCount inVal == 0 then bit gapIndex
                 else inVal

{-
-- | zero2GapWide converts a '0' or no bits set to gap (indel) value
zero2GapWide :: Word64 -> Word64 -> Word64
zero2GapWide gapChar inVal = if popCount inVal == 0  then bit gapIndex
                         else inVal

-- | zero2GapBV converts a '0' or no bits set to gap (indel) value
zero2GapBV :: BV.BitVector -> BV.BitVector -> BV.BitVector
zero2GapBV gapChar inVal = if popCount inVal == 0 then bit gapIndex
                         else inVal
-}

-- | maxIntervalDiff takes two ranges and gets the maximum difference between the two based on differences
-- in upp and lower ranges.
maxMinIntervalDiff :: (Int, Int)-> (Int, Int) -> (Int, Int)
maxMinIntervalDiff (a,b) (x,y) =
    let upper = max b y - min b y
        lower = max a x - min a x
    in
    (min upper lower, max  upper lower)

-- | getLowestCostMatrixStates takes a Vector Triple for matrix charxcter and returns lowest cost states as vector
-- of Ints
getLowestCostMatrixStates :: V.Vector MatrixTriple -> V.Vector Int
getLowestCostMatrixStates tripleVect =
    if V.null tripleVect then V.empty
    else
        let minCost = minimum $ fmap fst3 tripleVect
            stateCostPairList = V.zip (V.fromList [0..(V.length tripleVect - 1)]) (fmap fst3 tripleVect)
            (minStateVect, _) = V.unzip $ V.filter ((== minCost) . snd) stateCostPairList
        in
        minStateVect


-- | minMaxMatrixDiff takes twovetors of states and calculates the minimum and maximum state differnce cost
-- between the two
minMaxMatrixDiff :: S.Matrix Int -> V.Vector Int -> V.Vector Int -> (Int, Int)
minMaxMatrixDiff localCostMatrix uStatesV vStatesV =
    let statePairs = (V.toList uStatesV, V.toList vStatesV)
        cartesianPairs = cartProdPair statePairs
        costList = fmap (localCostMatrix S.!) cartesianPairs
    in
    {-THis ti check for errors
    if (not . null) costList then (minimum costList, maximum costList)
    else (-1, -1)
    -}
    -- trace ("MMD: " ++ (show (statePairs,cartesianPairs)))
    (minimum costList, maximum costList)
    

-- | createFinalAssignment takes vertex data (child or current vertex) and creates the final
-- assignment from parent (if not root or leaf) and 'child' ie current vertex
-- if root or leaf preliminary is assigned to final
   -- need to watch zipping for missing sequence data
-- this creates the IA during preorder from which final assignments are contructed
-- via addition post and preorder passes on IA fields.
createFinalAssignmentOverBlocks :: GlobalSettings
                                -> AssignmentMethod
                                -> Bool
                                -> NodeType
                                -> VertexBlockData
                                -> VertexBlockData
                                -> VertexBlockDataMaybe -- second parent if indegree 2 node
                                -> CharInfo
                                -> Bool
                                -> Bool
                                -> Bool
                                -> VertexBlockData
createFinalAssignmentOverBlocks inGS finalMethod staticIA childType childBlockData parentBlockData parent2BlockDataM charInfo isLeft isInOutDegree1 isIn2Out1 =
   -- if root or leaf final assignment <- preliminary asssignment
   V.zipWith3 (assignFinal inGS finalMethod staticIA childType isLeft charInfo isInOutDegree1 isIn2Out1) childBlockData parentBlockData parent2BlockDataM


-- | assignFinal takes a vertex type and single block of zip3 of child info, parent info, and character type
-- to create pre-order assignments
assignFinal :: GlobalSettings
            -> AssignmentMethod
            -> Bool
            -> NodeType
            -> Bool
            -> CharInfo
            -> Bool
            -> Bool
            -> V.Vector CharacterData
            -> V.Vector CharacterData
            -> V.Vector (Maybe CharacterData)
            -> V.Vector CharacterData
assignFinal inGS finalMethod staticIA childType isLeft charInfo isOutDegree1 isIn2Out1 = V.zipWith3 (setFinal inGS finalMethod staticIA childType isLeft charInfo isOutDegree1 isIn2Out1)

-- | setFinal takes a vertex type and single character of zip3 of child info, parent info, and character type
-- to create pre-order assignments
   -- | setFinalHTU takes a single character and its parent and sets the final state to prelim based
-- on character info.
-- non exact characters are vectors of characters of same type
-- this does the same things for sequence types, but also
-- performs preorder logic for exact characters
-- staticIA flage is for IA and static only optimization used in IA heuriastics for DO
-- no IA for networks--at least for now.Bool ->
setFinal :: GlobalSettings -> AssignmentMethod -> Bool -> NodeType -> Bool -> CharInfo -> Bool -> Bool -> CharacterData -> CharacterData -> Maybe CharacterData -> CharacterData
setFinal inGS finalMethod staticIA childType isLeft charInfo isIn1Out1 isIn2Out1 childChar parentChar parent2CharM =
   let localCharType = charType charInfo
       symbolCount = toEnum $ length $ costMatrix charInfo :: Int
       isTree = (graphType inGS) == Tree
   in
   -- Three cases, Root, leaf, HTU
   -- trace ("set final:" ++ (show (finalMethod, staticIA)) ++ " " ++ (show childType) ++ " " ++ (show isLeft) ++ " " ++ (show isIn1Out1) ++ " " ++ (show isIn2Out1)) (
   if childType == RootNode then

      if localCharType == Add then
         childChar {rangeFinal = snd3 $ rangePrelim childChar}

      else if localCharType == NonAdd then
        childChar {stateBVFinal = snd3 $ stateBVPrelim childChar}

      else if localCharType `elem` packedNonAddTypes then
        childChar {packedNonAddFinal = snd3 $ packedNonAddPrelim childChar}

      else if localCharType == Matrix then
        childChar {matrixStatesFinal = setMinCostStatesMatrix (fromEnum symbolCount) (matrixStatesPrelim childChar)}

      else if localCharType == AlignedSlim then
        childChar {alignedSlimFinal = snd3 $ alignedSlimPrelim childChar}

      else if localCharType == AlignedWide then
        childChar {alignedWideFinal = snd3 $ alignedWidePrelim childChar}

      else if localCharType == AlignedHuge then
        childChar {alignedHugeFinal = snd3 $ alignedHugePrelim childChar}

      -- need to set both final and alignment for sequence characters
      else if (localCharType == SlimSeq) || (localCharType == NucSeq) then
         let finalAssignment' = extractMedians $ slimGapped childChar
         in
         -- trace ("TNFinal-Root: " ++ (show finalAssignment') ++ " " ++ (show (GV.length finalAssignment', slimGapped childChar))) $
         --traceNoLF ("TNFinal-Root") $
         if staticIA then childChar {slimIAFinal = extractMediansGapped $ slimIAPrelim childChar}
         else childChar { slimFinal = finalAssignment'
                        , slimAlignment = slimGapped childChar
                                          --if isTree then slimGapped childChar
                                          --else mempty
                        }
         

      else if (localCharType == WideSeq) || (localCharType == AminoSeq) then
         let finalAssignment' = extractMedians $ wideGapped childChar
         in
         if staticIA then childChar {wideIAFinal = extractMediansGapped $ wideIAPrelim childChar}
         else childChar { wideFinal = finalAssignment'
                        , wideAlignment = wideGapped childChar
                                          -- if isTree then wideGapped childChar
                                          -- else mempty
                        }

      else if localCharType == HugeSeq then
         let finalAssignment' = extractMedians $ hugeGapped childChar
         in
         if staticIA then childChar {hugeIAFinal = extractMediansGapped $ hugeIAPrelim childChar}
         else childChar { hugeFinal = finalAssignment'
                        , hugeAlignment = hugeGapped childChar
                                          -- if isTree then hugeGapped childChar
                                          -- else mempty
                        }

      else error ("Unrecognized/implemented character type: " ++ show localCharType)

   else if childType == LeafNode then
      -- since leaf no neeed to precess final alignment fields for sequence characters
      if localCharType == Add then
        childChar {rangeFinal = snd3 $ rangePrelim childChar}

      else if localCharType == NonAdd then
        childChar {stateBVFinal = snd3 $ stateBVPrelim childChar}

      else if localCharType `elem` packedNonAddTypes then
        childChar {packedNonAddFinal = snd3 $ packedNonAddPrelim childChar}

      else if localCharType == Matrix then
         childChar {matrixStatesFinal = setMinCostStatesMatrix (fromEnum symbolCount) (matrixStatesPrelim childChar)}

      else if localCharType == AlignedSlim then
         childChar {alignedSlimFinal = extractMediansGapped $ alignedSlimPrelim childChar}

      else if localCharType == AlignedWide then
         childChar {alignedWideFinal = extractMediansGapped $ alignedWidePrelim childChar}

      else if localCharType == AlignedHuge then
        childChar {alignedHugeFinal = extractMediansGapped $ alignedHugePrelim childChar}

      -- need to set both final and alignment for sequence characters
      else if (localCharType == SlimSeq) || (localCharType == NucSeq) then
         let finalAlignment = doPreOrderWithParentCheck isLeft (slimAlignment parentChar) (slimGapped parentChar) (slimGapped childChar)
             -- finalAssignment' = extractMedians finalAlignment
         in
         -- traceNoLF ("TNFinal-Leaf") $
         -- trace ("TNFinal-Leaf:" ++ (show (GV.length $ fst3  (slimAlignment parentChar), GV.length $ fst3 finalAlignment, isLeft, (slimAlignment parentChar), (slimGapped parentChar) ,(slimGapped childChar))) ++ "\n->" ++ (show finalAlignment)) $
         if staticIA then childChar {slimIAFinal = extractMediansGapped $ slimIAPrelim childChar}
         else childChar { slimFinal = extractMedians $ slimGapped childChar -- finalAssignment'
                        , slimAlignment = finalAlignment
                                          --if isTree then finalAlignment
                                          --else mempty
                        , slimIAPrelim  = finalAlignment
                                          -- if isTree then finalAlignment
                                          -- else mempty
                        , slimIAFinal  =  extractMediansGapped finalAlignment
                                          --if isTree then extractMediansGapped $ finalAlignment
                                          --else mempty
                        }
         --)

      else if (localCharType == WideSeq) || (localCharType == AminoSeq) then
         let finalAlignment = doPreOrderWithParentCheck isLeft (wideAlignment parentChar) (wideGapped parentChar) (wideGapped childChar)
                              
             -- finalAssignment' = extractMedians finalAlignment
         in
         if staticIA then childChar {wideIAFinal = extractMediansGapped $ wideIAPrelim childChar}
         else childChar { wideFinal = extractMedians $ wideGapped childChar -- finalAssignment'
                        , wideAlignment = finalAlignment
                                          --if isTree then finalAlignment
                                          --else mempty
                        , wideIAPrelim =  finalAlignment
                                          --if isTree then finalAlignment
                                          --else mempty
                        , wideIAFinal =   extractMediansGapped finalAlignment
                                          --if isTree then extractMediansGapped $ finalAlignment
                                          --else mempty
                        }

      else if localCharType == HugeSeq then
         let finalAlignment = doPreOrderWithParentCheck isLeft (hugeAlignment parentChar) (hugeGapped parentChar) (hugeGapped childChar)

             -- finalAssignment' = extractMedians finalAlignment
         in
         if staticIA then childChar {hugeIAFinal = extractMediansGapped $ hugeIAPrelim childChar}
         else childChar { hugeFinal = extractMedians $ hugeGapped childChar -- finalAssignment'
                        , hugeAlignment = finalAlignment
                                          -- if isTree then finalAlignment
                                          -- else mempty
                        , hugeIAPrelim =  finalAlignment
                                          -- if isTree then finalAlignment
                                          -- else mempty
                        , hugeIAFinal =   extractMediansGapped finalAlignment
                                          -- if isTree then extractMediansGapped $ finalAlignment
                                          -- else mempty
                        }

      else error ("Unrecognized/implemented character type: " ++ show localCharType)

   else if childType == TreeNode && not isIn1Out1 then

      if localCharType == Add then
         -- add logic for pre-order
         let finalAssignment' = additivePreorder (rangePrelim childChar) (rangeFinal parentChar)
         in
         childChar {rangeFinal = finalAssignment'}

      else if localCharType == NonAdd then
         -- add logic for pre-order
         let finalAssignment' = nonAdditivePreorder (stateBVPrelim childChar) (stateBVFinal parentChar)
         in
         childChar {stateBVFinal = finalAssignment'}

      else if localCharType `elem` packedNonAddTypes then
         let finalAssignment' = BP.packedPreorder localCharType (packedNonAddPrelim childChar) (packedNonAddFinal parentChar)
         in
         childChar {packedNonAddFinal = finalAssignment'}

      else if localCharType == Matrix then
         -- add logic for pre-order
         let finalAssignment' = matrixPreorder isLeft (matrixStatesPrelim childChar) (matrixStatesFinal parentChar)
         in
         childChar {matrixStatesFinal = finalAssignment'}

      else if localCharType == AlignedSlim then
        let alignedFinal = M.getFinal3WaySlim (slimTCM charInfo) (alignedSlimFinal parentChar) (fst3 $ alignedSlimPrelim childChar) (thd3 $ alignedSlimPrelim childChar)
        in
        childChar {alignedSlimFinal = alignedFinal}

      else if localCharType == AlignedWide then
       let alignedFinal = M.getFinal3WayWideHuge (wideTCM charInfo) (alignedWideFinal parentChar) (fst3 $ alignedWidePrelim childChar) (thd3 $ alignedWidePrelim childChar)
        in
        childChar {alignedWideFinal = alignedFinal}

      else if localCharType == AlignedHuge then
        let alignedFinal = M.getFinal3WayWideHuge (hugeTCM charInfo) (alignedHugeFinal parentChar) (fst3 $ alignedHugePrelim childChar) (thd3 $ alignedHugePrelim childChar)
        in
        childChar {alignedHugeFinal = alignedFinal}

      -- need to set both final and alignment for sequence characters
      else if (localCharType == SlimSeq) || (localCharType == NucSeq) then
         -- traceNoLF ("TNFinal-Tree") $
         let finalGapped = doPreOrderWithParentCheck isLeft (slimAlignment parentChar) (slimGapped parentChar) (slimGapped childChar)
                              
             finalAssignmentDO = if finalMethod == DirectOptimization then
                                    let parentFinalDC = M.makeDynamicCharacterFromSingleVector (slimFinal parentChar)
                                        parentFinal = (parentFinalDC, mempty, mempty)
                                        -- parentGapped = (slimGapped parentChar, mempty, mempty)
                                        childGapped = (slimGapped childChar, mempty, mempty)
                                        finalAssignmentDOGapped = fst3 $ getDOFinal charInfo parentFinal childGapped
                                    in
                                    extractMedians finalAssignmentDOGapped
                                    -- really could/should be mempty since overwritten by IA later
                                 else extractMedians finalGapped
         in
         --trace ("TNFinal-Tree:" ++ (show (SV.length $ fst3  (slimAlignment parentChar), SV.length $ fst3 finalGapped,isLeft, (slimAlignment parentChar), (slimGapped parentChar) ,(slimGapped childChar))) ++ "->" ++ (show finalGapped)) (
         if staticIA then M.makeIAFinalCharacter finalMethod charInfo childChar parentChar
         else childChar { slimFinal = finalAssignmentDO
                        , slimAlignment = finalGapped
                                         --if isTree then finalGapped
                                         -- else mempty
                        }
         --)

      else if (localCharType == WideSeq) || (localCharType == AminoSeq) then
         let finalGapped = doPreOrderWithParentCheck isLeft (wideAlignment parentChar) (wideGapped parentChar) (wideGapped childChar)

             finalAssignmentDO = if finalMethod == DirectOptimization then
                                    let parentFinalDC = M.makeDynamicCharacterFromSingleVector (wideFinal parentChar)
                                        parentFinal = (mempty, parentFinalDC, mempty)
                                        --parentGapped = (mempty, wideGapped parentChar, mempty)
                                        childGapped = (mempty, wideGapped childChar, mempty)
                                        finalAssignmentDOGapped = snd3 $ getDOFinal charInfo parentFinal childGapped
                                    in
                                    extractMedians finalAssignmentDOGapped
                                 else extractMedians finalGapped
         in
         if staticIA then M.makeIAFinalCharacter finalMethod charInfo childChar parentChar
         else childChar { wideFinal = finalAssignmentDO
                        , wideAlignment = finalGapped
                                         --if isTree then finalGapped
                                          --else mempty
                        }

      else if localCharType == HugeSeq then
         let finalGapped = doPreOrderWithParentCheck isLeft (hugeAlignment parentChar) (hugeGapped parentChar) (hugeGapped childChar)

             finalAssignmentDO = if finalMethod == DirectOptimization then
                                    let parentFinalDC = M.makeDynamicCharacterFromSingleVector (hugeFinal parentChar)
                                        parentFinal = (mempty, mempty, parentFinalDC)
                                        -- parentGapped = (mempty, mempty, hugeGapped parentChar)
                                        childGapped = (mempty, mempty, hugeGapped childChar)
                                        finalAssignmentDOGapped = thd3 $ getDOFinal charInfo parentFinal childGapped
                                    in
                                     extractMedians finalAssignmentDOGapped
                                 else extractMedians finalGapped
         in
         if staticIA then M.makeIAFinalCharacter finalMethod charInfo childChar parentChar
         else childChar { hugeFinal = finalAssignmentDO
                        , hugeAlignment = finalGapped
                                         --if isTree then finalGapped
                                         -- else mempty
                        }

      else error ("Unrecognized/implemented character type: " ++ show localCharType)

   -- display tree indegree=outdegree=1
   -- since display trees here--indegree should be one as well
   -- this doens't work--need to redo pass logic--perhaps by doing "grandparent"
   else if isIn1Out1 then
      -- trace ("InOut1 preorder") (
      if localCharType == Add then
         childChar {rangeFinal = rangeFinal parentChar}

      else if localCharType == NonAdd then
         childChar {stateBVFinal = stateBVFinal parentChar}

      else if localCharType `elem` packedNonAddTypes then
         childChar {packedNonAddFinal = packedNonAddFinal parentChar}

      else if localCharType == Matrix then
         childChar {matrixStatesFinal = matrixStatesFinal parentChar}

      else if localCharType == AlignedSlim then
        childChar {alignedSlimFinal = alignedSlimFinal parentChar}

      else if localCharType == AlignedWide then
        childChar {alignedWideFinal = alignedWideFinal parentChar}

      else if localCharType == AlignedHuge then
        childChar {alignedHugeFinal = alignedHugeFinal parentChar}

      -- need to set both final and alignment for sequence characters
      else if (localCharType == SlimSeq) || (localCharType == NucSeq) then -- parentChar
         --traceNoLF ("TNFinal-1/1") $
         -- trace ("TNFinal-1/1:" ++ (show (isLeft, (slimAlignment parentChar), (slimGapped parentChar) ,(slimGapped childChar)))) $
         if staticIA then childChar { slimIAFinal = slimIAFinal parentChar}
         else childChar { slimFinal = slimFinal parentChar
                        , slimAlignment = slimAlignment parentChar 
                                          -- if isTree then slimAlignment parentChar -- finalGappedO -- slimAlignment parentChar -- finalGappedO-- slimAlignment parentChar
                                          -- else mempty
                        , slimGapped = slimGapped parentChar -- slimGapped' -- slimGapped parentChar -- finalGappedO --slimGapped parentChar
                        -- , slimIAPrelim = slimIAPrelim parentChar
                        , slimIAFinal =  slimFinal parentChar
                                        -- if isTree then slimFinal parentChar
                                        --else mempty
                        }
        -- )

    else if (localCharType == WideSeq) || (localCharType == AminoSeq) then -- parentChar
         -- trace ("TNFinal-1/1:" ++ (show (isLeft, (slimAlignment parentChar), (slimGapped parentChar) ,(slimGapped childChar)))) (
         if staticIA then childChar { wideIAFinal = wideIAFinal parentChar}
         else childChar { wideFinal = wideFinal parentChar
                        , wideAlignment = if isTree then wideAlignment parentChar -- finalGappedO -- wideAlignment parentChar -- finalGappedO-- wideAlignment parentChar
                                          else mempty
                        , wideGapped = wideGapped parentChar -- wideGapped' -- wideGapped parentChar -- finalGappedO --wideGapped parentChar
                        -- , wideIAPrelim = wideIAPrelim parentChar
                        , wideIAFinal = wideFinal parentChar
                                        -- if isTree then wideFinal parentChar
                                        -- else mempty
                        }
        -- )

    else if (localCharType == HugeSeq)  then -- parentChar
         -- trace ("TNFinal-1/1:" ++ (show (isLeft, (hugeAlignment parentChar), (hugeGapped parentChar) ,(hugeGapped childChar)))) (
         if staticIA then childChar { hugeIAFinal = hugeIAFinal parentChar}
         else childChar { hugeFinal = hugeFinal parentChar
                        , hugeAlignment = if isTree then hugeAlignment parentChar -- finalGappedO -- hugeAlignment parentChar -- finalGappedO-- hugeAlignment parentChar
                                          else mempty
                        , hugeGapped = hugeGapped parentChar -- hugeGapped' -- hugeGapped parentChar -- finalGappedO --hugeGapped parentChar
                        -- , hugeIAPrelim = hugeIAPrelim parentChar
                        , hugeIAFinal = hugeFinal parentChar
                                        --if isTree then hugeFinal parentChar
                                        --else mempty
                        }
        -- )

      else error ("Unrecognized/implemented character type: " ++ show localCharType)
      -- )

      --for Hardwired graphs
   else if isIn2Out1 then
        if (parent2CharM == Nothing) then error ("Nothing parent2char in setFinal")
        else 
            -- trace ("SF: " ++ "makeing 3-way final")
            TW.threeMedianFinal charInfo parentChar (fromJust parent2CharM) childChar

   else error ("Node type should not be here (pre-order on tree node only): " ++ show  childType)
   -- )

-- | doPreOrderWithParentCheck performs post order losig if parent non-zero--otherwise returns preliminary assignment
doPreOrderWithParentCheck :: (FiniteBits e, GV.Vector v e) => Bool -> (v e, v e, v e) -> (v e, v e, v e) -> (v e, v e, v e) -> (v e, v e, v e) 
doPreOrderWithParentCheck isLeft alignmentParent gappedParent gappedChild =
    if not $ GV.null $ extractMediansGapped alignmentParent then
            DOP.preOrderLogic isLeft alignmentParent gappedParent gappedChild 
    else gappedChild

-- | getDOFinal takes parent final, and node gapped (including its parent gapped) and performs a DO median
-- to get the final state.  This takes place in several steps
--    1) align (DOMedian) parent final with node gapped (ie node preliminary)
--    2) propagate new gaps in aligned node preliminary to child gapped in node triple (snd and thd)
--       creating a 3-way alignment with parent final and child preliminary
--    3) apply appropriate get3way for the structure
-- The final is then returned--with gaps to be filtered afterwards
-- getDOFinal :: (FiniteBits a, GV.Vector v a) => v a -> (v a, v a, v a) -> CharInfo -> v a
getDOFinal :: CharInfo
           -> (SlimDynamicCharacter, WideDynamicCharacter, HugeDynamicCharacter)
           -> (SlimDynamicCharacter, WideDynamicCharacter, HugeDynamicCharacter)
           -> (SlimDynamicCharacter, WideDynamicCharacter, HugeDynamicCharacter)
getDOFinal charInfo parentFinal nodeGapped =
   let (a,b,c,_) = M.pairwiseDO charInfo parentFinal nodeGapped
       parentNodeChar = (a,b,c)

       -- put "new" gaps into 2nd and thd gapped fields of appropriate seqeunce type
       gappedFinal = makeGappedLeftRight charInfo parentNodeChar nodeGapped
   in
   gappedFinal


-- | makeGappedLeftRight takes an alignment parent character and original node character and inserts "new" gaps into nodeCharcater
-- makeGappedLeftRight :: CharacterData -> CharacterData -> CharInfo -> CharacterData
-- makeGappedLeftRight gappedLeftRight nodeChar charInfo =
makeGappedLeftRight :: CharInfo
                    -> (SlimDynamicCharacter, WideDynamicCharacter, HugeDynamicCharacter)
                    -> (SlimDynamicCharacter, WideDynamicCharacter, HugeDynamicCharacter)
                    -> (SlimDynamicCharacter, WideDynamicCharacter, HugeDynamicCharacter)
makeGappedLeftRight charInfo gappedLeftRight nodeChar  =
   let localCharType = charType charInfo
   in
   if localCharType `elem` [SlimSeq, NucSeq] then
      let (parentGapped, leftChildGapped, rightChildGapped) = TW.addGapsToChildren  (fst3 gappedLeftRight) (fst3 nodeChar)
          newFinalGapped = M.getFinal3WaySlim (slimTCM charInfo) parentGapped leftChildGapped rightChildGapped
      in
      (M.makeDynamicCharacterFromSingleVector newFinalGapped, mempty, mempty)

   else if localCharType `elem` [AminoSeq, WideSeq] then
      let (parentGapped, leftChildGapped, rightChildGapped) = TW.addGapsToChildren  (snd3 gappedLeftRight) (snd3 nodeChar)
          newFinalGapped = M.getFinal3WayWideHuge (wideTCM charInfo) parentGapped leftChildGapped rightChildGapped
      in
      (mempty, M.makeDynamicCharacterFromSingleVector newFinalGapped, mempty)

   else if localCharType == HugeSeq then
      let (parentGapped, leftChildGapped, rightChildGapped) = TW.addGapsToChildren  (thd3 gappedLeftRight) (thd3 nodeChar)
          newFinalGapped = M.getFinal3WayWideHuge (hugeTCM charInfo)  parentGapped leftChildGapped rightChildGapped
      in
      (mempty, mempty, M.makeDynamicCharacterFromSingleVector newFinalGapped)

   else error ("Unrecognized character type: " ++ show localCharType)



-- |  additivePreorder assignment takes preliminary triple of child (= current vertex) and
-- final states of parent to create preorder final states of child
additivePreorder :: (V.Vector (Int, Int), V.Vector (Int, Int), V.Vector (Int, Int)) -> V.Vector (Int, Int) ->  V.Vector (Int, Int)
additivePreorder (leftChild, nodePrelim, rightChild) parentFinal =
   if null nodePrelim then mempty
   else
      V.zipWith4 makeAdditiveCharacterFinal nodePrelim leftChild rightChild parentFinal

-- |  nonAdditivePreorder assignment takes preliminary triple of child (= current vertex) and
-- final states of parent to create preorder final states of child
nonAdditivePreorder :: (V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector) -> V.Vector BV.BitVector -> V.Vector BV.BitVector
nonAdditivePreorder (leftChild, nodePrelim, rightChild) parentFinal =
   if null nodePrelim then mempty
   else
      V.zipWith4 makeNonAdditiveCharacterFinal nodePrelim leftChild rightChild parentFinal


-- | matrixPreorder assigment akes preliminary matrix states of child (= current vertex) and
-- final states of parent to create preorder final states of child
-- th eboolean says whether the node is a 'left' node or right based on bitvetor label
matrixPreorder :: Bool -> V.Vector (V.Vector MatrixTriple) -> V.Vector (V.Vector MatrixTriple) -> V.Vector (V.Vector MatrixTriple)
matrixPreorder isLeft nodePrelim parentFinal =
   if null nodePrelim then mempty
   else
      V.zipWith (makeMatrixCharacterFinal isLeft)  nodePrelim parentFinal


-- | makeAdditiveCharacterFinal takes vertex preliminary, and child preliminary states with well as parent final state
-- and constructs final state assignment
makeAdditiveCharacterFinal :: (Int, Int) -> (Int, Int) -> (Int, Int) -> (Int, Int) -> (Int, Int)
makeAdditiveCharacterFinal nodePrelim leftChild rightChild parentFinal =
   -- From Wheeler (20012) after Goloboff (1993)
   let interNodeParent = intervalIntersection nodePrelim parentFinal
   in
   -- trace (show inData) (
   -- Rule 1
   if interNodeParent == Just parentFinal then
      -- trace ("R1 " ++ show parentFinal)
      parentFinal
   -- Rule 2
   else if isJust ((leftChild `intervalUnion` rightChild) `intervalIntersection` parentFinal) then
      let xFactor = ((leftChild `intervalUnion` rightChild) `intervalUnion` nodePrelim) `intervalIntersection` parentFinal
      in
      if isNothing xFactor then error ("I don't think this should happen in makeAdditiveCharacterFinal" ++ show (nodePrelim, leftChild, rightChild, parentFinal))
      else
         if isJust (fromJust xFactor `intervalIntersection` nodePrelim) then
            -- trace ("R2a " ++ show (fromJust xFactor))
            fromJust xFactor
         else
            -- trace ("Rb " ++ show (lciClosest (fromJust xFactor) nodePrelim))
            lciClosest (fromJust xFactor) nodePrelim

   -- Rule 3
   else
      let unionLR = leftChild `intervalUnion` rightChild
          closestPtoA = stateFirstClosestToSecond nodePrelim parentFinal
          closestULRtoA = stateFirstClosestToSecond unionLR parentFinal
      in
      -- trace ("R3 " ++ show (min closestPtoA closestULRtoA, max closestPtoA closestULRtoA))
      (min closestPtoA closestULRtoA, max closestPtoA closestULRtoA)
   -- )

-- | stateFirstClosestToSecond takes teh states of the first interval and finds the state wiht smallest distance
-- to either state in the second
 -- assumes a <= b, x<= y
stateFirstClosestToSecond :: (Int, Int) -> (Int, Int) -> Int
stateFirstClosestToSecond (a,b) (x,y) =
   let distASecond
         | x > b = x - a
         | y < a = a - y
         | otherwise = error ("I don't think this should happen in makeAdditiveCharacterFinal" ++ show (a,b,x,y))
       distBSecond
         | x > b = x - b
         | y < a = b - y
         | otherwise = error ("I don't think this should happen in makeAdditiveCharacterFinal" ++ show (a,b,x,y))
   in
   if distASecond <= distBSecond then a
   else b

-- | lciClosest returns the "largest closed interval" between the first interval
-- and the closest state in the second interval
 -- assumes a <= b, x<= y
lciClosest :: (Int, Int) -> (Int, Int) -> (Int, Int)
lciClosest (a,b) (x,y)
  | x > b = (a,x)
  | y < a = (y,b)
  | otherwise = error ("I don't think this should happen in lciClosest" ++ show (a,b,x,y))

 -- | intervalIntersection is bit-analogue intersection for additive character operations
 -- takes two intervals and returnas range intersection
 -- Nothing signifies an empty intersection
 -- assumes a <= b, x<= y
intervalIntersection :: (Int, Int) -> (Int, Int) -> Maybe (Int, Int)
intervalIntersection (a,b) (x,y) =
   let newPair = (max a x, min b y)
   in
   if max a x > min b y then Nothing
   else Just newPair


-- | intervalUnion is bit-analogue union for additive character operations
-- takes two intervals and returnas union
 -- assumes a <= b, x<= y
intervalUnion :: (Int, Int) -> (Int, Int) -> (Int, Int)
intervalUnion (a,b) (x,y) = (min a x, max b y)


-- | makeNonAdditiveCharacterFinal takes vertex preliminary, and child preliminary states with well as parent final state
-- and constructs final state assignment
makeNonAdditiveCharacterFinal :: BV.BitVector -> BV.BitVector-> BV.BitVector-> BV.BitVector -> BV.BitVector
makeNonAdditiveCharacterFinal nodePrelim leftChild rightChild parentFinal =
   -- From Wheeler (2012) after Fitch (1971)
   -- trace (show inData) (
   if BV.isZeroVector ((complement nodePrelim) .&. parentFinal) then
      --trace ("R1 " ++ show parentFinal)
      parentFinal
   else if (BV.isZeroVector (leftChild .&. rightChild)) then
      --trace ("R2 " ++ show (nodePrelim .|. parentFinal))
      nodePrelim .|. parentFinal
   else
      -- trace ("R3 " ++ show (nodePrelim .|.  (leftChild .&. parentFinal) .|. (rightChild .&. parentFinal)))
      nodePrelim .|. (parentFinal .&. (leftChild .|. rightChild))
   -- )

-- | makeMatrixCharacterFinal vertex preliminary and parent final state
-- and constructs final state assignment
-- really just tracks the states on a traceback and sets the cost to maxBound :: Int for states not in the traceback
-- path
-- Bool for left right node
makeMatrixCharacterFinal :: Bool -> V.Vector MatrixTriple -> V.Vector MatrixTriple -> V.Vector MatrixTriple
makeMatrixCharacterFinal isLeft nodePrelim parentFinal =
   let numStates = length nodePrelim
       stateIndexList = V.fromList [0..(numStates - 1)]
       (stateCostList, stateLeftChildList, stateRightChildList) = V.unzip3 parentFinal
       (_, prelimStateLeftChildList, prelimStateRightChildList) = V.unzip3 nodePrelim
       allThree = if isLeft then V.zip3 stateCostList stateLeftChildList stateIndexList
                  else V.zip3 stateCostList stateRightChildList stateIndexList
       bestParentThree = V.filter ((/= (maxBound :: StateCost)). fst3) allThree
       bestPrelimStates = L.sort $ L.nub $ concatMap snd3 bestParentThree
       allFour = V.zipWith3 (setCostsAndStates bestPrelimStates)  prelimStateLeftChildList prelimStateRightChildList stateIndexList
       finalBestTriple = V.filter ((/= (maxBound :: StateCost)).fst3) allFour
   in
   finalBestTriple

-- | setCostsAndStates takes a list of states that are in teh set of 'best' and a four-tuple
-- of a matrix triple annd a fourth field of the state index
-- if the state is in the list of `best' indices it is kept and not if it isn't
setCostsAndStates :: [Int] -> [ChildStateIndex] -> [ChildStateIndex] -> Int -> (StateCost, [ChildStateIndex], [ChildStateIndex])
setCostsAndStates bestPrelimStates leftChildState rightChildStates stateIndex =
   if stateIndex `elem` bestPrelimStates then (stateIndex, leftChildState, rightChildStates)
   else (maxBound :: StateCost, leftChildState, rightChildStates)


-- | setMinCostStatesMatrix  sets the cost of non-minimal cost states to maxBounnd :: StateCost (Int)
setMinCostStatesMatrix ::  Int -> V.Vector (V.Vector MatrixTriple) ->  V.Vector (V.Vector MatrixTriple)
setMinCostStatesMatrix numStates inStateVect =
    let outStates = V.filter ((/= (maxBound :: StateCost)).fst3) <$> fmap (nonMinCostStatesToMaxCost (V.fromList [0.. (numStates - 1)])) inStateVect
    in
    outStates

-- | nonMinCostStatesToMaxCost takes an individual pair of minimum state cost and matrix character triple
-- returns a new character with the states cost either the minium value or maxBound if not
-- this only applied at root or leaf--other vertices minimum costs may not be part of the
-- miniumm cost assignment, but may be useful heuristically
nonMinCostStatesToMaxCost :: V.Vector StateCost -> V.Vector MatrixTriple -> V.Vector MatrixTriple
nonMinCostStatesToMaxCost stateIndexList tripleVect =
   let minStateCost = V.minimum $ fmap fst3 tripleVect
       result = V.zipWith (modifyStateCost minStateCost) tripleVect stateIndexList
   in
   -- trace ((show stateIndexList) ++ " " ++ (show $ V.zip tripleVect stateIndexList))
   result
      where
         modifyStateCost d (a,b,c) e = if a == d then (e,b,c)
                                       else (maxBound :: StateCost ,b,c)

-- | setFinalToPreliminaryStates takes VertexBlockData and sets the final values to Preliminary
setFinalToPreliminaryStates :: VertexBlockData -> VertexBlockData
setFinalToPreliminaryStates inVertBlockData =
    if V.null inVertBlockData then mempty
    else
        fmap setBlockFinalToPrelim inVertBlockData

-- | setBlockFinalToPrelim sets characters in a block final values to Preliminary
setBlockFinalToPrelim :: V.Vector CharacterData -> V.Vector CharacterData
setBlockFinalToPrelim inCharVect =
    if V.null inCharVect then mempty
    else fmap setFinalToPrelimCharacterData inCharVect

-- | setFinalFinalToPrelimCharacterData takes a single chartcater and sets final values to Preliminary
setFinalToPrelimCharacterData :: CharacterData -> CharacterData
setFinalToPrelimCharacterData inChar =
    inChar {                 stateBVFinal       = snd3 $ stateBVPrelim inChar
                           , rangeFinal         = snd3 $ rangePrelim inChar
                           , matrixStatesFinal  = matrixStatesPrelim inChar
                           , slimAlignment      = slimGapped inChar
                           , slimFinal          = slimPrelim inChar
                           , slimIAFinal        = snd3 $ slimIAPrelim inChar

                           , wideAlignment      = wideGapped inChar
                           , wideFinal          = widePrelim inChar
                           , wideIAFinal        = snd3 $ wideIAPrelim inChar

                           , hugeAlignment      = hugeGapped inChar
                           , hugeFinal          = hugePrelim inChar
                           , hugeIAFinal        = snd3 $ hugeIAPrelim inChar

                           , alignedSlimFinal  = snd3 $ alignedSlimPrelim inChar
                           , alignedWideFinal  = snd3 $ alignedWidePrelim inChar
                           , alignedHugeFinal  = snd3 $ alignedHugePrelim inChar

                           , packedNonAddFinal = snd3 $ packedNonAddPrelim inChar
            }

-- | setPreliminaryToFinalStates takes VertexBlockData and sets the Preliminary states to final values
setPreliminaryToFinalStates :: VertexBlockData -> VertexBlockData
setPreliminaryToFinalStates inVertBlockData =
    if V.null inVertBlockData then mempty
    else
        fmap setBlockPrelimToFinal inVertBlockData

-- | setBlockPrelimToFinal sets characters in a block preliminary data to final
setBlockPrelimToFinal :: V.Vector CharacterData -> V.Vector CharacterData
setBlockPrelimToFinal inCharVect =
    if V.null inCharVect then mempty
    else fmap setPrelimToFinalCharacterData inCharVect

-- | setPrelimToFinalCharacterData takes a single chartcater and sets preliminary values to final
setPrelimToFinalCharacterData :: CharacterData -> CharacterData
setPrelimToFinalCharacterData inChar =
    inChar {                 stateBVPrelim      = (stateBVFinal inChar, stateBVFinal inChar, stateBVFinal inChar)
                           , rangePrelim        = (rangeFinal inChar, rangeFinal inChar, rangeFinal inChar)
                           , matrixStatesPrelim  = matrixStatesFinal inChar
                           -- , slimAlignment      = slimGapped inChar
                           , slimGapped          = (slimFinal inChar, slimFinal inChar, slimFinal inChar)
                           , slimIAPrelim        = (slimIAFinal inChar,slimIAFinal inChar, slimIAFinal inChar)

                           -- , wideAlignment      = wideGapped inChar
                           , wideGapped          = (wideFinal inChar, wideFinal inChar, wideFinal inChar)
                           , wideIAPrelim        = (wideIAFinal inChar, wideIAFinal inChar, wideIAFinal inChar)

                           -- , hugeAlignment      = hugeGapped inChar
                           , hugeGapped          = (hugeFinal inChar, hugeFinal inChar, hugeFinal inChar)
                           , hugeIAPrelim        = (hugeIAFinal inChar, hugeIAFinal inChar, hugeIAFinal inChar)

                           , alignedSlimPrelim  = (alignedSlimFinal inChar, alignedSlimFinal inChar, alignedSlimFinal inChar)
                           , alignedWidePrelim  = (alignedWideFinal inChar, alignedWideFinal inChar, alignedWideFinal inChar)
                           , alignedHugePrelim  = (alignedHugeFinal inChar, alignedHugeFinal inChar, alignedHugeFinal inChar)

                           , packedNonAddPrelim = (packedNonAddFinal inChar, packedNonAddFinal inChar, packedNonAddFinal inChar)
            }

