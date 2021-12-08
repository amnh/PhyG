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

{-
ToDo:
   Add parallel optimization overblocks and characters?
-}


module GraphOptimization.PreOrderFunctions  ( createFinalAssignmentOverBlocks
                                            , get2WaySlim
                                            , get2WayWideHuge
                                            , getFinal3WaySlim
                                            , getFinal3WayWideHuge
                                            , preOrderTreeTraversal
                                            , getBlockCostPairsFinal
                                            , setFinalToPreliminaryStates
                                            ) where

import           Bio.DynamicCharacter
import Data.Alphabet
import           Types.Types
import GeneralUtilities
import qualified DirectOptimization.PreOrder as DOP
import qualified Data.Vector                                                 as V
import qualified Data.BitVector.LittleEndian as BV
import           Data.Maybe
import qualified Data.List as L
import Debug.Trace
import qualified Data.TCM.Dense as TCMD
import qualified Data.MetricRepresentation as MR
import qualified Data.Vector.Generic                                         as GV
import Data.Bits
import qualified Data.Vector.Storable         as SV
import           Foreign.C.Types             (CUInt)
import qualified GraphOptimization.Medians as M
import qualified Utilities.LocalGraph as LG
import qualified SymMatrix                   as S
import qualified Graphs.GraphOperations as GO


-- | preOrderTreeTraversal takes a preliminarily labelled PhylogeneticGraph
-- and returns a full labels with 'final' assignments based on character decorated graphs
-- created postorder (5th of 6 fields).  
-- the preorder states are creted by traversing the traversal DecoratedGraphs in the 5th filed of PhylogeneticGraphs
-- these are by block and character,  Exact charcaters are vectors of standard characters and each seqeunce (non-exact) 
-- has its own traversal graph. These should be treees here (could be forests later) and should all have same root (number taxa)
-- but worth checking to make sure.
-- these were creted by "splitting" them after preorder 

-- For sequence charcaters (slim/wide/huge) final states are created either by DirectOptimization or ImpliedAlignment
-- if DO--then does a  median between parent andchgaps wherre needed--then  doing a 3way state assignmeent filteringgaps out
-- if IA--then a separate post and preorder pass are donne on the slim/wide/huge/AI fields to crete full IA assignments 
-- that are then filtered of gaps and assigned to th efinal fields

-- The final states are propagated back to the second field
-- DecoratedGraph of the full Phylogenetic graph--which does NOT have the character based preliminary assignments
-- ie postorder--since those are traversal specific
-- the character specific decorated graphs have appropriate post and pre-order assignments
-- the traversal begins at the root (for a tree) and proceeds to leaves.
preOrderTreeTraversal :: GlobalSettings -> AssignmentMethod -> Bool -> Int -> PhylogeneticGraph -> PhylogeneticGraph
preOrderTreeTraversal inGS finalMethod hasNonExact rootIndex inPGraph@(inSimple, inCost, inDecorated, blockDisplayV, blockCharacterDecoratedVV, inCharInfoVV) =
    --trace ("PreO: " ++ (show finalMethod) ++ " " ++ (show $ fmap (fmap charType) inCharInfoVV)) (
    if LG.isEmpty (thd6 inPGraph) then emptyPhylogeneticGraph
    else
        -- trace ("In PreOrder\n" ++ "Simple:\n" ++ (LG.prettify inSimple) ++ "Decorated:\n" ++ (LG.prettify $ GO.convertDecoratedToSimpleGraph inDecorated) ++ "\n" ++ (GFU.showGraph inDecorated)) (
        -- mapped recursive call over blkocks, later characters
        let -- preOrderBlockVect = fmap doBlockTraversal $ Debug.debugVectorZip inCharInfoVV blockCharacterDecoratedVV
            preOrderBlockVect = V.zipWith (doBlockTraversal finalMethod rootIndex) inCharInfoVV blockCharacterDecoratedVV

            -- if final non-exact states determined by IA then perform passes and assignments of final and final IA fields
            -- always do IA pass--but only assign to final if finalMethod == ImpliedAlignment
            preOrderBlockVect' = -- if (finalMethod == ImpliedAlignment) && hasNonExact then V.zipWith makeIAAssignments finalMethod preOrderBlockVect inCharInfoVV
                                 if hasNonExact then V.zipWith (makeIAAssignments finalMethod rootIndex) preOrderBlockVect inCharInfoVV
                                 else preOrderBlockVect

            fullyDecoratedGraph = assignPreorderStatesAndEdges inGS finalMethod rootIndex preOrderBlockVect' inCharInfoVV inDecorated
        in
        {-
        let blockPost = GO.showDecGraphs blockCharacterDecoratedVV
            blockPre = GO.showDecGraphs preOrderBlockVect
        in
        trace ("BlockPost:\n" ++ blockPost ++ "BlockPre:\n" ++ blockPre ++ "After Preorder\n" ++  (LG.prettify $ GO.convertDecoratedToSimpleGraph fullyDecoratedGraph))
        -}
        (inSimple, inCost, fullyDecoratedGraph, blockDisplayV, preOrderBlockVect, inCharInfoVV)
        -- )

-- | makeIAAssignments takes the vector of vector of character trees and (if) slim/wide/huge
-- does an additional pot and pre order pass to assign IA fileds and final fields in slim/wide/huge
makeIAAssignments :: AssignmentMethod -> Int -> V.Vector DecoratedGraph -> V.Vector CharInfo -> V.Vector DecoratedGraph
makeIAAssignments finalMethod rootIndex = V.zipWith (makeCharacterIA finalMethod rootIndex)

-- | makeCharacterIA takes an individual character postorder tree and if non-exact perform post and preorder IA passes
-- and assignment to final field in slim/wide/huge
makeCharacterIA :: AssignmentMethod -> Int -> DecoratedGraph -> CharInfo -> DecoratedGraph
makeCharacterIA finalMethod rootIndex inGraph charInfo =
    if charType charInfo `notElem` nonExactCharacterTypes then inGraph
    else
        let postOrderIATree = postOrderIA inGraph charInfo [(rootIndex, fromJust $ LG.lab inGraph rootIndex)]
            preOrderIATree = preOrderIA postOrderIATree finalMethod charInfo $ zip [(rootIndex, fromJust $ LG.lab postOrderIATree rootIndex)] [(rootIndex, fromJust $ LG.lab postOrderIATree rootIndex)] 
        in
        preOrderIATree

-- | postOrderIA performs a post-order IA pass assigning leaf preliminary states
-- from the "alignment" fields and setting HTU preliminary by calling the apropriate 2-way
-- matrix
postOrderIA :: DecoratedGraph -> CharInfo -> [LG.LNode VertexInfo] -> DecoratedGraph
postOrderIA inGraph charInfo inNodeList =
    if null inNodeList then inGraph
    else
        let inNode@(nodeIndex, nodeLabel) = head inNodeList
            (inNodeEdges, outNodeEdges) = LG.getInOutEdges inGraph nodeIndex
            inCharacter = V.head $ V.head $ vertData nodeLabel
            nodeType' = GO.getNodeType inGraph nodeIndex
        in

        -- trace ("POIA Node: " ++ (show nodeIndex) ++ " " ++ (show $ nodeType nodeLabel) ++ " " ++ (show  $ fmap fst inNodeList)) (
        -- checking sanity of data
        if V.null $ vertData nodeLabel then error "Null vertData in postOrderIA"
        else if V.null $ V.head $ vertData nodeLabel then error "Null vertData data in postOrderIA"

        -- leaf take assignment from alignment field
        else if nodeType nodeLabel  == LeafNode then
            {-
            let newCharacter
                  | characterType `elem` [SlimSeq, NucSeq] =
                     inCharacter { slimIAPrelim = slimAlignment inCharacter'
                                 , slimIAFinal = extractMediansGapped $ slimAlignment inCharacter'
                                 , slimFinal = extractMedians $ slimAlignment inCharacter'
                                 }
                  | characterType `elem` [WideSeq, AminoSeq] =
                     inCharacter { wideIAPrelim = wideAlignment inCharacter'
                                 , wideIAFinal = extractMediansGapped $ wideAlignment inCharacter'
                                 , wideFinal = extractMedians $ wideAlignment inCharacter'
                                 }
                  | characterType == HugeSeq =
                     inCharacter { hugeIAPrelim = hugeAlignment inCharacter'
                                 , hugeIAFinal = extractMediansGapped $ hugeAlignment inCharacter'
                                 , hugeFinal = extractMedians $ hugeAlignment inCharacter'
                                 }
                  | otherwise = error ("Unrecognized character type " ++ show characterType)
                newLabel = nodeLabel  {vertData = V.singleton (V.singleton newCharacter)}
                newGraph = LG.insEdges (inNodeEdges ++ outNodeEdges) $ LG.insNode (nodeIndex, newLabel) $ LG.delNode nodeIndex inGraph
            in
            -}
            -- trace ("PostOLeaf: " ++ (show nodeIndex) ++ " " ++ (show $ slimFinal $ V.head $ V.head $ vertData nodeLabel))
            postOrderIA inGraph charInfo (tail inNodeList)
            -- postOrderIA newGraph charInfo (tail inNodeList)

        -- HTU take create assignment from children
        else
            let childNodes = LG.labDescendants inGraph inNode
                childTree = postOrderIA inGraph charInfo childNodes
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
                else if V.null $ V.head $ vertData childLabel then error "Null vertData data in postOrderIA"
                else
                    let newLabel = nodeLabel  {vertData = V.singleton (V.singleton childCharacter), nodeType = nodeType'}
                        newGraph = LG.insEdges (inNodeEdges ++ outNodeEdges) $ LG.insNode (nodeIndex, newLabel) $ LG.delNode nodeIndex childTree
                    in
                    -- trace ("PostO1Child: " ++ (show nodeIndex) ++ " " ++ (show $ slimFinal childCharacter))
                    postOrderIA newGraph charInfo (tail inNodeList)

            -- two children 
            else
                let childIndices = fmap fst childNodes
                    childlabels = fmap (fromJust . LG.lab childTree) childIndices
                    childCharacters = fmap vertData childlabels
                    leftChar = V.head $ V.head $ head childCharacters
                    rightChar = V.head $ V.head $ last childCharacters
                    newCharacter = makeIAPrelimCharacter charInfo inCharacter leftChar rightChar
                    newLabel = nodeLabel  {vertData = V.singleton (V.singleton newCharacter), nodeType = nodeType'}
                    newGraph = LG.insEdges (inNodeEdges ++ outNodeEdges) $ LG.insNode (nodeIndex, newLabel) $ LG.delNode nodeIndex childTree
                in
                -- trace ("PostO2hildren: " ++ (show nodeIndex) ++ " " ++ (show $ slimFinal newCharacter) ++ " " ++ (show $ nodeType newLabel)) -- ++ " From: " ++ (show childlabels))
                postOrderIA newGraph charInfo (tail inNodeList)
            -- )
    -- )

-- | makeIAPrelimCharacter takes two characters and performs 2-way assignment 
-- based on character type and nodeChar--only IA fields are modified
makeIAPrelimCharacter :: CharInfo -> CharacterData -> CharacterData -> CharacterData -> CharacterData
makeIAPrelimCharacter charInfo nodeChar leftChar rightChar =
     let characterType = charType charInfo
     in
     if characterType `elem` [SlimSeq, NucSeq] then
        let prelimChar = get2WaySlim (slimTCM charInfo) (extractMediansGapped $ slimIAPrelim leftChar) (extractMediansGapped $ slimIAPrelim rightChar)
        in
        -- trace ("MPC: " ++ (show prelimChar) ++ "\nleft: " ++ (show $ extractMediansGapped $ slimIAPrelim leftChar) ++ "\nright: " ++ (show $ extractMediansGapped $ slimIAPrelim rightChar))
        nodeChar {slimIAPrelim = (extractMediansGapped $ slimIAPrelim leftChar, prelimChar,  extractMediansGapped $ slimIAPrelim rightChar)}
     else if characterType `elem` [WideSeq, AminoSeq] then
        let prelimChar = get2WayWideHuge (wideTCM charInfo) (extractMediansGapped $ wideIAPrelim leftChar) (extractMediansGapped $ wideIAPrelim rightChar)
        in
        nodeChar {wideIAPrelim = (extractMediansGapped $ wideIAPrelim leftChar, prelimChar, extractMediansGapped $ wideIAPrelim rightChar)}
     else if characterType == HugeSeq then
        let prelimChar = get2WayWideHuge (hugeTCM charInfo) (extractMediansGapped $ hugeIAPrelim leftChar) (extractMediansGapped $ hugeIAPrelim rightChar)
        in
        nodeChar {hugeIAPrelim = (extractMediansGapped $ hugeIAPrelim leftChar, prelimChar, extractMediansGapped $ hugeIAPrelim rightChar)}
     else error ("Unrecognized character type " ++ show characterType)

-- | makeIAFinalharacter takes two characters and performs 2-way assignment 
-- based on character type and nodeChar--only IA fields are modified
makeIAFinalCharacter :: AssignmentMethod -> CharInfo -> CharacterData -> CharacterData -> CharacterData -> CharacterData -> CharacterData
makeIAFinalCharacter finalMethod charInfo nodeChar parentChar leftChar rightChar =
     let characterType = charType charInfo
     in
     if characterType `elem` [SlimSeq, NucSeq] then
        let finalIAChar = getFinal3WaySlim (slimTCM charInfo) (slimIAFinal parentChar) (extractMediansGapped $ slimIAPrelim leftChar) (extractMediansGapped $ slimIAPrelim rightChar)
            finalChar =  if finalMethod == ImpliedAlignment then extractMedians $ M.makeDynamicCharacterFromSingleVector finalIAChar
                         else slimFinal nodeChar 
        in
        nodeChar { slimIAFinal = finalIAChar
                 , slimFinal = finalChar
                 }
     else if characterType `elem` [WideSeq, AminoSeq] then
        let finalIAChar = getFinal3WayWideHuge (wideTCM charInfo) (wideIAFinal parentChar) (extractMediansGapped $ wideIAPrelim leftChar) (extractMediansGapped $ wideIAPrelim rightChar)
            finalChar = if finalMethod == ImpliedAlignment then extractMedians $ M.makeDynamicCharacterFromSingleVector finalIAChar
                        else wideFinal nodeChar 
        in
        nodeChar { wideIAFinal = finalIAChar
                 , wideFinal = finalChar
                 }
     else if characterType == HugeSeq then
        let finalIAChar = getFinal3WayWideHuge (hugeTCM charInfo) (hugeIAFinal parentChar) (extractMediansGapped $ hugeIAPrelim leftChar) (extractMediansGapped $ hugeIAPrelim rightChar)
            finalChar = if finalMethod == ImpliedAlignment then extractMedians $ M.makeDynamicCharacterFromSingleVector finalIAChar
                        else hugeFinal nodeChar 
        in
        nodeChar { hugeIAFinal = finalIAChar
                 , hugeFinal = finalChar
                 }
     else error ("Unrecognized character type " ++ show characterType)



-- | preOrderIA performs a pre-order IA pass assigning via the apropriate 3-way matrix
-- the "final" fields are also set by filtering out gaps and 0.
preOrderIA :: DecoratedGraph -> AssignmentMethod -> CharInfo -> [(LG.LNode VertexInfo, LG.LNode VertexInfo)] -> DecoratedGraph
preOrderIA inGraph finalMethod charInfo inNodePairList =
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
        else if nodeType nodeLabel == LeafNode then preOrderIA inGraph finalMethod charInfo (tail inNodePairList)

        else if nodeType nodeLabel == RootNode then
            let newCharacter
                  | characterType `elem` [SlimSeq, NucSeq] =
                     inCharacter { slimIAFinal = extractMediansGapped $ slimIAPrelim inCharacter'
                                 -- , slimFinal = extractMedians $ slimIAPrelim  inCharacter'
                                 }
                  | characterType `elem` [WideSeq, AminoSeq] =
                     inCharacter { wideIAFinal = extractMediansGapped $ wideIAPrelim inCharacter'
                                 -- , wideFinal = extractMedians $ wideIAPrelim inCharacter'
                                 }
                  | characterType == HugeSeq =
                     inCharacter { hugeIAFinal = extractMediansGapped $ hugeIAPrelim inCharacter'
                                 -- , hugeFinal = extractMedians $ hugeIAPrelim inCharacter'
                                 }
                  | otherwise = error ("Unrecognized character type " ++ show characterType)
                newLabel = nodeLabel  {vertData = V.singleton (V.singleton newCharacter)}
                newGraph = LG.insEdges (inNodeEdges ++ outNodeEdges) $ LG.insNode (nodeIndex, newLabel) $ LG.delNode nodeIndex inGraph
                parentNodeList = replicate (length childNodes) (nodeIndex, newLabel)
            in
            -- trace ("PreIARoot: " ++ (show nodeIndex) ++ " IAFinal: " ++ (show $ slimIAFinal newCharacter) ++ " Final: " ++ (show $ slimFinal newCharacter)) 
            preOrderIA newGraph finalMethod charInfo (tail inNodePairList ++ zip childNodes parentNodeList)

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
                      | otherwise = error ("Unrecognized character type " ++ show characterType)

                    newLabel = nodeLabel  {vertData = V.singleton (V.singleton newCharacter)}
                    newGraph = LG.insEdges (inNodeEdges ++ outNodeEdges) $ LG.insNode (nodeIndex, newLabel) $ LG.delNode nodeIndex inGraph
                    parentNodeList = replicate (length childNodes) (nodeIndex, newLabel)
                in
                -- trace ("PreIANet: " ++ (show nodeIndex) ++ " IAFinal: " ++ (show $ slimIAFinal newCharacter) ++ " Final: " ++ (show $ slimFinal newCharacter)) 
                preOrderIA newGraph finalMethod charInfo (tail inNodePairList ++ zip childNodes parentNodeList)

        -- 2 children, make 3-way 
        else
            let childLabels = fmap snd childNodes
                leftChar = V.head $ V.head $ vertData $ head childLabels
                rightChar = V.head $ V.head $ vertData $ last childLabels
                finalCharacter = makeIAFinalCharacter finalMethod charInfo inCharacter parentCharacter leftChar rightChar

                newLabel = nodeLabel  {vertData = V.singleton (V.singleton finalCharacter)}
                newGraph = LG.insEdges (inNodeEdges ++ outNodeEdges) $ LG.insNode (nodeIndex, newLabel) $ LG.delNode nodeIndex inGraph
                parentNodeList = replicate (length childNodes) (nodeIndex, newLabel)
            in
            -- trace ("PreIATree: " ++ (show nodeIndex) ++ " IAFinal: " ++ (show $ slimIAFinal finalCharacter) ++ " Final: " ++ (show $ slimFinal finalCharacter)) 
            preOrderIA newGraph finalMethod charInfo (tail inNodePairList ++ zip childNodes parentNodeList)

        -- )

-- | doBlockTraversal takes a block of postorder decorated character trees character info  
-- could be moved up preOrderTreeTraversal, but like this for legibility
doBlockTraversal :: AssignmentMethod -> Int -> V.Vector CharInfo -> V.Vector DecoratedGraph -> V.Vector DecoratedGraph
doBlockTraversal finalMethod rootIndex inCharInfoV traversalDecoratedVect =
    --trace ("BlockT:" ++ (show $ fmap charType inCharInfoV)) 
    V.zipWith (doCharacterTraversal finalMethod rootIndex) inCharInfoV traversalDecoratedVect

-- | doCharacterTraversal performs preorder traversal on single character tree
-- with single charInfo
-- this so each character can be independently "rooted" for optimal traversals.
doCharacterTraversal :: AssignmentMethod -> Int -> CharInfo -> DecoratedGraph -> DecoratedGraph
doCharacterTraversal finalMethod rootIndex inCharInfo inGraph =
    -- find root--index should = number of leaves 
    --trace ("charT:" ++ (show $ charType inCharInfo)) (
    let isolateNodeList = LG.getIsolatedNodes inGraph
        (_, leafVertexList, _, _)  = LG.splitVertexList inGraph
        inEdgeList = LG.labEdges inGraph
    in
    -- remove these two lines if working
    -- if rootIndex /=  length leafVertexList then error ("Root index not =  number leaves in doCharacterTraversal" ++ show (rootIndex, length leafVertexList))
    -- else
        -- root vertex, repeat of label info to avoid problem with zero length zip later, second info ignored for root
        let rootLabel = fromJust $ LG.lab inGraph rootIndex
            rootFinalVertData = createFinalAssignmentOverBlocks finalMethod RootNode (vertData rootLabel) (vertData rootLabel) inCharInfo True False
            rootChildren =LG.labDescendants inGraph (rootIndex, rootLabel)

            -- left / right to match post-order
            rootChildrenBV = fmap (bvLabel . snd) rootChildren
            rootChildrenIsLeft
              | length rootChildrenBV == 1 = [True]
              | head rootChildrenBV > (rootChildrenBV !! 1) = [False, True]
              | otherwise = [True, False]
            newRootNode = (rootIndex, rootLabel {vertData = rootFinalVertData})
            rootChildrenPairs = zip3 rootChildren (replicate (length rootChildren) newRootNode) rootChildrenIsLeft
            upDatedNodes = makeFinalAndChildren finalMethod inGraph rootChildrenPairs [newRootNode] inCharInfo

            -- update isolated nodes with final == preliminary as with root nodes (and leaves, but without postorder logic)
            updatedIsolateNodes = fmap (updateIsolatedNode finalMethod inCharInfo) isolateNodeList
        in
        -- hope this is the most efficient way since all nodes have been remade
        -- trace (U.prettyPrintVertexInfo $ snd newRootNode)
        LG.mkGraph (upDatedNodes ++ updatedIsolateNodes) inEdgeList
        --)

-- | updateIsolatedNode updates the final states of an isolated node as if it were a root with final=preliminary
-- states without preorder logic as in regular leaves
-- NB IA length won't match if compared since not in graph
updateIsolatedNode :: AssignmentMethod -> CharInfo -> LG.LNode VertexInfo -> LG.LNode VertexInfo
updateIsolatedNode finalMethod inCharInfo (inNodeIndex, inNodeLabel) = 
    let newVertData = createFinalAssignmentOverBlocks finalMethod RootNode (vertData inNodeLabel) (vertData inNodeLabel) inCharInfo True False
    in
    (inNodeIndex, inNodeLabel {vertData = newVertData})

-- | makeFinalAndChildren takes a graph, list of pairs of (labelled nodes,parent node) to make final assignment and a liss of updated nodes
-- the input nodes are relabelled by preorder functions and added to the list of processed nodes and recursed to their children
-- nodes are retuned in reverse order at they are made--need to check if this will affect graph identity or indexing in fgl
makeFinalAndChildren :: AssignmentMethod
                     -> DecoratedGraph
                     -> [(LG.LNode VertexInfo, LG.LNode VertexInfo, Bool)]
                     -> [LG.LNode VertexInfo]
                     -> CharInfo
                     -> [LG.LNode VertexInfo]
makeFinalAndChildren finalMethod inGraph nodesToUpdate updatedNodes inCharInfo =
    --trace ("mFAC:" ++ (show $ charType inCharInfo)) (
    if null nodesToUpdate then updatedNodes
    else
        let (firstNode, firstParent, isLeft) = head nodesToUpdate
            firstLabel = snd firstNode
            firstNodeType = nodeType firstLabel
            firstVertData = vertData firstLabel
            firstParentVertData = vertData $ snd firstParent
            firstChildren = LG.labDescendants inGraph firstNode

            -- this OK with one or two children
            firstChildrenBV = fmap (bvLabel . snd) firstChildren
            firstChildrenIsLeft
              | length firstChildrenBV == 1 = [True]
              | head firstChildrenBV > (firstChildrenBV !! 1) = [False, True]
              | otherwise = [True, False]
            firstFinalVertData = createFinalAssignmentOverBlocks finalMethod firstNodeType firstVertData firstParentVertData inCharInfo isLeft (length firstChildren == 1)
            newFirstNode = (fst firstNode, firstLabel {vertData = firstFinalVertData})
            childrenPairs = zip3 firstChildren (replicate (length firstChildren) newFirstNode) firstChildrenIsLeft
        in
        -- trace (U.prettyPrintVertexInfo $ snd newFirstNode)
        makeFinalAndChildren finalMethod inGraph (childrenPairs ++ tail nodesToUpdate) (newFirstNode : updatedNodes) inCharInfo
        --)

-- | assignPreorderStatesAndEdges takes a postorder decorated graph (should be but not required) and propagates 
-- preorder character states from individual character trees.  Exact characters (Add, nonAdd, matrix) postorder
-- states should be based on the outgroup rooted tree.  
-- root should be median of finals of two descendets--for non-exact based on final 'alignments' field with gaps filtered
-- postorder assignment and preorder will be out of whack--could change to update with correponding postorder
-- but that would not allow use of base decorated graph for incremental optimization (which relies on postorder assignments) in other areas
-- optyion code ikn there to set root final to outgropu final--but makes thigs scewey in matrix character and some pre-order assumptions
assignPreorderStatesAndEdges :: GlobalSettings -> AssignmentMethod -> Int -> V.Vector (V.Vector DecoratedGraph) -> V.Vector (V.Vector CharInfo) -> DecoratedGraph  -> DecoratedGraph
assignPreorderStatesAndEdges inGS finalMethd rootIndex preOrderBlockTreeVV inCharInfoVV inGraph =
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
            newEdgeList = if (graphType inGS == Tree || graphType inGS == HardWired) then fmap (updateEdgeInfoTree finalMethd inCharInfoVV (V.fromList newNodeList)) postOrderEdgeList
                          else fmap (updateEdgeInfoSoftWired finalMethd inCharInfoVV blockTreePairVV rootIndex) postOrderEdgeList
        in
        -- make new graph
        -- LG.mkGraph newNodeList' newEdgeList
        LG.mkGraph newNodeList newEdgeList
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
-- only updating preorder assignment--except for root, that is needed to draw final state for brnach lengths
updateCharacter :: CharacterData -> CharacterData -> CharType  -> CharacterData
updateCharacter postOrderCharacter preOrderCharacter localCharType
  | localCharType == Add =
    postOrderCharacter { rangeFinal = rangeFinal preOrderCharacter }
  | localCharType == NonAdd =
    postOrderCharacter { stateBVFinal = stateBVFinal preOrderCharacter }
  | localCharType == Matrix =
    postOrderCharacter { matrixStatesFinal = matrixStatesFinal preOrderCharacter }
  | localCharType == SlimSeq || localCharType == NucSeq =
    postOrderCharacter { slimAlignment = slimAlignment preOrderCharacter
                  , slimFinal = slimFinal preOrderCharacter
                  , slimIAFinal = slimIAFinal preOrderCharacter
              }
  | localCharType == WideSeq || localCharType == AminoSeq =
    postOrderCharacter { wideAlignment = wideAlignment preOrderCharacter
                  , wideFinal = wideFinal preOrderCharacter
                  , wideIAFinal = wideIAFinal preOrderCharacter
                  }
  | localCharType == HugeSeq =
    postOrderCharacter { hugeAlignment = hugeAlignment preOrderCharacter
                  , hugeFinal = hugeFinal preOrderCharacter
                  , hugeIAFinal = hugeIAFinal preOrderCharacter
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


-- | getEdgeBlockWeightSoftWired takes a block of charcter trees and maps character distances if edge exists in block tree
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


-- | getEdgeVerts retuns vertex labels if edge in vect or if a virtual edge including root
getEdgeVerts :: Int -> Int -> Int -> V.Vector (LG.LNode VertexInfo) -> V.Vector (LG.LEdge EdgeInfo) -> Maybe (VertexInfo, VertexInfo)
getEdgeVerts uNode vNode rootIndex nodeVect edgeVect =
    if edgeInVect (uNode, vNode) edgeVect then Just (snd $ nodeVect V.! uNode, snd $ nodeVect V.! vNode)
    else if (edgeInVect (rootIndex, uNode) edgeVect) && (edgeInVect (rootIndex, vNode) edgeVect) then Just (snd $ nodeVect V.! uNode, snd $ nodeVect V.! vNode)
    else Nothing

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

-- | getEdgeWeight takes a preorder decorated decorated graph and an edge and gets the weight information for that edge
-- basically a min/max distance between the two
getEdgeWeight :: AssignmentMethod -> V.Vector (V.Vector CharInfo) -> V.Vector (LG.LNode VertexInfo) -> (Int, Int) -> (VertexCost, VertexCost)
getEdgeWeight finalMethod inCharInfoVV nodeVector (uNode, vNode) =
    if V.null nodeVector then error "Empty node list in getEdgeWeight"
    else
        trace ("GEW: " ++ (show $ fmap fst nodeVector) ++ " " ++ (show (uNode, vNode))) (
        let uNodeInfo = vertData $ snd $ nodeVector V.! uNode
            vNodeInfo = vertData $ snd $ nodeVector V.! vNode
            blockCostPairs = V.zipWith3 (getBlockCostPairsFinal finalMethod) uNodeInfo vNodeInfo inCharInfoVV
            minCost = sum $ fmap fst blockCostPairs
            maxCost = sum $ fmap snd blockCostPairs
        in

        (minCost, maxCost)
        )

-- | getBlockCostPairsFinal takes a block of two nodes and character infomation and returns the min and max block branch costs
getBlockCostPairsFinal :: AssignmentMethod -> V.Vector CharacterData -> V.Vector CharacterData -> V.Vector CharInfo -> (VertexCost, VertexCost)
getBlockCostPairsFinal finalMethod uNodeCharDataV vNodeCharDataV charInfoV =
    let characterCostPairs = V.zipWith3 (getCharacterDistFinal finalMethod) uNodeCharDataV vNodeCharDataV charInfoV
        minCost = sum $ fmap fst characterCostPairs
        maxCost = sum $ fmap snd characterCostPairs
    in
    (minCost, maxCost)

-- | getCharacterDistFinal takes a pair of characters and character type, returning the minimum and maximum character distances
-- for sequence charcaters this is based on slim/wide/hugeAlignment field, hence all should be O(n) in num characters/sequence length
getCharacterDistFinal :: AssignmentMethod -> CharacterData -> CharacterData -> CharInfo -> (VertexCost, VertexCost)
getCharacterDistFinal finalMethod uCharacter vCharacter charInfo =
    let thisWeight = weight charInfo
        thisMatrix = costMatrix charInfo
        thisCharType = charType charInfo
    in
    if thisCharType == Add then
        let minCost = localCost (M.intervalAdd thisWeight uCharacter vCharacter)
            maxDiff = V.sum $ V.zipWith maxIntervalDiff  (rangeFinal uCharacter) (rangeFinal vCharacter)
            maxCost = thisWeight * fromIntegral maxDiff
        in
        (minCost, maxCost)


    else if thisCharType == NonAdd then
        let minCost = localCost (M.interUnion thisWeight uCharacter vCharacter)
            maxDiff = length $ V.filter (==False) $ V.zipWith (==) (stateBVFinal uCharacter) (stateBVFinal vCharacter)
            maxCost = thisWeight * fromIntegral maxDiff
        in
        (minCost, maxCost)

    else if thisCharType == Matrix then
        let minMaxListList= V.zipWith (minMaxMatrixDiff thisMatrix)  (fmap fst3 <$> matrixStatesFinal uCharacter) (fmap fst3 <$> matrixStatesFinal vCharacter)
            minDiff = V.sum $ fmap fst minMaxListList
            maxDiff = V.sum $ fmap snd minMaxListList
            minCost = thisWeight * fromIntegral minDiff
            maxCost = thisWeight * fromIntegral maxDiff
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
                                zipWith  (generalSequenceDiff thisMatrix (length thisMatrix)) (GV.toList newU) (GV.toList newV)
                             else zipWith  (generalSequenceDiff thisMatrix (length thisMatrix)) (GV.toList $ slimIAFinal uCharacter) (GV.toList $ slimIAFinal vCharacter)
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
                                zipWith  (generalSequenceDiff thisMatrix (length thisMatrix)) (GV.toList newU) (GV.toList newV)
                             else GV.toList $ GV.zipWith (generalSequenceDiff thisMatrix (length thisMatrix))  (wideIAFinal uCharacter) (wideIAFinal vCharacter)
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
                                zipWith  (generalSequenceDiff thisMatrix (length thisMatrix)) (GV.toList newU) (GV.toList newV)
                             else GV.toList $ GV.zipWith (generalSequenceDiff thisMatrix (length thisMatrix)) (hugeIAFinal uCharacter) (hugeIAFinal vCharacter)
            (minDiff, maxDiff) = unzip minMaxDiffList
            minCost = thisWeight * fromIntegral (sum minDiff)
            maxCost = thisWeight * fromIntegral (sum maxDiff)
        in
        (minCost, maxCost)

    else error ("Character type unimplemented : " ++ show thisCharType)
{-
-- | zero2Gap converts a '0' or no bits set to gap (indel) value 
zero2Gap :: (FiniteBits a) => a -> a
zero2Gap inVal = if inVal == zeroBits then bit gapIndex
                 else inVal

-- | zero2GapWide converts a '0' or no bits set to gap (indel) value 
zero2GapWide :: Word64 -> Word64 -> Word64
zero2GapWide gapChar inVal = if inVal == zeroBits  then bit gapIndex
                         else inVal

-- | zero2GapBV converts a '0' or no bits set to gap (indel) value 
zero2GapBV :: BV.BitVector -> BV.BitVector -> BV.BitVector
zero2GapBV gapChar inVal = if inVal == zeroBits then bit gapIndex
                         else inVal
-}

-- | maxIntervalDiff takes two ranges and gets the maximum difference between the two based on differences
-- in upp and lower ranges.
maxIntervalDiff :: (Int, Int)-> (Int, Int) -> Int
maxIntervalDiff (a,b) (x,y) =
    let upper = max b y - min b y
        lower = max a x - min a x
    in
    max upper lower

-- | minMaxMatrixDiff takes twovetors of states and calculates the minimum and maximum state differnce cost 
-- between the two
minMaxMatrixDiff :: S.Matrix Int -> V.Vector Int -> V.Vector Int -> (Int, Int)
minMaxMatrixDiff localCostMatrix uStatesV vStatesV =
    let statePairs = (V.toList uStatesV, V.toList vStatesV)
        cartesianPairs = cartProdPair statePairs
        costList = fmap (localCostMatrix S.!) cartesianPairs
    in
    --trace (show cartesianPairs  ++ " " ++ show costList) 
    (minimum costList, maximum costList)



-- | generalSequenceDiff  takes two sequnce elemental bit types and retuns min and max integer 
-- cost differences using matrix values
-- if value has no bits on--it is set to 0th bit on for GAP
generalSequenceDiff :: (FiniteBits a) => S.Matrix Int -> Int -> a -> a -> (Int, Int)
generalSequenceDiff thisMatrix numStates uState vState =
    -- trace ("GSD: " ++ (show (numStates, uState, vState))) (
    let uState' = if uState == zeroBits then bit gapIndex else uState
        vState' = if vState == zeroBits then bit gapIndex else vState
        uStateList = fmap snd $ filter ((== True).fst) $ zip (fmap (testBit uState') [0.. numStates - 1]) [0.. numStates - 1]
        vStateList = fmap snd $ filter ((== True).fst) $ zip (fmap (testBit vState') [0.. numStates - 1]) [0.. numStates - 1]
        uvCombinations = cartProd uStateList vStateList
        costOfPairs = fmap (thisMatrix S.!) uvCombinations
    in
    -- trace ("GSD: " ++ (show uStateList) ++ " " ++ (show vStateList) ++ " min " ++ (show $ minimum costOfPairs) ++ " max " ++ (show $  maximum costOfPairs))
    (minimum costOfPairs, maximum costOfPairs)
    -- )

-- | createFinalAssignment takes vertex data (child or current vertex) and creates the final 
-- assignment from parent (if not root or leaf) and 'child' ie current vertex
-- if root or leaf preliminary is assigned to final
   -- need to watch zipping for missing sequence data
-- this creates the IA during preorder from which final assignments are contructed
-- via addition post and preorder passes on IA fields.
createFinalAssignmentOverBlocks :: AssignmentMethod
                                -> NodeType
                                -> VertexBlockData
                                -> VertexBlockData
                                -> CharInfo
                                -> Bool
                                -> Bool
                                -> VertexBlockData
createFinalAssignmentOverBlocks finalMethod childType childBlockData parentBlockData charInfo isLeft isOutDegree1 =
   -- if root or leaf final assignment <- preliminary asssignment
   V.zipWith (assignFinal finalMethod childType isLeft charInfo isOutDegree1) childBlockData parentBlockData

  -- | assignFinal takes a vertex type and single block of zip3 of child info, parent info, and character type 
-- to create pre-order assignments
assignFinal :: AssignmentMethod -> NodeType -> Bool -> CharInfo -> Bool -> V.Vector CharacterData -> V.Vector CharacterData -> V.Vector CharacterData
assignFinal finalMethod childType isLeft charInfo isOutDegree1 = V.zipWith (setFinal finalMethod childType isLeft charInfo isOutDegree1)

-- | setFinal takes a vertex type and single character of zip3 of child info, parent info, and character type 
-- to create pre-order assignments
   -- | setFinalHTU takes a single character and its parent and sets the final state to prelim based 
-- on character info. 
-- non exact charcaters are vectors of characters of same type
-- this does the same things for sequence types, but also 
-- performs preorder logic for exact characters
setFinal :: AssignmentMethod -> NodeType -> Bool -> CharInfo -> Bool -> CharacterData-> CharacterData -> CharacterData
setFinal finalMethod childType isLeft charInfo isOutDegree1 childChar parentChar =
   let localCharType = charType charInfo
       symbolCount = toEnum $ length $ costMatrix charInfo :: Int
   in
   -- Three cases, Root, leaf, HTU
   --trace ("set final:" ++ (show isLeft) ++ " " ++ (show isOutDegree1) ++ " " ++ (show $ slimAlignment parentChar) ++ " " 
   --   ++ (show $ slimGapped parentChar) ++ " " ++ (show $ slimGapped childChar)) (
   if childType == RootNode then

      if localCharType == Add then
         childChar {rangeFinal = snd3 $ rangePrelim childChar}

      else if localCharType == NonAdd then childChar {stateBVFinal = snd3 $ stateBVPrelim childChar}

      else if localCharType == Matrix then childChar {matrixStatesFinal = setMinCostStatesMatrix (fromEnum symbolCount) (localCostVect childChar) (matrixStatesPrelim childChar)}

      -- need to set both final and alignment for sequence characters
      else if (localCharType == SlimSeq) || (localCharType == NucSeq) then
         let finalAssignment' = extractMedians $ slimGapped childChar
         in
         -- trace ("SF Root: " ++ (show finalAssignment') ++ "\n" ++ (show $ slimGapped childChar))
         childChar {slimFinal = finalAssignment', slimAlignment = slimGapped childChar}

      else if (localCharType == WideSeq) || (localCharType == AminoSeq) then
         let finalAssignment' = extractMedians $ wideGapped childChar
         in
         childChar {wideFinal = finalAssignment', wideAlignment = wideGapped childChar}

      else if localCharType == HugeSeq then
         let finalAssignment' = extractMedians $ hugeGapped childChar
         in
         childChar {hugeFinal = finalAssignment', hugeAlignment = hugeGapped childChar}

      else error ("Unrecognized/implemented character type: " ++ show localCharType)

   else if childType == LeafNode then
      -- since leaf no neeed to precess final alignment fields for sequence characters
      if localCharType == Add then childChar {rangeFinal = snd3 $ rangePrelim childChar}

      else if localCharType == NonAdd then childChar {stateBVFinal = snd3 $ stateBVPrelim childChar}

      else if localCharType == Matrix then
         childChar {matrixStatesFinal = setMinCostStatesMatrix (fromEnum symbolCount) (V.replicate  (fromEnum symbolCount) 0) (matrixStatesPrelim childChar)}

      -- need to set both final and alignment for sequence characters
      else if (localCharType == SlimSeq) || (localCharType == NucSeq) then
         let finalAlignment = DOP.preOrderLogic isLeft (slimAlignment parentChar) (slimGapped parentChar) (slimGapped childChar)
             finalAssignment' = extractMedians finalAlignment
         in
         --trace ("Leaf " ++ show (slimPrelim childChar, slimPrelim childChar, finalAlignment, slimGapped childChar, slimAlignment parentChar))
         -- trace ("SF Leaf: " ++ (show isLeft) ++ " FA': " ++  (show finalAssignment') ++ "\nFA: " ++ (show finalAlignment) ++ "\nSAP: " ++
         --    (show $ slimAlignment parentChar) ++ "\nSGP: " ++ (show $ slimGapped parentChar) ++ "\nSGC: " ++ (show $ slimGapped childChar))
         childChar {slimFinal = finalAssignment', slimAlignment = finalAlignment, slimIAPrelim = finalAlignment, slimIAFinal = extractMediansGapped $ finalAlignment}

      else if (localCharType == WideSeq) || (localCharType == AminoSeq) then
         let finalAlignment = DOP.preOrderLogic isLeft (wideAlignment parentChar) (wideGapped parentChar) (wideGapped childChar)
             finalAssignment' = extractMedians finalAlignment
         in
         childChar {wideFinal = finalAssignment', wideAlignment = finalAlignment, wideIAPrelim = finalAlignment, wideIAFinal = extractMediansGapped $ finalAlignment}

      else if localCharType == HugeSeq then
         let finalAlignment = DOP.preOrderLogic isLeft (hugeAlignment parentChar) (hugeGapped parentChar) (hugeGapped childChar)
             finalAssignment' = extractMedians finalAlignment
         in
         childChar {hugeFinal = finalAssignment', hugeAlignment = finalAlignment, hugeIAPrelim = finalAlignment, hugeIAFinal = extractMediansGapped $ finalAlignment}

      else error ("Unrecognized/implemented character type: " ++ show localCharType)

   else if childType == TreeNode && not isOutDegree1 then

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

      else if localCharType == Matrix then
         -- add logic for pre-order
         let finalAssignment' = matrixPreorder isLeft (matrixStatesPrelim childChar) (matrixStatesFinal parentChar)
         in
         childChar {matrixStatesFinal = finalAssignment'}

      -- need to set both final and alignment for sequence characters
      else if (localCharType == SlimSeq) || (localCharType == NucSeq) then
         let finalGapped = DOP.preOrderLogic isLeft (slimAlignment parentChar) (slimGapped parentChar) (slimGapped childChar)
             finalAssignmentDO = if finalMethod == DirectOptimization then
                                    let parentFinalDC = M.makeDynamicCharacterFromSingleVector (slimFinal parentChar)
                                        parentFinal = (parentFinalDC, mempty, mempty)
                                        -- parentGapped = (slimGapped parentChar, mempty, mempty)
                                        childGapped = (slimGapped childChar, mempty, mempty)
                                        finalAssignmentDOGapped = fst3 $ getDOFinal charInfo parentFinal  childGapped
                                    in
                                    extractMedians finalAssignmentDOGapped
                                    -- really could/should be mempty since overwritten by IA later
                                 else extractMedians finalGapped
         in
         -- trace ("SF Tree: " ++ (show isLeft) ++ " FA': " ++  (show finalAssignment) ++ "\nFA: " ++ (show finalGapped) ++ "\nSAP: " ++
         --   (show $ slimAlignment parentChar) ++ "\nSGP: " ++ (show $ slimGapped parentChar) ++ "\nSGC: " ++ (show $ slimGapped childChar))
         --trace ("SF Tree: " ++ (show finalGapped))
         childChar {slimFinal = finalAssignmentDO, slimAlignment = finalGapped}
         -- For debugging IA/DO bug 
         -- childChar {slimFinal = mempty, slimAlignment = finalGapped}

      else if (localCharType == WideSeq) || (localCharType == AminoSeq) then
         let finalGapped = DOP.preOrderLogic isLeft (wideAlignment parentChar) (wideGapped parentChar) (wideGapped childChar)
             finalAssignmentDO = if finalMethod == DirectOptimization then
                                    let parentFinalDC = M.makeDynamicCharacterFromSingleVector (wideFinal parentChar)
                                        parentFinal = (mempty, parentFinalDC, mempty)
                                        --parentGapped = (mempty, wideGapped parentChar, mempty)
                                        childGapped = (mempty, wideGapped childChar, mempty)
                                        finalAssignmentDOGapped = snd3 $ getDOFinal charInfo parentFinal  childGapped
                                    in
                                    extractMedians finalAssignmentDOGapped
                                 else extractMedians finalGapped
         in
         childChar {wideFinal = finalAssignmentDO, wideAlignment = finalGapped}

      else if localCharType == HugeSeq then
         let finalGapped = DOP.preOrderLogic isLeft (hugeAlignment parentChar) (hugeGapped parentChar) (hugeGapped childChar)
             finalAssignmentDO = if finalMethod == DirectOptimization then
                                    let parentFinalDC = M.makeDynamicCharacterFromSingleVector (hugeFinal parentChar)
                                        parentFinal = (mempty, mempty, parentFinalDC)
                                        -- parentGapped = (mempty, mempty, hugeGapped parentChar)
                                        childGapped = (mempty, mempty, hugeGapped childChar)
                                        finalAssignmentDOGapped = thd3 $ getDOFinal charInfo parentFinal  childGapped
                                    in
                                    extractMedians finalAssignmentDOGapped
                                 else extractMedians finalGapped
         in
         childChar {hugeFinal = finalAssignmentDO, hugeAlignment = finalGapped}

      else error ("Unrecognized/implemented character type: " ++ show localCharType)

   -- display tree indegree=outdegree=1
   -- since display trees here--indegree should be one as well
   else if isOutDegree1 then
      -- trace ("InOut1 preorder") (
      if localCharType == Add then
         -- add logic for pre-order
         let lFinalAssignment = rangeFinal parentChar
         in
         childChar {rangeFinal = lFinalAssignment}

      else if localCharType == NonAdd then
         -- add logic for pre-order
         let lFinalAssignment = stateBVFinal parentChar
         in
         childChar {stateBVFinal = lFinalAssignment}

      else if localCharType == Matrix then
         -- add logic for pre-order
         let lFinalAssignment = matrixStatesFinal parentChar
         in
         childChar {matrixStatesFinal = lFinalAssignment}

      -- need to set both final and alignment for sequence characters
      else if (localCharType == SlimSeq) || (localCharType == NucSeq) then
         -- trace ("SF Net: " ++ (show $ slimAlignment parentChar))
         childChar { slimFinal = slimFinal parentChar
                   , slimAlignment = slimAlignment parentChar
                   , slimIAPrelim = slimIAPrelim parentChar
                   , slimIAFinal = slimFinal parentChar}

      else if (localCharType == WideSeq) || (localCharType == AminoSeq) then
         childChar { wideFinal = wideFinal parentChar
                   , wideAlignment = wideAlignment parentChar
                   , wideIAPrelim = wideIAPrelim parentChar
                   , wideIAFinal = wideFinal parentChar}

      else if localCharType == HugeSeq then
         childChar { hugeFinal = hugeFinal parentChar
                   , hugeAlignment =  hugeAlignment parentChar
                   , hugeIAPrelim = hugeIAPrelim parentChar
                   , hugeIAFinal = hugeFinal parentChar}

      else error ("Unrecognized/implemented character type: " ++ show localCharType)
      -- )

   else error ("Node type should not be here (pre-order on tree node only): " ++ show  childType)
   -- )

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


-- | makeGappedLeftRight takes an alignment parent charcater and original node character and inserts "new" gaps into nodeCharcater
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
      let (parentGapped, leftChildGapped, rightChildGapped) = addGapsToChildren  (fst3 gappedLeftRight) (fst3 nodeChar)
          newFinalGapped = getFinal3WaySlim (slimTCM charInfo) parentGapped leftChildGapped rightChildGapped
      in
      (M.makeDynamicCharacterFromSingleVector newFinalGapped, mempty, mempty)

   else if localCharType `elem` [AminoSeq, WideSeq] then
      let (parentGapped, leftChildGapped, rightChildGapped) = addGapsToChildren  (snd3 gappedLeftRight) (snd3 nodeChar)
          newFinalGapped = getFinal3WayWideHuge (wideTCM charInfo) parentGapped leftChildGapped rightChildGapped
      in
      (mempty, M.makeDynamicCharacterFromSingleVector newFinalGapped, mempty)

   else if localCharType == HugeSeq then
      let (parentGapped, leftChildGapped, rightChildGapped) = addGapsToChildren  (thd3 gappedLeftRight) (thd3 nodeChar)
          newFinalGapped = getFinal3WayWideHuge (hugeTCM charInfo)  parentGapped leftChildGapped rightChildGapped
      in
      (mempty, mempty, M.makeDynamicCharacterFromSingleVector newFinalGapped)

   else error ("Unrecognized character type: " ++ show localCharType)

-- | addGapsToChildren pads out "new" gaps based on identity--if not identical--adds a gap based on cost matrix size
addGapsToChildren :: (FiniteBits a, GV.Vector v a) => (v a, v a, v a) -> (v a, v a, v a) -> (v a, v a, v a)
addGapsToChildren  (reGappedParentFinal, _, reGappedNodePrelim) (gappedLeftChild, gappedNodePrelim, gappedRightChild) =
   let (reGappedLeft, reGappedRight) = slideRegap reGappedNodePrelim gappedNodePrelim gappedLeftChild gappedRightChild mempty mempty
   in
   if (GV.length reGappedParentFinal /= GV.length reGappedLeft) || (GV.length reGappedParentFinal /= GV.length reGappedRight) then error ("Vectors not same length "
      ++ show (GV.length reGappedParentFinal, GV.length reGappedLeft, GV.length reGappedRight))
   else (reGappedParentFinal, reGappedLeft, reGappedRight)

-- | slideRegap takes two version of same vectors (1st and snd) one with additional gaps and if the two aren't equal then adds gaps 
-- to the 3rd and 4th input vectors
slideRegap :: (FiniteBits a, GV.Vector v a) => v a -> v a -> v a -> v a -> [a] -> [a] -> (v a, v a)
slideRegap reGappedNode gappedNode gappedLeft gappedRight newLeftList newRightList =
   -- trace ("SRG: " ++ (show (GV.length reGappedNode, GV.length gappedNode))) (
   if GV.null reGappedNode then (GV.fromList $ reverse newLeftList, GV.fromList $ reverse newRightList)
   else
      let firstRGN = GV.head reGappedNode
          firstGN = GV.head  gappedNode
      in

      -- gap in reGappedNode, null gappedNode is gap at end of reGappedNode
      -- can copmplete the remainder of the slide as gaps only
      if GV.null gappedNode then
        let gapList = replicate  (GV.length reGappedNode) (bit gapIndex)
        in
        (GV.fromList $ reverse (gapList ++ newLeftList), GV.fromList $ reverse (gapList ++ newRightList))

      else if firstRGN /= firstGN then
         let gap = bit gapIndex
         in
         slideRegap (GV.tail reGappedNode) gappedNode gappedLeft gappedRight (gap : newLeftList) (gap : newRightList)
      
      -- no "new gap"  
      else -- if firstRGN == firstGN then
         slideRegap (GV.tail reGappedNode) (GV.tail gappedNode) (GV.tail gappedLeft) (GV.tail gappedRight) (GV.head gappedLeft : newLeftList) (GV.head gappedRight : newRightList)
    -- )


-- | getFinal3Way takes parent final assignment (including indel characters) and descendent
-- preliminary gapped assingment from postorder and creates a gapped final assignment based on 
-- minimum cost median for the three inputs.  THis is done to preserve the ImpliedAlignment
-- information to create a final assingment with out an additional DO call to keep the 
-- creation linear in sequence length.  Since gaps remain--they must be filtered when output or 
-- used as true final sequence assignments using M.createUngappedMedianSequence
getFinal3WaySlim :: TCMD.DenseTransitionCostMatrix -> SV.Vector CUInt -> SV.Vector CUInt -> SV.Vector CUInt -> SV.Vector CUInt
getFinal3WaySlim lSlimTCM parentFinal descendantLeftPrelim descendantRightPrelim =
   let newFinal = SV.zipWith3 (local3WaySlim lSlimTCM) parentFinal descendantLeftPrelim descendantRightPrelim
   in
   newFinal

-- | getFinal3WayWideHuge like getFinal3WaySlim but for wide and huge characters
getFinal3WayWideHuge :: (FiniteBits a, GV.Vector v a) => MR.MetricRepresentation a -> v a -> v a -> v a -> v a
getFinal3WayWideHuge whTCM parentFinal descendantLeftPrelim descendantRightPrelim =
   let newFinal = GV.zipWith3 (local3WayWideHuge whTCM) parentFinal descendantLeftPrelim descendantRightPrelim
   in
   newFinal

-- | local3WayWideHuge takes tripples for wide and huge sequence types and returns median
local3WayWideHuge :: (FiniteBits a) => MR.MetricRepresentation a -> a -> a -> a -> a
local3WayWideHuge lWideTCM b c d =
   let  -- b' = if b == zeroBits then gap else b
        -- c' = if c == zeroBits then gap else c
        -- d' = if d == zeroBits then gap else d
        (median, _) = MR.retreiveThreewayTCM lWideTCM b c d
   in
   -- trace ((show b) ++ " " ++ (show c) ++ " " ++ (show d) ++ " => " ++ (show median))
   median

-- | local3WaySlim takes triple of CUInt and retuns median
local3WaySlim :: TCMD.DenseTransitionCostMatrix -> CUInt -> CUInt -> CUInt -> CUInt
local3WaySlim lSlimTCM b c d =
 -- trace ("L3WS: " ++ (show (b,c,d))) (
 let  -- b' = if b == zeroBits then gap else b
      -- c' = if c == zeroBits then gap else c
      -- d' = if d == zeroBits then gap else d
 in
 -- trace ("L3WS: " ++ (show (b',c',d'))) (
 let (median, _) = TCMD.lookupThreeway lSlimTCM b c d
 in
 -- trace ("3way: " ++ (show b) ++ " " ++ (show c) ++ " " ++ (show d) ++ " => " ++ (show (median, cost)))
 median
 -- )

-- | get2WaySlim takes two slim vectors an produces a preliminary median
get2WaySlim :: TCMD.DenseTransitionCostMatrix -> SV.Vector CUInt -> SV.Vector CUInt -> SV.Vector CUInt
get2WaySlim lSlimTCM descendantLeftPrelim descendantRightPrelim =
   let median = SV.zipWith (local2WaySlim lSlimTCM) descendantLeftPrelim descendantRightPrelim
   in
   median

-- | local2WaySlim takes pair of CUInt and retuns median
local2WaySlim :: TCMD.DenseTransitionCostMatrix -> CUInt -> CUInt -> CUInt
local2WaySlim lSlimTCM b c =
 let  -- b' = if b == zeroBits then gap else (b `shiftL` 1) -- 2 * b -- temp fix for 2-way oddness
      -- c' = if c == zeroBits then gap else c
      -- (median, _) = TCMD.lookupPairwise lSlimTCM (b `shiftL` 1) c -- shiflt is a temp fix for 2-way oddness
      (median, _) = TCMD.lookupPairwise lSlimTCM b c
 in
 
 -- if b == zeroBits || c == zeroBits then trace ("Uh oh zero bits on") median
 -- trace ("Hello 2Way: " ++ (show b) ++ " " ++ (show c) ++ " " ++ " => " ++ (show (median, cost)))
 -- else median

 median

-- | get2WayWideHuge like get2WaySlim but for wide and huge characters
get2WayWideHuge :: (FiniteBits a, GV.Vector v a) => MR.MetricRepresentation a -> v a -> v a -> v a
get2WayWideHuge whTCM  descendantLeftPrelim descendantRightPrelim =
   let median = GV.zipWith (local2WayWideHuge whTCM) descendantLeftPrelim descendantRightPrelim
   in
   median

-- | local3WayWideHuge takes tripples for wide and huge sequence types and returns median
local2WayWideHuge :: (FiniteBits a) => MR.MetricRepresentation a -> a -> a -> a
local2WayWideHuge lWideTCM b c =
   let  -- b' = if b == zeroBits then gap else b
        -- c' = if c == zeroBits then gap else c
        (median, _) = MR.retreivePairwiseTCM lWideTCM b c
   in
   --trace ((show b) ++ " " ++ (show c) ++ " " ++ (show d) ++ " => " ++ (show median))
   median


-- |  additivePreorder assignment takes preliminary triple of child (= current vertex) and
-- final states of parent to create preorder final states of child
additivePreorder :: (V.Vector (Int, Int), V.Vector (Int, Int), V.Vector (Int, Int)) -> V.Vector (Int, Int) ->  V.Vector (Int, Int)
additivePreorder (nodePrelim, leftChild, rightChild) parentFinal =
   if null nodePrelim then mempty
   else
      V.zipWith4 makeAdditiveCharacterFinal nodePrelim leftChild rightChild parentFinal

-- |  nonAdditivePreorder assignment takes preliminary triple of child (= current vertex) and
-- final states of parent to create preorder final states of child
nonAdditivePreorder :: (V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector) -> V.Vector BV.BitVector -> V.Vector BV.BitVector
nonAdditivePreorder (nodePrelim, leftChild, rightChild) parentFinal =
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

-- | lciClosest returns thhe "largest closed interval" between the first interval
-- and the closest state in the second interval
 -- assumes a <= b, x<= y
lciClosest :: (Int, Int) -> (Int, Int) -> (Int, Int)
lciClosest (a,b) (x,y)
  | x > b = (a,x)
  | y < a = (y,b)
  | otherwise = error ("I don't think this should happen in lciClosest" ++ show (a,b,x,y))

 -- | intervalIntersection is bit-analogue intersection for additive character operations
 -- takes two intervals and returnas intesection
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

{-
-- | largestClosedInterval is the maximum interval created from two intervals
largestClosedInterval :: (Int, Int) -> (Int, Int) -> (Int, Int)
largestClosedInterval = intervalUnion
-}

-- | makeNonAdditiveCharacterFinal takes vertex preliminary, and child preliminary states with well as parent final state
-- and constructs final state assignment
makeNonAdditiveCharacterFinal :: BV.BitVector -> BV.BitVector-> BV.BitVector-> BV.BitVector -> BV.BitVector
makeNonAdditiveCharacterFinal nodePrelim leftChild rightChild parentFinal =
   -- From Wheeler (2012) after Fitch (1971)
   -- trace (show inData) (
   if (nodePrelim .&. parentFinal) == parentFinal then
      --trace ("R1 " ++ show parentFinal) 
      parentFinal
   else if BV.isZeroVector (leftChild .&. rightChild) && (leftChild .|. rightChild) == nodePrelim then
      --trace ("R2 " ++ show (nodePrelim .|. parentFinal)) 
      nodePrelim .|. parentFinal
   else
      -- trace ("R3 " ++ show (nodePrelim .|.  (leftChild .&. parentFinal) .|. (rightChild .&. parentFinal))) 
      nodePrelim .|.  (leftChild .&. parentFinal) .|. (rightChild .&. parentFinal)
   -- )

-- | makeMatrixCharacterFinal vertex preliminaryavnd parent final state
-- and constructs final state assignment
-- really just tracks the states on a traceback and sets the cost to maxBound:: Int for states not in the traceback
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
setMinCostStatesMatrix ::  Int -> V.Vector StateCost -> V.Vector (V.Vector MatrixTriple) ->  V.Vector (V.Vector MatrixTriple)
setMinCostStatesMatrix numStates inCostVect inStateVect =
    V.filter ((/= (maxBound :: StateCost)).fst3) <$> V.zipWith (nonMinCostStatesToMaxCost (V.fromList [0.. (numStates - 1)])) inCostVect inStateVect

-- | nonMinCostStatesToMaxCost takes an individual pair of minimum state cost and matrix character triple 
-- retiurns a new character with the states cost either the minium value or maxBound iof not
-- this only really useful at root--other vertices minimu costs may not be paert of the
-- miniumm cost assignment, but may be useful heuristically
nonMinCostStatesToMaxCost :: V.Vector StateCost -> StateCost -> V.Vector MatrixTriple -> V.Vector MatrixTriple
nonMinCostStatesToMaxCost stateIndexList minStateCost tripleVect =
   let result = V.zipWith (modifyStateCost minStateCost) tripleVect stateIndexList
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
        fmap setBlockPrelimToFinal inVertBlockData

-- | setBlockPrelimToFinal sets characyters in a block preliminary data to final
setBlockPrelimToFinal :: V.Vector CharacterData -> V.Vector CharacterData
setBlockPrelimToFinal inCharVect = 
    if V.null inCharVect then mempty 
    else fmap setPrelimToFinalCharacterData inCharVect

-- | setPrelimToFinalCharacterData takes a single chartcater and sets preliminary values to final
setPrelimToFinalCharacterData :: CharacterData -> CharacterData
setPrelimToFinalCharacterData inChar =
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
            }
