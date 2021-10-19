{- |
Module      :  PostOrderFunctions.hs
Description :  Module specifying post-order graph functions
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


module GraphOptimization.PostOrderFunctions  ( rerootPhylogeneticTree
                                             , rerootPhylogeneticTree'
                                             , rerootPhylogeneticNetwork
                                             , rerootPhylogeneticNetwork'
                                             , createVertexDataOverBlocks
                                             , divideDecoratedGraphByBlockAndCharacterTree
                                             , divideDecoratedGraphByBlockAndCharacterSoftWired
                                             ) where

import           Data.Bits                 ((.|.))
import           Data.Maybe
import qualified Data.Vector               as V
import qualified GraphOptimization.Medians as M
import           Types.Types
import qualified Utilities.LocalGraph      as LG
import qualified Utilities.Utilities       as U
import qualified Graphs.GraphOperations    as GO
-- import Debug.Debug
import Debug.Trace

-- | reOptimizeNodes takes a decorated graph and a list of nodes and reoptimizes (relabels)
-- them based on children in input graph
-- simple recursive since each node depends on children
-- remove check for debugging after it works
-- check for out-degree 1, doesn't matter for trees however.
reOptimizeNodes :: GraphType -> V.Vector (V.Vector CharInfo) -> DecoratedGraph -> [LG.LNode VertexInfo] -> DecoratedGraph
reOptimizeNodes localGraphType charInfoVectVect inGraph oldNodeList =
  if null oldNodeList then inGraph
  else
    -- make sure that nodes are optimized in correct order so that nodes are only reoptimized using updated children
    -- this should really not have to happen--order should be determined a priori
    let curNode@(curNodeIndex, curNodeLabel) = head oldNodeList
        nodeChildren = LG.descendants inGraph curNodeIndex  -- should be 1 or 2, not zero since all leaves already in graph
        foundCurChildern = filter (`elem` nodeChildren) $ fmap fst (tail oldNodeList)
    in
    if (not $ null foundCurChildern) then
      -- trace ("Current node " ++ (show curNodeIndex) ++ " has children " ++ (show nodeChildren) ++ " in optimize list (optimization order error)" ++ (show $ fmap fst $ tail oldNodeList))
      reOptimizeNodes localGraphType charInfoVectVect inGraph ((tail oldNodeList) ++ [curNode])
    
    else
        let leftChild = (head nodeChildren)
            rightChild = (last nodeChildren)
            -- leftChildLabel = fromJust $ LG.lab inGraph leftChild
            -- rightChildLabel = fromJust $ LG.lab inGraph rightChild

            -- this ensures that left/right choices are based on leaf BV for consistency and label invariance
            (leftChildLabel, rightChildLabel) = U.leftRightChildLabelBV (fromJust $ LG.lab inGraph leftChild, fromJust $ LG.lab inGraph rightChild)
            newVertexData = createVertexDataOverBlocksNonExact (vertData leftChildLabel) (vertData  rightChildLabel) charInfoVectVect []
            --newVertexData = createVertexDataOverBlocks  (vertData leftChildLabel) (vertData  rightChildLabel) charInfoVectVect []
        in
        {-
        --debug remove when not needed--checking to see if node should not be re optimized
        if (sort nodeChildren) == (sort $ V.toList $ children curnodeLabel) then
          trace ("Children for vertex unchanged " ++ (show curNodeIndex)
          reOptimizeNodes localGraphType charInfoVectVect inGraph (tail oldNodeList)
        else
        -}
        if localGraphType == Tree then 
           let newCost =  if (length nodeChildren < 2) then 0
                          else V.sum $ V.map (V.sum) $ V.map (V.map snd) newVertexData
               newVertexLabel = VertexInfo {  index = curNodeIndex
                                            -- this bit labelling incorect for outdegree = 1, need to prepend bits
                                            , bvLabel = (bvLabel leftChildLabel) .|. (bvLabel rightChildLabel)
                                            , parents = V.fromList $ LG.parents inGraph curNodeIndex
                                            , children = V.fromList nodeChildren
                                            , nodeType = nodeType curNodeLabel
                                            , vertName = vertName curNodeLabel
                                            , vertexResolutionData = mempty
                                            , vertData = if (length nodeChildren < 2) then vertData leftChildLabel
                                                         else V.map (V.map fst) newVertexData
                                            , vertexCost = newCost
                                            , subGraphCost = if (length nodeChildren < 2) then subGraphCost leftChildLabel
                                                             else (subGraphCost leftChildLabel) + (subGraphCost rightChildLabel) + newCost
                                            }
               -- this to add back edges deleted with nodes (undocumented but sensible in fgl)
               replacementEdges = (fst $ LG.getInOutEdges inGraph curNodeIndex) ++ (snd $ LG.getInOutEdges inGraph curNodeIndex)
               newGraph = LG.insEdges replacementEdges $ LG.insNode (curNodeIndex, newVertexLabel) $ LG.delNode curNodeIndex inGraph
            in
            --trace ("New vertexCost " ++ show newCost) --  ++ " lcn " ++ (show (vertData leftChildLabel, vertData rightChildLabel, vertData curnodeLabel)))
            reOptimizeNodes localGraphType charInfoVectVect newGraph (tail oldNodeList)
        
        else if  localGraphType == SoftWired then
            errorWithoutStackTrace ("Graph type yet implemented: " ++ show localGraphType)

        else  errorWithoutStackTrace ("Graph type unrecognized/not yet implemented: " ++ show localGraphType)



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
            firstBlockMedian = if (leftBlockLength == 0) then V.zip (V.head rightBlockData) (V.replicate rightBlockLength 0)
                               else if (rightBlockLength == 0) then V.zip (V.head leftBlockData) (V.replicate leftBlockLength 0)
                               else medianFunction (V.head leftBlockData) (V.head rightBlockData) (V.head blockCharInfoVect)
        in
        generalCreateVertexDataOverBlocks medianFunction (V.tail leftBlockData) (V.tail rightBlockData) (V.tail blockCharInfoVect) (firstBlockMedian : curBlockData)



-- | rerootPhylogeneticNetwork take a vertex index and reroots phylogenetic network
-- roots as for tree, but when a netwrok edge is chosen as root edge, its 'sister' network edge
-- is deleted (leaving two in=out=1 nodes), th ereooting takes place,
-- and the edge is re-added in both directions.  Nodes may change between
-- tree and network.  This is updated 
rerootPhylogeneticNetwork :: Int -> PhylogeneticGraph -> PhylogeneticGraph
rerootPhylogeneticNetwork rerootIndex inGraph@(inSimple, _, inDecGraph, blockDisplayForestVV, _, charInfoVectVect) = 
    if LG.isEmpty inSimple then inGraph
    else
        let originalRoots = LG.getRoots inSimple
            originalRootIndex = fst $ head originalRoots
            originalRootChildren = LG.descendants inSimple originalRootIndex
            newRootChildren = LG.descendants inSimple rerootIndex
            newRootParents = LG.parents inSimple rerootIndex
        in

        -- sanity check on graph
        if length originalRoots /= 1 then errorWithoutStackTrace ("Input graph has <>1 roots: " ++ (show originalRoots) ++ "\n" ++ (LG.prettify inSimple))

        -- would produce same root position
        else if rerootIndex `elem` originalRootChildren then inGraph

        -- OK to reroot
        else 
            --if length newRootParents > 1 then trace ("Root on network edge") inGraph
            --else 
                let newSimple = GO.rerootTree rerootIndex inSimple
                in
                -- currently just bags out on network branch
                if (null $ LG.getRoots newSimple) then 
                    trace ("RRN No roots:" ++ (show originalRootIndex) ++ " children " ++ (show originalRootChildren) ++ "-> " ++ (show rerootIndex) ++ " parents " 
                        ++ (show newRootParents) ++ " children " ++ (show newRootChildren) ++ "\nNewGraph:\n" ++ (LG.prettify newSimple)) 
                    inGraph
                else rerootPhylogeneticTree rerootIndex inGraph
              

-- | rerootPhylogeneticNetwork' flipped version of rerootPhylogeneticNetwork
rerootPhylogeneticNetwork' :: PhylogeneticGraph -> Int -> PhylogeneticGraph
rerootPhylogeneticNetwork' inGraph rerootIndex = rerootPhylogeneticNetwork rerootIndex inGraph

-- | rerootPhylogeneticTree' flipped version of rerootPhylogeneticTree
rerootPhylogeneticTree' :: PhylogeneticGraph -> Int -> PhylogeneticGraph
rerootPhylogeneticTree'  inGraph rerootIndex = rerootPhylogeneticTree rerootIndex inGraph

-- | rerootGraph takes a phylogenetic graph and reroots based on a vertex index (usually leaf outgroup)
--   if input is a forest then only roots the component that contains the vertex wil be rerooted
--   unclear how will effect network edges--will need to verify that does not create cycles
--   multi-rooted components (as opposed to forests) are unaffected with trace warning thrown
--   after checking for existing root and multiroots, should be O(n) where 'n is the length
--   of the path between the old and new root
--   the PhyloGenetic graph is (minimally) reoptimized along the spine of edegs that are redirected
--   in order that the root costr be correct.  The block display forest (components are always trees--for soft-wired
--   graphs only) is also rerooted, and
--   Character foci set to the new root edge
--   NB--Uses calls to rerootGraph since traversals are for different graphs so wouldn't save
--   much time by consolidating--also since labels are all different--can't re-use alot of info
--   from graph to graph.
--   NNB only deals with post-order states
rerootPhylogeneticTree :: Int -> PhylogeneticGraph -> PhylogeneticGraph
rerootPhylogeneticTree  rerootIndex inPhyGraph@(inSimple, _, inDecGraph, blockDisplayForestVV, _, charInfoVectVect) =
  if LG.isEmpty inSimple || LG.isEmpty inDecGraph then error "Empty graph in rerootPhylogeneticGraph"
  --else if inCost == 0 then error ("Input graph with cost zero--likely non decorated input graph in rerootPhylogeneticGraph\n" ++ (LG.prettify $ convertDecoratedToSimpleGraph inDecGraph))
  else
    let -- simple graph rerooted Boolean to specify that non-exact characters need NOT be reoptimized if affected
        newSimpleGraph = GO.rerootTree rerootIndex inSimple

        -- decorated graph Boolean to specify that non-exact characters need to be reoptimized if affected
        -- could just update with needges? from simple graph rerooting
        newDecGraph = GO.rerootTree rerootIndex inDecGraph

        -- reoptimize nodes here
        -- nodes on spine from new root to old root that needs to be reoptimized
        (nodesToOptimize, _) = LG.pathToRoot inDecGraph (rerootIndex, fromJust $ LG.lab inDecGraph rerootIndex)

        -- reversed because ll these node edges are reversed so preorder would be in reverse orientation
        -- this only reoptimizes non-exact characters since rerooting doesn't affect 'exact" character optimization'
        newDecGraph' = reOptimizeNodes Tree charInfoVectVect newDecGraph (reverse nodesToOptimize)

        -- sum of root costs on Decorated graph
        newGraphCost = sum $ fmap subGraphCost $ fmap snd $ LG.getRoots newDecGraph'

        -- rerooted diplay forests--don't care about costs--I hope (hence Bool False)
        newBlockDisplayForestVV = if V.null blockDisplayForestVV then mempty
                                  else fmap (fmap (GO.rerootTree rerootIndex)) blockDisplayForestVV

        in
        --trace ("=")
        -- Same root, so no need to redo
        if (length nodesToOptimize == 1) then inPhyGraph
        else
          {-
          trace ("To optimize:" ++ (show nodesToOptimize) ++ "\nOG " ++ (show inCost) ++ " :"
          ++ (LG.prettify inDecGraph) ++ "\nRRG:" ++ ((LG.prettify newDecGraph)) ++ "\nNG " ++ (show newGraphCost) ++ " :" ++ (LG.prettify newDecGraph')
          ++ "\nSG:" ++ (LG.prettify newSimpleGraph))
          -}
          -- (newSimpleGraph, newGraphCost, newDecGraph', newBlockDisplayForestVect, V.replicate (length charInfoVectVect) (V.singleton newDecGraph'), charInfoVectVect)
          (newSimpleGraph, newGraphCost, newDecGraph', newBlockDisplayForestVV, divideDecoratedGraphByBlockAndCharacterTree newDecGraph', charInfoVectVect)

-- | divideDecoratedGraphByBlockAndCharacterSoftWired takes a Vector of a list of DecoratedGraph
-- continaing a list of decorated tryees that are the display trees for that block
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

-- | divideDecoratedGraphByBlockAndCharacterTree takes a DecoratedGraph with (potentially) multiple blocks
-- and (potentially) multiple character per block and creates a Vector of Vector of Decorated Graphs
-- over blocks and characyets with the same graph, but only a single block and character for each graph
-- this to be used to create the "best" cost over alternate graph traversals
-- vertexCost and subGraphCost will be taken from characterData localcost/localcostVect and globalCost
divideDecoratedGraphByBlockAndCharacterTree :: DecoratedGraph -> V.Vector (V.Vector DecoratedGraph)
divideDecoratedGraphByBlockAndCharacterTree inGraph = 
  if LG.isEmpty inGraph then V.empty
  else 
    let numBlocks = V.length $ vertData $ snd $ head $ LG.labNodes inGraph
        blockGraphList = fmap (pullBlock inGraph) [0.. (numBlocks - 1)] 
        characterGraphList = fmap makeCharacterGraph blockGraphList
    in
    -- trace ("Blocks " ++ (show numBlocks) ++ " Characters " ++ (show $ fmap length $ vertData $ snd $ head $ LG.labNodes inGraph)) 
    V.fromList characterGraphList

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
-- and cretes a new list of a singleton block from the input block index
makeBlockNodeLabels :: Int -> VertexInfo -> VertexInfo
makeBlockNodeLabels blockIndex inVertexInfo =
  let newVertexData = (vertData inVertexInfo) V.! blockIndex
      newVertexCost = V.sum $ fmap localCost newVertexData
      newsubGraphCost = V.sum $ fmap globalCost newVertexData
  in
  -- trace ("MBD " ++ (show $ length newVertexData) ++ " from " ++ (show $ length (vertData inVertexInfo)))
  inVertexInfo { vertData     = V.singleton newVertexData
               , vertexCost   = newVertexCost
               , subGraphCost = newsubGraphCost
               }

-- | makeCharacterGraph takes a blockGraph and creates a vector of character graphs
-- each with a single block and single character 
-- updating costs
makeCharacterGraph :: DecoratedGraph -> V.Vector DecoratedGraph
makeCharacterGraph inBlockGraph =
  if LG.isEmpty inBlockGraph then V.empty
  else 
    let numCharacters =  V.length $ V.head $ vertData $ snd $ head $ LG.labNodes inBlockGraph
        characterGraphList = if (numCharacters > 0) then fmap (pullCharacter False inBlockGraph) [0.. (numCharacters - 1)]
                             -- missing data
                             else [pullCharacter True inBlockGraph 0]
    in
    if (V.length $ vertData $ snd $ head $ LG.labNodes inBlockGraph) /= 1 then error "Number of blocks /= 1 in makeCharacterGraph"
    else 
      -- trace ("Chars: " ++ show numCharacters)
      V.fromList characterGraphList

-- | pullCharacter takes a DecoratedGraph with a single block and
-- creates a new DecoratedGraph with a single character form the input index
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
makeCharacterLabels :: Bool -> Int -> VertexInfo -> VertexInfo
makeCharacterLabels isMissing characterIndex inVertexInfo =
  let newVertexData = (V.head $ vertData inVertexInfo) V.! characterIndex
      newVertexCost = localCost newVertexData
      newSubGraphCost = globalCost newVertexData
  in
  -- trace ("MCL " ++ (show $ V.length $ vertData inVertexInfo) ++ " " ++ (show $ fmap  V.length $ vertData inVertexInfo) )
  {-
  if V.null $ V.head $ vertData inVertexInfo then 
      inVertexInfo { vertData    = V.singleton $ V.singleton emptyCharacter
                  , vertexCost   = 0
                  , subGraphCost = 0
                  }
  else
  -} 
     inVertexInfo { vertData     = if (not isMissing) then V.singleton $ V.singleton newVertexData
                                   else V.singleton V.empty
                  , vertexCost   = newVertexCost
                  , subGraphCost = newSubGraphCost
                  }
