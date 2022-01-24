{- |
Module      :  Wagner.hs
Description :  Module with distance tree construction methods Distance Wagner Farris 1972
              -- but with added refinement based on 4-point metric
Copyright   :  (c) 2020 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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
Need to integerize costs for swapping very slow on Double values
  do due spurious precision
-}

module Search.DistanceWagner (doWagnerS, performRefinement) where

import           Control.Parallel.Strategies
import           Data.Maybe
import qualified Data.Number.Transfinite     as NT
import qualified Data.Vector                 as V
import           Debug.Trace
import           GeneralUtilities
import           ParallelUtilities as PU
import qualified SymMatrix                   as M
import           Types.DistanceTypes
import           Utilities.DistanceUtilities
--import qualified LocalSequence as LS
import qualified Data.Vector as LS

-- | getStartingPair returns starying pair for Wagner build
--  closts mnimal cost pair
--  furthest maximal cost pair
--  random chooses uniformly at random from leaf set
getStartingPair :: String -> M.Matrix Double -> Edge
getStartingPair choiceOpt distMatrix
  | choiceOpt == "closest" = getMatrixMinPair distMatrix (-1, -1, NT.infinity) 0 0
  | choiceOpt == "furthest" = getMatrixMaxPair distMatrix (-1 , -1, 0 :: Double) 0 0
  | choiceOpt == "random" = errorWithoutStackTrace "Initial pair option 'random' not yet implemented"
  | otherwise = errorWithoutStackTrace ("Initial pair option " ++ choiceOpt ++ " unrecognized.  Must be 'closest', 'furthest', or 'random'")

-- | getBestEdgeTree take list of edge tuples and return trhe one with best addition cost
getBestEdgeTree :: V.Vector (Double, Tree, M.Matrix Double) -> Double -> (Double, Tree, M.Matrix Double) -> (Double, Tree, M.Matrix Double)
getBestEdgeTree edgeTreeList curBestCost curBestResult =
  if V.null edgeTreeList then curBestResult
  else
    let (firstAddCost, _, _) = V.head edgeTreeList
    in
    if firstAddCost < curBestCost then getBestEdgeTree (V.tail edgeTreeList) firstAddCost (V.head edgeTreeList)
    else getBestEdgeTree (V.tail edgeTreeList) curBestCost curBestResult

-- | add leaf to an edge creating new tree with distances and add cost, also augmented distance matrix
-- but this for swap so returns entire new3-edge cost  so not Farris triangle it is sum of three diveded by 2
addToEdgeSwap :: M.Matrix Double -> Int -> Tree -> Int -> Edge -> (Double, Tree, M.Matrix Double)
addToEdgeSwap distMatrix leaf initialTree newLeafIndex inEdge =
  let (eVertex, uVertex, inWeight) = inEdge
      (initialVertexVect, initialEdgeVect) = initialTree
      addCost = ((distMatrix M.! (leaf, eVertex)) + (distMatrix M.! (leaf, uVertex)) - (distMatrix M.! (eVertex, uVertex))) / 2.0
      eVertLeafDist = (distMatrix M.! (leaf, eVertex)) - addCost
      uVertLeafDist = (distMatrix M.! (leaf, uVertex)) - addCost
      newVertexVect = V.snoc initialVertexVect leaf
      newEdges = V.fromList [(leaf,newLeafIndex, addCost),(eVertex, newLeafIndex, eVertLeafDist),(uVertex, newLeafIndex, uVertLeafDist)]
      cleanupEdges = V.filter (/= inEdge) initialEdgeVect
      newEdgeVect = cleanupEdges V.++ newEdges
      newTree = (newVertexVect, newEdgeVect)
      -- add new costs from added vertex to each reamaining leaf
      augmentedDistMatrix = getNewDistMatrix distMatrix addCost eVertLeafDist uVertLeafDist eVertex uVertex leaf
  in
  (addCost + eVertLeafDist + uVertLeafDist - inWeight , newTree, augmentedDistMatrix)


-- | add leaf to an edge creating new tree with distances and add cost, also augmented distance matrix
addToEdge :: M.Matrix Double -> Int -> Tree -> Int -> Edge -> (Double, Tree, M.Matrix Double)
addToEdge distMatrix leaf initialTree newLeafIndex inEdge =
  -- trace ("In addToEdge with " ++ (show (leaf, initialTree, newLeafIndex, (M.rows distMatrix), inEdge))) (
  let (eVertex, uVertex, _) = inEdge
      (initialVertexVect, initialEdgeVect) = initialTree
      addCost = ((distMatrix M.! (leaf, eVertex)) + (distMatrix M.! (leaf, uVertex)) - (distMatrix M.! (eVertex, uVertex))) / 2.0
      eVertLeafDist = (distMatrix M.! (leaf, eVertex)) - addCost
      uVertLeafDist = (distMatrix M.! (leaf, uVertex)) - addCost
      newVertexVect = V.snoc initialVertexVect leaf
      newEdges = V.fromList [(leaf,newLeafIndex, addCost),(eVertex, newLeafIndex, eVertLeafDist),(uVertex, newLeafIndex, uVertLeafDist)]
      cleanupEdges = V.filter (/= inEdge) initialEdgeVect
      newEdgeVect = V.map orderEdge $ cleanupEdges V.++ newEdges
      newTree = (newVertexVect, newEdgeVect)
      -- add new costs from added vertex to each reamaining leaf
      augmentedDistMatrix = getNewDistMatrix distMatrix addCost eVertLeafDist uVertLeafDist eVertex uVertex leaf
  in
  (addCost, newTree, augmentedDistMatrix)


-- | addTaxonToTree takes distMatrix, an initialTree, Vector of leavesToAdd, and leaf index to add
-- and retursn a tuple wiht the addition cost, the new tree, the new leaves to add, and new distance matrix (enhanced)
addTaxonToTree :: M.Matrix Double -> Tree -> V.Vector Int -> Int -> Int -> (Double, Tree, V.Vector Int, M.Matrix Double)
addTaxonToTree distMatrix initialTree leavesToAdd newVertexIndex leaf =
  if V.null leavesToAdd then (0.0, initialTree, leavesToAdd, distMatrix)
  else
    let leavesRemaining = V.filter (/= leaf) leavesToAdd
        (_, edgesInitial) = initialTree

        -- Parallelize heretoo much and destroys lazy matrix update
        addEdgeList = V.map (addToEdge distMatrix leaf initialTree newVertexIndex) edgesInitial
        (firstAddCost, _, _) = V.head addEdgeList -- this to initialize getBestEdge below
    in
    --filter for best addition point
    let (addCost, newTree, augmentedDistMatrix) = getBestEdgeTree (V.tail addEdgeList) firstAddCost (V.head addEdgeList)
    in
    (addCost, newTree, leavesRemaining, augmentedDistMatrix)

-- | getBestLeafAdd chooses best leaf to add based on cost field
getBestLeafAdd ::  V.Vector (Double, Tree, V.Vector Int, M.Matrix Double) -> Double -> (Double, Tree, V.Vector Int, M.Matrix Double) -> (Double, Tree, V.Vector Int, M.Matrix Double)
getBestLeafAdd addPosVect curBestCost curBestLeaf =
  if V.null addPosVect then curBestLeaf
  else
    let (thisCost, _, _, _) = V.head addPosVect
    in
    if thisCost < curBestCost then getBestLeafAdd (V.tail addPosVect) thisCost (V.head addPosVect)
    else getBestLeafAdd (V.tail addPosVect) curBestCost curBestLeaf


-- | wagBest takes distMatrix, and intial tree of two leaves, a vector of leavesToAdd, the input nuber of leaves
-- and returns the Farris 1972 distance Wagner adding the "closest" leaf at each iteration
wagBest :: M.Matrix Double -> Tree -> V.Vector Int -> Int -> Int ->  V.Vector Int -> String -> (Tree, V.Vector Int, M.Matrix Double)
wagBest distMatrix inTree leavesToAdd nOTUs newVertexIndex leavesToMap choiceOpt
  | length leavesToAdd == nOTUs =
  let (eVertex, uVertex, edgeWeight) =  orderEdge $ getStartingPair choiceOpt distMatrix
      initialTree = (V.fromList[eVertex, uVertex],V.fromList [(eVertex, uVertex, edgeWeight)])
      leavesToAdd' = V.filter (/= eVertex) $ V.filter (/= uVertex) leavesToAdd
  in
  wagBest distMatrix initialTree leavesToAdd' nOTUs nOTUs leavesToAdd' choiceOpt
  | V.null leavesToAdd = (inTree, leavesToAdd, distMatrix)
  | otherwise =
  let addPosVect = V.map (addTaxonToTree distMatrix inTree leavesToAdd newVertexIndex) leavesToMap
      (firstLeafCost, _, _, _ ) = V.head addPosVect -- To initialize below
      (_, newTree, newLeavesToAdd, augmentedDistMatrix) = getBestLeafAdd (V.tail addPosVect) firstLeafCost (V.head addPosVect)
  in
  let progress = takeWhile (/='.') $ show  ((fromIntegral (100 * (newVertexIndex - nOTUs))/fromIntegral (nOTUs - 2)) :: Double)
  in
  if last progress == '0' then
      trace ("\t\t"++ progress ++ "%")
      wagBest augmentedDistMatrix newTree newLeavesToAdd nOTUs (newVertexIndex + 1) newLeavesToAdd choiceOpt
  else wagBest augmentedDistMatrix newTree newLeavesToAdd nOTUs (newVertexIndex + 1) newLeavesToAdd choiceOpt


-- | calculateWagnerTrees takes an input distance matrix (and options later) and returns
-- a tree (V,E) discription of Wagner tree with labelled internal veritices and branch lengths
calculateWagnerTrees :: M.Matrix Double -> String -> (Tree, M.Matrix Double)
calculateWagnerTrees distMatrix choiceOpt =
  if M.dim distMatrix == (0,0) then errorWithoutStackTrace "Null distance matrix"
  else
    -- get initial pair of leaves and create initial tree
    let nOTUs = M.cols distMatrix
        allLeaves = V.fromList [0..(nOTUs - 1)]
    in
    let (newTree, _, endMatrix) = wagBest distMatrix (V.empty, V.empty) allLeaves nOTUs nOTUs allLeaves choiceOpt
    in
    (newTree, endMatrix)

-- | makeTreeFromOrder  takes an input order and other arguemnts and cretes tree using a single additoin
-- seqeunce, best plaecment for the leaf each round
makeTreeFromOrder :: M.Matrix Double -> Tree -> Int -> Int -> V.Vector Int -> (Tree, M.Matrix Double)
makeTreeFromOrder distMatrix initialTree nOTUs vertexIndex leavesToAdd =
  if null leavesToAdd then (initialTree, distMatrix)
  else
      let leaf = V.head leavesToAdd
          (_, newTree, _, augmentedDistMatrix) = addTaxonToTree distMatrix initialTree leavesToAdd vertexIndex leaf
      in
      makeTreeFromOrder augmentedDistMatrix newTree nOTUs (vertexIndex + 1) (V.tail leavesToAdd)

-- | getRandomAdditionSequence initializes based on input sequence and adds in order from there
getRandomAdditionSequence :: V.Vector String -> M.Matrix Double -> Int -> V.Vector Int -> TreeWithData
getRandomAdditionSequence leafNames distMatrix outgroup initiaLeavesToAdd =
  let nOTUs = V.length leafNames
  in
  let eVertex = initiaLeavesToAdd V.! 0
      uVertex = initiaLeavesToAdd V.! 1
      edgeWeight = distMatrix M.! (eVertex, uVertex)
      initialTree = (V.fromList[eVertex, uVertex],V.fromList [(eVertex, uVertex, edgeWeight)])
      leavesToAdd = V.filter (/= eVertex) $ V.filter (/= uVertex) initiaLeavesToAdd
  in
  let thisTree = makeTreeFromOrder distMatrix initialTree nOTUs nOTUs leavesToAdd
      -- (_, edgeVect) = fst thisTree
      treeCost = getTreeCost $ fst thisTree -- V.sum $ V.map getEdgeCost edgeVect
      newickTree = convertToNewick leafNames outgroup (fst thisTree)
      newickTree' = take (length newickTree - 3) newickTree ++ "[" ++ showDouble precision treeCost ++ "]" ++ ";"
  in
  (newickTree', fst thisTree, treeCost, snd thisTree)

-- | doWagnerS takes user options and produces the Wagner tree methods desired (best, asis, or random)
-- outputs newick rep list
doWagnerS :: V.Vector String -> M.Matrix Double -> String -> Int -> String -> [V.Vector Int]-> [TreeWithData]
doWagnerS leafNames distMatrix firstPairMethod outgroup addSequence replicateSequences =
  let nOTUs = V.length leafNames
  in
  if addSequence == "best" then
     let wagnerResult = calculateWagnerTrees distMatrix firstPairMethod
         -- (_, edgeVect) = fst wagnerResult
         treeCost = getTreeCost $ fst wagnerResult --- V.sum $ V.map getEdgeCost edgeVect
         newickTree = convertToNewick leafNames outgroup (fst wagnerResult)
         newickTree' = take (length newickTree - 3) newickTree ++ "[" ++ showDouble precision treeCost ++ "]" ++ ";"
      in
      [(newickTree', fst wagnerResult, treeCost, snd wagnerResult)]
  else if addSequence == "asis" then
      let initialTree = (V.fromList[0, 1],V.fromList [(0, 1, distMatrix M.! (0,1))])
          leavesToAdd = V.fromList [2..(nOTUs-1)]
          asIsResult = makeTreeFromOrder distMatrix initialTree nOTUs nOTUs leavesToAdd
          treeCost = getTreeCost $ fst asIsResult -- V.sum $ V.map getEdgeCost asIsEdges
          newickTree = convertToNewick leafNames outgroup (fst asIsResult)
          newickTree' = take (length newickTree - 3) newickTree ++ "[" ++ showDouble precision treeCost ++ "]" ++ ";"
      in
      [(newickTree', fst asIsResult, treeCost, snd asIsResult)]
  else if head addSequence == 'r' then
      if null replicateSequences then errorWithoutStackTrace "Zero replicate additions specified--could be error in configuration file"
      else
        let randomAddTrees = fmap (getRandomAdditionSequence leafNames distMatrix outgroup) replicateSequences `using` myParListChunkRDS -- was rseq not sure whats better
            -- randomAddTrees = parmap rseq (getRandomAdditionSequence leafNames distMatrix outgroup) replicateSequences
        in
        randomAddTrees
  else errorWithoutStackTrace ("Addition sequence " ++ addSequence ++ " not implemented")

-- | edgeHasVertex takes an vertex and an edge and returns Maybe Int
-- of other vertex
edgeHasVertex :: Vertex -> Edge -> Maybe (Vertex, Edge)
edgeHasVertex inVert inEdge =
  let (a, b, _) = inEdge
  in
  if a == inVert then Just (b, inEdge)
  else if b == inVert then Just (a, inEdge)
  else Nothing

-- | getSubEdges take a Vector of edges and retuns a list of edges connected to input vertex
-- uses nOTUs to know when to stop recursing
getSubEdges :: [Vertex] -> Int -> V.Vector Edge -> V.Vector Edge -> String -> V.Vector Edge
getSubEdges inVertexList nOTUs edgeVect subEdgeVect howMany
  | V.null edgeVect = subEdgeVect
  | null inVertexList = subEdgeVect
  | otherwise =
    let startVertex = head inVertexList
        foundVect = V.filter (/= Nothing) $ V.map (edgeHasVertex startVertex) edgeVect
    in
    if V.null foundVect then getSubEdges (tail inVertexList) nOTUs edgeVect subEdgeVect howMany-- only terminals left
    else if V.length foundVect /= 2 then error ("Index (" ++ howMany ++ ")" ++ show startVertex ++ "->found " ++ show (V.length foundVect) ++ " but should be two edges in " ++ show foundVect ++ " in " ++ show edgeVect)
    else
      let thingsFound = V.map fromJust foundVect
          verticesFound = V.toList $ V.map fst thingsFound
          edgesFound = V.map snd thingsFound

          -- add found edges to subEdgeSet
          newSubEdgeVect =  subEdgeVect V.++ edgesFound
          -- delete delete edges from those to be searched
          newEdgeVect = V.filter (/= V.last edgesFound) $ V.filter (/= V.head edgesFound) edgeVect
      in
      -- recurse on vertices that were found
      if howMany == "first2" then edgesFound
      else getSubEdges (verticesFound ++ tail inVertexList) nOTUs newEdgeVect newSubEdgeVect howMany

-- | adjustInternalEdgeVertex adjusts vertex of an internal (ie non-pendent, none 2 terminal edge)
-- assumes hv > lv
adjustInternalEdgeVertex :: Vertex -> Vertex -> Vertex -> Int -> Int -> Vertex
adjustInternalEdgeVertex inV hV lV maxOffSet nOTUs
  | inV <= max nOTUs lV = inV
  | inV > hV = inV - maxOffSet
  | inV > lV = inV - 1
  | otherwise = error ("This can't happen " ++ show (inV, lV, hV))
  -- )

-- | adjustVertex reduces vertex (hV,lV,_) is the dge that was deleted
-- index by 2,1, or 0 if greater than hVert, lVert, or neither
-- asumes hv > lv; and inE > inU
-- if selfe edge  its a termin and only reduce by 1
-- assumes edge is ordered e > u
adjustVertex :: Edge -> Vertex -> Vertex -> Int -> Edge
adjustVertex (inE, inU, w) hV lV nOTUs
  | (inE <= max lV nOTUs) && (inU <= max lV nOTUs) = (inE, inU, w)
  | lV < nOTUs =  -- update pendant edge was deleted since hV > lV; hV must be > nOTUs
    if (inE > max hV nOTUs) && (inU > max hV nOTUs) then (inE - 1, inU - 1, w)
    else if (inE <= max hV nOTUs) && (inU > max hV nOTUs) then (inE, inU - 1, w)
    else if (inE > max hV nOTUs) && (inU <= max hV nOTUs) then (inE - 1, inU, w)
    else (inE, inU, w)
  | otherwise =                                                          -- internal edge was deleted                                                         -- both deleted verteces internal
    let newE = adjustInternalEdgeVertex inE hV lV 2 nOTUs
        newU = adjustInternalEdgeVertex inU hV lV 2 nOTUs
    in
    (newE, newU, w)

-- | updateVertexNUmbersOnEdges taked vertex numbers and updates HTU indices to reflect
-- the deletion of those vertices
-- ASSUMES edges are ordered (a,b,weight) a > b
-- subtracts 2 from index if > the bigger of teh two, subtract 1 if bigger than lower,
-- otherwise leaves unchanged
updateVertexNUmbersOnEdges  :: Vertex -> Vertex -> V.Vector Edge -> Int -> V.Vector Edge
updateVertexNUmbersOnEdges eVert uVert edgeList nOTUs =
  -- trace ("In updateVertexNUmbersOnEdges") (
  if V.null edgeList then V.empty
  else
      let hVert = max eVert uVert
          lVert = min eVert uVert
          newVertex = adjustVertex (V.head edgeList) hVert lVert nOTUs
      in
      V.cons newVertex (updateVertexNUmbersOnEdges eVert uVert (V.tail edgeList) nOTUs)


-- | updateDistMatrix updates distance matrix to remove eVertex and uVertex columns and rows as in updateVertexNUmbersOnEdges
-- update costs for two contracte edges c1Edge and c2Edge
-- the distance update takes place first and row/coumn removal second
-- this to keep the proper indices (since the non-deleted rows and columns indices are changed)
updateDistMatrix :: Vertex -> Vertex ->  M.Matrix Double -> Int -> Edge -> Edge -> M.Matrix Double
updateDistMatrix eVert uVert distMatrix nOTUs c1Edge c2Edge =
  let validEdgeList = filter ((>= 0).fst3) [c1Edge, c2Edge]
      newMatrix = M.unsafeUpdateMatrix distMatrix validEdgeList
      newMatrix' = M.deleteRowsAndColumns newMatrix (filter (> (nOTUs - 1)) [eVert, uVert])
  in
  newMatrix'



-- | getEndVertices takes a pair of edges and returns the non-index vertices
getEndVertices :: V.Vector Edge -> Vertex -> (Vertex, Vertex)
getEndVertices inEdges index =
  if V.length inEdges /= 2 then error ("Edge number should be 2 not " ++ show (V.length inEdges) ++ " in getEndVertices")
  else
    let (fEVertex, fUVertex, _) = V.head inEdges
        (sEVertex, sUVertex, _) = V.last inEdges
    in
    if fEVertex == index then
      if sEVertex == index then (fUVertex, sUVertex)
      else
          (fUVertex, sEVertex)
    else if fUVertex == index then
      if sEVertex == index then
         (fEVertex, sUVertex)
      else
          (fEVertex, sEVertex)
    else error ("Error finding ends of " ++ show (V.head inEdges) ++ " and " ++ show (V.last inEdges))

-- | getContractedWeight gets edge contractd edge weights ither form distMatrix for  OTUs or 1 OTU and 1
-- internal wertex or by 4 point metric (maximum of two estimations)
getContractedWeight :: Vertex -> Vertex -> M.Matrix Double -> Int -> V.Vector Edge -> Double
getContractedWeight aVert bVert distMatrix nOTUs edgeVect
  | aVert < nOTUs && bVert < nOTUs = distMatrix M.! (aVert, bVert)
  | aVert < nOTUs = distMatrix M.! (aVert, bVert)
  | bVert < nOTUs = distMatrix M.! (aVert, bVert)
  | otherwise =
  -- both internal 4-point mertic estimation
  -- get edges connected to the contracted edge
  let aEdgeVect = getSubEdges [aVert] nOTUs edgeVect V.empty "first2"
      bEdgeVect = getSubEdges [bVert] nOTUs edgeVect V.empty "first2"
      (a,b) = getEndVertices aEdgeVect aVert
      (c,d) = getEndVertices bEdgeVect bVert
      firstEstimate  = (distMatrix M.! (a,c)) + (distMatrix M.! (b,d)) - (distMatrix M.! (a,b)) - (distMatrix M.! (c,d))
      secondEstimate = (distMatrix M.! (a,b)) + (distMatrix M.! (b,c)) - (distMatrix M.! (a,b)) - (distMatrix M.! (c,d))
  in
  max firstEstimate secondEstimate / 2.0


-- | contractEdges takes two edges that share a vertex and fuses them
-- If terminal then return selfdge
-- update contracted edge weight
-- dummy edge (-1,-1,-1) is returned when no edge to contract so won't be updated in distance Matrix
contractEdges :: M.Matrix Double -> Int -> V.Vector Edge -> Vertex -> V.Vector Edge -> Edge
contractEdges distMatrix nOTUs edgeVect index allEdges
  | V.null edgeVect = (-1,-1,0)
  | V.length edgeVect == 1 =
    let (a,b,w) = V.head edgeVect
    in
    if (a == index) || (b == index) then (index, index, w)
    else error ("Contacting single edge: " ++ show (V.head edgeVect) ++ " with vertex " ++ show index ++ " not found")
  | otherwise =
    let (aVert,bVert) = getEndVertices edgeVect index
        newWeight = getContractedWeight aVert bVert distMatrix nOTUs (subtractVector edgeVect allEdges)
    in
    orderEdge (aVert, bVert, newWeight)

-- | getDistance creates distances from remaing OTUs to new vertex (rowID) via Farris 1972
-- rowID is the row (or column) of distance Matrix
getDistance :: M.Matrix Double -> Double -> Double -> Double -> Int -> Int -> Int -> Int -> Double
getDistance origDist addCost eVertLeafDist uVertLeafDist leafIndex eVertex uVertex rowID =
  let first  = (origDist M.! (rowID, leafIndex)) - addCost
      second = (origDist M.! (rowID, eVertex)) - eVertLeafDist
      third  = (origDist M.! (rowID, uVertex)) - uVertLeafDist
  in
  maximum [first, second, third]

-- | getNewDistMatrix takes distMatrix and adds Cost for new eVertLeafDist uVertLeafDist
-- created in build process (new HTUs)
-- should be complete (on input) for leaves already added (initail paiwise distances and HTUs added)
-- adds a single new row (minus last 0.0 as new row at end) which is appended
getNewDistMatrix :: M.Matrix Double -> Double -> Double -> Double -> Int -> Int -> Int -> M.Matrix Double
getNewDistMatrix origDist addCost eVertLeafDist uVertLeafDist eVertex uVertex leafIndex =
    let columnHolder = LS.fromList [0..(M.rows origDist - 1)] -- List of HTU and OTU indices in pairwise dist matrix
        newDistRow = LS.map (getDistance origDist addCost eVertLeafDist uVertLeafDist leafIndex eVertex uVertex) columnHolder
        newDistRow' = newDistRow `LS.snoc` (0.0 :: Double)
    in
    M.addMatrixRow origDist newDistRow'

-- | enterNewEdgeCost cretes new row of costs from edge vectors
-- infty NT.infinity if not there
-- this makes update n^2 which is dumb
enterNewEdgeCost :: Int -> V.Vector Edge -> Int -> Double
enterNewEdgeCost columnNumber edgeVect rowNumber =
  if V.null edgeVect then NT.infinity --not found so Infty
  else
    let (a, b, weight) = V.head edgeVect
    in
    if (columnNumber == a && rowNumber == b) || (columnNumber == b && rowNumber == a) then weight else enterNewEdgeCost columnNumber (V.tail edgeVect) rowNumber

-- | addEdgesToDistMatrix adds new edges from internal node join to existing matrix
-- 0 if not created so only the two new vertices and the edges they touch (3 each)
-- so need to add two columns and two rows, for new vertices I, and II (here numIn and numIn + 1)
getNewDistMatrixInternal :: M.Matrix Double -> V.Vector Edge -> M.Matrix Double
getNewDistMatrixInternal inMatrix newEdgeVect =
  -- trace ("In getNewDistMatrixInternal") (
  if V.length newEdgeVect /= 5 then error ("Wrong size edgeVector shoud be 5 and is " ++ show (V.length newEdgeVect))
  else
    let numIn = M.rows inMatrix
        columnHolder = [0..(numIn - 1)]
        newDistColumnI = fmap (enterNewEdgeCost numIn newEdgeVect) columnHolder ++ [0.0]
        newDistColumnII = fmap (enterNewEdgeCost (numIn + 1) newEdgeVect) columnHolder ++ [enterNewEdgeCost numIn newEdgeVect (numIn + 1), 0.0]
    in
    M.addMatrices inMatrix  (LS.fromList [LS.fromList newDistColumnI, LS.fromList newDistColumnII])

-- | connectEdges takes two vectors of edges and adds and edge between the two in edgesToConnect
-- this deletes the two old edges from the edge Vectors and creates a new tree with the five new edges
-- added in; the addition cost is for total of created edges minus edges that were destroyed
connectEdges :: M.Matrix Double -> V.Vector Edge -> V.Vector Edge -> V.Vector Edge -> (Double, Tree, M.Matrix Double)
connectEdges distMatrix eEdges uEdges edgesToConnect
  | V.null eEdges = error "Empty e Edge vector in connectEdges"
  | V.null uEdges = error "Empty u Edge vector in connectEdges"
  | V.length edgesToConnect /= 2 = error ("There shoule be 2 edges to connect and ther are " ++ show (V.length edgesToConnect))
  | otherwise =
  let edgesKept = subtractVector edgesToConnect eEdges V.++ subtractVector edgesToConnect uEdges
      numInTree = M.rows distMatrix
      -- order new edges
      (a, b, wAB) = V.head edgesToConnect
      (c, d, wCD) = V.last edgesToConnect
      firstEstimate  = (distMatrix M.! (a,c)) + (distMatrix M.! (b,d)) - (distMatrix M.! (a,b)) - (distMatrix M.! (c,d))
      secondEstimate = (distMatrix M.! (a,b)) + (distMatrix M.! (b,c)) - (distMatrix M.! (a,b)) - (distMatrix M.! (c,d))
      centralEdgeCost  = max firstEstimate secondEstimate / 2.0
      centralEdge = orderEdge (numInTree, numInTree + 1, centralEdgeCost)
      aCost = ((distMatrix M.! (c,a)) - (distMatrix M.! (c,b)) + (distMatrix M.! (a,b))) / 2.0
      bCost = (distMatrix M.! (a,b)) - aCost
      cCost = ((distMatrix M.! (a,c)) - (distMatrix M.! (a,d)) + (distMatrix M.! (c,d))) / 2.0
      dCost = (distMatrix M.! (c,d)) - cCost
      newAEdge = orderEdge (a, numInTree, aCost) -- edge cost not unique either (a or b, numInTree) could be max or min
      newBEdge = orderEdge (b, numInTree, bCost) -- edge cost not unique either (a or b, numInTree) could be max or min
      newCEdge = orderEdge (c, numInTree + 1, cCost) -- edge cost not unique either (a or b, numInTree + 1) could be max or min
      newDEdge = orderEdge (d, numInTree + 1, dCost) -- edge cost not unique either (a or b, numInTree + 1) could be max or min
      newEdgeVect = V.fromList [centralEdge, newAEdge, newBEdge, newCEdge, newDEdge]
      newDistMatrix = getNewDistMatrixInternal distMatrix newEdgeVect
  in
  (centralEdgeCost + aCost + bCost + cCost + dCost - wAB - wCD, (V.empty, newEdgeVect V.++ edgesKept), newDistMatrix)


-- | addEdgeToSplit adds a new edge to a specifc pair of input edges detemined addition cost
-- and creted new edge set with appropriate weights
addEdgeToSplit :: V.Vector Edge -> V.Vector Edge -> Edge -> Edge -> M.Matrix Double -> V.Vector Edge -> (Double, Tree, M.Matrix Double)
addEdgeToSplit eEdges uEdges eTerminal uTerminal distMatrix edgesToConnect
  | V.null eEdges &&  V.null uEdges = error "Empty e/u Edge vectors in addEdgeToSplit"
  | V.null edgesToConnect = error "Empty eEdge vector in addEdgeToSplit"
  | fst3 (V.head eEdges) == (-1) = -- adding eOTU, V.last since the (-1,-1,0) edge should always be first
      addToEdgeSwap distMatrix (fst3 eTerminal) (V.empty, uEdges) (M.rows distMatrix) (V.last edgesToConnect)
  | fst3 (V.head uEdges) == (-1) = -- adding uOTU
      addToEdgeSwap distMatrix (fst3 uTerminal) (V.empty, eEdges) (M.rows distMatrix) (V.last edgesToConnect)
  | otherwise = -- both internal edges
      connectEdges distMatrix eEdges uEdges edgesToConnect

-- | splitTree takes a tree description and its edgeList and return pairs of edge list
-- split at input edge in tree with "repaird"/contracted edges, delta, and original
-- edge (pairs of vertices and weight) for each split
splitTree :: M.Matrix Double -> Tree -> Double -> Edge -> SplitTreeData
splitTree distMatrix inTree inTreeCost edgeToRemove =
  -- check if proper tree--remove later
  let (_, edgeVect) = inTree
      (eVertex, uVertex, _) = edgeToRemove

      -- newEdgeSet = subtractVector (V.cons edgeToRemove $ eEdges V.++ uEdges) edgeVect
      newEdgeSet = V.filter (/= edgeToRemove) edgeVect
      nOTUs = div (3 + V.length edgeVect) 2
      eSubEdges = getSubEdges [eVertex] nOTUs newEdgeSet V.empty "all"
      uSubEdges = getSubEdges [uVertex] nOTUs newEdgeSet V.empty "all"

      -- get edges that need to be contracted and re-estimate weights
      eEdges = getSubEdges [eVertex] nOTUs eSubEdges V.empty "first2"
      uEdges = getSubEdges [uVertex] nOTUs uSubEdges V.empty "first2"
      eMergedEdge = contractEdges distMatrix nOTUs eEdges eVertex edgeVect
      uMergedEdge = contractEdges distMatrix nOTUs uEdges uVertex edgeVect

      -- remove non-contracted edges and add in contracted edges
      eSubEdges' = V.cons eMergedEdge (subtractVector eEdges eSubEdges)
      uSubEdges' = V.cons uMergedEdge (subtractVector uEdges uSubEdges)

      -- need to know this order for SPR/TBR so which edge was in which set
      previousEdges = V.fromList [eMergedEdge, uMergedEdge]

      -- map new HTU indices in edges and create new distance matrix
      -- HTU indices are updated first--then values updated in distMatrix as rows/columns are
      -- deleted to remove info from delted edge
      eSubEdges'' = updateVertexNUmbersOnEdges eVertex uVertex (V.map orderEdge eSubEdges') nOTUs
      uSubEdges'' = updateVertexNUmbersOnEdges eVertex uVertex (V.map orderEdge uSubEdges') nOTUs
      previousEdges'' = updateVertexNUmbersOnEdges eVertex uVertex (V.map orderEdge previousEdges) nOTUs

      -- Update with deleted node and reestimated contracted edges
      -- update matrix the costs (2 ij x 2 ji) of contracted edges
      distMatrix'' = updateDistMatrix eVertex uVertex distMatrix nOTUs eMergedEdge uMergedEdge -- (V.head previousEdges'') (V.last previousEdges'')

      -- Delta calcualted by differene in original tree cost and split and readdition to split edges (then tail splits later
      -- so not remake the original tree)
      splitCost = getTreeCost (V.empty, eSubEdges'') + getTreeCost (V.empty, uSubEdges'')

      -- Delta of tree length will be weights of the the edges removed - the weights of the two edges contracted and reestimated
      -- delta = weight + (V.sum $ V.map thd3 eEdges) +  (V.sum $ V.map thd3 uEdges) - (V.sum $ V.map thd3 previousEdges'')
      delta =  inTreeCost - splitCost   -- readditionCost
  in
  if eVertex < nOTUs then (V.singleton (eVertex, eVertex, 0.0), uSubEdges'', delta, previousEdges'', distMatrix'') --this so know a pendant edge
  else if uVertex < nOTUs then (V.singleton (uVertex, uVertex, 0.0), eSubEdges'', delta, previousEdges'', distMatrix'') --this so know a pendant edge
  else
      -- neded to have e first then u for SPR/TBR so can coordinate with previous edges
      (eSubEdges'', uSubEdges'', delta, previousEdges'', distMatrix'')


-- | sieveTrees takes a list of (addition cost, Tree, distnce matrix) and returns list
-- of better or equal trees to input delta
sieveTrees :: Double -> Double -> V.Vector (Double, Tree, M.Matrix Double) -> V.Vector String -> Int -> [TreeWithData] -> [TreeWithData]
sieveTrees inDelta curBestCost inAddList leafNames outgroup savedTrees =
  if null inAddList then savedTrees
  else
      let firstTuple = V.head inAddList
          (firstDelta, firstTree, firstMatrix) = firstTuple
          newCost = curBestCost - inDelta + firstDelta
          -- checkCost = getTreeCost firstTree
          newickTree = convertToNewick leafNames outgroup firstTree
          newickTree' = take (length newickTree - 3) newickTree ++ "[" ++ showDouble precision newCost ++ "]" ++ ";"
          newTuple = (newickTree', firstTree, newCost, firstMatrix)
      in
      if firstDelta > inDelta then sieveTrees inDelta curBestCost  (V.tail inAddList) leafNames outgroup savedTrees
      else
        if newCost < curBestCost then
            sieveTrees inDelta curBestCost (V.tail inAddList) leafNames outgroup [newTuple]
        else if withinEpsilon newCost curBestCost then sieveTrees inDelta curBestCost (V.tail inAddList) leafNames outgroup (newTuple : savedTrees)
        else sieveTrees inDelta curBestCost (V.tail inAddList) leafNames outgroup savedTrees


-- | reAddTerminals checks to see if split on pendant edge--if so reads terminal to each edge but saves equal cost if found
-- and saveMethod specifies it.  Idenitical to the wagner addition process with saveing equal as option
-- could use addEdgeToSplit for consistancey with SPR/TBR
reAddTerminals :: String -> Double -> V.Vector String -> Int -> SplitTreeData -> [TreeWithData]
reAddTerminals rejoinType curBestCost leafNames outGroup split =
  if rejoinType /= "otu" then error ("Incorrect swap function in reAddTerminals: " ++ rejoinType)
  else
    let (eEdgeVect, uEdgeVect, delta, _, distMatrix) = split
        nOTUs = V.length leafNames
    in
    if (V.length eEdgeVect > 1) || ((fst3 (V.head eEdgeVect) /= snd3 (V.head eEdgeVect)) && (fst3 (V.head eEdgeVect) < nOTUs) && (snd3 (V.head eEdgeVect) < nOTUs)) then []
    else (if M.rows distMatrix /= ((2 * nOTUs) - 3) then error ("Dist Matrix incorrect size " ++ show (M.dim distMatrix) ++ " should be " ++ show ((2 * nOTUs) - 3, (2 * nOTUs) - 3))
    else
      let newLeafIndex = M.rows distMatrix
      -- take tail of uEdgeVect so not regerate input tree
          additionList = V.map (addToEdgeSwap distMatrix (fst3 $ V.head eEdgeVect) (V.empty,uEdgeVect) newLeafIndex) uEdgeVect -- (V.tail uEdgeVect) -- tail so not hit original tree, leave all to reestimate if necesary
          minAdditionCost = V.minimum (V.map fst3 additionList)
      in
      if minAdditionCost > delta then []
      else sieveTrees delta curBestCost additionList leafNames outGroup [])


-- | add leaf to an edge creating new tree with distances and add cost, also augmented distance matrix
-- but this for swap so returns entire new3-edge cost  so not Farris triangle it is sum of three diveded by 2
addToEdgeSwapRecurse :: Double -> M.Matrix Double -> Int -> Tree -> Int -> V.Vector Edge -> (Double, Tree, M.Matrix Double)
addToEdgeSwapRecurse inDelta distMatrix leaf initialTree newLeafIndex inEdgeVect =
  if V.null inEdgeVect then (inDelta, initialTree, distMatrix)
  else
    let inEdge@(eVertex, uVertex, inWeight) = V.head inEdgeVect
        (initialVertexVect, initialEdgeVect) = initialTree
        addCost = ((distMatrix M.! (leaf, eVertex)) + (distMatrix M.! (leaf, uVertex)) - (distMatrix M.! (eVertex, uVertex))) / 2.0
        eVertLeafDist = (distMatrix M.! (leaf, eVertex)) - addCost
        uVertLeafDist = (distMatrix M.! (leaf, uVertex)) - addCost
        newVertexVect = V.snoc initialVertexVect leaf
        newEdges = V.fromList [(leaf,newLeafIndex, addCost),(eVertex, newLeafIndex, eVertLeafDist),(uVertex, newLeafIndex, uVertLeafDist)]
        cleanupEdges = V.filter (/= inEdge) initialEdgeVect
        newEdgeVect = cleanupEdges V.++ newEdges
        newTree = (newVertexVect, newEdgeVect)
        -- add new costs from added vertex to each reamaining leaf
        augmentedDistMatrix = getNewDistMatrix distMatrix addCost eVertLeafDist uVertLeafDist eVertex uVertex leaf
        newDelta = addCost + eVertLeafDist + uVertLeafDist - inWeight
    in
    if newDelta < inDelta then (newDelta, newTree, augmentedDistMatrix)
    else addToEdgeSwapRecurse inDelta distMatrix leaf initialTree newLeafIndex (V.tail inEdgeVect)

-- | getVectorAllVectorPairs takes two vectors and creates a vector of avector of two elements each for each
-- pairwise combinatrion of elements
getVectorAllVectorPairs :: V.Vector a -> V.Vector a -> V.Vector (V.Vector a)
getVectorAllVectorPairs firstVect secondVect =
  if V.null firstVect then V.empty
  else
    let firstElement = V.head firstVect
        firstPairs = V.map (V.cons firstElement) $ V.map V.singleton secondVect
    in
    firstPairs V.++ getVectorAllVectorPairs (V.tail firstVect) secondVect

-- | createVectorEdgePairs creates teh Vector of Vectors of edges (2 in each case) to connect
-- if SPR then takes the initial (previous e edge) and pairs with all in u edge vect
-- if TBR then all conbinations of pairs
createVectorEdgePairs :: String -> V.Vector Edge -> V.Vector Edge -> V.Vector Edge -> V.Vector (V.Vector Edge)
createVectorEdgePairs pairSet previousEdges eEdgeVect uEdgeVect
  | pairSet == "spr" =
    let eEdgePrev = V.head previousEdges
    in
    V.map (V.cons eEdgePrev) $ V.map V.singleton uEdgeVect
  | pairSet == "tbr" = getVectorAllVectorPairs eEdgeVect uEdgeVect
  | otherwise = errorWithoutStackTrace ("Pair set option " ++ pairSet ++ " not implemented")


-- | addEdgeToSplitRecurse like addToEdgeSplit but recursiblye yeilds a single best tree
addEdgeToSplitRecurse :: V.Vector Edge -> V.Vector Edge -> Edge -> Edge -> M.Matrix Double -> V.Vector (V.Vector Edge) -> (Double, Tree, M.Matrix Double) -> (Double, Tree, M.Matrix Double)
addEdgeToSplitRecurse eEdges uEdges eTerminal uTerminal distMatrix edgesToConnectVect origTriple@(inDelta, _, _) =
  if V.null edgesToConnectVect then origTriple
  else
    let edgesToConnect = V.head edgesToConnectVect
    in
    if V.null eEdges &&  V.null uEdges then error "Empty e/u Edge vectors in addEdgeToSplit"
    else if V.null edgesToConnect then error "Empty eEdge vector in addEdgeToSplit"
    else if fst3 (V.head eEdges) == (-1) then -- adding eOTU, V.last since the (-1,-1,0) edge should always be first
      let (newDelta, newTree, newMatrix) = addToEdgeSwap distMatrix (fst3 eTerminal) (V.empty, uEdges) (M.rows distMatrix) (V.last edgesToConnect)
      in
      if newDelta < inDelta then (newDelta, newTree, newMatrix)
      else addEdgeToSplitRecurse eEdges uEdges eTerminal uTerminal distMatrix (V.tail edgesToConnectVect) origTriple
    else if fst3 (V.head uEdges) == (-1) then -- adding uOTU
      let (newDelta, newTree, newMatrix) = addToEdgeSwap distMatrix (fst3 uTerminal) (V.empty, eEdges) (M.rows distMatrix) (V.last edgesToConnect)
      in
      if newDelta < inDelta then (newDelta, newTree, newMatrix)
      else addEdgeToSplitRecurse eEdges uEdges eTerminal uTerminal distMatrix (V.tail edgesToConnectVect) origTriple
    else -- both internal edges
      let (newDelta, newTree, newMatrix) = connectEdges distMatrix eEdges uEdges edgesToConnect
      in
      if newDelta < inDelta then (newDelta, newTree, newMatrix)
      else addEdgeToSplitRecurse eEdges uEdges eTerminal uTerminal distMatrix (V.tail edgesToConnectVect) origTriple


-- | doSPRTBR takes split tree and rejoins by creating edges from a single one of the first edge set to each of the mebers of the second
-- important that e and u edges come in correct order and previous edges are e and u sets in order as well
doSPRTBR :: String -> Double -> V.Vector String -> Int -> SplitTreeData -> [TreeWithData]
doSPRTBR rejoinType curBestCost leafNames outGroup  split =
  -- trace ("In doSPRTBR") (
  let (eEdgeVect, uEdgeVect, delta, previousEdges, distMatrix) = split
  in
  if V.null eEdgeVect || V.null uEdgeVect then error "Empty edge vectors in doSPRTBR"
  else
      if fst3 (V.head eEdgeVect) == snd3 (V.head eEdgeVect) then -- if an OTU call readdTerminal
        let newLeafIndex = M.rows distMatrix
            -- keep whole uEdgeVect so recheck input tree edges
            additionList = V.map (addToEdgeSwap distMatrix (fst3 $ V.head eEdgeVect) (V.empty,uEdgeVect) newLeafIndex) uEdgeVect
            minAdditionCost = V.minimum (V.map fst3 additionList)
        in
        if minAdditionCost > delta then []
        else sieveTrees delta curBestCost additionList leafNames outGroup []
      else -- internal edge or edge with two OTUs as vertices
          -- check to make sure edges are where they should be can remove later
          -- fix e edge and join to each u edge for SPR (n^2)
          let edgesToConnect = createVectorEdgePairs rejoinType previousEdges eEdgeVect uEdgeVect
              -- e and u terminal should not be used here since OTUs are shortcircuited above
              eTerminal = (-1,-1,0)
              uTerminal = (-1,-1,0)
              additionList = V.map (addEdgeToSplit eEdgeVect uEdgeVect eTerminal uTerminal distMatrix) edgesToConnect
              minAdditionCost = V.minimum (V.map fst3 additionList)
          in
          if minAdditionCost > delta then []
          else sieveTrees delta curBestCost additionList leafNames outGroup []


-- | doSPRTBRSteep like doSPRTBR but only saves a sinlge (and better) tree
doSPRTBRSteep :: String -> Double -> V.Vector String -> Int -> SplitTreeData -> TreeWithData -> TreeWithData
doSPRTBRSteep rejoinType curBestCost leafNames outGroup split origTree@(_, inTree, _, inMatrix) =
  -- trace ("In doSPRTBR") (
  let (eEdgeVect, uEdgeVect, delta, previousEdges, distMatrix) = split
  in
  if V.null eEdgeVect || V.null uEdgeVect then error "Empty edge vectors in doSPRTBR"
  else
      if fst3 (V.head eEdgeVect) == snd3 (V.head eEdgeVect) then -- if an OTU call readdTerminal
        let newLeafIndex = M.rows distMatrix
            -- keep whole uEdgeVect so recheck input tree edges
            (newDelta, newTree, newMatrix) = addToEdgeSwapRecurse delta distMatrix (fst3 $ V.head eEdgeVect) (V.empty,uEdgeVect) newLeafIndex uEdgeVect
            newCost = curBestCost - delta + newDelta
            newickTree = convertToNewick leafNames outGroup newTree
            newickTree' = take (length newickTree - 3) newickTree ++ "[" ++ showDouble precision newCost ++ "]" ++ ";"
        in
        if newCost < curBestCost then (newickTree', newTree, newCost, newMatrix)
        else origTree
      else -- internal edge or edge with two OTUs as vertices
          -- check to make sure edges are where they should be can remove later
          -- fix e edge and join to each u edge for SPR (n^2)
          let edgesToConnect = createVectorEdgePairs rejoinType previousEdges eEdgeVect uEdgeVect
              -- e and u terminal should not be used here since OTUs are shortcircuited above
              eTerminal = (-1,-1,0)
              uTerminal = (-1,-1,0)
              (newDelta, newTree, newMatrix) = addEdgeToSplitRecurse eEdgeVect uEdgeVect eTerminal uTerminal distMatrix edgesToConnect (delta, inTree, inMatrix)
              newCost = curBestCost - delta + newDelta
              newickTree = convertToNewick leafNames outGroup newTree
              newickTree' = take (length newickTree - 3) newickTree ++ "[" ++ showDouble precision newCost ++ "]" ++ ";"
          in
          if newCost < curBestCost then
            --trace ("->" ++ show newCost)
            (newickTree', newTree, newCost, newMatrix)
          else origTree


-- | reAddTerminalsSteep like readdTerminals but only returns one tree keeping better
reAddTerminalsSteep :: String -> Double -> V.Vector String -> Int -> SplitTreeData -> TreeWithData -> TreeWithData
reAddTerminalsSteep rejoinType curBestCost leafNames outGroup split origTree =
  if rejoinType /= "otu" then error ("Incorrect swap function in reAddTerminals: " ++ rejoinType)
  else
    let (eEdgeVect, uEdgeVect, delta, _, distMatrix) = split
        nOTUs = V.length leafNames
    in
    if (V.length eEdgeVect > 1) || ((fst3 (V.head eEdgeVect) /= snd3 (V.head eEdgeVect)) && (fst3 (V.head eEdgeVect) < nOTUs) && (snd3 (V.head eEdgeVect) < nOTUs)) then origTree
    else if M.rows distMatrix /= ((2 * nOTUs) - 3) then error ("Dist Matrix incorrect size " ++ show (M.dim distMatrix) ++ " should be " ++ show ((2 * nOTUs) - 3, (2 * nOTUs) - 3))
    else
      let newLeafIndex = M.rows distMatrix
      -- take tail of uEdgeVect so not regerate input tree
          (newDelta, newTree, newMatrix) = addToEdgeSwapRecurse delta distMatrix (fst3 $ V.head eEdgeVect) (V.empty,uEdgeVect) newLeafIndex uEdgeVect -- (V.tail uEdgeVect) -- tail so not hit original tree, leave all to reestimate if necesary
          newCost = curBestCost - delta + newDelta
          newickTree = convertToNewick leafNames outGroup newTree
          newickTree' = take (length newickTree - 3) newickTree ++ "[" ++ showDouble precision newCost ++ "]" ++ ";"
      in
      if newCost < curBestCost then
        --trace ("->" ++ show newCost)
        (newickTree', newTree, newCost, newMatrix)
      else origTree



-- | filterNewTreesOnCost returns list of all unique new best cost trees from list
-- assumes curBestCost = cost of sabed trees
filterNewTreesOnCost :: Double -> [TreeWithData] -> [TreeWithData] -> [TreeWithData]
filterNewTreesOnCost curBestCost firstTreeList savedTrees =
  if null firstTreeList then savedTrees
  else
      let firstTree = head firstTreeList
          (_, _, firstCost, _) = firstTree
      in
      if firstCost < curBestCost then filterNewTreesOnCost firstCost (tail firstTreeList) [firstTree]
      else if firstCost > curBestCost then filterNewTreesOnCost curBestCost (tail firstTreeList) savedTrees
      else
          let uniqueTree = filterNewTrees savedTrees firstTree
          in
          if isNothing uniqueTree then filterNewTreesOnCost curBestCost (tail firstTreeList) savedTrees
          else filterNewTreesOnCost curBestCost (tail firstTreeList) (fromJust uniqueTree : savedTrees )

-- | filterNewTrees takes the first tree and checks if in the second list
filterNewTrees :: [TreeWithData] -> TreeWithData -> Maybe TreeWithData
filterNewTrees secondTreeList firstTree =
  if null secondTreeList then Just firstTree
  else
    let (firstNewick, _, _, _) = firstTree
        (secondNewick, _, _, _) = head secondTreeList
    in
    if firstNewick == secondNewick then Nothing
    else filterNewTrees (tail secondTreeList) firstTree

-- | getSaveNumber returns number ot save or Infty
getSaveNumber :: String -> Int
getSaveNumber inString =
  if length inString == 4 then maxBound :: Int
  else (read $ drop 5 inString) :: Int

-- | splitJoin does both split and rejoin operations in a fashion that if a better (shorter) tree is found is shortcircuits and
-- begins again on the new tree, else proceeds untill all splits and joins are completed, but only on a single tree
splitJoin :: (String -> Double -> V.Vector String -> Int -> SplitTreeData -> TreeWithData -> TreeWithData) -> String -> V.Vector String -> Int -> V.Vector Edge -> TreeWithData -> TreeWithData
splitJoin swapFunction refineType leafNames outGroup edgeVect curTreeWithData@(_, curTree, curTreeCost, curTreeMatrix) =
  if V.null edgeVect then curTreeWithData -- All splits tested, nothing better found
  else
    let firstEdge = V.head edgeVect
        firstSplit = splitTree curTreeMatrix curTree curTreeCost firstEdge
        firstTree@(_, firstNewTree, firstTreeCost, _) = swapFunction refineType curTreeCost leafNames outGroup firstSplit curTreeWithData
    in
    if firstTreeCost < curTreeCost then splitJoin swapFunction refineType leafNames outGroup (snd firstNewTree) firstTree
    else splitJoin swapFunction refineType leafNames outGroup (V.tail edgeVect) curTreeWithData


-- | splitJoinWrapper wraps around splitJoin to allow parallel execution
-- the reason is to allow the consumption of "edgeVect" recursively within the same tree
splitJoinWrapper :: (String -> Double -> V.Vector String -> Int -> SplitTreeData -> TreeWithData -> TreeWithData) -> String -> V.Vector String -> Int -> TreeWithData -> TreeWithData
splitJoinWrapper swapFunction refineType leafNames outGroup curTreeWithData@(_, curTree, _, _) =
    let edgeVect = snd curTree
    in
    splitJoin swapFunction refineType leafNames outGroup edgeVect curTreeWithData

-- | getGeneralSwapSteepestOne performs refinement as in getGeneralSwap but saves on a single tree (per split/swap) and
-- immediately accepts a Better (shorter) tree and resumes the search on that new tree
-- relies heavily on laziness of splitTree so not parallel at this level
getGeneralSwapSteepestOne :: String -> (String -> Double -> V.Vector String -> Int -> SplitTreeData -> TreeWithData -> TreeWithData) -> V.Vector String -> Int -> [TreeWithData] -> [TreeWithData] -> [TreeWithData]
getGeneralSwapSteepestOne refineType swapFunction leafNames outGroup inTreeList savedTrees =
  if null inTreeList then savedTrees
  else
      trace ("In "++ refineType ++ " Swap (steepest) with " ++ show (length inTreeList) ++ " trees with minimum length " ++ show (minimum $ fmap thd4 inTreeList)) (
      let steepTreeList = fmap (splitJoinWrapper swapFunction refineType leafNames outGroup) inTreeList `using` myParListChunkRDS
          steepCost = minimum $ fmap thd4 steepTreeList
      in
      --this to maintina the trajectories untill final swap--otherwise could converge down to single tree prematurely
      keepTrees steepTreeList "unique" "first" steepCost
      )

-- | getGeneralSwap performs a "re-add" of terminal identical to wagner build addition to available edges
-- performed on all splits recursively until no more better/equal cost trees found
-- this won't work to save "all" just unique and best and unique of best
-- add "steep-est" descent
getGeneralSwap :: String -> (String -> Double -> V.Vector String -> Int -> SplitTreeData -> [TreeWithData]) -> String -> String -> V.Vector String -> Int -> [TreeWithData] -> [TreeWithData] -> [TreeWithData]
getGeneralSwap refineType swapFunction saveMethod keepMethod leafNames outGroup inTreeList savedTrees =
  let maxNumSave = getSaveNumber saveMethod
  in
  if null inTreeList then savedTrees
  else
      trace ("In "++ refineType ++ " Swap with " ++ show (length inTreeList) ++ " trees with minimum length " ++ show (minimum $ fmap thd4 inTreeList)) (
      let curFullTree = head inTreeList
          overallBestCost = minimum $ fmap thd4 savedTrees
          (_, curTree, curTreeCost, curTreeMatrix) = curFullTree
          -- parallelize here
          splitTreeList = fmap (splitTree curTreeMatrix curTree curTreeCost) (V.toList $ snd curTree)  `using` myParListChunkRDS
          firstTreeList = fmap (swapFunction refineType curTreeCost leafNames outGroup) splitTreeList  `using` myParListChunkRDS
          firstTreeList' = filterNewTreesOnCost overallBestCost  (curFullTree : concat firstTreeList) savedTrees -- keepTrees (concat $ V.toList firstTreeList) saveMethod overallBestCost
      in
      -- Work around for negative NT.infinity tree costs (could be dst matrix issue)
      if NT.isInfinite curTreeCost || null firstTreeList' then getGeneralSwap refineType swapFunction saveMethod keepMethod leafNames outGroup (tail inTreeList) savedTrees else (
   let (_, _, costOfFoundTrees, _) = head firstTreeList'
   in
   --workaround for negatrive NT.infinity trees
   if NT.isInfinite costOfFoundTrees then getGeneralSwap refineType swapFunction saveMethod keepMethod leafNames outGroup (tail inTreeList) savedTrees
   else if costOfFoundTrees < overallBestCost then
     let uniqueTreesToAdd = fmap fromJust $ filter (/= Nothing ) $ fmap (filterNewTrees inTreeList) firstTreeList'
         treesToSwap = keepTrees (tail inTreeList ++ uniqueTreesToAdd) saveMethod keepMethod costOfFoundTrees
     in
     getGeneralSwap refineType swapFunction saveMethod keepMethod leafNames outGroup treesToSwap (take maxNumSave firstTreeList')
   else if costOfFoundTrees == overallBestCost then
     if length savedTrees >= maxNumSave then getGeneralSwap refineType swapFunction saveMethod keepMethod leafNames outGroup (tail inTreeList) savedTrees
     else
      let uniqueTreesToAdd = fmap fromJust $ filter (/= Nothing ) $ fmap (filterNewTrees inTreeList) firstTreeList'
          treesToSwap = keepTrees (tail inTreeList ++ uniqueTreesToAdd) saveMethod keepMethod costOfFoundTrees
      in
      getGeneralSwap refineType swapFunction saveMethod keepMethod leafNames outGroup treesToSwap (take maxNumSave $ savedTrees ++ firstTreeList')
   else getGeneralSwap refineType swapFunction saveMethod keepMethod leafNames outGroup (tail inTreeList) savedTrees)
        ) -- )

-- | performRefinement takes input trees in TRE and Newick format and performs different forms of tree refinement
-- at present just OTU (remove leaves and re-add), SPR and TBR
performRefinement :: String -> String -> String -> V.Vector String -> Int -> TreeWithData -> [TreeWithData]
performRefinement refinement saveMethod keepMethod leafNames outGroup inTree
  | refinement == "none" = [inTree]
  | refinement == "otu" =
    let newTrees = getGeneralSwapSteepestOne "otu" reAddTerminalsSteep leafNames outGroup [inTree] [([],(V.empty,V.empty), NT.infinity, M.empty)]
        newTrees' = getGeneralSwap "otu" reAddTerminals saveMethod keepMethod leafNames outGroup newTrees [([],(V.empty,V.empty), NT.infinity, M.empty)]
    in
    if refinement == "best:1" then
      if not $ null newTrees then newTrees
      else [inTree]
    else if not (null newTrees') then newTrees'
    else
      trace "OTU swap did not find any new trees"
      [inTree]
  | refinement == "spr" =
    let newTrees = getGeneralSwapSteepestOne "otu" reAddTerminalsSteep leafNames outGroup [inTree] [([],(V.empty,V.empty), NT.infinity, M.empty)]
        newTrees' = getGeneralSwapSteepestOne "spr" doSPRTBRSteep leafNames outGroup newTrees [([],(V.empty,V.empty), NT.infinity, M.empty)]
        newTrees'' = getGeneralSwap "spr" doSPRTBR saveMethod keepMethod leafNames outGroup newTrees' [([],(V.empty,V.empty), NT.infinity, M.empty)]
    in
    if saveMethod == "best:1" then
      if not $ null newTrees' then newTrees'
      else [inTree]
    else if not (null newTrees'') then newTrees''
    else
      trace "SPR swap did not find any new trees"
      [inTree]
  | refinement == "tbr" =
    let newTrees = getGeneralSwapSteepestOne "otu" reAddTerminalsSteep leafNames outGroup [inTree] [([],(V.empty,V.empty), NT.infinity, M.empty)]
        newTrees' = getGeneralSwapSteepestOne "spr" doSPRTBRSteep leafNames outGroup newTrees [([],(V.empty,V.empty), NT.infinity, M.empty)]
        newTrees'' = getGeneralSwapSteepestOne "tbr" doSPRTBRSteep leafNames outGroup newTrees' [([],(V.empty,V.empty), NT.infinity, M.empty)]
        newTrees''' = getGeneralSwap "tbr" doSPRTBR saveMethod keepMethod leafNames outGroup newTrees'' [([],(V.empty,V.empty), NT.infinity, M.empty)]
    in
    if saveMethod == "best:1" then
      if not $ null newTrees' then newTrees''
      else [inTree]
    else if not (null newTrees''') then newTrees'''
    else
      trace "TBR swap did not find any new trees"
      [inTree]
  | otherwise = errorWithoutStackTrace ("Unrecognized refinement method: " ++ refinement)
