{- |
Module      :  Utilities.hs
Description :  Module with useful functionsfor  distance tree construction methods dWag, Neightbor-Joining, UPGMA, and WPGMA
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

module Utilities.DistanceUtilities where

import qualified Data.Graph.Inductive.Graph        as G
import qualified Data.Graph.Inductive.PatriciaTree as P
import           Data.Maybe
import qualified Data.Number.Transfinite           as NT
import qualified Data.Set                          as Set
import qualified Data.Text.Lazy                    as T
import qualified Data.Vector                       as V
import           GeneralUtilities
import qualified GraphFormatUtilities              as PP
import           ParallelUtilities
import qualified SymMatrix                         as M
import           System.IO.Unsafe
import qualified System.Random                     as Rand
import qualified System.Random.Shuffle             as RandS
import           Types.DistanceTypes
--import qualified LocalSequence as LS
import qualified Data.Vector as LS
import qualified Data.List as L



-- | localRoundtakes a double multiplies by 10^precisoin, rounds to integer then divides
-- by precision
localRound :: Double -> Int -> Double
localRound val places =
  let factor = 10.0 ^^ places
      newVal = val * factor
      roundedVal = round newVal :: Int
  in
  fromIntegral roundedVal / factor

-- | showDouble integer number of string postion )including sign and decimal) and returns string
-- of that length--basically truncation
showDouble :: Int -> Double -> String
showDouble places val = show $ localRound val places -- reverse $ dropWhile (== '0') $ reverse $ take places $ show val

-- | test for equality with epsilon
withinEpsilon :: Double -> Double -> Bool
withinEpsilon a b = a == b
  {-
  if abs (a - b) < epsilon then True
  else False
  -}

-- | vertex2FGLNode take vertex of Int  and a Vector of Strings (leaf names) and returns
-- fgl node with type T.Text
vertex2FGLNode :: V.Vector String -> Vertex -> (Int, T.Text)
vertex2FGLNode leafVect vertIndex=
  if vertIndex < V.length leafVect then (vertIndex,  T.pack (leafVect V.! vertIndex))
  else
    let vertexName = "HTU" ++ show (vertIndex - V.length leafVect)
    in
    (vertIndex, T.pack vertexName)

-- | treeFGL take a Treee type and converts to an fgl graph
-- to be used with PhylParsers module hence Text
tree2FGL :: Tree -> V.Vector String -> P.Gr T.Text Double
tree2FGL inTree@(inVertexVect, inEdgeVect) leafNameVect =
  if inTree == emptyTree then error "Empty tree in tree2FGL"
  else
    let fglNodes = V.map (vertex2FGLNode leafNameVect) inVertexVect
    in
    -- trace ("tree2FGL orig, vertices, edges =  " ++ (show $ length inVertexVect ) ++ " "  ++ (show $ length fglNodes ) ++ " " ++ (show $ length inEdgeVect ))
    G.mkGraph (V.toList fglNodes) (V.toList inEdgeVect)

-- reIndexEdges take a list of vertices as integer indices and an edges (Int, Int, Double)
-- and retuns a new vertex with the indives of the vertex labels
reIndexEdges :: [Int] -> Edge -> Edge
reIndexEdges inVertexList (e,u,w) =
  if null inVertexList then error "Null inVertexList in reIndexEdges"
  else
    let newE = L.elemIndex e inVertexList
        newU = L.elemIndex u inVertexList
    in
    --un safe not testing but list is created from edges first
    (fromJust newE, fromJust newU, w)

-- | makeVertexNames takes vertgex indices and returns leaf name if < nOTUs and "HTU" ++ show Index
-- if not
makeVertexNames :: [Vertex] -> Int -> V.Vector String -> Bool -> [String]
makeVertexNames vertList nOTUs leafNames nameHTUs =
  if null vertList then []
  else
      let firstVert = head vertList
      in
      if firstVert < nOTUs then (leafNames V.! firstVert) : makeVertexNames (tail vertList) nOTUs leafNames nameHTUs
      else if nameHTUs then ("HTU" ++ show firstVert) : makeVertexNames (tail vertList) nOTUs leafNames nameHTUs
      else "" : makeVertexNames (tail vertList) nOTUs leafNames nameHTUs

-- | directSingleEdge takes an Int and makes that 'e' and otehr vertex as 'u' in edge (e->u)
directSingleEdge :: Int -> Edge -> Edge
directSingleEdge index (a,b,w)
  | a == index = (a,b,w)
  | b == index = (b,a,w)
  | otherwise = error ("Index " ++ show index ++ " doesn't match edge " ++ show (a,b,w) ++ " in directSingleEdge")

-- | getChildEdges returns the two edges that are childre of a vertex
getChildEdges :: Int -> Int -> V.Vector Edge -> V.Vector Edge
getChildEdges vertIndex nLeaves inEdgeVect
  | V.null inEdgeVect = V.empty
  | vertIndex < nLeaves = error ("Looking for child of leaf " ++ show (vertIndex, nLeaves))
  | otherwise =
    let (a,b,w) = V.head inEdgeVect
    in
    if (a == vertIndex) || (b == vertIndex) then V.cons (a,b,w) (getChildEdges vertIndex nLeaves (V.tail inEdgeVect)) else getChildEdges vertIndex nLeaves (V.tail inEdgeVect)


-- | directexEdges takes a vector of edges and outgrop index and directs the edges (parent -> child vertices) based on that
directEdges :: Int -> Int -> Bool -> V.Vector Edge -> V.Vector Edge
directEdges vertIndex nLeaves isFirst inEdgeVect
  | V.null inEdgeVect = V.empty
  | isFirst = --to find out group edge order larger to smaller will have outgroup index second
    let outgroupEdge = getEdgeRoot vertIndex inEdgeVect
        remainingEdgeVect = subtractVector (V.singleton outgroupEdge) inEdgeVect
        (a,b,w) = orderEdge outgroupEdge
    in
    V.cons (a,b,w) (directEdges a nLeaves False remainingEdgeVect)
  | vertIndex < nLeaves = V.empty
  | otherwise = -- not outgroup but regular node, get two child edges
    let descdendantEdges = getChildEdges vertIndex nLeaves inEdgeVect
        remainingEdgeVect = subtractVector descdendantEdges inEdgeVect
        newDescEdges = V.map (directSingleEdge vertIndex) descdendantEdges
    in
    if V.length newDescEdges /= 2 then error ("There should be 2 child edges for index " ++ show vertIndex ++ " and there are(is) " ++ show (V.length newDescEdges) ++ " " ++ show newDescEdges)
    else
        let (_, bf, _) = V.head newDescEdges
            (_, bs, _) = V.last newDescEdges
            firstSubEdges = directEdges bf nLeaves False remainingEdgeVect
            remainingEdgeVect' = subtractVector firstSubEdges remainingEdgeVect
            secondSubEdges = directEdges bs nLeaves False remainingEdgeVect'
        in
        (newDescEdges V.++ (firstSubEdges V.++ secondSubEdges))


-- | convertToGraph takes Vertex of names and a tree and return inductive Graph format
convertToDirectedGraph :: V.Vector String -> Int -> Tree -> P.Gr String Double
convertToDirectedGraph leafList outgroupIndex inTree =
  let (_, edgeVect) = inTree
      nOTUs = length leafList
      -- should be stright 0->n-1 but in case some vertex number if missing
      vertexList = L.sort $ Set.toList $ getVertexSet edgeVect
      vertexNames = makeVertexNames vertexList nOTUs leafList True
      labelledVertexList = L.zip vertexList vertexNames
      edgeList = V.toList $ directEdges outgroupIndex nOTUs True edgeVect
  in
  G.mkGraph labelledVertexList edgeList

-- | convertToGraphText takes Vertex of names and a tree and return inductive Graph format
convertToDirectedGraphText :: V.Vector String -> Int -> Tree -> P.Gr T.Text Double
convertToDirectedGraphText leafList outgroupIndex inTree =
  let (_, edgeVect) = inTree
      nOTUs = length leafList
      -- should be stright 0->n-1 but in case some vertex number if missing
      vertexList = L.sort $ Set.toList $ getVertexSet edgeVect
      vertexNames = makeVertexNames vertexList nOTUs leafList False
      labelledVertexList = L.zip vertexList (fmap T.pack vertexNames)
      edgeList = V.toList $ directEdges outgroupIndex nOTUs True edgeVect
  in
  G.mkGraph labelledVertexList edgeList


-- | convertToNewick generates a newick file by converting Tree type to FGL Graph,
-- adds a root and two new edges (deleting root edge)
-- and calls convert function from PhyloParsers
convertToNewick :: V.Vector String -> Int -> Tree -> String
convertToNewick leafNames outGroup inTree
  | inTree == emptyTree = error "Empty tree in convertToNewick"
  | V.null leafNames = error "Empty leaf names in convertToNewick"
  | otherwise =
    let fglTree = convertToDirectedGraphText leafNames outGroup inTree
    in
    --trace ("ConverttoNewick in-vertices in-edges" ++ (show $ length inVertexVect ) ++ " "  ++ (show $ V.toList inVertexVect ) ++ "\n" ++ (show $ length vertexVect ) ++ " "  ++ (show vertexVect ) ++ "\n" ++ (show $ length edgeVect ) ++ " " ++ show edgeVect ++ "\n" ++ show fglTree)
    PP.fglList2ForestEnhancedNewickString [fglTree] True False


-- | getEdgeRootIndex takes edge Vect, Index, and determines edges from root
getEdgeRootIndex :: Int -> Int -> V.Vector Edge -> (Int, Edge)
getEdgeRootIndex edgeIndex outgroup edgeVect =
  if V.null edgeVect then error "Root edge not found"
  else
   let (eVect, uVect, _) = V.head edgeVect
   in
   if (eVect == outgroup) || (uVect == outgroup) then (edgeIndex, V.head edgeVect)
  else getEdgeRootIndex (edgeIndex + 1) outgroup (V.tail edgeVect)



-- | convertToNewick wrapper to remove double commas
convertToNewick' :: V.Vector String -> Int -> Tree -> String
convertToNewick' leafNames outGroup wagTree = removeCrap $ convertToNewickGuts leafNames outGroup wagTree

-- | removeDoubleCommas removes second of double comas ",," -> ","
-- this a hack to fix problem in convertToNewick
removeCrap :: String -> String
removeCrap inString =
  if length inString == 1 then inString
  else
    let firstChar = head inString
        secondChar = inString !! 1
    in
    if firstChar == ',' && secondChar == ',' then ',' : removeCrap (drop 2 inString)
    else if firstChar == ',' && secondChar == ')' then ')' : removeCrap (drop 2 inString)
    else firstChar : removeCrap (tail inString)

-- | convertToNewick converts Tree rep to Newick String
-- includes edge cost--splits root edge cost into halves to help
-- tree viewers like FigTree (so input is unrooted)
-- NEED TO ADD smaller group left larger group right for more legible trees
convertToNewickGuts :: V.Vector String ->Int ->  Tree -> String
convertToNewickGuts leafNames outGroup wagTree =
  let (inLeaves, inEdges) = wagTree
      newEdges = fmap orderEdge inEdges
      (_, edgeVect) = orderTree (inLeaves, newEdges)
      foundEdge = getEdgeRoot outGroup edgeVect
  in
  let (firstVert, secondVert, weight) = foundEdge
      remainderEdges = V.filter (/= foundEdge) edgeVect
  in
  -- this is embarassing bullshit  -- converting  ",,"  to ","
  if firstVert == outGroup then "(" ++ (leafNames V.! outGroup)  ++ ":" ++ showDouble 8 (weight/2.0) ++ "," ++ getEdgesNonRoot secondVert remainderEdges (V.length leafNames) leafNames ++ ":" ++ showDouble 8 (weight/2.0) ++ ")"
  else "(" ++ (leafNames V.! outGroup)  ++ ":" ++ showDouble 8 (weight/2.0) ++ "," ++ getEdgesNonRoot firstVert remainderEdges (V.length leafNames) leafNames ++ ":" ++ showDouble 8 (weight/2.0) ++ ")"

-- | orderEdge takes an Edge and puts high index first then lower
orderEdge :: Edge -> Edge
orderEdge (a,b,w) =
  if a > b then (a,b,w)
  else (b,a,w)

-- | orderTree puts Tree edges in order based on edges
orderTree :: Tree -> Tree
orderTree (leaves, edges) =
  let edgeList = L.sort $ V.toList edges
  in
  (leaves, V.fromList edgeList)

-- | getEdges takes root Index and determines edges from root
getEdgeRoot :: Int -> V.Vector Edge -> Edge
getEdgeRoot edgeIndex edgeVect =
  if V.null edgeVect then (-1,-1,-1.0)
  else
   let (eVect, uVect, _) = V.head edgeVect
   in
   if (eVect == edgeIndex) || (uVect == edgeIndex) then V.head edgeVect else getEdgeRoot edgeIndex (V.tail edgeVect)

-- | getEdgesNonRoot takes root Index and determines edges from root returns String of Taxa
-- alwasy getting ordered edges and ordered tree
-- so eVect > uVect always; eVect can never be < nOTUs
--Need to add smaller tree left, bigger right
getEdgesNonRoot :: Int -> V.Vector Edge -> Int -> V.Vector String -> String
getEdgesNonRoot edgeIndex edgeVect nOTUs leafNames =
  --trace (show edgeIndex) (
  let terminal = (-1,-1,-1.0)
  in
  if V.null edgeVect then "END"
  else
   let thisEdge = getEdgeRoot edgeIndex edgeVect
   in
   if thisEdge == terminal then "ERROR"
   else
      let (eVect, uVect, weight) = thisEdge
          remainderEdges = V.filter (/= thisEdge) edgeVect
          eDesc = getEdgeRoot eVect remainderEdges
          uDesc = getEdgeRoot uVect remainderEdges
          eSubTree = getEdgesNonRoot eVect remainderEdges nOTUs leafNames
          uSubTree = getEdgesNonRoot uVect remainderEdges nOTUs leafNames
      in
      if V.null remainderEdges then
        (leafNames V.! uVect) ++ ":" ++ showDouble precision weight  ++ ","

      else if eVect == edgeIndex then

        if eVect < nOTUs then
          if uDesc /= terminal then "(" ++ (leafNames V.! eVect) ++ ":" ++ showDouble precision weight ++ "," ++ uSubTree ++ ")"
          else (leafNames V.! eVect) ++ ":" ++ showDouble precision weight ++ ","
        else
        if uVect < nOTUs then
          if eDesc /= terminal then "(" ++ (leafNames V.! uVect) ++ ":" ++ showDouble precision weight ++ "," ++ eSubTree ++  ")"
          else (leafNames V.! uVect) ++ ":" ++ showDouble precision weight ++ ","
        else
            if (eDesc /= terminal) && (uDesc == terminal) then eSubTree  ++ ":" ++ showDouble precision weight  ++ ","
            else if (eDesc == terminal) && (uDesc /= terminal) then uSubTree ++ ":" ++ showDouble precision weight ++ ","
            else
              if length eSubTree < length uSubTree then "(" ++ eSubTree ++ "," ++ uSubTree ++ ":" ++ showDouble precision weight ++ ")"
              else "(" ++ uSubTree ++ "," ++ eSubTree ++ ":" ++ showDouble precision weight ++ ")"

      else if uVect == edgeIndex then
        if uVect < nOTUs then
          if eDesc /= terminal then "(" ++ (leafNames V.! uVect) ++ ":" ++ showDouble precision weight ++  "," ++ eSubTree ++ ")"
          else (leafNames V.!uVect) ++ ":" ++ showDouble precision weight ++ ","

        else if eVect < nOTUs then
          if uDesc /= terminal then "(" ++ (leafNames V.! eVect) ++ ":" ++ showDouble precision weight ++ "," ++ uSubTree ++ ")"
          else (leafNames V.!eVect) ++ ":" ++ showDouble precision weight  ++ ","

        else
          if (eDesc /= terminal) && (uDesc == terminal) then eSubTree  ++ ":" ++ showDouble precision weight ++ ","
          else if (eDesc == terminal) && (uDesc /= terminal) then uSubTree ++ ":" ++ showDouble precision weight ++ ","
          else
              if length eSubTree < length uSubTree then "(" ++ eSubTree ++ "," ++ uSubTree ++ ":" ++ showDouble precision weight ++ ")"
              else "(" ++ uSubTree ++ ":" ++ showDouble precision weight ++ "," ++ eSubTree ++ ")"

      else getEdgesNonRoot edgeIndex remainderEdges nOTUs leafNames ++ ":" ++ showDouble precision weight ++ ","

-- getBestTrees takes newick and L.sorts on comment at end with cost
getBestTrees :: String -> Int -> [TreeWithData] -> Double -> [TreeWithData] -> [TreeWithData]
getBestTrees keepMethod number inList curBestCost curBestTrees =
  if null inList then
      -- apply keep method if not keeping all
      if length curBestTrees <= number then curBestTrees
      else if keepMethod == "first" then take number curBestTrees
      else if keepMethod == "last" then reverse $ take number $ reverse curBestTrees
      else if keepMethod == "random" then
        -- not Fitch but shuffles elements and takes first n
        let randIntList = betterRandomList (length curBestTrees) (length curBestTrees - 1)
            newList =  RandS.shuffle curBestTrees randIntList
        in
        take number newList
      else error ("Keep method " ++ keepMethod ++ " not implemented")
  else
    let firstTree = head inList
        (_, _, firstCost, _) = firstTree
    in
    if firstCost < curBestCost then getBestTrees keepMethod number (tail inList) firstCost [firstTree]
    else if withinEpsilon firstCost curBestCost then getBestTrees keepMethod number (tail inList) firstCost (firstTree : curBestTrees)
    else getBestTrees keepMethod number (tail inList) curBestCost curBestTrees

-- | getUniqueTrees saves uniqe newick trees
-- different paths--ehnce different distMatrices coud result in same newick tree
getUniqueTrees :: [TreeWithData] -> [TreeWithData] -> [TreeWithData]
getUniqueTrees inList uniqueList =
  if null inList then uniqueList
  else
    let firstTree = head inList
        fstNewick = fst4 firstTree
    in
    if fstNewick `notElem` fmap fst4 uniqueList then getUniqueTrees (tail inList) (firstTree : uniqueList)
    else getUniqueTrees (tail inList) uniqueList

-- | keepTrees filters newick trees based on options
-- all keep all
-- best shortest (and unique) allows number of max to save
-- unique unique  representations irespective of length
-- keep metyhod for save first | last | atRandom if buffer full
keepTrees :: [TreeWithData] -> String -> String -> Double -> [TreeWithData]
keepTrees inList saveMethod keepMethod curBestCost
  | null inList = []
  | saveMethod == "all" = inList
  | take 6 saveMethod == "unique" =
    if length saveMethod == 6 then getUniqueTrees inList []
    else
      let number = read (drop 7 saveMethod) :: Int
      in
      take number $ L.sortOn thd4 $ getUniqueTrees inList []
  | take 4 saveMethod == "best" =
    if length saveMethod == 4 then getUniqueTrees (getBestTrees keepMethod (maxBound :: Int) inList curBestCost []) []
    else if (saveMethod !! 4) == ':' then
      let number =  read (drop 5 saveMethod) :: Int
          saveTrees = take number $ L.sortOn thd4 $ getUniqueTrees inList []
          (_, _, bestCost, _) = head saveTrees
      in
      getBestTrees keepMethod number saveTrees bestCost []
    else error ("Save method " ++ saveMethod ++ " improperly formatted")
  | otherwise = error ("Save method " ++ saveMethod ++ " not implemented")

  -- | ranList generates random list of positive integers
ranList :: Rand.StdGen -> Int -> Int -> [Int]
ranList sg n maxValue = take n $ Rand.randomRs (0,maxValue) sg

-- | driver function to generate list of positive integers
{-# NOINLINE betterRandomList #-}
betterRandomList :: Int -> Int -> [Int]
betterRandomList n maxValue= unsafePerformIO $ do
    sg <- Rand.getStdGen
    return $ ranList sg n maxValue

-- | Subtrace vector subtracts elements of vector a from vector b
-- is thins n^2 ?
-- edges are directed
subtractVector :: (Eq a) => V.Vector a -> V.Vector a -> V.Vector a
subtractVector a b
  | V.null a = b
  | V.null b = V.empty
  | otherwise =
    let firstB = V.head b
        notFound = V.notElem firstB a
    in
    if notFound then V.cons firstB (subtractVector a (V.tail b))
    else subtractVector a (V.tail b)

-- | getVertexSet take a vector of edges and creates the set of vertex numbers
getVertexSet :: V.Vector Edge -> Set.Set Vertex
getVertexSet edgeVect =
  if V.null edgeVect then Set.empty
  else
      let (a,b,_) = V.head edgeVect
          thisSet =  Set.fromList [a,b]
      in
      Set.union thisSet (getVertexSet $ V.tail edgeVect)

-- | getMinRowDistMatrix distMatrix tabuList
getMinRowDistMatrix :: M.Matrix Double -> [Int] -> (Int, Double) -> Int -> Int -> (Int, Int, Double)
getMinRowDistMatrix distMatrix tabuList minPair@(minCol, minVal) curColumn row
  | curColumn == LS.length (distMatrix LS.! row) = (row, minCol, minVal)
  | row `elem` tabuList = (-1,-1, NT.infinity)
  | curColumn == row = getMinRowDistMatrix distMatrix tabuList minPair (curColumn + 1) row
  | curColumn `elem` tabuList = getMinRowDistMatrix distMatrix tabuList minPair (curColumn + 1) row
  | otherwise =
    let firstVal = distMatrix M.! (row, curColumn)
    in
    if firstVal < minVal then getMinRowDistMatrix distMatrix tabuList (curColumn, firstVal) (curColumn + 1) row
    else getMinRowDistMatrix distMatrix tabuList minPair (curColumn + 1) row


-- | compareTriples take two triples and orders them based on the smaller thrid element
minTriples :: (Ord c) => (a,b,c) -> (a,b,c) -> Ordering
minTriples (_,__,c) (_,_,d)
  | c < d = LT
  | c > d = GT
  | otherwise = EQ

-- | getMatrixMinPairTabu' takes distMatrix initial integer pair and value
-- traverses the matrix (skippiong rows and columns in tabuList and return minimum distance and index pair
-- if tie takes first
-- gets minimum by row and parallelizes ovcer rows
  -- call with (-1, -1, NT.infinity) 0 0
getMatrixMinPairTabu :: M.Matrix Double -> [Int] -> (Int, Int, Double)
getMatrixMinPairTabu distMatrix tabuList =
  if M.null distMatrix then error "Empty matrix in getMatrixPairTabu"
  else
    let minValueList = seqParMap myStrategy (getMinRowDistMatrix distMatrix tabuList (-1, NT.infinity) 0) [0..(M.rows distMatrix - 1)]
    in
    L.minimumBy minTriples minValueList

-- | getMatrixMinPairTabu takes distMatrix initial integer pair and value
-- traverses the matrix (skippiong rows and columns in tabuList and return minimum distance and index pair
-- if tie takes first
  -- call with (-1, -1, NT.infinity) 0 0
getMatrixMinPairTabu' :: M.Matrix Double -> [Int] -> (Int, Int, Double) -> Int -> Int -> (Int, Int, Double)
getMatrixMinPairTabu' distMatrix tabuList curBest curRow curColumn
  | curRow == M.rows distMatrix = curBest
  | curColumn == M.cols distMatrix = getMatrixMinPairTabu' distMatrix tabuList curBest (curRow + 1) 0
  | curColumn == curRow = getMatrixMinPairTabu' distMatrix tabuList curBest curRow (curColumn + 1)
  | (curColumn `elem` tabuList) || (curRow `elem` tabuList) = getMatrixMinPairTabu' distMatrix tabuList curBest curRow (curColumn + 1)
  | otherwise =
    let (_, _, currentBestDistance) = curBest
    in
    if  distMatrix M.! (curRow, curColumn) < currentBestDistance then
      getMatrixMinPairTabu' distMatrix tabuList (curRow, curColumn, distMatrix M.! (curRow, curColumn)) curRow (curColumn + 1)
    else getMatrixMinPairTabu' distMatrix tabuList curBest curRow (curColumn + 1)



-- | getMatrixMinPair takes distMatrix initla pinteger pair and value
-- traverses teh matrix and return minimum distance and index pair
-- if tie takes first
  -- call with (-1, -1, NT.infinity) 0 0
getMatrixMinPair :: M.Matrix Double ->  (Int, Int, Double) -> Int -> Int -> (Int, Int, Double)
getMatrixMinPair distMatrix curBest curRow curColumn
  | curRow == M.rows distMatrix = curBest
  | curColumn == M.cols distMatrix = getMatrixMinPair distMatrix curBest (curRow + 1) 0
  | curColumn == curRow = getMatrixMinPair distMatrix curBest curRow (curColumn + 1)
  | otherwise =
    let (_, _, currentBestDistance) = curBest
    in
    if  distMatrix M.! (curRow, curColumn) < currentBestDistance then
      getMatrixMinPair distMatrix (curRow, curColumn, distMatrix M.! (curRow, curColumn)) curRow (curColumn + 1)
    else getMatrixMinPair distMatrix curBest curRow (curColumn + 1)


-- | getMatrixMaxPair takes distMatrix initla pinteger pair and value
-- traverses teh matrix and return maximum distance and index pair
-- if tie takes first
-- call with (-1 , -1, 0 :: Double) 0 0
getMatrixMaxPair :: M.Matrix Double ->  (Int, Int, Double) -> Int -> Int -> (Int, Int, Double)
getMatrixMaxPair distMatrix curBest curRow curColumn
  | curRow == M.rows distMatrix = curBest
  | curColumn == M.cols distMatrix = getMatrixMaxPair distMatrix curBest (curRow + 1) 0
  | curColumn == curRow = getMatrixMaxPair distMatrix curBest curRow (curColumn + 1)
  | otherwise =
    let (_, _, currentBestDistance) = curBest
    in
    if  distMatrix M.! (curRow, curColumn) > currentBestDistance then
      getMatrixMaxPair distMatrix (curRow, curColumn, distMatrix M.! (curRow, curColumn)) curRow (curColumn + 1)
    else getMatrixMaxPair distMatrix curBest curRow (curColumn + 1)

-- | getTreeCost takes Tree and returns cost based on sum of edge weights
getTreeCost :: Tree -> Double
getTreeCost inTree =
  V.sum $ V.map getEdgeCost $ snd inTree

-- | getEdgeCost returns weight form edge tuple
getEdgeCost :: (Vertex, Vertex, Weight) -> Double
getEdgeCost (_, _, edgeWeight) = edgeWeight


