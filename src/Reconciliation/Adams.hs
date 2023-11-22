{- |
Functions to create Adams II consensus trees (unlabelled internal veritices)
-}

module Reconciliation.Adams (makeAdamsII) where

import Control.Evaluation
import Control.Evaluation.Verbosity (Verbosity (..))
import Control.Monad.IO.Class (MonadIO (..))
import Data.Graph.Inductive.Graph qualified as G
import Data.Graph.Inductive.PatriciaTree qualified as P
import Data.List qualified as L
import Data.Maybe
import Data.Set qualified as Set
import Data.Text.Lazy qualified as T
import Data.Vector qualified as V
import GraphFormatUtilities qualified as PhyP
import Types.Types
-- import ParallelUtilities qualified as PU
-- import Debug.Trace

data VertexType = Root | Internal | Leaf | Network | Tree
    deriving stock (Eq, Ord, Read, Show) --NFData ?

type Vertex = (String, [Int], [Int], VertexType) -- name, child vertices, parent vertices, type (a bit redundant)
type Edge = (Int, Int, Maybe Double) -- terminal vertices (by numbers) and potential length
type PhyloGraphVect = (V.Vector Vertex, V.Vector Edge)
type GenPhyNetNode = (String, [String], [String]) --Make list for both so unresolved network (node name, [descendant name], [ancestor name])
type GenPhyNet = [GenPhyNetNode]

-- | null PhyloGraphVect
nullGraphVect :: PhyloGraphVect
nullGraphVect = (V.empty, V.empty)

-- | getAdamsIIPair inputs 2 PhyloGraphVects and returns AdamsII consensus
getAdamsIIPair ::  PhyloGraphVect -> PhyloGraphVect -> PhyloGraphVect
getAdamsIIPair inGraphVectA inGraphVectB =
        let inGraphVectList = [inGraphVectA, inGraphVectB]
            (sameLeafSet, leafSets) = getAndCheckLeafSets inGraphVectList
            curVertexSets = map fst inGraphVectList
            rootPairList = map (findRoot 0) inGraphVectList
            --rootIndexList = map fst rootPairList
            rootVertexList = map snd rootPairList
            rootSplits = map getSecond rootVertexList
            rootSplitLeafListList = getSplitLeafListList rootSplits inGraphVectList
            rootLUBPre = leastUpperBound rootSplitLeafListList
            rootLUB = [L.sort x | x <- rootLUBPre, not (null x)] --need map sort $

            --create nodes based on LUBs
            rootNode :: (String, [String], [a])
            rootNode = ("root", map lub2TreeRep rootLUB, [])
            leavesPlaced = concat [x | x <- rootLUB, length x < 3]
            vertexLeafSetList = map (map getLeafSetFromNodeName . V.toList) curVertexSets
            potentialVertexSets = map (map getSecond . V.toList) curVertexSets
        in
        if not sameLeafSet then errorWithoutStackTrace ("Leaf sets of input graphs do not match"  <> show leafSets)
        else
          --return processed when have all nodes
            let allAdamsNodes = makeAdamsNodes [rootNode] "root" rootLUB leavesPlaced (zip potentialVertexSets vertexLeafSetList) --curVertexSets vertexLeafSetList
            in
            genPhyNet2PhyloGraphVect allAdamsNodes

-- | mkGraphPair take LNode list lEdge List pair and reutns fgl graph
mkGraphPair :: ([G.LNode a], [G.LEdge b]) -> P.Gr a b
mkGraphPair (nodeList, edgeList) = G.mkGraph nodeList edgeList


-- | makeAdamsII takes a list of fgl graphs, convertes them to PhyloGraphVect
-- makes the Adamns consensus and then converts back to fgl for return to EUN code
makeAdamsII :: [G.LNode String] -> [P.Gr String Double] ->  PhyG (P.Gr String Double)
makeAdamsII leafNodeList inFGList
  | null leafNodeList = error "Null leaf node list in makeAdamsII"
  | null inFGList = pure G.empty
  | otherwise =
    let inGraphNodes = fmap G.labNodes inFGList
        inGraphEdges = fmap G.labEdges inFGList
        inGraphNonLeafNodes = fmap (drop $ length leafNodeList) inGraphNodes
        newNodeListList = fmap (leafNodeList <> ) inGraphNonLeafNodes
        inFGList' = mkGraphPair <$> zip  newNodeListList inGraphEdges

        -- parallel 
        action :: (Ord a) => P.Gr a b -> Bool
        action = isTree

    in do
        actionPar <- getParallelChunkMap
        let allTreesList = actionPar action inFGList'
        -- let allTreesList = PU.seqParMap PU.myStrategyRDS isTree inFGList' -- `using` PU.myParListChunkRDS
        let allTrees = L.foldl1' (&&) allTreesList
        
        if not allTrees then errorWithoutStackTrace ("Input graphs are not all trees in makeAdamsII: " <> show allTreesList)
        else if not (leafSetConstant [] inFGList) then errorWithoutStackTrace "Input leaf sets not constant in makeAdamsII"
        else
          let inPGVList = fmap fgl2PGV inFGList' -- paralle problem with NFData seqParMap myStrategy fgl2PGV inFGList
              adamsPGV = L.foldl1' getAdamsIIPair inPGVList
          in
          -- trace ("\nAdams: " <> show adamsPGV)
          pure $ pgv2FGL adamsPGV

-- | fgl2PGVEdge takes an fgl Labeled edge (e,u,label)
-- and returns a PGV edge  with no brnach length (e,u,Nothing)
fgl2PGVEdge :: G.LEdge b -> Edge
fgl2PGVEdge (e, u, _) = (e, u, Nothing)

-- | fgl2PGVNode takes a tripple of an fgl labelled node (Int, a), its
-- child verteces and parent versitices and returns the PGV Vertex including its type
fgl2PGVNode :: (Show b) => P.Gr T.Text b -> (G.LNode String, [Int], [Int]) -> Vertex
fgl2PGVNode inGraph ((index, inLabel), childList, parentList) =
    --if null label then error "Null node label in fgl2PGVNode"
    --else
    let guts = T.init $ PhyP.component2Newick inGraph False False (index, T.pack inLabel)
        label = if PhyP.checkIfLeaf guts then T.unpack $ T.tail $ T.init guts else T.unpack guts
    in
    if null parentList then (label, childList, parentList, Root)
    else if null childList then (label, childList, parentList, Leaf)
    else if length parentList == 1 then (label, childList, parentList, Reconciliation.Adams.Tree)
    else if length parentList > 1 then (label, childList, parentList, Network)
    else error "This can't happen in fgl2PGVNode"

-- | fgl2PGV takes an fgl (functional graph) and convertes to PhyloGraphVect
-- to use local (and old) Adams consensus functions
-- retuns "Nothing" for edge labels ( no need for branch lengths)
fgl2PGV :: P.Gr String Double -> PhyloGraphVect
fgl2PGV inGraph =
    if G.isEmpty inGraph then nullGraphVect
    else
        let fglNodeList = G.labNodes inGraph
            fglNodeParentList = fmap (G.pre inGraph . fst) fglNodeList
            fglNodeChildList = fmap (G.suc inGraph . fst) fglNodeList
            fglNodeInfoList = zip3 fglNodeList fglNodeChildList fglNodeParentList
            fglEdgeList = G.labEdges inGraph
            pgvNodeList = fmap (fgl2PGVNode (PhyP.stringGraph2TextGraph inGraph)) fglNodeInfoList
            pgvEdgeList = fmap fgl2PGVEdge fglEdgeList
        in
        (V.fromList pgvNodeList, V.fromList pgvEdgeList)

-- | vertex2Node take alist of vertices and returns a list of fgl Labelled nodes
vertex2Node :: Int -> V.Vector Vertex -> [G.LNode String]
vertex2Node counter inVertexVect =
    if V.null inVertexVect then []
    else
        let (label, _, _, _) = V.head inVertexVect
        in
        (counter, label) : vertex2Node (counter + 1) (V.tail inVertexVect)

-- | edge2FGLEdge take vertex of Int Int Maybe Double and returns
-- fgl node with type Double
edge2FGLEdge :: (Int, Int, Maybe Double) -> (Int, Int, Double)
edge2FGLEdge (e, u, _) = (e,u, 1.0 :: Double)

-- | pgv2FGL take a PhyloGraphVect and converts to an fgl graph
pgv2FGL :: PhyloGraphVect -> P.Gr String Double
pgv2FGL (inVertexVect, inEdgeVect) =
    let fglNodes = vertex2Node 0 inVertexVect
        fglEdges = V.map edge2FGLEdge inEdgeVect
    in
    G.mkGraph fglNodes (V.toList fglEdges)


-- | getLeafList returns sorted leaf complement of graph fgl
getLeafLabelListFGL ::  (Ord a) => P.Gr a b -> [a]
getLeafLabelListFGL inGraph =
  if G.isEmpty inGraph then error "Empty graph in getLeafLabelListFGL"
  else
    let degOutList = G.outdeg inGraph <$> G.nodes inGraph
        newNodePair = zip degOutList (G.labNodes inGraph)
        leafPairList = filter ((== 0).fst ) newNodePair
        (_, leafList) = unzip leafPairList
    in
    L.sort $ snd <$> leafList

-- | leafSetConstant takes a series of fgl graphs and checks if leaf sets are the same for
-- all of them
leafSetConstant :: (Ord a) => [a] -> [P.Gr a b] -> Bool
leafSetConstant leafList inFGLList
  | null inFGLList = True
  | null leafList =
    -- first graph
    let firstGraph = head inFGLList
        firstLeaves = getLeafLabelListFGL firstGraph
    in
    leafSetConstant firstLeaves (tail inFGLList)
      | otherwise =
    let thisGraph = head inFGLList
        theseLeaves = getLeafLabelListFGL thisGraph
    in
    theseLeaves == leafList && leafSetConstant leafList (tail inFGLList)

-- | isTree takes fgl graph and checks is conected, no self edges, single root (includes connected), no indegree
-- > 1 nodes, leaf labels appear only once
isTree :: (Ord a) => P.Gr a b -> Bool
isTree inGraph =
  not (G.isEmpty inGraph) && (
  let nodeIndegrees = G.indeg inGraph <$> G.nodes inGraph
      maxIndegree = maximum nodeIndegrees
      rootNodes =  filter (== 0) nodeIndegrees
      leafLabels = getLeafLabelListFGL inGraph
      uniqueLeafLabels = L.nub leafLabels
      eList = fst <$> G.edges inGraph
      uList = snd <$> G.edges inGraph
      selfList = filter id $ zipWith (==) eList uList
  in
  not ((((length rootNodes /= 1) || (maxIndegree > 1)) || (length leafLabels /= length uniqueLeafLabels)) || not (null selfList)))

-- | getRootNamesFromGenPhyNet extracts non-leaf-non-root
-- names from vertices in order found
getRootNamesFromGenPhyNet :: GenPhyNet  -> [String]
getRootNamesFromGenPhyNet inNet =
    if null inNet then []
    else
        let (name, desc, anc) = head inNet
            vertType = getVertType (length desc) (length anc)
        in
        if vertType == Root then
            name : getRootNamesFromGenPhyNet (tail inNet)
        else
            getRootNamesFromGenPhyNet (tail inNet)

-- | getNonLeafNonRootNamesFromGenPhyNet extracts non-leaf-non-root
-- names from vertices in order found
getNonLeafNonRootNamesFromGenPhyNet :: GenPhyNet  -> [String]
getNonLeafNonRootNamesFromGenPhyNet inNet =
    if null inNet then []
    else
        let (name, desc, anc) = head inNet
            vertType = getVertType (length desc) (length anc)
        in
        if (vertType /= Leaf) && (vertType /= Root) then
            name : getNonLeafNonRootNamesFromGenPhyNet (tail inNet)
        else
            getNonLeafNonRootNamesFromGenPhyNet (tail inNet)

-- | getLeafNamesFromGenPhyNet extracts leaf names from vertices in order found
getLeafNamesFromGenPhyNet :: GenPhyNet  -> [String]
getLeafNamesFromGenPhyNet inNet =
    if null inNet then []
    else
        let (name, desc, anc) = head inNet
        in
        if getVertType (length desc) (length anc) == Leaf then
            name : getLeafNamesFromGenPhyNet (tail inNet)
        else
            getLeafNamesFromGenPhyNet (tail inNet)

-- | getVertType takes list of desc and anc to determine type of vertex
getVertType :: Int -> Int -> VertexType
getVertType nDesc nAnc
  | nDesc == 0 && nAnc == 0 = error "Isolated node"
  | nAnc == 0 = Root
  | nDesc == 0 = Leaf
  | nAnc == 1 = Reconciliation.Adams.Tree
  | nAnc > 2 = Network
  | otherwise = error ("Screwey node: indegree " <> show nDesc <> " outdegree " <> show nAnc)

-- | getVertNum takes a list of vertex names and teh complete list and
-- returns a list of the indices (integers) of the names
getVertNum :: [String] -> [String] -> [Int]
getVertNum nameList vertexNameList =
    if null vertexNameList then []
    else
        let firstVertName = head vertexNameList
            vertNum = L.elemIndex firstVertName nameList
        in
        if isNothing vertNum then error ("Error in vertex name index: " <> show firstVertName <> " in " <> show nameList)
        else
            fromJust vertNum : getVertNum nameList (tail vertexNameList)

-- | oldIndex2New takes a PhyloGraphVect and creates a list of reorder based ordered name list
oldIndex2New :: V.Vector Vertex  -> [String] -> [(Int, Vertex)]
oldIndex2New inVertexVect nameList =
    if V.null inVertexVect then []
    else
        let curVert = V.head inVertexVect
            (vertName, _, _, _) = curVert
            vertNum = L.elemIndex vertName nameList
        in
        (fromJust vertNum, curVert) : oldIndex2New (V.tail inVertexVect) nameList

-- | genForestToPhyloGraph converts GenForest to PhyloGraph (so can use legacy
-- ENewick etc parsers
-- takes flattened vector of GenPhyNetNodes and builds vertices (leaf, internal,
-- and root) and edges.  Vertices and edges are added to input null
-- PhyloGraphVect
-- EDGES SEEM TO BE INCORRECT IN PLACES
genForestToPhyloGraphVect :: V.Vector GenPhyNetNode -> PhyloGraphVect -> [String] -> PhyloGraphVect
genForestToPhyloGraphVect inGen inPhyVect nameList =
    if V.null inGen then inPhyVect
    else
        let (inVertexName, inVertexDescNameList, inVertexAncNameList) = V.head inGen
            (curVertVect, curEdgeVect) = inPhyVect
            descNumList = getVertNum nameList inVertexDescNameList
            ancNumList = getVertNum nameList inVertexAncNameList
            vertType = getVertType (length descNumList) (length ancNumList)
            newEdgeVect ::  V.Vector (Int, Int, Maybe a)
            newEdgeVect = V.zip3 (V.fromList ancNumList) (V.replicate (length ancNumList)
                (head $ getVertNum nameList [inVertexName]))
                (V.replicate (length ancNumList) Nothing) --edge from anc to current, no weight info
        in
        genForestToPhyloGraphVect (V.tail inGen)
            (V.snoc curVertVect (inVertexName, descNumList, ancNumList, vertType),
            curEdgeVect <> newEdgeVect) nameList

-- | getNamesFromGenPhyNet extracts names from vertices in order found
-- leaves are first, then internal, root last
getNamesFromGenPhyNet :: GenPhyNet  -> [String]
getNamesFromGenPhyNet inNet =
    L.sort (getLeafNamesFromGenPhyNet inNet) <> getNonLeafNonRootNamesFromGenPhyNet inNet
        <> getRootNamesFromGenPhyNet inNet

-- | getShortestList takes list and length and list of lists and return
-- shortest list
getShortestList :: ([a], [b]) -> Int -> [([a],[b])] -> ([a],[b])
getShortestList bestList lengthBestList inListList =
    if null inListList then bestList
    else
        let curList = head inListList
            lengthCurList = length $ fst curList
        in
        if lengthCurList < lengthBestList then getShortestList curList lengthCurList (tail inListList)
        else  getShortestList bestList lengthBestList (tail inListList)

-- | getSplitList take an LUB, list of placed taxa, and vector of tree vertices
-- and returns a list of splits for each input tree (from tree vertices)
-- also filters out placed taxa
-- CLEANUP--many more ioperations than needed--should be passed as better
-- structure
getSplitList :: [String] -> [String] -> ([[Int]], [[String]]) -> [[String]]
getSplitList curLUB placedTaxa (potentialVerts, vertLeafSet) =
    if null curLUB then error "Null LUB in getSplitList"
    else
       let  vertList =  [(x, y) | (x, y)  <- zip vertLeafSet potentialVerts, L.intersect x curLUB == curLUB]
            smallestVert = snd $ getShortestList (head vertList) (length $ fst $ head vertList) (tail vertList)
            vectVertLeafSet = V.fromList vertLeafSet --adds factor of "n", could pass another variable?
            rawLUBs = map (vectVertLeafSet V.!) smallestVert
            newLUBs = map (L.\\ placedTaxa) rawLUBs
        in
        newLUBs

-- | replaceChar take set of charcters to be replaced by a char in a String
replaceChar :: String -> Char -> Char -> Char
replaceChar inSet2Replace replaceChar2 inChar =
    if inChar `elem` inSet2Replace then replaceChar2
    else inChar

-- | getVertsFromIndexList takes a list of vertex vector indices and returns a list of
-- vertices
getVertsFromIndexList :: [Int] -> PhyloGraphVect -> [Vertex]
getVertsFromIndexList indexList inGraphVect  =
    if null indexList then []
    else
        let (vertVect, _) = inGraphVect
        in
        (vertVect V.! head indexList) :  getVertsFromIndexList (tail indexList) inGraphVect

-- | ggenPhyNet2PhyloGraphVect takes as input GenPhyNet and return
-- PhyloGraphVect with root as last node
genPhyNet2PhyloGraphVect :: GenPhyNet -> PhyloGraphVect
genPhyNet2PhyloGraphVect inGenPhyNet =
    if null inGenPhyNet then error "Null GenPhyNet in genPhyNet2PhyloGraphVect"
    else
        let nameList = getNamesFromGenPhyNet inGenPhyNet
            (vertVect, edgeVect) = genForestToPhyloGraphVect (V.fromList inGenPhyNet) nullGraphVect nameList
            newVertVect = vertVect V.// oldIndex2New vertVect nameList
        in
        (newVertVect, edgeVect)

-- | makeAdamsNodes takes root Adams node, rootLUB, vertex sets of input
-- trees and placed leaf set and constructs each Adams node in turn.
makeAdamsNodes :: GenPhyNet -> String -> [[String]] -> [String] -> [([[Int]], [[String]])] -> GenPhyNet
makeAdamsNodes inAdamsTree parentName inLUBList placedTaxa bothLeafLists = --inTreeVertexLists vertexLeafSetList =
   if null inLUBList then inAdamsTree
   else
        let curLUB = head inLUBList
        in
        if length curLUB == 1 then --make nodes since done
            let newNode :: (String, [a], [String])
                newNode = (head curLUB, [], [parentName])
            in  makeAdamsNodes
                  (newNode : inAdamsTree)
                  parentName
                  (tail inLUBList)
                  (head curLUB : placedTaxa) bothLeafLists --inTreeVertexLists vertexLeafSetList
        else if length curLUB == 2 then
            let leftChild  = lub2TreeRep  [head curLUB]
                rightChild = lub2TreeRep  [last curLUB]
                newNode1, newNode2, newNode3 :: (String, [String], [String])
                newNode1 = (lub2TreeRep curLUB, [leftChild, rightChild], [parentName])
                newNode2 = ( leftChild, [], [lub2TreeRep curLUB])
                newNode3 = (rightChild, [], [lub2TreeRep curLUB])
                newGenPhyNet = newNode2 : (newNode3 : (newNode1 : inAdamsTree))
                newPlacedTaxa = lub2TreeRep curLUB : (leftChild : (rightChild : placedTaxa))
            in
            makeAdamsNodes newGenPhyNet parentName (tail inLUBList)
                newPlacedTaxa bothLeafLists --inTreeVertexLists vertexLeafSetList
        else --core case with LUB creation and taxon placementg
            let splitListList = map (getSplitList curLUB placedTaxa) bothLeafLists --(zip inTreeVertexLists vertexLeafSetList)
                newLUBpre = leastUpperBound splitListList
                newLUB =  [L.sort x | x <- newLUBpre, not (null x)] --had "map sort $" was this "sort" necessary?  for List intersection?
                newNode = (lub2TreeRep curLUB, map lub2TreeRep newLUB, [parentName])
            in
            --trace ("New LUBs " <> show newLUB <> " newNode " <> show newNode)
            makeAdamsNodes (newNode : inAdamsTree) (lub2TreeRep curLUB) (tail inLUBList)
                placedTaxa bothLeafLists <> --inTreeVertexLists vertexLeafSetList) <>
                makeAdamsNodes [] (lub2TreeRep curLUB) newLUB placedTaxa bothLeafLists --inTreeVertexLists vertexLeafSetList)

-- | getLeafSetFromNodeName takes String name of node and returns sorted list of leaf
-- names--ASSUMES node names are not given in input and are assigned as trees
-- are parsed
getLeafSetFromNodeName :: Vertex -> [String]
getLeafSetFromNodeName inVertex =
    let (nodeName, _, _, _) = inVertex
    in
    if null nodeName then error "Null node name in getLeafSetFromNodeName"
    else
        let rawList = map (replaceChar ['(', ')', ','] ' ') nodeName
        in
        L.sort $ words rawList --this sort required

-- | lub2TreeRep takes list of names and makes into unresolved subtree
-- in parens
lub2TreeRep :: [String] -> String
lub2TreeRep inStringList
  | null inStringList = error "Null input in lub2TreeRep"
  | length inStringList == 1 = head inStringList
  | otherwise =
    let inside = init $ concatMap (<> ",") inStringList
    in
    ( '(' : inside ) <> ")"

-- | getDecendantLeafList iputs a vertex and returns leaf set (as list of
-- leaf names as strings) descdended from
-- that vertex, if a leaf, returns that leaf
getDecendantLeafList :: [Vertex] -> PhyloGraphVect -> [String]
getDecendantLeafList inVertexList inGraphVect =
    if null inVertexList then []
    else
        let (curVertName, descList, _, vertType) = head inVertexList
            descVertList = getVertsFromIndexList  descList inGraphVect
        in
        if vertType == Leaf then
            curVertName : getDecendantLeafList (tail inVertexList) inGraphVect
        else
            getDecendantLeafList [head descVertList] inGraphVect
                <> getDecendantLeafList (tail descVertList) inGraphVect
                <> getDecendantLeafList (tail inVertexList) inGraphVect

-- | getSplitLeafList takes a node and returns a list of list of descendent leaves
getSplitLeafList :: [Int] -> PhyloGraphVect -> [[String]]
getSplitLeafList descList inGraphVect =
    if null descList then []
    else
        let curDesc = head descList
            (vertexVect, _)  = inGraphVect
            curLeaves = getDecendantLeafList [vertexVect V.! curDesc] inGraphVect
        in curLeaves : getSplitLeafList (tail descList) inGraphVect

-- | getSplitLeafListList takes list of descenndents for PhyloGraphVect and
-- returns a list of descendant list for each split of each tree
getSplitLeafListList :: [[Int]] -> [PhyloGraphVect] -> [[[String]]]
getSplitLeafListList descListList inGraphVectList
  | null descListList = []
  | null inGraphVectList = error "Diff numbers of descdent lists and graphs"
  | otherwise =
    let curIntList = head descListList
        curGraphVectList = head inGraphVectList
    in
    getSplitLeafList curIntList curGraphVectList :
        getSplitLeafListList (tail descListList) (tail inGraphVectList)

-- | lub2 takes two lists of lists of names and generates the pairswise set of
-- intersections
lub2 :: [[String]] -> [[String]] -> [[String]]
lub2 s1 s2
  | null s1 = []
  | null s2 = []
  | otherwise =
    let intersectFirst = L.intersect (head s1) (head s2) : lub2 [head s1] (tail s2)
    in
    intersectFirst <> lub2 (tail s1) s2

-- | leastUpperBound takes list of list vertex leaf descendants (as Strings)
-- and returns LUB of Adams II (1972) consensus
leastUpperBound :: [[[String]]] -> [[String]]
leastUpperBound inVertexListList
  | length inVertexListList < 2 =
    error "Too few name lists in leastUpperBound"
  | length inVertexListList == 2 =
    let x = head inVertexListList
        y = last inVertexListList
    in
    lub2 x y
      | otherwise =
    let x = head inVertexListList
        y = head $ tail inVertexListList
        z = tail $ tail inVertexListList
        t = lub2 x y
    in
    leastUpperBound (t : z)

-- | get second retriueves 2nd element of 4
getSecond :: (a, b, c, d) -> b
getSecond inTuple  =
    let (_, b2, _, _) = inTuple
    in
    b2

-- | leafSetFromVertexVect takes vector of veritces and returns set of leaf
-- names
leafSetFromVertexVect :: Set.Set String -> V.Vector Vertex -> Set.Set String
leafSetFromVertexVect inSet inVerts =
    if V.null inVerts then inSet
    else
        let (curName, _, _, curType) = V.head inVerts
        in
        if curType == Leaf then
            leafSetFromVertexVect (Set.insert curName inSet) (V.tail inVerts)
        else
            leafSetFromVertexVect inSet (V.tail inVerts)

-- | getLeafSet tgake a list pf PhyloGraphVect and returns a pair with
-- True if the leaf sets are identical, and a list of the sets
getLeafSet :: PhyloGraphVect -> Set.Set String
getLeafSet inGraphVect =
    let (inVerts, _) = inGraphVect
    in
    leafSetFromVertexVect Set.empty inVerts


-- | setEqual checks for set equality by difference between union and
-- intersection is empty
setEqual :: Ord a => Set.Set a -> Set.Set a-> Bool
setEqual firstSet secondSet =
    let combinedElem = Set.union firstSet secondSet
        sameElem = Set.intersection firstSet secondSet
    in
    Set.empty == Set.difference combinedElem sameElem


-- | getAndCheckLeafSets take graphs and checks that leaf sets are identical
getAndCheckLeafSets :: [PhyloGraphVect] -> (Bool, [Set.Set String])
getAndCheckLeafSets inGraphs =
    if null inGraphs then error "Empty graph list in getAndCheckLeafSets"
    else
        let leafSetList = map getLeafSet inGraphs
            firstSet = head leafSetList
            setDiffList = map (setEqual firstSet) (tail leafSetList)
            allEmpty = and setDiffList
        in
        (allEmpty, leafSetList)

-- | findRoot take PhyloGraphVect and return root index and Vertex
findRoot :: Int -> PhyloGraphVect -> (Int, Vertex)
findRoot index inGraph  =
    let (vertexVect, _) = inGraph
    in
    if index < V.length vertexVect then
        let (_, _, _, vertexType) = vertexVect V.! index
        in
        if vertexType == Root then
            (index, vertexVect V.! index)
        else
            findRoot (index + 1) inGraph
     else error "Index exceeeds vertex number in findRoot"

