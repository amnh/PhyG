{-# LANGUAGE BangPatterns #-}

{- |
Module to calculate distance tree construction methods Neightbor-Joining, WPGMA,
and WPGMAcbut with added refinement based on 4-point metric.
-}
module Search.DistanceMethods (neighborJoining, wPGMA, doWagnerS, performWagnerRefinement) where

import Control.Monad (when)
import Data.Number.Transfinite as NT
import Data.Vector qualified as V
import GeneralUtilities
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import PHANE.Evaluation.Verbosity (Verbosity (..))
import Search.DistanceWagner qualified as W
import SymMatrix qualified as M
import Types.DistanceTypes
import Types.Types
import Utilities.DistanceUtilities


{- | wPGMA takes a list of leaves and a distance matrixx and returns
an WGPGMA tree
WPGMA not UPGMA since the linkages are from the two descendent linkages and not weighted by
the number of taxa in each group as in https://en.wikipedia.org/wiki/WPGMA
-}
wPGMA ∷ V.Vector String → M.Matrix Double → Int → PhyG TreeWithData
wPGMA leafNames distMatrix outgroup =
    if M.null distMatrix
        then error "Null matrix in WPGMA"
        else
            let numLeaves = V.length leafNames
                leafVertexVect = V.fromList [0 .. (numLeaves - 1)]
            in  do
                    ((vertexVect, edgeVect), finalMatrix) ← addTaxaWPGMA distMatrix numLeaves (leafVertexVect, V.empty) []
                    let wPGMATree' = convertLinkagesToEdgeWeights vertexVect (V.toList edgeVect) (V.toList edgeVect) [] numLeaves
                    -- linkage leves are not same as edgeweights so need to be converted
                    let newickString = convertToNewick leafNames outgroup wPGMATree'
                    let treeCost = getTreeCost wPGMATree'

                    logWith LogInfo "\n\tBuilding WPGMA tree\n"
                    -- trace (show wPGMATree <> "\n" <> show wPGMATree')
                    return $ (take (length newickString - 3) newickString <> "[" <> show treeCost <> "];\n", wPGMATree', treeCost, finalMatrix)


-- | pulls dWagner function from module Wagner
doWagnerS
    ∷ GlobalSettings → V.Vector String → M.Matrix Double → String → Int → String → Int → [V.Vector Int] → PhyG [TreeWithData]
doWagnerS inGS leafNames distMatrix firstPairMethod outgroup addSequence numToKeep replicateSequences =
    do
        logWith LogInfo ("\tBuilding " <> show (length replicateSequences) <> " Wagner tree(s)" <> "\n")
        W.doWagnerS inGS leafNames distMatrix firstPairMethod outgroup addSequence numToKeep replicateSequences


-- | pulls Wagner refinement from Wagner module
performWagnerRefinement ∷ String → String → String → V.Vector String → Int → TreeWithData → PhyG [TreeWithData]
performWagnerRefinement = W.performRefinement


{- | neighborJoining takes a list of leaves and a distance matrixx and returns
an NJ tree
-}
neighborJoining ∷ V.Vector String → M.Matrix Double → Int → PhyG TreeWithData
neighborJoining leafNames distMatrix outgroup =
    if M.null distMatrix
        then error "Null matrix in neighborJoining"
        else do
            logWith LogInfo "\n\tBuilding NJ tree\n"
            -- get intial matrices
            let numLeaves = V.length leafNames
            let leafVertexVect = V.fromList [0 .. (numLeaves - 1)]
            (nJTree, finalLittleDMatrix) ← addTaxaNJ distMatrix numLeaves (leafVertexVect, V.empty) []
            let newickString = convertToNewick leafNames outgroup nJTree
            let treeCost = getTreeCost nJTree

            return $ (take (length newickString - 3) newickString <> "[" <> show treeCost <> "];\n", nJTree, treeCost, finalLittleDMatrix)


-- | sumAvail sums the row only thos values not already added to tree
sumAvail ∷ [Int] → Int → [Double] → Double
sumAvail vertInList index distList
    | null distList = 0.0
    | index `elem` vertInList = sumAvail vertInList (index + 1) (tail distList)
    | otherwise =
        let firstDist = head distList
        in  firstDist + sumAvail vertInList (index + 1) (tail distList)


-- | makeDMatrixRow make a single row of the bif D matrix
makeDMatrixRow ∷ M.Matrix Double → [Int] → Int → Int → V.Vector Double -- V.Seq Double
makeDMatrixRow inObsMatrix vertInList column row
    | M.null inObsMatrix = error "Null matrix in makeInitialDMatrix"
    | row `elem` vertInList = V.replicate (V.length (inObsMatrix V.! row)) NT.infinity
    | column == V.length (inObsMatrix V.! row) = V.empty
    | column == row = V.cons 0.0 (makeDMatrixRow inObsMatrix vertInList (column + 1) row)
    | column `notElem` vertInList =
        let dij = inObsMatrix M.! (row, column)
            divisor = (fromIntegral (M.rows inObsMatrix) - 2) - fromIntegral (length vertInList)
            ri = (sumAvail vertInList 0 $ M.getFullRow inObsMatrix row)
            rj = (sumAvail vertInList 0 $ M.getFullRow inObsMatrix column)
            bigDij = dij - ((ri + rj) / divisor)
        in  V.cons bigDij (makeDMatrixRow inObsMatrix vertInList (column + 1) row)
    | otherwise = V.cons NT.infinity (makeDMatrixRow inObsMatrix vertInList (column + 1) row)


{- | makeIDMatrix makes adjusted matrix (D) from observed (d) values
assumes matrix is square and symmetrical
makes values Infinity if already added
adjust ri and rj to bew based on on values not in termInList
does by row so can be parallelized call with column = 0 update list []
makes DMatrix direclty not via M.updateMatrix
-}
makeDMatrix' ∷ M.Matrix Double → [Int] → Int → Int → [(Int, Int, Double)] → M.Matrix Double
makeDMatrix' inObsMatrix vertInList row column updateList
    | M.null inObsMatrix = error "Null matrix in makeInitialDMatrix"
    | row == M.rows inObsMatrix = M.updateMatrix inObsMatrix updateList
    | column == M.cols inObsMatrix = makeDMatrix' inObsMatrix vertInList (row + 1) 0 updateList
    | column == row = makeDMatrix' inObsMatrix vertInList row (column + 1) ((row, column, 0.0) : updateList)
    | (column `elem` vertInList) || (row `elem` vertInList) =
        makeDMatrix' inObsMatrix vertInList row (column + 1) ((row, column, NT.infinity) : updateList)
    | otherwise =
        let dij = inObsMatrix M.! (row, column)
            divisor = (fromIntegral (M.rows inObsMatrix) - 2) - fromIntegral (length vertInList)
            ri = (sumAvail vertInList 0 $ M.getFullRow inObsMatrix row)
            rj = (sumAvail vertInList 0 $ M.getFullRow inObsMatrix column)
            bigDij = dij - ((ri + rj) / divisor)
        in  makeDMatrix' inObsMatrix vertInList row (column + 1) ((row, column, bigDij) : updateList)


{- | makeIDMatrix makes adjusted matrix (D) from observed (d) values
assumes matrix is square and symmetrical
makes values Infinity if already added
adjust ri and rj to be based on on values not in termInList
does by row so can be parallelized call with column = 0 update list []
makes DMatrix direclty not via M.updateMatrix
-}
makeDMatrix ∷ M.Matrix Double → [Int] → PhyG (M.Matrix Double)
makeDMatrix inObsMatrix vertInList =
    if M.null inObsMatrix
        then error "Null matrix in makeInitialDMatrix"
        else
            let makeRowAction ∷ Int → V.Vector Double
                makeRowAction = makeDMatrixRow inObsMatrix vertInList 0
            in  do
                    makeDRowPar ← getParallelChunkMap
                    let newMatrix = makeDRowPar makeRowAction [0 .. (M.rows inObsMatrix - 1)]

                    pure $ V.fromList newMatrix


{- | pickNearestUpdateMatrix takes d and D matrices, pickes nearest based on D
then updates d and D to reflect new node and distances created
updates the column/row for vertices that are joined to be infinity so
won't be chosen to join again
-}
pickNearestUpdateMatrixNJ ∷ M.Matrix Double → [Int] → PhyG (M.Matrix Double, Vertex, Edge, Edge, [Int])
pickNearestUpdateMatrixNJ littleDMatrix vertInList =
    if M.null littleDMatrix
        then error "Null d matrix in pickNearestUpdateMatrix"
        else do
            newMatrix ← makeDMatrix littleDMatrix vertInList
            minPairResult ← getMatrixMinPairTabu newMatrix vertInList
            let (iMin, jMin, distIJ) = minPairResult
            -- (iMin, jMin, distIJ) = getMatrixMinPairTabu (makeDMatrix' littleDMatrix vertInList 0 0 []) vertInList

            -- trace ("First pair " <> show (iMin, jMin, distIJ) <> " matrix: " <> (show (makeDMatrix littleDMatrix vertInList))
            --  <> "\nVertsIn: " <> show vertInList <> " matrix': " <> show (makeDMatrix' littleDMatrix vertInList 0 0 [])) (
            if distIJ == NT.infinity
                then error "No minimum found in pickNearestUpdateMatrix"
                else
                    let -- new vertex is size of distance matrix (0 indexed)
                        newVertIndex = M.rows littleDMatrix
                        dij = littleDMatrix M.! (iMin, jMin)
                        divisor = fromIntegral (M.rows littleDMatrix) - 2 - fromIntegral (length vertInList)
                        -- only count values of those in
                        ri = (sumAvail vertInList 0 $ M.getFullRow littleDMatrix iMin)
                        rj = (sumAvail vertInList 0 $ M.getFullRow littleDMatrix jMin)
                        -- seem reversed compared to examples, seems arbitrary to me (for leaf pairs at least)
                        -- diMinNewVert = (dij / 2.0) - ((ri - rj) / (2.0 * divisor))
                        -- djMinNewVert = dij - diMinNewVert
                        djMinNewVert = (dij / 2.0) - ((ri - rj) / (2.0 * divisor))
                        diMinNewVert = dij - djMinNewVert

                        newVertInList = vertInList <> [iMin, jMin]

                        -- get distances to existing vertices
                        otherVertList = [0 .. (M.rows littleDMatrix - 1)]

                        getNewDistAction ∷ Int → Double
                        getNewDistAction = getNewDist littleDMatrix dij iMin jMin diMinNewVert djMinNewVert
                    in  do
                            -- NOTE: This may be a problem for parallel calls but looks OK
                            newDistPar ← getParallelChunkMap
                            let newLittleDRow = newDistPar getNewDistAction otherVertList
                            -- newLittleDRow = PU.seqParMap PU.myStrategyR0 (getNewDist littleDMatrix dij iMin jMin diMinNewVert djMinNewVert) otherVertList -- `using` myParListChunkRDS
                            let newLittleDMatrix = M.addMatrixRow littleDMatrix (V.fromList $ newLittleDRow <> [0.0])
                            -- recalculate whole D matrix since new row affects all the original ones  (except those merged)
                            -- included vertex values set to infinity so won't be chosen later
                            -- newBigDMatrix = makeDMatrix newLittleDMatrix newVertInList -- 0 0 []

                            -- create new edges
                            let newEdgeI = (newVertIndex, iMin, diMinNewVert)
                            let newEdgeJ = (newVertIndex, jMin, djMinNewVert)

                            pure (newLittleDMatrix, newVertIndex, newEdgeI, newEdgeJ, newVertInList)


-- )

-- | getNewDist get ditance of new vertex to existing vertices
getNewDist ∷ M.Matrix Double → Double → Int → Int → Double → Double → Int → Double
getNewDist littleDMatrix dij iMin jMin diMinNewVert djMinNewVert otherVert
    | otherVert == iMin = diMinNewVert
    | otherVert == jMin = djMinNewVert
    | otherwise =
        let dik = littleDMatrix M.! (iMin, otherVert)
            djk = littleDMatrix M.! (jMin, otherVert)
        in  (dik + djk - dij) / 2.0


{- | addTaxaNJ recursively calls pickNearestUpdateMatrix untill all internal nodes are created
recursively called until all (n - 2) internal vertices are created.
-}
addTaxaNJ ∷ M.Matrix Double → Int → Tree → [Int] → PhyG (Tree, M.Matrix Double)
addTaxaNJ littleDMatrix numLeaves (vertexVect, edgeVect) vertInList =
    if V.length vertexVect == (2 * numLeaves) - 2
        then
            let -- (iMin, jMin, _) = getMatrixMinPairTabu  (makeDMatrix littleDMatrix vertInList) vertInList
                last2 = subtractVector (V.fromList vertInList) vertexVect
                iMin = last2 V.! 0
                jMin = last2 V.! 1
                {-This is wrong--not last two
                iMin = vertexVect V.! ((V.length vertexVect) - 1)
                jMin = vertexVect V.! ((V.length vertexVect) - 2)
                -}
                lastEdge = (iMin, jMin, littleDMatrix M.! (iMin, jMin))
            in  {-
                trace (show last2 <> " from " <> show vertexVect <> "\nlast edge: " <> " size " <> (show $ V.length vertexVect) <> " matrix: " <> (show littleDMatrix)
                  <> " edge: "
                  <> (show lastEdge) <> " vertInList: " <> show vertInList)
                  -}
                return ((vertexVect, edgeVect `V.snoc` lastEdge), littleDMatrix)
        else -- more to add
        do
            pickNearestResult ← pickNearestUpdateMatrixNJ littleDMatrix vertInList
            let (newLittleDMatrix, newVertIndex, newEdgeI, newEdgeJ, newVertInList) = pickNearestResult
            let newVertexVect = vertexVect `V.snoc` newVertIndex
            let newEdgeVect = edgeVect <> V.fromList [newEdgeI, newEdgeJ]

            -- trace (M.showMatrixNicely newLittleDMatrix <> "\n" <> M.showMatrixNicely bigDMatrix)
            let -- progress = takeWhile (/='.') $ show ((fromIntegral (100 * (V.length vertexVect - numLeaves))/fromIntegral (numLeaves - 2)) :: Double)
            let (percentAdded, _) = divMod (100 * (V.length vertexVect - numLeaves)) (numLeaves - 2)
            let (decileNumber, decileRemainder) = divMod percentAdded 10
            let (_, oddRemainder) = divMod (V.length vertexVect - numLeaves) 2

            if decileRemainder == 0 && oddRemainder == 0
                then do
                    logWith LogInfo ("\t" <> show (10 * decileNumber) <> "%")
                    addTaxaNJ newLittleDMatrix numLeaves (newVertexVect, newEdgeVect) newVertInList
                else addTaxaNJ newLittleDMatrix numLeaves (newVertexVect, newEdgeVect) newVertInList


-- | addTaxaWPGMA perfomrs recursive reduction of distance matrix until all internal vertices are created
addTaxaWPGMA ∷ M.Matrix Double → Int → Tree → [Int] → PhyG (Tree, M.Matrix Double)
addTaxaWPGMA distMatrix numLeaves (vertexVect, edgeVect) vertInList =
    if V.length vertexVect == (2 * numLeaves) - 2
        then
            let -- (iMin, jMin, _) = getMatrixMinPairTabu distMatrix vertInList -- (-1, -1, NT.infinity) 0 0
                last2 = subtractVector (V.fromList vertInList) vertexVect
                iMin = last2 V.! 0
                jMin = last2 V.! 1
                lastEdge = (iMin, jMin, distMatrix M.! (iMin, jMin))
            in  {-
                trace ((show last2) <> " from " <> show vertexVect <> "\nlast edge: " <> " size " <> (show $ V.length vertexVect) <> " matrix: " <> (show distMatrix)
                  <> " last edge: "
                  <> (show lastEdge) <> " vertInList: " <> (show vertInList)
                  <> " all edges " <> show (edgeVect `V.snoc` lastEdge))
                  -}
                pure ((vertexVect, edgeVect `V.snoc` lastEdge), distMatrix)
        else do
            -- building
            upfdateResult ← pickUpdateMatrixWPGMA distMatrix vertInList
            let (newDistMatrix, newVertIndex, newEdgeI, newEdgeJ, newVertInList) = upfdateResult
            let newVertexVect = vertexVect `V.snoc` newVertIndex
            let newEdgeVect = edgeVect <> V.fromList [newEdgeI, newEdgeJ]

            -- trace (M.showMatrixNicely distMatrix)
            let -- progress = takeWhile (/='.') $ show ((fromIntegral (100 * (V.length vertexVect - numLeaves))/fromIntegral (numLeaves - 2)) :: Double)
            let (percentAdded, _) = divMod (100 * (V.length vertexVect - numLeaves)) (numLeaves - 2)
            let (decileNumber, decileRemainder) = divMod percentAdded 10
            let (_, oddRemainder) = divMod (V.length vertexVect - numLeaves) 2

            if decileRemainder == 0 && oddRemainder == 0
                then do
                    logWith LogInfo ("\t" <> show (10 * decileNumber) <> "%")
                    addTaxaWPGMA newDistMatrix numLeaves (newVertexVect, newEdgeVect) newVertInList
                else do
                    addTaxaWPGMA newDistMatrix numLeaves (newVertexVect, newEdgeVect) newVertInList


{- | pickUpdateMatrixWPGMA takes d matrix, pickes closesst based on d
then updates d  to reflect new node and distances created
updates the column/row for vertices that are joined to be infinity so
won't be chosen to join again
-}
pickUpdateMatrixWPGMA ∷ M.Matrix Double → [Int] → PhyG (M.Matrix Double, Vertex, Edge, Edge, [Int])
pickUpdateMatrixWPGMA distMatrix vertInList =
    if M.null distMatrix
        then error "Null d matrix in pickNearestUpdateMatrix"
        else do
            minPairResult ← getMatrixMinPairTabu distMatrix vertInList -- (-1, -1, NT.infinity) 0 0
            let (iMin, jMin, dij) = minPairResult

            -- trace ("First pair " <> show (iMin, jMin, dij)) (
            if dij == NT.infinity
                then error "No minimum found in pickNearestUpdateMatrix"
                else
                    let -- new vertex is size of distance matrix (0 indexed)
                        newVertIndex = M.rows distMatrix

                        diMinNewVert = dij / 2.0
                        djMinNewVert = dij / 2.0

                        newVertInList = vertInList <> [iMin, jMin]

                        -- get distances to existing vertices
                        otherVertList = [0 .. (M.rows distMatrix - 1)]

                        wpgmaAction ∷ Int → Double
                        wpgmaAction = getNewDistWPGMA distMatrix iMin jMin diMinNewVert djMinNewVert
                    in  do
                            wpgmaPar ← getParallelChunkMap
                            let newDistRow = wpgmaPar wpgmaAction otherVertList
                            -- newDistRow = PU.seqParMap PU.myStrategyR0 (getNewDistWPGMA distMatrix iMin jMin diMinNewVert djMinNewVert) otherVertList -- `using` myParListChunkRDS
                            let newDistMatrix = M.addMatrixRow distMatrix (V.fromList $ newDistRow <> [0.0])

                            -- create new edges
                            let newEdgeI = (newVertIndex, iMin, diMinNewVert)
                            let newEdgeJ = (newVertIndex, jMin, djMinNewVert)

                            pure (newDistMatrix, newVertIndex, newEdgeI, newEdgeJ, newVertInList)


-- | getNewDistWPGMA get ditance of new vertex to existing vertices WPGMA--cluster levels
getNewDistWPGMA ∷ M.Matrix Double → Int → Int → Double → Double → Int → Double
getNewDistWPGMA distMatrix iMin jMin diMinNewVert djMinNewVert otherVert
    | otherVert == iMin = diMinNewVert
    | otherVert == jMin = djMinNewVert
    | otherwise =
        let dik = distMatrix M.! (iMin, otherVert)
            djk = distMatrix M.! (jMin, otherVert)
        in  (dik + djk) / 2.0


{- | convertLinkagesToEdgeWeights converts linake leevs to branch lengths
by subtracting descnedent linakge from edge weight
edges are created in linakege order so can use that direction, leaves are
always 2nd element in edge
-}
convertLinkagesToEdgeWeights ∷ V.Vector Vertex → [Edge] → [Edge] → [Edge] → Int → Tree
convertLinkagesToEdgeWeights vertexVect fullEdgeList inEdgeList curEdgeList numLeaves
    | null fullEdgeList = error "Null edge set in convertLinkagesToEdgeWeights"
    | V.null vertexVect = error "Null vertex set in convertLinkagesToEdgeWeights"
    | null inEdgeList = (vertexVect, V.fromList curEdgeList)
    | null curEdgeList =
        -- first case, take last edge (highest linkage) special case need
        -- to subtract both descendent linkages
        let (eVert, uVert, linkage) = last inEdgeList
            eDescLinkage = if eVert >= numLeaves then getWeightDescLink eVert fullEdgeList else 0.0
            uDescLinkage = if uVert >= numLeaves then getWeightDescLink uVert fullEdgeList else 0.0
            newEdgeList = [(eVert, uVert, linkage - eDescLinkage - uDescLinkage)]
        in  convertLinkagesToEdgeWeights vertexVect fullEdgeList (init inEdgeList) newEdgeList numLeaves
    | otherwise -- not first but still some edegs to go
        =
        let (eVert, uVert, linkage) = head inEdgeList
            uDescLinkage = if uVert >= numLeaves then getWeightDescLink uVert fullEdgeList else 0.0
            newEdgeList = (eVert, uVert, linkage - uDescLinkage) : curEdgeList
        in  convertLinkagesToEdgeWeights vertexVect fullEdgeList (tail inEdgeList) newEdgeList numLeaves


{- | getWeightDescLink takes the m,pore derived (smaller linkage, 2nd vertex of edge) and returns
the weight of that edge (assumes not leaf--checked earlier)
-}
getWeightDescLink ∷ Int → [Edge] → Double
getWeightDescLink uVert fullEdgeList =
    if null fullEdgeList
        then error "Edge not found in getWeightDescLink"
        else
            let (eVert, _, weight) = head fullEdgeList
            in  if eVert == uVert
                    then weight
                    else getWeightDescLink uVert (tail fullEdgeList)
