{- |
Module specifying data types
-}

module Utilities.Distances
  ( getPairwiseDistances
  , getBlockDistance
  , getPairwiseBlockDistance
  ) where

import Control.Parallel.Strategies
import Data.Vector qualified as V
import Debug.Trace
import GeneralUtilities
import GraphOptimization.Medians qualified as M
import ParallelUtilities qualified as P
import SymMatrix qualified as S
import Types.Types
import Data.List qualified as L
import Utilities.Utilities qualified as U

{- |
getPairwiseDistances takes Processed data
and retuns a matrix (list of lists of Double) of pairwise
distances among vertices in data set over blocks ans all character types
sums over blocks
-}
getPairwiseDistances :: ProcessedData ->  [[VertexCost]]
getPairwiseDistances (nameVect, _, blockDataVect)
  | V.null nameVect = error "Null name vector in getPairwiseDistances"
  | V.null blockDataVect = error "Null Block Data vector in getPairwiseDistances"
  | otherwise =
    let -- get maximum observations and pairwise max observations to normalize distances
        maxDistance = U.getMaxNumberObservations blockDataVect
        pairList = makeIndexPairs True (V.length nameVect) (V.length nameVect) 0 0
        pairListCosts = fmap (U.getPairwiseObservations blockDataVect) pairList `using` P.myParListChunkRDS
        normFactorList = fmap (maxDistance /) $ fmap (max 1.0) pairListCosts
        initialFactorMatrix = S.fromLists $ replicate (V.length nameVect) $ replicate (V.length nameVect) 0.0
        (iLst, jList) = unzip pairList
        threeList = zip3 iLst jList normFactorList
        factorMatrix = S.updateMatrix initialFactorMatrix threeList


        -- get pairwise distances
        blockDistancesList =  V.toList $ V.map (getPairwiseBlockDistance (V.length nameVect)) blockDataVect
        summedBlock = L.foldl1' (S.zipWith (+)) blockDistancesList

        -- rescaled pairwsie distances 
        rescaledDistanceMatrix = S.zipWith (*) factorMatrix summedBlock
    in
    -- trace ("Factor:" <> show maxDistance <> " : "  <> (show normFactorList))
    trace ("\tGenerating pairwise distances for " <> show (V.length blockDataVect) <> " character blocks") 
    S.toFullLists rescaledDistanceMatrix -- summedBlock



-- | getBlockDistance takes Block data and returns distance between
-- vertices based on block data
-- this can be done for leaves only or all via the input processed
-- data leaves are first--then HTUs follow
getBlockDistance :: BlockData -> (Int, Int) -> VertexCost
getBlockDistance (_, localVertData, blockCharInfo) (firstIndex, secondIndex)  =
    let pairCost = V.sum $ V.map snd $ M.median2 (localVertData V.! firstIndex) (localVertData V.! secondIndex) blockCharInfo
    in
    pairCost

-- | getPairwiseBlocDistance returns pairwsie distances among vertices for
-- a block of data
-- this can be done for ;leaves only or all via the input processed
-- data leaves are first--then HTUs follow
getPairwiseBlockDistance :: Int -> BlockData-> S.Matrix VertexCost
getPairwiseBlockDistance numVerts inData  =
    let pairList = makeIndexPairs True numVerts numVerts 0 0
        initialPairMatrix = S.fromLists $ replicate numVerts $ replicate numVerts 0.0
        pairListCosts = fmap (getBlockDistance inData) pairList `using` P.myParListChunkRDS
        (iLst, jList) = unzip pairList
        threeList = zip3 iLst jList pairListCosts
        newMatrix = S.updateMatrix initialPairMatrix threeList
    in
    --trace ("NM:\n" <> (show threeList) <> "\n" <>(show $ S.toFullLists newMatrix))
    newMatrix

