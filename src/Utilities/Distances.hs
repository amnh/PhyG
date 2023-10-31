{- |
Module specifying data types
-}

module Utilities.Distances
  ( getPairwiseDistances
  , getBlockDistance
  , getPairwiseBlockDistance
  ) where

import Control.Evaluation
import Control.Monad.Logger (LogLevel (..), Logger (..), Verbosity (..))
import Data.Vector qualified as V
import GeneralUtilities
import GraphOptimization.Medians qualified as M
import ParallelUtilities qualified as P
import SymMatrix qualified as S
import System.ErrorPhase (ErrorPhase (..))
import Types.Types
import Data.List qualified as L
import Utilities.Utilities qualified as U
-- import Debug.Trace

{- |
getPairwiseDistances takes Processed data
and retuns a matrix (list of lists of Double) of pairwise
distances among vertices in data set over blocks ans all character types
sums over blocks
-}

getPairwiseDistances :: ProcessedData ->  PhyG [[VertexCost]]
getPairwiseDistances (nameVect, _, blockDataVect)
  | V.null nameVect = error "Null name vector in getPairwiseDistances"
  | V.null blockDataVect = error "Null Block Data vector in getPairwiseDistances"
  | otherwise =
    let -- get maximum observations and pairwise max observations to normalize distances
        maxDistance = U.getMaxNumberObservations blockDataVect
        pairList = makeIndexPairs True (V.length nameVect) (V.length nameVect) 0 0
        --pairListCosts = fmap (U.getPairwiseObservations blockDataVect) pairList `using` P.myParListChunkRDS
        --TODO
        -- pairListCosts = P.seqParMap rdeepseq   (U.getPairwiseObservations blockDataVect) pairList 
        --pairListCosts = fmap  (U.getPairwiseObservations blockDataVect) pairList 

        -- parallel setup
        action ::  (Int, Int) -> VertexCost
        action = U.getPairwiseObservations blockDataVect

    in do 
        pTraverse <- getParallelChunkMap
        let pairListCosts = pTraverse action pairList
        

        let normFactorList = fmap (maxDistance /) $ fmap (max 1.0) pairListCosts
        let initialFactorMatrix = S.fromLists $ replicate (V.length nameVect) $ replicate (V.length nameVect) 0.0
        let (iLst, jList) = unzip pairList
        let threeList = zip3 iLst jList normFactorList
        let factorMatrix = S.updateMatrix initialFactorMatrix threeList


        -- get pairwise distances
        -- let blockDistancesList =  V.toList $ V.map (getPairwiseBlockDistance (V.length nameVect)) blockDataVect
        blockDistancesList' <- mapM (getPairwiseBlockDistance (V.length nameVect)) blockDataVect
        let blockDistancesList =  V.toList blockDistancesList'

        let summedBlock = L.foldl1' (S.zipWith (+)) blockDistancesList

        -- rescaled pairwsie distances 
        let rescaledDistanceMatrix = S.zipWith (*) factorMatrix summedBlock
    
        -- trace ("Factor:" <> show maxDistance <> " : "  <> (show normFactorList))
        logWith LogInfo  ("\tGenerating pairwise distances for " <> show (V.length blockDataVect) <> " character blocks\n") 
        pure $ S.toFullLists rescaledDistanceMatrix -- summedBlock



-- | getBlockDistance takes Block data and returns distance between
-- vertices based on block data
-- this can be done for leaves only or all via the input processed
-- data leaves are first--then HTUs follow
getBlockDistance :: BlockData -> (Int, Int) -> VertexCost
getBlockDistance (_, localVertData, blockCharInfo) (firstIndex, secondIndex)  =
    let pairCost = V.sum $ V.map snd $ M.median2 (localVertData V.! firstIndex) (localVertData V.! secondIndex) blockCharInfo
    in
    pairCost

-- | getPairwiseBlocDistance returns pairwisee distances among vertices for
-- a block of data
-- this can be done for ;leaves only or all via the input processed
-- data leaves are first--then HTUs follow
getPairwiseBlockDistance :: Int -> BlockData-> PhyG (S.Matrix VertexCost)
getPairwiseBlockDistance numVerts inData  =
    let pairList = makeIndexPairs True numVerts numVerts 0 0
        initialPairMatrix = S.fromLists $ replicate numVerts $ replicate numVerts 0.0
        --pairListCosts = fmap (getBlockDistance inData) pairList `using` P.myParListChunkRDS
        --pairListCosts = P.seqParMap rdeepseq  (getBlockDistance inData) pairList 
        -- TODO 
        -- pairListCosts = fmap  (getBlockDistance inData) pairList 
        action ::  (Int, Int) -> VertexCost
        action = getBlockDistance inData
    in do
        pTraverse <- getParallelChunkMap
        let pairListCosts = pTraverse action pairList
        
        let (iLst, jList) = unzip pairList
        let threeList = zip3 iLst jList pairListCosts
        let newMatrix = S.updateMatrix initialPairMatrix threeList
    
        --trace ("NM:\n" <> (show threeList) <> "\n" <>(show $ S.toFullLists newMatrix))
        pure newMatrix

