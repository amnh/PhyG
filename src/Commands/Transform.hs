{- |
Module      :  Transform.hs
Description :  Module to coordinate transform command execution
Copyright   :  (c) 2022 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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

module Commands.Transform
  ( transform
  ) where

import Types.Types
import qualified Search.Build as B
import qualified GraphOptimization.Traversals as T
import qualified Utilities.Utilities as U
import qualified ParallelUtilities            as PU
import Control.Parallel.Strategies
import           Data.Maybe
import           Text.Read
import           Data.Char
import qualified Graphs.GraphOperations  as GO
import GeneralUtilities
import qualified Data.List              as L
import qualified Utilities.LocalGraph    as LG
import qualified Data.Vector as V

import qualified SymMatrix                   as S


-- | transformArgList is the list of valid transform arguments
transformArgList :: [String]
transformArgList = ["totree", "tosoftwired", "tohardwired", "staticapprox", "dynamic", "atrandom", "first", "displaytrees"]


-- | transform changes aspects of data sande settings during execution
-- as opposed to Set with all happens at begginign of program execution
transform :: [Argument] -> GlobalSettings -> ProcessedData -> ProcessedData -> Int -> [PhylogeneticGraph] -> (GlobalSettings, ProcessedData, [PhylogeneticGraph])
transform inArgs inGS origData inData rSeed inGraphList = 
   let fstArgList = fmap (fmap toLower . fst) inArgs
       sndArgList = fmap (fmap toLower . snd) inArgs
       lcArgList = zip fstArgList sndArgList
       checkCommandList = U.checkCommandArgs "transform" fstArgList transformArgList
   in
   -- check for valid command options
   if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'transform': " ++ show inArgs)
   else 
        let displayBlock = filter ((=="displaytrees").fst) lcArgList
            numDisplayTrees
               | length displayBlock > 1 =
                  errorWithoutStackTrace ("Multiple displayTree number specifications in tansform--can have only one: " ++ show inArgs)
               | null displayBlock = Just 10
               | null (snd $ head displayBlock) = Just 10
               | otherwise = readMaybe (snd $ head displayBlock) :: Maybe Int

            toTree = any ((=="totree").fst) lcArgList
            toSoftWired = any ((=="tosoftwired").fst) lcArgList
            toHardWired = any ((=="tohardwired").fst) lcArgList
            toStaticApprox = any ((=="staticapprox").fst) lcArgList
            toDynamic = any ((=="todynamic").fst) lcArgList
            atRandom = any ((=="atrandom").fst) lcArgList
            chooseFirst = any ((=="first").fst) lcArgList


        in
        if (length $ filter (== True) [toTree, toSoftWired, toHardWired]) > 1 then 
            errorWithoutStackTrace ("Multiple graph transform commands--can only have one : " ++ (show inArgs))
        else if toStaticApprox && toDynamic then 
            errorWithoutStackTrace ("Multiple staticApprox/Dynamic transform commands--can only have one : " ++ (show inArgs))
        else if atRandom && chooseFirst then 
            errorWithoutStackTrace ("Multiple display tree choice commands in transform (first, atRandom)--can only have one : " ++ (show inArgs))
        else if (toTree || toSoftWired || toHardWired) && (toDynamic || toDynamic) then 
            errorWithoutStackTrace ("Multiple transform operations in transform (e.g. toTree, staticApprox)--can only have one at a time: " ++ (show inArgs))
        else
            let pruneEdges = False
                warnPruneEdges = False
                startVertex = Nothing
            in          
            -- transform nets to tree
            if toTree then
               -- already Tree return
               if (graphType inGS == Tree) then (inGS, inData, inGraphList)
               else 
                  let newGS = inGS {graphType = Tree}
                  
                      -- generate and return display trees-- displayTreNUm / graph
                      displayGraphList = if chooseFirst then fmap (take (fromJust numDisplayTrees) . GO.generateDisplayTrees) (fmap fst6 inGraphList)
                                         else fmap (GO.generateDisplayTreesRandom rSeed (fromJust numDisplayTrees)) (fmap fst6 inGraphList)

                      -- prob not required
                      displayGraphs = fmap GO.ladderizeGraph $ fmap GO.renameSimpleGraphNodes (concat displayGraphList)

                      -- reoptimize as Trees
                      newPhylogeneticGraphList = fmap (T.multiTraverseFullyLabelGraph newGS inData pruneEdges warnPruneEdges startVertex) displayGraphs `using` PU.myParListChunkRDS
                  in
                  (newGS, inData, newPhylogeneticGraphList)
            
            -- transform to softwired
            else if toSoftWired then 
               if (graphType inGS == SoftWired) then (inGS, inData, inGraphList)
               else 
                  let newGS = inGS {graphType = SoftWired}
                      newPhylogeneticGraphList = fmap (T.multiTraverseFullyLabelGraph newGS inData pruneEdges warnPruneEdges startVertex) (fmap fst6 inGraphList)  `using` PU.myParListChunkRDS
                  in
                  (newGS, inData, newPhylogeneticGraphList)

            -- transform to hardwired
            else if toHardWired then
               if (graphType inGS == HardWired) then (inGS, inData, inGraphList)
               else 
                  let newGS = inGS {graphType = HardWired}
                      
                      newPhylogeneticGraphList = fmap (T.multiTraverseFullyLabelGraph newGS inData pruneEdges warnPruneEdges startVertex) (fmap fst6 inGraphList)  `using` PU.myParListChunkRDS
                  in
                  (newGS, inData, newPhylogeneticGraphList)

            -- roll back to dynamic data from static approx      
            else if toDynamic then 
               if (origData /= inData) then
                  let newPhylogeneticGraphList = fmap (T.multiTraverseFullyLabelGraph inGS origData pruneEdges warnPruneEdges startVertex) (fmap fst6 inGraphList)  `using` PU.myParListChunkRDS
                  in
                  (inGS, origData, newPhylogeneticGraphList)

               else (inGS, inData, inGraphList)

            -- transform to static approx--using first Tree
            else if toStaticApprox then
               let newData = makeStaticApprox inGS inData (head $ L.sortOn snd6 inGraphList) 
                   newPhylogeneticGraphList = fmap (T.multiTraverseFullyLabelGraph inGS newData pruneEdges warnPruneEdges startVertex) (fmap fst6 inGraphList)  `using` PU.myParListChunkRDS

               in
               (inGS, newData, newPhylogeneticGraphList)


            else error ("Transform type not implemented/recognized" ++ (show inArgs))
            
  
-- | makeStaticApprox takes ProcessedData and returns static approx (implied alignment recoded) ProcessedData
-- if Tree take SA fields and recode appropriatrely given cost regeme of character
-- if Softwired--use display trees for SA
-- if hardWired--convert to softwired and use display trees for SA
-- since for heuristic searcing--uses additive weight for sequences and simple cost matrices, otherwise
-- matrix characters
makeStaticApprox :: GlobalSettings -> ProcessedData -> PhylogeneticGraph -> ProcessedData
makeStaticApprox inGS inData inGraph = 
   if LG.isEmpty (fst6 inGraph) then error "Empty graph in makeStaticApprox"

   -- tree type
   else if graphType inGS == Tree then
      let charInfoVV = six6 inGraph
          decGraph = thd6 inGraph
          (nameV, nameBVV, blockDataV) = inData

          -- do each block in turn pulling and transforming data from inGraph
          newBlockDataV = fmap (pullGraphBlockDataAndTransform decGraph  charInfoVV inData) [0..(length blockDataV - 1)] `using` PU.myParListChunkRDS
      in
      (nameV, nameBVV, V.fromList newBlockDataV)

   else error ("Static Approx not yet implemented for gaph type :" ++ (show $ graphType inGS))


-- | pullGraphBlockDataAndTransform takes a DecoratedGrpah and block index and pulls 
-- the character data of the block and transforms the leaf data by using implied alignment
-- feilds for dynamic characters
pullGraphBlockDataAndTransform :: DecoratedGraph -> V.Vector (V.Vector CharInfo) -> ProcessedData -> Int -> BlockData
pullGraphBlockDataAndTransform inDecGraph charInfoVV (blockName, _, blockCharInfoV) blockIndex =
   let (_, leafLabelList) = unzip $ LG.labNodes inDecGraph
       leafBlockData = fmap (V.! blockIndex) (fmap vertData leafLabelList)
       (transformedLeafBlockData, transformedBlockInfo) = unzip $ fmap (transformData (thd3 $ blockCharInfoV V.! blockIndex)) leafBlockData
   in
   (blockName V.! blockIndex, V.fromList transformedLeafBlockData, head transformedBlockInfo)


-- | transformData takes originalCharacter info and chracter data and transforms to static if dynamic noting chracter type
transformData :: V.Vector CharInfo -> V.Vector CharacterData -> (V.Vector CharacterData, V.Vector CharInfo) 
transformData inCharInfoV inCharDataV =
   if V.null inCharInfoV then 
      (V.empty, V.empty)
   else 
      let (inCharDataV, inCharInfoV) = V.unzip $ V.zipWith transformCharacter inCharDataV inCharInfoV 
      in
      (inCharDataV, inCharInfoV)

-- transformCharacter takes a single characer info and character and returns IA if statis as is if not
transformCharacter :: CharacterData -> CharInfo -> (CharacterData, CharInfo)
transformCharacter inCharData inCharInfo =
   let inCharType = charType inCharInfo
       inCostMatrix = costMatrix inCharInfo

       -- determine if matrix is all same costs => nonadditive
       --                        all same except fort single indel costs => non add with gap binary chars
       --                        not either => matrix char
       inCostMatrixType = getRecodingType inCostMatrix

   in
   if inCharType `elem` exactCharacterTypes then (inCharData, inCharInfo)
   else 
      -- different types--vector wrangling
      if inCharType `elem` [SlimSeq, NucSeq] then 
         (inCharData, inCharInfo)
      else if inCharType `elem` [WideSeq, AminoSeq] then
         (inCharData, inCharInfo)
      else if inCharType == HugeSeq then
         (inCharData, inCharInfo)
      else 
         error ("Unrecognized character type in transformCharacter: " ++ (show inCharType)) 
      

-- | getRecodingType takes a cost matrix and detemines if it can be recodes as non-additive, 
-- non-additive with gap chars, or matrix
-- assumes indel costs are in last row and column
getRecodingType :: S.Matrix Int -> String
getRecodingType inMatrix =
   if S.null inMatrix then error "Null matrix in getRecodingType"
   else
      if (not . S.isSymmetric) inMatrix then "matrix"
      else 
         let rows = S.rows inMatrix 
             matrixLL = S.toFullLists inMatrix
             lastRow = L.last matrixLL
             numUniqueCosts = length $ L.group $ (filter (/= 0) $ concat matrixLL) 

         in
         -- all same except for 0
         if numUniqueCosts == 1 then "nonAdd"

         -- all same except for gaps
         else if numUniqueCosts == 2 then
            if  (length $ L.group $ filter (/= 0) lastRow) == 1 then "nonAddGap"

            -- some no gaps different
            else "matrix"

         -- to many types for nonadd coding
         else "matrix"