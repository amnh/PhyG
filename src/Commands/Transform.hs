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
            -- transform nets to tree
            if toTree then
               -- already Tree return
               if (graphType inGS == Tree) then (inGS, inData, inGraphList)
               else 
                  let newGS = inGS {graphType == Tree}
                      pruneEdges = False
                      warnPruneEdges = False
                      startVertex = Nothing

                      -- generate and return display trees-- displayTreNUm / graph
                      displayGraphList = if chooseFirst then fmap (take (fromJust numDisplayTrees) $ GO.generateDisplayTrees) (fmap fst6 inGraphList)
                                         else GO.generateDisplayTreesRandom rSeed (fromJust numDisplayTrees) (fmap fst6 inGraphList)

                      -- prob not required
                      displayGraphs = fmap GO.ladderizeGraph $ fmap GO.renameSimpleGraphNodes (concat displayGraphList)

                      -- reoptimize as Trees
                      newPhylogeneticGraphList  = fmap (T.multiTraverseFullyLabelGraph newGS inData pruneEdges warnPruneEdges startVertex) displayGraphs `using` PU.myParListChunkRDS
                  in
                  (newGS, inData, newPhylogeneticGraphList)
            
            else if toSoftWired then 
               let newGS = inGS {graphType == SoftWired}
                   pruneEdges = False
                   warnPruneEdges = False
                   startVertex = Nothing
                   newPhylogeneticGraphList  = fmap (T.multiTraverseFullyLabelGraph newGS inData pruneEdges warnPruneEdges startVertex) (fmap fst6 inGraphList)  `using` PU.myParListChunkRDS
               in
               (newGS, inData, newPhylogeneticGraphList)

            else if toHardWired then
               let newGS = inGS {graphType == HardWired}
                   pruneEdges = False
                   warnPruneEdges = False
                   startVertex = Nothing
                   newPhylogeneticGraphList  = fmap (T.multiTraverseFullyLabelGraph newGS inData pruneEdges warnPruneEdges startVertex) (fmap fst6 inGraphList)  `using` PU.myParListChunkRDS
               in
               (newGS, inData, newPhylogeneticGraphList)


            else error "Graph transform no implemenmted/recognized" ++ (show inArgs)
            
  