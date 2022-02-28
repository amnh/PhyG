{- |
Module      :  Support.hs
Description :  Module containing support functions
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

module Support.Support  ( supportGraph
                        ) where

import Types.Types
import qualified ParallelUtilities       as PU
import Control.Parallel.Strategies
import Debug.Trace
import GeneralUtilities
import qualified Graphs.GraphOperations  as GO
import Utilities.Utilities               as U
import Data.Maybe
import           Data.Char
import           Text.Read
import qualified Data.Vector as V
import qualified Data.List as L
import qualified Search.Build as B
import qualified Search.Refinement as R
import qualified Utilities.Distances     as DD
import qualified Reconciliation.ReconcileGraphs as REC
import qualified Utilities.LocalGraph      as LG

-- | refinement arguments
supportArgList :: [String]
supportArgList = ["jackknife", "goodmanbremer", "gb", "gbsample", "replicates", "buildonly"] -- "bootstrap", 

-- | driver for overall support
supportGraph :: [Argument] -> GlobalSettings -> ProcessedData -> Int -> [PhylogeneticGraph] -> [PhylogeneticGraph]
supportGraph inArgs inGS inData rSeed inGraphList = 
   if null inGraphList then error ("No graphs input to calculate support") 
   else 
      let fstArgList = fmap (fmap toLower . fst) inArgs
          sndArgList = fmap (fmap toLower . snd) inArgs
          lcArgList = zip fstArgList sndArgList
          checkCommandList = U.checkCommandArgs "support" fstArgList supportArgList
     in
     -- check for valid command options
     if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'support': " ++ show inArgs)
     else 
         let doBootStrap = any ((=="bootstrap").fst) lcArgList
             onlyBuild = any ((=="buildonly").fst) lcArgList
             
             jackList   = filter ((=="jackknife").fst) lcArgList 
             jackFreq'   
              | length jackList > 1 =
                errorWithoutStackTrace ("Multiple jackknife sampling frequency specifications in support command--can have only one (e.g. jackknife:0.62): " ++ show inArgs)
              | null jackList = Just 0.6321 -- 1- 1/e
              | null (snd $ head jackList) = Just 0.6321 
              | otherwise = readMaybe (snd $ head jackList) :: Maybe Double

             replicateList   = filter ((=="replicates").fst) lcArgList 
             replicates'   
              | length replicateList > 1 =
                errorWithoutStackTrace ("Multiple resampling replicate specifications in support command--can have only one (e.g. replicates:100): " ++ show inArgs)
              | null replicateList = Just 100
              | otherwise = readMaybe (snd $ head replicateList) :: Maybe Int

             goodBremList   = filter ((`elem` ["goodmanbremer", "gb"]).fst) lcArgList 
             goodBremMethod   
              | length goodBremList > 1 =
                errorWithoutStackTrace ("Multiple Goodman-Bremer method specifications in support command--can have only one (e.g. gb:tbr): " ++ show inArgs)
              | null goodBremList = Just "tbr" 
              | otherwise = readMaybe (snd $ head goodBremList) :: Maybe String

             goodBremSampleList   = filter ((`elem` ["gbsample"]).fst) lcArgList 
             goodBremSample   
              | length goodBremSampleList > 1 =
                errorWithoutStackTrace ("Multiple Goodman-Bremer sample specifications in support command--can have only one (e.g. gbsample:1000): " ++ show inArgs)
              | null goodBremSampleList = Just (maxBound :: Int) 
              | otherwise = readMaybe (snd $ head goodBremSampleList) :: Maybe Int
  
         in
         if isNothing jackFreq' then errorWithoutStackTrace ("Jacknife frequency not a float (e.g. jackknife:0.5) in support: " ++ show (snd $ head jackList))
         else if isNothing replicates' then errorWithoutStackTrace ("Resampling replicates specification not a string (e.g. replicates:100) in support: " ++ show (snd $ head replicateList))
         else if isNothing goodBremMethod then errorWithoutStackTrace ("Goodman-Bremer method specification not a string (e.g. goodmanBremer:SPR) in support: " ++ show (snd $ head goodBremList))
         else if isNothing goodBremSample then errorWithoutStackTrace ("Goodman-Bremer sample specification not an integer (e.g. gbsample:1000) in support: " ++ show (snd $ head goodBremSampleList))
         else 
            let method = if doBootStrap && (not . null) jackList && (null goodBremList) then trace ("Bootstrap and Jackknife specified--defaulting to Jackknife") "jackknife"
                         else if (doBootStrap || (not . null) jackList) && (not . null) goodBremList then trace ("Resampling (Bootstrap or Jackknife) and Goodman-Bremer specified--defaulting to Goodman-Bremer") "goodBrem"
                         else if doBootStrap then 
                           trace ("Bootstrap not currently implemented--defaulting to Jackknife") 
                           -- "bootstrap"
                           "jackknife"
                         else if (not . null) jackList then "jackknife"
                         else "goodBrem"
                
                gbSampleSize = if goodBremSample == Just (maxBound :: Int)  then Nothing
                               else goodBremSample
                
                replicates = if fromJust replicates' < 0 then 
                                 trace ("Negative replicates number--defaulting to 100")
                                 100
                             else fromJust replicates'
                jackFreq = if fromJust jackFreq' <= 0 || fromJust jackFreq' >= 1.0 then
                              trace ("Jackknife frequency must be on (0.0, 1.0) defaulting to 0.6321")
                              0.6321
                           else fromJust jackFreq'

                buildOptions = [("distance",""), ("replicates", show 100), ("best", show 1), ("rdwag", ""), ("dWag", "")]
                swapOptions = if onlyBuild then []
                              else [("tbr", ""), ("steepest", ""), ("keep", show 1)]
                supportGraph = if method == "bootstrap" || method == "jackknife" then getResampleGraph inGS inData rSeed method replicates buildOptions swapOptions jackFreq inGraphList
                                   else getGoodBremGraphs inGS inData rSeed swapOptions gbSampleSize inGraphList
            in
            [supportGraph]
     
-- | getResampledGraphs performs resampling and search for bootstrap and jackknife support
getResampleGraph :: GlobalSettings -> ProcessedData -> Int -> String -> Int -> [(String, String)] -> [(String, String)] -> Double -> [PhylogeneticGraph] -> PhylogeneticGraph
getResampleGraph inGS inData rSeed resampleType replicates buildOptions swapOptions jackFreq inGraphList = 
   let resampledGraphList = fmap (makeResampledDataAndGraph inGS inData resampleType buildOptions swapOptions jackFreq) (take replicates $ randomIntList rSeed) `using` PU.myParListChunkRDS
       -- create appropriate support graph >50% ?
       -- need to add args
       reconcileArgs = if graphType inGS == Tree then [("method","majority"), ("compare","identity"), ("edgelabel","true"), ("vertexlabel","true"), ("connect","true"), ("threshold","51"), ("outformat", "dot")]
                       else [("method","eun"), ("compare","identity"), ("edgelabel","true"),  ("vertexlabel","true"), ("connect","true"), ("threshold","51"),("outformat", "dot")]
         -- majority ruke consensus if no args
       (reconciledGraphString, reconciledGraph) = REC.makeReconcileGraph REC.reconcileCommandList reconcileArgs (fmap fst6 resampledGraphList)
   in
   -- trace ("GRG: \n" ++ reconciledGraphString) (
   -- generate resampled graph
   if null inGraphList then (reconciledGraph, infinity, LG.empty, V.empty, V.empty, V.empty)

   -- label in graph with edge frequencies of resampled graph
   else (reconciledGraph, infinity, LG.empty, V.empty, V.empty, V.empty)
   -- )

-- | makeResampledDataAndGraph takes paramters, resmaples data and find a graph based on search parameters
-- returning the resampled graph
makeResampledDataAndGraph :: GlobalSettings -> ProcessedData -> String -> [(String, String)] -> [(String, String)] -> Double -> Int -> PhylogeneticGraph
makeResampledDataAndGraph inGS inData resampleType buildOptions swapOptions jackFreq rSeed =
   let randomIntegerList = randomIntList rSeed
       -- create resampled data
       newData = resampleData (randomIntegerList !! 0) resampleType jackFreq inData

       -- pairwise distances for distance analysis
       pairwiseDistances = DD.getPairwiseDistances newData

       -- build graphs
       buildGraphs = B.buildGraph buildOptions inGS newData pairwiseDistances (randomIntegerList !! 1)
       bestBuildGraphList = GO.selectPhylogeneticGraph [("best", "")] 0 ["best"] buildGraphs

       -- if not a tree then try to add net edges
       netAddArgs = [("netAdd", ""), ("keep", show 1), ("steepest", ""), ("atRandom", "")] 
       netGraphList = if (graphType inGS == Tree) then bestBuildGraphList
                      else R.netEdgeMaster netAddArgs inGS newData (randomIntegerList !! 2) bestBuildGraphList

       --simple swap refinement
       swapGraphList = if null swapOptions then netGraphList
                       else R.swapMaster swapOptions inGS newData (randomIntegerList !! 3) netGraphList
   in
   -- no data in there
   if (V.null . thd3) newData then emptyPhylogeneticGraph
   else head swapGraphList

-- | getGoodBremGraphs performs Goodman-Bremer support
-- can do sample of trees at random if specified
getGoodBremGraphs :: GlobalSettings -> ProcessedData -> Int -> [(String, String)] -> Maybe Int -> [PhylogeneticGraph] -> PhylogeneticGraph
getGoodBremGraphs inGS inData rSeed swapOptions sampleSize inGraphList = 
   emptyPhylogeneticGraph 

-- | resampleData perfoms a single randomized data resampling 
-- based on either with replacement (bootstrp) or without (jackknife)
-- jackknife moves through processed data and cretes a new data set
--    based on simple prob
-- bootStrap draws chars from input directly copying--currently disabled
-- if a block of data end up with zero resampled characters it is deleted
resampleData :: Int -> String -> Double -> ProcessedData -> ProcessedData
resampleData rSeed resampleType sampleFreq (nameVect, nameBVVect, blockDataVect) =
   if V.null blockDataVect then error "Null input data in resampleData"
   else
      let randomIntegerList = randomIntList rSeed
      in
      if resampleType == "bootstrap" then error "Bootstrap not currently implemented"
      else
         --Jackknife resampling
         let newBlockDataVect' = V.zipWith (resampleBlock resampleType sampleFreq) (V.fromList randomIntegerList) blockDataVect 
             -- filter any zero length blocks
             newBlockDataVect = V.filter ((not . V.null) . thd3) newBlockDataVect'
         in
         (nameVect, nameBVVect, newBlockDataVect)

-- | resampleBlock takes BlockData and a seed and creates a resampled BlockData
resampleBlock :: String -> Double -> Int -> BlockData -> BlockData
resampleBlock resampleType sampleFreq rSeed (nameText, charDataVV, charInfoV) =
   if resampleType == "bootstrap" then error "Bootsdtrap not currently implemented"
   else
      let randomIntegerList = randomIntList rSeed
          randomIntegerList2 = randomIntList (head randomIntegerList)
          acceptanceList = fmap (randAccept sampleFreq) randomIntegerList
          acceptanceList2 = fmap (randAccept sampleFreq) randomIntegerList2
          -- newCharInfoV = makeSampledVect acceptanceVect [] charInfoV
          -- newCharDataV = fmap (makeSampledVect acceptanceVect []) charDataVV
          (newCharDataVV, newCharInfoV) = V.unzip $ fmap (makeSampledPairVect acceptanceList acceptanceList2 [] [] charInfoV) charDataVV
          

      in
      trace ("RB length " ++ (show $ V.length charInfoV) ++ " -> " ++ (show $ V.length $ V.head newCharInfoV) )
      (nameText, newCharDataVV, V.head newCharInfoV)

   where randAccept b a = let (_, randVal) = divMod (abs a) 1000
                              critVal = floor (1000 * b)
                           in
                           -- trace ("RA : " ++ (show (b,a, randVal, critVal, randVal < critVal)))
                           randVal < critVal

-- | makeSampledCharCharInfoVect takes a vectot of Bool and a vector of cahrdata and a vector of charinfo 
-- if teh data type is not static--the character is returns if Bool is True not otherwise
-- if the char is static (add, non add, matrix) then the bool array is applied
-- across the vetor of those characters (since they re vectors of charcters themselves
-- returned as a pair of vectors (reversed--but shouldn't matter for resampling purposes) 
-- does not check if equal in length
makeSampledVect :: [Bool] -> [a] -> V.Vector a -> V.Vector a
makeSampledVect boolList accumList inVect  =
   if V.null inVect then 
    -- trace ("MSV R: " ++ (show $ length accumList)) 
    V.fromList accumList
   else
      -- trace ("MSV: " ++ (show $ head boolList)) (
      if head boolList then makeSampledVect (tail boolList) ((V.head inVect) : accumList) (V.tail inVect) 

      else makeSampledVect (tail boolList) accumList (V.tail inVect) 
      -- )
      
-- | makeSampledVect takes a liust of Bool and avector and returns those values 
-- with True as a vector (reversed--but shouldn't matter for resampling purposes) 
-- does not check if equal in length
makeSampledPairVect :: [Bool] -> [Bool] -> [CharacterData] -> [CharInfo] -> V.Vector CharInfo -> V.Vector CharacterData -> (V.Vector CharacterData, V.Vector CharInfo)
makeSampledPairVect fullBoolList boolList accumCharDataList accumCharInfoList inCharInfoVect inCharDataVect=
   if V.null inCharInfoVect then (V.fromList accumCharDataList, V.fromList accumCharInfoList)  
   else
      let firstCharInfo = V.head inCharInfoVect
          firstCharData = V.head inCharDataVect
          firstCharType = charType firstCharInfo
      in

      -- straight resample if dynamic
      if firstCharType `notElem` exactCharacterTypes
         then 
            if head boolList then makeSampledPairVect fullBoolList (tail boolList) (firstCharData : accumCharDataList) (firstCharInfo : accumCharInfoList) (V.tail inCharInfoVect) (V.tail inCharDataVect) 

            else makeSampledPairVect fullBoolList (tail boolList) accumCharDataList accumCharInfoList (V.tail inCharInfoVect) (V.tail inCharDataVect) 

      -- static character--keep in sample, but need to sample in the vector
      else
         let (a1, a2, a3) = rangePrelim firstCharData
             (na1, na2, na3) = stateBVPrelim firstCharData
             m1 = matrixStatesPrelim firstCharData
         in
         if firstCharType == Add then
            let newCharData = firstCharData {rangePrelim = (makeSampledVect fullBoolList [] a1, makeSampledVect fullBoolList [] a2, makeSampledVect fullBoolList [] a3)}
            in
            -- trace ("Length Add: " ++ (show $ V.length $ snd3 $ rangePrelim newCharData)) (
            if V.null (makeSampledVect fullBoolList [] a2) then makeSampledPairVect fullBoolList (tail boolList) accumCharDataList accumCharInfoList (V.tail inCharInfoVect) (V.tail inCharDataVect) 
            else makeSampledPairVect fullBoolList (tail boolList) (newCharData : accumCharDataList) (firstCharInfo : accumCharInfoList) (V.tail inCharInfoVect) (V.tail inCharDataVect) 
            -- )

         else if firstCharType == NonAdd then
            let newCharData = firstCharData {stateBVPrelim = (makeSampledVect fullBoolList [] na1, makeSampledVect fullBoolList [] na2, makeSampledVect fullBoolList [] na3)}
            in
            -- trace ("Length NonAdd: " ++ (show $ V.length $ snd3 $ stateBVPrelim newCharData))  (
            if V.null (makeSampledVect fullBoolList [] na2) then makeSampledPairVect fullBoolList (tail boolList) accumCharDataList accumCharInfoList(V.tail inCharInfoVect)  (V.tail inCharDataVect) 
            else makeSampledPairVect fullBoolList (tail boolList) (newCharData : accumCharDataList) (firstCharInfo : accumCharInfoList) (V.tail inCharInfoVect) (V.tail inCharDataVect)
            -- )

         else if firstCharType == Matrix then
            let  newCharData = firstCharData {matrixStatesPrelim = (makeSampledVect fullBoolList [] m1)}
            in
            if V.null (makeSampledVect fullBoolList [] m1) then makeSampledPairVect fullBoolList (tail boolList) accumCharDataList accumCharInfoList (V.tail inCharInfoVect) (V.tail inCharDataVect) 
            else makeSampledPairVect fullBoolList (tail boolList) (newCharData : accumCharDataList) (firstCharInfo : accumCharInfoList) (V.tail inCharInfoVect) (V.tail inCharDataVect)


         else error ("Incorrect character type in makeSampledPairVect: " ++ show firstCharType)

