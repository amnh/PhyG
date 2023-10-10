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

import Commands.Verify qualified as VER
import Data.Char
import Data.List qualified as L
import Data.Maybe
import Data.Vector qualified as V
import Data.Vector.Generic qualified as GV
import Data.Vector.Unboxed qualified as UV
import Debug.Trace
import GeneralUtilities
import GraphOptimization.Traversals qualified as T
import Graphs.GraphOperations qualified as GO
import ParallelUtilities qualified as PU
import Reconciliation.ReconcileGraphs qualified as REC
import Search.Build qualified as B
import Search.NetworkAddDelete qualified as N
import Search.Refinement qualified as R
import Text.Read
import Types.Types
import Utilities.Distances qualified as DD
import Utilities.LocalGraph qualified as LG

-- "High level" paralleization used for overall graphs contruction

-- | driver for overall support
supportGraph :: [Argument] -> GlobalSettings -> ProcessedData -> Int -> [ReducedPhylogeneticGraph] -> [ReducedPhylogeneticGraph]
supportGraph inArgs inGS inData rSeed inGraphList =
   if null inGraphList then error "No graphs input to calculate support"
   else
      let fstArgList = fmap (fmap toLower . fst) inArgs
          sndArgList = fmap (fmap toLower . snd) inArgs
          lcArgList = zip fstArgList sndArgList
          checkCommandList = checkCommandArgs "support" fstArgList VER.supportArgList
     in
     -- check for valid command options
     if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'support': " <> show inArgs)
     else
         let supportMeasure
               | any (( == "bootstrap").fst) lcArgList = Bootstrap
               | any (( == "jackknife").fst) lcArgList = Jackknife
               | any (( == "goodmanbremer").fst) lcArgList = GoodmanBremer
               | otherwise = GoodmanBremer

             useSPR = any (( == "spr").fst) lcArgList 
             useTBR = any (( == "spr").fst) lcArgList 


             onlyBuild = any (( == "buildonly").fst) lcArgList

             jackList   = filter (( == "jackknife").fst) lcArgList
             jackFreq'
              | length jackList > 1 =
                errorWithoutStackTrace ("Multiple jackknife sampling frequency specifications in support command--can have only one (e.g. jackknife:0.62): " <> show inArgs)
              | null jackList = Just 0.6321 -- 1- 1/e
              | null (snd $ head jackList) = Just 0.6321
              | otherwise = readMaybe (snd $ head jackList) :: Maybe Double

             replicateList   = filter (( == "replicates").fst) lcArgList
             replicates'
              | length replicateList > 1 =
                errorWithoutStackTrace ("Multiple resampling replicate specifications in support command--can have only one (e.g. replicates:100): " <> show inArgs)
              | null replicateList = Just 100
              | otherwise = readMaybe (snd $ head replicateList) :: Maybe Int

             goodBremList   = filter ((`elem` ["goodmanbremer", "gb"]).fst) lcArgList
             {-
             goodBremMethod
              | length goodBremList > 1 =
                errorWithoutStackTrace ("Multiple Goodman-Bremer method specifications in support command--can have only one (e.g. gb:tbr): " <> show inArgs)
              | null (snd $ head goodBremList) = Just "tbr"
              | otherwise = Just $ snd $ head goodBremList
              -}

             goodBremSampleList   = filter ((`elem` ["gbsample"]).fst) lcArgList
             goodBremSample
              | length goodBremSampleList > 1 =
                errorWithoutStackTrace ("Multiple Goodman-Bremer sample specifications in support command--can have only one (e.g. gbsample:1000): " <> show inArgs)
              | null goodBremSampleList = Just (maxBound :: Int)
              | otherwise = readMaybe (snd $ head goodBremSampleList) :: Maybe Int

         in
         if isNothing jackFreq' then errorWithoutStackTrace ("Jacknife frequency not a float (e.g. jackknife:0.5) in support: " <> show (snd $ head jackList))
         else if isNothing replicates' then errorWithoutStackTrace ("Resampling replicates specification not a string (e.g. replicates:100) in support: " <> show (snd $ head replicateList))
         --else if isNothing goodBremMethod then errorWithoutStackTrace ("Goodman-Bremer method specification not a string (e.g. goodmanBremer:SPR) in support: " <> (show (snd $ head goodBremList)) <> (show lcArgList))
         else if isNothing goodBremSample then errorWithoutStackTrace ("Goodman-Bremer sample specification not an integer (e.g. gbsample:1000) in support: " <> show (snd $ head goodBremSampleList))

         else
            let thisMethod
                  | (supportMeasure == Bootstrap) && ((not . null) jackList && null goodBremList) = trace "Bootstrap and Jackknife specified--defaulting to Jackknife"
                                Jackknife
                  | (supportMeasure == Bootstrap) || ((not . null) jackList && (not . null) goodBremList) = trace "Resampling (Bootstrap or Jackknife) and Goodman-Bremer specified--defaulting to Goodman-Bremer"
                                GoodmanBremer
                  | supportMeasure == Bootstrap = Bootstrap
                  | (not . null) jackList = Jackknife
                  | otherwise = GoodmanBremer

                gbSampleSize = if goodBremSample == Just (maxBound :: Int) then Nothing
                               else goodBremSample

                -- sample trees uniformly at random--or "nth"
                gbRandomSample = if isJust gbSampleSize then True -- any (( == "atrandom").fst) lcArgList
                                 else False

                replicates = if fromJust replicates' < 0 then
                                 trace "Negative replicates number--defaulting to 100"
                                 100
                             else fromJust replicates'
                jackFreq = if fromJust jackFreq' <= 0 || fromJust jackFreq' >= 1.0 then
                              trace "Jackknife frequency must be on (0.0, 1.0) defaulting to 0.6321"
                              0.6321
                           else fromJust jackFreq'

                buildOptions = [("distance",""), ("replicates", show (100 :: Int)), ("best", show (1 :: Int)), ("rdwag", "")] -- [("replicates", show 10), ("best", show 1)]
                swapOptions = if onlyBuild then []
                              else [("tbr", ""), ("steepest", ""), ("keep", show (1 :: Int))]
                supportGraphList = if thisMethod == Bootstrap || thisMethod == Jackknife then
                                       let extraString = if  thisMethod == Jackknife then " with delete fraction  " <> show (1 - jackFreq)
                                                         else ""
                                       in
                                       trace ("Generating " <> show thisMethod <> " resampling support with " <> show replicates <> " replicates" <> extraString)
                                       [getResampleGraph inGS inData rSeed thisMethod replicates buildOptions swapOptions jackFreq]
                                   else
                                       let neighborhood = if useTBR then "tbr"
                                                          else if useSPR then "spr"
                                                          else "tbr"
                                           extraString = if isJust gbSampleSize then 
                                                " using " <> neighborhood <> " based on " <> show (fromJust gbSampleSize) <> " samples at random"
                                                         else " using " <> neighborhood
                                       in
                                       trace ("Generating Goodman-Bremer support" <> extraString)
                                       fmap (getGoodBremGraphs inGS inData rSeed neighborhood gbSampleSize gbRandomSample) inGraphList
            in

            supportGraphList


-- | getResampledGraphs performs resampling and search for Bootstrap and jackknife support
getResampleGraph :: GlobalSettings
                 -> ProcessedData
                 -> Int
                 -> SupportMethod
                 -> Int
                 -> [(String, String)]
                 -> [(String, String)]
                 -> Double
                 -> ReducedPhylogeneticGraph
getResampleGraph inGS inData rSeed resampleType replicates buildOptions swapOptions jackFreq =
   let resampledGraphList = PU.seqParMap (parStrategy $ strictParStrat inGS) (makeResampledDataAndGraph inGS inData resampleType buildOptions swapOptions jackFreq) (take replicates $ randomIntList rSeed) -- `using` PU.myParListChunkRDS
       -- create appropriate support graph >50% ?
       -- need to add args
       reconcileArgs = if graphType inGS == Tree then [("method","majority"), ("compare","identity"), ("edgelabel","true"), ("vertexlabel","true"), ("connect","true"), ("threshold","51"), ("outformat", "dot")]
                       else [("method","eun"), ("compare","identity"), ("edgelabel","true"),  ("vertexlabel","true"), ("connect","true"), ("threshold","51"),("outformat", "dot")]
         -- majority ruke consensus if no args
       (_, reconciledGraph) = REC.makeReconcileGraph VER.reconcileArgList reconcileArgs (fmap fst5 resampledGraphList)
   in
   -- trace ("GRG: \n" <> reconciledGraphString) (
   -- generate resampled graph
   -- can't really relabel  easily wihtout bv and maybe not necessary anyway--node numebrs inconsistent
  (reconciledGraph, infinity, LG.empty, V.empty, V.empty)
   -- )

-- | makeResampledDataAndGraph takes paramters, resmaples data and find a graph based on search parameters
-- returning the resampled graph
makeResampledDataAndGraph :: GlobalSettings
                          -> ProcessedData
                          -> SupportMethod
                          -> [(String, String)]
                          -> [(String, String)]
                          -> Double
                          -> Int
                          -> ReducedPhylogeneticGraph
makeResampledDataAndGraph inGS inData resampleType buildOptions swapOptions jackFreq rSeed =
   let randomIntegerList1 = randomIntList rSeed
       -- create resampled data
       newData = resampleData (head randomIntegerList1) resampleType jackFreq inData

       -- pairwise distances for distance analysis
       pairwiseDistances = DD.getPairwiseDistances newData

       -- build graphs
       buildGraphs = B.buildGraph buildOptions inGS newData pairwiseDistances (randomIntegerList1 !! 1)
       bestBuildGraphList = GO.selectGraphs Best (maxBound::Int) 0.0 (-1) buildGraphs

       -- if not a tree then try to add net edges
       netAddArgs = [("netadd", ""), ("keep", show (1 :: Int)), ("steepest", ""), ("atrandom", ""), ("maxnetedges","5")]
       netGraphList = if graphType inGS == Tree then bestBuildGraphList
                      else R.netEdgeMaster netAddArgs inGS newData (randomIntegerList1 !! 2) bestBuildGraphList

       --simple swap refinement
       swapGraphList = if null swapOptions then netGraphList
                       else R.swapMaster swapOptions inGS newData (randomIntegerList1 !! 3) netGraphList
   in
   -- no data in there
   if (V.null . thd3) newData then emptyReducedPhylogeneticGraph
   else head swapGraphList

-- | resampleData perfoms a single randomized data resampling
-- based on either with replacement (bootstrp) or without (jackknife)
-- jackknife moves through processed data and creates a new data set
--    based on simple prob
-- Bootstrap draws chars from input directly copying--currently disabled
-- if a block of data end up with zero resampled characters it is deleted
resampleData :: Int -> SupportMethod -> Double -> ProcessedData -> ProcessedData
resampleData rSeed resampleType sampleFreq (nameVect, nameBVVect, blockDataVect) =
   if V.null blockDataVect then error "Null input data in resampleData"
   else
      let lRandomIntegerList = randomIntList rSeed
      in
         --Bootstrap  or Jackknife resampling
         let newBlockDataVect' = if resampleType == Bootstrap then V.zipWith resampleBlockBootstrap (V.fromList lRandomIntegerList) blockDataVect
                                 else V.zipWith (resampleBlockJackknife  sampleFreq) (V.fromList lRandomIntegerList) blockDataVect
             -- filter any zero length blocks
             newBlockDataVect = V.filter ((not . V.null) . thd3) newBlockDataVect'
         in
         (nameVect, nameBVVect, newBlockDataVect)


-- | resampleBlockBootstrap takes BlockData and a seed and creates a Bootstrap resampled BlockData
resampleBlockBootstrap :: Int -> BlockData -> BlockData
resampleBlockBootstrap rSeed (nameText, charDataVV, charInfoV) =
      let lRandomIntegerList = randomIntList rSeed
          randomIntegerList2 = randomIntList (head lRandomIntegerList)

          -- maps over taxa in data bLock
          (newCharDataVV, newCharInfoV) = V.unzip $ fmap (makeSampledPairVectBootstrap lRandomIntegerList randomIntegerList2 charInfoV) charDataVV
      in
      (nameText, newCharDataVV, V.head newCharInfoV)

-- | makeSampledPairVectBootstrap takes a list of Int and a vectors of charinfo and char data
-- and returns new vectors of chardata and charinfo based on randomly sampled character indices
-- this to create a Bootstrap replicate of equal size
-- this for a single taxon hecen pass teh random ints so same for each one
makeSampledPairVectBootstrap :: [Int] -> [Int] -> V.Vector CharInfo -> V.Vector CharacterData -> (V.Vector CharacterData, V.Vector CharInfo)
makeSampledPairVectBootstrap fullRandIntlList randIntList inCharInfoVect inCharDataVect =
      -- get character numbers and set resampling indices for dynamic characters--statric happen within those characters
      let -- zip so can filter static and dynamic characters
          dataInfoPairV = V.zip inCharDataVect inCharInfoVect

          -- filter static
          (staticCharsV, staticCharsInfoV)  = V.unzip $ V.filter ((`elem` exactCharacterTypes) . charType . snd ) dataInfoPairV

          -- filter dynamic
          (dynamicCharsV, dynamicCharsInfoV) = V.unzip $ V.filter ((`notElem` exactCharacterTypes) . charType .snd ) dataInfoPairV

          numDynamicChars   = V.length dynamicCharsV
          dynCharIndices   = fmap (randIndex numDynamicChars) (take numDynamicChars randIntList)

      in
      -- trace ("MSPVB: " <> (show $ length dynCharIndices) <> " " <> (show dynCharIndices)) (
      -- resample dynamic characters
         -- straight resample
      let resampleDynamicChars    = V.map (dynamicCharsV V.!) (V.fromList dynCharIndices)
          resmapleDynamicCharInfo = V.map (dynamicCharsInfoV V.!) (V.fromList dynCharIndices)

          -- static chars do each one mapping random choices within the character type
          -- but keeping each one--hence char info is staticCharsInfoV
          resampleStaticChars    = V.zipWith (subSampleStatic fullRandIntlList) staticCharsV staticCharsInfoV

      in
      -- cons the vectors for chrater data and character info
      (resampleStaticChars <> resampleDynamicChars, staticCharsInfoV <> resmapleDynamicCharInfo)
      -- )
      where randIndex :: forall {b}. Integral b => b -> b -> b
            randIndex a b  = snd $  divMod (abs b) a


-- | subSampleStatic takes a random int list and a static charcter
-- Bootstrap resamples that character based on ransom it list and number of "subcharacters" in character
subSampleStatic :: [Int] -> CharacterData -> CharInfo -> CharacterData
subSampleStatic randIntList inCharData inCharInfo =
   let (a1, a2, a3) = rangePrelim inCharData
       (na1, na2, na3) = stateBVPrelim inCharData
       (pa1, pa2, pa3) = packedNonAddPrelim inCharData
       m1 = matrixStatesPrelim inCharData
       inCharType = charType inCharInfo

       charLength
         | inCharType == Add = V.length a2
         | inCharType == NonAdd = V.length na2
         | inCharType == Matrix = V.length m1
         | otherwise = error ("Dynamic character in subSampleStatic: " <> show inCharType)

       -- get character indices based on number "subcharacters"
       staticCharIndices = V.fromList $ fmap (randIndex charLength) (take charLength randIntList)
       staticCharIndicesUV = UV.fromList $ fmap (randIndex charLength) (take charLength randIntList)


   in
   -- trace ("SSS:" <> (show $ V.length staticCharIndices) <> " " <> (show staticCharIndices)) (
   if inCharType == Add then
         inCharData {rangePrelim = (V.map (a1 V.!) staticCharIndices, V.map (a2 V.!) staticCharIndices, V.map (a3 V.!) staticCharIndices)}

   else if inCharType == NonAdd then
         inCharData {stateBVPrelim = (V.map (na1 V.!) staticCharIndices, V.map (na2 V.!) staticCharIndices, V.map (na3 V.!) staticCharIndices)}

   else if inCharType `elem` packedNonAddTypes then
         inCharData {packedNonAddPrelim = (UV.map (pa1 UV.!) staticCharIndicesUV, UV.map (pa2 UV.!) staticCharIndicesUV, UV.map (pa3 UV.!) staticCharIndicesUV)}

   else if inCharType == Matrix then
         inCharData {matrixStatesPrelim = V.map (m1 V.!) staticCharIndices}

   else error ("Incorrect character type in subSampleStatic: " <> show inCharType)
   -- )
   where randIndex :: forall {b}. Integral b => b -> b -> b
         randIndex a b  = snd $  divMod (abs b) a

-- | makeSampledCharCharInfoVect takes a vector of Int and a vector of charData and a vector of charinfo
-- if teh data type is not static--the character is returns if Bool is True not otherwise
-- if the char is static (add, non add, matrix) then the bool array is applied
-- across the vetor of those characters (since they re vectors of charcters themselves
-- returned as a pair of vectors (reversed--but shouldn't matter for resampling purposes)
-- does not check if equal in length
makeSampledVect ::  (GV.Vector v a) => [Bool] -> [a] -> v a -> v a
makeSampledVect boolList accumList inVect  =
   if GV.null inVect then
    -- trace ("MSV R: " <> (show $ length accumList))
    GV.fromList accumList
   else
      -- trace ("MSV: " <> (show $ head boolList)) (
      if head boolList then makeSampledVect (tail boolList) (GV.head inVect : accumList) (GV.tail inVect)

      else makeSampledVect (tail boolList) accumList (GV.tail inVect)
      -- )

-- | makeSampledVect takes a list of Bool and avector and returns those values
-- with True as a vector (reversed--but shouldn't matter for resampling purposes)
-- does not check if equal in length
makeSampledPairVect :: [Bool] -> [Bool] -> [CharacterData] -> [CharInfo] -> V.Vector CharInfo -> V.Vector CharacterData -> (V.Vector CharacterData, V.Vector CharInfo)
makeSampledPairVect fullBoolList boolList accumCharDataList accumCharInfoList inCharInfoVect inCharDataVect =
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
             (pa1, pa2, pa3) = packedNonAddPrelim firstCharData
             m1 = matrixStatesPrelim firstCharData
         in
         if firstCharType == Add then
            let newCharData = firstCharData {rangePrelim = (makeSampledVect fullBoolList [] a1, makeSampledVect fullBoolList [] a2, makeSampledVect fullBoolList [] a3)}
            in
            -- trace ("Length Add: " <> (show $ V.length $ snd3 $ rangePrelim newCharData)) (
            if V.null (makeSampledVect fullBoolList [] a2) then makeSampledPairVect fullBoolList (tail boolList) accumCharDataList accumCharInfoList (V.tail inCharInfoVect) (V.tail inCharDataVect)
            else makeSampledPairVect fullBoolList (tail boolList) (newCharData : accumCharDataList) (firstCharInfo : accumCharInfoList) (V.tail inCharInfoVect) (V.tail inCharDataVect)
            -- )

         else if firstCharType == NonAdd then
            let newCharData = firstCharData {stateBVPrelim = (makeSampledVect fullBoolList [] na1, makeSampledVect fullBoolList [] na2, makeSampledVect fullBoolList [] na3)}
            in
            -- trace ("Length NonAdd: " <> (show $ V.length $ snd3 $ stateBVPrelim newCharData))  (
            if V.null (makeSampledVect fullBoolList [] na2) then makeSampledPairVect fullBoolList (tail boolList) accumCharDataList accumCharInfoList (V.tail inCharInfoVect)  (V.tail inCharDataVect)
            else makeSampledPairVect fullBoolList (tail boolList) (newCharData : accumCharDataList) (firstCharInfo : accumCharInfoList) (V.tail inCharInfoVect) (V.tail inCharDataVect)
            -- )

         else if firstCharType `elem` packedNonAddTypes then
            let newCharData = firstCharData {packedNonAddPrelim = (makeSampledVect fullBoolList [] pa1, makeSampledVect fullBoolList [] pa2, makeSampledVect fullBoolList [] pa3)}
            in
            -- trace ("Length NonAdd: " <> (show $ V.length $ snd3 $ stateBVPrelim newCharData))  (
            if GV.null (makeSampledVect fullBoolList [] pa2) then makeSampledPairVect fullBoolList (tail boolList) accumCharDataList accumCharInfoList (GV.tail inCharInfoVect)  (GV.tail inCharDataVect)
            else makeSampledPairVect fullBoolList (tail boolList) (newCharData : accumCharDataList) (firstCharInfo : accumCharInfoList) (GV.tail inCharInfoVect) (GV.tail inCharDataVect)
            -- )

         else if firstCharType == Matrix then
            let  newCharData = firstCharData {matrixStatesPrelim = makeSampledVect fullBoolList [] m1}
            in
            if V.null (makeSampledVect fullBoolList [] m1) then makeSampledPairVect fullBoolList (tail boolList) accumCharDataList accumCharInfoList (V.tail inCharInfoVect) (V.tail inCharDataVect)
            else makeSampledPairVect fullBoolList (tail boolList) (newCharData : accumCharDataList) (firstCharInfo : accumCharInfoList) (V.tail inCharInfoVect) (V.tail inCharDataVect)


         else error ("Incorrect character type in makeSampledPairVect: " <> show firstCharType)

-- | resampleBlockJackknife takes BlockData and a seed and creates a jackknife resampled BlockData
resampleBlockJackknife :: Double -> Int -> BlockData -> BlockData
resampleBlockJackknife sampleFreq rSeed inData@(nameText, charDataVV, charInfoV) =
      let randomIntegerList1 = randomIntList rSeed
          randomIntegerList2 = randomIntList (head randomIntegerList1)
          acceptanceList = fmap (randAccept sampleFreq) randomIntegerList1
          acceptanceList2 = fmap (randAccept sampleFreq) randomIntegerList2
          -- newCharInfoV = makeSampledVect acceptanceVect [] charInfoV
          -- newCharDataV = fmap (makeSampledVect acceptanceVect []) charDataVV
          (newCharDataVV, newCharInfoV) = V.unzip $ fmap (makeSampledPairVect acceptanceList acceptanceList2 [] [] charInfoV) charDataVV


      in
      -- trace ("RB length " <> (show $ V.length charInfoV) <> " -> " <> (show $ V.length $ V.head newCharInfoV) ) (
      if V.length (V.head newCharInfoV) > 0 then (nameText, newCharDataVV, V.head newCharInfoV)
      else resampleBlockJackknife sampleFreq (head randomIntegerList1) inData
      -- )

   where randAccept :: forall {p} {a}.
                      (RealFrac p, Integral a) =>
                      p -> a -> Bool
         randAccept b a = let (_, randVal) = divMod (abs a) 1000
                              critVal = floor (1000 * b)
                           in
                           -- trace ("RA : " <> (show (b,a, randVal, critVal, randVal < critVal)))
                           randVal < critVal

-- | getGoodBremGraphs performs Goodman-Bremer support
-- examines complete SPR or TBR swap neighborhood chekcing for presence/absence of edges in input Graph List
-- can do sample of trees either "nth" or at random if specified
-- sample based on SPR-- 4n^2 - 26n - 42 for TBR 8n^3 for now
-- this will only examine bridge edges for networks, networkedge values willl be doen via net delete
-- MAPs for each graph?
getGoodBremGraphs :: GlobalSettings -> ProcessedData -> Int -> String -> Maybe Int -> Bool -> ReducedPhylogeneticGraph -> ReducedPhylogeneticGraph
getGoodBremGraphs inGS inData rSeed swapType sampleSize sampleAtRandom inGraph =
   if LG.isEmpty (fst5 inGraph) then error "Null graph in getGoodBremGraphs" -- maybe should be error?
   else
      -- create list of edges for input graph and a structure with egde node indices and bitvector values
      -- requires index BV of each node
      {-
      let egdeList = LG.edges (fst5 inGraph)

          -- graph node list
          nodeList = LG.labNodes (thd5 inGraph)
          nodeIndexBVPairList = fmap makeindexBVPair nodeList

          -- list of vectors for contant time access via index = fst (a, bv)
          nodeIndexBVPairVect= V.fromList nodeIndexBVPairList

          -- make tuple for each edge in each graph
          -- (uIndex,vINdex,uBV, vBV, graph cost)
          tupleList = makeGraphEdgeTuples nodeIndexBVPairVect infinity egdeList
      -}
      let tupleList = getGraphTupleList inGraph

          -- traverse neighborhood (and net edge removal) keeping min cost without edges
          supportEdgeTupleList = getGBTuples inGS inData rSeed swapType sampleSize sampleAtRandom tupleList inGraph

          simpleGBGraph = LG.mkGraph (LG.labNodes $ fst5 inGraph) (fmap (tupleToSimpleEdge (snd5 inGraph)) supportEdgeTupleList)
      in
      -- trace ("GGBG: " <> (show $ length tupleList) <> " -> " <> (show $ length supportEdgeTupleList))
      (simpleGBGraph, snd5 inGraph, thd5 inGraph, fth5 inGraph, fft5 inGraph)

      where tupleToSimpleEdge :: forall {c1} {a} {b} {c2} {d}.
                             Num c1 =>
                             c1 -> (a, b, c2, d, c1) -> (a, b, c1)
            tupleToSimpleEdge d (a,b, _, _, c) = (a, b, c - d)

-- | getGraphTupleList takes a graph and cost (maybe initialized to infinity) returns tuple list
getGraphTupleList :: ReducedPhylogeneticGraph -> [(Int, Int, NameBV, NameBV, VertexCost)]
getGraphTupleList inGraph =
   if LG.isEmpty (fst5 inGraph) then error "Null graph in getGraphTupleList"
   else
      let egdeList = LG.edges (fst5 inGraph)

          -- graph node list
          nodeList = LG.labNodes (thd5 inGraph)
          nodeIndexBVPairList = fmap makeindexBVPair nodeList

          -- list of vectors for contant time access via index = fst (a, bv)
          nodeIndexBVPairVect= V.fromList nodeIndexBVPairList

          -- make tuple for each edge in each graph
          -- (uIndex,vINdex,uBV, vBV, graph cost)
          tupleList = makeGraphEdgeTuples nodeIndexBVPairVect infinity egdeList
      in
      tupleList
      where makeindexBVPair :: forall {a}. (a, VertexInfo) -> (a, NameBV)
            makeindexBVPair (a,b) = (a, bvLabel b)

-- | getGBTuples takes a tuple list from graph containing initialized values and update those values based
-- on each graph in the inGraph neigborhood
-- first does this via swap--for network does edge net edge in turn by removing using netDel
getGBTuples :: GlobalSettings
            -> ProcessedData
            -> Int
            -> String
            -> Maybe Int
            -> Bool
            -> [(Int, Int, NameBV, NameBV, VertexCost)]
            -> ReducedPhylogeneticGraph
            -> [(Int, Int, NameBV, NameBV, VertexCost)]
getGBTuples inGS inData rSeed swapType sampleSize sampleAtRandom inTupleList inGraph =
    -- traverse swap (SPR/TBR) neighborhood optimizing each graph fully
    let swapTuples = performGBSwap inGS inData rSeed swapType sampleSize sampleAtRandom inTupleList inGraph

        -- network edge support if not Tree
        netTuples = if (graphType inGS == Tree) || LG.isTree (fst5 inGraph) then
                        -- swap only for Tree-do nothing
                        swapTuples

                    -- SoftWired => delete edge -- could add net move if needed
                     else if graphType inGS == SoftWired then
                        PU.seqParMap (parStrategy $ strictParStrat inGS) (updateDeleteTuple inGS inData inGraph) swapTuples -- `using` PU.myParListChunkRDS

                    -- HardWired => move edge
                    else
                        PU.seqParMap (parStrategy $ strictParStrat inGS) (updateMoveTuple inGS inData inGraph) swapTuples -- `using` PU.myParListChunkRDS
        in
        netTuples

-- | updateDeleteTuple take a graph and and edge and delete a network edge (or retunrs tuple if not network)
-- if this were a HardWWired graph--cost would always go down, so only applied to softwired graphs
updateDeleteTuple :: GlobalSettings 
                  -> ProcessedData 
                  -> ReducedPhylogeneticGraph 
                  -> (Int, Int, NameBV, NameBV, VertexCost) 
                  -> (Int, Int, NameBV, NameBV, VertexCost)
updateDeleteTuple inGS inData inGraph inTuple@(inE, inV, inEBV, inVBV, inCost) =
   let isNetworkEdge = LG.isNetworkEdge (fst5 inGraph) (inE, inV)
   in
   if not isNetworkEdge then inTuple
   else
      -- True to force full evalutation
      let deleteCost = snd5 $ N.deleteNetEdge inGS inData inGraph True (inE, inV)
      in
      (inE, inV, inEBV, inVBV, min inCost deleteCost)


-- | updateMoveTuple take a graph and and edge and moves a network edge (or returns tuple if not network)
-- if this were a HardWWired graph--cost would always go down, so only applied to softwired graphs
    -- max bound because its a place holder for max num net edges
updateMoveTuple :: GlobalSettings 
                -> ProcessedData
                ->  ReducedPhylogeneticGraph 
                -> (Int, Int, NameBV, NameBV, VertexCost) 
                -> (Int, Int, NameBV, NameBV, VertexCost)
updateMoveTuple inGS inData inGraph inTuple@(inE, inV, inEBV, inVBV, inCost) =
   let isNetworkEdge = LG.isNetworkEdge (fst5 inGraph) (inE, inV)
   in
   if not isNetworkEdge then inTuple
   else
      -- True to force full evalutation
      let steepest = False
          randomOrder = False
          keepNum = 10 -- really could be one since sorted by cost, but just to make sure)Order
          rSeed = 0
          
          saParams :: forall {a}. Maybe a
          saParams = Nothing
          
          moveCost = minimum (snd5 <$> N.deleteOneNetAddAll inGS inData (maxBound :: Int) keepNum steepest randomOrder inGraph [(inE, inV)] rSeed saParams)
      in
      (inE, inV, inEBV, inVBV, min inCost moveCost)



-- | performGBSwap takes parameters and  graphs and traverses swap neighborhood
-- examining each (or nth, or random) Graphs examining each ecah in each graph for Goodman-Bremer
-- optimality support
performGBSwap   :: GlobalSettings
                -> ProcessedData
                -> Int
                -> String
                -> Maybe Int
                -> Bool
                -> [(Int, Int, NameBV, NameBV, VertexCost)]
                -> ReducedPhylogeneticGraph
                -> [(Int, Int, NameBV, NameBV, VertexCost)]
performGBSwap inGS inData rSeed swapType sampleSize sampleAtRandom inTupleList inGraph =
    if LG.isEmpty (fst5 inGraph) then error "Null graph in performGBSwap"
    else
        let -- work with simple graph
            inSimple = fst5 inGraph
            (firstRootIndex, _) = head $ LG.getRoots inSimple

            -- determine edges to break on--'bridge' edges only for network
            -- filter out edges from root since no use--would just rejoin
            breakEdgeList = if graphType inGS == Tree then filter ((/= firstRootIndex) . fst3) $ LG.labEdges inSimple
                          else filter ((/= firstRootIndex) . fst3) $ LG.getEdgeSplitList inSimple

            -- get random integer lists for swap
            lRandomIntegerList = randomIntList rSeed
            randomIntegerListList = fmap randomIntList lRandomIntegerList

            -- integerized critical value for prob accept
            -- based on approx (leaves - netnodes)^2 or (leaves - netnodes)^3
            (_, leafList, _, netVertList) = LG.splitVertexList (fst5 inGraph)
            intProbAccept = if swapType == "spr" then floor ((1000.0 * fromIntegral (fromJust sampleSize)) / ((2.0 * fromIntegral (length leafList - length netVertList)) ** 2) :: Double)
                            else floor ((1000.0 * fromIntegral (fromJust sampleSize)) / ((2.0 * fromIntegral (length leafList - length netVertList)) ** 3) :: Double)


            -- generate tuple lists for each break edge parallelized at this level
            tupleListList = PU.seqParMap (parStrategy $ strictParStrat inGS) (splitRejoinGB' inGS inData swapType intProbAccept sampleAtRandom inTupleList inSimple breakEdgeList) (zip randomIntegerListList breakEdgeList)  -- `using` PU.myParListChunkRDS

            -- merge tuple lists--should all be in same order
            newTupleList = mergeTupleLists (filter (not . null) tupleListList) []
        in
        -- trace ("PGBS:" <> (show $ fmap length tupleListList) <> " -> " <> (show $ length newTupleList))
        newTupleList


-- | splitRejoinGB' is  wrapper for splitRejoinGB to allow for seqParMap
splitRejoinGB'   :: GlobalSettings
                -> ProcessedData
                -> String
                -> Int
                -> Bool
                -> [(Int, Int, NameBV, NameBV, VertexCost)]
                -> SimpleGraph
                -> [LG.LEdge Double]
                -> ([Int], LG.LEdge Double)
                -> [(Int, Int, NameBV, NameBV, VertexCost)]
splitRejoinGB' inGS inData swapType intProbAccept sampleAtRandom inTupleList inGraph originalBreakEdgeList (inRandomIntegerList, breakEdge) =
        splitRejoinGB inGS inData swapType intProbAccept sampleAtRandom inTupleList inGraph originalBreakEdgeList inRandomIntegerList breakEdge

-- | splitRejoinGB take parameters and splits input graph at specified edge and rejoins at all available edge
-- (reroots the pruned subgraph if TBR) and creates and gets cost of graph (lazy takes care of post order only)
-- with optimized graph, tuple list is creted and compared to input graph tuple list.
-- original edge was (parentPrunedGraphRoot, prunedGraphRootIndex)
-- working with SimpleGraph
splitRejoinGB   :: GlobalSettings
                -> ProcessedData
                -> String
                -> Int
                -> Bool
                -> [(Int, Int, NameBV, NameBV, VertexCost)]
                -> SimpleGraph
                -> [LG.LEdge Double]
                -> [Int]
                -> LG.LEdge Double
                -> [(Int, Int, NameBV, NameBV, VertexCost)]
splitRejoinGB inGS inData swapType intProbAccept sampleAtRandom inTupleList inGraph originalBreakEdgeList inRandomIntegerList breakEdge =

    let
      -- split graph on breakEdge
      (splitGraph, _, prunedGraphRootIndex,  _, _, edgeDeleteList) = LG.splitGraphOnEdge' inGraph breakEdge

      -- get edges in base graph to be invaded (ie not in pruned graph)
      prunedGraphRootNode = (prunedGraphRootIndex, fromJust $ LG.lab splitGraph prunedGraphRootIndex)
      (prunedSubTreeNodes, prunedSubTreeEdges) = LG.nodesAndEdgesAfter splitGraph [prunedGraphRootNode]
      edgesNotToInvade = (LG.toEdge breakEdge : edgeDeleteList) <> fmap LG.toEdge prunedSubTreeEdges
      edgesToInvade = filter (LG.notMatchEdgeIndices edgesNotToInvade) originalBreakEdgeList

      -- rejoin, evaluate, get better tuple
      -- check if there are tbr-type rearrangements to do (rerooting pruned graph)
      -- create TBR rerooot split graphs if required
      splitGraphList = if (length prunedSubTreeNodes < 3) || swapType == "spr" then [splitGraph]

                       -- generate "tbr" rerootings in split graph
                       else getTBRSplitGraphs inGS splitGraph breakEdge

      -- new random lists for rejoin
      randomIntegerListList = fmap randomIntList inRandomIntegerList

      -- parallel at break level above
      rejoinTupleListList = zipWith (rejoinGB inGS inData intProbAccept sampleAtRandom inTupleList splitGraphList breakEdge) randomIntegerListList edgesToInvade

      -- merge tuples
      newTupleList = mergeTupleLists rejoinTupleListList []
    in
    newTupleList

-- | rejoinGB rejoins split graph at specific edge, id SPR then that's it, if TBR reroot pruned subgraph
-- splitGraph is SimpleGraph
-- the rejoin is SPR type relying on teh list lengt of split graph to present the TBR reroots
rejoinGB :: GlobalSettings
         -> ProcessedData
         -> Int
         -> Bool
         -> [(Int, Int, NameBV, NameBV, VertexCost)]
         -> [SimpleGraph]
         -> LG.LEdge Double
         -> [Int]
         -> LG.LEdge Double
         -> [(Int, Int, NameBV, NameBV, VertexCost)]
rejoinGB inGS inData intProbAccept sampleAtRandom inTupleList splitGraphList originalBreakEdge@(eBreak, _, _) randIntList edgeToInvade =
   if null splitGraphList then inTupleList
   else
      let splitGraph = head  splitGraphList
          doGraph = (not sampleAtRandom || (let (_, intRandVal) = divMod (abs (head randIntList)) 1000
                                            in
                                            intRandVal < intProbAccept))
      in
      if doGraph then
         let newGraph = LG.joinGraphOnEdge splitGraph edgeToInvade eBreak
             pruneEdges = False
             warnPruneEdges = False
             
             startVertex :: forall {a}. Maybe a
             startVertex = Nothing
             
             newPhylogeneticGraph
               | (graphType inGS == Tree) || LG.isTree newGraph = T.multiTraverseFullyLabelGraphReduced inGS inData pruneEdges warnPruneEdges startVertex newGraph
               | (not . LG.cyclic) newGraph && (not . LG.parentInChain) newGraph = T.multiTraverseFullyLabelGraphReduced inGS inData pruneEdges warnPruneEdges startVertex newGraph
               | otherwise = emptyReducedPhylogeneticGraph
         in
         -- return original
         if newPhylogeneticGraph == emptyReducedPhylogeneticGraph then rejoinGB inGS inData intProbAccept sampleAtRandom inTupleList (tail splitGraphList) originalBreakEdge (tail randIntList) edgeToInvade

         -- update tuple list based on new graph
         else
            let updatedTupleList = getLowerGBEdgeCost inTupleList newPhylogeneticGraph -- ((2 * numTaxa) -1)
            in
            rejoinGB inGS inData intProbAccept sampleAtRandom updatedTupleList (tail splitGraphList) originalBreakEdge (tail randIntList) edgeToInvade


      -- return original
      else rejoinGB inGS inData intProbAccept sampleAtRandom inTupleList (tail splitGraphList) originalBreakEdge (tail randIntList) edgeToInvade

-- | mergeTupleLists takes a list of list of tuples and merges them choosing the better each recursive round
mergeTupleLists :: [[(Int, Int, NameBV, NameBV, VertexCost)]] -> [(Int, Int, NameBV, NameBV, VertexCost)] -> [(Int, Int, NameBV, NameBV, VertexCost)]
mergeTupleLists inTupleListList accumList
  | null inTupleListList = accumList
  | null accumList = mergeTupleLists (tail inTupleListList) (head inTupleListList)
  | otherwise = let firstTupleList = head inTupleListList
                    newTupleList = zipWith chooseBetterTuple firstTupleList accumList
                in
                mergeTupleLists (tail inTupleListList) newTupleList

-- | chooseBetterTuple takes two (Int, Int, NameBV, NameBV, VertexCost) and returns better cost
chooseBetterTuple :: (Int, Int, NameBV, NameBV, VertexCost) -> (Int, Int, NameBV, NameBV, VertexCost) -> (Int, Int, NameBV, NameBV, VertexCost)
chooseBetterTuple (aE, aV, aEBV, aVBV, aCost) (_, _, _, _, bCost) = (aE, aV, aEBV, aVBV, min aCost bCost)

-- | makeGraphEdgeTuples take node and edge,cost tuples from a graph and returns a list of tuples of the form
-- (uIndex,vINdex,uBV, vBV, graph cost)
-- this for edge comparisons for Goodman-Bremer and other optimality-type support
makeGraphEdgeTuples :: V.Vector (Int, NameBV) -> VertexCost -> [(Int, Int)] -> [(Int, Int, NameBV, NameBV, VertexCost)]
makeGraphEdgeTuples nodeBVVect graphCost edgeList =
   -- trace ("MET: " <> (show $ V.length nodeBVVect))
   fmap (make5Tuple nodeBVVect graphCost) edgeList
   where make5Tuple :: forall {a} {d} {e}.
                      V.Vector (a, d) -> e -> (Int, Int) -> (Int, Int, d, d, e)
         make5Tuple nv c (a, b) = (a, b, snd (nv V.! a), snd (nv V.! b), c)

-- | getLowerGBEdgeCost take a list of edge tuples of (uIndex,vINdex,uBV, vBV, graph cost) from the graph
-- whose supports are being calculated and a new graph and updates the edge cost (GB value) if that edge
-- is NOT present in the graph taking the minimum of the original GB value and the new graph cost
getLowerGBEdgeCost :: [(Int, Int, NameBV, NameBV, VertexCost)] -> ReducedPhylogeneticGraph -> [(Int, Int, NameBV, NameBV, VertexCost)]
getLowerGBEdgeCost edgeTupleList inGraph =
   if LG.isEmpty (fst5 inGraph) || null edgeTupleList then error "Empty graph or null edge tuple list in getLowerGBEdgeCost"
   else
      let inGraphTupleList = getGraphTupleList inGraph
      in
      fmap (updateEdgeTuple (snd5 inGraph) inGraphTupleList) edgeTupleList


-- | updateEdgeTuple checks is edge is NOT in input graph edge tuple list and if not takes minimum
-- of edge cost GB value and in graph cost, else returns unchanged
updateEdgeTuple :: VertexCost -> [(Int, Int, NameBV, NameBV, VertexCost)] -> (Int, Int, NameBV, NameBV, VertexCost) -> (Int, Int, NameBV, NameBV, VertexCost)
updateEdgeTuple inGraphCost inGraphTupleList (uIndex, vIndex, uBV, vBV, edgeGBValue) =
   let edgeNotFoundCost = getNotFoundCost uBV vBV inGraphCost inGraphTupleList
   in
   if isNothing edgeNotFoundCost then (uIndex, vIndex, uBV, vBV, edgeGBValue)
   else (uIndex, vIndex, uBV, vBV, min edgeGBValue (fromJust edgeNotFoundCost))

-- | getNotFoundCost take a pair of BitVectors (of vertices in graph) from an edge
-- and a list of  (Int, Int, NameBV, NameBV, VertexCost) tuples and returns
-- Nothing is the BVs of the two match (= signifying edge in graph) or
-- Just graph cost if not present for Goodman-Bremer calculations
getNotFoundCost :: NameBV -> NameBV -> VertexCost -> [(Int, Int, NameBV, NameBV, VertexCost)] -> Maybe VertexCost
getNotFoundCost uBV vBV inTupleCost inTupleList =
   if null inTupleList then Just inTupleCost
   else
      let (_, _, uInBV, vInBV, _) = head inTupleList
      in
      if uBV == uInBV && vBV == vInBV then Nothing
      else getNotFoundCost uBV vBV inTupleCost (tail inTupleList)


-- | getTBRSplitGraphs takes a split gaph and the original split edge and
-- returns a list of rerooted subgrahs split graphs suitable for rejoining
-- via SPR-type rejoin each to generate TBR neighborhood
-- much of this is modified from Swap.hs but removing data and delta portions
getTBRSplitGraphs :: GlobalSettings -> SimpleGraph -> LG.LEdge Double -> [SimpleGraph]
getTBRSplitGraphs inGS splitGraph splitEdge =
   if LG.isEmpty splitGraph then error "Empty graph in getTBRSplitGraphs"
   else
      -- get edges in pruned graph and reroot on those edges that are 1) not from original "root" of prune
      -- and 2) not network edges
      let prunedGraphRootNode = (snd3 splitEdge, fromJust $ LG.lab splitGraph $ snd3 splitEdge)
          edgesInPrunedSubGraph = snd $ LG.nodesAndEdgesAfter splitGraph [prunedGraphRootNode]

          nonNetWorkEdgeList = if graphType inGS /= Tree then filter (not . LG.isNetworkLabEdge splitGraph) edgesInPrunedSubGraph
                               else edgesInPrunedSubGraph

          -- original pruned root edges
          prunedRootEdges = LG.out splitGraph $ fst prunedGraphRootNode

          -- edges available for rerooting
          edgeAfterList = nonNetWorkEdgeList L.\\ prunedRootEdges

          -- get edges to add and delete for TBR rerooting
          tbrEdits = fmap (getTBREdits splitGraph prunedGraphRootNode edgesInPrunedSubGraph . LG.toEdge) edgeAfterList

          -- TBR split graph list
          tbrGraphList = fmap (LG.insertDeleteEdges splitGraph) tbrEdits

      in
      splitGraph : tbrGraphList

-- | getTBREdits takes and edge and returns the list of edits to pruned subgraph
-- as a pair of edges to add and those to delete
-- since reroot edge is directed (e,v), edges away from v will have correct
-- orientation. Edges between 'e' and the root will have to be flipped
-- original root edges and reroort edge are deleted and new root and edge spanning orginal root created
-- returns ([add], [delete])
-- modified from function in swap to be more general and operate on SimpleGraphs as are used here
getTBREdits :: (Eq a, Eq b) => LG.Gr a b -> LG.LNode a -> [LG.LEdge b] -> LG.Edge -> ([LG.LEdge b],[LG.Edge])
getTBREdits inGraph prunedGraphRootNode edgesInPrunedSubGraph rerootEdge =
   --trace ("Gettiung TBR Edits for " <> (show rerootEdge)) (
   let prunedGraphRootIndex = fst prunedGraphRootNode
       originalRootEdgeNodes = LG.descendants inGraph prunedGraphRootIndex
       originalRootEdges = LG.out inGraph prunedGraphRootIndex

       -- get path from new root edge fst vertex to orginal root and flip those edges
       closerToPrunedRootEdgeNode = (fst rerootEdge, fromJust $ LG.lab inGraph $ fst rerootEdge)
       (nodesInPath, edgesinPath) = LG.postOrderPathToNode inGraph closerToPrunedRootEdgeNode prunedGraphRootNode

       -- don't want original root edges to be flipped since deleted
       edgesToFlip = edgesinPath L.\\ originalRootEdges
       flippedEdges = fmap LG.flipLEdge edgesToFlip

       -- dummyEdgeLabel so can be type "b"
       dummyEdgeLabel =  thd3 $ head edgesInPrunedSubGraph

       -- new edges on new root position and spanning old root
       -- add in closer vertex to root to make sure direction of edge is correct
       newEdgeOnOldRoot = if snd3 (head originalRootEdges) `elem` (fst rerootEdge : fmap fst nodesInPath) then (snd3 $ head originalRootEdges, snd3 $ last originalRootEdges, dummyEdgeLabel)
                          else (snd3 $ last originalRootEdges, snd3 $ head originalRootEdges, dummyEdgeLabel)
       newRootEdges = [(prunedGraphRootIndex, fst rerootEdge, dummyEdgeLabel),(prunedGraphRootIndex, snd rerootEdge, dummyEdgeLabel)]


   in
   -- original root edge so no change
   if fst rerootEdge `elem` originalRootEdgeNodes &&  snd rerootEdge `elem` originalRootEdgeNodes then ([],[])

   -- rerooted
   else
      -- delete orignal root edges and rerootEdge
      -- add new root edges
      -- and new edge on old root--but need orientation
      -- flip edges from new root to old (delete and add list)
      --trace ("\n\nIn Graph:\n" <> (LG.prettify $ GO.convertDecoratedToSimpleGraph inGraph) <> "\nTBR Edits: " <> (show (rerootEdge, prunedGraphRootIndex, fmap LG.toEdge flippedEdges))
      --   <> "\nEdges to add: " <> (show $ fmap LG.toEdge $ newEdgeOnOldRoot : (flippedEdges <> newRootEdges)) <> "\nEdges to delete: " <> (show $ rerootEdge : (fmap LG.toEdge (edgesToFlip <> originalRootEdges))))
      (newEdgeOnOldRoot : (flippedEdges <> newRootEdges), rerootEdge : fmap LG.toEdge (edgesToFlip <> originalRootEdges))
      -- )
