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
  , makeStaticApprox
  ) where

import           Data.Alphabet
import qualified Data.BitVector.LittleEndian  as BV
import           Data.Bits
import           Data.Char
import qualified Data.List                    as L
import           Data.Maybe
import qualified Data.Vector                  as V
import qualified Data.Vector.Generic          as GV
import qualified Data.Vector.Storable         as SV
import qualified Data.Vector.Unboxed          as UV
import           Data.Word
import           Debug.Trace
import           Foreign.C.Types              (CUInt)
import           GeneralUtilities
import qualified GraphOptimization.Traversals as T
import qualified Graphs.GraphOperations       as GO
import qualified Input.Reorganize             as R
import qualified ParallelUtilities            as PU
import           Text.Read
import           Types.Types
import qualified Utilities.LocalGraph         as LG
import qualified Utilities.Utilities          as U
-- import qualified Input.DataTransformation    as TRANS
import qualified Commands.Verify              as VER
import qualified Data.Bits                    as B
import qualified Data.Char                    as C
import qualified Data.Text.Lazy               as TL
import qualified Input.BitPack                as BP



-- | transform changes aspects of data sande settings during execution
-- as opposed to Set with all happens at begginign of program execution
transform :: [Argument] -> GlobalSettings -> ProcessedData -> ProcessedData -> Int -> [ReducedPhylogeneticGraph] -> (GlobalSettings, ProcessedData, ProcessedData, [ReducedPhylogeneticGraph])
transform inArgs inGS origData inData rSeed inGraphList =
   let fstArgList = fmap (fmap toLower . fst) inArgs
       sndArgList = fmap (fmap toLower . snd) inArgs
       lcArgList = zip fstArgList sndArgList
       checkCommandList = checkCommandArgs "transform" fstArgList VER.transformArgList
   in
   -- check for valid command options
   if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'transform': " <> show inArgs)
   else
        let displayBlock = filter ((=="displaytrees").fst) lcArgList
            numDisplayTrees
               | length displayBlock > 1 =
                  errorWithoutStackTrace ("Multiple displayTree number specifications in transform--can have only one: " <> show inArgs)
               | null displayBlock = Just 10
               | null (snd $ head displayBlock) = Just 10
               | otherwise = readMaybe (snd $ head displayBlock) :: Maybe Int

            toTree = any ((=="totree").fst) lcArgList
            toSoftWired = any ((=="tosoftwired").fst) lcArgList
            toHardWired = any ((=="tohardwired").fst) lcArgList
            toStaticApprox = any ((=="staticapprox").fst) lcArgList
            toDynamic = any ((=="dynamic").fst) lcArgList
            atRandom = any ((=="atrandom").fst) lcArgList
            chooseFirst = any ((=="first").fst) lcArgList
            reWeight =  any ((=="weight").fst) lcArgList
            changeEpsilon = any ((=="dynamicepsilon").fst) lcArgList
            reRoot = any ((=="outgroup").fst) lcArgList
            changeGraphsSteepest = any ((=="graphssteepest").fst) lcArgList
            changeSoftwiredMethod = any ((=="softwiredmethod").fst) lcArgList
            changeGraphFactor = any ((=="graphfactor").fst) lcArgList
            changeCompressionResolutions = any ((=="compressresolutions").fst) lcArgList
            changeMultiTraverse = any ((=="multitraverse").fst) lcArgList
            changeUnionThreshold = any ((=="jointhreshold").fst) lcArgList


            reweightBlock = filter ((=="weight").fst) lcArgList
            weightValue
               | length reweightBlock > 1 =
                  errorWithoutStackTrace ("Multiple weight specifications in transform--can have only one: " <> show inArgs)
               | null reweightBlock = Just 1.0
               | null (snd $ head reweightBlock) = Just 1
               | otherwise = readMaybe (snd $ head reweightBlock) :: Maybe Double

            changeCompressionBlock = filter ((=="compressresolutions").fst) lcArgList
            compressionValue
               | length changeCompressionBlock > 1 =
                  errorWithoutStackTrace ("Multiple compressResolutions specifications in transform--can have only one: " <> show inArgs)
               | null changeCompressionBlock = Just $ fmap toLower $ show $ compressResolutions inGS
               | null (snd $ head changeCompressionBlock) = Just $ fmap toLower $ show $ compressResolutions inGS
               | otherwise = readMaybe (show $ snd $ head changeCompressionBlock) :: Maybe String

            changeEpsilonBlock = filter ((=="dynamicepsilon").fst) lcArgList
            epsilonValue
               | length changeEpsilonBlock > 1 =
                  errorWithoutStackTrace ("Multiple dynamicEpsilon specifications in transform--can have only one: " <> show inArgs)
               | null changeEpsilonBlock = Just $ dynamicEpsilon inGS
               | null (snd $ head changeEpsilonBlock) = Just $ dynamicEpsilon inGS
               | otherwise = readMaybe (snd $ head changeEpsilonBlock) :: Maybe Double

            changeGraphFactorBlock = filter ((=="graphfactor").fst) lcArgList
            newGraphFactor
               | length changeGraphFactorBlock > 1 =
                  errorWithoutStackTrace ("Multiple graphFactor specifications in transform--can have only one: " <> show inArgs)
               | null changeGraphFactorBlock = Just $ fmap toLower $ show $ graphFactor inGS
               | null (snd $ head changeGraphFactorBlock) = Just $ fmap toLower $ show $ graphFactor inGS
               | otherwise = readMaybe (show $ snd $ head changeGraphFactorBlock) :: Maybe String

            changeGraphsSteepestBlock = filter ((=="graphssteepest").fst) lcArgList
            newGraphsSteepest
               | length changeGraphsSteepestBlock > 1 =
                  errorWithoutStackTrace ("Multiple graphsSteepest specifications in transform--can have only one: " <> show inArgs)
               | null changeGraphsSteepestBlock = Just $ graphsSteepest inGS
               | null (snd $ head changeGraphsSteepestBlock) = Just $ graphsSteepest inGS
               | otherwise = readMaybe (snd $ head changeGraphsSteepestBlock) :: Maybe Int

            changeMultiTraverseBlock = filter ((=="multitraverse").fst) lcArgList
            multiTraverseValue
               | length changeMultiTraverseBlock > 1 =
                  errorWithoutStackTrace ("Multiple multiTraverse specifications in transform--can have only one: " <> show inArgs)
               | null changeMultiTraverseBlock = Just $ fmap toLower $ show $ multiTraverseCharacters inGS
               | null (snd $ head changeMultiTraverseBlock) = Just $ fmap toLower $ show $ multiTraverseCharacters inGS
               | otherwise = readMaybe (show $ snd $ head changeMultiTraverseBlock) :: Maybe String

            changeSoftwiredMethodBlock = filter ((=="softwiredmethod").fst) lcArgList
            newSoftwiredMethod
               | length changeSoftwiredMethodBlock > 1 =
                  errorWithoutStackTrace ("Multiple softwiredMethod specifications in transform--can have only one: " <> show inArgs)
               | null changeSoftwiredMethodBlock = Just $ fmap toLower $ show $ softWiredMethod inGS
               | null (snd $ head changeSoftwiredMethodBlock) = Just $ fmap toLower $ show $ softWiredMethod inGS
               | otherwise = readMaybe (show $ snd $ head changeSoftwiredMethodBlock) :: Maybe String

            reRootBlock = filter ((=="outgroup").fst) lcArgList
            outgroupValue
               | length reRootBlock > 1 =
                  errorWithoutStackTrace ("Multiple outgroup specifications in transform--can have only one: " <> show inArgs)
               | null reRootBlock = Just $ outGroupName inGS
               | null (snd $ head reRootBlock) = Just $ outGroupName inGS
               | otherwise = readMaybe (snd $ head reRootBlock) :: Maybe TL.Text

            changeUnionBlock = filter ((=="jointhreshold").fst) lcArgList
            unionValue
               | length changeUnionBlock > 1 =
                  errorWithoutStackTrace ("Multiple joinThreshold specifications in transform--can have only one: " <> show inArgs)
               | null changeUnionBlock = Just $ unionThreshold inGS
               | null (snd $ head changeUnionBlock) = Just $ unionThreshold inGS
               | otherwise = readMaybe (snd $ head changeUnionBlock) :: Maybe Double

            nameList = fmap TL.pack $ fmap (filter (/= '"')) $ fmap snd $ filter ((=="name").fst) lcArgList
            charTypeList = fmap snd $ filter ((=="type").fst) lcArgList
        

        in
        if (length $ filter (== True) [toTree, toSoftWired, toHardWired]) > 1 then
            errorWithoutStackTrace ("Multiple graph transform commands--can only have one : " <> (show inArgs))
        else if toStaticApprox && toDynamic then
            errorWithoutStackTrace ("Multiple staticApprox/Dynamic transform commands--can only have one : " <> (show inArgs))
        else if atRandom && chooseFirst then
            errorWithoutStackTrace ("Multiple display tree choice commands in transform (first, atRandom)--can only have one : " <> (show inArgs))
        else if (toTree || toSoftWired || toHardWired) && (toDynamic || toStaticApprox) then
            errorWithoutStackTrace ("Multiple transform operations in transform (e.g. toTree, staticApprox)--can only have one at a time: " <> (show inArgs))
        else
            let pruneEdges = False
                warnPruneEdges = False
                startVertex = Nothing
            in
            -- transform nets to tree
            if toTree then
               -- already Tree return
               if (graphType inGS == Tree) then (inGS, origData, inData, inGraphList)
               else
                  let newGS = inGS {graphType = Tree}

                      -- generate and return display trees-- displayTreNUm / graph
                      contractIn1Out1Nodes = True
                      displayGraphList = if chooseFirst then fmap (take (fromJust numDisplayTrees) . (LG.generateDisplayTrees) contractIn1Out1Nodes) (fmap fst5 inGraphList)
                                         else fmap (LG.generateDisplayTreesRandom rSeed (fromJust numDisplayTrees)) (fmap fst5 inGraphList)

                      -- prob not required
                      displayGraphs = fmap GO.ladderizeGraph $ fmap GO.renameSimpleGraphNodes (concat displayGraphList)

                      -- reoptimize as Trees
                      newPhylogeneticGraphList = PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced newGS inData pruneEdges warnPruneEdges startVertex) displayGraphs -- `using` PU.myParListChunkRDS
                  in
                  (newGS, origData, inData, newPhylogeneticGraphList)

            -- transform to softwired
            else if toSoftWired then
               if (graphType inGS == SoftWired) then (inGS, origData, inData, inGraphList)
               else
                  let newGS = inGS {graphType = SoftWired}
                      newPhylogeneticGraphList = PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced newGS inData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList)  -- `using` PU.myParListChunkRDS
                  in
                  (newGS, origData, inData, newPhylogeneticGraphList)

            -- transform to hardwired
            else if toHardWired then
               if (graphType inGS == HardWired) then (inGS, origData, inData, inGraphList)
               else
                  let newGS = inGS {graphType = HardWired, graphFactor = NoNetworkPenalty}

                      newPhylogeneticGraphList = PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced newGS inData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList)  -- `using` PU.myParListChunkRDS
                  in
                  trace ("Changing GraphFactor to NoNetworkPenalty for HardWired graphs")
                  (newGS, origData, inData, newPhylogeneticGraphList)

            -- roll back to dynamic data from static approx
            else if toDynamic then
               let newPhylogeneticGraphList = PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced inGS origData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList) -- `using` PU.myParListChunkRDS
               in
               trace ("Transforming data to dynamic: " <> (show $ minimum $ fmap snd5 inGraphList) <> " -> " <> (show $ minimum $ fmap snd5 newPhylogeneticGraphList))
               (inGS, origData, origData, newPhylogeneticGraphList)

            -- transform to static approx--using first Tree
            else if toStaticApprox then
               let newData = makeStaticApprox inGS False inData (head $ L.sortOn snd5 inGraphList)
                   newPhylogeneticGraphList = PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced inGS newData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList) -- `using` PU.myParListChunkRDS

               in
               trace ("Transforming data to staticApprox: " <> (show $ minimum $ fmap snd5 inGraphList) <> " -> " <> (show $ minimum $ fmap snd5 newPhylogeneticGraphList))
               (inGS, origData, newData, newPhylogeneticGraphList)

            -- change weight values in charInfo and reoptimize
            -- reweights both origData and inData so weighting doens't get undone by static approc to and from transfomrations
            else if reWeight then
               let newOrigData = reWeightData (fromJust weightValue) charTypeList nameList origData
                   newData = reWeightData (fromJust weightValue) charTypeList nameList inData
                   newPhylogeneticGraphList = PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced inGS newData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList) -- `using` PU.myParListChunkRDS
               in
               if isNothing weightValue then errorWithoutStackTrace ("Reweight value is not specified correcty. Must be a double (e.g. 1.2): " <> (show  (snd $ head reweightBlock)))
               else
                  trace ("Reweighting types " <> (show charTypeList) <> " and/or characters " <> (L.intercalate ", " $ fmap TL.unpack nameList) <> " to " <> (show $ fromJust weightValue)
                  <> "\n\tReoptimizing graphs")
                  (inGS, newOrigData, newData, newPhylogeneticGraphList)

            -- changes the softwired optimization algorithm--this really for experimental use
            else if changeCompressionResolutions then
               if isNothing compressionValue then errorWithoutStackTrace ("CompressResolutions value is not specified correcty. Must be either 'True' or 'False': " <> (show (snd $ head changeCompressionBlock)))
               else
                  let newMethod = if fromJust compressionValue == "true" then True
                                  else if fromJust compressionValue == "false" then False
                                  else errorWithoutStackTrace ("CompressResolutions value is not specified correcty. Must be either 'True' or 'False': " <> (show (snd $ head changeSoftwiredMethodBlock)))
                      newPhylogeneticGraphList = PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced (inGS  {compressResolutions = newMethod}) origData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList) -- `using` PU.myParListChunkRDS
                  in
                  if newMethod /=  compressResolutions inGS then
                     trace ("Changing compressResolutions method to " <> (show newMethod))
                     (inGS {compressResolutions = newMethod}, origData, inData, newPhylogeneticGraphList)
                  else (inGS {compressResolutions = newMethod}, origData, inData, inGraphList)

            -- changes dynamicEpsilon error check factor
            else if changeEpsilon then
               if isNothing epsilonValue then errorWithoutStackTrace ("DynamicEpsilon value is not specified correcty. Must be a double (e.g. 0.02): " <> (show (snd $ head changeEpsilonBlock)))
               else
                  trace ("Changing dynamicEpsilon factor to " <> (show $ fromJust epsilonValue))
                  (inGS {dynamicEpsilon = 1.0 + ((fromJust epsilonValue) * (fractionDynamic inGS))}, origData, inData, inGraphList)

            -- changes the softwired optimization algorithm--this really for experimental use
            else if changeGraphFactor then
               if isNothing newGraphFactor then errorWithoutStackTrace ("GraphFactor value is not specified correcty. Must be either 'NoPenalty', 'W15'. 'W23', or 'PMDL': " <> (show (snd $ head changeGraphFactorBlock)))
               else
                  let newMethod = if fromJust newGraphFactor == "nopenalty" then NoNetworkPenalty
                                  else if fromJust newGraphFactor == "w15" then Wheeler2015Network
                                  else if fromJust newGraphFactor == "w23" then Wheeler2023Network
                                  else if fromJust newGraphFactor == "pmdl" then PMDLGraph
                                  else errorWithoutStackTrace ("GraphFactor value is not specified correcty. Must be either 'NoPenalty', 'W15'. 'W23', or 'PMDL': " <> (show (snd $ head changeGraphFactorBlock)))
                      newPhylogeneticGraphList = PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced (inGS  {graphFactor = newMethod}) origData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList) -- `using` PU.myParListChunkRDS
                  in
                  if newMethod /=  graphFactor inGS then
                     trace ("Changing graphFactor method to " <> (show newMethod))
                     (inGS {graphFactor = newMethod}, origData, inData, newPhylogeneticGraphList)
                  else (inGS {graphFactor = newMethod}, origData, inData, inGraphList)

            -- changes graphsSteepest -- maximum number of graphs evaluated in paralell at each "steepest" phse in swpa dn netadd/delete
            else if changeGraphsSteepest then
               if isNothing newGraphsSteepest then errorWithoutStackTrace ("GraphsSteepest value is not specified correcty. Must be an Integer (e.g. 5): " <> (show (snd $ head changeGraphsSteepestBlock)))
               else
                  trace ("Changing GraphsSteepest factor to " <> (show $ fromJust newGraphsSteepest))
                  (inGS {graphsSteepest = fromJust newGraphsSteepest}, origData, inData, inGraphList)

            -- changes the multiTraverse behavior for all graphs
            else if changeMultiTraverse then
               if isNothing multiTraverseValue then errorWithoutStackTrace ("MultiTraverse value is not specified correcty. Must be either 'True' or 'False': " <> (show (snd $ head changeMultiTraverseBlock)))
               else
                  let newMethod = if fromJust multiTraverseValue == "true" then True
                                  else if fromJust multiTraverseValue == "false" then False
                                  else errorWithoutStackTrace ("MultiTraverse value is not specified correcty. Must be either 'True' or 'False': " <> (show (snd $ head changeMultiTraverseBlock)))
                      newPhylogeneticGraphList = PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced (inGS  {multiTraverseCharacters = newMethod}) origData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList) -- `using` PU.myParListChunkRDS
                  in
                  if newMethod /=  multiTraverseCharacters inGS then
                     let lengthChangeString = if null inGraphList then ""
                                              else (":" <> (show $ minimum $ fmap snd5 inGraphList) <> " -> " <> (show $ minimum $ fmap snd5 newPhylogeneticGraphList))
                     in
                     trace ("Changing multiTraverse to " <> (show newMethod) <> lengthChangeString)
                     (inGS {multiTraverseCharacters = newMethod}, origData, inData, newPhylogeneticGraphList)
                  else (inGS {multiTraverseCharacters = newMethod}, origData, inData, inGraphList)

            -- changes the softwired optimization algorithm--this really for experimental use
            else if changeSoftwiredMethod then
               if isNothing newSoftwiredMethod then errorWithoutStackTrace ("SoftwiredMethod value is not specified correcty. Must be either 'Exhaustive' or 'ResolutionCache': " <> (show (snd $ head changeSoftwiredMethodBlock)))
               else
                  let newMethod = if fromJust newSoftwiredMethod == "naive" then Naive
                                  else if fromJust newSoftwiredMethod == "exhaustive" then Naive
                                  else if fromJust newSoftwiredMethod == "resolutioncache" then ResolutionCache
                                  else errorWithoutStackTrace ("SoftwiredMethod value is not specified correcty. Must be either 'Naive' or 'ResolutionCache': " <> (show (snd $ head changeSoftwiredMethodBlock)))
                      newPhylogeneticGraphList = PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced (inGS  {softWiredMethod = newMethod}) origData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList) -- `using` PU.myParListChunkRDS
                      newMethodString = if newMethod == ResolutionCache then "ResolutionCache"
                                        else "Exhaustive"
                  in
                  if newMethod /=  softWiredMethod inGS then
                     trace ("Changing softwired optimization method to " <> newMethodString)
                     (inGS {softWiredMethod = newMethod}, origData, inData, newPhylogeneticGraphList)
                  else (inGS {softWiredMethod = newMethod}, origData, inData, inGraphList)

            -- changes outgroup
            else if reRoot then
               if isNothing outgroupValue then errorWithoutStackTrace ("Outgroup is not specified correctly. Must be a string (e.g. \"Name\"): " <> (snd $ head reRootBlock))
               else
                  let newOutgroupName = TL.filter (/= '"') $ fromJust outgroupValue
                      newOutgroupIndex =  V.elemIndex newOutgroupName (fst3 origData)
                      newPhylogeneticGraphList = PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced inGS origData pruneEdges warnPruneEdges startVertex) (fmap (LG.rerootTree (fromJust newOutgroupIndex)) $ fmap fst5 inGraphList)
                  in
                  if isNothing newOutgroupIndex then errorWithoutStackTrace ("Outgoup name not found: " <> (snd $ head reRootBlock))
                  else
                     trace ("Changing outgroup to " <> (TL.unpack newOutgroupName))
                     (inGS {outgroupIndex = fromJust newOutgroupIndex, outGroupName = newOutgroupName}, origData, inData, newPhylogeneticGraphList)

            -- changes unionThreshold error check factor
            else if changeUnionThreshold then
               if isNothing unionValue then errorWithoutStackTrace ("UninThreshold value is not specified correcty. Must be a double (e.g. 1.17): " <> (show (snd $ head changeUnionBlock)))
               else
                  trace ("Changing uninoTHreshold factor to " <> (show $ fromJust unionValue))
                  (inGS {dynamicEpsilon = (fromJust unionValue)}, origData, inData, inGraphList)

            -- modify the use of Network Add heurisitcs in network optimization
            else if (fst $ head lcArgList) == "usenetaddheuristic"  then
            let localCriterion
                  | ((snd $ head lcArgList) == "true") = True
                  | ((snd $ head lcArgList) == "false") = False
                  | otherwise = errorWithoutStackTrace ("Error in 'transform' command. UseNetAddHeuristic '" <> (snd $ head lcArgList) <> "' is not 'true' or 'false'")
            in
            trace ("UseNetAddHeuristic set to " <> (snd $ head lcArgList))
            (inGS {useNetAddHeuristic = localCriterion}, origData, inData, inGraphList)

            else error ("Transform type not implemented/recognized" <> (show inArgs))

-- | reWeightData sets weights to new values based on
reWeightData :: Double -> [String] -> [NameText] -> ProcessedData -> ProcessedData
reWeightData weightValue charTypeStringList charNameList (inName, inNameBV, inBlockDataV) =
   let charTypeList = concatMap stringToType charTypeStringList
       newBlockData = fmap (reweightBlockData weightValue charTypeList charNameList) inBlockDataV
   in
   (inName, inNameBV, newBlockData)

-- |  stringToType takes  String and returns typelist
stringToType :: String -> [CharType]
stringToType inString =
   if null inString then []
   else
      let inVal = fmap C.toLower inString
          typeList = if inVal == "all" then exactCharacterTypes <> sequenceCharacterTypes
                     else if inVal == "prealigned" then prealignedCharacterTypes
                     else if inVal `elem` ["nonexact", "dynamic"] then nonExactCharacterTypes
                     else if inVal == "nonadditive" then [NonAdd, Packed2, Packed4, Packed5, Packed8, Packed64]
                     else if inVal == "additive" then [Add]
                     else if inVal == "matrix" then [Matrix]
                     else if inVal == "sequence" then sequenceCharacterTypes
                     else if inVal == "packed" then   packedNonAddTypes
                     else if inVal == "packed2" then [Packed2]
                     else if inVal == "packed4" then [Packed4]
                     else if inVal == "packed5" then [Packed5]
                     else if inVal == "packed8" then [Packed8]
                     else if inVal == "packed64" then [Packed64]
                     else if inVal `elem` ["static", "exact", "qualitative"] then exactCharacterTypes

                     else errorWithoutStackTrace ("Error in transform : Unrecognized character type '" <> inString <> "'")
      in
      typeList

-- | reweightBlockData applies new weight to catagories of data
reweightBlockData :: Double -> [CharType] -> [NameText] -> BlockData -> BlockData
reweightBlockData  weightValue charTypeList charNameList (blockName, blockData, charInfoV) =
   let newCharacterInfoV = fmap (reweightCharacterData weightValue charTypeList charNameList) charInfoV
   in
   (blockName, blockData, newCharacterInfoV)

-- | reweightCharacterData changes weight in charInfo based on type or name
reweightCharacterData ::  Double -> [CharType] -> [NameText] -> CharInfo -> CharInfo
reweightCharacterData weightValue charTypeList charNameList charInfo =
   let wildCardMatchCharName = filter (== True) $ fmap (textMatchWildcards (name charInfo)) charNameList
   in
   -- trace ("RWC Wildcards: " <> (show $ fmap (textMatchWildcards (name charInfo)) charNameList)) (
   if null wildCardMatchCharName && (charType charInfo) `notElem` charTypeList then
      -- trace ("RWC not : " <> (show $ name charInfo) <> " of " <> (show charNameList) <> " " <> (show $ charType charInfo) <> " of " <> (show charTypeList))
      charInfo
   else
      -- trace ("RWC: " <> (show $ name charInfo) <> " " <> (show $ charType charInfo))
      charInfo {weight = weightValue}
   -- )

-- | makeStaticApprox takes ProcessedData and returns static approx (implied alignment recoded) ProcessedData
-- if Tree take SA fields and recode appropriatrely given cost regeme of character
-- if Softwired--use display trees for SA
-- if hardWired--convert to softwired and use display trees for SA
-- since for heuristic searcing--uses additive weight for sequences and simple cost matrices, otherwise
-- matrix characters
makeStaticApprox :: GlobalSettings -> Bool -> ProcessedData -> ReducedPhylogeneticGraph -> ProcessedData
makeStaticApprox inGS leavePrealigned inData inGraph =
   if LG.isEmpty (fst5 inGraph) then error "Empty graph in makeStaticApprox"

   -- tree type
   else if graphType inGS == Tree then
      let decGraph = thd5 inGraph
          (nameV, nameBVV, blockDataV) = inData

          -- do each block in turn pulling and transforming data from inGraph
          newBlockDataV = PU.seqParMap (parStrategy $ strictParStrat inGS) (pullGraphBlockDataAndTransform leavePrealigned decGraph inData) [0..(length blockDataV - 1)] -- `using` PU.myParListChunkRDS

          -- convert prealigned to non-additive if all 1's tcm

          -- remove constants from new prealigned
          newProcessedData  = R.removeConstantCharactersPrealigned (nameV, nameBVV, V.fromList newBlockDataV)

          -- bit pack any new non-additive characters
          newProcessedData' = BP.packNonAdditiveData inGS newProcessedData
      in
      -- trace ("MSA:" <> (show (fmap (V.length . thd3) blockDataV, fmap (V.length . thd3) newBlockDataV)))
      if leavePrealigned then (nameV, nameBVV, V.fromList newBlockDataV)
      else newProcessedData'

   else trace ("Static Approx not yet implemented for graph type :" <> (show $ graphType inGS) <> " skipping") inData


-- | pullGraphBlockDataAndTransform takes a DecoratedGrpah and block index and pulls
-- the character data of the block and transforms the leaf data by using implied alignment
-- feilds for dynamic characters
pullGraphBlockDataAndTransform :: Bool -> DecoratedGraph -> ProcessedData -> Int -> BlockData
pullGraphBlockDataAndTransform leavePrealigned inDecGraph  (_, _, blockCharInfoV) blockIndex =
   let (_, leafVerts, _, _) = LG.splitVertexList inDecGraph
       (_, leafLabelList) = unzip leafVerts
       leafBlockData = fmap (V.! blockIndex) (fmap vertData leafLabelList)

       -- new recoded data-- need to filter out constant chars after recoding
       -- nedd character legnth for missing values
       charLengthV = V.zipWith U.getMaxCharacterLength (thd3 $ blockCharInfoV V.! blockIndex) (V.fromList $ fmap V.toList leafBlockData)

       (transformedLeafBlockData, transformedBlockInfo) = unzip $ fmap (transformData leavePrealigned (thd3 $ blockCharInfoV V.! blockIndex) charLengthV) leafBlockData
   in
   -- trace ("PGDT: " <> show charLengthV)
   (fst3 $ blockCharInfoV V.! blockIndex, V.fromList transformedLeafBlockData, head transformedBlockInfo)


-- | transformData takes original Character info and character data and transforms to static if dynamic noting chracter type
transformData :: Bool -> V.Vector CharInfo -> V.Vector Int -> V.Vector CharacterData -> (V.Vector CharacterData, V.Vector CharInfo)
transformData leavePrealigned inCharInfoV inCharLengthV inCharDataV  =
   if V.null inCharInfoV then
      (V.empty, V.empty)
   else
      let (outCharDataV, outCharInfoV) = V.unzip $ V.zipWith3 (transformCharacter leavePrealigned) inCharDataV inCharInfoV inCharLengthV
      in
      (outCharDataV, outCharInfoV)

-- transformCharacter takes a single character info and character and returns IA if dynamic as is if not
-- checks if all gaps with the GV.filter.  If all gaps--it means the sequence char was missing and
-- implied alignment produced all gaps.  The passing of character length is not necessary when changed missing seq to empty
-- character--but leaving in case change back to []
-- "nonAddGap" not currently implemented
transformCharacter :: Bool -> CharacterData -> CharInfo -> Int -> (CharacterData, CharInfo)
transformCharacter leavePrealigned inCharData inCharInfo charLength =
   let inCharType = charType inCharInfo
       inCostMatrix = costMatrix inCharInfo
       alphSize = length $ alphabet inCharInfo

       -- determine if matrix is all same costs => nonadditive
       --                        all same except fort single indel costs => non add with gap binary chars
       --                        not either => matrix char
       (inCostMatrixType, gapCost) = R.getRecodingType inCostMatrix

   in
   -- trace ("TC:" <> (show alphSize) <> " " <> (show $ alphabet inCharInfo)) (
   --trace ("TC:" <> (show charLength) <> " " <> (show (GV.length $ snd3 $ slimAlignment inCharData, GV.length $ snd3 $ wideAlignment inCharData, GV.length $ snd3 $ hugeAlignment inCharData))) (
   if inCharType `elem` exactCharacterTypes then (inCharData, inCharInfo)

   else if inCharType `elem` prealignedCharacterTypes then (inCharData, inCharInfo)

   else
      -- trace ("TC: " <> inCostMatrixType) (
      -- different types--vector wrangling
      -- missing data fields set if no implied alignment ie missing data
      if inCharType `elem` [SlimSeq, NucSeq] then
         let gapChar = (0 :: CUInt) `setBit` fromEnum gapIndex
             impliedAlignChar = if (not . GV.null $ GV.filter (/= gapChar) $ snd3 $ slimAlignment inCharData) then slimAlignment inCharData
                                else
                                  let missingElement = SV.replicate charLength $ B.complement (0 :: CUInt) -- TRANS.setMissingBits (0 :: CUInt) 0 alphSize
                                  in
                                  (missingElement, missingElement, missingElement)

             newPrelimBV = R.convert2BV 32 impliedAlignChar
             newPrelimBVGaps = addGaps2BV gapCost newPrelimBV
         in
         -- trace ("TC-Slim:" <> (show $ GV.length $ snd3 $ slimAlignment inCharData) <> " " <> (show $ snd3 $ impliedAlignChar)) (

         if leavePrealigned then
            (inCharData {alignedSlimPrelim = impliedAlignChar}, inCharInfo {charType =  AlignedSlim})

         else if inCostMatrixType == "nonAdd" then
            (inCharData {stateBVPrelim = newPrelimBV}, inCharInfo {charType = NonAdd})

         else if inCostMatrixType == "nonAddGap" then
            (inCharData {stateBVPrelim = newPrelimBVGaps}, inCharInfo {charType = NonAdd})

         else -- matrix recoding
            (inCharData {alignedSlimPrelim = impliedAlignChar}, inCharInfo {charType =  AlignedSlim})

      else if inCharType `elem` [WideSeq, AminoSeq] then
         let gapChar = (0 :: Word64) `setBit` fromEnum gapIndex
             impliedAlignChar = if (not . GV.null $ GV.filter (/= gapChar) $ snd3 $ wideAlignment inCharData)  then wideAlignment inCharData
                                else
                                  let missingElement = UV.replicate charLength $ B.complement (0 :: Word64) -- TRANS.setMissingBits (0 :: Word64) 0 alphSize
                                  in (missingElement, missingElement, missingElement)

             newPrelimBV = R.convert2BV 64 impliedAlignChar
             newPrelimBVGaps = addGaps2BV gapCost newPrelimBV
         in
         if leavePrealigned then
            (inCharData {alignedWidePrelim = impliedAlignChar}, inCharInfo {charType =  AlignedWide})

         else if inCostMatrixType == "nonAdd" then
            (inCharData {stateBVPrelim = newPrelimBV}, inCharInfo {charType = NonAdd})

         else if inCostMatrixType == "nonAddGap" then
            (inCharData {stateBVPrelim = newPrelimBVGaps}, inCharInfo {charType = NonAdd})

         else -- matrix recoding
            (inCharData {alignedWidePrelim = impliedAlignChar}, inCharInfo {charType =  AlignedWide})

      else if inCharType == HugeSeq then
         let gapChar = (BV.fromBits $ replicate alphSize False) `setBit` fromEnum gapIndex
             impliedAlignChar = if (not . GV.null $ GV.filter (/= gapChar) $ snd3 $ hugeAlignment inCharData) then hugeAlignment inCharData
                                else
                                 let missingElement = V.replicate alphSize $ (BV.fromBits $ replicate alphSize True)
                                 in (missingElement, missingElement, missingElement)

             newPrelimBV = impliedAlignChar
             newPrelimBVGaps = addGaps2BV gapCost newPrelimBV
         in
         if leavePrealigned then
            (inCharData {alignedHugePrelim = impliedAlignChar}, inCharInfo {charType =  AlignedHuge})

         else if inCostMatrixType == "nonAdd" then
            (inCharData {stateBVPrelim = newPrelimBV}, inCharInfo {charType = NonAdd})

         else if inCostMatrixType == "nonAddGap" then
            (inCharData {stateBVPrelim = newPrelimBVGaps}, inCharInfo {charType = NonAdd})

         else -- matrix recoding
            (inCharData {alignedHugePrelim = impliedAlignChar}, inCharInfo {charType =  AlignedHuge})

      else
         error ("Unrecognized character type in transformCharacter: " <> (show inCharType))
      -- )

-- | addGaps2BV adds gap characters 0 = nonGap, 1 = Gap to Vector
-- of states to non-additive charcaters for static approx.  gapCost - 1 characters are added
-- sems wasteful, but comctant filtered out and recoded later when non-add/add charsa re optimized and bitpacked
-- since this only for leaves assume inM good for all
addGaps2BV :: Int -> (V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector) -> (V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector)
addGaps2BV gapCost (_, inM, _) =
   -- trace ("AG2BV: " <> (show inM)) (
   let gapChar = BV.fromNumber (BV.dimension $ V.head inM) (1 :: Int)
       noGap = L.replicate (gapCost - 1) $ BV.fromNumber (BV.dimension $ V.head inM) (1 :: Int)
       hasGap =  L.replicate (gapCost - 1) $ BV.fromNumber (BV.dimension $ V.head inM) (2 :: Int)
       gapCharV = createGapChars inM gapChar [] noGap hasGap
       outM = inM <> gapCharV
   in
   (outM, outM, outM)
   -- )

-- | createGapChars takes a vector of bitvector coded states and checks if first states == 1 (= gap)
-- if so a number based on gap cost are created.. Will create n * original klength so need to
-- filter out constant characters later
createGapChars :: V.Vector BV.BitVector -> BV.BitVector -> [BV.BitVector] -> [BV.BitVector] -> [BV.BitVector] -> V.Vector BV.BitVector
createGapChars origBVV gapCharacter newCharL noGapL hasGapL =
   if V.null origBVV then V.fromList newCharL
   else
      if V.head origBVV == gapCharacter then createGapChars (V.tail origBVV) gapCharacter (hasGapL <> newCharL) noGapL hasGapL
      else createGapChars (V.tail origBVV) gapCharacter (noGapL <> newCharL) noGapL hasGapL

