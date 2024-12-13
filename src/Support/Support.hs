{- |
Module containing support functions
-}
module Support.Support (
    supportGraph,
) where

import Commands.Verify qualified as VER
import Control.Monad (when)
-- import Control.Monad.IO.Class (MonadIO (..))
import Control.Monad.Random.Class
import Data.Char
import Data.List qualified as L
import Data.Maybe
import Data.Vector qualified as V
import Data.Vector.Generic qualified as GV
import Data.Vector.Primitive (convert)
import Data.Vector.Unboxed qualified as UV
import GeneralUtilities
import GraphOptimization.PostOrderSoftWiredFunctions qualified as POSW
import GraphOptimization.Traversals qualified as T
import Graphs.GraphOperations qualified as GO
import PHANE.Evaluation
-- import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
-- import PHANE.Evaluation.Verbosity (Verbosity (..))
import Reconciliation.ReconcileGraphs qualified as REC
import Search.Build qualified as B
import Search.NetworkAddDelete qualified as N
import Search.Refinement qualified as R
import Text.Read
import Types.Types
-- import Utilities.Distances qualified as DD
import Utilities.LocalGraph qualified as LG
import Utilities.Utilities qualified as U
--import Debug.Trace

-- | driver for overall support
supportGraph ∷ [Argument] → GlobalSettings → ProcessedData → [ReducedPhylogeneticGraph] → PhyG [ReducedPhylogeneticGraph]
supportGraph inArgs inGS inData inGraphList =
    if null inGraphList
        then error "No graphs input to calculate support"
        else
            let fstArgList = fmap (fmap toLower . fst) inArgs
                sndArgList = fmap (fmap toLower . snd) inArgs
                lcArgList = zip fstArgList sndArgList
                checkCommandList = checkCommandArgs "support" fstArgList VER.supportArgList
            in  -- check for valid command options
                if not checkCommandList
                    then errorWithoutStackTrace ("Unrecognized command in 'support': " <> show inArgs)
                    else
                        let supportMeasure
                                | any ((== "bootstrap") . fst) lcArgList = Bootstrap
                                | any ((== "jackknife") . fst) lcArgList = Jackknife
                                | any ((== "goodmanbremer") . fst) lcArgList = GoodmanBremer
                                | otherwise = GoodmanBremer

                            useSPR = any ((== "spr") . fst) lcArgList
                            useTBR = any ((== "tbr") . fst) lcArgList

                            onlyBuild = any ((== "buildonly") . fst) lcArgList

                            jackList = filter ((== "jackknife") . fst) lcArgList
                            jackFreq'
                                | length jackList > 1 =
                                    errorWithoutStackTrace
                                        ( "Multiple jackknife sampling frequency specifications in support command--can have only one (e.g. jackknife:0.62): "
                                            <> show inArgs
                                        )
                                | null jackList = Just 0.6321 -- 1- 1/e
                                | null (snd $ head jackList) = Just 0.6321
                                | otherwise = readMaybe (snd $ head jackList) ∷ Maybe Double

                            replicateList = filter ((== "replicates") . fst) lcArgList
                            replicates'
                                | length replicateList > 1 =
                                    errorWithoutStackTrace
                                        ("Multiple resampling replicate specifications in support command--can have only one (e.g. replicates:100): " <> show inArgs)
                                | null replicateList = Just 100
                                | otherwise = readMaybe (snd $ head replicateList) ∷ Maybe Int

                            goodBremList = filter ((`elem` ["goodmanbremer", "gb"]) . fst) lcArgList
                            {-
                            goodBremMethod
                             | length goodBremList > 1 =
                               errorWithoutStackTrace ("Multiple Goodman-Bremer method specifications in support command--can have only one (e.g. gb:tbr): " <> show inArgs)
                             | null (snd $ head goodBremList) = Just "tbr"
                             | otherwise = Just $ snd $ head goodBremList
                             -}

                            goodBremSampleList = filter ((`elem` ["gbsample"]) . fst) lcArgList
                            goodBremSample
                                | length goodBremSampleList > 1 =
                                    errorWithoutStackTrace
                                        ("Multiple Goodman-Bremer sample specifications in support command--can have only one (e.g. gbsample:1000): " <> show inArgs)
                                | null goodBremSampleList = Just (maxBound ∷ Int)
                                | otherwise = readMaybe (snd $ head goodBremSampleList) ∷ Maybe Int

                            maxParallelValue = filter ((== "maxparallel") . fst) lcArgList
                            maximizeParallel'  
                                | length maxParallelValue > 1 =
                                            errorWithoutStackTrace ("Multiple maxParallel specifications in support--can have only one: " <> show inArgs)
                                | null maxParallelValue = Just "false"
                                | null (snd $ head maxParallelValue) = errorWithoutStackTrace ("MaxParallel support option must be 'True' or 'False'" <> show inArgs)
                                | otherwise = readMaybe (show $ snd $ head maxParallelValue) ∷ Maybe String

                            maximizeParallel = if isNothing maximizeParallel' then errorWithoutStackTrace ("MaxParallel fuse option must be 'True' or 'False'" <> show inArgs)
                                               else if fromJust maximizeParallel'  == "true" then True
                                               else if fromJust maximizeParallel'  == "false" then False
                                               else errorWithoutStackTrace ("MaxParallel fuse option must be 'True' or 'False'" <> show inArgs)

                            -- set level of swap heristric intensity
                            levelList = filter ((== "level") . fst) lcArgList
                            levelNumber
                                | length levelList > 1 =
                                    errorWithoutStackTrace ("Multiple 'level' number specifications in swap command--can have only one: " <> show inArgs)
                                | null levelList = Just 2
                                | otherwise = readMaybe (snd $ head levelList) ∷ Maybe Int

                            swapLevel
                                | all ((/= "level") . fst) lcArgList = (-1)
                                | fromJust levelNumber < 0 = 0
                                | fromJust levelNumber > 3 = 3
                                | otherwise = fromJust levelNumber

                            levelParams 
                                | swapLevel == (-1) = [("joinall", ""), ("bettern",""), ("multitraverse","false")] --defaut really level 1
                                | swapLevel == (0) = [("joinall", ""), ("bestall",""), ("multitraverse","true")]
                                | swapLevel == (1) = [("joinall", ""), ("bettern",""), ("multitraverse","false")] 
                                | swapLevel == (2) = [("joinpruned", ""), ("bettern",""), ("multitraverse","false")] 
                                | swapLevel == (3) = [("joinpruned", ""), ("bestonly",""), ("multitraverse","false")] 
                                | otherwise = [("joinall", ""), ("bettern",""), ("multitraverse","false")] -- level 1

                            swapParams 
                                | useSPR = [("spr", "")]
                                | useTBR = [("tbr", "")]
                                | onlyBuild = []
                                | otherwise = [("tbr", "")] --default to TBR
                                                                                    
                        in  
                            --trace ("SG: " <> (show supportMeasure) <> " " <> (show lcArgList)) $
                            if isNothing levelNumber 
                                then errorWithoutStackTrace ("Support 'level' specification not an integer (e.g. level:2): " <> show (snd $ head jackList))
                            else if isNothing jackFreq'
                                then errorWithoutStackTrace ("Jacknife frequency not a float (e.g. jackknife:0.5) in support: " <> show (snd $ head jackList))
                                else
                                    if isNothing replicates'
                                        then
                                            errorWithoutStackTrace
                                                ("Resampling replicates specification not a string (e.g. replicates:100) in support: " <> show (snd $ head replicateList))
                                        else -- else if isNothing goodBremMethod then errorWithoutStackTrace ("Goodman-Bremer method specification not a string (e.g. goodmanBremer:SPR) in support: " <> (show (snd $ head goodBremList)) <> (show lcArgList))

                                            if isNothing goodBremSample
                                                then
                                                    errorWithoutStackTrace
                                                        ("Goodman-Bremer sample specification not an integer (e.g. gbsample:1000) in support: " <> show (snd $ head goodBremSampleList))
                                                else
                                                    let thisMethod
                                                            | (supportMeasure == Bootstrap) && ((not . null) jackList) && (null goodBremList) =
                                                                -- trace
                                                                --    "Bootstrap and Jackknife specified--defaulting to Jackknife"
                                                                Jackknife
                                                            | ((supportMeasure == Bootstrap) || ((not . null) jackList)) && ((not . null) goodBremList) =
                                                                -- trace
                                                                --    "Resampling (Bootstrap or Jackknife) and Goodman-Bremer specified--defaulting to Goodman-Bremer"
                                                                GoodmanBremer
                                                            | supportMeasure == Bootstrap = Bootstrap
                                                            | (not . null) jackList = Jackknife
                                                            | otherwise = GoodmanBremer

                                                        gbSampleSize =
                                                            if goodBremSample == Just (maxBound ∷ Int)
                                                                then Nothing
                                                                else goodBremSample

                                                        -- sample trees uniformly at random--or "nth"
                                                        gbRandomSample =
                                                            if isJust gbSampleSize
                                                                then True -- any (( == "atrandom").fst) lcArgList
                                                                else False

                                                        replicates =
                                                            if fromJust replicates' < 0
                                                                then -- logWith LogWarn "Negative replicates number--defaulting to 100"
                                                                    100
                                                                else fromJust replicates'
                                                        jackFreq =
                                                            if fromJust jackFreq' <= 0 || fromJust jackFreq' >= 1.0
                                                                then -- trace "Jackknife frequency must be on (0.0, 1.0) defaulting to 0.6321"
                                                                    0.6321
                                                                else fromJust jackFreq'

                                                        buildOptions = [("distance", ""), ("replicates", show (100 ∷ Int)), ("best", show (1 ∷ Int)), ("rdwag", "")] -- [("replicates", show 10), ("best", show 1)]
                                                        swapOptions =
                                                            if onlyBuild
                                                                then []
                                                                else swapParams <> levelParams <> [("support", ""), ("steepest", ""), ("keep", show (1 ∷ Int))]
                                                        supportGraphList =
                                                            if thisMethod == Bootstrap || thisMethod == Jackknife
                                                                then
                                                                    let extraString =
                                                                            if thisMethod == Jackknife
                                                                                then " with delete fraction  " <> show (1 - jackFreq)
                                                                                else ""
                                                                    in  do
                                                                            logWith LogInfo $
                                                                                unwords
                                                                                    [ "Generating"
                                                                                    , show thisMethod
                                                                                    , "resampling support with"
                                                                                    , show replicates
                                                                                    , "replicates"
                                                                                    , extraString
                                                                                    ]
                                                                            g ← getResampleGraph inGS inData maximizeParallel thisMethod replicates buildOptions swapOptions jackFreq
                                                                            pure [g]
                                                                else
                                                                    let neighborhood =
                                                                            if useTBR
                                                                                then "tbr"
                                                                                else
                                                                                    if useSPR
                                                                                        then "spr"
                                                                                        else "tbr"
                                                                        extraString =
                                                                            if isJust gbSampleSize
                                                                                then " using " <> neighborhood <> " based on " <> show (fromJust gbSampleSize) <> " samples at random"
                                                                                else " using " <> neighborhood
                                                                    in  do
                                                                            logWith LogInfo $ "Generating Goodman-Bremer support" <> extraString <> "\n"
                                                                            -- TODO
                                                                            mapM (getGoodBremGraphs inGS inData maximizeParallel neighborhood gbSampleSize gbRandomSample) inGraphList
                                                    in  do
                                                            -- Option warnings
                                                            when ((supportMeasure == Bootstrap) && ((not . null) jackList) && (null goodBremList)) $
                                                                logWith LogWarn "Bootstrap and Jackknife specified--defaulting to Jackknife"
                                                            when (((supportMeasure == Bootstrap) || ((not . null) jackList)) && ((not . null) goodBremList)) $
                                                                logWith LogWarn "Resampling (Bootstrap or Jackknife) and Goodman-Bremer specified--defaulting to Goodman-Bremer"
                                                            when (fromJust replicates' < 0) $
                                                                logWith LogWarn "Negative replicates number--defaulting to 100"
                                                            when (fromJust jackFreq' <= 0 || fromJust jackFreq' >= 1.0) $
                                                                logWith LogWarn "Jackknife frequency must be on (0.0, 1.0) defaulting to 0.6321"

                                                            supportGraphList


-- | getResampledGraphs performs resampling and search for Bootstrap and jackknife support
getResampleGraph
    ∷ GlobalSettings
    → ProcessedData
    -> Bool
    → SupportMethod
    → Int
    → [(String, String)]
    → [(String, String)]
    → Double
    → PhyG ReducedPhylogeneticGraph
getResampleGraph inGS inData maximizeParallel resampleType replicates buildOptions swapOptions jackFreq =
    let -- create appropriate support graph >50% ?
        -- need to add args
        reconcileArgs = case graphType inGS of
            Tree →
                [ ("method", "majority")
                , ("compare", "identity")
                --, ("edgelabel", "true")
                --, ("vertexlabel", "true")
                , ("connect", "true")
                , ("threshold", "51")
                --, ("outformat", "dot")
                ]
            _ →
                [ ("method", "eun")
                , ("compare", "identity")
                --, ("edgelabel", "true")
                --, ("vertexlabel", "true")
                , ("connect", "true")
                , ("threshold", "51")
                --, ("outformat", "dot")
                ]
        (numSets, leftOver) = divMod replicates (graphsSteepest inGS)
        -- parallel stuff
        action ∷ PhyG ReducedPhylogeneticGraph
        action = makeResampledDataAndGraph inGS inData resampleType buildOptions swapOptions jackFreq
    in  -- majority ruke consensus if no args
        do
            -- the replicate to performs number replicates
            resampledGraphList ← if maximizeParallel then
                                    getParallelChunkTraverse >>= \pTraverse →
                                        const action `pTraverse` replicate replicates ()
                                 else do
                                    firstSetList <- mapM (makeDataGraphReplicates inGS inData resampleType buildOptions swapOptions jackFreq (graphsSteepest inGS)) [0 .. numSets - 1]
                                    remainderSetList <- makeDataGraphReplicates inGS inData resampleType buildOptions swapOptions jackFreq leftOver 0
                                    pure  $ remainderSetList <> (concat firstSetList)

            recResult ← REC.makeReconcileGraph VER.reconcileArgList reconcileArgs $ fst5 <$> resampledGraphList
            let (_, reconciledGraph) = recResult

            -- generate resampled graph
            -- can't really relabel  easily without bv and maybe not necessary anyway--node numbers inconsistent
            pure (reconciledGraph, infinity, LG.empty, V.empty, V.empty)

{- | makeDataGraphReplicates is a wrapper around makeResampledDataAndGraph
    To allow for finner control of parallel execution-- ie reduce amount of prallelism due
    to memory footprint issues.
-}
makeDataGraphReplicates 
    ∷ GlobalSettings
    → ProcessedData
    → SupportMethod 
    → [(String, String)]
    → [(String, String)]
    → Double
    -> Int
    -> Int 
    -> PhyG [ReducedPhylogeneticGraph]
makeDataGraphReplicates inGS inData resampleType buildOptions swapOptions jackFreq replicates _ =
    let -- parallel stuff
        action ∷ PhyG ReducedPhylogeneticGraph
        action = makeResampledDataAndGraph inGS inData resampleType buildOptions swapOptions jackFreq
    in do
         getParallelChunkTraverse >>= \pTraverse →
                const action `pTraverse` replicate replicates ()

{- | makeResampledDataAndGraph takes paramters, resmaples data and find a graph based on search parameters
returning the resampled graph
-}
makeResampledDataAndGraph
    ∷ GlobalSettings
    → ProcessedData
    → SupportMethod
    → [(String, String)]
    → [(String, String)]
    → Double
    → PhyG ReducedPhylogeneticGraph
makeResampledDataAndGraph inGS inData resampleType buildOptions swapOptions jackFreq = do
    newData ← resampleData resampleType jackFreq inData
    -- pairwise distances for distance analysis
    -- pairwiseDistances ← DD.getPairwiseDistances newData
    -- let buildGraphs = B.buildGraph buildOptions inGS newData

    -- if not a tree then try to add net edges
    let netAddArgs = [("netadd", ""), ("keep", show (1 ∷ Int)), ("steepest", ""), ("atrandom", ""), ("maxnetedges", "5")]

    -- simple swap refinement
    if V.null $ thd3 newData
        then pure emptyReducedPhylogeneticGraph
        else do
            -- build graphs
            buildGraphs ← B.buildGraph buildOptions inGS newData
            bestBuildGraphList ← GO.selectGraphs Best (outgroupIndex inGS) (maxBound ∷ Int) 0.0 buildGraphs

            -- Do net edges if not tree
            netGraphList <- if Tree == graphType inGS then
                                    pure bestBuildGraphList
                            else do
                                    R.netEdgeMaster netAddArgs inGS newData bestBuildGraphList

            -- Do swap if specified
            swapGraphList <- if null swapOptions then
                                pure netGraphList
                             else do
                                R.swapMaster swapOptions inGS newData netGraphList

            pure $ head swapGraphList


{- | resampleData perfoms a single randomized data resampling
based on either with replacement (bootstrap) or without (jackknife)
jackknife moves through processed data and creates a new data set
   based on simple prob
Bootstrap draws chars from input directly copying--currently disabled
if a block of data end up with zero resampled characters it is deleted
-}
resampleData ∷ SupportMethod → Double → ProcessedData → PhyG ProcessedData
resampleData resampleType sampleFreq (nameVect, nameBVVect, blockDataVect)
    | V.null blockDataVect = error "Null input data in resampleData"
    | otherwise -- Bootstrap  or Jackknife resampling
        =
        let resampler = case resampleType of
                Bootstrap → resampleBlockBootstrap
                _ → resampleBlockJackknife sampleFreq
        in  do
                newBlockDataVect' ← traverse resampler blockDataVect
                -- filter any zero length blocks
                let newBlockDataVect = V.filter ((not . V.null) . thd3) newBlockDataVect'
                pure $ (nameVect, nameBVVect, newBlockDataVect)


{- |
Takes BlockData and a seed and creates a Bootstrap resampled BlockData
-}
resampleBlockBootstrap ∷ BlockData → PhyG BlockData
resampleBlockBootstrap (nameText, charDataVV, charInfoV) = do
    -- maps over taxa in data bLock
    (newCharDataVV, newCharInfoV) ← V.unzip <$> traverse (makeSampledPairVectBootstrap charInfoV) charDataVV
    pure (nameText, newCharDataVV, V.head newCharInfoV)


{- | makeSampledPairVectBootstrap takes a list of Int and a vectors of charinfo and char data
and returns new vectors of chardata and charinfo based on randomly sampled character indices
this to create a Bootstrap replicate of equal size
this for a single taxon hecen pass teh random ints so same for each one
-}
makeSampledPairVectBootstrap
    ∷ V.Vector CharInfo → V.Vector CharacterData → PhyG (V.Vector CharacterData, V.Vector CharInfo)
makeSampledPairVectBootstrap inCharInfoVect inCharDataVect =
    -- get character numbers and set resampling indices for dynamic characters--statric happen within those characters
    let -- zip so can filter static and dynamic characters
        dataInfoPairV = V.zip inCharDataVect inCharInfoVect

        -- filter static
        (staticCharsV, staticCharsInfoV) = V.unzip $ V.filter ((`elem` exactCharacterTypes) . charType . snd) dataInfoPairV

        -- filter dynamic
        (dynamicCharsV, dynamicCharsInfoV) = V.unzip $ V.filter ((`notElem` exactCharacterTypes) . charType . snd) dataInfoPairV

        numDynamicChars = V.length dynamicCharsV
    in  do
            dynCharIndices ← V.replicateM numDynamicChars $ getRandomR (0, numDynamicChars - 1)
            let resampleDynamicChars = V.map (dynamicCharsV V.!) dynCharIndices
            let resampleDynamicCharInfo = V.map (dynamicCharsInfoV V.!) dynCharIndices

            -- static chars do each one mapping random choices within the character type
            -- but keeping each one--hence char info is staticCharsInfoV
            resampleStaticChars ← V.zipWithM subSampleStatic staticCharsV staticCharsInfoV
            -- cons the vectors for chrater data and character info
            pure (resampleStaticChars <> resampleDynamicChars, staticCharsInfoV <> resampleDynamicCharInfo)


{- | subSampleStatic takes a random int list and a static charcter
Bootstrap resamples that character based on ransom it list and number of "subcharacters" in character
-}
subSampleStatic ∷ CharacterData → CharInfo → PhyG CharacterData
subSampleStatic inCharData inCharInfo =
    let (a1, a2, a3) = rangePrelim inCharData
        (na1, na2, na3) = stateBVPrelim inCharData
        (pa1, pa2, pa3) = packedNonAddPrelim inCharData
        m1 = matrixStatesPrelim inCharData
        inCharType = charType inCharInfo

        charLength
            | inCharType == Add = V.length a2
            | inCharType == NonAdd = V.length na2
            | inCharType == Matrix = V.length m1
            | inCharType `elem` packedNonAddTypes = UV.length pa2
            | otherwise = error ("Dynamic character in subSampleStatic: " <> show inCharType)
    in  do
            -- get character indices based on number "subcharacters"
            staticCharIndices ← V.generateM charLength . const $ randIndex charLength <$> getRandom
            let staticCharIndicesUV ∷ UV.Vector Int
                staticCharIndicesUV = convert staticCharIndices

            -- trace ("SSS:" <> (show $ V.length staticCharIndices) <> " " <> (show staticCharIndices)) (
            case inCharType of
                Add →
                    pure $
                        inCharData
                            { rangePrelim = (V.map (a1 V.!) staticCharIndices, V.map (a2 V.!) staticCharIndices, V.map (a3 V.!) staticCharIndices)
                            }
                NonAdd →
                    pure $
                        inCharData
                            { stateBVPrelim = (V.map (na1 V.!) staticCharIndices, V.map (na2 V.!) staticCharIndices, V.map (na3 V.!) staticCharIndices)
                            }
                val
                    | val `elem` packedNonAddTypes →
                        pure $
                            inCharData
                                { packedNonAddPrelim =
                                    (UV.map (pa1 UV.!) staticCharIndicesUV, UV.map (pa2 UV.!) staticCharIndicesUV, UV.map (pa3 UV.!) staticCharIndicesUV)
                                }
                Matrix → pure $ inCharData{matrixStatesPrelim = V.map (m1 V.!) staticCharIndices}
                _ → error ("Incorrect character type in subSampleStatic: " <> show inCharType)
    where
        randIndex ∷ ∀ {b}. (Integral b) ⇒ b → b → b
        randIndex a b = snd $ divMod (abs b) a


{- | makeSampledCharCharInfoVect takes a vector of Int and a vector of charData and a vector of charinfo
if teh data type is not static--the character is returns if Bool is True not otherwise
if the char is static (add, non add, matrix) then the bool array is applied
across the vetor of those characters (since they re vectors of charcters themselves
returned as a pair of vectors (reversed--but shouldn't matter for resampling purposes)
does not check if equal in length
-}
makeSampledVect ∷ (GV.Vector v a) ⇒ [Bool] → [a] → v a → v a
makeSampledVect boolList accumList inVect =
    if GV.null inVect
        then GV.fromList accumList
        else
            if head boolList
                then makeSampledVect (tail boolList) (GV.head inVect : accumList) (GV.tail inVect)
                else makeSampledVect (tail boolList) accumList (GV.tail inVect)


{- | makeSampledVect takes a list of Bool and avector and returns those values
with True as a vector (reversed--but shouldn't matter for resampling purposes)
does not check if equal in length
-}
makeSampledPairVect
    ∷ [Bool]
    → [Bool]
    → [CharacterData]
    → [CharInfo]
    → V.Vector CharInfo
    → V.Vector CharacterData
    → (V.Vector CharacterData, V.Vector CharInfo)
makeSampledPairVect fullBoolList boolList accumCharDataList accumCharInfoList inCharInfoVect inCharDataVect =
    if V.null inCharInfoVect
        then (V.fromList accumCharDataList, V.fromList accumCharInfoList)
        else
            let firstCharInfo = V.head inCharInfoVect
                firstCharData = V.head inCharDataVect
                firstCharType = charType firstCharInfo
            in  -- straight resample if dynamic
                if firstCharType `notElem` exactCharacterTypes
                    then
                        if head boolList
                            then
                                makeSampledPairVect
                                    fullBoolList
                                    (tail boolList)
                                    (firstCharData : accumCharDataList)
                                    (firstCharInfo : accumCharInfoList)
                                    (V.tail inCharInfoVect)
                                    (V.tail inCharDataVect)
                            else
                                makeSampledPairVect
                                    fullBoolList
                                    (tail boolList)
                                    accumCharDataList
                                    accumCharInfoList
                                    (V.tail inCharInfoVect)
                                    (V.tail inCharDataVect)
                    else -- static character--keep in sample, but need to sample in the vector

                        let (a1, a2, a3) = rangePrelim firstCharData
                            (na1, na2, na3) = stateBVPrelim firstCharData
                            (pa1, pa2, pa3) = packedNonAddPrelim firstCharData
                            m1 = matrixStatesPrelim firstCharData
                        in  if firstCharType == Add
                                then
                                    let newCharData =
                                            firstCharData
                                                { rangePrelim = (makeSampledVect fullBoolList [] a1, makeSampledVect fullBoolList [] a2, makeSampledVect fullBoolList [] a3)
                                                }
                                    in  -- trace ("Length Add: " <> (show $ V.length $ snd3 $ rangePrelim newCharData)) (
                                        if V.null (makeSampledVect fullBoolList [] a2)
                                            then
                                                makeSampledPairVect
                                                    fullBoolList
                                                    (tail boolList)
                                                    accumCharDataList
                                                    accumCharInfoList
                                                    (V.tail inCharInfoVect)
                                                    (V.tail inCharDataVect)
                                            else
                                                makeSampledPairVect
                                                    fullBoolList
                                                    (tail boolList)
                                                    (newCharData : accumCharDataList)
                                                    (firstCharInfo : accumCharInfoList)
                                                    (V.tail inCharInfoVect)
                                                    (V.tail inCharDataVect)
                                else -- )

                                    if firstCharType == NonAdd
                                        then
                                            let newCharData =
                                                    firstCharData
                                                        { stateBVPrelim = (makeSampledVect fullBoolList [] na1, makeSampledVect fullBoolList [] na2, makeSampledVect fullBoolList [] na3)
                                                        }
                                            in  -- trace ("Length NonAdd: " <> (show $ V.length $ snd3 $ stateBVPrelim newCharData))  (
                                                if V.null (makeSampledVect fullBoolList [] na2)
                                                    then
                                                        makeSampledPairVect
                                                            fullBoolList
                                                            (tail boolList)
                                                            accumCharDataList
                                                            accumCharInfoList
                                                            (V.tail inCharInfoVect)
                                                            (V.tail inCharDataVect)
                                                    else
                                                        makeSampledPairVect
                                                            fullBoolList
                                                            (tail boolList)
                                                            (newCharData : accumCharDataList)
                                                            (firstCharInfo : accumCharInfoList)
                                                            (V.tail inCharInfoVect)
                                                            (V.tail inCharDataVect)
                                        else -- )

                                            if firstCharType `elem` packedNonAddTypes
                                                then
                                                    let newCharData =
                                                            firstCharData
                                                                { packedNonAddPrelim =
                                                                    (makeSampledVect fullBoolList [] pa1, makeSampledVect fullBoolList [] pa2, makeSampledVect fullBoolList [] pa3)
                                                                }
                                                    in  -- trace ("Length NonAdd: " <> (show $ V.length $ snd3 $ stateBVPrelim newCharData))  (
                                                        if GV.null (makeSampledVect fullBoolList [] pa2)
                                                            then
                                                                makeSampledPairVect
                                                                    fullBoolList
                                                                    (tail boolList)
                                                                    accumCharDataList
                                                                    accumCharInfoList
                                                                    (GV.tail inCharInfoVect)
                                                                    (GV.tail inCharDataVect)
                                                            else
                                                                makeSampledPairVect
                                                                    fullBoolList
                                                                    (tail boolList)
                                                                    (newCharData : accumCharDataList)
                                                                    (firstCharInfo : accumCharInfoList)
                                                                    (GV.tail inCharInfoVect)
                                                                    (GV.tail inCharDataVect)
                                                else -- )

                                                    if firstCharType == Matrix
                                                        then
                                                            let newCharData = firstCharData{matrixStatesPrelim = makeSampledVect fullBoolList [] m1}
                                                            in  if V.null (makeSampledVect fullBoolList [] m1)
                                                                    then
                                                                        makeSampledPairVect
                                                                            fullBoolList
                                                                            (tail boolList)
                                                                            accumCharDataList
                                                                            accumCharInfoList
                                                                            (V.tail inCharInfoVect)
                                                                            (V.tail inCharDataVect)
                                                                    else
                                                                        makeSampledPairVect
                                                                            fullBoolList
                                                                            (tail boolList)
                                                                            (newCharData : accumCharDataList)
                                                                            (firstCharInfo : accumCharInfoList)
                                                                            (V.tail inCharInfoVect)
                                                                            (V.tail inCharDataVect)
                                                        else error ("Incorrect character type in makeSampledPairVect: " <> show firstCharType)


-- | resampleBlockJackknife takes BlockData and a seed and creates a jackknife resampled BlockData
resampleBlockJackknife ∷ Double → BlockData → PhyG BlockData
resampleBlockJackknife sampleFreq inData@(nameText, charDataVV, charInfoV) =
    let getRandomAcceptances ∷ PhyG [Bool]
        getRandomAcceptances = fmap (randAccept sampleFreq) <$> getRandoms

        randAccept ∷ Double → Word → Bool
        randAccept b a =
            let (_, randVal) = divMod (abs a) 1000
                critVal = floor (1000 * b)
            in  -- trace ("RA : " <> (show (b,a, randVal, critVal, randVal < critVal)))
                randVal < critVal

        jackknifeSampling ∷ PhyG (V.Vector (V.Vector CharacterData), V.Vector (V.Vector CharInfo))
        jackknifeSampling = do
            accepts1 ← getRandomAcceptances
            accepts2 ← getRandomAcceptances
            pure . V.unzip $ makeSampledPairVect accepts1 accepts2 [] [] charInfoV <$> charDataVV
    in  do
            (newCharDataVV, newCharInfoVV) ← jackknifeSampling
            let newCharInfoV ∷ V.Vector CharInfo
                newCharInfoV = V.head newCharInfoVV

            --logWith LogInfo $ "Jacknife sample size: " <> (show $ V.length charInfoV) <> " -> " <> (show $ V.length newCharInfoV)

            case V.length newCharInfoV of
                0 → resampleBlockJackknife sampleFreq inData
                _ → pure (nameText, newCharDataVV, newCharInfoV)


{- | getGoodBremGraphs performs Goodman-Bremer support
examines complete SPR or TBR swap neighborhood chekcing for presence/absence of edges in input Graph List
can do sample of trees either "nth" or at random if specified
sample based on SPR-- 4n^2 - 26n - 42 for TBR 8n^3 for now
this will only examine bridge edges for networks, networkedge values willl be doen via net delete
MAPs for each graph?
-}
getGoodBremGraphs
    ∷ GlobalSettings 
    → ProcessedData 
    → Bool 
    -> String 
    → Maybe Int 
    → Bool 
    → ReducedPhylogeneticGraph 
    → PhyG ReducedPhylogeneticGraph
getGoodBremGraphs inGS inData maximizeParallel swapType sampleSize sampleAtRandom inGraph =
    if LG.isEmpty (fst5 inGraph)
        then error "Null graph in getGoodBremGraphs" -- maybe should be error?
        else do
            -- create list of edges for input graph and a structure with egde node indices and bitvector values
            -- requires index BV of each node

            let tupleList = getGraphTupleList inGraph

            -- traverse neighborhood (and net edge removal) keeping min cost without edges
            supportEdgeTupleList ← getGBTuples inGS inData maximizeParallel swapType sampleSize sampleAtRandom tupleList inGraph

            let simpleGBGraph = LG.mkGraph (LG.labNodes $ fst5 inGraph) (fmap (tupleToSimpleEdge (snd5 inGraph)) supportEdgeTupleList)
            -- trace ("GGBG: " <> (show $ length tupleList) <> " -> " <> (show $ length supportEdgeTupleList))
            pure (simpleGBGraph, snd5 inGraph, thd5 inGraph, fth5 inGraph, fft5 inGraph)
    where
        tupleToSimpleEdge
            ∷ ∀ {c1} {a} {b} {c2} {d}
             . (Num c1)
            ⇒ c1
            → (a, b, c2, d, c1)
            → (a, b, c1)
        tupleToSimpleEdge d (a, b, _, _, c) = (a, b, c - d)


-- | getGraphTupleList takes a graph and cost (maybe initialized to infinity) returns tuple list
getGraphTupleList ∷ ReducedPhylogeneticGraph → [(Int, Int, NameBV, NameBV, VertexCost)]
getGraphTupleList inGraph =
    if LG.isEmpty (fst5 inGraph)
        then error "Null graph in getGraphTupleList"
        else
            let egdeList = LG.edges (fst5 inGraph)

                -- graph node list
                nodeList = LG.labNodes (thd5 inGraph)
                nodeIndexBVPairList = fmap makeindexBVPair nodeList

                -- list of vectors for contant time access via index = fst (a, bv)
                nodeIndexBVPairVect = V.fromList nodeIndexBVPairList

                -- make tuple for each edge in each graph
                -- (uIndex,vINdex,uBV, vBV, graph cost)
                tupleList = makeGraphEdgeTuples nodeIndexBVPairVect infinity egdeList
            in  tupleList
    where
        makeindexBVPair ∷ ∀ {a}. (a, VertexInfo) → (a, NameBV)
        makeindexBVPair (a, b) = (a, bvLabel b)


{- | getGBTuples takes a tuple list from graph containing initialized values and update those values based
on each graph in the inGraph neigborhood
first does this via swap--for network does edge net edge in turn by removing using netDel
-}
getGBTuples
    ∷ GlobalSettings
    → ProcessedData
    -> Bool 
    → String
    → Maybe Int
    → Bool
    → [(Int, Int, NameBV, NameBV, VertexCost)]
    → ReducedPhylogeneticGraph
    → PhyG [(Int, Int, NameBV, NameBV, VertexCost)]
getGBTuples inGS inData maximizeParallel swapType sampleSize sampleAtRandom inTupleList inGraph = do
    -- traverse swap (SPR/TBR) neighborhood optimizing each graph fully
    swapTuples ← performGBSwap inGS inData maximizeParallel swapType sampleSize sampleAtRandom inTupleList inGraph
    case graphType inGS of
        -- swap only for Tree-do nothing
        Tree → pure swapTuples
        _ | LG.isTree (fst5 inGraph) → pure swapTuples
        -- network edge support if not Tree
        -- SoftWired => delete edge -- could add net move if needed
        SoftWired →
            let deleteAction ∷ (Int, Int, NameBV, NameBV, VertexCost) → PhyG (Int, Int, NameBV, NameBV, VertexCost)
                deleteAction = updateDeleteTuple inGS inData inGraph
            in  getParallelChunkTraverse >>= \pTraverse →
                    deleteAction `pTraverse` swapTuples
        -- HardWired => move edge
        _ →
            let moveAction ∷ (Int, Int, NameBV, NameBV, VertexCost) → PhyG (Int, Int, NameBV, NameBV, VertexCost)
                moveAction = updateMoveTuple inGS inData inGraph
            in  getParallelChunkTraverse >>= \pTraverse →
                    moveAction `pTraverse` swapTuples


{- | updateDeleteTuple take a graph and and edge and delete a network edge (or retunrs tuple if not network)
if this were a HardWWired graph--cost would always go down, so only applied to softwired graphs
-}
updateDeleteTuple
    ∷ GlobalSettings
    → ProcessedData
    → ReducedPhylogeneticGraph
    → (Int, Int, NameBV, NameBV, VertexCost)
    → PhyG (Int, Int, NameBV, NameBV, VertexCost)
updateDeleteTuple inGS inData inGraph inTuple@(inE, inV, inEBV, inVBV, inCost) =
    let isNetworkEdge = LG.isNetworkEdge (fst5 inGraph) (inE, inV)
    in  if not isNetworkEdge
            then do
                pure inTuple
            else do
                -- True to force full evalutation
                deleteRedgeResults ← N.deleteNetEdge inGS inData inGraph True (inE, inV)
                let deleteCost = snd5 deleteRedgeResults
                pure (inE, inV, inEBV, inVBV, min inCost deleteCost)


{- | updateMoveTuple take a graph and and edge and moves a network edge (or returns tuple if not network)
if this were a HardWWired graph--cost would always go down, so only applied to softwired graphs
max bound because its a place holder for max num net edges
-}
updateMoveTuple
    ∷ GlobalSettings
    → ProcessedData
    → ReducedPhylogeneticGraph
    → (Int, Int, NameBV, NameBV, VertexCost)
    → PhyG (Int, Int, NameBV, NameBV, VertexCost)
updateMoveTuple inGS inData inGraph inTuple@(inE, inV, inEBV, inVBV, inCost) =
    let isNetworkEdge = LG.isNetworkEdge (fst5 inGraph) (inE, inV)
    in  if not isNetworkEdge
            then do
                pure inTuple
            else -- True to force full evalutation

                let steepest = False
                    randomOrder = False
                    keepNum = 10 -- really could be one since sorted by cost, but just to make sure)Order
                    saParams ∷ ∀ {a}. Maybe a
                    saParams = Nothing
                in  do
                        deleteAddGraphs ←
                            N.deleteOneNetAddAll inGS inData (maxBound ∷ Int) keepNum steepest randomOrder inGraph [(inE, inV)] saParams
                        let moveCost = minimum (snd5 <$> deleteAddGraphs)

                        pure (inE, inV, inEBV, inVBV, min inCost moveCost)


{- | performGBSwap takes parameters and  graphs and traverses swap neighborhood
examining each (or nth, or random) Graphs examining each ecah in each graph for Goodman-Bremer
optimality support
-}
performGBSwap
    ∷ GlobalSettings
    → ProcessedData
    -> Bool
    → String
    → Maybe Int
    → Bool
    → [(Int, Int, NameBV, NameBV, VertexCost)]
    → ReducedPhylogeneticGraph
    → PhyG [(Int, Int, NameBV, NameBV, VertexCost)]
performGBSwap inGS inData maximizeParallel swapType sampleSize sampleAtRandom inTupleList inGraph
    | LG.isEmpty (fst5 inGraph) = error "Null graph in performGBSwap"
    | otherwise =
        let -- work with simple graph
            inSimple = fst5 inGraph
            (firstRootIndex, _) = head $ LG.getRoots inSimple

            -- determine edges to break on--'bridge' edges only for network
            -- filter out edges from root since no use--would just rejoin
            breakEdgeList = case graphType inGS of
                Tree → filter ((/= firstRootIndex) . fst3) $ LG.labEdges inSimple
                _ → filter ((/= firstRootIndex) . fst3) $ LG.getEdgeSplitList inSimple
        in  do
                -- integerized critical value for prob accept
                -- based on approx (leaves - netnodes)^2 or (leaves - netnodes)^3
                let (_, leafList, _, netVertList) = LG.splitVertexList (fst5 inGraph)
                let intProbAccept = case swapType of
                        "spr" →
                            floor
                                ((1000.0 * fromIntegral (fromJust sampleSize)) / ((2.0 * fromIntegral (length leafList - length netVertList)) ** 2) ∷ Double)
                        _ →
                            floor
                                ((1000.0 * fromIntegral (fromJust sampleSize)) / ((2.0 * fromIntegral (length leafList - length netVertList)) ** 3) ∷ Double)

                -- logWith LogWarn $ "PGBS: " <> (show (intProbAccept, (fromJust sampleSize), ((2.0 * fromIntegral (length leafList - length netVertList)) ** 2) ∷ Double, ((2.0 * fromIntegral (length leafList - length netVertList)) ** 3) ∷ Double)) <> "\n"
                -- splitRejoinAction ∷ ([Int], LG.LEdge Double) → PhyG [(Int, Int, NameBV, NameBV, VertexCost)]
                let splitRejoinAction = splitRejoinGB inGS inData swapType intProbAccept sampleAtRandom inTupleList inSimple breakEdgeList

                -- generate tuple lists for each break edge parallelized at this level
                tupleListList ← if maximizeParallel then
                                    getParallelChunkTraverse >>= \pTraverse →
                                        splitRejoinAction `pTraverse` breakEdgeList
                                else mapM splitRejoinAction breakEdgeList

                -- merge tuple lists--should all be in same order
                let newTupleList = mergeTupleLists (filter (not . null) tupleListList) []
                -- trace ("PGBS:" <> (show $ fmap length tupleListList) <> " -> " <> (show $ length newTupleList))
                pure newTupleList


{- | splitRejoinGB take parameters and splits input graph at specified edge and rejoins at all available edge
(reroots the pruned subgraph if TBR) and creates and gets cost of graph (lazy takes care of post order only)
with optimized graph, tuple list is creted and compared to input graph tuple list.
original edge was (parentPrunedGraphRoot, prunedGraphRootIndex)
working with SimpleGraph
-}
splitRejoinGB
    ∷ GlobalSettings
    → ProcessedData
    → String
    → Int
    → Bool
    → [(Int, Int, NameBV, NameBV, VertexCost)]
    → SimpleGraph
    → [LG.LEdge Double]
    → LG.LEdge Double
    → PhyG [(Int, Int, NameBV, NameBV, VertexCost)]
splitRejoinGB inGS inData swapType intProbAccept sampleAtRandom inTupleList inGraph originalBreakEdgeList breakEdge =
    let -- split graph on breakEdge
        (splitGraph, _, prunedGraphRootIndex, _, _, edgeDeleteList) = LG.splitGraphOnEdge' inGraph breakEdge

        -- get edges in base graph to be invaded (ie not in pruned graph)
        prunedGraphRootNode = (prunedGraphRootIndex, fromJust $ LG.lab splitGraph prunedGraphRootIndex)
        (prunedSubTreeNodes, prunedSubTreeEdges) = LG.nodesAndEdgesAfter splitGraph [prunedGraphRootNode]
        edgesNotToInvade = (LG.toEdge breakEdge : edgeDeleteList) <> fmap LG.toEdge prunedSubTreeEdges
        edgesToInvade = filter (LG.notMatchEdgeIndices edgesNotToInvade) originalBreakEdgeList

        -- rejoin, evaluate, get better tuple
        -- check if there are tbr-type rearrangements to do (rerooting pruned graph)
        -- create TBR rerooot split graphs if required
        splitGraphList =
            if (length prunedSubTreeNodes < 3) || swapType == "spr"
                then [splitGraph]
                else -- generate "tbr" rerootings in split graph
                    getTBRSplitGraphs inGS splitGraph breakEdge

        action ∷ LG.LEdge Double → PhyG [(Int, Int, NameBV, NameBV, VertexCost)]
        action = rejoinGB inGS inData intProbAccept sampleAtRandom inTupleList splitGraphList breakEdge
    in  do
            -- parallel at break level above
            rejoinTupleListList ←
                getParallelChunkTraverse >>= \pTraverse →
                    action `pTraverse` edgesToInvade

            -- merge tuples
            pure $ mergeTupleLists rejoinTupleListList []


{- | rejoinGB rejoins split graph at specific edge, id SPR then that's it, if TBR reroot pruned subgraph
splitGraph is SimpleGraph
the rejoin is SPR type relying on teh list lengt of split graph to present the TBR reroots
-}
rejoinGB
    ∷ GlobalSettings
    → ProcessedData
    → Int
    → Bool
    → [(Int, Int, NameBV, NameBV, VertexCost)]
    → [SimpleGraph]
    → LG.LEdge Double
    → LG.LEdge Double
    → PhyG [(Int, Int, NameBV, NameBV, VertexCost)]
rejoinGB inGS inData intProbAccept sampleAtRandom inTupleList splitGraphList originalBreakEdge@(eBreak, _, _) edgeToInvade = case splitGraphList of
    [] → pure inTupleList
    splitGraph : otherGraphs →
        let proceedWithSampling
                | not sampleAtRandom = pure False
                | otherwise = do
                                -- getRandomR (0, 999) >>= \rVal → pure $ rVal <= intProbAccept
                                rVal <- getRandomR (0, 999)
                                --- logWith LogWarn $ "RV: " <> (show rVal)
                                pure $ rVal <= intProbAccept

            
            rejoinUsingTuples givenTuples =
                rejoinGB
                    inGS
                    inData
                    intProbAccept
                    sampleAtRandom
                    givenTuples
                    otherGraphs
                    originalBreakEdge
                    edgeToInvade
            

            resultOfSampling = rejoinUsingTuples inTupleList

        in do
            shouldSampleRandomly ← proceedWithSampling

            -- Shortcircuit if not to sample based on randval and critical value
            if sampleAtRandom && (not shouldSampleRandomly) then do
                --logWith LogInfo $ "-" <> (show (sampleAtRandom, not shouldSampleRandomly))
                resultOfSampling

            else
                let resultWithoutSampling =
                        let newGraph = LG.joinGraphOnEdge splitGraph edgeToInvade eBreak
                            --pruneEdges = False
                            --warnPruneEdges = False

                            startVertex ∷ ∀ {a}. Maybe a
                            startVertex = Nothing

                            nonExactCharacters = U.getNumberSequenceCharacters (thd3 inData)

                            leafGraph = if graphType inGS /= SoftWired then
                                            GO.makeLeafGraph inData
                                        else POSW.makeLeafGraphSoftWired inGS inData

                            {- Moved to monadic part so could use postoder function

                            generatedResult = T.multiTraverseFullyLabelGraphReduced inGS inData pruneEdges warnPruneEdges startVertex newGraph

                            generaterNewGraph
                                | graphType inGS == Tree || LG.isTree newGraph || ((not . LG.cyclic) newGraph && (not . LG.parentInChain) newGraph) =
                                    generatedResult
                                | otherwise = pure emptyReducedPhylogeneticGraph
                            -}
                        in  do
                                generatedResult <- T.generalizedGraphPostOrderTraversal 
                                                    inGS
                                                    nonExactCharacters
                                                    inData
                                                    Nothing
                                                    leafGraph
                                                    False
                                                    startVertex
                                                    newGraph

                                let newPhylogeneticGraph 
                                        | graphType inGS == Tree || LG.isTree newGraph || ((not . LG.cyclic) newGraph && (not . LG.parentInChain) newGraph) =
                                            GO.convertPhylogeneticGraph2Reduced $ fst generatedResult
                                        | otherwise = emptyReducedPhylogeneticGraph

                                -- newPhylogeneticGraph ← pure generaterNewGraph
                                
                                let tupleList
                                        | newPhylogeneticGraph == emptyReducedPhylogeneticGraph = inTupleList
                                        -- update tuple list based on new graph
                                        | otherwise = getLowerGBEdgeCost inTupleList newPhylogeneticGraph -- ((2 * numTaxa) -1)
                                rejoinUsingTuples tupleList
                in  do
                        -- there were issues with this logic and rVal above
                        --shouldSampleRandomly ← proceedWithSampling
                        --if shouldSampleRandomly
                        --    then resultOfSampling
                        --    else resultWithoutSampling
                        resultWithoutSampling

-- | mergeTupleLists takes a list of list of tuples and merges them choosing the better each recursive round
mergeTupleLists
    ∷ [[(Int, Int, NameBV, NameBV, VertexCost)]]
    → [(Int, Int, NameBV, NameBV, VertexCost)]
    → [(Int, Int, NameBV, NameBV, VertexCost)]
mergeTupleLists inTupleListList accumList
    | null inTupleListList = accumList
    | null accumList = mergeTupleLists (tail inTupleListList) (head inTupleListList)
    | otherwise =
        let firstTupleList = head inTupleListList
            newTupleList = zipWith chooseBetterTuple firstTupleList accumList
        in  mergeTupleLists (tail inTupleListList) newTupleList


-- | chooseBetterTuple takes two (Int, Int, NameBV, NameBV, VertexCost) and returns better cost
chooseBetterTuple
    ∷ (Int, Int, NameBV, NameBV, VertexCost) → (Int, Int, NameBV, NameBV, VertexCost) → (Int, Int, NameBV, NameBV, VertexCost)
chooseBetterTuple (aE, aV, aEBV, aVBV, aCost) (_, _, _, _, bCost) = (aE, aV, aEBV, aVBV, min aCost bCost)


{- | makeGraphEdgeTuples take node and edge,cost tuples from a graph and returns a list of tuples of the form
(uIndex,vINdex,uBV, vBV, graph cost)
this for edge comparisons for Goodman-Bremer and other optimality-type support
-}
makeGraphEdgeTuples ∷ V.Vector (Int, NameBV) → VertexCost → [(Int, Int)] → [(Int, Int, NameBV, NameBV, VertexCost)]
makeGraphEdgeTuples nodeBVVect graphCost edgeList =
    -- trace ("MET: " <> (show $ V.length nodeBVVect))
    fmap (make5Tuple nodeBVVect graphCost) edgeList
    where
        make5Tuple
            ∷ ∀ {a} {d} {e}
             . V.Vector (a, d)
            → e
            → (Int, Int)
            → (Int, Int, d, d, e)
        make5Tuple nv c (a, b) = (a, b, snd (nv V.! a), snd (nv V.! b), c)


{- | getLowerGBEdgeCost take a list of edge tuples of (uIndex,vINdex,uBV, vBV, graph cost) from the graph
whose supports are being calculated and a new graph and updates the edge cost (GB value) if that edge
is NOT present in the graph taking the minimum of the original GB value and the new graph cost
-}
getLowerGBEdgeCost
    ∷ [(Int, Int, NameBV, NameBV, VertexCost)] → ReducedPhylogeneticGraph → [(Int, Int, NameBV, NameBV, VertexCost)]
getLowerGBEdgeCost edgeTupleList inGraph =
    if LG.isEmpty (fst5 inGraph) || null edgeTupleList
        then error "Empty graph or null edge tuple list in getLowerGBEdgeCost"
        else
            let inGraphTupleList = getGraphTupleList inGraph
            in  fmap (updateEdgeTuple (snd5 inGraph) inGraphTupleList) edgeTupleList


{- | updateEdgeTuple checks is edge is NOT in input graph edge tuple list and if not takes minimum
of edge cost GB value and in graph cost, else returns unchanged
-}
updateEdgeTuple
    ∷ VertexCost
    → [(Int, Int, NameBV, NameBV, VertexCost)]
    → (Int, Int, NameBV, NameBV, VertexCost)
    → (Int, Int, NameBV, NameBV, VertexCost)
updateEdgeTuple inGraphCost inGraphTupleList (uIndex, vIndex, uBV, vBV, edgeGBValue) =
    let edgeNotFoundCost = getNotFoundCost uBV vBV inGraphCost inGraphTupleList
    in  if isNothing edgeNotFoundCost
            then (uIndex, vIndex, uBV, vBV, edgeGBValue)
            else (uIndex, vIndex, uBV, vBV, min edgeGBValue (fromJust edgeNotFoundCost))


{- | getNotFoundCost take a pair of BitVectors (of vertices in graph) from an edge
and a list of  (Int, Int, NameBV, NameBV, VertexCost) tuples and returns
Nothing is the BVs of the two match (= signifying edge in graph) or
Just graph cost if not present for Goodman-Bremer calculations
-}
getNotFoundCost ∷ NameBV → NameBV → VertexCost → [(Int, Int, NameBV, NameBV, VertexCost)] → Maybe VertexCost
getNotFoundCost uBV vBV inTupleCost inTupleList =
    if null inTupleList
        then Just inTupleCost
        else
            let (_, _, uInBV, vInBV, _) = head inTupleList
            in  if uBV == uInBV && vBV == vInBV
                    then Nothing
                    else getNotFoundCost uBV vBV inTupleCost (tail inTupleList)


{- | getTBRSplitGraphs takes a split gaph and the original split edge and
returns a list of rerooted subgrahs split graphs suitable for rejoining
via SPR-type rejoin each to generate TBR neighborhood
much of this is modified from Swap.hs but removing data and delta portions
-}
getTBRSplitGraphs ∷ GlobalSettings → SimpleGraph → LG.LEdge Double → [SimpleGraph]
getTBRSplitGraphs inGS splitGraph splitEdge =
    if LG.isEmpty splitGraph
        then error "Empty graph in getTBRSplitGraphs"
        else -- get edges in pruned graph and reroot on those edges that are 1) not from original "root" of prune
        -- and 2) not network edges

            let prunedGraphRootNode = (snd3 splitEdge, fromJust $ LG.lab splitGraph $ snd3 splitEdge)
                edgesInPrunedSubGraph = snd $ LG.nodesAndEdgesAfter splitGraph [prunedGraphRootNode]

                nonNetWorkEdgeList =
                    if graphType inGS /= Tree
                        then filter (not . LG.isNetworkLabEdge splitGraph) edgesInPrunedSubGraph
                        else edgesInPrunedSubGraph

                -- original pruned root edges
                prunedRootEdges = LG.out splitGraph $ fst prunedGraphRootNode

                -- edges available for rerooting
                edgeAfterList = nonNetWorkEdgeList L.\\ prunedRootEdges

                -- get edges to add and delete for TBR rerooting
                tbrEdits = fmap (getTBREdits splitGraph prunedGraphRootNode edgesInPrunedSubGraph . LG.toEdge) edgeAfterList

                -- TBR split graph list
                tbrGraphList = fmap (LG.insertDeleteEdges splitGraph) tbrEdits
            in  splitGraph : tbrGraphList


{- | getTBREdits takes and edge and returns the list of edits to pruned subgraph
as a pair of edges to add and those to delete
since reroot edge is directed (e,v), edges away from v will have correct
orientation. Edges between 'e' and the root will have to be flipped
original root edges and reroort edge are deleted and new root and edge spanning orginal root created
returns ([add], [delete])
modified from function in swap to be more general and operate on SimpleGraphs as are used here
-}
getTBREdits ∷ (Eq a, Eq b) ⇒ LG.Gr a b → LG.LNode a → [LG.LEdge b] → LG.Edge → ([LG.LEdge b], [LG.Edge])
getTBREdits inGraph prunedGraphRootNode edgesInPrunedSubGraph rerootEdge =
    -- trace ("Gettiung TBR Edits for " <> (show rerootEdge)) (
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
        dummyEdgeLabel = thd3 $ head edgesInPrunedSubGraph

        -- new edges on new root position and spanning old root
        -- add in closer vertex to root to make sure direction of edge is correct
        newEdgeOnOldRoot =
            if snd3 (head originalRootEdges) `elem` (fst rerootEdge : fmap fst nodesInPath)
                then (snd3 $ head originalRootEdges, snd3 $ last originalRootEdges, dummyEdgeLabel)
                else (snd3 $ last originalRootEdges, snd3 $ head originalRootEdges, dummyEdgeLabel)
        newRootEdges = [(prunedGraphRootIndex, fst rerootEdge, dummyEdgeLabel), (prunedGraphRootIndex, snd rerootEdge, dummyEdgeLabel)]
    in  -- original root edge so no change
        if fst rerootEdge `elem` originalRootEdgeNodes && snd rerootEdge `elem` originalRootEdgeNodes
            then ([], [])
            else -- rerooted

            -- delete orignal root edges and rerootEdge
            -- add new root edges
            -- and new edge on old root--but need orientation
            -- flip edges from new root to old (delete and add list)
            -- trace ("\n\nIn Graph:\n" <> (LG.prettify $ GO.convertDecoratedToSimpleGraph inGraph) <> "\nTBR Edits: " <> (show (rerootEdge, prunedGraphRootIndex, fmap LG.toEdge flippedEdges))
            --   <> "\nEdges to add: " <> (show $ fmap LG.toEdge $ newEdgeOnOldRoot : (flippedEdges <> newRootEdges)) <> "\nEdges to delete: " <> (show $ rerootEdge : (fmap LG.toEdge (edgesToFlip <> originalRootEdges))))
                (newEdgeOnOldRoot : (flippedEdges <> newRootEdges), rerootEdge : fmap LG.toEdge (edgesToFlip <> originalRootEdges))

-- )
