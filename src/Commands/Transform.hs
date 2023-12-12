{-# OPTIONS_GHC -Wno-missed-specialisations #-}

{- |
Module to coordinate transform command execution.
-}
module Commands.Transform (
    transform,
    makeStaticApprox,
) where

import Bio.DynamicCharacter.Element (SlimState, WideState)
import Commands.CommandUtilities qualified as CU
import Commands.Verify qualified as VER
import Control.Monad (when)
import Control.Monad.IO.Class (MonadIO (..))
import Control.Parallel.Strategies
import Data.Alphabet
import Data.BitVector.LittleEndian qualified as BV
import Data.Bits
import Data.Char
import Data.Char qualified as C
import Data.List qualified as L
import Data.Maybe
import Data.Text.Lazy qualified as TL
import Data.Vector qualified as V
import Data.Vector.Generic qualified as GV
import Data.Vector.Storable qualified as SV
import Data.Vector.Unboxed qualified as UV
import GeneralUtilities
import GraphOptimization.Traversals qualified as T
import Graphs.GraphOperations qualified as GO
import Input.BitPack qualified as BP
import Input.Reorganize qualified as R
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import PHANE.Evaluation.Verbosity (Verbosity (..))
import Text.Read
import Types.Types
import Utilities.LocalGraph qualified as LG
import Utilities.Utilities qualified as U


{- | transform changes aspects of data sande settings during execution
as opposed to Set with all happens at begginign of program execution
-}
transform
    ∷ [Argument]
    → GlobalSettings
    → ProcessedData
    → ProcessedData
    → Int
    → [ReducedPhylogeneticGraph]
    → PhyG (GlobalSettings, ProcessedData, ProcessedData, [ReducedPhylogeneticGraph])
transform inArgs inGS origData inData rSeed inGraphList =
    let fstArgList = fmap (fmap toLower . fst) inArgs
        sndArgList = fmap (fmap toLower . snd) inArgs
        lcArgList = zip fstArgList sndArgList
        checkCommandList = checkCommandArgs "transform" fstArgList VER.transformArgList
    in  -- check for valid command options
        if not checkCommandList
            then errorWithoutStackTrace ("Unrecognized command in 'transform': " <> show inArgs)
            else
                let displayBlock = filter ((== "displaytrees") . fst) lcArgList
                    numDisplayTrees
                        | length displayBlock > 1 =
                            errorWithoutStackTrace ("Multiple displayTree number specifications in transform--can have only one: " <> show inArgs)
                        | null displayBlock = Just 10
                        | null (snd $ head displayBlock) = Just 10
                        | otherwise = readMaybe (snd $ head displayBlock) ∷ Maybe Int

                    toTree = any ((== "totree") . fst) lcArgList
                    toSoftWired = any ((== "tosoftwired") . fst) lcArgList
                    toHardWired = any ((== "tohardwired") . fst) lcArgList
                    toStaticApprox = any ((== "staticapprox") . fst) lcArgList
                    toDynamic = any ((== "dynamic") . fst) lcArgList
                    atRandom = any ((== "atrandom") . fst) lcArgList
                    chooseFirst = any ((== "first") . fst) lcArgList
                    reWeight = any ((== "weight") . fst) lcArgList
                    changeEpsilon = any ((== "dynamicepsilon") . fst) lcArgList
                    reRoot = any ((== "outgroup") . fst) lcArgList
                    changeGraphsSteepest = any ((== "graphssteepest") . fst) lcArgList
                    changeSoftwiredMethod = any ((== "softwiredmethod") . fst) lcArgList
                    changeGraphFactor = any ((== "graphfactor") . fst) lcArgList
                    changeCompressionResolutions = any ((== "compressresolutions") . fst) lcArgList
                    changeMultiTraverse = any ((== "multitraverse") . fst) lcArgList
                    changeUnionThreshold = any ((== "jointhreshold") . fst) lcArgList

                    reweightBlock = filter ((== "weight") . fst) lcArgList
                    weightValue
                        | length reweightBlock > 1 =
                            errorWithoutStackTrace ("Multiple weight specifications in transform--can have only one: " <> show inArgs)
                        | null reweightBlock = Just 1.0
                        | null (snd $ head reweightBlock) = Just 1
                        | otherwise = readMaybe (snd $ head reweightBlock) ∷ Maybe Double

                    changeCompressionBlock = filter ((== "compressresolutions") . fst) lcArgList
                    compressionValue
                        | length changeCompressionBlock > 1 =
                            errorWithoutStackTrace ("Multiple compressResolutions specifications in transform--can have only one: " <> show inArgs)
                        | null changeCompressionBlock = Just $ fmap toLower $ show $ compressResolutions inGS
                        | null (snd $ head changeCompressionBlock) = Just $ fmap toLower $ show $ compressResolutions inGS
                        | otherwise = readMaybe (show $ snd $ head changeCompressionBlock) ∷ Maybe String

                    changeEpsilonBlock = filter ((== "dynamicepsilon") . fst) lcArgList
                    epsilonValue
                        | length changeEpsilonBlock > 1 =
                            errorWithoutStackTrace ("Multiple dynamicEpsilon specifications in transform--can have only one: " <> show inArgs)
                        | null changeEpsilonBlock = Just $ dynamicEpsilon inGS
                        | null (snd $ head changeEpsilonBlock) = Just $ dynamicEpsilon inGS
                        | otherwise = readMaybe (snd $ head changeEpsilonBlock) ∷ Maybe Double

                    changeGraphFactorBlock = filter ((== "graphfactor") . fst) lcArgList
                    newGraphFactor
                        | length changeGraphFactorBlock > 1 =
                            errorWithoutStackTrace ("Multiple graphFactor specifications in transform--can have only one: " <> show inArgs)
                        | null changeGraphFactorBlock = Just $ fmap toLower $ show $ graphFactor inGS
                        | null (snd $ head changeGraphFactorBlock) = Just $ fmap toLower $ show $ graphFactor inGS
                        | otherwise = readMaybe (show $ snd $ head changeGraphFactorBlock) ∷ Maybe String

                    changeGraphsSteepestBlock = filter ((== "graphssteepest") . fst) lcArgList
                    newGraphsSteepest
                        | length changeGraphsSteepestBlock > 1 =
                            errorWithoutStackTrace ("Multiple graphsSteepest specifications in transform--can have only one: " <> show inArgs)
                        | null changeGraphsSteepestBlock = Just $ graphsSteepest inGS
                        | null (snd $ head changeGraphsSteepestBlock) = Just $ graphsSteepest inGS
                        | otherwise = readMaybe (snd $ head changeGraphsSteepestBlock) ∷ Maybe Int

                    changeMultiTraverseBlock = filter ((== "multitraverse") . fst) lcArgList
                    multiTraverseValue
                        | length changeMultiTraverseBlock > 1 =
                            errorWithoutStackTrace ("Multiple multiTraverse specifications in transform--can have only one: " <> show inArgs)
                        | null changeMultiTraverseBlock = Just $ fmap toLower $ show $ multiTraverseCharacters inGS
                        | null (snd $ head changeMultiTraverseBlock) = Just $ fmap toLower $ show $ multiTraverseCharacters inGS
                        | otherwise = readMaybe (show $ snd $ head changeMultiTraverseBlock) ∷ Maybe String

                    changeSoftwiredMethodBlock = filter ((== "softwiredmethod") . fst) lcArgList
                    newSoftwiredMethod
                        | length changeSoftwiredMethodBlock > 1 =
                            errorWithoutStackTrace ("Multiple softwiredMethod specifications in transform--can have only one: " <> show inArgs)
                        | null changeSoftwiredMethodBlock = Just $ fmap toLower $ show $ softWiredMethod inGS
                        | null (snd $ head changeSoftwiredMethodBlock) = Just $ fmap toLower $ show $ softWiredMethod inGS
                        | otherwise = readMaybe (show $ snd $ head changeSoftwiredMethodBlock) ∷ Maybe String

                    reRootBlock = filter ((== "outgroup") . fst) lcArgList
                    outgroupValue
                        | length reRootBlock > 1 =
                            errorWithoutStackTrace ("Multiple outgroup specifications in transform--can have only one: " <> show inArgs)
                        | null reRootBlock = Just $ outGroupName inGS
                        | null (snd $ head reRootBlock) = Just $ outGroupName inGS
                        | otherwise = readMaybe (snd $ head reRootBlock) ∷ Maybe TL.Text

                    changeUnionBlock = filter ((== "jointhreshold") . fst) lcArgList
                    unionValue
                        | length changeUnionBlock > 1 =
                            errorWithoutStackTrace ("Multiple joinThreshold specifications in transform--can have only one: " <> show inArgs)
                        | null changeUnionBlock = Just $ unionThreshold inGS
                        | null (snd $ head changeUnionBlock) = Just $ unionThreshold inGS
                        | otherwise = readMaybe (snd $ head changeUnionBlock) ∷ Maybe Double

                    nameList = fmap TL.pack $ fmap (filter (/= '"')) $ fmap snd $ filter ((== "name") . fst) lcArgList
                    charTypeList = fmap snd $ filter ((== "type") . fst) lcArgList
                in  if (length $ filter (== True) [toTree, toSoftWired, toHardWired]) > 1
                        then errorWithoutStackTrace ("Multiple graph transform commands--can only have one : " <> (show inArgs))
                        else
                            if toStaticApprox && toDynamic
                                then errorWithoutStackTrace ("Multiple staticApprox/Dynamic transform commands--can only have one : " <> (show inArgs))
                                else
                                    if atRandom && chooseFirst
                                        then
                                            errorWithoutStackTrace
                                                ("Multiple display tree choice commands in transform (first, atRandom)--can only have one : " <> (show inArgs))
                                        else
                                            if (toTree || toSoftWired || toHardWired) && (toDynamic || toStaticApprox)
                                                then
                                                    errorWithoutStackTrace
                                                        ("Multiple transform operations in transform (e.g. toTree, staticApprox)--can only have one at a time: " <> (show inArgs))
                                                else
                                                    let pruneEdges = False
                                                        warnPruneEdges = False
                                                        startVertex ∷ Maybe a
                                                        startVertex = Nothing

                                                        -- setup parallel
                                                        reoptimizeAction ∷ GlobalSettings → ProcessedData → Bool → Bool → Maybe Int → SimpleGraph → PhyG ReducedPhylogeneticGraph
                                                        reoptimizeAction = T.multiTraverseFullyLabelGraphReduced
                                                    in  -- transform nets to tree
                                                        if toTree
                                                            then -- already Tree return

                                                                if (graphType inGS == Tree)
                                                                    then do
                                                                        pure (inGS, origData, inData, inGraphList)
                                                                    else
                                                                        let newGS = inGS{graphType = Tree}

                                                                            -- generate and return display trees-- displayTreNUm / graph
                                                                            contractIn1Out1Nodes = True
                                                                        in do
                                                                            displayGraphList <-
                                                                                if chooseFirst
                                                                                    then pure $ fmap (take (fromJust numDisplayTrees) . (LG.generateDisplayTrees) contractIn1Out1Nodes) (fmap fst5 inGraphList)
                                                                                    else mapM (LG.generateDisplayTreesRandom rSeed (fromJust numDisplayTrees)) (fmap fst5 inGraphList)

                                                                            -- prob not required
                                                                            let displayGraphs = fmap GO.ladderizeGraph $ fmap GO.renameSimpleGraphNodes (concat displayGraphList)
                                                                            -- reoptimize as Trees
                                                                            reoptimizePar ← getParallelChunkTraverse
                                                                            newPhylogeneticGraphList ← reoptimizePar (reoptimizeAction newGS inData pruneEdges warnPruneEdges startVertex) displayGraphs
                                                                            -- newPhylogeneticGraphList = PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced newGS inData pruneEdges warnPruneEdges startVertex) displayGraphs -- `using` PU.myParListChunkRDS
                                                                            pure (newGS, origData, inData, newPhylogeneticGraphList)
                                                            else -- transform to softwired

                                                                if toSoftWired
                                                                    then
                                                                        if (graphType inGS == SoftWired)
                                                                            then do
                                                                                pure (inGS, origData, inData, inGraphList)
                                                                            else
                                                                                let newGS = inGS{graphType = SoftWired}
                                                                                in  do
                                                                                        reoptimizePar ← getParallelChunkTraverse
                                                                                        newPhylogeneticGraphList ←
                                                                                            reoptimizePar (reoptimizeAction newGS inData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList)
                                                                                        -- PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced newGS inData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList)  -- `using` PU.myParListChunkRDS
                                                                                        pure (newGS, origData, inData, newPhylogeneticGraphList)
                                                                    else -- transform to hardwired

                                                                        if toHardWired
                                                                            then
                                                                                if (graphType inGS == HardWired)
                                                                                    then do
                                                                                        pure (inGS, origData, inData, inGraphList)
                                                                                    else
                                                                                        let newGS = inGS{graphType = HardWired, graphFactor = NoNetworkPenalty}
                                                                                        in  do
                                                                                                logWith LogInfo ("Changing GraphFactor to NoNetworkPenalty for HardWired graphs" <> "\n")
                                                                                                reoptimizePar ← getParallelChunkTraverse
                                                                                                newPhylogeneticGraphList ←
                                                                                                    reoptimizePar (reoptimizeAction newGS inData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList)
                                                                                                -- PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced newGS inData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList)  -- `using` PU.myParListChunkRDS
                                                                                                pure (newGS, origData, inData, newPhylogeneticGraphList)
                                                                            else -- roll back to dynamic data from static approx

                                                                                if toDynamic
                                                                                    then do
                                                                                        reoptimizePar ← getParallelChunkTraverse
                                                                                        newPhylogeneticGraphList ←
                                                                                            reoptimizePar (reoptimizeAction inGS origData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList)
                                                                                        -- PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced inGS origData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList) -- `using` PU.myParListChunkRDS
                                                                                        logWith
                                                                                            LogInfo
                                                                                            ( "Transforming data to dynamic: "
                                                                                                <> (show $ minimum $ fmap snd5 inGraphList)
                                                                                                <> " -> "
                                                                                                <> (show $ minimum $ fmap snd5 newPhylogeneticGraphList)
                                                                                                <> "\n"
                                                                                            )
                                                                                        pure (inGS, origData, origData, newPhylogeneticGraphList)
                                                                                    else -- transform to static approx--using first Tree

                                                                                        if toStaticApprox
                                                                                            then do
                                                                                                newData ← makeStaticApprox inGS False inData (head $ L.sortOn snd5 inGraphList)
                                                                                                reoptimizePar ← getParallelChunkTraverse
                                                                                                newPhylogeneticGraphList ←
                                                                                                    reoptimizePar (reoptimizeAction inGS newData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList)
                                                                                                -- PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced inGS newData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList) -- `using` PU.myParListChunkRDS

                                                                                                if null inGraphList
                                                                                                    then do
                                                                                                        logWith LogInfo ("No graphs to base static approximation on--skipping." <> "\n")
                                                                                                        pure (inGS, origData, origData, inGraphList)
                                                                                                    else do
                                                                                                        logWith
                                                                                                            LogInfo
                                                                                                            ( "Transforming data to staticApprox: "
                                                                                                                <> (show $ minimum $ fmap snd5 inGraphList)
                                                                                                                <> " -> "
                                                                                                                <> (show $ minimum $ fmap snd5 newPhylogeneticGraphList)
                                                                                                                <> "\n"
                                                                                                            )
                                                                                                        pure (inGS, origData, newData, newPhylogeneticGraphList)
                                                                                            else -- change weight values in charInfo and reoptimize
                                                                                            -- reweights both origData and inData so weighting doens't get undone by static approc to and from transfomrations

                                                                                                if reWeight
                                                                                                    then
                                                                                                        let newOrigData = reWeightData (fromJust weightValue) charTypeList nameList origData
                                                                                                            newData = reWeightData (fromJust weightValue) charTypeList nameList inData
                                                                                                        in  if isNothing weightValue
                                                                                                                then
                                                                                                                    errorWithoutStackTrace
                                                                                                                        ("Reweight value is not specified correcty. Must be a double (e.g. 1.2): " <> (show (snd $ head reweightBlock)))
                                                                                                                else do
                                                                                                                    logWith
                                                                                                                        LogInfo
                                                                                                                        ( "Reweighting types "
                                                                                                                            <> (show charTypeList)
                                                                                                                            <> " and/or characters "
                                                                                                                            <> (L.intercalate ", " $ fmap TL.unpack nameList)
                                                                                                                            <> " to "
                                                                                                                            <> (show $ fromJust weightValue)
                                                                                                                            <> "\n\tReoptimizing graphs"
                                                                                                                            <> "\n"
                                                                                                                        )
                                                                                                                    reoptimizePar ← getParallelChunkTraverse
                                                                                                                    newPhylogeneticGraphList ←
                                                                                                                        reoptimizePar (reoptimizeAction inGS newData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList)
                                                                                                                    -- PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced inGS newData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList) -- `using` PU.myParListChunkRDS
                                                                                                                    pure (inGS, newOrigData, newData, newPhylogeneticGraphList)
                                                                                                    else -- changes the softwired optimization algorithm--this really for experimental use

                                                                                                        if changeCompressionResolutions
                                                                                                            then
                                                                                                                if isNothing compressionValue
                                                                                                                    then
                                                                                                                        errorWithoutStackTrace
                                                                                                                            ( "CompressResolutions value is not specified correcty. Must be either 'True' or 'False': "
                                                                                                                                <> (show (snd $ head changeCompressionBlock))
                                                                                                                            )
                                                                                                                    else
                                                                                                                        let newMethod =
                                                                                                                                if fromJust compressionValue == "true"
                                                                                                                                    then True
                                                                                                                                    else
                                                                                                                                        if fromJust compressionValue == "false"
                                                                                                                                            then False
                                                                                                                                            else
                                                                                                                                                errorWithoutStackTrace
                                                                                                                                                    ( "CompressResolutions value is not specified correcty. Must be either 'True' or 'False': "
                                                                                                                                                        <> (show (snd $ head changeSoftwiredMethodBlock))
                                                                                                                                                    )
                                                                                                                        in  do
                                                                                                                                reoptimizePar ← getParallelChunkTraverse
                                                                                                                                newPhylogeneticGraphList ←
                                                                                                                                    reoptimizePar
                                                                                                                                        (reoptimizeAction (inGS{compressResolutions = newMethod}) origData pruneEdges warnPruneEdges startVertex)
                                                                                                                                        (fmap fst5 inGraphList)
                                                                                                                                -- PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced (inGS  {compressResolutions = newMethod}) origData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList)
                                                                                                                                if newMethod /= compressResolutions inGS
                                                                                                                                    then do
                                                                                                                                        logWith LogInfo ("Changing compressResolutions method to " <> (show newMethod) <> "\n")
                                                                                                                                        pure (inGS{compressResolutions = newMethod}, origData, inData, newPhylogeneticGraphList)
                                                                                                                                    else do
                                                                                                                                        pure (inGS{compressResolutions = newMethod}, origData, inData, inGraphList)
                                                                                                            else -- changes dynamicEpsilon error check factor

                                                                                                                if changeEpsilon
                                                                                                                    then
                                                                                                                        if isNothing epsilonValue
                                                                                                                            then
                                                                                                                                errorWithoutStackTrace
                                                                                                                                    ("DynamicEpsilon value is not specified correcty. Must be a double (e.g. 0.02): " <> (show (snd $ head changeEpsilonBlock)))
                                                                                                                            else do
                                                                                                                                logWith LogInfo ("Changing dynamicEpsilon factor to " <> (show $ fromJust epsilonValue) <> "\n")
                                                                                                                                pure (inGS{dynamicEpsilon = 1.0 + ((fromJust epsilonValue) * (fractionDynamic inGS))}, origData, inData, inGraphList)
                                                                                                                    else -- changes the softwired optimization algorithm--this really for experimental use

                                                                                                                        if changeGraphFactor
                                                                                                                            then
                                                                                                                                if isNothing newGraphFactor
                                                                                                                                    then
                                                                                                                                        errorWithoutStackTrace
                                                                                                                                            ( "GraphFactor value is not specified correcty. Must be either 'NoPenalty', 'W15'. 'W23', or 'PMDL': "
                                                                                                                                                <> (show (snd $ head changeGraphFactorBlock))
                                                                                                                                            )
                                                                                                                                    else
                                                                                                                                        let newMethod =
                                                                                                                                                if fromJust newGraphFactor == "nopenalty"
                                                                                                                                                    then NoNetworkPenalty
                                                                                                                                                    else
                                                                                                                                                        if fromJust newGraphFactor == "w15"
                                                                                                                                                            then Wheeler2015Network
                                                                                                                                                            else
                                                                                                                                                                if fromJust newGraphFactor == "w23"
                                                                                                                                                                    then Wheeler2023Network
                                                                                                                                                                    else
                                                                                                                                                                        if fromJust newGraphFactor == "pmdl"
                                                                                                                                                                            then PMDLGraph
                                                                                                                                                                            else
                                                                                                                                                                                errorWithoutStackTrace
                                                                                                                                                                                    ( "GraphFactor value is not specified correcty. Must be either 'NoPenalty', 'W15'. 'W23', or 'PMDL': "
                                                                                                                                                                                        <> (show (snd $ head changeGraphFactorBlock))
                                                                                                                                                                                    )
                                                                                                                                        in  do
                                                                                                                                                reoptimizePar ← getParallelChunkTraverse
                                                                                                                                                newPhylogeneticGraphList ←
                                                                                                                                                    reoptimizePar
                                                                                                                                                        (reoptimizeAction (inGS{graphFactor = newMethod}) origData pruneEdges warnPruneEdges startVertex)
                                                                                                                                                        (fmap fst5 inGraphList)
                                                                                                                                                -- PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced (inGS  {graphFactor = newMethod}) origData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList) -- `using` PU.myParListChunkRDS
                                                                                                                                                if newMethod /= graphFactor inGS
                                                                                                                                                    then do
                                                                                                                                                        logWith LogInfo ("Changing graphFactor method to " <> (show newMethod) <> "\n")
                                                                                                                                                        pure (inGS{graphFactor = newMethod}, origData, inData, newPhylogeneticGraphList)
                                                                                                                                                    else do
                                                                                                                                                        pure (inGS{graphFactor = newMethod}, origData, inData, inGraphList)
                                                                                                                            else -- changes graphsSteepest -- maximum number of graphs evaluated in paralell at each "steepest" phse in swpa dn netadd/delete

                                                                                                                                if changeGraphsSteepest
                                                                                                                                    then
                                                                                                                                        if isNothing newGraphsSteepest
                                                                                                                                            then
                                                                                                                                                errorWithoutStackTrace
                                                                                                                                                    ( "GraphsSteepest value is not specified correcty. Must be an Integer (e.g. 5): " <> (show (snd $ head changeGraphsSteepestBlock))
                                                                                                                                                    )
                                                                                                                                            else do
                                                                                                                                                logWith LogInfo ("Changing GraphsSteepest factor to " <> (show $ fromJust newGraphsSteepest) <> "\n")
                                                                                                                                                pure (inGS{graphsSteepest = fromJust newGraphsSteepest}, origData, inData, inGraphList)
                                                                                                                                    else -- changes the multiTraverse behavior for all graphs

                                                                                                                                        if changeMultiTraverse
                                                                                                                                            then
                                                                                                                                                if isNothing multiTraverseValue
                                                                                                                                                    then
                                                                                                                                                        errorWithoutStackTrace
                                                                                                                                                            ( "MultiTraverse value is not specified correcty. Must be either 'True' or 'False': "
                                                                                                                                                                <> (show (snd $ head changeMultiTraverseBlock))
                                                                                                                                                            )
                                                                                                                                                    else
                                                                                                                                                        let newMethod =
                                                                                                                                                                if fromJust multiTraverseValue == "true"
                                                                                                                                                                    then True
                                                                                                                                                                    else
                                                                                                                                                                        if fromJust multiTraverseValue == "false"
                                                                                                                                                                            then False
                                                                                                                                                                            else
                                                                                                                                                                                errorWithoutStackTrace
                                                                                                                                                                                    ( "MultiTraverse value is not specified correcty. Must be either 'True' or 'False': "
                                                                                                                                                                                        <> (show (snd $ head changeMultiTraverseBlock))
                                                                                                                                                                                    )
                                                                                                                                                        in  do
                                                                                                                                                                reoptimizePar ← getParallelChunkTraverse
                                                                                                                                                                newPhylogeneticGraphList ←
                                                                                                                                                                    reoptimizePar
                                                                                                                                                                        (reoptimizeAction (inGS{multiTraverseCharacters = newMethod}) origData pruneEdges warnPruneEdges startVertex)
                                                                                                                                                                        (fmap fst5 inGraphList)
                                                                                                                                                                -- PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced (inGS  {multiTraverseCharacters = newMethod}) origData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList) -- `using` PU.myParListChunkRDS
                                                                                                                                                                if newMethod /= multiTraverseCharacters inGS
                                                                                                                                                                    then
                                                                                                                                                                        let lengthChangeString =
                                                                                                                                                                                if null inGraphList
                                                                                                                                                                                    then ""
                                                                                                                                                                                    else (":" <> (show $ minimum $ fmap snd5 inGraphList) <> " -> " <> (show $ minimum $ fmap snd5 newPhylogeneticGraphList))
                                                                                                                                                                        in  do
                                                                                                                                                                                logWith LogInfo ("Changing multiTraverse to " <> (show newMethod) <> lengthChangeString <> "\n")
                                                                                                                                                                                pure (inGS{multiTraverseCharacters = newMethod}, origData, inData, newPhylogeneticGraphList)
                                                                                                                                                                    else do
                                                                                                                                                                        pure (inGS{multiTraverseCharacters = newMethod}, origData, inData, inGraphList)
                                                                                                                                            else -- changes the softwired optimization algorithm--this really for experimental use

                                                                                                                                                if changeSoftwiredMethod
                                                                                                                                                    then
                                                                                                                                                        if isNothing newSoftwiredMethod
                                                                                                                                                            then
                                                                                                                                                                errorWithoutStackTrace
                                                                                                                                                                    ( "SoftwiredMethod value is not specified correcty. Must be either 'Exhaustive' or 'ResolutionCache': "
                                                                                                                                                                        <> (show (snd $ head changeSoftwiredMethodBlock))
                                                                                                                                                                    )
                                                                                                                                                            else
                                                                                                                                                                let newMethod =
                                                                                                                                                                        if fromJust newSoftwiredMethod == "naive"
                                                                                                                                                                            then Naive
                                                                                                                                                                            else
                                                                                                                                                                                if fromJust newSoftwiredMethod == "exhaustive"
                                                                                                                                                                                    then Naive
                                                                                                                                                                                    else
                                                                                                                                                                                        if fromJust newSoftwiredMethod == "resolutioncache"
                                                                                                                                                                                            then ResolutionCache
                                                                                                                                                                                            else
                                                                                                                                                                                                errorWithoutStackTrace
                                                                                                                                                                                                    ( "SoftwiredMethod value is not specified correcty. Must be either 'Naive' or 'ResolutionCache': "
                                                                                                                                                                                                        <> (show (snd $ head changeSoftwiredMethodBlock))
                                                                                                                                                                                                    )
                                                                                                                                                                    newMethodString =
                                                                                                                                                                        if newMethod == ResolutionCache
                                                                                                                                                                            then "ResolutionCache"
                                                                                                                                                                            else "Exhaustive"
                                                                                                                                                                in  do
                                                                                                                                                                        reoptimizePar ← getParallelChunkTraverse

                                                                                                                                                                        newPhylogeneticGraphList ←
                                                                                                                                                                            reoptimizePar
                                                                                                                                                                                (reoptimizeAction (inGS{softWiredMethod = newMethod}) origData pruneEdges warnPruneEdges startVertex)
                                                                                                                                                                                (fmap fst5 inGraphList)
                                                                                                                                                                        -- PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced (inGS  {softWiredMethod = newMethod}) origData pruneEdges warnPruneEdges startVertex) (fmap fst5 inGraphList) -- `using` PU.myParListChunkRDS

                                                                                                                                                                        if newMethod /= softWiredMethod inGS
                                                                                                                                                                            then do
                                                                                                                                                                                logWith LogInfo ("Changing softwired optimization method to " <> newMethodString <> "\n")
                                                                                                                                                                                pure (inGS{softWiredMethod = newMethod}, origData, inData, newPhylogeneticGraphList)
                                                                                                                                                                            else do
                                                                                                                                                                                pure (inGS{softWiredMethod = newMethod}, origData, inData, inGraphList)
                                                                                                                                                    else -- changes outgroup

                                                                                                                                                        if reRoot
                                                                                                                                                            then
                                                                                                                                                                if isNothing outgroupValue
                                                                                                                                                                    then errorWithoutStackTrace ("Outgroup is not specified correctly. Must be a string (e.g. \"Name\"): " <> (snd $ head reRootBlock))
                                                                                                                                                                    else
                                                                                                                                                                        let newOutgroupName = TL.filter (/= '"') $ fromJust outgroupValue
                                                                                                                                                                            newOutgroupIndex = V.elemIndex newOutgroupName (fst3 origData)
                                                                                                                                                                        in  do
                                                                                                                                                                                reoptimizePar ← getParallelChunkTraverse
                                                                                                                                                                                newPhylogeneticGraphList ←
                                                                                                                                                                                    reoptimizePar
                                                                                                                                                                                        (reoptimizeAction inGS origData pruneEdges warnPruneEdges startVertex)
                                                                                                                                                                                        (fmap (LG.rerootTree (fromJust newOutgroupIndex)) $ fmap fst5 inGraphList)
                                                                                                                                                                                -- PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced inGS origData pruneEdges warnPruneEdges startVertex) (fmap (LG.rerootTree (fromJust newOutgroupIndex)) $ fmap fst5 inGraphList)
                                                                                                                                                                                if isNothing newOutgroupIndex
                                                                                                                                                                                    then errorWithoutStackTrace ("Outgoup name not found: " <> (snd $ head reRootBlock))
                                                                                                                                                                                    else do
                                                                                                                                                                                        reoptimizePar ← getParallelChunkTraverse
                                                                                                                                                                                        newPhylogeneticGraphList ←
                                                                                                                                                                                            reoptimizePar
                                                                                                                                                                                                (reoptimizeAction inGS origData pruneEdges warnPruneEdges startVertex)
                                                                                                                                                                                                (fmap (LG.rerootTree (fromJust newOutgroupIndex)) $ fmap fst5 inGraphList)
                                                                                                                                                                                        -- PU.seqParMap (parStrategy $ strictParStrat inGS) (T.multiTraverseFullyLabelGraphReduced inGS origData pruneEdges warnPruneEdges startVertex) (fmap (LG.rerootTree (fromJust newOutgroupIndex)) $ fmap fst5 inGraphList)
                                                                                                                                                                                        logWith LogInfo ("Changing outgroup to " <> (TL.unpack newOutgroupName) <> "\n")
                                                                                                                                                                                        pure
                                                                                                                                                                                            (inGS{outgroupIndex = fromJust newOutgroupIndex, outGroupName = newOutgroupName}, origData, inData, newPhylogeneticGraphList)
                                                                                                                                                            else -- changes unionThreshold error check factor

                                                                                                                                                                if changeUnionThreshold
                                                                                                                                                                    then
                                                                                                                                                                        if isNothing unionValue
                                                                                                                                                                            then
                                                                                                                                                                                errorWithoutStackTrace
                                                                                                                                                                                    ("UninThreshold value is not specified correcty. Must be a double (e.g. 1.17): " <> (show (snd $ head changeUnionBlock)))
                                                                                                                                                                            else do
                                                                                                                                                                                logWith LogInfo ("Changing uninoTHreshold factor to " <> (show $ fromJust unionValue) <> "\n")
                                                                                                                                                                                pure (inGS{dynamicEpsilon = (fromJust unionValue)}, origData, inData, inGraphList)
                                                                                                                                                                    else -- modify the use of Network Add heurisitcs in network optimization

                                                                                                                                                                        if (fst $ head lcArgList) == "usenetaddheuristic"
                                                                                                                                                                            then
                                                                                                                                                                                let localCriterion
                                                                                                                                                                                        | ((snd $ head lcArgList) == "true") = True
                                                                                                                                                                                        | ((snd $ head lcArgList) == "false") = False
                                                                                                                                                                                        | otherwise =
                                                                                                                                                                                            errorWithoutStackTrace
                                                                                                                                                                                                ("Error in 'transform' command. UseNetAddHeuristic '" <> (snd $ head lcArgList) <> "' is not 'true' or 'false'")
                                                                                                                                                                                in  do
                                                                                                                                                                                        logWith LogInfo ("UseNetAddHeuristic set to " <> (snd $ head lcArgList) <> "\n")
                                                                                                                                                                                        pure (inGS{useNetAddHeuristic = localCriterion}, origData, inData, inGraphList)
                                                                                                                                                                            else error ("Transform type not implemented/recognized" <> (show inArgs))


-- | reWeightData sets weights to new values based on
reWeightData ∷ Double → [String] → [NameText] → ProcessedData → ProcessedData
reWeightData weightValue charTypeStringList charNameList (inName, inNameBV, inBlockDataV) =
    let charTypeList = concatMap stringToType charTypeStringList
        newBlockData = fmap (reweightBlockData weightValue charTypeList charNameList) inBlockDataV
    in  (inName, inNameBV, newBlockData)


-- |  stringToType takes  String and returns typelist
stringToType ∷ String → [CharType]
stringToType inString =
    if null inString
        then []
        else
            let inVal = fmap C.toLower inString
                typeList =
                    if inVal == "all"
                        then exactCharacterTypes <> sequenceCharacterTypes
                        else
                            if inVal == "prealigned"
                                then prealignedCharacterTypes
                                else
                                    if inVal `elem` ["nonexact", "dynamic"]
                                        then nonExactCharacterTypes
                                        else
                                            if inVal == "nonadditive"
                                                then [NonAdd, Packed2, Packed4, Packed5, Packed8, Packed64]
                                                else
                                                    if inVal == "additive"
                                                        then [Add]
                                                        else
                                                            if inVal == "matrix"
                                                                then [Matrix]
                                                                else
                                                                    if inVal == "sequence"
                                                                        then sequenceCharacterTypes
                                                                        else
                                                                            if inVal == "packed"
                                                                                then packedNonAddTypes
                                                                                else
                                                                                    if inVal == "packed2"
                                                                                        then [Packed2]
                                                                                        else
                                                                                            if inVal == "packed4"
                                                                                                then [Packed4]
                                                                                                else
                                                                                                    if inVal == "packed5"
                                                                                                        then [Packed5]
                                                                                                        else
                                                                                                            if inVal == "packed8"
                                                                                                                then [Packed8]
                                                                                                                else
                                                                                                                    if inVal == "packed64"
                                                                                                                        then [Packed64]
                                                                                                                        else
                                                                                                                            if inVal `elem` ["static", "exact", "qualitative"]
                                                                                                                                then exactCharacterTypes
                                                                                                                                else errorWithoutStackTrace ("Error in transform : Unrecognized character type '" <> inString <> "'")
            in  typeList


-- | reweightBlockData applies new weight to catagories of data
reweightBlockData ∷ Double → [CharType] → [NameText] → BlockData → BlockData
reweightBlockData weightValue charTypeList charNameList (blockName, blockData, charInfoV) =
    let newCharacterInfoV = fmap (reweightCharacterData weightValue charTypeList charNameList) charInfoV
    in  (blockName, blockData, newCharacterInfoV)


-- | reweightCharacterData changes weight in charInfo based on type or name
reweightCharacterData ∷ Double → [CharType] → [NameText] → CharInfo → CharInfo
reweightCharacterData weightValue charTypeList charNameList charInfo =
    let wildCardMatchCharName = filter (== True) $ fmap (textMatchWildcards (name charInfo)) charNameList
    in  -- trace ("RWC Wildcards: " <> (show $ fmap (textMatchWildcards (name charInfo)) charNameList)) (
        if null wildCardMatchCharName && (charType charInfo) `notElem` charTypeList
            then -- trace ("RWC not : " <> (show $ name charInfo) <> " of " <> (show charNameList) <> " " <> (show $ charType charInfo) <> " of " <> (show charTypeList))
                charInfo
            else -- trace ("RWC: " <> (show $ name charInfo) <> " " <> (show $ charType charInfo))
                charInfo{weight = weightValue}


-- )

{- | makeStaticApprox takes ProcessedData and returns static approx (implied alignment recoded) ProcessedData
if Tree take SA fields and recode appropriatrely given cost regeme of character
if Softwired--use display trees for SA
if hardWired--convert to softwired and use display trees for SA
since for heuristic searcing--uses additive weight for sequences and simple cost matrices, otherwise
matrix characters
-}
makeStaticApprox ∷ GlobalSettings → Bool → ProcessedData → ReducedPhylogeneticGraph → PhyG ProcessedData
makeStaticApprox inGS leavePrealigned inData@(nameV, nameBVV, blockDataV) inGraph =
    if LG.isEmpty (fst5 inGraph)
        then error "Empty graph in makeStaticApprox"
        else -- tree type

            if graphType inGS == Tree
                then do
                    let decGraph = thd5 inGraph

                    -- parallel setup
                    -- action :: Int -> BlockData
                    let action = pullGraphBlockDataAndTransform leavePrealigned decGraph blockDataV

                    pTraverse ← getParallelChunkMap
                    -- do each block in turn pulling and transforming data from inGraph
                    let newBlockDataV = pTraverse action [0 .. (length blockDataV - 1)]
                    -- PU.seqParMap (parStrategy $ strictParStrat inGS) (pullGraphBlockDataAndTransform leavePrealigned decGraph blockDataV) [0..(length blockDataV - 1)] -- `using` PU.myParListChunkRDS

                    if leavePrealigned
                        then do
                            pure (nameV, nameBVV, V.fromList newBlockDataV)

                    else do -- convert prealigned to non-additive if all 1's tcm

                        -- remove constants from new prealigned  this may be redundant since bit packing also removes constants
                        -- error here in case where there is missing seqeunce data for all but one input block for a character
                        -- newProcessedData ← R.removeConstantCharactersPrealigned (nameV, nameBVV, V.fromList newBlockDataV)

                        -- bit pack any new non-additive characters
                        newProcessedData' ← BP.packNonAdditiveData inGS (nameV, nameBVV, V.fromList newBlockDataV) -- newProcessedData

                        -- trace ("MSA:" <> (show (fmap (V.length . thd3) blockDataV, fmap (V.length . thd3) newBlockDataV)))
                        -- issues if no variation in block reducing length to zero so need leave "prealigned" if so
                        pure newProcessedData'

                else -- network static approx relies on display tree implied alignments after contacting out in=out=1 vertices
                -- harwired based on softwired optimization

                    if graphType inGS `elem` [SoftWired, HardWired]
                        then
                            let -- parallel setup
                                action ∷ (ProcessedData, SimpleGraph) → PhyG PhylogeneticGraph
                                action = T.multiTraverseFullyLabelGraphPair (inGS{graphType = Tree}) False False Nothing
                            in  do
                                    -- get display trees for each data block-- takes first of potentially multiple
                                    inFullGraph ←
                                        if graphType inGS == SoftWired
                                            then pure $ GO.convertReduced2PhylogeneticGraph inGraph
                                            else -- rediagnose HardWired as softwired
                                                T.multiTraverseFullyLabelGraph (inGS{graphType = SoftWired}) inData False False Nothing (fst5 inGraph)

                                    let blockDisplayList = fmap (GO.contractIn1Out1EdgesRename . GO.convertDecoratedToSimpleGraph . head) (fth6 inFullGraph)

                                    -- create seprate processed data for each block
                                    let blockProcessedDataList = fmap (CU.makeBlockData (fst3 inData) (snd3 inData)) (thd3 inData)

                                    pTraverse ← getParallelChunkTraverse
                                    decoratedBlockTreeList' ← pTraverse action (zip (V.toList blockProcessedDataList) (V.toList blockDisplayList))

                                    -- Perform full optimizations on display trees (as trees) with single block data (blockProcessedDataList) to create IAs
                                    let decoratedBlockTreeList = V.fromList decoratedBlockTreeList'
                                    -- V.fromList (zipWith (T.multiTraverseFullyLabelGraph' (inGS {graphType = Tree}) False False Nothing) (V.toList blockProcessedDataList) (V.toList blockDisplayList) `using` PU.myParListChunkRDS)

                                    -- get new processed (leaf) data
                                    let newBlockDataV = V.zipWith (getBlockLeafDataFromDisplayTree leavePrealigned) (fmap thd6 decoratedBlockTreeList) blockDataV

                                    if leavePrealigned
                                        then do
                                            pure (nameV, nameBVV, newBlockDataV)

                                    else do
                                        -- remove constants from new prealigned-- this may be redundant since bit packing also removes constants
                                        -- error here in case where there is missing seqeunce data for all but one input block for a character
                                        -- newProcessedData ← R.removeConstantCharactersPrealigned (nameV, nameBVV, newBlockDataV)

                                        -- bit pack any new non-additive characters
                                        newProcessedData' ← BP.packNonAdditiveData inGS (nameV, nameBVV, newBlockDataV) -- newProcessedData

                                        pure newProcessedData'
                        else do
                            logWith LogWarn ("Static Approx not yet implemented for graph type : " <> (show $ graphType inGS) <> " skipping" <> "\n")
                            pure inData


{- | getBlockLeafDataFromDisplayTree take a dispay tree and the block dat for that tree
and returns leae data--prealiged IA data
like pullGraphBlockDataAndTransform but for a single block and display tree
-}
getBlockLeafDataFromDisplayTree ∷ Bool → DecoratedGraph → BlockData → BlockData
getBlockLeafDataFromDisplayTree leavePrealigned inDecGraph blockCharInfo =
    let (_, leafVerts, _, _) = LG.splitVertexList inDecGraph
        (_, leafLabelList) = unzip leafVerts
        leafBlockData = fmap V.toList $ V.fromList $ fmap V.head $ fmap vertData leafLabelList

        -- new recoded data-- need to filter out constant chars after recoding
        -- need character length for missing values
        charLengthV = V.zipWith U.getMaxCharacterLength (thd3 blockCharInfo) leafBlockData

        (transformedLeafBlockData, transformedBlockInfo) = unzip $ fmap (transformData leavePrealigned (thd3 blockCharInfo) charLengthV) (fmap V.fromList $ V.toList leafBlockData)
    in  (fst3 blockCharInfo, V.fromList transformedLeafBlockData, head transformedBlockInfo)


{- | pullGraphBlockDataAndTransform takes a DecoratedGraph and block index and pulls
the character data of the block and transforms the leaf data by using implied alignment
feilds for dynamic characters
-}
pullGraphBlockDataAndTransform ∷ Bool → DecoratedGraph → V.Vector BlockData → Int → BlockData
pullGraphBlockDataAndTransform leavePrealigned inDecGraph blockCharInfoV blockIndex =
    let (_, leafVerts, _, _) = LG.splitVertexList inDecGraph
        (_, leafLabelList) = unzip leafVerts
        leafBlockData = fmap (V.! blockIndex) (fmap vertData leafLabelList)

        -- new recoded data-- need to filter out constant chars after recoding
        -- need character length for missing values
        charLengthV = V.zipWith U.getMaxCharacterLength (thd3 $ blockCharInfoV V.! blockIndex) (V.fromList $ fmap V.toList leafBlockData)

        (transformedLeafBlockData, transformedBlockInfo) = unzip $ fmap (transformData leavePrealigned (thd3 $ blockCharInfoV V.! blockIndex) charLengthV) leafBlockData
    in  -- trace ("PGDT: " <> show charLengthV)
        (fst3 $ blockCharInfoV V.! blockIndex, V.fromList transformedLeafBlockData, head transformedBlockInfo)


-- | transformData takes original Character info and character data and transforms to static if dynamic noting chracter type
transformData
    ∷ Bool → V.Vector CharInfo → V.Vector Int → V.Vector CharacterData → (V.Vector CharacterData, V.Vector CharInfo)
transformData leavePrealigned inCharInfoV inCharLengthV inCharDataV =
    if V.null inCharInfoV
        then (V.empty, V.empty)
        else
            let (outCharDataV, outCharInfoV) = V.unzip $ V.zipWith3 (transformCharacter leavePrealigned) inCharDataV inCharInfoV inCharLengthV
            in  (outCharDataV, outCharInfoV)


-- transformCharacter takes a single character info and character and returns IA if dynamic as is if not
-- checks if all gaps with the GV.filter.  If all gaps--it means the sequence char was missing and
-- implied alignment produced all gaps.  The passing of character length is not necessary when changed missing seq to empty
-- character--but leaving in case change back to []
-- "nonAddGap" not currently implemented
transformCharacter ∷ Bool → CharacterData → CharInfo → Int → (CharacterData, CharInfo)
transformCharacter leavePrealigned inCharData inCharInfo charLength =
    let inCharType = charType inCharInfo
        inCostMatrix = costMatrix inCharInfo
        alphSize = length $ alphabet inCharInfo

        -- determine if matrix is all same costs => nonadditive
        --                        all same except fort single indel costs => non add with gap binary chars
        --                        not either => matrix char
        (inCostMatrixType, gapCost) = R.getRecodingType inCostMatrix
    in  -- trace ("TC:" <> (show alphSize) <> " " <> (show $ alphabet inCharInfo)) (
        -- trace ("TC:" <> (show charLength) <> " " <> (show (GV.length $ snd3 $ slimAlignment inCharData, GV.length $ snd3 $ wideAlignment inCharData, GV.length $ snd3 $ hugeAlignment inCharData))) (
        if inCharType `elem` exactCharacterTypes
            then (inCharData, inCharInfo)
            else
                if inCharType `elem` prealignedCharacterTypes
                    then (inCharData, inCharInfo)
                    else -- trace ("TC: " <> inCostMatrixType) (
                    -- different types--vector wrangling
                    -- missing data fields set if no implied alignment ie missing data

                        if inCharType `elem` [SlimSeq, NucSeq]
                            then
                                let gapChar = (0 ∷ SlimState) `setBit` fromEnum gapIndex
                                    missingState = L.foldl' (setBit) (0 ∷ SlimState) [0 .. alphSize - 1]
                                    impliedAlignChar =
                                        if (not . GV.null $ GV.filter (/= gapChar) $ snd3 $ slimAlignment inCharData)
                                            then slimAlignment inCharData
                                            else
                                                let missingElement = SV.replicate charLength missingState -- if simple all ON then segfault do to lookup outside of cost matrix
                                                {-
                                                             impliedAlignChar = if (not . GV.null $ GV.filter (/= gapChar) $ snd3 $ slimAlignment inCharData) then slimAlignment inCharData
                                                                                else
                                                                                  let missingElement = SV.replicate charLength $ B.complement (0 :: SlimState) -- TRANS.setMissingBits (0 :: SlimState) 0 alphSize
                                                -}
                                                in  (missingElement, missingElement, missingElement)

                                    newPrelimBV = R.convert2BV 32 impliedAlignChar
                                    newPrelimBVGaps = addGaps2BV gapCost newPrelimBV
                                in  -- trace ("TC-Slim:" <> (show $ GV.length $ snd3 $ slimAlignment inCharData) <> " " <> (show $ snd3 $ impliedAlignChar)) (

                                    if leavePrealigned
                                        then (inCharData{alignedSlimPrelim = impliedAlignChar}, inCharInfo{charType = AlignedSlim})
                                        else
                                            if inCostMatrixType == "nonAdd"
                                                then (inCharData{stateBVPrelim = newPrelimBV}, inCharInfo{charType = NonAdd})
                                                else
                                                    if inCostMatrixType == "nonAddGap"
                                                        then (inCharData{stateBVPrelim = newPrelimBVGaps}, inCharInfo{charType = NonAdd})
                                                        else -- matrix recoding
                                                            (inCharData{alignedSlimPrelim = impliedAlignChar}, inCharInfo{charType = AlignedSlim})
                            else
                                if inCharType `elem` [WideSeq, AminoSeq]
                                    then
                                        let gapChar = (0 ∷ WideState) `setBit` fromEnum gapIndex
                                            missingState = L.foldl' (setBit) (0 ∷ WideState) [0 .. alphSize - 1]
                                            impliedAlignChar =
                                                if (not . GV.null $ GV.filter (/= gapChar) $ snd3 $ wideAlignment inCharData)
                                                    then wideAlignment inCharData
                                                    else
                                                        let missingElement = UV.replicate charLength missingState -- if simple all ON then segfault do to lookup outside of cost matrix
                                                        in  {-
                                                                         impliedAlignChar = if (not . GV.null $ GV.filter (/= gapChar) $ snd3 $ wideAlignment inCharData)  then wideAlignment inCharData
                                                                                            else
                                                                                              let missingElement = UV.replicate charLength $ B.complement (0 :: WideState) -- TRANS.setMissingBits (0 :: WideState) 0 alphSize
                                                            -}
                                                            (missingElement, missingElement, missingElement)

                                            newPrelimBV = R.convert2BV 64 impliedAlignChar
                                            newPrelimBVGaps = addGaps2BV gapCost newPrelimBV
                                        in  if leavePrealigned
                                                then (inCharData{alignedWidePrelim = impliedAlignChar}, inCharInfo{charType = AlignedWide})
                                                else
                                                    if inCostMatrixType == "nonAdd"
                                                        then (inCharData{stateBVPrelim = newPrelimBV}, inCharInfo{charType = NonAdd})
                                                        else
                                                            if inCostMatrixType == "nonAddGap"
                                                                then (inCharData{stateBVPrelim = newPrelimBVGaps}, inCharInfo{charType = NonAdd})
                                                                else -- matrix recoding
                                                                    (inCharData{alignedWidePrelim = impliedAlignChar}, inCharInfo{charType = AlignedWide})
                                    else
                                        if inCharType == HugeSeq
                                            then
                                                let gapChar = (BV.fromBits $ replicate alphSize False) `setBit` fromEnum gapIndex
                                                    missingState = BV.fromBits $ replicate alphSize True
                                                    impliedAlignChar =
                                                        if (not . GV.null $ GV.filter (/= gapChar) $ snd3 $ hugeAlignment inCharData)
                                                            then hugeAlignment inCharData
                                                            else
                                                                let missingElement = V.replicate alphSize missingState
                                                                in  (missingElement, missingElement, missingElement)

                                                    newPrelimBV = impliedAlignChar
                                                    newPrelimBVGaps = addGaps2BV gapCost newPrelimBV
                                                in  if leavePrealigned
                                                        then (inCharData{alignedHugePrelim = impliedAlignChar}, inCharInfo{charType = AlignedHuge})
                                                        else
                                                            if inCostMatrixType == "nonAdd"
                                                                then (inCharData{stateBVPrelim = newPrelimBV}, inCharInfo{charType = NonAdd})
                                                                else
                                                                    if inCostMatrixType == "nonAddGap"
                                                                        then (inCharData{stateBVPrelim = newPrelimBVGaps}, inCharInfo{charType = NonAdd})
                                                                        else -- matrix recoding
                                                                            (inCharData{alignedHugePrelim = impliedAlignChar}, inCharInfo{charType = AlignedHuge})
                                            else error ("Unrecognized character type in transformCharacter: " <> (show inCharType))


-- )

{- | addGaps2BV adds gap characters 0 = nonGap, 1 = Gap to Vector
of states to non-additive charcaters for static approx.  gapCost - 1 characters are added
sems wasteful, but comctant filtered out and recoded later when non-add/add charsa re optimized and bitpacked
since this only for leaves assume inM good for all
-}
addGaps2BV
    ∷ Int
    → (V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector)
    → (V.Vector BV.BitVector, V.Vector BV.BitVector, V.Vector BV.BitVector)
addGaps2BV gapCost (_, inM, _) =
    -- trace ("AG2BV: " <> (show inM)) (
    let gapChar = BV.fromNumber (BV.dimension $ V.head inM) (1 ∷ Int)
        noGap = L.replicate (gapCost - 1) $ BV.fromNumber (BV.dimension $ V.head inM) (1 ∷ Int)
        hasGap = L.replicate (gapCost - 1) $ BV.fromNumber (BV.dimension $ V.head inM) (2 ∷ Int)
        gapCharV = createGapChars inM gapChar [] noGap hasGap
        outM = inM <> gapCharV
    in  (outM, outM, outM)


-- )

{- | createGapChars takes a vector of bitvector coded states and checks if first states == 1 (= gap)
if so a number based on gap cost are created.. Will create n * original klength so need to
filter out constant characters later
-}
createGapChars
    ∷ V.Vector BV.BitVector → BV.BitVector → [BV.BitVector] → [BV.BitVector] → [BV.BitVector] → V.Vector BV.BitVector
createGapChars origBVV gapCharacter newCharL noGapL hasGapL =
    if V.null origBVV
        then V.fromList newCharL
        else
            if V.head origBVV == gapCharacter
                then createGapChars (V.tail origBVV) gapCharacter (hasGapL <> newCharL) noGapL hasGapL
                else createGapChars (V.tail origBVV) gapCharacter (noGapL <> newCharL) noGapL hasGapL
