{- |
Module specifying graph building functions

Distance builds (@Wagner@, @NJ@, @WPGMA@) are imported from @Wag2020@:

https://www.github.com/wardwheeler/wag2020
-}
module Search.Build (
    buildGraph,
) where

import Commands.Verify qualified as VER
import Control.Evaluation
import Control.Monad (when)
import Control.Monad.Logger (LogLevel (..), Logger (..))
import Data.Char
import Data.Foldable (fold)
import Data.Functor (($>))
import Data.List qualified as L
import Data.List.NonEmpty (NonEmpty ((:|)))
import Data.Text.Lazy qualified as TL
import Data.Vector qualified as V
import GeneralUtilities
import GraphOptimization.Traversals qualified as T
import Graphs.GraphOperations qualified as GO
import ParallelUtilities qualified as PU
import Reconciliation.ReconcileGraphs qualified as R
import Search.DistanceMethods qualified as DM
import Search.DistanceWagner qualified as DW
import Search.WagnerBuild qualified as WB
import SymMatrix qualified as M
import System.ErrorPhase (ErrorPhase (..))
import Text.Read
import Types.Types
import Utilities.DistanceUtilities qualified as DU
import Utilities.Distances qualified as DD
import Utilities.LocalGraph qualified as LG
import Utilities.Utilities qualified as U


{- | buildGraph wraps around build tree--build trees and adds network edges after build if network
with appropriate options
transforms graph type to Tree for builds then back to initial graph type
-}
buildGraph ∷ [Argument] → GlobalSettings → ProcessedData → [[VertexCost]] → Int → PhyG [ReducedPhylogeneticGraph]
buildGraph inArgs inGS inData pairwiseDistances rSeed =
    let getKeyBy ∷ (Eq a) ⇒ (((a, b) → Bool) → [([Char], [Char])] → t) → a → t
        getKeyBy f key = f ((== key) . fst) lcArgList
        hasKey = getKeyBy any
        filterKey = fmap snd . getKeyBy filter

        fstArgList = fmap (fmap toLower . fst) inArgs
        sndArgList = fmap (fmap toLower . snd) inArgs
        lcArgList = zip fstArgList sndArgList
        checkCommandList = checkCommandArgs "build" fstArgList VER.buildArgList

        -- block build options including number of display trees to return
        doEUN' = hasKey "eun"
        doCUN' = hasKey "cun"
        doEUN = (not doEUN' && not doCUN') || doEUN'
        returnTrees' = hasKey "displaytrees"
        returnGraph' = hasKey "graph"
        returnRandomDisplayTrees' = hasKey "atrandom"
        returnFirst' = hasKey "first"
        buildDistance = hasKey "distance"

        -- temporary change (if needed) to build tree structures
        inputGraphType = graphType inGS
        treeGS = inGS{graphType = Tree}

        -- really only trees now--but maybe later if can ensure phylogenetic graph from recocnile
        (returnGraph, returnTrees)
            | (graphType inGS) == Tree = (False, True)
            | returnGraph' || returnTrees' = (returnGraph', returnTrees')
            | otherwise = (False, True)

        -- default to return reandom and overrides if both specified
        (returnRandomDisplayTrees, _)
            | returnRandomDisplayTrees' || returnFirst' = (returnRandomDisplayTrees', returnFirst')
            | otherwise = (True, False)
    in do
        --let processedDataList ∷ [ProcessedData]
        let processedDataList = U.getProcessDataByBlock True inData
        parwiseDistanceResult <- mapM DD.getPairwiseDistances processedDataList
    
        -- initial build of trees from combined data--or by blocks
        let initialBuild numDisplayTrees
                | hasKey "filter" = buildTree True inArgs treeGS inData pairwiseDistances rSeed
                | otherwise = do
                 -- removing taxa with missing data for block
                 logWith LogInfo "Block building initial graph(s)"
                 let distanceMatrixList ∷ [[[VertexCost]]]
                     distanceMatrixList
                       -- | buildDistance = PU.seqParMap PU.myStrategyHighLevel DD.getPairwiseDistances processedDataList
                       -- TODO
                         | buildDistance = parwiseDistanceResult
                         | otherwise = replicate (length processedDataList) []

                 blockTrees ← fmap fold . traverse (buildTree' True inArgs treeGS rSeed) $ zip distanceMatrixList processedDataList
                 -- blockTrees = concat (PU.myChunkParMapRDS (buildTree' True inArgs treeGS inputGraphType seed) (zip distanceMatrixList processedDataList))
                 -- reconcile trees and return graph and/or display trees (limited by numDisplayTrees) already re-optimized with full data set
                 returnGraphs ← reconcileBlockTrees rSeed blockTrees numDisplayTrees returnTrees returnGraph returnRandomDisplayTrees doEUN
                 -- seqParMap ∷ (Traversable t) ⇒ Strategy b → (a → b) → t a → t b
                 -- TODO
                 --pure $ PU.seqParMap PU.myStrategyHighLevel (T.multiTraverseFullyLabelGraphReduced inGS inData True True Nothing) returnGraphs
                 pure $ fmap (T.multiTraverseFullyLabelGraphReduced inGS inData True True Nothing) returnGraphs

        -- check for valid command options
        failWhen (not checkCommandList) $ "Unrecognized command in 'build': " <> show inArgs

        numDisplayTrees ← case filterKey "displaytrees" of
            [] → pure 10
            [x] → case readMaybe x ∷ Maybe Int of
                Just i → pure i
                Nothing → failParseKeyInteger "displayTree" x
            _ → failParseKeyDuplicates "displayTree" inArgs

        numReturnTrees ← case filterKey "return" of
            [] → pure 10
            [x] → case readMaybe x ∷ Maybe Int of
                Just i → pure i
                Nothing → failParseKeyInteger "return" x
            _ → failParseKeyDuplicates "return" inArgs

        -- initial build of trees from combined data--or by blocks
        firstGraphs' ← initialBuild numDisplayTrees

        -- this to allow 'best' to return more trees then later 'returned' and contains memory by letting other graphs go out of scope
        let firstGraphs
                | not $ hasKey "filter" = GO.selectGraphs Unique numReturnTrees 0.0 (-1) firstGraphs'
                | otherwise = firstGraphs'

        -- reporting info
        logWith LogMore $ getBuildLogMessage "Block" "returned" "graphs" firstGraphs
        if inputGraphType == Tree || hasKey "filter"
            then return firstGraphs
            else do
                logWith LogInfo $ unwords ["\tRediagnosing as", show $ graphType inGS]
                pure $
                    PU.seqParMap
                        PU.myStrategyHighLevel
                        (T.multiTraverseFullyLabelGraphReduced inGS inData False False Nothing)
                        (fmap fst5 firstGraphs) -- `using` PU.myParListChunkRDS


{- | reconcileBlockTrees takes a lists of trees (with potentially varying leave complement) and reconciled them
as per the arguments producing a set of displayTrees (ordered or resolved random), and/or the reconciled graph
all outputs are re-optimzed and ready to go
-}
reconcileBlockTrees ∷ Int → [ReducedPhylogeneticGraph] → Int → Bool → Bool → Bool → Bool → PhyG [SimpleGraph]
reconcileBlockTrees rSeed blockTrees numDisplayTrees returnTrees returnGraph returnRandomDisplayTrees doEUN =
    -- trace ("Reconcile producing " <> (show numDisplayTrees)) (
    let -- numLeaves = V.length $ fst3 inData
        -- fullLeafSet = zip [0..(numLeaves - 1)] (V.toList $ fst3 inData)
        simpleGraphList = fmap fst5 blockTrees
        -- fullLeafGraphList = fmap (E.makeProcessedGraph fullLeafSet) simpleGraphList

        reconcileArgList ∷ ∀ {a}. [(String, [a])]
        reconcileArgList =
            if doEUN
                then [("eun", []), ("vertexLabel:true", []), ("connect:True", [])]
                else [("cun", []), ("vertexLabel:true", []), ("connect:True", [])]

        -- create reconciled graph--NB may NOT be phylogenetic graph--time violations etc.
        reconciledGraphInitial = snd $ R.makeReconcileGraph VER.reconcileArgList reconcileArgList simpleGraphList

        -- ladderize, time consistent-ized, removed chained network edges, removed treenodes with all network edge children
        {-
        reconciledGraph' = GO.convertGeneralGraphToPhylogeneticGraph "correct" reconciledGraphInitial
        noChainedGraph = LG.removeChainedNetworkNodes False reconciledGraph'
        noTreeNodesWithAllNetChildren = LG.removeTreeEdgeFromTreeNodeWithAllNetworkChildren $ fromJust noChainedGraph
        contractedGraph = GO.contractIn1Out1EdgesRename noTreeNodesWithAllNetChildren
        reconciledGraph = GO.convertGeneralGraphToPhylogeneticGraph "correct" contractedGraph -}

        -- chained was separate in past, now in convertGeneralGraphToPhylogeneticGraph
        reconciledGraph = GO.convertGeneralGraphToPhylogeneticGraph True reconciledGraphInitial

        -- this for non-convertable graphs
        reconciledGraph'
            | not $ LG.isEmpty reconciledGraph = reconciledGraph
            | otherwise = reconciledGraphInitial

        displayGraphs'
            | not returnRandomDisplayTrees = take numDisplayTrees $ LG.generateDisplayTrees True reconciledGraph'
            | otherwise = LG.generateDisplayTreesRandom rSeed numDisplayTrees reconciledGraph'

        -- need this to fix up some graphs after other stuff changed
        displayGraphs = fmap (GO.convertGeneralGraphToPhylogeneticGraph True) displayGraphs'

        -- displayGraphs = fmap GO.ladderizeGraph $ fmap GO.renameSimpleGraphNodes displayGraphs'
        numNetNodes = length $ fth4 (LG.splitVertexList reconciledGraph)
        strNetNodes =
            fold
                [ "Reconciled graph has "
                , show numNetNodes
                , " network nodes hence up to 2^"
                , show numNetNodes
                , " display trees"
                ]
    in  do
            when (LG.isEmpty reconciledGraph && not returnTrees) $
                failWithPhase
                    Unifying
                    "\n\n\tError--reconciled graph could not be converted to phylogenetic graph.  Consider modifying block tree search options or returning display trees."

            if not (LG.isEmpty reconciledGraph) && not returnTrees
                then logWith LogMore (strNetNodes <> " for softwired network") $> [reconciledGraph]
                else
                    if not returnGraph && returnTrees
                        then pure displayGraphs
                        else
                            if not (LG.isEmpty reconciledGraph)
                                then logWith LogMore ("\n\t" <> strNetNodes) $> reconciledGraph : displayGraphs
                                else
                                    logWith
                                        LogWarn
                                        "\n\tWarning--reconciled graph could not be converted to phylogenetic graph.  Consider modifying block tree search options or performing standard builds."
                                        $> displayGraphs


-- ))

-- | buildTree' wraps build tree and changes order of arguments for mapping
buildTree' ∷ Bool → [Argument] → GlobalSettings → Int → ([[VertexCost]], ProcessedData) → PhyG [ReducedPhylogeneticGraph]
buildTree' simpleTreeOnly inArgs inGS rSeed (pairwiseDistances, inData) =
    buildTree simpleTreeOnly inArgs inGS inData pairwiseDistances rSeed


{- | buildTree takes build options and returns constructed graphList
simpleTreeOnly (for block build) returns a single best tree to reduce edges in
reconcile step
-}
buildTree ∷ Bool → [Argument] → GlobalSettings → ProcessedData → [[VertexCost]] → Int → PhyG [ReducedPhylogeneticGraph]
buildTree simpleTreeOnly inArgs inGS inData@(nameTextVect, _, _) pairwiseDistances rSeed =
    let getKeyBy ∷ (Eq a) ⇒ (((a, b) → Bool) → [([Char], [Char])] → t) → a → t
        getKeyBy f key = f ((== key) . fst) lcArgList
        hasKey = getKeyBy any
        filterKey = fmap snd . getKeyBy filter

        fstArgList = fmap (fmap toLower . fst) inArgs
        sndArgList = fmap (fmap toLower . snd) inArgs
        lcArgList = zip fstArgList sndArgList
        checkCommandList = checkCommandArgs "build" fstArgList VER.buildArgList

        buildDistance = hasKey "distance"
        buildCharacter = hasKey "character"

        -- character build
        performBuildCharacter ∷ Int → PhyG [ReducedPhylogeneticGraph]
        performBuildCharacter numReplicates =
            let treeList = WB.rasWagnerBuild inGS inData rSeed numReplicates
                treeList'
                    | simpleTreeOnly = GO.selectGraphs Best 1 0.0 (-1) treeList
                    | otherwise = treeList
            in  do
                    logWith LogMore $ getBuildLogMessage "Character" "yielded" "trees" treeList'
                    pure treeList'

        -- distance build
        performBuildDistance ∷ Int → Int → PhyG [ReducedPhylogeneticGraph]
        performBuildDistance numReplicates numToSave =
            -- do all options in line and add together for return tree list
            let outgroupElem = outgroupIndex inGS
                nameStringVect = fmap TL.unpack nameTextVect
                distMatrix = M.fromLists pairwiseDistances

                whenKey ∷ (Monoid p) ⇒ [Char] → p → p
                whenKey key value
                    | hasKey key = value
                    | otherwise = mempty

                refinement
                    | hasKey "tbr" = "tbr"
                    | hasKey "spr" = "spr"
                    | hasKey "otu" = "otu"
                    | otherwise = "none"

                {-
                treeList1 =
                    whenKey "rdwag" $
                        randomizedDistanceWagner
                            simpleTreeOnly
                            inGS
                            inData
                            nameStringVect
                            distMatrix
                            outgroupElem
                            numReplicates
                            rSeed
                            numToSave
                            refinement
                treeList2 = whenKey "dwag" . pure $ distanceWagner simpleTreeOnly inGS inData nameStringVect distMatrix outgroupElem refinement
                treeList3 = whenKey "nj" . pure $ neighborJoin simpleTreeOnly inGS inData nameStringVect distMatrix outgroupElem refinement
                treeList4 = whenKey "doWPGMA" . pure $ wPGMA simpleTreeOnly inGS inData nameStringVect distMatrix outgroupElem refinement
                -}

                
            in  do
                rdWagTree <- randomizedDistanceWagner simpleTreeOnly inGS inData nameStringVect distMatrix outgroupElem numReplicates rSeed numToSave refinement
                dWagTree <- distanceWagner simpleTreeOnly inGS inData nameStringVect distMatrix outgroupElem refinement
                njTree <- neighborJoin simpleTreeOnly inGS inData nameStringVect distMatrix outgroupElem refinement
                wPGMATree <-  wPGMA simpleTreeOnly inGS inData nameStringVect distMatrix outgroupElem refinement
                let treeList1 = if (hasKey "rdwag") then rdWagTree
                                else mempty
                let treeList2 = if (hasKey "dWag") then [dWagTree]
                                else mempty
                let treeList3 = if (hasKey "nj") then [njTree]
                                else mempty
                let treeList4 = if (hasKey "doWPGMA") then [wPGMATree]
                                else mempty
                let treeListFull = fold [treeList1, treeList2, treeList3, treeList4]
                logWith LogInfo "\tBuilding Distance Wagner"
                case treeListFull of
                    [] → failWithPhase Computing $ "Distance build is specified, but without any method: " <> show inArgs
                    xs → do
                        logWith LogMore $ getBuildLogMessage "Distance" "yielded" "trees" $ xs
                        if not simpleTreeOnly
                            then pure xs
                            else pure $ GO.selectGraphs Best 1 0.0 (-1) xs
    in  do
            failWhen (not checkCommandList) $ "Unrecognized command in 'build': " <> show inArgs
            failWhen (buildDistance && buildCharacter) $
                "Cannot specify both 'character' and 'distance' builds in same build command" <> show inArgs
            numReplicates ← case filterKey "replicates" of
                [] → pure 10
                [x] → case readMaybe x ∷ Maybe Int of
                    Just i → pure i
                    Nothing → failParseKeyInteger "replicates" x
                _ → failParseKeyDuplicates "replicate" inArgs
            numToSave ← case filterKey "best" of
                [] → pure numReplicates
                [x] → case readMaybe x ∷ Maybe Int of
                    Just i → pure i
                    Nothing → failParseKeyInteger "best" x
                _ → failParseKeyDuplicates "best" inArgs

            if buildDistance
                then performBuildDistance numReplicates numToSave
                else performBuildCharacter numReplicates


{- | distanceWagner takes Processed data and pairwise distance matrix and returns
'best' addition sequence Wagner (defined in Farris, 1972) as fully decorated tree (as Graph)
-}
distanceWagner
    ∷ Bool → GlobalSettings → ProcessedData → V.Vector String → M.Matrix Double → Int → String → PhyG ReducedPhylogeneticGraph
distanceWagner simpleTreeOnly inGS inData leafNames distMatrix outgroupValue refinement =
    do
    distWagTreeList <- DM.doWagnerS leafNames distMatrix "closest" outgroupValue "best" 1 []
    let distWagTree = head distWagTreeList
    let distWagTree' = head $ DW.performRefinement refinement "best:1" "first" leafNames outgroupValue distWagTree
    let distWagTreeSimpleGraph = DU.convertToDirectedGraphText leafNames outgroupValue (snd4 distWagTree')
    let charInfoVV = V.map thd3 $ thd3 inData
    if not simpleTreeOnly
        then
            return $ T.multiTraverseFullyLabelGraphReduced
                inGS
                inData
                False
                False
                Nothing
                (GO.renameSimpleGraphNodes $ GO.dichotomizeRoot outgroupValue $ LG.switchRootTree (length leafNames) distWagTreeSimpleGraph)
        else
            let simpleWag = GO.renameSimpleGraphNodes $ GO.dichotomizeRoot outgroupValue $ LG.switchRootTree (length leafNames) distWagTreeSimpleGraph
            in  return $ (simpleWag, 0.0, LG.empty, V.empty, charInfoVV)


{- | randomizedDistanceWagner takes Processed data and pairwise distance matrix and returns
random addition sequence Wagner trees fully decorated tree (as Graph)
-}
randomizedDistanceWagner
    ∷ Bool
    → GlobalSettings
    → ProcessedData
    → V.Vector String
    → M.Matrix Double
    → Int
    → Int
    → Int
    → Int
    → String
    → PhyG [ReducedPhylogeneticGraph]
randomizedDistanceWagner simpleTreeOnly inGS inData leafNames distMatrix outgroupValue numReplicates rSeed numToKeep refinement =
    do
    let randomizedAdditionSequences = V.fromList <$> shuffleInt rSeed numReplicates [0 .. (length leafNames - 1)]
    randomizedAdditionWagnerTreeList <- DM.doWagnerS leafNames distMatrix "random" outgroupValue "random" numToKeep randomizedAdditionSequences
    let randomizedAdditionWagnerTreeList' = take numToKeep $ L.sortOn thd4 randomizedAdditionWagnerTreeList
    let randomizedAdditionWagnerTreeList'' =
            head <$>
               {- PU.seqParMap
               --     PU.myStrategyHighLevel
               TODO
               -}
                    fmap
                    (DW.performRefinement refinement "best:1" "first" leafNames outgroupValue)
                    randomizedAdditionWagnerTreeList'
        randomizedAdditionWagnerSimpleGraphList = fmap (DU.convertToDirectedGraphText leafNames outgroupValue . snd4) randomizedAdditionWagnerTreeList''
    let charInfoVV = V.map thd3 $ thd3 inData
    if not simpleTreeOnly
            then -- fmap ((T.multiTraverseFullyLabelGraphReduced inGS inData False False Nothing . GO.renameSimpleGraphNodes . GO.dichotomizeRoot outgroupValue) . LG.switchRootTree (length leafNames)) randomizedAdditionWagnerSimpleGraphList `using` PU.myParListChunkRDS
                do
                return $ PU.seqParMap
                    PU.myStrategyHighLevel
                    ( ( T.multiTraverseFullyLabelGraphReduced inGS inData False False Nothing
                            . GO.renameSimpleGraphNodes
                            . GO.dichotomizeRoot outgroupValue
                      )
                        . LG.switchRootTree (length leafNames)
                    )
                    randomizedAdditionWagnerSimpleGraphList
            else
                let numTrees = length randomizedAdditionWagnerSimpleGraphList
                    -- simpleRDWagList = fmap (GO.dichotomizeRoot outgroupValue . LG.switchRootTree (length leafNames)) randomizedAdditionWagnerSimpleGraphList `using` PU.myParListChunkRDS
                    simpleRDWagList =
                        PU.seqParMap
                            PU.myStrategyHighLevel
                            (GO.dichotomizeRoot outgroupValue . LG.switchRootTree (length leafNames))
                            randomizedAdditionWagnerSimpleGraphList
                in  
                do 
                    return $ L.zip5
                        simpleRDWagList
                        (replicate numTrees 0.0)
                        (replicate numTrees LG.empty)
                        (replicate numTrees V.empty)
                        (replicate numTrees charInfoVV)


{- | neighborJoin takes Processed data and pairwise distance matrix and returns
Neighbor-Joining tree as fully decorated tree (as Graph)
-}
neighborJoin
    ∷ Bool → GlobalSettings → ProcessedData → V.Vector String → M.Matrix Double → Int → String → PhyG ReducedPhylogeneticGraph
neighborJoin simpleTreeOnly inGS inData leafNames distMatrix outgroupValue refinement =
    do
    njTree <- DM.neighborJoining leafNames distMatrix outgroupValue
    let njTree' = head $ DW.performRefinement refinement "best:1" "first" leafNames outgroupValue njTree
    let njSimpleGraph = DU.convertToDirectedGraphText leafNames outgroupValue (snd4 njTree')
    let charInfoVV = V.map thd3 $ thd3 inData
    if not simpleTreeOnly
        then 
            do
            return $ T.multiTraverseFullyLabelGraphReduced
                inGS
                inData
                False
                False
                Nothing
                (GO.renameSimpleGraphNodes $ GO.dichotomizeRoot outgroupValue $ LG.switchRootTree (length leafNames) njSimpleGraph)
        else
            do
            let simpleNJ = GO.dichotomizeRoot outgroupValue $ LG.switchRootTree (length leafNames) njSimpleGraph
            return $ (simpleNJ, 0.0, LG.empty, V.empty, charInfoVV)


{- | wPGMA takes Processed data and pairwise distance matrix and returns
WPGMA tree as fully decorated tree (as Graph)
since root index not nOTUs as with other tres--chanegd as with dWag and NJ to make consistent.
-}
wPGMA
    ∷ Bool → GlobalSettings → ProcessedData → V.Vector String → M.Matrix Double → Int → String → PhyG ReducedPhylogeneticGraph
wPGMA simpleTreeOnly inGS inData leafNames distMatrix outgroupValue refinement = 
    do
    wpgmaTree <- DM.wPGMA leafNames distMatrix outgroupValue
    let wpgmaTree' = head $ DW.performRefinement refinement "best:1" "first" leafNames outgroupValue wpgmaTree
    let wpgmaSimpleGraph = DU.convertToDirectedGraphText leafNames outgroupValue (snd4 wpgmaTree')
    let charInfoVV = V.map thd3 $ thd3 inData
    if not simpleTreeOnly
            then
                return $ T.multiTraverseFullyLabelGraphReduced
                    inGS
                    inData
                    False
                    False
                    Nothing
                    (GO.renameSimpleGraphNodes $ GO.dichotomizeRoot outgroupValue $ LG.switchRootTree (length leafNames) wpgmaSimpleGraph)
    else
        let simpleWPGMA = GO.dichotomizeRoot outgroupValue $ LG.switchRootTree (length leafNames) wpgmaSimpleGraph
        in  return (simpleWPGMA, 0.0, LG.empty, V.empty, charInfoVV)


{-
### ### ###
### Log Helper functions
### ### ###
-}

getBuildLogMessage ∷ (Show b, Ord b) ⇒ String → String → String → [(a, b, c, d, e)] → String
getBuildLogMessage pref verb obj list =
    let listLength = show $ length list
        costRangeStr = case list of
            [] → mempty
            x : xs →
                let costValues = snd5 <$> x :| xs
                in  "at cost range " <> show (minimum costValues, maximum costValues)
    in  unwords ["\t" <> pref, "build", verb, listLength, obj, costRangeStr]


failParse ∷ String → Evaluation env a
failParse = failWithPhase Parsing


failParseKeyDuplicates ∷ (Show a) ⇒ String → a → PhyG b
failParseKeyDuplicates key args = failParse $ fold ["Multiple '", key, "' number specifications in command--can have only one: ", show args]


failParseKeyInteger ∷ (Show a) ⇒ String → a → PhyG b
failParseKeyInteger key val = failParse $ fold ["Key '", key, "' specification in 'build' not an integer: ", show val]


failWhen ∷ Bool → String → PhyG ()
failWhen cond = when cond . failParse
