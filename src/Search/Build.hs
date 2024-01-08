{- |
Module specifying graph building functions

Distance builds (@Wagner@, @NJ@, @WPGMA@) are imported from @Wag2020@:

https://www.github.com/wardwheeler/wag2020
-}
module Search.Build (
    buildGraph,
) where

import Commands.Verify qualified as VER
import Control.Monad (replicateM, when)
import Control.Monad.Random.Class
import Data.Char
import Data.Foldable (fold)
import Data.Functor (($>), (<&>))
import Data.List qualified as L
import Data.List.NonEmpty (NonEmpty ((:|)))
import Data.Maybe
import Data.Text.Lazy qualified as TL
import Data.Vector qualified as V
import GeneralUtilities
import GraphOptimization.Traversals qualified as T
import Graphs.GraphOperations qualified as GO
import PHANE.Evaluation
import PHANE.Evaluation.ErrorPhase (ErrorPhase (..))
import PHANE.Evaluation.Logging (LogLevel (..), Logger (..))
import PHANE.Evaluation.Verbosity (Verbosity (..))
import Reconciliation.ReconcileGraphs qualified as R
import Search.DistanceMethods qualified as DM
import Search.DistanceWagner qualified as DW
import Search.WagnerBuild qualified as WB
import SymMatrix qualified as M
import Text.Read
import Types.DistanceTypes
import Types.Types
import Utilities.DistanceUtilities qualified as DU
import Utilities.Distances qualified as D
import Utilities.Distances qualified as DD
import Utilities.LocalGraph qualified as LG
import Utilities.Utilities qualified as U


{- | buildGraph wraps around build tree--build trees and adds network edges after build if network
with appropriate options
transforms graph type to Tree for builds then back to initial graph type
-}
buildGraph ∷ [Argument] → GlobalSettings → ProcessedData → PhyG [ReducedPhylogeneticGraph]
buildGraph inArgs inGS inData =
    let fstArgList = fmap (fmap toLower . fst) inArgs
        sndArgList = fmap (fmap toLower . snd) inArgs
        lcArgList = zip fstArgList sndArgList
        checkCommandList = checkCommandArgs "build" fstArgList VER.buildArgList
        failure = failWithPhase Parsing . unwords
    in  -- check for valid command options
        if not checkCommandList
            then failure [ "Unrecognized command in 'build':",  show inArgs ]
            else
                let -- block build options including number of display trees to return
                    buildBlock = filter ((== "block") . fst) lcArgList
                    displayBlock = filter ((== "displaytrees") . fst) lcArgList
                    numDisplayTrees
                        | length displayBlock > 1 =
                            errorWithoutStackTrace ("Multiple displayTree number specifications in command--can have only one: " <> show inArgs)
                        | null displayBlock = Just 10
                        | null (snd $ head displayBlock) = Just 10
                        | otherwise = readMaybe (snd $ head displayBlock) ∷ Maybe Int

                    returnList = filter ((== "return") . fst) lcArgList
                    numReturnTrees
                        | length returnList > 1 =
                            errorWithoutStackTrace ("Multiple 'return' number specifications in command--can have only one: " <> show inArgs)
                        | null returnList = Just (maxBound ∷ Int)
                        | null (snd $ head returnList) = Just (maxBound ∷ Int)
                        | otherwise = readMaybe (snd $ head returnList) ∷ Maybe Int

                    doEUN' = any ((== "eun") . fst) lcArgList
                    doCUN' = any ((== "cun") . fst) lcArgList
                    doEUN
                        | not doEUN' && not doCUN' = True
                        | otherwise = doEUN'
                        
                    returnTrees' = any ((== "displaytrees") . fst) lcArgList
                    returnGraph' = any ((== "graph") . fst) lcArgList
                    returnRandomDisplayTrees' = any ((== "atrandom") . fst) lcArgList
                    returnFirst' = any ((== "first") . fst) lcArgList
                    buildDistance = any ((== "distance") . fst) lcArgList

                    -- temprary change (if needed) to buyild tree structures
                    inputGraphType = graphType inGS
                    treeGS = inGS{graphType = Tree}

                    -- really only trees now--but maybe later if can ensure phylogenetic graph from recocnile
                    (returnGraph, returnTrees) = case graphType inGS of
                        Tree -> (False, True)
                        _ | returnGraph' || returnTrees' -> (returnGraph', returnTrees')
                        _ -> (False, True)

                    -- default to return reandom and overrides if both specified
                    (returnRandomDisplayTrees, _) =
                        if returnRandomDisplayTrees' || returnFirst'
                            then (returnRandomDisplayTrees', returnFirst')
                            else (True, False)

                    processedDataList = U.getProcessDataByBlock True inData

                    -- parallel setup
                    pairwiseAction ∷ ProcessedData → PhyG [[VertexCost]]
                    pairwiseAction = DD.getPairwiseDistances

                    buildAction ∷ (ProcessedData, [[VertexCost]]) → PhyG [ReducedPhylogeneticGraph]
                    buildAction = uncurry $ buildTree' True inArgs treeGS

                    traverseAction ∷ Bool → Bool → Maybe Int → SimpleGraph → PhyG ReducedPhylogeneticGraph
                    traverseAction = T.multiTraverseFullyLabelGraphReduced inGS inData -- False False Nothing
                in  do
                        when (not $ null buildBlock) $ logWith LogInfo "Block building initial graph(s)\n"
                        -- initial build of trees from combined data--or by blocks
                        -- distance calculation moved here to get out of main scope
                        firstGraphs' ←
                            if null buildBlock
                                then
                                    let simpleTreeOnly = False
                                    in  do
                                            pairwiseDistances' ← D.getPairwiseDistances inData
                                            buildTreeList ← buildTree simpleTreeOnly inArgs treeGS inData pairwiseDistances'
                                            -- logWith LogTech $ fold [ "BL:\t", show $ length buildTreeList, "\n" ]
                                            pure buildTreeList
                                else do
                                    -- removing taxa with missing data for block
                                    -- trace ("Block building initial graph(s)") $
                                    distanceMatrixList ←
                                        if buildDistance
                                            then getParallelChunkTraverse >>= \pTraverse -> pairwiseAction `pTraverse` processedDataList
                                            else pure $ replicate (length processedDataList) []
                                    blockList ← getParallelChunkTraverse >>= \pTraverse ->
                                          buildAction `pTraverse` zip processedDataList distanceMatrixList
                                    let blockTrees = concat blockList

                                    returnGraphs ←
                                        reconcileBlockTrees blockTrees (fromJust numDisplayTrees) returnTrees returnGraph returnRandomDisplayTrees doEUN

                                    getParallelChunkTraverse >>= \pTraverse -> traverseAction True True Nothing `pTraverse` returnGraphs

                        -- this to allow 'best' to return more trees then later 'returned' and contains memory by letting other graphs go out of scope
                        firstGraphs ← case buildBlock of
                            [] → GO.selectGraphs Unique (fromJust numReturnTrees) 0.0 firstGraphs'
                            _ → pure firstGraphs'

                        -- reporting info
                        let (costString, returnString) = case firstGraphs of
                                [] -> ("\t\tBlock build returned 0 graphs", "\t\tReturning 0 graphs")
                                g:gs ->
                                    let gCosts = snd5 <$> firstGraphs
                                        hiCost = maximum gCosts
                                        loCost = minimum gCosts
                                        suffix = unwords
                                            [ show $ length firstGraphs
                                            , "graphs at cost range"
                                            , fold [ "(", show loCost, ", ", show hiCost, ")" ]
                                            ]
                                    in  ("\tBlock build yielded" <> suffix, "\tReturning " <> suffix)


                        case numDisplayTrees of
                            Nothing -> failure [ "DisplayTrees specification in build not an integer:", show . snd $ head displayBlock ]
                            Just nDistplayTrees -> case numReturnTrees of
                                Nothing -> failure [ "Return number specifications in build not an integer:", show . snd $ head returnList ]
                                Just nReturnTrees -> do
                                    logWith LogInfo $ returnString <> "\n"
                                    case (inputGraphType, buildBlock) of
                                        (Tree, []) -> pure firstGraphs
                                        (_, x:xs) -> logWith LogInfo (costString <> "\n") $> firstGraphs
                                        (_, []) -> do
                                                logWith LogInfo $ unwords [ "\tRediagnosing as", show $ graphType inGS, "\n" ]
                                                getParallelChunkTraverse >>= \pTraverse ->
                                                    pTraverse (traverseAction False False Nothing) $ fst5 <$> firstGraphs


{- | reconcileBlockTrees takes a lists of trees (with potentially varying leave complement) and reconciled them
as per the arguments producing a set of displayTrees (ordered or resolved random), and/or the reconciled graph
all outputs are re-optimzed and ready to go
-}
reconcileBlockTrees ∷ [ReducedPhylogeneticGraph] → Int → Bool → Bool → Bool → Bool → PhyG [SimpleGraph]
reconcileBlockTrees blockTrees numDisplayTrees returnTrees returnGraph returnRandomDisplayTrees doEUN =
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

        -- parallel setup
        convertAction ∷ SimpleGraph → PhyG SimpleGraph
        convertAction = GO.convertGeneralGraphToPhylogeneticGraph True
    in  do
            -- create reconciled graph--NB may NOT be phylogenetic graph--time violations etc.
            reconciledGraphInitial' ← R.makeReconcileGraph VER.reconcileArgList reconcileArgList simpleGraphList
            let reconciledGraphInitial = snd reconciledGraphInitial'

            -- ladderize, time consistent-ized, removed chained network edges, removed treenodes with all network edge children
            reconciledGraph ← GO.convertGeneralGraphToPhylogeneticGraph True reconciledGraphInitial

            -- this for non-convertable graphs
            let reconciledGraph'
                    | not $ LG.isEmpty reconciledGraph = reconciledGraph
                    | otherwise = reconciledGraphInitial

            displayGraphs' <-
                    if not returnRandomDisplayTrees then pure $ take numDisplayTrees $ LG.generateDisplayTrees True reconciledGraph'
                    else LG.generateDisplayTreesRandom numDisplayTrees reconciledGraph'

            -- need this to fix up some graphs after other stuff changed
            -- displayGraphs <- mapM (GO.convertGeneralGraphToPhylogeneticGraph True) displayGraphs'
            displayGraphs ← getParallelChunkTraverse >>= \pTraverse ->
                convertAction `pTraverse` displayGraphs'

            -- displayGraphs = fmap GO.ladderizeGraph $ fmap GO.renameSimpleGraphNodes displayGraphs'
            let numNetNodes = length $ fth4 (LG.splitVertexList reconciledGraph)
            let strNetNodes =
                    fold
                        [ "Reconciled graph has "
                        , show numNetNodes
                        , " network nodes hence up to 2^"
                        , show numNetNodes
                        , " display trees"
                        ]

            when (LG.isEmpty reconciledGraph && not returnTrees) $
                error
                    "\n\n\tError--reconciled graph could not be converted to phylogenetic graph.  Consider modifying block tree search options or returning display trees."

            if not (LG.isEmpty reconciledGraph) && not returnTrees
                then logWith LogMore (strNetNodes <> " for softwired network\n") $> [reconciledGraph]
                else
                    if not returnGraph && returnTrees
                        then pure displayGraphs
                        else
                            if not (LG.isEmpty reconciledGraph)
                                then logWith LogMore ("\n\t" <> strNetNodes) $> reconciledGraph : displayGraphs
                                else
                                    logWith
                                        LogWarn
                                        "\n\tWarning--reconciled graph could not be converted to phylogenetic graph.  Consider modifying block tree search options or performing standard builds.\n"
                                        $> displayGraphs


{- |
'buildTree'' wraps build tree and changes order of arguments for mapping.
-}
buildTree' ∷ Bool → [Argument] → GlobalSettings → ProcessedData → [[VertexCost]] → PhyG [ReducedPhylogeneticGraph]
buildTree' simpleTreeOnly inArgs inGS inData pairwiseDistances =
    buildTree simpleTreeOnly inArgs inGS inData pairwiseDistances


{- |
'buildTree' takes build options and returns constructed graphList
simpleTreeOnly (for block build) returns a single best tree to reduce edges in
reconcile step.
-}
buildTree ∷ Bool → [Argument] → GlobalSettings → ProcessedData → [[VertexCost]] → PhyG [ReducedPhylogeneticGraph]
buildTree simpleTreeOnly inArgs inGS inData@(nameTextVect, _, _) pairwiseDistances =
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

        {-
        -- character build
        performBuildCharacter ∷ Int → PhyG [ReducedPhylogeneticGraph]
        performBuildCharacter numReplicates =
            let treeList = WB.rasWagnerBuild inGS inData numReplicates
                treeList'
                    | simpleTreeOnly = GO.selectGraphs Best 1 0.0 (-1) treeList
                    | otherwise = treeList
            in  do
                    logWith LogMore $ getBuildLogMessage "Character" "yielded" "trees" treeList'
                    pure treeList'
        -}

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
            in  {-
                treeList1 = if hasKey "rdwag" then randomizedDistanceWagner simpleTreeOnly inGS inData nameStringVect distMatrix outgroupElem numReplicates numToSave refinement
                            else pure []
                treeList2 = if hasKey  "dwag" then distanceWagner simpleTreeOnly inGS inData nameStringVect distMatrix outgroupElem refinement
                            else pure []
                treeList3 = if hasKey  "nj" then neighborJoin simpleTreeOnly inGS inData nameStringVect distMatrix outgroupElem refinement
                            else pure []
                treeList4 = if hasKey  "wpgma" then wPGMA simpleTreeOnly inGS inData nameStringVect distMatrix outgroupElem refinement
                            else pure []

                -- treeListFull = fold [treeList1, treeList2, treeList3, treeList4]
                treeListFull = treeList1 <> treeList2 <> treeList3 <> treeList4
                -}

                do
                    -- logWith LogInfo ("L455: " <> (show (numReplicates,numToSave)))
                    treeList1 ← {-# SCC buildTree_treeList1 #-}
                        if hasKey "rdwag"
                            then
                                randomizedDistanceWagner
                                    simpleTreeOnly
                                    inGS
                                    inData
                                    nameStringVect
                                    distMatrix
                                    outgroupElem
                                    numReplicates
                                    numToSave
                                    refinement
                            else pure []
                    treeList2 ←
                        if hasKey "dwag"
                            then distanceWagner simpleTreeOnly inGS inData nameStringVect distMatrix outgroupElem refinement
                            else pure []
                    treeList3 ←
                        if hasKey "nj"
                            then neighborJoin simpleTreeOnly inGS inData nameStringVect distMatrix outgroupElem refinement
                            else pure []
                    treeList4 ←
                        if hasKey "wpgma"
                            then wPGMA simpleTreeOnly inGS inData nameStringVect distMatrix outgroupElem refinement
                            else pure []

                    -- let treeListFull = fold [treeList1, treeList2, treeList3, treeList4]
                    let treeListFull = treeList1 <> treeList2 <> treeList3 <> treeList4
                    logWith LogInfo "\tBuilding Distance Tree\n"
                    case treeListFull of
                        [] → errorWithoutStackTrace $ "Distance build is specified, but without any method: " <> show inArgs
                        xs → do
                            logWith LogMore $ (getBuildLogMessage "Distance" "yielded" "trees" $ xs) <> "\n"
                            if not simpleTreeOnly
                                then pure xs
                                else GO.selectGraphs Best 1 0 xs
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
                else -- else performBuildCharacter numReplicates
                do
                    -- character build
                    treeList ← WB.rasWagnerBuild inGS inData numReplicates
                    treeList' ← GO.selectGraphs Best 1 0 treeList
                    if simpleTreeOnly
                        then do
                            logWith LogMore $ (getBuildLogMessage "Character" "yielded" "trees" treeList') <> "\n"
                            pure treeList'
                        else do
                            logWith LogMore $ (getBuildLogMessage "Character" "yielded" "trees" treeList) <> "\n"
                            pure treeList


{- | distanceWagner takes Processed data and pairwise distance matrix and returns
'best' addition sequence Wagner (defined in Farris, 1972) as fully decorated tree (as Graph)
-}
distanceWagner
    ∷ Bool → GlobalSettings → ProcessedData → V.Vector String → M.Matrix Double → Int → String → PhyG [ReducedPhylogeneticGraph]
distanceWagner simpleTreeOnly inGS inData leafNames distMatrix outgroupValue refinement =
    do
        distWagTreeList ← DM.doWagnerS inGS leafNames distMatrix "closest" outgroupValue "best" 1 []
        let distWagTree = head distWagTreeList
        dWagRefined ← DW.performRefinement refinement "best:1" "first" leafNames outgroupValue distWagTree
        let distWagTree' = head dWagRefined
        let distWagTreeSimpleGraph = DU.convertToDirectedGraphText leafNames outgroupValue (snd4 distWagTree')
        let charInfoVV = V.map thd3 $ thd3 inData
        if not simpleTreeOnly
            then do
                result ←
                    T.multiTraverseFullyLabelGraphReduced
                        inGS
                        inData
                        False
                        False
                        Nothing
                        (GO.renameSimpleGraphNodes $ GO.dichotomizeRoot outgroupValue $ LG.switchRootTree (length leafNames) distWagTreeSimpleGraph)
                return $ [result]
            else
                let simpleWag = GO.renameSimpleGraphNodes $ GO.dichotomizeRoot outgroupValue $ LG.switchRootTree (length leafNames) distWagTreeSimpleGraph
                in  return $ [(simpleWag, 0.0, LG.empty, V.empty, charInfoVV)]


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
    → String
    → PhyG [ReducedPhylogeneticGraph]
randomizedDistanceWagner simpleTreeOnly inGS inData leafNames distMatrix outgroupValue numReplicates numToKeep refinement = {-# SCC randomizedDistanceWagner_TOP_DEF #-}
    -- set up parallel structures
    let refineAction ∷ TreeWithData → PhyG [TreeWithData]
        refineAction = DW.performRefinement refinement "best:1" "first" leafNames outgroupValue

        traverseGraphAction ∷ SimpleGraph → PhyG ReducedPhylogeneticGraph
        traverseGraphAction g0 = {-# SCC randomizedDistanceWagner_traverseGraphAction #-} do
            let -- debugger :: (Logger m, Show a, Show b, Show c) => a -> LG.Gr b c -> m ()
                -- debugger n g = logWith LogTech $ "In: 'randomizedDistanceWagner.traverseGraphAction'\n  Graph [ G_" <> show n <> " ]:\n" <> LG.prettify g
            -- debugger 0 g0
            let g1 = LG.switchRootTree (length leafNames) g0
            -- debugger 1 g1
            let g2 = GO.dichotomizeRoot outgroupValue g1
            -- debugger 2 g2
            result@(g3,_,_,_,_) <- T.multiTraverseFullyLabelGraphReduced inGS inData False False Nothing g2
            -- debugger 3 g3
            pure result


        dichotomizeAction ∷ SimpleGraph → SimpleGraph
        dichotomizeAction = GO.dichotomizeRoot outgroupValue . (LG.switchRootTree (length leafNames))

        directedGraphAction ∷ TreeWithData → SimpleGraph
        directedGraphAction = DU.convertToDirectedGraphText leafNames outgroupValue . snd4

        leafIndexVec = V.generate (V.length leafNames) id
    in  do
            randomizedAdditionSequences ← replicateM numReplicates $ shuffleList leafIndexVec
            randomizedAdditionWagnerTreeList ←
                DM.doWagnerS inGS leafNames distMatrix "random" outgroupValue "random" numToKeep randomizedAdditionSequences

            let randomizedAdditionWagnerTreeList' = take numToKeep $ L.sortOn thd4 randomizedAdditionWagnerTreeList

            -- logWith LogInfo ("L567: " <> (show (numToKeep, length randomizedAdditionWagnerTreeList')) <> "\n")

            randomizedAdditionWagnerTreeList'' :: [TreeWithData] ← getParallelChunkTraverse >>= \pTraverse -> fmap fold $
                refineAction `pTraverse` randomizedAdditionWagnerTreeList'

            randomizedAdditionWagnerSimpleGraphList <- getParallelChunkMap <&> \pMap ->
                  directedGraphAction `pMap` randomizedAdditionWagnerTreeList''
            -- fmap (DU.convertToDirectedGraphText leafNames outgroupValue . snd4) randomizedAdditionWagnerTreeList''
            let charInfoVV = V.map thd3 $ thd3 inData

            if not simpleTreeOnly
                then getParallelChunkTraverse >>= \pTraverse -> traverseGraphAction `pTraverse` randomizedAdditionWagnerSimpleGraphList
                else {-
                     return $ PU.seqParMap
                         PU.myStrategyHighLevel
                         ( ( T.multiTraverseFullyLabelGraphReduced inGS inData False False Nothing
                                 . GO.renameSimpleGraphNodes
                                 . GO.dichotomizeRoot outgroupValue
                           )
                             . LG.switchRootTree (length leafNames)
                         )
                         randomizedAdditionWagnerSimpleGraphList
                     -}

                    let numTrees = length randomizedAdditionWagnerSimpleGraphList
                    in  -- simpleRDWagList = fmap (GO.dichotomizeRoot outgroupValue . LG.switchRootTree (length leafNames)) randomizedAdditionWagnerSimpleGraphList `using` PU.myParListChunkRDS
                        {-
                        simpleRDWagList =
                            PU.seqParMap
                                PU.myStrategyHighLevel
                                (GO.dichotomizeRoot outgroupValue . LG.switchRootTree (length leafNames))
                                randomizedAdditionWagnerSimpleGraphList
                        -}

                        do
                            simpleRDWagList <- getParallelChunkMap <&> \pMap ->
                                dichotomizeAction `pMap` randomizedAdditionWagnerSimpleGraphList
                            return $
                                L.zip5
                                    simpleRDWagList
                                    (replicate numTrees 0.0)
                                    (replicate numTrees LG.empty)
                                    (replicate numTrees V.empty)
                                    (replicate numTrees charInfoVV)


{- | neighborJoin takes Processed data and pairwise distance matrix and returns
Neighbor-Joining tree as fully decorated tree (as Graph)
-}
neighborJoin
    ∷ Bool → GlobalSettings → ProcessedData → V.Vector String → M.Matrix Double → Int → String → PhyG [ReducedPhylogeneticGraph]
neighborJoin simpleTreeOnly inGS inData leafNames distMatrix outgroupValue refinement =
    do
        njTree ← DM.neighborJoining leafNames distMatrix outgroupValue
        njRefined ← DW.performRefinement refinement "best:1" "first" leafNames outgroupValue njTree
        let njTree' = head njRefined
        let njSimpleGraph = DU.convertToDirectedGraphText leafNames outgroupValue (snd4 njTree')
        let charInfoVV = V.map thd3 $ thd3 inData
        if not simpleTreeOnly
            then do
                result ←
                    T.multiTraverseFullyLabelGraphReduced
                        inGS
                        inData
                        False
                        False
                        Nothing
                        (GO.renameSimpleGraphNodes $ GO.dichotomizeRoot outgroupValue $ LG.switchRootTree (length leafNames) njSimpleGraph)
                pure [result]
            else do
                let simpleNJ = GO.dichotomizeRoot outgroupValue $ LG.switchRootTree (length leafNames) njSimpleGraph
                pure [(simpleNJ, 0.0, LG.empty, V.empty, charInfoVV)]


{- | wPGMA takes Processed data and pairwise distance matrix and returns
WPGMA tree as fully decorated tree (as Graph)
since root index not nOTUs as with other tres--chanegd as with dWag and NJ to make consistent.
-}
wPGMA
    ∷ Bool → GlobalSettings → ProcessedData → V.Vector String → M.Matrix Double → Int → String → PhyG [ReducedPhylogeneticGraph]
wPGMA simpleTreeOnly inGS inData leafNames distMatrix outgroupValue refinement =
    do
        wpgmaTree ← DM.wPGMA leafNames distMatrix outgroupValue
        wpgmaRefined ← DW.performRefinement refinement "best:1" "first" leafNames outgroupValue wpgmaTree
        let wpgmaTree' = head wpgmaRefined
        let wpgmaSimpleGraph = DU.convertToDirectedGraphText leafNames outgroupValue (snd4 wpgmaTree')
        let charInfoVV = V.map thd3 $ thd3 inData
        if not simpleTreeOnly
            then do
                result ←
                    T.multiTraverseFullyLabelGraphReduced
                        inGS
                        inData
                        False
                        False
                        Nothing
                        (GO.renameSimpleGraphNodes $ GO.dichotomizeRoot outgroupValue $ LG.switchRootTree (length leafNames) wpgmaSimpleGraph)
                return $ [result]
            else
                let simpleWPGMA = GO.dichotomizeRoot outgroupValue $ LG.switchRootTree (length leafNames) wpgmaSimpleGraph
                in  return [(simpleWPGMA, 0.0, LG.empty, V.empty, charInfoVV)]


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
                in  "at cost range " <> show (minimum costValues, maximum costValues) <> "\n"
    in  unwords ["\t" <> pref, "build", verb, listLength, obj, costRangeStr]


failParse ∷ String → Evaluation env a
failParse = failWithPhase Parsing


failParseKeyDuplicates ∷ (Show a) ⇒ String → a → PhyG b
failParseKeyDuplicates key args = failParse $ fold ["Multiple '", key, "' number specifications in command--can have only one: ", show args]


failParseKeyInteger ∷ (Show a) ⇒ String → a → PhyG b
failParseKeyInteger key val = failParse $ fold ["Key '", key, "' specification in 'build' not an integer: ", show val]


failWhen ∷ Bool → String → PhyG ()
failWhen cond = when cond . failParse
