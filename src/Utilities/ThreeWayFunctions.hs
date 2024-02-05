{-
ToDo:
   Add parallel optimization overblocks and characters?
-}

{- |
Module specifying three way optimization functions for use in pre-order of
HardWired graphs and iterative pass-type optimization for Trees
-}
module Utilities.ThreeWayFunctions (
    threeMedianFinal,
    addGapsToChildren,
    threeWayGeneric,
) where

import Bio.DynamicCharacter (HugeState, extractMedians)
import Bio.DynamicCharacter.Element (SlimState, WideState)
import Data.Alphabet
import Data.BitVector.LittleEndian qualified as BV
import Data.Bits
import Data.List qualified as L
import Data.MetricRepresentation qualified as MR
import Data.TCM.Dense qualified as TCMD
import Data.Vector qualified as V
import Data.Vector.Generic qualified as GV
import Data.Vector.Storable qualified as SV
import Data.Vector.Unboxed qualified as UV
import GeneralUtilities
import GraphOptimization.Medians qualified as M
import Input.BitPack qualified as BP
import SymMatrix qualified as S
import Types.Types


{- | threeMedianFinal calculates a 3-median for data types in a single character
for dynamic characters this is done by 3 min-trees
taking best result.  Later true 3-way via ffi can be incorporated
first type of operation two parents and current node-- since prelim
has (left child preliminary, node preliminary, right child preliminary)
that information can be used if needed
since assumes in 2 out 1 only the node preliminary field is used
need to remove naked gaps from medians for dynamic characters --hence the extract medians call
-}
threeMedianFinal ∷ CharInfo → CharacterData → CharacterData → CharacterData → CharacterData
threeMedianFinal charInfo parent1 parent2 curNode =
    let localCharType = charType charInfo
    in  if localCharType == Add
            then
                let threeFinal = V.zipWith3 threeWayAdditive (rangeFinal parent1) (rangeFinal parent2) (snd3 $ rangePrelim curNode)
                in  curNode{rangeFinal = threeFinal}
            else
                if localCharType == NonAdd
                    then
                        let threeFinal = V.zipWith3 threeWayNonAdditive (stateBVFinal parent1) (stateBVFinal parent2) (snd3 $ stateBVPrelim curNode)
                        in  curNode{stateBVFinal = threeFinal}
                    else
                        if localCharType `elem` packedNonAddTypes
                            then
                                let threeFinal = BP.threeWayPacked localCharType (packedNonAddFinal parent1) (packedNonAddFinal parent2) (snd3 $ packedNonAddPrelim curNode)
                                in  curNode{packedNonAddFinal = threeFinal}
                            else
                                if localCharType == Matrix
                                    then
                                        let threeFinal =
                                                V.zipWith3
                                                    (threeWayMatrix (costMatrix charInfo))
                                                    (matrixStatesFinal parent1)
                                                    (matrixStatesFinal parent2)
                                                    (matrixStatesPrelim curNode)
                                        in  curNode{matrixStatesFinal = threeFinal}
                                    else
                                        if localCharType == AlignedSlim
                                            then
                                                let threeFinal =
                                                        M.getFinal3WaySlim (slimTCM charInfo) (alignedSlimFinal parent1) (alignedSlimFinal parent2) (snd3 $ alignedSlimPrelim curNode)
                                                in  curNode{alignedSlimFinal = threeFinal}
                                            else
                                                if localCharType == AlignedWide
                                                    then
                                                        let threeFinal =
                                                                M.getFinal3WayWideHuge
                                                                    (wideTCM charInfo)
                                                                    (alignedWideFinal parent1)
                                                                    (alignedWideFinal parent2)
                                                                    (snd3 $ alignedWidePrelim curNode)
                                                        in  curNode{alignedWideFinal = threeFinal}
                                                    else
                                                        if localCharType == AlignedHuge
                                                            then
                                                                let threeFinal =
                                                                        M.getFinal3WayWideHuge
                                                                            (hugeTCM charInfo)
                                                                            (alignedHugeFinal parent1)
                                                                            (alignedHugeFinal parent2)
                                                                            (snd3 $ alignedHugePrelim curNode)
                                                                in  curNode{alignedHugeFinal = threeFinal}
                                                            else
                                                                if (localCharType == SlimSeq) || (localCharType == NucSeq)
                                                                    then
                                                                        let threeFinal = threeWaySlim charInfo parent1 parent2 curNode
                                                                        in  curNode
                                                                                { slimFinal = extractMedians $ M.makeDynamicCharacterFromSingleVector threeFinal
                                                                                }
                                                                    else
                                                                        if (localCharType == WideSeq) || (localCharType == AminoSeq)
                                                                            then
                                                                                let threeFinal = threeWayWide charInfo parent1 parent2 curNode
                                                                                in  curNode
                                                                                        { wideFinal = extractMedians $ M.makeDynamicCharacterFromSingleVector threeFinal
                                                                                        }
                                                                            else
                                                                                if localCharType == HugeSeq
                                                                                    then
                                                                                        let threeFinal = threeWayHuge charInfo parent1 parent2 curNode
                                                                                        in  curNode
                                                                                                { hugeFinal = extractMedians $ M.makeDynamicCharacterFromSingleVector threeFinal
                                                                                                }
                                                                                    else error ("Unrecognized/implemented character type: " <> show localCharType)


-- | threeWayNonAdditive takes the union/intersection operation over 3 non additive states
threeWayNonAdditive ∷ BV.BitVector → BV.BitVector → BV.BitVector → BV.BitVector
threeWayNonAdditive inA inB inC =
    let intersection3 = (inA .&. inB) .&. inC
        intersectionAB = inA .&. inB
        intersectionAC = inA .&. inC
        intersectionBC = inB .&. inC
        union3 = (inA .|. inB) .|. inC
    in  if not (BV.isZeroVector intersection3)
            then intersection3
            else
                if not (BV.isZeroVector intersectionAB)
                    then intersectionAB .|. inC
                    else
                        if not (BV.isZeroVector intersectionAC)
                            then intersectionAC .|. inB
                            else
                                if not (BV.isZeroVector intersectionBC)
                                    then intersectionBC .|. inA
                                    else union3


{- | threeWayAdditive take three additive states and returns median
the idea is the interval between the minimum of all three maximum (min maxA, maxB, maxC)
and the maximum of all three minima (max minA, minB, minC)
ordered such that the fst of pair not greater than second
-}
threeWayAdditive ∷ (Int, Int) → (Int, Int) → (Int, Int) → (Int, Int)
threeWayAdditive (minA, maxA) (minB, maxB) (minC, maxC) =
    let minOfMaxs = minimum [maxA, maxB, maxC]
        maxOfMins = maximum [minA, minB, minC]
    in  if maxOfMins > minOfMaxs
            then (maxOfMins, minOfMaxs)
            else (minOfMaxs, maxOfMins)


{- | threeWayMatrix creates median best state vector from a traceback, since parents could conflict
on traceback does a minimization.
The final states of parents will have non-maximum costs and these compared
to the the child states with pointers to their children are set for
traceback from current node to child(ren) from preliminary assignment
since type assumes two children--they are both set to same value so if either left or right
is set later the process will be correct
-}
threeWayMatrix
    ∷ S.Matrix Int → V.Vector MatrixTriple → V.Vector MatrixTriple → V.Vector MatrixTriple → V.Vector MatrixTriple
threeWayMatrix inCostMatrix parent1 parent2 curNode =
    let numStates = S.rows inCostMatrix

        -- get the costs of each state for each node, for prents non-maximal cost will be final states
        parent1StatesCost = fmap fst3 parent1
        parent2StatesCost = fmap fst3 parent2
        curNodeStatesCost = fmap fst3 curNode

        -- get the minimum cost for each state given combinations of all three nodes and the min cost child state
        minCost3States = getMinStatePair inCostMatrix (maxBound ∷ StateCost) numStates parent1StatesCost parent2StatesCost curNodeStatesCost
    in  -- minStateCost = V.minimum $ fmap fst minCost3States
        -- finalStatesTriple = fmap (assignMinMaxCost minStateCost (maxBound :: StateCost)) minCost3States

        minCost3States


{- | getMinStatePair takes cost matrix and state costs (vector of Int) and returns best median cost state of child for that best cost
if either parent or child has maxbound cost then that state get max bound cost
-}
getMinStatePair
    ∷ S.Matrix Int
    → StateCost
    → Int
    → V.Vector StateCost
    → V.Vector StateCost
    → V.Vector StateCost
    → V.Vector (StateCost, [ChildStateIndex], [ChildStateIndex])
getMinStatePair inCostMatrix maxCost numStates p1CostV p2CostV curCostV =
    let f ∷ (Num a) ⇒ a → a → a → a
        f a b c = a + b + c

        range = [0 .. numStates - 1]

        bestMedianCost vec = getBestPairCost inCostMatrix maxCost numStates vec <$> range

        -- get costs to parents-- will assume parent costs are 0 or max
        bestMedianCostP1 = bestMedianCost p1CostV
        bestMedianCostP2 = bestMedianCost p2CostV

        -- get costs to single child via preliminary states
        medianChildCostPairVect =
            getBestPairCostAndState inCostMatrix maxCost numStates curCostV <$> range

        -- get 3 sum costs and best state value
        threeWayStateCostList = zipWith3 f bestMedianCostP1 bestMedianCostP2 (fmap fst medianChildCostPairVect)
        minThreeWayCost = minimum threeWayStateCostList

        finalStateCostL = zipWith (assignBestMax minThreeWayCost maxCost) threeWayStateCostList medianChildCostPairVect
    in  V.fromList finalStateCostL


{- |
assignBestMax checks 3-way median state cost and if minimum sets to that otherwise sets to max
double 2nd field for 2-child type asumption
-}
assignBestMax
    ∷ StateCost → StateCost → StateCost → (StateCost, [ChildStateIndex]) → (StateCost, [ChildStateIndex], [ChildStateIndex])
assignBestMax minCost maxCost stateCost (_, stateChildList)
    | stateCost == minCost = (minCost, stateChildList, stateChildList)
    | otherwise = (maxCost, stateChildList, stateChildList)


{- |
getBestPairCost gets the baest cost for a state to each of parent states--does not keep parent state
-}
getBestPairCost ∷ S.Matrix Int → StateCost → Int → V.Vector StateCost → Int → StateCost
getBestPairCost inCostMatrix maxCost numStates parentStateCostV medianStateIndex =
    let stateCost = V.minimum $ V.zipWith (g inCostMatrix maxCost medianStateIndex) parentStateCostV (V.fromList [0 .. (numStates - 1)])

        g ∷ (Eq a) ⇒ S.Matrix a → a → Int → a → Int → a
        g cM mC mS pC pS
            | pC == mC = mC
            | otherwise = cM S.! (mS, pS)
    in  stateCost


{- |
getBestPairCostAndState gets best pair of median state and chikd states based on preliminarr states of node
-}
getBestPairCostAndState ∷ S.Matrix Int → StateCost → Int → V.Vector StateCost → Int → (StateCost, [ChildStateIndex])
getBestPairCostAndState inCostMatrix maxCost numStates childStateCostV medianStateIndex =
    let g ∷ (Eq a) ⇒ S.Matrix a → a → Int → a → Int → (a, Int)
        g cM mC mS pC pS
            | pC == mC = (mC, pS)
            | otherwise = (cM S.! (mS, pS), pS)

        statecostV = V.zipWith (g inCostMatrix maxCost medianStateIndex) childStateCostV (V.fromList [0 .. (numStates - 1)])
        minStateCost = V.minimum $ fmap fst statecostV
        bestPairs = V.filter ((== minStateCost) . fst) statecostV
        bestChildStates = V.toList $ fmap snd bestPairs
    in  (minStateCost, L.sort bestChildStates)


{- | threeWaySlim take charInfo, 2 parents, and curNOde and creates 3 median via
1) 3 DO medians (choosing lowest cost median) ((p1,p2), cn), ((cn,p1), p2), and ((cn,p2), p1)
2) inserting gaps to make all 3 line up
3) creating 3-medians
4) choosing lowest cost median
No change adjust is False since the 3-way lookup shold include that factor when it returns cost
-}
threeWaySlim ∷ CharInfo → CharacterData → CharacterData → CharacterData → SV.Vector SlimState
threeWaySlim charInfo parent1 parent2 curNode =
    -- trace ("3WSlim: ") (
    let noChangeAdjust = False
        -- pairwise median structures
        p1p2 = M.getDOMedianCharInfo noChangeAdjust charInfo parent1 parent2
        p1cN = M.getDOMedianCharInfo noChangeAdjust charInfo parent1 curNode
        p2cN = M.getDOMedianCharInfo noChangeAdjust charInfo parent2 curNode

        -- get 3rd to pairwise
        p1p2cN = M.getDOMedianCharInfo noChangeAdjust charInfo p1p2 curNode
        p1cNp2 = M.getDOMedianCharInfo noChangeAdjust charInfo p1cN parent2
        p2cNp1 = M.getDOMedianCharInfo noChangeAdjust charInfo p2cN parent1

        (a1, b1, c1) = addGapsToChildren (slimGapped p1p2cN) (slimGapped p1p2)
        (median1, cost1) = get3WayGeneric (TCMD.lookupThreeway (slimTCM charInfo)) a1 b1 c1

        (a2, b2, c2) = addGapsToChildren (slimGapped p1cNp2) (slimGapped p1cN)
        (median2, cost2) = get3WayGeneric (TCMD.lookupThreeway (slimTCM charInfo)) a2 b2 c2

        (a3, b3, c3) = addGapsToChildren (slimGapped p2cNp1) (slimGapped p2cN)
        (median3, cost3) = get3WayGeneric (TCMD.lookupThreeway (slimTCM charInfo)) a3 b3 c3

        minCost = minimum [cost1, cost2, cost3]
    in  if cost1 == minCost
            then median1
            else
                if cost2 == minCost
                    then median2
                    else median3


-- )

{- | threeWayWide take charInfo, 2 parents, and curNOde and creates 3 median via
1) 3 DO medians (choosing lowest cost median) ((p1,p2), cn), ((cn,p1), p2), and ((cn,p2), p1)
2) inserting gaps to make all 3 line up
3) creating 3-medians
4) choosing lowest cost median
No change adjust is False since the 3-way lookup shold include that factor when it returns cost
-}
threeWayWide ∷ CharInfo → CharacterData → CharacterData → CharacterData → UV.Vector WideState
threeWayWide charInfo parent1 parent2 curNode =
    let noChangeAdjust = False
        -- pairwise median structures
        p1p2 = M.getDOMedianCharInfo noChangeAdjust charInfo parent1 parent2
        p1cN = M.getDOMedianCharInfo noChangeAdjust charInfo parent1 curNode
        p2cN = M.getDOMedianCharInfo noChangeAdjust charInfo parent2 curNode

        -- get 3rd to pairwise
        p1p2cN = M.getDOMedianCharInfo noChangeAdjust charInfo p1p2 curNode
        p1cNp2 = M.getDOMedianCharInfo noChangeAdjust charInfo p1cN parent2
        p2cNp1 = M.getDOMedianCharInfo noChangeAdjust charInfo p2cN parent1

        (a1, b1, c1) = addGapsToChildren (wideGapped p1p2cN) (wideGapped p1p2)
        (median1, cost1) = get3WayGeneric (MR.retreiveThreewayTCM (wideTCM charInfo)) a1 b1 c1

        (a2, b2, c2) = addGapsToChildren (wideGapped p1cNp2) (wideGapped p1cN)
        (median2, cost2) = get3WayGeneric (MR.retreiveThreewayTCM (wideTCM charInfo)) a2 b2 c2

        (a3, b3, c3) = addGapsToChildren (wideGapped p2cNp1) (wideGapped p2cN)
        (median3, cost3) = get3WayGeneric (MR.retreiveThreewayTCM (wideTCM charInfo)) a3 b3 c3

        minCost = minimum [cost1, cost2, cost3]
    in  if cost1 == minCost
            then median1
            else
                if cost2 == minCost
                    then median2
                    else median3


{- | threeWayHuge take charInfo, 2 parents, and curNOde and creates 3 median via
1) 3 DO medians (choosing lowest cost median) ((p1,p2), cn), ((cn,p1), p2), and ((cn,p2), p1)
2) inserting gaps to make all 3 line up
3) creating 3-medians
4) choosing lowest cost median
No change adjust is False since the 3-way lookup shold include that factor when it returns cost
-}
threeWayHuge ∷ CharInfo → CharacterData → CharacterData → CharacterData → V.Vector HugeState
threeWayHuge charInfo parent1 parent2 curNode =
    let noChangeAdjust = False
        -- pairwise median structures
        p1p2 = M.getDOMedianCharInfo noChangeAdjust charInfo parent1 parent2
        p1cN = M.getDOMedianCharInfo noChangeAdjust charInfo parent1 curNode
        p2cN = M.getDOMedianCharInfo noChangeAdjust charInfo parent2 curNode

        -- get 3rd to pairwise
        p1p2cN = M.getDOMedianCharInfo noChangeAdjust charInfo p1p2 curNode
        p1cNp2 = M.getDOMedianCharInfo noChangeAdjust charInfo p1cN parent2
        p2cNp1 = M.getDOMedianCharInfo noChangeAdjust charInfo p2cN parent1

        (a1, b1, c1) = addGapsToChildren (hugeGapped p1p2cN) (hugeGapped p1p2)
        (median1, cost1) = get3WayGeneric (MR.retreiveThreewayTCM (hugeTCM charInfo)) a1 b1 c1

        (a2, b2, c2) = addGapsToChildren (hugeGapped p1cNp2) (hugeGapped p1cN)
        (median2, cost2) = get3WayGeneric (MR.retreiveThreewayTCM (hugeTCM charInfo)) a2 b2 c2

        (a3, b3, c3) = addGapsToChildren (hugeGapped p2cNp1) (hugeGapped p2cN)
        (median3, cost3) = get3WayGeneric (MR.retreiveThreewayTCM (hugeTCM charInfo)) a3 b3 c3

        minCost = minimum [cost1, cost2, cost3]
    in  if cost1 == minCost
            then median1
            else
                if cost2 == minCost
                    then median2
                    else median3


{- | addGapsToChildren pads out "new" gaps based on identity--if not identical--adds a gap based on cost matrix size
importand node filed orders correct--has moved around
-}
addGapsToChildren ∷ (FiniteBits a, GV.Vector v a) ⇒ (v a, v a, v a) → (v a, v a, v a) → (v a, v a, v a)
addGapsToChildren (reGappedParentFinal, _, reGappedNodePrelim) (gappedLeftChild, gappedNodePrelim, gappedRightChild) =
    -- trace ("AG2C:") (
    let (reGappedLeft, reGappedRight) = slideRegap reGappedNodePrelim gappedNodePrelim gappedLeftChild gappedRightChild mempty mempty
    in  if (GV.length reGappedParentFinal /= GV.length reGappedLeft) || (GV.length reGappedParentFinal /= GV.length reGappedRight)
            then
                error
                    ( "Vectors not same length "
                        <> show (GV.length reGappedParentFinal, GV.length reGappedLeft, GV.length reGappedRight)
                    )
            else (reGappedParentFinal, reGappedLeft, reGappedRight)


-- )

{- | slideRegap takes two version of same vectors (1st and snd) one with additional gaps and if the two aren't equal then adds gaps
to the 3rd and 4th input vectors
-}
slideRegap ∷ (FiniteBits a, GV.Vector v a) ⇒ v a → v a → v a → v a → [a] → [a] → (v a, v a)
slideRegap reGappedNode gappedNode gappedLeft gappedRight newLeftList newRightList =
    let finalize = GV.fromList . reverse
    in  case GV.uncons reGappedNode of
            Nothing → (finalize newLeftList, finalize newRightList)
            Just (headRGN, tailRGN) →
                let gapState = (headRGN `xor` headRGN) `setBit` fromEnum gapIndex
                    nextCall = slideRegap tailRGN
                in  -- gap in reGappedNode, null gappedNode is gap at end of reGappedNode
                    -- can copmplete the remainder of the slide as gaps only
                    case GV.uncons gappedNode of
                        Nothing →
                            let gapList = replicate (GV.length reGappedNode) gapState
                                gapApplication = finalize . (gapList <>)
                            in  (gapApplication newLeftList, gapApplication newRightList)
                        Just (headGN, tailGN) →
                            if headRGN /= headGN
                                then nextCall gappedNode gappedLeft gappedRight (gapState : newLeftList) (gapState : newRightList)
                                else -- no "new gap"

                                    nextCall
                                        tailGN
                                        (GV.tail gappedLeft)
                                        (GV.tail gappedRight)
                                        (GV.head gappedLeft : newLeftList)
                                        (GV.head gappedRight : newRightList)


-- | get3WayGeneric takes thee vectors and produces a (median, cost) pair
get3WayGeneric ∷ (GV.Vector v e) ⇒ (e → e → e → (e, Word)) → v e → v e → v e → (v e, Word)
get3WayGeneric tcm in1 in2 in3 =
    let len = GV.length in1
        vt = V.generate len $ \i → tcm (in1 GV.! i) (in2 GV.! i) (in3 GV.! i)

        gen ∷ (GV.Vector v a) ⇒ V.Vector (a, b) → v a
        gen v = let med i = fst $ v V.! i in GV.generate len med

        add ∷ (Num b) ⇒ V.Vector (a, b) → b
        add = V.foldl' (\x e → x + snd e) 0
    in  (,) <$> gen <*> add $ vt


{-Not using this now--but could.  Would need to add Aligned Types-}

{- | threeWayGeneric take charInfo, 2 parents, and curNOde and creates 3 median via
1) 3 DO medians (choosing lowest cost median) ((p1,p2), cn), ((cn,p1), p2), and ((cn,p2), p1)
2) inserting gaps to make all 3 line up
3) creating 3-medians
4) choosing lowest cost median
No change adjust is False since the 3-way lookup shold include that factor when it returns cost
-}
threeWayGeneric ∷ CharInfo → CharacterData → CharacterData → CharacterData → CharacterData
threeWayGeneric charInfo parent1 parent2 curNode =
    let noChangeAdjust = False
        localCharType = charType charInfo
        -- pairwise medina structures
        p1p2 = M.getDOMedianCharInfo noChangeAdjust charInfo parent1 parent2
        p1cN = M.getDOMedianCharInfo noChangeAdjust charInfo parent1 curNode
        p2cN = M.getDOMedianCharInfo noChangeAdjust charInfo parent2 curNode

        -- get 3rd to pairwise
        p1p2cN = M.getDOMedianCharInfo noChangeAdjust charInfo p1p2 curNode
        p1cNp2 = M.getDOMedianCharInfo noChangeAdjust charInfo p1cN parent2
        p2cNp1 = M.getDOMedianCharInfo noChangeAdjust charInfo p2cN parent1

        (median1Slim, median1Wide, median1Huge, cost1) =
            if localCharType `elem` [SlimSeq, NucSeq]
                then
                    let (a, b, c) = addGapsToChildren (slimGapped p1p2cN) (slimGapped p1p2)
                        (median, cost) = get3WayGeneric (TCMD.lookupThreeway (slimTCM charInfo)) a b c
                    in  (median, mempty, mempty, cost)
                else
                    if localCharType `elem` [AminoSeq, WideSeq]
                        then
                            let (a, b, c) = addGapsToChildren (wideGapped p1p2cN) (wideGapped p1p2)
                                (median, cost) = get3WayGeneric (MR.retreiveThreewayTCM (wideTCM charInfo)) a b c
                            in  (mempty, median, mempty, cost)
                        else
                            if localCharType == HugeSeq
                                then
                                    let (a, b, c) = addGapsToChildren (hugeGapped p1p2cN) (hugeGapped p1p2)
                                        (median, cost) = get3WayGeneric (MR.retreiveThreewayTCM (hugeTCM charInfo)) a b c
                                    in  (mempty, mempty, median, cost)
                                else error ("Unrecognized character type: " <> show localCharType)

        (median2Slim, median2Wide, median2Huge, cost2) =
            if localCharType `elem` [SlimSeq, NucSeq]
                then
                    let (a, b, c) = addGapsToChildren (slimGapped p1cNp2) (slimGapped p1cN)
                        (median, cost) = get3WayGeneric (TCMD.lookupThreeway (slimTCM charInfo)) a b c
                    in  (median, mempty, mempty, cost)
                else
                    if localCharType `elem` [AminoSeq, WideSeq]
                        then
                            let (a, b, c) = addGapsToChildren (wideGapped p1cNp2) (wideGapped p1cN)
                                (median, cost) = get3WayGeneric (MR.retreiveThreewayTCM (wideTCM charInfo)) a b c
                            in  (mempty, median, mempty, cost)
                        else
                            if localCharType == HugeSeq
                                then
                                    let (a, b, c) = addGapsToChildren (hugeGapped p1cNp2) (hugeGapped p1cN)
                                        (median, cost) = get3WayGeneric (MR.retreiveThreewayTCM (hugeTCM charInfo)) a b c
                                    in  (mempty, mempty, median, cost)
                                else error ("Unrecognized character type: " <> show localCharType)

        (median3Slim, median3Wide, median3Huge, cost3) =
            if localCharType `elem` [SlimSeq, NucSeq]
                then
                    let (a, b, c) = addGapsToChildren (slimGapped p2cNp1) (slimGapped p2cN)
                        (median, cost) = get3WayGeneric (TCMD.lookupThreeway (slimTCM charInfo)) a b c
                    in  (median, mempty, mempty, cost)
                else
                    if localCharType `elem` [AminoSeq, WideSeq]
                        then
                            let (a, b, c) = addGapsToChildren (wideGapped p2cNp1) (wideGapped p2cN)
                                (median, cost) = get3WayGeneric (MR.retreiveThreewayTCM (wideTCM charInfo)) a b c
                            in  (mempty, median, mempty, cost)
                        else
                            if localCharType == HugeSeq
                                then
                                    let (a, b, c) = addGapsToChildren (hugeGapped p2cNp1) (hugeGapped p2cN)
                                        (median, cost) = get3WayGeneric (MR.retreiveThreewayTCM (hugeTCM charInfo)) a b c
                                    in  (mempty, mempty, median, cost)
                                else error ("Unrecognized character type: " <> show localCharType)

        minCost = minimum [cost1, cost2, cost3]
        (medianBestSlim, medianBestWide, medianBestHuge) =
            if cost1 == minCost
                then (median1Slim, median1Wide, median1Huge)
                else
                    if cost2 == minCost
                        then (median2Slim, median2Wide, median2Huge)
                        else (median3Slim, median3Wide, median3Huge)
    in  -- set for correct data type
        if localCharType `elem` [SlimSeq, NucSeq]
            then emptyCharacter{slimFinal = medianBestSlim}
            else
                if localCharType `elem` [AminoSeq, WideSeq]
                    then emptyCharacter{wideFinal = medianBestWide}
                    else
                        if localCharType == HugeSeq
                            then emptyCharacter{hugeFinal = medianBestHuge}
                            else error ("Unrecognized character type: " <> show localCharType)
