-----------------------------------------------------------------------------
-- |
-- Module      :  DirectOptimization.Pairwise.Test
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-- Test suite for dynamic characters
--
-----------------------------------------------------------------------------

module DirectOptimization.Pairwise.Test
  ( testSuite
  ) where

import           Bio.DynamicCharacter
import           Bio.DynamicCharacter.HandleGaps
import           Data.Foldable
import           Data.List                                  (intercalate, partition, sortBy)
import           Data.Map                                   (assocs, insertWith)
import           Data.Maybe
import           Data.MetricRepresentation
import           Data.Ord
import qualified Data.Vector.Generic                        as GV
import qualified Data.Vector.Storable                       as SV
import           Data.Word
import           Test.Aligners
import           Test.QuickCheck.Instances.DynamicCharacter
import           Test.Tasty
import           Test.Tasty.QuickCheck


-- |
-- The test-suite for pairwise string alignment algorithms.
testSuite :: TestTree
testSuite = testGroup "Pairwise alignment tests"
    [ knownErrorCases
    , characterSuite
    , gapSubstitutionSuite
    , validitySuite
    , consistencySuite
    ]


testPropertyRigorously :: Testable a => TestName -> a -> TestTree
testPropertyRigorously x = localOption (QuickCheckTests  1000) . testProperty x


testPropertyLeniently :: Testable a => TestName -> a -> TestTree
testPropertyLeniently  x = localOption (QuickCheckMaxRatio 32) . testProperty x


knownErrorCases :: TestTree
knownErrorCases = testGroup "Known Error Cases"
    [ specialCase1 "001" insertGapsDeleteGapsRightIdentity
        "H "
        "h-"
        " g"

    , specialCase2 "002" (commutivity $ withAlignmentλ 2)
        "T A y GR"
        "tcav-gGA"
        "C-V GNA "

        "T A y GR"
        "tcav-gGA"
        "C-V GNA "

    , specialCase2 "003" (commutivity $ withAlignmentλ 8)
        " AC"
        "-WY"
        "-TT"

        "VGT-C??-GC "
        "vgttcTV-SC-"
        " - T TV CC-"

    , specialCase2 "004" (commutivity $ withAlignmentλ 2)
        "AG"
        "Mv"
        "Cm"

        "Y-"
        "y-"
        "- "

    , specialCase2 "005" (commutivity $ withAlignmentλ 5)
        "R"
        "r"
        " "

        "C- "
        "Sc-"
        "GC-"

    , specialCase2 "006" (commutivity $ withAlignmentλ 5)
        "D?HvGGC?CY ACctAH-GC G"
        "NCT-ggSGMTtaSy-dNtgCtG"
        "CCT   GGATT GT kGT CTs"

        "-AAYG--TGGT TSGC  tATA "
        "aAaBGcaWkGWcTsVc--watMa"
        "AA G?CAAtGACw M d-A- CA"

    , specialCase2 "007" (commutivity $ withAlignmentλ 5)
        "H  T TCCssGC"
        "TgttgTcCsvgC"
        "TGT GT csA N"

        "AkkyN GmC-W?GTs?"
        "a--c?dg-c-WBgB-T"
        " - s D   ?wB S-T"

    , specialCase2 "008" (checkAlignmentResults $ withMetric 0)
        "GGG AGsAAH"
        "GDgtagvWAA"
        "GW T  ATaA"

        "TCN-GGTAAG"
        "tS?cSRTNAB"
        " G CCATBAY"

    , specialCase2 "009" (checkAlignmentResults $ withMetric 0)
        "GA?AHGGW T sA"
        "S??MhgBTcKsca"
        "Cb?C  YbCGSh-"

        "?c G Cd ACATG"
        "hCgRmC?cWMRYG"
        "hBGAMBCCTAGCG"
    
    , specialCase2 "010" (checkAlignmentResults $ withMetric 0)
        "TMCG ACATTCA"
        "WCCggNMaWTcR"
        "ACB GBA-AT G"

        "ACB GBA-AT G"
        "svagggb?aWTt"
        "Cr G- sA TT "

    , specialCase2 "011" (checkAlignmentResults $ withMetric 0)
        "TH TTTT ?"
        "TNtttKbtG"
        "TGT  GsTG"

        "-  A -TTK"
        "cgcaattWB"
        "CGC-AT AC"

    , specialCase2 "012" (checkAlignmentResults $ withMetric 0)
        "?G?TT Y  G"
        "TGNTYgytc?"
        "TGNyCG-TCh"

        "TCTTTTGMTG"
        "wcttTtRHtS"
        "a  -y AT-C"
    ]
  where
    withMetric     = (metricChoices !!)
    withAlignmentλ = snd . (alignmentMetricCombinations !!)


specialCase1
  :: String
  -> (DNA -> Property)
  -> String -> String -> String
  -> TestTree
specialCase1 str prop x y z =
    testProperty ("Special Case: " <> str) . prop $ buildDNA x y z


specialCase2
  :: String
  -> (DyadDNA -> Property)
  -> String -> String -> String
  -> String -> String -> String
  -> TestTree
specialCase2 str prop a b c x y z =
    testProperty ("Special Case: " <> str) . prop $ buildDNA a b c :×: buildDNA x y z


characterSuite :: TestTree
characterSuite = testGroup "Dynamic Character operations"
    [ contextQueries
    , missingQueries
    , transposeProperties
    ]


contextQueries :: TestTree
contextQueries = testGroup "Context Queries"
    [ testProperty "Alignment --> not Deletion " alignIsNotDelete
    , testProperty "Alignment --> not Insertion" alignIsNotInsert
    , testProperty "Alignment --> not Gapped   " alignIsNotGapped

    , testProperty "Deletion  --> not Alignment" deleteIsNotAlign
    , testProperty "Deletion  --> not Insertion" deleteIsNotInsert
    , testProperty "Deletion  --> not Gapped   " deleteIsNotGapped

    , testProperty "Insertion --> not Alignment" insertIsNotAlign
    , testProperty "Insertion --> not Deletion " insertIsNotDelete
    , testProperty "Insertion --> not Gapped   " insertIsNotGapped
    ]
  where
    alignIsNotDelete :: DNA -> Property
    alignIsNotDelete (DNA char) =
        isAlign char 0 ==> not (isDelete char 0)

    alignIsNotInsert :: DNA -> Property
    alignIsNotInsert (DNA char) =
        isAlign char 0 ==> not (isInsert char 0)

    alignIsNotGapped :: DNA -> Property
    alignIsNotGapped (DNA char) =
        isAlign char 0 ==> not (isGapped char 0)

    deleteIsNotAlign :: DNA -> Property
    deleteIsNotAlign (DNA char) =
        isDelete char 0 ==> not (isAlign char 0)

    deleteIsNotInsert :: DNA -> Property
    deleteIsNotInsert (DNA char) =
        isDelete char 0 ==> not (isInsert char 0)

    deleteIsNotGapped :: DNA -> Property
    deleteIsNotGapped (DNA char) =
        isDelete char 0 ==> not (isGapped char 0)

    insertIsNotAlign :: DNA -> Property
    insertIsNotAlign (DNA char) =
        isInsert char 0 ==> not (isAlign char 0)

    insertIsNotDelete :: DNA -> Property
    insertIsNotDelete (DNA char) =
        isInsert char 0 ==> not (isDelete char 0)

    insertIsNotGapped :: DNA -> Property
    insertIsNotGapped (DNA char) =
        isInsert char 0 ==> not (isGapped char 0)


missingQueries :: TestTree
missingQueries = testGroup "Missing Character Query"
    [ testPropertyLeniently "Missing --> length = 0" missingIsZeroLength
--    , testProperty "Missing --> transpose = id"
    ]
  where
    missingIsZeroLength :: DNA -> Property
    missingIsZeroLength (DNA char) =
        isMissing char ==> characterLength char === 0


transposeProperties :: TestTree
transposeProperties = testGroup "Transpose Properties"
    [ testProperty "Transpose = id --> Missing"  transposeIsIdIsMissing
    , testProperty "Length . Transpose = Length" transposeLengthInvariant
    , testProperty "Medians . Transpose = Medians" transposeMediansInvariant
    , testProperty "GappedMedians . Transpose = GappedMedians" transposeGappedMediansInvariant
    , testProperty "Left--Medians . Transpose = Right-Medians" transposeLeftMediansAreRightMedians
    , testProperty "Right-Medians . Transpose = Left--Medians" transposeRightMediansAreLeftMedians
    , testProperty "Alignment . Transpose = Alignment" transposeAlignIsAlign
    , testProperty "Deletion  . Transpose = Insertion" transposeDeleteIsInsert
    , testProperty "Insertion . Transpose = Deletion " transposeInsertIsDelete
    , testProperty "Gapped    . Transpose = Gapped   " transposeGappedIsGapped
    ]
  where
    transposeIsIdIsMissing :: DNA -> Property
    transposeIsIdIsMissing (DNA char) =
        transposeCharacter char == char ==> isMissing char

    transposeLengthInvariant :: DNA -> Property
    transposeLengthInvariant (DNA char) =
        (characterLength . transposeCharacter) char === characterLength char

    transposeMediansInvariant :: DNA -> Property
    transposeMediansInvariant (DNA char) =
        (extractMedians . transposeCharacter) char === extractMedians char

    transposeGappedMediansInvariant :: DNA -> Property
    transposeGappedMediansInvariant (DNA char) =
        (extractMediansGapped . transposeCharacter) char === extractMediansGapped char

    transposeLeftMediansAreRightMedians :: DNA -> Property
    transposeLeftMediansAreRightMedians (DNA char) =
        (extractMediansLeft . transposeCharacter) char === extractMediansRight char

    transposeRightMediansAreLeftMedians :: DNA -> Property
    transposeRightMediansAreLeftMedians (DNA char) =
        (extractMediansRight . transposeCharacter) char === extractMediansLeft char

    transposeAlignIsAlign :: DNA -> Property
    transposeAlignIsAlign (DNA char) =
        ((`isAlign` 0) . transposeCharacter) char === isAlign char 0

    transposeDeleteIsInsert :: DNA -> Property
    transposeDeleteIsInsert (DNA char) =
        ((`isDelete` 0) . transposeCharacter) char === isInsert char 0

    transposeInsertIsDelete :: DNA -> Property
    transposeInsertIsDelete (DNA char) =
        ((`isInsert` 0) . transposeCharacter) char === isDelete char 0

    transposeGappedIsGapped :: DNA -> Property
    transposeGappedIsGapped (DNA char) =
        ((`isGapped` 0) . transposeCharacter) char === isGapped char 0


gapSubstitutionSuite :: TestTree
gapSubstitutionSuite = testGroup "Gap Removal/Reinsertion Properties"
    [ testProperty "Medians . DeleteGaps = Medians" deleteGapsIsExtractMedians
    , testProperty "InsertGaps . DeleteGaps = Medians" insertGapsDeleteGapsLeftIdentity
    , testProperty "InsertGaps . DeleteGaps = Medians" insertGapsDeleteGapsRightIdentity
    ]
  where
    deleteGapsIsExtractMedians :: DNA -> Property
    deleteGapsIsExtractMedians (DNA char) =
        (extractMedians . snd . deleteGaps) char === extractMedians char



insertGapsDeleteGapsIdentity :: (GapSet -> DNA -> DNA) ->  DNA -> Property
insertGapsDeleteGapsIdentity applyGapsλ input@(DNA char) = counterexample shownMessage $
    extractMediansGapped char === extractMediansGapped char''
  where
    (gaps, char' ) = deleteGaps char
    (DNA   char'') = applyGapsλ gaps $ DNA char'
    shownMessage = unlines
        [ "Input:"
        , show input
        , "Output:"
        , show $ DNA char''
        , "Comparing gapped medians:"
        , show . Snip $ extractMediansGapped char
        , show . Snip $ extractMediansGapped char''
        ]


insertGapsDeleteGapsLeftIdentity :: DNA -> Property
insertGapsDeleteGapsLeftIdentity = insertGapsDeleteGapsIdentity f
  where
    f gaps (DNA char) =
        let gap   = getNucleotide nucleotideGap
            char' = ( extractMediansGapped char
                    , extractMediansGapped char
                    , SV.replicate (fromEnum $ characterLength char) 0
                    )
        in DNA $ insertGaps gap gaps nullGapSet char'

insertGapsDeleteGapsRightIdentity :: DNA -> Property
insertGapsDeleteGapsRightIdentity = insertGapsDeleteGapsIdentity f
  where
    f gaps (DNA char) =
        let gap   = getNucleotide nucleotideGap
            char' = ( SV.replicate (fromEnum $ characterLength char) 0
                    , extractMediansGapped char
                    , extractMediansGapped char
                    )
        in DNA $ insertGaps gap nullGapSet gaps char'


-- |
-- Test properties of pairwise string alignment algorithms.
validitySuite :: TestTree
validitySuite = testGroup "Validity of Alignments" $
    uncurry isValidPairwiseAlignment <$> alignmentMetricCombinations


isValidPairwiseAlignment
  :: String
  -> (SlimDynamicCharacter -> SlimDynamicCharacter -> (Word, SlimDynamicCharacter))
  -> TestTree
isValidPairwiseAlignment testLabel alignmentλ = testGroup ("Validity of " <> testLabel)
    [
       testPropertyRigorously "alignment function is commutative"         $ commutivity alignmentλ
     , testProperty           "output length is >= input length"            greaterThanOrEqualToInputLength
     , testProperty           "alignment length is =< sum of input lengths" alignmentLengthRelation
     , testPropertyLeniently  "alignment was not erroneously reversed"      isNotReversed
     , testPropertyLeniently  "alignment was not erroneously transposed"    isNotTransposed
     , testPropertyLeniently  "alignment without indels == input medians"   factorStatesEqualInputMedians
    ]
  where
    greaterThanOrEqualToInputLength :: DyadDNA -> Property
    greaterThanOrEqualToInputLength (DNA lhs :×: DNA rhs) =
        alnLen >= lhsLen .&&. alnLen >= rhsLen
      where
        (alnLen, lhsLen, rhsLen) = alignmentLengths lhs rhs

    alignmentLengthRelation :: DyadDNA -> Property
    alignmentLengthRelation (DNA lhs :×: DNA rhs) =
        counterexample shownCounterexample $ alnLen <= lhsLen + rhsLen
      where
        (alnLen, lhsLen, rhsLen) = alignmentLengths lhs rhs
        shownCounterexample = unwords [ show alnLen, ">", show lhsLen, "+", show rhsLen ]

    isNotReversed :: DyadDNA -> Property
    isNotReversed (DNA lhs :×: DNA rhs) =
        counterexample shownLHS (notPalindrome lhs' ==> notReversed lhs' lhs'') .&&.
        counterexample shownRHS (notPalindrome rhs' ==> notReversed rhs' rhs'')
      where
        (lhs' , rhs' ) = principalMedians lhs rhs
        (lhs'', rhs'') = alignmentFactors lhs rhs
        (_,   aligned) = alignmentλ lhs rhs

        shownLHS = fold [ "LHS:\n", show (DNA lhs), "\n", show (DNA aligned) ]
        shownRHS = fold [ "RHS:\n", show (DNA rhs), "\n", show (DNA aligned) ]

        notPalindrome x = notReversed x x
        notReversed x y = GV.reverse x /= y

    isNotTransposed :: DyadDNA -> Property
    isNotTransposed (DNA lhs :×: DNA rhs) =
        lhs' /= rhs' && lhs'' /= rhs'' ==>
            counterexample "lhs === rhs'" (lhs' /= rhs'') .&&.
            counterexample "rhs === lhs'" (rhs' /= lhs'')
      where
        (lhs' , rhs' ) = principalMedians lhs rhs
        (lhs'', rhs'') = alignmentFactors lhs rhs

    factorStatesEqualInputMedians :: DyadDNA -> Property
    factorStatesEqualInputMedians (DNA lhs :×: DNA rhs) =
        not (isMissing lhs) && not (isMissing rhs) ==>
            counterexample shownBoth (lhs' === lhs'' .||. rhs' === rhs'') .&&.
            counterexample shownLHS  (lhs' === lhs'') .&&.
            counterexample shownRHS  (rhs' === rhs'')
      where
        (lhs' , rhs' ) = principalMedians lhs rhs
        (lhs'', rhs'') = alignmentFactors lhs rhs
        (_,   aligned) = alignmentλ lhs rhs

        shownBoth = showAlignmentAndInputDNA "Both:" (DNA aligned) [DNA lhs, DNA rhs]
        shownLHS  = showAlignmentAndInputDNA "LHS:"  (DNA aligned) [DNA lhs]
        shownRHS  = showAlignmentAndInputDNA "RHS:"  (DNA aligned) [DNA rhs]

        showAlignmentAndInputDNA str alignment inDNAs =
            let sep = fold (replicate 30 "~=") <> "~\n"
            in  fold
                  [ str, "\n"
                  , sep
                  , intercalate sep $ show <$> inDNAs
                  , sep
                  , show alignment
                  , sep
                  ]

    -- Helpers

    alignmentLengths lhs rhs = ( alnLen, lhsLen, rhsLen )
      where
        alnLen = characterLength . snd $ alignmentλ lhs rhs
        lhsLen = characterLength lhs
        rhsLen = characterLength rhs

    alignmentFactors lhs rhs = (lhs', rhs')
      where
        (_, aligned) = alignmentλ lhs rhs
        lhs' = extractMediansLeft  aligned
        rhs' = extractMediansRight aligned

    principalMedians lhs rhs = (extractMediansGapped lhs, extractMediansGapped rhs)


commutivity
  :: (SlimDynamicCharacter -> SlimDynamicCharacter -> (Word, SlimDynamicCharacter))
  -> DyadDNA
  -> Property
commutivity alignmentλ (DNA lhs :×: DNA rhs) =
    counterexample shownCounterexample $
        x === y .||. (x === (c, transposeCharacter d) .&&. (a, transposeCharacter b) === y)
  where
    x@(a,b) = alignmentλ lhs rhs
    y@(c,d) = alignmentλ rhs lhs
    shownCounterexample = unlines
          [ "LHS ⊗ RHS:  ⇆  " <> show a
          , show (DNA b)
          , "RHS ⊗ LHS:  ⇆  " <> show c
          , show (DNA d)
          ]


-- |
-- Ensure that all three implementations return the same results
consistencySuite :: TestTree
consistencySuite = testGroup description $
    consistentResults <$> metricChoices
  where
    preamble = "All of these implementations return same alignments:"
    description = intercalate "\n        * " $ preamble : (fst <$> alignmentChoices)
      

consistentResults :: (String, MetricRepresentation Word32) -> TestTree
consistentResults param@(testLabel, _) =
    testPropertyRigorously ("Consistenty over " <> testLabel) $ checkAlignmentResults param


checkAlignmentResults :: (String, MetricRepresentation Word32) -> DyadDNA -> Property
checkAlignmentResults (testLabel, metric) (DNA lhs :×: DNA rhs) =
        allAreEqual alignmentResults
      where
        alignmentResults = getResults <$> alignmentChoices
        getResults = fmap (\f -> f metric lhs rhs)

        allAreEqual :: Foldable f => f (String, (Word, SlimDynamicCharacter)) -> Property
        allAreEqual input =
            case values of
              []   -> z
              x:xs -> counterexample errMsg $
                        foldr (\e a -> (e == x) .&&. a) z xs
          where
            z = property True
            p = intercalate "\n" . fmap ("    " <>) . lines
            renderToken (note, alignVal ) = note <> ":\n" <> renderAlign alignVal
            renderAlign (cost, alignment) = unlines [ p $ "Cost: " <> show cost, p . show $ DNA alignment ]
            values = snd <$> toList input
            errMsg = fromMaybe " " $
              case occurrences values of
                []  -> Nothing
                [_] -> Nothing
                (x,_):_ ->
                  let matchModeVal = (== x) . snd
                      (same, diff) = partition matchModeVal $ toList input
                  in  Just $ unlines
                        [ ""
                        , "The method is:"
                        , "  " <> testLabel
                        , "The mode is:"
                        , renderAlign x
                        , "Matched exactly by:"
                        , unlines $ p . fst <$> same
                        , "The following differed from the mode"
                        ] <> intercalate "\n" (p . renderToken <$> diff)


-- |
-- \( \mathcal{O} \left( n * \log_2 n \right) \)
--
-- Returns a mapping of each unique element in the list paired with how often
-- the element occurs in the list.
--
-- The elements are in descending order of occurrence.
--
-- ==_Example==
--
-- >>> occurrences "GATACACATCAGATT"
-- [('A',6),('T',4),('C',3),('G',2)]
--
-- >>> occurrences "AABCDDDEFGGT"
-- [('D',3),('A',2),('G',2),('B',1),('C',1),('E',1),('F',1),('T',1)]
occurrences :: (Foldable t, Ord a) => t a -> [(a, Word)]
occurrences = collateOccuranceMap . buildOccuranceMap
  where
    buildOccuranceMap = foldr occurrence mempty
      where
        occurrence e = insertWith (const succ) e 1
    collateOccuranceMap = sortBy comparator . assocs
      where
        comparator x y = descending $ comparing snd x y
        descending LT = GT
        descending GT = LT
        descending x  = x
