-----------------------------------------------------------------------------
-- |
-- Module      :  Test.Aligners
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

module Test.Aligners
  ( alignmentMetricCombinations
  , alignmentChoices
  , metricChoices
  ) where

import Bio.DynamicCharacter
import Data.Bifunctor
import Data.Bits
import Data.Coerce
import Data.Word
import DirectOptimization.Pairwise
import DirectOptimization.Pairwise.Swapping
import DirectOptimization.Pairwise.Ukkonen
import Foreign.C.Types                            (CUInt(..))
import Measure.Transition.Diagnosis
import Measure.Transition.Representation
import Layout.Compact.States
import Layout.Compact.Symbols                 (fromList)
import Measure.Unit.AlignmentCost
import Measure.Unit.SymbolCount
import Test.QuickCheck.Instances.DynamicCharacter


alignmentMetricCombinations
  :: [ (String, SlimDynamicCharacter -> SlimDynamicCharacter -> (AlignmentCost, SlimDynamicCharacter)) ]
alignmentMetricCombinations = do
    -- Generate all combinations
    (mLabel, metric ) <- metricChoices
    (aLabel, mkAlign) <- alignmentChoices
    let niceLabel = unwords [ mLabel, "＆", aLabel ]
    pure (niceLabel, mkAlign metric)


alignmentChoices
  :: [ ( String
       ,    TransitionMatrix Word32
         -> SlimDynamicCharacter
         -> SlimDynamicCharacter
         -> (AlignmentCost, SlimDynamicCharacter)
       )
     ]
alignmentChoices =
    [ ("FFI Code", slimPairwiseDO . metricRepresentationToCompactTCM )
    , ("Swapping", swappingDO . translateSlimStateTCM )
    , ("Ukkonens", \m -> ukkonenDO (maxEdit m) (translateSlimStateTCM m) )
    ]


metricRepresentationToCompactTCM :: TransitionMatrix a -> TCMρ
metricRepresentationToCompactTCM =
    let dim = SymbolCount . toEnum $ length nucleotideAlphabet
    in  fromSDMλ dim 0 . symbolDistances


translateSlimStateTCM :: TransitionMatrix Word32 -> TCM2Dλ SlimState
translateSlimStateTCM m =
    let tcm = stateTransitionPairwiseDispersion m
        c2w = coerce :: SlimState -> Word32
        w2c = coerce :: Word32 -> SlimState
    in  \x y -> second w2c $ tcm (c2w x) (c2w y)


metricChoices :: [(String, TransitionMatrix Word32)]
metricChoices =
    [ ("Discrete Metric", discreteMetric len)
    , ("1st Linear Norm", linearNorm     len)
    , ("Sub-InDel (1:2)", subInDel 1 2      )
    , ("Sub-InDel (2:1)", subInDel 2 1      )
    ]
  where
    len = SymbolCount . toEnum $ finiteBitSize nucleotideGap

    subInDel :: Word -> Word -> TransitionMatrix Word32
    subInDel x g = metricRepresentation . factoredSDM . fromList $
        let indices = [ 0 .. length nucleotideAlphabet - 1 ]
        in  [ if i == j then 0 else if i == 0 || j == 0 then g else x | i <- indices, j <- indices ]
