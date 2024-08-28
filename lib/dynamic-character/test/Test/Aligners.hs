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
import Data.MetricRepresentation
import Data.TCM                                   (fromList)
import Data.TCM.Dense
import Data.Word
import Foreign.C.Types (CUInt)
import DirectOptimization.Pairwise
import DirectOptimization.Pairwise.Swapping
import DirectOptimization.Pairwise.Ukkonen
import Test.QuickCheck.Instances.DynamicCharacter


alignmentMetricCombinations
  :: [ (String, SlimDynamicCharacter -> SlimDynamicCharacter -> (Word, SlimDynamicCharacter)) ]
alignmentMetricCombinations = do
    -- Generate all combinations
    (mLabel, metric ) <- metricChoices
    (aLabel, mkAlign) <- alignmentChoices
    let niceLabel = unwords [ mLabel, "＆", aLabel ]
    pure (niceLabel, mkAlign metric)


alignmentChoices
  :: [ ( String
       ,    MetricRepresentation Word32
         -> SlimDynamicCharacter
         -> SlimDynamicCharacter
         -> (Word, SlimDynamicCharacter)
       )
     ]
alignmentChoices =
      [ ("FFI Code", slimPairwiseDO . metricRepresentationToDenseTCM )
      , ("Swapping", swappingDO . translateSlimStateTCM )
      , ("Ukkonens", \m -> ukkonenDO (maxInDelCost m) (translateSlimStateTCM m) )
      ]


metricRepresentationToDenseTCM :: MetricRepresentation a -> DenseTransitionCostMatrix
metricRepresentationToDenseTCM =
      generateDenseTransitionCostMatrix 0 (toEnum (length nucleotideAlphabet)) . retreiveSCM


translateSlimStateTCM :: MetricRepresentation Word32 -> SlimState -> SlimState -> (SlimState, Word)
translateSlimStateTCM m =
    let tcm = retreivePairwiseTCM m
        c2w = (toEnum . fromEnum) :: SlimState -> Word32
        w2c = (toEnum . fromEnum) :: Word32 -> SlimState
    in  \x y -> first w2c $ tcm (c2w x) (c2w y)


metricChoices :: [(String, MetricRepresentation Word32)]
metricChoices =
      [ ("Discrete Metric", discreteMetric)
      , ("1st Linear Norm", linearNorm len)
{-
      , ("Sub-InDel (1:2)", subInDel 1 2  )
      , ("Sub-InDel (2:1)", subInDel 2 1  )
-}
      ]
  where
    len = toEnum $ finiteBitSize nucleotideGap

-- Comment out the monadic code
{-
    subInDel :: Word -> Word -> MetricRepresentation Word32
    subInDel x g = metricRepresentation . snd . fromList $
        let indices = [ 0 .. length nucleotideAlphabet - 1 ]
        in  [ if i == j then 0 else if i == 0 || j == 0 then g else x | i <- indices, j <- indices ]
-}
