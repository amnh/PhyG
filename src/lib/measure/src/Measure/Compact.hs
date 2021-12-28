-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Compact
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE DeriveAnyClass     #-}
{-# LANGUAGE DeriveGeneric      #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE FlexibleContexts   #-}
{-# LANGUAGE Strict             #-}
{-# LANGUAGE UnboxedSums        #-}
{-# LANGUAGE LambdaCase         #-}
{-# LANGUAGE FlexibleInstances  #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Measure.Compact
  ( -- * Measure.Matrixs
    SCMλ
  , TCM2Dλ
  , TCM3Dλ
    -- * Measure Components
  , Distance
  , SymbolCount
  , SymbolIndex
    -- * Smart Constructors
  , CompactMeasure()
  , discreteCrossGap
  , discreteMetric
  , linearNorm
  , metricRepresentation
    -- * Accessor Type-classes
  , HasSymbolCount(..)
  , HasEditCosts(..)
  , HasSymbolChangeMatrix(..)
  , HasTransitionCostMatrix(..)
  , HasDenseMatrix(..)
  ) where

import Control.DeepSeq
import Data.Bits
import Data.Hashable
import Data.Hashable.Memoize
import Measure.Compact.Overlap
import qualified Measure.Compact.Discrete         as Dis
import qualified Measure.Compact.DiscreteCrossGap as Gap
import qualified Measure.Compact.L1Norm           as L1N
import Measure.Symbols.Dense (SCMρ)
import Measure.States.Dense (TCMρ, fromSCMρ)
import Measure.Matrix
import Measure.Unit
import GHC.Generics


-- |
-- Represents the metric for some discrete characters or dynamic characters.
-- The representation is highly optimized for both time and space efficiency.
-- Because of this, the term "compact" is used to describe the representation.
--
-- If any measure has an 'iota' symbol count, i.e. a symbol count of less than
-- or equal to 'infimumSymbolLimit', the compact representation will pre-compute
-- the entire transition cost matrix in a structure which is 'Storable' and
-- inter-operable with the C FFI. If and only if the measure is specified for an
-- 'iota' symbol count will 'getTCMρ' return a @Just@ value.
--
-- Additionally, the compact representation notes if the discrete metric, the
-- discrete metric adjoined by the gap symbol, or the L1 norm are the specified
-- measure. If any of these metrics are specified, specialized functions which
-- are more efficient will be returned when calling 'getSCMλ', 'getTCM2Dλ', and
-- 'getTCM3Dλ'.
--
-- Notably, if it is the case that /both/ the measure has an 'iota' symbol count
-- /and/ the measure is a specialized metric described above, then calling 'getTCMρ'
-- returns a /compile-time/ pre-computed transition cost matrix. Multiple
-- "constructions" of the same metric will result in the same "singlton" compact
-- representation.
--
-- Finally, if a specified measure has more than an 'iota' symbol count /and/ is
-- not one of the specialized metrics described above, then a /sparse/ and
-- /memoized/ representation is stored.
--
-- It is important to use this type in the metadata decorations rather than store
-- a function because a function cannot be contained in a compact region.
--
-- Use the 'getSCMλ', 'getTCM2Dλ', 'getTCM3Dλ', and 'getTCMρ' to the retrieve the
-- desired functions.
data  CompactMeasure a
    = ExplicitLayoutρ   {-# UNPACK #-} SCMρ        {-# UNPACK #-} TCMρ
    | ExplicitLayoutλ   {-# UNPACK #-} SCMρ        {-# UNPACK #-} Distance {-# UNPACK #-} Distance (TCM2Dλ a) (TCM3Dλ a)
    | DiscreteCrossGapρ {-# UNPACK #-} TCMρ
    | DiscreteCrossGapλ {-# UNPACK #-} SymbolCount {-# UNPACK #-} Distance {-# UNPACK #-} Distance
    | DiscreteMetricρ   {-# UNPACK #-} TCMρ
    | DiscreteMetricλ   {-# UNPACK #-} SymbolCount
    | LinearNormρ       {-# UNPACK #-} TCMρ
    | LinearNormλ       {-# UNPACK #-} SymbolCount
    deriving stock    (Generic)
    deriving anyclass (NFData)


instance Eq (CompactMeasure a) where

    (==) (ExplicitLayoutρ   scm _      ) (ExplicitLayoutρ   scm' _      ) = scm == scm'
    (==) (ExplicitLayoutλ   scm _ _ _ _) (ExplicitLayoutλ   scm' _ _ _ _) = scm == scm'
    (==) (DiscreteCrossGapρ tcm        ) (DiscreteCrossGapρ tcm'        ) = tcm == tcm'
    (==) (DiscreteCrossGapλ n g s      ) (DiscreteCrossGapλ n' g' s'    ) = n == n' && g == g' && s == s'
    (==) (DiscreteMetricρ   tcm        ) (DiscreteMetricρ   tcm'        ) = symbolCount tcm == symbolCount tcm'
    (==) (DiscreteMetricλ   n          ) (DiscreteMetricλ   n'          ) = n == n'
    (==) (LinearNormρ       tcm        ) (LinearNormρ       tcm'        ) = symbolCount tcm == symbolCount tcm'
    (==) (LinearNormλ       n          ) (LinearNormλ       n'          ) = n == n'
    (==) _ _                                                              = False


instance HasDenseMatrix (CompactMeasure a) where

    getTCMρ (ExplicitLayoutρ   _ tcm) = Just tcm
    getTCMρ  ExplicitLayoutλ   {}     = Nothing
    getTCMρ (DiscreteCrossGapρ tcm  ) = Just tcm
    getTCMρ  DiscreteCrossGapλ {}     = Nothing
    getTCMρ (DiscreteMetricρ   tcm  ) = Just tcm
    getTCMρ  DiscreteMetricλ   {}     = Nothing
    getTCMρ (LinearNormρ       tcm  ) = Just tcm
    getTCMρ  LinearNormλ       {}     = Nothing


instance HasEditCosts (CompactMeasure a) where

    {-# INLINEABLE maxEdit #-}
    maxEdit (ExplicitLayoutρ   _ tcm    ) = maxEdit tcm
    maxEdit (ExplicitLayoutλ   _ _ e _ _) = e
    maxEdit (DiscreteCrossGapρ tcm      ) = maxEdit tcm
    maxEdit (DiscreteCrossGapλ _ g _    ) = g
    maxEdit (DiscreteMetricρ   tcm      ) = maxEdit tcm
    maxEdit  DiscreteMetricλ   {}         = 1
    maxEdit (LinearNormρ       tcm      ) = maxEdit tcm
    maxEdit  LinearNormλ       {}         = 1

    maxDeletion (ExplicitLayoutρ   _ tcm      ) = maxDeletion tcm
    maxDeletion (ExplicitLayoutλ   scm _ _ _ _) = maxDeletion scm
    maxDeletion (DiscreteCrossGapρ tcm        ) = maxDeletion tcm
    maxDeletion (DiscreteCrossGapλ _ g _      ) = g
    maxDeletion (DiscreteMetricρ   tcm        ) = maxDeletion tcm
    maxDeletion  DiscreteMetricλ   {}           = 1
    maxDeletion (LinearNormρ       tcm        ) = maxDeletion tcm
    maxDeletion  LinearNormλ       {}           = 1

    maxInsertion (ExplicitLayoutρ   _ tcm      ) = maxInsertion tcm
    maxInsertion (ExplicitLayoutλ   scm _ _ _ _) = maxInsertion scm
    maxInsertion (DiscreteCrossGapρ tcm        ) = maxInsertion tcm
    maxInsertion (DiscreteCrossGapλ _ g _      ) = g
    maxInsertion (DiscreteMetricρ   tcm        ) = maxInsertion tcm
    maxInsertion  DiscreteMetricλ   {}           = 1
    maxInsertion (LinearNormρ       tcm        ) = maxInsertion tcm
    maxInsertion  LinearNormλ       {}           = 1

    {-# INLINEABLE minEdit #-}
    minEdit (ExplicitLayoutρ   _ tcm    ) = minEdit tcm
    minEdit (ExplicitLayoutλ   _ e _ _ _) = e
    minEdit (DiscreteCrossGapρ tcm      ) = minEdit tcm
    minEdit (DiscreteCrossGapλ _ g _    ) = g
    minEdit (DiscreteMetricρ   tcm      ) = minEdit tcm
    minEdit  DiscreteMetricλ   {}         = 1
    minEdit (LinearNormρ       tcm      ) = minEdit tcm
    minEdit  LinearNormλ       {}         = 1

    minDeletion (ExplicitLayoutρ   _ tcm      ) = minDeletion tcm
    minDeletion (ExplicitLayoutλ   scm _ _ _ _) = minDeletion scm
    minDeletion (DiscreteCrossGapρ tcm        ) = minDeletion tcm
    minDeletion (DiscreteCrossGapλ _ g _      ) = g
    minDeletion (DiscreteMetricρ   tcm        ) = minDeletion tcm
    minDeletion  DiscreteMetricλ   {}           = 1
    minDeletion (LinearNormρ       tcm        ) = minDeletion tcm
    minDeletion  LinearNormλ       {}           = 1

    minInsertion (ExplicitLayoutρ   _ tcm      ) = minInsertion tcm
    minInsertion (ExplicitLayoutλ   scm _ _ _ _) = minInsertion scm
    minInsertion (DiscreteCrossGapρ tcm        ) = minInsertion tcm
    minInsertion (DiscreteCrossGapλ _ g _      ) = g
    minInsertion (DiscreteMetricρ   tcm        ) = minInsertion tcm
    minInsertion  DiscreteMetricλ   {}           = 1
    minInsertion (LinearNormρ       tcm        ) = minInsertion tcm
    minInsertion  LinearNormλ       {}           = 1


instance HasSymbolChangeMatrix (CompactMeasure a) where

    getSCMλ (ExplicitLayoutρ   _ tcm      ) = getSCMλ tcm
    getSCMλ (ExplicitLayoutλ   scm _ _ _ _) = getSCMλ scm
    getSCMλ (DiscreteCrossGapρ tcm        ) = getSCMλ tcm
    getSCMλ (DiscreteCrossGapλ _ g s      ) = Gap.scmλ g s
    getSCMλ (DiscreteMetricρ   tcm        ) = getSCMλ tcm
    getSCMλ  DiscreteMetricλ   {}           = Dis.scmλ
    getSCMλ (LinearNormρ       tcm        ) = getSCMλ tcm
    getSCMλ  LinearNormλ       {}           = L1N.scmλ


instance HasSymbolCount (CompactMeasure a) where

    symbolCount (ExplicitLayoutρ   _ tcm      ) = symbolCount tcm
    symbolCount (ExplicitLayoutλ   scm _ _ _ _) = symbolCount scm
    symbolCount (DiscreteCrossGapρ tcm        ) = symbolCount tcm
    symbolCount (DiscreteCrossGapλ n _ _      ) = n
    symbolCount (DiscreteMetricρ   tcm        ) = symbolCount tcm
    symbolCount (DiscreteMetricλ   n          ) = n
    symbolCount (LinearNormρ       tcm        ) = symbolCount tcm
    symbolCount (LinearNormλ       n          ) = n


instance (Enum a, FiniteBits a) => HasTransitionCostMatrix (CompactMeasure a) a where

    getTCM2Dλ (ExplicitLayoutρ   _ tcm       ) = getTCM2Dλ tcm
    getTCM2Dλ (ExplicitLayoutλ   _ _ _ tcmλ _) = tcmλ
    getTCM2Dλ (DiscreteCrossGapρ tcm         ) = getTCM2Dλ tcm
    getTCM2Dλ (DiscreteCrossGapλ n g s       ) = Gap.tcm2Dλ n g s
    getTCM2Dλ (DiscreteMetricρ   tcm         ) = getTCM2Dλ tcm
    getTCM2Dλ  DiscreteMetricλ   {}            = Dis.tcm2Dλ
    getTCM2Dλ (LinearNormρ       tcm         ) = getTCM2Dλ tcm
    getTCM2Dλ (LinearNormλ       n           ) = L1N.tcm2Dλ n

    getTCM3Dλ (ExplicitLayoutρ   _ tcm       ) = getTCM3Dλ tcm
    getTCM3Dλ (ExplicitLayoutλ   _ _ _ _ tcmλ) = tcmλ
    getTCM3Dλ (DiscreteCrossGapρ tcm         ) = getTCM3Dλ tcm
    getTCM3Dλ (DiscreteCrossGapλ n g s       ) = Gap.tcm3Dλ n g s
    getTCM3Dλ (DiscreteMetricρ   tcm         ) = getTCM3Dλ tcm
    getTCM3Dλ  DiscreteMetricλ   {}            = Dis.tcm3Dλ
    getTCM3Dλ (LinearNormρ       tcm         ) = getTCM3Dλ tcm
    getTCM3Dλ (LinearNormλ       n           ) = L1N.tcm3Dλ n


instance Show (CompactMeasure a) where

    show ExplicitLayoutρ   {} = "General Metric"
    show ExplicitLayoutλ   {} = "General Metric"
    show DiscreteCrossGapρ {} = "Discrete Metric ⨯ Gap Symbol"
    show DiscreteCrossGapλ {} = "Discrete Metric ⨯ Gap Symbol"
    show DiscreteMetricρ   {} = "Discrete Metric"
    show DiscreteMetricλ   {} = "Discrete Metric"
    show LinearNormρ       {} = "1st Linear Norm"
    show LinearNormλ       {} = "1st Linear Norm"


-- |
-- Nullary constructor for the <https://en.wikipedia.org/wiki/Discrete_space discrete metric>.
discreteMetric :: SymbolCount -> CompactMeasure a
discreteMetric =
    \case
      2 -> DiscreteMetricρ Dis.tcmρ2
      3 -> DiscreteMetricρ Dis.tcmρ3
      4 -> DiscreteMetricρ Dis.tcmρ4
      5 -> DiscreteMetricρ Dis.tcmρ5
      6 -> DiscreteMetricρ Dis.tcmρ6
      7 -> DiscreteMetricρ Dis.tcmρ7
      8 -> DiscreteMetricρ Dis.tcmρ8
      n -> DiscreteMetricλ n


-- |
-- Nullary constructor for the <https://en.wikipedia.org/wiki/Discrete_space discrete metric>.
discreteCrossGap :: Distance -> Distance -> SymbolCount -> CompactMeasure a
discreteCrossGap 1 2 =
    \case
      2 -> DiscreteCrossGapρ Gap.tcm12ρ2
      3 -> DiscreteCrossGapρ Gap.tcm12ρ3
      4 -> DiscreteCrossGapρ Gap.tcm12ρ4
      5 -> DiscreteCrossGapρ Gap.tcm12ρ5
      6 -> DiscreteCrossGapρ Gap.tcm12ρ6
      7 -> DiscreteCrossGapρ Gap.tcm12ρ7
      8 -> DiscreteCrossGapρ Gap.tcm12ρ8
      n -> DiscreteCrossGapλ n 1 2
discreteCrossGap gap sub =
    \case
      n | iota n -> DiscreteCrossGapρ $ Gap.tcmρ n gap sub
      n          -> DiscreteCrossGapλ n gap sub


-- |
-- Nullary constructor for the <https://en.wikipedia.org/wiki/Lp_space 1st linear norm>.
linearNorm :: SymbolCount -> CompactMeasure a
linearNorm =
    \case
      2 -> LinearNormρ L1N.tcmρ2
      3 -> LinearNormρ L1N.tcmρ3
      4 -> LinearNormρ L1N.tcmρ4
      5 -> LinearNormρ L1N.tcmρ5
      6 -> LinearNormρ L1N.tcmρ6
      7 -> LinearNormρ L1N.tcmρ7
      8 -> LinearNormρ L1N.tcmρ8
      n -> LinearNormλ n


-- |
-- General constructor for an arbitrary metric.
--
-- Performs memoization so repeated value queries are not recomputed.
metricRepresentation
  :: ( FiniteBits a
     , Hashable a
     , NFData a
     )
  => SCMρ
  -> CompactMeasure a
metricRepresentation scmρ
  | iota scmρ = ExplicitLayoutρ scmρ $ fromSCMρ 0 scmρ
  | otherwise = 
    let scmλ  = getSCMλ scmρ
        minΔ  = minEdit scmρ
        maxΔ  = maxEdit scmρ
        size  = symbolCount scmρ
        tcm2D = memoize2 $ overlap2 size scmλ
        tcm3D = memoize3 $ overlap3 size scmλ
    in  ExplicitLayoutλ scmρ minΔ maxΔ tcm2D tcm3D

