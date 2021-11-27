-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Diagnosis
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE DeriveAnyClass     #-}
{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE DeriveGeneric      #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE OverloadedLists    #-}
{-# LANGUAGE Strict             #-}
{-# LANGUAGE TypeFamilies       #-}

module Measure.Diagnosis
  ( -- * Diagnosis of structure
    DiagnosisOfSCM(..)
    -- * Diagnosing functions
  , diagnoseSCMρ
  , diagnoseSCMλ
  ) where

import           Control.DeepSeq
import           Data.Data
import qualified Data.Vector.Storable     as V
import           GHC.Generics
import           Measure.Metricity
import           Measure.Symbols.Internal
import           Measure.Matrix
import           Measure.Unit.SymbolCount


-- |
-- The result of a call to 'diagnoseSCM'.
data  DiagnosisOfSCM scm
    = DiagnosisOfSCM
    { metricity   :: !Metricity
      -- ^ The most restrictive present in the 'factoredSCM'.
    , coefficient :: {-# UNPACK #-} !Rational
      -- ^ The multiplicative constant factor of a symbol change matrix. Minimum
      -- value of the multiplicative identity /one/.
    , factoredSCM :: !scm
      -- ^ The new symbol change matrix with each value divided by the 'coefficient'.
    }
    deriving stock    (Data, Eq, Generic, Show, Typeable)
    deriving anyclass (NFData)


diagnoseSCMρ :: SCMρ -> DiagnosisOfSCM SCMρ
diagnoseSCMρ scmρ@(SCMρ n vec) =
    let value = V.foldl1' gcd vec
        scmρ' | value <= 1 = scmρ
              | otherwise  = SCMρ n $ V.map (`div` value) vec
    in  DiagnosisOfSCM
        { metricity   = metricityOfSCMρ scmρ
        , coefficient = fromIntegral value
        , factoredSCM = scmρ'
        }


diagnoseSCMλ :: SymbolCount -> SCMλ -> DiagnosisOfSCM SCMλ
diagnoseSCMλ size scmλ =
    let dim   = fromIntegral $ size * size
        idx i = fromIntegral . uncurry scmλ $ fromIntegral i `quotRem` fromIntegral size
        pairs = V.generate dim idx :: V.Vector Word
        value = fromIntegral $ V.foldl1' gcd pairs
        scmλ' | value <= 1 = scmλ
              | otherwise  = \i j -> scmλ i j `div` value
    in  DiagnosisOfSCM
        { metricity   = metricityOfSCMλ size scmλ'
        , coefficient = fromIntegral value
        , factoredSCM = scmλ'
        }
