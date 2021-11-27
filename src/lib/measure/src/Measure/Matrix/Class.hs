-----------------------------------------------------------------------------
-- |
-- Module      :  Measure.Matrix.Class
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE FlexibleContexts           #-}
{-# LANGUAGE FlexibleInstances          #-}
{-# LANGUAGE MultiParamTypeClasses      #-}
{-# LANGUAGE Strict                     #-}
{-# LANGUAGE TypeSynonymInstances       #-}

module Measure.Matrix.Class
  ( HasDenseMatrix(..)
  , HasSymbolChangeMatrix(..)
  , HasTransitionCostMatrix(..)
  ) where

import Data.Bits
import Foreign.C.Types (CUInt)
import Measure.Unit.SymbolIndex
import Measure.Matrix.Types
import Measure.States.Dense.Structure

-- |
-- Any structural representation which can produce a Transition Cost Matrix.
class HasDenseMatrix a where

    getTCMρ :: a -> Maybe TCMρ 


-- |
-- Any structural representation which can produce a Symbol Change Matrix.
class HasSymbolChangeMatrix a where

    getSCMλ :: a -> SCMλ


-- |
-- Any structural representation which can produce a Transition Cost Matrix.
class HasTransitionCostMatrix a e where

    getTCM2Dλ :: a -> TCM2Dλ e

    getTCM3Dλ :: a -> TCM3Dλ e


instance HasSymbolChangeMatrix SCMλ where

    {-# INLINE getSCMλ #-}
    getSCMλ = id


instance HasSymbolChangeMatrix TCMρ where

    getSCMλ tcmρ i =
        let f = bit . fromEnum :: SymbolIndex -> CUInt
        in  fst . getTCM2Dλ tcmρ (f i) . f


instance Enum a => HasTransitionCostMatrix TCMρ a where

    -- |
    -- /O(1)/
    {-# SPECIALISE getTCM2Dλ :: TCMρ -> TCM2Dλ CUInt #-}
    getTCM2Dλ tcmρ   = getTCMρ2D tcmρ

    -- |
    -- /O(1)/
    {-# SPECIALISE getTCM3Dλ :: TCMρ -> TCM3Dλ CUInt #-}
    getTCM3Dλ tcmρ i = getTCMρ3D tcmρ i


instance HasDenseMatrix TCMρ where

    getTCMρ = Just

