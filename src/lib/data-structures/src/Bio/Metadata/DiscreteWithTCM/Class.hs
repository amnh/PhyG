-----------------------------------------------------------------------------
-- |
-- Module      :  Bio.Metadata.DiscreteWithTCM.Class
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE FlexibleContexts       #-}
{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE MultiParamTypeClasses  #-}

module Bio.Metadata.DiscreteWithTCM.Class
  ( DiscreteWithTcmCharacterMetadata()
  , GetSymbolChangeMatrix(..)
  , GetPairwiseTransitionCostMatrix(..)
  , GetSparseTransitionCostMatrix(..)
  , HasCharacterAlphabet(..)
  , HasCharacterName(..)
  , HasCharacterWeight(..)
  ) where

import Bio.Character.Encodable
import Bio.Metadata.Discrete
import Bio.Metadata.Metric


-- |
-- A decoration of the encoding of the metadata for a metric or dynamic character
-- which has the appropriate 'Control.Lens.Type.Lens' and class constraints.
class ( DiscreteCharacterMetadata     s
      , EncodableStreamElement        c
      , GetSymbolChangeMatrix         s (Word -> Word -> Word)
      , GetPairwiseTransitionCostMatrix       s c Word
      ) => DiscreteWithTcmCharacterMetadata s c | s -> c where
