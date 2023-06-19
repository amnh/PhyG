------------------------------------------------------------------------------
-- |
-- Module      :  Software.Metadata
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# Language OverloadedStrings #-}
{-# Language TemplateHaskell   #-}

module Software.Metadata
  ( softwareName
  , shortVersionInformation
  , fullVersionInformation
  , timeOfCompilation
  ) where

import Data.Foldable (fold)
import Data.String
import Data.Version (showVersion)
import Development.GitRev (gitCommitCount, gitHash)
import Paths_PhyGraph (version)
import Software.Metadata.TimeStamp (compilationTimeStamp)


-- |
-- Name of the software package.
softwareName :: IsString s => s
softwareName = "Phylogenetic Graph"


-- |
-- Brief description of the software version.
shortVersionInformation :: (IsString s, Semigroup s) => s
shortVersionInformation = "(PhyG) β-version " <> fromString (showVersion version)


-- |
-- Full escription of the software version.
--
-- Uses @TemplateHaskell@ to splice in git hash and commit count information
-- from the compilation environment.
fullVersionInformation :: (IsString s, Monoid s) => s
fullVersionInformation = fold
    [ softwareName
    , " "
    , shortVersionInformation
    , " ["
    , fromString $ take 7 $(gitHash)
    , "] ("
    , fromString $(gitCommitCount)
    , " commits)"
    ]


{- |
The UTC system time at which (this module of) the binary was compiled.
-}
timeOfCompilation :: String
timeOfCompilation = $$(compilationTimeStamp)
