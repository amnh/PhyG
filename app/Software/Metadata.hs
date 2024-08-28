{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE Strict #-}
{-# LANGUAGE TemplateHaskell #-}

{- |
Important project metadata which is programmatically accessible.

All temporally variable information is generated at compile-time, reflecting the
state of the project when the program was built.
-}
module Software.Metadata (
    -- * Data-type
    AbreviatedName (),

    -- ** Accessors
    nameAbbreviation,
    nameExpansion,

    -- * Metadata Names
    projectName,
    softwareName,

    -- * Versioning
    shortVersionInformation,
    fullVersionInformation,

    -- * Time of build
    timeOfCompilation,
) where

import Data.Foldable (toList)
import Data.Functor ((<&>))
import Data.String
import Data.Version (showVersion)
import Development.GitRev (gitCommitCount, gitHash)
import PackageInfo_PhyG (name, version)
import Software.Metadata.TimeStamp (compilationTimeStamp, renderTimeStampAsLocalTime)


{- |
A name which also is an acronym.
-}
data AbreviatedName = AbbreviationOf
    { getNameAbreviated ∷ String
    , getNameExpansion ∷ String
    }


{- |
Get the abbreviation or acronym of an 'AbreviatedName'.
-}
nameAbbreviation ∷ (IsString s) ⇒ AbreviatedName → s
nameAbbreviation = fromString . getNameAbreviated


{- |
Get the full expansion on an 'AbreviatedName'.
-}
nameExpansion ∷ (IsString s) ⇒ AbreviatedName → s
nameExpansion = fromString . getNameExpansion


{- |
Name of the larger software project.
-}
projectName ∷ AbreviatedName
projectName = "PHANE" `AbbreviationOf` "Phylogenetic Haskell Analytic Network Engine"


{- |
Name of the larger software project.
-}
softwareName ∷ AbreviatedName
softwareName = name `AbbreviationOf` "Phylogenetic Graph"


{- |
Brief description of the software version.
-}
shortVersionInformation ∷ (IsString s, Semigroup s) ⇒ s
--shortVersionInformation = "β-version " <> fromString (showVersion version)
shortVersionInformation = "Version " <> fromString (showVersion version)


{- |
Full description of the software version.

Uses @TemplateHaskell@ to splice in git hash and commit count information
from the compilation environment.
-}
fullVersionInformation ∷ (IsString s, Monoid s) ⇒ IO s
fullVersionInformation =
    let intercalate' ∷ (Foldable f, Monoid s) ⇒ s → f s → s
        intercalate' sep =
            let consider [] = mempty
                consider [x] = x
                consider (x : xs) = x <> foldMap (sep <>) xs
            in  consider . toList

        commits = "(" <> fromString $(gitCommitCount) <> " commits)"
        hashInfo = "[" <> fromString (take 7 $(gitHash)) <> "]"
        lessName = "(" <> nameAbbreviation softwareName <> ")"
        longName = nameExpansion softwareName
    in  timeOfCompilation <&> \time →
            intercalate'
                " "
                [ longName
                , lessName
                , shortVersionInformation
                , hashInfo
                , commits
                , "on"
                , time
                ]


{- |
The UTC system time at which (this module of) the binary was compiled.
-}
timeOfCompilation ∷ (IsString s) ⇒ IO s
timeOfCompilation = renderTimeStampAsLocalTime $$(compilationTimeStamp)
