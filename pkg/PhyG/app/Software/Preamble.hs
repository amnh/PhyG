{- |
Human readable message describing the program being used.
-}

{-# Language LambdaCase #-}
{-# Language OverloadedStrings #-}
{-# Language Strict #-}
{-# Language TemplateHaskell #-}

module Software.Preamble
  ( preambleText
  ) where

import Data.Functor ((<&>))
import Data.List.NonEmpty (NonEmpty(..))
import Data.String (IsString(fromString))
import Data.Text.Builder.Linear
import Data.Text qualified as T
import Software.License
import Software.Metadata
import System.Environment (getProgName)


preambleText :: IO Builder
preambleText =
    let copywriting = getLicenseLine 3
        licenseName = getLicenseLine 1
    in  getProgName <&> \progName -> intercalate' "\n" $
            fullVersionInformation :|
            [ ""
            , copywriting
            , fromString progName <> " comes with ABSOLUTELY NO WARRANTY;"
            , "This is free software, and may be redistributed under the " <> licenseName
            ]


getLicenseLine :: Word -> Builder
getLicenseLine lineNum =
    let fullLicenseText = $(licenseText)
        n = fromEnum lineNum
    in  fromText . head . drop (n - 1) . take n $ T.lines fullLicenseText


intercalate' :: Builder -> NonEmpty Builder -> Builder
intercalate' sep = \case
    x :| [] -> x
    x :| xs -> x <> foldMap (sep <>) xs
