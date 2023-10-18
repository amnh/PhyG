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

import Data.Foldable (fold)
import Data.List.NonEmpty (NonEmpty(..))
import Data.String (IsString(fromString))
import Data.Text.Builder.Linear
import Data.Text (unpack)
import Data.Text qualified as T
import PackageInfo_PhyG (copyright, homepage, synopsis)
import Software.License
import Software.Metadata
import System.Environment (getProgName)


{- |
Text to be displayed at the beginning of each computation.
-}
preambleText :: IO Builder
preambleText =
    let copywriting = fromString copyright
        description = fromString synopsis
        homepageURL = fromString homepage
        licenseName = getLicenseLine 1

        above = encloseAbove :| []
        below = encloseBelow :| []
        inner progName =
            let makeMainBlock = flip (:|)
                    [ ""
                    , copywriting
                    , "The program '" <> fromString progName <> "' comes with ABSOLUTELY NO WARRANTY;"
                    , "This is free software, and may be redistributed under the " <> licenseName <> "."
                    , ""
                    , description
                    , "For more information, see: " <> homepageURL
                    ]
            in  fmap encloseLine . makeMainBlock <$> fullVersionInformation

        surround x = above <> x <> below

        topProgBlock = fmap surround . inner

    in  getProgName >>= fmap (intercalate' "\n") . topProgBlock


enclosing :: Char -> String -> String -> String -> Builder
enclosing pad strL strR strM =
    let widthPreamble :: Int
        widthPreamble = 100
        widthL = length strL
        widthR = length strR
        widthM = length strM
        padded = replicate (widthPreamble - widthL - widthM - widthR) pad
    in  fromString $ fold [ strL, strM, padded, strR ]


encloseAbove :: Builder
encloseAbove = enclosing '─' "┌" "┐" ""


encloseBelow :: Builder
encloseBelow = enclosing '─' "└" "┘" ""


encloseLine :: String -> Builder
encloseLine = enclosing ' ' "│ " " │"


getLicenseLine :: Word -> String
getLicenseLine lineNum =
    let n = fromEnum lineNum
    in  case drop (n - 1) . take n $ T.lines licenseText of
          [] -> mempty
          x:_ -> unpack x


intercalate' :: Builder -> NonEmpty Builder -> Builder
intercalate' sep = \case
    x :| [] -> x
    x :| xs -> x <> foldMap (sep <>) xs
