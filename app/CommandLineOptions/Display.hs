{- |
Human readable renderings for textual output of important program information.
-}

{-# Language LambdaCase #-}
{-# Language OverloadedStrings #-}
{-# Language Strict #-}
{-# Language TemplateHaskell #-}

module CommandLineOptions.Display
  ( printInformationDisplay
  ) where

import Data.Foldable
import Data.List.NonEmpty (NonEmpty(..), sort)
import Data.Text.Builder.Linear
import Data.Text.IO
import CommandLineOptions.Types
import Software.Credits
import Software.License
import Software.Metadata
import Software.SplashImage
import Prelude hiding (putStrLn)


{- |
Gather the information blocks specified from the CLI options and nicely render the blocks together.
-}
printInformationDisplay :: NonEmpty DisplayInfoBlock -> IO ()
printInformationDisplay =
    let bordering :: Builder
        bordering = "\n"

        delimiter :: Builder
        delimiter = "\n\n\n"

        encloseNonSingle :: Builder -> (NonEmpty Builder -> Builder) -> NonEmpty Builder -> Builder
        encloseNonSingle endcap f blocks@(_ :| bs) =
            let joined = f blocks
            in  case bs of
                  [] -> joined
                  _ -> endcap <> joined <> endcap

        getBuilder :: DisplayInfoBlock -> IO Builder
        getBuilder = \case
            DisplayVersion -> builderVersionInfo
            DisplayLicense -> pure builderLicenseText
            DisplaySplash  -> pure builderSplashImage
            DisplayCredits -> pure builderCreditsRoll

        joinBuilders :: NonEmpty Builder -> Builder
        joinBuilders = encloseNonSingle bordering (intercalate' delimiter)

        printBuilder :: IO Builder -> IO ()
        printBuilder b = fmap runBuilder b >>= putStrLn

    in  printBuilder . fmap joinBuilders . traverse getBuilder . sort


builderLicenseText :: Builder
builderLicenseText = fromText licenseText


builderCreditsRoll :: Builder
builderCreditsRoll = contributors


builderVersionInfo :: IO Builder
builderVersionInfo = fullVersionInformation


builderSplashImage :: Builder
builderSplashImage = fromText $$(splashImage)


intercalate' :: Builder -> NonEmpty Builder -> Builder
intercalate' sep = \case
    x :| [] -> x
    x :| xs -> x <> foldMap (sep <>) xs
