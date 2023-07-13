{- |
Compile-time embeddings of project contributions via Template Haskell.
-}

{-# Language OverloadedStrings #-}
{-# Language Strict #-}

module Software.License
  ( licenseText
  ) where

import Data.Text (intercalate, lines)
import Data.Text.IO (readFile)
import Instances.TH.Lift ()
import Language.Haskell.TH hiding (Inline)
import Language.Haskell.TH.Syntax hiding (Inline)
import Paths_PhyG (getDataFileName)
import Prelude hiding (lines, readFile)
import System.FilePath (normalise)


{- |
The text of the software license under which the PHANE tool is distributed.
-}
licenseText :: ExpQ
licenseText =
    let gatherData = getDataFileName "doc/LICENSE" >>= polishData . readFile . normalise
        polishData = fmap (intercalate "\n" . lines)
    in  runIO gatherData >>= lift
