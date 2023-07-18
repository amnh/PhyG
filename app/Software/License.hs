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
--import Paths_PhyG (getDataFileName)
import Paths_PhyG
import Prelude hiding (lines, readFile)
import System.FilePath (normalise)
import System.Directory


{- |
The text of the software license under which the PHANE tool is distributed.
-}
licenseText :: ExpQ
licenseText =
    let filePath = "doc/LICENSE"
        gatherData = polishData . readFile . normalise
        polishData = fmap (intercalate "\n" . lines)
        inlineData = gatherData filePath >>= lift
        noting key = putStrLn . ("\n>>> [TEMPLATE HASKELL NOTE]:\n\t" <>) . (key <>) . (":\t" <>)
        workingIO = do
            getDataFileName "LICENSE" >>= noting "getDataFileName"
            getBinDir     >>= noting "getBinDir"
            getLibDir     >>= noting "getLibDir"
            getDynLibDir  >>= noting "getDynLibDir"
            getDataDir    >>= noting "getDataDir"
            getLibexecDir >>= noting "getLibexecDir"
            getSysconfDir >>= noting "getSysconfDir"
            getCurrentDirectory >>= noting "CurrentDir"
            inlineData
    in  addDependentFile filePath *>
            (location >>= runIO . noting "location" . show)
            *> runIO workingIO
{-
licenseText :: ExpQ
licenseText =
    let locateData = runIO $ getDataFileName "LICENSE"
        gatherData = polishData . readFile . normalise
        polishData = fmap (intercalate "\n" . lines)
    in  do  filePath <- locateData
            addDependentFile filePath
            runIO $ gatherData filePath >>= lift 
-}
