{- |
Compile-time embeddings of project contributions via Template Haskell.
-}

{-# Language OverloadedStrings #-}
{-# Language Strict #-}

module Software.License
  ( licenseText
  ) where

import Control.Monad (filterM)
import Data.Foldable
import Data.Text (intercalate, lines)
import Data.Text.IO (readFile)
import Instances.TH.Lift ()
import Language.Haskell.TH hiding (Inline)
import Language.Haskell.TH.Syntax hiding (Inline)
--import Paths_PhyG (getDataFileName)
import Paths_PhyG
import Prelude hiding (lines, readFile)
import System.FilePath ((</>), normalise)
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
        findLicense = getDirFiltered (const (pure True))  -- (pure . (== "LICENSE"))
        
        workingIO = do
            getDataFileName "LICENSE" >>= noting "getDataFileName"
            getBinDir     >>= noting "getBinDir"
            getLibDir     >>= noting "getLibDir"
            getDynLibDir  >>= noting "getDynLibDir"
            getDataDir    >>= noting "getDataDir"
            getLibexecDir >>= noting "getLibexecDir"
            getSysconfDir >>= noting "getSysconfDir"
            cwd <- getCurrentDirectory
            noting "CurrentDir" cwd
            listDirectory cwd >>= noting "listDirectory" . show
            found <- findLicense cwd
            noting "findLicense" $ show found
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


{-# INLINE getDirFiltered #-}
-- | Recursively get all files and subdirectories in the given directory that
-- satisfy the given predicate. Note that the content of subdirectories not
-- matching the filter is ignored. In particular, that means something like
-- @getDirFiltered doesFileExist@ will /not/ recursively return all files.
--
-- @since 0.2.2.0
getDirFiltered :: (FilePath -> IO Bool) -- ^ Filepath filter
               -> FilePath
               -> IO [FilePath]
getDirFiltered check path =
    let mkRel = (path </>)
        foldMapA = (fmap fold .) . traverse
        canDecend fp = do
            isDir <- doesDirectoryExist fp
            isDig <- searchable <$> getPermissions fp
            pure $ isDir && isDig

    in  do  all' <- listDirectory path
            curr <- filterM check (mkRel <$> all')
            dirs <- filterM canDecend all'
            case dirs of
                [] -> pure curr
                ds -> do
                    next <- foldMapA (getDirFiltered check) ds
                    pure $ curr <> next


