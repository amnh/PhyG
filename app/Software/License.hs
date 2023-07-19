{- |
Compile-time embeddings of project contributions via Template Haskell.
-}

{-# Language LambdaCase #-}
{-# Language OverloadedStrings #-}
{-# Language Strict #-}

module Software.License
  ( licenseText
  ) where

import Data.Map ((!))
import Data.Text (Text)
import Software.Metadata.Embedded (embeddedDataFiles)
{-
import Control.Monad (filterM)
import Data.FileEmbed
import Data.Foldable
import Data.List (isPrefixOf, sort)
import Data.List.NonEmpty (init, last, nonEmpty)
import Data.Text (intercalate, lines)
import Data.Text.IO (readFile)
import Instances.TH.Lift ()
import Language.Haskell.TH hiding (Inline)
import Language.Haskell.TH.Syntax hiding (Inline)
--import Paths_PhyG (getDataFileName)
import Paths_PhyG
import Prelude hiding (init, last, lines, readFile)
import System.FilePath ((</>), normalise, splitDirectories, splitFileName, takeFileName)
import System.Directory
-}

licenseText :: Text
licenseText = embeddedDataFiles ! "LICENSE"

{-
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
        findLicense =
            let predicate path =
                    let (fPath, fName) = splitFileName path
                        allDirs = nonEmpty $ splitDirectories fPath
                        check dirs = last dirs == "doc" && ("PhyG" `isPrefixOf`) `any` init dirs
                    in  fName == "LICENSE" && maybe False check allDirs

            in  getFilesFilteredBy $ pure . predicate

        searchPaths =
            let pred = \case
                    '.':_ -> pure False
                    fName -> doesDirectoryExist fName
            in  getCurrentDirectory >>= listDirectory >>= filterM pred

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
            searchPaths >>= noting "searchPaths" . show
            found <- searchPaths >>= findLicense
            noting "findLicense" $ unlines found
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


{-# INLINE getFilesFilteredBy #-}
-- | Recursively get all files and subdirectories in the given directory that
-- satisfy the given predicate. Note that the content of subdirectories not
-- matching the filter is ignored. In particular, that means something like
-- @getDirFiltered doesFileExist@ will /not/ recursively return all files.
--
-- @since 0.2.2.0
getFilesFilteredBy
  :: (FilePath -> IO Bool) -- ^ File filter
  -> [FilePath] -- ^ Input paths
  -> IO [FilePath]
getFilesFilteredBy predicate = foldMapA (getFilesFilteredBy' predicate)


foldMapA :: (FilePath -> IO [FilePath]) -> [FilePath] -> IO [FilePath]
foldMapA = (fmap fold .) . traverse


{-# INLINE getFilesFilteredBy' #-}
getFilesFilteredBy'
  :: (FilePath -> IO Bool) -- ^ Filepath filter
  -> FilePath
  -> IO [FilePath]
getFilesFilteredBy' check path =
    let prepend = (path </>)
        canDecend fp = do
            isDir <- doesDirectoryExist fp
            perms <- getPermissions fp
            let result = isDir && readable perms && searchable perms
            let action
                    | result    = "SCAN"
                    | otherwise = "PASS"
            noting fp "permissions" $ show perms
            noting fp "Decendable" action
            pure result

        noting file key val = putStrLn $ unwords [ "Dir Scan @", file, key <> ":\t", val ]

        consider fp = do
            isFile <- doesFileExist fp
            if not isFile
            then noting fp "Consideration" "SKIP" *> pure False
            else do
                isGood <- check fp
                let action
                      | isGood    = "KEEP"
                      | otherwise = "SKIP"
                noting fp "Consideration" action
                pure $ isFile && isGood

    in  do  all' <- fmap prepend . sort <$> listDirectory path
            curr <- filterM consider  all'
            dirs <- filterM canDecend all'
            case dirs of
                [] -> pure curr
                ds -> do
                    next <- foldMapA (getFilesFilteredBy' check) ds
                    pure $ curr <> next
-}
