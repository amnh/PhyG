-----------------------------------------------------------------------------
-- |
-- Module      :  Data.FileSource.IO
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-- Exposes data input utility which more granualarly controls I/O.
--
-----------------------------------------------------------------------------

{-# Language LambdaCase #-}
{-# Language OverloadedStrings #-}

module Data.FileSource.IO
    ( -- * Inputing Text
      readFile
    , readFiles
    , readSTDIN
      -- * Output Streams
    , FileStream ()
    , appendFile
    , streamBytes
    , streamText
    , writeFile
    , writeFileWithMove
    , writeSTDOUT
      -- * Binary data I/O
    , deserializeBinary
    , serializeBinary
      -- * Error types of I/O
    , InputStreamError ()
    , OutputStreamError ()
    , ParseStreamError ()
    ) where

import Control.DeepSeq
import Control.Exception
import Control.Monad
import Control.Monad.IO.Class
import Control.Monad.Trans.Validation
import Data.Bifunctor
import Data.Binary (Binary, decodeFileOrFail, encode)
import Data.ByteString.Lazy (ByteString)
import Data.ByteString.Lazy qualified as BS
import Data.Char (isNumber)
import Data.FileSource
import Data.FileSource.InputStreamError
import Data.FileSource.OutputStreamError
import Data.FileSource.ParseStreamError
import Data.Foldable
import Data.List (isPrefixOf)
import Data.List.NonEmpty (NonEmpty(..))
import Data.MonoTraversable
import Data.String
import Data.Text.Lazy (Text)
import Data.Text.Lazy.IO qualified as T
import Data.Validation
import Pipes (for, runEffect, yield)
import Prelude hiding (appendFile, readFile, writeFile)
import System.Directory
import System.FilePath.Glob
import System.FilePath.Posix (takeDirectory, takeExtension)
import System.IO hiding (appendFile, readFile, writeFile)
import System.IO.Error


-- |
-- Read textual the contents of a file.
--
-- If the 'FileSource' exists verbatim, the contents will be read.
--
-- If the 'FileSource' does not exist verbatim, the 'FileSource' will be
-- interpreted as a "file-glob" pattern. If there exists a single, unabiguous
-- file matching the "file-glob" pattern, the contents will be read. If there
-- exist multiple files which match the "file-glob" pattern, an "ambiguous file"
-- failure state will be returned.
readFile :: FileSource -> ValidationT InputStreamError IO Text
readFile fSource = readFilesAndLocate readFileContent (invalid . makeAmbiguousFiles fSource) fSource


-- |
-- Read the textual contents of one or more files matching a "file globbing" pattern.
--
-- If the 'FileSource' exists verbatim, the contents will be read.
--
-- If the 'FileSource' does not exist verbatim, the 'FileSource' will be
-- interpreted as a "file-glob" pattern. The contents of each each file matching
-- the glob pattern will be read, and the name of the matching file will be
-- tagged to the corresponding file contents.
readFiles :: FileSource -> ValidationT InputStreamError IO (NonEmpty (FileSource, Text))
readFiles = readFilesAndLocate (fmap pure . readContentsAndTag) (traverse readContentsAndTag)
    where
        readContentsAndTag :: FileSource -> ValidationT InputStreamError IO (FileSource, Text)
        readContentsAndTag path = (\z -> (path, z)) <$> readFileContent path


-- |
-- Read textual the contents of a file.
--
-- Returns a 'Failure' if the STDIN stream is empty.
readSTDIN :: ValidationT InputStreamError IO Text
readSTDIN = do
    nonEmptyStream <- liftIO $ hReady stdin
    if nonEmptyStream then liftIO T.getContents else invalid $ makeEmptyFileStream "STDIN"


-- |
-- Represents a stream of data, either textual or of raw byte data.
--
-- A 'FileStream' will be lazily rendered to it's output source in constant memory.
--
-- Create a 'FileStream' with
--
--   * 'streamBytes'
--   * 'streamText'
--
-- Render a stream with
--
--   * 'appendFile'
--   * 'writeFile'
--   * 'writeFileWithMove'
--   * 'writeSTDOUT'
--
data  FileStream
    = T Text
    | B ByteString


-- |
-- Convert a /lazy/ 'Text' stream to be used in output streaming functions.
{-# INLINE streamText #-}
streamText :: Text -> FileStream
streamText = T


-- |
-- Convert a /lazy/ 'ByteString' stream to be used in output streaming functions.
{-# INLINE streamBytes #-}
streamBytes :: ByteString -> FileStream
streamBytes = B


-- |
-- Write textual stream to the /end/ of the file.
appendFile :: FileSource -> FileStream -> ValidationT OutputStreamError IO ()
appendFile fSource str = ValidationT $ catch
    (Success <$> streamToFile AppendMode fSource str)
    (runValidationT . outputErrorHandling fSource)


-- |
-- Write textual stream to the to the file, overwriting any existing file contents.
writeFile :: FileSource -> FileStream -> ValidationT OutputStreamError IO ()
writeFile fSource str = ValidationT $ catch
    (Success <$> streamToFile WriteMode fSource str)
    (runValidationT . outputErrorHandling fSource)


-- |
-- Write textual stream to the to the file.
--
-- If the specified file already exists, rename the existing file, so that the
-- specified file can be written to without overwriting existing data.
--
-- The existing file is renamed, adding a numeric suffix to the end. The
-- function will try to rename the existing file path by adding the suffix ".0",
-- however if that filepath also exists, it will add ".1", ".2", ".3", ",.4", etc.
-- The suffix added will be one greater than the highest existing numeric suffix.
writeFileWithMove :: FileSource -> FileStream -> ValidationT OutputStreamError IO ()
writeFileWithMove fSource str = liftIO (safelyMoveFile fSource) *> writeFile fSource str


-- |
-- Render the stream to STDIN.
writeSTDOUT :: FileStream -> ValidationT OutputStreamError IO ()
writeSTDOUT = liftIO . \case
    T s -> T.putStr s
    B s -> BS.putStr s


-- |
-- Deserialize binary encodable content from the specified file path.
--
-- Operational inverse of 'serializeBinary'.
deserializeBinary :: Binary a => FileSource -> ValidationT (Either InputStreamError ParseStreamError) IO a
deserializeBinary fSource =
    let deserialize :: Binary b => FileSource -> ValidationT (Either InputStreamError ParseStreamError) IO b
        deserialize fs = do
            res <- ValidationT $ catch
                (Success <$> decodeFileOrFail (otoList fs))
                (fmap (first Left) . runValidationT . inputErrorHandling fs)
            case res of
                Right val        -> pure val
                Left  (off, err) -> invalid . Right . makeDeserializeErrorInBinaryEncoding fs $ fold
                    ["At stream offset ", fromString $ show off, ", ", fromString err]
    in  readFilesAndLocate' Left deserialize (invalid . Left . makeAmbiguousFiles fSource) fSource


-- |
-- Serialize binary encodable content to the specified file path.
--
-- Operational inverse of 'deserializeBinary'.
serializeBinary :: Binary a => FileSource -> a -> ValidationT OutputStreamError IO ()
serializeBinary fSource val = ValidationT $ catch
    (fmap Success . streamToFile WriteMode fSource . streamBytes $ encode val)
    (runValidationT . outputErrorHandling fSource)


-- |
-- Read the textual contents of one or more files matching a "file globbing" pattern.
readFilesAndLocate
    :: (FileSource -> ValidationT InputStreamError IO a)          -- ^ What to do in the unambiguous case
    -> (NonEmpty FileSource -> ValidationT InputStreamError IO a) -- ^ What to do in the ambiguous case
    -> FileSource -- ^ Path description
    -> ValidationT InputStreamError IO a -- ^ Result
readFilesAndLocate = readFilesAndLocate' id


-- |
-- Read the textual contents of one or more files matching a "file globbing" pattern.
readFilesAndLocate'
    :: Semigroup e
    => (InputStreamError -> e)                     -- ^ How to project an input error to the type e
    -> (FileSource -> ValidationT e IO a)          -- ^ What to do in the unambiguous case
    -> (NonEmpty FileSource -> ValidationT e IO a) -- ^ What to do in the ambiguous case
    -> FileSource -- ^ Path description
    -> ValidationT e IO a -- ^ Result
readFilesAndLocate' e f g fSource = do
    -- Check if the file exists exactly as specified
    exists <- liftIO $ doesFileExist (otoList fSource)
    if exists
    -- If it exists exactly as specified, read it in
        then f fSource
        else do
        -- If the file does not exists exactly as specified
        -- try to match other files to the given path
        -- by interpreting the path as a 'glob'
            matches <- fmap (fmap fromString) . liftIO . glob $ otoList fSource
            case matches of
                []     -> invalid . e $ makeFileNotFound fSource
                [ x ]  -> f x
                x : xs -> g $ x :| xs


-- |
-- Shared utility between 'readFile' & 'readFiles'.
--
-- Checks if the file permissions allows the file contents to be read.
readFileContent :: FileSource -> ValidationT InputStreamError IO Text
readFileContent fSource =
    let path = force $ otoList fSource
    in  do
            canRead <- liftIO $ readable <$> getPermissions path
            if not canRead
                then invalid $ makeFileNoReadPermissions fSource
                else do
                    txt <- ValidationT $ catch
                        (Success <$> streamfromFile fSource)
                        (runValidationT . inputErrorHandling fSource)
                    if onull txt then invalid $ makeEmptyFileStream fSource else pure txt


-- |
-- Smartly handle certain I/O errors that can occur while inputing a data stream.
--
-- Re-throws errors not specially handled and reported by 'InputStreamError'.
inputErrorHandling :: FileSource -> IOError -> ValidationT InputStreamError IO a
inputErrorHandling fSource e
    | isAlreadyInUseError e = invalid $ makeFileInUseOnRead fSource
    | isPermissionError e   = invalid $ makeFileNoReadPermissions fSource
    | isDoesNotExistError e = invalid $ makeFileNotFound fSource
    |
  -- Re-throw if it is not an error we explicitly handle and report
      otherwise             = ValidationT $ ioError e


-- |
-- Smartly handle certain I/O errors that can occur while outputting a data stream.
--
-- Re-throws errors not specially handled and reported by 'OutputStreamError'.
outputErrorHandling :: FileSource -> IOError -> ValidationT OutputStreamError IO a
outputErrorHandling fSource e
    | isAlreadyInUseError e = invalid $ makeFileInUseOnWrite fSource
    | isFullError e         = invalid $ makeNotEnoughSpace fSource
    | isPermissionError e   = invalid $ makeFileNoWritePermissions fSource
    | isDoesNotExistError e = invalid $ makePathDoesNotExist fSource
    |
  -- Re-throw if it is not an error we explicitly handle and report
      otherwise             = ValidationT $ ioError e


-- |
-- Checks to see if the supplied file path exists.
--
-- If it does, it moves the existing file path, so that the supplied file path
-- can be written to without overwriting data.
--
-- The existing file path is renamed, adding a numeric suffix to the end. The
-- function will try to rename the existing file path by adding the suffix ".0",
-- however if that filepath also exists, it will add ".1", ".2", ".3", ",.4", etc.
-- The suffix added will be one greater than the highest existing numeric suffix.
safelyMoveFile :: FileSource -> IO ()
safelyMoveFile fSource =
    let filePath            = otoList fSource
        getFilePathPrefixes = fmap (drop (length filePath)) . filter (filePath `isPrefixOf`)
        getNumericSuffixes =
            let hasDotThenNumberSuffix ('.' : x : xs) = all isNumber $ x : xs
                hasDotThenNumberSuffix _              = False
            in  fmap tail . filter hasDotThenNumberSuffix . fmap takeExtension

        getLargestNumericSuffix :: (Num a, Ord a, Read a) => [String] -> a
        getLargestNumericSuffix []       = -1
        getLargestNumericSuffix (x : xs) = maximum . fmap read $ x :| xs
    in  do
        absPath <- makeAbsolute filePath
        exists  <- doesFileExist absPath
        when exists $ do
            allFiles <- getDirectoryContents $ takeDirectory absPath
            let prefixed = getFilePathPrefixes allFiles
            let numbers  = getNumericSuffixes prefixed
            let lastNum  = getLargestNumericSuffix numbers
            let nextNum  = lastNum + 1 :: Word
            let newName  = absPath <> "." <> show nextNum
            renameFile absPath newName


-- |
-- Streams text from a file in constant memory.
streamfromFile :: FileSource -> IO Text
streamfromFile = T.readFile . otoList
-- This streaming from file doesn't work...
{-
streamfromFile fSource = do
    h   <- openFile (otoList fSource) ReadMode
    txt <- runEffect $ (liftIO (T.hGetLine h)) >~ await
    hClose h
    pure txt
-}


-- |
-- Streams text to a file in constant memory.
streamToFile :: IOMode -> FileSource -> FileStream -> IO ()
streamToFile m fSource = \case
    T txt -> runStream T.hPutStr m fSource txt
    B bts -> runStream BS.hPutStr m fSource bts


-- |
-- Given a streaming function to a file handle, write out a data stream.
runStream :: (Handle -> v -> IO ()) -> IOMode -> FileSource -> v -> IO ()
runStream f mode fSource v = do
    h <- openFile (otoList fSource) mode
    runEffect $ for (yield v) (liftIO . f h)
    hClose h
