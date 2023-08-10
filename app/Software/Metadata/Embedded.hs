{- |
Compile-time embedding of the UTC time at which the program was built.
-}

{-# Language TemplateHaskell #-}

module Software.Metadata.Embedded
  ( embeddedDataFiles
  ) where

import Data.Bifunctor (bimap)
import Data.ByteString (ByteString)
import Data.Map (Map)
import Data.Text (Text)
import Data.Text.Encoding (decodeUtf8Lenient)
import Data.FileEmbed
import GHC.Exts (fromListN)
import System.FilePath (splitDirectories)
import System.FilePath.Posix (joinPath)


embeddedDataFiles :: Map FilePath Text
embeddedDataFiles =
    let -- We ensure that the data files are stored using Posix
        -- path separators (/), even on Windows.
        correctPath = joinPath . splitDirectories
        correctData = decodeUtf8Lenient
        alterations = bimap correctPath correctData
    in  fromListN 3 $ alterations <$> dataFiles'


dataFiles' :: [(FilePath, ByteString)]
dataFiles' =
    [ ( "Authors.md", $(embedFile "doc/Authors.md") )
    , ( "Funding.md", $(embedFile "doc/Funding.md") )
    , ( "LICENSE"   , $(embedFile "doc/LICENSE")    )
    ]
