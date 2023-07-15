{- |
Record which holds all the resulting preferences of the user as indicated by the
supplied command line options.
-}

{-# Language DeriveGeneric #-}
{-# Language GeneralizedNewtypeDeriving #-}
{-# Language Strict #-}

module CommandLineOptions.Types
  ( -- * Primary data-types
    CommandLineOptions(..)
  , DisplayInfoBlock(..)
    -- * Secondary data-types
    -- ** CommandLineOptionsInput
  , CommandLineOptionsInput(..)
    -- ** DisplayFlags
  , DisplayFlags(DisplayFlags)
    -- *** Accessor
  , displaySomeInformation
    -- ** InputFile
  , InputFile()
    -- *** Constructor
  , missingInputFile
    -- *** Accessor
  , maybeInputFile
  ) where

import Control.DeepSeq
import Data.Foldable (fold)
import Data.List.NonEmpty (NonEmpty, nonEmpty)
import Data.String (IsString(..))
import GHC.Generics
import Text.Read(Read(readPrec))


{- |
Valid command line options
-}
newtype CommandLineOptions = OptionsCLI { getArgsCLI :: Either (NonEmpty DisplayInfoBlock) FilePath }
    deriving stock (Generic)


{- |
Collection of flags gathered from command line interface which has yet to be validated.
-}
data  CommandLineOptionsInput
    = CommandLineOptionsInput
    { inputFile    :: InputFile
    , displayFlags :: DisplayFlags
    } deriving stock (Generic)


{- |
Collection of binary flags indicating requests for software metadata information to be displayed.
-}
data  DisplayFlags
    = DisplayFlags
    { printCredits :: Bool
    , printLicense :: Bool
    , printSplash  :: Bool
    , printVersion :: Bool
    } deriving stock (Generic)


{- |
Enumeration of metadata information blocks which can be displayed.

/Ordering of this type determines the ordering of information displayed!/
-}
data  DisplayInfoBlock
    = DisplayVersion
    | DisplayLicense
    | DisplaySplash
    | DisplayCredits
    deriving stock (Eq, Ord, Generic)


{- |
Possibly missing input file for the software.
-}
newtype InputFile = InputFile { maybeInputFile :: Maybe FilePath }
    deriving stock   (Generic)
    deriving newtype (Eq, Ord, Show)


instance IsString InputFile where

    fromString []  = InputFile Nothing
    fromString fp = InputFile $ Just fp


instance NFData CommandLineOptions


instance NFData CommandLineOptionsInput


instance NFData DisplayFlags


instance NFData DisplayInfoBlock


instance NFData InputFile


instance Read InputFile where

    readPrec = InputFile . Just <$> readPrec


{- |
Accessor for 'DisplayFlags'.
-}
displaySomeInformation :: DisplayFlags -> Maybe (NonEmpty DisplayInfoBlock)
displaySomeInformation =
    let coalesce
          :: ( Applicative f
             , Monoid (f DisplayInfoBlock)
             )
          => (DisplayFlags -> Bool) -> DisplayInfoBlock -> DisplayFlags -> f DisplayInfoBlock
        coalesce check val obj
            | check obj = pure val
            | otherwise = mempty

        gatherSpecifiedFlags = fold . (getFlags <*>) . pure

        getFlags = 
            [ coalesce printVersion DisplayVersion
            , coalesce printLicense DisplayLicense
            , coalesce printSplash  DisplaySplash
            , coalesce printCredits DisplayCredits
            ]

    in  nonEmpty . gatherSpecifiedFlags


{- |
Represents the default value for when no input file exists.
-}
missingInputFile :: InputFile
missingInputFile = InputFile Nothing
