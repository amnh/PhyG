{- |
Compile-time embedding of the UTC time at which the program was built.
-}

{-# Language Safe #-}
{-# Language TemplateHaskellQuotes #-}
{-# OPTIONS_GHC -fforce-recomp #-}
{-# OPTIONS_GHC -Wno-implicit-lift #-}

module Software.Metadata.TimeStamp
  ( compilationTimeStamp
  ) where

import Data.Time.Clock (getCurrentTime)
import Data.Time.Format
import Language.Haskell.TH.Syntax


{- |
The UTC system time at which (this module of) the binary was compiled.
-}
compilationTimeStamp :: Code Q String
compilationTimeStamp = (formatTime defaultTimeLocale "%Y-%m-%d %T UTC" <$> runIO getCurrentTime) `bindCode` (\time -> [|| time ||] )
