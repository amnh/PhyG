{-# LANGUAGE CPP #-}
{-# LANGUAGE DeriveLift #-}
{-# LANGUAGE Safe #-}
{-# LANGUAGE TemplateHaskellQuotes #-}
{-# OPTIONS_GHC -Wno-implicit-lift #-}


{- |
Compile-time embedding of the UTC time at which the program was built.
-}

#ifdef ENFORCE_TIMESTAMP
{-# OPTIONS_GHC -fforce-recomp #-}
#endif

module Software.Metadata.TimeStamp (
    compilationTimeStamp,
    renderTimeStampAsLocalTime,
) where

import Data.String (IsString (fromString))
import Data.Time.Clock (getCurrentTime)
import Data.Time.Format
import Data.Time.LocalTime (getTimeZone, utcToLocalTime)
import Language.Haskell.TH.Syntax


{- |
The UTC system time at which (this module of) the binary was compiled.
-}
compilationTimeStamp ∷ Code Q String
compilationTimeStamp = (show <$> runIO getCurrentTime) `bindCode` (\time → [||time||])


renderTimeStampAsLocalTime ∷ (IsString s) ⇒ String → IO s
renderTimeStampAsLocalTime timeStr =
    let timeStamp = read timeStr
        formatTimeLocal = formatTime defaultTimeLocale "%Y-%m-%d @ %T"
        formatTimeZoned = formatTime defaultTimeLocale "%EZ"
        formatTimeStamp zone =
            let time = utcToLocalTime zone timeStamp
            in  unwords [formatTimeLocal time, formatTimeZoned zone]
    in  fromString . formatTimeStamp <$> getTimeZone timeStamp
