{- |
Compile-time embeddings of the ASCII art "splash image" via Template Haskell.
-}

{-# Language OverloadedStrings #-}
{-# Language Strict #-}
{-# Language TemplateHaskellQuotes #-}

module Software.SplashImage
  ( splashImage
  ) where

import Data.Text (Text)
import Language.Haskell.TH (Code, Q)
import Prettyprinter (LayoutOptions(..), PageWidth(..), indent, layoutPretty, vsep)
import Prettyprinter.Render.Text (renderStrict)


{- |
The ASCII art "splash image" for the software.
-}
splashImage :: Code Q Text
splashImage =
  [|| renderStrict . layoutPretty (LayoutOptions Unbounded) . indent 1 $ vsep
        [ " _____ _          ___ _  _   _   ___ ___   ___          _        _"
        , "|_   _| |_  ___  | _ \\ || | /_\\ / __| __| | _ \\_ _ ___ (_)___ __| |_"
        , "  | | | ' \\/ -_) |  _/ __ |/ _ \\ (_ | _|  |  _/ '_/ _ \\| / -_) _|  _|"
        , "  |_| |_||_\\___| |_| |_||_/_/ \\_\\___|___| |_| |_| \\___// \\___\\__|\\__|"
        , "  ___                     _        _                 |__/"
        , " | _ \\_ _ ___ ___ ___ _ _| |_ ___ (_)"
        , " |  _/ '_/ -_|_-</ -_) ' \\  _(_-<  _"
        , " |_| |_| \\___/__/\\___|_||_\\__/__/ (_)"
        , ""
        , "     ____  __          __                            __  _"
        , "    / __ \\/ /_  __  __/ /___  ____ ____  ____  ___  / /_(_)____"
        , "   / /_/ / __ \\/ / / / / __ \\/ __ `/ _ \\/ __ \\/ _ \\/ __/ / ___/"
        , "  / ____/ / / / /_/ / / /_/ / /_/ /  __/ / / /  __/ /_/ / /__"
        , " /_/   /_/ /_/\\__, /_/\\____/\\__, /\\___/_/ /_/\\___/\\__/_/\\___/"
        , "             /____/        /____/"
        , "    ______                 _         __   ____  __          ______   _"
        , "   / ____/________ _____  / /_     _/_/  / __ \\/ /_  __  __/ ____/  | |"
        , "  / / __/ ___/ __ `/ __ \\/ __ \\   / /   / /_/ / __ \\/ / / / / __    / /"
        , " / /_/ / /  / /_/ / /_/ / / / /  / /   / ____/ / / / /_/ / /_/ /   / /"
        , " \\____/_/   \\__,_/ .___/_/ /_/  / /   /_/   /_/ /_/\\__, /\\____/  _/_/"
        , "                /_/             |_|               /____/        /_/"
        ]
    ||]
