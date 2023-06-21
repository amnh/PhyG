------------------------------------------------------------------------------
-- |
-- Module      :  Software.Metadata.TimeStamp
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# Language OverloadedStrings #-}
{-# Language TemplateHaskell   #-}

module Software.Metadata.TimeStamp
  ( compilationTimeStamp
  ) where

import Data.Time.Clock (getCurrentTime)
import Data.Time.Format.ISO8601
import Language.Haskell.TH.Syntax


{- |
The UTC system time at which (this module of) the binary was compiled.
-}
compilationTimeStamp :: Code Q String
compilationTimeStamp = (formatShow iso8601Format <$> runIO getCurrentTime) `bindCode` (\time -> [|| time ||])
