-----------------------------------------------------------------------------
-- |
-- Module      :  PCG.Syntax
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-- Exports all the PCG syntax stuff.
--
-- I think that eventually the Combinators module should *not* be export from
-- this module, just the stream parser and command types.
--
-----------------------------------------------------------------------------

module PCG.Syntax
    ( module PCG.Syntax.Combinators
    , module PCG.Syntax.Parser
    ) where

import PCG.Syntax.Combinators
import PCG.Syntax.Parser
