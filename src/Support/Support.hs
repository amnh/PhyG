{- |
Module      :  Support.hs
Description :  Module containg support functions
Copyright   :  (c) 2022 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
License     :

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

Maintainer  :  Ward Wheeler <wheeler@amnh.org>
Stability   :  unstable
Portability :  portable (I hope)

-}

module Support.Support  ( supportGraph
                        ) where

import Types.Types
import qualified ParallelUtilities       as PU
import Control.Parallel.Strategies
import Debug.Trace
import GeneralUtilities
import qualified Graphs.GraphOperations  as GO
import Utilities.Utilities               as U
import Data.Maybe
import           Data.Char
import           Text.Read

-- | refinement arguments
supportArgList :: [String]
supportArgList = ["netadd", "netdel", "netdelete", "netadddel", "netmove","geneticalgorithm", "ga"] ++ fuseArgList ++ netEdgeArgList ++ geneticAlgorithmArgList

-- | driver for overall support
supportGraph :: [Argument] -> GlobalSettings -> ProcessedData -> Int -> [PhylogeneticGraph] -> [PhylogeneticGraph]
supportGraph inArgs inGS inData rSeed inGraphList = 
   if null inGraphList then error ("No graphs input to calculate support") 
   else inGraphList
     