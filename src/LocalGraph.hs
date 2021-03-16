{- |
Module      :  LocalGraph.hs
Description :  Module specifying graph types and functionality
				This is for indirection so can change underlying graph library
				without  polutting the rest of the code
Copyright   :  (c) 2021 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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

{-# LANGUAGE ScopedTypeVariables #-}

module LocalGraph  where


import qualified Data.Graph.Inductive.PatriciaTree as P
import GraphFormatUtilities
import qualified Data.Text as T
import           Data.GraphViz                     as GV
import           Data.GraphViz.Attributes.Complete (Attribute (Label),
                                                    Label (..))
import           Data.GraphViz.Commands.IO         as GVIO
import qualified Data.Graph.Inductive.Graph        as G


-- | Gr local graph definition using FGL
type Gr a b = P.Gr a b
type Node = G.Node
type DotGraph = GV.DotGraph

-- | getFENLocal maps to forestEnhancedNewickStringList2FGLList in GraphFormatUtilities
-- to allow for potnetial swapping FGL graph backend
-- requires leading and trailing space and newlines to be removed
getFENLocal :: T.Text -> [Gr T.Text Double] 
getFENLocal inText = forestEnhancedNewickStringList2FGLList inText


readDotLocal :: String -> IO [DotGraph G.Node]
readDotLocal fileName = do
	(dotGraph :: [DotGraph G.Node]) <- GVIO.readDotFile fileName
    return dotGraph

-- let (inputGraphListDot :: [P.Gr Attributes Attributes]) = fmap dotToGraph  dotGraphList

-- | dotToGraph local mapo dor 
dotToGraph ::  [LocalGraph.DotGraph Node] -> [LocalGraph.Gr Attributes Attributes]
dotToGraph dotGraphList = 
	let (outGraph :: [P.Gr Attributes Attributes]) = GV.dotToGraph dotGraphList
	in outGraph
