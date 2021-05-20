{- |
Module      :  Traversals.hs
Description :  Module specifying graph traversal functions for PhyGraph
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

module Traversals  (  fullyLabelGraph
                    , postOrderTreeTraversal
                    , preOrderTreeTraversal
                   ) where

import           Types
import qualified Data.Vector as V
import qualified LocalGraph as LG
import GeneralUtilities
import Debug.Trace
import qualified GraphFormatUtilities as GFU

-- | emptyPhyloGeneticGraph for checking inputs
emptyPhyloGeneticGraph :: PhylogeneticGraph
emptyPhyloGeneticGraph = (LG.empty, 0.0, V.empty, V.empty, (V.empty, V.empty))

-- | fullyLabelGraph takes an unlabelled "simple' graph, performs post and preorder passes to 
-- fully label the graph and return a PhylogeenticGraph
fullyLabelGraph :: GlobalSettings -> ProcessedData -> SimpleGraph -> PhylogeneticGraph
fullyLabelGraph inGS inData inGraph = 
    if LG.isEmpty inGraph then (LG.empty, 0.0, V.empty, V.empty, inData)
    else 
        let postOrderTree = postOrderTreeTraversal inGS inData inGraph
            preOrderTree = preOrderTreeTraversal inGS inData postOrderTree
        in
        (inGraph, 0.0, V.empty, V.empty, inData)


-- | postOrderTreeTraversal takes a 'simple' graph and generates 'preliminary' assignments
-- vi post-order traversal, yields cost as well
-- for a binary tree only
postOrderTreeTraversal :: GlobalSettings -> ProcessedData -> SimpleGraph -> PhylogeneticGraph
postOrderTreeTraversal inGS inData inGraph = 
    if LG.isEmpty inGraph then (LG.empty, 0.0, V.empty, V.empty, inData)
    else
        -- Assumes root is Number of Leaves  
    	let rootIndex = V.length $ fst inData
    	in
    	trace ("It Begins at " ++ show rootIndex) (
    	if not $ LG.isRoot inGraph rootIndex then error ("Index "  ++ (show rootIndex) ++ " not root in graph:\n" ++ (GFU.showGraph inGraph))
    	else emptyPhyloGeneticGraph
    	)


-- | preOrderTreeTraversal takes a preliminarily labelled PhylogeneticGraph
-- and returns a full labbels with 'final' assignments
-- invafiant that root is HTU !! nLeaves?
preOrderTreeTraversal :: GlobalSettings -> ProcessedData -> PhylogeneticGraph -> PhylogeneticGraph
preOrderTreeTraversal inGS inData inPGraph = 
    if LG.isEmpty (fst5 inPGraph) then emptyPhyloGeneticGraph
    else 
    	inPGraph
    	

        



