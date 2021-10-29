{- |
Module      :  ReconcileGraphs.hs
Description :  Module to call graph reconciliation functions
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

module Reconciliation.ReconcileGraphs  ( makeReconcileGraph 
                                       ) where

import           Types.Types
import qualified Reconciliation.Eun    as E
import qualified Utilities.LocalGraph  as LG
import qualified GraphFormatUtilities  as GFU


-- | makeReconcileGraph is a wrapper around eun.hs functions to return String of reconciled graph
makeReconcileGraph :: [String] -> [SimpleGraph] -> (String, SimpleGraph)
makeReconcileGraph commandList inGraphList =
   if null inGraphList then "Error: No input graphs to reconcile"
   else   
      let -- convert SimpleGraph to String String from Text Double
          stringGraphs = fmap (GFU.modifyVertexEdgeLabels True True) fmap GFU.textGraph2StringGraph inGraphList
          method = "eun"
          compareMethod = "combinable"
          threshold = 100
          connectComponents = True
          edgeLabel = True
          vertexLabel = False
          outputFormat = "dot"
          outputFile = "bleh.dot"
          (reconcileString, reconcileGraph) = reconcile (method, compareMethod, threshold, connectComponents, edgeLabel, vertexLabel, outputFormat, outputFile, stringGraphs)
          reconcileSimpleGraph = GFU.stringGraph2TextGraphDouble reconcileGraph
      in
      (reconcileString, reconcileSimpleGraph)





