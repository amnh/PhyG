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

import qualified Data.List            as L
import qualified Data.Text.Lazy       as T
import           GeneralUtilities
import qualified GraphFormatUtilities as GFU
import qualified Reconciliation.Eun   as E
import           Types.Types
import qualified Utilities.LocalGraph as LG
import           Data.Maybe      
-- import           Debug.Trace

-- | makeReconcileGraph is a wrapper around eun.hs functions to return String of reconciled graph
makeReconcileGraph :: [String] -> [(String, String)] -> [SimpleGraph] -> (String, SimpleGraph)
makeReconcileGraph validCommandList commandPairList inGraphList =
   if null inGraphList then ("Error: No input graphs to reconcile", LG.empty)
   else
      let -- convert SimpleGraph to String String from Text Double
          stringGraphs = fmap (GFU.modifyVertexEdgeLabels True True . GFU.textGraph2StringGraph) inGraphList

          -- parse arguements
          commandList = (mergePair <$> filter (('"' `notElem`).snd) commandPairList)
          (localMethod, compareMethod, threshold, connectComponents, edgeLabel, vertexLabel, outputFormat) = processReconcileArgs validCommandList commandList

          -- call EUN/reconcile functions
          (reconcileString, reconcileGraph) = E.reconcile (localMethod, compareMethod, threshold, connectComponents, edgeLabel, vertexLabel, outputFormat,stringGraphs)

          -- convert eun format graph back to SimpleGraph
          reconcileSimpleGraph = GFU.stringGraph2TextGraphDouble reconcileGraph
      in
      --trace ("MRG :" <> (show (localMethod, compareMethod, threshold, connectComponents, edgeLabel, vertexLabel, outputFormat)) <> "\n" <> reconcileString
      --  <> "\n" <> (LG.prettyIndices reconcileSimpleGraph))
      (reconcileString, reconcileSimpleGraph)
      where mergePair (a,b) = if a /= [] && b /= [] then a <> (':' : b)
                              else a <> b


-- | processReconcileArgs takes a list of strings and returns values of commands for proram execution
-- including defaults
-- checks commands for misspellings
processReconcileArgs :: [String] -> [String] -> (String, String, Int, Bool, Bool, Bool, String)
processReconcileArgs validCommandList inList' =
    let inList = inList' L.\\ ["overwrite", "append", "reconcile"]
    in
    if null inList then
      let -- default values
          localMethod = "eun"
          compareMethod = "combinable"
          threshold = 0
          connectComponents = True
          edgeLabel = True
          vertexLabel = True
          outputFormat = "dot"
      in
      (localMethod, compareMethod, threshold, connectComponents, edgeLabel, vertexLabel, outputFormat)

    else
        -- trace ("Rec args: " <> (show inList)) (
        let inTextList = fmap T.pack inList
            inTextListLC = fmap T.toLower inTextList
            commandList = filter (T.any (== ':')) inTextListLC
            stringCommands = fmap (T.unpack . T.takeWhile (/= ':')) commandList
            (editCostList, matchList) = unzip $ fmap (getBestMatch (maxBound :: Int ,"no suggestion") validCommandList) stringCommands
            commandMatch = zip3 editCostList stringCommands matchList
            notMatchedList = filter ((>0).fst3) commandMatch
            localMethod = getMethod inTextListLC
            compareMethod = getCompareMethod inTextListLC
            connect = getConnect inTextListLC
            edgeLabel = getEdgeLabel inTextListLC
            vertexLabel = getVertexLabel inTextListLC
            threshold
              | localMethod == "cun" = 0
              | localMethod == "strict" = 100
              | otherwise = getThreshold inTextListLC
            outFormat = getOutputFormat inTextListLC
        in
        if null notMatchedList then
            (localMethod, compareMethod, threshold, connect, edgeLabel, vertexLabel, outFormat)
        else errorWithoutStackTrace ("\n\nError(s) in reconcile command specification (case insensitive):\n" <> getCommandErrorString notMatchedList)
        -- )

-- | getMethod returns method value or dedfault otherwise
-- assumes in lower case
getMethod :: [T.Text] -> String
getMethod inTextList =
    -- default
    if null inTextList then "eun"
    else
        let firstCommand = T.takeWhile (/= ':') $ head inTextList
            firstOption = T.tail $ T.dropWhile (/= ':') $ head inTextList
        in
        if isNothing (T.find (== ':') (head inTextList)) then getMethod (tail inTextList)
        else if firstCommand == T.pack "method" then
            let option = T.unpack firstOption
            in
            if option == "eun" then "eun"
            else if option == "cun" then "cun"
            else if option == "majority" then "majority"
            else if option == "strict" then "strict"
            else if option == "adams" then "adams"
            else errorWithoutStackTrace ("Reconcile option \'" <> option <> "\' not recognized (eun|cun|majority|strict)")
        else getMethod (tail inTextList)

-- | getCompareMethod returns compareMethod value or default otherwise
-- assumes in lower case
getCompareMethod :: [T.Text] -> String
getCompareMethod inTextList =
    -- default
    if null inTextList then "combinable"
    else
        let firstCommand = T.takeWhile (/= ':') $ head inTextList
            firstOption = T.tail $ T.dropWhile (/= ':') $ head inTextList
        in
        if isNothing (T.find (== ':') (head inTextList)) then getCompareMethod (tail inTextList)
        else if firstCommand == T.pack "compare" then
            let option = T.unpack firstOption
            in
            if option == "combinable" then "combinable"
            else if option == "identity" then "identity"
            else errorWithoutStackTrace ("Compare option \'" <> option <> "\' not recognized (combinable|identity)")
        else getCompareMethod (tail inTextList)

-- | getConect returns connect value or default otherwise (True|False)
-- assumes in lower case
getConnect :: [T.Text] -> Bool
getConnect inTextList =
    -- default
    not (null inTextList) && (let firstCommand = T.takeWhile (/= ':') $ head inTextList
                                  firstOption = T.tail $ T.dropWhile (/= ':') $ head inTextList
                              in
                              if isNothing (T.find (== ':') (head inTextList)) then getConnect (tail inTextList)
                              else if firstCommand == T.pack "connect" then
                                  let option = T.unpack firstOption
                                  in
                                  (option == "true") || (option /= "false" && errorWithoutStackTrace ("Connect option \'" <> option <> "\' not recognized (True|False)"))
                              else getConnect (tail inTextList))

-- | getEdgeLabel returns edgeLabel value or default otherwise (True|False)
-- assumes in lower case
getEdgeLabel :: [T.Text] -> Bool
getEdgeLabel inTextList =
    -- default
    null inTextList || (let firstCommand = T.takeWhile (/= ':') $ head inTextList
                            firstOption = T.tail $ T.dropWhile (/= ':') $ head inTextList
                        in
                        if isNothing (T.find (== ':') (head inTextList)) then getEdgeLabel (tail inTextList)
                        else if firstCommand == T.pack "edgelabel" then
                            let option = T.unpack firstOption
                            in
                            (option == "true") || (option /= "false" && errorWithoutStackTrace ("EdgeLAbel option \'" <> option <> "\' not recognized (True|False)"))
                        else getEdgeLabel (tail inTextList))

-- | getVertexLabel returns edgeLabel value or default otherwise (True|False)
-- assumes in lower case
getVertexLabel :: [T.Text] -> Bool
getVertexLabel inTextList =
    -- default
    not (null inTextList) && (let firstCommand = T.takeWhile (/= ':') $ head inTextList
                                  firstOption = T.tail $ T.dropWhile (/= ':') $ head inTextList
                              in
                              if isNothing (T.find (== ':') (head inTextList)) then getVertexLabel (tail inTextList)
                              else if firstCommand == T.pack "vertexlabel" then
                                  let option = T.unpack firstOption
                                  in
                                  (option == "true") || (option /= "false" && errorWithoutStackTrace ("VertexLabel option \'" <> option <> "\' not recognized (True|False)"))
                              else getVertexLabel (tail inTextList))


-- | getThreshold returns threshold value or default otherwise
-- assumes in lower case
getThreshold :: [T.Text] -> Int
getThreshold inTextList =
    -- default
    if null inTextList then 0 :: Int
    else
        let firstCommand = T.takeWhile (/= ':') $ head inTextList
            firstOption = T.tail $ T.dropWhile (/= ':') $ head inTextList
        in
        if isNothing (T.find (== ':') (head inTextList)) then getThreshold (tail inTextList)
        else if firstCommand == T.pack "threshold" then read (T.unpack firstOption) :: Int
        else getThreshold (tail inTextList)

-- | getOutputFormat returns output file format or default otherwise
-- assumes in lower case
getOutputFormat :: [T.Text] -> String
getOutputFormat inTextList =
    -- default
    if null inTextList then "dot"
    else
         --removed prefix for output graph format
        let firstOption = head inTextList
            outFormat = T.unpack firstOption
         in
         if outFormat == "dot" then "dot"
         else (if (outFormat == "fenewick") || (outFormat == "newick") then "fenewick" else getOutputFormat (tail inTextList))


-- |
