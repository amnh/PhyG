{- |
Module      :  Graphs
Description :  Functions to generate (algorithmic) complexity of graphs
               Generates a Haskell program (compiles)
				       with a description of a graph. The output program can be executed with
				       GHCi interpreter. ghci --:load progName
               also outputs Huffman binary code of program
Copyright   :  (c) 2018-2019 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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

{- |E| = |L| + (|L|-2) - 2(|R|-1) + 3|N| + |S|
       = 2|L| - 2|R| + 3|N| + |S|
   |V| = |L| + (|L|-1) - (|R|-1) + 2|N| + |S|
       = 2|L| - |R| + 2|N| + |S|
    E = edge set, L = leave set, R = root set, N = network edge set,
    S = "singleton" (single leaf + root) components.

    To Do:  read graph from file and do get numbers that way?
-}

module Complexity.Graphs
  (  makeProgramStringGraph
  ,  makeDisplayGraphString
  )  where

import Complexity.CodeStrings
--import Debug.Trace

mainStartString :: String
mainStartString = "main=do\n"

-- | makeProgramString is wrapper for makeGraphProgramString to simplify interface
makeProgramStringGraph :: Int -> Int -> Int -> Int -> String
makeProgramStringGraph  numLeaves numSingle numRoots numNetEdges  =
  let (outProgram, middleString, sumString) = makeGraphProgramString numLeaves numRoots numSingle programStartStringGraph mainStartString "" numLeaves numRoots numNetEdges numSingle
  in
  outProgram  ++ middleString ++ sumString

-- | makeGraphProgramString puts apropriate code and values into string to define a general graph
makeGraphProgramString :: Int -> Int -> Int -> String -> String -> String -> Int -> Int -> Int -> Int -> (String, String, String)
makeGraphProgramString numLeavesOrig numRootsOrig numSingleOrig beforeMain afterMain sumString numLeaves numRoots numNetEdges numSingle =
  let numLeaves2 =  numLeavesOrig - numSingleOrig
      numRoots2 = numRootsOrig - numSingleOrig
      numLeaves3 = numLeaves2 - max 0 (2 * (numRootsOrig - numSingleOrig - 1))
      maxVertex = (2 * numSingleOrig) + max 0 (3 * (numRootsOrig - numSingleOrig - 1))
      maxVertex2 = maxVertex + (2 * numLeaves3) - 1
      doSingle = (numSingleOrig > 0)
      doMinimal = (numRootsOrig > (numSingleOrig + 1))
      doTree = numLeavesOrig > max 0 (2*(numRootsOrig - numSingleOrig - 1)) + numSingleOrig
      doEdges = numNetEdges > 0
  in
  if numSingle > 0 then --make singleton trees one root -> one leaf each
    if doMinimal || doTree then -- More stuff to come
      makeGraphProgramString numLeavesOrig numRootsOrig numSingleOrig (beforeMain ++ getSingletonEdgesString ++ nullString) (afterMain ++ "  let s=aG 0 " ++ show numSingleOrig ++ "\n") "  p0 \"\" (s++"  (numLeaves - numSingle) (numRoots - numSingle) numNetEdges 0
    else --singles only
      makeGraphProgramString numLeavesOrig numRootsOrig numSingleOrig (beforeMain ++ getSingletonEdgesString ++ nullString) (afterMain ++ "  let s=aG 0 " ++ show numSingleOrig ++ "\n") "  p0 \"\" s"  (numLeaves - numSingle) (numRoots - numSingle) numNetEdges 0
  else if numRoots > (numSingle + 1) then  -- have trivial trees.  one root -> two leaves--always do a tree so always a tree follows if this is done. hence "m++""
    if not doSingle then
      makeGraphProgramString numLeavesOrig numRootsOrig numSingleOrig (beforeMain ++ minimalTreesString ++ nullString) (afterMain ++ "  let m=bG " ++ show (max 0 $ 2*numSingleOrig) ++ " " ++ show numLeaves2 ++ " " ++ show numRoots2 ++ "\n") (sumString ++ "  p0 \"\" (m++") (numLeaves - 2*(numRoots - 1)) 1 numNetEdges 0
    else
      makeGraphProgramString numLeavesOrig numRootsOrig numSingleOrig (beforeMain ++ minimalTreesString ++ nullString) (afterMain ++ "  let m=bG " ++ show (max 0 $ 2*numSingleOrig) ++ " " ++ show numLeaves2 ++ " " ++ show numRoots2 ++ "\n") (sumString ++ "m++") (numLeaves - 2*(numRoots - 1)) 1 numNetEdges 0
  else if numLeaves > 0 then --pectinate tree for remaining leaves one root -> all remining leaves
      if not doEdges then
        if not doSingle && not doMinimal then
          (beforeMain ++ fullTreeString ++ nullString, afterMain ++ "  let t=cG True " ++ show maxVertex ++ " " ++ show numLeaves3 ++ "\n", sumString ++ "  p0 \"\" t")
        else
          (beforeMain ++ fullTreeString, afterMain ++ "  let t=cG True " ++ show maxVertex ++ " " ++ show numLeaves3 ++ "\n", sumString ++ "t)")
      else
        if not doSingle && not doMinimal then
          (beforeMain ++ fullTreeString++ nullString ++ addEdgeString, afterMain ++ "  let t=cG True " ++ show maxVertex ++ " " ++ show numLeaves3 ++ "\n" ++ "  let n=dG "++ show maxVertex2 ++ " t " ++ show numNetEdges ++ "\n", sumString ++ "  p0 \"\" n")
        else
          (beforeMain ++ fullTreeString ++ nullString++ addEdgeString, afterMain ++ "  let t=cG True " ++ show maxVertex ++ " " ++ show numLeaves3 ++ "\n" ++ "  let n=dG " ++ show maxVertex2 ++ " t " ++ show numNetEdges ++ "\n", sumString ++ "n)")
  else (beforeMain, afterMain, sumString)

-- | makeBaseStringGraph creates graph code for a graph that can become a display tree
-- so conditions already checked--leaves > 0, numSunbgle == 0, numRoots == 1 
makeBaseStringGraph :: Int -> String -> String -> String -> Int -> Int -> (String, String, String)
makeBaseStringGraph numLeavesOrig beforeMain afterMain sumString numLeaves numNetEdges =
  let maxVertex2 = (2 * numLeavesOrig) - 1
      zero0 = 0 :: Int -- quiets warning
      necessaryFunctionStrings = beforeMain ++ fullTreeString ++ addEdgeString ++ fmapString ++ sndString ++ headString ++ tailString
                                 ++ elemString ++ notElemString ++ getRepeatedElementsString ++ childrenParentsOfNodeString ++ lastString ++ filterString ++ displayEdgesString ++ nullString
      bitList = replicate numNetEdges (0 :: Int)
  in
  
  (necessaryFunctionStrings, afterMain ++ "  let t=cG True " ++ show zero0 ++ " " ++ show numLeavesOrig ++ "\n" ++ "  let n=dG " ++ show maxVertex2 ++ " t " ++ show numNetEdges ++ "\n" ++ "  let b=" ++ (show bitList) ++ "\n"  ++ "  let d=dE b n\n", sumString ++ "  p0 \"\" d")
  

-- | makeDisplayGraphString cretges code to generate a general gaph and then output a display graph based on that graph
-- the ouput is the display graph
-- if input is tree then return graph string as is
makeDisplayGraphString :: Int -> Int -> Int -> Int -> String
makeDisplayGraphString numLeaves numSingle numRoots numNetEdges =
  if numNetEdges == 0 then makeProgramStringGraph numLeaves numSingle numRoots numNetEdges
  else 
    let (outProgram, middleString, sumString) = makeBaseStringGraph numLeaves programStartStringGraph mainStartString "" numLeaves numNetEdges 
    in 
    outProgram  ++ middleString ++ sumString

{-
Resolve general graph to display tree
(Assuming here can be ie 1 root, leaves >4, net edges > 1, no isolated nodes etc)

0) Get node list from edge list (don't care about root--only once there for sure
  so can just get all n odes that are terminations of edges
1) get indegree > 1 nodes
2) take first netNode
3) get edges and nodes that will be effected
4) delete edges from list that need to be
5) add new edges to list that need to be create
6) since only dealing with edges no need to update nodes
7) recurse till no net nodes
-}

{- -- | Functions that are needed
  need 
    n = list of edges from previous code
    gM = fmapString 
    iM = sndString  
    a5 = headString
    a6 = tailString
    eS = elemString
    nE = notElementString
    eR = getRepeatedElementsString
    cN = childrenParentsOfNode
    a3 = lastString
    g2 = filterString
    dE = displayEdgesString
    nU = nullString
-}

{- general flow
-- get list of nodes that terminate an edge
--  fmap snd [(u,v)]
"let v=gM iM n\n"

--for each member of list, is it repeated--ie found in tail
-- return list of repeated elements

"let r=rE v\n"

-- if r is empty then return current node list for final print and termination
"if null r then p0 \"\" n\n"

-- else resolve (arbitrarily first since only need cost of doing one resolution really) first vertex in list
"else "
  vertex = head r
"let h=a5 r\n"

-- get children and parents of network node from list of edges

-- call with empty lists for child and parent
childrenParentsOfNode :: Int -> [(Int,Int)] -> [Int] -> [Int] -> ([Int], [Int])
childrenParentsOfNode node edgeList childList parentList =
  if null edgeList then (childList, parentList)
  else 
    let (e,u) = head edgeList
        if node == e then childrenParentsOfNode node (tail edgeList) (u:childList) parentList
        else if node == u then childrenParentsOfNode node (tail edgeList) childList (e:parentList)
        else childrenParentsOfNode node (tail edgeList) childList parentList

"let (c,p)=cN h n [] []\n"

-- get grandparents of network edge (need for new edges)
"let (_,pL)=cN (a5 p) n [] []\n"
"    (_,pR)=cn (a3 p) n [] []\n"
--"    gL=a5 pL\n"
"    gR=a5 pR\n"

-- get left and right siblings of net node (h)
--"    (cL, _)=cN (a5 p) n [] []\n"
"    (cR, _)=cN (a3 p) n [] []\n"
--"    sL=a5 $ g2 (/=h) cL\n"
"    sR=a5 $ g2 (/=h) cR\n"


-- get list of edges to delete
  -- connecting child of network node (h) to "left" side
  -- net node to child of net node
  -- parent left to net node
  -- parent right to net node
  -- grand parent right to parent right
  -- parent right to sibling right
"let d=[(h,a5 c),(a5 p, h),(a3 p, h),(gR,a3 p),(a3 p,sR)]\n"

-- get edges to add
"    a=[(a5 p, a5 c),(gR, sR)]\n"

-- filter out nodes to delete and add new ones
--  displayNodes = filter (`notElem` d) n
"   w=(gM(`nE` d)h)++a \n"


-- print edge set and terminate

"   p0 "" d\n"


-- wrap all inside of a recursive function
displayEdges :: [(Int, Int)] -> [(Int, Int)] 
displayEdges n =
  if null n then []
  else 
      -- get node list
      "let v=gM iM n\n"
      -- get repeated node list
      "    r=rE v\n"
      "in\n"
      -- if no repeated nodes then return current list
      "if null r then n\n"
      -- else resolve first net (repeated) node by removing 5 edges and adding 2 new
      "else\n"
      "    let h=a5 r\n"
      -- get children and parents of net node
      "        (c,p)=cN h n [] []\n"
      -- get left and right parent and right grandparent (since adding net child to left parent of net node) 
      --       of network edge (need for new edges)
      "        (_,pL)=cN (a5 p) n [] []\n"
      "        (_,pR)=cn (a3 p) n [] []\n"
      "        gR=a5 pR\n"
      -- get right sibling of net node (h)
      "       (cR, _)=cN (a3 p) n [] []\n"
      "       sR=a5 $ g2 (/=h) cR\n"
      -- get list of edges to delete
      -- connecting child of network node (h) to "left" side
      -- net node to child of net node
      -- parent left to net node
      -- parent right to net node
      -- grand parent right to parent right
      -- parent right to sibling right
      "       d=[(h,a5 c),(a5 p, h),(a3 p, h),(gR,a3 p),(a3 p,sR)]\n"
      -- get edges to add
      "       a=[(a5 p, a5 c),(gR, sR)]\n"
      -- filter out nodes to delete and add new ones and return
      "       w=(gM(`nE` d)h)++a \n"
      "    in displayEdges w\n"
-}