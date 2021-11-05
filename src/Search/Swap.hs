{- |
Module      :  Swap.hs
Description :  Module specifying graph swapping rearrangement functions
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

module Search.Swap  ( swapMaster
                    ) where

import Types.Types
import qualified ParallelUtilities       as PU
import GeneralUtilities
import qualified Graphs.GraphOperations  as GO
import qualified Utilities.LocalGraph    as LG
import Utilities.Utilities               as U
import Debug.Trace
import           Data.Char
import           Text.Read
import           Data.Maybe

-- | buildArgList is the list of valid build arguments
swapArgList :: [String]
swapArgList = ["spr","tbr", "keep", "steepest", "all"]


-- | swapMaster processes and spawns the swap functions
swapMaster ::  [Argument] -> GlobalSettings -> Int -> [PhylogeneticGraph] -> [PhylogeneticGraph]
swapMaster inArgs inGS rSeed inGraphList =
   trace ("Swapping " ++ (show $ length inGraphList) ++ " graphs") (
   if null inGraphList then []
   else 
      let fstArgList = fmap (fmap toLower . fst) inArgs
          sndArgList = fmap (fmap toLower . snd) inArgs
          lcArgList = zip fstArgList sndArgList
          checkCommandList = U.checkCommandArgs "swap" fstArgList swapArgList
   in
   -- check for valid command options
   if not checkCommandList then errorWithoutStackTrace ("Unrecognized command in 'swap': " ++ show inArgs)
   else 
       let keepList = filter ((=="keep").fst) lcArgList
           keepNum
            | length keepList > 1 =
              errorWithoutStackTrace ("Multiple 'keep' number specifications in swap command--can have only one: " ++ show inArgs)
            | null keepList = Just 1
            | otherwise = readMaybe (snd $ head keepList) :: Maybe Int
      in
      if isNothing keepNum then errorWithoutStackTrace ("Keep specification not an integer: "  ++ show (snd $ head keepList))
      else 
         let doSPR' = any ((=="spr").fst) lcArgList
             doTBR = any ((=="tbr").fst) lcArgList
             doSteepest' = any ((=="steepest").fst) lcArgList
             doAll = any ((=="all").fst) lcArgList
             doSPR = if (not doSPR' && not doTBR) then True
                     else doSPR'
             doSteepest = if (not doSteepest' && not doAll) then True
                          else doSteepest'
             newGraphList  = if doSPR then concatMap (swapSPR inGS (fromJust keepNum) doSteepest) inGraphList
                             else inGraphList
             newGraphList' = if doTBR then concatMap (swapTBR inGS (fromJust keepNum) doSteepest) newGraphList
                             else newGraphList
            in
            newGraphList'
   )


-- | swapSPR perfomrs SPR branch (edge) swapping on graphs
swapSPR :: GlobalSettings -> Int -> Bool -> PhylogeneticGraph -> [PhylogeneticGraph]
swapSPR inGS numToKeep steepest inGraph = 
   if LG.isEmpty (fst6 inGraph) then []
   else 
      [inGraph]

-- | swapTBR performs TBR branch (edge) swapping on graphs
swapTBR :: GlobalSettings -> Int -> Bool -> PhylogeneticGraph -> [PhylogeneticGraph]
swapTBR inGS numToKeep steepest inGraph = 
   if LG.isEmpty (fst6 inGraph) then []
   else 
      [inGraph]