{- |
Module      :  Medians.hs
Description :  Module specifying data type medians
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

module Medians (  median2
                ) where

import           Types
import           Debug.Trace
import qualified Data.Vector as V
import qualified Data.BitVector as BV
import           GeneralUtilities

-- | median2 takes character data and returns median character and cost
median2 :: CharacterData -> CharacterData -> CharInfo -> (CharacterData, VertexCost)
median2 firstVertChar secondVertChar inCharInfo = 
    let thisType = charType inCharInfo
        thisWeight = weight inCharInfo
    in
    if thisType == Add then 
        let newCharVect = intervalAdd thisWeight firstVertChar secondVertChar 
        in
        (newCharVect, localCost  newCharVect)

    else if thisType == NonAdd then 
        let newCharVect = interUnion thisWeight firstVertChar secondVertChar 
        in
        (newCharVect, localCost  newCharVect)

    else if thisType == Matrix then (firstVertChar, 0.0)

    else if thisType `elem` [SmallAlphSeq, NucSeq] then (firstVertChar, 0.0)

    else if thisType `elem` [AminoSeq, GenSeq] then (firstVertChar, 0.0)

    else error ("Character type " ++ show thisType ++ " unrecongized/not implemented")

-- | localOr wrapper for BV.or for vector elements
localOr :: BV.BV -> BV.BV -> BV.BV
localOr lBV rBV = BV.or [lBV, rBV]

-- | localAnd wrapper for BV.and for vector elements
localAnd :: BV.BV -> BV.BV -> BV.BV
localAnd lBV rBV = BV.and [lBV, rBV]

-- | localAndOr takes the intesection vect and union vect elements
-- and return intersection is /= 0 otherwise union
localAndOr :: BV.BV -> BV.BV -> BV.BV -> BV.BV
localAndOr localZero interBV unionBV = if (interBV BV.==. localZero) then unionBV else interBV


-- | interUnion takes two non-additive chars and creates newCharcter as 2-median
-- in post-order pass to create preliminary states assignment
-- assumes a single weight for all
-- performs two passes though chars to get cost of assignments
interUnion :: Double -> CharacterData -> CharacterData -> CharacterData
interUnion thisWeight leftChar rightChar =
    let intersectVect =  V.zipWith localAnd (stateBVPrelim leftChar) (stateBVPrelim rightChar)
        unionVect = V.zipWith localOr (stateBVPrelim leftChar) (stateBVPrelim rightChar)
        localZero = BV.bitVec (BV.size $ V.head $ stateBVPrelim leftChar) (0 :: Integer)
        numUnions = V.length $ V.filter (BV.==. localZero) intersectVect
        newCost = thisWeight * (fromIntegral numUnions)
        newStateVect = V.zipWith (localAndOr localZero) intersectVect unionVect
        newCharcater = CharacterData {  stateBVPrelim = newStateVect
                                      , minRangePrelim = V.empty
                                      , maxRangePrelim = V.empty
                                      , matrixStatesPrelim = V.empty
                                      , stateBVFinal = V.singleton BV.nil
                                      , minRangeFinal = V.empty
                                      , maxRangeFinal = V.empty
                                      , matrixStatesFinal = V.empty
                                      , approxMatrixCost = V.singleton 0
                                      , localCostVect = V.singleton 0
                                      , localCost = newCost
                                      , globalCost = newCost + (globalCost leftChar) + (globalCost rightChar) 
                                      }
    in 
    trace ("NonAdditive: " ++ (show numUnions) ++ " " ++ (show newCost) ++ "\n\t" ++ (show $ stateBVPrelim leftChar) ++ "\n\t" ++ (show $ stateBVPrelim rightChar) ++ "\n\t"
        ++ (show intersectVect) ++ "\n\t" ++ (show unionVect) ++ "\n\t" ++ (show newStateVect))
    newCharcater

-- | getNewRange takes min and max range of two additive charcaters and returns 
-- a triple of (newMin, newMax, Cost)
getNewRange :: (Int, Int, Int, Int) -> (Int, Int, Int)
getNewRange inStuff@(lMin, lMax, rMin, rMax) = 
    -- subset
    if (rMin >= lMin) && (rMax <= lMax) then (rMin, rMax, 0)
    else if  (lMin >= rMin) && (lMax <= rMax) then (lMin, lMax, 0)
    -- overlaps
    else if (rMin >= lMin) && (rMax >= lMax) then (rMin, lMax,0)
    else if (lMin >= rMin) && (lMax >= rMax) then (lMin, rMax,0)
    -- newInterval
    else if (lMax <= rMin) then (lMax, rMin ,rMin - lMax)
    else if (rMax <= lMin) then (rMax, lMin, lMin - rMax)
    else error ("This can't happen " ++ show inStuff)

-- | intervalAdd takes two additive chars and creates newCharcter as 2-median
-- in post-order pass to create preliminary states assignment
-- assumes a single weight for all
intervalAdd :: Double -> CharacterData -> CharacterData -> CharacterData
intervalAdd thisWeight leftChar rightChar =
    let newRangeCosts = V.map getNewRange $ V.zip4 (minRangePrelim leftChar) (maxRangePrelim leftChar) (minRangePrelim rightChar) (maxRangePrelim rightChar)
        newMinRange = V.map fst3 newRangeCosts
        newMaxRange = V.map snd3 newRangeCosts
        newCost = thisWeight * (fromIntegral $ V.sum $ V.map thd3 newRangeCosts)
        newCharcater = CharacterData {  stateBVPrelim = V.empty
                                      , minRangePrelim = newMinRange
                                      , maxRangePrelim = newMaxRange
                                      , matrixStatesPrelim = V.empty
                                      , stateBVFinal = V.singleton BV.nil
                                      , minRangeFinal = V.empty
                                      , maxRangeFinal = V.empty
                                      , matrixStatesFinal = V.empty
                                      , approxMatrixCost = V.singleton 0
                                      , localCostVect = V.singleton 0
                                      , localCost = newCost
                                      , globalCost = newCost + (globalCost leftChar) + (globalCost rightChar) 
                                      }
    in 
    trace ("Additive: " ++ (show newCost) ++ "\n\t" ++ (show $ minRangePrelim leftChar) ++ "\n\t" ++ (show $ maxRangePrelim leftChar) ++ "\n\t"
        ++ (show $ minRangePrelim rightChar) ++ "\n\t" ++ (show $ maxRangePrelim rightChar) ++ "\n\t"  
        ++ (show newMinRange) ++ "\n\t" ++ (show newMaxRange) ++ "\n\t" 
        ++ (show newRangeCosts))
    newCharcater
