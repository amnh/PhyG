{- |
Module      :  Debug.hs
Description :  Module with Debug version of functions 
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



module Debug.Debug where

import           Types.Types
-- import Data.List
import qualified Data.Vector as V

debugZip :: [a] -> [b] -> [(a,b)]
debugZip la lb = 
  if not isDebug then zip la lb
  else
     if length la /= length lb then error ("Zip arguments not equal in length: " ++ show (length la, length lb))
     else if null la then error ("First list null in debugZip " ++ show (length la, length lb))
     else if null lb then error ("Second list null in debugZip " ++ show (length la, length lb))  
     else zip la lb 
    

debugZip3 :: [a] -> [b] -> [c] -> [(a,b,c)]
debugZip3 la lb lc = 
  if not isDebug then zip3 la lb lc
  else
     if (length la /= length lb) || (length la /= length lc) || (length lb /= length lc) then error ("Zip3 arguments not equal in length: " ++ show (length la, length lb, length lc))
     else if null la then error ("First list null in debugZip3 " ++ show (length la, length lb, length lc))
     else if null lb then error ("Second list null in debugZip3 " ++ show (length la, length lb, length lc))  
     else if null lc then error ("Third list null in debugZip3 " ++ show (length la, length lb, length lc))  
     else zip3 la lb  lc

debugVectorZip :: V.Vector a -> V.Vector b -> V.Vector (a,b)
debugVectorZip la lb = 
  if not isDebug then V.zip la lb
  else
     if V.length la /= V.length lb then error ("Zip arguments not equal in length: " ++ show (V.length la, V.length lb))
     else if V.null la then error ("First vector null in debugZip " ++ show (V.length la, V.length lb))
     else if V.null lb then error ("Second vector null in debugZip " ++ show (V.length la, V.length lb))  
     else V.zip la lb 
    

debugVectorZip3 :: V.Vector a -> V.Vector b -> V.Vector c -> V.Vector (a,b,c)
debugVectorZip3 la lb lc = 
  if not isDebug then V.zip3 la lb lc
  else
     if (V.length la /= V.length lb) || (V.length la /= V.length lc) || (V.length lb /= V.length lc) then error ("Zip3 arguments not equal in length: " ++ show (V.length la, V.length lb, V.length lc))
     else if V.null la then error ("First vector null in debugZip3 " ++ show (V.length la, V.length lb, V.length lc))
     else if V.null lb then error ("Second vector null in debugZip3 " ++ show (V.length la, V.length lb, V.length lc))  
     else if V.null lc then error ("Third vector null in debugZip3 " ++ show (V.length la, V.length lb, V.length lc))  
     else V.zip3 la lb  lc

debugVectorZip4 :: V.Vector a -> V.Vector b -> V.Vector c -> V.Vector d -> V.Vector (a,b,c,d)
debugVectorZip4 la lb lc ld = 
  if not isDebug then V.zip4 la lb lc ld
  else
     if (V.length la /= V.length lb) || (V.length la /= V.length lc) || (V.length la /= V.length ld) || (V.length lb /= V.length lc) || (V.length lb /= V.length ld) || (V.length lc /= V.length ld) then error ("Zip3 arguments not equal in length: " ++ show (V.length la, V.length lb, V.length lc, V.length ld))
     else if V.null la then error ("First vector null in debugZip4 " ++ show (V.length la, V.length lb, V.length lc, V.length ld))
     else if V.null lb then error ("Second vector null in debugZip4 " ++ show (V.length la, V.length lb, V.length lc, V.length ld))  
     else if V.null lc then error ("Third vector null in debugZip4 " ++ show (V.length la, V.length lb, V.length lc, V.length ld))  
     else if V.null ld then error ("Fourth vector null in debugZip4 " ++ show (V.length la, V.length lb, V.length lc, V.length ld))  
     else V.zip4 la lb  lc ld

