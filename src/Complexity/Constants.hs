{- |
Module      :  Constants
Description :  Constant values used by functions
Copyright   :  (c) 2019-2020 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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

module Complexity.Constants where

-- | epsilon value for floating point comparisons
epsilon :: Double
epsilon = 0.0001 :: Double

-- | maximum iterations for Newton's method
maxIterations :: Int
maxIterations = 100 :: Int

--maximum rate modifier for discrete Gamma
maxGammaRate :: Double
maxGammaRate = 10.0 :: Double

-- maximum time parameter for exponential distribution
maxTime :: Double
maxTime = 2.0 :: Double

-- fixed precision for some functins like expE
-- this becuase of precision issues
-- ned this high with facvtorial on Doubles
fixedPrecision :: Int
fixedPrecision = 100 :: Int

-- strings for use in code generation
maxTimeString :: String
maxTimeString= show maxTime

epsilonString :: String
epsilonString = show epsilon

maxIterationsString :: String
maxIterationsString = show maxIterations

maxGammaRateString :: String
maxGammaRateString = show maxGammaRate

fixedPrecisionString :: String
fixedPrecisionString = show fixedPrecision
