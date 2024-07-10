{-# LANGUAGE Safe #-}

{- |
Exposed constant values for complexity calculations.
-}
module Complexity.Constants
( module Complexity.Constants
) where


-- | epsilon value for floating point comparisons
epsilon ∷ Double
epsilon = 0.0001 ∷ Double


-- | maximum iterations for Newton's method
maxIterations ∷ Int
maxIterations = 100 ∷ Int


-- maximum rate modifier for discrete Gamma
maxGammaRate ∷ Double
maxGammaRate = 10.0 ∷ Double


-- maximum time parameter for exponential distribution
maxTime ∷ Double
maxTime = 2.0 ∷ Double


-- fixed precision for some functins like expE
-- this becuase of precision issues
-- ned this high with facvtorial on Doubles
fixedPrecision ∷ Int
fixedPrecision = 100 ∷ Int


-- strings for use in code generation
maxTimeString ∷ String
maxTimeString = show maxTime


epsilonString ∷ String
epsilonString = show epsilon


maxIterationsString ∷ String
maxIterationsString = show maxIterations


maxGammaRateString ∷ String
maxGammaRateString = show maxGammaRate


fixedPrecisionString ∷ String
fixedPrecisionString = show fixedPrecision
