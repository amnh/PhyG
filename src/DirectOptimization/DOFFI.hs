module DirectOptimization.DOFFI where

import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise
import Bio.Character.Encodable
import Data.Alphabet
import Data.Foldable
import Data.TCM.Dense
import qualified Data.BitVector as BV
import Data.List.NonEmpty (NonEmpty(..))
import qualified Data.List.NonEmpty as NE
import Data.Vector (Vector, (!))
import qualified Data.Vector as V
import qualified SymMatrix as SM

import Debug.Trace


wrapperPCG_DO_FFI :: Vector BV.BV -> Vector BV.BV -> Vector (Vector Int) -> (Vector BV.BV, Int)
-- wrapperPCG_DO_FFI lhs rhs tcm | trace (show tcm) False= undefined
wrapperPCG_DO_FFI lhs rhs tcm = (resultingMedians, fromEnum resultCost)
    where
        (resultCost, resultFFI) = foreignPairwiseDO tcmDense lhsDC rhsDC

        bitStreams = decodeStream specializedAlphabetToDNA resultFFI

        resultingMedians = V.fromList . toList $ fmap (BV.fromBits . g 0 . toList) bitStreams

        lhsDC = buildDC lhs 
        rhsDC = buildDC rhs

        tcmDense = generateDenseTransitionCostMatrix 0 5 getCost

        getCost i j = 
            let x = SM.getFullVects tcm
            in  toEnum $ (x ! fromEnum i) ! fromEnum j

        buildDC :: Vector BV.BV -> DynamicCharacter
        buildDC = encodeStream specializedAlphabetToDNA . fmap (NE.fromList . f 0 . BV.toBits) . NE.fromList  . toList

        f :: Word -> [Bool] -> [String]
        f _ [] = []
        f n (x:xs)
          | x = show n : f (n+1) xs
          | otherwise = f (n+1) xs

        g :: Word -> [String] -> [Bool]
        g 5 _ = []
        g n [] = False : g (n+1) []
        g n (x:xs)
          | read x == n = True  : g (n+1)    xs
          | otherwise   = False : g (n+1) (x:xs)


specializedAlphabetToDNA :: Alphabet String
specializedAlphabetToDNA = fromSymbols $ show <$> (0 :: Word) :| [1 .. 4]