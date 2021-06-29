module DOSlim where

import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise
import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Slim
import Bio.Character.Encodable
import Data.Alphabet
import Data.Bits
import Data.Foldable
import Data.TCM.Dense
import           Data.BitVector.LittleEndian (BitVector)
import qualified Data.BitVector.LittleEndian as BV
import Data.List.NonEmpty (NonEmpty(..))
import qualified Data.List.NonEmpty as NE
import           Data.Vector (Vector, (!))
import qualified Data.Vector  as V
import qualified Data.Vector.Storable as SV
import qualified Data.Vector as V
import           Data.Word
import           Foreign.C.Types (CUInt)
import qualified SymMatrix as SM

import Debug.Trace


wrapperSlimDO :: Vector BitVector -> Vector BitVector -> Vector (Vector Int) -> (Vector BitVector, Int)
-- wrapperPCG_DO_FFI lhs rhs tcm | trace (show tcm) False= undefined
wrapperSlimDO lhs rhs tcm = (resultMedians, fromEnum resultCost)
  where
    (resultCost, resultMedians) = slimDC2BVs 5 <$> smallAlphabetPairwiseDO tcmDense lhsDC rhsDC

    lhsDC = bvs2SlimDC lhs 
    rhsDC = bvs2SlimDC rhs
    tcmDense = generateDenseTransitionCostMatrix 0 5 getCost
    getCost i j = 
         let x = SM.getFullVects tcm
         in  toEnum $ (x ! fromEnum i) ! fromEnum j


--specializedAlphabetToDNA :: Alphabet String
--specializedAlphabetToDNA = fromSymbols $ show <$> (0 :: Word) :| [1 .. 4]


bvs2SlimDC :: V.Vector BitVector -> SlimDynamicCharacter
bvs2SlimDC v = (x,x,x)
  where
    x = SV.generate (V.length v) $ \i -> bv2w (v ! i)
    
    bv2w :: BitVector -> CUInt
    bv2w bv =
        let f i a
              | bv `testBit` i = a `setBit` i
              | otherwise      = a
        in  foldr f 0 [0 .. fromEnum $ BV.dimension bv - 1]


slimDC2BVs :: Int -> SlimDynamicCharacter -> V.Vector BitVector
slimDC2BVs n (x,_,_) = V.generate (SV.length x) $ \i -> w2bv (x SV.! i)
  where
    w2bv :: CUInt -> BitVector
    w2bv w =
        let f i a
              | w `testBit` i = a `setBit` i
              | otherwise     = a
        in  foldr f (BV.fromNumber (toEnum n) 0) [0 .. n - 1]
