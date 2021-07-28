module DirectOptimization.DOWide where

import           Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Wide
import           Control.Arrow (first)
import           Data.BitVector.LittleEndian                                 (BitVector)
import qualified Data.BitVector.LittleEndian                                 as BV
import           Data.Bits
import           Data.MetricRepresentation
import           Data.Vector                                                 (Vector, (!))
import qualified Data.Vector                                                 as V
import qualified Data.Vector.Unboxed                                         as UV
import           Data.Word


wrapperWideDO
  :: Vector BitVector
  -> Vector BitVector
  -> MetricRepresentation BV.BitVector  --Word64
  -> (Vector BitVector, Int)
wrapperWideDO lhs rhs metric = (wideDC2BVs (fromIntegral n) resultMedians, fromEnum resultCost)
  where
    (resultCost, resultMedians) = widePairwiseDO (fromIntegral n) tcm lhsDC rhsDC

    tcm :: Word64 -> Word64 -> (Word64, Word)
    tcm x y = first BV.toUnsignedNumber $ (retreivePairwiseTCM metric) (BV.fromNumber 64 x) (BV.fromNumber 64 y)

    n = case length lhs of
          0 -> case length rhs of
                 0 -> 64
                 _ -> BV.dimension $ V.head rhs
          _ -> BV.dimension $ V.head lhs

    lhsDC = bvs2WideDC lhs
    rhsDC = bvs2WideDC rhs

bvs2WideDC :: V.Vector BitVector -> WideDynamicCharacter
bvs2WideDC v = (x,x,x)
  where
    x = UV.generate (V.length v) $ \i -> bv2w (v ! i)

    bv2w :: BitVector -> Word64
    bv2w bv =
        let f i a
              | bv `testBit` i = a `setBit` i
              | otherwise      = a
        in  foldr f 0 [0 .. fromEnum $ BV.dimension bv - 1]


wideDC2BVs :: Int -> WideDynamicCharacter -> V.Vector BitVector
wideDC2BVs n (x,_,_) = V.generate (UV.length x) $ \i -> w2bv (x UV.! i)
  where
    w2bv :: Word64 -> BitVector
    w2bv w =
        let f i a
              | w `testBit` i = a `setBit` i
              | otherwise     = a
        in  foldr f (BV.fromNumber (toEnum n) (0 :: Word)) [0 .. n - 1]
