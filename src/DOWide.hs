module DOWide where

import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise
import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Wide
import Bio.Character.Encodable
import Data.Alphabet
import Data.Bits
import Data.Foldable
import           Data.BitVector.LittleEndian (BitVector)
import qualified Data.BitVector.LittleEndian as BV
import           Data.List.NonEmpty  (NonEmpty(..))
import qualified Data.List.NonEmpty  as NE
import           Data.MetricRepresentation
import           Data.Vector         (Vector, (!))
import qualified Data.Vector         as  V
import qualified Data.Vector.Unboxed as UV
import           Data.Word
import qualified SymMatrix           as SM


wrapperSlimDO
  :: Vector BitVector
  -> Vector BitVector
  -> MetricRepresentation (Word64 -> Word64 -> (Word64, Word))
  -> (Vector BitVector, Int)
-- wrapperPCG_DO_FFI lhs rhs tcm | trace (show tcm) False= undefined
wrapperSlimDO lhs rhs metric = (resultMedians, fromEnum resultCost)
  where
    (resultCost, resultMedians) = wideDC2BVs (fromEnum n) <$> widePairwiseDO (fromIntegral n) tcm lhsDC rhsDC

    tcm = retreivePairwiseTCM metric

    n = case length lhs of
          0 -> case length rhs of
                 0 -> 64
                 _ -> BV.dimension $ V.head rhs
          _ -> BV.dimension $ V.head lhs

    lhsDC = bvs2WideDC lhs 
    rhsDC = bvs2WideDC rhs

    
{-    
    getCost i j = 
         let x = SM.getFullVects tcm
         in  toEnum $ (x ! fromEnum i) ! fromEnum j
-}

--specializedAlphabetToDNA :: Alphabet String
--specializedAlphabetToDNA = fromSymbols $ show <$> (0 :: Word) :| [1 .. 4]


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
        in  foldr f (BV.fromNumber (toEnum n) 0) [0 .. n - 1]
