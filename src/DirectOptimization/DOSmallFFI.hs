module DirectOptimization.DOSmallFFI (wrapperPCG_DO_Small_FFI) where

import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise
import Bio.Character.Encodable  -- this has DynamicCharacter reference
import Data.Alphabet
import Data.Foldable
import Data.TCM.Dense
--import qualified Data.BitVector as BV
import qualified Data.BitVector.LittleEndian as BV
import Data.List.NonEmpty (NonEmpty(..))
import qualified Data.List.NonEmpty as NE
import Data.Vector (Vector, (!))
import qualified Data.Vector as V
import qualified SymMatrix as SM
import GeneralUtilities
import Utilities.Utilities
import Debug.Trace



wrapperPCG_DO_Small_FFI :: Vector BV.BitVector -> Vector BV.BitVector -> SM.Matrix Int -> DenseTransitionCostMatrix -> (Vector BV.BitVector, Int)
wrapperPCG_DO_Small_FFI lhs rhs tcmVect tcmDense = (resultingMedians, fromEnum resultCost)
    where
        --(resultCost, resultFFI) = foreignPairwiseDO tcmDense lhsDC rhsDC
        (resultCost, resultFFI) = foreignPairwiseDO tcmDense lhsDC rhsDC

        bitStreams = decodeStream arbitraryAlphabet resultFFI

        resultingMedians = V.fromList . toList $ fmap (BV.fromBits . g 0 . toList) bitStreams
        
        
        lhsDC = buildDC lhs 
        rhsDC = buildDC rhs

        buildDC :: Vector BV.BitVector -> DynamicCharacter
        -- buildDC = encodeStream specializedAlphabetToDNA . fmap (NE.fromList . f 0 . BV.toBits) . NE.fromList  . toList
        buildDC = encodeStream arbitraryAlphabet . fmap (NE.fromList . f 0 . BV.toBits) . NE.fromList  . toList

        f :: Word -> [Bool] -> [String]
        f _ [] = []
        f n (x:xs)
            | x = show n : f (n+1) xs
            | otherwise = f (n+1) xs
                
        -- modified to allow for variable size alphabets
        g :: Word -> [String] -> [Bool]
        g n thang
          | fromEnum n >= length tcmVect = []
          | null thang = False : g (n+1) []
          | read (head thang) == n = True  : g (n+1) (tail thang)
          | otherwise   = False : g (n+1) thang
        
        {-
        tcmDense = generateDenseTransitionCostMatrix 0 5 getCost

        getCost i j = 
            let x = SM.getFullVects tcm
            in  toEnum $ (x ! fromEnum i) ! fromEnum j
        -}
        arbitraryAlphabet :: Alphabet String
        arbitraryAlphabet = fromSymbols $ show <$> 0 :| [1 .. length tcmVect - 1]



{-
wrapperPCG_DO_Small_FFI :: Vector BV.BitVector -> Vector BV.BitVector -> SM.Matrix Int 
                         -> DenseTransitionCostMatrix -> (Vector BV.BitVector, Int)
wrapperPCG_DO_Small_FFI lhs rhs tcmVect tcmDense = 
    --trace ("FFI: \n" ++ (show $ snd4 $ convertDynamicCharacter resultFFI)) (resultingMedians, fromEnum resultCost)
    (snd4 $ dynamicCharacterTo3Vector resultFFI, fromEnum resultCost)
    where
        --(resultCost, resultFFI) = foreignPairwiseDO tcmDense lhsDC rhsDC
        (resultCost, resultFFI) = foreignPairwiseDO tcmDense lhsDC rhsDC

        -- bitStreams = decodeStream specializedAlphabetToDNA resultFFI
        bitStreams = decodeStream arbitraryAlphabet resultFFI

        resultingMedians = V.fromList . toList $ fmap (BV.fromBits . g 0 . toList) bitStreams
        --resultingMedians = V.fromList . toList $ fmap (BV.fromBits . g' (length tcmVect) 0 . toList) bitStreams
        
        {-
        lhsDC = buildDC lhs 
        rhsDC = buildDC rhs
        -}
        lhsDC = convertVectorToDynamicCharacter lhs
        rhsDC = convertVectorToDynamicCharacter rhs

        buildDC :: Vector BV.BitVector -> DynamicCharacter
        -- buildDC = encodeStream specializedAlphabetToDNA . fmap (NE.fromList . f 0 . BV.toBits) . NE.fromList  . toList
        buildDC = encodeStream arbitraryAlphabet . fmap (NE.fromList . f 0 . BV.toBits) . NE.fromList  . toList

        f :: Word -> [Bool] -> [String]
        f _ [] = []
        f n (x:xs)
            | x = show n : f (n+1) xs
            | otherwise = f (n+1) xs
                
        -- modified to allow for variable size alphabets
        g :: Word -> [String] -> [Bool]
        g n thang
          | fromEnum n >= length tcmVect = []
          | null thang = False : g (n+1) []
          | read (head thang) == n = True  : g (n+1) (tail thang)
          | otherwise   = False : g (n+1) thang
        
        {- never halts due to 1st line
        g :: Word -> [String] -> [Bool]
        g n [] = False : g (n+1) []
        g n (x:xs)
          -- | fromEnum n >= TCM.size tcm = []
          | fromEnum n >= length tcmVect = []
          | read x == n = True  : g (n+1)    xs
          | otherwise   = False : g (n+1) (x:xs)
        
        -}
        {-
        g :: Word -> [String] -> [Bool]
        g 5 _ = []
        g n [] = False : g (n+1) []
        g n (x:xs)
            | read x == n = True  : g (n+1)    xs
            | otherwise   = False : g (n+1) (x:xs)
        -}        
                
        arbitraryAlphabet :: Alphabet String
        arbitraryAlphabet = fromSymbols $ show <$> 0 :| [1 .. length tcmVect - 1]
-}
{-
g' :: Int -> Int -> [String] -> [Bool]
g' alphSize n inList = 
    if n == alphSize then []
    else if null inList then g' alphSize (n + 1) []
    else 
        let fOne = head inList
        in
        if read fOne == n then True : g' alphSize (n+1) (tail inList)
        else False : g' alphSize (n+1) inList
-}


       

    {-
    if V.null inVector then Missing (0 :: Word)
    else 
        let neV = NE.fromList $ toList $ V.zip3 inVector inVector inVector
        in
        DC neV
    -}