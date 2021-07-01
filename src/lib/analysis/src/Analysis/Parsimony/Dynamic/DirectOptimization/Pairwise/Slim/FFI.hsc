-----------------------------------------------------------------------------
-- |
-- an FFI interface for C code that efficiently aligns either two or three
-- sequences, using Ukkonen when appropriate, in both affine and non-affine
-- cases.
--
-- This example uses pointers, both to structs and to fields within the
-- structs. This is much easier to accomplish via .hsc rather than doing
-- straight FFI. A .hsc file are read by hsc2hs, which then creates a .c
-- file, which is compiled and run to create an .hs file, which is then
-- compiled for use in outside modules.
--
-- For notes on usage, data construction and external see referenced C
-- compilation units, and also driver.c, which is not imported, but is
-- included indirectory for reference.
--
-----------------------------------------------------------------------------

-- TODO: do I need this: https://hackage.haskell.org/package/base-4.9.0.0/docs/Foreign-StablePtr.html
{-# LANGUAGE BangPatterns             #-}
{-# LANGUAGE DeriveGeneric            #-}
{-# LANGUAGE FlexibleContexts         #-}
{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE StrictData               #-}

module Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Slim.FFI
  ( DenseTransitionCostMatrix
  , smallAlphabetPairwiseDO
--  , foreignThreeWayDO
  ) where

import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.DynamicCharacter2
import Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.Slim.Internal
import Data.Coerce
import Data.MonoTraversable
import Data.TCM.Dense
import           Data.Vector.Storable (Vector)
import qualified Data.Vector.Storable as V
import Foreign
import Foreign.C.Types
import GHC.ForeignPtr
import Prelude   hiding (sequence, tail)
import System.IO.Unsafe (unsafePerformIO)

#include "c_alignment_interface2.h"
#include "c_code_alloc_setup.h"
#include "costMatrix.h"
#include "alignmentMatrices.h"


-- |
-- Specify whether or not to compute median state values
data MedianContext = ComputeMedians | DoNotComputeMedians


-- |
-- Specify whether or not to compute union state values
data UnionContext  = ComputeUnions  | DoNotComputeUnions


instance Enum MedianContext where

    fromEnum      ComputeMedians = 1
    fromEnum DoNotComputeMedians = 0

    toEnum 0 = DoNotComputeMedians
    toEnum _ =      ComputeMedians


instance Enum UnionContext where

    fromEnum      ComputeUnions = 1
    fromEnum DoNotComputeUnions = 0

    toEnum 0 = DoNotComputeUnions
    toEnum _ =      ComputeUnions


foreign import ccall unsafe "c_alignment_interface2.h cAlign2D"

    align2dFn_c
      :: Ptr CUInt -- ^ character1, input & output (lesser)
      -> Ptr CUInt -- ^ character2, input & output (longer)
      -> Ptr CUInt -- ^ gapped median output
      -> Ptr CSize -- ^ length median output
      -> CSize     -- ^ size of each buffer
      -> CSize     -- ^ length of character1
      -> CSize     -- ^ length of character2
      -> Ptr CostMatrix2d
      -> CInt      -- ^ compute ungapped & not   gapped medians
      -> CInt      -- ^ compute   gapped & not ungapped medians
      -> CInt      -- ^ compute union
      -> CSize     -- ^ cost


{-
foreign import ccall unsafe "c_alignment_interface2.h cAlignAffine2D"

    align2dAffineFn_c
      :: Ptr CUInt -- ^ character1, input & output (lesser)
      -> Ptr CUInt -- ^ character2, input & output (longer)
      -> Ptr CUInt -- ^ gapped median output
      -> Ptr CSize -- ^ length median output
      -> CSize     -- ^ size of each buffer
      -> CSize     -- ^ length of character1
      -> CSize     -- ^ length of character2
      -> Ptr CostMatrix2d
      -> CInt      -- ^ compute medians
      -> CSize     -- ^ cost
-}


{-
-- | Create and allocate cost matrix
-- first argument, TCM, is only for non-ambiguous nucleotides, and it used to generate
-- the entire cost matrix, which includes ambiguous elements.
-- TCM is row-major, with each row being the left character element.
-- It is therefore indexed not by powers of two, but by cardinal integer.
foreign import ccall unsafe "c_alignment_interface.h align3d"

    align3dFn_c :: Ptr Align_io -- ^ character1, input
      -> Ptr Align_io -- ^ character2, input
      -> Ptr Align_io -- ^ character3, input
      -> Ptr Align_io -- ^ character1, output
      -> Ptr Align_io -- ^ character2, output
      -> Ptr Align_io -- ^ character3, output
      -> Ptr Align_io -- ^ gapped median output
      -> Ptr Align_io -- ^ ungapped median output
      -> Ptr CostMatrix3d
      -> CInt        -- ^ substitution cost
      -> CInt        -- ^ gap open cost
      -> CInt        -- ^ indel cost
      -> CInt        -- ^ alignment cost
-}


-- |
-- Align two dynamic characters using an FFI call for more efficient computation
-- on small alphabet sizes.
--
-- Requires a pre-generated 'DenseTransitionCostMatrix' from a call to
-- 'generateDenseTransitionCostMatrix' defining the alphabet and transition costs.
-- {-# INLINE foreignPairwiseDO #-}
{-# SCC smallAlphabetPairwiseDO #-}
smallAlphabetPairwiseDO
  :: DenseTransitionCostMatrix    -- ^ Structure defining the transition costs between character states
  -> SlimDynamicCharacter         -- ^ First  dynamic character
  -> SlimDynamicCharacter         -- ^ Second dynamic character
  -> (Word, SlimDynamicCharacter) -- ^ The /ungapped/ character derived from the the input characters' N-W-esque matrix traceback
smallAlphabetPairwiseDO = algn2d DoNotComputeUnions ComputeMedians


{-
-- |
-- Align three dynamic characters using an FFI call for more efficient computation
-- on small (or smallish) alphabet sizes.
--
-- Requires a pre-generated 'DenseTransitionCostMatrix' from a call to
-- 'generateDenseTransitionCostMatrix' defining the alphabet and transition costs.
foreignThreeWayDO :: ( EncodableDynamicCharacter s
                     , ExportableElements s
--                     , Show s
                     )
                  => s                         -- ^ First  dynamic character
        -> s                         -- ^ Second dynamic character
        -> s                         -- ^ Third  dynamic character
        -> Int                       -- ^ Mismatch cost
        -> Int                       -- ^ Gap open cost
        -> Int                       -- ^ Indel cost
        -> DenseTransitionCostMatrix -- ^ Structure defining the transition costs between character states
        -> (Word, s, s, s, s, s)     -- ^ The /ungapped/ character derived from the the input characters' N-W-esque matrix traceback
foreignThreeWayDO char1 char2 char3 costMatrix = algn3d char1 char2 char3 costMatrix
-}


-- |
-- Performs a naive direct optimization
-- Takes in two characters to run DO on and a metadata object
-- Returns an assignment character, the cost of that assignment, the assignment character with gaps included,
-- the aligned version of the first input character, and the aligned version of the second input character
-- The process for this algorithm is to generate a traversal matrix then perform a traceback.
-- {-# INLINE algn2d #-}
algn2d
  :: UnionContext
  -> MedianContext
  -> DenseTransitionCostMatrix    -- ^ Structure defining the transition costs between character states
  -> SlimDynamicCharacter         -- ^ First  dynamic character
  -> SlimDynamicCharacter         -- ^ Second dynamic character
  -> (Word, SlimDynamicCharacter) -- ^ The cost of the alignment
algn2d computeUnion computeMedians denseTCMs lhs rhs =
{-
    let (gapsChar1, ungappedChar1) = (\v@(y,x) -> trace ("CHAR 1: " <> show x <> "\ngaps: " <> show y) v) $ deleteGaps char1
        (gapsChar2, ungappedChar2) = (\v@(y,x) -> trace ("CHAR 2: " <> show x <> "\ngaps: " <> show y) v) $ deleteGaps char2
        (swapped, shorterChar, longerChar) = (\v@(s,_,_) -> trace ("SWAPPED: " <> show s) v) $ measureCharacters ungappedChar1 ungappedChar2
-}
    let (swapped, gapsLesser, gapsLonger, (shorterChar,_,_), (longerChar,_,_)) = measureAndUngapCharacters (coerceEnum symbolCount) lhs rhs
        (alignmentCost, ungappedAlignment) =
          if      olength shorterChar == 0
          then if olength  longerChar == 0
               -- Niether character was Missing, but both are empty when gaps are removed
               then (0, (mempty,mempty,mempty))
               -- Niether character was Missing, but one of them is empty when gaps are removed
               else let gap = bit $ fromEnum symbolCount :: CUInt
                        vec = V.generate (V.length longerChar) $ \i ->
                                  fst (lookupPairwise denseTCMs (longerChar V.! i) gap)
                    in  (0, (vec, V.replicate (V.length longerChar) 0, longerChar))
               -- Both have some non-gap elements, perform string alignment
          else f shorterChar longerChar
        transformation = if swapped then \(m,l,r) -> (m,r,l) else id
{-
        (gapsLesser, gapsLonger, transformation)
          | swapped   = (gapsChar2, gapsChar1, omap swapContext)
          | otherwise = (gapsChar1, gapsChar2, id)
-}
        regappedAlignment = insertGaps (coerceEnum symbolCount) gapsLesser gapsLonger ungappedAlignment
        alignmentContext  :: SlimDynamicCharacter 
        alignmentContext  = transformation regappedAlignment
    in  (alignmentCost, alignmentContext)
  where
    symbolCount = matrixDimension denseTCMs
    f :: Vector CUInt -> Vector CUInt -> (Word, SlimDynamicCharacter)
    f lesser longer = unsafePerformIO . V.unsafeWith lesser $ \lesserPtr -> V.unsafeWith longer $ \longerPtr -> do
        let lesserLength = V.length lesser
        let longerLength = V.length longer
        -- Add two because the C code needs stupid gap prepended to each character.
        -- Forgetting to do this will eventually corrupt the heap memory
        let bufferLength = lesserLength + longerLength + 2
        lesserBuffer <- allocCharacterBuffer bufferLength lesserLength lesserPtr
        longerBuffer <- allocCharacterBuffer bufferLength longerLength longerPtr
        medianBuffer <- allocCharacterBuffer bufferLength            0   nullPtr
        resultLength <- malloc :: IO (Ptr CSize)
        strategy  <- getAlignmentStrategy <$> peek costStruct
        let medianOpt = coerceEnum computeMedians
        let !cost = case strategy of
                      Affine -> {-# SCC affine_undefined #-}
                        undefined -- align2dAffineFn_c lesserBuffer longerBuffer medianBuffer resultLength (ics bufferLength) (ics lesserLength) (ics longerLength) costStruct medianOpt
                      _      -> {-# SCC align2dFn_c #-}
                        align2dFn_c lesserBuffer longerBuffer medianBuffer resultLength (ics bufferLength) (ics lesserLength) (ics longerLength) costStruct neverComputeOnlyGapped medianOpt (coerceEnum computeUnion)

        alignedLength <- coerce <$> peek resultLength
        let g = buildResult bufferLength (csi alignedLength)
        alignedLesser <- {-# SCC alignedLesser #-} g lesserBuffer
        alignedLonger <- {-# SCC alignedLonger #-} g longerBuffer
        alignedMedian <- {-# SCC alignedMedian #-} g medianBuffer

{--
        !_ <- trace ("\n Len: " <> show (length resultingAlignedChar1)) $ pure ()
        !_ <- trace ("  Gapped Char: " <> show       resultingGapped) $ pure ()
        !_ <- trace (" Aligned LHS : " <> show resultingAlignedChar1) $ pure ()
        !_ <- trace (" Aligned RHS : " <> show resultingAlignedChar2) $ pure ()
--}
        pure $ {-# SCC ffi_result #-} (fromIntegral cost, (alignedMedian,alignedLesser,alignedLonger))

      where
        costStruct  = costMatrix2D denseTCMs
        neverComputeOnlyGapped = 0
        ics :: Int -> CSize
        ics = coerce . (toEnum :: Int -> Word64)
        csi :: CSize -> Int
        csi = (fromEnum :: Word64 -> Int) . coerce 

{-
-- |
-- Performs a naive direct optimization
-- Takes in three characters to run DO on and a metadata object.
-- Returns an assignment character, the cost of that assignment, the assignment character with gaps included,
-- the aligned versions of the three input characters.
-- The process for this algorithm is to generate a traversal matrix, then perform a traceback.
algn3d
  :: ( EncodableDynamicCharacter s
     , ExportableElements s
     )
  => s                         -- ^ First  dynamic character
  -> s                         -- ^ Second dynamic character
  -> s                         -- ^ Third  dynamic character
  -> Int                       -- ^ Mismatch cost
  -> Int                       -- ^ Gap open cost
  -> Int                       -- ^ Indel cost
  -> DenseTransitionCostMatrix -- ^ Structure defining the transition costs between character states
  -> (Word, s, s, s, s, s)     -- ^ The cost of the alignment
                               --
                               --   The /ungapped/ character derived from the the input characters' N-W-esque matrix traceback
                               --
                               --   The /gapped/ character derived from the the input characters' N-W-esque matrix traceback
                               --
                               --   The gapped alignment of the /first/ input character when aligned with the second & third character
                               --
                               --   The gapped alignment of the /second/ input character when aligned with the first & third character
                               --
                               --   The gapped alignment of the /third/ input character when aligned with the first & second character
                               --
algn3d char1 char2 char3 mismatchCost openningGapCost indelCost denseTCMs = handleMissingCharacterThreeway someFun char1 char2 char3 $
    case (toExportableElements char1, toExportableElements char2, toExportableElements char3) of
      (Just x, Just y, Just z) -> f x y z
      (     _,      _,      _) -> error "3DO: There's a dynamic character missing!"
  where
    someFun = undefined
    f exportedChar1 exportedChar2 exportedChar3 = unsafePerformIO $ do
--        !_ <- trace ("char 1: " <> show char1) $ pure ()
--        !_ <- trace ("char 2: " <> show char2) $ pure ()
        char1ToSend <- allocInitAlign_io bufferLength exportedChar1Len . fmap coerceEnum $ exportedCharacterElements exportedChar1
        char2ToSend <- allocInitAlign_io bufferLength exportedChar2Len . fmap coerceEnum $ exportedCharacterElements exportedChar2
        char3ToSend <- allocInitAlign_io bufferLength exportedChar3Len . fmap coerceEnum $ exportedCharacterElements exportedChar3
        char1Return <- allocInitAlign_io bufferLength 0 []    -- Note that the next six can be empty as their C-side
        char2Return <- allocInitAlign_io bufferLength 0 []    -- internal arrays are alloc'ed
        char3Return <- allocInitAlign_io bufferLength 0 []
        retGapped   <- allocInitAlign_io bufferLength 0 []
        retUngapped <- allocInitAlign_io bufferLength 0 []
        -- retUnion    <- allocInitALignIO 0 []

        let !cost = align3dFn_c char1ToSend char2ToSend char3ToSend
                                char1Return char2Return char3Return
                                retGapped   retUngapped
                                costStruct
                                (coerceEnum mismatchCost)
                                (coerceEnum openningGapCost)
                                (coerceEnum indelCost)

        resultingAlignedChar1 <- extractFromAlign_io elemWidth char1Return
        resultingAlignedChar2 <- extractFromAlign_io elemWidth char2Return
        resultingAlignedChar3 <- extractFromAlign_io elemWidth char3Return
        resultingGapped       <- extractFromAlign_io elemWidth retGapped
        resultingUngapped     <- extractFromAlign_io elemWidth retUngapped

        pure ( fromIntegral cost
             , resultingUngapped
             , resultingGapped
             , resultingAlignedChar1
             , resultingAlignedChar2
             , resultingAlignedChar3
             )

      where
        costStruct       = costMatrix3D denseTCMs -- TODO: get memoized matrix wedged in here

        elemWidth        = exportedChar1 ^. exportedElementWidth

        exportedChar1Len = coerceEnum $ exportedChar1 ^. exportedElementCount
        exportedChar2Len = coerceEnum $ exportedChar2 ^. exportedElementCount
        exportedChar3Len = coerceEnum $ exportedChar3 ^. exportedElementCount

        bufferLength      = exportedChar1Len + exportedChar2Len + exportedChar3Len
-}


{- Generic helper functions -}


-- |
-- Allocates space for an align_io struct to be sent to C.
allocCharacterBuffer :: Int -> Int  -> Ptr CUInt -> IO (Ptr CUInt)
allocCharacterBuffer maxSize elemCount elements = do
    let e   = min maxSize elemCount
    buffer <- mallocBytes maxSize
    let ref = plusPtr buffer . fromEnum $ maxSize - e
    copyArray ref elements e
    pure $ coerce buffer


buildResult :: Int -> Int -> Ptr CUInt -> IO (Vector CUInt)
buildResult bufferLength alignedLength alignedBuffer = do
    let e   = min bufferLength alignedLength
    vector <- mallocBytes e
    let ref = plusPtr alignedBuffer . fromEnum $ bufferLength - e
    copyArray vector ref e
    free alignedBuffer
    fPtr <- newConcForeignPtr vector (free vector)
    let y = V.unsafeFromForeignPtr0 fPtr e :: Vector CUInt
    pure $ coerce y


{-
-- |
-- Converts the data behind an 'Align_io' pointer to an 'Exportable' type.
{-# INLINE extractFromAlign_io #-}
{-# SPECIALISE extractFromAlign_io :: Word -> Ptr Align_io -> IO DynamicCharacter #-}
extractFromAlign_io :: ExportableElements s => Word -> Ptr Align_io -> IO s
extractFromAlign_io elemWidth ptr = do
    Align_io bufferPtr charLenC bufferLenC <- peek ptr
    let    charLength = fromEnum   charLenC
    let  bufferLength = fromEnum bufferLenC
    buffer <- peekArray bufferLength bufferPtr
    let !charElems = drop (bufferLength - charLength) buffer
    let  exportVal = ReImportableCharacterElements (toEnum charLength) elemWidth charElems
    _ <- free bufferPtr
    _ <- free ptr
    pure $ fromExportableElements exportVal
-}


-- |
-- Coercing one 'Enum' to another through their corresponding 'Int' values.
{-# INLINE coerceEnum #-}
{-# SPECIALISE coerceEnum :: Word  -> Int   #-}
{-# SPECIALISE coerceEnum :: Int   -> Word  #-}
{-# SPECIALISE coerceEnum :: Int   -> CUInt #-}
{-# SPECIALISE coerceEnum :: CUInt -> Int   #-}
{-# SPECIALISE coerceEnum :: Word  -> CUInt #-}
{-# SPECIALISE coerceEnum :: CUInt -> Word  #-}
coerceEnum :: (Enum a, Enum b) => a -> b
coerceEnum = toEnum . fromEnum

