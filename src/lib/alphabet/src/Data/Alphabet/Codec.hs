-----------------------------------------------------------------------------
-- |
-- Module      :  Data.Alphabet.Codec
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-- Facilitates encoding and decoding symbols from an 'Alphabet' into a bit-state.
-- Works for any 'Bits' instance.
-----------------------------------------------------------------------------

{-# LANGUAGE GADTs  #-}
{-# LANGUAGE Strict #-}

module Data.Alphabet.Codec
  ( decodeState
  , encodeState
  ) where

import           Data.Alphabet.Internal
import           Data.Bits
import           Data.Foldable
import           Data.List.NonEmpty     (NonEmpty)
import           Data.Set               (Set)
import qualified Data.Set               as Set
import           Data.Text.Short        (ShortText)
import qualified Data.Vector            as V
import           Data.Word
import           Foreign.C.Types
import           GHC.Exts               (IsList(fromList), Item)


{-# INLINEABLE encodeState #-}
{-# SPECIALISE encodeState :: (Foldable f, Ord s) => Alphabet s -> (Word -> CUChar ) -> f s -> CUChar  #-}
{-# SPECIALISE encodeState :: (Foldable f, Ord s) => Alphabet s -> (Word -> CUShort) -> f s -> CUShort #-}
{-# SPECIALISE encodeState :: (Foldable f, Ord s) => Alphabet s -> (Word -> CUInt  ) -> f s -> CUInt   #-}
{-# SPECIALISE encodeState :: (Foldable f, Ord s) => Alphabet s -> (Word -> CULong ) -> f s -> CULong  #-}
{-# SPECIALISE encodeState :: (Foldable f, Ord s) => Alphabet s -> (Word -> Word   ) -> f s -> Word    #-}
{-# SPECIALISE encodeState :: (Foldable f, Ord s) => Alphabet s -> (Word -> Word8  ) -> f s -> Word8   #-}
{-# SPECIALISE encodeState :: (Foldable f, Ord s) => Alphabet s -> (Word -> Word16 ) -> f s -> Word16  #-}
{-# SPECIALISE encodeState :: (Foldable f, Ord s) => Alphabet s -> (Word -> Word32 ) -> f s -> Word32  #-}
{-# SPECIALISE encodeState :: (Foldable f, Ord s) => Alphabet s -> (Word -> Word64 ) -> f s -> Word64  #-}
{-# SPECIALISE encodeState :: Ord s => Alphabet s -> (Word -> CUChar ) -> Set s -> CUChar  #-}
{-# SPECIALISE encodeState :: Ord s => Alphabet s -> (Word -> CUShort) -> Set s -> CUShort #-}
{-# SPECIALISE encodeState :: Ord s => Alphabet s -> (Word -> CUInt  ) -> Set s -> CUInt   #-}
{-# SPECIALISE encodeState :: Ord s => Alphabet s -> (Word -> CULong ) -> Set s -> CULong  #-}
{-# SPECIALISE encodeState :: Ord s => Alphabet s -> (Word -> Word   ) -> Set s -> Word    #-}
{-# SPECIALISE encodeState :: Ord s => Alphabet s -> (Word -> Word8  ) -> Set s -> Word8   #-}
{-# SPECIALISE encodeState :: Ord s => Alphabet s -> (Word -> Word16 ) -> Set s -> Word16  #-}
{-# SPECIALISE encodeState :: Ord s => Alphabet s -> (Word -> Word32 ) -> Set s -> Word32  #-}
{-# SPECIALISE encodeState :: Ord s => Alphabet s -> (Word -> Word64 ) -> Set s -> Word64  #-}
{-# SPECIALISE encodeState :: Alphabet String -> (Word -> CUChar ) -> Set String -> CUChar  #-}
{-# SPECIALISE encodeState :: Alphabet String -> (Word -> CUShort) -> Set String -> CUShort #-}
{-# SPECIALISE encodeState :: Alphabet String -> (Word -> CUInt  ) -> Set String -> CUInt   #-}
{-# SPECIALISE encodeState :: Alphabet String -> (Word -> CULong ) -> Set String -> CULong  #-}
{-# SPECIALISE encodeState :: Alphabet String -> (Word -> Word   ) -> Set String -> Word    #-}
{-# SPECIALISE encodeState :: Alphabet String -> (Word -> Word8  ) -> Set String -> Word8   #-}
{-# SPECIALISE encodeState :: Alphabet String -> (Word -> Word16 ) -> Set String -> Word16  #-}
{-# SPECIALISE encodeState :: Alphabet String -> (Word -> Word32 ) -> Set String -> Word32  #-}
{-# SPECIALISE encodeState :: Alphabet String -> (Word -> Word64 ) -> Set String -> Word64  #-}
{-# SPECIALISE encodeState :: Alphabet ShortText -> (Word -> CUChar ) -> Set ShortText -> CUChar  #-}
{-# SPECIALISE encodeState :: Alphabet ShortText -> (Word -> CUShort) -> Set ShortText -> CUShort #-}
{-# SPECIALISE encodeState :: Alphabet ShortText -> (Word -> CUInt  ) -> Set ShortText -> CUInt   #-}
{-# SPECIALISE encodeState :: Alphabet ShortText -> (Word -> CULong ) -> Set ShortText -> CULong  #-}
{-# SPECIALISE encodeState :: Alphabet ShortText -> (Word -> Word   ) -> Set ShortText -> Word    #-}
{-# SPECIALISE encodeState :: Alphabet ShortText -> (Word -> Word8  ) -> Set ShortText -> Word8   #-}
{-# SPECIALISE encodeState :: Alphabet ShortText -> (Word -> Word16 ) -> Set ShortText -> Word16  #-}
{-# SPECIALISE encodeState :: Alphabet ShortText -> (Word -> Word32 ) -> Set ShortText -> Word32  #-}
{-# SPECIALISE encodeState :: Alphabet ShortText -> (Word -> Word64 ) -> Set ShortText -> Word64  #-}
{-# SPECIALISE encodeState :: Alphabet String -> (Word -> CUChar ) -> NonEmpty String -> CUChar  #-}
{-# SPECIALISE encodeState :: Alphabet String -> (Word -> CUShort) -> NonEmpty String -> CUShort #-}
{-# SPECIALISE encodeState :: Alphabet String -> (Word -> CUInt  ) -> NonEmpty String -> CUInt   #-}
{-# SPECIALISE encodeState :: Alphabet String -> (Word -> CULong ) -> NonEmpty String -> CULong  #-}
{-# SPECIALISE encodeState :: Alphabet String -> (Word -> Word   ) -> NonEmpty String -> Word    #-}
{-# SPECIALISE encodeState :: Alphabet String -> (Word -> Word8  ) -> NonEmpty String -> Word8   #-}
{-# SPECIALISE encodeState :: Alphabet String -> (Word -> Word16 ) -> NonEmpty String -> Word16  #-}
{-# SPECIALISE encodeState :: Alphabet String -> (Word -> Word32 ) -> NonEmpty String -> Word32  #-}
{-# SPECIALISE encodeState :: Alphabet String -> (Word -> Word64 ) -> NonEmpty String -> Word64  #-}
{-# SPECIALISE encodeState :: Alphabet ShortText -> (Word -> CUChar ) -> NonEmpty ShortText -> CUChar  #-}
{-# SPECIALISE encodeState :: Alphabet ShortText -> (Word -> CUShort) -> NonEmpty ShortText -> CUShort #-}
{-# SPECIALISE encodeState :: Alphabet ShortText -> (Word -> CUInt  ) -> NonEmpty ShortText -> CUInt   #-}
{-# SPECIALISE encodeState :: Alphabet ShortText -> (Word -> CULong ) -> NonEmpty ShortText -> CULong  #-}
{-# SPECIALISE encodeState :: Alphabet ShortText -> (Word -> Word   ) -> NonEmpty ShortText -> Word    #-}
{-# SPECIALISE encodeState :: Alphabet ShortText -> (Word -> Word8  ) -> NonEmpty ShortText -> Word8   #-}
{-# SPECIALISE encodeState :: Alphabet ShortText -> (Word -> Word16 ) -> NonEmpty ShortText -> Word16  #-}
{-# SPECIALISE encodeState :: Alphabet ShortText -> (Word -> Word32 ) -> NonEmpty ShortText -> Word32  #-}
{-# SPECIALISE encodeState :: Alphabet ShortText -> (Word -> Word64 ) -> NonEmpty ShortText -> Word64  #-}
encodeState
  :: ( Bits e
     , Foldable f
     , Ord s
     )
  => Alphabet s  -- ^ Alphabet of symbols
  -> (Word -> e) -- ^ Constructor for an empty element, taking the alphabet size
  -> f s         -- ^ ambiguity groups of symbols
  -> e           -- ^ Encoded dynamic character element
encodeState alphabet f symbols = getSubsetIndex alphabet symbolsSet emptyElement
  where
    emptyElement = f . toEnum $ length alphabet
    symbolsSet   = Set.fromList $ toList symbols


{-# INLINEABLE decodeState #-}
{-# SPECIALISE decodeState :: Alphabet s -> CUChar  -> [s] #-}
{-# SPECIALISE decodeState :: Alphabet s -> CUShort -> [s] #-}
{-# SPECIALISE decodeState :: Alphabet s -> CUInt   -> [s] #-}
{-# SPECIALISE decodeState :: Alphabet s -> CULong  -> [s] #-}
{-# SPECIALISE decodeState :: Alphabet s -> Word    -> [s] #-}
{-# SPECIALISE decodeState :: Alphabet s -> Word8   -> [s] #-}
{-# SPECIALISE decodeState :: Alphabet s -> Word16  -> [s] #-}
{-# SPECIALISE decodeState :: Alphabet s -> Word32  -> [s] #-}
{-# SPECIALISE decodeState :: Alphabet s -> Word64  -> [s] #-}
{-# SPECIALISE decodeState :: Alphabet String -> CUChar  -> [String] #-}
{-# SPECIALISE decodeState :: Alphabet String -> CUShort -> [String] #-}
{-# SPECIALISE decodeState :: Alphabet String -> CUInt   -> [String] #-}
{-# SPECIALISE decodeState :: Alphabet String -> CULong  -> [String] #-}
{-# SPECIALISE decodeState :: Alphabet String -> Word    -> [String] #-}
{-# SPECIALISE decodeState :: Alphabet String -> Word8   -> [String] #-}
{-# SPECIALISE decodeState :: Alphabet String -> Word16  -> [String] #-}
{-# SPECIALISE decodeState :: Alphabet String -> Word32  -> [String] #-}
{-# SPECIALISE decodeState :: Alphabet String -> Word64  -> [String] #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> CUChar  -> [ShortText] #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> CUShort -> [ShortText] #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> CUInt   -> [ShortText] #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> CULong  -> [ShortText] #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> Word    -> [ShortText] #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> Word8   -> [ShortText] #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> Word16  -> [ShortText] #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> Word32  -> [ShortText] #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> Word64  -> [ShortText] #-}
{-# SPECIALISE decodeState :: Ord s => Alphabet s -> CUChar  -> Set s #-}
{-# SPECIALISE decodeState :: Ord s => Alphabet s -> CUShort -> Set s #-}
{-# SPECIALISE decodeState :: Ord s => Alphabet s -> CUInt   -> Set s #-}
{-# SPECIALISE decodeState :: Ord s => Alphabet s -> CULong  -> Set s #-}
{-# SPECIALISE decodeState :: Ord s => Alphabet s -> Word    -> Set s #-}
{-# SPECIALISE decodeState :: Ord s => Alphabet s -> Word8   -> Set s #-}
{-# SPECIALISE decodeState :: Ord s => Alphabet s -> Word16  -> Set s #-}
{-# SPECIALISE decodeState :: Ord s => Alphabet s -> Word32  -> Set s #-}
{-# SPECIALISE decodeState :: Ord s => Alphabet s -> Word64  -> Set s #-}
{-# SPECIALISE decodeState :: Alphabet String -> CUChar  -> Set String #-}
{-# SPECIALISE decodeState :: Alphabet String -> CUShort -> Set String #-}
{-# SPECIALISE decodeState :: Alphabet String -> CUInt   -> Set String #-}
{-# SPECIALISE decodeState :: Alphabet String -> CULong  -> Set String #-}
{-# SPECIALISE decodeState :: Alphabet String -> Word    -> Set String #-}
{-# SPECIALISE decodeState :: Alphabet String -> Word8   -> Set String #-}
{-# SPECIALISE decodeState :: Alphabet String -> Word16  -> Set String #-}
{-# SPECIALISE decodeState :: Alphabet String -> Word32  -> Set String #-}
{-# SPECIALISE decodeState :: Alphabet String -> Word64  -> Set String #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> CUChar  -> Set ShortText #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> CUShort -> Set ShortText #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> CUInt   -> Set ShortText #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> CULong  -> Set ShortText #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> Word    -> Set ShortText #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> Word8   -> Set ShortText #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> Word16  -> Set ShortText #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> Word32  -> Set ShortText #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> Word64  -> Set ShortText #-}
{-# SPECIALISE decodeState :: Alphabet s -> CUChar  -> NonEmpty s #-}
{-# SPECIALISE decodeState :: Alphabet s -> CUShort -> NonEmpty s #-}
{-# SPECIALISE decodeState :: Alphabet s -> CUInt   -> NonEmpty s #-}
{-# SPECIALISE decodeState :: Alphabet s -> CULong  -> NonEmpty s #-}
{-# SPECIALISE decodeState :: Alphabet s -> Word    -> NonEmpty s #-}
{-# SPECIALISE decodeState :: Alphabet s -> Word8   -> NonEmpty s #-}
{-# SPECIALISE decodeState :: Alphabet s -> Word16  -> NonEmpty s #-}
{-# SPECIALISE decodeState :: Alphabet s -> Word32  -> NonEmpty s #-}
{-# SPECIALISE decodeState :: Alphabet s -> Word64  -> NonEmpty s #-}
{-# SPECIALISE decodeState :: Alphabet String -> CUChar  -> NonEmpty String #-}
{-# SPECIALISE decodeState :: Alphabet String -> CUShort -> NonEmpty String #-}
{-# SPECIALISE decodeState :: Alphabet String -> CUInt   -> NonEmpty String #-}
{-# SPECIALISE decodeState :: Alphabet String -> CULong  -> NonEmpty String #-}
{-# SPECIALISE decodeState :: Alphabet String -> Word    -> NonEmpty String #-}
{-# SPECIALISE decodeState :: Alphabet String -> Word8   -> NonEmpty String #-}
{-# SPECIALISE decodeState :: Alphabet String -> Word16  -> NonEmpty String #-}
{-# SPECIALISE decodeState :: Alphabet String -> Word32  -> NonEmpty String #-}
{-# SPECIALISE decodeState :: Alphabet String -> Word64  -> NonEmpty String #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> CUChar  -> NonEmpty ShortText #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> CUShort -> NonEmpty ShortText #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> CUInt   -> NonEmpty ShortText #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> CULong  -> NonEmpty ShortText #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> Word    -> NonEmpty ShortText #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> Word8   -> NonEmpty ShortText #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> Word16  -> NonEmpty ShortText #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> Word32  -> NonEmpty ShortText #-}
{-# SPECIALISE decodeState :: Alphabet ShortText -> Word64  -> NonEmpty ShortText #-}
decodeState
  :: ( Bits e
     , IsList (f s)
     , Item (f s) ~ s
     )
  => Alphabet s  -- ^ Alphabet of symbols
  -> e           -- ^ State to decode
  -> f s
decodeState alphabet state = fromList $ foldr pollSymbol mempty indices
  where
    indices = [ 0 .. len - 1 ]
    len = length vec
    vec = alphabetSymbols alphabet
    pollSymbol i polled
      | state `testBit` i = (vec V.! i) : polled
      | otherwise         = polled

