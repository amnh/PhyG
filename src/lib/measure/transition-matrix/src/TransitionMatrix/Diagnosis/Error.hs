-----------------------------------------------------------------------------
-- |
-- Module      :  TransitionMatrix.Diagnosis.Error
-- Copyright   :  (c) 2015-2021 Ward Wheeler
-- License     :  BSD-style
--
-- Maintainer  :  wheeler@amnh.org
-- Stability   :  provisional
-- Portability :  portable
--
-----------------------------------------------------------------------------

{-# LANGUAGE DeriveAnyClass        #-}
{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE DeriveGeneric         #-}
{-# LANGUAGE DerivingStrategies    #-}
{-# LANGUAGE FlexibleContexts      #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE LambdaCase            #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE Strict                #-}
{-# LANGUAGE UnboxedSums           #-}

module TransitionMatrix.Diagnosis.Error
  ( -- * Failure
    DiagnosisFailure()
    -- * Error types
  , DiagnosisError(..)
    -- * Smart constructor
  , makeDiagnosisFailure
    -- * Error extraction
  , getDiagnosisErrors
  ) where

import           Control.DeepSeq
import           Data.Char (toLower)
import           Data.Data
import           Data.Foldable
import           Data.IntSet                     (IntSet)
import qualified Data.IntSet                     as IS
import           Data.List          (intercalate, partition, sort)
import           Data.List.NonEmpty (NonEmpty(..))
import           Data.Set                        (Set)
import           Data.Word                       (Word16)
import           GHC.Generics
import           GHC.Exts (fromList)


data  DiagnosisError a
    = JaggedColumns !Word !IntSet
    | JaggedRows    !Word !IntSet
    | NotSquareGrid !Word !Word
    | NotSquareList !Word
    | MatrixDimension0
    | MatrixDimension1
    | ValueNegative !(Set a)
    | ValueOverflow !(Set a)
    deriving stock    (Data, Eq, Ord, Generic, Typeable)
    deriving anyclass (NFData)


newtype DiagnosisFailure a = DiagnosisFailure (NonEmpty (DiagnosisError a))
    deriving newtype  (NFData, Eq, Typeable)
    deriving stock    (Data, Generic)


instance (Ord a) => Semigroup (DiagnosisFailure a) where

    (<>) (DiagnosisFailure lhs@(x:|xs)) (DiagnosisFailure rhs) =
        let rhs' = filter (`notElem` lhs) $ toList rhs
            end  = sort $ xs <> rhs'
        in  DiagnosisFailure $ x :| end


instance Show a => Show (DiagnosisError a) where

    show =
        \case
         JaggedColumns n xs -> intercalate "\n"
             [ "Unequal column lengths in the proffered transition measure range."
             , fold [ "Expected all columns to be of length ", show n, " but found the following:"]
             , nicelyRenderList $ IS.toList xs
             ]
         JaggedRows    n xs -> intercalate "\n"
             [ "Unequal row lengths in the proffered transition measure range."
             , fold [ "Expected all rows to be of length ", show n, " but found the following:"]
             , nicelyRenderList $ IS.toList xs
             ]
         NotSquareGrid m n -> intercalate "\n"
             [ "Unequal dimensions proffered for transition measure range."
             , fold [ "Cannot construct a non-square, ", show m, "⨉", show n, " transition measure" ]
             ]
         NotSquareList n   -> intercalate "\n"
             [ "List of non-square length proffered as transition measure range."
             , let root  = sqrt $ fromIntegral n :: Double
                   lower = floor   root :: Word
                   upper = ceiling root :: Word
               in  unwords
                   [ "The number of elements", show n
                   , "lies between the valid square numbers", show lower
                   , "and", show upper <> "."
                   ]
             ]
         MatrixDimension0 -> intercalate "\n"
             [ "0-dimensional transition measure proffered."
             , "Cannot construct an empty transition measure"
             ]
         MatrixDimension1 -> intercalate "\n"
             [ "1-dimensional transition measure proffered."
             , "Cannot construct a transition measure of size 1, must be of size 2 or greater."
             ]
         ValueNegative xs -> intercalate "\n"
             [ "Negative values in the proffered transition measure range."
             , "Expected only non-negative values in range but found the following:"
             , nicelyRenderList $ toList xs
             ]
         ValueOverflow xs -> intercalate "\n"
             [ "Values exceeding " <> show (maxBound :: Word16) <> " in the transition measure range."
             , "Expected smaller values in range but found the following:"
             , nicelyRenderList $ toList xs
             ]


instance Show a => Show (DiagnosisFailure a) where

    show =
        let tenseSensitiveShow = \case
                x:|[] ->
                    case show x of
                        -- Empty case *should* never happen!
                       []    -> "Transition Measure diagnosis failed for an unknown reason"
                       c:str -> "Transition Measure diagnosis failed due to " <> (toLower c : str)
                xs -> fold [ "Transition Measure diagnosis failed for the following reasons:\n"
                           , prettifyErr xs
                           ]
            bullet    []  = []
            bullet (x:xs) = ("  • " <> drop 4 x):xs
            prettifyErr = intercalate "\n" . fmap renderError . toList
            indentError = intercalate "\n" . bullet . fmap ("    " <>) . lines
            renderError = indentError . show
        in  tenseSensitiveShow . getDiagnosisErrors


nicelyRenderList :: Show a => [a] -> String
nicelyRenderList =
    let limit = 75
        showTokens = words . intercalate ", " . fmap show
        revUnwords = unwords . reverse
        revUnlines = intercalate "\n" . reverse
        breakLines =
            let go acc = \case
                    [] -> revUnlines acc
                    ts ->
                        let (taken, leftover) = takeLine (0, []) ts
                        in  go (taken:acc) leftover
            in  go []
        takeLine (n, taken) = \case
            [] -> (revUnwords taken, [])
            toks@(t:ts) ->
                let n' = n + length t
                in  if   n' > limit
                    then (revUnwords taken, toks) 
                    else takeLine (n' + 1, t:taken) ts
    in  breakLines . showTokens


-- |
-- Ensures that ther is at most /one/ 'DiagnosisFailure' which satisfies the
-- predicate 'malformedInputShape' and that if such an error exists, that it is
-- the first error in the list.
makeDiagnosisFailure :: Ord a => NonEmpty (DiagnosisError a) -> DiagnosisFailure a
makeDiagnosisFailure xs =
    let sortedErrors = sort $ toList xs
        uniqueErrors = pruneDuplicates sortedErrors
        (badShapes, other) = partition malformedInputShape uniqueErrors
    in  DiagnosisFailure $ case badShapes of
            []  -> fromList other
            x:_ -> x :| other


-- |
-- Extract the list of diagnosed errors with the proffered transition measure.
getDiagnosisErrors :: DiagnosisFailure a -> NonEmpty (DiagnosisError a)
getDiagnosisErrors (DiagnosisFailure xs) = xs


pruneDuplicates :: [DiagnosisError a] -> [DiagnosisError a]
pruneDuplicates      []  = []
pruneDuplicates v@(x:xs) =
    let f = sameErrorType
    in  case xs of
          []           -> v
          y:ys | f x y -> x : pruneDuplicates ys
          _            -> x : pruneDuplicates xs


malformedInputShape :: DiagnosisError a -> Bool
malformedInputShape =
    \case
        JaggedColumns {} -> True
        JaggedRows    {} -> True
        NotSquareGrid {} -> True
        NotSquareList {} -> True
        _                -> False


sameErrorType :: DiagnosisError a -> DiagnosisError a -> Bool
sameErrorType JaggedColumns    {} JaggedColumns    {} = True
sameErrorType JaggedRows       {} JaggedRows       {} = True
sameErrorType NotSquareGrid    {} NotSquareGrid    {} = True
sameErrorType NotSquareList    {} NotSquareList    {} = True
sameErrorType MatrixDimension0 {} MatrixDimension0 {} = True
sameErrorType MatrixDimension1 {} MatrixDimension1 {} = True
sameErrorType ValueNegative    {} ValueNegative    {} = True
sameErrorType ValueOverflow    {} ValueOverflow    {} = True
sameErrorType _ _ = False
