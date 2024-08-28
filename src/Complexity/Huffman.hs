{-# LANGUAGE Safe #-}

{-
Module from:
http://hackage.haskell.org/package/huffman-1.0.1/docs/Data-Compression-Huffman.html
modified to work with newer versino of Containers
-}

module Complexity.Huffman (
    HuffmanTree (..),
    Bit (..),
    Code,
    huffman,
    huffmanSorted,
    codewords,
    ppCode,
) where

import Control.Arrow (first, second)
import Data.List (intercalate)
import Data.PriorityQueue.FingerTree qualified as PQ
import Data.Sequence (
    ViewL (EmptyL, (:<)),
    viewl,
    (|>),
 )
import Data.Sequence qualified as S


data Bit = Zero | One
    deriving stock (Eq)


instance Show Bit where
    show Zero = "0"
    show One = "1"


data HuffmanTree a
    = Empty
    | Node (HuffmanTree a) (HuffmanTree a)
    | Leaf a
    deriving stock (Show)


type role HuffmanTree representational


type Code a = [(a, [Bit])]


{- |
__Time:__ \( \mathcal{O}\left( n \log n \right) \)

Simple implementation.
-}
huffman ∷ ∀ a w. (Ord w, Num w) ⇒ [(a, w)] → HuffmanTree a
huffman =
    let build ∷ PQ.PQueue w (HuffmanTree a) → HuffmanTree a
        build pq = case PQ.minViewWithKey pq of
            Nothing → Empty
            Just ((w, x), pq') →
                case PQ.minViewWithKey pq' of
                    Nothing → x
                    Just ((w', y), pq'') → build $ PQ.insert (w + w') (Node x y) pq''

        prepare ∷ [(a, w)] → PQ.PQueue w (HuffmanTree a)
        prepare = PQ.fromList . map (\(x, w) → (w, Leaf x))
    in  build . prepare


{- |
__Time:__ \( \mathcal{O}\left( n \right) \)

More efficient implementation.
Requires that the input list of symbols and weight is sorted by increasing weight.
-}
huffmanSorted ∷ ∀ a w. (Ord w, Num w) ⇒ [(a, w)] → HuffmanTree a
huffmanSorted =
    let prepare ∷ [(a, w)] → S.Seq (HuffmanTree a, w)
        prepare = S.fromList . map (first Leaf)

        dequeue
            ∷ S.Seq (HuffmanTree a, w)
            → S.Seq (HuffmanTree a, w)
            → Maybe ((HuffmanTree a, w), S.Seq (HuffmanTree a, w), S.Seq (HuffmanTree a, w))
        dequeue s t =
            case (viewl s, viewl t) of
                (EmptyL, EmptyL) → Nothing
                (EmptyL, x :< ts) → Just (x, s, ts)
                (x :< ss, EmptyL) → Just (x, ss, t)
                ((x, w) :< ss, (y, w') :< ts)
                    | w < w' → Just ((x, w), ss, t)
                    | otherwise → Just ((y, w'), s, ts)

        build ∷ S.Seq (HuffmanTree a, w) → S.Seq (HuffmanTree a, w) → HuffmanTree a
        build s t =
            case dequeue s t of
                Nothing → Empty
                Just ((x, w), s', t') →
                    case dequeue s' t' of
                        Nothing → x
                        Just ((y, w'), s'', t'') → build (s'' |> (Node x y, w + w')) t''
    in  build S.empty . prepare


-- Derive the prefix-free binary code from a huffman tree.
codewords ∷ HuffmanTree a → Code a
codewords = code' []
    where
        code' ∷ [Bit] → HuffmanTree a → [(a, [Bit])]
        code' _ Empty = []
        code' bits (Leaf x) = [(x, bits)]
        code' bits (Node l r) =
            map (second (Zero :)) (code' bits l)
                ++ map (second (One :)) (code' bits r)


-- Pretty-print a binary code.  Mostly useful for debugging.
ppCode ∷ (Show a) ⇒ Code a → String
ppCode =
    intercalate "\n"
        . map (\(x, bits) → show x ++ ": " ++ concatMap show bits)
