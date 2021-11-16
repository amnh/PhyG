{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies     #-}

module Main (main) where

import Bio.DynamicCharacter
import Data.Foldable
import Data.List                                  (isPrefixOf, uncons)
import Data.MetricRepresentation
import Data.Word
import System.Environment                         (getArgs)
import Test.Aligners
import Test.QuickCheck.Instances.DynamicCharacter
import Test.Tasty
import Test.Tasty.QuickCheck


data OperationalMode a
   = Search
   | SearchAgainst    a
   | RenderComparison a a
   | TooManyParameters


main :: IO ()
main = do
    args <- getArgs
    case parseArgs args of
      Search                       -> performCounterExampleSearch' Nothing
      SearchAgainst    char1       -> performCounterExampleSearch' $  Just  char1
      RenderComparison char1 char2 -> performImplementationComparison char1 char2
      TooManyParameters            -> putStrLn "Expecting only two parameters!"


parseArgs :: [String] -> OperationalMode DNA
parseArgs args =
    case args of
      []          -> Search
      arg1:_ | "--quickcheck-replay" `isPrefixOf` arg1 -> Search
      arg1:xs     ->
        let char1 = readSequence arg1
        in  case xs of
               []     -> SearchAgainst    char1
               [arg2] -> RenderComparison char1 (readSequence arg2)
               _      -> TooManyParameters
  where
    readSequence xs =
       let (a,ys) = break (=='|') xs
           (b,zs) = break (=='|') $ tail ys
           c      = tail zs
       in  if length a == length b && length b == length c
           then buildDNA a b c
           else error $ unlines
                  [ show xs
                  , show a
                  , show b
                  , show c
                  ]


performCounterExampleSearch' :: Maybe DNA -> IO ()
performCounterExampleSearch' valueMay =
        case valueMay of
          Nothing  -> makeMain . testProperty preamble $ uncurry counterExampleCheck
          Just dna -> makeMain . testProperty preamble $ counterExampleCheck dna
  where
    preamble = "Performing stochastic counter-example search"
    makeMain = defaultMain . localOption (QuickCheckTests 1000000) . localOption (QuickCheckShowReplay True)


counterExampleCheck :: DNA -> DNA -> Property
counterExampleCheck lhs rhs = uncurry counterexample $ gatherContexts chosenMetric lhs rhs


performImplementationComparison :: DNA -> DNA -> IO ()
performImplementationComparison lhs rhs = putStrLn . fst $ gatherContexts chosenMetric lhs rhs


chosenMetric :: MetricRepresentation Word32
chosenMetric = snd . head $ metricChoices


gatherContexts
  :: MetricRepresentation Word32
  -> DNA
  -> DNA
  -> (String, Bool)
gatherContexts tcm (DNA lhs) (DNA rhs) = (contextRendering, contextSameness)
  where
    alignmentλs :: [(String, SlimDynamicCharacter -> SlimDynamicCharacter -> (Word, SlimDynamicCharacter))]
    alignmentλs = selectAlignmentλ tcm

    alignmentResults =
        let performAlignment (x,f) = (x, f lhs rhs)
        in  performAlignment <$> alignmentλs

    contextSameness  = sameAlignment $ snd <$> alignmentResults

    contextRendering = renderContexts lhs rhs alignmentResults


sameAlignment :: (Foldable t, Eq c) => t (c, SlimDynamicCharacter) -> Bool
sameAlignment v =
    case uncons $ toList v of
      Nothing -> True
      Just ((c,a),xs) ->
          let sameCost = all (== c) $ fst <$> xs
              sameStr  = all (== a) $ snd <$> xs
          in  sameCost && sameStr


renderContexts
 :: ( Eq c
    , Foldable f
    , Functor f
    , Show c
    )
 => SlimDynamicCharacter
 -> SlimDynamicCharacter
 -> f (String, (c, SlimDynamicCharacter))
 -> String
renderContexts m n xs = unlines . (\x -> [prefix] <> x <> [suffix]) . fmap f $ toList xs
  where
    f (s, c) = s <> "\n" <> renderResult c
    renderResult (cost, aligned) = unlines
        [ "  Cost: "  <> show cost
        , "  Alignment:\n" <> renderAlignment aligned
        ]

    renderAlignment = unlines . fmap ("    " <>) . lines . show . DNA

    noDiff = sameAlignment $ snd <$> xs
    prefix = unlines [ show m, show n ]
    suffix
      | noDiff    = "[!] Results MATCH"
      | otherwise = "[X] Results DO NOT MATCH"


selectAlignmentλ
  :: MetricRepresentation Word32
  -> [ (String, SlimDynamicCharacter -> SlimDynamicCharacter -> (Word, SlimDynamicCharacter)) ]
selectAlignmentλ metric = fmap ($ metric) <$> alignmentChoices
