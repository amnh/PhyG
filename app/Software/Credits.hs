{-# LANGUAGE ImportQualifiedPost #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE Strict #-}

{- |
Compile-time embeddings of project contributions via Template Haskell.
-}
module Software.Credits (
    contributors,
) where

import Control.Foldl qualified as L
import Control.Monad
import Data.Foldable
import Data.List.NonEmpty (NonEmpty (..))
import Data.Map ((!))
import Data.String (IsString)
import Data.Text hiding (filter, intersperse, replicate)
import Data.Text.Builder.Linear
import Instances.TH.Lift ()
import Software.Metadata.Embedded (embeddedDataFiles)
import Text.MMark
import Text.MMark.Extension
import Text.URI qualified as URI
import Prelude hiding (readFile)


{- |
Rendering of the financial and techncal contributors to the software.
-}
contributors ∷ Builder
contributors =
    let joinLists ∷ (Semigroup a, IsString a) ⇒ a → a → a
        joinLists x y = x <> "\n\n\n" <> y
    in  joinLists authorsList fundingList


{- |
List of authors who have contributed to the PHANE project.
-}
authorsList ∷ Builder
authorsList =
    let drawAuthorLines ∷ [Text] → Builder
        drawAuthorLines rawAuthorLines =
            let renderedAuthorLines = bulletPrefix 1 '•' <$> rawAuthorLines
                renderedHeaderLines =
                    unlines'
                        [ "  │ Project Contributors: │"
                        , "  ╘═══════════════════════╛"
                        ]
            in  seperateLines $ renderedHeaderLines : renderedAuthorLines

        readAuthorLines ∷ Text → [Text]
        readAuthorLines = fmap fst . fromMarkdown processMarkdown

        renderAuthorData ∷ Text → Builder
        renderAuthorData = drawAuthorLines . readAuthorLines
    in  renderAuthorData `fromEmbeddedFileData` "Authors.md"


{- |
List of funding sources which have contributed to PHANE project.
-}
fundingList ∷ Builder
fundingList =
    let drawFunderLines ∷ [(Text, Maybe Text)] → Builder
        drawFunderLines rawFundingSources =
            let listingPrefix = bulletPrefix 1 '•'
                linkingPrefix = bulletPrefix 2 '›'
                processFunder (x, y) = unlines' [listingPrefix x, maybe "" linkingPrefix y]
                renderedFunderLines = processFunder <$> rawFundingSources
                renderedHeaderLines =
                    unlines'
                        [ "  │ Funding Provided By:  │"
                        , "  ╘═══════════════════════╛"
                        ]
            in  seperateLines $ renderedHeaderLines : renderedFunderLines

        readFunderLines ∷ Text → [(Text, Maybe Text)]
        readFunderLines = fromMarkdown processMarkdown

        renderFunderData ∷ Text → Builder
        renderFunderData = drawFunderLines . readFunderLines
    in  renderFunderData `fromEmbeddedFileData` "Funding.md"


fromEmbeddedFileData ∷ (Text → Builder) → FilePath → Builder
fromEmbeddedFileData processor = processor . (embeddedDataFiles !)


processMarkdown ∷ MMark → [(Text, Maybe Text)]
processMarkdown = (`runScanner` L.foldMap g f)
    where
        f ∷ [Block (NonEmpty Inline)] → [(Text, Maybe Text)]
        f = foldMap renderItem

        g ∷ Block a → [Block a]
        g = getListBlocks


getListBlocks ∷ Block a → [Block a]
getListBlocks = fold . foldMap toList . toList . getListMay


getListMay ∷ Block a → Maybe (NonEmpty [Block a])
getListMay (OrderedList _ xs) = Just xs
getListMay (UnorderedList xs) = Just xs
getListMay _ = Nothing


fromMarkdown ∷ (Monoid a) ⇒ (MMark → a) → Text → a
fromMarkdown f = foldMap f . parse ""


renderItem ∷ Block (NonEmpty Inline) → [(Text, Maybe Text)]
renderItem (CodeBlock _ val) = [(val, Nothing)]
renderItem (Naked val) = foldMap renderInline val
renderItem (Paragraph val) = foldMap renderInline val
renderItem _ = []


renderInline ∷ Inline → [(Text, Maybe Text)]
renderInline (Plain txt) = [(txt, Nothing)]
renderInline (Emphasis val) = [(asPlainText val, Nothing)]
renderInline (Strong val) = [(asPlainText val, Nothing)]
renderInline (Strikeout val) = [(asPlainText val, Nothing)]
renderInline (Subscript val) = [(asPlainText val, Nothing)]
renderInline (Superscript val) = [(asPlainText val, Nothing)]
renderInline (CodeSpan txt) = [(txt, Nothing)]
renderInline (Link val uri _) = [(asPlainText val, Just $ URI.render uri)]
renderInline _ = []


bulletPrefix ∷ Word → Char → Text → Builder
bulletPrefix depth bullet =
    let padded = fold $ replicate (fromEnum depth) indentation
        prefix = padded <> fromChar bullet <> " "
    in  (prefix <>) . fromText


indentation ∷ Builder
indentation = "  "


unlines' ∷ (Foldable f) ⇒ f Builder → Builder
unlines' = intercalate' "\n"


intercalate' ∷ (Foldable f) ⇒ Builder → f Builder → Builder
intercalate' sep =
    let consider [] = mempty
        consider [x] = x
        consider (x : xs) = x <> foldMap (sep <>) xs
    in  consider . toList


seperateLines ∷ (Foldable f) ⇒ f Builder → Builder
seperateLines = intercalate' "\n\n"
