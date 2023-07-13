{- |
Command line options parser for the Phylogenetic Graph (PhyG) tool.
-}

{-# Language OverloadedStrings #-}

module CommandLineOptions
  ( -- * Data-types
    CommandLineOptions(getArgsCLI)
    -- * Parser
  , parseCommandLineOptions
    -- * Display information
  , printInformationDisplay
  ) where

import CommandLineOptions.Display
import CommandLineOptions.Types
import Data.Foldable
import Data.Functor ((<&>))
import Data.String
import Options.Applicative
import Options.Applicative.NonEmpty (some1)
import Software.Metadata
import System.Environment
import Prettyprinter ((<+>), align, fillSep, hardline, indent, parens, vsep)


{- |
Command to parse the command line options.
-}
parseCommandLineOptions :: IO CommandLineOptions
parseCommandLineOptions = parserInformation >>= customExecParser parserPreferences


{- |
The preferences for the CLI parser.
-}
parserPreferences :: ParserPrefs
parserPreferences = prefs $ fold
    [ columns 100
    , disambiguate
    , showHelpOnEmpty
    , showHelpOnError
    ]


{- |
Information regarding which command line options are valid and how they are
parsed and interpreted.
-}
parserInformation :: IO (ParserInfo CommandLineOptions)
parserInformation = parserDescription >>= \description -> 
    let displayFlag :: DisplayInfoBlock -> String -> String -> Parser DisplayInfoBlock
        displayFlag val name desc = flag' val $ fold [ long name, help desc ]

        pCommandLineOptions2 =
            let caseInput = Right <$> pInputFile
                casePrint = Left  <$> pDisplayInfos
            in  fmap OptionsCLI $ caseInput <|> casePrint

        pDisplayInfos = some1 pDisplayFlag

        pDisplayFlag = asum
                [ displayFlag DisplayCredits "credits" "Roll project contributions credits"
                , displayFlag DisplayLicense "license" "Emit the software license"
                , displayFlag DisplaySplash  "splash"  "Show splash image"
                , displayFlag DisplayVersion "version" "Show version information"
                ]

        pInputFile = strArgument $ fold
            [ help $ fold ["Filepath to PhyG command script (required)"]
            , metavar "FILE"
            , action "file"
            ]

    in  pure $ info (helper <*> pCommandLineOptions2) description


parserDescription :: IO (InfoMod a)
parserDescription = getProgName <&> \name ->
    let projectLess = fromString $ nameAbbreviation projectName
        projectMore = fromString $ nameExpansion   projectName
        programLess = fromString $ nameAbbreviation softwareName
        programLong = fromString $ nameExpansion   softwareName
        programTell = fromString name

        headInfo :: InfoMod a
        headInfo =
            let projectPreamble = fillSep ["The", projectMore, parens projectLess, "project presents:" ]
                programNameBlock = indent 2 . align $ vsep
                    [ programLong <+> parens programLess
                    , programTell <+> shortVersionInformation
                    ]
            in  headerDoc . Just $ vsep [ projectPreamble <> hardline, programNameBlock ]

        progInfo :: InfoMod a
        progInfo =
            let sentences = 
                    [ programLong : parens programLess : words' "is a multi-platform programram designed to produce phylogenetic graphs from input data and graphs via heuristic searching of general phylogenetic graph space."
                    , words' "The bio-informatics framework libraries of the broader" <> [ projectMore, parens projectLess ] <> words' "project are the foundation upon which" <> [programLess ] <> words' "is constructed."
                    , programLess : words' "offers vast functionality, including the optimization of unaligned sequences, and the ability to implement search strategies such as random addition sequence, swapping, and tree fusing."
                    , "Furthermore," : programLess : words' "can generate outputs in the form of implied alignments and graphical representations of cladograms and graphs."
                    , "What" : "sets" : programLess : words' "apart from other phylogenetic analysis programrams, is the extension to broader classes of input data and phylogenetic graphs."
                    , words' "The phylogenetic graph inputs and outputs of" <> [ programLess ] <> words' "include trees, as well as softwired and hardwired networks."
                    ]
                textStream = align . fillSep $ fold sentences
                words' = fmap fromString . words

            in  progDescDoc $ Just textStream 
        
    in  fold
            [ failureCode 2
            , fullDesc
            , headInfo
            , noIntersperse
            , progInfo
            ]

