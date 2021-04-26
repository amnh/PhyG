{- |
Module      :  ReadInputFiles.hs
Description :  Module to read input files for phylogenetic analysis
Copyright   :  (c) 2021 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
License     :

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

Maintainer  :  Ward Wheeler <wheeler@amnh.org>
Stability   :  unstable
Portability :  portable (I hope)

-}

{-
Functions to peform the input file reading for PhyG
-}

module ReadInputFiles  
(  executeReadCommands
 , getReadArgs
 , extractDataGraphPair
) where

import           Types
import           Debug.Trace
import           Data.Char
import           System.IO
import           Data.List
import qualified Data.Text.Lazy  as T
import qualified Data.Text.Short as ST
import qualified LocalGraph as LG
import qualified GraphFormatUtilities as GFU
import qualified TNTUtilities as TNT
import qualified DataTransformation as DT
import qualified Data.Graph.Inductive.Basic        as B
import qualified GeneralUtilities as GU
import qualified SymMatrix as SM


{-
Added strictness on read so could close files after reading to 
   allow for large number (1000's) of files to be input.
   this may not be the way to go, especially if files are big and 
   closing them is not an issue
-}

-- | extractDataGraphPair  takes the list of pairs from mapped executeReadCommands
-- and returns ([RawData], [SimpleGraph])
extractDataGraphPair :: [([RawData], [SimpleGraph])] -> ([RawData], [SimpleGraph])
extractDataGraphPair dataGraphList =  
    let (inDataList, inGraphList) = unzip dataGraphList
        rawData = concat inDataList
        rawGraphs = concat inGraphList
    in
    (rawData, rawGraphs)



-- | executeReadCommands reads iput files and returns raw data 
-- assumes that "prealigned" and "tcm:File" are tghe first arguments if they are specified 
-- so that they can apply to all teh files in the command without order depence
executeReadCommands :: [RawData] -> [SimpleGraph] -> Bool -> ([ST.ShortText], [[Int]]) -> [Argument] -> IO ([RawData], [SimpleGraph])
executeReadCommands curData curGraphs isPrealigned tcmPair argList = do
    if null argList then return (curData, curGraphs)
    else do
        let isPrealigned' = if isPrealigned then True
                            else if ("prealigned" `elem`  (fmap fst argList)) then True
                            else False
        let (firstOption, firstFile) = head argList
        -- Check for prealigned
        if firstOption == "prealigned" then
            executeReadCommands curData curGraphs isPrealigned' tcmPair (tail argList)
        else do
            fileHandle <- openFile  firstFile ReadMode
            canBeReadFrom <- hIsReadable fileHandle
            if not canBeReadFrom then errorWithoutStackTrace ("\n\n'Read' error: file " ++ firstFile ++ " cannot be read")
            else hPutStrLn stderr ("Reading " ++ firstFile ++ " with option " ++ firstOption)

            -- this is awkward but need to use dot utilities
            if firstOption == "dot" then do
                dotGraph <- LG.hGetDotLocal fileHandle
                let inputDot = GFU.relabelFGL $ LG.dotToGraph dotGraph
                let hasLoops = B.hasLoop inputDot 
                if hasLoops then errorWithoutStackTrace ("Input graph in " ++ firstFile ++ "  has loops/self-edges")
                else hPutStr stderr ""
                let hasCycles = GFU.cyclic inputDot
                if hasCycles then errorWithoutStackTrace ("Input graph in " ++ firstFile ++ " has at least one cycle")
                else executeReadCommands curData (inputDot : curGraphs) isPrealigned' tcmPair (tail argList)
            -- not "dot" files
            else do
                -- destroys lazyness but allows closing right away
                -- this so don't have huge numbers of open files for large data sets
                fileContents <- hGetContents' fileHandle
                hClose fileHandle
                -- try to figure out file type based on first and following characters
                if firstOption == "tcm" then
                    let newTCMPair = processTCMContents fileContents firstFile
                    in
                    executeReadCommands curData curGraphs isPrealigned' newTCMPair (tail argList)
                else if null firstOption then 
                    let firstChar = head $ dropWhile (== ' ') fileContents
                    in
                    if (toLower firstChar == '/') || (toLower firstChar == 'd') || (toLower firstChar == 'g')then do
                        hPutStrLn stderr ("\tTrying to parse " ++ firstFile ++ " as dot") 
                        -- destroys lazyness but allows closing right away
                        -- this so don't have huge numbers of open files for large data sets
                        fileHandle2 <- openFile  firstFile ReadMode
                        !dotGraph <- LG.hGetDotLocal fileHandle2
                        hClose fileHandle2
                        let inputDot = GFU.relabelFGL $ LG.dotToGraph dotGraph
                        let hasLoops = B.hasLoop inputDot 
                        if hasLoops then errorWithoutStackTrace ("Input graph in " ++ firstFile ++ "  has loops/self-edges")
                        else hPutStr stderr ""
                        let hasCycles = GFU.cyclic inputDot
                        if hasCycles then errorWithoutStackTrace ("Input graph in " ++ firstFile ++ " has at least one cycle")
                        else executeReadCommands curData (inputDot : curGraphs) isPrealigned' tcmPair (tail argList)
                    else if (toLower firstChar == '<') || (toLower firstChar == '(')  then
                        let thisGraphList = getFENewickGraph fileContents
                            hasCycles = filter (== True) $ fmap GFU.cyclic thisGraphList
                            hasLoops = filter (== True) $ fmap B.hasLoop thisGraphList 
                        in 
                        if (not $ null hasLoops) then errorWithoutStackTrace ("Input graph in " ++ firstFile ++ "  has loops/self-edges")
                        else if (not $ null hasCycles) then errorWithoutStackTrace ("Input graph in " ++ firstFile ++ " has at least one cycle")
                        else executeReadCommands curData (thisGraphList ++ curGraphs) isPrealigned' tcmPair (tail argList)

                    else if toLower firstChar == 'x' then 
                        let tntData = TNT.getTNTData fileContents firstFile
                        in
                        trace ("\tTrying to parse " ++ firstFile ++ " as TNT") 
                        executeReadCommands (tntData : curData) curGraphs isPrealigned' tcmPair (tail argList)
                    else 
                        let fileContents' =  unlines $ filter (not.null) $ fmap (takeWhile (/= ';')) $ lines fileContents
                        in 
                          if (head $ dropWhile (== ' ') fileContents') == '>' then 
                            let secondLine = head $ tail $ lines fileContents'
                                numSpaces = length $ elemIndices ' ' secondLine
                                -- numWords = length $ words secondLine
                            in
                            -- spaces between alphabet elements suggest fastc
                            if (numSpaces > 0) then 
                                let fastcData = getFastC firstOption fileContents' firstFile
                                    fastcCharInfo = getFastcCharInfo fastcData firstFile isPrealigned' tcmPair 
                                in
                                trace ("\tTrying to parse " ++ firstFile ++ " as fastc") 
                                executeReadCommands ((fastcData, [fastcCharInfo]) : curData) curGraphs isPrealigned' tcmPair (tail argList)
                            else 
                                let fastaData = getFastA firstOption fileContents' firstFile
                                    fastaCharInfo = getFastaCharInfo fastaData firstFile firstOption isPrealigned' tcmPair 
                                in
                                trace ("\tTrying to parse " ++ firstFile ++ " as fasta")
                                executeReadCommands ((fastaData, [fastaCharInfo]) : curData) curGraphs isPrealigned' tcmPair (tail argList)
                        
                    else errorWithoutStackTrace ("Cannot determine file type for " ++ firstOption ++ " need to prepend type")
                -- fasta
                else if (firstOption `elem` ["fasta", "nucleotide", "aminoacid"]) then 
                    let fastaData = getFastA firstOption fileContents firstFile
                        fastaCharInfo = getFastaCharInfo fastaData firstFile firstOption isPrealigned' tcmPair 
                    in
                    executeReadCommands ((fastaData, [fastaCharInfo]) : curData) curGraphs isPrealigned' tcmPair (tail argList)
                -- fastc
                else if (firstOption `elem` ["fastc", "custom_alphabet"])  then 
                    let fastcData = getFastC firstOption fileContents firstFile
                        fastcCharInfo = getFastcCharInfo fastcData firstFile isPrealigned' tcmPair 
                    in
                    executeReadCommands ((fastcData, [fastcCharInfo]) : curData) curGraphs isPrealigned' tcmPair (tail argList)
                -- tnt
                else if firstOption == "tnt" then
                    let tntData = TNT.getTNTData fileContents firstFile
                    in
                    executeReadCommands (tntData : curData) curGraphs isPrealigned' tcmPair (tail argList)
                else if firstOption == "prealigned" then executeReadCommands curData curGraphs isPrealigned' tcmPair (tail argList)
                -- FENEwick
                else if (firstOption `elem` ["newick" , "enewick", "fenewick"])  then 
                    let thisGraphList = getFENewickGraph fileContents
                        hasLoops = filter (== True) $ fmap B.hasLoop thisGraphList 
                        hasCycles = filter (== True) $ fmap GFU.cyclic thisGraphList
                    in 
                    if (not $ null hasLoops) then errorWithoutStackTrace ("Input graphin " ++ firstFile ++ "  has loops/self-edges")
                    else if (not $ null hasCycles) then errorWithoutStackTrace ("Input graph in " ++ firstFile ++ " has at least one cycle")
                    else executeReadCommands curData (thisGraphList ++ curGraphs) isPrealigned' tcmPair (tail argList)
                else errorWithoutStackTrace ("\n\n'Read' command error: option " ++ firstOption ++ " not recognized/implemented")
                

-- | Read arg list allowable modifiers in read
readArgList :: [String]
readArgList = ["tcm", "prealigned", "nucleotide", "aminoacid", "custom_alphabet", "fasta", "fastc", "tnt", "csv", 
    "dot", "newick" , "enewick", "fenewick"]

-- | getReadArgs processes arguments ofr the 'read' command
-- should allow mulitple files and gracefully error check
-- also allows tcm file specification (limit 1 tcm per command?)
-- as fasta, fastc, tnt, tcm, prealigned
getReadArgs :: String -> [(String, String)] -> [Argument]
getReadArgs fullCommand argList = 
    if null argList then []
    else 
        let (firstPart, secondPart) = head argList 
        in
        -- plain file name with no modifier
        if (null firstPart) then
            if (head secondPart == '"') || (last secondPart == '"') then (firstPart, (init $ tail secondPart)) : getReadArgs fullCommand (tail argList)
            else errorWithoutStackTrace ("\n\n'Read' command error (*) '" ++ (secondPart) ++"' : Need to specify filename in double quotes") 
        -- Change to allowed modifiers
        else if ((fmap toLower firstPart) `notElem` readArgList) then
            errorWithoutStackTrace ("\n\n'Read' command error (**): " ++ fullCommand ++ " contains unrecognized option '" ++ firstPart ++ "'")
        else if ((null secondPart) && (firstPart == "prealigned"))  then
            (firstPart, []) : getReadArgs fullCommand (tail argList)
        else (firstPart, (init $ tail secondPart)) : getReadArgs fullCommand (tail argList)

-- | getAlphabet takse a list of short-text lists and returns alphabet as list of short-text
-- filters out '?' '[' and ']' adds in '-' for indel Gap
getAlphabet :: [String] -> [ST.ShortText] -> [ST.ShortText] 
getAlphabet curList inList =
    let notAlphElement = fmap ST.fromString ["?", "[", "]"]
    in
    if null inList then 
        filter (`notElem` notAlphElement) $ fmap ST.fromString $ (sort curList) `union` ["-"]
    else 
        let firstChars = fmap (:[]) $ nub $ ST.toString $ head inList 
        in
        getAlphabet (firstChars `union` curList) (tail inList)


-- | generateDefaultMatrix takes an alphabet and generates cost matrix (assuming '-'
--   in already)
generateDefaultMatrix :: [ST.ShortText] -> Int -> [[Int]]
generateDefaultMatrix inAlph rowCount =
    if null inAlph then []
    else if rowCount == length inAlph then []
    else 
        let firstPart = replicate rowCount 1
            thirdPart = replicate ((length inAlph) - rowCount - 1) 1
        in
        (firstPart ++ [0] ++ thirdPart) : generateDefaultMatrix inAlph (rowCount + 1)

-- | getFastaCharInfo get alphabet , names etc from processed fasta data
-- this doesn't separate ambiguities from elements--processed later
-- need to read in TCM or default
getFastaCharInfo :: [TermData] -> String -> String -> Bool -> ([ST.ShortText], [[Int]]) -> CharInfo
getFastaCharInfo inData dataName dataType isPrealigned localTCM = 
    if null inData then error "Empty inData in getFastaCharInfo"
    else 
        let nucleotideAlphabet = fmap ST.fromString ["A","C","G","T","U","R","Y","S","W","K","M","B","D","H","V","N","?","-"]
            aminoAcidAlphabet  = fmap ST.fromString ["A","B","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","X","Y","Z", "-","?"]
            --onlyInNucleotides = [ST.fromString "U"]
            --onlyInAminoAcids = fmap ST.fromString ["E","F","I","L","P","Q","X","Z"]
            sequenceData = getAlphabet [] $ concat $ fmap snd inData
            seqType = if dataType == "nucleotide" then NucSeq
                      else if dataType == "aminoacid" then AminoSeq
                      else if dataType == "custom_alphabet" then GenSeq
                      else if (sequenceData `intersect` nucleotideAlphabet == sequenceData) then trace ("Assuming file " ++ dataName 
                          ++ " is nucleotide data. Specify `aminoacid' filetype if this is incorrect.") NucSeq
                      else if (sequenceData `intersect` aminoAcidAlphabet == sequenceData) then trace ("Assuming file " ++ dataName 
                          ++ " is amino acid data. Specify `nucleotide' filetype if this is incorrect.") AminoSeq
                      -- can fit in byte with reasonable pre-calculation costs
                      else if (length sequenceData) < 9 then SmallAlphSeq
                      else GenSeq
            seqAlphabet = if seqType == NucSeq then fmap ST.fromString ["A","C","G","T","-"]
                       else if seqType == AminoSeq then fmap ST.fromString ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y", "-"]
                       else if seqType == Binary then fmap ST.fromString ["0","1"]
                       else fst localTCM
            defaultGenSeqCharInfo = CharInfo {
                                       charType = seqType
                                     , activity = True
                                     , weight = 1.0
                                     , costMatrix = if localTCM == ([],[]) then SM.fromLists $ generateDefaultMatrix seqAlphabet 0
                                                    else SM.fromLists $ snd localTCM
                                     , name = T.pack ((filter (/= ' ') dataName) ++ ":0")
                                     , alphabet = if localTCM == ([],[]) then seqAlphabet
                                                  else fst localTCM
                                     , prealigned = isPrealigned
                                     }
        in
        trace ("Warning: no tcm file specified for use with fasta file : " ++ dataName ++ ". Using default, all 1 diagonal 0 cost matrix.")
        defaultGenSeqCharInfo

-- | getSequenceAphabet take a list of ShortText with inform ation and accumulatiors
-- For both nonadditive and additve looks for [] to denote ambiguity and splits states 
--  if splits on spaces if there are spaces (within []) (ala fastc or multicharacter states)
--  else if no spaces 
--    if non-additive then each symbol is split out as an alphabet element -- as in TNT
--    if is additive splits on '-' to denote range
-- rescales (integerizes later) additive characters with decimal places to an integer type rep
-- for additive charcaters if states are not nummerical then throuws an error
getSequenceAphabet :: [ST.ShortText]  -> [ST.ShortText] -> [ST.ShortText]
getSequenceAphabet newAlph inStates = 
    if null inStates then 
        -- removes indel gap from alphabet if present and then (re) adds at end
        (filter (/= (ST.singleton '-')) $ sort $ nub newAlph) ++ [ST.singleton '-']
    else
        let firstState = ST.toString $ head inStates
        in
        if (head firstState) /= '['  then 
            if (firstState `elem` ["?"]) then getSequenceAphabet  newAlph (tail inStates)
            else getSequenceAphabet ((head inStates) : newAlph) (tail inStates)
        else -- ambiguity
            let newAmbigStates  = fmap ST.fromString $ words $ filter (`notElem` ['[',']']) firstState
            in
            getSequenceAphabet (newAmbigStates ++ newAlph) (tail inStates) 




-- | getFastcCharInfo get alphabet , names etc from processed fasta data
-- this doesn't separate ambiguities from elements--processed later
-- need to read in TCM or default
getFastcCharInfo :: [TermData] -> String -> Bool -> ([ST.ShortText], [[Int]]) -> CharInfo
getFastcCharInfo inData dataName isPrealigned localTCM = 
    if null inData then error "Empty inData in getFastcCharInfo"
    else 
        --if null $ fst localTCM then errorWithoutStackTrace ("Must specify a tcm file with fastc data for fie : " ++ dataName)
        let thisAlphabet = if (not $ null $ fst localTCM) then fst localTCM
                           else getSequenceAphabet [] $ concat $ fmap snd inData
            inMatrix = if (not $ null $ fst localTCM) then SM.fromLists $ snd localTCM
                       else SM.fromLists $ generateDefaultMatrix thisAlphabet 0
            defaultGenSeqCharInfo = CharInfo {
                                       charType = if (length $ thisAlphabet) < 9 then SmallAlphSeq
                                                     else GenSeq
                                     , activity = True
                                     , weight = 1.0
                                     , costMatrix = inMatrix
                                     , name = T.pack ((filter (/= ' ') dataName) ++ ":0")
                                     , alphabet = thisAlphabet
                                     , prealigned = isPrealigned
                                     }
        in
        trace ("Warning: no tcm file specified for use with fastc file : " ++ dataName ++ ". Using default, all 1 diagonal 0 cost matrix.")
        defaultGenSeqCharInfo

-- | getSequenceAphabet takes a list of ShortText and returns the alp[habet and adds '-' if not present  

-- | getFastA processes fasta file 
-- assumes single character alphabet
-- deletes '-' (unless "prealigned"), and spaces
getFastA :: String -> String -> String -> [TermData] 
getFastA modifier fileContents' fileName  =
    if null fileContents' then errorWithoutStackTrace ("\n\n'Read' command error: empty file")
    else 
        -- removes ';' comments   
        let fileContents =  unlines $ filter (not.null) $ fmap (takeWhile (/= ';')) $ lines fileContents'
        in 
        if (head fileContents) /= '>' then errorWithoutStackTrace ("\n\n'Read' command error: fasta file must start with '>'")
        else 
            let terminalSplits = T.split (=='>') $ T.pack fileContents 
                pairData =  getRawDataPairsFastA modifier (tail terminalSplits)
                (hasDupTerminals, dupList) = DT.checkDuplicatedTerminals pairData
            in
            -- tail because initial split will an empty text
            if hasDupTerminals then errorWithoutStackTrace ("\tInput file " ++ fileName ++ " has duplicate terminals: " ++ show dupList)
            else pairData 
        
       
        

-- | getRawDataPairsFastA takes splits of Text and returns terminalName, Data pairs--minimal error checking
getRawDataPairsFastA :: String -> [T.Text] -> [TermData]
getRawDataPairsFastA modifier inTextList =
    if null inTextList then []
    else 
        let firstText = head inTextList
            firstName = T.takeWhile (/= '$') $ T.takeWhile (/= ';') $ head $ T.lines firstText
            firstData = T.filter (/= ' ') $ T.toUpper $ T.concat $ tail $ T.lines firstText
            firstDataNoGaps = T.filter (/= '-') firstData
            firtDataSTList = fmap ST.fromText $ fmap T.toStrict $ T.chunksOf 1 firstData
            firstDataNoGapsSTList = fmap ST.fromText $ fmap T.toStrict $ T.chunksOf 1 firstDataNoGaps
        in
        --trace (T.unpack firstName ++ "\n"  ++ T.unpack firstData) (
        --trace ("FA " ++ (show $ length firstDataNoGapsSTList)) (
        if modifier == "prealigned" then (firstName, firtDataSTList) : getRawDataPairsFastA modifier (tail inTextList)
        else (firstName, firstDataNoGapsSTList) : getRawDataPairsFastA modifier (tail inTextList)
        --)
        
-- | getFastC processes fasta file 
-- assumes spaces between alphabet elements
-- deletes '-' (unless "prealigned")
-- NEED TO ADD AMBIGUITY
getFastC :: String -> String -> String -> [TermData] 
getFastC modifier fileContents' fileName =
    if null fileContents' then errorWithoutStackTrace ("\n\n'Read' command error: empty file")
    else 
        -- removes ';' comments   
        let fileContents =  unlines $ filter (not.null) $ fmap (takeWhile (/= ';')) $ lines fileContents'
        in 
        if (head fileContents) /= '>' then errorWithoutStackTrace ("\n\n'Read' command error: fasta file must start with '>'")
        else 
            let terminalSplits = T.split (=='>') $ T.pack fileContents 
                pairData = recodeFASTCAmbiguities fileName $ getRawDataPairsFastC modifier (tail terminalSplits)
                (hasDupTerminals, dupList) = DT.checkDuplicatedTerminals pairData
            in
            -- tail because initial split will an empty text
            if hasDupTerminals then errorWithoutStackTrace ("\tInput file " ++ fileName ++ " has duplicate terminals: " ++ show dupList)
            else pairData
            
-- | recodeFASTCAmbiguities take list of TermData and scans for ambiguous groups staring with '['' and ending with ']
recodeFASTCAmbiguities :: String -> [TermData] -> [TermData] 
recodeFASTCAmbiguities fileName inData =
    if null inData then []
    else 
        let (firstName, firstData) = head inData
            newData = concatAmbig fileName firstData
        in 
        (firstName, newData) : recodeFASTCAmbiguities fileName (tail inData)

-- | concatAmbig takes a list of ShortText and concatanates ambiguyous states '['X Y Z...']' into a
-- single Short Tex for later processing
concatAmbig :: String -> [ST.ShortText] -> [ST.ShortText] 
concatAmbig fileName inList = 
    if null inList then []
    else 
        let firstGroup = ST.toString $ head inList 
        in
        -- not ambiguity group
        -- trace (firstGroup ++ show inList) (
        if null firstGroup then concatAmbig fileName (tail inList)
        else if head firstGroup /= '[' then (head inList) : concatAmbig fileName (tail inList)
        else 
            let ambiguityGroup = (head inList) : getRestAmbiguityGroup fileName (tail inList)
            in
            --trace (show ambiguityGroup) 
            (ST.concat ambiguityGroup) : concatAmbig fileName (drop (length ambiguityGroup) inList)
            --)

-- | getRestAmbiguityGroup takes a list of ShorText and keeps added them until one is found with ']'
getRestAmbiguityGroup :: String -> [ST.ShortText] -> [ST.ShortText]
getRestAmbiguityGroup fileName inList = 
    if null inList then errorWithoutStackTrace ("\n\n'Read' command error: fastc file " ++ fileName ++ " with unterminated ambiguity specification ']'")
    else 
        let firstGroup = ST.toString $ head inList
        in
        if ']' `notElem` firstGroup then (ST.cons ' ' $ head inList) : getRestAmbiguityGroup fileName (tail inList)
        else [ST.cons ' ' $ head inList]

-- | getRawDataPairsFastA takes splits of Text and returns terminalName, Data pairs--minimal error checking
-- this splits on spaces in sequences
getRawDataPairsFastC :: String -> [T.Text] -> [TermData]
getRawDataPairsFastC modifier inTextList =
    if null inTextList then []
    else 
        let firstText = head inTextList
            firstName = T.takeWhile (/= '$') $ T.takeWhile (/= ';') $ head $ T.lines firstText
            firstData = T.split (== ' ') $ T.concat $ tail $ T.lines firstText
            firstDataNoGaps = fmap (T.filter (/= '-')) firstData
        in
        -- trace (T.unpack firstName ++ "\n"  ++ (T.unpack $ T.intercalate (T.pack " ") firstData)) (
        if modifier == "prealigned" then (firstName, fmap ST.fromText  $ fmap T.toStrict firstData) : getRawDataPairsFastC modifier (tail inTextList)
        else (firstName, fmap ST.fromText $ fmap T.toStrict firstDataNoGaps) : getRawDataPairsFastC modifier (tail inTextList)
        --
        
-- | getFENewickGraph takes graph contents and returns local graph format
-- could be mulitple input graphs
getFENewickGraph :: String -> [LG.Gr T.Text Double] 
getFENewickGraph fileString = 
    -- trace (fileString)
    LG.getFENLocal (T.filter (/= '\n') $ T.strip $ T.pack fileString) 

-- | processTCMContents
processTCMContents :: String -> String -> ([ST.ShortText], [[Int]])
processTCMContents inContents fileName = 
    if null inContents then errorWithoutStackTrace ("\n\n'Read' 'tcm' command error: empty tcmfile `" ++ fileName)
    else 
        --First line is alphabet
        let tcmLines = lines inContents
            localAlphabet =  fmap ST.pack $ words $ head tcmLines
            -- to account for '-'
            numElements = 1 + (length localAlphabet)
            costMatrixStrings = fmap words $ tail tcmLines
            localCostMatrix = filter (/= []) $ fmap (fmap (GU.stringToInt fileName)) costMatrixStrings
            numLines = length localCostMatrix
            lengthRowList = fmap length localCostMatrix
            rowsCorrectLength = foldl' (&&) True $ fmap (==  numElements) lengthRowList
        in
        if (not $ null (filter (== ST.singleton '-') localAlphabet)) then errorWithoutStackTrace ("\n\n'Read' 'tcm' file format error: '-' (InDel/Gap) should not be specifid as an alphabet element.  It is added in automatically and assumed to be the last row and column of tcm matrix")
        else if numLines /= numElements then errorWithoutStackTrace ("\n\n'Read' 'tcm' file format error: incorrect number of lines in matrix from " 
            ++ fileName ++ " there should be (including '-') " ++ show numElements ++ " and there are " ++ show numLines)
        else if (not rowsCorrectLength) then errorWithoutStackTrace ("\n\n'Read' 'tcm' file format error: incorrect lines length(s) in matrix from " 
            ++ fileName ++ " there should be (including '-') " ++ show numElements)
        else
            -- trace (show localAlphabet ++ "\n" ++ show localCostMatrix) 
            (localAlphabet ++ [ST.singleton '-'], localCostMatrix)