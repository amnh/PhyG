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
import qualified Data.Graph.Inductive.Basic as B
import qualified GeneralUtilities as GU
import qualified FastAC as FAC


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
executeReadCommands :: [RawData] -> [SimpleGraph] -> Bool -> ([ST.ShortText], [[Int]], Double) -> [Argument] -> IO ([RawData], [SimpleGraph])
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
            else    
                if not $ null firstOption then hPutStrLn stderr ("Reading " ++ firstFile ++ " with option " ++ firstOption)
                else hPutStrLn stderr ("Reading " ++ firstFile ++ " with no options")

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
                if (null fileContents) then errorWithoutStackTrace ("Error: Input file " ++ firstFile ++ " is empty")
                else 
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
                                    let fastcData = FAC.getFastC firstOption fileContents firstFile
                                        fastcCharInfo = FAC.getFastcCharInfo fastcData firstFile isPrealigned' tcmPair 
                                    in
                                    trace ("\tTrying to parse " ++ firstFile ++ " as fastc") 
                                    executeReadCommands ((fastcData, [fastcCharInfo]) : curData) curGraphs isPrealigned' tcmPair (tail argList)
                                else 
                                    let fastaData = FAC.getFastA  firstOption fileContents firstFile
                                        fastaCharInfo = FAC.getFastaCharInfo fastaData firstFile firstOption isPrealigned' tcmPair 
                                    in
                                    trace ("\tTrying to parse " ++ firstFile ++ " as fasta")
                                    executeReadCommands ((fastaData, [fastaCharInfo]) : curData) curGraphs isPrealigned' tcmPair (tail argList)
                            
                        else errorWithoutStackTrace ("Cannot determine file type for " ++ firstOption ++ " need to prepend type")
                    -- fasta
                    else if (firstOption `elem` ["fasta", "nucleotide", "aminoacid"]) then 
                        let fastaData = FAC.getFastA  firstOption fileContents firstFile
                            fastaCharInfo = FAC.getFastaCharInfo fastaData firstFile firstOption isPrealigned' tcmPair 
                        in
                        executeReadCommands ((fastaData, [fastaCharInfo]) : curData) curGraphs isPrealigned' tcmPair (tail argList)
                    -- fastc
                    else if (firstOption `elem` ["fastc", "custom_alphabet"])  then 
                        let fastcData = FAC.getFastC firstOption fileContents firstFile
                            fastcCharInfo = FAC.getFastcCharInfo fastcData firstFile isPrealigned' tcmPair 
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

-- | getFENewickGraph takes graph contents and returns local graph format
-- could be mulitple input graphs
getFENewickGraph :: String -> [LG.Gr T.Text Double] 
getFENewickGraph fileString = 
    -- trace (fileString)
    LG.getFENLocal (T.filter (/= '\n') $ T.strip $ T.pack fileString) 

-- | processTCMContents
-- functionality added to integerize tcm and add weight factor to allow for decimal tcm values
processTCMContents :: String -> String -> ([ST.ShortText], [[Int]], Double)
processTCMContents inContents fileName = 
    if null inContents then errorWithoutStackTrace ("\n\n'Read' 'tcm' command error: empty tcmfile `" ++ fileName)
    else 
        --First line is alphabet
        let tcmLines = lines inContents
            localAlphabet =  fmap ST.pack $ words $ head tcmLines
            -- to account for '-'
            numElements = 1 + (length localAlphabet)
            costMatrixStrings = fmap words $ tail tcmLines
            -- localCostMatrix = filter (/= []) $ fmap (fmap (GU.stringToInt fileName)) costMatrixStrings
            (scaleFactor, localCostMatrix) = getCostMatrixAndScaleFactor fileName costMatrixStrings 
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
            --trace (show (scaleFactor, localCostMatrix)) 
            -- trace (show localAlphabet ++ "\n" ++ show localCostMatrix) 
            (localAlphabet ++ [ST.singleton '-'], localCostMatrix, scaleFactor)


-- | getCostMatrixAndScaleFactor takes [[String]] and returns cost matrix as
-- [[Int]] but if there are decimal values, a scalerFactor is determined and the 
-- cost matrix are integerized by multiplication by 1/scaleFactor
-- scaleFactor will alway be a factor of 10 to allow for easier 
-- compilation of tcm charcaters later (weights the same...)
getCostMatrixAndScaleFactor :: String -> [[String]] -> (Double, [[Int]])
getCostMatrixAndScaleFactor fileName inStringListList =
    if null inStringListList then error "Empty list in inStringListList"
    else 
        let maxDecimalPlaces = maximum $ fmap getDecimals $ concat inStringListList
            scaleFactor = if maxDecimalPlaces == 0 then 1.0
                          else (0.1) ** (fromIntegral maxDecimalPlaces)
            in
            if maxDecimalPlaces == 0 then (scaleFactor, filter (/= []) $ fmap (fmap (GU.stringToInt fileName)) inStringListList)
            else 
                let newStringListList = fmap (fmap (integerizeString maxDecimalPlaces)) inStringListList
                in (scaleFactor, filter (/= []) $ fmap (fmap (GU.stringToInt fileName)) newStringListList)
                

-- | getDecimals returns the number of decimal places in a string rep of a number
getDecimals :: String -> Int
getDecimals inString =
    if null inString then 0
    else 
        let decimalPart = dropWhile (/= '.') inString
        in
        if null decimalPart then 0
        else (length decimalPart) - 1

-- | integerizeString integerizes a String by removing the '.' and adding the number of '0's to pad out 
-- adds maxDecimalPlaces - the nymber in the String  
integerizeString :: Int -> String -> String
integerizeString maxDecimalPlaces inString = 
    if null inString then error "Null string in integerizeString"
    else 
        let decimalPart = dropWhile (/= '.') inString
        in
        if inString == "0" then inString
        else if null decimalPart then inString ++ (replicate maxDecimalPlaces '0')
        else (filter (/= '.') inString) ++ (replicate (maxDecimalPlaces - ((length decimalPart) - 1)) '0')
