{- |
Module      :  SymMatrix.hs
Description :  Progam to manipulate square symmetric lower diagonal matrices with diagnonal values
                as if they were normal matrices
Copyright   :  (c) 2020 Ward C. Wheeler, Division of Invertebrate Zoology, AMNH. All rights reserved.
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

{-# Language ImportQualifiedPost #-}

module SymMatrix
    ( empty
    , dim
    , fromLists
    , Matrix
    , SymMatrix.null
    , cols
    , rows
    , (!)
    , toLists
    , toRows
    , fromRows
    , toFullLists
    , getFullRow
    , isSymmetric
    , updateMatrix
    , unsafeUpdateMatrix
    , addMatrixRow
    , addMatrices
    , deleteRowsAndColumns
    , showMatrixNicely
    , SymMatrix.map
    , SymMatrix.flatten
    , getFullRowVect
    , SymMatrix.zipWith
    , SymMatrix.zip
    , combine
    , safeIndex
    , makeDefaultMatrix
    , getFullVects
    ) where

import Data.List           qualified as L
import Data.Sort           qualified as S
import Data.Vector         qualified as V
import Data.Vector.Generic qualified as G

-- | Matrix type as Vector of Vectors
type Matrix a = V.Vector (V.Vector a)

-- | functions for triples
fst3 :: (a,b,c) -> a
fst3 (d,_,_) = d

-- | empty matrix value from Vector
empty :: Matrix a
empty = G.empty

-- | dim returns dimmension (rows = cols)
dim :: (Eq a) => Matrix a -> (Int, Int)
dim inMatrix =
    if SymMatrix.null inMatrix then (0,0)
    else (V.length inMatrix, V.length inMatrix)

-- | rows returns number rows in Matrix (rows = cols)
rows :: (Eq a) => Matrix a -> Int
rows inM = fst $ dim inM

-- | cols returns number cols in Matrix (rows = cols)
cols :: (Eq a) => Matrix a -> Int
cols inM = fst $ dim inM

-- | null returns True of row number is 0
null :: (Eq a) => Matrix a -> Bool
null inMatrix = inMatrix == empty

-- | isSymmetric is true by defineition--when created error if not
isSymmetric :: (Eq a) => Matrix a -> Bool
isSymmetric inM =
    not (SymMatrix.null inM) || error "Null matrix in isSymmetric"

-- | fromRows creates a lower diagonal matrix (with diagonal)
fromRows :: (Eq a, Show a) => [V.Vector a] -> Matrix a
fromRows inVectList = fromLists $ fmap V.toList inVectList

-- | toRows converts a Matrix to a list of Vectors
-- unequal in length
toRows :: Matrix a -> [V.Vector a]
toRows = V.toList

-- | fromLists takes list of list of a and returns lower diagnoal (with diagonal)
-- matrix as Matrix (Vector of Vectors)
fromLists :: (Eq a, Show a) => [[a]] -> Matrix a
fromLists inListList =
    if L.null inListList then empty
    else
        let initialSquare = V.map V.fromList $ V.fromList inListList
            colsH = V.length $ V.head initialSquare
            rowsH = V.length initialSquare
        in
        if colsH /= rowsH then error ("Input matrix is not square " ++ (show (colsH, rowsH)) ++ " " ++ show inListList)
        else
            let indexPairs = cartProd [0..(rowsH - 1)] [0..(rowsH - 1)]
                sym = checkSymmetry initialSquare indexPairs
            in
            if not sym then error "Input matrix not symmetrical"
            else makeLowerDiag initialSquare 0 rowsH

-- | toLists takes a Matrix and returns a list of lists (not all same length)
toLists :: Eq a => Matrix a -> [[a]]
toLists inM =
    if SymMatrix.null inM then []
    else
        V.toList $ V.map V.toList inM

-- | toFullLists takes a Matrix and returns a list of lists of full length
-- square matrix
toFullLists :: Eq a => Matrix a -> [[a]]
toFullLists inM =
    if SymMatrix.null inM then []
    else
        fmap (getFullRow inM) [0..(rows inM - 1)]

-- | getFullRow returns a specific full row (is if matrix were square)
-- as a list
getFullRow :: Eq a => Matrix a -> Int -> [a]
getFullRow inM index =
    if SymMatrix.null inM then []
    else
        let firstPart = V.toList $ inM V.! index -- initial [0..index] of elements
            restMatrix = V.drop (index + 1) inM
            restByColumn = V.toList $ V.map (V.! index) restMatrix
        in
        firstPart ++ restByColumn

-- | getFullVects returns full size verctor of vector not lower triangular
getFullVects :: Eq a => Matrix a -> Matrix a
getFullVects inLV =
  if SymMatrix.null inLV then  SymMatrix.empty
  else
    V.fromList $ fmap (getFullRowVect inLV) [0..((SymMatrix.rows inLV) - 1)]     

-- | getFullRowVect reurns a specific full row (is if matrix were square)
-- as a Vector
getFullRowVect :: Eq a => Matrix a -> Int -> V.Vector a
getFullRowVect inM index =
    if SymMatrix.null inM then V.empty
    else
        let firstPart = inM V.! index -- initial [0..index] of elements
            restMatrix = V.drop (index + 1) inM
            restByColumn = V.map (V.! index) restMatrix
        in
        firstPart V.++ restByColumn

-- | indexing lower diag matrix
(!) :: Matrix a -> (Int, Int) -> a
(!) inM (iIndex,jIndex) =
    if iIndex > jIndex then
        (inM V.! iIndex) V.! jIndex
    else
        (inM V.! jIndex) V.! iIndex

-- | indexing lower diag matrix
safeIndex :: Matrix a -> (Int, Int) -> a
safeIndex inM (iIndex,jIndex) =
    if iIndex >= (length inM) then error ("First out of bounds " ++ show (iIndex, (length inM)))
    else if jIndex >= (length inM) then error ("Second out of bounds " ++ show (jIndex, (length inM)))
    else if iIndex > jIndex then
        (inM V.! iIndex) V.! jIndex
    else
        (inM V.! jIndex) V.! iIndex

-- | makeLowerDiag take a Vector of Vetors (Matrix) and returns a lower diagonal matrix
-- including diagonal
makeLowerDiag :: (Eq a) => Matrix a -> Int -> Int -> Matrix a
makeLowerDiag inM row numRows
  | SymMatrix.null inM = error "Input matrix is empty in makeLowerDiag"
  | row == numRows = V.empty
  | otherwise =
    let origRow = inM V.! row
        newRow = V.take (row + 1) origRow
    in  V.cons newRow (makeLowerDiag inM (row + 1) numRows)

-- | makeDefaultMatrix creates an all 1 wth diagonal 0 matrix of size n x n
makeDefaultMatrix :: Int -> Matrix Int 
makeDefaultMatrix n =
    let row = replicate n 1
        rowList = replicate n row
        initialMattrix = fromLists rowList
        updateIndexList  = [0..(n-1)]
        zeroList = replicate n 0
        updateList = zip3 updateIndexList updateIndexList zeroList
    in
    updateMatrix initialMattrix updateList


-- | cartesian product of two lists
cartProd :: [a] -> [a] -> [(a,a)]
cartProd xs ys = [(x,y) | x <- xs, y <- ys]

-- | checkSymmetry takes a Vector of VEctors of a
-- and checks for symmetry
checkSymmetry :: (Eq a, Show a) => Matrix a -> [(Int, Int)] -> Bool
checkSymmetry inVV pairList
  | SymMatrix.null inVV = error "Input matrix is empty"
  | L.null pairList = True
  | otherwise =
    let (iIndex, jIndex) = head pairList
        firstCheck = ((inVV V.! iIndex) V.! jIndex) == ((inVV V.! jIndex) V.! iIndex)
    in
    if firstCheck then checkSymmetry inVV (tail pairList)
    else error ("Matrix is not symmetrical:" ++ show (iIndex, jIndex) ++ "=>" ++ show ((inVV V.! iIndex) V.! jIndex) ++ " /= " ++ show ((inVV V.! jIndex) V.! iIndex))

-- | addMatrixRow add a row to existing matrix to extend Matrix dimension
-- used when adding HTU distances to existing distance matrix as Wagner tree is built
addMatrixRow :: Matrix a -> V.Vector a -> Matrix a
addMatrixRow inM newRow =
    if V.null newRow then inM
    else
        inM `V.snoc` newRow


-- | addMatrices  adds a Matrix to existing matrix to extend Matrix dimension
addMatrices :: (Eq a) => Matrix a -> Matrix a -> Matrix a
addMatrices inM newMatrix =
    if SymMatrix.null newMatrix then inM
    else
        inM V.++ newMatrix

-- | reIndexTriple taske (i,j,k) and returns (max i j, min i j, k)
reIndexTriple :: (Ord a) => (a, a, b) -> (a, a, b)
reIndexTriple trip@(iIndex, jIndex, value) =
    if iIndex > jIndex then trip
    else (jIndex, iIndex, value)

-- | updateMatrix takes a list of triples and update matrix
-- update all at once checking for bounds
-- could naively do each triple in turn, but would be alot of copying
updateMatrix :: (Show a, Ord a) => Matrix a -> [(Int, Int, a)] -> Matrix a
updateMatrix inM modList =
    if L.null modList then inM
    else
        let orderedTripleList = S.uniqueSort $ fmap reIndexTriple modList
            minRow = fst3 $ head orderedTripleList
            maxRow = fst3 $ last orderedTripleList
        in
        if minRow < 0 then error ("Update matrix out of bounds: " ++ show orderedTripleList)
        else if maxRow >= rows inM then error ("Update matrix out of bounds, row = " ++ show (rows inM) ++ " and trying to update row " ++ show maxRow)
        else
            let firstPart = V.unsafeTake minRow inM
                restPart  = V.unsafeDrop minRow inM
                modifiedRemainder = updateRows restPart orderedTripleList minRow
            in
            addMatrices firstPart modifiedRemainder

-- | unsafeUpdateMatrix unsafe version of updateMatrix
unsafeUpdateMatrix :: (Show a, Ord a) => Matrix a -> [(Int, Int, a)] -> Matrix a
unsafeUpdateMatrix inM modList =
    if L.null modList then inM
    else
        let orderedTripleList = S.uniqueSort $ fmap reIndexTriple modList
            minRow = fst3 $ head orderedTripleList
            firstPart = V.unsafeTake minRow inM
            restPart  = V.unsafeDrop minRow inM
            modifiedRemainder = updateRows restPart orderedTripleList minRow
        in
        addMatrices firstPart modifiedRemainder

-- | updateRows takes the section of the matrix containing rows that wil be modified
-- (some not) and modifes or copies rows and rerns a Matrix (vector of roow vectors)
updateRows :: (Show a, Eq a) => Matrix a -> [(Int, Int, a)] -> Int -> Matrix a
updateRows inM tripList currentRow
  | L.null tripList = inM
  | SymMatrix.null inM = error ("Matrix is empty and there are modifications that remain: " ++ show tripList)
  | otherwise =
    let (rowIndex, columnIndex, value) = L.head tripList
        firstOrigRow = V.head inM
    in
    if currentRow /= rowIndex then firstOrigRow `V.cons` updateRows (V.tail inM) tripList (currentRow + 1)
    else -- account for multiple modifications to same row
        let (newRow, newTripList) = modifyRow firstOrigRow columnIndex value currentRow (L.tail tripList)
        in
        -- This for debug--remove after test
        newRow `V.cons` updateRows (V.tail inM) newTripList (currentRow + 1)


-- | modifyRow takes an initial modification (column and value) and then checks to see if there are more modifications in that
-- row (rowNumber) in the remainder of the list of modifications, returning the new row and mod list as a pair
-- assumes that sorted triples sort by first, second, then third elements
modifyRow :: V.Vector a -> Int -> a -> Int -> [(Int, Int, a)] -> (V.Vector a, [(Int, Int, a)])
modifyRow inRow colIndex value rowNumber modList =
    if colIndex >= V.length inRow then error ("Column to modify is outside length of row " ++ show (rowNumber, colIndex))
    else
        let firstPart = V.unsafeTake colIndex inRow
            remainderPart = V.unsafeDrop (colIndex + 1) inRow
            newRow = firstPart V.++ (value `V.cons` remainderPart)
        in
        if L.null modList then (newRow, modList)
        else continueRow (firstPart `V.snoc` value) inRow (colIndex + 1) rowNumber modList


-- | continueRow continues to modify a row with multiple column modifcations
-- assumes that sorted triples sorted by first, second, then third elements
continueRow :: V.Vector a ->V.Vector a -> Int -> Int -> [(Int, Int, a)] -> (V.Vector a, [(Int, Int, a)])
continueRow partRow origRow colIndex rowNumber modList
  | colIndex == V.length origRow = (partRow, modList)
  | L.null modList = (partRow V.++ V.unsafeDrop colIndex origRow, modList)
  | otherwise =
    let (nextRowNumber, nextColIndex, nextValue) = L.head modList
    in
    if nextRowNumber /= rowNumber then (partRow V.++ V.unsafeDrop colIndex origRow, modList)
    else
        if nextColIndex /= colIndex then continueRow (partRow `V.snoc` (origRow V.! colIndex)) origRow (colIndex + 1) rowNumber modList
        else continueRow (partRow `V.snoc` nextValue) origRow (colIndex + 1) rowNumber (L.tail modList)

-- | makeNiceRow pretty preints a list
makeNiceRow :: (Show a) => V.Vector a -> String
makeNiceRow aVect =
  if V.null aVect then "\n"
  else
    show (V.head aVect) ++ " " ++  makeNiceRow (V.tail aVect)

-- | showNicely pretty prins matrix
showMatrixNicely :: (Show a, Eq a) => Matrix a -> String
showMatrixNicely inM =
  let mRows = rows inM
      mCols = cols inM
      niceRows = V.map makeNiceRow inM
  in
  ("Dimensions: :" ++ show mRows ++ " " ++ show mCols ++ "\n" ++ concat niceRows)

-- | deleteRowsAndColumns take a list of rows (and same index for columns)
-- to delete from Matrix. Uses lisyt to do in single pass
deleteRowsAndColumns :: (Show a, Eq a) => Matrix a -> [Int] -> Matrix a
deleteRowsAndColumns inM deleteList =
    if L.null deleteList then inM
    else deleteRC inM deleteList (rows inM) 0


-- | deleteRC takes matri delete list and counter to delte coumns and rows
deleteRC :: (Show a, Eq a) => Matrix a -> [Int] -> Int -> Int -> Matrix a
deleteRC inM deleteList origRows rowCounter =
    if rowCounter == origRows then empty
    else
        let firstRow = V.head inM
            toKeep = rowCounter `L.notElem` deleteList
            newRow = deleteColumn firstRow deleteList (rowCounter + 1) 0
        in
        if toKeep then newRow `V.cons` deleteRC (V.tail inM) deleteList origRows (rowCounter + 1)
        else deleteRC (V.tail inM) deleteList origRows (rowCounter + 1)

-- | deleteColumn takes a row of a matrix (lower diagnonal), its length,
-- a list of cilumns to delete and a column counter and creates a new row
deleteColumn :: (Show a, Eq a) => V.Vector a -> [Int] -> Int -> Int -> V.Vector a
deleteColumn origRow deleteList rowLength colCounter =
    if colCounter == rowLength then V.empty
    else
        let firstValue = V.head origRow
            toKeep = colCounter `L.notElem` deleteList
        in
        if toKeep then firstValue `V.cons` deleteColumn (V.tail origRow) deleteList rowLength (colCounter + 1)
        else deleteColumn (V.tail origRow) deleteList rowLength (colCounter + 1)

-- | map maps a function over the matrix returning a matrix of new type
map :: (Eq a) => (a->b) -> Matrix a -> Matrix b
map f m =
  if SymMatrix.null m then empty
  else V.map (V.map f) m

-- | flatten concats rows of matrix to make a single Vector
flatten :: Eq a => Matrix a -> V.Vector a
flatten m =
  if SymMatrix.null m then V.empty
  else
    let rowList = fmap (getFullRowVect m) [0..(rows m - 1)]
    in
    V.concat rowList

-- | zip takes two matrices and zips into a matrix of pairs
zip :: (Eq a, Eq b) => Matrix a -> Matrix b -> Matrix (a,b)
zip m1 m2 = 
  if dim m1 /= dim m2 then error ("Cannot zip matrices with unequal dimensions " ++ (show $ dim m1) ++ " " ++ (show $ dim m2))
  else if V.null m1 then V.empty 
  else 
    let m1r = V.head m1
        m2r = V.head m2
        newRow = V.zip m1r m2r
    in
    V.cons newRow (SymMatrix.zip (V.tail m1) (V.tail m2))

-- | zip takes two matrices and zips into a matrix using f
zipWith :: (Eq a, Eq b) => (a -> b -> c) -> Matrix a -> Matrix b -> Matrix c
zipWith f m1 m2 = 
  if dim m1 /= dim m2 then error ("Cannot zip matrices with unequal dimensions " ++ (show $ dim m1) ++ " " ++ (show $ dim m2))
  else if V.null m1 then V.empty 
  else 
    let m1r = V.head m1
        m2r = V.head m2
        newRow = V.zipWith f m1r m2r
    in
    V.cons newRow (SymMatrix.zipWith f (V.tail m1) (V.tail m2))
    --where g (a,b) = f a b

-- | combine takes an operator f (Enforcing Num as opposed to zipWith) and two matrices
-- applying f to each element of the two matrices M1 f M2
-- to create the output
combine :: Eq a => (a -> a -> a) -> Matrix a -> Matrix a -> Matrix a
combine f m1 m2 = 
  if SymMatrix.null m1 then error "Null matrix 1 in combine"
  else if SymMatrix.null m2 then error "Null matrix 2 in combine"
  else if dim m1 /= dim m2 then error ("Cannot combine matrices with unequal dimensions " ++ (show $ dim m1) ++ " " ++ (show $ dim m2))
  else 
    SymMatrix.map g $ SymMatrix.zip m1 m2  
    where g (a,b) = f a b
