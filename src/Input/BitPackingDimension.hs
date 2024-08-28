module BitPackingDimension where


data BitDimension
    = BitDim0
    | BitDim1
    | BitDim2
    | BitDim3
    | BitDim4
    | BitDim5
    | BitDim6
    | BitDimN
    deriving (Bounded, Eq, Ord)


instance Enum BitDimension where
    fromEnum BitDim0 = 0
    fromEnum BitDim1 = 1
    fromEnum BitDim2 = 2
    fromEnum BitDim3 = 3
    fromEnum BitDim4 = 4
    fromEnum BitDim5 = 5
    fromEnum BitDim6 = 6
    fromEnum BitDimN = maxBound


    toEnum i = case i `mod` 7 of
        0 -> BitDim0
        1 -> BitDim1
        2 -> BitDim2
        3 -> BitDim3
        4 -> BitDim4
        5 -> BitDim5
        6 -> BitDim6
        _ -> BitDimN


instance Show BitDimension where
    show d =
        let getExp BitDim0 = '⁰'
            getExp BitDim1 = '¹'
            getExp BitDim2 = '²'
            getExp BitDim3 = '³'
            getExp BitDim4 = '⁴'
            getExp BitDim5 = '⁵'
            getExp BitDim6 = '⁶'
            getExp BitDimN = 'ⁿ'
        in  '❬' : '2' : getExp d : "❭"


bitWidth :: BitDimension -> Maybe Word
bitWidth BitDimN = Nothing
bitWidth bitDimV = Just . (2 ^) $ fromEnum bitDimV


computeDimension :: Word -> BitDimension
computeDimension x
    | x > 2 ^ 6 = BitDimN
    | x > 2 ^ 5 = BitDim6
    | x > 2 ^ 4 = BitDim5
    | x > 2 ^ 3 = BitDim4
    | x > 2 ^ 2 = BitDim3
    | x > 2 ^ 1 = BitDim2
    | x > 2 ^ 0 = BitDim1
    | otherwise = BitDim0
