module DOFFI where

import Analysis.Parsimony.DirectOptimization.Pairwise
import Bio.Character.Encodable
import Data.Alphabet
import Data.TCM


wrapperPCG_DO_FFI :: V.Vector BV.BV -> V.Vector BV.BV -> V.Vector (V.Vector Int) -> (V.Vector BV.BV, Int)
wrapperPCG_DO_FFI = undefined