module Analysis.Parsimony.Dynamic.DirectOptimization.Pairwise.DynamicCharacter2 where

import           Data.BitVector.LittleEndian
import qualified Data.Vector          as  V
import qualified Data.Vector.Storable as SV
import qualified Data.Vector.Unboxed  as UV
import           Data.Word
import           Foreign.C.Types

type SlimDynamicCharacter = (SV.Vector CUInt    , SV.Vector CUInt    , SV.Vector CUInt    )
type WideDynamicCharacter = (UV.Vector Word64   , UV.Vector Word64   , UV.Vector Word64   )
type HugeDynamicCharacter = ( V.Vector BitVector,  V.Vector BitVector,  V.Vector BitVector)
