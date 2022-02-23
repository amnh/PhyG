{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE Strict       #-}

module System.Timing
  ( CPUTime()
  , fromPicoseconds
  , fromMilliseconds
  , fromMicroseconds
  , toPicoseconds
  , toMicroseconds
  , toMilliseconds
  , timeOp
  ) where


import           Control.DeepSeq
import           Control.Monad.IO.Class
import           Data.Foldable
import           Numeric.Natural
import           System.CPUTime


-- | CPU time with picosecond resolution
newtype CPUTime = CPUTime Natural


instance NFData CPUTime where

    rnf (CPUTime !_) = ()


instance Show CPUTime where

    show (CPUTime x)
      | x < nSecond = let (q,_) = x `quotRem` 1       in fold [show q, ".", "???"                       , "ps" ]
      | x < μSecond = let (q,r) = x `quotRem` nSecond in fold [show q, ".", zeroPad 3 (r `div` 1       ), "ns" ]
      | x < mSecond = let (q,r) = x `quotRem` μSecond in fold [show q, ".", zeroPad 3 (r `div` nSecond ), "μs" ]
      | x <  second = let (q,r) = x `quotRem` mSecond in fold [show q, ".", zeroPad 3 (r `div` μSecond ), "ms" ]
      | x <  minute = let (q,r) = x `quotRem`  second in fold [show q, ".", zeroPad 3 (r `div` mSecond ), "s " ]
      | x <    hour = let (q,r) = x `quotRem`  minute in fold [show q, "m", zeroPad 2 (r `div`  second ), "sec"]
      | x <     day = let (q,r) = x `quotRem`    hour in fold [show q, "h", zeroPad 2 (r `div`  minute ), "min"]
      | otherwise   = let (q,r) = x `quotRem`     day in fold [show q, "d", zeroPad 2 (r `div`    hour ), "hrs"]
      where
        nSecond = 1000
        μSecond = 1000 * nSecond
        mSecond = 1000 * μSecond
        second  = 1000 * mSecond
        minute  = 60   *  second
        hour    = 60   *  minute
        day     = 24   *  hour


zeroPad :: Int -> Natural -> String
zeroPad k i = replicate (k - length shown) '0' <> shown
  where
    shown = show i


timeOp :: (MonadIO m, NFData a) => m a -> m (CPUTime, a)
timeOp ioa = do
    t1 <- liftIO getCPUTime
    a  <- force <$> ioa
    t2 <- liftIO getCPUTime
    let t = CPUTime . fromIntegral $ t2 - t1
    pure (t, a)


fromPicoseconds :: Natural -> CPUTime
fromPicoseconds = CPUTime


fromMicroseconds :: Natural -> CPUTime
fromMicroseconds = CPUTime . (*1000000)


fromMilliseconds :: Natural -> CPUTime
fromMilliseconds = CPUTime . (*1000000000)


toPicoseconds :: CPUTime -> Natural
toPicoseconds (CPUTime x) = x


toMicroseconds :: CPUTime -> Natural
toMicroseconds (CPUTime x) = x `div` 1000000


toMilliseconds :: CPUTime -> Natural
toMilliseconds (CPUTime x) = x `div` 1000000000
