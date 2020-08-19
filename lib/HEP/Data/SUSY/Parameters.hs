module HEP.Data.SUSY.Parameters where

-- data ModelParams = ModelParams { _M0   :: Double
--                                , _c    :: ModularWeights
--                                , _tanb :: Double
--                                } deriving Show

data ModularWeights = ModularWeights { _cHu :: Double
                                     , _cHd :: Double
                                     , _cQ  :: Double
                                     , _ctR :: Double
                                     , _cbR :: Double
                                     , _cL  :: Double
                                     , _ceR :: Double
                                     } deriving Show

-- data HuHdLoopFactor = HuHdLoopFactor { _kHu :: Double
--                                      , _kHd :: Double
--                                      } deriving Show

mStar :: Double
mStar = 1e+3

getMu :: ModularWeights
      -> Double  -- ^ M_0
      -> Double
getMu ModularWeights { _cHd = cHd } m0 =
    mStar * (mStar / m0) ** (7.0 / 12) / cHd ** (1.0 / 3)

cos2Beta :: Double -> Double
cos2Beta tanb = (1.0 - tanbSq) / (1.0 + tanbSq)
  where tanbSq = tanb * tanb

cosBeta :: Double -> Double
cosBeta tanb = 1.0 / sqrt (1 + tanb * tanb)

sinBeta :: Double -> Double
sinBeta tanb = tanb / sqrt (1 + tanb * tanb)
