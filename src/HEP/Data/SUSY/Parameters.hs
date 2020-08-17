{-# LANGUAGE RecordWildCards #-}

module HEP.Data.SUSY.Parameters where

import HEP.Data.Constants
import HEP.Data.Kinematics (Mass (..))

data ModelParams = ModelParams { _M0   :: Double
                               , _c    :: ModularWeights
                               , _tanb :: Double
                               } deriving Show

data ModularWeights = ModularWeights { _Hu :: Double
                                     , _Hd :: Double
                                     , _Q  :: Double
                                     , _tR :: Double
                                     , _bR :: Double
                                     , _L  :: Double
                                     , _eR :: Double
                                     } deriving Show

mHiggs :: ModelParams -> Double -> Double -> Maybe Mass
mHiggs ModelParams {..} mu as | mhSq <= 0 = Nothing
                              | otherwise = Just $ Mass (sqrt mhSq)
  where
    mhSq = mZ2 * cos2b * cos2b + 3.0 / (4 * pi2) * mt2 * mt2 / vEW2 * termT

    cos2b = cos2Beta _tanb

    mSUSY2 = _M0 * _M0 / 2.0
    loopT = log (mSUSY2 / mt2)

    aT = _M0 - mu / _tanb
    aT2 = aT * aT
    xT = 2.0 * aT2 / mSUSY2 * (1.0 - aT2 / (12.0 * mSUSY2))

    termT = 0.5 * xT + loopT
            + (3.0 / 2 * mt2 / vEW2 - 32 * pi * as)
            * (xT * loopT + loopT * loopT) / (16.0 * pi2)

getMu :: ModelParams
      -> Double  -- ^ m_*
      -> Double
getMu ModelParams {..} mStar =
    mStar * (mStar / _M0) ** (7.0 / 12) / _Hd _c ** (1.0 / 3)

cos2Beta :: Double -> Double
cos2Beta tanb = (1.0 - tanbSq) / (1.0 + tanbSq)
  where tanbSq = tanb * tanb

cosBeta :: Double -> Double
cosBeta tanb = 1.0 / sqrt (1 + tanb * tanb)
