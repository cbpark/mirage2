{-# LANGUAGE RecordWildCards #-}

module HEP.Data.SUSY.Higgs where

import HEP.Data.Constants       (mZ2, mb, mt2, pi2, vEW, vEW2)
import HEP.Data.Kinematics      (Mass (..))
import HEP.Data.SUSY.Parameters
import HEP.Data.SUSY.Squark     (getMSUSY)

getMu :: ModelParams
      -> Double  -- ^ m_*
      -> Double
getMu ModelParams {..} mStar =
    mStar * (mStar / _M0) ** (7.0 / 12) / _Hd _c ** (1.0 / 3)

mHiggs :: ModelParams -> Double -> Double -> Maybe Mass
mHiggs m@ModelParams {..} mu as | mhSq <= 0 = Nothing
                                | otherwise = Just $ Mass (sqrt mhSq)
  where
    mhSq = mZ2 * cos2b * cos2b
           + 3.0 / (4 * pi2) * mt2 * mt2 / vEW2 * termT
           - yb2 * yb2 * vEW2 * loopFac * mu4 / (mSUSY2 * mSUSY2) * termB

    cos2b = cos2Beta _tanb

    mSUSY = getMSUSY m mu
    mSUSY2 = mSUSY * mSUSY
    loopT = log (mSUSY2 / mt2)

    aT = _M0 - mu / _tanb
    aT2 = aT * aT
    xT = 2.0 * aT2 / mSUSY2 * (1.0 - aT2 / (12.0 * mSUSY2))
    loopFac = 1.0 / (16 * pi2)

    termT = 0.5 * xT + loopT
            + loopFac * (3.0 / 2 * mt2 / vEW2 - 32 * pi * as)
            * (xT * loopT + loopT * loopT)

    cosb = cosBeta _tanb
    mu4 = mu ** 4

    yb = getMass mb / (vEW * cosb)
    yb2 = yb * yb
    termB = 1.0 + loopFac * loopT * (9.0 * yb2 - 5 * mt2 / vEW2 - 64 * pi * as)
