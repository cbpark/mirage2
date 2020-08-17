module HEP.Data.SUSY.Parameters where

import HEP.Data.Constants
import HEP.Data.Kinematics (Mass (..))

-- data ModelParams = ModelParam { _M0   :: Double
--                               , _cHu  :: Double
--                               , _cHd  :: Double
--                               , _cQ   :: Double
--                               , _ctR  :: Double
--                               , _cbR  :: Double
--                               , _cL   :: Double
--                               , _c
--                               , _tanb :: Double
--                               } deriving Show

-- data ModularWeight = ModularWeight

mHiggs :: Double -> Double -> Double -> Double -> Maybe Mass
mHiggs as tanb muV mSUSY | mhSq <= 0 = Nothing
                         | otherwise = Just $ Mass (sqrt mhSq)
  where
    mhSq = mZ2 * cos2b * cos2b + 3.0 / (4 * pi2) * mt2 * mt2 / vEW2 * termT

    cos2b = cos2Beta tanb
    loopT = log (mSUSY * mSUSY / mt2)
    aT2 = 2.0  -- aT = At / MSUSY, At ~ M0, MSUSY ~ M0 / sqrt(2), so aT^2 ~ 2.
    xT = 2.0 * aT2 * (1.0 - aT2 / 12.0)
    termT = 0.5 * xT + loopT
            + (3.0 / 2 * mt2 / vEW2 - 32 * pi * as)
            * (xT * loopT + loopT * loopT) / (16.0 * pi2)

cos2Beta :: Double -> Double
cos2Beta tanb = (1.0 - tanbSq) / (1.0 + tanbSq)
  where tanbSq = tanb * tanb

mu :: Double  -- ^ c_{H_d)
   -> Double  -- ^ m_*
   -> Double  -- ^ M_0
   -> Double
mu cHd mStar m0 = mStar * (mStar / m0) ** (7.0 / 12) / cHd ** (1.0 / 3)
