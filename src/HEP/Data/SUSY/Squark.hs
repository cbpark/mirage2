{-# LANGUAGE RecordWildCards #-}

module HEP.Data.SUSY.Squark where

import HEP.Data.Constants       (mZ2, mb, mt, sinThetaW2)
import HEP.Data.Kinematics      (Mass (..), massSq)
import HEP.Data.SUSY.Parameters

mSquark :: Mass    -- ^ m_{~q_R}
        -> Mass    -- ^ quark mass
        -> Double  -- ^ T3(q)
        -> Double  -- ^ R(q)
        -> Double  -- ^ electric charge
        -> ModelParams
        -> Double  -- ^ mu
        -> (Mass, Mass)
mSquark (Mass mqR) mq t3 rq qq ModelParams {..} mu =
    (Mass (sqrt m1Sq), Mass (sqrt m2Sq))
  where
    mQ3 = _M0 * _Q  _c
    mQ3Sq = mQ3 * mQ3
    mqRSq = mqR * mqR
    mqSq  = massSq mq
    cos2b = cos2Beta _tanb
    aq = _M0

    term1 = mQ3Sq + mqRSq + 2.0 * mqSq + mZ2 * cos2b * t3
    term2 = sqrt $ (mQ3Sq - mqRSq + mZ2 * cos2b * (t3 - 2.0 * qq * sinThetaW2)) ** 2
            + 4.0 * mqSq * (aq - mu * rq) ** 2
    m1Sq = abs $ 0.5 * (term1 - term2)
    m2Sq = abs $ 0.5 * (term1 + term2)

mStop :: ModelParams -> Double -> (Mass, Mass)
mStop m@ModelParams {..} = mSquark mtR mt 0.5 (1.0 / _tanb) (2.0 / 3) m
  where mtR = Mass $ _tR _c * _M0

mSbottom :: ModelParams -> Double -> (Mass, Mass)
mSbottom m@ModelParams {..} = mSquark mbR mb (-0.5) _tanb (-1.0 / 3) m
  where mbR = Mass $ _bR _c * _M0

getMSUSY :: ModelParams -> Double -> Double
getMSUSY mpoint mu = sqrt $ mt1 * mt2
  where (Mass mt1, Mass mt2) = mStop mpoint mu
