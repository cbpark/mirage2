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
        -> ModularWeights
        -> Double  -- ^ tan(beta)
        -> Double  -- ^ M_0
        -> Double  -- ^ mu
        -> (Mass, Mass)
mSquark (Mass mqR) mq t3 rq qq ModularWeights {..} tanb m0 mu =
    (Mass (sqrt m1Sq), Mass (sqrt m2Sq))
  where
    mQ3 = m0 * _cQ
    mQ3Sq = mQ3 * mQ3
    mqRSq = mqR * mqR
    mqSq  = massSq mq
    cos2b = cos2Beta tanb
    aq = m0

    term1 = mQ3Sq + mqRSq + 2.0 * mqSq + mZ2 * cos2b * t3
    term2 = sqrt $
            (mQ3Sq - mqRSq + mZ2 * cos2b * (t3 - 2.0 * qq * sinThetaW2)) ** 2
            + 4.0 * mqSq * (aq - mu * rq) ** 2
    m1Sq = abs $ 0.5 * (term1 - term2)
    m2Sq = abs $ 0.5 * (term1 + term2)

mStop :: ModularWeights -> Double -> Double -> Double -> (Mass, Mass)
mStop cs@ModularWeights {..} tanb m0 =
    mSquark mtR mt 0.5 (1.0 / tanb) (2.0 / 3) cs tanb m0
  where mtR = Mass $ _ctR * m0

mSbottom :: ModularWeights -> Double -> Double -> Double -> (Mass, Mass)
mSbottom cs@ModularWeights {..} tanb m0 =
    mSquark mbR mb (-0.5) tanb (-1.0 / 3) cs tanb m0
  where mbR = Mass $ _cbR * m0

getMSUSY :: ModularWeights -> Double -> Double -> Double -> Double
getMSUSY cs tanb m0 mu = sqrt $ mt1 * mt2
  where (Mass mt1, Mass mt2) = mStop cs tanb m0 mu
