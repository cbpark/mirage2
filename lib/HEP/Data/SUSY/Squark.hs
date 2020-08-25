module HEP.Data.SUSY.Squark
    (
      mStop
    , mSbottom
    , getMSUSY
    , getMSUSY2
    , getM0FromStop
    ) where

import HEP.Data.Constants       (mZ2, mb, mt, sinThetaW2)
import HEP.Data.Kinematics      (Mass (..), massSq)
import HEP.Data.SUSY.Parameters (ModularWeights (..), cos2Beta, getMu)
import HEP.Data.Util            (riddersSolver)

mSquark :: Mass            -- ^ m_{~q_R}
        -> Mass            -- ^ quark mass
        -> Double          -- ^ T3(q) : either 1/2 or -1/2
        -> Double          -- ^ R(q) : either cot(beta) or tan(beta)
        -> Double          -- ^ electric charge
        -> ModularWeights
        -> Double          -- ^ tan(beta)
        -> Double          -- ^ M_0
        -> Double          -- ^ mu
        -> (Mass, Mass)
mSquark (Mass mqR) mq t3 rq qq ModularWeights { _cQ = cQ } tanb m0 mu =
    (Mass (sqrt m1Sq), Mass (sqrt m2Sq))
  where
    mQ3Sq = m0 * m0 * cQ
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
mStop cs@ModularWeights { _ctR = ctR } tanb m0 =
    mSquark mtR mt 0.5 (1.0 / tanb) (2.0 / 3) cs tanb m0
  where mtR = Mass $ sqrt ctR * m0

mSbottom :: ModularWeights -> Double -> Double -> Double -> (Mass, Mass)
mSbottom cs@ModularWeights { _cbR = cbR } tanb m0 =
    mSquark mbR mb (-0.5) tanb (-1.0 / 3) cs tanb m0
  where mbR = Mass $ sqrt cbR * m0

getMSUSY :: ModularWeights -> Double -> Double -> Double -> Double
getMSUSY cs tanb m0 mu = sqrt $ mst1 * mst2
  where (Mass mst1, Mass mst2) = mStop cs tanb m0 mu

getMSUSY2 :: ModularWeights -> Double -> Double -> Double -> Double
getMSUSY2 cs tanb m0 mu = 0.5 * (mst1 * mst1 + mst2 * mst2)
  where (Mass mst1, Mass mst2) = mStop cs tanb m0 mu

getM0FromStop :: ModularWeights
              -> Double            -- ^ m_*
              -> (Mass, Double)    -- ^ (mtMS, alpha_s(mtMS))
              -> Double            -- ^ m_{stop}
              -> (Double, Double)  -- ^ (xlow, xup)
              -> Double            -- ^ tan(beta)
              -> Maybe Double
getM0FromStop cs mStar mta3 mstop range tanb = riddersSolver stopF range
  where
    stopF m0 = let Mass mst1 = fst $ mStop cs tanb m0 (getMu cs mStar mta3 m0)
               in mstop - mst1
