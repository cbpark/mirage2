{-# LANGUAGE RecordWildCards #-}

module HEP.Data.SUSY.Higgs where

import HEP.Data.Constants       (mZ2, mtau, pi2, vEW, vEW2)
import HEP.Data.Kinematics      (Mass (..), massSq)
import HEP.Data.SUSY.Parameters
import HEP.Data.SUSY.Squark     (getMSUSY)

import Numeric.RootFinding

mStar :: Double
mStar = 1e+3

getMu :: ModularWeights
      -> Double  -- ^ M_0
      -> Double
getMu ModularWeights { _cHd = cHd } m0 =
    mStar * (mStar / m0) ** (7.0 / 12) / cHd ** (1.0 / 3)

-- | From <https://arxiv.org/abs/1112.3336 arXiv:1112.3336>
--
--   valid for moderate or large tan(beta) and
--   (m_{~t_2}^2 - m_{~t_1}^2) / (m_{~t_2}^2 + m_{~t_1}^2) < 0.5.
mHiggs :: ModularWeights
       -> (Mass, Mass)  -- ^ (mtMS, mbMS)
       -> Double        -- ^ alpha_s
       -> Double        -- ^ tan(beta)
       -> Double        -- ^ M_0
       -> Double
mHiggs cs@ModularWeights {..} (mtMS, mbMS) as tanb m0
    | mhSq <= 0 = 0
    | otherwise = sqrt mhSq
  where
    mhSq = mZ2 * cos2b * cos2b
           + 3.0 / (4 * pi2) * mt2 * mt2 / vEW2 * termT
           - yb2 * yb2 * vEW2 * loopFac * termB
           - ytau2 * ytau2 * vEW2 * loopFac / 3 * termTau

    [cos2b, cosb] = fmap ($ tanb) [cos2Beta, cosBeta]
    mt2 = massSq mtMS
    mu = getMu cs m0
    mu4 = mu ** 4

    mSUSY = getMSUSY cs tanb m0 mu
    mSUSY2 = mSUSY * mSUSY
    loopT = log (mSUSY2 / mt2)

    aT = m0 - mu / tanb
    aT2 = aT * aT
    xT = 2.0 * aT2 / mSUSY2 * (1.0 - aT2 / (12.0 * mSUSY2))
    loopFac = 1.0 / (16 * pi2)

    termT = 0.5 * xT + loopT
            + loopFac * (1.5 * mt2 / vEW2 - 32 * pi * as)
            * (xT * loopT + loopT * loopT)

    -- the contribution from sbottom
    yb = getMass mbMS / (vEW * cosb)
    yb2 = yb * yb
    termB = mu4 / (mSUSY2 * mSUSY2)
            * (1.0 + loopFac
               * loopT * (9.0 * yb2 - 5 * mt2 / vEW2 - 64 * pi * as))

    -- the contribution from stau
    ytau = getMass mtau / (vEW * cosb)
    ytau2 = ytau * ytau
    mStau = m0 * _cL
    termTau = mu4 / mStau ** 4

mHiggsFunc :: ModularWeights
           -> Mass          -- ^ Higgs mass
           -> (Mass, Mass)  -- ^ (mtMS, mbMS)
           -> Double        -- ^ alpha_s
           -> Double        -- ^ tan(beta)
           -> (Double -> Double)
mHiggsFunc cs (Mass mh) (mtMS, mbMS) as tanb m0 =
    mHiggs cs (mtMS, mbMS) as tanb m0 - mh

getM0Sol :: ModularWeights
         -> Mass              -- ^ Higgs mass
         -> (Mass, Mass)      -- ^ (mtMS, mbMS)
         -> Double            -- ^ alpha_s
         -> (Double, Double)  -- ^ (low, upper)
         -> Double            -- ^ tan(beta)
         -> Maybe Double
getM0Sol cs mh (mtMS, mbMS) as (xlow, xupper) tanb = do
    let mhFunc = mHiggsFunc cs mh (mtMS, mbMS) as tanb
    if xupper <= xlow
        then Nothing
        else do let param = RiddersParam 1000 (AbsTol 1.0e-3)
                case ridders param (xlow, xupper) mhFunc of
                    Root m0      -> return m0
                    NotBracketed -> getM0Sol cs mh (mtMS, mbMS) as
                                    (xlow, xupper * 0.95) tanb
                    _            -> Nothing

getMHdSq :: ModularWeights
         -> Double  -- ^ kHd
         -> Double  -- ^ M0
         -> Double
getMHdSq ModularWeights { _cHd = cHd } kHd m0 =
    m0 * m0 * (cHd + kHd / (8 * pi2))

-- | (mHu^2, mHd^2, mu)
type HiggsMassParams = (Double, Double, Double)

-- | the mHd solution from the EWSB.
getMHParams :: ModularWeights
            -> Double  -- ^ kHd
            -> Double  -- ^ tan(beta)
            -> Double  -- ^ M0
            -> HiggsMassParams
getMHParams cs kHd tanb m0 =
    ( mHdSq / tanbSq - (0.5 * mZ2 + mu * mu) * (1.0 - 1 / tanbSq)
    , mHdSq
    , mu )
  where
    mHdSq = getMHdSq cs kHd m0
    mu = getMu cs m0
    tanbSq = tanb * tanb

-- | the B solution from the EWSB.
getBmu :: HiggsMassParams -> Double -> Double
getBmu (mHuSq, mHdSq, mu) tanb =
    ((mHuSq + mHdSq) / absmu + 2 * absmu) * tanb / (1 + tanb * tanb)
  where
    absmu = abs mu
