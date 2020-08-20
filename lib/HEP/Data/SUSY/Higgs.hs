{-# LANGUAGE RecordWildCards #-}

module HEP.Data.SUSY.Higgs where

import HEP.Data.Constants       (mZ2, mtau, pi2, vEW, vEW2)
import HEP.Data.Kinematics      (Mass (..), massSq)
import HEP.Data.SUSY.Parameters
import HEP.Data.SUSY.Squark     (getMSUSY2, mSbottom, mStop)
import HEP.Util                 (riddersSolver)

-- | From <https://arxiv.org/abs/1112.3336 arXiv:1112.3336> and
--   <https://arxiv.org/abs/hep-ph/0001002 arXiv:hep-ph/0001002>
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
    -- for mA >> mZ, it is mZ2 * cos2b * cos2b.
    mh0Sq = 0.5 * (mA2 + mZ2
                   - sqrt ((mA2 + mZ2) ** 2 - 4 * mA2 * mZ2 * cos2b * cos2b))
    mA2 = _cHd * m0 * m0

    mhSq = mh0Sq -- mZ2 * cos2b * cos2b
           + 3.0 / (4 * pi2) * mt2 * mt2 / vEW2 * termT
           - yb2 * yb2 * vEW2 * loopFac * termB
           - ytau2 * ytau2 * vEW2 * loopFac / 3 * termTau

    [cos2b, cosb] = fmap ($ tanb) [cos2Beta, cosBeta]
    mt2 = massSq mtMS
    mu = getMu cs m0

    -- mSUSY = getMSUSY cs tanb m0 mu
    -- mSUSY2 = mSUSY * mSUSY
    mSUSY2 = getMSUSY2 cs tanb m0 mu
    loopT = log (1 + mSUSY2 / mt2)

    aT = m0 - mu / tanb
    aT2 = aT * aT
    xT = 2.0 * aT2 / mSUSY2 * (1.0 - aT2 / (12.0 * mSUSY2))
    loopFac = 1.0 / (16 * pi2)
    yt = getMass mtMS / vEW

    termT = 0.5 * xT + loopT
            + loopFac * (1.5 * mt2 / vEW2 - 32 * pi * as)
            * (xT * loopT + loopT * loopT)
            + (4 * as / (3 * pi) - 5 * yt * yt * loopFac) * loopT

    -- the contribution from sbottom
    (mst1, mst2) = mStop    cs tanb m0 mu
    (msb1, msb2) = mSbottom cs tanb m0 mu
    mu2 = mu * mu
    -- from Eq.(40) of https://arxiv.org/pdf/hep-ph/9402253.pdf
    mgluino = m0
    dhb = (2 * as / 3 * mgluino
           * integralFunc (massSq msb1) (massSq msb2) (mgluino * mgluino)
           + yt / 4 * aT
           * integralFunc (massSq mst1) (massSq mst2) (mu * mu))
          * mu * tanb / pi
    yb = getMass mbMS / (vEW * cosb * dhb)
    yb2 = yb * yb
    mu4 = mu2 * mu2
    termB = mu4 / (mSUSY2 * mSUSY2)
            * (1.0 + loopFac
               * loopT * (9.0 * yb2 - 5 * mt2 / vEW2 - 64 * pi * as))

    -- the contribution from stau
    ytau = getMass mtau / (vEW * cosb)
    ytau2 = ytau * ytau
    mStau = m0 * sqrt _cL
    termTau = mu4 / mStau ** 4

integralFunc :: Double -> Double -> Double -> Double
integralFunc a b c = num / den
  where
    num = a * b * log (a / b) + b * c * log (b / c) + a * c * log (c / a)
    den = (a - b) * (b - c) * (a - c)

mHiggsFunc :: ModularWeights
           -> Mass          -- ^ Higgs mass
           -> (Mass, Mass)  -- ^ (mtMS, mbMS)
           -> Double        -- ^ alpha_s
           -> Double        -- ^ tan(beta)
           -> (Double -> Double)
mHiggsFunc cs (Mass mh) mqMS as tanb m0 = mHiggs cs mqMS as tanb m0 - mh

getM0FromHiggs :: ModularWeights
               -> Mass              -- ^ Higgs mass
               -> (Mass, Mass)      -- ^ (mtMS, mbMS)
               -> Double            -- ^ alpha_s
               -> (Double, Double)  -- ^ (low, upper)
               -> Double            -- ^ tan(beta)
               -> Maybe Double
getM0FromHiggs cs mh mqMS as range tanb = do
    let mhFunc = mHiggsFunc cs mh mqMS as tanb
    riddersSolver mhFunc range

getMHParams :: ModularWeights
            -> Double            -- ^ kHd
            -> Mass              -- ^ SM Higgs mass
            -> (Mass, Mass)      -- ^ (mt, mb)
            -> Double            -- ^ alpha_s
            -> (Double, Double)  -- ^ the range of M0
            -> Double            -- ^ tan(beta)
            -> Maybe HiggsParams
getMHParams cs kHd mh mqs a3 xrange tanb = do
    m0 <- getM0FromHiggs cs mh mqs a3 xrange tanb
    let mH2@(mHu2, mHd2, mu) = getMHParams' cs kHd tanb m0
        bH = getB mH2 tanb
    return $ HiggsParams { _M0      = m0
                         , _tanbeta = tanb
                         , _mHu2    = mHu2
                         , _mHd2    = mHd2
                         , _B       = bH
                         , _mu      = mu }

-- | (mHu^2, mHd^2, mu)
type HiggsMassParams = (Double, Double, Double)

-- | the mHd solution from the EWSB.
getMHParams' :: ModularWeights
             -> Double  -- ^ kHd
             -> Double  -- ^ tan(beta)
             -> Double  -- ^ M0
             -> HiggsMassParams
getMHParams' cs kHd tanb m0 =
    ( mHdSq / tanbSq - (0.5 * mZ2 + mu * mu) * (1.0 - 1 / tanbSq)
    , mHdSq
    , mu )
  where
    mHdSq = getMHdSq cs kHd m0
    mu = getMu cs m0
    tanbSq = tanb * tanb

getMHdSq :: ModularWeights
         -> Double  -- ^ kHd
         -> Double  -- ^ M0
         -> Double
getMHdSq ModularWeights { _cHd = cHd } kHd m0 =
    m0 * m0 * (cHd + kHd / (8 * pi2))

-- | the B solution from the EWSB.
getB :: HiggsMassParams -> Double -> Double
getB (mHuSq, mHdSq, mu) tanb =
    ((mHuSq + mHdSq) / absmu + 2 * absmu) * tanb / (1 + tanb * tanb)
  where
    absmu = abs mu

getM0FromEWSB :: ModularWeights
              -> Double            -- ^ kHd
              -> (Double, Double)  -- ^ (xlow, xup)
              -> Double            -- ^ tan(beta)
              -> Maybe Double
getM0FromEWSB cs kHd range tanb = do
    let ewsbF m0 = let mH@(mHu2, mHd2, mu) = getMHParams' cs kHd tanb m0
                       b = getB mH tanb
                       mu2 = mu * mu
                   in b * b * mu2 - (mHu2 + mu2) * (mHd2 + mu2)
    riddersSolver ewsbF range

getM0FromDBI :: ModularWeights
             -> Double            -- ^ kHd
             -> Double            -- ^ (Delta B)^{-1}
             -> (Double, Double)  -- ^ (xlow, xup)
             -> Double            -- ^ tan(beta)
             -> Maybe Double
getM0FromDBI cs kHd dBI range tanb = do
    let dBIF m0 = let mH@(mHu2, mHd2, mu) = getMHParams' cs kHd tanb m0
                      b = getB mH tanb
                      hparams = HiggsParams { _M0      = m0
                                            , _tanbeta = tanb
                                            , _mHu2    = mHu2
                                            , _mHd2    = mHd2
                                            , _B       = b
                                            , _mu      = mu }
                  in (1.0 / deltaB hparams - dBI)
    riddersSolver dBIF range

deltaB :: HiggsParams -> Double
deltaB HiggsParams {..} =
    4 * t2 / (t2 - 1) ** 2 * (1 + (t2 + 1) / _tanbeta * bMu / mZ2)
  where t2 = _tanbeta * _tanbeta
        bMu = abs (_B * _mu)
