module HEP.Data.SUSY.Higgs where

import HEP.Data.Constants       (loopFac, mZ2, mtau, pi2, vEW, vEW2)
import HEP.Data.Kinematics      (Mass (..), massSq)
import HEP.Data.SUSY.Parameters
import HEP.Data.SUSY.Squark     (getMSUSY2, mSbottom, mStop)
import HEP.Data.Util            (riddersSolver)

-- | From <https://arxiv.org/abs/1112.3336 arXiv:1112.3336> and
--   <https://arxiv.org/abs/hep-ph/0001002 arXiv:hep-ph/0001002>
--
--   valid for moderate or large tan(beta) and
--   (m_{~t_2}^2 - m_{~t_1}^2) / (m_{~t_2}^2 + m_{~t_1}^2) < 0.5.
mHiggs :: ModularWeights
       -> Double        -- ^ m_*
       -> (Mass, Mass)  -- ^ (mtMS, mbMS)
       -> Double        -- ^ alpha_s
       -> Double        -- ^ tan(beta)
       -> Double        -- ^ M_0
       -> Double
mHiggs cs@ModularWeights {_cL = cL} mStar (mtMS, mbMS) as tanb m0
    | mhSq <= 0 = 0
    | otherwise = sqrt mhSq
  where
    mhSq = mHiggsTreeSq cs tanb m0 -- mZ2 * cos2b * cos2b
           + 3.0 / (4 * pi2) * mt2 * mt2 / vEW2 * termT
           - yb2 * yb2 * vEW2 * loopFac * termB
           - ytau2 * ytau2 * vEW2 * loopFac / 3 * termTau

    mt2 = massSq mtMS
    mu = getMu cs mStar m0

    -- mSUSY = getMSUSY cs tanb m0 mu
    -- mSUSY2 = mSUSY * mSUSY
    mSUSY2 = getMSUSY2 cs tanb m0 mu
    loopT = log (1 + mSUSY2 / mt2)

    aT = m0 - mu / tanb
    aT2 = aT * aT / mSUSY2
    xT = 2 * aT2 * (1 - aT2 / 12)
    yt = getMass mtMS / vEW

    termT = 0.5 * xT + loopT
            + loopFac * (1.5 * mt2 / vEW2 - 32 * pi * as)
            * (xT * loopT + loopT * loopT)
            + (4 * as / (3 * pi) - 5 * yt * yt * loopFac) * loopT

    (mst1, mst2) = mStop    cs tanb m0 mu
    (msb1, msb2) = mSbottom cs tanb m0 mu
    mu2 = mu * mu
    mu4 = mu2 * mu2
    cosb = cosBeta tanb
    -- from Eq.(40) of https://arxiv.org/pdf/hep-ph/9402253.pdf
    mgluino = m0
    dhb = (2 * as / 3 * mgluino
           * integralFunc (massSq msb1) (massSq msb2) (mgluino * mgluino)
           + yt / 4 * aT
           * integralFunc (massSq mst1) (massSq mst2) mu2)
          * mu * tanb / pi
    yb = getMass mbMS / (vEW * cosb * (1 + dhb))
    yb2 = yb * yb
    termB = mu4 / (mSUSY2 * mSUSY2)
            * (1 + loopFac * loopT * (9 * yb2 - 5 * mt2 / vEW2 - 64 * pi * as))

    -- the contribution from stau
    ytau = getMass mtau / (vEW * cosb)
    ytau2 = ytau * ytau
    mStau2 = m0 * m0 * cL
    termTau = mu4 / (mStau2 * mStau2)

-- | Squared tree-level Higgs mass.
--
--   for mA >> mZ, it is mZ2 * cos2b * cos2b.
mHiggsTreeSq :: ModularWeights -> Double -> Double -> Double
mHiggsTreeSq ModularWeights {_cHd = cHd} tanb m0 =
    0.5 * (mA2mZ2
           - sqrt (mA2mZ2 * mA2mZ2 - 4 * mA2 * mZ2 * cos2b * cos2b))
  where
    mA2 = cHd * m0 * m0
    mA2mZ2 = mA2 + mZ2
    cos2b = cos2Beta tanb

integralFunc :: Double -> Double -> Double -> Double
integralFunc a b c = num / den
  where
    num = a * b * log (a / b) + b * c * log (b / c) + a * c * log (c / a)
    den = (a - b) * (b - c) * (a - c)

getM0FromHiggs :: ModularWeights
               -> Double            -- ^ m_*
               -> Mass              -- ^ Higgs mass
               -> (Mass, Mass)      -- ^ (mtMS, mbMS)
               -> Double            -- ^ alpha_s
               -> (Double, Double)  -- ^ (low, upper)
               -> Double            -- ^ tan(beta)
               -> Maybe Double
getM0FromHiggs cs mStar (Mass mh) mqMS as range tanb = riddersSolver mhFunc range
  where mhFunc m0 = mHiggs cs mStar mqMS as tanb m0 - mh

getMHParams :: ModularWeights
            -> Double            -- ^ m_*
            -> Double            -- ^ kHd
            -> Mass              -- ^ SM Higgs mass
            -> (Mass, Mass)      -- ^ (mt, mb)
            -> Double            -- ^ alpha_s
            -> (Double, Double)  -- ^ the range of M0
            -> Double            -- ^ tan(beta)
            -> Maybe HiggsParams
getMHParams cs mStar kHd mh mqs a3 xrange tanb = do
    m0 <- getM0FromHiggs cs mStar mh mqs a3 xrange tanb
    let mH2@(mHu2, mHd2, mu) = getMHParams' cs mStar kHd tanb m0
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
             -> Double  -- ^ m_*
             -> Double  -- ^ kHd
             -> Double  -- ^ tan(beta)
             -> Double  -- ^ M0
             -> HiggsMassParams
getMHParams' cs mStar kHd tanb m0
    | tanb <= 0 = (0, 0, 0)
    | otherwise = ( mHdSq / tanbSq - (0.5 * mZ2 + mu * mu) * (1 - 1 / tanbSq)
                  , mHdSq
                  , mu )
  where
    mHdSq = getMHdSq cs kHd m0
    mu = getMu cs mStar m0
    tanbSq = tanb * tanb

getMHdSq :: ModularWeights
         -> Double  -- ^ kHd
         -> Double  -- ^ M0
         -> Double
getMHdSq ModularWeights { _cHd = cHd } kHd m0 =
    m0 * m0 * (cHd + kHd / (8 * pi2))

-- | the B solution from the EWSB.
getB :: HiggsMassParams -> Double -> Double
getB (mHuSq, mHdSq, mu) tanb
    | mu == 0   = 0
    | otherwise = ((mHuSq + mHdSq) / absmu + 2 * absmu)
                  * tanb / (1 + tanb * tanb)
  where
    absmu = abs mu

getM0FromEWSB :: ModularWeights
              -> Double            -- ^ m_*
              -> Double            -- ^ kHd
              -> (Double, Double)  -- ^ (xlow, xup)
              -> Double            -- ^ tan(beta)
              -> Maybe Double
getM0FromEWSB cs mStar kHd range tanb = riddersSolver ewsbF range
  where
    ewsbF m0 = let mH@(mHu2, mHd2, mu) = getMHParams' cs mStar kHd tanb m0
                   b = getB mH tanb
                   mu2 = mu * mu
               in b * b * mu2 - (mHu2 + mu2) * (mHd2 + mu2)

getM0FromFT :: (HiggsParams -> Double)  -- ^ the function to calculate Delta(FT)
            -> ModularWeights           -- ^ modular weights
            -> Double                   -- ^ m_*
            -> Double                   -- ^ kHd
            -> Double                   -- ^ the Delta^{-1} value
            -> (Double, Double)         -- ^ (xlow, xup)
            -> Double                   -- ^ tan(beta)
            -> Maybe Double
getM0FromFT deltaF cs mStar kHd delta range tanb = riddersSolver f range
  where
    f m0 = let mH@(mHu2, mHd2, mu) = getMHParams' cs mStar kHd tanb m0
               b = getB mH tanb
               hparams = HiggsParams { _M0      = m0
                                     , _tanbeta = tanb
                                     , _mHu2    = mHu2
                                     , _mHd2    = mHd2
                                     , _B       = b
                                     , _mu      = mu }
           in (1 / deltaF hparams - delta)

getM0FromFTB :: ModularWeights
             -> Double            -- ^ m_*
             -> Double            -- ^ kHd
             -> Double            -- ^ the (Delta_B)^{-1} value
             -> (Double, Double)  -- ^ (xlow, xup)
             -> Double            -- ^ tan(beta)
             -> Maybe Double
getM0FromFTB = getM0FromFT deltaB

deltaB :: HiggsParams -> Double
deltaB HiggsParams {_tanbeta = tanb, _B = b, _mu = mu} =
    4 * t2 / (t2 - 1) ** 2 * (1 + (t2 + 1) / tanb * bMu / mZ2)
  where
    t2 = tanb * tanb
    bMu = abs (b * mu)

-- | obtain M0 from the condition that B = k * M0.
getM0FromB :: ModularWeights
           -> Double            -- ^ m_*
           -> Double            -- ^ kHd
           -> Double            -- ^ k of B = k * M0
           -> (Double, Double)  -- ^ (xlow, xup)
           -> Double            -- ^ tan(beta)
           -> Maybe Double
getM0FromB cs mStar kHd k range tanb = riddersSolver bF range
  where
    bF m0 = let b = getB (getMHParams' cs mStar kHd tanb m0) tanb
            in b - k * m0

getMuFromHiggs :: ModularWeights
               -> Double            -- ^ m_*
               -> Mass              -- ^ Higgs mass
               -> (Mass, Mass)      -- ^ (mtMS, mbMS)
               -> Double            -- ^ alpha_s
               -> (Double, Double)  -- ^ (low, upper)
               -> Double            -- ^ tan(beta)
               -> Maybe Double
getMuFromHiggs cs mStar mh mqMS as range tanb =
    getMu cs mStar <$> getM0FromHiggs cs mStar mh mqMS as range tanb

getM0FromFTMuSq :: ModularWeights
                -> Double            -- ^ m_*
                -> Double            -- ^ kHd
                -> Double            -- ^ the (Delta_Mu^2)^{-1} value
                -> (Double, Double)  -- ^ (xlow, xup)
                -> Double            -- ^ tan(beta)
                -> Maybe Double
getM0FromFTMuSq = getM0FromFT (abs . deltaMuSq)

deltaMuSq :: HiggsParams -> Double
deltaMuSq HiggsParams {_tanbeta = tanb, _B = b', _mu = mu'} =
    - 2 * mu * mu / mZ2
    + 2 * t2 / (t2 - 1) ** 2 * (1 - 4 * tanb / (t2 + 1) * mu / b)
    * (1 + (t2 + 1) / tanb * b * mu / mZ2)
  where
    t2 = tanb * tanb
    b = abs b'
    mu = abs mu'
