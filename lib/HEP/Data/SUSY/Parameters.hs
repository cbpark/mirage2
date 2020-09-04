module HEP.Data.SUSY.Parameters where

-- import HEP.Data.Constants  (b2, b3, gW2, pi2)
import HEP.Data.Kinematics (Mass (..))

-- data ModelParams = ModelParams { _M0   :: Double
--                                , _c    :: ModularWeights
--                                , _tanb :: Double
--                                } deriving Show

data ModularWeights = ModularWeights { _cHu :: !Double
                                     , _cHd :: !Double
                                     , _cQ  :: !Double
                                     , _ctR :: !Double
                                     , _cbR :: !Double
                                     , _cL  :: !Double
                                     , _ceR :: !Double
                                     } deriving Show

-- data HuHdLoopFactor = HuHdLoopFactor { _kHu :: Double
--                                      , _kHd :: Double
--                                      } deriving Show

data HiggsParams = HiggsParams { _M0      :: !Double
                               , _tanbeta :: !Double
                               , _mHu2    :: !Double
                               , _mHd2    :: !Double
                               , _B       :: !Double
                               , _mu      :: !Double
                               } deriving Show

foreign import ccall "math.h cbrt" cbrt :: Double -> Double

getMu :: ModularWeights
      -> Double          -- ^ m_*
      -> (Mass, Double)  -- ^ (mtMS, alpha_s(mtMS))
      -> Double          -- ^ M_0
      -> Double
getMu ModularWeights {_cHd = cHd} mStar _ m0 =
    -- mStar * (mStar / m0) ** (7.0 / 12) / cHd ** (1.0 / 8)
    mStar / cHd ** (1.0 / 8)
    * (sqrt . sqrt) (mStar / m0)  -- (mStar / m0) ** (1.0 / 4)
    * cbrt (mStar / m2)           -- (mStar / m2) ** (1.0 / 3)
    * cbrt (m3    / m2) ** 7      -- (m3    / m2) ** (7.0 / 3)
  where
    -- cf = 1 / (8 * pi2) * log (m0 / mtMS)
    m2 = m0 -- m0 * (1 - b2 * gW2 * cf)
    m3 = m0 -- m0 * (1 - b3 * (4 * pi * as) * cf)

cos2Beta :: Double -> Double
cos2Beta tanb = (1.0 - tanbSq) / (1.0 + tanbSq)
  where tanbSq = tanb * tanb

cosBeta :: Double -> Double
cosBeta tanb = 1.0 / sqrt (1 + tanb * tanb)

sinBeta :: Double -> Double
sinBeta tanb = tanb / sqrt (1 + tanb * tanb)
