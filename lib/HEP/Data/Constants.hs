module HEP.Data.Constants where

import HEP.Data.Kinematics             (Mass (..), massSq)

import Numeric.MathFunctions.Constants (m_1_sqrt_2)

mW, mZ, mhSM :: Mass
mW   = Mass 80.379
mZ   = Mass 91.1876
mhSM = Mass 125.09

mW2, mZ2 :: Double
mW2 = massSq mW
mZ2 = massSq mZ

-- | pole masses
mt, mb, mc, mtau :: Mass
mt   = Mass 173.0
mb   = Mass 4.78
mc   = Mass 1.67
mtau = Mass 1.777

mt2, mb2, mc2, mtau2 :: Double
mt2   = massSq mt
mb2   = massSq mb
mc2   = massSq mc
mtau2 = massSq mtau

gFermi, gW, gW2 :: Double
gFermi = 1.1663787e-5
gW2    = 8 * mW2 * gFermi * m_1_sqrt_2
gW     = sqrt gW2

-- | vEW = 174.1 GeV.
vEW, vEW2 :: Double
vEW2 = 2 * mW2 / gW2
vEW  = sqrt vEW2

alphasMZ, sinThetaW2 :: Double
alphasMZ = 0.118
sinThetaW2 = 0.23122

b2, b3 :: Double
b2 = 1.0
b3 = - 3.0

pi2 :: Double
pi2 = 9.869604401089358

-- | 1 / (16 * pi).
loopFac :: Double
loopFac = 0.0625 / pi2
