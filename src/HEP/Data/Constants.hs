module HEP.Data.Constants where

import HEP.Data.Kinematics (Mass (..), massSq)

mW, mZ, mhSM, mt, mb :: Mass
mW   = Mass 80.379
mZ   = Mass 91.1876
mhSM = Mass 125.10
mt   = Mass 173.0
mb   = Mass 4.78

mW2, mZ2, mt2 :: Double
mW2 = massSq mW
mZ2 = massSq mZ
mt2 = massSq mt

gFermi, gW, gW2 :: Double
gFermi = 1.1663787e-5
gW2    = 8 * mW2 * gFermi / sqrt2
gW     = sqrt gW2

-- | vEW = 174.1 GeV.
vEW, vEW2 :: Double
vEW2 = 2 * mW2 / gW2
vEW  = sqrt vEW2

alphasMZ, sinThetaW2 :: Double
alphasMZ = 0.118
sinThetaW2 = 0.23122

sqrt2 :: Double
sqrt2 = 1.4142135623730951

pi2 :: Double
pi2 = 9.869604401089358
