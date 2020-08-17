module HEP.Data.Constants where

import HEP.Data.Kinematics (Mass (..), massSq)

mW, mZ, mh, mt :: Mass
mW = Mass 80.379
mZ = Mass 91.1876
mh = Mass 125.10
mt = Mass 173.0

mW2 :: Double
mW2 = massSq mW

gFermi, gW, gW2 :: Double
gFermi = 1.1663787e-5
gW2    = 8 * mW2 * gFermi / sqrt2
gW     = sqrt gW2

vEW, vEW2 :: Double
vEW2 = 2 * mW2 / gW2
vEW  = sqrt vEW2

alphasMZ :: Double
alphasMZ = 0.118

sqrt2 :: Double
sqrt2 = 1.4142135623730951
