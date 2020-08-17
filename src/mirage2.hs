module Main where

import HEP.Data.AlphaS
import HEP.Data.Constants
import HEP.Data.SUSY.Parameters
import HEP.Data.SUSY.Squark

main :: IO ()
main = do
    as <- initAlphaS >>= alphasQ mhSM
    print as

    let mu = getMu point1 1000.0
    print mu

    let mtS = mStop point1 mu
    print mtS

    let mSUSY = getMSUSY point1 mu
    print mSUSY

    -- let mbS = mSbottom point1 mu
    -- print mbS

    -- let mh = mHiggs point1 mu as
    -- print mh

point1 :: ModelParams
point1 = ModelParams { _M0 = 1500.0
                     , _c  = ModularWeights { _Hu = 0.0
                                            , _Hd = 1.0
                                            , _Q  = 0.5
                                            , _tR = 0.5
                                            , _bR = 0.5
                                            , _L  = 0.5
                                            , _eR = 0.5 }
                     , _tanb = 10.0
                     }
