module Main where

import HEP.Data.AlphaS
import HEP.Data.Constants
import HEP.Data.Kinematics
import HEP.Data.Quark
import HEP.Data.SUSY

main :: IO ()
main = do
    alphaS <- initAlphaS
    a3 <- alphasQ mhSM alphaS
    print a3

    let mu = getMu point1 1000.0
    putStrLn $ "mu = " ++ show mu

    let mtS = mStop point1 mu
    putStrLn $ "mstop = " ++ show mtS

    let mSUSY = getMSUSY point1 mu
    putStrLn $ "mSUSY = " ++ show mSUSY

    let mbS = mSbottom point1 mu
    putStrLn $ "msbottom = " ++ show mbS

    (mtMS,    _, _) <- mMSbarHeavy alphaS (getMass mt)
    (   _, mbMS, _) <- mMSbarHeavy alphaS (getMass mtMS)
    print (mtMS, mbMS)
    let mh = mHiggs point1 mu (mtMS, mbMS) a3
    print mh

point1 :: ModelParams
point1 = ModelParams { _M0 = 2950.0
                     , _c  = ModularWeights { _Hu = 0.0
                                            , _Hd = 1.0
                                            , _Q  = 0.5
                                            , _tR = 0.5
                                            , _bR = 0.5
                                            , _L  = 0.5
                                            , _eR = 0.5 }
                     , _tanb = 10.0
                     }
