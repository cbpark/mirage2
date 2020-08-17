module Main where

import HEP.Data.AlphaS    (alphasQ, initAlphaS)
import HEP.Data.Constants (mt)
import HEP.Data.Quark     (mMSbarHeavy)
import HEP.Data.SUSY

main :: IO ()
main = do
    alphaS <- initAlphaS

    let tanb = 10.0
        m0 = 2230.0

    let mu = getMu point1 1000.0 m0
    putStrLn $ "mu = " ++ show mu

    let mtS = mStop point1 tanb m0 mu
        mbS = mSbottom point1 tanb m0 mu
    putStrLn $ "mstop = " ++ show mtS
    putStrLn $ "msbottom = " ++ show mbS

    let mSUSY = getMSUSY point1 tanb m0 mu
    putStrLn $ "mSUSY = " ++ show mSUSY

    (mtMS,    _, _) <- mMSbarHeavy alphaS mt
    (   _, mbMS, _) <- mMSbarHeavy alphaS mtMS
    print (mtMS, mbMS)

    a3 <- alphasQ mtMS alphaS
    print a3

    let mh = mHiggs point1 (mtMS, mbMS) a3 mu tanb m0
    print mh

    let m0sol = getM0Sol point1 (mtMS, mbMS) a3 mu (5.0e+2, 5.0e+4) 30
    putStrLn $ "m0sol = " ++ show m0sol

point1 :: ModularWeights
point1 = ModularWeights { _cHu = 0.0
                        , _cHd = 1.0
                        , _cQ  = 0.5
                        , _ctR = 0.5
                        , _cbR = 0.5
                        , _cL  = 0.5
                        , _ceR = 0.5 }
