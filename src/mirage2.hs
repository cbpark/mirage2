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

    let tanb = 10.0
        m0 = 3000.0

    let mu = getMu point1 1000.0 m0
    putStrLn $ "mu = " ++ show mu

    let mtS = mStop point1 tanb m0 mu
        mbS = mSbottom point1 tanb m0 mu
    putStrLn $ "mstop = " ++ show mtS
    putStrLn $ "msbottom = " ++ show mbS

    let mSUSY = getMSUSY point1 tanb m0 mu
    putStrLn $ "mSUSY = " ++ show mSUSY

    (mtMS,    _, _) <- mMSbarHeavy alphaS (getMass mt)
    (   _, mbMS, _) <- mMSbarHeavy alphaS (getMass mtMS)
    print (mtMS, mbMS)

    let mh = mHiggs point1 (mtMS, mbMS) a3 mu tanb m0
    print mh

point1 :: ModularWeights
point1 = ModularWeights { _cHu = 0.0
                        , _cHd = 1.0
                        , _cQ  = 0.5
                        , _ctR = 0.5
                        , _cbR = 0.5
                        , _cL  = 0.5
                        , _ceR = 0.5 }
