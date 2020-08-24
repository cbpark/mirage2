module Main where

import HEP.Data.AlphaS    (alphasQ, initAlphaS)
import HEP.Data.Constants (mhSM, mt)
import HEP.Data.Quark     (mMSbarHeavy)
import HEP.Data.SUSY

import Data.Maybe         (fromJust)

main :: IO ()
main = do
    alphaS <- initAlphaS

    let tanb = 10.0
        m0 = 4290
        mStar = 1000.0

    let mu = getMu point1 mStar m0
    putStrLn $ "mu = " ++ show mu

    let mtS = mStop point1 tanb m0 mu
        mbS = mSbottom point1 tanb m0 mu
    putStrLn $ "mstop = " ++ show mtS
    putStrLn $ "msbottom = " ++ show mbS

    let mSUSY = sqrt $ getMSUSY2 point1 tanb m0 mu
    putStrLn $ "mSUSY = " ++ show mSUSY

    (mtMS,    _, _) <- mMSbarHeavy alphaS mt
    (   _, mbMS, _) <- mMSbarHeavy alphaS mtMS
    print (mtMS, mbMS)

    a3 <- alphasQ mtMS alphaS
    print a3

    let mh = mHiggs point1 mStar (mtMS, mbMS) a3 tanb m0
    print mh

    let m0sol = getM0FromHiggs point1 mStar
                mhSM (mtMS, mbMS) a3 (mStar, 5e+4) 10
    putStrLn $ "m0sol (Higgs)= " ++ show m0sol

    let m0sol' = getM0FromStop point1 mStar 1000.0 (mStar, 1e+4) 10
    putStrLn $ "m0sol (stop) = " ++ show m0sol'

    let mHparams = getMHParams point1 mStar 0
                   mhSM (mtMS, mbMS) a3 (mStar, 1e+4) tanb
        dB = deltaB (fromJust mHparams)
        dMu2 = deltaMuSq (fromJust mHparams)
    print mHparams
    putStrLn $ "delta B^{-1} = " ++ show (1.0 / dB)
    putStrLn $ "delta Mu2^{-1} = " ++ show (1.0 / dMu2)

point1 :: ModularWeights
point1 = ModularWeights { _cHu = 0.0
                        , _cHd = 1.0
                        , _cQ  = 0.5
                        , _ctR = 0.5
                        , _cbR = 0.5
                        , _cL  = 0.5
                        , _ceR = 0.5 }
