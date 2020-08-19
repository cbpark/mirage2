module Main where

import HEP.Data.AlphaS    (alphasQ, initAlphaS)
import HEP.Data.Constants (mhSM, mt)
import HEP.Data.Quark     (mMSbarHeavy)
import HEP.Data.SUSY

import Data.Maybe

main :: IO ()
main = do
    -- let mH2@(mHu2, mHd2, mu) = getMHParams' point1 0 10 3000.0
    -- putStrLn $ "mHu2 = " <> show mHu2 <> ", mHd2 = " <> show mHd2
    --     <> ", mu = " <> show mu
    -- let bH = getBmu mH2 10
    -- putStrLn $ "B = " <> show bH

    alphaS <- initAlphaS
    (mtMS,    _, _) <- mMSbarHeavy alphaS mt
    (   _, mbMS, _) <- mMSbarHeavy alphaS mtMS
    a3 <- alphasQ mtMS alphaS

    let tanb = 11
        mHparams = getMHParams point1 0 mhSM (mtMS, mbMS) a3 (1e+3, 1e+4) tanb
        dB = deltaB (fromJust mHparams)
    print mHparams
    putStrLn $ "delta B = " ++ show (1.0 / dB)

point1 :: ModularWeights
point1 = ModularWeights { _cHu = 0.0
                        , _cHd = 1.0
                        , _cQ  = 0.5
                        , _ctR = 0.5
                        , _cbR = 0.5
                        , _cL  = 0.5
                        , _ceR = 0.5 }
