module Main where

import HEP.Data.SUSY

main :: IO ()
main = do
    let mH2@(mHu2, mHd2, mu) = getMHParams point1 0 10 3000.0
    putStrLn $ "mHu2 = " <> show mHu2 <> ", mHd2 = " <> show mHd2
        <> ", mu = " <> show mu
    let bH = getBmu mH2 10
    putStrLn $ "B = " <> show bH

point1 :: ModularWeights
point1 = ModularWeights { _cHu = 0.0
                        , _cHd = 1.0
                        , _cQ  = 0.5
                        , _ctR = 0.5
                        , _cbR = 0.5
                        , _cL  = 0.5
                        , _ceR = 0.5 }
