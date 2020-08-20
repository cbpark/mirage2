module Main where

import HEP.Data.SUSY

main :: IO ()
main = do
    let tanb = 10
        m0sol = getM0FromDBI point1 0 0.01 (1e+3, 1e+5) tanb
    putStrLn $ "m0sol = " ++ show m0sol

point1 :: ModularWeights
point1 = ModularWeights { _cHu = 0.0
                        , _cHd = 1.0
                        , _cQ  = 0.5
                        , _ctR = 0.5
                        , _cbR = 0.5
                        , _cL  = 0.5
                        , _ceR = 0.5 }
