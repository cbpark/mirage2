module Main where

import HEP.Data.AlphaS
import HEP.Data.Constants
import HEP.Data.SUSY.Parameters

main :: IO ()
main = do
    as <- initAlphaS >>= alphasQ mhSM

    let mh = mHiggs as 10 1000.0
    print mh
