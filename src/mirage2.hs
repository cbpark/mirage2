module Main where

import HEP.Data.AlphaS
import HEP.Data.Constants

main :: IO ()
main = do
    as <- initAlphaS >>= alphasQ mh
    print as
