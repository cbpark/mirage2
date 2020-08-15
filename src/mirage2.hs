module Main where

import HEP.Data.AlphaS

main :: IO ()
main = do
    as <- initAlphaS
    alphasQ as 125.0 >>= print
