{-# LANGUAGE Strict #-}

module HEP.Data.Util where

import Numeric.RootFinding

riddersSolver :: (Double -> Double) -> (Double, Double) -> Maybe Double
riddersSolver f (xlow, xup) =
    if xup <= xlow
    then Nothing
    else do let param = RiddersParam 1000 (AbsTol 1e-3)
            case ridders param (xlow, xup) f of
                Root x       -> return x
                NotBracketed -> riddersSolver f (xlow, xup * 0.95)
                _            -> Nothing
