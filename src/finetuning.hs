module Main where

import           HEP.Data.SUSY

import qualified Data.Vector.Unboxed as U

import           Data.Maybe          (fromMaybe)
import           System.Environment  (getArgs)
import           System.IO           (IOMode (..), hPutStrLn, withFile)
import           Text.Printf         (hPrintf)

main :: IO ()
main = do
    outfile <- head <$> getArgs

    -- let tanb = 10
    --     m0sol = getM0FromDBI point1 0 0.01 (1e+3, 1e+5) tanb
    -- putStrLn $ "m0sol = " ++ show m0sol

    let kHd = 0
        tanbs = U.enumFromStepN 6.0 0.2 200
        getM0 dBI = fromMaybe 0
                    . getM0FromDBI point1 kHd dBI (1e+3, 1e+5)

        m0s0 = U.map (getM0 0.005) tanbs
        m0s1 = U.map (getM0 0.01 ) tanbs
        m0s2 = U.map (getM0 0.05 ) tanbs
        m0s3 = U.map (getM0 0.1  ) tanbs
        result = U.zip5 tanbs m0s0 m0s1 m0s2 m0s3

    withFile outfile WriteMode $ \h -> do
        hPutStrLn h "# tan(beta), M0(1), M0(2), M0(3), M0(4)"
        U.mapM_ (\(tanb, m0, m1, m2, m3) ->
                     hPrintf h "%5.1f  %11.4f  %11.4f  %11.4f  %11.4f\n"
                     tanb m0 m1 m2 m3) result

    putStrLn $ "-- " <> outfile <> " generated."

point1 :: ModularWeights
point1 = ModularWeights { _cHu = 0.0
                        , _cHd = 1.0
                        , _cQ  = 0.5
                        , _ctR = 0.5
                        , _cbR = 0.5
                        , _cL  = 0.5
                        , _ceR = 0.5 }
