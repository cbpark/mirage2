module Main where

import           HEP.Data.SUSY       (ModularWeights (..), getM0FromStop)

import qualified Data.Vector.Unboxed as U

import           Data.Maybe          (fromMaybe)
import           System.Environment  (getArgs)
import           System.IO           (IOMode (..), hPutStrLn, withFile)
import           Text.Printf         (hPrintf)

main :: IO ()
main = do
    outfile <- head <$> getArgs

    let mstop = 1000.0
        tanbs = U.enumFromStepN 6.0 0.2 200
        getM0 = fromMaybe 0 . getM0FromStop point1 mstop (1e+2, 1e+4)

        m0s = U.map getM0 tanbs
        result = U.zip tanbs m0s

    withFile outfile WriteMode $ \h -> do
        hPutStrLn h "# tan(beta), M0"
        U.mapM_ (uncurry (hPrintf h "%5.1f  %11.4f\n")) result

    putStrLn $ "-- " <> outfile <> " generated."

point1 :: ModularWeights
point1 = ModularWeights { _cHu = 0.0
                        , _cHd = 1.0
                        , _cQ  = 0.5
                        , _ctR = 0.5
                        , _cbR = 0.5
                        , _cL  = 0.5
                        , _ceR = 0.5 }
