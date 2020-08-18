module Main where

import           HEP.Data.AlphaS     (alphasQ, initAlphaS)
import           HEP.Data.Constants  (mhSM, mt)
import           HEP.Data.Quark      (mMSbarHeavy)
import           HEP.Data.SUSY       (ModularWeights (..), getM0Sol)

import qualified Data.Vector.Unboxed as U

import           Data.Maybe          (fromMaybe)
import           System.Environment  (getArgs)
import           System.IO           (IOMode (..), hPutStrLn, withFile)
import           Text.Printf         (hPrintf)

main :: IO ()
main = do
    outfile <- head <$> getArgs

    alphaS <- initAlphaS
    (mtMS,    _, _) <- mMSbarHeavy alphaS mt
    (   _, mbMS, _) <- mMSbarHeavy alphaS mtMS
    a3 <- alphasQ mtMS alphaS

    let tanbs = U.enumFromStepN 6 0.2 10 -- 102
        m0s = U.map (fromMaybe 0
                     . getM0Sol point1 mhSM (mtMS, mbMS) a3 (1e+3, 1e+4)) tanbs
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
