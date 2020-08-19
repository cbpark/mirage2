module Main where

import           HEP.Data.AlphaS     (alphasQ, initAlphaS)
import           HEP.Data.Constants  (mhSM, mt)
import           HEP.Data.Kinematics (Mass (..))
import           HEP.Data.Quark      (mMSbarHeavy)
import           HEP.Data.SUSY       (ModularWeights (..), getM0FromHiggs)

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

    let tanbs = U.enumFromStepN 6.0 0.2 200
        getM0 mh = fromMaybe 0
                   . getM0FromHiggs point1 mh (mtMS, mbMS) a3 (1e+3, 1e+5)

        deltaMH = Mass 1.0
        m0s0 = U.map (getM0 mhSM) tanbs
        m0s1 = U.map (getM0 (mhSM - deltaMH)) tanbs
        m0s2 = U.map (getM0 (mhSM + deltaMH)) tanbs
        result = U.zip4 tanbs m0s0 m0s1 m0s2

    withFile outfile WriteMode $ \h -> do
        hPutStrLn h "# tan(beta), M0, M0(lower), M0(upper)"
        U.mapM_ (\(tanb, m0, m1, m2) ->
                     hPrintf h "%5.1f  %11.4f  %11.4f  %11.4f\n"
                     tanb m0 m1 m2) result

    putStrLn $ "-- " <> outfile <> " generated."

point1 :: ModularWeights
point1 = ModularWeights { _cHu = 0.0
                        , _cHd = 1.0
                        , _cQ  = 0.5
                        , _ctR = 0.5
                        , _cbR = 0.5
                        , _cL  = 0.5
                        , _ceR = 0.5 }
