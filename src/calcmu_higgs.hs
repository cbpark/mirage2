{-# LANGUAGE OverloadedStrings #-}

module Main where

import           HEP.Data.AlphaS     (alphasQ, initAlphaS)
import           HEP.Data.Constants  (mhSM, mt)
import           HEP.Data.Interface  (InputArgs (..))
import           HEP.Data.Kinematics (Mass (..))
import           HEP.Data.Quark      (mMSbarHeavy)
import           HEP.Data.SUSY       (ModularWeights (..), getMuFromHiggs)

import qualified Data.Vector.Unboxed as U
import           Options.Generic     (unwrapRecord)

import           Data.Maybe          (fromMaybe)
import           System.IO           (IOMode (..), hPutStrLn, withFile)
import           Text.Printf         (hPrintf)

main :: IO ()
main = do
    input <- unwrapRecord "Calculate mu from m_h (125 GeV)"

    alphaS <- initAlphaS
    (mtMS,    _, _) <- mMSbarHeavy alphaS mt
    (   _, mbMS, _) <- mMSbarHeavy alphaS mtMS
    a3 <- alphasQ mtMS alphaS

    let mStar   = msusy input
        outfile = output input

        tanbs = U.enumFromStepN 6.0 0.05 800
        muF mh = fromMaybe 0
                 . getMuFromHiggs point1 mStar mh (mtMS, mbMS) a3 (1e+2, 1e+5)

        deltaMH = Mass 1.0
        mus0 = U.map (muF mhSM) tanbs
        mus1 = U.map (muF (mhSM - deltaMH)) tanbs
        mus2 = U.map (muF (mhSM + deltaMH)) tanbs
        result = U.zip4 tanbs mus0 mus1 mus2

    withFile outfile WriteMode $ \h -> do
        hPutStrLn h "# tan(beta), mu, mu(lower), mu(upper)"
        U.mapM_ (\(tanb, m0, m1, m2) ->
                     hPrintf h "%6.2f  %11.4f  %11.4f  %11.4f\n"
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
