{-# LANGUAGE DataKinds          #-}
{-# LANGUAGE DeriveGeneric      #-}
{-# LANGUAGE FlexibleInstances  #-}
{-# LANGUAGE OverloadedStrings  #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE TypeOperators      #-}

module Main where

import           HEP.Data.AlphaS     (alphasQ, initAlphaS)
import           HEP.Data.Constants  (mt)
import           HEP.Data.Quark      (mMSbarHeavy)
import           HEP.Data.SUSY       (ModularWeights (..), mHiggs)

import qualified Data.Vector.Unboxed as U
import           Options.Generic
import           System.IO           (IOMode (..), hPutStrLn, withFile)

main :: IO ()
main = do
    input <- unwrapRecord "Calculate the Higgs mass for given tan(beta)"
    let mStar   = msusy input
        tanbeta = tanb input
        outfile = output input

    putStrLn $ "-- tan(beta) = " <> show tanbeta

    alphaS <- initAlphaS
    (mtMS,    _, _) <- mMSbarHeavy alphaS mt
    (   _, mbMS, _) <- mMSbarHeavy alphaS mtMS
    a3 <- alphasQ mtMS alphaS

    let m0s = U.enumFromStepN mStar 0.5 18000
        mhs = U.map (mHiggs point1 mStar (mtMS, mbMS) a3 tanbeta) m0s
        result = U.zip m0s mhs

    withFile outfile WriteMode $ \h -> do
        hPutStrLn h "# M0, tan(beta), mh"
        U.mapM_ (\(m0, mh) -> hPutStrLn h (show m0
                                           <> "\t" <> show tanbeta
                                           <> "\t" <> show mh)) result

    putStrLn $ "-- " <> outfile <> " generated."

data InputArgs w = InputArgs
    { msusy  :: w ::: Double <?> "m(sfermion)"
    , tanb   :: w ::: Double <?> "tan(beta)"
    , output :: w ::: String <?> "the name of the output file"
    } deriving Generic

instance ParseRecord (InputArgs Wrapped)
deriving instance Show (InputArgs Unwrapped)

point1 :: ModularWeights
point1 = ModularWeights { _cHu = 0.0
                        , _cHd = 1.0
                        , _cQ  = 0.5
                        , _ctR = 0.5
                        , _cbR = 0.5
                        , _cL  = 0.5
                        , _ceR = 0.5 }
