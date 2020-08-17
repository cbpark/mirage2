module HEP.Data.SUSY.Parameters where

-- data ModelParams = ModelParams { _M0   :: Double
--                                , _c    :: ModularWeights
--                                , _tanb :: Double
--                                } deriving Show

data ModularWeights = ModularWeights { _cHu :: Double
                                     , _cHd :: Double
                                     , _cQ  :: Double
                                     , _ctR :: Double
                                     , _cbR :: Double
                                     , _cL  :: Double
                                     , _ceR :: Double
                                     } deriving Show

cos2Beta :: Double -> Double
cos2Beta tanb = (1.0 - tanbSq) / (1.0 + tanbSq)
  where tanbSq = tanb * tanb

cosBeta :: Double -> Double
cosBeta tanb = 1.0 / sqrt (1 + tanb * tanb)
