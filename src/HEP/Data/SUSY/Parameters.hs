module HEP.Data.SUSY.Parameters where

data ModelParams = ModelParams { _M0   :: Double
                               , _c    :: ModularWeights
                               , _tanb :: Double
                               } deriving Show

data ModularWeights = ModularWeights { _Hu :: Double
                                     , _Hd :: Double
                                     , _Q  :: Double
                                     , _tR :: Double
                                     , _bR :: Double
                                     , _L  :: Double
                                     , _eR :: Double
                                     } deriving Show

cos2Beta :: Double -> Double
cos2Beta tanb = (1.0 - tanbSq) / (1.0 + tanbSq)
  where tanbSq = tanb * tanb

cosBeta :: Double -> Double
cosBeta tanb = 1.0 / sqrt (1 + tanb * tanb)
