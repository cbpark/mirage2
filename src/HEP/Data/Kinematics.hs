module HEP.Data.Kinematics where

newtype Mass = Mass { getMass :: Double } deriving (Eq, Ord)

instance Show Mass where
    show = show . getMass

massSq :: Mass -> Double
massSq (Mass m) = m * m
