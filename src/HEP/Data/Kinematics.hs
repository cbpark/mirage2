module HEP.Data.Kinematics where

newtype Mass = Mass { getMass :: Double }

instance Show Mass where
    show = show . getMass

massSq :: Mass -> Double
massSq (Mass m) = m * m
