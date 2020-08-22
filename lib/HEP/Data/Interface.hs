{-# LANGUAGE DataKinds          #-}
{-# LANGUAGE DeriveGeneric      #-}
{-# LANGUAGE FlexibleInstances  #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE TypeOperators      #-}

module HEP.Data.Interface (InputArgs (..)) where

import Options.Generic

data InputArgs w = InputArgs
    { msusy  :: w ::: Double <?> "m(sfermion)"
    , output :: w ::: String <?> "the name of the output file"
    } deriving Generic

instance ParseRecord (InputArgs Wrapped)
deriving instance Show (InputArgs Unwrapped)
