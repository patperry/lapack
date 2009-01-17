-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.House
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Householder reflections.
--

module Data.Matrix.House (
    -- * Elementary reflectors
    Reflector,
    reflector,
    getReflector,
    setReflector,
    
    ) where

import Data.Matrix.ReflectorBase
