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
    
    -- * Overloaded interface for matrices
    module Data.Matrix.Class,
    module Data.Matrix.Class.IMatrix,
    module Data.Matrix.Class.MMatrix,
    module Data.Matrix.Class.ISolve,
    module Data.Matrix.Class.MSolve,
    ) where

import Data.Matrix.ReflectorBase

import Data.Matrix.Class
import Data.Matrix.Class.IMatrix
import Data.Matrix.Class.MMatrix
import Data.Matrix.Class.ISolve
import Data.Matrix.Class.MSolve
