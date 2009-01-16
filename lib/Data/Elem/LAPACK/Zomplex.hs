{-# LANGUAGE  ForeignFunctionInterface #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Elem.LAPACK.Zomplex
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Elem.LAPACK.Zomplex
    where

import Data.Complex( Complex )
import Foreign( Ptr )

type Zomplex = Complex Double

foreign import ccall unsafe "BLAS.h lapack_zlarfg"
    zlarfg :: Int -> Ptr Zomplex -> Ptr Zomplex -> Int -> Ptr Zomplex -> IO ()
