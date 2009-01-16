{-# LANGUAGE ForeignFunctionInterface #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Elem.LAPACK.Double
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Elem.LAPACK.Double
    where

import Foreign( Ptr )

foreign import ccall unsafe "BLAS.h lapack_dlarfg"
    dlarfg :: Int -> Ptr Double -> Ptr Double -> Int -> Ptr Double -> IO ()
