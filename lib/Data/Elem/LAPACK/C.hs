{-# LANGUAGE FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Elem.LAPACK.C
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Low-level interface to LAPACK.
--

module Data.Elem.LAPACK.C
    where

import Data.Complex( Complex )
import Data.Elem.BLAS.Level3
import Foreign

import Data.Elem.LAPACK.Double
import Data.Elem.LAPACK.Zomplex

-- | The LAPACK typeclass.
class (BLAS3 e) => LAPACK e where
    larfg :: Int -> Ptr e -> Ptr e -> Int -> IO e

instance LAPACK Double where
    larfg n alpha x incx = with 0 $ \pBeta -> 
        dlarfg n alpha x incx pBeta >> peek pBeta

instance LAPACK (Complex Double) where
    larfg n alpha x incx = with 0 $ \pBeta -> 
        zlarfg n alpha x incx pBeta >> peek pBeta
