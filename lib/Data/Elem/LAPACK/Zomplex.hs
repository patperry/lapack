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
import LAPACK.CTypes

type Zomplex = Complex Double

foreign import ccall unsafe "LAPACK.h lapack_zgeqrf"
    zgeqrf :: Int -> Int -> Ptr Zomplex -> Int -> Ptr Zomplex -> Ptr Zomplex -> Int -> IO Int

foreign import ccall unsafe "LAPACK.h lapack_zgelqf"
    zgelqf :: Int -> Int -> Ptr Zomplex -> Int -> Ptr Zomplex -> Ptr Zomplex -> Int -> IO Int

foreign import ccall unsafe "LAPACK.h lapack_zunmqr"
    zunmqr :: CBLASSide -> CBLASTrans -> Int -> Int -> Int -> Ptr Zomplex -> Int -> Ptr Zomplex
           -> Ptr Zomplex -> Int -> Ptr Zomplex -> Int -> IO Int

foreign import ccall unsafe "LAPACK.h lapack_zunmlq"
    zunmlq :: CBLASSide -> CBLASTrans -> Int -> Int -> Int -> Ptr Zomplex -> Int -> Ptr Zomplex
           -> Ptr Zomplex -> Int -> Ptr Zomplex -> Int -> IO Int

foreign import ccall unsafe "LAPACK.h lapack_zlarfg"
    zlarfg :: Int -> Ptr Zomplex -> Ptr Zomplex -> Int -> Ptr Zomplex -> IO ()

