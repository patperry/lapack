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
import LAPACK.CTypes

foreign import ccall unsafe "LAPACK.h lapack_dgeqrf"
    dgeqrf :: Int -> Int -> Ptr Double -> Int -> Ptr Double -> Ptr Double -> Int -> IO Int

foreign import ccall unsafe "LAPACK.h lapack_dgelqf"
    dgelqf :: Int -> Int -> Ptr Double -> Int -> Ptr Double -> Ptr Double -> Int -> IO Int

foreign import ccall unsafe "LAPACK.h lapack_dormqr"
    dormqr :: CBLASSide -> CBLASTrans -> Int -> Int -> Int -> Ptr Double -> Int -> Ptr Double
           -> Ptr Double -> Int -> Ptr Double -> Int -> IO Int

foreign import ccall unsafe "LAPACK.h lapack_dormlq"
    dormlq :: CBLASSide -> CBLASTrans -> Int -> Int -> Int -> Ptr Double -> Int -> Ptr Double
           -> Ptr Double -> Int -> Ptr Double -> Int -> IO Int

foreign import ccall unsafe "LAPACK.h lapack_dlarfg"
    dlarfg :: Int -> Ptr Double -> Ptr Double -> Int -> Ptr Double -> IO ()
