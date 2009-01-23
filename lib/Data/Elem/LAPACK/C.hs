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

import Control.Exception( assert )
import Control.Monad
import Data.Complex( Complex )
import Data.Elem.BLAS.Level3
import Data.Matrix.Class( TransEnum(..), SideEnum(..) )
import Foreign
import LAPACK.CTypes

import Data.Elem.LAPACK.Double
import Data.Elem.LAPACK.Zomplex

-- | The LAPACK typeclass.
class (BLAS3 e) => LAPACK e where
    larfg :: Int -> Ptr e -> Ptr e -> Int -> IO e
    unmqr :: SideEnum -> TransEnum -> Int -> Int -> Int -> Ptr e -> Int -> Ptr e -> Ptr e -> Int -> IO ()
    unmlq :: SideEnum -> TransEnum -> Int -> Int -> Int -> Ptr e -> Int -> Ptr e -> Ptr e -> Int -> IO ()

callWithWork :: (Storable e) => (Ptr e -> Int -> IO a) -> IO a
callWithWork call =
    alloca $ \pQuery -> do
        call pQuery (-1)
        ldWork <- peek (castPtr pQuery) :: IO Double
        let lWork = ceiling ldWork
        allocaArray lWork $ \pWork -> do
            call pWork lWork

checkInfo :: Int -> IO ()
checkInfo info = assert (info == 0) $ return ()
        
instance LAPACK Double where
    larfg n alpha x incx = with 0 $ \pTau -> 
        dlarfg n alpha x incx pTau >> peek pTau
    unmqr s t m n k pA ldA pTau pC ldC =
        checkInfo =<< callWithWork (dormqr (cblasSide s) (cblasTrans t) m n k pA ldA pTau pC ldC)
    unmlq s t m n k pA ldA pTau pC ldC =
        checkInfo =<< callWithWork (dormlq (cblasSide s) (cblasTrans t) m n k pA ldA pTau pC ldC)

instance LAPACK (Complex Double) where
    larfg n alpha x incx = with 0 $ \pTau -> 
        zlarfg n alpha x incx pTau >> peek pTau
    unmqr s t m n k pA ldA pTau pC ldC =
        checkInfo =<< callWithWork (zunmqr (cblasSide s) (cblasTrans t) m n k pA ldA pTau pC ldC)
    unmlq s t m n k pA ldA pTau pC ldC =
        checkInfo =<< callWithWork (zunmlq (cblasSide s) (cblasTrans t) m n k pA ldA pTau pC ldC)
