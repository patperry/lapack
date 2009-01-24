{-# LANGUAGE FlexibleInstances, GADTs, MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.QR
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- QR factorization of a matrix
--

module Data.Matrix.QR (
    -- * The QR decomposition type
    QR(..),
    qrQ,
    qrR,
    
    -- * Computing QR factorizations
    qrFactor,
    getQRFactor,
    qrFactorize,

    -- * Overloaded interface for solving linear systems
    module Data.Matrix.Class.ISolve,
    module Data.Matrix.Class.MSolve,
    
    ) where

import Control.Monad.Interleave
import Control.Monad.ST( runST )
import Data.Maybe( fromJust )

import Data.Elem.BLAS
import Data.Elem.LAPACK.C

import Data.Matrix.Class.ISolve
import Data.Matrix.Class.MSolve
import Data.Matrix.Dense
import Data.Matrix.Dense.IO
import Data.Matrix.Dense.ST

import Data.Matrix.House
import Data.Matrix.Tri

import Data.Vector.Dense.IO
import Data.Vector.Dense.ST

import Unsafe.BLAS

-- | A QR factorization of a dense matrix.
data QR np e where
    QR :: House (n,n) e -> Tri Matrix (n,p) e -> QR (n,p) e

-- | Get the Q part of the factorization.    
qrQ :: QR (n,p) e -> House (n,n) e
qrQ (QR q _) = q
{-# INLINE qrQ #-}

-- | Get the R part of the factorization.
qrR :: QR (n,p) e -> (Tri Matrix (n,p) e)
qrR (QR _ r) = r
{-# INLINE qrR #-}

instance Shaped QR (Int,Int) where
    shape  (QR _ r) = shape r
    bounds (QR _ r) = bounds r

instance MatrixShaped QR

instance ISolve QR where
    unsafeSSolveVector k qr y = runSTVector $ unsafeGetSSolveVector k qr y
    unsafeSSolveMatrix k qr c = runSTMatrix $ unsafeGetSSolveMatrix k qr c

instance (MonadInterleave m) => MSolve QR m where
    unsafeDoSolveVector = unsafeDoQRSolveVector
    unsafeDoSolveMatrix = unsafeDoQRSolveMatrix

unsafeDoQRSolveVector :: (ReadVector y m, WriteVector x m, BLAS2 e)
                      => QR (n,p) e -> y n e -> x p e -> m ()
unsafeDoQRSolveVector (QR q r) y x = do
    qty <- unsafeGetSApplyVector 1 (herm q) y
    unsafeDoSolveVector r qty x
{-# INLINE unsafeDoQRSolveVector #-}

unsafeDoQRSolveMatrix :: (ReadMatrix c m, WriteMatrix b m, BLAS3 e)
                      => QR (n,p) e -> c (n,q) e -> b (p,q) e -> m ()
unsafeDoQRSolveMatrix (QR q r) c b = do
    qtc <- unsafeGetSApplyMatrix 1 (herm q) c
    unsafeDoSolveMatrix r qtc b
{-# INLINE unsafeDoQRSolveMatrix #-}

-- | Get the QR factorization of an immutable dense matrix.
qrFactor :: (LAPACK e) => Matrix (n,p) e -> QR (n,p) e
qrFactor a = runST $ getQRFactor a

-- | Get the QR factorization of a dense matrix.
getQRFactor :: (ReadMatrix a m, LAPACK e) => a (n,p) e -> m (QR (n,p) e)
getQRFactor a = unsafePerformIOWithMatrix a $ \a' -> do
    a'' <- newCopyMatrix a'
    qrFactorize a''
{-# INLINE getQRFactor #-}

-- | Compute the QR factorization of a matrix in-place and return the 
-- result.  Note that this destroys the values stored in the input matrix
-- and that the returned object assumes ownership of the memory.
qrFactorize :: (WriteMatrix a m, LAPACK e) => a (n,p) e -> m (QR (n,p) e)
qrFactorize a
    | isHermMatrix a =
        unsafePerformIOWithMatrix (herm a) $ \a' -> do
            tau <- newVector_ np
            withIOMatrix a' $ \pA ->
                withIOVector tau $ \pTau ->
                    gelqf n' p' pA ldA pTau
            a''  <- unsafeFreezeMatrix a'
            tau' <- unsafeFreezeVector tau
            if n' > p'
                then let q = fromJust $ flip maybeHouseFromRows tau' $ upperU $
                                 unsafeSubmatrix a'' (0,0) (p',p')
                         l = lower a''
                     in return $ QR (herm q) (herm l)
                else let q = fromJust $ flip maybeHouseFromRows tau' $ upperU a''
                         l = lower a''
                    in return $ QR (herm q) (herm l)
    | otherwise =
        unsafePerformIOWithMatrix a $ \a' -> do
            tau <- newVector_ np
            withIOMatrix a' $ \pA ->
                withIOVector tau $ \pTau ->
                    geqrf n p pA ldA pTau
            a''  <- unsafeFreezeMatrix a'
            tau' <- unsafeFreezeVector tau
            if n < p
                then let q = fromJust $ flip maybeHouseFromCols tau' $ lowerU $
                                 unsafeSubmatrix a'' (0,0) (n,n)
                         r = upper a''
                     in return $ QR q r
                else let q = fromJust $ flip maybeHouseFromCols tau' $ lowerU a''
                         r = upper a''
                     in return $ QR q r
  where (n,p)   = shape a
        (n',p') = if isHermMatrix a then (p,n) else (n,p)
        np      = min n p
        ldA     = ldaMatrix a
{-# INLINE qrFactorize #-}
