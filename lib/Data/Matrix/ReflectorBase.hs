{-# LANGUAGE FlexibleInstances, GADTs, MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.ReflectorBase
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--

module Data.Matrix.ReflectorBase
    where

import Control.Monad
import Control.Monad.Interleave
import Control.Monad.ST
import Foreign
import Unsafe.Coerce

import Data.Elem.BLAS
import Data.Elem.LAPACK.C

import Data.Vector.Dense
import Data.Vector.Dense.ST( runSTVector )
import Data.Vector.Dense.Class
import Data.Matrix.Dense.Class
import Data.Matrix.Dense.ST( runSTMatrix )
import Data.Matrix.Class
import Data.Tensor.Class
import Data.Tensor.Class.MTensor( unsafeReadElem, unsafeWriteElem, 
    unsafeModifyElem )
import Unsafe.BLAS( IOVector(..), IMatrix(..), MMatrix(..), ISolve(..),
    MSolve(..), unsafePerformIOWithVector, unsafeIOVectorToVector,
    unsafeSubvectorView, unsafeCopyVector, unsafeAxpyVector, unsafeGetDot, 
    unsafeSubmatrixView, unsafeRowView, unsafeCopyMatrix,
    unsafeRank1UpdateMatrix, unsafeGetSSolveVector, unsafeGetSSolveMatrix )

-- | An elementary Housholder reflector.
data Reflector nn e where
    Reflector :: (BLAS1 e) => Vector n e -> e -> Reflector (n,n) e

coerceReflector :: Reflector np e -> Reflector np' e
coerceReflector = unsafeCoerce
{-# INLINE coerceReflector #-}

hermReflector :: Reflector (n,p) e -> Reflector (p,n) e
hermReflector (Reflector v tau) = Reflector v (conjugate tau)

unsafeDoSApplyReflectorVector_ :: (WriteVector x m, BLAS1 e)
                               => e -> Reflector (n,n) e -> x n e -> m ()
unsafeDoSApplyReflectorVector_ k (Reflector v tau) x
    | n == 0    = return ()
    | n == 1    = unsafeModifyElem x 0 (\e -> (1-tau)*(k*e))
    | otherwise =
        let x2 = unsafeSubvectorView x 1 (n-1)
        in do
            scaleBy k x
            x1  <- unsafeReadElem x 0
            vx2 <- unsafeGetDot v x2
            let alpha = -tau*(x1 + vx2)
            unsafeWriteElem x 0 (x1 + alpha)
            unsafeAxpyVector alpha v x2
  where
    n  = dim x
{-# INLINE unsafeDoSApplyReflectorVector_ #-}
    
unsafeDoSApplyReflectorMatrix_ :: (WriteMatrix a m, BLAS3 e)
                               => e -> Reflector (n,n) e -> a (n,p) e -> m ()
unsafeDoSApplyReflectorMatrix_ k (Reflector v tau) a
    | n == 0    = return ()
    | n == 1    = scaleBy ((1-tau)*k) a
    | otherwise =
        let a1 = unsafeRowView a 0
            a2 = unsafeSubmatrixView a (1,0) (n-1,p)
        in do
            scaleBy k a
            z <- newCopyVector (conj a1)
            unsafeDoSApplyAddVector 1 (herm a2) v 1 z
            unsafeAxpyVector (-tau) (conj z) a1
            unsafeRank1UpdateMatrix a2 (-tau) v z
  where (n,p) = shape a
{-# INLINE unsafeDoSApplyReflectorMatrix_ #-}

unsafeGetSApplyReflectorVector :: (ReadVector x m, WriteVector y m, BLAS1 e)
                               => e -> Reflector (n,p) e -> x p e -> m (y n e)
unsafeGetSApplyReflectorVector k r x = do
        x' <- newCopyVector x
        unsafeDoSApplyReflectorVector_ k (unsafeCoerce r) x'
        return $ unsafeCoerce x'

unsafeGetSApplyReflectorMatrix :: (ReadMatrix a m, WriteMatrix b m, BLAS3 e)
                               => e 
                               -> Reflector (n,p) e 
                               -> a (p,q) e
                               -> m (b (n,q) e)
unsafeGetSApplyReflectorMatrix k r a = do
        a' <- newCopyMatrix a
        unsafeDoSApplyReflectorMatrix_ k (unsafeCoerce r) a'
        return $ unsafeCoerce a'
{-# INLINE unsafeGetSApplyReflectorMatrix #-}
    
unsafeGetColReflector :: (WriteVector x m)
             => Reflector (n,p) e -> Int -> m (x n e)
unsafeGetColReflector r@(Reflector _ _) j = do
    e <- newBasisVector (numRows r) j
    unsafeDoSApplyReflectorVector_ 1 (unsafeCoerce r) e
    return e
{-# INLINE unsafeGetColReflector #-}
    
-- | Compute an elementary reflector @H@ such that @H' x = beta e_1@.  The
-- reflector H is represented as @H = I - tau (1; v) (1; v)'@.  The function
-- sets the first element of @x@ to @beta@, sets the rest of the components
-- to @v@, and returns @tau@.  The function returns the resulting reflector,
-- along with the value @beta@.  Note that this is a destructive operation
-- and that the returned object assumes ownership of the memory in @x@.
setReflector :: (WriteVector x m, LAPACK e) 
             => x n e 
             -> m (Reflector (n,n) e, e)
setReflector x | dim x == 0 = 
    fail $ "setReflector <vector of dim 0>: dimension must be positive."
setReflector x = unsafePerformIOWithVector x $ 
    \(IOVector c n f p inc) -> 
        let p1 = p `advancePtr` inc 
            v  = unsafeIOVectorToVector $ IOVector NoConj (n-1) f p1 inc
        in do
            when (c == Conj) $ doConjVector (IOVector c n f p inc)
            tau  <- larfg n p p1 inc
            beta <- peek p
            touchForeignPtr f
            return $ (Reflector v tau, beta)
{-# INLINE setReflector #-}

-- | Get an elementary reflector @H@ that transforms the vector @x@ to be
-- parallel with the first basis vector, along with value @(H' x)_1@.
getReflector :: (ReadVector x m, LAPACK e)
             => x n e
             -> m (Reflector (n,n) e, e)
getReflector x = unsafePerformIOWithVector x $ \x' -> do
    y <- newCopyVector x'
    setReflector y
{-# INLINE getReflector #-}

-- | Get an elementary reflector @H@ that transforms the vector @x@ to be
-- parallel with the first basis vector, along with value @(H' x)_1@.
reflector :: (LAPACK e) => Vector n e -> (Reflector (n,n) e, e)
reflector x = runST $ getReflector x
{-# INLINE reflector #-}

instance Shaped Reflector (Int,Int) where
    shape (Reflector v _) = (n,n) where n = 1 + dim v
    {-# INLINE shape #-}
    bounds a = ((0,0),(n-1,n-1)) where (n,_) = shape a
    {-# INLINE bounds #-}

instance MatrixShaped Reflector where
    herm = hermReflector
    {-# INLINE herm #-}

instance IMatrix Reflector where
    unsafeSApplyVector alpha a x = runSTVector $ unsafeGetSApplyVector alpha a x
    {-# INLINE unsafeSApplyVector #-}
    unsafeSApplyMatrix alpha a b = runSTMatrix $ unsafeGetSApplyMatrix alpha a b
    {-# INLINE unsafeSApplyMatrix #-}
    unsafeCol a j = runSTVector $ unsafeGetCol a j
    {-# INLINE unsafeCol #-}
    unsafeRow a i = runSTVector $ unsafeGetRow a i
    {-# INLINE unsafeRow #-}

instance (MonadInterleave m) => MMatrix Reflector m where
    unsafeGetSApplyVector = unsafeGetSApplyReflectorVector
    {-# INLINE unsafeGetSApplyVector #-}
    unsafeGetSApplyMatrix = unsafeGetSApplyReflectorMatrix
    {-# INLINE unsafeGetSApplyMatrix #-}
    unsafeDoSApplyVector_ = unsafeDoSApplyReflectorVector_
    {-# INLINE unsafeDoSApplyVector_ #-}
    unsafeDoSApplyMatrix_ = unsafeDoSApplyReflectorMatrix_
    {-# INLINE unsafeDoSApplyMatrix_ #-}
    unsafeGetCol = unsafeGetColReflector
    {-# INLINE unsafeGetCol #-}

instance ISolve Reflector where
    unsafeSSolveVector alpha a x = runSTVector $ unsafeGetSSolveVector alpha a x
    {-# INLINE unsafeSSolveVector #-}
    unsafeSSolveMatrix alpha a b = runSTMatrix $ unsafeGetSSolveMatrix alpha a b
    {-# INLINE unsafeSSolveMatrix #-}

instance (MonadInterleave m) => MSolve Reflector m where
    unsafeDoSSolveVector alpha a x y = do
        unsafeCopyVector (coerceVector y) x
        unsafeDoSApplyVector_ alpha (coerceReflector $ herm a) (coerceVector y)
    {-# INLINE unsafeDoSSolveVector #-}
    unsafeDoSSolveMatrix alpha a b c = do
        unsafeCopyMatrix (coerceMatrix c) b
        unsafeDoSApplyMatrix_ alpha (coerceReflector $ herm a) (coerceMatrix c)
    {-# INLINE unsafeDoSSolveMatrix #-}
    