{-# LANGUAGE FlexibleInstances, GADTs, MultiParamTypeClasses #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.HouseBase
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Products of elementary Householder reflectors.
--

module Data.Matrix.HouseBase
    where

import Foreign
import Text.Printf

import Control.Monad.Interleave
import Data.Elem.BLAS   
import Data.Elem.LAPACK.C
import Data.Vector.Dense
import Data.Vector.Dense.IO
import Data.Vector.Dense.ST
import Data.Matrix.Dense
import Data.Matrix.Dense.IO
import Data.Matrix.Dense.ST
import Data.Matrix.Tri
import Data.Matrix.Class
import Data.Tensor.Class
import Unsafe.BLAS
import Unsafe.Coerce

-- | An @n@-by-@n@ House matrix is represented as a product of @p@ 
-- elementary Householder reflectors, where @p <= n@.  Specifying the 
-- reflectors requires a unit-diagonal triangular matrix, @V@, and a
-- vector of multipliers, @tau@.  
-- 
-- The @i@th reflector is given by
--
--     @P[i] = I - tau[i] v[i] v[i]^H@, 
--
-- where @I@ is the @n@-by-@n@ identity matrix and @tau[i]@ is the @i@th
-- element of @tau@.  When @V@ is lower-triangular, @v[i]@ is its @i@th
-- column of @V@; otherwise, (when @V@ is upper-triangular) @v[i]@ is its
-- @i@th row.
--
-- We denote by @Q@ the product @P[1] ... P[p]@.  When the 
-- 'TransEnum' in the data constructor is 'NoTrans', we are representing
-- @Q@.  When it is 'ConjTrans', we are representing @Q^H@.
--
-- The triangular matrix must have unit diagonal.  Also, the underlying matrix
-- must have a 'transEnumMatrix' of 'NoTrans'.  Moreover, the vector of
-- multipliers must have a 'conjEnum' of 'NoConj'.
--
data House nn e where 
    ColumnWise :: (LAPACK e) => TransEnum -> Tri Matrix (n,p) e -> Vector p e -> House (n,n) e
    RowWise    :: (LAPACK e) => TransEnum -> Tri Matrix (p,n) e -> Vector p e -> House (n,n) e

-- When the matrix is lower-triangular, we use the 'ColumnWise' data
-- constructor.  Otherwise, we use 'RowWise'.

coerceHouse :: House nn e -> House nn' e
coerceHouse = unsafeCoerce

-- | Given a unit lower-triangular matrix storing reflectors as columns
-- and a vector of multipliers, possibly return a Householder matrix.  This
-- will fail if the matrix is hermed or the vector is conjugated.
maybeHouseFromCols :: (LAPACK e) 
                   => Tri Matrix (n,p) e -> Vector p e -> Maybe (House (n,n) e)
maybeHouseFromCols v@(Tri l u a) tau
    | numCols a /= dim tau =
        err "dimension mismatch"
    | l /= Lower =
        err "triangular matrix must be lower-triangular"
    | u /= Unit =
        err "triangular matrix must have unit diagonal"
    | otherwise =
        case (transEnumMatrix a, conjEnum tau) of
            (NoTrans, NoConj) -> Just $ ColumnWise NoTrans v tau
            _                 -> Nothing
  where
    err = error $ printf "maybeHouseFromCols <Tri %s %s <matrix of shape %s>> <vector of dim %d>"
                         (show l) (show u) (show (shape a)) (dim tau)
                         
                         
-- | Given a unit upper-triangular matrix storing reflectors as rows and a
-- vector of multipliers, possibly return a Householder matrix.  This will fail
-- if the matrix is hermed or the vector is conjugated.
maybeHouseFromRows :: (LAPACK e)
                   => Tri Matrix (p,n) e -> Vector p e -> Maybe (House (n,n) e)
maybeHouseFromRows v@(Tri l u a) tau
    | numRows a /= dim tau =
        err "dimension mismatch"
    | l /= Upper =
        err "triangular matrix must be upper-triangular"
    | u /= Unit =
        err "triangular matrix must have unit diagonal"
    | otherwise =
        case (transEnumMatrix a, conjEnum tau) of
            (NoTrans, NoConj) -> Just $ RowWise NoTrans v tau
            _                 -> Nothing
  where
    err = error $ printf "maybeHouseFromRows <Tri %s %s <matrix of shape %s>> <vector of dim %d>"
                         (show l) (show u) (show (shape a)) (dim tau)


unsafeDoSApplyVectorHouse_ :: (WriteVector x m) 
                           => e -> House (n,n) e -> x n e -> m ()
unsafeDoSApplyVectorHouse_ alpha a x = 
    unsafePerformIOWithVector x $ \iox ->
        case maybeViewVectorAsCol iox of
            Just c  -> unsafeDoSApplyMatrixHouse_ alpha a c
            Nothing -> 
                case maybeViewVectorAsRow (conj iox) of
                    Just c' -> unsafeDoSApplyMatrixHouse_ alpha a (herm c')
                    Nothing -> do
                        x' <- newCopyVector' iox
                        unsafeDoSApplyVectorHouse_ alpha a x'
                        unsafeCopyVector iox x'
{-# INLINE unsafeDoSApplyVectorHouse_ #-}

unsafeDoSApplyMatrixHouse_ :: (WriteMatrix a m) 
                           => e -> House (n,n) e -> a (n,p) e -> m ()
unsafeDoSApplyMatrixHouse_ alpha (ColumnWise trans (Tri _ _ a) tau) c =
    let (m,n)   = shape c
        k       = numCols a
        trans'  = if isHermMatrix c then flipTrans trans else trans
        side    = if isHermMatrix c then RightSide else LeftSide
        (m',n') = if isHermMatrix c then (n,m) else (m,n)
        ldA     = ldaMatrix a
        ldC     = ldaMatrix c
    in unsafePerformIOWithMatrix c $ \c' ->
       withIOMatrix (unsafeMatrixToIOMatrix a) $ \pA ->
       withIOVector (unsafeVectorToIOVector tau) $ \pTau ->
       withIOMatrix c' $ \pC -> do
           scaleByMatrix alpha c'
           unmqr side trans' m' n' k pA ldA pTau pC ldC
unsafeDoSApplyMatrixHouse_ alpha (RowWise trans (Tri _ _ a) tau) c =
    let (m,n)   = shape c
        k       = numRows a
        trans'  = if isHermMatrix c then flipTrans trans else trans
        side    = if isHermMatrix c then RightSide else LeftSide
        (m',n') = if isHermMatrix c then (n,m) else (m,n)
        ldA     = ldaMatrix a
        ldC     = ldaMatrix c
    in unsafePerformIOWithMatrix c $ \c' ->
       withIOMatrix (unsafeMatrixToIOMatrix a) $ \pA ->
       withIOVector (unsafeVectorToIOVector tau) $ \pTau ->
       withIOMatrix c' $ \pC -> do
           scaleByMatrix alpha c'
           unmlq side trans' m' n' k pA ldA pTau pC ldC
{-# INLINE unsafeDoSApplyMatrixHouse_ #-}

unsafeGetSApplyVectorHouse :: (ReadVector x m, WriteVector y m, BLAS1 e)
                               => e -> House (n,p) e -> x p e -> m (y n e)
unsafeGetSApplyVectorHouse k a x = do
        x' <- newCopyVector x
        unsafeDoSApplyVectorHouse_ k (unsafeCoerce a) x'
        return $ unsafeCoerce x'

unsafeGetSApplyMatrixHouse :: (ReadMatrix a m, WriteMatrix b m, BLAS3 e)
                               => e 
                               -> House (n,p) e 
                               -> a (p,q) e
                               -> m (b (n,q) e)
unsafeGetSApplyMatrixHouse k a b = do
        b' <- newCopyMatrix b
        unsafeDoSApplyMatrixHouse_ k (unsafeCoerce a) b'
        return $ unsafeCoerce b'
{-# INLINE unsafeGetSApplyMatrixHouse #-}
    
unsafeGetColHouse :: (WriteVector x m, Elem e)
             => House (n,p) e -> Int -> m (x n e)
unsafeGetColHouse a j = do
    e <- newBasisVector (numRows a) j
    unsafeDoSApplyVectorHouse_ 1 (unsafeCoerce a) e
    return e
{-# INLINE unsafeGetColHouse #-}

instance Shaped House (Int,Int) where
    shape (ColumnWise _ a _) = (n,n) where n = numRows a
    shape (RowWise    _ a _) = (n,n) where n = numCols a
    {-# INLINE shape #-}
    bounds a = case shape a of (m,n) -> ((0,0),(m-1,n-1))
    {-# INLINE bounds #-}

instance MatrixShaped House where
    herm (ColumnWise h a tau) = ColumnWise (flipTrans h) a tau
    herm (RowWise    h a tau) = RowWise    (flipTrans h) a tau
    {-# INLINE herm #-}

instance IMatrix House where
    unsafeSApplyVector alpha a x = runSTVector $ unsafeGetSApplyVector alpha a x
    {-# INLINE unsafeSApplyVector #-}
    unsafeSApplyMatrix alpha a b = runSTMatrix $ unsafeGetSApplyMatrix alpha a b
    {-# INLINE unsafeSApplyMatrix #-}
    unsafeCol a j = runSTVector $ unsafeGetCol a j
    {-# INLINE unsafeCol #-}
    unsafeRow a i = runSTVector $ unsafeGetRow a i
    {-# INLINE unsafeRow #-}

instance (MonadInterleave m) => MMatrix House m where
    unsafeGetSApplyVector = unsafeGetSApplyVectorHouse
    {-# INLINE unsafeGetSApplyVector #-}
    unsafeGetSApplyMatrix = unsafeGetSApplyMatrixHouse
    {-# INLINE unsafeGetSApplyMatrix #-}
    unsafeDoSApplyVector_ = unsafeDoSApplyVectorHouse_
    {-# INLINE unsafeDoSApplyVector_ #-}
    unsafeDoSApplyMatrix_ = unsafeDoSApplyMatrixHouse_
    {-# INLINE unsafeDoSApplyMatrix_ #-}
    unsafeGetCol = unsafeGetColHouse
    {-# INLINE unsafeGetCol #-}

instance ISolve House where
    unsafeSSolveVector alpha a x = runSTVector $ unsafeGetSSolveVector alpha a x
    {-# INLINE unsafeSSolveVector #-}
    unsafeSSolveMatrix alpha a b = runSTMatrix $ unsafeGetSSolveMatrix alpha a b
    {-# INLINE unsafeSSolveMatrix #-}

instance (MonadInterleave m) => MSolve House m where
    unsafeDoSSolveVector alpha a x y = do
        unsafeCopyVector (coerceVector y) x
        unsafeDoSApplyVector_ alpha (coerceHouse $ herm a) (coerceVector y)
    {-# INLINE unsafeDoSSolveVector #-}
    unsafeDoSSolveMatrix alpha a b c = do
        unsafeCopyMatrix (coerceMatrix c) b
        unsafeDoSApplyMatrix_ alpha (coerceHouse $ herm a) (coerceMatrix c)
    {-# INLINE unsafeDoSSolveMatrix #-}
