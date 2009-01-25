{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module     : Data.Matrix.Perm
-- Copyright  : Copyright (c) , Patrick Perry <patperry@stanford.edu>
-- License    : BSD3
-- Maintainer : Patrick Perry <patperry@stanford.edu>
-- Stability  : experimental
--
-- Permutation matrices.
--

module Data.Matrix.Perm (
    -- * Permutation matrices
    Perm,
    identityPerm,
    
    -- * Converting to/from @Permute@s
    permFromPermute,
    permuteFromPerm,

    -- * Overloaded matrix interface
    module Data.Matrix.Class.IMatrix,
    module Data.Matrix.Class.MMatrix,
    module Data.Matrix.Class.ISolve,
    module Data.Matrix.Class.MSolve,
    ) where

import Control.Monad ( forM_ )
import Control.Monad.Interleave
import Data.AEq

import Data.Elem.BLAS
import Data.Matrix.Class.IMatrix
import Data.Matrix.Class.MMatrix
import Data.Matrix.Class.ISolve
import Data.Matrix.Class.MSolve

import Data.Matrix.Dense.ST
import Data.Vector.Dense.ST
import Unsafe.BLAS

import Data.Permute ( Permute )
import qualified Data.Permute as P

import Unsafe.Coerce

-- | The permutation matrix data type.
data Perm nn e = 
      Perm !TransEnum !Permute
    | Identity !Int

-- | Get an identity permutaiton matrix of the given size.
identityPerm :: Int -> Perm (n,n) e
identityPerm = Identity

-- | Get a matrix from a permutation.
permFromPermute :: Permute -> Perm (n,n) e
permFromPermute = Perm NoTrans

-- | Get a permutaiton from a matrix.
permuteFromPerm :: Perm (n,n) e -> Permute
permuteFromPerm (Identity n)   = P.permute n
permuteFromPerm (Perm h sigma) = 
    if h == ConjTrans then P.inverse sigma else sigma

coercePerm :: Perm np e -> Perm np' e
coercePerm = unsafeCoerce

instance Shaped Perm (Int,Int) where
    shape (Perm _ sigma) = (n,n) where n = P.size sigma
    shape (Identity n)   = (n,n)
    
    bounds p = ((0,0), (m-1,n-1)) where (m,n) = shape p
    
instance MatrixShaped Perm where
    
instance HasHerm Perm where
    herm (Perm h sigma) = Perm (flipTrans h) sigma
    herm a@(Identity _) = coercePerm a

instance IMatrix Perm where
    unsafeSApplyVector alpha a x = 
        runSTVector $ unsafeGetSApplyVector alpha a x
    unsafeSApplyMatrix alpha a b = 
        runSTMatrix $ unsafeGetSApplyMatrix alpha a b
    unsafeCol a j = runSTVector $ unsafeGetCol a j
    unsafeRow a i = runSTVector $ unsafeGetRow a i

instance (MonadInterleave m) => MMatrix Perm m where
    unsafeDoSApplyAddVector = unsafeDoSApplyAddPerm
    {-# INLINE unsafeDoSApplyAddVector #-}
    unsafeDoSApplyAddMatrix = unsafeDoSApplyAddMatPerm
    {-# INLINE unsafeDoSApplyAddMatrix #-}    
    unsafeDoSApplyVector_ = unsafeDoSApplyPerm_
    {-# INLINE unsafeDoSApplyVector_ #-}    
    unsafeDoSApplyMatrix_ = unsafeDoSApplyMatPerm_
    {-# INLINE unsafeDoSApplyMatrix_ #-}  
    unsafeGetCol = unsafeGetColPerm
    {-# INLINE unsafeGetCol #-}  

unsafeGetColPerm :: (WriteVector x m, Elem e) 
                 => Perm (n,p) e -> Int -> m (x n e)
unsafeGetColPerm (Identity n)   j = newBasisVector n j
unsafeGetColPerm (Perm h sigma) j
    | h == NoTrans = newBasisVector n (P.indexOf sigma j)
    | otherwise    = newBasisVector n (P.unsafeAt sigma j)
  where
    n = P.size sigma
{-# INLINE unsafeGetColPerm #-}

unsafeDoSApplyPerm_ :: (WriteVector y m, BLAS1 e) => 
    e -> Perm (k,k) e -> y k e -> m ()
unsafeDoSApplyPerm_ alpha (Identity _)   x = scaleBy alpha x
unsafeDoSApplyPerm_ alpha (Perm h sigma) x
    | h == ConjTrans = do
        scaleBy alpha x
        sequence_ [ swap i j | (i,j) <- P.invSwaps sigma ]
    | otherwise = do
        scaleBy alpha x
        sequence_ [ swap i j | (i,j) <- P.swaps sigma ]        
  where
    swap  = unsafeSwapElems x
{-# INLINE unsafeDoSApplyPerm_ #-}

unsafeDoSApplyMatPerm_ :: (WriteMatrix c m, BLAS1 e) => 
    e -> Perm (k,k) e -> c (k,l) e -> m ()
unsafeDoSApplyMatPerm_ alpha (Identity _)   a = scaleBy alpha a
unsafeDoSApplyMatPerm_ alpha (Perm h sigma) a
    | h == ConjTrans = do
        scaleBy alpha a
        sequence_ [ swap i j | (i,j) <- P.invSwaps sigma ]        
    | otherwise = do
        scaleBy alpha a
        sequence_ [ swap i j | (i,j) <- P.swaps sigma ]        
  where
    swap i j = unsafeSwapVector (unsafeRowView a i) (unsafeRowView a j)
{-# INLINE unsafeDoSApplyMatPerm_ #-}

unsafeDoSApplyAddPerm :: (ReadVector x m, WriteVector y m, BLAS1 e) =>
    e -> Perm (k,l) e -> x l e -> e -> y k e -> m ()
unsafeDoSApplyAddPerm alpha (Identity _) x beta y = do
    scaleBy beta y
    unsafeAxpyVector alpha (coerceVector x) y
unsafeDoSApplyAddPerm alpha p@(Perm h sigma) x beta y
    | isConj x =
        unsafeDoSApplyAddPerm (conjugate alpha) p (conj x) 
                              (conjugate beta) (conj y)
    | otherwise =
        let n     = dim x
        in do
            scaleBy beta y
            forM_ [0..(n-1)] $ \i ->
                let i' = P.unsafeAt sigma i in
                case h of
                   ConjTrans -> do
                        e <- unsafeReadElem x i  
                        f <- unsafeReadElem y i'
                        unsafeWriteElem y i' (alpha*e + f)
                   NoTrans  -> do
                        e <- unsafeReadElem x i'
                        f <- unsafeReadElem y i
                        unsafeWriteElem y i (alpha*e + f)
{-# INLINE unsafeDoSApplyAddPerm #-}

unsafeDoSApplyAddMatPerm :: (ReadMatrix b m, WriteMatrix c m, BLAS1 e) =>
    e -> Perm (r,s) e -> b (s,t) e -> e -> c (r,t) e -> m ()
unsafeDoSApplyAddMatPerm alpha (Identity _) b beta c = do
    scaleBy beta c
    unsafeAxpyMatrix alpha b (coerceMatrix c)
unsafeDoSApplyAddMatPerm alpha p@(Perm h sigma) b beta c =
    let m     = numCols p
    in do
        scaleBy beta c
        forM_ [0..(m-1)] $ \i ->
            let i' = P.unsafeAt sigma i in
            case h of
                ConjTrans -> unsafeAxpyVector alpha (unsafeRowView b i)
                                                    (unsafeRowView c i') 
                NoTrans   -> unsafeAxpyVector alpha (unsafeRowView b i') 
                                                    (unsafeRowView c i)
{-# INLINE unsafeDoSApplyAddMatPerm #-}

instance ISolve Perm where
    unsafeSSolveVector alpha a y = 
        runSTVector $ unsafeGetSSolveVector alpha a y
    unsafeSSolveMatrix alpha a c = 
        runSTMatrix $ unsafeGetSSolveMatrix alpha a c
    
instance (MonadInterleave m) => MSolve Perm m where    
    unsafeDoSSolveVector_ alpha p = 
        unsafeDoSApplyPerm_ alpha (herm p)
    {-# INLINE unsafeDoSSolveVector_ #-}
    unsafeDoSSolveMatrix_ alpha p = 
        unsafeDoSApplyMatPerm_ alpha (herm p)
    {-# INLINE unsafeDoSSolveMatrix_ #-}
    unsafeDoSSolveVector alpha p x y =
        unsafeDoSApplyAddPerm alpha (coercePerm $ herm p) x 0 y
    {-# INLINE unsafeDoSSolveVector #-}
    unsafeDoSSolveMatrix alpha p a b = 
        unsafeDoSApplyAddMatPerm alpha (coercePerm $ herm p) a 0 b
    {-# INLINE unsafeDoSSolveMatrix #-}

instance Show (Perm (n,n) e) where
    show   (Identity n) = "identityPerm " ++ show n
    show p@(Perm h sigma) 
        | h == ConjTrans = "herm (" ++ show (herm p) ++ ")"
        | otherwise      = "permFromPermute (" ++ show sigma ++ ")"
    
instance Eq (Perm (n,n) e) where
    (==) (Identity n) (Identity n') = n == n'
    (==) (Identity n) p@(Perm h _)
        | h == ConjTrans = (==) (Identity n) (herm p)
        | otherwise      = (==) (permFromPermute $ P.permute n) p
    (==) p (Identity n)  = (==) (Identity n) p

    (==) (Perm h sigma) (Perm h' sigma') 
        | h == h'   = sigma == sigma'
        | otherwise = P.size sigma == P.size sigma'
                      && sigma == (P.inverse sigma')

instance AEq (Perm (n,n) e) where
    (===) = (==)
    (~==) = (==)
    