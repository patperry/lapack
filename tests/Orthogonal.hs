{-# LANGUAGE ScopedTypeVariables #-}
module Orthogonal
    where

import Driver
import Monadic
import Test.QuickCheck hiding ( vector )
import qualified Test.QuickCheck as QC
import Test.QuickCheck.BLAS( Index(..), Pos(..), Nat(..), Nat2(..) )
import qualified Test.QuickCheck.BLAS as Test

import Control.Monad

import Data.Permute( Permute )
import qualified Data.Permute as P

import Data.Elem.BLAS
import Data.Vector.Dense
import Data.Vector.Dense.ST
import Data.Matrix.Dense
import Data.Matrix.Dense.ST
import Data.Matrix.House
import Data.Matrix.Perm
import Data.Matrix.QR


testPermute :: Int -> Gen Permute
testPermute n = do
    es <- QC.vector n :: Gen [Int]
    return $ P.order n es

testPerm :: Int -> Gen (Perm (n,n) e)
testPerm n = do
    p <- liftM permFromPermute $ testPermute n
    elements [ identityPerm n, p, herm p ]
    
prop_perm_herm (Nat n) =
    forAll (testPerm n) $ \p ->
        permuteFromPerm (herm p) == P.inverse (permuteFromPerm p)
        
prop_perm_col (Index n i) =
    forAll (testPerm n) $ \p ->
        col p i
        === 
        p <*> (basisVector n i :: V) 

prop_perm_apply_basis (Index n i) =
    forAll (testPerm n) $ \p ->
        p <*> (basisVector n i :: V) 
        === 
        basisVector n (P.indexOf (permuteFromPerm p) i)

prop_perm_herm_apply (Nat n) =
    forAll (testPerm n)    $ \p ->
    forAll (Test.vector n) $ \(x :: V) ->
        p <*> herm p <*> x === x

prop_herm_perm_apply (Nat n) =
    forAll (testPerm n)    $ \p ->
    forAll (Test.vector n) $ \(x :: V) ->
        herm p <*> p <*> x === x

prop_perm_solve (Nat n) =
    forAll (testPerm n)    $ \p ->
    forAll (Test.vector n) $ \(x :: V) ->
        p <\> x === herm p <*> x

prop_perm_applyMat_cols (Nat2 (m,n)) =
    forAll (testPerm m) $ \p ->
    forAll (Test.matrix (m,n)) $ \(a :: M) ->
        cols (p <**> a) === map (p <*>) (cols a)

prop_perm_herm_applyMat (Nat2 (m,n)) =
    forAll (testPerm m) $ \p ->
    forAll (Test.matrix (m,n)) $ \(a :: M) ->
        p <**> herm p <**> a === a

prop_herm_perm_applyMat (Nat2 (m,n)) =
    forAll (testPerm m) $ \p ->
    forAll (Test.matrix (m,n)) $ \(a :: M) ->
        herm p <**> p <**> a === a

prop_perm_solveMat_cols (Nat2 (m,n)) =
    forAll (testPerm m) $ \p ->
    forAll (Test.matrix (m,n)) $ \(a :: M) ->
        cols (p <\\> a) === map (p <\>) (cols a)

prop_perm_solveMat (Nat2 (m,n)) =
    forAll (testPerm m) $ \p ->
    forAll (Test.matrix (m,n)) $ \(a :: M) ->
        p <\\> a === herm p <**> a

prop_perm_doSApplyVector_ alpha (Nat n) =
    forAll (testPerm n) $ \p ->
    forAll (Test.vector n) $ \(x :: V) ->
        monadicST $ do
            x'  <- run $ unsafeThawVector x
            x'' <- run $ freezeVector x'
            run $ doSApplyVector_ alpha p x'
            assert $ x ~== p <*> (alpha *> x'')

prop_perm_doSApplyMatrix_ alpha (Nat2 (m,n)) =
    forAll (testPerm m) $ \p ->
    forAll (Test.matrix (m,n)) $ \(b :: M) ->
        monadicST $ do
            b'  <- run $ unsafeThawMatrix b
            b'' <- run $ freezeMatrix b'
            run $ doSApplyMatrix_ alpha p b'
            assert $ b ~== p <**> (alpha *> b'')

prop_perm_doSSolveVector alpha (Nat n) =
    forAll (testPerm n) $ \p ->
    forAll (Test.vector n) $ \(x :: V) ->
    forAll (Test.vector n) $ \y -> 
        monadicST $ do
            x' <- run $ unsafeThawVector x
            run $ doSSolveVector alpha p y x'
            assert $ x ~== p <\> (alpha *> y)

prop_perm_doSSolveMatrix alpha (Nat2 (m,n)) =
    forAll (testPerm m) $ \p ->
    forAll (Test.matrix (m,n)) $ \(b :: M) ->
    forAll (Test.matrix (m,n)) $ \c -> 
        monadicST $ do
            b' <- run $ unsafeThawMatrix b
            run $ doSSolveMatrix alpha p c b'
            assert $ b ~== p <\\> (alpha *> c)

prop_setReflector_snd (Pos n) = 
    monadicST $ do
        (x :: V) <- pick (Test.vector n)
        x' <- run $ freezeVector x
        (_,beta) <- run $ setReflector =<< unsafeThawVector x
        assert $ norm beta ~== norm2 x'
        
prop_reflector_fst (Pos n) =
    forAll (Test.vector n) $ \(x :: V) ->
        let (r,beta) = reflector x
        in herm r <*> x ~== beta *> basisVector n 0

prop_reflector_vector (Pos n) =
    forAll (Test.vector n) $ \(x :: V) ->
        let (_,beta) = reflector x
        in norm beta ~== norm2 x

prop_reflector_matrix (Pos n) =
    forAll (Test.vector n) $ \(x :: V) ->
    forAll (liftM snd $ Test.shape) $ \p ->
    forAll (Test.elems p) $ \es ->
        let (r,beta) = reflector x
            a        = colsMatrix (n,p) [ k *> x | k <- es ]
            ra       = colsMatrix (n,p) [ (k*beta) *> basisVector n 0 | k <- es ]
        in herm r <**> a ~== ra

prop_qrFactor (Nat2 mn) =
    forAll (Test.matrix mn) $ \(a :: M) ->
        let qr    = qrFactor a
            (q,r) = (qrQ qr, qrR qr)
            i     = identityMatrix (numCols r, numCols r)
            a'    = q <**> r <**> i
        in a' ~== a

prop_qrFactor_solveVector (Nat2 mn) =
    forAll (Test.matrix mn) $ \(a :: M) ->
        let a' = runSTMatrix (do
                     ma <- unsafeThawMatrix a
                     setConstant 1 (diagView ma 0)
                     return ma) in
        forAll (Test.vector $ snd mn) $ \x ->
            let y  = a' <*> x
                x' = qrFactor a' <\> y
                y' = a' <*> x'
            in y' ~== y

prop_qrFactor_doSSolveVector alpha (Nat2 (n,p)) =
    forAll (Test.matrix (n,p)) $ \(a :: M) ->
        let a' = runSTMatrix (do
                     ma <- unsafeThawMatrix a
                     setConstant 1 (diagView ma 0)
                     return ma) in
        forAll (Test.vector p) $ \x ->
        forAll (Test.vector n) $ \y -> 
            let qr = qrFactor a' in 
            monadicST $ do
                x' <- run $ unsafeThawVector x
                run $ doSSolveVector alpha qr y x'
                assert $ x ~== qr <\> (alpha *> y)

prop_qrFactor_doSSolveMatrix alpha (Nat2 (n,p)) =
    forAll (Test.matrix (n,p)) $ \(a :: M) ->
        let a' = runSTMatrix (do
                     ma <- unsafeThawMatrix a
                     setConstant 1 (diagView ma 0)
                     return ma) in
        forAll (liftM snd Test.shape) $ \q ->
        forAll (Test.matrix (p,q)) $ \b ->
        forAll (Test.matrix (n,q)) $ \c -> 
            let qr = qrFactor a' in 
            monadicST $ do
                b' <- run $ unsafeThawMatrix b
                run $ doSSolveMatrix alpha qr c b'
                assert $ b ~== qr <\\> (alpha *> c)

tests_Orthogonal =
    [ 
      ("perm herm", mytest prop_perm_herm)
    , ("perm col", mytest prop_perm_col)
    , ("perm apply basis", mytest prop_perm_apply_basis)
    , ("perm herm apply", mytest prop_perm_herm_apply)
    , ("herm perm apply", mytest prop_herm_perm_apply)
    , ("perm solve", mytest prop_perm_solve)
    , ("perm applyMat cols", mytest prop_perm_applyMat_cols)
    , ("perm herm applyMat", mytest prop_perm_herm_applyMat)
    , ("herm perm applyMat", mytest prop_herm_perm_applyMat)
    , ("perm solveMat cols", mytest prop_perm_solveMat_cols)
    , ("perm solveMat", mytest prop_perm_solveMat)
    , ("perm doApplyVector_", mytest prop_perm_doSApplyVector_)
    , ("perm doApplyMatrix_", mytest prop_perm_doSApplyMatrix_)
    , ("perm doSolveVector", mytest prop_perm_doSSolveVector)
    , ("perm doSolveMatrix", mytest prop_perm_doSSolveMatrix)
    
    , ("snd . setReflector", mytest prop_setReflector_snd)
    , ("fst . reflector", mytest prop_reflector_fst)
    , ("reflector <*>", mytest prop_reflector_vector)
    , ("reflector <**>", mytest prop_reflector_matrix)
    , ("qrFactor", mytest prop_qrFactor)
    , ("qrFactor solveVector", mytest prop_qrFactor_solveVector)
    , ("qrFactor doSolveVector", mytest prop_qrFactor_doSSolveVector)    
    , ("qrFactor doSolveMatrix", mytest prop_qrFactor_doSSolveMatrix)
    ]
