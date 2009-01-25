{-# LANGUAGE ScopedTypeVariables #-}
module Orthogonal
    where

import Driver
import Monadic
import Test.QuickCheck hiding ( vector )
import Test.QuickCheck.BLAS( Pos(..), Nat2(..) )
import qualified Test.QuickCheck.BLAS as Test

import Control.Monad

import Data.Elem.BLAS
import Data.Vector.Dense
import Data.Vector.Dense.ST
import Data.Matrix.Dense
import Data.Matrix.Dense.ST
import Data.Matrix.House
import Data.Matrix.QR

import Debug.Trace

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
    [ ("snd . setReflector", mytest prop_setReflector_snd)
    , ("fst . reflector", mytest prop_reflector_fst)
    , ("reflector <*>", mytest prop_reflector_vector)
    , ("reflector <**>", mytest prop_reflector_matrix)
    , ("qrFactor", mytest prop_qrFactor)
    , ("qrFactor solveVector", mytest prop_qrFactor_solveVector)
    , ("qrFactor doSolveVector", mytest prop_qrFactor_doSSolveVector)    
    , ("qrFactor doSolveMatrix", mytest prop_qrFactor_doSSolveMatrix)
    ]
