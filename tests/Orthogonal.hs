{-# LANGUAGE ScopedTypeVariables #-}
module Orthogonal
    where

import Driver
import Monadic
import Test.QuickCheck hiding ( vector )
import Test.QuickCheck.BLAS( Pos(..) )
import qualified Test.QuickCheck.BLAS as Test

import Control.Monad

import Data.Elem.BLAS
import Data.Vector.Dense
import Data.Vector.Dense.ST
import Data.Matrix.Dense
import Data.Matrix.House

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

tests_Orthogonal =
    [ ("snd . setReflector", mytest prop_setReflector_snd)
    , ("fst . reflector", mytest prop_reflector_fst)
    , ("reflector <*>", mytest prop_reflector_vector)
    , ("reflector <**>", mytest prop_reflector_matrix)
    ]