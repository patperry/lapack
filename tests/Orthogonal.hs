{-# LANGUAGE ScopedTypeVariables #-}
module Orthogonal
    where

import Driver
import Test.QuickCheck hiding ( vector )
import Test.QuickCheck.BLAS( Pos(..) )
import qualified Test.QuickCheck.BLAS as Test

import Data.Vector.Dense
import Data.Matrix.House

prop_reflector (Pos n) =
    forAll (Test.vector n) $ \(x :: V) ->
        let (r,alpha) = reflector x
        in r <*> x ~== alpha *> basisVector n 0

tests_Orthogonal =
    [ ("reflector", mytest prop_reflector)
    ]