{-# LANGUAGE CPP #-}
module Driver (
    E,
    V,
    M,
    
    AEq,
    (===),
    (~==),

    field,

    mytest,
    mycheck,
    mytests,
    done,

    ) where

import Data.AEq( AEq )
import qualified Data.AEq as AEq
import Data.Complex
import Data.List
import Debug.Trace
import System.IO
import System.Random
import Test.QuickCheck
import Text.Printf
import Text.Show.Functions

import Data.Vector.Dense( Vector )
import Data.Matrix.Dense( Matrix )

#ifdef COMPLEX
field = "Complex Double"
type E = Complex Double 
#else
field = "Double"
type E = Double
#endif

type V = Vector Int E
type M = Matrix (Int,Int) E

infix 4 ===, ~==

x === y | (AEq.===) x y = True
        | otherwise = trace (printf "expected `%s', but got `%s'" (show y) (show x)) False

x ~== y | (AEq.~==) x y = True
        | otherwise = trace (printf "expected `%s', but got `%s'" (show y) (show x)) False


------------------------------------------------------------------------
--
-- QC driver ( taken from xmonad-0.6 )
--

debug = False

mytest :: Testable a => a -> Int -> IO (Bool, Int)
mytest a n = mycheck defaultConfig
    { configMaxTest=n
    , configEvery   = \n args -> let s = show n in s ++ [ '\b' | _ <- s ] } a
 -- , configEvery= \n args -> if debug then show n ++ ":\n" ++ unlines args else [] } a

mycheck :: Testable a => Config -> a -> IO (Bool, Int)
mycheck config a = do
    rnd <- newStdGen
    mytests config (evaluate a) rnd 0 0 []

mytests :: Config -> Gen Result -> StdGen -> Int -> Int -> [[String]] -> IO (Bool, Int)
mytests config gen rnd0 ntest nfail stamps
    | ntest == configMaxTest config = done "OK," ntest stamps >> return (True, ntest)
    | nfail == configMaxFail config = done "Arguments exhausted after" ntest stamps >> return (True, ntest)
    | otherwise               =
      do putStr (configEvery config ntest (arguments result)) >> hFlush stdout
         case ok result of
           Nothing    ->
             mytests config gen rnd1 ntest (nfail+1) stamps
           Just True  ->
             mytests config gen rnd1 (ntest+1) nfail (stamp result:stamps)
           Just False ->
             putStr ( "Falsifiable after "
                   ++ show ntest
                   ++ " tests:\n"
                   ++ unlines (arguments result)
                    ) >> hFlush stdout >> return (False, ntest)
     where
      result      = generate (configSize config ntest) rnd2 gen
      (rnd1,rnd2) = split rnd0

done :: String -> Int -> [[String]] -> IO ()
done mesg ntest stamps = putStr ( mesg ++ " " ++ show ntest ++ " tests" ++ table )
  where
    table = display
            . map entry
            . reverse
            . sort
            . map pairLength
            . group
            . sort
            . filter (not . null)
            $ stamps

    display []  = ".\n"
    display [x] = " (" ++ x ++ ").\n"
    display xs  = ".\n" ++ unlines (map (++ ".") xs)

    pairLength xss@(xs:_) = (length xss, xs)
    entry (n, xs)         = percentage n ntest
                       ++ " "
                       ++ concat (intersperse ", " xs)

    percentage n m        = show ((100 * n) `div` m) ++ "%"

------------------------------------------------------------------------

