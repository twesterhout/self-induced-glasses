module Main (main) where

import qualified Data.Vector.Generic as G
import SelfInducedGlasses.Core
import SelfInducedGlasses.Metropolis
import System.Random

main :: IO ()
main = do
  putStrLn ("Test suite is not implemented" :: String)
  let g = mkStdGen 42
      simplify (Configuration v) = G.map ((\x -> (x + 1) `div` 2) . (round :: ℝ -> Int)) v
      (x, g') = randomConfiguration 100 g
  print (simplify x)
  let y = packConfiguration x
  print y
  let z = unpackConfiguration 100 y
  print (simplify z)
  print (z == x)
