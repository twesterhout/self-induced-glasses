module Main (main) where

import Control.Monad
import SelfInducedGlasses
import System.IO

main :: IO ()
main = do
  withFile "resummation.dat" WriteMode $ \h ->
    forM_ (take 2500 $ effectiveInteractionDebug testModel1 (P 2 1)) $ \(r, δj, j) ->
      when (r `mod` 100 == 0) $
        hPutStrLn h (show r <> "\t" <> show δj <> "\t" <> show j)
