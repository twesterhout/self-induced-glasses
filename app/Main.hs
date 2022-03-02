module Main (main) where

import Control.Monad
import Control.Monad.State.Strict
import Control.Monad.Trans
import Data.Vector (Vector)
import qualified Data.Vector.Generic as G
import SelfInducedGlasses
import SelfInducedGlasses.Core
import SelfInducedGlasses.Metropolis
import System.IO
import System.Random.Stateful

generateInitialState :: (G.Vector v ℝ, StatefulGen g m) => Int -> g -> m (v ℝ)
generateInitialState n g = undefined

experiment1 :: IO ()
experiment1 = do
  let n = 10
      lattice = Lattice (n, n) squareLatticeVectors
      λ = 40.0
      model = Model lattice λ
      couplings = buildInteractionMatrix lattice (effectiveInteraction model 100)
      g = mkStdGen 42
      sweepSize = n * n
      numberSweeps = 1
      sample :: RandomGen g => ℝ -> MetropolisT (StateGenM g) (StateT g IO) ()
      sample β = do
        rs <- manySweeps β numberSweeps sweepSize
        forM_ rs $ \(MetropolisStats acceptance σ) ->
          let e = computeEnergyPerSite couplings σ
              m = computeMagnetizationPerSite σ
           in lift . lift $ do
                putStrLn $ show acceptance <> "\t" <> show e <> "\t" <> show m
                saveForGnuplot lattice "picture.dat" σ
      run :: RandomGen g => StateGenM g -> StateT g IO ()
      run g' = runMetropolisT options (sample 0.001)
        where
          options = SamplingOptions couplings g' Nothing
  _ <- runStateGenT g run
  pure ()

main :: IO ()
main = do
  experiment1

-- let lattice = Lattice (5, 5) squareLatticeVectors
--     model = Model lattice 3
--     couplings = buildInteractionMatrix lattice (effectiveInteraction model 100)
--     β = 0.5
--     g = mkStdGen 42
-- let run g' =
--       let options = SamplingOptions couplings β g' Nothing
--        in runMetropolisT options (sweep 100)
-- (MetropolisStats acceptance spins', _) <- runStateGenT g run
-- print spins'
-- print acceptance

-- withFile "resummation.dat" WriteMode $ \h ->
--   forM_ (take 2500 $ effectiveInteractionDebug testModel1 (P 2 1)) $ \(r, δj, j) ->
--     when (r `mod` 100 == 0) $
--       hPutStrLn h (show r <> "\t" <> show δj <> "\t" <> show j)
