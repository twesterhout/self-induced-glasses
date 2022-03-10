{-# LANGUAGE OverloadedStrings #-}

module Main (main) where

import Control.Monad
import Control.Monad.Primitive
import Control.Monad.State.Strict
import Control.Monad.Trans
import qualified Data.HDF5 as H5
import Data.Text (Text, pack)
import Data.Vector (Vector)
import qualified Data.Vector.Generic as G
import SelfInducedGlasses
import SelfInducedGlasses.Analysis
import SelfInducedGlasses.Core
import SelfInducedGlasses.Metropolis
import SelfInducedGlasses.Random
import System.IO
import System.Random.Stateful
import Text.Printf (printf)

generateInitialState :: (G.Vector v ℝ, StatefulGen g m) => Int -> g -> m (v ℝ)
generateInitialState n g = undefined

β :: ℝ
β = 0.12

couplingsRenormalization :: IO ()
couplingsRenormalization = do
  let n = 10
      lattice = Lattice (n, n) squareLatticeVectors
  forM_ [2, 5] $ \y ->
    forM_ [2, 5, 20] $ \λ ->
      let model = Model lattice λ
          pⱼ = P 3 y
          filename = pack $ printf "couplings_renormalization_%f_%d.csv" λ y
       in dumpCouplingEvolutionToFile filename $ take 1000 (effectiveInteractionDebug model pⱼ)

experiment2 :: IO ()
experiment2 = do
  let inputPath = pack $ "experiment1_" <> show β <> ".h5"
  states <- loadStates inputPath
  couplings@(DenseMatrix n _ _) <- loadCouplings inputPath
  let configurations = fmap (unpackConfiguration n) (denseMatrixRows states)
      energy = fmap (energyPerSite couplings) configurations
      magnetization = fmap magnetizationPerSite configurations
      autocorrFunction = autocorrStates n states
  -- autocorrFunction100 = autocorr 50000 n states
  -- withFile ("energy_" <> show β <> ".dat") WriteMode $ \h ->
  --   forM_ energy (hPutStrLn h . show)
  -- withFile ("magnetization_" <> show β <> ".dat") WriteMode $ \h ->
  --   forM_ magnetization (hPutStrLn h . show)
  print $ integratedAutocorrTime (G.take 1000 autocorrFunction)
  print $ integratedAutocorrTime (G.take 10000 autocorrFunction)
  print $ integratedAutocorrTime (G.take 20000 autocorrFunction)
  print $ integratedAutocorrTime (G.take 30000 autocorrFunction)
  withFile ("autocorr_" <> show β <> ".dat") WriteMode $ \h ->
    G.forM_ autocorrFunction (hPutStrLn h . show)
-- withFile ("autocorr50000_" <> show β <> ".dat") WriteMode $ \h ->
--   G.forM_ autocorrFunction100 (hPutStrLn h . show)
{-# SCC experiment2 #-}

experiment1 :: IO ()
experiment1 = do
  let n = 8
      lattice = Lattice (n, n) squareLatticeVectors
      λ = 20.0
      model = Model lattice λ
      couplings = buildInteractionMatrix lattice (effectiveInteraction model 100)
      sweepSize = n * n
      numberSweeps = 100000
      thermalizationSweeps = 50000
      filename = pack $ "experiment1_" <> show β <> ".h5"
      options = SamplingOptions couplings Nothing
      sample β g = do
        lift $ putStrLn "Running Monte Carlo ..."
        _ <- thermalize β thermalizationSweeps sweepSize g
        (states, acceptance) <- manySweeps β numberSweeps sweepSize g
        lift $ putStrLn "Saving results ..."
        lift $ print acceptance
        lift $
          H5.withFile filename H5.WriteTruncate $ \h -> do
            H5.createDataset h "/couplings" couplings
            H5.createDataset h "/states" states

  -- forM_ rs $ \(MetropolisStats acceptance σ) ->
  --   let e = computeEnergyPerSite couplings σ
  --       m = computeMagnetizationPerSite σ
  --    in lift . lift $ do
  --         putStrLn $ show acceptance <> "\t" <> show e <> "\t" <> show m
  --         saveForGnuplot lattice "picture.dat" σ
  -- g = mkStdGen 42
  (g :: Xoshiro256PlusPlus (PrimState IO)) <- mkXoshiro256PlusPlus 42
  _ <- runMetropolisT (sample β) options g
  -- :: Xoshiro256PlusPlus (PrimState IO) -> IO ()) g
  -- _ <- runStateGenT g (runMetropolisT (sample β) options)
  pure ()
{-# SCC experiment1 #-}

main :: IO ()
main = do
  -- couplingsRenormalization
  experiment1

-- experiment2

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
