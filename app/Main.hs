{-# LANGUAGE LambdaCase #-}
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
import qualified Data.Vector.Storable as S
import qualified ListT
import SelfInducedGlasses
import SelfInducedGlasses.Analysis
import SelfInducedGlasses.Core
import SelfInducedGlasses.Interaction
import SelfInducedGlasses.Metropolis
import SelfInducedGlasses.Random
import System.IO
import System.Random.Stateful
import Text.Printf (printf)

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
  let n = (10 :: Int)
      λ = (20 :: ℝ)
      filename = pack $ printf "data/annealing_result_n=%d_λ=%f.h5" n λ
  H5.withFile filename H5.ReadOnly $ \h -> do
    (couplings :: DenseMatrix S.Vector ℝ) <- H5.open h "couplings" >>= H5.readDataset
    ListT.traverse_ pure $
      H5.forGroupM h $ \case
        g@H5.Group -> do
          states <- H5.open g "states" >>= H5.readDataset
          (H5.Scalar (β :: ℝ)) <- H5.readAttribute g "β"
          liftIO $ computeLocalObservables couplings states (pack $ printf "data/observables_n=%d_λ=%f_β=%f" n λ β)
        _ -> pure ()
-- computeLocalObservables :: DenseMatrix S.Vector ℝ -> DenseMatrix S.Vector Word64 -> Text -> IO ()

-- let β = (0.12 :: ℝ)
--     inputPath = pack $ "experiment1_" <> show β <> ".h5"
-- states <- loadStates inputPath
-- couplings@(DenseMatrix n _ _) <- loadCouplings inputPath
-- let configurations = fmap (unpackConfiguration n) (denseMatrixRows states)
--     energy = fmap (energyPerSite couplings) configurations
--     magnetization = fmap magnetizationPerSite configurations
--     autocorrFunction = autocorrStates n states
-- -- autocorrFunction100 = autocorr 50000 n states
-- -- withFile ("energy_" <> show β <> ".dat") WriteMode $ \h ->
-- --   forM_ energy (hPutStrLn h . show)
-- -- withFile ("magnetization_" <> show β <> ".dat") WriteMode $ \h ->
-- --   forM_ magnetization (hPutStrLn h . show)
-- print $ integratedAutocorrTime (G.take 1000 autocorrFunction)
-- print $ integratedAutocorrTime (G.take 10000 autocorrFunction)
-- print $ integratedAutocorrTime (G.take 20000 autocorrFunction)
-- print $ integratedAutocorrTime (G.take 30000 autocorrFunction)
-- withFile ("autocorr_" <> show β <> ".dat") WriteMode $ \h ->
--   G.forM_ autocorrFunction (hPutStrLn h . show)
-- withFile ("autocorr50000_" <> show β <> ".dat") WriteMode $ \h ->
--   G.forM_ autocorrFunction100 (hPutStrLn h . show)
{-# SCC experiment2 #-}

experiment1 :: IO ()
experiment1 = do
  let n = 10
      lattice = Lattice (n, n) squareLatticeVectors
      λ = 20.0
      model = Model lattice λ
      !couplings = buildCouplings model
      options = SamplingOptions couplings Nothing
      annealingSteps =
        [ (0.10, 10000, 10000),
          (0.11, 10000, 10000),
          (0.12, 10000, 10000),
          (0.13, 10000, 10000),
          (0.14, 10000, 10000),
          (0.15, 10000, 10000)
        ]
      filename = pack $ printf "data/annealing_result_n=%d_λ=%f.h5" n λ
      sample g = do
        lift $ putStrLn "[*] Running Monte Carlo ..."
        results <- anneal annealingSteps g
        -- _ <- thermalize β thermalizationSweeps sweepSize g
        -- (states, acceptance) <- manySweeps β numberSweeps sweepSize g
        lift $ putStrLn "[*] Saving results ..."
        -- lift $ print acceptance
        lift $
          H5.withFile filename H5.WriteTruncate $ \h -> do
            H5.createDataset h "couplings" couplings
            H5.open h "couplings" >>= \(d :: H5.Dataset) -> H5.writeAttribute d "λ" (H5.Scalar λ)
            forM_ (zip annealingSteps results) $ \((β, numberThermalization, _), (states, acceptance)) -> do
              group <- H5.createGroup h (pack $ printf "%f" β)
              H5.createDataset group "states" states
              H5.writeAttribute group "acceptance" (H5.Scalar acceptance)
              H5.writeAttribute group "thermalization" (H5.Scalar numberThermalization)
              H5.writeAttribute group "β" (H5.Scalar β)
  -- H5.createDataset group "acceptance" (H5.Scalar acceptance)

  -- forM_ rs $ \(MetropolisStats acceptance σ) ->
  --   let e = computeEnergyPerSite couplings σ
  --       m = computeMagnetizationPerSite σ
  --    in lift . lift $ do
  --         putStrLn $ show acceptance <> "\t" <> show e <> "\t" <> show m
  --         saveForGnuplot lattice "picture.dat" σ
  -- g = mkStdGen 42
  (g :: Xoshiro256PlusPlus (PrimState IO)) <- mkXoshiro256PlusPlus 42
  _ <- runMetropolisT sample options g
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
