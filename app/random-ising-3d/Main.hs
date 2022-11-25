{-# LANGUAGE DuplicateRecordFields #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE OverloadedStrings #-}

module Main (main) where

import qualified Control.Foldl as Foldl
import Control.Monad (forM_)
import Control.Monad.ST (RealWorld)
import Data.List (intersperse)
import Data.Stream.Monadic (Step (..), Stream (..))
import qualified Data.Stream.Monadic as Stream
import Data.Text (pack)
import Data.Vector (Vector)
import Data.Vector.Generic ((!))
import qualified Data.Vector.Generic as G
import qualified ListT
import Options.Applicative
import SelfInducedGlasses.Random
import SelfInducedGlasses.Sampling
import System.Random.Stateful (freezeGen, thawGen, uniformM)
import Text.Printf (printf)

data MainSettings = MainSettings {size :: !Int, seed :: !Int, sweeps :: !Int}

pMainSettings :: Parser MainSettings
pMainSettings =
  MainSettings
    <$> option auto (long "size" <> help "Side length")
    <*> option auto (long "seed" <> help "Random number generator seed")
    <*> option auto (long "sweeps" <> help "Number sweeps between measurements")

run :: MainSettings -> IO ()
run (MainSettings n seed numSweeps) = do
  let -- ts = G.fromList [0.84999996, 0.95815915, 0.96502686, 0.97176474, 0.978413, 0.9850579, 0.9917415, 0.99841505, 1.0049784, 1.0113064, 1.017554, 1.0241934, 1.030724, 1.0371666, 1.0437058, 1.0504944, 1.0576115, 1.0650607, 1.0727555, 1.0805116, 1.2]
      -- ts = G.fromList [0.6, 0.7249834, 0.8357433, 0.9394234, 1.0331062, 1.1203706, 1.2068709, 1.2935127, 1.3753805, 1.4549646, 1.536332, 1.6221249, 1.7133656, 1.8060983, 1.9043283, 2.0104682, 2.1263938, 2.2658665, 2.43726, 2.6608124, 3.0]
      -- 2D Ising model
      -- ts = G.fromList [2.0e-1, 2.9529545e-1, 4.4920352e-1, 6.530335e-1, 8.7971103e-1, 1.0804988e0, 1.2323471e0, 1.3524361e0, 1.4534732e0, 1.5459211e0, 1.6282915e0, 1.7041376e0, 1.7745655e0, 1.8454263e0, 1.916044e0, 1.9817193e0, 2.0396886e0, 2.0944474e0, 2.1460757e0, 2.1929586e0, 2.2347448e0, 2.2718368e0, 2.3065984e0, 2.3419597e0, 2.381181e0, 2.4252942e0, 2.4755163e0, 2.5366852e0, 2.6128287e0, 2.7022533e0, 2.8031707e0, 2.9164104e0, 3.0429513e0, 3.183308e0, 3.3447561e0, 3.5340526e0, 3.755327e0, 4.0079594e0, 4.2939982e0, 4.614977e0, 5.0e0]
      -- 3D random Ising model
      -- ts = G.fromList [0.2, 0.25510567, 0.30581513, 0.3524917, 0.39541757, 0.43577713, 0.4751937, 0.5144822, 0.55365485, 0.5926007, 0.6310563, 0.668978, 0.7065402, 0.7439146, 0.7813722, 0.8203625, 0.85991955, 0.90033203, 0.9439652, 0.988771, 1.0348577, 1.0828927, 1.1327728, 1.1871444, 1.2472969, 1.3144703, 1.3901908, 1.4752965, 1.5737754, 1.6874822, 1.8182957, 1.9684793, 2.1413357, 2.3426905, 2.5771008, 2.845036, 3.1539936, 3.514956, 3.9378889, 4.4272513, 5.0]
      -- ts = [0.1,0.12819228,0.18958376,0.28366432,0.38042215,0.4772853,0.5505419,0.6222352,0.6933297,0.76180136,0.82273734,0.88004243,0.94099456,1.0189248,1.0968534,1.1777678,1.2509673,1.319273,1.3786286,1.4336286,1.4838532,1.5284123,1.5686824,1.6088673,1.6481054,1.6878262,1.7298468,1.7803172,1.8372859,1.9006613,2.0]
      ts = G.map (\β -> 1 / β) βs
      -- βs = G.map (\t -> 1 / t) ts
      βs = betasGeometric 17 0.2 2
      numReplicas = G.length βs

  g <- mkXoshiro256PlusPlus seed
  h <- randomIsingModel3D n <$> uniformM g
  -- let h = ferromagneticIsingModelSquare2D n 0
  -- print h

  let k = 10
      go (i :: Int) ts
        | i >= k = do
            writeVectorToCsv (pack $ printf "temperatures_%02d.csv" i) ts
            pure ts
        | otherwise = do
            writeVectorToCsv (pack $ printf "temperatures_%02d.csv" i) ts
            (ObservableState replicas) <-
              Stream.last . Stream.take 2 $
                monteCarloSampling
                  {- sweepSize -} (n * n * n)
                  {- numberSweeps -} ((2 ^ i) * numSweeps)
                  h
                  (G.map (\t -> 1 / t) ts)
                  g
            let stats = G.map (\r -> r.hist) replicas
                ts' = updateTemperatures ts stats
                accept = G.map (\r -> realToFrac r.stats.ssAcceptProb) replicas
                getFlow r = fromIntegral r.hist.numBlue / fromIntegral (r.hist.numRed + r.hist.numBlue)
                flow = G.map getFlow replicas
            print stats
            print ts'
            writeVectorToCsv (pack $ printf "acceptance_%02d.csv" (i + 1)) accept
            writeVectorToCsv (pack $ printf "flow_%02d.csv" i) flow
            go (i + 1) ts'
  ts' <- go 0 (G.map (\β -> 1 / β) βs)
  -- let ts' = ts

  -- g1 <- mkXoshiro256PlusPlus (seed + 1)
  -- g2 <- mkXoshiro256PlusPlus (seed + 2)

  -- let numMeasure = 10000
  --     stream gen =
  --       Stream.take numMeasure $
  --         monteCarloSampling
  --           {- sweepSize -} (n * n)
  --           {- numberSweeps -} (numSweeps `div` 100)
  --           h
  --           (G.map (\t -> 1 / t) ts')
  --           gen
  --     fold =
  --       (,,,)
  --         <$> pairFolds (energyFold numReplicas) (energyFold numReplicas)
  --         <*> pairFolds (magnetizationFold numReplicas) (magnetizationFold numReplicas)
  --         <*> overlapFold numReplicas
  --         <*> overlapSquaredFold numReplicas
  -- ((es1, es2), (ms1, ms2), qs, qSqrs) <-
  --   streamFoldM fold $
  --     -- ((,) <$> overlapFold numReplicas <*> overlapSquaredFold numReplicas) $
  --     Stream.zip (stream g1) (stream g2)

  -- forM_ ([0 .. numReplicas - 1] :: [Int]) $ \i -> do
  --   let β = realToFrac $ βs ! i
  --       t = 1 / β
  --       MeanVariance meanE1 varE1 = es1 ! i
  --       MeanVariance meanM1 varM1 = ms1 ! i
  --       MeanVariance meanE2 varE2 = es2 ! i
  --       MeanVariance meanM2 varM2 = ms2 ! i
  --       MeanVariance meanQ varQ = qs ! i
  --       MeanVariance meanQSqr varQSqr = qSqrs ! i
  --   putStrLn . mconcat . intersperse "\t" $
  --     fmap show [t, meanE1, varE1, meanE2, varE2, meanM1, varM1, meanM2, varM2, meanQ, varQ, meanQSqr, varQSqr]

  pure ()

main :: IO ()
main = run =<< execParser opts
  where
    opts = info (pMainSettings <**> helper) fullDesc
