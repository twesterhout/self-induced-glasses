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

data MainSettings = MainSettings {size :: !Int, seed :: !Int}

pMainSettings :: Parser MainSettings
pMainSettings =
  MainSettings
    <$> option auto (long "size" <> help "Side length")
    <*> option auto (long "seed" <> help "Random number generator seed")

run :: MainSettings -> IO ()
run (MainSettings n seed) = do
  let -- ts = G.fromList [0.84999996, 0.95815915, 0.96502686, 0.97176474, 0.978413, 0.9850579, 0.9917415, 0.99841505, 1.0049784, 1.0113064, 1.017554, 1.0241934, 1.030724, 1.0371666, 1.0437058, 1.0504944, 1.0576115, 1.0650607, 1.0727555, 1.0805116, 1.2]
      -- ts = G.fromList [0.6, 0.7249834, 0.8357433, 0.9394234, 1.0331062, 1.1203706, 1.2068709, 1.2935127, 1.3753805, 1.4549646, 1.536332, 1.6221249, 1.7133656, 1.8060983, 1.9043283, 2.0104682, 2.1263938, 2.2658665, 2.43726, 2.6608124, 3.0]
      ts = G.map (\β -> 1 / β) βs
      -- βs = G.map (\t -> 1 / t) ts
      βs = betasGeometric 41 0.2 5
      numReplicas = G.length βs

  g <- mkXoshiro256PlusPlus seed
  h <- randomIsingModel3D n <$> uniformM g
  -- print h

  let k = 5
      go (i :: Int) ts
        | i >= k = do
            writeVectorToCsv (pack $ printf "temperatures_%02d.csv" i) ts
            pure ts
        | otherwise = do
            writeVectorToCsv (pack $ printf "temperatures_%02d.csv" i) ts
            (ObservableState replicas) <-
              Stream.last . Stream.take 10 $
                monteCarloSampling
                  {- sweepSize -} n
                  {- numberSweeps -} ((2 ^ i) * 10000)
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
            writeVectorToCsv (pack $ printf "flow_%02d.csv" (i + 1)) flow
            go (i + 1) ts'
  ts' <- go 0 (G.map (\β -> 1 / β) βs)

  g1 <- mkXoshiro256PlusPlus (seed + 1)
  g2 <- mkXoshiro256PlusPlus (seed + 2)

  let stream gen =
        Stream.take 10000 $
          monteCarloSampling
            {- sweepSize -} n
            {- numberSweeps -} 100
            h
            (G.map (\t -> 1 / t) ts')
            gen
  results <- streamFoldM (overlapSquaredFold numReplicas) $ Stream.zip (stream g1) (stream g2)
  print results

  forM_ ([0 .. numReplicas - 1] :: [Int]) $ \i -> do
    let β = realToFrac $ βs ! i
        t = 1 / β
        MeanVariance meanq2 varq2 = results ! i
    putStrLn . mconcat . intersperse "\t" $
      fmap show [t, meanq2, varq2]

  pure ()

main :: IO ()
main = run =<< execParser opts
  where
    opts = info (pMainSettings <**> helper) fullDesc
