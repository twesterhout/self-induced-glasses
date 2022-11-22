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
  let ts = G.fromList [0.84999996, 0.92132324, 0.93324476, 0.9447699, 0.9559653, 0.96667635, 0.9768916, 0.98701656, 0.9973598, 1.0075397, 1.0173488, 1.0275521, 1.0375121, 1.0476379, 1.0583582, 1.0696511, 1.0812087, 1.0928099, 1.1044731, 1.1163883, 1.2]
      βs = G.map (\t -> 1 / t) ts
      -- betasGeometric 21 0.85 1.2
      numReplicas = G.length βs

  g <- mkXoshiro256PlusPlus seed
  h <- randomIsingModel3D n <$> uniformM g
  -- print h

  let k = 5
      go (i :: Int) ts
        | i >= k = do
            writeVectorToCsv (pack $ printf "temperatures_%02d.csv" i) ts
        | otherwise = do
            writeVectorToCsv (pack $ printf "temperatures_%02d.csv" i) ts
            (ObservableState replicas) <-
              Stream.last . Stream.take 10 $
                monteCarloSampling
                  {- sweepSize -} n
                  {- numberSweeps -} ((2 ^ i) * 100000)
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
  go 0 (G.map (\β -> 1 / β) βs)
  pure ()

main :: IO ()
main = run =<< execParser opts
  where
    opts = info (pMainSettings <**> helper) fullDesc
