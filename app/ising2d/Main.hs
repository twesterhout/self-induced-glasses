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
import System.Random.Stateful (freezeGen, thawGen)
import Text.Printf (printf)

data MainSettings = MainSettings
  { settingsSideLength :: !Int,
    settingsSeed :: !Int,
    -- settingsTemperature :: !ℝ,
    settingsField :: !ℝ
  }

pMainSettings :: Parser MainSettings
pMainSettings =
  MainSettings
    <$> option auto (short 'L' <> help "System side length")
    <*> option auto (long "seed" <> help "Random number generator seed")
    <*> option auto (short 'B' <> help "Magnetic field")

data Acc = Acc
  { accEnergy :: !Double,
    accEnergySquared :: !Double,
    accMagnetization :: !Double,
    accMagnetizationSquared :: !Double,
    accProbability :: !Double
  }

instance Semigroup Acc where
  (<>) a b =
    Acc
      { accEnergy = accEnergy a + accEnergy b,
        accEnergySquared = accEnergySquared a + accEnergySquared b,
        accMagnetization = accMagnetization a + accMagnetization b,
        accMagnetizationSquared = accMagnetizationSquared a + accMagnetizationSquared b,
        accProbability = accProbability a + accProbability b
      }

instance Monoid Acc where
  mempty =
    Acc
      { accEnergy = 0,
        accEnergySquared = 0,
        accMagnetization = 0,
        accMagnetizationSquared = 0,
        accProbability = 0
      }

data SecondMomentAcc a = SecondMomentAcc {smaCount :: !Int, smaMean :: !a, smaM2 :: !a}
  deriving stock (Show, Eq)

zeroAcc :: Num a => SecondMomentAcc a
zeroAcc = SecondMomentAcc 0 0 0

updateAcc :: Fractional a => SecondMomentAcc a -> a -> SecondMomentAcc a
updateAcc (SecondMomentAcc n μ m2) !x = SecondMomentAcc n' μ' m2'
  where
    !n' = n + 1
    !μ' = μ + δx / fromIntegral n'
    !m2' = m2 + δx * (x - μ')
    !δx = x - μ

extractMoments :: Fractional a => SecondMomentAcc a -> (a, a)
extractMoments (SecondMomentAcc n μ m2)
  | n == 0 = (0 / 0, 0 / 0)
  | n == 1 = (μ, 0 / 0)
  | otherwise = (μ, m2 / fromIntegral n)

run :: MainSettings -> IO ()
run settings = do
  let n = settingsSideLength settings
      b = settingsField settings
      βs = betasGeometric 21 0.1 5.0
      -- ts = G.map (\β -> 1 / β) βs
      numReplicas = G.length βs
      numSpins = n * n
      h = ferromagneticIsingModelSquare2D n b

  -- print βs

  g <- mkXoshiro256PlusPlus (settingsSeed settings)

  -- [0.1,0.19870201,0.31110653,0.45704943,0.7657603,1.1811346,1.4555979,1.6799563,1.8841829,1.9659275,2.0476718,2.1294162,2.2111607,2.2945392,2.3972168,2.4998949,2.6025727,2.7052505,2.8919132,3.309102,5.0]
  -- [0.1,0.45570597,1.1840377,1.448896,1.6580452,1.852513,1.9770592,2.0619533,2.1342323,2.2020564,2.2522564,2.2995596,2.3455229,2.3914862,2.4514763,2.5184853,2.603465,2.6941645,2.8551204,3.1773748,5.0]
  --
  --
  -- [0.1,0.17065893,0.33317283,0.4913308,0.9386934,1.4166639,1.6537647,1.8354808,1.9360626,2.0100803,2.084098,2.1581159,2.2321339,2.3155596,2.4250574,2.5345552,2.644053,2.753551,3.0828357,3.5608215,5.0]
  -- [0.1,0.36153185,0.74150276,1.1706158,1.4656792,1.6417818,1.7887982,1.9057987,1.9985542,2.0870912,2.1582294,2.2238133,2.2897656,2.3617556,2.4409685,2.5375166,2.642117,2.7515118,2.9799993,3.265915,5.0]

  let k = 2
      go (i :: Int) ts
        | i >= k = do
            writeVectorToCsv (pack $ printf "temperatures_%02d.csv" i) ts
        | otherwise = do
            writeVectorToCsv (pack $ printf "temperatures_%02d.csv" i) ts
            (ObservableState replicas) <-
              Stream.last . Stream.take 10 $
                monteCarloSampling
                  {- sweepSize -} 400
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
  go 0 (G.map (\β -> 1 / β) βs)
  pure ()

-- forM_ ([0 .. G.length βs - 1] :: [Int]) $ \i -> do
--   let t = ts ! i
--       r = replicas ! i
--       f = (fromIntegral r.hist.numBlue :: Double) / fromIntegral (r.hist.numRed + r.hist.numBlue)
--   printf "%f\t%f\n" t f

{-
g <- freezeGen =<< mkXoshiro256PlusPlus (settingsSeed settings)
(es, ms, spinSpin, perSiteM) <-
  monteCarloSampling'
    {- sweepSize -} 100
    {- numberSweeps -} 50
    {- numberSkip -} 100
    {- numberMeasure -} 10000
    h
    βs
    ((,,,) <$> energyFold <*> magnetizationFold <*> structureFactorFold <*> perSiteMagnetizationFold)
    g

forM_ ([0 .. G.length βs - 1] :: [Int]) $ \i -> do
  let β = realToFrac $ βs ! i
      t = 1 / β
      MeanVariance μE varE = es ! i
      MeanVariance μM varM = ms ! i
      cV = β * β * varE
      χM = β * varM
  putStrLn . mconcat . intersperse "\t" $
    fmap show [t, μE / fromIntegral numSpins, cV, μM / fromIntegral numSpins, χM]

forM_ ([0 .. G.length βs - 1] :: [Int]) $ \i -> do
  let β = βs ! i
      t = 1 / β
  writeMatrixToCsv (pack $ printf "structure_factor_t=%.2f.svg" t) $
    fourierTransformStructureFactorSquare n (spinSpin ! i) (perSiteM ! i)
-}

-- let initial₂ = pure (0, 0)
--     step₂ (eAcc, pAcc) state@(Configuration numSpins@9 _) = do
--       let e' = totalEnergy h state / fromIntegral numSpins
--           -- m' = fromIntegral (totalMagnetization state) / fromIntegral numSpins
--           p' = exp (-β * e')
--       pure (eAcc + e' * p', pAcc + p')
--     extract₂ (eAcc, pAcc) = pure $ eAcc / pAcc
-- μE₂ <- ListT.applyFoldM (Foldl.FoldM step₂ initial₂ extract₂) (allConfigurations (n * n))

-- forM_ ([0 .. G.length βs - 1] :: [Int]) $ \i -> do
--   let β = realToFrac $ βs ! i
--       t = 1 / β
--       μE = μEs ! i
--       cV = β ^ (2 :: Int) * (varEs ! i)
--       μM = μMs ! i
--       χM = β * (varMs ! i)
--   putStrLn $ show t <> "\t" <> show μE <> "\t" <> show cV <> "\t" <> show μM <> "\t" <> show χM

-- putStrLn $
--   show (settingsTemperature settings)
--     <> "\t"
--     <> show μE
--     <> "\t"
--     -- <> show (β ^ (2 :: Int) * varE)
--     -- <> "\t"
--     <> show μM

-- <> "\t"
-- <> show (β * varM)
-- <> "\t"
-- <> show μE₂

-- (energy, heatCapacity, magnetization, susceptibility) <-
--   ListT.applyFoldM (Foldl.FoldM step initial extract) (allConfigurations (n * n))
-- putStrLn $
--   show (settingsTemperature settings)
--     <> "\t"
--     <> show energy
--     <> "\t"
--     <> show heatCapacity
--     <> "\t"
--     <> show magnetization
--     <> "\t"
--     <> show susceptibility

main :: IO ()
main = run =<< execParser opts
  where
    opts = info (pMainSettings <**> helper) fullDesc
