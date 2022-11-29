{-# LANGUAGE DuplicateRecordFields #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TypeApplications #-}

module Main (main) where

import Data.ByteString.Builder (Builder)
import qualified Data.ByteString.Builder as Builder
import qualified Data.List
import Data.Strict.Tuple (Pair (..))
import Data.Text (Text, pack, unpack)
import qualified Data.Vector.Fusion.Bundle.Monadic as Bundle
import Data.Vector.Fusion.Bundle.Size (Size (..))
import qualified Data.Vector.Fusion.Stream.Monadic as Stream
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U
import Options.Applicative
import SelfInducedGlasses.Sampling
import System.IO
import Text.Printf (printf)
import qualified UnliftIO.Async as UnliftIO

data MainSettings = MainSettings
  { sideLength :: !Int,
    lambda :: !Double,
    order :: !Int
  }

pMainSettings :: Parser MainSettings
pMainSettings =
  MainSettings
    <$> option auto (long "side-length" <> help "Width and height of the lattice")
    <*> option auto (long "lambda" <> help "Parameter λ of the Hamiltonian")
    <*> option auto (long "order" <> value 10000 <> help "Resummation order")

writeToCsv :: Text -> [(Int, Int, Double)] -> IO ()
writeToCsv filename table =
  withFile (unpack filename) WriteMode $ \h ->
    Builder.hPutBuilder h $
      mconcat . Data.List.intersperse (Builder.charUtf8 '\n') $
        fmap renderRow table
  where
    renderRow (x, y, coupling) =
      mconcat $
        Data.List.intersperse (Builder.charUtf8 ',') $
          [Builder.intDec x, Builder.intDec y, Builder.doubleDec coupling]

run :: MainSettings -> IO ()
run settings = do
  let n = settings.sideLength
      λ = settings.lambda
      points :: [Pair Int Int]
      points = do
        x <- [1 .. n - 1]
        y <- [0 .. x]
        pure (x :!: y)
  couplings <-
    UnliftIO.pooledForConcurrently points $ \(x :!: y) -> do
      c <-
        Stream.last . Stream.take settings.order $
          resummedRkkyInteraction λ (n :!: n) (x :!: y)
      pure (x, y, c)
  writeToCsv (pack $ printf "kolmus_model_n=%d_λ=%.4f.csv" n λ) couplings

  pure ()

main :: IO ()
main = run =<< execParser opts
  where
    opts = info (pMainSettings <**> helper) fullDesc
