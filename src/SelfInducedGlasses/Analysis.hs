{-# LANGUAGE OverloadedStrings #-}

module SelfInducedGlasses.Analysis where

import Data.Bits
import qualified Data.HDF5 as H5
import Data.Text (Text)
import Data.Vector.Generic (Vector, (!))
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Storable as S
import Data.Word
import GHC.Stack
import SelfInducedGlasses.Core

loadStates :: Text -> IO (DenseMatrix S.Vector Word64)
loadStates filename = H5.withFile filename H5.ReadOnly $ \f ->
  H5.open f "states" >>= H5.readDataset

wordOverlap :: HasCallStack => Int -> Word64 -> Word64 -> Int
wordOverlap n a b = abs $ 2 * k - 64
  where
    mask = complement 0 `shiftR` (64 - n)
    k = popCount $
      case compare n 64 of
        LT -> a .&. b .&. mask
        EQ -> a .&. b
        _ -> error "this should never happen"

computeOverlap :: Vector v Word64 => Int -> v Word64 -> v Word64 -> ℝ
computeOverlap n a b = fromIntegral (go 0 0) / fromIntegral n
  where
    (blocks, rest) = n `divMod` 64
    go !acc !i
      | i < blocks = go (acc + wordOverlap 64 (a ! i) (b ! i)) (i + 1)
      | rest /= 0 = acc + wordOverlap rest (a ! i) (b ! i)
      | otherwise = acc

autocorr :: (Vector v Word64, Vector v ℝ) => Int -> Int -> Int -> DenseMatrix v Word64 -> v ℝ
autocorr offset δt numberBits states = G.fromList $ zipWith (computeOverlap numberBits) a b
  where
    a = drop offset (denseMatrixRows states)
    b = drop (offset + δt) (denseMatrixRows states)

magnetizationPerSite :: Configuration -> ℝ
magnetizationPerSite (Configuration v) = G.sum v / fromIntegral (G.length v)

energyPerSite :: G.Vector v ℝ => DenseMatrix v ℝ -> Configuration -> ℝ
energyPerSite couplings (Configuration v) =
  totalEnergy couplings (G.convert v) / fromIntegral (G.length v)

saveForGnuplot :: Lattice -> FilePath -> Configuration -> IO ()
saveForGnuplot (Lattice (width, height) _) filepath (Configuration v) = do
  let m = DenseMatrix height width v
  withFile filepath WriteMode $ \h ->
    forM_ [0 .. height - 1] $ \i ->
      hPutStrLn h
        . intercalate "\t"
        . fmap show
        . G.toList
        . getRow (height - 1 - i)
        $ m
