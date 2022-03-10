{-# LANGUAGE CApiFFI #-}
{-# LANGUAGE OverloadedStrings #-}

module SelfInducedGlasses.Analysis where

import Control.Monad (forM_, unless)
import Control.Monad.ST (runST)
import Data.Bits
import qualified Data.ByteString.Builder as Builder
import qualified Data.HDF5 as H5
import Data.List (intercalate)
import Data.Text (Text, pack, unpack)
import Data.Vector.Generic (Vector, (!))
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Storable as S
import qualified Data.Vector.Storable.Mutable as SM
import Data.Word
import Foreign.C.Types (CFloat (..), CInt (..), CPtrdiff (..))
import Foreign.Ptr (Ptr, castPtr)
import GHC.Stack
import SelfInducedGlasses.Core
import System.IO
import qualified System.IO.Unsafe

loadCouplings :: Text -> IO (DenseMatrix S.Vector ℝ)
loadCouplings filename = H5.withFile filename H5.ReadOnly $ \f ->
  H5.open f "couplings" >>= H5.readDataset

loadStates :: Text -> IO (DenseMatrix S.Vector Word64)
loadStates filename = H5.withFile filename H5.ReadOnly $ \f ->
  H5.open f "states" >>= H5.readDataset

wordOverlap :: Int -> Word64 -> Word64 -> Int
wordOverlap n a b = 64 - 2 * k
  where
    mask = complement 0 `shiftR` (64 - n)
    k = popCount $
      case compare n 64 of
        LT -> (a `xor` b) .&. mask
        EQ -> a `xor` b
        _ -> error "this should never happen"
{-# INLINE wordOverlap #-}

computeOverlap :: Vector v Word64 => Int -> v Word64 -> v Word64 -> ℝ
computeOverlap n a b = fromIntegral (go 0 0) / fromIntegral n
  where
    (!blocks, !rest) = n `divMod` 64
    go !acc !i
      | i < blocks = go (acc + wordOverlap 64 (G.unsafeIndex a i) (G.unsafeIndex b i)) (i + 1)
      | rest /= 0 = acc + wordOverlap rest (G.unsafeIndex a i) (G.unsafeIndex b i)
      | otherwise = acc
{-# INLINE computeOverlap #-}

-- autocorr :: (HasCallStack, Vector v Word64, Vector v ℝ) => Int -> Int -> DenseMatrix v Word64 -> v ℝ
-- autocorr offset numberBits states = G.fromList $ fmap (computeOverlap numberBits s₀) states'
--   where
--     (s₀, states') = case drop offset (denseMatrixRows states) of
--       xs@(x : _) -> (x, xs)
--       _ -> error "too few states"
-- {-# SCC autocorr #-}

magnetizationPerSite :: Configuration -> ℝ
magnetizationPerSite (Configuration v) = G.sum v / fromIntegral (G.length v)
{-# SCC magnetizationPerSite #-}

energyPerSite :: DenseMatrix S.Vector ℝ -> Configuration -> ℝ
energyPerSite couplings@(DenseMatrix n _ _) x = totalEnergy couplings x / fromIntegral n

computeLocalObservables :: DenseMatrix S.Vector ℝ -> DenseMatrix S.Vector Word64 -> Text -> IO ()
computeLocalObservables couplings@(DenseMatrix n _ _) states filename = do
  let configurations = fmap (unpackConfiguration n) (denseMatrixRows states)
      observables = fmap (\ !σ -> (energyPerSite couplings σ, magnetizationPerSite σ)) configurations
      renderRow (e, m) = mconcat $ [Builder.floatDec e, Builder.charUtf8 ',', Builder.floatDec m]
      renderTable rows = mconcat [renderRow r <> Builder.charUtf8 '\n' | r <- rows]
  withFile (unpack filename) WriteMode $ \h ->
    Builder.hPutBuilder h (renderTable observables)

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

foreign import capi unsafe "helpers.h autocorrelation"
  c_autocorrelation :: CPtrdiff -> Ptr CFloat -> Ptr CFloat -> IO CInt

autocorr :: HasCallStack => S.Vector ℝ -> S.Vector ℝ
autocorr signal = System.IO.Unsafe.unsafePerformIO $ do
  let n = G.length signal
  out <- SM.unsafeNew n
  e <- S.unsafeWith signal $ \(inPtr :: Ptr Float) ->
    SM.unsafeWith out $ \(outPtr :: Ptr Float) ->
      c_autocorrelation (fromIntegral n) (castPtr inPtr) (castPtr outPtr)
  unless (e == 0) $ error "c_autocorrelation ran out of memory"
  S.unsafeFreeze out

cumsum :: (G.Vector v a, Num a) => v a -> v a
cumsum xs = G.scanl1' (+) xs

autoWindow :: HasCallStack => ℝ -> S.Vector ℝ -> Int
autoWindow !c τs
  | G.length τs /= 0 = go 0
  | otherwise = error "τs must be non-empty"
  where
    go !i
      | i < G.length τs =
        if fromIntegral i >= c * (G.unsafeIndex τs i)
          then i
          else go (i + 1)
      | otherwise = G.length τs - 1

integratedAutocorrTime :: HasCallStack => S.Vector ℝ -> ℝ
integratedAutocorrTime ρs = G.unsafeIndex τs window
  where
    τs = G.map (\ρ -> 2 * ρ - 1) (cumsum ρs)
    window = autoWindow 5 τs

extractSingleSpinEvolution :: DenseMatrix S.Vector Word64 -> Int -> S.Vector ℝ
extractSingleSpinEvolution (DenseMatrix n numberWords v) i =
  G.generate n (\ !k -> toFloat $ v ! (wordIndex + k * numberWords) `unsafeShiftR` bitIndex)
  where
    (!wordIndex, !bitIndex) = i `divMod` 64
    toFloat 0 = -1
    toFloat _ = 1

autocorrStates :: Int -> DenseMatrix S.Vector Word64 -> S.Vector ℝ
autocorrStates numberBits m@(DenseMatrix n _ _) = runST $ do
  buffer <- G.unsafeThaw (G.replicate n 0)
  let addToBuffer !v = addToBuffer' 0 v
      addToBuffer' !i !v
        | i < n = GM.modify buffer (+ v ! i) i >> addToBuffer' (i + 1) v
        | otherwise = pure ()
  forM_ [0 .. numberBits - 1] $ \i ->
    addToBuffer $ autocorr (extractSingleSpinEvolution m i)
  G.map (/ fromIntegral numberBits) <$> G.unsafeFreeze buffer

-- autocorr :: (HasCallStack, Vector v Word64, Vector v ℝ) => Int -> Int -> DenseMatrix v Word64 -> v ℝ
-- autocorr offset numberBits states = G.fromList $ fmap (computeOverlap numberBits s₀) states'
--   where
--     (s₀, states') = case drop offset (denseMatrixRows states) of
--       xs@(x : _) -> (x, xs)
--       _ -> error "too few states"
-- {-# SCC autocorr #-}

-- data AnalysisOptions = AnalysisOptions
--   {
--
--   }
