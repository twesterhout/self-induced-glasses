{-# LANGUAGE CApiFFI #-}
{-# LANGUAGE OverloadedStrings #-}

module SelfInducedGlasses.Analysis where

-- import Control.DeepSeq (force)
import Control.Monad (forM_, unless)
-- import Control.Monad.ST (runST)
import Control.Scheduler
import Data.Coerce
import Data.Bits
import qualified Data.ByteString.Builder as Builder
import qualified Data.HDF5 as H5
-- import Data.List (intercalate)
import Data.Text (Text, unpack)
-- import qualified Data.Vector as B
import Data.Vector.Generic (Vector)
import qualified Data.Vector.Generic as G
-- import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Storable as S
import qualified Data.Vector.Storable.Mutable as SM
import qualified Data.Vector.Unboxed as U
import Data.Word
import Foreign.C.Types (CFloat (..), CInt (..), CPtrdiff (..))
import Foreign.Ptr (Ptr, castPtr)
-- import GHC.IO.Handle (BufferMode (BlockBuffering), hSetBuffering)
import GHC.Stack
import SelfInducedGlasses.Core
-- import SelfInducedGlasses.Interaction (Lattice (..))
import System.IO
import qualified System.IO.Unsafe

loadCouplings :: Text -> IO (DenseMatrix S.Vector ℝ)
loadCouplings filename = H5.withFile filename H5.ReadOnly $ \f ->
  H5.open f "couplings" >>= H5.readDataset

loadStates :: Text -> IO (DenseMatrix S.Vector Word64)
loadStates filename = H5.withFile filename H5.ReadOnly $ \f ->
  H5.open f "states" >>= H5.readDataset

-- wordOverlap :: Int -> Word64 -> Word64 -> Int
-- wordOverlap n a b = 64 - 2 * k
--   where
--     mask = complement 0 `shiftR` (64 - n)
--     k = popCount $
--       case compare n 64 of
--         LT -> (a `xor` b) .&. mask
--         EQ -> a `xor` b
--         _ -> error "this should never happen"
-- {-# INLINE wordOverlap #-}

-- computeOverlap :: Vector v Word64 => Int -> v Word64 -> v Word64 -> ℝ
-- computeOverlap n a b = fromIntegral (go 0 0) / fromIntegral n
--   where
--     (!blocks, !rest) = n `divMod` 64
--     go !acc !i
--       | i < blocks = go (acc + wordOverlap 64 (G.unsafeIndex a i) (G.unsafeIndex b i)) (i + 1)
--       | rest /= 0 = acc + wordOverlap rest (G.unsafeIndex a i) (G.unsafeIndex b i)
--       | otherwise = acc
-- {-# INLINE computeOverlap #-}

-- autocorr :: (HasCallStack, Vector v Word64, Vector v ℝ) => Int -> Int -> DenseMatrix v Word64 -> v ℝ
-- autocorr offset numberBits states = G.fromList $ fmap (computeOverlap numberBits s₀) states'
--   where
--     (s₀, states') = case drop offset (denseMatrixRows states) of
--       xs@(x : _) -> (x, xs)
--       _ -> error "too few states"
-- {-# SCC autocorr #-}

foreign import capi unsafe "hamming_weight"
  hamming_weight :: CPtrdiff -> Ptr Word64 -> IO CPtrdiff

magnetizationPerSite :: Configuration -> ℝ
magnetizationPerSite (Configuration n v) = m / fromIntegral n
  where
    !m =
      fromIntegral $
        System.IO.Unsafe.unsafePerformIO $
          S.unsafeWith v $ hamming_weight (fromIntegral n)
{-# SCC magnetizationPerSite #-}

energyPerSite :: Couplings -> Configuration -> ℝ
energyPerSite couplings x@(Configuration n _) = totalEnergy couplings x / fromIntegral n

localObservables :: Couplings -> ConfigurationBatch -> U.Vector (ℝ, ℝ)
localObservables couplings states =
  System.IO.Unsafe.unsafePerformIO $ do
    putStrLn "[*] Computing local observables ..."
    !r <- G.concat <$> traverseConcurrently Par process chunks
    putStrLn "[+] Done!"
    pure r
  where
    chunks = batchChunksOf chunkSize states
    compute σ =
      let !e = energyPerSite couplings σ
          !m = magnetizationPerSite σ
       in (e, m)
    -- We go through unboxed vector to force all results to be fully computed
    process = pure . U.fromList . fmap compute . batchToList
    chunkSize = 1000

computeLocalObservables :: Couplings -> ConfigurationBatch -> Text -> IO ()
computeLocalObservables couplings states filename = do
  let observables = localObservables couplings states
      renderRow (e, m) = Builder.floatDec e <> Builder.charUtf8 ',' <> Builder.floatDec m
      renderTable rows = G.foldl' (\ !b !r -> b <> renderRow r <> Builder.charUtf8 '\n') mempty rows
  {-# SCC saving #-} withBinaryFile (unpack filename) WriteMode $ \h -> do
    hSetBuffering h (BlockBuffering Nothing) -- (Just 65536))
    Builder.hPutBuilder h (renderTable observables)

computeAutocorrFunction :: ConfigurationBatch -> Text -> IO ()
computeAutocorrFunction states filename = do
  let table =
        mconcat $
          fmap (\f -> Builder.floatDec f <> Builder.charUtf8 '\n') $
            G.toList $
              autocorrStates states
  withFile (unpack filename) WriteMode $ \h ->
    Builder.hPutBuilder h table

foreign import capi unsafe "two_point_autocorr_function"
  two_point_autocorr_function :: Int -> Int -> Ptr Word64 -> Int -> Ptr CFloat -> IO ()

twoPointAutocorrFunction :: Int -> ConfigurationBatch -> S.Vector ℝ
twoPointAutocorrFunction t_w states@(ConfigurationBatch numberBits (DenseMatrix n _ bits))
  | t_w < n = System.IO.Unsafe.unsafePerformIO $ do
    let size = n - t_w
    buffer <- G.unsafeThaw (G.replicate size 0)
    S.unsafeWith bits $ \(bitsPtr :: Ptr Word64) ->
      SM.unsafeWith buffer $ \(outPtr :: Ptr Float) ->
        two_point_autocorr_function n t_w bitsPtr numberBits (castPtr outPtr)
    G.unsafeFreeze buffer
-- twoPointAutocorrFunction :: Int -> ConfigurationBatch -> S.Vector ℝ
-- twoPointAutocorrFunction t_w states@(ConfigurationBatch numberBits (DenseMatrix n _ _))
--   | t_w < n = System.IO.Unsafe.unsafePerformIO $ do
--     let size = n - t_w
--     buffer <- G.unsafeThaw (G.replicate size 0)
--     let process i =
--           let !σ = G.drop t_w $ extractSingleSpinEvolution states i
--               !σ₀ = coerce (G.unsafeHead σ) :: CFloat
--            in S.unsafeWith σ $ \(xPtr :: Ptr Float) ->
--                 SM.unsafeWith buffer $ \(yPtr :: Ptr Float) ->
--                   contiguous_axpy size (σ₀ / fromIntegral numberBits) (castPtr xPtr) (castPtr yPtr)
--     forM_ [0 .. numberBits - 1] process
--     G.unsafeFreeze buffer

computeTwoPointAutocorrFunction :: Int -> ConfigurationBatch -> Text -> IO ()
computeTwoPointAutocorrFunction t_w states filename = do
  let table =
        mconcat $
          fmap (\f -> Builder.floatDec f <> Builder.charUtf8 '\n') $
            G.toList $
              twoPointAutocorrFunction t_w states
  withFile (unpack filename) WriteMode $ \h ->
    Builder.hPutBuilder h table
-- autocorrStates :: ConfigurationBatch -> S.Vector ℝ
-- autocorrStates states@(ConfigurationBatch numberBits (DenseMatrix n _ _)) =
--   System.IO.Unsafe.unsafePerformIO $ do
--     let process = pure . autocorr . extractSingleSpinEvolution states
--     putStrLn "[*] Running process ... "
--     autocorrs <- traverseConcurrently Seq process [0 .. numberBits - 1]
--     putStrLn "[*] Reducing ..."
--     buffer <- G.unsafeThaw (G.replicate n 0)
--     forM_ autocorrs $ \x ->
--       S.unsafeWith x $ \(xPtr :: Ptr Float) ->
--         SM.unsafeWith buffer $ \(yPtr :: Ptr Float) ->
--           contiguous_axpy n (1 / fromIntegral numberBits) (castPtr xPtr) (castPtr yPtr)
--     G.unsafeFreeze buffer

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

-- saveForGnuplot :: Lattice -> FilePath -> Configuration -> IO ()
-- saveForGnuplot (Lattice (width, height) _) filepath (Configuration v) = do
--   let m = DenseMatrix height width v
--   withFile filepath WriteMode $ \h ->
--     forM_ [0 .. height - 1] $ \i ->
--       hPutStrLn h
--         . intercalate "\t"
--         . fmap show
--         . G.toList
--         . getRow (height - 1 - i)
--         $ m

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

foreign import capi unsafe "helpers.h extract_evolution"
  extract_evolution :: Int -> Int -> Ptr Word64 -> Int -> Ptr CFloat -> IO ()

extractSingleSpinEvolution :: ConfigurationBatch -> Int -> S.Vector ℝ
extractSingleSpinEvolution (ConfigurationBatch _ (DenseMatrix n numberWords v)) i =
  -- System.IO.Unsafe.unsafePerformIO $ do
  --   dest <- SM.unsafeNew n
  --   S.unsafeWith v $ \inputPtr ->
  --     SM.unsafeWith dest $ \outPtr ->
  --       extract_evolution n numberWords inputPtr i (castPtr outPtr)
  --   S.unsafeFreeze dest
  G.generate n (toFloat . extractSpin)
  where
    (!wordIndex, !bitIndex) = i `divMod` 64
    extractSpin !k = (.&. 1) . (`unsafeShiftR` bitIndex) $ G.unsafeIndex v (wordIndex + k * numberWords)
    toFloat !x = 2 * fromIntegral x - 1
{-# SCC extractSingleSpinEvolution #-}

foreign import capi unsafe "helpers.h contiguous_axpy"
  contiguous_axpy :: Int -> CFloat -> Ptr CFloat -> Ptr CFloat -> IO ()

autocorrStates :: ConfigurationBatch -> S.Vector ℝ
autocorrStates states@(ConfigurationBatch numberBits (DenseMatrix n _ _)) =
  System.IO.Unsafe.unsafePerformIO $ do
    let process = pure . autocorr . extractSingleSpinEvolution states
    putStrLn "[*] Running process ... "
    autocorrs <- traverseConcurrently Seq process [0 .. numberBits - 1]
    putStrLn "[*] Reducing ..."
    buffer <- G.unsafeThaw (G.replicate n 0)
    forM_ autocorrs $ \x ->
      S.unsafeWith x $ \(xPtr :: Ptr Float) ->
        SM.unsafeWith buffer $ \(yPtr :: Ptr Float) ->
          contiguous_axpy n (1 / fromIntegral numberBits) (castPtr xPtr) (castPtr yPtr)
    G.unsafeFreeze buffer

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
