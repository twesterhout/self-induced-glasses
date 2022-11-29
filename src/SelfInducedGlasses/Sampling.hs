-- {-# LANGUAGE CApiFFI #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DuplicateRecordFields #-}
{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE TypeFamilies #-}

-- |
-- Module      : SelfInducedGlasses.Sampling
-- Description : Monte Carlo sampling for spin glasses
-- Copyright   : (c) Tom Westerhout, 2022
--
-- Here is a longer description of this module, containing some
-- commentary with @some markup@.
module SelfInducedGlasses.Sampling
  ( -- * General types
    ℝ,
    DenseMatrix (..),
    Hamiltonian (..),
    ferromagneticIsingModelSquare2D,
    randomIsingModel3D,
    sherringtonKirkpatrickModel,
    resummedRkkyInteraction,
    kolmusRkkyModel,
    Configuration (..),
    ConfigurationBatch (..),
    MutableConfiguration (..),
    SweepStats (..),
    ProbDist (..),
    MeanVariance (..),
    ReplicaExchangeSchedule (..),
    ReplicaExchangeState (..),
    ObservableState (..),
    ColorStats (..),
    indexConfiguration,
    totalEnergy,
    replicaOverlap,
    replicaOverlapSquared,
    totalMagnetization,
    energyFold,
    overlapFold,
    overlapSquaredFold,
    magnetizationFold,
    structureFactorFold,
    perSiteMagnetizationFold,
    pairFolds,
    stackFolds,
    streamFoldM,
    autocorrelation,
    autocorrTime,
    fourierTransformStructureFactorSquare,
    monteCarloSampling',
    monteCarloSampling,
    energyChanges,
    maybeExchangeReplicas,
    prettySchedule,
    writeMatrixToCsv,
    writeVectorToCsv,
    measureOne,
    measureTwo,
    singleMeasurementToCsv,
    pairMeasurementToCsv,
    writeMeasurementsToCsv,
    mkIntervals,
    mkSchedule,
    withVector,
    withMutableConfiguration,
    -- doSweep,
    doManySweeps,
    allConfigurations,
    randomConfigurationM,
    updateTemperatures,
    betasGeometric,
  )
where

import Control.DeepSeq
import Control.Foldl (FoldM (..))
import qualified Control.Foldl as Foldl
import Control.Monad (forM_, unless)
import Control.Monad.IO.Unlift
import Control.Monad.Primitive
import Data.Bits
import Data.ByteString.Builder (Builder)
import qualified Data.ByteString.Builder as Builder
import Data.Complex
import qualified Data.List
import Data.Primitive.ByteArray (mutableByteArrayContents)
import Data.Primitive.PVar
import Data.Strict.Tuple (Pair (..))
import qualified Data.Strict.Tuple as Strict
import Data.Text (Text)
import qualified Data.Text as Text
import qualified Data.Text.IO as Text
import Data.Vector (Vector)
import qualified Data.Vector as B
import Data.Vector.Fusion.Stream.Monadic (Step (..), Stream (..))
import qualified Data.Vector.Fusion.Stream.Monadic as Stream
import Data.Vector.Fusion.Util (Id (..))
import Data.Vector.Generic ((!))
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Storable as S
import qualified Data.Vector.Storable.Mutable as SM
import qualified Data.Vector.Unboxed as U
import Data.Word
import Foreign.Ptr
import qualified GHC.Exts as GHC (IsList (..))
import GHC.Generics
import GHC.Stack (HasCallStack)
import SelfInducedGlasses.Random
import System.IO
import System.IO.Unsafe (unsafePerformIO)
import qualified System.Random.MWC.Distributions as MWC
import System.Random.Stateful
import Text.PrettyPrint.ANSI.Leijen (Doc, Pretty (..))
import qualified Text.PrettyPrint.ANSI.Leijen as Pretty
import UnliftIO.Async
import qualified UnliftIO.Foreign as UnliftIO

-- | Real number
type ℝ = Float

-- | Complex number
type ℂ = Complex Float

-- | Dense matrix in row-major order
--
-- We use dense matrices for two purposes:
--
--   * storing the matrix of couplings \(J_{ij}\)
--   * storing bitstrings produced by Monte Carlo (such that each row is a bitstring)
data DenseMatrix v a = DenseMatrix
  { numRows :: {-# UNPACK #-} !Int,
    numCols :: {-# UNPACK #-} !Int,
    elements :: !(v a)
  }
  deriving stock (Show, Eq, Generic)
  deriving anyclass (NFData)

-- | Hamiltonian with Ising-type interaction
data Hamiltonian = IsingLike
  { -- | Matrix of couplings \(J_{ij}\). It is assumed to be symmetric.
    hInteraction :: {-# UNPACK #-} !(DenseMatrix S.Vector ℝ),
    -- | External magnetic field \(B_i\).
    hMagneticField :: {-# UNPACK #-} !(S.Vector ℝ)
  }
  deriving stock (Show, Eq, Generic)
  deriving anyclass (NFData)

-- | Bitstring representing a spin configuration.
data Configuration = Configuration
  { -- | Number of spins
    numSpins :: {-# UNPACK #-} !Int,
    -- | Bit string
    bits :: {-# UNPACK #-} !(S.Vector Word64)
  }
  deriving stock (Show, Eq, Generic)
  deriving anyclass (NFData)

-- | A batch of spin configurations.
data ConfigurationBatch
  = ConfigurationBatch
      {-# UNPACK #-} !Int
      -- ^ Number of spins
      {-# UNPACK #-} !(DenseMatrix S.Vector Word64)
      -- ^ Bitstrings -- each row of the matrix is a spin configuration
  deriving stock (Show, Eq, Generic)
  deriving anyclass (NFData)

-- | Mutable counterpart of 'Configuration'.
data MutableConfiguration = MutableConfiguration
  { -- | Number of spins
    numSpins :: {-# UNPACK #-} !Int,
    -- | Bit string
    bits :: {-# UNPACK #-} !(S.MVector RealWorld Word64)
  }
  deriving stock (Generic)
  deriving anyclass (NFData)

-- | Mutable state that Metropolis-Hastings algorithm keeps around.
data MetropolisState = MetropolisState
  { -- | Current spin configuration
    configuration :: {-# UNPACK #-} !MutableConfiguration,
    -- | Current energy
    energy :: {-# UNPACK #-} !(PVar Double RealWorld),
    -- | Potential energy changes upon spin flips
    deltaEnergies :: {-# UNPACK #-} !(S.MVector RealWorld ℝ)
  }
  deriving stock (Generic)
  deriving anyclass (NFData)

-- | General information about the sampling that is gathered during a sweep.
--
-- 'SweepStats' obtained from multiple sweeps can be combined using the 'Monoid' instance.
data SweepStats = SweepStats {ssTime :: {-# UNPACK #-} !Int, ssAcceptProb :: {-# UNPACK #-} !Double}
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

-- | Probability distribution \(\exp^{-\beta H}\)
data ProbDist = ProbDist {pdBeta :: {-# UNPACK #-} !ℝ, pdHamiltonian :: !Hamiltonian}
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

data ReplicaExchangeSchedule = ReplicaExchangeSchedule
  { resIntervals :: {-# UNPACK #-} !(Vector [Pair Int Int]),
    resSwaps :: {-# UNPACK #-} !(Vector (Pair Int ℝ))
  }
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

data ReplicaColor = ReplicaBlue | ReplicaRed
  deriving stock (Eq, Show, Generic)
  deriving anyclass (NFData)

data ColorStats = ColorStats {numBlue :: {-# UNPACK #-} !Int, numRed :: {-# UNPACK #-} !Int}
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

data ReplicaExchangeState g = ReplicaExchangeState
  { prob :: !ProbDist,
    state :: !MetropolisState,
    color :: !(Maybe ReplicaColor),
    hist :: !ColorStats,
    stats :: !SweepStats,
    gen :: !g
  }
  deriving stock (Generic)
  deriving anyclass (NFData)

newtype ObservableState g = ObservableState {unObservableState :: (Vector (ReplicaExchangeState g))}
  deriving stock (Generic)
  deriving anyclass (NFData)

data SingleMeasurement = SingleMeasurement
  { energy :: !Double,
    magnetization :: !Double,
    flow :: !Double
  }
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

data PairMeasurement = PairMeasurement
  { energy :: !(Pair Double Double),
    magnetization :: !(Pair Double Double),
    flow :: !(Pair Double Double),
    overlap :: !Double
  }
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

updateVarianceAcc :: Fractional a => (Int, a, a) -> a -> (Int, a, a)
updateVarianceAcc (!n, !μ, !m2) !x = (n', μ', m2')
  where
    !n' = n + 1
    !μ' = μ + δx / fromIntegral n'
    !m2' = m2 + δx * (x - μ')
    !δx = x - μ
{-# INLINE updateVarianceAcc #-}

data MeanVariance a = MeanVariance {mean :: !a, variance :: !a}
  deriving stock (Show, Eq, Generic)
  deriving anyclass (NFData)

extractMoments :: Fractional a => (Int, a, a) -> MeanVariance a
extractMoments (!n, !μ, !m2)
  | n == 0 = MeanVariance (0 / 0) (0 / 0)
  | n == 1 = MeanVariance μ (0 / 0)
  | otherwise = MeanVariance μ (m2 / fromIntegral n)

instance Semigroup ColorStats where
  (ColorStats a b) <> (ColorStats c d) = ColorStats (a + c) (b + d)

instance Monoid ColorStats where
  mempty = ColorStats 0 0

updateColorStats :: ReplicaExchangeState g -> ReplicaExchangeState g
updateColorStats g = g {hist = hist'}
  where
    (ColorStats blue red) = g.hist
    hist' = case g.color of
      Just ReplicaBlue -> ColorStats (blue + 1) red
      Just ReplicaRed -> ColorStats blue (red + 1)
      Nothing -> ColorStats blue red
{-# INLINE updateColorStats #-}

updateTemperatures :: Vector ℝ -> Vector ColorStats -> Vector ℝ
updateTemperatures ts colorStats = G.generate (G.length ts) getTemperature
  where
    m = G.length ts
    -- Calculate f(Tᵢ)
    fs = G.map (\c -> fromIntegral c.numBlue / fromIntegral (c.numBlue + c.numRed)) colorStats
    -- Calculate √|f(Tᵢ₊₁) - f(Tᵢ)|
    factors = G.zipWith (\fnext f -> (sqrt . abs) (f - fnext)) (G.drop 1 fs) fs
    -- Inverse normalization constant C'
    c' = G.sum factors
    getTemperature !k
      | k == 0 = G.head ts
      | k == G.length ts - 1 = G.last ts
      | otherwise = go 0 0 (c' * fromIntegral k / fromIntegral (m - 1))
      where
        go !s !i !target
          | s' < target = go s' (i + 1) target
          | otherwise =
              let t = ts ! i
                  δt = ts ! (i + 1) - t
               in t + δt * (target - s) / (factors ! i)
          where
            s' = s + (factors ! i)

betasGeometric ::
  -- | Number temperatures
  Int ->
  -- | Minimal temperature
  Double ->
  -- | Maximal temperature
  Double ->
  -- | βs
  Vector ℝ
betasGeometric n tₘᵢₙ tₘₐₓ = G.generate n (\i -> realToFrac $ 1 / (tₘᵢₙ * (rate ^ i)))
  where
    rate = (tₘₐₓ / tₘᵢₙ) ** ((1 :: Double) / fromIntegral (n - 1))

secondMomentFold ::
  forall m a r.
  (MonadUnliftIO m, Fractional a) =>
  (r -> m a) ->
  FoldM m r (MeanVariance a)
secondMomentFold f = FoldM step begin extract
  where
    begin = pure (0, 0, 0)
    step !acc !x = updateVarianceAcc acc <$> f x
    extract = pure . extractMoments

stackFolds ::
  MonadUnliftIO m =>
  Int ->
  FoldM m a b ->
  FoldM m (Vector a) (Vector b)
stackFolds count (FoldM step begin extract) = FoldM step' begin' extract'
  where
    begin' = withRunInIO $ \u -> GM.replicateM count (u begin)
    step' !accs !xs = withRunInIO $ \u -> do
      GM.iforM_ accs $ \i acc ->
        GM.write accs i =<< u (step acc (xs ! i))
      pure accs
    extract' !accs = G.mapM extract =<< liftIO (G.unsafeFreeze accs)

pairFolds ::
  Monad m =>
  FoldM m a b ->
  FoldM m c d ->
  FoldM m (a, c) (b, d)
pairFolds fold1 fold2 =
  (,) <$> Foldl.premapM (pure . fst) fold1 <*> Foldl.premapM (pure . snd) fold2

energyFold ::
  MonadUnliftIO m =>
  Int ->
  FoldM m (ObservableState g) (Vector (MeanVariance Double))
energyFold numReplicas =
  Foldl.premapM (pure . unObservableState) $
    stackFolds numReplicas $
      secondMomentFold (\r -> liftIO $ readPVar r.state.energy)

magnetizationFold ::
  MonadUnliftIO m =>
  Int ->
  FoldM m (ObservableState g) (Vector (MeanVariance Double))
magnetizationFold numReplicas =
  Foldl.premapM (pure . unObservableState) $
    stackFolds numReplicas $
      secondMomentFold compute
  where
    compute !r =
      (abs . totalMagnetization)
        <$> unsafeFreezeConfiguration r.state.configuration

genericOverlapFold ::
  MonadUnliftIO m =>
  (Configuration -> Configuration -> Double) ->
  Int ->
  FoldM m (ObservableState g, ObservableState g) (Vector (MeanVariance Double))
genericOverlapFold f numReplicas =
  Foldl.premapM preprocess $
    stackFolds numReplicas $
      secondMomentFold compute
  where
    preprocess (a, b) = pure $ G.zip (unObservableState a) (unObservableState b)
    compute !(r1, r2) =
      f
        <$> unsafeFreezeConfiguration r1.state.configuration
        <*> unsafeFreezeConfiguration r2.state.configuration

overlapFold ::
  MonadUnliftIO m =>
  Int ->
  FoldM m (ObservableState g, ObservableState g) (Vector (MeanVariance Double))
overlapFold = genericOverlapFold (\a b -> abs $ replicaOverlap a b)

overlapSquaredFold ::
  MonadUnliftIO m =>
  Int ->
  FoldM m (ObservableState g, ObservableState g) (Vector (MeanVariance Double))
overlapSquaredFold = genericOverlapFold replicaOverlapSquared

perSiteMagnetizationFold ::
  forall m g.
  (MonadUnliftIO m) =>
  FoldM m (ObservableState g) (B.Vector (S.Vector ℝ))
perSiteMagnetizationFold = Foldl.FoldM step begin extract
  where
    begin :: m (Int, Maybe (B.Vector (S.MVector RealWorld ℝ)))
    begin = pure (0, Nothing)
    go ::
      B.Vector (S.MVector RealWorld ℝ) ->
      ObservableState g ->
      Int ->
      m ()
    go !accs !s !i
      | i < G.length accs = do
          x@(Configuration n _) <-
            unsafeFreezeConfiguration (s.unObservableState ! i).state.configuration
          withConfiguration x $ \xPtr ->
            withMutableVector (accs ! i) $ \outPtr ->
              liftIO $ c_addMagnetization n xPtr outPtr
          go accs s (i + 1)
      | otherwise = pure ()
    step (!k, !maybeAccs) !s = do
      accs <- case maybeAccs of
        Just accs -> pure accs
        Nothing ->
          let numReplicas = G.length . unObservableState $ s
              (MutableConfiguration n _) =
                (G.head s.unObservableState).state.configuration
           in G.replicateM numReplicas $ liftIO (GM.new n)
      go accs s 0 >> pure (k + 1, Just accs)
    extract (!k, !maybeAccs) = do
      case maybeAccs of
        Just accs ->
          G.forM accs $ \v ->
            G.map (/ fromIntegral k) <$> liftIO (G.unsafeFreeze v)
        Nothing -> pure G.empty

structureFactorFold ::
  forall m g.
  (MonadUnliftIO m) =>
  FoldM m (ObservableState g) (B.Vector (DenseMatrix S.Vector ℝ))
structureFactorFold = Foldl.FoldM step begin extract
  where
    begin :: m (Int, Maybe (B.Vector (DenseMatrix (S.MVector RealWorld) ℝ)))
    begin = pure (0, Nothing)
    go ::
      B.Vector (DenseMatrix (S.MVector RealWorld) ℝ) ->
      ObservableState g ->
      Int ->
      m ()
    go !accs !s !i
      | i < G.length accs = do
          x <- unsafeFreezeConfiguration (s.unObservableState ! i).state.configuration
          let (DenseMatrix n _ out) = accs ! i
          withConfiguration x $ \xPtr ->
            withMutableVector out $ \outPtr ->
              liftIO $ c_addSpinSpinCorrelation n xPtr outPtr
          go accs s (i + 1)
      | otherwise = pure ()
    step (!k, !maybeAccs) !s = do
      accs <- case maybeAccs of
        Just accs -> pure accs
        Nothing ->
          let numReplicas = G.length . unObservableState $ s
              (MutableConfiguration n _) =
                (G.head s.unObservableState).state.configuration
           in G.replicateM numReplicas $
                DenseMatrix n n <$> liftIO (GM.new (n * n))
      go accs s 0 >> pure (k + 1, Just accs)
    extract (!k, !maybeAccs) = do
      case maybeAccs of
        Just accs ->
          G.forM accs $ \(DenseMatrix nRows nCols v) ->
            DenseMatrix nRows nCols . G.map (/ fromIntegral k)
              <$> liftIO (G.unsafeFreeze v)
        Nothing -> pure G.empty

measureOne :: MonadUnliftIO m => ReplicaExchangeState g -> m SingleMeasurement
measureOne r = do
  e <- liftIO $ readPVar r.state.energy
  m <- totalMagnetization <$> unsafeFreezeConfiguration r.state.configuration
  let c = r.hist
      f = fromIntegral c.numBlue / fromIntegral (c.numBlue + c.numRed)
  pure $ SingleMeasurement e m f
{-# INLINE measureOne #-}

measureTwo :: MonadUnliftIO m => ReplicaExchangeState g -> ReplicaExchangeState g -> m PairMeasurement
measureTwo r1 r2 = do
  (SingleMeasurement e1 m1 f1) <- measureOne r1
  (SingleMeasurement e2 m2 f2) <- measureOne r2
  q <-
    replicaOverlap
      <$> unsafeFreezeConfiguration r1.state.configuration
      <*> unsafeFreezeConfiguration r2.state.configuration
  pure $ PairMeasurement (e1 :!: e2) (m1 :!: m2) (f1 :!: f2) q
{-# INLINE measureTwo #-}

fourierTransformStructureFactorSquare ::
  Int ->
  DenseMatrix S.Vector ℝ ->
  S.Vector ℝ ->
  DenseMatrix S.Vector ℝ
fourierTransformStructureFactorSquare numPoints spinSpin spin =
  generateDenseMatrix numPoints numPoints $ \ !i !j ->
    let !qx = prefactor * fromIntegral j
        !qy = prefactor * fromIntegral (numPoints - 1 - i)
     in computeOne qx qy
  where
    !prefactor = 2 * pi / fromIntegral numPoints
    !n = G.length spin
    !sideLength = round $ sqrt (fromIntegral n :: Double)
    !xs = U.generate n (`mod` sideLength)
    !ys = U.generate n (`div` sideLength)
    term :: ℝ -> ℝ -> Int -> Int -> ℝ
    term !qx !qy !i !j = (1 / fromIntegral n) * cos (qx * δx + qy * δy) * (sᵢsⱼ - sᵢ * sⱼ)
      where
        !δx = fromIntegral $ G.unsafeIndex xs i - G.unsafeIndex xs j
        !δy = fromIntegral $ G.unsafeIndex ys i - G.unsafeIndex ys j
        !sᵢsⱼ = indexDenseMatrix spinSpin i j
        !sᵢ = G.unsafeIndex spin i
        !sⱼ = G.unsafeIndex spin j
    computeOne !qx !qy =
      iFold 0 (< n) (+ 1) 0 $ \ !accᵢ !i ->
        (accᵢ +) $
          iFold 0 (< n) (+ 1) 0 $ \ !accⱼ !j ->
            accⱼ + term qx qy i j

loopM :: Monad m => i -> (i -> Bool) -> (i -> i) -> (i -> m ()) -> m ()
loopM i₀ cond inc action = go i₀
  where
    go !i
      | cond i = do () <- action i; go (inc i)
      | otherwise = pure ()
{-# INLINE loopM #-}

{-
iFoldM :: Monad m => i -> (i -> Bool) -> (i -> i) -> a -> (a -> i -> m a) -> m a
iFoldM i₀ cond inc x₀ action = go x₀ i₀
  where
    go !x !i
      | cond i = do !x' <- action x i; go x' (inc i)
      | otherwise = pure x
{-# INLINE iFoldM #-}
-}

iFold :: i -> (i -> Bool) -> (i -> i) -> a -> (a -> i -> a) -> a
iFold i₀ cond inc x₀ action = go x₀ i₀
  where
    go !x !i
      | cond i = let !x' = action x i in go x' (inc i)
      | otherwise = x
{-# INLINE iFold #-}

-- | Get matrix shape
-- dmShape :: DenseMatrix v a -> (Int, Int)
-- dmShape m = (m.numRows, m.numCols)
denseMatrixFromList :: G.Vector v a => [[a]] -> DenseMatrix v a
denseMatrixFromList rs
  | G.length elts == nRows * nCols = DenseMatrix nRows nCols elts
  | otherwise = error "nested list has irregular shape"
  where
    !nRows = length rs
    !nCols = case rs of
      [] -> 0
      (r : _) -> length r
    elts = G.fromListN (nRows * nCols) $ mconcat rs

denseMatrixToList :: G.Vector v a => DenseMatrix v a -> [[a]]
denseMatrixToList (DenseMatrix _ nCols v) = go (G.toList v)
  where
    go elts = case splitAt nCols elts of
      (row, []) -> [row]
      (row, rest) -> row : go rest

instance (G.Vector v a) => GHC.IsList (DenseMatrix v a) where
  type Item (DenseMatrix v a) = [a]
  fromList = denseMatrixFromList
  toList = denseMatrixToList

indexDenseMatrix :: G.Vector v a => DenseMatrix v a -> Int -> Int -> a
indexDenseMatrix (DenseMatrix nRows nCols v) i j
  | 0 <= i && i < nRows && 0 <= j && j < nCols = G.unsafeIndex v (i * nCols + j)
  | otherwise = error $ "invalid index: " <> show (i, j)

generateDenseMatrix :: G.Vector v a => Int -> Int -> (Int -> Int -> a) -> DenseMatrix v a
generateDenseMatrix nRows nCols f = runST $ do
  v <- GM.unsafeNew (nRows * nCols)
  loopM 0 (< nRows) (+ 1) $ \ !i ->
    loopM 0 (< nCols) (+ 1) $ \ !j ->
      GM.unsafeWrite v (i * nRows + j) $ f i j
  DenseMatrix nRows nCols <$> G.unsafeFreeze v
{-# INLINE generateDenseMatrix #-}

writeMatrixToCsv :: G.Vector v ℝ => Text -> DenseMatrix v ℝ -> IO ()
writeMatrixToCsv filename matrix@(DenseMatrix nRows nCols _) =
  withFile (Text.unpack filename) WriteMode $ \h ->
    Builder.hPutBuilder h renderMatrix
  where
    renderElement !i !j = Builder.formatFloat Builder.scientific (indexDenseMatrix matrix i j)
    renderRow !i =
      mconcat $
        Data.List.intersperse (Builder.charUtf8 ',') $
          [renderElement i j | j <- [0 .. nCols - 1]]
    renderMatrix =
      mconcat $
        Data.List.intersperse (Builder.charUtf8 '\n') $
          [renderRow i | i <- [0 .. nRows - 1]]

renderFloat :: Float -> Builder
renderFloat = Builder.formatFloat Builder.scientific

writeVectorToCsv :: G.Vector v ℝ => Text -> v ℝ -> IO ()
writeVectorToCsv filename v =
  withFile (Text.unpack filename) WriteMode $ \h ->
    Builder.hPutBuilder h renderVector
  where
    renderElement !i = Builder.formatFloat Builder.scientific (v ! i)
    renderVector =
      mconcat $
        Data.List.intersperse (Builder.charUtf8 '\n') $
          [renderElement i | i <- [0 .. G.length v - 1]]

singleMeasurementToCsv :: SingleMeasurement -> Builder
singleMeasurementToCsv (SingleMeasurement e m f) =
  mconcat $
    Data.List.intersperse (Builder.charUtf8 ',') $
      (renderFloat . realToFrac) <$> [e, m, f]

pairMeasurementToCsv :: PairMeasurement -> Builder
pairMeasurementToCsv (PairMeasurement (e1 :!: e2) (m1 :!: m2) (f1 :!: f2) q) =
  mconcat $
    Data.List.intersperse (Builder.charUtf8 ',') $
      (renderFloat . realToFrac) <$> [e1, e2, m1, m2, f1, f2, q]

writeMeasurementsToCsv ::
  MonadUnliftIO m =>
  (a -> Builder) ->
  Vector Text ->
  Stream m (Vector a) ->
  Stream m (Vector a)
writeMeasurementsToCsv convert paths stream = Stream.mapM go stream
  where
    writeOne filename x = liftIO $ do
      withFile (Text.unpack filename) AppendMode $ \h ->
        Builder.hPutBuilder h $ convert x <> Builder.charUtf8 '\n'
    go xs = G.zipWithM_ writeOne paths xs >> pure xs

indexConfiguration :: ConfigurationBatch -> Int -> Configuration
indexConfiguration (ConfigurationBatch n (DenseMatrix _ numWords v)) i =
  Configuration n (G.slice (i * numWords) numWords v)

maybeUpdateColor ::
  Int ->
  Int ->
  ReplicaExchangeState g ->
  ReplicaExchangeState g
maybeUpdateColor n i r@(ReplicaExchangeState _ _ c _ _ _)
  | i == 0 && (c == Nothing || c == Just ReplicaRed) = r {color = Just ReplicaBlue}
  | i == n - 1 && (c == Nothing || c == Just ReplicaBlue) = r {color = Just ReplicaRed}
  | otherwise = r

-- | Try exchanging two replicas
maybeExchangeReplicas ::
  MonadUnliftIO m =>
  -- | First replica
  ReplicaExchangeState g ->
  -- | Second replica
  ReplicaExchangeState g ->
  -- | A uniform random number in @[0, 1)@. Required to make the decision stochastic.
  ℝ ->
  -- | Potentially updated replicas
  m (ReplicaExchangeState g, ReplicaExchangeState g, Bool)
maybeExchangeReplicas
  !r1@(ReplicaExchangeState (ProbDist β1 _) s1 c1 _ _ _)
  !r2@(ReplicaExchangeState (ProbDist β2 _) s2 c2 _ _ _)
  !u = do
    !δe <-
      liftIO . fmap realToFrac $
        (-) <$> readPVar s2.energy <*> readPVar s1.energy
    let !δβ = β2 - β1
    if u <= exp (δβ * δe)
      then
        let -- We swap the spin configurations only, leaving the βs and random number generators
            -- unchanged.
            !r1' = r1 {state = s2, color = c2}
            !r2' = r2 {state = s1, color = c1}
         in pure (r1', r2', True)
      else pure (r1, r2, False)
{-# INLINE maybeExchangeReplicas #-}

runReplicaExchangeSchedule ::
  forall g m.
  (MonadUnliftIO m, PrimMonad m, g ~ Xoshiro256PlusPlus (PrimState m)) =>
  Int ->
  ReplicaExchangeSchedule ->
  B.Vector (ReplicaExchangeState g) ->
  m (B.Vector (ReplicaExchangeState g))
runReplicaExchangeSchedule sweepSize schedule₀ states₀ = do
  -- liftIO $ do
  --   putStrLn "Running schedule ..."
  --   print schedule₀
  --   Pretty.putDoc $ pretty schedule₀
  --   putStrLn ""
  -- liftIO $ print (resSwaps schedule₀)
  -- liftIO $ do
  --   Pretty.putDoc $ pretty schedule₀
  --   putStrLn ""
  Foldl.foldM (FoldM step initial extract) (resSwaps schedule₀)
  where
    -- spawn ::
    --   ReplicaExchangeState g ->
    --   [(Int, Int)] ->
    --   m (Async (ReplicaExchangeState g), [(Int, Int)])
    spawn !replica@(ReplicaExchangeState probDist _ _ _ _ g) !intervals =
      case intervals of
        ((start :!: end) : others) -> do
          future <- async $ do
            stats' <- doManySweeps probDist (end - start) sweepSize replica.state g
            -- liftIO . putStrLn $ "Result of (" <> show probDist.pdBeta <> "): " <> show stats'
            pure $ replica {stats = replica.stats <> stats'}
          pure (future, others)
        [] -> do
          !future <- async $ pure replica
          pure (future, [])
    -- initial :: m (B.Vector (Async (ReplicaExchangeState g), [(Int, Int)]))
    initial = G.zipWithM spawn states₀ (resIntervals schedule₀)
    -- step ::
    --   B.Vector (Async (ReplicaExchangeState g), [(Int, Int)]) ->
    --   (Int, ℝ) ->
    --   m (B.Vector (Async (ReplicaExchangeState g), [(Int, Int)]))
    step accumulators (i :!: u)
      | G.length accumulators > 1 = do
          let n = G.length accumulators
              (future1, intervals1) = accumulators ! i
              (future2, intervals2) = accumulators ! (i + 1)
          (state1, state2) <- waitBoth future1 future2
          -- TODO: how to update the colors
          (state1', state2', _) <- maybeExchangeReplicas state1 state2 u
          acc1' <- flip spawn intervals1 $ updateColorStats $ maybeUpdateColor n i state1'
          acc2' <- flip spawn intervals2 $ updateColorStats $ maybeUpdateColor n (i + 1) state2'
          pure $ accumulators G.// [(i, acc1'), (i + 1, acc2')]
      | otherwise = do
          let (future1, intervals1) = accumulators ! i
          state1 <- wait future1
          acc1' <- spawn state1 intervals1
          pure $ accumulators G.// [(i, acc1')]
    -- extract ::
    --   B.Vector (Async (ReplicaExchangeState g), [(Int, Int)]) ->
    --   m (B.Vector (ReplicaExchangeState g))
    extract accumulators =
      G.forM accumulators $ \(future, intervals) -> do
        unless (null intervals) $ error "not all intervals have been processed"
        wait future
{-# SCC runReplicaExchangeSchedule #-}

{-
initSweepTaskInputOutputs ::
  MonadUnliftIO m =>
  Vector (ReplicaExchangeState g) ->
  m (Vector (MVar (ReplicaExchangeState g)))
initSweepTaskInputOutputs rs = G.mapM UnliftIO.newMVar rs

scheduleToSweepTasks :: ReplicaExchangeSchedule -> [SweepTask]
scheduleToSweepTasks schedule =
  Data.List.sort $
    G.toList (sweepTasks G.++ swapTasks)
  where
    mkTasksOneReplica i intervals =
      G.fromList $
        fmap (\(start, stop) -> SweepTask i start stop) intervals
    sweepTasks = G.concat . G.toList $ G.imap mkTasksOneReplica schedule.resIntervals
    swapTasks = G.imap (\t (i :!: u) -> SwapTask (i :!: i + 1) (t + 1) u) schedule.resSwaps

instance Ord SweepTask where
  compare (SweepTask _ start1 stop1) (SweepTask _ start2 stop2) =
    compare (start1, stop1) (start2, stop2)
  compare (SwapTask _ t1 _) (SwapTask _ t2 _) = compare t1 t2
  compare (SweepTask _ start1 _) (SwapTask _ t2 _)
    | start1 < t2 = LT
    | otherwise = GT
  compare x@(SwapTask _ _ _) y@(SweepTask _ _ _) =
    case compare y x of
      LT -> GT
      EQ -> EQ
      GT -> LT

runSweepTasks ::
  forall g m.
  (MonadUnliftIO m, PrimMonad m, g ~ Xoshiro256PlusPlus (PrimState m)) =>
  Int ->
  Vector (ReplicaExchangeState g) ->
  [SweepTask] ->
  m (Vector (ReplicaExchangeState g))
runSweepTasks sweepSize replicas₀ tasks = do
  -- liftIO . putStrLn $ "==== runSweepTasks ===="
  -- liftIO . print $ tasks
  -- inputOutputs <- initSweepTaskInputOutputs replicas₀
  let n = G.length replicas₀
  -- capabilityRef <- UnliftIO.newEmptyMVar
  -- currentNumThreads <- UnliftIO.newMVar 0
  maxNumThreads <- min n <$> UnliftIO.getNumCapabilities
  -- numThreads <- newPVar (0 :: Int)
  -- nextTask <- UnliftIO.newEmptyMVar

  let -- spawnNextThread capability task@(SweepTask i start stop) = do
      --   replica <- UnliftIO.takeMVar (inputOutputs ! i)
      --   -- (UnliftIO.mkWeakThreadId =<<) $
      --   -- liftIO . putStrLn $ "Scheduling (" <> show task <> ") on " <> show capability
      --   -- UnliftIO.forkOn capability $ do
      --   UnliftIO.forkIO $ do
      --     !sweepStats <-
      --       doManySweeps
      --         replica.prob
      --         (stop - start)
      --         sweepSize
      --         replica.state
      --         replica.gen
      --     -- liftIO . putStrLn $ "Result of (" <> show task <> ", " <> show replica.prob.pdBeta <> "): " <> show sweepStats
      --     UnliftIO.putMVar (inputOutputs ! i) $
      --       replica {stats = replica.stats <> sweepStats}
      --     UnliftIO.putMVar capabilityRef capability

      -- go ::
      --   Vector (Either (ReplicaExchangeState g) (Async (ReplicaExchangeState g))) ->
      --   [SweepTask] ->
      --   m (Vector (ReplicaExchangeState g))
      go !rs [] = do
        -- liftIO . putStrLn $ "====== waiting for all tasks to complete ====="
        -- forM_ [0 .. min maxNumThreads numThreads - 1] $ \_ -> do
        --   liftIO . putStrLn $ "* takeMVar"
        --   UnliftIO.takeMVar nextTask
        rs' <- G.freeze rs
        G.forM rs' $ \x -> case x of
          Left r -> pure r
          Right f -> wait f
      -- forM_ [0 .. (min numThreads maxNumThreads) - 1] $ \_ -> do
      --   _ <- UnliftIO.takeMVar capabilityRef
      --   pure ()
      -- forM_ tids UnliftIO.killThread
      -- forM_ tids $ \weakThreadId -> do
      --   liftIO (deRefWeak weakThreadId) >>= maybe (pure ()) UnliftIO.killThread
      go !rs (task@(SwapTask (i :!: j) _ u) : otherTasks) = do
        -- liftIO . putStrLn $ "SwapTask: waiting for " <> show i <> " and " <> show j
        r1 <-
          GM.read rs i >>= \x -> case x of
            Right f -> wait f
            Left r -> error $ "in " <> show task
        r2 <-
          GM.read rs i >>= \x -> case x of
            Right f -> wait f
            Left r -> error $ "in " <> show task

        -- liftIO . putStrLn $ "SwapTask: running the swap of " <> show i <> " and " <> show j
        (r1', r2', _) <- maybeExchangeReplicas r1 r2 u
        let !r1'' = updateColorStats $ maybeUpdateColor n i r1'
            !r2'' = updateColorStats $ maybeUpdateColor n j r2'

        -- liftIO . putStrLn $ "SwapTask: donw with swap of " <> show i <> " and " <> show j
        GM.write rs i (Left r1'')
        GM.write rs j (Left r2'')
        go rs otherTasks
      -- go (rs G.// [(i, Left r1''), (j, Left r2'')]) otherTasks
      go !rs (task@(SweepTask i start stop) : otherTasks) = do
        -- let waitForTasks !t = do
        --       k <- atomicReadIntPVar numThreads
        --       if k >= maxNumThreads
        --         then UnliftIO.threadDelay (if t < 10 then 100 else 1000) >> waitForTasks (t + 1)
        --         else atomicAddIntPVar numThreads 1 >> pure ()
        -- {-# SCC waitForTasks #-} waitForTasks 0
        -- oldNumThreads <- atomicAddIntPVar numThreads 1
        -- when (oldNumThreads == maxNumThreads + 1) $ UnliftIO.takeMVar nextTask
        -- numThreads <- UnliftIO.takeMVar currentNumThreads
        -- if numThreads >= maxNumThreads
        --   then do
        --     UnliftIO.putMVar currentNumThreads numThreads
        --     UnliftIO.threadDelay 1000
        --     go rs (task : otherTasks)
        --   else do
        -- UnliftIO.putMVar currentNumThreads (numThreads + 1)
        -- liftIO . putStrLn $ "SweepTask checking futures[" <> show i <> "] for duplicates"
        x <- GM.read rs i
        let (Left replica) = x
        -- liftIO . putStrLn $ "SweepTask spawning a task on " <> show i
        !f <- asyncOn (i `mod` maxNumThreads) $ do
          -- liftIO . putStrLn $ "doManySweeps on " <> show i
          sweepStats <-
            doManySweeps
              replica.prob
              (stop - start)
              sweepSize
              replica.state
              replica.gen
          -- liftIO . putStrLn $ "doManySweeps done on " <> show i
          let !r' = replica {stats = replica.stats <> sweepStats}
          -- numThreads <- UnliftIO.takeMVar currentNumThreads
          -- UnliftIO.putMVar currentNumThreads (numThreads - 1)
          -- liftIO . putStrLn $ "putMVar"
          -- UnliftIO.putMVar nextTask ()
          -- oldNumThreads' <- atomicSubIntPVar numThreads 1
          -- when (oldNumThreads' <= maxNumThreads) $ do
          --   _ <- UnliftIO.tryPutMVar nextTask ()
          --   pure ()
          pure r'
        -- liftIO . putStrLn $ "SweepTask is done on " <> show i
        GM.write rs i (Right f)
        -- go (rs G.// [(i, Right f)]) otherTasks
        -- atomicAddIntPVar numThreads 1
        go rs otherTasks

  rs <- G.thaw $ G.map Left replicas₀
  go rs tasks
-- liftIO . putStrLn $ "====== waiting for all tasks to complete ====="
-- G.generateM n $ \i -> do
--   maybeReplica <- UnliftIO.tryTakeMVar (replicas ! i)
--   case maybeReplica of
--     Just r -> pure r
--     Nothing -> do
--       UnliftIO.takeMVar (futures ! i)
-- G.mapM UnliftIO.takeMVar inputOutputs
-- liftIO . putStrLn $ "returning replicas"
-- G.freeze replicas
{-# SCC runSweepTasks #-}

runReplicaExchangeSchedule' ::
  (MonadUnliftIO m, PrimMonad m, g ~ Xoshiro256PlusPlus (PrimState m)) =>
  Int ->
  ReplicaExchangeSchedule ->
  Vector (ReplicaExchangeState g) ->
  m (Vector (ReplicaExchangeState g))
runReplicaExchangeSchedule' sweepSize schedule states = do
  -- liftIO $ print schedule
  -- liftIO $ do
  --   putStrLn "Running schedule ..."
  --   print schedule
  --   Pretty.putDoc $ pretty schedule
  --   putStrLn ""
  runSweepTasks sweepSize states (scheduleToSweepTasks schedule)
-}

-- data SweepTask
--   = SweepTask
--       { replicaIdx :: {-# UNPACK #-} !Int,
--         startTime :: {-# UNPACK #-} !Int,
--         stopTime :: {-# UNPACK #-} !Int
--       }
--   | SwapTask
--       { replicaIdices :: {-# UNPACK #-} !(Pair Int Int),
--         timePoint :: {-# UNPACK #-} !Int
--       }
-- data ReplicaExchangeSchedule = ReplicaExchangeSchedule
--   { resIntervals :: !(B.Vector [(Int, Int)]),
--     resSwaps :: !(B.Vector (Int, ℝ))
--   }
--   deriving stock (Show, Generic)
--   deriving anyclass (NFData)

-- 🡼 🡽 🡾 🡿
-- "╸╺"
--
-- rowIdx: 0   ━━🡾🡽━━╸╺━━╸╺━━🡾🡽  [(0, 1), (1, 4)]
--         1   ━━🡽🡾━━🡾🡽━━🡾🡽━━🡽🡾  [(0, 1), (1, 2), (2, 3), (3, 4)]
--         2   ━━╸╺━━🡽🡾━━🡽🡾━━    [(0, 2), (2, 3), (3, 4)]
--              [0,  1,  1,  0]

prettyInterval :: G.Vector v (Pair Int ℝ) => Int -> Pair Int Int -> v (Pair Int ℝ) -> Doc
prettyInterval rowIdx (start :!: end) swaps =
  Pretty.text left <> Pretty.hcat middle <> Pretty.text right
  where
    middle = Pretty.punctuate (Pretty.text "╸╺") (replicate (end - start) (Pretty.text "━━"))
    left
      | start <= 0 = "╺"
      | i == rowIdx = "🡽"
      | i + 1 == rowIdx = "🡾"
      | otherwise = "╺"
      where
        ~(i :!: _) = swaps ! (start - 1)
    right
      | end <= 0 = "╸"
      | i == rowIdx = "🡾"
      | i + 1 == rowIdx = "🡽"
      | otherwise = "╸"
      where
        ~(i :!: _) = swaps ! (end - 1)

prettyRow :: ReplicaExchangeSchedule -> Int -> Doc
prettyRow schedule rowIdx =
  mconcat $ fmap (\i -> prettyInterval rowIdx i swaps) intervals
  where
    intervals = (resIntervals schedule) ! rowIdx
    swaps = resSwaps schedule

prettySchedule :: ReplicaExchangeSchedule -> Doc
prettySchedule schedule =
  Pretty.vcat $
    fmap (prettyRow schedule) [0 .. G.length (resIntervals schedule) - 1]

instance Pretty ReplicaExchangeSchedule where
  pretty = prettySchedule

{-
mkIntervals :: G.Vector v (Pair Int ℝ) => v (Pair Int ℝ) -> Int -> [Pair Int Int]
mkIntervals swaps replicaIdx = go [0 :!: G.length swaps] (G.length swaps - 1)
  where
    go acc@((start :!: end) : rest) !i
      | i <= 0 = acc
      | Strict.fst (swaps ! (i - 1)) == replicaIdx
          || Strict.fst (swaps ! (i - 1)) == replicaIdx - 1 =
          go ((start :!: i) : (i :!: end) : rest) (i - 1)
      | otherwise = go acc (i - 1)
    go _ _ = error "this should never happen by construction"
-}

mkIntervals :: G.Vector v (Pair Int ℝ) => v (Pair Int ℝ) -> Int -> [Pair Int Int]
mkIntervals swaps replicaIdx = go 0 relevantSwaps
  where
    n = G.length swaps
    go !s (t : others) = (s :!: t) : go t others
    go !s [] = if s < n then [(s :!: n)] else []
    relevantSwaps =
      fmap Strict.fst $
        filter (\(_ :!: (i :!: _)) -> i == replicaIdx || i + 1 == replicaIdx) $
          zipWith (:!:) [1 ..] (G.toList swaps)

mkSchedule :: G.Vector v (Pair Int ℝ) => Int -> v (Pair Int ℝ) -> ReplicaExchangeSchedule
mkSchedule numReplicas swaps =
  ReplicaExchangeSchedule
    (G.generate numReplicas (mkIntervals swaps))
    (G.convert swaps)

randomReplicaExchangeScheduleM :: StatefulGen g m => Int -> Int -> g -> m ReplicaExchangeSchedule
randomReplicaExchangeScheduleM numReplicas numSwaps g = do
  swaps <-
    B.replicateM numSwaps $
      (:!:)
        <$> uniformRM (0, max 0 (numReplicas - 2)) g
        <*> uniformRM (0, 1) g
  pure $ mkSchedule numReplicas swaps

-- | Generate a random spin configuration of given length.
randomConfigurationM ::
  StatefulGen g m =>
  -- | Desired length
  Int ->
  -- | Random number generator
  g ->
  m Configuration
randomConfigurationM n g
  | n < 0 = error "negative length"
  | otherwise = Configuration n <$> G.generateM numWords mkWord
  where
    numWords = (n + 63) `div` 64
    rest = n `mod` 64
    mkWord !i
      | rest /= 0 && i == numWords - 1 =
          let mask = ((1 :: Word64) `unsafeShiftL` rest) - 1
           in (mask .&.) <$> uniformM g
      | i < numWords = uniformM g
      | otherwise = error $ "should never happen; n=" <> show n <> ", i=" <> show i

-- | Generate a random 'MetropolisState'.
randomMetropolisStateM ::
  (MonadUnliftIO m, StatefulGen g m) =>
  -- | The Hamiltonian
  Hamiltonian ->
  -- | Random number generator
  g ->
  m MetropolisState
randomMetropolisStateM hamiltonian g = do
  spins <- randomConfigurationM n g
  MetropolisState <$> thawConfiguration spins <*> mkCurrentEnergy spins <*> mkDeltaEnergies spins
  where
    n = G.length . hMagneticField $ hamiltonian
    mkCurrentEnergy spins = liftIO . newPVar $ totalEnergy hamiltonian spins
    mkDeltaEnergies spins = liftIO . G.thaw $ energyChanges hamiltonian spins

randomReplicaExchangeStateM ::
  (MonadUnliftIO m, StatefulGen g m) =>
  -- | The Hamiltonian
  Hamiltonian ->
  -- | Inverse temperature β
  ℝ ->
  -- | Random number generator
  g ->
  m (ReplicaExchangeState g)
randomReplicaExchangeStateM hamiltonian β g = do
  s <- randomMetropolisStateM hamiltonian g
  pure $ ReplicaExchangeState (ProbDist β hamiltonian) s Nothing mempty mempty g

-- | Compute energy of a spin configuration
totalEnergy :: Hamiltonian -> Configuration -> Double
totalEnergy hamiltonian s@(Configuration n _) =
  unsafePerformIO $!
    withDenseMatrix (hInteraction hamiltonian) $
      \interactionPtr ->
        withVector (hMagneticField hamiltonian) $ \fieldPtr ->
          withConfiguration s $ \statePtr ->
            pure $ realToFrac (c_totalEnergy n interactionPtr fieldPtr statePtr)

-- | Compute the total magnetization of a spin configuration
totalMagnetization :: Configuration -> Double
totalMagnetization s@(Configuration n _) =
  unsafePerformIO $!
    withConfiguration s $ \statePtr ->
      pure . (/ fromIntegral n) . fromIntegral $ c_totalMagnetization n statePtr

replicaOverlap :: HasCallStack => Configuration -> Configuration -> Double
replicaOverlap state1@(Configuration n _) state2@(Configuration n' _)
  | n == n' =
      unsafePerformIO $!
        withConfiguration state1 $ \statePtr1 ->
          withConfiguration state2 $ \statePtr2 ->
            pure $ (/ fromIntegral n) . realToFrac $ c_replicaOverlap n statePtr1 statePtr2
  | otherwise = error $ "spin configurations have incompatible length: " <> show n <> " != " <> show n'

replicaOverlapSquared :: HasCallStack => Configuration -> Configuration -> Double
replicaOverlapSquared state1@(Configuration n _) state2@(Configuration n' _)
  | n == n' =
      unsafePerformIO $!
        withConfiguration state1 $ \statePtr1 ->
          withConfiguration state2 $ \statePtr2 ->
            pure $
              (/ fromIntegral (n * n)) . realToFrac $
                c_replicaOverlapSquared n statePtr1 statePtr2
  | otherwise = error $ "spin configurations have incompatible length: " <> show n <> " != " <> show n'

-- | Compute a vector of @ΔE@ that one would get by flipping each spin.
energyChanges :: Hamiltonian -> Configuration -> S.Vector ℝ
energyChanges hamiltonian s@(Configuration n _) =
  G.generate n $ \i ->
    unsafePerformIO $!
      withDenseMatrix (hInteraction hamiltonian) $ \couplingsPtr ->
        withVector (hMagneticField hamiltonian) $ \fieldPtr ->
          withConfiguration s $ \statePtr ->
            pure $ c_energyChangeUponFlip n couplingsPtr fieldPtr statePtr i

monteCarloSampling' ::
  forall m a.
  (MonadUnliftIO m, PrimMonad m) =>
  -- | Sweep size
  Int ->
  -- | Number of sweeps between measurements
  Int ->
  -- | Number thermalization steps
  Int ->
  -- | Number measurements
  Int ->
  -- | The Hamiltonian
  Hamiltonian ->
  -- | Inverse temperatures βs
  B.Vector ℝ ->
  -- | What to measure
  FoldM m (ObservableState (Xoshiro256PlusPlus (PrimState m))) a ->
  -- | Random number generator
  Xoshiro256PlusPlusState ->
  m a
monteCarloSampling' sweepSize numberSweeps numberSkip numberMeasure hamiltonian βs fold gₘₐᵢₙ = do
  -- let gs = splitForParallel (G.length βs) gₘₐᵢₙ
  gₘₐᵢₙ' <- thawGen gₘₐᵢₙ
  -- (gs' :: B.Vector (Xoshiro256PlusPlus (PrimState m))) <- G.mapM thawGen gs
  streamFoldM fold . Stream.take numberMeasure . Stream.drop numberSkip $
    monteCarloSampling sweepSize numberSweeps hamiltonian βs gₘₐᵢₙ'

streamFoldM :: Monad m => FoldM m a b -> Stream m a -> m b
streamFoldM (FoldM step begin extract) stream = do
  x0 <- begin
  r <- Stream.foldlM step x0 stream
  extract r

customIterateM :: Monad m => (a -> m a) -> m a -> Stream m a
customIterateM f s₀ = Stream step (Left s₀)
  where
    step (Left !sM) = do
      !s <- sM
      pure $ Yield s (Right s)
    step (Right !s) = do
      !s' <- f s
      pure $ Yield s' (Right s')

monteCarloSampling ::
  forall g m.
  (MonadUnliftIO m, PrimMonad m, g ~ Xoshiro256PlusPlus (PrimState m)) =>
  -- | Sweep size
  Int ->
  -- | Number of sweeps between measurements
  Int ->
  -- | The Hamiltonian
  Hamiltonian ->
  -- | Inverse temperatures βs
  B.Vector ℝ ->
  -- | Random number generator
  g ->
  Stream m (ObservableState g)
monteCarloSampling sweepSize numberSweeps hamiltonian βs gₘₐᵢₙ =
  customIterateM unfoldStep initUnfoldState
  where
    !numReplicas = G.length βs
    initUnfoldState = do
      gₘₐᵢₙ' <- freezeGen gₘₐᵢₙ
      let gs = splitForParallel ShortJump (G.length βs) gₘₐᵢₙ'
      (gs' :: B.Vector (Xoshiro256PlusPlus (PrimState m))) <- G.mapM thawGen gs
      ObservableState
        <$> G.zipWithM (\β g -> randomReplicaExchangeStateM hamiltonian β g) βs gs'
    unfoldStep (ObservableState replicas) = do
      !schedule <- randomReplicaExchangeScheduleM numReplicas numberSweeps gₘₐᵢₙ
      !replicas' <- ObservableState <$> runReplicaExchangeSchedule sweepSize schedule replicas
      pure replicas'

{-
freezeReplicas :: MonadUnliftIO m => B.Vector (ReplicaExchangeState g) -> m ConfigurationBatch
freezeReplicas replicas
  | G.null replicas = error "replicas array should not be empty"
  | otherwise = do
      let (MutableConfiguration numSpins _) = (G.head replicas).state.configuration
          numWords = (numSpins + 63) `div` 64
      buf <- liftIO $ SM.new (G.length replicas * numWords)
      withMutableVector buf $ \bufPtr ->
        G.iforM_ replicas $ \i r ->
          withMutableConfiguration r.state.configuration $ \statePtr ->
            liftIO $
              Ptr.copyPtr
                (bufPtr `Ptr.advancePtr` (i * numWords))
                statePtr
                numWords
      ConfigurationBatch numSpins
        <$> DenseMatrix (G.length replicas) numWords
        <$> liftIO (S.unsafeFreeze buf)
-}

rkkyInteraction :: Double -> Pair Int Int -> Double
rkkyInteraction λ (x :!: y) = 1 / r2 * sin (2 * pi / λ * r)
  where
    r2 = fromIntegral $ x * x + y * y
    r = sqrt r2
{-# INLINE rkkyInteraction #-}

resummedSquare :: Double -> Pair Int Int -> Pair Int Int -> Int -> Double
resummedSquare λ (width :!: height) (x₀ :!: y₀) a =
  side top + side right + side bottom + side left
  where
    top k = (-a + k) :!: a
    bottom k = (a - k) :!: (-a)
    right k = a :!: (a - k)
    left k = (-a) :!: (-a + k)

    shiftedPoint (x :!: y) = x₀ + width * x :!: y₀ + height * y
    side f =
      unId $
        Stream.foldl' (+) 0 $
          Stream.map (rkkyInteraction λ . shiftedPoint) $
            Stream.generate (2 * a) f
{-# SCC resummedSquare #-}

resummedRkkyInteraction :: Monad m => Double -> Pair Int Int -> Pair Int Int -> Stream m Double
resummedRkkyInteraction λ shape p₀ =
  Stream.scanl' (+) (rkkyInteraction λ p₀) resumAllSquares
  where
    square a = resummedSquare λ shape p₀ a
    resumAllSquares =
      let f !a = let !x = square a in Just (x, a + 1)
       in Stream.unfoldr f 1

kolmusRkkyModel :: Text -> IO Hamiltonian
kolmusRkkyModel filename = do
  let parse [i, j, c] = (read (Text.unpack i), read (Text.unpack j), read (Text.unpack c))
      parse x = error $ "expected i,j,coupling; got " <> show x
  tuples <-
    fmap parse
      <$> fmap (Text.split (== ','))
      <$> Text.lines
      <$> Text.readFile (Text.unpack filename)
  let sideLength = (+ 1) . maximum $ fmap (\(i, j, _) -> max i j) tuples
      n = sideLength * sideLength
      magneticField = G.replicate n 0
  buf <- SM.new (n * n)
  loopM 0 (< sideLength) (+ 1) $ \ !x₀ ->
    loopM 0 (< sideLength) (+ 1) $ \ !y₀ -> do
      let !i = y₀ * sideLength + x₀
          write δx δy c =
            let x = (x₀ + δx) `mod` sideLength
                y = (y₀ + δy) `mod` sideLength
                j = y * sideLength + x
             in SM.write buf (i * n + j) c
      forM_ tuples $ \(δx, δy, c) -> do
        write δx δy c
        write δy δx c
  interactionMatrix <- DenseMatrix n n <$> S.unsafeFreeze buf
  pure $ IsingLike interactionMatrix magneticField

ferromagneticIsingModelSquare2D ::
  -- | Length of the lattice
  Int ->
  -- | Magnetic field
  ℝ ->
  -- | The Hamiltonian
  Hamiltonian
ferromagneticIsingModelSquare2D n b = IsingLike interactionMatrix magneticField
  where
    nSpins = n * n
    magneticField = G.replicate nSpins b
    indexToCoord !i = (i `mod` n, i `div` n)
    coordToIndex (!x, !y) = y * n + x
    interactionMatrix = unsafePerformIO $ do
      buf <- SM.new (nSpins * nSpins)
      forM_ [0 .. nSpins - 1] $ \i -> do
        let (x, y) = indexToCoord i
        forM_ (coordToIndex <$> [((x + 1) `mod` n, y), (x, (y + 1) `mod` n)]) $ \j -> do
          SM.write buf (i * nSpins + j) (-0.5)
          SM.write buf (j * nSpins + i) (-0.5)
      DenseMatrix nSpins nSpins <$> S.unsafeFreeze buf

sherringtonKirkpatrickModel ::
  -- | Size of the system
  Int ->
  -- | Seed
  Int ->
  -- | The Hamiltonian
  Hamiltonian
sherringtonKirkpatrickModel n seed = IsingLike interactionMatrix magneticField
  where
    magneticField = G.replicate n 0
    interactionMatrix = unsafePerformIO $ do
      g <- mkXoshiro256PlusPlus seed
      buf <- SM.new (n * n)
      forM_ [0 .. n - 1] $ \i ->
        forM_ [i + 1 .. n - 1] $ \j -> do
          h <- realToFrac <$> MWC.standard g
          SM.write buf (i * n + j) h
          SM.write buf (j * n + i) h
      DenseMatrix n n <$> S.unsafeFreeze buf

{-
shuffleVectorM :: (PrimMonad m, StatefulGen g m, GM.MVector v a) => v (PrimState m) a -> g -> m ()
shuffleVectorM !v !g = go (GM.length v - 1)
  where
    swap !i !j = do
      a <- GM.read v i
      b <- GM.read v j
      GM.write v i b
      GM.write v j a
    go !j
      | j > 0 = do
          swap j =<< uniformRM (0, j - 1) g
          go (j - 1)
      | otherwise = pure ()
-}

randomIsingModel3D ::
  -- | Side length of the lattice
  Int ->
  -- | Seed
  Int ->
  -- | The Hamiltonian
  Hamiltonian
randomIsingModel3D n seed = IsingLike interactionMatrix magneticField
  where
    nSpins = n * n * n
    magneticField = G.replicate nSpins 0
    indexToCoord !i =
      let z = i `div` (n * n)
          y = (i `mod` (n * n)) `div` n
          x = (i `mod` (n * n)) `mod` n
       in (x, y, z)
    coordToIndex (!x, !y, !z) = z * n * n + y * n + x
    interactionMatrix = unsafePerformIO $ do
      g <- mkXoshiro256PlusPlus seed
      -- couplings <- G.unsafeThaw $ U.generate (3 * numSpins) (\i -> fromIntegral $ 2 * (i `mod` 2) - 1)
      -- shuffleVectorM couplings g
      couplings <- U.replicateM (3 * nSpins) (realToFrac . (/ 2) <$> MWC.standard g)

      buf <- SM.new (nSpins * nSpins)
      forM_ [0 .. nSpins - 1] $ \i -> do
        let (x, y, z) = indexToCoord i
            neighbors =
              G.map coordToIndex $
                B.fromList
                  [ ((x + 1) `mod` n, y, z),
                    (x, (y + 1) `mod` n, z),
                    (x, y, (z + 1) `mod` n)
                  ]
        G.iforM_ neighbors $ \k j -> do
          let h = couplings ! (3 * i + k)
          GM.write buf (i * nSpins + j) h
          GM.write buf (j * nSpins + i) h

      DenseMatrix nSpins nSpins <$> G.unsafeFreeze buf

allConfigurations :: MonadUnliftIO m => Int -> Stream m Configuration
allConfigurations n
  | n >= 32 = error "too many spins"
  | otherwise = Stream.unfoldr step (0 :: Word64)
  where
    toConfiguration !i = Configuration n (G.singleton i)
    step !i
      | i < bit n = Just (toConfiguration i, i + 1)
      | otherwise = Nothing

-- generateSchedule :: StategulGen g m => Int -> Int ->

--
-- data MutableConfigurationBatch s
--   = MutableConfigurationBatch
--       {-# UNPACK #-} !Int
--       {-# UNPACK #-} !(DenseMatrix (S.MVector s) Word64)
--
withVector ::
  (MonadUnliftIO m, Storable a) =>
  S.Vector a ->
  (Ptr a -> m b) ->
  m b
withVector v action = do
  let (fp, _) = S.unsafeToForeignPtr0 v
  !r <- action (UnliftIO.unsafeForeignPtrToPtr fp)
  UnliftIO.touchForeignPtr fp
  pure r
{-# INLINE withVector #-}

withMutableVector ::
  MonadUnliftIO m =>
  S.MVector (PrimState IO) a ->
  (Ptr a -> m b) ->
  m b
withMutableVector (S.MVector _ fp) action = do
  !r <- action (UnliftIO.unsafeForeignPtrToPtr fp)
  UnliftIO.touchForeignPtr fp
  pure r
{-# INLINE withMutableVector #-}

-- withCouplings :: MonadUnliftIO m => Couplings -> (Ptr ℝ -> m a) -> m a
-- withCouplings (Couplings (DenseMatrix _ _ v)) = withVector v

withDenseMatrix :: (MonadUnliftIO m, Storable a) => DenseMatrix S.Vector a -> (Ptr a -> m b) -> m b
withDenseMatrix (DenseMatrix _ _ v) = withVector v
{-# INLINE withDenseMatrix #-}

withConfiguration ::
  MonadUnliftIO m =>
  Configuration ->
  (Ptr Word64 -> m a) ->
  m a
withConfiguration (Configuration _ v) = withVector v
{-# INLINE withConfiguration #-}

withMutableConfiguration ::
  MonadUnliftIO m =>
  MutableConfiguration ->
  (Ptr Word64 -> m a) ->
  m a
withMutableConfiguration (MutableConfiguration _ (S.MVector _ fp)) action = do
  !r <- action (UnliftIO.unsafeForeignPtrToPtr fp)
  UnliftIO.touchForeignPtr fp
  pure r
{-# INLINE withMutableConfiguration #-}

instance Semigroup SweepStats where
  (<>) (SweepStats t1 r1) (SweepStats t2 r2) = SweepStats t (s / fromIntegral t)
    where
      !t = t1 + t2
      !s = fromIntegral t1 * r1 + fromIntegral t2 * r2
  {-# INLINE (<>) #-}

instance Monoid SweepStats where
  mempty = SweepStats 0 0
  {-# INLINE mempty #-}

-- doSweep ::
--   MonadUnliftIO m =>
--   ProbDist ->
--   Int ->
--   Ptr Int ->
--   Ptr ℝ ->
--   MetropolisState ->
--   m SweepStats
-- doSweep (ProbDist β hamiltonian) !numberSteps !randIntsPtr !randFloatsPtr !state = do
--   -- liftIO $ putStrLn "Calling doSweep..."
--   let !numberBits = G.length (hMagneticField hamiltonian)
--   withDenseMatrix (hInteraction hamiltonian) $ \couplingsPtr ->
--     withVector (hMagneticField hamiltonian) $ \fieldPtr ->
--       withMutableConfiguration state.configuration $ \statePtr ->
--         withMutableVector state.deltaEnergies $ \deltaEnergiesPtr ->
--           liftIO $ do
--             -- x <- unsafeFreezeConfiguration state.configuration
--             -- putStrLn $ "Starting sweep from " <> show x
--             -- print =<< peek statePtr
--             e₀ <- readIORef state.energy
--             with e₀ $ \energyPtr -> do
--               acceptance <-
--                 {-# SCC c_runOneSweep #-}
--                 c_runOneSweep
--                   numberBits
--                   numberSteps
--                   β
--                   couplingsPtr
--                   randIntsPtr
--                   randFloatsPtr
--                   statePtr
--                   energyPtr
--                   deltaEnergiesPtr
--               !e <- peek energyPtr
--               writeIORef state.energy e
--               -- !e' <- totalEnergy hamiltonian <$> unsafeFreezeConfiguration state.configuration
--               -- print acceptance
--               -- assert (e == e') $
--               pure $ SweepStats numberSteps acceptance
-- {-# INLINE doSweep #-}

doManySweeps ::
  (PrimMonad m, MonadUnliftIO m) =>
  ProbDist ->
  Int ->
  Int ->
  MetropolisState ->
  Xoshiro256PlusPlus (PrimState m) ->
  m SweepStats
doManySweeps (ProbDist β hamiltonian) numberSweeps sweepSize s (Xoshiro256PlusPlus g) = do
  -- k <- unsafeFreezeConfiguration s.configuration
  -- g' <- freezeGen gen
  -- liftIO . putStrLn $ "doManySweeps " <> show (β, numberSweeps, sweepSize, k, g')
  let numberSpins = G.length hamiltonian.hMagneticField
      totalSteps = numberSweeps * sweepSize
      gPtr = castPtr (mutableByteArrayContents g)
  withDenseMatrix hamiltonian.hInteraction $ \couplingsPtr ->
    withMutableConfiguration s.configuration $ \statePtr ->
      withMutableVector s.deltaEnergies $ \deltaEnergiesPtr -> do
        e₀ <- liftIO $ readPVar s.energy
        UnliftIO.with e₀ $ \energyPtr -> do
          !acceptance <-
            liftIO $ do
              r <-
                {-# SCC c_runOneSweep #-}
                c_runOneSweep
                  numberSpins
                  totalSteps
                  β
                  couplingsPtr
                  statePtr
                  deltaEnergiesPtr
                  gPtr
                  energyPtr
              touch g
              pure r
          liftIO $ peek energyPtr >>= writePVar s.energy
          pure $ SweepStats totalSteps acceptance

-- | Clone a spin configuration into a mutable one
thawConfiguration :: MonadUnliftIO m => Configuration -> m MutableConfiguration
thawConfiguration (Configuration n v) = withRunInIO $ \_ ->
  MutableConfiguration n <$> G.thaw v
{-# INLINE thawConfiguration #-}

{-
-- | Freeze a spin configuration into an immutable one
freezeConfiguration :: MonadUnliftIO m => MutableConfiguration -> m Configuration
freezeConfiguration (MutableConfiguration n v) = withRunInIO $ \_ ->
  Configuration n <$> G.freeze v
{-# INLINE freezeConfiguration #-}
-}

{-
-- | Unsafe version of 'thawConfiguration' which avoids copying. Using this function is only legal
-- if the original spin configuration will never be used again.
unsafeThawConfiguration :: MonadUnliftIO m => Configuration -> m MutableConfiguration
unsafeThawConfiguration (Configuration n v) = withRunInIO $ \_ ->
  MutableConfiguration n <$> G.unsafeThaw v
{-# INLINE unsafeThawConfiguration #-}
-}

-- | Unsafe version of 'freezeConfiguration' which avoids copying. Using this function is only legal
-- if the original spin configuration will never be modified again.
unsafeFreezeConfiguration :: MonadUnliftIO m => MutableConfiguration -> m Configuration
unsafeFreezeConfiguration (MutableConfiguration n v) = withRunInIO $ \_ ->
  Configuration n <$> G.unsafeFreeze v
{-# INLINE unsafeFreezeConfiguration #-}

autocorrelation :: (G.Vector v a, G.Vector v ℝ, G.Vector v ℂ, Real a) => v a -> v ℝ
autocorrelation v = G.map (/ G.head autocorr) autocorr
  where
    size = G.length v
    n = (1 `shiftL`) . ceiling . logBase 2 . (fromIntegral :: Int -> Double) $ G.length v
    preprocessed = G.map (\x -> realToFrac x :+ 0) v G.++ G.replicate (2 * n - size) 0
    f = fft preprocessed
    autocorr = G.map realPart . G.take size . ifft $ G.map (\x -> x * conjugate x) f

autocorrTime :: (G.Vector v ℝ, G.Vector v (Int, ℝ)) => v ℝ -> ℝ
autocorrTime f =
  case G.find (\(i, τ) -> fromIntegral i < c * τ) $ G.indexed τs of
    Just (_, τ) -> τ
    Nothing -> G.last τs
  where
    τs = G.map (\x -> 2 * x - 1) $ G.scanl1' (+) f
    c = 5

-- r"""Estimate the normalized autocorrelation function of a 1D array."""
-- if x.ndim != 1:
--     raise ValueError("x has wrong shape: {}; expected a 1D array".format(x.shape))
-- n = 1 << math.ceil(math.log2(len(x)))
-- f = np.fft.fft(x - np.mean(x), n=2 * n)
-- autocorr = np.fft.ifft(f * np.conj(f))[: len(x)].real
-- autocorr /= autocorr[0]
-- return autocorr

genericFFT ::
  G.Vector v (Complex Float) =>
  (Int -> Ptr (Complex Float) -> Ptr (Complex Float) -> IO ()) ->
  v (Complex Float) ->
  v (Complex Float)
genericFFT f v =
  unsafePerformIO $ do
    signal <- S.unsafeThaw (G.convert v)
    output <- SM.unsafeNew n
    SM.unsafeWith signal $ \signalPtr ->
      SM.unsafeWith output $ \outputPtr ->
        f n signalPtr outputPtr
    G.convert <$> S.unsafeFreeze output
  where
    n = G.length v

fft :: G.Vector v (Complex Float) => v (Complex Float) -> v (Complex Float)
fft = genericFFT c_fft

ifft :: G.Vector v (Complex Float) => v (Complex Float) -> v (Complex Float)
ifft = genericFFT c_ifft

{-
unsafeFreezeDenseMatrix ::
  (MonadUnliftIO m, GM.MVector (G.Mutable v) a, G.Vector v a) =>
  DenseMatrix (G.Mutable v RealWorld) a ->
  m (DenseMatrix v a)
unsafeFreezeDenseMatrix (DenseMatrix nRows nCols v) =
  DenseMatrix nRows nCols <$> liftIO (G.unsafeFreeze v)
{-# INLINE unsafeFreezeDenseMatrix #-}
-}
foreign import ccall unsafe "fft.h fft"
  c_fft :: Int -> Ptr (Complex Float) -> Ptr (Complex Float) -> IO ()

foreign import ccall unsafe "fft.h ifft"
  c_ifft :: Int -> Ptr (Complex Float) -> Ptr (Complex Float) -> IO ()

foreign import ccall unsafe "metropolis.h energy_change_upon_flip"
  c_energyChangeUponFlip :: Int -> Ptr Float -> Ptr Float -> Ptr Word64 -> Int -> Float

foreign import ccall unsafe "metropolis.h total_energy"
  c_totalEnergy :: Int -> Ptr Float -> Ptr Float -> Ptr Word64 -> Float

foreign import ccall unsafe "metropolis.h total_magnetization"
  c_totalMagnetization :: Int -> Ptr Word64 -> Int

foreign import ccall unsafe "metropolis.h replica_overlap"
  c_replicaOverlap :: Int -> Ptr Word64 -> Ptr Word64 -> Float

foreign import ccall unsafe "metropolis.h replica_overlap_squared"
  c_replicaOverlapSquared :: Int -> Ptr Word64 -> Ptr Word64 -> Float

foreign import ccall unsafe "metropolis.h add_spin_spin_correlation"
  c_addSpinSpinCorrelation :: Int -> Ptr Word64 -> Ptr Float -> IO ()

foreign import ccall unsafe "metropolis.h add_magnetization"
  c_addMagnetization :: Int -> Ptr Word64 -> Ptr Float -> IO ()

foreign import ccall unsafe "metropolis.h run_one_sweep"
  c_runOneSweep :: Int -> Int -> Float -> Ptr Float -> Ptr Word64 -> Ptr Float -> Ptr Word64 -> Ptr Double -> IO Double
