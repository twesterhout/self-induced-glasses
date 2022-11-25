{-# LANGUAGE CApiFFI #-}
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
    fourierTransformStructureFactorSquare,
    monteCarloSampling',
    monteCarloSampling,
    energyChanges,
    maybeExchangeReplicas,
    prettySchedule,
    writeMatrixToCsv,
    writeVectorToCsv,
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
import Control.Monad (forM, forM_, unless, (<=<))
import Control.Monad.IO.Unlift
import Control.Monad.Primitive
import Control.Monad.ST.Strict (runST)
import Data.Bits
import qualified Data.ByteString.Builder as Builder
import qualified Data.ByteString.Builder.RealFloat as Builder
import Data.IORef
import qualified Data.List
import Data.Primitive.ByteArray (mutableByteArrayContents)
import qualified Data.Primitive.Ptr as Ptr
import Data.Stream.Monadic (Step (..), Stream (..))
import qualified Data.Stream.Monadic as Stream
import Data.Text (Text)
import qualified Data.Text as Text
import Data.Vector (Vector)
import qualified Data.Vector as B
import Data.Vector.Generic ((!))
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Storable as S
import qualified Data.Vector.Storable.Mutable as SM
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import Data.Word
import Debug.Trace
import Foreign.Marshal.Utils
import Foreign.Ptr
import Foreign.Storable
import qualified GHC.Exts as GHC (IsList (..))
import GHC.Generics
import GHC.Stack (HasCallStack)
import ListT (ListT)
import qualified ListT
import SelfInducedGlasses.Random
import System.IO
import System.IO.Unsafe (unsafePerformIO)
import qualified System.Random.MWC.Distributions as MWC
import System.Random.Stateful
import Text.PrettyPrint.ANSI.Leijen (Doc, Pretty (..))
import qualified Text.PrettyPrint.ANSI.Leijen as Pretty
import UnliftIO.Async
import UnliftIO.Exception (assert)
import qualified UnliftIO.Foreign as UnliftIO
import qualified UnliftIO.IORef as UnliftIO

-- | Real number
type ℝ = Float

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
    hInteraction :: !(DenseMatrix S.Vector ℝ),
    -- | External magnetic field \(B_i\).
    hMagneticField :: !(S.Vector ℝ)
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
    energy :: {-# UNPACK #-} !(IORef Double),
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
data ProbDist = ProbDist {pdBeta :: !ℝ, pdHamiltonian :: !Hamiltonian}
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

data ReplicaExchangeSchedule = ReplicaExchangeSchedule
  { resIntervals :: !(B.Vector [(Int, Int)]),
    resSwaps :: !(B.Vector (Int, ℝ))
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

newtype ObservableState g = ObservableState {unObservableState :: (B.Vector (ReplicaExchangeState g))}
  deriving stock (Generic)
  deriving anyclass (NFData)

updateVarianceAcc :: Fractional a => (Int, a, a) -> a -> (Int, a, a)
updateVarianceAcc (!n, !μ, !m2) !x = (n', μ', m2')
  where
    !n' = n + 1
    !μ' = μ + δx / fromIntegral n'
    !m2' = m2 + δx * (x - μ')
    !δx = x - μ
{-# INLINE updateVarianceAcc #-}

data MeanVariance a = MeanVariance
  { mean :: {-# UNPACK #-} !a,
    variance :: {-# UNPACK #-} !a
  }
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
    (ColorStats numBlue numRed) = g.hist
    hist' = case g.color of
      Just ReplicaBlue -> ColorStats (numBlue + 1) numRed
      Just ReplicaRed -> ColorStats numBlue (numRed + 1)
      Nothing -> ColorStats numBlue numRed
{-# INLINE updateColorStats #-}

updateTemperatures :: Vector ℝ -> Vector ColorStats -> Vector ℝ
updateTemperatures ts stats = G.generate (G.length ts) getTemperature
  where
    m = fromIntegral (G.length ts)
    -- Calculate f(Tᵢ)
    fs = G.map (\(ColorStats numBlue numRed) -> fromIntegral numBlue / fromIntegral (numBlue + numRed)) stats
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

data PiecewiseLinear = PiecewiseLinear [(ℝ, ℝ)]

data PiecewiseConstant = PiecewiseConstant !(Vector ℝ) !(Vector ℝ)

integratePiecewiseConst :: PiecewiseConstant -> PiecewiseLinear
integratePiecewiseConst (PiecewiseConstant xs₀ ys₀) = PiecewiseLinear $ go 0 (G.toList xs₀) (G.toList ys₀)
  where
    go !acc (!x : xnext : xrest) (y : yrest) =
      let area = y * (xnext - x)
       in (x, acc) : go (acc + area) (xnext : xrest) yrest
    go acc (x : []) [] = [(x, acc)]
    go _ _ _ = error "invalid PiecewiseConstant"

invertPiecewiseLinear :: PiecewiseLinear -> PiecewiseLinear
invertPiecewiseLinear (PiecewiseLinear points) = PiecewiseLinear $ fmap (\(x, y) -> (y, x)) points

-- evaluatePiecewiseLinear :: PiecewiseLinear -> ℝ -> ℝ
-- evaluatePiecewiseLinear (PiecewiseLinear points) x = undefined
--   where
--     go

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
      secondMomentFold (\r -> liftIO $ readIORef r.state.energy)

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
      (fromIntegral . abs . totalMagnetization)
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
  forall m a g.
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
    go !accs !state !i
      | i < G.length accs = do
          x@(Configuration n _) <-
            unsafeFreezeConfiguration (state.unObservableState ! i).state.configuration
          withConfiguration x $ \xPtr ->
            withMutableVector (accs ! i) $ \outPtr ->
              liftIO $ c_addMagnetization n xPtr outPtr
          go accs state (i + 1)
      | otherwise = pure ()
    step (!k, !maybeAccs) !state = do
      accs <- case maybeAccs of
        Just accs -> pure accs
        Nothing ->
          let numReplicas = G.length . unObservableState $ state
              (MutableConfiguration numSpins _) =
                (G.head state.unObservableState).state.configuration
           in G.replicateM numReplicas $ liftIO (GM.new numSpins)
      go accs state 0 >> pure (k + 1, Just accs)
    extract (!k, !maybeAccs) = do
      case maybeAccs of
        Just accs ->
          G.forM accs $ \v ->
            G.map (/ fromIntegral k) <$> liftIO (G.unsafeFreeze v)
        Nothing -> pure G.empty

structureFactorFold ::
  forall m a g.
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
    go !accs !state !i
      | i < G.length accs = do
          x <- unsafeFreezeConfiguration (state.unObservableState ! i).state.configuration
          let (DenseMatrix n _ out) = accs ! i
          withConfiguration x $ \xPtr ->
            withMutableVector out $ \outPtr ->
              liftIO $ c_addSpinSpinCorrelation n xPtr outPtr
          go accs state (i + 1)
      | otherwise = pure ()
    step (!k, !maybeAccs) !state = do
      accs <- case maybeAccs of
        Just accs -> pure accs
        Nothing ->
          let numReplicas = G.length . unObservableState $ state
              (MutableConfiguration numSpins _) =
                (G.head state.unObservableState).state.configuration
           in G.replicateM numReplicas $
                DenseMatrix numSpins numSpins <$> liftIO (GM.new (numSpins * numSpins))
      go accs state 0 >> pure (k + 1, Just accs)
    extract (!k, !maybeAccs) = do
      case maybeAccs of
        Just accs ->
          G.forM accs $ \(DenseMatrix nRows nCols v) ->
            DenseMatrix nRows nCols . G.map (/ fromIntegral k)
              <$> liftIO (G.unsafeFreeze v)
        Nothing -> pure G.empty

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
    !numSpins = G.length spin
    !sideLength = round $ sqrt (fromIntegral numSpins)
    !xs = U.generate numSpins (`mod` sideLength)
    !ys = U.generate numSpins (`div` sideLength)
    term :: ℝ -> ℝ -> Int -> Int -> ℝ
    term !qx !qy !i !j = (1 / fromIntegral numSpins) * cos (qx * δx + qy * δy) * (sᵢsⱼ - sᵢ * sⱼ)
      where
        !δx = fromIntegral $ G.unsafeIndex xs i - G.unsafeIndex xs j
        !δy = fromIntegral $ G.unsafeIndex ys i - G.unsafeIndex ys j
        !sᵢsⱼ = indexDenseMatrix spinSpin i j
        !sᵢ = G.unsafeIndex spin i
        !sⱼ = G.unsafeIndex spin j
    computeOne !qx !qy =
      iFold 0 (< numSpins) (+ 1) 0 $ \ !accᵢ !i ->
        (accᵢ +) $
          iFold 0 (< numSpins) (+ 1) 0 $ \ !accⱼ !j ->
            accⱼ + term qx qy i j

loopM :: Monad m => i -> (i -> Bool) -> (i -> i) -> (i -> m ()) -> m ()
loopM i₀ cond inc action = go i₀
  where
    go !i
      | cond i = do () <- action i; go (inc i)
      | otherwise = pure ()
{-# INLINE loopM #-}

iFoldM :: Monad m => i -> (i -> Bool) -> (i -> i) -> a -> (a -> i -> m a) -> m a
iFoldM i₀ cond inc x₀ action = go x₀ i₀
  where
    go !x !i
      | cond i = do !x' <- action x i; go x' (inc i)
      | otherwise = pure x
{-# INLINE iFoldM #-}

iFold :: i -> (i -> Bool) -> (i -> i) -> a -> (a -> i -> a) -> a
iFold i₀ cond inc x₀ action = go x₀ i₀
  where
    go !x !i
      | cond i = let !x' = action x i in go x' (inc i)
      | otherwise = x
{-# INLINE iFold #-}

-- | Get matrix shape
dmShape :: DenseMatrix v a -> (Int, Int)
dmShape m = (m.numRows, m.numCols)

denseMatrixFromList :: G.Vector v a => [[a]] -> DenseMatrix v a
denseMatrixFromList rs
  | G.length elements == nRows * nCols = DenseMatrix nRows nCols elements
  | otherwise = error "nested list has irregular shape"
  where
    !nRows = length rs
    !nCols = case rs of
      [] -> 0
      (r : _) -> length r
    elements = G.fromListN (nRows * nCols) $ mconcat rs

denseMatrixToList :: G.Vector v a => DenseMatrix v a -> [[a]]
denseMatrixToList (DenseMatrix _ nCols v) = go (G.toList v)
  where
    go elements = case splitAt nCols elements of
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
writeMatrixToCsv filename matrix@(DenseMatrix nRows nCols v) =
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

indexConfiguration :: ConfigurationBatch -> Int -> Configuration
indexConfiguration (ConfigurationBatch numSpins (DenseMatrix numReplicas numWords v)) i =
  Configuration numSpins (G.slice (i * numWords) numWords v)

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
    !δe <- liftIO $ do
      !e1 <- readIORef s1.energy
      !e2 <- readIORef s2.energy
      pure $ realToFrac (e2 - e1)
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
  (MonadUnliftIO m, g ~ Xoshiro256PlusPlus (PrimState m)) =>
  Int ->
  ReplicaExchangeSchedule ->
  B.Vector (ReplicaExchangeState g) ->
  m (B.Vector (ReplicaExchangeState g))
runReplicaExchangeSchedule sweepSize schedule₀ states₀ = do
  -- liftIO $ do
  --   putStrLn "Running schedule ..."
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
    spawn !replica@(ReplicaExchangeState probDist state _ _ stats g) !intervals =
      case intervals of
        ((!start, !end) : others) -> do
          !future <- async $ do
            !stats' <- doManySweeps probDist (end - start) sweepSize state g
            pure . force $ replica {stats = stats <> stats'}
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
    step accumulators (i, u)
      | G.length accumulators > 1 = do
          let n = G.length accumulators
              (future1, intervals1) = accumulators ! i
              (future2, intervals2) = accumulators ! (i + 1)
          state1 <- wait future1
          state2 <- wait future2
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
    extract ::
      B.Vector (Async (ReplicaExchangeState g), [(Int, Int)]) ->
      m (B.Vector (ReplicaExchangeState g))
    extract accumulators =
      G.forM accumulators $ \(future, intervals) -> do
        unless (null intervals) $ error "not all intervals have been processed"
        wait future

-- 🡼 🡽 🡾 🡿
-- "╸╺"
--
-- rowIdx: 0   ━━🡾🡽━━╸╺━━╸╺━━🡾🡽  [(0, 1), (1, 4)]
--         1   ━━🡽🡾━━🡾🡽━━🡾🡽━━🡽🡾  [(0, 1), (1, 2), (2, 3), (3, 4)]
--         2   ━━╸╺━━🡽🡾━━🡽🡾━━    [(0, 2), (2, 3), (3, 4)]
--              [0,  1,  1,  0]

prettyInterval :: G.Vector v (Int, ℝ) => Int -> (Int, Int) -> v (Int, ℝ) -> Doc
prettyInterval rowIdx (start, end) swaps =
  Pretty.text left <> Pretty.hcat middle <> Pretty.text right
  where
    middle = Pretty.punctuate (Pretty.text "╸╺") (replicate (end - start) (Pretty.text "━━"))
    left
      | start <= 0 = "╺"
      | fst (swaps ! (start - 1)) == rowIdx = "🡽"
      | fst (swaps ! (start - 1)) == rowIdx - 1 = "🡾"
      | otherwise = "╺"
    right
      | end <= 0 = "╸"
      | fst (swaps ! (end - 1)) == rowIdx = "🡾"
      | fst (swaps ! (end - 1)) == rowIdx - 1 = "🡽"
      | otherwise = "╸"

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

mkIntervals :: G.Vector v (Int, ℝ) => v (Int, ℝ) -> Int -> [(Int, Int)]
mkIntervals swaps replicaIdx = go [(0, G.length swaps)] (G.length swaps - 2)
  where
    go acc@((!start, !end) : rest) !i
      | i <= 0 = acc
      | fst (swaps ! (i - 1)) == replicaIdx
          || fst (swaps ! (i - 1)) == replicaIdx - 1 =
          go ((start, i) : (i, end) : rest) (i - 1)
      | otherwise = go acc (i - 1)
    go _ _ = error "this should never happen by construction"

mkSchedule :: G.Vector v (Int, ℝ) => Int -> v (Int, ℝ) -> ReplicaExchangeSchedule
mkSchedule numReplicas swaps =
  ReplicaExchangeSchedule
    (G.generate numReplicas (mkIntervals swaps))
    (G.convert swaps)

randomReplicaExchangeScheduleM :: StatefulGen g m => Int -> Int -> g -> m ReplicaExchangeSchedule
randomReplicaExchangeScheduleM numReplicas numSwaps g = do
  swaps <-
    U.replicateM numSwaps $
      (,)
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

-- | Generate a random 'MetropolisState'.
randomMetropolisStateM ::
  (MonadUnliftIO m, StatefulGen g m) =>
  -- | The Hamiltonian
  Hamiltonian ->
  -- | Random number generator
  g ->
  m MetropolisState
randomMetropolisStateM hamiltonian g = do
  spins <- randomConfigurationM numSpins g
  MetropolisState
    <$> thawConfiguration spins
    <*> mkCurrentEnergy spins
    <*> mkDeltaEnergies spins
  where
    numSpins = G.length . hMagneticField $ hamiltonian
    mkCurrentEnergy spins = liftIO $ newIORef (totalEnergy hamiltonian spins)
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
  state <- randomMetropolisStateM hamiltonian g
  pure $ ReplicaExchangeState (ProbDist β hamiltonian) state Nothing mempty mempty g

-- | Compute energy of a spin configuration
totalEnergy :: Hamiltonian -> Configuration -> Double
totalEnergy hamiltonian state@(Configuration n _) =
  unsafePerformIO $!
    withDenseMatrix (hInteraction hamiltonian) $
      \interactionPtr ->
        withVector (hMagneticField hamiltonian) $ \fieldPtr ->
          withConfiguration state $ \statePtr ->
            pure $ realToFrac (c_totalEnergy n interactionPtr fieldPtr statePtr)

-- | Compute the total magnetization of a spin configuration
totalMagnetization :: Configuration -> Int
totalMagnetization state@(Configuration n _) =
  unsafePerformIO $!
    withConfiguration state $
      \statePtr ->
        pure $ c_totalMagnetization n statePtr

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
energyChanges hamiltonian state@(Configuration n _) =
  G.generate n $ \i ->
    unsafePerformIO $!
      withDenseMatrix (hInteraction hamiltonian) $
        \couplingsPtr ->
          withVector (hMagneticField hamiltonian) $ \fieldPtr ->
            withConfiguration state $ \statePtr ->
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
      let gs = splitForParallel (G.length βs) gₘₐᵢₙ'
      (gs' :: B.Vector (Xoshiro256PlusPlus (PrimState m))) <- G.mapM thawGen gs
      ObservableState
        <$> G.zipWithM (\β g -> randomReplicaExchangeStateM hamiltonian β g) βs gs'
    unfoldStep (ObservableState replicas) = do
      !schedule <- randomReplicaExchangeScheduleM numReplicas numberSweeps gₘₐᵢₙ
      !replicas' <- ObservableState <$> runReplicaExchangeSchedule sweepSize schedule replicas
      pure replicas'

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

ferromagneticIsingModelSquare2D ::
  HasCallStack =>
  -- | Length of the lattice
  Int ->
  -- | Magnetic field
  ℝ ->
  -- | The Hamiltonian
  Hamiltonian
ferromagneticIsingModelSquare2D n b = IsingLike interactionMatrix magneticField
  where
    numSpins = n * n
    magneticField = G.replicate numSpins b
    indexToCoord !i = (i `mod` n, i `div` n)
    coordToIndex (!x, !y) = y * n + x
    interactionMatrix = unsafePerformIO $ do
      buf <- SM.new (numSpins * numSpins)
      forM_ [0 .. numSpins - 1] $ \i -> do
        let (x, y) = indexToCoord i
        forM_ (coordToIndex <$> [((x + 1) `mod` n, y), (x, (y + 1) `mod` n)]) $ \j -> do
          SM.write buf (i * numSpins + j) (-0.5)
          SM.write buf (j * numSpins + i) (-0.5)
      DenseMatrix numSpins numSpins <$> S.unsafeFreeze buf

sherringtonKirkpatrickModel ::
  HasCallStack =>
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

randomIsingModel3D ::
  HasCallStack =>
  -- | Side length of the lattice
  Int ->
  -- | Seed
  Int ->
  -- | The Hamiltonian
  Hamiltonian
randomIsingModel3D n seed = IsingLike interactionMatrix magneticField
  where
    numSpins = n * n * n
    magneticField = G.replicate numSpins 0
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
      couplings <- U.replicateM (3 * numSpins) (realToFrac . (/ 2) <$> MWC.standard g)

      buf <- SM.new (numSpins * numSpins)
      forM_ [0 .. numSpins - 1] $ \i -> do
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
          GM.write buf (i * numSpins + j) h
          GM.write buf (j * numSpins + i) h

      DenseMatrix numSpins numSpins <$> G.unsafeFreeze buf

allConfigurations :: MonadUnliftIO m => Int -> ListT m Configuration
allConfigurations n
  | n >= 32 = error "too many spins"
  | otherwise = ListT.unfold step (0 :: Word64)
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
  (MonadUnliftIO m, Storable a) =>
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
  MonadUnliftIO m =>
  ProbDist ->
  Int ->
  Int ->
  MetropolisState ->
  Xoshiro256PlusPlus (PrimState m) ->
  m SweepStats
doManySweeps (ProbDist β hamiltonian) numberSweeps sweepSize state (Xoshiro256PlusPlus g) = do
  let numberSpins = G.length hamiltonian.hMagneticField
      totalSteps = numberSweeps * sweepSize
      gPtr = castPtr (mutableByteArrayContents g)
  withDenseMatrix hamiltonian.hInteraction $ \couplingsPtr ->
    withMutableConfiguration state.configuration $ \statePtr ->
      withMutableVector state.deltaEnergies $ \deltaEnergiesPtr -> do
        e₀ <- UnliftIO.readIORef state.energy
        UnliftIO.with e₀ $ \energyPtr -> do
          !acceptance <-
            liftIO $
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
          !e <- liftIO $ peek energyPtr
          UnliftIO.writeIORef state.energy e
          pure $ SweepStats totalSteps acceptance

-- | Clone a spin configuration into a mutable one
thawConfiguration :: MonadUnliftIO m => Configuration -> m MutableConfiguration
thawConfiguration (Configuration n v) = withRunInIO $ \_ ->
  MutableConfiguration n <$> G.thaw v
{-# INLINE thawConfiguration #-}

-- | Freeze a spin configuration into an immutable one
freezeConfiguration :: MonadUnliftIO m => MutableConfiguration -> m Configuration
freezeConfiguration (MutableConfiguration n v) = withRunInIO $ \_ ->
  Configuration n <$> G.freeze v
{-# INLINE freezeConfiguration #-}

-- | Unsafe version of 'thawConfiguration' which avoids copying. Using this function is only legal
-- if the original spin configuration will never be used again.
unsafeThawConfiguration :: MonadUnliftIO m => Configuration -> m MutableConfiguration
unsafeThawConfiguration (Configuration n v) = withRunInIO $ \_ ->
  MutableConfiguration n <$> G.unsafeThaw v
{-# INLINE unsafeThawConfiguration #-}

-- | Unsafe version of 'freezeConfiguration' which avoids copying. Using this function is only legal
-- if the original spin configuration will never be modified again.
unsafeFreezeConfiguration :: MonadUnliftIO m => MutableConfiguration -> m Configuration
unsafeFreezeConfiguration (MutableConfiguration n v) = withRunInIO $ \_ ->
  Configuration n <$> G.unsafeFreeze v
{-# INLINE unsafeFreezeConfiguration #-}

unsafeFreezeDenseMatrix ::
  (MonadUnliftIO m, GM.MVector (G.Mutable v) a, G.Vector v a) =>
  DenseMatrix (G.Mutable v RealWorld) a ->
  m (DenseMatrix v a)
unsafeFreezeDenseMatrix (DenseMatrix nRows nCols v) =
  DenseMatrix nRows nCols <$> liftIO (G.unsafeFreeze v)
{-# INLINE unsafeFreezeDenseMatrix #-}

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
