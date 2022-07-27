{-# LANGUAGE CApiFFI #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveGeneric #-}
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
    Configuration (..),
    ConfigurationBatch (..),
    MutableConfiguration (..),
    SweepStats (..),
    ProbDist (..),
    ReplicaExchangeSchedule (..),
    ReplicaExchangeState (..),
    indexConfiguration,
    totalEnergy,
    totalMagnetization,
    energyFold,
    magnetizationFold,
    structureFactorFold,
    monteCarloSampling',
    monteCarloSampling,
    energyChanges,
    maybeExchangeReplicas,
    prettySchedule,
    mkIntervals,
    mkSchedule,
    withVector,
    withMutableConfiguration,
    doSweep,
    doManySweeps,
    allConfigurations,
    randomConfigurationM,
  )
where

import Control.DeepSeq
import Control.Foldl (FoldM (..))
import qualified Control.Foldl as Foldl
import Control.Monad (forM, forM_, unless)
import Control.Monad.IO.Unlift
import Control.Monad.Primitive
import Data.Bits
import Data.IORef
import qualified Data.Primitive.Ptr as Ptr
import qualified Data.Vector as B
import Data.Vector.Generic ((!))
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Storable as S
import qualified Data.Vector.Storable.Mutable as SM
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import Data.Word
import Foreign.Marshal.Utils
import Foreign.Ptr
import Foreign.Storable
import qualified GHC.Exts as GHC (IsList (..))
import GHC.Generics
import GHC.Stack (HasCallStack)
import ListT (ListT)
import qualified ListT
import SelfInducedGlasses.Random
import System.IO.Unsafe (unsafePerformIO)
import System.Random.Stateful
import Text.PrettyPrint.ANSI.Leijen (Doc, Pretty (..))
import qualified Text.PrettyPrint.ANSI.Leijen as Pretty
import UnliftIO.Async
import UnliftIO.Exception (assert)
import qualified UnliftIO.Foreign as UnliftIO

type ℝ = Float

-- | Dense matrix in row-major order
--
-- We use dense matrices for two purposes:
--
--   * storing the matrix of couplings \(J_{ij}\)
--   * storing bitstrings produced by Monte Carlo (such that each row is a bitstring)
data DenseMatrix v a = DenseMatrix
  { dmNumRows :: {-# UNPACK #-} !Int,
    dmNumCols :: {-# UNPACK #-} !Int,
    dmData :: !(v a)
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
data Configuration
  = Configuration
      {-# UNPACK #-} !Int
      -- ^ Number of spins
      {-# UNPACK #-} !(S.Vector Word64)
      -- ^ Bitstring
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
data MutableConfiguration
  = MutableConfiguration
      {-# UNPACK #-} !Int
      {-# UNPACK #-} !(S.MVector RealWorld Word64)
  deriving stock (Generic)
  deriving anyclass (NFData)

-- | Mutable state that Metropolis-Hastings algorithm keeps around.
data MetropolisState = MetropolisState
  { -- | Current spin configuration
    msConfiguration :: {-# UNPACK #-} !MutableConfiguration,
    -- | Current energy
    msCurrentEnergy :: {-# UNPACK #-} !(IORef Double),
    -- | Potential energy changes upon spin flips
    msDeltaEnergies :: {-# UNPACK #-} !(S.MVector RealWorld ℝ)
  }
  deriving stock (Generic)
  deriving anyclass (NFData)

-- | General information about the sampling that is gathered during a sweep.
--
-- 'SweepStats' obtained from multiple sweeps can be combined using the 'Monoid' instance.
data SweepStats = SweepStats {ssTime :: !Int, ssAcceptProb :: !Double}
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

data ReplicaExchangeState g = ReplicaExchangeState
  { resProb :: !ProbDist,
    resState :: !MetropolisState,
    resColor :: !(Maybe ReplicaColor),
    resStats :: !SweepStats,
    resGen :: !g
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

extractMoments :: Fractional a => (Int, a, a) -> (a, a)
extractMoments (!n, !μ, !m2)
  | n == 0 = (0 / 0, 0 / 0)
  | n == 1 = (μ, 0 / 0)
  | otherwise = (μ, m2 / fromIntegral n)

secondMomentFold ::
  forall m a g.
  (MonadUnliftIO m, Fractional a, U.Unbox a) =>
  (ReplicaExchangeState g -> m a) ->
  FoldM m (ObservableState g) (U.Vector (a, a))
secondMomentFold f = Foldl.FoldM step begin extract
  where
    begin :: m (Maybe (U.MVector RealWorld (Int, a, a)))
    begin = pure Nothing
    go !accs !state !i
      | i < GM.length accs = do
        y <- f . (! i) . unObservableState $ state
        liftIO $ GM.modify accs (\acc -> updateVarianceAcc acc y) i
        go accs state (i + 1)
      | otherwise = pure ()
    step !maybeAccs !state = do
      accs <- case maybeAccs of
        Just accs -> pure accs
        Nothing -> liftIO $ GM.replicate (G.length (unObservableState state)) (0, 0, 0)
      go accs state 0 >> pure (Just accs)
    extract !maybeAccs = do
      case maybeAccs of
        Just accs -> liftIO $ G.map extractMoments <$> (G.unsafeFreeze accs)
        Nothing -> pure G.empty

energyFold ::
  forall m g.
  MonadUnliftIO m =>
  FoldM m (ObservableState g) (U.Vector (Double, Double))
energyFold = secondMomentFold (liftIO . readIORef . msCurrentEnergy . resState)

magnetizationFold ::
  forall m g.
  MonadUnliftIO m =>
  FoldM m (ObservableState g) (U.Vector (Double, Double))
magnetizationFold = secondMomentFold compute
  where
    compute !state = do
      x <- unsafeFreezeConfiguration . msConfiguration . resState $ state
      pure (fromIntegral $ totalMagnetization x)

structureFactorFold ::
  forall m a g.
  (MonadUnliftIO m) =>
  FoldM m (ObservableState g) (B.Vector (DenseMatrix S.Vector ℝ))
structureFactorFold = Foldl.FoldM step begin extract
  where
    begin :: m (Maybe (B.Vector (DenseMatrix (S.MVector RealWorld) ℝ)))
    begin = pure Nothing
    go ::
      B.Vector (DenseMatrix (S.MVector RealWorld) ℝ) ->
      ObservableState g ->
      Int ->
      m ()
    go !accs !state !i
      | i < G.length accs = do
        x <-
          unsafeFreezeConfiguration
            . msConfiguration
            . resState
            . (! i)
            . unObservableState
            $ state
        let (DenseMatrix n _ out) = accs ! i
        withConfiguration x $ \xPtr ->
          withMutableVector out $ \outPtr ->
            liftIO $ c_addSpinSpinCorrelation n xPtr outPtr
        go accs state (i + 1)
      | otherwise = pure ()
    step !maybeAccs !state = do
      accs <- case maybeAccs of
        Just accs -> pure accs
        Nothing ->
          let numReplicas = G.length . unObservableState $ state
              (MutableConfiguration numSpins _) =
                msConfiguration . resState . G.head . unObservableState $ state
           in G.replicateM numReplicas $
                DenseMatrix numSpins numSpins <$> liftIO (GM.new (numSpins * numSpins))
      go accs state 0 >> pure (Just accs)
    extract !maybeAccs = do
      case maybeAccs of
        Just accs -> liftIO $ G.mapM unsafeFreezeDenseMatrix accs
        Nothing -> pure G.empty

-- | Get matrix shape
dmShape :: DenseMatrix v a -> (Int, Int)
dmShape m = (dmNumRows m, dmNumCols m)

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

indexConfiguration :: ConfigurationBatch -> Int -> Configuration
indexConfiguration (ConfigurationBatch numSpins (DenseMatrix numReplicas numWords v)) i =
  Configuration numSpins (G.slice (i * numWords) numWords v)

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
  !r1@(ReplicaExchangeState (ProbDist β1 _) s1 c1 _ _)
  !r2@(ReplicaExchangeState (ProbDist β2 _) s2 c2 _ _)
  !u = do
    !δe <- liftIO $ do
      !e1 <- readIORef (msCurrentEnergy s1)
      !e2 <- readIORef (msCurrentEnergy s2)
      pure $ realToFrac (e2 - e1)
    let !δβ = β2 - β1
    if u <= exp (δβ * δe)
      then
        let -- We swap the spin configurations only, leaving the βs and random number generators
            -- unchanged.
            !r1' = r1 {resState = s2, resColor = c2}
            !r2' = r2 {resState = s1, resColor = c1}
         in pure (r1', r2', True)
      else pure (r1, r2, False)
{-# INLINE maybeExchangeReplicas #-}

runReplicaExchangeSchedule ::
  forall g m.
  (MonadUnliftIO m, NFData g, VectorizedUniformRange g m Int, VectorizedUniformRange g m Float) =>
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
    spawn !replica@(ReplicaExchangeState probDist state _ stats g) !intervals =
      case intervals of
        ((!start, !end) : others) -> do
          !future <- async $ do
            !stats' <- doManySweeps probDist (end - start) sweepSize state g
            pure . force $ replica {resStats = stats <> stats'}
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
        let (future1, intervals1) = accumulators ! i
            (future2, intervals2) = accumulators ! (i + 1)
        state1 <- wait future1
        state2 <- wait future2
        -- TODO: how to update the colors
        (state1', state2', _) <- maybeExchangeReplicas state1 state2 u
        acc1' <- spawn state1' intervals1
        acc2' <- spawn state2' intervals2
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
      | i < numWords - 1 = uniformM g

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
  pure $ ReplicaExchangeState (ProbDist β hamiltonian) state Nothing mempty g

-- | Compute energy of a spin configuration
totalEnergy :: Hamiltonian -> Configuration -> Double
totalEnergy hamiltonian state@(Configuration n _) =
  unsafePerformIO
    $! withDenseMatrix (hInteraction hamiltonian)
    $ \interactionPtr ->
      withVector (hMagneticField hamiltonian) $ \fieldPtr ->
        withConfiguration state $ \statePtr ->
          pure $ realToFrac (c_totalEnergy n interactionPtr fieldPtr statePtr)

-- | Compute the total magnetization of a spin configuration
totalMagnetization :: Configuration -> Int
totalMagnetization state@(Configuration n _) =
  unsafePerformIO
    $! withConfiguration state
    $ \statePtr ->
      pure $ c_totalMagnetization n statePtr

-- | Compute a vector of @ΔE@ that one would get by flipping each spin.
energyChanges :: Hamiltonian -> Configuration -> S.Vector ℝ
energyChanges hamiltonian state@(Configuration n _) =
  G.generate n $ \i ->
    unsafePerformIO
      $! withDenseMatrix (hInteraction hamiltonian)
      $ \couplingsPtr ->
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
  let gs = splitForParallel (G.length βs) gₘₐᵢₙ
  gₘₐᵢₙ' <- thawGen gₘₐᵢₙ
  (gs' :: B.Vector (Xoshiro256PlusPlus (PrimState m))) <- G.mapM thawGen gs
  l <- monteCarloSampling sweepSize numberSweeps hamiltonian βs gₘₐᵢₙ' gs'
  ListT.applyFoldM fold . ListT.take numberMeasure . ListT.drop numberSkip $ l

monteCarloSampling ::
  forall g m.
  (MonadUnliftIO m, NFData g, VectorizedUniformRange g m Int, VectorizedUniformRange g m Float) =>
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
  -- | More random number generators
  B.Vector g ->
  m (ListT m (ObservableState g))
monteCarloSampling sweepSize numberSweeps hamiltonian βs gₘₐᵢₙ gs = do
  s₀ <- initUnfoldState
  pure $ ListT.unfoldM unfoldStep s₀
  where
    !numReplicas = G.length βs
    initUnfoldState =
      ObservableState
        <$> G.zipWithM (\β g -> randomReplicaExchangeStateM hamiltonian β g) βs gs
    unfoldStep (ObservableState replicas) = do
      !schedule <- randomReplicaExchangeScheduleM numReplicas numberSweeps gₘₐᵢₙ
      !replicas' <- ObservableState <$> runReplicaExchangeSchedule sweepSize schedule replicas
      pure $ Just (replicas', replicas')

freezeReplicas :: MonadUnliftIO m => B.Vector (ReplicaExchangeState g) -> m ConfigurationBatch
freezeReplicas replicas
  | G.null replicas = error "replicas array should not be empty"
  | otherwise = do
    let (MutableConfiguration numSpins _) = msConfiguration . resState $ G.head replicas
        numWords = (numSpins + 63) `div` 64
    buf <- liftIO $ SM.new (G.length replicas * numWords)
    withMutableVector buf $ \bufPtr ->
      G.iforM_ replicas $ \i r ->
        withMutableConfiguration (msConfiguration (resState r)) $ \statePtr ->
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

withConfiguration ::
  MonadUnliftIO m =>
  Configuration ->
  (Ptr Word64 -> m a) ->
  m a
withConfiguration (Configuration _ v) = withVector v

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
      t = t1 + t2
      s = fromIntegral t1 * r1 + fromIntegral t2 * r2

instance Monoid SweepStats where
  mempty = SweepStats 0 (0 / 0)

doSweep ::
  MonadUnliftIO m =>
  ProbDist ->
  Int ->
  Ptr Int ->
  Ptr ℝ ->
  MetropolisState ->
  m SweepStats
doSweep (ProbDist β hamiltonian) numberSteps randIntsPtr randFloatsPtr state = do
  -- liftIO $ putStrLn "Calling doSweep..."
  let numberBits = G.length (hMagneticField hamiltonian)
  withDenseMatrix (hInteraction hamiltonian) $ \couplingsPtr ->
    withVector (hMagneticField hamiltonian) $ \fieldPtr ->
      withMutableConfiguration (msConfiguration state) $ \statePtr ->
        withMutableVector (msDeltaEnergies state) $ \deltaEnergiesPtr ->
          liftIO $ do
            x <- unsafeFreezeConfiguration (msConfiguration state)
            -- putStrLn $ "Starting sweep from " <> show x
            -- print =<< peek statePtr
            e₀ <- readIORef (msCurrentEnergy state)
            with e₀ $ \energyPtr -> do
              acceptance <-
                {-# SCC c_runOneSweep #-}
                c_runOneSweep
                  numberBits
                  numberSteps
                  β
                  couplingsPtr
                  fieldPtr
                  randIntsPtr
                  randFloatsPtr
                  statePtr
                  energyPtr
                  deltaEnergiesPtr
              e <- peek energyPtr
              writeIORef (msCurrentEnergy state) e
              e' <- totalEnergy hamiltonian <$> unsafeFreezeConfiguration (msConfiguration state)
              -- print acceptance
              assert (e == e') $
                pure $ SweepStats numberSteps acceptance

doManySweeps ::
  forall g m.
  (MonadUnliftIO m, VectorizedUniformRange g m Float, VectorizedUniformRange g m Int) =>
  ProbDist ->
  Int ->
  Int ->
  MetropolisState ->
  g ->
  m SweepStats
doManySweeps probDist numberSweeps sweepSize state g = do
  randInts <- {-# SCC allocRandInts #-} liftIO $ GM.unsafeNew sweepSize
  randFloats <- {-# SCC allocRandFloats #-} liftIO $ GM.unsafeNew sweepSize
  let !numberSpins = G.length . hMagneticField . pdHamiltonian $ probDist
      go :: SweepStats -> Int -> m SweepStats
      go !stats !i
        | i < numberSweeps = do
          stats' <-
            withMutableVector randInts $ \randIntsPtr ->
              withMutableVector randFloats $ \randFloatsPtr -> do
                {-# SCC fillRandInts #-} vectorizedUniformRM (0, numberSpins) sweepSize randIntsPtr g
                {-# SCC fillRandFloats #-} vectorizedUniformRM (0, 1) sweepSize randFloatsPtr g
                doSweep probDist sweepSize randIntsPtr randFloatsPtr state
          go (stats <> stats') (i + 1)
        | otherwise = pure stats
  go mempty 0

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

foreign import ccall unsafe "metropolis.h add_spin_spin_correlation"
  c_addSpinSpinCorrelation :: Int -> Ptr Word64 -> Ptr Float -> IO ()

foreign import ccall unsafe "metropolis.h run_one_sweep"
  c_runOneSweep :: Int -> Int -> Float -> Ptr Float -> Ptr Float -> Ptr Int -> Ptr Float -> Ptr Word64 -> Ptr Double -> Ptr Float -> IO Double
