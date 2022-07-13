{-# LANGUAGE CApiFFI #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveGeneric #-}

-- |
-- Module      : SelfInducedGlasses.Sampling
-- Description : Monte Carlo sampling for spin glasses
-- Copyright   : (c) Tom Westerhout, 2022
--
-- Here is a longer description of this module, containing some
-- commentary with @some markup@.
module SelfInducedGlasses.Sampling
  ( -- * General types
    DenseMatrix (..),
    Couplings (..),
    Configuration (..),
    MutableConfiguration (..),
    SweepStats (..),
    ProbDist (..),
    ReplicaExchangeSchedule (..),
    ReplicaExchangeState (..),
    maybeExchangeReplicas,
    prettySchedule,
    mkIntervals,
    mkSchedule,
    withCouplings,
    withVector,
    withMutableConfiguration,
    doSweep,
    doManySweeps,
  )
where

import Control.DeepSeq
import Control.Foldl (FoldM (..))
import qualified Control.Foldl as Foldl
import Control.Monad (forM, unless)
import Control.Monad.IO.Unlift
import Control.Monad.Primitive
import Data.IORef
import qualified Data.Vector as B
import Data.Vector.Generic ((!))
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Storable as S
import qualified Data.Vector.Storable.Mutable as SM
import Data.Word
import Foreign.Marshal.Utils
import Foreign.Ptr
import Foreign.Storable
import GHC.Generics
import System.Random.Stateful
import Text.PrettyPrint.ANSI.Leijen (Doc, Pretty (..))
import qualified Text.PrettyPrint.ANSI.Leijen as Pretty
import UnliftIO.Async

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

-- | Matrix of couplings \(J_{ij}\). It is assumed to be symmetric.
newtype Couplings = Couplings {unCouplings :: DenseMatrix S.Vector ℝ}
  deriving stock (Show, Eq, Generic)
  deriving anyclass (NFData)

-- | Bitstring representing a spin configuration
data Configuration
  = Configuration {-# UNPACK #-} !Int {-# UNPACK #-} !(S.Vector Word64)
  deriving stock (Show, Eq, Generic)
  deriving anyclass (NFData)

-- | Mutable counterpart of 'Configuration'.
data MutableConfiguration
  = MutableConfiguration
      {-# UNPACK #-} !Int
      {-# UNPACK #-} !(S.MVector (PrimState IO) Word64)
  deriving stock (Generic)
  deriving anyclass (NFData)

-- | Mutable state that Metropolis-Hastings algorithm keeps around.
data MetropolisState = MetropolisState
  { msConfiguration :: {-# UNPACK #-} !MutableConfiguration,
    msCurrentEnergy :: {-# UNPACK #-} !(IORef Double),
    msDeltaEnergies :: {-# UNPACK #-} !(S.MVector (PrimState IO) ℝ)
  }

-- | General information about the sampling that is gathered during a sweep.
--
-- 'SweepStats' obtained from multiple sweeps can be combined using the 'Monoid' instance.
data SweepStats = SweepStats {ssTime :: !Int, ssAcceptProb :: !Double}
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

-- | Probability distribution \(\exp^{-\beta H}\)
data ProbDist = ProbDist {pdBeta :: !ℝ, pdHamiltonian :: !Couplings}
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

data ReplicaExchangeSchedule = ReplicaExchangeSchedule
  { resIntervals :: B.Vector [(Int, Int)],
    resSwaps :: B.Vector (Int, ℝ)
  }
  deriving stock (Show)

data ReplicaColor = ReplicaBlue | ReplicaRed

data ReplicaExchangeState g = ReplicaExchangeState
  { resProb :: !ProbDist,
    resState :: !MetropolisState,
    resColor :: !ReplicaColor,
    resStats :: !SweepStats,
    resGen :: !g
  }

-- | Try exchanging two replicas at indices @i@ and @i+1@.
maybeExchangeReplicas ::
  MonadUnliftIO m =>
  ReplicaExchangeState g ->
  ReplicaExchangeState g ->
  -- | A uniform random number in @[0, 1)@. Required to make the decision stochastic.
  ℝ ->
  -- | Final state where replicas @i@ and @i+1@ could have been swapped.
  m (ReplicaExchangeState g, ReplicaExchangeState g)
maybeExchangeReplicas
  r1@(ReplicaExchangeState (ProbDist β1 _) s1 c1 _ _)
  r2@(ReplicaExchangeState (ProbDist β2 _) s2 c2 _ _)
  u = do
    δe <- withRunInIO $ \_ -> do
      e1 <- readIORef (msCurrentEnergy s1)
      e2 <- readIORef (msCurrentEnergy s2)
      pure $ realToFrac (e2 - e1)
    if u <= exp (δβ * δe)
      then
        let -- We swap the spin configurations only, leaving the βs and random number generators
            -- unchanged.
            r1' = r1 {resState = s2, resColor = c2}
            r2' = r2 {resState = s1, resColor = c1}
         in pure (r1', r2')
      else pure (r1, r2)
    where
      δβ = β2 - β1

runReplicaExchangeSchedule ::
  forall g m.
  (MonadUnliftIO m, StatefulGen g m) =>
  Int ->
  ReplicaExchangeSchedule ->
  B.Vector (ReplicaExchangeState g) ->
  m (B.Vector (ReplicaExchangeState g))
runReplicaExchangeSchedule sweepSize schedule₀ states₀ =
  Foldl.foldM (FoldM step initial extract) (resSwaps schedule₀)
  where
    spawn ::
      ReplicaExchangeState g ->
      [(Int, Int)] ->
      m (Async (ReplicaExchangeState g), [(Int, Int)])
    spawn replica@(ReplicaExchangeState probDist state _ stats g) intervals =
      case intervals of
        ((start, end) : others) -> do
          future <- async $ do
            stats' <- doManySweeps probDist (end - start) sweepSize state g
            pure $ replica {resStats = stats <> stats'}
          pure (future, others)
        [] -> do
          future <- async $ pure replica
          pure (future, [])
    initial :: m (B.Vector (Async (ReplicaExchangeState g), [(Int, Int)]))
    initial = G.zipWithM spawn states₀ (resIntervals schedule₀)
    step ::
      B.Vector (Async (ReplicaExchangeState g), [(Int, Int)]) ->
      (Int, ℝ) ->
      m (B.Vector (Async (ReplicaExchangeState g), [(Int, Int)]))
    step accumulators (i, u) = do
      let (future1, intervals1) = accumulators ! i
          (future2, intervals2) = accumulators ! (i + 1)
      state1 <- wait future1
      state2 <- wait future2
      -- TODO: how to update the colors
      (state1', state2') <- maybeExchangeReplicas state1 state2 u
      acc1' <- spawn state1' intervals1
      acc2' <- spawn state2' intervals2
      pure $ accumulators G.// [(i, acc1'), (i + 1, acc2')]
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
mkIntervals swaps replicaIdx = go [(0, G.length swaps - 1)] (G.length swaps - 2)
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

-- generateSchedule :: StategulGen g m => Int -> Int ->

--
-- data MutableConfigurationBatch s
--   = MutableConfigurationBatch
--       {-# UNPACK #-} !Int
--       {-# UNPACK #-} !(DenseMatrix (S.MVector s) Word64)
--
-- data ConfigurationBatch
--   = ConfigurationBatch
--       {-# UNPACK #-} !Int
--       {-# UNPACK #-} !(DenseMatrix S.Vector Word64)

withVector ::
  (MonadUnliftIO m, Storable a) =>
  S.Vector a ->
  (Ptr a -> m b) ->
  m b
withVector v f = withRunInIO $ \runInIO -> S.unsafeWith v (runInIO . f)

withMutableVector ::
  (MonadUnliftIO m, Storable a) =>
  S.MVector (PrimState IO) a ->
  (Ptr a -> m b) ->
  m b
withMutableVector v f = withRunInIO $ \runInIO -> SM.unsafeWith v (runInIO . f)

withCouplings :: MonadUnliftIO m => Couplings -> (Ptr ℝ -> m a) -> m a
withCouplings (Couplings (DenseMatrix _ _ v)) = withVector v

withMutableConfiguration ::
  MonadUnliftIO m =>
  MutableConfiguration ->
  (Ptr Word64 -> m a) ->
  m a
withMutableConfiguration (MutableConfiguration _ v) = withMutableVector v

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
  S.Vector Int ->
  S.Vector ℝ ->
  MetropolisState ->
  m SweepStats
doSweep (ProbDist β couplings) randInts randFloats state = do
  let numberSteps = G.length randInts
      numberBits = dmNumRows (unCouplings couplings)
  withCouplings couplings $ \couplingsPtr ->
    withVector randInts $ \randIntsPtr ->
      withVector randFloats $ \randFloatsPtr ->
        withMutableConfiguration (msConfiguration state) $ \statePtr ->
          withMutableVector (msDeltaEnergies state) $ \deltaEnergiesPtr ->
            withRunInIO $ \_ ->
              with (0 :: Double) $ \energyPtr -> do
                acceptance <-
                  c_runOneSweep
                    numberBits
                    numberSteps
                    β
                    couplingsPtr
                    randIntsPtr
                    randFloatsPtr
                    statePtr
                    energyPtr
                    deltaEnergiesPtr
                writeIORef (msCurrentEnergy state) =<< peek energyPtr
                pure $ SweepStats numberSteps acceptance

doManySweeps ::
  forall g m.
  (MonadUnliftIO m, StatefulGen g m) =>
  ProbDist ->
  Int ->
  Int ->
  MetropolisState ->
  g ->
  m SweepStats
doManySweeps probDist numberSweeps sweepSize state g = go mempty 0
  where
    !numberSpins = dmNumRows . unCouplings . pdHamiltonian $ probDist
    go :: SweepStats -> Int -> m SweepStats
    go !stats !i
      | i < numberSweeps = do
        randInts <- G.replicateM sweepSize (uniformRM (0, numberSpins - 1) g)
        randFloats <- G.replicateM sweepSize (uniformRM (0, 1) g)
        stats' <- doSweep probDist randInts randFloats state
        go (stats <> stats') (i + 1)
      | otherwise = pure stats

foreign import capi unsafe "metropolis.h energy_change_upon_flip"
  c_energyChangeUponFlip :: Int -> Ptr Float -> Ptr Word64 -> Int -> Float

foreign import capi unsafe "metropolis.h run_one_sweep"
  c_runOneSweep :: Int -> Int -> Float -> Ptr Float -> Ptr Int -> Ptr Float -> Ptr Word64 -> Ptr Double -> Ptr Float -> IO Double
