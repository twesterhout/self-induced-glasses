{-# LANGUAGE DataKinds #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE MagicHash #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UnboxedTuples #-}
{-# LANGUAGE UndecidableInstances #-}

module SelfInducedGlasses.Metropolis where

import Control.Monad.Primitive
import Control.Monad.Reader
import Control.Monad.State.Strict
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Storable as S
import qualified Data.Vector.Storable.Mutable as SM
import Data.Word
import Foreign.Storable (Storable)
import SelfInducedGlasses.Core
import System.Random.Stateful

data SamplingOptions = SamplingOptions
  { soCoupling :: {-# UNPACK #-} !(DenseMatrix S.Vector ℝ),
    soInitialConfiguration :: !(Maybe Configuration)
  }

data MetropolisState s = MetropolisState
  { msCoupling :: {-# UNPACK #-} !(DenseMatrix S.Vector ℝ),
    msConfiguration :: {-# UNPACK #-} !(MutableConfiguration s),
    msDeltaEnergies :: {-# UNPACK #-} !(S.MVector s ℝ)
  }

-- data MetropolisStats = MetropolisStats
--   { msAcceptance :: {-# UNPACK #-} !ℝ,
--     msState :: {-# UNPACK #-} !Configuration
--   }

newtype MetropolisT m a = MetropolisT
  { unMetropolis :: ReaderT (MetropolisState (PrimState m)) m a
  }
  deriving newtype (Functor, Applicative, Monad, PrimMonad)

instance MonadTrans MetropolisT where
  lift f = MetropolisT (lift f)

deriving newtype instance
  (PrimMonad m, s ~ PrimState m) =>
  MonadReader (MetropolisState s) (MetropolisT m)

-- deriving newtype instance (PrimMonad m, MonadState g m) => MonadState g (MetropolisT m)

-- instance (PrimMonad m, StatefulGen g m) => StatefulGen g (MetropolisT m) where
--   uniformWord32R r g = lift $ uniformWord32R r g
--   uniformWord64R r g = lift $ uniformWord64R r g
--   uniformWord8 g = lift $ uniformWord8 g
--   uniformWord16 g = lift $ uniformWord16 g
--   uniformWord32 g = lift $ uniformWord32 g
--   uniformWord64 g = lift $ uniformWord64 g
--   uniformShortByteString n g = lift $ uniformShortByteString n g

randomConfigurationM :: (Monad m, StatefulGen g m) => Int -> g -> m Configuration
randomConfigurationM n g = Configuration . G.map fromBool <$> G.replicateM n (uniformM g)
  where
    fromBool True = 1
    fromBool False = -1

randomConfiguration :: forall g. RandomGen g => Int -> g -> (Configuration, g)
randomConfiguration n g = runStateGen g (randomConfigurationM n)

prepareInitialState ::
  (PrimMonad m, StatefulGen g m) =>
  SamplingOptions ->
  g ->
  m (MetropolisState (PrimState m))
prepareInitialState options g = do
  let coupling@(DenseMatrix n _ _) = soCoupling options
  v₀ <- case soInitialConfiguration options of
    Nothing -> randomConfigurationM n g
    Just v -> pure v
  MetropolisState coupling <$> thaw v₀ <*> GM.new n

runMetropolisT ::
  (PrimMonad m, StatefulGen g m) =>
  (g -> MetropolisT m a) ->
  SamplingOptions ->
  g ->
  m a
runMetropolisT m options g = do
  s <- prepareInitialState options g
  runReaderT (unMetropolis (m g)) s

flipSpin :: PrimMonad m => Int -> MutableConfiguration (PrimState m) -> m ()
flipSpin i (MutableConfiguration v) = GM.unsafeModify v (* (-1)) i
{-# INLINE flipSpin #-}

numberSpins :: MutableConfiguration s -> Int
numberSpins (MutableConfiguration v) = GM.length v
{-# INLINE numberSpins #-}

step :: forall m g. (PrimMonad m, StatefulGen g m) => ℝ -> g -> MetropolisT m Bool
step !β !g = do
  x <- asks msConfiguration
  i <- lift $! uniformRM (0, numberSpins x - 1) g
  coupling <- asks msCoupling
  σ <- unsafeFreeze x
  let δe = energyDifferenceUponFlip coupling i σ
  if δe > 0
    then do
      u <- lift $! uniformRM (0, 1) g
      if u < exp (-β * δe)
        then flipSpin i x >> pure True
        else pure False
    else flipSpin i x >> pure True
{-# INLINE step #-}

sweep :: forall m g v. (PrimMonad m, StatefulGen g m, GM.MVector v Word64) => ℝ -> Int -> Maybe (v (PrimState m) Word64) -> g -> MetropolisT m ℝ
sweep !β !n !buffer !g = go 0 n
  where
    go :: Int -> Int -> MetropolisT m ℝ
    go !acc !i
      | i > 0 = do
        accepted <- step β g
        if accepted
          then go (acc + 1) (i - 1)
          else go acc (i - 1)
      | otherwise = do
        case buffer of
          Just dest -> do
            σ <- unsafeFreeze =<< asks msConfiguration
            packConfiguration σ dest
          Nothing -> pure ()
        pure $ fromIntegral acc / fromIntegral n
{-# SCC sweep #-}

manySweeps ::
  (PrimMonad m, StatefulGen g m) =>
  ℝ ->
  Int ->
  Int ->
  g ->
  MetropolisT m (DenseMatrix S.Vector Word64, ℝ)
manySweeps !β !numberSweeps !sweepSize !g = do
  (DenseMatrix n _ _) <- asks msCoupling
  let numberWords = (n + 63) `div` 64
  buffer <- DenseMatrix numberSweeps numberWords <$> GM.new (numberSweeps * numberWords)
  let go !i !acc
        | i < numberSweeps = do
          acc' <- sweep β sweepSize (Just (getMutableRow i buffer)) g
          go (i + 1) (acc + acc')
        | otherwise = pure (acc / fromIntegral numberSweeps)
  acceptance <- go 0 0
  states <- unsafeFreezeDenseMatrix buffer
  pure (states, acceptance)
{-# SCC manySweeps #-}

thermalize ::
  (PrimMonad m, StatefulGen g m) =>
  ℝ ->
  Int ->
  Int ->
  g ->
  MetropolisT m ℝ
thermalize β numberSweeps sweepSize g = snd <$> manySweeps β 1 (numberSweeps * sweepSize) g
{-# SCC thermalize #-}

anneal ::
  (PrimMonad m, StatefulGen g m) =>
  [(ℝ, Int, Int)] ->
  g ->
  MetropolisT m [(DenseMatrix S.Vector Word64, ℝ)]
anneal steps g = do
  (DenseMatrix n _ _) <- asks msCoupling
  let sweepSize = n
  forM steps $ \(β, numberThermalization, numberGathering) -> do
    thermalize β numberThermalization sweepSize g
    manySweeps β numberGathering sweepSize g
