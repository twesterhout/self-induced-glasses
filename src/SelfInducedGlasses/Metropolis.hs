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

data MetropolisStats = MetropolisStats
  { msAcceptance :: {-# UNPACK #-} !ℝ,
    msState :: {-# UNPACK #-} !Configuration
  }

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
flipSpin i (MutableConfiguration v) = GM.modify v (* (-1)) i

step :: forall m g. (PrimMonad m, StatefulGen g m) => ℝ -> g -> MetropolisT m Bool
step β g = do
  !x <- asks msConfiguration
  let !numberSpins = case x of (MutableConfiguration v) -> GM.length v
  !i <- lift $ uniformRM (0, numberSpins - 1) g
  !δe <- energyDifferenceUponFlip <$> asks msCoupling <*> pure i <*> unsafeFreeze x
  if δe > 0
    then do
      !u <- lift $ uniformRM (0, 1) g
      if u < exp (-β * δe)
        then flipSpin i x >> pure True
        else pure False
    else flipSpin i x >> pure True

sweep :: forall m g. (PrimMonad m, StatefulGen g m) => ℝ -> Int -> g -> MetropolisT m MetropolisStats
sweep β n g = go 0 n
  where
    go :: Int -> Int -> MetropolisT m MetropolisStats
    go !acc !i
      | i > 0 = do
        accepted <- step β g
        if accepted
          then go (acc + 1) (i - 1)
          else go acc (i - 1)
      | otherwise =
        MetropolisStats (fromIntegral acc / fromIntegral n) <$> (freeze =<< asks msConfiguration)

manySweeps ::
  (PrimMonad m, StatefulGen g m) =>
  ℝ ->
  Int ->
  Int ->
  g ->
  MetropolisT m (DenseMatrix S.Vector Word64, ℝ)
manySweeps β numberSweeps sweepSize g = do
  (DenseMatrix n _ _) <- asks msCoupling
  let numberWords = (n + 63) `div` 64
  buffer <- DenseMatrix numberSweeps numberWords <$> GM.new (numberSweeps * numberWords)
  let go !i !acc
        | i < numberSweeps = do
          (MetropolisStats acc' σ) <- sweep β sweepSize g
          packConfiguration σ (getMutableRow i buffer)
          go (i + 1) (acc + acc')
        | otherwise = pure (acc / fromIntegral numberSweeps)
  acceptance <- go 0 0
  states <- unsafeFreezeDenseMatrix buffer
  pure (states, acceptance)
