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
import Foreign.Storable (Storable)
import SelfInducedGlasses.Core
import System.Random.Stateful

data SamplingOptions g = SamplingOptions
  { soCoupling :: {-# UNPACK #-} !(DenseMatrix S.Vector ℝ),
    -- soBeta :: {-# UNPACK #-} !ℝ,
    soGenerator :: {-# UNPACK #-} !g,
    soInitialConfiguration :: !(Maybe Configuration)
  }

data MetropolisState g s = MetropolisState
  { msOptions :: !(SamplingOptions g),
    msConfiguration :: {-# UNPACK #-} !(MutableConfiguration s)
  }

data MetropolisStats = MetropolisStats
  { msAcceptance :: {-# UNPACK #-} !ℝ,
    msState :: {-# UNPACK #-} !Configuration
  }

newtype MetropolisT g m a = MetropolisT
  { unMetropolis :: ReaderT (MetropolisState g (PrimState m)) m a
  }
  deriving newtype (Functor, Applicative, Monad, PrimMonad)

instance MonadTrans (MetropolisT g) where
  lift f = MetropolisT (lift f)

deriving newtype instance (PrimMonad m, s ~ PrimState m) => MonadReader (MetropolisState g s) (MetropolisT g m)

deriving newtype instance (PrimMonad m, MonadState g m) => MonadState g (MetropolisT (StateGenM g) m)

-- instance (PrimMonad m, s ~ PrimState m, StatefulGen g m) => StatefulGen g (MetropolisT g m) where
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

prepareInitialState :: (PrimMonad m, StatefulGen g m) => SamplingOptions g -> m (MetropolisState g (PrimState m))
prepareInitialState options = do
  v₀ <- case soInitialConfiguration options of
    Nothing ->
      let (DenseMatrix n _ _) = soCoupling options
       in randomConfigurationM n (soGenerator options)
    Just v -> pure v
  MetropolisState options <$> thaw v₀

runMetropolisT :: (PrimMonad m, StatefulGen g m) => SamplingOptions g -> MetropolisT g m a -> m a
runMetropolisT options m = do
  s <- prepareInitialState options
  runReaderT (unMetropolis m) s

flipSpin :: PrimMonad m => Int -> MutableConfiguration (PrimState m) -> m ()
flipSpin i (MutableConfiguration v) = GM.modify v (* (-1)) i

step :: forall m g. (PrimMonad m, StatefulGen g (MetropolisT g m)) => ℝ -> MetropolisT g m Bool
step β = do
  !x <- asks msConfiguration
  let !numberSpins = case x of (MutableConfiguration v) -> GM.length v
  !i <- uniformRM (0, numberSpins - 1) =<< asks (soGenerator . msOptions)
  !δe <-
    energyDifferenceUponFlip
      <$> asks (soCoupling . msOptions)
      <*> pure i
      <*> unsafeFreeze x
  if (δe > 0)
    then do
      !u <- uniformRM (0, 1) =<< asks (soGenerator . msOptions)
      if u < exp (-β * δe)
        then flipSpin i x >> pure True
        else pure False
    else flipSpin i x >> pure True

sweep :: forall m g. (PrimMonad m, StatefulGen g (MetropolisT g m)) => ℝ -> Int -> MetropolisT g m MetropolisStats
sweep β n = go 0 n
  where
    go :: Int -> Int -> MetropolisT g m MetropolisStats
    go !acc !i
      | i > 0 = do
        accepted <- step β
        if accepted
          then go (acc + 1) (i - 1)
          else go acc (i - 1)
      | otherwise =
        MetropolisStats (fromIntegral acc / fromIntegral n) <$> (freeze =<< asks msConfiguration)

manySweeps :: (PrimMonad m, StatefulGen g (MetropolisT g m)) => ℝ -> Int -> Int -> MetropolisT g m [MetropolisStats]
manySweeps β numberSweeps sweepSize = replicateM numberSweeps (sweep β sweepSize)
