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
import Control.Monad.ST
-- import Control.Monad.State.Strict
import Data.Bits
import Data.Primitive.Ptr
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Storable as S
import qualified Data.Vector.Storable.Mutable as SM
import Data.Word
import Foreign.ForeignPtr (withForeignPtr)
-- import Foreign.Storable (Storable)
import SelfInducedGlasses.Core
import SelfInducedGlasses.Random
import System.Random.Stateful

data SamplingOptions = SamplingOptions
  { soCoupling :: {-# UNPACK #-} !Couplings,
    soInitialConfiguration :: !(Maybe Configuration)
  }

data MetropolisState s = MetropolisState
  { msCoupling :: {-# UNPACK #-} !Couplings,
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

randomConfigurationM :: (PrimMonad m, StatefulGen g m) => Int -> g -> m Configuration
randomConfigurationM n g = do
  let numberWords = (n + 63) `div` 64
      rest = n `mod` 64
  buffer <- G.unsafeThaw =<< G.replicateM numberWords (uniformM g)
  when (rest /= 0) $ do
    let mask = ((1 :: Word64) `unsafeShiftL` rest) - 1
    GM.unsafeModify buffer (.&. mask) (numberWords - 1)
  Configuration n <$> G.unsafeFreeze buffer

randomConfiguration :: forall g. RandomGen g => Int -> g -> (Configuration, g)
randomConfiguration n g = runST $ runStateGenT g (randomConfigurationM n)

prepareInitialState ::
  (PrimMonad m, StatefulGen g m) =>
  SamplingOptions ->
  g ->
  m (MetropolisState (PrimState m))
prepareInitialState options g = do
  let coupling@(Couplings (DenseMatrix n _ _)) = soCoupling options
  v₀ <- case soInitialConfiguration options of
    Nothing -> randomConfigurationM n g
    Just v -> pure v
  MetropolisState coupling <$> thawConfiguration v₀ <*> GM.new n

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
flipSpin !i (MutableConfiguration _ (SM.MVector _ fp)) = unsafeIOToPrim $
  withForeignPtr fp $ \p -> do
    w <- readOffPtr p wordIndex
    writeOffPtr p wordIndex (complementBit w bitIndex)
  where
    (!wordIndex, !bitIndex) = i `divMod` 64
-- {-# SCC flipSpin #-}
{-# INLINE flipSpin #-}

numberSpins :: MutableConfiguration s -> Int
numberSpins (MutableConfiguration n _) = n
{-# INLINE numberSpins #-}

step :: forall m g. (PrimMonad m, StatefulGen g m) => ℝ -> g -> MetropolisT m Bool
step !β !g = do
  (MetropolisState !coupling !x _) <- ask
  !i <- lift $ uniformRM (0, numberSpins x - 1) g
  !δe <- energyChangeUponFlip coupling i x
  if δe > 0
    then do
      u <- lift $ uniformRM (0, 1) g
      if u < exp (-β * δe)
        then flipSpin i x >> pure True
        else pure False
    else flipSpin i x >> pure True
{-# SCC step #-}

sweep :: forall m g. (PrimMonad m, StatefulGen g m) => ℝ -> Int -> g -> MetropolisT m ℝ
sweep !β !n !g = go 0 n
  where
    go :: Int -> Int -> MetropolisT m ℝ
    go !acc !i
      | i > 0 = do
        accepted <- step β g
        if accepted
          then go (acc + 1) (i - 1)
          else go acc (i - 1)
      | otherwise =
        let !acceptance = fromIntegral acc / fromIntegral n
         in pure acceptance
{-# SCC sweep #-}

manySweeps ::
  (PrimMonad m, StatefulGen g m) =>
  ℝ ->
  Int ->
  Int ->
  g ->
  MetropolisT m (ConfigurationBatch, ℝ)
manySweeps !β !numberSweeps !sweepSize !g = do
  (Couplings (DenseMatrix n _ _)) <- asks msCoupling
  let numberWords = (n + 63) `div` 64
  buffer <-
    MutableConfigurationBatch n
      <$> DenseMatrix numberSweeps numberWords
      <$> GM.new (numberSweeps * numberWords)
  let go !i !acc
        | i < numberSweeps = do
          acc' <- sweep β sweepSize g
          asks msConfiguration
            >>= unsafeFreezeConfiguration
            >>= insertConfiguration buffer i
          go (i + 1) (acc + acc')
        | otherwise = pure (acc / fromIntegral numberSweeps)
  acceptance <- go 0 0
  states <- unsafeFreezeBatch buffer
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
  Int ->
  [(ℝ, Int, Int)] ->
  g ->
  MetropolisT m [(ConfigurationBatch, ℝ)]
anneal sweepSize steps g = do
  -- (Couplings (DenseMatrix sweepSize _ _)) <- asks msCoupling
  forM steps $ \(β, numberThermalization, numberGathering) -> do
    _ <- thermalize β numberThermalization sweepSize g
    manySweeps β numberGathering sweepSize g

annealIO ::
  Int ->
  Couplings ->
  [(ℝ, Int, Int)] ->
  Int ->
  IO [(ConfigurationBatch, ℝ)]
annealIO sweepSize couplings steps seed = do
  let run g =
        forM steps $ \(β, numberThermalization, numberGathering) -> do
          _ <- thermalize β numberThermalization sweepSize g
          manySweeps β numberGathering sweepSize g
  (g :: Xoshiro256PlusPlus (PrimState IO)) <- mkXoshiro256PlusPlus seed
  runMetropolisT run (SamplingOptions couplings Nothing) g
