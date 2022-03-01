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

type ℝ = Float

data MetropolisState s = MetropolisState
  { msCoupling :: {-# UNPACK #-} !(DenseMatrix S.Vector ℝ),
    msBeta :: {-# UNPACK #-} !ℝ,
    msSpins :: {-# UNPACK #-} !(S.MVector s ℝ)
  }

data MetropolisStats = MetropolisStats
  { msAcceptance :: {-# UNPACK #-} !ℝ,
    msState :: {-# UNPACK #-} !(S.Vector ℝ)
  }

newtype MetropolisT g m a = MetropolisT
  { unMetropolis :: StateT g (ReaderT (MetropolisState (PrimState m)) m) a
  }
  deriving newtype (Functor, Applicative, Monad, PrimMonad)

deriving newtype instance (PrimMonad m, s ~ PrimState m) => MonadReader (MetropolisState s) (MetropolisT g m)

deriving newtype instance Monad m => MonadState g (MetropolisT g m)

flipSpin :: (PrimMonad m, GM.MVector v a, Num a) => Int -> v (PrimState m) a -> m ()
flipSpin i s = GM.modify s (* (-1)) i

step :: forall m g. (PrimMonad m, RandomGen g) => MetropolisT g m Bool
step = do
  !x <- asks msSpins
  let !numberSpins = GM.length x
  !i <- uniformRM (0, numberSpins - 1) (StateGenM :: StateGenM g)
  !δe <-
    energyDifferenceUponFlip
      <$> asks msCoupling
      <*> pure i
      <*> S.unsafeFreeze x
  if (δe > 0)
    then do
      !u <- uniformRM (0, 1) (StateGenM :: StateGenM g)
      !β <- asks msBeta
      if u < exp (-β * δe)
        then flipSpin i x >> pure True
        else pure False
    else flipSpin i x >> pure True

sweep :: forall m g. (PrimMonad m, RandomGen g) => Int -> MetropolisT g m MetropolisStats
sweep n = go 0 n
  where
    go :: Int -> Int -> MetropolisT g m MetropolisStats
    go !acc !i
      | i > 0 = do
        accepted <- step
        if accepted
          then go (acc + 1) (i - 1)
          else go acc (i - 1)
      | otherwise =
        MetropolisStats (fromIntegral acc / fromIntegral n) <$> (S.freeze =<< asks msSpins)
