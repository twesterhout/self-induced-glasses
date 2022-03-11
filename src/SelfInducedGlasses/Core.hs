{-# LANGUAGE CApiFFI #-}
{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE TypeFamilies #-}

module SelfInducedGlasses.Core where

import Control.DeepSeq
-- import Control.Monad (forM_)
import Control.Monad.Primitive
-- import Control.Monad.ST
-- import Data.Bits
-- import qualified Data.ByteString.Builder as Builder
import Data.Coerce (coerce)
import qualified Data.HDF5 as H5
-- import Data.List (foldl', intercalate)
-- import Data.MemoTrie
-- import Data.Text (Text, pack, unpack)
import Data.Vector.Generic (Mutable, Vector)
import qualified Data.Vector.Generic as G
import Data.Vector.Generic.Mutable (MVector)
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Storable as S
import qualified Data.Vector.Storable.Mutable as SM
import Data.Word
-- import Debug.Trace
import Foreign.C.Types (CFloat (..), CPtrdiff (..))
import Foreign.ForeignPtr (withForeignPtr)
import Foreign.Marshal.Utils (copyBytes)
import Foreign.Ptr (Ptr, castPtr, plusPtr)
import Foreign.Storable (Storable, sizeOf)
import GHC.Generics (Generic)
import GHC.Stack (HasCallStack)
-- import System.IO
import qualified System.IO.Unsafe

type ℝ = Float

-- | Dense matrix in row-major order
--
-- We use dense matrices for two purposes:
--
--   * storing the matrix of couplings Jᵢⱼ
--   * storing bitstrings produced by Monte Carlo (such that each row is a bitstring)
data DenseMatrix v a = DenseMatrix {-# UNPACK #-} !Int {-# UNPACK #-} !Int !(v a)
  deriving stock (Show, Eq, Generic)
  deriving anyclass (NFData)

type instance H5.ElementOf (DenseMatrix v a) = a

instance (H5.KnownDatatype a, Storable a) => H5.KnownDataset' (DenseMatrix S.Vector a) where
  withArrayView' (DenseMatrix r c v) action = action (H5.ArrayView' fp [r, c] [c, 1])
    where
      (fp, _) = S.unsafeToForeignPtr0 v
  fromArrayView' (H5.ArrayView' fp [r, c] [c', 1])
    | c == c' = pure $ DenseMatrix r c (S.unsafeFromForeignPtr0 fp (r * c))
  fromArrayView' _ = error "invalid shapes or strides"

sliceDenseMatrix :: (HasCallStack, Vector v a) => Int -> Int -> Int -> DenseMatrix v a -> DenseMatrix v a
sliceDenseMatrix 0 offset size (DenseMatrix r c v) = DenseMatrix size' c (G.slice (offset' * c) (size' * c) v)
  where
    offset'
      | offset >= 0 = min offset r
      | otherwise = error "invalid offset"
    size'
      | size >= 0 = min size (r - offset')
      | size == -1 = r - offset'
      | otherwise = error "invalid size"
sliceDenseMatrix 1 _ _ _ = error "not implemented"
sliceDenseMatrix _ _ _ _ = error "invalid dimension"

-- | Interpret the matrix as a list of rows
denseMatrixRows :: Vector v a => DenseMatrix v a -> [v a]
denseMatrixRows (DenseMatrix r c v) = go 0
  where
    go !i
      | i < r = G.slice (i * c) c v : go (i + 1)
      | otherwise = []
{-# SCC denseMatrixRows #-}

unsafeFreezeDenseMatrix :: (PrimMonad m, Vector v a) => DenseMatrix (Mutable v (PrimState m)) a -> m (DenseMatrix v a)
unsafeFreezeDenseMatrix (DenseMatrix r c v) = DenseMatrix r c <$> G.unsafeFreeze v
{-# INLINE unsafeFreezeDenseMatrix #-}

newtype Couplings = Couplings (DenseMatrix S.Vector ℝ)

getRow :: (HasCallStack, Vector v a) => Int -> DenseMatrix v a -> v a
getRow i (DenseMatrix nRows nCols v)
  | i < nRows = G.slice (i * nCols) nCols v
  | otherwise = error "index out of bounds"

getMutableRow :: (HasCallStack, MVector v a) => Int -> DenseMatrix (v s) a -> v s a
getMutableRow i (DenseMatrix nRows nCols v)
  | i < nRows = GM.slice (i * nCols) nCols v
  | otherwise = error "index out of bounds"

-- dotProduct :: (Vector v a, Num a) => v a -> v a -> a
-- dotProduct a b = G.sum $ G.zipWith (*) a b

foreign import capi unsafe "helpers.h total_energy"
  total_energy :: CPtrdiff -> Ptr CFloat -> Ptr Word64 -> IO CFloat

totalEnergy :: Couplings -> Configuration -> ℝ
totalEnergy (Couplings (DenseMatrix _ _ v)) (Configuration n x) =
  (coerce :: CFloat -> Float) . System.IO.Unsafe.unsafePerformIO $
    S.unsafeWith v $ \(couplingsPtr :: Ptr Float) ->
      S.unsafeWith x $ \xPtr ->
        total_energy (fromIntegral n) (castPtr couplingsPtr) xPtr
{-# SCC totalEnergy #-}

foreign import ccall unsafe "helpers.h energy_change_upon_flip"
  energy_change_upon_flip :: Int -> Ptr Float -> Ptr Word64 -> Int -> IO Float

energyChangeUponFlip :: PrimMonad m => Couplings -> Int -> MutableConfiguration (PrimState m) -> m ℝ
energyChangeUponFlip (Couplings (DenseMatrix _ _ v)) i (MutableConfiguration n (SM.MVector _ x)) =
  unsafeIOToPrim $
    S.unsafeWith v $ \couplingsPtr ->
      withForeignPtr x $ \xPtr ->
        energy_change_upon_flip n couplingsPtr xPtr i
{-# INLINE energyChangeUponFlip #-}

-- {-# SCC energyChangeUponFlip #-}

data Configuration
  = Configuration
      {-# UNPACK #-} !Int
      {-# UNPACK #-} !(S.Vector Word64)
  deriving stock (Eq)

data MutableConfiguration s
  = MutableConfiguration
      {-# UNPACK #-} !Int
      {-# UNPACK #-} !(S.MVector s Word64)

data MutableConfigurationBatch s
  = MutableConfigurationBatch
      {-# UNPACK #-} !Int
      {-# UNPACK #-} !(DenseMatrix (S.MVector s) Word64)

data ConfigurationBatch
  = ConfigurationBatch
      {-# UNPACK #-} !Int
      {-# UNPACK #-} !(DenseMatrix S.Vector Word64)

insertConfiguration :: PrimMonad m => MutableConfigurationBatch (PrimState m) -> Int -> Configuration -> m ()
insertConfiguration
  (MutableConfigurationBatch _ (DenseMatrix _ numberWords (SM.MVector _ fp)))
  i
  (Configuration _ x) =
    unsafeIOToPrim $
      withForeignPtr fp $ \bufferPtr ->
        S.unsafeWith x $ \src ->
          let bytesWidth = numberWords * sizeOf (undefined :: Word64)
              dest = bufferPtr `plusPtr` (i * bytesWidth)
           in copyBytes dest src bytesWidth

batchToList :: ConfigurationBatch -> [Configuration]
batchToList (ConfigurationBatch n matrix) = Configuration n <$> denseMatrixRows matrix

unsafeFreezeBatch :: PrimMonad m => MutableConfigurationBatch (PrimState m) -> m ConfigurationBatch
unsafeFreezeBatch (MutableConfigurationBatch n matrix) = ConfigurationBatch n <$> unsafeFreezeDenseMatrix matrix

batchChunksOf :: Int -> ConfigurationBatch -> [ConfigurationBatch]
batchChunksOf chunkSize (ConfigurationBatch n matrix@(DenseMatrix totalSize _ _)) = go 0
  where
    getChunk !i = ConfigurationBatch n (sliceDenseMatrix 0 i chunkSize matrix)
    go !i
      | i + chunkSize <= totalSize = getChunk i : go (i + chunkSize)
      | i < totalSize = [getChunk i]
      | otherwise = []

-- newtype Configuration = Configuration (S.Vector ℝ)
--   deriving stock (Show, Eq)
--
-- newtype MutableConfiguration s = MutableConfiguration (S.MVector s ℝ)

-- loadWord :: Vector v ℝ => Int -> Int -> v ℝ -> Word64
-- loadWord offset n v = go 0 0
--   where
--     go :: Word64 -> Int -> Word64
--     go !acc !i
--       | i < n =
--         let x = G.unsafeIndex v (offset + i)
--             acc' =
--               if x == 1
--                 then acc .|. ((1 :: Word64) `unsafeShiftL` i)
--                 else acc
--          in go acc' (i + 1)
--       | otherwise = acc

-- storeWord :: forall m v. (PrimMonad m, MVector v ℝ) => Int -> v (PrimState m) ℝ -> Int -> Word64 -> m ()
-- storeWord n v offset w = go w 0
--   where
--     go :: Word64 -> Int -> m ()
--     go !acc !i
--       | i < n = do
--         let x = case acc .&. 1 of
--               0 -> -1
--               1 -> 1
--               _ -> error "this should never happen"
--         GM.unsafeWrite v (offset + i) x
--         go (acc `shiftR` 1) (i + 1)
--       | otherwise = pure ()

-- packConfiguration :: (PrimMonad m, MVector v Word64) => Configuration -> v (PrimState m) Word64 -> m ()
-- packConfiguration (Configuration v) buffer = do
--   let numberBits = G.length v
--       numberWords = (numberBits + 63) `div` 64
--   let go !i
--         | i + 64 <= numberBits = do
--           GM.unsafeWrite buffer (i `div` 64) (loadWord i 64 v)
--           go (i + 64)
--         | i < numberBits = do
--           GM.unsafeWrite buffer (i `div` 64) (loadWord i (numberBits - i) v)
--           pure ()
--         | otherwise = pure ()
--   go 0
-- {-# SCC packConfiguration #-}

-- unpackConfiguration :: Int -> S.Vector Word64 -> Configuration
-- unpackConfiguration numberBits v = runST $ do
--   buffer <- GM.new numberBits
--   let go !i
--         | i + 64 <= numberBits = do
--           storeWord 64 buffer i (v ! (i `div` 64))
--           go (i + 64)
--         | i < numberBits = do
--           storeWord (numberBits - i) buffer i (v ! (i `div` 64))
--           pure ()
--         | otherwise = pure ()
--   go 0
--   Configuration <$> G.unsafeFreeze buffer

thawConfiguration :: PrimMonad m => Configuration -> m (MutableConfiguration (PrimState m))
thawConfiguration (Configuration n v) = MutableConfiguration n <$> G.thaw v
{-# INLINE thawConfiguration #-}

freezeConfiguration :: PrimMonad m => MutableConfiguration (PrimState m) -> m Configuration
freezeConfiguration (MutableConfiguration n v) = Configuration n <$> G.freeze v
{-# INLINE freezeConfiguration #-}

unsafeFreezeConfiguration :: PrimMonad m => MutableConfiguration (PrimState m) -> m Configuration
unsafeFreezeConfiguration (MutableConfiguration n v) = Configuration n <$> G.unsafeFreeze v
{-# INLINE unsafeFreezeConfiguration #-}
