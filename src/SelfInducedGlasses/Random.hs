{-# LANGUAGE MagicHash #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UnboxedTuples #-}

module SelfInducedGlasses.Random (Xoshiro256PlusPlus, mkXoshiro256PlusPlus, splitForParallel) where

import Control.Monad.Primitive
import Control.Monad.ST
-- import Data.Bits
import Data.Primitive.ByteArray
import qualified Data.Vector as B
import qualified Data.Vector.Generic as G
import Data.Word
import Foreign.Ptr (Ptr, castPtr)
-- import GHC.Prim
-- import GHC.Word (Word64 (..))
import System.Random.Stateful hiding (next)

newtype Xoshiro256PlusPlus s = Xoshiro256PlusPlus (MutableByteArray s)

newtype Xoshiro256PlusPlusState = Xoshiro256PlusPlusState ByteArray

mkXoshiro256PlusPlus :: PrimMonad m => Int -> m (Xoshiro256PlusPlus (PrimState m))
mkXoshiro256PlusPlus seed = do
  let g0 = mkStdGen seed
      (w0, g1) = genWord64 g0
      (w1, g2) = genWord64 g1
      (w2, g3) = genWord64 g2
      (w3, _) = genWord64 g3
  b <- newAlignedPinnedByteArray {-size-} 32 {-alignment-} 64
  writeByteArray b 0 w0
  writeByteArray b 1 w1
  writeByteArray b 2 w2
  writeByteArray b 3 w3
  pure $! Xoshiro256PlusPlus b

foreign import ccall unsafe "xoshiro256plusplus_next"
  nextXoshiro256PlusPlus :: Ptr Word64 -> IO Word64

foreign import ccall unsafe "xoshiro256plusplus_jump"
  jumpXoshiro256PlusPlus :: Ptr Word64 -> IO ()

-- unsafeRotateL :: Word64 -> Int -> Word64
-- unsafeRotateL w k = (unsafeShiftL w k) .|. (unsafeShiftR w (64 - k))
-- {-# INLINE unsafeRotateL #-}

-- nextImpl :: PrimMonad m => MutableByteArray (PrimState m) -> m Word64
-- nextImpl !s = do
--   !s0 <- readByteArray s 0
--   !s1 <- readByteArray s 1
--   !s2 <- readByteArray s 1
--   !s3 <- readByteArray s 3
--   let !result = unsafeRotateL (s0 + s3) 23 + s0
--       !t = s1 `unsafeShiftL` 17
--       !s2' = s0 `xor` s2
--       !s3' = s1 `xor` s3
--       !s1' = s2' `xor` s1
--       !s0' = s3' `xor` s0
--       !s2'' = t `xor` s2'
--       !s3'' = rotateL s3' 45
--   writeByteArray s 0 s0'
--   writeByteArray s 1 s1'
--   writeByteArray s 2 s2''
--   writeByteArray s 3 s3''
--   pure result
-- {-# INLINE nextImpl #-}

-- next :: PrimMonad m => Xoshiro256PlusPlus (PrimState m) -> m Word64
-- next (Xoshiro256PlusPlus !s) = nextImpl s
-- {-# INLINE next #-}

-- nextIO :: Xoshiro256PlusPlus (PrimState IO) -> IO Word64
-- nextIO s = next s -- (Xoshiro256PlusPlus s) = nextXoshiro256PlusPlus (castPtr (mutableByteArrayContents s))

instance (PrimMonad m, PrimState m ~ s) => StatefulGen (Xoshiro256PlusPlus s) m where
  uniformWord64 g = next g
  uniformShortByteString _ _ = undefined

instance (PrimMonad m) => FrozenGen Xoshiro256PlusPlusState m where
  type MutableGen Xoshiro256PlusPlusState m = Xoshiro256PlusPlus (PrimState m)
  freezeGen (Xoshiro256PlusPlus s) = Xoshiro256PlusPlusState <$> freezeByteArray s 0 32
  thawGen (Xoshiro256PlusPlusState s) = Xoshiro256PlusPlus <$> thawByteArray s 0 32

-- instance (PrimMonad m, StatefulGen g m) => StatefulGen g (MetropolisT m) where
--   uniformWord32R r g = lift $ uniformWord32R r g
--   uniformWord64R r g = lift $ uniformWord64R r g
--   uniformWord8 g = lift $ uniformWord8 g
--   uniformWord16 g = lift $ uniformWord16 g
--   uniformWord32 g = lift $ uniformWord32 g
--   uniformWord64 g = lift $ uniformWord64 g
--   uniformShortByteString n g = lift $ uniformShortByteString n g

next :: PrimMonad m => Xoshiro256PlusPlus (PrimState m) -> m Word64
next (Xoshiro256PlusPlus v) =
  unsafeIOToPrim $
    nextXoshiro256PlusPlus (castPtr (mutableByteArrayContents v))
-- primitive $ \s0 ->
-- case nextImpl# v s0 of
--   (# s1, r #) -> (# s1, (W64# r) #)
{-# SCC next #-}

jump :: PrimMonad m => Xoshiro256PlusPlus (PrimState m) -> m ()
jump (Xoshiro256PlusPlus v) =
  unsafeIOToPrim $
    jumpXoshiro256PlusPlus (castPtr (mutableByteArrayContents v))
-- primitive $ \s0 ->
-- case nextImpl# v s0 of
--   (# s1, r #) -> (# s1, (W64# r) #)
{-# SCC jump #-}

jump' :: Xoshiro256PlusPlusState -> Xoshiro256PlusPlusState
jump' (Xoshiro256PlusPlusState v) = runST $ do
  v' <- thawByteArray v 0 (sizeofByteArray v)
  jump (Xoshiro256PlusPlus v')
  Xoshiro256PlusPlusState <$> unsafeFreezeByteArray v'

splitForParallel :: Int -> Xoshiro256PlusPlusState -> B.Vector Xoshiro256PlusPlusState
splitForParallel n s₀ = G.iterateN n jump' s₀

-- unsafeRotateL# :: Word# -> Int# -> Word#
-- unsafeRotateL# w k = (uncheckedShiftL# w k) `or#` (uncheckedShiftRL# w (64# -# k))

-- nextImpl# :: MutableByteArray# d -> State# d -> (# State# d, Word# #)
-- nextImpl# v s0 =
--   case readWordArray# v 0# s0 of
--     (# s1, v0 #) ->
--       case readWordArray# v 1# s1 of
--         (# s2, v1 #) ->
--           case readWordArray# v 2# s2 of
--             (# s3, v2 #) ->
--               case readWordArray# v 3# s3 of
--                 (# s4, v3 #) ->
--                   let result = (unsafeRotateL# (v0 `plusWord#` v3) 23#) `plusWord#` v0
--                       t = v1 `uncheckedShiftL#` 17#
--                       v2' = v0 `xor#` v2
--                       v3' = v1 `xor#` v3
--                       v1' = v2' `xor#` v1
--                       v0' = v3' `xor#` v0
--                       v2'' = t `xor#` v2'
--                       v3'' = unsafeRotateL# v3' 45#
--                    in case writeWordArray# v 0# v0' s4 of
--                         s5 -> case writeWordArray# v 1# v1' s5 of
--                           s6 -> case writeWordArray# v 2# v2'' s6 of
--                             s7 -> case writeWordArray# v 3# v3'' s7 of
--                               s8 -> (# s8, result #)

-- void jump(uint64_t s[4]) {
--   static const uint64_t JUMP[] = {0x180ec6d33cfd0aba, 0xd5a61266f0c9392c,
--                                   0xa9582618e03fc9aa, 0x39abdc4529b1661c};
--
--   uint64_t s0 = 0;
--   uint64_t s1 = 0;
--   uint64_t s2 = 0;
--   uint64_t s3 = 0;
--   for (int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
--     for (int b = 0; b < 64; b++) {
--       if (JUMP[i] & UINT64_C(1) << b) {
--         s0 ^= s[0];
--         s1 ^= s[1];
--         s2 ^= s[2];
--         s3 ^= s[3];
--       }
--       next(s);
--     }
--
--   s[0] = s0;
--   s[1] = s1;
--   s[2] = s2;
--   s[3] = s3;
-- }
