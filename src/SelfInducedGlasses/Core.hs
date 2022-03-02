{-# LANGUAGE DeriveAnyClass #-}

module SelfInducedGlasses.Core where

import Control.Monad (forM_)
import Control.Monad.Primitive
import Control.Monad.ST
import Data.Bits
import Data.List (foldl', intercalate)
import Data.MemoTrie
import Data.Vector.Generic (Vector, (!))
import qualified Data.Vector.Generic as G
import Data.Vector.Generic.Mutable (MVector)
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Storable as S
import Data.Word
import Debug.Trace
import GHC.Stack (HasCallStack)
import System.IO

-- data System = System
type ℝ = Float

data DenseMatrix v a = DenseMatrix {-# UNPACK #-} !Int {-# UNPACK #-} !Int {-# UNPACK #-} !(v a)
  deriving (Show, Eq)

getRow :: (HasCallStack, Vector v a) => Int -> DenseMatrix v a -> v a
getRow i (DenseMatrix nRows nCols v)
  | i < nRows = G.slice (i * nCols) nCols v
  | otherwise = error "index out of bounds"

dotProduct :: (Vector v a, Num a) => v a -> v a -> a
dotProduct a b = G.sum $ G.zipWith (*) a b

totalEnergy :: (HasCallStack, Vector v a, Num a) => DenseMatrix v a -> v a -> a
totalEnergy m@(DenseMatrix nRows nCols v) x = -dotProduct x matvec
  where
    matvec = G.generate nRows (\i -> dotProduct (getRow i m) x)

energyDifferenceUponFlip :: (HasCallStack, Vector v ℝ) => DenseMatrix v ℝ -> Int -> Configuration -> ℝ
energyDifferenceUponFlip m i (Configuration x) = 2 * dotProduct (getRow i m) (G.convert x) * (x ! i)

someFunc :: IO ()
someFunc = putStrLn ("someFunc" :: String)

data Point a = P {-# UNPACK #-} !a {-# UNPACK #-} !a
  deriving stock (Show, Eq)

instance Num a => Num (Point a) where
  (+) (P a b) (P c d) = P (a + c) (b + d)
  (-) (P a b) (P c d) = P (a - c) (b - d)
  (*) _ _ = error "Num instance of Point does not implement (*)"
  negate (P a b) = P (-a) (-b)
  abs _ = error "Num instance of Point does not implement abs"
  signum _ = error "Num instance of Point does not implement signum"
  fromInteger _ = error "Num instance of Point does not implement fromInteger"

scale :: Num a => a -> Point a -> Point a
scale c (P a b) = P (c * a) (c * b)

data LatticeVectors = LatticeVectors !(Point ℝ) !(Point ℝ)
  deriving stock (Show, Eq)

data Lattice = Lattice !(Int, Int) !LatticeVectors
  deriving stock (Show, Eq)

data Model = Model
  { modelLattice :: !Lattice,
    modelLambda :: !ℝ
  }

testLattice1 :: Lattice
testLattice1 = Lattice (10, 10) squareLatticeVectors

testModel1 :: Model
testModel1 = Model testLattice1 7

norm :: Floating a => Point a -> a
norm (P x y) = sqrt (x * x + y * y)

toCartesian :: Integral a => LatticeVectors -> Point a -> Point ℝ
toCartesian (LatticeVectors r₁ r₂) (P a b) = (fromIntegral a) `scale` r₁ + (fromIntegral b) `scale` r₂

indexToPoint :: Lattice -> Int -> Point Int
indexToPoint (Lattice (width, _) _) i = P x y
  where
    x = i `mod` width
    y = i `div` width

pointToIndex :: Lattice -> Point Int -> Int
pointToIndex (Lattice (width, _) _) (P x y) = y * width + x

squareLatticeVectors :: LatticeVectors
squareLatticeVectors = LatticeVectors (P 1 0) (P 0 1)

rawInteraction :: LatticeVectors -> ℝ -> Point Int -> Point Int -> ℝ
rawInteraction latticeVectors λ pᵢ pⱼ = 1 / (r ^ 2) * sin (2 * pi / λ * r)
  where
    !r = norm $ toCartesian latticeVectors (pᵢ - pⱼ)
{-# INLINE rawInteraction #-}

effectiveInteractionLoop :: Model -> Point Int -> Int -> ℝ
effectiveInteractionLoop (Model (Lattice (!width, !height) latticeVectors) λ) !pⱼ !radius
  | radius > 0 = top + right + bottom + left
  | radius == 0 = rawInteraction latticeVectors λ (P 0 0) pⱼ
  | otherwise = error "invalid radius"
  where
    combine !acc !pᵢ = acc + rawInteraction latticeVectors λ pᵢ pⱼ
    rescale !x !y = P (x * width) (y * height)
    partialSum = foldl' combine 0.0
    !top = partialSum [rescale x radius | x <- [-radius .. (radius - 1)]]
    !right = partialSum [rescale radius y | y <- [(-radius + 1) .. radius]]
    !bottom = partialSum [rescale x (-radius) | x <- [(-radius + 1) .. radius]]
    !left = partialSum [rescale (-radius) y | y <- [-radius .. radius - 1]]
{-# INLINEABLE effectiveInteractionLoop #-}

effectiveInteraction :: Model -> Int -> Point Int -> ℝ
effectiveInteraction model rₘₐₓ pⱼ
  | pⱼ == (P 0 0) = 0
  | otherwise = trace ("Computing " <> show pⱼ <> " ...") $ go 0 rₘₐₓ
  where
    go !acc !r
      | r < 0 = acc
      | otherwise = go (acc + effectiveInteractionLoop model pⱼ r) (r - 1)

effectiveInteractionDebug :: Model -> Point Int -> [(Int, ℝ, ℝ)]
effectiveInteractionDebug model@(Model (Lattice _ latticeVectors) λ) pⱼ = go 0 0
  where
    go !r !acc =
      let !δj = effectiveInteractionLoop model pⱼ r
          !acc' = acc + δj
       in (r, δj, acc') : go (r + 1) acc'

buildMatrix ::
  Vector v a =>
  (Int, Int) ->
  (Int -> Int -> a) ->
  DenseMatrix v a
buildMatrix (rows, cols) f = runST $ do
  v <- GM.new (rows * cols)
  let go2 !i !j
        | j < cols = do
          let !z = f i j
          GM.write v (i * cols + j) z >> go2 i (j + 1)
        | otherwise = pure ()
      go1 !i
        | i < rows = go2 i 0 >> go1 (i + 1)
        | otherwise = pure ()
  go1 0
  DenseMatrix rows cols <$> G.unsafeFreeze v

closestToOrigin :: Int -> Int -> Int
closestToOrigin n i = if abs i'' <= abs i' then i'' else i'
  where
    i' = i - signum i * (abs i `div` n) * n
    i'' = i' - signum i * n

projectUsingTranslations :: Lattice -> Point Int -> Point Int
projectUsingTranslations (Lattice (width, height) _) (P x y)
  | y' < 0 = -p
  | y' == 0 && x' < 0 = -p
  | otherwise = p
  where
    p@(P x' y') = P (closestToOrigin width x) (closestToOrigin height y)

data MemoTree a = MemoTree (MemoTree a) a (MemoTree a)

instance Functor MemoTree where
  fmap f (MemoTree l m r) = MemoTree (fmap f l) (f m) (fmap f r)

buildInteractionMatrix ::
  forall a v.
  Vector v a =>
  Lattice ->
  (Point Int -> a) ->
  DenseMatrix v a
buildInteractionMatrix lattice@(Lattice (width, height) _) f =
  buildMatrix (width * height, width * height) f'
  where
    f' :: Int -> Int -> a
    f' !i !j = compute (x, y)
      where
        pᵢ = indexToPoint lattice i
        pⱼ = indexToPoint lattice j
        (P x y) = projectUsingTranslations lattice (pⱼ - pᵢ)
    compute :: (Int, Int) -> a
    compute = memo $ \(x, y) -> f (P x y)

newtype Configuration = Configuration (S.Vector ℝ)
  deriving stock (Show, Eq)

newtype MutableConfiguration s = MutableConfiguration (S.MVector s ℝ)

loadWord :: Vector v ℝ => Int -> Int -> v ℝ -> Word64
loadWord offset n v = go 0 0
  where
    go :: Word64 -> Int -> Word64
    go !acc !i
      | i < n =
        let x = G.unsafeIndex v (offset + i)
            acc' =
              if x == 1
                then acc .|. ((1 :: Word64) `shiftL` i)
                else acc
         in go acc' (i + 1)
      | otherwise = acc

storeWord :: forall m v. (PrimMonad m, MVector v ℝ) => Int -> v (PrimState m) ℝ -> Int -> Word64 -> m ()
storeWord n v offset w = go w 0
  where
    go :: Word64 -> Int -> m ()
    go !acc !i
      | i < n = do
        let x = case acc .&. 1 of
              0 -> -1
              1 -> 1
              _ -> error "this should never happen"
        GM.unsafeWrite v (offset + i) x
        go (acc `shiftR` 1) (i + 1)
      | otherwise = pure ()

packConfiguration :: Configuration -> S.Vector Word64
packConfiguration (Configuration v) = runST $ do
  let numberBits = G.length v
      numberWords = (numberBits + 63) `div` 64
  buffer <- GM.new numberWords
  let go !i
        | i + 64 <= numberBits = do
          GM.write buffer (i `div` 64) (loadWord i 64 v)
          go (i + 64)
        | i < numberBits = do
          GM.write buffer (i `div` 64) (loadWord i (numberBits - i) v)
          pure ()
        | otherwise = pure ()
  go 0
  G.unsafeFreeze buffer

unpackConfiguration :: Int -> S.Vector Word64 -> Configuration
unpackConfiguration numberBits v = runST $ do
  buffer <- GM.new numberBits
  let go !i
        | i + 64 <= numberBits = do
          storeWord 64 buffer i (v ! (i `div` 64))
          go (i + 64)
        | i < numberBits = do
          storeWord (numberBits - i) buffer i (v ! (i `div` 64))
          pure ()
        | otherwise = pure ()
  go 0
  Configuration <$> G.unsafeFreeze buffer

thaw :: PrimMonad m => Configuration -> m (MutableConfiguration (PrimState m))
thaw (Configuration v) = MutableConfiguration <$> G.thaw v

freeze :: PrimMonad m => MutableConfiguration (PrimState m) -> m Configuration
freeze (MutableConfiguration v) = Configuration <$> G.freeze v

unsafeFreeze :: PrimMonad m => MutableConfiguration (PrimState m) -> m Configuration
unsafeFreeze (MutableConfiguration v) = Configuration <$> G.unsafeFreeze v

computeMagnetizationPerSite :: Configuration -> ℝ
computeMagnetizationPerSite (Configuration v) = G.sum v / fromIntegral (G.length v)

computeEnergyPerSite :: G.Vector v ℝ => DenseMatrix v ℝ -> Configuration -> ℝ
computeEnergyPerSite couplings (Configuration v) =
  totalEnergy couplings (G.convert v) / fromIntegral (G.length v)

saveForGnuplot :: Lattice -> FilePath -> Configuration -> IO ()
saveForGnuplot (Lattice (width, height) _) filepath (Configuration v) = do
  let m = DenseMatrix height width v
  withFile filepath WriteMode $ \h ->
    forM_ [0 .. height - 1] $ \i ->
      hPutStrLn h
        . intercalate "\t"
        . fmap show
        . G.toList
        . getRow (height - 1 - i)
        $ m
