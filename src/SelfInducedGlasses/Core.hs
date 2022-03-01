module SelfInducedGlasses.Core where

import Control.Monad.ST
import Data.List (foldl')
import Data.Vector.Generic (Vector, (!))
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import GHC.Stack (HasCallStack)

-- data System = System

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

energyDifferenceUponFlip :: (Vector v a, Num a) => DenseMatrix v a -> Int -> v a -> a
energyDifferenceUponFlip m i x = 2 * dotProduct (getRow i m) x * (x ! i)

someFunc :: IO ()
someFunc = putStrLn ("someFunc" :: String)

data Point a = P {-# UNPACK #-} !a {-# UNPACK #-} !a
  deriving stock (Show, Eq)

instance Num a => Num (Point a) where
  (+) (P a b) (P c d) = P (a + c) (b + d)
  (-) (P a b) (P c d) = P (a - c) (b - d)
  (*) _ _ = error "Num instance of Point does not implement (*)"
  abs _ = error "Num instance of Point does not implement abs"
  signum _ = error "Num instance of Point does not implement signum"

scale :: Num a => a -> Point a -> Point a
scale c (P a b) = P (c * a) (c * b)

data LatticeVectors = LatticeVectors !(Point Double) !(Point Double)
  deriving stock (Show, Eq)

data Lattice = Lattice !(Int, Int) !LatticeVectors
  deriving stock (Show, Eq)

data Model = Model
  { modelLattice :: !Lattice,
    modelLambda :: !Double
  }

testLattice1 :: Lattice
testLattice1 = Lattice (10, 10) squareLatticeVectors

testModel1 :: Model
testModel1 = Model testLattice1 7

norm :: Floating a => Point a -> a
norm (P x y) = sqrt (x * x + y * y)

toCartesian :: Integral a => LatticeVectors -> Point a -> Point Double
toCartesian (LatticeVectors r₁ r₂) (P a b) = (fromIntegral a) `scale` r₁ + (fromIntegral b) `scale` r₂

squareLatticeVectors :: LatticeVectors
squareLatticeVectors = LatticeVectors (P 1 0) (P 0 1)

rawInteraction :: LatticeVectors -> Double -> Point Int -> Point Int -> Double
rawInteraction latticeVectors λ pᵢ pⱼ = 1 / (r ^ 2) * sin (2 * pi / λ * r)
  where
    !r = norm $ toCartesian latticeVectors (pᵢ - pⱼ)
{-# INLINE rawInteraction #-}

effectiveInteractionLoop :: Model -> Point Int -> Int -> Double
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

effectiveInteraction :: Model -> Int -> Point Int -> Double
effectiveInteraction model rₘₐₓ pⱼ = go 0 rₘₐₓ
  where
    go !acc !r
      | r < 0 = acc
      | otherwise = go (acc + effectiveInteractionLoop model pⱼ r) (r - 1)

effectiveInteractionDebug :: Model -> Point Int -> [(Int, Double, Double)]
effectiveInteractionDebug model@(Model (Lattice _ latticeVectors) λ) pⱼ = go 0 0
  where
    go !r !acc =
      let !δj = effectiveInteractionLoop model pⱼ r
          !acc' = acc + δj
       in (r, δj, acc') : go (r + 1) acc'

buildMatrixCaching ::
  Vector v a =>
  (Int -> Int -> Maybe (Int, Int)) ->
  (Int, Int) ->
  (Int -> Int -> a) ->
  DenseMatrix v a
buildMatrixCaching isCached (rows, cols) f = runST $ do
  v <- GM.new (rows * cols)
  let go2 !i !j
        | j < cols = do
          !x <- case isCached i j of
            Just (i', j')
              | (i', j') < (i, j) -> GM.read v (i' * cols + j')
              | otherwise -> error "entry not yet computed"
            Nothing -> pure (f i j)
          GM.write v (i * cols + j) x
          go2 i (j + 1)
        | otherwise = pure ()
      go1 !i
        | i < rows = go2 i 0 >> go1 (i + 1)
        | otherwise = pure ()
  go1 0
  DenseMatrix rows cols <$> G.unsafeFreeze v

buildMatrix :: Vector v a => (Int, Int) -> (Int -> Int -> a) -> DenseMatrix v a
buildMatrix = buildMatrixCaching (\_ _ -> Nothing)

buildSymmetricMatrix :: Vector v a => (Int, Int) -> (Int -> Int -> a) -> DenseMatrix v a
buildSymmetricMatrix = buildMatrixCaching (\i j -> if i <= j then Nothing else Just (j, i))

closestToOrigin :: Int -> Int -> Int
closestToOrigin n i = if abs i'' <= abs i' then i'' else i'
  where
    i' = i - signum i * (abs i `div` n) * n
    i'' = i' - signum i * n

projectUsingTranslations :: Lattice -> Point Int -> Point Int
projectUsingTranslations (Lattice (width, height) _) (P x y) =
  P (closestToOrigin width x) (closestToOrigin height y)

projectUsingReflection :: Point Int -> Point Int
projectUsingReflection (P x y)
  | x <= y = P x y
  | otherwise = P y x
