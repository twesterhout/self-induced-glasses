module SelfInducedGlasses.Interaction
  ( Model (..),
    Lattice (..),
    LatticeVectors (..),
    Point (..),
    squareLatticeVectors,
    effectiveInteractionDebug,
    dumpCouplingEvolutionToFile,
    buildCouplings,
    buildSKModel,
  )
where

import Control.Monad.Primitive (PrimMonad)
import Control.Monad.ST
import Control.Scheduler
import qualified Data.ByteString.Builder as Builder
import Data.List (foldl')
import qualified Data.Map.Strict as Map
import qualified Data.Set as Set
import Data.Text (Text, unpack)
import Data.Vector.Generic (Vector)
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Storable as S
import SelfInducedGlasses.Core
import System.IO (IOMode (..), withFile)
import qualified System.IO.Unsafe
import qualified System.Random.MWC.Distributions as Distributions
import System.Random.Stateful

-- | Point on a 2-dimensional lattice
data Point a = P !a !a
  deriving stock (Show, Eq, Ord)

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
{-# INLINE scale #-}

data LatticeVectors = LatticeVectors {-# UNPACK #-} !(Point ℝ) {-# UNPACK #-} !(Point ℝ)
  deriving stock (Show, Eq)

data Lattice = Lattice !(Int, Int) !LatticeVectors
  deriving stock (Show, Eq)

data Model = Model
  { modelLattice :: !Lattice,
    modelLambda :: !ℝ
  }
  deriving stock (Show, Eq)

norm :: Floating a => Point a -> a
norm (P x y) = sqrt (x * x + y * y)
{-# INLINE norm #-}

toCartesian :: Integral a => LatticeVectors -> Point a -> Point ℝ
toCartesian (LatticeVectors r₁ r₂) (P a b) = (fromIntegral a) `scale` r₁ + (fromIntegral b) `scale` r₂
{-# INLINE toCartesian #-}

indexToPoint :: Lattice -> Int -> Point Int
indexToPoint (Lattice (width, _) _) i = P x y
  where
    (y, x) = i `divMod` width
{-# INLINE indexToPoint #-}

-- pointToIndex :: Lattice -> Point Int -> Int
-- pointToIndex (Lattice (width, _) _) (P x y) = y * width + x
-- {-# INLINE pointToIndex #-}

squareLatticeVectors :: LatticeVectors
squareLatticeVectors = LatticeVectors (P 1 0) (P 0 1)

rawInteraction :: LatticeVectors -> ℝ -> Point Int -> Point Int -> ℝ
rawInteraction latticeVectors λ pᵢ pⱼ = 1 / (r * r) * sin (2 * pi / λ * r)
  where
    !r = norm $ toCartesian latticeVectors (pᵢ - pⱼ)
{-# INLINE rawInteraction #-}

-- foreign import capi "effective_interaction_loop"
--   effective_interaction_loop :: Cmodel -> CInt -> CInt -> CInt -> CFloat

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

effectiveInteraction :: Model -> Int -> Point Int -> ℝ
effectiveInteraction model rₘₐₓ pⱼ
  | pⱼ == (P 0 0) = 0
  | otherwise = go 0 rₘₐₓ -- trace ("Computing " <> show pⱼ <> " ...") $
  where
    go !acc !r
      | r < 0 = acc
      | otherwise = go (acc + effectiveInteractionLoop model pⱼ r) (r - 1)

effectiveInteractionDebug :: Model -> Point Int -> [(Int, ℝ, ℝ)]
effectiveInteractionDebug model pⱼ = go 0 0
  where
    go !r !acc =
      let !δj = effectiveInteractionLoop model pⱼ r
          !acc' = acc + δj
       in (r, δj, acc') : go (r + 1) acc'

dumpCouplingEvolutionToFile :: Text -> [(Int, ℝ, ℝ)] -> IO ()
dumpCouplingEvolutionToFile filename v =
  withFile (unpack filename) WriteMode $ \h ->
    Builder.hPutBuilder h (renderTable v)
  where
    renderTable rs = mconcat [renderRow r <> Builder.charUtf8 '\n' | r <- rs]
    renderRow (r, δj, j) =
      Builder.intDec r
        <> Builder.charUtf8 ','
        <> Builder.floatDec δj
        <> Builder.charUtf8 ','
        <> Builder.floatDec j

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

requiredPoints ::
  Model ->
  [Point Int]
requiredPoints (Model lattice@(Lattice (width, height) _) _) = go1 (Set.empty :: Set.Set (Point Int)) 0
  where
    !n = width * height
    go1 !acc !i
      | i < n = go2 acc i i
      | otherwise = Set.toAscList acc
    go2 !acc !i !j
      | j < n =
        let !pᵢ = indexToPoint lattice i
            !pⱼ = indexToPoint lattice j
            !δp = projectUsingTranslations lattice (pⱼ - pᵢ)
         in go2 (Set.insert δp acc) i (j + 1)
      | otherwise = go1 acc (i + 1)

evaluateInteraction ::
  Model ->
  Int ->
  [Point Int] ->
  [(Point Int, ℝ)]
evaluateInteraction model rₘₐₓ points =
  System.IO.Unsafe.unsafePerformIO $
    traverseConcurrently Par f points
  where
    f !p = let !e = effectiveInteraction model rₘₐₓ p in pure (p, e)

buildInteractionMatrix ::
  Model ->
  Int ->
  DenseMatrix S.Vector ℝ
buildInteractionMatrix model@(Model lattice@(Lattice (width, height) _) _) rₘₐₓ = buildMatrix (n, n) f
  where
    !n = width * height
    !cache = Map.fromAscList $ evaluateInteraction model rₘₐₓ (requiredPoints model)
    f :: Int -> Int -> ℝ
    f !i !j = case Map.lookup δp cache of
      Just e -> e
      Nothing -> error "this should never happen"
      where
        !pᵢ = indexToPoint lattice i
        !pⱼ = indexToPoint lattice j
        !δp = projectUsingTranslations lattice (pⱼ - pᵢ)

buildCouplings :: Model -> Couplings
buildCouplings model = Couplings $ buildInteractionMatrix model rₘₐₓ
  where
    rₘₐₓ = 1000
{-# SCC buildCouplings #-}

buildSKModel :: (PrimMonad m, StatefulGen g m) => Int -> g -> m Couplings
buildSKModel n g = do
  v <- GM.new (n * n)
  let f !i !j
        | i <= j = realToFrac <$> Distributions.normal 0 1 g
        | otherwise = error "nope, use symmetry"
      go2 !i !j
        | j < n = do
          z <- f i j
          GM.write v (i * n + j) z
          GM.write v (j * n + i) z
          go2 i (j + 1)
        | otherwise = pure ()
      go1 !i
        | i < n = go2 i (i + 1) >> go1 (i + 1)
        | otherwise = pure ()
  go1 0
  Couplings <$> DenseMatrix n n <$> G.unsafeFreeze v
