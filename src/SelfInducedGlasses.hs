-- |
-- Copyright: (c) 2022 Tom Westerhout
-- SPDX-License-Identifier: BSD-3-Clause
-- Maintainer: Tom Westerhout <14264576+twesterhout@users.noreply.github.com>
--
-- See README for more info
module SelfInducedGlasses
  ( someFunc,
    Point (..),
    reSumBoundaries1D,
    effectiveInteractionDebug,
    effectiveInteraction,
    testLattice1,
    testModel1,
  )
where

import Data.List (foldl')

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
effectiveInteraction model pⱼ rₘₐₓ = go 0 rₘₐₓ
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
