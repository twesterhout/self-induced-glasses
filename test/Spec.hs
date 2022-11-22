{-# LANGUAGE BinaryLiterals #-}
{-# LANGUAGE OverloadedLists #-}

module Main (main) where

import Data.Bits (bit)
import qualified Data.Vector.Generic as G
import qualified ListT
import SelfInducedGlasses.Random
import SelfInducedGlasses.Sampling
import System.Random
import Test.Hspec

main :: IO ()
main = hspec $ do
  describe "allConfigurations" $ do
    it "builds all possible spin configurations" $ do
      states <- ListT.toList $ allConfigurations 2
      states
        `shouldBe` fmap (Configuration 2 . G.singleton) [0, 1, 2, 3]
  describe "ising2d" $ do
    it "builds the Hamiltonian" $ do
      let !h2 = ferromagneticIsingModelSquare2D 2 0.5
      hMagneticField h2 `shouldBe` [0.5, 0.5, 0.5, 0.5]
      hInteraction h2
        `shouldBe` [ [0, -0.5, -0.5, 0],
                     [-0.5, 0, 0, -0.5],
                     [-0.5, 0, 0, -0.5],
                     [0, -0.5, -0.5, 0]
                   ]
      let !h3 = ferromagneticIsingModelSquare2D 3 (-0.1)
      hMagneticField h3 `shouldSatisfy` G.all (== (-0.1))
      hInteraction h3
        `shouldBe` [ [0, -0.5, -0.5, -0.5, 0, 0, -0.5, 0, 0],
                     [-0.5, 0, -0.5, 0, -0.5, 0, 0, -0.5, 0],
                     [-0.5, -0.5, 0, 0, 0, -0.5, 0, 0, -0.5],
                     [-0.5, 0, 0, 0, -0.5, -0.5, -0.5, 0, 0],
                     [0, -0.5, 0, -0.5, 0, -0.5, 0, -0.5, 0],
                     [0, 0, -0.5, -0.5, -0.5, 0, 0, 0, -0.5],
                     [-0.5, 0, 0, -0.5, 0, 0, 0, -0.5, -0.5],
                     [0, -0.5, 0, 0, -0.5, 0, -0.5, 0, -0.5],
                     [0, 0, -0.5, 0, 0, -0.5, -0.5, -0.5, 0]
                   ]
  describe "totalEnergy & totalMagnetization" $ do
    it "correctly calculates total energy & magnetization" $ do
      let !h₁ = ferromagneticIsingModelSquare2D 3 0
      totalEnergy h₁ (Configuration 9 [0b100110011]) `shouldBe` 6
      totalEnergy h₁ (Configuration 9 [0b000000111]) `shouldBe` (-6)
      totalEnergy h₁ (Configuration 9 [0b110000001]) `shouldBe` 2
  describe "randomConfigurationM" $ do
    it "generates random bit strings" $ do
      g <- mkXoshiro256PlusPlus 42
      state <- randomConfigurationM 9 g
      print state
      state `shouldSatisfy` (\(Configuration 9 x) -> G.length x == 1)
      state `shouldSatisfy` (\(Configuration 9 x) -> 0 < G.head x && G.head x < bit 9)
