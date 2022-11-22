let
  nixpkgs = import (import nix/sources.nix { }).nixpkgs { };
  # hp = nixpkgs.pkgsStatic.haskell.compiler.ghc924 { };
  # hp = nixpkgs.pkgsStatic.haskellPackages.ghc9_2_4;
  hp = nixpkgs.haskell.packages.ghc924;
  #   overrides = self: super: {
  #     vty = (super.vty.override {
  #       terminfo = super.terminfo_0_4_1_5;
  #     });
  #   };
  # };
in
hp.callCabal2nix "random-ising-3d" ./. { }
