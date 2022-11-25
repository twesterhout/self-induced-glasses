{ enableStatic ? false,
  enableProfiling ? false
}:
let
  pkgs = import (import nix/sources.nix { }).nixpkgs { };
  # hp = nixpkgs.pkgsStatic.haskellPackages.ghc9_2_4;
  hpRaw = if enableStatic then pkgs.pkgsStatic.haskell.packages.native-bignum.ghc924
                          else pkgs.haskell.packages.ghc924;

  cabal2nixOptions = "" + (pkgs.lib.optionalString enableProfiling "--enable-profiling");
  cabalOptions = "--ghc-options=-fexpose-all-unfoldings -fspecialise-aggressively";
  
  hp = hpRaw.override {
    overrides = self: super: {
      unliftio = super.unliftio.overrideAttrs (oldAttrs: {
        configureFlags = oldAttrs.configureFlags ++ [ cabalOptions ];
      });
      vector = super.vector.overrideAttrs (oldAttrs: {
        configureFlags = oldAttrs.configureFlags ++ [ cabalOptions ];
      });
      primitive = super.primitive.overrideAttrs (oldAttrs: {
        configureFlags = oldAttrs.configureFlags ++ [ cabalOptions ];
      });
    };
  };

  self-induced-glasses = 
    (hp.callCabal2nixWithOptions "self-induced-glasses" ./. cabal2nixOptions {}).overrideAttrs (oldAttrs: {
      configureFlags = oldAttrs.configureFlags ++ [ cabalOptions ];
    });

  pythonEnv = pkgs.python3.withPackages(ps: [ ps.numpy ps.scipy ]);

  projectShell = hp.shellFor {
    packages = p: [ self-induced-glasses ];
    buildInputs = [
      hp.cabal-install
      pythonEnv
    ];
  };
in
  if pkgs.lib.inNixShell then projectShell else self-induced-glasses
  # hp.callCabal2nixWithOptions "self-induced-glasses" (./.) "--shell" { }
  # hp.developPackage {
  #   root = ./.;
  #   name = "self-induced-glasses";
  #   withHoogle = false;
  #   modifier = drv:
  #     nixpkgs.pkgsStatic.haskell.lib.addBuildTools drv [
  #       hp.cabal-install
  #     ];
  # }
