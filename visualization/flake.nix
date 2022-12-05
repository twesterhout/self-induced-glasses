{
  description = "Your jupyterWith project";

  nixConfig.extra-substituters = [
    "https://tweag-jupyter.cachix.org"
  ];
  nixConfig.extra-trusted-public-keys = [
    "tweag-jupyter.cachix.org-1:UtNH4Zs6hVUFpFBTLaA4ejYavPo5EFFqgd7G7FxGW9g="
  ];

  inputs.flake-compat.url = "github:edolstra/flake-compat";
  inputs.flake-compat.flake = false;
  inputs.flake-utils.url = "github:numtide/flake-utils";
  # inputs.nixpkgs.url = "github:nixos/nixpkgs/bc73da29ee0978fa544c3733348b172bceafa549";
  inputs.nixpkgs.url = "github:nixos/nixpkgs/e0990771a605839bc63eaef5edadffb48a040664";
  inputs.jupyterWith.url = "github:tweag/jupyterWith";
  inputs.jupyterWith.inputs.nixpkgs.follows = "nixpkgs";

  outputs = {
    self,
    flake-compat,
    flake-utils,
    nixpkgs,
    jupyterWith,
  }:
  flake-utils.lib.eachSystem [ flake-utils.lib.system.x86_64-linux ] (system:
    let
      inherit (jupyterWith.lib.${system}) mkJupyterlab mkJupyterlabFromPath;
      # mkJupyterlabFromPath ./kernels {inherit system;};
      pkgs = nixpkgs.legacyPackages.${system};
      haskellCompiler = "ghc924";

      cabalOptions = "--ghc-options=-fexpose-all-unfoldings -fspecialise-aggressively";
      cabal2nixOptions = enableProfiling:
        "" + (pkgs.lib.optionalString enableProfiling "--enable-profiling");

      mkHaskellPackages = { enableStatic, enableProfiling }:
        let
          hpRaw = if enableStatic
            then pkgs.pkgsStatic.haskell.packages.native-bignum.${haskellCompiler}
            else pkgs.haskell.packages.${haskellCompiler};
        in
          hpRaw.override {
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

      mkSelfInducedGlasses = { enableStatic, enableProfiling }:
        let
          hp = mkHaskellPackages { inherit enableStatic enableProfiling; };
          cabal2nixOptions = (pkgs.lib.optionalString enableProfiling "--enable-profiling");
        in
          (hp.callCabal2nixWithOptions "self-induced-glasses" ./.. cabal2nixOptions {}
          ).overrideAttrs (oldAttrs: {
            configureFlags = oldAttrs.configureFlags ++ [ cabalOptions ];
          });

      haskellPackagesNormal = mkHaskellPackages { enableStatic = false; enableProfiling = false; };
      haskellPackagesProfiling = mkHaskellPackages { enableStatic = false; enableProfiling = true; };
      haskellPackagesStatic = mkHaskellPackages { enableStatic = true; enableProfiling = false; };

      haskellKernel = k: k.haskell {
        inherit system pkgs;
        name = "haskell";
        displayName = "Haskell";
        # runtimePackages = [ haskellPackagesNormal ];
        haskellCompiler = "ghc924";
        extraHaskellPackages = ps: [
          (ps.callCabal2nix "self-induced-glasses" ./.. { })
          ps.hvega
          ps.ihaskell-hvega
          # (pkgs.haskell.lib.dontCheck ps.hspec-contrib)
          ps.text
          ps.bytestring
          ps.vector
          ps.binary
          ps.aeson
          ps.aeson-pretty
        ];
      };
      jupyterlab = mkJupyterlab {
        kernels = k: [ (haskellKernel k) ];
      };

      projectShell = haskellPackagesProfiling.shellFor {
        packages = p: [
          (mkSelfInducedGlasses { enableStatic = false; enableProfiling = true; })
        ];
        buildInputs = [
          haskellPackagesProfiling.cabal-install
        ];
      };

    in rec {
      # packages = {inherit jupyterlab;};
      # packages.jupyterlab = jupyterlab;
      packages.normal =
        mkSelfInducedGlasses { enableStatic = false; enableProfiling = false; };
      packages.profiling =
        mkSelfInducedGlasses { enableStatic = false; enableProfiling = true; };
      packages.static =
        mkSelfInducedGlasses { enableStatic = true; enableProfiling = false; };
      # packages.profiling = haskellPackagesProfiling.self-induced-glasses;
      # packages.static = haskellPackagesStatic.self-induced-glasses;
      apps.visualization.program = "${jupyterlab}/bin/jupyter-notebook";
      apps.visualization.type = "app";

      devShells.default = projectShell;
    }
  );
}
