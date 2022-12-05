{
  name,
  availableKernels,
  extraArgs,
}:
availableKernels.haskell {
  name = "${name}";
  inherit (extraArgs) system;
  haskellCompiler = "ghc924";
  extraHaskellPackages = ps: [
    ps.hvega
    ps.ihaskell-hvega
    ps.text
    ps.bytestring
    ps.vector
    ps.aeson
    ps.aeson-pretty
  ];
}
