{
  name,
  availableKernels,
  extraArgs,
}:
availableKernels.python {
  name = "${name}";
  inherit (extraArgs) system;
}
