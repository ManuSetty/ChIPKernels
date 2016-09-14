# ChIPKernels

ChIPKernels is an R package for building different string kernels used for DNA Sequence analysis. A dictionary of the desired kernel must be built and this dictionary can be used for determining kernels for DNA Sequences. This package has been tested with R 3.2

The following representations are supported.
  1. Mismatch/Dimismatch kernel: ```build.mismatch.dictionary```
  2. Positional wildcard kernel: ```build.wildcard.dictionary```

After the dictionary is built, use the ```build.features.kernels``` function to determine the feature matrix and kernel for DNA Sequences.
