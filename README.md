# MolecularFingerprints

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://molecularfingerprints.lukaszsztukiewicz.com/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://molecularfingerprints.lukaszsztukiewicz.com/dev/)
[![Build Status](https://github.com/LukaszSztukiewicz/MolecularFingerprints.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/LukaszSztukiewicz/MolecularFingerprints.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/LukaszSztukiewicz/MolecularFingerprints.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/LukaszSztukiewicz/MolecularFingerprints.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

MolecularFingerprints.jl is a Julia package for calculating molecular fingerprints using various algorithms. It provides an easy-to-use interface for generating fingerprints from molecular structures, enabling efficient similarity searches, clustering, and machine learning applications in cheminformatics.

## Quick Start
```julia
(@v1.12) pkg> add https://github.com/LukaszSztukiewicz/MolecularFingerprints.jl
julia> using MolecularFingerprints
julia> using MolecularGraph  # Required to represent molecular structures
julia> fp = ECFP{2048}(3)  # Create an ECFP fingerprint generator of size 2048 with radius 3
julia> mol = smilestomol("CCO")  # Ethanol
julia> result = fingerprint(mol, fp)  # Calculate the fingerprint and output the resulting sparse BitVector
0
0
0
⋮
0
0
0
julia> findall(result)  # To see which bits are set, use this command
6-element Vector{Int64}:
   81
  223
  295
  808
 1058
 1411
```

## Documentation

The documentation for MolecularFingerprints.jl can be found at [https://molecularfingerprints.lukaszsztukiewicz.com/stable](https://molecularfingerprints.lukaszsztukiewicz.com/stable/).

To build the documentation locally, you can use the following commands:

```julia
activate docs/
add Documenter
include("docs/make.jl")
```

This will generate the documentation in the `docs/build/` directory.

## Testing

To run the tests for MolecularFingerprints.jl, you can use the following commands:

```julia
] activate .
test
```

This will execute the test suite and report any failures or errors.

To set up the test environment and add necessary dependencies, you can use the following commands: 

```julia
activate test/
add YourPackageName
```



