# MolecularFingerprints

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://molecularfingerprints.lukaszsztukiewicz.com/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://molecularfingerprints.lukaszsztukiewicz.com/dev/)
[![Build Status](https://github.com/LukaszSztukiewicz/MolecularFingerprints.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/LukaszSztukiewicz/MolecularFingerprints.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/LukaszSztukiewicz/MolecularFingerprints.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/LukaszSztukiewicz/MolecularFingerprints.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

MolecularFingerprints.jl is a Julia package for calculating molecular fingerprints using various algorithms. It provides an easy-to-use interface for generating fingerprints from molecular structures, enabling efficient similarity searches, clustering, and machine learning applications in cheminformatics.


# Quick Start

## Installation

There are two ways of using the MolecularFingerprints.jl package: 
1. Using it with temporary environment (best for trying out the package).
2. Using it with a persistent environment (best for using the package in your own projects).

### Using with Temporary Environment
You can try out the MolecularFingerprints.jl package without installing it permanently by using a temporary environment. Open a Julia REPL and run the following commands:
```julia
using Pkg
Pkg.activate(temp=true)
Pkg.add(url="https://github.com/LukaszSztukiewicz/MolecularFingerprints.jl")
using MolecularFingerprints
```

## Usage
Once you have installed the MolecularFingerprints.jl package, you can start using it to calculate molecular fingerprints. Here is a simple example:

```julia
using MolecularFingerprints
# Load a molecule from a SMILES string
molecule = "C1=CC=CC=C1"  # Benzene
# Choose the fingerprint calculator
fp = ECFP{1024}(2)  # Create an ECFP fingerprint generator of size 1024 with radius 2
# Calculate the fingerprint
fingerprint_vector = fingerprint(molecule, fp)
println(fingerprint_vector)
findall(fingerprint_vector)  # Indices of bits set to 1
```


# Documentation

The documentation for MolecularFingerprints.jl can be found at [https://molecularfingerprints.lukaszsztukiewicz.com/stable](https://molecularfingerprints.lukaszsztukiewicz.com/stable/).
