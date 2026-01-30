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

## Minimalistic Usage
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

## Usage of all Fingerprint Types
The MolecularFingerprints.jl package supports several types of molecular fingerprints. Here is an example of how to use all available fingerprint types:

```julia
smiles = "C1=CC=CC=C1"

# 2. This package implements 4 types of fingerprints. 
# All of them could be customized with parameters, but here we use default settings.
ecfp_calc = ECFP() # Extended Connectivity Fingerprints
mhfp_calc = MHFP() # MinHash Fingerprints
torsion_calc = TopologicalTorsion() # Topological Torsion Fingerprints
maccs_calc = MACCS() # MACCS Keys

# 3. Execution: Compute the fingerprint for each type
ecfp_vector = fingerprint(smiles, ecfp_calc)
mhfp_vector = fingerprint(smiles, mhfp_calc)
torsion_vector = fingerprint(smiles, torsion_calc)
maccs_vector = fingerprint(smiles, maccs_calc)

# 4. Analysis: Find indices of active features

# ECFP returns BitVector to see active bits, we can use findall
println("ECFP active bits: ", findall(ecfp_vector))

# MACCS returns BitVector to see active bits, we can use findall
println("MACCS active bits: ", findall(maccs_vector))

# MHFP returns Vector{Int64} with each non-zero entry, so all bits are active
# You will see that are of the 2048 bits are being listed
println("MHFP active bits: ", findall(mhfp_vector .!= 0))

# TopologicalTorsion returns SparseArrays.SparseVector{Int32, Int64} so it is easy to find non-zero entries
using SparseArrays
println("Topological Torsion active bits: ", SparseArrays.findnz(torsion_vector)[1])
```


# Documentation

The documentation for MolecularFingerprints.jl can be found at [https://molecularfingerprints.lukaszsztukiewicz.com/stable](https://molecularfingerprints.lukaszsztukiewicz.com/stable/).
