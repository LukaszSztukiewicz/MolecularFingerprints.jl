# Getting Started with MolecularFingerprints.jl

This guide will help you set up your environment and compute your first molecular representations.

## Installation

We recommend using Julia's built-in package manager (`Pkg`) to manage dependencies. Choose the method that best fits your workflow:

#### Option 1: Sandbox (Trial)

Best for a quick "Hello World" or testing a specific feature without modifying your global state.

```julia
using Pkg
Pkg.activate(temp=true)
Pkg.add(url="https://github.com/LukaszSztukiewicz/MolecularFingerprints.jl")
using MolecularFingerprints

```

#### Option 2: Project-Specific (Recommended)

Best for building reproducible research or production pipelines. This ensures your project’s dependencies are locked in a `Project.toml` file.

```julia
using Pkg
Pkg.activate(".") 
Pkg.add(url="https://github.com/LukaszSztukiewicz/MolecularFingerprints.jl")

```

---

## Usage

Molecular fingerprints are essentially feature extraction steps in a pipeline. The API is designed to be functional: you define a Calculator (the model) and apply it to your Data.

### Basic Pipeline

```julia
using MolecularFingerprints

# 1. Input: SMILES string (Benzene)
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
println("ECFP active bits: ", ecfp_vector)
println("MHFP active bits: ", mhfp_vector)
println("Topological Torsion active bits: ", torsion_vector)
println("MACCS active bits: ", maccs_vector)

```

### High-Throughput Processing

For large datasets, the package provides a vectorized implementation that leverages multithreading.

```julia
# A list of SMILES (e.g., from a CSV)
dataset = ["CCO", "C1=CC=CC=C1", "CC(=O)O"]

# The vectorized call automatically parallelizes over available threads
batch_vectors = fingerprint(dataset, calc)

```

If you have never used molecular fingerprints before, see [Explanation](explanation.md) for an introduction to the concept. 

For more detailed examples and advanced usage, please refer to the [API Reference](api_reference.md) and tutorials on [Solubility Prediction](tutorial_solubility_prediction.md) and [Similarity Search](tutorial_similarity_search.md).
