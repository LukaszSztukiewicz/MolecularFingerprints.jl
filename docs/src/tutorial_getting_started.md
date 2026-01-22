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

```@example main
using MolecularFingerprints

# 1. Input: SMILES string (Benzene)
smiles = "C1=CC=CC=C1"

# 2. Configuration: ECFP (Extended Connectivity Fingerprints)
# Parameters: <Bit-length>(Radius)
calc = ECFP{1024}(2) 

# 3. Execution: Compute the fingerprint
vector = fingerprint(smiles, calc)

# 4. Analysis: Find indices of active features
active_bits = findall(vector)
println("Active bit indices: ", active_bits)

```

### High-Throughput Processing

For large datasets, the package provides a vectorized implementation that leverages multithreading.

```@example hightroughput

using MolecularFingerprints

calc = ECFP{1024}(2) 

# A list of SMILES (e.g., from a CSV)
dataset = ["CCO", "C1=CC=CC=C1", "CC(=O)O"]

# The vectorized call automatically parallelizes over available threads
batch_vectors = fingerprint(dataset, calc)

println("Active bit indices in 1st: ", findall(batch_vectors[1]))
println("Active bit indices in 2nd: ", findall(batch_vectors[2]))
println("Active bit indices in 3rd: ", findall(batch_vectors[3]))

```

If you have never used molecular fingerprints before, see [Explanation: What are Fingerprints?](explanation_fingerprints.md) for a detailed introduction.

For more detailed examples and advanced usage, please refer to the [API Reference: Public API](public_api_reference.md) and tutorials on [Solubility Prediction](tutorial_solubility_prediction.md) and [Similarity Search](tutorial_similarity_search.md).
