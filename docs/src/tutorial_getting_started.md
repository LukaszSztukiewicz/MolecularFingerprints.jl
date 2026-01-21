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

### Using with Persistent Environment
To use the MolecularFingerprints.jl package in your own projects, you can add it to your project's environment. Open a terminal, navigate to your project directory, and run the following commands:
```bash
git clone https://github.com/LukaszSztukiewicz/MolecularFingerprints.jl
cd MolecularFingerprints.jl
julia --project=.
```

Then, in the Julia REPL, run:
```julia
using Pkg
Pkg.instantiate()
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