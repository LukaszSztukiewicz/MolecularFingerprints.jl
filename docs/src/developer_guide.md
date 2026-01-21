# MolecularFingerprints.jl Developers Guide

## Contributing to MolecularFingerprints.jl

We have complete CONTRIBUTING guidelines in the [CONTRIBUTING.md](https://github.com/LukaszSztukiewicz/MolecularFingerprints.jl/blob/main/CONTRIBUTING.md) file. Please read it if you are interested in contributing to the project.

This document provides a brief overview of the miscellaneous aspects of contributing to the MolecularFingerprints.jl package.

## Documentation
The documentation for MolecularFingerprints.jl is built using Documenter.jl. To build the documentation locally, follow these steps:
```julia
activate docs/
add Documenter
include("docs/make.jl")
```
This will generate the documentation in the `docs/build/` directory.


## Testing
To run the tests for MolecularFingerprints.jl, you can use the following commands:
```julia
] activate test
test MolecularFingerprints
```
or
```julia
] activate .
test
```
or
```julia
julia --project=test test/runtests.jl
```
This will execute the test suite and report any failures or errors.
To set up the test environment and add necessary dependencies, you can use the following commands: 
```julia
activate test/
add YourPackageName
```
### Update Test Dependencies
To update the test dependencies, you can use the following commands from the root of the repository:
```julia
using Pkg
Pkg.activate("test")
Pkg.update()
Pkg.develop(path=".")
```
This will ensure that your test environment is up to date with the latest dependencies and changes in the main package.
