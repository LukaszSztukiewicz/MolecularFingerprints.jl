Since our native Julia implementation uses Julia's hash() function, the specific bits set will differ from RDKit's C++ implementation (which uses specific random seeds). Therefore, strict bitwise equality (==) is impossible to test.

Instead, we validate by Ranking Correlation: We verify that MolecularFingerprints.jl identifies the same molecules as "similar" that RDKit does. If the chemical logic is correct, the Tanimoto similarity scores for a set of molecules should be highly correlated between the two libraries.

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