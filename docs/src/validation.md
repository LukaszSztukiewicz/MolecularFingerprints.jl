# Validation & Testing

How do we know if our implementation is good?

## 1. Validation Strategy -> Cross-Library Testing

When developing fingerprinting algorithms, we compare our results against industry standards like RDKit (C++). However, direct bit-vector comparison is not a valid test for correctness due to implementation details:

* Hashing Variance: Our native Julia implementation utilizes the internal `hash()` function, whereas RDKit uses specific PRNG (Pseudo-Random Number Generator) seeds in C++.
* Result: The specific indices of bits set will differ between libraries.

### Verification via Ranking Correlation

Instead of bitwise equality, we validate using Statistical Correlation. If the chemical logic (subgraph extraction) is identical, the Tanimoto similarity between pairs of molecules should be highly correlated across both libraries.

We verify that if Molecule A is "most similar" to Molecule B in RDKit, `MolecularFingerprints.jl` should produce the same ranking order, even if the underlying bit-vectors are different.


## 2. Testing Framework

We use the standard Julia `Test` module to ensure high code coverage and functional correctness.

### Running Tests

The most efficient way to run tests is via the Julia REPL. From the root of the repository:

```julia
# Method 1: The standard Pkg way (Recommended)
pkg> activate .
pkg> test

# Method 2: Running the test script directly from terminal
# julia --project=test test/runtests.jl

```

### Managing the Test Environment

The tests reside in an independent environment located in the `/test` directory. This keeps the main package dependencies lightweight by excluding testing-only packages (like `RDKitMinimalLib`) from the production environment.

Adding a new test-only dependency:

```julia
pkg> activate test
pkg> add Statistics  # Example: adding a stats package for validation

```

Syncing your local changes for testing:
If you make changes to the source code in `/src`, ensure the test environment is tracking your local version:

```julia
pkg> activate test
pkg> dev . 

```

### CI Integration

Every Pull Request is automatically tested against multiple Julia versions and operating systems via GitHub Actions. We also track Code Coverage; please ensure that any new fingerprint types added include corresponding tests in `test/runtests.jl`.