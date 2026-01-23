## Best Practices for Using MolecularFingerprints.jl

To get the most out of `MolecularFingerprints.jl`, keep the following engineering principles in mind:

## 1. Leverage Type Specialization

Julia's compiler specializes code based on types. By defining your calculator once (e.g., `calc = ECFP{2048}(3)`), you allow the compiler to optimize the internal loops for that specific bit-length. Avoid redefining the calculator inside loops.

## 2. Threading Performance

The batch `fingerprint` function uses `Threads.@threads`. To ensure you are actually using multiple cores, check your environment:

```julia
using Base.Threads
println(nthreads()) # Should be > 1

```

If it returns `1`, start Julia with `julia -t auto`.

## 3. Bit-Vector vs. Dense Vector

Note that fingerprints are returned as `BitVector` objects. This is memory-efficient (1 bit per element) and allows for high-speed logical operations (AND, OR, XOR) which are critical for calculating Tanimoto similarities.

## 4. Canonicalization

While the library handles SMILES parsing, remember that different SMILES strings can represent the same molecular graph (e.g., `C1=CC=CC=C1` vs `c1ccccc1`). For consistency in machine learning tasks, ensure your input data is canonicalized.