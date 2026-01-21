```@meta
CurrentModule = MolecularFingerprints
```

# MolecularFingerprints.jl API Reference

## Index of Available Functions and Types

```@index

```

## Core Interface

```@docs
MolecularFingerprints.fingerprint

```

## Type Hierarchy

### Abstract Types

```@docs
MolecularFingerprints.AbstractCalculator
MolecularFingerprints.AbstractFingerprint
MolecularFingerprints.AbstractDescriptor

```

### Concrete Types

```@docs
MolecularFingerprints.ECFP
MolecularFingerprints.MHFP
MolecularFingerprints.MACCS
MolecularFingerprints.TopologicalTorsion
```

## Miscellaneous

```@autodocs
Modules = [MolecularFingerprints]
Filter = t -> !(t in [
    MolecularFingerprints.fingerprint,
    MolecularFingerprints.AbstractCalculator,
    MolecularFingerprints.AbstractFingerprint,
    MolecularFingerprints.AbstractDescriptor,
    MolecularFingerprints.ECFP,
    MolecularFingerprints.MHFP,
    MolecularFingerprints.MACCS,
    MolecularFingerprints.TopologicalTorsion,
])

```
