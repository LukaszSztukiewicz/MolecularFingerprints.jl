# --- Abstract Interfaces for Fingerprint Calculators ---

"""
    AbstractCalculator

Supertype for all molecular property calculators. Subtypes should implement 
specific calculation logic for molecular properties.
"""
abstract type AbstractCalculator end

"""
    AbstractFingerprint <: AbstractCalculator

Abstract type for calculators that produce representations 
of molecular features (e.g., MACCS, ECFP).

Unlike descriptors, fingerprints typically represent the presence or absence 
of specific substructures or patterns within a molecule.
"""
abstract type AbstractFingerprint <: AbstractCalculator end

"""
    AbstractDescriptor <: AbstractCalculator

Abstract type for calculators that produce scalar or numerical molecular properties 
(e.g., LogP, Molecular Weight, TPSA).
"""
abstract type AbstractDescriptor <: AbstractCalculator end

# --- Fingerprint calculation functions ---

"""
    fingerprint(smiles::String, calc::AbstractCalculator)

Calculate the fingerprint for a single SMILES string using the provided `calc`.

# Arguments
- `smiles`: A string representing the molecule in SMILES format.
- `calc`: A subtype of `AbstractCalculator` defining the fingerprint type.

# Returns
- A fingerprint representation (the specific type depends on `calc`).
"""
function fingerprint(smiles::String, calc::AbstractCalculator)
    mol = smilestomol(smiles)
    return fingerprint(mol, calc)
end

"""
    fingerprint(smiles_list::Vector{String}, calc::AbstractCalculator)

Calculate fingerprints for a collection of SMILES strings.

This method uses multithreading to process the list. Ensure that `JULIA_NUM_THREADS` 
is set appropriately in your environment to see performance gains.

# Arguments
- `smiles_list`: A vector of SMILES strings.
- `calc`: The calculator instance to apply to each molecule.

# Returns
- `Vector`: A collection of fingerprints, typed according to the first successful calculation.

!!! note
    This function is thread-parallelized using `Threads.@threads`.
"""
function fingerprint(smiles_list::Vector{String}, calc::AbstractCalculator) 
    # Determine return type from the first element
    first_fp = fingerprint(smiles_list[1], calc)
    ResultType = typeof(first_fp)
    
    fingerprints = Vector{ResultType}(undef, length(smiles_list))
    fingerprints[1] = first_fp

    Threads.@threads for i in 2:length(smiles_list)
        # Assuming smilestomol is thread-safe
        fingerprints[i] = fingerprint(smilestomol(smiles_list[i]), calc)
    end
    
    return fingerprints
end