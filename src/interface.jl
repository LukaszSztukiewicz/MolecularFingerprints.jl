# --- Abstract Interfaces for Fingerprint Calculators ---
abstract type AbstractCalculator end
abstract type AbstractFingerprint <: AbstractCalculator end
abstract type AbstractDescriptor <: AbstractCalculator end

# --- Fingerprint calculation functions ---
"""
    fingerprint(smiles_list::Vector{String}, calc::AbstractCalculator)
Calculates the fingerprint for a list of SMILES strings using the specified calculator.

Uses parallelization over input vector.

# Returns a vector of fingerprints.
"""
function fingerprint(smiles_list::Vector{String}, calc::AbstractCalculator) 
    fingerprints = Vector{Any}(undef, length(smiles_list))

    # Use Threads for parallel computation
    Threads.@threads for i in eachindex(smiles_list)
        fingerprints[i] = fingerprint(smilestomol(smiles_list[i]), calc)
    end
    
    return fingerprints
 end 

"""
    fingerprint(smiles::String, calc::AbstractCalculator)
Calculates the fingerprint for a single SMILES string using the specified calculator.

# Returns the fingerprint.
"""
function fingerprint(smiles::String, calc::AbstractCalculator)
    
    # Convert SMILES to molecular graph
    mol = smilestomol(smiles)

    return fingerprint(mol, calc)
end

function fingerprint(mol::MolGraph, calc::ECFP{N}) where N
    return _fingerprint(mol, calc)
end

function fingerprint(mol::MolGraph, calc::MHFP)
    return _fingerprint(mol, calc)
end

function fingerprint(mol::MolGraph, calc::MACCSFingerprint)
    return _fingerprint(mol, calc)
end

function fingerprint(mol::MolGraph, calc::TopologicalTorsion)
    return _fingerprint(mol, calc)
end