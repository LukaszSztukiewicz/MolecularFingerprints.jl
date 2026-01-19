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
    
    # Convert first element to determine return type safely
    first_fp = fingerprint(smiles_list[1], calc)
    ResultType = typeof(first_fp)
    
    fingerprints = Vector{ResultType}(undef, length(smiles_list))
    fingerprints[1] = first_fp

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

# # --- Specific fingerprint calculation methods ---
# function fingerprint(mol::MolGraph, calc::ECFP{N}) where N
#     return _fingerprint(mol, calc::ECFP{N})
# end

# function fingerprint(mol::MolGraph, calc::MHFP)
#     return _fingerprint(mol, calc::MHFP)
# end

# function fingerprint(mol::MolGraph, calc::MACCSFingerprint)
#     return _fingerprint(mol, calc::MACCSFingerprint)
# end

# function fingerprint(mol::MolGraph, calc::TopologicalTorsion)
#     return _fingerprint(mol, calc::TopologicalTorsion)
# end