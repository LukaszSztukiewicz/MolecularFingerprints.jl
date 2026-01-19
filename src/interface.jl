# --- Abstract Interfaces for Fingerprint Calculators ---
abstract type AbstractCalculator end
abstract type AbstractFingerprint <: AbstractCalculator end
abstract type AbstractDescriptor <: AbstractCalculator end

# --- Fingerprint calculation functions ---

"""
    fingerprint(smiles::String, calc::AbstractCalculator)
Calculates the fingerprint for a single SMILES string using the specified calculator.

# Returns the fingerprint.
"""
function fingerprint(smiles::String, calc::AbstractCalculator)
    
    # Convert SMILES to molecular graph #RDKitMinimalLib (MolecularGraph)
    mol = smilestomol(smiles)

    return fingerprint(mol, calc)
end

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

