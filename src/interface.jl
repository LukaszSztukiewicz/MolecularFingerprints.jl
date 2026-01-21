# --- Abstract Interfaces for Fingerprint Calculators ---

"""
    AbstractCalculator
An abstract type representing a general calculator for molecular properties.
We differentiate between calculators for fingerprints and descriptors.
"""
abstract type AbstractCalculator end


"""
    AbstractFingerprint
An abstract type representing a calculator for molecular fingerprints.

Fingerprint is a binary representation of molecular features.
It captures the presence or absence of specific substructures or properties in a molecule.
It is often used in cheminformatics for similarity searching, clustering, and machine learning tasks.

We differentiate it from descriptors, which are typically numerical values.
"""
abstract type AbstractFingerprint <: AbstractCalculator end

"""
    AbstractDescriptor
An abstract type representing a calculator for molecular descriptors.

Descriptors are numerical values that describe molecular properties. 
For example, molecular weight, logP, etc.
"""
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

