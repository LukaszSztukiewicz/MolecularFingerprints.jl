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

# --- Helper Functions ---
"""
    safe_smilestomol(smiles::String)

Attempts to parse a SMILES string. Returns nothing if it fails
instead of crashing the entire thread.
"""
function safe_smilestomol(smiles::String)
    try
        return smilestomol(smiles)
    catch e
        @warn "Parsing failed for: $smiles - Error: $e"
        return nothing
    end
end

"""
    smiles_to_neutralized_mol(smiles_string::String)
Convert a SMILES string to a neutralized `MolGraph` instance.
This function identifies the largest fragment in the SMILES string,
removes charges from common organic elements, and returns the corresponding `MolGraph`.
# Arguments
- `smiles_string`: A string representing the molecule in SMILES format.
# Returns
- A `MolGraph` instance of the neutralized largest fragment.
"""
function smiles_to_neutralized_mol(smiles_string::String)
    parts = split(smiles_string, ".")
    fragment = String(parts[argmax(length.(parts))])
    
    organic_subset = Set(["B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"])
    neutral_str = replace(fragment, r"\[([A-Z][a-z]?)[\+\-]*\d*H*\d*\]" => function (m)
        el = match(r"([A-Z][a-z]?)", m).captures[1]
        return el in organic_subset ? el : "[$el]"
    end)

    mol = safe_smilestomol(neutral_str)
    
    return mol
end



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
    mol = smiles_to_neutralized_mol(smiles)
    return fingerprint(mol, calc)
end

#FIXME should be implemented in each calculator
"""
    fingerprint(mol::nothing, calc::AbstractCalculator)
Calculate the fingerprint for a single `MolGraph` instance using the provided `calc`.
# Arguments
- `mol`: A `MolGraph` instance representing the molecule.
- `calc`: A subtype of `AbstractCalculator` defining the fingerprint type.
# Returns
- A fingerprint representation (the specific type depends on `calc`).
"""
function fingerprint(mol::Nothing, calc::AbstractCalculator)
    @warn "Molecule is invalid or could not be parsed. Returning default empty fingerprint."
    # Return a default empty fingerprint based on calculator type
    if calc isa AbstractFingerprint
        return calc isa MACCS ? BitVector(undef, 166) :
               calc isa ECFP ? BitVector(undef, 1024) : #make it calc.nbits
               calc isa MHFP ? Vector{Int64}(undef, 2048) : #make it calc.nbits
               calc isa TopologicalTorsion ? spzeros(Int32, 68719476736) : #68719476736-element SparseVector{Int32, Int64}
               error("Unknown fingerprint type for empty molecule.")
    elseif calc isa AbstractDescriptor
        return NaN  # or some other sentinel value for descriptors
    else
        error("Unknown calculator type.")
    end
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
        fingerprints[i] = fingerprint(smiles_list[i], calc)
    end
    
    return fingerprints
end

