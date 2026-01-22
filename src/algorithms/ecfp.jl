"""
    ECFP{N} <: AbstractFingerprint

Extended-Connectivity Fingerprint (ECFP) calculator.

ECFPs are circular fingerprints encoding a local molecular environment around each atom
up to a specified radius. This implementation closely follows the RDKit algorithm.

# Fields
- `radius::Int`: The maximum number of bonds to traverse from each atom (default: 2)

# Type Parameters
- `N`: The size of the fingerprint bit vector

# Examples
```julia
# Create ECFP4 (radius=2) with 2048 bits
fp_calc = ECFP{2048}(2)

# Create ECFP6 (radius=3) with 1024 bits
fp_calc = ECFP{1024}(3)
```

# References
Rogers, D., & Hahn, M. (2010). Extended-connectivity fingerprints.
Journal of Chemical Information and Modeling, 50(5), 742-754.
"""
struct ECFP{N} <: AbstractFingerprint
    radius::Int

    function ECFP{N}(radius::Int = 2) where N
        radius >= 0 || throw(ArgumentError("radius must be non-negative"))
        N > 0 || throw(ArgumentError("fingerprint size must be positive"))
        new{N}(radius)
    end
end

"""
    ecfp_atom_invariant(smiles::AbstractString)
    ecfp_atom_invariant(mol::AbstractMolGraph)
    ecfp_atom_invariant(mol::AbstractMolGraph, atom_index)

Calculate atomic invariants for ECFP fingerprint generation.

The atomic invariants are properties of an atom that don't depend on initial atom numbering,
based on the Daylight atomic invariants. This implementation follows the RDKit approach.

# Arguments
- `smiles::AbstractString`: SMILES string representation of a molecule
- `mol::AbstractMolGraph`: Molecular graph structure
- `atom_index`: Index of the specific atom to compute invariants for. If not specified, invariants for all atoms are computed and returned

# Returns
- For single atom: `Vector{UInt32}` containing the invariant components
- For all atoms: `Vector{Vector{UInt32}}` with invariants for each atom

# Invariant Components
The computed invariants include (in order):
1. Atomic number
2. Total degree (number of neighbors including implicit hydrogens)
3. Total number of hydrogens (implicit + explicit)
4. Atomic charge
5. Delta mass (difference from standard isotope mass)
6. Ring membership indicator (1 if atom is in a ring, omitted otherwise)

# References
RDKit implementation: https://github.com/rdkit/rdkit/blob/Release_2025_09_4/Code/GraphMol/Fingerprints/FingerprintUtil.cpp#L244
"""
ecfp_atom_invariant(smiles::AbstractString) = ecfp_atom_invariant(smilestomol(smiles))

ecfp_atom_invariant(mol::AbstractMolGraph) = [ecfp_atom_invariant(mol, i) for i in 1:nv(mol)]

function ecfp_atom_invariant(mol::AbstractMolGraph, atom_index)
    atom = mol.vprops[atom_index]

    # Get number of implicit and explicit hydrogens
    implicit_hs = implicit_hydrogens(mol)[atom_index]
    explicit_hs = explicit_hydrogens(mol)[atom_index]
    total_hs = implicit_hs + explicit_hs

    # Get number of neighbor atoms (including the implicit hydrogens)
    total_degree = degree(mol.graph, atom_index) + implicit_hs

    # Get valence
    valence_info = valence(mol)
    v = valence_info[atom_index]

    # Get atomic number
    at_number = atom_number(atom)

    # Get difference from atomic mass to average standard weight
    atom_mass = exact_mass(atom)
    standard_mass = monoiso_mass(atom)
    delta_mass = trunc(Int, atom_mass - standard_mass)

    # Get formal charge
    at_charge = atom_charge(atom)

    # Check if atom is in a ring
    ring_info = is_in_ring(mol)
    in_ring = ring_info[atom_index]

    # Return all invariants as an UInt32 Vector. Add an extra component if we are in a ring.
    components = UInt32[
        at_number, # atomic number
        total_degree, # number of neighbors (including implicit hydrogens)
        total_hs, # total number of hydrogens
        at_charge, # atomic charge
        delta_mass, # difference between atom and standard mass
    ]

    if in_ring
        push!(components, UInt32(1))
    end

    return components
end

"""
    ecfp_hash_combine(seed::UInt32, value::UInt32)

Combine two hash values using the boost hash_combine algorithm.

This function implements the hash combining strategy used in RDKit's ECFP implementation,
which is based on the boost C++ library's hash_combine function.

# Arguments
- `seed::UInt32`: Current hash seed value
- `value::UInt32`: New value to combine into the hash

# Returns
- `UInt32`: Combined hash value

# References
Boost hash implementation, as provided by RDKit: https://github.com/rdkit/rdkit/blob/Release_2025_09_4/Code/RDGeneral/hash/hash.hpp
"""
function ecfp_hash_combine(seed::UInt32, value::UInt32)
    return seed ⊻ (value + UInt32(0x9e3779b9) + (seed << 6) + (seed >> 2))
end

"""
    ecfp_hash(v::Vector{UInt32})

Generate a hash value from a vector of UInt32 values.

Iteratively combines all values in the vector using the ECFP hash combining algorithm
to produce a single hash value representing the entire vector.

# Arguments
- `v::Vector{UInt32}`: Vector of values to hash

# Returns
- `UInt32`: Hash value representing the input vector

# References
Boost hash implementation, as provided by RDKit: https://github.com/rdkit/rdkit/blob/Release_2025_09_4/Code/RDGeneral/hash/hash.hpp
"""
function ecfp_hash(v::Vector{UInt32})
    seed = UInt32(0)
    for value in v
        seed = ecfp_hash_combine(seed, value)
    end
    return seed
end

"""
    MorganAtomEnv

Internal structure representing a Morgan atom environment.

Stores the hash code, atom identifier (index), and layer/radius for each atomic environment
encountered during ECFP fingerprint generation.

# Fields
- `code::UInt32`: Hash code representing the atomic environment
- `atom_id::Int`: Identifier of the central atom
- `layer::Int`: Radius/layer at which this environment was computed
"""
struct MorganAtomEnv
    code::UInt32
    atom_id::Int
    layer::Int

    MorganAtomEnv(code::UInt32, atom_id::Int, layer::Int) = new(code, atom_id, layer)
end

"""
    AccumTuple

Internal structure for tracking and comparing atomic neighborhoods during ECFP generation.

Used to detect duplicate neighborhoods and maintain consistency with RDKit's algorithm
by storing bond connectivity patterns along with invariant hashes.

# Fields
- `bits::BitVector`: Bit representation of the bond neighborhood
- `invariant::UInt32`: Hash invariant for this neighborhood
- `atom_index::Int`: Index of the central atom
"""
struct AccumTuple
    bits::BitVector
    invariant::UInt32
    atom_index::Int

    AccumTuple(bits::BitVector, invariant::UInt32, atom_index::Int) = new(bits, invariant, atom_index)
end


function Base.isless(a::AccumTuple, b::AccumTuple)
    # Compare ra and rb in reverse order, because this is how boost::dynamic_bitset does it
    for (abit, bbit) in zip(Iterators.reverse(a.bits), Iterators.reverse(b.bits))
        abit != bbit && return abit < bbit
    end
    length(a.bits) != length(b.bits) && return length(a.bits) < length(b.bits) # Should not occur in our case, but let's put it here for completeness

    a.invariant != b.invariant && return a.invariant < b.invariant
    a.atom_index < b.atom_index
end

"""
    get_atom_invariants(mol::MolGraph)

Compute hashed atom invariants for all atoms in a molecule.

# Arguments
- `mol::MolGraph`: Input molecular graph

# Returns
- `Vector{UInt32}`: Hashed invariant values for each atom

# References
RDKit implementation: https://github.com/rdkit/rdkit/blob/Release_2025_09_4/Code/GraphMol/Fingerprints/MorganGenerator.cpp#L42
"""
function get_atom_invariants(mol::MolGraph)
    return [ecfp_hash(x) for x in ecfp_atom_invariant(mol)]
end

"""
    rdkit_bond_type(bond::SMILESBond)

Convert a SMILES bond to RDKit's bond type encoding.

Maps bond properties to integer codes matching RDKit's bond type enumeration.

# Arguments
- `bond::SMILESBond`: Input bond object

# Returns
- `Int`: Bond type code (1-6 for single to hextuple, 12 for aromatic, 20 for other, 21 for zero)

# Known Issues
Due to differences in the internal representation of bonds within MolecularGraph.jl,
we currently only support the most common bond types (1 to 6).

# References
RDKit bond types: https://github.com/rdkit/rdkit/blob/Release_2025_09_4/Code/GraphMol/Bond.h#L55
"""
function rdkit_bond_type(bond::SMILESBond)
    if !bond.isaromatic && bond.order in 1:6
        return bond.order # SINGLE..HEXTUPLE
    else
        error("Unsupported bond type")
    end

    # This is not well tested on out of scope for this project
    # if bond.isaromatic
    #     return 12 # AROMATIC
    # elseif bond.order == 0
    #     return 21 # ZERO
    # elseif bond.order in 1:6
    #     return bond.order # SINGLE..HEXTUPLE
    # else
    #     return 20 # OTHER
    # end
end

"""
    get_bond_invariants(mol::MolGraph)

Compute bond type invariants for all bonds in a molecule.

# Arguments
- `mol::MolGraph`: Input molecular graph

# Returns
- `Vector{UInt32}`: Bond type codes for each bond in the molecule

# Known Issue
The edge properties provided by [MolecularGraph.jl](https://github.com/mojaie/MolecularGraph.jl) are not in the same order as in RDKit.
This results in different hashes and, ultimately, in different fingerprints for larger
molecules compared to RDKit. As this would require rework on the smilestomol algorithm
provided by [MolecularGraph.jl](https://github.com/mojaie/MolecularGraph.jl), a fix for this issue is currently not in scope of this project.

# References
RDKit implementation: https://github.com/rdkit/rdkit/blob/Release_2025_09_4/Code/GraphMol/Fingerprints/MorganGenerator.cpp#L126
"""
function get_bond_invariants(mol::MolGraph)
    return [UInt32(rdkit_bond_type(bond)) for (_, bond) in mol.eprops]
end

"""
    fingerprint(mol::MolGraph, calc::ECFP{N}) where N

Generate an ECFP (Extended-Connectivity Fingerprint) for a molecule.

This function implements the Morgan/ECFP algorithm as described in the original paper
and matching the RDKit implementation. It generates circular fingerprints by iteratively
expanding atomic neighborhoods up to the specified radius.

# Algorithm Overview
1. Compute initial atom invariants (layer 0)
2. For each layer up to the specified radius:
   - Expand atomic neighborhoods by one bond
   - Hash neighborhood information to create new invariants
   - Detect and eliminate duplicate neighborhoods
   - Store unique atomic environments
3. Map all environment hashes to bit positions in the fingerprint

# Arguments
- `mol::MolGraph`: Input molecular graph
- `calc::ECFP{N}`: ECFP calculator specifying radius and fingerprint size

# Returns
- `BitVector`: Binary fingerprint of length N with bits set for detected molecular features

# Examples
```julia
mol = smilestomol("CCO")  # Ethanol
fp_calc = ECFP{2048}(2)   # ECFP4 with 2048 bits
fp = fingerprint(mol, fp_calc)
```

# References
- Rogers, D., & Hahn, M. (2010). Extended-connectivity fingerprints. J. Chem. Inf. Model., 50(5), 742-754.
- RDKit implementation: https://github.com/rdkit/rdkit/blob/Release_2025_09_4/Code/GraphMol/Fingerprints/MorganGenerator.cpp#L257
"""
function fingerprint(mol::MolGraph, calc::ECFP{N}) where N
    num_atoms = nv(mol)
    num_bonds = ne(mol)
    ernk = edge_rank(mol)

    # Generate atom hashes
    atom_invariants::Vector{UInt32} = get_atom_invariants(mol)
    bond_invariants::Vector{UInt32} = get_bond_invariants(mol)

    result::Vector{MorganAtomEnv} = []

    current_invariants::Vector{UInt32} = copy(atom_invariants)
    next_layer_invariants::Vector{UInt32} = zeros(num_atoms)

    neighborhood_invariants::Vector{Pair{UInt32, UInt32}} = []

    neighborhoods::Set{BitVector} = Set()
    atom_neighborhoods::Vector{BitVector} = [falses(num_bonds) for _ in 1:num_atoms]
    round_atom_neighborhoods::Vector{BitVector} = [falses(num_bonds) for _ in 1:num_atoms]

    # These atoms are skipped
    dead_atoms = falses(num_atoms)

    # Add round 0 invariants
    for (i, hash) in enumerate(current_invariants)
        push!(result, MorganAtomEnv(hash, i, 0)) # Note: atom id should start at 0
    end

    # Do subsequent rounds
    for layer in 1:calc.radius
        all_neighborhoods_this_round::Vector{AccumTuple} = []

        for atom_index in 1:num_atoms
            # If the atom is marked as dead, skip
            if dead_atoms[atom_index]
                continue
            end

            # If the atom has no neighbors, also skip
            neighbor_indices = mol.graph.fadjlist[atom_index]
            if isempty(neighbor_indices)
                dead_atoms[atom_index] = true
                continue
            end

            # Add up-to-date variants of neighbors
            empty!(neighborhood_invariants)

            for neighbor_index in neighbor_indices
                bond_index = edge_rank(ernk, atom_index, neighbor_index)
                round_atom_neighborhoods[atom_index][bond_index] = true;

                round_atom_neighborhoods[atom_index] .|= atom_neighborhoods[neighbor_index]
                push!(neighborhood_invariants, Pair(
                    bond_invariants[bond_index],
                    current_invariants[neighbor_index],
                ))
            end

            # Sort neighbor list
            sort!(neighborhood_invariants)

            # Calculate the new atom invariant by combining the neighbors'
            invar::UInt32 = layer - 1
            invar = ecfp_hash_combine(invar, current_invariants[atom_index])

            for (bond_invar, neighbor_invar) in neighborhood_invariants
                # hash the pair separately first, then combine with the total invariant (this is how rdkit does it)
                nb_invar = ecfp_hash([bond_invar, neighbor_invar])
                invar = ecfp_hash_combine(invar, nb_invar)
            end

            next_layer_invariants[atom_index] = invar;
            push!(all_neighborhoods_this_round, AccumTuple(
                round_atom_neighborhoods[atom_index],
                invar,
                atom_index,
            ))
        end

        # Sort the boolean array in reverse order, the other tuple members in forward order.
        # This is consistent with the way C++ sorts the Bitstream member
        sort!(all_neighborhoods_this_round)

        for nbh in all_neighborhoods_this_round
            if nbh.bits in neighborhoods
                # If this bit sequence has been seen before, mark the atom as dead
                dead_atoms[nbh.atom_index] = true
            else
                # Add the neighborhood invariant to the result list
                push!(result, MorganAtomEnv(
                    nbh.invariant,
                    nbh.atom_index,
                    layer,
                ))
                push!(neighborhoods, nbh.bits)
            end
        end

        # Swap current and next layer invariant vectors
        current_invariants = copy(next_layer_invariants)
        next_layer_invariants = zeros(num_atoms)

        # Start with this round's neighbors on the next rounds
        atom_neighborhoods = round_atom_neighborhoods
    end

    fp::BitVector = falses(N)

    for env in result
        seed::UInt32 = env.code

        bit_id = seed
        if N != 0
            bit_id = mod(bit_id, N)
        end

        # We could also use an integer vector for the fingerprint here, in that case we would
        # add one and get a count instead of just a bit map (rdkit provides a configuration option for this)
        # Note: use 1-based indexing
        fp[bit_id + 1] = true
    end

    return fp
end

