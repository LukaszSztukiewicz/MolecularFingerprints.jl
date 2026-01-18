using MolecularGraph: AbstractMolGraph, edge_rank

struct ECFP{N} <: AbstractFingerprint
    radius::Int

    function ECFP{N}(radius::Int = 2) where N
        radius >= 0 || throw(ArgumentError("radius must be non-negative"))
        N > 0 || throw(ArgumentError("fingerprint size must be positive"))
        new{N}(radius)
    end
end

ecfp_atom_invariant(smiles::AbstractString) = ecfp_atom_invariant(smilestomol(smiles))

ecfp_atom_invariant(mol::AbstractMolGraph) = [ecfp_atom_invariant(mol, i) for i in 1:nv(mol)]

function ecfp_atom_invariant(mol::AbstractMolGraph, atom_index)
    """
    From the ECFP Paper:

    The Daylight atomic invariants are six properties of an atom in a molecule that do not depend
    on initial atom numbering. These properties are:
      - The number of immediate neighbors who are "heavy" (non-hydrogen) atoms
      - The valence minus the number of hydrogens
      - The atomic number
      - The atomic mass
      - The atomic charge
      - The number of attached hydrogens (both implicit and explicit)
      - Whether the atom is contained in at least one ring

    RDKit uses specific combinations and orders of these variants, the exact implementation can be found here:
      https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Fingerprints/FingerprintUtil.cpp#L244
    """

    atom = mol.vprops[atom_index]

    # Get number of implicit and explicit hydrogens
    implicit_hs = implicit_hydrogens(mol, atom_index)
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

function ecfp_hash_combine(seed::UInt32, value::UInt32)
    return seed ⊻ (value + UInt32(0x9e3779b9) + (seed << 6) + (seed >> 2))
end

function ecfp_hash(v::Vector{UInt32})
    # return hash(invariant, UInt(0xECFECF00)) % UInt32

    # Reference: https://github.com/rdkit/rdkit/blob/master/Code/RDGeneral/hash/hash.hpp
    seed = UInt32(0)
    for value in v
        seed = ecfp_hash_combine(seed, value)
    end
    return seed
end

struct MorganAtomEnv 
    code::UInt32
    atom_id::Int
    layer::Int

    MorganAtomEnv(code::UInt32, atom_id::Int, layer::Int) = new(code, atom_id, layer)
end

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

function get_atom_invariants(mol::SMILESMolGraph)
    # Reference: https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Fingerprints/MorganGenerator.cpp#L42
    return [ecfp_hash(x) for x in ecfp_atom_invariant(mol)]
end

function rdkit_bond_type(bond::SMILESBond)
    # Reference: https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Bond.h#L55

    if bond.isaromatic
        return 12 # AROMATIC
    elseif bond.order == 0
        return 21 # ZERO
    elseif bond.order in 1:6
        return bond.order # SINGLE..HEXTUPLE
    else
        return 20 # OTHER
    end
end

function get_bond_invariants(mol::SMILESMolGraph)
    # Reference: https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Fingerprints/MorganGenerator.cpp#L126
    return [UInt32(rdkit_bond_type(bond)) for (_, bond) in mol.eprops]
end

function fingerprint(mol::SMILESMolGraph, calc::ECFP{N}) where N
    # Reference: https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Fingerprints/MorganGenerator.cpp#L257

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

export ecfp_atom_invariant, ecfp_hash, ECFP, fingerprint
