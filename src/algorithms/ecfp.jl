using MolecularGraph: AbstractMolGraph, edge_rank

struct ECFP{N} <: AbstractFingerprint
    radius::Int

    function ECFP{N}(radius::Int = 2) where N
        radius >= 0 || throw(ArgumentError("radius must be non-negative"))
        N > 0 || throw(ArgumentError("fingerprint size must be positive"))
        new{N}(radius)
    end
end

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

    # Return all invariants as an UInt32 Vector
    return convert(Vector{UInt32}, [
        at_number, # atom number
        total_degree, # number of neighbors (including implicit hydrogens)
        total_hs, # total number of hydrogens
        at_charge, # atomic charge
        delta_mass, # difference between atom and standard mass
        in_ring, # if the atom is part of a ring
    ])
end

function ecfp_hash(invariant)
    return hash(invariant, UInt(0xECFECF00)) % UInt32
end

function fingerprint(mol::SMILESMolGraph, calc::ECFP{N}) where N
    n_atoms = nv(mol)
    n_bonds = ne(mol)
    ernk = edge_rank(mol)

    # Handle empty molecule
    if n_atoms == 0
        return falses(N)
    end

    # Initialize atom identifiers
    atom_hashes = Vector{UInt32}(undef, n_atoms)
    for i in 1:n_atoms
        invariant = ecfp_atom_invariant(mol, i)
        atom_hashes[i] = ecfp_hash(invariant)
    end

    # Collect all features (hashes at each iteration)
    features = Set{UInt32}()

    # Track which bonds each atom's neighborhood includes
    atom_neighborhoods = [falses(n_bonds) for _ in 1:n_atoms]

    # Track neighborhoods we've already seen (as Sets or BitVectors)
    seen_neighborhoods = Set{BitVector}()

    # Track "dead" atoms that won't produce unique environments anymore
    dead_atoms = falses(n_atoms)

    # Add initial atom hashes (iteration 0)
    for h in atom_hashes
        push!(features, h)
    end

    # Iterate for the specified radius
    for _ in 1:calc.radius
        # Store new hashes for all atoms
        new_hashes = Vector{UInt32}(undef, n_atoms)
        round_neighborhoods = deepcopy(atom_neighborhoods)
        round_results = []  # (neighborhood, hash, atom_index)

        for atom_index in 1:n_atoms
            # Skip if the atom is marked as dead
            if dead_atoms[atom_index]
                continue
            end

            # If the atom has no neighbors, skip immediately
            neighbor_indices = mol.graph.fadjlist[atom_index]
            if isempty(neighbor_indices)
                dead_atoms[atom_index] = true
                continue
            end

            # Get neighbor's hashes
            neighbor_hashes = [atom_hashes[neighbor_index] for neighbor_index in neighbor_indices]

            for neighbor_index in neighbor_indices
                # Get bond index and mark it in this atom's neighborhood
                bond_idx = edge_rank(ernk, atom_index, neighbor_index)
                round_neighborhoods[atom_index][bond_idx] = true

                # Union with neighbor's previous neighborhood
                round_neighborhoods[atom_index] .|= atom_neighborhoods[neighbor_index]

                push!(neighbor_hashes, atom_hashes[neighbor_index])
            end

            # Sort neighbor hashes for canonical ordering
            sorted_neighbors = sort(neighbor_hashes)

            # Combine current hash with sorted neighbor hashes and store it as the new hash for this atom
            combined = (atom_hashes[atom_index], sorted_neighbors...)
            new_hashes[atom_index] = ecfp_hash(combined)

            # Add to features
            push!(round_results, (round_neighborhoods[atom_index], new_hashes[atom_index], atom_index))
        end

        # Sort and deduplicate
        sort!(round_results)
        for (neighborhood, hash, atom_idx) in round_results
            if neighborhood ∉ seen_neighborhoods
                push!(features, hash)
                push!(seen_neighborhoods, copy(neighborhood))
            else
                dead_atoms[atom_idx] = true
            end
        end

        # Update hashes for next iteration
        atom_hashes = new_hashes
        atom_neighborhoods = round_neighborhoods
    end

    # Create bit vector by folding features
    fp = falses(N)
    for feature in features
        # Use module to fold feature into bit vector
        bit_index = (feature % N) + 1 # +1 for 1-based indexing
        fp[bit_index] = true
    end

    return fp
end

export ecfp_atom_invariant, ecfp_hash, ECFP, fingerprint
