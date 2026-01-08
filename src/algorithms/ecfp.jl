using MolecularGraph

struct ECFP{N} <: AbstractFingerprint
    radius::Int

    function ECFP{N}(radius::Int = 2) where N
        radius >= 0 || throw(ArgumentError("radius must be non-negative"))
        N > 0 || throw(ArgumentError("fingerprint size must be positive"))
        new{N}(radius)
    end
end

function get_atom_invariant(mol, atom_index)
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
    delta_mass = round(Int, atom_mass - standard_mass)

    # Get formal charge
    at_charge = atom_charge(atom)

    # Check if atom is in a ring
    ring_info = is_in_ring(mol)
    in_ring = ring_info[atom_index]

    # Combine into a tuple (will be hashed)
    return (
        at_number,
        total_degree,
        total_hs,
        at_charge,
        delta_mass,
        in_ring,
    )
end

function hash_invariant(invariant)
    return hash(invariant, UInt(0xECFECF00)) % UInt32
end

function get_neighborhood_hash(current_hash, neighbor_hashes)
    # Sort neighbor hashes for canonical ordering
    sorted_neighbors = sort(neighbor_hashes)

    # Combine current hash with sorted neighbor hashes
    combined = (current_hash, sorted_neighbors...)

    return hash(combined, UInt(0xECFECF00)) % UInt32
end

function fingerprint(mol::SMILESMolGraph, calc::ECFP{N}) where N
    n_atoms = length(mol.vprops)

    # Handle empty molecule
    if n_atoms == 0
        return falses(N)
    end

    # Initialize atom identifiers
    atom_hashes = Vector{UInt32}(undef, n_atoms)
    for i in 1:n_atoms
        invariant = get_atom_invariant(mol, i)
        atom_hashes[i] = hash_invariant(invariant)
    end

    # Collect all features (hashes at each iteration)
    features = Set{UInt32}()

    # Add initial atom hashes (iteration 0)
    for h in atom_hashes
        push!(features, h)
    end

    # Iterate for the specified radius
    for _ in 1:calc.radius
        # Store new hashes for all atoms
        new_hashes = Vector{UInt32}(undef, n_atoms)

        for atom_index in 1:n_atoms
            # Get neighbor's hashes
            neighbor_indices = mol.graph.fadjlist[atom_index]
            neighbor_hashes = [atom_hashes[neighbor_index] for neighbor_index in neighbor_indices]

            # Compute new hash for this atom
            new_hashes[atom_index] = get_neighborhood_hash(atom_hashes[atom_index], neighbor_hashes)

            # Add to features
            push!(features, new_hashes[atom_index])
        end

        # Update hashes for next iteration
        atom_hashes = new_hashes
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

export ECFP, fingerprint
