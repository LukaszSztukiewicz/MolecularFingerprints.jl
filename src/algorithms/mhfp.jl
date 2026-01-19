"""
    MHFP

Class for MHFP fingerprint generator/featurizer. Contains settings and parameters for 
MHFP fingerprint generation.
"""
struct MHFP <: AbstractFingerprint
    radius::Int
    min_radius::Int
    rings::Bool
    n_permutations::Int
    max_hash::Int
    seed::Int
    _mersenne_prime::Int
    _permutations_a::Vector{UInt32}
    _permutations_b::Vector{UInt32}

    function MHFP(
        radius::Int,
        min_radius::Int, 
        rings::Bool, 
        n_permutations::Int, 
        max_hash::Int, 
        seed::Int,
        _mersenne_prime::Int,
        _permutations_a::Vector{UInt32},
        _permutations_b::Vector{UInt32})
        
        # Inner constructor to ensure mersenne prime is not given manually.
        _mersenne_prime == (1 << 61) - 1 || error(
            """Incorrect value set for the mersenne prime, must be $((1 << 61) - 1), 
            got $_mersenne_prime. Don't give this value explicitly.""")
        
        return new(
            radius, 
            min_radius,
            rings,
            n_permutations,
            max_hash,
            seed,
            _mersenne_prime,
            _permutations_a,
            _permutations_b)
    end

end

##### Outer constructur for the MHFP encoder ##############################################

# The most detailed call a user should ever do manually. All other attributes are calculated
# automatically.
function MHFP(
    radius::Int = 3, 
    min_radius::Int = 1,
    rings::Bool = true;
    n_permutations::Int = 2048,  # keyword arguments
    max_hash::Int = (1 << 32) - 1, 
    seed::Int = 42)
    
    _mersenne_prime = (1 << 61) -1

    _permutations_a = Vector{UInt32}()
    _permutations_b = Vector{UInt32}()

    Random.seed!(seed)

    for i in 1:n_permutations
        a = rand(UInt32(1):UInt32(max_hash))
        b = rand(UInt32(0):UInt32(max_hash))

        # redraw values if already present in _permutations_a
        while a in _permutations_a
            a = rand(UInt32(1):UInt32(max_hash))
        end

        # redraw values if already present in _permutations_b
        while b in _permutations_b
            b = rand(UInt32(0):UInt32(max_hash))
        end

        push!(_permutations_a, a)
        push!(_permutations_b, b)
    end
    
    return MHFP(
        radius,
        min_radius,
        rings,
        n_permutations,
        max_hash,
        seed,
        _mersenne_prime,
        _permutations_a,
        _permutations_b
    )
end



"""
    fingerprint(mol::SMILESMolGraph, calc::MHFP{N})

Calculates the MHFP fingerprint of the given molecule and returns it as a bit vector
"""
function fingerprint(mol::MolGraph, calc::MHFP)
    return mhfp_hash_from_molecular_shingling(  # TODO make mhfp_shingling_from_mol accept the encoder "calc" as a parameter instead of separate arguments, and change call here correspondingly
        mhfp_shingling_from_mol(mol, calc.radius, calc.rings, calc.min_radius),
        calc)  # FIXME This doesn't return a BitVector yet, but a vector of UInt32's
end

"""
    mhfp_shingling_from_mol(
        mol::MolGraph,
        radius::Int = 3,
        rings::Bool = true,
        min_radius::Int = 1)

Calculate the "molecular shingling" of a given molecule.

A molecular shingling is a vector of "SMILES"-strings, calculated from the ring 
structures and atom types of the molecule (optional), and the circular substructures 
around each heavy (=non-hydrogen) atom of the molecule.

# Arguments
- `mol::MolGraph`: the molecule for which to calculate the shingling.
- `radius::Int=3`: the radius up to which the substructures around each atom should 
    be considered
- `rings::Bool=true`: if true (default), extract the smallest set of smallest rings
    (sssr) of the molecule, and include the SMILES-string of each contained ring in 
    the molecular shingling
- `min_radius::Int=1`: minimum radius of the substructures around each atom to be 
    considered. Default is 1. If min_radius=0, include the SMILES-string of each 
    heavy (i.e., non-hydrogen) atom of the molecule is included in the molecular 
    shingling.
"""
function mhfp_shingling_from_mol(
    mol::MolGraph,
    radius::Int = 3,
    rings::Bool = true,
    min_radius::Int = 1,
)

    # remove all hydrogens in the molecule, as we only want to consider heavy atoms in all 
    # following steps
    remove_all_hydrogens!(mol)

    shingling = []

    # Consider rings of the molecule, if corresponding parameter is set
    if rings
        append!(shingling, smiles_from_rings(mol))
    end

    # If min_radius == 0, add SMILES string of all heavy atoms to shingling
    if min_radius == 0
        append!(shingling, smiles_from_atoms(mol))
        # Increase min_radius, as the case of radius=0 has now been dealt with.
        min_radius += 1
    end


    # Add SMILES strings of circular substructures around heavy atoms to the shingling
    append!(shingling, smiles_from_circular_substructures(mol, radius, min_radius))

    # Make shingling entries unique

    sort!(shingling)   # unique! is much faster on sorted arrays
    # (see https://docs.julialang.org/en/v1/base/collections/)
    
    unique!(shingling)

    return shingling
end

"""
    smiles_from_rings(mol::MolGraph)

Return vector containing SMILES strings of all rings in the SSSR of the given molecule.

SSSR stands for the smallest set of smallest rings of the molecule.

Note: This function uses the function sssr from MolecularGraph.jl, which returns a "true"
smallest set of smallest rings of the given molecule. However, in the original
implementation of the mhfp algorithm, the "symmetrisized sssr" is used, which in some cases
is non-minimal, i.e., contains an additional ring.
The rdkit function to get the symmetrisized sssr is not available in MolecularGraph.jl or in
RDKitMinimalLib, which is why the standard sssr is used.
In most cases, this will not have any effect, but for some molecules, such as cubane, it 
will.
"""
function smiles_from_rings(mol::MolGraph)
    shingling_snippet = []


    # Go through all rings in the sssr
    for ring in MolecularGraph.sssr(mol)
        # For each ring in the sssr, create smiles string of the submolecule corresponding
        # to the ring and add to the shingling.
        push!(shingling_snippet, smiles(induced_subgraph(mol, ring)[1]))
    end

    return shingling_snippet
end

"""
    smiles_from_atoms(mol::MolGraph)

Return vector containing SMILES strings of all atoms of the given molecule.
"""
function smiles_from_atoms(mol::MolGraph)
    shingling_snippet = []

    aromatic_atoms = is_aromatic(mol)

    for atom in vertices(mol)
        # create "molecule" containing only the current atom
        atom_as_mol = induced_subgraph(mol, [atom])[1]

        # Note: The original authors add the SMARTS string (not SMILES string) of the atoms
        # to the shingling. However, in the rdkit implementation, they use SMILES strings, 
        # see https://github.com/rdkit/rdkit/blob/da08e8d954c0a923300a43d24faa37ad224b320d/Code/GraphMol/Fingerprints/MHFP.cpp#L124
        
        # Get smiles string of atom
        smiles_string_of_atom = smiles(atom_as_mol)

        # the smiles string calculated in the line above is always a capital letter.
        # However, aromatic atoms are supposed to be written with a lowercase letter.
        # (The aromaticity of an atom is lost in the process above, as we are creating
        # an artificial molecule with only one atom (which cannot be aromatic).)
        # Therefore, the strings are manually converted to lowercase if the atom was
        # aromatic in the original molecule.
        if aromatic_atoms[atom] == 1
            smiles_string_of_atom = lowercase(smiles_string_of_atom)
        end
        
        push!(shingling_snippet, smiles_string_of_atom)
    end
    return shingling_snippet
end

"""
    smiles_from_circular_substructures(
        mol::MolGraph,
        radius::Int,
        min_radius::Int)

Return vector of SMILES strings of circular substructures around all atoms of a molecule.

For each atom of the given molecule, extract the substructures of radii min_radius to
radius, and generate their corresponding SMILES strings.
"""
function smiles_from_circular_substructures(mol::MolGraph, radius::Int, min_radius::Int)
    shingling_snippet = []

    for atom_index in vertices(mol)  # go through all atoms
        
        for i = min_radius:radius  # go through al selected radii
            
            atoms_in_substructure_of_radius_i = neighborhood(mol, atom_index, i)

            submol, atom_map = induced_subgraph(mol, atoms_in_substructure_of_radius_i)

            # TODO: This test is copied from the original authors.
            # I don't know what this test is for, as I don't see why it could be that 
            # atom_index is not contained in the atom map.
            # I guess this is just to make sure that we don't get an error but simply
            # continue. 
            # Maybe exclude this test for speedup in the future.
            if atom_index ∉ atom_map
                continue
                @warn "Current atom not found in atom map, skipping" atom_index atom_map
            end
            
            # Find index of current atom (atom_index) in the submolecule
            pos_of_atom_index_in_submol = findfirst(atom_map .== atom_index)

            smiles_of_substructure = smiles(
                submol,
                Dict{String,Any}("rootedAtAtom" => pos_of_atom_index_in_submol),
            )

            if smiles_of_substructure != ""
                # Add smiles of substructure to shingling.
                push!(shingling_snippet, smiles_of_substructure)
            end

        end
    end

    return shingling_snippet
end

"""
    mhfp_hash_from_molecular_shingling(shingling::Vector{String}, Encoder::MHFP)

Calculate the MinHash values from a given Molecular shingling.

The given encoder contains parameters such as the length of the hashed vector or the 
maximum permitted hash number, but also the seed of the rng, and the random numbers that
were generated using this seed, that are used in the hashing process.

The algorithm is described in more detail in the original authors paper.
"""
function mhfp_hash_from_molecular_shingling(shingling::Vector{String}, Encoder::MHFP)
    hash_values = zeros(UInt32, (Encoder.n_permutations))
    fill!(hash_values, Encoder.max_hash)

    for s in shingling
        sha_from_string = sha1(s)[begin:4]  # create sha1 hash from the string
        # Note: we are only using the first 4 sha1 bytes, as we want a 32-bit hash

        #let 
        buf = IOBuffer(sha_from_string)  # make the sha bytes a buffer
        s_h = Int(  # read bytes from sha hash into integer
            htol(  # ensure little-endian format of the integer
                 read(buf, UInt8)))  # read the buffer bytes as unsigned integer
        # end

        hashes = mod.(
            mod.(
                Encoder._permutations_a * s_h + Encoder._permutations_b,
                Encoder._mersenne_prime
            ), Encoder.max_hash
        )
        
        hash_values = min.(hash_values, hashes)

        
    end


    return hash_values
end
