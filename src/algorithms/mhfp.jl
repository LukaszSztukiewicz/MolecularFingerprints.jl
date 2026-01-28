"""
    MHFP(
        radius::Int = 3,
        min_radius::Int = 1,
        rings::Bool = true,
        fp_size::Int = 2048,
        seed::Int = 42
    )

Type for MHFP fingerprint calculators. Contains settings and parameters for 
MHFP fingerprint generation.

The MHFP fingerprint is a vector of UInt32's, calculated for a given molecule by:

1. generating the "molecular shingling" of the molecule, which is a set of strings, 
    containing:
    1. The SMILES strings of all rings in the smallest set of smallest rings (sssr) 
        of the molecule (optional, corresponds to setting `rings=true` in the MHFP 
        calculator object),
    2. The SMILES strings of the circular substructures of radii `min_radius:radius` around
        each heavy atom of the molecule. Note: if `min_radius=0`, the corresponding 
        substructures are just the atoms themselves. 
2. Hashing the molecular shingling, which consists of:
    1. Converting each string to a 32-bit integer using SHA1 (and only using the first 32
        bits of the hashed result)
    2. Applying the MinHash scheme to the set of 32-bit integers in order to generate the 
        final fingerprint. The exact formula is given in the original authors paper, but
        we note here that it takes a vector of 32-bit integers as input, and is furthermore
        dependent on two vectors a and b, each of a given length k, which is also the length
        of the resulting fingerprint vector. The two vectors are sampled at random, but
        must be the same for comparable fingerprints. Note: in the fields of MHFP objects,
        the vectors a, b and their length k are named `_a`, `_b` 
        and `fp_size`, respectively.


The avaliable parameters of the calculator object are:
# Given as arguments to the constructor:
- `radius::Int`: The *maximum* radius of circular substructures around each heavy atom 
    of a molecule that are to be included in the fingerprint. Recommended values are 2 or 3
    according to the original authors, with 3 (default) giving best results.
- `min_radius::Int`: The *minimum* radius of circular substructures around each heavy 
    atom of a molecule that are to be considered. Will be 1 (default) in most cases, however
    0 is also valid; in this case information about the heavy atoms of the molecules is 
    included explicitly in the fingerprints. The original paper only considers the case
    `min_radius=1`.
- `rings::Bool`: If true (default), information about rings in the molecules is included in 
    the fingerprints explicitly. This matches the original authors description of the 
    fingerprint in their paper.

# Given as keyword arguments to the constructor:
- `fp_size::Int`: length of the fingerprint. Also means that this is the length of the 
    random vectors a and b which are used in the hashing process. Default is 2048, as 
    recommended by the original authors in their paper.
- `seed::Int`: seed for the generation of the random vectors `a` and `b` which are used
     in the hashing process. Must be the same for comparable fingerprints. Default is 42.

Also contains the fields `_mersenne_prime`, `_max_hash`, `_a` and 
`_b`, which are internal and cannot be set explicitly.
The first two are constants, and the second two are random vectors which are generated 
automatically based on the given `seed`.

# References
- [Probst D., Reymond, JL. 2018](https://doi.org/10.1186/s13321-018-0321-8)
"""
struct MHFP <: AbstractFingerprint
    radius::Int
    min_radius::Int
    rings::Bool
    fp_size::Int
    seed::Int
    _mersenne_prime::UInt64  # strict types below, as fields are set automatically anyway.
    _max_hash::UInt32
    _a::Vector{UInt32}
    _b::Vector{UInt32}

    # Inner constructor. Checks for invalid given parameters, initializes the constants
    # _max_hash and _mersenne_prime and generates the random vectors a and b based on the 
    # given seed.
    function MHFP(
        radius::Integer = 3, 
        min_radius::Integer = 1,
        rings::Bool = true;
        fp_size::Integer = 2048,  # keyword arguments
        seed::Integer = 42
        )
        
        ### Ensure given values are valid
        # Ensure radius and min_radius are non-negative
        radius ≥ 0 || error("""Given radius must be non-negative, got radius=$radius.""")
        min_radius ≥ 0 || error(
            """Given min_radius must be non-negative, got min_radius=$min_radius.""")

        # Ensure radius is at least as large as min_radius
        radius ≥ min_radius || error(
            """Given radius must be larger or equal to given min_radius.
            Got radius=$radius but min_radius=$min_radius.""")

        # Ensure fp_size is positive
        fp_size > 0 || error("""fp_size must be strictly positive. Got 
            fp_size=$fp_size.""")

        ### set fixed values
        _mersenne_prime = (1 << 61) -1
        _max_hash = (1 << 32) - 1

        ### generate vectors a, b
        # initialize vectors
        _a = Vector{UInt32}()
        _b = Vector{UInt32}()

        # set seed
        seed!(seed)

        # fill vectors entry by entry, to ensure pairwise unique entries within the vectors
        for i in 1:fp_size
            a = rand(UInt32(1):UInt32(_max_hash))
            b = rand(UInt32(0):UInt32(_max_hash))

            # redraw values if already present in _a
            while a in _a
                a = rand(UInt32(1):UInt32(_max_hash))
            end

            # redraw values if already present in _b
            while b in _b
                b = rand(UInt32(0):UInt32(_max_hash))
            end

            push!(_a, a)
            push!(_b, b)
        end
        
        return new(
            radius, 
            min_radius,
            rings,
            fp_size,
            seed,
            _mersenne_prime,
            _max_hash,
            _a,
            _b)
    end

end

"""
    helper_custom_on_update!(mol::SimpleMolGraph)

Custom function for which actions are to be performed on a MolGraph when properties change.
In particular, we skip kekulization, as it will not always be possible on our substructures
(since they may have invalid molecular properties). Other than that, the function is copied 
from the default `smiles_on_update!` from 
https://github.com/mojaie/MolecularGraph.jl/blob/1c4498363381cdfd6162368f33d54d67dd3f1e04/src/smarts/base.jl#L66C1-L77C4.
Note that this is the function in the 0.21.1 release of MolecularGraph (which is what 
MolecularFingerprints uses as a dependency).
"""
function helper_custom_on_update!(mol::SimpleMolGraph)    
    # preprocessing
    default_atom_charge!(mol)
    default_bond_order!(mol)
    # kekulize!(mol)  # skip kekulization!

    # recalculate bottleneck descriptors
    sssr!(mol)
    apparent_valence!(mol)
    valence!(mol)
    lone_pair!(mol)
    is_ring_aromatic!(mol)
end

"""
    helper_custom_on_init!(mol::SimpleMolGraph)

Custom function for which actions are to be performed on a MolGraph on initialization.
In particular, as we skip kekulization in our custom on_update! function, we add
kekulization on initialization, so that it has been done once at least.
Other than that, the function is copied from the default `smiles_on_init!` from
https://github.com/mojaie/MolecularGraph.jl/blob/1c4498363381cdfd6162368f33d54d67dd3f1e04/src/smarts/base.jl#L61C1-L64C4.
Note that this is the function in the 0.21.1 release of MolecularGraph (which is what
MolecularFingerprints is using as a dependency), and newer releases include an additional
step `check_valence`.
"""
function helper_custom_on_init!(mol::SimpleMolGraph)
    stereocenter_from_smiles!(mol)
    stereobond_from_smiles!(mol)

    kekulize!(mol)  # add kekulization to init step
end

"""
    helper_copy_mol(mol::MolGraph)

Creates a copy of the given molecule, and disables automatic kekulization on update.

This serves two goals:
- User-given molecules are not modified. This is important, as hydrogens are removed before
    creating the MHFP fingerprint, which is a side effect the user may not want on their 
    given molecule. Furthermore, the molecule is equipped with custom on_init and on_update
    functions (see second point below), which is not something the user should have to 
    consider.
- Automatic kekulization of the molecule upon modification (such as creating a subgraph) is
    disabled. This is important as kekulization is not always successful on the 
    substructures we generate, as some of their properties are considered invalid.
    This is done as follows: we give the MolGraph object custom on_init and on_update 
    functions, where in the latter, the kekulization is not included (which it otherwise is 
    by default), while in the former, we added kekulization (so that it is at least 
    performed once, upon initialization)

# Remark
For some reason, MolecularGraph.jl has defined MolState{T, F1, F2} as a generic type struct, 
where F1 and F2 are the types of the two functions on_init and on_update, respectively.
This means that we cannot modify the fields on_init and on_update in-place, since once a 
molecule is initialized, F1 and F2 have fixed types, namely, the types 
`typeof(<the current on_init function>)` and analogously for `on_update`. Since the new 
functions have types like `typeof(<the new on_init function>)`, replacing in-place is not 
possible.
Instead, we create a new MolGraph object.
In MolecularGraph v0.22.0, this has been changed, and F1, F2 are no longer generic types.
If MolecularFingerprints is ported to use this version in the future, the function below 
could be simplified into `mol=copy(mol); mol.state.on_init = <new_on_init_function>` etc.
"""
function helper_copy_mol(mol::MolGraph)
    # new molecule will be of this type
    current_mol_type = typeof(mol)

    # state of newe molecule will be of type MolState{T, F1, F2}, where
    # F1, F2 are typeof(helper_custom_on_init!), typeof(helper_custom_on_update!) and
    # T is as determined in the following line:
    T = typeof(mol.state).parameters[1]

    mol = current_mol_type(
        mol.graph, # graph, vertex and edge properties all stay the same
        mol.vprops,
        mol.eprops,
        mol.gprops,
        MolState{  # mol.state is changed
            T, # first parameter of MolState unchanged
            typeof(helper_custom_on_init!),  # second and third type parameters changed
            typeof(helper_custom_on_update!)
            }(
            mol.state.initialized,  # MolState stays mostly the same
            mol.state.has_updates,
            mol.state.has_new_edges,
            mol.state.force_calculate,
            helper_custom_on_init!,  # add custom on_init and on_update functions
            helper_custom_on_update!)
    )
end

"""
    fingerprint(mol::MolGraph, calc::MHFP)

Calculates the MHFP fingerprint of the given molecule and returns it as a vector of UInt32's
"""
function fingerprint(mol::MolGraph, calc::MHFP)
    return mhfp_hash_from_molecular_shingling(  # calculate hash
        mhfp_shingling_from_mol(mol, calc),     # given the shingling of the molecule
        calc)                                   # using the parameters stored in calc
end

"""
    mhfp_shingling_from_mol(
        mol::MolGraph,
        calc::MHFP)

Calculate the "molecular shingling" of a given molecule.

A molecular shingling is a vector of "SMILES"-strings, calculated from the ring 
structures and atom types of the molecule (optional), and the circular substructures 
around each heavy (=non-hydrogen) atom of the molecule.

# Arguments
- `mol::MolGraph`: the molecule for which to calculate the shingling.
- `calc::MHFP`: fingeprint "calculator" object, containing the relevant parameters for the 
    fingerprint calculation, e.g., the radii of the circular substructures to be considered
    and whether to include ring information explicitly in the fingerprints
"""
function mhfp_shingling_from_mol(
    mol::MolGraph,
    calc::MHFP
)   

    # make copy of molecule, so modifications, such as removing Hs or changing the on_update
    # function of the MolGraph, are not affecting the given molecule
    mol = helper_copy_mol(mol)

    # read parameters from calculator object
    radius = calc.radius
    rings = calc.rings
    min_radius = calc.min_radius

    # remove all hydrogens in the molecule, as we only want to consider heavy atoms in all 
    # following steps
    remove_all_hydrogens!(mol)

    shingling::Vector{String} = []

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
    shingling_snippet::Vector{String} = []


    # Go through all rings in the sssr
    for ring in sssr(mol)
        # For each ring in the sssr, create smiles string of the submolecule corresponding
        # to the ring and add to the shingling.
        smiles_from_ring = smiles(induced_subgraph(mol, ring)[1])

        isnothing(smiles_from_ring) || push!(shingling_snippet, smiles_from_ring)
    end

    return shingling_snippet
end

"""
    smiles_from_atoms(mol::MolGraph)

Return vector containing SMILES strings of all atoms of the given molecule.
"""
function smiles_from_atoms(mol::MolGraph)
    shingling_snippet::Vector{String} = []

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
    shingling_snippet::Vector{String} = []

    # Ensure radius >= 1
    radius >=   0 || error("""radius must be strictly positive in this function.\nGot
        radius=$radius.\nTo generate the SMILES strings for individual atoms, call 
        the function smiles_from_atoms instead.""")
    # Ensure min_radius >= 1
    min_radius > 0 || error("""min_radius must be strictly positive in this function.\nGot
        min_radius=$min_radius.\nTo generate the SMILES strings for individual atoms, call 
        the function smiles_from_atoms instead.""")

    # # Ensure radius >= min_radius
    # min_radius ≤ radius || error("""radius must be larger or equal to min_radius.\nGot 
    # radius=$radius but min_radius=$min_radius.""")

    for atom_index in vertices(mol)  # go through all atoms
        
        for i = min_radius:radius  # go through al selected radii
            
            atoms_in_substructure_of_radius_i = neighborhood(mol, atom_index, i)

            submol, atom_map = induced_subgraph(mol, atoms_in_substructure_of_radius_i)
            
            # NOTE: This test is copied from the original authors.
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
            # need to subtract one, as the rootedAtAtom keyword is passed directly to C++ code, which uses 0-based indexing
            pos_of_atom_index_in_submol = findfirst(atom_map .== atom_index) - 1  
            
            smiles_of_substructure = smiles(
                submol,
                Dict{String,Any}("rootedAtAtom" => pos_of_atom_index_in_submol),
            )

            if !isnothing(smiles_of_substructure) && smiles_of_substructure != ""
                # Add smiles of substructure to shingling.
                push!(shingling_snippet, smiles_of_substructure)
            end

        end
    end

    return shingling_snippet
end

"""
    mhfp_hash_from_molecular_shingling(shingling::Vector{String}, calc::MHFP)

Calculate the MinHash values from a given Molecular shingling.

The given calculator contains parameters such as the length of the random vectors a , b that
are used in the hashing scheme, as well as the seed used when generating them.
The algorithm is described in more detail in the original authors paper.
"""
function mhfp_hash_from_molecular_shingling(shingling::Vector{String}, calc::MHFP)
    hash_values = zeros(UInt32, (calc.fp_size))
    fill!(hash_values, calc._max_hash)
    
    for s in shingling
        # create sha1 hash from the string
        sha_from_string = sha1(s)[begin:4]
        
        # Note: we are only using the first 4 sha1 bytes, as we want a 32-bit hash

        buf = IOBuffer(sha_from_string)  # make the sha bytes a buffer
        
        s_h = UInt32(  # read bytes from sha hash into integer
            htol(  # ensure little-endian format of the integer
                read(buf, UInt32))
                )  # read the buffer bytes as unsigned integer
        
        # apply equation 2 from the original authors paper
        hashes = Vector{UInt32}(mod.(
            mod.(
                calc._a * s_h + calc._b,
                calc._mersenne_prime
            ), calc._max_hash
        ))
        
        hash_values = min.(hash_values, hashes)
                
    end

    return hash_values
end
