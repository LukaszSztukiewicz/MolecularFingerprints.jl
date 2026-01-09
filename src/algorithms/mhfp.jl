using MolecularGraph
using Random
using Graphs
using RDKitMinimalLib

"""
    MHFP{N}

Class for MHFP fingerprint generator/featurizer. Contains settings and parameters for 
MHFP fingerprint generation.
"""
struct MHFP{N} <: AbstractFingerprint
    radius::Int
end


"""
    fingerprint(mol::SMILESMolGraph, calc::MHFP{N})

Calculates the MHFP fingerprint of the given molecule and returns it as a bit vector
"""
function fingerprint(mol::SMILESMolGraph, calc::MHFP{N}) where {N}
    # Placeholder implementation for MHFP fingerprint calculation
    # In a real implementation, this would compute the MHFP fingerprint
    # based on the circular substructures around each atom of the molecule,
    # depending on the specified radius.
    return BitVector(rand(Bool, 1024))
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
around each atom of the molecule.

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

    shingling = []

    # Consider rings of the molecule, if corresponding parameter is set
    if rings
        append!(shingling, smiles_from_rings(mol))
    end

    # If min_radius == 0, add SMILES string of all atoms to shingling
    if min_radius == 0
        append!(shingling, smiles_from_atoms(mol))
        min_radius += 1  # Increase min_radius, as the case of radius=0 has now been dealt with.
    end


    # Add SMILES strings of all circular substructures around all heavy atoms to the shingling
    append!(shingling, smiles_from_circular_substructures(mol, radius, min_radius))

    # Make shingling entries unique
    sort!(shingling)   # unique! is much faster on sorted arrays (see https://docs.julialang.org/en/v1/base/collections/)
    unique!(shingling)

    return shingling
end

"""
    smiles_from_rings(mol::MolGraph)

Return vector containing SMILES strings of all rings in the SSSR of the given molecule.

SSSR stands for the smallest set of smallest rings of the molecule.
"""
function smiles_from_rings(mol::MolGraph)
    shingling_snippet = []

    # NOTE:
    # In the original implementation by the authors, they proceed like this:
    # Go through all rings, go through all pairs of contained vertices, for pairs that have a bond, add the bond to "bonds".
    # Then, create a submolecule from these bonds using AllChem.PathToSubmol(mol, list_of_bonds) from rdkit,
    # and finally, they write the smiles string for that submolecule with AllChem.MolToSmiles.

    # We instead directly take the induced subgraph of the molecule (aka the molecular graph), induced by the atoms in the sssr. Then we write the
    # smiles string of that using smiles() (which comes from rdkitminimallib.jl, but is supplied by MolecularGraph.)

    # TODO: verify that these approaches have the same result.

    for ring in MolecularGraph.sssr(mol)  # TODO: verify if sssr from MolecularGraph.jl is the "symmetrisized" one from rdkit. Probably it isn't. Discuss/decide if this is a problem.
        # For each ring in the sssr, create smiles string of the submolecule corresponding to the ring and add to the shingling.
        push!(shingling_snippet, smiles(induced_subgraph(mol, ring)[1]))
    end

    return shingling_snippet
end

"""
    smiles_from_atoms(mol::MolGraph)

Return vector containing SMILES strings of all heavy atoms of the given molecule.

Heavy atoms means all non-hydrogen atoms.
"""
function smiles_from_atoms(mol::MolGraph)
    shingling_snippet = []

    for atom in vertices(mol)
        atom_as_mol = induced_subgraph(mol, [atom])[1]  # create "molecule" containing only the atom
        # remove_all_hydrogens!(atom_as_mol)  # remove all hydrogens

        # TODO: the smiles function below returns capital C for the carbon in biphenyl (c1ccc(-c2ccccc2)cc1), which are supposed to be lower-case. rdkit in python returns lower-case.
        # This very likely comes from the fact that in python, we can iterate through the atom objects and extract smarts strings/smiles strings, while here, I create molecules containing only the atoms, which then lose there aromaticity information (which is what decides between c and C)
        # Worst case, use mol[atom].is_aromatic to shift manually.
        push!(shingling_snippet, smiles(atom_as_mol))  # Note: The original authors add the SMARTS string (not SMILES string) of the atoms to the shingling. However, in the rdkit implementation, they use SMILES strings, see https://github.com/rdkit/rdkit/blob/da08e8d954c0a923300a43d24faa37ad224b320d/Code/GraphMol/Fingerprints/MHFP.cpp#L124
    end

    return shingling_snippet
end

"""
    smiles_from_circular_substructures(
        mol::MolGraph,
        radius::Int,
        min_radius::Int)

Return vector of SMILES strings of circular substructures around all atoms of a molecule.

For each heavy (i.e., non-hydrogen) atom of the given molecule, extract the substructures
of radii min_radius to radius, and generate their corresponding SMILES strings.
"""
function smiles_from_circular_substructures(mol::MolGraph, radius::Int, min_radius::Int)
    shingling_snippet = []

    for atom_index in vertices(mol)  # go through all (heavy) atoms
        
        for i = min_radius:radius  # go through al selected radii
            
            atoms_in_substructure_of_radius_i = neighborhood(mol, atom_index, i)

            submol, atom_map = induced_subgraph(mol, atoms_in_substructure_of_radius_i)

            # TODO: This test is copied from the original authors.
            # I don't know what this test is for, as I don't see why it could be that 
            # atom_index is not contained in the atom map.
            # I guess this is just to make sure that we don't get an error but simply continue. 
            # Maybe exclude this test for speedup in the future.
            if atom_index ∉ atom_map
                continue
                @warn "Current atom was not found in the atom map, skipping" atom_index atom_map
            end
            
            # Find index of current atom (atom_index) in the submolecule
            pos_of_atom_index_in_submol = findfirst(atom_map .== atom_index)

            smiles_of_substructure = smiles(
                submol,
                Dict{String,Any}("rootedAtAtom" => pos_of_atom_index_in_submol),
            )

            if smiles_of_substructure != ""
                push!(shingling_snippet, smiles_of_substructure)  # Add smiles of substructure to shingling.
            end

        end
    end

    return shingling_snippet
end

export MHFP, mhfp_shingling_from_mol, fingerprint