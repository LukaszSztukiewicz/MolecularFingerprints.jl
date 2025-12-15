module MACCS

using MolecularGraph
using SparseArrays
using Graphs

export MACCSFingerprint, fingerprint

include("../interface.jl")

abstract type AbstractFingerprint end
nbits(::AbstractFingerprint) = error("nbits not implemented")

struct MACCSFingerprint <: AbstractFingerprint
    count::Bool        # false = bit, true = count-based
    sparse::Bool       # false = dense, true = sparse
end

# 166 bits
nbits(::MACCSFingerprint) = 166

# check whether atom X is contained in molecule
function has_atom(mol, sym::Symbol)
    for v in vertices(mol.graph)
        atom = mol.vprops[v]
        if atom_symbol(atom) == sym
            return true
        end
    end
    return false
end

# count how many atoms are in molecule
function count_atom(mol, sym::Symbol)
    c = 0
    for v in vertices(mol.graph) # iteration over atoms in a graph of a molecule
        atom = mol.vprops[v]  # v = the number of the atom we are checking to see if it is the atom we are looking for, mol.vprops[v] object of the atom 
        if atom_symbol(atom) == sym
            c += 1
        end
    end
    return c
end

# check whether an bond exists between two atoms, order=1/2 single/double binding
function has_bond(mol, s1::Symbol, s2::Symbol, order::Int)
    for e in edges(mol.graph) # bindings in molecule
        v1 = src(e)
        v2 = dst(e)

        a1 = atom_symbol(mol.vprops[v1])    # atoms on both sides of the bond 
        a2 = atom_symbol(mol.vprops[v2])

        bond = mol.eprops[e] # information about bond

        if ((a1 == s1 && a2 == s2) ||
            (a1 == s2 && a2 == s1)) &&
        bond.order == order
            return true
        end
    end
    return false
end

# ------------------------------------------------------------------------------
# MACCS rules 
# ------------------------------------------------------------------------------
const MACCS_RULES = Dict{Int, Function}(
    # 1: contains nitrogen
    1 => mol -> count_atom(mol, :N),

    # 2: two or more oxygens
    2 => mol -> count_atom(mol, :O),

    # 3: carbonyl C=O
    3 => mol -> (has_bond(mol, :C, :O, 2) ? 1 : 0),
    
    # 4: contains sulfur
    4 => mol -> (has_atom(mol, :S) ? 1 : 0),

    # 5: contains halogen
    5 => mol -> (
        has_atom(mol, :F) ||
        has_atom(mol, :Cl) ||
        has_atom(mol, :Br) ||
        has_atom(mol, :I)
    ) ? 1 : 0
)

function compute_maccs(mol, fp::MACCSFingerprint)
    vec = zeros(Int, nbits(fp)) # [0, 0, 0, 0, 0, ..., 0]  (166 elemetns)

    for (idx, rule) in MACCS_RULES
        # hit = rule(mol)
        # vec[idx] = fp.count ? Int(hit) : (hit ? 1 : 0)
        val = rule(mol)
        vec[idx] = fp.count ? val : (val > 0 ? 1 : 0)
    end

    return fp.sparse ? sparse(vec) : vec
    # fingerprint looks often like [0,0,1,0,0,0,0,1,0,0,0,...] - mostly zeros → waste of memory - use sparse vector
end

# get molecule from SMILES
function mol_from_smiles(smiles::AbstractString)
    return MolecularGraph.smilestomol(smiles)
end

function fingerprint(fp::MACCSFingerprint, smiles::AbstractString)
    mol = mol_from_smiles(smiles)  # molecular structure with atoms, bonds
    return compute_maccs(mol, fp)
end

end

# fp = MACCSFingerprint(false, false)
# #fp = MACCSFingerprint(true, false)
# println("Test: acetic acid CC(=O)O")
# v = fingerprint(fp, "CC(=O)O")
# # smiles CC(=O)O - molecule in a text form
# println("Vector: ", v)
# println("Bits set: ", findall(!=(0), v))




# using MolecularGraph

# include("../interface.jl")

# struct MACCS{N} <: AbstractFingerprint
#     radius::Int
# end

# function fingerprint(mol::MolecularGraph.Mol, calc::MACCS{N}) where N
#     # Placeholder implementation for MACCS fingerprint calculation
#     # In a real implementation, this would compute the MACCS fingerprint
#     # based on the molecular structure and the specified radius.
#     return BitVector(rand(Bool, 1024))  # Example: return a random 1024-bit fingerprint
# end