module MolecularFingerprints

using Graphs:
    all_simple_paths,
    cycle_basis,
    degree,
    dst,
    edges,
    fadjlist,
    induced_subgraph,
    ne,
    neighborhood,
    neighbors,
    nv,
    src,
    vertices
using MolecularGraph:
    AbstractMolGraph,
    MolGraph,
    MolState,
    SimpleMolGraph,
    SMILESBond,
    apparent_valence!,
    atom_charge,
    atom_number,
    atom_symbol,
    default_atom_charge!,
    default_bond_order!,
    edge_rank,
    exact_mass,
    explicit_hydrogens,
    get_prop,
    implicit_hydrogens,
    is_aromatic,
    is_in_ring,
    is_ring_aromatic!,
    lone_pair!,
    monoiso_mass,
    pi_electron,
    remove_all_hydrogens!,
    set_prop!,
    smiles,
    smilestomol,
    sssr,
    sssr!,
    subset,
    valence,
    valence!
using Random: rand, randstring, seed!
using RDKitMinimalLib: smiles
using SHA: sha1
using SparseArrays: sparse, spzeros, SparseVector

# NOTE: In Julia, order of includes matters for dependencies

# Abstract Interfaces
include("interface.jl")

# Utility Functions
include("utils/tanimoto_similarity.jl")
include("utils/cosine_similarity.jl")

# Fingerprint Algorithms
include("algorithms/mhfp.jl")
include("algorithms/ecfp.jl")
include("algorithms/maccs.jl")
include("algorithms/torsions.jl")

export safe_smilestomol, smiles_to_neutralized_mol

export AbstractCalculator, AbstractFingerprint, AbstractDescriptor
export MACCS, TopologicalTorsion, MHFP, ECFP
export tanimoto_similarity, cosine_similarity
export fingerprint

end
