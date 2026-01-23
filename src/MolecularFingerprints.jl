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
    SMILESBond,
    atom_charge,
    atom_number,
    atom_symbol,
    edge_rank,
    exact_mass,
    explicit_hydrogens,
    implicit_hydrogens,
    is_aromatic,
    is_in_ring,
    monoiso_mass,
    pi_electron,
    remove_all_hydrogens!,
    smilestomol,
    sssr,
    subset,
    valence
using Random: rand, randstring, seed!
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

export AbstractCalculator, AbstractFingerprint, AbstractDescriptor
export MACCS, TopologicalTorsion, MHFP, ECFP
export tanimoto_similarity, cosine_similarity
export fingerprint

end
