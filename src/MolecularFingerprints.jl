module MolecularFingerprints

using MolecularGraph
using Random
using Graphs
using RDKitMinimalLib
using SHA


using MolecularGraph: AbstractMolGraph, edge_rank

using PythonCall: Py, pyimport, pyconvert
using MolecularGraph: SMILESMolGraph, smilestomol
using Graphs: vertices, edges, neighbors, src, dst, degree, cycle_basis
using SparseArrays: sparse


using MolecularGraph
using SHA
using Graphs: induced_subgraph, nv, vertices, all_simple_paths
using SparseArrays

"""
    MolecularFingerprints

    A Julia package for computing various molecular fingerprints used in cheminformatics.

    # Modules Included
    - Abstract Interfaces
    - Utility Functions
    - Fingerprint Algorithms

    # Exported Functions
    - `tanimoto`
"""

# In julia inclusion of files make the further iclusions to see previously defined symbols,
# so the order of includes matters.

# Abstract Interfaces
include("interface.jl")

# Utility Functions
include("utils/tanimoto.jl")

# Fingerprint Algorithms
include("algorithms/mhfp.jl")
include("algorithms/ecfp.jl")
include("algorithms/maccs.jl")
include("algorithms/torsions.jl")

#FIXME in the end, probably we shouldn't export mhfp_shingling_from_mol
export tanimoto, mhfp_shingling_from_mol, fingerprint
# FIXME remove mhfp_hash_from_molecular_shingling and mhfp_shingling_from_mol from export
# later (it is useful right now for (manual))
export MHFP, mhfp_shingling_from_mol, fingerprint, mhfp_hash_from_molecular_shingling 
export ecfp_atom_invariant, ecfp_hash, ECFP, fingerprint
export MACCSFingerprint, fingerprint, fingerprint_rdkit
export TopologicalTorsion, fingerprint, getAtomCode, get4paths, getTopologicalTorsionFP, getPathsOfLengthN # FIXME remove unnecessary exports

end
