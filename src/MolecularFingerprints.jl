module MolecularFingerprints

using MolecularGraph
using Random
using Graphs
using RDKitMinimalLib
using SHA
using MolecularGraph
using SparseArrays

using PythonCall: Py, pyimport, pyconvert
using MolecularGraph: SMILESMolGraph, smilestomol, AbstractMolGraph, edge_rank
using Graphs: vertices, edges, neighbors, src, dst, degree, cycle_basis
using Graphs: induced_subgraph, nv, vertices, all_simple_paths
using SparseArrays: sparse

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

export AbstractCalculator, AbstractFingerprint, AbstractDescriptor
export MACCS, TopologicalTorsion, MHFP, ECFP
export tanimoto
export fingerprint

end
