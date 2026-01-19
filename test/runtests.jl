using Test
using MolecularGraph
using MolecularFingerprints
using Random
using SHA
using Graphs: vertices, induced_subgraph, neighborhood

using Test
using MolecularGraph
using MolecularFingerprints

using PythonCall: Py, pyimport, pyconvert, pybuiltins
using Graphs: nv

using Test
using PythonCall
using MolecularGraph
using MolecularFingerprints
using SparseArrays

using Test
using MolecularGraph
using MolecularFingerprints

@testset "MolecularFingerprints.jl" begin
    @testset "Unit Tests" begin
        include("unit/utils/tanimoto_tests.jl")
        include("unit/algorithms/ecfp_tests.jl")
        include("unit/algorithms/mhfp_tests.jl")
        include("unit/algorithms/maccs_tests.jl")
        include("unit/algorithms/torsions_tests.jl")
    end
end
