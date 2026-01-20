using Test
using MolecularGraph
using MolecularGraph: smilestomol, MolGraph
using MolecularFingerprints
using Random
using SHA
using SparseArrays
using Graphs: nv, all_simple_paths, degree
using Graphs: vertices, induced_subgraph, neighborhood
using PythonCall: Py, pyimport, pyconvert, pybuiltins


@testset "MolecularFingerprints.jl" begin

    @testset "1. Unit Tests" begin

        @testset "1.1 Interface Tests" begin
            @info "Running Interface Tests..."
            include("unit/interface/mock_interface_tests.jl")
        end

        @testset "1.2 Utility Function Tests" begin
            @info "Running Utility Function Tests..."
            include("unit/utils/tanimoto_tests.jl")
        end

        @testset "1.3 Fingerprint Algorithm Tests" begin
            @info "Running Fingerprint Algorithm Tests..."

            @testset "1.3.1 ECFP Tests" begin
                @info "Running ECFP Tests..."
                include("unit/algorithms/ecfp_tests.jl")
            end

            @testset "1.3.2 MHFP Tests" begin
                @info "Running MHFP Tests..."
                include("unit/algorithms/mhfp_tests.jl")
            end

            @testset "1.3.3 MACCS Tests" begin
                @info "Running MACCS Tests..."
                include("unit/algorithms/maccs_tests.jl")
            end

            @testset "1.3.4 Topological Torsion Tests" begin
                @info "Running Topological Torsion Tests..."
                include("unit/algorithms/torsions_tests.jl")
            end
        end
    end

    @testset "2. Reference Tests (Validation)" begin
        @info "Running Validation Reference Tests..."
        include("validation/rdkit_comparison_tests.jl")
    end
end
