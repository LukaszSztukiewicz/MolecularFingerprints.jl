using Test
using MolecularGraph
using MolecularFingerprints
using Random
using SHA
using SparseArrays
using Graphs: nv, all_simple_paths, degree
using Graphs: vertices, induced_subgraph, neighborhood
using PythonCall: Py, pyimport, pyconvert, pybuiltins


@testset "MolecularFingerprints.jl" begin

    @testset "Unit Tests" begin

        @testset "Interface Tests" begin
            @info "Running Interface Tests..."
            include("unit/interface/mock_interface_tests.jl")
        end

        @testset "Utility Function Tests" begin
            @info "Running Utility Function Tests..."
            include("unit/utils/tanimoto_tests.jl")
        end

        @testset "Fingerprint Algorithm Tests" begin
            @info "Running Fingerprint Algorithm Tests..."

            @testset "ECFP Tests" begin
                @info "Running ECFP Tests..."
                include("unit/algorithms/ecfp_tests.jl")
            end

            @testset "MHFP Tests" begin
                @info "Running MHFP Tests..."
                include("unit/algorithms/mhfp_tests.jl")
            end

            @testset "MACCS Tests" begin
                @info "Running MACCS Tests..."
                include("unit/algorithms/maccs_tests.jl")
            end

            @testset "Topological Torsion Tests" begin
                @info "Running Topological Torsion Tests..."
                include("unit/algorithms/torsions_tests.jl")
            end
        end
    end

    @testset "Reference Tests (Validation)" begin
        @info "Running Validation Reference Tests..."
        include("validation/rdkit_comparison_tests.jl")
    end
end
