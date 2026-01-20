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

# Set seed for reproducibility across all tests
Random.seed!(42)

@testset "MolecularFingerprints.jl" begin

    @testset "Unit Tests" begin
        
        @testset "Interface" begin
            @info "Running Interface Tests..."
            include("unit/interface/mock_interface_tests.jl")
        end

        @testset "Utilities" begin
            @info "Running Utility Function Tests..."
            include("unit/utils/tanimoto_tests.jl")
        end

        @testset "Algorithms" begin
            # Using a loop for repetitive algorithm test sets
            algorithms = [
                ("ECFP", "unit/algorithms/ecfp_tests.jl"),
                ("MHFP", "unit/algorithms/mhfp_tests.jl"),
                ("MACCS", "unit/algorithms/maccs_tests.jl"),
                ("Topological Torsion", "unit/algorithms/torsions_tests.jl")
            ]

            for (name, path) in algorithms
                @testset "$name" begin
                    @info "$(findfirst(x -> x[1] == name, algorithms)). Running $name Tests..."
                    include(path)
                end
            end
        end
    end

    @testset "Reference Validation (RDKit)" begin
        @info "Running Reference Validation Tests..."
        include("validation/rdkit_comparison_tests.jl")
    end

end