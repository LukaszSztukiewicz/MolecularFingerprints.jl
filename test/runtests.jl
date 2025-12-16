using Test

@testset "MolecularFingerprints.jl" begin
    @testset "Unit Tests" begin
        include("unit/utils/tanimoto_tests.jl")
        include("unit/algorithms/ecfp_tests.jl")
        include("unit/algorithms/mhfp_tests.jl")
        include("unit/algorithms/maccs_tests.jl")
        include("unit/algorithms/torsions_tests.jl")
    end

    # @testset "RDKit Validation" begin
    #     include("validation/rdkit_comparison.jl")
    # end
end
