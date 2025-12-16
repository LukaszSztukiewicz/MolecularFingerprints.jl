using MolecularFingerprints
using Test

@testset "MolecularFingerprints.jl" begin
    @testset "Unit Tests" begin
        include("unit/utils/tanimoto_tests.jl")
        include("unit/algorithms/ecfp_tests.jl")
        include("unit/algorithms/mhfp_tests.jl")
    end
end
