using MolecularFingerprints
using Test

@testset "Tanimoto Similarity Tests" begin
    
    # edge case: all bits set to zero
    a = falses(10)
    b = falses(10)
    @test isapprox(tanimoto(a, b), 0.0)


end
