using MolecularFingerprints
using Test

@testset "Tanimoto Similarity Tests" begin
    
    # edge case: all bits set to zero
    a = falses(10)
    b = falses(10)
    @test isapprox(tanimoto(a, b), 0.0)

    # full overlap
    a = BitVector([1, 0, 1, 1, 0])
    b = BitVector([1, 0, 1, 1, 0])
    @test isapprox(tanimoto(a, b), 1.0)

    # no overlap
    a = BitVector([1, 0, 0, 0, 1])
    b = BitVector([0, 1, 1, 1, 0])
    @test isapprox(tanimoto(a, b), 0.0)


    # partial overlap
    # A: 1 1 0 1 0
    # B: 1 0 1 1 0
    # intersection = 2 (positions 1 and 4)
    # union = 4 (positions 1,2,3,4)
    # Tanimoto = 2 / 4 = 0.5
    a = BitVector([1, 1, 0, 1, 0])
    b = BitVector([1, 0, 1, 1, 0])
    @test isapprox(tanimoto(a, b), 0.5)

    # partial overlap 2
    # A: 1 0 1 1 1 0
    # B: 1 1 0 1 0 1
    # intersection = 2 (positions 1 and 4)
    # union = 6 (positions 1,2,3,4,5,6)
    # Tanimoto = 2 / 6 = 0.3333...
    a = BitVector([1, 0, 1, 1, 1, 0])
    b = BitVector([1, 1, 0, 1, 0, 1])
    @test isapprox(tanimoto(a, b), 1/3)

    # different lengths should throw an error
    a = BitVector([1, 0, 1])
    b = BitVector([1, 0])
    @test_throws AssertionError tanimoto(a, b)

    # larger fingerprints
    a = BitVector(rand(Bool, 1024))
    b = BitVector(rand(Bool, 1024))
    sim = tanimoto(a, b)
    @test 0.0 <= sim <= 1.0

    # identical large fingerprints
    a = BitVector(rand(Bool, 2048))
    b = copy(a)
    @test isapprox(tanimoto(a, b), 1.0)



end
