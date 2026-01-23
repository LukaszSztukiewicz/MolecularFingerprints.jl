using MolecularFingerprints
using Test

@testset "tanimoto_similarity Similarity Tests" begin
    
    # edge case: all bits set to zero
    a = falses(10)
    b = falses(10)
    @test isapprox(tanimoto_similarity(a, b), 0.0)

    # full overlap
    a = BitVector([1, 0, 1, 1, 0])
    b = BitVector([1, 0, 1, 1, 0])
    @test isapprox(tanimoto_similarity(a, b), 1.0)

    # no overlap
    a = BitVector([1, 0, 0, 0, 1])
    b = BitVector([0, 1, 1, 1, 0])
    @test isapprox(tanimoto_similarity(a, b), 0.0)


    # partial overlap
    # A: 1 1 0 1 0
    # B: 1 0 1 1 0
    # intersection = 2 (positions 1 and 4)
    # union = 4 (positions 1,2,3,4)
    # tanimoto_similarity = 2 / 4 = 0.5
    a = BitVector([1, 1, 0, 1, 0])
    b = BitVector([1, 0, 1, 1, 0])
    @test isapprox(tanimoto_similarity(a, b), 0.5)

    # partial overlap 2
    # A: 1 0 1 1 1 0
    # B: 1 1 0 1 0 1
    # intersection = 2 (positions 1 and 4)
    # union = 6 (positions 1,2,3,4,5,6)
    # tanimoto_similarity = 2 / 6 = 0.3333...
    a = BitVector([1, 0, 1, 1, 1, 0])
    b = BitVector([1, 1, 0, 1, 0, 1])
    @test isapprox(tanimoto_similarity(a, b), 1/3)

    # different lengths should throw an error
    a = BitVector([1, 0, 1])
    b = BitVector([1, 0])
    @test_throws ArgumentError tanimoto_similarity(a, b)

    # larger fingerprints
    a = BitVector(rand(Bool, 1024))
    b = BitVector(rand(Bool, 1024))
    sim = tanimoto_similarity(a, b)
    @test 0.0 <= sim <= 1.0

    # identical large fingerprints
    a = BitVector(rand(Bool, 2048))
    b = copy(a)
    @test isapprox(tanimoto_similarity(a, b), 1.0)

    # different lengths should throw an error
    a = [1, 2, 3]
    b = [1, 2]
    @test_throws ArgumentError tanimoto_similarity(a, b)

    # integer vector fingerprints
    a = [1, 2, 3, 4]
    b = [3, 4, 5, 6]
    # intersection = {3,4} -> 2
    # union = {1,2,3,4,5,6} -> 6
    # tanimoto_similarity = 2 / 6 = 0.3333...
    @test isapprox(tanimoto_similarity(a, b), 1/3)

    # SparseVector fingerprints
    a = sparsevec([1, 3, 5], [1, 2, 3], 5)
    b = sparsevec([2, 3, 4], [4, 5, 6], 5)
    # intersection = min(0,4) + min(2,5) + min(0,6) = 2
    # union = max(1,0) + max(0,4) + max(2,5) + max(0,6) + max(3,0) = 1 + 4 + 5 + 6 + 3 = 19
    # tanimoto_similarity = 2 / 19
    @test isapprox(tanimoto_similarity(a, b), 2 / 19)

end
