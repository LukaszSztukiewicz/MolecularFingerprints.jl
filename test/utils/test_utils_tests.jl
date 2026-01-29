@testset "Cosine Similarity Tests" begin
    
    # Test case 1: identical fingerprints
    a = Vector{Int}([1, 2, 3, 4, 5])
    b = Vector{Int}([1, 2, 3, 4, 5])
    @test isapprox(cosine_similarity(a, b), 1.0)

    # Test case 2: orthogonal fingerprints
    a = Vector{Int}([1, 0, 0])
    b = Vector{Int}([0, 1, 0])
    @test isapprox(cosine_similarity(a, b), 0.0)

    # Test case 3: partial overlap
    a = Vector{Int}([1, 2, 0, 0])
    b = Vector{Int}([0, 2, 3, 0])
    # cosine similarity = (2*2) / (sqrt(1^2 + 2^2) * sqrt(2^2 + 3^2)) = 4 / (sqrt(5) * sqrt(13)) = 4 / sqrt(65)
    expected_sim = 4 / sqrt(65)
    @test isapprox(cosine_similarity(a, b), expected_sim)

    # Test case 4: zero vector fingerprint
    a = Vector{Int}([0, 0, 0])
    b = Vector{Int}([1, 2, 3])
    @test isapprox(cosine_similarity(a, b), 0.0)

    # Test case 5: different lengths should throw an error
    a = Vector{Int}([1, 2, 3])
    b = Vector{Int}([1, 2])
    @test_throws ArgumentError cosine_similarity(a, b)

    # Test case 6: BitVector fingerprints
    a = BitVector([1, 0, 1, 1])
    b = BitVector([1, 1, 0, 1])
    # cosine similarity = (1*1 + 0*1 + 1*0 + 1*1) / (sqrt(1^2 + 0^2 + 1^2 + 1^2) * sqrt(1^2 + 1^2 + 0^2 + 1^2)) = 2 / (sqrt(3) * sqrt(3))
    expected_sim = 2 / (sqrt(3) * sqrt(3))
    @test isapprox(cosine_similarity(a, b), expected_sim)

    # Test case 7: SparseVector fingerprints
    a = sparsevec([1, 3, 5], [1, 2, 3], 5)
    b = sparsevec([2, 3, 4], [4, 2, 1], 5)
    # cosine similarity = (0*4 + 2*2 + 0*1) / (sqrt(1^2 + 2^2 + 3^2) * sqrt(4^2 + 2^2 + 1^2)) = 4 / (sqrt(14) * sqrt(21))
    expected_sim = 4 / (sqrt(14) * sqrt(21))
    @test isapprox(cosine_similarity(a, b), expected_sim)

    # Test case 8: both zero SparseVectors
    a = spzeros(Int, 5)
    b = spzeros(Int, 5)
    @test isapprox(cosine_similarity(a, b), 0.0)

    # Test case 9: SparseVector edge case with one zero vector
    a = sparsevec([1, 2], [1, 2], 5)
    b = spzeros(Int, 5)
    @test isapprox(cosine_similarity(a, b), 0.0)

end