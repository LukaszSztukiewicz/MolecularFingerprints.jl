using MolecularFingerprints
using Test

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

end
