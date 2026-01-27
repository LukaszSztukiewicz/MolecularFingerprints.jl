""" 
    cosine_similarity(fp1::Vector{T}, fp2::Vector{T}) where T<:Integer
Calculate the cosine similarity between two fingerprints represented as integer vectors.
# Arguments
- `fp1`: First fingerprint as a vector of integers.
- `fp2`: Second fingerprint as a vector of integers.
# Returns
- Cosine similarity value between 0.0 and 1.0.
"""
function cosine_similarity(fp1::Vector{T}, fp2::Vector{T}) where T<:Integer
    length(fp1) != length(fp2) && throw(ArgumentError("Fingerprints must be of the same length"))
    
    if all(iszero, fp1) || all(iszero, fp2)
        return 0.0
    end

    # cast to Float64 to prevent integer overflow during internal squaring/summing
    # we use cosine_dist from Distances.jl which computes 1 - cosine similarity
    calculated_dist = cosine_dist(Float64.(fp1), Float64.(fp2))
    
    sim =  isnan(calculated_dist) ? 0.0 : 1.0 - calculated_dist
    return clamp(sim, 0.0, 1.0)
end

function cosine_similarity(fp1::BitVector, fp2::BitVector)
    return cosine_similarity(Vector{Int}(fp1), Vector{Int}(fp2))
end