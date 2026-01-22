"""
    cosine_similarity(fp1::Vector{T}, fp2::Vector{T}) where T<:Integer

Calculate the cosine similarity between two integer vector fingerprints.
Formula: (A · B) / (||A|| * ||B||)
"""
function cosine_similarity(fp1::Vector{T}, fp2::Vector{T}) where T<:Integer
    if length(fp1) != length(fp2)
        throw(ArgumentError("Fingerprints must be of the same length"))
    end

    dot_product = sum(x * y for (x, y) in zip(fp1, fp2))
    norm_fp1 = sqrt(sum(x^2 for x in fp1))
    norm_fp2 = sqrt(sum(y^2 for y in fp2))
    if norm_fp1 == 0.0 || norm_fp2 == 0.0
        return 0.0 # Handle edge case of zero-vector fingerprints
    end
    return dot_product / (norm_fp1 * norm_fp2)
end