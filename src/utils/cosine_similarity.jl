"""
    cosine_similarity(fp1::Vector{T}, fp2::Vector{T}) where T<:Integer

Calculate the cosine similarity between two integer vector fingerprints.
Formula: (A · B) / (||A|| * ||B||)
"""
function cosine_similarity(fp1::Vector{T}, fp2::Vector{T}) where T<:Integer
    if length(fp1) != length(fp2)
        throw(ArgumentError("Fingerprints must be of the same length"))
    end

    maxval1 = maximum(abs, fp1)
    maxval2 = maximum(abs, fp2)

    if maxval1 == 0 || maxval2 == 0  # Handle edge case of zero-vector fingerprints
        return 0.0
    end

    # normalize vectors to avoid overflow
    fp1 /= maxval1
    fp2 /= maxval2

    dot_product = sum(x * y for (x, y) in zip(fp1, fp2))
    norm_fp1 = sqrt(sum(x^2 for x in fp1))
    norm_fp2 = sqrt(sum(y^2 for y in fp2))
    
    return dot_product / (norm_fp1 * norm_fp2)
end

"""
    cosine_similarity(fp1::BitVector, fp2::BitVector)
Calculate the cosine similarity between two BitVector fingerprints.
Formula: (A · B) / (||A|| * ||B||)
"""
function cosine_similarity(fp1::BitVector, fp2::BitVector)
    if length(fp1) != length(fp2)
        throw(ArgumentError("Fingerprints must be of the same length"))
    end
    dot_product = count(fp1 .& fp2)
    norm_fp1 = sqrt(count(fp1))
    norm_fp2 = sqrt(count(fp2))
    if norm_fp1 == 0.0 || norm_fp2 == 0.0
        return 0.0 # Handle edge case of zero-vector fingerprints
    end
    return dot_product / (norm_fp1 * norm_fp2)
end

function cosine_similarity(fp1::SparseVector, fp2::SparseVector)
    # Get indices where either vector has a non-zero value
    common_indices = union(fp1.nzind, fp2.nzind)
    
    numerator = 0.0
    sum_sq_fp1 = 0.0
    sum_sq_fp2 = 0.0
    
    for i in common_indices
        v1 = fp1[i]
        v2 = fp2[i]
        numerator += v1 * v2
        sum_sq_fp1 += v1^2
        sum_sq_fp2 += v2^2
    end
    
    denominator = sqrt(sum_sq_fp1) * sqrt(sum_sq_fp2)
    return denominator == 0.0 ? 0.0 : numerator / denominator
end