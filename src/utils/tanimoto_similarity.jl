"""
    tanimoto_similarity(a::BitVector, b::BitVector)

Calculate the tanimoto_similarity similarity coefficient (Jaccard Index) between two fingerprints.
Formula: c / (a + b - c) where c is intersection count.
"""
function tanimoto_similarity(a::BitVector, b::BitVector)
    if length(a) != length(b)
        throw(ArgumentError("Fingerprints must be of the same length"))
    end
    
    # High-performance bit counting
    intersection = count(a .& b)
    union_sum = count(a) + count(b)
    
    if union_sum == 0
        return 0.0 # Handle edge case of two empty fingerprints
    end
    
    return intersection / (union_sum - intersection)
end

function tanimoto_similarity(a::Vector{Int}, b::Vector{Int})
    throw(DomainError("You seem to be using integer vector fingerprints. Please use cosine_similarity function for integer vectors."))
end