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
    if length(a) != length(b)
        throw(ArgumentError("Fingerprints must be of the same length"))
    end

    set_a = Set(a)
    set_b = Set(b)
    intersection_sets = length(intersect(set_a, set_b))
    union_sets = length(union(set_a, set_b))
    if union_sets == 0
        return 0.0
    end
    return intersection_sets / union_sets
end

function tanimoto_similarity(fp1::SparseVector, fp2::SparseVector)
    # Get indices where either vector has a non-zero value
    common_indices = union(fp1.nzind, fp2.nzind)
    
    numerator = 0
    denominator = 0
    
    for i in common_indices
        v1 = fp1[i]
        v2 = fp2[i]
        numerator += min(v1, v2)
        denominator += max(v1, v2)
    end
    
    return denominator == 0 ? 0.0 : numerator / denominator
end