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
