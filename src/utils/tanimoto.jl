"""
    tanimoto(a::BitVector, b::BitVector)

Calculate the Tanimoto similarity coefficient (Jaccard Index) between two fingerprints.
Formula: c / (a + b - c) where c is intersection count.
"""
function tanimoto(a::BitVector, b::BitVector)
    @assert length(a) == length(b) "Fingerprints must be of the same length"
    
    # High-performance bit counting
    intersection = count(a .& b)
    union_sum = count(a) + count(b)
    
    if union_sum == 0
        return 0.0 # Handle edge case of two empty fingerprints
    end
    
    return intersection / (union_sum - intersection)
end