"""
    cosine_similarity(fp1::AbstractVector, fp2::AbstractVector) -> Float64
Calculate the cosine similarity between two fingerprint vectors. Using Distances.jl for efficient computation.
"""
function cosine_similarity(fp1::AbstractVector, fp2::AbstractVector)
    if length(fp1) != length(fp2)
        throw(ArgumentError("Fingerprints must be of the same length"))
    end
    
    #handle all zero vectors
    if allzeros(fp1) || allzeros(fp2)
        return 0.0
    end

    dist = cosine_dist(fp1, fp2)
    return 1.0 - dist
end