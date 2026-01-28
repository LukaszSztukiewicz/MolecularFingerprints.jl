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

function cosine_similarity(fp1::SparseVector, fp2::SparseVector)
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

"""
    rdkit_to_julia_sparse(rdkit_sparse_vect::PyObject, vector_length::Union{Nothing, Int}=nothing)
Convert an RDKit sparse fingerprint (Python object) to a native Julia `SparseVector`.
# Arguments
- `rdkit_sparse_vect`: The RDKit sparse fingerprint object.
- `vector_length`: Optional length of the resulting vector. If not provided,
    it will be determined from the RDKit object.
# Returns
- A `SparseVector` representing the fingerprint in Julia.
"""
function rdkit_to_julia_sparse(rdkit_sparse_vect, vector_length=nothing)

    py_dict = rdkit_sparse_vect.GetNonzeroElements()
    
    idxs = Int64[]
    vals = Int64[]
    
    for (k, v) in py_dict.items()
        push!(idxs, pyconvert(Int64, k))
        push!(vals, pyconvert(Int64, v))
    end

    len = isnothing(vector_length) ? pyconvert(Int, rdkit_sparse_vect.GetLength()) : vector_length

    return sparsevec(idxs, vals, len)
end
