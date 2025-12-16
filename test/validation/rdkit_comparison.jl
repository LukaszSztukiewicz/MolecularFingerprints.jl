using Statistics: cor # For correlation
using Random
using ReferenceTests
using MolecularGraph
using RDKitMinimalLib
using MolecularFingerprints


function tanimoto_rdkit(fp1_str::AbstractString, fp2_str::AbstractString)
    # RDKitMinimalLib returns a string of '0' and '1'
    # Convert to BitVector: '1' becomes true, '0' becomes false
    b1 = collect(fp1_str) .== '1'
    b2 = collect(fp2_str) .== '1'
    return MolecularFingerprints.tanimoto(b1, b2)
end

@testset "ECFP4 RDKit Correlation" begin
    smiles = [
        "CC(=O)Oc1ccccc1C(=O)O", # Aspirin
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O", # Ibuprofen
        "Cc1ccccc1", # Toluene
        "c1ccccc1", # Benzene
        "CCO", # Ethanol
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" # Caffeine
    ]

    mols_julia = [smilestomol(s) for s in smiles]
    mols_rdkit = [get_mol(s) for s in smiles]

    rdkit_sims = Float64[]
    julia_sims = Float64[]

    n = length(smiles)
    for i in 1:n
        for j in (i+1):n
            fp_i_rd_str = get_morgan_fp(mols_rdkit[i]; radius=2, nBits=2048)
            fp_j_rd_str = get_morgan_fp(mols_rdkit[j]; radius=2, nBits=2048)
            push!(rdkit_sims, tanimoto_rdkit(fp_i_rd_str, fp_j_rd_str))

            fp_i_jl = fingerprint(mols_julia[i], ECFP{2}(2))
            fp_j_jl = fingerprint(mols_julia[j], ECFP{2}(2))
            push!(julia_sims, tanimoto(fp_i_jl, fp_j_jl))
        end
    end

    # 3. Statistical Validation
    # We check if the similarity scores are correlated.
    # A high correlation (> 0.85) means we are capturing the same chemical features.
    correlation = cor(rdkit_sims, julia_sims)
    
    println("ECFP4 Correlation with RDKit: $correlation")
    
    # Assert that our implementation broadly agrees with RDKit
    @test correlation > 0.85
end

@testset "Consistency (ReferenceTests)" begin
    Random.seed!(42)
    # Regression testing: Ensure the hashing doesn't drift in future versions
    mol = smilestomol("CC(=O)Oc1ccccc1C(=O)O") # Aspirin
    fp = fingerprint(mol, ECFP{2}(2))
    @test_reference "references/aspirin_ecfp4.txt" BitVector(fp)
end