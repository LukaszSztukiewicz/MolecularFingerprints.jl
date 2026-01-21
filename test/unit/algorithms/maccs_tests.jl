@testset "MACCS - Julia vs RDKit consistency" begin

    # smiles = "C=C(C)C"    
    # smiles = "COCC" 
    # smiles = "C=CN"  

    # smiles = "OS(=O)O"  
    # !!! Bit 81: julia[81] = 1 i rdkit[81] = 0

    smiles = "C1=COC=C1" 
    # !!! Bit 99: julia[99] = 1 i rdkit[99] = 0  
    # !!! Bit 157: julia[157] = 1 i rdkit[157] = 0  

    #smiles = "OS(=O)O"
    mol = MolecularGraph.smilestomol(smiles)

    fp_julia = MACCS(false, false)
    julia_fp = fingerprint(smiles, fp_julia)
    # julia_fp = fingerprint(mol, fp_julia)

    rdkit_fp = MolecularFingerprints.fingerprint_rdkit(smiles)

    @test length(julia_fp) == 166
    @test length(rdkit_fp) == 166

    @test julia_fp[58] == rdkit_fp[58]
    @debug "Bit 58" julia=julia_fp[58] rdkit=rdkit_fp[58]

    @test julia_fp[166] == rdkit_fp[166]
    @debug "Bit 166" julia=julia_fp[166] rdkit=rdkit_fp[166]

    @test julia_fp[1] == rdkit_fp[1]
    @debug "Bit 1" julia=julia_fp[1] rdkit=rdkit_fp[1]

    @debug("\n=== MACCS bit-by-bit comparison for SMILES: $smiles ===")

    for i in 1:166
        julia_fp[i] == 9 && continue

        if julia_fp[i] != rdkit_fp[i]
            @debug "Bit $i" julia=julia_fp[i] rdkit=rdkit_fp[i]
        end
    end

    @debug "=============================================="
end

