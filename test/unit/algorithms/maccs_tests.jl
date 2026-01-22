@testset "MACCS - Julia vs RDKit consistency" begin

    _Chem  = Ref{Py}()
    _MACCS = Ref{Py}()

    function __init__()
        _Chem[]  = pyimport("rdkit.Chem")
        _MACCS[] = pyimport("rdkit.Chem.MACCSkeys")
    end

    function fingerprint_rdkit(smiles::AbstractString)
        mol = _Chem[].MolFromSmiles(smiles)

        fp = _MACCS[].GenMACCSKeys(mol)

        nbits = pyconvert(Int, fp.GetNumBits())

        fp_array = Vector{Int}(undef, nbits - 1) # 166 bits

        for i in 1:nbits-1
            fp_array[i] = pyconvert(Bool, fp.GetBit(i)) ? 1 : 0
        end

        return fp_array
    end

    __init__()

    @testset "test 1" begin 
        # smiles = "C=C(C)C"    
        # smiles = "COCC" 
        # smiles = "C=CN"  

        # smiles = "OS(=O)O"  
        # !!! Bit 81: julia[81] = 1 i rdkit[81] = 0

        smiles = "C1=COC=C1" 
        # !!! Bit 99: julia[99] = 1 i rdkit[99] = 0  
        # !!! Bit 157: julia[157] = 1 i rdkit[157] = 0  

        #smiles = "OS(=O)O"
        # mol = smilestomol(smiles)
        fp_julia = MACCS(false, false)
        julia_fp = fingerprint(smiles, fp_julia)
        # julia_fp = fingerprint(mol, fp_julia)

        rdkit_fp = fingerprint_rdkit(smiles)

        @test length(julia_fp) == 166
        @test length(rdkit_fp) == 166

        # println("Julia MACCS: ", julia_fp)
        # println("RDKit MACCS: ", rdkit_fp)

        # @test julia_fp[58] == rdkit_fp[58]
        # @debug "Bit 58" julia=julia_fp[58] rdkit=rdkit_fp[58]

        # @test julia_fp[166] == rdkit_fp[166]
        # @debug "Bit 166" julia=julia_fp[166] rdkit=rdkit_fp[166]

        # @test julia_fp[1] == rdkit_fp[1]
        # @debug "Bit 1" julia=julia_fp[1] rdkit=rdkit_fp[1]

        @debug("\n=== MACCS bit-by-bit comparison for SMILES: $smiles ===")

        for i in 1:166
            julia_fp[i] == -1 && continue
            #println("Bit $i: Julia=$(julia_fp[i]) RDKit=$(rdkit_fp[i])")

            if julia_fp[i] != rdkit_fp[i]
                @debug "Bit $i" julia=julia_fp[i] rdkit=rdkit_fp[i]
            end
        end
        @debug "=============================================="
    end

    @testset "rule_56 - exact N neighbors (2O + 1C)" begin

        smiles = "ON(O)C" 
        # smiles = "ON(O)CO" 

        # mol = smilestomol(smiles)
        fp_julia = MACCS(false, false)
        julia_fp = fingerprint(smiles, fp_julia)
        rdkit_fp = fingerprint_rdkit(smiles)

        @test length(julia_fp) == 166
        @test length(rdkit_fp) == 166

        @debug("\n=== MACCS bit-by-bit comparison for SMILES rule_56: $smiles ===")

        for i in 1:166
            julia_fp[i] == 9 && continue

            if julia_fp[i] != rdkit_fp[i]
                @debug "Bit $i" julia=julia_fp[i] rdkit=rdkit_fp[i]
            end
        end
        @debug "=============================================="

    end

    @testset "rule_74 - CH3-X-CH3 path" begin

        smiles = "CC(C)" 
        # smiles = "CCC" 

        # mol = smilestomol(smiles)
        fp_julia = MACCS(false, false)
        julia_fp = fingerprint(smiles, fp_julia)
        rdkit_fp = fingerprint_rdkit(smiles)

        @test length(julia_fp) == 166
        @test length(rdkit_fp) == 166

        @debug("\n=== MACCS bit-by-bit comparison for SMILES rule_56: $smiles ===")

        for i in 1:166
            julia_fp[i] == 9 && continue

            if julia_fp[i] != rdkit_fp[i]
                @debug "Bit $i" julia=julia_fp[i] rdkit=rdkit_fp[i]
            end
        end
        @debug "=============================================="

    end

    @testset "rule_79 - N-X-X-N path" begin

        smiles = "NCCN" 
        # smiles = "NCN" 

        fp_julia = MACCS(false, false)
        julia_fp = fingerprint(smiles, fp_julia)
        rdkit_fp = fingerprint_rdkit(smiles)

        @test length(julia_fp) == 166
        @test length(rdkit_fp) == 166

        @debug("\n=== MACCS bit-by-bit comparison for SMILES rule_56: $smiles ===")

        for i in 1:166
            julia_fp[i] == 9 && continue

            if julia_fp[i] != rdkit_fp[i]
                @debug "Bit $i" julia=julia_fp[i] rdkit=rdkit_fp[i]
            end
        end
        @debug "=============================================="

    end

    @testset "rule_80 - N-X-X-X-N path" begin

        smiles = "NCCCN" 
        # smiles = "NCCN" 

        fp_julia = MACCS(false, false)
        julia_fp = fingerprint(smiles, fp_julia)
        rdkit_fp = fingerprint_rdkit(smiles)

        @test length(julia_fp) == 166
        @test length(rdkit_fp) == 166

        @debug("\n=== MACCS bit-by-bit comparison for SMILES rule_56: $smiles ===")

        for i in 1:166
            julia_fp[i] == 9 && continue

            if julia_fp[i] != rdkit_fp[i]
                @debug "Bit $i" julia=julia_fp[i] rdkit=rdkit_fp[i]
            end
        end
        @debug "=============================================="

    end

end


