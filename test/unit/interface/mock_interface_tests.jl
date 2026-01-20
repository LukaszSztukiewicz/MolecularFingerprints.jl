@testset "MolecularFingerprints Interface" begin

    # --- Mocks ---
    struct MockFingerprintCalc <: AbstractFingerprint end
    struct MockDescriptorCalc <: AbstractDescriptor end

    function Main.fingerprint(mol::MolGraph, calc::MockFingerprintCalc)
        BitVector([nv(mol) % 2 == 0]) 
    end
    
    function Main.fingerprint(mol::MolGraph, ::MockDescriptorCalc)
        Float64[nv(mol) * 1.5]
    end

    fp_calc = MockFingerprintCalc()
    desc_calc = MockDescriptorCalc()

    @testset "Interface Dispatch" begin
        # 'CCO' has 3 non-H atoms. 3%2 is not 0, so should be [0]
        @test fingerprint("CCO", fp_calc) == BitVector([0])
        @test fingerprint("CCO", desc_calc) == [4.5]
    end

    @testset "Batch Processing & Parallelization" begin
        smiles_list = ["CCO", "c1ccccc1", "CC(=O)NCCC1=CNc2c1cc(OC)cc2"]
        fp_results = [fingerprint(smiles, fp_calc) for smiles in smiles_list]
        desc_results = [fingerprint(smiles, desc_calc) for smiles in smiles_list]
        @test length(fp_results) == 3
        @test length(desc_results) == 3
        @test all(fp isa BitVector for fp in fp_results)
        @test all(desc isa Vector{Float64} for desc in desc_results)        
    end

    @testset "Error Handling" begin
        @test_throws BoundsError fingerprint(String[], fp_calc)
    end
end