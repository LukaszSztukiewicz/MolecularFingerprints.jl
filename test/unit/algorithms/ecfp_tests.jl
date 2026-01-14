using Test
using MolecularGraph
using MolecularFingerprints

bitset_to_string(bitset) = join(Int.(bitset), "")

@testset "ECFP Algorithm Tests" begin

    @testset "Constructor" begin
        # Test valid construction
        ecfp2 = ECFP{1024}(2)
        @test ecfp2.radius == 2

        ecfp3 = ECFP{2048}(3)
        @test ecfp3.radius == 3

        # Test with radius 0
        ecfp0 = ECFP{512}(0)
        @test ecfp0.radius == 0

         # Test invalid radius
        @test_throws ArgumentError ECFP{1024}(-1)

        # Test invalid size
        @test_throws ArgumentError ECFP{0}(2)
        @test_throws ArgumentError ECFP{-1}(2)
    end

    @testset "Basic Properties" begin
        mol = smilestomol("CCO")  # Ethanol
        ecfp = ECFP{1024}(2)

        fp = fingerprint(mol, ecfp)

        # Check output type and size
        @test fp isa BitVector
        @test length(fp) == 1024

        # Check that some bits are set
        @test any(fp)

        # Check that not all bits are set
        @test !all(fp)
    end

     @testset "Different Sizes" begin
        mol = smilestomol("CCO")

        # Test different fingerprint sizes
        for size in [512, 1024, 2048, 4096]
            ecfp = ECFP{size}(2)
            fp = fingerprint(mol, ecfp)
            @test length(fp) == size
            @test fp isa BitVector
        end
    end

    @testset "Different Radii" begin
        mol = smilestomol("CC(=O)NCCC1=CNc2c1cc(OC)cc2")

        # Test radius 2 (ECFP4, diameter 4)
        ecfp2 = ECFP{2048}(2)
        fp2 = fingerprint(mol, ecfp2)
        @test length(fp2) == 2048
        @test any(fp2)

        # Test radius 3 (ECFP6, diameter 6)
        ecfp3 = ECFP{2048}(3)
        fp3 = fingerprint(mol, ecfp3)
        @test length(fp3) == 2048
        @test any(fp3)

        # Different radii should generally produce different fingerprints
        # (though this is not guaranteed for all molecules)
        @test fp2 != fp3
    end

    @testset "Determinism" begin
        mol = smilestomol("CCO")
        ecfp = ECFP{2048}(2)

        # Generate fingerprint multiple times
        fp1 = fingerprint(mol, ecfp)
        fp2 = fingerprint(mol, ecfp)
        fp3 = fingerprint(mol, ecfp)

        # All should be identical
        @test fp1 == fp2
        @test fp2 == fp3
    end

    @testset "Different Molecules" begin
        ecfp = ECFP{2048}(2)

        # Test different molecules
        mol_ethanol = smilestomol("CCO")
        mol_methanol = smilestomol("CO")
        mol_propanol = smilestomol("CCCO")
        mol_benzene = smilestomol("c1ccccc1")

        fp_ethanol = fingerprint(mol_ethanol, ecfp)
        fp_methanol = fingerprint(mol_methanol, ecfp)
        fp_propanol = fingerprint(mol_propanol, ecfp)
        fp_benzene = fingerprint(mol_benzene, ecfp)

        # Different molecules should have different fingerprints
        @test fp_ethanol != fp_methanol
        @test fp_ethanol != fp_benzene
        @test fp_methanol != fp_benzene

        # Methanol should be a substructure-like of ethanol
        # (some overlap expected but not identical)
        @test fp_ethanol != fp_methanol
    end

    @testset "Similarity" begin
        ecfp = ECFP{2048}(2)

        # Similar molecules should have similar fingerprints
        mol1 = smilestomol("CCO")  # Ethanol
        mol2 = smilestomol("CCCO")  # Propanol
        mol3 = smilestomol("c1ccccc1")  # Benzene

        fp1 = fingerprint(mol1, ecfp)
        fp2 = fingerprint(mol2, ecfp)
        fp3 = fingerprint(mol3, ecfp)

        # Ethanol and propanol should have more overlap than ethanol and benzene
        @test tanimoto(fp1, fp2) > tanimoto(fp1, fp3)
    end

    @testset "Ring Structures" begin
        ecfp = ECFP{2048}(2)

        # Test molecules with rings
        benzene = smilestomol("c1ccccc1")
        cyclohexane = smilestomol("C1CCCCC1")
        toluene = smilestomol("Cc1ccccc1")

        fp_benzene = fingerprint(benzene, ecfp)
        fp_cyclohexane = fingerprint(cyclohexane, ecfp)
        fp_toluene = fingerprint(toluene, ecfp)

        @test any(fp_benzene)
        @test any(fp_cyclohexane)
        @test any(fp_toluene)

        # Benzene and toluene should be at least as similar as benzene and cyclohexane
        # (toluene contains benzene as a substructure)
        @test tanimoto(fp_benzene, fp_toluene) >= tanimoto(fp_benzene, fp_cyclohexane)
    end

end

# @testset "ECFP Fingerprint Tests" begin

#     # dummy test to ensure the test file is being executed
#     # Replace with actual ECFP fingerprint tests when implemented

#     mol = MolecularGraph.smilestomol("CC(=O)O")
#     ecfp_calc = ECFP{2}(2)  # Example type ECFP{2} with radius 2
#     fp = fingerprint(mol, ecfp_calc)
#     @test length(fp) == 2048 #FIXME probably should be 1024, now 2048 for the demonstration purpose, change if needed

# end
