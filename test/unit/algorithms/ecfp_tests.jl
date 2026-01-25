const rdkitChem = pyimport("rdkit.Chem")
const rdkitMolDescriptors = pyimport("rdkit.Chem.rdMolDescriptors")
const skfpFingerprints = pyimport("skfp.fingerprints")

function rdkit_get_atom_invariants(smiles::AbstractString)
    mol::Py = rdkitChem.MolFromSmiles(smiles)
    result = rdkitMolDescriptors.GetConnectivityInvariants(mol; includeRingMembership = true)
    return pyconvert(Vector{UInt32}, result)
end

function rdkit_fingerprint(smiles::AbstractString, fp_size=64, radius=2)
    fp = skfpFingerprints.ECFPFingerprint(fp_size, radius)
    return pyconvert(Vector{Bool}, fp.transform([smiles])[0])
end

bitset_to_string(bitset) = join(Int.(bitset), "")



@testset "ECFP Algorithm Tests" begin
    @testset "Atom invariant hashes: Comparison to RDKit" begin
        test_cases = [
            "C",
            "CC",
            "CCC",
            "CC[2H]",
            "CC(=O)OCC[N+](C)(C)C",
            "CC(C[N+](C)(C)C)OC(=O)C",
            "O=C1CCCN1CC#CC[N+](C)(C)C",
            "NC(=O)OCC[N+](C)(C)C",
            "CC(C[N+](C)(C)C)OC(=O)N",
            "COC(=O)C1=CCCN(C1)C",
            "O=C1CCCN1CC#CCN1CCCC1",
            "CON=CC1=CCCN(C1)C",
            "CN1CC(=CCC1)C(=O)OCC#C",
        ]

        for tc in test_cases
            actual_hashes = MolecularFingerprints.get_atom_invariants(tc)
            expected_hashes = rdkit_get_atom_invariants(tc)

            @testset "Molecule $tc, atom index $i" for (i, (actual, expected)) in enumerate(zip(actual_hashes, expected_hashes))
                # Check hash result
                @test actual == expected
            end
        end
    end

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
        @test_throws DomainError ECFP{1024}(-1)

        # Test invalid size
        @test_throws DomainError ECFP{0}(2)
        @test_throws DomainError ECFP{-1}(2)
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
        @test tanimoto_similarity(fp1, fp2) > tanimoto_similarity(fp1, fp3)
    end

    @testset "Fingerprint: Comparison to RDKit" begin
        test_cases = [
            ("C", 2),
            ("CC", 2),
            ("CCC", 2),
            ("C", 3),
            ("CC", 3),
            ("CCC", 3),
            ("C", 4),
            ("CC", 4),
            ("CCC", 4),
        ]

        @testset "Molecule $smiles, radius=$radius" for (smiles, radius) in test_cases
            expected = rdkit_fingerprint(smiles, 16, radius)

            mol = smilestomol(smiles)
            fp = ECFP{16}(radius)
            actual = fingerprint(mol, fp)

            @test actual == expected
        end
    end

    @testset "Unsupported Bond Type" begin
        mol = smilestomol("C1C:C:C:C:C1")
        fp = ECFP{16}(3)
        @test_throws "Unsupported bond type" fingerprint(mol, fp)
    end
end
