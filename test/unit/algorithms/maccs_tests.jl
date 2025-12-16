using Test
using MolecularGraph
using MolecularFingerprints

@testset "MACCS Fingerprint - acetic acid CC(=O)O" begin

    @testset "Bit-based MACCS (count = false)" begin
        smiles = "CC(=O)O"
        fp = MACCSFingerprint(false, false)
        println("Created fingerprint: MACCSFingerprint(false, false) for smiles: ", smiles)

        mol = MolecularGraph.smilestomol(smiles)
        println("Molecule from SMILES: ", mol)
        
        v = fingerprint(mol, fp)
        println("  Fingerprint vector:")
        println("  ", v)

        bits = findall(!=(0), v)
        println("  Values at these positions: ", v[bits])

        @test length(v) == 166
        @test v[2] == 1   # two oxygens - true = 1 
        @test v[3] == 1   # C=O - true = 1 
        @test sum(v) == 2
    end

    @testset "Count-based MACCS (count = true)" begin
        smiles = "CC(=O)O"
        fp = MACCSFingerprint(true, false)
        println("Created fingerprint: MACCSFingerprint(true, false) for smiles: ", smiles)
        
        mol = MolecularGraph.smilestomol(smiles)
        println("Molecule from SMILES: ", mol)
        
        v = fingerprint(mol, fp)
        println("  Fingerprint vector:")
        println("  ", v)

        bits = findall(!=(0), v)
        println("  Values at these positions: ", v[bits])

        @test length(v) == 166
        @test v[2] == 2   # two oxygens - 2
        @test v[3] == 1   # C=O - 1
        @test sum(v) == 3
        
    end

end
