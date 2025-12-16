using Test
using MolecularGraph
using MolecularFingerprints

@testset "ECFP Fingerprint Tests" begin

    # dummy test to ensure the test file is being executed
    # Replace with actual ECFP fingerprint tests when implemented

    mol = MolecularGraph.smilestomol("CC(=O)O")
    ecfp_calc = ECFP{2}(2)  # Example type ECFP{2} with radius 2
    fp = fingerprint(mol, ecfp_calc)
    @test length(fp) == 2048 #FIXME probably should be 1024, now 2048 for the demonstration purpose, change if needed

end
