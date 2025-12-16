using Test
using MolecularGraph
using MolecularFingerprints

@testset "Topological Torsion Fingerprint Tests" begin

    # dummy test to ensure the test file is being executed
    # Replace with actual Topological Torsion fingerprint tests when implemented

    mol = MolecularGraph.smilestomol("CC(=O)O")
    torsion_calc = TopologicalTorsion{2}(2)  # Example type TopologicalTorsion{2} with radius 2
    fp = fingerprint(mol, torsion_calc)
    @test length(fp) == 1024

end
