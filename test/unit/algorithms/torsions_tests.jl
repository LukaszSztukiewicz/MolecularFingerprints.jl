using Test
using MolecularGraph
using MolecularFingerprints

@testset "Topological Torsion Fingerprint Tests" begin
    @testset "testPathFinder" begin
            data = "c1ccccc1"
            mol = smilestomol(data)
            pathList = [[1,2,3,4], [2,3,4,5], [3,4,5,6], [4,5,6,1], [5,6,1,2], [6,1,2,3]]
            
            paths = MolecularFingerprints.getPathsOfLengthN(mol, 4)
            @info paths
            @test length(pathList) == length(paths)
            #@test pathList == paths
        end

    
    
        # torsion_calc = TopologicalTorsion(4)  # Example type TopologicalTorsion with path length 4
        # mol = MolecularGraph.smilestomol(data)
        # fp = fingerprint(mol, torsion_calc)

end
