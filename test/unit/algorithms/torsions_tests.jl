@testset "Topological Torsion Fingerprint Tests" begin
    @testset "testPathFinder" begin
            # this molecule is a ring made from 6 carbon atoms 
            data = "c1ccccc1"
            mol = smilestomol(data)
            # all paths of length 4
            pathList = [[1,2,3,4], [2,3,4,5], [3,4,5,6], [4,5,6,1], [5,6,1,2], [6,1,2,3]]
            paths = MolecularFingerprints.getPathsOfLengthN(mol, 4)
            @test length(pathList) == length(paths)

            # the direction in which the paths are found is arbitrary, to obtain a unique fingerprint the paths are MolecularFingerprints.canonicalized
            # later in TTFP code, this is why we have to check if either the found path or the reversed path is in our list
            for path in paths
                @test path in pathList || reverse(path) in pathList
            end
            # all paths of length 6
            pathList = [[1,2,3,4,5,6],[2,3,4,5,6,1], [3,4,5,6,1,2], [4,5,6,1,2,3], [5,6,1,2,3,4], [6,1,2,3,4,5]]
            paths = MolecularFingerprints.getPathsOfLengthN(mol, 6)
            @test length(pathList) == length(paths)
            for path in paths
                @test path in pathList || reverse(path) in pathList
            end

            # there is only one 6-cycle, but the same cycle will be found 6 times, but we will only use one of them in getTopologicalTorsionFP
            pathList = [[1,2,3,4,5,6,1], [2,3,4,5,6,1,2], [3,4,5,6,1,2,3], [4,5,6,1,2,3,4], [5,6,1,2,3,4,5], [6,1,2,3,4,5,6]]
            paths = MolecularFingerprints.getPathsOfLengthN(mol, 7)
            @test length(pathList) == length(paths)

            for path in paths
                @test path in pathList || reverse(path) in pathList
            end
        end

    @testset "testCanonicalization" begin
        # As described in https://depth-first.com/articles/2021/10/06/molecular-graph-canonicalization/,
        # SMILES representations are not unique. We have to use canonicalization to obtain unique fingerprints.
        # For testing, we calculate the fingerprints from SMILES strings representing the same molecule and check 
        # if the fingerprints are identical.

        data = ["OC(C)C", "CC(C)O", "CC(O)C", "C(C)(C)O", "C(O)(C)C"]
        fingerprints = []
        torsion_calc = TopologicalTorsion()
        for d in data
            mol = MolecularGraph.smilestomol(d)
            fp = fingerprint(mol, torsion_calc)
            push!(fingerprints, fp)
        end
        @test fingerprints[2] == fingerprints[1]
        @test fingerprints[2] == fingerprints[3]
        @test fingerprints[3] == fingerprints[4]
        @test fingerprints[4] == fingerprints[5]

        # this example was taken from https://www.leskoff.com/s01812-0
        rep1 = "O=C(C)Oc1ccccc1C(=O)O"
        rep2 =  "CC(=O)OC1=C(C=CC=C1)C(=O)O"
        mol1 = MolecularGraph.smilestomol(rep1)
        fp1 = fingerprint(mol1, torsion_calc)
        mol2 = MolecularGraph.smilestomol(rep2)
        fp2 = fingerprint(mol2, torsion_calc)

        @test fp1 == fp2

        # this example was taken from https://ctr.fandom.com/wiki/Convert_a_SMILES_string_to_canonical_SMILES
        rep1 = "CN2C(=O)N(C)C(=O)C1=C2N=CN1C"
        rep2 =  "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"
        mol1 = MolecularGraph.smilestomol(rep1)
        fp1 = fingerprint(mol1, torsion_calc)
        mol2 = MolecularGraph.smilestomol(rep2)
        fp2 = fingerprint(mol2, torsion_calc)

        @test fp1 == fp2

        # this example was taken from https://en.wikipedia.org/wiki/Simplified_Molecular_Input_Line_Entry_System 
        # under section "Examples"
        rep1 = "O=Cc1ccc(O)c(OC)c1"
        rep2 = "COc1cc(C=O)ccc1O"
        mol1 = MolecularGraph.smilestomol(rep1)
        fp1 = fingerprint(mol1, torsion_calc)
        mol2 = MolecularGraph.smilestomol(rep2)
        fp2 = fingerprint(mol2, torsion_calc)

        @test fp1 == fp2


        # rings are found multiple times, we only keep the one starting at the lowest numbered vertex
        # test if we keep the correct ring
        paths = [[5,1,3,4,5], [1,3,4,5,1], [3,4,5,1,3], [4,5,1,3,4]]
        keepRing = [false, true, false, false]

        for i = 1:length(paths)
            @test MolecularFingerprints.handleRings(paths[i]) == keepRing[i]
        end

    end
    
    @testset "testAtomCodes" begin
        # check if atom code is the same as in the rdkit implementation
        data = "CCOC"
        mol = MolecularGraph.smilestomol(data)

        piBonds = pi_electron(mol) 
	    atomicNumber = atom_number(mol)
	    deg = degree(mol)

        rdkitAtomCode = [33, 34, 98, 33]
        atomCode = Int[]

        for v in vertices(mol)
            push!(atomCode, MolecularFingerprints.getAtomCode(deg[v], piBonds[v], atomicNumber[v]))
        end

        @assert atomCode == rdkitAtomCode
    end

    @testset "testNumPiBonds" begin
        # check if number of pi bonds is the same as in the rdkit implementation
        rdkitNumPi = [ 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0]
        data = "CN2C(=O)N(C)C(=O)C1=C2N=CN1C"
        mol = MolecularGraph.smilestomol(data)
        numPi = MolecularFingerprints.numPiBonds(mol)
        @assert rdkitNumPi == numPi

        rdkitNumPi = [ 0, 0, 0, 0]
        data = "CCOC"
        mol = MolecularGraph.smilestomol(data)
        numPi = MolecularFingerprints.numPiBonds(mol)
        @assert rdkitNumPi == numPi

        rdkitNumPi = [1, 1, 1, 1, 1, 1]
        data = "c1ccccc1"
        mol = MolecularGraph.smilestomol(data)
        numPi = MolecularFingerprints.numPiBonds(mol)
        @assert rdkitNumPi == numPi
      end

    @testset "testFingerprintResult" begin
        # check if fingerprint is the same as in the rdkit implementation
        torsion_calc = TopologicalTorsion()
        data = "CCOC"
        mol = MolecularGraph.smilestomol(data)
        fp = fingerprint(mol, torsion_calc)
        rdkit_ind = 4320149536
        rdkit_entry = 1
        rdkit_length = 68719476735

        @assert length(fp) == rdkit_length
        @assert findnz(fp)[1][1] == rdkit_ind + 1 # account for zero-based indexing in C++
        @assert findnz(fp)[2][1] == rdkit_entry 

        data = "CCCCC"
        mol = MolecularGraph.smilestomol(data)
        fp = fingerprint(mol, torsion_calc)
        rdkit_ind = 4437590048
        rdkit_entry = 2
        rdkit_length = 68719476735

        @assert length(fp) == rdkit_length
        @assert findnz(fp)[1][1] == rdkit_ind + 1 # account for zero-based indexing in C++
        @assert findnz(fp)[2][1] == rdkit_entry 

        data = "c1ccccc1"
        mol = MolecularGraph.smilestomol(data)
        fp = fingerprint(mol, torsion_calc)
        rdkit_ind = 5513433129
        rdkit_entry = 6
        rdkit_length = 68719476735

        @assert length(fp) == rdkit_length
        @assert findnz(fp)[1][1] == rdkit_ind + 1 # account for zero-based indexing in C++
        @assert findnz(fp)[2][1] == rdkit_entry 

        data = "CN2C(=O)N(C)C(=O)C1=C2N=CN1C"
        mol = MolecularGraph.smilestomol(data)
        fp = fingerprint(mol, torsion_calc)
        rdkit_ind = [5647929385, 5647929888, 5647929897, 5656302122, 9808417312, 9808417322, 9808663082, 9808679456, 9808679466, 9942634538, 9942880810, 9942880841, 9942880842, 9942896681, 9942897184, 9942897194, 13969412650, 13969412682, 13969429024, 13969429034]
        rdkit_entry = [1, 4, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2, 2, 1, 1, 3, 3]
        rdkit_length = 68719476735

        @assert length(fp) == rdkit_length
        @assert findnz(fp)[1] == rdkit_ind .+ 1 # account for zero-based indexing in C++
        @assert findnz(fp)[2] == rdkit_entry 
    end
    
end
