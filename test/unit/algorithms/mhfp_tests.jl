using Test
using MolecularGraph
using MolecularFingerprints
import MolecularFingerprints: fingerprint

# must be included, otherwise mhfp.jl cannot be included
include("../../../src/interface.jl")

include("../../../src/algorithms/mhfp.jl")

@testset "MHFP Fingerprint Tests" begin
    @testset "MHFP Molecular shingling tests" begin
        ##### Note on how the testing is done #############################################
        # TODO move this text into the documentation
        
        # In general, our implmentation cannot be expected to yield identical shinglings
        # as the fct written by the original authors in python, due to several reasons:
        
        # When writing smiles strings using MolecularGraph.jl, we can only pass a molecule
        # object to the corresponding function. To generate the smiles string of a
        # substructure of a molecule, we thus have to create a new molecule containing only
        # these atoms, which we do using the function "induced_subgraph" from Graphs.jl
        # This however leads to the issue that information about bonds between the atoms in
        # the substructure and atoms outside the substructure gets lost. For instance, this
        # can mean that information about the aromaticity of atoms is lost, which in that 
        # case decides whether a, say, carbon atom is written "c" (aromatic carbon) or "C"
        # (non-aromatic) in the smiles string of the substructure. In another case, the 
        # python rdkit function returned "[nH]" for an atom, which denotes an aromatic N
        # atom with an additional hydrogen; while our implementation simply returns "n".
        # This as well could be due information lost due to our approach.
        # The rdkit implementation used by the original authors seems to allow to write 
        # smiels strings of substructures without losing this kind of information. Thus,
        # our shinglings don't necessary match.

        # Furthermore, there can exist different smiles strings that are both a valid 
        # description of the same molecule or substructure. It seems like the functions in
        # MolecularGraph.jl do not necessary return the same variant as rdkit in python.
        # (although it is hard to rule out that this may also be due to how we create
        # the substructure as an induced subgraph)

        
        # In conclusion, it is not always possible to test our shingling generating 
        # functions by comparing their results to the results from the original authors
        # function.
        # Instead, what we do to test is the following:
        # - We verify for a larger molecule (which was used by the original authors
        #   to test their code) that the number of strings found by our implementation
        #   matches their implementation.
        #   - Note that we cannot compare the actual strings themselves, as they differ
        #   - Note that we cannot verify the number of unique strings, only the number
        #       of non-unique strings found. This is because our strings differ from
        #       those of the original authors, so we may also have a different number
        #       of strings that overlap, and thus a different number of unique strings
        # - We verify for a small, simple molecule that the actual strings also match
        #   - This is only possible for the strings generated from rings and individual
        #       atoms of the molecule. The strings from the circular substructures of
        #       the molecule differ already for this smaller molecule.
        # - Lastly, we verify for the large atom that the shingling returned by 
        #   mhfp_shingling_from_mol is the union of the strings returned by
        #   smiles_from_rings, smiles_from_atoms and smiles_from_circular_substructures.
        
        # A last note: As the function implemented by the original authors returns the
        # entire shingling at once, the reference strings that we test against were 
        # generated manually by running parts of their function mhfp_shingling_from_mol
        # separately.

        ##### Definition of the molecules used for testing ################################
        
        # Larger molecule, taken from mhfp python implementation, 
        # written by the original authors
        # https://github.com/reymond-group/mhfp/blob/master/test/test_encoder.py
        mol = MolecularGraph.smilestomol(
            "CCCC1=NN(C2=C1NC(=NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C")

        # Simpler molecule (Benzoic acid), taken from 
        # https://chemicbook.com/2021/02/13/smiles-strings-explained-for-beginners-part-1.html
        simpler_mol = smilestomol("c1cc(C(O)=O)ccc1")

        @testset "Shingling snippet from rings test" begin
            ##### Testing if number of strings is correct, on the larger molecule #########

            # this was manually calculated using the original author's code applied to the
            # above molecule
            ref_shingling_snippet_rings = ["c1cnnc1", "c1cnc[nH]c1", "c1ccccc1", "C1CNCCN1"]
            calculated_shingling_snippet_rings = smiles_from_rings(mol)
            
            # test if number of returned strings is correct
            @test length(ref_shingling_snippet_rings) == length(
                calculated_shingling_snippet_rings)
            
            # Note: In fact, there are even molecules for which even the number of rings 
            # found in our implementation is not the same as the number of rings found 
            # in the original authors implementation. The reason for this is that the 
            # original implementation uses rdkits "symmetrisized sssr", which sometimes
            # adds an extra ring to the sssr to make it symmetric in some sense. This
            # function is not available in MolecularGraph (nor in RDKitMinimalLib), which
            # is why we can only use the standard sssr.
            # For example for the molecule cubane (smilestomol("C12C3C4C1C5C2C3C45"))
            # has 5 rings in the stanrard sssr, while the symmetrisized sssr contains 6.


            ##### Testing if the actual string generated matches, on the simpler molecule #
            ref_shingling_snippet_rings_simplemol = ["c1ccccc1"]
            @test ref_shingling_snippet_rings_simplemol == smiles_from_rings(simpler_mol)
            
        end

        @testset "Shingling snippet from atoms tests" begin
            ##### Testing if the number of strings is correct, on the larger molecule #####
            
            # this was manually calculated using the original author's code applied to the
            # above molecule
            ref_shingling_snippet_atoms = ["C", "C", "C", "c", "n", "n", "c", "c", "[nH]", 
            "c", "n", "c", "O", "c", "c", "c", "c", "c", "c", "S", "O", "O", "N", "C", "C",
            "N", "C", "C", "C", "O", "C","C", "C"]

            @test length(ref_shingling_snippet_atoms) == length(smiles_from_atoms(mol))

            ##### Testing if the actual strings are correct, on the smaller molecule ######

            ref_shingling_snippet_atoms_simplemol = [
                "c", "c", "c", "C", "O", "O", "c", "c", "c"]

            @test ref_shingling_snippet_atoms_simplemol == smiles_from_atoms(simpler_mol)
        end

        @testset "Shingling snippet from circular substructures tests" begin
            ##### Testing if the number of strings is correct, on the larger molecule #####
            
            ref_shingling_snippet_circular_substructures = [
                "CC",
                "CCC",
                "CCCc",
                "C(C)C",
                "C(C)Cc",
                "C(C)Cc(c)n",
                "C(C)c",
                "C(CC)c(c)n",
                "C(CC)c(nn)c(c)[nH]",
                "c(c)(C)n",
                "c(CC)(nn)c(c)[nH]",
                "c1(CCC)nn(C)c(c)c1[nH]c",
                "n(c)n",
                "n(c(c)C)n(c)C",
                "n1c(CC)c([nH])c(c)n1C",
                "n(c)(C)n",
                "n(C)(nc)c(c)c",
                "n1(C)nc(C)c([nH])c1c(n)=O",
                "c(c)(c)n",
                "c(c(n)=O)(c(c)[nH])n(C)n",
                "c12c(=O)nc[nH]c1c(C)nn2C",
                "c(c)(c)[nH]",
                "c([nH]c)(c(C)n)c(c)n",
                "c12[nH]c(-c)nc(=O)c1n(C)nc2CC",
                "[nH](c)c",
                "[nH](c(-c)n)c(c)c",
                "[nH]1c(-c(c)c)ncc(n)c1c(C)n",
                "c(-c)([nH])n",
                "c(nc)([nH]c)-c(c)c",
                "c1(-c(cc)c(c)O)nc(=O)cc(c)[nH]1",
                "n(c)c",
                "n(c(-c)[nH])c(c)=O",
                "n1c(-c(c)c)[nH]cc(n)c1=O",
                "c(c)(n)=O",
                "c(=O)(nc)c(c)n",
                "c1(=O)nc(-c)[nH]c(c)c1n(C)n",
                "O=c",
                "O=c(c)n",
                "O=c(nc)c(c)n",
                "c(c)(c)-c",
                "c(cc)(-c([nH])n)c(c)O",
                "c1(-c(nc)[nH]c)cc(S)ccc1OC",
                "c(c)(c)O",
                "c(cc)(OC)c(c)-c",
                "c1(OCC)ccccc1-c([nH])n",
                "c(c)c",
                "c(cc)c(c)O",
                "c1cc(S)cc(-c)c1OC",
                "c(c)c",
                "c(cc)c(c)S",
                "c1cc(O)ccc1S(N)(=O)=O",
                "c(c)(c)S",
                "c(cc)(cc)S(N)(=O)=O",
                "c1(S(=O)(=O)N(C)C)cccc(-c)c1",
                "c(c)c",
                "c(c(c)-c)c(c)S",
                "c1c(S(N)(=O)=O)ccc(O)c1-c([nH])n",
                "S(c)(N)(=O)=O",
                "S(=O)(=O)(c(c)c)N(C)C",
                "S(=O)(=O)(c(cc)cc)N(CC)CC",
                "O=S",
                "O=S(c)(N)=O",
                "O=S(=O)(c(c)c)N(C)C",
                "O=S",
                "O=S(c)(N)=O",
                "O=S(=O)(c(c)c)N(C)C",
                "N(C)(C)S",
                "N(CC)(CC)S(c)(=O)=O",
                "N1(S(=O)(=O)c(c)c)CCNCC1",
                "C(C)N",
                "C(CN)N(C)S",
                "C1CN(C)CCN1S(c)(=O)=O",
                "C(C)N",
                "C(CN)N(C)C",
                "C1CN(S)CCN1C",
                "N(C)(C)C",
                "N(C)(CC)CC",
                "N1(C)CCNCC1",
                "C(C)N",
                "C(CN)N(C)C",
                "C1CN(S)CCN1C",
                "C(C)N",
                "C(CN)N(C)S",
                "C1CN(C)CCN1S(c)(=O)=O",
                "CN",
                "CN(C)C",
                "CN(CC)CC",
                "O(c)C",
                "O(CC)c(c)c",
                "O(CC)c(cc)c(c)-c",
                "C(C)O",
                "C(C)Oc",
                "C(C)Oc(c)c",
                "CC",
                "CCO",
                "CCOc",
                "Cn",
                "Cn(c)n",
                "Cn(nc)c(c)c"]

            @test length(ref_shingling_snippet_circular_substructures) == length(
                smiles_from_circular_substructures(mol, 3, 1))
        end

        @testset "Complete shingling tests" begin
            ##### Testing mhfp_shingling_from_mol #########################################
            calculated_shingling = mhfp_shingling_from_mol(mol, 3, true, 0)

            # Test that the shingling returned from mhfp_shingling_from_mol is the union
            # of smiles_from_rings, smiles_from_atoms & smiles_from_circular_substructures
            @test symdiff(calculated_shingling, union(smiles_from_rings(mol), 
                smiles_from_atoms(mol), 
                smiles_from_circular_substructures(mol, 3, 1))) == []
        end
    end
    
end
