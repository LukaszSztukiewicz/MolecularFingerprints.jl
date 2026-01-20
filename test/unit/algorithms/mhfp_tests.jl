import MolecularFingerprints: fingerprint

# must be included, otherwise mhfp.jl cannot be included
include("../../../src/interface.jl")

include("../../../src/algorithms/mhfp.jl")

@testset "MHFP Fingerprint Tests" begin

    ##### Definition of the molecules used for testing ################################
        
        # Larger molecule, taken from mhfp python implementation, 
        # written by the original authors
        # https://github.com/reymond-group/mhfp/blob/master/test/test_encoder.py
        mol = MolecularGraph.smilestomol(
            "CCCC1=NN(C2=C1NC(=NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C")

        # Simpler molecule (Benzoic acid), taken from 
        # https://chemicbook.com/2021/02/13/smiles-strings-explained-for-beginners-part-1.html
        simpler_mol = smilestomol("c1cc(C(O)=O)ccc1")

    
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
            
            # this was manually calculated using the original author's code applied to the
            # above molecule
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

            ## calculating with radius 3, min_radius 1 
            calc_shingling_snippet_circ_subst_r_3_mr_1 = smiles_from_circular_substructures(
                mol, 3, 1)
            # ensure length of shingling is the same as the reference shingling
            @test length(ref_shingling_snippet_circular_substructures) == length(
                calc_shingling_snippet_circ_subst_r_3_mr_1
                )

            ## calculating with radius 2, min_radius 1
            calc_shingling_snippet_circ_subst_r_2_mr_1 = smiles_from_circular_substructures(
                mol, 2, 1)
            # ensure shingling is smaller as a smaller radius was used
            @test length(ref_shingling_snippet_circular_substructures) > length(
                calc_shingling_snippet_circ_subst_r_2_mr_1
                )
            # ensure the smaller shingling is contained in the larger shingling
            @test issubset(calc_shingling_snippet_circ_subst_r_2_mr_1, 
                calc_shingling_snippet_circ_subst_r_3_mr_1)

            ## calculating with radius 3, min_radius 2
            calc_shingling_snippet_circ_subst_r_3_mr_2 = smiles_from_circular_substructures(
                mol, 3, 2)
            # ensure shingling is smaller as a smaller radius was used
            @test length(ref_shingling_snippet_circular_substructures) > length(
                calc_shingling_snippet_circ_subst_r_3_mr_2
                )
            # ensure the smaller shingling is contained in the larger shingling
            @test issubset(calc_shingling_snippet_circ_subst_r_3_mr_2, 
                calc_shingling_snippet_circ_subst_r_3_mr_1)

        end

        @testset "Complete shingling tests" begin
            ##### Testing mhfp_shingling_from_mol #########################################
            
            # set up calculator with parameters
            calculator = MHFP(3, 0, true)  # radius, min_radius, rings

            # calculate shingling
            calculated_shingling = mhfp_shingling_from_mol(mol, calculator)

            # Test that the shingling returned from mhfp_shingling_from_mol is the union
            # of smiles_from_rings, smiles_from_atoms & smiles_from_circular_substructures
            @test symdiff(calculated_shingling, union(smiles_from_rings(mol), 
                smiles_from_atoms(mol), 
                # Note: min_radius is now 1, even though we set it to 0 in the calculator 
                # above. This is because the mhfp_shingling_from_mol function increases the
                # min_radius to at least 1 before calling
                # smiles_from_circular_substructures, as the special case of radius 0 is 
                # already taken care of by the function smiles_from_atoms
                smiles_from_circular_substructures(mol, 3, 1))  # radius, min_radius
                ) == []
            
            # test that explicitly given hydrogens are removed and not included in smiles
            simple_mol_with_hydrogen = smilestomol("c1cc(C(O)=O)ccc1")
            add_hydrogens!(simple_mol_with_hydrogen)
            # the molecule now has explicit H in its graph

            # initialize calculator, with min_radius 0, such that if hydrogen removal 
            # doesn't work correctly, "[H]" will be part of the shingling
            calc_with_atoms = MHFP(3, 0)  # radius, min_radius

            # test that "[H]" is not in the shingling
            @test "[H]" ∉ mhfp_shingling_from_mol(simple_mol_with_hydrogen, MHFP(3, 0))

        end
    end

    @testset "MHFP Hashing function tests" begin
        ##### Testing mhfp_hash_from_molecular_shingling ##################################

        for n_permutations in [512, 2048]  # test default 2048, and 2048/4 = 512
            for seed in [42, (1 << 31)]  # test default value 42 and some high number

                # the values radius, min_radius, rins can be left at default because they
                # are not used in hashing process
                mhfp_calc = MHFP(n_permutations=n_permutations, seed=seed)

                ## Testing length of generated shingling ##################################

                # We can only test the length of the hashed shingling, as the actual numbers
                # in our implementation differ from those gotten from the original authors 
                # code. However, the rng used by numpy seems to differ from the one used in 
                # julia, hence, this is not surprising.

                # Example shingling taken from the original authors tests.
                # (See https://github.com/reymond-group/mhfp/blob/ea514f8fd4b21b0d0d732452cf7062c282edfbde/test/test_encoder.py#L12)
                ref_shingling = sort(["c1cnnc1", "Cn(nc)c(c)c", "C(C)Oc", "c1(OCC)ccccc1-c([nH])n", "S(c)(N)(=O)=O", "c(nc)([nH]c)-c(c)c", "CN(C)C", "n(c(-c)[nH])c(c)=O", "C(C)O", "N(C)(C)S", "S(=O)(=O)(c(cc)cc)N(CC)CC", "CC", "c1(CCC)nn(C)c(c)c1[nH]c", "c1(S(=O)(=O)N(C)C)cccc(-c)c1", "N(C)(C)C", "c(cc)c(c)S", "N1(S(=O)(=O)c(c)c)CCNCC1", "c(c)(c)[nH]", "O=c(nc)c(c)n", "N(C)(CC)CC", "n(c(c)C)n(c)C", "C1CNCCN1", "c1ccccc1", "C(C)N", "n(c)(C)n", "C(CC)c(c)n", "c(c(c)-c)c(c)S", "O(c)C", "CCO", "CN(CC)CC", "[nH](c)c", "n(c)c", "n1c(-c(c)c)[nH]cc(n)c1=O", "N(CC)(CC)S(c)(=O)=O", "C(CN)N(C)C", "S(=O)(=O)(c(c)c)N(C)C", "O(CC)c(cc)c(c)-c", "C(C)Oc(c)c", "C(C)Cc", "C(CN)N(C)S", "c1(-c(nc)[nH]c)cc(S)ccc1OC", "c(cc)(OC)c(c)-c", "n1c(CC)c([nH])c(c)n1C", "c(c)(c)-c", "c([nH]c)(c(C)n)c(c)n", "c(c)(c)O", "c1cc(O)ccc1S(N)(=O)=O", "O=S(c)(N)=O", "c12[nH]c(-c)nc(=O)c1n(C)nc2CC", "c(-c)([nH])n", "c1(=O)nc(-c)[nH]c(c)c1n(C)n", "c1cc(S)cc(-c)c1OC", "O(CC)c(c)c", "c(c)(c)n", "C(C)C", "n(c)n", "CCOc", "c(cc)(-c([nH])n)c(c)O", "c(CC)(nn)c(c)[nH]", "c(c)(c)S", "CCC", "N1(C)CCNCC1", "c(c(n)=O)(c(c)[nH])n(C)n", "n(C)(nc)c(c)c", "c(c)(n)=O", "Cn(c)n", "c(cc)(cc)S(N)(=O)=O", "O=c(c)n", "c(c)c", "n1(C)nc(C)c([nH])c1c(n)=O", "c1c(S(N)(=O)=O)ccc(O)c1-c([nH])n", "C(C)Cc(c)n", "C(C)c", "c(=O)(nc)c(c)n","[nH](c(-c)n)c(c)c", "C(CC)c(nn)c(c)[nH]", "O=c", "c1(-c(cc)c(c)O)nc(=O)cc(c)[nH]1", "c(cc)c(c)O", "c(c)(C)n", "O=S", "c1cnc[nH]c1", "O=S(=O)(c(c)c)N(C)C", "C1CN(C)CCN1S(c)(=O)=O", "c12c(=O)nc[nH]c1c(C)nn2C", "C1CN(S)CCN1C", "[nH]1c(-c(c)c)ncc(n)c1c(C)n", "Cn", "CCCc", "CN"])
                
                @test n_permutations == length(mhfp_hash_from_molecular_shingling(
                    ref_shingling, 
                    mhfp_calc))

                
                ## Testing similarity of hashed vectors depending on the similarity of the 
                ## original vectors #######################################################
                
                # The original authors claim that the MinHash can be used to estimate the 
                # tanimoto similarity of two sets.
                # However, testing the claim with their own implementation yields that the 
                # tanimoto similarity of the hashed vectors indeed relates to the similarity
                # of the sets that were hashed, but it's not a 1:1 relation.
                # Instead, for sets with low similarity, the corresponding hashed vectors 
                # have a similarity which is only around half as large.
                # For sets with higher similarity, the similarity of the corresponding 
                # hashed vectors is still smaller than the similarity of the original sets, 
                # but the difference is less pronounced.
                
                # To verify that our hashing implementation behaves similarly, we generate 
                # pairs of sets that share a certain amount of entries, respectively, which 
                # implies a certain tanimoto similarity of the pair of sets.
                # Then their MinHash vectors are generated, and then the tanimoto similarity
                # is calculated for these vectors.
                # This is repeated 1000 times and the average calculated, to avoid random 
                # effects.

                # This process was performed in python with the original authors MinHash
                # implementation (specifically, their function from_molecular_shingling), to 
                # generate the reference values we test against.
                
                # Below, the same process is performed with our implementation and the 
                # results compared to those test values.

                # We will generate 2 sets that share 2, 12, 30, 50 and 80 elements, 
                # respectively
                overlap_radii = [1, 6, 15, 25, 40]  # = [2, 12, 30, 50, 80] / 2
                
                # reference values calculated with the original python implementation
                ref_minhash_tanimoto_values = [0.010, 0.064, 0.177, 0.333, 0.667]
                
                for (i, overlap_radius) in enumerate(overlap_radii)
                    original_tanimoto_values = []
                    minhash_tanimoto_values = []
                    for j in 1:500  # repeat 500 times to avoid random factors
                        # Random.seed!(j)
                        
                        # random set of strings of size 100
                        test_set = [randstring(25) for k in 1:100]
                        
                        # split the 100 strings into two sets, sharing 2*overlap_radius 
                        # elements
                        test_set_1 = test_set[begin:50 + overlap_radius]
                        test_set_2 = test_set[50 + 1 - overlap_radius:end]

                        # As probably almost all non-overlapping strings are pairwise 
                        # different, we should get a tanimoto (jaccard) similarity of 
                        # approximately:
                        # 2 / 100 = 0.02,
                        # 12 / 100 = 0.12,
                        # 30 / 100 = 0.3,
                        # 50 / 100 = 0.5 and
                        # 80 / 100 = 0.8, respectively.
                        original_tanimoto_similarity = length(
                            intersect(test_set_1, test_set_2)) / length(
                                union(test_set_1, test_set_2))

                        push!(original_tanimoto_values, original_tanimoto_similarity)
                        
                        # Calculate hash vectors from the given test sets
                        minhash_1 = mhfp_hash_from_molecular_shingling(
                            test_set_1, mhfp_calc)
                        minhash_2 = mhfp_hash_from_molecular_shingling(
                            test_set_2, mhfp_calc)

                        # Calculating tanimoto (jaccard) similarity of the hash vectors.
                        tanimoto_of_minhash =  length(
                            intersect(minhash_1, minhash_2)) / length(
                                union(minhash_1, minhash_2))
                        
                            push!(minhash_tanimoto_values, tanimoto_of_minhash)
            
                    end

                    # Calculate average tanimoto similarity for original sets
                    avg_original_tanimoto_value = sum(original_tanimoto_values) / length(
                        original_tanimoto_values)
                    # Calculate average tanimoto similarity for the hashed vectors
                    avg_minhash_tanimoto_value = sum(minhash_tanimoto_values) / length(
                        minhash_tanimoto_values)

                    # test whether the average tanimoto similarity is off less than 5 % 
                    # compared to the reference values
                    @test avg_minhash_tanimoto_value ≈ ref_minhash_tanimoto_values[i] rtol=0.05
                    
                end
            end
        end
        
    end

    @testset "fingerprint method tests" begin
    
        # We go through combinations of different input parameters and verify the correct
        # length of the fingerprint, and that the resulting fingerprints are not identical

        for n_permutations in [512, 2048]
            fingerprint_results = []
            for seed in [42, (1 << 31)]
                for radius in 1:4
                    for min_radius in 0:min(2, radius)
                        for rings in [true, false]

                            mhfp_calc = MHFP(
                                radius, 
                                min_radius, rings, 
                                n_permutations=n_permutations, 
                                seed=seed)

                            fp = fingerprint(mol, mhfp_calc)
                            push!(fingerprint_results, fp)

                            # verify correct length
                            @test length(fp) == n_permutations
                        end
                    end
                end
            end
            # verify that the fingerprints differ with different parameters
            @test unique(fingerprint_results) == fingerprint_results
        end
    end

    @testset "Giving invalid input tests" begin
        ##### Testing MHFP constructor ####################################################
        @test_throws "must be non-negative" MHFP(-3)  # giving negative radius
        @test_throws "must be non-negative" MHFP(2, -1)  # giving negative min_radius
        @test_throws "must be larger or equal" MHFP(2, 3)  # giving min_radius >= radius
        # giving non-positive n_permutations
        @test_throws "must be strictly positive" MHFP(n_permutations = 0)

        ##### Testing smiles_from_circular_substructures ##################################
        test_mol = smilestomol("c1cc(C(O)=O)ccc1")
        
        # giving non-positive radius
        @test_throws "must be strictly positive" smiles_from_circular_substructures(
            test_mol,
            0,  # radius
            0   # min_radius
        )

        # giving non-positive min_radius
        @test_throws "must be strictly positive" smiles_from_circular_substructures(
            test_mol,
            3,  # radius
            0   # min_radius
        )


    end
    
end
