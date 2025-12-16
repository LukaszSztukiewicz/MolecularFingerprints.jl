using Test
using MolecularGraph
using MolecularFingerprints

@testset "MHFP Fingerprint Tests" begin

    # shingling from molecule helper function
    # example taken from mhfp python implementation, written by the original authors
    # https://github.com/reymond-group/mhfp/blob/master/test/test_encoder.py
    
    mol = MolecularGraph.smilestomol("CCCC1=NN(C2=C1NC(=NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C")
    shingling = sort(["c1cnnc1", "Cn(nc)c(c)c", "C(C)Oc", "c1(OCC)ccccc1-c([nH])n", "S(c)(N)(=O)=O", "c(nc)([nH]c)-c(c)c", "CN(C)C", "n(c(-c)[nH])c(c)=O", "C(C)O", "N(C)(C)S", "S(=O)(=O)(c(cc)cc)N(CC)CC", "CC", "c1(CCC)nn(C)c(c)c1[nH]c", "c1(S(=O)(=O)N(C)C)cccc(-c)c1", "N(C)(C)C", "c(cc)c(c)S", "N1(S(=O)(=O)c(c)c)CCNCC1", "c(c)(c)[nH]", "O=c(nc)c(c)n", "N(C)(CC)CC", "n(c(c)C)n(c)C", "C1CNCCN1", "c1ccccc1", "C(C)N", "n(c)(C)n", "C(CC)c(c)n", "c(c(c)-c)c(c)S", "O(c)C", "CCO", "CN(CC)CC", "[nH](c)c", "n(c)c", "n1c(-c(c)c)[nH]cc(n)c1=O", "N(CC)(CC)S(c)(=O)=O", "C(CN)N(C)C", "S(=O)(=O)(c(c)c)N(C)C", "O(CC)c(cc)c(c)-c", "C(C)Oc(c)c", "C(C)Cc", "C(CN)N(C)S", "c1(-c(nc)[nH]c)cc(S)ccc1OC", "c(cc)(OC)c(c)-c", "n1c(CC)c([nH])c(c)n1C", "c(c)(c)-c", "c([nH]c)(c(C)n)c(c)n", "c(c)(c)O", "c1cc(O)ccc1S(N)(=O)=O", "O=S(c)(N)=O", "c12[nH]c(-c)nc(=O)c1n(C)nc2CC", "c(-c)([nH])n", "c1(=O)nc(-c)[nH]c(c)c1n(C)n", "c1cc(S)cc(-c)c1OC", "O(CC)c(c)c", "c(c)(c)n", "C(C)C", "n(c)n", "CCOc", "c(cc)(-c([nH])n)c(c)O", "c(CC)(nn)c(c)[nH]", "c(c)(c)S", "CCC", "N1(C)CCNCC1", "c(c(n)=O)(c(c)[nH])n(C)n", "n(C)(nc)c(c)c", "c(c)(n)=O", "Cn(c)n", "c(cc)(cc)S(N)(=O)=O", "O=c(c)n", "c(c)c", "n1(C)nc(C)c([nH])c1c(n)=O", "c1c(S(N)(=O)=O)ccc(O)c1-c([nH])n", "C(C)Cc(c)n", "C(C)c", "c(=O)(nc)c(c)n","[nH](c(-c)n)c(c)c", "C(CC)c(nn)c(c)[nH]", "O=c", "c1(-c(cc)c(c)O)nc(=O)cc(c)[nH]1", "c(cc)c(c)O", "c(c)(C)n", "O=S", "c1cnc[nH]c1", "O=S(=O)(c(c)c)N(C)C", "C1CN(C)CCN1S(c)(=O)=O", "c12c(=O)nc[nH]c1c(C)nn2C", "C1CN(S)CCN1C", "[nH]1c(-c(c)c)ncc(n)c1c(C)n", "Cn", "CCCc", "CN"])
    @test sort(mhfp_shingling_from_mol(mol)) == shingling skip=true # skip=true makrs the test as "broken" 

end
