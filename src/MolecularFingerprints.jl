module MolecularFingerprints

include("utils/tanimoto.jl")
include("algorithms/mhfp.jl")

export tanimoto, mhfp_shingling_from_mol, fingerprint

end
