module MolecularFingerprints

include("interface.jl")

include("utils/tanimoto.jl")
include("algorithms/mhfp.jl")
include("algorithms/ecfp.jl")
include("algorithms/maccs.jl")
include("algorithms/torsions.jl")

#FIXME in the end, probably we shouldn't export mhfp_shingling_from_mol
export tanimoto, mhfp_shingling_from_mol, fingerprint

end
