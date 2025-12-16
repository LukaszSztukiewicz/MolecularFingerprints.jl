module MolecularFingerprints

include("interface.jl")
include("utils/tanimoto.jl")
include("algorithms/mhfp.jl")
include("algorithms/ecfp.jl")
include("algorithms/maccs.jl")
include("algorithms/torsions.jl")

export tanimoto, mhfp_shingling_from_mol, fingerprint

end
