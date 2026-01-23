using CSV
using DataFrames
using BenchmarkTools
using Random
using PythonCall
using LinearAlgebra
using Statistics

# --- RDKIT SETUP ---

const Chem = pyimport("rdkit.Chem")
const AllChem = pyimport("rdkit.Chem.AllChem")
const DataStructs = pyimport("rdkit.DataStructs")
const MACCSkeys = pyimport("rdkit.Chem.MACCSkeys")
const rdFingerprintGenerator = pyimport("rdkit.Chem.rdFingerprintGenerator")
const rdMHFPFingerprint = pyimport("rdkit.Chem.rdMHFPFingerprint")


# --- GENERATORS (Lazy Initialization) ---
const MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
const TT_GEN = rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=2048)
const MHFP_ENCODER = rdMHFPFingerprint.MHFPEncoder()

# --- RDKIT INTERFACE (Multiple Dispatch) ---

"""
    fingerprint_rdkit(smiles::String, calc::AbstractCalculator)

Internal dispatch mechanism to handle different RDKit fingerprint types.
"""
function fingerprint_rdkit(smiles::String, calc::AbstractCalculator)
    mol = Chem.MolFromSmiles(smiles)
    pyisnone(mol) && return nothing
    return _compute_rdkit_fp(mol, calc)
end

function _compute_rdkit_fp(mol, ::MACCS)
    fp = MACCSkeys.GenMACCSKeys(mol)
    # Convert RDKit bitstring to BitVector, skipping bit 0 (RDKit specific)
    bit_str = pyconvert(String, fp.ToBitString())
    return BitVector(c == '1' for c in bit_str[2:end])
end

function _compute_rdkit_fp(mol, ::ECFP)
    fp = MORGAN_GEN.GetFingerprint(mol)
    return _py_bitvect_to_julia(fp)
end

function _compute_rdkit_fp(mol, ::TopologicalTorsion)
    fp = TT_GEN.GetFingerprint(mol)
    return _py_bitvect_to_julia(fp)
end

function _compute_rdkit_fp(mol, ::MHFP)
    fp = MHFP_ENCODER.EncodeMol(mol)
    return pyconvert(Vector{Int64}, fp)
end

# Helper to convert RDKit ExplicitBitVect to Julia BitVector
function _py_bitvect_to_julia(py_fp)
    bit_str = pyconvert(String, py_fp.ToBitString())
    return BitVector(c == '1' for c in bit_str)
end

# --- DATA UTILITIES ---

const DATASET_CONFIG = Dict(
    "bace" => :mol,
    # "BBBP" => :smiles,
    # "qm8"  => :smiles,
    # "esol" => :smiles
)

function load_test_data(dataset::String; folder="./validation")
    haskey(DATASET_CONFIG, dataset) || throw(ArgumentError("Unknown dataset: $dataset"))
    
    path = joinpath(folder, "$dataset.csv")
    df = CSV.read(path, DataFrame)
    return Vector{String}(df[!, DATASET_CONFIG[dataset]])
end

# --- ANALYSIS ENGINE ---

function run_similarity_comparison(smiles_list::Vector{String}, calc::AbstractCalculator)
    # Filter out invalid molecules first to avoid mid-loop failures
    valid_smiles = filter(s -> !pyisnone(Chem.MolFromSmiles(s)), smiles_list)
    
    results = map(valid_smiles) do smi
        rd_fp = fingerprint_rdkit(smi, calc)
        jl_fp = fingerprint(smi, calc) # Assumed native Julia implementation
        
        (
            tanimoto = tanimoto_similarity(jl_fp, rd_fp),
            cosine   = cosine_similarity(jl_fp, rd_fp)
        )
    end

    t_scores = [r.tanimoto for r in results]
    c_scores = [r.cosine for r in results]

    print_stats(calc, t_scores, c_scores)
    return t_scores, c_scores
end

function print_stats(calc, t_scores, c_scores)
    @info "Results for $(typeof(calc))"
    println("Tanimoto | Mean: $(round(mean(t_scores), digits=4)) | Min: $(round(minimum(t_scores), digits=4))")
    println("Cosine   | Mean: $(round(mean(c_scores), digits=4)) | Min: $(round(minimum(c_scores), digits=4))")
end

# --- TEST SUITE ---

function run_validation_suite()
    calculators = [ECFP(2), MACCS()]
    datasets = ["bace"]
    limit = 100

    for calc in calculators
        @info "=== Testing Calculator: $(typeof(calc)) ==="
        
        for ds_name in datasets
            smiles = load_test_data(ds_name)[1:limit]
            t_scores, c_scores = run_similarity_comparison(smiles, calc)
            
            # Use @assert for logic, but prefer explicit checks for validation
            if !all(0 .<= t_scores .<= 1)
                @error "Tanimoto scores out of bounds for $ds_name"
            end
        end
    end

    # --- Benchmarking vs RDKit ---
    @info "=== Benchmarking vs RDKit ==="
    sample_smiles = load_test_data("esol")[1:100]
    calc = ECFP(2)

    jl_time = @belapsed fingerprint($sample_smiles, $calc)
    rd_time = @belapsed begin
        for smi in $sample_smiles
            fingerprint_rdkit(smi, $calc)
        end
    end

    println("Julia is $(round(rd_time/jl_time, digits=2))x faster than RDKit wrapper.")
end
