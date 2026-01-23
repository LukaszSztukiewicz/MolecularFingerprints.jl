using CSV
using DataFrames
using BenchmarkTools
using Random
using PythonCall
using LinearAlgebra

# --- 1. SETUP RDKIT VIA PYTHONCALL ---
const Chem = pyimport("rdkit.Chem")
const AllChem = pyimport("rdkit.Chem.AllChem")
const DataStructs = pyimport("rdkit.DataStructs")

# --- 3. DATA LOADING & EDGE CASES ---
const DATASET_CONFIG = Dict(
    "bace" => :mol,
    "BBBP" => :smiles,
    "qm8"  => :smiles,
    "esol" => :smiles
)

const SMILES_EDGE_CASES = [
    "C[C@H](F)Cl", "C[C@@H](F)Cl",        # Chirality
    "C1CCCCCCCCCCCCCCCCCCCC1",            # Macrocycle
    "CC(C)(C)C(C)(C)C(C)(C)C",            # Branching
    "[O-]S(=O)(=O)[O-].[Mg+2]",           # Ions
    "c1ccccc1", "C1=CC=CC=C1",            # Aromaticity
    "CC[Se]CC", "B1OC(C)CC1"              # Heteroatoms
]

"""
    load_test_data(dataset::AbstractString; folder_path="./validation")

Loads SMILES data from a CSV and appends standard edge cases.
"""
function load_test_data(dataset::AbstractString; folder_path::AbstractString="./validation")
    # 1. Validation using the Dict keys
    if !haskey(DATASET_CONFIG, dataset)
        throw(ArgumentError("Dataset '$dataset' not supported. Use: $(keys(DATASET_CONFIG))"))
    end

    # 2. Use joinpath for OS-agnostic path handling
    file_path = joinpath(folder_path, "$dataset.csv")

    # 3. Read and extract column
    df = CSV.read(file_path, DataFrame)
    smiles_col = DATASET_CONFIG[dataset]
    
    # 4. Convert and combine
    # Using [a; b] is a concise way to call vcat()
    return [Vector{String}(df[!, smiles_col]); SMILES_EDGE_CASES]
end

# --- 4. RDKIT REFERENCE FUNCTION ---
function get_rdkit_scores(query::String, db::Vector{String}, r, n)
    q_mol = Chem.MolFromSmiles(query)
    db_mols = [Chem.MolFromSmiles(s) for s in db]
    
    q_fp = AllChem.GetMorganFingerprintAsBitVect(q_mol, r, nBits=n)
    db_fps = [AllChem.GetMorganFingerprintAsBitVect(m, r, nBits=n) for m in db_mols]
    
    py_scores = DataStructs.BulkTanimotoSimilarity(q_fp, db_fps)
    return pyconvert(Vector{Float64}, py_scores)
end

# --- 5. MAIN BENCHMARK & TEST SUITE ---
function run_all_tests()
    #print current working directory
    #println("Current working directory: ", pwd())

    println("=== Chemical Fingerprint Test Suite ===")
    
    # Configuration
    csv_path = "./validation/bace.csv"
    query_smiles = "O=C(O)C1=CC=CC=C1" # Benzoic Acid
    radius = 2
    nbits = 2048
    calc = ECFP{nbits}(radius)

    # 1. Load Data
    if !isfile(csv_path)
        error("bace.csv not found in folder. Please ensure the file is present.")
    end
    database = load_test_data("bace")
    #take first 10 molecules for performance
    database = database[1:10]
    println("Database loaded with $(length(database)) molecules (BACE + Edge Cases).")

    # 2. Performance Benchmark (Julia)
    println("\n--- Performance: Julia Fingerprinting ---")
    # Warmup
    fingerprint(database, calc)
    # Benchmark
    t_julia = @benchmark fingerprint($database[1:10], $calc)
    display(t_julia)

    # 3. Accuracy Check
    println("\n--- Accuracy: Julia vs RDKit ---")
    
    # Julia Results
    jl_fps = fingerprint(database, calc)
    q_fp_jl = fingerprint(query_smiles, calc)
    julia_scores = [tanimoto_similarity(q_fp_jl, fp) for fp in jl_fps]

    # RDKit Results
    rdkit_scores = get_rdkit_scores(query_smiles, database, radius, nbits)

    # Comparison Logic
    k = 10
    jl_top_idx = sortperm(julia_scores, rev=true)[1:k]
    rd_top_idx = sortperm(rdkit_scores, rev=true)[1:k]

    common_hits = intersect(jl_top_idx, rd_top_idx)
    println("Similarity in Top $k rankings: $(length(common_hits)) / $k")
    
    if length(common_hits) > 0
        println("Sample common hit: ", database[common_hits[1]])
    end

    # Comparison Summary Table
    println("\nTop 3 Comparison:")
    println("-"^50)
    printf_fmt = "%-20s | %-15s | %-15s\n"
    @info "Rank" "Julia Index" "RDKit Index"
    for i in 1:3
        println("Rank $i: Julia #$(jl_top_idx[i]) vs RDKit #$(rd_top_idx[i])")
    end
end

# EXECUTE
run_all_tests()