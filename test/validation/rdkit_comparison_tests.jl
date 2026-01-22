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
function load_test_data(file_path::String)
    # Load BACE from local CSV
    df = CSV.read(file_path, DataFrame)
    db = Vector{String}(df.mol)

    edge_cases = [
        "C[C@H](F)Cl", "C[C@@H](F)Cl",        # Chirality
        "C1CCCCCCCCCCCCCCCCCCCC1",            # Macrocycle
        "CC(C)(C)C(C)(C)C(C)(C)C",            # Branching
        "[O-]S(=O)(=O)[O-].[Mg+2]",           # Ions
        "c1ccccc1", "C1=CC=CC=C1",            # Aromaticity
        "CC[Se]CC", "B1OC(C)CC1"              # Heteroatoms
    ]
    append!(db, edge_cases)
    return db
end

# --- 4. RDKIT REFERENCE FUNCTION ---
const rdFPGen = pyimport("rdkit.Chem.rdFingerprintGenerator")
function get_rdkit_scores(query::String, db::Vector{String}, r, n)
    q_mol = Chem.MolFromSmiles(query)
    # Filter out invalid SMILES (which return None in Python) to prevent crashes
    db_mols = [m for m in (Chem.MolFromSmiles(s) for s in db) if !pyisnone(m)]
    
    # Use the explicit rdFingerprintGenerator module
    gen = rdFPGen.GetMorganGenerator(radius=r, fpSize=n)
    
    # CHANGE: GetFingerprintAsBitVect -> GetFingerprint
    q_fp = gen.GetFingerprint(q_mol)
    db_fps = [gen.GetFingerprint(m) for m in db_mols]
    
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
    database = load_test_data(csv_path)
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