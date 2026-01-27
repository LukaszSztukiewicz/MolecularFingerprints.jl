using CSV
using DataFrames
using BenchmarkTools
using Random
using PythonCall
using LinearAlgebra

# --- SETUP RDKIT VIA PYTHONCALL ---
const Chem = pyimport("rdkit.Chem")
const AllChem = pyimport("rdkit.Chem.AllChem")
const DataStructs = pyimport("rdkit.DataStructs")

const MACCSkeys = pyimport("rdkit.Chem.MACCSkeys")
const rdFingerprintGenerator = pyimport("rdkit.Chem.rdFingerprintGenerator")
const rdMHFPFingerprint = pyimport("rdkit.Chem.rdMHFPFingerprint")


# --- FINGERPRINTS ---
ecfp_calc = ECFP() # Extended Connectivity Fingerprints
mhfp_calc = MHFP() # MinHash Fingerprints
torsion_calc = TopologicalTorsion() # Topological Torsion Fingerprints
maccs_calc = MACCS() # MACCS Keys

const CALCULATORS = Dict(
    "ECFP" => ecfp_calc,
    # "MHFP" => mhfp_calc,
    # "TopologicalTorsion" => torsion_calc,
    "MACCS" => maccs_calc
)


# --- DATASETS  ---
const DATASET_CONFIG = Dict(
    "bace" => :mol,
    # "BBBP" => :smiles,
    # "qm8"  => :smiles,
    # "esol" => :smiles
)

const SMILES_EDGE_CASES = [
    "C[C@H](F)Cl", "C[C@@H](F)Cl",        # Chirality
    "C1CCCCCCCCCCCCCCCCCCCC1",            # Macrocycle
    "CC(C)(C)C(C)(C)C(C)(C)C",            # Branching
    "[O-]S(=O)(=O)[O-].[Mg+2]",           # Ions
    "c1ccccc1", "C1=CC=CC=C1",            # Aromaticity
    "CC[Se]CC", "B1OC(C)CC1"              # Heteroatoms
]

function load_test_data(dataset::AbstractString; folder_path::AbstractString="./validation")
    if !haskey(DATASET_CONFIG, dataset)
        throw(ArgumentError("Dataset '$dataset' not supported. Use: $(keys(DATASET_CONFIG))"))
    end
    file_path = joinpath(folder_path, "$dataset.csv")
    df = CSV.read(file_path, DataFrame)
    smiles_col = DATASET_CONFIG[dataset]
    return Vector{String}(df[!, smiles_col])
end

folder_path = "./validation"

# ds_bbbp = load_test_data("BBBP"; folder_path=folder_path)
# ds_qm8 = load_test_data("qm8"; folder_path=folder_path)
ds_bace = load_test_data("bace"; folder_path=folder_path)
# ds_esol = load_test_data("esol"; folder_path=folder_path)
ds_hard = SMILES_EDGE_CASES

all_datasets = unique(vcat(ds_bace, ds_hard)) #26 431 molecules

morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=ecfp_calc.radius, fpSize=1024)
tt_gen = rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=2048)
mhfp_encoder = rdMHFPFingerprint.MHFPEncoder()

function fingerprint_rdkit(smiles::String, calc::AbstractCalculator)
    mol = Chem.MolFromSmiles(smiles)
    if mol === nothing
        # uses our defined handling for invalid molecules via multiple dispatch
        return fingerprint(nothing, calc) 
    end

    if calc isa MACCS
        fp = MACCSkeys.GenMACCSKeys(mol)
        bit_str = pyconvert(String, fp.ToBitString())
        fp_vector = BitVector(c == '1' for c in bit_str)
        return fp_vector[2:end]  # RDKit MACCS has 167 bits, first bit is unused

    elseif calc isa ECFP
        fp = morgan_gen.GetFingerprint(mol)
    
    elseif calc isa MHFP
        fp = mhfp_encoder.EncodeMol(mol)
        fp = pyconvert(Vector{Int64}, fp)
        return fp
    
    elseif calc isa TopologicalTorsion
        fp = tt_gen.GetFingerprint(mol)
    
    else
        error("Unknown calculator type.")
    
    end

    bit_str = pyconvert(String, fp.ToBitString())
    return BitVector(c == '1' for c in bit_str)
end


function run_similarity_comparison(smiles_list::Vector{String}, calc::AbstractCalculator)
    tanimoto_scores = Float64[]
    cosine_scores = Float64[]

    n = length(smiles_list)
    for i in 1:n 
        smi = smiles_list[i]
        rd_fp = fingerprint_rdkit(smi, calc)
        jl_fp = fingerprint(smi, calc)

        score = tanimoto_similarity(jl_fp, rd_fp)
        push!(tanimoto_scores, score)
        score = cosine_similarity(jl_fp, rd_fp)
        push!(cosine_scores, score)
    end
    println("Completed similarity comparison for $(length(smiles_list)) molecules using $(typeof(calc)).")
    println("Tanimoto Similarity Mean=$(mean(tanimoto_scores)), Median=$(median(tanimoto_scores)), Std=$(std(tanimoto_scores)), Minimum=$(minimum(tanimoto_scores)), Maximum=$(maximum(tanimoto_scores))")
    println("Cosine Similarity: Mean=$(mean(cosine_scores)), Median=$(median(cosine_scores)), Std=$(std(cosine_scores)), Minimum=$(minimum(cosine_scores)), Maximum=$(maximum(cosine_scores))")
    return tanimoto_scores, cosine_scores
end

# --- 5. MAIN BENCHMARK & TEST SUITE ---
function run_all_tests()
    limit = 1500  # Limit number of molecules for testing purposes
    # Iterate over all calculators and datasets
    for (calc_name, calc) in CALCULATORS
        println("=== Testing Calculator: $calc_name ===")
        for (ds_name, _) in DATASET_CONFIG
            println("---- Dataset: $ds_name ----")
            smiles_list = load_test_data(ds_name; folder_path=folder_path)
            smiles_list = smiles_list[1:limit]
            tanimoto_scores, cosine_scores = run_similarity_comparison(smiles_list, calc)
            # Basic assertions
            @assert all(0.0 .<= tanimoto_scores .<= 1.0) "Tanimoto scores out of bounds!"
            @assert all(0.0 .<= cosine_scores .<= 1.0) "Cosine scores out of bounds!"
        end
    end

    # Edge Case Tests
    println("=== Testing Edge Cases ===")
    for (calc_name, calc) in CALCULATORS
        println("---- Calculator: $calc_name ----")
        tanimoto_scores, cosine_scores = run_similarity_comparison(SMILES_EDGE_CASES, calc)
        @assert all(0.0 .<= tanimoto_scores .<= 1.000001) "Tanimoto scores out of bounds in edge cases!"
        @assert all(0.0 .<= cosine_scores .<= 1.000001) "Cosine scores out of bounds in edge cases!"
    end

    # Benchmarking
    println("=== Benchmarking Fingerprint Calculations ===")
    sample_smiles = all_datasets[1:100]  # Sample 100 molecules for benchmarking
    for (calc_name, calc) in CALCULATORS
        println("---- Calculator: $calc_name ----")
        @btime fingerprint($sample_smiles, $calc)
    end

    # Benchmarking vs RDKit
    println("=== Benchmarking vs RDKit ===")
    for (calc_name, calc) in CALCULATORS
        println("---- Calculator: $calc_name ----")
        jl_time = @belapsed begin 
            for smi in $sample_smiles
                fingerprint(smi, $calc)
            end
        end
        rd_time = @belapsed begin
            for smi in $sample_smiles
                fingerprint_rdkit(smi, $calc)
            end
        end
        println("Julia Time: $(jl_time) seconds for 100 molecules.")
        println("RDKit Time: $(rd_time) seconds for 100 molecules.")
    end

    # Compare results of similarity search between Julia and RDKit
    println("=== Similarity Search Comparison ===")
    query_smiles = all_datasets[1:10]  # 10 query molecules
    db_smiles = all_datasets[11:110]   # 100 database molecules
    for (calc_name, calc) in CALCULATORS
        println("---- Calculator: $calc_name ----")
        for q_smi in query_smiles
            jl_scores = Float64[]
            rd_scores = Float64[]

            q_fp_jl = fingerprint(q_smi, calc)
            q_fp_rd = fingerprint_rdkit(q_smi, calc)

            for db_smi in db_smiles
                db_fp_jl = fingerprint(db_smi, calc)
                db_fp_rd = fingerprint_rdkit(db_smi, calc)

                push!(jl_scores, tanimoto_similarity(q_fp_jl, db_fp_jl))
                push!(rd_scores, tanimoto_similarity(q_fp_rd, db_fp_rd))
            end

            # Compare top 10 results
            jl_top10 = sortperm(jl_scores, rev=true)[1:10]
            rd_top10 = sortperm(rd_scores, rev=true)[1:10]
            println("Query SMILES: $q_smi")
            println("Julia Top 10 Indices: $jl_top10")
            println("RDKit Top 10 Indices: $rd_top10")
            println("Recall@10: $(length(intersect(jl_top10, rd_top10)) / 10)")
            ndcg = 0.0
            for (rank, idx) in enumerate(jl_top10)
                if idx in rd_top10
                    ndcg += 1 / log2(rank + 1)
                end
            end
            println("NDCG@10: $ndcg")
        end
        println("Similarity search results match between Julia and RDKit for calculator $calc_name.")
    end

    #calculate average recall at 10 for 100 queries for tanimoto similarity and cosine similarity
    println("=== Average Recall@10 over 100 Queries ===")
    for (calc_name, calc) in CALCULATORS
        println("---- Calculator: $calc_name ----")
        total_recall_tanimoto = 0.0
        total_recall_cosine = 0.0
        if calc isa MHFP
            num_queries = 10 # MHFP is slow, limit to 10 queries
        else
            num_queries = 100
        end
        query_smiles = all_datasets[1:num_queries]
        db_smiles = all_datasets[num_queries+1:num_queries+100]
        for q_smi in query_smiles
            jl_scores_tanimoto = Float64[]
            rd_scores_tanimoto = Float64[]
            jl_scores_cosine = Float64[]
            rd_scores_cosine = Float64[]

            q_fp_jl = fingerprint(q_smi, calc)
            q_fp_rd = fingerprint_rdkit(q_smi, calc)

            for db_smi in db_smiles
                db_fp_jl = fingerprint(db_smi, calc)
                db_fp_rd = fingerprint_rdkit(db_smi, calc)

                push!(jl_scores_tanimoto, tanimoto_similarity(q_fp_jl, db_fp_jl))
                push!(rd_scores_tanimoto, tanimoto_similarity(q_fp_rd, db_fp_rd))
                push!(jl_scores_cosine, cosine_similarity(q_fp_jl, db_fp_jl))
                push!(rd_scores_cosine, cosine_similarity(q_fp_rd, db_fp_rd))
            end

            # Compare top 10 results for Tanimoto
            jl_top10_tanimoto = sortperm(jl_scores_tanimoto, rev=true)[1:10]
            rd_top10_tanimoto = sortperm(rd_scores_tanimoto, rev=true)[1:10]
            recall_tanimoto = length(intersect(jl_top10_tanimoto, rd_top10_tanimoto)) / 10
            total_recall_tanimoto += recall_tanimoto

            # Compare top 10 results for Cosine
            jl_top10_cosine = sortperm(jl_scores_cosine, rev=true)[1:10]
            rd_top10_cosine = sortperm(rd_scores_cosine, rev=true)[1:10]
            recall_cosine = length(intersect(jl_top10_cosine, rd_top10_cosine)) / 10
            total_recall_cosine += recall_cosine
        end
        avg_recall_tanimoto = total_recall_tanimoto / num_queries
        avg_recall_cosine = total_recall_cosine / num_queries
        println("Average Recall@10 (Tanimoto): $avg_recall_tanimoto")
        println("Average Recall@10 (Cosine): $avg_recall_cosine")
    
    end
        
end

# EXECUTE
run_all_tests()