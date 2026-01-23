# Similarity Search: Finding Chemically Related Molecules

In this tutorial, we will find the most chemically similar molecules in a database relative to a "query" molecule.

---

## 1. Setup and Library Loading

We will use the same core libraries: `MolecularFingerprints.jl` for generating descriptors and `MolecularGraph.jl` for parsing molecules.

```julia
using LinearAlgebra
using MolecularFingerprints
using MolecularGraph

# A small "database" of molecules
database_smiles = [
    "CCO",                # Ethanol
    "CC(=O)Oc1ccccc1C(=O)O", # Aspirin
    "c1ccccc1",           # Benzene
    "CCCCCCCC",           # Octane
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", # Caffeine
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"  # Ibuprofen
]

# Our Query Molecule: Acetaminophen (Paracetamol)
query_smiles = "CC(=O)Nc1ccc(O)cc1"

```

---

## 2. Fingerprinting the Database

To compare molecules, we must first map them into the same vector space. We'll use **Morgan Fingerprints (ECFP4)**, which capture the local neighborhood of atoms.

```julia
# Convert database to molecules and then to fingerprints
db_mols = [smilestomol(s) for s in database_smiles]
featurizer = ECFP{2}(2) # Radius 2 (ECFP4 equivalent)

# Generate a list of BitVectors
db_fps = [fingerprint(m, featurizer) for m in db_mols]

# Process the Query Molecule
query_mol = smilestomol(query_smiles)
query_fp = fingerprint(query_mol, featurizer)

```

---

## 3. Calculating Tanimoto Similarity

In cheminformatics, the **tanimoto_similarity Coefficient** (or Jaccard Index) is the industry standard for comparing bit-vector fingerprints. It measures the intersection of features divided by the union.

A score of **1.0** means the fingerprints are identical; **0.0** means they share no common features.

```julia
# Function to calculate tanimoto_similarity between two BitVectors
function tanimoto_similarity_similarity(fp1, fp2)
    intersection = sum(fp1 .& fp2)
    union = sum(fp1 .| fp2)
    return intersection / union
end

# Calculate similarity between query and every item in the database
similarities = [tanimoto_similarity(query_fp, db_fp) for db_fp in db_fps]

# Pair the SMILES with their scores for easy reading
results = collect(zip(database_smiles, similarities))

```

---

## 4. Ranking and Results

Finally, we sort our database based on the similarity score to find the closest matches.

```julia
# Sort by similarity score in descending order
sort!(results, by = x -> x[2], rev = true)

println("Similarity Search Results for Query: $query_smiles")
println("-"^45)
for (smiles, score) in results
    println("Score: $(round(score, digits=3)) | Molecule: $smiles")
end

```