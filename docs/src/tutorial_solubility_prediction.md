# QSAR Pipeline: Solubility Prediction

In this tutorial, we will build a Random Forest model to predict the aqueous solubility of molecules using `MolecularFingerprints.jl` and `DecisionTree.jl`.

## 1. Setup and Data Loading

We will use a small subset of the **Delaney (ESOL)** dataset. For this demo, we mock the data loading, but you would normally use `CSV.read`.

```julia
using Random
using Statistics
using DecisionTree
using MolecularFingerprints
using MolecularGraph

# For reproducibility
Random.seed!(42)

# Small subset of the Delaney Solubility Dataset (SMILES, Measured LogS)
data = [
    ("CCO", 0.8),               # Ethanol (High solubility)
    ("CC(=O)Oc1ccccc1C(=O)O", -2.1), # Aspirin (Moderate)
    ("c1ccccc1", -2.0),         # Benzene (Low)
    ("CCCCCCCC", -4.5),         # Octane (Very Low)
    ("O=C(C)Oc1ccccc1C(=O)O", -2.2), # Aspirin analog
    ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", -0.9) # Caffeine
]

# Separate X (SMILES) and y (Labels)
smiles_list = [d[1] for d in data]
y = [d[2] for d in data]
```


## 2. Feature Generation (Fingerprinting)
We convert the raw SMILES strings into fixed-length numeric vectors using ECFP4.

```julia
# Parse SMILES to GraphMol objects
mols = [smilestomol(s) for s in smiles_list]

# Define the featurizer
featurizer_ecfp = ECFP{2}(2) # 2048-bit ECFP with radius 2
featurizer_maccs = MACCS(true, false) # MACCS keys

# Generate BitVectors
# Note: DecisionTree.jl expects a standard Matrix{Float64} or Matrix{Int}
fingerprints = [fingerprint(m, featurizer_ecfp) for m in mols]
# Convert Vector of BitVectors to a Matrix (Samples x Features)
# We transpose (') because hcat stacks them as columns
X = hcat(fingerprints...)'
X = Matrix(X) # Convert to standard dense matrix for ML
```

## 3. Model Training

```julia
# Train/Test Split (Simple 80/20 manual split for demo)
n_samples = length(y)
train_idx = shuffle(1:n_samples)[1:floor(Int, 0.8*n_samples)]
test_idx = setdiff(1:n_samples, train_idx)

X_train, y_train = X[train_idx, :], y[train_idx]
X_test, y_test = X[test_idx, :], y[test_idx]

# Initialize and Train Model
# n_trees=10, min_samples=2
model = build_forest(y_train, X_train, 2, 10, 0.7, -1)

println("Model Trained successfully!")
```

## 4. Model Evaluation

```julia
# Predict on Test Set
preds = apply_forest(model, X_test)
rmse = sqrt(mean((preds .- y_test).^2))

println("Test RMSE: $rmse")

# Inference on New Molecule
new_mol_smiles = "CC(=O)N" # Acetamide
new_mol = smilestomol(new_mol_smiles)
new_fp = fingerprint(new_mol, featurizer_ecfp)

# Reshape to 1xN Matrix for prediction
new_X = reshape(new_fp, 1, length(new_fp))

predicted_solubility = apply_forest(model, new_X)[1]
println("Predicted Solubility for $new_mol_smiles: $predicted_solubility")
```


