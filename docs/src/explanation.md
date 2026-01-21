# MolecularFingerprints.jl

## Explanation

### Why we need fingerprints

Most chemical data is stored as SMILES strings, like `CCO` for ethanol or `C1=CC=CC=C1` for benzene. While these are great for humans to read and for databases to store, they are pretty difficult to use for actual computation. You cannot perform vector math on a string, and standard string distance metrics like Levenshtein distance do not correspond to chemical similarity. A tiny change in a molecule’s structure can result in a completely different string, which makes searching and machine learning nearly impossible (see the cliassic example of `CCO` vs `COC` - ethanol and dimethyl ether—which are structurally quite different but only one Levenshtein edit apart).

This is why we use molecular fingerprints. We essentially take the molecular graph and transform it into a high-dimensional vector. This moves the chemistry from the world of linguistics into the world of linear algebra.

### How they work

A molecular fingerprint is typically a fixed-length vector where each bit represents the presence or absence of a specific structural feature. We generate these through a few different strategies:

* Path-based methods traverse the molecular graph and identify all possible paths of a certain length, hashing those paths into the bit vector.
* Circular methods look at the local neighborhood around each individual atom, iteratively expanding outward to capture the local environment.
* Substructure methods check the molecule against a predefined library of chemical motifs, like searching for a specific regex in a block of text.

### Practical applications

#### Similarity searching
Once you have converted a molecule into a bit vector, you can use concepts you are already familiar with in Computer Science. For example, we use the Tanimoto similarity (which is just another name for the Jaccard similarity) to calculate how similar two molecules are. This allows us to perform lightning-fast similarity searches across databases containing millions of compounds.

#### Machine learning
These vectors also serve as the standard input for machine learning models. If you want to predict whether a molecule is toxic or if it will bind to a specific protein, these fingerprints provide the fixed-length feature set you need for a random forest, a support vector machine, or a neural network.

### The Julia advantage

We built this package in Julia because fingerprinting often involves processing massive datasets where performance is critical. By using Julia, we get the execution speed of C++ while maintaining the high-level syntax needed for complex data science workflows. This library is designed to be a fast, type-safe way to integrate chemical data into your computational pipelines.