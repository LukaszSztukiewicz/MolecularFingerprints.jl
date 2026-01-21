# Why MolecularFingerprints.jl?

If you are coming from a computer science background, you are likely used to thinking about data in terms of arrays, hashes, or objects. In chemoinformatics, our primary challenge is representing a real chemical molecule in a format that a computer can actually work with efficiently.

## Why we need fingerprints

Most chemical data is stored as Simplified Molecular Input Line Entry System [(SMILES)](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) strings, like `CCO` for ethanol or `C1=CC=CC=C1` for benzene. While these are great for humans to read and for databases to store, they are pretty difficult to use for actual computation. 

You cannot perform vector math on a string, and standard string distance metrics like Levenshtein distance do not correspond to chemical similarity. A small change in a SMILES string can represent a huge change in the actual molecule, and vice versa. For example, the SMILES strings `CCO` (ethanol) and `CC=O` (acetaldehyde) differ by just one character, but the molecules have very different properties. On the other hand, `CCO` (ethanol) and `C(C)O` (also ethanol) look quite different as strings but represent the same molecule.

This is why we use molecular fingerprints. We essentially take the molecular graph and transform it into a high-dimensional vector.

## How they work

A molecular fingerprint is typically a fixed-length vector where each vector entry represents the presence or absence of a specific structural feature. We generate these through a few different strategies:

* Path-based methods traverse the molecular graph and identify all possible paths of a certain length, hashing those paths into the fingerprint vector.
* Circular methods look at the local neighborhood around each individual atom, iteratively expanding outward to capture the local environment.
* Substructure methods check the molecule against a predefined library of chemical substructures [(motifs)](https://en.wikipedia.org/wiki/Chemical_substructure), like searching for a specific regex in a block of text.

## Practical applications

### Similarity searching
Once you have converted a molecule into a fingerprint vector, you can use concepts you are already familiar with in Computer Science. For example, we use the Tanimoto similarity (which is just another name for the Jaccard similarity) to calculate how similar two molecules are. This allows us to perform lightning-fast similarity searches across databases containing millions of compounds.

### Machine learning
These vectors also serve as the standard input for machine learning models. If you want to predict whether a molecule is toxic or if it will bind to a specific protein, these fingerprints provide the fixed-length feature set you need for a random forest, a support vector machine, or a neural network.

## Why MolecularFingerprints.jl in Julia?

There exist several cheminformatics libraries in other programming languages, such as RDKit in Python or CDK in Java. However, we wanted to create a native Julia implementation to leverage Julia's strengths in scientific computing, performance, and ease of use. Julia's multiple dispatch system allows us to create flexible and extensible fingerprinting algorithms that can be easily integrated into larger cheminformatics workflows. 

MolecularFingerprints.jl is designed to have minimal dependencies, making it lightweight and easy to install. This allows users to quickly get started with molecular fingerprinting without the overhead of large external libraries.

It is also designed to interoperate seamlessly with other Julia packages in the cheminformatics ecosystem, such as MolecularGraph.jl for molecular representation and manipulation. 