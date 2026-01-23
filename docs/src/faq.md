# Frequently Asked Questions (FAQ)

### What is a molecular fingerprint?
A molecular fingerprint is a fixed-length vector representation of a molecule, where each bit indicates the presence or absence of specific structural features or patterns within the molecule. This representation allows for efficient computation, similarity searching, and machine learning applications in cheminformatics. More extensive explanations can be found in the [Explanation](explanation.md) section of the documentation.

### How do I choose the right fingerprint type for my application?
The choice of fingerprint type depends on your specific application and the characteristics of the molecules you are working with. Common types include ECFP (Extended Connectivity Fingerprints) for capturing local atomic environments, MACCS keys for substructure presence, and path-based fingerprints for capturing linear paths in the molecular graph. Consider the size of your dataset, the nature of your molecules, and the computational resources available when selecting a fingerprint type.

### Can I customize the size and parameters of the fingerprints?
Yes, most fingerprint types in MolecularFingerprints.jl allow you to customize parameters such as the size of the fingerprint vector and specific algorithm settings (e.g., radius for ECFP). Refer to the documentation for details on how to configure these parameters. See the [API Reference](api_reference.md) for more information.

### Will this package provide descriptors in addition to fingerprints?
MolecularFingerprints.jl is primarily focused on generating molecular fingerprints. However, it is designed to be compatible with other Julia packages that provide molecular descriptors, allowing you to use both fingerprints and descriptors in your analyses.