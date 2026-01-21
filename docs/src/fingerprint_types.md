# MolecularFingerprints.jl - Fingerprint Types

This section breaks down the specific algorithms we have implemented in this package. As a computer scientist, it helps to think of these as different ways to hash a graph. Depending on which structural features you want to emphasize—local neighborhoods, specific functional groups, or long-range connectivity—you might choose one over the other.

## Extended Connectivity Fingerprints (ECFP)

ECFPs are the industry standard for most tasks today. In the CS world, you can think of these as a variation of the Weisfeiler-Lehman graph isomorphism test. The algorithm works by looking at each atom and its immediate neighbors, assigning them an initial integer identifier. It then iteratively updates these identifiers by looking at neighbors at increasing distances—usually referred to as the radius. Once the iterations are complete, all the unique identifiers are hashed into a bit vector of a fixed size, like 1024 or 2048. Because they capture local circular environments, they are incredibly effective at identifying similar molecular "building blocks" regardless of where they appear in the molecule.

## MACCS Keys

If ECFPs are like a dynamic hashing algorithm, MACCS keys are more like a checklist. This is a dictionary-based approach where each bit in the vector corresponds to a specific, predefined chemical question. For instance, bit 160 might ask: "Is there at least one oxygen atom present?" and bit 161 might ask: "Is there a nitrogen-hydrogen bond?". There are 166 of these keys in total. While this is a older and much more rigid approach than modern hashing methods, it remains very popular because the results are highly interpretable. You know exactly what each bit represents, which is not usually the case with hashed fingerprints.

## MHFP (Molecular Hash Fingerprints)

MHFP is a more recent development designed specifically for very large-scale datasets. It borrows heavily from the MinHash technique used in natural language processing and document deduplication. Instead of just hashing structural fragments into a bit vector, MHFP uses a MinHash-based scheme to create a signature that is very efficient for locality-sensitive hashing (LSH). This makes it possible to perform similarity searches across billions of molecules in sub-second time. If you are building a search engine for chemical space, this is likely the algorithm you would reach for.

## Topological Torsion Fingerprints

Topological torsion fingerprints take a linear approach rather than a circular one. Instead of looking at neighborhoods, they identify all sequences of four bonded atoms—which we call torsions—and record the types of atoms, their bond orders, and the number of non-hydrogen neighbors. These paths are then hashed into the fingerprint. This method is particularly sensitive to the overall shape and "skeleton" of the molecule. It is often used as a complementary tool to ECFP because it captures the long-range connectivity of the molecular graph that circular methods might overlook.
