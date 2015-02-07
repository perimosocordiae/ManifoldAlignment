# Manifold Alignment

We have two slightly different versions of alignment approaches.
The input to represent the weight matrix for each domain can be given in two forms
1. a sparse matrix. This version is the one that I am currently using.
2. an array modeling the k-nearest neighbor for each instance. The files are in /knn/* 
	EU test is based on this version.

## Folder: Alignment

### Affine Matching

 * AffineMatching.m

### Procrustes Alignment

 * Procrustes.m

### Canonical Correlation Analysis (CCA)

 * CCATwo.m
 * CCAThree.m

### Feature-level Manifold Projections

 * wmapGeneralThree.m
 * wmapGeneralTwo.m

### Instance-level Manifold Projections

 * wmapGeneralThreeInstance.m
 * wmapGeneralTwoInstance.m

### Unsupervised Alignment

The code is the same as feature-level manifold projections.
The only difference is how to create the correspondence matrix.

 * generateWeight3.m, used to generate weight matrix, calls:
   * computeOptimalMatch.m
   * decompose3.m 

## Folder: Utility

cmpEmbedding.m - compares different embedding results.
Examples are represented as columns.

knnsearch.m - k-nearest neighbor search

LaplacianEigenmaps.m

LPP.m - Locality Preserving Projections

Showtopics.m

createAllConnectedGraph.m

createKnnGraph.m

L2_distance.m - All pairs Euclidean distance

