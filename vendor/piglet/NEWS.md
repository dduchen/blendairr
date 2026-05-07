Version 1.3.0:  Feb 2026
-------------------------------------------------------------------------------

API CHANGES:

+ Column renames in `AlleleClusterTable` with backward-compatible deprecation warnings:
  - `imgt_allele` → `iuis_allele` (IUIS is the correct naming authority; IMGT is a database)
  - `Family` → `family` (lowercase snake_case for consistency)
  - `Allele_Cluster` → `allele_cluster` (lowercase snake_case for consistency)

+ Backward compatibility: accessing old column names via `$` or `[[` on an
  `AlleleClusterTable` object emits a deprecation warning and redirects to the
  new name. Functions accepting an `alleleClusterTable` as input also silently
  rename old columns via an internal `.compat_allele_table()` shim.

+ `assignAlleleClusters()`: default `from_col` changed from `"imgt_allele"` to
  `"iuis_allele"`.

+ Output columns of `inferGenotypeAllele_asc()` renamed: `imgt_alleles` →
  `iuis_alleles`, `genotyped_imgt_alleles` → `genotyped_iuis_alleles`.

NEW FEATURES:

+ `GermlineCluster` objects now store the raw family cut data.frame as the
  hidden slot `.familiesCut` (columns `iuis_allele` and `family`), produced
  at the hierarchical-clustering step before merging with allele-cluster
  assignments.

Version 1.2.1:  Feb 2026
-------------------------------------------------------------------------------

BUG FIXES:

+ Leiden clustering now uses hierarchical clustering for Family assignment,
  matching the behavior of the reference script. Previously, Family was set
  equal to Allele_Cluster (both from Leiden), losing the hierarchical family
  structure. Now Family comes from hierarchical clustering and Allele_Cluster
  from Leiden community detection.

+ Leiden results now populate `hclustAlleleCluster` and `family_threshold`,
  enabling `plotTruncatedTree()` to work with Leiden-based results.

Version 1.2.0:  Feb 2026
-------------------------------------------------------------------------------

MAJOR CHANGES:

+ Leiden Community Detection: Added support for Leiden algorithm as an alternative
  to hierarchical clustering. Use `clustering_method = "leiden"` in `inferAlleleClusters()`.
  This is particularly useful for D and J segments with variable-length sequences.

+ Generalized Distance Functions: The new `igDistance()` function supports multiple
  distance calculation methods:
  - "decipher": DECIPHER-based distance (default, best for aligned V segments)
  - "hamming": Hamming distance (requires equal-length sequences)
  - "lv": Levenshtein distance (handles variable-length sequences, good for D/J)

+ Multi-locus Support: `inferAlleleClusters()` now explicitly supports all V/D/J
  segments with the `locus` parameter.

+ New Visualization Functions:
  - `plotTruncatedTree()`: Publication-quality circular tree with collapsed subgroups
  - `plotCommunityNetwork()`: Network visualization of cluster relationships
  - `plotSilhouetteOptimization()`: Resolution parameter optimization plot
  - `plotClusterComparison()`: Compare hierarchical and Leiden clustering results

CLASS SYSTEM CHANGES:

+ S4 to S3 Migration: `GermlineCluster` is now an S3 class instead of S4.
  - Access slots with `$` instead of `@` (e.g., `result$alleleClusterTable`)
  - Added `print()`, `summary()`, and `plot()` S3 methods
  - New slots: `clusteringMethod`, `communityObject`, `graphObject`, `distanceMatrix`,
    `silhouetteScore`, `resolutionParameter`, `locus`

API CHANGES:

+ `ighvDistance()` is deprecated in favor of `igDistance()` (backward compatible wrapper provided)
+ `ighvClust()` is deprecated in favor of `igClust()` (backward compatible wrapper provided)
+ `inferAlleleClusters()` has new parameters but maintains backward compatibility

DEPENDENCIES:

+ Added: `igraph`, `stringdist`, `cluster`, `ape`
+ Moved to Suggests: `ComplexHeatmap`, `dplyr`, `ggtree`
+ Removed from Imports: `dplyr` (replaced with `data.table` operations)

MINOR CHANGES:

+ Improved CRAN compliance with proper `Author` and `Maintainer` fields
+ Added test suite using `testthat` (edition 3)
+ Redesigned Shiny app with modular architecture
+ Added Conda recipe for bioconda distribution

BUG FIXES:

+ Fixed deprecated `dplyr::between()` calls
+ Fixed deprecated `dplyr::bind_rows()` calls
+ Replaced deprecated S4 slot access patterns


Version 1.0.7:  Mar 19, 2025
-------------------------------------------------------------------------------

+ Working version for CRAN.

Version 1.0.6:  Nov 7, 2024
-------------------------------------------------------------------------------

+ An asc_annotation flag was added to `inferGenotypeAllele` to simplify genotype by asc annotation.


Version 1.0.5:  Oct 1, 2024
-------------------------------------------------------------------------------

+ New Feature! We added a confidence level to genotype inference `genotype_confidence`.

The confidence is calculated as the Z score of the allele count, repertoire depth, and the allele absolute threshold.

Version 1.0.4:  April 24, 2024
-------------------------------------------------------------------------------

+ Fixed bug in `download_zenodo_files`.

Version 1.0.3:  December 13, 2023
-------------------------------------------------------------------------------

+ Fixed bug in `download_zenodo_files`.

Version 1.0.1:  September 13, 2023
-------------------------------------------------------------------------------

+ Fixed bug in `download_zenodo_files`.

Version 1.0.1:  August 16, 2023
-------------------------------------------------------------------------------

+ First CRAN release.

Version 1.0.0:  June 9, 2023
-------------------------------------------------------------------------------

+ First release.

Version 0.0.1:  December 7, 2022
-------------------------------------------------------------------------------

+ Stable version.
