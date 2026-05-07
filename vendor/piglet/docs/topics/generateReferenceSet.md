**generateReferenceSet** - *Generate allele similarity reference set*

Description
--------------------

Generates the allele clusters reference set based on the clustering from [ighvClust](ighvClust.md). The function collapse
similar alleles and assign them into their respective allele clusters and family clusters. See details for naming scheme


Usage
--------------------
```
generateReferenceSet(germline_distance, germline_set, alleleClusterTable)
```

Arguments
-------------------

germline_distance
:   A germline set distance matrix created by [ighvDistance](ighvDistance.md).

germline_set
:   A character list of the IMGT aligned IGHV allele sequences. See details for curating options.

alleleClusterTable
:   A data.frame of the alleles and their clusters created by [ighvClust](ighvClust.md).




Value
-------------------

A `list` with the re-named germline set, and a table of the allele clusters and thresholds.


Details
-------------------

Each allele is named by this scheme:
IGHVF1-G1*01 - IGH = chain, V = region, F1 = family cluster numbering,
G1 - allele cluster numbering, and 01 = allele numbering (given by clustering order, no connection to the expression)









