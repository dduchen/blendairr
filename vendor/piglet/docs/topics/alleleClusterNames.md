**alleleClusterNames** - *Allele similarity cluster naming scheme*

Description
--------------------

For a given cluster the function collapse similar sequences and renames the sequences based on the ASC name scheme


Usage
--------------------
```
alleleClusterNames(cluster, allele.cluster.table, germ.dist, chain, segment)
```

Arguments
-------------------

cluster
:   A vector with the cluster identifier - the family and allele cluster number.

allele.cluster.table
:   A data.frame with the list of all germline sequences and their clusters.

germ.dist
:   A matrix with the germline distance between the germline set sequences.

chain
:   A character with the chain identifier: IGH/IGL/IGK/TRB/TRA... (Currently only IGH is supported)

segment
:   A character with the segment identifier: IGHV/IGHD/IGHJ.... (Currently only IGHV is supported)




Value
-------------------

A data.frame with the clusters renamed alleles based on the ASC scheme.









